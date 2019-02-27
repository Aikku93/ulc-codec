/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
#include <stdint.h>
/**************************************/

//! Analysis key structure
//! NOTE: ceil(log2(BlockSize * nChan)) <= 32
#define ANALYSIS_KEY_MAX_BITS 32
struct AnalysisKey_t {
	uint32_t Key; //! Band | Chan<<log2(BlockSize)
	float    Val; //! Analysis value (eg. for biased preferencing)
};

/**************************************/

//! Insert key to analysis buffers
//! Sorted by highest power (passed in 'v') first
//! NOTE:
//!  -This implies that nKeys must increase in the caller
static void Analysis_KeyInsert(const struct AnalysisKey_t *Key, struct AnalysisKey_t *Keys, size_t nKeys) {
	//! Get insertion index (binary search for yeet points)
	size_t InsertIdx; {
		int Lo = 0, Hi = nKeys-1;
		while(Lo <= Hi) {
			size_t Mid = (Lo+Hi) / 2u;
			if(Key->Val > Keys[Mid].Val)
				Hi = Mid-1;
			else
				Lo = Mid+1;
		}
		InsertIdx = Lo;
	}

	//! Shift and insert
	//! PONDER: Faster way of doing the shifting?
	size_t i;
	for(i=nKeys;i>InsertIdx;i--) Keys[i] = Keys[i-1];
	Keys[InsertIdx] = *Key;
}

//! Remove key from analysis buffer
//! NOTE:
//!  -This implies that nKeys must decrease in the caller
static void Analysis_KeyRemove(size_t Key, struct AnalysisKey_t *Keys, size_t nKeys) {
	//! PONDER: Faster way of doing the shifting?
	size_t i;
	for(i=Key+1;i<nKeys;i++) Keys[i-1] = Keys[i];
}

/**************************************/

//! Sort keys in increasing Chan>Band order
//! eg. Ch=0,Band=0..n, Ch=1,Band=0..m, etc.
//! NOTE: Due to combining the Band,Chan values as compactly as
//!       possible (in terms of bits), we can use a radix sort
static void Analysis_KeysSort(struct AnalysisKey_t *Beg, struct AnalysisKey_t *End, size_t Bit) {
	//! If size <= 1, then already sorted
	if(Beg+1 >= End) return;

	//! Swap unsorted keys until L/R meet
	struct AnalysisKey_t *l = Beg, *r = End-1, TempKey;
	for(;;) {
		//! Find swap positions
		while(l < r && (l->Key & Bit) == 0) l++; //! Find LHS.Bit==1
		while(l < r && (r->Key & Bit) != 0) r--; //! Find RHS.Bit==0

		//! Found swap point?
		if(l < r) {
			TempKey = *l;
			*l++ = *r;
			*r-- = TempKey;
		}

		//! L/R met or crossed over?
		if(l >= r) break;
	}

	//! If L/R meet (odd-sized search area), then fix up LHS pointer as needed
	if((l->Key & Bit) == 0 && l < End) l++;

	//! Sort LHS (Bit==0), RHS (Bit==1)
	if(Bit >>= 1) {
		Analysis_KeysSort(Beg, l, Bit);
		Analysis_KeysSort(l, End, Bit);
	}
}

/**************************************/
//! EOF
/**************************************/
