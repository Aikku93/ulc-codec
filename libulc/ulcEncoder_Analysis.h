/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
/**************************************/

//! Analysis key structure
#define MIN_BANDS 8
#define MAX_BANDS 65536
#define MIN_CHANS 1
#define MAX_CHANS 65536
struct AnalysisKey_t {
	//! Awkward structuring for faster sorting comparisons
	union {
		struct {
			uint16_t Band; //! Band index
			uint16_t Chan; //! Channel index
		};
		uint32_t SortKey;
	};
	float Val; //! Analysis value (eg. for biased preferencing)
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
//! PONDER:
//!  If we combine the sort key just right (eg. to compact all bits
//!  and leave the topmost bits as 0), then a radix sort might be
//!  better here.
static int Analysis_KeysSort_Comparator(const void *a, const void *b) {
	return ((struct AnalysisKey_t*)a)->SortKey - ((struct AnalysisKey_t*)b)->SortKey;
}
static void Analysis_KeysSort(struct AnalysisKey_t *Keys, size_t nKeys) {
	qsort(Keys, nKeys, sizeof(struct AnalysisKey_t), Analysis_KeysSort_Comparator);
}

/**************************************/
//! EOF
/**************************************/
