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
#include <string.h>
/**************************************/

//! Analysis key structure
#define ANALYSIS_KEY_UNUSED 0x7FFFFFFF
struct AnalysisKey_t {
	union {
		struct {
			uint16_t Band;
			uint16_t Chan;
		};
		uint32_t Key; //! Band | Chan<<16
	};
	union {
		float   Val;   //! Sorting value (ie. Coef^2 * Importance; for biased preferencing)
		int16_t Quant; //! Quantizer (set after calling Block_Encode_BuildQuants() in ulcEncoder_Quantizer.h)
	};
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
	memmove(Keys + InsertIdx+1, Keys + InsertIdx, (nKeys - InsertIdx)*sizeof(struct AnalysisKey_t));
	Keys[InsertIdx] = *Key;
}

/**************************************/

//! Sort keys in increasing Chan>Band order
//! TODO: Replace this sort with something else
//!       Perhaps radix sort on channels, followed by counting sort on bands.
static int Analysis_KeysSort_Comparator(const void *_a, const void *_b) {
	const struct AnalysisKey_t *a = _a;
	const struct AnalysisKey_t *b = _b;
	return a->Key - b->Key;
}
static void Analysis_KeysSort(struct AnalysisKey_t *Keys, size_t nKeys) {
	qsort(Keys, nKeys, sizeof(struct AnalysisKey_t), Analysis_KeysSort_Comparator);
}

/**************************************/
//! EOF
/**************************************/
