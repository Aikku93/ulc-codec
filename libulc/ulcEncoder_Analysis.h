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

//! Sort keys in increasing Chan>Band order
static int Analysis_KeysSort_Comparator(const void *_a, const void *_b) {
	const struct AnalysisKey_t *a = _a;
	const struct AnalysisKey_t *b = _b;
	return a->Key - b->Key;
}
static void Analysis_KeysSort(struct AnalysisKey_t *Keys, size_t nKeys) {
	qsort(Keys, nKeys, sizeof(struct AnalysisKey_t), Analysis_KeysSort_Comparator);
}

//! Sort keys in decreasing Val order
static int Analysis_KeysValSort_Comparator(const void *_a, const void *_b) {
	const struct AnalysisKey_t *a = _a;
	const struct AnalysisKey_t *b = _b;
	return a->Val < b->Val ? (+1) : (-1);
}
static void Analysis_KeysValSort(struct AnalysisKey_t *Keys, size_t nKeys) {
	qsort(Keys, nKeys, sizeof(struct AnalysisKey_t), Analysis_KeysValSort_Comparator);
}

/**************************************/
//! EOF
/**************************************/
