/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
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
			uint8_t  Chan;
			uint8_t  QBand;
		};
		uint32_t Key; //! Band | Chan<<16 | QBand<<24
	};
	union {
		float Val;   //! Sorting value (higher value = greater importance)
		float Quant; //! Quantizer scale (set after calling Block_Encode_BuildQuants() in ulcEncoder_Quantizer.h)
	};
};

/**************************************/

//! Sort keys in increasing Chan>Band order
//! NOTE: All keys' QBand members must be set to 0 before calling this
static int Analysis_KeysSort_Comparator(const void *_a, const void *_b) {
	const struct AnalysisKey_t *a = _a;
	const struct AnalysisKey_t *b = _b;
	return a->Key - b->Key;
}
static void Analysis_KeysSort(struct AnalysisKey_t *Keys, int nKeys) {
	qsort(Keys, nKeys, sizeof(struct AnalysisKey_t), Analysis_KeysSort_Comparator);
}

//! Sort keys in decreasing Val order
static int Analysis_KeysValSort_Comparator(const void *_a, const void *_b) {
	const struct AnalysisKey_t *a = _a;
	const struct AnalysisKey_t *b = _b;
	return a->Val < b->Val ? (+1) : (-1);
}
static void Analysis_KeysValSort(struct AnalysisKey_t *Keys, int nKeys) {
	qsort(Keys, nKeys, sizeof(struct AnalysisKey_t), Analysis_KeysValSort_Comparator);
}

/**************************************/
//! EOF
/**************************************/
