/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <math.h>
#include <stddef.h>
#include <stdint.h>
/**************************************/
#include "ulcEncoder.h"
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_Helper.h"
/**************************************/

//! Maximum allowed quantizer optimization passes
//! More passes should give better results at the cost of increased complexity
#define MAX_QUANTIZER_PASSES 16

/**************************************/

//! Get quantizer band from band index
static size_t Block_Encode_BuildQuants_GetQBand(size_t Band, const uint16_t *QuantBw) {
	size_t QBand;
	for(QBand=0;;QBand++) {
		size_t Bw = QuantBw[QBand];
		if(Band < Bw) return QBand;
		Band -= Bw;
	}
}

//! Build quantizer from sum of raised-power values and sum of absolutes
//! NOTE: Biased by the median of possible quantized values (ie. x^2, with x = {1..7}),
//!       and the maximum range is reduced by the same amount for compactness.
static float Block_Encode_BuildQuantizer(float Sum, float Weight) {
	if(Weight == 0.0f) return 0.0f;

	//! `sd` will always be greater than 4 (due to the bias)
	//! NOTE: Adding 0.5 and then truncating to avoid lrint()
	int sd = (int)(4.5f - 0x1.715476p0f*logf(Sum / Weight)); //! 0x1.715476p0 == 1/Ln[2]
	if(sd < 0) sd = 0; //! Sometimes happens because of overflow?
	if(sd > 4 + 0xE + 15) sd = 4 + 0xE + 15; //! 4+Eh+15 = Maximum extended-precision quantizer value (including a bias of 4)
	return exp2f(sd);
}

/**************************************/

//! Process keys to work with and build quantizers
//! Returns the number of keys to encode
static inline void Block_Encode_ProcessKeys_QuantAdd(float *Sum, float *Weight, float Val) {
	float w = Val;
	*Sum    += w*Val;
	*Weight += w;
}
static inline void Block_Encode_ProcessKeys_QuantRemove(float *Sum, float *Weight, float Val) {
	float w = Val;
	*Sum    -= w*Val;
	*Weight -= w;
}
static size_t Block_Encode_ProcessKeys(const struct ULC_EncoderState_t *State, size_t nNzMax, size_t nKeys) {
	size_t Key;

	//! Key data fetcher
	size_t Band, Chan, QBand;
	float  Val;
#define FETCH_KEY_DATA(KeyIdx)     \
	Band  = Keys[KeyIdx].Band, \
	Chan  = Keys[KeyIdx].Chan, \
	QBand = Block_Encode_BuildQuants_GetQBand(Band, QuantsBw[Chan]), \
	Val   = ABS(CoefBuffer[Chan][Band])
#define MARK_KEY_UNUSED(KeyIdx) Keys[KeyIdx].Key = ANALYSIS_KEY_UNUSED

	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	size_t     nChan        = State->nChan;
	float    **CoefBuffer   = State->TransformBuffer;
	uint16_t **QuantsBw     = State->QuantsBw;
	float    **Quants       = State->Quants;
	float    **QuantsSum    = State->QuantsSum;
	float    **QuantsWeight = State->QuantsWeight;
	struct AnalysisKey_t *Keys = State->AnalysisKeys;

	//! Clear quantizer state
	for(Chan=0;Chan<nChan;Chan++) {
		for(QBand=0;QBand<ULC_MAX_QBANDS;QBand++) {
			QuantsSum   [Chan][QBand] = 0.0f;
			QuantsWeight[Chan][QBand] = 0.0f;
			//Quants      [Chan][QBand] = 0.0f; //! <- Will be set during first pass below
		}
	}

	//! Clip to maximum available keys
	//! This can happen if rate control says we can fit more
	//! non-zero bands than are actually present for this block
	if(nNzMax > nKeys) nNzMax = nKeys;

	//! Build initial quantizers by considering all available keys (ie. the 'first pass')
	//! Note that this is a weighted average, with the weights being the absolute of the value
	for(Key=0;Key<nNzMax;Key++) {
		FETCH_KEY_DATA(Key);
		Block_Encode_ProcessKeys_QuantAdd(&QuantsSum[Chan][QBand], &QuantsWeight[Chan][QBand], Val);
		Keys[Key].Val = 1.0f; //! Anything non-negative (see below)
	}
	for(;Key<nKeys;Key++) Keys[Key].Val = -1.0f; //! Anything negative (see below)
	for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<ULC_MAX_QBANDS;QBand++) {
		Quants[Chan][QBand] = Block_Encode_BuildQuantizer(QuantsSum[Chan][QBand], QuantsWeight[Chan][QBand]);
	}

	//! Fit quantizers to keys until they stabilize
	//! NOTE: Reusing Keys[].Val as a check for 'key is in analysis'
	size_t nPass;
	for(nPass=0;nPass<MAX_QUANTIZER_PASSES;nPass++) {
		size_t nKeysEncode = 0;
		size_t nNzMaxCur = nNzMax; //! nNzMax for this pass

		//! Check for any collapsed keys
		if(nKeysEncode < nNzMaxCur) do {
			FETCH_KEY_DATA(nKeysEncode);

			//! Collapse?
			float Qnt = Quants[Chan][QBand];
			if(Val*Qnt < 0.5f) {
				//! Key still present in quantizer?
				if(Keys[nKeysEncode].Val >= 0.0f) {
					//! Remove key from analysis (DON'T rebuild the quantizer)
					Block_Encode_ProcessKeys_QuantRemove(&QuantsSum[Chan][QBand], &QuantsWeight[Chan][QBand], Val);
					Keys[nKeysEncode].Val = -1.0f; //! Anything negative
				}

				//! Increase number of keys to analyze
				if(nNzMaxCur < nKeys) nNzMaxCur++;
			} else {
				//! Key not present in quantizer?
				if(Keys[nKeysEncode].Val < 0.0f) {
					//! Re-add key to analysis (DON'T rebuild the quantizer)
					Block_Encode_ProcessKeys_QuantAdd(&QuantsSum[Chan][QBand], &QuantsWeight[Chan][QBand], Val);
					Keys[nKeysEncode].Val = +1.0f; //! Anything non-negative
				}
			}
		} while(++nKeysEncode < nNzMaxCur);

		//! Rebuild the quantizers
		size_t nQuantChanges = 0;
		for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<ULC_MAX_QBANDS;QBand++) {
			float QntOld = Quants[Chan][QBand];
			float QntNew = Block_Encode_BuildQuantizer(QuantsSum[Chan][QBand], QuantsWeight[Chan][QBand]);
			if(QntNew != QntOld) {
				//! Don't count creating or destroying a quantizer
				Quants[Chan][QBand] = QntNew;
				if(QntOld != 0.0f && QntNew != 0.0f) nQuantChanges++;
			}
		}

		//! If converged, early exit
		if(!nQuantChanges) break;
	}

	//! Final pass to remove 'dead' keys
	//! NOTE: Saving quantizers to key value; this avoids a lookup inside the coding loop
	size_t nKeysEncode;
	size_t nDeadKeys = 0;
	for(nKeysEncode=0;nKeysEncode<nNzMax;nKeysEncode++) {
		FETCH_KEY_DATA(nKeysEncode);

		//! Collapse?
		float Qnt = Quants[Chan][QBand];
		if(Val*Qnt < 0.5f) {
			//! Mark key unused and increase number of keys to analyze
			MARK_KEY_UNUSED(nKeysEncode);
			nDeadKeys++;
			if(nNzMax < nKeys) nNzMax++;
		} else Keys[nKeysEncode].Quant = Qnt;
	}

	//! Sort the keys and return how many we work with
	Analysis_KeysSort(Keys, nKeysEncode);
	return nKeysEncode - nDeadKeys;
#undef MARK_KEY_UNUSED
#undef FETCH_KEY_DATA
}

/**************************************/
//! EOF
/**************************************/
