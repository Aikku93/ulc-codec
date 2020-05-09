/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <math.h>
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

//! Build quantizer from weighted average
//! NOTE: Biased by the log base-2 of the maximum possible quantized coefficient
//!       (ie. 7^2). This allows coding up to the maximum range in a quantizer.
//! NOTE: The average is performed over the Neper-domain coefficients, so there is
//!       no need to apply a further logarithm here to get the base-2 logarithm.
static float Block_Encode_BuildQuantizer(float Sum, float Weight) {
	if(Weight == 0.0f) return 0.0f;

	//! NOTE: `sd` will always be greater than 5 due to the bias
	//! NOTE: Truncate (ie. round down); do NOT round off, as if
	//! we were to round up, then too many coefficients would be
	//! clipped and this distorts the output; it's better to
	//! have reduced precision and avoid clipping.
	int sd = (int)(0x1.675768p2f - 0x1.715476p0f*Sum / Weight); //! 0x1.715476p0 == 1/Ln[2], as input is in natural log
	if(sd < 5) sd = 5; //! Sometimes happens because of overflow?
	if(sd > 5 + 0xE + 0xC) sd = 5 + 0xE + 0xC; //! 5+Eh+Ch = Maximum extended-precision quantizer value (including a bias of 5)
	return (float)(1u << sd);
}

/**************************************/

//! Process keys to work with and build quantizers
//! Returns the number of keys to encode
static inline void Block_Encode_ProcessKeys_QuantAdd(float *Sum, float *Weight, float Val, float w) {
	*Sum    += w*Val;
	*Weight += w;
}
static int Block_Encode_ProcessKeys(const struct ULC_EncoderState_t *State, int nNzMax, int nKeys) {
	int Key;

	//! Key data fetcher
	int Band, Chan, QBand;
	float  Val, wVal;
#define FETCH_KEY_DATA(KeyIdx)      \
	Band  = Keys[KeyIdx].Band,  \
	Chan  = Keys[KeyIdx].Chan,  \
	QBand = Keys[KeyIdx].QBand, \
	wVal  = Keys[KeyIdx].Val,   \
	Val   = ABS(CoefBuffer[Chan][Band])
#define MARK_KEY_UNUSED(KeyIdx) Keys[KeyIdx].Key = ANALYSIS_KEY_UNUSED

	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	int        nChan        = State->nChan;
	float    **CoefBuffer   = State->TransformBuffer, **CoefNepers = State->TransformNepers;
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
	//! Note that this is a weighted average, with the weights being the 'audible' level after masking
	for(Key=0;Key<nNzMax;Key++) {
		//! NOTE: wVal must be > 0 for this to function properly, which it always will be at this point
		FETCH_KEY_DATA(Key);
		Block_Encode_ProcessKeys_QuantAdd(&QuantsSum[Chan][QBand], &QuantsWeight[Chan][QBand], CoefNepers[Chan][Band], wVal);
	}
	for(;Key<nKeys;Key++) Keys[Key].Val *= -1.0f; //! Turn negative ('unused')
	for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<ULC_MAX_QBANDS;QBand++) {
		Quants[Chan][QBand] = Block_Encode_BuildQuantizer(QuantsSum[Chan][QBand], QuantsWeight[Chan][QBand]);
	}

	//! Fit quantizers to keys until they stabilize
	//! NOTE: Reusing Keys[].Val as a check for 'key is in analysis'
	int nPass;
	for(nPass=0;nPass<MAX_QUANTIZER_PASSES;nPass++) {
		int nKeysEncode = 0;
		int nNzMaxCur = nNzMax; //! nNzMax for this pass

		//! Check for any collapsed keys
		if(nKeysEncode < nNzMaxCur) do {
			FETCH_KEY_DATA(nKeysEncode);

			//! Collapse?
			float Qnt = Quants[Chan][QBand];
			if(Val*Qnt < 0.5f) {
				//! Key still present in quantizer?
				if(wVal > 0.0f) {
					//! Remove key from analysis (DON'T rebuild the quantizer)
					Block_Encode_ProcessKeys_QuantAdd(&QuantsSum[Chan][QBand], &QuantsWeight[Chan][QBand], CoefNepers[Chan][Band], -wVal);
					Keys[nKeysEncode].Val = -wVal; //! Turn negative ('unused')
				}

				//! Increase number of keys to analyze
				if(nNzMaxCur < nKeys) nNzMaxCur++;
			} else {
				//! Key not present in quantizer?
				if(wVal < 0.0f) {
					//! Re-add key to analysis (DON'T rebuild the quantizer)
					Block_Encode_ProcessKeys_QuantAdd(&QuantsSum[Chan][QBand], &QuantsWeight[Chan][QBand], CoefNepers[Chan][Band], -wVal);
					Keys[nKeysEncode].Val = -wVal; //! Turn positive ('used')
				}
			}
		} while(++nKeysEncode < nNzMaxCur);

		//! Rebuild the quantizers
		int nQuantChanges = 0;
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
	int nKeysEncode;
	int nDeadKeys = 0;
	for(nKeysEncode=0;nKeysEncode<nNzMax;nKeysEncode++) {
		FETCH_KEY_DATA(nKeysEncode);

		//! Collapse?
		float Qnt = Quants[Chan][QBand];
		if(Val*Qnt < 0.5f) {
			//! Mark key unused and increase number of keys to analyze
			MARK_KEY_UNUSED(nKeysEncode);
			nDeadKeys++;
			if(nNzMax < nKeys) nNzMax++;
		} else {
			Keys[nKeysEncode].QBand = 0; //! Must be cleared before sorting
			Keys[nKeysEncode].Quant = Qnt;
		}
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
