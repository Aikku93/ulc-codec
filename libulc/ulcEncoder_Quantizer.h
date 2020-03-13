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

//! Maximum allowed quantizers
//! NOTE: NO GREATER THAN 48
#define MAX_QUANTS 48

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
//! NOTE: Biased by the mean of possible quantized values (ie. x^2, with x = {1..7}),
//!       and the maximum range is extended by 4 bits for the same reason.
static float Block_Encode_BuildQuantizer(double Pow, double Abs) {
	if(Abs == 0.0) return 0.0f;

	//! `sd` will always be greater than 4 (due to the bias)
	long int sd = lrint(0x1.149A784BCD1B9p2 - log2(Pow / Abs));
	if(sd > 4 + 0xE + 15) sd = 4 + 0xE + 15; //! 4+Eh+15 = Maximum extended-precision quantizer value (including a bias of 4)
	return exp2f(sd);
}

//! Build quantizer bands
static void Block_Encode_BuildQBands(const float *CoefsNp, uint16_t *QuantsBw, size_t BlockSize) {
	//! Average step size of quantized values, in Nepers (ie. -Ln[7])
	const float AvgMaxRange = 0x1.F2272Bp0;

	//! Split into quantizer bands based on range distortion
	//! For some reason, it works better to average in the log domain
	//! here, while just use a weighted average in the linear domain
	//! when calculating the actual quantizer values
	float  Avg   = 0.0f, nAvg   = 0.0f; //! Initializing all these shuts gcc up
	float  AvgLo = 0.0f, nAvgLo = 0.0f;
	float  AvgHi = 0.0f, nAvgHi = 0.0f;
	size_t Band, nQBands = 0, LastQBandOfs = 0;
	for(Band=0;Band<BlockSize;Band++) {
		float vNp = CoefsNp[Band]; if(vNp == ULC_COEF_NEPER_OUT_OF_RANGE) continue;

		//! New quantizer band?
		if(Avg == 0.0f) {
			Avg   = vNp, nAvg   = 1.0f;
			AvgLo = Avg, nAvgLo = 1.0f;
			AvgHi = Avg, nAvgHi = 1.0f;
		} else {
			//! Sum to lower-/higher-than-average parts
#if 0
			if(vNp < Avg/nAvg) AvgLo += vNp, nAvgLo += 1.0f;
			else               AvgHi += vNp, nAvgHi += 1.0f;
#else //! Math optimization (avoid division)
			if(vNp*nAvg < Avg) AvgLo += vNp, nAvgLo += 1.0f;
			else               AvgHi += vNp, nAvgHi += 1.0f;
#endif
			Avg += vNp, nAvg += 1.0f;

			//! Check range against threshold
#if 0
			float r = ABS(AvgHi/nAvgHi - AvgLo/nAvgLo);
			if(r > AvgMaxRange) {
#else //! Math optimization (avoid division)
			float r = ABS(AvgHi*nAvgLo - AvgLo*nAvgHi);
			if(r > AvgMaxRange*nAvgLo*nAvgHi) {
#endif
				//! A very sharp transient inside a small bandwidth
				//! is likely to be a statistical anomaly, so ignore it
				//! and let it stabilize on its own
				size_t Bw = Band-1 - LastQBandOfs;
				if(Bw > 4) {
					//! Create quantizer band
					LastQBandOfs = Band;
					QuantsBw[nQBands] = Bw;
					if(++nQBands == MAX_QUANTS-1) break;

					//! Reset for a new band starting from the transient band
					Avg = 0.0f;
					Band--;
				}
			}
		}
	}
	QuantsBw[nQBands++] = BlockSize - LastQBandOfs;
}

/**************************************/

//! Process keys to work with and build quantizers
//! Returns the number of keys to encode
static size_t Block_Encode_ProcessKeys(const struct ULC_EncoderState_t *State, size_t nNzMax, size_t nKeys) {
	size_t Key;
	size_t BlockSize = State->BlockSize;

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
	size_t nChan          = State->nChan;
	float    **CoefBuffer = State->TransformBuffer;
	uint16_t **QuantsBw   = State->QuantsBw;
	float    **Quants     = State->Quants;
	double   **QuantsPow  = State->QuantsPow;
	double   **QuantsAbs  = State->QuantsAbs;
	struct AnalysisKey_t *Keys = State->AnalysisKeys;

	//! Clear quantizer state and build quantizer bands
	for(Chan=0;Chan<nChan;Chan++) {
		for(QBand=0;QBand<MAX_QUANTS;QBand++) {
			QuantsPow[Chan][QBand] = 0.0;
			QuantsAbs[Chan][QBand] = 0.0;
			//Quants   [Chan][QBand] = 0.0f; //! <- Will be set during first pass below
		}
		Block_Encode_BuildQBands(State->TransformNepers[Chan], QuantsBw[Chan], BlockSize);
	}

	//! Clip to maximum available keys
	//! This can happen if rate control says we can fit more
	//! non-zero bands than are actually present for this block
	if(nNzMax > nKeys) nNzMax = nKeys;

	//! Build initial quantizers by considering all available keys (ie. the 'first pass')
	//! Note that this is a weighted average, with the weights being the absolute of the value
	for(Key=0;Key<nNzMax;Key++) {
		FETCH_KEY_DATA(Key);

		QuantsPow[Chan][QBand] += SQR((double)Val);
		QuantsAbs[Chan][QBand] += Val;
		Keys[Key].Val = 1.0f; //! Anything non-negative (see below)
	}
	for(;Key<nKeys;Key++) Keys[Key].Val = -1.0f; //! Anything negative (see below)
	for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<MAX_QUANTS;QBand++) {
		Quants[Chan][QBand] = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsAbs[Chan][QBand]);
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
					QuantsPow[Chan][QBand] -= SQR((double)Val);
					QuantsAbs[Chan][QBand] -= Val;
					Keys[nKeysEncode].Val = -1.0f; //! Anything negative
				}

				//! Increase number of keys to analyze
				if(nNzMaxCur < nKeys) nNzMaxCur++;
			} else {
				//! Key not present in quantizer?
				if(Keys[nKeysEncode].Val < 0.0f) {
					//! Re-add key to analysis (DON'T rebuild the quantizer)
					QuantsPow[Chan][QBand] += SQR((double)Val);
					QuantsAbs[Chan][QBand] += Val;
					Keys[nKeysEncode].Val = +1.0f; //! Anything non-negative
				}
			}
		} while(++nKeysEncode < nNzMaxCur);

		//! Rebuild the quantizers
		size_t nQuantChanges = 0;
		for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<MAX_QUANTS;QBand++) {
			float QntOld = Quants[Chan][QBand];
			float QntNew = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsAbs[Chan][QBand]);
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
