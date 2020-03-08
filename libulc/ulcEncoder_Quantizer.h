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
//! NOTE: Currently using Sum[x^2]/Sum[x]. This somewhat favours larger
//!       values which mask lower-power values anyway so it works out better.
//! NOTE: Biased by the mean of possible quantized values (ie. x^2, with x = {1..7}),
//!       and the maximum range is extended by 4 bits for the same reason.
static float Block_Encode_BuildQuantizer(double Pow, double Abs) {
	if(Abs == 0.0) return 0.0f;

	//! `sd` will always be greater than 4 (due to the bias)
	long int sd = lrint(0x1.149A784BCD1B9p2 - log2(Pow / Abs));
	if(sd > 4 + 0xE + 15) sd = 4 + 0xE + 15; //! 4+Eh+15 = Maximum extended-precision quantizer value (plus a bias of 4)
	return exp2f(sd);
}

//! Build quantizer bands
static void Block_Encode_BuildQBands(const float *Coefs, uint16_t *QuantsBw, const float *Flatness, size_t BlockSize, float NyquistHz, float RateKbps) {
	size_t Band;
	size_t BandsRem = BlockSize;
	size_t nQBands = 0;
	size_t QBandBw = 0, QBandNzBw = 0;
	double SumSqr  = 0.0;

	//! RateScale adjusts things based on the target bitrate
	//! so as to avoid creating too many quantizer bands at
	//! very low bitrates (eg. 32kbps, etc.)
	//! The expression is compacted from:
	//!  0.5 * 5*(1 - (RateKbps - 64)/(320*NyquistHz/22050))
	//! Key features:
	//!  - Scale by 0.5: There are up to 48 quantizer bands,
	//!    and there are approximately 24 critical bands in
	//!    human hearing.
	//!  - At 64kbps, we target 5x scaling relative to the
	//!    'peak' rate of 320kbps. At 320kbps, RateScale is
	//!    1.0 (or rather, 0.5 from scaling), corresponding
	//!    to the rate at which the quantizer was tuned.
	//!  - Below 64kbps, slowly ramp towards even larger
	//!    quantizer bandwidths without exploding as would
	//!    be the case if we used eg. 320/Rate.
	//!  - Past 320kbps, scaling gets smaller than 1.0 so
	//!    as to provide finer quantization at high rates.
	//!  - Scale the 'peak' rate of 320kbps according to
	//!    the sampling rate of the audio (eg. 320kbps does
	//!    not make sense for this codec at eg. 8000Hz).
	float RateScale = 2.5f - 172.265625f*(RateKbps-64.0f)/NyquistHz; if(RateScale < 0.5f) RateScale = 0.25f*exp2f(RateScale);
	float Flat_mu   = 0.0f;
	float Flat_Step = (2.0f*FLATNESS_COUNT) / BlockSize;
	float Flat_Cur  = *Flatness++;
	float Flat_Nxt  = *Flatness++;
	for(Band=0;Band<BlockSize;Band++) {
		//! Codeable?
		double vNew = SQR((double)Coefs[Band]);
		if(vNew >= SQR(0.5*ULC_COEF_EPS)) {
			//! Enough bands to decide on a split?
			//! NOTE: Adjust for flatness; 'flatter' bands don't need fine quantization
			float Flat = Flat_Cur*(1.0f - Flat_mu) + Flat_Nxt*Flat_mu;
			size_t QBandBwThres = (size_t)((0.25f + 0.75f*Flat) * RateScale * MaskingBandwidth(Band*NyquistHz/BlockSize)*BlockSize/NyquistHz);
			if(QBandNzBw > QBandBwThres) {
				//! Coefficient not in range?
				//! NOTE: Somewhat arbitrary (though tuned) thresholds (derived from 1/(2*RMS[{1..7}]) and RMS[{1..7}])
				double t = vNew*QBandNzBw;
				if(t < 0x1.3CF53F3A97312p-6*SumSqr || t > 0x1.9D87F87E71422p4*SumSqr) {
					//! Create a split point
					//! NOTE: Last band is built from remaining coefficients
					BandsRem -= QBandBw;
					QuantsBw[nQBands++] = QBandBw;
					if(nQBands == MAX_QUANTS-1) break;

					//! Reset state for a new band
					QBandBw = QBandNzBw = 0;
					SumSqr  = 0.0;
				}
			}

			//! Add to quantizer band
			SumSqr += vNew;
			QBandNzBw++;
		}

		//! Increase bandwidth
		QBandBw++;

		//! Step through flatness
		if(Band%2 == 1) {
			Flat_mu += Flat_Step;
			if(Flat_mu >= 1.0f) {
				Flat_mu -= 1.0f;
				Flat_Cur = Flat_Nxt;
				Flat_Nxt = *Flatness++;
			}
		}
	}
	QuantsBw[nQBands++] = BandsRem;
}

/**************************************/

//! Process keys to work with and build quantizers
//! Returns the number of keys to encode
static size_t Block_Encode_ProcessKeys(const struct ULC_EncoderState_t *State, size_t nNzMax, size_t nKeys, float RateKbps) {
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
		Block_Encode_BuildQBands(CoefBuffer[Chan], QuantsBw[Chan], State->TransformFlatness[Chan], BlockSize, State->RateHz * 0.5f, RateKbps);
	}

	//! Clip to maximum available keys
	//! This can happen if rate control says we can fit more
	//! non-zero bands than are actually present for this block
	if(nNzMax > nKeys) nNzMax = nKeys;

	//! Build initial quantizers by considering
	//! all available keys (ie. the 'first pass')
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
