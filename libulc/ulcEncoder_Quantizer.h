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

//! Build quantizer from weighted average
//! NOTE: Biased by the log base-2 of the maximum possible quantized coefficient
//!       (ie. 7^2). This allows coding up to the maximum range in a quantizer.
//! NOTE: The average is performed over the Neper-domain coefficients, so there is
//!       no need to apply a further logarithm here to get the base-2 logarithm.
static float Block_Encode_BuildQuantizer(float Sum, float Weight) {
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

//! Group coefficients into quantizer zones based on their dynamic range
static void Block_Encode_ProcessQuantizerZones(const struct ULC_EncoderState_t *State, int nKeys) {
	//! Analyze the bands and build quantizer bands as we move along
	int   Key;
	int   Chan = -1;         //! Always start on a dummy/invalid channel to force a channel change
	int   QuantStartKey = 0; //! Always store to a valid key on the first [dummy] coefficient
	float QuantMin    =  100.0f;
	float QuantMax    = -100.0f;
	float QuantSum    = 0.0f;
	float QuantWeight = 0.0f;
	float **Coef   = State->TransformBuffer;
	float **CoefNp = State->TransformNepers;
	struct AnalysisKey_t *Keys = State->AnalysisKeys;

	//! Sort keys in band/channel order before processing
	Analysis_KeysSort(Keys, nKeys);

	//! Begin analysis of quantizer zones
#define STORE_QUANTIZER() \
	Keys[QuantStartKey].Quant = Block_Encode_BuildQuantizer(QuantSum, QuantWeight), \
	QuantStartKey = Key, \
	QuantSum    = 0.0f, \
	QuantWeight = 0.0f
	for(Key=0;Key<nKeys;Key++) {
		//! Moved to another channel?
		struct AnalysisKey_t *CurKey = &Keys[Key];
		if(CurKey->Chan != Chan) {
			STORE_QUANTIZER();
			QuantMin =  100.0;
			QuantMax = -100.0;
			Chan  = CurKey->Chan;
		}

		//! Get this key's data
		int   Band    = CurKey->Band;
		float BandPow = SQR(Coef[Chan][Band]);

		//! Level out of range in this quantizer zone?
		if(BandPow < QuantMin) QuantMin = BandPow;
		if(BandPow > QuantMax) QuantMax = BandPow;
		if(QuantMax/QuantMin > SQR(8)) { //! Maximum quantized range before needing a split
			//! Set the new key at which the next quantizer will be stored
			STORE_QUANTIZER();
			QuantMin = BandPow;
			QuantMax = BandPow;
		}

		//! Accumulate to the quantizer
		QuantSum    += BandPow * CoefNp[Chan][Band];
		QuantWeight += BandPow;
		CurKey->Quant = 0.0f; //! Will be set later if there's a quantizer change starting here
	}
	STORE_QUANTIZER();
#undef STORE_QUANTIZER
}

/**************************************/
//! EOF
/**************************************/
