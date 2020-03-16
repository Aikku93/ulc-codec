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
#include "ulcEncoder_Helper.h"
/**************************************/

//! How many masking bands are used per actual band
//! 1 corresponds to no decimation. This would ideally give
//! the most correct results (while being the slowest option),
//! but in testing it was found that 16 or 32 works better in
//! practice, while also being significantly faster
#define MASKING_BAND_DECIMATION_FACTOR 32

/**************************************/

//! Compute masking energy estimation
//! NOTE: This uses an inlined version of SpectralFlatness(), as we
//! are already iterating over those bands so it makes no sense to
//! use multiple loops to do different tasks
static inline float Block_Transform_ComputeMaskingPower(
	const float *Coef,
	const float *CoefNp,
	size_t Band,
	size_t BlockSize,
	float NyquistHz
) {
#if ULC_USE_PSYCHOACOUSTICS
	size_t n;

	//! These settings are mostly based on trial and error
	float Bw = MaskingBandwidth((Band+MASKING_BAND_DECIMATION_FACTOR*0.5f) * NyquistHz/BlockSize) * BlockSize/NyquistHz;
	float BandPow  = 0.0f; for(n=0;n<MASKING_BAND_DECIMATION_FACTOR;n++) BandPow += ABS(Coef[Band+n]);
	      BandPow  = BandPow / MASKING_BAND_DECIMATION_FACTOR;
	float MaskFW   = 1.0f - BandPow;
	float MaskBW   = sqrtf(1.0f - MaskFW);
	float fBwFW    = Bw * MaskFW;
	float fBwBW    = Bw * MaskBW;
	size_t BwFW    = (size_t)(fBwFW);
	size_t BwBW    = (size_t)(fBwBW);
	if(Band+MASKING_BAND_DECIMATION_FACTOR+BwFW > BlockSize) BwFW = BlockSize - (Band+MASKING_BAND_DECIMATION_FACTOR);
	if(Band                                     < BwBW)      BwBW = Band;
	if(!BwFW && !BwBW) return 0.0f;

	//! Get the energy of the masked bandwidths
	float SumPow = 0.0f, SumFlat = 0.0f; {
		float ScaleLin, ScaleExp;
		float DecayLin, DecayExp;

		//! Forwards
		//! Target decay of -1.6Np (~0.2 linear scaling)
		ScaleLin = 0.0f;
		ScaleExp = 1.0f;
		DecayLin = -1.6f / fBwFW;
		DecayExp = expf(DecayLin);
		for(n=0;n<BwFW;n++) {
			float Val2  = ABS(Coef  [Band+n+MASKING_BAND_DECIMATION_FACTOR]) * (ScaleExp *= DecayExp);
			float ValNp =    (CoefNp[Band+n+MASKING_BAND_DECIMATION_FACTOR]) + (ScaleLin += DecayLin);
			SumPow  += Val2;
			SumFlat += Val2 * ValNp;
		}

		//! Backwards
		//! Target decay of -0.7Np (~0.5 linear scaling)
		ScaleLin = 0.0f;
		ScaleExp = 1.0f;
		DecayLin = -1.6f / fBwBW;
		DecayExp = expf(DecayLin);
		for(n=0;n<BwBW;n++) {
			float Val2  = ABS(Coef  [Band-n-1]) * (ScaleExp *= DecayExp);
			float ValNp =    (CoefNp[Band-n-1]) + (ScaleLin += DecayLin);
			SumPow  += Val2;
			SumFlat += Val2*ValNp;
		}
	}

	//! Add the energy of the 'center' bandwidths to get the
	//! final masking power and flatness
	for(n=0;n<MASKING_BAND_DECIMATION_FACTOR;n++) {
		float Val2  = ABS(Coef  [Band+n]); SumPow  += Val2;
		float ValNp =    (CoefNp[Band+n]); SumFlat += Val2*ValNp;
	}
	if(SumPow == 0.0f) return 0.0f;
	float Flat = exp2f((logf(SumPow) - SumFlat/SumPow) / logf(BwBW + MASKING_BAND_DECIMATION_FACTOR + BwFW)) - 1.0f;

	//! Return the final masked power
	//! Scaling by a flatness curve seems to give better results,
	//! as flatter signals mask much more than tonal ones
	return SumFlat/SumPow + 2.5f-4.5f*(1.0f-Flat);
#else
	(void)Coef;
	(void)CoefNp;
	(void)Band;
	(void)BlockSize;
	(void)NyquistHz;
	return 0.0f;
#endif
}

/**************************************/
//! EOF
/**************************************/
