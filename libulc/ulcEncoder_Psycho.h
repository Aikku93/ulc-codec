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
#define MASKING_BAND_DECIMATION_FACTOR 8

/**************************************/

//! Compute masking energy estimation
static inline float Block_Transform_ComputeMaskingPower(
	const float *Coef,
	const float *CoefNp,
	size_t Band,
	size_t BlockSize,
	float NyquistHz
) {
#if ULC_USE_PSYHOACOUSTICS
	size_t n;

	//! These settings are mostly based on trial and error
	float Bw = MaskingBandwidth((Band+MASKING_BAND_DECIMATION_FACTOR*0.5f) * NyquistHz/BlockSize) * BlockSize/NyquistHz;
	float BandPow  = 0.0f; for(n=0;n<MASKING_BAND_DECIMATION_FACTOR;n++) BandPow += ABS(Coef[Band+n]);
	      BandPow  = SQR(BandPow / MASKING_BAND_DECIMATION_FACTOR);
	float MaskFW   = 1.0f - BandPow;
	float MaskBW   = sqrtf(1.0f - MaskFW);
	size_t BwFW    = (size_t)(Bw * MaskFW);
	size_t BwBW    = (size_t)(Bw * MaskBW);
	if(Band+MASKING_BAND_DECIMATION_FACTOR+BwFW > BlockSize) BwFW = BlockSize - (Band+MASKING_BAND_DECIMATION_FACTOR);
	if(Band                                     < BwBW)      BwBW = Band;
	if(!BwFW && !BwBW) return 0.0f;

	//! Perform convolution
	//! Scaling by a flatness curve seems to give better results,
	//! where this flatness curve doubles the masking power at
	//! "medium" flatness (corresponding to complex signals that
	//! can have MDCT coefficients discarded without too much
	//! loss of fidelity to the original signal after IMDCT), but
	//! also drops the masking power to 0 when very tonal
	//! NOTE: Flatness calculation is inlined to speed things up
	float SumPow  = 0.0f, Scale = 1.0f / (BwFW + BwBW);
	float SumFlat = BandPow * CoefNp[Band];
	for(n=0;n<BwFW;n++) {
		float Val2  = SQR(Coef  [Band+n+MASKING_BAND_DECIMATION_FACTOR]); SumPow += Val2;
		float ValNp =    (CoefNp[Band+n+MASKING_BAND_DECIMATION_FACTOR]);
		SumFlat += Val2 * ValNp;
	}
	for(n=0;n<BwBW;n++) {
		float Val2  = SQR(Coef  [Band-n-1]); SumPow += Val2;
		float ValNp =    (CoefNp[Band-n-1]);
		SumFlat += Val2 * ValNp;
	}
	float Nrm  = BandPow + SumPow;
	float Flat = exp2f((Nrm*logf(Nrm) - 2.0f*SumFlat) / (Nrm * logf(BwBW+1+BwFW))) - 1.0f;
	return SumPow*Scale * exp2f(1.0f - SQR(5.0f*(Flat - 0.6f)));
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
