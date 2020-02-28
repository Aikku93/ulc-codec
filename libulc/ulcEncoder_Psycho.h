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
#include "ulcEncoder_Helper.h"
/**************************************/

//! Number of spectral flatness bands
#define FLATNESS_COUNT 32

/**************************************/

//! Compute flatness bands
static void Block_Transform_ComputeFlatness(const float *Coef, float *Flatness, size_t BlockSize) {
	size_t i, Width = BlockSize / FLATNESS_COUNT;
	for(i=0;i<BlockSize;i+=Width) *Flatness++ = SpectralFlatness(Coef + i, Width);
	*Flatness = Flatness[-1]; //! Interpolation
}

/**************************************/

static inline __attribute__((always_inline)) float Block_Transform_ComputeMaskingPower_Convolve(
	size_t Band,
	size_t BlockSize,
	const float *Coef,
	float NyquistHz
) {
	//! These settings are mostly based on trial and error
	float Bw = MaskingBandwidth((Band+1.0f) * NyquistHz/BlockSize) * BlockSize/NyquistHz;
	float BandPow = sqrtf(SQR(Coef[Band]) + SQR(Coef[Band+1])) / (float)(M_SQRT2 * 32768.0);
	float MaskFW  = 1.0f - BandPow;
	float MaskBW  = sqrtf(1.0f - MaskFW);
	size_t BwFW   = (size_t)(Bw * MaskFW); if(Band+2+BwFW > BlockSize) BwFW = BlockSize - (Band+2);
	size_t BwBW   = (size_t)(Bw * MaskBW); if(Band        < BwBW)      BwBW = Band;
	if(!BwFW && !BwBW) return 0.0f;

	//! At this point, the coefficients are scaled by 0.5 in preparation for
	//! the upcoming sum/difference transform, so as to normalize the final
	//! coefficients. However, masking energy is subtracted /after/ that
	//! transform, so it must be scaled by 2.0 here since it was calculated
	//! when the coefficients were still scaled by 0.5.
	float Sum = 0.0f, Scale = 2.0f / (BwFW + BwBW);
#if 0 //! No intercarrier interference compensation; can result in rumbling and birdies
	if(BwFW) do Sum += SQR(Coef[Band+BwFW+1]); while(--BwFW);
	if(BwBW) do Sum += SQR(Coef[Band-BwBW  ]); while(--BwBW);
#else //! Rudimentary protection
	size_t i;
	static const size_t N_PHASECOMPENSATION = 8;
	static const float PhaseCompensation[] = { //! 4^-(n+1)
		0.25, 0.0625, 0.015625, 0.00390625, 0.000976563, 0.000244141, 0.0000610352, 0.0000152588
	};
	for(i=0;i<BwFW;i++) Sum += SQR(Coef[Band+i+2]) - ((i < N_PHASECOMPENSATION) ? BandPow*PhaseCompensation[i] : 0.0f);
	for(i=0;i<BwBW;i++) Sum += SQR(Coef[Band-i-1]) - ((i < N_PHASECOMPENSATION) ? BandPow*PhaseCompensation[i] : 0.0f);
#endif
	return Sum*Scale;
}
static void Block_Transform_ComputeMaskingPower(
	const float *Coef,
	float *MaskingPower,
	const float *Flatness,
	size_t BlockSize,
	float NyquistHz
) {
	size_t i;
	float Flat_mu   = 0.0f;
	float Flat_Step = (2.0f*FLATNESS_COUNT) / BlockSize;
	float Flat_Cur  = *Flatness++;
	float Flat_Nxt  = *Flatness++;
	for(i=0;i<BlockSize;i+=2) {
		//! Get flatness and step through
		//! Mapping flatness through an inverse exponential curve seems to give better results
		float Flat = Flat_Cur*(1.0f - Flat_mu) + Flat_Nxt*Flat_mu;
		      Flat = 1.0f - expf(-2.0f*(float)M_PI * Flat);
		Flat_mu += Flat_Step;
		if(Flat_mu >= 1.0f) {
			Flat_mu -= 1.0f;
			Flat_Cur = Flat_Nxt;
			Flat_Nxt = *Flatness++;
		}

		//! Convolve masking power
		//! Flatter (noisy) bands contribute more masking than sharper (tonal) ones
		MaskingPower[i/2] = Block_Transform_ComputeMaskingPower_Convolve(i, BlockSize, Coef, NyquistHz)*Flat;
	}
}

/**************************************/
//! EOF
/**************************************/
