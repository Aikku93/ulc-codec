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

//! Compute masking energy estimation
//! NOTE: Masking is only computed for every 2 coefficients to save processing time
static inline __attribute__((always_inline)) float Block_Transform_ComputeMaskingPower_Convolve(
	size_t Band,
	size_t BlockSize,
	const float *Coef,
	float NyquistHz,
	float Flat
) {
	//! These settings are mostly based on trial and error
	//! Flatter bands have a wider masking bandwidth than diffuse ones
	float Bw = Flat*MaskingBandwidth((Band+1.0f) * NyquistHz/BlockSize) * BlockSize/NyquistHz;
	float BandPow = sqrtf(SQR(Coef[Band]) + SQR(Coef[Band+1])) * (float)M_SQRT1_2;
	float MaskFW  = 1.0f - BandPow;
	float MaskBW  = sqrtf(1.0f - MaskFW);
	size_t BwFW   = (size_t)(Bw * MaskFW); if(Band+2+BwFW > BlockSize) BwFW = BlockSize - (Band+2);
	size_t BwBW   = (size_t)(Bw * MaskBW); if(Band        < BwBW)      BwBW = Band;
	if(!BwFW && !BwBW) return 0.0f;

	//! Perform convolution
	//! Scaling by 1+Flatness seems to give better results
	float Sum = 0.0f, Scale = (1.0f+Flat) / (BwFW + BwBW);
	if(BwFW) do Sum += SQR(Coef[Band+BwFW+1]); while(--BwFW);
	if(BwBW) do Sum += SQR(Coef[Band-BwBW  ]); while(--BwBW);
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
		float Flat = Flat_Cur*(1.0f - Flat_mu) + Flat_Nxt*Flat_mu;
		Flat_mu += Flat_Step;
		if(Flat_mu >= 1.0f) {
			Flat_mu -= 1.0f;
			Flat_Cur = Flat_Nxt;
			Flat_Nxt = *Flatness++;
		}

		//! Convolve masking power
		MaskingPower[i/2] = Block_Transform_ComputeMaskingPower_Convolve(i, BlockSize, Coef, NyquistHz, Flat);
	}
}

/**************************************/
//! EOF
/**************************************/
