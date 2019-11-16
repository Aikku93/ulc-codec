/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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
#define FLATNESS_COUNT 64

/**************************************/

//! Compute flatness bands
static void Block_Transform_ComputeFlatness(const float *Coef, float *Flatness, size_t BlockSize) {
	size_t i, Width = BlockSize / FLATNESS_COUNT;
	for(i=0;i<BlockSize;i+=Width) *Flatness++ = SpectralFlatness(Coef + i, Width);
	*Flatness = Flatness[-1]; //! Interpolation
}

//! Estimate masking for each band
static inline __attribute__((always_inline)) float Block_Transform_ComputeMaskingPower_Bandwidth(float Fc, float BandsPerHz, float CurveSpread) {
	//! Based on ERB scale (Moore and Glasberg, 1990; with the assumption of being symmetric), then tapered at the end
	float k;
	     if(Fc < 10000.0f) k = 12.35f + 0.0539695f*Fc;
	else if(Fc < 20000.0f) k = 1091.74f - Fc*0.0539695f;
	else                   k = 12.35f;
	return k*CurveSpread*BandsPerHz;
}
static inline __attribute__((always_inline)) float Block_Transform_ComputeMaskingPower_Convolve(
	size_t Band,
	size_t BlockSize,
	float  Nyquist_Hz,
	float  PowSum,
	const float *Coef,
	int Direction
) {
	size_t N;
	if(Direction == -1) { Coef += Band-1; N = Band-1; }
	if(Direction == +1) { Coef += Band+2; N = BlockSize - (Band+2); } //! {Coef[n..n+1]} are the 'center of interest' bands, so skip them
	if(N > 0 && N < BlockSize) { //! NOTE: (N < BlockSize) relies on unsigned overflow behaviour
		float Fc     = (Band+1.0f) * Nyquist_Hz / BlockSize;
		float Spread = 1.0f;
		if(Direction == -1) {
			//! Back-masking becomes about equal to forward-masking at ~16kHz
			Spread = 0.3f + (0.7f/16000.0f)*Fc;
			if(Spread > 1.0f) Spread = 1.0f;
		}

		//! Determine if it's worth computing masking
		float MaskingBw = Block_Transform_ComputeMaskingPower_Bandwidth(Fc, BlockSize / Nyquist_Hz, Spread);
		if(MaskingBw > 1.0f) {
			//! Get oscillator parameters (using linear prediction)
			//! NOTE: Only approximating the area under the curve, but it
			//! is extremely close for MaskingBw < 8. The square root is
			//! to account for squaring in the sum inside the loop
			float InvMaskingBw = 1.0f / MaskingBw;
			float CurveMul = sqrtf(InvMaskingBw);
			float CurveOmg = 2.0f*cosf((float)M_PI_2 * InvMaskingBw);
			float CurveOld = 1.0f * CurveMul;
			float Curve    = 0.5f * CurveMul * CurveOmg;

			//! Begin convolution
			//! Most of the processing time for this computation will
			//! be spent in this loop, so tried to optimize as best
			//! as possible for the compiler
			for(Coef -= Direction;;) {
				Coef   += Direction;
				PowSum += SQR(Curve * (*Coef));
				if(!--N) break;

				float t = Curve;
				Curve = Curve*CurveOmg - CurveOld, CurveOld = t;
				if(Curve < 0.005f) break; //! Sqrt[2^-15]
			}
		}
	}
	return PowSum;
}
static void Block_Transform_ComputeMaskingPower(const float *Coef, float *MaskingPower, const float *Flatness, size_t BlockSize, float Nyquist_Hz) {
	size_t i;
	float Flat_mu   = 0.0f;
	float Flat_Step = (2.0f*FLATNESS_COUNT) / BlockSize;
	float Flat_Cur  = *Flatness++;
	float Flat_Nxt  = *Flatness++;
	for(i=0;i<BlockSize;i+=2) {
		//! Get flatness
		float Flat = SmoothStep(Flat_mu);
		      Flat = Flat_Cur*(1.0f - Flat) + Flat_Nxt*Flat;

		//! Convolve masking power
		float PowSum = 0.0f;
		PowSum = Block_Transform_ComputeMaskingPower_Convolve(i, BlockSize, Nyquist_Hz, PowSum, Coef, -1);
		PowSum = Block_Transform_ComputeMaskingPower_Convolve(i, BlockSize, Nyquist_Hz, PowSum, Coef, +1);
		MaskingPower[i/2] = PowSum*Flat;

		//! Step flatness
		Flat_mu += Flat_Step;
		if(Flat_mu >= 1.0f) {
			Flat_mu -= 1.0f;
			Flat_Cur = Flat_Nxt;
			Flat_Nxt = *Flatness++;
		}
	}
}

/**************************************/
//! EOF
/**************************************/
