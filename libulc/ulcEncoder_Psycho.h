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
static inline __attribute__((always_inline)) float Block_Transform_ComputeMaskingPower_CurveParam(float BandsPerHz, float Fc, float CurveSpread) {
	float k;
	     if(Fc <  1000.0f) k =  50.0f + 200.0f*SmoothStep( Fc          * (1.0f/ 1000.0f)); //! 50..250Hz
	else if(Fc < 20000.0f) k = 250.0f - 225.0f*SmoothStep((Fc-1000.0f) * (1.0f/19000.0f)); //! -250..25Hz
	else                   k =  25.0f;
	k *= CurveSpread;

	//! Rescale to bands. If outside of limits, just return 0, as the
	//! cosine approximation breaks outside of the x=[0,2] range
	k *= BandsPerHz;
	if(k <= 1.0f) return 0.0f;

	//! Return oscillator parameter for differential equation
	return 2.0f*cosf((float)M_PI_2 / k);
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
		float Band_Norm = (Band+1.0f) / BlockSize;
		float Spread = (Direction == +1) ? 1.0f : (0.3f + 0.7f*Band_Norm);

		//! Get oscillator parameters (for linear prediction)
		//! PONDER: Not sure why scaling is needed here;
		//! it becomes approximately 1/4 after squaring.
		//! Maybe it's got to do with combining two
		//! spectral bands per masking band?
		const float CURVE_SCALE = 1.0f/2.0f;
		float CurveOmg = Block_Transform_ComputeMaskingPower_CurveParam(BlockSize / Nyquist_Hz, Band_Norm * Nyquist_Hz, Spread);
		float CurveOld = 1.0f * CURVE_SCALE;
		float Curve    = 0.5f * CURVE_SCALE * CurveOmg;

		//! Begin convolution
		//! Most of the processing time for this computation will
		//! be spent in this loop, so tried to optimize as best
		//! as possible for the compiler
		if(Curve > 0.0f) for(Coef -= Direction;;) {
			Coef   += Direction;
			PowSum += SQR(Curve * (*Coef));
			if(!--N) break;

			float t = Curve;
			Curve = Curve*CurveOmg - CurveOld, CurveOld = t;
			if(Curve < 0.0f) break;
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
		MaskingPower[i/2] = Flat*PowSum;

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
