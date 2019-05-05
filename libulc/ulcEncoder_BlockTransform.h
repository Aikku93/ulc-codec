/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#if defined(__AVX__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include <math.h>
#include <stddef.h>
#include <stdint.h>
/**************************************/
#include "ulcUtility.h"
/**************************************/
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_Helper.h"
/**************************************/

//! Use psychoacoustic model
#define USE_PSYHOACOUSTICS 1

//! Number of spectral flatness bands
#define FLATNESS_COUNT 64

/**************************************/

//! Copy scaled samples to buffer
static void Block_Transform_CopySamples(float *DataDst, const float *DataSrc, size_t N, float Scale) {
	size_t i;
#if defined(__AVX__)
	__m256 mScale = _mm256_set1_ps(Scale);
	for(i=0;i<N;i+=64) {
		__m256 v0 = _mm256_load_ps(DataSrc + i+ 0);
		__m256 v1 = _mm256_load_ps(DataSrc + i+ 8);
		__m256 v2 = _mm256_load_ps(DataSrc + i+16);
		__m256 v3 = _mm256_load_ps(DataSrc + i+24);
		__m256 v4 = _mm256_load_ps(DataSrc + i+32);
		__m256 v5 = _mm256_load_ps(DataSrc + i+40);
		__m256 v6 = _mm256_load_ps(DataSrc + i+48);
		__m256 v7 = _mm256_load_ps(DataSrc + i+56);
		v0 = _mm256_mul_ps(v0, mScale);
		v1 = _mm256_mul_ps(v1, mScale);
		v2 = _mm256_mul_ps(v2, mScale);
		v3 = _mm256_mul_ps(v3, mScale);
		v4 = _mm256_mul_ps(v4, mScale);
		v5 = _mm256_mul_ps(v5, mScale);
		v6 = _mm256_mul_ps(v6, mScale);
		v7 = _mm256_mul_ps(v7, mScale);
		_mm256_store_ps(DataDst + i+ 0, v0);
		_mm256_store_ps(DataDst + i+ 8, v1);
		_mm256_store_ps(DataDst + i+16, v2);
		_mm256_store_ps(DataDst + i+24, v3);
		_mm256_store_ps(DataDst + i+32, v4);
		_mm256_store_ps(DataDst + i+40, v5);
		_mm256_store_ps(DataDst + i+48, v6);
		_mm256_store_ps(DataDst + i+56, v7);
	}
#elif defined(__SSE__)
	__m128 mScale = _mm_set1_ps(Scale);
	for(i=0;i<N;i+=16) {
		__m128 v0 = _mm_load_ps(DataSrc + i+ 0);
		__m128 v1 = _mm_load_ps(DataSrc + i+ 4);
		__m128 v2 = _mm_load_ps(DataSrc + i+ 8);
		__m128 v3 = _mm_load_ps(DataSrc + i+12);
		__m128 v4 = _mm_load_ps(DataSrc + i+16);
		__m128 v5 = _mm_load_ps(DataSrc + i+20);
		__m128 v6 = _mm_load_ps(DataSrc + i+24);
		__m128 v7 = _mm_load_ps(DataSrc + i+28);
		v0 = _mm_mul_ps(v0, mScale);
		v1 = _mm_mul_ps(v1, mScale);
		v2 = _mm_mul_ps(v2, mScale);
		v3 = _mm_mul_ps(v3, mScale);
		v4 = _mm_mul_ps(v4, mScale);
		v5 = _mm_mul_ps(v5, mScale);
		v6 = _mm_mul_ps(v6, mScale);
		v7 = _mm_mul_ps(v7, mScale);
		_mm_store_ps(DataDst + i+ 0, v0);
		_mm_store_ps(DataDst + i+ 4, v1);
		_mm_store_ps(DataDst + i+ 8, v2);
		_mm_store_ps(DataDst + i+12, v3);
		_mm_store_ps(DataDst + i+16, v4);
		_mm_store_ps(DataDst + i+20, v5);
		_mm_store_ps(DataDst + i+24, v6);
		_mm_store_ps(DataDst + i+28, v7);
	}
#else
	for(i=0;i<N;i++) DataDst[i] = DataSrc[i] * Scale;
#endif
}

/**************************************/
#if USE_PSYHOACOUSTICS
/**************************************/

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
	float s, c;
	Fourier_SinCos(1.0f / k, &s, &c);
	return 2.0f*c;
}
static inline __attribute__((always_inline)) float Block_Transform_ComputeMaskingPower_Convolve(size_t Band, size_t BlockSize, float Nyquist_Hz, float PowSum, const float *Coef, int Direction) {
	size_t N;
	if(Direction == -1) { Coef += Band-1; N = Band-1; }
	if(Direction == +1) { Coef += Band+2; N = BlockSize - (Band+2); } //! {Coef[n..n+1]} are the 'center of interest' bands, so skip them
	if(N > 0 && N < BlockSize) { //! NOTE: (N < BlockSize) relies on unsigned overflow behaviour
		float Band_Norm = (Band+1.0f) / BlockSize;
		float Spread = (Direction == +1) ? 1.0f : (0.3f + 0.7f*Band_Norm);

		//! Get oscillator parameters (for linear prediction)
		//! PONDER: Not sure why scaling is needed here;
		//! it becomes approximately 1/6 after squaring
		const float CURVE_SCALE = 1.0f/2.5f;
		float CurveOmg = Block_Transform_ComputeMaskingPower_CurveParam(BlockSize / Nyquist_Hz, Band_Norm * Nyquist_Hz, Spread);
		float CurveOld = 1.0f * CURVE_SCALE;
		float Curve    = 0.5f * CURVE_SCALE * CurveOmg;

		//! Begin convolution
		//! Most of the processing time for this computation will
		//! be spent in this loop, so tried to optimize as best
		//! as possible for the compiler
		if(Curve > 0.0f) for(Coef -= Direction;;) {
			Coef += Direction;
			PowSum += SQR(Curve * (*Coef));
			if(!--N) break;

			float t = Curve;
			Curve = Curve*CurveOmg - CurveOld, CurveOld = t;
			if(Curve < 0.0f) break;
		}
	}
	return PowSum;
}
static void Block_Transform_ComputeMaskingPower(const float *Coef, float *MaskingPower, size_t BlockSize, float Nyquist_Hz) {
	size_t i;
	for(i=0;i<BlockSize;i+=2) {
		float PowSum = 0.0f;
		PowSum = Block_Transform_ComputeMaskingPower_Convolve(i, BlockSize, Nyquist_Hz, PowSum, Coef, -1);
		PowSum = Block_Transform_ComputeMaskingPower_Convolve(i, BlockSize, Nyquist_Hz, PowSum, Coef, +1);
		MaskingPower[i/2] = PowSum;
	}
}

//! Compute flatness
static void Block_Transform_ComputeFlatness(const float *Coef, float *Flatness, size_t BlockSize) {
	size_t i, Width = BlockSize / FLATNESS_COUNT;
	for(i=0;i<BlockSize;i+=Width) *Flatness++ = SpectralFlatness(Coef + i, Width);
	*Flatness = Flatness[-1]; //! Interpolation
}

/**************************************/
#endif
/**************************************/

//! Insert keys for block coefficients
//! Returns updated number of keys in list
//! NOTE:
//!  AnalysisPower is used to alter the preference for the currently-being-analyzed channel
static size_t Block_Transform_InsertKeys(const float *Coef, size_t BlockSize, size_t Chan, struct AnalysisKey_t *Keys, size_t nKeys, float AnalysisPower, const float *MaskingPower, const float *Flatness) {
	size_t i;
#if !USE_PSYCHOACOUSTICS
	(void)MaskingPower;
	(void)Flatness;
#endif
	//! Start inserting keys
#if USE_PSYHOACOUSTICS
	float Flat_mu   = 0.0f;
	float Flat_Step = (float)FLATNESS_COUNT / BlockSize;
	float Flat_Cur  = *Flatness++;
	float Flat_Nxt  = *Flatness++;
#endif
	for(i=0;i<BlockSize;i++) {
		//! Check that value doesn't collapse to 0 under the smallest quantizer (1.0)
		float v2 = SQR(Coef[i]);
		if(v2 >= SQR(0.5f)) {
#if USE_PSYHOACOUSTICS
			//! Get flatness
			float Flat = SmoothStep(Flat_mu);
			      Flat = Flat_Cur*(1.0f - Flat) + Flat_Nxt*Flat;
#endif
			//! Build and insert key
			Keys[nKeys].Band = i;
			Keys[nKeys].Chan = Chan;
#if USE_PSYHOACOUSTICS
			Keys[nKeys].Val  = v2*AnalysisPower - Flat*MaskingPower[i/2];
#else
			Keys[nKeys].Val  = v2*AnalysisPower;
#endif
			nKeys++;
		}
#if USE_PSYHOACOUSTICS
		//! Step flatness
		Flat_mu += Flat_Step;
		if(Flat_mu >= 1.0f) {
			Flat_mu -= 1.0f;
			Flat_Cur = Flat_Nxt;
			Flat_Nxt = *Flatness++;
		}
#endif
	}
	return nKeys;
}

/**************************************/

//! Apply block transform
//!  -Fetches data
//!  -Applies MDCT
//!  -Applies anti-pre-echo formula
//!  -Stores keys for block coefficients
//! Returns the number of keys stored
static size_t Block_Transform(const struct ULC_EncoderState_t *State, const float *Data, float PowerDecay) {
	size_t nChan     = State->nChan;
	size_t BlockSize = State->BlockSize;

	size_t Chan;
	size_t nKeys = 0;
	float  AnalysisPower = 1.0f;
	for(Chan=0;Chan<nChan;Chan++) {
		//! Get buffer pointers
		//! NOTE:
		//!  TransformTemp has 2*BlockSize elements
		//!  The first half will contain our (scaled) sample data
		//!  The second half will be for processing
		float *BufferSample    = State->TransformTemp;
		float *BufferTemp      = State->TransformTemp + BlockSize;
		float *BufferMasking   = BufferTemp;
		float *BufferFlatness  = BufferTemp + BlockSize/2;
		float *BufferTransform = State->TransformBuffer[Chan];
		float *BufferFwdLap    = State->TransformFwdLap[Chan];

		//! Fetch sample data
		//! Pre-scale for scaled IMDCT(*2.0/BlockSize) and SumDif transform(*0.5)
		Block_Transform_CopySamples(BufferSample, Data + Chan*BlockSize, BlockSize, 2.0f/BlockSize * 0.5f);

		//! Apply transforms
		//! NOTE: Masking power stored to BufferTemp
		Fourier_MDCT(BufferTransform, BufferSample, BufferFwdLap, BufferTemp, BlockSize, State->BlockOverlap);
#if USE_PSYHOACOUSTICS
		Block_Transform_ComputeMaskingPower(BufferTransform, BufferMasking, BlockSize, State->RateHz*0.5f);
		Block_Transform_ComputeFlatness(BufferTransform, BufferFlatness, BlockSize);
#endif
		ULC_Transform_AntiPreEcho(BufferTransform, BlockSize);

		//! Insert coefficient keys
		nKeys = Block_Transform_InsertKeys(BufferTransform, BlockSize, Chan, State->AnalysisKeys, nKeys, AnalysisPower, BufferMasking, BufferFlatness);
		AnalysisPower *= PowerDecay;
	}
	Analysis_KeysValSort(State->AnalysisKeys, nKeys);
	return nKeys;
}

/**************************************/
//! EOF
/**************************************/
