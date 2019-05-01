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

//! Bands per spectral flatness measure
#define FLATNESS_WIDTH 32

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
	     if(Fc <=   500.0f) k = 100.0f + 100.0f*SplineCurve( Fc         * (1.0f/  500.0f));
	else if(Fc <= 20000.0f) k = 250.0f - 200.0f*SplineCurve((Fc-500.0f) * (1.0f/19500.0f));
	else                    k =  50.0f;
	k *= CurveSpread;

	//! Rescale to bands. If outside of limits, just return 0, as the
	//! cosine approximation breaks outside of the x=[0,1] range
	k *= BandsPerHz;
	if(k <= 1.0f) return 0.0f;

	//! Return oscillator parameter for differential equation
	float s, c;
	Fourier_SinCos(1.0f / k, &s, &c);
	return 2.0f*c;
}
static void Block_Transform_ComputeMaskingPower(const float *Coef, float *MaskingPower, size_t BlockSize, float Nyquist_Hz) {
	size_t i, j;

	//! Convolve with spreading function
	//! The point of this particular psychoacoustic model
	//! is that neighbouring bands 'support' a given band
	for(i=0;i<BlockSize;i+=2) {
		float Fc = (i+1.0f) * Nyquist_Hz / BlockSize;

		//! Convolve LHS, RHS
		//! Based on a raised cosine curve
		float PowSum = 0.0f; {
			float CurveOmg;
			float CurveOld, Curve;

			//! Using linear prediction for the cosine curve
			CurveOmg = Block_Transform_ComputeMaskingPower_CurveParam(BlockSize / Nyquist_Hz, Fc, 0.3f);
			CurveOld = 1.0f, Curve = CurveOmg*0.5f;
			for(j=1;(Curve > 0.0f && j <= i);j++) {
				PowSum += SQR(Curve * Coef[i-j]);

				float t = Curve;
				Curve = Curve*CurveOmg - CurveOld, CurveOld = t;
			}
			CurveOmg = Block_Transform_ComputeMaskingPower_CurveParam(BlockSize / Nyquist_Hz, Fc, 1.0f);
			CurveOld = 1.0f, Curve = CurveOmg*0.5f;
			for(j=2;(Curve > 0.0f && i+j < BlockSize);j++) {
				PowSum += SQR(Curve * Coef[i+j]);

				float t = Curve;
				Curve = Curve*CurveOmg - CurveOld, CurveOld = t;
			}
		}
		MaskingPower[i/2] = PowSum;
	}
}

//! Compute flatness
static void Block_Transform_ComputeFlatness(const float *Coef, float *Flatness, size_t BlockSize) {
	size_t i;
	for(i=0;i<BlockSize;i+=FLATNESS_WIDTH) *Flatness++ = SpectralFlatness(Coef + i, FLATNESS_WIDTH);
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

	//! Start inserting keys
	float Flat_mu   = 0.0f;
	float Flat_Step = 1.0f / FLATNESS_WIDTH;
	float Flat_Cur  = *Flatness++;
	float Flat_Nxt  = *Flatness++;
	for(i=0;i<BlockSize;i++) {
		//! Check that value doesn't collapse to 0 under the smallest quantizer (1.0)
		float v2 = SQR(Coef[i]);
		if(v2 >= SQR(0.5f)) {
			//! Get flatness
			//! NOTE: mu^2 to get a sharper curve
			float Flat = SQR(Flat_mu);
			      Flat = Flat_Cur*(1.0f - Flat) + Flat_Nxt*Flat;

			//! Build and insert key
			Keys[nKeys].Band = i;
			Keys[nKeys].Chan = Chan;
#if USE_PSYHOACOUSTICS
			//! PONDER: Not sure why this scaling is needed
			Keys[nKeys].Val  = v2*AnalysisPower - (1.0f/6.0f)*Flat*MaskingPower[i/2];
#else
			Keys[nKeys].Val  = v2*AnalysisPower;
#endif
			nKeys++;
		}

		//! Step flatness
		Flat_mu += Flat_Step;
		if(Flat_mu >= 1.0f) {
			Flat_mu -= 1.0f;
			Flat_Cur = Flat_Nxt;
			Flat_Nxt = *Flatness++;
		}
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
