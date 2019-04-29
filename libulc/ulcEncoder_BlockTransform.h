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
#if defined(__SSE3__)
# include <pmmintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include <stddef.h>
#include <stdint.h>
/**************************************/
#include "ulcUtility.h"
/**************************************/
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_Helper.h"
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

//! Compute adjusted power for bands
//! NOTE: N = BlockSize/2
static inline __attribute__((always_inline)) float Block_Transform_ComputeFinalPower_SupportCurveParam(float BandsPerHz, float Fc, float CurveSpread) {
	float k = (500.0f + 0.01f*Fc) * CurveSpread;
	if(k <  25.0f) k =  25.0f;
	if(k > 750.0f) k = 750.0f;

	//! Rescale to bands and return oscillator parameter
	k *= BandsPerHz;
	float s, c; Fourier_SinCos(1.0f / k, &s, &c);
	return 2.0f*c;
}
void Block_Transform_ComputeFinalPower(const float *Coef, float *FinalPower, size_t BlockSize, float Nyquist_Hz) {
	size_t i, j;
	size_t N = BlockSize/2;

	//! Convolve with spreading function
	//! The point of this particular psychoacoustic model
	//! is that neighbouring bands 'support' a given band
	//! NOTE: There are twice as many bands in Coef than we store
	for(i=0;i<N;i++) {
		float Fc = (i+0.5f) * Nyquist_Hz / N;
		const float *CoefSrc = Coef + i*2;

		//! Convolve LHS, RHS
		//! Based on a squared cosine curve
		float BandPow = SQR(CoefSrc[0]) + SQR(CoefSrc[1]);
		float PowSum = 0.125f * SQR(BandPow); {
			float Curve, CurveOld;
			float CurveOmg = Block_Transform_ComputeFinalPower_SupportCurveParam(N / Nyquist_Hz, Fc, 1.0f);

			//! Using linear prediction for the cosine curve
			CurveOld = 1.0f, Curve = CurveOmg*0.5f;
			for(j=1;j<=i*2;j++) {
				PowSum += BandPow*ABS(CoefSrc[-j]) * SQR(Curve);

				float t = Curve;
				Curve = Curve*CurveOmg - CurveOld, CurveOld = t;
				if(Curve < 0.0f) break;
			}
			CurveOld = 1.0f, Curve = CurveOmg*0.5f;
			for(j=2;j<(N-i)*2;j++) {
				PowSum += BandPow*ABS(CoefSrc[+j]) * SQR(Curve);

				float t = Curve;
				Curve = Curve*CurveOmg - CurveOld, CurveOld = t;
				if(Curve < 0.0f) break;
			}
		}
		FinalPower[i] = PowSum;
	}
}

/**************************************/

//! Insert keys for block coefficients
//! Returns updated number of keys in list
//! NOTE:
//!  AnalysisPower is used to alter the preference for the currently-being-analyzed channel
static size_t Block_Transform_InsertKeys(const float *Coef, size_t BlockSize, size_t Chan, struct AnalysisKey_t *Keys, size_t nKeys, float AnalysisPower, const float *FinalPower) {
	size_t i;
	for(i=0;i<BlockSize;i++) {
		//! Check that value doesn't collapse to 0 under the smallest quantizer (1.0)
		float v2 = SQR(Coef[i]);
		if(v2 < SQR(0.5f)) continue;

		//! Build and insert key
		Keys[nKeys].Band = i;
		Keys[nKeys].Chan = Chan;
		Keys[nKeys].Val  = v2 * FinalPower[i/2] * AnalysisPower;
		nKeys++;
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
		//!  TempBuffer has 2*BlockSize elements
		//!  The first half will contain our (scaled) sample data
		//!  The second half will be for MDCT processing
		float *BufferSample    = State->TransformTemp;
		float *BufferTemp      = State->TransformTemp + BlockSize;
		float *BufferTransform = State->TransformBuffer[Chan];
		float *BufferFwdLap    = State->TransformFwdLap[Chan];

		//! Fetch sample data
		//! Pre-scale for scaled IMDCT(*2.0/BlockSize) and SumDif transform(*0.5)
		Block_Transform_CopySamples(BufferSample, Data + Chan*BlockSize, BlockSize, 2.0f/BlockSize * 0.5f);

		//! Apply transforms
		//! NOTE: Final power stored to BufferTemp
		Fourier_MDCT(BufferTransform, BufferSample, BufferFwdLap, BufferTemp, BlockSize, State->BlockOverlap);
		Block_Transform_ComputeFinalPower(BufferTransform, BufferTemp, BlockSize, State->RateHz*0.5f);
		ULC_Transform_AntiPreEcho(BufferTransform, BlockSize);

		//! Insert coefficient keys
		nKeys = Block_Transform_InsertKeys(BufferTransform, BlockSize, Chan, State->AnalysisKeys, nKeys, AnalysisPower, BufferTemp);
		AnalysisPower *= PowerDecay;
	}
	Analysis_KeysValSort(State->AnalysisKeys, nKeys);
	return nKeys;
}

/**************************************/
//! EOF
/**************************************/
