/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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
#include <stddef.h>
#include <stdint.h>
/**************************************/
#include "ulcUtility.h"
/**************************************/
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_Helper.h"
/**************************************/
#include "ulcEncoder_Psycho.h"
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

//! Insert keys for block coefficients
//! Returns updated number of keys in list
//! NOTE:
//!  AnalysisPower is used to alter the preference for the currently-being-analyzed channel
static size_t Block_Transform_InsertKeys(const float *Coef, size_t BlockSize, size_t Chan, struct AnalysisKey_t *Keys, size_t nKeys, float AnalysisPower, const float *MaskingPower) {
	size_t i;
#if !USE_PSYCHOACOUSTICS
	(void)MaskingPower;
#endif
	//! Start inserting keys
	for(i=0;i<BlockSize;i++) {
		//! Check that value doesn't collapse to 0 under the smallest quantizer (1.0)
		float v2 = SQR(Coef[i]);
		if(v2 >= SQR(0.5f)) {
			//! Build and insert key
			Keys[nKeys].Band = i;
			Keys[nKeys].Chan = Chan;
#if ULC_USE_PSYHOACOUSTICS
			Keys[nKeys].Val  = v2*AnalysisPower - MaskingPower[i/2];
#else
			Keys[nKeys].Val  = v2*AnalysisPower;
#endif
			nKeys++;
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
		float *BufferFlatness  = State->TransformFlatness[Chan];
		float *BufferTransform = State->TransformBuffer[Chan];
		float *BufferFwdLap    = State->TransformFwdLap[Chan];

		//! Fetch sample data
		//! Pre-scale for scaled IMDCT(*2.0/BlockSize) and SumDif transform(*0.5)
		Block_Transform_CopySamples(BufferSample, Data + Chan*BlockSize, BlockSize, 2.0f/BlockSize * 0.5f);

		//! Apply transforms
		//! NOTE: Masking power stored to BufferTemp
		Fourier_MDCT(BufferTransform, BufferSample, BufferFwdLap, BufferTemp, BlockSize, State->BlockOverlap);
		Block_Transform_ComputeFlatness(BufferTransform, BufferFlatness, BlockSize);
#if ULC_USE_PSYHOACOUSTICS
		Block_Transform_ComputeMaskingPower(BufferTransform, BufferMasking, BufferFlatness, BlockSize, State->RateHz * 0.5f);
#endif
		ULC_Transform_AntiPreEcho(BufferTransform, BlockSize);

		//! Insert coefficient keys
		nKeys = Block_Transform_InsertKeys(BufferTransform, BlockSize, Chan, State->AnalysisKeys, nKeys, AnalysisPower, BufferMasking);
		AnalysisPower *= PowerDecay;
	}
	Analysis_KeysValSort(State->AnalysisKeys, nKeys);
	return nKeys;
}

/**************************************/
//! EOF
/**************************************/
