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
static void Block_Transform_InsertKeys(
	struct AnalysisKey_t *Keys,
	const float *Coef,
	const float *CoefNp,
	size_t BlockSize,
	size_t *nKeys,
	size_t Chan,
	float AnalysisPower,
	float NyquistHz
) {
	size_t Band;
#if MASKING_BAND_DECIMATION_FACTOR > 1
	float MaskCur = 0.0f; //! <- Initializing shuts gcc up
	float MaskNxt = Block_Transform_ComputeMaskingPower(Coef, CoefNp, 0, BlockSize, NyquistHz);
#endif
	for(Band=0;Band<BlockSize;Band++) {
#if MASKING_BAND_DECIMATION_FACTOR > 1
		//! Update masking
		size_t MaskMu = Band%MASKING_BAND_DECIMATION_FACTOR;
		if(MaskMu == 0) {
			MaskCur = MaskNxt;
			size_t NextBand = Band + MASKING_BAND_DECIMATION_FACTOR;
			if(NextBand < BlockSize) MaskNxt = Block_Transform_ComputeMaskingPower(Coef, CoefNp, NextBand, BlockSize, NyquistHz);
		}
		float Mask = (MaskCur*(MASKING_BAND_DECIMATION_FACTOR - MaskMu) + MaskNxt*MaskMu) * (1.0f / MASKING_BAND_DECIMATION_FACTOR);
#endif
		//! Check that the value is in range of the smallest quantization
		float ValNp = SQR(Coef[Band]);
		if(ValNp >= 0.5f*SQR(ULC_COEF_EPS)) {
			Keys[*nKeys].Band = Band;
			Keys[*nKeys].Chan = Chan;
#if MASKING_BAND_DECIMATION_FACTOR > 1
			Keys[*nKeys].Val  = ValNp*AnalysisPower - Mask;
#else
			Keys[*nKeys].Val  = ValNp*AnalysisPower - Block_Transform_ComputeMaskingPower(Coef, CoefNp, Band, BlockSize, NyquistHz);
#endif
			(*nKeys)++;
		}
	}
}

/**************************************/

//! Apply block transform
//!  -Fetches data
//!  -Applies MDCT
//!  -Stores keys for block coefficients
//! Returns the number of keys stored
static inline void Block_Transform_ToNepers(float *Dst, float *Src, size_t BlockSize) {
	size_t Band;
	for(Band=0;Band<BlockSize;Band++) {
		float v = ABS(Src[Band]);
		      v = (v < 0.5f*ULC_COEF_EPS) ? ULC_COEF_NEPER_OUT_OF_RANGE : logf(v);
		Dst[Band] = v - ULC_COEF_EPS*(v == 0.0f);
	}
}
static size_t Block_Transform(const struct ULC_EncoderState_t *State, const float *Data, float PowerDecay) {
	size_t nChan     = State->nChan;
	size_t BlockSize = State->BlockSize;

	size_t Chan;
	size_t nKeys = 0;
	float  AnalysisPower = 1.0f;
	for(Chan=0;Chan<nChan;Chan++) {
		//! Get buffer pointers
		float *BufferTransform = State->TransformBuffer[Chan];
		float *BufferNepers    = State->TransformNepers[Chan];
		float *BufferFwdLap    = State->TransformFwdLap[Chan];
		float *BufferTemp      = State->TransformTemp;

		//! Fetch sample data, pre-scaled for scaled IMDCT(*2.0/BlockSize)
		//! NOTE: Sample data temporarily stored to BufferNepers
		Block_Transform_CopySamples(BufferNepers, Data + Chan*BlockSize, BlockSize, 2.0f/BlockSize);

		//! Apply transforms and insert keys
		Fourier_MDCT(BufferTransform, BufferNepers, BufferFwdLap, BufferTemp, BlockSize, State->BlockOverlap);
		Block_Transform_ToNepers(BufferNepers, BufferTransform, BlockSize);
		Block_Transform_InsertKeys(State->AnalysisKeys, BufferTransform, BufferNepers, BlockSize, &nKeys, Chan, AnalysisPower, State->RateHz * 0.5f);
		AnalysisPower *= PowerDecay;
	}
	Analysis_KeysValSort(State->AnalysisKeys, nKeys);
	return nKeys;
}

/**************************************/
//! EOF
/**************************************/
