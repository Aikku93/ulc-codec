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

//! Copy scaled samples to buffer
static void Block_Transform_CopySamples(float *DataDst, const float *DataSrc, size_t N, float Scale) {
	size_t i;
#if defined(__AVX__)
	__m256 mScale = _mm256_set1_ps(Scale);
	for(i=0;i<N;i+=16) {
		__m256 v0 = _mm256_load_ps(DataSrc + i+0);
		__m256 v1 = _mm256_load_ps(DataSrc + i+8);
		v0 = _mm256_mul_ps(v0, mScale);
		v1 = _mm256_mul_ps(v1, mScale);
		_mm256_store_ps(DataDst + i+0, v0);
		_mm256_store_ps(DataDst + i+8, v1);
	}
#elif defined(__SSE__)
	__m128 mScale = _mm_set1_ps(Scale);
	for(i=0;i<N;i+=16) {
		__m128 v0 = _mm_load_ps(DataSrc + i+ 0);
		__m128 v1 = _mm_load_ps(DataSrc + i+ 4);
		__m128 v2 = _mm_load_ps(DataSrc + i+ 8);
		__m128 v3 = _mm_load_ps(DataSrc + i+12);
		v0 = _mm_mul_ps(v0, mScale);
		v1 = _mm_mul_ps(v1, mScale);
		v2 = _mm_mul_ps(v2, mScale);
		v3 = _mm_mul_ps(v3, mScale);
		_mm_store_ps(DataDst + i+ 0, v0);
		_mm_store_ps(DataDst + i+ 4, v1);
		_mm_store_ps(DataDst + i+ 8, v2);
		_mm_store_ps(DataDst + i+12, v3);
	}
#else
	for(i=0;i<N;i++) DataDst[i] = DataSrc[i] * Scale;
#endif
}

/**************************************/

//! Insert keys for block coefficients
//! Returns updated number of keys in list
//! NOTE:
//!  AnalysisPower is used to alter the preference for
//!  the currently-being-analyzed channel
static size_t Block_Transform_InsertKeys(const float *Coef, size_t BlockSize, size_t Chan,  struct AnalysisKey_t *Keys, size_t nKeys, float AnalysisPower) {
	size_t i;
	struct AnalysisKey_t Key = {.Chan = Chan};
	for(i=0;i<BlockSize;i++) {
		//! Check that value doesn't collapse to 0 under the smallest quantizer (1.0)
		float v = Coef[i]*Coef[i];
		if(v < 0.5f*0.5f) continue;

		//! To avoid excessive muffling, we slightly decrease
		//! the importance of the low frequency bands using
		//! a pseudo-logarithmic curve (inverted parabola)
		float lfScale = (BlockSize-i - 0.5f)/BlockSize;
		      lfScale = 1.0f - 0.5f*lfScale*lfScale;

		//! Build and insert key
		Key.Band = i;
		//Key.Chan = Chan;
		Key.Val  = v * AnalysisPower * lfScale;
		Analysis_KeyInsert(&Key, Keys, nKeys);
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
static size_t Block_Transform(const struct ULC_EncoderState_t *State, const float *Data) {
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
		Fourier_MDCT(BufferTransform, BufferSample, BufferFwdLap, BufferTemp, BlockSize);
		ULC_Transform_AntiPreEcho(BufferTransform, BlockSize);

		//! Insert coefficient keys
		//! NOTE:
		//!  Channels get increasingly lower analysis power,
		//!  as the input is expected to have some form of
		//!  decorrelation transform applied, such as M/S,
		//!  Hadamard, DCT, etc.
		nKeys = Block_Transform_InsertKeys(BufferTransform, BlockSize, Chan, State->AnalysisKeys, nKeys, AnalysisPower);
		AnalysisPower *= 0x1.6A09E6p-1; //! 1/sqrt[2]
	}
	return nKeys;
}

/**************************************/
//! EOF
/**************************************/