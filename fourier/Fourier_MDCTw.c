/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#if defined(__AVX__) || defined(__FMA__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include <stddef.h>
/**************************************/
#include "Fourier.h"
/**************************************/

void Fourier_MDCTw(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, size_t N, const float *Window) {
	size_t i;

	//! Copy state from old block
	for(i=0;i<N/2;i++) BufOut[N/2+i] = BufLap[i];

	//! Perform windowed lapping
#if defined(__AVX__)
	__m256 a, b;
	__m256 s, c;
	__m256 t0, t1;

	//! First quarter
	for(i=0;i<N/4;i+=8) {
		b = _mm256_load_ps(BufIn - i + N-8);
		a = _mm256_load_ps(BufIn + i);
		s = _mm256_load_ps(Window + i);
		c = _mm256_load_ps(Window - i + N-8);
		b = _mm256_shuffle_ps(b, b, 0x1B);
		b = _mm256_permute2f128_ps(b, b, 0x01);
		c = _mm256_shuffle_ps(c, c, 0x1B);
		c = _mm256_permute2f128_ps(c, c, 0x01);
#if defined(__FMA__)
		t0 = _mm256_mul_ps(s, b);
		t1 = _mm256_mul_ps(c, b);
		t0 = _mm256_fmadd_ps (c, a, t0);
		t1 = _mm256_fnmadd_ps(s, a, t1);
#else
		t0 = _mm256_add_ps(_mm256_mul_ps(s, b), _mm256_mul_ps(c, a));
		t1 = _mm256_sub_ps(_mm256_mul_ps(c, b), _mm256_mul_ps(s, a));
#endif
		t0 = _mm256_shuffle_ps(t0, t0, 0x1B);
		t0 = _mm256_permute2f128_ps(t0, t0, 0x01);
		_mm256_store_ps(BufOut - i + N/2-8, t0);
		_mm256_store_ps(BufLap + i,         t1);
	}

	//! Second quarter
	for(i=0;i<N/4;i+=8) {
		b = _mm256_load_ps(BufIn - i + N/2-8);
		a = _mm256_load_ps(BufIn + i + N/2);
		s = _mm256_load_ps(Window + i + N/2);
		c = _mm256_load_ps(Window - i + N/2-8);
		b = _mm256_shuffle_ps(b, b, 0x1B);
		b = _mm256_permute2f128_ps(b, b, 0x01);
		c = _mm256_shuffle_ps(c, c, 0x1B);
		c = _mm256_permute2f128_ps(c, c, 0x01);
#if defined(__FMA__)
		t0 = _mm256_mul_ps(c, a);
		t1 = _mm256_mul_ps(s, a);
		t0 = _mm256_fmadd_ps (s, b, t0);
		t1 = _mm256_fnmadd_ps(c, b, t1);
#else
		t0 = _mm256_add_ps(_mm256_mul_ps(c, a), _mm256_mul_ps(s, b));
		t1 = _mm256_sub_ps(_mm256_mul_ps(s, a), _mm256_mul_ps(c, b));
#endif
		t1 = _mm256_shuffle_ps(t1, t1, 0x1B);
		t1 = _mm256_permute2f128_ps(t1, t1, 0x01);
		_mm256_store_ps(BufOut + i,         t0);
		_mm256_store_ps(BufLap - i + N/2-8, t1);
	}
#elif defined(__SSE__)
	__m128 a, b;
	__m128 s, c;
	__m128 t0, t1;

	//! First quarter
	for(i=0;i<N/4;i+=4) {
		b = _mm_loadr_ps(BufIn - i + N-4);
		a = _mm_load_ps (BufIn + i);
		s = _mm_load_ps (Window + i);
		c = _mm_loadr_ps(Window - i + N-4);
#if defined(__FMA__)
		t0 = _mm_mul_ps(s, b);
		t1 = _mm_mul_ps(c, b);
		t0 = _mm_fmadd_ps (c, a, t0);
		t1 = _mm_fnmadd_ps(s, a, t1);
#else
		t0 = _mm_add_ps(_mm_mul_ps(s, b), _mm_mul_ps(c, a));
		t1 = _mm_sub_ps(_mm_mul_ps(c, b), _mm_mul_ps(s, a));
#endif
		_mm_storer_ps(BufOut - i + N/2-4, t0);
		_mm_store_ps (BufLap + i,         t1);
	}

	//! Second quarter
	for(i=0;i<N/4;i+=4) {
		b = _mm_loadr_ps(BufIn - i + N/2-4);
		a = _mm_load_ps (BufIn + i + N/2);
		s = _mm_load_ps (Window + i + N/2);
		c = _mm_loadr_ps(Window - i + N/2-4);
#if defined(__FMA__)
		t0 = _mm_mul_ps(c, a);
		t1 = _mm_mul_ps(s, a);
		t0 = _mm_fmadd_ps (s, b, t0);
		t1 = _mm_fnmadd_ps(c, b, t1);
#else
		t0 = _mm_add_ps(_mm_mul_ps(c, a), _mm_mul_ps(s, b));
		t1 = _mm_sub_ps(_mm_mul_ps(s, a), _mm_mul_ps(c, b));
#endif
		_mm_store_ps (BufOut + i,         t0);
		_mm_storer_ps(BufLap - i + N/2-4, t1);
	}
#else
	//! First quarter
	for(i=0;i<N/4;i++) {
		BufOut[N/2-1-i] =  Window[N-1-i]*BufIn[i] + Window[    i]*BufIn[N-1-i];
		BufLap[i]       = -Window[    i]*BufIn[i] + Window[N-1-i]*BufIn[N-1-i];
	}

	//! Second quarter
	for(i=0;i<N/4;i++) {
		BufOut[i]       =  Window[N/2  +i]*BufIn[N/2-1-i] + Window[N/2-1-i]*BufIn[N/2+i];
		BufLap[N/2-1-i] = -Window[N/2-1-i]*BufIn[N/2-1-i] + Window[N/2  +i]*BufIn[N/2+i];
	}
#endif
	//! Do actual transform
	Fourier_DCT4(BufOut, BufTmp, N);
}

/**************************************/
//! EOF
/**************************************/
