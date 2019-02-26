/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#if defined(__AVX__)
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

void Fourier_MDCT(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, size_t N) {
	size_t i;
	float Ns = 1.0f / N;

	//! Copy state from old block
	for(i=0;i<N/2;i++) BufOut[N/2+i] = BufLap[i];

	//! Perform windowed lapping
#if defined(__AVX__)
	__m256 c, s;
	__m256 wc, ws;

	//! Get oscillator parameter
	{
		float c, s;
		Fourier_SinCos(8.0f * Ns, &s, &c);
		wc = _mm256_set1_ps(c);
		ws = _mm256_set1_ps(s);
	}

	//! First quarter
	c = _mm256_setr_ps(0.5f, 1.5f, 2.5f, 3.5f, 4.5f, 5.5f, 6.5f, 7.5f);
	c = _mm256_mul_ps(c, _mm256_set1_ps(Ns));
	Fourier_SinCosAVX(c, &s, &c);
	c = _mm256_xor_ps(c, _mm256_set1_ps(-0.0f));
	s = _mm256_xor_ps(s, _mm256_set1_ps(-0.0f));
	for(i=0;i<N/4;i+=8) {
		__m256 b = _mm256_load_ps(BufIn - i + N-8);
		__m256 a = _mm256_load_ps(BufIn + i);
		b = _mm256_shuffle_ps(b, b, 0x1B);
		b = _mm256_permute2f128_ps(b, b, 0x01);
		__m256 t0 = _mm256_add_ps(_mm256_mul_ps(b, s), _mm256_mul_ps(a, c));
		__m256 t1 = _mm256_sub_ps(_mm256_mul_ps(b, c), _mm256_mul_ps(a, s));
		t0 = _mm256_shuffle_ps(t0, t0, 0x1B);
		t0 = _mm256_permute2f128_ps(t0, t0, 0x01);
		_mm256_store_ps(BufOut - i + N/2-8, t0);
		_mm256_store_ps(BufLap + i,         t1);
		t0 = _mm256_sub_ps(_mm256_mul_ps(wc, c), _mm256_mul_ps(ws, s));
		t1 = _mm256_add_ps(_mm256_mul_ps(ws, c), _mm256_mul_ps(wc, s));
		c = t0;
		s = t1;
	}

	//! Second quarter
	c = _mm256_setr_ps(0.5f, 1.5f, 2.5f, 3.5f, 4.5f, 5.5f, 6.5f, 7.5f);
	c = _mm256_mul_ps(c, _mm256_set1_ps(Ns));
	c = _mm256_add_ps(c, _mm256_set1_ps(0.5f));
	Fourier_SinCosAVX(c, &s, &c);
	c = _mm256_xor_ps(c, _mm256_set1_ps(-0.0f));
	s = _mm256_xor_ps(s, _mm256_set1_ps(-0.0f));
	for(i=0;i<N/4;i+=8) {
		__m256 b = _mm256_load_ps(BufIn - i + N/2-8);
		__m256 a = _mm256_load_ps(BufIn + i + N/2);
		b = _mm256_shuffle_ps(b, b, 0x1B);
		b = _mm256_permute2f128_ps(b, b, 0x01);
		__m256 t0 = _mm256_add_ps(_mm256_mul_ps(a, c), _mm256_mul_ps(b, s));
		__m256 t1 = _mm256_sub_ps(_mm256_mul_ps(a, s), _mm256_mul_ps(b, c));
		t1 = _mm256_shuffle_ps(t1, t1, 0x1B);
		t1 = _mm256_permute2f128_ps(t1, t1, 0x01);
		_mm256_store_ps(BufOut + i,         t0);
		_mm256_store_ps(BufLap - i + N/2-8, t1);
		t0 = _mm256_sub_ps(_mm256_mul_ps(wc, c), _mm256_mul_ps(ws, s));
		t1 = _mm256_add_ps(_mm256_mul_ps(ws, c), _mm256_mul_ps(wc, s));
		c = t0;
		s = t1;
	}
#elif defined(__SSE__)
	__m128 c, s;
	__m128 wc, ws;

	//! Get oscillator parameter
	{
		float c, s;
		Fourier_SinCos(4.0f * Ns, &s, &c);
		wc = _mm_set1_ps(c);
		ws = _mm_set1_ps(s);
	}

	//! First quarter
	c = _mm_setr_ps(0.5f, 1.5f, 2.5f, 3.5f);
	c = _mm_mul_ps(c, _mm_set1_ps(Ns));
	Fourier_SinCosSSE(c, &s, &c);
	c = _mm_xor_ps(c, _mm_set1_ps(-0.0f));
	s = _mm_xor_ps(s, _mm_set1_ps(-0.0f));
	for(i=0;i<N/4;i+=4) {
		__m128 b = _mm_loadr_ps(BufIn - i + N-4);
		__m128 a = _mm_load_ps (BufIn + i);
		__m128 t0 = _mm_add_ps(_mm_mul_ps(b, s), _mm_mul_ps(a, c));
		__m128 t1 = _mm_sub_ps(_mm_mul_ps(b, c), _mm_mul_ps(a, s));
		_mm_storer_ps(BufOut - i + N/2-4, t0);
		_mm_store_ps (BufLap + i,         t1);
		t0 = _mm_sub_ps(_mm_mul_ps(wc, c), _mm_mul_ps(ws, s));
		t1 = _mm_add_ps(_mm_mul_ps(ws, c), _mm_mul_ps(wc, s));
		c = t0;
		s = t1;
	}

	//! Second quarter
	c = _mm_setr_ps(0.5f, 1.5f, 2.5f, 3.5f);
	c = _mm_mul_ps(c, _mm_set1_ps(Ns));
	c = _mm_add_ps(c, _mm_set1_ps(0.5f));
	Fourier_SinCosSSE(c, &s, &c);
	c = _mm_xor_ps(c, _mm_set1_ps(-0.0f));
	s = _mm_xor_ps(s, _mm_set1_ps(-0.0f));
	for(i=0;i<N/4;i+=4) {
		__m128 b = _mm_loadr_ps(BufIn - i + N/2-4);
		__m128 a = _mm_load_ps (BufIn + i + N/2);
		__m128 t0 = _mm_add_ps(_mm_mul_ps(a, c), _mm_mul_ps(b, s));
		__m128 t1 = _mm_sub_ps(_mm_mul_ps(a, s), _mm_mul_ps(b, c));
		_mm_store_ps (BufOut + i,         t0);
		_mm_storer_ps(BufLap - i + N/2-4, t1);
		t0 = _mm_sub_ps(_mm_mul_ps(wc, c), _mm_mul_ps(ws, s));
		t1 = _mm_add_ps(_mm_mul_ps(ws, c), _mm_mul_ps(wc, s));
		c = t0;
		s = t1;
	}
#else
	float c, s;
	float wc, ws;
	Fourier_SinCos(Ns, &wc, &ws);

	//! First quarter
	Fourier_SinCos(0.5f*Ns, &s, &c);
	c = -c;
	s = -s;
	for(i=0;i<N/4;i++) {
		BufOut[N/2-1-i] =  BufIn[i]*c + BufIn[N-1-i]*s;
		BufLap[i]       = -BufIn[i]*s + BufIn[N-1-i]*c;
		float _c = wc*c - ws*s;
		float _s = ws*c + wc*s;
		c = _c, s = _s;
	}

	//! Second quarter
	Fourier_SinCos(0.5f*Ns + 0.5f, &s, &c);
	c = -c;
	s = -s;
	for(i=0;i<N/4;i++) {
		BufOut[i]       =  BufIn[N/2-1-i]*s + BufIn[N/2+i]*c;
		BufLap[N/2-1-i] = -BufIn[N/2-1-i]*c + BufIn[N/2+i]*s;
		float _c = wc*c - ws*s;
		float _s = ws*c + wc*s;
		c = _c, s = _s;
	}
#endif
	//! Do actual transform
	Fourier_DCT4(BufOut, BufTmp, N);
}

/**************************************/
//! EOF
/**************************************/
