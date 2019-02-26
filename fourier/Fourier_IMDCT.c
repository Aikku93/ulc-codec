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
/**************************************/
#include "Fourier.h"
/**************************************/

void Fourier_IMDCT(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, size_t N) {
	size_t i;
	float Ns = M_PI_2 / N;

	//! Undo transform
	for(i=0;i<N;i++) BufTmp[i] = BufIn[i];
	Fourier_DCT4(BufTmp, BufOut, N);

	//! Undo lapping
#if defined(__AVX__)
	__m256 c, s;
	__m256 wc = _mm256_set1_ps(cos(8.0f*Ns));
	__m256 ws = _mm256_set1_ps(sin(8.0f*Ns));

	//! First quarter
	c = _mm256_setr_ps(
		cos((0 + 0.5f)*Ns),
		cos((1 + 0.5f)*Ns),
		cos((2 + 0.5f)*Ns),
		cos((3 + 0.5f)*Ns),
		cos((4 + 0.5f)*Ns),
		cos((5 + 0.5f)*Ns),
		cos((6 + 0.5f)*Ns),
		cos((7 + 0.5f)*Ns)
	);
	s = _mm256_setr_ps(
		sin((0 + 0.5f)*Ns),
		sin((1 + 0.5f)*Ns),
		sin((2 + 0.5f)*Ns),
		sin((3 + 0.5f)*Ns),
		sin((4 + 0.5f)*Ns),
		sin((5 + 0.5f)*Ns),
		sin((6 + 0.5f)*Ns),
		sin((7 + 0.5f)*Ns)
	);
	c = _mm256_xor_ps(c, _mm256_set1_ps(-0.0f));
	s = _mm256_xor_ps(s, _mm256_set1_ps(-0.0f));
	for(i=0;i<N/4;i+=8) {
		__m256 b = _mm256_load_ps(BufLap - i + N/2-8);
		__m256 a = _mm256_load_ps(BufTmp + i + N/2);
		b = _mm256_shuffle_ps(b, b, 0x1B);
		b = _mm256_permute2f128_ps(b, b, 0x01);
		__m256 t0 = _mm256_sub_ps(_mm256_mul_ps(b, c), _mm256_mul_ps(a, s));
		__m256 t1 = _mm256_add_ps(_mm256_mul_ps(b, s), _mm256_mul_ps(a, c));
		t1 = _mm256_shuffle_ps(t1, t1, 0x1B);
		t1 = _mm256_permute2f128_ps(t1, t1, 0x01);
		_mm256_store_ps(BufOut + i,       t0);
		_mm256_store_ps(BufOut - i + N-8, t1);
		t0 = _mm256_sub_ps(_mm256_mul_ps(wc, c), _mm256_mul_ps(ws, s));
		t1 = _mm256_add_ps(_mm256_mul_ps(ws, c), _mm256_mul_ps(wc, s));
		c = t0;
		s = t1;
	}

	//! Second quarter
	c = _mm256_setr_ps(
		cos((N/2+0 + 0.5f)*Ns),
		cos((N/2+1 + 0.5f)*Ns),
		cos((N/2+2 + 0.5f)*Ns),
		cos((N/2+3 + 0.5f)*Ns),
		cos((N/2+4 + 0.5f)*Ns),
		cos((N/2+5 + 0.5f)*Ns),
		cos((N/2+6 + 0.5f)*Ns),
		cos((N/2+7 + 0.5f)*Ns)
	);
	s = _mm256_setr_ps(
		sin((N/2+0 + 0.5f)*Ns),
		sin((N/2+1 + 0.5f)*Ns),
		sin((N/2+2 + 0.5f)*Ns),
		sin((N/2+3 + 0.5f)*Ns),
		sin((N/2+4 + 0.5f)*Ns),
		sin((N/2+5 + 0.5f)*Ns),
		sin((N/2+6 + 0.5f)*Ns),
		sin((N/2+7 + 0.5f)*Ns)
	);
	c = _mm256_xor_ps(c, _mm256_set1_ps(-0.0f));
	s = _mm256_xor_ps(s, _mm256_set1_ps(-0.0f));
	for(i=0;i<N/4;i+=8) {
		__m256 b = _mm256_load_ps(BufTmp - i + N-8);
		__m256 a = _mm256_load_ps(BufLap + i);
		b = _mm256_shuffle_ps(b, b, 0x1B);
		b = _mm256_permute2f128_ps(b, b, 0x01);
		__m256 t0 = _mm256_sub_ps(_mm256_mul_ps(a, s), _mm256_mul_ps(b, c));
		__m256 t1 = _mm256_add_ps(_mm256_mul_ps(a, c), _mm256_mul_ps(b, s));
		t0 = _mm256_shuffle_ps(t0, t0, 0x1B);
		t0 = _mm256_permute2f128_ps(t0, t0, 0x01);
		_mm256_store_ps(BufOut - i + N/2-8, t0);
		_mm256_store_ps(BufOut + i + N/2,   t1);
		t0 = _mm256_sub_ps(_mm256_mul_ps(wc, c), _mm256_mul_ps(ws, s));
		t1 = _mm256_add_ps(_mm256_mul_ps(ws, c), _mm256_mul_ps(wc, s));
		c = t0;
		s = t1;
	}
#elif defined(__SSE__)
	__m128 c, s;
	__m128 wc = _mm_set1_ps(cos(4.0f*Ns));
	__m128 ws = _mm_set1_ps(sin(4.0f*Ns));

	//! First quarter
	c = _mm_setr_ps(
		cos((0 + 0.5f)*Ns),
		cos((1 + 0.5f)*Ns),
		cos((2 + 0.5f)*Ns),
		cos((3 + 0.5f)*Ns)
	);
	s = _mm_setr_ps(
		sin((0 + 0.5f)*Ns),
		sin((1 + 0.5f)*Ns),
		sin((2 + 0.5f)*Ns),
		sin((3 + 0.5f)*Ns)
	);
	c = _mm_xor_ps(c, _mm_set1_ps(-0.0f));
	s = _mm_xor_ps(s, _mm_set1_ps(-0.0f));
	for(i=0;i<N/4;i+=4) {
		__m128 b = _mm_loadr_ps(BufLap - i + N/2-4);
		__m128 a = _mm_load_ps (BufTmp + i + N/2);
		__m128 t0 = _mm_sub_ps(_mm_mul_ps(b, c), _mm_mul_ps(a, s));
		__m128 t1 = _mm_add_ps(_mm_mul_ps(b, s), _mm_mul_ps(a, c));
		_mm_store_ps (BufOut + i,       t0);
		_mm_storer_ps(BufOut - i + N-4, t1);
		t0 = _mm_sub_ps(_mm_mul_ps(wc, c), _mm_mul_ps(ws, s));
		t1 = _mm_add_ps(_mm_mul_ps(ws, c), _mm_mul_ps(wc, s));
		c = t0;
		s = t1;
	}

	//! Second quarter
	c = _mm_setr_ps(
		cos((N/2+0 + 0.5f)*Ns),
		cos((N/2+1 + 0.5f)*Ns),
		cos((N/2+2 + 0.5f)*Ns),
		cos((N/2+3 + 0.5f)*Ns)
	);
	s = _mm_setr_ps(
		sin((N/2+0 + 0.5f)*Ns),
		sin((N/2+1 + 0.5f)*Ns),
		sin((N/2+2 + 0.5f)*Ns),
		sin((N/2+3 + 0.5f)*Ns)
	);
	c = _mm_xor_ps(c, _mm_set1_ps(-0.0f));
	s = _mm_xor_ps(s, _mm_set1_ps(-0.0f));
	for(i=0;i<N/4;i+=4) {
		__m128 b = _mm_loadr_ps(BufTmp - i + N-4);
		__m128 a = _mm_load_ps (BufLap + i);
		__m128 t0 = _mm_sub_ps(_mm_mul_ps(a, s), _mm_mul_ps(b, c));
		__m128 t1 = _mm_add_ps(_mm_mul_ps(a, c), _mm_mul_ps(b, s));
		_mm_storer_ps(BufOut - i + N/2-4, t0);
		_mm_store_ps (BufOut + i + N/2,   t1);
		t0 = _mm_sub_ps(_mm_mul_ps(wc, c), _mm_mul_ps(ws, s));
		t1 = _mm_add_ps(_mm_mul_ps(ws, c), _mm_mul_ps(wc, s));
		c = t0;
		s = t1;
	}
#else
	float c, s, wc = cos(1.0f*Ns), ws = sin(1.0f*Ns);

	//! First quarter
	c = -cos((0 + 0.5f)*Ns), s = -sin((0 + 0.5f)*Ns);
	for(i=0;i<N/4;i++) {
		BufOut[i]       = BufLap[N/2-1-i]*c - BufTmp[N/2+i]*s;
		BufOut[N-1-i]   = BufLap[N/2-1-i]*s + BufTmp[N/2+i]*c;
		float _c = wc*c - ws*s;
		float _s = ws*c + wc*s;
		c = _c, s = _s;
	}

	//! Second quarter
	c = -cos((N/2 + 0.5f)*Ns), s = -sin((N/2 + 0.5f)*Ns);
	for(i=0;i<N/4;i++) {
		BufOut[N/2-1-i] = BufLap[i]*s - BufTmp[N-1-i]*c;
		BufOut[N/2+i]   = BufLap[i]*c + BufTmp[N-1-i]*s;
		float _c = wc*c - ws*s;
		float _s = ws*c + wc*s;
		c = _c, s = _s;
	}
#endif
	//! Copy state to old block
	for(i=0;i<N/2;i++) BufLap[i] = BufTmp[i];
}

/**************************************/
//! EOF
/**************************************/
