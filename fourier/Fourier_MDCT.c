/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#if defined(__AVX__) || defined(__FMA__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include "Fourier.h"
/**************************************/

void Fourier_MDCT(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, int N, int Overlap, float *BufMDST) {
	int i;
	const float *WinS = Fourier_SinTableN(Overlap);
	const float *WinC = WinS + Overlap;
	const float *InLo = BufIn;
	const float *InHi = BufIn  + N;
	      float *Out  = BufOut + N/2;
	      float *Lap  = BufLap;

	//! Copy state from old block
	for(i=0;i<N/2;i++) Out[i] = Lap[i];

	//! Perform windowed lapping
#if defined(__AVX__)
	__m256 a, b;
	__m256 t0, t1;
	__m256 c, s;

	for(i=0;i<(N-Overlap)/2;i+=8) {
		InHi -= 8; b = _mm256_load_ps(InHi);
		a = _mm256_load_ps(InLo); InLo += 8;
		b = _mm256_shuffle_ps(b, b, 0x1B);
		b = _mm256_permute2f128_ps(b, b, 0x01);
		a = _mm256_shuffle_ps(a, a, 0x1B);
		a = _mm256_permute2f128_ps(a, a, 0x01);
		_mm256_store_ps(Lap, b); Lap += 8;
		Out -= 8; _mm256_store_ps(Out, a);
	}
	for(;i<N/2;i+=8) {
		InHi -= 8; b = _mm256_load_ps(InHi);
		a = _mm256_load_ps(InLo); InLo += 8;
		b = _mm256_shuffle_ps(b, b, 0x1B);
		b = _mm256_permute2f128_ps(b, b, 0x01);
		WinC -= 8; c = _mm256_load_ps(WinC);
		s = _mm256_load_ps(WinS); WinS += 8;
		c = _mm256_shuffle_ps(c, c, 0x1B);
		c = _mm256_permute2f128_ps(c, c, 0x01);
#if defined(__FMA__)
		t0 = _mm256_mul_ps(s, b);
		t1 = _mm256_mul_ps(s, a);
		t0 = _mm256_fmadd_ps(c, a, t0);
		t1 = _mm256_fmsub_ps(c, b, t1);
#else
		t0 = _mm256_add_ps(_mm256_mul_ps(c, a), _mm256_mul_ps(s, b));
		t1 = _mm256_sub_ps(_mm256_mul_ps(c, b), _mm256_mul_ps(s, a));
#endif
		t0 = _mm256_shuffle_ps(t0, t0, 0x1B);
		t0 = _mm256_permute2f128_ps(t0, t0, 0x01);
		Out -= 8; _mm256_store_ps(Out, t0);
		_mm256_store_ps(Lap, t1); Lap += 8;
	}
#elif defined(__SSE__)
	__m128 a, b;
	__m128 t0, t1;
	__m128 c, s;

	for(i=0;i<(N-Overlap)/2;i+=4) {
		InHi -= 4; b = _mm_loadr_ps(InHi);
		a = _mm_load_ps(InLo); InLo += 4;
		_mm_store_ps(Lap, b); Lap += 4;
		Out -= 4; _mm_storer_ps(Out, a);
	}
	for(;i<N/2;i+=4) {
		InHi -= 4; b = _mm_loadr_ps(InHi);
		a = _mm_load_ps(InLo); InLo += 4;
		s = _mm_load_ps(WinS); WinS += 4;
		WinC -= 4; c = _mm_loadr_ps(WinC);
#if defined(__FMA__)
		t0 = _mm_mul_ps(s, b);
		t1 = _mm_mul_ps(s, a);
		t0 = _mm_fmadd_ps(c, a, t0);
		t1 = _mm_fmsub_ps(c, b, t1);
#else
		t0 = _mm_add_ps(_mm_mul_ps(c, a), _mm_mul_ps(s, b));
		t1 = _mm_sub_ps(_mm_mul_ps(c, b), _mm_mul_ps(s, a));
#endif
		Out -= 4; _mm_storer_ps(Out, t0);
		_mm_store_ps(Lap, t1); Lap += 4;
	}
#else
	for(i=0;i<(N-Overlap)/2;i++) {
		float a = *InLo++;
		float b = *--InHi;
		*--Out = a;
		*Lap++ = b;
	}
	for(;i<N/2;i++) {
		float a = *InLo++;
		float b = *--InHi;
		float c = *--WinC;
		float s = *WinS++;
		*--Out =  c*a + s*b;
		*Lap++ = -s*a + c*b;
	}
#endif
	//! Do actual transform
	if(BufMDST) {
		//! Use aliased data to compute the DST via
		//! a DCT using trigonometric relations
		{
#if defined(__AVX__)
			__m256 v, NegMask = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f);
			for(i=0;i<N;i+=8) {
				v = _mm256_load_ps(BufOut + i);
				v = _mm256_xor_ps(v, NegMask);
				_mm256_store_ps(BufMDST + i, v);
			}
#elif defined(__SSE__)
			__m128 v, NegMask = _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f);
			for(i=0;i<N;i+=4) {
				v = _mm_load_ps(BufOut + i);
				v = _mm_xor_ps(v, NegMask);
				_mm_store_ps(BufMDST + i, v);
			}
#else
			for(i=0;i<N;i+=2) {
				BufMDST[i  ] =  BufOut[i  ];
				BufMDST[i+1] = -BufOut[i+1];
			}
#endif
		}
		Fourier_DCT4(BufMDST, BufTmp, N);
		{
			float *BufLo = BufMDST;
			float *BufHi = BufMDST + N;
#if defined(__AVX_)
			__m256 v0, v1;
			for(i=0;i<N/2;i+=8) {
				BufHi -= 8;
				v0 = _mm256_load_ps(BufLo);
				v1 = _mm256_load_ps(BufHi);
				v0 = _mm256_shuffle_ps(v0, v0, 0x1B);
				v1 = _mm256_shuffle_ps(v1, v1, 0x1B);
				v0 = _mm256_permute2f128_ps(v0, v0, 0x01);
				v1 = _mm256_permute2f128_ps(v1, v1, 0x01);
				_mm256_store_ps(BufHi, v0);
				_mm256_store_ps(BufLo, v1);
				BufLo += 8;
			}
#elif defined(__SSE__)
			__m128 v0, v1;
			for(i=0;i<N/2;i+=4) {
				BufHi -= 4;
				v0 = _mm_load_ps(BufLo);
				v1 = _mm_load_ps(BufHi);
				_mm_storer_ps(BufHi, v0);
				_mm_storer_ps(BufLo, v1);
				BufLo += 4;
			}
#else
			for(i=0;i<N/2;i++) {
				BufHi--;
				float t = *BufLo;
				*BufLo = *BufHi;
				*BufHi = t;
				BufLo++;
			}
#endif
		}
	}
	Fourier_DCT4(BufOut, BufTmp, N);
}

/**************************************/
//! EOF
/**************************************/
