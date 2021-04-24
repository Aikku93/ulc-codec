/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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

//! DCT-III (N=8)
static void DCT3_8(float *x) {
	const float sqrt1_2 = 0x1.6A09E6p-1f;
	const float c1_4 = 0x1.F6297Dp-1f, s1_4 = 0x1.8F8B84p-3f;
	const float c3_4 = 0x1.A9B663p-1f, s3_4 = 0x1.1C73B4p-1f;
	const float c6_4 = 0x1.87DE2Ap-2f, s6_4 = 0x1.D906BDp-1f;

	float a0 = x[0];
	float b0 = x[4] * sqrt1_2;
	float c0 = x[2];
	float d0 = x[6];
	float a1 = (x[1] + x[7]) * sqrt1_2;
	float d1 = (x[1] - x[7]) * sqrt1_2;
	float b1 = x[5];
	float c1 = x[3];

	float ss07s34 = a0*0.5f + b0;
	float ss16s25 = a0*0.5f - b0;
	float ds16s25 = c6_4*c0 - s6_4*d0;
	float ds07s34 = s6_4*c0 + c6_4*d0;
	float d34d07y = a1 + c1;
	float d25d16x = a1 - c1;
	float d34d07x = d1 + b1;
	float d25d16y = d1 - b1;

	float s07 = ss07s34 + ds07s34;
	float s34 = ss07s34 - ds07s34;
	float s16 = ss16s25 + ds16s25;
	float s25 = ss16s25 - ds16s25;
	float d34 = c3_4*d34d07x - s3_4*d34d07y;
	float d07 = s3_4*d34d07x + c3_4*d34d07y;
	float d25 = c1_4*d25d16x - s1_4*d25d16y;
	float d16 = s1_4*d25d16x + c1_4*d25d16y;

	x[0] = s07 + d07;
	x[7] = s07 - d07;
	x[1] = s16 + d16;
	x[6] = s16 - d16;
	x[2] = s25 + d25;
	x[5] = s25 - d25;
	x[3] = s34 + d34;
	x[4] = s34 - d34;
}

/**************************************/

void Fourier_DCT3(float *Buf, float *Tmp, int N) {
	int i;

	//! Stop condition
	if(N == 8) {
		DCT3_8(Buf);
		return;
	}

	{
		      float *DstLo = Tmp;
		      float *DstHi = Tmp + N/2;
		const float *Src   = Buf;
#if defined(__AVX__)
		__m256 a0, b0;
		__m256 a1, b1;
		for(i=0;i<N/2;i+=8) {
			a0 = _mm256_load_ps(Src); Src += 8;
			b0 = _mm256_load_ps(Src); Src += 8;
			a1 = _mm256_permute2f128_ps(a0, b0, 0x20);
			b1 = _mm256_permute2f128_ps(a0, b0, 0x31);
			a0 = _mm256_shuffle_ps(a1, b1, 0x88);
			b0 = _mm256_shuffle_ps(a1, b1, 0xDD);
			_mm256_store_ps(DstLo, a0); DstLo += 8;
			_mm256_store_ps(DstHi, b0); DstHi += 8;
		}
#elif defined(__SSE__)
		__m128 a0, b0;
		__m128 a1, b1;
		for(i=0;i<N/2;i+=4) {
			a0 = _mm_load_ps(Src); Src += 4;
			b0 = _mm_load_ps(Src); Src += 4;
			a1 = _mm_shuffle_ps(a0, b0, 0x88);
			b1 = _mm_shuffle_ps(a0, b0, 0xDD);
			_mm_store_ps(DstLo, a1); DstLo += 4;
			_mm_store_ps(DstHi, b1); DstHi += 4;
		}
#else
		float a, b;
		for(i=0;i<N/2;i++) {
			a = *Src++;
			b = *Src++;
			*DstLo++ = a;
			*DstHi++ = b;
		}
#endif
	}

	Fourier_DCT3 (Tmp,       Buf,       N/2);
	Fourier_DCT4T(Tmp + N/2, Buf + N/2, N/2);

	{
		const float *SrcLo = Tmp;
		const float *SrcHi = Tmp + N/2;
		      float *DstLo = Buf;
		      float *DstHi = Buf + N;
#if defined(__AVX__)
		__m256 a, b;
		__m256 s, d;
		for(i=0;i<N/2;i+=8) {
			a = _mm256_load_ps(SrcLo); SrcLo += 8;
			b = _mm256_load_ps(SrcHi); SrcHi += 8;
			s = _mm256_add_ps(a, b);
			d = _mm256_sub_ps(a, b);
			d = _mm256_shuffle_ps(d, d, 0x1B);
			d = _mm256_permute2f128_ps(d, d, 0x01);
			_mm256_store_ps(DstLo, s); DstLo += 8;
			DstHi -= 8; _mm256_store_ps(DstHi, d);
		}
#elif defined(__SSE__)
		__m128 a, b;
		__m128 s, d;
		for(i=0;i<N/2;i+=4) {
			a = _mm_load_ps(SrcLo); SrcLo += 4;
			b = _mm_load_ps(SrcHi); SrcHi += 4;
			s = _mm_add_ps(a, b);
			d = _mm_sub_ps(a, b);
			_mm_store_ps(DstLo, s); DstLo += 4;
			DstHi -= 4; _mm_storer_ps(DstHi, d);
		}
#else
		float a, b;
		for(i=0;i<N/2;i++) {
			a = *SrcLo++;
			b = *SrcHi++;
			*DstLo++ = a + b;
			*--DstHi = a - b;
		}
#endif
	}
}

/**************************************/
//! EOF
/**************************************/
