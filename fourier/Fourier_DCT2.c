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

//! DCT-II (N=8)
static void DCT2_8(float *x) {
	const float sqrt1_2 = 0x1.6A09E6p-1f;
	const float c1_4 = 0x1.F6297Dp-1f, s1_4 = 0x1.8F8B84p-3f;
	const float c3_4 = 0x1.A9B663p-1f, s3_4 = 0x1.1C73B4p-1f;
	const float c6_4 = 0x1.87DE2Ap-2f, s6_4 = 0x1.D906BDp-1f;

	//! First stage butterflies (DCT2_8)
	float s07 = x[0]+x[7];
	float d07 = x[0]-x[7];
	float s16 = x[1]+x[6];
	float d16 = x[1]-x[6];
	float s25 = x[2]+x[5];
	float d25 = x[2]-x[5];
	float s34 = x[3]+x[4];
	float d34 = x[3]-x[4];

	//! Second stage (DCT2_4, DCT4_4)
	float ss07s34 = s07+s34;
	float ds07s34 = s07-s34;
	float ss16s25 = s16+s25;
	float ds16s25 = s16-s25;
	float d34d07x =  c3_4*d34 + s3_4*d07;
	float d34d07y = -s3_4*d34 + c3_4*d07;
	float d25d16x =  c1_4*d25 + s1_4*d16;
	float d25d16y = -s1_4*d25 + c1_4*d16;

	//! Third stage (rotation butterflies; DCT2_2, DCT4_2, DCT2_2, DCT2_2)
	float a0 =       ss07s34 +      ss16s25;
	float b0 =       ss07s34 -      ss16s25;
	float c0 =  c6_4*ds16s25 + s6_4*ds07s34;
	float d0 = -s6_4*ds16s25 + c6_4*ds07s34;
	float a1 =       d34d07y +      d25d16x;
	float c1 =       d34d07y -      d25d16x;
	float d1 =       d34d07x +      d25d16y;
	float b1 =       d34d07x -      d25d16y;

	//! Permute and final DCT4 stage
	x[0] = a0;
	x[4] = b0 * sqrt1_2;
	x[2] = c0;
	x[6] = d0;
	x[1] = (a1 + d1) * sqrt1_2;
	x[5] = b1;
	x[3] = c1;
	x[7] = (a1 - d1) * sqrt1_2;
}

/**************************************/

void Fourier_DCT2(float *Buf, float *Tmp, int N) {
	int i;

	//! Stop condition
	if(N == 8) {
		DCT2_8(Buf);
		return;
	}

	//! Perform butterflies
	//!  u = H_n.x
	{
		const float *SrcLo = Buf;
		const float *SrcHi = Buf + N;
		      float *DstLo = Tmp;
		      float *DstHi = Tmp + N/2;
#if defined(__AVX__)
		__m256 a, b;
		__m256 s, d;
		for(i=0;i<N/2;i+=8) {
			SrcHi -= 8; b = _mm256_load_ps(SrcHi);
			a = _mm256_load_ps(SrcLo); SrcLo += 8;
			b = _mm256_shuffle_ps(b, b, 0x1B);
			b = _mm256_permute2f128_ps(b, b, 0x01);
			s = _mm256_add_ps(a, b);
			d = _mm256_sub_ps(a, b);
			_mm256_store_ps(DstLo, s); DstLo += 8;
			_mm256_store_ps(DstHi, d); DstHi += 8;
		}
#elif defined(__SSE__)
		__m128 a0, b0;
		__m128 a1, b1;
		__m128 s0, d0;
		__m128 s1, d1;
		for(i=0;i<N/2;i+=8) {
			SrcHi -= 8;
			b0 = _mm_loadr_ps(SrcHi + 4);
			b1 = _mm_loadr_ps(SrcHi + 0);
			a0 = _mm_load_ps (SrcLo + 0);
			a1 = _mm_load_ps (SrcLo + 4);
			SrcLo += 8;
			s0 = _mm_add_ps(a0, b0);
			d0 = _mm_sub_ps(a0, b0);
			s1 = _mm_add_ps(a1, b1);
			d1 = _mm_sub_ps(a1, b1);
			_mm_store_ps(DstLo + 0, s0);
			_mm_store_ps(DstLo + 4, s1); DstLo += 8;
			_mm_store_ps(DstHi + 0, d0);
			_mm_store_ps(DstHi + 4, d1); DstHi += 8;
		}
#else
		float a, b;
		for(i=0;i<N/2;i++) {
			a = *SrcLo++;
			b = *--SrcHi;
			*DstLo++ = a + b;
			*DstHi++ = a - b;
		}
#endif
	}

	//! Perform recursion
	//!  z1 = cos2([u_j][j=0..n/2-1],n/2)
	//!  z2 = cos4([u_j][j=n/2..n-1],n/2)
	Fourier_DCT2(Tmp,       Buf,       N/2);
	Fourier_DCT4(Tmp + N/2, Buf + N/2, N/2);

	//! Combine
	//!  y = (P_n)^T.(z1^T, z2^T)^T
	{
		const float *SrcLo = Tmp;
		const float *SrcHi = Tmp + N/2;
		      float *Dst   = Buf;
#if defined(__AVX__)
		__m256 a, b;
		__m256 a0, b0;
		__m256 t0, t1;
		for(i=0;i<N/2;i+=8) {
			a0 = _mm256_load_ps(SrcLo); SrcLo += 8;
			b0 = _mm256_load_ps(SrcHi); SrcHi += 8;
			t0 = _mm256_unpacklo_ps(a0, b0);
			t1 = _mm256_unpackhi_ps(a0, b0);
			a  = _mm256_permute2f128_ps(t0, t1, 0x20);
			b  = _mm256_permute2f128_ps(t0, t1, 0x31);
			_mm256_store_ps(Dst + 0, a);
			_mm256_store_ps(Dst + 8, b); Dst += 16;
		}
#elif defined(__SSE__)
		__m256 a, b, c, d;
		__m128 a0, a1;
		__m128 b0, b1;
		for(i=0;i<N/2;i+=8) {
			a0 = _mm_load_ps(SrcLo + 0);
			b0 = _mm_load_ps(SrcHi + 0);
			a1 = _mm_load_ps(SrcLo + 4); SrcLo += 8;
			b1 = _mm_load_ps(SrcHi + 4); SrcHi += 8;
			a = _mm_unpacklo_ps(a0, b0);
			b = _mm_unpackhi_ps(a0, b0);
			c = _mm_unpacklo_ps(a1, b1);
			d = _mm_unpackhi_ps(a1, b1);
			_mm_store_ps(Dst +  0, a);
			_mm_store_ps(Dst +  4, b);
			_mm_store_ps(Dst +  8, c);
			_mm_store_ps(Dst + 12, d); Dst += 16;
		}
#else
		float a, b;
		for(i=0;i<N/2;i++) {
			a = *SrcLo++;
			b = *SrcHi++;
			*Dst++ = a;
			*Dst++ = b;
		}
#endif
	}
}

/**************************************/
//! EOF
/**************************************/
