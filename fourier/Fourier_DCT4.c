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

//! DCT-IV (N=8)
static void DCT4_8(float *x) {
	const float sqrt1_2 = 0x1.6A09E6p-1f;
	const float c1_3 = 0x1.D906BDp-1f, s1_3 = 0x1.87DE2Ap-2f;
	const float c1_5 = 0x1.FD88DAp-1f, s1_5 = 0x1.917A6Cp-4f;
	const float c3_5 = 0x1.E9F415p-1f, s3_5 = 0x1.294063p-2f;
	const float c5_5 = 0x1.C38B2Fp-1f, s5_5 = 0x1.E2B5D4p-2f;
	const float c7_5 = 0x1.8BC807p-1f, s7_5 = 0x1.44CF32p-1f;
#if defined(__SSE__)
	__m128 x0, x1;
	__m128 a0, a1;

	x0 = _mm_load_ps (x+0);
	x1 = _mm_loadr_ps(x+4);
	a0 = _mm_mul_ps(x0, _mm_setr_ps( c1_5,  c3_5,  c5_5,  c7_5));
	a1 = _mm_mul_ps(x0, _mm_setr_ps( s1_5, -s3_5,  s5_5, -s7_5));
#if defined(__FMA__)
	x0 = _mm_fmadd_ps(x1, _mm_setr_ps( s1_5,  s3_5,  s5_5,  s7_5), a0);
	x1 = _mm_fmadd_ps(x1, _mm_setr_ps(-c1_5,  c3_5, -c5_5,  c7_5), a1);
#else
	x0 = _mm_mul_ps(x1, _mm_setr_ps( s1_5,  s3_5,  s5_5,  s7_5));
	x1 = _mm_mul_ps(x1, _mm_setr_ps(-c1_5,  c3_5, -c5_5,  c7_5));
	x0 = _mm_add_ps(x0, a0);
	x1 = _mm_add_ps(x1, a1);
#endif
	a0 = _mm_shuffle_ps(x0, x1, 0xE4);
	a1 = _mm_shuffle_ps(x0, x1, 0x1B);
	x0 = _mm_add_ps(a0, a1);
	x1 = _mm_sub_ps(a0, a1);

	a0 = _mm_shuffle_ps(x0, x0, 0xB1);
	a0 = _mm_xor_ps(a0, _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f));
	x0 = _mm_add_ps(x0, a0);
	a0 = _mm_shuffle_ps(x1, x1, 0xCC);
	a1 = _mm_shuffle_ps(x1, x1, 0x99);
#if defined(__FMA__)
	x1 = _mm_mul_ps  (a0, _mm_setr_ps(c1_3, c1_3,  s1_3,  s1_3));
	x1 = _mm_fmadd_ps(a1, _mm_setr_ps(s1_3, s1_3, -c1_3, -c1_3), x1);
#else
	a0 = _mm_mul_ps(a0, _mm_setr_ps(c1_3, c1_3,  s1_3,  s1_3));
	a1 = _mm_mul_ps(a1, _mm_setr_ps(s1_3, s1_3, -c1_3, -c1_3));
	x1 = _mm_add_ps(a0, a1);
#endif
	a0 = _mm_shuffle_ps(x0, x0, 0x6C);
	a0 = _mm_xor_ps(a0, _mm_setr_ps( 0.0f, 0.0f,  0.0f, -0.0f));
	x0 = _mm_add_ps(x0, a0);
	x0 = _mm_mul_ps(x0, _mm_setr_ps(0.5f, -sqrt1_2, 0.5f, sqrt1_2));
	a1 = _mm_shuffle_ps(x1, x1, 0x1B);
	a1 = _mm_xor_ps(a1, _mm_setr_ps(-0.0f, 0.0f, -0.0f,  0.0f));
	x1 = _mm_add_ps(x1, a1);
	a0 = _mm_shuffle_ps(x0, x1, 0xCC);
	a1 = _mm_shuffle_ps(x0, x1, 0x99);
	x0 = _mm_shuffle_ps(a0, a0, 0x78);
	x1 = _mm_shuffle_ps(a1, a1, 0x6C);
	_mm_store_ps(x+0, x0);
	_mm_store_ps(x+4, x1);
#else
	//! First stage (rotation butterflies; DCT4_8)
	float ax =  c1_5*x[0] + s1_5*x[7];
	float ay =  s1_5*x[0] - c1_5*x[7];
	float bx =  c3_5*x[1] + s3_5*x[6];
	float by = -s3_5*x[1] + c3_5*x[6];
	float cx =  c5_5*x[2] + s5_5*x[5];
	float cy =  s5_5*x[2] - c5_5*x[5];
	float dx =  c7_5*x[3] + s7_5*x[4];
	float dy = -s7_5*x[3] + c7_5*x[4];

	//! Second stage (butterflies; DCT2_4)
	float saxdx = ax + dx;
	float daxdx = ax - dx;
	float sbxcx = bx + cx;
	float dbxcx = bx - cx;
	float sdyay = dy + ay;
	float ddyay = dy - ay;
	float scyby = cy + by;
	float dcyby = cy - by;

	//! Third stage (rotation butterflies; DCT2_2, DCT4_2)
	float sx =      saxdx +      sbxcx;
	float sy =      saxdx -      sbxcx;
	float tx = c1_3*daxdx + s1_3*dbxcx;
	float ty = s1_3*daxdx - c1_3*dbxcx;
	float ux =      sdyay +      scyby;
	float uy =      sdyay -      scyby;
	float vx = c1_3*ddyay + s1_3*dcyby;
	float vy = s1_3*ddyay - c1_3*dcyby;

	//! Permute and final DCT4 stage
	x[0] = sx;
	x[1] = (tx - vy);
	x[2] = (tx + vy);
	x[3] = (sy + uy) * sqrt1_2;
	x[4] = (sy - uy) * sqrt1_2;
	x[5] = (ty - vx);
	x[6] = (ty + vx);
	x[7] = ux;
#endif
}

/**************************************/

void Fourier_DCT4(float *Buf, float *Tmp, int N) {
	int i;

	//! Stop condition
	if(N == 8) {
		DCT4_8(Buf);
		return;
	}

	//! Perform rotation butterflies
	//!  u = R_n.x
	{
		float Ns = 1.0f / N;
		const float *SrcLo = Buf;
		const float *SrcHi = Buf + N;
		      float *DstLo = Tmp;
		      float *DstHi = Tmp + N/2;
#if defined(__AVX__)
		__m256 a, b;
		__m256 t0, t1;
		__m256 wc, ws; {
			float c, s;
			Fourier_SinCos(8.0f * Ns, &s, &c);
			wc = _mm256_set1_ps(c);
			ws = _mm256_set1_ps(s);
		}
		__m256 c, s; {
			c = _mm256_setr_ps(0.5f, 1.5f, 2.5f, 3.5f, 4.5f, 5.5f, 6.5f, 7.5f);
			c = _mm256_mul_ps(c, _mm256_set1_ps(Ns));
			Fourier_SinCosAVX(c, &s, &c);
		}

		for(i=0;i<N/2;i+=8) {
			SrcHi -= 8;
			b = _mm256_load_ps(SrcHi);
			a = _mm256_load_ps(SrcLo);
			SrcLo += 8;
			b = _mm256_shuffle_ps(b, b, 0x1B);
			b = _mm256_permute2f128_ps(b, b, 0x01);
#if defined(__FMA__)
			t1 = _mm256_mul_ps(s, a);
			t0 = _mm256_mul_ps(c, a);
			t1 = _mm256_fnmadd_ps(c, b, t1);
			t0 = _mm256_fmadd_ps (s, b, t0);
#else
			t1 = _mm256_sub_ps(_mm256_mul_ps(s, a), _mm256_mul_ps(c, b));
			t0 = _mm256_add_ps(_mm256_mul_ps(c, a), _mm256_mul_ps(s, b));
#endif
			t1 = _mm256_xor_ps(t1, _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f));
			_mm256_store_ps(DstLo, t0); DstLo += 8;
			_mm256_store_ps(DstHi, t1); DstHi += 8;
#if defined(__FMA__)
			t0 = _mm256_mul_ps(wc, c);
			t1 = _mm256_mul_ps(ws, c);
			t0 = _mm256_fnmadd_ps(ws, s, t0);
			t1 = _mm256_fmadd_ps (wc, s, t1);
#else
			t0 = _mm256_sub_ps(_mm256_mul_ps(wc, c), _mm256_mul_ps(ws, s));
			t1 = _mm256_add_ps(_mm256_mul_ps(ws, c), _mm256_mul_ps(wc, s));
#endif
			c = t0;
			s = t1;
		}
#elif defined(__SSE__)
		__m128 a, b;
		__m128 t0, t1;
		__m128 wc, ws; {
			float c, s;
			Fourier_SinCos(4.0f * Ns, &s, &c);
			wc = _mm_set1_ps(c);
			ws = _mm_set1_ps(s);
		}
		__m128 c, s; {
			c = _mm_setr_ps(0.5f, 1.5f, 2.5f, 3.5f);
			c = _mm_mul_ps(c, _mm_set1_ps(Ns));
			Fourier_SinCosSSE(c, &s, &c);
		}

		for(i=0;i<N/2;i+=4) {
			SrcHi -= 4;
			b = _mm_loadr_ps(SrcHi);
			a = _mm_load_ps(SrcLo);
			SrcLo += 4;
#if defined(__FMA__)
			t1 = _mm_mul_ps(s, a);
			t0 = _mm_mul_ps(c, a);
			t1 = _mm_fnmadd_ps(c, b, t1);
			t0 = _mm_fmadd_ps (s, b, t0);
#else
			t1 = _mm_sub_ps(_mm_mul_ps(s, a), _mm_mul_ps(c, b));
			t0 = _mm_add_ps(_mm_mul_ps(c, a), _mm_mul_ps(s, b));
#endif
			t1 = _mm_xor_ps(t1, _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f));

			_mm_store_ps(DstLo, t0); DstLo += 4;
			_mm_store_ps(DstHi, t1); DstHi += 4;
#if defined(__FMA__)
			t0 = _mm256_mul_ps(wc, c);
			t1 = _mm256_mul_ps(ws, c);
			t0 = _mm256_fnmadd_ps(ws, s, t0);
			t1 = _mm256_fmadd_ps (wc, s, t1);
#else
			t0 = _mm_sub_ps(_mm_mul_ps(wc, c), _mm_mul_ps(ws, s));
			t1 = _mm_add_ps(_mm_mul_ps(ws, c), _mm_mul_ps(wc, s));
#endif
			c = t0;
			s = t1;
		}
#else
#define UPDATE_SINCOS() \
	do { \
		_c = wc*c - ws*s; \
		_s = ws*c + wc*s; \
		c = _c; \
		s = _s; \
	} while(0)
		float a, b;
		float c, s, _c, _s;
		float wc, ws;
		Fourier_SinCos(0.5f * Ns, &s,  &c);
		Fourier_SinCos(1.0f * Ns, &ws, &wc);
		for(i=0;i<N/2;i+=2) {
			a = *SrcLo++;
			b = *--SrcHi;
			*DstLo++ =  c*a + s*b;
			*DstHi++ =  s*a - c*b;
			UPDATE_SINCOS();

			a = *SrcLo++;
			b = *--SrcHi;
			*DstLo++ =  c*a + s*b;
			*DstHi++ = -s*a + c*b;
			UPDATE_SINCOS();
		}
#undef UPDATE_SINCOS
#endif
	}

	//! Perform recursion
	//!  z1 = cos2([u_j][j=0..n/2-1],n/2)
	//!  z2 = cos2([u_j][j=n/2..n-1],n/2)
	Fourier_DCT2(Tmp,       Buf,       N/2);
	Fourier_DCT2(Tmp + N/2, Buf + N/2, N/2);

	//! Combine
	//!  w = U_n.(z1^T, z2^T)^T
	//!  y = (P_n)^T.w
	//! TODO: Avoid unaligned load/store in SSE/AVX code?
	{
		const float *TmpLo = Tmp;
		const float *TmpHi = Tmp + N;
		      float *Dst   = Buf;
#if defined(__AVX__)
		__m256 a, b;
		__m256 t0, t1;

		*Dst++ = *TmpLo++;
		for(i=0;i<N/2-8;i+=8) {
			TmpHi -= 8; b = _mm256_load_ps(TmpHi);
			a = _mm256_loadu_ps(TmpLo); TmpLo += 8;
			b = _mm256_shuffle_ps(b, b, 0x1B);
			b = _mm256_permute2f128_ps(b, b, 0x01);

			t0 = _mm256_add_ps(a, b);
			t1 = _mm256_sub_ps(a, b);
			a  = _mm256_unpacklo_ps(t0, t1);
			b  = _mm256_unpackhi_ps(t0, t1);
			t0 = _mm256_permute2f128_ps(a, b, 0x20);
			t1 = _mm256_permute2f128_ps(a, b, 0x31);
			_mm256_storeu_ps(Dst, t0); Dst += 8;
			_mm256_storeu_ps(Dst, t1); Dst += 8;
		}
		{
			TmpHi -= 8; b = _mm256_load_ps(TmpHi);
			a = _mm256_loadu_ps(TmpLo); TmpLo += 8;
			_mm_store_ss(Dst + 14, _mm256_castps256_ps128(b));
			b = _mm256_shuffle_ps(b, b, 0x1B);
			b = _mm256_permute2f128_ps(b, b, 0x01);

			t0 = _mm256_add_ps(a, b);
			t1 = _mm256_sub_ps(a, b);
			a = _mm256_unpacklo_ps(t0, t1);
			b = _mm256_unpackhi_ps(t0, t1);
			t0 = _mm256_permute2f128_ps(a, b, 0x20);
			t1 = _mm256_permute2f128_ps(a, b, 0x31);
			_mm256_storeu_ps(Dst, t0);                              Dst += 8;
			_mm_storeu_ps(Dst, _mm256_castps256_ps128(t1));         Dst += 4;
			t1 = _mm256_permute2f128_ps(t1, t1, 0x01);
			_mm_storel_pi((__m64*)Dst, _mm256_castps256_ps128(t1)); Dst += 2;
		}
#elif defined(__SSE__)
		__m128 a, b;
		__m128 t0, t1;

		*Dst++ = *TmpLo++;
		for(i=0;i<N/2-4;i+=4) {
			TmpHi -= 4; b = _mm_loadr_ps(TmpHi);
			a = _mm_loadu_ps(TmpLo); TmpLo += 4;

			t0 = _mm_add_ps(a, b);
			t1 = _mm_sub_ps(a, b);
			a = _mm_unpacklo_ps(t0, t1);
			b = _mm_unpackhi_ps(t0, t1);
			_mm_storeu_ps(Dst, a); Dst += 4;
			_mm_storeu_ps(Dst, b); Dst += 4;
		}
		{
			TmpHi -= 4; b = _mm_load_ps(TmpHi);
			a = _mm_loadu_ps(TmpLo); TmpLo += 4;
			_mm_store_ss(Dst + 6, b);
			b = _mm_shuffle_ps(b, b, 0x1B);

			t0 = _mm_add_ps(a, b);
			t1 = _mm_sub_ps(a, b);
			a = _mm_unpacklo_ps(t0, t1);
			b = _mm_unpackhi_ps(t0, t1);
			_mm_storeu_ps(Dst, a);         Dst += 4;
			_mm_storel_pi((__m64*)Dst, b); Dst += 2;
		}
#else
		float a, b;

		*Dst++ = *TmpLo++;
		for(i=0;i<N/2-1;i++) {
			a = *TmpLo++;
			b = *--TmpHi;
			*Dst++ = a + b;
			*Dst++ = a - b;
		}
		*Dst++ = *--TmpHi;
#endif
	}
}

/**************************************/
//! EOF
/**************************************/
