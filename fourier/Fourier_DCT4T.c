/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <stddef.h>
#if defined(__AVX__) || defined(__FMA__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include "Fourier.h"
#include "FourierHelper.h"
/**************************************/

//! DCT-IV (N=8)
static void DCT4T_8(float *x) {
	FOURIER_ASSUME_ALIGNED(x, 32);

	const float sqrt1_2 = 0x1.6A09E6p-1f;
	const float c1_3 = 0x1.D906BCp-1f, s1_3 = 0x1.87DE2Ap-2f;
	const float c1_5 = 0x1.FD88DAp-1f, s1_5 = 0x1.917A6Cp-4f;
	const float c3_5 = 0x1.E9F416p-1f, s3_5 = 0x1.294062p-2f;
	const float c5_5 = 0x1.C38B30p-1f, s5_5 = 0x1.E2B5D4p-2f;
	const float c7_5 = 0x1.8BC806p-1f, s7_5 = 0x1.44CF32p-1f;

	float sx = x[0];
	float tx = (x[2] + x[1]);
	float vy = (x[2] - x[1]);
	float sy = (x[3] + x[4]) * sqrt1_2;
	float uy = (x[3] - x[4]) * sqrt1_2;
	float ty = (x[6] + x[5]);
	float vx = (x[6] - x[5]);
	float ux = x[7];

	float saxdx = sx + sy;
	float sbxcx = sx - sy;
	float daxdx = c1_3*tx + s1_3*ty;
	float dbxcx = s1_3*tx - c1_3*ty;
	float sdyay = ux + uy;
	float scyby = ux - uy;
	float ddyay = c1_3*vx + s1_3*vy;
	float dcyby = s1_3*vx - c1_3*vy;

	float ax = saxdx + daxdx;
	float dx = saxdx - daxdx;
	float bx = sbxcx + dbxcx;
	float cx = sbxcx - dbxcx;
	float dy = sdyay + ddyay;
	float ay = sdyay - ddyay;
	float cy = scyby + dcyby;
	float by = scyby - dcyby;

	x[0] = c1_5*ax + s1_5*ay;
	x[7] = s1_5*ax - c1_5*ay;
	x[1] = c3_5*bx - s3_5*by;
	x[6] = s3_5*bx + c3_5*by;
	x[2] = c5_5*cx + s5_5*cy;
	x[5] = s5_5*cx - c5_5*cy;
	x[3] = c7_5*dx - s7_5*dy;
	x[4] = s7_5*dx + c7_5*dy;
}

/**************************************/

void Fourier_DCT4T(float *Buf, float *Tmp, int N) {
	int i;
	FOURIER_ASSUME_ALIGNED(Buf, 32);
	FOURIER_ASSUME_ALIGNED(Tmp, 32);
	FOURIER_ASSUME(N >= 8 && N <= 8192);

	//! Stop condition
	if(N == 8) {
		DCT4T_8(Buf);
		return;
	}

	{
		      float *DstLo = Tmp;
		      float *DstHi = Tmp + N;
		const float *Src   = Buf;

		*DstLo++ = *Src++ * 2.0f;
#if defined(__AVX__)
		__m256 a, b;
		__m256 t0, t1;
		for(i=0;i<N/2-8;i+=8) {
			a  = _mm256_loadu_ps(Src); Src += 8;
			b  = _mm256_loadu_ps(Src); Src += 8;
			t0 = _mm256_permute2f128_ps(a, b, 0x20);
			t1 = _mm256_permute2f128_ps(a, b, 0x31);
			a  = _mm256_shuffle_ps(t0, t1, 0x88);
			b  = _mm256_shuffle_ps(t0, t1, 0xDD);
			t0 = _mm256_add_ps(a, b);
			t1 = _mm256_sub_ps(a, b);
			t1 = _mm256_shuffle_ps(t1, t1, 0x1B);
			t1 = _mm256_permute2f128_ps(t1, t1, 0x01);
			_mm256_storeu_ps(DstLo, t0); DstLo += 8;
			DstHi -= 8; _mm256_store_ps(DstHi, t1);
		}
		{
			a  = _mm256_loadu_ps(Src); Src += 8;     //! {1,2,3,4,5,6,7,8}
			b  = _mm256_loadu_ps(Src); Src += 8;     //! {9,A,B,C,D,E,F,X}
			t0 = _mm256_permute2f128_ps(a, b, 0x20); //! {1,2,3,4,9,A,B,C}
			t1 = _mm256_permute2f128_ps(a, b, 0x31); //! {5,6,7,8,D,E,F,X}
			a  = _mm256_shuffle_ps(t0, t1, 0x88);    //! {1,3,5,7,9,B,D,F}
			b  = _mm256_shuffle_ps(t0, t1, 0xDD);    //! {2,4,6,8,A,C,E,X} <- Note the X here
			b  = _mm256_blend_ps(b, a, 0x80);        //! {2,4,6,8,A,C,E,F} <- Note the F here now
			t0 = _mm256_add_ps(a, b);                //! a[7] = 2.0*F
			t1 = _mm256_sub_ps(a, b);                //! a[7] = 0; "garbage"
			t1 = _mm256_shuffle_ps(t1, t1, 0x1B);
			t1 = _mm256_permute2f128_ps(t1, t1, 0x01);
			DstHi -= 8; _mm256_store_ps(DstHi, t1);  //! Stores "garbage" for final DstHi
			_mm256_storeu_ps(DstLo, t0); DstLo += 8; //! Stores correct final DstHi (2.0*F)
		}
#elif defined(__SSE__)
		__m128 a, b;
		__m128 t0, t1;
		for(i=0;i<N/2-4;i+=4) {
			a  = _mm_loadu_ps(Src); Src += 4;
			b  = _mm_loadu_ps(Src); Src += 4;
			t0 = _mm_shuffle_ps(a, b, 0x88);
			t1 = _mm_shuffle_ps(a, b, 0xDD);
			a = _mm_add_ps(t0, t1);
			b = _mm_sub_ps(t0, t1);
			_mm_storeu_ps(DstLo, a); DstLo += 4;
			DstHi -= 4; _mm_storer_ps(DstHi, b);
		}
		{
			a = _mm_loadu_ps(Src); Src += 4;     //! {9,A,B,C}
			b = _mm_loadu_ps(Src); Src += 4;     //! {D,E,F,X}
			t0 = _mm_shuffle_ps(a, b, 0x88);     //! {9,B,D,F}
			t1 = _mm_shuffle_ps(a, b, 0x9D);     //! {A,C,E,F} <- Note the F here
			a = _mm_add_ps(t0, t1);              //! a[3] = 2.0*F
			b = _mm_sub_ps(t0, t1);              //! a[3] = 0; "garbage"
			DstHi -= 4; _mm_storer_ps(DstHi, b); //! Stores "garbage" for final DstHi
			_mm_storeu_ps(DstLo, a); DstLo += 4; //! Stores correct final DstHi (2.0*F)
		}
#else
		float a, b;
		for(i=0;i<N/2-1;i++) {
			a = *Src++;
			b = *Src++;
			*DstLo++ = a + b;
			*--DstHi = a - b;
		}
		*--DstHi = *Src++ * 2.0f;
#endif
	}

	Fourier_DCT3(Tmp,       Buf,       N/2);
	Fourier_DCT3(Tmp + N/2, Buf + N/2, N/2);

	{
		const float *Win   = Fourier_SinTableN(N, NULL);
		const float *SrcLo = Tmp;
		const float *SrcHi = Tmp + N/2;
		      float *DstLo = Buf;
		      float *DstHi = Buf + N;
#if defined(__AVX__)
		__m256 a, b;
		__m256 t0, t1;
		__m256 c, s;
		for(i=0;i<N/2;i+=8) {
			a = _mm256_load_ps(SrcLo); SrcLo += 8;
			b = _mm256_load_ps(SrcHi); SrcHi += 8;
			c  = _mm256_load_ps(Win); Win += 8;
			s  = _mm256_load_ps(Win); Win += 8;
			t0 = _mm256_permute2f128_ps(c, s, 0x20);
			t1 = _mm256_permute2f128_ps(c, s, 0x31);
			c  = _mm256_shuffle_ps(t0, t1, 0x88);
			s  = _mm256_shuffle_ps(t0, t1, 0xDD);
			t0 = _mm256_mul_ps(s, b);
			t1 = _mm256_mul_ps(c, b);
			t0 = _mm256_xor_ps(t0, _mm256_set1_ps(-0.0f));
#if defined(__FMA__)
			t0 = _mm256_fmaddsub_ps(c, a, t0);
			t1 = _mm256_fmaddsub_ps(s, a, t1);
#else
			t0 = _mm256_addsub_ps(_mm256_mul_ps(c, a), t0);
			t1 = _mm256_addsub_ps(_mm256_mul_ps(s, a), t1);
#endif
			t1 = _mm256_shuffle_ps(t1, t1, 0x1B);
			t1 = _mm256_permute2f128_ps(t1, t1, 0x01);
			_mm256_store_ps(DstLo, t0); DstLo += 8;
			DstHi -= 8; _mm256_store_ps(DstHi, t1);
		}
#elif defined(__SSE__)
		__m128 a, b;
		__m128 t0, t1;
		__m128 c, s;
		for(i=0;i<N/2;i+=4) {
			a  = _mm_load_ps(SrcLo); SrcLo += 4;
			b  = _mm_load_ps(SrcHi); SrcHi += 4;
			t0 = _mm_load_ps(Win); Win += 4;
			t1 = _mm_load_ps(Win); Win += 4;
			c  = _mm_shuffle_ps(t0, t1, 0x88);
			s  = _mm_shuffle_ps(t0, t1, 0xDD);
			t0 = _mm_mul_ps(s, b);
			t1 = _mm_mul_ps(c, b);
#if defined(__FMA__)
			t0 = _mm_xor_ps(t0, _mm_set1_ps(-0.0f));
			t0 = _mm_fmaddsub_ps(c, a, t0);
			t1 = _mm_fmaddsub_ps(s, a, t1);
#else
			t0 = _mm_xor_ps(t0, _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f));
			t1 = _mm_xor_ps(t1, _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f));
			t0 = _mm_add_ps(_mm_mul_ps(c, a), t0);
			t1 = _mm_sub_ps(_mm_mul_ps(s, a), t1);
#endif
			_mm_store_ps(DstLo, t0); DstLo += 4;
			DstHi -= 4; _mm_storer_ps(DstHi, t1);
		}
#else
		float a, b;
		float c, s;
		for(i=0;i<N/2;i+=2) {
			a = *SrcLo++;
			b = *SrcHi++;
			c = *Win++;
			s = *Win++;
			*DstLo++ =  c*a + s*b;
			*--DstHi =  s*a - c*b;

			a = *SrcLo++;
			b = *SrcHi++;
			c = *Win++;
			s = *Win++;
			*DstLo++ =  c*a - s*b;
			*--DstHi =  s*a + c*b;
		}
#endif
	}
}

/**************************************/
//! EOF
/**************************************/
