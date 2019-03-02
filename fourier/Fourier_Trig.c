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
/*!

6th order cosine approximation
Estimates Cos[x*Pi/2], x=[0,1]
RMSE=0.0000127281 (PSNR=97.9dB)
------------------------------

a[x_] := Cos[x*Pi/2];
b[x_] := (x^2 - 1)*((c1 + c2*x^2)*x^2 - 1);
Solve[{
  a'[1] == b'[1],
  Integrate[(a[x] - b[x]), {x, 0, 1}] == 0
  }, {c1, c2}]

This gives:
  c1 = (420 - 128*Pi - 3*Pi^2) / (-16*Pi)
  c2 = 105/(4*Pi) - 7*Pi/16 - 7

!*/
/**************************************/

#define C1 ( 0x1.DE08378EAE872p-3) //! (420 - 128*Pi - 3*Pi^2) / (-16*Pi)
#define C2 (-0x1.343864FDCE740p-6) //! 105/(4*Pi) - 7*Pi/16 - 7

/**************************************/
#if defined(__AVX__)
void Fourier_SinCosAVX(__m256 x, __m256 *Sin, __m256 *Cos) {
	__m256 sx = _mm256_sub_ps(x, _mm256_set1_ps(1.0f));
	__m256 cx = x;
	__m256 s, c;
	sx = _mm256_mul_ps(sx, sx);
	cx = _mm256_mul_ps(cx, cx);
	s  = _mm256_mul_ps(sx, _mm256_set1_ps(C2));
	c  = _mm256_mul_ps(cx, _mm256_set1_ps(C2));
	s  = _mm256_add_ps(s,  _mm256_set1_ps(C1));
	c  = _mm256_add_ps(c,  _mm256_set1_ps(C1));
	s  = _mm256_mul_ps(s,  sx);
	c  = _mm256_mul_ps(c,  cx);
	s  = _mm256_sub_ps(s,  _mm256_set1_ps(1.0f));
	c  = _mm256_sub_ps(c,  _mm256_set1_ps(1.0f));
	sx = _mm256_sub_ps(sx, _mm256_set1_ps(1.0f));
	cx = _mm256_sub_ps(cx, _mm256_set1_ps(1.0f));
	s  = _mm256_mul_ps(s,  sx);
	c  = _mm256_mul_ps(c,  cx);
	*Sin = s;
	*Cos = c;
}
#endif
/**************************************/
#if defined(__SSE__)
void Fourier_SinCosSSE(__m128 x, __m128 *Sin, __m128 *Cos) {
	__m128 sx = _mm_sub_ps(x, _mm_set1_ps(1.0f));
	__m128 cx = x;
	__m128 s, c;
	sx = _mm_mul_ps(sx, sx);
	cx = _mm_mul_ps(cx, cx);
	s  = _mm_mul_ps(sx, _mm_set1_ps(C2));
	c  = _mm_mul_ps(cx, _mm_set1_ps(C2));
	s  = _mm_add_ps(s,  _mm_set1_ps(C1));
	c  = _mm_add_ps(c,  _mm_set1_ps(C1));
	s  = _mm_mul_ps(s,  sx);
	c  = _mm_mul_ps(c,  cx);
	s  = _mm_sub_ps(s,  _mm_set1_ps(1.0f));
	c  = _mm_sub_ps(c,  _mm_set1_ps(1.0f));
	sx = _mm_sub_ps(sx, _mm_set1_ps(1.0f));
	cx = _mm_sub_ps(cx, _mm_set1_ps(1.0f));
	s  = _mm_mul_ps(s,  sx);
	c  = _mm_mul_ps(c,  cx);
	*Sin = s;
	*Cos = c;
}
#endif
/**************************************/

void Fourier_SinCos(float x, float *Sin, float *Cos) {
	float sx = x - 1.0f, cx = x; //! Phase-change sine to cosine
	sx *= sx; *Sin = (sx - 1.0f)*((C1 + C2*sx)*sx - 1.0f);
	cx *= cx; *Cos = (cx - 1.0f)*((C1 + C2*cx)*cx - 1.0f);
}

/**************************************/
//! EOF
/**************************************/
