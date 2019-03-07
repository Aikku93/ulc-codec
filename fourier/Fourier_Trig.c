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

8th order cosine approximation
Estimates Cos[x*Pi/2], x=[0,1]
RMSE=0.000000114459 (PSNR=138.8dB)
------------------------------

a[x_] := Cos[x*Pi/2];
b[x_] := (x^2 - 1)*((c1 + (c2 + c3*x^2)*x^2)*x^2 - 1);
Solve[{
  a[1/4] == b[1/4],
  a[1/2] == b[1/2],
  Integrate[a[x] - b[x], {x, 0, 7/8}] == 0
}, {c1, c2, c3}]

This gives:
  c1 ~= 0.23
  c2 ~= -0.02
  c3 ~= 0.00087

We generally only call SinCos for x/N with small x,
so we minimize error in those regions rather than
aiming for a best fit over the whole range.
However, the MDCT/IMDCT routines also use x/N+0.5,
so we also need some optimization for that area.

!*/
/**************************************/

#define C1 ( 0x1.DE9E45DE9FBFFp-3)
#define C2 (-0x1.4716e4380e769p-6)
#define C3 ( 0x1.C9344264DB172p-11)

/**************************************/
#if defined(__AVX__)
void Fourier_SinCosAVX(__m256 x, __m256 *Sin, __m256 *Cos) {
	__m256 sx = _mm256_sub_ps(x, _mm256_set1_ps(1.0f));
	__m256 cx = x;
	__m256 s, c;
	sx = _mm256_mul_ps(sx, sx);
	cx = _mm256_mul_ps(cx, cx);
	s  = _mm256_mul_ps(sx, _mm256_set1_ps(C3));
	c  = _mm256_mul_ps(cx, _mm256_set1_ps(C3));
	s  = _mm256_add_ps(s,  _mm256_set1_ps(C2));
	c  = _mm256_add_ps(c,  _mm256_set1_ps(C2));
	s  = _mm256_mul_ps(s,  sx);
	c  = _mm256_mul_ps(c,  cx);
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
	s  = _mm_mul_ps(sx, _mm_set1_ps(C3));
	c  = _mm_mul_ps(cx, _mm_set1_ps(C3));
	s  = _mm_add_ps(s,  _mm_set1_ps(C2));
	c  = _mm_add_ps(c,  _mm_set1_ps(C2));
	s  = _mm_mul_ps(s,  sx);
	c  = _mm_mul_ps(c,  cx);
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
	sx *= sx; *Sin = (sx - 1.0f)*((C1 + (C2 + C3*sx)*sx)*sx - 1.0f);
	cx *= cx; *Cos = (cx - 1.0f)*((C1 + (C2 + C3*cx)*cx)*cx - 1.0f);
}

/**************************************/
//! EOF
/**************************************/
