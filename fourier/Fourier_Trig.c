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
/*!

8th order cosine approximation
Estimates Cos[x*Pi/2], x=[0,1]
RMSE=0.000000040504495 (PSNR=147.9dB)
------------------------------

a[x_] := Cos[x*Pi/2];
b[x_] := (x^2 - 1)*((c1 + (c2 + c3*x^2)*x^2)*x^2 - 1);
Solve[{
  a[430/1000] == b[430/1000],
  a[710/1000] == b[710/1000],
  a[895/1000] == b[895/1000]
}, {c1, c2, c3}]

This gives:
  c1 ~= 0.23
  c2 ~= -0.02
  c3 ~= 0.00086

Although SinCos is generally called as x/N with small
x and large N, the same equation is used for BOTH
sine /and/ cosine, so we must optimize at both ends.
And since DCT4 gradually closes in on x=0.5, we must
optimize the whole space anyway.

!*/
/**************************************/

#define C1 ( 0x1.DE9D73D230B71p-3)
#define C2 (-0x1.46EAA650565C1p-6)
#define C3 ( 0x1.C20AB73E45F6Bp-11)

/**************************************/
#if defined(__AVX__)
void Fourier_SinCosAVX(__m256 x, __m256 *Sin, __m256 *Cos) {
	__m256 sx = _mm256_sub_ps(x, _mm256_set1_ps(1.0f));
	__m256 cx = x;
	__m256 s, c;
	sx = _mm256_mul_ps(sx, sx);
	cx = _mm256_mul_ps(cx, cx);
#if defined(__FMA__)
	s  = _mm256_set1_ps(C3);
	c  = _mm256_set1_ps(C3);
	s  = _mm256_fmadd_ps(s, sx, _mm256_set1_ps(C2));
	c  = _mm256_fmadd_ps(c, cx, _mm256_set1_ps(C2));
	s  = _mm256_fmadd_ps(s, sx, _mm256_set1_ps(C1));
	c  = _mm256_fmadd_ps(c, cx, _mm256_set1_ps(C1));
	s  = _mm256_fmadd_ps(s, sx, _mm256_set1_ps(-1.0f));
	c  = _mm256_fmadd_ps(c, cx, _mm256_set1_ps(-1.0f));
	s  = _mm256_fmsub_ps(s, sx, s);
	c  = _mm256_fmsub_ps(c, cx, c);
#else
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
	s  = _mm256_add_ps(s,  _mm256_set1_ps(-1.0f));
	c  = _mm256_add_ps(c,  _mm256_set1_ps(-1.0f));
	sx = _mm256_add_ps(sx, _mm256_set1_ps(-1.0f));
	cx = _mm256_add_ps(cx, _mm256_set1_ps(-1.0f));
	s  = _mm256_mul_ps(s,  sx);
	c  = _mm256_mul_ps(c,  cx);
#endif
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
#if defined(__FMA__)
	s  = _mm_set1_ps(C3);
	c  = _mm_set1_ps(C3);
	s  = _mm_fmadd_ps(s, sx, _mm_set1_ps(C2));
	c  = _mm_fmadd_ps(c, cx, _mm_set1_ps(C2));
	s  = _mm_fmadd_ps(s, sx, _mm_set1_ps(C1));
	c  = _mm_fmadd_ps(c, cx, _mm_set1_ps(C1));
	s  = _mm_fmadd_ps(s, sx, _mm_set1_ps(-1.0f));
	c  = _mm_fmadd_ps(c, cx, _mm_set1_ps(-1.0f));
	s  = _mm_fmsub_ps(s, sx, s);
	c  = _mm_fmsub_ps(c, cx, c);
#else
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
	s  = _mm_add_ps(s,  _mm_set1_ps(-1.0f));
	c  = _mm_add_ps(c,  _mm_set1_ps(-1.0f));
	sx = _mm_add_ps(sx, _mm_set1_ps(-1.0f));
	cx = _mm_add_ps(cx, _mm_set1_ps(-1.0f));
	s  = _mm_mul_ps(s,  sx);
	c  = _mm_mul_ps(c,  cx);
#endif
	*Sin = s;
	*Cos = c;
}
#endif
/**************************************/

void Fourier_SinCos(float x, float *Sin, float *Cos) {
	float sx = x + -1.0f, cx = x; //! Phase-change sine to cosine
	sx *= sx; *Sin = (sx + -1.0f)*((C1 + (C2 + C3*sx)*sx)*sx + -1.0f);
	cx *= cx; *Cos = (cx + -1.0f)*((C1 + (C2 + C3*cx)*cx)*cx + -1.0f);
}

/**************************************/
//! EOF
/**************************************/
