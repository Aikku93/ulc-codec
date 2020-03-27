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
/*!

8th order cosine approximation
Estimates Cos[x*Pi/2] (with x=[-1,1]) using Legendre polynomials,
using the polynomial form of:
 y = c1 + (c2 + (c3 + (c4 + c5*x^2)*x^2)*x^2)*x^2;

NOTE: This has bad behaviour at x={0,-1,+1} but these values are
unused in all calls to these functions, so this is still good.

!*/
/**************************************/

#define C1 ( 0x1.FFFFFFp-1f)
#define C2 (-0x1.3BD3AEp0f)
#define C3 ( 0x1.03BDD4p-2f)
#define C4 (-0x1.550D85p-6f)
#define C5 ( 0x1.C390DDp-11f)

/**************************************/
#if defined(__AVX__)
void Fourier_SinCosAVX(__m256 x, __m256 *Sin, __m256 *Cos) {
	__m256 sx = _mm256_sub_ps(x, _mm256_set1_ps(1.0f));
	__m256 cx = x;
	__m256 s, c;
	sx = _mm256_mul_ps(sx, sx);
	cx = _mm256_mul_ps(cx, cx);
#if defined(__FMA__)
	s  = _mm256_set1_ps(C5);
	c  = _mm256_set1_ps(C5);
	s  = _mm256_fmadd_ps(s, sx, _mm256_set1_ps(C4));
	c  = _mm256_fmadd_ps(c, cx, _mm256_set1_ps(C4));
	s  = _mm256_fmadd_ps(s, sx, _mm256_set1_ps(C3));
	c  = _mm256_fmadd_ps(c, cx, _mm256_set1_ps(C3));
	s  = _mm256_fmadd_ps(s, sx, _mm256_set1_ps(C2));
	c  = _mm256_fmadd_ps(c, cx, _mm256_set1_ps(C2));
	s  = _mm256_fmadd_ps(s, sx, _mm256_set1_ps(C1));
	c  = _mm256_fmadd_ps(c, cx, _mm256_set1_ps(C1));
#else
	s  = _mm256_mul_ps(sx, _mm256_set1_ps(C5));
	c  = _mm256_mul_ps(cx, _mm256_set1_ps(C5));
	s  = _mm256_add_ps(s,  _mm256_set1_ps(C4));
	c  = _mm256_add_ps(c,  _mm256_set1_ps(C4));
	s  = _mm256_mul_ps(s,  sx);
	c  = _mm256_mul_ps(c,  cx);
	s  = _mm256_add_ps(s,  _mm256_set1_ps(C3));
	c  = _mm256_add_ps(c,  _mm256_set1_ps(C3));
	s  = _mm256_mul_ps(s,  sx);
	c  = _mm256_mul_ps(c,  cx);
	s  = _mm256_add_ps(s,  _mm256_set1_ps(C2));
	c  = _mm256_add_ps(c,  _mm256_set1_ps(C2));
	s  = _mm256_mul_ps(s,  sx);
	c  = _mm256_mul_ps(c,  cx);
	s  = _mm256_add_ps(s,  _mm256_set1_ps(C1));
	c  = _mm256_add_ps(c,  _mm256_set1_ps(C1));
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
	s  = _mm_set1_ps(C5);
	c  = _mm_set1_ps(C5);
	s  = _mm_fmadd_ps(s, sx, _mm_set1_ps(C4));
	c  = _mm_fmadd_ps(c, cx, _mm_set1_ps(C4));
	s  = _mm_fmadd_ps(s, sx, _mm_set1_ps(C3));
	c  = _mm_fmadd_ps(c, cx, _mm_set1_ps(C3));
	s  = _mm_fmadd_ps(s, sx, _mm_set1_ps(C2));
	c  = _mm_fmadd_ps(c, cx, _mm_set1_ps(C2));
	s  = _mm_fmadd_ps(s, sx, _mm_set1_ps(C1));
	c  = _mm_fmadd_ps(c, cx, _mm_set1_ps(C1));
#else
	s  = _mm_mul_ps(sx, _mm_set1_ps(C5));
	c  = _mm_mul_ps(cx, _mm_set1_ps(C5));
	s  = _mm_add_ps(s,  _mm_set1_ps(C4));
	c  = _mm_add_ps(c,  _mm_set1_ps(C4));
	s  = _mm_mul_ps(s,  sx);
	c  = _mm_mul_ps(c,  cx);
	s  = _mm_add_ps(s,  _mm_set1_ps(C3));
	c  = _mm_add_ps(c,  _mm_set1_ps(C3));
	s  = _mm_mul_ps(s,  sx);
	c  = _mm_mul_ps(c,  cx);
	s  = _mm_add_ps(s,  _mm_set1_ps(C2));
	c  = _mm_add_ps(c,  _mm_set1_ps(C2));
	s  = _mm_mul_ps(s,  sx);
	c  = _mm_mul_ps(c,  cx);
	s  = _mm_add_ps(s,  _mm_set1_ps(C1));
	c  = _mm_add_ps(c,  _mm_set1_ps(C1));
#endif
	*Sin = s;
	*Cos = c;
}
#endif
/**************************************/

void Fourier_SinCos(float x, float *Sin, float *Cos) {
	float sx = x + -1.0f, cx = x; //! Phase-change sine to cosine
	sx *= sx; *Sin = C1 + (C2 + (C3 + (C4 + C5*sx)*sx)*sx)*sx;
	cx *= cx; *Cos = C1 + (C2 + (C3 + (C4 + C5*cx)*cx)*cx)*cx;
}

/**************************************/
//! EOF
/**************************************/
