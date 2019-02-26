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

7th order sine approximation
Estimates Sin[x*Pi/2], x=[0,1]
RMSE=0.00000128434 (PSNR=117.83dB)
------------------------------

a[x_] := Sin[x*Pi/2];
b[x_] := u*x + v*x^3 + w*x^5 + y*x^7;
Solve[{
  a[0] == b[0],
  a'[0] == b'[0],
  a[1] == b[1],
  a'[1] == b'[1],
  Integrate[a[x], {x, 0, 1}] == Integrate[b[x], {x, 0, 1}]
  }, {u, v, w, y}]

After a lot of wrangling, we arrive at:
  eq[x_] := (-x/(2*Pi)) * (-Pi^2 + x^2*(-96 + 13*Pi + 6*Pi^2) + x^2*(x^2*((192 - 33*Pi - 9*Pi^2) + x^2*(-96 + 18*Pi + 4*Pi^2))))
Pulling out constants:
  s0 = -Pi^2
  s1 = -96 + 13*Pi + 6*Pi^2
  s2 = 192 - 33*Pi - 9*Pi^2
  s3 = -96 + 18*Pi + 4*Pi^2
  s4 = -1/(2*Pi)
  eq[x_] := x*s4 * (s0 + x^2*(s1 + x^2*(s2 + x^2*s3)))

6th order cosine approximation
Estimates Cos[x*Pi/2], x=[0,1]
RMSE=0.0000127222 (PSNR=97.91dB)
------------------------------

a[x_] := Cos[x*Pi/2];
b[x_] := u + v*x^2 + w*x^4 + y*x^6;
Solve[{
  a[0] == b[0],
  a'[0] == b'[0],
  a[1] == b[1],
  a'[1] == b'[1],
  Integrate[a[x], {x, 0, 1}] == Integrate[b[x], {x, 0, 1}]
  }, {u, v, w, y}]

After a lot of wrangling, we arrive at:
  eq[x_] := (x^2-1) * (16*Pi + x^2((420 - 128*Pi - 3*Pi^2) - x^2(420 - 112*Pi - 7*Pi^2))) * (-1/(16*Pi))
Pulling out constants:
  c0 = 16*Pi
  c1 = 420 - 128*Pi - 3*Pi^2
  c2 = (-420 - 112*Pi - 7*Pi^2)
  c3 = -1/(16*Pi)
  eq[x_] := (x^2-1) * (c0 + x^2(c1 + c2*x^2)) * c3

!*/
/**************************************/

#define S0 (-0x1.3BD3CC9BE45DEp+3) //! -Pi^2
#define S1 ( 0x1.03BBB18A66B88p+2) //! -96 + 13*Pi + 6*Pi^2
#define S2 (-0x1.FEF91DC5D1B00p-2) //! 192 - 33*Pi - 9*Pi^2
#define S3 ( 0x1.BBC4457C4C800p-6) //! -96 + 18*Pi + 4*Pi^2
#define S4 (-0x1.45F306DC9C883p-3) //! -1/(2*Pi)

#define C0 ( 0x1.921FB54442D18p+5) //! 16*Pi
#define C1 (-0x1.77720E5C0749Ap+3) //! 420 - 128*Pi - 3*Pi^2
#define C2 ( 0x1.E426BBA8D7B00p-1) //! -(420 - 112*Pi - 7*Pi^2)
#define C3 (-0x1.45F306DC9C883p-6) //! -1/(16*Pi)

/**************************************/
#if defined(__AVX__)
void Fourier_SinCosAVX(__m256 x, __m256 *Sin, __m256 *Cos) {
	__m256 c, s;
	__m256 x2 = _mm256_mul_ps(x, x);
	s = _mm256_add_ps(_mm256_set1_ps(S2), _mm256_mul_ps(x2, _mm256_set1_ps(S3)));
	c = _mm256_add_ps(_mm256_set1_ps(C1), _mm256_mul_ps(x2, _mm256_set1_ps(C2)));
	s = _mm256_add_ps(_mm256_set1_ps(S1), _mm256_mul_ps(x2, s));
	c = _mm256_add_ps(_mm256_set1_ps(C0), _mm256_mul_ps(x2, c));
	s = _mm256_add_ps(_mm256_set1_ps(S0), _mm256_mul_ps(x2, s));
	c = _mm256_mul_ps(_mm256_set1_ps(C3), c);
	s = _mm256_mul_ps(_mm256_set1_ps(S4), s);
	c = _mm256_mul_ps(c, _mm256_sub_ps(x2, _mm256_set1_ps(1.0f)));
	s = _mm256_mul_ps(s, x);
	*Cos = c;
	*Sin = s;
}
#endif
/**************************************/
#if defined(__SSE__)
void Fourier_SinCosSSE(__m128 x, __m128 *Sin, __m128 *Cos) {
	__m128 c, s;
	__m128 x2 = _mm_mul_ps(x, x);
	s = _mm_add_ps(_mm_set1_ps(S2), _mm_mul_ps(x2, _mm_set1_ps(S3)));
	c = _mm_add_ps(_mm_set1_ps(C1), _mm_mul_ps(x2, _mm_set1_ps(C2)));
	s = _mm_add_ps(_mm_set1_ps(S1), _mm_mul_ps(x2, s));
	c = _mm_add_ps(_mm_set1_ps(C0), _mm_mul_ps(x2, c));
	s = _mm_add_ps(_mm_set1_ps(S0), _mm_mul_ps(x2, s));
	c = _mm_mul_ps(_mm_set1_ps(C3), c);
	s = _mm_mul_ps(_mm_set1_ps(S4), s);
	c = _mm_mul_ps(c, _mm_sub_ps(x2, _mm_set1_ps(1.0f)));
	s = _mm_mul_ps(s, x);
	*Cos = c;
	*Sin = s;
}
#endif
/**************************************/

void Fourier_SinCos(float x, float *Sin, float *Cos) {
	float c, s;
	float x2 = x*x;
	s  = S2 + S3*x2;
	c  = C1 + C2*x2;
	s  = S1 +  s*x2;
	c  = C0 +  c*x2;
	s  = S0 +  s*x2;
	c *= C3;
	s *= S4;
	c *= x2 - 1.0f;
	s *= x;
	*Cos = c;
	*Sin = s;
}

/**************************************/
//! EOF
/**************************************/
