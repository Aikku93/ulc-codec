/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2022, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#if defined(__AVX__) || defined(__FMA__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#define FOURIER_FORCED_INLINE static inline __attribute__((always_inline))
#define FOURIER_ASSUME(Cond) (Cond) ? ((void)0) : __builtin_unreachable()
#define FOURIER_ASSUME_ALIGNED(x,Align) x = __builtin_assume_aligned(x,Align)
/**************************************/
#if defined(__AVX__)
  typedef __m256 Fourier_Vec_t;
# define FOURIER_VSTRIDE            8
# define FOURIER_VLOAD(Src)         _mm256_load_ps(Src)
# define FOURIER_VLOADU(Src)        _mm256_loadu_ps(Src)
# define FOURIER_VSTORE(Dst, x)     _mm256_store_ps(Dst, x)
# define FOURIER_VSTOREU(Dst, x)    _mm256_storeu_ps(Dst, x)
# define FOURIER_VSET1(x)           _mm256_set1_ps(x)
# define FOURIER_VSET_LINEAR_RAMP() _mm256_setr_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f)
# define FOURIER_VADD(x, y)         _mm256_add_ps(x, y)
# define FOURIER_VSUB(x, y)         _mm256_sub_ps(x, y)
# define FOURIER_VMUL(x, y)         _mm256_mul_ps(x, y)
# define FOURIER_VREVERSE_LANE(x)   _mm256_shuffle_ps(x, x, 0x1B)
# define FOURIER_VREVERSE(x)        _mm256_permute2f128_ps(FOURIER_VREVERSE_LANE(x), FOURIER_VREVERSE_LANE(x), 0x01)
# define FOURIER_VNEGATE_ODD(x)     _mm256_xor_ps(x, _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f))
# if defined(__FMA__)
#  define FOURIER_VFMA(x, y, a)     _mm256_fmadd_ps(x, y, a)
#  define FOURIER_VFMS(x, y, a)     _mm256_fmsub_ps(x, y, a)
#  define FOURIER_VNFMA(x, y, a)    _mm256_fnmadd_ps(x, y, a)
# else
#  define FOURIER_VFMA(x, y, a)     _mm256_add_ps(_mm256_mul_ps(x, y), a)
#  define FOURIER_VFMS(x, y, a)     _mm256_sub_ps(_mm256_mul_ps(x, y), a)
#  define FOURIER_VNFMA(x, y, a)    _mm256_sub_ps(a, _mm256_mul_ps(x, y))
# endif
  FOURIER_FORCED_INLINE void FOURIER_VINTERLEAVE(Fourier_Vec_t a, Fourier_Vec_t b, Fourier_Vec_t *Lo, Fourier_Vec_t *Hi) {
	__m256 t0 = _mm256_unpacklo_ps(a, b);
	__m256 t1 = _mm256_unpackhi_ps(a, b);
	*Lo = _mm256_permute2f128_ps(t0, t1, 0x20);
	*Hi = _mm256_permute2f128_ps(t0, t1, 0x31);
  }
  FOURIER_FORCED_INLINE void FOURIER_VSPLIT_EVEN_ODD(Fourier_Vec_t l, Fourier_Vec_t h, Fourier_Vec_t *Even, Fourier_Vec_t *Odd) {
	__m256 a = _mm256_permute2f128_ps(l, h, 0x20);
	__m256 b = _mm256_permute2f128_ps(l, h, 0x31);
	*Even = _mm256_shuffle_ps(a, b, 0x88);
	*Odd  = _mm256_shuffle_ps(a, b, 0xDD);
  }
  FOURIER_FORCED_INLINE void FOURIER_VSPLIT_EVEN_ODDREV(Fourier_Vec_t l, Fourier_Vec_t h, Fourier_Vec_t *Even, Fourier_Vec_t *Odd) {
	__m256 a = _mm256_permute2f128_ps(l, h, 0x20);
	__m256 b = _mm256_permute2f128_ps(l, h, 0x31);
	*Even = _mm256_shuffle_ps(a, b, 0x88);
	b     = _mm256_shuffle_ps(b, a, 0x77);
	*Odd  = _mm256_permute2f128_ps(b, b, 0x01);
  }
#elif defined(__SSE__)
  typedef __m128 Fourier_Vec_t;
# define FOURIER_VSTRIDE            4
# define FOURIER_VLOAD(Src)         _mm_load_ps(Src)
# define FOURIER_VLOADU(Src)        _mm_loadu_ps(Src)
# define FOURIER_VSTORE(Dst, x)     _mm_store_ps(Dst, x)
# define FOURIER_VSTOREU(Dst, x)    _mm_storeu_ps(Dst, x)
# define FOURIER_VSET1(x)           _mm_set1_ps(x)
# define FOURIER_VSET_LINEAR_RAMP() _mm_setr_ps(0.0f, 1.0f, 2.0f, 3.0f)
# define FOURIER_VADD(x, y)         _mm_add_ps(x, y)
# define FOURIER_VSUB(x, y)         _mm_sub_ps(x, y)
# define FOURIER_VMUL(x, y)         _mm_mul_ps(x, y)
# define FOURIER_VREVERSE(x)        _mm_shuffle_ps(x, x, 0x1B)
# define FOURIER_VNEGATE_ODD(x)     _mm_xor_ps(x, _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f))
# if defined(__FMA__)
#  define FOURIER_VFMA(x, y, a)     _mm_fmadd_ps(x, y, a)
#  define FOURIER_VFMS(x, y, a)     _mm_fmsub_ps(x, y, a)
#  define FOURIER_VNFMA(x, y, a)    _mm_fnmadd_ps(x, y, a)
# else
#  define FOURIER_VFMA(x, y, a)     _mm_add_ps(_mm_mul_ps(x, y), a)
#  define FOURIER_VFMS(x, y, a)     _mm_sub_ps(_mm_mul_ps(x, y), a)
#  define FOURIER_VNFMA(x, y, a)    _mm_sub_ps(a, _mm_mul_ps(x, y))
# endif
  FOURIER_FORCED_INLINE void FOURIER_VINTERLEAVE(Fourier_Vec_t a, Fourier_Vec_t b, Fourier_Vec_t *Lo, Fourier_Vec_t *Hi) {
	*Lo = _mm_unpacklo_ps(a, b);
	*Hi = _mm_unpackhi_ps(a, b);
  }
  FOURIER_FORCED_INLINE void FOURIER_VSPLIT_EVEN_ODD(Fourier_Vec_t l, Fourier_Vec_t h, Fourier_Vec_t *Even, Fourier_Vec_t *Odd) {
	*Even = _mm_shuffle_ps(l, h, 0x88);
	*Odd  = _mm_shuffle_ps(l, h, 0xDD);
  }
  FOURIER_FORCED_INLINE void FOURIER_VSPLIT_EVEN_ODDREV(Fourier_Vec_t l, Fourier_Vec_t h, Fourier_Vec_t *Even, Fourier_Vec_t *Odd) {
	*Even = _mm_shuffle_ps(l, h, 0x88);
	*Odd  = _mm_shuffle_ps(h, l, 0x77);
  }
#else
  typedef float Fourier_Vec_t;
# define FOURIER_VSTRIDE        1
# define FOURIER_VSET1(x)      ((x))
# define FOURIER_VADD(x, y)    ((x) + (y))
# define FOURIER_VSUB(x, y)    ((x) - (y))
# define FOURIER_VMUL(x, y)    ((x) * (y))
# define FOURIER_VFMA(x, y, a) ((x) * (y) + (a))
#endif
/**************************************/

//! Sin[z] approximation (where z = x*Pi/2, x < Pi/2)
//! Coefficients derived by RMSE minimization
FOURIER_FORCED_INLINE
Fourier_Vec_t Fourier_Sin(Fourier_Vec_t x) {
	Fourier_Vec_t x2  = FOURIER_VMUL(x, x);
	Fourier_Vec_t Res = FOURIER_VFMA(x2, FOURIER_VSET1(+0x1.3C8B08p-13f), FOURIER_VSET1(-0x1.3237ECp-8f));
	              Res = FOURIER_VFMA(x2, Res, FOURIER_VSET1(+0x1.4667ACp-4f));
	              Res = FOURIER_VFMA(x2, Res, FOURIER_VSET1(-0x1.4ABBB8p-1f));
	              Res = FOURIER_VFMA(x2, Res, FOURIER_VSET1(+0x1.921FB4p0f));
	return FOURIER_VMUL(x, Res);
}

/**************************************/

//! Cos[z] approximation (where z = x*Pi/2, x < Pi/2)
//! Coefficients derived by RMSE minimization
FOURIER_FORCED_INLINE
Fourier_Vec_t Fourier_Cos(Fourier_Vec_t x) {
	Fourier_Vec_t x2  = FOURIER_VMUL(x, x);
	Fourier_Vec_t Res = FOURIER_VFMA(x2, FOURIER_VSET1(+0x1.C2EA62p-11f), FOURIER_VSET1(-0x1.550810p-6f));
	              Res = FOURIER_VFMA(x2, Res, FOURIER_VSET1(+0x1.03BDD6p-2f));
	              Res = FOURIER_VFMA(x2, Res, FOURIER_VSET1(-0x1.3BD3B2p0f));
	              Res = FOURIER_VFMA(x2, Res, FOURIER_VSET1(1.0f));
	return Res;
}

/**************************************/
//! EOF
/**************************************/
