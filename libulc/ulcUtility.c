/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#if defined(__AVX__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include <stddef.h>
/**************************************/

//! Apply pre-echo reduction formula
void ULC_Transform_AntiPreEcho(float *Buf, size_t N) {
	size_t i;
#if defined(__AVX__) //! 16x256bit, but only use 8x for smaller loop code
	__m256 a, va;
	__m256 b, vb;
	__m256 c, vc;
	__m256 d, vd;
	for(i=0;i<N;i+=32) {
		 a = _mm256_load_ps(Buf + i+ 0);
		 b = _mm256_load_ps(Buf + i+ 8);
		 c = _mm256_load_ps(Buf + i+16);
		 d = _mm256_load_ps(Buf + i+24);
		va = _mm256_permute2f128_ps(a, b, 0x20);
		vb = _mm256_permute2f128_ps(a, b, 0x31);
		vc = _mm256_permute2f128_ps(c, d, 0x20);
		vd = _mm256_permute2f128_ps(c, d, 0x31);
		 a = _mm256_shuffle_ps(va, vb, 0x88);
		 b = _mm256_shuffle_ps(va, vb, 0xDD);
		 c = _mm256_shuffle_ps(vc, vd, 0x88);
		 d = _mm256_shuffle_ps(vc, vd, 0xDD);
		va = _mm256_add_ps(a, b);
		vb = _mm256_sub_ps(a, b);
		vc = _mm256_add_ps(c, d);
		vd = _mm256_sub_ps(c, d);
		 a = _mm256_unpacklo_ps(va, vb);
		 b = _mm256_unpackhi_ps(va, vb);
		 c = _mm256_unpacklo_ps(vc, vd);
		 d = _mm256_unpackhi_ps(vc, vd);
		va = _mm256_permute2f128_ps(a, b, 0x20);
		vb = _mm256_permute2f128_ps(a, b, 0x31);
		vc = _mm256_permute2f128_ps(c, d, 0x20);
		vd = _mm256_permute2f128_ps(c, d, 0x31);
		_mm256_store_ps(Buf + i+ 0, va);
		_mm256_store_ps(Buf + i+ 8, vb);
		_mm256_store_ps(Buf + i+16, vc);
		_mm256_store_ps(Buf + i+24, vd);
	}
#elif defined(__SSE__) //! 8x128bit
	__m128 a, va;
	__m128 b, vb;
	__m128 c, vc;
	__m128 d, vd;
	for(i=0;i<N;i+=16) {
		 a = _mm_load_ps(Buf + i+ 0);
		 b = _mm_load_ps(Buf + i+ 4);
		 c = _mm_load_ps(Buf + i+ 8);
		 d = _mm_load_ps(Buf + i+12);
		va = _mm_shuffle_ps(a, b, 0x88);
		vb = _mm_shuffle_ps(a, b, 0xDD);
		vc = _mm_shuffle_ps(c, d, 0x88);
		vd = _mm_shuffle_ps(c, d, 0xDD);
		 a = _mm_add_ps(va, vb);
		 b = _mm_sub_ps(va, vb);
		 c = _mm_add_ps(vc, vd);
		 d = _mm_sub_ps(vc, vd);
		va = _mm_unpacklo_ps(a, b);
		vb = _mm_unpackhi_ps(a, b);
		vc = _mm_unpacklo_ps(c, d);
		vd = _mm_unpackhi_ps(c, d);
		_mm_store_ps(Buf + i+ 0, va);
		_mm_store_ps(Buf + i+ 4, vb);
		_mm_store_ps(Buf + i+ 8, vc);
		_mm_store_ps(Buf + i+12, vd);
	}
#else
	for(i=0;i<N;i+=2) {
		float *a = &Buf[i+0], va = *a;
		float *b = &Buf[i+1], vb = *b;
		*a = va + vb;
		*b = va - vb;
	}
#endif
}

/**************************************/
//! EOF
/**************************************/
