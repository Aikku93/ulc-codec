/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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

//! Implementation notes for IMDCT:
//!  IMDCT is implemented via DCT-IV, which can be thought of
//!  as splitting the MDCT inputs into four regions:
//!   {A,B,C,D}
//!  and then taking the DCT-IV of:
//!   {C_r + D, B_r - A}
//!  On IMDCT, we get back these latter values following an
//!  inverse DCT-IV (which is itself a DCT-IV due to its
//!  involutive nature).
//!  From a prior call, we keep C_r+D buffered, which becomes
//!  A_r+B after accounting for movement to the next block.
//!  We can then state:
//!   Reverse(A_r + B) -        (B_r - A) = (A + B_r) - (B_r - A) = 2A
//!          (A_r + B) + Reverse(B_r - A) = (A_r + B) + (B - A_r) = 2B
//!           ^ Buffered         ^ New input data
//!  Allowing us to reconstruct the inputs A,B.
void Fourier_IMDCT(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, int N, int Overlap) {
	int i;
	const float *WinS = Fourier_SinTableN(Overlap);
	const float *WinC = WinS + Overlap;
	const float *Lap   = BufLap + N/2;
	const float *Tmp   = BufTmp + N/2;
	      float *OutLo = BufOut;
	      float *OutHi = BufOut + N;

	//! Undo transform
	for(i=0;i<N;i++) BufTmp[i] = BufIn[i];
	Fourier_DCT4(BufTmp, BufOut, N);

	//! Undo lapping
#if defined(__AVX__)
	__m256 a, b;
	__m256 t0, t1;
	__m256 c, s;

	for(i=0;i<(N-Overlap)/2;i+=8) {
		Lap -= 8; a = _mm256_load_ps(Lap);
		b = _mm256_load_ps(Tmp); Tmp += 8;
		a = _mm256_shuffle_ps(a, a, 0x1B);
		a = _mm256_permute2f128_ps(a, a, 0x01);
		b = _mm256_shuffle_ps(b, b, 0x1B);
		b = _mm256_permute2f128_ps(b, b, 0x01);
		_mm256_store_ps(OutLo, a); OutLo += 8;
		OutHi -= 8; _mm256_store_ps(OutHi, b);
	}
	for(;i<N/2;i+=8) {
		Lap -= 8; a = _mm256_load_ps(Lap);
		b = _mm256_load_ps(Tmp); Tmp += 8;
		a = _mm256_shuffle_ps(a, a, 0x1B);
		a = _mm256_permute2f128_ps(a, a, 0x01);
		WinC -= 8; c = _mm256_load_ps(WinC);
		s = _mm256_load_ps(WinS); WinS += 8;
		c = _mm256_shuffle_ps(c, c, 0x1B);
		c = _mm256_permute2f128_ps(c, c, 0x01);
#if defined(__FMA__)
		t0 = _mm256_mul_ps(s, b);
		t1 = _mm256_mul_ps(c, b);
		t0 = _mm256_fmsub_ps(c, a, t0);
		t1 = _mm256_fmadd_ps(s, a, t1);
#else
		t0 = _mm256_sub_ps(_mm256_mul_ps(c, a), _mm256_mul_ps(s, b));
		t1 = _mm256_add_ps(_mm256_mul_ps(s, a), _mm256_mul_ps(c, b));
#endif
		t1 = _mm256_shuffle_ps(t1, t1, 0x1B);
		t1 = _mm256_permute2f128_ps(t1, t1, 0x01);
		_mm256_store_ps(OutLo, t0); OutLo += 8;
		OutHi -= 8; _mm256_store_ps(OutHi, t1);
	}
#elif defined(__SSE__)
	__m128 a, b;
	__m128 t0, t1;
	__m128 c, s;

	for(i=0;i<(N-Overlap)/2;i+=4) {
		Lap -= 4; a = _mm_loadr_ps(Lap);
		b = _mm_load_ps(Tmp); Tmp += 4;
		_mm_store_ps(OutLo, a); OutLo += 4;
		OutHi -= 4; _mm_storer_ps(OutHi, b);
	}
	for(;i<N/2;i+=4) {
		Lap -= 4; a = _mm_loadr_ps(Lap);
		b = _mm_load_ps(Tmp); Tmp += 4;
		s = _mm_load_ps(WinS); WinS += 4;
		WinC -= 4; c = _mm_loadr_ps(WinC);
#if defined(__FMA__)
		t0 = _mm_mul_ps(s, b);
		t1 = _mm_mul_ps(c, b);
		t0 = _mm_fmsub_ps(c, a, t0);
		t1 = _mm_fmadd_ps(s, a, t1);
#else
		t0 = _mm_sub_ps(_mm_mul_ps(c, a), _mm_mul_ps(s, b));
		t1 = _mm_add_ps(_mm_mul_ps(s, a), _mm_mul_ps(c, b));
#endif
		_mm_store_ps(OutLo, t0); OutLo += 4;
		OutHi -= 4; _mm_storer_ps(OutHi, t1);
	}
#else
	for(i=0;i<(N-Overlap)/2;i++) {
		float a = *--Lap;
		float b = *Tmp++;
		*OutLo++ = a;
		*--OutHi = b;
	}
	for(;i<N/2;i++) {
		float a = *--Lap;
		float b = *Tmp++;
		float c = *--WinC;
		float s = *WinS++;
		*OutLo++ = c*a - s*b;
		*--OutHi = s*a + c*b;
	}
#endif
	//! Copy state to old block
	for(i=0;i<N/2;i++) BufLap[i] = BufTmp[i];
}

/**************************************/
//! EOF
/**************************************/
