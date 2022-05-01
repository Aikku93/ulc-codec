/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2022, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#if defined(__AVX__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include "Fourier.h"
#include "FourierHelper.h"
/**************************************/

//! Implementation notes for MDCT:
//!  MDCT is implemented via DCT-IV, which can be thought of
//!  as splitting the MDCT inputs into four regions:
//!   {A,B,C,D}
//!  and then taking the DCT-IV of:
//!   {C_r + D, B_r - A}
//!  During non-overlap regions, A,B scales by 0.0 and C,D by 1.0.
//!  As an optimization, the `C_r + D` region is written to
//!  in a backwards direction, leading to `D_r + C`; this
//!  leads to reuse of the trigonometric constants, avoiding
//!  reading the array of sin/cos twice per call.
//! Pseudocode:
//!   for(n=0;n<(N-Overlap)/2;n++) {
//!     MDCT[N/2-1-n] =  New[n];
//!     MDST[N/2-1-n] = -New[n];
//!     MDCT[N/2+n]   =  Old[N-1-n];
//!     MDST[N/2+n]   = -Old[N-1-n];
//!   }
//!   for(;n<N/2;n++) {
//!     MDCT[N/2-1-n] =  c*New[n] + s*New[N-1-n];
//!     MDST[N/2-1-n] = -c*New[n] + s*New[N-1-n];
//!     MDCT[N/2+n]   = -s*Old[n] + c*Old[N-1-n];
//!     MDST[N/2+n]   = -s*Old[n] - c*Old[N-1-n];
//!   }
//!   DCT4(MDCT);
//!   DCT4(MDST);
//!   for(n=0;n<N/2;n++) Swap(MDST[n], MDST[N-1-n])
//! If we don't care about MDST, this can be refactored:
//!   for(n=0;n<(N-Overlap)/2;n++) {
//!     MDCT[N/2-1-n] = New[n];
//!     MDCT[N/2+n]   = Old[n];
//!     Old[n]        = New[N-1-n];
//!   }
//!   for(;n<N/2;n++) {
//!     MDCT[N/2-1-n] =  c*New[n] + s*New[N-1-n];
//!     MDCT[N/2+n]   = Old[n];
//!     Old[n]        = -s*New[n] + c*New[N-1-n];
//!   }
//!   DCT4(MDCT);
void Fourier_MDCT_MDST(float *MDCT, float *MDST, const float *New, float *Lap, float *BufTmp, int N, int Overlap, const float *ModulationWindow) {
	int n;
	FOURIER_ASSUME_ALIGNED(MDCT,   32);
	FOURIER_ASSUME_ALIGNED(MDST,   32);
	FOURIER_ASSUME_ALIGNED(New,    32);
	FOURIER_ASSUME_ALIGNED(Lap,    32);
	FOURIER_ASSUME_ALIGNED(BufTmp, 32);
	FOURIER_ASSUME_ALIGNED(ModulationWindow, 32);
	FOURIER_ASSUME(N >= 16 && N <= 8192);
	FOURIER_ASSUME(Overlap >= 0 && Overlap <= N);

	const float *Win   = Fourier_SinTableN(Overlap, ModulationWindow);
	      float *LapLo = Lap;
	      float *LapHi = Lap + N;
	const float *NewLo = New;
	const float *NewHi = New + N;
	      float *MDCTMid = MDCT + N/2;
	      float *MDSTMid = MDST + N/2;

	//! Perform windowed lapping
	//! NOTE: We compute MDST using DCT-IV, so sign-flip every other
	//! value here, and then reverse the whole array after DCT-IV.
	{
#if defined(__AVX__)
		__m256 c, s;
		__m256 A, Br, C, Dr;
		__m256 Zero = _mm256_setzero_ps();
		__m256 XorMask = _mm256_setr_ps(0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f, 0.0f, -0.0f);
		for(n=0;n<(N-Overlap)/2;n+=8) {
			A  = _mm256_load_ps(LapLo   + n);
			Br = _mm256_load_ps(LapHi-8 - n);
			C  = _mm256_load_ps(NewLo   + n);
			Dr = _mm256_load_ps(NewHi-8 - n);
			Br = _mm256_permute2f128_ps(Br, Br, 0x01);
			Br = _mm256_shuffle_ps     (Br, Br, 0x1B);
			C  = _mm256_permute2f128_ps(C, C, 0x01);
			C  = _mm256_shuffle_ps     (C, C, 0x1B);
			_mm256_store_ps(LapLo     + n, Zero);
			_mm256_store_ps(LapHi  -8 - n, Dr);
			_mm256_store_ps(MDCTMid-8 - n, C);
			_mm256_store_ps(MDSTMid-8 - n, _mm256_xor_ps(C, XorMask));
			c = _mm256_sub_ps(Br, A);
			s = _mm256_add_ps(Br, A);
			_mm256_store_ps(MDCTMid   + n, c);
			_mm256_store_ps(MDSTMid   + n, _mm256_xor_ps(s, XorMask));
		}
		for(;n<N/2;n+=8) {
			c  = _mm256_load_ps(Win); Win += 8;
			s  = _mm256_load_ps(Win); Win += 8;
			A  = _mm256_permute2f128_ps(c, s, 0x20);
			C  = _mm256_permute2f128_ps(c, s, 0x31);
			c  = _mm256_shuffle_ps(A, C, 0x88);
			s  = _mm256_shuffle_ps(A, C, 0xDD);
			A  = _mm256_load_ps(LapLo   + n);
			Br = _mm256_load_ps(LapHi-8 - n);
			C  = _mm256_load_ps(NewLo   + n);
			Dr = _mm256_load_ps(NewHi-8 - n);
			Zero = _mm256_permute2f128_ps(c, c, 0x01); //! <- Reverse for Dr*c :/
			Zero = _mm256_shuffle_ps(Zero, Zero, 0x1B);
			_mm256_store_ps(LapLo   + n, _mm256_mul_ps(s, C));
			_mm256_store_ps(LapHi-8 - n, _mm256_mul_ps(Zero, Dr));
			Br = _mm256_permute2f128_ps(Br, Br, 0x01);
			Br = _mm256_shuffle_ps     (Br, Br, 0x1B);
			Dr = _mm256_permute2f128_ps(Dr, Dr, 0x01);
			Dr = _mm256_shuffle_ps     (Dr, Dr, 0x1B);
			C  = _mm256_mul_ps(C,  c);
			Dr = _mm256_mul_ps(Dr, s);
			c  = _mm256_add_ps(C,  Dr);
			s  = _mm256_sub_ps(C,  Dr);
			c  = _mm256_permute2f128_ps(c, c, 0x01);
			c  = _mm256_shuffle_ps     (c, c, 0x1B);
			s  = _mm256_permute2f128_ps(s, s, 0x01);
			s  = _mm256_shuffle_ps     (s, s, 0x1B);
			_mm256_store_ps(MDCTMid-8 - n, c);
			_mm256_store_ps(MDSTMid-8 - n, _mm256_xor_ps(s, XorMask));
			c  = _mm256_sub_ps(Br, A);
			s  = _mm256_add_ps(Br, A);
			_mm256_store_ps(MDCTMid   + n, c);
			_mm256_store_ps(MDSTMid   + n, _mm256_xor_ps(s, XorMask));
		}
#elif defined(__SSE__)
		__m128 c, s;
		__m128 A, Br, C, Dr;
		__m128 Zero = _mm_setzero_ps();
		__m128 XorMask = _mm_setr_ps(0.0f, -0.0f, 0.0f, -0.0f);
		for(n=0;n<(N-Overlap)/2;n+=4) {
			A  = _mm_load_ps(LapLo   + n);
			Br = _mm_load_ps(LapHi-4 - n);
			C  = _mm_load_ps(NewLo   + n);
			Dr = _mm_load_ps(NewHi-4 - n);
			Br = _mm_shuffle_ps(Br, Br, 0x1B);
			C  = _mm_shuffle_ps(C,  C,  0x1B);
			_mm_store_ps(LapLo     + n, Zero);
			_mm_store_ps(LapHi  -4 - n, Dr);
			_mm_store_ps(MDCTMid-4 - n, C);
			_mm_store_ps(MDSTMid-4 - n, _mm_xor_ps(C, XorMask));
			c = _mm_sub_ps(Br, A);
			s = _mm_add_ps(Br, A);
			_mm_store_ps(MDCTMid   + n, c);
			_mm_store_ps(MDSTMid   + n, _mm_xor_ps(s, XorMask));
		}
		for(;n<N/2;n+=4) {
			A  = _mm_load_ps(Win); Win += 4;
			C  = _mm_load_ps(Win); Win += 4;
			c  = _mm_shuffle_ps(C, A, 0x22); //! <- Reversed
			s  = _mm_shuffle_ps(A, C, 0xDD);
			A  = _mm_load_ps(LapLo   + n);
			Br = _mm_load_ps(LapHi-4 - n);
			C  = _mm_load_ps(NewLo   + n);
			Dr = _mm_load_ps(NewHi-4 - n);
			_mm_store_ps(LapLo   + n, _mm_mul_ps(s, C));
			_mm_store_ps(LapHi-4 - n, _mm_mul_ps(c, Dr));
			c  = _mm_shuffle_ps(c,  c,  0x1B); //! <- Restore
			Br = _mm_shuffle_ps(Br, Br, 0x1B);
			Dr = _mm_shuffle_ps(Dr, Dr, 0x1B);
			C  = _mm_mul_ps(C,  c);
			Dr = _mm_mul_ps(Dr, s);
			c  = _mm_add_ps(C, Dr);
			s  = _mm_sub_ps(C, Dr);
			c  = _mm_shuffle_ps(c, c, 0x1B);
			s  = _mm_shuffle_ps(s, s, 0x1B);
			_mm_store_ps(MDCTMid-4 - n, c);
			_mm_store_ps(MDSTMid-4 - n, _mm_xor_ps(s, XorMask));
			c  = _mm_sub_ps(Br, A);
			s  = _mm_add_ps(Br, A);
			_mm_store_ps (MDCTMid   + n, c);
			_mm_store_ps (MDSTMid   + n, _mm_xor_ps(s, XorMask));
		}
#else
		float c, s;
		float A, Br, C, Dr;
		for(n=0;n<(N-Overlap)/2;n++) {
			A  = LapLo[   n];
			Br = LapHi[-1-n];
			C  = NewLo[   n];
			Dr = NewHi[-1-n];
			LapLo  [   n] =  0.0f;
			LapHi  [-1-n] =  Dr;
			MDCTMid[-1-n] =  C;
			MDSTMid[-1-n] =  C;
			MDCTMid[   n] = -A + Br;
			MDSTMid[   n] =  A + Br;
		}
		for(;n<N/2;n++) {
			c = *Win++;
			s = *Win++;
			A  = LapLo[   n];
			Br = LapHi[-1-n];
			C  = NewLo[   n];
			Dr = NewHi[-1-n];
			LapLo  [   n] = s*C;
			LapHi  [-1-n] = c*Dr;
			C *= c, Dr *= s;
			MDCTMid[-1-n] =  C + Dr;
			MDSTMid[-1-n] =  C - Dr;
			MDCTMid[   n] = -A + Br;
			MDSTMid[   n] =  A + Br;

		}
		for(n=1;n<N;n+=2) MDST[n] = -MDST[n];
#endif
	}

	//! Do actual transforms
	Fourier_DCT4T(MDCT, BufTmp, N);
	Fourier_DCT4T(MDST, BufTmp, N);

	//! Reverse array for MDST
	{
		float *BufLo = MDST;
		float *BufHi = MDST + N;
#if defined(__AVX_)
		__m256 v0, v1;
		for(n=0;n<N/2;n+=8) {
			BufHi -= 8;
			v0 = _mm256_load_ps(BufLo);
			v1 = _mm256_load_ps(BufHi);
			v0 = _mm256_shuffle_ps(v0, v0, 0x1B);
			v1 = _mm256_shuffle_ps(v1, v1, 0x1B);
			v0 = _mm256_permute2f128_ps(v0, v0, 0x01);
			v1 = _mm256_permute2f128_ps(v1, v1, 0x01);
			_mm256_store_ps(BufHi, v0);
			_mm256_store_ps(BufLo, v1);
			BufLo += 8;
		}
#elif defined(__SSE__)
		__m128 v0, v1;
		for(n=0;n<N/2;n+=4) {
			BufHi -= 4;
			v0 = _mm_load_ps(BufLo);
			v1 = _mm_load_ps(BufHi);
			_mm_storer_ps(BufHi, v0);
			_mm_storer_ps(BufLo, v1);
			BufLo += 4;
		}
#else
		for(n=0;n<N/2;n++) {
			BufHi--;
			float t = *BufLo;
			*BufLo = *BufHi;
			*BufHi = t;
			BufLo++;
		}
#endif
	}
}
void Fourier_MDCT(float *MDCT, const float *New, float *Lap, float *BufTmp, int N, int Overlap, const float *ModulationWindow) {
	int n;
	FOURIER_ASSUME_ALIGNED(MDCT,   32);
	FOURIER_ASSUME_ALIGNED(New,    32);
	FOURIER_ASSUME_ALIGNED(Lap,    32);
	FOURIER_ASSUME_ALIGNED(BufTmp, 32);
	FOURIER_ASSUME_ALIGNED(ModulationWindow, 32);
	FOURIER_ASSUME(N >= 16 && N <= 8192);
	FOURIER_ASSUME(Overlap >= 0 && Overlap <= N);

	const float *Win   = Fourier_SinTableN(Overlap, ModulationWindow);
	const float *NewLo = New;
	const float *NewHi = New + N;
	      float *MDCTMid = MDCT + N/2;

	//! Copy state from old block
	for(n=0;n<N/2;n++) MDCTMid[n] = Lap[n];

	//! Perform windowed lapping
	{
#if defined(__AVX__)
		__m256 a, b;
		__m256 t0, t1;
		__m256 c, s;

		for(n=0;n<(N-Overlap)/2;n+=8) {
			NewHi -= 8; b = _mm256_load_ps(NewHi);
			a = _mm256_load_ps(NewLo); NewLo += 8;
			b = _mm256_shuffle_ps(b, b, 0x1B);
			b = _mm256_permute2f128_ps(b, b, 0x01);
			a = _mm256_shuffle_ps(a, a, 0x1B);
			a = _mm256_permute2f128_ps(a, a, 0x01);
			_mm256_store_ps(Lap, b); Lap += 8;
			MDCTMid -= 8; _mm256_store_ps(MDCTMid, a);
		}
		for(;n<N/2;n+=8) {
			NewHi -= 8; b = _mm256_load_ps(NewHi);
			a = _mm256_load_ps(NewLo); NewLo += 8;
			b = _mm256_shuffle_ps(b, b, 0x1B);
			b = _mm256_permute2f128_ps(b, b, 0x01);
			c  = _mm256_load_ps(Win); Win += 8;
			s  = _mm256_load_ps(Win); Win += 8;
			t0 = _mm256_permute2f128_ps(c, s, 0x20);
			t1 = _mm256_permute2f128_ps(c, s, 0x31);
			c  = _mm256_shuffle_ps(t0, t1, 0x88);
			s  = _mm256_shuffle_ps(t0, t1, 0xDD);
#if defined(__FMA__)
			t0 = _mm256_mul_ps(s, b);
			t1 = _mm256_mul_ps(s, a);
			t0 = _mm256_fmadd_ps(c, a, t0);
			t1 = _mm256_fmsub_ps(c, b, t1);
#else
			t0 = _mm256_add_ps(_mm256_mul_ps(c, a), _mm256_mul_ps(s, b));
			t1 = _mm256_sub_ps(_mm256_mul_ps(c, b), _mm256_mul_ps(s, a));
#endif
			t0 = _mm256_shuffle_ps(t0, t0, 0x1B);
			t0 = _mm256_permute2f128_ps(t0, t0, 0x01);
			MDCTMid -= 8; _mm256_store_ps(MDCTMid, t0);
			_mm256_store_ps(Lap, t1); Lap += 8;
		}
#elif defined(__SSE__)
		__m128 a, b;
		__m128 t0, t1;
		__m128 c, s;

		for(n=0;n<(N-Overlap)/2;n+=4) {
			NewHi -= 4; b = _mm_loadr_ps(NewHi);
			a = _mm_load_ps(NewLo); NewLo += 4;
			_mm_store_ps(Lap, b); Lap += 4;
			MDCTMid -= 4; _mm_storer_ps(MDCTMid, a);
		}
		for(;n<N/2;n+=4) {
			NewHi -= 4; b = _mm_loadr_ps(NewHi);
			a = _mm_load_ps(NewLo); NewLo += 4;
			t0 = _mm_load_ps(Win); Win += 4;
			t1 = _mm_load_ps(Win); Win += 4;
			c  = _mm_shuffle_ps(t0, t1, 0x88);
			s  = _mm_shuffle_ps(t0, t1, 0xDD);
#if defined(__FMA__)
			t0 = _mm_mul_ps(s, b);
			t1 = _mm_mul_ps(s, a);
			t0 = _mm_fmadd_ps(c, a, t0);
			t1 = _mm_fmsub_ps(c, b, t1);
#else
			t0 = _mm_add_ps(_mm_mul_ps(c, a), _mm_mul_ps(s, b));
			t1 = _mm_sub_ps(_mm_mul_ps(c, b), _mm_mul_ps(s, a));
#endif
			MDCTMid -= 4; _mm_storer_ps(MDCTMid, t0);
			_mm_store_ps(Lap, t1); Lap += 4;
		}
#else
		for(n=0;n<(N-Overlap)/2;n++) {
			float a = *NewLo++;
			float b = *--NewHi;
			*--MDCTMid = a;
			*Lap++     = b;
		}
		for(;n<N/2;n++) {
			float a = *NewLo++;
			float b = *--NewHi;
			float c = *Win++;
			float s = *Win++;
			*--MDCTMid =  c*a + s*b;
			*Lap++     = -s*a + c*b;
		}
#endif
	}

	//! Do actual transforms
	Fourier_DCT4T(MDCT, BufTmp, N);
}

/**************************************/
//! EOF
/**************************************/
