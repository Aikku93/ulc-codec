/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#if defined(__AVX__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/

//! Sine table for DCT analysis
//! Contains Table[Sin[(n+0.5)*(Pi/2)/N], {n,0,N-1}]
//! for N={16,32,64,128,256,512,1024,2048,4096,8192}
extern const float Fourier_SinTable[];
static inline __attribute__((always_inline)) const float *Fourier_SinTableN(int N) {
	//! NOTE: N must be > 8
	return Fourier_SinTable + (N-16);
}

/**************************************/

//! DCT-II/DCT-IV (scaled)
//! Arguments:
//!  Buf[N]
//!  Tmp[N]
//! Implemented transforms (matrix form):
//!  mtxDCTII = Table[Cos[(n-1/2)(k-1  )Pi/N], {k,N}, {n,N}]
//!  mtxDCTIV = Table[Cos[(n-1/2)(k-1/2)Pi/N], {k,N}, {n,N}]
//! Implementations from:
//!  "Signal Processing based on Stable radix-2 DCT I-IV Algorithms having Orthogonal Factors"
//!  DOI: 10.13001/1081-3810.3207
//! NOTE:
//!  -N must be a power of two, and >= 8
void Fourier_DCT2(float *Buf, float *Tmp, int N);
void Fourier_DCT4(float *Buf, float *Tmp, int N);

//! MDCT/IMDCT (based on DCT-IV; scaled)
//! Arguments:
//!  BufOut[N]
//!  BufIn[N]
//!  BufLap[N/2]
//!  BufTmp[N]
//! Implemented transforms (matrix form):
//!  mtxMDCT  = Table[Cos[(n-1/2 + N/2 + N*2)(k-1/2)Pi/N], {k, N}, {n,2N}]
//!  mtxIMDCT = Table[Cos[(n-1/2 + N/2 + N*2)(k-1/2)Pi/N], {k,2N}, {n, N}]
//! NOTE:
//!  -N must be a power of two, and >= 8
//!  -Overlap must be a multiple of 16
//!  -BufOut must not be the same as BufIn
//!  -Sine window (modulated lapped transform)
//!  -Shifted basis that results in phase inversion
//!   relative to 'normal' I/MDCT calculations.
//!   Negate input (MDCT) and output (IMDCT) if
//!   'correct' I/MDCT coefficients are needed
void Fourier_MDCT (float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, int N, int Overlap);
void Fourier_IMDCT(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, int N, int Overlap);

/**************************************/
//! EOF
/**************************************/
