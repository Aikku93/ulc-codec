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
#include <stddef.h>
/**************************************/

//! Sine/Cosine approximations
//! Approximates {Sin,Cos}[x*Pi/2] (with x in range [-1,1])
#if defined(__AVX__)
void Fourier_SinCosAVX(__m256 x, __m256 *Sin, __m256 *Cos);
#endif
#if defined(__SSE__)
void Fourier_SinCosSSE(__m128 x, __m128 *Sin, __m128 *Cos);
#endif
float Fourier_SinCos(float x, float *Sin, float *Cos);

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
void Fourier_DCT2(float *Buf, float *Tmp, size_t N);
void Fourier_DCT4(float *Buf, float *Tmp, size_t N);

//! MDCT/IMDCT (based on DCT-IV; scaled)
//! Arguments:
//!  BufOut[N]
//!  BufIn[N]
//!  BufLap[N/2]
//!  BufTmp[N]
//! Implemented transforms (matrix form):
//!  mtxMDCT  = Table[Cos[(n-1/2 + N/2)(k-1/2)Pi/N], {k, N}, {n,2N}]
//!  mtxIMDCT = Table[Cos[(n-1/2 + N/2)(k-1/2)Pi/N], {k,2N}, {n, N}]
//! NOTE:
//!  -N must be a power of two, and >= 8
//!  -Overlap must be a multiple of 16
//!  -BufOut must not be the same as BufIn
//!  -Sine window (modulated lapped transform)
void Fourier_MDCT (float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, size_t N, size_t Overlap);
void Fourier_IMDCT(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, size_t N, size_t Overlap);

/**************************************/
//! EOF
/**************************************/
