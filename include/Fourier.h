/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/

//! Sine table for DCT analysis
//! Contains Table[Sin[(n+0.5)*(Pi/2)/N], {n,0,N-1}]
//! for N={16,32,64,128,256,512,1024,2048,4096,8192}
static inline __attribute__((always_inline)) const float *Fourier_SinTableN(int N) {
	extern const float Fourier_SinTable[];

	//! NOTE: N must be > 8
	return Fourier_SinTable + (N-16);
}

/**************************************/

//! DCT-II/DCT-III/DCT-IV (scaled)
//! Arguments:
//!  Buf[N]
//!  Tmp[N]
//! Implemented transforms (matrix form):
//!  mtxDCTII  = Table[Cos[(n-1/2)(k-1  )Pi/N], {k,N}, {n,N}]
//!  mtxDCTIII = Table[Cos[(n-1  )(k-1/2)Pi/N], {k,N}, {n,N}]
//!  mtxDCTIV  = Table[Cos[(n-1/2)(k-1/2)Pi/N], {k,N}, {n,N}]
//! Implementations from:
//!  "Signal Processing based on Stable radix-2 DCT I-IV Algorithms having Orthogonal Factors"
//!  DOI: 10.13001/1081-3810.3207
//! NOTE:
//!  -N must be a power of two, and >= 8
//!  -DCT3() is the transposed version of DCT2(); all its code
//!   is based on direct reversal of the steps in DCT2().
//!  -DCT4T() is the transposed version of DCT4(); all its code
//!   is based on direct reversal of the steps in DCT4().
void Fourier_DCT2 (float *Buf, float *Tmp, int N);
void Fourier_DCT3 (float *Buf, float *Tmp, int N);
void Fourier_DCT4 (float *Buf, float *Tmp, int N);
void Fourier_DCT4T(float *Buf, float *Tmp, int N);

//! MDCT+MDST/IMDCT (based on DCT-IV; scaled)
//! Arguments:
//!  MDCT[N]
//!  MDST[n]
//!  New[N]
//!  Lap[N]
//!  BufOut[N]
//!  BufIn[N]
//!  BufLap[N/2]
//!  BufTmp[N]
//! Implemented transforms (matrix form):
//!  mtxMDCT  = Table[Cos[(n-1/2 + N/2 + N*2)(k-1/2)Pi/N], {k, N}, {n,2N}]
//!  mtxMDST  = Table[Sin[(n-1/2 + N/2 + N*2)(k-1/2)Pi/N], {k, N}, {n,2N}]
//!  mtxIMDCT = Table[Cos[(n-1/2 + N/2 + N*2)(k-1/2)Pi/N], {k,2N}, {n, N}]
//! NOTE:
//!  -N must be a power of two, and >= 16
//!  -Overlap must be a multiple of 16
//!  -BufOut must not be the same as BufIn
//!  -Sine window (modulated lapped transform)
//!  -Shifted basis that results in phase inversion
//!   relative to 'normal' I/MDCT calculations.
//!   Negate input (MDCT) and output (IMDCT) if
//!   'correct' I/MDCT coefficients are needed.
//!   MDST coefficients are NOT phase inverted.
//!  -BufIn can be the same as BufTmp. However, this
//!   implies trashing of the buffer contents.
//!  -MDCT uses DCT4T(), IMDCT uses DCT4()
void Fourier_MDCT_MDST(float *MDCT, float *MDST, const float *New, float *Lap, float *BufTmp, int N, int Overlap);
void Fourier_IMDCT(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, int N, int Overlap);

/**************************************/
//! EOF
/**************************************/
