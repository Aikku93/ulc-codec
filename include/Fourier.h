/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
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
//!   is based on direct reversal of the steps in DCT4(). This
//!   is mathematically equivalent to DCT4(), but may result in
//!   better/worse round-off error, depending on the input.
//!   Recommend using DCT4T() for converting time-domain signals
//!   to the frequency domain, and DCT4() for the inverse.
void Fourier_DCT2 (float *Buf, float *Tmp, int N);
void Fourier_DCT3 (float *Buf, float *Tmp, int N);
void Fourier_DCT4 (float *Buf, float *Tmp, int N);
void Fourier_DCT4T(float *Buf, float *Tmp, int N);

//! MDCT+MDST/IMDCT (based on DCT-IV; scaled)
//! Arguments:
//!  MDCT[N]
//!  MDST[N] (for Fourier_MDCT_MDST() only)
//!  New[N]
//!  Lap[N]  (for Fourier_MDCT() only: Lap[N/2])
//!  BufTmp[N]
//! Implemented transforms (matrix form):
//!  mtxMDCT = Table[Cos[(n-1/2 - N/2 + N*2)(k-1/2)Pi/N], {k, N}, {n,2N}]
//!  mtxMDST = Table[Sin[(n-1/2 + N/2 + N*2)(k-1/2)Pi/N], {k, N}, {n,2N}]
//! NOTE:
//!  -N must be a power of two, and >= 16
//!  -Overlap must be a power of two, and >= 16
//!  -Shifted basis (note the signs in the matrices)
//!  -Sine window (modulated lapped transform) is
//!   used with ModulationWindow == NULL. To use
//!   custom windows, they must match the format of
//!   Fourier_SinTable[] (that is, sequential arrays
//!   for different block sizes 16..N, in powers of 2)
//!  -New can be the same as BufTmp. However, this
//!   implies trashing of the buffer contents.
//!  -MDCT uses Fourier_DCT4T() internally
void Fourier_MDCT_MDST(float *MDCT, float *MDST, const float *New, float *Lap, float *BufTmp, int N, int Overlap, const float *ModulationWindow);
void Fourier_MDCT     (float *MDCT,              const float *New, float *Lap, float *BufTmp, int N, int Overlap, const float *ModulationWindow);

//! IMDCT (based on DCT-IV; scaled)
//! Arguments:
//!  BufOut[N]
//!  BufIn[N]
//!  BufLap[N/2]
//!  BufTmp[N]
//! Implemented transforms (matrix form):
//!  mtxIMDCT = Table[Cos[(n-1/2 - N/2 + N*2)(k-1/2)Pi/N], {k,2N}, {n, N}]
//! NOTE:
//!  -N must be a power of two, and >= 16
//!  -Overlap must be a power of two, and >= 16
//!  -BufOut must not be the same as BufIn
//!  -Shifted basis (note the sign in the matrix).
//!   This re-inverts (ie. removes) the phase shift
//!   in the MDCT implementation from above.
//!  -ModulationWindow operates same as Fourier_MDCT()
//!  -BufIn can be the same as BufTmp. However, this
//!   implies trashing of the buffer contents.
//!  -IMDCT uses Fourier_DCT4() internally
void Fourier_IMDCT(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, int N, int Overlap, const float *ModulationWindow);

/**************************************/
//! EOF
/**************************************/
