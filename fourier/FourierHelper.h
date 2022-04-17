/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2022, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#define FOURIER_FORCED_INLINE static inline __attribute__((always_inline))
#define FOURIER_ASSUME(Cond) (Cond) ? ((void)0) : __builtin_unreachable()
#define FOURIER_ASSUME_ALIGNED(x,Align) x = __builtin_assume_aligned(x,Align)
/**************************************/

//! Sine table for DCT analysis
//! Contains E^(I*(n+0.5)*Pi/N) for n=0..N/2-1 and
//! N={16,32,64,128,256,512,1024,2048,4096,8192}.
//! Arranged as Re,Im per element.
FOURIER_FORCED_INLINE const float *Fourier_SinTableN(int N, const float *Table) {
	extern const float Fourier_SinTable[];

	//! NOTE: N must be > 8
	if(!Table) Table = Fourier_SinTable;
	return Table + (N-16);
}

/**************************************/
//! EOF
/**************************************/
