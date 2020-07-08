/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <math.h>
/**************************************/
#define ABS(x) ((x) < 0 ? (-(x)) : (x))
#define SQR(x) ((x)*(x))
/**************************************/

//! Pattern for decimated subblocks
//! Usage:
//!  PatternBase = &PatternTable[DecimationNybble >> 1];
//!  PatternNext = PatternBase;
//!  for(Coef in Coefs) {
//!   *Buffer = Coef; Buffer += BlockSize >> (*PatternNext++);
//!   if(Buffer >= BufferEnd) Buffer -= BufferSize-1, PatternNext = PatternBase;
//!  }
//! The decimation nybble comes from Block_Transform_GetLogOverlapScale(),
//! corresponding to the upper/second nybble for decimated blocks.
static inline const int8_t *ULC_Helper_SubBlockDecimationPattern(int Decimation) {
	static const int8_t Pattern[8][4] = {
		{0},          //! 0001: N/1
		{1, 1},       //! 001x: N/2,N/2
		{2, 2, 1},    //! 010x: N/4,N/4,N/2
		{1, 2, 2},    //! 011x: N/2,N/4,N/4
		{3, 3, 2, 1}, //! 100x: N/8,N/8,N/4,N/2
		{2, 3, 3, 1}, //! 101x: N/4,N/8,N/8,N/2
		{1, 3, 3, 2}, //! 110x: N/2,N/8,N/8,N/4
		{1, 2, 3, 3}, //! 111x: N/2,N/4,N/8,N/8
	};
	return Pattern[Decimation >> 1];
}

//! Transient subblock index from decimation pattern
//! Due to the way the decimation pattern is built, the
//! subblock index is conveniently obtained via a
//! POPCNT (minus one for the unary 'stop' bit).
static inline int ULC_Helper_TransientSubBlockIndex(int Decimation) {
	return __builtin_popcount(Decimation) - 1;
}

/**************************************/

//! Spectral flatness measure
//! Adapted from "Note on measures for spectral flatness"
//! DOI: 10.1049/el.2009.1977
//! Optimized derivation (for a single loop):
//!  Given a sequence x[0..N-1]:
//!  Let `f` be the output value, and `b` an arbitrary base:
//!   f = 2^(-g / log_b(N)) - 1
//!   where `g` = Sum(y[n]*log_b(y[n]))
//!   and y = x[n] / Sum(x[0..N])
//!  Letting `a` = Sum(x[0..N]) and replacing `y` by its
//!  definition, we then get:
//!   g = Sum(x[n]/a * log_b(x[n]/a))
//!  Use logarithm rules to expand:
//!   g = Sum(x[n]/a * (log_b(x[n]) - log_b(a)))
//!   g = Sum(x[n] * log_b(x[n]) / a) - Sum(x[n] * log_b(a) / a)
//!  Factor out `a` and `log_b(a)`:
//!   g = (1/a)*Sum(x[n] * log_b(x[n])) - (1/a)*log_b(a)*Sum(x[n])
//!  And since `a` = Sum(x[0..N]):
//!   g = (1/a)*Sum(x[n] * log_b(x[n])) - log_b(Sum(x[n]))
//!  This means that we only need to compute two sums:
//!   a = Sum(x[n])
//!   f = Sum(x[n] * log_b(x[n]))
//!  and combine as
//!   g = (1/a)*f - log_b(a)
//!  As the sums are independent, this can be performed in a single step.
//! NOTE: Unused, and only here for reference
static inline __attribute__((always_inline)) float ULC_Helper_SpectralFlatness(const float *Buf, int N) {
	int i;
	float a = 0.0f, f = 0.0f;
	for(i=0;i<N;i++) {
		float Val = SQR(Buf[i]);
		a += Val;
		if(Val != 0.0f) f += Val*logf(Val);
	}
	if(a == 0.0f) return 1.0f;
	return exp2f((logf(a) - f/a) / logf(N)) - 1.0f;
}

//! Masking bandwidth estimation
//! NOTE: Unused, and only here for reference
static inline __attribute__((always_inline)) float ULC_Helper_MaskingBandwidth(float Fc) {
#if 0 //! Bark scale
	return 50.21f + Fc*(1.0f/8.73f + Fc*(1.0f/93945.23f));
#else
	//! ERB scale
	return 24.7f + 0.107939f*Fc;
#endif
}

/**************************************/
//! EOF
/**************************************/
