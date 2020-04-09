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

//! Maximum possible coding rate
//! NOTE: Does NOT account for quantizers, as these mess up the calculations
static inline __attribute__((always_inline)) float MaxCodingKbps(int BlockSize, int nChan, int RateHz) {
	return nChan*(4*BlockSize) * (float)RateHz/BlockSize * (1.0f/1000.0f);
}

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
static inline __attribute__((always_inline)) float SpectralFlatness(const float *Buf, int N) {
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
static inline __attribute__((always_inline)) float MaskingBandwidth(float Fc) {
#if 0 //! Bark scale
	return 50.21f + Fc*(1.0f/8.73f + Fc*(1.0f/93945.23f));
#else
	//! ERB scale
	return 24.7f + 0.107939f*Fc;
#endif
}

//! A-Weight approximation, in Nepers
static inline __attribute__((always_inline)) float AWeightNp(float Fc) {
	if(Fc < 1000.0f) return 0x1.95BC61p-2f - 420.0f/(60.0f + Fc);
	if(Fc < 6000.0f) {
		float f = logf((Fc - 1000.0f) * (1.0f/5000.0f));
		return -0.35f * expf(0.87f*f) * f;
	}
	return (Fc - 6000.0f) * (-1.0f/13000.0f);
}

/**************************************/
//! EOF
/**************************************/
