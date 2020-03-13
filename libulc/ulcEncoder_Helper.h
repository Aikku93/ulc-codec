/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
/**************************************/
#define ABS(x) ((x) < 0 ? (-(x)) : (x))
#define SQR(x) ((x)*(x))
/**************************************/

//! Integer floor(log2(x))
//! Using unsigned int because that's what the compiler expects
static inline size_t IntLog2(unsigned int x) {
	return sizeof(unsigned int)*8-1 - __builtin_clz(x);
}

//! Spectral flatness measure
//! Adapted from "Note on measures for spectral flatness"
//! and then optimized, and slightly tweaked (specifically,
//! normalization was factored out to perform a single pass;
//! this is not technically the same, but the results are
//! extremely close, and so it's worth optimizing like this)
//! DOI: 10.1049/el.2009.1977
static inline float SpectralFlatness(const float *Buf, size_t N) {
	size_t i;
	float Sum = 0.0f;
	float Nrm = 0.0f;
	for(i=0;i<N;i++) {
		float v = SQR(Buf[i]);
		Nrm += v; if(v != 0.0f) Sum += v*logf(v);
	}
	return exp2f((Nrm*logf(Nrm) - Sum) / (Nrm * logf(N))) - 1.0f;
}

//! Masking bandwidth estimation
//! The curve peaks at around 3000Hz, corresponding to human hearing
static inline __attribute__((always_inline)) float MaskingBandwidth(float Fc) {
	return 1000.0f * expf(-2.0f*(float)M_PI * SQR(Fc/16000.0f)) * (1.0f - expf(-2.0f*(float)M_PI * SQR((Fc + 80.0f)/5000.0f)));
}

//! NOTE: Unused, and mostly for reference
//! ERB scale-based (Moore and Glasberg, 1990) masking bandwidth estimation, but made into an exponential curve
static inline __attribute__((always_inline)) float MaskingBandwidth_ERB(float Fc) {
	return 1200.0f * (1.0f - expf((float)(-2.0*M_PI/SQR(16000.0f)) * SQR(Fc)));
}

/**************************************/
//! EOF
/**************************************/
