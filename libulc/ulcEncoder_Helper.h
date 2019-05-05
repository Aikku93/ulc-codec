/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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

//! Smoothstep curve
static inline float SmoothStep(float x) {
	return (3.0f - 2.0f*x)*x*x;
}

//! Spectral flatness measure
//! Adapted from "Note on measures for spectral flatness"
//! DOI: 10.1049/el.2009.1977
static inline double SpectralFlatness(const float *Buf, size_t N) {
	size_t i;

	//! Get normalization factor
	double Nrm = 0.0;
	for(i=0;i<N;i++) Nrm += SQR((double)Buf[i]);
	if(Nrm == 0.0) return 0.0;
	Nrm = 1.0 / Nrm;

	//! Measure
	double Sum = 0.0;
	for(i=0;i<N;i++) {
		double v = SQR((double)Buf[i]) * Nrm;
		if(v != 0.0) Sum += v*log(v);
	}
	return pow(2.0, -Sum/log(N)) - 1.0;
}

/**************************************/
//! EOF
/**************************************/
