/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
/**************************************/

//! Integer floor(log2(x))
//! Using unsigned int because that's what the compiler expects
static inline size_t IntLog2(unsigned int x) {
	return sizeof(unsigned int)*8-1 - __builtin_clz(x);
}

/**************************************/
//! EOF
/**************************************/
