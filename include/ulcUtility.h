/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
/**************************************/

//! Apply pre-echo reduction formula
//! Based on the Haar wavelet, just use a sum/difference
//! This will make low-frequency transients compact into
//! the even elements, and high-frequency transients will
//! go into the odd frequency components, resulting in
//! better transient coding at lower bitrates
void ULC_Transform_AntiPreEcho(float *Buf, size_t N);

/**************************************/
//! EOF
/**************************************/
