/************************************************/
#ifndef __AVX__
# define __AVX__
#endif
/************************************************/
#include "Fourier.h"
#include "../Fourier_Helper.h"
/************************************************/
#if (defined(FOURIER_IS_X86) && defined(FOURIER_ALLOW_AVX))
/************************************************/

#include "../Fourier_IMDCT_Template.h"
void Fourier_IMDCT_AVX(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, int N, int Overlap) {
	Fourier_IMDCT_Template(BufOut, BufIn, BufLap, BufTmp, N, Overlap);
}

/************************************************/
#endif
/************************************************/
// EOF
/************************************************/
