/************************************************/
#ifndef __SSE__
# define __SSE__
#endif
#ifndef __FMA__
# define __FMA__
#endif
/************************************************/
#include "Fourier.h"
#include "../Fourier_Helper.h"
/************************************************/
#if (defined(FOURIER_IS_X86) && defined(FOURIER_ALLOW_SSE) && defined(FOURIER_ALLOW_FMA))
/************************************************/

#include "../Fourier_IMDCT_Template.h"
void Fourier_IMDCT_SSE_FMA(float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, int N, int Overlap) {
	Fourier_IMDCT_Template(BufOut, BufIn, BufLap, BufTmp, N, Overlap);
}

/************************************************/
#endif
/************************************************/
// EOF
/************************************************/
