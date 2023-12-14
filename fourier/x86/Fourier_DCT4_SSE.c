/************************************************/
#ifndef __SSE__
# define __SSE__
#endif
/************************************************/
#include "Fourier.h"
#include "../Fourier_Helper.h"
/************************************************/
#if (defined(FOURIER_IS_X86) && defined(FOURIER_ALLOW_SSE))
/************************************************/

#include "../Fourier_DCT4_Template.h"
void Fourier_DCT4_SSE(float *Buf, float *Tmp, int N) {
	Fourier_DCT4_Template(Buf, Tmp, N);
}

/************************************************/
#endif
/************************************************/
// EOF
/************************************************/
