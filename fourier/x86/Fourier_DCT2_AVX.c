/************************************************/
#ifndef __AVX__
# define __AVX__
#endif
/************************************************/
#include "Fourier.h"
#include "../Fourier_Helper.h"
/************************************************/
#if (defined(FOURIER_IS_X86) && defined(FOURIER_ALLOW_SSE))
/************************************************/

#include "../Fourier_DCT2_Template.h"
void Fourier_DCT2_AVX(float *Buf, float *Tmp, int N) {
	Fourier_DCT2_Template(Buf, Tmp, N);
}

/************************************************/
#endif
/************************************************/
// EOF
/************************************************/
