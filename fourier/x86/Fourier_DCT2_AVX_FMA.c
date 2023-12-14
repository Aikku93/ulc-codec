/************************************************/
#ifndef __AVX__
# define __AVX__
#endif
#ifndef __FMA__
# define __FMA__
#endif
/************************************************/
#include "Fourier.h"
#include "../Fourier_Helper.h"
/************************************************/
#if (defined(FOURIER_IS_X86) && defined(FOURIER_ALLOW_AVX) && defined(FOURIER_ALLOW_FMA))
/************************************************/

#include "../Fourier_DCT2_Template.h"
void Fourier_DCT2_AVX_FMA(float *Buf, float *Tmp, int N) {
	Fourier_DCT2_Template(Buf, Tmp, N);
}

/************************************************/
#endif
/************************************************/
// EOF
/************************************************/
