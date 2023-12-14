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

#include "../Fourier_MDCT_Template.h"
void Fourier_MDCT_MDST_AVX(float *MDCT, float *MDST, const float *New, float *Lap, float *BufTmp, int N, int Overlap) {
	Fourier_MDCT_MDST_Template(MDCT, MDST, New, Lap, BufTmp, N, Overlap);
}

/************************************************/
#endif
/************************************************/
// EOF
/************************************************/
