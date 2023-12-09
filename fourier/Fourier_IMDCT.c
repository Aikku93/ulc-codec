/************************************************/
#include "Fourier.h"
#include "Fourier_Helper.h"
/************************************************/

#define FUNC_ARGS_LIST float *BufOut, const float *BufIn, float *BufLap, float *BufTmp, int N, int Overlap

/************************************************/

// Declare specialized routines
void Fourier_IMDCT_Generic(FUNC_ARGS_LIST);
#ifdef FOURIER_IS_X86
# ifdef FOURIER_ALLOW_SSE
void Fourier_IMDCT_SSE    (FUNC_ARGS_LIST);
#  ifdef FOURIER_ALLOW_FMA
void Fourier_IMDCT_SSE_FMA(FUNC_ARGS_LIST);
#  endif
# endif
# ifdef FOURIER_ALLOW_AVX
void Fourier_IMDCT_AVX    (FUNC_ARGS_LIST);
#  ifdef FOURIER_ALLOW_FMA
void Fourier_IMDCT_AVX_FMA(FUNC_ARGS_LIST);
#  endif
# endif
#endif

/************************************************/

static const struct Fourier_FuncTbl_t DispatchTbl = {
	.Generic = Fourier_IMDCT_Generic,
#ifdef FOURIER_IS_X86
# ifdef FOURIER_ALLOW_SSE
	.SSE = Fourier_IMDCT_SSE,
#  ifdef FOURIER_ALLOW_FMA
	.SSE_FMA = Fourier_IMDCT_SSE_FMA,
#  endif
# endif
# ifdef FOURIER_ALLOW_AVX
	.AVX = Fourier_IMDCT_AVX,
#  ifdef FOURIER_ALLOW_FMA
	.AVX_FMA = Fourier_IMDCT_AVX_FMA,
#  endif
# endif
#endif
};

/************************************************/

typedef void (*Dispatcher_t)(FUNC_ARGS_LIST);
static void InitDispatcher(FUNC_ARGS_LIST);

static Dispatcher_t Dispatcher = InitDispatcher;

static void InitDispatcher(FUNC_ARGS_LIST) {
	Dispatcher = (Dispatcher_t)Fourier_GetDispatchFnc(&DispatchTbl);
	Dispatcher(BufOut, BufIn, BufLap, BufTmp, N, Overlap);
}

/************************************************/

void Fourier_IMDCT(FUNC_ARGS_LIST) {
	Dispatcher(BufOut, BufIn, BufLap, BufTmp, N, Overlap);
}

/************************************************/

#include "Fourier_IMDCT_Template.h"
void Fourier_IMDCT_Generic(FUNC_ARGS_LIST) {
	Fourier_IMDCT_Template(BufOut, BufIn, BufLap, BufTmp, N, Overlap);
}

/************************************************/
// EOF
/************************************************/
