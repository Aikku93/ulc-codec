/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <math.h>
#include <stdint.h>
/**************************************/
#include "ulcEncoder.h"
#include "ulcEncoder_Helper.h"
/**************************************/

//! Masking threshold estimation
//! NOTE: ERB calculations are inlined as an optimization
struct Block_Transform_MaskingState_t {
	double SumLin;
	double SumLog;
	float BwBase;
	int   BandBeg;
	int   BandEnd;
};
static inline void Block_Transform_MaskingState_Init(
	struct Block_Transform_MaskingState_t *State,
	const float *Coef,
	const float *CoefNp,
	int   BlockSize,
	float NyquistHz
) {
	State->SumLin = SQR(Coef[0]);
	State->SumLog = SQR(Coef[0]) * CoefNp[0];
	State->BwBase = 24.7f*BlockSize/NyquistHz;
	State->BandBeg = 0;
	State->BandEnd = 0;
}
static inline float Block_Transform_GetMaskedLevel(
	struct Block_Transform_MaskingState_t *State,
	const float *Coef,
	const float *CoefNp,
	int   Band,
	int   BlockSize
) {
	//! These settings are mostly based on trial and error
	float Bw  = State->BwBase + 0.107939f*Band;
	int   BandBeg = Band - (int)(Bw * 0.25f);
	int   BandEnd = Band + (int)(Bw * 0.75f);
	if(BandBeg <          0) BandBeg = 0;
	if(BandEnd >= BlockSize) BandEnd = BlockSize-1;

	//! Update masking calculations
	//! NOTE: Double-precision is necessary here to reduce
	//! catastrophic cancellation. Integer maths would be
	//! ideal (as all subtractions would cancel exactly),
	//! but may give suboptimal results. Alternatively, a
	//! convolution may be performed per call to this
	//! function, but this is much more costly.
	double SumLin = State->SumLin;
	double SumLog = State->SumLog;
	{
#define ADD_TO_STATE(Val, ValNp) SumLin += SQR(Val), SumLog += SQR(Val)*(ValNp)
#define SUB_TO_STATE(Val, ValNp) SumLin -= SQR(Val), SumLog -= SQR(Val)*(ValNp)
		int Beg = State->BandBeg, End = State->BandEnd;
		while(BandBeg < Beg) --Beg, ADD_TO_STATE(Coef[Beg], CoefNp[Beg]); //! Expand
		while(BandEnd > End) ++End, ADD_TO_STATE(Coef[End], CoefNp[End]);
		while(BandBeg > Beg) SUB_TO_STATE(Coef[Beg], CoefNp[Beg]), Beg++; //! Contract
		while(BandEnd < End) SUB_TO_STATE(Coef[End], CoefNp[End]), End--;
		State->BandBeg = BandBeg, State->BandEnd = BandEnd;
#undef SUB_TO_STATE
#undef ADD_TO_STATE
	}
	State->SumLin = SumLin;
	State->SumLog = SumLog;

	//! Get final masked amplitude after mapping from the psychoacoustic
	//! domain to the linear domain
	//! NOTE: SumLin can never be 0, as this update routine is only
	//! called upon encountering a valid non-zero coefficient, which is
	//! included in the calculations. Thus, division by 0 cannot occur.
	//! NOTE: The idea behind these equations is to translate and scale
	//! the coefficient power using the classical formula of:
	//!  x' = NewOrigin + Scale*(x - OldOrigin)
	//! OldOrigin is the psychoacoustic domain, and NewOrigin is the
	//! linear domain. Scale is set in such a way that the power of
	//! the coefficients scales proportionally with the psychoacoustic
	//! power such that it matches overall.
	//! NOTE: Using squared inputs seems to be critical for maximizing
	//! brightness at low bitrates.
	float nSumLog = 2.0f * SumLog / SumLin;
	float LogSumLin = logf(SumLin / (BandEnd - BandBeg + 1));
	float Curve = expf(nSumLog - LogSumLin);
	return LogSumLin + Curve*(2.0f*CoefNp[Band] - nSumLog);
}

/**************************************/
//! EOF
/**************************************/
