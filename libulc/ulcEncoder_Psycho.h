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
//! NOTE: ERB and flatness calculations are inlined as an optimization
struct Block_Transform_MaskingState_t {
	float SumLin;
	float SumLog;
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
	State->BwBase = 24.7f*BlockSize/NyquistHz + 0.0539695f;
	State->BandBeg = 0;
	State->BandEnd = 0;
}
static inline float Block_Transform_UpdateMaskingThreshold(
	struct Block_Transform_MaskingState_t *State,
	const float *Coef,
	const float *CoefNp,
	int   Band,
	int   BlockSize,
	float *_Flat
) {
	//! These settings are mostly based on trial and error
	float Bw = State->BwBase + 0.107939f*Band;
	float MaskFW   = 1.0f - ABS(Coef[Band]);
	float MaskBW   = sqrtf(1.0f - MaskFW);
	int   BandEnd  = Band + (int)(Bw * MaskFW + 0.5f);
	int   BandBeg  = Band - (int)(Bw * MaskBW + 0.5f);
	if(BandEnd >= BlockSize) BandEnd = BlockSize-1;
	if(BandBeg <          0) BandBeg = 0;

	//! Update masking calculations
	{
#define ADD_TO_STATE(Val, ValNp) State->SumLin += SQR(Val), State->SumLog += SQR(Val)*(ValNp)
#define SUB_TO_STATE(Val, ValNp) State->SumLin -= SQR(Val), State->SumLog -= SQR(Val)*(ValNp)
		int Beg = State->BandBeg, End = State->BandEnd;
		while(BandBeg < Beg) Beg--, ADD_TO_STATE(Coef[Beg], CoefNp[Beg]); //! Expand
		while(BandEnd > End) End++, ADD_TO_STATE(Coef[End], CoefNp[End]);
		while(BandBeg > Beg) SUB_TO_STATE(Coef[Beg], CoefNp[Beg]), Beg++; //! Contract
		while(BandEnd < End) SUB_TO_STATE(Coef[End], CoefNp[End]), End--;
		State->BandBeg = BandBeg, State->BandEnd = BandEnd;
#undef SUB_TO_STATE
#undef ADD_TO_STATE
	}

	//! Get final masking threshold and flatness from the intermediates
	//! NOTE: Limit at (0.5*eps)^2; this would correspond to one single
	//! coefficient in the sum, so anything smaller than that would be
	//! rounding error from limited precision representation.
	//! NOTE: Multiply by two in the flatness equation ONLY. This will
	//! return the correct masking value, while still using the squared
	//! coefficients necessary for the calculation, as that equation
	//! relies on kind of 'comparing' the log/linear values so they must
	//! match. The reason for this is that log(x^2) = 2*log(x)
	//! NOTE: 0x1.62E430p-1f = Log[2]
	float  SumLin = State->SumLin; if(SumLin < SQR(0.5f*ULC_COEF_EPS)) { *_Flat = 1.0f; return 0.0f; }
	float nSumLog = State->SumLog / SumLin;
	*_Flat        = expf(0x1.62E430p-1f * ((logf(SumLin) - 2.0f*nSumLog) / logf(BandEnd - BandBeg + 1))) - 1.0f;
	return nSumLog;
}

/**************************************/
//! EOF
/**************************************/
