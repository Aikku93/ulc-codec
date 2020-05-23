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
//! NOTE: Flatness is returned as its reciprocal, as this seems to be
//! more numerically stable.
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
	State->SumLin = SQR((double)Coef[0]);
	State->SumLog = SQR((double)Coef[0]) * CoefNp[0];
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
	float *_InvFlat
) {
	//! These settings are mostly based on trial and error
	float Bw = State->BwBase + 0.107939f*Band;
	float MaskBW   = 0.75f;
	float MaskFW   = 0.25f;
	int   BandBeg  = Band - (int)ceilf(Bw * MaskBW);
	int   BandEnd  = Band + (int)ceilf(Bw * MaskFW);
	if(BandBeg <          0) BandBeg = 0;
	if(BandEnd >= BlockSize) BandEnd = BlockSize-1;

	//! Update masking calculations
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

	//! Get final masking threshold and flatness from the intermediates
	//! NOTE: Multiply by two in the flatness equation ONLY. This will
	//! return the correct masking value, while still using the squared
	//! coefficients necessary for the calculation, as that equation
	//! relies on kind of 'comparing' the log/linear values so they must
	//! match. The reason for this is that log(x^2) = 2*log(x)
	//! NOTE: SumLin can never be 0, as this update routine is only called
	//! upon encountering a valid non-zero coefficient. Therefore, this
	//! coefficient will always allow the routine to succeed
	double LogSumLin = log(SumLin);
	double LogMaskBw = log(BandEnd - BandBeg + 1);
	double nSumLog   = SumLog/SumLin;
	*_InvFlat = (float)(LogMaskBw / (LogSumLin - 2.0*nSumLog));
	return nSumLog;
}

/**************************************/
//! EOF
/**************************************/
