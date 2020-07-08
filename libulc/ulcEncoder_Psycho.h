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
#include "ulcHelper.h"
/**************************************/

//! Simultaneous masking-compensated energy estimation
//! NOTE: ERB calculations are inlined as an optimization
struct Block_Transform_MaskingState_t {
	double Energy;
	double Nepers;
	float  BwBase;
	int    BandBeg;
	int    BandEnd;
};
static inline void Block_Transform_MaskingState_Init(
	struct Block_Transform_MaskingState_t *State,
	const float *Energy,
	const float *EnergyNp,
	int   BlockSize,
	float NyquistHz
) {
	State->Energy = Energy[0];
	State->Nepers = Energy[0] * EnergyNp[0];
	State->BwBase  = 24.7f*BlockSize/NyquistHz + 0.0539695f; //! Offset at Band+0.5
	State->BandBeg = 0;
	State->BandEnd = 0;
}
static inline float Block_Transform_GetMaskedLevel(
	struct Block_Transform_MaskingState_t *State,
	const float *Energy,
	const float *EnergyNp,
	int   Band,
	int   BlockSize,
	float Gamma
) {
	//! These settings are mostly based on trial and error
	int BandBeg, BandEnd; {
		//! Because lower frequencies mask higher frequencies more
		//! than higher frequencies mask lower ones, the frequency
		//! bands below the current band are critical for analysis
		//! and the masking bandwidth is thus much larger there
		float Bw = State->BwBase + 0.107939f*Band;
		BandBeg = Band - (int)(Bw * 0.5f + 0x1.FFFFFFp-1f); //! Add 0.99999... for ceiling
		BandEnd = Band + (int)(Bw * 0.5f + 0x1.FFFFFFp-1f);
		if(BandBeg <          0) BandBeg = 0;
		if(BandEnd >= BlockSize) BandEnd = BlockSize-1;
	}

	//! Update energy for this critical band
	//! NOTE: Double-precision is necessary here to reduce
	//! catastrophic cancellation. Integer maths would be
	//! ideal (as all subtractions would cancel exactly),
	//! but the weighted geometric mean we are computing
	//! here is hard to get right. Alternatively, a
	//! convolution may be performed per call to this
	//! function, but this is much more costly.
	double EnergySum = State->Energy;
	double EnergyLog = State->Nepers;
	{
#define ADD_TO_STATE(Band) EnergySum += Energy[Band], EnergyLog += Energy[Band] * EnergyNp[Band]
#define SUB_TO_STATE(Band) EnergySum -= Energy[Band], EnergyLog -= Energy[Band] * EnergyNp[Band]
		int Beg = State->BandBeg, End = State->BandEnd;
		while(BandBeg < Beg) --Beg, ADD_TO_STATE(Beg); //! Expand
		while(BandEnd > End) ++End, ADD_TO_STATE(End);
		while(BandBeg > Beg) SUB_TO_STATE(Beg), Beg++; //! Contract
		while(BandEnd < End) SUB_TO_STATE(End), End--;
		State->BandBeg = BandBeg, State->BandEnd = BandEnd;
#undef SUB_TO_STATE
#undef ADD_TO_STATE
	}
	State->Energy = EnergySum;
	State->Nepers = EnergyLog;

	//! De/emphasize band energy relative to the background
	//! level. This is similar in effect to a gamma correction
	//! filter in image processing, where we wish to use it to
	//! extract 'important' information from a noisy background.
	float BandNp = EnergyNp[Band];
	return BandNp + Gamma*(BandNp - EnergyLog/EnergySum);
}

/**************************************/
//! EOF
/**************************************/
