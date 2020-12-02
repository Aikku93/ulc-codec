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
	float  BwBase; //! Scaled by 0.5
	int    SumShift;
	int    BandBeg;
	int    BandEnd;
	uint64_t Energy;
	uint64_t Nepers;
};
static inline void Block_Transform_MaskingState_Init(
	struct Block_Transform_MaskingState_t *State,
	const uint32_t *Energy,
	const uint32_t *EnergyNp,
	int   BlockSize,
	float NyquistHz
) {
	State->BwBase   = 12.35f*BlockSize/NyquistHz + 0.02698475f; //! Offset at Band+0.5
	State->SumShift = 31 - __builtin_clz(BlockSize); //! Log2[BlockSize]
	State->BandBeg  = 0;
	State->BandEnd  = 0;
	State->Energy   = Energy[0];
	State->Nepers   = Energy[0] * (uint64_t)EnergyNp[0] >> State->SumShift;
}
static inline float Block_Transform_GetMaskedLevel(
	struct Block_Transform_MaskingState_t *State,
	const uint32_t *Energy,
	const uint32_t *EnergyNp,
	int   Band,
	int   BlockSize
) {
	int BandBeg, BandEnd; {
		float Bw = State->BwBase + 0.0539695f*Band;
		BandBeg = Band - (int)(Bw + 0x1.FFFFFFp-1f); //! Add 0.99999... for ceiling
		BandEnd = Band + (int)(Bw + 0x1.FFFFFFp-1f);
		if(BandBeg <          0) BandBeg = 0;
		if(BandEnd >= BlockSize) BandEnd = BlockSize-1;
	}

	//! Update energy for this critical band
	//! NOTE: The maximum value for Energy*EnergyNp correspond to the values
	//! Energy == FFFFFFFFh, EnergyNp == FFFFFFFFh. Being pessimistic, we
	//! can assume that the bandwidth of the sum is the maximum allowable
	//! block size (8192), thus we must scale down the EnergyLog sum by
	//! BlockSize, which conveniently can be performed with a LSR.
	uint64_t EnergySum = State->Energy;
	uint64_t EnergyLog = State->Nepers;
	{
#define ADD_TO_STATE(Band) EnergySum += Energy[Band], EnergyLog += Energy[Band] * (uint64_t)EnergyNp[Band] >> State->SumShift
#define SUB_TO_STATE(Band) EnergySum -= Energy[Band], EnergyLog -= Energy[Band] * (uint64_t)EnergyNp[Band] >> State->SumShift
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
	//! This gamma value was derived experimentally.
	//! NOTE: Do not return the 'gamma-corrected' value, but
	//! instead return the correction factor.
	//! NOTE: EnergySum must scale down, as EnergyLog cannot
	//! scale up the necessary bits without potential overflow.
	const float InvLogScale = 0x1.0p-27f; //! Reciprocal of LogScale in Block_Transform_ComputePowerSpectrum()
	EnergySum = (EnergySum >> State->SumShift) + ((EnergySum << (64-State->SumShift)) != 0); //! Round up
	return logf(BandEnd-BandBeg+1) * InvLogScale*(int64_t)(EnergyNp[Band] - EnergyLog/EnergySum);
}

/**************************************/
//! EOF
/**************************************/
