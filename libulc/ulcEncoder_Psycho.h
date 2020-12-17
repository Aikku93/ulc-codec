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
	float  BwBase;
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
	State->BwBase   = 24.7f*BlockSize/NyquistHz + 0.0539695f; //! Offset at Band+0.5
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
		//! NOTE: Bw is doubled relative to actual ERB as this
		//! seems to give slightly improved results.
		float Bw = State->BwBase + 0.107939f*Band;
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

	//! Return the entropy of this frequency bin's critical band,
	//! normalized about its L1 norm.
	//! NOTE: EnergySum must scale down, as EnergyLog cannot
	//! scale up the necessary bits without potential overflow.
	//! NOTE: This value will be subtracted from the Neper-domain
	//! amplitude of the MDCT coefficient. I'm not entirely
	//! sure what the idea is, but it seems that the implicit
	//! "exponents" (ie. scale, in the log domain) must add up
	//! to 1.0 in the end. The values we get here have an
	//! exponent of 2.0 (due to working with X^2 values), but
	//! the Neper MDCT coefficients have an exponent of 1.0,
	//! giving us the expression:
	//!  x - 2 == 1
	//! which obviously gives x=3. However, the actual scale
	//! doesn't matter for our purposes, as we are only using
	//! these values for their rank/importance, so we can
	//! simply divide everything by 3 - giving the output here
	//! an exponent of 1/3 - and fold this into the scaling
	//! for getting out of fixed-point representation, which
	//! should (hopefully) be more precise than scaling the
	//! output /after/ casting to float.
	//! NOTE: F3FCE0F5h == Floor[(1/LogScale)*(1/3) * 2^61 + 0.5]
	//! LogScale ((2^32) / Log[2*2^32]) is defined in Block_Transform().
	EnergySum = (EnergySum >> State->SumShift) + ((EnergySum << (64-State->SumShift)) != 0);
	uint64_t r = (EnergyLog/EnergySum) * 0xF3FCE0F5ull;
	return r * 0x1.0p-61f;
}

/**************************************/
//! EOF
/**************************************/
