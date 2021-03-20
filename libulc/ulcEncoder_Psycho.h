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
struct Block_Transform_MaskingState_t {
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
	int   BlockSize
) {
	//! The maximum bandwidth of a masking band is 0.25*BlockSize,
	//! meaning that we can include an extra 2 bits of precision.
	//! 0.25*BlockSize is never achieved in practice because the
	//! block is only so big, but this gives us an upper limit.
	State->SumShift = 31-2 - __builtin_clz(BlockSize); //! Log2[BlockSize] - UnusedBits
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
		//! These curves are similar to the Bark-scale bandwidths,
		//! with the assumption that the Bark bands are not discrete
		//! NOTE: Offset at Band+0.5.
		//! NOTE: Round up the end band so that we always have at
		//! least two bands to analyze for spectral flatness.
		BandBeg = (int)(0.90f*Band);
		BandEnd = (int)(1.15f*Band);
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
	//! NOTE: This value will be subtracted from the Neper-domain
	//! amplitude of the MDCT coefficient. The overall idea of
	//! the implemented model is to form this equation:
	//!  ImportanceLevel = CoefficientAmplitude * BandPower/BackgroundPower
	//! In the log domain, this becomes:
	//!  Log[ImportanceLevel] = Log[CoefficientAmplitude] + Log[BandPower]-Log[BackgroundPower]
	//! In this function, we are calculating Log[BackgroundPower],
	//! and thus we need Log[CoefficientAmplitude] as well as
	//! Log[BandPower]. From testing, using the MDCT+MDST
	//! squared amplitude as BandPower results in inferior
	//! results compared to the squared amplitude of
	//! CoefficientAmplitude. So we can then state:
	//!  Log[ImportanceLevel] = Log[CoefficientAmplitude] + Log[CoefficientAmplitude^2]-Log[BackgroundPower]
	//!                       = Log[CoefficientAmplitude] + 2*Log[CoefficientAmplitude]-Log[BackgroundPower]
	//!                       = 3*Log[CoefficientAmplitude] - Log[BackgroundPower]
	//! And finally, since ImportanceLevel is scale-invariant
	//! for our purposes, we divide by 3 on both sides:
	//!  Log[ImportanceLevel]/3 = Log[CoefficientAmplitude] - Log[BackgroundPower]/3
	//! NOTE: EnergySum must scale down, as EnergyLog cannot
	//! scale up the necessary bits without potential overflow.
	//! NOTE: F3FCE0F5h == Floor[0.5 + 2^61*(1/LogScale)/3]
	//! LogScale ((2^32) / Log[2*2^32]) is defined in Block_Transform().
	//! NOTE: Energy[] was calculated without raising to any power,
	//! whereas these calculations depend on it being raised to the
	//! 2nd power, so we compensate by scaling by 2.0 here.
	EnergySum = (EnergySum >> State->SumShift) + ((EnergySum << (64-State->SumShift)) != 0);
	if(EnergySum == 0) return 0.0f; //! <- Doesn't matter; no coefficient will use this
	return (float)(EnergyLog/EnergySum * 0xF3FCE0F5ull) * (0x1.0p-61f * 2.0f);
}

/**************************************/

static inline void Block_Transform_CalculatePsychoacoustics(float *MaskingNp, const uint32_t *Energy, const uint32_t *EnergyNp, int BlockSize, const int8_t *DecimationPattern) {
	struct Block_Transform_MaskingState_t MaskingState;

	//! Masking is performed for each sub-block, independent of each other
	const float *MaskingNpEnd = MaskingNp + BlockSize;
	int SubBlockIdx;
	for(SubBlockIdx=0;(MaskingNp < MaskingNpEnd);SubBlockIdx++) {
		int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];
		Block_Transform_MaskingState_Init(&MaskingState, Energy, EnergyNp, SubBlockSize);

		//! Perform simultaneous masking for each band
		int Band;
		for(Band=0;Band<SubBlockSize;Band++) MaskingNp[Band] = Block_Transform_GetMaskedLevel(&MaskingState, Energy, EnergyNp, Band, SubBlockSize);

		//! Advance buffers
		MaskingNp += SubBlockSize;
		Energy    += SubBlockSize;
		EnergyNp  += SubBlockSize;
	}
}

/**************************************/
//! EOF
/**************************************/
