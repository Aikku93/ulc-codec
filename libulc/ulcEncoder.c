/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
/**************************************/
#include "Fourier.h"
#include "ulcEncoder.h"
/**************************************/
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_BlockTransform.h"
#include "ulcEncoder_Encode.h"
#include "ulcEncoder_Helper.h"
/**************************************/
#if defined(__AVX__)
# define BUFFER_ALIGNMENT 32u //! __mm256
#elif defined(__SSE__)
# define BUFFER_ALIGNMENT 16u //! __mm128
#else
# define BUFFER_ALIGNMENT 4u //! float
#endif
/**************************************/

#define MIN_CHANS     1
#define MIN_BANDS    64 //! Mostly depends on the SIMD routines (currently limited by Block_Transform_CopySamples())
#define MAX_CHANS 65535
#define MAX_BANDS 65535
#define MIN_OVERLAP  16 //! Depends on SIMD routines; setting as 16 arbitrarily

//! Initialize encoder state
int ULC_EncoderState_Init(struct ULC_EncoderState_t *State) {
	//! Clear anything that is needed for EncoderState_Destroy()
	State->BufferData = NULL;

	//! Verify parameters
	int nChan      = State->nChan;
	int BlockSize  = State->BlockSize;
	int MinOverlap = State->MinOverlap;
	int MaxOverlap = State->MaxOverlap;
	if(nChan     < MIN_CHANS || nChan     > MAX_CHANS) return -1;
	if(BlockSize < MIN_BANDS || BlockSize > MAX_BANDS) return -1;
	if((BlockSize & (-BlockSize)) != BlockSize)        return -1;
	if(MinOverlap < MIN_OVERLAP) MinOverlap = State->MinOverlap = MIN_OVERLAP;
	if(MinOverlap > BlockSize)                         return -1;
	if((MinOverlap & (-MinOverlap)) != MinOverlap)     return -1;
	if(MaxOverlap < MIN_OVERLAP) MaxOverlap = State->MaxOverlap = MIN_OVERLAP;
	if(MaxOverlap > BlockSize)                         return -1;
	if((MaxOverlap & (-MaxOverlap)) != MaxOverlap)     return -1;

	//! Get buffer offsets+sizes
	//! PONDER: This... is probably not ideal
	int TransformBuffer_Size  = sizeof(float)                * (nChan* BlockSize   );
	int TransformNepers_Size  = sizeof(float)                * (nChan* BlockSize   );
	int TransformFwdLap_Size  = sizeof(float)                * (nChan*(BlockSize/2));
	int TransformTemp_Size    = sizeof(float)                * (       BlockSize   );
	int AnalysisKeys_Size     = sizeof(struct AnalysisKey_t) * (nChan* BlockSize   );
	int QuantsSum_Size        = sizeof(float)                * (nChan* ULC_MAX_QBANDS);
	int QuantsWeight_Size     = sizeof(float)                * (nChan* ULC_MAX_QBANDS);
	int Quants_Size           = sizeof(float)                * (nChan* ULC_MAX_QBANDS);
	int QuantsBw_Size         = sizeof(uint16_t)             * (nChan* ULC_MAX_QBANDS);
	int _TransformBuffer_Size = sizeof(float*)               * (nChan              );
	int _TransformNepers_Size = sizeof(float*)               * (nChan              );
	int _TransformFwdLap_Size = sizeof(float*)               * (nChan              );
	int _QuantsSum_Size       = sizeof(float*)               * (nChan              );
	int _QuantsWeight_Size    = sizeof(float*)               * (nChan              );
	int _Quants_Size          = sizeof(float*)               * (nChan              );
	int _QuantsBw_Size        = sizeof(uint16_t*)            * (nChan              );
	int TransformBuffer_Offs  = 0;
	int TransformNepers_Offs  = TransformBuffer_Offs    + TransformBuffer_Size;
	int TransformFwdLap_Offs  = TransformNepers_Offs    + TransformNepers_Size;
	int TransformTemp_Offs    = TransformFwdLap_Offs    + TransformFwdLap_Size;
	int AnalysisKeys_Offs     = TransformTemp_Offs      + TransformTemp_Size;
	int QuantsSum_Offs        = AnalysisKeys_Offs       + AnalysisKeys_Size;
	int QuantsWeight_Offs     = QuantsSum_Offs          + QuantsSum_Size;
	int Quants_Offs           = QuantsWeight_Offs       + QuantsWeight_Size;
	int QuantsBw_Offs         = Quants_Offs             + Quants_Size;
	int _TransformBuffer_Offs = QuantsBw_Offs           + QuantsBw_Size;
	int _TransformNepers_Offs = _TransformBuffer_Offs   + _TransformBuffer_Size;
	int _TransformFwdLap_Offs = _TransformNepers_Offs   + _TransformNepers_Size;
	int _QuantsSum_Offs       = _TransformFwdLap_Offs   + _TransformFwdLap_Size;
	int _QuantsWeight_Offs    = _QuantsSum_Offs         + _QuantsSum_Size;
	int _Quants_Offs          = _QuantsWeight_Offs      + _QuantsWeight_Size;
	int _QuantsBw_Offs        = _Quants_Offs            + _Quants_Size;
	int AllocSize = _QuantsBw_Offs + _QuantsBw_Size;

	//! Allocate buffer space
	char *Buf = State->BufferData = malloc(BUFFER_ALIGNMENT-1 + AllocSize);
	if(!Buf) return -1;

	//! Set initial state
	State->BitBudget   = 0.0;
	State->CoefBitRate = -1.0f; //! Will be set on first run
	State->LastTrackedRMSNp = -1.0e30;

	//! Initialize pointers
	int i, Chan;
	Buf += (-(uintptr_t)Buf) & (BUFFER_ALIGNMENT-1);
	State->TransformBuffer = (float   **)(Buf + _TransformBuffer_Offs);
	State->TransformNepers = (float   **)(Buf + _TransformNepers_Offs);
	State->TransformFwdLap = (float   **)(Buf + _TransformFwdLap_Offs);
	State->TransformTemp   = (float    *)(Buf + TransformTemp_Offs);
	State->AnalysisKeys    = (void     *)(Buf + AnalysisKeys_Offs);
	State->QuantsSum       = (float   **)(Buf + _QuantsSum_Offs);
	State->QuantsWeight    = (float   **)(Buf + _QuantsWeight_Offs);
	State->Quants          = (float   **)(Buf + _Quants_Offs);
	State->QuantsBw        = (uint16_t**)(Buf + _QuantsBw_Offs);
	for(Chan=0;Chan<nChan;Chan++) {
		State->TransformBuffer[Chan] = (float   *)(Buf + TransformBuffer_Offs  ) + Chan*BlockSize;
		State->TransformNepers[Chan] = (float   *)(Buf + TransformNepers_Offs  ) + Chan*BlockSize;
		State->TransformFwdLap[Chan] = (float   *)(Buf + TransformFwdLap_Offs  ) + Chan*(BlockSize/2);
		State->QuantsSum      [Chan] = (float   *)(Buf + QuantsSum_Offs        ) + Chan*ULC_MAX_QBANDS;
		State->QuantsWeight   [Chan] = (float   *)(Buf + QuantsWeight_Offs     ) + Chan*ULC_MAX_QBANDS;
		State->Quants         [Chan] = (float   *)(Buf + Quants_Offs           ) + Chan*ULC_MAX_QBANDS;
		State->QuantsBw       [Chan] = (uint16_t*)(Buf + QuantsBw_Offs         ) + Chan*ULC_MAX_QBANDS;

		//! Everything can remain uninitialized except for the lapping buffer
		for(i=0;i<BlockSize/2;i++) State->TransformFwdLap[Chan][i] = 0.0f;
	}

	//! Success
	return 1;
}

/**************************************/

//! Destroy encoder state
void ULC_EncoderState_Destroy(struct ULC_EncoderState_t *State) {
	//! Free buffer space
	free(State->BufferData);
}

/**************************************/

//! Encode block
int ULC_EncodeBlock(struct ULC_EncoderState_t *State, uint8_t *DstBuffer, const float *SrcData, float RateKbps, float PowerDecay) {
	//! Refill bit budget
	float AvgBitBudget = RateKbps*1000.0f/State->RateHz * State->BlockSize;
	State->BitBudget += AvgBitBudget;

	//! Transform input, build keys, and get maximum number of non-zero bands
	int nKeys = Block_Transform(State, SrcData, RateKbps, PowerDecay);
	int nNzMax; {
		//! First run?
		if(State->CoefBitRate < 0.0f) {
			//! Experimentally-derived average bitrate
			float MaxRate = MaxCodingKbps(State->BlockSize, State->nChan, State->RateHz);
			float CoefRate = (State->RateHz / 1000.0f / SQR(MaxRate)) * (MaxRate + RateKbps);
			State->CoefBitRate = CoefRate < 1.0f ? CoefRate : 1.0f;
		}

		//! Limit to 1.25x target bitrate to avoid quality spikes,
		//! as these have a tendency to sound bad in context
		float MinBudget = AvgBitBudget*0.75f;
		float MaxBudget = AvgBitBudget*1.25f;

		//! Clip budget to limit, derive nNzMax from that
		float Budget = State->BitBudget;
		if(Budget < MinBudget) Budget = MinBudget;
		if(Budget > MaxBudget) Budget = MaxBudget;
		nNzMax = lrintf(Budget * State->CoefBitRate);

		//! Sometimes we may get less non-zero bands than we're
		//! able to code (eg. during silence), so just clip here
		if(nNzMax > nKeys) nNzMax = nKeys;
	}

	//! Encode block
	int nNzCoded;
	int BlockBits = Block_Encode(State, DstBuffer, nNzMax, nKeys, &nNzCoded);

	//! Update state
	State->BitBudget -= BlockBits;
	if(nNzCoded) {
		//! RateDecayFactor specifies how much of the coefficient
		//! bit rate of the previous block we should take into
		//! account for future blocks:
		//!  0.0: Use old cost exclusively
		//!  1.0: Use new cost exclusively
		//! With the following formulation, we use the new cost
		//! exclusively if we hit our nNzMax target exactly,
		//! but if we didn't hit this target, then we rely more
		//! on the last-known 'goal' coefficient rate.
		//! NOTE:
		//!  -nNzCoded can be greater than nNzMax, depending on
		//!   the statistics of the zero runs (since we might
		//!   decide to code a coefficient that we didn't think
		//!   to include in our calculations), so we need to
		//!   clip to 1.0 here (since we met our goal)
		//!  -Slightly biased towards the new cost regardless of
		//!   meeting the target; this allows some adaptation to
		//!   take place, hopefully giving us a more precise cost
		//!   measurement in future blocks
		float RateDecayFactor = 0.25f + nNzCoded*0.75f/nNzMax;
		if(RateDecayFactor > 1.0f) RateDecayFactor = 1.0f;

		float nZCoefCost = (float)nNzCoded / BlockBits;
		State->CoefBitRate = State->CoefBitRate*(1.0f-RateDecayFactor) + nZCoefCost*RateDecayFactor;
	}

	//! Return block size
	return BlockBits;
}

/**************************************/
//! EOF
/**************************************/
