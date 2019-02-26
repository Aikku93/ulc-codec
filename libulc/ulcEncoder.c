/**************************************/
#include <stddef.h>
#include <stdint.h>
#include <math.h>
/**************************************/
#include "Fourier.h"
#include "ulcEncoder.h"
#include "ulcUtility.h"
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

//! Initialize encoder state
int ULC_EncoderState_Init(struct ULC_EncoderState_t *State) {
	//! Clear anything that is needed for EncoderState_Destroy()
	State->BufferData = NULL;

	//! Verify parameters
	size_t nChan     = State->nChan;
	size_t BlockSize = State->BlockSize;
	size_t nQuants   = State->nQuants;
	if(nChan     < MIN_CHANS || nChan     > MAX_CHANS) return -1;
	if(BlockSize < MIN_BANDS || BlockSize > MAX_BANDS) return -1;

	//! Get buffer offsets+sizes
	//! PONDER: This... is probably not ideal
	size_t TransformBuffer_Size  = sizeof(float)                * (nChan* BlockSize   );
	size_t TransformTemp_Size    = sizeof(float)                * (2    * BlockSize   );
	size_t TransformFwdLap_Size  = sizeof(float)                * (nChan*(BlockSize/2));
	size_t AnalysisKeys_Size     = sizeof(struct AnalysisKey_t) * (nChan* BlockSize   );
	size_t Quants_Size           = sizeof(int16_t)              * (nChan* nQuants     );
	size_t QuantsCnt_Size        = sizeof(uint16_t)             * (nChan* nQuants     );
	size_t QuantsPow_Size        = sizeof(double)               * (nChan* nQuants     );
	size_t _TransformBuffer_Size = sizeof(float*)               * (nChan              );
	size_t _TransformFwdLap_Size = sizeof(float*)               * (nChan              );
	size_t _Quants_Size          = sizeof(int16_t*)             * (nChan              );
	size_t _QuantsCnt_Size       = sizeof(uint16_t*)            * (nChan              );
	size_t _QuantsPow_Size       = sizeof(double*)              * (nChan              );
	size_t TransformBuffer_Offs  = 0;
	size_t TransformTemp_Offs    = TransformBuffer_Offs  + TransformBuffer_Size;
	size_t TransformFwdLap_Offs  = TransformTemp_Offs    + TransformTemp_Size;
	size_t AnalysisKeys_Offs     = TransformFwdLap_Offs  + TransformFwdLap_Size;
	size_t Quants_Offs           = AnalysisKeys_Offs     + AnalysisKeys_Size;
	size_t QuantsCnt_Offs        = Quants_Offs           + Quants_Size;
	size_t QuantsPow_Offs        = QuantsCnt_Offs        + QuantsCnt_Size;
	size_t _TransformBuffer_Offs = QuantsPow_Offs        + QuantsPow_Size;
	size_t _TransformFwdLap_Offs = _TransformBuffer_Offs + _TransformBuffer_Size;
	size_t _Quants_Offs          = _TransformFwdLap_Offs + _TransformFwdLap_Size;
	size_t _QuantsCnt_Offs       = _Quants_Offs          + _Quants_Size;
	size_t _QuantsPow_Offs       = _QuantsCnt_Offs       + _QuantsCnt_Size;
	size_t AllocSize = _QuantsPow_Offs + _QuantsPow_Size;

	//! Allocate buffer space
	char *Buf = State->BufferData = malloc(BUFFER_ALIGNMENT-1 + AllocSize);
	if(!Buf) return -1;

	//! Set initial state
	State->BitBudget   = 0.0;
	State->CoefBitRate = 1.0 / 7.5; //! Initial approximation (reciprocal of bits per nZ coefficient)

	//! Initialize pointers
	size_t i, Chan;
	Buf += -(uintptr_t)Buf % BUFFER_ALIGNMENT;
	State->TransformBuffer = (float   **)(Buf + _TransformBuffer_Offs);
	State->TransformTemp   = (float    *)(Buf + TransformTemp_Offs);
	State->TransformFwdLap = (float   **)(Buf + _TransformFwdLap_Offs);
	State->AnalysisKeys    = (void     *)(Buf + AnalysisKeys_Offs);
	State->Quants          = (int16_t **)(Buf + _Quants_Offs);
	State->QuantsCnt       = (uint16_t**)(Buf + _QuantsCnt_Offs);
	State->QuantsPow       = (double  **)(Buf + _QuantsPow_Offs);
	for(Chan=0;Chan<nChan;Chan++) {
		State->TransformBuffer[Chan] = (float   *)(Buf + TransformBuffer_Offs) + Chan*BlockSize;
		State->TransformFwdLap[Chan] = (float   *)(Buf + TransformFwdLap_Offs) + Chan*(BlockSize/2);
		State->Quants         [Chan] = (int16_t *)(Buf + Quants_Offs         ) + Chan*nQuants;
		State->QuantsCnt      [Chan] = (uint16_t*)(Buf + QuantsCnt_Offs      ) + Chan*nQuants;
		State->QuantsPow      [Chan] = (double  *)(Buf + QuantsPow_Offs      ) + Chan*nQuants;

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
size_t ULC_EncodeBlock(struct ULC_EncoderState_t *State, uint8_t *DstBuffer, const float *SrcData, double RateKbps) {
	//! Refill bit budget
	double AvgBitBudget = RateKbps*1000.0/State->RateHz * State->BlockSize;
	State->BitBudget += AvgBitBudget;

	//! Transform input, build keys, and get maximum number of non-zero bands
	size_t nKeys = Block_Transform(State, SrcData);
	size_t nNzMax; {
		//! Limit to 1.25x target bitrate to avoid quality spikes,
		//! as these have a tendency to sound bad in context
		double MinBudget = AvgBitBudget*0.75;
		double MaxBudget = AvgBitBudget*1.25;

		//! Clip budget to limit, derive nNzMax from that
		double Budget = State->BitBudget;
		if(Budget < MinBudget) Budget = MinBudget;
		if(Budget > MaxBudget) Budget = MaxBudget;
		nNzMax = round(Budget * State->CoefBitRate);

		//! Sometimes we may get less non-zero bands than we're
		//! able to code (eg. during silence), so just clip here
		if(nNzMax > nKeys) nNzMax = nKeys;
	}

	//! Encode block
	size_t nNzCoded;
	size_t BlockBits = Block_Encode(State, DstBuffer, nNzMax, nKeys, &nNzCoded);

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
		double RateDecayFactor = 0.25 + nNzCoded*0.75/nNzMax;
		if(RateDecayFactor > 1.0) RateDecayFactor = 1.0;

		double nZCoefCost = (double)nNzCoded / BlockBits;
		State->CoefBitRate = State->CoefBitRate*(1.0-RateDecayFactor) + nZCoefCost*RateDecayFactor;
	}

	//! Return block size
	return BlockBits;
}

/**************************************/
//! EOF
/**************************************/
