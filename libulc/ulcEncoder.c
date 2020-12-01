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
#include "ulcHelper.h"
/**************************************/
#include "ulcEncoder_BlockTransform.h"
#include "ulcEncoder_Encode.h"
/**************************************/
#define BUFFER_ALIGNMENT 64u //! Always align memory to 64-byte boundaries (preparation for AVX-512)
/**************************************/

#define MIN_CHANS    1
#define MAX_CHANS  255
#define MIN_BANDS   64 //! 64-point MDCT is pretty extreme, so setting this as the limit
#define MAX_BANDS 8192
#define MIN_OVERLAP 16 //! Depends on SIMD routines; setting as 16 arbitrarily

//! Initialize encoder state
int ULC_EncoderState_Init(struct ULC_EncoderState_t *State) {
	//! Clear anything that is needed for EncoderState_Destroy()
	State->BufferData = NULL;

	//! Verify parameters
	int nChan      = State->nChan;
	int BlockSize  = State->BlockSize;
	if(nChan     < MIN_CHANS || nChan     > MAX_CHANS) return -1;
	if(BlockSize < MIN_BANDS || BlockSize > MAX_BANDS) return -1;
	if((BlockSize & (-BlockSize)) != BlockSize)        return -1;

	//! Get buffer offsets+sizes
	//! PONDER: This... is probably not ideal
	//! NOTE: Psychoacoustics needs at least two BlockSize temporary buffers
	//! to store the energy information (raw/linear and log-domain)
	int SampleBuffer_Size    = sizeof(float) * (nChan* BlockSize   );
	int TransformBuffer_Size = sizeof(float) * (nChan* BlockSize   );
	int TransformNepers_Size = sizeof(float) * (nChan* BlockSize   );
	int TransientEnergy_Size = sizeof(float) * (       BlockSize   );
	int TransformFwdLap_Size = sizeof(float) * (nChan*(BlockSize/2));
	int TransformTemp_Size   = sizeof(float) * ((nChan + (nChan < 2)) * BlockSize);
	int TransformIndex_Size  = sizeof(int)   * (nChan* BlockSize   );
	int SampleBuffer_Offs    = 0;
	int TransformBuffer_Offs = SampleBuffer_Offs    + SampleBuffer_Size;
	int TransformNepers_Offs = TransformBuffer_Offs + TransformBuffer_Size;
	int TransientEnergy_Offs = TransformNepers_Offs + TransformNepers_Size;
	int TransformFwdLap_Offs = TransientEnergy_Offs + TransientEnergy_Size;
	int TransformTemp_Offs   = TransformFwdLap_Offs + TransformFwdLap_Size;
	int TransformIndex_Offs  = TransformTemp_Offs   + TransformTemp_Size;
	int AllocSize            = TransformIndex_Offs  + TransformIndex_Size;

	//! Allocate buffer space
	char *Buf = State->BufferData = malloc(BUFFER_ALIGNMENT-1 + AllocSize);
	if(!Buf) return -1;

	//! Initialize pointers
	Buf += (-(uintptr_t)Buf) & (BUFFER_ALIGNMENT-1);
	State->SampleBuffer    = (float*)(Buf + SampleBuffer_Offs);
	State->TransformBuffer = (float*)(Buf + TransformBuffer_Offs);
	State->TransformNepers = (float*)(Buf + TransformNepers_Offs);
	State->TransientEnergy = (float*)(Buf + TransientEnergy_Offs);
	State->TransformFwdLap = (float*)(Buf + TransformFwdLap_Offs);
	State->TransformTemp   = (float*)(Buf + TransformTemp_Offs);
	State->TransformIndex  = (int  *)(Buf + TransformIndex_Offs);

	//! Set initial state
	int i;
	State->NextWindowCtrl = 0x10; //! No decimation, full overlap
	for(i=0;i<nChan*(BlockSize  );i++) State->SampleBuffer   [i] = 0.0f;
	for(i=0;i<      (BlockSize  );i++) State->TransientEnergy[i] = 0.0f;
	for(i=0;i<nChan*(BlockSize/2);i++) State->TransformFwdLap[i] = 0.0f;

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

//! Encode block (CBR mode)
int ULC_EncodeBlock_CBR(struct ULC_EncoderState_t *State, uint8_t *DstBuffer, const float *SrcData, float RateKbps) {
	int Size;
	int nOutCoef  = -1;
	int BitBudget = (int)((State->BlockSize * RateKbps) * 1000.0f/State->RateHz); //! NOTE: Truncate

	//! Perform a binary search for the optimal nOutCoef
	int Lo = 0, Hi = Block_Transform(State, SrcData);
	if(Lo < Hi) {
		do {
			nOutCoef = (Lo + Hi) / 2u;
			Size = Block_Encode_EncodePass(State, DstBuffer, nOutCoef);
			     if(Size < BitBudget) Lo = nOutCoef;
			else if(Size > BitBudget) Hi = nOutCoef-1;
			else {
				//! Should very, VERY rarely happen, but just in case
				Lo = nOutCoef;
				break;
			}
		} while(Lo < Hi-1);
	}

	//! Avoid going over budget
	int nOutCoefFinal = Lo;
	if(nOutCoefFinal != nOutCoef) Size = Block_Encode_EncodePass(State, DstBuffer, nOutCoef = nOutCoefFinal);
	return Size;
}

/**************************************/

//! Encode block (VBR mode)
int ULC_EncodeBlock_VBR(struct ULC_EncoderState_t *State, uint8_t *DstBuffer, const float *SrcData, float Quality) {
	int nOutCoef = (int)(Quality*(1.0f/100.0f) * State->nChan * State->BlockSize + 0x1.FFFFFFp-1f); //! NOTE: Ceiling
	int MaxCoef  = Block_Transform(State, SrcData);
	if(nOutCoef > MaxCoef) nOutCoef = MaxCoef;
	return Block_Encode_EncodePass(State, DstBuffer, nOutCoef);
}

/**************************************/
//! EOF
/**************************************/
