/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stdint.h>
/**************************************/

//! Decoder state structure
//! NOTE:
//!  -The global state data must be set before calling ULC_EncoderState_Init()
//!  -{nChan, BlockSize} must not change after calling ULC_EncoderState_Init()
struct ULC_DecoderState_t {
	//! Global state
	int nChan;       //! Channels in encoding scheme
	int BlockSize;   //! Transform block size
	int OverlapSize; //! Cached overlap from last block

	//! Encoding state
	//! Buffer memory layout:
	//!  Data:
	//!   char    _Padding[];
	//!   float   TransformBuffer[BlockSize]
	//!   float   TransformTemp  [BlockSize]
	//!   float   TransformInvLap[nChan][BlockSize/2]
	//!  MD-array pointers:
	//!   float    *_TransformInvLap[nChan]
	//! BufferData contains the pointer returned by malloc()
	void   *BufferData;
	float  *TransformBuffer;
	float  *TransformTemp;
	float **TransformInvLap;
};

/**************************************/

//! Initialize decoder state
//! On success, returns a non-negative value
//! On failure, returns a negative value
int ULC_DecoderState_Init(struct ULC_DecoderState_t *State);

//! Destroy decoder state
void ULC_DecoderState_Destroy(struct ULC_DecoderState_t *State);

/**************************************/

//! Decode block
//! NOTE:
//!  -Output data will have its channels arranged sequentially;
//!   For example:
//!   {
//!    0,1,2,3...BlockSize-1, //! Chan0
//!    0,1,2,3...BlockSize-1, //! Chan1
//!   }
//! Returns the number of bits read.
int ULC_DecodeBlock(struct ULC_DecoderState_t *State, float *DstData, const uint8_t *SrcBuffer);

/**************************************/
//! EOF
/**************************************/
