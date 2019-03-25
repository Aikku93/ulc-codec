/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
#include <stdint.h>
/**************************************/

//! Decoder state structure
//! NOTE:
//!  -The global state data must be set before calling ULC_EncoderState_Init()
//!  -{nChan, BlockSize, BlockOverlap} must not change after calling ULC_EncoderState_Init()
struct ULC_DecoderState_t {
	//! Global state
	size_t nChan;        //! Channels in encoding scheme
	size_t BlockSize;    //! Transform block size
	size_t BlockOverlap; //! Block overlap

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
	void *BufferData;
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
//!  -Maximum size (in bits) for each block is:
//!    (4 + 12*(MAX_QUANTS-1) + 4*BlockSize)*nChan
//!   So input buffer should generally contain at least
//!   that much cached.
//!  -MAX_QUANTS is an internal value from the encoder, no larger than 48
//!  -Output data will have its channels arranged sequentially;
//!   For example:
//!   {
//!    0,1,2,3...BlockSize-1, //! Chan0
//!    0,1,2,3...BlockSize-1, //! Chan1
//!   }
//! Returns the number of bits read
size_t ULC_DecodeBlock(const struct ULC_DecoderState_t *State, float *DstData, const uint8_t *SrcBuffer);

/**************************************/
//! EOF
/**************************************/
