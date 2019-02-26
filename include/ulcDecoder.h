/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
#include <stdint.h>
/**************************************/

//! Decoder state structure
//! NOTE:
//!  -The global state data must be set before calling ULC_EncoderState_Init()
//!  -{nChan, BlockSize} must not change after calling ULC_EncoderState_Init()
//!  -{nQuants,QuantsBw} can be changed on a block-to-block basis
//! Example for QuantsBw (BlockSize=512, nQuants=8):
//!  {4,8,12,24,48,96,144,176}
//! NOTE: UNDEFINED BEHAVIOUR MAY OCCUR IF THE TOTAL OF THE
//!       QUANTIZER BANDWIDTHS IS NOT EQUAL TO THE BLOCK SIZE.
struct ULC_DecoderState_t {
	//! Global state
	      size_t    nChan;     //! Channels in encoding scheme
	      size_t    BlockSize; //! Transform block size
	      size_t    nQuants;   //! Number of quantizer bands
	const uint16_t *QuantsBw;  //! Quantizer band bandwidths

	//! Encoding state
	//! Buffer memory layout:
	//!  Data:
	//!   char    _Padding[];
	//!   float   TransformBuffer[BlockSize]
	//!   float   TransformTemp  [BlockSize]
	//!   float   TransformInvLap[nChan][BlockSize/2]
	//!   uint8_t Quants         [nQuants]
	//!  MD-array pointers:
	//!   float    *_TransformInvLap[nChan]
	//! BufferData contains the pointer returned by malloc()
	void *BufferData;
	float   *TransformBuffer;
	float   *TransformTemp;
	float  **TransformInvLap;
	uint8_t *Quants;
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
//!    (4*nQuant + 4*BlockSize)*nChan
//!   So input buffer should generally contain at least
//!   that much cached.
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
