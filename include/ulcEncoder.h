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

//! Encoder state structure
//! NOTE:
//!  -The global state data must be set before calling ULC_EncoderState_Init()
//!  -{RateHz, nChan, BlockSize} must not change after calling ULC_EncoderState_Init()
//!  -{nQuants,QuantsBw} can be changed on a block-to-block basis
//! Example for QuantsBw (BlockSize=512, nQuants=8):
//!  {4,8,12,24,48,96,144,176}
//! NOTE: UNDEFINED BEHAVIOUR MAY OCCUR IF THE TOTAL OF THE
//!       QUANTIZER BANDWIDTHS IS NOT EQUAL TO THE BLOCK SIZE.
struct ULC_EncoderState_t {
	//! Global state
	      size_t    RateHz;    //! Playback rate (used for rate control)
	      size_t    nChan;     //! Channels in encoding scheme
	      size_t    BlockSize; //! Transform block size
	      size_t    nQuants;   //! Number of quantizer bands
	const uint16_t *QuantsBw;  //! Quantizer band bandwidths

	//! Rate control state
	double BitBudget;   //! Bit budget left over from previous block (similar to a bit reservoir)
	double CoefBitRate; //! Reciprocal of average bits per non-zero coefficient

	//! Encoding state
	//! Buffer memory layout:
	//!  Data:
	//!   char          _Padding[];
	//!   float         TransformBuffer[nChan][BlockSize]
	//!   float         TransformTemp  [2*BlockSize]
	//!   float         TransformFwdLap[nChan][BlockSize/2]
	//!   AnalysisKey_t AnalysisKeys   [nChan*BlockSize]
	//!   int16_t       Quants         [nChan][nQuants]
	//!   double        QuantsPow      [nChan][nQuants]
	//!   double        QuantsAvg      [nChan][nQuants]
	//!  MD-array pointers:
	//!   float    *_TransformBuffer[nChan]
	//!   float    *_TransformFwdLap[nChan]
	//!   int16_t  *_Quants         [nChan]
	//!   double   *_QuantsPow      [nChan]
	//!   double   *_QuantsAvg      [nChan]
	//! BufferData contains the pointer returned by malloc()
	void *BufferData;
	float    **TransformBuffer;
	float     *TransformTemp;
	float    **TransformFwdLap;
	void      *AnalysisKeys;
	int16_t  **Quants;
	double   **QuantsPow;
	double   **QuantsAvg;
};

/**************************************/

//! Initialize encoder state
//! On success, returns a non-negative value
//! On failure, returns a negative value
int ULC_EncoderState_Init(struct ULC_EncoderState_t *State);

//! Destroy encoder state
void ULC_EncoderState_Destroy(struct ULC_EncoderState_t *State);

/**************************************/

//! Encode block
//! NOTE:
//!  -Maximum size (in bits) for each block is:
//!    (4*nQuant + 4*BlockSize)*nChan
//!   So output buffer size must be at least that size
//!  -Input data must have its channels arranged sequentially;
//!   For example:
//!   {
//!    0,1,2,3...BlockSize-1, //! Chan0
//!    0,1,2,3...BlockSize-1, //! Chan1
//!   }
//! Returns the block size in bits
size_t ULC_EncodeBlock(struct ULC_EncoderState_t *State, uint8_t *DstBuffer, const float *SrcData, double RateKbps);

/**************************************/
//! EOF
/**************************************/
