/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
#include <stdint.h>
/**************************************/

//! 0 == No psychoacoustic optimizations
//! 1 == Use psychoacoustic model
#define ULC_USE_PSYHOACOUSTICS 1

/**************************************/

//! Encoder state structure
//! NOTE:
//!  -The global state data must be set before calling ULC_EncoderState_Init()
//!  -{RateHz, nChan, BlockSize, BlockOverlap} must not change after calling ULC_EncoderState_Init()
//!  -MAX_QUANTS is an internal value, no larger than 48
struct ULC_EncoderState_t {
	//! Global state
	size_t RateHz;       //! Playback rate (used for rate control)
	size_t nChan;        //! Channels in encoding scheme
	size_t BlockSize;    //! Transform block size
	size_t BlockOverlap; //! Block overlap

	//! Rate control state
	double BitBudget;   //! Bit budget left over from previous block (similar to a bit reservoir)
	double CoefBitRate; //! Reciprocal of average bits per non-zero coefficient

	//! Encoding state
	//! Buffer memory layout:
	//!  Data:
	//!   char          _Padding[];
	//!   float         TransformBuffer  [nChan][BlockSize]
	//!   float         TransformTemp    [2*BlockSize]
	//!   float         TransformFwdLap  [nChan][BlockSize/2]
	//!   float         TransformFlatness[nChan][N_FLATNESS]
	//!   AnalysisKey_t AnalysisKeys     [nChan*BlockSize]
	//!   double        QuantsPow        [nChan][MAX_QUANTS]
	//!   double        QuantsAbs        [nChan][MAX_QUANTS]
	//!   uint16_t      QuantsBw         [nChan][MAX_QUANTS]
	//!   int16_t       Quants           [nChan][MAX_QUANTS]
	//!  Followed by MD-array pointers:
	//!   float    *_TransformBuffer  [nChan]
	//!   float    *_TransformFwdLap  [nChan]
	//!   float    *_TransformFlatness[nChan]
	//!   double   *_QuantsPow        [nChan]
	//!   double   *_QuantsAbs        [nChan]
	//!   uint16_t *_QuantsBw         [nChan]
	//!   int16_t  *_Quants           [nChan]
	//! BufferData contains the pointer returned by malloc()
	void *BufferData;
	float    **TransformBuffer;
	float     *TransformTemp;
	float    **TransformFwdLap;
	float    **TransformFlatness;
	void      *AnalysisKeys;
	double   **QuantsPow;
	double   **QuantsAbs;
	uint16_t **QuantsBw;
	int16_t  **Quants;
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
//!    (4 + 12*(MAX_QUANTS-1) + 4*BlockSize)*nChan
//!   So output buffer size must be at least that size
//!  -Input data must have its channels arranged sequentially;
//!   For example:
//!   {
//!    0,1,2,3...BlockSize-1, //! Chan0
//!    0,1,2,3...BlockSize-1, //! Chan1
//!   }
//!  -PowerDecay is meant for M/S transform and the like, where
//!   channels are increasingly 'less important'. Pass 1.0 if
//!   all channels are equally important.
//! Returns the block size in bits
size_t ULC_EncodeBlock(struct ULC_EncoderState_t *State, uint8_t *DstBuffer, const float *SrcData, double RateKbps, float PowerDecay);

/**************************************/
//! EOF
/**************************************/
