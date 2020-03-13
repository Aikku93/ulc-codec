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

//! Lowest possible coefficient value
#define ULC_COEF_EPS (0x1.0p-33) //! 4+0xE+15 = Maximum extended-precision quantizer

//! Used in Neper-scale coefficients
//! dB calculations would add computational cost for the exact same results,
//! as logf() is faster than log2f() which is faster than log10f()... somehow
#define ULC_COEF_NEPER_OUT_OF_RANGE 0.0f

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
	//!   float         TransformBuffer[nChan][BlockSize]
	//!   float         TransformNepers[nChan][BlockSize]
	//!   float         TransformFwdLap[nChan][BlockSize/2]
	//!   float         TransformTemp  [BlockSize]
	//!   AnalysisKey_t AnalysisKeys   [nChan*BlockSize]
	//!   double        QuantsPow      [nChan][MAX_QUANTS]
	//!   double        QuantsAbs      [nChan][MAX_QUANTS]
	//!   float         Quants         [nChan][MAX_QUANTS]
	//!   uint16_t      QuantsBw       [nChan][MAX_QUANTS]
	//!  Followed by MD-array pointers:
	//!   float    *_TransformBuffer  [nChan]
	//!   float    *_TransformNepers  [nChan]
	//!   float    *_TransformFwdLap  [nChan]
	//!   double   *_QuantsPow        [nChan]
	//!   double   *_QuantsAbs        [nChan]
	//!   float    *_Quants           [nChan]
	//!   uint16_t *_QuantsBw         [nChan]
	//! BufferData contains the pointer returned by malloc()
	void *BufferData;
	float    **TransformBuffer;
	float    **TransformNepers;
	float    **TransformFwdLap;
	float     *TransformTemp;
	void      *AnalysisKeys;
	double   **QuantsPow;
	double   **QuantsAbs;
	float    **Quants;
	uint16_t **QuantsBw;
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
