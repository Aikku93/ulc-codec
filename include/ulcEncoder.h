/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stdint.h>
/**************************************/

//! 0 == No psychoacoustic optimizations
//! 1 == Use psychoacoustic model
#define ULC_USE_PSYCHOACOUSTICS 1

//! 0 == No noise-fill coding
//! 1 == Use noise-fill where useful
#define ULC_USE_NOISE_CODING 1

//! 0 == No window switching
//! 1 == Use window switching
#define ULC_USE_WINDOW_SWITCHING 1

//! Maximum number of subblocks present in a block
#define ULC_MAX_SUBBLOCKS 4

//! Lowest possible coefficient value
#define ULC_COEF_EPS (0x1.0p-31f) //! 4+0xE+0xC = Maximum extended-precision quantizer

//! Used in Neper-scale coefficients
//! dB calculations would add computational cost for the exact same results,
//! as logf() is faster than log2f() which is faster than log10f()... somehow.
//! This value is set to 0.0, as the only times that out-of-range coefficients
//! are used are during MDCT/MDST (and psychoacoustics) calculations, where
//! these log-domain values are used as part of a weighted geometric mean
#define ULC_COEF_NEPER_OUT_OF_RANGE 0.0f

/**************************************/

//! Encoder state structure
//! NOTE:
//!  -The global state data must be set before calling ULC_EncoderState_Init()
//!  -{RateHz, nChan, BlockSize, BlockOverlap} must not change after calling ULC_EncoderState_Init()
struct ULC_EncoderState_t {
	//! Global state
	int RateHz;     //! Playback rate (used for rate control)
	int nChan;      //! Channels in encoding scheme
	int BlockSize;  //! Transform block size

	//! Encoding state
	//! Buffer memory layout:
	//!  Data:
	//!   char  _Padding[];
	//!   float SampleBuffer   [nChan*BlockSize]
	//!   float TransformBuffer[nChan*BlockSize]
	//!   float TransformNoise [nChan*BlockSize] <- With ULC_USE_NOISE_CODING only
	//!   float TransformFwdLap[nChan*BlockSize/2]
	//!   float TransformTemp  [MAX(2,nChan)*BlockSize]
	//!   int   TransformIndex [nChan*BlockSize]
	//! BufferData contains the original pointer returned by malloc()
	int    WindowCtrl;        //! Window control parameter (for last coded block)
	int    NextWindowCtrl;    //! Window control parameter (for data in SampleBuffer)
	float  BlockComplexity;   //! Coefficient distribution complexity (0 = Highly tonal, 1 = Highly noisy)
	void  *BufferData;
	float *SampleBuffer;
	float *TransformBuffer;
#if ULC_USE_NOISE_CODING
	float *TransformNoise;
#endif
	float *TransformFwdLap;
	float *TransformTemp;
	int   *TransformIndex;
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
//!    8 + nChan*(8+4 + (16+4)*(BlockSize-1))
//!     8    = Window shape[s] selection
//!     8+4  = Initial quantizer ([8h,0h,]Eh,Xh) and first coefficient (Xh)
//!     16+4 = Quantizer (8h,0h,Eh,Xh) + coefficient (Xh)
//!   So output buffer size should be at least that size
//!  -Input data must have its channels arranged sequentially;
//!   For example:
//!   {
//!    0,1,2,3...BlockSize-1, //! Chan0
//!    0,1,2,3...BlockSize-1, //! Chan1
//!   }
//! Notes regarding coding modes:
//!  -CBR will encode as many coefficients as possible for a given
//!   RateKbps, never exceeding this value (for example, if 128.0kbps
//!   is the target and the encoding choice is between 127.0kbps and
//!   128.01kbps, 127.0kbps will always be chosen).
//!   The rate is matched via binary search, and so this encoding
//!   mode (and ABR, which uses the same mechanism) is the slowest.
//!  -ABR mode tries to balance the number of coefficients in each
//!   block based on their complexity. It will achieve an average
//!   bitrate very close to the target, but may be slightly off due
//!   to rounding errors.
//!   To run in ABR mode, the file must first be analyzed to get the
//!   average complexity, and then this is passed to the routine
//!   alongside the desired RateKbps (note that AvgComplexity can be
//!   passed arbitrarily without a pre-pass, but the target bitrate
//!   might not be achieved).
//!   This encoding mode is just as slow as CBR mode, as the algorithm
//!   chooses a target bitrate for each block, which must go through
//!   the CBR encoding routine for each block.
//!  -VBR mode attempts to estimate the average complexity of the file
//!   based on a Quality parameter, and then proceeds to encode the
//!   file based on this metric alone. This mode is the fastest, but
//!   there is no guarantee for the bitrate.
//!   Generally, though:
//!     0 < Quality <= 10 = Average  <30kbps
//!    10 < Quality <= 20 = Average  <40kbps
//!    20 < Quality <= 30 = Average  <50kbps
//!    30 < Quality <= 40 = Average  <60kbps
//!    40 < Quality <= 50 = Average  <75kbps
//!    50 < Quality <= 60 = Average  <95kbps
//!    60 < Quality <= 70 = Average <125kbps
//!    70 < Quality <= 80 = Average <175kbps
//!    80 < Quality <= 90 = Average <300kbps
//! Returns the block size in bits
int ULC_EncodeBlock_CBR(struct ULC_EncoderState_t *State, uint8_t *DstBuffer, const float *SrcData, float RateKbps);
int ULC_EncodeBlock_ABR(struct ULC_EncoderState_t *State, uint8_t *DstBuffer, const float *SrcData, float RateKbps, float AvgComplexity);
int ULC_EncodeBlock_VBR(struct ULC_EncoderState_t *State, uint8_t *DstBuffer, const float *SrcData, float Quality);

/**************************************/
//! EOF
/**************************************/
