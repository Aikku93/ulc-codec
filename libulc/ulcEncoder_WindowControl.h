/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <math.h>
/**************************************/
#include "ulcEncoder.h"
#include "ulcHelper.h"
/**************************************/

//! Get optimal log base-2 overlap and window scalings for transients
//! The idea is that if a transient is relatively centered with the
//! transition region of a subblock, then we can just set the overlap
//! amount to account for it and avoid reducing the window size too
//! much, preserving the quality gains of a larger transform. At the
//! same time, we also need to /make/ the transient sit within a
//! transition region to take advantage of this, and so we combine
//! the overlap scaling and window-switching strategies.
//! NOTE: Bit codes for transient region coding, and their window sizes:
//!  First nybble:
//!   0xxx: No decimation. xxx = Overlap scaling
//!   1xxx: Decimate. xxx = Overlap scaling for the transient subblock
//!  Second nybble (when first nybble is 1xxx; otherwise, this is implicitly 0001):
//!   1xxx: Decimation by 1/8: Position = 0~7
//!    1000: N/8*,N/8,N/4,N/2
//!    1001: N/8,N/8*,N/4,N/2
//!    1010: N/4,N/8*,N/8,N/2
//!    1011: N/4,N/8,N/8*,N/2
//!    1100: N/2,N/8*,N/8,N/4
//!    1101: N/2,N/8,N/8*,N/4
//!    1110: N/2,N/4,N/8*,N/8
//!    1111: N/2,N/4,N/8,N/8*
//!   01xx: Decimation by 1/4: Position = 0~3
//!    0100: N/4*,N/4,N/2
//!    0101: N/4,N/4*,N/2
//!    0110: N/2,N/4*,N/4
//!    0111: N/2,N/4,N/4*
//!   001x: Decimation by 1/2: Position = 0~1
//!    0010: N/2*,N/2
//!    0011: N/2,N/2*
//!   0001: No decimation (not coded in the bitstream)
//!    0001: N/1*
//!  The starred subblocks set overlap scaling with regards to the last [sub]block.
//! NOTE:
//!  -TransientBuffer[] must be ULC_MAX_BLOCK_DECIMATION_FACTOR*4 in size (LL,L,M,R at maximum decimation)
//!   and will be updated as (LL=M_Old, L=R_Old, M=M_New, R=R_New)
//!  -TmpBuffer[] must be BlockSize in size
//!  -SmoothingTaps[] must be 2 elements in size
#pragma GCC push_options
#pragma GCC optimize("fast-math") //! Should improve things, hopefully, maybe
static inline void Block_Transform_GetWindowCtrl_TransientFiltering(
	const float *ThisBlockData,
	const float *LastBlockData,
	      float *TransientBuffer,
	      float *TmpBuffer,
	      float *SmoothingTaps,
	      int    BlockSize,
	      int    nChan
) {
	int n, Chan;

	//! Perform a bandpass filter to isolate the energy that is
	//! important to transient detection. Generally, LF energy
	//! and HF energy are 'unimportant', and it's the MF energy
	//! that has most of the information we're interested in.
	//! Transfer function:
	//!  H(z) = z^1 - z^-1
	//! NOTE: This filter does not have unity gain, as doing so
	//! would add some multiplications that reduce performance.
	//! NOTE: We end up losing the last sample of the new block,
	//! but this shouldn't affect things; we instead compensate
	//! with a 1 sample delay to preserve accuracy.
	//! NOTE: DOFILTER() accepts z^-1,1,z^1 for flexibility if we
	//! ever need to change the filter formula.
	for(n=0;n<BlockSize;n++) TmpBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
#define DOFILTER(zM1, z0, z1) SQR((zM1) - (z1))
		float *Dst = TmpBuffer;
		const float *SrcOld = LastBlockData + Chan*BlockSize + BlockSize;
		const float *SrcNew = ThisBlockData + Chan*BlockSize;
		{
			*Dst++ += DOFILTER(SrcOld[-2], SrcOld[-1], SrcNew[ 0]); //! <- Belongs to the last block
			*Dst++ += DOFILTER(SrcOld[-1], SrcNew[ 0], SrcNew[+1]), SrcNew++;
		}
		for(n=1;n<BlockSize-1;n++) {
			*Dst++ += DOFILTER(SrcNew[-1], SrcNew[ 0], SrcNew[+1]), SrcNew++;
		}
#undef DOFILTER
	}

	//! Apply two lowpass filters and take their difference to
	//! obtain a low-frequency bandpass filter, and accumulate
	//! energy to the M/R bins while simultaneously restoring
	//! the LL/L bins from te cache and swapping in new data.
	//! Theory:
	//!  Transients result in pulses close to DC, so we try to
	//!  remove harmonic reflections in the higher freqs. We
	//!  then apply another filter to remove DC content, as this
	//!  causes biasing of the signal analysis.
	//! NOTE: It's important to keep the smoothing taps accurate,
	//! so we save it across blocks and buffer the filtered data.
	//! NOTE: We perform the filtering in a companded domain,
	//! as this emphasizes the transient structure far better.
	//! NOTE: Slightly refactored to remove a multiplication; a
	//! gain of 1/(1-LPDecay) is applied to LPTap, and a gain of
	//! (1/(1-DCDecay))/(1/(1-LPDecay))=(1-LPDecay)/(1-DCDecay)
	//! to DCTap, which then cancels out the normalization gain
	//! on DCTap, resulting in an overall gain of 1/(1-LPDecay)
	//! in the output, plus whatever gain we had as input.
	{
#define DOFILTER()                        \
	v      = sqrtf(sqrtf(*Src++)),    \
	LPTap *= LPDecay,                 \
	DCTap *= DCDecay,                 \
	LPTap += v,                       \
	DCTap += v * DCGainCompensation,  \
	v      = SQR(SQR(LPTap - DCTap))
		//! BlockSize(={LL,L}) + BlockSize(={M,R}) = BlockSize*2 = {LL,L,M,R}
		//! -> BlockSize*2 / (MAX_DECIMATION*4)
		//! =  BlockSize / (MAX_DECIMATION*2)
		int i, BinSize = BlockSize/(ULC_MAX_BLOCK_DECIMATION_FACTOR*2);
		float v, Sum;
		float LPTap = SmoothingTaps[0], LPDecay = 252/256.0f, OneMinusLPDecay = 1.0f - LPDecay;
		float DCTap = SmoothingTaps[1], DCDecay = 255/256.0f, OneMinusDCDecay = 1.0f - DCDecay;
		float DCGainCompensation = OneMinusDCDecay / OneMinusLPDecay;
		      float *Dst = TransientBuffer + ULC_MAX_BLOCK_DECIMATION_FACTOR*2; //! Align to M,R segment
		const float *Src = TmpBuffer;
		{
			//! Update the filter taps, but do NOT fix up the
			//! values in the previous {M,R} segment, or this
			//! will cause issues with transients from silence.
			DOFILTER(); //Dst[+ULC_MAX_BLOCK_DECIMATION_FACTOR*2-1] += v;
		}
		n = ULC_MAX_BLOCK_DECIMATION_FACTOR*2; do {
			//! Filter and accumulate for this bin
			Sum = 0.0f;
			i = BinSize - (n == 1); do { //! Last sample unavailable, so skip it in the last bin
				DOFILTER(); Sum += v;
			} while(--i);

			//! {LL,L} = {M,R}, then set {M,R} with new data
			Dst[-ULC_MAX_BLOCK_DECIMATION_FACTOR*2] = *Dst;
			*Dst++ = Sum;
		} while(--n);
		SmoothingTaps[0] = LPTap;
		SmoothingTaps[1] = DCTap;
#undef DOFILTER
	}
}
#pragma GCC pop_options
static inline float Block_Transform_GetWindowCtrl_Log2DecimationRatio(float Ratio2, int Log2SubBlockSize) {
	//! Full, unsimplified expression:
	//!  LogRatio2         = Log[Ratio^2] = 2*Log[Ratio]
	//!  OverlapSamples    = E^(-2*Log[Ratio]) * 10000; experimentally determined
	//!  OverlapDecimation = Log2[SubBlockSize / OverlapSamples] = Log2[SubBlockSize] - Log2[10000] + Log2[Ratio^2]
	return Log2SubBlockSize - 0x1.A934F1p3f + 0x1.715476p0f*logf(Ratio2);
}
static inline int Block_Transform_GetWindowCtrl(
	const float *ThisBlockData,
	const float *LastBlockData,
	      float *TransientBuffer,
	      float *TmpBuffer,
	      float *SmoothingTaps,
	      int    BlockSize,
	      int    nChan
) {
	int n;

	//! Perform filtering to obtain pre-echo analysis
	Block_Transform_GetWindowCtrl_TransientFiltering(ThisBlockData, LastBlockData, TransientBuffer, TmpBuffer, SmoothingTaps, BlockSize, nChan);

	//! Split the block until the transient signal stabilizes,
	//! then take the ratio on the L side for overlap scaling.
	float DecimationRatio;
	int Decimation  = 0b0001;
	int AnalysisLen = ULC_MAX_BLOCK_DECIMATION_FACTOR;
	int Log2SubBlockSize = 31 - __builtin_clz(BlockSize);
	for(TransientBuffer += AnalysisLen*3;;) { //! Pivot from the next block's center (R segment)
		//! Find the peak ratio for the L/M segments
		enum { POS_L, POS_M, POS_N};
		int PeakPos;
		float Ratio[POS_N];  {
			//! Get the energy of each segment (LL/L/M)
			//! NOTE: Do not use FLT_MIN as the bias, as we need
			//! some room for the ratio to grow into upon division.
			//! NOTE: The R segment doesn't matter, as it becomes the
			//! L segment in the next block, so we will analyze then.
			float LL = 0x1.0p-64f;
			float L  = 0x1.0p-64f;
			float M  = 0x1.0p-64f;
			for(n=0;n<AnalysisLen;n++) {
				LL += TransientBuffer[-3*AnalysisLen + n];
				L  += TransientBuffer[-2*AnalysisLen + n];
				M  += TransientBuffer[-1*AnalysisLen + n];
			}

			//! Get the ratios between the segments and select the largest
			   (Ratio[POS_L] = L / LL);                 PeakPos = POS_L;
			if((Ratio[POS_M] = M / L) > Ratio[PeakPos]) PeakPos = POS_M;
		}

		//! If the transient merits decimation, keep going.
		//! NOTE: Minimum subblock size of 64 samples.
		//! NOTE: Checking AnalysisLen should be better than checking
		//! Decimation directly, as then we can change the maximum allowed
		//! decimation without changing this code.
		if(ULC_USE_WINDOW_SWITCHING && AnalysisLen > 1 && Log2SubBlockSize > 6) {
			float r = Block_Transform_GetWindowCtrl_Log2DecimationRatio(Ratio[PeakPos], Log2SubBlockSize);
			if(r > 0.0f) {
				//! Update the decimation pattern and continue
				//! NOTE: When PeakPos==L, we've simply shifted the
				//! next subblock's center point, hence pivoting about
				//! the R segment instead of anywhere else.
				Decimation = (Decimation<<1) | (PeakPos != POS_L);
				TransientBuffer -= AnalysisLen * (PeakPos == POS_L);
				AnalysisLen /= 2;
				Log2SubBlockSize--;
				continue;
			}
		}

		//! No more decimation - Store the L ratio and break out
		DecimationRatio = Block_Transform_GetWindowCtrl_Log2DecimationRatio(Ratio[POS_L], Log2SubBlockSize);
		break;
	}

	//! Determine overlap size from the ratio
	int OverlapScale = (DecimationRatio < 0.5f) ? 0 : (DecimationRatio >= 6.5f) ? 7 : (int)(0.5f + DecimationRatio);
	if(Log2SubBlockSize-OverlapScale < 4) OverlapScale = Log2SubBlockSize-4; //! Minimum 16-sample overlap

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
