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
//!  -TransientFilter[] must be 2 elements in size
#pragma GCC push_options
#pragma GCC optimize("fast-math") //! Should improve things, hopefully, maybe
static inline void Block_Transform_GetWindowCtrl_TransientFiltering(
	const float *ThisBlockData,
	const float *LastBlockData,
	      float *TransientBuffer,
	      float *TransientFilter,
	      float *TmpBuffer,
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
	//! This is compensated for later.
	//! NOTE: We introduce a 1-sample delay here to compensate
	//! for incomplete data at n=BlockSize-1.
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

	//! Pass the signal power through two smoothing filters,
	//! and then take their difference as the transient signal.
	//! One filter operates in the Power domain, and the other
	//! in a companded domain, which appears to capture the
	//! intent of transients fairly cleanly.
	{
		//! BlockSize(={LL,L}) + BlockSize(={M,R}) = BlockSize*2 = {LL,L,M,R}
		//! -> BlockSize*2 / (MAX_DECIMATION*4)
		//! =  BlockSize / (MAX_DECIMATION*2)
		//! NOTE: The values are renormalized to avoid issues
		//! with subnormal collapse. The output of the filter
		//! steps is generally a low-amplitude signal, so it
		//! stands to reason that it can collapse as such.
		//! It should be noted that the gain will be squared
		//! during the summation (eg. 2^32 becomes 2^64) for
		//! the weighted sum.
		int i, BinSize = BlockSize / (ULC_MAX_BLOCK_DECIMATION_FACTOR*2);
		float Gain = 0x1.0p48f / (4 * nChan); //! The filter from earlier has a gain of 2.0, and is then squared, giving 4.0
		float PowerTap = TransientFilter[0], PowerDecay = 255/256.0f;
		float FloorTap = TransientFilter[1], FloorDecay = 252/256.0f;
		      float *Dst = TransientBuffer + ULC_MAX_BLOCK_DECIMATION_FACTOR*2; //! Align to M,R segment
		const float *Src = TmpBuffer;
		i = ULC_MAX_BLOCK_DECIMATION_FACTOR*2; do {
			float v, Sum = 0.0f, SumW = 0.0f;
			n = BinSize; do {
				//! Extract the transient energy signal
				v         = *Src++ * Gain;
				FloorTap += (v - FloorTap)*(1.0f-FloorDecay);
				v         = sqrtf(sqrtf(v));
				PowerTap += (v - PowerTap)*(1.0f-PowerDecay);
				v        = ABS(SQR(SQR(PowerTap)) - FloorTap);
				Sum      += v*v, SumW += v;
			} while(--n);

			//! {LL,L} = {M,R}, then set {M,R} with new data
			Dst[-ULC_MAX_BLOCK_DECIMATION_FACTOR*2] = *Dst;
			*Dst++ = Sum ? (Sum/SumW) : 0.0f;
		} while(--i);
		TransientFilter[0] = PowerTap;
		TransientFilter[1] = FloorTap;
	}
}
#pragma GCC pop_options
static inline float Block_Transform_GetWindowCtrl_Log2DecimationRatio(float LogRatio, int Log2SubBlockSize) {
	//! I have no idea what is going on here; this is from trial and error.
	//! 0x1.715476p0 = 1/Log[2] for change of base
	return (Log2SubBlockSize - 12) + 0x1.715476p0f*LogRatio;
}
static inline int Block_Transform_GetWindowCtrl(
	const float *ThisBlockData,
	const float *LastBlockData,
	      float *TransientBuffer,
	      float *TransientFilter,
	      float *TmpBuffer,
	      int    BlockSize,
	      int    nChan
) {
	int n;

	//! Perform filtering to obtain transient analysis
	Block_Transform_GetWindowCtrl_TransientFiltering(ThisBlockData, LastBlockData, TransientBuffer, TransientFilter, TmpBuffer, BlockSize, nChan);

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
			//! NOTE: The R segment doesn't matter, as it becomes the
			//! L segment in the next block, so we will analyze then.
			//! NOTE: Use logarithms to avoid overflow on division.
			float v;
			float LL = 0.0f, LLw = 0.0f;
			float L  = 0.0f, Lw  = 0.0f;
			float M  = 0.0f, Mw  = 0.0f;
			const float *Src = TransientBuffer;
			const float MIN_LOG = -100.0f;
			n = AnalysisLen; do v = *--Src, M  += v*v, Mw  += v; while(--n);
			n = AnalysisLen; do v = *--Src, L  += v*v, Lw  += v; while(--n);
			n = AnalysisLen; do v = *--Src, LL += v*v, LLw += v; while(--n);
			LL = LL ? logf(LL / LLw) : MIN_LOG;
			L  = L  ? logf(L  / Lw)  : MIN_LOG;
			M  = M  ? logf(M  / Mw)  : MIN_LOG;

			//! Get the ratios between the segments and select the largest
			Ratio[POS_L] = L - LL;
			Ratio[POS_M] = M - L;
			PeakPos = (Ratio[POS_L] > Ratio[POS_M]) ? POS_L : POS_M;
		}

		//! If there is a transient in the M segment (M > L)
		//! or there is post-echo (L > M), keep splitting.
		//! NOTE: Minimum subblock size of 64 samples.
		//! NOTE: Checking AnalysisLen should be better than checking
		//! Decimation directly, as then we can change the maximum allowed
		//! decimation without changing this code.
		if(ULC_USE_WINDOW_SWITCHING && AnalysisLen > 1 && Log2SubBlockSize > 6) {
			if(ABS(Ratio[POS_M] - Ratio[POS_L]) > 0x1.62E430p0f) { //! Log[2]
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
		DecimationRatio = Block_Transform_GetWindowCtrl_Log2DecimationRatio(ABS(Ratio[POS_L]), Log2SubBlockSize);
		break;
	}

	//! Determine overlap size from the ratio
	int OverlapScale = (DecimationRatio <= 0.0f) ? 0 : (DecimationRatio >= 7.0f) ? 7 : (int)DecimationRatio;
	if(Log2SubBlockSize-OverlapScale < 4) OverlapScale = Log2SubBlockSize-4; //! Minimum 16-sample overlap

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
