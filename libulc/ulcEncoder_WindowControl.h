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
//!  -TransientBuffer[] must be ULC_MAX_BLOCK_DECIMATION_FACTOR*2 in size
//!   and will be updated as (L = R_Old, R = New). This array will store
//!   logarithmic values, so must be initialized as an approximation to
//!   Log[0] (-100 should be fine for all cases).
//!  -TmpBuffer[] must be BlockSize*2 in size
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

	//! Perform a highpass and bandpass filter over the signal.
	//! This removes unwanted low frequency oscillations (eg.
	//! beating) that can cause issues on transient detection,
	//! and gives us separate spectra to analyze, as in some
	//! signals, transients only have very short-lived HF data
	//! (with correspondingly longer-lived MF/LF data), while
	//! more common transients simply have a lot of HF data.
	//! Transfer functions:
	//!  Highpass: H(z) = -z^-1 + 2 - z^1 (Gain: 4.0)
	//!  Bandpass: H(z) = -z^-1     + z^1 (Gain: 2.0)
	//! NOTE: We need to center the buffer between the two
	//! blocks (Old/New) for correct MDCT window alignment.
	//! With window switching enabled, we also want to fix
	//! edge cases where we need to apply jitter correction,
	//! which means the "true" center is off by one segment.
	//! Without window switching, we present the data fully
	//! centered, but cannot perform jitter correction. It
	//! really should not be necessary, though, as it won't
	//! improve the output much, if at all.
	//! NOTE: The filters used here is noncausal, so we
	//! apply a 1-sample delay to compensate.
	for(n=0;n<BlockSize*2;n++) TmpBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
#define DOHP(zM1, z0, z1) *Dst++ += SQR(-(zM1) + 2*(z0) - (z1))
#define DOBP(zM1, z0, z1) *Dst++ += SQR(-(zM1)          + (z1))
#define DOFILTER(zM1, z0, z1) DOHP(zM1, z0, z1), DOBP(zM1, z0, z1)
		int Lag = BlockSize/2; if(ULC_USE_WINDOW_SWITCHING) Lag -= BlockSize/ULC_MAX_BLOCK_DECIMATION_FACTOR;
		Lag += 1; //! Noncausal filter compensation
		float *Dst = TmpBuffer;
		const float *SrcOld = LastBlockData + Chan*BlockSize + BlockSize-Lag;
		const float *SrcNew = ThisBlockData + Chan*BlockSize;
		for(n=0;n<Lag-1;n++) {
			DOFILTER(SrcOld[-1], SrcOld[ 0], SrcOld[+1]), SrcOld++;
		}
		{
			DOFILTER(SrcOld[-1], SrcOld[ 0], SrcNew[ 0]), SrcOld++, n++;
			DOFILTER(SrcOld[-1], SrcNew[ 0], SrcNew[+1]), SrcNew++, n++;
		}
		for(;n<BlockSize;n++) {
			DOFILTER(SrcNew[-1], SrcNew[ 0], SrcNew[+1]), SrcNew++;
		}
#undef DOFILTER
#undef DOBP
#undef DOHP
	}

	//! Pass the signal through a compressor envelope
	//! generator, and then divide the signal by this
	//! envelope (in effect producing an expander),
	//! which then emphasizes transients.
	{
		//! NOTE: (UNITY_GAIN/SILENCE_BIAS)^2 must not overflow.
		const float UNITY_GAIN   = 0x1.0p31f;
		const float SILENCE_BIAS = 0x1.0p-32f;
		int i, BinSize = BlockSize / ULC_MAX_BLOCK_DECIMATION_FACTOR;
		float HPTap = TransientFilter[0], HPDecay = 254/256.0f, HPGain = (UNITY_GAIN / SQR(4.0f)) / nChan;
		float BPTap = TransientFilter[1], BPDecay = 254/256.0f, BPGain = (UNITY_GAIN / SQR(2.0f)) / nChan;
		      float *Dst = TransientBuffer + ULC_MAX_BLOCK_DECIMATION_FACTOR; //! Align to new block
		const float *Src = TmpBuffer;
		i = ULC_MAX_BLOCK_DECIMATION_FACTOR; do {
			float HP, SumHP = 0.0f;
			float BP, SumBP = 0.0f;
			n = BinSize; do {
				HP     = *Src++ * HPGain;
				BP     = *Src++ * BPGain;
				HPTap += (UNITY_GAIN - HP)*(1.0f-HPDecay);
				BPTap += (UNITY_GAIN - BP)*(1.0f-BPDecay);
				SumHP += HP / (HPTap + SILENCE_BIAS);
				SumBP += BP / (BPTap + SILENCE_BIAS);
			} while(--n);

			//! Swap out the old "new" data, and replace with new "new" data
			float Sum = SumHP*SumBP;
			Dst[-ULC_MAX_BLOCK_DECIMATION_FACTOR] = *Dst;
			*Dst++ = Sum ? logf(Sum) : (-100.0f); //! -100 = Placeholder for Log[0]
		} while(--i);
		TransientFilter[0] = HPTap;
		TransientFilter[1] = BPTap;
	}
}
#pragma GCC pop_options
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
	int Log2BlockSize = 31 - __builtin_clz(BlockSize);
	float L, R, r;

	//! Control thresholds
	//! NOTE: The "raw" thresholds (without THRES_CORRECTION) were
	//! tuned for BlockSize=2048, so we compensate by multiplying
	//! with 2048/BlockSize:
	//!  Log[CorrectedThreshold] = Log[Threshold * 2048/BlockSize]
	//!  =Log[Threshold] + Log[2048/BlockSize]
	//!  =Log[Threshold] + (1/Log2[E])*Log2[2048/BlockSize]
	//!  =Log[Threshold] + (1/Log2[E])*(11 - Log2[BlockSize])
	const float THRES_CORRECTION = 0x1.62E430p-1f*(11 - Log2BlockSize); //! 0x1.62E430p-1 = 1/Log2[E], for change of base
	const float ATT_THRES  = THRES_CORRECTION + 0x1.62E430p0f;  //! Log[ 4]; attack threshold
	const float DEC_THRES  = THRES_CORRECTION + 0x1.62E430p1f;  //! Log[16]; decay threshold
	const float LEAK_THRES = THRES_CORRECTION + 0x1.62E430p-1f; //! Log[ 2]; post-echo leakage threshold
	const float EDGE_THRES = THRES_CORRECTION + 0x1.62E430p-1f; //! Log[ 2]; edge-case threshold (x[n+1]/x[n-1] must be this much higher/lower than x[n]/x[n-1])

	//! Perform filtering to obtain transient analysis
	//! then seek to this "new" block's transient data
	//! (which is centered between the two blocks).
	//! NOTE: Without window switching, the output of
	//! the filter has already been centered and we
	//! simply move to the next block at once.
	//! With window switching, we're shifted by one
	//! segment relative to the center, so we adjust.
	Block_Transform_GetWindowCtrl_TransientFiltering(ThisBlockData, LastBlockData, TransientBuffer, TransientFilter, TmpBuffer, BlockSize, nChan);
	TransientBuffer += ULC_MAX_BLOCK_DECIMATION_FACTOR - ULC_USE_WINDOW_SWITCHING;

	//! Find the most relevant transient
	int   TransientIndex;
	float TransientRatio; {
		//! Find the segment with the largest attack/release
		int MaxIndex = 0; float MaxRatio = -1000.0f;
		int MinIndex = 0; float MinRatio = +1000.0f;
		for(n=0;n<ULC_MAX_BLOCK_DECIMATION_FACTOR;n++) {
			L = TransientBuffer[n-1];
			R = TransientBuffer[n];
			r = R-L;
			if(r > MaxRatio) MaxIndex = n, MaxRatio = r;
			if(r < MinRatio) MinIndex = n, MinRatio = r;
		}

		//! In some cases, jitter happens and our transient ratio
		//! becomes skewed due to a higher-than-expected value in
		//! the previous segment. This attempts to compensate by
		//! comparing the next segment against the previous, and
		//! deciding whether to use that ratio instead of
		//! the 'canonical' ratio of x[n]/x[n-1].
		//! NOTE: Despite what intuition tells us, we dot NOT
		//! alter the transient subblock index. This can cause
		//! serious artifacts when transients "only just" begin
		//! on that subblock, as then that large spike in energy
		//! will pre-echo very audibly. A similar issue happens
		//! with decay transients and post-echo.
		if(ULC_USE_WINDOW_SWITCHING) {
			//! Attack
			L = TransientBuffer[MaxIndex-1];
			R = TransientBuffer[MaxIndex+1];
			r = R-L;
			if(r > MaxRatio+EDGE_THRES) MaxRatio = r;

			//! Decay
			L = TransientBuffer[MinIndex-2];
			R = TransientBuffer[MinIndex];
			r = R-L;
			if(r < MinRatio-EDGE_THRES) MinRatio = r;
		}

		//! If we don't have a significant attack, use release as transient marker
		     if(MaxRatio > +ATT_THRES) TransientIndex = MaxIndex, TransientRatio =  MaxRatio - ATT_THRES;
		else if(MinRatio < -DEC_THRES) TransientIndex = MinIndex, TransientRatio = -MinRatio - DEC_THRES;
		else                           TransientIndex = -1,       TransientRatio =  0.0f;
	}

	//! Attempt to enlarge the window to the right (if there is no post-echo)
	int Decimation       = 0b0001; //! Default to SubBlocks={N/1}
	int Log2SubBlockSize = Log2BlockSize;
	if(TransientIndex != -1) {
		int AnalysisLen = 1;
		int PostEchoIsHandled = 0;
		Decimation        = ULC_MAX_BLOCK_DECIMATION_FACTOR + TransientIndex;
		Log2SubBlockSize -= 31 - __builtin_clz(ULC_MAX_BLOCK_DECIMATION_FACTOR); //! Log2[BlockSize/MAX_DECIMATION_FACTOR]
		while((Decimation & 1) == 0) {
			//! Check for decaying segments
			//! NOTE: We already checked AnalysisLen/2 segments,
			//! so we can start at AnalysisLen/2.
			for(n=AnalysisLen/2;n<AnalysisLen;n++) {
				L = TransientBuffer[TransientIndex+n];
				R = TransientBuffer[TransientIndex+n+1];
				r = R-L;
				if(r < -LEAK_THRES) { PostEchoIsHandled = 1; break; }
			}
			if(PostEchoIsHandled) break;

			//! No decay transient (yet) - Keep enlarging
			Log2SubBlockSize++;
			AnalysisLen *= 2;
			Decimation >>= 1;
		}

		//! When leakage happens on the edge of a block,
		//! ensure that the next segment triggers a decay
		//! transient event.
		if(!PostEchoIsHandled) {
			L = TransientBuffer[ULC_MAX_BLOCK_DECIMATION_FACTOR-1];
			R = TransientBuffer[ULC_MAX_BLOCK_DECIMATION_FACTOR];
			r = R-L;
			if(r < -LEAK_THRES) {
				TransientBuffer[ULC_MAX_BLOCK_DECIMATION_FACTOR-1] = L + (DEC_THRES - LEAK_THRES);
			}
		}
	}

	//! Determine overlap size from the ratio
	float DecimationRatio = (Log2SubBlockSize - 12) + 0x1.715476p1f*TransientRatio; //! 0x1.715476p0 = 1/Log[2] for change of base
	int OverlapScale = (DecimationRatio <= 0.0f) ? 0 : (DecimationRatio >= 7.0f) ? 7 : (int)DecimationRatio;
	if(Log2SubBlockSize-OverlapScale < 4) OverlapScale = Log2SubBlockSize-4; //! Minimum 16-sample overlap

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
