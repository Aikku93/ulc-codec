/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2023, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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
//!   and will be updated as (L = R_Old, R = New). Initialize with 0.
//!   Internal layout is:
//!   {
//!     struct ULC_TransientData_t L[ULC_MAX_BLOCK_DECIMATION_FACTOR];
//!     struct ULC_TransientData_t R[ULC_MAX_BLOCK_DECIMATION_FACTOR];
//!   }
//!  -TransientFilter[] must be 3 elements in size. Initialize with 0.
//!   Internal layout is:
//!   {
//!     float EnvGainHP; //! Gain envelope tap for HP filter
//!     float EnvGainBP; //! Gain envelope tap for BP filter
//!     float EnvPost;   //! Post-echo-compensated gain tap
//!   }
//!  -TmpBuffer[] must be BlockSize*2 elements in size.
#pragma GCC push_options
#pragma GCC optimize("fast-math") //! Should improve things, hopefully, maybe
static inline void Block_Transform_GetWindowCtrl_TransientFiltering(
	const float *BlockData,
	struct ULC_TransientData_t *TransientBuffer,
	      float *TransientFilter,
	      float *TmpBuffer,
	      int    BlockSize,
	      int    nChan,
	      int    RateHz
) {
	int n, Chan;

	//! Extract the energy of a highpass and bandpass filter
	//! Transfer functions:
	//!  H(z) = -z^-1 + 2 - z^1 (Highpass; Gain = 4.0)
	//!  H(z) = z^1 - z^-1      (Bandpass; Gain = 2.0)
	//! This section outputs interleaved {Highpass,Bandpass}
	//! into BufEnergy[].
	float *BufEnergy = TmpBuffer; {
		for(n=0;n<BlockSize*2;n++) BufEnergy[n] = 0.0f;
		for(Chan=0;Chan<nChan;Chan++) {
#define DOFILTER(XzM1,Xz0,Xz1) *Dst++ += SQR(-XzM1+2*Xz0-Xz1), *Dst++ += SQR(Xz1-XzM1)
			      int    Lag    = BlockSize/2; //! MDCT alignment
			      float *Dst    = BufEnergy;
			const float *SrcOld = BlockData + Chan*BlockSize + BlockSize-Lag;
			const float *SrcNew = BlockData + Chan*BlockSize + nChan*BlockSize;
			n = Lag-1; do {
				DOFILTER(SrcOld[-1], SrcOld[0], SrcOld[+1]), SrcOld++;
			} while(--n);
			{
				DOFILTER(SrcOld[-1], SrcOld[0], SrcNew[ 0]), SrcOld++;
				DOFILTER(SrcOld[-1], SrcNew[0], SrcNew[+1]), SrcNew++;
			}
			n = BlockSize - (Lag-1) - 2; do {
				DOFILTER(SrcNew[-1], SrcNew[0], SrcNew[+1]), SrcNew++;
			} while(--n);
#undef DOFILTER
		}
	}

	//! Model the energy curve and integrate it over each segment
	{
		int i, BinSize = BlockSize / ULC_MAX_BLOCK_DECIMATION_FACTOR;
		float EnvGainHP = TransientFilter[0], GainHPRate = expf(-0x1.5A92D6p2f / RateHz); //! -0.047dB/ms (1000 * Log[2^-(1/128)])
		float EnvGainBP = TransientFilter[1], GainBPRate = expf(-0x1.5A92D6p2f / RateHz); //! -0.047dB/ms (1000 * Log[2^-(1/128)])
		float EnvPost   = TransientFilter[2], PostRate   = expf(-0x1.5A92D6p6f / RateHz); //! -0.753dB/ms (1000 * Log[2^-(1/8)])
		struct ULC_TransientData_t *Dst = TransientBuffer + ULC_MAX_BLOCK_DECIMATION_FACTOR; //! Align to new block
		const float *Src = BufEnergy;
		i = ULC_MAX_BLOCK_DECIMATION_FACTOR; do {
			//! Swap out the old "new" data, and clear the new "new" data
			Dst[-ULC_MAX_BLOCK_DECIMATION_FACTOR] = *Dst;
			*Dst = (struct ULC_TransientData_t){.Sum = 0.0f, .SumW = 0.0f};
			n = BinSize; do {
				//! Accumulate and subtract 'drift' from the HP and BP
				//! signals to "unbias" the data for the next steps.
				//! This essentially 'enhances' energy differences.
				//! NOTE: This calculation must be done in the amplitude
				//! domain, as the power domain behaves too erratically.
				float vHP = sqrtf(Src[0]) - EnvGainHP;
				float vBP = sqrtf(Src[1]) - EnvGainBP;
				EnvGainHP += vHP * (1.0f-GainHPRate);
				EnvGainBP += vBP * (1.0f-GainBPRate);
				Src += 2;

				//! Update the post-echo-compensating energy curve, and
				//! store the updated sum for this segment.
				//! NOTE: Delta HP and delta BP are cross-multiplied by
				//! their respective gains to sort of 'normalize' with
				//! respect to one another. This balances out their
				//! weaknesses somewhat.
				float vPost = SQR(vHP*EnvGainBP) + SQR(vBP*EnvGainHP) - EnvPost;
				EnvPost += vPost * (1.0f-PostRate);
				Dst->Sum += EnvPost, Dst->SumW += (SQR(EnvGainBP) + SQR(EnvGainHP));
			} while(--n);
		} while(Dst++, --i);
		TransientFilter[0] = EnvGainHP;
		TransientFilter[1] = EnvGainBP;
		TransientFilter[2] = EnvPost;
	}
}
#pragma GCC pop_options
static inline int Block_Transform_GetWindowCtrl(
	const float *BlockData,
	struct ULC_TransientData_t *TransientBuffer,
	      float *TransientFilter,
	      float *TmpBuffer,
	      int    BlockSize,
	      int    nChan,
	      int    RateHz
) {
	int n;

	//! Perform filtering to obtain transient analysis
	//! then seek to this "new" block's transient data
	Block_Transform_GetWindowCtrl_TransientFiltering(BlockData, TransientBuffer, TransientFilter, TmpBuffer, BlockSize, nChan, RateHz);
	TransientBuffer += ULC_MAX_BLOCK_DECIMATION_FACTOR;

	//! Keep trying to increase the window size until the
	//! transient drops off compared to the last window
	//! NOTE: Default to maximum overlap, no decimation
	int   Log2SubBlockSize = 31 - __builtin_clz(BlockSize/ULC_MAX_BLOCK_DECIMATION_FACTOR);
	int   Decimation     = 0b0001;
	float TransientRatio = 0.0f; {
		//! First, enforce a minimum SubBlockSize of 64
		int nSegments   = ULC_MAX_BLOCK_DECIMATION_FACTOR;
		int SegmentSize = 1;
		if(Log2SubBlockSize < 6) {
			int Shift = 6 - Log2SubBlockSize;
			nSegments   >>= Shift;
			SegmentSize <<= Shift;
			Log2SubBlockSize = 6;
		}

		//! Begin searching for the right windows
		for(;;) {
			//! For this iteration, we are preparing for a larger window.
			//! Putting this here avoids awkward formulations later on.
			Log2SubBlockSize++;

			//! Get the maximum Attack/Release ratio for all segments
			//! by removing the masking (Release) level and then comparing
			//! the remaining energy to that of the last segment
			int Segment;
			int MaxSegment = 0;
			float AvgRatio = 0.0f, MaxRatio = -1000.0f;
			for(Segment=0;Segment<nSegments;Segment++) {
				struct ULC_TransientData_t L = {.Sum = 0.0f, .SumW = 0.0f};
				struct ULC_TransientData_t R = {.Sum = 0.0f, .SumW = 0.0f};
				const struct ULC_TransientData_t *Src = TransientBuffer + Segment*SegmentSize;
				for(n=0;n<SegmentSize;n++) {
#define ADDSEGMENT(Dst, Src) Dst.Sum += Src.Sum, Dst.SumW += Src.SumW
					ADDSEGMENT(L, Src[n-SegmentSize]);
					ADDSEGMENT(R, Src[n]);
#undef ADDSEGMENT
				}
				L.Sum = L.Sum ? logf(L.Sum / L.SumW) : (-100.0f); //! -100 = Placeholder for Log[0]
				R.Sum = R.Sum ? logf(R.Sum / R.SumW) : (-100.0f);

				//! Get final energy ratio
				float Ratio = ABS(R.Sum - L.Sum);
				AvgRatio += Ratio;
				if(Ratio > MaxRatio) MaxSegment = Segment, MaxRatio = Ratio;
			}
			AvgRatio /= nSegments;

			//! If this window causes the transient ratio to drop too much,
			//! stop enlarging and use the windows from the last iteration.
			//! On the other hand: If the ratio /increases/ as the window
			//! size increases, then the transient is likely close to the
			//! center of two subblocks; in this case, we /can/ continue
			//! increasing the window size, but we must be very careful to
			//! only do this if the transient is significant enough - when
			//! it's not, increasing the window size will cause 'metallic'
			//! artifacts around the transient, as well as post-echo.
			if(MaxRatio-TransientRatio < 0x1.62E430p-1f) break; //! 0x1.62E430p-1 = Log[2]

			//! Set new decimation pattern + transient ratio, and continue
			//! NOTE: Only try larger subblocks if the transient isn't too powerful,
			//! otherwise we risk smearing it, which would defeat the purpose.
			Decimation     = nSegments + MaxSegment;
			TransientRatio = MaxRatio;
			if(nSegments > 1 && TransientRatio < 0x1.62E430p-1f) { //! 0x1.62E430p-1 = Log[2]
				//! Increase spectral resolution
				nSegments   /= 2;
				SegmentSize *= 2;
			} else break;
		}
	}

	//! If the transient isn't powerful enough, just return a full-overlap block
	if(TransientRatio < 0x1.62E430p-2f) return 0x10;

	//! Determine overlap size from the ratio
	//! NOTE: Log2SubBlockSize was pre-incremented earlier, so account for this here.
	TransientRatio *= 0x1.715476p0f; //! 0x1.715476p0 = 1/Log[2] for change of base
	int OverlapScale = (TransientRatio < 0.5f) ? 0 : (TransientRatio >= 6.5f) ? 7 : lrintf(TransientRatio);
	if(Log2SubBlockSize-OverlapScale < 5+1) OverlapScale = Log2SubBlockSize - (5+1); //! Minimum 32-sample overlap

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
