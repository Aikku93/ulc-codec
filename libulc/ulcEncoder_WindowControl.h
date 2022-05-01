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
//!   and will be updated as (L = R_Old, R = New). Initialize with 0.
//!   Internal layout is:
//!   {
//!     struct ULC_TransientData_t L[ULC_MAX_BLOCK_DECIMATION_FACTOR];
//!     struct ULC_TransientData_t R[ULC_MAX_BLOCK_DECIMATION_FACTOR];
//!   }
//!  -TransientFilter[] must be 3 elements in size. Initialize with 0.
//!   Internal layout is:
//!   {
//!     float EnvGain; //! Gain envelope tap
//!     float EnvAtt;  //! Attack envelope tap
//!     float EnvRel;  //! Release envelope tap
//!   }
//!  -TmpBuffer[] must be BlockSize*2 elements in size.
#pragma GCC push_options
#pragma GCC optimize("fast-math") //! Should improve things, hopefully, maybe
static inline void Block_Transform_GetWindowCtrl_TransientFiltering(
	const float *ThisBlockData,
	const float *LastBlockData,
	struct ULC_TransientData_t *TransientBuffer,
	      float *TransientFilter,
	      float *TmpBuffer,
	      int    BlockSize,
	      int    nChan,
	      int    RateHz
) {
	int n, Chan;

	//! Extract the energy of a highpass and bandpass filter,
	//! which will later be combined as a geometric mean.
	//! The filters are 90deg out of phase w.r.t. one another
	//! such that multiplying them together should not result
	//! in beating artifacts, allowing cleaner extraction of
	//! the transient signal ("mid" band modulated by the
	//! higher-energy/high-frequency band).
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
			const float *SrcOld = LastBlockData + Chan*BlockSize + BlockSize-Lag;
			const float *SrcNew = ThisBlockData + Chan*BlockSize;
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

	//! Model the attack and release envelopes, then
	//! integrate them separately over each segment.
	{
		//! The AttRate parameter refers to how long a signal takes to reach unity.
		//! The RelRate parameter refers to how long a signal takes to decay to 0.
		int i, BinSize = BlockSize / ULC_MAX_BLOCK_DECIMATION_FACTOR;
		float EnvGain = TransientFilter[0], GainRate = expf(-0x1.0A2B24p3f / BlockSize); //! -36dB/block (Log[2^-6])
		float EnvAtt  = TransientFilter[1], AttRate  = expf(-0x1.5A92D7p6f / RateHz); //! -0.75dB/ms (1000 * Log[2^-0.125])
		float EnvRel  = TransientFilter[2], RelRate  = expf(-0x1.5A92D7p5f / RateHz); //! -0.38dB/ms (1000 * Log[2^-0.0625])
		struct ULC_TransientData_t *Dst = TransientBuffer + ULC_MAX_BLOCK_DECIMATION_FACTOR; //! Align to new block
		const float *Src = BufEnergy;
		i = ULC_MAX_BLOCK_DECIMATION_FACTOR; do {
			//! Swap out the old "new" data, and clear the new "new" data
			Dst[-ULC_MAX_BLOCK_DECIMATION_FACTOR] = *Dst;
			*Dst = (struct ULC_TransientData_t){.Att = 0.0f, .Rel = 0.0f};
			n = BinSize; do {
				//! Accumulate and subtract the "drift" bias
				//! (ie. the amplitude at which everything is
				//! happening around) to "unbias" the signal
				//! and center it about 0 for attack/release.
				//! NOTE: Square root to extract geometric mean
				//! of combined highpass/bandpass energy.
				float v  = sqrtf(Src[0]*Src[1]) - EnvGain; Src += 2;
				EnvGain += v*(1.0f-GainRate);

				//! Update attack/release envelopes and integrate
				//! NOTE: EnvRel is sign-inverted.
				//! NOTE: Envelope attack (fade-in) is proportional to
				//! the signal level, whereas envelope decay (fade-out)
				//! follows a constant falloff.
				if(v >= EnvAtt) EnvAtt += (v-EnvAtt)*(1.0f-AttRate); else EnvAtt *= AttRate;
				if(v <= EnvRel) EnvRel += (v-EnvRel)*(1.0f-RelRate); else EnvRel *= RelRate;
				Dst->Att += EnvAtt;
				Dst->Rel -= EnvRel; //! EnvRel is sign-inverted, so flip again
			} while(--n);
		} while(Dst++, --i);
		TransientFilter[0] = EnvGain;
		TransientFilter[1] = EnvAtt;
		TransientFilter[2] = EnvRel;
	}
}
#pragma GCC pop_options
static inline int Block_Transform_GetWindowCtrl(
	const float *ThisBlockData,
	const float *LastBlockData,
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
	Block_Transform_GetWindowCtrl_TransientFiltering(ThisBlockData, LastBlockData, TransientBuffer, TransientFilter, TmpBuffer, BlockSize, nChan, RateHz);
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
			int MaxSegment = 0; float MaxRatio = -1000.0f;
			for(Segment=0;Segment<nSegments;Segment++) {
				struct ULC_TransientData_t L = {.Att = 0.0f, .Rel = 0.0f};
				struct ULC_TransientData_t R = {.Att = 0.0f, .Rel = 0.0f};
				const struct ULC_TransientData_t *Src = TransientBuffer + Segment*SegmentSize;
				for(n=0;n<SegmentSize;n++) {
#define ADDSEGMENT(Dst, Src) Dst.Att += Src.Att, Dst.Rel += Src.Rel
					ADDSEGMENT(L, Src[n-SegmentSize]);
					ADDSEGMENT(R, Src[n]);
#undef ADDSEGMENT
				}
				L.Att = L.Att ? logf(L.Att) : (-100.0f); //! -100 = Placeholder for Log[0]
				L.Rel = L.Rel ? logf(L.Rel) : (-100.0f);
				R.Att = R.Att ? logf(R.Att) : (-100.0f);
				R.Rel = R.Rel ? logf(R.Rel) : (-100.0f);
				float Ratio = ABS((R.Att - R.Rel) - (L.Att - L.Rel));
				if(Ratio > MaxRatio) MaxSegment = Segment, MaxRatio = Ratio;
			}

			//! If this window causes the transient ratio to drop too much,
			//! stop enlarging and use the windows from the last iteration.
			//! On the other hand: If the ratio /increases/ as the window
			//! size increases, then the transient is likely close to the
			//! center of two subblocks; in this case, we /can/ continue
			//! increasing the window size, but we must be very careful to
			//! only do this if the transient is significant enough - when
			//! it's not, increasing the window size will cause 'metallic'
			//! artifacts around the transient, as well as post-echo.
			if(MaxRatio-TransientRatio < 0x1.62E430p-1f) break; //! 0x1.62E430p-1 = Log[2.0]

			//! Set new decimation pattern + transient ratio, and continue
			Decimation     = nSegments + MaxSegment;
			TransientRatio = MaxRatio;
			if(nSegments > 1) {
				//! Increase spectral resolution
				nSegments   /= 2;
				SegmentSize *= 2;
			} else break;
		}
	}
	if(TransientRatio == 0.0f) {
		//! When we don't detect a transient during the first pass,
		//! Log2SubBlockSize is incorrect, so fix it up here. This
		//! is still cleaner code than the alternative :/
		Log2SubBlockSize = 31 - __builtin_clz(BlockSize);
		Log2SubBlockSize++; //! <- Account for pre-increment
	}

	//! Determine overlap size from the ratio
	//! NOTE: Log2SubBlockSize was pre-incremented earlier, so account for this here.
	TransientRatio *= 0x1.715476p0f; //! 0x1.715476p0 = 1/Log[2] for change of base
	int OverlapScale = (TransientRatio < 0.5f) ? 0 : (TransientRatio >= 6.5f) ? 7 : (int)(TransientRatio+0.5f);
	if(Log2SubBlockSize-OverlapScale < 5+1) OverlapScale = Log2SubBlockSize - (5+1); //! Minimum 32-sample overlap

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
