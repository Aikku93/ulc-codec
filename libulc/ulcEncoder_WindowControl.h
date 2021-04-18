/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <math.h>
/**************************************/
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
//! NOTE: StepBuffer must be 2*BlockSize in size.
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
//!  Transient subblocks are thus conveniently indexed via
//!  POPCNT (minus 1 to remove the unary count 'stop' bit)
struct Block_Transform_GetWindowCtrl_TransientFiltering_Sum_t {
	float Sum, SumW;
};
static inline void Block_Transform_GetWindowCtrl_TransientFiltering(
	const float *Data,
	const float *LastBlockData,
	      float *StepBuffer,
	int BlockSize,
	int nChan
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
	//! NOTE: We end up missing the first sample of the old
	//! block, and the last sample of the new block; this really
	//! shouldn't affect things, though.
	//! NOTE: We decimate by a factor of 2 to save a bit of memory
	//! and because the results should be almost the same anyway.
	//! NOTE: BPFILT() accepts z^-1,1,z^1 for flexibility if we
	//! ever need to change the filter formula.
	for(n=0;n<BlockSize;n++) StepBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
#define BPFILT(zM1, z0, z1) ((z1) - (zM1))
		float *Dst = StepBuffer;
		const float *SrcOld = LastBlockData + Chan*BlockSize;
		const float *SrcNew = Data          + Chan*BlockSize;
		{
			*Dst   += 0.0f; //! H(z) = (z^1 - z^-1) becomes 0 with even symmetry about n=1/2
			*Dst++ += SQR(BPFILT(SrcOld[0], SrcOld[1], SrcOld[2])), SrcOld++;
		}
		for(n=1;n<BlockSize/2-1;n++) {
			*Dst   += SQR(BPFILT(SrcOld[0], SrcOld[1], SrcOld[2])), SrcOld++;
			*Dst++ += SQR(BPFILT(SrcOld[0], SrcOld[1], SrcOld[2])), SrcOld++;
		}
		{
			*Dst   += SQR(BPFILT(SrcOld[0], SrcOld[1], SrcOld[2])), SrcOld++;
			*Dst++ += SQR(BPFILT(SrcOld[0], SrcOld[1], SrcNew[0])), SrcOld++;
			*Dst   += SQR(BPFILT(SrcOld[0], SrcNew[0], SrcNew[1]));
			*Dst++ += SQR(BPFILT(SrcNew[0], SrcNew[1], SrcNew[2])), SrcNew++;
		}
		for(n=1;n<BlockSize/2-1;n++) {
			*Dst   += SQR(BPFILT(SrcNew[0], SrcNew[1], SrcNew[2])), SrcNew++;
			*Dst++ += SQR(BPFILT(SrcNew[0], SrcNew[1], SrcNew[2])), SrcNew++;
		}
		{
			*Dst   += SQR(BPFILT(SrcNew[0], SrcNew[1], SrcNew[2])), SrcNew++;
			*Dst++ += 0.0f; //! H(z) = (z^1 - z^-1) becomes 0 with even symmetry about n=N-1/2
		}
#undef BPFILT
	}

	//! Compute back-propagated energy (ie. how much energy we
	//! can tolerate as part of pre-echo).
	//! NOTE: We increase the signal gain here so that we are
	//! less likely to run into underflow isuses with exp-decay.
	//! To avoid adding an extra multiplication just to increase
	//! the gain, we fold it into OneMinusDecay for both backwards
	//! and forwards spreading calculations.
	//! Note that if RescaleGain is changed, the bias in the
	//! FINALIZE_DATA() macro must be updated (we rely on the
	//! step data being normalized at that stage).
	//! NOTE: While we're at it, we also apply the square root
	//! needed to get the RMS energy over all channels.
	float RescaleGain = 0x1.0p24f / (2.0f * 2.0f * nChan); //! *0.5 for decimation, *0.5 for BP gain
	float *BackMasking = StepBuffer + BlockSize; {
		//! Decay curve: -0.5dB/sample
		//! Calculation:
		//!  10^(-dBPerSample/20 * 2) = E^(Log[10]*(-dBPerSample/20) * 2)
		//! NOTE: Exponent multiplied by 2 to account for decimation by 2.
		float *Src = StepBuffer  + BlockSize; //! Iterated backwards
		float *Dst = BackMasking + BlockSize; //! Iterated backwards
		float  Tap = 0.0f, Decay = 0x1.C8520Bp-1f, OneMinusDecay = 0x1.BD6FA8p-4f * RescaleGain;
		for(n=0;n<BlockSize;n++) {
			float v = sqrtf(*--Src); *Src = v;
			Tap += v * OneMinusDecay;
			*--Dst = Tap;
			Tap *= Decay;
		}
	}

	//! Compute forward-propagated energy (ie. how much energy we
	//! can tolerate as part of post-echo), and differentiate with
	//! respect to backwards-propagated energy to capture changes
	//! in the signal envelope for analysis. This differentiated
	//! energy is then put into a smooth-max/entropy sum for
	//! accelerating analysis in window-switch/overlap-scale.
	{
		//! Decay curve: -0.2dB/sample
		int AnalysisIntervalMask = BlockSize/(ULC_HELPER_SUBBLOCK_INTERLEAVE_MODULO*4) - 1; //! Break up into LL/L/M/R (*4)
		const float *Src = StepBuffer;
		      float  Tap = 0.0f, Decay = 0x1.E8F4CAp-1f, OneMinusDecay = 0x1.70B363p-5f * RescaleGain;
		struct Block_Transform_GetWindowCtrl_TransientFiltering_Sum_t *Dst = (void*)StepBuffer; //! void* cast avoids a nasty, long typecast
		struct Block_Transform_GetWindowCtrl_TransientFiltering_Sum_t  Tmp = {.Sum = 0.0f, .SumW = 0.0f};
		for(n=0;n<BlockSize;n++) {
			Tap += *Src++ * OneMinusDecay;
			float d = *BackMasking++ - Tap;
			Tap *= Decay;

			//! Accumulate to sums.
			//! Because everything would be summed up in the search loop
			//! of Block_Transform_GetWindowCtrl(), we sum as much as we
			//! can here to reuse as many computations as possible.
			//! NOTE: These sums bear some similarity to Shannon entropy,
			//! but I'm not sure of the connection. This formulation seems
			//! to give the best results at any block size, though.
			//! NOTE: Even though we are using non-linear weights, the
			//! sum normalizes to a linear logarithm (this is important
			//! in the threshold analyses). Also, taking the logarithm of
			//! the original value should preserve accuracy a bit better.
			//! NOTE: We shouldn't need very precise logarithms here, so
			//! just use a cheaper approximation.
			float w = SQR(SQR(d));
			Tmp.SumW += w, Tmp.Sum += w*ULC_FastLnApprox(ABS(d));

			//! Wrapping around to next segment?
			if(((n+1) & AnalysisIntervalMask) == 0) *Dst++ = Tmp, Tmp.Sum = Tmp.SumW = 0.0f;
		}
	}
}
static inline int Block_Transform_GetWindowCtrl(
	const float *Data,
	const float *LastBlockData,
	      float *StepBuffer,
	int BlockSize,
	int nChan
) {
	int n;

	//! Perform filtering to obtain pre-echo analysis
	//! NOTE: Output data is stored as a struct to improve performance (SIMD
	//! optimizable, in theory). The void* cast avoids a nasty, long typecast.
	Block_Transform_GetWindowCtrl_TransientFiltering(Data, LastBlockData, StepBuffer, BlockSize, nChan);
	struct Block_Transform_GetWindowCtrl_TransientFiltering_Sum_t *TransientData = (void*)StepBuffer;

	//! Begin binary search for transient segment until it stops on the R side,
	//! at which point the ratio for the transition region is stored
	float Ratio;
	int Decimation    = 0b0001;
	int SubBlockSize  = BlockSize;
	int AnalysisLen  = ULC_HELPER_SUBBLOCK_INTERLEAVE_MODULO;
	for(TransientData += AnalysisLen;;) { //! MDCT transition region begins -BlockSize/2 samples from the new block (ie. L segment, in LL/L/M/R notation)
		//! Find the peak ratio within each segment (L/M/R),
		//! making sure to store the transition region ratio
		enum { POS_L, POS_M, POS_R};
		int   RatioPos;
		float RatioR; {
			//! This is used as a placeholder for Log[0]
			const float MIN_LOG = -100.0f;

			//! Get the energy of each segment (LL/L/M/R)
			struct Block_Transform_GetWindowCtrl_TransientFiltering_Sum_t LL, L, M, R;
			LL.Sum = LL.SumW = 0.0f;
			L .Sum = L .SumW = 0.0f;
			M .Sum = M .SumW = 0.0f;
			R .Sum = R .SumW = 0.0f;
			for(n=0;n<AnalysisLen;n++) {
#define SUM_DATA(Dst, Src) Dst.Sum  += Src.Sum, Dst.SumW += Src.SumW
				SUM_DATA(LL, TransientData[-1*AnalysisLen + n]);
				SUM_DATA(L,  TransientData[ 0*AnalysisLen + n]);
				SUM_DATA(M,  TransientData[+1*AnalysisLen + n]);
				SUM_DATA(R,  TransientData[+2*AnalysisLen + n]);
#undef SUM_DATA
			}
#define FINALIZE_DATA(x) x.Sum = (x.Sum != 0.0f) ? (x.Sum / x.SumW - 0x1.0A2B24p4f) : MIN_LOG //! -0x1.0A2B24p4 = Log[2^-24]; undo rescaling
			FINALIZE_DATA(LL);
			FINALIZE_DATA(L);
			FINALIZE_DATA(M);
			FINALIZE_DATA(R);
#undef FINALIZE_DATA
			//! Get the ratios between the segments
			float RatioL =   L.Sum - LL.Sum;
			float RatioM =   M.Sum - L .Sum;
			      RatioR = 2*R.Sum - M .Sum;

			//! Select the largest ratio of L/M/R
			                   RatioPos = POS_L, Ratio = RatioL;
			if(RatioM > Ratio) RatioPos = POS_M, Ratio = RatioM;
			if(RatioR > Ratio) RatioPos = POS_R, Ratio = RatioR;
		}

		//! Can we decimate?
		//! NOTE: Minimum subblock size of 64 samples.
		//! NOTE: Checking AnalysisLen should be better than checking
		//! Decimation directly, as then we can change the maximum allowed
		//! decimation without changing this code.
		if(ULC_USE_WINDOW_SWITCHING && AnalysisLen > 1 && SubBlockSize > 64) {
			//! If the transient is not in the transition region and
			//! is still significant, decimate the subblock further
			if(RatioPos != POS_R && Ratio >= 0x1.62E430p-2f) { //! Log[Sqrt[2]]
				//! Update the decimation pattern and continue
				if(RatioPos == POS_L)
					Decimation     = (Decimation<<1) | 0;
				else
					Decimation     = (Decimation<<1) | 1,
					TransientData += AnalysisLen;
				AnalysisLen  /= 2;
				SubBlockSize /= 2;
				continue;
			}
		}

		//! No more decimation - break out of the decimation loop
		Ratio = RatioR;
		break;
	}

	//! Determine the overlap scaling for the transition region
	int OverlapScale = 0;
	if(Ratio >= 0x1.62E430p-2f) {     //! Ratio >= Log[2^0.5]
		if(Ratio < 0x1.205967p2f) //! Ratio <  Log[2^6.5]
			OverlapScale = (int)(0x1.715476p0f*Ratio + 0.5f); //! 1/Log[2]
		else
			OverlapScale = 7;
		while((SubBlockSize >> OverlapScale) < 16) OverlapScale--; //! Minimum 16-sample overlap
	}

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
