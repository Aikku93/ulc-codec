/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
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
static inline int Block_Transform_GetWindowCtrl(
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
	//! Assumes an even-symmetric structure at the boundaries
	//! NOTE: This filter does not have unity gain, as doing so
	//! would add some multiplications that reduce performance.
	for(n=0;n<2*BlockSize;n++) StepBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
#define BPFILT(zM1, z0, z1) ((z1) - (zM1))
		float *Dst = StepBuffer;
		const float *SrcOld = LastBlockData + Chan*BlockSize;
		const float *SrcNew = Data          + Chan*BlockSize;

		//! Get step energy for last block
		*Dst++ += SQR(BPFILT(SrcOld[1], SrcOld[0], SrcOld[1]));
		for(n=1;n<BlockSize-1;n++) {
			*Dst++ += SQR(BPFILT(SrcOld[n-1], SrcOld[n], SrcOld[n+1]));
		}
		*Dst++ += SQR(BPFILT(SrcOld[n-1], SrcOld[n], SrcNew[0]));

		//! Get step energy for this block
		*Dst++ += SQR(BPFILT(SrcOld[n], SrcNew[0], SrcNew[1]));
		for(n=1;n<BlockSize-1;n++) {
			*Dst++ += SQR(BPFILT(SrcNew[n-1], SrcNew[n], SrcNew[n+1]));
		}
		*Dst++ += SQR(BPFILT(SrcNew[n-1], SrcNew[n], SrcNew[n-1]));
#undef BPFILT
	}

	//! Filter the BP energy to isolate energy peaks.
	//! Doing this tends to isolate energy spikes a
	//! lot more cleanly (less glitching) and also
	//! takes care of decay transients simultaneously,
	//! avoiding the need for special handling later.
	//! NOTE: The first output of the bandpass filter
	//! is always 0 due to the even-symmetric structure
	//! used. This is incorrect, but is convenient in
	//! that it avoids having to store even more data
	//! about the previous block. However, this also
	//! means that the first two outputs of this delta
	//! filter's output are "wrong" and must be set to
	//! 0 so that they don't mess up the analysis.
	{
		float *Buf = StepBuffer;
		float  Tap = sqrtf(Buf[1]);
		Buf[0] = 0.0f;
		Buf[1] = 0.0f;
		for(n=2;n<BlockSize*2;n++) {
			float v = sqrtf(Buf[n]);
			Buf[n] = SQR(v - Tap), Tap = v;
		}
	}

	//! Begin binary search for transient segment until it stops on the R side,
	//! at which point the ratio for the transition region is stored
	float RatioL;
	float RatioM;
	float RatioR;
	int   Decimation     = 0b0001;
	int   SubBlockSize_2 = BlockSize/2;
	for(StepBuffer += SubBlockSize_2;;) { //! MDCT transition region begins -BlockSize/2 samples from the new block
		//! Find the peak ratio within each segment (L/M/R)
		enum { POS_L, POS_M, POS_R};
		int   RatioPos;
		float Ratio; {
			//! Get the smooth-max of each segment (LL/L/M/R)
			const float Bias = 0x1.0p-63f, Bias2 = SQR(Bias);
			float LL = Bias2, LLw = Bias;
			float L  = Bias2, Lw  = Bias;
			float M  = Bias2, Mw  = Bias;
			float R  = Bias2, Rw  = Bias;
			for(n=0;n<SubBlockSize_2;n++) {
				float ll = StepBuffer[n - SubBlockSize_2];
				float l  = StepBuffer[n];
				float m  = StepBuffer[n + SubBlockSize_2];
				float r  = StepBuffer[n + SubBlockSize_2*2];
				LL += SQR(ll), LLw += ll;
				L  += SQR(l),  Lw  += l;
				M  += SQR(m),  Mw  += m;
				R  += SQR(r),  Rw  += r;
			}
			LL /= LLw;
			L  /= Lw;
			M  /= Mw;
			R  /= Rw;

			//! Get the ratios between segments
			RatioL = L / LL;
			RatioM = M / L;
			RatioR = R / M;

			//! Select the largest ratio of L/M/R
			                   RatioPos = POS_L, Ratio = RatioL;
			if(RatioM > Ratio) RatioPos = POS_M, Ratio = RatioM;
			if(RatioR > Ratio) RatioPos = POS_R, Ratio = RatioR;
		}

		//! Can we decimate?
		//! NOTE: Minimum subblock size of 64 samples.
		if(ULC_USE_WINDOW_SWITCHING && Decimation < 0x8 && SubBlockSize_2 > 64/2) {
			//! If the transient is not in the transition region and
			//! is still significant, decimate the subblock further
			if(RatioPos != POS_R && Ratio >= 1.5f) {
				//! Update the decimation pattern and continue
				if(RatioPos == POS_L)
					Decimation  = (Decimation<<1) | 0;
				else
					Decimation  = (Decimation<<1) | 1,
					StepBuffer += SubBlockSize_2;
				SubBlockSize_2 /= 2;
				continue;
			}
		}

		//! No more decimation - break out of the decimation loop
		break;
	}

	//! Determine the overlap scaling for the transition region
	int OverlapScale = 0;
	if(RatioR >= 0x1.6A09E6p0f) { //! Log2[RatioR] >= 0.5
		if(RatioR < 0x1.6A09E6p6f) //! Log2[RatioR] < 6.5
			OverlapScale = (int)(0x1.715476p0f*logf(RatioR) + 0.5f); //! 1/Log[2] for change-of-base
		else
			OverlapScale = 7;
		while((SubBlockSize_2 >> OverlapScale) < 16/2) OverlapScale--; //! Minimum 16-sample overlap
	}

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
