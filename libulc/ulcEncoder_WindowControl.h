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
	float *LastTransientEnergy,
	float *StepBuffer,
	int BlockSize,
	int nChan
) {
	int n, Chan;

	//! Copy the transient energy from the previous block
	for(n=0;n<BlockSize;n++) StepBuffer[n] = LastTransientEnergy[n];

	//! Perform a bandpass filter to isolate the energy that is
	//! important to transient detection. Generally, LF energy
	//! and HF energy are 'unimportant', and it's the MF energy
	//! that has most of the information we're interested in.
	//! Transfer function:
	//!  H(z) = z^1 - z^-1
	//! NOTE: This filter does not have unity gain, as doing so
	//! would add some multiplications that reduce performance.
	//! NOTE: Recompute the sample at the boundary between the
	//! previoos block and this one, as it was "wrong" last time.
	for(n=BlockSize-1;n<BlockSize*2;n++) StepBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
#define BPFILT(zM1, z1) ((z1) - (zM1))
		float *Dst = StepBuffer + BlockSize;
		const float *SrcOld = LastBlockData + Chan*BlockSize;
		const float *SrcNew = Data          + Chan*BlockSize;
		n = BlockSize-1;
		Dst[-1] += SQR(BPFILT(SrcOld[n-1], SrcNew[0])); //! Fix up last sample of previous block
		*Dst++  += SQR(BPFILT(SrcOld[n],   SrcNew[1]));
		for(n=1;n<BlockSize-1;n++) {
			*Dst++ += SQR(BPFILT(SrcNew[n-1], SrcNew[n+1]));
		}
		//*Dst++ += 0.0f; //! H(z) = (z^1 - z^-1) becomes 0 with even symmetry
#undef BPFILT
	}

	//! Filter the BP energy to isolate energy peaks.
	//! Doing this tends to isolate energy spikes a lot more cleanly
	//! (less glitching from false positives) and also takes care of
	//! decay transients simultaneously, avoiding the need for
	//! special handling later.
	//! NOTE: The outputs of the filtered energy are passed through
	//! a very janky expander. Despite how janky it looks, it
	//! actually improves coding quality tremendously by being very
	//! sensitive to real transients, while ignoring non-transients.
	//! The "Norm" term is just a normalization factor to avoid
	//! having the numbers blowing up to +inf; I'm not sure what it
	//! really should be, but the value it's set to was derived
	//! from experimenting and should mostly stay within [0,10].
	//! Hopefully this will be improved in future.
	{
		float Gain = 0.0f;
		float Norm = SQR(SQR(0.25f / nChan));
		float *Buf = StepBuffer + BlockSize;
		float  Tap = sqrtf(Buf[-1]);
		Buf[-1] = 0.0f; //! We don't know what Gain or Buf[-2] was, so set the value to 0
		for(n=0;n<BlockSize-1;n++) {
			float v = sqrtf(Buf[n]);
			Buf[n] = SQR(v - Tap) * Gain;
			Tap = v;
			Gain = 0.75f*Gain + Buf[n] + Norm;
		}
		Buf[n] = 0.0f; //! Final BP sample is always 0, so this sample is unreliable - set it to 0
	}

	//! Copy the new transient energy back to the intermediate buffer
	for(n=0;n<BlockSize;n++) LastTransientEnergy[n] = StepBuffer[BlockSize + n];

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
			//! Get the energy of each segment (LL/L/M/R), making
			//! sure to window it in such a way that most the most
			//! importance is given to values /farthest/ from the
			//! corresponding transition region. This is because
			//! values close to the transition region "don't matter"
			//! as they can be handled by changing the overlap amount,
			//! whereas the farther they are, the more pre-echo will
			//! become audible.
			const float *Sin = Fourier_SinTableN(SubBlockSize_2);
			const float Bias = 0x1.0p-63f;
			float LL = Bias;
			float L  = Bias;
			float M  = Bias;
			float R  = Bias;
			for(n=0;n<SubBlockSize_2;n++) {
				float s  = SQR(Sin[n]);
				float ll = StepBuffer[n - SubBlockSize_2];
				float l  = StepBuffer[n];
				float m  = StepBuffer[n + SubBlockSize_2];
				float r  = StepBuffer[n + SubBlockSize_2*2];
				LL += s*ll;
				L  += s*l;
				M  += s*m;
				R  += s*r;
			}

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
			if(RatioPos != POS_R && Ratio >= 2.0f) {
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

	//! Get the final overlapping ratio. Unlike the decimation
	//! window from earlier, this part places importance on the
	//! energy closest to the transition region, as anything too
	//! far will pre-echo no matter what.
	{
		const float *Sin = Fourier_SinTableN(SubBlockSize_2);
		const float *Cos = Sin + SubBlockSize_2;
		const float Bias = 0x1.0p-63f;
		float L = Bias;
		float R = Bias;
		for(n=0;n<SubBlockSize_2;n++) {
			float s = SQR(Sin[n]);
			float c = SQR(Cos[-1-n]);
			float l = StepBuffer[n];
			float m = StepBuffer[n + SubBlockSize_2];
			float r = StepBuffer[n + SubBlockSize_2*2];
			L += s*l + c*m;
			R += s*m + c*r;
		}
		RatioR = R / L;
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
