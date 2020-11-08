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
static inline float Block_Transform_GetWindowCtrl_TransientRatio(float R, float L, float DecayScale) {
	if(DecayScale*R < L) { float t = L; L = DecayScale*R, R = t; }
	return R / L;
}
static inline int Block_Transform_GetWindowCtrl(
	const float *Data,
	const float *LastBlockData,
	float *StepBuffer,
	int BlockSize,
	int MinOverlap,
	int MaxOverlap,
	int nChan,
	int RateHz
) {
	int n, Chan;

	//! Perform a highpass filter to get the step/transient energies
	//! Transfer function:
	//!  H(z) = 1 - z^-1
	//! Assumes an even-symmetric structure at the boundaries
	//! NOTE: This filter does not have unity gain, as doing so
	//! would add some multiplications that would reduce performance.
	//! NOTE: The following filter was attempted:
	//!  H(z) = -0.5z^1 + 1 - 0.5z^-1
	//! However, this gave very, very poor results with the method
	//! of transient detection used at present.
	for(n=0;n<2*BlockSize;n++) StepBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
#define HPFILT(zM1, z0, z1) (-(zM1) + (z0))
		float *Dst = StepBuffer;
		const float *SrcOld = LastBlockData + Chan*BlockSize;
		const float *SrcNew = Data          + Chan*BlockSize;

		//! Get step energy for last block
		*Dst++ += SQR(HPFILT(SrcOld[1], SrcOld[0], SrcOld[1]));
		for(n=1;n<BlockSize-1;n++) {
			*Dst++ += SQR(HPFILT(SrcOld[n-1], SrcOld[n], SrcOld[n+1]));
		}
		*Dst++ += SQR(HPFILT(SrcOld[n-1], SrcOld[n], SrcNew[0]));

		//! Get step energy for this block
		*Dst++ += SQR(HPFILT(SrcOld[n], SrcNew[0], SrcNew[1]));
		for(n=1;n<BlockSize-1;n++) {
			*Dst++ += SQR(HPFILT(SrcNew[n-1], SrcNew[n], SrcNew[n+1]));
		}
		*Dst++ += SQR(HPFILT(SrcNew[n-1], SrcNew[n], SrcNew[n-1]));
#undef HPFILT
	}

	//! Spread the energy forwards a bit, so that transients at the
	//! edge of a transition zone are more likely to be detected, as
	//! well as masking subtle jitter in the data.
	//! Note that overdoing this (ie. setting the decay curve too
	//! low) will miss some transients, and not doing it enough (ie.
	//! setting it too high) will be overly sensitive.
	//! TL;DR: This gives transients extra oomph, smooths out
	//! irregularities in the signal, and does post-masking.
	{
		//! Transient decay curve: 3dB/ms (in X^2 domain):
		//!  E^(-Log[10]*1000*(dBDecayPerMs/10) / RateHz)
		float *Buf = StepBuffer;
		float a = expf(-0x1.596344p9f / RateHz);
		float SampleLast = *Buf++;
		for(n=1;n<BlockSize*2;n++) {
			float v = *Buf;
			*Buf++ = SampleLast = a*SampleLast + v;
		}
	}

	//! Window off the tail end of the sample data
	//! This avoids doubling up on the transient by overlap
	//! adjustment on one block, and decimation on the next.
	{
		const float *Sin = Fourier_SinTableN(BlockSize/4);
		float *Buf = StepBuffer + 2*BlockSize;
		for(n=0;n<BlockSize/4;n++) *--Buf *= Sin[n];
	}

	//! Finally, convert the filtered energy into transient
	//! ratio segments for detecting spikes within the subblocks.
	//! From thorough experimentation, subdividing into small
	//! segments from which we extract ratios seems to result
	//! in the best transient detection.
	//! NOTE: We don't necessarily use all of SegmentBuffer later
	//! (meaning there is unnecessary calculations here), but it
	//! appears to not be a bottleneck at present, so this has
	//! not been optimized out.
	int nRatioSegments;
	float *SegmentBuffer = StepBuffer; {
		const int RatioSegmentSamples = 64;
		nRatioSegments = 2*BlockSize / RatioSegmentSamples;

		//! Get the ratios for each segment
		//! NOTE: Adding a bias avoids issues with x/0.
		//! NOTE: Everything here works much better with squared ratios.
		//! NOTE: Cap the ratio at sane limits.
		int Segment;
		const float Bias = 0x1.0p-48f;
		float L, R = Bias;
		const float *Src = StepBuffer;
		for(Segment=0;Segment<nRatioSegments;Segment++) {
			L = R;
			R = Bias;
			for(n=0;n<RatioSegmentSamples;n++) R += *Src++;
			float Ratio = Block_Transform_GetWindowCtrl_TransientRatio(R, L, SQR(4.0f));
			Ratio *= Ratio;
			SegmentBuffer[Segment] = (Ratio < 1.0f) ? 1.0f : (Ratio > 0x1.0p31f) ? 0x1.0p31f : Ratio;
		}

		//! Offset by -BlockSize/2 relative to the start of
		//! the data given for the block being analyzed. This
		//! is because MDCT places the overlap region here:
		//!  +nRatioSegments/2 places us on the next block
		//!  -nRatioSegments/4 rewinds to the transition region
		//!  Adding together gives the offset +nRatioSegments/4
		SegmentBuffer += nRatioSegments/4;
	}

	//! Begin binary search for transient segment until it stops on the R side,
	//! at which point the ratio for the transition region is stored
	int Decimation = 1, SegmentSize = nRatioSegments / 4, SubBlockSize = BlockSize;
	float Ratio;
	for(;;) {
		//! Find the peak ratio within each segment (L/M/R)
		enum { POS_L, POS_M, POS_R};
		int   RatioPos;
		float RatioL;
		float RatioM;
		float RatioR; {
			//! Find the smooth-max of each segment
			float L = 0.0f, LW = 0.0f;
			float M = 0.0f, MW = 0.0f;
			float R = 0.0f, RW = 0.0f;
			for(n=0;n<SegmentSize;n++) {
				float l = SegmentBuffer[n];
				float m = SegmentBuffer[n + SegmentSize];
				float r = SegmentBuffer[n + SegmentSize*2];
				L += SQR(l), LW += l;
				M += SQR(m), MW += m;
				R += SQR(r), RW += r;
			}

			//! Determine the segment's ratio by dividing
			//! the smooth-max ratio by the average ratio,
			//! and then select the largest of L/M/R.
			//! dividing by the average ratio helps to
			//! smooth out waveform irregularities (eg.
			//! false-positives on a saw wave, etc).
			//! NOTE: Remove the 'peak' from the average
			//! before dividing by it. This gives us an
			//! 'unbiased' average, by ignoring the peak.
			//! NOTE: This formula can be simplified, but
			//! is left as-is to avoid over/underflow.
			RatioL = L / LW;
			RatioM = M / MW;
			RatioR = R / RW;
			if(SegmentSize > 1) { //! With SegmentSize==1, this would cause division by 0
				RatioL = RatioL*SegmentSize / (LW - RatioL);
				RatioM = RatioM*SegmentSize / (MW - RatioM);
				RatioR = RatioR*SegmentSize / (RW - RatioR);
			}
			                   RatioPos = POS_L, Ratio = RatioL;
			if(RatioM > Ratio) RatioPos = POS_M, Ratio = RatioM;
			if(RatioR > Ratio) RatioPos = POS_R, Ratio = RatioR;
		}

		//! Can we decimate?
		if(ULC_USE_WINDOW_SWITCHING && SegmentSize >= 1 && Decimation < 0x8 && SubBlockSize > MinOverlap) {
			//! If the transient is not in the transition region and
			//! is still significant, decimate the subblock further
			if(RatioPos != POS_R && Ratio > 1.25f) {
				//! Update the decimation pattern and continue
				if(RatioPos == POS_L)
					Decimation     = (Decimation<<1) | 0;
				else
					Decimation     = (Decimation<<1) | 1,
					SegmentBuffer += SegmentSize;
				SegmentSize  /= 2;
				SubBlockSize /= 2;
				continue;
			}
		}

		//! No more decimation - set the ratio to that on the R side and break out
		Ratio = RatioR;
		break;
	}

	//! Determine the overlap amount for the transition region
	int OverlapScale;
	if(Ratio >= 0x1.6A09E6p0f) { //! Ratio==2^0.5 is the first ratio to result in OverlapScale > 0
		//! Ratio >= 2^6.5 gives the upper limit for OverlapScale = 7
		if(Ratio < 0x1.6A09E6p6f) {
			//! 0x1.715476p0 = 1/Log[2], to get the log base-2
			OverlapScale = (int)(0x1.715476p0f*logf(Ratio) + 0.5f);
		} else OverlapScale = 7;
		while(OverlapScale > 0 && (SubBlockSize >> OverlapScale) < MinOverlap) OverlapScale--;
		while(OverlapScale < 7 && (SubBlockSize >> OverlapScale) > MaxOverlap) OverlapScale++;
	} else OverlapScale = 0;

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
