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
struct Block_Transform_GetWindowCtrl_TransientRatio_t {
	//! The weighted average is performed in the log domain for
	//! better numerical stability, and the weights are computed
	//! directly from the L and R values for the same reason.
	float wLog2Ratio; //! Weight*Log2[Ratio]
	float Weight;     //! Ratio-1.0
};
static inline struct Block_Transform_GetWindowCtrl_TransientRatio_t Block_Transform_GetWindowCtrl_TransientRatio(
	float  R,
	float  L,
	float  DecayScale,
	double RatioExponent
) {
	struct Block_Transform_GetWindowCtrl_TransientRatio_t Ret;
	if(DecayScale*R < L) { float t = L; L = DecayScale*R, R = t; }
	L += 0x1.0p-31f; //! NOTE: Adding a bias avoids issues with x/0

	//! NOTE: Cap the ratio at sane limits.
	if(R < L) {
		//! Soft decay - Ignore, and use Ratio=1.0
		Ret.wLog2Ratio = 0.0f;
		Ret.Weight     = 0.0f;
	} else {
		double r = pow((double)R / L, RatioExponent);
		if(r > 0x1.0p25) r = 0x1.0p25; //! Cap at a sane limit
		double w = r - 1.0;
		Ret.wLog2Ratio = (float)(w*log(r));
		Ret.Weight     = (float)w;
	}
	return Ret;
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

	//! Perform a bandpass filter to isolate the energy that is
	//! important to transient detection. Generally, LF energy
	//! and HF energy are 'unimportant', and it's the MF energy
	//! that has most of the information we're interested in.
	//! Transfer function:
	//!  H(z) = z^1 - z^-1
	//! Assumes an even-symmetric structure at the boundaries
	//! NOTE: This filter does not have unity gain, as doing so
	//! would add some multiplications that would reduce performance.
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

	//! Spread the energy forwards a bit, so that transients at the
	//! edge of a transition zone are more likely to be detected, as
	//! well as masking subtle jitter in the data.
	//! Note that overdoing this (ie. setting the decay curve too
	//! low) will miss some transients, and not doing it enough (ie.
	//! setting it too high) will be overly sensitive.
	//! TL;DR: This gives transients extra oomph, smooths out
	//! irregularities in the signal, and does post-masking.
	{
		//! Transient decay curve: 1.5dB/ms:
		//!  E^(-Log[10]*1000*(dBDecayPerMs/10) / RateHz)
		//! 1.5dB/ms is /massive/ overkill, but we really
		//! want to isolate just the first transient. This
		//! makes the transient threshold /very/ touchy,
		//! however, and will be extremely close to 1.0
		//! without exponentiating the ratios later.
		float *Buf = StepBuffer;
		float a = expf(-0x1.596344p8f / RateHz);
		float SampleLast = *Buf++;
		for(n=1;n<BlockSize*2;n++) {
			float v = *Buf;
			*Buf++ = SampleLast = a*SampleLast + v;
		}
	}

	//! Finally, convert the filtered energy into transient
	//! ratio segments for detecting spikes within the subblocks.
	//! From thorough experimentation, subdividing into small
	//! segments from which we extract ratios seems to result
	//! in the best transient detection.
	//! NOTE: We don't necessarily use all of RatioBuffer later
	//! (meaning there is unnecessary calculations here), but it
	//! appears to not be a bottleneck at present, so this has
	//! not been optimized out.
	//! TODO: Attempt using the old strategy of analyzing the
	//! total energy in each segment (LL/L/M/R) and dividing,
	//! but this time exponentiate the ratio using, eg.
	//! Log[nSamplesInSegment] (or since we're using log ratios,
	//! multiply it by Log[nSamplesInSegment]).
	int nRatioSegments;
	struct Block_Transform_GetWindowCtrl_TransientRatio_t *RatioBuffer = (struct Block_Transform_GetWindowCtrl_TransientRatio_t*)StepBuffer; {
		const int RatioSegmentSamples = 64;
		nRatioSegments = 2*BlockSize / RatioSegmentSamples;

		//! Get the ratios for each segment
		//! NOTE: The exponent was determined experimentally via
		//! testing and may not be the best choice.
		int Segment;
		float L, R = 0.0f;
		const float *Src = StepBuffer;
		for(Segment=0;Segment<nRatioSegments;Segment++) {
			L = R;
			R = 0.0f;
			for(n=0;n<RatioSegmentSamples;n++) R += *Src++;
			RatioBuffer[Segment] = Block_Transform_GetWindowCtrl_TransientRatio(R, L, SQR(4.0f), 20.0/6);
		}

		//! Offset by -BlockSize/2 relative to the start of
		//! the data given for the block being analyzed. This
		//! is because MDCT places the overlap region here:
		//!  +nRatioSegments/2 places us on the next block
		//!  -nRatioSegments/4 rewinds to the transition region
		//!  Adding together gives the offset +nRatioSegments/4
		RatioBuffer += nRatioSegments/4;
	}

	//! Begin binary search for transient segment until it stops on the R side,
	//! at which point the ratio for the transition region is stored
	int Decimation = 1, SegmentSize = nRatioSegments / 4, SubBlockSize = BlockSize;
	float LogRatio;
	for(;;) {
		//! Find the peak ratio within each segment (L/M/R)
		enum { POS_L, POS_M, POS_R};
		int   RatioPos;
		float RatioL;
		float RatioM;
		float RatioR; {
			//! Find the pseudo-peak ratio of each segment
			float L = 0.0f, LW = 0.0f;
			float M = 0.0f, MW = 0.0f;
			float R = 0.0f, RW = 0.0f;
			for(n=0;n<SegmentSize;n++) {
				struct Block_Transform_GetWindowCtrl_TransientRatio_t l, m, r;
				l = RatioBuffer[n];
				m = RatioBuffer[n + SegmentSize];
				r = RatioBuffer[n + SegmentSize*2];
				L += l.wLog2Ratio, LW += l.Weight;
				M += m.wLog2Ratio, MW += m.Weight;
				R += r.wLog2Ratio, RW += r.Weight;
			}
			RatioL = (LW > 0.0f) ? (L / LW) : 1.0f;
			RatioM = (MW > 0.0f) ? (M / MW) : 1.0f;
			RatioR = (RW > 0.0f) ? (R / RW) : 1.0f;

			//! Select the largest ratio of L/M/R
			                      RatioPos = POS_L, LogRatio = RatioL;
			if(RatioM > LogRatio) RatioPos = POS_M, LogRatio = RatioM;
			if(RatioR > LogRatio) RatioPos = POS_R, LogRatio = RatioR;
		}

		//! Can we decimate?
		if(ULC_USE_WINDOW_SWITCHING && SegmentSize > 1 && Decimation < 0x8 && SubBlockSize > MinOverlap) {
			//! If the transient is not in the transition region and
			//! is still significant, decimate the subblock further
			if(RatioPos != POS_R && LogRatio >= 1.0f) { //! Log2[2.0] = 1.0
				//! Update the decimation pattern and continue
				if(RatioPos == POS_L)
					Decimation   = (Decimation<<1) | 0;
				else
					Decimation   = (Decimation<<1) | 1,
					RatioBuffer += SegmentSize;
				SegmentSize  /= 2;
				SubBlockSize /= 2;
				continue;
			}
		}

		//! No more decimation - set the ratio to that on the R side and break out
		LogRatio = RatioR;
		break;
	}

	//! Determine the overlap scaling for the transition region
	int OverlapScale = 0;
	if(LogRatio >= 1.0f) {
		OverlapScale = (LogRatio < 6.5f) ? (int)(LogRatio + 0.5f) : 7;
		while(OverlapScale > 0 && (SubBlockSize >> OverlapScale) < MinOverlap) OverlapScale--;
		while(OverlapScale < 7 && (SubBlockSize >> OverlapScale) > MaxOverlap) OverlapScale++;
	}

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
