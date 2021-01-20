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
static inline void Block_Transform_GetWindowCtrl_TransientFiltering(
	const float *Data,
	const float *LastBlockData,
	float *LastTransientEnergy,
	float *StepBuffer,
	float *CompressorGain,
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
	//! This will be compensated for later.
	//! NOTE: Recompute the sample at the boundary between the
	//! previous block and this one, as it wasn't available in
	//! the previous block (and was set to 0).
	for(n=BlockSize;n<BlockSize*2;n++) StepBuffer[n] = 0.0f;
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

	//! Filter the BP energy to extract transient spikes.
	//! The idea is to take the difference between the
	//! instantaneous energy, and the integrated energy
	//! (using a leaky integrator). The main point of
	//! interest here is the final step:
	//!  Buf[n] = d^8
	//! This exponentiation corresponds to a dynamic-range
	//! expander (1:8 ratio), and what it achieves is
	//! to lower the noise floor (corresponding to minor
	//! variations in energy and 'soft' transients) so that
	//! non-transient events are less likely to trigger at
	//! a given threshold, while improving the resolution
	//! of 'real' transients by separating them cleanly
	//! from the noise floor, allowing us to use a much
	//! larger threshold than expected. This large threshold
	//! is key to avoiding triggering during non-transient
	//! events, which maintains high-quality output.
	//! NOTE: I'm really not sure what to set the gain to.
	//! Technically, there's no need for any renormalization,
	//! but this filter can generate huge values that could
	//! overflow during subsequent calculations.
	{
		const float Decay = 0x1.C8520Bp-1f;                 //! Smoothness of SmoothedEnergy; -1.0dB/sample (10^(-1.0/20))
		      float Norm  = 0x1.8386AFp-7f * 0.25f / nChan; //! (1-Decay)^2 * BandpassEnergyGain^-1. Squaring cancels after the square root
		float *Buf = StepBuffer;
		float SmoothGain = *CompressorGain;
		for(n=BlockSize-1;n<BlockSize*2-1;n++) {
			float v = sqrtf(Buf[n] * Norm);
			float d = v - SmoothGain;
			SmoothGain = SmoothGain*Decay + v;
			Buf[n] = SQR(SQR(SQR(d)));
		}
		*CompressorGain = SmoothGain;
		Buf[n] = 0.0f; //! Last sample is unavailable - Exclude from analysis
	}

	//! Save the new transient energy back to the caching buffer
	for(n=0;n<BlockSize;n++) LastTransientEnergy[n] = StepBuffer[BlockSize + n];
}
static inline int Block_Transform_GetWindowCtrl(
	const float *Data,
	const float *LastBlockData,
	float *LastTransientEnergy,
	float *StepBuffer,
	float *CompressorGain,
	int BlockSize,
	int nChan
) {
	int n;

	//! Perform filtering to enhance transient analysis
	Block_Transform_GetWindowCtrl_TransientFiltering(Data, LastBlockData, LastTransientEnergy, StepBuffer, CompressorGain, BlockSize, nChan);

	//! Begin binary search for transient segment until it stops on the R side,
	//! at which point the ratio for the transition region is stored
	int Decimation     = 0b0001;
	int SubBlockSize_2 = BlockSize/2;
	for(StepBuffer += SubBlockSize_2;;) { //! MDCT transition region begins -BlockSize/2 samples from the new block
		//! Find the peak ratio within each segment (L/M/R)
		enum { POS_L, POS_M, POS_R};
		int   RatioPos;
		float Ratio; {
			//! This is used as a placeholder for Log[0]
			const float MIN_LOG = -100.0f;

			//! Get the energy of each segment (LL/L/M/R)
			//! NOTE: We are computing a log-domain smooth-max here,
			//! which has some relation to Shannon entropy. I'm not
			//! entirely sure what the connection is, but this seems
			//! to give the best results even for large block sizes.
			float LL = 0.0f, LLw = 0.0f;
			float L  = 0.0f, Lw  = 0.0f;
			float M  = 0.0f, Mw  = 0.0f;
			float R  = 0.0f, Rw  = 0.0f;
			for(n=0;n<SubBlockSize_2;n++) {
				float ll = StepBuffer[n - SubBlockSize_2];
				float l  = StepBuffer[n];
				float m  = StepBuffer[n + SubBlockSize_2];
				float r  = StepBuffer[n + SubBlockSize_2*2];
				if(ll != 0.0f) LLw += ll, LL += ll * logf(ll);
				if(l  != 0.0f) Lw  += l,  L  += l  * logf(l);
				if(m  != 0.0f) Mw  += m,  M  += m  * logf(m);
				if(r  != 0.0f) Rw  += r,  R  += r  * logf(r);
			}
			LL = (LL != 0.0f) ? (LL / LLw) : MIN_LOG;
			L  = (L  != 0.0f) ? (L  / Lw)  : MIN_LOG;
			M  = (M  != 0.0f) ? (M  / Mw)  : MIN_LOG;
			R  = (R  != 0.0f) ? (R  / Rw)  : MIN_LOG;

			//! Get the ratios between segments
			float RatioL = L - LL;
			float RatioM = M - L;
			float RatioR = R - M;

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
			if(RatioPos != POS_R && Ratio >= 0x1.62E430p1f) { //! Log[16.0]
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
	int OverlapScale = 0; {
		//! Perform a final pass over the transition region,
		//! this time with the aim of finding how concentrated
		//! the transient energy is about the transition region.
		//! NOTE: We really do care more about what's in the
		//! transition region itself, so we use these curves:
		//!  In  = {1-Cos,Cos}
		//!  Out = {Cos,1-Cos}
		//! This sometimes results in popping, in which case
		//! the following curves may be used:
		//!  In  = {Sin^2,Cos}
		//!  Out = {Cos,Sin^2}
		//! However, this variation results in less sharpness.
		float In  = 0.0f;
		float Out = 0.0f;
		const float *Sin = Fourier_SinTableN(SubBlockSize_2);
		const float *Cos = Sin + SubBlockSize_2;
		for(n=0;n<SubBlockSize_2;n++) {
			float c = *--Cos;
			float s = 1.0f - c;
			float l = StepBuffer[n + SubBlockSize_2];
			float r = StepBuffer[n + SubBlockSize_2*2];
			In  += s*l + c*r;
			Out += c*l + s*r;
		}

		//! Finally determine the scaling from the ratio of In:Out
		if(In*0x1.6A09E6p-1f >= Out) {      //! In/Out >= 2^0.5
			if(In*0x1.6A09E6p-7f < Out) //! In/Out <  2^6.5
				OverlapScale = (int)(0x1.715476p0f*logf(In/Out) + 0.5f); //! 1/Log[2] for change-of-base
			else
				OverlapScale = 7;
			while((SubBlockSize_2 >> OverlapScale) < 16/2) OverlapScale--; //! Minimum 16-sample overlap
		}
	}

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(Decimation != 1) + 0x10*Decimation;
}

/**************************************/
//! EOF
/**************************************/
