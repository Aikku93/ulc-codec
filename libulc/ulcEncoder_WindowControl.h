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
	//! previous block and this one, as it was "wrong" last time.
	//! Additionally, get the sample before that for the z^-1
	//! tap of the differentiator in the next section.
	float EnergyTap = 0.0f;
	StepBuffer[BlockSize-1] = 0.0f; //! Separating this one out might give better-optimized code
	for(n=BlockSize;n<BlockSize*2;n++) StepBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
#define BPFILT(zM1, z1) ((z1) - (zM1))
		float *Dst = StepBuffer + BlockSize;
		const float *SrcOld = LastBlockData + Chan*BlockSize;
		const float *SrcNew = Data          + Chan*BlockSize;
		n = BlockSize-1;
		EnergyTap += SQR(BPFILT(SrcOld[n-2], SrcOld[n])); //! Get energy of sample before last for a z^-1 tap
		Dst[-1]   += SQR(BPFILT(SrcOld[n-1], SrcNew[0])); //! Fix up last sample of previous block
		*Dst++    += SQR(BPFILT(SrcOld[n],   SrcNew[1]));
		for(n=1;n<BlockSize-1;n++) {
			*Dst++ += SQR(BPFILT(SrcNew[n-1], SrcNew[n+1]));
		}
		//*Dst++ += 0.0f; //! H(z) = (z^1 - z^-1) becomes 0 with even symmetry
#undef BPFILT
	}

	//! Filter the BP energy to isolate energy peaks, and then apply
	//! a strong compressor to the signal.
	//! Doing this tends to isolate energy spikes a lot more cleanly
	//! (less glitching from false positives) and also takes care of
	//! decay transients simultaneously, avoiding the need for
	//! special handling later.
	//! NOTE: Pushing the signal through a strong compressor with no
	//! release time and a short attack causes transients to pop out
	//! a lot more clearly than even just differentiation of the
	//! signal energy. However, this can also result in odd artifacts
	//! during non-transients. I've tried to tune the terms in the
	//! compressor to give the best tradeoff between sensitivity and
	//! accurate (true positive) detection, but there will always be
	//! cases where this will not work right. Further improvements
	//! are very much welcomed.
	{
		//! Attack coefficient:
		//!  Sensitivity^(-1/AttackSamples)
		//! ie. This solves with the number of samples it takes
		//! for Gain to return to unity after a transient pulse.
		//! Default: Sensitivity = 64, AttackSamples = 24 samples.
		//! NOTE: Sensitivity is kept non-const because it will be
		//! modified later for a more optimized processing loop.
		//! NOTE: AttackSamples essentially controls post-masking.
		//! However, this isn't "traditional" post-masking that
		//! lasts a while, but rather more of a jitter-prevention.
		//! Setting it too low will make the detector pick up too
		//! many false-positives, and setting it too high will
		//! make it harder to pick up transients (as well as cause
		//! the signal gain to explode a lot easier).
		      float Sensitivity        = 64.0f;
		const float AttackCoef         = 0x1.AE89FAp-1f;
		const float OneMinusAttackCoef = 0x1.45D81Ap-3f;

		//! NOTE: The output of the bandpass filter has a gain of
		//! 2^2*nChan. However, that same output is also strictly
		//! non-negative. So when we apply a square root and then
		//! differentiate before squaring, we are not changing the
		//! gain term at all, as the output of the differentiator
		//! can never be larger than its inputs. This means that
		//! the unity gain factor is still just 1/(2^2*nChan) at
		//! that point. However, this signal is being driven by
		//! Gain which is itself derived from the same inputs,
		//! and thus squares the values again, and so the final
		//! normalization factor becomes 1/(2^2*nChan)^2.
		//! The purpose of normalization is to avoid having the
		//! values exploding out of numerical represensation.
		//! NOTE: The above normalization is folded into UnityGain,
		//! but it must also be folded into Sensitivity to account
		//! for the signal feeding back into itself. This allows
		//! us to remove a multiplication from the processing loop.
		//! NOTE: The gain control expression can be expressed as:
		//!  Gain' = Gain*AttackCoef + UnityGain*(1.0 - AttackCoef) <- Form A
		//!        = Gain + (UnityGain - Gain)*(1.0 - AttackCoef)   <- Form B
		//! Form B makes it clear as to what we are doing here,
		//! but Form A is more numerically stable and efficient.
		//! NOTE: Because Gain can easily explode to infinity
		//! regardless of our safeguards in normalization, we must
		//! still limit it. Even though the range should ideally
		//! be 0.0 .. 1.0, we are enforcing a limit of 256.0, as
		//! this should be more than enough for any case.
		float GainNorm  = SQR(0.25f) / SQR(nChan);
		float GainLimit = GainNorm * 256.0f;
		float UnityGain = GainNorm * OneMinusAttackCoef;
		float Gain = *CompressorGain;
		float *Buf = StepBuffer;
		float  Tap = sqrtf(EnergyTap);
		Sensitivity *= GainNorm;
		for(n=BlockSize-1;n<BlockSize*2-1;n++) {
			float v = sqrtf(Buf[n]);
			float d = SQR(v - Tap);
			Buf[n] = d * Gain, Tap = v;
			Gain  = Gain*AttackCoef + UnityGain;
			Gain += d*Sensitivity;

			//! Gain can easily explode to infinity, so limit things here
			if(Gain > GainLimit) Gain = GainLimit;
		}
		Buf[n] = ABS(2*Buf[n-1] - Buf[n-2]); //! Last sample is unavailable - linearly extrapolate from its neighbours
		*CompressorGain = Gain;
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
			//! Minimum value for Log[x]. This is used as a
			//! placeholder for Log[0] in the analysis.
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
				LLw += ll; if(ll != 0.0f) LL += ll * logf(ll);
				Lw  += l;  if(l  != 0.0f) L  += l  * logf(l);
				Mw  += m;  if(m  != 0.0f) M  += m  * logf(m);
				Rw  += r;  if(r  != 0.0f) R  += r  * logf(r);
			}
			LL = (LL != 0.0f) ? (LL / LLw) : MIN_LOG;
			L  = (L  != 0.0f) ? (L  / Lw)  : MIN_LOG;
			M  = (M  != 0.0f) ? (M  / Mw)  : MIN_LOG;
			R  = (R  != 0.0f) ? (R  / Rw)  : MIN_LOG;

			//! Get the ratios between segments
			RatioL = L - LL;
			RatioM = M - L;
			RatioR = R - M;

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
			if(RatioPos != POS_R && Ratio >= 0x1.62E430p-1f) { //! Log[2.0]
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
	if(RatioR >= 0x1.62E430p-2f) { //! Log[2^0.5]
		if(RatioR < 0x1.205967p2f) //! Log[2^6.5]
			OverlapScale = (int)(0x1.715476p0f*RatioR + 0.5f); //! 1/Log[2] for change-of-base
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
