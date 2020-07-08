/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#if defined(__AVX__) || defined(__FMA__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include <stdint.h>
/**************************************/
#include "Fourier.h"
#include "ulcEncoder_Psycho.h"
#include "ulcHelper.h"
/**************************************/

//! Write the sort values for all coefficients in a block
//! and return the number of codeable non-zero coefficients
//! NOTE:
//!  AnalysisPower is used to alter the preference for the currently-being-analyzed channel
static inline void Block_Transform_WriteSortValues(
	float *CoefIdx,
	const float *Energy,
	const float *EnergyNp,
	const float *CoefNp,
	int   *nNzCoef,
	int   BlockSize,
	float AnalysisPowerNp,
	float NyquistHz
) {
	int Band;
#if ULC_USE_PSYCHOACOUSTICS
	float Gamma = logf(2*BlockSize);
	struct Block_Transform_MaskingState_t MaskingState;
	Block_Transform_MaskingState_Init(&MaskingState, Energy, EnergyNp, BlockSize, NyquistHz);
#else
	(void)Energy;
	(void)NyquistHz;
#endif
	for(Band=0;Band<BlockSize;Band++) {
		//! Inside the codeable range?
		float ValNp = CoefNp[Band];
		if(ValNp != ULC_COEF_NEPER_OUT_OF_RANGE) {
#if ULC_USE_PSYCHOACOUSTICS
			//! Apply psychoacoustic corrections to this band energy
			ValNp += Block_Transform_GetMaskedLevel(&MaskingState, Energy, EnergyNp, Band, BlockSize, Gamma);
#else
			//! Not so much a psychoacoustic optimization as an MDCT
			//! energy correction factor
			ValNp += EnergyNp[Band];
#endif
			//! Store the sort value for this coefficient
			CoefIdx[Band] = ValNp + AnalysisPowerNp;
			(*nNzCoef)++;
		} else CoefIdx[Band] = -1.0e30f; //! Unusable coefficient; map to the end of the list
	}
}

/**************************************/

//! Get optimal log base-2 overlap and window scalings for transients
//! The idea is that overlap occurs roughly at the center of a block,
//! and if we have a transient placed right smack in the middle, we
//! can simply reduce the overlap to control its pre-echo. This is
//! combined with window switching so that the transient lies mostly
//! in the center of a subblock, at which point overlap scaling does
//! its job to reduce the pre-echo.
//! The original transient detection idea was based off this paper:
//!  "Codierung von Audiosignalen mit uberlappender Transformation und adaptiven Fensterfunktionen"
//!  (Coding of Audio Signals with Overlapping Block Transform and Adaptive Window Functions)
//!  DOI: 10.1515/FREQ.1989.43.9.252
//! However, it has been modified so extensively that the only
//! takeaway from the paper was the idea of dividing the
//! averaged high-passed energy as a measure of transientness.
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
	float *StepBuffer,
	float *LastBlockEnergy,
	float *LastBlockSample,
	int BlockSize,
	int MinOverlap,
	int MaxOverlap,
	int nChan
) {
	int n, Chan;

	//! Perform a highpass filter on this block to get the step/transient energies
	//! PONDER: Use a proper halfband filter with cutoff at Fs/4, similar to LAME?
	for(n=0;n<BlockSize;n++) StepBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
		float SampleLast = LastBlockSample[Chan];
		const float *Src = Data + Chan*BlockSize;
		for(n=0;n<BlockSize;n++) {
			float v = Src[n];
			StepBuffer[n] += ABS(v - SampleLast);
			SampleLast = v;
		}
		LastBlockSample[Chan] = SampleLast;
	}

	//! Begin binary search for transient region until it centers within a subblock
	int   Decimation = 1, SubBlockSize = BlockSize;
	float StepL, StepLW;
	float StepM, StepMW;
	float StepR, StepRW;
	float StepOld = *LastBlockEnergy, StepOldW = 1.0f;
	for(;;) {
		//! Analyze the step energy for this subblock to find where
		//! the peak transient energy lies (left, middle, right).
		//! This is a weighted average, similar to a smooth max function
		StepM = 0.0f, StepMW = 0.0f;
		StepL = 0.0f, StepLW = 0.0f;
		StepR = 0.0f, StepRW = 0.0f; {
			const float *SrcL = StepBuffer;
			const float *SrcR = StepBuffer + SubBlockSize;
			const float *Cos = Fourier_SinTableN(SubBlockSize/2) + SubBlockSize/2;
			const float *Sin = Fourier_SinTableN(SubBlockSize/4);
			for(n=0;n<SubBlockSize/4;n++) { //! 0..Pi/4 (L/R)
				float c  = *--Cos;
				float dL = *SrcL++;
				float dR = *--SrcR;
				StepL += c*SQR(dL), StepLW += c*dL;
				StepR += c*SQR(dR), StepRW += c*dR;
			}
			for(n=0;n<SubBlockSize/4;n++) { //! Pi/4..Pi/2 (L/R + Mid)
				float s  = *Sin++;
				float c  = *--Cos;
				float dL = *SrcL++;
				float dR = *--SrcR;
				StepL += c*SQR(dL), StepLW += c*dL;
				StepR += c*SQR(dR), StepRW += c*dR;
				StepM += s*(SQR(dL) + SQR(dR)), StepMW += s*(dL + dR);
			}
		}

		//! Can we use window switching at all?
		//! NOTE: Limited to a minimum subblock size of 128 samples, and
		//! maximum decimation of 1/8 (Decimation > 8h: 1yyy = Decimation by 1/8)
		//! NOTE: Avoid division (and divide-by-0 errors) by cross-multiplying.
		if(ULC_USE_WINDOW_SWITCHING && SubBlockSize > 128 && Decimation < 0x8) {
			//! Check to see if the peak is at the center yet (MaxRatio(StepM/Step[L/R]) >= 0.5)
			float tL = StepM*StepLW, tLW = StepL*StepMW;
			float tR = StepM*StepRW, tRW = StepR*StepMW;
			if(tL < 0.5f*tLW || tLW < 0.5f*tL || tR < 0.5f*tRW || tRW < 0.5f*tR) {
				//! Decimate further and update the decimation pattern
				SubBlockSize /= 2;
				if(StepL*StepRW > StepR*StepLW) Decimation  = (Decimation<<1) | 0;
				else                            Decimation  = (Decimation<<1) | 1,
				                                StepBuffer += SubBlockSize,
				                                StepOld     = StepL,
				                                StepOldW    = StepLW;
				continue;
			}
		}

		//! If we reach here, then the transient either lies in between the left/right
		//! boundaries, or exists on both sides equally. Either way, we can't improve
		//! it any further, so we stop here and allow overlap switching to take over.
		break;
	}

	//! Use the jump to the middle from the last [sub]block to judge the overlap
	//! NOTE: As above, minimize division here by cross-multiplying instead.
	int OverlapScale;
	float Ra = StepM  * StepOldW;
	float Rb = StepMW * StepOld;
	if(Ra*0x1.6A09E6p-1f >= Rb) { //! a/b==Sqrt[2] is the first ratio to result in a ratio > 0
		//! a/b >= 2^6.5 gives the upper limit ratio of 7.0
		if(Ra*0x1.6A09E6p-7f < Rb) {
			//! 0x1.715476p0 = 1/Log[2], to get the log base-2
			float r = Ra / Rb;
			OverlapScale = (int)(0x1.715476p0f*logf(r) + 0.5f);
		} else OverlapScale = 7;
		while(OverlapScale > 0 && (SubBlockSize >> OverlapScale) < MinOverlap) OverlapScale--;
		while(OverlapScale < 7 && (SubBlockSize >> OverlapScale) > MaxOverlap) OverlapScale++;
	} else OverlapScale = 0;

	//! Return the combined overlap+window switching parameters
	*LastBlockEnergy = (StepRW > 0.0f) ? (StepR / StepRW) : 0.0f;
	return OverlapScale + 0x8*(SubBlockSize != BlockSize) + 0x10*Decimation;
}

/**************************************/

//! Transform a block and prepare its coefficients
static float *Block_Transform_SortComparator_KeysArray;
static int Block_Transform_SortComparator(const void *a, const void *b) {
	const float *BufferIdx = Block_Transform_SortComparator_KeysArray;
	return (BufferIdx[*(int*)a] < BufferIdx[*(int*)b]) ? (+1) : (-1);
}
static inline void Block_Transform_ScaleAndToNepers(float *Coef, float *CoefNp, int N) {
	int n;
	for(n=0;n<N;n++) {
		float v = Coef[n] * 2.0f/N;
		Coef  [n] = v;
		CoefNp[n] = (ABS(v) < 0.5f*ULC_COEF_EPS) ? ULC_COEF_NEPER_OUT_OF_RANGE : logf(ABS(v));
	}
}
static inline void Block_Transform_ComputePowerSpectrum(float *Power, float *PowerNp, const float *Re, const float *Im, int N) {
	int i;
	for(i=0;i<N;i++) {
		//! `v` should technically be scaled by (2/N)^2, but this power
		//! spectrum is only used for psychoacoustic analysis on this
		//! buffer /only/, and so scaling really doesn't matter here
		float v = SQR(Re[i]) + SQR(Im[i]);
		Power  [i] = v;
		PowerNp[i] = (ABS(v) < 0x1.0p-126f) ? ULC_COEF_NEPER_OUT_OF_RANGE : logf(ABS(v)); //! Do not allow subnormals
	}
}
static void Block_Transform_BufferInterleave(float *Buf, float *Tmp, int BlockSize, int Decimation) {
	//! The interleaving patterns here were chosen to try and
	//! optimize coefficient clustering across block sizes
	int n; for(n=0;n<BlockSize;n++) Tmp[n] = Buf[n];
	int n2 = BlockSize/2;
	int n4 = BlockSize/4;
	int n8 = BlockSize/8;
	switch(Decimation >> 1) { //! Lowermost bit only controls which subblock gets overlap scaling, so ignore it
		//! 001x: a=N/2, b=N/2
		case 0b001: {
			for(n=0;n<n2;n++) {
				*Buf++ = Tmp[  +n]; //! a
				*Buf++ = Tmp[n2+n]; //! b
			}
		} break;

		//! 010x: a=N/4, b=N/4, c=N/2
		case 0b010: {
			for(n=0;n<n4;n++) {
				*Buf++ = Tmp[  +n*1+0]; //! a
				*Buf++ = Tmp[n4+n*1+0]; //! b
				*Buf++ = Tmp[n2+n*2+0]; //! c0
				*Buf++ = Tmp[n2+n*2+1]; //! c1
			}
		} break;

		//! 011x: a=N/2, b=N/4, c=N/4
		case 0b011: {
			for(n=0;n<n4;n++) {
				*Buf++ = Tmp[     +n*2+0]; //! a0
				*Buf++ = Tmp[     +n*2+1]; //! a1
				*Buf++ = Tmp[n2   +n*1+0]; //! b
				*Buf++ = Tmp[n2+n4+n*1+0]; //! c
			}
		} break;

		//! 100x: a=N/8, b=N/8, c=N/4, d=N/2
		case 0b100: {
			for(n=0;n<n8;n++) {
				*Buf++ = Tmp[  +n*1+0]; //! a
				*Buf++ = Tmp[n8+n*1+0]; //! b
				*Buf++ = Tmp[n4+n*2+0]; //! c0
				*Buf++ = Tmp[n4+n*2+1]; //! c1
				*Buf++ = Tmp[n2+n*4+0]; //! d0
				*Buf++ = Tmp[n2+n*4+1]; //! d1
				*Buf++ = Tmp[n2+n*4+2]; //! d2
				*Buf++ = Tmp[n2+n*4+3]; //! d3
			}
		} break;

		//! 101x: a=N/4, b=N/8, c=N/8, d=N/2
		case 0b101: {
			for(n=0;n<n8;n++) {
				*Buf++ = Tmp[     +n*2+0]; //! a0
				*Buf++ = Tmp[     +n*2+1]; //! a1
				*Buf++ = Tmp[n4   +n*1+0]; //! b
				*Buf++ = Tmp[n4+n8+n*1+0]; //! c
				*Buf++ = Tmp[n2   +n*4+0]; //! d0
				*Buf++ = Tmp[n2   +n*4+1]; //! d1
				*Buf++ = Tmp[n2   +n*4+2]; //! d2
				*Buf++ = Tmp[n2   +n*4+3]; //! d3
			}
		} break;

		//! 110x: a=N/2, b=N/8, c=N/8, d=N/4
		case 0b110: {
			for(n=0;n<n8;n++) {
				*Buf++ = Tmp[     +n*4+0]; //! a0
				*Buf++ = Tmp[     +n*4+1]; //! a1
				*Buf++ = Tmp[     +n*4+2]; //! a2
				*Buf++ = Tmp[     +n*4+3]; //! a3
				*Buf++ = Tmp[n2   +n*1+0]; //! b
				*Buf++ = Tmp[n2+n8+n*1+0]; //! c
				*Buf++ = Tmp[n2+n4+n*2+0]; //! d0
				*Buf++ = Tmp[n2+n4+n*2+1]; //! d1
			}
		} break;

		//! 111x: a=N/2, b=N/4, c=N/8, d=N/8
		case 0b111: {
			for(n=0;n<n8;n++) {
				*Buf++ = Tmp[        +n*4+0]; //! a0
				*Buf++ = Tmp[        +n*4+1]; //! a1
				*Buf++ = Tmp[        +n*4+2]; //! a2
				*Buf++ = Tmp[        +n*4+3]; //! a3
				*Buf++ = Tmp[n2      +n*2+0]; //! b0
				*Buf++ = Tmp[n2      +n*2+1]; //! b1
				*Buf++ = Tmp[n2+n4   +n*1+0]; //! c
				*Buf++ = Tmp[n2+n4+n8+n*1+0]; //! d
			}
		} break;
	}
}
static int Block_Transform(struct ULC_EncoderState_t *State, const float *Data, float PowerDecay) {
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;

	//! Get the window control parameters for this block and the next
	int WindowCtrl     = State->WindowCtrl     = State->NextWindowCtrl;
	int NextWindowCtrl = State->NextWindowCtrl = Block_Transform_GetWindowCtrl(
		Data,
		State->TransformTemp,
		&State->LastBlockEnergy,
		State->LastBlockSample,
		BlockSize,
		State->MinOverlap,
		State->MaxOverlap,
		nChan
	);
	int NextSubBlockOverlap = BlockSize; {
		//! Adjust for the first [sub]block's size
		int NextDecimation = NextWindowCtrl >> 4;
		NextSubBlockOverlap >>= ULC_Helper_SubBlockDecimationPattern(NextDecimation)[0];

		//! Adjust for 'first [sub]block is transient'
		if(ULC_Helper_TransientSubBlockIndex(NextDecimation) == 0) NextSubBlockOverlap >>= (NextWindowCtrl & 0x7);
	}

	//! Transform channels and insert keys for each codeable coefficient
	//! It's not /strictly/ required to calculate nNzCoef, but it can
	//! speed things up in the rate-control step
	int nNzCoef = 0; {
		int n, Chan;
		float  AnalysisPowerNp = 0.0f; PowerDecay = logf(PowerDecay);
		float *BufferSamples   = State->SampleBuffer;
		float *BufferTransform = State->TransformBuffer;
		float *BufferNepers    = State->TransformNepers;
		float *BufferIndex     = (float*)State->TransformIndex;
		float *BufferFwdLap    = State->TransformFwdLap;
		float *BufferTemp      = State->TransformTemp;

		//! Transform the input data
		int ThisTransientSubBlockIdx = ULC_Helper_TransientSubBlockIndex(WindowCtrl >> 4);
		const int8_t *DecimationPattern = ULC_Helper_SubBlockDecimationPattern(WindowCtrl >> 4);
		for(Chan=0;Chan<nChan;Chan++) {
			//! Process each subblock sequentially
			float *BufferTransformEnd = BufferTransform + BlockSize;
			int SubBlockIdx;
			for(SubBlockIdx=0;(BufferTransform < BufferTransformEnd);SubBlockIdx++) {
				//! Get the size of this subblock and its overlap
				int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];
				int OverlapSize = BlockSize >> DecimationPattern[SubBlockIdx];
				if(SubBlockIdx == ThisTransientSubBlockIdx) OverlapSize >>= (WindowCtrl & 0x7);

				//! If the next block's first subblock is smaller than
				//! this one, set overlap to the smaller of the two.
				//! NOTE: This is literally THE ONLY reason that there is a
				//! coding delay in this codec (and most MDCT codecs).
				if(BufferTransform+SubBlockSize == BufferTransformEnd) {
					if(NextSubBlockOverlap < OverlapSize) OverlapSize = NextSubBlockOverlap;
				}

				//! If we have a long block, all the data is prepared and
				//! we can move straight on to the transform
				const float *SmpBuf = BufferSamples;
				if(SubBlockSize != BlockSize) {
					//! Use a scratch buffer for sample data
					SmpBuf = BufferTemp;

					//! With decimation, we cycle data through the lapping
					//! buffer, so that it will be ready for future calls
					int LapExtra = (BlockSize - SubBlockSize) / 2;
					float *LapBuf = BufferFwdLap + LapExtra;

					//! Get samples from the lapping buffer
					int LapCopy = (LapExtra < SubBlockSize) ? LapExtra : SubBlockSize;
					for(n=0;n<LapCopy;n++)   BufferTemp[n] = LapBuf[-1-n];
					for(;n<SubBlockSize;n++) BufferTemp[n] = *BufferSamples++;

					//! Shuffle new samples into it
					int NewCopy = LapExtra - LapCopy;
					for(n=0;n<NewCopy;n++) LapBuf[-1-n] = LapBuf[-1-n-LapCopy];
					for(;n<LapExtra;n++) LapBuf[-1-n] = *BufferSamples++;
				} else BufferSamples += SubBlockSize;

				//! Perform the actual analyses
				float *BufferEnergy   = BufferTemp;
				float *BufferEnergyNp = BufferTemp + SubBlockSize;
				Fourier_MDCT(
					BufferTransform,
					SmpBuf,
					BufferFwdLap + (BlockSize-SubBlockSize)/2,
					BufferTemp,
					SubBlockSize,
					OverlapSize,
					BufferNepers //! MDST coefficients stored here temporarily
				);
				Block_Transform_ComputePowerSpectrum(BufferEnergy, BufferEnergyNp, BufferTransform, BufferNepers, SubBlockSize);
				Block_Transform_ScaleAndToNepers(BufferTransform, BufferNepers, SubBlockSize);
				Block_Transform_WriteSortValues(
					BufferIndex,
					BufferEnergy,
					BufferEnergyNp,
					BufferNepers,
					&nNzCoef,
					SubBlockSize,
					AnalysisPowerNp,
					State->RateHz * 0.5f
				);

				//! Move to the next subblock
				BufferTransform += SubBlockSize;
				BufferNepers    += SubBlockSize;
				BufferIndex     += SubBlockSize;
			}

			//! Cache the sample data for the next block
			for(n=0;n<BlockSize;n++) BufferSamples[n-BlockSize] = *Data++;

			//! Move to the next channel
			BufferFwdLap    += BlockSize/2u;
			AnalysisPowerNp += PowerDecay;
		}

		//! Interleave the transform data for coding
		if(WindowCtrl & 0x8) {
			int Decimation = WindowCtrl >> 4;
			BufferTransform = State->TransformBuffer;
			BufferNepers    = State->TransformNepers;
			BufferIndex     = (float*)State->TransformIndex;
			for(Chan=0;Chan<nChan;Chan++) {
				Block_Transform_BufferInterleave(BufferTransform, BufferTemp, BlockSize, Decimation);
				Block_Transform_BufferInterleave(BufferNepers,    BufferTemp, BlockSize, Decimation);
				Block_Transform_BufferInterleave(BufferIndex,     BufferTemp, BlockSize, Decimation);
				BufferTransform += BlockSize;
				BufferNepers    += BlockSize;
				BufferIndex     += BlockSize;
			}
		}
	}

	//! Create the coefficient sorting indices
	//! TODO: Be more efficient; this is suboptimal, but
	//! fixing would require implementing a customized
	//! sorting routine
	{
		int Idx, nIdx = nChan * BlockSize;
		int *BufferTmp = (int*)State->TransformTemp;
		int *BufferIdx = State->TransformIndex;
		for(Idx=0;Idx<nIdx;Idx++) BufferTmp[Idx] = Idx;
		Block_Transform_SortComparator_KeysArray = (float*)State->TransformIndex;
		qsort(BufferTmp, nIdx, sizeof(int), Block_Transform_SortComparator);
		for(Idx=0;Idx<nIdx;Idx++) BufferIdx[BufferTmp[Idx]] = Idx;
	}
	return nNzCoef;
}

/**************************************/
//! EOF
/**************************************/
