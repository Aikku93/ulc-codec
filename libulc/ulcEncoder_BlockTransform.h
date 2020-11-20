/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
#include <stdint.h>
/**************************************/
#include "Fourier.h"
#include "ulcEncoder_Psycho.h"
#include "ulcEncoder_WindowControl.h"
#include "ulcHelper.h"
/**************************************/

//! Write the sort values for all coefficients in a block
//! and return the number of codeable non-zero coefficients
//! NOTE:
//!  AnalysisPower is used to alter the preference for the currently-being-analyzed channel
static inline void Block_Transform_WriteSortValues(
	float *CoefIdx,
	const uint32_t *Energy,
	const uint32_t *EnergyNp,
	const float *CoefNp,
	int   *nNzCoef,
	int   BlockSize,
	float AnalysisPowerNp,
	float NyquistHz
) {
	int Band;
#if ULC_USE_PSYCHOACOUSTICS
	struct Block_Transform_MaskingState_t MaskingState;
	Block_Transform_MaskingState_Init(&MaskingState, Energy, EnergyNp, BlockSize, NyquistHz);
#else
	(void)Energy;
	(void)EnergyNp;
	(void)NyquistHz;
	AnalysisPowerNp *= 0.5f; //! Assumed to operate in the energy (X^2) domain, so account for operating in a linear domain
#endif
	for(Band=0;Band<BlockSize;Band++) {
		//! Inside the codeable range?
		float ValNp = CoefNp[Band];
		if(ValNp != ULC_COEF_NEPER_OUT_OF_RANGE) {
#if ULC_USE_PSYCHOACOUSTICS
			//! Apply psychoacoustic corrections to this band energy
			//! NOTE: Operate in the energy (X^2) domain, as this
			//! seems to result in slightly improved quality.
			ValNp += ValNp + Block_Transform_GetMaskedLevel(&MaskingState, Energy, EnergyNp, Band, BlockSize);
#endif
			//! Store the sort value for this coefficient
			CoefIdx[Band] = ValNp + AnalysisPowerNp;
			(*nNzCoef)++;
		} else CoefIdx[Band] = -0x1.0p126f; //! Unusable coefficient; map to the end of the list
	}
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
static inline void Block_Transform_ComputePowerSpectrum(uint32_t *Power, uint32_t *PowerNp, const float *Re, const float *Im, int N) {
	int i;

	//! Find the maximum amplitude in the analysis
	//! This is used to scale the buffer and avoid precision issues
	double Scale = 0.0; {
		double Max = 0.0;
		for(i=0;i<N;i++) {
			double v = SQR((double)Re[i]) + SQR((double)Im[i]);
			if(v > Max) Max = v;
		}
		if(Max) Scale = 0x1.0p32 / Max;
	}

	//! Rescale and convert data to fixed-point
	//! NOTE: Maximum value for PowerNp is Log[2 * 2^32] == 22.87,
	//! allowing us to scale this by (2^32 - 1)/Log[2 * 2^32]
	const double LogScale = 0x1.66235B759A0B1p27; //! (2^32 - 1)/Log[2 * 2^32]
	const double CeilBias = 0x1.FFFFFFFFFFFFFp-1; //! 0.99999... for ceiling
	for(i=0;i<N;i++) {
		double v  = SQR((double)Re[i]) + SQR((double)Im[i]);
		       v *= Scale;
		Power  [i] = (uint32_t)(v + CeilBias);
		PowerNp[i] = (v < 0.5) ? 0 : (uint32_t)(LogScale*log(2.0*v) + CeilBias); //! Scale by 2 to keep values strictly non-negative
	}
}
static void Block_Transform_BufferInterleave(float *Buf, float *Tmp, int BlockSize, int Decimation) {
	//! The interleaving patterns here were chosen to try and
	//! optimize coefficient clustering across block sizes
	int n; for(n=0;n<BlockSize;n++) Tmp[n] = Buf[n];
	switch(Decimation >> 1) { //! Lowermost bit only controls which subblock gets overlap scaling, so ignore it
		//! 001x: a=N/2, b=N/2
		case 0b001: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/2;
			for(n=0;n<BlockSize/2;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
			}
		} break;

		//! 010x: a=N/4, b=N/4, c=N/2
		case 0b010: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/4;
			float *SrcC = SrcB + BlockSize/4;
			for(n=0;n<BlockSize/4;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcC++;
			}
		} break;

		//! 011x: a=N/2, b=N/4, c=N/4
		case 0b011: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/2;
			float *SrcC = SrcB + BlockSize/4;
			for(n=0;n<BlockSize/4;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
			}
		} break;

		//! 100x: a=N/8, b=N/8, c=N/4, d=N/2
		case 0b100: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/8;
			float *SrcC = SrcB + BlockSize/8;
			float *SrcD = SrcC + BlockSize/4;
			for(n=0;n<BlockSize/8;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
			}
		} break;

		//! 101x: a=N/4, b=N/8, c=N/8, d=N/2
		case 0b101: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/4;
			float *SrcC = SrcB + BlockSize/8;
			float *SrcD = SrcC + BlockSize/8;
			for(n=0;n<BlockSize/8;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
			}
		} break;

		//! 110x: a=N/2, b=N/8, c=N/8, d=N/4
		case 0b110: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/2;
			float *SrcC = SrcB + BlockSize/8;
			float *SrcD = SrcC + BlockSize/8;
			for(n=0;n<BlockSize/8;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
			}
		} break;

		//! 111x: a=N/2, b=N/4, c=N/8, d=N/8
		case 0b111: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/2;
			float *SrcC = SrcB + BlockSize/4;
			float *SrcD = SrcC + BlockSize/8;
			for(n=0;n<BlockSize/8;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcD++;
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
		State->SampleBuffer,
		State->TransformTemp,
		BlockSize,
		State->MinOverlap,
		State->MaxOverlap,
		nChan,
		State->RateHz
	);
	int NextBlockSubBlockSize = BlockSize; {
		//! Adjust for the first [sub]block's size
		int NextDecimation = NextWindowCtrl >> 4;
		NextBlockSubBlockSize >>= ULC_Helper_SubBlockDecimationPattern(NextDecimation)[0];
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

		//! Apply M/S transform to the data
		//! NOTE: Fully normalized; not orthogonal.
		if(nChan == 2) for(n=0;n<BlockSize;n++) {
			float L = BufferSamples[n];
			float R = BufferSamples[n + BlockSize];
			BufferSamples[n]             = (L+R) * 0.5f;
			BufferSamples[n + BlockSize] = (L-R) * 0.5f;
		}

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
				if(BufferTransform+SubBlockSize < BufferTransformEnd) {
					int NextSubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx+1];
					if(OverlapSize > NextSubBlockSize) OverlapSize = NextSubBlockSize;
				} else {
					if(OverlapSize > NextBlockSubBlockSize) OverlapSize = NextBlockSubBlockSize;
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
				uint32_t *BufferEnergy   = (uint32_t*)BufferTemp;
				uint32_t *BufferEnergyNp = (uint32_t*)BufferTemp + SubBlockSize;
				Fourier_MDCT(
					BufferTransform,
					SmpBuf,
					BufferFwdLap + (BlockSize-SubBlockSize)/2,
					BufferTemp,
					SubBlockSize,
					OverlapSize,
#if ULC_USE_PSYCHOACOUSTICS
					BufferNepers //! MDST coefficients stored here temporarily
#else
					NULL
#endif
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
