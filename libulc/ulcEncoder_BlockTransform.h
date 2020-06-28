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
#include "ulcEncoder_Helper.h"
#include "ulcEncoder_Psycho.h"
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
			//! Convert to the energy domain
			ValNp *= 2.0f;
#if ULC_USE_PSYCHOACOUSTICS
			//! Apply psychoacoustic corrections to this band energy
			ValNp += Block_Transform_GetMaskedLevel(&MaskingState, Energy, EnergyNp, Band, BlockSize);
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

//! Get optimal log base-2 overlap scaling for transients
//! The idea is that with reduced overlap, some transients
//! (especially 'pops') come for 'free' and additionally
//! reduce the amount of spectral distortion that may affect
//! psychoacoustics calculations.
//! The original transient detection idea was based off this paper:
//!  "Codierung von Audiosignalen mit uberlappender Transformation und adaptiven Fensterfunktionen"
//!  (Coding of Audio Signals with Overlapping Block Transform and Adaptive Window Functions)
//!  DOI: 10.1515/FREQ.1989.43.9.252
//! However, it has been modified so extensively that the only
//! takeaway from the paper was the idea of dividing the
//! averaged high-passed energy as a measure of transientness.
//! NOTE: When window switching is used, the function will return
//! a negative value corresponding to the transform size decimation.
static inline int Block_Transform_GetLogOverlapScale(
	const float *Data,
	float *LastBlockEnergy,
	float *LastBlockSample,
	int BlockSize,
	int MinOverlap,
	int MaxOverlap,
	int nChan
) {
	int n, Chan;
	int OverlapScale = 0; //! Default to long block

	//! Get step energy for this block (first and second half)
	float aIn = 0.0f, aOut = 0.0f;
	float bIn = 0.0f, bOut = 0.0f;
	const float *Sin = Fourier_SinTableN(BlockSize/2);
	const float *Cos = Sin + BlockSize/2; //! Accessed backwards
	for(Chan=0;Chan<nChan;Chan++) {
		const float *Src = Data + Chan*BlockSize;
		float v, d, SampleLast = LastBlockSample[Chan];
		for(n=0;n<BlockSize/2;n++) {
			v = Src[n];
			d = v - SampleLast;
			aIn  += SQR(d*Sin[n]);
			aOut += SQR(d*Cos[-1-n]);
			SampleLast = v;
		}
		for(n=0;n<BlockSize/2;n++) {
			v = Src[BlockSize/2+n];
			d = v - SampleLast;
			bIn  += SQR(d*Sin[n]);
			bOut += SQR(d*Cos[-1-n]);
			SampleLast = v;
		}
		LastBlockSample[Chan] = SampleLast;
	}

	//! Overlap switching
	float OverlapSwitchRatio = 1.0f;
	if(OverlapScale == 0) {
		//! Relate the average step energy in this block to that of the last block
		//! EnergyLapped contains the contributions of the second half of
		//! the last block PLUS the first half of this block, with a peak
		//! at the midpoint. EnergyCenter contains the contributions of
		//! both halves of this block, with a peak at the midpoint. This
		//! allows better detection of transients whilst preserving (as
		//! much as we can) smoothness in the pre-transient section; any
		//! block distortion before the transient tends to result in very
		//! audible artifacts, so these must be minimized.
		float Ra = 2.0f*aIn + bOut;
		float Rb = *LastBlockEnergy + 2.0f*aOut;
		if(Ra*0x1.6A09E6p-1f >= Rb) { //! Ra/Rb==Sqrt[2] is the first ratio to result in a ratio > 0
			if(Ra*0x1.6A09E6p-7f < Rb) { //! Ra/Rb >= 2^6.5 gives the upper limit ratio of 7.0
				//! 0x1.715476p0 = 1/Log[2], to get the log base-2
				//! NOTE: It works better to use the squared energy ratio here.
				OverlapSwitchRatio = Ra / Rb;
				OverlapScale = (int)(0x1.715476p0f*logf(OverlapSwitchRatio) + 0.5f);
			} else OverlapScale = 0xF;
			while((BlockSize >> OverlapScale) < MinOverlap) OverlapScale--;
			while((BlockSize >> OverlapScale) > MaxOverlap) OverlapScale++;
		}
	}

	//! Window switching (limited to a minimum of 128-sample subblocks)
	if(ULC_USE_WINDOW_SWITCHING && BlockSize > 128) {
#if 0 //! The window switching criteria is too poor to be usable at present
		//! Relate the average step energy in this block to that of the last block,
		//! but unlike the below overlap switching, use the energy that is present
		//! in the whole block instead of just the transition region
		float Ra = aIn + 2.0f*aOut;
		float Rb = (2.0f*(*LastBlockEnergy) + aIn) * OverlapSwitchRatio;
		if(Ra*0x1.6A09E6p-1f >= Rb) { //! Ra/Rb==Sqrt[2] is the first ratio to result in a ratio > 0
			int Scale;
			if(Ra*0x1.6A09E6p-7f < Rb) { //! Ra/Rb >= 2^6.5 gives the upper limit ratio of 7.0
				//! 0x1.715476p0 = 1/Log[2], to get the log base-2
				//! NOTE: As above, using the squared ratio works better.
				float r = Ra / Rb;
				Scale = (int)(0x1.715476p0f*logf(r) + 0.5f);
			} else Scale = 7;
			while((BlockSize >> Scale) < 128) Scale--;
			OverlapScale = -Scale;
		}
#endif
	}

	//! Save state for next block and return overlap/window switching
	*LastBlockEnergy = bIn;
	return OverlapScale;
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
static void Block_Transform_BufferInterleave(float *Buf, float *Tmp, int BlockSize, int nSubBlocks) {
	int n, SubBlock;
	int SubBlockSize = BlockSize / nSubBlocks;
	for(n=0;n<BlockSize;n++) Tmp[n] = Buf[n];
	for(SubBlock=0;SubBlock<nSubBlocks;SubBlock++) {
		for(n=0;n<SubBlockSize;n++) Buf[n*nSubBlocks+SubBlock] = Tmp[SubBlock*SubBlockSize+n];
	}
}
static int Block_Transform(struct ULC_EncoderState_t *State, const float *Data, float PowerDecay) {
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;

	//! Get the overlap scaling for the next block's transition
	int ThisOverlapScale = State->ThisOverlap = State->NextOverlap;
	int NextOverlapScale = State->NextOverlap = Block_Transform_GetLogOverlapScale(
		Data,
		&State->LastBlockEnergy,
		State->LastBlockSample,
		BlockSize,
		State->MinOverlap,
		State->MaxOverlap,
		nChan
	);
#if 0 //! Testing
	if(ThisOverlapScale < NextOverlapScale) NextOverlapScale = State->NextOverlap = ThisOverlapScale;
#endif
	//! Transform channels and insert keys for each codeable coefficient
	//! It's not /strictly/ required to calculate nNzCoef, but it can
	//! speed things up in the rate-control step
	int nNzCoef = 0; {
		int n, Chan, SubBlock;
		float  AnalysisPowerNp = 0.0f; PowerDecay = logf(PowerDecay);
		float *BufferTransform = State->TransformBuffer;
		float *BufferNepers    = State->TransformNepers;
		float *BufferIndex     = (float*)State->TransformIndex;
		float *BufferFwdLap    = State->TransformFwdLap;
		float *BufferTemp      = State->TransformTemp;
		float *BufferEnergy    = BufferTemp;
		float *BufferEnergyNp  = BufferTemp + BlockSize;
		int    OverlapSize     = BlockSize;
		int    SubBlockSize    = BlockSize;
		int    nSubBlocks      = 1;

		//! Adjust for window switching (or overlap switching)
		if(ThisOverlapScale < 0) {
			SubBlockSize >>= -ThisOverlapScale;
			nSubBlocks   <<= -ThisOverlapScale;
			OverlapSize    =  SubBlockSize;
		} else OverlapSize >>= ABS(NextOverlapScale);

		//! Transform the input data
		for(Chan=0;Chan<nChan;Chan++) {
			//! Long block?
			if(nSubBlocks == 1) {
				Fourier_MDCT(
					BufferTransform,
					Data,
					BufferFwdLap,
					BufferTemp,
					SubBlockSize,
					OverlapSize,
					BufferNepers //! MDST coefficients stored here temporarily
				);
			} else {
				//! Transform half the subblocks from the lapping buffer
				//! and the other half from the first half of this block
				float *DataTemp = BufferTemp + SubBlockSize;
				float *LapBuf = BufferFwdLap + (BlockSize-SubBlockSize)/2;
				for(SubBlock=0;SubBlock<nSubBlocks/2-1;SubBlock++) {
					for(n=0;n<SubBlockSize;n++) DataTemp[n] = *--LapBuf;
					Fourier_MDCT(
						BufferTransform,
						DataTemp,
						BufferFwdLap + (BlockSize-SubBlockSize)/2,
						BufferTemp,
						SubBlockSize,
						OverlapSize,
						BufferNepers //! MDST coefficients stored here temporarily
					);
					BufferTransform += SubBlockSize;
					BufferNepers    += SubBlockSize;
				}
				{
					for(n=0;n<SubBlockSize/2;n++) {
						DataTemp[n] = *--LapBuf;
						DataTemp[SubBlockSize/2+n] = Data[n];
					}
					Fourier_MDCT(
						BufferTransform,
						DataTemp,
						BufferFwdLap + (BlockSize-SubBlockSize)/2,
						BufferTemp,
						SubBlockSize,
						OverlapSize,
						BufferNepers //! MDST coefficients stored here temporarily
					);
					BufferTransform += SubBlockSize;
					BufferNepers    += SubBlockSize;
				}
				for(SubBlock++;SubBlock<nSubBlocks;SubBlock++) {
					Fourier_MDCT(
						BufferTransform,
						Data + SubBlock*SubBlockSize + SubBlockSize/2 - BlockSize/2,
						BufferFwdLap + (BlockSize-SubBlockSize)/2,
						BufferTemp,
						SubBlockSize,
						OverlapSize,
						BufferNepers //! MDST coefficients stored here temporarily
					);
					BufferTransform += SubBlockSize;
					BufferNepers    += SubBlockSize;
				}
				BufferTransform -= BlockSize;
				BufferNepers    -= BlockSize;

				//! Refill the lapping buffer with samples from this block
				for(n=0;n<(BlockSize-SubBlockSize)/2;n++) BufferFwdLap[n] = Data[BlockSize-1-n];
			}

			//! Perform analysis of the subblocks
			for(SubBlock=0;SubBlock<nSubBlocks;SubBlock++) {
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

			//! Move to the next channel
			Data            += BlockSize;
			BufferFwdLap    += BlockSize/2u;
			AnalysisPowerNp += PowerDecay;
		}

		//! Interleave the transform data for coding
		BufferTransform = State->TransformBuffer;
		BufferNepers    = State->TransformNepers;
		BufferIndex     = (float*)State->TransformIndex;
		if(nSubBlocks > 1) for(Chan=0;Chan<nChan;Chan++) {
			Block_Transform_BufferInterleave(BufferTransform, BufferTemp, BlockSize, nSubBlocks);
			Block_Transform_BufferInterleave(BufferNepers,    BufferTemp, BlockSize, nSubBlocks);
			Block_Transform_BufferInterleave(BufferIndex,     BufferTemp, BlockSize, nSubBlocks);
			BufferTransform += BlockSize;
			BufferNepers    += BlockSize;
			BufferIndex     += BlockSize;
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
