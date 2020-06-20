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
//! The idea is that with reduced overlap, transients need fewer
//! coefficients to sound correct (at the cost of distortion)
//! Transient detection is loosely based on ideas found in:
//!  "Codierung von Audiosignalen mit uberlappender Transformation und adaptiven Fensterfunktionen"
//!  (Coding of Audio Signals with Overlapping Block Transform and Adaptive Window Functions)
//!  DOI: 10.1515/FREQ.1989.43.9.252
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

	//! Get the weighted step energy, modulating for fade-in
	//! and fade-out. As far as MDCT is concerned, this current
	//! block will be fading in, and so the important transient
	//! energy is towards the end of the block. The previous
	//! block is simultaneously fading out, so the transient
	//! energy to analyze against should use a reversed window.
	//! NOTE: There is no point in normalizing, as these energy
	//! values will be divided by the previous block's, and as
	//! they have the same normalization factor, it cancels out.
	float EnergyCenter = 0.0f;
	float EnergyEdge   = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
		float SampleLast = LastBlockSample[Chan];
		const float *Src  = Data + Chan*BlockSize;
		const float *WinS = Fourier_SinTableN(BlockSize);
		const float *WinC = WinS + BlockSize;
		for(n=0;n<BlockSize;n++) {
			//! NOTE: From testing, it appears that taking the difference
			//! between non-linearly-mapped samples works much better than
			//! using the samples directly.
			//! NOTE: It is important to remove the sign from the input
			//! values prior to doing anything with them. As an extreme
			//! example of what happens when this is not done: Consider a
			//! sine wave at Nyquist frequency; this results in a very
			//! large step energy accumulation that may cause selection of
			//! a very narrow overlap, despite audible distortion.
			float v = ABS(Src[n]); v *= SQR(v);
			float d = SQR(v - SampleLast);
			float s = WinS[n];
			float c = WinC[-1-n];
			EnergyCenter += s*d;
			EnergyEdge   += c*d;
			SampleLast    = v;
		}
		LastBlockSample[Chan] = SampleLast;
	}

	//! Relate the average step size of this block to that of the last block
	int OverlapScale; {
		float Ra = EnergyCenter;
		float Rb = *LastBlockEnergy;
		if(Ra >= Rb*2.0f) { //! Ra/Rb==2.0 is the first ratio to result in OverlapScale > 0
			if(Ra*0x1.0p-29f < Rb) { //! Ra/Rb >= 2^29 gives the upper limit of OverlapScale == 15
				//! 0x1.715476p0 = 1/Log[2], to get the log base-2
				float r = Ra / Rb;
				OverlapScale = (int)(0x1.715476p0f*0.5f*logf(r) + 0.5f);
			} else OverlapScale = 0xF;
		} else OverlapScale = 0;
		*LastBlockEnergy = EnergyEdge;
	}

	//! Clip to the minimum/maximum allowed overlap
	while((BlockSize >> OverlapScale) < MinOverlap) OverlapScale--;
	while((BlockSize >> OverlapScale) > MaxOverlap) OverlapScale++;
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
	float Scale2 = SQR(2.0f/N);
	for(i=0;i<N;i++) {
		float v = Scale2 * (SQR(Re[i]) + SQR(Im[i]));
		Power  [i] = v;
		PowerNp[i] = (ABS(v) < SQR(0.5f*ULC_COEF_EPS)) ? ULC_COEF_NEPER_OUT_OF_RANGE : logf(ABS(v));
	}
}
static int Block_Transform(struct ULC_EncoderState_t *State, const float *Data, float PowerDecay) {
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;

	//! Get the overlap scaling for this block
	int OverlapScale = State->ThisOverlap = Block_Transform_GetLogOverlapScale(
		Data,
		&State->LastBlockEnergy,
		State->LastBlockSample,
		BlockSize,
		State->MinOverlap,
		State->MaxOverlap,
		nChan
	);

	//! Transform channels and insert keys for each codeable coefficient
	//! It's not /strictly/ required to calculate nNzCoef, but it can
	//! speed things up in the rate-control step
	int nNzCoef = 0; {
		int Chan;
		float  AnalysisPowerNp = 0.0f; PowerDecay = logf(PowerDecay);
		float *BufferTransform = State->TransformBuffer;
		float *BufferNepers    = State->TransformNepers;
		float *BufferIndex     = (float*)State->TransformIndex;
		float *BufferFwdLap    = State->TransformFwdLap;
		float *BufferTemp      = State->TransformTemp;
		float *BufferEnergy    = BufferTemp;
		float *BufferEnergyNp  = BufferTemp + BlockSize;
		for(Chan=0;Chan<nChan;Chan++) {
			//! Do transform processing and write the sort values
			Fourier_MDCT(
				BufferTransform,
				Data,
				BufferFwdLap,
				BufferTemp,
				BlockSize,
				BlockSize >> OverlapScale,
				BufferNepers //! MDST coefficients stored here temporarily
			);
			Block_Transform_ComputePowerSpectrum(BufferEnergy, BufferEnergyNp, BufferTransform, BufferNepers, BlockSize);
			Block_Transform_ScaleAndToNepers(BufferTransform, BufferNepers, BlockSize);
			Block_Transform_WriteSortValues(
				BufferIndex,
				BufferEnergy,
				BufferEnergyNp,
				BufferNepers,
				&nNzCoef,
				BlockSize,
				AnalysisPowerNp,
				State->RateHz * 0.5f
			);

			//! Move to the next channel
			Data            += BlockSize;
			BufferTransform += BlockSize;
			BufferNepers    += BlockSize;
			BufferIndex     += BlockSize;
			BufferFwdLap    += BlockSize/2u;
			AnalysisPowerNp += PowerDecay;
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
