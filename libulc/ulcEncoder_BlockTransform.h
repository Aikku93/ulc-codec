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
	float *Coef,
	float *CoefNp,
	int   *nNzCoef,
	int   BlockSize,
	float AnalysisPowerNp,
	float NyquistHz
) {
	int Band;
#if ULC_USE_PSYCHOACOUSTICS
	struct Block_Transform_MaskingState_t MaskingState;
	Block_Transform_MaskingState_Init(&MaskingState, Coef, CoefNp, BlockSize, NyquistHz);
#else
	(void)NyquistHz;
#endif
	for(Band=0;Band<BlockSize;Band++) {
		//! Inside the codeable range?
		float ValNp = CoefNp[Band];
		if(ValNp != ULC_COEF_NEPER_OUT_OF_RANGE) {
#if ULC_USE_PSYCHOACOUSTICS
			//! Get the effective energy of this band by removing the background contribution
			//! NOTE: Not sure why this equation is the way it is. It was developed
			//! over 8 hours of experimentation to find what worked best in all cases
			float InvFlat, Mask = Block_Transform_UpdateMaskingThreshold(&MaskingState, Coef, CoefNp, Band, BlockSize, &InvFlat);
			ValNp += (expf(0x1.62E430p-1f*InvFlat))*(ValNp - Mask); //! 0x1.62E430p-1 = Log[2], to get 2^x
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
	float *EnergyBuffer,
	float *LastBlockEnergy,
	float *LastSampleEnergy,
	int BlockSize,
	int MinOverlap,
	int MaxOverlap,
	int nChan
) {
	int n, Chan;

	//! Combine all channel energy into a single buffer
	for(n=0;n<BlockSize;n++) EnergyBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
#if defined(__AVX__)
		for(n=0;n<BlockSize;n+=8) {
			__m256 v = _mm256_load_ps(Data + Chan*BlockSize + n);
			__m256 b = _mm256_load_ps(EnergyBuffer + n);
#if defined(__FMA__)
			b = _mm256_fmadd_ps(v, v, b);
#else
			v = _mm256_mul_ps(v, v);
			b = _mm256_add_ps(b, v);
#endif
			_mm256_store_ps(EnergyBuffer + n, b);
		}
#elif defined(__SSE__)
		for(n=0;n<BlockSize;n+=4) {
			__m128 v = _mm_load_ps(Data + Chan*BlockSize + n);
			__m128 b = _mm_load_ps(EnergyBuffer + n);
#if defined(__FMA__)
			b = _mm_fmadd_ps(v, v, b);
#else
			v = _mm_mul_ps(v, v);
			b = _mm_add_ps(b, v);
#endif
			_mm_store_ps(EnergyBuffer + n, b);
		}
#else
		for(n=0;n<BlockSize;n++) EnergyBuffer[n] += SQR(Data[Chan*BlockSize + n]);
#endif
	}

	//! Analyze sample energy to build a log-domain ratio
	int OverlapScale; {
		//! Get the squared step energy for each sample, and then
		//! apply a window to the step energy to get a weighted sum
		//! of squares. As far as MDCT is concerned, this current
		//! block will be fading in, and so the important transient
		//! energy is towards the end of the block. The previous
		//! block is simultaneously fading out, so the transient
		//! energy to analyze against should use a reversed window.
		//! NOTE: There is no point in normalizing, as these energy
		//! values will be divided by the previous block's, and as
		//! they have the same normalization factor, it cancels out.
		//! TODO: Try using different highpass filters; this filter
		//! (using the step energy) does not give very good bandpass
		//! selectivity, which may be affecting results; LAME uses a
		//! filter that rejects all frequencies below Fs/2 @ -50dB.
		float LastEnergy   = *LastSampleEnergy;
		float EnergyCenter = 1.0e-30f; //! Add a small bias to avoid division by zero
		float EnergyEdge   = 1.0e-30f;
		const float *WinS = Fourier_SinTableN(BlockSize);
		const float *WinC = WinS + BlockSize;
		for(n=0;n<BlockSize;n++) {
			float v = sqrtf(*EnergyBuffer++);
			float s = *WinS++;
			float c = *--WinC;
			float d = v - LastEnergy;
			LastEnergy = v;
			EnergyCenter += s*SQR(d);
			EnergyEdge   += c*SQR(d);
		}
		*LastSampleEnergy = LastEnergy;

		//! Relate the average step size of this block to that of the last block
		//! NOTE: If the energy decayed in this block, don't shrink this block's
		//! overlap; release transients are easily masked, and if it wasn't really
		//! a transient, then it was very likely a simple volume envelope.
		//! NOTE: Scaling by EnergyCenter/EnergyEdge gives better results, as this
		//! basically controls just how far we can narrow the overlap without any
		//! click/pop artifacts from discontinuities at the overlap boundaries.
		//! When a transient is extremely sharp/poppy, then we can decrease the
		//! overlap a lot further, as the transient itself will be masking any
		//! discontinuity artifacts from too narrow of an overlap.
		//! NOTE: Correctly map EnergyCenter/EnergyEdge via the square-root. Not
		//! doing this tends to overdo poppy transients in smooth regions.
		if(EnergyCenter > *LastBlockEnergy) {
			float LogRatio = sqrtf(EnergyCenter / EnergyEdge) * logf(EnergyCenter / *LastBlockEnergy);
			OverlapScale = (int)(0x1.715476p0f*LogRatio + 0.5f); //! 0x1.715476p0 = 1/Log[2], to get the log base-2
			if(OverlapScale > 0xF) OverlapScale = 0xF;
		} else OverlapScale = 0;
		*LastBlockEnergy = EnergyEdge;
	}

	//! Use the above-derived ratio to set an overlap amount for this block
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
static int Block_Transform(struct ULC_EncoderState_t *State, const float *Data, float PowerDecay) {
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;

	//! Get the overlap scaling for this block
	int OverlapScale = State->ThisOverlap = Block_Transform_GetLogOverlapScale(
		Data,
		State->TransformTemp,
		&State->LastBlockEnergy,
		&State->LastSampleEnergy,
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
		for(Chan=0;Chan<nChan;Chan++) {
			//! Apply transform and write the sort values
			Fourier_MDCT(
				BufferTransform,
				Data,
				BufferFwdLap,
				BufferTemp,
				BlockSize,
				BlockSize >> OverlapScale
			);
			Block_Transform_ScaleAndToNepers(BufferTransform, BufferNepers, BlockSize);
			Block_Transform_WriteSortValues(
				BufferIndex,
				BufferTransform,
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
