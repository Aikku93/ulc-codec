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
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_Helper.h"
#include "ulcEncoder_Psycho.h"
/**************************************/

//! Insert keys for block coefficients
//! Returns updated number of keys in list
//! NOTE:
//!  AnalysisPower is used to alter the preference for the currently-being-analyzed channel
static inline void Block_Transform_InsertKeys(
	struct AnalysisKey_t *Keys,
	const float *Coef,
	const float *CoefNp,
	int   BlockSize,
	int  *nKeys,
	int   Chan,
	float AnalysisPowerNp,
	float NyquistHz
) {
	int Band;
#if ULC_USE_PSYCHOACOUSTICS
	struct Block_Transform_MaskingState_t MaskingState;
	Block_Transform_MaskingState_Init(&MaskingState, Coef, CoefNp, BlockSize, NyquistHz);
#else
	(void)Coef;
	(void)NyquistHz;
#endif
	//! Analyze to create the quantizer zones
	for(Band=0;Band<BlockSize;Band++) {
		//! Check that the value is in range of the smallest quantization
		float ValNp = CoefNp[Band]; if(ValNp == ULC_COEF_NEPER_OUT_OF_RANGE) continue;
#if ULC_USE_PSYCHOACOUSTICS
		//! Get the effective energy of this band by removing the background contribution
		//! NOTE: Not sure why this masking equation is the way it is.
		//! From experiments, this seems to give the best tradeoff between
		//! brightness/clarity and tonal stability
		float Flat, Mask = Block_Transform_UpdateMaskingThreshold(&MaskingState, Coef, CoefNp, Band, BlockSize, &Flat);
		ValNp = (2.0f+Flat)*ValNp - (1.0f+Flat)*Mask;
#endif
		//! Insert key for this band
		//! NOTE: It is not necessary to re-map back to the linear
		//! domain, as this value is only used for sorting
		Keys[*nKeys].Band  = Band;
		Keys[*nKeys].Chan  = Chan;
		Keys[*nKeys].Val   = 2.0f*ValNp + AnalysisPowerNp;
		(*nKeys)++;
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
		//! Calculate the energy leading up to the transient point

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
		//! NOTE: Do NOT count 'release' jumps; these are very easily
		//! masked and a larger overlap generally provides better
		//! perceived audio quality.
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
			if(d > 0.0f) {
				EnergyCenter += s*SQR(d);
				EnergyEdge   += c*SQR(d);
			}
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
		if(EnergyCenter > *LastBlockEnergy) {
			float LogRatio = (EnergyCenter / EnergyEdge) * logf(EnergyCenter / *LastBlockEnergy);
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

//! Apply block transform
//!  -Fetches data
//!  -Applies MDCT
//!  -Stores keys for block coefficients
//! Returns the number of keys stored
static inline void Block_Transform_ScaleAndToNepers(float *Dst, float *Src, int BlockSize) {
	int Band;
	for(Band=0;Band<BlockSize;Band++) {
		float v = Src[Band] * 2.0f/BlockSize;
		Dst[Band] = (ABS(v) < 0.5f*ULC_COEF_EPS) ? ULC_COEF_NEPER_OUT_OF_RANGE : logf(ABS(v));
		Src[Band] = v;
	}
}
static int Block_Transform(struct ULC_EncoderState_t *State, const float *Data, float PowerDecay) {
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;
	float AnalysisPowerNp = 0.0f; PowerDecay = logf(PowerDecay);

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
	int Chan;
	int nKeys = 0;
	for(Chan=0;Chan<nChan;Chan++) {
		float *BufferTransform = State->TransformBuffer[Chan];
		float *BufferNepers    = State->TransformNepers[Chan];
		float *BufferFwdLap    = State->TransformFwdLap[Chan];
		float *BufferTemp      = State->TransformTemp;

		//! Apply transforms and insert keys
		Fourier_MDCT(BufferTransform, Data + Chan*BlockSize, BufferFwdLap, BufferTemp, BlockSize, BlockSize >> OverlapScale);
		Block_Transform_ScaleAndToNepers(BufferNepers, BufferTransform, BlockSize);
		Block_Transform_InsertKeys(
			State->AnalysisKeys,
			BufferTransform,
			BufferNepers,
			BlockSize,
			&nKeys,
			Chan,
			AnalysisPowerNp,
			State->RateHz * 0.5f
		);
		AnalysisPowerNp += PowerDecay;
	}
	Analysis_KeysValSort(State->AnalysisKeys, nKeys);
	return nKeys;
}

/**************************************/
//! EOF
/**************************************/
