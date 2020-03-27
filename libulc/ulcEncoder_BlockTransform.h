/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#if defined(__AVX__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include <stdint.h>
/**************************************/
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_Helper.h"
/**************************************/
#include "ulcEncoder_Psycho.h"
/**************************************/

//! Insert keys for block coefficients
//! Returns updated number of keys in list
//! NOTE:
//!  AnalysisPower is used to alter the preference for the currently-being-analyzed channel
struct Block_Transform_InsertKeys_UpdateQuantizers_QuantState_t {
	float Avg,   nAvg;
	float AvgLo, nAvgLo;
	float AvgHi, nAvgHi;
	int   nQBands, LastQBandOfs;
};
static inline int Block_Transform_InsertKeys_UpdateQuantizers(
	struct Block_Transform_InsertKeys_UpdateQuantizers_QuantState_t *State,
	int Band,
	float vNp,
	float QuantRange
) {
	if(State->nQBands == ULC_MAX_QBANDS-1) return 0;

	//! Creating the first quantizer band?
	if(State->Avg == 0.0f) {
		State->Avg   = vNp, State->nAvg   = 1.0f;
		State->AvgLo = vNp, State->nAvgLo = 1.0f;
		State->AvgHi = vNp, State->nAvgHi = 1.0f;
		return 0;
	}

	//! Add to the lower/upper average
#if 0
	if(vNp < State->Avg/State->nAvg) State->AvgLo += vNp, State->nAvgLo += 1.0f;
	else                             State->AvgHi += vNp, State->nAvgHi += 1.0f;
#else //! Math optimization (avoid division)
	if(vNp*State->nAvg < State->Avg) State->AvgLo += vNp, State->nAvgLo += 1.0f;
	else                             State->AvgHi += vNp, State->nAvgHi += 1.0f;
#endif
	State->Avg += vNp, State->nAvg += 1.0f;

	//! Check range against threshold
#if 0
	float r = ABS(4.0f*State->AvgHi/State->nAvgHi - 3.0f*State->AvgLo/State->nAvgLo);
	if(r > QuantRange) {
#else //! Math optimization (avoid division)
	float r = ABS(4.0f*State->AvgHi*State->nAvgLo - 3.0f*State->AvgLo*State->nAvgHi);
	if(r > QuantRange*State->nAvgLo*State->nAvgHi) {
#endif
		//! A very sharp transient inside a small bandwidth
		//! is likely to be a statistical anomaly, so ignore it
		//! and let it stabilize on its own
		int Bw = Band-1 - State->LastQBandOfs;
		if(Bw > 4) {
			//! Reset for a new band starting from the transient band
			State->Avg   = vNp, State->nAvg   = 1.0f;
			State->AvgLo = vNp, State->nAvgLo = 1.0f;
			State->AvgHi = vNp, State->nAvgHi = 1.0f;
			return Bw;
		}
	}
	return 0;
}
static inline void Block_Transform_InsertKeys(
	struct AnalysisKey_t *Keys,
	const float *Coef,
	const float *CoefNp,
	int   BlockSize,
	int  *nKeys,
	uint16_t *QuantsBw,
	int   Chan,
	float AnalysisPowerNp,
	float NyquistHz,
	float QuantRange
) {
#if ULC_USE_PSYCHOACOUSTICS
	struct Block_Transform_MaskingState_t MaskingState;
	Block_Transform_MaskingState_Init(&MaskingState, Coef, CoefNp, BlockSize, NyquistHz);
#else
	(void)Coef;
	(void)NyquistHz;
#endif
	int Band;
	struct Block_Transform_InsertKeys_UpdateQuantizers_QuantState_t QuantState = {
		.Avg   = 0, .nAvg   = 0,
		.AvgLo = 0, .nAvgLo = 0,
		.AvgHi = 0, .nAvgHi = 0,
		.nQBands = 0, .LastQBandOfs = 0,
	};
	for(Band=0;Band<BlockSize;Band++) {
		//! Check that the value is in range of the smallest quantization
		float ValNp = CoefNp[Band]; if(ValNp == ULC_COEF_NEPER_OUT_OF_RANGE) continue;
#if ULC_USE_PSYCHOACOUSTICS
		//! Get masking threshold and spectral flatness
		float Flat, Mask = Block_Transform_UpdateMaskingThreshold(&MaskingState, Coef, CoefNp, Band, BlockSize, &Flat);

		//! Get final channel value, taking masking and flatness into account
		//! Scaling by a flatness curve seems to give better results,
		//! as flatter signals mask much more than tonal ones
		//! TMN ratio: Around -6dB
		ValNp = 4.0f*ValNp - 3.0f*(Mask + 0.7f*(Flat - 1.0f));
#endif
		//! Update quantizer bands
		//! NOTE: Pass the /perceived/ level here
#if ULC_USE_PSYCHOACOUSTICS
		int QuantBw = Block_Transform_InsertKeys_UpdateQuantizers(&QuantState, Band, ValNp, QuantRange);
#else
		int QuantBw = Block_Transform_InsertKeys_UpdateQuantizers(&QuantState, Band, ValNp, QuantRange);
#endif
		if(QuantBw != 0) {
			QuantState.LastQBandOfs = Band;
			QuantsBw[QuantState.nQBands++] = QuantBw;
		}

		//! Insert key for this band
		//! NOTE: Store the 'audible' level to the key value. This
		//! will be used later as a weight for the quantizers.
		Keys[*nKeys].Band = Band;
		Keys[*nKeys].Chan = Chan;
		Keys[*nKeys].Val  = expf(ValNp + AnalysisPowerNp);
		(*nKeys)++;
	}

	//! Create the final quantizer band
	QuantsBw[QuantState.nQBands++] = BlockSize - QuantState.LastQBandOfs;
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
static int Block_Transform(struct ULC_EncoderState_t *State, const float *Data, float RateKbps, float PowerDecay) {
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;

	int Chan;
	int nKeys = 0;
	float AnalysisPowerNp = 0.0f; PowerDecay = logf(PowerDecay);

	//! Attempt to find transients in order to decide on an overlap amount
	int OverlapScale; {
		//! Divide block into 64-sample subblocks and get their RMS
		int n;
		float *SubBlockRMS = State->TransformTemp;
		for(n=0;n<64;n++) SubBlockRMS[n] = 0.0f;
		for(Chan=0;Chan<nChan;Chan++) {
			for(n=0;n<BlockSize;n++) SubBlockRMS[n/64] += SQR(Data[Chan*BlockSize+n]);
		}

		//! Find peak difference between subblocks
		//! NOTE: PeakDif is scaled by 2.0 because we used squared inputs,
		//! and biased by Log[64] because it wasn't normalized. The bias
		//! is unimportant, however, as it cancels out in the subtraction.
		float PeakDif = 0.0f;
		float RMSaNp = State->LastTrackedRMSNp;
		for(n=0;n<BlockSize/64;n++) {
			const float Silence = 1.0e-30f;
			float RMS = SubBlockRMS[n]; if(RMS < Silence) RMS = Silence;
			float RMSbNp = logf(RMS);
			float Dif = RMSbNp - RMSaNp;
			if(Dif < 0.0f) Dif *= -0.5f; //! Release transients are less important
			if(Dif > PeakDif) PeakDif = Dif;

			//! NOTE: Remove part of the linear term from the last block.
			//! This reduces overdetection during volume envelopes
			const float mu = 0.3f;
			RMSaNp = mu*RMSaNp + (1.0f-mu)*RMSbNp;
		}
		State->LastTrackedRMSNp = RMSaNp;

		//! Set overlap size from peak difference
		//! Target values:
		//!  3dB (~0.3Np) = 50.0ms
		//!  6dB (~0.7Np) = 25.0ms
		//!  9dB (~1.0Np) = 12.5ms, etc.
		//! ie.
		//!  OverlapMs = 50 * E^(Log[2] - Np*(Log[2] * (20/Log[10])/3))
		//! Note that because PeakDif is scaled by 2.0, we must scale
		//! by 0.5 to account for that (0x1.00E...p0 instead of 0x1.00E...p1).
		float OverlapSec = (50.0f/1000.0f) * expf(0x1.62E430p-1f - PeakDif*0x1.00E102p0f);
		float OverlapSmp = OverlapSec*State->RateHz;
		OverlapScale = lrintf(log2f(BlockSize / (1.0e-30f + OverlapSmp)));
		if(OverlapScale < 0x0) OverlapScale = 0x0;
		if(OverlapScale > 0xF) OverlapScale = 0xF;
		while((BlockSize >> OverlapScale) < State->MinOverlap) OverlapScale--;
	}
	State->ThisOverlap = OverlapScale;

	//! Neper range before a quantizer change should happen
	//! 0x1.51CCA1p2  = Log[(2*7)^2]; Dynamic range of quantized coefficients, in Nepers
	//! 0x1.62E430p-1 = Log[2]
	float QuantRange = 0x1.51CCA1p2f * (1.0f - 4.0f*RateKbps/MaxCodingKbps(BlockSize, nChan, State->RateHz));
	if(QuantRange < 0x1.62E430p-1f) QuantRange = 0x1.62E430p-1f; //! Limit at dynamic range of 2 units
	for(Chan=0;Chan<nChan;Chan++) {
		//! Get buffer pointers
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
			State->QuantsBw[Chan],
			Chan,
			AnalysisPowerNp,
			State->RateHz * 0.5f,
			QuantRange
		);
		AnalysisPowerNp += PowerDecay;
	}
	Analysis_KeysValSort(State->AnalysisKeys, nKeys);
	return nKeys;
}

/**************************************/
//! EOF
/**************************************/
