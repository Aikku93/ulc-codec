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
	float NyquistHz,
	float QuantRange
) {
	struct Block_Transform_MaskingState_t MaskingState;
	Block_Transform_MaskingState_Init(&MaskingState, Coef, CoefNp, BlockSize, NyquistHz);

	int   Band;
	int   nQBands   = 0;
	float QBandAvg  = 0.0f;
	float QBandAvgW = 0.0f;
	for(Band=0;Band<BlockSize;Band++) {
		//! Check that the value is in range of the smallest quantization
		float ValNp = CoefNp[Band]; if(ValNp == ULC_COEF_NEPER_OUT_OF_RANGE) continue;

		//! Get masking and equal-loudness parameters to get the final value
		float Flat, Mask = Block_Transform_UpdateMaskingThreshold(&MaskingState, Coef, CoefNp, Band, BlockSize, &Flat);
		float InvAWeight = (Flat-1.0f)*AWeightNp((Band+0.5f)*NyquistHz / BlockSize);
		ValNp += InvAWeight;
		Mask  += InvAWeight;
		ValNp  = (2.0f*ValNp - Mask);
		float ValMasked = expf(ValNp + AnalysisPowerNp);

		//! Check the 'background' level of this quantizer band against the masking threshold
		//! NOTE: Accumulation is done with weighting of the band power to simulate
		//! the fact that low-power bands are masked by loud bands
		if(QBandAvg > (Mask+QuantRange)*QBandAvgW) {
			if(nQBands < ULC_MAX_QBANDS-1) {
				QBandAvg = 0.0f, QBandAvgW = 0.0f;
				nQBands++;
			}
		}
		QBandAvg  += ValMasked*ValNp;
		QBandAvgW += ValMasked;

		//! Insert key for this band
		//! NOTE: Store the 'audible' level to the key value. This
		//! will be used later as a weight for the quantizers.
		Keys[*nKeys].Band  = Band;
		Keys[*nKeys].Chan  = Chan;
		Keys[*nKeys].QBand = nQBands;
		Keys[*nKeys].Val   = ValMasked;
		(*nKeys)++;
	}
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
		//! Divide block into 64-sample subblocks and get their sum of squares
		int n;
		float *SubBlockRMS = State->TransformTemp;
		for(n=0;n<64;n++) SubBlockRMS[n] = 1.0e-30f;
		for(Chan=0;Chan<nChan;Chan++) {
			for(n=0;n<BlockSize;n++) SubBlockRMS[n/64] += SQR(Data[Chan*BlockSize+n]);
		}

		//! Find minimum (or maximum) ratio between subblocks
		//! NOTE: The ratios are squared by design; this seems to
		//! result in greater transient selectivity
		float MinRatio = 1.0f;
		float RMSa = State->LastTrackedRMS;
		for(n=0;n<BlockSize/64;n++) {
			//! Get the ratio of A:B (or B:A) and save the smallest
			//! NOTE: Even though the original RMS values were scaled
			//! by 64*nChan (from the analysis), this cancels in the
			//! division applied, only needing the square root for RMS.
			//! NOTE: When A:B > 1.0, then we detected a release/decay
			//! transient; these are less important than attack ones,
			//! so we apply a curve to reduce their importance a bit.
			float RMSb = SubBlockRMS[n];
			float Ratio = RMSa/RMSb;
			if(Ratio > 1.0f) {
				Ratio = 1.0f - RMSb/RMSa;
				Ratio = 1.0f - SQR(SQR(Ratio));
			}
			if(Ratio < MinRatio) MinRatio = Ratio;

			//! Remove part of the linear term from the last block,
			//! based on how large of a jump there was between them.
			//! This reduces overdetection during volume envelopes
			RMSa = Ratio*RMSa + (1.0f-Ratio)*RMSb;
		}
		State->LastTrackedRMS = RMSa;

		//! Set overlap size from the smallest (or largest) ratio
		//! NOTE: Maximum overlap time of 50ms. The rounding point
		//! is also at 0.75, and NOT 0.5 as this would be too much.
		float OverlapSec = (50.0f/1000.0f) * MinRatio;
		float OverlapSmp = OverlapSec*State->RateHz;
		OverlapScale = (int)(log2f(BlockSize / (1.0e-30f + OverlapSmp)) + 0.25f);
		if(OverlapScale < 0x0) OverlapScale = 0x0;
		if(OverlapScale > 0xF) OverlapScale = 0xF;
		while((BlockSize >> OverlapScale) < State->MinOverlap) OverlapScale--;
		while((BlockSize >> OverlapScale) > State->MaxOverlap) OverlapScale++;
	}
	State->ThisOverlap = OverlapScale;

	//! Neper range before a quantizer change should happen
	//! 0x1.51CCA1p2 = Log[(2*7)^2]; Range of quantized coefficients, in Nepers (45.8dB)
	float QuantRange = 0x1.51CCA1p2f * (1.0f - RateKbps/MaxCodingKbps(BlockSize, nChan, State->RateHz));
	if(QuantRange < 0.0f) QuantRange = 0.0f;
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
