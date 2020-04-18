/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
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
#if ULC_USE_PSYCHOACOUSTICS
	struct Block_Transform_MaskingState_t MaskingState;
	Block_Transform_MaskingState_Init(&MaskingState, Coef, CoefNp, BlockSize, NyquistHz);
#else
	(void)NyquistHz;
#endif
	int   Band;
	int   QBand     = 0;
	float QBandAvg  = 0.0f;
	float QBandAvgW = 0.0f;
	for(Band=0;Band<BlockSize;Band++) {
		//! Check that the value is in range of the smallest quantization
		float ValNp = CoefNp[Band]; if(ValNp == ULC_COEF_NEPER_OUT_OF_RANGE) continue;

		//! Check the 'background' level of this quantizer band against the current value
		if((ValNp + QuantRange)*QBandAvgW < QBandAvg || (ValNp - QuantRange)*QBandAvgW > QBandAvg) {
			//! Out of range - split a new quantizer band
			if(QBand < ULC_MAX_QBANDS-1) {
				QBandAvg = 0.0f, QBandAvgW = 0.0f;
				QBand++;
			}
		}
		QBandAvg  += SQR(Coef[Band])*ValNp;
		QBandAvgW += SQR(Coef[Band]);

		//! Insert key for this band
#if ULC_USE_PSYCHOACOUSTICS
		//! NOTE: Not sure why this masking equation is the way it is.
		//! Using 2*ValNp-Mask does not give very impressive results
		//! whereas this trial-and-error form gives substantially
		//! better results (values correspond to 30dB and 22dB in Np).
		//! NOTE: Reduce importance of non-tonal/non-noise bands by 17.37dB.
		float Flat, Mask = Block_Transform_UpdateMaskingThreshold(&MaskingState, Coef, CoefNp, Band, BlockSize, &Flat);
		ValNp  = 0x1.BA18AAp1f*ValNp - 0x1.443438p1f*Mask;
		ValNp += 2.0f * 4.0f*SQR(Flat)*(SQR(Flat) - 1.0f);
#endif
		//! NOTE: Store the SQUARED post-masking energy as weights.
		Keys[*nKeys].Band  = Band;
		Keys[*nKeys].Chan  = Chan;
		Keys[*nKeys].QBand = QBand;
		Keys[*nKeys].Val   = expf(2.0f*ValNp + AnalysisPowerNp);
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
		const unsigned int SubBlockLength = 32u;
		int nSubBlocks = BlockSize / SubBlockLength;

		//! Divide block into short subblocks and get their sum of squares
		int n;
		float *SubBlockRMS = State->TransformTemp;
		for(n=0;n<nSubBlocks;n++) SubBlockRMS[n] = 1.0e-30f;
		for(Chan=0;Chan<nChan;Chan++) for(n=0;n<nSubBlocks;n++) {
			unsigned int i;
			float Sum = 0.0f;
			const float *Src = Data + Chan*BlockSize + n*SubBlockLength;
			for(i=0;i<SubBlockLength;i++) Sum += SQR(Src[i]);
			SubBlockRMS[n] += Sum;
		}

		//! Find minimum (or maximum) ratio between subblocks
		//! NOTE: The ratios are squared by design; this seems to
		//! result in greater transient selectivity and is faster
		//! than using a square root everywhere
		float RMSa = State->LastTrackedRMS;
		float RatioMin = 1.0f, RatioSum = 1.0e-30f;
		for(n=0;n<nSubBlocks;n++) {
			//! Get the ratio of A:B (or B:A) and save the smallest
			//! NOTE: Even though the original RMS values were scaled
			//! by SubBlockLength*nChan (from the analysis), this will
			//! cancel in the division applied.
			//! NOTE: When A:B > 1.0, then we detected a release/decay
			//! transient; these are less important than attack ones,
			//! so we apply a curve to reduce their importance a bit.
			float Ratio, RMSb = SubBlockRMS[n];
			if(RMSa < RMSb) Ratio = RMSa/RMSb;
			else {
				Ratio = 1.0f - RMSb/RMSa;
				Ratio = 1.0f - SQR(Ratio);
			}
			if(Ratio < RatioMin) RatioMin = Ratio;
			RatioSum += SQR(Ratio); //! Squaring here seems to make it less oversensitive

			//! Remove part of the linear term from the last block,
			//! based on how large of a jump there was between them.
			//! This reduces overdetection during volume envelopes
			RMSa = Ratio*RMSa + (1.0f-Ratio)*RMSb;
		}
		State->LastTrackedRMS = RMSa;
		float Ratio = RatioMin*nSubBlocks / RatioSum;

		//! Set overlap size from the smallest (or largest) ratio,
		//! taking into account its step behaviour
		//! NOTE: The rounding point is at 0.75, and NOT 0.5 as this
		//! would result in too much unnecessary narrowing.
		OverlapScale = (int)(-0x1.715476p0f*logf(Ratio) + 0.25f); //! 0x1.715476p0 = 1/Log[2], to get the log base-2
		if(OverlapScale < 0x0) OverlapScale = 0x0;
		if(OverlapScale > 0xF) OverlapScale = 0xF;
		while((BlockSize >> OverlapScale) < State->MinOverlap) OverlapScale--;
		while((BlockSize >> OverlapScale) > State->MaxOverlap) OverlapScale++;
	}
	State->ThisOverlap = OverlapScale;

	//! Neper range before a quantizer change should happen
	//! 0x1.25701Bp2 = Log[(2*7)^2 / 2]; Half the range of quantized coefficients, in Nepers (39.8dB)
	float QuantRange = 0x1.25701Bp2f * (2.0f - RateKbps/MaxCodingKbps(BlockSize, nChan, State->RateHz));
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
