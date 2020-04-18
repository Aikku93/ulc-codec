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

//! Get optimal log base-2 overlap scaling for transients
//! The idea is that with reduced overlap, transients need fewer
//! coefficients to sound correct (at the cost of distortion)
static inline int Block_Transform_GetLogOverlapScale(const float *Data, float *SubBlockRMSBuffer, float *LastTrackedRMS, int BlockSize, int MinOverlap, int MaxOverlap, int nChan) {
	const unsigned int SubBlockLength = 32u;
	int nSubBlocks = BlockSize / SubBlockLength;

	//! Divide block into short subblocks and get their sum of squares
	int n, Chan; {
		for(n=0;n<nSubBlocks;n++) SubBlockRMSBuffer[n] = 1.0e-30f;
		for(Chan=0;Chan<nChan;Chan++) {
			for(n=0;n<nSubBlocks;n++) {
				unsigned int i;
				float Sum;
#if defined(__AVX__)
				__m256 vSum = _mm256_setzero_ps();
				for(i=0;i<SubBlockLength;i+=8) {
					__m256 v = _mm256_load_ps(Data); Data += 8;
#if defined(__FMA__)
					vSum = _mm256_fmadd_ps(v, v, vSum);
#else
					v    = _mm256_mul_ps(v, v);
					vSum = _mm256_add_ps(vSum, v);
#endif
				}
				__m128 vSumL = _mm256_extractf128_ps(vSum, 0);
				__m128 vSumH = _mm256_extractf128_ps(vSum, 1);
				vSumL = _mm_add_ps    (vSumL, vSumH);
				vSumH = _mm_movehl_ps (vSumL, vSumL);
				vSumL = _mm_add_ps    (vSumL, vSumH);
				vSumH = _mm_shuffle_ps(vSumL, vSumL, 0x01);
				vSumL = _mm_add_ss    (vSumL, vSumH);
				Sum   = _mm_cvtss_f32 (vSumL);
#elif defined(__SSE__)
				__m128 vSum = _mm_setzero_ps();
				for(i=0;i<SubBlockLength;i+=4) {
					__m128 v = _mm_load_ps(Data); Data += 4;
#if defined(__FMA__)
					vSum = _mm_fmadd_ps(v, v, vSum);
#else
					v    = _mm_mul_ps(v, v);
					vSum = _mm_add_ps(vSum, v);
#endif
				}
				__m128 vSumL = vSum;
				__m128 vSumH;
				vSumH = _mm_movehl_ps (vSumL, vSumL);
				vSumL = _mm_add_ps    (vSumL, vSumH);
				vSumH = _mm_shuffle_ps(vSumL, vSumL, 0x01);
				vSumL = _mm_add_ss    (vSumL, vSumH);
				Sum   = _mm_cvtss_f32 (vSumL);
#else
				Sum = 0.0f;
				for(i=0;i<SubBlockLength;i++) {
					float v = *Data++;
					Sum += SQR(v);
				}
#endif
				SubBlockRMSBuffer[n] += Sum;
			}
		}
	}

	//! Find minimum (or maximum) ratio between subblocks
	//! NOTE: The ratios are squared by design; this seems to
	//! result in greater transient selectivity and is faster
	//! than using a square root everywhere
	float RMSa = *LastTrackedRMS;
	float RatioMin = 1.0f, RatioSum = 1.0e-30f;
	for(n=0;n<nSubBlocks;n++) {
		//! Get the ratio of A:B (or B:A) and save the smallest
		//! NOTE: Even though the original RMS values were scaled
		//! by SubBlockLength*nChan (from the analysis), this will
		//! cancel in the division applied.
		//! NOTE: When A:B > 1.0, then we detected a release/decay
		//! transient; these are less important than attack ones,
		//! so we apply a curve to reduce their importance a bit.
		float Ratio, RMSb = SubBlockRMSBuffer[n];
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
	*LastTrackedRMS = RMSa;
	float Ratio = RatioMin*nSubBlocks / RatioSum;

	//! Set overlap size from the smallest (or largest) ratio,
	//! taking into account its step behaviour
	//! NOTE: The rounding point is at 0.75, and NOT 0.5 as this
	//! would result in too much unnecessary narrowing.
	int OverlapScale = (int)(-0x1.715476p0f*logf(Ratio) + 0.25f); //! 0x1.715476p0 = 1/Log[2], to get the log base-2
	if(OverlapScale < 0x0) OverlapScale = 0x0;
	if(OverlapScale > 0xF) OverlapScale = 0xF;
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
static int Block_Transform(struct ULC_EncoderState_t *State, const float *Data, float RateKbps, float PowerDecay) {
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;
	float AnalysisPowerNp = 0.0f; PowerDecay = logf(PowerDecay);

	//! Get the overlap scaling for this block
	int OverlapScale = State->ThisOverlap = Block_Transform_GetLogOverlapScale(
		Data,
		State->TransformTemp,
		&State->LastTrackedRMS,
		BlockSize,
		State->MinOverlap,
		State->MaxOverlap,
		nChan
	);

	//! Get the allowed dynamic range in a quantizer zone
	//! 0x1.25701Bp2 = Log[(2*7)^2 / 2]; Half the range of quantized coefficients, in Nepers (39.8dB)
	float QuantRangeScale = 2.0f - RateKbps/MaxCodingKbps(BlockSize, nChan, State->RateHz);
	if(QuantRangeScale < 1.0f) QuantRangeScale = 1.0f; //! Avoid creating too many quantizer zones
	float QuantRange = 0x1.25701Bp2f * QuantRangeScale;

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
