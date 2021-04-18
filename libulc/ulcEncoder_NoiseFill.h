/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <math.h>
/**************************************/
#include "ulcHelper.h"
/**************************************/

//! Get the quantized noise amplitude for encoding
static int Block_Encode_EncodePass_GetNoiseQ(const float *Coef, float q, int Band, int N, int WindowCtrl) {
	//! Analyze for the noise amplitude
	float Amplitude; {
		//! NOTE: The purpose of putting the values into bins is
		//! to avoid noise-fill pre-echo by getting the noise
		//! amplitudes at all subblock positions and then using
		//! an average that favours low amplitudes (geometric,
		//! harmonic, etc). Comparing against the previous block
		//! might also be helpful here.
		int n;
		const int8_t *Mapping = ULC_Helper_SubBlockInterleavePattern(WindowCtrl >> 4);

		//! Perform the actual analysis (geometric mean in each bin)
		float Sum [ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumW[ULC_MAX_SUBBLOCKS] = {0.0f};
		for(n=0;n<N;n++) {
			int Bin  = Mapping[(Band+n) % ULC_HELPER_SUBBLOCK_INTERLEAVE_MODULO];
			Sum [Bin] += ULC_FastLnApprox(Coef[Band+n]);
			SumW[Bin] += 1.0f;
		}

		//! Get harmonic mean of all bins
		float Total  = 0.0f;
		int   TotalN = 0;
		int   nSubBlocks = ULC_Helper_SubBlockCount(WindowCtrl >> 4);
		for(n=0;n<nSubBlocks;n++) {
			float s  = Sum [n];
			float sW = SumW[n]; if(sW == 0.0f) return 0;
			Total  += expf(-s/sW);
			TotalN += 1;
		}
		Amplitude = Total ? (TotalN/Total) : 0.0f;
	}

	//! Quantize the noise amplitude into final code
	//! NOTE: This is encoded at higher precision, because it spans
	//! the full 4bit range, meaning we have an extra Log2[16^2 / 7^2]
	//! bits to play with (2.385 bits, so use 3.0 for simplicity,
	//! especially since noise should be lower than the maximum value).
	int NoiseQ = (int)sqrtf(Amplitude * 8.0f*q); //! <- Round down
	if(NoiseQ > 0xF+1) NoiseQ = 0xF+1;
	return NoiseQ;
}

//! Compute quantized HF extension parameters for encoding
static void Block_Encode_EncodePass_GetHFExtParams(const float *Coef, float q, int Band, int N, int WindowCtrl, int *_NoiseQ, int *_NoiseDecay) {
	//! Solve for weighted least-squares (in the log domain, for exponential fitting)
	float Amplitude = 0.0f;
	float Decay     = 0.0f; {
		//! As with normal noise-fill, group everything into bins
		int n;
		const int8_t *Mapping = ULC_Helper_SubBlockInterleavePattern(WindowCtrl >> 4);

		//! Next, accumulate all the data into their respective bins.
		//! NOTE: The analysis is the same as in normal noise-fill.
		float SumX [ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumX2[ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumXY[ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumY [ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumW [ULC_MAX_SUBBLOCKS] = {0.0f};
		int   SumN [ULC_MAX_SUBBLOCKS] = {0};
		for(n=0;n<N;n++) {
			int   Bin = Mapping[(Band+n) % ULC_HELPER_SUBBLOCK_INTERLEAVE_MODULO];
			float x = (float)n;
			float w = SQR(N-n);
			float yLog = ULC_FastLnApprox(Coef[Band+n]);
			SumX [Bin] += w*x;
			SumX2[Bin] += w*x*x;
			SumXY[Bin] += w*x*yLog;
			SumY [Bin] += w  *yLog;
			SumW [Bin] += w;
			SumN [Bin] += 1;
		}

		//! Solve for each bin and accumulate for harmonic mean
		int Total = 0;
		int nSubBlocks = ULC_Helper_SubBlockCount(WindowCtrl >> 4);
		for(n=0;n<nSubBlocks;n++) {
			float Det = SQR(SumX[n]) - SumW[n]*SumX2[n]; //! NOTE: Negated to get E^(-x) for harmonic mean
			if(Det != 0.0f) {
				Amplitude += expf((SumX2[n]*SumY [n] - SumX[n]*SumXY[n]) / Det);
				Decay     += expf((SumW [n]*SumXY[n] - SumX[n]*SumY [n]) / Det);
				Total++;
			} else {
				//! Play it safe and disable HF extension
				Amplitude = 0.0f;
				Decay     = 0.0f;
				Total     = 0;
				break;
			}
		}

		//! Get final values
		if(Total) {
			Amplitude = Total/Amplitude;
			Decay     = Total/Decay;
		}
	}

	//! Quantize amplitude and decay
	Decay = (Decay > 1.0f) ? 0.0f : sqrtf(1.0f - Decay); //! <- Re-map
	int NoiseQ     = (int)sqrtf(Amplitude * q); //! <- Round down
	int NoiseDecay = (int)ceilf(Decay*32 - 1);  //! <- Round up (rounds down after inverse mapping)
	if(NoiseQ     > 0xF) NoiseQ     = 0xF;
	if(NoiseDecay < 0x0) NoiseDecay = 0x0;
	if(NoiseDecay > 0xF) NoiseDecay = 0xF;
	*_NoiseQ     = NoiseQ;
	*_NoiseDecay = NoiseDecay;
}

/**************************************/
//! EOF
/**************************************/
