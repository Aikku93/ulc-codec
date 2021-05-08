/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <math.h>
/**************************************/
#include "ulcHelper.h"
/**************************************/

//! Get the quantized noise amplitude for encoding
static int Block_Encode_EncodePass_GetNoiseQ(const float *LogCoef, int Band, int N, int WindowCtrl, float q) {
	//! Analyze for the noise amplitude
	float Amplitude; {
		//! NOTE: The purpose of putting the values into bins is
		//! to avoid noise-fill pre-echo by getting the noise
		//! amplitudes at all subblock positions and then using
		//! a measure that favours low amplitudes (eg. minimum).
		int n;
		const int8_t *Mapping = ULC_Helper_SubBlockInterleavePattern(WindowCtrl >> 4);

		//! Perform the actual analysis (geometric mean in each bin)
		float Sum [ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumW[ULC_MAX_SUBBLOCKS] = {0.0f};
		for(n=0;n<N;n++) {
			int Bin = Mapping[(Band+n) % ULC_HELPER_SUBBLOCK_INTERLEAVE_MODULO];
			Sum [Bin] += LogCoef[Band+n];
			SumW[Bin] += 1.0f;
		}

		//! Find the minimum bin amplitude
		Amplitude = 1.0f; //! <- Arbitrarily large, but shouldn't need to be larger than 1.0
		int nSubBlocks = ULC_Helper_SubBlockCount(WindowCtrl >> 4);
		for(n=0;n<nSubBlocks;n++) {
			float s  = Sum [n];
			float sW = SumW[n]; if(sW == 0.0f) continue; //! Bin not affected by fill - ignore it
			float a  = s/sW;
			if(a < Amplitude) Amplitude = a;
		}
		Amplitude = expf(Amplitude);
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
static void Block_Encode_EncodePass_GetHFExtParams(const float *LogCoef, int Band, int N, int WindowCtrl, float q, int *_NoiseQ, int *_NoiseDecay) {
	//! Solve for weighted least-squares (in the log domain, for exponential fitting)
	float Amplitude, Decay; {
		//! As with normal noise-fill, group everything into bins
		int n;
		const int8_t *Mapping = ULC_Helper_SubBlockInterleavePattern(WindowCtrl >> 4);

		//! Next, accumulate all the data into their respective bins.
		//! NOTE: The analysis is the same as in normal noise-fill.
		//! But this time carrying a weight that decays the higher
		//! in frequency we go (with the idea being that these
		//! reconstructed frequencies will matter less and less,
		//! so long as the start is close enough to the original).
		float w = 1.0f;
		float SumX [ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumX2[ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumXY[ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumY [ULC_MAX_SUBBLOCKS] = {0.0f};
		float SumW [ULC_MAX_SUBBLOCKS] = {0.0f};
		for(n=0;n<N;n++) {
			int   Bin = Mapping[(Band+n) % ULC_HELPER_SUBBLOCK_INTERLEAVE_MODULO];
			float x = (float)n;
			float yLog = LogCoef[Band+n];
			SumX [Bin] += w*x;
			SumX2[Bin] += w*x*x;
			SumXY[Bin] += w*x*yLog;
			SumY [Bin] += w  *yLog;
			SumW [Bin] += w;
			w *= 0.99f;
		}

		//! Solve for each bin, then take the minimum amplitude
		//! of each bin, alongside the geometric mean decay
		Amplitude = 1.0f; //! <- Arbitrarily large, but shouldn't need to be larger than 1.0
		Decay     = 0.0f;
		int nSubBlocks = ULC_Helper_SubBlockCount(WindowCtrl >> 4);
		for(n=0;n<nSubBlocks;n++) {
			if(SumW[n] == 0.0f) continue; //! Bin not affected by fill - ignore it
			float Det = SumW[n]*SumX2[n] - SQR(SumX[n]);
			if(Det != 0.0f) {
				float a = (SumX2[n]*SumY [n] - SumX[n]*SumXY[n]) / Det;
				float d = (SumW [n]*SumXY[n] - SumX[n]*SumY [n]) / Det;
				if(a < Amplitude) Amplitude = a;
				Decay += d;
			} else {
				//! Play it safe and disable HF extension
				*_NoiseQ = *_NoiseDecay = 0;
				return;
			}
		}
		Amplitude = expf(Amplitude);
		Decay     = expf(Decay/nSubBlocks);
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
