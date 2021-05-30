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
//! NOTE: N must be at least ULC_HELPER_SUBBLOCK_INTERLEAVE_MODULO.
static int Block_Encode_EncodePass_GetNoiseQ(const float *LogCoef, int Band, int N, float q, int WindowCtrl) {
	//! Analyze for the noise amplitude
	float Amplitude; {
		//! NOTE: The purpose of putting the values into bins is
		//! to avoid noise-fill pre-echo by getting the noise
		//! amplitudes at all subblock positions and then using
		//! a measure that favours low amplitudes (eg. minimum).
		//! NOTE: This mapping process is especially suited to
		//! removing transients in the squared domain. When the
		//! transients are removed in the linear domain, this
		//! tends to give very low estimates for noise. Without
		//! the mapping, and removing bias in the squared domain,
		//! noise leaks into prior bins, causing pre-echo leakage.
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
			float sW = SumW[n];
			float a  = s/sW;
			if(a < Amplitude) Amplitude = a;
		}
		Amplitude = expf(Amplitude);
	}

	//! Quantize the noise amplitude into final code
	int NoiseQ = Block_Encode_Quantize(Amplitude*q * 8.0f);
	if(NoiseQ > 0x7) NoiseQ = 0x7;
	return NoiseQ;
}

//! Compute quantized HF extension parameters for encoding
//! NOTE: N must be at least ULC_HELPER_SUBBLOCK_INTERLEAVE_MODULO.
static void Block_Encode_EncodePass_GetHFExtParams(const float *LogCoef, int Band, int N, float q, int WindowCtrl, int *_NoiseQ, int *_NoiseDecay) {
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
		float w = 1.0f, SpectralDecayRate = expf(-0x1.7069E3p3f / N); //! E^(-100/20*Log[10]/n) (ie. decay by 100dB upon reaching the end)
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
			w *= SpectralDecayRate;
		}

		//! Solve for each bin, then take the minimum amplitude
		//! of each bin, alongside the geometric mean decay
		Amplitude = 1.0f; //! <- Arbitrarily large, but shouldn't need to be larger than 1.0
		Decay     = 0.0f;
		int nSubBlocks = ULC_Helper_SubBlockCount(WindowCtrl >> 4);
		for(n=0;n<nSubBlocks;n++) {
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
	int NoiseQ     = Block_Encode_Quantize(Amplitude*q * 8.0f);
	int NoiseDecay = Block_Encode_Quantize((1.0f-Decay) * 0x1.0p16f);
	if(NoiseDecay > 50) NoiseQ = 0; //! When decay is too steep (half of max decay), disable fill
	else {
		if(NoiseQ     >  0x7+1) NoiseQ     =  0x7+1;
		if(NoiseDecay < 0x00+1) NoiseDecay = 0x00+1;
		if(NoiseDecay > 0x1F+1) NoiseDecay = 0x1F+1;
	}
	*_NoiseQ     = NoiseQ;
	*_NoiseDecay = NoiseDecay;
}

/**************************************/
//! EOF
/**************************************/
