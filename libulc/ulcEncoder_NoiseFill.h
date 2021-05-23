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
static int Block_Encode_EncodePass_GetNoiseQ(const float *LogCoef, int Band, int N, float q) {
	//! Analyze for the noise amplitude (geometric mean)
	float Amplitude = 0.0f; {
		int n;
		for(n=0;n<N;n++) Amplitude += LogCoef[Band+n];
		Amplitude = expf(Amplitude/N);
	}

	//! Quantize the noise amplitude into final code
	int NoiseQ = Block_Encode_Quantize(Amplitude, q, 1);
	if(NoiseQ > 0x7) NoiseQ = 0x7;
	return NoiseQ;
}

//! Compute quantized HF extension parameters for encoding
static void Block_Encode_EncodePass_GetHFExtParams(const float *LogCoef, int Band, int N, float q, int *_NoiseQ, int *_NoiseDecay) {
	//! Solve for weighted least-squares (in the log domain, for exponential fitting)
	float Amplitude, Decay; {
		int n;
		float SumX  = 0.0f;
		float SumX2 = 0.0f;
		float SumXY = 0.0f;
		float SumY  = 0.0f;
		float SumW  = 0.0f;
		float w = 1.0f;
		float SpectralDecayRate = expf(-0x1.7069E3p3f / N); //! E^(-100/20*Log[10]/n) (ie. decay by 100dB upon reaching the end)
		for(n=0;n<N;n++) {
			float x = (float)n;
			float yLog = LogCoef[Band+n];
			SumX  += w*x;
			SumX2 += w*x*x;
			SumXY += w*x*yLog;
			SumY  += w*  yLog;
			SumW  += w;
			w *= SpectralDecayRate;
		}
		float Det = SumW*SumX2 - SQR(SumX);
		if(Det == 0.0f) { *_NoiseQ = *_NoiseDecay = 0; return; }
		Amplitude = expf((SumX2*SumY  - SumX*SumXY) / Det);
		Decay     = expf((SumW *SumXY - SumX*SumY ) / Det);
		if(Decay > 1.0f) Decay = 1.0f;
	}

	//! Quantize amplitude and decay
	int NoiseQ     = Block_Encode_Quantize(Amplitude, q, 1);
	int NoiseDecay = Block_Encode_Quantize(1.0f-Decay, SQR(128.0f), 0);
	if(NoiseDecay > 50) NoiseQ = 0; //! When decay is too steep (half of max decay), disable fill
	else {
		if(NoiseQ     >  0x7+1) NoiseQ     =  0x7+1;
		if(NoiseDecay > 0x1F+1) NoiseDecay = 0x1F+1;
	}
	*_NoiseQ     = NoiseQ;
	*_NoiseDecay = NoiseDecay;
}

/**************************************/
//! EOF
/**************************************/
