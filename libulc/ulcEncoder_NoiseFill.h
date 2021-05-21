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
		for(n=0;n<N;n++) {
			float x = (float)n;
			float yLog = LogCoef[Band+n];
			SumX  += x;
			SumX2 += x*x;
			SumXY += x*yLog;
			SumY  +=   yLog;
		}
		float Det = N*SumX2 - SQR(SumX);
		if(Det == 0.0f) { *_NoiseQ = *_NoiseDecay = 0; }
		Amplitude = expf((SumX2*SumY  - SumX*SumXY) / Det);
		Decay     = expf((    N*SumXY - SumX*SumY ) / Det);
	}

	//! Quantize amplitude and decay
	Decay = (Decay > 1.0f) ? 0.0f : sqrtf(1.0f - Decay); //! <- Re-map
	int NoiseQ     = Block_Encode_Quantize(Amplitude, q, 1);
	int NoiseDecay = (int)ceilf(Decay*32 - 1); //! <- Round up (rounds down after inverse mapping)
	if(NoiseQ     > 0xF) NoiseQ     = 0xF;
	if(NoiseDecay < 0x0) NoiseDecay = 0x0;
	if(NoiseDecay > 0xF) NoiseDecay = 0xF;
	*_NoiseQ     = NoiseQ;
	*_NoiseDecay = NoiseDecay;
}

/**************************************/
//! EOF
/**************************************/
