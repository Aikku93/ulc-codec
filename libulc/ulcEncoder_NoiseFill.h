/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <math.h>
/**************************************/
#include "ulcHelper.h"
/**************************************/

//! Compute noise spectrum (logarithmic output)
//! The code here is very similar to the one used in
//! psychoacoustics (see ulcEncoder_Psycho.h for details).
static inline void Block_Transform_CalculateNoiseLogSpectrum(float *Data, void *Temp, int N) {
	int n;
	float v;

	//! Find the subblock's normalization factor
	float Norm = 0.0f;
	for(n=0;n<N;n++) if((v = Data[n]) > Norm) Norm = v;
	if(Norm == 0.0f) return;

	//! Normalize the logarithmic energy and convert to fixed-point
	Norm = (Norm > 0x1.0p-96f) ? (0x1.FFFFFCp31f / Norm) : 0x1.FFFFFCp127f;
	float LogScale = 0x1.EC709Cp30f / N; //! (2^32/Log[2^32]) / (N * (1-29/32)) = (2^32/Log[2^32] / (1-29/32)) / N
	uint32_t *LogPower = (uint32_t*)Temp;
	for(n=0;n<N;n++) {
		v = Data[n] * Norm;
		LogPower[n] = (v <= 1.0f) ? 0 : (uint32_t)(logf(v) * LogScale);
	}
	float LogNorm     = 0x1.62E430p-1f - logf(Norm); //! Pre-scale by Log[2] to account for noise mean at 0.5
	float InvLogScale = N * 0x1.0A2B24p-31f;

	//! Thoroughly smooth/flatten out the spectrum for noise analysis.
	//! This is achieved by using a geometric mean over each frequency
	//! band's critical bandwidth.
	int BandBeg = 0, BandEnd = 0;
	uint32_t Sum = 0;
	for(n=0;n<N;n++) {
		//! Re-focus analysis window
		const int RangeScaleFxp = 5;
		int BandLen; {
			int Old, New;
			const int LoRangeScale = 29; //! Beg = 0.90625*Band
			const int HiRangeScale = 32; //! End = 1.00000*Band

			//! Remove samples that went out of focus
			Old = BandBeg >> RangeScaleFxp, BandBeg += LoRangeScale;
			New = BandBeg >> RangeScaleFxp;
			if(Old < New) {
				Sum -= LogPower[Old];
			}
			BandLen = New;

			//! Add samples that came into focus
			Old = BandEnd >> RangeScaleFxp, BandEnd += HiRangeScale;
			New = BandEnd >> RangeScaleFxp; if(New > N) New = N;
			if(Old < New) do {
				Sum += LogPower[Old];
			} while(++Old < New);
			BandLen = New - BandLen;
		}

		//! Store the geometric mean
		Data[n] = (Sum/BandLen)*InvLogScale + LogNorm;
	}
}

/**************************************/

//! Get the quantized noise amplitude for encoding
static int Block_Encode_EncodePass_GetNoiseQ(const float *LogCoef, int Band, int N, float q) {
	//! Analyze for the noise amplitude (geometric mean over N coefficients)
	float Amplitude; {
		int n;
		Amplitude = 0.0f;
		for(n=0;n<N;n++) Amplitude += LogCoef[Band+n];
		Amplitude = expf(Amplitude/N);
	}

	//! Quantize the noise amplitude into final code
	int NoiseQ = ULC_CompandedQuantizeUnsigned(Amplitude*q*2.0f);
	if(NoiseQ > 0x7) NoiseQ = 0x7;
	return NoiseQ;
}

//! Compute quantized HF extension parameters for encoding
#pragma GCC push_options
#pragma GCC optimize("fast-math") //! Should improve things, hopefully, maybe... please.
static void Block_Encode_EncodePass_GetHFExtParams(const float *LogCoef, int Band, int N, float q, int *_NoiseQ, int *_NoiseDecay) {
	//! Solve for least-squares (in the log domain, for exponential fitting)
	float Amplitude, Decay; {
		//! NOTE: The analysis is the same as in normal noise-fill,
		//! but with a decaying weight parameter. This appears to
		//! be necessary to avoid overfitting to -inf dB (with a
		//! corresponding spike about x=0), but is not perfect.
		int n;
		float SumX  = 0.0f;
		float SumX2 = 0.0f;
		float SumXY = 0.0f;
		float SumY  = 0.0f;
		float SumW  = 0.0f;
		float w = 15.0f;
		for(n=0;n<N;n++) {
			float tw = 1.0f + w;
			float x = (float)n;
			float yLog = LogCoef[Band+n];
			SumX  += tw*x;
			SumX2 += tw*x*x;
			SumXY += tw*x*yLog;
			SumY  += tw  *yLog;
			SumW  += tw;
			w *= 0.99f;
		}

		//! Solve for amplitude and decay
		float Det = SumW*SumX2 - SQR(SumX);
		if(Det == 0.0f) {
			//! Play it safe and disable HF extension
			*_NoiseQ = *_NoiseDecay = 0;
			return;
		}
		Amplitude = (SumX2*SumY  - SumX*SumXY) / Det;
		Decay     = (SumW *SumXY - SumX*SumY ) / Det;

		//! Convert to linear units
		Amplitude = expf(Amplitude);
		Decay     = expf(Decay);
	}

	//! Quantize amplitude and decay
	int NoiseQ     = ULC_CompandedQuantizeUnsigned(Amplitude*q*8.0f); //! <- Account for dynamic range of NoiseQ vs. coefficients
	int NoiseDecay = ULC_CompandedQuantizeUnsigned((Decay-1.0f) * -0x1.0p16f); //! (1-Decay) * 2^16
	if(NoiseQ     > 1 + 0xF) NoiseQ     = 1 + 0xF;
	if(NoiseDecay <    0x00) NoiseDecay =    0x00;
	if(NoiseDecay >    0xFF) NoiseDecay =    0xFF;
	*_NoiseQ     = NoiseQ;
	*_NoiseDecay = NoiseDecay;
}
#pragma GCC pop_options

/**************************************/
//! EOF
/**************************************/
