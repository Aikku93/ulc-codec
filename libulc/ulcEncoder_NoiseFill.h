/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2023, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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
//! The main difference is that we're extracting the noise
//! level after masking with the tone level, rather than
//! the other way around.
static inline void Block_Transform_CalculateNoiseLogSpectrum(float *Data, void *Temp, int N, int RateHz) {
	const int ULC_N_BARK_BANDS = 25;
	float NyquistHz = (float)RateHz * 0.5f;

	//! DCT+DST -> Pseudo-DFT
	N /= 2;

	//! Compute logarithm for all lines to speed up calculations
	float *LogData = (float*)Temp; {
		int Line;
		for(Line=0;Line<N;Line++) {
			LogData[Line] = logf(0x1.0p-127f + Data[Line]);
		}
	}

	//! Iterate over all Bark bands
	int BarkBand;
	float *BarkMask = LogData + N;
	for(BarkBand=0;BarkBand<ULC_N_BARK_BANDS;BarkBand++) {
		//! Get the lines corresponding to this Bark band
		float FreqBeg = BarkToFreq(BarkBand+0);
		float FreqEnd = BarkToFreq(BarkBand+1);
		int   LineBeg = (int)floorf(FreqToLine(FreqBeg, NyquistHz, N));
		int   LineEnd = (int)ceilf (FreqToLine(FreqEnd, NyquistHz, N));
		if(LineBeg < 0) LineBeg = 0;
		if(LineEnd < 0) LineEnd = 0;
		if(LineBeg > N-1) LineBeg = N-1;
		if(LineEnd > N)   LineEnd = N;

		//! Sum levels for this band
		double SumFloor = 0.0;
		double SumPeak  = 0.0;
		double SumPeakW = 0.0;
		int nLines = LineEnd - LineBeg;
		if(nLines > 0) {
			int Line;
			const float *Src    = Data    + LineBeg;
			const float *SrcLog = LogData + LineBeg;
			for(Line=0;Line<nLines;Line++) {
				double v    = (double)Src   [Line];
				double vLog = (double)SrcLog[Line];
				SumFloor += vLog;
				SumPeak  += vLog * v;
				SumPeakW += v;
			}
		}

		//! Get the final noise ratio for this band
		float MaskRatio = 0.0f;
		if(SumPeakW != 0.0) {
			SumPeak   = SumPeak  / SumPeakW;
			SumFloor  = SumFloor / (float)nLines;
			MaskRatio = (float)(SumFloor - SumPeak);
		}
		BarkMask[BarkBand] = MaskRatio;
	}

	//! Now generate noise level for each frequency line
	//! NOTE: If the tone-to-noise ratio is too high, we assume
	//! that this line must NOT be noise-coded. The higher weight
	//! given to these lines ensures that the large negative log
	//! value we feed it will collapse any noise fill that tries
	//! to use it.
	int Line;
	for(Line=0;Line<N;Line++) {
		float BarkBand = FreqToBark(LineToFreq(Line, NyquistHz, N));
		int   BandIdx  = (int)BarkBand;
		      BandIdx  = (BandIdx >=               0) ? BandIdx : 0;
		      BandIdx  = (BandIdx < ULC_N_BARK_BANDS) ? BandIdx : (ULC_N_BARK_BANDS-1);
		float w     = expf(BarkMask[BandIdx]);
		float Noise = ((w < (float)M_SQRT1_2) ? (LogData[Line] + BarkMask[BandIdx]) : -100.0f);
		Data[Line*2+0] = w;
		Data[Line*2+1] = w * (Noise*0.5f + 0x1.62E430p-1f); //! Pre-scale by Scale=4.0/2 for noise quantizer (by adding Log[Scale]));
	}

}

/**************************************/

//! Get the quantized noise amplitude for encoding
static int Block_Encode_EncodePass_GetNoiseQ(const float *Data, int Band, int N, float q) {
	//! Fixup for DCT+DST -> Pseudo-DFT
	Data += Band / 2 * 2;
	N = (N + (Band & 1) + 1) / 2;

	//! Analyze for the noise amplitude (geometric mean over N coefficients)
	float Amplitude; {
		int n;
		float Sum = 0.0f, SumW = 0.0f;
		for(n=0;n<N;n++) {
			float w  = Data[n*2+0];
			float wy = Data[n*2+1]; //! wy = w * y
			Sum += wy, SumW += w;
		}
		if(Sum == 0.0f) return 0;
		Amplitude = expf(Sum/SumW);
	}

	//! Quantize the noise amplitude into final code
	int NoiseQ = ULC_CompandedQuantizeCoefficientUnsigned(Amplitude*q, 1 + 0x7);
	return NoiseQ;
}

//! Compute quantized HF extension parameters for encoding
static int Block_Encode_EncodePass_GetHFExtParams_LeastSquares(const float *Data, int N, float *Amplitude, float *Decay) {
	//! Assumes Data[] is {Weight,Value} pairs
	int n;
	float SumX  = 0.0f;
	float SumX2 = 0.0f;
	float SumXY = 0.0f;
	float SumY  = 0.0f;
	float SumW  = 0.0f;
	for(n=0;n<N;n++) {
		float x  = n * 2.0f;
		float w  = Data[n*2+0];
		float wy = Data[n*2+1]; //! wy = w * y
		SumX  += w*x;
		SumX2 += w*x*x;
		SumXY +=   x*wy;
		SumY  +=     wy;
		SumW  += w;
	}

	//! Solve for amplitude and decay
	float Det = SumW*SumX2 - SQR(SumX); if(Det == 0.0f) return 0;
	*Amplitude = (SumX2*SumY  - SumX*SumXY) / Det;
	*Decay     = (SumW *SumXY - SumX*SumY ) / Det;
	return 1;
}
static void Block_Encode_EncodePass_GetHFExtParams(const float *Data, int Band, int N, float q, int *_NoiseQ, int *_NoiseDecay) {
	//! Fixup for DCT+DST -> Pseudo-DFT
	Data += Band / 2 * 2;
	N = (N + (Band & 1) + 1) / 2;

	//! Solve for least-squares (in the log domain, for exponential fitting)
	float Amplitude, Decay; {
		//! If couldn't solve, disable noise fill
		if(!Block_Encode_EncodePass_GetHFExtParams_LeastSquares(Data, N, &Amplitude, &Decay)) {
			*_NoiseQ = *_NoiseDecay = 0;
			return;
		}

		//! Convert to linear units
		Amplitude = expf(Amplitude);
		Decay     = (Decay < 0.0f) ? expf(Decay) : 1.0f; //! <- Ensure E^LogDecay is <= 1.0
	}

	//! Quantize amplitude and decay
	//! Amplitude has already been scaled by 4.0 (plus normalization),
	//! but we need to scale to 16.0 here because HF extension uses
	//! a 4bit amplitude instead of 3bit like "normal" noise fill does
	int NoiseQ     = ULC_CompandedQuantizeCoefficientUnsigned(Amplitude*q*4.0f, 1 + 0xF);
	int NoiseDecay = ULC_CompandedQuantizeUnsigned((Decay-1.0f) * -0x1.0p19f); //! (1-Decay) * 2^19
	if(!NoiseDecay) return; //! <- On immediate falloff, disable fill
	if(NoiseDecay > 0xFF) NoiseDecay = 0xFF;
	*_NoiseQ     = NoiseQ;
	*_NoiseDecay = NoiseDecay;
}

/**************************************/
//! EOF
/**************************************/
