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
static inline void Block_Transform_CalculateNoiseLogSpectrum(float *Data, void *Temp1, void *Temp2, int N) {
	int n;
	float v;

	//! Find the subblock's normalization factor
	float Norm = 0.0f;
	for(n=0;n<N;n++) if((v = Data[n]) > Norm) Norm = v;
	if(Norm == 0.0f) return;

	//! Normalize the energy and convert to fixed-point
	//! NOTE: The weights are derived from the coefficient
	//! power, not amplitude. This favours the tone noise
	//! levels (which should be lower) and avoid overdoing
	//! the noise content in tonal passages.
	Norm = (Norm > 0x1.0p-96f) ? (0x1.FFFFFCp31f / Norm) : 0x1.FFFFFCp127f;
	float LogScale = 0x1.83CBE3p28f / N; //! (2^32/Log[2^32]) / (N * (1-22/42)) = (2^32/Log[2^32] / (1-22/42)) / N
	uint32_t *Weight   = (uint32_t*)Temp1;
	uint32_t *LogAmp   = (uint32_t*)Data; //! <- Aliasing of Data[]
	uint32_t *LogFloor = (uint32_t*)Temp2;
	for(n=0;n<N;n++) {
		v = Data[n] * Norm;
		LogAmp[n] = (v <= 1.0f) ? 0 : (uint32_t)(logf(v) * LogScale);
		Weight[n] = (v <= 1.0f) ? 1 : (uint32_t)v;
	}
	float LogNorm     = 0x1.3687AAp1f - 0.5f*logf(Norm); //! Pre-scale by Scale=16.0/Sqrt[2] for noise quantizer (by adding Log[Scale]=0x1.3687AAp1)
	float InvLogScale = N * 0x1.51FDE5p-30f;

	//! Thoroughly smooth/flatten out the spectrum for noise analysis.
	//! This is achieved by using a geometric mean over each frequency
	//! band's critical bandwidth, giving us a rough "noise floor"
	int FloorBeg = 0, FloorEnd = 0; uint32_t FloorSum = 0;
	for(n=0;n<N;n++) {
		//! Re-focus analysis window
		const int RangeScaleFxp = 5;
		int BandLen; {
			//! NOTE: This curve is overkill, but is very resistant
			//! to amplitude spikes/dips across the analysis region,
			//! making things sound a lot nicer/consistent.
			//! NOTE: This has a known problem with harmonic sounds
			//! combined with noise, where the noise amplitude is
			//! overestimated due to not enough contrast between
			//! the tonal spikes and the background noise.
			//! I've tried to compensate for this by using the power
			//! spectrum as the weight (instead of the amplitude),
			//! but it still shows up sometimes.
			int Old, New;
			const int LoRangeScale = 22; //! Beg = 1 - 0.3125*Band
			const int HiRangeScale = 42; //! End = 1 + 0.3125*Band

			//! Remove samples that went out of focus
			Old = FloorBeg >> RangeScaleFxp, FloorBeg += LoRangeScale;
			New = FloorBeg >> RangeScaleFxp;
			if(Old < New) {
				FloorSum -= LogAmp[Old];
			}
			BandLen = New;

			//! Add samples that came into focus
			Old = FloorEnd >> RangeScaleFxp, FloorEnd += HiRangeScale;
			New = FloorEnd >> RangeScaleFxp; if(New > N) New = N;
			if(Old < New) do {
				FloorSum += LogAmp[Old];
			} while(++Old < New);
			BandLen = New - BandLen;
		}
		LogFloor[n] = FloorSum / BandLen;
	}

	//! Finally, take a weighted sum of the above-smoothed "noise floor".
	//! In essence, this forces tones (which get a large weight at their
	//! peak, but should have a low "noise floor" due to neighbouring bins
	//! biasing the geometric mean towards -inf) to spread their low noise
	//! amplitude out across more bins, but noise is preserved (as the
	//! weights are somewhat evenly distributed amongst the peak levels).
	int NoiseBeg = 0, NoiseEnd = 0; uint64_t NoiseSum = 0, NoiseSumW = 0;
	for(n=0;n<N;n++) {
		//! Re-focus analysis window
		const int RangeScaleFxp = 5;
		{
			int Old, New;
			const int LoRangeScale = 27; //! Beg = 1 - 0.15625*Band
			const int HiRangeScale = 37; //! End = 1 + 0.15625*Band

			//! Remove samples that went out of focus
			Old = NoiseBeg >> RangeScaleFxp, NoiseBeg += LoRangeScale;
			New = NoiseBeg >> RangeScaleFxp;
			if(Old < New) {
				NoiseSumW -= Weight[Old];
				NoiseSum  -= Weight[Old] * (uint64_t)LogFloor[Old];
			}

			//! Add samples that came into focus
			Old = NoiseEnd >> RangeScaleFxp, NoiseEnd += HiRangeScale;
			New = NoiseEnd >> RangeScaleFxp; if(New > N) New = N;
			if(Old < New) do {
				NoiseSumW += Weight[Old];
				NoiseSum  += Weight[Old] * (uint64_t)LogFloor[Old];
			} while(++Old < New);
		}
		Data[n] = (NoiseSum/NoiseSumW)*InvLogScale + LogNorm;
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
	int NoiseQ = ULC_CompandedQuantizeCoefficientUnsigned(Amplitude*q, 0xF);
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
	int NoiseQ     = ULC_CompandedQuantizeCoefficientUnsigned(Amplitude*q, 1 + 0xF);
	int NoiseDecay = ULC_CompandedQuantizeUnsigned((Decay-1.0f) * -0x1.0p16f); //! (1-Decay) * 2^16
	if(NoiseDecay < 0x00) NoiseDecay = 0x00;
	if(NoiseDecay > 0xFF) NoiseDecay = 0xFF;
	*_NoiseQ     = NoiseQ;
	*_NoiseDecay = NoiseDecay;
}
#pragma GCC pop_options

/**************************************/
//! EOF
/**************************************/
