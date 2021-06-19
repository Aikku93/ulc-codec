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

//! In calculating Log[Re^2+Im^2], we used a filter to remove
//! transient bias, which resulted in an overall gain of Sqrt[2].
//! Sqrt[Re^2+Im^2] gives amplitudes that are generally scaled
//! by Sqrt[2] for noise (relative to Re, which is where we're
//! substituting noise), giving a scale of 2.0. Finally, we will
//! use a geometric mean for determining noise fill levels, which
//! converges to 1/E for random noise, meaning we should scale by
//! E/2.0 (0x1.5BF0A9p0) for deciding on the noise fill amplitude.

/**************************************/

//! Compute noise spectrum (logarithmic output)
static inline void Block_Transform_CalculateNoiseLogSpectrum(float *LogNoise, float *Power, int N) {
	int n;
	float v;

	//! Find the subblock's normalization factor
	float LogNorm = 0.0f;
	for(n=0;n<N;n++) if((v = Power[n]) > LogNorm) LogNorm = v;
	if(LogNorm == 0.0f) {
		//! Empty spectrum - fill with -100.0Np data just in case
		for(n=0;n<N;n++) LogNoise[n] = -100.0f;
		return;
	}

	//! Normalize the logarithmic energy and convert to fixed-point
	//! NOTE: Strictly speaking, we could normalize by 1/LogNorm here.
	//! However, normalizing to 2^32 should reduce subnormal collapse.
	LogNorm = 0x1.0p32f / LogNorm;
	uint32_t *LogPower = (uint32_t*)Power;
	for(n=0;n<N;n++) {
		v = Power[n] * LogNorm;
		v = (v > 1.0f) ? ceilf(logf(v) * 0x1.715476p27f) : 0.0f;
		LogPower[n] = (v >= 0x1.0p32f) ? 0xFFFFFFFFu : (uint32_t)v;
	}
	LogNorm = -0.5f*logf(LogNorm); //! Scale by 1/2 to convert Power to Amplitude

	//! Thoroughly smooth/flatten out the spectrum for noise analysis.
	//! This is achieved by using a geometric mean over each frequency
	//! band's critical bandwidth. The code is very similar to the one
	//! used in psychoacoustics (see ulcEncoder_Psycho.h for details).
	int BandBeg = 0, BandEnd = 0;
	uint64_t Sum = 0ull;
	for(n=0;n<N;n++) {
		//! Re-focus analysis window
		const int RangeScaleFxp = 2;
		int BandLen; {
			int Old, New;
			const int LoRangeScale = 3; //! Beg = 0.75*Band
			const int HiRangeScale = 5; //! End = 1.25*Band

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

		//! Store the geometric mean for this band
		LogNoise[n] = (Sum / BandLen)*0x1.62E430p-29f + LogNorm; //! 0x1.62E430p-29 = 1/LogScale / 2 (/2 to convert Power to Amplitude)
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
	int NoiseQ = (int)(Amplitude*q * 0x1.5BF0A9p0f);
	if(NoiseQ > 0x7) NoiseQ = 0x7;
	return NoiseQ;
}

//! Compute quantized HF extension parameters for encoding
static void Block_Encode_EncodePass_GetHFExtParams(const float *LogCoef, int Band, int N, float q, int *_NoiseQ, int *_NoiseDecay) {
	//! Solve for least-squares (in the log domain, for exponential fitting)
	float Amplitude, Decay; {
		//! NOTE: The analysis is the same as in normal noise-fill,
		//! but with a decaying weight parameter. This appears to
		//! be necessary to avoid overfitting to -inf dB, but is
		//! still not perfect.
		//! NOTE: We always have at least 16 coefficients before
		//! Band, so we use these to regularize the fit.
		int n;
		float SumX  = 0.0f;
		float SumX2 = 0.0f;
		float SumXY = 0.0f;
		float SumY  = 0.0f;
		float SumW  = 0.0f;
		float w = 1.0f;
		for(n=-16;n<N;n++) {
			float x = (float)n;
			float yLog = LogCoef[Band+n];
			SumX  += w*x;
			SumX2 += w*x*x;
			SumXY += w*x*yLog;
			SumY  += w  *yLog;
			SumW  += w;
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
	int NoiseQ     = (int)(Amplitude*q * 0x1.5BF0A9p0f);
	int NoiseDecay = ULC_CompandedQuantizeUnsigned((Decay-1.0f) * -0x1.0p16f); //! (1-Decay) * 2^16
	if(NoiseQ     > 1 + 0x7) NoiseQ     = 1 + 0x7;
	if(NoiseDecay <     0x0) NoiseDecay =     0x0;
	if(NoiseDecay >    0xFF) NoiseDecay =    0xFF;
	*_NoiseQ     = NoiseQ;
	*_NoiseDecay = NoiseDecay;
}

/**************************************/
//! EOF
/**************************************/
