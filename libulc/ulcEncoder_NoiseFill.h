/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <math.h>
/**************************************/
#include "ulcHelper.h"
/**************************************/

//! Fit an exponential curve using linear least-squares
static int Block_Encode_EncodePass_FitExpCurve(const float *X, const float *Y, int N, float *m, float *c) {
	int n;
	float SumX  = 0.0f, SumX2 = 0.0f;
	float SumY  = 0.0f, SumY2 = 0.0f;
	float SumXY = 0.0f;
	int   SumN  = 0;
	for(n=0;n<N;n++) if(Y[n] != 0.0f) {
		float x = X[n];
		float y = logf(Y[n]);
		SumX  += x, SumX2 += SQR(x);
		SumY  += y, SumY2 += SQR(y);
		SumXY += x*y;
		SumN++;
	}
	float Det = SumN*SumX2 - SQR(SumX); if(Det == 0.0f) return 0;
	*m = expf((SumN*SumXY - SumX*SumY)  / Det);
	*c = expf((SumY*SumX2 - SumX*SumXY) / Det);
	return 1;
}

//! Get the amplitude of noise in a segment (via mean and noise-to-signal ratio)
static float Block_Encode_EncodePass_GetNoiseAmplitude(const float *Coef, int Band, int N, int WindowCtrl) {
#define MAX_BINS 4
	static const int8_t BinMapping[8][8] = {
		{0,0,0,0,0,0,0,0}, //! 000x: N/1
		{0,1,0,1,0,1,0,1}, //! 001x: N/2,N/2
		{0,1,2,2,0,1,2,2}, //! 010x: N/4,N/4,N/2
		{0,0,1,2,0,0,1,2}, //! 011x: N/2,N/4,N/4
		{0,1,2,2,3,3,3,3}, //! 100x: N/8,N/8,N/4,N/2
		{0,0,1,2,3,3,3,3}, //! 101x: N/4,N/8,N/8,N/2
		{0,0,0,0,1,2,3,3}, //! 110x: N/2,N/8,N/8,N/4
		{0,0,0,0,1,1,2,3}, //! 111x: N/2,N/4,N/8,N/8
	};

	//! Analyze the values at all subblocks
	//! NOTE: We really want to emphasize tone-masks-noise
	//! here (ie. in the presence of tones, we want noise to
	//! be extremely quiet), so we analyze in an expanded domain.
	//! NOTE: The purpose of putting the values into bins is
	//! to avoid noise-fill pre-echo by getting the noise
	//! amplitudes at all subblock positions and then using
	//! an average that favours low amplitudes (geometric,
	//! harmonic, etc). Comparing against the previous block
	//! might also be helpful here.
	int n;
	int   SumN[MAX_BINS] = {0};
	float Sum1[MAX_BINS] = {0.0f};
	float Sum2[MAX_BINS] = {0.0f};
	const int8_t *Mapping = BinMapping[WindowCtrl >> (4+1)]; //! Low bit of high nybble only selects L/R overlap scaling
	for(n=0;n<N;n++) {
		int Bin = Mapping[(Band+n) % 8u];
		float c = ABS(Coef[Band+n]);
		Sum1[Bin] += sqrtf(c);
		Sum2[Bin] += c;
		SumN[Bin]++;
	}

	//! Finally take the harmonic mean of the noise amplitudes
	//! of all subblocks and combine it into the final fill level
	float Total  = 0.0f;
	int   TotalW = 0;
	for(n=0;n<MAX_BINS;n++) {
		int   sN = SumN[n]; if(!sN) continue;
		float s1 = Sum1[n]; if(s1 == 0.0f) return 0.0f;
		float s2 = Sum2[n];
		Total += sqrtf(s2 * sN*SQR((float)sN)) / SQR(s1); //! 1/((s1/sN)^2 / Sqrt[s2/sN])
		TotalW++;
	}
	return TotalW ? 2.0f*SQR(TotalW / Total) : 0.0f; //! 0x1.6A09E6p1 = 2*Sqrt[2]
#undef MAX_BINS
}

/**************************************/

//! Get the quantized noise amplitude for encoding
static int Block_Encode_EncodePass_GetNoiseQ(float q, const float *Coef, int Band, int N, int WindowCtrl) {
	//! Quantize the noise amplitude into final code
	//! NOTE: This is encoded at higher precision, because it spans
	//! the full 4bit range, meaning we have an extra Log2[16^2 / 7^2]
	//! bits to play with (2.385 bits, so use 3.0 for simplicity,
	//! especially since noise should be lower than the maximum value).
	float Amplitude = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band, N, WindowCtrl);
	int NoiseQ = (int)sqrtf(Amplitude * 8.0f*q); //! <- Round down
	if(NoiseQ > 0xF+1) NoiseQ = 0xF+1;
	return NoiseQ;
}

//! Compute quantized HF extension parameters for encoding
static void Block_Encode_EncodePass_GetHFExtParams(const float *Coef, float q, int Band, int N, int FirstBand, int WindowCtrl, int *_NoiseQ, int *_NoiseDecay) {
	//! Estimate the amplitude and decay parameters
	float Amplitude = 0.0f;
	float Decay     = 0.0f; {
		int x0End = 0;
		int x1End = N*1/4u;
		int x2End = N*2/4u;
		int x3End = N*3/4u;
		int x4End = N*4/4u;
		int x0Beg = -(N*4/4u); if(Band+x0Beg < FirstBand) x0Beg = FirstBand - Band;
		int x1Beg = -(N*3/4u); if(Band+x1Beg < FirstBand) x1Beg = FirstBand - Band;
		int x2Beg = -(N*2/4u); if(Band+x2Beg < FirstBand) x2Beg = FirstBand - Band;
		int x3Beg = -(N*1/4u); if(Band+x3Beg < FirstBand) x3Beg = FirstBand - Band;
		int x4Beg = -(N*0/4u); if(Band+x4Beg < FirstBand) x4Beg = FirstBand - Band;
		float x[5], y[5];
		x[0] = (x0Beg + x0End) * 0.5f;
		x[1] = (x1Beg + x1End) * 0.5f;
		x[2] = (x2Beg + x2End) * 0.5f;
		x[3] = (x3Beg + x3End) * 0.5f;
		x[4] = (x4Beg + x4End) * 0.5f;
		y[0] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band+x0Beg, x0End-x0Beg, WindowCtrl) + 0x1.0p-32f;
		y[1] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band+x1Beg, x1End-x1Beg, WindowCtrl) + 0x1.0p-32f;
		y[2] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band+x2Beg, x2End-x2Beg, WindowCtrl) + 0x1.0p-32f;
		y[3] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band+x3Beg, x3End-x3Beg, WindowCtrl) + 0x1.0p-32f;
		y[4] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band+x4Beg, x4End-x4Beg, WindowCtrl) + 0x1.0p-32f;
		Block_Encode_EncodePass_FitExpCurve(x, y, 5, &Decay, &Amplitude);
	}

	//! Quantize amplitude and decay
	     if(Decay < SQR(1/32.0f)*0.5f) Amplitude = 0.0f; //! Disable fill with too fast decay
	else if(Decay > 1.0f) Decay = 0.0f;
	else                  Decay = sqrtf(1.0f - Decay); //! <- Re-map
	int NoiseQ     = (int)sqrtf(Amplitude * q);        //! <- Round down
	int NoiseDecay = (int)(Decay*32 - 1 + 0.5f);
	if(NoiseQ     > 0xF) NoiseQ     = 0xF;
	if(NoiseDecay < 0x0) NoiseDecay = 0x0;
	if(NoiseDecay > 0xF) NoiseDecay = 0xF;
	*_NoiseQ     = NoiseQ;
	*_NoiseDecay = NoiseDecay;
}

/**************************************/
//! EOF
/**************************************/
