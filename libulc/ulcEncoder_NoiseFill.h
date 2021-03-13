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
#define MAX_BINS 8
	//! Analyze the values at all subblocks
	//! NOTE: The purpose of putting the values into bins is
	//! to avoid noise-fill pre-echo by getting the noise
	//! amplitudes at all subblock positions and then using
	//! an average that favours low amplitudes (geometric,
	//! harmonic, etc). Comparing against the previous block
	//! might also be helpful here.
	//! Notes:
	//!  The expression:
	//!    Sqrt[2] * (Total[v]^2 / (N*Total[v^2]))
	//!  approximates the noise-to-signal ratio. Therefore,
	//!  multiplying this by the mean (Total[v]/N) gives the
	//!  average level of the noise. Multiplying by 2.0 then
	//!  gives the noise amplitude, resulting in:
	//!    2*Sqrt[2] * Total[v]^3 / (N^2*Total[v^2])
	int n;
	int BinN[MAX_BINS] = {0};
	float Sum[MAX_BINS] = {0.0f}, Sum2[MAX_BINS] = {0.0f};
	for(n=0;n<N;n++) {
		int Bin = (Band+n) % (unsigned int)MAX_BINS;
		float c = ABS(Coef[Band+n]);
		Sum[Bin] += c, Sum2[Bin] += SQR(c), BinN[Bin]++;
	}

	//! Combine the bins into subblocks
	int nBins = 0; //! <- Shuts gcc up
	switch(WindowCtrl >> (4+1)) { //! Low bit of high nybble only selects L/R overlap scaling
		//! N/1
		case 0b000: {
#define COMBINE_A(x) x[0] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7]
			nBins = 1;
			COMBINE_A(Sum), COMBINE_A(Sum2), COMBINE_A(BinN);
#undef COMBINE_A
		} break;

		//! N/2,N/2
		case 0b001: {
#define COMBINE_A(x) x[0] = x[0] + x[2] + x[4] + x[6]
#define COMBINE_B(x) x[1] = x[1] + x[3] + x[5] + x[7]
			nBins = 2;
			COMBINE_A(Sum), COMBINE_A(Sum2), COMBINE_A(BinN);
			COMBINE_B(Sum), COMBINE_B(Sum2), COMBINE_B(BinN);
#undef COMBINE_B
#undef COMBINE_A
		} break;

		//! N/4,N/4,N/2
		case 0b010: {
#define COMBINE_A(x) x[0] = x[0] + x[4]
#define COMBINE_B(x) x[1] = x[1] + x[5]
#define COMBINE_C(x) x[2] = x[2] + x[3] + x[6] + x[7]
			nBins = 3;
			COMBINE_A(Sum), COMBINE_A(Sum2), COMBINE_A(BinN);
			COMBINE_B(Sum), COMBINE_B(Sum2), COMBINE_B(BinN);
			COMBINE_C(Sum), COMBINE_C(Sum2), COMBINE_C(BinN);
#undef COMBINE_C
#undef COMBINE_B
#undef COMBINE_A
		} break;

		//! N/2,N/4,N/4
		case 0b011: {
#define COMBINE_A(x) x[0] = x[0] + x[1] + x[4] + x[5]
#define COMBINE_B(x) x[1] = x[2] + x[6]
#define COMBINE_C(x) x[2] = x[3] + x[7]
			nBins = 3;
			COMBINE_A(Sum), COMBINE_A(Sum2), COMBINE_A(BinN);
			COMBINE_B(Sum), COMBINE_B(Sum2), COMBINE_B(BinN);
			COMBINE_C(Sum), COMBINE_C(Sum2), COMBINE_C(BinN);
#undef COMBINE_C
#undef COMBINE_B
#undef COMBINE_A
		} break;

		//! N/8,N/8,N/4,N/2
		case 0b100: {
#define COMBINE_A(x) x[0] = x[0]
#define COMBINE_B(x) x[1] = x[1]
#define COMBINE_C(x) x[2] = x[2] + x[3]
#define COMBINE_D(x) x[3] = x[4] + x[5] + x[6] + x[7]
			nBins = 4;
			COMBINE_A(Sum), COMBINE_A(Sum2), COMBINE_A(BinN);
			COMBINE_B(Sum), COMBINE_B(Sum2), COMBINE_B(BinN);
			COMBINE_C(Sum), COMBINE_C(Sum2), COMBINE_C(BinN);
			COMBINE_D(Sum), COMBINE_D(Sum2), COMBINE_D(BinN);
#undef COMBINE_D
#undef COMBINE_C
#undef COMBINE_B
#undef COMBINE_A
		} break;

		//! N/4,N/8,N/8,N/2
		case 0b101: {
#define COMBINE_A(x) x[0] = x[0] + x[1]
#define COMBINE_B(x) x[1] = x[2]
#define COMBINE_C(x) x[2] = x[3]
#define COMBINE_D(x) x[3] = x[4] + x[5] + x[6] + x[7]
			nBins = 4;
			COMBINE_A(Sum), COMBINE_A(Sum2), COMBINE_A(BinN);
			COMBINE_B(Sum), COMBINE_B(Sum2), COMBINE_B(BinN);
			COMBINE_C(Sum), COMBINE_C(Sum2), COMBINE_C(BinN);
			COMBINE_D(Sum), COMBINE_D(Sum2), COMBINE_D(BinN);
#undef COMBINE_D
#undef COMBINE_C
#undef COMBINE_B
#undef COMBINE_A
		} break;

		//! N/2,N/8,N/8,N/4
		case 0b110: {
#define COMBINE_A(x) x[0] = x[0] + x[1] + x[2] + x[3]
#define COMBINE_B(x) x[1] = x[4]
#define COMBINE_C(x) x[2] = x[5]
#define COMBINE_D(x) x[3] = x[6] + x[7]
			nBins = 4;
			COMBINE_A(Sum), COMBINE_A(Sum2), COMBINE_A(BinN);
			COMBINE_B(Sum), COMBINE_B(Sum2), COMBINE_B(BinN);
			COMBINE_C(Sum), COMBINE_C(Sum2), COMBINE_C(BinN);
			COMBINE_D(Sum), COMBINE_D(Sum2), COMBINE_D(BinN);
#undef COMBINE_D
#undef COMBINE_C
#undef COMBINE_B
#undef COMBINE_A
		} break;

		//! N/2,N/4,N/8,N/8
		case 0b111: {
#define COMBINE_A(x) x[0] = x[0] + x[1] + x[2] + x[3]
#define COMBINE_B(x) x[1] = x[4] + x[5]
#define COMBINE_C(x) x[2] = x[6]
#define COMBINE_D(x) x[3] = x[7]
			nBins = 4;
			COMBINE_A(Sum), COMBINE_A(Sum2), COMBINE_A(BinN);
			COMBINE_B(Sum), COMBINE_B(Sum2), COMBINE_B(BinN);
			COMBINE_C(Sum), COMBINE_C(Sum2), COMBINE_C(BinN);
			COMBINE_D(Sum), COMBINE_D(Sum2), COMBINE_D(BinN);
#undef COMBINE_D
#undef COMBINE_C
#undef COMBINE_B
#undef COMBINE_A
		} break;
	}

	//! Finally take the harmonic mean of the noise amplitudes
	//! of all subblocks and combine it into the final fill level
	float Total = 0.0f;
	for(n=0;n<nBins;n++) {
		float s = Sum[n], s2 = Sum2[n];
		if(s2 == 0.0f) return 0.0f;
		Total += (SQR(BinN[n]) * s2) / (s*SQR(s));
	}
	return 0x1.6A09E6p1f * (nBins / Total); //! 0x1.6A09E6p1 = 2*Sqrt[2]
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
static void Block_Encode_EncodePass_GetHFExtParams(const float *Coef, float q, int Band, int FirstBand, int N, int WindowCtrl, int *_NoiseQ, int *_NoiseDecay) {
	//! Estimate the amplitude and decay parameters
	float Amplitude = 0.0f;
	float Decay     = 0.0f; {
		//! Find the noise amplitude by using a 4-point estimate
		//! across the entire segment to be analyzed, together
		//! with a constraint at X=0 to regularize the estimate.
		//! TODO: Improve this. There has to be a better solution
		//! than to simply do a point-wise approximation.
		int InitBand = Band - N/8u, EndBand = Band + N/8u;
		if(InitBand < FirstBand) InitBand = FirstBand;
		float x[5], y[5];
		int N0 = EndBand-InitBand;
		int N1 = N/8u;
		int N2 = N/4u - N/8u;
		int N3 = N/2u - N/4u;
		int N4 = N - N/2u;
		x[0] = 0.0f;
		x[1] = N1*0.5f;
		x[2] = N2*0.5f + N1;
		x[3] = N3*0.5f + N2 + N1;
		x[4] = N4*0.5f + N3 + N2 + N1;
		y[0] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, InitBand,      N0, WindowCtrl);
		y[1] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band,          N1, WindowCtrl);
		y[2] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band+N1,       N2, WindowCtrl);
		y[3] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band+N1+N2,    N3, WindowCtrl);
		y[4] = Block_Encode_EncodePass_GetNoiseAmplitude(Coef, Band+N1+N2+N3, N4, WindowCtrl);
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
