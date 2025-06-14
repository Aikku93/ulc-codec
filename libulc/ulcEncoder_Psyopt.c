/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2025, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <math.h>
#include <stdint.h>
/**************************************/
#include "ulcEncoder.h"
#include "ulcHelper.h"
#include "ulcEncoder_Internals.h"
/**************************************/
#if (ULC_USE_PSYCHOACOUSTICS || ULC_USE_NOISE_CODING)
/**************************************/

struct LineSum_t {
	int LineEnd;
	double Floor;
	double Peak;
	double PeakW;
};

static void InitLineSum(struct LineSum_t *LineSum) {
	LineSum->LineEnd = 0;
	LineSum->Floor = 0.0;
	LineSum->Peak  = 0.0;
	LineSum->PeakW = 0.0;
}

static void UpdateLineSum(
	const float *Src,
	const float *SrcLog,
	struct LineSum_t *LineSum,
	int LineEnd //! <- Range /exclusive/
) {
	int Line;
	double SumFloor = LineSum->Floor;
	double SumPeak  = LineSum->Peak;
	double SumPeakW = LineSum->PeakW;
	for(Line=LineSum->LineEnd;Line<LineEnd;Line++) {
		double v    = (double)Src   [Line];
		double vLog = (double)SrcLog[Line];
		SumFloor += vLog;
		SumPeak  += vLog * v;
		SumPeakW += v;
	}
	LineSum->LineEnd = LineEnd;
	LineSum->Floor   = SumFloor;
	LineSum->Peak    = SumPeak;
	LineSum->PeakW   = SumPeakW;
}

/**************************************/
#endif
/**************************************/
#if ULC_USE_PSYCHOACOUSTICS
/**************************************/

//! Calculate masking levels for each frequency line
void ULCi_CalculatePsychoacoustics(
	float *MaskingNp,
	float *BufferAmp2,
	void  *BufferTemp,
	int    BlockSize,
	int    RateHz,
	uint32_t WindowCtrl
) {
	const int ULC_N_BARK_BANDS = 25;
	float NyquistHz = (float)RateHz * 0.5f;

	//! DCT+DST -> Pseudo-DFT
	BlockSize /= 2;

	//! Compute logarithm for all lines to speed up calculations
	{
		int Line;
		for(Line=0;Line<BlockSize;Line++) {
			MaskingNp[Line] = FastLog(0x1.0p-126f + BufferAmp2[Line]);
		}
	}

	//! Compute masking levels for each [sub-]block
	ULC_SubBlockDecimationPattern_t DecimationPattern = ULCi_SubBlockDecimationPattern(WindowCtrl);
	do {
		int SubBlockSize = BlockSize >> (DecimationPattern&0x7);

		//! Iterate over all Bark bands
		//! If any band is silent, then we use the ratio of the
		//! prior band, since we will later interpolate using
		//! these values.
		//! Because the Bark bands (and thus frequency lines)
		//! are being scanned in increasing order, we optimize
		//! by accumulating a "low" and "high" sum, which will
		//! correspond to the start and end of the analysis
		//! region, and subtract low from high to get the sum.
		int BarkBand;
		float MaskRatio = 0.0f;
		float *BarkMask = (float*)BufferTemp;
		struct LineSum_t LineSumLo, LineSumHi;
		InitLineSum(&LineSumLo);
		InitLineSum(&LineSumHi);
		for(BarkBand=0;BarkBand<ULC_N_BARK_BANDS;BarkBand++) {
			//! Get the lines corresponding to this Bark band
			float FreqBeg = ULCi_BarkToFreq(BarkBand+0);
			float FreqEnd = ULCi_BarkToFreq(BarkBand+1);
			int   LineBeg = (int)floorf(ULCi_FreqToLine(FreqBeg, NyquistHz, SubBlockSize));
			int   LineEnd = (int)ceilf (ULCi_FreqToLine(FreqEnd, NyquistHz, SubBlockSize));
			if(LineBeg < 0) LineBeg = 0;
			if(LineEnd < 0) LineEnd = 0;
			if(LineBeg > SubBlockSize-1) LineBeg = SubBlockSize-1;
			if(LineEnd > SubBlockSize)   LineEnd = SubBlockSize;

			//! Update levels for this band
			UpdateLineSum(BufferAmp2, MaskingNp, &LineSumLo, LineBeg);
			UpdateLineSum(BufferAmp2, MaskingNp, &LineSumHi, LineEnd);
			double SumFloor = LineSumHi.Floor - LineSumLo.Floor;
			double SumPeak  = LineSumHi.Peak  - LineSumLo.Peak;
			double SumPeakW = LineSumHi.PeakW - LineSumLo.PeakW;

			//! Get the final masking ratio for this band
			//! This basically amounts to how much CANNOT be
			//! masked by this band, and then applies a
			//! scaling to this Bark band so as to normalize
			//! by the amount of energy present.
			if(SumPeakW > 0.0) {
				SumPeak   = SumPeak  / SumPeakW;
				SumFloor  = SumFloor / (double)(LineEnd - LineBeg);
				MaskRatio = (float)(SumFloor - SumPeak - log(SumPeakW));
			}
			BarkMask[BarkBand] = MaskRatio;
		}

		//! Now generate masking level for each frequency line
		int Line;
		for(Line=0;Line<SubBlockSize;Line++) {
			float BarkBand = ULCi_FreqToBark(ULCi_LineToFreq(Line, NyquistHz, SubBlockSize));
			int   BandIdx  = (int)BarkBand;
			float BarkFrac = BarkBand - (float)BandIdx;
			      BandIdx  = (BandIdx >=               0) ? BandIdx : 0;
			      BandIdx  = (BandIdx < ULC_N_BARK_BANDS) ? BandIdx : (ULC_N_BARK_BANDS-1);
			float BarkL = BarkMask[BandIdx];
			float BarkR = (BandIdx+1 < ULC_N_BARK_BANDS) ? BarkMask[BandIdx+1] : BarkL;
			MaskingNp[Line] = BarkL*(1.0f-BarkFrac) + BarkR*BarkFrac;
		}

		//! Move to next subblock
		MaskingNp  += SubBlockSize;
		BufferAmp2 += SubBlockSize;
	} while(DecimationPattern >>= 4);
}

/**************************************/
#endif
/**************************************/
#if ULC_USE_NOISE_CODING
/**************************************/

//! Compute noise spectrum (logarithmic output)
//! The code here is very similar to the one used in
//! psychoacoustics (see ulcEncoder_Psyopt.c for details).
//! The main difference is that we're extracting the noise
//! level after masking with the tone level, rather than
//! the other way around.
void ULCi_CalculateNoiseLogSpectrum(float *Data, void *Temp, int N, int RateHz) {
	const int ULC_N_BARK_BANDS = 25;
	float NyquistHz = (float)RateHz * 0.5f;

	//! DCT+DST -> Pseudo-DFT
	N /= 2;

	//! Compute logarithm for all lines to speed up calculations
	float *LogData = (float*)Temp; {
		int Line;
		for(Line=0;Line<N;Line++) {
			LogData[Line] = FastLog(0x1.0p-126f + Data[Line]);
		}
	}

	//! Iterate over all Bark bands
	int BarkBand;
	float *BarkMask = LogData + N;
	struct LineSum_t LineSumLo, LineSumHi;
	InitLineSum(&LineSumLo);
	InitLineSum(&LineSumHi);
	for(BarkBand=0;BarkBand<ULC_N_BARK_BANDS;BarkBand++) {
		//! Get the lines corresponding to this Bark band
		float FreqBeg = ULCi_BarkToFreq(BarkBand+0);
		float FreqEnd = ULCi_BarkToFreq(BarkBand+1);
		int   LineBeg = (int)floorf(ULCi_FreqToLine(FreqBeg, NyquistHz, N));
		int   LineEnd = (int)ceilf (ULCi_FreqToLine(FreqEnd, NyquistHz, N));
		if(LineBeg < 0) LineBeg = 0;
		if(LineEnd < 0) LineEnd = 0;
		if(LineBeg > N-1) LineBeg = N-1;
		if(LineEnd > N)   LineEnd = N;

		//! Sum levels for this band
		UpdateLineSum(Data, LogData, &LineSumLo, LineBeg);
		UpdateLineSum(Data, LogData, &LineSumHi, LineEnd);
		double SumFloor = LineSumHi.Floor - LineSumLo.Floor;
		double SumPeak  = LineSumHi.Peak  - LineSumLo.Peak;
		double SumPeakW = LineSumHi.PeakW - LineSumLo.PeakW;

		//! Get the final noise ratio for this band
		float MaskRatio = 0.0f;
		if(SumPeakW > 0.0) {
			SumPeak   = SumPeak  / SumPeakW;
			SumFloor  = SumFloor / (double)(LineEnd - LineBeg);
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
		float BarkBand = ULCi_FreqToBark(ULCi_LineToFreq(Line, NyquistHz, N));
		int   BandIdx  = (int)BarkBand;
		      BandIdx  = (BandIdx >=               0) ? BandIdx : 0;
		      BandIdx  = (BandIdx < ULC_N_BARK_BANDS) ? BandIdx : (ULC_N_BARK_BANDS-1);
		float w     = expf(BarkMask[BandIdx]);
		float Noise = LogData[Line] + BarkMask[BandIdx];
		Data[Line*2+0] = w;
		Data[Line*2+1] = w * (Noise*0.5f + 0x1.62E430p-1f); //! Pre-scale by Scale=4.0/2 for noise quantizer (by adding Log[Scale]));
	}

}

/**************************************/
#endif
/**************************************/
//! EOF
/**************************************/
