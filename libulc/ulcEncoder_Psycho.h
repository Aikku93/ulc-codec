/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2023, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <math.h>
#include <stdint.h>
/**************************************/
#include "ulcEncoder.h"
#include "ulcHelper.h"
/**************************************/

static inline void Block_Transform_CalculatePsychoacoustics(
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
			MaskingNp[Line] = logf(0x1.0p-127f + BufferAmp2[Line]);
		}
	}

	//! Compute masking levels for each [sub-]block
	ULC_SubBlockDecimationPattern_t DecimationPattern = ULC_SubBlockDecimationPattern(WindowCtrl);
	do {
		int SubBlockSize = BlockSize >> (DecimationPattern&0x7);

		//! Iterate over all Bark bands
		int BarkBand;
		float *BarkMask = (float*)BufferTemp;
		for(BarkBand=0;BarkBand<ULC_N_BARK_BANDS;BarkBand++) {
			//! Get the lines corresponding to this Bark band
			float FreqBeg = BarkToFreq(BarkBand+0);
			float FreqEnd = BarkToFreq(BarkBand+1);
			int   LineBeg = (int)floorf(FreqToLine(FreqBeg, NyquistHz, SubBlockSize));
			int   LineEnd = (int)ceilf (FreqToLine(FreqEnd, NyquistHz, SubBlockSize));
			if(LineBeg < 0) LineBeg = 0;
			if(LineEnd < 0) LineEnd = 0;
			if(LineBeg > SubBlockSize-1) LineBeg = SubBlockSize-1;
			if(LineEnd > SubBlockSize)   LineEnd = SubBlockSize;

			//! Sum levels for this band
			double SumFloor = 0.0;
			double SumPeak  = 0.0;
			double SumPeakW = 0.0;
			int nLines = LineEnd - LineBeg;
			if(nLines > 0) {
				int Line;
				const float *Src    = BufferAmp2 + LineBeg;
				const float *SrcLog = MaskingNp  + LineBeg;
				for(Line=0;Line<nLines;Line++) {
					double v    = (double)Src   [Line];
					double vLog = (double)SrcLog[Line];
					SumFloor += vLog;
					SumPeak  += vLog * v;
					SumPeakW += v;
				}
			}

			//! Get the final masking ratio for this band
			float MaskRatio = 0.0f;
			if(SumPeakW != 0.0) {
				SumPeak   = SumPeak  / SumPeakW;
				SumFloor  = SumFloor / (float)nLines;
				MaskRatio = (float)(SumPeak + SumFloor);
			}
			BarkMask[BarkBand] = MaskRatio;
		}

		//! Now generate masking level for each frequency line
		int Line;
		for(Line=0;Line<SubBlockSize;Line++) {
			float BarkBand = FreqToBark(LineToFreq(Line, NyquistHz, SubBlockSize));
			int   BandIdx  = (int)BarkBand;
			      BandIdx  = (BandIdx >=               0) ? BandIdx : 0;
			      BandIdx  = (BandIdx < ULC_N_BARK_BANDS) ? BandIdx : (ULC_N_BARK_BANDS-1);
			MaskingNp[Line] = BarkMask[BandIdx] - MaskingNp[Line];
		}

		//! Move to next subblock
		MaskingNp  += SubBlockSize;
		BufferAmp2 += SubBlockSize;
	} while(DecimationPattern >>= 4);
}

/**************************************/
//! EOF
/**************************************/
