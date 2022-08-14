/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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

static inline void Block_Transform_CalculatePsychoacoustics_CalcFreqWeightTable(float *Dst, int BlockSize, float NyquistHz) {
	//! DCT+DST -> Pseudo-DFT
	BlockSize /= 2;

	//! Compute window for all subblock sizes in a sequential window
	//! This should improve cache locality vs a single large window.
	int n, SubBlockSize = BlockSize / ULC_MAX_BLOCK_DECIMATION_FACTOR;
	float LogFreqStep = logf(NyquistHz / SubBlockSize) + -0x1.BA18AAp2f; //! -0x1.BA18AAp2 = Log[1/1000]
	do {
		for(n=0;n<SubBlockSize;n++) {
			//! The weight function is a log-normal distribution function,
			//! which is used to change the amplitude of the analyzed
			//! coefficients. The idea is to drop the amplitude of the
			//! frequencies that the ear is most sensitive to (~1kHz)
			//! so that these cannot take part in masking calculations.
			//! This greatly improves stability of tones during transient
			//! passages, as well as with smaller BlockSize.
			float x = logf(n+0.5f) + LogFreqStep; //! Log[(n+0.5)*NyquistHz/SubBlockSize / 1000]
			*Dst++ = 1.0f - expf(-0.5f*SQR(x));
		}
	} while(LogFreqStep += -0x1.62E430p-1f, (SubBlockSize *= 2) <= BlockSize); //! -0x1.62E430p-1 = Log[0.5], ie. FreqStep *= 0.5
}
static inline void Block_Transform_CalculatePsychoacoustics(
	float *MaskingNp,
	float *BufferAmp2,
	void  *BufferTemp,
	int    BlockSize,
	const float *FreqWeightTable,
	uint32_t WindowCtrl
) {
	int n;
	float v;

	//! DCT+DST -> Pseudo-DFT
	BlockSize /= 2;

	//! Compute masking levels for each [sub-]block
	ULC_SubBlockDecimationPattern_t DecimationPattern = ULC_SubBlockDecimationPattern(WindowCtrl);
	do {
		int SubBlockSize = BlockSize >> (DecimationPattern&0x7);

		//! Find the subblock's normalization factor
		float Norm = 0.0f;
		for(n=0;n<SubBlockSize;n++) if((v = BufferAmp2[n]) > Norm) Norm = v;
		if(Norm != 0.0f) {
			//! Normalize the energy and convert to fixed-point
			//! This normalization step forces the sums to be as precise as
			//! possible without overflowing.
                        //! NOTE: Renormalization of the log values is derived as follows:
                        //!  -The bandwidth of a given line is given by:
                        //!    Bw(f) = HiRangeScale*f - LoRangeScale*f
                        //!  -The maximum number of summation terms occurs when the high
                        //!   range is at MaskEnd=N, as when we try to go higher, we do
                        //!   not have any more terms that we can add. Thus:
                        //!    HiRangeScale*f_Max = N
                        //!    f_Max = N/HiRangeScale
                        //!  -Combining the two equations, we get:
                        //!    Bw(f_Max) = HiRangeScale*f_Max - LoRangeScale*f_Max
                        //!              = HiRangeScale*(N/HiRangeScale) - LoRangeScale*(N/HiRangeScale)
                        //!              = N * (1 - LoRangeScale/HiRangeScale)
                        //!   Thus, the renormalization factor becomes 1/Bw(f_Max):
                        //!    Renormalization = 1 / (N * (1 - LoRangeScale/HiRangeScale))
			//! NOTE: Ensure that Weight[] is not zero or division by 0 may
			//! occur if the accumulated sums are all zeros.
			//! NOTE: If the maximum value will cause overflow to +inf,
			//! just use the largest possible scaling. This really should
			//! never happen in practice, but it's here just to be sure.
			uint32_t *Weight   = (uint32_t*)BufferTemp;
			uint32_t *EnergyNp = (uint32_t*)BufferAmp2;
			Norm = (Norm > 0x1.0p-96f) ? (0x1.FFFFFCp31f / Norm) : 0x1.FFFFFCp127f; //! NOTE: 2^32-eps*2 to ensure we don't overflow
			float LogScale = 0x1.14FF58p29f / SubBlockSize; //! (2^32/Log[2^32] / (1-2/3)) / N (round down)
			const float *ThisFreqWeightTable = FreqWeightTable + SubBlockSize-BlockSize/ULC_MAX_BLOCK_DECIMATION_FACTOR;
			for(n=0;n<SubBlockSize;n++) {
				v = BufferAmp2[n] * Norm;
				float ve = v*ThisFreqWeightTable[n];
				float vw = v;
				EnergyNp[n] = (ve <= 1.0f) ? 0 : (uint32_t)(logf(ve) * LogScale);
				Weight  [n] = (vw <= 1.0f) ? 1 : (uint32_t)vw;
			}
			float LogNorm     = -logf(Norm);                  //! Log[1/Norm]
			float InvLogScale = 0x1.D93040p-30f*SubBlockSize; //! Inverse (round up)

			//! Extract the masking levels for each line
			int      MaskBeg = 0, MaskEnd  = 0;
			uint64_t MaskSum = 0, MaskSumW = 0;
			for(n=0;n<SubBlockSize;n++) {
				//! Re-focus the analysis window
				int Old, New;
				const int RangeScaleFxp = 2;
				const int LoRangeScale = 3; //! Beg = (1-0.25)*Band
				const int HiRangeScale = 6; //! End = (1+0.50)*Band

				//! Remove samples that went out of focus
				//! NOTE: We skip /at most/ one sample, so don't loop.
				Old = MaskBeg >> RangeScaleFxp, MaskBeg += LoRangeScale;
				New = MaskBeg >> RangeScaleFxp;
				if(Old < New) {
					MaskSumW -= Weight[Old];
					MaskSum  -= Weight[Old] * (uint64_t)EnergyNp[Old];
				}

				//! Add samples that came into focus
				//! NOTE: We usually skip /at least/ one sample, but when we
				//! reach the end of the buffer (ie. x=xMax), we stop adding
				//! samples, so we can't go straight into a do-while loop.
				Old = MaskEnd >> RangeScaleFxp, MaskEnd += HiRangeScale;
				New = MaskEnd >> RangeScaleFxp; if(New > SubBlockSize) New = SubBlockSize;
				if(Old < New) do {
					MaskSumW += Weight[Old];
					MaskSum  += Weight[Old] * (uint64_t)EnergyNp[Old];
				} while(++Old < New);

				//! Extract level
				MaskingNp[n] = (MaskSum / MaskSumW)*InvLogScale + LogNorm;
			}
		}

		//! Move to next subblock
		MaskingNp  += SubBlockSize;
		BufferAmp2 += SubBlockSize;
	} while(DecimationPattern >>= 4);
}

/**************************************/
//! EOF
/**************************************/
