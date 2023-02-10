/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2022, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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
			//! which is used to basically apply a perceptual-amplitude
			//! re-normalization, by protecting sensitive bands <1kHz.
			//! This greatly improves stability of tones during transient
			//! passages, as well as with smaller BlockSize.
			//! NOTE: We "protect" everything below 1kHz by forcing the
			//! masking calculations to rely only on the floor level.
			float x = logf(n+0.5f) + LogFreqStep; //! Log[(n+0.5)*NyquistHz/SubBlockSize / 1000]
			*Dst++ = (x > 0.0f) ? expf(-2.0f*SQR(x)) : 1.0f; //! If below 1kHz (x < 0), clip (E^0 == 1.0)
		}
	} while(LogFreqStep += -0x1.62E430p-1f, (SubBlockSize *= 2) <= BlockSize); //! -0x1.62E430p-1 = Log[0.5], ie. FreqStep *= 0.5
}
static inline void Block_Transform_CalculatePsychoacoustics(
	float *MaskingNp,
	float *BufferAmp2,
	void  *BufferTemp,
	int    BlockSize,
	int    RateHz,
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
			//! Get the window bandwidth scaling constants
			int RangeScaleFxp = 16;
			int LoRangeScale; {
				float s = (2*20000.0f) / RateHz;
				if(s >= 1.0f) s = 0x1.FFFFFEp-1f; //! <- Ensure this is always < 1.0
				LoRangeScale = (int)floorf((1<<RangeScaleFxp) * s);
			}
			int HiRangeScale; {
				float s = RateHz / (2*15000.0f);
				if(s < 1.0f) s = 1.0f; //! <- Ensure this is always >= 1.0
				HiRangeScale = (int)ceilf((1<<RangeScaleFxp) * s);
			}

			//! Normalize the energy and convert to fixed-point
			//! This normalization step forces the sums to be as precise as
			//! possible without overflowing.
                        //! NOTE: We can't renormalize by taking into account the number
                        //! of summation terms, as this can vary depending on RateHz. In
                        //! theory, it's possible, but would require adjusting rounding
                        //! behaviour so that 1-Lo/Hi rounds up, which makes it nasty.
			//! NOTE: Ensure that Weight[] is not zero or division by 0 may
			//! occur if the accumulated sums are all zeros.
			//! NOTE: If the maximum value will cause overflow to +inf,
			//! just use the largest possible scaling. This really should
			//! never happen in practice, but it's here just to be sure.
			uint32_t *Weight   = (uint32_t*)BufferTemp;
			uint32_t *EnergyNp = (uint32_t*)BufferAmp2;
			Norm = (Norm > 0x1.0p-96f) ? (0x1.FFFFFCp31f / Norm) : 0x1.FFFFFCp127f; //! NOTE: 2^32-eps*2 to ensure we don't overflow
			float LogScale = 0x1.715476p27f / SubBlockSize; //! 2^32/Log[2^32] / N
			const float *ThisFreqWeightTable = FreqWeightTable + SubBlockSize-BlockSize/ULC_MAX_BLOCK_DECIMATION_FACTOR;
			for(n=0;n<SubBlockSize;n++) {
				v = BufferAmp2[n] * Norm;
				float ve = v;
				float vw = 0x1.0p16f * sqrtf(v);
				      vw = vw + (0x1.FFFFFCp31f-vw)*ThisFreqWeightTable[n];
				EnergyNp[n] = (ve <= 1.0f) ? 0 : (uint32_t)(logf(ve) * LogScale);
				Weight  [n] = (vw <= 1.0f) ? 1 : (uint32_t)vw;
			}
			float LogNorm     = -logf(Norm); //! Log[1/Norm]
			float InvLogScale = 0x1.62E430p-28f*SubBlockSize; //! Inverse (round up)

			//! Extract the masking levels for each line
			int      MaskBeg = 0, MaskEnd  = 0;
			uint64_t MaskSum = 0, MaskSumW = 0;
			uint32_t FloorSum = 0;
			for(n=0;n<SubBlockSize;n++) {
				int Old, New;

				//! Remove samples that went out of focus
				//! NOTE: We skip /at most/ one sample, so don't loop.
				Old = MaskBeg >> RangeScaleFxp, MaskBeg += LoRangeScale;
				New = MaskBeg >> RangeScaleFxp;
				if(Old < New) {
					MaskSumW -= Weight[Old];
					MaskSum  -= Weight[Old] * (uint64_t)EnergyNp[Old];
					FloorSum -= EnergyNp[Old];
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
					FloorSum += EnergyNp[Old];
				} while(++Old < New);

				//! Extract level
				//! NOTE: FloorBw's computation must be done exactly as written;
				//! attempting to simplify the maths will go wrong. Additionally,
				//! using the actual number of bands in the sum of FloorSum can
				//! cause strange artifacts at higher frequencies, so we use the
				//! theoretical number of bands that would be in that bandwidth.
				int64_t Mask  = MaskSum / MaskSumW;
				int64_t Floor = ((int64_t)FloorSum << RangeScaleFxp) / (MaskEnd - MaskBeg);
				MaskingNp[n] = (2*Mask - Floor)*InvLogScale + LogNorm;
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
