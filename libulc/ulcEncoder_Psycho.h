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

static inline void Block_Transform_CalculatePsychoacoustics(float *MaskingNp, const float *BufferAmp2, uint32_t *BufferTemp, int BlockSize, uint32_t WindowCtrl) {
	int n;
	float v;

	//! Compute masking levels for each [sub-]block
	uint32_t *Energy   = (uint32_t*)(BufferTemp);
	uint32_t *EnergyNp = (uint32_t*)(BufferTemp + BlockSize);
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
			//! NOTE: The normalization is based on the widest bandwidth we
			//! will encounter in the loop; the tone analysis is bound by
			//! 1-LoScale/HiScale (see below for derivation), and the noise
			//! analysis is bound by (HiScale-LowScale), so we use whichever
			//! of these is wider.
			//! NOTE: Ensure that Energy[] is not zero or division by 0 may
			//! occur if the accumulated sums are all zeros.
			//! NOTE: Do NOT remove the square root in computing Energy[].
			//! Although the difference is marginal, it works better.
			//! NOTE: If the maximum value will cause overflow to +inf,
			//! just use the largest possible scaling. This really should
			//! never happen in practice, but it's here just to be sure.
			Norm = (Norm > 0x1.0p-96f) ? (0x1.FFFFFCp31f / Norm) : 0x1.FFFFFCp127f; //! NOTE: 2^32-eps*2 to ensure we don't overflow
			float LogScale = 0x1.FBD422p28f / SubBlockSize; //! (2^32/Log[2^32]) / (N * (1-14/22)) = (2^32/Log[2^32] / (1-14/22)) / N (round down)
			for(n=0;n<SubBlockSize;n++) {
				v = BufferAmp2[n] * Norm;
				EnergyNp[n] = (v <= 1.0f) ? 0 : (uint32_t)(logf(v) * LogScale);
				v = 0x1.0p16f * sqrtf(v);
				Energy  [n] = (v <= 1.0f) ? 1 : (uint32_t)v;
			}
			float LogNorm     = 0x1.555555p-2f*logf(Norm); //! Log[Norm]/3
			float InvLogScale = SubBlockSize * -0x1.582318p-31f; //! Inverse, scaled by -1/3 (round up)

			//! Compute expected level of each band's critical bandwidth
			//! NOTE: We can solve for the maximum bandwidth used in
			//! practice (given the limited range of the block size) by
			//! finding the intersection at yMax=SubBlockSize for xMax:
			//!  yMax = SubBlockSize = xMax*HiRangeScale
			//!  xMax = SubBlockSize/HiRangeScale
			//! Then we plug xMax into the bandwidth:
			//!  yBw = xMax*HiRangeScale - xMax*LoRangeScale
			//!      = SubBlockSize * (1 - LoRangeScale/HiRangeScale)
			//! Setting SubBlockSize=1 gives us the normalized bandwidth:
			//!  MaxBandwidth = 1 - LoRangeScale/HiRangeScale
			int ToneBeg  = 0, ToneEnd  = 0; uint64_t ToneSum  = 0, ToneSumW = 0;
			int NoiseBeg = 0, NoiseEnd = 0; uint32_t NoiseSum = 0;
			int NoiseNormShift = 31 - __builtin_clz(SubBlockSize); //! 1/SubBlockSize. Not sure why it normalizes like this
			for(n=0;n<SubBlockSize;n++) {
				//! Re-focus the tone analysis window
				{
					int Old, New;
					const int RangeScaleFxp = 4;
					const int LoRangeScale = 14; //! Beg = 0.875*Band
					const int HiRangeScale = 22; //! End = 1.375*Band

					//! Remove samples that went out of focus
					//! NOTE: We skip /at most/ one sample, so don't loop.
					Old = ToneBeg >> RangeScaleFxp, ToneBeg += LoRangeScale;
					New = ToneBeg >> RangeScaleFxp;
					if(Old < New) {
						ToneSumW -= Energy[Old];
						ToneSum  -= Energy[Old] * (uint64_t)EnergyNp[Old];
					}

					//! Add samples that came into focus
					//! NOTE: We usually skip /at least/ one sample, but when we
					//! reach the end of the buffer (ie. x=xMax), we stop adding
					//! samples, so we can't go straight into a do-while loop.
					Old = ToneEnd >> RangeScaleFxp, ToneEnd += HiRangeScale;
					New = ToneEnd >> RangeScaleFxp; if(New > SubBlockSize) New = SubBlockSize;
					if(Old < New) do {
						ToneSumW += Energy[Old];
						ToneSum  += Energy[Old] * (uint64_t)EnergyNp[Old];
					} while(++Old < New);
				}

				//! Re-focus the noise analysis window
				//! Same idea as above, except only summing the log values
				//! NOTE: These values were chosen so that they are exactly
				//! half the bandwidth of the 'main' analysis window. This
				//! appears to give the best results.
				//! NOTE: The range will NOT be bound by the rules set out
				//! earlier, as we re-scale our sum when we hit the end.
				//! The maximum bandwidth is thus HiScale-LoScale.
				{
					int Old, New;
					const int RangeScaleFxp = 4;
					const int LoRangeScale = 15; //! Beg = 0.9375*Band
					const int HiRangeScale = 19; //! End = 1.1875*Band

					//! Remove samples that went out of focus
					Old = NoiseBeg >> RangeScaleFxp, NoiseBeg += LoRangeScale;
					New = NoiseBeg >> RangeScaleFxp;
					if(Old < New) {
						NoiseSum -= EnergyNp[Old];
					}

					//! Add samples that came into focus
					int NewBeg = New;
					Old = NoiseEnd >> RangeScaleFxp, NoiseEnd += HiRangeScale;
					New = NoiseEnd >> RangeScaleFxp;
					if(New > SubBlockSize) {
						//! When we've gone out of range, re-scale the sum
						//! NOTE: (NewEnd-NewBeg) / (OldEnd-NewBeg).
						//! I haven't gone crazy. We use NewBeg in the divisor
						//! because we've already removed samples from OldBeg,
						//! so we must use NewBeg as a reference point.
						NoiseSum = NoiseSum * (uint64_t)(New-NewBeg) / (Old-NewBeg);
					} else do {
						NoiseSum += EnergyNp[Old];
					} while(++Old < New);
				}

				//! Store the expected value for this band.
				//! This is essentially a contraharmonic mean in the log domain
				//! The overall idea is to implement this equation:
				//!  ImportanceLevel = CoefAmplitude * CoefPower/CritBandPower
				//! Since we're working in the log domain, and the values
				//! are scale-invariant (only used for comparing):
				//!  LogImportanceLevel = 3*Log[CoefAmplitude] - Log[CritBandPower]
				//!                     = Log[CoefAmplitude] - Log[CritBandPower]/3
				//! NOTE: Cast the tone sum to 32bit, on the assumption that
				//! the compiler won't do this prior to adding the noise sum.
				uint32_t x = (uint32_t)(ToneSum / ToneSumW) + (NoiseSum >> NoiseNormShift);
				MaskingNp[n] = (float)x*InvLogScale + LogNorm;
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
