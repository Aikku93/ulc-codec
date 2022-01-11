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

static inline void Block_Transform_CalculatePsychoacoustics(float *MaskingNp, float *BufferAmp2, void *BufferTemp, int BlockSize, uint32_t WindowCtrl) {
	int n;
	float v;

	//! Compute masking levels for each [sub-]block
	ULC_SubBlockDecimationPattern_t DecimationPattern = ULC_SubBlockDecimationPattern(WindowCtrl);
	do {
		int SubBlockSize = BlockSize >> (DecimationPattern&0x7);

		//! Find the subblock's normalization factor
		float Norm = 0.0f;
		for(n=0;n<SubBlockSize;n++) if((v = BufferAmp2[n]) > Norm) Norm = v;
		if(Norm != 0.0f) {
			int Log2SubBlockSize = 31 - __builtin_clz(SubBlockSize);

			//! Normalize the energy and convert to fixed-point
			//! This normalization step forces the sums to be as precise as
			//! possible without overflowing.
			//! NOTE: The normalization of the log values is based on the
			//! final masking equation, where two values are interpolated:
			//!  (NoiseMask*(N-n) + ToneMask*n)/N
			//! 1/N is then factored out, leaving us with this sum:
			//!  NoiseMask*(N-n) + ToneMask*n
			//! which has a maxima of Max(Max(NoiseMask), Max(ToneMask))*N.
			//! NoiseMask and ToneMask are constrained to have a maximum
			//! value of 2^32-1, meaning that the normalization that will
			//! avoid overflow but still have maximum precision is 1/N.
			//! NOTE: Ensure that Energy[] is not zero or division by 0 may
			//! occur if the accumulated sums are all zeros.
			//! NOTE: If the maximum value will cause overflow to +inf,
			//! just use the largest possible scaling. This really should
			//! never happen in practice, but it's here just to be sure.
			//! NOTE: The output into MaskingNp[] is scaled by 0.5 so that
			//! we can just subtract this from the log coefficient amplitude
			//! to get the final level. The reason being that the importance
			//! level is guided as:
			//!  Importance = CoefAmplitude + CoefAmplitude-Masking
			//!             = 2*CoefAmplitude - Masking
			//! and because we are scale invariant (only need to sort):
			//!  Importance = CoefAmplitude - 0.5*Masking
			uint32_t *Energy   = (uint32_t*)BufferTemp;
			uint32_t *EnergyNp = (uint32_t*)BufferAmp2;
			Norm = (Norm > 0x1.0p-96f) ? (0x1.FFFFFCp31f / Norm) : 0x1.FFFFFCp127f; //! NOTE: 2^32-eps*2 to ensure we don't overflow
			float LogScale = 0x1.715476p27f / SubBlockSize; //! (2^32/Log[2^32]) / N (round down)
			for(n=0;n<SubBlockSize;n++) {
				v = BufferAmp2[n] * Norm;
				EnergyNp[n] = (v <= 1.0f) ? 0 : (uint32_t)(logf(v) * LogScale);
				Energy  [n] = (v <= 1.0f) ? 1 : (uint32_t)v;
			}
			float LogNorm     = -0.25f*logf(Norm); //! Log[1/Norm] * 1/2 * 1/2 (for converting Power to Amplitude plus weight for coefficient amplitude)
			float InvLogScale = 0x1.62E430p-30f;   //! Inverse, scaled by (1/2 * 1/2)/SubBlockSize (round up)

			//! Extract the tone masking levels (spectral peaks) and
			//! make the level proportional to the masking bandwidth,
			//! then extract the noise masking level (geometric mean)
			//! /without/ making these levels proportional, and then
			//! add them together for the final output.
			int ToneBeg  = 0, ToneEnd  = 0; uint64_t ToneSum  = 0, ToneSumW = 0;
			int NoiseBeg = 0, NoiseEnd = 0; uint32_t NoiseSum = 0;
			for(n=0;n<SubBlockSize;n++) {
				//! Re-focus the tone analysis window
				uint32_t ToneMask; {
					int Old, New;
					const int RangeScaleFxp = 4;
					const int LoRangeScale = 14; //! Beg = (1-0.125)*Band
					const int HiRangeScale = 22; //! End = (1+0.375)*Band

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

					//! Extract level
					ToneMask = ToneSum / ToneSumW;
				}

				//! Re-focus the noise analysis window
				uint32_t NoiseMask; {
					int Old, New;
					uint32_t MaskBw;
					const int RangeScaleFxp = 4;
					const int LoRangeScale = 15; //! Beg = (1-0.0625)*Band
					const int HiRangeScale = 21; //! End = (1+0.3125)*Band

					//! Remove samples that went out of focus
					Old = NoiseBeg >> RangeScaleFxp, NoiseBeg += LoRangeScale;
					New = NoiseBeg >> RangeScaleFxp;
					if(Old < New) {
						NoiseSum -= EnergyNp[Old];
					}
					MaskBw = New;

					//! Add samples that came into focus
					//! NOTE: Clip New /before/ updating MaskBw or the
					//! geometric mean will turn out wrong by addition
					//! of implicit 0 padding
					Old = NoiseEnd >> RangeScaleFxp, NoiseEnd += HiRangeScale;
					New = NoiseEnd >> RangeScaleFxp; if(New > SubBlockSize) New = SubBlockSize;
					if(Old < New) do {
						NoiseSum += EnergyNp[Old];
					} while(++Old < New);
					MaskBw = New - MaskBw;

					//! Extract level
					NoiseMask = NoiseSum / MaskBw;
				}

				//! Final masking level calculation:
				//!  NoiseMask + (ToneMask-NoiseMask)*n/N
				//!  =(NoiseMask*N + (ToneMask-NoiseMask)*n) / N
				//! Then factor out 1/N and combine into InvLogScale.
				//! The idea is that lower frequencies will rely more
				//! on the noise masking level (which is computed to
				//! be close to the noise floor), while frequencies
				//! that are much higher in the spectrum will rely on
				//! the pseudo-peak "tone" masking level. This tends
				//! to stabilize low/mid-freq tones while preserving
				//! brightness in the higher frequencies.
				//! NOTE: This equation can be refactored as:
				//!  (NoiseMask*(N-n) + ToneMask*n)/N
				//! Meaning that as long as we don't drop bits, the output
				//! is always UNSIGNED and we don't need to do anything
				//! special, as long as wraparound is observed in the sum;
				//! this form just trades a multiply for a shift.
				MaskingNp[n] = ((NoiseMask<<Log2SubBlockSize) + (ToneMask-NoiseMask)*n)*InvLogScale + LogNorm;
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
