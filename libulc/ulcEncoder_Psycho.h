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

//! Ultra-stable psychoacoustics toggle
//!  0 = Weigh tones and noise equally (can be unstable in tones)
//!  1 = Weight out noise (can be muffled)
//! Default: Weight out the noise when using noise coding
//! (the idea being to synthesize it and use more bits to
//! code tone signals instead).
#define PSYCHO_ULTRASTABLE ULC_USE_NOISE_CODING

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
			//! NOTE: We re-scale the logarithm by (2^32)/Log[2^32], which
			//! we term LogScale (0x1.715476p27).
			Norm = 0x1.0p32f / Norm;
				for(n=0;n<SubBlockSize;n++) {
					float p   = BufferAmp2[n] * Norm;
					float pNp = (p > 1.0f) ? logf(p) : 0.0f;
					p   = ceilf(sqrtf(p) * 0x1.0p16f); //! Re-normalize to .32fxp
					pNp = ceilf(pNp*0x1.715476p27f);
					uint32_t ip   = (p   >= 0x1.0p32f) ? 0xFFFFFFFFu : (uint32_t)p;
					uint32_t ipNp = (pNp >= 0x1.0p32f) ? 0xFFFFFFFFu : (uint32_t)pNp;
					Energy  [n] = ip;
					EnergyNp[n] = ipNp;
				}
			Norm = logf(Norm)*0x1.555555p-2f; //! Log[LogNorm]/3

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
			int Log2MaxBandwidth = 1;
			int SumShift = 31-__builtin_clz(SubBlockSize) - Log2MaxBandwidth; //! Log2[SubBlockSize * MaxBandwidth]
			int BandBeg = 0, BandEnd = 0;
			uint64_t Sum = 0ull, SumW = 0ull;
#if PSYCHO_ULTRASTABLE
			int NoiseBeg = 0, NoiseEnd = 0;
			uint64_t NoiseSum = 0ull;
#endif
			for(n=0;n<SubBlockSize;n++) {
				//! Re-focus the main analysis window
				{
					int Old, New;
					const int RangeScaleFxp = 5;
					const int LoRangeScale = 29; //! Beg = 0.90625*Band
					const int HiRangeScale = 45; //! End = 1.40625*Band

					//! Remove samples that went out of focus
					//! NOTE: We skip /at most/ one sample, so don't loop.
					Old = BandBeg >> RangeScaleFxp, BandBeg += LoRangeScale;
					New = BandBeg >> RangeScaleFxp;
					if(Old < New) {
						SumW -= Energy[Old];
						Sum  -= Energy[Old] * (uint64_t)EnergyNp[Old] >> SumShift;
					}

					//! Add samples that came into focus
					//! NOTE: We usually skip /at least/ one sample, but when we
					//! reach the end of the buffer (ie. x=xMax), we stop adding
					//! samples, so we can't go straight into a do-while loop.
					Old = BandEnd >> RangeScaleFxp, BandEnd += HiRangeScale;
					New = BandEnd >> RangeScaleFxp; if(New > SubBlockSize) New = SubBlockSize;
					if(Old < New) do {
						SumW += Energy[Old];
						Sum  += Energy[Old] * (uint64_t)EnergyNp[Old] >> SumShift;
					} while(++Old < New);
				}
#if PSYCHO_ULTRASTABLE
				//! Re-focus the noise analysis window
				//! Same idea as above, except only summing the log values
				{
					int Old, New;
					const int RangeScaleFxp = 4;
					const int LoRangeScale = 15; //! Beg = 0.9375*Band
					const int HiRangeScale = 20; //! End = 1.2500*Band

					//! Remove samples that went out of focus
					Old = NoiseBeg >> RangeScaleFxp, NoiseBeg += LoRangeScale;
					New = NoiseBeg >> RangeScaleFxp;
					if(Old < New) {
						NoiseSum -= EnergyNp[Old];
					}

					//! Add samples that came into focus
					Old = NoiseEnd >> RangeScaleFxp, NoiseEnd += HiRangeScale;
					New = NoiseEnd >> RangeScaleFxp; if(New > SubBlockSize) New = SubBlockSize;
					if(Old < New) do {
						NoiseSum += EnergyNp[Old];
					} while(++Old < New);
				}
#endif
				//! Store the expected value for this band.
				//! This is essentially a contraharmonic mean in the log domain
				//! The overall idea is to implement this equation:
				//!  ImportanceLevel = CoefRe * CoefRe^2 / BandAbs^2
				//! Since we're working in the log domain, and the values
				//! are scale-invariant (only used for comparing):
				//!  LogImportanceLevel = Log[CoefRe^3] - Log[BandAbs^2]
				//!                     = Log[CoefRe] - Log[BandAbs^2]/3
				//! NOTE: Ideally, we would shift up as Sum<<SumShift prior
				//! to dividing, but Sum is already full-width 64bit and
				//! cannot shift up without a larger integer type.
				uint32_t w = (SumW >> SumShift) + ((SumW << SumShift) != 0); //! Ceiling division appears to be important sometimes
				uint64_t x = Sum/w;
#if PSYCHO_ULTRASTABLE
				x += NoiseSum >> SumShift >> Log2MaxBandwidth; //! NoiseSum/SubBlockSize. Not sure why it normalizes like this
#endif
				MaskingNp[n] = x*-0x1.D93040p-30f + Norm; //! 0x1.D93040p-30 = 1/LogScale / 3
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
