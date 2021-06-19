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
		float LogNorm = 0.0f;
		for(n=0;n<SubBlockSize;n++) if((v = BufferAmp2[n]) > LogNorm) LogNorm = v;
		if(LogNorm != 0.0f) {
			//! Normalize the energy and convert to fixed-point
			//! NOTE: We re-scale the logarithm by (2^32)/Log[2^32], which
			//! we term LogScale (0x1.715476p27).
			//! NOTE: Because we normalized the coefficients here, we must
			//! undo the normalization upon storing to EnergyNp[]. This adds
			//! a small amount of extra work, but should improve precision.
			//! Note the scaling; this is to match the scaling of EnergyNp[].
			LogNorm = 0x1.0p32f / LogNorm;
				for(n=0;n<SubBlockSize;n++) {
					float p   = BufferAmp2[n] * LogNorm;
					float pNp = (p > 1.0f) ? logf(p) : 0.0f;
					p   = ceilf(sqrtf(p) * 0x1.0p16f); //! Re-normalize to .32fxp
					pNp = ceilf(pNp*0x1.715476p27f);
					uint32_t ip   = (p   >= 0x1.0p32f) ? 0xFFFFFFFFu : (uint32_t)p;
					uint32_t ipNp = (pNp >= 0x1.0p32f) ? 0xFFFFFFFFu : (uint32_t)pNp;
					Energy  [n] = ip;
					EnergyNp[n] = ipNp;
				}
			LogNorm = logf(LogNorm)*0x1.555555p-2f; //! Log[LogNorm]/3

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
			int SumShift = 31-__builtin_clz(SubBlockSize) - 2; //! Log2[SubBlockSize * MaxBandwidth]
			int BandBeg = 0, BandEnd = 0;
			uint64_t Sum = 0ull, SumW = 0ull;
			for(n=0;n<SubBlockSize;n++) {
				//! Re-focus analysis window
				const int RangeScaleFxp = 7;
				{
					int Old, New;
					const int LoRangeScale = 121; //! Beg = 0.9453125*Band
					const int HiRangeScale = 161; //! End = 1.2578125*Band

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
				uint32_t w = SumW >> SumShift; w += (w == 0); //! Avoid division by 0
				MaskingNp[n] = (Sum/w)*-0x1.D93040p-30f + LogNorm; //! 0x1.D93040p-30 = 1/LogScale / 3
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
