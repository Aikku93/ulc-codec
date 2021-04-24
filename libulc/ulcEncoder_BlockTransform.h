/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
#include <stdint.h>
/**************************************/
#include "Fourier.h"
#include "ulcEncoder_Psycho.h"
#include "ulcEncoder_WindowControl.h"
#include "ulcHelper.h"
/**************************************/

//! Store the Neper-domain coefficients and sorting (importance) indices of
//! a block, and update the number of codeable non-zero coefficients
//! NOTE:
//!  AnalysisPower is used to alter the preference for the currently-being-analyzed channel
static inline void Block_Transform_WriteNepersAndIndices(
	      float *CoefIdx,
#if ULC_USE_PSYCHOACOUSTICS
	const float *MaskingNp,
#endif
	const float *Coef,
	      int   *nNzCoef,
	      int    BlockSize
) {
	int Band;
	for(Band=0;Band<BlockSize;Band++) {
		//! Coefficient inside codeable range?
		float Val = ABS(Coef[Band]);
		if(Val < 0.5f*ULC_COEF_EPS) {
			CoefIdx[Band] = -0x1.0p126f; //! Unusable coefficient; map to the end of the list
		} else {
			float ValNp = logf(Val);
			float MaskedValNp = ValNp;
#if ULC_USE_PSYCHOACOUSTICS
			//! Apply psychoacoustic corrections to this band energy
			MaskedValNp -= MaskingNp[Band];
#endif
			//! Store the sort value for this coefficient
			CoefIdx[Band] = MaskedValNp;
			(*nNzCoef)++;
		}
	}
}

/**************************************/

//! Transform a block and prepare its coefficients
static void Block_Transform_BufferInterleave(float *Buf, float *Tmp, int BlockSize, int Decimation) {
	//! The interleaving patterns here were chosen to try and
	//! optimize coefficient clustering across block sizes
	int n; for(n=0;n<BlockSize;n++) Tmp[n] = Buf[n];
	switch(Decimation >> 1) { //! Lowermost bit only controls which subblock gets overlap scaling, so ignore it
		//! 001x: a=N/2, b=N/2
		case 0b001: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/2;
			for(n=0;n<BlockSize/2;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
			}
		} break;

		//! 010x: a=N/4, b=N/4, c=N/2
		case 0b010: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/4;
			float *SrcC = SrcB + BlockSize/4;
			for(n=0;n<BlockSize/4;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcC++;
			}
		} break;

		//! 011x: a=N/2, b=N/4, c=N/4
		case 0b011: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/2;
			float *SrcC = SrcB + BlockSize/4;
			for(n=0;n<BlockSize/4;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
			}
		} break;

		//! 100x: a=N/8, b=N/8, c=N/4, d=N/2
		case 0b100: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/8;
			float *SrcC = SrcB + BlockSize/8;
			float *SrcD = SrcC + BlockSize/4;
			for(n=0;n<BlockSize/8;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
			}
		} break;

		//! 101x: a=N/4, b=N/8, c=N/8, d=N/2
		case 0b101: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/4;
			float *SrcC = SrcB + BlockSize/8;
			float *SrcD = SrcC + BlockSize/8;
			for(n=0;n<BlockSize/8;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
			}
		} break;

		//! 110x: a=N/2, b=N/8, c=N/8, d=N/4
		case 0b110: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/2;
			float *SrcC = SrcB + BlockSize/8;
			float *SrcD = SrcC + BlockSize/8;
			for(n=0;n<BlockSize/8;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcD++;
				*Buf++ = *SrcD++;
			}
		} break;

		//! 111x: a=N/2, b=N/4, c=N/8, d=N/8
		case 0b111: {
			float *SrcA = Tmp;
			float *SrcB = SrcA + BlockSize/2;
			float *SrcC = SrcB + BlockSize/4;
			float *SrcD = SrcC + BlockSize/8;
			for(n=0;n<BlockSize/8;n++) {
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcA++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcB++;
				*Buf++ = *SrcC++;
				*Buf++ = *SrcD++;
			}
		} break;
	}
}
static inline void Block_Transform_SortIndices_SiftDown(const float *SortValues, int *Order, int Root, int N) {
	//! NOTE: Most of the sorting time is spent in this function,
	//! so I've tried to optimize it as best as I could. The
	//! code below is what performed best on my Intel i7 CPU,
	//! so mileage may vary. There may exist other variations
	//! that have even better performance, but this is good
	//! enough that it probably doesn't need it anyway.
	int Child = 2*Root+1;
	if(Child < N) for(;;) {
		//! Get the values of the root and child
		//! NOTE: If a sibling exists, check it also.
		//! NOTE: Using pre-increment (combined with post-decrement
		//! for the FALSE case) appears to result in better performance?
		float RootVal  = SortValues[Order[Root]];
		float ChildVal = SortValues[Order[Child]], v;
		if(++Child < N && (v = SortValues[Order[Child]]) < ChildVal)
			ChildVal = v;
		else Child--;

		//! Already in order? Exit
		if(ChildVal > RootVal) return;

		//! Swap to restore heap order, and travel further down the heap
		//! NOTE: Stashing Order[Root],Order[Child] in temps appears to
		//! give slightly worse performance. This may be an x86 quirk
		//! due to its low register count, though.
		int t = Order[Root];
		Order[Root]  = Order[Child];
		Order[Child] = t;
		Root = Child, Child = 2*Root+1;
		if(Child >= N) return;
	}
}
static inline void Block_Transform_SortIndices(int *SortedIndices, const float *SortValues, int *Temp, int N) {
	int n;

	//! Start with mapping the indices directly
	int *Order = Temp;
	for(n=0;n<N;n++) Order[n] = n;

	//! Begin sorting
	//! NOTE: This was heavily borrowed from Rosetta Code (Heapsort).
	//! qsort() has some edge cases that degrade performance enough
	//! to roll our own variant here.
	{
		//! Heapify the array
		n = N/2u - 1;
		do Block_Transform_SortIndices_SiftDown(SortValues, Order, n, N); while(--n >= 0);

		//! Pop all elements off the heap one at a time
		n = N-1;
		do {
			int t = Order[n];
			Order[n] = Order[0];
			Order[0] = t;
			Block_Transform_SortIndices_SiftDown(SortValues, Order, 0, n);
		} while(--n);
	}

	//! Remap indices based on their sort order
	for(n=0;n<N;n++) SortedIndices[Order[n]] = n;
}
static int Block_Transform(struct ULC_EncoderState_t *State, const float *Data) {
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;

	//! Get the window control parameters for this block and the next
	int WindowCtrl     = State->WindowCtrl     = State->NextWindowCtrl;
	int NextWindowCtrl = State->NextWindowCtrl = Block_Transform_GetWindowCtrl(
		Data,
		State->SampleBuffer,
		State->TransformTemp,
		BlockSize,
		nChan
	);
	int NextBlockSubBlockSize = BlockSize; {
		//! Adjust for the first [sub]block's size
		int NextDecimation = NextWindowCtrl >> 4;
		NextBlockSubBlockSize >>= ULC_Helper_SubBlockDecimationPattern(NextDecimation)[0];
	}

	//! Transform channels and insert keys for each codeable coefficient
	//! It's not /strictly/ required to calculate nNzCoef, but it can
	//! speed things up in the rate-control step
	int nNzCoef = 0; {
		int n, Chan;
		float *BufferSamples = State->SampleBuffer;
		float *BufferMDCT    = State->TransformBuffer;
		float *BufferIndex   = (float*)State->TransformIndex;
		float *BufferFwdLap  = State->TransformFwdLap;
#if ULC_USE_NOISE_CODING
		float *BufferNoise   = State->TransformNoise;
#endif
		float *BufferTemp    = State->TransformTemp;
		float *BufferMDST    = BufferIndex;                       //! NOTE: Aliasing of BufferIndex
#if ULC_USE_PSYCHOACOUSTICS
		float *MaskingNp     = BufferIndex + (nChan-1)*BlockSize; //! NOTE: Aliasing of BufferIndex in last channel
		float *BufferAmp2    = BufferTemp + BlockSize;            //! NOTE: Using upper half of BufferTemp

		//! Clear the amplitude buffer; we'll be accumulating all channels here
		for(n=0;n<BlockSize;n++) BufferAmp2[n] = 0.0f;
#endif
		//! Apply M/S transform to the data
		//! NOTE: Fully normalized; not orthogonal.
		if(nChan == 2) for(n=0;n<BlockSize;n++) {
			float L = BufferSamples[n];
			float R = BufferSamples[n + BlockSize];
			BufferSamples[n]             = (L+R) * 0.5f;
			BufferSamples[n + BlockSize] = (L-R) * 0.5f;
		}

		//! Transform the input data and get complexity measure (ABR, VBR modes)
		int ThisTransientSubBlockIdx = ULC_Helper_TransientSubBlockIndex(WindowCtrl >> 4);
		const int8_t *DecimationPattern = ULC_Helper_SubBlockDecimationPattern(WindowCtrl >> 4);
		float Complexity = 0.0f, ComplexityW = 0.0f;
		for(Chan=0;Chan<nChan;Chan++) {
			//! Process each subblock sequentially
			const float *BufferMDCTEnd = BufferMDCT + BlockSize;
			int SubBlockIdx;
			for(SubBlockIdx=0;(BufferMDCT < BufferMDCTEnd);SubBlockIdx++) {
				//! Get the size of this subblock and its overlap
				int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];
				int OverlapSize  = SubBlockSize;
				if(SubBlockIdx == ThisTransientSubBlockIdx) OverlapSize >>= (WindowCtrl & 0x7);

				//! If the next block's first subblock is smaller than
				//! this one, set overlap to the smaller of the two.
				//! NOTE: This is literally THE ONLY reason that there is a
				//! coding delay in this codec (and most MDCT codecs).
				if(BufferMDCT+SubBlockSize < BufferMDCTEnd) {
					int NextSubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx+1];
					if(OverlapSize > NextSubBlockSize) OverlapSize = NextSubBlockSize;
				} else {
					if(OverlapSize > NextBlockSubBlockSize) OverlapSize = NextBlockSubBlockSize;
				}

				//! Cycle data through the lapping buffer
				float *SmpBuf = BufferTemp; {
					float *SmpDst = SmpBuf;

					/*!   |            . | .____________|
					      |            . |/.            |
					      |            . | .            |
					      |            . | .            |
					      |____________./| .            |
					      |            . | .            |
					      <-    L    -><-M-><-    R    ->
					      <-BlockSize/2->|<-BlockSize/2->
					      <-         BlockSize         ->

					    L_Size = (BlockSize - SubBlockSize)/2
					    M_Size = SubBlockSize
					    R_Size = (BlockSize - SubBlockSize)/2
					    L_Offs = 0
					    M_Offs = (BlockSize - SubBlockSize)/2
					    R_Offs = (BlockSize + SubBlockSize)/2

					    The L segment contains all 0s.
					    The M segment contains our lapping data.
					    The R segment contains data that we must transform.

					    Therefore:
					    When using Fourier_MDCT_MDST(), we must align BufferFwdLap
					    with the M segment, and cycle data through the R segment.
					!*/

					//! Point to the R segment and get its size
					float *LapBuf = BufferFwdLap + (BlockSize+SubBlockSize)/2;
					int LapRem = (BlockSize-SubBlockSize)/2;

					//! Do we have the full [sub]block in the lapping buffer?
					if(LapRem < SubBlockSize) {
						//! Don't have enough samples in lapping buffer
						//! for the full block - stream new data in and
						//! refill the lapping buffer
						for(n=0;n<LapRem;n++) *SmpDst++ = LapBuf[n];
						for(n=LapRem;n<SubBlockSize;n++) *SmpDst++ = *BufferSamples++;
						for(n=0;n<LapRem;n++) LapBuf[n] = *BufferSamples++;
					} else {
						//! We got a full [sub]block, and we might have data
						//! remaining in the lapping buffer. This data must
						//! now be shifted down and then the lapping buffer
						//! refilled with new data
						float *Dst = LapBuf;
						float *Src = LapBuf + SubBlockSize;
						for(n=0;n<SubBlockSize;n++) *SmpDst++ = LapBuf[n];
						for(n=SubBlockSize;n<LapRem;n++) *Dst++ = *Src++;
						for(n=0;n<SubBlockSize;n++) *Dst++ = *BufferSamples++;
					}
				}

				//! Perform the actual MDCT+MDST
				Fourier_MDCT_MDST(
					BufferMDCT,
					BufferMDST,
					SmpBuf,
					BufferFwdLap + (BlockSize-SubBlockSize)/2,
					BufferTemp,
					SubBlockSize,
					OverlapSize
				);

				//! Normalize spectrum, and accumulate amplitude by
				//! treating MDCT as Re and MDST as Im (akin to DFT)
				float Norm = 2.0f / SubBlockSize;
				for(n=0;n<SubBlockSize;n++) {
					float Re = (BufferMDCT[n] *= Norm);
					float Im = (BufferMDST[n] *= Norm);
					float Abs2 = SQR(Re) + SQR(Im);
#if ULC_USE_NOISE_CODING
					BufferTemp[n]  = Abs2;
#endif
#if ULC_USE_PSYCHOACOUSTICS
					BufferAmp2[n] += Abs2;
#endif
				}

				//! Remove transient bias for noise amplitude spectrum
				//! and add this subblock to the complexity analysis.
				//! Transfer function: H(z) = -z^-1 + 2 - z^1 (zero-phase highpass)
				//! NOTE: In theory, this has non-unity-gain of 4.0. In practice,
				//! however, the gain is 2.0 because of the non-negative inputs.
				//! NOTE: 'Transient biasing' refers to the fact that transients
				//! form a sinusoid in the Re and Im parts, but they are exactly
				//! Pi/2 radians apart. By adding Re^2+Im^2, they form a constant
				//! amplitude throughout the spectrum, which we then subtract to
				//! avoid biasing issues in noise and complexity analysis.
				//! NOTE: For noise analysis, we use a geometric mean to determine
				//! the correct noise amplitude. However, the geometric mean tends
				//! towards 1/E rather than 0.5 for white noise. So we scale by E
				//! to compensate and get unity amplitude, but then scale by 0.5
				//! to account for this filter's gain. 0x1.5BF0A9p0 = E/2.
				//! NOTE: We apply the filter over the squared samples, because
				//! this appears to give better results for some reason. The
				//! scaling then works out to have a gain of Sqrt[2], but this
				//! appears to sound better than unity gain for some reason.
				{
					float v, v2;
#if ULC_USE_NOISE_CODING
# define STORE_VALUE(n, Expr) \
	v2 = (Expr), v = sqrtf(v2), \
	Complexity += v2, ComplexityW += v, \
	BufferNoise[n] = 0x1.5BF0A9p0f * v
#else
# define STORE_VALUE(n, Expr) \
	v2 = (Expr), v = sqrtf(v2), \
	Complexity += v2, ComplexityW += v
#endif
					STORE_VALUE(0, 2.0f * ABS(BufferTemp[0] - BufferTemp[1])); //! H(z) = -z^1 + 2 - z^1 = 2 - 2z^1
					for(n=1;n<SubBlockSize-1;n++) {
						STORE_VALUE(n, ABS(-BufferTemp[n-1] + 2.0f*BufferTemp[n] - BufferTemp[n+1]));
					}
					STORE_VALUE(n, 2.0f * ABS(BufferTemp[n] - BufferTemp[n-1])); //! H(z) = -z^1 + 2 - z^1 = 2 - 2z^1
#undef STORE_VALUE
				}

				//! Move to the next subblock
				BufferMDCT  += SubBlockSize;
				BufferMDST  += SubBlockSize;
#if ULC_USE_PSYCHOACOUSTICS
				BufferAmp2  += SubBlockSize;
#endif
#if ULC_USE_NOISE_CODING
				BufferNoise += SubBlockSize;
#endif
			}

			//! Cache the sample data for the next block
			for(n=0;n<BlockSize;n++) BufferSamples[n-BlockSize] = *Data++;

			//! Move to the next channel
			BufferFwdLap += BlockSize;
#if ULC_USE_PSYCHOACOUSTICS
			BufferAmp2   -= BlockSize; //! <- Accumulated across all channels - rewind
#endif
		}
		BufferMDCT -= BlockSize*nChan; //! Rewind to start of buffer
		BufferMDST -= BlockSize*nChan;

		//! Perform encoding analysis
		{
			//! Finalize and store block complexity
			if(Complexity) {
				//! Based off the same principles of normalized entropy:
				//!  Entropy = (Log[Total[x]] - Total[x*Log[x]]/Total[x]) / Log[N]
				//! Instead of accumulating log values, we accumulate
				//! raw values, meaning we need to take the log:
				//!  Total[x*Log[x]]/Total[x] -> Log[Total[x*x]/Total[x]]
				//! Simplifying:
				//!   (Log[Total[x]] - Log[Total[x*x]/Total[x]]) / Log[N]
				//!  =(Log[Total[x] / (Total[x*x]/Total[x])) / Log[N]
				//!  =Log[Total[x]^2 / Total[x^2]] / Log[N]
				float ComplexityScale = 0x1.62E430p-1f*(31 - __builtin_clz(BlockSize)); //! 0x1.62E430p-1 = 1/Log2[E] for change-of-base
				Complexity = logf(SQR(ComplexityW) / Complexity) / ComplexityScale;
				if(Complexity < 0.0f) Complexity = 0.0f; //! In case of round-off error
				if(Complexity > 1.0f) Complexity = 1.0f;
			}
			State->BlockComplexity = Complexity;

			//! Perform psychoacoustic analysis
#if ULC_USE_PSYCHOACOUSTICS
			{
				//! Combine the energy of all channels' MDCT+MDST coefficients
				//! into normalized (by maximum) fixed-point integer values.
				//! For some reason, it seems to work better to combine the
				//! channels for analysis, rather than use each one separately.
				uint32_t *BufferEnergy   = (uint32_t*)(BufferTemp);
				uint32_t *BufferEnergyNp = (uint32_t*)(BufferTemp + BlockSize);

				//! Find the maximum value for normalization
				//! NOTE: Ideally, we'd take the maximum of each subblock, but
				//! this should work well enough for our purposes, as all
				//! coefficients are "normalized" to begin with, and this is
				//! just a slight nudge to improve integer precision.
				//! NOTE: Using M/S makes no difference to using L/R; the
				//! equations cancel out the inner terms, leaving 2L^2+2R^2.
				float v, Norm = 0.0f;
				for(n=0;n<BlockSize;n++) if((v = BufferAmp2[n]) > Norm) Norm = v;

				//! Find the integer normalization factor and convert to integers
				//! NOTE: Maximum value for BufferEnergyNp is Log[2 * 2^32], so rescale by
				//! (2^32)/Log[2 * 2^32] (and clip to prevent overflow issues on some CPUs).
				if(Norm != 0.0f) {
					const float LogScale = 0x1.66235Bp27f; //! (2^32) / Log[2 * 2^32]
					Norm = 0x1.0p32f / Norm;
					for(n=0;n<BlockSize;n++) {
						float p   = BufferAmp2[n] * Norm;
						float pNp = (p < 0.5f) ? 0.0f : logf(2.0f*p); //! Scale by 2 to keep values strictly non-negative
						p   = ceilf(p);
						pNp = ceilf(pNp*LogScale);
						uint32_t ip   = (p   >= 0x1.0p32f) ? 0xFFFFFFFFu : (uint32_t)p;
						uint32_t ipNp = (pNp >= 0x1.0p32f) ? 0xFFFFFFFFu : (uint32_t)pNp;
						BufferEnergy  [n] = ip;
						BufferEnergyNp[n] = ipNp;
					}

					//! Buffer the psychoacoustics analysis into the last channel's buffer
					//! so that we can save a bit of memory to store these values into and
					//! avoid re-calculating for each channel, as the analysis isn't cheap
					Block_Transform_CalculatePsychoacoustics(MaskingNp, BufferEnergy, BufferEnergyNp, BlockSize, DecimationPattern);
				}
			}
#endif
			//! Analyze each channel
			for(Chan=0;Chan<nChan;Chan++) {
				//! Analyze each subblock separately
				const float *BufferMDCTEnd = BufferMDCT + BlockSize;
				int SubBlockIdx;
				for(SubBlockIdx=0;(BufferMDCT < BufferMDCTEnd);SubBlockIdx++) {
					int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];

					//! Store the sorting (importance) indices
					Block_Transform_WriteNepersAndIndices(
						BufferIndex,
#if ULC_USE_PSYCHOACOUSTICS
						MaskingNp,
#endif
						BufferMDCT,
						&nNzCoef,
						SubBlockSize
					);

					//! Move to the next subblock
					BufferMDCT  += SubBlockSize;
					BufferIndex += SubBlockSize;
#if ULC_USE_PSYCHOACOUSTICS
					MaskingNp   += SubBlockSize;
#endif
				}
#if ULC_USE_PSYCHOACOUSTICS
				//! Psychoacoustic analysis is re-used across channels - rewind
				MaskingNp -= BlockSize;
#endif
			}
		}

		//! Interleave the transform data for coding
		if(WindowCtrl & 0x8) {
			int Decimation = WindowCtrl >> 4;
			BufferMDCT  = State->TransformBuffer;
#if ULC_USE_NOISE_CODING
			BufferNoise = State->TransformNoise;
#endif
			BufferIndex = (float*)State->TransformIndex;
			for(Chan=0;Chan<nChan;Chan++) {
#define INTERLEAVE(Buf) Block_Transform_BufferInterleave(Buf, BufferTemp, BlockSize, Decimation); Buf += BlockSize
				INTERLEAVE(BufferMDCT);
#if ULC_USE_NOISE_CODING
				INTERLEAVE(BufferNoise);
#endif
				INTERLEAVE(BufferIndex);
#undef INTERLEAVE
			}
		}
	}

	//! Create the coefficient sorting indices
	{
		int *BufferTmp = (int*)State->TransformTemp;
		int *BufferIdx = State->TransformIndex;
		Block_Transform_SortIndices(BufferIdx, (float*)BufferIdx, BufferTmp, nChan * BlockSize);
	}
	return nNzCoef;
}

/**************************************/
//! EOF
/**************************************/
