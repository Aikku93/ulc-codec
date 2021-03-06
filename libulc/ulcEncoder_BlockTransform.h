/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
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
	      float *CoefNp,
	      int   *nNzCoef,
	      int    BlockSize
) {
	int Band;
	for(Band=0;Band<BlockSize;Band++) {
		float Val = Coef[Band];

		//! Coefficient inside codeable range?
		if(ABS(Val) < 0.5f*ULC_COEF_EPS) {
			CoefNp [Band] = ULC_COEF_NEPER_OUT_OF_RANGE;
			CoefIdx[Band] = -0x1.0p126f; //! Unusable coefficient; map to the end of the list
		} else {
			float ValNp = logf(ABS(Val));
			float MaskedValNp = ValNp;
#if ULC_USE_PSYCHOACOUSTICS
			//! Apply psychoacoustic corrections to this band energy
			MaskedValNp -= MaskingNp[Band];
#endif
			//! Store the sort value for this coefficient
			CoefNp [Band] = ValNp;
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
		State->TransientEnergy,
		State->TransformTemp,
		&State->TransientCompressorGain,
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
		float *BufferNepers  = State->TransformNepers;
		float *BufferIndex   = (float*)State->TransformIndex;
		float *BufferFwdLap  = State->TransformFwdLap;
		float *BufferTemp    = State->TransformTemp;
		float *BufferMDST    = BufferNepers;                       //! NOTE: Aliasing of BufferNepers
#if ULC_USE_PSYCHOACOUSTICS
		float *MaskingNp     = BufferNepers + (nChan-1)*BlockSize; //! NOTE: Aliasing of BufferNepers in last channel
#endif
		//! Apply M/S transform to the data
		//! NOTE: Fully normalized; not orthogonal.
		if(nChan == 2) for(n=0;n<BlockSize;n++) {
			float L = BufferSamples[n];
			float R = BufferSamples[n + BlockSize];
			BufferSamples[n]             = (L+R) * 0.5f;
			BufferSamples[n + BlockSize] = (L-R) * 0.5f;
		}

		//! Transform the input data
		int ThisTransientSubBlockIdx = ULC_Helper_TransientSubBlockIndex(WindowCtrl >> 4);
		const int8_t *DecimationPattern = ULC_Helper_SubBlockDecimationPattern(WindowCtrl >> 4);
		for(Chan=0;Chan<nChan;Chan++) {
			//! Process each subblock sequentially
			const float *BufferMDCTEnd = BufferMDCT + BlockSize;
			int SubBlockIdx;
			for(SubBlockIdx=0;(BufferMDCT < BufferMDCTEnd);SubBlockIdx++) {
				//! Get the size of this subblock and its overlap
				int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];
				int OverlapSize = BlockSize >> DecimationPattern[SubBlockIdx];
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

				//! If we have a long block, all the data is prepared and
				//! we can move straight on to the transform
				const float *SmpBuf = BufferSamples;
				if(SubBlockSize != BlockSize) {
					//! Use a scratch buffer for sample data
					SmpBuf = BufferTemp;

					//! With decimation, we cycle data through the lapping
					//! buffer, so that it will be ready for future calls
					int LapExtra = (BlockSize - SubBlockSize) / 2;
					float *LapBuf = BufferFwdLap + LapExtra;

					//! Get samples from the lapping buffer
					int LapCopy = (LapExtra < SubBlockSize) ? LapExtra : SubBlockSize;
					for(n=0;n<LapCopy;n++)   BufferTemp[n] = LapBuf[-1-n];
					for(;n<SubBlockSize;n++) BufferTemp[n] = *BufferSamples++;

					//! Shuffle new samples into it
					int NewCopy = LapExtra - LapCopy;
					for(n=0;n<NewCopy;n++) LapBuf[-1-n] = LapBuf[-1-n-LapCopy];
					for(;n<LapExtra;n++) LapBuf[-1-n] = *BufferSamples++;
				} else BufferSamples += SubBlockSize;

				//! Perform the actual MDCT
				Fourier_MDCT(
					BufferMDCT,
					SmpBuf,
					BufferFwdLap + (BlockSize-SubBlockSize)/2,
					BufferTemp,
					SubBlockSize,
					OverlapSize,
					BufferMDST
				);

				//! Normalize the MDCT (and MDST) coefficients
				//! This simplifies analysis later on and is 'more correct'
				//! in the first place, as the output of the MDCT itself
				//! should've been normalized already.
				//! PONDER: Modify Fourier_MDCT() to do this automatically?
				float Norm = 2.0f / SubBlockSize;
				for(n=0;n<SubBlockSize;n++) {
					BufferMDCT[n] *= Norm;
					BufferMDST[n] *= Norm;
				}

				//! Move to the next subblock
				BufferMDCT += SubBlockSize;
				BufferMDST += SubBlockSize;
			}

			//! Cache the sample data for the next block
			for(n=0;n<BlockSize;n++) BufferSamples[n-BlockSize] = *Data++;

			//! Move to the next channel
			BufferFwdLap += BlockSize/2u;
		}
		BufferMDCT -= BlockSize*nChan; //! Rewind to start of buffer
		BufferMDST -= BlockSize*nChan;

		//! Perform importance analysis
		//! NOTE: Without psychoacoustics, we still perform some
		//! analysis for ABR coding modes (complexity analysis).
		{
			//! Combine the energy of all channels' MDCT+MDST coefficients
			//! into normalized (by maximum) fixed-point integer values.
			//! For some reason, it seems to work better to combine the
			//! channels for analysis, rather than use each one separately.
#if ULC_USE_PSYCHOACOUSTICS
			uint32_t *BufferEnergy   = (uint32_t*)(BufferTemp);
			uint32_t *BufferEnergyNp = (uint32_t*)(BufferTemp + BlockSize);
#endif
			{
				float *fBufferEnergy = BufferTemp;
				for(n=0;n<BlockSize;n++) fBufferEnergy[n] = 0.0f;

				//! Sum energy and get the maximum value
				//! NOTE: Ideally, we'd take the maximum of each subblock, but
				//! this should work well enough for our purposes, as all
				//! coefficients are "normalized" to begin with, and this is
				//! just a slight nudge to improve integer precision.
				//! NOTE: Using M/S makes no difference to using L/R; the
				//! equations cancel out the inner terms, leaving 2L^2+2R^2.
				//! NOTE: Using only the MDST coefficients here appears to
				//! improve results. Intuitively, this makes some sense (in
				//! particular: transients in MDST are rotated by Pi/2 radians
				//! and so when the final masked levels are computed, the
				//! Re^2+Im^2 cancels out to 1.0), but this probably bears more
				//! investigation before this comment can be removed/corrected.
				float Norm = 0.0f;
				const float *MDSTSrc = BufferMDST;
				for(Chan=0;Chan<nChan;Chan++) for(n=0;n<BlockSize;n++) {
					float Im = *MDSTSrc++;
					float v = (fBufferEnergy[n] += SQR(Im));
					if(v > Norm) Norm = v;
				}

				//! Find the integer normalization factor and convert
				//! NOTE: Maximum value for BufferEnergyNp is Log[2 * 2^32], so rescale by
				//! (2^32)/Log[2 * 2^32] (and clip to prevent overflow issues on some CPUs).
#if ULC_USE_PSYCHOACOUSTICS
				const float LogScale = 0x1.66235Bp27f; //! (2^32) / Log[2 * 2^32]
#endif
				if(Norm != 0.0f) Norm = 0x1.0p32f / Norm;
				float Complexity = 0.0f, ComplexityW = 0.0f;
				int ComplexityScale = 31 - __builtin_clz(BlockSize);
				for(n=0;n<BlockSize;n++) {
					float p   = fBufferEnergy[n] * Norm;
					float pNp = (p < 0.5f) ? 0.0f : logf(2.0f*p); //! Scale by 2 to keep values strictly non-negative
					Complexity  += p * pNp;
					ComplexityW += p;
#if ULC_USE_PSYCHOACOUSTICS
					p   = ceilf(p);
					pNp = ceilf(pNp*LogScale);
					uint32_t ip   = (p   >= 0x1.0p32f) ? 0xFFFFFFFFu : (uint32_t)p;
					uint32_t ipNp = (pNp >= 0x1.0p32f) ? 0xFFFFFFFFu : (uint32_t)pNp;
					BufferEnergy  [n] = ip;
					BufferEnergyNp[n] = ipNp;
#endif
				}

				//! Use the log and linear sums to get the normalized entropy
				//! which will be used as a measure of the block's complexity
				if(ComplexityW) {
					//! Entropy doesn't directly translate to quality
					//! very well, so map it through a square-root curve
					float LogComplexityW = logf(ComplexityW * 0x1.6A09E6p0f); //! * Sqrt[2] to account for scaling in Log[x]
					float LogComplexity = Complexity / ComplexityW;
					Complexity = (LogComplexityW - LogComplexity) / (0x1.62E430p-1f*ComplexityScale); //! 0x1.62E430p-1 = 1/Log2[E] for change-of-base
					State->BlockComplexity = sqrtf(Complexity);
				} else State->BlockComplexity = 0.0f;
#if ULC_USE_PSYCHOACOUSTICS
				//! Buffer the psychoacoustics analysis into the last channel's buffer
				//! so that we can save a bit of memory to store these values into and
				//! avoid re-calculating for each channel, as the analysis isn't cheap
				Block_Transform_CalculatePsychoacoustics(MaskingNp, BufferEnergy, BufferEnergyNp, BlockSize, DecimationPattern);
#endif
			}

			//! Analyze each channel
			for(Chan=0;Chan<nChan;Chan++) {
				//! Analyze each subblock separately
				const float *BufferMDCTEnd = BufferMDCT + BlockSize;
				int SubBlockIdx;
				for(SubBlockIdx=0;(BufferMDCT < BufferMDCTEnd);SubBlockIdx++) {
					int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];

					//! Store the Neper-domain coefficients and sorting (importance) indices
					Block_Transform_WriteNepersAndIndices(
						BufferIndex,
#if ULC_USE_PSYCHOACOUSTICS
						MaskingNp,
#endif
						BufferMDCT,
						BufferNepers,
						&nNzCoef,
						SubBlockSize
					);

					//! Move to the next subblock
					BufferMDCT   += SubBlockSize;
					BufferNepers += SubBlockSize;
					BufferIndex  += SubBlockSize;
#if ULC_USE_PSYCHOACOUSTICS
					MaskingNp    += SubBlockSize;
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
			BufferMDCT   = State->TransformBuffer;
			BufferNepers = State->TransformNepers;
			BufferIndex  = (float*)State->TransformIndex;
			for(Chan=0;Chan<nChan;Chan++) {
				Block_Transform_BufferInterleave(BufferMDCT,   BufferTemp, BlockSize, Decimation);
				Block_Transform_BufferInterleave(BufferNepers, BufferTemp, BlockSize, Decimation);
				Block_Transform_BufferInterleave(BufferIndex,  BufferTemp, BlockSize, Decimation);
				BufferMDCT   += BlockSize;
				BufferNepers += BlockSize;
				BufferIndex  += BlockSize;
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
