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

//! Transform a block and prepare its coefficients
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
		State->TransientWindow,
		State->TransformTemp,
		State->WindowCtrlTaps,
		BlockSize,
		nChan,
		State->RateHz
	);
	int NextBlockOverlap; {
		int NextDecimation   = NextWindowCtrl >> 4;
		    NextBlockOverlap = BlockSize >> ULC_Helper_SubBlockDecimationPattern(NextDecimation)[0];
		if(ULC_Helper_TransientSubBlockIndex(NextDecimation) == 0) NextBlockOverlap >>= (NextWindowCtrl & 0x7);
	}

	//! Transform channels and insert keys for each codeable coefficient
	//! It's not /strictly/ required to calculate nNzCoef, but it can
	//! speed things up in the rate-control step
	int nNzCoef = 0; {
		int n, Chan;
		const float *ModulationWindow = State->ModulationWindow;
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
				//! Get the size of this subblock and the overlap at the next
				int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];
				int OverlapSize;
				if(BufferMDCT+SubBlockSize < BufferMDCTEnd) {
					OverlapSize = BlockSize >> DecimationPattern[SubBlockIdx+1];
					if(SubBlockIdx+1 == ThisTransientSubBlockIdx) OverlapSize >>= (WindowCtrl & 0x7);
				} else OverlapSize = NextBlockOverlap;

				//! Limit overlap to the maximum allowed by this subblock
				if(OverlapSize > SubBlockSize) OverlapSize = SubBlockSize;

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
					OverlapSize,
					ModulationWindow
				);

				//! Normalize spectrum, and accumulate amplitude by
				//! treating MDCT as Re and MDST as Im (akin to DFT)
				float Norm = 2.0f / SubBlockSize;
				for(n=0;n<SubBlockSize;n++) {
					float Re = (BufferMDCT[n] *= Norm);
					float Im = (BufferMDST[n] *= Norm);
					float Abs2 = SQR(Re) + SQR(Im);
#if ULC_USE_NOISE_CODING
					BufferTemp[n] = Abs2;
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
				//! This isn't perfect (adding even one more transient event results
				//! in the original problem once again), but is much better than
				//! doing nothing at all.
				//! NOTE: For noise analysis, we use a geometric mean to determine
				//! the correct noise amplitude. However, the geometric mean tends
				//! towards 1/E rather than 0.5 for white noise. So we scale by E
				//! to compensate and get unity amplitude, but then scale by 0.5
				//! to account for this filter's gain. As the computations work
				//! in the log domain, we take the logarithm here and then add
				//! Log[E/2] (0x1.3A37A0p-2) to change a multiply into an add.
				//! NOTE: We apply the filter over the squared samples, because
				//! this appears to give consistent results for some reason. The
				//! scaling then works out to have a gain of Sqrt[2], but this
				//! appears to sound better than unity gain for some reason.
				{
					float v, v2;
#if ULC_USE_NOISE_CODING
# define STORE_VALUE(n, Expr) \
	v2 = (Expr), v = sqrtf(v2), \
	Complexity += v2, ComplexityW += v, \
	BufferNoise[n] = ULC_FastLnApprox(v) + 0x1.3A37A0p-2f
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

#if ULC_USE_PSYCHOACOUSTICS
		//! Perform psychoacoustics analysis
		Block_Transform_CalculatePsychoacoustics(MaskingNp, BufferAmp2, (uint32_t*)BufferTemp, BlockSize, DecimationPattern);
#endif

		//! Perform importance analysis for all coefficients
		for(Chan=0;Chan<nChan;Chan++) {
			//! Analyze each subblock separately
			const float *BufferMDCTEnd = BufferMDCT + BlockSize;
			int SubBlockIdx;
			for(SubBlockIdx=0;(BufferMDCT < BufferMDCTEnd);SubBlockIdx++) {
				//! Store the sorting (importance) indices of a block's coefficients,
				//! and update the number of codeable non-zero coefficients
				int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];
				for(n=0;n<SubBlockSize;n++) {
					//! Coefficient inside codeable range?
					float Val = ABS(BufferMDCT[n]);
					if(Val < 0.5f*ULC_COEF_EPS) {
						BufferIndex[n] = -0x1.0p126f; //! Unusable coefficient; map to the end of the list
					} else {
						float ValNp = logf(Val);
						float MaskedValNp = ValNp;
#if ULC_USE_PSYCHOACOUSTICS
						//! Apply psychoacoustic corrections to this band energy
						MaskedValNp -= MaskingNp[n];
#endif
						//! Store the sort value for this coefficient
						BufferIndex[n] = MaskedValNp;
						nNzCoef++;
					}
				}

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
