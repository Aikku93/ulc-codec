/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#if defined(__AVX__) || defined(__FMA__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include <stdint.h>
/**************************************/
#include "Fourier.h"
#include "ulcEncoder_Psycho.h"
#include "ulcHelper.h"
/**************************************/

//! Write the sort values for all coefficients in a block
//! and return the number of codeable non-zero coefficients
//! NOTE:
//!  AnalysisPower is used to alter the preference for the currently-being-analyzed channel
static inline void Block_Transform_WriteSortValues(
	float *CoefIdx,
	const float *Energy,
	const float *EnergyNp,
	const float *CoefNp,
	int   *nNzCoef,
	int   BlockSize,
	float AnalysisPowerNp,
	float NyquistHz
) {
	int Band;
#if ULC_USE_PSYCHOACOUSTICS
#if 0 //! Log2[2*BlockSize] can be done in integer math
	float Gamma = 0x1.715476p0f*logf(2*BlockSize); //! Log2[2*BlockSize] == Log[2*BlockSize]/Log[2]
#else
	float Gamma = (float)(32 - __builtin_clz(BlockSize)); //! Log2[2*x] == 1+Log2[x], Log2[x] == 31-Clz_32bit[x]
#endif
	struct Block_Transform_MaskingState_t MaskingState;
	Block_Transform_MaskingState_Init(&MaskingState, Energy, EnergyNp, BlockSize, NyquistHz);
#else
	(void)Energy;
	(void)NyquistHz;
#endif
	for(Band=0;Band<BlockSize;Band++) {
		//! Inside the codeable range?
		float ValNp = CoefNp[Band];
		if(ValNp != ULC_COEF_NEPER_OUT_OF_RANGE) {
#if ULC_USE_PSYCHOACOUSTICS
			//! Apply psychoacoustic corrections to this band energy
			//! NOTE: Operate in the energy (X^2) domain, as this
			//! seems to result in slightly improved quality.
			ValNp += ValNp + Block_Transform_GetMaskedLevel(&MaskingState, Energy, EnergyNp, Band, BlockSize, Gamma);
#else
			//! Not so much a psychoacoustic optimization as an MDCT
			//! energy correction factor
			ValNp += ValNp + EnergyNp[Band];
#endif
			//! Store the sort value for this coefficient
			CoefIdx[Band] = ValNp + AnalysisPowerNp;
			(*nNzCoef)++;
		} else CoefIdx[Band] = -0x1.0p126f; //! Unusable coefficient; map to the end of the list
	}
}

/**************************************/

//! Get optimal log base-2 overlap and window scalings for transients
//! The idea is that overlap occurs roughly at the center of a block,
//! and if we have a transient placed right smack in the middle, we
//! can simply reduce the overlap to control its pre-echo. This is
//! combined with window switching so that the transient lies mostly
//! in the center of a subblock, at which point overlap scaling does
//! its job to reduce the pre-echo.
//! NOTE: StepBuffer must be 2*BlockSize in size.
//! NOTE: Bit codes for transient region coding, and their window sizes:
//!  First nybble:
//!   0xxx: No decimation. xxx = Overlap scaling
//!   1xxx: Decimate. xxx = Overlap scaling for the transient subblock
//!  Second nybble (when first nybble is 1xxx; otherwise, this is implicitly 0001):
//!   1xxx: Decimation by 1/8: Position = 0~7
//!    1000: N/8*,N/8,N/4,N/2
//!    1001: N/8,N/8*,N/4,N/2
//!    1010: N/4,N/8*,N/8,N/2
//!    1011: N/4,N/8,N/8*,N/2
//!    1100: N/2,N/8*,N/8,N/4
//!    1101: N/2,N/8,N/8*,N/4
//!    1110: N/2,N/4,N/8*,N/8
//!    1111: N/2,N/4,N/8,N/8*
//!   01xx: Decimation by 1/4: Position = 0~3
//!    0100: N/4*,N/4,N/2
//!    0101: N/4,N/4*,N/2
//!    0110: N/2,N/4*,N/4
//!    0111: N/2,N/4,N/4*
//!   001x: Decimation by 1/2: Position = 0~1
//!    0010: N/2*,N/2
//!    0011: N/2,N/2*
//!   0001: No decimation (not coded in the bitstream)
//!    0001: N/1*
//!  Transient subblocks are thus conveniently indexed via
//!  POPCNT (minus 1 to remove the unary count 'stop' bit)
//! NOTE: There is a lot of cross-multiplication in this function, the purpose
//! of which is to avoid division and its associated issues (eg. divide-by-0).
static inline int Block_Transform_GetWindowCtrl(
	const float *Data,
	const float *LastBlockData,
	float *StepBuffer,
	int BlockSize,
	int MinOverlap,
	int MaxOverlap,
	int nChan,
	int RateHz
) {
	int n, Chan;
	float StepData[8]; //! StepLL[W], StepL[W], StepM[W], StepR[W]. Janky setup for vectorization
#define STEP_LL  StepData[0]
#define STEP_LLW StepData[1]
#define STEP_L   StepData[2]
#define STEP_LW  StepData[3]
#define STEP_M   StepData[4]
#define STEP_MW  StepData[5]
#define STEP_R   StepData[6]
#define STEP_RW  StepData[7]

	//! Perform a highpass filter to get the step/transient energies
	for(n=0;n<BlockSize;n++) StepBuffer[n] = 0.0f;
	for(Chan=0;Chan<nChan;Chan++) {
		const float *Src;
		float v, d, SampleLast, *Dst = StepBuffer;

		//! Get step energy for last block
		//! NOTE: The first step sample will always be zero, but this is ok
		//! as the formulation is finding a smooth-max rather than average.
		Src = LastBlockData + Chan*BlockSize;
		SampleLast = Src[0];
		for(n=0;n<BlockSize;n++) {
			v = *Src++, d = v - SampleLast, SampleLast = v, *Dst++ += SQR(d);
		}

		//! Get step energy for this block
		Src = Data + Chan*BlockSize;
		for(n=0;n<BlockSize;n++) {
			v = *Src++, d = v - SampleLast, SampleLast = v, *Dst++ += SQR(d);
		}
	}

	//! Spread the energy backwards a bit to account for pre-echo
	{
		//! Pre-echo decay curve: 24dB/ms (in X^2 domain):
		//!  E^(-Log[10]*1000*(dBDecayPerMs/10) / RateHz)
		float *Buf = StepBuffer + BlockSize*2;
		float a = expf(-0x1.596344p12f / RateHz);
		float b = 1.0f - a;
		float SampleLast = 0.0f;
		for(n=0;n<BlockSize*2;n++) {
			float v = *--Buf;
			*Buf = SampleLast = a*SampleLast + b*v;
		}
	}

	//! NOTE: Offset by -BlockSize/2 relative to the data pointed to by
	//! `StepBuffer`. This is due to MDCT overlap placing the start of
	//! the lapping region at this point in time.
	StepBuffer += /*BlockSize - */BlockSize/2; //! x - x/2 == x/2

	//! Begin binary search for transient region until it centers within a subblock
	//! NOTE: Operating with SubBlockSize/2 should result in more optimized code.
	int Decimation = 1, SubBlockSize_2 = BlockSize / 2;
	for(;;) {
		//! Calculate the smooth-max energy of the left/middle/right sections
		//! (ie. compute a sparseness-compensating average energy)
		{
#if defined(__AVX__)
			__m256 vStepLL = _mm256_setzero_ps(), vStepLLW = _mm256_setzero_ps();
			__m256 vStepL  = _mm256_setzero_ps(), vStepLW  = _mm256_setzero_ps();
			__m256 vStepM  = _mm256_setzero_ps(), vStepMW  = _mm256_setzero_ps();
			__m256 vStepR  = _mm256_setzero_ps(), vStepRW  = _mm256_setzero_ps();
			for(n=0;n<SubBlockSize_2;n+=8) {
				__m256 dLL, dL, dM, dR;

				dLL = _mm256_load_ps(StepBuffer + n - SubBlockSize_2);
				dL  = _mm256_load_ps(StepBuffer + n);
				dM  = _mm256_load_ps(StepBuffer + n + SubBlockSize_2);
				dR  = _mm256_load_ps(StepBuffer + n + SubBlockSize_2*2);
				vStepLLW = _mm256_add_ps(vStepLLW, dLL);
				vStepLW  = _mm256_add_ps(vStepLW,  dL);
				vStepMW  = _mm256_add_ps(vStepMW,  dM);
				vStepRW  = _mm256_add_ps(vStepRW,  dR);
#if defined(__FMA__)
				vStepLL = _mm256_fmadd_ps(dLL, dLL, vStepLL);
				vStepL  = _mm256_fmadd_ps(dL,  dL,  vStepL);
				vStepM  = _mm256_fmadd_ps(dM,  dM,  vStepM);
				vStepR  = _mm256_fmadd_ps(dR,  dR,  vStepR);
#else
				dLL = _mm256_mul_ps(dLL, dLL);
				dL  = _mm256_mul_ps(dL,  dL);
				dM  = _mm256_mul_ps(dM,  dM);
				dR  = _mm256_mul_ps(dR,  dR);
				vStepLL = _mm256_add_ps(vStepLL, dLL);
				vStepL  = _mm256_add_ps(vStepL,  dL);
				vStepM  = _mm256_add_ps(vStepM,  dM);
				vStepR  = _mm256_add_ps(vStepR,  dR);
#endif
			}
			{
				__m256 a, b, c, d;

				//! a = L: {StepA[0+1],StepA[2+3],StepAW[0+1],StepAW[2+3]}
				//!     H: {StepA[4+5],StepA[6+7],StepAW[4+5],StepAW[6+7]}
				//! b = L: {StepB[0+1],StepB[2+3],StepBW[0+1],StepBW[2+3]}
				//!     H: {StepB[4+5],StepB[6+7],StepBW[4+5],StepBW[6+7]}
				a = _mm256_hadd_ps(vStepLL, vStepLLW);
				b = _mm256_hadd_ps(vStepL,  vStepLW);
				c = _mm256_hadd_ps(vStepM,  vStepMW);
				d = _mm256_hadd_ps(vStepR,  vStepRW);
				//! a = L: {StepA[0+1+2+3],StepAW[0+1+2+3],StepB[0+1+2+3],StepBW[0+1+2+3]}
				//!     H: {StepA[4+5+6+7],StepAW[4+5+6+7],StepB[4+5+6+7],StepBW[4+5+6+7]}
				a = _mm256_hadd_ps(a, b);
				b = _mm256_hadd_ps(c, d);
				//! c = L: {StepA[0+1+2+3],StepAW[0+1+2+3],StepB[0+1+2+3],StepBW[0+1+2+3]}
				//!     H: {StepC[0+1+2+3],StepCW[0+1+2+3],StepD[0+1+2+3],StepDW[0+1+2+3]}
				//! d = L: {StepA[4+5+6+7],StepAW[4+5+6+7],StepB[4+5+6+7],StepBW[4+5+6+7]}
				//!     H: {StepC[4+5+6+7],StepCW[4+5+6+7],StepD[4+5+6+7],StepDW[4+5+6+7]}
				c = _mm256_permute2f128_ps(a, b, 0x20);
				d = _mm256_permute2f128_ps(a, b, 0x31);
				c = _mm256_add_ps(c, d);
				_mm256_store_ps(StepData, c);
			}
#elif defined(__SSE__)
			__m128 vStepLL = _mm_setzero_ps(), vStepLLW = _mm_setzero_ps();
			__m128 vStepL  = _mm_setzero_ps(), vStepLW  = _mm_setzero_ps();
			__m128 vStepM  = _mm_setzero_ps(), vStepMW  = _mm_setzero_ps();
			__m128 vStepR  = _mm_setzero_ps(), vStepRW  = _mm_setzero_ps();
			for(n=0;n<SubBlockSize_2;n+=4) {
				__m128 dLL, dL, dM, dR;

				dLL = _mm_load_ps(StepBuffer + n - SubBlockSize_2);
				dL  = _mm_load_ps(StepBuffer + n);
				dM  = _mm_load_ps(StepBuffer + n + SubBlockSize_2);
				dR  = _mm_load_ps(StepBuffer + n + SubBlockSize_2*2);
				vStepLLW = _mm_add_ps(vStepLLW, dLL);
				vStepLW  = _mm_add_ps(vStepLW,  dL);
				vStepMW  = _mm_add_ps(vStepMW,  dM);
				vStepRW  = _mm_add_ps(vStepRW,  dR);
#if defined(__FMA__)
				vStepLL = _mm_fmadd_ps(dLL, dLL, vStepLL);
				vStepL  = _mm_fmadd_ps(dL,  dL,  vStepL);
				vStepM  = _mm_fmadd_ps(dM,  dM,  vStepM);
				vStepR  = _mm_fmadd_ps(dR,  dR,  vStepR);
#else
				dLL = _mm_mul_ps(dLL, dLL);
				dL  = _mm_mul_ps(dL,  dL);
				dM  = _mm_mul_ps(dM,  dM);
				dR  = _mm_mul_ps(dR,  dR);
				vStepLL = _mm_add_ps(vStepLL, dLL);
				vStepL  = _mm_add_ps(vStepL,  dL);
				vStepM  = _mm_add_ps(vStepM,  dM);
				vStepR  = _mm_add_ps(vStepR,  dR);
#endif
			}
			{
				__m128 a, b, c, d;
				a = _mm_shuffle_ps(vStepLL, vStepLLW, 0x88); //! {StepX[0,2],StepXW[0,2]}
				b = _mm_shuffle_ps(vStepLL, vStepLLW, 0xDD); //! {StepX[1,3],StepXW[1,3]}
				c = _mm_add_ps(a, b);                        //! {StepX[0+1],StepX[2+3],StepXW[0+1],StepXW[2+3]}
				a = _mm_shuffle_ps(vStepL, vStepLW, 0x88);
				b = _mm_shuffle_ps(vStepL, vStepLW, 0xDD);
				d = _mm_add_ps(a, b);
				a = _mm_shuffle_ps(c, d, 0x88);              //! {StepX[0+1],StepXW[0+1],StepY[0+1],StepYW[0+1]}
				b = _mm_shuffle_ps(c, d, 0xDD);              //! {StepX[2+3],StepXW[2+3],StepY[2+3],StepYW[2+3]}
				c = _mm_add_ps(a, b);                        //! {StepX,StepXW,StepY,StepYW}
				_mm_store_ps(StepData, c);
				a = _mm_shuffle_ps(vStepM, vStepMW, 0x88);
				b = _mm_shuffle_ps(vStepM, vStepMW, 0xDD);
				c = _mm_add_ps(a, b);
				a = _mm_shuffle_ps(vStepR, vStepRW, 0x88);
				b = _mm_shuffle_ps(vStepR, vStepRW, 0xDD);
				d = _mm_add_ps(a, b);
				a = _mm_shuffle_ps(c, d, 0x88);
				b = _mm_shuffle_ps(c, d, 0xDD);
				c = _mm_add_ps(a, b);
				_mm_store_ps(StepData+4, c);
			}
#else
			float StepLL = 0.0f, StepLLW = 0.0f;
			float StepL  = 0.0f, StepLW  = 0.0f;
			float StepM  = 0.0f, StepMW  = 0.0f;
			float StepR  = 0.0f, StepRW  = 0.0f;
			for(n=0;n<SubBlockSize_2;n++) {
				float dLL = StepBuffer[n-SubBlockSize_2];
				float dL  = StepBuffer[n];
				float dM  = StepBuffer[n+SubBlockSize_2];
				float dR  = StepBuffer[n+SubBlockSize_2*2];
				StepLL += SQR(dLL), StepLLW += dLL;
				StepL  += SQR(dL),  StepLW  += dL;
				StepM  += SQR(dM),  StepMW  += dM;
				StepR  += SQR(dR),  StepRW  += dR;
			}
			STEP_LL = StepLL, STEP_LLW = StepLLW;
			STEP_L  = StepL,  STEP_LW  = StepLW;
			STEP_M  = StepM,  STEP_MW  = StepMW;
			STEP_R  = StepR,  STEP_RW  = StepRW;
#endif
		}

		//! Normalize the step values
		{
			//! NOTE: Avoid division-by-zero issues by adding a tiny bias.
			//! NOTE: Keep the bias tiny by scaling up the values prior to biasing.
			//! The maximum value that STEP_X can take on is (2^2)^2 * (MaxBlockSize/2),
			//! giving us a large margin to scale into.
			const float NormScale = 0x1.0p+64f;
			const float StepBias  = 0x1.0p-63f, StepBias2 = StepBias*StepBias;
			STEP_LL = (NormScale*STEP_LL + StepBias2) / (NormScale*STEP_LLW + StepBias);
			STEP_L  = (NormScale*STEP_L  + StepBias2) / (NormScale*STEP_LW  + StepBias);
			STEP_M  = (NormScale*STEP_M  + StepBias2) / (NormScale*STEP_MW  + StepBias);
			STEP_R  = (NormScale*STEP_R  + StepBias2) / (NormScale*STEP_RW  + StepBias);
		}

		//! Can we use window switching at all?
		//! NOTE: Limited to a minimum subblock size of 128 samples, and
		//! maximum decimation of 1/8 (Decimation > 8h: 1yyy = Decimation by 1/8)
		if(ULC_USE_WINDOW_SWITCHING && SubBlockSize_2 > 128/2 && Decimation < 0x8) {
			enum { POS_L, POS_M, POS_R};

			//! Determine the transient ratios (ie. how much energy jumps between subblocks)
			//! NOTE: Decays are scaled by half the amplitude.
			float RatioL = (2.0f*STEP_L > STEP_LL) ? (STEP_L / STEP_LL) : (STEP_LL / (2.0f*STEP_L));
			float RatioM = (2.0f*STEP_M > STEP_L)  ? (STEP_M / STEP_L)  : (STEP_L  / (2.0f*STEP_M));
			float RatioR = (2.0f*STEP_R > STEP_M)  ? (STEP_R / STEP_M)  : (STEP_M  / (2.0f*STEP_R));

			//! Determine which segment (L/M/R) is most transient
			//! NOTE: R corresponds to the end of the block pointed to by
			//! `Data`, which can use overlap switching to avoid decimation.
			int   TransientPos   = POS_L;
			float TransientRatio = RatioL;
			if(RatioM > TransientRatio) {
				TransientPos   = POS_M;
				TransientRatio = RatioM;
			}
			if(RatioR > TransientRatio) {
				TransientPos   = POS_R;
				TransientRatio = RatioR;
			}

			//! If the transient is not on the overlap boundary and
			//! is still significant, decimate the subblock further
			if(TransientPos != POS_R && TransientRatio > 2.0f) {
				//! Update the decimation pattern and continue
				if(TransientPos == POS_L)
					Decimation  = (Decimation<<1) | 0;
				else
					Decimation  = (Decimation<<1) | 1,
					StepBuffer += SubBlockSize_2;
				SubBlockSize_2 /= 2;
				continue;
			}
		}

		//! If we reach here, then the transient either lies in between the left/right
		//! boundaries, or exists on both sides equally. Either way, we can't improve
		//! it any further, so we stop here and allow overlap switching to take over.
		break;
	}

	//! Determine the overlap amount for the transient subblock
	int OverlapScale, SubBlockSize = SubBlockSize_2*2; {
		//! R/M or M/L; when R/M isn't a transient, then odds are that M/L is
		float Ra = STEP_R;
		float Rb = STEP_M;
		if(Ra < Rb) Ra = STEP_M, Rb = STEP_L;
		if(Ra*0x1.6A09E6p-1f >= Rb) { //! a/b==2^0.5 is the first ratio to result in a ratio > 0
			//! a/b >= 2^6.5 gives the upper limit ratio of 7.0
			if(Ra*0x1.6A09E6p-7f < Rb) {
				//! 0x1.715476p0 = 1/Log[2], to get the log base-2
				float r = Ra / Rb;
				OverlapScale = (int)(0x1.715476p0f*logf(r) + 0.5f);
			} else OverlapScale = 7;
			while(OverlapScale > 0 && (SubBlockSize >> OverlapScale) < MinOverlap) OverlapScale--;
			while(OverlapScale < 7 && (SubBlockSize >> OverlapScale) > MaxOverlap) OverlapScale++;
		} else OverlapScale = 0;
	}

	//! Return the combined overlap+window switching parameters
	return OverlapScale + 0x8*(SubBlockSize != BlockSize) + 0x10*Decimation;
#undef STEP_LL
#undef STEP_LLW
#undef STEP_L
#undef STEP_LW
#undef STEP_M
#undef STEP_MW
#undef STEP_R
#undef STEP_RW
}

/**************************************/

//! Transform a block and prepare its coefficients
static float *Block_Transform_SortComparator_KeysArray;
static int Block_Transform_SortComparator(const void *a, const void *b) {
	const float *BufferIdx = Block_Transform_SortComparator_KeysArray;
	return (BufferIdx[*(int*)a] < BufferIdx[*(int*)b]) ? (+1) : (-1);
}
static inline void Block_Transform_ScaleAndToNepers(float *Coef, float *CoefNp, int N) {
	int n;
	for(n=0;n<N;n++) {
		float v = Coef[n] * 2.0f/N;
		Coef  [n] = v;
		CoefNp[n] = (ABS(v) < 0.5f*ULC_COEF_EPS) ? ULC_COEF_NEPER_OUT_OF_RANGE : logf(ABS(v));
	}
}
static inline void Block_Transform_ComputePowerSpectrum(float *Power, float *PowerNp, const float *Re, const float *Im, int N) {
	int i;
	for(i=0;i<N;i++) {
		//! `v` should technically be scaled by (2/N)^2, but this power
		//! spectrum is only used for psychoacoustic analysis on this
		//! buffer /only/, and so scaling really doesn't matter here
		float v = SQR(Re[i]) + SQR(Im[i]);
		Power  [i] = v;
		PowerNp[i] = (ABS(v) < 0x1.0p-126f) ? ULC_COEF_NEPER_OUT_OF_RANGE : logf(ABS(v)); //! Do not allow subnormals
	}
}
static void Block_Transform_BufferInterleave(float *Buf, float *Tmp, int BlockSize, int Decimation) {
	//! The interleaving patterns here were chosen to try and
	//! optimize coefficient clustering across block sizes
	int n; for(n=0;n<BlockSize;n++) Tmp[n] = Buf[n];
	int n2 = BlockSize/2;
	int n4 = BlockSize/4;
	int n8 = BlockSize/8;
	switch(Decimation >> 1) { //! Lowermost bit only controls which subblock gets overlap scaling, so ignore it
		//! 001x: a=N/2, b=N/2
		case 0b001: {
			for(n=0;n<n2;n++) {
				*Buf++ = Tmp[  +n]; //! a
				*Buf++ = Tmp[n2+n]; //! b
			}
		} break;

		//! 010x: a=N/4, b=N/4, c=N/2
		case 0b010: {
			for(n=0;n<n4;n++) {
				*Buf++ = Tmp[  +n*1+0]; //! a
				*Buf++ = Tmp[n4+n*1+0]; //! b
				*Buf++ = Tmp[n2+n*2+0]; //! c0
				*Buf++ = Tmp[n2+n*2+1]; //! c1
			}
		} break;

		//! 011x: a=N/2, b=N/4, c=N/4
		case 0b011: {
			for(n=0;n<n4;n++) {
				*Buf++ = Tmp[     +n*2+0]; //! a0
				*Buf++ = Tmp[     +n*2+1]; //! a1
				*Buf++ = Tmp[n2   +n*1+0]; //! b
				*Buf++ = Tmp[n2+n4+n*1+0]; //! c
			}
		} break;

		//! 100x: a=N/8, b=N/8, c=N/4, d=N/2
		case 0b100: {
			for(n=0;n<n8;n++) {
				*Buf++ = Tmp[  +n*1+0]; //! a
				*Buf++ = Tmp[n8+n*1+0]; //! b
				*Buf++ = Tmp[n4+n*2+0]; //! c0
				*Buf++ = Tmp[n4+n*2+1]; //! c1
				*Buf++ = Tmp[n2+n*4+0]; //! d0
				*Buf++ = Tmp[n2+n*4+1]; //! d1
				*Buf++ = Tmp[n2+n*4+2]; //! d2
				*Buf++ = Tmp[n2+n*4+3]; //! d3
			}
		} break;

		//! 101x: a=N/4, b=N/8, c=N/8, d=N/2
		case 0b101: {
			for(n=0;n<n8;n++) {
				*Buf++ = Tmp[     +n*2+0]; //! a0
				*Buf++ = Tmp[     +n*2+1]; //! a1
				*Buf++ = Tmp[n4   +n*1+0]; //! b
				*Buf++ = Tmp[n4+n8+n*1+0]; //! c
				*Buf++ = Tmp[n2   +n*4+0]; //! d0
				*Buf++ = Tmp[n2   +n*4+1]; //! d1
				*Buf++ = Tmp[n2   +n*4+2]; //! d2
				*Buf++ = Tmp[n2   +n*4+3]; //! d3
			}
		} break;

		//! 110x: a=N/2, b=N/8, c=N/8, d=N/4
		case 0b110: {
			for(n=0;n<n8;n++) {
				*Buf++ = Tmp[     +n*4+0]; //! a0
				*Buf++ = Tmp[     +n*4+1]; //! a1
				*Buf++ = Tmp[     +n*4+2]; //! a2
				*Buf++ = Tmp[     +n*4+3]; //! a3
				*Buf++ = Tmp[n2   +n*1+0]; //! b
				*Buf++ = Tmp[n2+n8+n*1+0]; //! c
				*Buf++ = Tmp[n2+n4+n*2+0]; //! d0
				*Buf++ = Tmp[n2+n4+n*2+1]; //! d1
			}
		} break;

		//! 111x: a=N/2, b=N/4, c=N/8, d=N/8
		case 0b111: {
			for(n=0;n<n8;n++) {
				*Buf++ = Tmp[        +n*4+0]; //! a0
				*Buf++ = Tmp[        +n*4+1]; //! a1
				*Buf++ = Tmp[        +n*4+2]; //! a2
				*Buf++ = Tmp[        +n*4+3]; //! a3
				*Buf++ = Tmp[n2      +n*2+0]; //! b0
				*Buf++ = Tmp[n2      +n*2+1]; //! b1
				*Buf++ = Tmp[n2+n4   +n*1+0]; //! c
				*Buf++ = Tmp[n2+n4+n8+n*1+0]; //! d
			}
		} break;
	}
}
static int Block_Transform(struct ULC_EncoderState_t *State, const float *Data, float PowerDecay) {
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;

	//! Get the window control parameters for this block and the next
	int WindowCtrl     = State->WindowCtrl     = State->NextWindowCtrl;
	int NextWindowCtrl = State->NextWindowCtrl = Block_Transform_GetWindowCtrl(
		Data,
		State->SampleBuffer,
		State->TransformTemp,
		BlockSize,
		State->MinOverlap,
		State->MaxOverlap,
		nChan,
		State->RateHz
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
		float  AnalysisPowerNp = 0.0f; PowerDecay = logf(PowerDecay);
		float *BufferSamples   = State->SampleBuffer;
		float *BufferTransform = State->TransformBuffer;
		float *BufferNepers    = State->TransformNepers;
		float *BufferIndex     = (float*)State->TransformIndex;
		float *BufferFwdLap    = State->TransformFwdLap;
		float *BufferTemp      = State->TransformTemp;

		//! Transform the input data
		int ThisTransientSubBlockIdx = ULC_Helper_TransientSubBlockIndex(WindowCtrl >> 4);
		const int8_t *DecimationPattern = ULC_Helper_SubBlockDecimationPattern(WindowCtrl >> 4);
		for(Chan=0;Chan<nChan;Chan++) {
			//! Process each subblock sequentially
			float *BufferTransformEnd = BufferTransform + BlockSize;
			int SubBlockIdx;
			for(SubBlockIdx=0;(BufferTransform < BufferTransformEnd);SubBlockIdx++) {
				//! Get the size of this subblock and its overlap
				int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];
				int OverlapSize = BlockSize >> DecimationPattern[SubBlockIdx];
				if(SubBlockIdx == ThisTransientSubBlockIdx) OverlapSize >>= (WindowCtrl & 0x7);

				//! If the next block's first subblock is smaller than
				//! this one, set overlap to the smaller of the two.
				//! NOTE: This is literally THE ONLY reason that there is a
				//! coding delay in this codec (and most MDCT codecs).
				if(BufferTransform+SubBlockSize < BufferTransformEnd) {
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

				//! Perform the actual analyses
				float *BufferEnergy   = BufferTemp;
				float *BufferEnergyNp = BufferTemp + SubBlockSize;
				Fourier_MDCT(
					BufferTransform,
					SmpBuf,
					BufferFwdLap + (BlockSize-SubBlockSize)/2,
					BufferTemp,
					SubBlockSize,
					OverlapSize,
					BufferNepers //! MDST coefficients stored here temporarily
				);
				Block_Transform_ComputePowerSpectrum(BufferEnergy, BufferEnergyNp, BufferTransform, BufferNepers, SubBlockSize);
				Block_Transform_ScaleAndToNepers(BufferTransform, BufferNepers, SubBlockSize);
				Block_Transform_WriteSortValues(
					BufferIndex,
					BufferEnergy,
					BufferEnergyNp,
					BufferNepers,
					&nNzCoef,
					SubBlockSize,
					AnalysisPowerNp,
					State->RateHz * 0.5f
				);

				//! Move to the next subblock
				BufferTransform += SubBlockSize;
				BufferNepers    += SubBlockSize;
				BufferIndex     += SubBlockSize;
			}

			//! Cache the sample data for the next block
			for(n=0;n<BlockSize;n++) BufferSamples[n-BlockSize] = *Data++;

			//! Move to the next channel
			BufferFwdLap    += BlockSize/2u;
			AnalysisPowerNp += PowerDecay;
		}

		//! Interleave the transform data for coding
		if(WindowCtrl & 0x8) {
			int Decimation = WindowCtrl >> 4;
			BufferTransform = State->TransformBuffer;
			BufferNepers    = State->TransformNepers;
			BufferIndex     = (float*)State->TransformIndex;
			for(Chan=0;Chan<nChan;Chan++) {
				Block_Transform_BufferInterleave(BufferTransform, BufferTemp, BlockSize, Decimation);
				Block_Transform_BufferInterleave(BufferNepers,    BufferTemp, BlockSize, Decimation);
				Block_Transform_BufferInterleave(BufferIndex,     BufferTemp, BlockSize, Decimation);
				BufferTransform += BlockSize;
				BufferNepers    += BlockSize;
				BufferIndex     += BlockSize;
			}
		}
	}

	//! Create the coefficient sorting indices
	//! TODO: Be more efficient; this is suboptimal, but
	//! fixing would require implementing a customized
	//! sorting routine
	{
		int Idx, nIdx = nChan * BlockSize;
		int *BufferTmp = (int*)State->TransformTemp;
		int *BufferIdx = State->TransformIndex;
		for(Idx=0;Idx<nIdx;Idx++) BufferTmp[Idx] = Idx;
		Block_Transform_SortComparator_KeysArray = (float*)State->TransformIndex;
		qsort(BufferTmp, nIdx, sizeof(int), Block_Transform_SortComparator);
		for(Idx=0;Idx<nIdx;Idx++) BufferIdx[BufferTmp[Idx]] = Idx;
	}
	return nNzCoef;
}

/**************************************/
//! EOF
/**************************************/
