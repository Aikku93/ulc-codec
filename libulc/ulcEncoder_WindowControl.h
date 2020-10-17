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
static inline float Block_Transform_GetWindowCtrl_TransientRatio(float R, float L) {
	//! NOTE: Decays are scaled by half the amplitude.
	if(R < 0.5f*L) { float t = 0.5f*L; L = R, R = t; }
	return R / L;
}
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
	for(n=0;n<2*BlockSize;n++) StepBuffer[n] = 0.0f;
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

	//! Spread the energy backwards a bit to account for pre-echo.
	//! The idea is to have the transient detector begin decimating
	//! at the onset of the transient, rather than immediately when
	//! the transient strikes, giving a less jarring transition.
	{
		//! Pre-echo decay curve: 6dB/ms (in X^2 domain):
		//!  E^(-Log[10]*1000*(dBDecayPerMs/10) / RateHz)
		float *Buf = StepBuffer + BlockSize*2;
		float a = expf(-0x1.596344p10f / RateHz);
		float b = 1.0f - a;
		float SampleLast = 0.0f;
		for(n=0;n<BlockSize*2;n++) {
			float v = *--Buf;
			*Buf = SampleLast = a*SampleLast + b*v;
		}
	}

	//! Offset by -BlockSize/2 relative to the data pointed to by
	//! `StepBuffer`. This is due to MDCT overlap placing the start
	//! of the lapping region at this point in time.
	StepBuffer += /*BlockSize - */BlockSize/2; //! x - x/2 == x/2

	//! Begin binary search for transient region until it centers between subblocks
	//! NOTE: Operating with SubBlockSize/2 should result in more optimized code.
	int Decimation = 1, SubBlockSize_2 = BlockSize / 2;
	for(;;) {
		//! Calculate the smooth-max energy of the left/middle/right sections
		//! (ie. compute a sparseness-compensating average energy)
		{
			//! NOTE: Adding a small bias avoids issues with division-by-zero
			const float Bias = 0x1.0p-31f, Bias2 = SQR(Bias);
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
				const __m256 vBias = _mm256_setr_ps(Bias2, Bias, Bias2, Bias, Bias2, Bias, Bias2, Bias);

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
				c = _mm256_add_ps(c, vBias);
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
				const __m128 vBias = _mm_setr_ps(Bias2, Bias, Bias2, Bias);
				a = _mm_shuffle_ps(vStepLL, vStepLLW, 0x88); //! {StepX[0,2],StepXW[0,2]}
				b = _mm_shuffle_ps(vStepLL, vStepLLW, 0xDD); //! {StepX[1,3],StepXW[1,3]}
				c = _mm_add_ps(a, b);                        //! {StepX[0+1],StepX[2+3],StepXW[0+1],StepXW[2+3]}
				a = _mm_shuffle_ps(vStepL, vStepLW, 0x88);
				b = _mm_shuffle_ps(vStepL, vStepLW, 0xDD);
				d = _mm_add_ps(a, b);
				a = _mm_shuffle_ps(c, d, 0x88);              //! {StepX[0+1],StepXW[0+1],StepY[0+1],StepYW[0+1]}
				b = _mm_shuffle_ps(c, d, 0xDD);              //! {StepX[2+3],StepXW[2+3],StepY[2+3],StepYW[2+3]}
				c = _mm_add_ps(a, b);                        //! {StepX,StepXW,StepY,StepYW}
				c = _mm_add_ps(c, vBias);
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
				c = _mm_add_ps(c, vBias);
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
			STEP_LL = StepLL + Bias2, STEP_LLW = StepLLW + Bias;
			STEP_L  = StepL  + Bias2, STEP_LW  = StepLW  + Bias;
			STEP_M  = StepM  + Bias2, STEP_MW  = StepMW  + Bias;
			STEP_R  = StepR  + Bias2, STEP_RW  = StepRW  + Bias;
#endif
		}

		//! Can we use window switching at all?
		//! NOTE: Limited to a minimum subblock size of 128 samples, and
		//! maximum decimation of 1/8 (Decimation > 8h: 1yyy = Decimation by 1/8)
		if(ULC_USE_WINDOW_SWITCHING && SubBlockSize_2 > 128/2 && Decimation < 0x8) {
			enum { POS_L, POS_M, POS_R};

			//! Determine the transient ratios (ie. how much energy jumps between subblocks)
			//! and then scale by the (squared) L2:L1 ratio to improve sensitivity
			float MaxLL  = STEP_LL / STEP_LLW;
			float MaxL   = STEP_L  / STEP_LW;
			float MaxM   = STEP_M  / STEP_MW;
			float MaxR   = STEP_R  / STEP_RW;
			float RatioL = Block_Transform_GetWindowCtrl_TransientRatio(MaxL, MaxLL);
			float RatioM = Block_Transform_GetWindowCtrl_TransientRatio(MaxM, MaxL);
			float RatioR = Block_Transform_GetWindowCtrl_TransientRatio(MaxR, MaxM);
			RatioL *= (STEP_L*SubBlockSize_2) / SQR(STEP_LW);
			RatioM *= (STEP_M*SubBlockSize_2) / SQR(STEP_MW);
			RatioR *= (STEP_R*SubBlockSize_2) / SQR(STEP_RW);

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
			if(TransientPos != POS_R) {
				if(TransientRatio > 2.0f) {
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
		}

		//! If we reach here, then the transient either lies in between the left/right
		//! boundaries, or exists on both sides equally. Either way, we can't improve
		//! it any further, so we stop here and allow overlap switching to take over.
		break;
	}

	//! Determine the overlap amount for the transient subblock
	int OverlapScale, SubBlockSize = SubBlockSize_2*2; {
		//! Find the maximum in the M+R region (ie. the transition region)
		//! and then divide this by the energy in the M+L region.
		float Ra = SQR(STEP_M  + STEP_R ) * SubBlockSize;
		float Rb = SQR(STEP_MW + STEP_RW) * (STEP_M + STEP_L);
		if(Ra < 0.5f*Rb) { float t = 0.5f*Rb; Rb = Ra, Ra = t; }
		if(Ra*0x1.0p-1f >= Rb) { //! a/b==2^1 is the first ratio to result in a ratio > 0
			//! a/b >= 2^13 gives the upper limit ratio of 7.0
			if(Ra*0x1.0p-13f < Rb) {
				//! 0x1.715476p0 = 1/Log[2], to get the log base-2
				//! Scale by 0.5 to apply a square root to the result
				float r = Ra / Rb;
				OverlapScale = (int)(0x1.715476p-1f*logf(r) + 0.5f);
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
//! EOF
/**************************************/
