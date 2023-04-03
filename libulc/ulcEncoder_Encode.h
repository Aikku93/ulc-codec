/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2023, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <math.h>
#include <stdint.h>
/**************************************/
#include "Fourier.h"
#include "ulcEncoder.h"
/**************************************/
#if ULC_USE_NOISE_CODING
# include "ulcEncoder_NoiseFill.h"
#endif
/**************************************/

//! The target memory is always aligned to 64 bytes, so just
//! use whatever is most performant on the target architecture
typedef uint32_t BitStream_t; //! <- MUST BE UNSIGNED
#define BISTREAM_NBITS (8u*sizeof(BitStream_t))

/**************************************/

ULC_FORCED_INLINE void Block_Encode_WriteNybble(BitStream_t x, BitStream_t **Dst, int *Size) {
	//! Push nybble
	*(*Dst) >>= 4;
	*(*Dst)  |= x << (BISTREAM_NBITS - 4);
	*Size += 4;

	//! Next byte?
	if((*Size)%BISTREAM_NBITS == 0) (*Dst)++;
}
ULC_FORCED_INLINE void Block_Encode_WriteQuantizer(int qi, BitStream_t **DstBuffer, int *Size, int Lead) {
	int s = qi - 5;
	if(Lead) {
		Block_Encode_WriteNybble(0xF, DstBuffer, Size);
	}
	if(s < 0xE) {
		//! Fh,0h..Dh: Quantizer change
		Block_Encode_WriteNybble(s, DstBuffer, Size);
	} else {
		//! Fh,Eh,0h..Ch: Quantizer change (extended precision)
		Block_Encode_WriteNybble(  0xE, DstBuffer, Size);
		Block_Encode_WriteNybble(s-0xE, DstBuffer, Size);
	}
}

/**************************************/

//! Build quantizer from weighted average
ULC_FORCED_INLINE int Block_Encode_BuildQuantizer(float MaxVal) {
	//! NOTE: `MaxVal` approaches (4/Pi)^2 as BlockSize
	//! approaches infinity, due to the p=1 norm of
	//! the MDCT matrix. Therefore, we use a bias of 5
	//! in the syntax to allow for a full range (that is
	//! to say: 7^2 * 2^-5 = 1.53125, 1.53125 >= 4/Pi).
	//! We then round this up to the nearest integer.
	int q = (int)ceilf(5.0f - 0x1.715476p0f*logf(MaxVal)); //! 0x1.715476p0 == 1/Ln[2] for change of base
	if(q < 5) q = 5;
	if(q > 5 + 0xE + 0xC) q = 5 + 0xE + 0xC; //! 5+Eh+Ch = Maximum extended-precision quantizer value (including a bias of 5)
	return q;
}

/**************************************/

//! Encode a range of coefficients
static inline int Block_Encode_EncodePass_WriteQuantizerZone(
	int           CurIdx,
	int           EndIdx,
	float         Quant,
	const float  *Coef,
#if ULC_USE_NOISE_CODING
	const float  *CoefNoise,
#endif
	const int    *CoefIdx,
	int           NextCodedIdx,
	int           nOutCoef,
	BitStream_t **DstBuffer,
	int          *Size
) {
	//! Write the coefficients
	do {
		//! Seek the next viable coefficient
		int Qn;
#if 1
		do if(CoefIdx[CurIdx] < nOutCoef) {
			//! We can only code +/-2..+/-7, so we
			//! check that Qn is inside this range
			Qn = ULC_CompandedQuantizeCoefficient(Coef[CurIdx]*Quant, 0x7);
			if(ABS(Qn) > 1) break;
		} while(++CurIdx < EndIdx);
		if(CurIdx >= EndIdx) break;
#else //! This code path can result in really bad noise when quantization gets extreme
		//! If the psymodel says we /need/ this coefficient,
		//! make sure we encode it rather than skip it, even
		//! if we have more distortion than if we skipped it
		while(CurIdx < EndIdx && CoefIdx[CurIdx] >= nOutCoef) CurIdx++;
		if(CurIdx >= EndIdx) break;
		Qn = ULC_CompandedQuantizeCoefficient(Coef[CurIdx]*Quant, 0x7);
		if(ABS(Qn) < 2) Qn = (Coef[CurIdx] < 0) ? (-2) : (+2);
#endif

		//! Code the zero runs, and get coefficient energy loss.
		//! The re-shaping basically accounts for how much energy we
		//! lose when coding zeros and/or noise, and adjusts the
		//! coefficient we code accordingly.
		float LossSum  = 0.0f;
		float LossSumW = 0.0f;
		int n, v, zR = CurIdx - NextCodedIdx;
		while(zR) {
			//! If we plan to skip two or less coefficients, try to encode
			//! them instead, as this is cheaper+better than filling with zero
			if(zR == 1) {
				int Qn1 = ULC_CompandedQuantizeCoefficient(Coef[NextCodedIdx]*Quant, 0x7);
				if(ABS(Qn1) > 1) {
					//! -7h..-2h, +2h..+7h: Normal coefficient
					Block_Encode_WriteNybble(Qn1, DstBuffer, Size);
					NextCodedIdx += 1;
					break;
				}
			}
			if(zR == 2) {
				int Qn1 = ULC_CompandedQuantizeCoefficient(Coef[NextCodedIdx+0]*Quant, 0x7);
				int Qn2 = ULC_CompandedQuantizeCoefficient(Coef[NextCodedIdx+1]*Quant, 0x7);
				if(ABS(Qn1) > 1 && ABS(Qn2) > 1) {
					//! -7h..-2h, +2h..+7h: Normal coefficient
					Block_Encode_WriteNybble(Qn1, DstBuffer, Size);
					Block_Encode_WriteNybble(Qn2, DstBuffer, Size);
					NextCodedIdx += 2;
					break;
				}
			}
#if ULC_USE_NOISE_CODING
			//! Determine the quantized coefficient for noise-fill mode
			//! NOTE: The range of noise samples is coded as a trade-off
			//! between noise-fill commands and coded coefficients, as
			//! noise-fill competes with coefficients we /want/ to code
			//! for the bit bandwidth available. Setting the start range
			//! too low favours noise-fill over psychoacoustically
			//! significant coefficients, but setting it too high can
			//! cause degradation in small block sizes as well as when
			//! we get many short runs.
			//! TODO: Is there a better min-number-of-coefficients?
			//! Too low and this fights too much with the coefficients
			//! we are actually trying to code (reducing the coding
			//! rate of them, and increasing the coding rate of noise).
			//! However, with small BlockSize (or SubBlockSize when
			//! decimating), smaller runs are useful.
			int NoiseQ = 0;
			if(zR >= 16) {
				v = zR - 16; if(v > 0x01FF) v = 0x01FF;
				n = v  + 16;
				NoiseQ = Block_Encode_EncodePass_GetNoiseQ(CoefNoise, NextCodedIdx, n, Quant);
			}
			if(NoiseQ) {
				//! 8h,Zh,Yh,Xh: 16 .. 527 noise fill (Xh != 0)
				Block_Encode_WriteNybble(0x8,  DstBuffer, Size);
				Block_Encode_WriteNybble(v>>5, DstBuffer, Size);
				Block_Encode_WriteNybble(v>>1, DstBuffer, Size);
				Block_Encode_WriteNybble((v&1) | ((NoiseQ-1)<<1), DstBuffer, Size);
			} else {
#endif
				//! Determine which run type to use and get the number of zeros coded
				//! NOTE: A short run takes 2 nybbles, and a long run takes 4 nybbles.
				//! So two short runs of maximum length code up to 32 zeros with the
				//! same efficiency as a long run, meaning that long runs start with
				//! 33 zeros.
				if(zR < 33) {
					//! 0h,0h..Fh: Zeros fill (1 .. 16 coefficients)
					v = zR - 1; if(v > 0xF) v = 0xF;
					n = v  + 1;
					Block_Encode_WriteNybble(0x0, DstBuffer, Size);
					Block_Encode_WriteNybble(v,   DstBuffer, Size);
				} else {
					//! 1h,Yh,Xh: 33 .. 542 zeros fill (Xh == 0)
					v = zR - 33; if(v > 0xFF) v = 0xFF;
					n = v  + 33;
					Block_Encode_WriteNybble(0x1,  DstBuffer, Size);
					Block_Encode_WriteNybble(v>>4, DstBuffer, Size);
					Block_Encode_WriteNybble(v>>0, DstBuffer, Size);
				}

				//! Finally, accumulate this run's energy loss
				int k;
				for(k=0;k<n;k++) {
					float v = Coef[NextCodedIdx+k];
					LossSum  += SQR(v);
					LossSumW += ABS(v);
				}
#if ULC_USE_NOISE_CODING
			}
#endif
			//! Skip the zeros
			NextCodedIdx += n;
			zR           -= n;
		}

		//! -7h..-2h, +2h..+7h: Normal coefficient
		{
			//! Sometimes, the math fails and we get an invalid coefficient.
			//! In these cases, just go with the original value.
			float Value = Coef[CurIdx];
			if(LossSum) Value += (LossSum/LossSumW) * ((Value < 0.0f) ? (-1.0f) : (+1.0f));
			int TestQn = ULC_CompandedQuantizeCoefficient(Value*Quant, 0x7);
			if(ABS(TestQn) > 1) Qn = TestQn;
		}
		Block_Encode_WriteNybble(Qn, DstBuffer, Size);
		NextCodedIdx++;

		//! Move to the next coefficient
		do CurIdx++; while(CurIdx < EndIdx && CoefIdx[CurIdx] >= nOutCoef);
	} while(CurIdx < EndIdx);
	return NextCodedIdx;
}

//! Encode a [sub]block
static inline void Block_Encode_EncodePass_WriteSubBlock(
	int           Idx,
	int           SubBlockSize,
	const float  *Coef,
#if ULC_USE_NOISE_CODING
	const float  *CoefNoise,
#endif
	const int    *CoefIdx,
	int           nOutCoef,
	BitStream_t **DstBuffer,
	int          *Size
) {
	//! Encode direct coefficients
	int EndIdx        = Idx+SubBlockSize;
	int NextCodedIdx  = Idx;
	int PrevQuant     = -1;
	int QuantStartIdx = -1;
	float QuantMin = 1000.0f, QuantMax = -1000.0f;
	do {
		//! Seek the next coefficient
		while(Idx < EndIdx && CoefIdx[Idx] >= nOutCoef) Idx++;

		//! Read coefficient and set the first quantizer's first coefficient index
		//! NOTE: Set NewMin=0.0 upon reaching the end. This causes the range
		//! check to fail in the next step, causing the final quantizer zone to dump
		//! if we had any data (if we don't, the check "passes" because NewMax==0).
		float NewMin = 0.0f;
		float NewMax = QuantMax;
		float CurLevel = 0.0f;
		if(Idx < EndIdx) {
			CurLevel = ABS(Coef[Idx]);
			NewMin = (CurLevel < QuantMin) ? CurLevel : QuantMin;
			NewMax = (CurLevel > QuantMax) ? CurLevel : QuantMax;
			if(QuantStartIdx == -1) QuantStartIdx = Idx;
		}

		//! Level out of range in this quantizer zone?
		const float MaxRange = 7.0f;
		if(NewMax > NewMin*MaxRange) {
			//! Write/update the quantizer
			int qi = Block_Encode_BuildQuantizer(QuantMax);
			if(qi != PrevQuant) {
				Block_Encode_WriteQuantizer(qi, DstBuffer, Size, PrevQuant != -1);
				PrevQuant = qi;
			}

			//! Write the quantizer zone we just searched through
			//! and start a new one from this coefficient
			NextCodedIdx = Block_Encode_EncodePass_WriteQuantizerZone(
				QuantStartIdx,
				Idx,
				(float)(1u << qi),
				Coef,
#if ULC_USE_NOISE_CODING
				CoefNoise,
#endif
				CoefIdx,
				NextCodedIdx,
				nOutCoef,
				DstBuffer,
				Size
			);
			QuantStartIdx = Idx;
			QuantMin = QuantMax = CurLevel;
		} else {
			//! Update ranges
			QuantMin = NewMin;
			QuantMax = NewMax;
		}
	} while(++Idx <= EndIdx);

	//! Decide what to do about the tail coefficients
	//! If we're at the edge of the block, it might work better to just fill with 0h
	int n = EndIdx - NextCodedIdx;
	if(n > 4) {
		//! If we coded anything, then we must specify the lead sequence
		if(PrevQuant != -1) {
			Block_Encode_WriteNybble(0xF, DstBuffer, Size);
		}

		//! Analyze the remaining data for noise-fill mode
#if ULC_USE_NOISE_CODING
		int NoiseQ = 0, NoiseDecay = 0;
		if(PrevQuant != -1 && n >= 16) { //! Don't use noise-fill for ultra-short tails
			Block_Encode_EncodePass_GetHFExtParams(
				CoefNoise,
				NextCodedIdx,
				n,
				(float)(1u << PrevQuant),
				&NoiseQ,
				&NoiseDecay
			);
		}
		if(NoiseQ) {
			//! Fh,Fh,Zh,Yh,Xh: Noise fill (to end; exp-decay)
			Block_Encode_WriteNybble(0xF,           DstBuffer, Size);
			Block_Encode_WriteNybble(NoiseQ-1,      DstBuffer, Size);
			Block_Encode_WriteNybble(NoiseDecay>>4, DstBuffer, Size);
			Block_Encode_WriteNybble(NoiseDecay,    DstBuffer, Size);
		} else {
#endif
			//! Fh,Eh,Fh: Stop
			Block_Encode_WriteNybble(0xE, DstBuffer, Size);
			Block_Encode_WriteNybble(0xF, DstBuffer, Size);
#if ULC_USE_NOISE_CODING
		}
#endif
	} else if(n > 0) {
		//! If we have less than 4 coefficients, it's cheaper to
		//! store a zero run than to do anything else.
		Block_Encode_WriteNybble(0x0, DstBuffer, Size);
		Block_Encode_WriteNybble(n-1, DstBuffer, Size);
	}
}

/**************************************/

//! Returns the block size (in bits) and the number of coded (non-zero) coefficients
static inline int Block_Encode_EncodePass(const struct ULC_EncoderState_t *State, void *_DstBuffer, int nOutCoef) {
	int BlockSize   = State->BlockSize;
	int Chan, nChan = State->nChan;
	const float *Coef      = State->TransformBuffer;
#if ULC_USE_NOISE_CODING
	const float *CoefNoise = State->TransformNoise;
#endif
	const int   *CoefIdx   = State->TransformIndex;
	BitStream_t *DstBuffer = _DstBuffer;

	//! Begin coding
	int Idx  = 0;
	int Size = 0; //! Block size (in bits)
	int WindowCtrl = State->WindowCtrl; {
		Block_Encode_WriteNybble(WindowCtrl, &DstBuffer, &Size);
		if(WindowCtrl & 0x8) Block_Encode_WriteNybble(WindowCtrl >> 4, &DstBuffer, &Size);
	}
	for(Chan=0;Chan<nChan;Chan++) {
		ULC_SubBlockDecimationPattern_t DecimationPattern = ULC_SubBlockDecimationPattern(WindowCtrl);
		do {
			int SubBlockSize = BlockSize >> (DecimationPattern&0x7);
			Block_Encode_EncodePass_WriteSubBlock(
				Idx,
				SubBlockSize,
				Coef,
#if ULC_USE_NOISE_CODING
				CoefNoise,
#endif
				CoefIdx,
				nOutCoef,
				&DstBuffer,
				&Size
			);
			Idx += SubBlockSize;
		} while(DecimationPattern >>= 4);
	}

	//! Align the output stream and pad size to bytes
	*DstBuffer >>= (-Size) % BISTREAM_NBITS;
	Size = (Size+7) &~ 7;
	return Size;
}

/**************************************/

#undef BISTREAM_NBITS

/**************************************/
//! EOF
/**************************************/
