/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <math.h>
#include <stdint.h>
/**************************************/
#include "Fourier.h"
#include "ulcEncoder.h"
/**************************************/

static inline __attribute__((always_inline)) void Block_Encode_WriteNybble(uint8_t x, uint8_t **Dst, int *Size) {
	//! Push nybble
	*(*Dst) >>= 4;
	*(*Dst)  |= x << 4;
	*Size += 4;

	//! Next byte?
	if((*Size)%8u == 0) (*Dst)++;
}
static inline __attribute__((always_inline)) void Block_Encode_WriteQuantizer(int qi, uint8_t **DstBuffer, int *Size, int Lead) {
	//! 8h,0h,0h..Eh[,0h..Ch]: Quantizer change
	int s = qi - 5;
	if(Lead) {
		Block_Encode_WriteNybble(0x8, DstBuffer, Size);
		Block_Encode_WriteNybble(0x0, DstBuffer, Size);
	}
	if(s < 0xE) Block_Encode_WriteNybble(s, DstBuffer, Size);
	else {
		//! 8h,0h,Eh,0h..Ch: Extended-precision quantizer
		Block_Encode_WriteNybble(  0xE, DstBuffer, Size);
		Block_Encode_WriteNybble(s-0xE, DstBuffer, Size);
	}
}

/**************************************/

//! Build quantizer from weighted average
//! NOTE: The quantizer scale is 2^-x, so when this is set to the minimum scaling
//!       (ie. x==5, setting q=1/32), then the maximum codeable value becomes
//!       7^2/32 == 1.53125. This is higher than 1.0, so we use this extra space
//!       as a sort of headroom, and shift the rounding point accordingly.
//! NOTE: The average is performed over the Neper-domain coefficients, so there is
//!       no need to apply a further logarithm here to get the base-2 logarithm.
static inline int Block_Encode_BuildQuantizer(float Sum, float Weight) {
	//! NOTE: `q` will always be greater than 5 due to the bias
	int q = (int)(0x1.675768p2f - 0x1.715476p0f*Sum / Weight); //! 0x1.715476p0 == 1/Ln[2], as input is in natural log
	if(q < 5) q = 5; //! Sometimes happens because of overflow?
	if(q > 5 + 0xE + 0xC) q = 5 + 0xE + 0xC; //! 5+Eh+Ch = Maximum extended-precision quantizer value (including a bias of 5)
	return q;
}

//! Quantize coefficient
static inline int Block_Encode_Quantize(float v, float q) {
	float av = ABS(v);
	int vq = (int)(sqrtf(av*q) + 0.5f);
	if(vq == 0) return vq;
#if 1 //! Optimal rounding; the above only ever goes /over/ by +1, and never under
	float dl = ABS(av*q - SQR(vq-1));
	float d  = ABS(av*q - SQR(vq  ));
	vq -= (dl < d);
#endif
	return (v < 0.0f) ? (-vq) : (+vq);
}

/**************************************/

//! Returns the block size (in bits) and the number of coded (non-zero) coefficients
static inline __attribute__((always_inline)) int Block_Encode_EncodePass_WriteQuantizerZone(
	int       CurIdx,
	int       EndIdx,
	float     QuantSum,
	float     QuantWeight,
	const float *Coef,
	const int   *CoefIdx,
	int       NextCodedIdx,
	int      *PrevQuant,
	uint8_t **DstBuffer,
	int      *Size,
	int       nOutCoef
) {
	//! Write/update the quantizer
	float q; {
		int qi = Block_Encode_BuildQuantizer(QuantSum, QuantWeight);
		if(qi != *PrevQuant) {
			Block_Encode_WriteQuantizer(qi, DstBuffer, Size, *PrevQuant != -1);
			*PrevQuant = qi;
		}
		q = (float)(1u << qi);
	}

	//! Write the coefficients
	do {
		//! Target coefficient doesn't collapse to 0?
		if(ABS(Coef[CurIdx]*q) >= SQR(0.5f)) { //! Block_Encode_Quantize(Coef[CurIdx], q) != 0
			//! Code the zero runs
			int zR = CurIdx - NextCodedIdx;
			while(zR >= 3) {
				Block_Encode_WriteNybble(0x8, DstBuffer, Size);

				//! Determine which run type to use and  get the number of zeros coded
				int n;
				if(zR < 17) {
					//! 8h,1h..Eh: 3 .. 16 zeros
					n = zR;
					Block_Encode_WriteNybble(n-2, DstBuffer, Size);
				} else {
					//! 8h,Fh,Yh,Xh: 17 .. 272 zeros
					int v = zR-17; if(v > 0xFF) v = 0xFF;
					n = v + 17;
					Block_Encode_WriteNybble(0xF,  DstBuffer, Size);
					Block_Encode_WriteNybble(v>>4, DstBuffer, Size);
					Block_Encode_WriteNybble(v,    DstBuffer, Size);
				}

				//! Skip the zeros
				NextCodedIdx += n;
				zR           -= n;
			}

			//! Insert coded coefficients
			//! NOTE:
			//!  We might still have more coefficients marked for skipping,
			//!  but this didn't take into account the actual statistics of
			//!  the coded zero runs. This means that the coefficients might
			//!  actually not collapse to 0, so we may as well code them anyway
			//!  as it would cost the same either way (though they might quantize
			//!  sub-optimally from not being considered originally)
			do {
				//! Get quantized coefficient
				//! -7h..+7h
				int Qn = Block_Encode_Quantize(Coef[NextCodedIdx], q);
				if(Qn < -7) Qn = -7;
				if(Qn > +7) Qn = +7;

				//! Write to output
				Block_Encode_WriteNybble(Qn, DstBuffer, Size);
			} while(++NextCodedIdx <= CurIdx);
		}

		//! Move to the next coefficient
		do CurIdx++; while(CurIdx < EndIdx && CoefIdx[CurIdx] >= nOutCoef);
	} while(CurIdx < EndIdx);
	return NextCodedIdx;
}
static inline int Block_Encode_EncodePass(const struct ULC_EncoderState_t *State, uint8_t *DstBuffer, int nOutCoef) {
	int BlockSize   = State->BlockSize;
	int Chan, nChan = State->nChan;
	const float *Coef    = State->TransformBuffer;
	const float *CoefNp  = State->TransformNepers;
	const int   *CoefIdx = State->TransformIndex;

	//! Begin coding
	int Idx  = 0;
	int Size = 0; //! Block size (in bits)
	Block_Encode_WriteNybble(State->ThisOverlap, &DstBuffer, &Size);
	for(Chan=0;Chan<nChan;Chan++) {
		int   NextCodedIdx  = Idx;
		int   PrevQuant     = -1;
		int   ChanLastIdx   = Idx + BlockSize;
		int   QuantStartIdx = -1;
		float QuantMin      =  1000.0f; //! Start out of range to force a reset
		float QuantMax      = -1000.0f;
		float QuantSum      = 0.0f;
		float QuantWeight   = 0.0f;
#define WRITE_QUANT_ZONE() \
	Block_Encode_EncodePass_WriteQuantizerZone(QuantStartIdx, Idx, QuantSum, QuantWeight, Coef, CoefIdx, NextCodedIdx, &PrevQuant, &DstBuffer, &Size, nOutCoef)
		for(;;Idx++) {
			//! Seek the next coefficient
			while(Idx < ChanLastIdx && CoefIdx[Idx] >= nOutCoef) Idx++;
			if(Idx >= ChanLastIdx) break;

			//! Level out of range in this quantizer zone?
			const float MaxRange2 = SQR(SQR(7.5f)); //! Maximum quantized range before needing a split, squared
			float BandCoef2  = SQR(Coef  [Idx]);
			float BandCoefNp =    (CoefNp[Idx]);
			if(QuantStartIdx == -1) QuantStartIdx = Idx;
			if(BandCoef2 < QuantMin) QuantMin = BandCoef2;
			if(BandCoef2 > QuantMax) QuantMax = BandCoef2;
			if(QuantMax > QuantMin*MaxRange2) { //! sqrtf(QuantMax/QuantMin) > MaxRange
				//! Write the quantizer zone we just searched through
				//! and start a new one from this coefficient
				NextCodedIdx  = WRITE_QUANT_ZONE();
				QuantStartIdx = Idx;
				QuantMin    = BandCoef2;
				QuantMax    = BandCoef2;
				QuantSum    = BandCoef2 * BandCoefNp;
				QuantWeight = BandCoef2;
			} else {
				//! Accumulate to the current quantizer zone
				QuantSum    += BandCoef2 * BandCoefNp;
				QuantWeight += BandCoef2;
			}
		}

		//! If there's anything in the last quantizer, we must write that quantizer zone
		if(QuantWeight != 0.0f) NextCodedIdx = WRITE_QUANT_ZONE();

		//! 8h,0h,Fh: Stop
		//! If we're at the edge of the block, it might work better to just fill with 0h
		int n = ChanLastIdx - NextCodedIdx;
		if(n > 3) {
			//! If we coded anything, then we must specify the lead sequence
			if(PrevQuant != -1) {
				Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
				Block_Encode_WriteNybble(0x0, &DstBuffer, &Size);
			}
			Block_Encode_WriteNybble(0xF, &DstBuffer, &Size);
		} else if(n > 0) {
			do Block_Encode_WriteNybble(0x0, &DstBuffer, &Size); while(--n);
		}
#undef WRITE_QUANT_ZONE
	}

	//! Pad final byte as needed
	if(Size % 8u) Block_Encode_WriteNybble(0x0, &DstBuffer, &Size);
	return Size;
}
static int Block_Encode(const struct ULC_EncoderState_t *State, uint8_t *DstBuffer, int nNzCoef, float RateKbps) {
	int Size;
	int nOutCoef  = -1;
	int BitBudget = (int)((State->BlockSize * RateKbps) * 1000.0f/State->RateHz); //! NOTE: Truncate

	//! Perform a binary search for the optimal nOutCoef
	int Lo = 0, Hi = nNzCoef;
	if(Lo < Hi) {
		do {
			nOutCoef = (Lo + Hi) / 2u;
			Size = Block_Encode_EncodePass(State, DstBuffer, nOutCoef);
			     if(Size < BitBudget) Lo = nOutCoef;
			else if(Size > BitBudget) Hi = nOutCoef-1;
			else {
				//! Should very, VERY rarely happen, but just in case
				Lo = nOutCoef;
				break;
			}
		} while(Lo < Hi-1);
	}

	//! Avoid going over budget
	int nOutCoefFinal = Lo;
	if(nOutCoefFinal != nOutCoef) Size = Block_Encode_EncodePass(State, DstBuffer, nOutCoef = nOutCoefFinal);
	return Size;
}

/**************************************/
//! EOF
/**************************************/
