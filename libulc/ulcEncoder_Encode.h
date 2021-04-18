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
#if ULC_USE_NOISE_CODING
# include "ulcEncoder_NoiseFill.h"
#endif
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
static inline int Block_Encode_BuildQuantizer(float Sum, float Weight) {
	//! NOTE: No rounding (ie. this truncates); this favours louder bands.
	//! NOTE: `q` will always be greater than 5 due to the bias so
	//! the quantizer code syntax is biased accordingly.
	int q = (int)(5.0f - 0x1.715476p0f*logf(Sum/Weight)); //! 0x1.715476p0 == 1/Ln[2] for change of base
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
	const float *CoefNoise,
	const int   *CoefIdx,
	int       NextCodedIdx,
	int      *PrevQuant,
	uint8_t **DstBuffer,
	int      *Size,
	int       nOutCoef,
	int       WindowCtrl
) {
	(void)CoefNoise;  //! CoefNoise is only used with ULC_USE_NOISE_CODING
	(void)WindowCtrl; //! Ditto WindowCtrl

	//! Write/update the quantizer
	float q; {
		int qi = Block_Encode_BuildQuantizer(QuantSum, QuantWeight);
		q = (float)(1u << qi);
		if(qi != *PrevQuant) {
			Block_Encode_WriteQuantizer(qi, DstBuffer, Size, *PrevQuant != -1);
			*PrevQuant = qi;
		}
	}

	//! Write the coefficients
	do {
		//! Target coefficient doesn't collapse to 0?
		if(ABS(Coef[CurIdx]*q) >= SQR(0.5f)) { //! Block_Encode_Quantize(Coef[CurIdx], q) != 0
			//! Code the zero runs
			int n, v, zR = CurIdx - NextCodedIdx;
			while(zR >= 3) {
				Block_Encode_WriteNybble(0x8, DstBuffer, Size);

				//! Determine the quantized coefficient for noise-fill mode
#if ULC_USE_NOISE_CODING
				//! 8h,Eh,Zh,Yh,Xh: 31 .. 286 noise samples
				//! NOTE: The range of noise samples is coded as a trade-off
				//! between noise-fill commands and coded coefficients, as
				//! noise-fill competes with coefficients we /want/ to code
				//! for the bit bandwidth available. Setting the start range
				//! too low favours noise-fill over psychoacoustically
				//! significant coefficients, but setting it too high can
				//! cause degradation in small block sizes as well as when
				//! we get many short runs.
				//! TODO: Check the target coefficient's value; if it's
				//! very close to the noise amplitude, combine it into
				//! the noise run?
				//! TODO: Is there a better min-number-of-coefficients?
				//! Too low and this fights too much with the coefficients
				//! we are actually trying to code (reducing the coding
				//! rate of them, and increasing the coding rate of noise).
				//! However, with small BlockSize (or SubBlockSize when
				//! decimating), smaller runs are useful.
				int NoiseQ = 0;
				if(zR >= 31) {
					v = zR - 31; if(v > 0xFF) v = 0xFF;
					n = v + 31;
					NoiseQ = Block_Encode_EncodePass_GetNoiseQ(CoefNoise, q, NextCodedIdx, n, WindowCtrl);
				}
				if(NoiseQ) {
					Block_Encode_WriteNybble(0xE,      DstBuffer, Size);
					Block_Encode_WriteNybble(v>>4,     DstBuffer, Size);
					Block_Encode_WriteNybble(v,        DstBuffer, Size);
					Block_Encode_WriteNybble(NoiseQ-1, DstBuffer, Size);
				} else {
#endif
					//! Determine which run type to use and get the number of zeros coded
					//! NOTE: A short run takes 2 nybbles, and a long run takes 4 nybbles.
					//! So two short runs of maximum length code up to 30 zeros with the
					//! same efficiency as a long run, meaning that long runs start with
					//! 31 zeros.
					if(zR < 31) {
						//! 8h,1h..Dh: 3 .. 15 zeros
						n = zR; if(n > 15) n = 15;
						Block_Encode_WriteNybble(n-2, DstBuffer, Size);
					} else {
						//! 8h,Fh,Yh,Xh: 31 .. 286 zeros
						v = zR-31; if(v > 0xFF) v = 0xFF;
						n = v + 31;
						Block_Encode_WriteNybble(0xF,  DstBuffer, Size);
						Block_Encode_WriteNybble(v>>4, DstBuffer, Size);
						Block_Encode_WriteNybble(v,    DstBuffer, Size);
					}
#if ULC_USE_NOISE_CODING
				}
#endif
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
	const float *Coef      = State->TransformBuffer;
	const float *CoefNoise = State->TransformBuffer;
	const int   *CoefIdx   = State->TransformIndex;
#if ULC_USE_NOISE_CODING
	CoefNoise = State->TransformNoise;
#endif

	//! Begin coding
	int Idx  = 0;
	int Size = 0; //! Block size (in bits)
	int WindowCtrl = State->WindowCtrl;
	Block_Encode_WriteNybble(WindowCtrl, &DstBuffer, &Size);
	if(WindowCtrl & 0x8) Block_Encode_WriteNybble(WindowCtrl >> 4, &DstBuffer, &Size);
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
	Block_Encode_EncodePass_WriteQuantizerZone( \
		QuantStartIdx, \
		Idx, \
		QuantSum, \
		QuantWeight, \
		Coef, \
		CoefNoise, \
		CoefIdx, \
		NextCodedIdx, \
		&PrevQuant, \
		&DstBuffer, \
		&Size, \
		nOutCoef, \
		WindowCtrl \
	)
		for(;;Idx++) {
			//! Seek the next coefficient
			while(Idx < ChanLastIdx && CoefIdx[Idx] >= nOutCoef) Idx++;
			if(Idx >= ChanLastIdx) break;

			//! Level out of range in this quantizer zone?
			//! NOTE: There is a certain amount of overlap
			//! between quantizers:
			//!  7^2 * 2^-1 == x^2 * 2^0
			//!  x == Sqrt[7^2 * 2^-1]
			//! Since 7.0 is the maximum value that goes into a
			//! quantizer, x^2 therefore becomes the usable range, ie.:
			//!  Range: 7^2 * 2^-1
			//! This is approximately 27.8dB.
			//! NOTE: Using 8^2*2^-1 (30.1dB) was tested and proved inferior.
			const float MaxRange = 24.5f;
			float BandCoef   = ABS(Coef[Idx]);
			float BandCoef2  = SQR(Coef[Idx]);
			if(QuantStartIdx == -1) QuantStartIdx = Idx;
			if(BandCoef < QuantMin) QuantMin = BandCoef;
			if(BandCoef > QuantMax) QuantMax = BandCoef;
			if(QuantMax > MaxRange*QuantMin) {
				//! Write the quantizer zone we just searched through
				//! and start a new one from this coefficient
				NextCodedIdx  = WRITE_QUANT_ZONE();
				QuantStartIdx = Idx;
				QuantMin    = BandCoef;
				QuantMax    = BandCoef;
				QuantSum    = BandCoef2;
				QuantWeight = BandCoef;
			} else {
				//! Accumulate to the current quantizer zone
				QuantSum    += BandCoef2;
				QuantWeight += BandCoef;
			}
		}

		//! If there's anything in the last quantizer, we must write that quantizer zone
		if(QuantWeight != 0.0f) NextCodedIdx = WRITE_QUANT_ZONE();

		//! 8h,0h,Fh,0h:        Stop
		//! 8h,0h,Fh,1h..Fh,Xh: Noise-fill
		//! If we're at the edge of the block, it might work better to just fill with 0h
		int n = ChanLastIdx - NextCodedIdx;
		if(n > 4) {
			//! If we coded anything, then we must specify the lead sequence
			if(PrevQuant != -1) {
				Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
				Block_Encode_WriteNybble(0x0, &DstBuffer, &Size);
			}
			Block_Encode_WriteNybble(0xF, &DstBuffer, &Size);

			//! Analyze the remaining data for noise-fill mode
			if(PrevQuant != -1) {
				int NoiseQ = 0, NoiseDecay = 0;
#if ULC_USE_NOISE_CODING
				Block_Encode_EncodePass_GetHFExtParams(
					CoefNoise,
					(float)(8u << PrevQuant),
					NextCodedIdx,
					n,
					WindowCtrl,
					&NoiseQ,
					&NoiseDecay
				);
#endif
				Block_Encode_WriteNybble(NoiseQ, &DstBuffer, &Size);
				if(NoiseQ) Block_Encode_WriteNybble(NoiseDecay, &DstBuffer, &Size);
			}
		} else if(n > 0) {
			//! If n < 4, then PrevQuant can't be invalid because
			//! we must have coded /something/ to get to this point.
			float q = (float)(1u << PrevQuant);
			do {
				//! At this point, it's more efficient to write coefficients
				//! than a stop code. It's unlikely we'll have any actually
				//! meaningful data here, but we may as well code it anyway.
				int Qn = Block_Encode_Quantize(Coef[NextCodedIdx], q);
				if(Qn < -7) Qn = -7;
				if(Qn > +7) Qn = +7;
				Block_Encode_WriteNybble(Qn, &DstBuffer, &Size);
			} while(++NextCodedIdx < ChanLastIdx);
		}
#undef WRITE_QUANT_ZONE
	}

	//! Pad final byte as needed
	if(Size % 8u) Block_Encode_WriteNybble(0x0, &DstBuffer, &Size);
	return Size;
}

/**************************************/
//! EOF
/**************************************/
