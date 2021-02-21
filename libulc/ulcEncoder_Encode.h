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
#if ULC_USE_NOISE_CODING
static int Block_Encode_EncodePass_GetNoiseQ(float q, const float *Coef, const float *CoefNp, int Band, int n) {
	//! Get the energy and spectral flatness of the run
	//! NOTE: In testing, it was found necessary to map
	//! Flatness through a fairly steep curve, and then
	//! adjust the overall gain to compensate. Without
	//! this, noise-fill is too loud in sections where
	//! it should be a lot more subtle (and conversely,
	//! too quiet when it should've been louder).
	float Sum = 0.0f, SumW = 0.0f;
	int BandEnd = Band+n;
	do {
		float c   = SQR(Coef  [Band]);
		float cNp =    (CoefNp[Band]); //! <- Scaled by 0.5 relative to Log[Coef^2]
		Sum  += c * cNp;
		SumW += c;
	} while(++Band < BandEnd);
	if(Sum == 0.0f) return 0; //! Empty
	float Flatness  = (logf(SumW) - 2.0f*Sum/SumW) / logf(n);
	float Amplitude = sqrtf(SumW / n) * SQR(SQR(Flatness)) * 2.0f;

	//! Quantize into final code
	//! NOTE: This is encoded at higher precision, because it spans
	//! the full 4bit range, meaning we have an extra Log2[16^2 / 7^2]
	//! bits to play with (2.385 bits, so use 3.0 for simplicity,
	//! especially since noise should be lower than the maximum value).
	int NoiseQ = Block_Encode_Quantize(Amplitude, q*8.0f);
	if(NoiseQ > 0xF+1) NoiseQ = 0xF+1;
	return NoiseQ;
}
#endif
static inline __attribute__((always_inline)) int Block_Encode_EncodePass_WriteQuantizerZone(
	int       CurIdx,
	int       EndIdx,
	float     QuantSum,
	float     QuantWeight,
	const float *Coef,
	const float *CoefNp,
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
			int zR = CurIdx - NextCodedIdx;
			while(zR >= 3) {
				Block_Encode_WriteNybble(0x8, DstBuffer, Size);

				//! Determine which run type to use and  get the number of zeros coded
				//! NOTE: A short run takes 2 nybbles, and a long run takes 4 nybbles.
				//! So two short runs of maximum length code up to 30 zeros with the
				//! same efficiency as a long run, meaning that long runs start with
				//! 31 zeros.
				int n;
				if(zR < 31) {
					//! 8h,1h..Dh: 3 .. 15 zeros
					n = zR; if(n > 15) n = 15;
					Block_Encode_WriteNybble(n-2, DstBuffer, Size);
				} else {
					//! 8h,Eh,Zh,Yh,Xh: 31 .. 286 noise samples
					//! 8h,Fh,Yh,Xh:    31 .. 286 zeros
					int v = zR-31; if(v > 0xFF) v = 0xFF;
					n = v + 31;

					//! Determine the quantized coefficient for noise-fill mode
					int NoiseQ = 0;
#if ULC_USE_NOISE_CODING
					NoiseQ = Block_Encode_EncodePass_GetNoiseQ(q, Coef, CoefNp, NextCodedIdx, n);
#endif
					if(NoiseQ) {
						Block_Encode_WriteNybble(0xE,      DstBuffer, Size);
						Block_Encode_WriteNybble(v>>4,     DstBuffer, Size);
						Block_Encode_WriteNybble(v,        DstBuffer, Size);
						Block_Encode_WriteNybble(NoiseQ-1, DstBuffer, Size);
					} else {
						Block_Encode_WriteNybble(0xF,  DstBuffer, Size);
						Block_Encode_WriteNybble(v>>4, DstBuffer, Size);
						Block_Encode_WriteNybble(v,    DstBuffer, Size);
					}
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
#if ULC_USE_NOISE_CODING
	const float *CoefNp  = State->TransformNepers;
#endif
	const int   *CoefIdx = State->TransformIndex;

	//! Begin coding
	int Idx  = 0;
	int Size = 0; //! Block size (in bits)
	Block_Encode_WriteNybble(State->WindowCtrl, &DstBuffer, &Size);
	if(State->WindowCtrl & 0x8) Block_Encode_WriteNybble(State->WindowCtrl >> 4, &DstBuffer, &Size);
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
		CoefNp, \
		CoefIdx, \
		NextCodedIdx, \
		&PrevQuant, \
		&DstBuffer, \
		&Size, \
		nOutCoef \
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
		if(n > 3) {
			//! If we coded anything, then we must specify the lead sequence
			if(PrevQuant != -1) {
				Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
				Block_Encode_WriteNybble(0x0, &DstBuffer, &Size);
			}
			Block_Encode_WriteNybble(0xF, &DstBuffer, &Size);

			//! Analyze the remaining data for noise-fill mode
			int NoiseQ = 0, NoiseDecay = 0;
			if(PrevQuant != -1) {
#if ULC_USE_NOISE_CODING
				//! Analyze the energy to form the noise-fill energy
				{
					int   Band;
					float Sum  = 0.0f, SumW = 0.0f;
					for(Band=NextCodedIdx;Band<ChanLastIdx;Band++) {
						float c2  = SQR(Coef  [Band]);
						float cNp =    (CoefNp[Band]);
						Sum += c2 * cNp, SumW += c2;
					}
					if(Sum) {
						//! This uses a similar method to noise-fill for selecting the gain,
						//! but using the smooth-max peak rather than the RMS amplitude,
						//! and a higher compensation gain when accounting for flatness.
						float Flatness = (logf(SumW) - 2.0f*Sum/SumW) / logf(n);
						float Amplitude = expf(Sum/SumW) * SQR(SQR(Flatness)) * 4.0f;
						NoiseQ = (int)sqrtf(Amplitude * (float)(8u << PrevQuant)); //! <- Round down
						if(NoiseQ > 0xF) NoiseQ = 0xF;
					}
				}

				//! If we are going to noise-code, analyze the decay
				//! rate between the start of the fill and the end
				if(NoiseQ) {
					int Band;

					//! Analyze the start and end simultaneously, with a
					//! tapered, overlapping window to capture the envelope.
					//! NOTE: Using Decay=End/Beg follows from the following:
					//!  If we take the average decay between each coefficient
					//!  in the log domain, we end up with:
					//!   Sum[Log[c[n+1]] - Log[c[n]], {n,N}] / N
					//!  However, the sum cancels out completely in the middle
					//!  sections, leading to:
					//!    (Log[c[N]] - Log[c[1]]) / N
					float Beg = 0.0f, End = 0.0f;
					for(Band=NextCodedIdx;Band<ChanLastIdx;Band++) {
						float c2 = SQR(Coef[Band]);
						Beg += c2 * (ChanLastIdx - Band);
						End += c2 * (Band - NextCodedIdx + 1);
					}
					float Decay = powf(End / Beg, 1.0f / (n*2)); //! *2 for square root in the RMS sum

					//! Decay is mapped through an inverted Sqrt curve for higher accuracy
					Decay = (Decay > 1.0f) ? 0.0f : sqrtf(1.0f - Decay);
					NoiseDecay = (int)(Decay*17 - 1 + 0.5f);
					if(NoiseDecay < 0x0) NoiseDecay = 0x0;
					if(NoiseDecay > 0xF) NoiseDecay = 0xF;
				}
#endif
				Block_Encode_WriteNybble(NoiseQ, &DstBuffer, &Size);
				if(NoiseQ) Block_Encode_WriteNybble(NoiseDecay, &DstBuffer, &Size);
			}
		} else if(n > 0) {
			do Block_Encode_WriteNybble(0x0, &DstBuffer, &Size); while(--n);
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
