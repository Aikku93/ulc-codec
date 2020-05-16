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
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_Quantizer.h"
/**************************************/

//! Returns the block size (in bits) and the number of coded (non-zero) coefficients
static inline int Block_Encode_Quantize(float v, float q) {
	float av = ABS(v);
	int vq = (int)(sqrtf(av*q) + 0.5f);
#if 1 //! Optimal rounding
	//! NOTE: Handle vq+1 first in case vq==0, as dl and dr would be equal, and
	//! so if vq-1 was handled first, it would set vq=-1, which would sign-flip
	float dl = ABS(av*q - SQR(vq-1)); int vql = vq-1;
	float d  = ABS(av*q - SQR(vq  ));
	float dr = ABS(av*q - SQR(vq+1)); int vqr = vq+1;
	if(dr < d) vq = vqr, d = dr;
	if(dl < d) vq = vql;
#endif
	return (v < 0.0f) ? (-vq) : (+vq);
}
static inline void Block_Encode_WriteNybble(uint8_t x, uint8_t **Dst, int *Size) {
	//! Push nybble
	*(*Dst) >>= 4;
	*(*Dst)  |= x << 4;
	*Size += 4;

	//! Next byte?
	if((*Size)%8u == 0) (*Dst)++;
}
static inline void Block_Encode_WriteQuantizer(float Quant, uint8_t **DstBuffer, int *Size, int Lead) {
	//! 8h,0h,0h..Eh[,0h..Ch]: Quantizer change
	int s = (int)log2f(Quant) - 5;
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
static int Block_Encode(const struct ULC_EncoderState_t *State, uint8_t *DstBuffer, int nKeys, int *_nNzCoded) {
	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	int     nChan           = State->nChan;
	int     BlockSize       = State->BlockSize;
	float **TransformBuffer = State->TransformBuffer;
	struct AnalysisKey_t *Keys = State->AnalysisKeys;

	//! Group the coefficients into quantizer zones
	Block_Encode_ProcessQuantizerZones(State, nKeys);

	//! Start coding
	int Chan;
	int Key      = 0;
	int Size     = 0; //! Block size (in bits)
	int nNzCoded = 0; //! Coded non-zero coefficients
	Block_Encode_WriteNybble(State->ThisOverlap, &DstBuffer, &Size);
	for(Chan=0;Chan<nChan;Chan++) {
		//! Code the coefficients
		//! NOTE:
		//!  Check Key < nKeys: If the last channel is silent, this avoids trying to code it by accident
		//!  Check Key.Chan==Chan: This correctly codes silent channels
		if(Key < nKeys && Keys[Key].Chan == Chan) {
			//! Code the first quantizer ([8h,0h,]0h..Eh) and start coding
			int   NextBand  = 0;
			float LastQuant = Keys[Key].Quant;
			Block_Encode_WriteQuantizer(LastQuant, &DstBuffer, &Size, 0);
			do {
				//! Unpack key data
				//! If we cross to the next channel, break out
				int tChan = Keys[Key].Chan; if(tChan != Chan) break;
				int tBand = Keys[Key].Band;

				//! Update quantizer?
				if(Keys[Key].Quant != 0.0f && Keys[Key].Quant != LastQuant) {
					LastQuant = Keys[Key].Quant;
					Block_Encode_WriteQuantizer(LastQuant, &DstBuffer, &Size, 1);
				}

				//! As we don't do an optimization pass over the coefficients, this
				//! key might collapse to 0. So check for this first, as if it does
				//! collapse, we can extend the zero run further for bit savings
				if(ABS(Block_Encode_Quantize(TransformBuffer[Chan][tBand], LastQuant)) == 0) continue;

				//! Code the zero runs
				//! NOTE: Escape-code-coded zero runs have a minimum size of 4 coefficients
				//!       This is because two zero coefficients can be coded as 0h,0h, so
				//!       we instead use 8h,0h for the 'stop' code. This also allows coding
				//!       of some coefficients we may have missed (see below)
				int zR = tBand - NextBand;
				while(zR >= 4) {
					//! Encode optimal run
					int n = zR;
					if(n < 28+2) {
						//! 8h,1h..Dh:  4.. 28 zeros (Step: 2)
						n = (n-2)/2u;
						Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
						Block_Encode_WriteNybble(n,   &DstBuffer, &Size);
						n = n*2+2;
					} else if(n < 90+4) {
						//! 8h,Eh,Xh:  30.. 90 zeros (Step: 4)
						n = (n-30)/4u;
						Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
						Block_Encode_WriteNybble(0xE, &DstBuffer, &Size);
						Block_Encode_WriteNybble(n,   &DstBuffer, &Size);
						n = n*4+30;
					} else {
						//! 8h,Fh,Xh:  94..214 zeros (Step: 8)
						n = (n-94)/8u; if(n > 0xF) n = 0xF;
						Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
						Block_Encode_WriteNybble(0xF, &DstBuffer, &Size);
						Block_Encode_WriteNybble(n,   &DstBuffer, &Size);
						n = n*8+94;
					}

					//! Insert zeros
					NextBand += n;
					zR       -= n;
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
					int Qn = Block_Encode_Quantize(TransformBuffer[Chan][NextBand], LastQuant);
					if(Qn < -7) Qn = -7;
					if(Qn > +7) Qn = +7;

					//! Write to output
					Block_Encode_WriteNybble(Qn, &DstBuffer, &Size);
				} while(++NextBand <= tBand);
				nNzCoded++;
			} while(++Key < nKeys);

			//! 8h,0h,Fh: Stop
			//! If we're at the edge of the block, it might work better to just fill with 0h
			int n = BlockSize - NextBand;
			if(n > 3) {
				Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
				Block_Encode_WriteNybble(0x0, &DstBuffer, &Size);
				Block_Encode_WriteNybble(0xF, &DstBuffer, &Size);
			} else if(n > 0) {
				do Block_Encode_WriteNybble(0x0, &DstBuffer, &Size); while(--n);
			}
		} else {
			//! [8h,0h,]Fh: Stop
			Block_Encode_WriteNybble(0xF, &DstBuffer, &Size);
		}
	}

	//! Shift down final byte if needed
	if(Size % 8u) *DstBuffer >>= 4;

	//! Return output size
	if(_nNzCoded) *_nNzCoded = nNzCoded;
	return Size;
}

/**************************************/
//! EOF
/**************************************/
