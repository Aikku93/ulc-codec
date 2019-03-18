/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <stddef.h>
#include <stdint.h>
#include <math.h>
/**************************************/
#include "Fourier.h"
#include "ulcEncoder.h"
#include "ulcUtility.h"
/**************************************/
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_Quantizer.h"
/**************************************/

//! Returns the block size (in bits) and the number of coded (non-zero) coefficients
static inline void Block_Encode_WriteNybble(uint8_t x, uint8_t **Dst, size_t *Size) {
	//! Push nybble
	*(*Dst) >>= 4;
	*(*Dst)  |= x << 4;
	*Size += 4;

	//! Next byte?
	if((*Size)%8u == 0) (*Dst)++;
}
static size_t Block_Encode(const struct ULC_EncoderState_t *State, uint8_t *DstBuffer, size_t nNzMax, size_t nKeys, size_t *_nNzCoded) {
	//! Generate quantizers and get number of non-zero bands
	size_t nNzBands = Block_Encode_BuildQuants(State, nNzMax, nKeys);

	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	size_t nChan            = State->nChan;
	size_t BlockSize        = State->BlockSize;
	float **TransformBuffer = State->TransformBuffer;
	struct AnalysisKey_t *Keys = State->AnalysisKeys;

	//! Sort keys by band index
	//! This avoids a search for the next non-zero band
	//! Also because the channel is coded in the high bits, we can
	//! code one channel at a time, making things easier as well
	size_t nChanLog2     = IntLog2(nChan);
	size_t BlockSizeLog2 = IntLog2(BlockSize);
	if(nChan > (1u<<nChanLog2)) nChanLog2++; //! Round to next power of two
	Analysis_KeysSort(Keys, Keys+nNzBands, 1u << (nChanLog2 + BlockSizeLog2 - 1));

	//! Start coding
	size_t Chan;
	size_t Key      = 0;
	size_t Size     = 0; //! Block size (in bits)
	size_t nNzCoded = 0; //! Coded non-zero coefficients
	for(Chan=0;Chan<nChan;Chan++) {
		//! Code the coefficients
		if(nNzBands) {
			size_t  NextBand  = 0;
			int16_t LastQuant = Keys[Key].Quant;

			//! Code the first quantizer ([8h,0h,]0h..Eh) and start coding
			Block_Encode_WriteNybble(IntLog2(LastQuant), &DstBuffer, &Size);
			do {
				//! Unpack key data
				//! If we cross to the next channel, break out
				size_t tChan = Keys[Key].Key >> BlockSizeLog2; if(tChan != Chan) break;
				size_t tBand = Keys[Key].Key & (BlockSize-1);

				//! Code the zero runs
				//! NOTE: Escape-code-coded zero runs have a minimum size of 4 coefficients
				//!       This is because two zero coefficients can be coded as 0h,0h, so
				//!       we instead use 8h,0h for the 'stop' code. This also allows coding
				//!       of some coefficients we may have missed (see below)
				size_t zR = tBand - NextBand;
				while(zR >= 4) {
					//! Small run?
					size_t n = zR;
					if(n < 26) {
						//! 8h,1h..Bh: 4..24 zeros
						n = (n-2)/2;
						Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
						Block_Encode_WriteNybble(n,   &DstBuffer, &Size);
						n = n*2+2;
					} else {
						//! 8h,Ch..Fh,Xh: 26..152 zeros (Ch + n>>4, n&Fh)
						n = (n-26)/2; if(n > 0x3F) n = 0x3F;
						Block_Encode_WriteNybble(0x8,          &DstBuffer, &Size);
						Block_Encode_WriteNybble(0xC + (n>>4), &DstBuffer, &Size);
						Block_Encode_WriteNybble(n&0xF,        &DstBuffer, &Size);
						n = n*2+26;
					}

					//! Insert zeros
					NextBand += n;
					zR       -= n;
				}

				//! Update quantizer?
				if(Keys[Key].Quant != LastQuant) {
					LastQuant = Keys[Key].Quant;

					//! 8h,0h,0h..Eh: Quantizer change
					size_t s = IntLog2(LastQuant);
					Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
					Block_Encode_WriteNybble(0x0, &DstBuffer, &Size);
					Block_Encode_WriteNybble(s,   &DstBuffer, &Size);
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
					int32_t Qn = (int32_t)round(TransformBuffer[Chan][NextBand] / LastQuant);
					if(Qn < -7) Qn = -7;
					if(Qn > +7) Qn = +7;

					//! Write to output
					Block_Encode_WriteNybble(Qn, &DstBuffer, &Size);
					if(Qn != 0) nNzCoded++;
				} while(++NextBand <= tBand);
			} while(++Key < nNzBands);

			//! 8h,0h,Fh: Stop
			//! If we're at the edge of the block, it might work better to just fill with 0h
			size_t n = BlockSize - NextBand;
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
