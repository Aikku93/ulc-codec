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
	size_t nChan             = State->nChan;
	size_t nQuants           = State->nQuants;
	float **TransformBuffer  = State->TransformBuffer;
	const uint16_t *QuantsBw = State->QuantsBw;
	int16_t **Quants         = State->Quants;
	struct AnalysisKey_t *Keys = State->AnalysisKeys;

	//! Sort keys by band index
	//! This avoids a search for the next non-zero band
	//! Also because the channel is coded in the high bits, we can
	//! code one channel at a time, making things easier as well
	Analysis_KeysSort(Keys, nNzBands);

	//! Start coding
	size_t Chan, QBand;
	size_t Key      = 0;
	size_t Size     = 0; //! Block size (in bits)
	size_t nNzCoded = 0; //! Coded non-zero coefficients
	for(Chan=0;Chan<nChan;Chan++) {
		//! Code the quantizer values (in log2 form)
		for(QBand=0;QBand<nQuants;QBand++) {
			size_t s = 31 - __builtin_clz(Quants[Chan][QBand]);
			Block_Encode_WriteNybble(s, &DstBuffer, &Size);
		}

		//! Start coding coefficients
		size_t NextNz, LastNz = 0;
		size_t NxtQuantBand = 0;
		for(;;) {
			//! Skip unused quantizer bands
			while(NxtQuantBand < nQuants && Quants[Chan][NxtQuantBand] == QUANTIZER_UNUSED) {
				LastNz += QuantsBw[NxtQuantBand++];
			}
			if(NxtQuantBand >= nQuants) break;

			//! Set limit for the /current/ quantizer band
			NextNz  = LastNz;
			LastNz += QuantsBw[NxtQuantBand];
			size_t CurQuantBand = NxtQuantBand;
			size_t CurQuantEnd  = LastNz;

			//! Set limit for coefficients, taking into account consecutive quantizer bands
			while(++NxtQuantBand < nQuants && Quants[Chan][NxtQuantBand] != QUANTIZER_UNUSED) LastNz += QuantsBw[NxtQuantBand];

			//! Code the coefficients
			do {
				//! Unpack key data
				//! If we cross to the next [coded] quantizer band or channel, break out
				size_t tBand = Keys[Key].Band; if(tBand >= LastNz) break;
				size_t tChan = Keys[Key].Chan; if(tChan != Chan)   break;

				//! Code the zero runs
				//! NOTE: Always in multiples of two
				size_t zR = (tBand - NextNz) / 2;
				if(zR) do {
					//! Small run?
					size_t n = zR;
					if(n < 26/2) {
						//! 8h,0h..Bh: 2..24 zeros
						n = n-2/2;
						Block_Encode_WriteNybble(0x8, &DstBuffer, &Size);
						Block_Encode_WriteNybble(n,   &DstBuffer, &Size);
						n = n+2/2;
					} else {
						//! 8h,Ch..Fh,Xh: 26..152 zeros (Ch + n>>4, n&Fh)
						n = n-26/2; if(n > 0x3F) n = 0x3F;
						Block_Encode_WriteNybble(0x8,          &DstBuffer, &Size);
						Block_Encode_WriteNybble(0xC + (n>>4), &DstBuffer, &Size);
						Block_Encode_WriteNybble(n&0xF,        &DstBuffer, &Size);
						n = n+26/2;
					}

					//! Insert zeros
					NextNz += n*2;
					zR     -= n;
				} while(zR);

				//! Insert coded coefficients
				//! NOTE:
				//!  We might still have one more coefficient marked for skipping
				//!  but this didn't take into account the actual statistics of
				//!  the coded zero runs. This means that the coefficient might
				//!  actually not collapse to 0, so we may as well code it anyway
				//!  as it would cost the same either way (though it might quantize
				//!  sub-optimally from not being considered originally)
				do {
					//! Crossed to the next quantizer band?
					//! NOTE: Can only cross one quantizer band at a time, or
					//!       that quantizer band would've been disabled; this
					//!       simplifies the following into an if() statement
					//!       rather than a while() loop
					if(NextNz >= CurQuantEnd) CurQuantEnd += QuantsBw[++CurQuantBand];

					//! Get quantized coefficient
					//! -7h..+7h
					int32_t Qn = (int32_t)round(TransformBuffer[Chan][NextNz] / Quants[Chan][CurQuantBand]);
					if(Qn < -7) Qn = -7;
					if(Qn > +7) Qn = +7;

					//! Write to output
					Block_Encode_WriteNybble(Qn, &DstBuffer, &Size);
					if(Qn != 0) nNzCoded++;
				} while(++NextNz <= tBand);
			} while(++Key < nNzBands);

			//! Finalize the block (0h,0h: Stop)
			//! If we're at the edge of the block, we may as well terminate
			//! normally, though, so just cap to max 2x 0h values
			size_t n = LastNz - NextNz; if(n > 2) n = 2;
			while(n--) Block_Encode_WriteNybble(0x0, &DstBuffer, &Size);
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
