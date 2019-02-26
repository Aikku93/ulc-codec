/**************************************/
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
/**************************************/
#include "Fourier.h"
#include "ulcDecoder.h"
#include "ulcUtility.h"
/**************************************/
#if defined(__AVX__)
# define BUFFER_ALIGNMENT 32u //! __mm256
#elif defined(__SSE__)
# define BUFFER_ALIGNMENT 16u //! __mm128
#else
# define BUFFER_ALIGNMENT 4u //! float
#endif
/**************************************/

//! Just for consistency
#define MIN_BANDS 8
#define MAX_BANDS 65536
#define MIN_CHANS 1
#define MAX_CHANS 65536

/**************************************/

//! Initialize decoder state
int ULC_DecoderState_Init(struct ULC_DecoderState_t *State) {
	//! Clear anything that is needed for EncoderState_Destroy()
	State->BufferData = NULL;

	//! Verify parameters
	size_t nChan     = State->nChan;
	size_t BlockSize = State->BlockSize;
	size_t nQuants   = State->nQuants;
	if(nChan     < MIN_CHANS || nChan     > MAX_CHANS) return -1;
	if(BlockSize < MIN_BANDS || BlockSize > MAX_BANDS) return -1;

	//! Get buffer offsets+sizes
	//! PONDER: As with the encoder, this... is probably not ideal
	size_t TransformBuffer_Size  = sizeof(float)    * (nChan* BlockSize   );
	size_t TransformTemp_Size    = sizeof(float)    * (       BlockSize   );
	size_t TransformInvLap_Size  = sizeof(float)    * (nChan*(BlockSize/2));
	size_t Quants_Size           = sizeof(uint8_t)  * (nQuants            );
	size_t _TransformInvLap_Size = sizeof(float*)   * (nChan              );
	size_t TransformBuffer_Offs  = 0;
	size_t TransformTemp_Offs    = TransformBuffer_Offs  + TransformBuffer_Size;
	size_t TransformInvLap_Offs  = TransformTemp_Offs    + TransformTemp_Size;
	size_t Quants_Offs           = TransformInvLap_Offs  + TransformInvLap_Size;
	size_t _TransformInvLap_Offs = Quants_Offs           + Quants_Size;
	size_t AllocSize = _TransformInvLap_Offs + _TransformInvLap_Size;

	//! Allocate buffer space
	char *Buf = State->BufferData = malloc(BUFFER_ALIGNMENT-1 + AllocSize);
	if(!Buf) return -1;

	//! Initialize pointers
	size_t i, Chan;
	Buf += -(uintptr_t)Buf % BUFFER_ALIGNMENT;
	State->TransformBuffer = (float  *)(Buf + TransformBuffer_Offs);
	State->TransformTemp   = (float  *)(Buf + TransformTemp_Offs);
	State->TransformInvLap = (float **)(Buf + _TransformInvLap_Offs);
	State->Quants          = (uint8_t*)(Buf + Quants_Offs);
	for(Chan=0;Chan<nChan;Chan++) {
		State->TransformInvLap[Chan] = (float*)(Buf + TransformInvLap_Offs) + Chan*(BlockSize/2);

		//! Everything can remain uninitialized except for the lapping buffer
		for(i=0;i<BlockSize/2;i++) State->TransformInvLap[Chan][i] = 0.0f;
	}

	//! Success
	return 1;
}

/**************************************/

//! Destroy decoder state
void ULC_DecoderState_Destroy(struct ULC_DecoderState_t *State) {
	//! Free buffer space
	free(State->BufferData);
}

/**************************************/

//! Decode block
static inline uint8_t Block_Decode_ReadNybble(const uint8_t **Src, size_t *Size) {
	//! Fetch and shift nybble
	uint8_t x = *(*Src);
	*Size += 4;
	if((*Size)%8u == 0) x >>= 4, (*Src)++;
	return x;
}
size_t ULC_DecodeBlock(const struct ULC_DecoderState_t *State, float *DstData, const uint8_t *SrcBuffer) {
	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	size_t nChan             = State->nChan;
	size_t BlockSize         = State->BlockSize;
	size_t nQuants           = State->nQuants;
	const uint16_t *QuantsBw = State->QuantsBw;
	float *TransformBuffer   = State->TransformBuffer;
	float *TransformTemp     = State->TransformTemp;
	float **TransformInvLap  = State->TransformInvLap;
	uint8_t *Quants          = State->Quants;

	size_t Chan, QBand;
	size_t Size = 0;
	for(Chan=0;Chan<nChan;Chan++) {
		//! Read quantizers
		for(QBand=0;QBand<nQuants;QBand++) Quants[QBand] = Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF;

		//! Start decoding coefficients
		size_t NextQuantBand = 0;
		float *CoefDst = TransformBuffer;
		for(;;) {
			size_t nZ;

			//! Insert zeros in unused quantizer bands
			nZ = 0;
			while(NextQuantBand < nQuants && Quants[NextQuantBand] == 0xF) nZ += QuantsBw[NextQuantBand++];
			if(nZ) do *CoefDst++ = 0.0f; while(--nZ);
			if(NextQuantBand >= nQuants) break;

			//! Get bandwidth for current quantizer band (to keep track of when
			//! it changes, for consecutive enabled quantizer bands) as well as
			//! the bandwidth until the last consecutive enabled quantizer band
			size_t CurQuantBand    = NextQuantBand;
			size_t MaxQuantBandBw  = QuantsBw[CurQuantBand];
			float *CurQuantBandEnd = CoefDst + MaxQuantBandBw;
			while(++NextQuantBand < nQuants && Quants[NextQuantBand] != 0xF) MaxQuantBandBw += QuantsBw[NextQuantBand];

			//! Read the coefficients
			int32_t LastV = 1; //! Anything non-zero
			for(;;) {
				int32_t v = ((int32_t)Block_Decode_ReadNybble(&SrcBuffer, &Size) << 28) >> 28;

				//! Normal coefficient? (-7..+7)
				if(v != -0x8) {
					//! If we have two zeros in a row, this is the stop code 0h,0h
					if((v|LastV) == 0) {
						do *CoefDst++ = 0.0f; while(--MaxQuantBandBw);
						break;
					}
					LastV = v;

					//! Crossed to the next quantizer band?
					//! NOTE: Can only cross one quantizer band at a time, or
					//!       that quantizer band would've been disabled; this
					//!       simplifies the following into an if() statement
					//!       rather than a while() loop
					if(CoefDst >= CurQuantBandEnd) CurQuantBandEnd += QuantsBw[++CurQuantBand];

					//! Store dequantized
					*CoefDst++ = (float)(v << Quants[CurQuantBand]);
					if(--MaxQuantBandBw == 0) break;
					continue;
				} else {
					//! Unpack zero run
					nZ = Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF;
					if(nZ < 0xC) {
						//! Small run
						//! 0h..Bh: 2..24 zeros
						nZ = nZ*2 + 2;
					} else {
						//! Long run
						//! Ch..Fh,Xh: 26..152 zeros (Ch + n>>4, n&Fh)
						nZ = (nZ-0xC)<<4 | (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF);
						nZ = nZ*2 + 26;
					}

					//! Prevent buffer uverflow
					//! PONDER: Maybe just report the block as broken?
					if(nZ > MaxQuantBandBw) {
						nZ = MaxQuantBandBw;
					}

					//! Insert zeros
					MaxQuantBandBw -= nZ;
					do *CoefDst++ = 0.0f; while(--nZ);
					if(!MaxQuantBandBw) break;
				}
			}
		}

		//! Inverse transform block
		ULC_Transform_AntiPreEcho(TransformBuffer, BlockSize);
		Fourier_IMDCT(DstData + Chan*BlockSize, TransformBuffer, TransformInvLap[Chan], TransformTemp, BlockSize);
	}
	return Size;
}

/**************************************/
//! EOF
/**************************************/
