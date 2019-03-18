/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
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
	if(nChan     < MIN_CHANS || nChan     > MAX_CHANS) return -1;
	if(BlockSize < MIN_BANDS || BlockSize > MAX_BANDS) return -1;

	//! Get buffer offsets+sizes
	//! PONDER: As with the encoder, this... is probably not ideal
	size_t TransformBuffer_Size  = sizeof(float)  * (nChan* BlockSize   );
	size_t TransformTemp_Size    = sizeof(float)  * (       BlockSize   );
	size_t TransformInvLap_Size  = sizeof(float)  * (nChan*(BlockSize/2));
	size_t _TransformInvLap_Size = sizeof(float*) * (nChan              );
	size_t TransformBuffer_Offs  = 0;
	size_t TransformTemp_Offs    = TransformBuffer_Offs + TransformBuffer_Size;
	size_t TransformInvLap_Offs  = TransformTemp_Offs   + TransformTemp_Size;
	size_t _TransformInvLap_Offs = TransformInvLap_Offs + TransformInvLap_Size;
	size_t AllocSize             = _TransformInvLap_Offs + _TransformInvLap_Size;

	//! Allocate buffer space
	char *Buf = State->BufferData = malloc(BUFFER_ALIGNMENT-1 + AllocSize);
	if(!Buf) return -1;

	//! Initialize pointers
	size_t i, Chan;
	Buf += -(uintptr_t)Buf % BUFFER_ALIGNMENT;
	State->TransformBuffer = (float  *)(Buf + TransformBuffer_Offs);
	State->TransformTemp   = (float  *)(Buf + TransformTemp_Offs);
	State->TransformInvLap = (float **)(Buf + _TransformInvLap_Offs);
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
	size_t nChan            = State->nChan;
	size_t BlockSize        = State->BlockSize;
	float *TransformBuffer  = State->TransformBuffer;
	float *TransformTemp    = State->TransformTemp;
	float **TransformInvLap = State->TransformInvLap;

	size_t Chan;
	size_t Size = 0;
	for(Chan=0;Chan<nChan;Chan++) {
		//! Start decoding coefficients
		int32_t v;
		float *CoefDst = TransformBuffer;
		size_t CoefRem = BlockSize;
		uint8_t Quant  = Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF;
		if(Quant == 0xF) {
			//! [8h,0h,]Fh: Stop
			do *CoefDst++ = 0.0f; while(--CoefRem);
		} else for(;;) {
			//! -7h..+7h: Normal
			v = ((int32_t)Block_Decode_ReadNybble(&SrcBuffer, &Size) << 28) >> 28;
			if(v != -0x8) {
				//! Store dequantized
				*CoefDst++ = (float)(v << Quant);
				if(--CoefRem == 0) break;
				continue;
			}

			//! Unpack escape code
			v = Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF;
			if(v != 0x0) {
				//! 8h,1h..Bh:     4.. 24 zeros
				//! 8h,Ch..Fh,Xh: 26..152 zeros
				size_t nZ = v;
				if(nZ < 0xC) nZ = nZ*2 + 2;
				else {
					nZ = (nZ-0xC)<<4 | (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF);
					nZ = nZ*2 + 26;
				}

				//! Clipping to avoid buffer overflow on corrupted blocks
				if(nZ > CoefRem) nZ = CoefRem;

				//! Insert zeros
				CoefRem -= nZ;
				do *CoefDst++ = 0.0f; while(--nZ);
				if(CoefRem == 0) break;
			} else {
				//! 8h,0h,0h..Eh: Quantizer change
				//! 8h,0h,Fh:     Stop
				Quant = Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF;
				if(Quant == 0xF) {
					do *CoefDst++ = 0.0f; while(--CoefRem);
					break;
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
