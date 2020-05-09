/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2020, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
/**************************************/
#include "Fourier.h"
#include "ulcDecoder.h"
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
#define MAX_BANDS 65535
#define MIN_CHANS 1
#define MAX_CHANS 255

/**************************************/

//! Initialize decoder state
int ULC_DecoderState_Init(struct ULC_DecoderState_t *State) {
	//! Clear anything that is needed for EncoderState_Destroy()
	State->BufferData  = NULL;
	State->NextOverlap = 0;

	//! Verify parameters
	int nChan     = State->nChan;
	int BlockSize = State->BlockSize;
	if(nChan     < MIN_CHANS || nChan     > MAX_CHANS) return -1;
	if(BlockSize < MIN_BANDS || BlockSize > MAX_BANDS) return -1;
	if((BlockSize & (-BlockSize)) != BlockSize)        return -1;

	//! Get buffer offsets+sizes
	//! PONDER: As with the encoder, this... is probably not ideal
	int TransformBuffer_Size  = sizeof(float)  * (nChan* BlockSize   );
	int TransformTemp_Size    = sizeof(float)  * (       BlockSize   );
	int TransformInvLap_Size  = sizeof(float)  * (nChan*(BlockSize/2));
	int _TransformInvLap_Size = sizeof(float*) * (nChan              );
	int TransformBuffer_Offs  = 0;
	int TransformTemp_Offs    = TransformBuffer_Offs + TransformBuffer_Size;
	int TransformInvLap_Offs  = TransformTemp_Offs   + TransformTemp_Size;
	int _TransformInvLap_Offs = TransformInvLap_Offs + TransformInvLap_Size;
	int AllocSize             = _TransformInvLap_Offs + _TransformInvLap_Size;

	//! Allocate buffer space
	char *Buf = State->BufferData = malloc(BUFFER_ALIGNMENT-1 + AllocSize);
	if(!Buf) return -1;

	//! Initialize pointers
	int i, Chan;
	Buf += (-(uintptr_t)Buf) & (BUFFER_ALIGNMENT-1);
	State->TransformBuffer = (float *)(Buf + TransformBuffer_Offs);
	State->TransformTemp   = (float *)(Buf + TransformTemp_Offs);
	State->TransformInvLap = (float**)(Buf + _TransformInvLap_Offs);
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
static inline uint8_t Block_Decode_ReadNybble(const uint8_t **Src, int *Size) {
	//! Fetch and shift nybble
	uint8_t x = *(*Src);
	*Size += 4;
	if((*Size)%8u == 0) x >>= 4, (*Src)++;
	return x; //! NOTE: Unmasked return value
}
static inline float Block_Decode_DecodeQuantizer(const uint8_t **Src, int *Size) {
	int8_t        qi  = Block_Decode_ReadNybble(Src, Size) & 0xF; if(qi == 0xF) return 0.0f;
	if(qi == 0xE) qi += Block_Decode_ReadNybble(Src, Size) & 0xF; //! 8h,0h,Eh,0h..Ch: Extended-precision quantizer
	return 1.0f / ((1u<<5) << qi);
}
int ULC_DecodeBlock(struct ULC_DecoderState_t *State, float *DstData, const uint8_t *SrcBuffer) {
	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	int    nChan            = State->nChan;
	int    BlockSize        = State->BlockSize;
	float *TransformBuffer  = State->TransformBuffer;
	float *TransformTemp    = State->TransformTemp;
	float **TransformInvLap = State->TransformInvLap;

	//! Begin decoding
	int Chan;
	int Size = 0;
	int BlockOverlap = State->NextOverlap;
	State->NextOverlap = BlockSize >> (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF);
	for(Chan=0;Chan<nChan;Chan++) {
		//! Start decoding coefficients
		int32_t v;
		float *CoefDst = TransformBuffer;
		int    CoefRem = BlockSize;
		float  Quant   = Block_Decode_DecodeQuantizer(&SrcBuffer, &Size);
		if(Quant == 0.0f) {
			//! [8h,0h,]Fh: Stop
			do *CoefDst++ = 0.0f; while(--CoefRem);
		} else for(;;) {
			//! -7h..+7h: Normal
			v = ((int32_t)Block_Decode_ReadNybble(&SrcBuffer, &Size) << 28) >> 28;
			if(v != -0x8) {
				//! Store linearized, dequantized coefficient
				v = (v < 0) ? (-v*v) : (+v*v);
				*CoefDst++ = v * Quant;
				if(--CoefRem == 0) break;
				continue;
			}

			//! Unpack escape code
			v = Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF;
			if(v != 0x0) {
				//! 8h,1h..Dh:  4.. 28 zeros (Step: 2)
				//! 8h,Eh,Xh:  30.. 90 zeros (Step: 4)
				//! 8h,Fh,Xh:  94..214 zeros (Step: 8)
				int nZ;
				     if(v == 0xF) nZ = 94 + 8*(Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF);
				else if(v == 0xE) nZ = 30 + 4*(Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF);
				else              nZ =  2 + 2*v;

				//! Clipping to avoid buffer overflow on corrupted blocks
				if(nZ > CoefRem) nZ = CoefRem;

				//! Insert zeros
				CoefRem -= nZ;
				do *CoefDst++ = 0.0f; while(--nZ);
				if(CoefRem == 0) break;
			} else {
				//! 8h,0h,0h..Eh[,Xh]: Quantizer change
				//! 8h,0h,Fh:          Stop
				Quant = Block_Decode_DecodeQuantizer(&SrcBuffer, &Size);
				if(Quant == 0.0f) {
					do *CoefDst++ = 0.0f; while(--CoefRem);
					break;
				}
			}
		}

		//! Inverse transform block
		Fourier_IMDCT(DstData + Chan*BlockSize, TransformBuffer, TransformInvLap[Chan], TransformTemp, BlockSize, BlockOverlap);
	}
	return Size;
}

/**************************************/
//! EOF
/**************************************/
