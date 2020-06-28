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
static void DecodeBlock_BufferDeinterleave(float *Buf, float *Tmp, int BlockSize, int nSubBlocks) {
	int n, SubBlock;
	int SubBlockSize = BlockSize / nSubBlocks;
	for(n=0;n<BlockSize;n++) Tmp[n] = Buf[n];
	for(SubBlock=0;SubBlock<nSubBlocks;SubBlock++) {
		for(n=0;n<SubBlockSize;n++) Buf[SubBlock*SubBlockSize+n] = Tmp[n*nSubBlocks+SubBlock];
	}
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
	int OverlapSize  = BlockSize;
	int SubBlockSize = BlockSize;
	int nSubBlocks   = 1; {
		int OverlapScale = ((int32_t)Block_Decode_ReadNybble(&SrcBuffer, &Size) << 28) >> 28;
		if(OverlapScale < 0) {
			SubBlockSize >>= -OverlapScale;
			nSubBlocks   <<= -OverlapScale;
			OverlapSize    =  SubBlockSize;
		} else OverlapSize >>= OverlapScale;
	}
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
				//! Short run?
				int nZ;
				if(v != 0xF) nZ = v + 2; //! 8h,1h..Eh: 3 .. 16 zeros
				else {
					//! 8h,Fh,Yh,Xh: 17 .. 272 zeros
					//! Read 16-zeros chunks followed by the tail
					nZ  = (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF);
					nZ  = (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF) | (nZ<<4);
					nZ += 17;
				}

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

		//! Long block?
		      float *Dst = DstData + Chan*BlockSize;
		const float *Src = TransformBuffer;
		      float *Lap = TransformInvLap[Chan];
		if(nSubBlocks == 1) {
			Fourier_IMDCT(Dst, Src, Lap, TransformTemp, SubBlockSize, OverlapSize);
		} else {
			//! Deinterleave coefficients
			DecodeBlock_BufferDeinterleave(TransformBuffer, TransformTemp, BlockSize, nSubBlocks);

			//! Output the non-overlap data we had from the previous block
			int n;
			for(n=0;n<(BlockSize-SubBlockSize)/2;n++) *Dst++ = Lap[BlockSize/2-1-n];

			//! The lapping buffer is now primed for short transform blocks, so decode
			//! half of the subblocks out to the stream, and the other half back to the
			//! lapping buffer
			float *DstLap = Lap + BlockSize/2;
			float *LapTemp = TransformTemp + SubBlockSize;
			for(n=0;n<nSubBlocks/2;n++) {
				Fourier_IMDCT(Dst, Src, Lap, TransformTemp, SubBlockSize, OverlapSize);
				Dst += SubBlockSize;
				Src += SubBlockSize;
			}
			{
				//! We now append the first half of this middle block to the
				//! end of the output block and begin writing back to the
				//! lapping buffer
				Fourier_IMDCT(LapTemp, Src, Lap, TransformTemp, SubBlockSize, OverlapSize);
				Src    += SubBlockSize;
				int i;
				for(i=0;i<SubBlockSize/2;i++) {
					*Dst++ = LapTemp[i];
					*--DstLap = LapTemp[i+SubBlockSize/2];
				}
			}
			for(n++;n<nSubBlocks;n++) {
				Fourier_IMDCT(LapTemp, Src, Lap, TransformTemp, SubBlockSize, OverlapSize);
				Src    += SubBlockSize;
				int i;
				for(i=0;i<SubBlockSize;i++) {
					*--DstLap = LapTemp[i];
				}
			}
		}
	}
	return Size;
}

/**************************************/
//! EOF
/**************************************/
