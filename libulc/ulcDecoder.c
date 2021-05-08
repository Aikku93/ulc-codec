/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
/**************************************/
#include "Fourier.h"
#include "ulcDecoder.h"
#include "ulcHelper.h"
/**************************************/
#define BUFFER_ALIGNMENT 64u //! Always align memory to 64-byte boundaries (preparation for AVX-512)
/**************************************/

//! Just for consistency
#define MIN_CHANS    1
#define MAX_CHANS  255
#define MIN_BANDS  256
#define MAX_BANDS 8192

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

	//! Get buffer offsets and allocation size
	int AllocSize = 0;
#define CREATE_BUFFER(Name, Sz) int Name##_Offs = AllocSize; AllocSize += Sz
	CREATE_BUFFER(TransformBuffer, sizeof(float) * (nChan* BlockSize   ));
	CREATE_BUFFER(TransformTemp,   sizeof(float) * (       BlockSize   ));
	CREATE_BUFFER(TransformInvLap, sizeof(float) * (nChan*(BlockSize/2)));
#undef CREATE_BUFFER

	//! Allocate buffer space
	char *Buf = State->BufferData = malloc(BUFFER_ALIGNMENT-1 + AllocSize);
	if(!Buf) return -1;

	//! Initialize state
	int i;
	Buf += (-(uintptr_t)Buf) & (BUFFER_ALIGNMENT-1);
	State->OverlapSize     = 0;
	State->TransformBuffer = (float*)(Buf + TransformBuffer_Offs);
	State->TransformTemp   = (float*)(Buf + TransformTemp_Offs);
	State->TransformInvLap = (float*)(Buf + TransformInvLap_Offs);
	for(i=0;i<nChan*(BlockSize/2);i++) State->TransformInvLap[i] = 0.0f;

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
static inline float Block_Decode_RandomCoef(void) {
	static uint32_t Seed = 1234567;
	Seed ^= Seed << 13; //! Xorshift
	Seed ^= Seed >> 17;
	Seed ^= Seed <<  5;
	return (int32_t)Seed * 0x1.0p-31f;
}
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
static void DecodeBlock_BufferDeinterleave(float *Dst, const float *Src, int BlockSize, int Decimation) {
	int n;
	switch(Decimation >> 1) { //! Lowermost bit only controls which subblock gets overlap scaling, so ignore it
		//! 001x: a=N/2, b=N/2
		case 0b001: {
			float *DstA = Dst;
			float *DstB = DstA + BlockSize/2;
			for(n=0;n<BlockSize/2;n++) {
				*DstA++ = *Src++;
				*DstB++ = *Src++;
			}
		} break;

		//! 010x: a=N/4, b=N/4, c=N/2
		case 0b010: {
			float *DstA = Dst;
			float *DstB = DstA + BlockSize/4;
			float *DstC = DstB + BlockSize/4;
			for(n=0;n<BlockSize/4;n++) {
				*DstA++ = *Src++;
				*DstB++ = *Src++;
				*DstC++ = *Src++;
				*DstC++ = *Src++;
			}
		} break;

		//! 011x: a=N/2, b=N/4, c=N/4
		case 0b011: {
			float *DstA = Dst;
			float *DstB = DstA + BlockSize/2;
			float *DstC = DstB + BlockSize/4;
			for(n=0;n<BlockSize/4;n++) {
				*DstA++ = *Src++;
				*DstA++ = *Src++;
				*DstB++ = *Src++;
				*DstC++ = *Src++;
			}
		} break;

		//! 100x: a=N/8, b=N/8, c=N/4, d=N/2
		case 0b100: {
			float *DstA = Dst;
			float *DstB = DstA + BlockSize/8;
			float *DstC = DstB + BlockSize/8;
			float *DstD = DstC + BlockSize/4;
			for(n=0;n<BlockSize/8;n++) {
				*DstA++ = *Src++;
				*DstB++ = *Src++;
				*DstC++ = *Src++;
				*DstC++ = *Src++;
				*DstD++ = *Src++;
				*DstD++ = *Src++;
				*DstD++ = *Src++;
				*DstD++ = *Src++;
			}
		} break;

		//! 101x: a=N/4, b=N/8, c=N/8, d=N/2
		case 0b101: {
			float *DstA = Dst;
			float *DstB = DstA + BlockSize/4;
			float *DstC = DstB + BlockSize/8;
			float *DstD = DstC + BlockSize/8;
			for(n=0;n<BlockSize/8;n++) {
				*DstA++ = *Src++;
				*DstA++ = *Src++;
				*DstB++ = *Src++;
				*DstC++ = *Src++;
				*DstD++ = *Src++;
				*DstD++ = *Src++;
				*DstD++ = *Src++;
				*DstD++ = *Src++;
			}
		} break;

		//! 110x: a=N/2, b=N/8, c=N/8, d=N/4
		case 0b110: {
			float *DstA = Dst;
			float *DstB = DstA + BlockSize/2;
			float *DstC = DstB + BlockSize/8;
			float *DstD = DstC + BlockSize/8;
			for(n=0;n<BlockSize/8;n++) {
				*DstA++ = *Src++;
				*DstA++ = *Src++;
				*DstA++ = *Src++;
				*DstA++ = *Src++;
				*DstB++ = *Src++;
				*DstC++ = *Src++;
				*DstD++ = *Src++;
				*DstD++ = *Src++;
			}
		} break;

		//! 111x: a=N/2, b=N/4, c=N/8, d=N/8
		case 0b111: {
			float *DstA = Dst;
			float *DstB = DstA + BlockSize/2;
			float *DstC = DstB + BlockSize/4;
			float *DstD = DstC + BlockSize/8;
			for(n=0;n<BlockSize/8;n++) {
				*DstA++ = *Src++;
				*DstA++ = *Src++;
				*DstA++ = *Src++;
				*DstA++ = *Src++;
				*DstB++ = *Src++;
				*DstB++ = *Src++;
				*DstC++ = *Src++;
				*DstD++ = *Src++;
			}
		} break;
	}
}
int ULC_DecodeBlock(struct ULC_DecoderState_t *State, float *DstData, const uint8_t *SrcBuffer) {
	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	int    n;
	int    nChan           = State->nChan;
	int    BlockSize       = State->BlockSize;
	float *TransformBuffer = State->TransformBuffer;
	float *TransformTemp   = State->TransformTemp;
	float *TransformInvLap = State->TransformInvLap;
	const float *ModulationWindow = State->ModulationWindow;

	//! Begin decoding
	int Chan, Size = 0;
	int LastOverlapSize = 0; //! <- Shuts gcc up
	int WindowCtrl; {
		//! Read window control information
		WindowCtrl = Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF;
		if(WindowCtrl & 0x8) WindowCtrl |= (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF) << 4;
		else                 WindowCtrl |= 1 << 4;
	}
	for(Chan=0;Chan<nChan;Chan++) {
		//! Reset overlap scaling for this channel
		LastOverlapSize = State->OverlapSize;

		//! Start decoding coefficients
		//! NOTE: If the block is decimated, then output coefficients to the
		//! temp buffer, so that they interleave to the transform buffer.
		int32_t v;
		float *CoefDst = (WindowCtrl & 0x8) ? TransformTemp : TransformBuffer;
		int    CoefRem = BlockSize;
		float  Quant   = Block_Decode_DecodeQuantizer(&SrcBuffer, &Size);
		if(Quant == 0.0f) {
			//! [8h,0h,]Fh: Stop
			do *CoefDst++ = 0.0f; while(--CoefRem);
		} else for(;;) {
			//! -7h..-1h, +1..+7h: Normal
			v = ((int32_t)Block_Decode_ReadNybble(&SrcBuffer, &Size) << 28) >> 28;
			if(v != -0x8 && v != 0x0) {
				//! Store linearized, dequantized coefficient
				v = (v < 0) ? (-v*v) : (+v*v);
				*CoefDst++ = v * Quant;
				if(--CoefRem == 0) break;
				continue;
			}

			//! 0h,Zh,Yh,Xh: Noise fill (16 .. 271 coefficients)
			if(v == 0x0) {
				n  = (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF);
				n  = (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF) | (n<<4);
				n += 16;
				if(n > CoefRem) n = CoefRem; //! <- Clip on corrupt blocks
				CoefRem -= n;

				//! NOTE: The scale is quantized in higher precision. See
				//! ulcEncoder_Encode.h for details.
				v = (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF) + 1;
				float p = (v*v) * Quant * (1 / 8.0f);
				do *CoefDst++ = p * Block_Decode_RandomCoef(); while(--n);
				if(CoefRem == 0) break;
				continue;
			}

			//! 8h,1h..Eh:   Zero run ( 1 ..  14 coefficients)
			//! 8h,Fh,Yh,Xh: Zero run (29 .. 284 coefficients)
			v = (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF);
			if(v != 0x0) {
				if(v < 0xF) n = v;
				else {
					n  = (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF);
					n  = (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF) | (n<<4);
					n += 29;
				}
				if(n > CoefRem) n = CoefRem; //! <- Clip on corrupt blocks
				CoefRem -= n;
				do *CoefDst++ = 0.0f; while(--n);
				if(CoefRem == 0) break;
				continue;
			}

			//! 8h,0h,0h..Eh[,0h..Ch]: Quantizer change
			//! 8h,0h,Fh,0h:           Stop
			//! 8h,0h,Fh,1h..Fh,Xh:    Noise fill (exp-decay to end)
			float q = Block_Decode_DecodeQuantizer(&SrcBuffer, &Size);
			if(q != 0.0f) Quant = q;
			else {
				v = Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF;
				if(v) {
					float p = (v*v) * Quant * (1 / 8.0f);
					v = (Block_Decode_ReadNybble(&SrcBuffer, &Size) & 0xF) + 1;
					float Decay = 1.0f - (v*v)*0x1.0p-10f; //! (1/32)^2
					do {
						*CoefDst++ = p * Block_Decode_RandomCoef();
						p *= Decay;
					} while(--CoefRem);
				} else do *CoefDst++ = 0.0f; while(--CoefRem);
				break;
			}
		}

		//! Deinterleave coefficients as needed
		if(WindowCtrl & 0x8) DecodeBlock_BufferDeinterleave(TransformBuffer, TransformTemp, BlockSize, WindowCtrl >> 4);

		//! Process subblocks
		int SubBlockIdx, TransientSubBlockIdx = ULC_Helper_TransientSubBlockIndex(WindowCtrl >> 4);
		const int8_t *DecimationPattern = ULC_Helper_SubBlockDecimationPattern(WindowCtrl >> 4);
		      float *Dst = DstData + Chan*BlockSize;
		const float *Src = TransformBuffer, *SrcEnd = TransformBuffer + BlockSize;
		      float *Lap = TransformInvLap;
		for(SubBlockIdx=0;(Src < SrcEnd);SubBlockIdx++) {
			//! Get the subblock size and overlap
			int SubBlockSize = BlockSize >> DecimationPattern[SubBlockIdx];
			int OverlapSize  = SubBlockSize;
			if(SubBlockIdx == TransientSubBlockIdx)
				OverlapSize >>= (WindowCtrl & 0x7);

			//! Limit overlap to that of the last subblock
			if(OverlapSize > LastOverlapSize) OverlapSize = LastOverlapSize;
			LastOverlapSize = SubBlockSize;

			//! A single long block can be read straight into the output buffer
			if(SubBlockSize == BlockSize) {
				Fourier_IMDCT(Dst, Src, Lap, TransformTemp, SubBlockSize, OverlapSize, ModulationWindow);
				break;
			}

			//! For small blocks, we store the decoded data to a scratch buffer
			float *DecBuf = TransformTemp + SubBlockSize;
			Fourier_IMDCT(DecBuf, Src, Lap, TransformTemp, SubBlockSize, OverlapSize, ModulationWindow);
			Src += SubBlockSize;

			//! Output samples from the lapping buffer, and cycle
			//! the new samples through it for the next call
			float *LapBuf = Lap + BlockSize/2;
			int LapExtra = (BlockSize - SubBlockSize) / 2;
			int LapCopy = (LapExtra < SubBlockSize) ? LapExtra : SubBlockSize;
			for(n=0;n<LapCopy;n++)   *Dst++ = LapBuf[-1-n];
			for(;n<SubBlockSize;n++) *Dst++ = *DecBuf++;
			int NewCopy = LapExtra - LapCopy;
			for(n=0;n<NewCopy;n++) LapBuf[-1-n] = LapBuf[-1-n-LapCopy];
			for(;n<LapExtra;n++) LapBuf[-1-n] = *DecBuf++;
		}

		//! Move to next channel
		TransformInvLap += BlockSize/2;
	}

	//! Undo M/S transform
	//! NOTE: Not orthogonal; must be fully normalized on the encoder side.
	if(nChan == 2) for(n=0;n<BlockSize;n++) {
		float M = DstData[n];
		float S = DstData[n + BlockSize];
		DstData[n]             = M+S;
		DstData[n + BlockSize] = M-S;
	}

	//! Store the last overlap, and return the number of bits read
	State->OverlapSize = LastOverlapSize;
	return Size;
}

/**************************************/
//! EOF
/**************************************/
