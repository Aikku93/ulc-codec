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
#define ESCAPE_SEQUENCE_STOP           (-1)
#define ESCAPE_SEQUENCE_STOP_NOISEFILL (-2)
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
	return x&0xF;
}
static inline int Block_Decode_ReadQuantizer(const uint8_t **Src, int *Size) {
	int           qi  = Block_Decode_ReadNybble(Src, Size); //! 8h,0h,0h..Dh:      Quantizer change
	if(qi == 0xF) return ESCAPE_SEQUENCE_STOP_NOISEFILL;    //! 8h,0h,Fh,Zh,Yh,Xh: Noise fill (to end; exp-decay)
	if(qi == 0xE) qi += Block_Decode_ReadNybble(Src, Size); //! 8h,0h,Eh,0h..Ch:   Quantizer change (extended precision)
	if(qi == 0xE + 0xF) return ESCAPE_SEQUENCE_STOP;        //! 8h,0h,Eh,Fh:       Zeros fill (to end)
	return qi;
}
static inline float Block_Decode_ExpandQuantizer(int qi) {
	return 1.0f / ((1u<<5) << qi);
}
static inline void Block_Decode_DecodeSubBlockCoefs(float *CoefDst, int N, const uint8_t **Src, int *Size) {
	int32_t n, v;

	//! Check first quantizer for Stop code
	v = Block_Decode_ReadQuantizer(Src, Size);
	if(v == ESCAPE_SEQUENCE_STOP) {
		//! [8h,0h,]Eh,Fh: Stop
		do *CoefDst++ = 0.0f; while(--N);
		return;
	}

	//! Unpack the [sub]block's coefficients
	float Quant = Block_Decode_ExpandQuantizer(v);
	for(;;) {
		//! -7h..-1h, +1..+7h: Normal
		v = Block_Decode_ReadNybble(Src, Size);
		if(v != 0x8 && v != 0x0) {
			//! Store linearized, dequantized coefficient
			v = (v^0x8) - 0x8; //! Sign extension
			v = (v < 0) ? (-v*v) : (+v*v);
			*CoefDst++ = v * Quant;
			if(--N == 0) break;
			continue;
		}

		//! 0h,Zh,Yh,Xh: 16 .. 527 noise fill (Xh.bit[1..3] != 0)
		//! 0h,Zh,Yh,Xh: 31 .. 542 zeros fill (Xh.bit[1..3] == 0)
		if(v == 0x0) {
			n  = Block_Decode_ReadNybble(Src, Size);
			n  = Block_Decode_ReadNybble(Src, Size) | (n<<4);
			v  = Block_Decode_ReadNybble(Src, Size);
			n  = (v&1) + (n<<1), v >>= 1;
			n += (v == 0) ? 31 : 16;
			if(n > N) n = N; //! <- Clip on corrupt blocks
			N -= n;
			if(v) {
				float p = v * Quant;
				do *CoefDst++ = p * Block_Decode_RandomCoef(); while(--n);
			} else do *CoefDst++ = 0.0f; while(--n);
			if(N == 0) break;
			continue;
		}

		//! 8h,1h..Fh: Zeros fill (1 .. 15 coefficients)
		v = Block_Decode_ReadNybble(Src, Size);
		if(v != 0x0) {
			n = v;
			if(n > N) n = N; //! <- Clip on corrupt blocks
			N -= n;
			do *CoefDst++ = 0.0f; while(--n);
			if(N == 0) break;
			continue;
		}

		//! 8h,0h,0h..Dh:    Quantizer change
		//! 8h,0h,Eh,0h..Ch: Quantizer change (extended precision)
		v = Block_Decode_ReadQuantizer(Src, Size);
		if(v >= 0) {
			Quant = Block_Decode_ExpandQuantizer(v);
			continue;
		}

		//! 8h,0h,Fh,Zh,Yh[,Xh]: Noise fill (to end; exp-decay)
		if(v == ESCAPE_SEQUENCE_STOP_NOISEFILL) {
			v = Block_Decode_ReadNybble(Src, Size);
			n = Block_Decode_ReadNybble(Src, Size);
			if(v&1) n = Block_Decode_ReadNybble(Src, Size) | (n<<4);
			v = 1 + (v>>1);
			n = n;
			float p = v * Quant;
			float r = 1.0f + (n*n)*-0x1.0p-16f;
			do {
				*CoefDst++ = p * Block_Decode_RandomCoef();
				p *= r;
			} while(--N);
			break;
		}

		//! 8h,0h,Eh,Dh: Unused
		//! 8h,0h,Eh,Eh: Unused
		//! 8h,0h,Eh,Fh: Zeros fill (to end)
		if(v == ESCAPE_SEQUENCE_STOP) {
			do *CoefDst++ = 0.0f; while(--N);
			break;
		}
	}
}
int ULC_DecodeBlock(struct ULC_DecoderState_t *State, float *DstData, const void *_SrcBuffer) {
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
	const uint8_t *SrcBuffer = _SrcBuffer;

	//! Begin decoding
	int Chan, Size = 0;
	int LastOverlapSize = 0; //! <- Shuts gcc up
	int WindowCtrl; {
		//! Read window control information
		WindowCtrl = Block_Decode_ReadNybble(&SrcBuffer, &Size);
		if(WindowCtrl & 0x8) WindowCtrl |= Block_Decode_ReadNybble(&SrcBuffer, &Size) << 4;
		else                 WindowCtrl |= 1 << 4;
	}
	for(Chan=0;Chan<nChan;Chan++) {
		//! Reset overlap scaling for this channel
		LastOverlapSize = State->OverlapSize;

		//! Process subblocks
		float *Dst = DstData + Chan*BlockSize;
		float *Src = TransformBuffer;
		float *Lap = TransformInvLap;
		ULC_SubBlockDecimationPattern_t DecimationPattern = ULC_SubBlockDecimationPattern(WindowCtrl);
		do {
			int SubBlockSize = BlockSize >> (DecimationPattern&0x7);
			Block_Decode_DecodeSubBlockCoefs(Src, SubBlockSize, &SrcBuffer, &Size);

			//! Get+update overlap size and limit to that of the last subblock
			int OverlapSize = SubBlockSize;
			if(DecimationPattern&0x8)
				OverlapSize >>= (WindowCtrl & 0x7);
			if(OverlapSize > LastOverlapSize)
				OverlapSize = LastOverlapSize;
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
		} while(DecimationPattern >>= 4);

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
