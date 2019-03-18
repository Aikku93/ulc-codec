/**************************************/
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**************************************/
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

#define MAX_QUANTS 48

/**************************************/

//! Clip to 16bit range
static inline float Clip16(float x) {
	if(x < -0x8000) return -0x8000;
	if(x > +0x7FFF) return +0x7FFF;
	return x;
}

/**************************************/

//! File header
struct FileHeader_t {
	uint32_t Magic[2];  //! [00h] Magic values
	uint32_t nSamp;     //! [08h] Number of samples
	uint32_t RateHz;    //! [0Ch] Playback rate
	uint32_t BlockSize; //! [10h] Transform block size
	uint16_t nChan;     //! [14h] Channels in stream
	uint16_t RateKbps;  //! [16h] Nominal coding rate
};

//! Decoding state
#define MAX_BLOCK_SIZE 4096
#define MAX_CHANS         4
static const size_t CacheMinBlockSize = ((4 + 4*MAX_BLOCK_SIZE + 12*(MAX_QUANTS-1))*MAX_CHANS + 7) / 8;
static const size_t CacheSize = 256*1024;
struct DecodeState_t {
	//! These need to be cleared to NULL
	FILE *FileIn;
	FILE *FileOut;
	char *AllocBuffer;

	float   *BlockBuffer;
	int16_t *BlockOutput;
	uint8_t *CacheBuffer;
	uint8_t *CacheNext;
};

//! Clean up decode state and exit
static void StateCleanupExit(const struct DecodeState_t *State, int ExitCode) {
	free(State->AllocBuffer);
	fclose(State->FileOut);
	fclose(State->FileIn);
	exit(ExitCode);
}

//! Initialize state
static void StateInit(struct DecodeState_t *State, const struct FileHeader_t *Header) {
	//! Allocate memory
	size_t BlockBuffer_Size = sizeof(float)    * Header->nChan*Header->BlockSize;
	size_t BlockOutput_Size = sizeof(int16_t)  * Header->nChan*Header->BlockSize;
	size_t CacheBuffer_Size = CacheSize;
	size_t BlockBuffer_Offs = 0;
	size_t BlockOutput_Offs = BlockBuffer_Offs + BlockBuffer_Size;
	size_t CacheBuffer_Offs = BlockOutput_Offs + BlockOutput_Size;
	size_t AllocSize = CacheBuffer_Offs + CacheBuffer_Size;
	char *Buf = State->AllocBuffer = malloc(BUFFER_ALIGNMENT-1 + AllocSize);
	if(!Buf) {
		printf("ERROR: Out of memory.\n");
		StateCleanupExit(State, -1);
	}

	//! Set pointers
	Buf += -(uintptr_t)Buf % BUFFER_ALIGNMENT;
	State->BlockBuffer = (float   *)(Buf + BlockBuffer_Offs);
	State->BlockOutput = (int16_t *)(Buf + BlockOutput_Offs);
	State->CacheBuffer = (uint8_t *)(Buf + CacheBuffer_Offs);
	State->CacheNext   = State->CacheBuffer;

	//! Fill cache
	fread(State->CacheNext, sizeof(uint8_t), CacheSize, State->FileIn);
}

//! Advance cached data
static void StateCacheAdvance(struct DecodeState_t *State, size_t nBytes) {
	State->CacheNext += nBytes;
	size_t Rem = State->CacheBuffer + CacheSize - State->CacheNext;
	if(Rem < CacheMinBlockSize) {
		memcpy(State->CacheBuffer, State->CacheNext, Rem);
		fread(State->CacheBuffer + Rem, sizeof(uint8_t), CacheSize-Rem, State->FileIn);
		State->CacheNext = State->CacheBuffer;
	}
}

/**************************************/

int main(int argc, const char *argv[]) {
	//! Check arguments
	if(argc != 3) {
		printf(
			"ulcDecodeTool - Ultra-Low Complexity Codec Decoding Tool\n"
			"Usage: ulcdecodetool Input.ulc Output.sw\n"
			"Multi-channel data will be interleaved.\n"
		);
		return 1;
	}

	//! Create decoding state
	struct DecodeState_t State = {
		.FileIn      = NULL,
		.FileOut     = NULL,
		.AllocBuffer = NULL,
	};

	//! Open input file
	State.FileIn = fopen(argv[1], "rb");
	if(!State.FileIn) {
		printf("ERROR: Unable to open input file.\n");
		StateCleanupExit(&State, -1);
	}

	//! Open output file
	State.FileOut = fopen(argv[2], "wb");
	if(!State.FileOut) {
		printf("ERROR: Unable to open output file.\n");
		StateCleanupExit(&State, -1);
	}

	//! Read header
	struct FileHeader_t Header; fread(&Header, sizeof(Header), 1, State.FileIn);
	if(Header.Magic[0] != (uint32_t)('U' | 'L'<<8 | 'C'<<16 | 'b'<<24) || Header.Magic[1] != (uint32_t)(0x01 | 0xFF<<8 | 0x02<<16 | 0xFE<<24)) {
		printf("ERROR: Invalid file.\n");
		StateCleanupExit(&State, -1);
	}
	if(Header.BlockSize > MAX_BLOCK_SIZE || Header.nChan > MAX_CHANS) {
		printf("ERROR: Unsupported specification.\n");
		StateCleanupExit(&State, -1);
	}

	//! Initialize state
	StateInit(&State, &Header);

	//! Create decoder
	struct ULC_DecoderState_t Decoder = {
		.nChan     = Header.nChan,
		.BlockSize = Header.BlockSize,
	};
	if(ULC_DecoderState_Init(&Decoder) > 0) {
		//! Process blocks
		float   *BlockBuffer = State.BlockBuffer;
		int16_t *BlockOutput = State.BlockOutput;
		size_t BlockSize = Header.BlockSize;
		size_t nChan = Header.nChan;
		size_t nBlk = (Header.nSamp + BlockSize-1) / BlockSize;
		for(size_t Blk=0;Blk<nBlk+1;Blk++) { //! +1 to account for coding delay
			//! Show progress
			printf("\rBlock %u/%u (%.2f%%)...", Blk, nBlk, Blk*100.0f/nBlk);
			fflush(stdout);

			//! Decode block
			size_t Size = ULC_DecodeBlock(&Decoder, BlockBuffer, State.CacheNext);
			StateCacheAdvance(&State, (Size + 7) / 8);

			//! Apply M/S transform
			if(nChan == 2) for(size_t n=0;n<BlockSize;n++) {
				float *a = &BlockBuffer[0*BlockSize+n], va = *a;
				float *b = &BlockBuffer[1*BlockSize+n], vb = *b;
				*a = va + vb;
				*b = va - vb;
			}

			//! Interleave to output buffer
			for(size_t Chan=0;Chan<nChan;Chan++) for(size_t n=0;n<BlockSize;n++) {
				BlockOutput[n*nChan+Chan] = (int16_t)Clip16(round(BlockBuffer[Chan*BlockSize+n]));
			}

			//! Write to file
			fwrite(BlockOutput, nChan*sizeof(int16_t), BlockSize, State.FileOut);
		}
	} else printf("ERROR: Unable to initialize decoder.\n");

	//! Done
	printf("\e[2K\rOk\n");
	StateCleanupExit(&State, 0);
	return 0;
}

/**************************************/
//! EOF
/**************************************/
