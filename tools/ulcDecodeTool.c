/**************************************/
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/**************************************/
#include "ulcDecoder.h"
/**************************************/
#define BUFFER_ALIGNMENT 64u //! __mm256
/**************************************/

#define HEADER_MAGIC (uint32_t)('U' | 'L'<<8 | 'C'<<16 | '2'<<24)

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
	uint32_t Magic;        //! [00h] Magic value/signature
	uint16_t BlockSize   ; //! [04h] Transform block size
	uint16_t MaxBlockSize; //! [06h] Largest block size (in bytes; 0 = Unknown)
	uint32_t nBlocks;      //! [08h] Number of blocks
	uint32_t RateHz;       //! [0Ch] Playback rate
	uint16_t nChan;        //! [10h] Channels in stream
	uint16_t RateKbps;     //! [12h] Nominal coding rate
	uint32_t StreamOffs;   //! [14h] Offset of data stream
};

//! Decoding state
#define MAX_BLOCK_SIZE 8192
#define MAX_CHANS         4
static const int CacheSize = 512*1024;
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
	int BlockBuffer_Size = sizeof(float)   * Header->nChan*Header->BlockSize;
	int BlockOutput_Size = sizeof(int16_t) * Header->nChan*Header->BlockSize;
	int CacheBuffer_Size = CacheSize;
	int BlockBuffer_Offs = 0;
	int BlockOutput_Offs = BlockBuffer_Offs + BlockBuffer_Size;
	int CacheBuffer_Offs = BlockOutput_Offs + BlockOutput_Size;
	int AllocSize = CacheBuffer_Offs + CacheBuffer_Size;
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

	//! Seek to start of stream and fill cache
	fseek(State->FileIn, Header->StreamOffs, SEEK_SET);
	fread(State->CacheNext, sizeof(uint8_t), CacheSize, State->FileIn);
}

//! Advance cached data
static void StateCacheAdvance(struct DecodeState_t *State, int nBytes, int MinSize) {
	State->CacheNext += nBytes;
	int Rem = State->CacheBuffer + CacheSize - State->CacheNext;
	if(Rem < MinSize) {
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
	if(Header.Magic != HEADER_MAGIC) {
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
		.nChan      = Header.nChan,
		.BlockSize  = Header.BlockSize,
	};
	if(ULC_DecoderState_Init(&Decoder) > 0) {
		const clock_t DISPLAY_UPDATE_RATE = CLOCKS_PER_SEC/2; //! Update every 0.5 seconds

		//! Process blocks
		int      n, Chan;
		int      nChan       = Header.nChan;
		int      BlockSize   = Header.BlockSize;
		float   *BlockBuffer = State.BlockBuffer;
		int16_t *BlockOutput = State.BlockOutput;
		uint32_t Blk, nBlk = Header.nBlocks;
		size_t BlkLastUpdate = 0;
		clock_t LastUpdateTime = clock() - DISPLAY_UPDATE_RATE;
		for(Blk=0;Blk<nBlk;Blk++) {
			//! Show progress
			//! NOTE: Take difference and use unsigned comparison to
			//! get correct results in the comparison on signed overflows.
			//! uint64_t might be overkill, depending on the implementation.
			if((uint64_t)(clock()-LastUpdateTime) >= DISPLAY_UPDATE_RATE) {
				size_t nBlkProcessed = 2 * (Blk-BlkLastUpdate); //! Updated every 0.5s, displayed as X*s^-1
				printf(
					"\rBlock %u/%u (%.2f%% | %.2f X rt)",
					Blk, nBlk, Blk*100.0/nBlk,
					nBlkProcessed*BlockSize / (double)Header.RateHz
				);
				fflush(stdout);
				LastUpdateTime += DISPLAY_UPDATE_RATE;
				BlkLastUpdate   = Blk;
			}

			//! Decode block
			int Size = ULC_DecodeBlock(&Decoder, BlockBuffer, State.CacheNext);
			StateCacheAdvance(&State, (Size + 7) / 8u, Header.MaxBlockSize);

			//! Interleave to output buffer
			for(Chan=0;Chan<nChan;Chan++) for(n=0;n<BlockSize;n++) {
				BlockOutput[n*nChan+Chan] = (int16_t)Clip16(lrintf(32768.0f * BlockBuffer[Chan*BlockSize+n]));
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
