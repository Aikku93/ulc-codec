/**************************************/
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/**************************************/
#include "ulcEncoder.h"
/**************************************/
#if defined(__AVX__)
# define BUFFER_ALIGNMENT 32u //! __mm256
#elif defined(__SSE__)
# define BUFFER_ALIGNMENT 16u //! __mm128
#else
# define BUFFER_ALIGNMENT 4u //! float
#endif
/**************************************/

//! Header magic value
#define HEADER_MAGIC (uint32_t)('U' | 'L'<<8 | 'C'<<16 | 'c'<<24)

//! Transform block size
//! Feel free to change this; the decoder doesn't care
#define BLOCK_SIZE    2048
#define BLOCK_OVERLAP 1536

/**************************************/

//! Cache memory
//! This avoids too many calls to fwrite(), which
//! hopefully avoids excessive system calls
#define CACHE_SIZE (256 * 1024) //! 256KiB
static uint8_t CacheMem[CACHE_SIZE];

/**************************************/

int main(int argc, const char *argv[]) {
	//! Check arguments
	int MidSideXfm = 1;
	if(argc > 6) MidSideXfm = (!strcmp(argv[6], "-nomidside")) ? 0 : (-1);
	if(argc < 5 || argc > 7) {
		printf(
			"ulcEncodeTool - Ultra-Low Complexity Codec Encoding Tool\n"
			"Usage: ulcencodetool Input.sw Output.ulc RateHz RateKbps [nChan=1 [-nomidside]]\n"
			"Multi-channel data must be interleaved.\n"
			"-nomidside disables M/S stereo.\n"
		);
		return 1;
	}

	//! Parse parameters
	size_t RateHz   = atoi(argv[3]);
	size_t RateKbps = atoi(argv[4]);
	size_t nChan    = (argc > 5) ? atoi(argv[5]) : 1;

	//! Verify parameters
	if(RateHz < 1 || RateHz > 0xFFFFFFFFu) {
		printf("ERROR: Invalid playback rate.\n");
		return -1;
	}
	if(RateKbps < 1 || RateKbps > 0xFFFFu) {
		printf("ERROR: Invalid coding rate.\n");
		return -1;
	}
	if(nChan < 1 || nChan > 0xFFFFu) {
		printf("ERROR: Invalid number of channels.\n");
		return -1;
	}

	//! Allocate buffers
	int16_t *BlockFetch   = malloc(sizeof(int16_t) * nChan*BLOCK_SIZE);
	char    *_BlockBuffer = malloc(sizeof(float)   * nChan*BLOCK_SIZE + BUFFER_ALIGNMENT-1);
	if(!BlockFetch || !_BlockBuffer) {
		printf("ERROR: Out of memory.\n");
		free(_BlockBuffer);
		free(BlockFetch);
		return -1;
	}
	float *BlockBuffer = (float*)(_BlockBuffer + (-(uintptr_t)_BlockBuffer % BUFFER_ALIGNMENT));

	//! Open input file
	size_t nSamp;
	FILE *InFile = fopen(argv[1], "rb");
	if(!InFile) {
		printf("ERROR: Unable to open input file.\n");
		free(_BlockBuffer);
		free(BlockFetch);
		return -1;
	} else {
		fseek(InFile, 0, SEEK_END);
		nSamp = ftell(InFile) / sizeof(int16_t) / nChan;
		rewind(InFile);
	}

	//! Open output file
	FILE *OutFile = fopen(argv[2], "wb");
	if(!OutFile) {
		printf("ERROR: Unable to open output file.\n");
		fclose(InFile);
		free(_BlockBuffer);
		free(BlockFetch);
		return -1;
	}

	//! Create file header and skip for now; written later
	struct {
		uint32_t Magic;        //! [00h] Magic value/signature
		uint32_t MaxBlockSize; //! [04h] Largest block size (in bytes; 0 = Unknown)
		uint32_t nSamp;        //! [08h] Number of samples
		uint32_t RateHz;       //! [0Ch] Playback rate
		uint16_t BlockSize;    //! [10h] Transform block size
		uint16_t BlockOverlap; //! [12h] Block overlap
		uint16_t nChan;        //! [14h] Channels in stream
		uint16_t RateKbps;     //! [16h] Nominal coding rate
	} FileHeader = {
		HEADER_MAGIC,
		0, //! MaxBlockSize filled later
		nSamp,
		RateHz,
		BLOCK_SIZE,
		BLOCK_OVERLAP,
		(uint16_t)nChan,
		(uint16_t)RateKbps,
	};
	long int FileHeaderOffs = ftell(OutFile);
	fseek(OutFile, +sizeof(FileHeader), SEEK_CUR);

	//! Create encoder
	struct ULC_EncoderState_t Encoder = {
		.RateHz       = RateHz,
		.nChan        = nChan,
		.BlockSize    = BLOCK_SIZE,
		.BlockOverlap = BLOCK_OVERLAP,
	};
	if(ULC_EncoderState_Init(&Encoder) > 0) {
		//! Process blocks
		size_t n, Chan;
		size_t Blk, nBlk = (nSamp + BLOCK_SIZE-1) / BLOCK_SIZE;
		uint64_t TotalSize = 0;
		size_t CacheIdx = 0;
		for(Blk=0;Blk<nBlk+1;Blk++) { //! +1 to account for coding delay
			//! Show progress
			printf("\rBlock %u/%u (%.2f%%)...", Blk, nBlk, Blk*100.0f/nBlk);
			fflush(stdout);

			//! Fill buffer data
			//! BlockFetch[] is free after this
			size_t nMax = fread(BlockFetch, nChan*sizeof(int16_t), BLOCK_SIZE, InFile);
			for(Chan=0;Chan<nChan;Chan++) for(n=0;n<BLOCK_SIZE;n++) {
				BlockBuffer[Chan*BLOCK_SIZE+n] = (n < nMax) ? BlockFetch[n*nChan+Chan] : 0.0f;
			}

			//! Apply M/S transform
			if(MidSideXfm && nChan == 2) for(n=0;n<BLOCK_SIZE;n++) {
				float *a = &BlockBuffer[0*BLOCK_SIZE+n], va = *a;
				float *b = &BlockBuffer[1*BLOCK_SIZE+n], vb = *b;
				*a = (va + vb) * 0.5f;
				*b = (va - vb) * 0.5f;
			}

			//! Encode block
			//! Reuse BlockFetch[] to avoid more memory allocation
			uint8_t *EncData = (uint8_t*)BlockFetch;
			size_t Size = ULC_EncodeBlock(&Encoder, EncData, BlockBuffer, RateKbps, MidSideXfm ? 0x1.6A09E6p-1f : 1.0f); //! 1/sqrt[2]);
			TotalSize += Size;

			//! Copy what we can into the cache
			Size = (Size+7) / 8;
			if(Size > FileHeader.MaxBlockSize) FileHeader.MaxBlockSize = Size;
			while(Size) {
				//! Copy up to the limits of the cache area
				size_t n = CACHE_SIZE - CacheIdx; //! =CacheRem
				if(Size < n) n = Size;
				Size -= n;

				memcpy(CacheMem+CacheIdx, EncData, n);
				EncData += n;
				CacheIdx += n;
				if(CacheIdx == CACHE_SIZE) {
					//! Flush to file
					fwrite(CacheMem, sizeof(uint8_t), CacheIdx, OutFile);
					CacheIdx = 0;
				}
			}
		}

		//! Flush cache
		fwrite(CacheMem, sizeof(uint8_t), CacheIdx, OutFile);

		//! Show statistics
		printf(
			"\e[2K\r"
			"Total size = %.2fKiB\n"
			"Avg rate = %.5fkbps (%.5f bits/sample)\n",
			TotalSize/8.0 / 1024,
			TotalSize*1.0 * RateHz/1000.0 / nSamp,
			TotalSize*1.0 / nSamp
		);

		//! Destroy encoder
		ULC_EncoderState_Destroy(&Encoder);
	} else printf("ERROR: Unable to initialize encoder.\n");

	//! Write file header
	fseek(OutFile, FileHeaderOffs, SEEK_SET);
	fwrite(&FileHeader, sizeof(FileHeader), 1, OutFile);

	//! Clean up
	fclose(OutFile);
	fclose(InFile);
	free(_BlockBuffer);
	free(BlockFetch);
	return 0;
}

/**************************************/
//! EOF
/**************************************/
