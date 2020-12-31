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
#define HEADER_MAGIC (uint32_t)('U' | 'L'<<8 | 'C'<<16 | 'j'<<24)

/**************************************/

//! Cache memory
//! This avoids too many calls to fwrite(), which
//! hopefully avoids excessive system calls
#define CACHE_SIZE (512 * 1024)
static uint8_t CacheMem[CACHE_SIZE];

/**************************************/

int main(int argc, const char *argv[]) {
	//! Check arguments
	if(argc < 5) {
		printf(
			"ulcEncodeTool - Ultra-Low Complexity Codec Encoding Tool\n"
			"Usage: ulcencodetool Input.sw Output.ulc RateHz RateKbps [Options]\n"
			"Options:\n"
			" -nc:1           - Set number of channels.\n"
			" -blocksize:2048 - Set number of coefficients per block (must be a power of 2).\n"
			"Multi-channel data must be interleaved.\n"
		);
		return 1;
	}

	//! Parse parameters
	int   BlockSize  = 2048;
	int   nChan      = 1;
	int   RateHz     = atoi(argv[3]);
	int   RateKbps   = atoi(argv[4]); {
		int n;
		for(n=5;n<argc;n++) {
			if(!memcmp(argv[n], "-nc:", 4)) {
				int x = atoi(argv[n] + 4);
				if(n > 0 && n < 65535) nChan = x;
				else printf("WARNING: Ignoring invalid parameter to number of channels (%d)\n", x);
			}

			else if(!memcmp(argv[n], "-blocksize:", 11)) {
				int x = atoi(argv[n] + 11);
				if(x >= 64 && x <= 8192 && (x & (-x)) == x) BlockSize = x;
				else printf("WARNING: Ignoring invalid parameter to block size (%d)\n", x);
			}

			else printf("WARNING: Ignoring unknown argument (%s)\n", argv[n]);
		}
	}

	//! Verify parameters
	if(RateHz < 1 || RateHz > 0x7FFFFFFF) {
		printf("ERROR: Invalid playback rate.\n");
		return -1;
	}
	if(RateKbps == 0 || RateKbps > 0xFFFF) {
		printf("ERROR: Invalid coding rate.\n");
		return -1;
	}
	if(nChan < 1 || nChan > 0xFFFF) {
		printf("ERROR: Invalid number of channels.\n");
		return -1;
	}

	//! Allocate buffers
	int16_t *BlockFetch   = malloc(sizeof(int16_t) * nChan*BlockSize);
	char    *_BlockBuffer = malloc(sizeof(float)   * nChan*BlockSize + BUFFER_ALIGNMENT-1);
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
		uint32_t BlockSize;    //! [10h] Transform block size
		uint16_t nChan;        //! [14h] Channels in stream
		uint16_t RateKbps;     //! [16h] Nominal coding rate
	} FileHeader = {
		HEADER_MAGIC,
		0, //! MaxBlockSize filled later
		nSamp,
		RateHz,
		BlockSize,
		nChan,
		RateKbps,
	};
	size_t FileHeaderOffs = ftell(OutFile);
	fseek(OutFile, +sizeof(FileHeader), SEEK_CUR);

	//! Create encoder
	struct ULC_EncoderState_t Encoder = {
		.RateHz     = RateHz,
		.nChan      = nChan,
		.BlockSize  = BlockSize,
	};
	if(ULC_EncoderState_Init(&Encoder) > 0) {
		//! Process blocks
		int n, Chan;
		int CacheIdx = 0;
		size_t Blk, nBlk = (nSamp + BlockSize-1) / BlockSize + 2; //! +1 to account for coding delay, +1 to account for MDCT delay
		uint64_t TotalSize = 0;
		for(Blk=0;Blk<nBlk;Blk++) {
			//! Show progress
			//! NOTE: Only displaying every 4th block; slowdown occurs otherwise
			if(Blk%4u == 0) {
				printf("\rBlock %u/%u (%.2f%%)...", Blk, nBlk, Blk*100.0f/nBlk);
				fflush(stdout);
			}

			//! Fill buffer data
			//! BlockFetch[] is free after this
			size_t nMax = fread(BlockFetch, nChan*sizeof(int16_t), BlockSize, InFile);
			for(Chan=0;Chan<nChan;Chan++) for(n=0;n<BlockSize;n++) {
				BlockBuffer[Chan*BlockSize+n] = ((size_t)n < nMax) ? (BlockFetch[n*nChan+Chan] * (1.0f/32768.0f)) : 0.0f;
			}

			//! Encode block
			//! Reuse BlockBuffer[] to avoid more memory allocation
			uint8_t *EncData = (uint8_t*)BlockBuffer;
			int Size;
			if(RateKbps > 0) Size = ULC_EncodeBlock_CBR(&Encoder, EncData, BlockBuffer,  RateKbps);
			else             Size = ULC_EncodeBlock_VBR(&Encoder, EncData, BlockBuffer, -RateKbps);
			TotalSize += Size;

			//! Copy what we can into the cache
			Size = (Size+7) / 8u;
			if((size_t)Size > FileHeader.MaxBlockSize) FileHeader.MaxBlockSize = Size;
			while(Size) {
				//! Copy up to the limits of the cache area
				int n = CACHE_SIZE - CacheIdx; //! =CacheRem
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
		size_t nEncodedSamples = BlockSize * nBlk;
		printf(
			"\e[2K\r" //! Clear line before CR
			"Total size = %.2fKiB\n"
			"Avg rate = %.5fkbps (%.5f bits/sample)\n"
			"Max rate = %.5fkbps (%.5f bits/sample)\n",
			TotalSize/8.0 / 1024,
			TotalSize               * 1.0 * RateHz/1000.0 / nEncodedSamples,
			TotalSize               * 1.0 / nEncodedSamples,
			FileHeader.MaxBlockSize * 8.0 * RateHz/1000.0 / BlockSize,
			FileHeader.MaxBlockSize * 8.0 / BlockSize
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
