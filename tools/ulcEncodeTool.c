/**************************************/
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/**************************************/
#include "ulcEncoder.h"
/**************************************/
#define BUFFER_ALIGNMENT 64u //! __mm256
/**************************************/

//! Header magic value
#define HEADER_MAGIC (uint32_t)('U' | 'L'<<8 | 'C'<<16 | '2'<<24)

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
			"Usage:\n"
			" ulcencodetool Input Output.ulc RateHz RateKbps[,AvgComplexity]|-Quality [Opt]\n"
			"Options:\n"
			" -nc:1           - Set number of channels.\n"
			" -blocksize:2048 - Set number of coefficients per block (must be a power of 2).\n"
			"Multi-channel data must be interleaved (packed).\n"
			"Passing AvgComplexity uses ABR mode.\n"
			"Passing negative RateKbps (-Quality) uses VBR mode.\n"
		);
		return 1;
	}

	//! Parse parameters
	int BlockSize = 2048;
	int nChan     = 1;
	int RateHz    = atoi(argv[3]);
	float RateKbps, AvgComplexity = 0.0f; sscanf(argv[4], "%f,%f", &RateKbps, &AvgComplexity);
	{
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

	//! Determine encoding mode (CBR/ABR/VBR) and set appropriate block encoding routine
	//! NOTE: This is kinda janky, but works fine. The main problem is passing AvgComplexity,
	//! which requires casting all the function pointers to a single form.
	typedef const uint8_t* (*BlockEncodeFnc_t)(struct ULC_EncoderState_t *State, const float *SrcData, int *Size, float Rate, float AvgComplexity);
	BlockEncodeFnc_t BlockEncodeFnc;
	                          BlockEncodeFnc = (BlockEncodeFnc_t)ULC_EncodeBlock_CBR;
	if(AvgComplexity > 0.0f)  BlockEncodeFnc = (BlockEncodeFnc_t)ULC_EncodeBlock_ABR;
	if(RateKbps < 0.0f)       BlockEncodeFnc = (BlockEncodeFnc_t)ULC_EncodeBlock_VBR, RateKbps = -RateKbps;

	//! Verify parameters
	if(RateHz < 1 || RateHz > 0x7FFFFFFF) {
		printf("ERROR: Invalid playback rate.\n");
		return -1;
	}
	if(RateKbps == 0.0f) {
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
		uint16_t BlockSize;    //! [04h] Transform block size
		uint16_t MaxBlockSize; //! [06h] Largest block size (in bytes; 0 = Unknown)
		uint32_t nBlocks;      //! [08h] Number of blocks
		uint32_t RateHz;       //! [0Ch] Playback rate
		uint16_t nChan;        //! [10h] Channels in stream
		uint16_t RateKbps;     //! [12h] Nominal coding rate
		uint32_t StreamOffs;   //! [14h] Offset of data stream
	} FileHeader = {
		.Magic      = HEADER_MAGIC,
		.BlockSize  = BlockSize,
		.nBlocks    = (nSamp + BlockSize-1) / BlockSize + 2, //! +1 to account for coding delay, +1 to account for MDCT delay
		.RateHz     = RateHz,
		.nChan      = nChan,
		.RateKbps   = (uint16_t)RateKbps,
	};
	size_t FileHeaderOffs = ftell(OutFile);
	fseek(OutFile, +sizeof(FileHeader), SEEK_CUR);

	//! Create encoder
	struct ULC_EncoderState_t Encoder = {
		.RateHz     = RateHz,
		.nChan      = nChan,
		.BlockSize  = BlockSize,
		.ModulationWindow = NULL,
	};
	if(ULC_EncoderState_Init(&Encoder) > 0) {
		const clock_t DISPLAY_UPDATE_RATE = CLOCKS_PER_SEC/2; //! Update every 0.5 seconds

		//! Store stream offset
		FileHeader.StreamOffs = ftell(OutFile);

		//! Process blocks
		int n, Chan;
		int CacheIdx = 0;
		size_t Blk, nBlk = FileHeader.nBlocks;
		uint64_t TotalSize = 0;
		double ActualAvgComplexity = 0.0;
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
					"\rBlock %u/%u (%.2f%% | %.2f X rt) | Average: %.2fkbps",
					Blk, nBlk, Blk*100.0/nBlk,
					nBlkProcessed*BlockSize / (double)RateHz,
					Blk ? (TotalSize * RateHz/1000.0 / (Blk * BlockSize)) : 0.0f
				);
				fflush(stdout);
				LastUpdateTime += DISPLAY_UPDATE_RATE;
				BlkLastUpdate   = Blk;
			}

			//! Fill buffer data
			size_t nMax = fread(BlockFetch, nChan*sizeof(int16_t), BlockSize, InFile);
			for(Chan=0;Chan<nChan;Chan++) for(n=0;n<BlockSize;n++) {
				BlockBuffer[Chan*BlockSize+n] = ((size_t)n < nMax) ? (BlockFetch[n*nChan+Chan] * (1.0f/32768.0f)) : 0.0f;
			}

			//! Encode block
			//! Reuse BlockBuffer[] to avoid more memory allocation
			int Size;
			const uint8_t *EncData = BlockEncodeFnc(&Encoder, BlockBuffer, &Size, RateKbps, AvgComplexity);
			TotalSize += Size;
			ActualAvgComplexity += Encoder.BlockComplexity;

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
			"Max rate = %.5fkbps (%.5f bits/sample)\n"
			"Avg complexity = %.5f\n",
			TotalSize/8.0 / 1024,
			TotalSize               * 1.0 * RateHz/1000.0 / nEncodedSamples,
			TotalSize               * 1.0 / nEncodedSamples,
			FileHeader.MaxBlockSize * 8.0 * RateHz/1000.0 / BlockSize,
			FileHeader.MaxBlockSize * 8.0 / BlockSize,
			ActualAvgComplexity/nBlk
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
