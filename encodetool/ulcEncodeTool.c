/**************************************/
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
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

//! Coding mode
//! This can be changed without rebuilding the decode tool
static const size_t BlockSize = 2048;
static const size_t nQuants = 8;
static const uint16_t QuantsBw[] = {16,32,48,96,192,384,576,704};

/**************************************/

int main(int argc, const char *argv[]) {
	//! Check arguments
	if(argc < 5 || argc > 6) {
		printf(
			"ulcEncodeTool - Ultra-Low Complexity Codec Encoding Tool\n"
			"Usage: ulcencodetool Input.sw Output.ulc RateHz RateKbps [nChan=1]\n"
			"Multi-channel data must be interleaved.\n"
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

	/*! Write file header !*/ {
		//! Core header
		struct {
			uint32_t Magic[2];  //! [00h] Magic values
			uint32_t nSamp;     //! [08h] Number of samples
			uint32_t RateHz;    //! [0Ch] Playback rate
			uint16_t BlockSize; //! [10h] Transform block size
			uint16_t nQuants;   //! [12h] Quantizer bands
			uint16_t nChan;     //! [14h] Channels in stream
			uint16_t RateKbps;  //! [16h] Nominal coding rate
		} Header = {
			{(uint32_t)('U' | 'L'<<8 | 'C'<<16 | 'a'<<24), (uint32_t)(0x01 | 0xFF<<8 | 0x02<<16 | 0xFE<<24)},
			nSamp,
			RateHz,
			(uint16_t)BlockSize,
			(uint16_t)nQuants,
			(uint16_t)nChan,
			(uint16_t)RateKbps,
		};
		fwrite(&Header, sizeof(Header), 1, OutFile);

		//! Quantizer bands
		fwrite(QuantsBw, sizeof(uint16_t), nQuants, OutFile);
	}

	//! Create encoder
	struct ULC_EncoderState_t Encoder = {
		.RateHz    = RateHz,
		.nChan     = nChan,
		.BlockSize = BlockSize,
		.nQuants   = nQuants,
		.QuantsBw  = QuantsBw,
	};
	if(ULC_EncoderState_Init(&Encoder) > 0) {
		//! Process blocks
		size_t n, Chan;
		size_t Blk, nBlk = (nSamp + BlockSize-1) / BlockSize;
		uint64_t TotalSize = 0;
		for(Blk=0;Blk<nBlk+1;Blk++) { //! +1 to account for coding delay
			//! Show progress
			printf("\rBlock %u/%u (%.2f%%)...", Blk, nBlk, Blk*100.0f/nBlk);
			fflush(stdout);

			//! Fill buffer data
			//! BlockFetch[] is free after this
			size_t nMax = fread(BlockFetch, nChan*sizeof(int16_t), BlockSize, InFile);
			for(Chan=0;Chan<nChan;Chan++) for(n=0;n<BlockSize;n++) {
				BlockBuffer[Chan*BlockSize+n] = (n < nMax) ? BlockFetch[n*nChan+Chan] : 0.0f;
			}

			//! Apply M/S transform
			if(nChan == 2) for(n=0;n<BlockSize;n++) {
				float *a = &BlockBuffer[0*BlockSize+n], va = *a;
				float *b = &BlockBuffer[1*BlockSize+n], vb = *b;
				*a = (va + vb) * 0.5f;
				*b = (va - vb) * 0.5f;
			}

			//! Encode block
			//! Reuse BlockFetch[] to avoid more memory allocation
			uint8_t *EncData = (uint8_t*)BlockFetch;
			size_t Size = ULC_EncodeBlock(&Encoder, EncData, BlockBuffer, RateKbps);
			fwrite(BlockFetch, 1, (Size+7)/8, OutFile);
			TotalSize += Size;
		}

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
