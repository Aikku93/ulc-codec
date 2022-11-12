/**************************************/
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/**************************************/
#include "ulc_Helper.h"
#include "ulcEncoder.h"
#include "WavIO.h"
/**************************************/

int main(int argc, const char *argv[]) {
	int   ExitCode = 0;
	FILE *FileOut;
	char *AllocBuffer;
	struct WAV_State_t FileIn;
	struct ULC_EncoderState_t Encoder;
	struct FileHeader_t FileHeader;

	//! Check arguments
	if(argc < 4) {
		printf(
			"ulcEncodeTool - Ultra-Low Complexity Codec Encoding Tool\n"
			"Usage:\n"
			" ulcencodetool Input.wav Output.ulc RateKbps[,AvgComplexity]|-Quality [Opt]\n"
			"Options:\n"
			" -blocksize:2048 - Set number of coefficients per block (must be a power of 2).\n"
			"Passing AvgComplexity uses ABR mode.\n"
			"Passing negative RateKbps (-Quality) uses VBR mode.\n"
			"Input file must be 8-bit, 16-bit, 24-bit, 32-bit, or 32-bit float.\n"
		);
		return 1;
	}

	//! Parse arguments
	int   BlockSize = 2048;
	float RateKbps;
	float AvgComplexity = 0.0f;
	sscanf(argv[3], "%f,%f", &RateKbps, &AvgComplexity);
	if(RateKbps == 0.0f) {
		printf("ERROR: Invalid coding rate (%.2f).\n", RateKbps);
		ExitCode = -1; goto Exit_BadArgs;
	}
	if(AvgComplexity < 0.0f) {
		printf("ERROR: Invalid AvgComplexity parameter (%.2f).\n", AvgComplexity);
		ExitCode = -1; goto Exit_BadArgs;
	}
	{
		int n;
		for(n=4;n<argc;n++) {
			if(!memcmp(argv[n], "-blocksize:", 11)) {
				int x = atoi(argv[n] + 11);
				if(x >= 256 && x <= 32768 && (x & (-x)) == x) BlockSize = x;
				else {
					printf("ERROR: Unsupported block size (%d).\n", x);
					ExitCode = -1; goto Exit_BadArgs;
				}
			}

			else printf("WARNING: Ignoring unknown argument (%s).\n", argv[n]);
		}
	}

	//! Open input file and verify
	{
		int Error = WAV_OpenR(&FileIn, argv[1]);
		if(Error < 0) {
			printf("ERROR: Unable to open input file (%s); error %s.\n", argv[1], WAV_ErrorCodeToString(Error));
			ExitCode = -1; goto Exit_FailOpenInFile;
		}
	}
	if(FileIn.fmt->nSamplesPerSec < 1) {
		printf("ERROR: Unsupported playback rate (%u).\n", FileIn.fmt->nSamplesPerSec);
		ExitCode = -1; goto Exit_FailInFileValidation;
	}
	if(FileIn.fmt->nChannels < 1/* || FileIn.fmt->nChannels > 0xFFFF*/) {
		printf("ERROR: Unsupported number of channels (%u).\n", FileIn.fmt->nChannels);
		ExitCode = -1; goto Exit_FailInFileValidation;
	}

	//! Allocate reading buffer
	AllocBuffer = malloc(BUFFER_ALIGNMENT-1 + sizeof(float)*BlockSize*FileIn.fmt->nChannels);
	if(!AllocBuffer) {
		printf("ERROR: Couldn't allocate reading buffer.\n");
		ExitCode = -1; goto Exit_FailCreateAllocBuffer;
	}
	float *ReadBuffer = (float*)(AllocBuffer + (-(uintptr_t)AllocBuffer % BUFFER_ALIGNMENT));

	//! Create file header
	//! nBlocks is +1 to account for coding delay, +1 to account for MDCT delay
	//! ::RateKbps and ::StreamOffs are written later
	FileHeader.Magic        = HEADER_MAGIC;
	FileHeader.BlockSize    = BlockSize;
	FileHeader.MaxBlockSize = 0;
	FileHeader.nBlocks      = (FileIn.nSamplePoints + BlockSize-1) / BlockSize + 2;
	FileHeader.RateHz       = FileIn.fmt->nSamplesPerSec;
	FileHeader.nChan        = FileIn.fmt->nChannels;

	//! Create encoder
	Encoder.RateHz    = FileHeader.RateHz;
	Encoder.nChan     = FileHeader.nChan;
	Encoder.BlockSize = FileHeader.BlockSize;
	if(ULC_EncoderState_Init(&Encoder) <= 0) {
		printf("ERROR: Unable to initialize encoder.\n");
		ExitCode = -1; goto Exit_FailCreateEncoder;
	}

	//! Open output file and skip header
	FileOut = fopen(argv[2], "wb");
	if(!FileOut) {
		printf("ERROR: Unable to open output file (%s).\n", argv[2]);
		ExitCode = -1; goto Exit_FailOpenFileOut;
	}
	size_t FileHeaderOffs = ftell(FileOut);
	fseek(FileOut, +sizeof(FileHeader), SEEK_CUR);

	//! Begin encoding
	{
		const clock_t DISPLAY_UPDATE_RATE = (clock_t)(CLOCKS_PER_SEC * 0.5); //! Update every 0.5 seconds

		//! Store stream offset
		FileHeader.StreamOffs = ftell(FileOut);

		//! Process blocks
		size_t Blk, nBlk = FileHeader.nBlocks;
		uint64_t TotalSize = 0;
		double ComplexitySum = 0.0;
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
					"\rBlock %zu/%zu (%.2f%% | %.2f X rt) | Average: %.2fkbps",
					Blk, nBlk, Blk*100.0/nBlk,
					nBlkProcessed*BlockSize / (double)FileHeader.RateHz,
					Blk ? (TotalSize * 8.0 * FileHeader.RateHz/1000.0 / (Blk * BlockSize)) : 0.0f
				);
				fflush(stdout);
				LastUpdateTime += DISPLAY_UPDATE_RATE;
				BlkLastUpdate   = Blk;
			}

			//! Read samples
			WAV_ReadAsFloat(&FileIn, ReadBuffer, BlockSize);

			//! Encode block
			int Size;
			const uint8_t *EncData;
			     if(RateKbps      < 0.0f) EncData = ULC_EncodeBlock_VBR(&Encoder, ReadBuffer, &Size, -RateKbps);
			else if(AvgComplexity > 0.0f) EncData = ULC_EncodeBlock_ABR(&Encoder, ReadBuffer, &Size,  RateKbps, AvgComplexity);
			else                          EncData = ULC_EncodeBlock_CBR(&Encoder, ReadBuffer, &Size,  RateKbps);

			//! Convert size to bytes and accumulate statistics
			Size = (Size+7) / 8u;
			TotalSize     += Size;
			ComplexitySum += Encoder.BlockComplexity;
			if((size_t)Size > FileHeader.MaxBlockSize) FileHeader.MaxBlockSize = Size;

			//! Write block to file
			fwrite(EncData, sizeof(uint8_t), Size, FileOut);
		}

		//! Show statistics and store RateKbps to header
		size_t nEncodedSamples = BlockSize * nBlk;
		double TotalSizeKiB  = TotalSize               * 1.0 / 1024;
		double AvgKbps       = TotalSize               * 8.0 * FileHeader.RateHz/1000.0 / nEncodedSamples;
		double AvgBitsPerSmp = TotalSize               * 8.0 / nEncodedSamples;
		double MaxKbps       = FileHeader.MaxBlockSize * 8.0 * FileHeader.RateHz/1000.0 / BlockSize;
		double MaxBitsPerSmp = FileHeader.MaxBlockSize * 8.0 / BlockSize;
		double Complexity    = ComplexitySum / nBlk;
		printf(
			"\n"
			"Total size = %.2fKiB\n"
			"Avg rate = %.5fkbps (%.5f bits/sample)\n"
			"Max rate = %.5fkbps (%.5f bits/sample)\n"
			"Avg complexity = %.5f\n",
			TotalSizeKiB,
			AvgKbps, AvgBitsPerSmp,
			MaxKbps, MaxBitsPerSmp,
			Complexity
		);
		FileHeader.RateKbps = lrint(AvgKbps);
	}

	//! Write file header
	fseek(FileOut, FileHeaderOffs, SEEK_SET);
	fwrite(&FileHeader, sizeof(FileHeader), 1, FileOut);

	//! Exit points
	fclose(FileOut);
Exit_FailOpenFileOut:
	ULC_EncoderState_Destroy(&Encoder);
Exit_FailCreateEncoder:
	free(AllocBuffer);
Exit_FailCreateAllocBuffer:
Exit_FailInFileValidation:
	WAV_Close(&FileIn);
Exit_FailOpenInFile:
Exit_BadArgs:
	return ExitCode;
}

/**************************************/
//! EOF
/**************************************/
