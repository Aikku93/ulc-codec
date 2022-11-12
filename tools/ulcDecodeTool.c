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
#include "ulcDecoder.h"
#include "WavIO.h"
/**************************************/

//! Possible output formats
#define FORMAT_PCM8    0
#define FORMAT_PCM16   1
#define FORMAT_PCM24   2
#define FORMAT_FLOAT32 3

/**************************************/

int main(int argc, const char *argv[]) {
	int   ExitCode = 0;
	FILE *FileIn;
	char *AllocBuffer;
	struct WAV_State_t FileOut;
	struct ULC_DecoderState_t Decoder;
	struct FileHeader_t FileHeader;

	//! Check arguments
	if(argc < 3) {
		printf(
			"ulcDecodeTool - Ultra-Low Complexity Codec Decoding Tool\n"
			"Usage: ulcdecodetool Input.ulc Output.wav [Opt]\n"
			"Options:\n"
			" -format:PCM16 - Set output format (PCM8, PCM16, PCM24, FLOAT32).\n"
		);
		return 1;
	}

	//! Parse arguments
	int FormatType = FORMAT_PCM16;
	{
		int n;
		for(n=3;n<argc;n++) {
			if(!memcmp(argv[n], "-format:", 8)) {
				const char *FmtStr = argv[n] + 8;
				     if(!strcmp(FmtStr, "PCM8")    || !strcmp(FmtStr, "pcm8"))
					FormatType = FORMAT_PCM8;
				else if(!strcmp(FmtStr, "PCM16")   || !strcmp(FmtStr, "pcm16"))
					FormatType = FORMAT_PCM16;
				else if(!strcmp(FmtStr, "PCM24")   || !strcmp(FmtStr, "pcm24"))
					FormatType = FORMAT_PCM24;
				else if(!strcmp(FmtStr, "FLOAT32") || !strcmp(FmtStr, "float32"))
					FormatType = FORMAT_FLOAT32;
				else {
					printf("ERROR: Ignoring invalid output format (%s).\n", FmtStr);
					ExitCode = -1; goto Exit_BadArgs;
				}
			}

			else printf("WARNING: Ignoring unknown argument (%s).\n", argv[n]);
		}
	}

	//! Open input file and verify
	FileIn = fopen(argv[1], "rb");
	if(!FileIn) {
		printf("ERROR: Unable to open input file (%s).\n", argv[1]);
		ExitCode = -1; goto Exit_FailOpenInFile;
	}
	if(fread(&FileHeader, sizeof(FileHeader), 1, FileIn) != 1 || FileHeader.Magic != HEADER_MAGIC) {
		printf("ERROR: Input file is not a valid ULC container.\n");
		ExitCode = -1; goto Exit_FailVerifyInFile;
	}

	//! Define the stream buffer size
	int StreamBufferSize = (16*1024);
	if((int)FileHeader.MaxBlockSize > StreamBufferSize) StreamBufferSize = FileHeader.MaxBlockSize;

	//! Allocate decoding buffer and stream buffer
	AllocBuffer = malloc(BUFFER_ALIGNMENT-1 + sizeof(float)*2*FileHeader.BlockSize*FileHeader.nChan + StreamBufferSize);
	if(!AllocBuffer) {
		printf("ERROR: Couldn't allocate decoding buffer.\n");
		ExitCode = -1; goto Exit_FailCreateAllocBuffer;
	}
	float   *DecodeBuffer = (float  *)(AllocBuffer + (-(uintptr_t)AllocBuffer % BUFFER_ALIGNMENT));
	uint8_t *StreamBuffer = (uint8_t*)(DecodeBuffer + FileHeader.BlockSize*FileHeader.nChan);

	//! Create decoder
	Decoder.nChan      = FileHeader.nChan;
	Decoder.BlockSize  = FileHeader.BlockSize;
	if(ULC_DecoderState_Init(&Decoder) <= 0) {
		printf("ERROR: Unable to initialize decoder.\n");
		ExitCode = -1; goto Exit_FailCreateDecoder;
	}

	//! Create output file
	{
		int BytesPerSmp = 0;
		switch(FormatType) {
			case FORMAT_PCM8:    BytesPerSmp =  8 / 8; break;
			case FORMAT_PCM16:   BytesPerSmp = 16 / 8; break;
			case FORMAT_PCM24:   BytesPerSmp = 24 / 8; break;
			case FORMAT_FLOAT32: BytesPerSmp = 32 / 8; break;
		}
		struct WAVE_fmt_t fmt;
		fmt.wFormatTag      = (FormatType == FORMAT_FLOAT32) ? WAVE_FORMAT_IEEE_FLOAT : WAVE_FORMAT_PCM;
		fmt.nChannels       = FileHeader.nChan;
		fmt.nSamplesPerSec  = FileHeader.RateHz;
		fmt.nAvgBytesPerSec = BytesPerSmp * FileHeader.nChan * FileHeader.RateHz;
		fmt.nBlockAlign     = BytesPerSmp * FileHeader.nChan;
		fmt.wBitsPerSample  = BytesPerSmp * 8;
		int Error = WAV_OpenW(&FileOut, argv[2], &fmt);
		if(Error < 0) {
			printf("ERROR: Unable to create output file (%s); error %s.", argv[1], WAV_ErrorCodeToString(Error));
			ExitCode = -1; goto Exit_FailCreateOutFile;
		}
	}

	//! Begin decoding
	{
		const clock_t DISPLAY_UPDATE_RATE = (clock_t)(CLOCKS_PER_SEC * 0.5); //! Update every 0.5 seconds

		//! Pre-fill the streaming buffer
		fseek(FileIn, FileHeader.StreamOffs, SEEK_SET);
		fread(StreamBuffer, StreamBufferSize, 1, FileIn);

		//! Process blocks
		int      BlockSize   = FileHeader.BlockSize;
		uint32_t Blk, nBlk = FileHeader.nBlocks;
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
					nBlkProcessed*BlockSize / (double)FileHeader.RateHz
				);
				fflush(stdout);
				LastUpdateTime += DISPLAY_UPDATE_RATE;
				BlkLastUpdate   = Blk;
			}

			//! Decode block
			int Size = (ULC_DecodeBlock(&Decoder, DecodeBuffer, StreamBuffer) + 7) / 8u;
			if(!Size) {
				printf("ERROR: Corrupted stream.\n");
				ExitCode = -1; goto Exit_FailCorruptStream;
			}

			//! Write samples
			WAV_WriteFromFloat(&FileOut, DecodeBuffer, BlockSize);

			//! Slide stream buffer
			memcpy(StreamBuffer, StreamBuffer+Size, StreamBufferSize-Size);
			fread(StreamBuffer + StreamBufferSize-Size, Size, 1, FileIn);
		}
	}

	//! Exit points
	printf("\nOk\n");
Exit_FailCorruptStream:
	WAV_Close(&FileOut);
Exit_FailCreateOutFile:
	ULC_DecoderState_Destroy(&Decoder);
Exit_FailCreateDecoder:
	free(AllocBuffer);
Exit_FailCreateAllocBuffer:
Exit_FailVerifyInFile:
	fclose(FileIn);
Exit_FailOpenInFile:
Exit_BadArgs:
	return ExitCode;
}

/**************************************/
//! EOF
/**************************************/
