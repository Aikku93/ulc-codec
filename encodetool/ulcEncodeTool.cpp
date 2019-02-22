/**************************************/
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
/**************************************/
#include "ulcCoder.h"
/**************************************/
using namespace std;
using namespace ULC;
/**************************************/

//! Block data buffer
__attribute__((aligned(32))) static int16_t BlockFetch [MAX_CHANS*BLOCK_SIZE];
__attribute__((aligned(32))) static float   BlockBuffer[MAX_CHANS*BLOCK_SIZE];

/**************************************/

//! Nybble coding class
class NybCoder_c {
	using _DatT = uint32_t;

	static constexpr size_t DUMP_SIZE = (8*1024) / sizeof(_DatT); //! 8KiB
	static constexpr size_t NYBS_PER_DATA = sizeof(_DatT)*8 / 4;

	size_t DataIdx;
	FILE *DataOutput;
	_DatT DataStore[DUMP_SIZE];

	void ClearData() {
		DataIdx = 0;
		//for(size_t i=0;i<DUMP_SIZE;i++) DataStore[i] = 0;
	}
public:
	void FinalizeBlock() {
		//! Shift down last element if needed
		if(size_t n = DataIdx%NYBS_PER_DATA) DataStore[DataIdx/NYBS_PER_DATA] >>= (NYBS_PER_DATA - n)*4;

		//! Dump to file
		if(size_t n = (DataIdx+NYBS_PER_DATA-1)/NYBS_PER_DATA) {
			fwrite(DataStore, sizeof(_DatT), n, DataOutput);
		}

		//! Reset
		ClearData();
	}

	void Write(int x) {
		_DatT &v = DataStore[DataIdx/NYBS_PER_DATA];
		v = (v>>4) | (x<<(sizeof(_DatT)*8 - 4));
		if(++DataIdx == DUMP_SIZE*NYBS_PER_DATA) FinalizeBlock();
	}

	NybCoder_c(FILE *Out) : DataOutput(Out) {
		ClearData();
	}
	~NybCoder_c() {
		FinalizeBlock();
	}
};

/**************************************/

//! Nybble-writing callback
static void NybWriteCbFnc(int x, void *User) {
	reinterpret_cast<NybCoder_c*>(User)->Write(x);
}

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
	if(nChan < 1 || nChan > MAX_CHANS) {
		printf("ERROR: Invalid number of channels.\n");
		return -1;
	}

	//! Open input file
	size_t nSmp;
	FILE *InFile = fopen(argv[1], "rb");
	if(!InFile) {
		printf("ERROR: Unable to open input file.\n");
		return -1;
	} else {
		fseek(InFile, 0, SEEK_END);
		nSmp = ftell(InFile) / sizeof(int16_t) / nChan;
		rewind(InFile);
	}

	//! Open output file
	FILE *OutFile = fopen(argv[2], "wb");
	if(!OutFile) {
		printf("ERROR: Unable to open output file.\n");
		fclose(InFile);
		return -1;
	}

	/*! Write file header !*/ {
		struct {
			uint32_t Magic;
			uint32_t nSamp;
			uint32_t RateHz;
			uint16_t RateKbps;
			uint16_t nChan;
		} Header = {
			uint32_t('U' | 'L'<<8 | 'C'<<16 | 'x'<<24),
			nSmp,
			RateHz,
			uint16_t(RateKbps),
			uint16_t(nChan)
		};
		fwrite(&Header, sizeof(Header), 1, OutFile);
	}

	//! Create encoder class, nybble writer
	Encoder_c  Encoder(RateHz, nChan);
	NybCoder_c NybWriter(OutFile);

	//! Process blocks
	size_t nBlk = (nSmp + BLOCK_SIZE-1) / BLOCK_SIZE;
	uint64_t TotalSize = 0;
	for(size_t Blk=0;Blk<nBlk+1;Blk++) { //! +1 to account for coding delay
		//! Show progress
		printf("\rBlock %u/%u (%.2f%%)...", Blk, nBlk, Blk*100.0f/nBlk);
		fflush(stdout);

		//! Fill buffer data
		fread(BlockFetch, nChan*sizeof(int16_t), BLOCK_SIZE, InFile);
		for(size_t Chan=0;Chan<nChan;Chan++) for(size_t n=0;n<BLOCK_SIZE;n++) {
			BlockBuffer[Chan*BLOCK_SIZE+n] = (Blk*BLOCK_SIZE+n < nSmp) ? BlockFetch[n*nChan+Chan] : 0.0f;
		}

		//! Apply M/S transform
		if(nChan == 2) for(size_t n=0;n<BLOCK_SIZE;n++) {
			float &a = BlockBuffer[0*BLOCK_SIZE+n], va = a;
			float &b = BlockBuffer[1*BLOCK_SIZE+n], vb = b;
			a = (va + vb) * 0.5f;
			b = (va - vb) * 0.5f;
		}

		//! Encode block
		TotalSize += Encoder.EncodeBlock(BlockBuffer, RateKbps, NybWriteCbFnc, &NybWriter);
	}

	//! Explicitly finalize the last block
	//! This is needed because by the time the destructor is called,
	//! the file handle will have been closed
	NybWriter.FinalizeBlock();

	//! Show statistics
	printf(
		"\e[2K\r"
		"Total size = %.2fKiB\n"
		"Avg rate = %.5fkbps (%.5f bits/sample)\n",
		TotalSize/8.0 / 1024,
		TotalSize*1.0 * RateHz/1000.0 / nSmp,
		TotalSize*1.0 / nSmp
	);

	//! Clean up
	fclose(OutFile);
	fclose(InFile);
	return 0;
}

/**************************************/
//! EOF
/**************************************/
