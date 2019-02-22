/**************************************/
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
__attribute__((aligned(32))) static int16_t BlockOutput[MAX_CHANS*BLOCK_SIZE];
__attribute__((aligned(32))) static float   BlockBuffer[MAX_CHANS*BLOCK_SIZE];

/**************************************/

//! Clip to 16bit range
template<typename T> static inline T Clip16(T x) {
	if(x < T(-0x8000)) return T(-0x8000);
	if(x > T(+0x7FFF)) return T(+0x7FFF);
	return x;
}

/**************************************/

//! Nybble reading class
//! TODO: Not exactly ideal. Will replace later with more streamlined format
class NybReader_c {
	using _DatT = uint32_t;

	static constexpr size_t CACHE_SIZE = (8*1024) / sizeof(_DatT); //! 8KiB
	static constexpr size_t NYBS_PER_DATA = sizeof(_DatT)*8 / 4;

	size_t DataIdx;
	FILE *DataInput;
	_DatT DataCache[CACHE_SIZE];
public:
	int Read() {
		if(DataIdx >= CACHE_SIZE*NYBS_PER_DATA) {
			DataIdx = 0;
			fread(DataCache, sizeof(_DatT), CACHE_SIZE, DataInput);
		}

		int v = DataCache[DataIdx/NYBS_PER_DATA] >> ((DataIdx%NYBS_PER_DATA) * 4);
		DataIdx++;
		return v;
	}

	NybReader_c(FILE *In) : DataIdx(~0U), DataInput(In) {}
};

/**************************************/

//! Nybble-reading callback
static int NybReadCbFnc(void *User) {
	return reinterpret_cast<NybReader_c*>(User)->Read();
}

int main(int argc, const char *argv[]) {
	struct {
		uint32_t Magic;
		uint32_t nSamp;
		uint32_t RateHz;
		uint16_t RateKbps;
		uint16_t nChan;
	} Header;

	//! Check arguments
	if(argc != 3) {
		printf(
			"ulcDecodeTool - Ultra-Low Complexity Codec Decoding Tool\n"
			"Usage: ulcdecodetool Input.ulc Output.sw\n"
			"Multi-channel data will be interleaved.\n"
		);
		return 1;
	}

	//! Open input file
	FILE *InFile = fopen(argv[1], "rb");
	if(!InFile) {
		printf("ERROR: Unable to open input file.\n");
		return -1;
	}

	//! Read header
	fread(&Header, sizeof(Header), 1, InFile);
	if(Header.Magic != uint32_t('U'  | 'L' <<8 | 'C' <<16 | 'x' <<24)) {
		printf("ERROR: Invalid file.\n");
		fclose(InFile);
		return -1;
	}

	//! Open output file
	FILE *OutFile = fopen(argv[2], "wb");
	if(!OutFile) {
		printf("ERROR: Unable to open output file.\n");
		fclose(InFile);
		return -1;
	}

	//! Create decoder class, nybble reader
	size_t nChan = Header.nChan;
	Decoder_c   Decoder(nChan);
	NybReader_c NybReader(InFile);

	//! Process blocks
	size_t nBlk = (Header.nSamp + BLOCK_SIZE-1) / BLOCK_SIZE;
	for(size_t Blk=0;Blk<nBlk+1;Blk++) { //! +1 to account for coding delay
		//! Show progress
		printf("\rBlock %u/%u (%.2f%%)...", Blk, nBlk, Blk*100.0f/nBlk);
		fflush(stdout);

		//! Decode block
		Decoder.DecodeBlock(BlockBuffer, NybReadCbFnc, &NybReader);

		//! Apply M/S transform
		if(nChan == 2) for(size_t n=0;n<BLOCK_SIZE;n++) {
			float &a = BlockBuffer[0*BLOCK_SIZE+n], va = a;
			float &b = BlockBuffer[1*BLOCK_SIZE+n], vb = b;
			a = va + vb;
			b = va - vb;
		}

		//! Interleave to output buffer
		for(size_t Chan=0;Chan<nChan;Chan++) for(size_t n=0;n<BLOCK_SIZE;n++) {
			BlockOutput[n*nChan+Chan] = Clip16(BlockBuffer[Chan*BLOCK_SIZE+n]);
		}

		//! Write to file
		fwrite(BlockOutput, nChan*sizeof(int16_t), BLOCK_SIZE, OutFile);
	}

	//! Done
	printf("\e[2K\rOk\n");

	//! Clean up
	fclose(OutFile);
	fclose(InFile);
	return 0;
}

/**************************************/
//! EOF
/**************************************/
