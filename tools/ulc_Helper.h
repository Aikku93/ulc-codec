/**************************************/
#pragma once
/**************************************/
#include <stdint.h>
/**************************************/
#define BUFFER_ALIGNMENT 64u //! __mm512
/**************************************/

//! File header
#define HEADER_MAGIC (uint32_t)('U' | 'L'<<8 | 'C'<<16 | '2'<<24)
struct FileHeader_t {
	uint32_t Magic;        //! [00h] Magic value/signature
	uint16_t BlockSize;    //! [04h] Transform block size
	uint16_t MaxBlockSize; //! [06h] Largest block size (in bytes; 0 = Unknown)
	uint32_t nBlocks;      //! [08h] Number of blocks
	uint32_t RateHz;       //! [0Ch] Playback rate
	uint16_t nChan;        //! [10h] Channels in stream
	uint16_t RateKbps;     //! [12h] Nominal coding rate
	uint32_t StreamOffs;   //! [14h] Offset of data stream
};

/**************************************/
//! EOF
/**************************************/
