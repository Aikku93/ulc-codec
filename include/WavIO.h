/**************************************/
#pragma once
/**************************************/
#include <stdint.h>
#include <stdio.h>
/**************************************/

//! WAV_Open() error codes
#define WAV_ENOFILE      (-1)
#define WAV_ENOMEM       (-2)
#define WAV_ENOTWAV      (-3)
#define WAV_EINVALID     (-4)
#define WAV_EUNSUPPORTED (-5)
#define WAV_EIO          (-6)

/**************************************/

//! WAVE -> fmt
#define WAVE_FORMAT_PCM        0x0001
#define WAVE_FORMAT_IEEE_FLOAT 0x0003
struct WAVE_fmt_t {
	uint16_t wFormatTag;
	uint16_t nChannels;
	uint32_t nSamplesPerSec;
	uint32_t nAvgBytesPerSec;
	uint16_t nBlockAlign;
	uint16_t wBitsPerSample;
};

/**************************************/

//! WAV chunk descriptor type
struct WAV_Chunk_t {
	uint32_t CkType;
	uint32_t CkSize;
	size_t   FileOffs;
	struct WAV_Chunk_t *Prev, *Next;
};

/**************************************/

//! Internal state type
struct WAV_State_t {
	FILE *File;
	uint8_t  Mode;
	uint32_t SamplePosition;
	uint32_t nSamplePoints;
	struct WAVE_fmt_t  *fmt;
	struct WAV_Chunk_t *dataCk;
	struct WAV_Chunk_t *Chunks;
};

/**************************************/

//! WAV_ErrorCodeToString(e)
//! Description: Convert error code to string.
//! Arguments:
//!  e: Error code.
//! Returns: Error code name (eg. "WAV_ENOFILE").
//! Notes:
//!  -Returns NULL on invalid error code.
const char *WAV_ErrorCodeToString(int e);

/**************************************/

//! WAV_OpenR(WavState, Filename)
//! Description: Open WAV file for reading.
//! Arguments:
//!   WavState: Structure to store internal state in.
//!   Filename: File to open.
//! Returns:
//!   On success, returns 0. On failure, returns a value < 0, corresponding to
//!   the error codes at the start of this file.
//! Notes:
//!  -Only the `fmt` and `data` chunks are parsed.
//!  -Only PCM8, PCM16, PCM24, PCM32, and FLOAT32 are supported formats.
int WAV_OpenR(struct WAV_State_t *WavState, const char *Filename);

//! WAV_ReadAsFloat(WavState, Dst, nSmpPoints)
//! Description: Read samples from file into float-type buffer.
//! Arguments:
//!   WavState:   Structure holding the internal state.
//!   Dst:        Memory into which to store samples.
//!   nSmpPoints: Number of sample points to read.
//! Returns: The number of sample points read from file.
//! Notes:
//!  -If EOF occurs before reading all requested samples, any remaining samples
//!   will be filled with silence rather than being left untouched.
//!  -Channels are interleaved as output.
uint32_t WAV_ReadAsFloat(struct WAV_State_t *WavState, float *Dst, uint32_t nSmpPoints);

/**************************************/

//! WAV_OpenW(WavState, Filename)
//! Description: Open WAV file for reading.
//! Arguments:
//!   WavState: Structure to store internal state in.
//!   Filename: File to open.
//!   fmt:      WAV format.
//! Returns:
//!   On success, returns 0. On failure, returns a value < 0, corresponding to
//!   the error codes at the start of this file.
//! Notes:
//!  -Only PCM8, PCM16, PCM24, PCM32, and FLOAT32 are supported formats.
//!  -The `fmt` header is copied locally.
int WAV_OpenW(struct WAV_State_t *WavState, const char *Filename, const struct WAVE_fmt_t *fmt);

//! WAV_WriteFromFloat(WavState, Src, nSmpPoints)
//! Description: Write samples to file from float-type buffer.
//! Arguments:
//!   WavState:   Structure holding the internal state.
//!   Src:        Memory from which to load samples.
//!   nSmpPoints: Number of sample points to write.
//! Returns: The number of sample points written to file.
//! Notes:
//!  -Internally, samples may be buffered; if the samples MUST be output, use
//!   the WAV_Flush() function following this call.
//!  -Channels must be interleaved as input.
int WAV_WriteFromFloat(struct WAV_State_t *WavState, const float *Src, uint32_t nSmpPoints);

//! WAV_Flush(WavState)
//! Description: Flush any buffered samples to file.
//! Arguments:
//!   WavState: Structure holding the internal state.
//! Returns: Nothing; samples flushed from buffer.
inline void WAV_Flush(struct WAV_State_t *WavState) {
	fflush(WavState->File);
}

/**************************************/

//! WAV_Close(WavState)
//! Description: Close WAV file.
//! Arguments:
//!   WavState: Structure holding the internal state.
//! Returns: Nothing; file is closed.
void WAV_Close(struct WAV_State_t *WavState);

/**************************************/
//! EOF
/**************************************/
