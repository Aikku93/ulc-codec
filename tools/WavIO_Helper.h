/**************************************/
#pragma once
/**************************************/
#include <stdint.h>
/**************************************/
#include "WavIO.h"
/**************************************/

//! WAV_State_t::Mode
#define WAV_STATE_MODE_READ  0
#define WAV_STATE_MODE_WRITE 1

/**************************************/

//! Convert data to normalized float
void WAV_ConvertToFloat_PCM8u(float *Dst, const void *RawMem, uint32_t N);
void WAV_ConvertToFloat_PCM16(float *Dst, const void *RawMem, uint32_t N);
void WAV_ConvertToFloat_PCM24(float *Dst, const void *RawMem, uint32_t N);

//! Convert data from normalized float
void WAV_ConvertFromFloat_PCM8u(void *RawDst, const float *Src, uint32_t N);
void WAV_ConvertFromFloat_PCM16(void *RawDst, const float *Src, uint32_t N);
void WAV_ConvertFromFloat_PCM24(void *RawDst, const float *Src, uint32_t N);

/**************************************/

//! Append new chunk to list
//! Returns pointer to allocated chunk
struct WAV_Chunk_t *WAV_AppendCkHeader(struct WAV_State_t *WavState, size_t ExtraData);

/**************************************/
//! EOF
/**************************************/
