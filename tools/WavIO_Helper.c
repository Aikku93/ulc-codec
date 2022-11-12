/**************************************/
#include <math.h>
#include <stdlib.h>
/**************************************/
#include "WavIO.h"
#include "WavIO_Helper.h"
/**************************************/

static inline float Clamp(float x, float Min, float Max) {
	return (x < Min) ? Min : (x > Max) ? Max : x;
}

/**************************************/

const char *WAV_ErrorCodeToString(int e) {
	switch(e) {
#define CREATECASE(x) case x: return #x;
		CREATECASE(WAV_ENOFILE);
		CREATECASE(WAV_ENOMEM);
		CREATECASE(WAV_ENOTWAV);
		CREATECASE(WAV_EINVALID);
		CREATECASE(WAV_EUNSUPPORTED);
		CREATECASE(WAV_EIO);
#undef CREATECASE
	}
	return NULL;
}

/**************************************/

void WAV_ConvertToFloat_PCM8u(float *Dst, const void *RawMem, uint32_t N) {
	uint32_t n;
	const int8_t *Src = (const int8_t*)RawMem;
	for(n=0;n<N;n++) {
		*Dst++ = (float)(*Src++ ^ 0x80) * 0x1.0p-7f;
	}
}

void WAV_ConvertFromFloat_PCM8u(void *RawDst, const float *Src, uint32_t N) {
	uint32_t n;
	int8_t *Dst = (int8_t*)RawDst;
	for(n=0;n<N;n++) {
		*Dst++ = (int8_t)lrintf(Clamp(*Src++ * 0x1.0p+7f, (float)-0x80, (float)+0x7F)) ^ 0x80;
	}
}

/**************************************/

void WAV_ConvertToFloat_PCM16(float *Dst, const void *RawMem, uint32_t N) {
	uint32_t n;
	const int16_t *Src = (const int16_t*)RawMem;
	for(n=0;n<N;n++) {
		*Dst++ = (float)(*Src++) * 0x1.0p-15f;
	}
}

void WAV_ConvertFromFloat_PCM16(void *RawDst, const float *Src, uint32_t N) {
	uint32_t n;
	int16_t *Dst = (int16_t*)RawDst;
	for(n=0;n<N;n++) {
		*Dst++ = (int16_t)lrintf(Clamp(*Src++ * 0x1.0p+15f, (float)-0x8000, (float)+0x7FFF));
	}
}

/**************************************/

void WAV_ConvertToFloat_PCM24(float *Dst, const void *RawMem, uint32_t N) {
	uint32_t n;
	const uint8_t *Src = (const uint8_t*)RawMem;
	for(n=0;n<N;n++) {
		int32_t x  = (int32_t)(*Src++) <<  8;
		        x |= (int32_t)(*Src++) << 16;
		        x |= (int32_t)(*Src++) << 24;
		*Dst++ = (float)x * 0x1.0p-31f;
	}
}

void WAV_ConvertFromFloat_PCM24(void *RawDst, const float *Src, uint32_t N) {
	uint32_t n;
	uint8_t *Dst = (uint8_t*)RawDst;
	for(n=0;n<N;n++) {
		uint32_t x = (uint32_t)lrintf(Clamp(*Src++ * 0x1.0p+23f, (float)-0x800000, (float)+0x7FFFFF));
		*Dst++ = (int8_t)(x >> 0);
		*Dst++ = (int8_t)(x >> 8);
		*Dst++ = (int8_t)(x >> 16);
	}
}

/**************************************/

//! Append new chunk to list
struct WAV_Chunk_t *WAV_AppendCkHeader(struct WAV_State_t *WavState, size_t ExtraData) {
	//! Attempt to allocate chunk
	struct WAV_Chunk_t *Ck = malloc(sizeof(struct WAV_Chunk_t) + ExtraData);
	if(!Ck) return NULL;

	//! Append to chunks list
	struct WAV_Chunk_t *Prev = NULL;
	struct WAV_Chunk_t *Next = WavState->Chunks;
	while(Next) {
		Prev = Next;
		Next = Next->Next;
	}
	Ck->Prev = Prev;
	Ck->Next = Next;
	if(Prev) Prev->Next = Ck; else WavState->Chunks = Ck;
	if(Next) Next->Prev = Ck;

	//! Return pointer to chunk
	return Ck;
}

/**************************************/

void WAV_Close(struct WAV_State_t *WavState) {
	FILE *f = WavState->File;
	if(WavState->Mode == WAV_STATE_MODE_READ) {
		//! Clean up all chunks
		struct WAV_Chunk_t *Ck = WavState->Chunks;
		while(Ck) {
			struct WAV_Chunk_t *Next = Ck->Next;
			free(Ck);
			Ck = Next;
		}
	} else {
		//! Finish up the data chunk
		uint32_t dataSize = ftell(f) - (12 + 8+16 + 8); //! 12 bytes for RIFF-WAVE, 8+16 bytes for fmt, 8 bytes for data header
		fseek(f, 12 + 8+16 + 4, SEEK_SET); //! dataCk.Size
		fwrite(&dataSize, sizeof(uint32_t), 1, f);
		fseek(f, 0, SEEK_END);
		if(ftell(f) & 1) fputc(0, f); //! <- Align chunk to 2 bytes

		//! Append any additional chunks
		struct WAV_Chunk_t *Ck = WavState->Chunks;
		while(Ck) {
			fwrite(Ck, sizeof(uint32_t)*2, 1, f); //! Write CkType,CkSize
			fwrite(Ck+1, Ck->CkSize, 1, f);
			if(ftell(f) & 1) fputc(0, f); //! <- Align chunk to 2 bytes
			Ck = Ck->Next;
		}

		//! Write final RIFF size
		uint32_t RIFFSize = ftell(f) - 8;
		fseek(f, 0+4, SEEK_SET);
		fwrite(&RIFFSize, sizeof(uint32_t), 1, f);
	}
	fclose(f);
}

/**************************************/
//! EOF
/**************************************/
