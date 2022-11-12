/**************************************/
#include <stdlib.h>
/**************************************/
#include "MiniRIFF.h"
#include "WavIO.h"
#include "WavIO_Helper.h"
/**************************************/

//! RIFF(WAVE) -> data
static int RIFF_WAVE_data(FILE *f, void *User, const struct RIFF_CkHeader_t *Ck) {
	struct WAV_State_t *WavState = (struct WAV_State_t*)User;

	//! Append chunk
	struct WAV_Chunk_t *WavCk = WAV_AppendCkHeader(WavState, sizeof(struct WAVE_fmt_t));
	if(!WavCk) return WAV_ENOMEM;
	WavCk->CkType    = Ck->Type;
	WavCk->CkSize    = Ck->Size;
	WavCk->FileOffs  = ftell(f);
	WavState->dataCk = WavCk;
	return 0;
}

/**************************************/

//! RIFF(WAVE) -> fmt
static int RIFF_WAVE_fmt(FILE *f, void *User, const struct RIFF_CkHeader_t *Ck) {
	struct WAV_State_t *WavState = (struct WAV_State_t*)User;

	//! Append a header plus the size of the `fmt` chunk
	struct WAV_Chunk_t *WavCk = WAV_AppendCkHeader(WavState, sizeof(struct WAVE_fmt_t));
	if(!WavCk) return WAV_ENOMEM;

	//! Read the chunk
	WavCk->CkType   = Ck->Type;
	WavCk->CkSize   = Ck->Size;
	WavCk->FileOffs = ftell(f);
	struct WAVE_fmt_t *fmt = (struct WAVE_fmt_t*)(WavCk+1);
	if(fread(fmt, sizeof(struct WAVE_fmt_t), 1, f) != 1) return WAV_EIO;
	WavState->fmt = fmt;
	return 0;
}

/**************************************/
//! WAVE chunk hierarchy (read bottom to top)
/**************************************/

//! RIFF(WAVE) -> Chunks
static const struct RIFF_CkHdl_t RIFF_WAVE_Ck[] = {
	{RIFF_FOURCC("fmt "),RIFF_WAVE_fmt},
	{RIFF_FOURCC("data"),RIFF_WAVE_data},
	{0},
};

//! RIFF(WAVE)
static const struct RIFF_CkListHdl_t RIFF_WAVE[] = {
	{RIFF_FOURCC("WAVE"),RIFF_WAVE_Ck,NULL,NULL,NULL},
	{0},
};

/**************************************/

int WAV_OpenR(struct WAV_State_t *WavState, const char *Filename) {
	//! Attempt to open file
	FILE *f = fopen(Filename, "rb");
	if(!f) return WAV_ENOFILE;

	//! Map out the RIFF structure
	WavState->Chunks = NULL;
	int RetVal = RIFF_CkRead(f, WavState, NULL, RIFF_WAVE);
	if(RetVal < 0) {
		fclose(f);
		return RetVal;
	}

	//! Check to see if the format is supported
	//! Check to see if this format is supported
	struct WAVE_fmt_t *fmt = WavState->fmt;
	if(!fmt->nChannels) return WAV_EINVALID;
	WavState->nSamplePoints = WavState->dataCk->CkSize / fmt->nChannels;
	switch(fmt->wFormatTag) {
		case WAVE_FORMAT_PCM: {
			     if(fmt->wBitsPerSample == 8) {
				WavState->nSamplePoints /= sizeof(uint8_t);
			}
			else if(fmt->wBitsPerSample == 16) {
				WavState->nSamplePoints /= sizeof(int16_t);
			}
			else if(fmt->wBitsPerSample == 24) {
				WavState->nSamplePoints /= 3; //! Yeah, looks ugly and out of place
			}
			else {
				fclose(f);
				return WAV_EUNSUPPORTED;
			}
		} break;
		case WAVE_FORMAT_IEEE_FLOAT: {
			if(fmt->wBitsPerSample == 32) {
				WavState->nSamplePoints /= sizeof(float);
			}
			else {
				fclose(f);
				return WAV_EUNSUPPORTED;
			}
		} break;
	}

	WavState->File           = f;
	WavState->Mode           = WAV_STATE_MODE_READ;
	WavState->SamplePosition = 0;
	return 0;
}

/**************************************/

uint32_t WAV_ReadAsFloat(struct WAV_State_t *WavState, float *Dst, uint32_t nSmpPoints) {
	struct WAVE_fmt_t  *fmt    = WavState->fmt;
	struct WAV_Chunk_t *dataCk = WavState->dataCk;
	uint32_t SmpPointSize = (fmt->wBitsPerSample/8)*fmt->nChannels;

	//! Seek to next sample
	size_t Pos = dataCk->FileOffs + WavState->SamplePosition*SmpPointSize;
	fseek(WavState->File, Pos, SEEK_SET);

	//! Read data to the end of the target memory to allow unpacking
	void *RawMem = (void*)((uintptr_t)(Dst + nSmpPoints*fmt->nChannels) - nSmpPoints*SmpPointSize);
	uint32_t nSmpPointsRead = fread(RawMem, SmpPointSize, nSmpPoints, WavState->File);
	if(WavState->SamplePosition+nSmpPointsRead > WavState->nSamplePoints) {
		nSmpPointsRead = WavState->nSamplePoints - WavState->SamplePosition;
	}

	//! Convert data as needed
	if(fmt->wFormatTag != WAVE_FORMAT_IEEE_FLOAT) {
		       if(fmt->wFormatTag == WAVE_FORMAT_PCM && fmt->wBitsPerSample == 8) {
			WAV_ConvertToFloat_PCM8u(Dst, RawMem, nSmpPointsRead*fmt->nChannels);
		} else if(fmt->wFormatTag == WAVE_FORMAT_PCM && fmt->wBitsPerSample == 16) {
			WAV_ConvertToFloat_PCM16(Dst, RawMem, nSmpPointsRead*fmt->nChannels);
		} else if(fmt->wFormatTag == WAVE_FORMAT_PCM && fmt->wBitsPerSample == 24) {
			WAV_ConvertToFloat_PCM24(Dst, RawMem, nSmpPointsRead*fmt->nChannels);
		}
	}

	//! Pad with 0
	uint32_t n, N = (nSmpPoints - nSmpPointsRead)*fmt->nChannels;
	float *ZeroPadTarget = Dst + nSmpPointsRead*fmt->nChannels;
	for(n=0;n<N;n++) *ZeroPadTarget++ = 0.0f;

	//! Advance and return number of samples read
	WavState->SamplePosition += nSmpPointsRead;
	return nSmpPointsRead;
}

/**************************************/
//! EOF
/**************************************/
