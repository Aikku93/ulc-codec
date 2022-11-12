/**************************************/
#include <stdlib.h>
/**************************************/
#include "MiniRIFF.h"
#include "WavIO.h"
#include "WavIO_Helper.h"
/**************************************/
#define PACK_BUFFER_SIZE (64*1024)
/**************************************/

static uint8_t PackBuffer[PACK_BUFFER_SIZE];

/**************************************/

int WAV_OpenW(struct WAV_State_t *WavState, const char *Filename, const struct WAVE_fmt_t *fmt) {
	//! Create a local copy of the format
	struct WAVE_fmt_t *fmtCopy = malloc(sizeof(struct WAVE_fmt_t));
	if(!fmtCopy) return WAV_ENOMEM;
	*fmtCopy = *fmt;
	WavState->fmt = fmtCopy;

	//! Attempt to open file
	FILE *f = fopen(Filename, "wb");
	if(!f) {
		free(fmtCopy);
		return WAV_ENOFILE;
	}

	//! Write the RIFF WAVE header, with Size=0
	static const struct {
		uint32_t CkType;
		uint32_t CkSize;
		uint32_t Type;
	} RIFF_WAVE = {
		RIFF_FOURCC("RIFF"),
		0,
		RIFF_FOURCC("WAVE"),
	};
	fwrite(&RIFF_WAVE, sizeof(RIFF_WAVE), 1, f);

	//! Write the fmt chunk
	static const struct RIFF_CkHeader_t fmt_Header = {
		RIFF_FOURCC("fmt "),
		sizeof(struct WAVE_fmt_t),
	};
	fwrite(&fmt_Header, sizeof(fmt_Header), 1, f);
	fwrite(fmt, sizeof(struct WAVE_fmt_t), 1, f);

	//! Write the data chunk header, with Size=0
	static const struct RIFF_CkHeader_t data_Header = {
		RIFF_FOURCC("data "),
		0,
	};
	fwrite(&data_Header, sizeof(data_Header), 1, f);

	//! Set the initial state
	WavState->File           = f;
	WavState->Mode           = WAV_STATE_MODE_WRITE;
	WavState->SamplePosition = 0;
	WavState->Chunks         = NULL;
	return 0;
}

/**************************************/

int WAV_WriteFromFloat(struct WAV_State_t *WavState, const float *Src, uint32_t nSmpPoints) {
	struct WAVE_fmt_t *fmt = WavState->fmt;

	//! If we're already targetting 32bit float, just write directly
	uint32_t nTotalWriteSmp = 0;
	uint32_t SmpSize = fmt->wBitsPerSample/8;
	if(fmt->wFormatTag == WAVE_FORMAT_IEEE_FLOAT) {
		nTotalWriteSmp = fwrite(Src, sizeof(float)*fmt->nChannels, nSmpPoints, WavState->File);
	} else {
		//! While we have data to process, keep packing
		uint32_t nSmpRem = nSmpPoints*fmt->nChannels;
		uint32_t MaxProcessSmp = PACK_BUFFER_SIZE / SmpSize;
		while(nSmpRem) {
			//! Get number of samples to process
			uint32_t nProcessSmp = nSmpRem;
			if(nProcessSmp > MaxProcessSmp) nProcessSmp = MaxProcessSmp;
			nSmpRem -= nProcessSmp;

			//! Perform conversion
			       if(fmt->wFormatTag == WAVE_FORMAT_PCM && fmt->wBitsPerSample == 8) {
				WAV_ConvertFromFloat_PCM8u(PackBuffer, Src, nProcessSmp);
			} else if(fmt->wFormatTag == WAVE_FORMAT_PCM && fmt->wBitsPerSample == 16) {
				WAV_ConvertFromFloat_PCM16(PackBuffer, Src, nProcessSmp);
			} else if(fmt->wFormatTag == WAVE_FORMAT_PCM && fmt->wBitsPerSample == 24) {
				WAV_ConvertFromFloat_PCM24(PackBuffer, Src, nProcessSmp);
			}
			Src += nProcessSmp;

			//! Dump data to file
			uint32_t nWriteSmp = fwrite(PackBuffer, SmpSize, nProcessSmp, WavState->File);
			nTotalWriteSmp += nWriteSmp;
			if(nWriteSmp != nProcessSmp) break;
		}
	}

	//! Advance and return number of sample points written
	WavState->SamplePosition += nTotalWriteSmp;
	return nTotalWriteSmp / SmpSize;
}

/**************************************/
//! EOF
/**************************************/
