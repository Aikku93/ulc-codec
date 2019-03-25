/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#include <math.h>
#include <stddef.h>
#include <stdint.h>
/**************************************/
#include "ulcEncoder.h"
#include "ulcEncoder_Analysis.h"
#include "ulcEncoder_Helper.h"
/**************************************/
#define SQR(x) ((x)*(x))
/**************************************/

//! Maximum allowed quantizers
//! NOTE: NO GREATER THAN 48
#define MAX_QUANTS 48

/**************************************/

//! Get quantizer band from band index
static size_t Block_Encode_BuildQuants_GetQBand(size_t Band, const uint16_t *QuantBw) {
	size_t QBand;
	for(QBand=0;;QBand++) {
		size_t Bw = QuantBw[QBand];
		if(Band < Bw) return QBand;
		Band -= Bw;
	}
}

//! Build quantizer from sum of raised-power values and sum of absolutes
//! NOTE: Currently using Sum[x^2]/Sum[x]. This somewhat favours larger
//!       values which mask lower-power values anyway so it works out better.
static inline int16_t Block_Encode_BuildQuantizer(double Pow, double Abs) {
	if(Abs == 0.0) return 0;
	double sd = round(log(Pow / Abs) * 0x1.71547652B82FEp0 - 2.0); //! Log2[x^(1/m)]=Log[x]/Log[2^m]
	if(sd <  0.0) sd =  0.0;
	if(sd > 14.0) sd = 14.0; //! Fh is reserved
	return 1 << (size_t)sd;
}

/**************************************/

//! Build quantizers
//! Returns the number of non-zero bands kept
//! TODO:
//!  -Optimization: Don't /actually/ remove keys, just skip them (this saves shifting the array)
static size_t Block_Encode_BuildQuants(const struct ULC_EncoderState_t *State, size_t nNzMax, size_t nKeys) {
	size_t Key;
	size_t BlockSize     = State->BlockSize;
	size_t BlockSizeLog2 = IntLog2(BlockSize);

	//! Key data fetcher
	size_t Band, Chan, QBand;
	double Val;
#define FETCH_KEY_DATA(KeyIdx) \
	Band  = Keys[KeyIdx].Key & (BlockSize-1),  \
	Chan  = Keys[KeyIdx].Key >> BlockSizeLog2, \
	QBand = Block_Encode_BuildQuants_GetQBand(Band, QuantsBw[Chan]), \
	Val   = CoefBuffer[Chan][Band]

	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	size_t nChan          = State->nChan;
	float    **CoefBuffer = State->TransformBuffer;
	uint16_t **QuantsBw   = State->QuantsBw;
	int16_t  **Quants     = State->Quants;
	double   **QuantsPow  = State->QuantsPow;
	double   **QuantsAbs  = State->QuantsAbs;
	struct AnalysisKey_t *Keys = State->AnalysisKeys;

	//! Build quantizer bands
	for(Chan=0;Chan<nChan;Chan++) {
		const float *Coefs = CoefBuffer[Chan];
		size_t BandsRem = BlockSize;
		size_t nQBands = 0;
		size_t QBandBw = 0, QBandNzBw = 0;
		double SumSqr  = 0.0;
		for(Band=0;Band<BlockSize;Band++) {
			//! Codeable?
			double vNew = Coefs[Band]; vNew = SQR(vNew);
			if(vNew >= SQR(0.5)) {
				//! Enough bands to decide on a split?
				//! NOTE: Somewhat arbitrary and less sensitive at high freq
				size_t QBandBwThres = 1 + Band/32;
				if(QBandNzBw > QBandBwThres) {
					//! Coefficient not in range?
					//! NOTE: Somewhat arbitrary (though tuned) thresholds
					double t = vNew*QBandBw;
					if(t < SQR(1.0/8.0)*SumSqr || t > SQR(4.0)*SumSqr) {
						//! Create a split point
						//! NOTE: Last band is built from remaining coefficients
						BandsRem -= QBandBw;
						QuantsBw[Chan][nQBands++] = QBandBw;
						if(nQBands == MAX_QUANTS-1) break;

						//! Reset state for a new band
						QBandBw = QBandNzBw = 0;
						SumSqr  = 0.0;
					}
				}

				//! Add to quantizer band
				SumSqr += vNew;
				QBandNzBw++;
			}

			//! Increase bandwidth
			QBandBw++;
		}
		QuantsBw[Chan][nQBands++] = BandsRem;
	}

	//! Clip to maximum available keys
	//! This can happen if rate control says we can fit more
	//! non-zero bands than are actually present for this block
	if(nNzMax > nKeys) nNzMax = nKeys;

	//! Clear quantizer state
	for(Chan=0;Chan<nChan;Chan++) {
		for(QBand=0;QBand<MAX_QUANTS;QBand++) {
			QuantsPow[Chan][QBand] = 0.0;
			QuantsAbs[Chan][QBand] = 0.0;
		}
	}

	//! Build initial quantizers by considering all the keys available
	for(Key=0;Key<nNzMax;Key++) {
		FETCH_KEY_DATA(Key);
		if(Val < 0.0) Val = -Val;

		QuantsPow[Chan][QBand] += Val*Val;
		QuantsAbs[Chan][QBand] += Val;
	}
	for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<MAX_QUANTS;QBand++) {
		Quants[Chan][QBand] = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsAbs[Chan][QBand]);
	}

	//! Remove any collapsed keys
	size_t nNzBands = 0, nNzBandsOrig = nNzMax;
	while(nNzBands < nNzBandsOrig) {
		FETCH_KEY_DATA(nNzBands);
		if(Val < 0.0) Val = -Val;

		//! Collapse?
		if(Val < 0.5*Quants[Chan][QBand]) {
			//! Remove key from analysis and rebuild quantizer
			QuantsPow[Chan][QBand] -= Val*Val;
			QuantsAbs[Chan][QBand] -= Val;
			Quants   [Chan][QBand] = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsAbs[Chan][QBand]);

			//! Remove key from list
			Analysis_KeyRemove(nNzBands, Keys, nKeys); nKeys--;
			nNzBandsOrig--;
		} else nNzBands++;
	}

	//! Try to add more keys if any collapsed above
	//! Note clipping again: nKeys might've dropped enough that nNzMax needs to change again
	if(nNzMax > nKeys) nNzMax = nKeys;
	while(nNzBands < nNzMax) {
		FETCH_KEY_DATA(nNzBands);
		if(Val < 0.0) Val = -Val;

		//! Does this key collapse if we try to fit it to the analysis?
		double  Pow = QuantsPow[Chan][QBand] + Val*Val;
		double  Avg = QuantsAbs[Chan][QBand] + Val;
		int32_t Qnt = Block_Encode_BuildQuantizer(Pow, Avg);
		if(Val < 0.5*Qnt) {
			//! Remove key from list
			Analysis_KeyRemove(nNzBands, Keys, nKeys); nKeys--;
			if(nNzMax > nKeys) nNzMax = nKeys;
		} else {
			//! No collapse - save new quantizer state
			QuantsPow[Chan][QBand] = Pow;
			QuantsAbs[Chan][QBand] = Avg;
			Quants   [Chan][QBand] = Qnt;
			nNzBands++;
		}
	}

	//! Double-check that no keys we kept collapse
	//! This is needed due to the psy-opt that amplified high frequency
	//! data, since now we might have smaller values stuck in the middle
	//! Additionally save the quantizer directly to the key structure,
	//! as this avoids a lookup inside the coding loop
	//! Finally, also try to keep any extra keys if we remove some, even
	//! if they will quantize sub-optimally
	for(Key=0;Key<nNzBands;) {
		FETCH_KEY_DATA(Key);
		if(Val < 0.0) Val = -Val;

		int16_t q = Quants[Chan][QBand];
		if(Val < 0.5*q) {
			Analysis_KeyRemove(Key, Keys, nKeys); nKeys--;
			if(nNzBands > nKeys) nNzBands = nKeys;
		} else Keys[Key++].Quant = q;
	}

	return nNzBands;
#undef FETCH_KEY_DATA
}

/**************************************/
//! EOF
/**************************************/
