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
#define ABS(x) ((x) < 0 ? (-(x)) : (x))
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
//! NOTE: Biased by -2 because 2^2=4 which is the 'middle' value between 1..7
//!       for quantized coefficients
static int16_t Block_Encode_BuildQuantizer(double Pow, double Abs) {
	if(Abs == 0.0) return 0;

	long int sd = lrint(log(Pow / Abs) * 0x1.71547652B82FEp0) - 2; //! Log2[x^(1/m)]=Log[x]/Log[2^m]
	if(sd <  0) sd =  0;
	if(sd > 14) sd = 14; //! Fh is reserved
	return 1 << sd;
}

//! Build quantizer bands
static void Block_Encode_BuildQBands(const float *Coefs, uint16_t *QuantsBw, size_t BlockSize) {
	size_t Band;
	size_t BandsRem = BlockSize;
	size_t nQBands = 0;
	size_t QBandBw = 0, QBandNzBw = 0;
	double SumSqr  = 0.0;
	for(Band=0;Band<BlockSize;Band++) {
		//! Codeable?
		double vNew = SQR((double)Coefs[Band]);
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
					QuantsBw[nQBands++] = QBandBw;
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
	QuantsBw[nQBands++] = BandsRem;
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
	float  Val;
#define FETCH_KEY_DATA(KeyIdx) \
	Band  = Keys[KeyIdx].Key & (BlockSize-1),  \
	Chan  = Keys[KeyIdx].Key >> BlockSizeLog2, \
	QBand = Block_Encode_BuildQuants_GetQBand(Band, QuantsBw[Chan]), \
	Val   = ABS(CoefBuffer[Chan][Band])

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

	//! Clear quantizer state
	for(Chan=0;Chan<nChan;Chan++) {
		for(QBand=0;QBand<MAX_QUANTS;QBand++) {
			QuantsPow[Chan][QBand] = 0.0;
			QuantsAbs[Chan][QBand] = 0.0;
			//Quants   [Chan][QBand] = 0; //! <- Will be set during first pass below
		}
	}

	//! Build quantizer bands
	for(Chan=0;Chan<nChan;Chan++) Block_Encode_BuildQBands(CoefBuffer[Chan], QuantsBw[Chan], BlockSize);

	//! Clip to maximum available keys
	//! This can happen if rate control says we can fit more
	//! non-zero bands than are actually present for this block
	if(nNzMax > nKeys) nNzMax = nKeys;

	//! Build initial quantizers by considering
	//! all available keys (ie. the 'first pass')
	for(Key=0;Key<nNzMax;Key++) {
		FETCH_KEY_DATA(Key);

		QuantsPow[Chan][QBand] += SQR((double)Val);
		QuantsAbs[Chan][QBand] += Val;
	}
	for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<MAX_QUANTS;QBand++) {
		Quants[Chan][QBand] = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsAbs[Chan][QBand]);
	}

	//! Loop until quantizers stabilize
	//! NOTE: Saving quantizers to key value; this avoids
	//! a lookup inside the coding loop
	size_t nNzBands;
	for(nNzBands = 0;;) {
		//! Remove any collapsed keys
		size_t nNzMaxOrig = nNzMax;
		while(nNzBands < nNzMaxOrig) {
			FETCH_KEY_DATA(nNzBands);

			//! Collapse?
			int16_t Qnt = Quants[Chan][QBand];
			if(Val < 0.5f*Qnt) {
				//! Remove key from analysis and rebuild quantizer
				QuantsPow[Chan][QBand] -= SQR((double)Val);
				QuantsAbs[Chan][QBand] -= Val;
				Quants   [Chan][QBand] = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsAbs[Chan][QBand]);

				//! Remove key from list
				Analysis_KeyRemove(nNzBands, Keys, nKeys); nKeys--;
				nNzMaxOrig--;
			} else Keys[nNzBands++].Quant = Qnt;
		}

		//! Clip nNzMax again: nKeys might've dropped enough that nNzMax needs to change again
		if(nNzMax > nKeys) nNzMax = nKeys;

		//! Try to add more keys if any collapsed above
		//! NOTE: Normally, no keys collapse above, so
		//! this section isn't used and we break out of
		//! the rate/quantizer optimization loop during
		//! the first iteration; this is 'just in case'
		size_t QuantsUpdated = 0;
		while(nNzBands < nNzMax) {
			FETCH_KEY_DATA(nNzBands);

			//! Does this key collapse if we try to fit it to the analysis?
			double  Pow = QuantsPow[Chan][QBand] + SQR((double)Val);
			double  Abs = QuantsAbs[Chan][QBand] + Val;
			int16_t Qnt = Block_Encode_BuildQuantizer(Pow, Abs);
			if(Val < 0.5f*Qnt) {
				//! Remove key from list
				Analysis_KeyRemove(nNzBands, Keys, nKeys); nKeys--;
				if(nNzMax > nKeys) nNzMax = nKeys;
			} else {
				//! No collapse - save new quantizer state
				QuantsPow[Chan][QBand] = Pow;
				QuantsAbs[Chan][QBand] = Abs;
				if(Quants[Chan][QBand] != Qnt) {
					//! Only count /modified/ quantizers,
					//! not the ones we just created
					if(Quants[Chan][QBand] != 0) QuantsUpdated++;
					Quants[Chan][QBand] = Qnt;
				}
				Keys[nNzBands++].Quant = Qnt;
			}
		}

		//! If quantizers aren't modified, we're done
		if(!QuantsUpdated) break;
	}
	return nNzBands;
#undef FETCH_KEY_DATA
}

/**************************************/
//! EOF
/**************************************/
