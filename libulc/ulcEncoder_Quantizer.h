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

//! Quantizers are coded in log2 form, and Fh is reserved for 'unused quantizer band'
#define QUANTIZER_UNUSED (-(1 << 0xF))

/**************************************/

//! Build quantizer from RMS of coefficients
static inline int16_t Block_Encode_BuildQuantizer(double QuantsPow, double QuantsAvg) {
	if(QuantsAvg == 0.0) return QUANTIZER_UNUSED;
	double sd = round(log(QuantsPow / QuantsAvg) * 0x1.71547652B82FEp0 - 2.0); //! Log2[x]=Log[x]/Log[2]
	if(sd <  0.0) sd =  0.0;
	if(sd > 14.0) sd = 14.0; //! Fh is reserved
	return 1 << (size_t)sd;
}

/**************************************/

//! Build quantizers for quantizer bands
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
	QBand = GetQuantBand(Band, QuantsBw),   \
	Val   = CoefBuffer[Chan][Band]

	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	size_t nChan             = State->nChan;
	size_t nQuants           = State->nQuants;
	float    **CoefBuffer    = State->TransformBuffer;
	const uint16_t *QuantsBw = State->QuantsBw;
	int16_t  **Quants        = State->Quants;
	double   **QuantsPow     = State->QuantsPow;
	double   **QuantsAvg     = State->QuantsAvg;
	struct AnalysisKey_t *Keys = State->AnalysisKeys;

	//! Clip to maximum available keys
	//! This can happen if rate control says we can fit more
	//! non-zero bands than are actually present for this block
	if(nNzMax > nKeys) nNzMax = nKeys;

	//! Clear quantizer state
	for(Chan=0;Chan<nChan;Chan++) {
		for(QBand=0;QBand<State->nQuants;QBand++) {
			QuantsPow[Chan][QBand] = 0.0;
			QuantsAvg[Chan][QBand] = 0.0;
		}
	}

	//! Build initial quantizers by considering all the keys available
	for(Key=0;Key<nNzMax;Key++) {
		FETCH_KEY_DATA(Key);
		if(Val < 0.0) Val = -Val;

		QuantsPow[Chan][QBand] += Val*Val;
		QuantsAvg[Chan][QBand] += Val;
	}
	for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<nQuants;QBand++) {
		Quants[Chan][QBand] = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsAvg[Chan][QBand]);
	}

	//! Remove any collapsed keys
	size_t nNzBands = 0, nNzBandsOrig = nNzMax;
	while(nNzBands < nNzBandsOrig) {
		FETCH_KEY_DATA(nNzBands);
		if(Val < 0.0) Val = -Val;

		//! Collapse?
		//! Compare with 0.5*Quant, because:
		//!  round(v / Quant) == 0 <-> abs(v) < 0.5*Quant
		//! eg. Quant == 1:
		//!  v == 0.5: round(v / Quant) -> 1 (not collapsed)
		//!  v == 0.4: round(v / Quant) -> 0 (collapsed)
		if(Val < 0.5*Quants[Chan][QBand]) {
			//! Remove key from analysis and rebuild quantizer
			QuantsPow[Chan][QBand] -= Val*Val;
			QuantsAvg[Chan][QBand] -= Val;
			Quants   [Chan][QBand] = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsAvg[Chan][QBand]);

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
		double  Avg = QuantsAvg[Chan][QBand] + Val;
		int32_t Qnt = Block_Encode_BuildQuantizer(Pow, Avg);
		if(Val < 0.5*Qnt) {
			//! Remove key from list
			Analysis_KeyRemove(nNzBands, Keys, nKeys); nKeys--;
			if(nNzMax > nKeys) nNzMax = nKeys;
		} else {
			//! No collapse - save new quantizer state
			QuantsPow[Chan][QBand] = Pow;
			QuantsAvg[Chan][QBand] = Avg;
			Quants   [Chan][QBand] = Qnt;
			nNzBands++;
		}
	}

	//! Double-check that no keys we kept collapse
	//! This is needed due to the psy-opt that amplified high frequency
	//! data, since now we might have smaller values stuck in the middle
	for(Key=0;Key<nNzBands;) {
		FETCH_KEY_DATA(Key);
		if(Val < 0.0) Val = -Val;

		if(Val < 0.5*Quants[Chan][QBand]) {
			Analysis_KeyRemove(Key, Keys, nKeys); nKeys--;
			nNzBands--;
		} else Key++;
	}
	return nNzBands;
#undef FETCH_KEY_DATA
}

/**************************************/
//! EOF
/**************************************/
