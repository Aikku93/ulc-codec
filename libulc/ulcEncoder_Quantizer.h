/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
#include <stdint.h>
/**************************************/
#include "ulcEncoder.h"
#include "ulcEncoder_Analysis.h"
/**************************************/

//! Quantizers are coded in log2 form, and Fh is reserved for 'unused quantizer band'
#define QUANTIZER_UNUSED (-(1 << 0xF))

/**************************************/

//! Build quantizer from RMS of coefficients
static inline int16_t Block_Encode_BuildQuantizer(double QuantsPow, size_t QuantsCnt) {
	if(!QuantsCnt) return QUANTIZER_UNUSED;
#if 0
	double sd = sqrt(QuantsPow / QuantsCnt); //! RMS
	       sd = log2(sd / 4.0);              //! Log2 form (division by 4.0 to set that as the 'middle' value for the quantized range 0..7)
#else //! Possibly faster
	double sd = log(QuantsPow / QuantsCnt)*0x1.71547652B82FEp-1 - 2.0; //! Log[QuantsPow/QuantsCnt]/Log[4] - 2
#endif
	       sd = round(sd);                   //! Rounding
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
	size_t Key, Band, Chan, QBand;
	double Val;
#define FETCH_KEY_DATA(Key) \
	Band  = Keys[Key].Band, \
	Chan  = Keys[Key].Chan, \
	QBand = GetQuantBand(Band, QuantsBw), \
	Val   = CoefBuffer[Chan][Band]

	//! Spill state to local variables to make things easier to read
	//! PONDER: Hopefully the compiler realizes that State is const and
	//!         doesn't just copy the whole thing out to the stack :/
	size_t nChan             = State->nChan;
	size_t nQuants           = State->nQuants;
	float    **CoefBuffer    = State->TransformBuffer;
	const uint16_t *QuantsBw = State->QuantsBw;
	int16_t  **Quants        = State->Quants;
	uint16_t **QuantsCnt     = State->QuantsCnt;
	double   **QuantsPow     = State->QuantsPow;
	struct AnalysisKey_t *Keys = State->AnalysisKeys;

	//! Clip to maximum available keys
	//! This can happen if rate control says we can fit more
	//! non-zero bands than are actually present for this block
	if(nNzMax > nKeys) nNzMax = nKeys;

	//! Clear quantizer state
	for(Chan=0;Chan<nChan;Chan++) {
		for(QBand=0;QBand<State->nQuants;QBand++) {
			QuantsCnt[Chan][QBand] = 0;
			QuantsPow[Chan][QBand] = 0.0;
		}
	}

	//! Build initial quantizers by considering all the keys available
	for(Key=0;Key<nNzMax;Key++) {
		FETCH_KEY_DATA(Key);
		QuantsCnt[Chan][QBand]++;
		QuantsPow[Chan][QBand] += Val*Val;
	}
	for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<nQuants;QBand++) {
		Quants[Chan][QBand] = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsCnt[Chan][QBand]);
	}

	//! Remove any collapsed keys
	size_t nNzBands = 0, nNzBandsOrig = nNzMax;
	while(nNzBands < nNzBandsOrig) {
		FETCH_KEY_DATA(nNzBands);

		//! Collapse?
		//! Compare with 0.5*Quant, because:
		//!  round(v / Quant) == 0 <-> abs(v) < 0.5*Quant
		//! eg. Quant == 1:
		//!  v == 0.5: round(v / Quant) -> 1 (not collapsed)
		//!  v == 0.4: round(v / Quant) -> 0 (collapsed)
		if(Val < 0.0) Val = -Val;
		if(Val < 0.5*Quants[Chan][QBand]) {
			//! Remove key from analysis and rebuild quantizer
			QuantsCnt[Chan][QBand]--;
			QuantsPow[Chan][QBand] -= Val*Val;
			Quants   [Chan][QBand] = Block_Encode_BuildQuantizer(QuantsPow[Chan][QBand], QuantsCnt[Chan][QBand]);

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

		//! Does this key collapse if we try to fit it to the analysis?
		size_t  Cnt = QuantsCnt[Chan][QBand] + 1;
		double  Pow = QuantsPow[Chan][QBand] + Val*Val;
		int32_t Qnt = Block_Encode_BuildQuantizer(Pow, Cnt);
		if(Val < 0.0) Val = -Val;
		if(Val < 0.5*Qnt) {
			//! Remove key from list
			Analysis_KeyRemove(nNzBands, Keys, nKeys); nKeys--;
			if(nNzMax > nKeys) nNzMax = nKeys;
		} else {
			//! No collapse - save new quantizer state
			QuantsCnt[Chan][QBand] = Cnt;
			QuantsPow[Chan][QBand] = Pow;
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

	//! If a quantizer has no keys, it needs to be cleared
	for(Chan=0;Chan<nChan;Chan++) for(QBand=0;QBand<nQuants;QBand++) {
		if(!QuantsCnt[Chan][QBand]) Quants[Chan][QBand] = QUANTIZER_UNUSED;
	}
	return nNzBands;
#undef FETCH_KEY_DATA
}

/**************************************/
//! EOF
/**************************************/