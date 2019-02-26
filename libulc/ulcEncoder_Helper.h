/**************************************/
#pragma once
/**************************************/
#include <stddef.h>
#include <stdint.h>
/**************************************/

//! Get quantizer band from band index
static size_t GetQuantBand(size_t Band, const uint16_t *QuantBw) {
	size_t QBand;
	for(QBand=0;;QBand++) {
		size_t Bw = QuantBw[QBand];
		if(Band < Bw) return QBand;
		Band -= Bw;
	}
}

/**************************************/
//! EOF
/**************************************/
