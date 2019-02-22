/**************************************/
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
/**************************************/
#include "Fourier.h"
#include "ulcCoder.h"
#include "ulcUtility.hpp"
/**************************************/
using namespace std;
using namespace ULC;
/**************************************/

//! Quantizer bandwidths and bin ranges
//! TODO: Avoid hardcoding this?
static const size_t QuantBinBw[MAX_QUANTS] = {16,32,48,96,192,384,768,512};
static const size_t QuantBinRng[MAX_QUANTS][2] = {{0,15},{16,47},{48,95},{96,191},{192,383},{384,767},{768,1535},{1536,2047}};

//! Transform buffers
__attribute__((aligned(32))) static float TransformBuffer[BLOCK_SIZE];
__attribute__((aligned(32))) static float TransformTemp  [BLOCK_SIZE];
__attribute__((aligned(32))) static float TransformInvLap[MAX_CHANS][BLOCK_SIZE/2] = {{0}};

/**************************************/

//! Get quantizer band from band index
static size_t GetQuantBand(size_t Band) {
	//! Perform binary search over bin ranges
	//! PONDER: Might be able to speed this up inside the
	//!         loops that use this, with linear search
	int Lo = 0, Hi = MAX_QUANTS-1;
	while(Lo <= Hi) {
		size_t Mid = (Lo + Hi) / 2;
		     if(Band < QuantBinRng[Mid][0]) Hi = Mid-1;
		else if(Band > QuantBinRng[Mid][1]) Lo = Mid+1;
		else return Mid;
	}
	return ~0U; //! Not used
}

/**************************************/

void Decoder_c::DecodeBlock(float *DstData, NybReadCb_t NybReadCbFnc, void *NybReadCbUsr) {
	for(size_t Chan=0;Chan<nChan;Chan++) {
		//! Read quantizers
		int32_t Quants[MAX_QUANTS];
		for(size_t i=0;i<MAX_QUANTS;i++) {
			Quants[i] = 1 << (NybReadCbFnc(NybReadCbUsr) & 0xF);
		}

		//! Start decoding coefficients
		size_t NextQuantBin = 0;
		for(;;) {
			//! Skip unused quantizer bands
			while(NextQuantBin < MAX_QUANTS && Quants[NextQuantBin] == 1<<15) {
				size_t Lo = QuantBinRng[NextQuantBin][0];
				size_t Hi = QuantBinRng[NextQuantBin][1];
				for(size_t i=Lo;i<=Hi;i++) TransformBuffer[i] = 0.0f;
				NextQuantBin++;
			}
			if(NextQuantBin >= MAX_QUANTS) break;

			//! Find range of bands we can read
			size_t NextNz = QuantBinRng[NextQuantBin][0];
			size_t LastNz;
			do LastNz = QuantBinRng[NextQuantBin][1]; while(++NextQuantBin < MAX_QUANTS && Quants[NextQuantBin] != 1<<15);

			//! Read the coefficients
			//! NOTE: Odd number of zeros run are coded as long run followed by coded zeros
			//!       ie. always {Run(n),0h}, never {0h,Run(n)} or {0h,Run(n),0h}
			int32_t LastV = 1; //! Anything non-zero
			while(NextNz <= LastNz) {
				//! Normal coefficient?
				int32_t v = (NybReadCbFnc(NybReadCbUsr) << 28) >> 28;
				if(v != -0x8) {
					//! Store dequantized
					TransformBuffer[NextNz] = v * Quants[GetQuantBand(NextNz)];
					if(++NextNz > LastNz) break;

					//! If we have two zeros in a row, this is the stop code 0h,0h
					if((v|LastV) == 0) {
						do TransformBuffer[NextNz] = 0.0f; while(++NextNz <= LastNz);
						break;
					}

					LastV = v;
					continue;
				}

				//! Unpack zero run
				LastV = v;
				size_t nZ = NybReadCbFnc(NybReadCbUsr) & 0xF;
				if(nZ < 0xC) {
					//! Small run
					//! 0h..Bh: 2..24 zeros
					nZ = nZ*2 + 2;
				} else {
					//! Long run
					//! Ch..Fh,Xh: 26..152 zeros (Ch + n>>4, n&Fh)
					nZ = (nZ-0xC)<<4 | (NybReadCbFnc(NybReadCbUsr) & 0xF);
					nZ = nZ*2 + 26;
				}

				//! Insert zeros
				do TransformBuffer[NextNz++] = 0.0f; while(--nZ);
			}
		}

		//! Inverse transform block
		Transform_SumDif(TransformBuffer, BLOCK_SIZE);
		Fourier::IMDCT(DstData + Chan*BLOCK_SIZE, TransformBuffer, TransformInvLap[Chan], TransformTemp, BLOCK_SIZE);
	}
}

/**************************************/
//! EOF
/**************************************/
