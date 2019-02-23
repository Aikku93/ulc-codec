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

//! Transform buffers
__attribute__((aligned(32))) static float  TransformBuffer[MAX_CHANS][BLOCK_SIZE];
__attribute__((aligned(32))) static float  TransformFwdLap[MAX_CHANS][BLOCK_SIZE/2] = {{0}};
__attribute__((aligned(32))) static float  TransformTemp  [MAX_CHANS][BLOCK_SIZE];
__attribute__((aligned(32))) static size_t AnalysisKeys   [MAX_CHANS*BLOCK_SIZE];

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

//! Encode analysis key data
static inline size_t Analysis_KeyEncode(size_t Band, size_t Chan) {
	return Band | Chan<<24;
}

//! Decode analysis key data
static inline void Analysis_KeyDecode(size_t &Band, size_t &Chan, size_t Key) {
	Chan = Key >> 24;
	Band = Key &~ (Chan<<24);
}

//! Insert key to analysis buffer
//! Sorted by highest power (passed in 'v') first
//! NOTE: Analysis values are stored in TransformTemp[][]
static void Analysis_KeyInsert(float v, size_t &nKeys, size_t Band, size_t Chan) {
	//! Get insertion index (binary search for yeet points)
	size_t InsertIdx; {
		int Lo = 0, Hi = nKeys-1;
		while(Lo <= Hi) {
			size_t Mid = (Lo + Hi) / 2;
			size_t MidBand, MidChan; Analysis_KeyDecode(MidBand, MidChan, AnalysisKeys[Mid]);
			float  vMid = TransformTemp[MidChan][MidBand];
			if(v > vMid) Hi = Mid-1;
			else         Lo = Mid+1;
		}
		InsertIdx = Lo;
	}

	//! Shift and insert
	//! PONDER: Faster way of doing the shifting?
	for(size_t i=nKeys;i>InsertIdx;i--) AnalysisKeys[i] = AnalysisKeys[i-1];
	AnalysisKeys[InsertIdx] = Analysis_KeyEncode(Band, Chan);
	TransformTemp[Chan][Band] = v;
	nKeys++;
}

//! Remove key from analysis buffer
static void Analysis_KeyRemove(size_t Key, size_t &nKeys) {
	//! PONDER: Faster way of doing the shifting?
	for(size_t i=Key+1;i<nKeys;i++) AnalysisKeys[i-1] = AnalysisKeys[i];
	nKeys--;
}

/**************************************/

//! Apply block transform
//!  -Fetches data
//!  -Applies MDCT
//!  -Applies anti-pre-echo formula
//!  -Stores keys for block coefficients
//! Returns the number of keys stored
static size_t Block_Transform(const float *Data, size_t nChan) {
	//! Do transforms
	for(size_t Chan=0;Chan<nChan;Chan++) {
		//! Fetch sample data
		//! Pre-scale for scaled IMDCT(*2.0) and SumDif transform(*0.5)
		constexpr float CoefScale = 2.0f*0.5f / BLOCK_SIZE;
		const float *Src = Data + Chan*BLOCK_SIZE;
#if defined(__AVX__)
		__m256 Scale = _mm256_set1_ps(CoefScale);
		for(size_t i=0;i<BLOCK_SIZE;i+=8) {
			__m256 v = _mm256_load_ps(Src + i);
			v = _mm256_mul_ps(v, Scale);
			_mm256_store_ps(TransformTemp[0] + i, v);
		}
#elif defined(__SSE__)
		__m128 Scale = _mm_set1_ps(CoefScale);
		for(size_t i=0;i<BLOCK_SIZE;i+=4) {
			__m128 v = _mm_load_ps(Src + i);
			v = _mm_mul_ps(v, Scale);
			_mm_store_ps(TransformTemp[0] + i, v);
		}
#else
		for(size_t i=0;i<BLOCK_SIZE;i++) TransformTemp[0][i] = Src[i] * CoefScale;
#endif
		//! Apply transforms
		Fourier::MDCT(TransformBuffer[Chan], TransformTemp[0], TransformFwdLap[Chan], TransformTemp[1], BLOCK_SIZE);
		Transform_SumDif(TransformBuffer[Chan], BLOCK_SIZE);
	}

	//! Build keys
	size_t nKeys = 0;
	float ChanPow = 1.0f;
	for(size_t Chan=0;Chan<nChan;Chan++,ChanPow *= float(M_SQRT1_2)) {
		for(size_t i=0;i<BLOCK_SIZE;i++) {
			//! Check that value doesn't collapse to 0 under the smallest quantizer (1.0)
			//! Check notes below for explanation of comparison against Quant*0.5
			float v = pow(TransformBuffer[Chan][i], 2.0f);
			if(v < pow(1.0f*0.5f, 2.0f)) continue;

			//! Do psy-opts and insert:
			//!  -Decrease side-chan importance
			//!  -Increase high-freq importance
			//! TODO: Accelerator table
			v *= ChanPow;
			v *= 1.0f - 0.5f*pow(1.0f - (i+0.5f)/BLOCK_SIZE, 2.0f);
			Analysis_KeyInsert(v, nKeys, i, Chan);
		}
	}
	return nKeys;
}

//! Build RMS-fit quantizer
//! NOTE:
//!  -Division by 4.0 to set that as the 'mid' quantized value
//!  -Fh = Unused quantizer band
//! PONDER: Shouldn't need to clip to <=14 here under normal circumstances:
//!  round(log2((32768*Pi/4) / 4)) = 13
//! 32768*Pi/4 comes from:
//!  32768: Signed 16bit data range
//!  Pi/4:  Infinity norm of DCT-IV (as N approaches infinity)
static int32_t BuildQuantizer(double BandPow, size_t n) {
	size_t s = 0xF;
	if(n) s = max(0.0, min(14.0, round(log2(sqrt(BandPow / n) / 4.0))));
	return 1 << s;
}

//! Build quantizers for quantizer bands
//! Returns the number of non-zero bands kept
static size_t Block_GetQuants(int32_t Quants[MAX_CHANS][MAX_QUANTS], size_t nNzBandsMax, size_t nKeys, size_t nChan) {
	//! LAMBDA EXPRESSION
	//! Fetch data associated with key
	float v;
	size_t Band, Chan, QBand;
	auto KeyDataFetch = [&](size_t Key) {
		Analysis_KeyDecode(Band, Chan, AnalysisKeys[Key]);
		QBand = GetQuantBand(Band);
		v = TransformBuffer[Chan][Band];
	};

	//! Clip to maximum available keys
	//! This can happen if rate control says we can fit more
	//! non-zero bands than are actually present for this block
	nNzBandsMax = min(nNzBandsMax, nKeys);

	//! Put all initial keys into analysis
	double BandPow[MAX_CHANS][MAX_QUANTS] = {{0.0}};
	size_t BandCnt[MAX_CHANS][MAX_QUANTS] = {{0}};
	for(size_t Key=0;Key<nNzBandsMax;Key++) {
		KeyDataFetch(Key);

		BandPow[Chan][QBand] += pow(v, 2.0);
		BandCnt[Chan][QBand]++;
	}
	for(size_t Chan=0;Chan<nChan;Chan++) for(size_t QBand=0;QBand<MAX_QUANTS;QBand++) {
		Quants[Chan][QBand] = BuildQuantizer(BandPow[Chan][QBand], BandCnt[Chan][QBand]);
	}

	//! Get number of [analyzed] keys we kept, and remove any collapsed keys
	size_t nNzBands = 0, nNzBandsOrig = nNzBandsMax;
	while(nNzBands < nNzBandsOrig) {
		KeyDataFetch(nNzBands);

		//! Collapse?
		//! Compare with 0.5*Quant, because:
		//!  round(v / Quant) == 0 <-> abs(v) < 0.5*Quant
		//! eg. Quant == 1:
		//!  v == 0.5: round(v / Quant) -> 1 (not collapsed)
		//!  v == 0.4: round(v / Quant) -> 0 (collapsed)
		int32_t &Quant = Quants[Chan][QBand];
		if(abs(v) < 0.5f*Quant) {
			//! Remove key from analysis and rebuild quantizer
			BandPow[Chan][QBand] -= pow(v, 2.0);
			BandCnt[Chan][QBand]--;
			Quant = BuildQuantizer(BandPow[Chan][QBand], BandCnt[Chan][QBand]);

			//! Remove key from list
			Analysis_KeyRemove(nNzBands, nKeys);
			nNzBandsOrig--;
		} else nNzBands++;
	}

	//! Try to add more keys if any collapsed above
	//! Note clipping again: nKeys might've dropped enough that nNzBandsMax needs to change again
	nNzBandsMax = min(nNzBandsMax, nKeys);
	while(nNzBands < nNzBandsMax) {
		KeyDataFetch(nNzBands);

		//! Does this key collapse if we try to fit it to the analysis?
		double  Pow = BandPow[Chan][QBand] + pow(v, 2.0);
		size_t  Cnt = BandCnt[Chan][QBand] + 1;
		int32_t q   = BuildQuantizer(Pow, Cnt);
		if(abs(v) < 0.5f*q) {
			//! Remove key from list
			Analysis_KeyRemove(nNzBands, nKeys);
			nNzBandsMax = min(nNzBandsMax, nKeys);
		} else {
			//! No collapse - save new quantizer state
			BandPow[Chan][QBand] = Pow;
			BandCnt[Chan][QBand] = Cnt;
			Quants [Chan][QBand] = q;
			nNzBands++;
		}
	}

	//! Double-check that no keys we kept collapse
	//! This is needed due to the psy-opt that amplified high frequency
	//! data, since now we might have smaller values stuck in the middle
	for(size_t Key=0;Key<nNzBands;) {
		KeyDataFetch(Key);

		if(abs(v) < 0.5f*Quants[Chan][QBand]) {
			Analysis_KeyRemove(Key, nKeys);
			nNzBands--;
		} else Key++;
	}

	//! If a quantizer has no keys, it needs to be cleared
	for(size_t Chan=0;Chan<nChan;Chan++) for(size_t QBand=0;QBand<MAX_QUANTS;QBand++) {
		if(!BandCnt[Chan][QBand]) Quants[Chan][QBand] = 1<<15;
	}
	return nNzBands;
}

//! Format:
//! -Pass quantizers (Fh = unused)
//! -Pass quantized values where quantizer != Fh
//! Returns the block size (in bits) and the number of coded (non-zero) coefficients
static size_t Block_Encode(size_t nNzBandsMax, size_t nKeys, size_t nChan, size_t *_nCoded, NybWriteCb_t NybWriteCbFnc, void *NybWriteCbUsr) {
	//! Get quantizers
	int32_t Quants[MAX_CHANS][MAX_QUANTS];
	size_t nNzBands = Block_GetQuants(Quants, nNzBandsMax, nKeys, nChan);

	//! Sort keys by band index
	//! This avoids a search for the next non-zero band
	//! Also because the channel is coded in the high bits, we can
	//! code one channel at a time, making things easier as well
	sort(AnalysisKeys, AnalysisKeys + nNzBands);

	//! Start coding
	size_t Size    = 0; //! Block size (in bits)
	size_t nCoded  = 0; //! Coded (non-zero) coefficients
	size_t NextKey = 0;
	for(size_t Chan=0;Chan<nChan;Chan++) {
		//! Code the quantizer values (in log2 form)
		for(size_t i=0;i<MAX_QUANTS;i++) {
			size_t s = 31 - __builtin_clz(Quants[Chan][i]);
			NybWriteCbFnc(s, NybWriteCbUsr); Size += 4;
		}

		//! Start coding coefficients
		size_t NextQuantBin = 0;
		for(;;) {
			//! Skip unused quantizer bands
			while(NextQuantBin < MAX_QUANTS && Quants[Chan][NextQuantBin] == 1<<15) NextQuantBin++;
			if(NextQuantBin >= MAX_QUANTS) break;

			//! Find range of bands we can store
			size_t NextNz = QuantBinRng[NextQuantBin][0];
			size_t LastNz;
			do LastNz = QuantBinRng[NextQuantBin][1]; while(++NextQuantBin < MAX_QUANTS && Quants[Chan][NextQuantBin] != 1<<15);

			//! Code the coefficients
			for(;;) {
				//! Unpack key data
				size_t Band; {
					size_t tChan;
					Analysis_KeyDecode(Band, tChan, AnalysisKeys[NextKey]);

					//! If we cross to the next quantizer band or
					//! next channel, stop coding.
					//! However, if the quantizer band was used, then
					//! we're guaranteed at least one coefficient here
					if(Band > LastNz) break;
					if(tChan != Chan) break;
				}

				//! Insert zeros run[s]
				//! NOTE: Always in multiples of two
				size_t zR = Band - NextNz;
				while(zR >= 2) {
					//! 8h
					NybWriteCbFnc(0x8, NybWriteCbUsr); Size += 4;

					//! Small run?
					size_t n;
					if(zR < 26) {
						//! 0h..Bh: 2..24 zeros
						size_t v = (zR-2) / 2;
						NybWriteCbFnc(v, NybWriteCbUsr); Size += 4;
						n = v*2 + 2;
					} else {
						//! Ch..Fh,Xh: 26..152 zeros (Ch + n>>4, n&Fh)
						size_t v = min(0x3Fu, (zR-26) / 2);
						NybWriteCbFnc(0xC + (v>>4), NybWriteCbUsr); Size += 4;
						NybWriteCbFnc(v&0xF,        NybWriteCbUsr); Size += 4;
						n = v*2 + 26;
					}

					//! Insert zeros
					zR -= n;
					NextNz += n;
				}

				//! Insert coded coefficients
				//! NOTE:
				//!  zR might still have one more coefficient marked for skipping
				//!  but this didn't take into account the actual statistics of
				//!  the coded zero runs. This means that the coefficient might
				//!  actually not collapse to 0, so we may as well code it anyway
				//!  as it would cost the same either way (though it might quantize
				//!  sub-optimally from not being considered originally)
				while(NextNz <= Band) {
					size_t  Band  = NextNz++;
					size_t  QBand = GetQuantBand(Band);
					int32_t Quant = Quants[Chan][QBand];

					//! -7h..+7h
					int32_t Qn = int32_t(max(-7.0f, min(+7.0f, round(TransformBuffer[Chan][Band] / Quant))));
					NybWriteCbFnc(Qn&0xF, NybWriteCbUsr); Size += 4;
					if(Qn != 0) nCoded++;
				}

				//! Coded all coefficients?
				if(++NextKey >= nNzBands) break;
			}

			//! Finalize the block (0h,0h: Stop)
			size_t n = min(2u, LastNz+1 - NextNz);
			for(size_t i=0;i<n;i++) {
				NybWriteCbFnc(0x0, NybWriteCbUsr); Size += 4;
			}
		}
	}

	//! Return output size
	if(_nCoded) *_nCoded = nCoded;
	return Size;
}

/**************************************/

size_t Encoder_c::EncodeBlock(const float *SrcData, size_t RateKbps, NybWriteCb_t NybWriteCbFnc, void *NybWriteCbUsr) {
	//! Refill bit budget
	double AvgBitBudget = RateKbps*1000.0/RateHz * BLOCK_SIZE;
	BitBudget += AvgBitBudget;

	//! Transform input, build keys, and get maximum number of non-zero bands
	//! NOTE: Limit to 1.25x target bitrate to avoid quality spikes
	//!       as these have a tendency to sound bad in context
	size_t nKeys = Block_Transform(SrcData, nChan);
	size_t nNzBandsMax = min(nKeys, size_t(max(0.0, round(min(BitBudget, AvgBitBudget*1.25) * CoefBitRate))));

	//! Encode block
	size_t nCoded;
	size_t BlockBits = Block_Encode(nNzBandsMax, nKeys, nChan, &nCoded, NybWriteCbFnc, NybWriteCbUsr);

	//! Update state
	double RateDecayFactor = (nCoded == 0) ? 0.0 : min(1.0, 0.25 + nCoded*0.75/nNzBandsMax); //! Bias slightly towards previous block
	double BlockCoefCost = double(nCoded) / BlockBits;
	BitBudget  -= BlockBits;
	CoefBitRate = CoefBitRate*(1.0-RateDecayFactor) + BlockCoefCost*RateDecayFactor;

	//! Return block size
	return BlockBits;
}

/**************************************/
//! EOF
/**************************************/
