/**************************************/
#pragma once
/**************************************/
#include <cstddef>
#include <cstdint>
#include <cstdio>
/**************************************/

namespace ULC {
	static constexpr size_t MAX_CHANS  = 2;
	static constexpr size_t MAX_QUANTS = 8;
	static constexpr size_t BLOCK_SIZE = 2048;

	//! Nybble-writing callback
	//! NOTE: x is already masked to 0h..Fh
	typedef void (*NybWriteCb_t)(int x, void *User);

	//! Nybble-reading callback
	//! NOTE: Nybbles are not expected to be masked
	typedef int (*NybReadCb_t)(void *User);

	//! Encoder class
	class Encoder_c {
		size_t RateHz;
		size_t nChan;

		double BitBudget;
		double CoefBitRate;
	public:
		//! Returns block size in bits
		//! Data is arranged sequentially for each channel
		//! eg. {L0,L1,L2,L3...Ln},{R0,R1,R2,R3...Rn}, etc
		size_t EncodeBlock(const float *SrcData, size_t RateKbps, NybWriteCb_t NybWriteCbFnc, void *NybWriteCbUsr = nullptr);

		Encoder_c(size_t _RateHz, size_t _nChan) {
			RateHz      = _RateHz;
			nChan       = _nChan;
			BitBudget   = 0.0;
			CoefBitRate = 1.0 / 7.5; //! Initial approximation (reciprocal of bits per nZ coefficient)
		}
	};

	//! Decoder class
	class Decoder_c {
		size_t nChan;
	public:
		//! Data is arranged sequentially for each channel
		//! eg. {L0,L1,L2,L3...Ln},{R0,R1,R2,R3...Rn}, etc
		void DecodeBlock(float *DstData, NybReadCb_t NybReadCbFnc, void *NybReadCbUsr = nullptr);

		Decoder_c(size_t _nChan) : nChan(_nChan) {}
	};
}

/**************************************/
//! EOF
/**************************************/
