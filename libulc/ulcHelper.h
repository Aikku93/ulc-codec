/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2021, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#pragma once
/**************************************/
#define ABS(x) ((x) < 0 ? (-(x)) : (x))
#define SQR(x) ((x)*(x))
/**************************************/

//! Subblock count for a block
static inline int ULC_Helper_SubBlockCount(int Decimation) {
	static const int8_t Count[8] = {
		1, //! 0001: N/1
		2, //! 001x: N/2,N/2
		3, //! 010x: N/4,N/4,N/2
		3, //! 011x: N/2,N/4,N/4
		4, //! 100x: N/8,N/8,N/4,N/2
		4, //! 101x: N/4,N/8,N/8,N/2
		4, //! 110x: N/2,N/8,N/8,N/4
		4, //! 111x: N/2,N/4,N/8,N/8
	};
	return Count[Decimation >> 1];
}

//! Subblock decimation pattern
//! Usage:
//!  PatternBase = &PatternTable[DecimationNybble >> 1];
//!  PatternNext = PatternBase;
//!  for(Coef in Coefs) {
//!   *Buffer = Coef; Buffer += BlockSize >> (*PatternNext++);
//!   if(Buffer >= BufferEnd) Buffer -= BufferSize-1, PatternNext = PatternBase;
//!  }
//! The decimation nybble comes from Block_Transform_GetLogOverlapScale(),
//! corresponding to the upper/second nybble for decimated blocks.
static inline const int8_t *ULC_Helper_SubBlockDecimationPattern(int Decimation) {
	static const int8_t Pattern[8][4] = {
		{0},          //! 0001: N/1
		{1, 1},       //! 001x: N/2,N/2
		{2, 2, 1},    //! 010x: N/4,N/4,N/2
		{1, 2, 2},    //! 011x: N/2,N/4,N/4
		{3, 3, 2, 1}, //! 100x: N/8,N/8,N/4,N/2
		{2, 3, 3, 1}, //! 101x: N/4,N/8,N/8,N/2
		{1, 3, 3, 2}, //! 110x: N/2,N/8,N/8,N/4
		{1, 2, 3, 3}, //! 111x: N/2,N/4,N/8,N/8
	};
	return Pattern[Decimation >> 1];
}

//! Subblock coefficient interleaving pattern
//! The decimation nybble comes from Block_Transform_GetLogOverlapScale(),
//! corresponding to the upper/second nybble for decimated blocks.
#define ULC_HELPER_SUBBLOCK_INTERLEAVE_MODULO 8u
static inline const int8_t *ULC_Helper_SubBlockInterleavePattern(int Decimation) {
	static const int8_t Mapping[8][8] = {
		{0,0,0,0,0,0,0,0}, //! 0001: N/1
		{0,1,0,1,0,1,0,1}, //! 001x: N/2,N/2
		{0,1,2,2,0,1,2,2}, //! 010x: N/4,N/4,N/2
		{0,0,1,2,0,0,1,2}, //! 011x: N/2,N/4,N/4
		{0,1,2,2,3,3,3,3}, //! 100x: N/8,N/8,N/4,N/2
		{0,0,1,2,3,3,3,3}, //! 101x: N/4,N/8,N/8,N/2
		{0,0,0,0,1,2,3,3}, //! 110x: N/2,N/8,N/8,N/4
		{0,0,0,0,1,1,2,3}, //! 111x: N/2,N/4,N/8,N/8
	};
	return Mapping[Decimation >> 1];
}

//! Transient subblock index from decimation pattern
//! Due to the way the decimation pattern is built, the
//! subblock index is conveniently obtained via a
//! POPCNT (minus one for the unary 'stop' bit).
static inline int ULC_Helper_TransientSubBlockIndex(int Decimation) {
	return __builtin_popcount(Decimation) - 1;
}

/**************************************/

//! Quick and dirty logarithm
//! Polynomial (with x = 0.0 .. 1.0):
//!  a*x + b*x^2 + c*x^3
//!  a = 1
//!  b = -15 - 3*Log[2] + 12*Log[4]
//!  c =  14 + 4*Log[2] - 12*Log[4]
//! This gives an approximation to Log[1+x], with PSNR = 59.3dB.
//! The exponent part is then trivial, forming:
//!  Log[2]*Exponent + Log[1+Mantissa]
//! NOTE: Returns Log[2^-126] for Log[0].
//! NOTE: Assumes x is normal (not NaN, +/-inf) and non-negative.
//! NOTE: Biased by 2^-126 to avoid dealing with subnormals.
static inline __attribute__((always_inline)) float ULC_FastLnApprox(float x) {
	union {
		float f;
		struct {
			uint32_t m:23;
			uint32_t e:9; //! Technically, {e:8,s:1}, but s=0 always
		};
	} v = {.f = x + 0x1.0p-126f};
	uint64_t m = 0xE348115793F6000ull - v.m*0x8C58828E7ull;      //! .61
		 m = ((uint64_t)v.m*v.m >> 14) * (m >> 29);          //! .64
		 m = ((uint64_t)v.m     << 41) - m;                  //! .64; final value for mantissa part
		 m = ((int)v.e - 127)*0xB17217F7D1CF78ll + (m >> 8); //! .56 (SIGNED!)
	return (int64_t)m * 0x1.0p-56f;
}

/**************************************/
//! EOF
/**************************************/
