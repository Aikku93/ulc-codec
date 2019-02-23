/**************************************/
#pragma once
/**************************************/
#include <cstddef>
/**************************************/
//! DCT-II/DCT-IV algorithms from:
//!  "Signal Processing based on Stable radix-2 DCT I-IV Algorithms having Orthogonal Factors"
//!  DOI: 10.13001/1081-3810.3207
/**************************************/

namespace Fourier {
	//! Buf[N]
	//! Tmp[N]
	void DCT2(float Buf[], float Tmp[], size_t N);
	void DCT4(float Buf[], float Tmp[], size_t N);

	//! BufOut[N]
	//! BufIn[N]
	//! BufLap[N/2]
	//! BufTmp[N]
	void  MDCT(float BufOut[], const float BufIn[], float BufLap[], float BufTmp[], size_t N);
	void IMDCT(float BufOut[], const float BufIn[], float BufLap[], float BufTmp[], size_t N);
}

/**************************************/
//! EOF
/**************************************/
