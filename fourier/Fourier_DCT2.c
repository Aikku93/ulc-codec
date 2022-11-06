/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2022, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include "Fourier.h"
#include "FourierHelper.h"
/**************************************/

//! DCT-II (N=8)
static void DCT2_8(float *x) {
	FOURIER_ASSUME_ALIGNED(x, 32);

	const float sqrt1_2 = 0x1.6A09E6p-1f;
	const float c1_4 = 0x1.F6297Cp-1f, s1_4 = 0x1.8F8B84p-3f;
	const float c3_4 = 0x1.A9B662p-1f, s3_4 = 0x1.1C73B4p-1f;
	const float c6_4 = 0x1.87DE2Ap-2f, s6_4 = 0x1.D906BCp-1f;

	//! First stage butterflies (DCT2_8)
	float s07 = x[0]+x[7];
	float d07 = x[0]-x[7];
	float s16 = x[1]+x[6];
	float d16 = x[1]-x[6];
	float s25 = x[2]+x[5];
	float d25 = x[2]-x[5];
	float s34 = x[3]+x[4];
	float d34 = x[3]-x[4];

	//! Second stage (DCT2_4, DCT4_4)
	float ss07s34 = s07+s34;
	float ds07s34 = s07-s34;
	float ss16s25 = s16+s25;
	float ds16s25 = s16-s25;
	float d34d07x =  c3_4*d34 + s3_4*d07;
	float d34d07y = -s3_4*d34 + c3_4*d07;
	float d25d16x =  c1_4*d25 + s1_4*d16;
	float d25d16y = -s1_4*d25 + c1_4*d16;

	//! Third stage (rotation butterflies; DCT2_2, DCT4_2, DCT2_2, DCT2_2)
	float a0 =       ss07s34 +      ss16s25;
	float b0 =       ss07s34 -      ss16s25;
	float c0 =  c6_4*ds16s25 + s6_4*ds07s34;
	float d0 = -s6_4*ds16s25 + c6_4*ds07s34;
	float a1 =       d34d07y +      d25d16x;
	float c1 =       d34d07y -      d25d16x;
	float d1 =       d34d07x +      d25d16y;
	float b1 =       d34d07x -      d25d16y;

	//! Permute and final DCT4 stage
	x[0] = a0;
	x[4] = b0 * sqrt1_2;
	x[2] = c0;
	x[6] = d0;
	x[1] = (a1 + d1) * sqrt1_2;
	x[5] = b1;
	x[3] = c1;
	x[7] = (a1 - d1) * sqrt1_2;
}

/**************************************/

void Fourier_DCT2(float *Buf, float *Tmp, int N) {
	int i;
	FOURIER_ASSUME_ALIGNED(Buf, 32);
	FOURIER_ASSUME_ALIGNED(Tmp, 32);
	FOURIER_ASSUME(N >= 8);

	//! Stop condition
	if(N == 8) {
		DCT2_8(Buf);
		return;
	}

	//! Perform butterflies
	//!  u = H_n.x
	{
		const float *SrcLo = Buf;
		const float *SrcHi = Buf + N;
		      float *DstLo = Tmp;
		      float *DstHi = Tmp + N/2;
#if FOURIER_VSTRIDE > 1
		Fourier_Vec_t a, b;
		Fourier_Vec_t s, d;
		for(i=0;i<N/2;i+=FOURIER_VSTRIDE) {
			SrcHi -= FOURIER_VSTRIDE; b = FOURIER_VREVERSE(FOURIER_VLOAD(SrcHi));
			a = FOURIER_VLOAD(SrcLo); SrcLo += FOURIER_VSTRIDE;
			s = FOURIER_VADD(a, b);
			d = FOURIER_VSUB(a, b);
			FOURIER_VSTORE(DstLo, s); DstLo += FOURIER_VSTRIDE;
			FOURIER_VSTORE(DstHi, d); DstHi += FOURIER_VSTRIDE;
		}
#else
		float a, b;
		for(i=0;i<N/2;i++) {
			a = *SrcLo++;
			b = *--SrcHi;
			*DstLo++ = a + b;
			*DstHi++ = a - b;
		}
#endif
	}

	//! Perform recursion
	//!  z1 = cos2([u_j][j=0..n/2-1],n/2)
	//!  z2 = cos4([u_j][j=n/2..n-1],n/2)
	Fourier_DCT2(Tmp,       Buf,       N/2);
	Fourier_DCT4(Tmp + N/2, Buf + N/2, N/2);

	//! Combine
	//!  y = (P_n)^T.(z1^T, z2^T)^T
	{
		const float *SrcLo = Tmp;
		const float *SrcHi = Tmp + N/2;
		      float *Dst   = Buf;
#if FOURIER_VSTRIDE > 1
		Fourier_Vec_t a, b;
		for(i=0;i<N/2;i+=FOURIER_VSTRIDE) {
			a = FOURIER_VLOAD(SrcLo); SrcLo += FOURIER_VSTRIDE;
			b = FOURIER_VLOAD(SrcHi); SrcHi += FOURIER_VSTRIDE;
			FOURIER_VINTERLEAVE(a, b, &a, &b);
			FOURIER_VSTORE(Dst+0,               a);
			FOURIER_VSTORE(Dst+FOURIER_VSTRIDE, b); Dst += 2*FOURIER_VSTRIDE;
		}
#else
		float a, b;
		for(i=0;i<N/2;i++) {
			a = *SrcLo++;
			b = *SrcHi++;
			*Dst++ = a;
			*Dst++ = b;
		}
#endif
	}
}

/**************************************/
//! EOF
/**************************************/
