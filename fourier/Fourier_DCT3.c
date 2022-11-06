/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2022, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include "Fourier.h"
#include "FourierHelper.h"
/**************************************/

//! DCT-III (N=8)
static void DCT3_8(float *x) {
	FOURIER_ASSUME_ALIGNED(x, 32);

	const float sqrt1_2 = 0x1.6A09E6p-1f;
	const float c1_4 = 0x1.F6297Cp-1f, s1_4 = 0x1.8F8B84p-3f;
	const float c3_4 = 0x1.A9B662p-1f, s3_4 = 0x1.1C73B4p-1f;
	const float c6_4 = 0x1.87DE2Ap-2f, s6_4 = 0x1.D906BCp-1f;

	float a0 = x[0] * 0.5f;
	float b0 = x[4] * sqrt1_2;
	float c0 = x[2];
	float d0 = x[6];
	float a1 = (x[1] + x[7]) * sqrt1_2;
	float d1 = (x[1] - x[7]) * sqrt1_2;
	float b1 = x[5];
	float c1 = x[3];

	float ss07s34 = a0 + b0;
	float ss16s25 = a0 - b0;
	float ds16s25 = c6_4*c0 - s6_4*d0;
	float ds07s34 = s6_4*c0 + c6_4*d0;
	float d34d07y = a1 + c1;
	float d25d16x = a1 - c1;
	float d34d07x = d1 + b1;
	float d25d16y = d1 - b1;

	float s07 = ss07s34 + ds07s34;
	float s34 = ss07s34 - ds07s34;
	float s16 = ss16s25 + ds16s25;
	float s25 = ss16s25 - ds16s25;
	float d34 = c3_4*d34d07x - s3_4*d34d07y;
	float d07 = s3_4*d34d07x + c3_4*d34d07y;
	float d25 = c1_4*d25d16x - s1_4*d25d16y;
	float d16 = s1_4*d25d16x + c1_4*d25d16y;

	x[0] = s07 + d07;
	x[7] = s07 - d07;
	x[1] = s16 + d16;
	x[6] = s16 - d16;
	x[2] = s25 + d25;
	x[5] = s25 - d25;
	x[3] = s34 + d34;
	x[4] = s34 - d34;
}

/**************************************/

void Fourier_DCT3(float *Buf, float *Tmp, int N) {
	int i;
	FOURIER_ASSUME_ALIGNED(Buf, 32);
	FOURIER_ASSUME_ALIGNED(Tmp, 32);
	FOURIER_ASSUME(N >= 8);

	//! Stop condition
	if(N == 8) {
		DCT3_8(Buf);
		return;
	}

	{
		      float *Dst = Tmp;
		const float *Src = Buf;
#if FOURIER_VSTRIDE > 1
		Fourier_Vec_t a, b;
		for(i=0;i<N/2;i+=FOURIER_VSTRIDE) {
			a = FOURIER_VLOAD(Src); Src += FOURIER_VSTRIDE;
			b = FOURIER_VLOAD(Src); Src += FOURIER_VSTRIDE;
			FOURIER_VSPLIT_EVEN_ODD(a, b, &a, &b);
			FOURIER_VSTORE(Dst,     a);
			FOURIER_VSTORE(Dst+N/2, b); Dst += FOURIER_VSTRIDE;
		}
#else
		float a, b;
		for(i=0;i<N/2;i++) {
			a = *Src++;
			b = *Src++;
			Dst[N/2] = b;
			*Dst++   = a;
		}
#endif
	}

	Fourier_DCT3 (Tmp,       Buf,       N/2);
	Fourier_DCT4T(Tmp + N/2, Buf + N/2, N/2);

	{
		const float *SrcLo = Tmp;
		const float *SrcHi = Tmp + N/2;
		      float *DstLo = Buf;
		      float *DstHi = Buf + N;
#if FOURIER_VSTRIDE > 1
		Fourier_Vec_t a, b;
		Fourier_Vec_t s, d;
		for(i=0;i<N/2;i+=FOURIER_VSTRIDE) {
			a = FOURIER_VLOAD(SrcLo); SrcLo += FOURIER_VSTRIDE;
			b = FOURIER_VLOAD(SrcHi); SrcHi += FOURIER_VSTRIDE;
			s = FOURIER_VADD(a, b);
			d = FOURIER_VSUB(a, b);
			d = FOURIER_VREVERSE(d);
			FOURIER_VSTORE(DstLo, s); DstLo += FOURIER_VSTRIDE;
			DstHi -= FOURIER_VSTRIDE; FOURIER_VSTORE(DstHi, d);
		}
#else
		float a, b;
		for(i=0;i<N/2;i++) {
			a = *SrcLo++;
			b = *SrcHi++;
			*DstLo++ = a + b;
			*--DstHi = a - b;
		}
#endif
	}
}

/**************************************/
//! EOF
/**************************************/
