/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2022, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include "Fourier.h"
#include "FourierHelper.h"
/**************************************/

//! DCT-IV (N=8)
static void DCT4T_8(float *x) {
	FOURIER_ASSUME_ALIGNED(x, 32);

	const float sqrt1_2 = 0x1.6A09E6p-1f;
	const float c1_3 = 0x1.D906BCp-1f, s1_3 = 0x1.87DE2Ap-2f;
	const float c1_5 = 0x1.FD88DAp-1f, s1_5 = 0x1.917A6Cp-4f;
	const float c3_5 = 0x1.E9F416p-1f, s3_5 = 0x1.294062p-2f;
	const float c5_5 = 0x1.C38B30p-1f, s5_5 = 0x1.E2B5D4p-2f;
	const float c7_5 = 0x1.8BC806p-1f, s7_5 = 0x1.44CF32p-1f;

	float sx = x[0];
	float tx = (x[2] + x[1]);
	float vy = (x[2] - x[1]);
	float sy = (x[3] + x[4]) * sqrt1_2;
	float uy = (x[3] - x[4]) * sqrt1_2;
	float ty = (x[6] + x[5]);
	float vx = (x[6] - x[5]);
	float ux = x[7];

	float saxdx = sx + sy;
	float sbxcx = sx - sy;
	float daxdx = c1_3*tx + s1_3*ty;
	float dbxcx = s1_3*tx - c1_3*ty;
	float sdyay = ux + uy;
	float scyby = ux - uy;
	float ddyay = c1_3*vx + s1_3*vy;
	float dcyby = s1_3*vx - c1_3*vy;

	float ax = saxdx + daxdx;
	float dx = saxdx - daxdx;
	float bx = sbxcx + dbxcx;
	float cx = sbxcx - dbxcx;
	float dy = sdyay + ddyay;
	float ay = sdyay - ddyay;
	float cy = scyby + dcyby;
	float by = scyby - dcyby;

	x[0] = c1_5*ax + s1_5*ay;
	x[7] = s1_5*ax - c1_5*ay;
	x[1] = c3_5*bx - s3_5*by;
	x[6] = s3_5*bx + c3_5*by;
	x[2] = c5_5*cx + s5_5*cy;
	x[5] = s5_5*cx - c5_5*cy;
	x[3] = c7_5*dx - s7_5*dy;
	x[4] = s7_5*dx + c7_5*dy;
}

/**************************************/

void Fourier_DCT4T(float *Buf, float *Tmp, int N) {
	int i;
	FOURIER_ASSUME_ALIGNED(Buf, 32);
	FOURIER_ASSUME_ALIGNED(Tmp, 32);
	FOURIER_ASSUME(N >= 8);

	//! Stop condition
	if(N == 8) {
		DCT4T_8(Buf);
		return;
	}

	{
		      float *DstLo = Tmp;
		      float *DstHi = Tmp + N;
		const float *Src   = Buf;
		*DstLo++ = *Src++ * 2.0f;
#if FOURIER_VSTRIDE > 1
		{
			Fourier_Vec_t a,  b;
			Fourier_Vec_t t0, t1;
			for(i=0;i<N/2-FOURIER_VSTRIDE;i+=FOURIER_VSTRIDE) {
				a  = FOURIER_VLOADU(Src); Src += FOURIER_VSTRIDE;
				b  = FOURIER_VLOADU(Src); Src += FOURIER_VSTRIDE;
				FOURIER_VSPLIT_EVEN_ODD(a, b, &a, &b);
				t0 = FOURIER_VADD(a, b);
				t1 = FOURIER_VSUB(a, b);
				t1 = FOURIER_VREVERSE(t1);
				FOURIER_VSTOREU(DstLo, t0); DstLo += FOURIER_VSTRIDE;
				DstHi -= FOURIER_VSTRIDE; FOURIER_VSTORE(DstHi, t1);
			}
		}
#else
		i = 0;
#endif
		float a, b;
		for(;i<N/2-1;i++) {
			a = *Src++;
			b = *Src++;
			*DstLo++ = a + b;
			*--DstHi = a - b;
		}
		*--DstHi = *Src++ * 2.0f;
	}

	Fourier_DCT3(Tmp,       Buf,       N/2);
	Fourier_DCT3(Tmp + N/2, Buf + N/2, N/2);

	{
		const float *SrcLo = Tmp;
		const float *SrcHi = Tmp + N/2;
		      float *DstLo = Buf;
		      float *DstHi = Buf + N;
#if FOURIER_VSTRIDE > 1
		Fourier_Vec_t a, b;
		Fourier_Vec_t t0, t1 = FOURIER_VMUL(FOURIER_VSET1(1.0f/N), FOURIER_VADD(FOURIER_VSET_LINEAR_RAMP(), FOURIER_VSET1(0.5f)));
		Fourier_Vec_t c  = Fourier_Cos(t1);
		Fourier_Vec_t s  = Fourier_Sin(t1);
		Fourier_Vec_t wc = Fourier_Cos(FOURIER_VSET1((float)FOURIER_VSTRIDE / N));
		Fourier_Vec_t ws = Fourier_Sin(FOURIER_VSET1((float)FOURIER_VSTRIDE / N));
		for(i=0;i<N/2;i+=FOURIER_VSTRIDE) {
			a  = FOURIER_VLOAD(SrcLo); SrcLo += FOURIER_VSTRIDE;
			b  = FOURIER_VLOAD(SrcHi); SrcHi += FOURIER_VSTRIDE;
			t0 = FOURIER_VMUL(s, b);
			t1 = FOURIER_VMUL(c, b);
			t0 = FOURIER_VNEGATE_ODD(t0);
			t1 = FOURIER_VNEGATE_ODD(t1);
			t0 = FOURIER_VFMA(c, a, t0);
			t1 = FOURIER_VFMS(s, a, t1);
			t1 = FOURIER_VREVERSE(t1);
			FOURIER_VSTORE(DstLo, t0); DstLo += FOURIER_VSTRIDE;
			DstHi -= FOURIER_VSTRIDE; FOURIER_VSTORE(DstHi, t1);
			t0 = c;
			t1 = s;
			c = FOURIER_VNFMA(t1, ws, FOURIER_VMUL(t0, wc));
			s = FOURIER_VFMA (t1, wc, FOURIER_VMUL(t0, ws));
		}
#else
		float a, b;
		float c  = Fourier_Cos(0.5f / N);
		float s  = Fourier_Sin(0.5f / N);
		float wc = Fourier_Cos(1.0f / N);
		float ws = Fourier_Sin(1.0f / N);
		for(i=0;i<N/2;i+=2) {
			a = *SrcLo++;
			b = *SrcHi++;
			*DstLo++ =  c*a + s*b;
			*--DstHi =  s*a - c*b;
			a = c;
			b = s;
			c = wc*a - ws*b;
			s = ws*a + wc*b;

			a = *SrcLo++;
			b = *SrcHi++;
			*DstLo++ =  c*a - s*b;
			*--DstHi =  s*a + c*b;
			a = c;
			b = s;
			c = wc*a - ws*b;
			s = ws*a + wc*b;
		}
#endif
	}
}

/**************************************/
//! EOF
/**************************************/
