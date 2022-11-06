/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2022, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#include "Fourier.h"
#include "FourierHelper.h"
/**************************************/

//! DCT-IV (N=8)
static void DCT4_8(float *x) {
	FOURIER_ASSUME_ALIGNED(x, 32);

	const float sqrt1_2 = 0x1.6A09E6p-1f;
	const float c1_3 = 0x1.D906BCp-1f, s1_3 = 0x1.87DE2Ap-2f;
	const float c1_5 = 0x1.FD88DAp-1f, s1_5 = 0x1.917A6Cp-4f;
	const float c3_5 = 0x1.E9F416p-1f, s3_5 = 0x1.294062p-2f;
	const float c5_5 = 0x1.C38B30p-1f, s5_5 = 0x1.E2B5D4p-2f;
	const float c7_5 = 0x1.8BC806p-1f, s7_5 = 0x1.44CF32p-1f;

	//! First stage (rotation butterflies; DCT4_8)
	float ax =  c1_5*x[0] + s1_5*x[7];
	float ay =  s1_5*x[0] - c1_5*x[7];
	float bx =  c3_5*x[1] + s3_5*x[6];
	float by = -s3_5*x[1] + c3_5*x[6];
	float cx =  c5_5*x[2] + s5_5*x[5];
	float cy =  s5_5*x[2] - c5_5*x[5];
	float dx =  c7_5*x[3] + s7_5*x[4];
	float dy = -s7_5*x[3] + c7_5*x[4];

	//! Second stage (butterflies; DCT2_4)
	float saxdx = ax + dx;
	float daxdx = ax - dx;
	float sbxcx = bx + cx;
	float dbxcx = bx - cx;
	float sdyay = dy + ay;
	float ddyay = dy - ay;
	float scyby = cy + by;
	float dcyby = cy - by;

	//! Third stage (rotation butterflies; DCT2_2, DCT4_2)
	float sx =      saxdx +      sbxcx;
	float sy =      saxdx -      sbxcx;
	float tx = c1_3*daxdx + s1_3*dbxcx;
	float ty = s1_3*daxdx - c1_3*dbxcx;
	float ux =      sdyay +      scyby;
	float uy =      sdyay -      scyby;
	float vx = c1_3*ddyay + s1_3*dcyby;
	float vy = s1_3*ddyay - c1_3*dcyby;

	//! Permute and final DCT4 stage
	x[0] = sx;
	x[1] = (tx - vy);
	x[2] = (tx + vy);
	x[3] = (sy + uy) * sqrt1_2;
	x[4] = (sy - uy) * sqrt1_2;
	x[5] = (ty - vx);
	x[6] = (ty + vx);
	x[7] = ux;
}

/**************************************/

void Fourier_DCT4(float *Buf, float *Tmp, int N) {
	int i;
	FOURIER_ASSUME_ALIGNED(Buf, 32);
	FOURIER_ASSUME_ALIGNED(Tmp, 32);
	FOURIER_ASSUME(N >= 8);

	//! Stop condition
	if(N == 8) {
		DCT4_8(Buf);
		return;
	}

	//! Perform rotation butterflies
	//!  u = R_n.x
	{
		const float *SrcLo = Buf;
		const float *SrcHi = Buf + N;
		      float *DstLo = Tmp;
		      float *DstHi = Tmp + N/2;
#if FOURIER_VSTRIDE > 1
		Fourier_Vec_t a, b;
		Fourier_Vec_t t0, t1 = FOURIER_VMUL(FOURIER_VSET1(1.0f/N), FOURIER_VADD(FOURIER_VSET_LINEAR_RAMP(), FOURIER_VSET1(0.5f)));
		Fourier_Vec_t c  = Fourier_Cos(t1);
		Fourier_Vec_t s  = Fourier_Sin(t1);
		Fourier_Vec_t wc = Fourier_Cos(FOURIER_VSET1((float)FOURIER_VSTRIDE / N));
		Fourier_Vec_t ws = Fourier_Sin(FOURIER_VSET1((float)FOURIER_VSTRIDE / N));
		for(i=0;i<N/2;i+=FOURIER_VSTRIDE) {
			SrcHi -= FOURIER_VSTRIDE; b = FOURIER_VREVERSE(FOURIER_VLOAD(SrcHi));
			a = FOURIER_VLOAD(SrcLo); SrcLo += FOURIER_VSTRIDE;
			t1 = FOURIER_VMUL(s, a);
			t0 = FOURIER_VMUL(c, a);
			t1 = FOURIER_VNFMA(c, b, t1);
			t0 = FOURIER_VFMA (s, b, t0);
			t1 = FOURIER_VNEGATE_ODD(t1);
			FOURIER_VSTORE(DstLo, t0); DstLo += FOURIER_VSTRIDE;
			FOURIER_VSTORE(DstHi, t1); DstHi += FOURIER_VSTRIDE;
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
			b = *--SrcHi;
			*DstLo++ =  c*a + s*b;
			*DstHi++ =  s*a - c*b;
			a = c;
			b = s;
			c = wc*a - ws*b;
			s = ws*a + wc*b;

			a = *SrcLo++;
			b = *--SrcHi;
			*DstLo++ =  c*a + s*b;
			*DstHi++ = -s*a + c*b; //! <- Sign-flip for DST
			a = c;
			b = s;
			c = wc*a - ws*b;
			s = ws*a + wc*b;
		}
#endif
	}

	//! Perform recursion
	//!  z1 = cos2([u_j][j=0..n/2-1],n/2)
	//!  z2 = cos2([u_j][j=n/2..n-1],n/2)
	Fourier_DCT2(Tmp,       Buf,       N/2);
	Fourier_DCT2(Tmp + N/2, Buf + N/2, N/2);

	//! Combine
	//!  w = U_n.(z1^T, z2^T)^T
	//!  y = (P_n)^T.w
	{
		const float *TmpLo = Tmp;
		const float *TmpHi = Tmp + N;
		      float *Dst   = Buf;
		*Dst++ = *TmpLo++;
#if FOURIER_VSTRIDE > 1
		{
			Fourier_Vec_t a, b;
			Fourier_Vec_t t0, t1;
			for(i=0;i<N/2-FOURIER_VSTRIDE;i+=FOURIER_VSTRIDE) {
				TmpHi -= FOURIER_VSTRIDE; b = FOURIER_VREVERSE(FOURIER_VLOAD(TmpHi));
				a = FOURIER_VLOADU(TmpLo); TmpLo += FOURIER_VSTRIDE;
				t0 = FOURIER_VADD(a, b);
				t1 = FOURIER_VSUB(a, b);
				FOURIER_VINTERLEAVE(t0, t1, &a, &b);
				FOURIER_VSTOREU(Dst, a); Dst += FOURIER_VSTRIDE;
				FOURIER_VSTOREU(Dst, b); Dst += FOURIER_VSTRIDE;
			}
		}
#else
		i = 0;
#endif
		float a, b;
		for(;i<N/2-1;i++) {
			a = *TmpLo++;
			b = *--TmpHi;
			*Dst++ = a + b;
			*Dst++ = a - b;
		}
		*Dst++ = *--TmpHi;
	}
}

/**************************************/
//! EOF
/**************************************/
