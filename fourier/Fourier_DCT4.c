/**************************************/
//! ulc-codec: Ultra-Low-Complexity Audio Codec
//! Copyright (C) 2019, Ruben Nunez (Aikku; aik AT aol DOT com DOT au)
//! Refer to the project README file for license terms.
/**************************************/
#if defined(__AVX__)
# include <immintrin.h>
#endif
#if defined(__SSE__)
# include <xmmintrin.h>
#endif
/**************************************/
#include <math.h>
#include <stddef.h>
/**************************************/
#include "Fourier.h"
/**************************************/

//! DCT-IV (N=8)
static void DCT4_8(float *x) {
	const float sqrt1_2 = 0x1.6A09E667F3BCD0p-1;
	const float c1_3 = 0x1.D906BCF328D460p-1, s1_3 = 0x1.87DE2A6AEA9630p-2;
	const float c1_5 = 0x1.FD88DA3D125260p-1, s1_5 = 0x1.917A6BC29B42C0p-4;
	const float c3_5 = 0x1.E9F4156C62DDA0p-1, s3_5 = 0x1.294062ED59F050p-2;
	const float c5_5 = 0x1.C38B2F180BDB10p-1, s5_5 = 0x1.E2B5D3806F63B0p-2;
	const float c7_5 = 0x1.8BC806B1517410p-1, s7_5 = 0x1.44CF325091DD60p-1;

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

void Fourier_DCT4(float *Buf, float *Tmp, size_t N) {
	size_t i;

	//! Stop condition
	if(N == 8) {
		DCT4_8(Buf);
		return;
	}

	//! Perform rotation butterflies
	//!  u = R_n.x
	//! TODO: SIMD opt
	float c  = cos(M_PI_4 / N), s  = sin(M_PI_4 / N);
	float wc = cos(M_PI_2 / N), ws = sin(M_PI_2 / N);
	for(i=0;i<N/2;i++) {
		float a = Buf[i];
		float b = Buf[N-1-i];
		Tmp[i]     = ( c*a + s*b);
		Tmp[N-1-i] = (-s*a + c*b) * (i%2 ? (+1) : (-1));
		float _c = wc*c - ws*s;
		float _s = ws*c + wc*s;
		c = _c;
		s = _s;
	}

	//! Perform recursion
	//!  z1 = cos2([u_j][j=0..n/2-1],n/2)
	//!  z2 = cos2([u_j][j=n/2..n-1],n/2)
	Fourier_DCT2(Tmp,       Buf,       N/2);
	Fourier_DCT2(Tmp + N/2, Buf + N/2, N/2);

	//! Combine
	//!  w = U_n.(z1^T, z2^T)^T
	//!  y = (P_n)^T.w
	//! TODO: SIMD opt? How to even structure this loop?
	Buf[0] = Tmp[0];
	for(i=1;i<N/2;i++) {
		float a = Tmp[i];
		float b = Tmp[N-i] * (i%2 ? (-1) : (+1));
		Buf[i*2-1] = a + b;
		Buf[i*2+0] = a - b;
	}
	Buf[N-1] = Tmp[N/2];
}

/**************************************/
//! EOF
/**************************************/
