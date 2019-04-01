/**************************************/
.section .iwram, "ax", %progbits
.balign 4
/**************************************/

@ r0: &Buf
@ r1: &Tmp
@ r2:  N

.arm
Fourier_DCT2:
	CMP	r2, #0x08
	BEQ	.LDCT2_8

.LButterflies:
	STMFD	sp!, {r4-sl,lr}
0:	SUB	sl, r2, r2, lsl #0x10
	ADD	r9, r0, r2, lsl #0x02   @ SrcHi = Buf+N
	ADD	r8, r1, r2, lsl #0x02-1 @ DstHi = Tmp+N/2
1:	LDMIA	r0!, {r4-r7}
	LDMDB	r9!, {r2,r3,ip,lr}
	ADD	r4, r4, lr
	ADD	r5, r5, ip
	ADD	r6, r6, r3
	ADD	r7, r7, r2
	STMIA	r1!, {r4-r7}
	SUB	r4, r4, lr, lsl #0x01
	SUB	r5, r5, ip, lsl #0x01
	SUB	r6, r6, r3, lsl #0x01
	SUB	r7, r7, r2, lsl #0x01
	STMIA	r8!, {r4-r7}
	ADDS	sl, sl, #0x08<<16
	BCC	1b
2:	LDMFD	sp!, {r4-r7}

@ r8: Tmp+N
@ r9: Buf+N/2
@ sl: N

.LRecurse:
	SUB	r0, r8, sl, lsl #0x02   @ Buf=Tmp
	SUB	r1, r9, sl, lsl #0x02-1 @ Tmp=Buf
	MOV	r2, sl, lsr #0x01       @ N=N/2
	BL	Fourier_DCT2
	SUB	r0, r8, sl, lsl #0x02-1 @ Buf=Tmp+N/2
	MOV	r1, r9                  @ Tmp=Buf+N/2
	MOV	r2, sl, lsr #0x01       @ N=N/2
	BL	Fourier_DCT4

.LMerge:
	STMFD	sp!, {r4-r7}
0:	SUB	r0, r9, sl, lsl #0x02-1 @ Dst=Buf
	SUB	r1, r8, sl, lsl #0x02   @ SrcLo=Tmp
	SUB	r2, r8, sl, lsl #0x02-1 @ SrcHi=Tmp+N/2
1:	LDMIA	r1!, {r3,r5,r7,r9}
	LDMIA	r2!, {r4,r6,r8,lr}
	STMIA	r0!, {r3-r9,lr}
	LDMIA	r1!, {r3,r5,r7,r9}
	LDMIA	r2!, {r4,r6,r8,lr}
	STMIA	r0!, {r3-r9,lr}
	SUBS	sl, sl, #0x10
	BNE	1b
2:	LDMFD	sp!, {r4-sl,lr}
	BX	lr

/**************************************/

@ r0: &Buf

.LDCT2_8:
	STMFD	sp!, {r4-fp,lr}
	LDMIA	r0, {r1-r8}
0:	ADD	r1, r1, r8            @ s07 -> r1
	ADD	r2, r2, r7            @ s16 -> r2
	ADD	r3, r3, r6            @ s25 -> r3
	ADD	r4, r4, r5            @ s34 -> r4
	SUB	r5, r4, r5, lsl #0x01 @ d34 -> r5
	SUB	r6, r3, r6, lsl #0x01 @ d25 -> r6
	SUB	r7, r2, r7, lsl #0x01 @ d16 -> r7
	SUB	r8, r1, r8, lsl #0x01 @ d07 -> r8
1:	ADD	r1, r1, r4            @ ss07s34 -> r1
	ADD	r2, r2, r3            @ ss16s25 -> r2
	SUB	r3, r2, r3, lsl #0x01 @ ds16s25 -> r3
	SUB	r4, r1, r4, lsl #0x01 @ ds07s34 -> r4
	MOV	fp, #0x3500           @ c3_4[.14] -> fp
	MOV	lr, #0x2300           @ s3_4[.14] -> lr
	ORR	fp, fp, #0x37
	ORR	lr, lr, #0x8E
	MUL	r9, r5, fp            @ d34d07x =  c3_4*d34 + s3_4*d07 -> r9 [.14]
	MUL	sl, r8, fp            @ d34d07y = -s3_4*d34 + c3_4*d07 -> sl [.14]
	MLA	r9, r8, lr, r9
	MUL	lr, r5, lr
	SUB	sl, sl, lr
	MOV	fp, #0x7D00           @ c1_4[.15] -> fp
	MOV	lr, #0x6300           @ s1_4[.17] -> lr
	ORR	fp, fp, #0x8A
	ORR	lr, lr, #0xE3
	MUL	r5, r6, fp            @ d25d16x =  c1_4*d25 + s1_4*d16 -> r5 [.15]
	MUL	r8, r7, fp            @ d25d16y = -s1_4*d25 + c1_4*d16 -> r8 [.15]
	MUL	fp, r6, lr
	MUL	lr, r7, lr
	ADD	r5, r5, lr, asr #0x02
	SUB	r8, r8, fp, asr #0x02
2:	MOV	lr, #0x2D00           @ sqrt1_2[.14] -> lr
	ORR	lr, lr, #0x41
	ADD	r1, r1, r2            @ a0 =       ss07s34 +      ss16s25 -> r1 = X0
	SUB	r7, r1, r2, lsl #0x01 @ b0 =       ss07s34 -      ss16s25 -> r7 = X4/sqrt1_2
	MUL	r2, r7, lr
	MOV	r7, r2, asr #0x0E     @ [X4 -> r7]
	MOV	sl, sl, asr #0x0E
	ADD	sl, sl, r5, asr #0x0F @ a1 =       d34d07y +      d25d16x -> sl
	SUB	r5, sl, r5, asr #0x0E @ c1 =       d34d07y -      d25d16x -> r5 = X3
	MOV	r9, r9, asr #0x0E
	ADD	r9, r9, r8, asr #0x0F @ d1 =       d34d07x +      d25d16y -> r9
	SUB	r8, r9, r8, asr #0x0E @ b1 =       d34d07x -      d25d16y -> r8 = X5
	ADD	r2, sl, r9            @ (a1+d1)*sqrt1_2 = X1 -> r2
	SUB	ip, sl, r9            @ (a1-d1)*sqrt1_2 = X7 -> ip
	MUL	r9, r2, lr
	MUL	sl, ip, lr
	MOV	r2, r9, asr #0x0E
	MOV	ip, sl, asr #0x0E
	MOV	fp, #0x3000           @ c6_4[.15] -> fp
	MOV	lr, #0x3B00           @ s6_4[.14] -> lr
	ORR	fp, fp, #0xFB
	ORR	lr, lr, #0x21
	MUL	r9, r3, fp            @ c0 =  c6_4*ds16s25 + s6_4*ds07s34 -> r4 = X2
	MUL	sl, r4, fp            @ d0 = -s6_4*ds16s25 + c6_4*ds07s34 -> r9 = X6
	MUL	r6, r4, lr
	MUL	lr, r3, lr
	ADD	r4, r9, r6, lsl #0x01
	MOV	r4, r4, asr #0x0F
	SUB	r9, sl, lr, lsl #0x01
	MOV	r9, r9, asr #0x0F
	STMIA	r0, {r1,r2,r4,r5,r7,r8,r9,ip}
2:	LDMFD	sp!, {r4-fp,lr}
	BX	lr

/**************************************/
.size   Fourier_DCT2, .-Fourier_DCT2
.global Fourier_DCT2
/**************************************/
/* EOF                                */
/**************************************/
