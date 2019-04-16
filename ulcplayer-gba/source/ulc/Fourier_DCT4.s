/**************************************/
.section .iwram, "ax", %progbits
.balign 4
/**************************************/

@ r0: &Buf
@ r1: &Tmp
@ r2:  N (<=800h)
@ NOTE: Caller must be ARM code

.arm
Fourier_DCT4:
	CMP	r2, #0x08
	BEQ	.LDCT4_8

.LButterflies:
	STMFD	sp!, {r4-fp,lr}
0:	ADD	r8, r0, r2, lsl #0x02   @ SrcHi = Tmp+N
	ADD	r9, r1, r2, lsl #0x02-1 @ DstHi = Tmp+N/2
	SUB	sl, r2, r2, lsl #0x10
	LDR	fp, =Fourier_DCT4_CosSin
	CMP	r2, #0x10
	ADDHI	fp, fp, #0x10*2
	CMP	r2, #0x20
	ADDHI	fp, fp, #0x20*2
	CMP	r2, #0x40
	ADDHI	fp, fp, #0x40*2
	CMP	r2, #0x80
	ADDHI	fp, fp, #0x80*2
	CMP	r2, #0x0100
	ADDHI	fp, fp, #0x0100*2
	CMP	r2, #0x0200
	ADDHI	fp, fp, #0x0200*2
	CMP	r2, #0x0400
	ADDHI	fp, fp, #0x0400*2
1:	LDMIA	fp!, {ip,lr}          @ cs -> ip,lr
	LDMIA	r0!, {r2-r3}          @ a = *SrcLo++
	LDMDB	r8!, {r4-r5}          @ b = *--SrcHi
	MOV	r7, ip, lsr #0x10     @ s -> r7
	BIC	ip, ip, r7, lsl #0x10 @ c -> ip
	MUL	r6, ip, r2            @ *DstLo++ =  c*a + s*b -> r6
	MUL	r2, r7, r2            @ *DstHi++ =  s*a - c*b -> r2
	MLA	r6, r7, r5, r6
	MUL	r7, ip, r5
	SUB	r2, r2, r7
	MOV	ip, lr, lsr #0x10     @ s -> ip
	BIC	lr, lr, ip, lsl #0x10 @ c -> lr
	MUL	r7, lr, r3            @ *DstLo++ =  c*a + s*b -> r7
	MUL	r3, ip, r3            @ *DstHi++ = -s*a + c*b -> r3
	MLA	r7, ip, r4, r7
	MUL	ip, lr, r4
	RSB	r3, r3, ip
	MOV	r2, r2, asr #0x0F
	MOV	r3, r3, asr #0x0F
	MOV	r6, r6, asr #0x0F
	MOV	r7, r7, asr #0x0F
	STMIA	r1!, {r6-r7}
	STMIA	r9!, {r2-r3}
2:	LDMIA	fp!, {ip,lr}          @ cs -> ip,lr
	LDMIA	r0!, {r2-r3}          @ a = *SrcLo++
	LDMDB	r8!, {r4-r5}          @ b = *--SrcHi
	MOV	r7, ip, lsr #0x10     @ s -> r7
	BIC	ip, ip, r7, lsl #0x10 @ c -> ip
	MUL	r6, ip, r2            @ *DstLo++ =  c*a + s*b -> r6
	MUL	r2, r7, r2            @ *DstHi++ =  s*a - c*b -> r2
	MLA	r6, r7, r5, r6
	MUL	r7, ip, r5
	SUB	r2, r2, r7
	MOV	ip, lr, lsr #0x10     @ s -> ip
	BIC	lr, lr, ip, lsl #0x10 @ c -> lr
	MUL	r7, lr, r3            @ *DstLo++ =  c*a + s*b -> r7
	MUL	r3, ip, r3            @ *DstHi++ = -s*a + c*b -> r3
	MLA	r7, ip, r4, r7
	MUL	ip, lr, r4
	RSB	r3, r3, ip
	MOV	r2, r2, asr #0x0F
	MOV	r3, r3, asr #0x0F
	MOV	r6, r6, asr #0x0F
	MOV	r7, r7, asr #0x0F
	STMIA	r1!, {r6-r7}
	STMIA	r9!, {r2-r3}
3:	ADDS	sl, sl, #0x08<<16
	BCC	1b
0:	LDMFD	sp!, {r4-r7}

@ r8: Buf+N/2
@ r9: Tmp+N
@ sl: N

.LRecurse:
	SUB	r0, r9, sl, lsl #0x02   @ Buf=Tmp
	SUB	r1, r8, sl, lsl #0x02-1 @ Tmp=Buf
	MOV	r2, sl, lsr #0x01       @ N=N/2
	BL	Fourier_DCT2
	SUB	r0, r9, sl, lsl #0x02-1 @ Buf=Tmp+N/2
	MOV	r1, r8                  @ Tmp=Buf+N/2
	MOV	r2, sl, lsr #0x01       @ N=N/2
	BL	Fourier_DCT2

.LMerge:
	STMFD	sp!, {r4-r7}
	SUB	ip, r8, sl, lsl #0x02-1 @ Dst=Buf -> ip
	SUB	r8, r9, sl, lsl #0x02   @ SrcLo=Tmp -> r8, SrcHi=Tmp+N -> r9
0:	LDR	r0, [r8], #0x04
	STR	r0, [ip], #0x04
	SUB	sl, sl, #0x08
1:	LDMIA	r8!, {r0-r3}
	LDMDB	r9!, {r4-r7}
	ADD	r0, r0, r7
	SUB	r7, r0, r7, lsl #0x01
	ADD	r1, r1, r6
	SUB	r6, r1, r6, lsl #0x01
	ADD	r2, r2, r5
	SUB	r5, r2, r5, lsl #0x01
	ADD	r3, r3, r4
	SUB	r4, r3, r4, lsl #0x01
	STMIA	ip!, {r0,r7}
	STMIA	ip!, {r1,r6}
	STMIA	ip!, {r2,r5}
	STMIA	ip!, {r3,r4}
	SUBS	sl, sl, #0x08
	BNE	1b
2:	LDMIA	r8!, {r0-r2}
	LDMDB	r9!, {r4-r7}
	ADD	r0, r0, r7
	SUB	r7, r0, r7, lsl #0x01
	ADD	r1, r1, r6
	SUB	r6, r1, r6, lsl #0x01
	ADD	r2, r2, r5
	SUB	r3, r2, r5, lsl #0x01
	STMIA	ip!, {r0,r7}
	STMIA	ip!, {r1,r6}
	STMIA	ip!, {r2-r4}
3:	LDMFD	sp!, {r4-fp,pc}

/**************************************/

@ r0: &Buf
@ c1_5: 3FB1h [.14]
@ s1_5: 645Fh [.18]
@ c3_5: 7A7Dh [.15]
@ s3_5: 04A5h [.12]
@ c5_5: 70E3h [.15]
@ s5_5: 78ADh [.16]
@ c7_5: 3179h [.14]
@ s7_5: 144Dh [.13]

.LDCT4_8:
	STMFD	sp!, {r4-fp,lr}
	LDMIA	r0, {r1-r8}
0:	MOV	ip, #0x6400           @ s1_5[.18] -> ip
	ORR	ip, ip, #0x5F
	RSB	r9, r1, r1, lsl #0x04 @ ax =  c1_5*x[0] + s1_5*x[7] -> r9 [.14]
	ADD	r9, r9, r1, lsl #0x06
	RSB	r9, r9, r1, lsl #0x0E
	MUL	sl, ip, r1            @ ay =  s1_5*x[0] - c1_5*x[7] -> r1 [.14]
	MUL	r1, ip, r8
	ADD	r9, r9, r1, asr #0x04
	RSB	r1, r8, r8, lsl #0x04
	ADD	r1, r1, r8, lsl #0x06
	RSB	r1, r1, r8, lsl #0x0E
	RSB	r1, r1, sl, asr #0x04
0:	MOV	ip, #0x7A00           @ c3_5[.15] -> ip
	ORR	ip, ip, #0x7D
	MUL	r8, ip, r2            @ bx =  c3_5*x[1] + s3_5*x[6] -> r8 [.15]
	ADD	sl, r2, r2, lsl #0x02 @ by = -s3_5*x[1] + c3_5*x[6] -> r2 [.15]
	ADD	sl, sl, sl, lsl #0x05
	ADD	sl, sl, r2, lsl #0x0A
	ADD	r2, r7, r7, lsl #0x02
	ADD	r2, r2, r2, lsl #0x05
	ADD	r2, r2, r7, lsl #0x0A
	ADD	r8, r8, r2, lsl #0x03
	MUL	r2, ip, r7
	SUB	r2, r2, sl, lsl #0x03
0:	MOV	ip, #0x7000 @ c5_5[.15] -> ip
	MOV	lr, #0x7800 @ s5_5[.16] -> lr
	ORR	ip, ip, #0xE3
	ORR	lr, lr, #0xAD
	MUL	r7, ip, r3  @ cx =  c5_5*x[2] + s5_5*x[5] -> r7 [.15]
	MUL	sl, lr, r3  @ cy =  s5_5*x[2] - c5_5*x[5] -> r3 [.15]
	MUL	r3, lr, r6
	ADD	r7, r7, r3, asr #0x01
	MUL	r3, ip, r6
	RSB	r3, r3, sl, asr #0x01
0:	MOV	ip, #0x3100 @ c7_5[.14] -> ip
	MOV	lr, #0x2800 @ s7_5[.14] -> lr
	ORR	ip, ip, #0x79
	ORR	lr, lr, #0x9A
	MUL	r6, ip, r4  @ dx =  c7_5*x[3] + s7_5*x[4] -> r6 [.14]
	MUL	sl, lr, r4  @ dy = -s7_5*x[3] + c7_5*x[4] -> r4 [.14]
	MLA	r6, lr, r5, r6
	MUL	r4, ip, r5
	SUB	r4, r4, sl
1:	MOV	r9, r9, asr #0x0E
	ADD	r9, r9, r6, asr #0x0E @ saxdx = ax+dx -> r9
	SUB	r6, r9, r6, asr #0x0D @ daxdx = ax-dx -> r6
	MOV	r8, r8, asr #0x0F
	ADD	r8, r8, r7, asr #0x0F @ sbxcx = bx+cx -> r8
	SUB	r7, r8, r7, asr #0x0E @ dbxcx = bx-cx -> r7
	MOV	r4, r4, asr #0x0E
	ADD	r4, r4, r1, asr #0x0E @ sdyay = dy+ay -> r4
	SUB	r1, r4, r1, asr #0x0D @ ddyay = dy-ay -> r1
	MOV	r3, r3, asr #0x0F
	ADD	r3, r3, r2, asr #0x0F @ scyby = cy+by -> r3
	SUB	r2, r3, r2, asr #0x0E @ dcyby = cy-by -> r2
2:	MOV	ip, #0x3B00 @ c1_3[.14] -> ip
	MOV	lr, #0x1800 @ s1_3[.14] -> lr
	ORR	ip, ip, #0x21
	ORR	lr, lr, #0x7E
	MUL	sl, ip, r6  @ tx = c1_3*daxdx + s1_3*dbxcx -> sl [.14]
	MUL	fp, lr, r6  @ ty = s1_3*daxdx - c1_3*dbxcx -> fp [.14]
	MLA	sl, lr, r7, sl
	MUL	r5, ip, r7
	SUB	fp, fp, r5
	MUL	r6, ip, r1  @ vx = c1_3*ddyay + s1_3*dcyby -> r6 [.14]
	MUL	r7, lr, r1  @ vy = s1_3*ddyay - c1_3*dcyby -> r7 [.14]
	MLA	r6, lr, r2, r6
	MUL	r5, ip, r2
	SUB	r7, r7, r5
	ADD	r1, r9, r8  @ sx = saxdx + sbxcx -> r1 = X0
	ADD	lr, r4, r3  @ ux = sdyay + scyby -> lr = X7
	SUB	r4, r4, r3  @ uy = sdyay - scyby -> r4
	SUB	r9, r9, r8  @ sy = saxdx - sbxcx -> r9
3:	MOV	sl, sl, asr #0x0E
	ADD	r3, sl, r7, asr #0x0E @ X2 = tx+vy -> r3
	SUB	r2, r3, r7, asr #0x0D @ X1 = tx-vy -> r2
	ADD	r5, r9, r4            @ X3 = (sy+uy)*sqrt1_2 -> r5 -> r6
	SUB	r8, r9, r4            @ X4 = (sy-uy)*sqrt1_2 -> r8 -> r9
	MOV	fp, fp, asr #0x0E
	ADD	ip, fp, r6, asr #0x0E @ X6 = ty+vx -> ip
	SUB	sl, ip, r6, asr #0x0D @ X5 = ty-vx -> sl
	MOV	r7, #0x2D00           @ sqrt1_2[.14] -> r7
	ORR	r7, r7, #0x41
	MUL	r6, r5, r7
	MUL	r9, r8, r7
	MOV	r6, r6, asr #0x0E
	MOV	r9, r9, asr #0x0E
	STMIA	r0, {r1,r2,r3,r6,r9,sl,ip,lr}
2:	LDMFD	sp!, {r4-fp,pc}

/**************************************/
.size   Fourier_DCT4, .-Fourier_DCT4
.global Fourier_DCT4
/**************************************/
/* EOF                                */
/**************************************/
