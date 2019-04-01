/**************************************/
.include "source/ulc/ulc_Specs.inc"
/**************************************/
.if BLOCK_OVERLAP == 2048
	.equ IMDCT_WS_BITS, 22
	.equ IMDCT_WS, (0x0C91) @ Sin[1.0*(1/Overlap)*Pi/2] * 2^22
	.equ IMDCT_C,  (0x8000) @ Cos[0.5*(1/Overlap)*Pi/2] * 2^15
	.equ IMDCT_S,  (0x000D) @ Sin[0.5*(1/Overlap)*Pi/2] * 2^15
.elseif BLOCK_OVERLAP == 1536
	.equ IMDCT_WS_BITS, 26
	.equ IMDCT_WS, (0x010C15)
	.equ IMDCT_C,  (0x8000)
	.equ IMDCT_S,  (0x0011)
.elseif BLOCK_OVERLAP == 768
	.equ IMDCT_WS_BITS, 25
	.equ IMDCT_WS, (0x010C15)
	.equ IMDCT_C,  (0x8000)
	.equ IMDCT_S,  (0x0022)
.elseif BLOCK_OVERLAP == 512
	.equ IMDCT_WS_BITS, 20
	.equ IMDCT_WS, (0x0C91)
	.equ IMDCT_C,  (0x8000)
	.equ IMDCT_S,  (0x0032)
.endif
/**************************************/
.section .iwram, "ax", %progbits
.balign 4
/**************************************/

ulc_BlockProcess:
	STMFD	sp!, {r4-fp,lr}
	LDR	r4, =ulc_State
0:	LDRH	r0, [r4, #0x00] @ --nBufProc?
	SUBS	r0, r0, #0x0100
	STRCSH	r0, [r4, #0x00]
	BCC	.LExit
0:	TST	r0, #0x01
	LDMIB	r4, {r0,r5,r6}  @ nBlkRem -> r0, &SoundFile -> r5, &NextByte -> r6
	LDR	r5, =ulc_OutputBuffer
	ADDEQ	r5, r5, #BLOCK_SIZE
	SUBS	r0, r0, #0x01   @ --nBlkRem?
	STRCS	r0, [r4, #0x04]
	BCC	.LNoBlocksRem
1:	LDRB	r7, [r6], #0x01 @ *Stream++ -> r7
	MOV	r8, #0x00       @ 0 -> r8
	MOV	r9, #0x0E40     @ "SUB r0, r8, r0, asr #0x1C", lower hword -> r9
	LDR	sl, =ulc_TransformBuffer

/**************************************/

@ r4: &State
@ r5: &OutputBuf | Chan<<31
@ r6: &NextByte
@ r7:  StreamByte
@ r8:  0
@ r9:  0x0E40
@ sl: &CoefDst
@ fp:  CoefRem

.macro NextNybble
	ADDS	r6, r6, #0x80000000 @ Odd nybble?
	MOVCC	r7, r7, lsr #0x04   @  Even: Move to next nybble
	LDRCSB	r7, [r6], #0x01     @  Odd:  Move to next byte
.endm

.LChannels_Loop:
	MOV	fp, #BLOCK_SIZE
	B	.LDecodeCoefs_Start

.LDecodeCoefs_ChangeQuant:
	NextNybble

.LDecodeCoefs_Start:
	AND	r0, r7, #0x0F
	NextNybble
	CMP	r0, #0x0F @ Stop? (8h,0h,Fh)
	BEQ	.LDecodeCoefs_Stop
0:	SUB	r0, r9, r0, lsl #0x07
	STRH	r0, .LDecodeCoefs_Normal_Shifter @ Form "SUB r0, r8, r0, asr #0x1C-n", lower hword (0E40h - n<<7)

.LDecodeCoefs_DecodeLoop:
	SUBS	r0, r8, r7, lsl #0x1C    @ -QCoef -> r0?
	BVS	.LDecodeCoefs_EscapeCode @ Escape code? (8h)

.LDecodeCoefs_Normal:
	NextNybble
.LDecodeCoefs_Normal_Shifter:
	SUB	r0, r8, r0, asr #0x00 @ Coef=QCoef<<Quant -> r0 (NOTE: Self-modified dequantization)
	STR	r0, [sl], #0x04       @ Coefs[n++] = Coef
	SUBS	fp, fp, #0x01         @ --CoefRem?
	BNE	.LDecodeCoefs_DecodeLoop
1:	B	.LDecodeCoefs_NoMoreCoefs

.LDecodeCoefs_EscapeCode:
	NextNybble
	ANDS	r0, r7, #0x0F @ Quantizer change? (8h,0h,Xh)
	BEQ	.LDecodeCoefs_ChangeQuant
	NextNybble
1:	CMP	r0, #0x0C     @ Short?
	BCC	2f
1:	ORR	r0, r0, r7, lsl #0x1C @ Long: 8h,Ch..Fh,Xh -> 26..152 zeros
	MOV	r0, r0, ror #0x1C
	ADD	r0, r0, #(26-2)/2 - (0xC<<4)
	NextNybble
2:	SUB	fp, fp, r0, lsl #0x01 @ CoefRem -= (zR-1)*2
	MOV	r1, #0x00
20:	STMIA	sl!, {r1,r8}
	SUBS	r0, r0, #0x01
	BCS	20b
3:	SUBS	fp, fp, #0x02 @ CoefRem -= 2? (because biased earlier)
	BNE	.LDecodeCoefs_DecodeLoop
	B	.LDecodeCoefs_NoMoreCoefs

.LDecodeCoefs_Stop:
	MOV	r0, #0x00
	MOV	r1, #0x00
	MOV	r2, #0x00
	MOVS	fp, fp, lsr #0x01
	STRCS	r0, [sl], #0x04
	MOVS	fp, fp, lsr #0x01
	STMCSIA	sl!, {r0-r1}
1:	STMNEIA	sl!, {r0-r2,r8}
	SUBNES	fp, fp, #0x01
	BNE	1b

.LDecodeCoefs_NoMoreCoefs:
	SUB	sl, sl, #BLOCK_SIZE*4 @ Rewind buffer

/**************************************/

@ r0:
@ r1:
@ r2:
@ r3:
@ r4: &State
@ r5: &OutputBuf | Chan<<31
@ r6: &NextByte
@ r7:  StreamByte
@ r8:  0
@ r9:  0x0E40
@ sl: &CoefBuf
@ fp:
@ ip:
@ lr:

.LDecodeCoefs_BlockUnpack:
	STMFD	sp!, {r4-r9}

.LDecodeCoefs_AntiPreEcho:
	MOV	lr, #BLOCK_SIZE
1:	LDMIA	sl, {r0-r7}
	ADD	r0, r0, r1
	SUB	r1, r0, r1, lsl #0x01
	ADD	r2, r2, r3
	SUB	r3, r2, r3, lsl #0x01
	ADD	r4, r4, r5
	SUB	r5, r4, r5, lsl #0x01
	ADD	r6, r6, r7
	SUB	r7, r6, r7, lsl #0x01
	STMIA	sl!, {r0-r7}
	LDMIA	sl, {r0-r7}
	ADD	r0, r0, r1
	SUB	r1, r0, r1, lsl #0x01
	ADD	r2, r2, r3
	SUB	r3, r2, r3, lsl #0x01
	ADD	r4, r4, r5
	SUB	r5, r4, r5, lsl #0x01
	ADD	r6, r6, r7
	SUB	r7, r6, r7, lsl #0x01
	STMIA	sl!, {r0-r7}
	SUBS	lr, lr, #0x10
	BNE	1b
2:	SUB	sl, sl, #BLOCK_SIZE*4 @ Rewind buffer

.LDecodeCoefs_IMDCT:
	MOV	r0, sl
	ADD	r1, sl, #BLOCK_SIZE*4
	MOV	r2, #BLOCK_SIZE
	BL	Fourier_DCT4
0:	MOV	r9, #0x7F
	ADD	sl, sl, #(BLOCK_SIZE/2)*4                 @ Tmp = TmpBuf+N/2 (forwards)
	LDR	fp, =ulc_LappingBuffer + (BLOCK_SIZE/2)*4 @ Lap = LapBuf+N/2 (backwards)
	LDR	ip, [sp, #0x04]     @ OutLo
.if ULC_STEREO
	TST	ip, #0x80000000
	ADDNE	ip, ip, #BLOCK_SIZE*2
	ADDNE	fp, fp, #(BLOCK_SIZE/2)*4
	BIC	ip, ip, #0x80000000
.endif
	ADD	lr, ip, #BLOCK_SIZE @ OutHi
.if BLOCK_OVERLAP < BLOCK_SIZE
0:	SUB	r9, r9, #(BLOCK_SIZE-BLOCK_OVERLAP) << 8
1:	LDMDB	fp!, {r0-r7}
	MOV	r0, r0, asr #0x08
	MOV	r1, r1, asr #0x08
	MOV	r2, r2, asr #0x08
	MOV	r3, r3, asr #0x08
	MOV	r4, r4, asr #0x08
	MOV	r5, r5, asr #0x08
	MOV	r6, r6, asr #0x08
	MOV	r7, r7, asr #0x08
	EORS	r0, r0, r0, lsl #0x18
	ADCMI	r0, r9, #0x00
	EORS	r1, r1, r1, lsl #0x18
	ADCMI	r1, r9, #0x00
	EORS	r2, r2, r2, lsl #0x18
	ADCMI	r2, r9, #0x00
	EORS	r3, r3, r3, lsl #0x18
	ADCMI	r3, r9, #0x00
	EORS	r4, r4, r4, lsl #0x18
	ADCMI	r4, r9, #0x00
	EORS	r5, r5, r5, lsl #0x18
	ADCMI	r5, r9, #0x00
	EORS	r6, r6, r6, lsl #0x18
	ADCMI	r6, r9, #0x00
	EORS	r7, r7, r7, lsl #0x18
	ADCMI	r7, r9, #0x00
	AND	r7, r7, #0xFF
	AND	r6, r6, #0xFF
	ORR	r7, r7, r6, lsl #0x08
	AND	r5, r5, #0xFF
	ORR	r7, r7, r5, lsl #0x10
	ORR	r7, r7, r4, lsl #0x18
	AND	r8, r3, #0xFF
	AND	r2, r2, #0xFF
	ORR	r8, r8, r2, lsl #0x08
	AND	r1, r1, #0xFF
	ORR	r8, r8, r1, lsl #0x10
	ORR	r8, r8, r0, lsl #0x18
	STMIA	ip!, {r7-r8}
	LDMIA	sl!, {r0-r7}
	MOV	r0, r0, asr #0x08
	MOV	r1, r1, asr #0x08
	MOV	r2, r2, asr #0x08
	MOV	r3, r3, asr #0x08
	MOV	r4, r4, asr #0x08
	MOV	r5, r5, asr #0x08
	MOV	r6, r6, asr #0x08
	MOV	r7, r7, asr #0x08
	EORS	r0, r0, r0, lsl #0x18
	ADCMI	r0, r9, #0x00
	EORS	r1, r1, r1, lsl #0x18
	ADCMI	r1, r9, #0x00
	EORS	r2, r2, r2, lsl #0x18
	ADCMI	r2, r9, #0x00
	EORS	r3, r3, r3, lsl #0x18
	ADCMI	r3, r9, #0x00
	EORS	r4, r4, r4, lsl #0x18
	ADCMI	r4, r9, #0x00
	EORS	r5, r5, r5, lsl #0x18
	ADCMI	r5, r9, #0x00
	EORS	r6, r6, r6, lsl #0x18
	ADCMI	r6, r9, #0x00
	EORS	r7, r7, r7, lsl #0x18
	ADCMI	r7, r9, #0x00
	AND	r3, r3, #0xFF
	AND	r2, r2, #0xFF
	ORR	r3, r3, r2, lsl #0x08
	AND	r1, r1, #0xFF
	ORR	r3, r3, r1, lsl #0x10
	ORR	r3, r3, r0, lsl #0x18
	AND	r2, r7, #0xFF
	AND	r6, r6, #0xFF
	ORR	r2, r2, r6, lsl #0x08
	AND	r5, r5, #0xFF
	ORR	r2, r2, r5, lsl #0x10
	ORR	r2, r2, r4, lsl #0x18
	STMDB	lr!, {r2-r3}
2:	ADDS	r9, r9, #0x10<<8
	BCC	1b
.endif
.if BLOCK_OVERLAP
0:	LDR	r6, =IMDCT_WS
	LDR	r7, =IMDCT_C
	LDR	r8, =IMDCT_S
1:	LDR	r0, [fp, #-0x04]!   @ a = *--Lap
	LDR	r1, [sl], #0x04     @ b = *Tmp++
	MUL	r4, r8, r0          @ *--OutHi = s*a + c*b
	MLA	r4, r7, r1, r4
	MOV	r4, r4, asr #0x0F+8
	EORS	r4, r4, r4, lsl #0x18
	ADCMI	r4, r9, #0x00
	AND	r4, r4, #0xFF
	ORR	r3, r4, r3, lsl #0x08
	MUL	r4, r7, r0          @ *OutLo++ = c*a - s*b
	MUL	r1, r8, r1
	SUB	r4, r4, r1
	MOV	r4, r4, asr #0x0F+8
	EORS	r4, r4, r4, lsl #0x18
	ADCMI	r4, r9, #0x00
	MOV	r2, r2, lsr #0x08
	ORR	r2, r2, r4, lsl #0x18
0:	MUL	r0, r6, r8          @ _c = wc*c - ws*s
	MUL	r1, r6, r7          @ _s = ws*c + wc*s
	SUB	r7, r7, r0, lsr #IMDCT_WS_BITS
	ADD	r8, r8, r1, lsr #IMDCT_WS_BITS
0:	ADDS	ip, ip, #0x40000000 @ Repeat until have 4 samples in each register
	BCC	1b
	STR	r2, [ip], #0x04
	STR	r3, [lr, #-0x04]!
2:	CMP	ip, lr
	BNE	1b
.endif
0:	MOV	r9, #BLOCK_SIZE/2
	SUB	sl, sl, #BLOCK_SIZE*4
1:	LDMIA	sl!, {r0-r7}
	STMIA	fp!, {r0-r7}
	LDMIA	sl!, {r0-r7}
	STMIA	fp!, {r0-r7}
	SUBS	r9, r9, #0x10
	BNE	1b
0:	LDMFD	sp!, {r4-r9}
	SUB	sl, sl, #(BLOCK_SIZE/2)*4 @ Rewind buffer

.LDecodeCoefs_NextChan:
.if ULC_STEREO
	ADDS	r5, r5, #0x80000000
	BCC	.LChannels_Loop
.endif

/**************************************/

@ r4: &State
@ r6: &NextByte

.LSaveState_Exit:
.if ULC_STEREO
	TST	r6, #0x80000000
	BIC	r6, r6, #0x80000000 @ Clear nybble mask (always move to next byte on next frame)
	SUBEQ	r6, r6, #0x01       @ We are always one nybble ahead, so rewind the byte if on first nybble
.else
	SUBS	r6, r6, #0x80000000
	ADDCC	r6, r6, #0x80000000-1
.endif
	STR	r6, [r4, #0x0C]

.LExit:
	LDMFD	sp!, {r4-fp,lr}
	BX	lr

/**************************************/

.LNoBlocksRem:
	MOV	r0, #0x00
	MOV	r1, #0x00
	MOV	r2, #0x00
	MOV	r3, #0x00
	MOV	r4, #0x00
	MOV	r6, #0x00
	MOV	r7, #0x00
	MOV	r8, #0x00
	MOV	r9, #BLOCK_SIZE
.if ULC_STEREO
	ADD	sl, r5, #BLOCK_SIZE*2
.endif
1:	STMIA	r5!, {r0-r4,r6-r8}
	STMIA	r5!, {r0-r4,r6-r8}
	STMIA	r5!, {r0-r4,r6-r8}
	STMIA	r5!, {r0-r4,r6-r8}
.if ULC_STEREO
	STMIA	sl!, {r0-r4,r6-r8}
	STMIA	sl!, {r0-r4,r6-r8}
	STMIA	sl!, {r0-r4,r6-r8}
	STMIA	sl!, {r0-r4,r6-r8}
.endif
	SUBS	r9, r9, #0x08*16
	BNE	1b
2:	B	.LExit

/**************************************/
.size   ulc_BlockProcess, .-ulc_BlockProcess
.global ulc_BlockProcess
/**************************************/
.section .bss
.balign 4
/**************************************/

ulc_TransformBuffer:
	.space 4*BLOCK_SIZE
.size ulc_TransformBuffer, .-ulc_TransformBuffer

ulc_TransformTemp:
	.space 4*BLOCK_SIZE
.size ulc_TransformTemp, .-ulc_TransformTemp

ulc_LappingBuffer:
	.space 4*(BLOCK_SIZE/2)
.if ULC_STEREO
	.space 4*(BLOCK_SIZE/2)
.endif
.size ulc_LappingBuffer, .-ulc_LappingBuffer

/**************************************/
/* EOF                                */
/**************************************/
