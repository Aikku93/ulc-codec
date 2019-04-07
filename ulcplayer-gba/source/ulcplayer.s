/**************************************/
.include "source/ulc/ulc_Specs.inc"
/**************************************/
.text
.balign 2
/**************************************/
@ BgDesign to BG3
@ Glyphs to BG2
@ GraphR to BG1
@ GraphL to BG0
/**************************************/
.equ GRAPH_W,  112 @ Pixels
.equ GRAPH_H,   64 @ Pixels (NOTE: lower half reflected from top)
.equ GRAPH_TAIL, 8 @ Pixels (NOTE: fixed 8px in code, just for reference)
.equ GRAPHL_X,  56 @ Pixels
.equ GRAPHL_Y,  48 @ Pixels
.equ GRAPHR_X,  56 @ Pixels
.equ GRAPHR_Y,  48 @ Pixels
.equ ARTIST_X,  48 @ Pixels (NOTE: nominal, whole row is used)
.equ ARTIST_Y,  16 @ Pixels
.equ ARTIST_W, 144 @ Pixels (NOTE: nominal, whole row is used)
.equ TITLE_X,   48 @ Pixels (NOTE: nominal, whole row is used)
.equ TITLE_Y,   32 @ Pixels
.equ TITLE_W,  144 @ Pixels (NOTE: nominal, whole row is used)
.equ SPEAKER_LT_X,   8
.equ SPEAKER_LT_Y,  24
.equ SPEAKER_LB_X,   8
.equ SPEAKER_LB_Y,  72
.equ SPEAKER_RT_X, 200
.equ SPEAKER_RT_Y,  24
.equ SPEAKER_RB_X, 200
.equ SPEAKER_RB_Y,  72
.equ GRAPH_SMPSTRIDE_RCP, 627 @ (1<<22) / (VBlankRate * GRAPH_W)
/**************************************/
.equ BGDESIGN_NTILES, 211
.equ BGDESIGN_TILEMAP, 31
.equ BGDESIGN_NPAL16,   5

.equ GLYPHS_NTILES,  95
.equ GLYPHS_TILEMAP, 30
.equ GLYPHS_PAL16,    4
.equ GLYPHS_TILEOFS, BGDESIGN_NTILES
.equ GLYPHS_TILEADR, (0x06000000 + GLYPHS_TILEMAP*0x0800)

.equ GRAPH_NTILES, ((GRAPH_W+GRAPH_TAIL*2)*GRAPH_H/2) / (8*8)
.equ GRAPHL_TILEMAP, 29
.equ GRAPHL_TILEOFS, (GLYPHS_TILEOFS + GLYPHS_NTILES)
.equ GRAPHL_TILEADR, (0x06000000 + GRAPHL_TILEOFS*32)
.equ GRAPHL_PAL16,    2
.equ GRAPHR_TILEMAP, 28
.equ GRAPHR_TILEOFS, (GRAPHL_TILEOFS + GRAPH_NTILES)
.equ GRAPHR_TILEADR, (0x06000000 + GRAPHR_TILEOFS*32)
.equ GRAPHR_PAL16,    3

.equ SPEAKERBASS_NTILES, 128
/**************************************/

.thumb
.thumb_func
main:
.Lmain_LoadDesign:
0:	LDR	r0, =0x06000000
	LDR	r1, =BgDesign_Gfx
	LDR	r2, =(BGDESIGN_NTILES + GLYPHS_NTILES) * 32
	BL	.Lmain_Copy32
0:	LDR	r0, =0x06000000 + BGDESIGN_TILEMAP*0x0800
	LDR	r1, =BgDesign_Map
	LDR	r2, =0x0500
	BL	.Lmain_Copy32
0:	LDR	r0, =0x05000000
	LDR	r1, =BgDesign_Pal
	LDR	r2, =BGDESIGN_NPAL16 * 16*2
	BL	.Lmain_Copy32
0:	LDR	r0, =0x06010000
	LDR	r1, =BgDesignSpeakerBass_Gfx
	LDR	r2, =SPEAKERBASS_NTILES * 32
	BL	.Lmain_Copy32
0:	LDR	r0, =0x05000200
	LDR	r1, =BgDesignSpeakerBass_Pal
	LDR	r2, =16*2
	BL	.Lmain_Copy32
0:	LDR	r0, =0x06000000 + GRAPHL_TILEMAP*0x0800
	LDR	r1, =0
	LDR	r2, =0x0500
	BL	.Lmain_Set32
0:	LDR	r0, =0x06000000 + GRAPHR_TILEMAP*0x0800
	LDR	r1, =0
	LDR	r2, =0x0500
	BL	.Lmain_Set32
0:	BL	.Lmain_InitGraphTiles
0:	LDR	r0, =0x07000000
	LDR	r2, =SPEAKER_LT_Y|(1<<10) | (SPEAKER_LT_X | 2<<14)<<16
	LDR	r3, =0x0000
	LDR	r4, =SPEAKER_LB_Y|(1<<10) | (SPEAKER_LB_X | 2<<14)<<16
	LDR	r5, =0x0000
	STMIA	r0!, {r2-r5}
	LDR	r2, =SPEAKER_RT_Y|(1<<10) | (SPEAKER_RT_X | 2<<14)<<16
	LDR	r3, =0x0000
	LDR	r4, =SPEAKER_RB_Y|(1<<10) | (SPEAKER_RB_X | 2<<14)<<16
	LDR	r5, =0x0000
	STMIA	r0!, {r2-r5}
	MOV	r2, #0x80-4-1 @ Clear remaining OAMs
	LSL	r2, #0x03
	MOV	r1, #0x02
	LSL	r1, #0x08
1:	STRH	r1, [r0, r2]
	SUB	r2, #0x08
	BCS	1b

.LMain:
	LDR	r1, =_IRQTable
	LDR	r0, =UpdateGfx
	STR	r0, [r1, #0x04*0]
	LDR	r1, =0x04000004
	LDR	r0, =1<<3
	STRH	r0, [r1]
	LDR	r1, =0x04000200
	LDRH	r0, [r1]
	ADD	r0, #0x01
	STRH	r0, [r1]
0:	LDR	r0, =SoundFile
	BL	ulc_Init
1:	BL	ulc_BlockProcess
	MOV	r0, #0x00 @ IntrWait (return if already set)
	MVN	r1, r0    @ Any interrupt
	SWI	0x04
	B	1b
.Lbxr5:	BX	r5
.Lbxr6:	BX	r6

.Lmain_Copy32:
1:	LDMIA	r1!, {r3}
	STMIA	r0!, {r3}
	SUB	r2, #0x04
	BNE	1b
2:	BX	lr

.Lmain_Set32:
1:	STMIA	r0!, {r1}
	SUB	r2, #0x04
	BNE	1b
2:	BX	lr

.balign 4
.Lmain_InitGraphTiles:
	LDR	r0, =0x06000000 + GRAPHL_TILEMAP*0x0800
	LDR	r1, =0x06000000 + GRAPHR_TILEMAP*0x0800
	LDR	r2, =0x06000000 + GRAPHL_TILEMAP*0x0800 + ((GRAPH_H/8-1)*32)*2
	LDR	r3, =0x06000000 + GRAPHR_TILEMAP*0x0800 + ((GRAPH_H/8-1)*32)*2
	LDR	r4, =GRAPHL_TILEOFS | GRAPHL_PAL16<<12
	LDR	r5, =GRAPHR_TILEOFS | GRAPHR_PAL16<<12
	LDR	r6, =GRAPH_H/2/8 @ Reflected
	BX	pc
	NOP
.arm
1:	SUB	r6, r6, #((GRAPH_W+GRAPH_TAIL*2)/8)<<16
0:	STRH	r4, [r0], #0x02
	STRH	r5, [r1], #0x02
	ORR	ip, r4, #0x01<<11
	STRH	ip, [r2], #0x02
	ORR	ip, r5, #0x01<<11
	STRH	ip, [r3], #0x02
	ADD	r4, r4, #(GRAPH_H/2)/8
	ADD	r5, r5, #(GRAPH_H/2)/8
	ADDS	r6, r6, #0x01<<16
	BCC	0b
0:	ADD	r0, r0, #(32-(GRAPH_W+GRAPH_TAIL*2)/8)*2 @ Next row
	ADD	r1, r1, #(32-(GRAPH_W+GRAPH_TAIL*2)/8)*2
	SUB	r2, r2, #(32+(GRAPH_W+GRAPH_TAIL*2)/8)*2 @ Previous row
	SUB	r3, r3, #(32+(GRAPH_W+GRAPH_TAIL*2)/8)*2
	SUB	r4, r4, #(GRAPH_H/2/8)*((GRAPH_W+GRAPH_TAIL*2)/8)-1
	SUB	r5, r5, #(GRAPH_H/2/8)*((GRAPH_W+GRAPH_TAIL*2)/8)-1
	SUBS	r6, r6, #0x01
	BNE	1b
2:	BX	lr

/**************************************/
.size   main, .-main
.global main
/**************************************/
.section .iwram, "ax", %progbits
.balign 4
/**************************************/

.arm
UpdateGfx:
	STMFD	sp!, {r4-fp,lr}
	LDR	r4, =ulc_State
0:	MOV	r0, #0x04000000
	MOV	r1, #0x1F00
	ORR	r1, r1, #0x40
	STRH	r1, [r0]
	MOV	r1, #GRAPHL_TILEMAP<<8
	ORR	r1, r1, #(GRAPHR_TILEMAP<<8)<<16
	STR	r1, [r0, #0x08]
	MOV	r1, #GLYPHS_TILEMAP<<8
	ORR	r1, r1, #(BGDESIGN_TILEMAP<<8)<<16
	STR	r1, [r0, #0x0C]
	LDR	r1, =((-GRAPHL_X)&0xFFFF) | ((-GRAPHL_Y)&0xFFFF)<<16
	STR	r1, [r0, #0x10]
	LDR	r1, =((-GRAPHR_X)&0xFFFF) | ((-GRAPHR_Y)&0xFFFF)<<16
	STR	r1, [r0, #0x14]
	MOV	r1, #0x00
	STR	r1, [r0, #0x18]
	STR	r1, [r0, #0x1C]
	LDR	r1, =0x10102F53
	STR	r1, [r0, #0x50]!       @ Layer BG0,BG1,OBJ over BG0,BG1,BG2,BG3, additive blend
0:	LDRH	r5, [r0, #0x0104-0x50] @ Get SmpPos from timer -> r5
	LDRB	ip, [r4, #0x02]
	SUB	r5, r5, #0x010000+BLOCK_SIZE @ Adjust for double buffer
	ADD	r5, r5, ip, lsl #BLOCK_SIZE_LOG2

.LRedraw_Clear:
1:	LDR	r0, =GLYPHS_TILEADR + ((ARTIST_Y/8)*32)*2
	MOV	r1, #0x00
	MOV	r2, #0x20*2 @ Clear whole row
	BL	.LRedraw_Set32
1:	LDR	r0, =GLYPHS_TILEADR + ((TITLE_Y/8)*32)*2
	MOV	r1, #0x00
	MOV	r2, #0x20*2 @ Clear whole row
	BL	.LRedraw_Set32
1:	LDR	r0, =GRAPHL_TILEADR
	MOV	r1, #0x00
	MOV	r2, #GRAPH_NTILES * 32 * 2 @ Clear L+R
	BL	.LRedraw_Set32

.LRedraw_DrawTitle:
1:	LDR	r0, =GLYPHS_TILEADR + ((ARTIST_Y/8)*32)*2
	LDR	r1, =SoundFile_Artist
	BL	.LRedraw_DrawString
1:	LDR	r0, =GLYPHS_TILEADR + ((TITLE_Y/8)*32)*2
	LDR	r1, =SoundFile_Title
	BL	.LRedraw_DrawString

.LRedraw_GetSamples:
	LDR	r0, =.LRedraw_GraphDataL
.if ULC_STEREO
	MOV	r1, #(BLOCK_SIZE*2) >> 8             @ Distance to right channel
	SUB	r1, r1, #GRAPH_W<<24
.else
	MOV	r1, #(-GRAPH_W)<<24
.endif
	LDR	r3, =ulc_OutputBuffer + BLOCK_SIZE*2 @ End -> r3
	ADD	r2, r3, r5                           @ Src -> r2
	LDR	r4, [r4, #0x0C]
	LDR	r4, [r4, #0x0C]
	LDR	r5, =GRAPH_SMPSTRIDE_RCP
	MUL	r4, r5, r4 @ Step[.22fxp]
	MOV	r5, #0x00  @ PosMu (not important to track accurately across frames)
	LDR	r6, .LRedraw_LowPassL
	LDR	r7, .LRedraw_LowPassR
	MOV	r8, #0x00 @ PowL=0
	MOV	r9, #0x00 @ PowR=0
1:	ADD	r5, r5, r4              @ [PosMu += Step]
0:	LDRB	ip, [r2, r1, lsl #0x08] @ Abs[xR] -> ip
	MOVS	ip, ip, lsl #0x18
	RSB	lr, r7, ip, asr #0x18-5 @ LP_R = LP_R + (xR - LP_R)*1/32 (NOTE: 8.5fxp)
	ADD	r7, r7, lr, asr #0x05
	MLA	r9, r7, r7, r9          @ PowR += LP_R (16.10fxp)
	RSBMI	ip, ip, #0x00
	LDRB	lr, [r0, #GRAPH_W]      @ Combine with old (nicer effect)
	RSB	ip, lr, ip, lsr #0x18
	ADD	ip, lr, ip, asr #0x03
	STRB	ip, [r0, #GRAPH_W]
0:	LDRB	ip, [r2], r5, lsr #0x16 @ Abs[xL] -> ip, update position
	MOVS	ip, ip, lsl #0x18
	RSB	lr, r6, ip, asr #0x18-5 @ LP_L = LP_L + (xL - LP_L)*1/32
	ADD	r6, r6, lr, asr #0x05
	MLA	r8, r6, r6, r8          @ PowL += LP_L
	RSBMI	ip, ip, #0x00
	LDRB	lr, [r0]                @ Combine with old
	RSB	ip, lr, ip, lsr #0x18
	ADD	ip, lr, ip, asr #0x03
	STRB	ip, [r0], #0x01
0:	BIC	r5, r5, #0xFF<<22       @ Clear integer part (step is less than 8 at sane rates, so clearing only up to 255 is more than fine)
	CMP	r2, r3                  @ Wrap
	SUBCS	r2, r2, #BLOCK_SIZE*2
	ADDS	r1, r1, #0x01<<24
	BCC	1b
2:	STR	r6, .LRedraw_LowPassL
	STR	r7, .LRedraw_LowPassR

.LRedraw_DrawSpeakers:
	LDR	r6, .LRedraw_SpeakerLowPassL
	LDR	r7, .LRedraw_SpeakerLowPassL
	SUBS	r8, r8, r6
	ADDCC	r8, r6, r8, asr #0x02   @ Dampen decay  (heavy)
	ADDCS	r8, r6, r8, asr #0x01   @ Dampen attack (slight)
	SUBS	r9, r9, r7
	ADDCC	r9, r7, r9, asr #0x02
	ADDCS	r9, r7, r9, asr #0x01
	STR	r8, .LRedraw_SpeakerLowPassL
	STR	r9, .LRedraw_SpeakerLowPassL
	MVN	ip, #0x0F
	AND	r8, ip, r8, lsr #0x18-4 @ 16x 8x8 tiles per 32x32 area (arbitrary scaling)
	AND	r9, ip, r9, lsr #0x18-4
	CMP	r8, #0x07<<4            @ Clip animation frames
	MOVHI	r8, #0x07<<4
	CMP	r9, #0x07<<4
	MOVHI	r9, #0x07<<4
0:	MOV	ip, #0x07000000
	STRH	r8, [ip, #0x08*0+0x04] @ L-T (tile in Attr2 bit 0..9)
	STRH	r8, [ip, #0x08*1+0x04] @ L-B
	STRH	r9, [ip, #0x08*2+0x04] @ R-T
	STRH	r9, [ip, #0x08*3+0x04] @ R-B

.LRedraw_DrawGraphs:
	LDR	r2, =GRAPHL_TILEADR + (GRAPH_H/2-1)*4 + (GRAPH_H/2)*(GRAPH_TAIL/8)*4
1:	MOV	r5, #0x01            @ PixelStep
	LDRB	r6, [r0, #-GRAPH_W]! @ L-side (L)
	ADR	lr, 0f
0:	MOV	r5, r5, ror #0x04
	SUB	r7, r2, #(GRAPH_H/2)*(GRAPH_TAIL/8)*4
	MOVS	r6, r6, lsr #0x01
	BNE	.LRedraw_DrawGraphBar
2:	MOV	r5, #0x01
	LDRB	r6, [r0, #GRAPH_W]   @ L-side (R)
	ADR	lr, 0f
0:	MOV	r5, r5, ror #0x04
	SUB	r7, r2, #(GRAPH_H/2)*(GRAPH_TAIL/8)*4 - (GRAPHR_TILEADR-GRAPHL_TILEADR)
	MOVS	r6, r6, lsr #0x01
	BNE	.LRedraw_DrawGraphBar
1:	MOV	r5, #0x10000000
	LDRB	r6, [r0, #GRAPH_W-1] @ R-side (L)
	ADR	lr, 0f
0:	MOV	r5, r5, ror #0x1C
	ADD	r7, r2, #(GRAPH_H/2)*(GRAPH_W/8)*4
	MOVS	r6, r6, lsr #0x01
	BNE	.LRedraw_DrawGraphBar
2:	MOV	r5, #0x10000000
	LDRB	r6, [r0, #GRAPH_W*2-1] @ R-side (R)
	ADR	lr, 0f
0:	MOV	r5, r5, ror #0x1C
	ADD	r7, r2, #(GRAPH_H/2)*(GRAPH_W/8)*4 + (GRAPHR_TILEADR-GRAPHL_TILEADR)
	MOVS	r6, r6, lsr #0x01
	BNE	.LRedraw_DrawGraphBar
3:	LDR	r4, =GRAPH_W
	MOV	r5, #0x01 @ PixelStep (will be rotated every pixel)
1:	LDRB	r6, [r0, #GRAPH_W] @ ValueR -> r6
	ADD	r7, r2, #GRAPHR_TILEADR-GRAPHL_TILEADR
	MOVS	r6, r6, lsr #0x01  @ Rescale
	BLNE	.LRedraw_DrawGraphBar
2:	LDRB	r6, [r0], #0x01    @ ValueL -> r6
	MOV	r7, r2
	MOVS	r6, r6, lsr #0x01  @ Rescale
	BLNE	.LRedraw_DrawGraphBar
3:	MOV	r5, r5, ror #0x1C  @ Rotate PixelStep
	CMP	r5, #0x01          @ Wrapped around? Move to next tile
	ADDEQ	r2, r2, #(GRAPH_H/2)*4
	SUBS	r4, r4, #0x01
	BHI	1b

.LRedraw_Exit:
	LDMFD	sp!, {r4-fp,lr}
	BX	lr

@ NOTE: Copying 128 bits at a time
.LRedraw_Set32:
	STR	lr, [sp, #-0x04]!
	MOV	r3, r1
	MOV	ip, r1
	MOV	lr, r1
1:	STMIA	r0!, {r1,r3,ip,lr}
	SUBS	r2, r2, #0x10
	BNE	1b
2:	LDR	pc, [sp], #0x04

@ r0: &Dst
@ r1: &String
.LRedraw_DrawString:
	MOV	r2, r1
1:	LDRB	r3, [r2], #0x01
	CMP	r3, #0x01
	BCS	1b
2:	SBC	r2, r2, r1
	RSB	r2, r2, #0x1E
	MOV	r2, r2, lsr #0x01
	ADD	r0, r0, r2, lsl #0x01
	LDR	r2, =GLYPHS_TILEOFS + (GLYPHS_PAL16<<12)
0:	LDRB	r3, [r1], #0x01
1:	SUBS	r3, r3, #0x21
	ADDCS	r3, r3, r2
	MOVCC	r3, #0x00
	STRH	r3, [r0], #0x02
	LDRB	r3, [r1], #0x01
	CMP	r3, #0x00
	BNE	1b
2:	BX	lr

@ r5:  PixelStep
@ r6:  Value (guaranteed not zero, will NOT be destroyed (only clipped))
@ r7: &Target
@ r8:  [Scratch]
@ r9:  [Scratch]
@ ip:  [Scratch]
.LRedraw_DrawGraphBar:
	CMP	r6, #GRAPH_H/2
	MOVHI	r6, #GRAPH_H/2
	RSB	r8, r5, r5, lsl #0x04 @ Start at Fh (full opacity/brightness)
0:	SUBS	r9, r6, #0x0F         @ Handle small bars
	MULCC	r8, r5, r6
	BLS	2f
1:	LDR	ip, [r7]
	ORR	ip, ip, r8
	STR	ip, [r7], #-0x04
	SUBS	r9, r9, #0x01
	BNE	1b
2:	LDR	r9, [r7]
	ORR	r9, r9, r8
	STR	r9, [r7], #-0x04
	SUBS	r8, r8, r5
	BHI	2b
4:	BX	lr

.LRedraw_LowPassL: .word 0
.LRedraw_LowPassR: .word 0
.LRedraw_SpeakerLowPassL: .word 0
.LRedraw_SpeakerLowPassR: .word 0
.LRedraw_GraphDataL: .space GRAPH_W
.LRedraw_GraphDataR: .space GRAPH_W

/**************************************/
.size UpdateGfx, .-UpdateGfx
/**************************************/

.section .rodata
.balign 4

SoundFile:
	.incbin "source/res/SoundData.ulc"
.size SoundFile, .-SoundFile

SoundFile_Artist:
	.asciz "Artist"
.size SoundFile_Artist, .-SoundFile_Artist

SoundFile_Title:
	.asciz "Title..."
.size SoundFile_Title, .-SoundFile_Title

/**************************************/
/* EOF                                */
/**************************************/
