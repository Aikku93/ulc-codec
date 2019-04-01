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
.equ GRAPHL_X,  64 @ Pixels
.equ GRAPHL_Y,  48 @ Pixels
.equ GRAPHR_X,  64 @ Pixels
.equ GRAPHR_Y,  48 @ Pixels
.equ ARTIST_X,  64 @ Pixels
.equ ARTIST_Y,  16 @ Pixels
.equ ARTIST_W, 112 @ Pixels
.equ TITLE_X,   64 @ Pixels
.equ TITLE_Y,   32 @ Pixels
.equ TITLE_W,  112 @ Pixels
.equ SPEAKER_LT_X,   8
.equ SPEAKER_LT_Y,  24
.equ SPEAKER_LB_X,   8
.equ SPEAKER_LB_Y,  72
.equ SPEAKER_RT_X, 200
.equ SPEAKER_RT_Y,  24
.equ SPEAKER_RB_X, 200
.equ SPEAKER_RB_Y,  72
.equ VBLANK_READY_ADR, 0x0203FFFC
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

.equ GRAPH_NTILES, (GRAPH_W*GRAPH_H) / (8*8)
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
	LDR	r6, =GRAPH_H/8/2 @ Reflected
	BX	pc
	NOP
.arm
1:	SUB	r6, r6, #(GRAPH_W/8)<<16
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
0:	ADD	r0, r0, #(32-GRAPH_W/8)*2 @ Next row
	ADD	r1, r1, #(32-GRAPH_W/8)*2
	SUB	r2, r2, #(32+GRAPH_W/8)*2 @ Previous row
	SUB	r3, r3, #(32+GRAPH_W/8)*2
	SUB	r4, r4, #(GRAPH_H/2/8)*(GRAPH_W/8)-1
	SUB	r5, r5, #(GRAPH_H/2/8)*(GRAPH_W/8)-1
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
	LDR	r0, =0x04000000
	LDR	r1, =0x1F40
	STRH	r1, [r0]
	LDR	r1, =GRAPHL_TILEMAP<<8 | (GRAPHR_TILEMAP<<8)<<16
	STR	r1, [r0, #0x08]
	LDR	r1, =GLYPHS_TILEMAP<<8 | (BGDESIGN_TILEMAP<<8)<<16
	STR	r1, [r0, #0x0C]
	LDR	r1, =((-GRAPHL_X)&0xFFFF) | ((-GRAPHL_Y)&0xFFFF)<<16
	STR	r1, [r0, #0x10]
	LDR	r1, =((-GRAPHR_X)&0xFFFF) | ((-GRAPHR_Y)&0xFFFF)<<16
	STR	r1, [r0, #0x14]
	LDR	r1, =0
	STR	r1, [r0, #0x18]
	STR	r1, [r0, #0x1C]
	LDR	r1, =0x10102F53
	STR	r1, [r0, #0x50] @ Layer BG0,BG1,OBJ over BG0,BG1,BG2,BG3, additive blend
0:	LDR	r0, .LRedraw_BufOfs
	LDR	r1, .LRedraw_PosMu
	LDR	r2, =ulc_State
	LDR	r2, [r2, #0x08]
	LDR	r2, [r2, #0x0C]
	LDR	r3, =GRAPH_SMPSTRIDE_RCP
	MUL	r3, r2, r3
	LDR	r2, =GRAPH_W
	UMULL	r3, ip, r2, r3 @ GRAPH_W * Rate/Stride + Mu
	ADDS	r3, r3, r1
	ADC	ip, ip, #0x00
	MOV	r2, r3, lsr #0x16
	BIC	r3, r3, r2, lsl #0x16
	STR	r3, .LRedraw_PosMu
	ORR	r2, r2, ip, lsl #0x20-22
	ADDS	r0, r0, r2
	SUBCS	r0, r0, #BLOCK_SIZE*2
	STR	r0, .LRedraw_BufOfs
0:	STMFD	sp!, {r4-fp,lr}

.LRedraw_Clear:
1:	LDR	r0, =GLYPHS_TILEADR + ((ARTIST_Y/8)*32 + ARTIST_X/8)*2
	LDR	r1, =0
	LDR	r2, =ARTIST_W*2
	BL	.LRedraw_Set32
1:	LDR	r0, =GLYPHS_TILEADR + ((TITLE_Y/8)*32 + TITLE_X/8)*2
	LDR	r1, =0
	LDR	r2, =TITLE_W*2
	BL	.LRedraw_Set32
1:	LDR	r0, =GRAPHL_TILEADR
	LDR	r1, =0
	LDR	r2, =GRAPH_NTILES * 32
	BL	.LRedraw_Set32
1:	LDR	r0, =GRAPHR_TILEADR
	LDR	r1, =0
	LDR	r2, =GRAPH_NTILES * 32
	BL	.LRedraw_Set32

.LRedraw_DrawTitle:
1:	LDR	r0, =GLYPHS_TILEADR + ((ARTIST_Y/8)*32 + ARTIST_X/8)*2
	LDR	r1, =SoundFile_Artist
	LDR	r2, =ARTIST_W
	BL	.LRedraw_DrawString
1:	LDR	r0, =GLYPHS_TILEADR + ((TITLE_Y/8)*32 + TITLE_X/8)*2
	LDR	r1, =SoundFile_Title
	LDR	r2, =TITLE_W
	BL	.LRedraw_DrawString

.LRedraw_GetSamples:
	LDR	r0, =.LRedraw_GraphDataL
	LDR	r1, =.LRedraw_GraphDataR
	LDR	r2, =ulc_OutputBuffer + BLOCK_SIZE*2
	LDR	r3, =ulc_OutputBuffer + BLOCK_SIZE*2 + BLOCK_SIZE*2
	LDR	r4, =ulc_State
	LDR	r4, [r4, #0x08]
	LDR	r4, [r4, #0x0C]
	LDR	ip, =GRAPH_SMPSTRIDE_RCP
	MUL	r4, ip, r4 @ Step[.22fxp]
	LDR	r5, .LRedraw_PosMu
	LDR	r6, .LRedraw_LowPassL
	LDR	r7, .LRedraw_LowPassR
	MOV	r8, #0x00 @ PowL = 0 (for speaker effect)
	MOV	r9, #0x00 @ PowR = 0
	LDR	sl, .LRedraw_BufOfs
	LDR	fp, =GRAPH_W
1:	LDRB	ip, [r2, sl]            @ Abs[xL] -> ip
	CMP	ip, #0x80
	RSBCS	ip, ip, #0x0100
	RSB	r6, r6, r6, lsr #0x04   @ LP_L = LP_L*15/16 + xL*1/16
	RSB	r6, r6, ip, lsl #0x10-4 @ NOTE: LP_L in 8.16
	ADD	r8, r8, r6
	LDRB	lr, [r0]
	RSB	lr, lr, lr, lsr #0x02   @ Combine with old (nicer effect)
	RSB	ip, lr, ip, lsr #0x02
	STRB	ip, [r0], #0x01
	LDRB	ip, [r3, sl]            @ Abs[xR] -> ip
	CMP	ip, #0x80
	RSBCS	ip, ip, #0x0100
	RSB	r7, r7, r7, lsr #0x04   @ LP_R = LP_R*15/16 + xR*1/16
	RSB	r7, r7, ip, lsl #0x10-4 @ NOTE: LP_R in 8.16
	ADD	r9, r9, r7
	LDRB	lr, [r1]
	RSB	lr, lr, lr, lsr #0x02   @ Combine with old
	RSB	ip, lr, ip, lsr #0x02
	STRB	ip, [r1], #0x01
	ADD	r5, r5, r4            @ PosMu += Step
	ADDS	sl, sl, r5, lsr #0x16 @ Update position, wrap, clear integer part of PosMu
	BIC	r5, r5, #0xFF<<22     @ (step is less than 8 at sane rates, so clearing only up to 255 is more than fine)
	SUBCS	sl, sl, #BLOCK_SIZE*2
	SUBS	fp, fp, #0x01
	BNE	1b
2:	@STR	sl, .LRedraw_BufOfs @ Updated in VBlank or can go out of sync from long decode times
	@STR	r5, .LRedraw_PosMu
	STR	r6, .LRedraw_LowPassL
	STR	r7, .LRedraw_LowPassR

.LRedraw_DrawSpeakers:
	LDR	r5, =0x07000000
	LDR	r2, .LRedraw_PowerOldL
	LDR	r3, .LRedraw_PowerOldR
	RSB	r8, r8, r2
	RSB	r9, r9, r3
	SUB	r2, r2, r8, asr #0x01 @ Dampen attack/decay
	SUB	r3, r3, r9, asr #0x01
	STR	r2, .LRedraw_PowerOldL
	STR	r3, .LRedraw_PowerOldR
	UMULL	r8, r9, r2, r2 @ Arbitrary scaling
	UMULL	r2, ip, r3, r3
	MOV	r2, r9, lsr #0x17
	MOV	r3, ip, lsr #0x17
	CMP	r2, #0x03
	MOVHI	r2, #0x03
	CMP	r3, #0x03
	MOVHI	r3, #0x03
	AND	ip, r8, #0x03 @ Using r8 as source of entropy
	CMP	ip, #0x02
	MOVHI	ip, #0x02
	CMP	r2, #0x02
	ADDCS	r2, r2, ip, lsl #0x01
	CMP	r3, #0x02
	ADDCS	r3, r3, ip, lsl #0x01
	MOV	r2, r2, lsl #0x04 @ 16x 8x8 tiles per 32x32 area
	MOV	r3, r3, lsl #0x04 @ 16x 8x8 tiles per 32x32 area
	STRH	r2, [r5, #0x08*0+0x04] @ L-T (tile in Attr2 bit 0..9)
	STRH	r2, [r5, #0x08*1+0x04] @ L-B
	STRH	r3, [r5, #0x08*2+0x04] @ R-T
	STRH	r3, [r5, #0x08*3+0x04] @ R-B

.LRedraw_DrawGraphs:
	SUB	r0, r0, #GRAPH_W
	SUB	r1, r1, #GRAPH_W
	LDR	r2, =GRAPHL_TILEADR + (GRAPH_H/2-1)*4
	LDR	r3, =GRAPHR_TILEADR + (GRAPH_H/2-1)*4
	LDR	r4, =GRAPH_W
	MOV	r5, #0x01         @ PixelStep (will be rotated every pixel)
1:	LDRB	r6, [r0], #0x01   @ Value -> r6
	MOV	r7, r2
	MOVS	r6, r6, lsr #0x01 @ Rescale
	BLNE	.LRedraw_DrawGraphBar
2:	LDRB	r6, [r1], #0x01   @ Value -> r6
	MOV	r7, r3
	MOVS	r6, r6, lsr #0x01 @ Rescale
	BLNE	.LRedraw_DrawGraphBar
3:	MOV	r5, r5, ror #0x1C @ Rotate PixelStep
	CMP	r5, #0x01         @ Wrapped around? Move to next tile
	ADDEQ	r2, r2, #(GRAPH_H/2)*4
	ADDEQ	r3, r3, #(GRAPH_H/2)*4
	SUBS	r4, r4, #0x01
	BHI	1b

.LRedraw_Exit:
	LDMFD	sp!, {r4-fp,lr}
	BX	lr

.LRedraw_Set32:
1:	STR	r1, [r0], #0x04
	SUBS	r2, r2, #0x04
	BHI	1b
2:	BX	lr

@ r0: &Dst
@ r1: &String
@ r2:  WidthPx
.LRedraw_DrawString:
	MOV	r3, r1
1:	LDRB	ip, [r3], #0x01
	CMP	ip, #0x01
	BCS	1b
2:	SBC	r3, r3, r1
	RSB	r3, r3, r2, lsr #0x03
	MOV	r3, r3, lsr #0x01
	ADD	r0, r0, r3, lsl #0x01
	LDR	r3, =GLYPHS_TILEOFS + (GLYPHS_PAL16<<12)
0:	LDRB	ip, [r1], #0x01
1:	SUBS	ip, ip, #0x21
	ADDCS	ip, ip, r3
	MOVCC	ip, #0x00
	STRH	ip, [r0], #0x02
	LDRB	ip, [r1], #0x01
	CMP	ip, #0x00
	BNE	1b
2:	BX	lr

@ r5:  PixelStep
@ r6:  Value (guaranteed not zero, will be destroyed)
@ r7: &Target
@ r8:  [Scratch]
@ r9:  [Scratch]
.LRedraw_DrawGraphBar:
	CMP	r6, #GRAPH_H/2
	MOVHI	r6, #GRAPH_H/2
	RSB	r8, r5, r5, lsl #0x04 @ Start at Fh (full opacity/brightness)
0:	SUBS	r9, r6, #0x0F         @ Handle small bars
	MULCC	r8, r5, r6
	BLS	2f
1:	LDR	r6, [r7]
	ORR	r6, r6, r8
	STR	r6, [r7], #-0x04
	SUBS	r9, r9, #0x01
	BNE	1b
2:	LDR	r9, [r7]
	ORR	r9, r9, r8
	STR	r9, [r7], #-0x04
	SUBS	r8, r8, r5
	BHI	2b
4:	BX	lr

.LRedraw_BufOfs: .word -BLOCK_SIZE*2
.LRedraw_PosMu:  .word 0
.LRedraw_LowPassL: .word 0
.LRedraw_LowPassR: .word 0
.LRedraw_PowerOldL: .word 0
.LRedraw_PowerOldR: .word 0

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
	.asciz "Some Artist"
.size SoundFile_Artist, .-SoundFile_Artist

SoundFile_Title:
	.asciz "Song Name"
.size SoundFile_Title, .-SoundFile_Title

/**************************************/
/* EOF                                */
/**************************************/
