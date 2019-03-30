/**************************************/
.text
.balign 2
/**************************************/

.thumb
.thumb_func
main:
	LDR	r0, =SoundFile
	BL	ulc_Init
1:	BL	ulc_BlockProcess
	MOV	r0, #0x00 @ IntrWait (return if already set)
	MVN	r1, r0    @ Any interrupt
	SWI	0x04
	B	1b

/**************************************/
.size   main, .-main
.global main
/**************************************/
.section .rodata
.balign 4
/**************************************/

SoundFile: .incbin "source/SoundData.ulc"

/**************************************/
.size SoundFile, .-SoundFile
/**************************************/
/* EOF                                */
/**************************************/
