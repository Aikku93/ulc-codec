/**************************************/
.section .entry, "ax", %progbits
.balign 4
/**************************************/

.arm
_entry:
	B	.LEntry

.LHeader:
	.fill   156,1,0   @ [Nintendo logo]
	.fill	16,1,0    @ [Game title]
	.byte   32,32     @ [Maker code]
	.byte   0x96      @ Fixed value
	.byte   0x00      @ Unit code
	.byte   0x00      @ Device type
	.fill	7,1,0     @ [Unused]
	.byte	0x00      @ Software version
	.byte	0x11      @ Complement
	.byte	0x00,0x00 @ Checksum

.LEntry:
	MSR	cpsr, #0x12|0x80  @ Setup stacks [IRQ off]
	LDR	sp, =__sp_irq
	MSR	cpsr, #0x13|0x80
	LDR	sp, =__sp_svc
	MOV	r4, #0x04000000   @ &REG_BASE -> r4
	STR	r4, [r4, #0x200]  @ REG_IE = 0
	MOV	ip, #0x01
	STR	ip, [r4, #0x208]  @ IME = 1
	MSR	cpsr, #0x1F       @ [IRQ on]
	LDR	sp, =__sp_usr
0:	ADD	ip, pc, #0x01     @ Switch THUMB
	BX	ip

.thumb
.LCopyLMAVMA:
0:	LDR	r0, =__ewram_beg__     @ Copy EWRAM
	LDR	r1, =__ewram_lma_beg__
	LDR	r2, =__ewram_lma_end__
	BL	.LCopyLMAVMA_CopyLoop
1:	LDR	r1, =__ewram_end__     @ Clear SBSS
	BL	.LCopyLMAVMA_ClearLoop
0:	LDR	r0, =__iwram_beg__     @ Copy IWRAM
	LDR	r1, =__iwram_lma_beg__
	LDR	r2, =__iwram_lma_end__
	BL	.LCopyLMAVMA_CopyLoop
1:	LDR	r1, =__iwram_end__+4   @ Clear BSS (and BIOS.IRQFlag. Will return &BIOS.IRQProc in r0)
	BL	.LCopyLMAVMA_ClearLoop

.LConsoleInit:
1:	LDR	r3, =_IRQProc @ Set BIOS.IRQProc
	STR	r3, [r0]
2:	LDR	r0, =main
	BX	r0

/**************************************/

@ r0: &VMA
@ r1: &LMA_Beg
@ r2: &LMA_End
@ Returns end of copied VMA area

.LCopyLMAVMA_CopyLoop:
	CMP	r1, r2
	BCS	2f
1:	LDMIA	r1!, {r3}
	STMIA	r0!, {r3}
	CMP	r1, r2
	BCC	1b
2:	BX	lr

@ r0: &VMA
@ r1: &VMA_End

.LCopyLMAVMA_ClearLoop:
	CMP	r0, r1
	BCS	2f
0:	MOV	r3, #0x00
1:	STMIA	r0!, {r3}
	CMP	r0, r1
	BCC	1b
2:	BX	lr

/**************************************/
.size   _entry, .-_entry
.global _entry
/**************************************/
/* EOF                                */
/**************************************/
