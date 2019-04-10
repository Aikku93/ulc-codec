/**************************************/
.section .iwram, "ax", %progbits
.balign 4
/**************************************/

.arm
_IRQProc:
	MOV	r0, #0x04000000
	LDR	r1, [r0, #0x0208]     @ IME -> r1
	LDR	r2, [r0, #0x0200]!    @ IE,IF -> r2, &IE -> r0
	ANDS	r1, r1, r1            @ IME?
	ANDNES	r1, r2, r2, lsr #0x10 @ IE&IF -> r1?
1:	RSBNE	r2, r1, #0x00       @ Mask to lowest bit
	ANDNE	r1, r1, r2
	STRNEH	r1, [r0, #0x02]     @ Set IF
	LDRNE	r2, [r0, #-0x0208]! @ BIOS.IRQFlags -> r2, &BIOS.IRQFlags -> r0
	ORRNE	r2, r2, r1          @ BIOS.IRQFlags |= Bit
	STRNEH	r2, [r0]
2:	LDRNE	r0, =_IRQTable @ Call handler
	LDRNE	r2, =0x077CB531
	MULNE	r3, r2, r1
	ADRNE	r1, .LIRQLog2
	LDRNEB	r1, [r1, r3, lsr #0x20-5]
	LDRNE	r0, [r0, r1, lsl #0x02]
	CMPNE	r0, #0x00
3:	MRSNE	ip, spsr
	STRNE	lr, [sp, #-0x0C]! @ Save lr_irq
	STMNEIB	sp, {ip,lr}^      @ Save spsr,lr_sys
	MSRNE	cpsr, #0x1F       @ SYS mode, free to interrupt
	ADRNE	lr, 0f
	BXNE	r0
	BX	lr
0:	MSR	cpsr, #0x92       @ IRQ mode, IRQ-block
	LDMIB	sp, {ip,lr}^      @ Restore spsr,lr_sys
	MSR	spsr, ip
	LDR	pc, [sp], #0x0C   @ Return to BIOS

.LIRQLog2:
	.byte  0, 1,28, 2,29,14,24,3
	.byte 30,22,20,15,25,17, 4,8
	.byte 31,27,13,23,21,19,16,7
	.byte 26,12,18, 6,11, 5,10,9

/**************************************/
.size   _IRQProc, .-_IRQProc
.global _IRQProc
/**************************************/
.section .bss
.balign 4
/**************************************/

_IRQTable:
	.space 16*4

/**************************************/
.size   _IRQTable, .-_IRQTable
.global _IRQTable
/**************************************/
/* EOF                                */
/**************************************/
