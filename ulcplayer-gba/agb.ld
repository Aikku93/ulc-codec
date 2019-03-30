/****************************/
OUTPUT_FORMAT("elf32-littlearm", "elf32-bigarm", "elf32-littlearm")
OUTPUT_ARCH(arm)
ENTRY(_entry)
/****************************/

/* Memory Layout */
MEMORY {
	ewram : ORIGIN = 0x02000000, LENGTH = 256K
	iwram : ORIGIN = 0x03000000, LENGTH =  32K - 8 /* Exclude BIOS.IRQFlag,BIOS.IRQProc */
	rom   : ORIGIN = 0x08000000, LENGTH = 128M
}
__ewram_beg__ = ORIGIN(ewram);
__ewram_end__ = ORIGIN(ewram) + LENGTH(ewram);
__iwram_beg__ = ORIGIN(iwram);
__iwram_end__ = ORIGIN(iwram) + LENGTH(iwram);
__rom_beg__   = ORIGIN(rom);
__rom_end__   = ORIGIN(rom) + LENGTH(rom);

/* Stack Areas */
__sp_base = __iwram_end__;
__sp_irq_size = 0x00C0;
__sp_svc_size = 0x0040;
__sp_usr_size = 0x0300;
__sp_irq = __sp_base;
__sp_svc = __sp_base - __sp_irq_size;
__sp_usr = __sp_base - __sp_irq_size - __sp_svc_size;

/****************************/

SECTIONS {
	/* ROM area */
	.text : {
		KEEP (*(.entry))
		*(.text .text.*)
		*(.stub)
		*(.glue_7 .glue_7t)
	} >rom
	.rodata : {
		*(.rodata .rodata.*)
		*(.roda)
	} >rom
	__lma_area_beg__ = . ;

	/* EWRAM area */
	. = __ewram_beg__;
	__ewram_lma_beg__ = __lma_area_beg__ ;
		.ewram : AT(__ewram_lma_beg__) {
			*(.data .data.*)
			*(.ewram .ewram.*)
		} >ewram
	__ewram_lma_end__ = __ewram_lma_beg__ + SIZEOF(.ewram);
		.sbss (NOLOAD) : {
			*(.sbss .sbss.*)
		} >ewram
	__ewram_heap_beg__ = . ;
	__ewram_heap_end__ = __ewram_end__;

	/* IWRAM area */
	. = __iwram_beg__;
	__iwram_lma_beg__ = __ewram_lma_end__;
		.iwram : AT(__iwram_lma_beg__) {
			*(.iwram .iwram.*)
		} >iwram
	__iwram_lma_end__ = __iwram_lma_beg__ + SIZEOF(.iwram);
		.bss (NOLOAD) : {
			*(.bss .bss.*)
		} >iwram
	__iwram_heap_beg__ = . ;
	__iwram_heap_end__ = __sp_usr - __sp_usr_size;
}

/****************************/
