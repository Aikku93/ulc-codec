#----------------------------#
PATH := /x/Tools/GCCARM/bin:$(PATH)
#----------------------------#

SOURCES  := source source/res source/sys source/ulc
INCLUDES := include
BUILD    := build
TARGET   := ulcplayer

#----------------------------#

ARCHCROSS := arm-none-eabi-
ARCHFLAGS := -mcpu=arm7tdmi -mtune=arm7tdmi -mthumb-interwork -mlong-calls

INCLUDEFLAGS := $(foreach dir, $(INCLUDES), -I$(dir))

CCFLAGS := $(ARCHFLAGS) $(INCLUDEFLAGS) -O2 -Wall -Wextra
ASFLAGS := $(ARCHFLAGS) $(INCLUDEFLAGS) -x assembler-with-cpp
LDFLAGS := $(ARCHFLAGS)

CFILES := $(foreach dir, $(SOURCES), $(notdir $(wildcard $(dir)/*.c)))
SFILES := $(foreach dir, $(SOURCES), $(notdir $(wildcard $(dir)/*.s)))
OFILES := $(addprefix $(BUILD)/, $(CFILES:.c=.c.o) $(SFILES:.s=.s.o))

VPATH := $(SOURCES)

#----------------------------#

CC := $(ARCHCROSS)gcc
AS := $(ARCHCROSS)gcc
LD := $(ARCHCROSS)gcc
NM := $(ARCHCROSS)nm
OBJCOPY := $(ARCHCROSS)objcopy

#----------------------------#

all : $(TARGET).gba

$(TARGET).gba : $(BUILD)/$(TARGET).elf | $(BUILD)
	$(NM) -S -n $< > MemoryMap.txt
	$(OBJCOPY) -O binary $< $@

$(BUILD)/$(TARGET).elf : $(OFILES) | $(BUILD)
	$(LD) -nostdlib -T agb.ld -o $@ $^

$(BUILD)/%.c.o : %.c | $(BUILD)
	$(CC) $(CCFLAGS) -c -o $@ $<

$(BUILD)/%.s.o : %.s | $(BUILD)
	$(AS) $(ASFLAGS) -c -o $@ $<

$(BUILD):; mkdir -p $@

#----------------------------#
