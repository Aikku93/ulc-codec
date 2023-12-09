.phony: common
.phony: encodetool
.phony: decodetool
.phony: clean

#----------------------------#
# Directories
#----------------------------#

OBJDIR := build
RELDIR := release

INCDIR := include
COMMON_SRCDIR := fourier libulc
ENCODETOOL_SRCDIR := tools
DECODETOOL_SRCDIR := tools

#----------------------------#
# Cross-compilation, compile flags
#----------------------------#

ARCHCROSS :=
ARCHFLAGS :=

CCFLAGS := $(ARCHFLAGS) -fno-math-errno -ffast-math -O2 -Wall -Wextra $(foreach dir, $(INCDIR), -I$(dir))
LDFLAGS := -static -s

#----------------------------#
# Tools
#----------------------------#

CC := $(ARCHCROSS)gcc
LD := $(ARCHCROSS)gcc

#----------------------------#
# Files
#----------------------------#

COMMON_SRC     := $(foreach dir, $(COMMON_SRCDIR), $(wildcard $(dir)/*.c))
ENCODETOOL_SRC := $(filter-out $(ENCODETOOL_SRCDIR)/ulcDecodeTool.c, $(wildcard $(ENCODETOOL_SRCDIR)/*.c))
DECODETOOL_SRC := $(filter-out $(DECODETOOL_SRCDIR)/ulcEncodeTool.c, $(wildcard $(DECODETOOL_SRCDIR)/*.c))
COMMON_OBJ     := $(addprefix $(OBJDIR)/, $(addsuffix .o, $(COMMON_SRC)))
ENCODETOOL_OBJ := $(addprefix $(OBJDIR)/, $(addsuffix .o, $(ENCODETOOL_SRC)))
DECODETOOL_OBJ := $(addprefix $(OBJDIR)/, $(addsuffix .o, $(DECODETOOL_SRC)))
ENCODETOOL_EXE := $(RELDIR)/ulcencodetool
DECODETOOL_EXE := $(RELDIR)/ulcdecodetool

DFILES := $(COMMON_OBJ:.o=.d) $(ENCODETOOL_OBJ:.o=.d) $(DECODETOOL_OBJ:.o=.d)

#----------------------------#
# make all
#----------------------------#

all : encodetool decodetool common

$(OBJDIR) $(RELDIR) :; mkdir -p $@

#----------------------------#
# make common
#----------------------------#

common : $(COMMON_OBJ)

$(COMMON_OBJ) : $(COMMON_SRC) | $(OBJDIR)

#----------------------------#
# make encodetool
#----------------------------#

encodetool : $(ENCODETOOL_EXE)

$(ENCODETOOL_EXE) : $(COMMON_OBJ) $(ENCODETOOL_OBJ) | $(RELDIR)
	@echo Building encode tool $@...
	@$(LD) -o $@ $^ $(LDFLAGS)

$(ENCODETOOL_OBJ) : $(ENCODETOOL_SRC) | $(OBJDIR)

#----------------------------#
# make decodetool
#----------------------------#

decodetool : $(DECODETOOL_EXE)

$(DECODETOOL_EXE) : $(COMMON_OBJ) $(DECODETOOL_OBJ) | $(RELDIR)
	@echo Building decode tool $@...
	@$(LD) -o $@ $^ $(LDFLAGS)

$(DECODETOOL_OBJ) : $(DECODETOOL_SRC) | $(OBJDIR)

#----------------------------#
# Rules
#----------------------------#

$(OBJDIR)/%_AVX_FMA.c.o : %_AVX_FMA.c | $(OBJDIR)
	@echo $(notdir $<)
	@mkdir -p $(dir $@)
	@$(CC) $(CCFLAGS) -c -MD -MP -MF $(OBJDIR)/$<.d -o $@ $< -mavx -mfma

$(OBJDIR)/%_AVX.c.o : %_AVX.c | $(OBJDIR)
	@echo $(notdir $<)
	@mkdir -p $(dir $@)
	@$(CC) $(CCFLAGS) -c -MD -MP -MF $(OBJDIR)/$<.d -o $@ $< -mavx

$(OBJDIR)/%_SSE_FMA.c.o : %_SSE_FMA.c | $(OBJDIR)
	@echo $(notdir $<)
	@mkdir -p $(dir $@)
	@$(CC) $(CCFLAGS) -c -MD -MP -MF $(OBJDIR)/$<.d -o $@ $< -msse -mfma

$(OBJDIR)/%_SSE.c.o : %_SSE.c | $(OBJDIR)
	@echo $(notdir $<)
	@mkdir -p $(dir $@)
	@$(CC) $(CCFLAGS) -c -MD -MP -MF $(OBJDIR)/$<.d -o $@ $< -msse

$(OBJDIR)/%.c.o : %.c | $(OBJDIR)
	@echo $(notdir $<)
	@mkdir -p $(dir $@)
	@$(CC) $(CCFLAGS) -c -MD -MP -MF $(OBJDIR)/$<.d -o $@ $<

-include $(DFILES)

#----------------------------#
# make clean
#----------------------------#

clean :; rm -rf $(OBJDIR) $(RELDIR)

#----------------------------#
