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

# Alternatively, try "-march=native" for ARCHFLAGS
ARCHCROSS :=
ARCHFLAGS := -msse -msse2 -mavx -mavx2 -mfma

CCFLAGS := $(ARCHFLAGS) -fno-math-errno -O2 -Wall -Wextra $(foreach dir, $(INCDIR), -I$(dir))
LDFLAGS := -static

#----------------------------#
# Tools
#----------------------------#

CC := $(ARCHCROSS)gcc
LD := $(ARCHCROSS)gcc

#----------------------------#
# Files
#----------------------------#

COMMON_SRC     := $(foreach dir, $(COMMON_SRCDIR), $(wildcard $(dir)/*.c))
ENCODETOOL_SRC := $(ENCODETOOL_SRCDIR)/ulcEncodeTool.c
DECODETOOL_SRC := $(ENCODETOOL_SRCDIR)/ulcDecodeTool.c
COMMON_OBJ     := $(addprefix $(OBJDIR)/, $(notdir $(COMMON_SRC:.c=.o)))
ENCODETOOL_OBJ := $(addprefix $(OBJDIR)/, $(notdir $(ENCODETOOL_SRC:.c=.o)))
DECODETOOL_OBJ := $(addprefix $(OBJDIR)/, $(notdir $(DECODETOOL_SRC:.c=.o)))
ENCODETOOL_EXE := ulcencodetool.exe # Change this for other platforms
DECODETOOL_EXE := ulcdecodetool.exe # Change this for other platforms

DFILES := $(wildcard $(OBJDIR)/*.d)

VPATH := $(COMMON_SRCDIR) $(ENCODETOOL_SRCDIR) $(DECODETOOL_SRCDIR)

#----------------------------#
# General rules
#----------------------------#

$(OBJDIR)/%.o : %.c
	@echo $(notdir $<)
	@$(CC) $(CCFLAGS) -c -o $@ $< -MMD -MP -MF $(OBJDIR)/$*.d

#----------------------------#
# make all
#----------------------------#

all : common encodetool decodetool

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

$(ENCODETOOL_OBJ) : $(ENCODETOOL_SRC) | $(OBJDIR)

$(ENCODETOOL_EXE) : $(COMMON_OBJ) $(ENCODETOOL_OBJ) | $(RELDIR)
	$(LD) -o $(RELDIR)/$@ $^ $(LDFLAGS)

#----------------------------#
# make decodetool
#----------------------------#

decodetool : $(DECODETOOL_EXE)

$(DECODETOOL_OBJ) : $(DECODETOOL_SRC) | $(OBJDIR)

$(DECODETOOL_EXE) : $(COMMON_OBJ) $(DECODETOOL_OBJ) | $(RELDIR)
	$(LD) -o $(RELDIR)/$@ $^ $(LDFLAGS)

#----------------------------#
# make clean
#----------------------------#

clean :; rm -rf $(OBJDIR) $(RELDIR)

#----------------------------#
# Dependencies
#----------------------------#

include $(wildcard $(DFILES))

#----------------------------#
