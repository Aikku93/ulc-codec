#----------------------------#
# Directories
#----------------------------#

OBJDIR := build
RELDIR := release

INCDIR := include
COMMON_SRCDIR := fourier libulc
ENCODETOOL_SRCDIR := encodetool
DECODETOOL_SRCDIR := decodetool

#----------------------------#
# Cross-compilation, compile flags
#----------------------------#

ARCHCROSS :=
ARCHFLAGS := -msse -msse2 -mavx -mavx2

CCFLAGS := -O2 -Wall -Wextra $(foreach dir, $(INCDIR), -I$(dir))
LDFLAGS := -static

#----------------------------#
# Tools
#----------------------------#

CC := $(ARCHCROSS)g++
LD := $(ARCHCROSS)g++

#----------------------------#
# Files
#----------------------------#

COMMON_SRC     := $(foreach dir, $(COMMON_SRCDIR), $(wildcard $(dir)/*.cpp))
ENCODETOOL_SRC := $(foreach dir, $(ENCODETOOL_SRCDIR), $(wildcard $(dir)/*.cpp))
DECODETOOL_SRC := $(foreach dir, $(DECODETOOL_SRCDIR), $(wildcard $(dir)/*.cpp))
COMMON_OBJ     := $(addprefix $(OBJDIR)/, $(notdir $(COMMON_SRC:.cpp=.o)))
ENCODETOOL_OBJ := $(addprefix $(OBJDIR)/, $(notdir $(ENCODETOOL_SRC:.cpp=.o)))
DECODETOOL_OBJ := $(addprefix $(OBJDIR)/, $(notdir $(DECODETOOL_SRC:.cpp=.o)))
ENCODETOOL_EXE := ulcencodetool.exe # Change this for other platforms
DECODETOOL_EXE := ulcdecodetool.exe # Change this for other platforms

VPATH := $(COMMON_SRCDIR) $(ENCODETOOL_SRCDIR) $(DECODETOOL_SRCDIR)

#----------------------------#
# General rules
#----------------------------#

$(OBJDIR)/%.o : %.cpp
	@echo $(notdir $<)
	@$(CC) $(CCFLAGS) -c -o $@ $<

#----------------------------#
# make common
#----------------------------#

.phony: common
common : $(COMMON_OBJ)

$(COMMON_OBJ) : $(COMMON_SRC) | $(OBJDIR)

#----------------------------#
# make encodetool
#----------------------------#

.phony: encodetool
encodetool : $(ENCODETOOL_EXE)

$(ENCODETOOL_OBJ) : $(ENCODETOOL_SRC) | $(OBJDIR)

$(ENCODETOOL_EXE) : $(COMMON_OBJ) $(ENCODETOOL_OBJ) | $(RELDIR)
	$(LD) -o $(RELDIR)/$@ $^ $(LDFLAGS)

#----------------------------#
# make decodetool
#----------------------------#

.phony: decodetool
decodetool : $(DECODETOOL_EXE)

$(DECODETOOL_OBJ) : $(DECODETOOL_SRC) | $(OBJDIR)

$(DECODETOOL_EXE) : $(COMMON_OBJ) $(DECODETOOL_OBJ) | $(RELDIR)
	$(LD) -o $(RELDIR)/$@ $^ $(LDFLAGS)

#----------------------------#
# make all
#----------------------------#

all : common encodetool decodetool

$(OBJDIR) $(RELDIR) :; mkdir -p $@

#----------------------------#
# make clean
#----------------------------#

.phony: clean
clean :; rm -rf $(OBJDIR) $(RELDIR)

#----------------------------#
