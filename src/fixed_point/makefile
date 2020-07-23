# Options
AFL         = 0
CLANG       = 0
GCOV        = 0
NO_POST_REL = 0
OPTIM       = 0
SUBSET      = ALL
WMOPS       = 1

# Paths
VPATH  = . basic_op
BUILD  = build
CC     = gcc
LINK   = $(CC)

# Binary Name
NAME_LC3   = LC3plus

# Default tool settings
RM        = rm -f

ifndef VERBOSE
QUIET_CC  = @echo '   ' Compiling $<;
QUIET_LINK= @echo '   ' Linking $@;
QUIET     = @
endif

CFLAGS += -std=c99

# C compiler flags
# Preprocessor(-I/-D) / Compiler / Linker flags
CPPFLAGS += -Ibasic_op -DSUBSET_$(SUBSET)
CFLAGS   += -pedantic -Wcast-qual -Wall -W -Wextra -Wno-long-long     \
            -Wpointer-arith -Wstrict-prototypes -Wmissing-prototypes  \
            -Werror-implicit-function-declaration

ifneq "$(DEBUG)" "0"
CFLAGS   += -g3
LDFLAGS  += -g3
endif

LDFLAGS += -lm

DEPFLAGS = -MT $@ -MMD -MP -MF $(BUILD)/$*.Td

ifeq "$(GCOV)" "1"
CFLAGS  += -fprofile-arcs -ftest-coverage
LDFLAGS += -fprofile-arcs -ftest-coverage
endif

OPTIM    ?= 0
CFLAGS   += -O$(OPTIM)

CFLAGS   += $(foreach DIR,$(SRC_DIRS),-I$(DIR))

ifeq "$(NO_POST_REL)" "1"
CFLAGS   += -DNO_POST_REL_CHANGES
endif

# disable wmops instrumentation
ifeq "$(WMOPS)" "0"
    CPPFLAGS += -DWMOPS=0 -DDONT_COUNT_MEM
endif

# dependency magic
CC_FLAGS    = '$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS)'
POSTCOMPILE = mv -f $(BUILD)/$*.Td $(BUILD)/$*.d && touch $@

###############################################################################

SRCS := $(notdir $(foreach DIR, $(VPATH), $(wildcard $(DIR)/*.c)))

EXCL := ccConvert.c
SRCS := $(notdir $(foreach DIR, $(VPATH), $(wildcard $(DIR)/*.c)))
SRCS := $(filter-out $(EXCL), $(SRCS))

EXCL_CCC := $(BUILD)/codec_exe.o
OBJS_CCC_UNF := $(addprefix $(BUILD)/, $(SRCS:.c=.o))
OBJS_CCC := $(filter-out $(EXCL_CCC), $(OBJS_CCC_UNF))

OBJS := $(addprefix $(BUILD)/, $(SRCS:.c=.o))

###############################################################################

.PHONY: all clean help force

all: $(NAME_LC3)

help:
	@echo 'Targets:'
	@echo '    $(NAME_LC3) (default)'
	@echo '    ccConvert'
	@echo 'Syntax: make [OPTION=VALUE ...]'
	@echo 'Build options:'
	@echo '    NO_POST_REL $(NO_POST_REL) [0,1]'
	@echo '    OPTIM       $(OPTIM) [0-3]'
	@echo '    SUBSET      $(SUBSET) [NB,SQ,HQ,SWB,FB,ALL]'
	@echo '    WMOPS       $(WMOPS) [0,1]'
	@echo 'Debug options:'
	@echo '    AFL         $(AFL) [0,1]'
	@echo '    CLANG       $(CLANG) [0-3]'
	@echo '    GCOV        $(GCOV) [0,1]'

$(NAME_LC3): $(OBJS)
	@echo 'Linking' $@
	$(QUIET) $(LINK)  $(OBJS) -o $@ $(LDFLAGS)
    
ccConvert: $(BUILD)/ccConvert.o $(OBJS_CCC)
	@echo 'Linking' $@
	$(QUIET) $(LINK) $? -o $@ $(LDFLAGS)

clean:
	$(QUIET) rm -rf $(NAME_LC3) $(BUILD) ccConvert

$(BUILD)/%.o : %.c $(BUILD)/cc_flags
	@echo 'Compiling' $<
	$(QUIET) $(CC) $(DEPFLAGS) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<
	$(QUIET) $(POSTCOMPILE)

# force rebuild if compilation flags changed
$(BUILD)/cc_flags: force
	$(QUIET) mkdir -p $(BUILD)
	$(QUIET) echo $(CC_FLAGS) | cmp -s - $@ || echo $(CC_FLAGS) > $@

# force rebuild if include dependency changed
$(BUILD)/%.d: ;
include $(wildcard $(patsubst %, $(BUILD)/%.d, $(basename $(SRCS))))
