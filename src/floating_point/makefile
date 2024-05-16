#******************************************************************************
#                        ETSI TS 103 634 V1.5.1                               *
#              Low Complexity Communication Codec Plus (LC3plus)              *
#                                                                             *
# Copyright licence is solely granted through ETSI Intellectual Property      *
# Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# estoppel or otherwise.                                                      *
#*****************************************************************************/
# GNU Makefile

# Options
AFL                      = 0
CLANG                    = 0
GCOV                     = 0
NO_POST_REL_CHANGES_TEST = 0
OPTIM                    = 0
SHORT_PLC_FADEOUT        = 0

# Paths
VPATH  = .
BUILD  = build
CC     = gcc
LINK   = $(CC)

# Binary Name
NAME_LC3   = LC3plus
# Shared Library Name
LIB_LC3    = libLC3plus.so

# Default tool settings
CC        = gcc
RM        = rm -f

# Preprocessor(-I/-D) / Compiler / Linker flags
CFLAGS   += -std=c99 -fPIC                                            \
             -pedantic -Wcast-qual -Wall -W -Wextra -Wno-long-long     \
             -Wpointer-arith -Wstrict-prototypes -Wmissing-prototypes  \
             -Werror-implicit-function-declaration

LDFLAGS += -lm -g
# Include dependency flags
DEPFLAGS = -MT $@ -MMD -MP -MF $(BUILD)/$*.Td

ifneq "$(DEBUG)" "0"
    CFLAGS   += -g3
    LDFLAGS  += -g3
endif

ifeq "$(GCOV)" "1"
    CFLAGS  += -fprofile-arcs -ftest-coverage
    LDFLAGS += -fprofile-arcs -ftest-coverage
endif

OPTIM    ?= 0
CFLAGS   += -O$(OPTIM)
CFLAGS   += $(foreach DIR,$(SRC_DIRS),-I$(DIR))

ifeq "$(NO_POST_REL_CHANGES_TEST)" "1"
CFLAGS   += -DNO_POST_REL_CHANGES
endif

ifeq "$(SHORT_PLC_FADEOUT)" "1"
CFLAGS   += -DPLC_TUNING_SHORT_FADEOUT
endif

# memory sanitizer, find use of uninitialized memory
ifeq "$(CLANG)" "1"
    CC        = clang
    CFLAGS   += -fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer
    LDFLAGS  += -fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer
    OPTIM     = 2
endif

# address sanitizer, find buffer overflows
ifeq "$(CLANG)" "2"
    CC        = clang
    CFLAGS   += -fsanitize=address -fno-omit-frame-pointer
    LDFLAGS  += -fsanitize=address -fno-omit-frame-pointer
    OPTIM     = 2
endif

# undefined behavior sanitizer, find bugs like integer overflows
ifeq "$(CLANG)" "3"
    CC       = clang
    CFLAGS  += -fsanitize=undefined
    LDFLAGS += -fsanitize=undefined
    OPTIM    = 2
endif

# for code coverate test
ifeq "$(GCOV)" "1"
    CFLAGS  += -fprofile-arcs -ftest-coverage
    LDFLAGS += -fprofile-arcs -ftest-coverage
endif

# verbose output
ifneq "$(VERBOSE)" "1"
    QUIET = @
endif

# dependency magic
CC_FLAGS    = '$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS)'
POSTCOMPILE = mv -f $(BUILD)/$*.Td $(BUILD)/$*.d && touch $@

######## Generate source / object lists ########

SRCS := $(notdir $(foreach DIR, $(VPATH), $(wildcard $(DIR)/*.c)))
SRCS_LIB := $(notdir $(foreach DIR, $(VPATH), $(wildcard $(DIR)/*.c)))
OBJS := $(addprefix $(BUILD)/, $(SRCS:.c=.o))

LIBSRCS := $(filter-out $(DIR)/ccConvert.c $(DIR)/codec_exe.c, $(SRCS))
LIBOBJS := $(addprefix $(BUILD)/, $(LIBSRCS:.c=.o))

.PHONY: all clean help force

.PRECIOUS: $(BUILD)/%.d

all: $(NAME_LC3)

help:
	@echo 'Targets:'
	@echo '    $(NAME_LC3) (default)'
	@echo '    $(LIB_LC3)'
	@echo 'Syntax: make [OPTION=VALUE ...]'
	@echo 'Build options:'
	@echo '    NO_POST_REL_CHANGES_TEST $(NO_POST_REL_CHANGES_TEST) [0,1]'
	@echo '    OPTIM                    $(OPTIM) [0-3]'
	@echo '    SHORT_PLC_FADEOUT        $(SHORT_PLC_FADEOUT) [0,1]'
	@echo 'Debug options:'
	@echo '    AFL                      $(AFL) [0,1]'
	@echo '    CLANG                    $(CLANG) [0-3]'
	@echo '    GCOV                     $(GCOV) [0,1]'

$(NAME_LC3): $(OBJS)
	@echo 'Linking' $@
	$(QUIET) $(LINK)  $(OBJS) -o $@ $(LDFLAGS)

$(LIB_LC3): $(LIBOBJS)
	@echo 'Linking' $@
	$(QUIET) $(LINK) --shared $(OBJS) -o $@ $(LDFLAGS)

clean:
	$(QUIET) rm -rf $(NAME_LC3) $(LIB_LC3) $(BUILD) ccConvert

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
