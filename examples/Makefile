# platform
ifeq ($(OS),Windows_NT)
  RM       := del /Q
  EXE_EXT  := .exe
  OBJ_EXT  := .obj
  NULLDEV  := NUL
  SEP      := ;
else
  RM       := rm -f
  EXE_EXT  :=
  OBJ_EXT  := .o
  NULLDEV  := /dev/null
  SEP      := :
endif

CC       ?= gcc
CFLAGS   := -std=c11 -O2 -Wall -Wextra -Werror -pedantic -I../ -I/opt/homebrew/opt/openblas/include/
LDFLAG   := -L/opt/homebrew/opt/openblas/lib/ -lopenblas

PROG     := ex1

SRCS     := ex1.c
OBJS     := $(SRCS:.c=$(OBJ_EXT))
DEPS     := $(OBJS:.o=.d)
TARGET   := $(PROG)$(EXE_EXT)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAG)

# Rule to generate dependency files
%$(OBJ_EXT): %.c
	$(CC) -c $(CFLAGS) -MMD -MP -MF $*.d -o $@ $<

# Include dependency files
-include $(DEPS)

clean:
	-$(RM) $(OBJS) $(DEPS) $(TARGET) 2>$(NULLDEV)
