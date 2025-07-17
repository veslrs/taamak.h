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
CFLAGS   := -std=c11 -Wall -Wextra -Werror -Wno-unused-function -I../

PROG     := ex1

SRCS     := $(wildcard *.c)
OBJS     := $(SRCS:.c=$(OBJ_EXT))
TARGET   := $(PROG)$(EXE_EXT)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo Linking $@ ...
	$(CC) $(CFLAGS) -o $@ $^

%$(OBJ_EXT): %.c
	@echo Compiling $< ...
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@echo Cleaning up...
	-$(RM) $(OBJS) $(TARGET) 2>$(NULLDEV)
