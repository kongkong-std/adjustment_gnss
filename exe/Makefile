CC = gcc
CFLAGS = -Wall -g -c
C_SRC_PATH = ../src
OBJECT := $(patsubst %.c,%.o,$(wildcard $(C_SRC_PATH)/*.c))

.PHONY: all clean
all: adjust_exe

adjust_exe: $(OBJECT)
	$(CC) $(OBJECT) -o adjust_exe

%.o: %.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm adjust_exe $(OBJECT)
