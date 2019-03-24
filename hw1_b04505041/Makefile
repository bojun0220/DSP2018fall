.PHONY: all clean
CC=gcc
CFLAGS+=

all: train test
	$(CC) -o train train.cpp
	$(CC) -o test test.cpp 
clean:
	rm -f *.o train test   # type make clean to remove the compiled file
