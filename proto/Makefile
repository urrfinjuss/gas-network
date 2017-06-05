#Makefile for building C stuff with GSL

C_SRCS= input.c output.c message.c init.c main.c
C_OBJS= $(C_SRCS:.c=.o)

OS := $(shell uname)
ifeq ($(OS), Linux)
CFLAGS= -I/home/sdyachen/usr/include -Wall -Wno-unused-variable -Wno-unused-result -std=gnu99 -Ofast -flto
LFLAGS= -L/home/sdyachen/usr/lib -Wall 
else
CFLAGS= -I"/home/Orange/usr/include" -Wall -Wno-unused-variable -Wno-unused-result -std=gnu99 -Ofast -flto
LFLAGS= -L"/home/Orange/usr/lib/" -Wall
endif

LDFLAGS= -lgsl -lgslcblas -lm -lfftw3
CC=gcc

init.o: init.c header.h
	$(CC) $(CFLAGS) -c init.c $(LDFLAGS)

message.o: message.c header.h
	$(CC) $(CFLAGS) -c message.c $(LDFLAGS)

output.o: output.c header.h
	$(CC) $(CFLAGS) -c output.c $(LDFLAGS)

input.o: input.c header.h
	$(CC) $(CFLAGS) -c input.c $(LDFLAGS)

main.o: main.c header.h
	$(CC) $(CFLAGS) -c main.c $(LDFLAGS)

proto: $(C_OBJS)
	echo "$(OS)"
	$(CC) -o proto.x $(C_OBJS) $(LFLAGS) $(LDFLAGS)
	mkdir -p ./out
	mkdir -p ./ic

clean:
	rm -f *~ *.o
