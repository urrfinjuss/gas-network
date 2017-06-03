#Makefile for building C stuff with GSL


C_SRCS= adiabatic.c waux.c noise.c init.c evolve.c main.c bc.c
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

adiabatic.o: adiabatic.c network.h
	$(CC) $(CFLAGS) -c adiabatic.c $(LDFLAGS)

waux.o: waux.c network.h
	$(CC) $(CFLAGS) -c waux.c $(LDFLAGS)

noise.o: noise.c network.h
	$(CC) $(CFLAGS) -c noise.c $(LFLAGS) $(LDFLAGS)

init.o: init.c network.h
	$(CC) $(CFLAGS) -c init.c $(LDFLAGS)

evolve.o: evolve.c network.h
	$(CC) $(CFLAGS) -c evolve.c $(LDFLAGS)

main.o: main.c network.h
	$(CC) $(CFLAGS) -c main.c $(LDFLAGS)

bc.o:	bc.c network.h
	$(CC) $(CFLAGS) -c bc.c $(LDFLAGS)

simulate: $(C_OBJS)
	echo "$(OS)"
	$(CC) -o simulate.x $(C_OBJS) $(LFLAGS) $(LDFLAGS)
	cp simulate.x ./demo/simulation_1/
	mkdir -p ./demo/simulation_1/network
	mkdir -p ./demo/simulation_1/pipe_000
	mkdir -p ./demo/simulation_1/pipe_001
	mkdir -p ./demo/simulation_1/pipe_002
	mkdir -p ./demo/simulation_1/pipe_003
clean:
	rm -f *~ *.o simulate.x
