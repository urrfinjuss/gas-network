#Makefile for building C stuff with GSL


C_SRCS= adiabatic.c waux.c noise.c init.c evolve.c main.c bc.c
C_OBJS= $(C_SRCS:.c=.o)

OS := $(shell uname)
ifeq ($(OS), Linux)
CFLAGS= -L/usr/local/lib/ -Wall -Wno-unused-variable -Wno-unused-result -std=gnu99 -Ofast -flto
ICFLAGS= -L/usr/local/include
else
CFLAGS=  -L"C:/home/gas-tools/libraries/fftw3/" -I"C:/home/gas-tools/libraries/gsl-2.3/build" -Wall -Wno-unused-variable -Wno-unused-result -std=gnu99 -Ofast -flto
ICFLAGS= -Wall
endif

LDFLAGS= -lgsl -lgslcblas -lm -lfftw3
CC=gcc

adiabatic.o: adiabatic.c network.h
	$(CC) $(CFLAGS) -c adiabatic.c $(LDFLAGS)

waux.o: waux.c network.h
	$(CC) $(CFLAGS) -c waux.c $(LDFLAGS)

noise.o: noise.c network.h
	$(CC) $(CFLAGS) -c noise.c $(LDFLAGS)

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
	$(CC) $(CFLAGS) -o simulate.x $(C_OBJS) $(LDFLAGS)
	cp simulate.x ../debug/

clean:
	rm -f *~ *.o simulate.x
	rm -f ../debug/simulate.x
