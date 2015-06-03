# Makefile for building C stuff with GSL

C_SRCS= aux.c gsl_noise.c init.c main.c
C_OBJS= $(C_SRCS:.c=.o)
CFLAGS= -L/usr/local/lib/ -std=gnu99 -Ofast -flto
LDFLAGS= -lgsl -lgslcblas -lm
CC=gcc

aux.o: aux.c network.h
	$(CC) -c aux.c

gsl_noise.o: gsl_noise.c network.h
	$(CC) $(CFLAGS) -c gsl_noise.c $(LDFLAGS)

init.o: init.c network.h
	$(CC) $(CFLAGS) -c init.c $(LDFLAGS)

main.o: main.c network.h
	$(CC) $(CFLAGS) -c main.c $(LDFLAGS)

gsl_noise: aux.o gsl_noise.o init.o main.o
	$(CC) $(CFLAGS) -o gsl_noise.x aux.o gsl_noise.o init.o main.o $(LDFLAGS)
	cp gsl_noise.x ../debug/

clean:
	rm -f *~ *.o gsl_noise.x
	rm -f ../debug/gsl_noise.x
