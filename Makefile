# Makefile for building C stuff with GSL

C_SRCS= aux.c gsl_noise.c
C_OBJS= $(C_SRCS:.c=.o)
CFLAGS= -L/usr/local/lib/ -Ofast -flto
LDFLAGS= -lgsl -lgslcblas -lm
CC=gcc

aux.o: aux.c network.h
	$(CC) -c aux.c

gsl_noise.o: gsl_noise.c network.h
	$(CC) $(CFLAGS) -c gsl_noise.c $(LDFLAGS)

gsl_noise: aux.o gsl_noise.o
	$(CC) $(CFLAGS) -o gsl_noise.x aux.o gsl_noise.o $(LDFLAGS)
	cp gsl_noise.x ../debug/

clean:
	rm -f *~ *.o gsl_noise.x
	rm -f ../debug/gsl_noise.x
