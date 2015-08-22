# Makefile for building C stuff with GSL

C_SRCS= drawing.c aux.c noise.c init.c evolve.c main.c bc.c
C_OBJS= $(C_SRCS:.c=.o)
CFLAGS= -L/usr/local/lib/ -Wall -Wno-unused-variable -Wno-unused-result -std=gnu99 -Ofast -flto
ICFLAGS= -L/usr/local/include
LDFLAGS= -lgsl -lgslcblas -lmgl -lm
CC=gcc

drawing.o: drawing.c network.h
	$(CC) $(CFLAGS) -c drawing.c $(LDFLAGS)

aux.o: aux.c network.h
	$(CC) $(CFLAGS) -c aux.c $(LDFLAGS)

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
	$(CC) $(CFLAGS) -o simulate.x $(C_OBJS) $(LDFLAGS)
	cp simulate.x ../debug/

clean:
	rm -f *~ *.o simulate.x
	rm -f ../debug/simulate.x
