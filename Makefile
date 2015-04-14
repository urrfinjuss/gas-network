# Makefile for building C stuff with GSL

CFLAGS= -L/usr/local/lib/
LDFLAGS= -lgsl -lgslcblas -lm
CC=gcc

%: %.c
	$(CC) $(CFLAGS) $<  $(LDFLAGS) -o $@

# eg. do "make gsl_test" to make gsl_test from gsl_test.c
# then run with "gsl_test 10"

clean:
	rm -f *~ *.o core a.out
