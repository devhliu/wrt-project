CC=gcc
CFLAGS=-I.
LIBS=-lm -fopenmp -lgsl -lgslcblas
DEPS=init.h interpolation.h slicing.h
OBJ=reduction.o init.o interpolation.o slicing.o
ONAME=radon_reduction

%.o: %.c
	$(CC) -c -o $@ 

$(ONAME): $(OBJ)
	$(CC) -o $@ $(LIBS)

.PHONY: clean

clean:
rm -f *.o
