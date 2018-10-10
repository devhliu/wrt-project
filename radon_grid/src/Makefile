CC=gcc
CFLAGS=-I. -fopenmp
LIBS=-lm -lgsl -lgslcblas
DEPS=init.h integral.h interp.h
OBJ=radon3d_grid.o init.o integral.o interp.o
ODIR=../binary

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

$(ODIR)/radon3d_grid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
