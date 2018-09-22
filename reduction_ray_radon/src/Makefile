CC=gcc
CFLAGS=-I. 
LIBS=-lm -fopenmp -lgsl -lgslcblas
DEPS=init.h interpolation.h slicing.h
OBJ=reduction.o init.o interpolation.o slicing.o
ODIR=../binary

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

$(ODIR)/reduction_rayradon3d: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
