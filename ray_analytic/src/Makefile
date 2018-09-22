CC=gcc
CFLAGS=-I.
LIBS=-lm -fopenmp
DEPS=init.h integral.h test_function.h
TFUNC=
OBJ=ray3d_analytic.o init.o integral.o ${TFUNC%.*}.o
ODIR=../binary

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

$(ODIR)/ray3d_grid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o
