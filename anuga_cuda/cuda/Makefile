# makefile for for COMP8320 Lab 07, 2011
# Peter Strazdins, RSCS ANU, Sep 2011

CC=nvcc  -arch=compute_20 -code=sm_20
CCFLAGS=-Xptxas -v  --compiler-options -fno-strict-aliasing -O3
LDFLAGS=-lcuda -lcudart -lm

OBJS=gravity.o

EXE=main

all: $(EXE)

$(EXE): $(EXE).o $(OBJS)
	$(CC) $(CCFLAGS) -o $(EXE) $(EXE).o $(OBJS) $(LDFLAGS)

%.o: %.cu
	$(CC) $(CCFLAGS) -c $*.cu

clean:
	rm -f $(OBJS) $(EXE) $(EXE).o *~
