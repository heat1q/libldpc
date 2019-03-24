CFLAGS=-x cu -dc

INCLUDEPATH= -I.

LIBRARYPATH=

LIBS=-lm

SRC=sim_cuda.cu simulation.cu ldpc/ldpc.cu ldpc/decoder.cu

OBJ=sim_cuda.o simulation.o ldpc.o decoder.o

all: $(OBJ)
	nvcc $(OBJ) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) -o sim_cuda

$(OBJ): $(SRC)
	nvcc -std=c++11 $(CFLAGS) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) $(SRC)

clean:
	rm *.o
