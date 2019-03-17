CFLAGS=-c

INCLUDEPATH= -I.

LIBRARYPATH=

LIBS=-lm

SRC=sim_cuda.cpp simulation.cpp ldpc/ldpc.cpp ldpc/decoder_layered.cpp

OBJ=sim_cuda.o simulation.o ldpc.o decoder_layered.o

all: $(OBJ)
	g++ $(OBJ) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) -o sim_cuda

$(OBJ): $(SRC)
	g++ -std=c++11 $(CFLAGS) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) $(SRC)

clean:
	rm *.o
