CFLAGS=-c -pthread -Wall

INCLUDEPATH= -I.

LIBRARYPATH=

LIBS=-lrt -lm -fopenmp

SRC=main.c ../function/ldpc_functions.c ../LDPC_SSET/ldpc_cycles.c ldpc_decoder.c ldpc_decoder_layered.c ../function/scm_functions.c ../LDPC_SSET/ldpc_stoppingsets.c

OBJ=main.o ldpc_functions.o ldpc_cycles.o ldpc_decoder.o ldpc_decoder_layered.o scm_functions.o ldpc_stoppingsets.o

all: $(OBJ)
	gcc $(OBJ) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) -o decode

$(OBJ): $(SRC)
	gcc -std=c11 $(CFLAGS) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) $(SRC)

clean:
	rm -f $(OBJ)
