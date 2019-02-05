CFLAGS=-c -pthread -Wall

INCLUDEPATH= -I.

LIBRARYPATH=

LIBS=-lrt -lm -fopenmp

SRC=analysis_ldpc.c function/ldpc_functions.c LDPC_SSET/ldpc_cycles.c decoder/ldpc_decoder.c function/scm_functions.c LDPC_SSET/ldpc_stoppingsets.c

OBJ=analysis_ldpc.o ldpc_functions.o ldpc_cycles.o ldpc_decoder.o scm_functions.o ldpc_stoppingsets.o

all: $(OBJ)
	gcc $(OBJ) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) -o analysis_ldpc

$(OBJ): $(SRC)
	gcc -std=c11 $(CFLAGS) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) $(SRC)

clean:
	rm -f $(OBJ)
