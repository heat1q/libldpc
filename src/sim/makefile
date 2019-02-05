CFLAGS=-c -pthread -Wall

INCLUDEPATH= -I.

LIBRARYPATH=

LIBS=-lrt -lm -fopenmp

SRC=sim_ldpc.c ldpc_functions.c ldpc_decoder.c scm_functions.c

OBJ=sim_ldpc.o ldpc_functions.o ldpc_decoder.o scm_functions.o

all: $(OBJ)
	gcc $(OBJ) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) -o sim_ldpc

$(OBJ): $(SRC)
	gcc -std=c11 $(CFLAGS) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) $(SRC)

clean:
	rm -f $(OBJ)
