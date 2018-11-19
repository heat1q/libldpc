CFLAGS=-c -pthread -Wall

INCLUDEPATH= -I.

LIBRARYPATH=

LIBS=-lm -fopenmp

SRC=ldpc_main.c functions.c ldpc_cycles.c ldpc_decoder.c scm_functions.c ldpc_stoppingsets.c

OBJ=ldpc_main.o functions.o ldpc_cycles.o ldpc_decoder.o scm_functions.o ldpc_stoppingsets.o

all: $(OBJ)
	gcc $(OBJ) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) -o st_ldpc

$(OBJ): $(SRC)
	gcc -std=c11 $(CFLAGS) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) $(SRC)

clean:
	rm -f $(OBJ)
