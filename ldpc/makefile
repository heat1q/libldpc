CFLAGS=-c -pthread -Wall

INCLUDEPATH= -I.

LIBRARYPATH=

LIBS=

SRC=ldpc_main.c functions.c ldpc_cycles.c

OBJ=ldpc_main.o functions.o ldpc_cycles.o

all: $(OBJ)
	gcc $(OBJ) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) -o run

$(OBJ): $(SRC)
	gcc -std=c11 $(CFLAGS) $(INCLUDEPATH) $(LIBRARYPATH) $(LIBS) $(SRC)

clean:
	rm -f $(OBJ)
