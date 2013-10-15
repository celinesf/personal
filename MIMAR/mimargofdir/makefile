### Makefile for mimargof ###

CC     = gcc
FLAGS  = -c #-g -Wall
LIBS   = -lm
LINKTO = -o mimargof

all: mimargof.o params.o make_gametes.o rand1.o streec.o tajd.o
	$(CC) $(LINKTO) mimargof.o params.o make_gametes.o rand1.o streec.o tajd.o $(LIBS)

mimargof.o: mimargof.c mimargof.h
	$(CC) $(FLAGS) mimargof.c

params.o: params.c mimargof.h
	$(CC) $(FLAGS) params.c

make_gametes.o: make_gametes.c
	$(CC) $(FLAGS) make_gametes.c

rand1.o: rand1.c
	$(CC) $(FLAGS) rand1.c

streec.o: streec.c
	$(CC) $(FLAGS) streec.c

tajd.o: tajd.c
	$(CC) $(FLAGS) tajd.c

clean:
	rm -rf *.o *~
