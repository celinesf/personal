### Makefile for mimarsim ###

CC     = gcc
FLAGS  = -c #-g -Wall
LIBS   = -lm
LINKTO = -o mimarsim

# Comment celine 10/09/2007 - Uncomment for Celine's use. #
#all: mimarsim.o params.o make_gametes.o rand1.o streec.o tajd.o
#	$(CC) $(LINKTO) -O3 mimarsim.o params.o make_gametes.o rand1.o streec.o tajd.o $(LIBS)
##

# Comment for Celine's use
all: mimarsim.o params.o make_gametes.o rand1.o streec.o
	$(CC) $(LINKTO) -O3 mimarsim.o params.o make_gametes.o rand1.o streec.o $(LIBS)
##

mimarsim.o: mimarsim.c mimarsim.h
	$(CC) $(FLAGS) mimarsim.c

params.o: params.c mimarsim.h
	$(CC) $(FLAGS) params.c

make_gametes.o: make_gametes.c
	$(CC) $(FLAGS) make_gametes.c

rand1.o: rand1.c
	$(CC) $(FLAGS) rand1.c

streec.o: streec.c
	$(CC) $(FLAGS) streec.c

# Comment celine 10/09/2007 - Uncomment for Celine's use. #
#tajd.o: tajd.c	
#	$(CC) $(FLAGS) tajd.c
##

clean:
	rm -rf *.o *~
 
