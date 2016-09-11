CC = gcc
FLAGS = -c # -g  -Wall 
LIBS   = -lm
LINKTO = -o mimar


all: mimar.o params.o rand1.o streec.o
	$(CC) $(LINKTO) -O3 mimar.o  params.o rand1.o streec.o $(LIBS)

mimar.o: mimar.c mimar.h
	$(CC) $(FLAGS) mimar.c

params.o: params.c mimar.h 
	$(CC) $(FLAGS) params.c

rand1.o: rand1.c 
	$(CC) $(FLAGS) rand1.c 

streec.o: streec.c 
	$(CC) $(FLAGS) streec.c

clean:
	rm -rf *.o *~
 
