
CC = gcc
LIBS   = -lm 
CFLAGS =  -pthread -fPIC -Wall -g -O3  $(LIBS) 

all: stats_pop.o ms.o #EstL.o    output_cmd.o 

EstL.o:  EstL.c streec.c rand1.c msfun.c statsfun.c 
	CC ${CFLAGS} -o EstL  EstL.c  streec.c rand1.c msfun.c statsfun.c 

ms.o:  ms.c streec.c rand1.c msfun.c
	CC ${CFLAGS} -o ms ms.c streec.c rand1.c msfun.c

stats_pop.o: stats_pop.c statsfun.c 
	CC ${CFLAGS} -o stats_pop stats_pop.c  statsfun.c 

output_cmd.o: output_cmd.c 
	CC ${CFLAGS} -o output_cmd output_cmd.c 


clean:
	rm -f EstL  
	rm -f ms  stats_pop output_cmd 
	rm -r *~	
	rm -f *.dSYM

listing:	
	nm -g ms.o
	nm -g stats_pop.o
	nm -g output_cmd.o


