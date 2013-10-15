

INCLUDES= -I. -I$(RDIST)/include -I/usr/lib64/R/include -I/usr/share/R/include/ 
CC = gcc
LIBS   = -lm 

CFLAGS = -pthread -fPIC -Wall -g -O3   $(LIBS)  $(INCLUDES) 
#CFLAGS = -O3  $(LIBS)
RFLAGS = -pthread -fPIC -Wall -c -g -03 $(LIBS) -L/usr/X11/lib $(INCLUDES)

all: estlikC.o msR.so  msC.o stats_popC.o output_cmdC.o  output_cmdR.so  stats_popR.so  

estlikC.o:  estlikC.c streecC.c rand1.c msfun.c statsfun.c 
	CC ${CFLAGS} -o estlikC  estlikC.c  streecC.c rand1.c msfun.c statsfun.c 

msC.o:  msC.c streecC.c rand1.c msfun.c
	CC ${CFLAGS} -o msC msC.c streecC.c rand1.c msfun.c

stats_popC.o: stats_popC.c statsfun.c 
	CC ${CFLAGS} -o stats_popC stats_popC.c  statsfun.c 

output_cmdC.o: output_cmdC.c 
	CC ${CFLAGS} -o output_cmdC output_cmdC.c 

msR.so:  msR.c streecR.c rand1.c msfun.c 
	R CMD SHLIB ${RFLAGS} -o msR.so msR.c streecR.c rand1.c msfun.c 

stats_popR.so: stats_popR.c statsfun.c 
	R CMD SHLIB ${RFLAGS} -o stats_popR.so stats_popR.c  statsfun.c 

output_cmdR.so: output_cmdR.c 
	R CMD SHLIB ${RFLAGS} -o output_cmdR.so output_cmdR.c 

clean:
	rm -f estlikC  
	rm -f msC  msR.so
	rm -f  stats_popC  stats_popR.so 	
	rm -f output_cmdC  output_cmdR.so 
	rm -r *.dSYM 
	rm -r *~	

listing:	
	nm -g msR.so
	nm -g stats_popR.so
	nm -g output_cmdR.so


