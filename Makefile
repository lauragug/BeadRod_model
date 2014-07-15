PROG = beadRod
CC     = mpicc
#CFLAGS = -fopenmp -O2 -W -Wall 
#CFLAGS = -g3 -W -Wall 
CFLAGS = -O3 -W -Wall -I ~/Libraries/GSL/include -L ~/Libraries/lib/  
LIBS   = -lgsl -lgslcblas -llapack -lblas -lm -lgfortran

CFILES = $(wildcard *.c) 
OBJS    = $(CFILES:.c=.o)


all:   $(PROG)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

%.o:%.c 
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(OBJS) $(PROG) 1

