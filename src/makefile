all: main

Compiler	= gcc
FLAGS		= -g -Wall -O0
LIB		= -lz -lm -lpthread
HEADER		= main.h
SOURCE		= main.c index.c bwa_index.c seed_align.c seed.c generate_cigar.c DP.c
PROG		= DEEPL
OBJECT		= $(SOURCE:%.c=%.o) 

all:		$(PROG)

main:		$(OBJECT)
			$(Compiler) $(FLAGS) $(OBJECT) -o $(PROG) $(LIB)
                        
%.o:		%.c $(HEADER)
			$(Compiler) $(FLAGS) -c $<

clean:
		rm -f *.o
