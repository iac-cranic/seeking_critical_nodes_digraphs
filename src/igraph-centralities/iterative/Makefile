SDIR =.
IDIR=../include

GCC_NAME    =  gcc
GCC_FLAGS   =  -pedantic  -Wall -g -std=c99 -I$(IDIR) -ligraph
GCC_FLAGS 	+= -O3 
TARGET = igraph_iterative

CCC      = $(GCC_NAME)
FLAGS    = $(GCC_FLAGS)

vpath %.c $(SDIR)

.PHONY: all clean install distclean 

%.o: %.c
	$(CCC) -I${IDIR} -g -o $@ -c $< 

all: $(TARGET) 


igraph_iterative:igraph_iterative.o ../common/utils.o ../common/wrapper_igraph.o 
	$(CCC) -o $@ $^ $(FLAGS) 

../common/utils.o:
	$(CCC) -g -I${IDIR} -o ../common/utils.o -c ../common/utils.c
	
../common/wrapper_igraph.o:
	$(CCC) -g -I${IDIR} -o ../common/wrapper_igraph.o -c ../common/wrapper_igraph.c

### CLEAN

clean:
	rm -v *.o ../common/*.o
	
distclean:
	rm -v *.o $(TARGET) ../common/*.o  




