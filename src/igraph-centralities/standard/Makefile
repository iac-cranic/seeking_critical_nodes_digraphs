SDIR =.
IDIR=../include

GCC_NAME    =  gcc
GCC_FLAGS   =  -pedantic  -Wall  -std=c99  -I$(IDIR) -I$(IGRAPH_HEADERS) -ligraph
#GCC_FLAGS 	+= -DDEBUG 
GCC_FLAGS 	+= -g -O3 
TARGET = igraph_standard igraph_random 

# CHANGE THESE LINES TO USE YOUR FAVORITE COMPILER
CCC      = $(GCC_NAME)
FLAGS    = $(GCC_FLAGS)

vpath %.c $(SDIR)

.PHONY: all clean install distclean 

%.o: %.c
	$(CCC) $(FLAGS) -o $@ -c $< 


### MAKE ALL

all: $(TARGET) 




igraph_standard: igraph_standard.o ../common/utils.o ../common/wrapper_igraph.o 
	$(CCC) -o $@ $^ $(FLAGS) 

igraph_random: igraph_random.o ../common/utils.o ../common/wrapper_igraph.o 
	$(CCC) -o $@ $^ $(FLAGS) 

../common/utils.o:
	$(CCC) $(FLAGS)  -o ../common/utils.o -c ../common/utils.c
	
../common/wrapper_igraph.o:
	$(CCC) $(FLAGS)  -o ../common/wrapper_igraph.o -c ../common/wrapper_igraph.c

### CLEAN

clean:
	rm -v *.o ../common/*.o
	
distclean:
	rm -v *.o $(TARGET) ../common/*.o




