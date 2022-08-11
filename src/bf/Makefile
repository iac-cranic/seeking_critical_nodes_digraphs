SDIR=src
ODIR=obj
IDIR=include
BDIR=bin

CC=gcc
CFLAGS=-pedantic -Wall -g -I$(IDIR) 
CFLAGS+=-Xpreprocessor -fopenmp -O3
#CFLAGS+=-Xpreprocessor -fopenmp
LDFLAGS=-lm
MKDIR_P=mkdir -p
SH=/bin/bash
UNAME_S=$(shell uname -s)

#export OMP_NUM_THREADS=3
ifeq ($(UNAME_S),Darwin)
	LDFLAGS+=-lomp
else
	LDFLAGS+=-lgomp -lbsd
endif


vpath %.c $(TDIR)
vpath %.c $(SDIR) 

.PHONY: all clean  

TARGET= $(BDIR)/digraphCriticalNodesBF 


# Compile C files
$(ODIR)/%.o: %.c
	${MKDIR_P} ${ODIR}
	$(CC) $(CFLAGS) -o $@ -c $<


# Toolkit files
all: $(TARGET) 
	
$(BDIR)/digraphCriticalNodesBF: $(ODIR)/digraphCriticalNodesBF.o $(ODIR)/graph.o $(ODIR)/utility.o $(IDIR)/graph.h $(IDIR)/utility.h 
	${MKDIR_P} ${BDIR}
	$(CC) $(LDFLAGS) -o $@ $^	

# Clean files
clean:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*
	
	

