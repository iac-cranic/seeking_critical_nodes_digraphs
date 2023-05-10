SDIR=src
ODIR=obj
IDIR=include
BDIR=bin

CC=gcc
CFLAGS=-pedantic -Wall -g -I$(IDIR) 
CFLAGS+=-Xpreprocessor -fopenmp -O3
LDFLAGS=-lm
MKDIR_P=mkdir -p
SH=/bin/bash
UNAME_S=$(shell uname -s)

#export OMP_NUM_THREADS=3
ifeq ($(UNAME_S),Darwin)
	LDFLAGS+=-lomp
	CFLAGS+=-DDARWIN
else
	LDFLAGS+=-lgomp -lbsd
endif


vpath %.c $(TDIR)
vpath %.c $(SDIR) 

.PHONY: all clean createDir 

TARGET= $(BDIR)/digraphCriticalNodesBF 


# Compile C files
$(ODIR)/%.o: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ -c $<


# Toolkit files
all: createDir $(TARGET) 

# Create Directories
createDir:
	${MKDIR_P} ${ODIR}
	${MKDIR_P} ${BDIR}

$(BDIR)/digraphCriticalNodesBF: $(ODIR)/digraphCriticalNodesBF.o $(ODIR)/graph.o $(ODIR)/utility.o  
	$(CC) $(LDFLAGS) -o $@ $^	

# Clean files
clean:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*
	
	

