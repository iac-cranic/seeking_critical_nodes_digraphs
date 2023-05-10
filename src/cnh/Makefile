SDIR = src
ODIR = obj
TDIR = test
IDIR = include
BDIR = bin

GCC_NAME    =  gcc
GCC_FLAGS   =  -pedantic  -Wall  -std=c99  -I$(IDIR)
GCC_FLAGS   += -g
GCC_FLAGS   += -O3

TARGET = $(BDIR)/cnh
TARGET_LARGE = $(TARGET:=_large)
# CHANGE THESE LINES TO USE YOUR FAVORITE COMPILER
CCC      = $(GCC_NAME)
FLAGS    = $(GCC_FLAGS)
MKDIR_P  = mkdir -p

vpath %.c $(TDIR) 
vpath %.c $(SDIR)

.PHONY: clean install createDir distclean 

#.SUFFIXES: .c

$(ODIR)/%.o: %.c
	$(CCC) $(FLAGS) -o $@ -c $<

$(ODIR)/%.large.o: %.c 
	$(CCC) $(FLAGS) -DLARGE -c $^ -o $@

### MAKE ALL
all: createDir $(TARGET) $(TARGET_LARGE)  

# Create Directories
createDir:
	${MKDIR_P} ${ODIR}
	${MKDIR_P} ${BDIR}

# VTYPE (unsigned int)
#
$(BDIR)/cnh: $(ODIR)/cnh.o $(ODIR)/utils.o $(ODIR)/graph_functions.o 
	$(CC) $(LDFLAGS) -o $@ $^

# VTYPE (unsigned long long) 
#
$(BDIR)/cnh_large: $(ODIR)/cnh.large.o $(ODIR)/utils.large.o $(ODIR)/graph_functions.large.o  
	$(CCC) $(FLAGS) -DLARGE -o  $@ $^ 

### CLEAN

clean:
	rm -vf $(ODIR)/*.o
	
distclean:
	rm -rf $(ODIR) $(BDIR) 

