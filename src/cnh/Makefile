SDIR = src
ODIR = obj
TDIR = test
IDIR = include
BDIR = bin

GCC_NAME    =  gcc
GCC_FLAGS   =  -pedantic  -Wall  -std=c99  -I$(IDIR)
GCC_FLAGS   += -g
GCC_FLAGS   += -O3

TARGET = cnh 
TARGET_LARGE = $(TARGET:=_large)
# CHANGE THESE LINES TO USE YOUR FAVORITE COMPILER
CCC      = $(GCC_NAME)
FLAGS    = $(GCC_FLAGS)

vpath %.c $(TDIR) 
vpath %.c $(SDIR)

.PHONY: compile clean install distclean 

#.SUFFIXES: .c

$(ODIR)/%.o: %.c
	mkdir -p $(ODIR)
	$(CCC) $(FLAGS) -o $@ -c $<

$(ODIR)/%.large.o: %.c 
	mkdir -p $(ODIR)
	$(CCC) $(FLAGS) -DLARGE -c $^ -o $@


### MAKE ALL
compile:  $(TARGET) $(TARGET_LARGE)  


# VTYPE (unsigned int)
#
cnh: $(ODIR)/cnh.o $(ODIR)/utils.o $(ODIR)/graph_functions.o $(IDIR)/define.h $(IDIR)/graph_functions.h $(IDIR)/timer.h 
	$(CCC) $(FLAGS) -o $@ $^ 

#
# VTYPE (unsigned long long) 
#
cnh_large: $(ODIR)/cnh.large.o $(ODIR)/utils.large.o $(ODIR)/graph_functions.large.o $(IDIR)/define.h $(IDIR)/graph_functions.h $(IDIR)/timer.h 
	$(CCC) $(FLAGS) -DLARGE -o  $@ $^ 

### CLEAN

clean:
	rm -vf $(ODIR)/*.o
	
distclean:
	rm -vf $(ODIR)/*.o $(TARGET) $(TARGET_LARGE) 


install:
	make compile
	mkdir -p ${BDIR}
	mv $(TARGET) $(TARGET_LARGE) ${BDIR}

