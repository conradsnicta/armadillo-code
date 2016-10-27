# You may need to edit this file to reflect the type and capabilities of your system.
# The defaults are for a Linux system and may need to be changed for other systems (eg. Mac OS X).


CXX=g++

#CXX=CC
## When using the Sun Studio compiler


ARMA_INCLUDE_FLAG = -I ../include
## If you've installed Armadillo's headers manually, you may need to tell the compiler where they are.
## For example, change ../include to /usr/local/include


LIB_FLAGS = -lblas -llapack 
#LIB_FLAGS = -lopenblas -llapack 
#LIB_FLAGS = -framework Accelerate
#LIB_FLAGS = -library=sunperf

## NOTE: on Ubuntu and Debian based systems you may need to add -lgfortran to LIB_FLAGS
## NOTE: if you're using Mac OS, use the line with -framework Accelerate 
## NOTE: if you're using the Sun Studio compiler, use the line with -library=sunperf


OPT = -O2
## As the Armadillo library uses recursive templates, compilation times depend on the level of optimisation:
##
## -O0: quick compilation, but the resulting program will be slow
## -O1: good trade-off between compilation time and execution speed
## -O2: produces programs which have almost all possible speedups, but compilation takes longer
## -O3: enables auto vectorisation when using gcc

#OPT = -xO4 -xannotate=no
## When using the Sun Studio compiler


#EXTRA_OPT = -fwhole-program
## Uncomment the above line if you're compiling all source files into one program in a single hit


#DEBUG = -DARMA_EXTRA_DEBUG
## Uncomment the above line to enable low-level debugging.
## Lots of debugging information will be printed when a compiled program is run.
## Please enable this option when reporting bugs.


#FINAL = -DARMA_NO_DEBUG
## Uncomment the above line to disable Armadillo's checks.
## Not recommended unless your code has been first thoroughly tested!


CXXFLAGS = $(ARMA_INCLUDE_FLAG) $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT) 

all: example1

example1: example1.cpp
	$(CXX) $(CXXFLAGS)  -o $@  $<  $(LIB_FLAGS)


.PHONY: clean

clean:
	rm -f example1

