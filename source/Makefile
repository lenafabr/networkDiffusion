OBJS = keys.o readkey.o genutil.o mt19937.o inputparams.o networkutil.o diffutil.o testmod.o

FC = gfortran
FLAGS = -O3 
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -fcheck=all
LIBS =  -lblas -llapack

# executable name
EXE = ../networkdiffusion.exe
TESTEXE = ../testrun.exe

.SUFFIXES: .o .f90 

# general suffix rules
%.o: %.f90 
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

%.o: tools/%.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

clean: 
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
readkey.o: inputparams.o keys.o genutil.o mt19937.o
genutil.o: mt19937.o
