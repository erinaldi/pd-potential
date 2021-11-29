FC=gfortran
# use this on Hokusai after "module load intel"
#FC=ifort
FCFLAGS=-O2
LDFLAGS=-llapack -lblas
DEPENDALL=staticparameters.f90

all: ek.exe

ek: EK_strong_coupling_v13.o
	${FC} ${FCFLAGS} EK_strong_coupling_v13.o -o $@ ${LDFLAGS}

%.o: %.f90 $(DEPENDALL)
	${FC} -c ${FCFLAGS} $<

mtmod.mod: mt19937.f90
	${FC} -c ${FCFLAGS} $<

clean:
	echo cleaning up in .
	$(RM) -f *.o *.mod
