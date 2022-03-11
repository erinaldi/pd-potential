FC=gfortran
# use this on Hokusai after "module load intel" and use -qopenmp
#FC=ifort
FCFLAGS=-O2 -fopenmp
LDFLAGS=-llapack -lblas
DEPENDALL=size_EK_v13.inc

all: ek

ek: EK_strong_coupling_v13.o
	${FC} ${FCFLAGS} EK_strong_coupling_v13.o -o $@ ${LDFLAGS}

ekomp: EK_strong_coupling_v14.o
	${FC} ${FCFLAGS} EK_strong_coupling_v14.o -o $@ ${LDFLAGS}

%.o: %.f90 $(DEPENDALL)
	${FC} -c ${FCFLAGS} $<

mtmod.mod: mt19937.f90
	${FC} -c ${FCFLAGS} $<

clean:
	echo cleaning up in .
	$(RM) -f *.o *.mod
