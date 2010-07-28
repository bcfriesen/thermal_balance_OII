#-------------------------------darkstar------------------------------------
#F77=gfortran
#FFLAGS=
#LFLAGS=-llapack -lblas

#-------------------------------duke----------------------------------------
#F77=xlf
#FFLAGS=
#LFLAGS=-llapack_atl_xlf -lptcblas_xlf -lptf77blas_xlf -latlas_xlf
#---------------------------------------------------------------------------

#-------------------------------mactel--------------------------------------
F77=ifort
FFLAGS=
LFLAGS=-lblas -llapack
#---------------------------------------------------------------------------

SOURCES=\
therm_bal.f

all: therm_bal

therm_bal: $(SOURCES)
	$(F77) $(SOURCES) $(FFLAGS) $(LFLAGS) -o $@

%.o: %.f
	$(F77) $(FFLAGS) $<

clean:
	rm -rf therm_bal *.out *.o
