FC=ifort
FCFLAGS=-mkl

SOURCES=\
therm_bal.f90

all: therm_bal

therm_bal: $(SOURCES)
	$(FC) $(SOURCES) $(FCFLAGS) $(LDFLAGS) -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) $<

clean:
	rm -rf therm_bal *.out *.o
