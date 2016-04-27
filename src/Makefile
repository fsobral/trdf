all: lib

install:
	mkdir -p $(LIB)
	mv libtrdf.a $(LIB)

lib: cinterface.o trdf.o TRDF.o
	ar crv libtrdf.a $^

%.o: %.f tr_params.par
	$(FC) $(FCC) -c $< -o $@ 

trdf.o: trdf.f90
	gfortran -c $< -o $@ 

cinterface.o: cinterface.f90
	gfortran -c $< -o $@ 

clean:
	rm -vf *.o *.mod *~

.PHONY: lib clean all install
