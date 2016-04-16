all: lib

install:
	mkdir -p $(LIB)
	mv libtrdf.a $(LIB)

lib: cinterface.o TRDF.o functions.o
	ar crv libtrdf.a $^

%.o: %.f tr_params.par
	$(FC) $(FCC) -c $< -o $@ 

cinterface.o: cinterface.f90
	gfortran -c $< -o $@ 

clean:
	rm -vf *.o *.mod *~

.PHONY: lib clean all install
