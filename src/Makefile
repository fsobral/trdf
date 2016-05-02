all: lib

install:
	mkdir -p $(LIB)
	mv libtrdf.a $(LIB)

lib: cinterface.o trdf.o TRDF.o
	ar crv libtrdf.a $^

trdf.o: trdf.f90
	gfortran -c $< -o $@ 

cinterface.o: cinterface.f90
	gfortran -c $< -o $@ 

clean:
	rm -vf *.o *.mod *~

.PHONY: lib clean all install
