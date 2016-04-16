install:
	mkdir -p $(LIB)
	mv libhs.a $(LIB)

lib: hs_problems/prob.o hs_problems/conv.o
	ar cr libhs.a $^

%.o: %.f
	gfortran -c $< -o $@

clean:
	rm -vf *.o *.mod *.a *~ hs_problems/*.o

.PHONY: lib clean
