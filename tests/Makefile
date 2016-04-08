install:

hs:
	$(MAKE) -C hs all

%:
	$(MAKE) -C examples $@.o install

clean:
	rm -vf *~
	$(foreach i, $(shell ls -d */), $(MAKE) -C $(i) clean;)

# all: lib

# install:
# 	mkdir -p $(LIB)
# 	mv libtrdf.a $(LIB)

# lib: cinterface.o TRDF.o functions.o
# 	ar crv libtrdf.a $^

# %.o: %.f tr_params.par
# 	$(FC) $(FCC) -c $< -o $@ 

# cinterface.o: cinterface.f90
# 	gfortran -c $< -o $@ 

# clean:
# 	rm -vf *.o *.mod *~

# libtrdf.a: cinterface.o $(SOLVER).o TRDF.o functions.o
# 	ar crv libtrdf.a $^

# libhs.a: hs_problems/prob.o hs_problems/conv.o
# 	ar cr $@ $^

# %.o: %.f tr_params.par
# 	gfortran -c $(FCC) $< -o $@ 

# hs_problems/prob.o: hs_problems/prob.f
# 	gfortran -c $< -o $@

# hs_problems/conv.o: hs_problems/conv.f
# 	gfortran -c $< -o $@

# cinterface.o: cinterface.f90
# 	gfortran -c $< -o $@ 

# clean:
# 	rm -vf *.o *.mod *~ libtrdf.a libhs.a trdf c_trdf hstests \
# 	hs_problems/*.o

