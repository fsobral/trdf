install: $(SOLVER_INTERFACE).o
	mkdir -p $(OBJ)
	cp $< $(OBJ)/solver.o

%.o: %.f $(SRC)/tr_params.par
	$(FC) $(FCC) -I$(SRC) -c $< -o $@ 

clean:
	rm -vf *.o *.mod *~
