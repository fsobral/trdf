install:
	mkdir -p $(OBJ)
	mv $(USERMAIN).o $(OBJ)/prob.o

%.o: %.f $(SRC)/tr_params.par
	$(FC) $(FCC) -I$(SRC) -c $< -o $@

clean:
	rm -vf *.o *.mod *~
