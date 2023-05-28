COMPILER= gfortran

FLAGS = 

EXEC = derivatives

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

analytic_functions.o: analytic_functions.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

euler_formulas.o : euler_formulas.f90 analytic_functions.o types.o
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 euler_formulas.o analytic_functions.o types.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 types.o read_write.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
