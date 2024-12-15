
# run
# mpirun -np 10 ./mmBp.exe

# make

# Compiler and flags for NVidia HPC SDK OpenMPI GPU
FC = mpifort
FFLAGS = -O3 -cpp -acc  

# Source files
SRC = random.f90 normals.f90 sorts.f90 global_data.f90 timers.f90 functions.f90 mmBp.f90

# Object files (generated in the current directory)
OBJ = $(SRC:%.f90=%.o)

# Module files (generated during compilation)
MOD = $(SRC:%.f90=%.mod)

# Executable file (generated in the current directory)
EXEC = mmBp.exe

# Default target
all: $(EXEC)

# Rule to create the executable
$(EXEC): $(OBJ)
	$(FC) $(OBJ) -o $(EXEC) $(FFLAGS)  

# Rule to compile Fortran source files into object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# make clean
clean:
	rm -f $(OBJ) $(MOD) $(EXEC)

# Phony targets (not actual files)
.PHONY: all clean