# mmBp

This repository contains the implementation code of the paper "[Multiscale Bernstein Copula with Bayesian Model Averaging](https://www.economics.utoronto.ca/mburda/papers/mmBp.pdf)" by [Martin Burda](https://www.economics.utoronto.ca/mburda/) and [Artem Prokhorov](https://sites.google.com/site/artembprokhorov) (2024).

The code is written in [Modern Fortran](https://fortran-lang.org/), designed for scalable native [parallelization](https://wvuhpc.github.io/Modern-Fortran/20-Parallel-Programming/index.html) for distributed nodes of Unix High Performance Computing ([HPC](https://scinethpc.ca/)) systems with across-node communication via Message Passing Interface ([MPI](https://github.com/open-mpi/ompi)) and within-node acceleration on Graphical Processing Units ([GPU](https://www.nvidia.com/en-sg/data-center/hpc/)) with [OpenACC](https://github.com/OpenACC) offload directives.


## Code Files

1. *mmBp.f90* - the main program
2. *global_data.f90* - routines for data allocation, deallocation, offload
3. *functions.f90* - implementation functions and routines
4. *timers.f90* - routines for wallclock, CPU, and MPI run time

## External Dependencies
Additionally, the code also uses the following external dependencies:

1. Module `random` for generating random draws from the Beta and Gaussian densities, available [here](https://www.netlib.org/random/random.f90), save to file *random.f90* in the project directory.

2. Functions for efficient evaluation of :
- `normal_01_cdf` - standard Gaussian cdf
- `normal_cdf_inv` - inverse Gaussian cdf
- `normal_01_cdf_inv` - inverse standard Gaussian cdf
- `r8poly_value_horner` - a polynomial using Horner's method   
<br>The functions are available [here](https://people.math.sc.edu/Burkardt/f_src/prob/prob.f90), save in module `normals` in file *normals.f90* in the project directory. 

3. Module `sorts` containing sorting routines available [here](
https://www.mjr19.org.uk/IT/sorts/sorts.f90), save in file *sorts.f90* in the project directory.

## Data Processing
1. *_Data download.Rmd* - code for data download from Yahoo Finance
2. *_Plots.Rmd* - code for generating the figures used in the paper

## Compilation
If using a Unix/Linux HPC system, load the necessary modules and initialize environmental variables. For example, the [Mist](https://docs.scinet.utoronto.ca/index.php/Mist) cluster at [Scinet](https://scinethpc.ca/) is set up with the bash shell instructions in the file *setup_env.sh*, invoked by:

<span style="font-family: monospace;">bash setup_env.sh</span>

Compiler directives with flags for for NVidia HPC SDK are provided in *Makefile*, invoked by:

<span style="font-family: monospace;">make</span>

Object, module, and executable files are removed with:

<span style="font-family: monospace;">make clean</span>
