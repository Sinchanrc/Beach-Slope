#!/usr/bin/bash

#EXPORTING OPENMP ENVIRONMENT VARIABLES
export OMP_NUM_THREADS=12
export OMP_SCHEDULE="dynamic"
# export OMP_STACKSIZE=16M
# export OMP_PROC_BIND="close"


# -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all \
# -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan \
# -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all \
# -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan \-fno-trapping-math -fno-signed-zeros 

#COMPILING SOURCE FILES#

#GNU COMPILERS
gfortran-12 -ffree-line-length-512 -Ofast -fopenmp -march=native -mtune=native -std=f2008 -c -flto=auto \
-fimplicit-none -funroll-loops -frepack-arrays -fno-realloc-lhs -mavx2 \
particle.f90 interactions.f90 transfer.f90 domain.f90 initialize.f90 kernel.f90 functions.f90 \
search.f90 gradient.f90 setup.f90 solver.f90 boundary.f90 output.f90 internal.f90 \
memory.f90 isph.f90 integrator.f90 turbulence.f90 part_shift.f90 mls.f90 \
scalar.f90 probe.f90

# ifort  -Ofast -qopenmp -std03 -ipo -xHost -c  \
# particle.f90 interactions.f90 transfer.f90 domain.f90 initialize.f90 kernel.f90 functions.f90 \
# search.f90 gradient.f90 setup.f90 solver.f90 boundary.f90 output.f90 internal.f90 \
# memory.f90 isph.f90 integrator.f90 turbulence.f90 part_shift.f90 mls.f90 \
# scalar.f90 probe.f90

#COMPILING EXECUTABLE#

#GNU COMPILERS
gfortran-12 -ffree-line-length-512 -Ofast -fopenmp -march=native -mtune=native -std=f2008 -flto=auto \
-fimplicit-none -funroll-loops -frepack-arrays -fno-realloc-lhs -mavx2 \
particle.o interactions.o transfer.o domain.o initialize.o kernel.o functions.o \
search.o gradient.o setup.o solver.o boundary.o output.o internal.o \
memory.o isph.o integrator.o turbulence.o part_shift.o scalar.o mls.o \
probe.o main.f90 -o main

# ifort  -Ofast -qopenmp -std03 -ipo -xHost \
# particle.o interactions.o transfer.o domain.o initialize.o kernel.o functions.o \
# search.o gradient.o setup.o solver.o boundary.o output.o internal.o \
# memory.o isph.o integrator.o turbulence.o part_shift.o scalar.o mls.o \
# probe.o main.f90 -o main
