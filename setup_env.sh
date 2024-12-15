#!/bin/bash

# Load the necessary module
module load nvhpc/21.3

# Set the environment variables
export PATH=$PATH:$MODULE_NVHPC_PREFIX/Linux_ppc64le/21.3/comm_libs/mpi/bin
export CPATH=$CPATH:$MODULE_NVHPC_PREFIX/Linux_ppc64le/21.3/comm_libs/mpi/include
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MODULE_NVHPC_PREFIX/Linux_ppc64le/21.3/comm_libs/mpi/lib

echo "Environment variables set successfully"