#!/bin/bash -x

#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=00:10:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin "/home/hp120263/k00634/hodgkinhuxley/simple/hh.out ./"
#PJM -s

. /work/system/Env_base

export OMP_NUM_THREADS=1

MPIEXEC="mpiexec"

EXEC="./hh.out"

time  ${MPIEXEC} ${EXEC}

