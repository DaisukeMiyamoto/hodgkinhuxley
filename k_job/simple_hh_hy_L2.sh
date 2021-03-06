#!/bin/bash -x

#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=00:10:00"
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"
#PJM --stgin "/home/hp120263/k00634/hodgkinhuxley/simple/hh.out ./"
#PJM --stgout "rank=* %r:./prof_cache/* /data/hp120263/k00634/result/prof/simple/%j_cache/"
#PJM --stgout "rank=* %r:./prof_perf/* /data/hp120263/k00634/result/prof/simple/%j_perf/"
#PJM --stgout "rank=* %r:./prof_stat/* /data/hp120263/k00634/result/prof/simple/%j_stat/"
#PJM --stgout "rank=* %r:./prof_inst/* /data/hp120263/k00634/result/prof/simple/%j_inst/"
#PJM --stgout "rank=* %r:./prof_mem/* /data/hp120263/k00634/result/prof/simple/%j_mem/"
#PJM -s

. /work/system/Env_base

export OMP_NUM_THREADS=8

PROF_CACHE="fapp -C -d ./prof_cache -Ihwm -Hevent=Cache -L2"
PROF_PERF="fapp -C -d ./prof_perf -Ihwm -Hevent=Performance -L2"
PROF_STAT="fapp -C -d ./prof_stat -Ihwm -Hevent=Statistics -L2"
PROF_INST="fapp -C -d ./prof_inst -Ihwm -Hevent=Instructions -L2"
PROF_MEM="fapp -C -d ./prof_mem -Ihwm -Hevent=MEM_access -L2"

MPIEXEC="mpiexec"

EXEC="./hh.out"

time ${PROF_CACHE} ${MPIEXEC} ${EXEC}
time ${PROF_PERF} ${MPIEXEC} ${EXEC}
time ${PROF_STAT} ${MPIEXEC} ${EXEC}
time ${PROF_INST} ${MPIEXEC} ${EXEC}
time ${PROF_MEM} ${MPIEXEC} ${EXEC}


