#!/bin/bash -x
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=00:10:00"
#PJM --rsc-list "rscgrp=small"
#PJM --mpi "proc=8"
#PJM -s

# staging
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"

# CHANGE TO YOUR OWN DIR
#PJM --stgin-basedir /home/hp120263/k00634

#PJM --stgout "rank=* %r:./prof_cache/* /data/hp120263/k00634/result/hh/%j_cache/"
#PJM --stgout "rank=* %r:./prof_perf/* /data/hp120263/k00634/result/hh/%j_perf/"
#PJM --stgout "rank=* %r:./prof_stat/* /data/hp120263/k00634/result/hh/%j_stat/"
#PJM --stgout "rank=* %r:./prof_inst/* /data/hp120263/k00634/result/hh/%j_inst/"
#PJM --stgout "rank=* %r:./prof_mem/* /data/hp120263/k00634/result/hh/%j_mem/"

#PJM --stgin "rank=* ./neuron_kplus/stgin/* %r:./"
#PJM --stgin "rank=* ./neuron_kplus/specials/sparc64/special %r:./"

#PJM --stgin "rank=* ./hodgkinhuxley/hoc/* %r:./"


. /work/system/Env_base

#export OMP_NUM_THREADS=8

NRNIV="./special"
HOC_NAME="./hh_k.hoc"
NRNOPT=""

LPG="lpgparm -t 4MB -s 4MB -d 4MB -h 4MB -p 4MB"
#MPIEXEC="mpiexec -mca mpi_print_stats 1"
MPIEXEX=""
#MPIEXEC="mpiexec -mca mpi_print_stats 2 -mca mpi_print_stats_ranks 0"

PROF_CACHE="fapp -C -d ./prof_cache -Ihwm -Hevent=Cache -L1"
PROF_PERF="fapp -C -d ./prof_perf -Ihwm -Hevent=Performance -L1"
PROF_STAT="fapp -C -d ./prof_stat -Ihwm -Hevent=Statistics -L1"
PROF_INST="fapp -C -d ./prof_inst -Ihwm -Hevent=Instructions -L1"
PROF_MEM="fapp -C -d ./prof_mem -Ihwm -Hevent=MEM_access -L1"

#echo "${PROF} ${MPIEXEC} ${LPG} ${NRNIV} ${NRNOPT} ${HOC_NAME}"
#time ${PROF} ${MPIEXEC} ${LPG} ${NRNIV} ${NRNOPT} ${HOC_NAME}
time ${PROF_CACHE} ${MPIEXEC} ${LPG} ${NRNIV} ${NRNOPT} ${HOC_NAME}
time ${PROF_PERF}  ${MPIEXEC} ${LPG} ${NRNIV} ${NRNOPT} ${HOC_NAME}
time ${PROF_STAT}  ${MPIEXEC} ${LPG} ${NRNIV} ${NRNOPT} ${HOC_NAME}
time ${PROF_INST}  ${MPIEXEC} ${LPG} ${NRNIV} ${NRNOPT} ${HOC_NAME}
time ${PROF_MEM}   ${MPIEXEC} ${LPG} ${NRNIV} ${NRNOPT} ${HOC_NAME}

sync


