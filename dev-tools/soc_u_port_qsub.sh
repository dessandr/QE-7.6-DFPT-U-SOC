#!/usr/bin/env bash
# Submit a port validation batch script to a PBS compute node.
# Usage: bash dev-tools/soc_u_port_qsub.sh BATCH_SCRIPT [JOB_NAME]
#   SOC_U_QUEUE   (default RTX3090)   SOC_U_HOST (default pcs-gpu06)
#   SOC_U_NCPUS   (default 32)        SOC_U_WALLTIME (default 06:00:00)
#   SOC_U_NP / SOC_U_NK are exported to the batch (defaults 32 / 4).
# The job script and its output go next to BATCH_SCRIPT; --dry-run only
# writes the job script.
set -euo pipefail

dry=0
[[ ${1:-} == --dry-run ]] && { dry=1; shift; }
[[ $# -ge 1 ]] || { echo "Usage: bash $0 [--dry-run] BATCH_SCRIPT [JOB_NAME]" >&2; exit 2; }
batch=$(readlink -f "$1"); name=${2:-$(basename "${batch%.sh}")}
dir=$(dirname "$batch")
queue=${SOC_U_QUEUE:-RTX3090}; host=${SOC_U_HOST:-pcs-gpu06}
ncpus=${SOC_U_NCPUS:-32}; wall=${SOC_U_WALLTIME:-06:00:00}
np=${SOC_U_NP:-32}; nk=${SOC_U_NK:-4}
job=$dir/$name.qjob

cat > "$job" <<EOF
#PBS -N $name
#PBS -q $queue
#PBS -l select=1:ncpus=$ncpus:mpiprocs=$ncpus:mem=64gb:host=$host
#PBS -l walltime=$wall
#PBS -r n
#PBS -j oe
#PBS -o $dir/$name.o
cd "$dir"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset QE_SOC_U_DFPT_FIRST_ITER_THRESH QE_SOC_U_DFPT_ALL_ITER_THRESH
export I_MPI_PIN=1 I_MPI_PIN_DOMAIN=core
export PATH=/home/app/intel/oneapi/mpi/2021.18/bin:/home/app/intel/oneapi/compiler/2026.1/bin:\$PATH
export LD_LIBRARY_PATH=/home/grad/jmkang/lib/libxc/lib64:/home/app/intel/oneapi/mkl/2026.1/lib:/home/app/intel/oneapi/mpi/2021.18/lib:/home/app/intel/oneapi/compiler/2026.1/lib:/home/app/hdf5/2.2.0/intel_oneapi-2026.1.0/impi/lib64\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}
export SOC_U_NP=$np SOC_U_NK=$nk
echo "host \$(hostname)  start \$(date)"
bash "$batch" > "$dir/$name.log" 2>&1
echo "exit \$?  end \$(date)"
EOF
echo "job script: $job"
if [[ $dry == 0 ]]; then qsub "$job"; fi
