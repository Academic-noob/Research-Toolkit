#!/bin/bash
#SBATCH --job-name=nonrad_wswq
#SBATCH --partition=cpu1_q
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36


cd ${SLURM_SUBMIT_DIR}

VASP="${VASP:-mpirun /opt/ohpc/pub/apps/vasp6.3.2/bin/vasp_std}"
#VASP="${VASP:-mpirun /opt/ohpc/pub/apps/vasp5.4.4_nonrad/bin/vasp_std}"
# export OMP_NUM_THREADS=1
# export MKL_NUM_THREADS=1
# export OMP_STACKSIZE=1024M
# ulimit -s unlimited

CASE_ROOT="${CASE_ROOT:-$(pwd)/nonrad_case}"

run_dir () {
    local d="$1"
    echo "===== RUN WSWQ: $d ====="
    cd "$d"
    $VASP > vasp.out
    cd "$SLURM_SUBMIT_DIR"
}

cd "$SLURM_SUBMIT_DIR"
#excited
for branch in ground ; do
    for d in "${CASE_ROOT}/04_wswq/${branch}"/*; do
        [[ -d "$d" ]] || continue
        run_dir "$d"
    done
done
