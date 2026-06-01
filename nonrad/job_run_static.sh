#!/bin/bash
#SBATCH --job-name=nonrad_static
#SBATCH --partition=cpu2_q
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40

set -euo pipefail

VASP="${VASP:-mpirun /opt/ohpc/pub/apps/vasp6.3.2/bin/vasp_std}"
CASE_ROOT="${CASE_ROOT:-$(pwd)/nonrad_case}"

run_dir () {
    local d="$1"
    echo "===== RUN STATIC: $d ====="
    cd "$d"
    $VASP > vasp.out 2>&1
    cd "$SLURM_SUBMIT_DIR"
}

cd "$SLURM_SUBMIT_DIR"

for branch in ground excited; do
    for d in "${CASE_ROOT}/03_static/${branch}"/*; do
        [[ -d "$d" ]] || continue
        run_dir "$d"
    done
done
