#!/bin/bash
# Submit every fully prepared SnB directory exactly once.
# Run from workflow root after reviewing 00_reports/finalize_vasp_inputs_audit.csv.
set -euo pipefail

ROOT="${1:-02_snb_structures}"

find "${ROOT}" -type f -name job.sh -print0 | while IFS= read -r -d '' job; do
  dir=$(dirname "${job}")
  if [[ -f "${dir}/.submitted" ]]; then
    echo "SKIP already submitted: ${dir}"
    continue
  fi
  [[ -s "${dir}/POSCAR" ]] || { echo "ERROR missing POSCAR: ${dir}"; exit 1; }
  [[ -s "${dir}/INCAR" ]] || { echo "ERROR missing INCAR: ${dir}"; exit 1; }
  [[ -s "${dir}/KPOINTS" ]] || { echo "ERROR missing KPOINTS: ${dir}"; exit 1; }
  [[ -s "${dir}/POTCAR" ]] || { echo "ERROR missing POTCAR: ${dir}"; exit 1; }
  (
    cd "${dir}"
    sbatch job.sh | tee .submitted
  )
done
