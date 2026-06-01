#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# USER-EDITABLE SECTION
# Two path modes are supported:
#   1. Current behavior: leave PROJECT_ROOT empty and run inside the defect root.
#   2. Internal path: fill PROJECT_ROOT with an absolute path and run from anywhere.
# A second positional argument overrides PROJECT_ROOT for one run.
# =============================================================================
PROJECT_ROOT=""                       # example: /home/student/works/ljd/nonrad/Int_I
GROUND_CHARGE="charge_state_0"        # relative to PROJECT_ROOT, or use an absolute path
EXCITED_CHARGE="charge_state_-1"      # relative to PROJECT_ROOT, or use an absolute path
CASE_DIR="nonrad_case"                # relative to PROJECT_ROOT, or use an absolute path
NPOINTS=9
QMIN="-0.5"
QMAX="0.5"
REF_INDEX=4                           # for 9 points: 0..8, index 4 is Q=0
PYTHON_BIN="${PYTHON_BIN:-python}"
# =============================================================================

mode="${1:-help}"
root_override="${2:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -n "$root_override" ]]; then
    ROOT="$root_override"
elif [[ -n "$PROJECT_ROOT" ]]; then
    ROOT="$PROJECT_ROOT"
else
    ROOT="$(pwd)"
fi
ROOT="$(cd "$ROOT" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR"

resolve_path() {
    local base="$1"
    local value="$2"
    if [[ "$value" = /* ]]; then
        printf '%s\n' "$value"
    else
        printf '%s/%s\n' "$base" "$value"
    fi
}

CASE_ROOT="$(resolve_path "$ROOT" "$CASE_DIR")"
GROUND_SRC="$(resolve_path "$ROOT" "$GROUND_CHARGE")"
EXCITED_SRC="$(resolve_path "$ROOT" "$EXCITED_CHARGE")"

die() { echo "ERROR: $*" >&2; exit 1; }
need_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }

read_nelect() {
    local f="$1"
    need_file "$f"
    local val
    val="$(grep -Eo '[-+]?[0-9]+([.][0-9]+)?' "$f" | head -1 || true)"
    [[ -n "$val" ]] || die "Cannot parse numeric NELECT from file: $f"
    echo "$val"
}

inject_nelect() {
    local src="$1"
    local dst="$2"
    local nelect="$3"
    need_file "$src"
    cp "$src" "$dst"
    if grep -qiE '^[[:space:]]*NELECT[[:space:]]*=' "$dst"; then
        sed -i -E "s/^[[:space:]]*NELECT[[:space:]]*=.*/NELECT = ${nelect}/I" "$dst"
    else
        printf "\nNELECT = %s\n" "$nelect" >> "$dst"
    fi
}

copy_common() {
    local charge_dir="$1"
    local target_dir="$2"
    local incar_name="$3"
    local kpoints_name="$4"
    local nelect
    nelect="$(read_nelect "${charge_dir}/Total_charge")"
    need_file "${charge_dir}/POTCAR"
    need_file "${charge_dir}/${incar_name}"
    need_file "${charge_dir}/${kpoints_name}"
    cp "${charge_dir}/POTCAR" "${target_dir}/POTCAR"
    cp "${charge_dir}/${kpoints_name}" "${target_dir}/KPOINTS"
    inject_nelect "${charge_dir}/${incar_name}" "${target_dir}/INCAR" "$nelect"
}

setup_relax() {
    echo "[MODE] setup endpoint relax inputs"
    need_file "${GROUND_SRC}/CONTCAR"
    need_file "${EXCITED_SRC}/CONTCAR"
    mkdir -p "${CASE_ROOT}/01_relax/ground" "${CASE_ROOT}/01_relax/excited"
    cp "${GROUND_SRC}/CONTCAR" "${CASE_ROOT}/01_relax/ground/POSCAR"
    cp "${EXCITED_SRC}/CONTCAR" "${CASE_ROOT}/01_relax/excited/POSCAR"
    copy_common "${GROUND_SRC}" "${CASE_ROOT}/01_relax/ground" "INCAR.relax" "KPOINTS.relax"
    copy_common "${EXCITED_SRC}" "${CASE_ROOT}/01_relax/excited" "INCAR.relax" "KPOINTS.relax"
    echo "[OK] Relax inputs created:"
    echo "  ${CASE_ROOT}/01_relax/ground"
    echo "  ${CASE_ROOT}/01_relax/excited"
}

setup_ccd() {
    echo "[MODE] generate CCD structures from relaxed endpoints"
    need_file "${CASE_ROOT}/01_relax/ground/CONTCAR"
    need_file "${CASE_ROOT}/01_relax/excited/CONTCAR"
    mkdir -p "${CASE_ROOT}/02_ccd_structures"
    "${PYTHON_BIN}" "${SCRIPTS_DIR}/01_generate_ccd_structures.py" \
        --ground-contcar "${CASE_ROOT}/01_relax/ground/CONTCAR" \
        --excited-contcar "${CASE_ROOT}/01_relax/excited/CONTCAR" \
        --outdir "${CASE_ROOT}/02_ccd_structures" \
        --npoints "${NPOINTS}" \
        --qmin "${QMIN}" \
        --qmax "${QMAX}"
    echo "[OK] CCD structures created in ${CASE_ROOT}/02_ccd_structures"
}

setup_static() {
    echo "[MODE] prepare CCD static inputs"
    need_file "${CASE_ROOT}/02_ccd_structures/ccd_points.tsv"
    for branch in ground excited; do
        src_charge="${GROUND_SRC}"
        [[ "$branch" == "excited" ]] && src_charge="${EXCITED_SRC}"
        for poscar in "${CASE_ROOT}/02_ccd_structures/${branch}"/*/POSCAR; do
            [[ -f "$poscar" ]] || continue
            idx="$(basename "$(dirname "$poscar")")"
            d="${CASE_ROOT}/03_static/${branch}/${idx}"
            mkdir -p "$d"
            cp "$poscar" "$d/POSCAR"
            copy_common "$src_charge" "$d" "INCAR.static" "KPOINTS.static"
            echo "[OK] static ${branch}/${idx}"
        done
    done
    echo "[OK] Static inputs are in ${CASE_ROOT}/03_static"
}

setup_wswq() {
    echo "[MODE] prepare WSWQ inputs"
    need_file "${CASE_ROOT}/02_ccd_structures/ccd_points.tsv"
    for branch in ground excited; do
        src_charge="${GROUND_SRC}"
        [[ "$branch" == "excited" ]] && src_charge="${EXCITED_SRC}"
        ref="${CASE_ROOT}/03_static/${branch}/${REF_INDEX}"
        need_file "${ref}/WAVECAR"
        need_file "${ref}/CHGCAR"
        need_file "${ref}/vasprun.xml"
        for static_dir in "${CASE_ROOT}/03_static/${branch}"/*; do
            [[ -d "$static_dir" ]] || continue
            idx="$(basename "$static_dir")"
            d="${CASE_ROOT}/04_wswq/${branch}/${idx}"
            mkdir -p "$d"
            need_file "${static_dir}/WAVECAR"
            cp "${ref}/POSCAR" "$d/POSCAR"
            cp "${src_charge}/POTCAR" "$d/POTCAR"
            cp "${src_charge}/KPOINTS.wswq" "$d/KPOINTS"
            inject_nelect "${src_charge}/INCAR.wswq" "$d/INCAR" "$(read_nelect "${src_charge}/Total_charge")"
            cp "${ref}/WAVECAR" "$d/WAVECAR"
            cp "${static_dir}/WAVECAR" "$d/WAVECAR.qqq"
            cp "${ref}/CHGCAR" "$d/CHGCAR"
            echo "[OK] wswq ${branch}/${idx}"
        done
    done
    echo "[OK] WSWQ inputs are in ${CASE_ROOT}/04_wswq"
    echo "[NOTE] initial_vasprun for Wif should be:"
    echo "  ${CASE_ROOT}/03_static/<branch>/${REF_INDEX}/vasprun.xml"
}

print_help() {
    cat <<EOF
Usage:
  bash scripts/00_setup_nonrad_case.sh relax
  bash scripts/00_setup_nonrad_case.sh ccd
  bash scripts/00_setup_nonrad_case.sh static
  bash scripts/00_setup_nonrad_case.sh wswq
  bash scripts/00_setup_nonrad_case.sh all_after_relax

Optional one-run project-root override:
  bash scripts/00_setup_nonrad_case.sh relax /absolute/path/to/defect_root

Resolved paths:
  ROOT        = ${ROOT}
  CASE_ROOT   = ${CASE_ROOT}
  GROUND_SRC  = ${GROUND_SRC}
  EXCITED_SRC = ${EXCITED_SRC}
EOF
}

case "$mode" in
    relax) setup_relax ;;
    ccd) setup_ccd ;;
    static) setup_static ;;
    wswq) setup_wswq ;;
    all_after_relax)
        setup_ccd
        setup_static
        ;;
    help|-h|--help) print_help ;;
    *) print_help; die "Unknown mode: $mode" ;;
esac
