#!/usr/bin/env bash
set -euo pipefail

# USER-EDITABLE INTERNAL PATH
# Keep "nonrad_case" for the original relative-path usage, or fill an absolute path.
CASE_ROOT_INTERNAL="nonrad_case"

# Optional positional argument overrides the internal path for one run.
CASE_ROOT="${1:-${CASE_ROOT_INTERNAL}}"

echo "===== CASE ROOT ====="
echo "$CASE_ROOT"

echo "===== OUTCAR identity ====="
find "$CASE_ROOT" -type f -name OUTCAR | while read -r outcar; do
    d="$(dirname "$outcar")"
    echo "---- $d ----"
    grep -E "NKPTS|NBANDS|NELECT|ISTART|ICHARG|IALGO|LHFCALC|AEXX|HFSCREEN|LWSWQ|free  energy   TOTEN" "$outcar" | head -40 || true
    ls -lh "$d"/WAVECAR "$d"/WAVECAR.qqq "$d"/CHGCAR "$d"/WSWQ "$d"/vasprun.xml 2>/dev/null || true
done

echo "===== possible bad zero-size files ====="
find "$CASE_ROOT" -type f \( -name WAVECAR -o -name CHGCAR -o -name WSWQ -o -name WAVECAR.qqq \) -size 0 -print
