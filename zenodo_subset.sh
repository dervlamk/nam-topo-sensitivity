#!/usr/bin/env bash
# Subset + compress the Zenodo files with NCO.
#
#   FLOR monthly time series (cm2pt5-flor.*.monthly.nc):
#       domain subset, JAS-only, optional year range, netCDF4 deflate (+ optional bit rounding)
#   everything else (zsurf, land_mask, obs, MERRA-2):
#       netCDF4 deflate only (no spatial/time subsetting; grids and lon conventions differ)
#
# Usage:   ./zenodo_subset.sh IN_DIR OUT_DIR
# Env:     DRY_RUN=1            print the NCO commands, run nothing
#          LATMIN LATMAX        default 5 70
#          LONMIN LONMAX        default 200 300 (FLOR lon is 0-360)
#          YR_MIN YR_MAX        optional model years as integers, e.g. YR_MIN=101 YR_MAX=200
#          LEVEL                deflate level 1-9, default 5
#          PPC                  bit-rounding spec for --ppc, default "default=5"; set PPC= for lossless
#
# Assumes FLOR files are monthly, contiguous, start in January, and use dims named lat/lon/time.
# FLOR files that do not start in January are skipped with a message.
# JAS is selected with ncrcat's subcycle option (ncks does not support it), which needs time to be the
# record dimension; if it is not, a temporary copy with time as the record dimension is made first
# (already domain-subset, so small).

set -euo pipefail

IN_DIR=${1:?usage: $0 IN_DIR OUT_DIR}
OUT_DIR=${2:?usage: $0 IN_DIR OUT_DIR}
LATMIN=${LATMIN:-5};   LATMAX=${LATMAX:-70}
LONMIN=${LONMIN:-200}; LONMAX=${LONMAX:-300}
YR_MIN=${YR_MIN:-};    YR_MAX=${YR_MAX:-}
LEVEL=${LEVEL:-5}
PPC=${PPC-default=5}
DRY_RUN=${DRY_RUN:-0}

for c in ncks ncrcat ncdump; do
    command -v "$c" >/dev/null || { echo "$c not found (NCO / netCDF); try: module load nco" >&2; exit 1; }
done
mkdir -p "$OUT_DIR"

run() { if [[ $DRY_RUN == 1 ]]; then echo "+ $*"; else "$@"; fi; }
ppc_args=(); [[ -n $PPC ]] && ppc_args=(--ppc "$PPC")
tmp=""; trap '[[ -n $tmp ]] && rm -f "$tmp"' EXIT

# echo "YEAR MONTH NRECORDS" of the time axis (ncdump -t prints calendar dates for any CF calendar)
first_date() {
    local d
    d=$(ncdump -t -v time "$1" | sed -n '/^ *time *=/,$p') || return 1
    [[ $d =~ \"0*([0-9]+)-([0-9]{1,2})- ]] || return 1
    echo "${BASH_REMATCH[1]} $((10#${BASH_REMATCH[2]})) $(( $(tr -cd '"' <<<"$d" | wc -c) / 2 ))"
}

fsize() { stat -c %s "$1" 2>/dev/null || stat -f %z "$1"; }

tot_in=0; tot_out=0
shopt -s nullglob
for f in "$IN_DIR"/*.nc; do
    b=$(basename "$f"); out="$OUT_DIR/$b"
    if [[ -e $out ]]; then echo "skip (exists): $b"; continue; fi

    if [[ $b == cm2pt5-flor.*.monthly.nc ]]; then
        read -r y0 m0 nrec < <(first_date "$f") || { echo "SKIP $b: cannot read time axis"; continue; }
        if [[ $m0 != 1 ]]; then echo "SKIP $b: first month is $m0, not January (handle by hand)"; continue; fi
        ymin=${YR_MIN:-$y0}
        start=$(( (ymin - y0) * 12 + 6 ))                    # July of ymin (0-based record index)
        if [[ -n $YR_MAX ]]; then end=$(( (YR_MAX - y0) * 12 + 11 )); else end=$(( nrec - 1 )); fi
        if (( start < 0 || end < start || end >= nrec )); then
            echo "SKIP $b: year range outside file (first year $y0, $nrec records)"; continue
        fi
        echo "FLOR  $b  (first year $y0, $nrec records -> records $start..$end, JAS only)"
        dom=(-d lat,"$LATMIN.","$LATMAX." -d lon,"$LONMIN.","$LONMAX.")
        src=$f
        if ! ncdump -h "$f" | grep -q 'time = UNLIMITED'; then
            tmp=$(mktemp "$OUT_DIR/.tmp.XXXXXX"); mv "$tmp" "$tmp.nc"; tmp="$tmp.nc"
            run ncks -O --mk_rec_dmn time ${dom[@]+"${dom[@]}"} "$f" "$tmp"
            src=$tmp; dom=()                                   # domain already applied
        fi
        run ncrcat -O -4 -L "$LEVEL" ${ppc_args[@]+"${ppc_args[@]}"} ${dom[@]+"${dom[@]}"} \
            --mro -d time,"$start","$end",12,3 "$src" "$out"
        [[ -n $tmp ]] && { rm -f "$tmp"; tmp=""; }
    else
        echo "other $b  (compress only)"
        run ncks -O -4 -L "$LEVEL" ${ppc_args[@]+"${ppc_args[@]}"} "$f" "$out"
    fi

    if [[ $DRY_RUN != 1 && -e $out ]]; then
        si=$(fsize "$f"); so=$(fsize "$out")
        tot_in=$((tot_in + si)); tot_out=$((tot_out + so))
        printf '      %.2f GB -> %.2f GB\n' "$(bc -l <<<"$si/1e9")" "$(bc -l <<<"$so/1e9")"
    fi
done

[[ $DRY_RUN == 1 ]] || printf 'TOTAL  in %.1f GB  out %.1f GB\n' "$(bc -l <<<"$tot_in/1e9")" "$(bc -l <<<"$tot_out/1e9")"
