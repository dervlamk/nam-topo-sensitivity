#!/usr/bin/env bash
# Step 1 for zenodo_subset.sh: make `time` the unlimited (record) dimension of the FLOR time-series files.
# ncrcat's JAS subcycle selection only works on a record dimension. The lat/lon box is applied here too so the
# intermediate copies stay small; values are untouched (lossless deflate only).
#
# Everything else in IN_DIR (zsurf, land_mask, obs, MERRA-2) is symlinked into REC_DIR unchanged, so
# zenodo_subset.sh can then be run on REC_DIR alone.
#
# Usage: ./zenodo_make_time_record.sh IN_DIR REC_DIR      then:  ./zenodo_subset.sh REC_DIR OUT_DIR
# Env:   LATMIN LATMAX LONMIN LONMAX (defaults 5 70 200 300), LEVEL (deflate level, default 1), DRY_RUN=1

set -euo pipefail
IN_DIR=${1:?usage: $0 IN_DIR REC_DIR}; REC_DIR=${2:?usage: $0 IN_DIR REC_DIR}
LATMIN=${LATMIN:-5}; LATMAX=${LATMAX:-70}; LONMIN=${LONMIN:-200}; LONMAX=${LONMAX:-300}
LEVEL=${LEVEL:-1}; DRY_RUN=${DRY_RUN:-0}
command -v ncks >/dev/null || { echo "ncks not found; try: module load nco" >&2; exit 1; }
mkdir -p "$REC_DIR"

shopt -s nullglob
for f in "$IN_DIR"/cm2pt5-flor.*.monthly.nc; do
    b=$(basename "$f"); out="$REC_DIR/$b"
    if [[ -e $out ]]; then echo "skip (exists): $b"; continue; fi
    echo "$b"
    cmd=(ncks -O -4 -L "$LEVEL" --mk_rec_dmn time -d lat,"$LATMIN.","$LATMAX." -d lon,"$LONMIN.","$LONMAX." "$f" "$out.part")
    if [[ $DRY_RUN == 1 ]]; then echo "+ ${cmd[*]} && mv $out.part $out"; else "${cmd[@]}" && mv "$out.part" "$out"; fi
done

# symlink the non-time-series files so REC_DIR holds the full set
for f in "$IN_DIR"/*.nc; do
    b=$(basename "$f"); [[ $b == cm2pt5-flor.*.monthly.nc || -e $REC_DIR/$b ]] && continue
    if [[ $DRY_RUN == 1 ]]; then echo "+ ln -s $f $REC_DIR/$b"; else ln -s "$(cd "$IN_DIR" && pwd)/$b" "$REC_DIR/$b"; fi
done
