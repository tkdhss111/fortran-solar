#!/usr/bin/env bash
# Engerer2 vs Erbs (-> Perez POA) A/B over LFM GPV GHI.
#   usage: run_ab.sh '<parquet-glob>' [tilt_deg=30] [panel_azimuth_deg=180] [albedo=0.2]
# Library is zero-dep; this harness uses duckdb only to read the LFM parquet -> CSV.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
GLOB="${1:?parquet glob required}"
TILT="${2:-30}" ; AZI="${3:-180}" ; ALB="${4:-0.2}"
BIN=/tmp/ab_poa ; CSV=/tmp/lfm_ghi.csv ; OUT=/tmp/ab_poa_out.csv

echo "[1/4] build ab_poa (+ library)" >&2
gfortran -O2 -ffree-line-length-none -o "$BIN" "$ROOT/src/solar_geometry_mo.f90" "$HERE/ab_poa.f90"

echo "[2/4] LFM GHI -> CSV (daytime only)" >&2
duckdb -csv -c "
  SELECT site_id AS site, jst, lat, lon, solar_radiation_wm2 AS ghi
  FROM read_parquet('$GLOB', union_by_name=true)
  WHERE solar_radiation_wm2 IS NOT NULL AND solar_radiation_wm2 > 0
  ORDER BY jst, site_id" > "$CSV"

echo "[3/4] run A/B (tilt=$TILT az=$AZI albedo=$ALB)" >&2
"$BIN" "$TILT" "$AZI" "$ALB" < "$CSV" > "$OUT"

echo "[4/4] breakdown by GHI level" >&2
duckdb -c "
  WITH d AS (SELECT *, poa_engerer2 - poa_erbs AS df FROM read_csv_auto('$OUT'))
  SELECT CASE WHEN ghi<200 THEN '1: <200' WHEN ghi<500 THEN '2: 200-500'
              WHEN ghi<800 THEN '3: 500-800' ELSE '4: >800' END AS ghi_bin,
         count(*) n, round(avg(ghi),0) ghi, round(avg(poa_engerer2),0) poa_e2,
         round(avg(poa_erbs),0) poa_erbs, round(avg(df),1) mean_diff,
         round(avg(abs(df)),1) mean_absdiff, round(max(abs(df)),1) max_absdiff
  FROM d GROUP BY 1 ORDER BY 1;"
echo "per-row CSV: $OUT" >&2
