#!/usr/bin/env bash
set -euo pipefail
OUT="$(./build/bin/iv_export_and_qc -f data/market_samples/sample_market.csv)"
echo "$OUT"

# Calendar must be 0 violations
CAL_VIO=$(echo "$OUT" | sed -n 's/.*Calendar QC: checks=[0-9]\+ violations=\([0-9]\+\).*/\1/p')
[ -z "$CAL_VIO" ] && { echo "Parse error (calendar)"; exit 2; }
[ "$CAL_VIO" -eq 0 ] || { echo "Calendar violations: $CAL_VIO"; exit 3; }

# Mean rel.err (arb-OK) <= 0.02
MEAN=$(echo "$OUT" | sed -n 's/.*QC repricing (arb-OK only).*mean=\([0-9.]\+\).*/\1/p')
[ -z "$MEAN" ] && { echo "Parse error (mean)"; exit 2; }
awk "BEGIN {exit !($MEAN <= 0.02)}" || { echo "Mean rel.err too high: $MEAN"; exit 4; }

echo "IV surface QC thresholds: OK"
