#!/usr/bin/env bash
set -euo pipefail

BIN="build/bin/greeks_runner"
N=400000
BATCH=100000
SEED=42
TOL=-1
STEPS=1
ANTI=1
CRN=1

ZERO_SE_ABS=1e-12
DIFF_ABS_PASS=1e-8
DIFF_REL_PASS=1e-9

CASES=(
  "ATM_1Y|call|100|100|0.02|0.00|0.20|1.00"
  "OTM_3M|call|100|110|0.02|0.00|0.20|0.25"
  "ITM_2Y_PUT|put|100|120|0.02|0.00|0.20|2.00"
)

# --------- opts ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin) BIN="$2"; shift 2;;
    --n) N="$2"; shift 2;;
    --batch) BATCH="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --tol) TOL="$2"; shift 2;;
    --steps) STEPS="$2"; shift 2;;
    --no-anti) ANTI=0; shift;;
    --no-crn) CRN=0; shift;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done

if [[ ! -x "$BIN" ]]; then echo "ERROR: $BIN not found/executable"; exit 1; fi
ANTI_FLAG=$([[ $ANTI -eq 1 ]] && echo "--antithetic" || echo "")
CRN_FLAG=$([[ $CRN -eq 1 ]] && echo "--crn" || echo "--no-crn")

parse_field(){ echo "$1" | sed -E 's/^[^:]+:[[:space:]]*//'; }

run_greek(){
  local greek="$1" type="$2" S0="$3" K="$4" r="$5" q="$6" sigma="$7" T="$8"
  shift 8
  local PUT=""; [[ "$type" == "put" ]] && PUT="--put"
  local out
  out=$("$BIN" "$S0" "$K" "$r" "$q" "$sigma" "$T" "$SEED" "$N" "$BATCH" "$TOL" "$STEPS" \
        --greek "$greek" $PUT $CRN_FLAG $ANTI_FLAG "$@")
  local est se
  est=$(parse_field "$(echo "$out" | grep -E '^estimate')"); se=$(parse_field "$(echo "$out" | grep -E '^std_error')")
  echo "$est $se"
}

bs_vega(){
  python3 - "$@" <<'PY'
import sys, math
typ, S0, K, r, q, sigma, T = sys.argv[1], *map(float, sys.argv[2:])
SQRT2PI = (2.0*math.pi)**0.5
def Phi(x): return 0.5*math.erfc(-x/(2.0**0.5))
def pdf(x): return math.exp(-0.5*x*x)/SQRT2PI
if T<=0 or sigma<=0:
    print("0.000000000000"); sys.exit(0)
d1 = (math.log(S0/K) + (r - q + 0.5*sigma*sigma)*T)/(sigma*math.sqrt(T))
vega = S0*math.exp(-q*T)*math.sqrt(T)*pdf(d1)
print(f"{vega:.12f}")
PY
}

z_and_pass(){
  awk -v e="$1" -v s="$2" -v a="$3" -v zthr="$ZERO_SE_ABS" -v absthr="$DIFF_ABS_PASS" -v relthr="$DIFF_REL_PASS" '
    function abs(x){return x<0?-x:x}
    BEGIN{
      if (s<zthr){
        diff=abs(e-a); thr=absthr; absa=abs(a); if (relthr*absa>thr) thr=relthr*absa;
        if (diff<=thr) {print "0.000 YES"; exit} else {print "inf NO"; exit}
      } else {
        z=(e-a)/s; az=abs(z); pass=(az<=2.0)?"YES":"NO"; printf("%.3f %s\n", z, pass);
      }
    }'
}

halfwidth95(){ awk -v se="$1" 'BEGIN{z=1.959963984540054; printf("%.12f", z*se)}'; }
overlap_ci(){ awk -v l1="$2" -v h1="$3" -v l2="$5" -v h2="$6" 'BEGIN{print ((h1>=l2 && h2>=l1)?"YES":"NO")}'; }
hr(){ printf '%s\n' "-----------------------------------------------------------------------------------"; }

printf "\n== Vega LRM Bench (N=%d, batch=%d, steps=%d, seed=%d, anti=%s) ==\n\n" "$N" "$BATCH" "$STEPS" "$SEED" "$ANTI"

# 1) LRM vs BRV vs BS
echo "==== [1] Vega: LRM vs BRV vs BS ===="
for case in "${CASES[@]}"; do
  IFS='|' read -r NAME TYPE S0 K R Q SIGMA T <<<"$case"
  echo
  printf "Case: %s (type=%s, S0=%g K=%g r=%g q=%g sigma=%g T=%g)\n" \
    "$NAME" "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T"
  hr
  printf "%-8s %-10s %-14s %-14s %-12s %-6s\n" "Greek" "Mode" "MC_estimate" "MC_SE" "z_vs_BS" "PASS"

  ANA="$(bs_vega "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"

  read -r EST_LRM SE_LRM <<<"$(run_greek vega_lrm "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
  read -r Z_LRM PASS_LRM <<<"$(z_and_pass "$EST_LRM" "$SE_LRM" "$ANA")"
  printf "%-8s %-10s %-14.6f %-14.9f %-12s %-6s\n" "Vega" "LRM" "$EST_LRM" "$SE_LRM" "$Z_LRM" "$PASS_LRM"

  read -r EST_BRV SE_BRV <<<"$(run_greek vega "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T" --bump-abs-sigma 0.01)"
  read -r Z_BRV PASS_BRV <<<"$(z_and_pass "$EST_BRV" "$SE_BRV" "$ANA")"
  printf "%-8s %-10s %-14.6f %-14.9f %-12s %-6s\n" "Vega" "BRV" "$EST_BRV" "$SE_BRV" "$Z_BRV" "$PASS_BRV"

  RAT_STD=$(awk -v a="$SE_BRV" -v b="$SE_LRM" 'BEGIN{ if(b==0){print "inf"} else printf("%.3f", a/b) }')
  RAT_VAR=$(awk -v a="$SE_BRV" -v b="$SE_LRM" 'BEGIN{ if(b==0){print "inf"} else printf("%.3f", (a/b)*(a/b)) }')
  echo "Std ratio (BRV/LRM): $RAT_STD ; Var ratio: $RAT_VAR"
done
echo

# 2) Antithétiques ON/OFF (LRM)
echo "==== [2] Antithetic impact (LRM) ===="
for case in "${CASES[@]}"; do
  IFS='|' read -r NAME TYPE S0 K R Q SIGMA T <<<"$case"
  echo
  printf "Case: %s (type=%s)\n" "$NAME" "$TYPE"
  hr
  printf "%-12s %-14s %-14s\n" "Mode" "MC_estimate" "MC_SE"

  # anti OFF
  ANTISAVE="$ANTI"; ANTI=0; ANTI_FLAG=""
  read -r E0 S0se <<<"$(run_greek vega_lrm "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
  printf "%-12s %-14.6f %-14.9f\n" "LRM(noAnti)" "$E0" "$S0se"

  # anti ON
  ANTI=1; ANTI_FLAG="--antithetic"
  read -r E1 S1se <<<"$(run_greek vega_lrm "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
  printf "%-12s %-14.6f %-14.9f\n" "LRM(anti)" "$E1" "$S1se"

  RAT=$(awk -v a="$S0se" -v b="$S1se" 'BEGIN{ if(b==0){print "inf"} else printf("%.3f", a/b) }')
  echo "Std ratio (noAnti/anti): $RAT"
  ANTI="$ANTISAVE"; ANTI_FLAG=$([[ $ANTI -eq 1 ]] && echo "--antithetic" || echo "")
done
echo

# 3) Invariance aux pas (LRM)
echo "==== [3] Steps invariance (LRM) ===="
CASE0="${CASES[0]}"; IFS='|' read -r NAME TYPE S0 K R Q SIGMA T <<<"$CASE0"
printf "Case: %s (ATM_1Y)\n" "$NAME"; hr
printf "%-7s %-14s %-14s %-14s\n" "steps" "MC_estimate" "MC_SE" "halfwidth95"

STEPS_SAVE="$STEPS"

STEPS=1
read -r E1 S1 <<<"$(run_greek vega_lrm "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
H1=$(halfwidth95 "$S1")
printf "%-7d %-14.6f %-14.9f %-14.9f\n" 1 "$E1" "$S1" "$H1"

STEPS=10
read -r E10 S10 <<<"$(run_greek vega_lrm "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
H10=$(halfwidth95 "$S10")
printf "%-7d %-14.6f %-14.9f %-14.9f\n" 10 "$E10" "$S10" "$H10"

STEPS=50
read -r E50 S50 <<<"$(run_greek vega_lrm "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
H50=$(halfwidth95 "$S50")
printf "%-7d %-14.6f %-14.9f %-14.9f\n" 50 "$E50" "$S50" "$H50"

STEPS="$STEPS_SAVE"

echo
echo "Done."
