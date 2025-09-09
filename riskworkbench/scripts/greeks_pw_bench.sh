#!/usr/bin/env bash
set -euo pipefail

# ================= Greeks PW Bench =================
# - Compare Delta PW / Delta BRV / BS Delta (z-scores, SE, ratio)
# - Compare Gamma via Delta-PW (CRN) vs BS Gamma, eps = 1% et 0.5%
# - Invariance aux pas (Delta PW): steps = 1,10,50
#
# Prérequis: build/bin/greeks_runner (avec --greek delta|delta_pw|gamma_pw|vega|rho|theta)
# Par défaut : N=400000, batch=100000, seed=42, tol=-1, steps=1,
#              antithetic=ON, CRN=ON (utile surtout pour gamma_pw)
# ====================================================

BIN="build/bin/greeks_runner"
N=400000
BATCH=100000
SEED=42
TOL=-1
STEPS=1
ANTI=1
CRN=1

# Seuils robustesse (SE quasi-nul)
ZERO_SE_ABS=1e-12
DIFF_ABS_PASS=1e-8
DIFF_REL_PASS=1e-9

# Cas à tester: NAME|type|S0|K|r|q|sigma|T
CASES=(
  "ATM_1Y|call|100|100|0.02|0.00|0.20|1.00"
  "OTM_3M|call|100|110|0.02|0.00|0.20|0.25"
  "ITM_2Y_PUT|put|100|120|0.02|0.00|0.20|2.00"
)

# Bumps pour gamma via ΔPW
EPS_LIST=("0.010" "0.005")  # 1% et 0.5%

# ---------- options ----------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin) BIN="$2"; shift 2;;
    --n) N="$2"; shift 2;;
    --batch) BATCH="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --tol) TOL="$2"; shift 2;;
    --steps) STEPS="$2"; shift 2;;
    --anti) ANTI=1; shift;;
    --no-anti) ANTI=0; shift;;
    --crn) CRN=1; shift;;
    --no-crn) CRN=0; shift;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: binary not found/executable at: $BIN" >&2
  exit 1
fi

ANTI_FLAG=$([[ $ANTI -eq 1 ]] && echo "--antithetic" || echo "")
CRN_FLAG=$([[ $CRN -eq 1 ]] && echo "--crn" || echo "--no-crn")

# ---------- helpers ----------
parse_field() { echo "$1" | sed -E 's/^[^:]+:[[:space:]]*//'; } # "key : value" -> value

run_greek() {
  # args: greek type S0 K r q sigma T extra_flags...
  local greek="$1"; local type="$2"; local S0="$3"; local K="$4"; local r="$5"; local q="$6"; local sigma="$7"; local T="$8"
  shift 8
  local PUT=""
  [[ "$type" == "put" ]] && PUT="--put"
  local out
  out=$("$BIN" "$S0" "$K" "$r" "$q" "$sigma" "$T" "$SEED" "$N" "$BATCH" "$TOL" "$STEPS" \
        --greek "$greek" $PUT $CRN_FLAG $ANTI_FLAG "$@")
  local est se ci_low ci_high
  est=$(parse_field "$(echo "$out" | grep -E '^estimate')" )
  se=$( parse_field "$(echo "$out" | grep -E '^std_error')" )
  ci_low=$( parse_field "$(echo "$out" | grep -E '^ci_low')" )
  ci_high=$(parse_field "$(echo "$out" | grep -E '^ci_high')")
  echo "$est $se $ci_low $ci_high"
}

bs_greek() {
  # args: greek type S0 K r q sigma T
  python3 - "$@" <<'PY'
import sys, math
greek, typ, S0, K, r, q, sigma, T = sys.argv[1], sys.argv[2], *map(float, sys.argv[3:])
SQRT2PI = (2.0*math.pi)**0.5
def Phi(x): return 0.5*math.erfc(-x/(2.0**0.5))
def pdf(x): return math.exp(-0.5*x*x)/SQRT2PI

if T<=0 or sigma<=0:
    d1 = float('inf') if S0> K else float('-inf')
else:
    d1 = (math.log(S0/K) + (r - q + 0.5*sigma*sigma)*T)/(sigma*math.sqrt(T))
d2 = d1 - sigma*math.sqrt(T) if T>0 else d1

disc_q = math.exp(-q*T)
disc_r = math.exp(-r*T)

if typ=="call":
    Delta = disc_q*Phi(d1) if T>0 and sigma>0 else (1.0 if S0>K else 0.0)
else:
    Delta = -disc_q*Phi(-d1) if T>0 and sigma>0 else (-1.0 if S0<K else 0.0)

Vega  = S0*disc_q*math.sqrt(T)*pdf(d1) if T>0 and sigma>0 else 0.0
Rho   = ( T*K*disc_r*Phi(d2) if typ=="call" else -T*K*disc_r*Phi(-d2) ) if T>0 and sigma>0 else 0.0
Theta = (-(S0*disc_q*pdf(d1)*sigma)/(2*math.sqrt(T)) - r*K*disc_r*Phi(d2) + q*S0*disc_q*Phi(d1)) if (typ=="call" and T>0 and sigma>0) else \
        (-(S0*disc_q*pdf(d1)*sigma)/(2*math.sqrt(T)) + r*K*disc_r*Phi(-d2) - q*S0*disc_q*Phi(-d1)) if (typ!="call" and T>0 and sigma>0) else 0.0
Gamma = (disc_q*pdf(d1)/(S0*sigma*math.sqrt(T))) if T>0 and sigma>0 else 0.0

vals = {"delta":Delta, "gamma":Gamma, "vega":Vega, "rho":Rho, "theta":Theta}
print(f"{vals[greek]:.12f}")
PY
}

z_and_pass() {
  # args: EST SE ANA -> prints: "Z PASS"
  awk -v e="$1" -v s="$2" -v a="$3" -v zthr="$ZERO_SE_ABS" -v absthr="$DIFF_ABS_PASS" -v relthr="$DIFF_REL_PASS" '
    function abs(x){return x<0?-x:x}
    BEGIN{
      if (s<zthr){
        diff=abs(e-a); thr=absthr; absa=abs(a); if (relthr*absa>thr) thr=relthr*absa;
        if (diff<=thr) {print "0.000 YES"; exit}
        else {print "inf NO"; exit}
      } else {
        z=(e-a)/s; az=abs(z);
        pass=(az<=2.0)?"YES":"NO";
        printf("%.3f %s\n", z, pass);
      }
    }'
}

halfwidth95() {
  awk -v se="$1" 'BEGIN{z=1.959963984540054; printf("%.12f", z*se)}'
}

overlap_ci() {
  # args: mean1 low1 high1  mean2 low2 high2 -> YES/NO
  awk -v l1="$2" -v h1="$3" -v l2="$5" -v h2="$6" '
    BEGIN{ print ( (h1>=l2 && h2>=l1) ? "YES":"NO" ) }'
}

hr(){ printf '%s\n' "------------------------------------------------------------------------------------------"; }

printf "\n== Greeks PW Bench (N=%d, batch=%d, steps=%d, seed=%d, anti=%s, crn=%s) ==\n\n" \
  "$N" "$BATCH" "$STEPS" "$SEED" "$ANTI" "$CRN"

# =========================================================
# 6.1 Delta PW vs BRV vs BS
# =========================================================
echo "==== [6.1] Delta: PW vs BRV vs BS ===="
for case in "${CASES[@]}"; do
  IFS='|' read -r NAME TYPE S0 K R Q SIGMA T <<<"$case"
  echo
  printf "Case: %s (type=%s, S0=%g K=%g r=%g q=%g sigma=%g T=%g)\n" \
    "$NAME" "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T"
  hr
  printf "%-8s %-10s %-14s %-14s %-14s %-8s\n" "Greek" "Mode" "MC_estimate" "MC_SE" "z_vs_BS" "PASS"

  # BS Delta
  ANA_DELTA="$(bs_greek delta "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"

  # Delta PW
  read -r EST_PW SE_PW L1 L2 <<<"$(run_greek delta_pw "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
  read -r Z_PW PASS_PW <<<"$(z_and_pass "$EST_PW" "$SE_PW" "$ANA_DELTA")"

  # Delta BRV
  read -r EST_BRV SE_BRV L1b L2b <<<"$(run_greek delta "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
  read -r Z_BRV PASS_BRV <<<"$(z_and_pass "$EST_BRV" "$SE_BRV" "$ANA_DELTA")"

  printf "%-8s %-10s %-14.6f %-14.9f %-14s %-8s\n" "Delta" "PW"  "$EST_PW"  "$SE_PW"  "$Z_PW"  "$PASS_PW"
  printf "%-8s %-10s %-14.6f %-14.9f %-14s %-8s\n" "Delta" "BRV" "$EST_BRV" "$SE_BRV" "$Z_BRV" "$PASS_BRV"

  # Ratio de variance: R = (SE_BRV / SE_PW)^2  (affiche aussi ratio d'écart-types)
  RAT_STD=$(awk -v a="$SE_BRV" -v b="$SE_PW" 'BEGIN{ if(b==0){print "inf"} else printf("%.3f", a/b) }')
  RAT_VAR=$(awk -v a="$SE_BRV" -v b="$SE_PW" 'BEGIN{ if(b==0){print "inf"} else printf("%.3f", (a/b)*(a/b)) }')
  echo "Std ratio (BRV/PW): $RAT_STD ; Var ratio: $RAT_VAR"
done
echo

# =========================================================
# 6.2 Gamma via ΔPW (CRN) vs BS Gamma (eps=1%, 0.5%)
# =========================================================
echo "==== [6.2] Gamma via Delta-PW (CRN) vs BS Gamma ===="
for case in "${CASES[@]}"; do
  IFS='|' read -r NAME TYPE S0 K R Q SIGMA T <<<"$case"
  echo
  printf "Case: %s (type=%s, S0=%g K=%g r=%g q=%g sigma=%g T=%g)\n" \
    "$NAME" "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T"
  hr
  printf "%-8s %-8s %-14s %-14s %-14s %-8s\n" "Greek" "eps" "MC_estimate" "MC_SE" "z_vs_BS" "PASS"

  ANA_GAMMA="$(bs_greek gamma "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"

  for EPS in "${EPS_LIST[@]}"; do
    read -r EST_G SE_G Lh1 Lh2 <<<"$(run_greek gamma_pw "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T" --bump-rel-s0 "$EPS")"
    read -r ZG PASSG <<<"$(z_and_pass "$EST_G" "$SE_G" "$ANA_GAMMA")"
    printf "%-8s %-8s %-14.6f %-14.9f %-14s %-8s\n" "Gamma" "$EPS" "$EST_G" "$SE_G" "$ZG" "$PASSG"
  done
done
echo

# =========================================================
# 6.3 Invariance aux pas (Delta PW) : steps=1,10,50
# =========================================================
echo "==== [6.3] Invariance aux pas (Delta PW) ===="
# On fait l'exercice sur le cas ATM_1Y
CASE_INVAR="${CASES[0]}"
IFS='|' read -r NAME TYPE S0 K R Q SIGMA T <<<"$CASE_INVAR"
echo
printf "Case: %s (type=%s, S0=%g K=%g r=%g q=%g sigma=%g T=%g)\n" \
  "$NAME" "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T"
hr
printf "%-8s %-8s %-14s %-14s %-14s\n" "Greek" "steps" "MC_estimate" "MC_SE" "halfwidth95"

# steps=1
STEPS_SAVE="$STEPS"
STEPS=1
read -r E1 S1 C1L C1H <<<"$(run_greek delta_pw "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
H1=$(halfwidth95 "$S1")
printf "%-8s %-8d %-14.6f %-14.9f %-14.9f\n" "DeltaPW" 1 "$E1" "$S1" "$H1"

# steps=10
STEPS=10
read -r E10 S10 C10L C10H <<<"$(run_greek delta_pw "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
H10=$(halfwidth95 "$S10")
printf "%-8s %-8d %-14.6f %-14.9f %-14.9f\n" "DeltaPW" 10 "$E10" "$S10" "$H10"

# steps=50
STEPS=50
read -r E50 S50 C50L C50H <<<"$(run_greek delta_pw "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
H50=$(halfwidth95 "$S50")
printf "%-8s %-8d %-14.6f %-14.9f %-14.9f\n" "DeltaPW" 50 "$E50" "$S50" "$H50"

# check overlap des IC (1 vs 10) et (1 vs 50)
OV10=$(overlap_ci "$E1" "$(awk -v e="$E1" -v h="$H1" 'BEGIN{printf("%.12f", e-h)}')" "$(awk -v e="$E1" -v h="$H1" 'BEGIN{printf("%.12f", e+h)}')" \
                 "$E10" "$(awk -v e="$E10" -v h="$H10" 'BEGIN{printf("%.12f", e-h)}')" "$(awk -v e="$E10" -v h="$H10" 'BEGIN{printf("%.12f", e+h)}')" )
OV50=$(overlap_ci "$E1" "$(awk -v e="$E1" -v h="$H1" 'BEGIN{printf("%.12f", e-h)}')" "$(awk -v e="$E1" -v h="$H1" 'BEGIN{printf("%.12f", e+h)}')" \
                 "$E50" "$(awk -v e="$E50" -v h="$H50" 'BEGIN{printf("%.12f", e-h)}')" "$(awk -v e="$E50" -v h="$H50" 'BEGIN{printf("%.12f", e+h)}')" )

echo
echo "IC overlap (steps 1 vs 10): $OV10"
echo "IC overlap (steps 1 vs 50): $OV50"

# restore steps
STEPS="$STEPS_SAVE"

echo
echo "Done."
