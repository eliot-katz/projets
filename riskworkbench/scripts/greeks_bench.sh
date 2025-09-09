#!/usr/bin/env bash
set -euo pipefail

# ===========================================
# Greeks MC vs BS (z-scores) + extra checks
# ===========================================

BIN="build/bin/greeks_runner"
N=400000
BATCH=100000
STEPS=1
SEED=42
TOL=-1
PLOT=1
Z95=1.959963984540054

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin)     BIN="$2"; shift 2;;
    --n)       N="$2"; shift 2;;
    --batch)   BATCH="$2"; shift 2;;
    --steps)   STEPS="$2"; shift 2;;
    --seed)    SEED="$2"; shift 2;;
    --tol)     TOL="$2"; shift 2;;
    --no-plot) PLOT=0; shift;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: MC binary not found/executable at: $BIN" >&2
  exit 1
fi


# --------- Cas à tester (nom|type|S0|K|r|q|sigma|T) ----------
CASES=(
  "ATM_1Y|call|100|100|0.02|0.00|0.20|1.00"
  "OTM_3M|call|100|110|0.02|0.00|0.20|0.25"
  "ITM_2Y_PUT|put|100|120|0.02|0.00|0.20|2.00"
)

# Bumps par défaut
BUMP_REL_S0=0.01
BUMP_ABS_SIGMA=0.01
BUMP_ABS_R=0.0001
BUMP_ABS_T=0.0027397260273972603

# Seuils de robustesse pour SE quasi-nul et diff MC vs BS
ZERO_SE_ABS=1e-12
DIFF_ABS_PASS=1e-8
DIFF_REL_PASS=1e-9

# CRN & antithétiques (ON ici pour le bench principal)
CRN_FLAG="--crn"
ANTI_FLAG="--antithetic"

CSV_OUT="greeks_bench.csv"
: > "$CSV_OUT"
echo "case,greek,estimate,std_error,analytic,zscore,pass" >> "$CSV_OUT"

# --------- Utilitaires ---------
parse_field() { echo "$1" | sed -E 's/^[^:]+:[[:space:]]*//'; }
hr() { printf '%s\n' "----------------------------------------------------------------------------------------"; }

run_mc() {
  # Args: greek type S0 K r q sigma T
  local greek="$1" type="$2" S0="$3" K="$4" r="$5" q="$6" sigma="$7" T="$8"
  local PUTFLAG=""; [[ "$type" == "put" ]] && PUTFLAG="--put"
  "$BIN" "$S0" "$K" "$r" "$q" "$sigma" "$T" "$SEED" "$N" "$BATCH" "$TOL" "$STEPS" \
        --greek "$greek" $PUTFLAG $CRN_FLAG $ANTI_FLAG \
        --bump-rel-s0 "$BUMP_REL_S0" \
        --bump-abs-sigma "$BUMP_ABS_SIGMA" \
        --bump-abs-r "$BUMP_ABS_R" \
        --bump-abs-T "$BUMP_ABS_T"
}

run_mc_custom() {
  # Args: greek type S0 K r q sigma T  crn(0/1) anti(0/1) bump_rel_s0 bump_abs_sigma bump_abs_r bump_abs_T
  local greek="$1" type="$2" S0="$3" K="$4" r="$5" q="$6" sigma="$7" T="$8"
  local crn="$9" anti="${10}" bS0="${11}" bSig="${12}" bR="${13}" bT="${14}"
  local PUTFLAG=""; [[ "$type" == "put" ]] && PUTFLAG="--put"
  local CRN_SW="--no-crn"; [[ "$crn" -eq 1 ]] && CRN_SW="--crn"
  local ANTI_SW="";        [[ "$anti" -eq 1 ]] && ANTI_SW="--antithetic"
  "$BIN" "$S0" "$K" "$r" "$q" "$sigma" "$T" "$SEED" "$N" "$BATCH" "$TOL" "$STEPS" \
        --greek "$greek" $PUTFLAG $CRN_SW $ANTI_SW \
        --bump-rel-s0 "$bS0" \
        --bump-abs-sigma "$bSig" \
        --bump-abs-r "$bR" \
        --bump-abs-T "$bT"
}

get_est_se() {
  # lit estimate + std_error dans la sortie de run_mc(_custom)
  local out; out="$1"
  local est se
  est=$(parse_field "$(echo "$out" | grep -E '^estimate')" || true)
  se=$(parse_field  "$(echo "$out" | grep -E '^std_error')" || true)
  echo "$est $se"
}

bs_greek() {
  # Args: greek type S0 K r q sigma T
  python3 - "$@" <<'PY'
import sys, math
greek, opttype, S0, K, r, q, sigma, T = sys.argv[1], sys.argv[2], *map(float, sys.argv[3:])
def Phi(x): return 0.5*math.erfc(-x/math.sqrt(2.0))
def pdf(x): return math.exp(-0.5*x*x)/math.sqrt(2.0*math.pi)
if T<=0 or sigma<=0:
    d1 = float('inf') if S0>K else float('-inf')
else:
    d1 = (math.log(S0/K) + (r - q + 0.5*sigma*sigma)*T)/(sigma*math.sqrt(T))
d2 = d1 - sigma*math.sqrt(T) if T>0 else d1
disc_q = math.exp(-q*T); disc_r = math.exp(-r*T)
if opttype == "call":
    Delta = disc_q*Phi(d1)
    Rho   =  T * K * disc_r * Phi(d2)
    Theta = (-(S0*disc_q*pdf(d1)*sigma)/(2*math.sqrt(T)) - r*K*disc_r*Phi(d2) + q*S0*disc_q*Phi(d1)) if T>0 else 0.0
else:
    Delta = -disc_q*Phi(-d1)
    Rho   = -T * K * disc_r * Phi(-d2)
    Theta = (-(S0*disc_q*pdf(d1)*sigma)/(2*math.sqrt(T)) + r*K*disc_r*Phi(-d2) - q*S0*disc_q*Phi(-d1)) if T>0 else 0.0
Vega = S0*disc_q*math.sqrt(T)*pdf(d1) if T>0 else 0.0
val = {"delta":Delta, "vega":Vega, "rho":Rho, "theta":Theta}[greek]
print(f"{val:.12f}")
PY
}

# ==================== Bench principal: MC vs BS (z-scores) ====================
printf "\n== Greeks MC vs BS (N=%d, batch=%d, steps=%d, seed=%d) ==\n\n" "$N" "$BATCH" "$STEPS" "$SEED"

for case in "${CASES[@]}"; do
  IFS='|' read -r NAME TYPE S0 K R Q SIGMA T <<<"$case"

  printf "Case: %s (type=%s, S0=%.4f K=%.4f r=%.4f q=%.4f sigma=%.4f T=%.4f)\n" \
    "$NAME" "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T"
  printf "%-8s  %-12s  %-12s  %-12s  %-10s  %-4s\n" "Greek" "MC_estimate" "MC_SE" "BS_analytic" "z_score" "PASS"
  hr

  for G in delta vega rho theta; do
    OUT="$(run_mc "$G" "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"
    read -r EST SE <<<"$(get_est_se "$OUT")"
    ANA="$(bs_greek "$G" "$TYPE" "$S0" "$K" "$R" "$Q" "$SIGMA" "$T")"

    # z = (EST - ANA)/SE avec gestion SE ~ 0
    # Règle: si SE < ZERO_SE_ABS, on "passe" si |EST-ANA| <= max(DIFF_ABS_PASS, DIFF_REL_PASS*|ANA|)
    if awk -v se="$SE" -v zthr="$ZERO_SE_ABS" 'BEGIN{exit !(se<zthr)}'; then
      DIFF=$(awk -v e="$EST" -v a="$ANA" 'BEGIN{d=e-a; if(d<0)d=-d; printf("%.12g", d)}')
      THR=$(awk -v a="$ANA" -v thr_abs="$DIFF_ABS_PASS" -v thr_rel="$DIFF_REL_PASS" \
               'BEGIN{absa=a; if(absa<0)absa=-absa; t=thr_abs; rel=thr_rel*absa; if(rel>t)t=rel; printf("%.12g", t)}')
      if awk -v d="$DIFF" -v t="$THR" 'BEGIN{exit !(d<=t)}'; then
        Z="0.000"; PASS="YES"
      else
        Z="inf";   PASS="NO"
      fi
    else
      Z=$(awk -v e="$EST" -v a="$ANA" -v s="$SE" 'BEGIN{printf("%.3f",(e-a)/s)}')
      ABSZ=$(awk -v z="$Z" 'BEGIN{z=(z<0?-z:z); printf("%.3f",z)}')
      PASS=$(awk -v a="$ABSZ" 'BEGIN{print (a<=2.0)?"YES":"NO"}')
    fi


    printf "%-8s  %-12.6f  %-12.9f  %-12.6f  %-10s  %-4s\n" "$G" "$EST" "$SE" "$ANA" "$Z" "$PASS"
    echo "$NAME,$G,$EST,$SE,$ANA,$Z,$PASS" >> "$CSV_OUT"
  done
  echo
done

echo "CSV -> $CSV_OUT"

# ------- Plot z-scores (optionnel via gnuplot) -------
if [[ "$PLOT" -eq 1 ]] && command -v gnuplot >/dev/null 2>&1; then
  GNUPLOT_SCRIPT="$(mktemp)"
  cat > "$GNUPLOT_SCRIPT" <<'GP'
set datafile separator comma
set terminal pngcairo size 1000,600
set output 'greeks_zscores.png'
set title "Greeks MC vs BS — z-scores ((estimate - analytic)/SE)"
set key outside
set grid ytics
set yrange [-5:5]
set style data points
set pointsize 1.2
set xlabel "Row index"
set ylabel "z-score"
plot \
  'greeks_bench.csv' using (stringcolumn(2) eq "delta" ? $0 : 1/0):6 with points pt 7  title "delta", \
  'greeks_bench.csv' using (stringcolumn(2) eq "vega"  ? $0 : 1/0):6 with points pt 5  title "vega",  \
  'greeks_bench.csv' using (stringcolumn(2) eq "rho"   ? $0 : 1/0):6 with points pt 9  title "rho",   \
  'greeks_bench.csv' using (stringcolumn(2) eq "theta" ? $0 : 1/0):6 with points pt 13 title "theta"
unset output
GP
  gnuplot "$GNUPLOT_SCRIPT" && rm -f "$GNUPLOT_SCRIPT"
  echo "Plot -> greeks_zscores.png"
else
  echo "(Plot skipped: gnuplot not found or --no-plot)"
fi

# ====================== EXTRA CHECKS (exigences 7) ===========================

echo
echo "== Extra checks =="

# ---------- 1) Effet CRN (Delta, ATM 1Y call) ----------
echo
echo "-- CRN effect on SE (Delta, ATM_1Y call) --"
OUT_NOCRN="$(run_mc_custom delta call 100 100 0.02 0.00 0.20 1.00  0 1  0.01 0.01 0.0001 0.0027397260273972603)"
OUT_CRN="$(  run_mc_custom delta call 100 100 0.02 0.00 0.20 1.00  1 1  0.01 0.01 0.0001 0.0027397260273972603)"
read -r _ SE_no <<<"$(get_est_se "$OUT_NOCRN")"
read -r _ SE_crn <<<"$(get_est_se "$OUT_CRN")"
printf "SE(no-CRN) = %.6f ; SE(CRN) = %.6f ; gain = %.3fx on stdev\n" "$SE_no" "$SE_crn" "$(awk -v a="$SE_no" -v b="$SE_crn" 'BEGIN{print a/b}')"
PASS_CRN=$(awk -v a="$SE_no" -v b="$SE_crn" 'BEGIN{print (b<a)?"YES":"NO"}')
echo "PASS (SE drops with CRN) = $PASS_CRN"

# ---------- 2) Stabilité bump Delta (eps = 1%, 0.5%, 0.1%) ----------
echo
echo "-- Bump stability (Delta, eps = 1%, 0.5%, 0.1%) --"
declare -A EST_D SE_D
for EPS in 0.01 0.005 0.001; do
  OUT="$(run_mc_custom delta call 100 100 0.02 0.00 0.20 1.00  1 1  $EPS 0.01 0.0001 0.0027397260273972603)"
  read -r EST SE <<<"$(get_est_se "$OUT")"
  EST_D["$EPS"]="$EST"; SE_D["$EPS"]="$SE"
  printf "eps=%-6s  estimate=%.6f  SE=%.6f  CI=[%.6f, %.6f]\n" \
    "$EPS" "$EST" "$SE" "$(awk -v m="$EST" -v s="$SE" -v z="$Z95" 'BEGIN{printf("%.6f", m-z*s)}')" \
                         "$(awk -v m="$EST" -v s="$SE" -v z="$Z95" 'BEGIN{printf("%.6f", m+z*s)}')"
done
# Pairwise z on diff of means (conservatif)
z01=$(awk -v m1="${EST_D[0.01]}"  -v s1="${SE_D[0.01]}"  -v m2="${EST_D[0.005]}" -v s2="${SE_D[0.005]}" 'BEGIN{printf("%.3f", ((m1-m2)>=0?(m1-m2):-(m1-m2))/sqrt(s1*s1+s2*s2))}')
z02=$(awk -v m1="${EST_D[0.01]}"  -v s1="${SE_D[0.01]}"  -v m2="${EST_D[0.001]}" -v s2="${SE_D[0.001]}" 'BEGIN{printf("%.3f", ((m1-m2)>=0?(m1-m2):-(m1-m2))/sqrt(s1*s1+s2*s2))}')
z12=$(awk -v m1="${EST_D[0.005]}" -v s1="${SE_D[0.005]}" -v m2="${EST_D[0.001]}" -v s2="${SE_D[0.001]}" 'BEGIN{printf("%.3f", ((m1-m2)>=0?(m1-m2):-(m1-m2))/sqrt(s1*s1+s2*s2))}')
echo "pairwise |z|: 1% vs 0.5% = $z01 ; 1% vs 0.1% = $z02 ; 0.5% vs 0.1% = $z12"
PASS_BUMP_DELTA=$(awk -v a="$z01" -v b="$z02" -v c="$z12" 'BEGIN{print (a<=2 && b<=2 && c<=2)?"YES":"NO"}')
echo "PASS (all |z|<=2) = $PASS_BUMP_DELTA"

# ---------- 3) Stabilité bump Vega (h = 0.01, 0.005) ----------
echo
echo "-- Bump stability (Vega, h = 0.01, 0.005) --"
declare -A EST_V SE_V
for H in 0.01 0.005; do
  OUT="$(run_mc_custom vega call 100 100 0.02 0.00 0.20 1.00  1 1  0.01 $H 0.0001 0.0027397260273972603)"
  read -r EST SE <<<"$(get_est_se "$OUT")"
  EST_V["$H"]="$EST"; SE_V["$H"]="$SE"
  printf "h=%-6s   estimate=%.6f  SE=%.6f  CI=[%.6f, %.6f]\n" \
    "$H" "$EST" "$SE" "$(awk -v m="$EST" -v s="$SE" -v z="$Z95" 'BEGIN{printf("%.6f", m-z*s)}')" \
                       "$(awk -v m="$EST" -v s="$SE" -v z="$Z95" 'BEGIN{printf("%.6f", m+z*s)}')"
done
zv=$(awk -v m1="${EST_V[0.01]}" -v s1="${SE_V[0.01]}" -v m2="${EST_V[0.005]}" -v s2="${SE_V[0.005]}" 'BEGIN{printf("%.3f", ((m1-m2)>=0?(m1-m2):-(m1-m2))/sqrt(s1*s1+s2*s2))}')
echo "pairwise |z|: 0.01 vs 0.005 = $zv"
PASS_BUMP_VEGA=$(awk -v a="$zv" 'BEGIN{print (a<=2)?"YES":"NO"}')
echo "PASS (|z|<=2) = $PASS_BUMP_VEGA"

# ---------- 4) Antithétiques vs BRV (Delta, ATM 1Y) ----------
echo
echo "-- Antithetic vs plain (Delta, ATM_1Y call, CRN=ON) --"
SE_plain=$(run_mc_custom delta call 100 100 0.02 0.00 0.20 1.00  1 0  0.01 0.01 0.0001 0.0027397260273972603 | awk '/^std_error/ {print $3}')
SE_anti=$( run_mc_custom delta call 100 100 0.02 0.00 0.20 1.00  1 1  0.01 0.01 0.0001 0.0027397260273972603 | awk '/^std_error/ {print $3}')
R=$(awk -v sp="$SE_plain" -v sa="$SE_anti" 'BEGIN{printf("%.3f",(sp/sa)^2)}')
printf "SE_plain=%.6f  SE_anti=%.6f  R=(sp/sa)^2 = %s\n" "$SE_plain" "$SE_anti" "$R"
PASS_ANTI=$(awk -v R="$R" 'BEGIN{print (R>=1.1)?"YES":"NO"}')
echo "PASS (R >= 1.1 expected) = $PASS_ANTI"

echo
echo "== Summary extra checks =="
echo "CRN drop         : $PASS_CRN"
echo "Delta bump stable: $PASS_BUMP_DELTA"
echo "Vega  bump stable: $PASS_BUMP_VEGA"
echo "Anti vs plain    : $PASS_ANTI"
