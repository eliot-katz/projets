#!/usr/bin/env bash
# Variance Reduction Integration Tests (Antithetic & Control Variates)
# Requires:
#   - build/bin/bs_ref
#   - build/bin/mc_runner (supports --antithetic, --cv, --cv-pilot, --cv-custom-K, --cv-type)
set -euo pipefail
export LC_ALL=C

# ---------- Config (overridable via env) -------------------------------------
BUILD_DIR="${BUILD_DIR:-build}"
BIN_DIR="${BUILD_DIR}/bin"
BS_REF="${BIN_DIR}/bs_ref"
MC_RUNNER="${BIN_DIR}/mc_runner"

SEED=${SEED:-42}
N=${N:-400000}          # physical paths (same budget for all modes)
BATCH=${BATCH:-100000}  # physical paths per batch
STEPS=${STEPS:-1}
PILOT=${PILOT:-20000}   # pilot physical paths for beta
TOL=${TOL:--1}          # disable tol-based early stop by default

# Test set: S0 K r q sigma T type name
CASES=(
  "100 100 0.02 0.0 0.20 1.00 call ATM_1Y"
  "100 110 0.02 0.0 0.20 0.25 call OTM_3M"
  "100 120 0.02 0.0 0.20 2.00  put  ITM_2Y_PUT"
)

# CV voisin à tester (en plus du “same instrument”)
# ex: 0.95*K et 1.05*K (forte corrélation mais SE > 0)
KCV_FACTORS=(0.95 1.05)

# ---------- Build if needed ---------------------------------------------------
need_build=false
[[ ! -x "$BS_REF" ]] && need_build=true
[[ ! -x "$MC_RUNNER" ]] && need_build=true

if $need_build; then
  echo "[build] Binaries not found. Building (Release)..."
  cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
  cmake --build "${BUILD_DIR}" -j
  echo "[build] Done."
fi

# ---------- Helpers -----------------------------------------------------------
extract_val () { printf "%s\n" "$1" | grep -E "^$2[[:space:]]*:" | sed -E "s/^$2[[:space:]]*:[[:space:]]*//"; }
in_ci () { awk -v x="$1" -v lo="$2" -v hi="$3" 'BEGIN{print (x>=lo && x<=hi)?"YES":"NO"}'; }

bs_price () {
  local out price
  out="$("$BS_REF" "$1" "$2" "$3" "$4" "$5" "$6")"
  if [[ "$7" == "call" ]]; then
    price=$(printf "%s\n" "$out" | awk 'NF>0 && $1 ~ /^[[:space:]]*[0-9]/ {print $(NF-2)}' | tail -n1)
  else
    price=$(printf "%s\n" "$out" | awk 'NF>0 && $1 ~ /^[[:space:]]*[0-9]/ {print $(NF-1)}' | tail -n1)
  fi
  printf "%s" "$price"
}

# run_mc mode S0 K r q sigma T type [steps_override]
run_mc () {
  local mode="$1"; shift
  local S0="$1" K="$2" r="$3" q="$4" sig="$5" T="$6" typ="$7"
  local steps_override="${8:-}"
  local steps="$STEPS"; [[ -n "$steps_override" ]] && steps="$steps_override"

  local args=( "$S0" "$K" "$r" "$q" "$sig" "$T" "$SEED" "$N" "$BATCH" "$TOL" "$steps" )
  [[ "$typ" == "put" ]] && args+=( "--put" )
  case "$mode" in
    plain)    ;;
    anti)     args+=( "--antithetic" ) ;;
    cv)       args+=( "--cv" "--cv-pilot" "$PILOT" ) ;;
    anti_cv)  args+=( "--antithetic" "--cv" "--cv-pilot" "$PILOT" ) ;;
    *) echo "Unknown mode: $mode" >&2; exit 2;;
  esac
  "$MC_RUNNER" "${args[@]}"
}

# run_mc_cv_custom mode_with_cv S0 K ... type Kcv [steps_override]
run_mc_cv_custom () {
  local mode="$1"; shift
  local S0="$1" K="$2" r="$3" q="$4" sig="$5" T="$6" typ="$7" Kcv="$8"
  local steps_override="${9:-}"
  local steps="$STEPS"; [[ -n "$steps_override" ]] && steps="$steps_override"

  local args=( "$S0" "$K" "$r" "$q" "$sig" "$T" "$SEED" "$N" "$BATCH" "$TOL" "$steps" )
  [[ "$typ" == "put" ]] && args+=( "--put" )
  case "$mode" in
    cv)      args+=( "--cv" "--cv-custom-K" "$Kcv" "--cv-type" "$typ" "--cv-pilot" "$PILOT" ) ;;
    anti_cv) args+=( "--antithetic" "--cv" "--cv-custom-K" "$Kcv" "--cv-type" "$typ" "--cv-pilot" "$PILOT" ) ;;
    *) echo "Unknown custom-cv mode: $mode" >&2; exit 2;;
  esac
  "$MC_RUNNER" "${args[@]}"
}

# R = (SE_plain/SE_mode)^2 ; returns "inf" si SE_mode==0
compute_R () {
  local se_plain="$1" se_mode="$2"
  if awk -v b="$se_mode" 'BEGIN{exit (b==0)?0:1}'; then
    echo "inf"
  else
    awk -v a="$se_plain" -v b="$se_mode" 'BEGIN{print (a>0&&b>0)?(a/b)^2:0}'
  fi
}

# returns success (0) si R>=1, ou si R=="inf"
ratio_ok () {
  local R="$1"
  if [[ "$R" == "inf" ]]; then return 0; fi
  awk -v R="$R" 'BEGIN{exit (R>=1)?0:1}'
}

fmt_row () { printf "%-12s %-12s %-12s %-12s %-12s %-6s %-8s\n" "$@"; }

# ---------- Tests -------------------------------------------------------------
total_fail=0

echo
echo "==== Variance Reduction Integration Tests (N=${N}, batch=${BATCH}, steps=${STEPS}, seed=${SEED}) ===="
echo

for case_line in "${CASES[@]}"; do
  read -r S0 K r q sig T typ name <<<"$case_line"

  BS=$(bs_price "$S0" "$K" "$r" "$q" "$sig" "$T" "$typ")

  # Run 4 modes (CV same-instrument)
  out_plain="$(run_mc plain    "$S0" "$K" "$r" "$q" "$sig" "$T" "$typ")"
  out_anti="$( run_mc anti     "$S0" "$K" "$r" "$q" "$sig" "$T" "$typ")"
  out_cv="$(   run_mc cv       "$S0" "$K" "$r" "$q" "$sig" "$T" "$typ")"
  out_acv="$(  run_mc anti_cv  "$S0" "$K" "$r" "$q" "$sig" "$T" "$typ")"

  # Extract
  Pp=$(extract_val "$out_plain" price);  SEp=$(extract_val "$out_plain" std_error);  LOp=$(extract_val "$out_plain" ci_low);  HIp=$(extract_val "$out_plain" ci_high)
  Pa=$(extract_val "$out_anti"  price);  SEa=$(extract_val "$out_anti"  std_error);  LOa=$(extract_val "$out_anti"  ci_low);  HIa=$(extract_val "$out_anti"  ci_high)
  Pc=$(extract_val "$out_cv"    price);  SEc=$(extract_val "$out_cv"    std_error);  LOc=$(extract_val "$out_cv"    ci_low);  HIc=$(extract_val "$out_cv"    ci_high)
  Px=$(extract_val "$out_acv"   price);  SEx=$(extract_val "$out_acv"   std_error);  LOx=$(extract_val "$out_acv"   ci_low);  HIx=$(extract_val "$out_acv"   ci_high)

  # Ratios vs plain (gèrent SE_mode==0 -> inf)
  R_anti=$(compute_R "$SEp" "$SEa")
  R_cv=$(  compute_R "$SEp" "$SEc")
  R_acv=$( compute_R "$SEp" "$SEx")

  # CI checks
  CIp=$(in_ci "$BS" "$LOp" "$HIp")
  CIa=$(in_ci "$BS" "$LOa" "$HIa")
  CIc=$(in_ci "$BS" "$LOc" "$HIc")
  CIx=$(in_ci "$BS" "$LOx" "$HIx")

  # 1) All prices’ CI cover BS
  if [[ "$CIp" == "YES" && "$CIa" == "YES" && "$CIc" == "YES" && "$CIx" == "YES" ]]; then
    pass_price=1
  else
    pass_price=0
  fi

  # 2) Variance ratios >= 1 (compte "inf" comme OK)
  pass_var=1
  ratio_ok "$R_anti" || pass_var=0
  ratio_ok "$R_cv"   || pass_var=0
  ratio_ok "$R_acv"  || pass_var=0

  # 3) SE(anti+cv) <= SE(cv) — marche aussi si SE=0 dans les deux cas
  eps="${ACV_CV_REL_EPS:-0.02}"
  if awk -v c="$SEc" -v x="$SEx" -v eps="$eps" 'BEGIN{exit (x<=c*(1+eps))?0:1}'; then
    pass_acv_le_cv=1
  else
    pass_acv_le_cv=0
  fi

  # Print table
  echo "Case: $name (type=$typ, S0=$S0 K=$K r=$r q=$q sigma=$sig T=$T) | BS_ref=$BS"
  printf "%-12s %-12s %-12s %-12s %-12s %-6s %-8s\n" "mode" "price" "std_error" "ci_low" "ci_high" "in_CI" "R(plain)"
  printf "%s\n" "-----------------------------------------------------------------------------------------------"
  fmt_row plain    "$Pp" "$SEp" "$LOp" "$HIp" "$CIp" "1.00"
  fmt_row anti     "$Pa" "$SEa" "$LOa" "$HIa" "$CIa" "$(printf '%s' "$R_anti")"
  fmt_row "cv(same)"      "$Pc" "$SEc" "$LOc" "$HIc" "$CIc" "$(printf '%s' "$R_cv")"
  fmt_row "anti+cv(same)" "$Px" "$SEx" "$LOx" "$HIx" "$CIx" "$(printf '%s' "$R_acv")"

  echo "Checks(same): prices_in_CI=${pass_price}  ratios>=1=${pass_var}  anti+cv<=cv=${pass_acv_le_cv}"
  if (( pass_price != 1 || pass_var != 1 || pass_acv_le_cv != 1 )); then
    total_fail=$((total_fail+1))
  fi

  # -------- CV voisin (Custom Kcv) -------------------------------------------
  for f in "${KCV_FACTORS[@]}"; do
    Kcv=$(awk -v k="$K" -v f="$f" 'BEGIN{print k*f}')
    out_cv_c="$(  run_mc_cv_custom cv      "$S0" "$K" "$r" "$q" "$sig" "$T" "$typ" "$Kcv")"
    out_acv_c="$( run_mc_cv_custom anti_cv "$S0" "$K" "$r" "$q" "$sig" "$T" "$typ" "$Kcv")"

    Pc2=$(extract_val "$out_cv_c"   price);  SEc2=$(extract_val "$out_cv_c"   std_error);  LOc2=$(extract_val "$out_cv_c"   ci_low);  HIc2=$(extract_val "$out_cv_c"   ci_high)
    Px2=$(extract_val "$out_acv_c"  price);  SEx2=$(extract_val "$out_acv_c"  std_error);  LOx2=$(extract_val "$out_acv_c"  ci_low);  HIx2=$(extract_val "$out_acv_c"  ci_high)

    R_cv2=$(  compute_R "$SEp" "$SEc2")
    R_acv2=$( compute_R "$SEp" "$SEx2")

    CIc2=$(in_ci "$BS" "$LOc2" "$HIc2")
    CIx2=$(in_ci "$BS" "$LOx2" "$HIx2")

    pass_price2=1; [[ "$CIp" != "YES" || "$CIa" != "YES" || "$CIc2" != "YES" || "$CIx2" != "YES" ]] && pass_price2=0
    pass_var2=1; ratio_ok "$R_cv2" || pass_var2=0; ratio_ok "$R_acv2" || pass_var2=0
    pass_acv_le_cv2=1
    eps="${ACV_CV_REL_EPS:-0.02}"
    if awk -v c="$SEc2" -v x="$SEx2" -v eps="$eps" 'BEGIN{exit (x<=c*(1+eps))?0:1}'; then :; else pass_acv_le_cv2=0; fi

    echo
    echo "CV custom (Kcv=$(printf '%.4f' "$Kcv")):"
    printf "%-12s %-12s %-12s %-12s %-12s %-6s %-8s\n" "mode" "price" "std_error" "ci_low" "ci_high" "in_CI" "R(plain)"
    printf "%s\n" "-----------------------------------------------------------------------------------------------"
    fmt_row "cv(custom)"     "$Pc2" "$SEc2" "$LOc2" "$HIc2" "$CIc2" "$(printf '%s' "$R_cv2")"
    fmt_row "anti+cv(custom)" "$Px2" "$SEx2" "$LOx2" "$HIx2" "$CIx2" "$(printf '%s' "$R_acv2")"
    echo "Checks(custom): prices_in_CI=${pass_price2}  ratios>=1=${pass_var2}  anti+cv<=cv=${pass_acv_le_cv2}"
    if (( pass_price2 != 1 || pass_var2 != 1 || pass_acv_le_cv2 != 1 )); then
      total_fail=$((total_fail+1))
    fi
  done

  # -------- Extra: steps invariance (anti+cv same, steps=1 vs 10) ------------
  out_x_s1="$out_acv"
  out_x_s10="$(run_mc anti_cv "$S0" "$K" "$r" "$q" "$sig" "$T" "$typ" 10)"
  Px1=$(extract_val "$out_x_s1" price); LOx1=$(extract_val "$out_x_s1" ci_low); HIx1=$(extract_val "$out_x_s1" ci_high)
  Px10=$(extract_val "$out_x_s10" price); LOx10=$(extract_val "$out_x_s10" ci_low); HIx10=$(extract_val "$out_x_s10" ci_high)
  overlap=$(awk -v a1="$LOx1" -v b1="$HIx1" -v a2="$LOx10" -v b2="$HIx10" 'BEGIN{lo=(a1>a2?a1:a2); hi=(b1<b2?b1:b2); print (hi>=lo)?"1":"0";}')
  echo
  echo "Extra: anti+cv(same) steps invariance -> steps=1 price=${Px1} ; steps=10 price=${Px10} ; IC overlap=${overlap}"
  if [[ "$overlap" -ne 1 ]]; then total_fail=$((total_fail+1)); fi

  # -------- Extra: reproducibility (same seed/config) ------------------------
  out_rep="$(run_mc anti_cv "$S0" "$K" "$r" "$q" "$sig" "$T" "$typ")"
  Px_rep=$(extract_val "$out_rep" price)
  if awk -v a="$Px" -v b="$Px_rep" 'BEGIN{da=a-b; if (da<0) da=-da; exit (da<=1e-12)?0:1}'; then same=1; else same=0; fi
  echo "Extra: reproducibility (same seed & config): equal_price=${same}"
  if [[ "$same" -ne 1 ]]; then total_fail=$((total_fail+1)); fi

  echo
done

# ---------- Summary -----------------------------------------------------------
if [[ "$total_fail" -eq 0 ]]; then
  echo "✓ All variance reduction tests PASSED."
  exit 0
else
  echo "✗ Some variance reduction tests FAILED (count=${total_fail})."
  exit 1
fi
