#!/usr/bin/env bash
# Micro-expériences Phase 1 (RW)
# - Convergence en N (half_width_95 ~ 1/sqrt(N))
# - Invariance au pas (n_steps = 1,10,50)
# - Cas “or” : BS ∈ IC(MC) (ITM/OTM, T=0.25 et 2.0)
# - CSV de convergence (pour plot half_width vs 1/sqrt(N))
set -euo pipefail
export LC_ALL=C

# --- Config par défaut (modifiable) -------------------------------------------
BUILD_DIR="build"
BIN_DIR="${BUILD_DIR}/bin"
BS_REF="${BIN_DIR}/bs_ref"
MC_RUNNER="${BIN_DIR}/mc_runner"
GBM_MOMENTS="${BIN_DIR}/gbm_moments"   # optionnel (si tu l'as ajouté)

# Cas ATM 1Y
S0=100; K=100; r=0.02; q=0.0; sigma=0.2; T=1.0
SEED=42
BATCH=100000

# N pour la convergence
Ns=(50000 100000 200000 500000 1000000)

# Steps pour l'invariance
STEPS=(1 10 50)

# CSV convergence
CSV_OUT="conv.csv"

# --- Helpers -----------------------------------------------------------------
need_build=false
[[ ! -x "$BS_REF" ]] && need_build=true
[[ ! -x "$MC_RUNNER" ]] && need_build=true

if $need_build; then
  echo "[build] Binaries not found in ${BIN_DIR}. Building (Release)..."
  cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
  cmake --build "${BUILD_DIR}" -j
  echo "[build] Done."
fi

extract_val () {
  # $1: output, $2: key (e.g. price|std_error|ci_low|ci_high|elapsed_ms)
  printf "%s\n" "$1" | grep -E "^$2[[:space:]]*:" | sed -E "s/^$2[[:space:]]*:[[:space:]]*//"
}

run_mc () {
  # $1=N, $2=n_steps
  "$MC_RUNNER" "$S0" "$K" "$r" "$q" "$sigma" "$T" "$SEED" "$1" "$BATCH" "-1" "$2"
}

half_width () {
  # $1=lo, $2=hi
  awk -v lo="$1" -v hi="$2" 'BEGIN{print ((hi-lo)/2.0)}'
}

sqrtN () {
  awk -v n="$1" 'BEGIN{print (sqrt(n))}'
}

# --- 1) Convergence en N -----------------------------------------------------
echo
echo "=== [1] Convergence en N (ATM 1Y, seed=$SEED, steps=1) ==="
printf "%-10s %-14s %-14s %-14s %-12s\n" "N" "price" "std_error" "half_width" "c=half*sqrt(N)"
prev_c=""
for N in "${Ns[@]}"; do
  out="$(run_mc "$N" 1)"
  price="$(extract_val "$out" price)"
  se="$(extract_val "$out" std_error)"
  lo="$(extract_val "$out" ci_low)"
  hi="$(extract_val "$out" ci_high)"
  half="$(half_width "$lo" "$hi")"
  c=$(awk -v h="$half" -v n="$N" 'BEGIN{print h*sqrt(n)}')
  printf "%-10s %-14.8f %-14.8f %-14.8f %-12.8f\n" "$N" "$price" "$se" "$half" "$c"
  prev_c="$c"
done

# --- 2) Invariance au pas ----------------------------------------------------
echo
echo "=== [2] Invariance au pas (N=300000) ==="
Nsteps=300000
printf "%-8s %-14s %-14s %-14s\n" "steps" "price" "std_error" "half_width"
for s in "${STEPS[@]}"; do
  out="$(run_mc "$Nsteps" "$s")"
  price="$(extract_val "$out" price)"
  se="$(extract_val "$out" std_error)"
  lo="$(extract_val "$out" ci_low)"
  hi="$(extract_val "$out" ci_high)"
  half="$(half_width "$lo" "$hi")"
  printf "%-8s %-14.8f %-14.8f %-14.8f\n" "$s" "$price" "$se" "$half"
done

# --- 3) Cas “or” : BS ∈ IC(MC) ----------------------------------------------
echo
echo "=== [3] Cas “or” (BS dans l'IC MC) ==="
# Définition des cas : S0 K r q sigma T
cases=(
  "100  90  0.02 0 0.2 0.25"
  "100 110  0.02 0 0.2 0.25"
  "100  90  0.02 0 0.2 2.0"
  "100 110  0.02 0 0.2 2.0"
)
printf "%-5s %-7s %-7s %-6s %-6s %-6s %-6s  %-11s %-11s %-9s\n" "Type" "S0" "K" "r" "q" "sigma" "T" "BS_price" "MC_price" "PASS?"

for c in "${cases[@]}"; do
  read -r s0 k rr qq ss tt <<<"$c"
  # BS ref appelle un seul cas à la fois
  bs_out="$("$BS_REF" "$s0" "$k" "$rr" "$qq" "$ss" "$tt")"
  # bs_ref imprime une ligne de tableau ; on récupère call/put par position
  bs_call=$(printf "%s\n" "$bs_out" | awk 'NF>0 && $1 ~ /^[[:space:]]*[0-9]/ {print $(NF-2)}' | tail -n1)
  bs_put=$(printf "%s\n" "$bs_out"  | awk 'NF>0 && $1 ~ /^[[:space:]]*[0-9]/ {print $(NF-1)}' | tail -n1)

  # MC call
  mc_out="$("$MC_RUNNER" "$s0" "$k" "$rr" "$qq" "$ss" "$tt" "$SEED" 300000 "$BATCH" -1 1)"
  mc_price_call="$(extract_val "$mc_out" price)"
  lo="$(extract_val "$mc_out" ci_low)"; hi="$(extract_val "$mc_out" ci_high)"
  pass_call=$(awk -v p="$bs_call" -v lo="$lo" -v hi="$hi" 'BEGIN{print (p>=lo && p<=hi) ? "YES":"NO"}')

  printf "%-5s %-7s %-7s %-6s %-6s %-6s %-6s  %-11.6f %-11.6f %-9s\n" \
    "CALL" "$s0" "$k" "$rr" "$qq" "$ss" "$tt" "$bs_call" "$mc_price_call" "$pass_call"

  # MC put
  mc_out="$("$MC_RUNNER" "$s0" "$k" "$rr" "$qq" "$ss" "$tt" "$SEED" 300000 "$BATCH" -1 1 --put)"
  mc_price_put="$(extract_val "$mc_out" price)"
  lo="$(extract_val "$mc_out" ci_low)"; hi="$(extract_val "$mc_out" ci_high)"
  pass_put=$(awk -v p="$bs_put" -v lo="$lo" -v hi="$hi" 'BEGIN{print (p>=lo && p<=hi) ? "YES":"NO"}')

  printf "%-5s %-7s %-7s %-6s %-6s %-6s %-6s  %-11.6f %-11.6f %-9s\n" \
    "PUT"  "$s0" "$k" "$rr" "$qq" "$ss" "$tt" "$bs_put"  "$mc_price_put"  "$pass_put"
done

# --- 4) CSV de convergence ---------------------------------------------------
echo
echo "=== [4] CSV de convergence (pour tracé half_width vs 1/sqrt(N)) ==="
"$MC_RUNNER" "$S0" "$K" "$r" "$q" "$sigma" "$T" "$SEED" 1000000 50000 -1 1 --log --csv "$CSV_OUT"
echo "CSV écrit -> $CSV_OUT"
# Petit résumé: moyenne de c = half*sqrt(N) sur les points log
awk -F, 'NR>1 {c=$3*sqrt($1); s+=c; n++} END{ if(n>0) printf("c_avg (half*sqrt(N)) = %.8f over %d points\n", s/n, n); else print "no points"; }' "$CSV_OUT"

# --- 5) Moments du GBM -------------------------------------------------------
if [[ -x "$GBM_MOMENTS" ]]; then
  echo
  echo "=== Moments du GBM (N=200000) ==="
  "$GBM_MOMENTS" "$S0" "$r" "$q" "$sigma" "$T" "$SEED" 200000
else
  echo
  echo "[info] gbm_moments non trouvé. Ajoute l'exécutable si tu veux ce sanity check."
fi

echo
echo " Micro-expériences terminées."
