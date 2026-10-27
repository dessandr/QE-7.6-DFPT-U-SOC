#!/usr/bin/env bash
# Regression set for the SOC+U port: paths that must stay identical to the
# pristine QE 7.6 build (collinear +U, HP, U=0, noncollinear magnetic U=0).
# Usage: bash dev-tools/soc_u_port_regress.sh BIN_DIR OUT_DIR [CASE...]
#   BIN_DIR holds pw.x, ph.x, hp.x (e.g. ../qe-7.6-pert-buildlogs/pristine-bin
#   or the tree's bin/). Compare two OUT_DIRs with soc_u_port_compare.py.
# MPI ranks/pools: SOC_U_NP (default 4, local five-core limit), SOC_U_NK
# (default 2); one thread each. Compare only runs made with the same layout.
# Exit status 1 if any case fails.
set -euo pipefail

[[ $# -ge 2 ]] || { echo "Usage: bash $0 BIN_DIR OUT_DIR [CASE...]" >&2; exit 2; }
bin=$(cd -- "$1" && pwd); out=$2; shift 2
tree=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
ts=$tree/test-suite
sr_pp=${SOC_U_REGRESS_SR_PP:-$HOME/PP/ONCVPseudoPack/PseudoDojo/SR_v0.4.1/PBE_standard}
fr_pp=${SOC_U_REGRESS_FR_PP:-$HOME/PP/ONCVPseudoPack/PseudoDojo/FR_v0.4/PBE_standard}
# Pseudopotentials the QE test-suite normally downloads.
ts_pp=${SOC_U_REGRESS_TS_PP:-$tree/../qe-7.6-pert-buildlogs/testsuite-pseudo}
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
np=${SOC_U_NP:-4}; nk=${SOC_U_NK:-2}
run() { mpirun -np $np "$bin/$1" -nk $nk -in "$2" > "${2%.in}.out" 2>&1; }

# Copy test-suite inputs with local outdir and absolute pseudo_dir.
stage() {  # case_dir files...
  local src=$ts/$1; shift
  for f in "$@"; do
    sed -e "s#outdir *= *'[^']*'#outdir='./out'#" \
        -e "s#pseudo_dir *= *'[^']*'#pseudo_dir='$ts_pp/'#" "$src/$f" > "$f"
  done
}

hbn_inputs() {  # q vector
  cat > pw.in <<PWIN
&control
  calculation='scf', prefix='bn', pseudo_dir='$sr_pp', outdir='./out'
/
&system
  ibrav=4, celldm(1)=4.7419, celldm(3)=4.3, nat=2, ntyp=2, ecutwfc=40
/
&electrons
  conv_thr=1d-12
/
ATOMIC_SPECIES
N 14.00674 N.upf
B 10.81    B.upf
ATOMIC_POSITIONS {alat}
N 0.000000000 -0.288675135 0.000000000
B 0.000000000  0.288675135 0.000000000
K_POINTS {automatic}
4 4 1 0 0 0
HUBBARD {atomic}
U N-2p 2.0
PWIN
  printf "bn\n&inputph\n  tr2_ph=1d-14, prefix='bn', outdir='./out', fildyn='bn.dyn', verbosity='high'\n/\n%s\n" "$1" > ph.in
}

case_hbn_sr_gamma() { hbn_inputs "0.0 0.0 0.0"; run pw.x pw.in; run ph.x ph.in; }
case_hbn_sr_q() { hbn_inputs "0.25 0.0 0.0"; run pw.x pw.in; run ph.x ph.in; }

case_hp_hbn_sr() {  # collinear HP; exercises mix_potential through hp.x
  hbn_inputs "0.0 0.0 0.0"
  printf "&inputhp\n  prefix='bn', outdir='./out', nq1=1, nq2=1, nq3=1,\n  conv_thr_chi=1d-8, iverbosity=2\n/\n" > hp.in
  run pw.x pw.in; run hp.x hp.in
}

case_fe_nc_soc_mag() {  # noncollinear magnetic SOC phonon, U=0 (antiunitary q-star)
  cat > pw.in <<PWIN
&control
  calculation='scf', prefix='fe', pseudo_dir='$fr_pp', outdir='./out'
/
&system
  ibrav=3, celldm(1)=5.42, nat=1, ntyp=1, ecutwfc=35,
  occupations='smearing', smearing='mv', degauss=0.02,
  noncolin=.true., lspinorb=.true., starting_magnetization(1)=0.5
/
&electrons
  conv_thr=1d-11
/
ATOMIC_SPECIES
Fe 55.845 Fe.upf
ATOMIC_POSITIONS {alat}
Fe 0.0 0.0 0.0
K_POINTS {automatic}
4 4 4 1 1 1
PWIN
  printf "fe\n&inputph\n  tr2_ph=1d-12, prefix='fe', outdir='./out', fildyn='fe.dyn'\n/\n0.0 0.0 0.5\n" > ph.in
  run pw.x pw.in; run ph.x ph.in
}

case_ph_U_insulator_us() {
  stage ph_U_insulator_us BN.scf.in BN.phG.in BN.phq.in
  run pw.x BN.scf.in; run ph.x BN.phG.in; run ph.x BN.phq.in
}
case_ph_base_si() {
  stage ph_base si.scf.in si.phG.in si.phX.in
  run pw.x si.scf.in; run ph.x si.phG.in; run ph.x si.phX.in
}
# The next two need test-suite pseudopotentials that may be missing locally;
# run them by name when ts_pp has them.
case_hp_insulator_us() {
  stage hp_insulator_us LiCoO2.scf.in LiCoO2.hp.in
  run pw.x LiCoO2.scf.in; run hp.x LiCoO2.hp.in
}
case_ph_Ni_nc_spinorbit_mag() {
  stage ph_Ni_nc_spinorbit_mag Ni.scf.in Ni.ph.in
  run pw.x Ni.scf.in; run ph.x Ni.ph.in
}

all_cases=(hbn_sr_gamma hbn_sr_q hp_hbn_sr ph_U_insulator_us ph_base_si
           fe_nc_soc_mag)
[[ $# -gt 0 ]] && all_cases=("$@")
mkdir -p "$out"
sha256sum "$bin"/{pw,ph,hp}.x > "$out/BINARIES.sha256"
failed=0
for c in "${all_cases[@]}"; do
  rm -rf "${out:?}/${c:?}"; mkdir -p "$out/$c"
  if ( cd "$out/$c" && "case_$c" ); then echo "done: $c"
  else echo "CASE FAILED: $c" >&2; failed=1; fi
  rm -rf "${out:?}/${c:?}/out"
done
exit $failed
