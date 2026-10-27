#!/usr/bin/env bash
# Small NiO AFM SOC+U fixture for A/B checks between qe-pert (QE 7.3) and this
# tree. Reduced cutoff and k mesh: a comparison fixture, not converged physics.
# Usage: bash dev-tools/soc_u_port_nio_fixture.sh OUT_DIR PW_X PH_X MODE [Q]
#   MODE = bare76  : 7.6 ph.x, lhubbard_soc=.true., only_init
#          bare73  : qe-pert ph_bare_validate.x, only_init
#          full76  : 7.6 ph.x, lhubbard_soc=.true., full phonon at Q
#          full73  : qe-pert ph_soc_u_q_operator_validate.x, full phonon at Q
#          efield76/efield73 : static epsilon and Born charges only (Gamma)
#          plain   : any ph.x, no DFPT+U input (use with SOC_U_NIO_NOU=1)
#          efieldc/fullc : collinear (SOC_U_NIO_SPIN=lsda) epsilon/Z* and phonon
#   SOC_U_NIO_SPIN = soc  (default: noncollinear + SOC, in-plane AFM, FR PPs)
#                    ncz  (noncollinear, no SOC, AFM moments along +z/-z, SR PPs)
#                    lsda (collinear nspin=2 AFM, SR PPs) -- collinear limit of ncz
#   Q defaults to "0.0 0.0 0.0". tr2_ph = $SOC_U_NIO_TR2 (default 1d-12): with
#   1d-14 the Hubbard relative-residual gate (sqrt(tr2_ph)) cannot be met in
#   irreps whose occupation response vanishes by symmetry at Gamma.
# MPI ranks/pools: SOC_U_NP (default 4), SOC_U_NK (default 2); one thread each.
set -euo pipefail

[[ $# -ge 4 ]] || { echo "Usage: bash $0 OUT_DIR PW_X PH_X MODE [Q]" >&2; exit 2; }
out=$1; pw=$2; ph=$3; mode=$4; q=${5:-0.0 0.0 0.0}
pp=${SOC_U_NIO_PP:-$HOME/PP/ONCVPseudoPack/Abinit_v0.4/UPF/PBE}
tr2=${SOC_U_NIO_TR2:-1d-12}
# Hubbard convergence test of the 7.6 runs: SOC_U_NIO_HUBCONV = relative
# (default, same gate as the qe-pert operator build), potential (7.6 default)
# or none
conv76="hubbard_conv='${SOC_U_NIO_HUBCONV:-relative}', "
# SOC_U_NIO_HUBTHR: hubbard_conv_thr of the 7.6 runs (default: not set, sqrt(tr2_ph))
[[ -n ${SOC_U_NIO_HUBTHR:-} ]] && conv76="${conv76}hubbard_conv_thr=${SOC_U_NIO_HUBTHR}, "
# SOC_U_NIO_NOSYM=1: no symmetry (isolates symmetrization differences)
# SOC_U_NIO_SHIFT=dx: move the first O by dx (alat) along x, removing the
#   inversion centre at Ni so that the bare occupation response is nonzero
nosym_line=''; [[ ${SOC_U_NIO_NOSYM:-0} == 1 ]] && nosym_line='nosym=.true.,'
ox=$(python3 -c "print(0.5 + float('${SOC_U_NIO_SHIFT:-0}'))")
# SOC_U_NIO_NOU=1: no Hubbard U (separates version effects from DFPT+U)
# SOC_U_NIO_DFT: input_dft override (e.g. PZ); default: functional of the PPs
# SOC_U_NIO_ECUTRHO: density cutoff in Ry (default 240 = 4 x ecutwfc)
ecutrho=${SOC_U_NIO_ECUTRHO:-240}
# SOC_U_NIO_ECUTWFC: wavefunction cutoff in Ry (default 60); SOC_U_NIO_K: k mesh (default '2 2 2')
ecutwfc=${SOC_U_NIO_ECUTWFC:-60}; kmesh=${SOC_U_NIO_K:-2 2 2}
dft_line=''; [[ -n ${SOC_U_NIO_DFT:-} ]] && dft_line="input_dft='${SOC_U_NIO_DFT}',"
# SOC_U_NIO_U: Hubbard U of Ni-3d in eV (default 5.7)
uval=${SOC_U_NIO_U:-5.7}
hubbard_card=$'HUBBARD {atomic}\n'"U Ni1-3d $uval"$'\n'"U Ni2-3d $uval"
case ${SOC_U_NIO_SPIN:-soc} in
  soc)  spin_lines="noncolin=.true., lspinorb=.true.,
  starting_magnetization(1)=0.5, starting_magnetization(2)=0.5,
  angle1(1)=90.0, angle2(1)=0.0, angle1(2)=90.0, angle2(2)=180.0,"; rel=FR ;;
  ncz)  spin_lines="noncolin=.true., lspinorb=.false.,
  starting_magnetization(1)=0.5, starting_magnetization(2)=0.5,
  angle1(1)=0.0, angle2(1)=0.0, angle1(2)=180.0, angle2(2)=0.0,"; rel=SR ;;
  lsda) spin_lines="nspin=2, tot_magnetization=0,
  starting_magnetization(1)=0.5, starting_magnetization(2)=-0.5,"; rel=SR ;;
  *) echo "unknown SOC_U_NIO_SPIN" >&2; exit 2 ;;
esac
[[ ${SOC_U_NIO_NOU:-0} == 1 ]] && hubbard_card=''
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
np=${SOC_U_NP:-4}; nk=${SOC_U_NK:-2}
mkdir -p "$out"; cd "$out"

cat > pw.in <<PWIN
&control
  calculation='scf', prefix='NiO', pseudo_dir='$pp', outdir='./out'
/
&system
  ibrav=0, celldm(1)=7.88, nat=4, ntyp=3, ecutwfc=$ecutwfc, ecutrho=$ecutrho,
  occupations='fixed', $spin_lines
  $nosym_line $dft_line
/
&electrons
  conv_thr=1.0d-10, mixing_beta=0.3, electron_maxstep=200
/
ATOMIC_SPECIES
Ni1 58.693 Ni_PBE_$rel.ABIv0.4.UPF
Ni2 58.693 Ni_PBE_$rel.ABIv0.4.UPF
O   15.999 O_PBE_$rel.ABIv0.4.UPF
ATOMIC_POSITIONS {alat}
Ni1 0.0 0.0 0.0
Ni2 0.5 0.5 0.0
O   $ox 0.0 0.0
O   1.0 0.5 0.0
CELL_PARAMETERS {alat}
1.0 0.5 0.5
0.5 1.0 0.5
0.5 0.5 1.0
K_POINTS {automatic}
$kmesh 0 0 0
$hubbard_card
PWIN

case $mode in
  bare76) extra="lhubbard_soc=.true., only_init=.true.," ;;
  bare73) extra="only_init=.true.," ;;
  full76) extra="lhubbard_soc=.true., ${conv76}tr2_ph=$tr2, alpha_mix(1)=0.3," ;;
  full73) extra="tr2_ph=$tr2, alpha_mix(1)=0.3," ;;
  efield76) extra="lhubbard_soc=.true., ${conv76}trans=.false., epsil=.true., zeu=.true., tr2_ph=$tr2," ;;
  efield73) extra="trans=.false., epsil=.true., zeu=.true., tr2_ph=$tr2," ;;
  plain) extra="tr2_ph=$tr2, alpha_mix(1)=0.3," ;;
  efieldc) extra="trans=.false., epsil=.true., zeu=.true., tr2_ph=$tr2," ;;
  fullc) extra="tr2_ph=$tr2, alpha_mix(1)=0.3," ;;
  *) echo "unknown MODE $mode" >&2; exit 2 ;;
esac
# SOC_U_NIO_PHEXTRA: extra &inputph text (e.g. "max_seconds=120, " or "recover=.true., ")
extra="$extra ${SOC_U_NIO_PHEXTRA:-}"
printf "NiO\n&inputph\n  prefix='NiO', outdir='./out', fildyn='NiO.dyn', verbosity='high',\n  %s\n/\n%s\n" \
  "$extra" "$q" > ph.in

[[ -f out/NiO.save/data-file-schema.xml ]] || \
  mpirun -np $np "$pw" -nk $nk -in pw.in > pw.out 2>&1
mpirun -np $np "$ph" -nk $nk -in ph.in > ph.out 2>&1
