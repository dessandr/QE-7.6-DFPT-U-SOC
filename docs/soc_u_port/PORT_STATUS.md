# Port status

Newest entry last. Plan: [PORT_PLAN.md](PORT_PLAN.md).

## 2026-10-23: phase 0 — preparation

- Source fixed: qe-pert branch `ph-acc` `1d326ca` (main `c00d268` + `ph_acc.x`).
- Tree: copied `~/bin/qe-7.6` without `build/` and without
  `external/{mbd,wannier90,devxlib}/.git`; `diff -rq` against the source is
  empty. Baseline commit `5e4b172`, all 8806 entries tracked.
- `make.inc`: `./configure --with-scalapack=intel --with-hdf5=<impi hdf5>`
  then HDF5_LIBS, libxc include/libs and `-D__LIBXC` copied from qe-pert's
  make.inc. Generated and used files in `../qe-7.6-pert-buildlogs/`.
- Baseline build `make -j4 pw ph`: success, 2 min 49 s, Intel ifx/mpiifx
  2026.1. Pristine executables now present in `bin/`.
- Cross-version fixture (`../qe-7.6-pert-buildlogs/xver_hbn_sr/`): hBN,
  PseudoDojo SR NC, U(N-2p)=2 eV, 4 ranks/2 pools. 7.3 vs 7.6: energy
  -26.67577207 vs -26.67576952 Ry; max frequency difference 0.034 cm^-1
  (q=(0.25,0,0)), 0.013 cm^-1 at Gamma. The FR-PP collinear fixture used for
  `ph_acc.x` does not run on 7.6 (`average_pp` error).
- 3-way dry run: 18 clean / 15 conflicting files. Per-file analysis in
  FILE_MAP.md; key 7.6 facts re-checked in source.
- Not done: no feature code ported, no SOC+U fixture prepared, no commit of
  the prep files (AGENTS.md, docs/soc_u_port, dev-tools script).
- Next: phase 1.

## 2026-10-23: decisions D1, D2

- D1 (user): the SOC+U path is switched on in `ph.x` by `lhubbard_soc=.true.`
  in `&INPUTPH`; default `.false.` keeps pristine 7.6 behavior. Semantics and
  guards in PORT_PLAN 4.1.
- D2 (user): magnetic static epsil/BEC is in scope. It works in qe-pert (one
  branch plus adjoint term in `dnsq_scf`; 7.3 had no magnetic E-field guard),
  but 7.6 rejects magnetic E-field, uses `nsolv=2` and writes one `dvpsi_e`
  branch, so it is re-implemented as phase 3b (PORT_PLAN 4.2).

## 2026-10-23: phase 1 — infrastructure

- Applied (qe-pert hunks onto the 7.6 owners): `ldaU_lr` kernels
  (`pack_hubbard_blocks_nc`, `build_dvhub_nc`, time-reversal helpers,
  `apply_hubbard_matrix_nc`) plus `accumulate_ns_nc_transition` and
  `ns_nc_hermiticity_residual` moved from `ldaU_ph`; `mix_pot.f90` (7.6 file
  identical to 7.3, so the merge equals the qe-pert file); `ph_restart`
  `HUBBARD_RESPONSE` checkpoint; `allocate_phq`, `close_phq`, `print_clock_ph`,
  `phq_setup` (`comp_dspinldau`; USE conflict resolved on the 7.6 list plus
  `d_spin_ldau`); `dynmatrix`/`q2qstar_ph`/`rotate_and_add_dyn` (antiunitary
  q-star; `q2qstar_ph` keeps a compatibility wrapper for elphon,
  elph_tetra_mod and q2qstar.x); `lr_orthoUwfc` (`ikk` records);
  `commutator_Vhubx_psi` (`intra_pool_comm` -> `intra_bgrp_comm`);
  `force_hub` oracle dump (USE conflict resolved on the 7.6 list).
- `mix_potential` default path (`nconv=ndim`, auxiliary flag true) is
  arithmetically unchanged; new behavior only for NaN, an unreleased history
  from a non-converged irrep, or a zero-norm secant.
- Build: `make -j4 pw ph hp` OK (39 s).
- Regression (`dev-tools/soc_u_port_regress.sh`, pristine binaries in
  `../qe-7.6-pert-buildlogs/pristine-bin`, compare with
  `soc_u_port_compare.py`): hBN SR+U Gamma and q, HP on hBN, test-suite
  `ph_U_insulator_us` (USPP+U, epsil, Gamma and q), `ph_base` Si (U=0):
  every compared value bitwise identical. Test-suite pseudopotentials for
  `hp_insulator_us` and `ph_Ni_nc_spinorbit_mag` could not be downloaded
  (network ~1-3 kB/s); local substitutes are the hBN HP case and a bcc Fe
  noncollinear SOC magnetic phonon (U=0).
- Found: qe-pert's base `347e370` is QE 7.3 plus the perturbo tutorial
  `dnsbare_pattern` write in `dnsq_bare.f90` (PH otherwise identical to
  `~/bin/q-e-qe-7.3`; that tree carries a separate ehub_uv patch in PW, so
  the PW comparison is inconclusive).

## 2026-10-23: phase 2 — bare response and `lhubbard_soc`

- `lhubbard_soc` (`ldaU_lr`, `&INPUTPH`, `bcast_ph_input`, `INPUT_PH.def`):
  noncollinear + U now stops unless it is set; with it, magnetic NCPP
  insulators only (atomic projectors, no orbital-resolved U, no shifted q
  mesh, phonons; epsil/zue blocked until phase 3b). Deliberate change from
  pristine 7.6: nonmagnetic noncollinear + U is rejected (7.6 Hubbard
  response routines have no spinor code). The 7.6 `magnetic_sym` + U guard
  is subsumed by this block.
- Macro policy: `dev-tools/soc_u_port_resolve_macros.py` resolves qe-pert
  macros (q-operator path defined, control twins and staged guard modes
  undefined); checked against Intel fpp on 9 qe-pert files (identical).
- Ported: `dnsq_bare` (qe-pert version + 7.6 changes: `proj1/proj2` from
  `ldaU_lr`, `dot_product` projections, `ikmks` from `qpoint`; qe2pert
  `dnsbare_pattern` write kept); `dvqhub_barepsi_nc`; `dvqpsi_us` dispatch
  (noncollinear -> `dvqhub_barepsi_nc`, branch and oracle flags
  `hub_bare_time_reversed`/`hub_bare_oracle` in `ldaU_lr`, set by
  `solve_linter` around its two `dvqpsi_us` calls; avoids the collinear term
  on spinors and double counting); `sym_dns_nc` in `LR_Modules/sym_dns.f90`
  on `upert`/`upert_mq`; `syme_dns` dropped; `dynmat_hub_bare_nc` and the
  `__PH_ACC` kernel (`intra_bgrp_comm`); `ph_acc.x` Make target.
- Build OK. Collinear regression (6 cases) bitwise identical to pristine.
  `ph_acc.x` vs 7.6 `ph.x` (hBN SR+U): d2nsbare max|diff| 9.8e-15 (Gamma),
  1.0e-15 (q); dynamical matrices and frequencies identical.
- NiO SOC+U bare A/B: first attempt stopped because the reduced-cutoff PW
  SCF stalled at 2.5e-11 Ry with conv_thr 1e-11; fixture now uses 1e-10.

## 2026-10-23/24: phase 3a (SCF response) and 3b (static E-field) — ported, validation in progress

- 3a ported: `LR_Modules/dnsq_scf.f90` noncollinear branch (two Sternheimer
  records, spinor `dpsi`, no triangle/1/2 for nc; `lflag` replaced by
  `ldaU_lr::hub_dnsq_single_branch`); `adddvhubscf_nc` + `response_kernels`
  dispatch (frequency and two-chem copies stop for noncollinear+U);
  `dynmat_hub_scf_nc`; `phqscf` export/restore; `dfpt_kernel`: joint
  Hxc+Hubbard Broyden (`mix_hubbard_response`) for noncolin+domag+U,
  restart rebuild, `mix_iter`, post-solve `dnsq_scf` refresh, time reversal
  of `v%ns_nc` around the -B solve (module procedure), optional
  `QE_SOC_U_DFPT_{FIRST,ALL}_ITER_THRESH` overrides applied on top of the
  7.6 adaptive threshold (`thresh_adaptive` keeps the 7.6 history).
  `make.depend`: `dfpt_kernels.o : mix_pot.o`. `only_init` restriction lifted.
- 3b ported: `dfpt_kernel(..., single_branch)` (nsolv=1 and the adjoint
  dN_scf term) passed by `solve_e` as `lhubbard_soc`; `phq_readin` magnetic
  E-field guard relaxed only for `lhubbard_soc`.
- Collinear/U=0/HP/Fe noncollinear-magnetic regression after 3a and after
  3b: bitwise identical to pristine 7.6.
- NiO SOC+U A/B (7.3 qe-pert vs this tree; reduced fixture
  `dev-tools/soc_u_port_nio_fixture.sh`): ground states differ by 4.9e-6 Ry
  (same |M|). Bare Hubbard dynamical matrix (symmetrized, Cartesian), Gamma
  and q=(0,0,0.5): max|diff| 3.2e-6 (rel 3.7e-7). Displacement patterns
  differ between versions, so pattern-basis quantities are not compared. In
  NiO the bare occupation response vanishes at the centrosymmetric Ni site
  for every q; a fixture with one O shifted is used to test it.
- Static E-field (tr2_ph 1e-12): epsilon max|diff| 1.0e-3 (rel 1.3e-4); the
  yz element differs by 0.2%, more than the other elements. Being checked
  with a nosym A/B (symmetrization differs: 7.3 psyme/syme_dns vs 7.6
  psymdvscf/sym_dns_nc with upert=s).
- Gamma phonon: with tr2_ph 1e-14 qe-pert itself stops unconverged (irrep 11;
  92 min). With 1e-12 qe-pert converges all 12 irreps (59 min) while the 7.6
  run stops at irrep 6 (Hubbard relative residual 2.6e-6 vs gate 1e-6 after
  150 iterations; absolute residual floor ~1e-7 in both versions). Irreps
  with zero Hubbard response converge in the same iteration counts in both
  versions. The gate sqrt(tr2_ph) on the relative Hubbard residual sits at
  the noise floor of this fixture — a qe-pert convergence-criterion issue
  carried by the port. Rerunning both with tr2_ph 1e-10 to compare
  Cartesian dynamical matrices.

## 2026-10-24: Hubbard convergence gate made optional (user decision)

- Default exit of the DFPT loop for `lhubbard_soc` is now the Hxc residual
  only, as in collinear DFPT+U (7.6 `dfpt_kernel`: dnsscf recomputed from the
  previous dpsi, never gated). The joint Hxc+Hubbard Broyden vector is kept;
  the Hubbard residual (`SOC_U_DNSSCF_CONVERGENCE_CHECK`) is still printed
  every iteration. `lhubbard_conv=.true.` (&INPUTPH, only with
  `lhubbard_soc`) restores the qe-pert default gate (relative Hubbard
  residual < sqrt(tr2_ph)). Documented in `INPUT_PH.def`. Reason: the gate
  sits at the noise floor (NiO irreps 3/6/11, absolute floor ~1e-7) and is
  undefined when the response vanishes by symmetry.
- The NiO fixture passes `lhubbard_conv=.true.` to 7.6 runs for parity with
  qe-pert (`SOC_U_NIO_CONV76=0` for the 7.6 default).
- Finding under investigation: qe-pert static E-field, nosym NiO: Hxc dr2
  reaches 7.7e-22 but the Hubbard residual stalls at ~2e-6 with a response
  scale of 5e-3..1.8e-2, although inversion symmetry at Ni requires the
  E-field occupation response to vanish; the symmetrized run gives ~1e-16.
  The 7.6 counterpart is running.

## 2026-10-24: second NiO batch (interrupted by the session end)

- Shifted-O NiO (first O +0.02 alat, inversion at Ni removed), bare, Gamma:
  dnsbare now nonzero (max 8.6e-3); 7.3 vs 7.6 max|diff| 8.9e-7 (rel
  1.0e-4, ground states differ by 4.9e-6 Ry); symmetrized bare Hubbard
  dynamical matrix rel 2.0e-7; Hermiticity residuals at round-off. Bare
  occupation response path agrees within version noise.
- nosym static E-field, both versions: Hxc dr2 reaches 1e-17..1e-20 while the
  E-field occupation response has scale ~5e-3 and its residual stalls
  (7.3: ~2e-6 after 150 iterations; 7.6: 1.3e-4 at iteration 61 when the
  batch was killed). Inversion at Ni requires this response to vanish; the
  symmetrized runs give ~1e-16. The anomaly is in the qe-pert E-field
  Hubbard path (single branch + adjoint term) and is reproduced by the port;
  symmetrization removes it. Not investigated further.
- Not run: Gamma phonon at tr2_ph 1e-10 (both versions).
- Relaunched (nohup, `../qe-7.6-pert-buildlogs/phase3-checks3.sh`):
  collinear regression after the `lhubbard_conv` build; 7.6 nosym E-field
  and Gamma phonon (tr2_ph 1e-12) with the new Hxc-only default, the latter
  compared with the converged qe-pert run at the same tr2_ph.

## 2026-10-24: Hxc-only exit rejected; potential-unit Hubbard test as default

- Evidence (NiO Gamma, tr2_ph 1e-12): with the Hxc-only exit the irreps with
  a Hubbard response (3, 6, 12) stop after 6-8 iterations with a relative
  Hubbard residual of 6-16% (U*|dN residual| ~ 1e-3 Ry); frequencies differ
  from the converged qe-pert run by up to 23 cm^-1. Unlike collinear
  DFPT+U, where dN_scf is recomputed unmixed from the previous dpsi every
  iteration, the joint mixing feeds a Broyden-mixed dN that lags the Hxc
  convergence. The 2026-10-24 Hxc-only default was therefore wrong.
- Retro-evaluation on the logged histories: the qe-pert relative test stops
  irreps 3/6/11 after 68-107 iterations or never (7.6 irrep 6, noise floor);
  max_I U_I*|dN_out-dN_in| <= sqrt(tr2_ph) is reached after 30-39 iterations
  (7.6 irrep 6: 71) in every run.
- Implemented (user decision): `hubbard_conv` in &INPUTPH (replaces the
  short-lived `lhubbard_conv`): 'potential' (default), 'relative' (qe-pert),
  'none' (diagnostic). `mix_hubbard_response` gained optional
  `hubbard_criterion` and `hub_pot_residual` (absent: unchanged relative
  test); new output line `SOC_U_HUBBARD_POTENTIAL_RESIDUAL`; the existing
  `SOC_U_DNSSCF_CONVERGENCE_CHECK` line is unchanged. NiO fixture: 7.6 runs
  use `SOC_U_NIO_HUBCONV` (default 'relative' for parity with qe-pert).
- 7.6 runs with the Hxc-only exit (now 'none'): nosym E-field converged in
  191 s, epsilon close to the symmetrized run (max|diff| 3.4e-3); Gamma
  phonon 851 s but inaccurate as above.
- Born effective charges of this reduced fixture violate the charge sum rule
  by ~0.78 in every run, so the 7.3/7.6 Z* difference (0.015) is below the
  fixture's resolution.
- Running: regression after the change; 7.6 E-field and Gamma phonon with the
  'potential' default vs the converged qe-pert runs
  (`../qe-7.6-pert-buildlogs/phase3-checks4.sh`).

## 2026-10-24: 'potential' default checked on NiO Gamma

- Regression after `hubbard_conv`: bitwise identical to pristine 7.6.
- 7.6 Gamma phonon, 'potential', tr2_ph 1e-12: all 12 irreps converge in
  36 min (Hubbard irreps 3/6/12: 35/71/37 iterations; qe-pert relative gate:
  59 min). Modes carried by the Hubbard irreps agree with the converged
  qe-pert run (320.408/320.665 vs 320.403/320.702 cm^-1; Hxc-only gave
  321.9/343.7). Static E-field: 'potential' and 'relative' give the same
  epsilon.
- Open: modes from irreps without Hubbard response still differ from qe-pert
  by up to ~5 cm^-1 (e.g. 56.9 vs 51.8, 353.9 vs 358.2), far above the hBN
  version difference (0.03 cm^-1). `lr_apply_time_reversal` flips only vrs
  and deeq_nc, so the explicit v%ns_nc flip in dfpt_kernel is not a double
  flip. Running the same NiO fixture with U=0 on qe-pert 7.3, pristine 7.6
  and the ported 7.6 (`phase3-checks5.sh`) to separate version/fixture
  effects from the port.

## 2026-10-24: U=0 control — the remaining 7.3/7.6 phonon difference is not the port

- Same NiO fixture with U=0 (Gamma, tr2_ph 1e-12; `phase3-checks5.sh`, run
  on pcs-gpu): pristine 7.6 ph.x vs ported ph.x on the same ground state:
  frequencies and dynamical matrix bitwise identical. qe-pert 7.3 vs
  pristine 7.6: max |dfreq| 5.57 cm^-1 (dyn rel 1.8e-3), ground states
  -736.31086915 vs -736.31087414 Ry, |M| 2.83 in both.
- With U=5.7 the qe-pert vs ported-7.6 difference is 5.05 cm^-1 (dyn rel
  5.0e-3), the same size. The ~5 cm^-1 is a 7.3 -> 7.6 noncollinear magnetic
  DFPT/fixture effect (acoustic modes down to -80 cm^-1: large ASR
  violation), not a porting error. Hubbard-carried modes agree with qe-pert
  to 0.04 cm^-1 with the 'potential' test.
- From here on, validation runs go to PBS compute nodes
  (`dev-tools/soc_u_port_qsub.sh`, MPI layout via SOC_U_NP/SOC_U_NK).
- Still open: the static E-field occupation response anomaly inherited from
  qe-pert (nosym), a q/=0 phonon A/B, restart and MPI-layout checks, qe2pert.

## 2026-10-24: E-field / Hubbard response investigation (user request)

Runs on pcs-gpu10 (PBS jobs 852, 853), `../qe-7.6-pert-buildlogs/efield-investigation/`.

- Parity: in the AFM NiO cell the three Gamma irreps with a Hubbard response
  (3, 6, 11) are the O1-O2 antiphase modes, even under inversion at Ni; all
  other Gamma modes are odd and give dN = 0. The uniform E-field is odd, so
  its dN must vanish; with nosym it appears at ~1e-5 in iteration 1
  (ground-state asymmetry) and grows to ~5e-3 while the Hxc residual
  converges.
- Conventions checked: PW `v_hubbard_nc` V(m1s1,m2s2) = -U N(m2s2,m1s1)
  (orbital and spin transpose) equals `build_dvhub_nc`; energy
  -U/2 sum N_ab N_ba consistent. No storage-convention mismatch found.
- Collinear limit (SOC off, SR Abinit v0.4 PPs, tr2_ph 1e-14, 'potential'):
  lsda (upstream collinear DFPT+U) and noncollinear moments +z/-z (ported
  path) give identical ground states (-736.18842988 Ry). E-field: lsda
  converges in 8 iterations with or without symmetry (epsilon sym/nosym
  diff 3.0e-5, Z* 1.3e-4); ncz converges in 7 with symmetry but fails
  after 150 without. ncz vs lsda (symmetric): epsilon diff 5.0e-3, Z* 3.1e-2,
  Gamma phonon up to 27 cm^-1 (acoustic) and 3-6 cm^-1 (optical), and the
  lsda doublet 253.300 splits into 252.1/254.5. The noncollinear path does
  not reproduce the collinear limit (pw: 12 symmetries lsda, 4 ncz).
- Step 2 submitted (job 853): U=0 lsda vs ncz phonon (pristine-7.6
  noncollinear magnetic DFPT, no DFPT+U code) and U=1e-5 eV with
  lhubbard_soc (E-field and phonon), to locate the deviation.
- Step 2 (job 853, PBE): U=0 lsda vs ncz Gamma phonon (pristine-7.6
  noncollinear magnetic DFPT, no DFPT+U code) differs by up to 11.5 cm^-1
  with broken degeneracies; U=1e-5 eV through lhubbard_soc: 9.7 cm^-1,
  E-field epsilon 23.69 vs 23.81 (0.5%). The symmetric-run discrepancy with
  the collinear limit therefore exists without the port.
- nosym ncz E-field (PBE, U=5.7): the spurious dN grows mainly in the
  transverse spin blocks (ud/du 4.2e-6 -> 7.4e-3 over 150 iterations; uu
  1.1e-4 -> 7.7e-4, dd ~2e-4). Interpretation: without SOC the AFM has an
  exact spin-rotation (Goldstone) mode, a zero mode of the noncollinear
  linear-response operator; ground-state asymmetry noise accumulates along
  it. lsda has no transverse channel; with symmetry the mode (inversion-even)
  is projected out of the odd E-field response; with SOC it is a near-zero
  mode (small anisotropy), consistent with the residual floors seen there.
  Not an implementation error of the E-field Hubbard path.
- Step 3 submitted (job 854): the same collinear-limit comparison with LDA
  (input_dft='PZ') to test whether the symmetric-run lsda/ncz differences come
  from the noncollinear GGA response.
- Step 3 (job 854, LDA input_dft='PZ', same PPs): lsda vs ncz epsilon diff
  1.6e-4 (PBE 5.0e-3); nosym ncz E-field converges in 15 iterations (PBE:
  fails after 150); Gamma phonon U=0 2.1 cm^-1 (PBE 11.5), U=5.7 3.0 cm^-1
  (PBE 27); Z* 0.022 (PBE 0.031). After ASR (dynmat.x asr='crystal'):
  U=0 2.56 cm^-1, U=5.7 1.99 cm^-1, with the lsda doublet 156.32 split into
  154.38/156.60 at U=0 already.
- Conclusions: (1) the nosym E-field non-convergence is a property of the
  noncollinear GGA transverse response (amplification along the spin-rotation
  direction), not of the E-field Hubbard path; with LDA it converges.
  (2) Most of the symmetric-run deviation from the collinear limit is the
  noncollinear GGA response; the remaining ~2 cm^-1 (LDA, after ASR) is
  present in the pristine-7.6 noncollinear magnetic phonon path at U=0, so
  it is not introduced by the port, and the DFPT+U part adds no further
  deviation. (3) The Z* difference (0.02) is not resolved; this fixture's Z*
  sum-rule violation is ~0.8.
- Implication: qe-pert production runs use PBEsol (GGA) noncollinear AFM;
  their Hubbard residual floors and the open "occupation FD mismatch" item
  may share this GGA transverse-response origin (untested).
- GGA transverse kernel (code reading, LR_Modules): `setup_dgc` builds the
  local-frame densities with PW's `compute_rho` (same `ux`/`segni`), and
  `dgradcorr` adds 0.5*segni*vsgga*(dm - m^(m^.dm))/|m| (plus 1/|m|^3
  gradient terms), cut only at |m| > 1e-12. vsgga = V_up - V_dn of the GGA
  part. This is the exact second derivative of the local-frame functional
  E[n,|m|]; near AFM nodes the GGA spin potential need not follow m, so the
  transverse curvature can be large or negative (hypothesis, untested).
- SOC (final goal) check by frozen-phonon FD (PBE, SOC, U=5.7, Gamma,
  ecutrho 240, d=0.001/0.002 alat, Richardson): inversion-even O-antiphase
  pattern: DFPT vs FD rel 1.6% (port) / 2.2% (qe-pert); odd O-in-phase
  control: 7.9% / 6.4%, where the DFPT forces violate the sum rule (-0.017
  Ry/bohr^2) while FD forces sum to ~1e-4. At this cutoff the comparison is
  limited by the DFPT acoustic-sum-rule error, not by the transverse issue.
- Submitted (job 857): the same FD vs DFPT comparison at ecutrho 480 Ry for
  PBE and LDA.
- FD vs DFPT at ecutrho 480 (job 857, SOC, U=5.7, Gamma, ported ph.x):
  inversion-even (Hubbard-active) pattern: PBE 1.7% (240: 1.6%), LDA 0.14%.
  This measures the noncollinear-GGA effect on the Hubbard-active modes of
  this fixture (~1.5%). Odd O-in-phase pattern: O rows PBE 0.3-1.6%, LDA 16%;
  Ni rows far off in both (PBE 18%, LDA 66% max), DFPT sum-rule violation
  0.030 (PBE) / 0.069 (LDA) Ry/bohr^2 not reduced by ecutrho 240->480, and
  Ni1/Ni2 rows not equivalent (LDA -0.046 vs -0.057) while FD rows are.
  The collinear lsda LDA U=0 phonon of the same fixture also has acoustic
  modes at -44 cm^-1, so this odd-mode/sum-rule problem is present in
  pristine QE paths; both PPs have NLCC (core_correction=T). Not resolved.
- Energy-curvature cross-check discarded: u.C.u from total energies changes
  0.42 -> 0.35 between ecutrho 240 and 480 while the force-based value stays
  0.363 (energy differences ~1e-5 Ry are dominated by egg-box/SCF noise).

## 2026-10-25: noncollinear GGA xc kernel tested directly; two upstream bugs fixed

Tool `dev-tools/xckfd/` (`xckfd.x`, test use only): reads a pw.x save, sets
up the kernel as ph.x does (compute_ux, setup_dmuxc, setup_dgc) and compares
`dv_of_drho_xc` with the central difference of the pw.x `v_xc` (Richardson,
h = 1e-4) for q = 0 perturbations: charge dn = n, longitudinal dm = m, global
rotation dm = a x m (exact answer a x B_xc), modulated rotation
cos(2 pi i1/nr1) a x m. Runs and outputs:
`../qe-7.6-pert-buildlogs/efield-investigation/xckfd/` (NiO fixture, SOC PBE
`nio/pw76`, ncz PBE `pwz`, ncz LDA `lda/pwz`; 1 core, seconds).

- LDA (ncz): kernel = FD to 1e-10..1e-13 for every perturbation (control).
- Bug A, ph.x ground-state B_xc: `post_xml_init` (read_file) runs `v_of_rho`
  before anything sets the fixed GGA axis (`lsign`, `ux`, set only by
  `compute_ux` in PW `setup` and later in `phq_setup`). ph.x therefore keeps
  in `v%of_r`/`vrs` a B_xc built from |m| (lsign = .false.), while the pw.x
  SCF and the DFPT kernel use sign(m.ux)|m|. In the continuum the two agree
  (up/down swap); on the grid |m| has a kink at the AFM nodal surfaces and the
  FFT gradient differs: B_xc(lsign=F) - B_xc(lsign=T) = 57% (L2) of B_xc,
  max 0.17 Ry, located at the nodes (|m| ~ 1e-4), for both SOC and ncz PBE;
  E_xc differs by 2e-4 Ry. So at Gamma ph.x solves the Sternheimer equations
  with an H0 that is not the pw.x Hamiltonian of its wavefunctions, and the
  kernel's rotation response (a x B_xc(T), exact to 3e-13) does not match
  the H0 field B_xc(F) (50% off): the spin-rotation (Goldstone) condition is
  broken. With LDA lsign is irrelevant. The q /= 0 nscf path inherits
  whatever lsign is set at that point. Present in pristine 7.6.
- Bug B, `dgradcorr` noncollinear branch: the gradient of
  delta(rho +- s|m|)/2 was built by the chain rule from FFT gradients of m
  and dm, with `gc_lr%gmag` (grad m) allocated, zeroed and never computed
  (also in pristine 7.6). For collinear m (ncz) the missing terms cancel
  (errors 1e-5); with SOC canting the kernel errs by 6% (longitudinal,
  rotation, modulated rotation). Filling gmag fixes the rotation but not the
  longitudinal (6%) and modulated (1.6e-3) cases: near nodes s|m| is not
  smooth under canting, so the chain rule with FFT gradients is not the
  derivative of the discretised v_xc. Taking the FFT gradient of
  delta(rho +- s|m|)/2 itself (as PW gradcorr does for rho +- s|m|) matches
  FD to 6e-10..9e-10 (rotations) and 1.5e-5 (longitudinal, as in ncz).
- Remaining common GGA difference: charge channel 6.4e-4 (L2, spikes to
  5% of max), identical in SOC and ncz and for all variants, effect on the
  quadratic form drho.dv 1e-5. Not investigated further.
- SOC FD at h = 1e-3 fails for rotations (h vs 2h 46%) because sign(m.ux)
  flips at ~1e-3-rad rotations on near-node points (33892 sign-change grid
  links, 4244 with |m_perp| > |m.ux|); at h = 1e-4 FD converges (1e-9).
- Fixes (port tree, uncommitted): A: `PW/src/read_file_new.f90`
  `post_xml_init` computes ux/lsign from the saved starting magnetization
  and angles (as `setup`) before `v_of_rho`, only for noncolin + domag + GGA
  (local array, `m_loc` untouched). pw.x binary unchanged (SHA256 equal;
  pw.x does not link it). B: `LR_Modules/dgradcorr.f90` noncollinear branch
  takes `fft_qgradient` of drhoout; `gdmag`/`gmag` no longer used there
  (still allocated by setup_dgc). After the fixes: B in `v%of_r` from
  read_file = B_xc(lsign=T) to 3e-16; library kernel = FD as above.
  Pre-fix binaries kept in `../qe-7.6-pert-buildlogs/prefix-bin/` (SHA256).
- Validation submitted: PBS 860 (pcs-gpu06; ncz PBE collinear limit:
  nosym/sym E-field, U=0 and U=5.7 Gamma phonons vs the earlier lsda runs)
  and 861 (pcs-gpu07; SOC nosym E-field fixed vs pre-fix, Gamma phonon at
  ecutrho 240/480 vs FD). Regression set run locally (4 ranks, 2 pools).
- Regression set after the fixes (local, 4 ranks, 2 pools,
  `regress/ncgga-fix`): collinear +U, HP, U=0 cases bitwise identical to the
  previous run; `fe_nc_soc_mag` (Fe SOC GGA, q = (0,0,1/2), nscf path)
  converges (3 irreps, 133 iterations) and changes as expected:
  415.087 -> 414.212 and 453.631 -> 453.375 cm^-1 (the first attempt was
  killed at the session end; rerun complete).
- ncz PBE collinear limit, static E-field with symmetry (job 860): epsilon
  vs lsda 5.0e-3 -> 3.2e-4 (LDA 1.6e-4); Z* 0.022 -> 0.023 (LDA 0.013,
  the separate Z* sum-rule issue). nosym E-field: still not converged after
  150 iterations, but the spurious dN ends at 5.2e-3 (pre-fix 3.5e-2, LDA
  converged run 3.2e-3) and the Hubbard potential residual plateaus at 3e-5
  instead of growing (pre-fix 1.1e-4).
- Transverse kernel of the local-frame GGA (xckfd `kperp_stats`): with
  E_xc[n, s|m|] the transverse response is dB_perp = K dm_perp with
  K = (B_xc.m)/|m|^2. LDA (ncz): K in [-7.6, -0.11] Ry, Stoner sign
  everywhere. PBE (ncz / SOC): K in [-3.5e3, +1.3e4] / [-3.5e3, +3.3e3] Ry
  against ~2.8 Ry inside the moments; B.m > 0 on 13% of the grid (16% of
  the charge) and |K| > 10x the moment scale on 13% (3% of the charge). The
  GGA spin field vsgga does not vanish on the nodal surfaces, so K ~ 1/|m|
  there with either sign. This is a property of the pw.x noncollinear GGA
  functional, not of the DFPT code: the kernel is its exact derivative (FD
  above). It explains why the nosym E-field (exact spin-rotation zero mode
  without SOC) still does not converge with PBE.
- Collinear limit after the fixes (job 860, PBE, Gamma, crystal ASR via
  dynmat.x, max |ncz - lsda|): U=0 4.62 -> 2.61 cm^-1 (LDA 2.56);
  U=5.7 (lhubbard_soc) 5.66 -> 2.18 cm^-1 (LDA 1.99), and the lsda doublet
  327.01 (pre-fix 321.35/321.99) is degenerate again (325.72). The
  GGA-specific part of the deviation is gone; the ~2 cm^-1 left is the
  pristine noncollinear-path residue already seen with LDA.
- SOC FD vs DFPT after the fixes (job 861, U=5.7, Gamma): even
  (Hubbard-active) pattern ecutrho 240: 1.64% -> 0.53%; 480: 1.71% -> 1.06%
  (LDA 480: 0.14%). Odd pattern 240: 7.9% -> 6.3%; 480: 18.4% -> 11.6%
  (the separate sum-rule/Ni-row problem, also with LDA).
- SOC nosym static E-field, same settings fixed vs pre-fix ph.x: neither
  converges in 150 iterations (target 1e-6); Hubbard potential residual
  plateau ~1e-3 (fixed) vs ~4e-5 (pre-fix), dN 3.2e-2 vs 1.4e-2. With SOC
  canting s|m| jumps at the nodes, and the exact kernel of that functional
  (singular, sign-indefinite K_perp) conditions the nosym problem worse than
  the inconsistent pre-fix kernel did. Symmetric runs converge.
- Open decision (user): whether to regularise the noncollinear GGA
  functional itself (pw.x and ph.x together) near |m| -> 0.
- Decision D6 (user, 2026-10-26): the noncollinear GGA functional is left as
  is (no regularisation); DFPT is run with symmetry, nosym is not a target.
  The fixes A and B stay in the tree (uncommitted).

## 2026-10-26: odd-mode / acoustic-sum-rule error isolated (in progress)

Runs: `../qe-7.6-pert-buildlogs/efield-investigation/nlcc/` (PBS 865-873).
Tools: `compare_fd_dfpt2.py` (FD from raw pw.x forces, i.e. the sum of the
verbosity='high' contributions, because pw.x subtracts the mean force and
symmetrises; `--asr` prints max|sum_b C_ab|), `fd-gen.sh` (FD for any spin
mode, patterns even/odd/ni, FD_CONV/FD_ECUTWFC/FD_K), NLCC-free copies of the
PPs in `../qe-7.6-pert-buildlogs/pp-nonlcc` (core_correction="F", test only).
Fixture: new optional SOC_U_NIO_ECUTWFC and SOC_U_NIO_K (defaults unchanged).

- NLCC is not the cause: lsda LDA U=0 ASR violation 3.6% (NLCC) vs 2.7%
  (none); lsda PBE 0.43% vs 6.4% (without NLCC the GGA egg-box is large: SOC
  PBE no-NLCC FD raw net force 0.029). Not the grid: lsda LDA 240 -> 480 Ry
  2.7 -> 3.1%. Not DFPT+U: SOC PBE U=0 7.5%.
- The pw.x energy is translation invariant to ~1e-3 (derivative of the raw
  net force) in the LDA and NLCC-PBE cases, so a DFPT net force of 0.07
  (SOC LDA) is a DFPT error.
- Inversion-even pattern: DFPT = FD to 0.08% (SOC LDA no NLCC), 0.27% (SOC
  PBE no NLCC), 0.5% (SOC PBE NLCC). The errors are in the Ni-displacement
  response (Ni-Ni blocks, and C(Ni,O) through hermiticity): lsda LDA U=0 Ni1
  self term 0.22585 vs FD 0.23007 (1.8%; FD stable to 1e-4 between
  conv_thr 1e-10 and 1e-13); ncz LDA O patterns fine (0.3-0.5%).
- lsda: sym vs nosym DFPT (identical ground states) differ by up to 2.9e-3
  in the Ni-Ni blocks; SOC by up to 4.5e-2 (Ni1-Ni1). SOC nosym fixes the
  O rows of the odd pattern (0.1%) and the Ni1/Ni2 asymmetry, Ni rows still
  8% off.
- Cause found for lsda: incomplete DFPT SCF convergence. tr2_ph 1d-14 reports
  convergence, but the Ni self term is 1.8% off; with tr2_ph 1d-18 (same
  symmetric run, 5-6 more iterations per irrep) it is 0.2% and the ASR
  violation 0.5%. The force-constant formula is not variational, so the
  residual dvscf error enters C linearly. The sym/nosym differences are
  consistent with different convergence paths.
- Submitted: SOC U=5.7 LDA (no NLCC) and PBE (NLCC), and SOC PBE U=0, all
  at tr2_ph 1d-18, against the raw-force FD.
- lsda checks at tr2_ph 1d-14 (job 870): ecutwfc 90 (2.2%) and k 4x4x4
  (3.1%) do not remove the Ni self-term error; only tr2_ph does.
- SOC PBE U=0 (pristine noncollinear path, NLCC) at tr2_ph 1d-18 (job 873):
  ASR violation 7.5% -> 0.6%; DFPT net force of the odd pattern 0.0021 =
  FD raw 0.0021; Ni1/Ni2 rows equal; FD agreement even 0.8%, odd 1.9% (a
  uniform ~2% scale difference; the SOC+GGA functional is not smooth near
  canted nodes, see 2026-10-25).
- SOC PBE U=5.7 (NLCC) at tr2_ph 1d-18: with hubbard_conv='potential' the
  Hubbard-active irreps cannot converge (job 872: potential residual
  plateaus at 1.46e-8 against sqrt(tr2_ph) = 1e-9). With 'none' (job 874):
  ASR 9.5% -> 1.6%, odd 6.4% -> 2.9%, but even 0.5% -> 2.0% (Hxc-only exit,
  dN lags as before).
- SOC LDA U=5.7 (no NLCC) at 1d-18, 'potential' (job 871): 11 of 12 irreps
  converge; irreps 3 and 6 just pass 1e-9 after ~100 iterations, irrep 12
  stays at 1.09e-8 with the Hxc residual at 9.7e-25. The Hubbard response
  has a noise floor ~1e-8 independent of the Hxc convergence, so tying the
  Hubbard threshold to sqrt(tr2_ph) fails at tight tr2_ph.
- Conclusions: (1) most of the odd-mode / sum-rule error was incomplete DFPT
  convergence at tr2_ph 1d-12..1d-14, also in pristine collinear QE; this
  fixture needs tr2_ph ~1d-18 for Ni-displacement irreps. (2) The Hubbard
  convergence threshold needs to be separable from tr2_ph (open decision).
- New input `hubbard_conv_thr` (ldaU_lr, phq_readin, bcast_ph_input,
  dfpt_kernels, INPUT_PH.def; user decision 2026-10-26: separate Hubbard
  convergence as in the original design): threshold of the Hubbard test,
  default (<= 0) sqrt(tr2_ph). Regression set with the default: all 25
  compared files identical (`regress/hubthr`). Fixture: SOC_U_NIO_HUBTHR.
  With 1d-8 the Hubbard-active irrep 12 still stalls (LDA 1.09e-8, PBE
  8.3e-8, slowly creeping); 1d-7 converges every irrep. The stall looks like
  joint-Broyden stagnation of the small Hubbard block once the Hxc part is
  far below its threshold (CG threshold tracks sqrt(dr2) ~ 1e-13), not CG
  noise; mixer improvement left for later.
- Converged references (tr2_ph 1d-18, hubbard_conv_thr 1d-7 where U /= 0),
  FD from raw forces:
    lsda PBE NLCC U=0 : even 2.0e-4, odd 6.2e-4, Ni 2.8e-4; ASR 0.2%
    ncz  PBE NLCC U=0 : even 9.4e-4, odd 6.9e-4, Ni 4.2e-4; ASR 0.3%
    SOC  LDA U=5.7    : even 8.2e-4, odd 1.1e-3 (mean-removed 6.4e-4); ASR 0.46%
    SOC  PBE NLCC U=0 : even 0.8%, odd 1.9%; ASR 0.6%
    SOC  PBE NLCC U=5.7: even 2.0%, odd 2.9% (Ni rows 0.5%, O rows 2.9%); ASR 1.6%
  So the collinear and noncollinear GGA paths (with fixes A and B) and the
  ported SOC+U path are exact once converged; the ~2% left is SOC + GGA
  only. FD smoothness (d vs 2d) for SOC PBE 1.4-3.7e-4 against 5-8e-5 for
  the others, ten times below the difference. Interpretation: with SOC
  canting, sign(m.ux) flips on near-node grid points under a displacement;
  FD includes the resulting jumps, the DFPT linearisation (segni fixed)
  cannot. A property of the pw.x noncollinear GGA functional (decision D6:
  left as is).
- Practical settings for this kind of system: tr2_ph ~1d-18 (the default
  1d-12 leaves percent-level errors in the Ni-displacement irreps, also in
  pristine collinear QE) and hubbard_conv_thr ~1d-7 (or 1d-6; the even modes
  were already accurate at 1d-6).

## 2026-10-26: general-q validation (in progress)

- User decision: Born effective charges are left for now (to be handled by
  zasr post-processing later); next is general q.
- Method: supercell frozen-phonon FD (`../qe-7.6-pert-buildlogs/efield-
  investigation/qpt/fd-sc.sh`, `compare_sc.py`). Supercell N x a1, one atom of
  cell 0 displaced; C_{b,ax}(q = m b1/N) = sum_c X_b(c) exp(-2 pi i m c/N)
  with X = -dF/du from raw forces; the supercell k mesh 1 x K x K is the
  primitive K x K x K mesh folded, so the comparison is exact.
  Case: SOC, LDA (PZ), no NLCC, U=5.7, tr2_ph 1d-18, hubbard_conv_thr 1d-7.
  q = b1/2 = (0.75,-0.25,-0.25) (N=2, k 2x2x2; m=0 also re-checks Gamma)
  and the general q = b1/3 (N=3, k 3x3x3; complex dynamical matrix).
  PBS 880-883.
- Supercell FD method check at Gamma (N=2, m=0) against the converged Gamma
  DFPT (nlcc/tr18_lda_off_thr7): Ni1 column 1.8e-3, O1 column 7.6e-4.
- q = b1/2 (job 881, 504 iterations, 56 min): DFPT vs supercell FD, Ni1
  column 1.3e-3, O1 column 1.0e-3 (relative max). Zone-boundary q path (nscf
  at k+q, time-reversed -k/-k-q branches, Hubbard at q) agrees with FD as
  well as Gamma does.
- Restart (job 884): the q = b1/2 run stopped by max_seconds three times
  inside irreps 4, 8 and 12 and continued with recover=.true.; final dyn vs
  the uninterrupted run max|dC| 5.4e-7 (rel 1.4e-6), frequencies within
  4.3e-4 cm^-1. Fixture: new SOC_U_NIO_PHEXTRA (extra &inputph text).
- MPI layout (job 885, Gamma, tr2_ph 1d-14): ranks/pools 32/4, 16/2, 32/2
  give optical modes within 0.03 cm^-1; differences up to 1.6 cm^-1 only in
  the acoustic modes (-31..-5 cm^-1, ASR 1.5% at this tr2_ph) and
  max|dC| 1.8-2.7e-4 between any two runs, including the two nk=2 runs:
  convergence noise at 1d-14, not a layout dependence.
- General q = b1/3 = (0.5,-1/6,-1/6) (job 883, k 3x3x3, 72 k points in the
  ph.x run; supercell N=3, k 1x3x3, 12 atoms, job 882): complex dynamical
  matrix, DFPT vs FD Ni1 column 1.2e-3, O1 column 4.1e-3 (relative max);
  with the opposite phase convention 23% / 48%, so the phases are tested.
  The -q block of the dyn file (star by the antiunitary operation) matches
  the FD at -q to the same accuracy.
- Summary of the SOC+U DFPT validation (NiO fixture, SOC LDA no NLCC, U=5.7,
  tr2_ph 1d-18, hubbard_conv_thr 1d-7): Gamma, q = b1/2 and q = b1/3 agree
  with frozen-phonon FD to ~1e-3; restart and MPI layout verified. Open:
  Born charges (to be handled with zasr post-processing), qe2pert port,
  joint-mixer stagnation of the Hubbard block (workaround hubbard_conv_thr).
- Gamma at k 3x3x3 (job 883) vs the N=3 supercell FD (m=0): Ni1 column
  2.5e-3, O1 column 6.9e-3, the same level as q = b1/3.
- 2026-10-27: committed on branch soc-u-port (9ad5313 port, 20a3e47
  noncollinear GGA fixes, fc3a8b4 tools/docs) and pushed with the pristine
  baseline main to github.com:dessandr/QE-7.6-DFPT-U-SOC (user request).
