# Porting qe-pert DFPT+SOC+U to QE 7.6 — plan

Prepared 2026-10-23. Nothing has been ported yet; this tree is still pristine
QE 7.6 (`5e4b172`) plus port tooling. Per-file hunk mapping:
[FILE_MAP.md](FILE_MAP.md). Progress log: [PORT_STATUS.md](PORT_STATUS.md).

## 1. What is being ported

- Source: `~/bin/qe-pert`, branch `ph-acc` (`1d326ca`), diff from `347e370`
  (pristine QE 7.3). QE-core scope only: 33 files in `LR_Modules`,
  `PHonon/PH`, `PW/src`. `perturbo/qe2pert` is a later, separate step.
- The working SOC+U path in qe-pert is **not** the ordinary `ph.x`: it is the
  `__SOC_U_DFPT_Q_OPERATOR_VALIDATE` build (`ph_soc_u_q_operator_validate.x`)
  plus guard relaxations (`__SOC_U_DFPT_{FULL,BARE,Q_BARE}_VALIDATE`). The
  other `__SOC_U_*_VALIDATE` macros are control/diagnostic twins.
- qe-pert's open issues move with the code (independent occupation-FD
  mismatch, general-q SOC-off gate failure, q-star covariance residuals; see
  `~/bin/qe-pert/docs/soc_u_dfpt/AGENT_HANDOFF_KNOWN_ISSUES.md` section 3).
  Agreement with qe-pert is a porting check, not physics validation.

## 2. Mechanical picture

`dev-tools/soc_u_port_merge_dryrun.sh` (3-way: 7.3 / qe-pert / 7.6):
18 of 33 files merge without textual conflict, 15 conflict
(`solve_linter` 13, `solve_e` 5, `dnsq_bare` 4, `dnsq_scf` 4, others 1-2).
A clean textual merge is not a correct port: several clean files hit 7.6
API changes (communicators, relocated callers).

## 3. QE 7.3 -> 7.6 changes that shape the port (all verified in 7.6 source)

| 7.6 change | Consequence for the port |
| --- | --- |
| PH SCF loop moved to `MODULE dfpt_kernels` (`LR_Modules/dfpt_kernels.f90`), shared by `solve_linter`, `solve_e`, TDDFPT | qe-pert's `solve_linter`/`solve_e` hunks are re-implemented in `dfpt_kernel`. Per-executable `-D` recompilation cannot reach a module in `libqelr`: control must become runtime flags (decision D1). |
| `dvqpsi_us.f90:191` calls the collinear `dvqhub_barepsi_us` for any `lda_plus_u` | Without a noncolin dispatch there, spinor runs get the collinear bare term, and adding qe-pert's explicit `_nc` calls double-counts it. |
| `dnsq_scf`, `sym_dns` moved to `LR_Modules`; `sym_dns` uses `upert`/`upert_mq`; `syme_dns` removed; `dnsscf_all_modes` handled by `dnsq_store` | nc helpers must move from `ldaU_ph` to `ldaU_lr`; `sym_dns_nc` rewritten on `upert`; `syme_dns_nc` dropped (E-field handled through `upert=s`). |
| `dnsq_scf` lost `lflag`; allocates `dpsi(npwx,nbnd)` | Need a new phonon/E-field branch signal; `dpsi` must be `npwx*npol`. |
| `dfpt_kernel` CG threshold `min(0.1*sqrt(dr2), thresh)` keeps memory | qe-pert threshold overrides need a separate adaptive variable. |
| `intra_pool_comm` -> `intra_bgrp_comm` in Hubbard routines | `commutator_Vhubx_psi` merge does not compile as is; `ph_acc` kernel must switch (silently wrong with `-nb>1` otherwise). |
| `phq_readin.f90:548` forbids magnetic + epsil/zue; `:779` forbids `magnetic_sym` + U | Feature guards must be placed around these; magnetic E-field needs new work (decision D2). |
| magnetic E-field: `dfpt_kernel` uses `nsolv=2` but `dvpsi_e` writes only one branch | qe-pert's single-record-plus-adjoint E-field scheme has no 7.6 home yet. |
| `dnsbare_pattern` write (present in pristine 7.3) removed | Keep qe-pert's write: qe2pert reads it. |
| Orbital-resolved U (`Hubbard_Um_nc`, `orbital_resolved`) added | nc routines use `effU`; add a guard. |
| `average_pp` now errors for FR PPs with `lspinorb=.false.` | SOC-off control runs (qe-pert B9/Gate 5 style) need scalar-relativistic PPs on 7.6. |
| HP has `revert_mag_u`/`calc_vh_u` (ground Hubbard potential for -B) | Candidate reuse/cross-check for qe-pert `time_reverse_ground_hubbard_nc`; equivalence unverified. Do not flip `v%ns_nc` inside `lr_apply_time_reversal` (would double-flip HP). |

## 4. Decisions

| ID | Question | Decision |
| --- | --- | --- |
| D1 | How to select the SOC+U path in 7.6 | **Decided (user, 2026-10-23):** `ph.x` input `lhubbard_soc` in `&INPUTPH`, default `.false.`. See 4.1. |
| D2 | Magnetic static epsil/BEC with SOC+U | **Decided (user, 2026-10-23):** in scope. Already working in qe-pert, but must be re-implemented for 7.6 (see 4.2). Ported after the phonon SCF path (phase 3b). |
| D3 | Control twins (FIXED_CG, LEGACY_MIX, ORIGINAL_CONV, JOINT_HXC_ONLY, NOHUB) and debug paths | Defer until phase 3 passes; port only on request. |
| D4 | `ph_acc` collinear NCPP kernel | Port with `intra_bgrp_comm` as `ph_acc.x` (`dynmat_hub_bare.f90` stays a plain PH file, so `-D__PH_ACC` still works). Making it the default in `ph.x` is a later option. |
| D5 | `dnsbare_pattern` output | Keep (qe2pert). |
| D6 | Noncollinear GGA near abs(m) -> 0 (singular, sign-indefinite transverse kernel of the pw.x local-frame functional) | **Decided (user, 2026-10-26):** leave the functional as is. Production DFPT runs keep symmetry, where the runs converge; nosym noncollinear-GGA runs (static E-field) are not a target. The two consistency fixes (`post_xml_init` lsign/ux, `dgradcorr` FFT gradient) stay. |

### 4.1 `lhubbard_soc` (replaces the qe-pert `__SOC_U_DFPT_*_VALIDATE` executables)

- Declared in an LR module (`ldaU_lr` in `LR_Modules/lrcom.f90`), because
  `dfpt_kernel`, `dnsq_scf` and `sym_dns` live in `LR_Modules` and cannot use
  `control_ph`. Read in `&INPUTPH` by `phq_readin`, broadcast with the other
  input variables, documented in `PHonon/Doc/INPUT_PH.def`.
- `.false.` (default): ph.x behaves exactly as pristine 7.6, including the
  noncollinear+U guards (`:779` for `magnetic_sym`), except that qe-pert's
  stricter production guard (reject any noncolin+U) is kept.
- `.true.`: enables the qe-pert operator path (former
  `__SOC_U_DFPT_Q_OPERATOR_VALIDATE`) for noncollinear magnetic or
  non-magnetic spinor states, with or without SOC. Still rejected with a clear
  message: USPP/PAW, non-`atomic` projectors, `lda_plus_u_kind/=0`, background
  manifolds, `orbital_resolved`, Hubbard J0/alpha/beta, smearing/tetrahedra,
  `elph`, Raman, `elop`, `fpol`, `read_dns_bare`. `epsil`/`zeu` stay rejected
  until phase 3b is done.
- Setting it with `noncolin=.false.` is an error, not a silent no-op.
  `ph_acc.x` (collinear) is independent of it.
- The qe-pert staged modes (`BARE`, `Q_BARE`, `FULL`) were narrower guards for
  validation; with `only_init=.true.` the same bare-only runs remain possible.
- The control twins (fixed CG, legacy mix, original convergence, Hxc-only exit,
  no Hubbard feedback) become further runtime options only on request (D3).

### 4.2 Magnetic E-field in 7.6

qe-pert (7.3) solves the static E-field with one Sternheimer branch and adds
the time-reversed partner inside `dnsq_scf` as an adjoint term
(`lflag=.FALSE.`); pristine 7.3 had no magnetic E-field guard. 7.6 instead
(a) rejects magnetic + epsil/zue for all runs (`phq_readin.f90:548`), (b) sets
`nsolv=2` in `dfpt_kernel` for any `noncolin.AND.domag`
(`dfpt_kernels.f90:196`), and (c) `dvpsi_e` writes only the first-branch
record. Two ways to port, to be chosen when phase 3a is done:
1. Keep the qe-pert scheme: `nsolv=1` for the static E-field when
   `lhubbard_soc`, and the adjoint term in `dnsq_scf` (needs the phonon/E-field
   branch signal that replaced `lflag`).
2. Follow 7.6: write time-reversed `dvpsi_e` records and solve both branches.
Option 1 reproduces the validated qe-pert results directly; option 2 is the
7.6-native form. Either way, compare epsilon and Born charges with qe-pert.

## 5. Phases

Each step: name the qe-pert hunks and the 7.6 owner, edit, build
(`make -j4 pw ph`), run the step's check, log in PORT_STATUS.

| Phase | Content | Exit check |
| --- | --- | --- |
| 0 (done) | Baseline tree, make.inc, baseline build, cross-version fixture, dry-run tool | Pristine 7.6 `pw.x`/`ph.x` built; hBN reference outputs stored |
| 1 Infrastructure | `lrcom` kernels (+ helpers moved from `phcom`), `mix_pot` (`hubbard_response_mixing`, `mix_potential` wrapper), `ph_restart` checkpoint, `allocate_phq`, `close_phq`, `lr_orthoUwfc`, `phq_setup` (`comp_dspinldau`), `print_clock_ph`, `dynmatrix`/`q2qstar_ph`/`rotate_and_add_dyn` (t_rev), `commutator_Vhubx_psi` (comm fix), `force_hub` diagnostics | Build. Collinear regression identical to the pristine build: hBN fixtures and test-suite `ph_U_*`, `hp_*`, `ph_base`-type cases. `mix_potential` changes every caller, so HP/KCW cases must also match. |
| 2 Bare response | `lhubbard_soc` (4.1) and its guards; `dnsq_bare` nc; `dvqhub_barepsi_nc` + `dvqpsi_us` dispatch (no double count) + `zstar_eu`; `sym_dns_nc` on `upert`; `dynmat_hub_bare_nc`; `ph_acc` | Collinear unchanged; `ph_acc.x` vs 7.6 `ph.x` exact as in qe-pert; SOC+U `only_init` on a small fixture: `dnsbare`, `dyn_hub_bare` vs qe-pert within cross-version noise; Hermiticity/ASR diagnostics. |
| 3a SCF response (phonons) | `dnsq_scf` nc (two branches, spinor `dpsi`, branch signal); `adddvhubscf_nc` + `response_kernels` dispatch (+ errors in freq/twochem copies); `dfpt_kernel` hooks: ground Hubbard time reversal for isolv=2, `mix_iter`, joint mixing, threshold handling, restart rebuild; `phqscf` export/restore; `dynmat_hub_scf_nc`; final phonon guards | Small SOC+U Gamma and one general q: dynamical matrix and frequencies vs qe-pert `ph_soc_u_q_operator_validate.x`; restart roundtrip; MPI layout (pools) independence. |
| 3b Static E-field | 4.2: `solve_e`/`dfpt_kernel` branch handling, `dnsq_scf` adjoint or reversed `dvpsi_e` records, `commutator_Vhubx_psi`, `zstar_eu`, lift `:548` only for `lhubbard_soc` | epsilon and Born effective charges (and E-U charges) vs qe-pert on the same small fixture; ASR/charge neutrality residual reported. |
| 4 Optional | D3 twins, diagnostics (`c_bands_ph`, `mix_rho` residual, `phq_init` covariance) | Per item. |
| 5 Separate | qe2pert on 7.6 | Needs a 7.6-compatible perturbo first. |

## 6. Fixtures and tolerances

- Collinear: `../qe-7.6-pert-buildlogs/xver_hbn_sr/` (hBN, PseudoDojo SR NC,
  U(N-2p)=2 eV, Gamma and q=(0.25,0,0)). Pristine 7.3 vs 7.6 already differ
  by 2.6e-6 Ry and at most 0.034 cm^-1: cross-version comparisons cannot be
  exact. Within 7.6, the ported collinear path must reproduce the pristine
  build exactly.
- SOC+U: a small magnetic SOC insulator fixture (FR NC PPs, reduced cutoff
  and k mesh) still has to be prepared and run on both qe-pert and this tree.
  Candidate inputs: NiO under `~/work/NiO` (qe-pert runs) and
  `~/work/NiO/7.6/pw` (used by the separate 7.6 effort). Material-size runs
  need user authorization.
- Reusable analysis scripts in `~/bin/qe-pert/dev-tools/`
  (`check_soc_u_*.py`, `compare_soc_u_fd_*.py`); native fixtures
  (`test_soc_u_hubbard_restart.f90`, `test_soc_u_joint_mixing*.f90`) need
  adaptation to the 7.6 module layout.

## 7. Build notes

- Make build from the top directory; local `make.inc` = 7.6 configure output
  with HDF5 and libxc lines taken from qe-pert (7.6 configure misdetects
  `lib64` HDF5 and fails the libxc link test with MKL). Files:
  `../qe-7.6-pert-buildlogs/`.
- New source files must be registered in the Makefiles and the CMakeLists.
- Build products not in the upstream `.gitignore` are excluded locally in
  `.git/info/exclude`.
