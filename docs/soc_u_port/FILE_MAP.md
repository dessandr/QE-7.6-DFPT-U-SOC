# qe-pert -> QE 7.6 file map

Read-only source analysis of 2026-10-23 (qe-pert `ph-acc` vs pristine 7.3
`347e370` vs this tree's `5e4b172`). Line numbers: `qp:` = qe-pert, otherwise
7.6. Items marked (v) were re-checked directly in the 7.6 source; the rest
are analysis to confirm at the time of the edit.

Actions: **APPLY** textual merge is correct; **ADAPT** re-implement at the 7.6
owner; **DEFER** control/diagnostic, port later on request; **DROP**
superseded in 7.6. Classes: **core** (feature), **OPV** (feature-path oracle
under `__SOC_U_DFPT_Q_OPERATOR_VALIDATE`, no effect on results), **CTRL**
(control twin), **DBG** (debug/iverbosity diagnostic).

## A. Hubbard occupation and operators

| qe-pert file | What | Class | 7.6 target | Action | Notes |
| --- | --- | --- | --- | --- | --- |
| `LR_Modules/lrcom.f90` | `pack_hubbard_blocks_nc`, `build_dvhub_nc`, time-reverse helpers, `apply_hubbard_matrix_nc` | core | `ldaU_lr` CONTAINS (before `END MODULE`, ~240) | APPLY | 7.6 adds `lr_has_dnsorth`, `lr_dnsorth`, `proj1/proj2`, `vh_u_save` declarations: keep them before CONTAINS. |
| `PHonon/PH/phcom.f90` | `accumulate_ns_nc_transition`, `ns_nc_hermiticity_residual` | core | move to `ldaU_lr` | ADAPT | Consumer `dnsq_scf` is now in LR and cannot use `ldaU_ph`. `nc_mag_aux` is `lr_nc_mag` in 7.6. |
| `LR_Modules/adddvhubscf.f90` | `adddvhubscf_nc` (delta V_U = -U_eff delta N_scf, both branches) | core+OPV | append after 7.6 `calc_vh_uv` | APPLY | Conflict = two appends, keep both. |
| `LR_Modules/lr_orthoUwfc.f90` | save `iunhub`/`iunhub_noS` at `ikk` | core | :138, :148 | APPLY | Clean. |
| `PHonon/PH/dnsq_bare.f90` | spinor buffers, nc ZGEMM path, q-partner (ikmks) branch, no triangle fill / 1/2 for nc | core | same file | ADAPT | Conflicts 1-3; rebuild the collinear ELSE from 7.6 text. `proj1/proj2` from `ldaU_lr`. |
| ″ | `dndtau_k_nc` cross-check, branch dumps | DBG | same | APPLY | GPU build would need `wfcU` device updates. |
| ″ | `dnsbare_pattern` write | core (qe2pert) | conflict 4 (~:300) | KEEP | Present in pristine 7.3, removed in 7.6 (v). |
| `PHonon/PH/dnsq_scf.f90` | nc branch: two Sternheimer records, `accumulate_ns_nc_transition`, pool sum, Hermiticity print | core | `LR_Modules/dnsq_scf.f90` `(npe,lmetq0)` | ADAPT | `dpsi` must be `npwx*npol` (7.6: `npwx,nbnd`) (v). `lflag` gone: new phonon/E-field branch signal needed. |
| ″ | `dnsscf_all_modes` store, `dnsorth` | — | `dnsq_store`, `lr_dnsorth` | DROP | Already in 7.6. Post-solve refresh for joint mixing goes into `dfpt_kernel` (section B). |
| `PHonon/PH/sym_dns.f90` | `sym_dns_nc` (magnetic little group, `d_spin_ldau`, t_rev conj), `rotate_dns_nc_symmetry` | core | `LR_Modules/sym_dns.f90` | ADAPT | Replace `t(:,:,isym,irr)` with `upert(:,:,isym)`, `tmq` with `upert_mq`; no `modes`/`irr` in LR. |
| `PHonon/PH/syme_dns.f90` | `syme_dns_nc` | core | — | DROP | `phescf` sets `upert=s`; generic `sym_dns_nc` covers it. Remove from Makefile. |
| `PHonon/PH/dvqhub_barepsi_us.f90` | `dvqhub_barepsi_nc` (normal/time-reversed bare term) | core+OPV | append after :390 | APPLY | Callers: see next row. |
| (7.6) `PHonon/PH/dvqpsi_us.f90` | — | core | :191 (v) | ADAPT | Add noncolin dispatch with time-reversed/oracle args; remove qe-pert explicit calls in `solve_linter`/`zstar_eu` to avoid double counting. |
| `PHonon/PH/zstar_eu.f90` | nc bare term | core | — | ADAPT | Take 7.6 side; term now via `dvqpsi_us`. |
| `PHonon/PH/dynmat_hub_bare.f90` | `dynmat_hub_bare_nc` | core | same | APPLY | |
| ″ | `__PH_ACC` onsite kernel | ph_acc | same | ADAPT | 7.6 has no `#ifdef __MPI` there; use `intra_bgrp_comm` (v). |
| `PHonon/PH/dynmat_hub_scf.f90` | `dynmat_hub_scf_nc` | core | same | ADAPT | `dnsq_scf(nper,lmetq0)`; `rec_code_read` now in `control_lr`. |
| `PHonon/PH/phq_setup.f90` | `comp_dspinldau` after `d_matrix` | core | :445 | APPLY | USE from `lr_nc_mag`. |
| `PHonon/PH/allocate_phq.f90`, `close_phq.f90`, `print_clock_ph.f90` | spinor sizes, close `iunhub`, `mix_hubbard` clock | core | same | APPLY | |
| `PHonon/PH/phq_readin.f90` | guard relaxations: OPERATOR (magnetic phonons + static epsil/zeu), Q_BARE, FULL, BARE | core | same | ADAPT | Becomes the runtime switch (D1). 7.6 guards :548 (magnetic E-field) and :779 (`magnetic_sym`+U) (v) must be placed accordingly. Keep qe-pert's stricter production guard (7.6 lets nonmagnetic noncolin+U through). Add an `orbital_resolved` guard. |
| `PHonon/PH/phq_init.f90` | `validate_ground_hamiltonian_covariance_nc` | DBG (Q_BARE) | same | DEFER | |
| `PHonon/PH/{dynmatrix,q2qstar_ph,rotate_and_add_dyn}.f90` | antiunitary q-star rotation | core | same | APPLY | 7.6 code there identical to 7.3. |
| `PW/src/commutator_Vhubx_psi.f90` | noncollinear [V_U, r] for `dvpsi_e` | core | same | APPLY+fix | Clean merge does not compile: use `intra_bgrp_comm` (v). |
| `PW/src/force_hub.f90` | SOC_U oracle dump | DBG | :193, :335 | APPLY | USE list: 7.6 removed `nsg`. |

## B. SCF loop, mixing, restart

| qe-pert file | What | Class | 7.6 target | Action | Notes |
| --- | --- | --- | --- | --- | --- |
| `PHonon/PH/solve_linter.f90` qp:36-97 | USEs | core | `dfpt_kernels.f90::dfpt_kernel` USE block | ADAPT | Take 7.6 side of all 13 conflicts, then re-implement. |
| ″ qp:177-299 | FIRST_ITER/ALL_ITER threshold override (env, broadcast) | core (FIXED_CG nested = CTRL) | `dfpt_kernel` before loop | ADAPT | Runtime flag; separate variable because 7.6 threshold has memory (v). |
| ″ qp:329-356 | joint-mixing setup (`dnsscf_input`, `mixing_ueff`, `mixing_ldim`) | core | after `apply_dpot_allocate` | ADAPT | Needs phonon vs E-field flag; locals or a `dfpt_type` field. |
| ″ qp:445-461 | `dvqhub_barepsi_nc` for isolv 1/2 | core | `dvqpsi_us` dispatch | ADAPT | See section A. |
| ″ qp:477-511 | `mix_iter`, restart rebuild of `dnsscf_input`, Broyden cold start, use `dnsscf_input` | core | `dfpt_kernel` before `do kter`, :272 | ADAPT | 7.6 condition `.NOT.first_iter`. |
| ″ qp:520-608 | `time_reverse_ground_hubbard_nc` around the isolv=2 solve | core | around the isolv=2 `sternheimer_kernel` call | ADAPT | Not inside `lr_apply_time_reversal` (HP double flip). Host flip suffices (`vhpsi_nc_acc` copies `v%ns_nc`). |
| ″ qp:538-575 | CG threshold override + log | core | :262-267 | ADAPT | |
| ″ qp:717-758 | joint Hxc+Hubbard Broyden (`dnsq_scf` refresh, `mix_hubbard_response`), `mix_iter` in mixer calls | core (JOINT_HXC_ONLY = CTRL) | :395-406 | ADAPT | Keep `dfpt_dvscfp_to_dvscfs` after it; leave raw `dnsscf` for `dnsq_store`. |
| ″ qp:767-816 | final-refresh Hubbard gate (non-joint) | CTRL | after `check_all_convt` | DEFER | |
| ″ qp:951-1004 | `SUBROUTINE time_reverse_ground_hubbard_nc` | core | LR_Modules (`ldaU_lr`/`lr_nc_mag`) | ADAPT | Compare with HP `revert_mag_u`/`calc_vh_u` before choosing. |
| `PHonon/PH/solve_e.f90` | joint mixing with 3 fields, restart rebuild, logs | core (LEGACY_MIX = CTRL) | `dfpt_kernel` | ADAPT / DEFER | Magnetic E-field blocked in 7.6 (D2). |
| `PHonon/PH/phqscf.f90` | invalidate/restore per-irrep `data_hub`, gated flat export, `prepare_dnsscf_export` | core | after :69, :130, :171-176 | APPLY | USE conflict: `convt`, `rec_code` now from `control_lr`. |
| `PHonon/PH/ph_restart.f90` | `HUBBARD_RESPONSE` checkpoint write/read | core | same | APPLY | Clean; `dnsq_store` fills `dnsscf_all_modes` before `write_rec`. |
| `LR_Modules/mix_pot.f90` | `hubbard_response_mixing` module; `mix_potential` -> wrapper over `mix_potential_core` | core | same | APPLY | Identical in 7.3/7.6. Changes behavior for all callers (HP, KCW, `solve_e2`, `solve_e_fpol`, `dfpt_kernel`): regression-test them. |
| `LR_Modules/response_kernels.f90` | `adddvhubscf_nc` dispatch | core | :232 | APPLY | Also error out in the freq copy (:530) and `lr_two_chem.f90:421` for noncolin+U. |
| `LR_Modules/c_bands_ph.f90` | `SOC_U_WFCU_RECORD_CHECK` | DBG | same | DEFER | GPU false-pass risk. |
| `PW/src/mix_rho.f90` | residual decomposition | DBG | `MODULE mix` procedure | DEFER | |
| `PW/src/Makefile`, `PHonon/PH/Makefile` | debug/validation executables | — | — | per D1/D3 | Per-executable `-D` cannot reach `dfpt_kernels`. |

## C. Fixes made in the 7.6 port (not in qe-pert)

| 7.6 file | What | Notes |
| --- | --- | --- |
| `PW/src/read_file_new.f90` (`post_xml_init`) | set ux/lsign (`compute_ux`) before `v_of_rho` for noncolin + domag + GGA | Upstream bug: ph.x kept a lsign=.false. B_xc in `vrs`. See PORT_STATUS 2026-10-25. |
| `LR_Modules/dgradcorr.f90` | noncollinear branch: FFT gradient of delta(rho +- s abs(m))/2 | Upstream bug: `gc_lr%gmag` never computed; chain rule not the derivative of the discrete v_xc. |
| `dev-tools/xckfd/` | `xckfd.x`: FD check of the noncollinear xc kernel | Test tool; `make -C dev-tools/xckfd` after `make ph`. |
