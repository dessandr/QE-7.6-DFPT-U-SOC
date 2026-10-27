# QE 7.6 port of the qe-pert DFPT+SOC+U implementation

This tree is pristine QE 7.6 (commit `5e4b172`) into which the noncollinear/SOC
fixed-Dudarev DFPT+U of `~/bin/qe-pert` (QE 7.3) is being ported. The porting
plan, file map and open decisions are in
[docs/soc_u_port/PORT_PLAN.md](docs/soc_u_port/PORT_PLAN.md); record progress in
[docs/soc_u_port/PORT_STATUS.md](docs/soc_u_port/PORT_STATUS.md).

## Sources of truth

| Source | Role |
| --- | --- |
| This tree's QE 7.6 code | Interfaces, call flow, array shapes. 7.6 moved the PH SCF loop into `LR_Modules/dfpt_kernels.f90` and moved `dnsq_scf`/`sym_dns` to `LR_Modules`. |
| `~/bin/qe-pert`, branch `ph-acc`, diff from `347e370` (pristine 7.3) | The behavior to port. Never copy a whole qe-pert file over a 7.6 file; port hunks onto the 7.6 owner. |
| `~/bin/qe-pert/docs/soc_u_dfpt/` | Derivations, conventions, validation history and known issues. Its results are 7.3 results, not 7.6 evidence. |
| `~/bin/q-e-qe-7.6` | A separate, independent 7.6 implementation. Do not modify it; it may be in use by another session. |

In qe-pert the ordinary `ph.x` keeps the noncollinear+U guard. The working
SOC+U path is the `__SOC_U_DFPT_Q_OPERATOR_VALIDATE` build
(`ph_soc_u_q_operator_validate.x`); the other `__SOC_U_*_VALIDATE` macros are
diagnostic variants of it. In this tree that path is selected at run time by
the `&INPUTPH` switch `lhubbard_soc` (PORT_PLAN 4.1), not by separate builds.

## Physics scope (inherited, unchanged by the port)

Fixed U_eff with dU/du=0, Dudarev onsite U, atomic projectors, insulating
noncollinear/SOC states, norm-conserving pseudopotentials, one Hubbard manifold
per atom. Keep all four complex spin blocks, delta N = delta N_bare + delta
N_scf, delta V_U = -U_eff delta N, both magnetic Sternheimer branches, and the
existing guards. Preserve 7.6 collinear+U and U=0 behavior. The port moves
known qe-pert open issues along with the code (see the known-issues handoff);
matching qe-pert is a porting check, not a physics validation.

## Working rules

- Port in the order of PORT_PLAN. Each step: state which qe-pert hunks and
  which 7.6 owner, make the edit, build, and run that step's check.
- `dev-tools/soc_u_port_merge_dryrun.sh OUT_DIR` regenerates the 3-way
  base(7.3)/qe-pert/7.6 comparison without touching the tree.
- Build with Make from the top directory (`make -j4 pw ph`); `make.inc` is a
  local, untracked file (provenance in `../qe-7.6-pert-buildlogs/`). Register
  any new source file in both the Makefiles and the CMakeLists.
- Locally (pcs-gpu, the shared login node) keep total CPU use at or below
  five cores, counting the user's other jobs; numerical analysis
  single-threaded. Builds and short checks run locally. Validation runs that
  take more than a few minutes (NiO fixture phonons/E-field, regression sets)
  go to a PBS compute node with `dev-tools/soc_u_port_qsub.sh`; submitting
  these validation jobs is authorized. Larger material runs need explicit
  authorization.
- A/B checks compare this tree's executables with the matching qe-pert
  executables on identical inputs and report the observable, the max
  difference and the tolerance.
- Commit only when asked. Respond to the user concisely in Korean.
