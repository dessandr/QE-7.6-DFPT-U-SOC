!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  adapted from PW/src/atomic_wfc.f90
!
!  QE 7.6: the atomic wavefunctions are built by atomic_wfc_acc
!  (Modules/atomic_wfc_mod.f90), as in PW/src/atomic_wfc.f90. These wrappers
!  keep the Perturbo interface (explicit xk, npw, igk_k instead of the k-point
!  index) and pass the same options as the PW wrappers, so the Hubbard
!  projectors are those of pw.x and ph.x of this tree.
!-----------------------------------------------------------------------
SUBROUTINE pw_atomic_wfc(xk, npw, igk_k, wfcatom )
  !-----------------------------------------------------------------------
  !! This routine computes the superposition of atomic wavefunctions
  !! for k-point "xk" - output in "wfcatom".
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, tau, nsp, ityp
  USE basis,            ONLY : natomwfc
  USE wvfct,            ONLY : npwx
  USE noncollin_module, ONLY : noncolin, domag, npol, &
                               angle1, angle2, starting_spin_angle
  USE atomic_wfc_mod,   ONLY : atomic_wfc_acc
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: xk(3)
  !! xk in cartesian coordinate
  INTEGER, INTENT(IN) :: npw
  INTEGER, INTENT(IN) :: igk_k(npw)
  COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx,npol,natomwfc)
  !! Superposition of atomic wavefunctions
  !
  LOGICAL :: starting_spin_angle_, updown, spinorbit
  !
  ! as atomic_wfc in PW/src/atomic_wfc.f90
  starting_spin_angle_ = starting_spin_angle .OR. .NOT.domag
  updown = .FALSE.
  spinorbit = .FALSE.
  CALL atomic_wfc_acc( xk, npw, igk_k, nat, nsp, ityp, tau, &
       noncolin, spinorbit, updown, angle1, angle2, starting_spin_angle_, &
       npwx, npol, natomwfc, wfcatom )
  !
END SUBROUTINE pw_atomic_wfc
!
!-----------------------------------------------------------------------
SUBROUTINE pw_atomic_wfc_nc_updown (xk, npw, igk_k, wfcatom)
  !-----------------------------------------------------------------------
  !! For noncollinear case: builds up the superposition (for a k-point "xk")
  !! of pure spin up or spin down atomic wavefunctions; with fully
  !! relativistic pseudopotentials the radial functions are j-averaged.
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, tau, nsp, ityp
  USE basis,            ONLY : natomwfc
  USE wvfct,            ONLY : npwx
  USE noncollin_module, ONLY : noncolin, npol, angle1, angle2
  USE atomic_wfc_mod,   ONLY : atomic_wfc_acc
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: xk(3)
  !! xk in cartesian coordinate
  INTEGER, INTENT(IN) :: npw
  INTEGER, INTENT(IN) :: igk_k(npw)
  COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx,npol,natomwfc)
  !
  LOGICAL :: starting_spin_angle_, updown, spinorbit
  !
  ! as atomic_wfc_nc_updown in PW/src/atomic_wfc.f90
  starting_spin_angle_ = .FALSE.
  updown = .TRUE.
  spinorbit = .FALSE.
  CALL atomic_wfc_acc( xk, npw, igk_k, nat, nsp, ityp, tau, &
       noncolin, spinorbit, updown, angle1, angle2, starting_spin_angle_, &
       npwx, npol, natomwfc, wfcatom )
  !
END SUBROUTINE pw_atomic_wfc_nc_updown
