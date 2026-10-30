!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! Copyright (C) 2026 Perturbo developers
! This file is distributed under the terms of the GNU General Public License.
!
! Noncollinear, norm-conserving qe2pert adapter for the self-consistent
! derivative of the Hubbard potential.
!-----------------------------------------------------------------------------
SUBROUTINE ph_adddvhubscf_nc(ikk, swfcatomk, swfcatomkpq, dnsscf, dvhubscf)
  !---------------------------------------------------------------------------
  !! Apply
  !!
  !!   P(k+q) [-Ueff dN_scf(q)] P(k)^dagger |psi(k)>
  !!
  !! in the full orbital-spin representation used by PW.  PH response data
  !! are supplied by qe2pert instead of being taken from ldaU_lr globals.
  !
  USE kinds,              ONLY : DP
  USE ions_base,          ONLY : nat, ityp
  USE ldaU,               ONLY : Hubbard_lmax, Hubbard_l, offsetU, &
                                 is_hubbard, nwfcU
  USE ldaU_lr,            ONLY : build_dvhub_nc, apply_hubbard_matrix_nc
  USE noncollin_module,   ONLY : npol
  USE wvfct,              ONLY : npwx
  USE dvhub_mod,          ONLY : ldim, effU
  USE electronic_data,    ONLY : nb_sub, evc_sub
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ikk
  COMPLEX(DP), INTENT(IN) :: swfcatomk(npwx*npol,nwfcU)
  COMPLEX(DP), INTENT(IN) :: swfcatomkpq(npwx*npol,nwfcU)
  COMPLEX(DP), INTENT(IN) :: dnsscf(ldim,ldim,npol*npol,nat)
  COMPLEX(DP), INTENT(OUT) :: dvhubscf(npwx*npol,nb_sub)
  !
  INTEGER :: na, nt, ldim_atom, max_matrix_size
  COMPLEX(DP), ALLOCATABLE :: projection(:,:), dvhub(:,:)
  !
  IF (npol /= 2) CALL errore('ph_adddvhubscf_nc', &
                             'noncollinear Hubbard response requires npol=2', 1)
  !
  max_matrix_size = (2*Hubbard_lmax + 1)*npol
  ALLOCATE(projection(nwfcU,nb_sub))
  ALLOCATE(dvhub(max_matrix_size,max_matrix_size))
  !
  CALL ZGEMM('C', 'N', nwfcU, nb_sub, npwx*npol, &
             (1.0_DP,0.0_DP), swfcatomk, npwx*npol, &
             evc_sub(1,1,ikk), npwx*npol, (0.0_DP,0.0_DP), &
             projection, nwfcU)
  !
  dvhubscf = (0.0_DP,0.0_DP)
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        ldim_atom = 2*Hubbard_l(nt) + 1
        CALL build_dvhub_nc(ldim_atom, npol, effU(nt), &
                            dnsscf(:,:,:,na), dvhub)
        CALL apply_hubbard_matrix_nc(npwx*npol, nb_sub, nwfcU, &
             ldim_atom, npol, offsetU(na), swfcatomkpq, projection, &
             dvhub, dvhubscf)
     ENDIF
  ENDDO
  !
  DEALLOCATE(projection, dvhub)
  !
  RETURN
END SUBROUTINE ph_adddvhubscf_nc
