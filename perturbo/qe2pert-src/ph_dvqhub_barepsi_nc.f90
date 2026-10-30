!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! Copyright (C) 2026 Perturbo developers
! This file is distributed under the terms of the GNU General Public License.
!
! Noncollinear, norm-conserving qe2pert adapter for the bare derivative of
! the Hubbard potential.
!-----------------------------------------------------------------------------
SUBROUTINE ph_dvqhub_barepsi_nc(ikk, xk, xkq, npwq, igk_kq, &
     wfcatomk, swfcatomk, wfcatomkpq, swfcatomkpq, dnsbare, &
     icart_disp, na_disp, dvqhbar)
  !---------------------------------------------------------------------------
  !! Apply the three fixed-Dudarev bare terms
  !!
  !!   dP(k+q) V0 P(k)^dagger
  !! + P(k+q) V0 dP(k)^dagger
  !! + P(k+q) [-Ueff dN_bare(q)] P(k)^dagger
  !!
  !! to the qe2pert Bloch/Wannier-gauge states at k.
  !
  USE kinds,              ONLY : DP
  USE ions_base,          ONLY : nat, ityp
  USE ldaU,               ONLY : Hubbard_lmax, Hubbard_l, offsetU, &
                                 is_hubbard, nwfcU
  USE ldaU_lr,            ONLY : pack_hubbard_blocks_nc, &
                                 build_dvhub_nc, apply_hubbard_matrix_nc
  USE noncollin_module,   ONLY : npol
  USE scf,                ONLY : v
  USE wvfct,              ONLY : npwx
  USE dvhub_mod,          ONLY : ldim, effU
  USE electronic_data,    ONLY : ngk_tot, igk_k_tot, nb_sub, evc_sub
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ikk, npwq, icart_disp, na_disp
  INTEGER, INTENT(IN) :: igk_kq(npwx)
  REAL(DP), INTENT(IN) :: xk(3), xkq(3)
  COMPLEX(DP), INTENT(IN) :: wfcatomk(npwx*npol,nwfcU)
  COMPLEX(DP), INTENT(IN) :: swfcatomk(npwx*npol,nwfcU)
  COMPLEX(DP), INTENT(IN) :: wfcatomkpq(npwx*npol,nwfcU)
  COMPLEX(DP), INTENT(IN) :: swfcatomkpq(npwx*npol,nwfcU)
  COMPLEX(DP), INTENT(IN) :: dnsbare(ldim,ldim,npol*npol,nat)
  COMPLEX(DP), INTENT(OUT) :: dvqhbar(npwx*npol,nb_sub)
  !
  INTEGER :: na, nt, projector_spin, row_spin, m, iwfc, row_first
  INTEGER :: ldim_atom, max_matrix_size
  COMPLEX(DP), ALLOCATABLE :: dprojector_k(:,:), dprojector_kpq(:,:)
  COMPLEX(DP), ALLOCATABLE :: dcolumn(:), projection(:,:), dprojection(:,:)
  COMPLEX(DP), ALLOCATABLE :: ground_matrix(:,:), bare_matrix(:,:)
  !
  IF (npol /= 2) CALL errore('ph_dvqhub_barepsi_nc', &
                             'noncollinear Hubbard response requires npol=2', 1)
  IF (.NOT.ALLOCATED(v%ns_nc)) CALL errore('ph_dvqhub_barepsi_nc', &
                                           'PW Hubbard potential is unavailable', 1)
  IF (na_disp < 1 .OR. na_disp > nat .OR. &
      icart_disp < 1 .OR. icart_disp > 3) &
     CALL errore('ph_dvqhub_barepsi_nc', 'invalid atomic displacement', 1)
  !
  max_matrix_size = (2*Hubbard_lmax + 1)*npol
  ALLOCATE(dprojector_k(npwx*npol,nwfcU))
  ALLOCATE(dprojector_kpq(npwx*npol,nwfcU))
  ALLOCATE(dcolumn(npwx*npol))
  ALLOCATE(projection(nwfcU,nb_sub), dprojection(nwfcU,nb_sub))
  ALLOCATE(ground_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(bare_matrix(max_matrix_size,max_matrix_size))
  !
  dprojector_k = (0.0_DP,0.0_DP)
  dprojector_kpq = (0.0_DP,0.0_DP)
  nt = ityp(na_disp)
  IF (is_hubbard(nt)) THEN
     ldim_atom = 2*Hubbard_l(nt) + 1
     DO projector_spin = 1, npol
        DO m = 1, ldim_atom
           iwfc = offsetU(na_disp) + m + &
              ldim_atom*(projector_spin-1)
           dcolumn = (0.0_DP,0.0_DP)
           DO row_spin = 1, npol
              row_first = 1 + npwx*(row_spin-1)
              CALL ph_dwfc(ngk_tot(ikk), igk_k_tot(1,ikk), xk, &
                           icart_disp, wfcatomk(row_first,iwfc), &
                           dcolumn(row_first))
           ENDDO
           dprojector_k(:,iwfc) = dcolumn
           dcolumn = (0.0_DP,0.0_DP)
           DO row_spin = 1, npol
              row_first = 1 + npwx*(row_spin-1)
              CALL ph_dwfc(npwq, igk_kq, xkq, icart_disp, &
                           wfcatomkpq(row_first,iwfc), &
                           dcolumn(row_first))
           ENDDO
           dprojector_kpq(:,iwfc) = dcolumn
        ENDDO
     ENDDO
  ENDIF
  !
  CALL ZGEMM('C', 'N', nwfcU, nb_sub, npwx*npol, &
             (1.0_DP,0.0_DP), swfcatomk, npwx*npol, &
             evc_sub(1,1,ikk), npwx*npol, (0.0_DP,0.0_DP), &
             projection, nwfcU)
  CALL ZGEMM('C', 'N', nwfcU, nb_sub, npwx*npol, &
             (1.0_DP,0.0_DP), dprojector_k, npwx*npol, &
             evc_sub(1,1,ikk), npwx*npol, (0.0_DP,0.0_DP), &
             dprojection, nwfcU)
  !
  dvqhbar = (0.0_DP,0.0_DP)
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        ldim_atom = 2*Hubbard_l(nt) + 1
        CALL pack_hubbard_blocks_nc(ldim_atom, npol, &
                                    v%ns_nc(:,:,:,na), ground_matrix)
        CALL build_dvhub_nc(ldim_atom, npol, effU(nt), &
                            dnsbare(:,:,:,na), bare_matrix)
        CALL apply_hubbard_matrix_nc(npwx*npol, nb_sub, nwfcU, &
             ldim_atom, npol, offsetU(na), dprojector_kpq, projection, &
             ground_matrix, dvqhbar)
        CALL apply_hubbard_matrix_nc(npwx*npol, nb_sub, nwfcU, &
             ldim_atom, npol, offsetU(na), swfcatomkpq, dprojection, &
             ground_matrix, dvqhbar)
        CALL apply_hubbard_matrix_nc(npwx*npol, nb_sub, nwfcU, &
             ldim_atom, npol, offsetU(na), swfcatomkpq, projection, &
             bare_matrix, dvqhbar)
     ENDIF
  ENDDO
  !
  DEALLOCATE(dprojector_k, dprojector_kpq, dcolumn)
  DEALLOCATE(projection, dprojection, ground_matrix, bare_matrix)
  !
  RETURN
END SUBROUTINE ph_dvqhub_barepsi_nc
