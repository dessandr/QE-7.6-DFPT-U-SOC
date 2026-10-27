!                                         
! Copyright (C) 2001-2015 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------------------
SUBROUTINE adddvhubscf (ipert, ik)
  !--------------------------------------------------------------------------------------
  ! 
  !! DFPT+U  
  !! This routine calculates the SCF derivative of the Hubbard potential times \(\psi\):
  !! \begin{equation}\notag
  !! \begin{split}
  !!   |\Delta V_{SCF}(k+q,is) \psi(\text{ibnd},k,is)\rangle = 
  !!   - &\sum_{I,m1,m2} \text{Hubbard_U}(I)\cdot\text{dnsscf}(m1,m2,is,I,\text{imode})\cdot \\
  !!     &|S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2|\psi(\text{ibnd},k,is)\rangle
  !! \end{split}
  !! \end{equation}
  !
  !! Addition of the \(\text{J0}\) terms:
  !
  !! $$ + \sum_{I,m1,m2} \text{Hubbard_J0}(I)\cdot \text{dnsscf}(m1,m2,\text{isi},I,\text{imode})\cdot 
  !! |S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle $$
  !
  !! Where:  
  !! \(\text{is}\)  = current_spin;  
  !! \(\text{isi}\) = opposite of the current_spin.
  !
  !! Written by A. Floris. Modified by I. Timrov (01.10.2018).
  !
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : nwordwfcU
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, &
                            Hubbard_J0, nwfcU
  USE ldaU_lr,       ONLY : effU, swfcatomk, swfcatomkpq, dnsscf
  USE wvfct,         ONLY : npwx, nbnd 
  USE control_lr,    ONLY : lgamma, nbnd_occ
  USE lsda_mod,      ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions, ONLY : evc
  USE eqv,           ONLY : dvpsi  
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  USE klist,         ONLY : ngk
  USE buffers,       ONLY : get_buffer
  USE qpoint,        ONLY : xq, ikks, ikqs
  USE units_lr,      ONLY : iuatswfc
  USE noncollin_module,     ONLY : npol
  !  
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! the k point under consideration
  INTEGER, INTENT(IN) :: ipert  
  !! the index of perturbation
  !  
  ! Local variables
  !
  INTEGER :: npw, npwq, ikk, ikq, op_spin
  INTEGER ::  i, j, k, nt, l, ih, n, ig, ihubst, ihubst1, ihubst2, &
              nah, m, m1, m2, ibnd, ldim
  COMPLEX(DP), ALLOCATABLE :: dvhubscf(:,:), dvqi(:,:), proj1(:, :)
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  CALL start_clock ( 'adddvhubscf' )
  !
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (dvhubscf(npwx,nbnd))
  ALLOCATE (dvqi(npwx,nbnd))
  !  
  ldim = 2 * Hubbard_lmax + 1
  !
  ! At each iteration dvhubscf is set to zero, 
  ! because it is recomputed through dnsscf 
  !
  dvhubscf = (0.d0, 0.d0) 
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  IF (lsda) THEN
     current_spin = isk (ikk)
     IF (current_spin==1) THEN   
        op_spin = 2
     ELSE
        op_spin = 1
     END IF
  ELSE        
     op_spin = 1
  ENDIF
  !
  ! Read S * atomic wfc's at k and k+q on unit iuatswfc
  !  
  CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
  IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
  !
  DO nah = 1, nat
     !
     nt = ityp(nah)
     !
     IF (is_hubbard(nt)) THEN        
        !   
        ! Calculate proj1(ibnd,ihubst) = < psi(inbd,k) | S_{k}\phi_(k,I,m) >
        ldim = 2*Hubbard_l(nt) + 1
        ihubst = offsetU(nah) + 1   ! I m index
        CALL ZGEMM('C', 'N', nbnd, ldim, npw, (1.0d0, 0.0d0), &
                evc, npwx*npol, swfcatomk(1, ihubst), npwx*npol, &
                (0.0d0, 0.0d0), proj1(1, ihubst), nbnd)
        !
     ENDIF
     !
  ENDDO
  ! 
  CALL mp_sum (proj1, intra_bgrp_comm)
  !
  DO nah = 1, nat
     !
     nt = ityp(nah)
     !
     dvqi = (0.d0, 0.d0)
     !
     IF (is_hubbard(nt)) THEN
        !
        DO m1 = 1, 2*Hubbard_l(nt)+1
           !
           ihubst1 = offsetU(nah) + m1
           !
           DO m2 = 1, 2*Hubbard_l(nt)+1
              ! 
              ihubst2 = offsetU(nah) + m2
              ! 
              DO ibnd = 1, nbnd
                 DO ig = 1, npwq
                    dvqi(ig,ibnd) = dvqi(ig,ibnd) - effU(nt) * &
                                    dnsscf(m1,m2,current_spin,nah,ipert) * &
                                    swfcatomkpq(ig,ihubst1) * CONJG(proj1(ibnd,ihubst2))
                 ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
        DO ibnd = 1, nbnd
           DO ig = 1, npwq
              dvhubscf(ig,ibnd) = dvhubscf(ig,ibnd) + dvqi(ig,ibnd)
           ENDDO
        ENDDO
        !
     ENDIF
     !
     !  Hubbard_J0
     !
     dvqi = (0.d0, 0.d0)
     !
     IF (Hubbard_J0(nt).NE.0.d0) THEN
        !
        DO m1 = 1, 2*Hubbard_l(nt)+1
           !
           ihubst1 = offsetU(nah) + m1
           !
           DO m2 = 1, 2*Hubbard_l(nt)+1
              !
              ihubst2 = offsetU(nah) + m2
              ! 
              DO ibnd = 1, nbnd
                 DO ig = 1, npwq
                    ! Notice that sign change here
                    dvqi(ig,ibnd) = dvqi(ig,ibnd) + Hubbard_J0(nt) * &
                                    dnsscf(m1,m2,op_spin,nah,ipert) * &
                                    swfcatomkpq(ig,ihubst1) * CONJG(proj1(ibnd,ihubst2))
                 ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
        DO ibnd = 1, nbnd
           DO ig = 1, npwq
              dvhubscf(ig,ibnd) = dvhubscf(ig,ibnd) + dvqi(ig,ibnd)
           ENDDO
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !  
  ! Add the result to dvpsi
  !
  DO ibnd = 1, nbnd
     DO ig = 1, npwq
        dvpsi(ig,ibnd) = dvpsi(ig,ibnd) + dvhubscf(ig,ibnd)  
     ENDDO
  ENDDO
  !
  DEALLOCATE (proj1)
  DEALLOCATE (dvhubscf)
  DEALLOCATE (dvqi)
  !
  CALL stop_clock ( 'adddvhubscf' )
  !  
  RETURN
  ! 
END SUBROUTINE adddvhubscf

SUBROUTINE revert_mag_u (ns_nc)
   !
   ! This routine reverts the sign of the Hubbard magnetization,
   ! to be used in the noncollinear magnetic case.
   ! Written by L. Binci (2023)
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE lsda_mod,             ONLY : nspin
   USE ldaU,                 ONLY : lda_plus_u_kind, Hubbard_l, is_hubbard, &
                                    Hubbard_lmax
   USE noncollin_module,     ONLY : npol
   !
   IMPLICIT NONE
   INTEGER :: na, nt, m1, m2, is, is1, is2, ldim
   COMPLEX(DP) :: ns, ms(3)
   COMPLEX(DP), INTENT(INOUT) :: ns_nc(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat)
   !
   DO na = 1, nat
      nt = ityp(na)
      IF (is_hubbard(nt)) THEN
          ldim = 2 * Hubbard_l(nt) + 1
          DO m1 = 1, ldim
              DO m2 = 1, ldim
                   !
                   ! charge
                   ns = (ns_nc(m1,m2,1,na) + ns_nc(m1,m2,4,na))
                   ! magnetization
                   ms(1) =  (ns_nc(m1,m2,2,na)  + ns_nc(m1,m2,3,na))
                   ms(2) =  -(0.d0,1.0d0)*(ns_nc(m1,m2,2,na)  - ns_nc(m1,m2,3,na))
                   ms(3) =  (ns_nc(m1,m2,1,na)  - ns_nc(m1,m2,4,na))
                   !
                   ms(:) = -ms(:)
                   !
                   ns_nc(m1,m2,1,na) = 0.5d0*( ns + ms(3))
                   ns_nc(m1,m2,2,na) = 0.5d0*( ms(1) + (0.d0,1.0d0)*ms(2) )
                   ns_nc(m1,m2,3,na) = 0.5d0*( ms(1) - (0.d0,1.0d0)*ms(2) )
                   ns_nc(m1,m2,4,na) = 0.5d0*( ns - ms(3) )
              ENDDO
          ENDDO
      ENDIF
   ENDDO
   !
   RETURN
   !
 END SUBROUTINE revert_mag_u

 SUBROUTINE calc_vh_u( ns, v_hub)
   !
   !! Recomputes the noncollinear Hubbard potential.
   !! Used for a (spin-)Hubbard matrix with reversed magnetization
   !! Similar to v_hubbard_nc within PW/src/v_of_rho.f90
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U
   USE lsda_mod,             ONLY : nspin
   USE control_flags,        ONLY : iverbosity
   USE io_global,            ONLY : stdout
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
   !! Occupation matrix
   COMPLEX(DP), INTENT(INOUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
   !! Hubbard potential
   INTEGER  :: is, is1, na, nt, m1, m2
   ! ... local variables
   !
   v_hub(:,:,:,:) = 0.d0
   !
   DO na = 1, nat
      nt = ityp (na)
      IF (Hubbard_U(nt) /= 0.d0) THEN
         !
         ! is=1 and is=4 are diagonal components
         ! is=2 and is=3 are off-diagonal components
         DO is = 1, nspin
            IF (is == 2) THEN
             is1 = 3
            ELSEIF (is == 3) THEN
             is1 = 2
            ELSE
             is1 = is ! is=1 or is=4
            ENDIF
            !
            IF (is1 == is) THEN
               !
               ! Non spin-flip contribution (is=1 and is=4)
               ! (diagonal [spin indexes] occupancy matrices)
               DO m1 = 1, 2*Hubbard_l(nt) + 1
                  v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + 0.5D0*Hubbard_U(nt)
                  DO m2 = 1, 2 * Hubbard_l(nt) + 1
                     v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - Hubbard_U(nt) * ns(m2,m1,is,na)
                  ENDDO
               ENDDO
            ELSE
               !
               ! Spin-flip contribution (is=2 and is=3)
               ! (NON-diagonal [spin indexes] occupancy matrices)
               DO m1 = 1, 2*Hubbard_l(nt) + 1
                  DO m2 = 1, 2 * Hubbard_l(nt) + 1
                     v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - Hubbard_U(nt)*ns(m2,m1,is1,na)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO !is
      ENDIF
   ENDDO
   !
   RETURN
   !
 END SUBROUTINE calc_vh_u
 !--------------------------------------------------------------------------

 SUBROUTINE revert_mag_uv (nsg_nc)
   !
   ! This routine reverts the sign of the Hubbard magnetization,
   ! to be used in the noncollinear magnetic case.
   ! Written by L. Binci (2023)
   !
   USE kinds,                ONLY : DP
   USE ions_base,            ONLY : nat, ityp
   USE lsda_mod,             ONLY : nspin
   USE ldaU,                 ONLY : lda_plus_u_kind, Hubbard_l, is_hubbard, &
                                    Hubbard_lmax, ldmx_tot, max_num_neighbors,ldim_u,&
                                    neighood, at_sc
   USE noncollin_module,     ONLY : npol
   USE control_lr,           ONLY : nbnd_occ
   !
   IMPLICIT NONE
   INTEGER :: na, nt, m1, m2, is, is1, is2, ldim, ibnd, nt1, nt2, na1, na2, viz,&
              ldim1, ldim2, eq_na2, isp1, isp2
   INTEGER, EXTERNAL :: find_viz
   COMPLEX(DP) :: ns, ms(3)
   COMPLEX(DP), INTENT(INOUT) :: nsg_nc(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   !
   DO na1 = 1, nat
      nt1 = ityp (na1)
      IF ( ldim_u(nt1).GT.0 ) THEN
         ldim1 = ldim_u(nt1)
         DO viz = 1, neighood(na1)%num_neigh
            na2 = neighood(na1)%neigh(viz)
            eq_na2 = at_sc(na2)%at
            nt2 = ityp (eq_na2)
            ldim2 = ldim_u(nt2)
            !IF (na1.GT.na2) THEN
            !   DO m1 = 1, ldim1
            !      DO m2 = 1, ldim2
            !         DO is1 = 1, npol
            !            DO is2 = 1, npol
            !               isp1 = is2 + npol*(is1-1)
            !               isp2 = is1 + npol*(is2-1)
            !               nsg_nc(m2,m1,viz,na1,isp1) = &
            !                  CONJG(nsg_nc(m1,m2,find_viz(na2,na1),na2,isp2))
            !            ENDDO
            !         ENDDO
            !      ENDDO
            !   ENDDO
            !ELSE
               DO m1 = 1, ldim1
                  DO m2 = 1, ldim2
                     ! charge
                     ns = (nsg_nc(m2,m1,viz,na1,1) + nsg_nc(m2,m1,viz,na1,4))
                     ! magnetization
                     ms(1) = (nsg_nc(m2,m1,viz,na1,2) + nsg_nc(m2,m1,viz,na1,3))
                     ms(2) = -(0.d0,1.0d0)*(nsg_nc(m2,m1,viz,na1,2) - nsg_nc(m2,m1,viz,na1,3))
                     ms(3) = (nsg_nc(m2,m1,viz,na1,1) - nsg_nc(m2,m1,viz,na1,4))
                     !
                     ms(:) = -ms(:)
                     !
                     nsg_nc(m2,m1,viz,na1,1) = 0.5d0*( ns + ms(3))
                     nsg_nc(m2,m1,viz,na1,2) = 0.5d0*( ms(1) + (0.d0,1.0d0)*ms(2) )
                     nsg_nc(m2,m1,viz,na1,3) = 0.5d0*( ms(1) - (0.d0,1.0d0)*ms(2) )
                     nsg_nc(m2,m1,viz,na1,4) = 0.5d0*( ns - ms(3) )
                  ENDDO
               ENDDO
            !ENDIF
         ENDDO
      ENDIF
   ENDDO
   !
   RETURN
   !
 END SUBROUTINE revert_mag_uv

 SUBROUTINE calc_vh_uv( nsg, v_hub)
   !
   !! Recomputes the noncollinear Hubbard potential.
   !! Used for a (spin-)Hubbard matrix with reversed magnetization
   !! Similar to v_hubbard_nc within PW/src/v_of_rho.f90
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE ldaU,              ONLY : Hubbard_l, Hubbard_alpha, Hubbard_J0, Hubbard_beta,   &
                                 ldim_u, ldmx_tot, max_num_neighbors, at_sc, neighood, &
                                 Hubbard_V, Hubbard_alpha_back, is_hubbard, is_hubbard_back
   USE lsda_mod,          ONLY : nspin
   USE control_flags,     ONLY : iverbosity
   USE io_global,         ONLY : stdout
   USE noncollin_module,  ONLY : npol
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN)  :: nsg  (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   COMPLEX(DP), INTENT(OUT) :: v_hub(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   !
   ! Local variables
   !
   INTEGER :: is, is1, isop, na, na1, na2, nt, nt1, nt2, m1, m2, viz, equiv_na2
   INTEGER, EXTERNAL :: find_viz
   !
   v_hub(:,:,:,:,:) = (0.d0, 0.d0)
   !
   DO na1 = 1, nat
      nt1 = ityp(na1)
      IF ( is_hubbard(nt1) ) THEN
         DO is = 1, nspin
            IF (is == 2) THEN
               is1 = 3
            ELSEIF (is == 3) THEN
               is1 = 2
            ELSE
               is1 = is
            ENDIF
            DO viz = 1, neighood(na1)%num_neigh
               na2 = neighood(na1)%neigh(viz)
               equiv_na2 = at_sc(na2)%at
               nt2 = ityp(equiv_na2)
               !
               IF (is_hubbard(nt2) .AND. &
                   (Hubbard_V(na1,na2,1).NE.0.d0) ) THEN
                   !
                   ! Here no need to use is1: complex conjugation is enough
                   ! For both standard and background states of a center atom
                   DO m1 = 1, ldim_u(nt1)
                      ! For both standard and background states of the neighbor atom
                      DO m2 = 1, ldim_u(nt2)
                         v_hub(m2,m1,viz,na1,is) = - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,1)
                      ENDDO
                   ENDDO
                   !
                   IF ( na1.EQ.na2 .AND. is1.EQ.is) THEN
                      !
                      na = find_viz(na1,na1)
                      ! This is the diagonal term (like in the DFT+U only case)
                      DO m1 = 1, ldim_u(nt1)
                         v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) &
                                        + Hubbard_V(na1,na1,1) * 0.5d0
                      ENDDO
                      !
                   ENDIF
                   !
               ENDIF
               !
            ENDDO ! viz
            !
         ENDDO ! is
         !
      ENDIF
      !
   ENDDO ! na1
   !
   !
   RETURN
   !
 END SUBROUTINE calc_vh_uv
 !--------------------------------------------------------------------------

!----------------------------------------------------------------------------------------
SUBROUTINE adddvhubscf_nc (ipert, ik, time_reversed)
  !--------------------------------------------------------------------------------------
  !! Norm-conserving, fixed-Dudarev noncollinear response:
  !! apply delta V_U = -U_eff delta N to an ordinary or time-reversed
  !! Sternheimer branch using the PW full orbital-spin layout.
  !
  USE kinds,             ONLY : DP
  USE io_files,          ONLY : nwordwfcU
  USE ions_base,         ONLY : nat, ityp, ntyp => nsp
  USE ldaU,              ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, &
                                Hubbard_J0, Hubbard_alpha, Hubbard_beta, &
                                Hubbard_projectors, lda_plus_u_kind, nwfcU
  USE ldaU_lr,           ONLY : effU, swfcatomk, swfcatomkpq, dnsscf, &
                                pack_hubbard_blocks_nc, &
                                build_dvhub_nc, &
                                time_reverse_hubbard_matrix_nc, &
                                time_reverse_q_response_hubbard_matrix_nc, &
                                apply_hubbard_matrix_nc
  USE wvfct,             ONLY : npwx, nbnd
  USE control_lr,        ONLY : lgamma
  USE noncollin_module,  ONLY : npol
  USE wavefunctions,     ONLY : evc
  USE eqv,               ONLY : dvpsi
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_sum, mp_max
  USE qpoint,            ONLY : ikks, ikqs
  USE units_lr,          ONLY : iuatswfc
  USE buffers,           ONLY : get_buffer
  USE uspp,              ONLY : okvan
  USE control_flags,     ONLY : iverbosity
  USE io_global,         ONLY : stdout, ionode
  USE lsda_mod,          ONLY : nspin
  USE scf,               ONLY : rho
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ipert, ik
  LOGICAL, INTENT(IN) :: time_reversed
  !
  INTEGER :: ikk, ikq, na, nt, ldim_atom, max_matrix_size
  COMPLEX(DP), ALLOCATABLE :: projection(:,:), dvhub(:,:), dvhub_tr(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvhubpsi(:,:)
  COMPLEX(DP), POINTER :: output_projector(:,:)
  INTEGER :: is_bra, is_ket, iblock, jblock, m1, m2, matrix_size
  INTEGER :: saved_iverbosity
  REAL(DP), PARAMETER :: oracle_eta = 1.0e-3_DP
  COMPLEX(DP), PARAMETER :: ci = (0.0_DP, 1.0_DP)
  REAL(DP) :: eth_oracle
  REAL(DP) :: n_scale, mx_scale, my_scale, mz_scale
  REAL(DP) :: matrix_scale, matrix_residual, matrix_relative
  REAL(DP) :: action_scale, action_residual, action_relative
  LOGICAL :: oracle_active
  COMPLEX(DP) :: dns_value, dns_adjoint
  COMPLEX(DP), ALLOCATABLE :: dns_h0(:,:,:,:), dns_h90(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: ns_plus(:,:,:,:), ns_minus(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: v_plus(:,:,:,:), v_minus(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dv0(:,:,:,:), dv90(:,:,:,:), dvq_pw(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dvhub_pw(:,:), dvhub_minus_q(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvhub_pw_tr(:,:), dvhub_pw_used(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvhubpsi_pw(:,:)
  EXTERNAL :: v_hubbard_nc
  !
  CALL start_clock('adddvhubscf_nc')
  !
  IF (okvan) CALL errore('adddvhubscf_nc', &
                         'noncollinear Hubbard response requires NCPP', 1)
  IF (Hubbard_projectors /= 'atomic') CALL errore('adddvhubscf_nc', &
       'noncollinear Hubbard response requires atomic projectors', 1)
  IF (lda_plus_u_kind /= 0) CALL errore('adddvhubscf_nc', &
                                        'only Dudarev DFT+U is supported', 1)
  IF (npol /= 2) CALL errore('adddvhubscf_nc', &
                             'noncollinear Hubbard response requires npol=2', 1)
  DO nt = 1, ntyp
     IF (is_hubbard(nt) .AND. Hubbard_J0(nt) /= 0.0_DP) &
        CALL errore('adddvhubscf_nc', 'Hubbard J0 is outside the MVP', 1)
     IF (is_hubbard(nt) .AND. Hubbard_alpha(nt) /= 0.0_DP) &
        CALL errore('adddvhubscf_nc', 'Hubbard alpha is outside the MVP', 1)
     IF (is_hubbard(nt) .AND. Hubbard_beta(nt) /= 0.0_DP) &
        CALL errore('adddvhubscf_nc', 'Hubbard beta is outside the MVP', 1)
  ENDDO
  !
  max_matrix_size = (2 * Hubbard_lmax + 1) * npol
  ALLOCATE(projection(nwfcU,nbnd))
  ALLOCATE(dvhub(max_matrix_size,max_matrix_size))
  ALLOCATE(dvhub_tr(max_matrix_size,max_matrix_size))
  ALLOCATE(dvhubpsi(npwx*npol,nbnd))
  oracle_active = iverbosity == 1
  IF (oracle_active) THEN
     IF (nspin /= npol*npol) CALL errore('adddvhubscf_nc', &
          'PW Jacobian oracle requires four occupation spin blocks', 1)
     IF (.NOT.ALLOCATED(rho%ns_nc)) CALL errore('adddvhubscf_nc', &
          'PW Jacobian oracle requires ground-state ns_nc', 1)
     ALLOCATE(dns_h0(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
     ALLOCATE(dns_h90(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
     ALLOCATE(ns_plus(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
     ALLOCATE(ns_minus(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
     ALLOCATE(v_plus(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
     ALLOCATE(v_minus(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
     ALLOCATE(dv0(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
     ALLOCATE(dv90(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
     ALLOCATE(dvq_pw(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat))
     ALLOCATE(dvhub_pw(max_matrix_size,max_matrix_size))
     ALLOCATE(dvhub_minus_q(max_matrix_size,max_matrix_size))
     ALLOCATE(dvhub_pw_tr(max_matrix_size,max_matrix_size))
     ALLOCATE(dvhub_pw_used(max_matrix_size,max_matrix_size))
     ALLOCATE(dvhubpsi_pw(npwx*npol,nbnd))
     !
     ! Embed the complex q response into two Hermitian directions accepted
     ! by the ground-state PW Hubbard potential builder.
     !
     dns_h0 = (0.0_DP, 0.0_DP)
     dns_h90 = (0.0_DP, 0.0_DP)
     n_scale = 0.0_DP
     mx_scale = 0.0_DP
     my_scale = 0.0_DP
     mz_scale = 0.0_DP
     DO na = 1, nat
        nt = ityp(na)
        IF (is_hubbard(nt)) THEN
           ldim_atom = 2 * Hubbard_l(nt) + 1
           DO is_bra = 1, npol
              DO is_ket = 1, npol
                 iblock = npol * (is_bra - 1) + is_ket
                 jblock = npol * (is_ket - 1) + is_bra
                 DO m2 = 1, ldim_atom
                    DO m1 = 1, ldim_atom
                       dns_value = dnsscf(m1,m2,iblock,na,ipert)
                       dns_adjoint = CONJG( &
                            dnsscf(m2,m1,jblock,na,ipert))
                       dns_h0(m1,m2,iblock,na) = &
                            dns_value + dns_adjoint
                       dns_h90(m1,m2,iblock,na) = &
                            ci * (dns_value - dns_adjoint)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           n_scale = MAX(n_scale, MAXVAL(ABS( &
                dnsscf(1:ldim_atom,1:ldim_atom,1,na,ipert) + &
                dnsscf(1:ldim_atom,1:ldim_atom,4,na,ipert))))
           mx_scale = MAX(mx_scale, MAXVAL(ABS( &
                dnsscf(1:ldim_atom,1:ldim_atom,2,na,ipert) + &
                dnsscf(1:ldim_atom,1:ldim_atom,3,na,ipert))))
           my_scale = MAX(my_scale, MAXVAL(ABS(ci * ( &
                dnsscf(1:ldim_atom,1:ldim_atom,2,na,ipert) - &
                dnsscf(1:ldim_atom,1:ldim_atom,3,na,ipert)))))
           mz_scale = MAX(mz_scale, MAXVAL(ABS( &
                dnsscf(1:ldim_atom,1:ldim_atom,1,na,ipert) - &
                dnsscf(1:ldim_atom,1:ldim_atom,4,na,ipert))))
        ENDIF
     ENDDO
     !
     ns_plus = rho%ns_nc + oracle_eta * dns_h0
     ns_minus = rho%ns_nc - oracle_eta * dns_h0
     saved_iverbosity = iverbosity
     iverbosity = 0
     CALL v_hubbard_nc(ns_plus, v_plus, eth_oracle)
     CALL v_hubbard_nc(ns_minus, v_minus, eth_oracle)
     iverbosity = saved_iverbosity
     dv0 = (v_plus - v_minus) / (2.0_DP * oracle_eta)
     !
     ns_plus = rho%ns_nc + oracle_eta * dns_h90
     ns_minus = rho%ns_nc - oracle_eta * dns_h90
     saved_iverbosity = iverbosity
     iverbosity = 0
     CALL v_hubbard_nc(ns_plus, v_plus, eth_oracle)
     CALL v_hubbard_nc(ns_minus, v_minus, eth_oracle)
     iverbosity = saved_iverbosity
     dv90 = (v_plus - v_minus) / (2.0_DP * oracle_eta)
     dvq_pw = 0.5_DP * (dv0 - ci * dv90)
     dvhubpsi_pw = (0.0_DP, 0.0_DP)
     matrix_scale = 0.0_DP
     matrix_residual = 0.0_DP
  ENDIF
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  CALL get_buffer(swfcatomk, nwordwfcU, iuatswfc, ikk)
  IF (lgamma) THEN
     output_projector => swfcatomk
  ELSE
     CALL get_buffer(swfcatomkpq, nwordwfcU, iuatswfc, ikq)
     output_projector => swfcatomkpq
  ENDIF
  CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
             (1.0_DP,0.0_DP), swfcatomk, npwx*npol, evc, npwx*npol, &
             (0.0_DP,0.0_DP), projection, nwfcU)
  CALL mp_sum(projection, intra_bgrp_comm)
  !
  dvhubpsi = (0.0_DP, 0.0_DP)
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        ldim_atom = 2 * Hubbard_l(nt) + 1
        CALL build_dvhub_nc(ldim_atom, npol, effU(nt), &
                            dnsscf(:,:,:,na,ipert), dvhub)
        IF (time_reversed) THEN
           IF (lgamma) THEN
              CALL time_reverse_hubbard_matrix_nc(ldim_atom, npol, &
                                                   dvhub, dvhub_tr)
           ELSE
              CALL time_reverse_q_response_hubbard_matrix_nc( &
                   ldim_atom, npol, dvhub, dvhub_tr)
           ENDIF
           CALL apply_hubbard_matrix_nc(npwx*npol, nbnd, nwfcU, &
                ldim_atom, npol, offsetU(na), output_projector, projection, &
                dvhub_tr, dvhubpsi)
        ELSE
           CALL apply_hubbard_matrix_nc(npwx*npol, nbnd, nwfcU, &
                ldim_atom, npol, offsetU(na), output_projector, projection, &
                dvhub, dvhubpsi)
        ENDIF
        IF (oracle_active) THEN
           matrix_size = ldim_atom * npol
           CALL pack_hubbard_blocks_nc(ldim_atom, npol, &
                                       dvq_pw(:,:,:,na), dvhub_pw)
           IF (time_reversed) THEN
              IF (lgamma) THEN
                 CALL time_reverse_hubbard_matrix_nc(ldim_atom, npol, &
                                                      dvhub_pw, dvhub_pw_tr)
              ELSE
                 dvhub_minus_q = (0.0_DP, 0.0_DP)
                 dvhub_minus_q(1:matrix_size,1:matrix_size) = &
                      TRANSPOSE(CONJG( &
                           dvhub_pw(1:matrix_size,1:matrix_size)))
                 CALL time_reverse_hubbard_matrix_nc(ldim_atom, npol, &
                      dvhub_minus_q, dvhub_pw_tr)
              ENDIF
              dvhub_pw_used = dvhub_pw_tr
              matrix_scale = MAX(matrix_scale, MAXVAL(ABS( &
                   dvhub_pw_tr(1:matrix_size,1:matrix_size))))
              matrix_residual = MAX(matrix_residual, MAXVAL(ABS( &
                   dvhub_tr(1:matrix_size,1:matrix_size) - &
                   dvhub_pw_tr(1:matrix_size,1:matrix_size))))
           ELSE
              dvhub_pw_used = dvhub_pw
              matrix_scale = MAX(matrix_scale, MAXVAL(ABS( &
                   dvhub_pw(1:matrix_size,1:matrix_size))))
              matrix_residual = MAX(matrix_residual, MAXVAL(ABS( &
                   dvhub(1:matrix_size,1:matrix_size) - &
                   dvhub_pw(1:matrix_size,1:matrix_size))))
           ENDIF
           !
           ! Use intrinsic MATMUL as an implementation-independent check of
           ! the BLAS offset/layout in apply_hubbard_matrix_nc.
           !
           dvhubpsi_pw = dvhubpsi_pw + MATMUL( &
                output_projector(:,offsetU(na)+1:offsetU(na)+matrix_size), &
                MATMUL(dvhub_pw_used(1:matrix_size,1:matrix_size), &
                     projection(offsetU(na)+1:offsetU(na)+matrix_size,:)))
        ENDIF
     ENDIF
  ENDDO
  IF (oracle_active) THEN
     action_scale = MAXVAL(ABS(dvhubpsi_pw))
     action_residual = MAXVAL(ABS(dvhubpsi - dvhubpsi_pw))
     CALL mp_max(n_scale, intra_bgrp_comm)
     CALL mp_max(mx_scale, intra_bgrp_comm)
     CALL mp_max(my_scale, intra_bgrp_comm)
     CALL mp_max(mz_scale, intra_bgrp_comm)
     CALL mp_max(matrix_scale, intra_bgrp_comm)
     CALL mp_max(matrix_residual, intra_bgrp_comm)
     CALL mp_max(action_scale, intra_bgrp_comm)
     CALL mp_max(action_residual, intra_bgrp_comm)
     matrix_relative = matrix_residual / &
          MAX(matrix_scale, TINY(1.0_DP))
     action_relative = action_residual / &
          MAX(action_scale, TINY(1.0_DP))
     IF (ionode) WRITE(stdout, &
          '(5x,"SOC_U_PW_JACOBIAN_OPERATOR_CHECK",3(1x,i0),' // &
          '11(1x,es24.16))') ik, ipert, &
          MERGE(1,0,time_reversed), oracle_eta, n_scale, mx_scale, &
          my_scale, mz_scale, matrix_scale, matrix_residual, &
          matrix_relative, action_scale, action_residual, action_relative
  ENDIF
  !
  dvpsi(1:npwx*npol,1:nbnd) = dvpsi(1:npwx*npol,1:nbnd) + dvhubpsi
  !
  DEALLOCATE(projection)
  DEALLOCATE(dvhub)
  DEALLOCATE(dvhub_tr)
  DEALLOCATE(dvhubpsi)
  IF (oracle_active) THEN
     DEALLOCATE(dns_h0)
     DEALLOCATE(dns_h90)
     DEALLOCATE(ns_plus)
     DEALLOCATE(ns_minus)
     DEALLOCATE(v_plus)
     DEALLOCATE(v_minus)
     DEALLOCATE(dv0)
     DEALLOCATE(dv90)
     DEALLOCATE(dvq_pw)
     DEALLOCATE(dvhub_pw)
     DEALLOCATE(dvhub_minus_q)
     DEALLOCATE(dvhub_pw_tr)
     DEALLOCATE(dvhub_pw_used)
     DEALLOCATE(dvhubpsi_pw)
  ENDIF
  NULLIFY(output_projector)
  !
  CALL stop_clock('adddvhubscf_nc')
  !
  RETURN
END SUBROUTINE adddvhubscf_nc
