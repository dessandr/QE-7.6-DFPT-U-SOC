!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dvqhub_barepsi_us (ik, uact)
  !-----------------------------------------------------------------------
  !
  !! DFPT+U 
  !! This routines calculates the BARE derivative of the Hubbard potential
  !! times \(\text{psi}\).
  !
  ! is  = current_spin
  ! isi = opposite of the current_spin 
  !
  !! $$ |\Delta V_{BARE}(k+q,is) \psi(\text{ibnd},k,is)\rangle = 
  !!     + \sum_{I,m1,m2} \text{Hubbard_U}(I) \cdot [0.5\delta(m1,m2) - \text{ns}(m1,m2,is,I)]\cdot
  !!        [ |\text{dqsphi}(imode,I,k+q,m1)\rangle \langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle + 
  !!          |S\phi(I,k+q,m1)\rangle\langle\text{dmqsphi}(\text{imode},I,k,m2)|\psi(\text{ibnd},k,is)\rangle ]
  !!     - \sum_{I,m1,m2} \text{Hubbard_U}(I) \cdot \text{dnsbare}(m1,m2,is,I,\text{imode}) \cdot
  !!          |S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle $$
  !
  !! Addition of the J0 terms:
  !
  !! $$ + \sum_{I,m1,m2} \text{Hubbard_J0}(I)\cdot \text{ns}(m1,m2,isi,I)\cdot
  !!       [ |\text{dqsphi}(\text{imode},I,k+q,m1)\rangle\langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle + 
  !!         |S\phi(I,k+q,m1)\langle\rangle\text{dmqsphi}(\text{imode},I,k,m2)|\psi(\text{ibnd},k,is)\rangle ]
  !!    + \sum_{I,m1,m2} \text{Hubbard_J0}(I) \cdot \text{dnsbare}(m1,m2,isi,I,\text{imode}) \cdot
  !!                          |S\phi(I,k+q,m1)\rangle\langle S\phi(I,k,m2)|\psi(\text{ibnd},k,is)\rangle $$
  !
  !! Important: in this routine \(\text{vkb}\) is a beta function at k+q, and \(\text{vkb_}\) is beta at k.
  !! This is done so because \(\text{vkb}\) is calculated at k+q in solve_linter (i.e. before calling
  !! this routine), so we keep the same attribution here.
  ! 
  !! Written by A. Floris.  
  !! Modified by I. Timrov (01.10.2018).
  !                 
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode
  USE io_files,      ONLY : nwordwfcU
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE klist,         ONLY : xk, ngk, igk_k
  USE ldaU,          ONLY : Hubbard_l, is_hubbard, Hubbard_J0, offsetU, nwfcU
  USE ldaU_ph,       ONLY : wfcatomk, wfcatomkpq, dwfcatomkpq, &
                            sdwfcatomk, sdwfcatomkpq, dvkb, vkbkpq, dvkbkpq, dnsbare
  USE ldaU_lr,       ONLY : effU, swfcatomk, swfcatomkpq, proj1, proj2
  USE wvfct,         ONLY : npwx, nbnd
  USE uspp,          ONLY : vkb, nkb, ofsbeta
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE control_lr,    ONLY : lgamma
  USE units_lr,      ONLY : iuatwfc, iuatswfc
  USE uspp_param,    ONLY : nh
  USE lsda_mod,      ONLY : lsda, current_spin, isk
  USE wavefunctions, ONLY : evc
  USE eqv,           ONLY : dvpsi
  USE scf,           ONLY : rho
  USE mp_bands,      ONLY : intra_bgrp_comm, use_bgrp_in_hpsi
  USE mp,            ONLY : mp_sum
  USE buffers,       ONLY : get_buffer
  USE uspp_init,        ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! the k point under consideration
  COMPLEX(DP), INTENT(IN) :: uact(3*nat)
  !! the pattern of displacements
  !
  ! ... local variables
  !
  LOGICAL :: save_flag
  INTEGER :: i, j, k, icart, counter, na, nt, l, ih, n, mu, ig, &
             ihubst, ihubst1, ihubst2, nah, m, m1, m2, ibnd, op_spin, &
             ikk, ikq, npw, npwq, ibeta
  COMPLEX(DP), ALLOCATABLE :: aux1(:), aux2(:), aux3(:), aux4(:), aux5(:), &
                              dqsphi(:,:), dmqsphi(:,:), dvqi(:,:), dvqhbar(:,:,:,:), &
                              vkb_(:,:), dwfcatom_(:)
  !  
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (proj2(nbnd,nwfcU))
  ALLOCATE (aux1(npwx))
  ALLOCATE (aux2(npwx))
  ALLOCATE (aux3(npwx))
  ALLOCATE (aux4(npwx)) 
  ALLOCATE (aux5(npwx))
  ALLOCATE (dqsphi(npwx,nwfcU))
  ALLOCATE (dmqsphi(npwx,nwfcU))
  ALLOCATE (dvqi(npwx,nbnd))
  ALLOCATE (dvqhbar(npwx,nbnd,3,nat))
  ALLOCATE (vkb_(npwx,nkb))
  ALLOCATE (dwfcatom_(npwx))
  !
  save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi=.false.
  !
  proj1 = (0.d0, 0.d0)
  proj2 = (0.d0, 0.d0)
  ! 
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  !
  IF (lsda) THEN 
     current_spin = isk(ikk)
     IF (current_spin==1) THEN
        op_spin = 2
     ELSE
        op_spin = 1
     ENDIF
  ELSE        
     op_spin = 1
  ENDIF
  !
  ! Compute the beta function at k and put the result in vkb_
  !
  IF (.NOT.lgamma) THEN   
     CALL init_us_2 (npw, igk_k(1,ikk), xk(:,ikk), vkb_)
  ELSE
     vkb_ = vkb
  ENDIF
  !
  ! The beta function at k+q. Let us put it in the proper array vkbkpq,
  ! because in the following of this routine the array vkb will be 
  ! used as a workspace in the routine swfc. 
  !
  vkbkpq = vkb
  !
  ! Calculate the derivatives of beta functions 
  ! d^{icart}beta at k and k+q for all the bands and for 
  ! the 3 cartesian directions
  !
  DO icart = 1, 3
     DO na = 1, nat
        nt = ityp(na) 
        DO ih = 1, nh(nt)
           !
           ibeta = ofsbeta(na) + ih
           !
           CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                      vkb_(:,ibeta), dvkb(:,ibeta,icart))
           IF (.NOT.lgamma) &
           CALL dwfc (npwq, igk_k(1,ikq), ikq, icart, &
                      vkbkpq(:,ibeta), dvkbkpq(:,ibeta,icart))
           !
        ENDDO
     ENDDO
  ENDDO
  !
  ! Read \phi at k and k+q from file (unit iuatwfc)
  ! 
  CALL get_buffer (wfcatomk, nwordwfcU, iuatwfc, ikk)
  IF (.NOT.lgamma) CALL get_buffer (wfcatomkpq, nwordwfcU, iuatwfc, ikq)
  !
  ! Read S*\phi at k and k+q from file (unit iuatswfc)
  !
  CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
  IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
  ! 
  dvqhbar = (0.d0, 0.d0)
  !
  DO na = 1, nat
     !
     DO icart = 1, 3 
        !  
        dqsphi  = (0.d0, 0.d0)   
        dmqsphi = (0.d0, 0.d0)
        !
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           ! For Hubbard_U - Hubbard_J0
           ! 
           IF (is_hubbard(nt)) THEN
              !
              DO m = 1, 2*Hubbard_l(nt)+1
                 !
                 ihubst = offsetU(nah) + m   ! I m index
                 !
                 IF (nah==na) THEN
                    !
                    ! Calculate | d_icart\phi_(k,I,m)) >
                    !
                    CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                               wfcatomk(:,ihubst), dwfcatom_) 
                    !
                    ! Apply the S matrix: | S d_^(I,icart)\phi_(k,I,m) >
                    !
                    CALL swfc (npw, 1, vkb_, dwfcatom_, sdwfcatomk(:,ihubst))
                    ! 
                    IF (.NOT.lgamma) THEN
                       !
                       ! Calculate |d_icart\phi_(k+q,I,m)) >
                       !
                       CALL dwfc (npwq, igk_k(1,ikq), ikq, icart, &
                                  wfcatomkpq(:,ihubst), dwfcatom_)
                       ! 
                       ! Calculate | S d_^(I,icart)\phi_(k+q,I,m) >
                       !
                       CALL swfc (npwq, 1, vkbkpq, dwfcatom_, sdwfcatomkpq(:,ihubst))
                       ! 
                    ENDIF
                    ! 
                 ENDIF
                 !
                 ! Calculate |\Delta_q(S_k \phi_(k,I,m)) >  
                 ! and |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) >
                 !
                 CALL delta_sphi (ikk, ikq, na, icart, nah, ihubst, wfcatomk, wfcatomkpq,  &
                                  sdwfcatomk, sdwfcatomkpq, vkb_, vkbkpq, dvkb(:,:,icart), &
                                  dvkbkpq(:,:,icart), dqsphi, dmqsphi, 1)  
                 !
                 ! Calculate:
                 ! proj1 (ihubst,ibnd) = < S_{k}\phi_(k,I,m)| psi(ibnd,k) >
                 ! proj2 (ihubst,ibnd) = < \Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) | psi(ibnd,k) > 
                 !
                 !FIXME these are zgemms 
                 DO ibnd = 1, nbnd
                 proj1(ibnd,ihubst) = dot_product (swfcatomk(1:npw,ihubst), evc(1:npw,ibnd))
                 proj2(ibnd,ihubst) = dot_product(dmqsphi(1:npw,ihubst), evc(1:npw,ibnd))
                 ENDDO
                 !
              ENDDO ! m
              !
           ENDIF
           ! 
        ENDDO ! nah
        ! 
        CALL mp_sum (proj1, intra_bgrp_comm) 
        CALL mp_sum (proj2, intra_bgrp_comm)
        !
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           dvqi = (0.d0, 0.d0) 
           !
           IF (is_hubbard(nt)) THEN
              !
              DO ibnd = 1, nbnd
                 !
                 DO m1 = 1, 2*Hubbard_l(nt)+1
                    !
                    ihubst1 = offsetU(nah) + m1 
                    !
                    DO ig = 1, npwq
                       !
                       aux1(ig) = dqsphi(ig,ihubst1) * proj1(ibnd,ihubst1) 
                       !
                       aux3(ig) = swfcatomkpq(ig,ihubst1) * proj2(ibnd,ihubst1)  
                       !
                       dvqi(ig,ibnd) = dvqi(ig,ibnd) + 0.5d0 * (aux1(ig)+aux3(ig))
                       !
                    ENDDO
                    !
                    DO m2 = 1, 2*Hubbard_l(nt)+1
                       !
                       ihubst2 = offsetU(nah) + m2
                       !
                       DO ig = 1, npwq 
                          !                         
                          aux2(ig) = dqsphi(ig,ihubst1) * rho%ns(m1,m2,current_spin,nah) &
                                     * proj1(ibnd, ihubst2)
                          aux4(ig) = swfcatomkpq(ig,ihubst1) * rho%ns(m1,m2,current_spin,nah) &
                                     * proj2(ibnd, ihubst2)
                          aux5(ig) = swfcatomkpq(ig,ihubst1) &
                                     * dnsbare(m1,m2,current_spin,nah,icart,na) &
                                     * proj1(ibnd, ihubst2)
                          !
                          dvqi(ig, ibnd) = dvqi(ig, ibnd) - aux2(ig) - aux4(ig) - aux5(ig) 
                          ! 
                       ENDDO
                       ! 
                    ENDDO ! m2
                    !
                 ENDDO ! m1
                 !
              ENDDO ! ibnd
              !
              ! effU = Hubbard_U - Hubbard_J0
              !
              dvqi = dvqi * effU(nt)
              ! 
              DO ig = 1, npwq
                 dvqhbar(ig,:,icart,na) = dvqhbar(ig,:,icart,na) + dvqi(ig,:)
              ENDDO
              ! 
           ENDIF
           !
           ! Hubbard_J0 \= 0 case 
           !
           dvqi = (0.d0, 0.d0) 
           ! 
           IF (Hubbard_J0(nt).NE.0.d0) THEN
              !              
              DO ibnd = 1, nbnd
                 ! 
                 DO m1 = 1, 2*Hubbard_l(nt)+1
                    ! 
                    ihubst1 = offsetU(nah) + m1 
                    !
                    ! No diagonal term for J0
                    !                          
                    DO m2 = 1, 2*Hubbard_l(nt)+1
                       ! 
                       ihubst2 = offsetU(nah) + m2
                       ! 
                       DO ig = 1, npwq                          
                          aux2(ig) = dqsphi(ig, ihubst1) * rho%ns(m1,m2,op_spin,nah) &
                                     * proj1(ibnd, ihubst2)
                          aux4(ig) = swfcatomkpq(ig,ihubst1) * rho%ns(m1,m2,op_spin,nah) &
                                     * proj2(ibnd, ihubst2)
                          aux5(ig) = swfcatomkpq(ig,ihubst1) &
                                     * dnsbare (m1,m2,op_spin,nah,icart,na) &
                                     * proj1 (ibnd, ihubst2)
                          !
                          ! Note the sign change w.r.t. the case above
                          !
                          dvqi(ig, ibnd) = dvqi(ig, ibnd) + aux2(ig) + aux4(ig) + aux5(ig)
                          ! 
                       ENDDO
                       !
                    ENDDO ! m2
                    !
                 ENDDO ! m1
                 !
              ENDDO ! ibnd
              !
              dvqi = dvqi * Hubbard_J0(nt)
              !   
              DO ig = 1, npwq
                 dvqhbar(ig,:,icart,na) = dvqhbar(ig,:,icart,na) + dvqi(ig,:)
              ENDDO
              !
           ENDIF
           !
        ENDDO ! nah
        !
     ENDDO ! icart
     ! 
  ENDDO ! na
  ! 
  ! Compute the displacement along the pattern uact (for each band).
  ! The result is stored in aux1.
  !
  DO ibnd = 1, nbnd    
     !   
     aux1 = (0.d0, 0.d0)
     !
     DO na = 1, nat
        mu = 3 * (na - 1)
        ! Here is the basis transformation from cartesian to pattern uact
        DO icart = 1, 3 
           DO ig = 1, npwq
              aux1(ig) = aux1(ig) + dvqhbar(ig,ibnd,icart,na) * uact(mu+icart) 
           ENDDO
        ENDDO
     ENDDO
     !
     ! Add the result to dvpsi
     !
     DO ig = 1, npwq
        dvpsi(ig,ibnd) = dvpsi(ig,ibnd) + aux1(ig)
     ENDDO
     !
  ENDDO
  !
  use_bgrp_in_hpsi = save_flag
  !
  DEALLOCATE (proj1)
  DEALLOCATE (proj2)  
  DEALLOCATE (aux1) 
  DEALLOCATE (aux2)
  DEALLOCATE (aux3) 
  DEALLOCATE (aux4)
  DEALLOCATE (aux5)
  DEALLOCATE (dqsphi)
  DEALLOCATE (dmqsphi)
  DEALLOCATE (dvqi)
  DEALLOCATE (dvqhbar)
  DEALLOCATE (vkb_)
  DEALLOCATE (dwfcatom_) 
  !
  RETURN
  !
END SUBROUTINE dvqhub_barepsi_us

!-----------------------------------------------------------------------
SUBROUTINE dvqhub_barepsi_nc (ik, uact, time_reversed, run_bare_rhs_oracle)
  !-----------------------------------------------------------------------
  !! Norm-conserving, fixed-Dudarev noncollinear bare Hubbard
  !! perturbation:
  !!
  !!   delta(P V P^dagger) psi =
  !!     deltaP V (P^dagger psi)
  !!     + P V (deltaP^dagger psi)
  !!     + P deltaV_bare (P^dagger psi).
  !!
  !! The ground-state V is reused from PW, while deltaV_bare contains only
  !! the derivative map -U_eff deltaN_bare.
  !
  USE kinds,             ONLY : DP
  USE io_files,          ONLY : nwordwfcU
  USE ions_base,         ONLY : nat, ityp, ntyp => nsp
  USE klist,             ONLY : ngk, igk_k, xk
  USE ldaU,              ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, &
                                Hubbard_J0, Hubbard_alpha, Hubbard_beta, &
                                Hubbard_projectors, lda_plus_u_kind, nwfcU
  USE ldaU_ph,           ONLY : wfcatomk, wfcatomkpq, dnsbare
  USE ldaU_lr,           ONLY : effU, swfcatomk, swfcatomkpq, &
                                pack_hubbard_blocks_nc, build_dvhub_nc, &
                                time_reverse_hubbard_matrix_nc, &
                                time_reverse_q_response_hubbard_matrix_nc, &
                                apply_hubbard_matrix_nc
  USE wvfct,             ONLY : npwx, nbnd
  USE qpoint,            ONLY : ikks, ikqs
  USE control_lr,        ONLY : lgamma, nbnd_occ
  USE control_flags,     ONLY : iverbosity
  USE units_lr,          ONLY : iuatwfc, iuatswfc
  USE wavefunctions,     ONLY : evc
  USE eqv,               ONLY : dvpsi
  USE scf,               ONLY : v
  USE mp_bands,          ONLY : intra_bgrp_comm, nbgrp
  USE mp,                ONLY : mp_sum, mp_max
  USE buffers,           ONLY : get_buffer
  USE noncollin_module,  ONLY : npol
  USE uspp,              ONLY : okvan
  USE io_global,         ONLY : stdout, ionode
  USE cell_base,         ONLY : tpiba
  USE gvect,             ONLY : g, mill
  USE lsda_mod,          ONLY : nspin
  USE scf,               ONLY : rho
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  COMPLEX(DP), INTENT(IN) :: uact(3*nat)
  LOGICAL, INTENT(IN) :: time_reversed
  LOGICAL, INTENT(IN) :: run_bare_rhs_oracle
  !
  INTEGER :: ikk, ikq, npw, npwq, na, na_disp, nt, icart, mu, is
  INTEGER :: m, iwfc, is_bra
  INTEGER :: ldim_atom, max_ldim, max_matrix_size
  INTEGER :: row_first
  COMPLEX(DP), ALLOCATABLE, TARGET :: dprojector_k(:,:), &
                                      dprojector_kpq(:,:)
  COMPLEX(DP), ALLOCATABLE :: dcolumn(:)
  COMPLEX(DP), POINTER :: output_projector(:,:), doutput_projector(:,:)
  COMPLEX(DP), ALLOCATABLE :: projection(:,:), dprojection(:,:)
  COMPLEX(DP), ALLOCATABLE :: ground_matrix(:,:), ground_matrix_tr(:,:)
  COMPLEX(DP), ALLOCATABLE :: bare_matrix(:,:), bare_matrix_tr(:,:)
  COMPLEX(DP), ALLOCATABLE :: dns_pattern(:,:,:), dvhubpsi(:,:)
  REAL(DP) :: projection_scale, dprojection_scale
  REAL(DP) :: output_projector_scale, doutput_projector_scale
  REAL(DP) :: result_scale
  INTEGER :: is_ket, iblock, jblock, m1, m2
  INTEGER :: saved_iverbosity
  REAL(DP), PARAMETER :: bare_oracle_eta = 1.0e-4_DP
  REAL(DP), PARAMETER :: bare_oracle_pattern_tol = 1.0e-12_DP
  REAL(DP), PARAMETER :: potential_oracle_eta = 1.0e-3_DP
  COMPLEX(DP), PARAMETER :: ci = (0.0_DP, 1.0_DP)
  REAL(DP) :: eth_oracle, pattern_imag_scale
  REAL(DP) :: projector_k_scale, projector_k_residual
  REAL(DP) :: projector_k_relative, projector_kpq_scale
  REAL(DP) :: projector_kpq_residual, projector_kpq_relative
  REAL(DP) :: analytic_scale, coarse_residual, coarse_relative
  REAL(DP) :: richardson_residual, richardson_relative
  REAL(DP) :: phase_argument
  COMPLEX(DP) :: dns_value, dns_adjoint, numerical_derivative
  LOGICAL :: oracle_active, oracle_requested
  COMPLEX(DP), ALLOCATABLE :: dns_q_oracle(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dns_h0(:,:,:,:), dns_h90(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: ns_plus(:,:,:,:), ns_minus(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: v_plus(:,:,:,:), v_minus(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: v_ground_pw(:,:,:,:), dv0(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dv90(:,:,:,:), dvq_pw(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: projector_k_fd(:,:)
  COMPLEX(DP), ALLOCATABLE :: projector_kpq_fd(:,:)
  COMPLEX(DP), ALLOCATABLE :: projection_fd(:,:)
  COMPLEX(DP), ALLOCATABLE :: action_plus(:,:), action_minus(:,:)
  COMPLEX(DP), ALLOCATABLE :: action_half_plus(:,:)
  COMPLEX(DP), ALLOCATABLE :: action_half_minus(:,:)
  COMPLEX(DP), ALLOCATABLE :: derivative_coarse(:,:)
  COMPLEX(DP), ALLOCATABLE :: derivative_half(:,:)
  COMPLEX(DP), ALLOCATABLE :: derivative_richardson(:,:)
  COMPLEX(DP), ALLOCATABLE :: oracle_ground_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: oracle_ground_matrix_tr(:,:)
  COMPLEX(DP), ALLOCATABLE :: oracle_dv_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: oracle_dv_minus_q(:,:)
  COMPLEX(DP), ALLOCATABLE :: oracle_dv_matrix_tr(:,:)
  COMPLEX(DP), ALLOCATABLE :: oracle_matrix_used(:,:)
  COMPLEX(DP) :: term_checksum_previous(2,2)
  EXTERNAL :: v_hubbard_nc
  !
  CALL start_clock('dvqhub_barepsi_nc')
  !
  IF (okvan) CALL errore('dvqhub_barepsi_nc', &
                         'noncollinear Hubbard response requires NCPP', 1)
  IF (Hubbard_projectors /= 'atomic') CALL errore('dvqhub_barepsi_nc', &
       'noncollinear Hubbard response requires atomic projectors', 1)
  IF (lda_plus_u_kind /= 0) CALL errore('dvqhub_barepsi_nc', &
                                        'only Dudarev DFT+U is supported', 1)
  IF (npol /= 2) CALL errore('dvqhub_barepsi_nc', &
                             'noncollinear Hubbard response requires npol=2', 1)
  IF (.NOT.ALLOCATED(v%ns_nc)) CALL errore('dvqhub_barepsi_nc', &
                                           'PW Hubbard potential is unavailable', 1)
  DO nt = 1, ntyp
     IF (is_hubbard(nt) .AND. Hubbard_J0(nt) /= 0.0_DP) &
        CALL errore('dvqhub_barepsi_nc', 'Hubbard J0 is outside the MVP', 1)
     IF (is_hubbard(nt) .AND. Hubbard_alpha(nt) /= 0.0_DP) &
        CALL errore('dvqhub_barepsi_nc', 'Hubbard alpha is outside the MVP', 1)
     IF (is_hubbard(nt) .AND. Hubbard_beta(nt) /= 0.0_DP) &
        CALL errore('dvqhub_barepsi_nc', 'Hubbard beta is outside the MVP', 1)
  ENDDO
  !
  max_ldim = 2 * Hubbard_lmax + 1
  max_matrix_size = max_ldim * npol
  ALLOCATE(dprojector_k(npwx*npol,nwfcU))
  ALLOCATE(dprojector_kpq(npwx*npol,nwfcU))
  ALLOCATE(dcolumn(npwx*npol))
  ALLOCATE(projection(nwfcU,nbnd))
  ALLOCATE(dprojection(nwfcU,nbnd))
  ALLOCATE(ground_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(ground_matrix_tr(max_matrix_size,max_matrix_size))
  ALLOCATE(bare_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(bare_matrix_tr(max_matrix_size,max_matrix_size))
  ALLOCATE(dns_pattern(max_ldim,max_ldim,npol*npol))
  ALLOCATE(dvhubpsi(npwx*npol,nbnd))
  pattern_imag_scale = MAXVAL(ABS(AIMAG(uact)))
  oracle_requested = iverbosity == 1 .AND. run_bare_rhs_oracle
  IF (oracle_requested .AND. nbgrp /= 1) &
     CALL errore('dvqhub_barepsi_nc', &
          'bare term checksum requires one band group', 1)
  oracle_active = oracle_requested .AND. &
                  pattern_imag_scale <= bare_oracle_pattern_tol
  IF (oracle_requested .AND. .NOT.oracle_active .AND. ionode .AND. &
      ik == 1 .AND. .NOT.time_reversed) THEN
     WRITE(stdout, &
          '(5x,"SOC_U_PW_BARE_RHS_ORACLE_SKIPPED",2(1x,es24.16))') &
          pattern_imag_scale, bare_oracle_pattern_tol
  ENDIF
  IF (oracle_active) THEN
     IF (nspin /= npol*npol) CALL errore('dvqhub_barepsi_nc', &
          'bare RHS oracle requires four occupation spin blocks', 1)
     IF (.NOT.ALLOCATED(rho%ns_nc)) CALL errore('dvqhub_barepsi_nc', &
          'bare RHS oracle requires ground-state ns_nc', 1)
     ALLOCATE(dns_q_oracle(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(dns_h0(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(dns_h90(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(ns_plus(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(ns_minus(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(v_plus(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(v_minus(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(v_ground_pw(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(dv0(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(dv90(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(dvq_pw(max_ldim,max_ldim,nspin,nat))
     ALLOCATE(projector_k_fd(npwx*npol,nwfcU))
     ALLOCATE(projector_kpq_fd(npwx*npol,nwfcU))
     ALLOCATE(projection_fd(nwfcU,nbnd))
     ALLOCATE(action_plus(npwx*npol,nbnd))
     ALLOCATE(action_minus(npwx*npol,nbnd))
     ALLOCATE(action_half_plus(npwx*npol,nbnd))
     ALLOCATE(action_half_minus(npwx*npol,nbnd))
     ALLOCATE(derivative_coarse(npwx*npol,nbnd))
     ALLOCATE(derivative_half(npwx*npol,nbnd))
     ALLOCATE(derivative_richardson(npwx*npol,nbnd))
     ALLOCATE(oracle_ground_matrix(max_matrix_size,max_matrix_size))
     ALLOCATE(oracle_ground_matrix_tr(max_matrix_size,max_matrix_size))
     ALLOCATE(oracle_dv_matrix(max_matrix_size,max_matrix_size))
     ALLOCATE(oracle_dv_minus_q(max_matrix_size,max_matrix_size))
     ALLOCATE(oracle_dv_matrix_tr(max_matrix_size,max_matrix_size))
     ALLOCATE(oracle_matrix_used(max_matrix_size,max_matrix_size))
  ENDIF
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq = ngk(ikq)
  CALL get_buffer(wfcatomk, nwordwfcU, iuatwfc, ikk)
  CALL get_buffer(swfcatomk, nwordwfcU, iuatswfc, ikk)
  IF (.NOT.lgamma) THEN
     CALL get_buffer(wfcatomkpq, nwordwfcU, iuatwfc, ikq)
     CALL get_buffer(swfcatomkpq, nwordwfcU, iuatswfc, ikq)
  ENDIF
  !
  dprojector_k = (0.0_DP, 0.0_DP)
  dprojector_kpq = (0.0_DP, 0.0_DP)
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        ldim_atom = 2 * Hubbard_l(nt) + 1
        mu = 3 * (na - 1)
        DO is_bra = 1, npol
           DO m = 1, ldim_atom
              iwfc = offsetU(na) + m + ldim_atom * (is_bra - 1)
              DO icart = 1, 3
                 dcolumn = (0.0_DP, 0.0_DP)
                 DO is = 1, npol
                    row_first = 1 + npwx * (is - 1)
                    CALL dwfc(npw, igk_k(1,ikk), ikk, icart, &
                              wfcatomk(row_first,iwfc), &
                              dcolumn(row_first))
                 ENDDO
                 ! dprojector_k is conjugate-transposed below.  Store the
                 ! conjugated pattern coefficient here so that the resulting
                 ! bra projection remains linear in uact.
                 dprojector_k(:,iwfc) = dprojector_k(:,iwfc) + &
                                        CONJG(uact(mu+icart)) * dcolumn
                 IF (lgamma) THEN
                    ! The ket-side derivative is linear in uact and must
                    ! remain distinct from the bra-side auxiliary even when
                    ! k+q and k coincide.
                    dprojector_kpq(:,iwfc) = &
                         dprojector_kpq(:,iwfc) + &
                         uact(mu+icart) * dcolumn
                 ELSE
                    dcolumn = (0.0_DP, 0.0_DP)
                    DO is = 1, npol
                       row_first = 1 + npwx * (is - 1)
                       CALL dwfc(npwq, igk_k(1,ikq), ikq, icart, &
                                 wfcatomkpq(row_first,iwfc), &
                                 dcolumn(row_first))
                    ENDDO
                    dprojector_kpq(:,iwfc) = &
                         dprojector_kpq(:,iwfc) + &
                         uact(mu+icart) * dcolumn
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  IF (lgamma) THEN
     output_projector => swfcatomk
  ELSE
     output_projector => swfcatomkpq
  ENDIF
  doutput_projector => dprojector_kpq
  !
  CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
             (1.0_DP,0.0_DP), swfcatomk, npwx*npol, evc, npwx*npol, &
             (0.0_DP,0.0_DP), projection, nwfcU)
  CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
             (1.0_DP,0.0_DP), dprojector_k, npwx*npol, evc, npwx*npol, &
             (0.0_DP,0.0_DP), dprojection, nwfcU)
  CALL mp_sum(projection, intra_bgrp_comm)
  CALL mp_sum(dprojection, intra_bgrp_comm)
  IF (oracle_active) THEN
     !
     ! Verify dwfc independently against the exact translation phase of
     ! the NCPP atomic projectors.
     !
     projector_k_scale = 0.0_DP
     projector_k_residual = 0.0_DP
     projector_kpq_scale = 0.0_DP
     projector_kpq_residual = 0.0_DP
     DO na = 1, nat
        nt = ityp(na)
        IF (is_hubbard(nt)) THEN
           ldim_atom = 2 * Hubbard_l(nt) + 1
           mu = 3 * (na - 1)
           DO is_bra = 1, npol
              row_first = 1 + npwx * (is_bra - 1)
              DO m = 1, ldim_atom
                 iwfc = offsetU(na) + m + ldim_atom * (is_bra - 1)
                 DO is = 1, npw
                    phase_argument = 0.0_DP
                    DO icart = 1, 3
                       phase_argument = phase_argument + &
                            (g(icart,igk_k(is,ikk)) + xk(icart,ikk)) * &
                            tpiba * DBLE(uact(mu+icart))
                    ENDDO
                    numerical_derivative = -ci * &
                         SIN(bare_oracle_eta * phase_argument) / &
                         bare_oracle_eta * swfcatomk(row_first+is-1,iwfc)
                    projector_k_scale = MAX(projector_k_scale, &
                         ABS(dprojector_k(row_first+is-1,iwfc)))
                    projector_k_residual = MAX(projector_k_residual, &
                         ABS(dprojector_k(row_first+is-1,iwfc) - &
                              numerical_derivative))
                 ENDDO
                 DO is = 1, npwq
                    phase_argument = 0.0_DP
                    DO icart = 1, 3
                       phase_argument = phase_argument + &
                            (g(icart,igk_k(is,ikq)) + xk(icart,ikq)) * &
                            tpiba * DBLE(uact(mu+icart))
                    ENDDO
                    numerical_derivative = -ci * &
                         SIN(bare_oracle_eta * phase_argument) / &
                         bare_oracle_eta * &
                         output_projector(row_first+is-1,iwfc)
                    projector_kpq_scale = MAX(projector_kpq_scale, &
                         ABS(dprojector_kpq(row_first+is-1,iwfc)))
                    projector_kpq_residual = MAX( &
                         projector_kpq_residual, &
                         ABS(dprojector_kpq(row_first+is-1,iwfc) - &
                              numerical_derivative))
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     !
     ! Reconstruct the complex-q bare occupation response through two
     ! Hermitian directions accepted by the PW ground-state potential.
     !
     dns_q_oracle = (0.0_DP, 0.0_DP)
     dns_h0 = (0.0_DP, 0.0_DP)
     dns_h90 = (0.0_DP, 0.0_DP)
     DO na = 1, nat
        nt = ityp(na)
        IF (is_hubbard(nt)) THEN
           ldim_atom = 2 * Hubbard_l(nt) + 1
           DO na_disp = 1, nat
              DO icart = 1, 3
                 dns_q_oracle(:,:,1:nspin,na) = &
                      dns_q_oracle(:,:,1:nspin,na) + &
                      uact(3*(na_disp-1)+icart) * &
                      dnsbare(:,:,1:nspin,na,icart,na_disp)
              ENDDO
           ENDDO
           DO is_bra = 1, npol
              DO is_ket = 1, npol
                 iblock = npol * (is_bra - 1) + is_ket
                 jblock = npol * (is_ket - 1) + is_bra
                 DO m2 = 1, ldim_atom
                    DO m1 = 1, ldim_atom
                       dns_value = dns_q_oracle(m1,m2,iblock,na)
                       dns_adjoint = CONJG( &
                            dns_q_oracle(m2,m1,jblock,na))
                       dns_h0(m1,m2,iblock,na) = &
                            dns_value + dns_adjoint
                       dns_h90(m1,m2,iblock,na) = &
                            ci * (dns_value - dns_adjoint)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     saved_iverbosity = iverbosity
     iverbosity = 0
     CALL v_hubbard_nc(rho%ns_nc, v_ground_pw, eth_oracle)
     ns_plus = rho%ns_nc + potential_oracle_eta * dns_h0
     ns_minus = rho%ns_nc - potential_oracle_eta * dns_h0
     CALL v_hubbard_nc(ns_plus, v_plus, eth_oracle)
     CALL v_hubbard_nc(ns_minus, v_minus, eth_oracle)
     dv0 = (v_plus - v_minus) / (2.0_DP * potential_oracle_eta)
     ns_plus = rho%ns_nc + potential_oracle_eta * dns_h90
     ns_minus = rho%ns_nc - potential_oracle_eta * dns_h90
     CALL v_hubbard_nc(ns_plus, v_plus, eth_oracle)
     CALL v_hubbard_nc(ns_minus, v_minus, eth_oracle)
     dv90 = (v_plus - v_minus) / (2.0_DP * potential_oracle_eta)
     iverbosity = saved_iverbosity
     dvq_pw = 0.5_DP * (dv0 - ci * dv90)
  ENDIF
  !
  dvhubpsi = (0.0_DP, 0.0_DP)
  term_checksum_previous = (0.0_DP, 0.0_DP)
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        ldim_atom = 2 * Hubbard_l(nt) + 1
        CALL pack_hubbard_blocks_nc(ldim_atom, npol, &
                                     v%ns_nc(:,:,:,na), ground_matrix)
        !
        dns_pattern = (0.0_DP, 0.0_DP)
        DO na_disp = 1, nat
           DO icart = 1, 3
              dns_pattern = dns_pattern + &
                   uact(3*(na_disp-1)+icart) * &
                   dnsbare(:,:,:,na,icart,na_disp)
           ENDDO
        ENDDO
        CALL build_dvhub_nc(ldim_atom, npol, effU(nt), &
                            dns_pattern, bare_matrix)
        !
        IF (time_reversed) THEN
           CALL time_reverse_hubbard_matrix_nc(ldim_atom, npol, &
                                                ground_matrix, ground_matrix_tr)
           IF (lgamma) THEN
              CALL time_reverse_hubbard_matrix_nc(ldim_atom, npol, &
                                                   bare_matrix, bare_matrix_tr)
           ELSE
              CALL time_reverse_q_response_hubbard_matrix_nc( &
                   ldim_atom, npol, bare_matrix, bare_matrix_tr)
           ENDIF
           CALL apply_hubbard_matrix_nc(npwx*npol, nbnd, nwfcU, &
                ldim_atom, npol, offsetU(na), doutput_projector, projection, &
                ground_matrix_tr, dvhubpsi)
           IF (oracle_requested) CALL report_bare_term_checksum(na, 1)
           CALL apply_hubbard_matrix_nc(npwx*npol, nbnd, nwfcU, &
                ldim_atom, npol, offsetU(na), output_projector, dprojection, &
                ground_matrix_tr, dvhubpsi)
           IF (oracle_requested) CALL report_bare_term_checksum(na, 2)
           CALL apply_hubbard_matrix_nc(npwx*npol, nbnd, nwfcU, &
                ldim_atom, npol, offsetU(na), output_projector, projection, &
                bare_matrix_tr, dvhubpsi)
           IF (oracle_requested) CALL report_bare_term_checksum(na, 3)
        ELSE
           CALL apply_hubbard_matrix_nc(npwx*npol, nbnd, nwfcU, &
                ldim_atom, npol, offsetU(na), doutput_projector, projection, &
                ground_matrix, dvhubpsi)
           IF (oracle_requested) CALL report_bare_term_checksum(na, 1)
           CALL apply_hubbard_matrix_nc(npwx*npol, nbnd, nwfcU, &
                ldim_atom, npol, offsetU(na), output_projector, dprojection, &
                ground_matrix, dvhubpsi)
           IF (oracle_requested) CALL report_bare_term_checksum(na, 2)
           CALL apply_hubbard_matrix_nc(npwx*npol, nbnd, nwfcU, &
                ldim_atom, npol, offsetU(na), output_projector, projection, &
                bare_matrix, dvhubpsi)
           IF (oracle_requested) CALL report_bare_term_checksum(na, 3)
        ENDIF
     ENDIF
  ENDDO
  !
  IF (oracle_active) THEN
     CALL build_bare_fd_action(bare_oracle_eta, action_plus)
     CALL build_bare_fd_action(-bare_oracle_eta, action_minus)
     CALL build_bare_fd_action(0.5_DP * bare_oracle_eta, &
                               action_half_plus)
     CALL build_bare_fd_action(-0.5_DP * bare_oracle_eta, &
                               action_half_minus)
     derivative_coarse = (action_plus - action_minus) / &
                         (2.0_DP * bare_oracle_eta)
     derivative_half = (action_half_plus - action_half_minus) / &
                       bare_oracle_eta
     derivative_richardson = (4.0_DP * derivative_half - &
                              derivative_coarse) / 3.0_DP
     analytic_scale = MAXVAL(ABS(dvhubpsi))
     coarse_residual = MAXVAL(ABS(dvhubpsi - derivative_coarse))
     richardson_residual = MAXVAL(ABS( &
          dvhubpsi - derivative_richardson))
     CALL mp_max(projector_k_scale, intra_bgrp_comm)
     CALL mp_max(projector_k_residual, intra_bgrp_comm)
     CALL mp_max(projector_kpq_scale, intra_bgrp_comm)
     CALL mp_max(projector_kpq_residual, intra_bgrp_comm)
     CALL mp_max(analytic_scale, intra_bgrp_comm)
     CALL mp_max(coarse_residual, intra_bgrp_comm)
     CALL mp_max(richardson_residual, intra_bgrp_comm)
     projector_k_relative = projector_k_residual / &
          MAX(projector_k_scale, TINY(1.0_DP))
     projector_kpq_relative = projector_kpq_residual / &
          MAX(projector_kpq_scale, TINY(1.0_DP))
     coarse_relative = coarse_residual / &
          MAX(analytic_scale, TINY(1.0_DP))
     richardson_relative = richardson_residual / &
          MAX(analytic_scale, TINY(1.0_DP))
     IF (ionode) WRITE(stdout, &
          '(5x,"SOC_U_PW_BARE_RHS_ORACLE_CHECK",2(1x,i0),' // &
          '12(1x,es24.16))') ik, MERGE(1,0,time_reversed), &
          bare_oracle_eta, projector_k_scale, projector_k_residual, &
          projector_k_relative, projector_kpq_scale, &
          projector_kpq_residual, projector_kpq_relative, &
          analytic_scale, coarse_residual, coarse_relative, &
          richardson_residual, richardson_relative
  ENDIF
  !
  IF ((.NOT.lgamma) .AND. iverbosity == 1) THEN
     projection_scale = MAXVAL(ABS(projection))
     dprojection_scale = MAXVAL(ABS(dprojection))
     output_projector_scale = 0.0_DP
     doutput_projector_scale = 0.0_DP
     result_scale = 0.0_DP
     DO is = 1, npol
        row_first = 1 + npwx * (is - 1)
        output_projector_scale = MAX(output_projector_scale, &
             MAXVAL(ABS(output_projector( &
                  row_first:row_first+npwq-1,1:nwfcU))))
        doutput_projector_scale = MAX(doutput_projector_scale, &
             MAXVAL(ABS(doutput_projector( &
                  row_first:row_first+npwq-1,1:nwfcU))))
        result_scale = MAX(result_scale, MAXVAL(ABS(dvhubpsi( &
             row_first:row_first+npwq-1,1:nbnd))))
     ENDDO
     CALL mp_max(projection_scale, intra_bgrp_comm)
     CALL mp_max(dprojection_scale, intra_bgrp_comm)
     CALL mp_max(output_projector_scale, intra_bgrp_comm)
     CALL mp_max(doutput_projector_scale, intra_bgrp_comm)
     CALL mp_max(result_scale, intra_bgrp_comm)
     IF (ionode) WRITE(stdout, &
          '(5x,"SOC_U_Q_BARE_OPERATOR_CHECK",2(1x,i0),5(1x,es24.16))') &
          ik, MERGE(1,0,time_reversed), projection_scale, &
          dprojection_scale, output_projector_scale, &
          doutput_projector_scale, result_scale
  ENDIF
  !
  dvpsi(1:npwx*npol,1:nbnd) = dvpsi(1:npwx*npol,1:nbnd) + dvhubpsi
  !
  NULLIFY(output_projector)
  NULLIFY(doutput_projector)
  DEALLOCATE(dprojector_k)
  DEALLOCATE(dprojector_kpq)
  DEALLOCATE(dcolumn)
  DEALLOCATE(projection)
  DEALLOCATE(dprojection)
  DEALLOCATE(ground_matrix)
  DEALLOCATE(ground_matrix_tr)
  DEALLOCATE(bare_matrix)
  DEALLOCATE(bare_matrix_tr)
  DEALLOCATE(dns_pattern)
  DEALLOCATE(dvhubpsi)
  IF (oracle_active) THEN
     DEALLOCATE(dns_q_oracle)
     DEALLOCATE(dns_h0)
     DEALLOCATE(dns_h90)
     DEALLOCATE(ns_plus)
     DEALLOCATE(ns_minus)
     DEALLOCATE(v_plus)
     DEALLOCATE(v_minus)
     DEALLOCATE(v_ground_pw)
     DEALLOCATE(dv0)
     DEALLOCATE(dv90)
     DEALLOCATE(dvq_pw)
     DEALLOCATE(projector_k_fd)
     DEALLOCATE(projector_kpq_fd)
     DEALLOCATE(projection_fd)
     DEALLOCATE(action_plus)
     DEALLOCATE(action_minus)
     DEALLOCATE(action_half_plus)
     DEALLOCATE(action_half_minus)
     DEALLOCATE(derivative_coarse)
     DEALLOCATE(derivative_half)
     DEALLOCATE(derivative_richardson)
     DEALLOCATE(oracle_ground_matrix)
     DEALLOCATE(oracle_ground_matrix_tr)
     DEALLOCATE(oracle_dv_matrix)
     DEALLOCATE(oracle_dv_minus_q)
     DEALLOCATE(oracle_dv_matrix_tr)
     DEALLOCATE(oracle_matrix_used)
  ENDIF
  !
  CALL stop_clock('dvqhub_barepsi_nc')
  !
  RETURN
  !
  CONTAINS
  !
  SUBROUTINE report_bare_term_checksum(na_term, term_id)
    !---------------------------------------------------------------------
    !! Observe one already accumulated bare Hubbard term with two fixed
    !! complex-linear probes.  Only active k+q rows enter the checksum.
    !
    INTEGER, INTENT(IN) :: na_term, term_id
    INTEGER :: ibnd_term, ig_term, ig_global, is_out, row_first_term
    INTEGER :: hash_re, hash_im
    COMPLEX(DP) :: checksum_current(2,2), checksum_delta(2,2)
    COMPLEX(DP) :: probe
    !
    checksum_current = (0.0_DP, 0.0_DP)
    DO is_out = 1, npol
       row_first_term = 1 + npwx * (is_out - 1)
       DO ibnd_term = 1, nbnd_occ(ikk)
          DO ig_term = 1, npwq
             ig_global = igk_k(ig_term,ikq)
             hash_re = MODULO(17*mill(1,ig_global) + &
                  31*mill(2,ig_global) + 43*mill(3,ig_global) + &
                  47*ibnd_term + 11, 251) - 125
             hash_im = MODULO(59*mill(1,ig_global) + &
                  67*mill(2,ig_global) + 73*mill(3,ig_global) + &
                  79*ibnd_term + 13, 241) - 120
             probe = CMPLX(REAL(hash_re,DP)/128.0_DP, &
                           REAL(hash_im,DP)/128.0_DP, KIND=DP)
             checksum_current(1,is_out) = &
                  checksum_current(1,is_out) + probe * &
                  dvhubpsi(row_first_term+ig_term-1,ibnd_term)
             hash_re = MODULO(83*mill(1,ig_global) + &
                  89*mill(2,ig_global) + 97*mill(3,ig_global) + &
                  101*ibnd_term + 17, 239) - 119
             hash_im = MODULO(103*mill(1,ig_global) + &
                  107*mill(2,ig_global) + 109*mill(3,ig_global) + &
                  113*ibnd_term + 19, 233) - 116
             probe = CMPLX(REAL(hash_re,DP)/128.0_DP, &
                           REAL(hash_im,DP)/128.0_DP, KIND=DP)
             checksum_current(2,is_out) = &
                  checksum_current(2,is_out) + probe * &
                  dvhubpsi(row_first_term+ig_term-1,ibnd_term)
          ENDDO
       ENDDO
    ENDDO
    CALL mp_sum(checksum_current, intra_bgrp_comm)
    checksum_delta = checksum_current - term_checksum_previous
    IF (ionode) THEN
       DO is_out = 1, npol
          WRITE(stdout, &
               '(5x,"SOC_U_Q_BARE_TERM_CHECKSUM",6(1x,i0),' // &
               '4(1x,es24.16))') ik, ikk, MERGE(1,0,time_reversed), &
               na_term, term_id, is_out, &
               DBLE(checksum_delta(1,is_out)), &
               AIMAG(checksum_delta(1,is_out)), &
               DBLE(checksum_delta(2,is_out)), &
               AIMAG(checksum_delta(2,is_out))
       ENDDO
    ENDIF
    term_checksum_previous = checksum_current
    !
    RETURN
  END SUBROUTINE report_bare_term_checksum
  !
  SUBROUTINE build_bare_fd_action(amplitude, action)
    !
    REAL(DP), INTENT(IN) :: amplitude
    COMPLEX(DP), INTENT(OUT) :: action(npwx*npol,nbnd)
    INTEGER :: na_fd, nt_fd, ldim_fd, matrix_size_fd
    INTEGER :: ispin_fd, m_fd, iwfc_fd, ig_fd, icart_fd
    INTEGER :: row_first_fd, mu_fd
    REAL(DP) :: argument_fd
    COMPLEX(DP) :: phase_fd
    !
    projector_k_fd = swfcatomk
    projector_kpq_fd = output_projector
    DO na_fd = 1, nat
       nt_fd = ityp(na_fd)
       IF (is_hubbard(nt_fd)) THEN
          ldim_fd = 2 * Hubbard_l(nt_fd) + 1
          mu_fd = 3 * (na_fd - 1)
          DO ispin_fd = 1, npol
             row_first_fd = 1 + npwx * (ispin_fd - 1)
             DO m_fd = 1, ldim_fd
                iwfc_fd = offsetU(na_fd) + m_fd + &
                           ldim_fd * (ispin_fd - 1)
                DO ig_fd = 1, npw
                   argument_fd = 0.0_DP
                   DO icart_fd = 1, 3
                      argument_fd = argument_fd + &
                           (g(icart_fd,igk_k(ig_fd,ikk)) + &
                            xk(icart_fd,ikk)) * tpiba * &
                           DBLE(uact(mu_fd+icart_fd))
                   ENDDO
                   phase_fd = CMPLX(COS(amplitude * argument_fd), &
                                    -SIN(amplitude * argument_fd), &
                                    KIND=DP)
                   projector_k_fd(row_first_fd+ig_fd-1,iwfc_fd) = &
                        swfcatomk(row_first_fd+ig_fd-1,iwfc_fd) * phase_fd
                ENDDO
                DO ig_fd = 1, npwq
                   argument_fd = 0.0_DP
                   DO icart_fd = 1, 3
                      argument_fd = argument_fd + &
                           (g(icart_fd,igk_k(ig_fd,ikq)) + &
                            xk(icart_fd,ikq)) * tpiba * &
                           DBLE(uact(mu_fd+icart_fd))
                   ENDDO
                   phase_fd = CMPLX(COS(amplitude * argument_fd), &
                                    -SIN(amplitude * argument_fd), &
                                    KIND=DP)
                   projector_kpq_fd(row_first_fd+ig_fd-1,iwfc_fd) = &
                        output_projector(row_first_fd+ig_fd-1,iwfc_fd) * &
                        phase_fd
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    !
    CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
         (1.0_DP,0.0_DP), projector_k_fd, npwx*npol, evc, npwx*npol, &
         (0.0_DP,0.0_DP), projection_fd, nwfcU)
    CALL mp_sum(projection_fd, intra_bgrp_comm)
    action = (0.0_DP, 0.0_DP)
    DO na_fd = 1, nat
       nt_fd = ityp(na_fd)
       IF (is_hubbard(nt_fd)) THEN
          ldim_fd = 2 * Hubbard_l(nt_fd) + 1
          matrix_size_fd = ldim_fd * npol
          oracle_ground_matrix = (0.0_DP, 0.0_DP)
          oracle_ground_matrix_tr = (0.0_DP, 0.0_DP)
          oracle_dv_matrix = (0.0_DP, 0.0_DP)
          oracle_dv_minus_q = (0.0_DP, 0.0_DP)
          oracle_dv_matrix_tr = (0.0_DP, 0.0_DP)
          oracle_matrix_used = (0.0_DP, 0.0_DP)
          CALL pack_hubbard_blocks_nc(ldim_fd, npol, &
               v_ground_pw(:,:,:,na_fd), oracle_ground_matrix)
          CALL pack_hubbard_blocks_nc(ldim_fd, npol, &
               dvq_pw(:,:,:,na_fd), oracle_dv_matrix)
          IF (time_reversed) THEN
             CALL time_reverse_hubbard_matrix_nc(ldim_fd, npol, &
                  oracle_ground_matrix, oracle_ground_matrix_tr)
             oracle_dv_minus_q(1:matrix_size_fd,1:matrix_size_fd) = &
                  TRANSPOSE(CONJG(oracle_dv_matrix( &
                       1:matrix_size_fd,1:matrix_size_fd)))
             CALL time_reverse_hubbard_matrix_nc(ldim_fd, npol, &
                  oracle_dv_minus_q, oracle_dv_matrix_tr)
             oracle_matrix_used = oracle_ground_matrix_tr + &
                                  amplitude * oracle_dv_matrix_tr
          ELSE
             oracle_matrix_used = oracle_ground_matrix + &
                                  amplitude * oracle_dv_matrix
          ENDIF
          action = action + MATMUL( &
               projector_kpq_fd(:,offsetU(na_fd)+1: &
                    offsetU(na_fd)+matrix_size_fd), &
               MATMUL(oracle_matrix_used(1:matrix_size_fd, &
                                        1:matrix_size_fd), &
                    projection_fd(offsetU(na_fd)+1: &
                         offsetU(na_fd)+matrix_size_fd,:)))
       ENDIF
    ENDDO
    !
    RETURN
  END SUBROUTINE build_bare_fd_action
END SUBROUTINE dvqhub_barepsi_nc
