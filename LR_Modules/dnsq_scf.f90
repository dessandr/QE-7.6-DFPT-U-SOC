!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dnsq_scf (npe, lmetq0)
  !-----------------------------------------------------------------------
  !! DFPT+U: This routine calculates, for each SCF iteration, 
  !! the SCF variation of the occupation matrix ns, for npe perturbations.
  !! The result is stored in the variable dnsscf:
  !
  !! $$ \text{dnsscf}(m1,m2,\text{ispin},I,\text{ipert}) = 
  !!   \sum_{k,n} [ \langle \psi(n,k,\text{ispin})| S_{k}\phi(k,I,m1)\rangle  
  !!                \cdot \langle S_{k+q}\phi(k+q,I,m2)|
  !!                       d\psi(\text{ipert},n,k+q,\text{ispin})\rangle
  !!                + \langle\psi(n,k,\text{ispin})| S_{k}\phi(k,I,m2)\rangle  
  !!                \cdot \langle S_{k+q}\phi(k+q,I,m1)|
  !!                 d\psi(\text{ipert},n,k+q,\text{ispin})\rangle ] $$
  !
  !! Written  by A. Floris.  
  !! Modified by I. Timrov (01.10.2018).
  !
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : nwordwfcU
  USE units_lr,      ONLY : iuwfc, lrwfc, iudwf, lrdwf
  USE ions_base,     ONLY : nat, ityp
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_hubbard, offsetU, nwfcU, &
                            Hubbard_projectors
  USE ldaU_lr,       ONLY : swfcatomk, swfcatomkpq, dnsscf, lr_has_dnsorth, lr_dnsorth, &
                            proj1, proj2, accumulate_ns_nc_transition, &
                            ns_nc_hermiticity_residual, hub_dnsq_single_branch
  USE klist,         ONLY : wk, degauss, ngauss, ngk, degauss_cond
  USE wvfct,         ONLY : npwx, nbnd, et, nbnd_cond
  USE qpoint,        ONLY : nksq, ikks, ikqs, ikmks
  USE control_lr,    ONLY : lgamma, nbnd_occ
  USE units_lr,      ONLY : iuatswfc
  USE lsda_mod,      ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions, ONLY : evc
  USE ener,          ONLY : ef, ef_cond
  USE uspp,          ONLY : okvan 
  USE mp_pools,      ONLY : inter_pool_comm 
  USE mp_bands,      ONLY : intra_bgrp_comm   
  USE mp,            ONLY : mp_sum, mp_max
  USE io_global,     ONLY : stdout
  USE buffers,       ONLY : get_buffer
  USE efermi_shift,  ONLY : def
  USE lr_two_chem,   ONLY : def_val,def_cond
  USE two_chem,      ONLY : twochem
  USE control_flags, ONLY : iverbosity
  USE noncollin_module, ONLY : noncolin, domag, npol
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: npe 
  !! the number of perturbations
  LOGICAL,  INTENT(IN) :: lmetq0
  !! TRUE if \(xq=(0,0,0)\) in a metal
  !
  ! ... local variables
  !
  INTEGER  :: nt, ihubst, ihubst1, ihubst2, nah, m, m1, m2, ibnd, is, &
              ipert, nrec, ldim, npw, npwq, ik , ikk, ikq
  INTEGER  :: i, isolv, ikmk, ldim_atom
  INTEGER  :: nsolv_nc
  COMPLEX(DP), ALLOCATABLE :: proj_nc_left(:,:), proj_nc_right(:,:)
  COMPLEX(DP), ALLOCATABLE :: dnsscf_nc_branch(:,:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: replay_proj_left(:,:), replay_proj_right(:,:)
  COMPLEX(DP), ALLOCATABLE :: dnsscf_nc_replay_raw(:,:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dnsscf_nc_replay(:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: band_weight(:)
  REAL(DP) :: herm_res, branch_scale(2), branch_sum_res
  INTEGER :: is_bra, is_ket, iblock
  COMPLEX(DP) :: replay_transition
  REAL(DP) :: epsilon_spin(2,2)
  REAL(DP) :: replay_prod_norm, replay_ref_norm
  REAL(DP) :: replay_abs_res, replay_rel_res, replay_scale
  REAL(DP) :: replay_left_abs_max, replay_left_rel_max
  REAL(DP) :: replay_right_abs_max, replay_right_rel_max
  REAL(DP) :: replay_spin_coeff
  INTEGER :: iwfc1_replay, iwfc2_replay
  INTEGER :: out_bra, out_ket, out_block
  LOGICAL :: replay_active
  COMPLEX(DP), ALLOCATABLE :: dpsi(:,:)
  REAL(DP) :: weight, wdelta, w1
  REAL(DP), EXTERNAL :: w0gauss 
  !
  CALL start_clock( 'dnsq_scf' )
  ! 
  replay_active = .FALSE.
  ALLOCATE (dpsi(npwx*npol,nbnd))
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (proj2(nbnd,nwfcU))
  IF (noncolin) THEN
     ALLOCATE (proj_nc_left(nwfcU,nbnd))
     ALLOCATE (proj_nc_right(nwfcU,nbnd))
     ALLOCATE (band_weight(nbnd))
     ALLOCATE (dnsscf_nc_branch(2*Hubbard_lmax+1,2*Hubbard_lmax+1, &
                                nspin,nat,npe,2))
     replay_active = iverbosity == 1
     IF (replay_active) THEN
        ALLOCATE (replay_proj_left(nwfcU,nbnd))
        ALLOCATE (replay_proj_right(nwfcU,nbnd))
        ALLOCATE (dnsscf_nc_replay_raw(2*Hubbard_lmax+1, &
                                       2*Hubbard_lmax+1,nspin,nat,npe,2))
        ALLOCATE (dnsscf_nc_replay(2*Hubbard_lmax+1, &
                                   2*Hubbard_lmax+1,nspin,nat,npe,2))
        epsilon_spin = 0.0_DP
        epsilon_spin(1,2) = 1.0_DP
        epsilon_spin(2,1) = -1.0_DP
        replay_left_abs_max = 0.0_DP
        replay_left_rel_max = 0.0_DP
        replay_right_abs_max = 0.0_DP
        replay_right_rel_max = 0.0_DP
     ENDIF
  ENDIF
  !
  ldim = 2 * Hubbard_lmax + 1
  !   
  ! At each iteration dnsscf is set to zero, as it is recomputed
  ! 
  dnsscf = (0.d0, 0.d0) 
  IF (noncolin) dnsscf_nc_branch = (0.d0, 0.d0)
  IF (noncolin .AND. replay_active) THEN
     dnsscf_nc_replay_raw = (0.0_DP,0.0_DP)
     dnsscf_nc_replay = (0.0_DP,0.0_DP)
  ENDIF
  !
  IF (noncolin) THEN
     !
     ! MVP path: the magnetic noncollinear phonon response is
     ! reconstructed from the two Sternheimer records already owned by
     ! solve_linter.  The electric-field response owned by solve_e has one
     ! ordinary record; its adjoint transition is accumulated explicitly.
     !
     IF (.NOT.domag) CALL errore('dnsq_scf', &
                         'noncollinear Hubbard response requires a magnetic state', 1)
     IF (okvan) CALL errore('dnsq_scf', &
                    'noncollinear Hubbard response requires NCPP', 1)
     IF (Hubbard_projectors /= 'atomic') CALL errore('dnsq_scf', &
          'noncollinear Hubbard response requires atomic projectors', 1)
     nsolv_nc = 2
     IF (hub_dnsq_single_branch) nsolv_nc = 1
     !
     DO ik = 1, nksq
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq = ngk(ikq)
        !
        ! Both response branches are represented in the ordinary k+q G order:
        ! c_bands_nscf_ph applies apply_trev to the second ground-state
        ! record before storing it.
        !
        CALL get_buffer(swfcatomk, nwordwfcU, iuatswfc, ikk)
        IF (.NOT.lgamma) &
           CALL get_buffer(swfcatomkpq, nwordwfcU, iuatswfc, ikq)
        band_weight = 0.0_DP
        band_weight(1:nbnd_occ(ikk)) = wk(ikk)
        !
        DO isolv = 1, nsolv_nc
           IF (isolv == 1) THEN
              ikmk = ikk
           ELSE
              ikmk = ikmks(ik)
           ENDIF
           CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
           !
           CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                      (1.0_DP,0.0_DP), swfcatomk, npwx*npol, evc, &
                      npwx*npol, (0.0_DP,0.0_DP), proj_nc_left, nwfcU)
           CALL mp_sum(proj_nc_left, intra_bgrp_comm)
           IF (replay_active) THEN
              replay_proj_left = (0.0_DP,0.0_DP)
              DO ibnd = 1, nbnd_occ(ikk)
                 DO i = 1, nwfcU
                    replay_proj_left(i,ibnd) = DOT_PRODUCT( &
                         swfcatomk(1:npwx*npol,i), &
                         evc(1:npwx*npol,ibnd))
                 ENDDO
              ENDDO
              CALL mp_sum(replay_proj_left, intra_bgrp_comm)
              replay_abs_res = MAXVAL(ABS( &
                   proj_nc_left(:,1:nbnd_occ(ikk)) - &
                   replay_proj_left(:,1:nbnd_occ(ikk))))
              replay_scale = MAX( &
                   MAXVAL(ABS(proj_nc_left(:,1:nbnd_occ(ikk)))), &
                   MAXVAL(ABS(replay_proj_left(:,1:nbnd_occ(ikk)))))
              replay_left_abs_max = MAX(replay_left_abs_max, &
                                        replay_abs_res)
              replay_left_rel_max = MAX(replay_left_rel_max, &
                   replay_abs_res / MAX(replay_scale,TINY(1.0_DP)))
           ENDIF
           !
           DO ipert = 1, npe
              nrec = (isolv - 1) * npe * nksq + (ipert - 1) * nksq + ik
              CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
              CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                         (1.0_DP,0.0_DP), swfcatomkpq, npwx*npol, dpsi, &
                         npwx*npol, (0.0_DP,0.0_DP), proj_nc_right, nwfcU)
              CALL mp_sum(proj_nc_right, intra_bgrp_comm)
              IF (replay_active) THEN
                 replay_proj_right = (0.0_DP,0.0_DP)
                 DO ibnd = 1, nbnd_occ(ikk)
                    DO i = 1, nwfcU
                       replay_proj_right(i,ibnd) = DOT_PRODUCT( &
                            swfcatomkpq(1:npwx*npol,i), &
                            dpsi(1:npwx*npol,ibnd))
                    ENDDO
                 ENDDO
                 CALL mp_sum(replay_proj_right, intra_bgrp_comm)
                 replay_abs_res = MAXVAL(ABS( &
                      proj_nc_right(:,1:nbnd_occ(ikk)) - &
                      replay_proj_right(:,1:nbnd_occ(ikk))))
                 replay_scale = MAX( &
                      MAXVAL(ABS(proj_nc_right(:,1:nbnd_occ(ikk)))), &
                      MAXVAL(ABS(replay_proj_right(:,1:nbnd_occ(ikk)))))
                 replay_right_abs_max = MAX(replay_right_abs_max, &
                                            replay_abs_res)
                 replay_right_rel_max = MAX(replay_right_rel_max, &
                      replay_abs_res / MAX(replay_scale,TINY(1.0_DP)))
              ENDIF
              !
              DO nah = 1, nat
                 nt = ityp(nah)
                 IF (is_hubbard(nt)) THEN
                    ldim_atom = 2 * Hubbard_l(nt) + 1
                    CALL accumulate_ns_nc_transition(ldim_atom, offsetU(nah), &
                         npol, nbnd_occ(ikk), band_weight, proj_nc_left, &
                         proj_nc_right, &
                         dnsscf_nc_branch(:,:,:,nah,ipert,isolv), isolv == 2)
                    IF (hub_dnsq_single_branch) &
                       CALL accumulate_ns_nc_transition(ldim_atom, offsetU(nah), &
                            npol, nbnd_occ(ikk), band_weight, proj_nc_right, &
                            proj_nc_left, &
                            dnsscf_nc_branch(:,:,:,nah,ipert,isolv), .FALSE.)
                    IF (replay_active) THEN
                       DO is_bra = 1, npol
                          DO is_ket = 1, npol
                             iblock = npol * (is_bra - 1) + is_ket
                             DO m2 = 1, ldim_atom
                                iwfc2_replay = offsetU(nah) + m2 + &
                                     ldim_atom * (is_ket - 1)
                                DO m1 = 1, ldim_atom
                                   iwfc1_replay = offsetU(nah) + m1 + &
                                        ldim_atom * (is_bra - 1)
                                   replay_transition = (0.0_DP,0.0_DP)
                                   DO ibnd = 1, nbnd_occ(ikk)
                                      replay_transition = &
                                           replay_transition + &
                                           band_weight(ibnd) * CONJG( &
                                           replay_proj_left( &
                                           iwfc1_replay,ibnd)) * &
                                           replay_proj_right( &
                                           iwfc2_replay,ibnd)
                                   ENDDO
                                   dnsscf_nc_replay_raw( &
                                        m1,m2,iblock,nah,ipert,isolv) = &
                                        dnsscf_nc_replay_raw( &
                                        m1,m2,iblock,nah,ipert,isolv) + &
                                        replay_transition
                                   IF (hub_dnsq_single_branch) THEN
                                      replay_transition = (0.0_DP,0.0_DP)
                                      DO ibnd = 1, nbnd_occ(ikk)
                                         replay_transition = &
                                              replay_transition + &
                                              band_weight(ibnd) * CONJG( &
                                              replay_proj_right( &
                                              iwfc1_replay,ibnd)) * &
                                              replay_proj_left( &
                                              iwfc2_replay,ibnd)
                                      ENDDO
                                      dnsscf_nc_replay_raw( &
                                           m1,m2,iblock,nah,ipert,isolv) = &
                                           dnsscf_nc_replay_raw( &
                                           m1,m2,iblock,nah,ipert,isolv) + &
                                           replay_transition
                                   ENDIF
                                ENDDO
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
  ELSE
  !
  DO ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     !
     IF (lsda) current_spin = isk(ikk)
     !
     ! Read unperturbed KS wavefuctions psi(k)
     !
     IF (nksq.GT.1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
     !
     ! Read the atomic orbitals S*\phi at k and k+q from file (unit iuatswfc)
     !  
     CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
     IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
     !
     DO ipert = 1, npe
        !   
        nrec = (ipert - 1) * nksq + ik
        !
        ! At each SCF iteration for each ik and ipert read dpsi on iunit iudwf
        ! iudwf contains data for the actual perturbation ipert.
        !
        CALL get_buffer (dpsi, lrdwf, iudwf, nrec)
        !
        ! Calculate:
        ! proj1 (ibnd, ihubst) = < S_{k}\phi_(k,I,m) | psi(ibnd,k) >
        ! proj2 (ibnd, ihubst) = < S_{k+q}\phi_(k+q,I,m)) | dpsi(ipert,ibnd,k+q) >
        ! 
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           IF (is_hubbard(nt)) THEN
              !   
              DO m = 1, 2*Hubbard_l(nt)+1
                 !
                 ihubst = offsetU(nah) + m   ! I m index
                 !
                 DO ibnd = 1, nbnd_occ(ikk)
                    !
                    proj1(ibnd,ihubst) = dot_product (swfcatomk(1:npw,ihubst), evc(1:npw,ibnd))
                    proj2(ibnd,ihubst) = dot_product (swfcatomkpq(1:npwq,ihubst), dpsi(1:npwq,ibnd))
                    !
                 ENDDO
                 ! 
              ENDDO ! m
              !
           ENDIF
           !
        ENDDO
        !
        CALL mp_sum (proj1, intra_bgrp_comm)  
        CALL mp_sum (proj2, intra_bgrp_comm)
        !
        ! Calculate dnsscf. It is in the pattern basis, because dpsi is.                       
        ! The weights (theta functions) are already contained in dpsi, 
        ! according to Eq. (75) in Rev. Mod. Phys. 73, 515, (2001).
        !
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           IF (is_hubbard(nt)) THEN
              !  
              DO m1 = 1, 2*Hubbard_l(nt)+1
                 !
                 ihubst1 = offsetU(nah) + m1
                 !
                 DO m2 = m1, 2*Hubbard_l(nt)+1
                    !  
                    ihubst2 = offsetU(nah) + m2
                    !
                    DO ibnd = 1, nbnd_occ(ikk)
                       !    
                       dnsscf(m1,m2,current_spin,nah,ipert) = dnsscf(m1,m2,current_spin,nah,ipert) + &
                                        wk(ikk) * ( CONJG(proj1(ibnd,ihubst1)) * proj2(ibnd,ihubst2) &
                                                  + CONJG(proj1(ibnd,ihubst2)) * proj2(ibnd,ihubst1) )
                       !
                       ! Correction for metals at q=0 to be added when dpsi is NOT converged.
                       ! After the convergence the correction is built into 
                       ! dpsi by ef_fermi.f90 with flag==.true. (called in solve_linter)
                       ! wdelta is a smeared delta function at the Fermi level.  
                       !                       
                       IF (lmetq0) THEN
                          !
                          wdelta = w0gauss ( (ef-et(ibnd,ikk))/degauss, ngauss) / degauss
                          weight = wk(ikk)  
                          w1 = weight * wdelta   
                          !
                          IF (twochem) THEN
                           !two chem case, valence and condcution must be treated separately for the two fermi shifts
                             IF (ibnd.le.nbnd-nbnd_cond) then
                                  dnsscf(m1,m2,current_spin,nah,ipert) = &
                                  dnsscf(m1,m2,current_spin,nah,ipert) + w1 * def_val(ipert) * &
                                         CONJG(proj1(ibnd,ihubst1)) * proj1(ibnd,ihubst2)
                             ELSE
                                        wdelta = w0gauss ( (ef_cond-et(ibnd,ikk))/degauss_cond, ngauss) / degauss_cond
                                        dnsscf(m1,m2,current_spin,nah,ipert) = &
                                        dnsscf(m1,m2,current_spin,nah,ipert) + w1 * def_cond(ipert) * &
                                         CONJG(proj1(ibnd,ihubst1)) * proj1(ibnd,ihubst2)
                             !  
                             END IF
                          !
                          ELSE
                                          dnsscf(m1,m2,current_spin,nah,ipert) = &
                               dnsscf(m1,m2,current_spin,nah,ipert) + w1 * def(ipert) * &
                                         CONJG(proj1(ibnd,ihubst1)) * proj1(ibnd,ihubst2)
                          END IF
                          !
                       ENDIF
                       !
                    ENDDO
                    !
                 ENDDO
                 !
              ENDDO
              !
           ENDIF
           !
        ENDDO
        !
     ENDDO ! ipert
     !
  ENDDO ! ik 
  !
  ENDIF
  !
  ! The noncollinear response keeps all four complex spin blocks: no
  ! orbital-triangle fill and no nspin=1 factor.
  !
  IF (noncolin) THEN
     CALL mp_sum(dnsscf_nc_branch, inter_pool_comm)
     dnsscf = dnsscf_nc_branch(:,:,:,:,:,1) + &
              dnsscf_nc_branch(:,:,:,:,:,2)
     IF (replay_active) THEN
        CALL mp_sum(dnsscf_nc_replay_raw, inter_pool_comm)
        dnsscf_nc_replay(:,:,:,:,:,1) = &
             dnsscf_nc_replay_raw(:,:,:,:,:,1)
        !
        ! Reconstruct the partner without the production helper:
        ! epsilon * transpose(T) * transpose(epsilon), where the
        ! transpose exchanges the combined orbital-spin indices.
        !
        DO ipert = 1, npe
           DO nah = 1, nat
              nt = ityp(nah)
              IF (is_hubbard(nt)) THEN
                 ldim_atom = 2 * Hubbard_l(nt) + 1
                 DO out_bra = 1, npol
                    DO out_ket = 1, npol
                       out_block = npol * (out_bra - 1) + out_ket
                       DO is_bra = 1, npol
                          DO is_ket = 1, npol
                             iblock = npol * (is_bra - 1) + is_ket
                             replay_spin_coeff = &
                                  epsilon_spin(out_bra,is_ket) * &
                                  epsilon_spin(out_ket,is_bra)
                             IF (replay_spin_coeff /= 0.0_DP) THEN
                                DO m2 = 1, ldim_atom
                                   DO m1 = 1, ldim_atom
                                      dnsscf_nc_replay( &
                                           m2,m1,out_block,nah,ipert,2) = &
                                           dnsscf_nc_replay( &
                                           m2,m1,out_block,nah,ipert,2) + &
                                           replay_spin_coeff * &
                                           dnsscf_nc_replay_raw( &
                                           m1,m2,iblock,nah,ipert,2)
                                   ENDDO
                                ENDDO
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
        !
        CALL mp_max(replay_left_abs_max, inter_pool_comm)
        CALL mp_max(replay_left_rel_max, inter_pool_comm)
        CALL mp_max(replay_right_abs_max, inter_pool_comm)
        CALL mp_max(replay_right_rel_max, inter_pool_comm)
        WRITE(stdout,'(5x,"SOC_U_DNSSCF_REPLAY_BEGIN")')
        WRITE(stdout, &
             '(5x,"SOC_U_DNSSCF_REPLAY_PROJECTION_CHECK",' // &
             '4(1x,es24.16))') &
             replay_left_abs_max, replay_left_rel_max, &
             replay_right_abs_max, replay_right_rel_max
        DO isolv = 1, nsolv_nc
           DO ipert = 1, npe
              DO nah = 1, nat
                 nt = ityp(nah)
                 IF (is_hubbard(nt)) THEN
                    ldim_atom = 2 * Hubbard_l(nt) + 1
                    DO is_bra = 1, npol
                       DO is_ket = 1, npol
                          iblock = npol * (is_bra - 1) + is_ket
                          replay_prod_norm = SQRT(SUM(ABS( &
                               dnsscf_nc_branch(1:ldim_atom, &
                               1:ldim_atom,iblock,nah,ipert,isolv))**2))
                          replay_ref_norm = SQRT(SUM(ABS( &
                               dnsscf_nc_replay(1:ldim_atom, &
                               1:ldim_atom,iblock,nah,ipert,isolv))**2))
                          replay_abs_res = SQRT(SUM(ABS( &
                               dnsscf_nc_branch(1:ldim_atom, &
                               1:ldim_atom,iblock,nah,ipert,isolv) - &
                               dnsscf_nc_replay(1:ldim_atom, &
                               1:ldim_atom,iblock,nah,ipert,isolv))**2))
                          replay_rel_res = replay_abs_res / MAX( &
                               replay_prod_norm,replay_ref_norm, &
                               TINY(1.0_DP))
                          WRITE(stdout, &
                               '(5x,"SOC_U_DNSSCF_REPLAY_BLOCK_CHECK",' // &
                               '5(1x,i0),4(1x,es24.16))') &
                               isolv, ipert, nah, is_bra, is_ket, &
                               replay_prod_norm, replay_ref_norm, &
                               replay_abs_res, replay_rel_res
                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
           replay_prod_norm = SQRT(SUM(ABS( &
                dnsscf_nc_branch(:,:,:,:,:,isolv))**2))
           replay_ref_norm = SQRT(SUM(ABS( &
                dnsscf_nc_replay(:,:,:,:,:,isolv))**2))
           replay_abs_res = SQRT(SUM(ABS( &
                dnsscf_nc_branch(:,:,:,:,:,isolv) - &
                dnsscf_nc_replay(:,:,:,:,:,isolv))**2))
           replay_rel_res = replay_abs_res / MAX( &
                replay_prod_norm,replay_ref_norm,TINY(1.0_DP))
           WRITE(stdout, &
                '(5x,"SOC_U_DNSSCF_REPLAY_BRANCH_CHECK",' // &
                '1x,i0,4(1x,es24.16))') &
                isolv, replay_prod_norm, replay_ref_norm, &
                replay_abs_res, replay_rel_res
        ENDDO
        replay_prod_norm = SQRT(SUM(ABS(dnsscf)**2))
        replay_ref_norm = SQRT(SUM(ABS( &
             dnsscf_nc_replay(:,:,:,:,:,1) + &
             dnsscf_nc_replay(:,:,:,:,:,2))**2))
        replay_abs_res = SQRT(SUM(ABS( &
             dnsscf - dnsscf_nc_replay(:,:,:,:,:,1) - &
             dnsscf_nc_replay(:,:,:,:,:,2))**2))
        replay_rel_res = replay_abs_res / MAX( &
             replay_prod_norm,replay_ref_norm,TINY(1.0_DP))
        WRITE(stdout, &
             '(5x,"SOC_U_DNSSCF_REPLAY_SUM_CHECK",' // &
             '4(1x,es24.16))') replay_prod_norm, replay_ref_norm, &
             replay_abs_res, replay_rel_res
        WRITE(stdout,'(5x,"SOC_U_DNSSCF_REPLAY_END")')
     ENDIF
     !
     IF (iverbosity == 1) THEN
        branch_scale(1) = MAXVAL(ABS(dnsscf_nc_branch(:,:,:,:,:,1)))
        branch_scale(2) = MAXVAL(ABS(dnsscf_nc_branch(:,:,:,:,:,2)))
        branch_sum_res = MAXVAL(ABS(dnsscf - &
                         dnsscf_nc_branch(:,:,:,:,:,1) - &
                         dnsscf_nc_branch(:,:,:,:,:,2)))
        WRITE(stdout,'(5x,"SOC_U_DNSSCF_BRANCH_BEGIN")')
        DO isolv = 1, nsolv_nc
           WRITE(stdout,'(5x,"SOC_U_DNSSCF_BRANCH_SCALE",1x,i0,1x,es24.16)') &
                isolv, branch_scale(isolv)
           DO ipert = 1, npe
              DO nah = 1, nat
                 nt = ityp(nah)
                 IF (is_hubbard(nt)) THEN
                    ldim_atom = 2 * Hubbard_l(nt) + 1
                    DO is_bra = 1, npol
                       DO is_ket = 1, npol
                          iblock = npol * (is_bra - 1) + is_ket
                          DO m2 = 1, ldim_atom
                             DO m1 = 1, ldim_atom
                                WRITE(stdout, &
                                     '(5x,"SOC_U_DNSSCF_BRANCH_CHECK",' // &
                                     '7(1x,i0),2(1x,es24.16))') &
                                     isolv, ipert, nah, is_bra, is_ket, &
                                     m1, m2, &
                                     REAL(dnsscf_nc_branch( &
                                     m1,m2,iblock,nah,ipert,isolv),DP), &
                                     AIMAG(dnsscf_nc_branch( &
                                     m1,m2,iblock,nah,ipert,isolv))
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        WRITE(stdout, &
             '(5x,"SOC_U_DNSSCF_BRANCH_SUM_CHECK",1x,es24.16)') &
             branch_sum_res
        WRITE(stdout,'(5x,"SOC_U_DNSSCF_BRANCH_END")')
     ENDIF
  ELSE
  !
  CALL mp_sum (dnsscf, inter_pool_comm) 
  !
  ! Filling the m1 m2 lower triangular part 
  !
  DO m1 = 2, ldim
     DO m2 = 1, m1-1
        dnsscf(m1,m2,:,:,:) = dnsscf(m2,m1,:,:,:)
     ENDDO
  ENDDO
  !
  ! In nspin.eq.1 k-point weight wk is normalized to 2 el/band 
  ! in the whole BZ but we are interested in dns of one spin component
  !
  IF (nspin.EQ.1) dnsscf = 0.5d0 * dnsscf
  !
  ENDIF
  !
  ! USPP case: add to the dnsscf calculated with the P_c^+dpsi 
  ! the non-scf part coming from the orthogonality constraints. 
  ! The orthogonality correction comes only for an atomic displacement. 
  ! In the case of the electric field calculation, lr_has_dnsorth is false.
  !
  IF (okvan .AND. lr_has_dnsorth) THEN
     dnsscf = dnsscf + lr_dnsorth
  ENDIF
  !
  ! Symmetrize dnsscf
  !
  CALL sym_dns(ldim, npe, dnsscf)
  !
  ! Write symmetrized dnsscf in the pattern basis 
  ! to the standard output
  !
  IF (iverbosity==1) THEN
     WRITE(stdout,*) 'DNSSCF SYMMETRIZED IN THE PATTERN BASIS'
     DO ipert = 1, npe
        WRITE(stdout,'(a,1x,i2)') ' partner # ', ipert
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO is = 1, nspin
                 WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                 DO m1 = 1, ldim
                    WRITE(stdout,'(10(f15.10,1x))') dnsscf(m1,:,is,nah,ipert)
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO     
     IF (noncolin .AND. lgamma) THEN
        herm_res = 0.0_DP
        DO ipert = 1, npe
           DO nah = 1, nat
              nt = ityp(nah)
              IF (is_hubbard(nt)) THEN
                 ldim_atom = 2 * Hubbard_l(nt) + 1
                 herm_res = MAX(herm_res, ns_nc_hermiticity_residual( &
                            ldim_atom, dnsscf(:,:,:,nah,ipert)))
              ENDIF
           ENDDO
        ENDDO
        WRITE(stdout,'(5x,"dnsscf nc Hermiticity residual =",1x,es12.5)') herm_res
     ENDIF
  ENDIF
  !
  IF (noncolin) THEN
     DEALLOCATE (proj_nc_left)
     DEALLOCATE (proj_nc_right)
     DEALLOCATE (band_weight)
     DEALLOCATE (dnsscf_nc_branch)
     IF (replay_active) THEN
        DEALLOCATE (replay_proj_left)
        DEALLOCATE (replay_proj_right)
        DEALLOCATE (dnsscf_nc_replay_raw)
        DEALLOCATE (dnsscf_nc_replay)
     ENDIF
  ENDIF
  DEALLOCATE (proj1)
  DEALLOCATE (proj2)
  DEALLOCATE (dpsi)
  !
  CALL stop_clock( 'dnsq_scf' )
  ! 
  RETURN
  ! 
END SUBROUTINE dnsq_scf
!----------------------------------------------------------------------------
