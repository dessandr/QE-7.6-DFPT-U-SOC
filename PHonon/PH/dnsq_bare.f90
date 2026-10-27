!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dnsq_bare
  !-----------------------------------------------------------------------
  !! DFPT+U : This routine calculates dnsbare, i.e. the bare variation
  !! of the occupation matrix ns under the perturbation q, due to a
  !! displacement of the atom L in the direction icart.
  !
  !! $$ \text{dnsbare}(m1,m2,\text{ispin},I,\text{icart},L) =
  !!    \sum_{k,n} [ \langle\psi(n,k,\text{ispin})| S_{k}\phi_(k,I,m1)\rangle \cdot
  !!                 \langle\Delta^{L \text{icart}}_{-q}(S_{k+q}\phi(k+q,I,m2))
  !!                  |\psi(n,k,\text{ispin})\rangle
  !!               + \langle \psi(n,k,\text{ispin})| S_{k}\phi(k,I,m2)> \cdot
  !!                 \langle\Delta^{L \text{icart}}_{-q}(S_{k+q}\phi(k+q,I,m1))
  !!                  |\psi(n,k,\text{ispin})\rangle ] $$
  !
  !! Written  by A. Floris
  !! Modified by I. Timrov (01.10.2018)
  !
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : nwordwfcU
  USE units_lr,      ONLY : iuwfc, lrwfc
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE klist,         ONLY : xk, ngk, igk_k
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, nwfcU, &
                            Hubbard_projectors, wfcU
  USE ldaU_ph,       ONLY : wfcatomk, wfcatomkpq, dwfcatomkpq, &
                            sdwfcatomk, sdwfcatomkpq, dvkb, vkbkpq, dvkbkpq, &
                            dnsbare, dnsbare_all_modes, read_dns_bare
  USE ldaU_lr,       ONLY : swfcatomk, swfcatomkpq, proj1, proj2, &
                            accumulate_ns_nc_transition, ns_nc_hermiticity_residual
  USE wvfct,         ONLY : npwx, wg, nbnd
  USE uspp,          ONLY : vkb, nkb, ofsbeta, okvan
  USE qpoint,        ONLY : nksq, ikks, ikqs, xq, ikmks
  USE control_lr,    ONLY : lgamma
  USE units_lr,      ONLY : iuatwfc, iuatswfc
  USE uspp_param,    ONLY : nh, nhm
  USE lsda_mod,      ONLY : lsda, nspin, current_spin, isk
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp_pools,      ONLY : inter_pool_comm, me_pool, nproc_pool
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum, mp_bcast
  USE mp_world,      ONLY : world_comm
  USE io_files,      ONLY : seqopn
  USE buffers,       ONLY : get_buffer
  USE control_flags, ONLY : iverbosity
  USE wavefunctions, ONLY : evc
  USE noncollin_module, ONLY : noncolin, npol, domag
  USE scf,            ONLY : rho
  USE mp_images,     ONLY : intra_image_comm
  USE uspp_init,        ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: i, j, k, icart, ios, is, na, nt, l, ih,        &
             ibeta, n, ihubst, ihubst1, ihubst2, nah, m, m1, m2,  &
             ibnd, ldim, npw, npwq, ik, ikk, ikq, iundnsbare, is_bra, is_ket, &
             ldim_atom, iblock, jblock, nb_s, nb_e, mykey
  COMPLEX(DP), ALLOCATABLE :: dqsphi(:,:), dmqsphi(:,:), dwfcatom_(:)
  COMPLEX(DP), ALLOCATABLE :: proj_nc_left(:,:), proj_nc_right(:,:)
  COMPLEX(DP), ALLOCATABLE :: proj_nc_partner_left(:,:), &
                              proj_nc_partner_right(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_partner(:,:)
  COMPLEX(DP), ALLOCATABLE :: ns_nc_ground(:,:,:,:), ns_nc_ground_partner(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dnsbare_nc_branch_ordinary(:,:,:,:,:,:), &
                              dnsbare_nc_branch_partner(:,:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dnsbare_pw_kernel(:,:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dns_nc_pw(:,:,:,:), wfcu_saved(:,:)
  COMPLEX(DP) :: ground_average_value
  REAL(DP), ALLOCATABLE :: band_weight(:)
  REAL(DP) :: ground_res, ground_partner_res, ground_average_res
  REAL(DP) :: ground_scale_ph, ground_partner_scale, ground_average_scale
  REAL(DP) :: ground_scale_pw, herm_res
  REAL(DP) :: pw_kernel_scale_ph, pw_kernel_scale_pw
  REAL(DP) :: pw_kernel_res
  REAL(DP) :: q_branch_scale(2), q_branch_sum_res
  LOGICAL :: exst, do_ground_check, do_partner_check, do_pw_kernel_check
  LOGICAL :: use_q_partner, need_partner_record
  !
  CALL start_clock( 'dnsq_bare' )
  !
  ios = 0
  ldim = 2 * Hubbard_lmax + 1
  do_ground_check = noncolin .AND. (iverbosity == 1)
  do_partner_check = do_ground_check .AND. domag
  do_pw_kernel_check = noncolin .AND. lgamma .AND. (iverbosity == 1)
  use_q_partner = noncolin .AND. (.NOT.lgamma)
  need_partner_record = do_partner_check .OR. use_q_partner
  !
  ALLOCATE (dqsphi(npwx*npol,nwfcU))
  ALLOCATE (dmqsphi(npwx*npol,nwfcU))
  ALLOCATE (dwfcatom_(npwx*npol))
  IF (noncolin) THEN
     ALLOCATE (proj_nc_left(nwfcU,nbnd))
     ALLOCATE (proj_nc_right(nwfcU,nbnd))
     ALLOCATE (band_weight(nbnd))
     IF (do_ground_check) &
        ALLOCATE (ns_nc_ground(ldim,ldim,nspin,nat))
     IF (need_partner_record) THEN
        ALLOCATE (evc_partner(npwx*npol,nbnd))
     ENDIF
     IF (do_partner_check) THEN
        ALLOCATE (ns_nc_ground_partner(ldim,ldim,nspin,nat))
     ENDIF
     IF (use_q_partner) THEN
        ALLOCATE (proj_nc_partner_left(nwfcU,nbnd))
        ALLOCATE (proj_nc_partner_right(nwfcU,nbnd))
        ALLOCATE (dnsbare_nc_branch_ordinary(ldim,ldim,nspin,nat,3,nat))
        ALLOCATE (dnsbare_nc_branch_partner(ldim,ldim,nspin,nat,3,nat))
     ENDIF
     IF (do_pw_kernel_check) THEN
        ALLOCATE (dnsbare_pw_kernel(ldim,ldim,nspin,nat,3,nat))
        ALLOCATE (dns_nc_pw(ldim,ldim,nspin,nat))
        ALLOCATE (wfcu_saved(npwx*npol,nwfcU))
     ENDIF
  ELSE
     ALLOCATE (proj1(nbnd,nwfcU))
     ALLOCATE (proj2(nbnd,nwfcU))
  ENDIF
  !
  IF (noncolin) THEN
     proj_nc_left  = (0.d0, 0.d0)
     proj_nc_right = (0.d0, 0.d0)
     band_weight   = 0.d0
     IF (do_ground_check) ns_nc_ground = (0.d0, 0.d0)
     IF (need_partner_record) evc_partner = (0.d0, 0.d0)
     IF (do_partner_check) THEN
        ns_nc_ground_partner = (0.d0, 0.d0)
     ENDIF
     IF (use_q_partner) THEN
        proj_nc_partner_left = (0.d0, 0.d0)
        proj_nc_partner_right = (0.d0, 0.d0)
        dnsbare_nc_branch_ordinary = (0.d0, 0.d0)
        dnsbare_nc_branch_partner = (0.d0, 0.d0)
     ENDIF
     IF (do_pw_kernel_check) THEN
        dnsbare_pw_kernel = (0.d0, 0.d0)
        dns_nc_pw = (0.d0, 0.d0)
        wfcu_saved = (0.d0, 0.d0)
     ENDIF
  ELSE
     proj1 = (0.d0, 0.d0)
     proj2 = (0.d0, 0.d0)
  ENDIF
  dnsbare = (0.d0, 0.d0)
  !
  ! Note: sdwfcatomkpq is not needed in this routine when we call delta_sphi,
  ! so we set it to zero and do not calculate.
  !
  sdwfcatomk   = (0.d0, 0.d0)
  sdwfcatomkpq = (0.d0, 0.d0)
  !
  exst = .FALSE.
  !
  ! The unit number
  iundnsbare = 37
  !
  IF (ionode) CALL seqopn (iundnsbare, 'dnsbare', 'formatted', exst)
  !
  ! Read (unsymmetrized) dnsbare from file (if it was already computed)
  !
  IF (read_dns_bare) THEN
     !
     WRITE(stdout,*) 'Does the file dnsbare exist in tmp_dir and is it read correctly?'
     !
     IF (ionode) THEN
        IF (exst) THEN
           READ(iundnsbare,*,iostat=ios) dnsbare
           IF (ios==0) THEN
              WRITE(stdout,*) '...yes. the dnsbare matrix was read correctly from file'
           ELSE
              WRITE(stdout,*) '...no. the file exists but it seems to be corrupted'
           ENDIF
        ELSE
           WRITE(stdout,*) '...no. the file does not exist. '
        ENDIF
     ENDIF
     !
     CALL mp_bcast(exst, ionode_id, world_comm)
     CALL mp_bcast(ios,  ionode_id, world_comm)
     IF (exst .AND. ios==0) CALL mp_bcast(dnsbare, ionode_id, world_comm)
     !
     ! IT: Is it needed to broadcast for intra_image_comm?
     CALL mp_bcast(ios,  ionode_id, intra_image_comm)
     IF (exst .AND. ios==0) CALL mp_bcast(dnsbare, ionode_id, intra_image_comm)
     !
  ENDIF
  !
  ! Compute dnsbare (if it was not already done)
  !
  IF ((.NOT.exst) .OR. (ios/=0) .OR. (.NOT.read_dns_bare)) THEN
     !
     WRITE(stdout,'(/5x,"Calculating the dnsbare matrix...")')
     dnsbare = (0.d0, 0.d0)
     IF (noncolin) THEN
        IF ((.NOT.lgamma) .AND. (.NOT.domag)) CALL errore('dnsq_bare', &
             'non-Gamma noncollinear Hubbard response requires a magnetic state', 1)
        IF (okvan) CALL errore('dnsq_bare', &
                       'noncollinear Hubbard response requires NCPP', 1)
        IF (Hubbard_projectors /= 'atomic') CALL errore('dnsq_bare', &
             'noncollinear Hubbard response requires atomic projectors', 1)
     ENDIF
     IF (do_pw_kernel_check) THEN
        wfcu_saved = wfcU
        CALL block_distribute(nbnd, me_pool, nproc_pool, nb_s, nb_e, mykey)
     ENDIF
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
        ! Compute beta functions at k (vkb) and at k+q (vkbkpq)
        !
        CALL init_us_2 (npw, igk_k(1,ikk), xk(:,ikk), vkb)
        IF (.NOT.lgamma) CALL init_us_2 (npwq, igk_k(1,ikq), xk(:,ikq), vkbkpq)
        !
        ! Compute the derivatives of beta functions d^{icart}beta at k and k+q
        ! for all bands and for the 3 cartesian directions
        !
        DO icart = 1, 3
           DO na = 1, nat
              nt = ityp(na)
              DO ih = 1, nh(nt)
                 !
                 ibeta = ofsbeta(na) + ih
                 !
                 CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                            vkb(:,ibeta), dvkb(:,ibeta,icart))
                 !
                 IF (.NOT.lgamma) &
                 CALL dwfc (npwq, igk_k(1,ikq), ikq, icart, &
                            vkbkpq(:,ibeta), dvkbkpq(:,ibeta,icart))
                 !
              ENDDO
           ENDDO
        ENDDO
        !
        ! Read the atomic orbitals \phi at k and k+q from file (unit iuatwfc)
        !
        CALL get_buffer (wfcatomk, nwordwfcU, iuatwfc, ikk)
        IF (.NOT.lgamma) CALL get_buffer (wfcatomkpq, nwordwfcU, iuatwfc, ikq)
        !
        ! Read S*\phi at k from file (unit iuatswfc)
        !
        CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
        IF (do_pw_kernel_check) wfcU = wfcatomk
        !
        IF (noncolin) THEN
           !
           ! The ground-state projection does not depend on the displaced
           ! atom or Cartesian direction.  Build it once per k point.
           !
           CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                      (1.0_DP,0.0_DP), swfcatomk, npwx*npol, evc, &
                      npwx*npol, (0.0_DP,0.0_DP), proj_nc_left, nwfcU)
           CALL mp_sum(proj_nc_left, intra_bgrp_comm)
           band_weight(1:nbnd) = wg(1:nbnd,ikk)
           !
           IF (do_ground_check) THEN
              !
              ! N_ab = sum_kn w_kn * CONJG(P_a,kn) * P_b,kn.
              !
              DO nah = 1, nat
                 nt = ityp(nah)
                 IF (is_hubbard(nt)) THEN
                    ldim_atom = 2 * Hubbard_l(nt) + 1
                    CALL accumulate_ns_nc_transition(ldim_atom, offsetU(nah), &
                         npol, nbnd, band_weight, proj_nc_left, proj_nc_left, &
                         ns_nc_ground(:,:,:,nah), .FALSE.)
                 ENDIF
              ENDDO
           ENDIF
           !
           IF (do_partner_check) THEN
              !
              ! c_bands_nscf_ph has already applied time reversal to this
              ! record and rewritten it in the odd-record G-vector order.
              !
              CALL get_buffer(evc_partner, lrwfc, iuwfc, ikmks(ik))
              CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                         (1.0_DP,0.0_DP), swfcatomk, npwx*npol, evc_partner, &
                         npwx*npol, (0.0_DP,0.0_DP), proj_nc_right, nwfcU)
              CALL mp_sum(proj_nc_right, intra_bgrp_comm)
              DO nah = 1, nat
                 nt = ityp(nah)
                 IF (is_hubbard(nt)) THEN
                    ldim_atom = 2 * Hubbard_l(nt) + 1
                    CALL accumulate_ns_nc_transition(ldim_atom, offsetU(nah), &
                         npol, nbnd, band_weight, proj_nc_right, proj_nc_right, &
                         ns_nc_ground_partner(:,:,:,nah), .TRUE.)
                 ENDIF
              ENDDO
           ENDIF
           !
           IF (use_q_partner) THEN
              !
              ! For a non-Gamma magnetic response, the second bare
              ! transition is supplied by the existing time-reversed
              ! (-k,-k-q) ground-state branch.  c_bands_nscf_ph has
              ! already rewritten this record in the ordinary-k G order.
              !
              IF (do_partner_check) THEN
                 proj_nc_partner_left = proj_nc_right
              ELSE
                 CALL get_buffer(evc_partner, lrwfc, iuwfc, ikmks(ik))
                 CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                            (1.0_DP,0.0_DP), swfcatomk, npwx*npol, &
                            evc_partner, npwx*npol, (0.0_DP,0.0_DP), &
                            proj_nc_partner_left, nwfcU)
                 CALL mp_sum(proj_nc_partner_left, intra_bgrp_comm)
              ENDIF
           ENDIF
        ENDIF
        !
        DO na = 1, nat
           !
           DO icart = 1, 3
              !
              dqsphi  = (0.d0, 0.d0)
              dmqsphi = (0.d0, 0.d0)
              !
              IF (noncolin) THEN
                 !
                 ! MVP path: Gamma, norm-conserving, atomic projectors.
                 ! Keep all spinor rows and build dN in the new_ns_nc
                 ! orbital-spin convention.
                 !
                 DO nah = 1, nat
                    nt = ityp(nah)
                    IF (is_hubbard(nt)) THEN
                       ldim_atom = 2 * Hubbard_l(nt) + 1
                       DO is_bra = 1, npol
                          DO m = 1, ldim_atom
                             ihubst = offsetU(nah) + m + &
                                       ldim_atom * (is_bra - 1)
                             IF (nah == na) THEN
                                dwfcatom_ = (0.d0, 0.d0)
                                DO is = 1, npol
                                   CALL dwfc(npw, igk_k(1,ikk), ikk, icart, &
                                        wfcatomk(1+npwx*(is-1),ihubst), &
                                        dwfcatom_(1+npwx*(is-1)))
                                ENDDO
                                dmqsphi(:,ihubst) = dwfcatom_
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDIF
                 ENDDO
                 !
                 CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                            (1.0_DP,0.0_DP), dmqsphi, npwx*npol, evc, &
                            npwx*npol, (0.0_DP,0.0_DP), proj_nc_right, nwfcU)
                 CALL mp_sum(proj_nc_right, intra_bgrp_comm)
                 IF (use_q_partner) THEN
                    CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                         (1.0_DP,0.0_DP), dmqsphi, npwx*npol, evc_partner, &
                         npwx*npol, (0.0_DP,0.0_DP), &
                         proj_nc_partner_right, nwfcU)
                    CALL mp_sum(proj_nc_partner_right, intra_bgrp_comm)
                 ENDIF
                 !
                 DO nah = 1, nat
                    nt = ityp(nah)
                    IF (is_hubbard(nt)) THEN
                       ldim_atom = 2 * Hubbard_l(nt) + 1
                       IF (use_q_partner) THEN
                          CALL accumulate_ns_nc_transition(ldim_atom, &
                               offsetU(nah), npol, nbnd, band_weight, &
                               proj_nc_left, proj_nc_right, &
                               dnsbare_nc_branch_ordinary( &
                                    :,:,:,nah,icart,na), &
                               .FALSE.)
                          CALL accumulate_ns_nc_transition(ldim_atom, &
                               offsetU(nah), npol, nbnd, band_weight, &
                               proj_nc_partner_left, &
                               proj_nc_partner_right, &
                               dnsbare_nc_branch_partner( &
                                    :,:,:,nah,icart,na), &
                               .TRUE.)
                       ELSE
                          CALL accumulate_ns_nc_transition(ldim_atom, &
                               offsetU(nah), npol, nbnd, band_weight, &
                               proj_nc_left, proj_nc_right, &
                               dnsbare(:,:,:,nah,icart,na), .FALSE.)
                          CALL accumulate_ns_nc_transition(ldim_atom, &
                               offsetU(nah), npol, nbnd, band_weight, &
                               proj_nc_right, proj_nc_left, &
                               dnsbare(:,:,:,nah,icart,na), .FALSE.)
                       ENDIF
                    ENDIF
                 ENDDO
                 !
                 IF (do_pw_kernel_check) THEN
                    !
                    ! Reuse the PW force owner on the identical PH state.
                    ! dndtau_k_nc returns the transposed force convention;
                    ! map it back to the new_ns_nc orbital-spin layout.
                    !
                    CALL dndtau_k_nc(ldim, proj_nc_left, evc, na, &
                         ofsbeta(na), icart, ikk, nb_s, nb_e, mykey, 1, &
                         dns_nc_pw)
                    DO nah = 1, nat
                       nt = ityp(nah)
                       IF (is_hubbard(nt)) THEN
                          ldim_atom = 2 * Hubbard_l(nt) + 1
                          DO is_bra = 1, npol
                             DO is_ket = 1, npol
                                iblock = npol * (is_bra - 1) + is_ket
                                jblock = npol * (is_ket - 1) + is_bra
                                dnsbare_pw_kernel(1:ldim_atom,1:ldim_atom, &
                                     iblock,nah,icart,na) = &
                                     dnsbare_pw_kernel(1:ldim_atom,1:ldim_atom, &
                                     iblock,nah,icart,na) + TRANSPOSE( &
                                     dns_nc_pw(1:ldim_atom,1:ldim_atom, &
                                     jblock,nah))
                             ENDDO
                          ENDDO
                       ENDIF
                    ENDDO
                    !
                 ENDIF
                 !
              ELSE
                 !
                 DO nah = 1, nat
                    !
                    nt = ityp(nah)
                    !
                    IF (is_hubbard(nt)) THEN
                       !
                       DO m = 1, 2*Hubbard_l(nt)+1
                          !
                          ihubst = offsetU(nah) + m ! I m index
                          !
                          IF (nah==na) THEN
                             !
                             ! Calculate |d_icart\phi_(k,I,m)) >
                             !
                             CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                                        wfcatomk(:,ihubst), dwfcatom_)
                             !
                             ! Apply the S matrix: | S_{k} d_^(I,icart)\phi_(k,I,m) >
                             !
                             CALL swfc (npw, 1, vkb, dwfcatom_, sdwfcatomk(:,ihubst))
                             !
                          ENDIF
                          !
                          ! Calculate dmqsphi = |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) >
                          !
                          CALL delta_sphi (ikk, ikq, na, icart, nah, ihubst, wfcatomk, wfcatomkpq, &
                                           sdwfcatomk, sdwfcatomkpq, vkb, vkbkpq, dvkb(:,:,icart), &
                                           dvkbkpq(:,:,icart), dqsphi, dmqsphi, 0)
                          !
                          ! Calculate:
                          ! proj1 (ibnd,ihubst) = < S_{k}\phi_(k,I,m) | psi(ibnd,k) >
                          ! proj2 (ibnd,ihubst) = < \Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) | psi(ibnd,k)>
                          !
                          !FIXME these  are  two zgemms  that could  be done outside the loop 
                          DO ibnd = 1, nbnd
                            proj1 (ibnd,ihubst) = dot_product (swfcatomk(1:npw,ihubst), evc(1:npw,ibnd))
                            proj2 (ibnd,ihubst) = dot_product ( dmqsphi(1:npw ,ihubst ),evc(1:npw,ibnd))
                          ENDDO
                          !
                       ENDDO
                       !
                    ENDIF
                    !
                 ENDDO
                 !
                 CALL mp_sum (proj1, intra_bgrp_comm)
                 CALL mp_sum (proj2, intra_bgrp_comm)
                 !
                 ! Finally calculate dnsbare.
                 ! dnsbare is symmetric: dnsbare(m1,m2) = dnsbare(m2,m1)
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
                             DO ibnd = 1, nbnd
                                dnsbare (m1,m2,current_spin,nah,icart,na) = &
                                dnsbare (m1,m2,current_spin,nah,icart,na) + wg(ibnd,ikk) * &
                                     ( CONJG(proj1(ibnd,ihubst1)) * proj2(ibnd,ihubst2)  + &
                                       CONJG(proj1(ibnd,ihubst2)) * proj2(ibnd,ihubst1) )
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
              ENDIF
              !
           ENDDO ! icart
           !
        ENDDO ! na
        !
     ENDDO ! ik
     !
     IF (do_pw_kernel_check) wfcU = wfcu_saved
     IF (use_q_partner) THEN
        CALL mp_sum(dnsbare_nc_branch_ordinary, inter_pool_comm)
        CALL mp_sum(dnsbare_nc_branch_partner, inter_pool_comm)
        dnsbare = dnsbare_nc_branch_ordinary + dnsbare_nc_branch_partner
     ELSE
        CALL mp_sum(dnsbare, inter_pool_comm)
     ENDIF
     IF (use_q_partner .AND. iverbosity == 1) THEN
        q_branch_scale(1) = MAXVAL(ABS(dnsbare_nc_branch_ordinary))
        q_branch_scale(2) = MAXVAL(ABS(dnsbare_nc_branch_partner))
        q_branch_sum_res = MAXVAL(ABS(dnsbare - &
             dnsbare_nc_branch_ordinary - dnsbare_nc_branch_partner))
        IF (ionode) THEN
           WRITE(stdout,'(5x,"SOC_U_Q_BARE_BEGIN",3(1x,es24.16))') xq
           DO is = 1, 2
              WRITE(stdout, &
                   '(5x,"SOC_U_Q_BARE_BRANCH_SCALE",1x,i0,1x,es24.16)') &
                   is, q_branch_scale(is)
           ENDDO
           WRITE(stdout, &
                '(5x,"SOC_U_Q_BARE_BRANCH_SUM_CHECK",1x,es24.16)') &
                q_branch_sum_res
           DO na = 1, nat
              DO icart = 1, 3
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
                                        '(5x,"SOC_U_Q_BARE_COMPONENT",' // &
                                        '7(1x,i0),2(1x,es24.16))') &
                                        na, icart, nah, is_bra, is_ket, &
                                        m1, m2, &
                                        REAL(dnsbare(m1,m2,iblock,nah, &
                                             icart,na),DP), &
                                        AIMAG(dnsbare(m1,m2,iblock,nah, &
                                             icart,na))
                                ENDDO
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
           WRITE(stdout,'(5x,"SOC_U_Q_BARE_END")')
        ENDIF
     ENDIF
     IF (do_pw_kernel_check) THEN
        CALL mp_sum(dnsbare_pw_kernel, inter_pool_comm)
        pw_kernel_scale_ph = 0.0_DP
        pw_kernel_scale_pw = 0.0_DP
        pw_kernel_res = 0.0_DP
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              ldim_atom = 2 * Hubbard_l(nt) + 1
              pw_kernel_scale_ph = MAX(pw_kernel_scale_ph, MAXVAL(ABS( &
                   dnsbare(1:ldim_atom,1:ldim_atom,1:nspin,nah,:,:))))
              pw_kernel_scale_pw = MAX(pw_kernel_scale_pw, MAXVAL(ABS( &
                   dnsbare_pw_kernel(1:ldim_atom,1:ldim_atom,1:nspin, &
                   nah,:,:))))
              pw_kernel_res = MAX(pw_kernel_res, MAXVAL(ABS( &
                   dnsbare(1:ldim_atom,1:ldim_atom,1:nspin,nah,:,:) - &
                   dnsbare_pw_kernel(1:ldim_atom,1:ldim_atom,1:nspin, &
                   nah,:,:))))
           ENDIF
        ENDDO
        IF (ionode) WRITE(stdout, &
             '(5x,"SOC_U_BARE_PW_KERNEL_CHECK",4(1x,es24.16))') &
             pw_kernel_scale_ph, pw_kernel_scale_pw, pw_kernel_res, &
             pw_kernel_res / MAX(pw_kernel_scale_pw, TINY(1.0_DP))
     ENDIF
     IF (do_ground_check) THEN
        CALL mp_sum(ns_nc_ground, inter_pool_comm)
        IF (do_partner_check) &
           CALL mp_sum(ns_nc_ground_partner, inter_pool_comm)
        ground_scale_ph = 0.0_DP
        ground_partner_scale = 0.0_DP
        ground_average_scale = 0.0_DP
        ground_scale_pw = 0.0_DP
        ground_res = 0.0_DP
        ground_partner_res = 0.0_DP
        ground_average_res = 0.0_DP
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              ldim_atom = 2 * Hubbard_l(nt) + 1
              DO is = 1, nspin
                 DO m2 = 1, ldim_atom
                    DO m1 = 1, ldim_atom
                       ground_scale_ph = MAX(ground_scale_ph, &
                            ABS(ns_nc_ground(m1,m2,is,nah)))
                       ground_scale_pw = MAX(ground_scale_pw, &
                            ABS(rho%ns_nc(m1,m2,is,nah)))
                       ground_res = MAX(ground_res, &
                            ABS(ns_nc_ground(m1,m2,is,nah) - &
                                rho%ns_nc(m1,m2,is,nah)))
                       IF (do_partner_check) THEN
                          ground_average_value = 0.5_DP * ( &
                               ns_nc_ground(m1,m2,is,nah) + &
                               ns_nc_ground_partner(m1,m2,is,nah))
                          ground_partner_scale = MAX(ground_partner_scale, &
                               ABS(ns_nc_ground_partner(m1,m2,is,nah)))
                          ground_average_scale = MAX(ground_average_scale, &
                               ABS(ground_average_value))
                          ground_partner_res = MAX(ground_partner_res, &
                               ABS(ns_nc_ground_partner(m1,m2,is,nah) - &
                                   rho%ns_nc(m1,m2,is,nah)))
                          ground_average_res = MAX(ground_average_res, &
                               ABS(ground_average_value - &
                                   rho%ns_nc(m1,m2,is,nah)))
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
        IF (ionode) WRITE(stdout, &
             '(5x,"SOC_U_BARE_GROUND_CHECK",4(1x,es24.16))') &
             ground_scale_ph, ground_scale_pw, ground_res, &
             ground_res / MAX(ground_scale_pw, TINY(1.0_DP))
        IF (ionode .AND. do_partner_check) THEN
           WRITE(stdout, &
                '(5x,"SOC_U_BARE_GROUND_PARTNER_CHECK",4(1x,es24.16))') &
                ground_partner_scale, ground_scale_pw, ground_partner_res, &
                ground_partner_res / MAX(ground_scale_pw, TINY(1.0_DP))
           WRITE(stdout, &
                '(5x,"SOC_U_BARE_GROUND_AVERAGE_CHECK",4(1x,es24.16))') &
                ground_average_scale, ground_scale_pw, ground_average_res, &
                ground_average_res / MAX(ground_scale_pw, TINY(1.0_DP))
        ENDIF
     ENDIF
     !
     ! Filling the m1 m2 lower triangular part
     !
     IF (.NOT.noncolin) THEN
        DO m1 = 2, ldim
           DO m2 = 1, m1-1
              dnsbare(m1,m2,:,:,:,:) = dnsbare(m2,m1,:,:,:,:)
           ENDDO
        ENDDO
     ENDIF
     !
     ! In nspin.eq.1 k-point weight is normalized to 2 el/band
     ! in the whole BZ but we are interested in dns of one spin component
     !
     IF ((.NOT.noncolin) .AND. nspin.EQ.1) dnsbare = 0.5d0 * dnsbare
     !
     ! Write (unsymmetrized) dnsbare to file
     !
     IF (ionode) WRITE(iundnsbare,*) dnsbare
     !
  ENDIF
  !
  ! Close the unit iundnsbare
  !
  IF (ionode) CLOSE (unit=iundnsbare,status='keep')
  !
  ! Symmetrize dnsbare, and also save the result
  ! in the pattern basis dnsbare_all_modes
  !
  CALL sym_dns_wrapper (ldim, dnsbare, dnsbare_all_modes)
  !
  ! Write symmetryzed dnsbare in cartesian coordinates
  ! to the standard output
  !
  !added to output symmetrized dnsbare in the pattern basis
  IF(ionode) THEN
     CALL seqopn (iundnsbare, 'dnsbare_pattern', 'formatted', exst)
     WRITE(iundnsbare,*) dnsbare_all_modes
     CLOSE (unit=iundnsbare,status='keep')
  ENDIF
  ! Added from https://perturbo-code.github.io/mydoc_tutorial_hubbard.html

  IF (iverbosity==1) THEN
     WRITE(stdout,*) 'DNS_BARE SYMMETRIZED IN CARTESIAN COORDINATES'
     DO na = 1, nat
        DO icart = 1, 3
           WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') 'displaced atom L =', na, 'ipol=', icart
           DO nah = 1, nat
              nt = ityp(nah)
              IF (is_hubbard(nt)) THEN
                 DO is = 1, nspin
                    WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                    DO m1 = 1, ldim
                       WRITE(stdout,'(10(f15.10,1x))') dnsbare (m1,:,is,nah,icart,na)
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     IF (noncolin .AND. lgamma) THEN
        herm_res = 0.0_DP
        DO na = 1, nat
           DO icart = 1, 3
              DO nah = 1, nat
                 nt = ityp(nah)
                 IF (is_hubbard(nt)) THEN
                    ldim_atom = 2 * Hubbard_l(nt) + 1
                    herm_res = MAX(herm_res, ns_nc_hermiticity_residual( &
                               ldim_atom, dnsbare(:,:,:,nah,icart,na)))
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        WRITE(stdout,'(5x,"dnsbare nc Hermiticity residual =",1x,es12.5)') herm_res
     ENDIF
     WRITE(stdout,*)
  ENDIF
  !
  DEALLOCATE (dqsphi)
  DEALLOCATE (dmqsphi)
  DEALLOCATE (dwfcatom_)
  IF (noncolin) THEN
     DEALLOCATE (proj_nc_left)
     DEALLOCATE (proj_nc_right)
     DEALLOCATE (band_weight)
     IF (do_ground_check) DEALLOCATE (ns_nc_ground)
     IF (need_partner_record) DEALLOCATE (evc_partner)
     IF (do_partner_check) THEN
        DEALLOCATE (ns_nc_ground_partner)
     ENDIF
     IF (use_q_partner) THEN
        DEALLOCATE (proj_nc_partner_left)
        DEALLOCATE (proj_nc_partner_right)
        DEALLOCATE (dnsbare_nc_branch_ordinary)
        DEALLOCATE (dnsbare_nc_branch_partner)
     ENDIF
     IF (do_pw_kernel_check) THEN
        DEALLOCATE (dnsbare_pw_kernel)
        DEALLOCATE (dns_nc_pw)
        DEALLOCATE (wfcu_saved)
     ENDIF
  ELSE
     DEALLOCATE (proj1)
     DEALLOCATE (proj2)
  ENDIF
  !
  CALL stop_clock ('dnsq_bare')
  !
  RETURN
  !
END SUBROUTINE dnsq_bare
!----------------------------------------------------------------------------------------
