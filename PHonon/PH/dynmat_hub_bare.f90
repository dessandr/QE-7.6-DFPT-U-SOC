!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE dynmat_hub_bare
  !---------------------------------------------------------------------
  !! DFPT+U: This routine does two tasks:  
  !! 1) Computes \(\text{d2ns_bare}\) (very slow) or reads it;  
  !! 2) Adds to the dynamical matrix the bare part due the Hubbard term:
  !! \begin{equation}\notag
  !! \begin{split}
  !!  \sum_\text{nah} \text{Hubbard_U}&(\text{nah}) \sum_{is, m1, m2} [ 
  !!      (0.5\delta_{m1m2} - \text{ns}(m1,m2,is,\text{nah})) \cdot \text{d2ns_bare}\\
  !!      &- \text{conjg}(\text{dnsbare}( m1, m2, is, \text{nah}, \text{icart},\text{na})) 
  !!            \cdot \text{dnsbare}( m1, m2, is, \text{nah}, \text{jcart},\text{nap}) ]
  !! \end{split}
  !! \end{equation}
  !
  !!  Addition of the J0 term:
  !!  $$ \sum_\text{nah} \text{Hubbard_J0}(\text{nah}) \sum_{is, m1, m2} [ 
  !!      \text{ns}(m1,m2,-is,\text{nah}) \cdot \text{d2ns_bare}(m1, m2, is,
  !!      \text{nah}, \text{nap_jcart}, \text{na_icart})  
  !!      + \text{conjg}(\text{dnsbare}(m1,m2,is,\text{nah}, \text{icart},\text{na}))\cdot
  !!          \text{dnsbare}(m1,m2,-is, \text{nah}, \text{jcart},\text{nap}) ] $$
  !
  !! Written  by A. Floris.  
  !! Modified by I. Timrov (01.10.2018).
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ldaU,          ONLY : Hubbard_l, is_hubbard, Hubbard_J0, &
                            Hubbard_lmax, offsetU, nwfcU
  USE ldaU_ph,       ONLY : dnsbare, wfcatomk, wfcatomkpq, dvkb, vkbkpq, &
                            dvkbkpq, projpb, projpdb, read_dns_bare, d2ns_type
  USE ldaU_lr,       ONLY : effU, swfcatomk, proj1
  USE wavefunctions, ONLY : evc
  USE units_lr,      ONLY : iuwfc, lrwfc, iuatwfc, iuatswfc
  USE uspp,          ONLY : vkb, nkb, ofsbeta
  USE uspp_param,    ONLY : nh
  USE klist,         ONLY : xk, ngk, igk_k
  USE control_lr,    ONLY : lgamma
  USE io_files,      ONLY : nwordwfcU, tmp_dir 
  USE lsda_mod,      ONLY : lsda, current_spin, isk, nspin
  USE modes,         ONLY : u, nmodes
  USE dynmat,        ONLY : dyn
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE wvfct,         ONLY : npwx, nbnd
  USE control_flags, ONLY : iverbosity
  USE d2nsq_bare_module
  USE scf,           ONLY : rho
  USE mp,            ONLY : mp_sum, mp_bcast
  USE mp_pools,      ONLY : inter_pool_comm
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp_world,      ONLY : world_comm
  USE io_files,      ONLY : seqopn
  USE buffers,       ONLY : get_buffer
  ! test: variables useful to symmetrize dynmat_hub_bare, 
  ! before (and only for) printing it
  USE lr_symm_base,  ONLY : rtau, nsymq, irotmq, minus_q
  USE modes,         ONLY : u
  USE qpoint,        ONLY : xq
  USE symm_base,     ONLY : irt, s, invs
  USE cell_base,     ONLY : at, bg
  USE dynmat,        ONLY : dyn_hub_bare 
  USE io_global,     ONLY : ionode, ionode_id, stdout
  USE uspp_init,        ONLY : init_us_2
  USE noncollin_module, ONLY : noncolin
#if defined(__PH_ACC)
  USE uspp,          ONLY : okvan
#endif
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER :: icart, jcart, na, nap, nah, &
             ihubst1, ihubst2, ibeta, nt, n, l, is, ina, ih, & 
             na_icart, nap_jcart, na_icar, m1, m2, op_spin, isi, &
             imode, jmode, ldim, ik, ikk, ikq, icar, ibnd, &
             ios, iund2nsbare, npw, npwq
  REAL(DP)    :: nsaux
  COMPLEX(DP) :: dnsaux1, dnsaux2, d2ns_bare_aux, d2ns_bare_k, work
  LOGICAL :: exst
  COMPLEX(DP), ALLOCATABLE :: d2ns_bare(:,:,:,:,:,:), dynwrk(:,:)  
  !
  IF (noncolin) THEN
     CALL dynmat_hub_bare_nc()
     RETURN
  ENDIF
  !
  CALL start_clock ( 'dynmat_hub_bare' )
  !
  ios = 0
  !
  ldim = 2*Hubbard_lmax + 1
  !
  ALLOCATE (dyn_hub_bare(3*nat,3*nat))  
  ALLOCATE (d2ns_bare(ldim,ldim,nspin,nat,3*nat,3*nat) )
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (projpb(nbnd,nkb))
  ALLOCATE (projpdb(nbnd,nkb,3))
  ALLOCATE (dynwrk(3*nat,3*nat))
  !
  dyn_hub_bare = (0.d0, 0.d0)
  dynwrk       = (0.d0, 0.d0) 
  d2ns_bare    = (0.d0, 0.d0)
  !
  exst = .false.
  iund2nsbare = 39
  !
  IF (ionode) CALL seqopn(iund2nsbare, 'd2nsbare', 'formatted', exst)
  !
  ! Read dnsbare
  !
  IF (read_dns_bare) THEN 
     !
     WRITE( stdout,*) 'D2NS_BARE FILE EXISTS IN TMP_DIR AND IS READ CORRECTLY?'
     !   
     IF (ionode) THEN
        !
        IF (exst) THEN
           !
           ios=0
           READ(iund2nsbare,*,iostat=ios) d2ns_bare
           !
           IF (ios==0) THEN 
              WRITE( stdout,*) '...YES. THE D2NS_BARE MATRIX WAS CORRECTLY READ FROM FILE'
           ELSE
              WRITE( stdout,*) '...NO. THE FILE EXISTS BUT IT SEEMS TO BE CORRUPTED'
           ENDIF
           !
        ELSE
           WRITE( stdout,*) 'THE FILE DOES NOT EXIST'
        ENDIF
        !
     ENDIF
     !
     CALL mp_bcast(exst, ionode_id, world_comm)
     CALL mp_bcast(ios, ionode_id, world_comm)
     !
     IF (exst .AND. ios==0) CALL mp_bcast(d2ns_bare, ionode_id, world_comm)
     !           
  ENDIF
  !
  ! If dnsbare was not read, then compute it
  !
  IF ((.NOT.exst) .OR. (ios/=0) .OR. (.NOT.read_dns_bare)) THEN
     !     
     d2ns_bare = (0.d0, 0.d0)
     ! 
     WRITE(stdout,'(/5x,"Calculating the d2ns_bare matrix. It might take a while!")')
     !   
#if defined(__PH_ACC)
     ! ph_acc.x: with NCPP only the onsite blocks are nonzero and no
     ! beta-function derivative is needed (see d2ns_bare_ncpp_onsite).
     IF (.NOT.okvan) THEN
        CALL d2ns_bare_ncpp_onsite (ldim, d2ns_bare)
     ELSE
#endif
     DO ik = 1, nksq
        !
        ! Print information about how much the calculation has advanced
        ! because it may take really long...
        !
        WRITE(stdout,'(5x,"k point #",1x,i5,3x,"out of",1x,i5)') ik, nksq
        !
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq= ngk(ikq) 
        !
        IF (lsda) THEN 
           current_spin = isk(ikk)
           IF ((current_spin==1).and.(nspin==2)) THEN
              op_spin = 2
           ELSE
              op_spin = 1
           ENDIF
        ELSE        
           op_spin = 1     
        ENDIF
        !
        ! Read unperturbed wavefuctions psi(k)
        !
        IF (nksq.GT.1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
        !
        ! Compute the beta function vkb at k, and vkbkpq at k+q
        !
        CALL init_us_2 (npw, igk_k(1,ikk), xk(1,ikk), vkb)
        IF (.NOT.lgamma) CALL init_us_2 (npwq, igk_k(1,ikq), xk(1,ikq), vkbkpq)
        !
        ! Read the atomic orbitals \phi at k and k+q from file (unit iuatwfc)
        ! 
        CALL get_buffer (wfcatomk, nwordwfcU, iuatwfc, ikk)
        IF (.NOT.lgamma) CALL get_buffer (wfcatomkpq, nwordwfcU, iuatwfc, ikq)
        !
        ! Read S*\phi at k from file (unit iuatswfc)
        !
        CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
        !
        ! Calculate proj1(ibnd) = <S(k)*phi(k,I,m')|psi(ibnd,k)> 
        !
        proj1 = (0.d0, 0.d0)
        !
        !FIXME use zgemms here 
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO m1 = 1, 2*Hubbard_l(nt)+1
                 ihubst1 = offsetU(nah) + m1   ! I m index
                 DO ibnd = 1, nbnd
                 proj1(ibnd, ihubst1) = dot_product (swfcatomk(1:npw,ihubst1), evc(1:npw,ibnd))
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
        !
        ! Calculate: projpb(nbnd,nkb) = <psi|vkb> for all nbnd and vkb states, 
        ! useful in d2ns_bare_k 
        !
        projpb = (0.d0, 0.d0)
        !
        !FIXME zgemms of even calbec could be use here 
        DO na = 1, nat    
           nt = ityp(na)  
           DO ih = 1, nh(nt)
              ibeta = ofsbeta(na) + ih
              DO ibnd = 1, nbnd
                 projpb(ibnd, ibeta) = dot_product (evc(1:npw,ibnd), vkb(1:npw,ibeta))
              ENDDO
           ENDDO
        ENDDO
        !
        CALL mp_sum(proj1, intra_bgrp_comm)
        CALL mp_sum(projpb, intra_bgrp_comm)
        !
        ! Calculate the m1 m2 upper triangular part (UPT) of the UPT matrices
        !
        projpdb = (0.d0, 0.d0)
        !
        DO na = 1, nat  ! the displaced atom
           !
           DO icart = 1, 3
              !
              na_icart = 3*(na-1) + icart
              !
              DO nap = 1, nat ! the displaced atom
                 !
                 DO jcart = 1, 3
                    !
                    nap_jcart = 3*(nap-1) + jcart
                    !
                    IF (nap_jcart.GE.na_icart) THEN
                       !
                       ! Calculate the derivatives of the beta functions  
                       ! at k and k+q only at icart, jcart
                       ! and only afor atoms na, nap (for all the states)
                       !                  
                       DO icar = 1, 3
                          !
                          IF ((icar==icart) .OR. (icar==jcart)) THEN 
                             ! we want only icart jcart
                             !
                             DO nt = 1, ntyp
                                !
                                DO ina = 1, nat 
                                   !
                                   IF (ityp(ina).EQ.nt) THEN
                                      !
                                      IF ((ina==na) .OR. (ina==nap)) THEN 
                                         ! we want only dbeta at na or nap
                                         !
                                         DO ih = 1, nh(nt)
                                            !
                                            ibeta = ofsbeta(ina) + ih
                                            !
                                            CALL dwfc(npw, igk_k(:,ikk), ikk, icar, &
                                                      vkb(:,ibeta), dvkb(:,ibeta,icar))
                                            IF (.NOT.lgamma) &
                                            CALL dwfc(npwq, igk_k(:,ikq), ikq, icar, &
                                                      vkbkpq(:,ibeta), dvkbkpq(:,ibeta,icar))
                                            ! 
                                            DO ibnd = 1, nbnd
                                               projpdb(ibnd, ibeta, icar) = &
                                                       dot_product(evc(1:npw,ibnd), dvkb(1:npw,ibeta,icar))
                                            ENDDO
                                            !
                                         ENDDO ! ih
                                         !
                                      ENDIF
                                      !
                                   ENDIF
                                   !
                                ENDDO ! ina
                                !
                             ENDDO ! nt
                             !
                          ENDIF
                          !
                       ENDDO ! icar
                       !
                       CALL mp_sum(projpdb, intra_bgrp_comm)
                       !
                       DO nah = 1, nat ! the Hubbard atom
                          !
                          nt = ityp(nah)
                          !
                          IF (is_hubbard(nt)) THEN
                             !                          
                             DO m1 = 1, 2*Hubbard_l(nt)+1
                                !
                                ihubst1 = offsetU(nah) + m1  ! I m index
                                !
                                DO m2 = m1, 2*Hubbard_l(nt)+1
                                   !
                                   ihubst2 = offsetU(nah) + m2 ! I m index
                                   !
                                   ! Calculate d2ns_bare_k  = the bare part of 
                                   ! the 2nd derivative of the occupation matrix ns 
                                   ! it is an object calculated at 
                                   ! (ik,na,icart,nap,jcart,ihubst1,ihubst2)
                                   !
                                   IF ((d2ns_type=='full') .OR. (d2ns_type=='fmmp')) THEN
                                      !           
                                      CALL d2nsq_bare_k (ik, icart, jcart, na, nap, nah, &
                                                         ihubst1, ihubst2, d2ns_bare_k)
                                      ! 
                                   ELSEIF ((d2ns_type=='diag') .OR. (d2ns_type=='dmmp')) THEN
                                      !
                                      CALL d2nsq_bare_k_diag (ik, icart, jcart, na, nap, nah, &
                                                              ihubst1, ihubst2, d2ns_bare_k)
                                      !
                                   ELSE
                                      CALL errore( 'dynmat_hub_bare', &
                                                    & 'wrong d2ns_type value', 1 )
                                   ENDIF
                                   !
                                   d2ns_bare (m1, m2, current_spin, nah, nap_jcart, na_icart) = &
                                   d2ns_bare (m1, m2, current_spin, nah, nap_jcart, na_icart) + d2ns_bare_k    
                                   ! 
                                ENDDO ! m2
                                !
                             ENDDO ! m1
                             !
                          ENDIF
                          !
                       ENDDO ! nah
                       !
                    ENDIF  ! IF  (nap_jcart >= na_icart) 
                    ! 
                 ENDDO ! jcart
                 !
              ENDDO ! nap
              !
           ENDDO ! icart
           !
        ENDDO ! na
        !
     ENDDO ! ik
     !
     CALL mp_sum(d2ns_bare, inter_pool_comm)
#if defined(__PH_ACC)
     ENDIF
#endif
     !
     ! If nspin=1, k point weight is normalized to 2 el/band 
     !
     IF (nspin.EQ.1) d2ns_bare = 0.5d0 * d2ns_bare  
     !   
     ! Filling the m1 m2 lower triangular part (LTP) of the calculated UPT matrices
     !   
     DO na_icart = 1, 3*nat
        DO nap_jcart = na_icart, 3*nat
           DO m1 = 2, ldim
              DO m2 = 1, m1-1
                 DO nah = 1, nat ! the Hubbard atom
                    nt = ityp(nah)
                    IF (is_hubbard(nt)) THEN
                       DO is = 1, nspin
                          d2ns_bare (m1, m2, is, nah, nap_jcart, na_icart) = &
                          d2ns_bare (m2, m1, is, nah, nap_jcart, na_icart)   
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !   
     ! Filling the LTP matrices
     !   
     DO na_icart = 2, 3*nat
        DO nap_jcart = 1, na_icart-1
           d2ns_bare (:,:,:,:, nap_jcart, na_icart) = &
           d2ns_bare (:,:,:,:, na_icart, nap_jcart)   
        ENDDO
     ENDDO
     !
     ! Write d2ns_bare to file
     !
     IF (ionode) WRITE(iund2nsbare,*) d2ns_bare
     !
  ENDIF
  !
  IF (ionode) CLOSE(unit=iund2nsbare,status='KEEP')
  !   
  ! Writing the full d2ns_bare matrix to the standard output file
  !
  IF (iverbosity==1) THEN 
     WRITE(stdout,*)
     WRITE(stdout,*) 'd2ns_bare NOT SYMMETRIZED IN CARTESIAN COORDINATES' 
     DO na = 1, nat  ! the displaced atom
        DO icart = 1, 3
           na_icart = 3*(na-1) + icart
           DO nap = 1, nat ! the displaced atom
              DO jcart = 1, 3
                 nap_jcart = 3*(nap-1) + jcart
                 DO nah = 1, nat ! the Hubbard atom
                    nt = ityp(nah)
                    IF (is_hubbard(nt)) THEN
                       DO is = 1, nspin
                          WRITE(stdout,'(2(a,1x,i2,2x,a,1x,i2,2x))') &
                                  'displaced atom L =', na, 'ipol=', &
                                   icart,'displaced atom Lp =', nap, 'ipol=', jcart
                          WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                          DO m1 = 1, ldim
                             WRITE(stdout,'(14(f15.10,1x))') &
                                  (d2ns_bare(m1, m2, is, nah, nap_jcart, na_icart), m2 = 1, ldim)
                          ENDDO
                       ENDDO
                    ENDIF 
                 ENDDO 
              ENDDO 
           ENDDO  
        ENDDO 
     ENDDO
     WRITE( stdout,*)
  ENDIF
  !
  ! And now we compute the bare part of the dynamical matrix and add it to dyn
  !      
  DO na = 1, nat  ! the displaced atom
     !
     DO icart = 1, 3
        !
        na_icart = 3*(na-1) + icart
        !      
        DO nap = 1, nat ! the displaced atom
           !
           DO jcart = 1, 3
              !
              nap_jcart = 3*(nap-1) + jcart
              !      
              DO nah = 1, nat ! the Hubbard atom
                 !
                 nt = ityp(nah)
                 !
                 ! The effective Hubbard_U part
                 ! 
                 IF (is_hubbard(nt)) THEN
                    !
                    ! work = \sum_{is, m1, m2} [ 
                    ! (0.5*delta_m1m2- ns(m1,m2,is,nah))*d2ns_bare - 
                    ! - conjg(dnsbare ( m1, m2, is, nah, icart,na))*
                    !         dnsbare ( m1, m2, is, nah, jcart,nap) ]
                    !      
                    work = (0.d0, 0.d0)
                    !       
                    DO is = 1, nspin
                       !
                       DO m1 = 1, 2*Hubbard_l(nt) + 1
                          !
                          DO m2 = 1, 2*Hubbard_l(nt) + 1
                             !      
                             nsaux = rho%ns(m1, m2, is, nah)
                             !      
                             ! dnsbare matrix is symmetric i.e. dnsbare(m1,m2) = dnsbare(m2,m1) 
                             ! (when keeping only the terms in u and neglecting the ones in u*)
                             !      
                             dnsaux1 = dnsbare(m1, m2, is, nah, icart, na)
                             dnsaux2 = dnsbare(m1, m2, is, nah, jcart, nap)
                             d2ns_bare_aux = d2ns_bare(m1, m2, is, nah, nap_jcart, na_icart) 
                             !                              
                             ! Include the delta_m1m2 contribution              
                             !
                             IF (m1==m2) work = work + 0.5d0*d2ns_bare_aux
                             !      
                             work = work - nsaux   * d2ns_bare_aux &
                                         - dnsaux1 * CONJG(dnsaux2)
                             ! 
                          ENDDO ! m2
                          !
                       ENDDO ! m1
                       !      
                    ENDDO ! is
                    ! 
                    ! Summing the contributions for each Hubbard atom                    
                    !
                    dynwrk(na_icart, nap_jcart) = dynwrk(na_icart, nap_jcart) & 
                                                  + work * effU(nt)
                    !      
                 ENDIF
                 !
                 ! The Hubbard_J0 part (if present)
                 !
                 IF (Hubbard_J0(nt).NE.0.d0) THEN
                    !                     
                    ! work = \sum_{is, m1, m2} [ 
                    !    ns(m1,m2,-is,nah) * d2ns_bare - 
                    !    dnsbare( m1, m2, is, nah, icart,na) *
                    !    conjg(dnsbare( m1, m2, -is, nah, jcart,nap)) ]
                    !      
                    work = (0.d0, 0.d0)
                    !      
                    DO is = 1, nspin
                       !      
                       isi = nspin/is
                       !
                       DO m1 = 1, 2*Hubbard_l(nt) + 1
                          !
                          DO m2 = 1, 2*Hubbard_l(nt) + 1
                             !      
                             nsaux = rho%ns(m1, m2, isi, nah)
                             !      
                             ! dnsbare matrix is symmetric i.e. dnsbare(m1,m2) = dnsbare(m2,m1) 
                             ! (when keeping only the terms in u and neglecting the ones in u*)
                             !      
                             dnsaux1 = dnsbare(m1, m2, is,  nah, icart, na)
                             dnsaux2 = dnsbare(m1, m2, isi, nah, jcart, nap) 
                             d2ns_bare_aux = d2ns_bare(m1, m2, is, nah, nap_jcart, na_icart) 
                             !                              
                             ! DO NOT include the delta_m1m2 contribution  
                             ! Note the sign change            
                             !
                             work = work + nsaux *  d2ns_bare_aux + dnsaux1 * CONJG(dnsaux2) 
                             ! 
                          ENDDO ! m2
                          !
                       ENDDO ! m1
                       !      
                    ENDDO ! is
                    ! 
                    ! Summing the contributions for each Hubbard atom                    
                    !
                    dynwrk(na_icart, nap_jcart) = dynwrk(na_icart, nap_jcart) & 
                                                  + work * Hubbard_J0(nt)
                    !
                 ENDIF
                 !
              ENDDO ! nah
              !
           ENDDO ! jcart
           !
        ENDDO  ! nap
        !
     ENDDO ! icart
     !
  ENDDO ! na
  !      
  IF (nspin.EQ.1) dynwrk = 2.d0 * dynwrk  
  !      
  ! Transform from the Cartesian to the pattern basis u. 
  ! We rotate the dynamical matrix on the basis of patterns u.
  !
  DO imode = 1, nmodes
     DO jmode = 1, nmodes
        !
        work = (0.d0, 0.d0)
        !
        DO nap_jcart = 1, 3*nat
           DO na_icart = 1, 3*nat
              work = work + CONJG (u(na_icart, imode) ) * &
                            dynwrk(na_icart, nap_jcart) * &
                            u(nap_jcart, jmode)
           ENDDO
        ENDDO
        !
        ! Adding the contribution to the total dynamical matrix dyn
        !        
        dyn(imode, jmode) = dyn(imode, jmode) + work
        !     
        ! Keep separately the bare part of the dynamical matrix
        ! 
        dyn_hub_bare(imode, jmode) = dyn_hub_bare(imode, jmode) + work
        !
     ENDDO
  ENDDO
  !
  IF (iverbosity==1) THEN
     !
     ! Write the UNSYMMETRIZED bare Hubbard dynamical matrix 
     ! in the pattern basis
     !
     CALL tra_write_matrix_no_sym ('dyn_hub_bare NOT SYMMETRIZED',dyn_hub_bare,nat) 
     !
     ! Write the SYMMETRIZED bare Hubbard dynamical matrix 
     ! in carthesian coordinates
     !
     CALL tra_write_matrix ('dyn_hub_bare SYMMETRIZED',dyn_hub_bare,u,nat)                   
     !
  ENDIF
  !
  DEALLOCATE (proj1)
  DEALLOCATE (projpb)
  DEALLOCATE (projpdb)
  DEALLOCATE (d2ns_bare)
  DEALLOCATE (dynwrk)     
  !
  CALL stop_clock ( 'dynmat_hub_bare' )
  ! 
  RETURN
  !      
END SUBROUTINE dynmat_hub_bare

#if defined(__PH_ACC)
!-----------------------------------------------------------------------
SUBROUTINE d2ns_bare_ncpp_onsite (ldim, d2ns_bare)
  !---------------------------------------------------------------------
  !! ph_acc.x only: collinear NCPP d2ns_bare without the USPP machinery.
  !! With S=1 only term_one_1 and term_three_1 of d2nsq_bare_k(_diag)
  !! survive, and only for na = nap = nah:
  !!   term1(m1,m2) = \sum_n wg <psi_n|d2_{ij} phi_m1> <S phi_m2|psi_n>
  !!   term3(m1,m2) = \sum_n wg <psi_n|d_i phi_m1> <d_j phi_m2|psi_n>
  !! term2/term4 are term1/term3 with m1 <=> m2 ('full', 'diag') or
  !! repeated ('fmmp', 'dmmp'), and are added in the original order.
  !! The projections are computed once per k point instead of once per
  !! (na,icart,nap,jcart,nah,m1,m2), and only the onsite blocks are
  !! pool-summed. The same m1<=m2, jcart>=icart entries as the original
  !! loop are filled; dynmat_hub_bare completes the triangles.
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp
  USE ldaU,          ONLY : Hubbard_l, is_hubbard, offsetU, nwfcU
  USE ldaU_ph,       ONLY : wfcatomk, d2ns_type
  USE ldaU_lr,       ONLY : swfcatomk
  USE wavefunctions, ONLY : evc
  USE units_lr,      ONLY : iuwfc, lrwfc, iuatwfc, iuatswfc
  USE klist,         ONLY : ngk, igk_k
  USE io_files,      ONLY : nwordwfcU
  USE lsda_mod,      ONLY : lsda, current_spin, isk, nspin
  USE qpoint,        ONLY : nksq, ikks
  USE wvfct,         ONLY : npwx, nbnd, wg
  USE buffers,       ONLY : get_buffer
  USE mp,            ONLY : mp_sum
  USE mp_pools,      ONLY : inter_pool_comm
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE io_global,     ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldim
  COMPLEX(DP), INTENT(INOUT) :: d2ns_bare(ldim,ldim,nspin,nat,3*nat,3*nat)
  !
  INTEGER :: ik, ikk, npw, nah, nt, m, m1, m2, ihubst, ihubst1, ihubst2, &
             icart, jcart, ipair, ibnd
  INTEGER :: pair(3,3)
  LOGICAL :: mmp_sym
  COMPLEX(DP) :: res_one, res_two, res_three, res_four, d2ns_bare_k
  COMPLEX(DP), ALLOCATABLE :: dphi(:,:,:), d2phi(:,:,:), proj0(:,:), &
                              projd1(:,:,:), projd2(:,:,:), onsite(:,:,:,:,:,:)
  !
  IF ((d2ns_type/='full') .AND. (d2ns_type/='fmmp') .AND. &
      (d2ns_type/='diag') .AND. (d2ns_type/='dmmp')) &
     CALL errore('dynmat_hub_bare', 'wrong d2ns_type value', 1)
  ! With NCPP 'diag'/'dmmp' coincide with 'full'/'fmmp'.
  mmp_sym = (d2ns_type=='fmmp') .OR. (d2ns_type=='dmmp')
  !
  ! Index of the symmetric second derivative d2_{icart,jcart}, icart<=jcart
  !
  ipair = 0
  DO icart = 1, 3
     DO jcart = icart, 3
        ipair = ipair + 1
        pair(icart,jcart) = ipair
     ENDDO
  ENDDO
  !
  ALLOCATE (dphi(npwx,nwfcU,3))
  ALLOCATE (d2phi(npwx,nwfcU,6))
  ALLOCATE (proj0(nwfcU,nbnd))
  ALLOCATE (projd1(nbnd,nwfcU,3))
  ALLOCATE (projd2(nbnd,nwfcU,6))
  ALLOCATE (onsite(ldim,ldim,nspin,3,3,nat))
  !
  ! Columns outside the Hubbard manifolds are never written and stay zero
  !
  dphi   = (0.d0, 0.d0)
  d2phi  = (0.d0, 0.d0)
  onsite = (0.d0, 0.d0)
  !
  DO ik = 1, nksq
     !
     WRITE(stdout,'(5x,"k point #",1x,i5,3x,"out of",1x,i5)') ik, nksq
     !
     ikk = ikks(ik)
     npw = ngk(ikk)
     IF (lsda) current_spin = isk(ikk)
     !
     IF (nksq.GT.1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
     CALL get_buffer (wfcatomk, nwordwfcU, iuatwfc, ikk)
     CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
     !
     ! First and second derivatives of the Hubbard orbitals at k
     !
     DO nah = 1, nat
        nt = ityp(nah)
        IF (is_hubbard(nt)) THEN
           DO m = 1, 2*Hubbard_l(nt)+1
              ihubst = offsetU(nah) + m
              DO icart = 1, 3
                 CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                            wfcatomk(:,ihubst), dphi(:,ihubst,icart))
                 DO jcart = icart, 3
                    CALL d2wfc (npw, igk_k(1,ikk), ikk, icart, jcart, &
                                wfcatomk(:,ihubst), d2phi(:,ihubst,pair(icart,jcart)))
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     !
     ! proj0(I,ibnd) = <S phi_I|psi_ibnd>, the original proj1(ibnd,I)
     ! projd1(ibnd,I,icart) = <psi_ibnd|d_icart phi_I>, S=1
     ! projd2(ibnd,I,ipair) = <psi_ibnd|d2_ipair phi_I>, S=1
     !
     CALL ZGEMM ('C', 'N', nwfcU, nbnd, npw, (1.d0,0.d0), swfcatomk, npwx, &
                 evc, npwx, (0.d0,0.d0), proj0, nwfcU)
     DO icart = 1, 3
        CALL ZGEMM ('C', 'N', nbnd, nwfcU, npw, (1.d0,0.d0), evc, npwx, &
                    dphi(:,:,icart), npwx, (0.d0,0.d0), projd1(:,:,icart), nbnd)
     ENDDO
     DO ipair = 1, 6
        CALL ZGEMM ('C', 'N', nbnd, nwfcU, npw, (1.d0,0.d0), evc, npwx, &
                    d2phi(:,:,ipair), npwx, (0.d0,0.d0), projd2(:,:,ipair), nbnd)
     ENDDO
     CALL mp_sum(proj0,  intra_bgrp_comm)
     CALL mp_sum(projd1, intra_bgrp_comm)
     CALL mp_sum(projd2, intra_bgrp_comm)
     !
     DO nah = 1, nat
        nt = ityp(nah)
        IF (.NOT.is_hubbard(nt)) CYCLE
        DO icart = 1, 3
           DO jcart = icart, 3
              ipair = pair(icart,jcart)
              DO m1 = 1, 2*Hubbard_l(nt)+1
                 ihubst1 = offsetU(nah) + m1
                 DO m2 = m1, 2*Hubbard_l(nt)+1
                    ihubst2 = offsetU(nah) + m2
                    !
                    res_one = (0.d0, 0.d0)
                    res_three = (0.d0, 0.d0)
                    DO ibnd = 1, nbnd
                       res_one = res_one + wg(ibnd,ikk) * projd2(ibnd,ihubst1,ipair) * &
                                           proj0(ihubst2,ibnd)
                       res_three = res_three + wg(ibnd,ikk) * projd1(ibnd,ihubst1,icart) * &
                                               CONJG(projd1(ibnd,ihubst2,jcart))
                    ENDDO
                    !
                    IF (mmp_sym .OR. (ihubst1==ihubst2)) THEN
                       res_two = res_one
                       res_four = res_three
                    ELSE
                       res_two = (0.d0, 0.d0)
                       res_four = (0.d0, 0.d0)
                       DO ibnd = 1, nbnd
                          res_two = res_two + wg(ibnd,ikk) * projd2(ibnd,ihubst2,ipair) * &
                                              proj0(ihubst1,ibnd)
                          res_four = res_four + wg(ibnd,ikk) * projd1(ibnd,ihubst2,icart) * &
                                                CONJG(projd1(ibnd,ihubst1,jcart))
                       ENDDO
                    ENDIF
                    !
                    d2ns_bare_k = res_one
                    d2ns_bare_k = d2ns_bare_k + res_two
                    d2ns_bare_k = d2ns_bare_k + res_three
                    d2ns_bare_k = d2ns_bare_k + res_four
                    !
                    onsite(m1,m2,current_spin,jcart,icart,nah) = &
                    onsite(m1,m2,current_spin,jcart,icart,nah) + d2ns_bare_k
                    !
                 ENDDO ! m2
              ENDDO ! m1
           ENDDO ! jcart
        ENDDO ! icart
     ENDDO ! nah
     !
  ENDDO ! ik
  !
  CALL mp_sum(onsite, inter_pool_comm)
  !
  ! Same storage as the original loop: (nap_jcart, na_icart) with na=nap=nah
  !
  DO nah = 1, nat
     IF (.NOT.is_hubbard(ityp(nah))) CYCLE
     DO icart = 1, 3
        DO jcart = icart, 3
           d2ns_bare(:,:,:,nah,3*(nah-1)+jcart,3*(nah-1)+icart) = &
                onsite(:,:,:,jcart,icart,nah)
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE (dphi)
  DEALLOCATE (d2phi)
  DEALLOCATE (proj0)
  DEALLOCATE (projd1)
  DEALLOCATE (projd2)
  DEALLOCATE (onsite)
  !
  RETURN
  !
END SUBROUTINE d2ns_bare_ncpp_onsite
#endif

!-----------------------------------------------------------------------
SUBROUTINE dynmat_hub_bare_nc
  !---------------------------------------------------------------------
  !! NCPP fixed-Dudarev bare Hubbard dynamical matrix in the complete
  !! orbital-spin space.  At a general q the two first-order occupation
  !! responses are contracted through their q/-q adjoint relation; they
  !! are not made Hermitian at the same q.
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ityp, ntyp => nsp
  USE ldaU,              ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, &
                                Hubbard_J0, Hubbard_alpha, Hubbard_beta, &
                                Hubbard_projectors, lda_plus_u_kind, nwfcU
  USE ldaU_ph,           ONLY : dnsbare, wfcatomk
  USE ldaU_lr,           ONLY : effU, swfcatomk, pack_hubbard_blocks_nc, &
                                build_dvhub_nc, accumulate_ns_nc_transition, &
                                ns_nc_hermiticity_residual
  USE wavefunctions,     ONLY : evc
  USE units_lr,          ONLY : iuwfc, lrwfc, iuatwfc, iuatswfc
  USE klist,             ONLY : ngk, igk_k
  USE qpoint,            ONLY : nksq, ikks
  USE wvfct,             ONLY : npwx, nbnd, wg
  USE control_lr,        ONLY : lgamma
  USE io_files,          ONLY : nwordwfcU
  USE buffers,           ONLY : get_buffer
  USE noncollin_module,  ONLY : npol
  USE uspp,              ONLY : okvan
  USE scf,               ONLY : v
  USE modes,             ONLY : u, nmodes
  USE dynmat,            ONLY : dyn, dyn_hub_bare
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp_pools,          ONLY : inter_pool_comm
  USE mp,                ONLY : mp_sum
  USE control_flags,     ONLY : iverbosity
  USE io_global,         ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ikk, npw, na, nap, nt, nah, is_col, is_row, m, iwfc
  INTEGER :: icart, jcart, na_icart, nap_jcart
  INTEGER :: ldim_atom, max_ldim, matrix_size, max_matrix_size
  INTEGER :: row_first, imode, jmode
  REAL(DP), ALLOCATABLE :: band_weight(:)
  REAL(DP) :: d2n_herm_res, dyn_herm_res
  COMPLEX(DP) :: work
  COMPLEX(DP), EXTERNAL :: ZDOTU, ZDOTC
  COMPLEX(DP), ALLOCATABLE :: projector_d1(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: projector_d2(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: proj0(:,:), proj_d1(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: proj_d2(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: d2ns_bare(:,:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dynwrk(:,:)
  COMPLEX(DP), ALLOCATABLE :: ground_matrix(:,:), d2n_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: dn_j_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: dv_i_matrix(:,:)
  !
  CALL start_clock('dynmat_hub_bare')
  !
  IF (okvan) CALL errore('dynmat_hub_bare_nc', &
                         'noncollinear Hubbard response requires NCPP', 1)
  IF (Hubbard_projectors /= 'atomic') CALL errore('dynmat_hub_bare_nc', &
       'noncollinear Hubbard response requires atomic projectors', 1)
  IF (lda_plus_u_kind /= 0) CALL errore('dynmat_hub_bare_nc', &
                                        'only Dudarev DFT+U is supported', 1)
  IF (npol /= 2) CALL errore('dynmat_hub_bare_nc', &
                             'noncollinear Hubbard response requires npol=2', 1)
  IF (.NOT.ALLOCATED(v%ns_nc)) CALL errore('dynmat_hub_bare_nc', &
                                           'PW Hubbard potential is unavailable', 1)
  DO nt = 1, ntyp
     IF (is_hubbard(nt) .AND. Hubbard_J0(nt) /= 0.0_DP) &
        CALL errore('dynmat_hub_bare_nc', 'Hubbard J0 is outside the MVP', 1)
     IF (is_hubbard(nt) .AND. Hubbard_alpha(nt) /= 0.0_DP) &
        CALL errore('dynmat_hub_bare_nc', 'Hubbard alpha is outside the MVP', 1)
     IF (is_hubbard(nt) .AND. Hubbard_beta(nt) /= 0.0_DP) &
        CALL errore('dynmat_hub_bare_nc', 'Hubbard beta is outside the MVP', 1)
  ENDDO
  !
  max_ldim = 2 * Hubbard_lmax + 1
  max_matrix_size = max_ldim * npol
  ALLOCATE(projector_d1(npwx*npol,nwfcU,3))
  ALLOCATE(projector_d2(npwx*npol,nwfcU,3,3))
  ALLOCATE(proj0(nwfcU,nbnd))
  ALLOCATE(proj_d1(nwfcU,nbnd,3))
  ALLOCATE(proj_d2(nwfcU,nbnd,3,3))
  ALLOCATE(band_weight(nbnd))
  ALLOCATE(d2ns_bare(max_ldim,max_ldim,npol*npol,nat,3*nat,3*nat))
  ALLOCATE(dynwrk(3*nat,3*nat))
  ALLOCATE(ground_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(d2n_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(dn_j_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(dv_i_matrix(max_matrix_size,max_matrix_size))
  !
  d2ns_bare = (0.0_DP, 0.0_DP)
  dynwrk = (0.0_DP, 0.0_DP)
  !
  DO ik = 1, nksq
     ikk = ikks(ik)
     npw = ngk(ikk)
     IF (nksq > 1) CALL get_buffer(evc, lrwfc, iuwfc, ikk)
     CALL get_buffer(wfcatomk, nwordwfcU, iuatwfc, ikk)
     CALL get_buffer(swfcatomk, nwordwfcU, iuatswfc, ikk)
     !
     projector_d1 = (0.0_DP, 0.0_DP)
     projector_d2 = (0.0_DP, 0.0_DP)
     DO nah = 1, nat
        nt = ityp(nah)
        IF (is_hubbard(nt)) THEN
           ldim_atom = 2 * Hubbard_l(nt) + 1
           DO is_col = 1, npol
              DO m = 1, ldim_atom
                 iwfc = offsetU(nah) + m + ldim_atom * (is_col - 1)
                 DO icart = 1, 3
                    DO is_row = 1, npol
                       row_first = 1 + npwx * (is_row - 1)
                       CALL dwfc(npw, igk_k(1,ikk), ikk, icart, &
                                 wfcatomk(row_first,iwfc), &
                                 projector_d1(row_first,iwfc,icart))
                    ENDDO
                    DO jcart = icart, 3
                       DO is_row = 1, npol
                          row_first = 1 + npwx * (is_row - 1)
                          CALL d2wfc(npw, igk_k(1,ikk), ikk, icart, jcart, &
                                    wfcatomk(row_first,iwfc), &
                                    projector_d2(row_first,iwfc,icart,jcart))
                          IF (jcart /= icart) &
                             projector_d2(row_first:row_first+npwx-1,iwfc, &
                                          jcart,icart) = &
                             projector_d2(row_first:row_first+npwx-1,iwfc, &
                                          icart,jcart)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     !
     CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                (1.0_DP,0.0_DP), swfcatomk, npwx*npol, evc, npwx*npol, &
                (0.0_DP,0.0_DP), proj0, nwfcU)
     DO icart = 1, 3
        CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                   (1.0_DP,0.0_DP), projector_d1(:,:,icart), npwx*npol, &
                   evc, npwx*npol, (0.0_DP,0.0_DP), &
                   proj_d1(:,:,icart), nwfcU)
        DO jcart = 1, 3
           CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                      (1.0_DP,0.0_DP), projector_d2(:,:,icart,jcart), &
                      npwx*npol, evc, npwx*npol, (0.0_DP,0.0_DP), &
                      proj_d2(:,:,icart,jcart), nwfcU)
        ENDDO
     ENDDO
     CALL mp_sum(proj0, intra_bgrp_comm)
     CALL mp_sum(proj_d1, intra_bgrp_comm)
     CALL mp_sum(proj_d2, intra_bgrp_comm)
     band_weight(1:nbnd) = wg(1:nbnd,ikk)
     !
     DO nah = 1, nat
        nt = ityp(nah)
        IF (is_hubbard(nt)) THEN
           ldim_atom = 2 * Hubbard_l(nt) + 1
           DO icart = 1, 3
              na_icart = 3 * (nah - 1) + icart
              DO jcart = 1, 3
                 nap_jcart = 3 * (nah - 1) + jcart
                 CALL accumulate_ns_nc_transition(ldim_atom, offsetU(nah), &
                      npol, nbnd, band_weight, proj_d2(:,:,icart,jcart), &
                      proj0, d2ns_bare(:,:,:,nah,nap_jcart,na_icart), .FALSE.)
                 CALL accumulate_ns_nc_transition(ldim_atom, offsetU(nah), &
                      npol, nbnd, band_weight, proj_d1(:,:,icart), &
                      proj_d1(:,:,jcart), &
                      d2ns_bare(:,:,:,nah,nap_jcart,na_icart), .FALSE.)
                 CALL accumulate_ns_nc_transition(ldim_atom, offsetU(nah), &
                      npol, nbnd, band_weight, proj_d1(:,:,jcart), &
                      proj_d1(:,:,icart), &
                      d2ns_bare(:,:,:,nah,nap_jcart,na_icart), .FALSE.)
                 CALL accumulate_ns_nc_transition(ldim_atom, offsetU(nah), &
                      npol, nbnd, band_weight, proj0, &
                      proj_d2(:,:,icart,jcart), &
                      d2ns_bare(:,:,:,nah,nap_jcart,na_icart), .FALSE.)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  CALL mp_sum(d2ns_bare, inter_pool_comm)
  !
  d2n_herm_res = 0.0_DP
  DO nah = 1, nat
     nt = ityp(nah)
     IF (is_hubbard(nt)) THEN
        ldim_atom = 2 * Hubbard_l(nt) + 1
        DO icart = 1, 3
           na_icart = 3 * (nah - 1) + icart
           DO jcart = 1, 3
              nap_jcart = 3 * (nah - 1) + jcart
              d2n_herm_res = MAX(d2n_herm_res, &
                   ns_nc_hermiticity_residual(ldim_atom, &
                   d2ns_bare(:,:,:,nah,nap_jcart,na_icart)))
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  DO na = 1, nat
     DO icart = 1, 3
        na_icart = 3 * (na - 1) + icart
        DO nap = 1, nat
           DO jcart = 1, 3
              nap_jcart = 3 * (nap - 1) + jcart
              work = (0.0_DP, 0.0_DP)
              DO nah = 1, nat
                 nt = ityp(nah)
                 IF (is_hubbard(nt) .AND. nah == na .AND. nah == nap) THEN
                    ldim_atom = 2 * Hubbard_l(nt) + 1
                    matrix_size = ldim_atom * npol
                    CALL pack_hubbard_blocks_nc(ldim_atom, npol, &
                         v%ns_nc(:,:,:,nah), ground_matrix)
                    CALL pack_hubbard_blocks_nc(ldim_atom, npol, &
                         d2ns_bare(:,:,:,nah,nap_jcart,na_icart), d2n_matrix)
                    CALL pack_hubbard_blocks_nc(ldim_atom, npol, &
                         dnsbare(:,:,:,nah,jcart,nap), dn_j_matrix)
                    CALL build_dvhub_nc(ldim_atom, npol, effU(nt), &
                         dnsbare(:,:,:,nah,icart,na), dv_i_matrix)
                    work = work + ZDOTU(matrix_size*matrix_size, &
                                        ground_matrix, 1, d2n_matrix, 1)
                    IF (lgamma) THEN
                       work = work + ZDOTU(matrix_size*matrix_size, &
                                           dv_i_matrix, 1, dn_j_matrix, 1)
                    ELSE
                       ! dN_j(-q)=dN_j(q)^dagger.  ZDOTC supplies the
                       ! conjugation order used by the collinear
                       ! general-q contraction without modifying dN(q).
                       work = work + ZDOTC(matrix_size*matrix_size, &
                                           dn_j_matrix, 1, dv_i_matrix, 1)
                    ENDIF
                 ENDIF
              ENDDO
              dynwrk(na_icart,nap_jcart) = work
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  dyn_herm_res = 0.0_DP
  DO nap_jcart = 1, 3*nat
     DO na_icart = 1, 3*nat
        dyn_herm_res = MAX(dyn_herm_res, &
             ABS(dynwrk(na_icart,nap_jcart) - &
                 CONJG(dynwrk(nap_jcart,na_icart))))
     ENDDO
  ENDDO
  !
  ALLOCATE(dyn_hub_bare(3*nat,3*nat))
  dyn_hub_bare = (0.0_DP, 0.0_DP)
  DO imode = 1, nmodes
     DO jmode = 1, nmodes
        work = (0.0_DP, 0.0_DP)
        DO nap_jcart = 1, 3*nat
           DO na_icart = 1, 3*nat
              work = work + CONJG(u(na_icart,imode)) * &
                              dynwrk(na_icart,nap_jcart) * &
                              u(nap_jcart,jmode)
           ENDDO
        ENDDO
        dyn(imode,jmode) = dyn(imode,jmode) + work
        dyn_hub_bare(imode,jmode) = work
     ENDDO
  ENDDO
  !
  IF (iverbosity == 1) THEN
     IF (lgamma) THEN
        WRITE(stdout,'(5x,"d2N_bare nc Hermiticity residual =",1x,es12.5)') &
             d2n_herm_res
        WRITE(stdout,'(5x,"dyn_hub_bare nc Hermiticity residual =",1x,es12.5)') &
             dyn_herm_res
     ELSE
        DO jmode = 1, nmodes
           DO imode = 1, nmodes
              WRITE(stdout, &
                   '(5x,"SOC_U_Q_DYNMAT_BARE_CHECK",2(1x,i0),2(1x,es24.16))') &
                   imode, jmode, REAL(dyn_hub_bare(imode,jmode),DP), &
                   AIMAG(dyn_hub_bare(imode,jmode))
           ENDDO
        ENDDO
     ENDIF
     CALL tra_write_matrix_no_sym('dyn_hub_bare NOT SYMMETRIZED', &
                                  dyn_hub_bare, nat)
     CALL tra_write_matrix('dyn_hub_bare SYMMETRIZED', &
                           dyn_hub_bare, u, nat)
  ENDIF
  !
  DEALLOCATE(projector_d1)
  DEALLOCATE(projector_d2)
  DEALLOCATE(proj0)
  DEALLOCATE(proj_d1)
  DEALLOCATE(proj_d2)
  DEALLOCATE(band_weight)
  DEALLOCATE(d2ns_bare)
  DEALLOCATE(dynwrk)
  DEALLOCATE(ground_matrix)
  DEALLOCATE(d2n_matrix)
  DEALLOCATE(dn_j_matrix)
  DEALLOCATE(dv_i_matrix)
  !
  CALL stop_clock('dynmat_hub_bare')
  !
  RETURN
END SUBROUTINE dynmat_hub_bare_nc
