!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dynmat_hub_scf (irr, nu_i0, nper)
  !---------------------------------------------------------------------
  !! DFPT+U: This routine adds to the dynamical matrix the scf+orth Hubbard terms.   
  !! It adds the scf derivative of the Hubbard energy (terms 1 and 2). Moreover,
  !! it adds 3 terms due to the orthogonality constraints in the USPP formalism (terms 3,4,5). 
  !! Note, the orthogonality terms 4 and 5 DO NOT follow simply from the 2nd derivative of 
  !! the Hubbard energy \(\text{E_U}\); they stem from the coupling of the \(-\epsilon\ dS\ d\lambda\) in 
  !! the USPP forces with the part of \(d\psi\ d\mu\) projected onto the occupied manifold 
  !! (the variable dpsi stems from the scf linear system and lives in the unoccupied manifold).  
  !! See comments in the source code for other details on the implemented terms.
  !
  ! Terms implemented:
  !   dyn_hub_scf (ipert, imode)  
  !  1) =  +\sum_{I,m1,m2,is} Hubbard_U(I) [0.5\delta_m1m2-ns(m1, m2, is, I)]* 
  !        \sum_{ibnd, k} wg(ibnd,k)[ <dpsi(ipert,ibnd,k+q, is)|dqsphi(imode,I,k+q,m1)><S\phi(I,k,m2)| psi(ibnd,k,is> + 
  !                                   <dpsi(ipert,ibnd,k+q, is)|S\phi(I,k+q,m1)><dmqsphi(imode,I,k,m2)| psi(ibnd,k,is> +  m1<=>m2]   
  !  2) = -\sum_{I,m1,m2,is} Hubbard_U(I) * CONJG(dnsscf(m1,m2,is,I,ipert)) * dnsbare(m2,m1,is,I,imode)  
  ! 
  ! Orthogonality terms present only in USPP:
  !  3) = -\sum_{I,m1,m2,is} Hubbard_U(I) [0.5\delta_m1m2-ns(m1, m2, is, I)]* 
  !        \sum_{ibnd, k} wg(ibnd,k)[ <dpsi_orth(ipert,ibnd,k+q, is)|dqsphi(imode,I,k+q,m2)><S\phi(I,k,m1)| psi(ibnd,k,is> + 
  !                                   <dpsi_orth(ipert,ibnd,k+q, is)|S\phi(I,k+q,m2)><dmqsphi(imode,I,k,m1)| psi(ibnd,k,is>] 
  !  4) = -\sum_{I,m1,m2,is} Hubbard_U(I) [0.5\delta_m1m2-ns(m1, m2, is, I)]* 
  !        \sum_{ibnd, k} wg(ibnd,k)[ <dpsi_orth(imode,ibnd,k+q, is)|dqsphi(ipert,I,k+q,m1)><S\phi(I,k,m2)| psi(ibnd,k,is> + 
  !                                   <dpsi_orth(imode,ibnd,k+q, is)|S\phi(I,k+q,m1)><dmqsphi(ipert,I,k,m2)| psi(ibnd,k,is>]  
  !  5) = -\sum_{I,m1,m2,is} Hubbard_U(I) * [CONJG(dnsscf(m1,m2,is,I,ipert) + CONJG(dnsbare(m2,m1,is,I,ipert))] 
  !                                       * dnsorth(m1,m2,is,I,imode) 
  !
  ! Note: dnsscf includes the orthogonality part (dnsorth)
  !
  !! Written  by A. Floris.  
  !! Modified by I. Timrov (01.10.2018).
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ldaU,          ONLY : Hubbard_l, is_hubbard, Hubbard_J0
  USE ldaU_lr,       ONLY : dnsscf, effU
  USE ldaU_ph,       ONLY : dnsbare, dnsbare_all_modes, dnsorth_cart
  USE lsda_mod,      ONLY : lsda, current_spin, isk, nspin
  USE modes,         ONLY : u, nmodes
  USE dynmat,        ONLY : dyn, dyn_rec, dyn_hub_scf
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE eqv,           ONLY : dpsi
  USE wvfct,         ONLY : npwx, nbnd
  USE control_lr,    ONLY : lgamma, rec_code_read
  USE units_lr,      ONLY : iuwfc, lrwfc, iudwf, lrdwf
  USE wavefunctions, ONLY : evc
  USE klist,         ONLY : wk, lgauss, ltetra, ngk, igk_k
  USE uspp,          ONLY : okvan
  USE control_flags, ONLY : iverbosity
  USE io_global,     ONLY : stdout
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp_pools,      ONLY : inter_pool_comm       
  USE mp,            ONLY : mp_sum
  USE buffers,       ONLY : get_buffer
  USE noncollin_module, ONLY : noncolin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: irr
  !! the irreducible representation
  INTEGER, INTENT(IN) :: nu_i0
  !! the initial position of the mode
  INTEGER, INTENT(IN) :: nper
  !! number of perturbations in the current irreducible representation
  !
  ! ... local variables
  !
  INTEGER :: icar, jcar, na, nap, nah, ihubst1, ihubst2, nt, counter, n, l, is,    & 
             na_icar, nap_jcar, m1, m2, imode, jmode, ik, ikk, ikq, ipert, nrec1,  &
             ibnd, ig, isi, op_spin, npw, npwq
  COMPLEX(DP), ALLOCATABLE :: dyn1 (:,:), dyn_orth(:,:), term(:,:), term_aux(:,:), &
                              dpsi_orth_cart(:,:,:,:), dyn_orth_cart(:,:),         &
                              dvqhbar(:,:,:,:), dvqhbar_orth(:,:,:,:),             &
                              dvqhbar_orth_lm(:,:,:,:), aux1(:,:,:), dyn1_test(:,:) 
  REAL(DP), ALLOCATABLE :: wgg(:,:,:)
  COMPLEX(DP) :: dvi, prj, prj_orth
  LOGICAL :: lmetq0   ! .true. if q=0 for a metal
  ! 
  IF (noncolin) THEN
     CALL dynmat_hub_scf_nc(irr, nu_i0, nper)
     RETURN
  ENDIF
  !
  CALL start_clock( 'dynmat_hub_scf' )
  !
  ALLOCATE (dyn1(nper,nmodes))
  ALLOCATE (dyn_orth(nper,nmodes))
  ALLOCATE (term(nat,3))
  ALLOCATE (term_aux(nper,nmodes))
  ALLOCATE (dpsi_orth_cart(npwx,nbnd,3,nat))
  ALLOCATE (dyn_orth_cart(3*nat,3*nat))
  ALLOCATE (dvqhbar(npwx,nbnd,3,nat))
  ALLOCATE (aux1(npwx,nbnd,nmodes))
  ALLOCATE (dyn1_test(nmodes,nmodes))
  ALLOCATE (wgg(nbnd,nbnd,nksq)) 
  ALLOCATE (dvqhbar_orth(npwx,nbnd,3,nat))
  ALLOCATE (dvqhbar_orth_lm(npwx,nbnd,3,nat))
  ! 
  dyn1          = (0.d0, 0.d0)  
  dyn_orth_cart = (0.d0, 0.d0)  
  dyn1_test     = (0.d0, 0.d0)  
  ! 
  lmetq0 = (lgauss .OR. ltetra) .AND. lgamma
  !
  ! USPP: compute the weights as in square bracket 
  ! of Eq. (27) of A. Dal Corso PRB 64, 235118 (2001).
  ! Or see Eq. (B19) in the same reference.
  !
  IF (okvan) CALL compute_weight (wgg)
  ! 
  ! If recover=.true. recompute dnsscf for the current irr 
  ! (e.g. when convt=.true. in solve_linter but the code is 
  ! interrupted before the call of this routine)
  !
  IF (rec_code_read==10) THEN
     !WRITE(stdout,*) 'rec_code_read', rec_code_read
     CALL dnsq_scf (nper, lmetq0)
  ENDIF
  !
  ! We need a sum over all k points 
  !
  DO ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     !
     IF (lsda) THEN
        current_spin = isk(ikk)
        IF (current_spin .eq. 1) THEN
           op_spin = 2
        ELSE
           op_spin = 1
        END IF
     ELSE        
        op_spin = 1
     ENDIF
     !
     ! Read unperturbed KS wavefuctions psi(k)     
     !
     IF (nksq.GT.1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
     !
     ! Calculate d^{na,icar}V_tilde(q) |psi(ibnd,k)> =  
     ! = |ket> Bloch vector at k+q for all cartesian coordinates (na,icar)
     !
     CALL dvqhub_barepsi_us2 (ik, dvqhbar, dvqhbar_orth, dvqhbar_orth_lm)
     !
     ! Compute terms 3 and 4.
     ! Compute the orthogonality term in cartesian coordinates.
     ! Hubbard-potential analogue to the two terms in 
     ! Eq. (42) of A. Dal Corso PRB 64, 235118 (2001).
     !
     IF (okvan)  THEN
        !
        ! Calculate dpsi_orth_cart which is the valence component 
        ! of the response KS wavefunction (non-zero only in the USPP case).
        !
        CALL dpsi_orth (ik, wgg, dpsi_orth_cart)
        !
        DO na = 1, nat
           DO icar = 1, 3
              na_icar = 3*(na-1)+icar
              DO nap = 1, nat
                 DO jcar = 1, 3
                    nap_jcar = 3*(nap-1)+jcar
                    DO ibnd = 1, nbnd
                       !
                       !          from E_Hub alone
                       prj_orth = dot_product (dpsi_orth_cart(1:npwq,ibnd,jcar,nap), &
                                                 dvqhbar_orth(1:npwq,ibnd,icar,na))+ & 
                       !          extra term      
                                  dot_product (dvqhbar_orth_lm(1:npwq,ibnd,jcar,nap),&
                                               dpsi_orth_cart( 1:npwq,ibnd,icar,na))
                       ! 
                       CALL mp_sum(prj_orth, intra_bgrp_comm) 
                       ! 
                       dyn_orth_cart (nap_jcar,na_icar) = dyn_orth_cart (nap_jcar,na_icar)  &
                                                          - wk(ikk) * prj_orth 
                       !
                    ENDDO 
                 ENDDO 
              ENDDO 
           ENDDO
        ENDDO
        ! 
     ENDIF
     !
     ! Compute term 1.
     ! Transform dvqhbar from the Cartesian to the pattern basis u. 
     ! aux1 is dvqhbar in the pattern basis.
     ! 
     aux1 = (0.d0, 0.d0)
     !
     DO imode = 1, nmodes 
        DO na = 1, nat
           DO icar = 1, 3 
              na_icar = 3*(na-1)+icar
              DO ibnd = 1, nbnd           
                 DO ig = 1, npwq
                    aux1(ig,ibnd,imode) = aux1(ig,ibnd,imode) + & 
                         dvqhbar(ig,ibnd,icar,na) * u(na_icar,imode)
                 ENDDO 
              ENDDO 
           ENDDO 
        ENDDO 
     ENDDO 
     !
     DO ipert = 1, nper
        !                
        nrec1 = (ipert-1)*nksq+ik
        !
        ! Read the response KS wave functions dpsi from file 
        !
        CALL get_buffer (dpsi, lrdwf, iudwf, nrec1)
        ! 
        DO imode = 1, nmodes 
           ! 
           ! Calculate prj = <d^(ipert)dpsi_k+q | d^{imode}V_tilde(q) | psi(ibnd,k)>
           !
           DO ibnd = 1, nbnd
              !
              prj = dot_product  (dpsi(1:npwq,ibnd), aux1(1:npwq,ibnd,imode))
                          
              CALL mp_sum (prj, intra_bgrp_comm) 
              !
              ! dyn1 is a sum over all the k and ibnd
              !
              dyn1(ipert,imode) =  dyn1(ipert,imode) + wk(ikk)*prj 
              !
           ENDDO 
           !
        ENDDO 
        ! 
     ENDDO 
     !
  ENDDO ! ik
  !
  CALL mp_sum (dyn1, inter_pool_comm) 
  CALL mp_sum (dyn_orth_cart, inter_pool_comm)
  !
  ! Compute terms 2 and 5 (and the equivalent for Hubbard_J0)
  !
  DO ipert = 1, nper
     !
     term = (0.d0, 0.d0)   
     !    
     DO na = 1, nat
        !
        DO icar = 1, 3
           !
           DO nah = 1, nat
              !
              nt = ityp(nah)
              !
              ! For effU = Hubbard_U - Hubbard_J0
              !
              IF (is_hubbard(nt)) THEN
                 !
                 dvi = (0.d0, 0.d0)
                 !           
                 DO is = 1, nspin
                    !
                    DO m1 = 1, 2*Hubbard_l(nt)+1
                       !
                       DO m2 = 1, 2*Hubbard_l(nt)+1
                          ! 
                          ! Term 2
                          !
                          dvi = dvi + CONJG(dnsscf(m1,m2,is,nah,ipert)) * &
                                      dnsbare(m2,m1,is,nah,icar,na)   
                          ! 
                          ! Term 5
                          !
                          IF (okvan) dvi = dvi + CONJG( dnsscf(m1,m2,is,nah,ipert) + &
                                                        dnsbare_all_modes(m1,m2,is,nah,nu_i0+ipert) ) &
                                                      * dnsorth_cart(m1,m2,is,nah,icar,na)   
                          !
                       ENDDO
                       !
                    ENDDO
                    !
                 ENDDO
                 !                    
                 ! term is implicitly defined for each ipert
                 !
                 term(na,icar) = term(na,icar) - dvi*effU(nt)  
                 ! 
              ENDIF
              !
              ! For Hubbard_J0
              !
              IF (Hubbard_J0(nt).NE.0.d0) THEN
                 !
                 dvi = (0.d0, 0.d0)   
                 !        
                 DO is = 1, nspin
                    !                    
                    IF ((is==1).AND.(nspin==2)) THEN
                       isi = 2
                    ELSE 
                       isi = 1
                    ENDIF
                    !
                    DO m1 = 1, 2*Hubbard_l(nt)+1
                       !
                       DO m2 = 1, 2*Hubbard_l(nt)+1
                          !
                          ! Term 2
                          !
                          dvi = dvi - CONJG(dnsscf(m1,m2,isi,nah,ipert)) * & ! sign change
                                      dnsbare(m2,m1,is,nah,icar,na)   
                          !
                          ! Term 5
                          !
                          IF (okvan) dvi = dvi - CONJG( dnsscf(m1,m2,isi,nah,ipert) + & ! sign change
                                                        dnsbare_all_modes(m1,m2,isi,nah,nu_i0+ipert) ) &
                                                      * dnsorth_cart(m1,m2,is,nah,icar,na)  
                          ! 
                       ENDDO
                       !
                    ENDDO
                    !
                 ENDDO
                 !                    
                 ! term is implicitely defined for each ipert
                 !
                 term(na,icar) = term(na,icar) - dvi*Hubbard_J0(nt)  
                 !
              ENDIF
              !
           ENDDO ! nah
           !
        ENDDO ! icar
        !
     ENDDO ! na
     !
     ! transform term from the Cartesian to the pattern basis u. 
     ! The result is stored in term_aux.
     ! 
     DO imode = 1, nmodes
        term_aux(ipert,imode) = (0.d0, 0.d0)
        DO na = 1, nat
           DO icar = 1, 3
              na_icar = 3*(na-1)+icar
              term_aux(ipert,imode) = term_aux(ipert,imode) + &
                                      term(na,icar) * u(na_icar,imode)         
           ENDDO
        ENDDO
     ENDDO
     !
  ENDDO ! ipert
  !
  ! Factor of 2 taking into account that
  ! if nspin=1 the spin sum above does not occur   
  !
  IF (nspin==1) term_aux = 2.d0 * term_aux
  !
  ! Add the contribution to the dyn_hub_scf matrix
  !
  dyn1 = dyn1 + term_aux    
  !  
  dyn_orth = (0.d0, 0.d0)
  !
  ! Adding the orthogonality terms 3 and 4, 
  ! which were computed above (and summed up over k), to dyn1.  
  !
  IF (okvan) THEN
     DO ipert = 1, nper
        DO imode = 1, nmodes
           DO na_icar = 1, 3*nat
              DO nap_jcar = 1, 3*nat 
                 dyn_orth(ipert,imode) = dyn_orth(ipert,imode) &
                                         + u(na_icar,imode)    &  
                                         * dyn_orth_cart(nap_jcar,na_icar) &
                                         * CONJG(u(nap_jcar,nu_i0+ipert)) 
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     dyn1 = dyn1 + dyn_orth
  ENDIF
  ! 
  DO imode = 1, nmodes
     DO ipert = 1, nper
        !
        jmode = nu_i0 + ipert
        !
        dyn_hub_scf(jmode,imode) = dyn1(ipert,imode)
        !
        dyn_rec(jmode,imode) = dyn_rec(jmode,imode) + dyn1(ipert,imode)
        ! 
        dyn1_test(jmode,imode) = dyn1 (ipert, imode)
        !
     ENDDO
  ENDDO
  !
  ! Write the symmetrized dyn_hub_scf (upgraded for the new ipert) and total 
  ! dynamical matrix dyn in  cartesian coordinates.
  ! dyn1 is instead the contribution of the particular irrep only.
  !
  IF (iverbosity==1) THEN
     CALL tra_write_matrix_no_sym('dyn1 NOT SYMMETRIZED',dyn1_test,nat)
     CALL tra_write_matrix('dyn1 SYMMETRIZED',dyn1_test,u,nat)
     CALL tra_write_matrix('dyn',dyn,u,nat)
  ENDIF
  !
  ! Add Hubbard terms in the dynamical matrix to the total
  ! dynamical matrix.
  !
  DO imode = 1, nmodes
     DO ipert = 1, nper
        jmode = nu_i0 + ipert
        dyn(jmode,imode) = dyn(jmode,imode) + dyn1(ipert,imode)
     ENDDO
  ENDDO
  !
  DEALLOCATE (dyn1)
  DEALLOCATE (dyn_orth)
  DEALLOCATE (term)
  DEALLOCATE (term_aux)
  DEALLOCATE (dpsi_orth_cart)
  DEALLOCATE (dyn_orth_cart)
  DEALLOCATE (dvqhbar)
  DEALLOCATE (aux1)
  DEALLOCATE (dyn1_test)
  DEALLOCATE (wgg)
  DEALLOCATE (dvqhbar_orth)
  DEALLOCATE (dvqhbar_orth_lm)
  !   
  CALL stop_clock ( 'dynmat_hub_scf' )
  !
  RETURN
  !
END SUBROUTINE dynmat_hub_scf

!-----------------------------------------------------------------------
SUBROUTINE dynmat_hub_scf_nc (irr, nu_i0, nper)
  !---------------------------------------------------------------------
  !! Norm-conserving, fixed-Dudarev noncollinear Hubbard contribution
  !! to the self-consistent dynamical matrix:
  !!
  !!   term 1 = V_U : delta_i delta_j N_cross
  !!   term 2 = delta_i V_U_bare : delta_j N_scf .
  !!
  !! The cross occupation response is reconstructed from both magnetic
  !! Sternheimer records.  All contractions retain the complete
  !! orbital-spin matrix.
  !
  USE kinds,             ONLY : DP
  USE ions_base,         ONLY : nat, ityp, ntyp => nsp
  USE ldaU,              ONLY : Hubbard_lmax, Hubbard_l, offsetU, &
                                is_hubbard, Hubbard_J0, Hubbard_alpha, &
                                Hubbard_beta, Hubbard_projectors, &
                                lda_plus_u_kind, nwfcU
  USE ldaU_lr,           ONLY : dnsscf, effU, swfcatomk, swfcatomkpq, &
                                pack_hubbard_blocks_nc, build_dvhub_nc, &
                                accumulate_ns_nc_transition, &
                                ns_nc_hermiticity_residual
  USE ldaU_ph,           ONLY : dnsbare, wfcatomk, wfcatomkpq
  USE noncollin_module,  ONLY : domag, npol
  USE control_lr,        ONLY : lgamma, nbnd_occ, rec_code_read
  USE uspp,              ONLY : okvan
  USE qpoint,            ONLY : nksq, ikks, ikqs, ikmks
  USE wvfct,             ONLY : npwx, nbnd
  USE klist,             ONLY : wk, ngk, igk_k, lgauss, ltetra
  USE wavefunctions,     ONLY : evc
  USE eqv,               ONLY : dpsi
  USE units_lr,          ONLY : iuwfc, lrwfc, iudwf, lrdwf, &
                                iuatwfc, iuatswfc
  USE io_files,          ONLY : nwordwfcU
  USE buffers,           ONLY : get_buffer
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp_pools,          ONLY : inter_pool_comm
  USE mp,                ONLY : mp_sum
  USE scf,               ONLY : v
  USE modes,             ONLY : u, nmodes
  USE dynmat,            ONLY : dyn, dyn_rec, dyn_hub_scf
  USE control_flags,     ONLY : iverbosity
  USE io_global,         ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: irr
  INTEGER, INTENT(IN) :: nu_i0
  INTEGER, INTENT(IN) :: nper
  !
  INTEGER :: ik, ikk, ikq, ikmk, npw, npwq, isolv, ipert, nrec
  INTEGER :: na, nah, nt, icart, is_row, is_col, m, iwfc
  INTEGER :: row_first, na_icart, ldim_atom, max_ldim
  INTEGER :: matrix_size, max_matrix_size, imode, jmode
  REAL(DP) :: cross_herm_res
  REAL(DP), ALLOCATABLE :: band_weight(:)
  COMPLEX(DP) :: work
  COMPLEX(DP), EXTERNAL :: ZDOTU, ZDOTC
  COMPLEX(DP), ALLOCATABLE :: dprojector_k(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dprojector_kpq(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: proj0(:,:), proj_d1(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: proj_resp(:,:), proj_d1_resp(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: cross_ns(:,:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: ground_matrix(:,:), cross_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: scf_matrix(:,:), dv_bare_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: cart_terms(:,:), dyn1(:,:)
  COMPLEX(DP), ALLOCATABLE :: dyn1_test(:,:)
  LOGICAL :: lmetq0
  !
  CALL start_clock('dynmat_hub_scf')
  !
  IF (.NOT.domag) CALL errore('dynmat_hub_scf_nc', &
                         'noncollinear Hubbard response requires two magnetic solutions', 1)
  IF (okvan) CALL errore('dynmat_hub_scf_nc', &
                         'noncollinear Hubbard response requires NCPP', 1)
  IF (Hubbard_projectors /= 'atomic') CALL errore('dynmat_hub_scf_nc', &
       'noncollinear Hubbard response requires atomic projectors', 1)
  IF (lda_plus_u_kind /= 0) CALL errore('dynmat_hub_scf_nc', &
                                        'only Dudarev DFT+U is supported', 1)
  IF (npol /= 2) CALL errore('dynmat_hub_scf_nc', &
                             'noncollinear Hubbard response requires npol=2', 1)
  IF (.NOT.ALLOCATED(v%ns_nc)) CALL errore('dynmat_hub_scf_nc', &
                                           'PW Hubbard potential is unavailable', 1)
  DO nt = 1, ntyp
     IF (is_hubbard(nt) .AND. Hubbard_J0(nt) /= 0.0_DP) &
        CALL errore('dynmat_hub_scf_nc', 'Hubbard J0 is outside the MVP', 1)
     IF (is_hubbard(nt) .AND. Hubbard_alpha(nt) /= 0.0_DP) &
        CALL errore('dynmat_hub_scf_nc', 'Hubbard alpha is outside the MVP', 1)
     IF (is_hubbard(nt) .AND. Hubbard_beta(nt) /= 0.0_DP) &
        CALL errore('dynmat_hub_scf_nc', 'Hubbard beta is outside the MVP', 1)
  ENDDO
  !
  lmetq0 = (lgauss .OR. ltetra) .AND. lgamma
  IF (rec_code_read == 10) &
     CALL dnsq_scf(nper, lmetq0)
  !
  max_ldim = 2 * Hubbard_lmax + 1
  max_matrix_size = max_ldim * npol
  ALLOCATE(dprojector_k(npwx*npol,nwfcU,3))
  ALLOCATE(dprojector_kpq(npwx*npol,nwfcU,3))
  ALLOCATE(proj0(nwfcU,nbnd))
  ALLOCATE(proj_d1(nwfcU,nbnd,3))
  ALLOCATE(proj_resp(nwfcU,nbnd))
  ALLOCATE(proj_d1_resp(nwfcU,nbnd,3))
  ALLOCATE(band_weight(nbnd))
  ALLOCATE(cross_ns(max_ldim,max_ldim,npol*npol,nat,3,nper))
  ALLOCATE(ground_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(cross_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(scf_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(dv_bare_matrix(max_matrix_size,max_matrix_size))
  ALLOCATE(cart_terms(nper,3*nat))
  ALLOCATE(dyn1(nper,nmodes))
  ALLOCATE(dyn1_test(nmodes,nmodes))
  !
  cross_ns = (0.0_DP, 0.0_DP)
  !
  DO ik = 1, nksq
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
                                 dprojector_k(row_first,iwfc,icart))
                       IF (.NOT.lgamma) &
                          CALL dwfc(npwq, igk_k(1,ikq), ikq, icart, &
                                    wfcatomkpq(row_first,iwfc), &
                                    dprojector_kpq(row_first,iwfc,icart))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     IF (lgamma) dprojector_kpq = dprojector_k
     !
     band_weight = 0.0_DP
     band_weight(1:nbnd_occ(ikk)) = wk(ikk)
     !
     DO isolv = 1, 2
        IF (isolv == 1) THEN
           ikmk = ikk
        ELSE
           ikmk = ikmks(ik)
        ENDIF
        CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
        !
        CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                   (1.0_DP,0.0_DP), swfcatomk, npwx*npol, evc, &
                   npwx*npol, (0.0_DP,0.0_DP), proj0, nwfcU)
        DO icart = 1, 3
           CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                      (1.0_DP,0.0_DP), dprojector_k(:,:,icart), npwx*npol, &
                      evc, npwx*npol, (0.0_DP,0.0_DP), &
                      proj_d1(:,:,icart), nwfcU)
        ENDDO
        CALL mp_sum(proj0, intra_bgrp_comm)
        CALL mp_sum(proj_d1, intra_bgrp_comm)
        !
        DO ipert = 1, nper
           nrec = (isolv - 1) * nper * nksq + &
                  (ipert - 1) * nksq + ik
           CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
           CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                      (1.0_DP,0.0_DP), swfcatomkpq, npwx*npol, dpsi, &
                      npwx*npol, (0.0_DP,0.0_DP), proj_resp, nwfcU)
           DO icart = 1, 3
              CALL ZGEMM('C', 'N', nwfcU, nbnd, npwx*npol, &
                         (1.0_DP,0.0_DP), dprojector_kpq(:,:,icart), &
                         npwx*npol, dpsi, npwx*npol, &
                         (0.0_DP,0.0_DP), proj_d1_resp(:,:,icart), nwfcU)
           ENDDO
           CALL mp_sum(proj_resp, intra_bgrp_comm)
           CALL mp_sum(proj_d1_resp, intra_bgrp_comm)
           !
           DO nah = 1, nat
              nt = ityp(nah)
              IF (is_hubbard(nt)) THEN
                 ldim_atom = 2 * Hubbard_l(nt) + 1
                 DO icart = 1, 3
                    CALL accumulate_ns_nc_transition(ldim_atom, &
                         offsetU(nah), npol, nbnd_occ(ikk), band_weight, &
                         proj_d1(:,:,icart), proj_resp, &
                         cross_ns(:,:,:,nah,icart,ipert), isolv == 2)
                    CALL accumulate_ns_nc_transition(ldim_atom, &
                         offsetU(nah), npol, nbnd_occ(ikk), band_weight, &
                         proj0, proj_d1_resp(:,:,icart), &
                         cross_ns(:,:,:,nah,icart,ipert), isolv == 2)
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  CALL mp_sum(cross_ns, inter_pool_comm)
  !
  cross_herm_res = 0.0_DP
  DO ipert = 1, nper
     DO nah = 1, nat
        nt = ityp(nah)
        IF (is_hubbard(nt)) THEN
           ldim_atom = 2 * Hubbard_l(nt) + 1
           DO icart = 1, 3
              cross_herm_res = MAX(cross_herm_res, &
                   ns_nc_hermiticity_residual(ldim_atom, &
                   cross_ns(:,:,:,nah,icart,ipert)))
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  cart_terms = (0.0_DP, 0.0_DP)
  DO ipert = 1, nper
     DO nah = 1, nat
        nt = ityp(nah)
        IF (is_hubbard(nt)) THEN
           ldim_atom = 2 * Hubbard_l(nt) + 1
           matrix_size = ldim_atom * npol
           CALL pack_hubbard_blocks_nc(ldim_atom, npol, &
                                        v%ns_nc(:,:,:,nah), ground_matrix)
           CALL pack_hubbard_blocks_nc(ldim_atom, npol, &
                                        dnsscf(:,:,:,nah,ipert), scf_matrix)
           !
           DO icart = 1, 3
              na_icart = 3 * (nah - 1) + icart
              CALL pack_hubbard_blocks_nc(ldim_atom, npol, &
                   cross_ns(:,:,:,nah,icart,ipert), cross_matrix)
              ! cross_matrix is linear in the response ket.  Conjugate the
              ! complete contraction to recover <dpsi|dV_U|psi>.
              work = CONJG(ZDOTU(matrix_size*matrix_size, ground_matrix, 1, &
                                 cross_matrix, 1))
              cart_terms(ipert,na_icart) = &
                   cart_terms(ipert,na_icart) + work
           ENDDO
           !
           DO na = 1, nat
              DO icart = 1, 3
                 na_icart = 3 * (na - 1) + icart
                 CALL build_dvhub_nc(ldim_atom, npol, effU(nt), &
                      dnsbare(:,:,:,nah,icart,na), dv_bare_matrix)
                 IF (lgamma) THEN
                    work = ZDOTU(matrix_size*matrix_size, &
                                 dv_bare_matrix, 1, scf_matrix, 1)
                 ELSE
                    ! Reproduce the collinear general-q order
                    ! CONJG(dN_scf(q))*dN_bare(q) in the full
                    ! orbital-spin representation.
                    work = ZDOTC(matrix_size*matrix_size, scf_matrix, 1, &
                                 dv_bare_matrix, 1)
                 ENDIF
                 cart_terms(ipert,na_icart) = &
                      cart_terms(ipert,na_icart) + work
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  CALL ZGEMM('N', 'N', nper, nmodes, 3*nat, &
             (1.0_DP,0.0_DP), cart_terms, nper, u, 3*nat, &
             (0.0_DP,0.0_DP), dyn1, nper)
  !
  dyn1_test = (0.0_DP, 0.0_DP)
  DO imode = 1, nmodes
     DO ipert = 1, nper
        jmode = nu_i0 + ipert
        dyn_hub_scf(jmode,imode) = dyn1(ipert,imode)
        dyn_rec(jmode,imode) = dyn_rec(jmode,imode) + dyn1(ipert,imode)
        dyn(jmode,imode) = dyn(jmode,imode) + dyn1(ipert,imode)
        dyn1_test(jmode,imode) = dyn1(ipert,imode)
     ENDDO
  ENDDO
  !
  IF (iverbosity == 1) THEN
     IF (lgamma) THEN
        WRITE(stdout,'(5x,"cross dnsscf nc Hermiticity residual =",1x,es12.5)') &
             cross_herm_res
     ELSE
        DO imode = 1, nmodes
           DO ipert = 1, nper
              jmode = nu_i0 + ipert
              WRITE(stdout, &
                   '(5x,"SOC_U_Q_DYNMAT_SCF_CHECK",2(1x,i0),2(1x,es24.16))') &
                   jmode, imode, REAL(dyn1(ipert,imode),DP), &
                   AIMAG(dyn1(ipert,imode))
           ENDDO
        ENDDO
     ENDIF
     CALL tra_write_matrix_no_sym('dyn1 NC NOT SYMMETRIZED',dyn1_test,nat)
     CALL tra_write_matrix('dyn1 NC SYMMETRIZED',dyn1_test,u,nat)
     CALL tra_write_matrix('dyn',dyn,u,nat)
  ENDIF
  !
  DEALLOCATE(dprojector_k)
  DEALLOCATE(dprojector_kpq)
  DEALLOCATE(proj0)
  DEALLOCATE(proj_d1)
  DEALLOCATE(proj_resp)
  DEALLOCATE(proj_d1_resp)
  DEALLOCATE(band_weight)
  DEALLOCATE(cross_ns)
  DEALLOCATE(ground_matrix)
  DEALLOCATE(cross_matrix)
  DEALLOCATE(scf_matrix)
  DEALLOCATE(dv_bare_matrix)
  DEALLOCATE(cart_terms)
  DEALLOCATE(dyn1)
  DEALLOCATE(dyn1_test)
  !
  CALL stop_clock('dynmat_hub_scf')
  !
  RETURN
  !
END SUBROUTINE dynmat_hub_scf_nc
