!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE sym_dns (ldim, npe, dns)
  !-----------------------------------------------------------------------
  !! DFPT+U: This routine symmetrizes the first order variation of 
  !! the occupation matrices dns due to the perturbation caused
  !! by the displacement of atoms.  
  !
  !! Written  by S. de Gironcoli and A. Floris.
  !! Modified by I. Timrov (01.10.2018).
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE ions_base,        ONLY : nat, ityp
  USE ldaU,             ONLY : Hubbard_l, is_hubbard
  USE lsda_mod,         ONLY : nspin
  USE lr_symm_base,     ONLY : nsymq, minus_q, irotmq, rtau, upert, upert_mq
  USE qpoint,           ONLY : xq
  USE symm_base,        ONLY : d1, d2, d3, irt
  USE noncollin_module, ONLY : noncolin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldim, npe
  COMPLEX(DP), INTENT(INOUT) :: dns(ldim,ldim,nspin,nat,npe)
  !
  ! ... local variables
  !
  INTEGER :: nt, ip, jp, na, nb, is, m1, m2, &
             m0, m00, isym, irot
  COMPLEX(DP), ALLOCATABLE :: dnr(:,:,:,:,:), dnraux(:,:,:,:,:)
  COMPLEX(DP) :: phase
  REAL(DP) :: arg
  !
  IF (noncolin) THEN
     CALL sym_dns_nc (ldim, npe, dns)
     RETURN
  ENDIF
  IF ((nsymq==1) .AND. (.NOT.minus_q)) RETURN
  !
  ! Initialization
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  ALLOCATE (dnraux(ldim,ldim,nspin,nat,npe))
  ALLOCATE (dnr(ldim,ldim,nspin,nat,npe))
  !  
  dnraux = (0.d0, 0.d0)
  dnr    = (0.d0, 0.d0)
  !
  ! Impose hermiticity of dns_{m1,m2, is,na,ip} for m1<->m2
  ! and put it in dnr zeroing dns
  !
  DO ip = 1, npe
     DO na = 1, nat  
        nt = ityp(na)
        DO is = 1, nspin  
           DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1  
                 dnr(m1, m2, is, na,ip) = 0.5d0 * ( dns(m1, m2, is, na, ip) &
                                                  + dns(m2, m1, is, na, ip) )
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  dns = (0.d0, 0.d0)
  !
  ! Symmetrize with -q if present (output overwritten on dnr)
  !
  IF (minus_q) THEN
     DO ip = 1, npe
        DO na = 1, nat  
           nt = ityp(na)  
           IF (is_hubbard(nt)) THEN 
              ! 
              arg = ( xq(1) * rtau(1,irotmq,na) + &
                      xq(2) * rtau(2,irotmq,na) + &
                      xq(3) * rtau(3,irotmq,na) ) * tpi
              phase = CMPLX (COS(arg), SIN(arg), kind=DP)
              !
              DO is = 1, nspin  
                 DO m1 = 1, 2 * Hubbard_l(nt) + 1  
                    DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                       nb = irt(irotmq, na)  
                       DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                          DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                             DO jp = 1, npe
                                IF (Hubbard_l(nt).EQ.0) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   dnr(m0,m00,is,nb,jp) * upert_mq(jp,ip) *    &
                                   phase
                                ELSE IF (Hubbard_l(nt).EQ.1) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   d1(m0 ,m1,irotmq) * d1(m00,m2,irotmq) *     &
                                   dnr(m0,m00,is,nb,jp) * upert_mq(jp,ip) *    &
                                   phase
                                ELSE IF (Hubbard_l(nt).EQ.2) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   d2(m0 ,m1,irotmq) * d2(m00,m2,irotmq) *     &
                                   dnr(m0,m00,is,nb,jp) * upert_mq(jp,ip) *    &
                                   phase
                                ELSE IF (Hubbard_l(nt).EQ.3) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   d3(m0 ,m1,irotmq) * d3(m00,m2,irotmq) *     &
                                   dnr(m0,m00,is,nb,jp) * upert_mq(jp,ip) *    &
                                   phase
                                ELSE
                                   CALL errore ('sym_dns', &
                                        'angular momentum not implemented', &
                                        ABS(Hubbard_l(nt)) )
                                ENDIF
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     dnr(:,:,:,:,:) = 0.5d0 * ( dnr(:,:,:,:,:) + CONJG(dns(:,:,:,:,:)) )
     dns = (0.d0, 0.d0)
  ENDIF
  ! 
  dnraux = (0.d0, 0.d0)
  !
  ! Symmetryze dnr -> dns
  !
  DO isym = 1, nsymq
     !
     dns = (0.d0, 0.d0)  
     !
     irot = isym
     !
     DO ip = 1, npe
        DO na = 1, nat  
           nt = ityp(na)  
           IF (is_hubbard(nt)) THEN  
              ! 
              arg = ( xq(1) * rtau(1,irot,na) + &
                      xq(2) * rtau(2,irot,na) + &
                      xq(3) * rtau(3,irot,na) ) * tpi
              phase = CMPLX (COS(arg), SIN(arg), kind=DP )
              !
              DO m1 = 1, 2 * Hubbard_l(nt) + 1  
                 DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                    nb = irt (irot, na)  
                    DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                       DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                          do jp=1, npe
                             IF (Hubbard_l(nt).EQ.0) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) +  &
                                dnr(m0,m00,:,nb,jp) * upert(jp,ip,irot) * &
                                phase 
                             ELSE IF (Hubbard_l(nt).EQ.1) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                d1(m0 ,m1,irot) * d1(m00,m2,irot) *       &
                                dnr(m0,m00,:,nb,jp) * upert(jp,ip,irot) * &
                                phase 
                             ELSE IF (Hubbard_l(nt).EQ.2) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                d2(m0 ,m1,irot) * d2(m00,m2,irot) *       &
                                dnr(m0,m00,:,nb,jp) * upert(jp,ip,irot) * &
                                phase 
                             ELSE IF (Hubbard_l(nt).EQ.3) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                d3(m0 ,m1,irot) * d3(m00,m2,irot) *       &
                                dnr(m0,m00,:,nb,jp) * upert(jp,ip,irot) * &
                                phase 
                             ELSE
                                CALL errore ('new_ns', &
                                     'angular momentum not implemented', &
                                     ABS(Hubbard_l(nt)) )
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     dnraux = dnraux + dns / nsymq
  ENDDO
  dns = dnraux
  ! 
  DEALLOCATE (dnr)
  DEALLOCATE (dnraux)
  ! 
  RETURN
  ! 
END SUBROUTINE sym_dns
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE sym_dns_nc (ldim, npe, dns)
  !-------------------------------------------------------------------
  !! Symmetrize the noncollinear Hubbard occupation response in the
  !! explicit PW ns_nc(:,:,npol*(is_bra-1)+is_ket,...) convention.
  !!
  !! The orbital-spin contraction follows new_ns_nc component by component:
  !!
  !!   dN'_(m1 m2,s1 s2) = D_(m1 m3) D_(m2 m4)
  !!                         U^*_(s1 s3) dN_(m3 m4,s3 s4) U_(s2 s4),
  !!
  !! with the atom, perturbation, and q-dependent phase transformations
  !! already used by the collinear sym_dns path (upert/upert_mq).  For an antiunitary
  !! magnetic operation the complete q-space response coefficient is
  !! complex conjugated and d_spin_ldau supplies the existing PW
  !! time-reversal spin mapping.
  !! At general q no Hermitian projection is made.  If a symmetry maps
  !! q to -q, only the physical combined-index adjoint relation is used.
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE ions_base,        ONLY : nat, ityp
  USE ldaU,             ONLY : Hubbard_l, is_hubbard, d_spin_ldau
  USE lsda_mod,         ONLY : nspin
  USE lr_symm_base,     ONLY : nsymq, minus_q, irotmq, rtau, upert, upert_mq
  USE qpoint,           ONLY : xq
  USE symm_base,        ONLY : t_rev
  USE noncollin_module, ONLY : npol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldim, npe
  COMPLEX(DP), INTENT(INOUT) :: dns(ldim,ldim,nspin,nat,npe)
  !
  INTEGER :: isym, ip, jp, na, nt, m1, m2
  INTEGER :: is_bra, is_ket, iblock, jblock
  COMPLEX(DP), ALLOCATABLE :: dnr(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dnrot(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dnraux(:,:,:,:,:)
  COMPLEX(DP) :: perturb_transform(npe,npe)
  COMPLEX(DP) :: phase_by_atom(nat)
  REAL(DP) :: arg
  !
  IF ((nsymq == 1) .AND. (.NOT.minus_q)) RETURN
  IF (npol /= 2 .OR. nspin /= npol*npol) &
       CALL errore('sym_dns_nc', 'spinor Hubbard response requires four blocks', 1)
  IF (.NOT.ALLOCATED(d_spin_ldau)) &
       CALL errore('sym_dns_nc', 'spin rotation matrices are not initialized', 1)
  !
  ALLOCATE(dnr(ldim,ldim,nspin,nat,npe))
  ALLOCATE(dnrot(ldim,ldim,nspin,nat,npe))
  ALLOCATE(dnraux(ldim,ldim,nspin,nat,npe))
  !
  ! Unlike the collinear response, dN(q) is not Hermitian at a general q.
  dnr = dns
  !
  ! If an operation sends q to -q, rotate first and then use
  ! dN(-q) = dN(q)^dagger in the complete orbital-spin space.
  IF (minus_q) THEN
     DO ip = 1, npe
        DO jp = 1, npe
           perturb_transform(jp,ip) = upert_mq(jp,ip)
        ENDDO
     ENDDO
     DO na = 1, nat
        arg = tpi * (xq(1)*rtau(1,irotmq,na) + &
                     xq(2)*rtau(2,irotmq,na) + &
                     xq(3)*rtau(3,irotmq,na))
        phase_by_atom(na) = CMPLX(COS(arg), SIN(arg), KIND=DP)
     ENDDO
     CALL rotate_dns_nc_symmetry(ldim, npe, irotmq, .FALSE., &
          perturb_transform, phase_by_atom, dnr, dnrot)
     DO ip = 1, npe
        DO na = 1, nat
           nt = ityp(na)
           IF (is_hubbard(nt)) THEN
              DO is_bra = 1, npol
                 DO is_ket = 1, npol
                    iblock = npol * (is_bra - 1) + is_ket
                    jblock = npol * (is_ket - 1) + is_bra
                    DO m1 = 1, 2*Hubbard_l(nt)+1
                       DO m2 = 1, 2*Hubbard_l(nt)+1
                          dnr(m1,m2,iblock,na,ip) = 0.5_DP * ( &
                               dnr(m1,m2,iblock,na,ip) + &
                               CONJG(dnrot(m2,m1,jblock,na,ip)) )
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !
  ! Project onto the invariant subspace of the magnetic little group of q.
  dnraux = (0.0_DP, 0.0_DP)
  DO isym = 1, nsymq
     DO ip = 1, npe
        DO jp = 1, npe
           perturb_transform(jp,ip) = upert(jp,ip,isym)
        ENDDO
     ENDDO
     DO na = 1, nat
        arg = tpi * (xq(1)*rtau(1,isym,na) + &
                     xq(2)*rtau(2,isym,na) + &
                     xq(3)*rtau(3,isym,na))
        phase_by_atom(na) = CMPLX(COS(arg), SIN(arg), KIND=DP)
     ENDDO
     CALL rotate_dns_nc_symmetry(ldim, npe, isym, t_rev(isym) == 1, &
          perturb_transform, phase_by_atom, dnr, dnrot)
     dnraux = dnraux + dnrot / DBLE(nsymq)
  ENDDO
  dns = dnraux
  !
  DEALLOCATE(dnr, dnrot, dnraux)
  RETURN
END SUBROUTINE sym_dns_nc
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE rotate_dns_nc_symmetry(ldim, npe, irot, antiunitary, &
                                  perturb_transform, phase_by_atom, &
                                  dns_in, dns_out)
  !-------------------------------------------------------------------
  !! Apply one spatial-orbital-spin symmetry operation to a set of
  !! noncollinear Hubbard occupation responses.  The caller owns the
  !! symmetry list, perturbation representation, q phase, and averaging.
  !! This routine owns only the algebra shared by phonon and electric-field
  !! response symmetrization.
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ityp
  USE ldaU,             ONLY : Hubbard_l, is_hubbard, d_spin_ldau
  USE lsda_mod,         ONLY : nspin
  USE symm_base,        ONLY : d1, d2, d3, irt
  USE noncollin_module, ONLY : npol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldim, npe, irot
  LOGICAL, INTENT(IN) :: antiunitary
  COMPLEX(DP), INTENT(IN) :: perturb_transform(npe,npe)
  COMPLEX(DP), INTENT(IN) :: phase_by_atom(nat)
  COMPLEX(DP), INTENT(IN) :: dns_in(ldim,ldim,nspin,nat,npe)
  COMPLEX(DP), INTENT(OUT) :: dns_out(ldim,ldim,nspin,nat,npe)
  !
  INTEGER :: ip, jp, na, nb, nt, m1, m2, m0, m00
  INTEGER :: out_bra, out_ket, in_bra, in_ket
  INTEGER :: out_block, in_block
  REAL(DP), TARGET :: d0(1,1)
  REAL(DP), POINTER :: dl(:,:)
  COMPLEX(DP) :: response_value, spin_factor
  !
  IF (npol /= 2 .OR. nspin /= npol*npol) &
       CALL errore('rotate_dns_nc_symmetry', &
                   'spinor Hubbard response requires four blocks', 1)
  IF (.NOT.ALLOCATED(d_spin_ldau)) &
       CALL errore('rotate_dns_nc_symmetry', &
                   'spin rotation matrices are not initialized', 1)
  !
  dns_out = (0.0_DP, 0.0_DP)
  d0(1,1) = 1.0_DP
  !
  DO na = 1, nat
     nt = ityp(na)
     IF (.NOT.is_hubbard(nt)) CYCLE
     !
     SELECT CASE (Hubbard_l(nt))
     CASE (0)
        dl => d0
     CASE (1)
        dl => d1(:,:,irot)
     CASE (2)
        dl => d2(:,:,irot)
     CASE (3)
        dl => d3(:,:,irot)
     CASE DEFAULT
        CALL errore('rotate_dns_nc_symmetry', &
                    'angular momentum not implemented', ABS(Hubbard_l(nt)))
     END SELECT
     !
     nb = irt(irot,na)
     DO ip = 1, npe
        DO jp = 1, npe
           DO out_bra = 1, npol
              DO out_ket = 1, npol
                 out_block = npol * (out_bra - 1) + out_ket
                 DO in_bra = 1, npol
                    DO in_ket = 1, npol
                       in_block = npol * (in_bra - 1) + in_ket
                       spin_factor = &
                            CONJG(d_spin_ldau(out_bra,in_bra,irot)) * &
                                  d_spin_ldau(out_ket,in_ket,irot)
                       DO m1 = 1, 2*Hubbard_l(nt)+1
                          DO m2 = 1, 2*Hubbard_l(nt)+1
                             DO m0 = 1, 2*Hubbard_l(nt)+1
                                DO m00 = 1, 2*Hubbard_l(nt)+1
                                   response_value = &
                                        dns_in(m0,m00,in_block,nb,jp) * &
                                        perturb_transform(jp,ip) * &
                                        phase_by_atom(na)
                                   IF (antiunitary) &
                                        response_value = CONJG(response_value)
                                   dns_out(m1,m2,out_block,na,ip) = &
                                        dns_out(m1,m2,out_block,na,ip) + &
                                        dl(m1,m0) * dl(m2,m00) * &
                                        spin_factor * response_value
                                ENDDO
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE rotate_dns_nc_symmetry
!---------------------------------------------------------------------
