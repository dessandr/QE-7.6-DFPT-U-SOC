!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from PW/src/s_psi.f90/s_psi_
!----------------------------------------------------------------------------
SUBROUTINE pw_s_psi( lda, n, m, vkb, becp, psi, spsi )
  !----------------------------------------------------------------------------
  !! This routine applies the S matrix to m wavefunctions psi and puts 
  !! the results in spsi.
  !! Requires the products of psi with all beta functions in array 
  !! becp(nkb,m) (calculated in h_psi or by calbec).
  !
  USE kinds,            ONLY: DP
  USE becmod,           ONLY: bec_type !becp
  USE uspp,             ONLY: nkb, okvan, qq_at, qq_so, ofsbeta !vkb, 
  USE uspp_param,       ONLY: upf, nh, nhm
  USE ions_base,        ONLY: nat, nsp, ityp
  USE control_flags,    ONLY: gamma_only 
  USE noncollin_module, ONLY: npol, noncolin, lspinorb
  USE realus,           ONLY: real_space
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  !! the m wavefunctions
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  !! S matrix dot wavefunctions psi
  !
  !jjz
  COMPLEX(DP), INTENT(IN) :: vkb(lda, nkb)
  TYPE(bec_type), INTENT(IN) :: becp !becp(nkb, m)
  !
  !jjz: drop support for gamma_only and real_space for simplicity
  if(gamma_only) call errore('pw_s_psi',"Gamma-only calculation not supported",1)
  if(real_space) call errore('pw_s_psi',"real_space not supported",1)
  !
  ! ... initialize  spsi
  !
  CALL threaded_memcpy( spsi, psi, lda*npol*m*2 )
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
  !
  CALL start_clock( 's_psi' )  
  !
  ! ... The product with the beta functions
  !
  ! QE 7.6: the gamma_only and real_space branches of PW s_psi, rejected
  ! above, are removed (bec_type has no band distribution any more)
  IF ( noncolin ) THEN
     !
     CALL s_psi_nc()
     !
  ELSE
     !
     CALL s_psi_k()
     !
  ENDIF
  !
  CALL stop_clock( 's_psi' )
  !
  !
  RETURN
  !
  CONTAINS
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_k()
       !-----------------------------------------------------------------------
       !! k-points version of \(\textrm{s_psi}\) routine.
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ierr
       ! counters
       COMPLEX(DP), ALLOCATABLE :: ps(:,:), qqc(:,:)
       ! ps = product vkb and psi ; qqc = complex version of qq
       !
       ALLOCATE( ps( nkb, m ), STAT=ierr )
       !
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_k ', ' cannot allocate memory (ps) ', ABS(ierr) )
       !
       DO nt = 1, nsp
          !
          IF ( upf(nt)%tvanp ) THEN
             ! qq is real:  copy it into a complex variable to perform
             ! a zgemm - simple but sub-optimal solution
             ALLOCATE( qqc(nh(nt),nh(nt)) )
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   qqc(:,:) = CMPLX ( qq_at(1:nh(nt),1:nh(nt),na), 0.0_DP, KIND=DP )
                   CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_DP,0.0_DP), &
                        qqc, nh(nt), becp%k(ofsbeta(na)+1,1), nkb, &
                        (0.0_DP,0.0_DP), ps(ofsbeta(na)+1,1), nkb )
                ENDIF
             ENDDO
             DEALLOCATE( qqc )
             !
          ELSE
             !
             IF (nh(nt)>0) THEN
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      ps(ofsbeta(na)+1:ofsbeta(na)+nh(nt),1:m) = (0.0_DP,0.0_DP)
                   ENDIF
                ENDDO
             ENDIF
             !
          ENDIF
          !
       ENDDO
       !
       IF ( m == 1 ) THEN
          !
          CALL ZGEMV( 'N', n, nkb, ( 1.D0, 0.D0 ), vkb, &
                      lda, ps, 1, ( 1.D0, 0.D0 ), spsi, 1 )
          !
       ELSE
          !
          CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ), vkb, &
                      lda, ps, nkb, ( 1.D0, 0.D0 ), spsi, lda )
          !
       ENDIF
       !
       DEALLOCATE( ps )
       !
       !
       RETURN
       !
     END SUBROUTINE s_psi_k     
     !
     !
     !-----------------------------------------------------------------------
       SUBROUTINE s_psi_nc ( )
       !-----------------------------------------------------------------------
       !! k-points noncolinear/spinorbit version of \(\textrm{s_psi}\) routine.
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ibnd, ipol, ierr
       ! counters
       COMPLEX (DP), ALLOCATABLE :: ps(:,:,:)
       ! the product vkb and psi
       !
       ALLOCATE( ps(nkb,npol,m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' s_psi_nc ', ' cannot allocate memory (ps) ', ABS(ierr) )
       ps(:,:,:) = (0.D0,0.D0)
       !
       DO nt = 1, nsp
          !
          IF ( upf(nt)%tvanp ) THEN
             !
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ofsbeta(na) + ih
                      DO jh = 1, nh(nt)
                         jkb = ofsbeta(na) + jh
                         IF ( .NOT. lspinorb ) THEN
                            DO ipol = 1, npol
                               DO ibnd = 1, m
                                  ps(ikb,ipol,ibnd) = ps(ikb,ipol,ibnd) + &
                                       qq_at(ih,jh,na)*becp%nc(jkb,ipol,ibnd)
                               ENDDO
                            ENDDO
                         ELSE
                            DO ibnd = 1, m
                               ps(ikb,1,ibnd) = ps(ikb,1,ibnd) + &
                                    qq_so(ih,jh,1,nt)*becp%nc(jkb,1,ibnd)+ &
                                    qq_so(ih,jh,2,nt)*becp%nc(jkb,2,ibnd)
                               ps(ikb,2,ibnd) = ps(ikb,2,ibnd) + &
                                    qq_so(ih,jh,3,nt)*becp%nc(jkb,1,ibnd)+ &
                                    qq_so(ih,jh,4,nt)*becp%nc(jkb,2,ibnd)
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
             !
          ENDIF
          !
       ENDDO
       !
       CALL ZGEMM ( 'N', 'N', n, m*npol, nkb, (1.d0,0.d0) , vkb, &
                    lda, ps, nkb, (1.d0,0.d0) , spsi(1,1), lda )
       !
       DEALLOCATE( ps )
       !
       !
       RETURN
       !
    END SUBROUTINE s_psi_nc
    !
END SUBROUTINE pw_s_psi
