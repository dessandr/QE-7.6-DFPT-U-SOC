!
! Copyright (C) 2001-2018 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE phqscf
  !-----------------------------------------------------------------------
  !! This subroutine is the main driver of the self consistent cycle
  !! which gives as output the change of the wavefunctions and the
  !! change of the self-consistent potential due to a phonon of
  !! fixed q.
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat
  USE lsda_mod,         ONLY : nspin
  USE io_global,        ONLY : stdout, ionode
  USE fft_base,         ONLY : dffts, dfftp
  USE uspp,             ONLY : okvan
  USE efield_mod,       ONLY : zstarue0, zstarue0_rec
  USE control_ph,       ONLY : zue, current_iq
  USE partial,          ONLY : done_irr, comp_irr
  USE modes,            ONLY : nirr, npert
  USE uspp_param,       ONLY : nhm
  USE paw_variables,    ONLY : okpaw
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE recover_mod,      ONLY : write_rec
  USE mp_pools,         ONLY : inter_pool_comm
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp_images,        ONLY : intra_image_comm
  USE mp,               ONLY : mp_sum, mp_bcast
  USE io_global,        ONLY : ionode_id
  USE ph_restart,       ONLY : ph_readfile, ph_hubbard_checkpoint_begin
  USE lrus,             ONLY : int3, int3_nc, int3_paw
  USE dynmat,           ONLY : dyn_hub_scf
  USE ldaU,             ONLY : lda_plus_u, Hubbard_lmax
  USE ldaU_lr,          ONLY : dnsscf
  USE ldaU_ph,          ONLY : dnsscf_all_modes
  USE units_ph,         ONLY : iundnsscf
  USE control_flags,    ONLY : iverbosity
  USE control_lr,       ONLY : convt, rec_code
  USE dfpt_type,        ONLY : dfpt_data_type, allocate_dfpt_data, deallocate_dfpt_data
  USE write_hub

  IMPLICIT NONE

  INTEGER :: irr, irr1, imode0, npe, ierr
  LOGICAL, ALLOCATABLE :: dnsscf_present(:)
  ! counter on the representations
  ! counter on the representations
  ! counter on the modes
  ! npert(irr)
  REAL(DP) :: tcpu, get_clock
  ! timing variables
  TYPE(dfpt_data_type) :: dfpt_data
  !! Data that describes linear response quantities
  !
  EXTERNAL get_clock
  ! the change of density due to perturbations

  CALL start_clock ('phqscf')
  !
  ! DFPT+U
  !
  IF (lda_plus_u) THEN
     ALLOCATE (dnsscf_all_modes(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,3*nat))
     ALLOCATE (dyn_hub_scf(3*nat,3*nat))
     dnsscf_all_modes = (0.d0, 0.d0)
     dyn_hub_scf      = (0.d0, 0.d0)
     IF (noncolin) THEN
        ! Completion of dyn_rec does not establish availability of dnsscf.
        ! Preserve an old flat export, and keep the canonical file explicitly
        ! unreadable as response data until every irrep is present.
        CALL prepare_dnsscf_export()
        CALL ph_hubbard_checkpoint_begin()
        ALLOCATE(dnsscf_present(nirr))
        dnsscf_present = .FALSE.
        DO irr = 1, nirr
           IF (.NOT.done_irr(irr)) CYCLE
           CALL ph_readfile('data_hub', current_iq, irr, ierr)
           dnsscf_present(irr) = ierr == 0
           IF (dnsscf_present(irr)) THEN
              WRITE(stdout,'(5x,"SOC_U_DNSSCF_RESTORED",2(1x,i0))') current_iq, irr
           ELSE
              WRITE(stdout,'(5x,"SOC_U_DNSSCF_RESTORE_MISSING",2(1x,i0))') current_iq, irr
           ENDIF
        ENDDO
     ENDIF
  ENDIF
  !
  ! For each irreducible representation we compute the change
  ! of the wavefunctions
  !
  DO irr = 1, nirr
     IF ( (comp_irr (irr)) .AND. (.NOT.done_irr (irr)) ) THEN
        npe=npert(irr)
        !
        CALL allocate_dfpt_data(dfpt_data, npe)
        imode0 = 0
        DO irr1 = 1, irr - 1
           imode0 = imode0 + npert (irr1)
        ENDDO
        IF (npe == 1) THEN
           WRITE( stdout, '(//,5x,"Representation #",i4," mode #",i4)') &
                              irr, imode0 + 1
        ELSE
           WRITE( stdout, '(//,5x,"Representation #",i4," modes #",8i4)') &
                              irr, (imode0+irr1, irr1=1,npe)
        ENDIF
        !
        ! Initialize LR_Modules variables
        !
        CALL ph_set_upert_phonon(irr)
        IF (okvan) CALL dnsq_orth_set_irr(irr, imode0)
        !
        !    then for this irreducible representation we solve the linear system
        !
        IF (okvan) THEN
           ALLOCATE (int3 ( nhm, nhm, nat, nspin_mag, npe))
           IF (okpaw) ALLOCATE (int3_paw (nhm, nhm, nat, nspin_mag, npe))
           IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, nat, nspin, npe))
        ENDIF
        !
        ! DFPT+U: dnsscf in the phonon calculation
        ! is the scf change of atomic occupations ns 
        ! induced by atomic displacements.
        !
        IF (lda_plus_u) THEN
           ALLOCATE (dnsscf(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,npe))
           dnsscf = (0.d0, 0.d0)
        ENDIF
        !
        WRITE( stdout, '(/,5x,"Self-consistent Calculation")')
        CALL solve_linter (irr, imode0, dfpt_data)
        WRITE( stdout, '(/,5x,"End of self-consistent calculation")')
        !
        !   Add the contribution of this mode to the dynamical matrix
        !
        IF (convt) THEN
           CALL drhodv (imode0, npe, dfpt_data%drhos)
           !
           !   add the contribution of the modes imode0+1 -> imode+npe
           !   to the effective charges Z(Us,E) (Us=scf,E=bare)
           !
           IF (zue) CALL add_zstar_ue (imode0, npe )
           IF (zue.AND. okvan) CALL add_zstar_ue_us(imode0, npe )
           !
           ! DFPT+U: calculate the scf part of the Hubbard dynamical matrix
           IF (lda_plus_u) THEN
              CALL dynmat_hub_scf (irr, imode0, npe)
              IF (noncolin) dnsscf_present(irr) = .TRUE.
              DEALLOCATE (dnsscf)
           ENDIF
           !
           IF (zue) THEN
              CALL mp_sum ( zstarue0_rec, intra_bgrp_comm )
              CALL mp_sum ( zstarue0_rec, inter_pool_comm )
              zstarue0(:,:)=zstarue0(:,:)+zstarue0_rec(:,:)
           ENDIF
           !
           WRITE( stdout, '(/,5x,"Convergence has been achieved ")')
           done_irr (irr) = .TRUE.
        ELSE
           WRITE( stdout, '(/,5x,"No convergence has been achieved ")')
           CALL stop_smoothly_ph (.FALSE.)
        ENDIF
        rec_code=20
        CALL write_rec('done_drhod',irr,0.0_DP,-1000,.false., dfpt_data)
        !
        IF (okvan) THEN
           DEALLOCATE (int3)
           IF (okpaw) DEALLOCATE (int3_paw)
           IF (noncolin) DEALLOCATE(int3_nc)
        ENDIF
        CALL deallocate_dnsorth()
        CALL ph_deallocate_upert()
        CALL deallocate_dfpt_data(dfpt_data)
        !
        tcpu = get_clock ('PHONON')
        !
     ENDIF

  ENDDO
  !
  ! DFPT+U
  !
  IF (lda_plus_u) THEN
     !
     ! Write dnsscf_all_modes in the pattern basis u to file,
     ! because it is needed for el-ph calculations
     !
     IF (noncolin) THEN
        IF (ALL(dnsscf_present)) THEN
           IF (ionode) THEN
              REWIND(iundnsscf)
              WRITE(iundnsscf,*) dnsscf_all_modes
              ENDFILE(iundnsscf)
              FLUSH(iundnsscf)
           ENDIF
           WRITE(stdout,'(5x,"SOC_U_DNSSCF_EXPORT_COMPLETE",2(1x,i0))') current_iq, nirr
        ELSE
           ! Partial/image calculations remain valid PH calculations. Their
           ! per-irrep XML checkpoints can be collected, but this image must
           ! not publish a flat qe2pert input with zero-filled missing modes.
           WRITE(stdout,'(5x,"SOC_U_DNSSCF_EXPORT_SKIPPED",3(1x,i0))') &
                current_iq, COUNT(dnsscf_present), nirr
           WRITE(stdout,'(5x,a)') 'Incomplete Hubbard response: collect compatible per-irrep checkpoints '// &
                'and recover all modes on one image before using qe2pert.'
        ENDIF
     ELSE
        IF (ionode) WRITE(iundnsscf,*) dnsscf_all_modes
     ENDIF
     !
     ! Write dnsscf_all_modes in the cartesian coordinates
     ! to the standard output
     !
     IF (iverbosity==1) THEN
        IF (noncolin) THEN
           IF (ALL(dnsscf_present)) CALL write_dnsscf_ph()
        ELSE
           CALL write_dnsscf_ph()
        ENDIF
     ENDIF
     !
     ! Write the SCF and total Hubbard dynamical matrix
     ! to the standard output
     !
     IF (iverbosity==1) CALL write_dynmat_hub()
     !
     DEALLOCATE (dnsscf_all_modes)
     DEALLOCATE (dyn_hub_scf)
     IF (ALLOCATED(dnsscf_present)) DEALLOCATE(dnsscf_present)
     !
  ENDIF
  !
  CALL stop_clock ('phqscf')
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE prepare_dnsscf_export()
    USE clib_wrappers, ONLY : f_rename
    IMPLICIT NONE
    CHARACTER(len=1024) :: filename, backup
    CHARACTER(len=24) :: suffix
    INTEGER :: status, file_size, backup_index
    LOGICAL :: opened, exists
    status = 0
    IF (ionode) THEN
       opened = .FALSE.
       file_size = -1
       INQUIRE(UNIT=iundnsscf, OPENED=opened, NAME=filename, SIZE=file_size, IOSTAT=status)
       IF (status == 0 .AND. .NOT.opened) status = 1
       IF (status == 0) THEN
          IF (file_size > 0) THEN
             backup_index = 0
             DO
                backup_index = backup_index + 1
                WRITE(suffix,'(i0)') backup_index
                backup = TRIM(filename)//'.previous.'//TRIM(suffix)
                INQUIRE(FILE=TRIM(backup), EXIST=exists)
                IF (.NOT.exists) EXIT
             ENDDO
             CLOSE(iundnsscf, STATUS='keep', IOSTAT=status)
             IF (status == 0) status = f_rename(TRIM(filename), TRIM(backup))
             IF (status == 0) THEN
                WRITE(stdout,'(5x,a,1x,a)') 'SOC_U_DNSSCF_PREVIOUS_EXPORT', TRIM(backup)
                OPEN(UNIT=iundnsscf, FILE=TRIM(filename), FORM='formatted', STATUS='new', IOSTAT=status)
             ENDIF
          ENDIF
          IF (status == 0) THEN
             REWIND(iundnsscf)
             WRITE(iundnsscf,'(a,1x,i0)',IOSTAT=status) 'SOC_U_DNSSCF_INCOMPLETE q', current_iq
             IF (status == 0) THEN
                ENDFILE(iundnsscf, IOSTAT=status)
                FLUSH(iundnsscf)
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    CALL mp_bcast(status, ionode_id, intra_image_comm)
    IF (status /= 0) CALL errore('phqscf', 'cannot preserve or invalidate the old Hubbard response export', 1)
  END SUBROUTINE prepare_dnsscf_export
  !
END SUBROUTINE phqscf
