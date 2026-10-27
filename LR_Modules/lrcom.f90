!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for LR_Modules routines
!
MODULE control_lr
  !
  USE kinds,      ONLY : DP
  !
  ! ... The variables controlling the run of linear response codes
  !
  SAVE
  !
  INTEGER, ALLOCATABLE :: nbnd_occ(:)  ! occupied bands in metals
  REAL(DP) :: alpha_pv       ! the alpha value for shifting the bands
  LOGICAL  :: lgamma         ! if .TRUE. this is a q=0 computation
  LOGICAL  :: lrpa           ! if .TRUE. uses the Random Phace Approximation
  REAL(DP) :: ethr_nscf      ! convergence threshol for KS eigenvalues in the
                             ! NSCF calculation
  ! Sternheimer case 
  LOGICAL :: lgamma_gamma
  !! if TRUE this is a q=0 computation with k=0 only
  LOGICAL :: ext_recover
  !! if .TRUE. there is a recover file
  LOGICAL :: lnoloc
  !! if .TRUE. neglect local field effects, i.e. ignore the induced potential
  !
  ! Variables for recover
  !
  CHARACTER(LEN=10) :: where_rec = 'no_recover'
  !! where the ph run recovered
  INTEGER :: rec_code = -1000
  !! code for recover
  INTEGER :: rec_code_read = -1000  ! code for recover. Not changed during the run
  !
  INTEGER :: nbnd_occx        ! maximun value of nbnd_occ(:)
  LOGICAL :: reduce_io
  !! if TRUE reduces needed I/O
  !
  ! Parameters controlling DFPT self-consistent iteration
  !
  INTEGER, PARAMETER :: maxter = 150
  !! maximum number of iterations
  LOGICAL :: convt
  !! if TRUE the DFPT has converged
  CHARACTER(LEN=256) :: flmixdpot
  !! File for writing history for potential mixing
  INTEGER :: niter_ph
  !! maximum number of iterations (read from input)
  INTEGER :: nmix_ph
  !! mixing type
  REAL(DP) :: tr2_ph
  !! threshold for DFPT calculation
  REAL(DP) :: alpha_mix(maxter)
  !! the mixing parameter
  !
  LOGICAL :: lnolr = .FALSE.
  !! if TRUE G=0 component of the Hartree term is not added in dv_of_drho
  !
END MODULE control_lr
!
MODULE eqv
  !
  USE kinds,  ONLY : DP
  !
  ! ... The variables describing the linear response problem
  !
  SAVE
  !
  COMPLEX (DP), POINTER :: evq(:,:)
  ! the wavefunctions at point k+q
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:), dpsi(:,:)
  ! the product of dV psi
  ! the change of the wavefunctions
  REAL (DP), ALLOCATABLE :: dmuxc(:,:,:)        ! nrxx, nspin, nspin)
  ! the derivative of the xc potential
  REAL (DP), ALLOCATABLE, TARGET :: vlocq(:,:)  ! ngm, ntyp)
  ! the local potential at q+G
  !
END MODULE eqv
!
MODULE gc_lr
  !
  USE kinds, ONLY : DP
  !
  ! ... The variables needed for gradient corrected calculations
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       grho(:,:,:),    &! gradient of the unperturbed density  (3,nrxx,nspin)
       gmag(:,:,:),    &! 3, nrxx, nspin)
       vsgga(:),       &! nrxx)
       segni(:),       &! nrxx)
       dvxc_rr(:,:,:), &! derivatives of the E_xc functional w.r.t. r and s  
       dvxc_sr(:,:,:), &! r=rho and s=|grad(rho)|
       dvxc_ss(:,:,:), &! dimensions: (nrxx, nspin, nspin)
       dvxc_s(:,:,:)
  !
  ! in the noncollinear case gmag contains the gradient of the magnetization
  ! grho the gradient of rho+ and of rho-, the eigenvalues of the spin density
  ! vsgga= 0.5* (V_up-V_down) to be used in the calculation of the change
  ! of the exchange and correlation magnetic field.
  !
END MODULE gc_lr
!
MODULE lr_symm_base
  !
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the modes and the small group of q
  !
  SAVE
  !
  INTEGER :: irgq(48), nsymq=0, irotmq
  ! selects the operations of the small group
  ! the number of symmetry of the small group
  ! selects the symmetry sending q <-> -q+G
  REAL (DP), ALLOCATABLE :: rtau(:,:,:) !3, 48, nat)
  ! coordinates of direct translations
  REAL (DP) :: gi(3,48), gimq(3)
  ! the possible G associated to each symmetry
  ! the G associated to the symmetry q<->-q+G
  LOGICAL :: minus_q, & ! if .TRUE. there is the symmetry sending q<->-q
             invsymq    ! if .TRUE. the small group of q has inversion
  !
  ! Symmetry representation of the perturbations
  !
  INTEGER :: lr_npert
  !! Number of perturbations considered at the same time.
  !! e.g., for phonons: dimension of the irreducible representation
  !! e.g., for electric fields: 3
  COMPLEX(DP), ALLOCATABLE :: upert(:, :, :)
  !! Representation of the symmetry in the perturbation basis. Size (lr_npert, lr_npert, nsymq)
  !! e.g., for phonons: transformation matrix of the patterns
  !! e.g., for electric fields: transformation matrix of Cartesian vectors
  COMPLEX(DP), ALLOCATABLE :: upert_mq(:, :)
  !! Representation of the symmetry that transforms q to -q. Size (lr_npert, lr_npert)
  !
END MODULE lr_symm_base
!
MODULE lrus
  !
  USE kinds,  ONLY : DP
  USE becmod, ONLY : bec_type
  !
  ! ... These are additional variables needed for the linear response
  ! ... with US pseudopotentials and a generic perturbation Delta Vscf
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE :: &
       int3(:,:,:,:,:),     &! nhm, nhm, nat, nspin, npert)
       int3_paw(:,:,:,:,:), &! nhm, nhm, nat, nspin, npert)
       int3_nc(:,:,:,:,:),  &! nhm, nhm, nat, nspin, npert)
       intq(:,:,:),         &! nhm, nhm, nat)
       intq_nc(:,:,:,:)      ! nhm, nhm, nat, nspin)
  ! int3 -> \int (Delta V_Hxc) Q d^3r
  ! similarly for int_nc while
  ! int3_paw contains Delta (D^1-\tilde D^1)
  ! intq integral of e^iqr Q
  ! intq_nc integral of e^iqr Q in the noncollinear case
  !
  REAL (DP), ALLOCATABLE ::    dpqq(:,:,:,:)       ! nhm, nhm, 3, ntyp)
  COMPLEX (DP), ALLOCATABLE :: dpqq_so(:,:,:,:,:)  ! nhm, nhm, nspin, 3, ntyp)
  ! dpqq and dpqq_so: dipole moment of each Q multiplied by the fcoef factors 
  !
  type (bec_type), ALLOCATABLE, TARGET :: becp1(:) ! nksq)
  ! becp1 contains < beta_n | psi_i >
  !
  REAL (DP),    ALLOCATABLE :: bbg(:,:)      ! nkb, nkb)
  ! for gamma_only     
  COMPLEX (DP), ALLOCATABLE :: bbk(:,:,:)    ! nkb, nkb, nks)
  ! for k points
  COMPLEX (DP), ALLOCATABLE :: bbnc(:,:,:) ! nkb*npol, nkb*npol, nks)
  ! for the noncollinear case
  ! bbg = < beta^N_i | beta^P_j > 
  ! bbg/bbk/bbnc are the scalar products of beta functions 
  ! localized on atoms N and P.
  !
  ! Variables needed for magnetic calculations
  TYPE(bec_type), ALLOCATABLE :: becpt(:)
  TYPE(bec_type), ALLOCATABLE :: alphapt(:,:)
  !
END MODULE lrus
!
MODULE units_lr
  !
  USE kinds,  ONLY : DP
  !
  ! ... These are the units used in the linear response calculations
  !
  SAVE
  !
  INTEGER :: iuwfc,   & ! unit for wavefunctions
             lrwfc,   & ! the length of wavefunction record
             iuatwfc, & ! unit for atomic wavefunctions
             iuatswfc,& ! unit for atomic wavefunctions * S
             iudwf,   & ! unit with D psi
             lrdwf      ! length of D psi record
  !
END MODULE units_lr

MODULE ldaU_lr
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
  !
  REAL(DP) :: effU(ntypx)
  ! effective Hubbard parameter: effU = Hubbard_U - Hubbard_J0
  ! TODO: Can be moved to PW/ldaU
  !
  COMPLEX(DP), ALLOCATABLE :: dnsscf(:,:,:,:,:)
  !! SCF derivative of ns
  !
  COMPLEX(DP), ALLOCATABLE :: vh_u_save(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: vh_uv_save(:,:,:,:,:,:)
  ! to save the two unperturbed Hubbard U and V potentials;
  ! one normal, and the other with the time-reversed m_hubb
  !
  COMPLEX(DP), ALLOCATABLE, TARGET :: swfcatomk(:,:)
  !! S * atomic wfc at k
  COMPLEX(DP), POINTER :: swfcatomkpq(:,:)
  !! S * atomic wfc at k+q
  !
  LOGICAL :: lr_has_dnsorth = .FALSE.
  !! If true, add lr_dnsorth to dnsscf.
  COMPLEX(DP), ALLOCATABLE :: lr_dnsorth(:, :, :, :, :)
  !! Fixed term to be added to dnsscf. Size (ldim, ldim, nspin, nat, 3*nat)
  !
  COMPLEX (DP), ALLOCATABLE :: proj1(:,:), proj2(:,:)
  !! Arrays to store scalar products between vectors
  !
  LOGICAL :: lhubbard_soc = .FALSE.
  !! Input (&INPUTPH): enable noncollinear/SOC fixed-Dudarev DFPT+U
  CHARACTER(LEN=16) :: hubbard_conv = 'potential'
  !! Input (&INPUTPH), with lhubbard_soc: convergence test of the Hubbard
  !! occupation response, required together with the Hxc test.
  !! 'potential': max U*|dN_out-dN_in| <= sqrt(tr2_ph) (Ry, as the Hxc test);
  !! 'relative' : max|dN_out-dN_in|/max|dN| <= sqrt(tr2_ph) (QE 7.3 qe-pert);
  !! 'none'     : Hxc test only (diagnostic).
  REAL(DP) :: hubbard_conv_thr = -1.0_DP
  !! Input (&INPUTPH), with lhubbard_soc: threshold of the Hubbard test
  !! (Ry for 'potential'); <= 0: sqrt(tr2_ph).
  LOGICAL :: hub_bare_time_reversed = .FALSE.
  !! Set by the caller of dvqpsi_us for the second (time-reversed)
  !! magnetic Sternheimer branch; selects the reversed bare Hubbard term
  LOGICAL :: hub_bare_oracle = .FALSE.
  !! Set by the caller of dvqpsi_us to request the iverbosity=1 checks of
  !! the noncollinear bare Hubbard term
  LOGICAL :: hub_dnsq_single_branch = .FALSE.
  !! Noncollinear dnsq_scf: .FALSE. reconstructs dN_scf from both magnetic
  !! Sternheimer records (phonons); .TRUE. uses one record plus its adjoint
  !! transition (static electric field, QE 7.3 implementation)
  !
CONTAINS
  !
  SUBROUTINE pack_hubbard_blocks_nc (ldim_atom, npol_, blocks, matrix)
    !-----------------------------------------------------------------------
    !! Pack PW ns_nc(:,:,npol*(is_bra-1)+is_ket) blocks into the full
    !! orbital-spin matrix layout consumed by vhpsi_U_nc.
    !
    INTEGER, INTENT(IN) :: ldim_atom, npol_
    COMPLEX(DP), INTENT(IN) :: blocks(:,:,:)
    COMPLEX(DP), INTENT(OUT) :: matrix(:,:)
    !
    INTEGER :: is_bra, is_ket, iblock
    INTEGER :: row_first, row_last, col_first, col_last
    !
    IF (SIZE(blocks,1) < ldim_atom .OR. &
        SIZE(blocks,2) < ldim_atom .OR. &
        SIZE(blocks,3) < npol_*npol_ .OR. &
        SIZE(matrix,1) < ldim_atom*npol_ .OR. &
        SIZE(matrix,2) < ldim_atom*npol_) &
       CALL errore('pack_hubbard_blocks_nc', &
                   'invalid Hubbard block shape', 1)
    !
    matrix = (0.0_DP, 0.0_DP)
    DO is_bra = 1, npol_
       row_first = 1 + ldim_atom * (is_bra - 1)
       row_last = ldim_atom * is_bra
       DO is_ket = 1, npol_
          col_first = 1 + ldim_atom * (is_ket - 1)
          col_last = ldim_atom * is_ket
          iblock = npol_ * (is_bra - 1) + is_ket
          matrix(row_first:row_last,col_first:col_last) = &
               blocks(1:ldim_atom,1:ldim_atom,iblock)
       ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE pack_hubbard_blocks_nc
  !
  SUBROUTINE build_dvhub_nc (ldim_atom, npol_, ueff, dns_atom, dvhub)
    !-----------------------------------------------------------------------
    !! Build only the fixed-Dudarev response potential in the same full
    !! orbital-spin layout used by PW vhpsi_U_nc.  The ground-state U/2
    !! identity contribution is intentionally absent.
    !
    INTEGER, INTENT(IN) :: ldim_atom, npol_
    REAL(DP), INTENT(IN) :: ueff
    COMPLEX(DP), INTENT(IN) :: dns_atom(:,:,:)
    COMPLEX(DP), INTENT(OUT) :: dvhub(:,:)
    !
    INTEGER :: is_bra, is_ket, source_block
    INTEGER :: m1, m2, irow, icol
    !
    IF (SIZE(dns_atom,1) < ldim_atom .OR. &
        SIZE(dns_atom,2) < ldim_atom .OR. &
        SIZE(dns_atom,3) < npol_*npol_) &
       CALL errore('build_dvhub_nc', 'invalid occupation-response shape', 1)
    IF (SIZE(dvhub,1) < ldim_atom*npol_ .OR. &
        SIZE(dvhub,2) < ldim_atom*npol_) &
       CALL errore('build_dvhub_nc', 'invalid potential-response shape', 1)
    !
    dvhub = (0.0_DP, 0.0_DP)
    DO is_bra = 1, npol_
       DO is_ket = 1, npol_
          source_block = npol_ * (is_ket - 1) + is_bra
          DO m2 = 1, ldim_atom
             icol = m2 + ldim_atom * (is_ket - 1)
             DO m1 = 1, ldim_atom
                irow = m1 + ldim_atom * (is_bra - 1)
                dvhub(irow,icol) = -ueff * &
                   dns_atom(m2,m1,source_block)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE build_dvhub_nc
  !
  SUBROUTINE time_reverse_hubbard_matrix_nc (ldim_atom, npol_, &
                                              matrix_in, matrix_out)
    !-----------------------------------------------------------------------
    !! Apply epsilon K to an operator represented in the PW orbital-spin
    !! basis: V_TR = epsilon V* epsilon^dagger.
    !
    INTEGER, INTENT(IN) :: ldim_atom, npol_
    COMPLEX(DP), INTENT(IN) :: matrix_in(:,:)
    COMPLEX(DP), INTENT(OUT) :: matrix_out(:,:)
    !
    INTEGER :: is_bra, is_ket, in_bra, in_ket
    INTEGER :: m1, m2, irow, icol, in_row, in_col
    REAL(DP) :: spin_sign
    !
    IF (npol_ /= 2) CALL errore('time_reverse_hubbard_matrix_nc', &
                                'time reversal requires npol=2', 1)
    IF (SIZE(matrix_in,1) < ldim_atom*npol_ .OR. &
        SIZE(matrix_in,2) < ldim_atom*npol_ .OR. &
        SIZE(matrix_out,1) < ldim_atom*npol_ .OR. &
        SIZE(matrix_out,2) < ldim_atom*npol_) &
       CALL errore('time_reverse_hubbard_matrix_nc', &
                   'invalid Hubbard matrix shape', 1)
    !
    matrix_out = (0.0_DP, 0.0_DP)
    DO is_bra = 1, npol_
       DO is_ket = 1, npol_
          in_bra = npol_ + 1 - is_bra
          in_ket = npol_ + 1 - is_ket
          IF (is_bra == is_ket) THEN
             spin_sign = 1.0_DP
          ELSE
             spin_sign = -1.0_DP
          ENDIF
          DO m2 = 1, ldim_atom
             icol = m2 + ldim_atom * (is_ket - 1)
             in_col = m2 + ldim_atom * (in_ket - 1)
             DO m1 = 1, ldim_atom
                irow = m1 + ldim_atom * (is_bra - 1)
                in_row = m1 + ldim_atom * (in_bra - 1)
                matrix_out(irow,icol) = spin_sign * &
                                        CONJG(matrix_in(in_row,in_col))
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE time_reverse_hubbard_matrix_nc
  !
  SUBROUTINE time_reverse_q_response_hubbard_matrix_nc (ldim_atom, npol_, &
                                                         matrix_q, matrix_out)
    !-----------------------------------------------------------------------
    !! Build the Hubbard operator used by the second magnetic branch for a
    !! single nonzero-q response:
    !!
    !!   V_partner(q) = epsilon [V(q)^dagger]* epsilon^dagger
    !!                = epsilon V(q)^T epsilon^dagger.
    !!
    !! The existing time_reverse_hubbard_matrix_nc remains the owner for
    !! Hermitian ground-state and Gamma operators.
    !
    INTEGER, INTENT(IN) :: ldim_atom, npol_
    COMPLEX(DP), INTENT(IN) :: matrix_q(:,:)
    COMPLEX(DP), INTENT(OUT) :: matrix_out(:,:)
    !
    INTEGER :: matrix_size
    COMPLEX(DP), ALLOCATABLE :: matrix_minus_q(:,:)
    !
    matrix_size = ldim_atom * npol_
    IF (SIZE(matrix_q,1) < matrix_size .OR. &
        SIZE(matrix_q,2) < matrix_size .OR. &
        SIZE(matrix_out,1) < matrix_size .OR. &
        SIZE(matrix_out,2) < matrix_size) &
       CALL errore('time_reverse_q_response_hubbard_matrix_nc', &
                   'invalid Hubbard matrix shape', 1)
    !
    ALLOCATE(matrix_minus_q(matrix_size,matrix_size))
    matrix_minus_q = TRANSPOSE(CONJG( &
         matrix_q(1:matrix_size,1:matrix_size)))
    CALL time_reverse_hubbard_matrix_nc(ldim_atom, npol_, &
                                         matrix_minus_q, matrix_out)
    DEALLOCATE(matrix_minus_q)
    !
    RETURN
  END SUBROUTINE time_reverse_q_response_hubbard_matrix_nc
  !
  SUBROUTINE apply_hubbard_matrix_nc (nrow, nbnd_use, nwfc_total, &
                                       ldim_atom, npol_, offset_atom, &
                                       projector, projection, hub_matrix, &
                                       result)
    !-----------------------------------------------------------------------
    !! Accumulate P_I V_I (P_I^dagger psi) with the two BLAS contractions
    !! used by PW vhpsi_U_nc.
    !
    INTEGER, INTENT(IN) :: nrow, nbnd_use, nwfc_total
    INTEGER, INTENT(IN) :: ldim_atom, npol_, offset_atom
    COMPLEX(DP), INTENT(IN) :: projector(:,:), projection(:,:)
    COMPLEX(DP), INTENT(IN) :: hub_matrix(:,:)
    COMPLEX(DP), INTENT(INOUT) :: result(:,:)
    !
    INTEGER :: matrix_size
    COMPLEX(DP), ALLOCATABLE :: coefficients(:,:)
    !
    matrix_size = ldim_atom * npol_
    IF (offset_atom + matrix_size > nwfc_total .OR. &
        SIZE(projector,1) < nrow .OR. &
        SIZE(projector,2) < nwfc_total .OR. &
        SIZE(projection,1) < nwfc_total .OR. &
        SIZE(projection,2) < nbnd_use .OR. &
        SIZE(hub_matrix,1) < matrix_size .OR. &
        SIZE(hub_matrix,2) < matrix_size .OR. &
        SIZE(result,1) < nrow .OR. SIZE(result,2) < nbnd_use) &
       CALL errore('apply_hubbard_matrix_nc', &
                   'invalid projected Hubbard matrix shape', 1)
    !
    ALLOCATE(coefficients(matrix_size,nbnd_use))
    CALL ZGEMM('N', 'N', matrix_size, nbnd_use, matrix_size, &
               (1.0_DP,0.0_DP), hub_matrix, SIZE(hub_matrix,1), &
               projection(offset_atom+1,1), SIZE(projection,1), &
               (0.0_DP,0.0_DP), coefficients, matrix_size)
    CALL ZGEMM('N', 'N', nrow, nbnd_use, matrix_size, &
               (1.0_DP,0.0_DP), projector(1,offset_atom+1), &
               SIZE(projector,1), coefficients, matrix_size, &
               (1.0_DP,0.0_DP), result, SIZE(result,1))
    DEALLOCATE(coefficients)
    !
    RETURN
  END SUBROUTINE apply_hubbard_matrix_nc
  !
  ! accumulate_ns_nc_transition and ns_nc_hermiticity_residual were in
  ! ldaU_ph in the QE 7.3 implementation; their consumers are now in LR_Modules.
  !
  SUBROUTINE accumulate_ns_nc_transition (ldim_atom, offset_atom, npol_, &
                                           nbnd_use, band_weight, proj_left, &
                                           proj_right, dns_atom, time_reversed)
    !-----------------------------------------------------------------------
    !! Accumulate one projected transition matrix in the PW new_ns_nc
    !! convention.  For the second magnetic Sternheimer branch, apply the
    !! charge-even/magnetization-odd map used by incdrhoscf_nc and transpose
    !! the orbital pair, as in the second term of the collinear dnsq_scf
    !! reconstruction.
    !
    INTEGER, INTENT(IN) :: ldim_atom, offset_atom, npol_, nbnd_use
    REAL(DP), INTENT(IN) :: band_weight(:)
    COMPLEX(DP), INTENT(IN) :: proj_left(:,:), proj_right(:,:)
    COMPLEX(DP), INTENT(INOUT) :: dns_atom(:,:,:)
    LOGICAL, INTENT(IN) :: time_reversed
    !
    INTEGER :: is_bra, is_ket, out_bra, out_ket, iblock, out_block
    INTEGER :: m1, m2, out_m1, out_m2, iwfc1, iwfc2
    REAL(DP) :: spin_sign
    COMPLEX(DP) :: transition
    !
    IF (time_reversed .AND. npol_ /= 2) &
       CALL errore('accumulate_ns_nc_transition', &
                   'time-reversal map requires npol=2', 1)
    !
    DO is_bra = 1, npol_
       DO is_ket = 1, npol_
          iblock = npol_ * (is_bra - 1) + is_ket
          IF (time_reversed) THEN
             out_bra = npol_ + 1 - is_ket
             out_ket = npol_ + 1 - is_bra
             out_block = npol_ * (out_bra - 1) + out_ket
             IF (is_bra == is_ket) THEN
                spin_sign = 1.0_DP
             ELSE
                spin_sign = -1.0_DP
             ENDIF
          ELSE
             out_block = iblock
             spin_sign = 1.0_DP
          ENDIF
          !
          DO m2 = 1, ldim_atom
             iwfc2 = offset_atom + m2 + ldim_atom * (is_ket - 1)
             DO m1 = 1, ldim_atom
                iwfc1 = offset_atom + m1 + ldim_atom * (is_bra - 1)
                IF (time_reversed) THEN
                   out_m1 = m2
                   out_m2 = m1
                ELSE
                   out_m1 = m1
                   out_m2 = m2
                ENDIF
                transition = SUM( band_weight(1:nbnd_use) * &
                                  CONJG(proj_left(iwfc1,1:nbnd_use)) * &
                                  proj_right(iwfc2,1:nbnd_use) )
                dns_atom(out_m1,out_m2,out_block) = &
                   dns_atom(out_m1,out_m2,out_block) + spin_sign * transition
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE accumulate_ns_nc_transition
  !
  REAL(DP) FUNCTION ns_nc_hermiticity_residual (ldim_atom, dns_atom)
    !-----------------------------------------------------------------------
    !! Maximum residual of N_ab = CONJG(N_ba) in orbital-spin space.
    !
    INTEGER, INTENT(IN) :: ldim_atom
    COMPLEX(DP), INTENT(IN) :: dns_atom(:,:,:)
    INTEGER :: is_bra, is_ket, iblock, jblock, m1, m2
    !
    ns_nc_hermiticity_residual = 0.0_DP
    DO is_bra = 1, 2
       DO is_ket = 1, 2
          iblock = 2 * (is_bra - 1) + is_ket
          jblock = is_bra + 2 * (is_ket - 1)
          DO m2 = 1, ldim_atom
             DO m1 = 1, ldim_atom
                ns_nc_hermiticity_residual = MAX( &
                   ns_nc_hermiticity_residual, &
                   ABS(dns_atom(m1,m2,iblock) - &
                       CONJG(dns_atom(m2,m1,jblock))) )
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    RETURN
  END FUNCTION ns_nc_hermiticity_residual
  !
END MODULE ldaU_lr
