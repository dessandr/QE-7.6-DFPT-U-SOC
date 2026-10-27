!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE hubbard_response_mixing
  ! Keep the joint adapter with the existing Broyden/packing owner. Both
  ! phonon and electric drivers supply their current, unmixed response.
  USE kinds, ONLY : DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: mix_hubbard_response
CONTAINS
  SUBROUTINE mix_hubbard_response(hxc_out, hxc_in, hub_out, hub_in, ueff, &
       ldim, alphamix, dr2, tr2, rtol, iter, n_iter, file_extension, conv, &
       hub_residual, hub_scale, hub_relative, hub_convt, hxc_grid_size, metric_info, &
       require_hubbard_convergence, hubbard_criterion, hub_pot_residual, &
       hubbard_mix, alpha_hub)
    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp, ONLY : mp_sum
    IMPLICIT NONE
    COMPLEX(DP), INTENT(IN) :: hxc_out(:,:,:), hub_out(:,:,:,:,:)
    COMPLEX(DP), INTENT(INOUT) :: hxc_in(:,:,:), hub_in(:,:,:,:,:)
    REAL(DP), INTENT(IN) :: ueff(:), alphamix, tr2, rtol
    INTEGER, INTENT(IN) :: ldim(:), n_iter
    INTEGER, INTENT(INOUT) :: iter
    CHARACTER(LEN=256), INTENT(IN) :: file_extension
    REAL(DP), INTENT(OUT) :: dr2, hub_residual, hub_scale, hub_relative
    LOGICAL, INTENT(OUT) :: conv, hub_convt
    REAL(DP), OPTIONAL, INTENT(IN) :: hxc_grid_size
    ! Global physical complex field count, including spin/perturbations but
    ! excluding FFT padding. The fallback is intended for unpadded fixtures.
    REAL(DP), OPTIONAL, INTENT(OUT) :: metric_info(3)
    ! Raw Hxc RMS, raw induced-Hubbard-potential RMS, fixed onsite scale.
    LOGICAL, OPTIONAL, INTENT(IN) :: require_hubbard_convergence
    ! Validation policy only.  The default preserves production's joint gate;
    ! false keeps the joint coordinates but lets the Hxc prefix decide exit.
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: hubbard_criterion
    ! 'relative' (default): max|dN_out-dN_in| / max|dN| <= rtol.
    ! 'potential': max_I U_I*max|dN_out-dN_in| <= rtol, i.e. the induced
    ! Hubbard potential residual in Ry, as the Hxc residual; well defined
    ! when the response vanishes by symmetry.
    REAL(DP), OPTIONAL, INTENT(OUT) :: hub_pot_residual
    ! max_I U_I*max|dN_out-dN_in| (Ry)
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: hubbard_mix
    ! 'joint' (default): U*dN in the Hxc Broyden vector, with initial step
    ! alpha_hub; 'linear': Broyden on the Hxc potential alone and
    ! dN_in <- dN_in + alpha_hub*(dN_out-dN_in).
    REAL(DP), OPTIONAL, INTENT(IN) :: alpha_hub
    ! <= 0 or absent: alphamix ('joint') or 1 ('linear')
    INTEGER :: nhxc, nhub, nlocal, na, ld, nreal, nactive
    REAL(DP) :: roundoff_floor, field_count, allocated_count, hub_weight
    REAL(DP) :: hxc_r2, hub_pot_r2, hub_pot_max, alpha_onsite
    LOGICAL :: auxiliary_convt, potential_criterion, linear_mixing
    COMPLEX(DP), ALLOCATABLE :: pot_in(:,:,:,:,:), pot_out(:,:,:,:,:), hub_prev(:,:,:,:,:)
    COMPLEX(DP), ALLOCATABLE :: flat_in(:), flat_out(:), mixin(:), mixout(:)
    COMPLEX(DP), ALLOCATABLE :: field_in(:), field_out(:)
    REAL(DP), ALLOCATABLE :: real_in(:), real_out(:)
    !
    CALL start_clock('mix_hubbard')
    IF (ANY(SHAPE(hxc_in) /= SHAPE(hxc_out)) .OR. &
        ANY(SHAPE(hub_in) /= SHAPE(hub_out))) &
       CALL errore('mix_hubbard_response', 'inconsistent response shapes', 1)
    IF (SIZE(hub_out,3) /= 4 .OR. SIZE(ueff) /= SIZE(hub_out,4) .OR. &
        SIZE(ldim) /= SIZE(hub_out,4)) &
       CALL errore('mix_hubbard_response', 'expected four Hubbard spin blocks', 1)
    IF (.NOT.ieee_is_finite(rtol) .OR. rtol <= 0.0_DP .OR. &
        .NOT.ALL(ieee_is_finite(ueff))) &
       CALL errore('mix_hubbard_response', 'invalid tolerance or Hubbard U', 1)
    IF (.NOT.ALL(ieee_is_finite(REAL(hub_out,DP))) .OR. &
        .NOT.ALL(ieee_is_finite(AIMAG(hub_out))) .OR. &
        .NOT.ALL(ieee_is_finite(REAL(hub_in,DP))) .OR. &
        .NOT.ALL(ieee_is_finite(AIMAG(hub_in)))) &
       CALL errore('mix_hubbard_response', 'nonfinite Hubbard response', 1)
    !
    nhxc = SIZE(hxc_in)
    nhub = SIZE(hub_in)
    allocated_count = REAL(nhxc,DP)
    CALL mp_sum(allocated_count,intra_bgrp_comm)
    field_count = allocated_count
    IF (PRESENT(hxc_grid_size)) field_count = hxc_grid_size
    IF (.NOT.ieee_is_finite(field_count) .OR. field_count <= 0.0_DP .OR. &
        field_count > allocated_count) &
       CALL errore('mix_hubbard_response', 'invalid global Hxc grid size', 1)
    hxc_r2 = SUM(ABS(hxc_out-hxc_in)**2)
    CALL mp_sum(hxc_r2,intra_bgrp_comm)
    ALLOCATE(pot_in, MOLD=hub_in)
    ALLOCATE(pot_out, MOLD=hub_out)
    pot_in = (0.0_DP,0.0_DP)
    pot_out = (0.0_DP,0.0_DP)
    hub_residual = 0.0_DP
    hub_scale = 0.0_DP
    hub_pot_r2 = 0.0_DP
    hub_pot_max = 0.0_DP
    potential_criterion = .FALSE.
    IF (PRESENT(hubbard_criterion)) THEN
       SELECT CASE (TRIM(hubbard_criterion))
       CASE ('relative')
       CASE ('potential')
          potential_criterion = .TRUE.
       CASE DEFAULT
          CALL errore('mix_hubbard_response', 'unknown Hubbard convergence criterion', 1)
       END SELECT
    ENDIF
    linear_mixing = .FALSE.
    IF (PRESENT(hubbard_mix)) THEN
       SELECT CASE (TRIM(hubbard_mix))
       CASE ('joint')
       CASE ('linear')
          linear_mixing = .TRUE.
       CASE DEFAULT
          CALL errore('mix_hubbard_response', 'unknown Hubbard mixing mode', 1)
       END SELECT
    ENDIF
    alpha_onsite = MERGE(1.0_DP, alphamix, linear_mixing)
    IF (PRESENT(alpha_hub)) THEN
       IF (alpha_hub > 0.0_DP) alpha_onsite = alpha_hub
    ENDIF
    IF (.NOT.ieee_is_finite(alpha_onsite)) &
       CALL errore('mix_hubbard_response', 'nonfinite Hubbard mixing step', 1)
    nactive = 0
    DO na = 1, SIZE(ueff)
       ld = ldim(na)
       IF (ld < 0 .OR. ld > MIN(SIZE(hub_out,1),SIZE(hub_out,2))) &
          CALL errore('mix_hubbard_response', 'invalid onsite manifold size', 1)
       IF (ld == 0 .OR. ueff(na) == 0.0_DP) CYCLE
       ! Metadata/onsite arrays are replicated: do NOT MPI-sum this count or
       ! the onsite norm. setmixout distributes their contribution later.
       nactive = nactive + 4*ld*ld*SIZE(hub_in,5)
       ! -U times the stored occupation blocks is an isometric index
       ! permutation of build_dvhub_nc's potential. No new operator or
       ! Hermiticity projection is introduced by the mixing metric.
       pot_in(1:ld,1:ld,:,na,:) = -ueff(na)*hub_in(1:ld,1:ld,:,na,:)
       pot_out(1:ld,1:ld,:,na,:) = -ueff(na)*hub_out(1:ld,1:ld,:,na,:)
       hub_pot_r2 = hub_pot_r2 + &
            SUM(ABS(pot_out(1:ld,1:ld,:,na,:)-pot_in(1:ld,1:ld,:,na,:))**2)
       hub_residual = MAX(hub_residual, &
            MAXVAL(ABS(hub_out(1:ld,1:ld,:,na,:)-hub_in(1:ld,1:ld,:,na,:))))
       hub_pot_max = MAX(hub_pot_max, ABS(ueff(na)) * &
            MAXVAL(ABS(hub_out(1:ld,1:ld,:,na,:)-hub_in(1:ld,1:ld,:,na,:))))
       hub_scale = MAX(hub_scale, MAXVAL(ABS(hub_out(1:ld,1:ld,:,na,:))), &
                                 MAXVAL(ABS(hub_in(1:ld,1:ld,:,na,:))))
    ENDDO
    ! Retain the relative occupation criterion for resolved responses; a
    ! roundoff-only absolute floor avoids dividing a vanishing response by
    ! TINY. Zero-U observables remain raw outputs, not feedback variables.
    roundoff_floor = 32.0_DP*EPSILON(1.0_DP)*MAX(1.0_DP,hub_scale)
    hub_relative = hub_residual / MAX(hub_scale,roundoff_floor/rtol)
    hub_convt = hub_relative <= rtol
    IF (potential_criterion) hub_convt = hub_pot_max <= rtol
    IF (PRESENT(hub_pot_residual)) hub_pot_residual = hub_pot_max
    IF (.NOT.ieee_is_finite(hxc_r2) .OR. .NOT.ieee_is_finite(hub_pot_r2)) &
       CALL errore('mix_hubbard_response', 'nonfinite block residual norm', 1)
    IF (PRESENT(metric_info)) metric_info = [SQRT(hxc_r2/field_count),0.0_DP,0.0_DP]
    ! With no Hubbard feedback use the unmodified Hxc coordinates and core.
    ! This also avoids dividing by zero U or adding inactive history entries.
    IF (nactive == 0) THEN
       nreal = 2*nhxc
       ALLOCATE(real_in(nreal),real_out(nreal))
       real_in = TRANSFER(hxc_in,real_in)
       real_out = TRANSFER(hxc_out,real_out)
       CALL mix_potential_core(nreal,real_out,real_in,alphamix,dr2,tr2,iter, &
            n_iter,file_extension,conv,nreal,.TRUE.,alphamix)
       hxc_in = RESHAPE(TRANSFER(real_in,hxc_in),SHAPE(hxc_in))
       hub_in = hub_out
       DEALLOCATE(real_in,real_out,pot_in,pot_out)
       CALL stop_clock('mix_hubbard')
       RETURN
    ENDIF
    auxiliary_convt = hub_convt
    IF (PRESENT(require_hubbard_convergence)) THEN
       IF (.NOT.require_hubbard_convergence) auxiliary_convt = .TRUE.
    ENDIF
    IF (linear_mixing) THEN
       ! Broyden on the Hxc potential alone; the Hubbard test still gates the
       ! exit. alpha_onsite = 1 uses the unmixed dN_out, as collinear DFPT+U.
       IF (PRESENT(metric_info)) metric_info(2:3) = &
            [SQRT(hub_pot_r2/REAL(nactive,DP)),0.0_DP]
       nreal = 2*nhxc
       ALLOCATE(real_in(nreal),real_out(nreal))
       real_in = TRANSFER(hxc_in,real_in)
       real_out = TRANSFER(hxc_out,real_out)
       CALL mix_potential_core(nreal,real_out,real_in,alphamix,dr2,tr2,iter, &
            n_iter,file_extension,conv,nreal,auxiliary_convt,alphamix)
       hxc_in = RESHAPE(TRANSFER(real_in,hxc_in),SHAPE(hxc_in))
       ALLOCATE(hub_prev, SOURCE=hub_in)
       hub_in = hub_out
       DO na = 1, SIZE(ueff)
          ld = ldim(na)
          IF (ld == 0 .OR. ueff(na) == 0.0_DP) CYCLE
          hub_in(1:ld,1:ld,:,na,:) = hub_prev(1:ld,1:ld,:,na,:) + alpha_onsite* &
               (hub_out(1:ld,1:ld,:,na,:)-hub_prev(1:ld,1:ld,:,na,:))
       ENDDO
       DEALLOCATE(real_in,real_out,pot_in,pot_out,hub_prev)
       CALL stop_clock('mix_hubbard')
       RETURN
    ENDIF
    ! Block-RMS Broyden metric, up to a harmless common factor field_count:
    ! ||r_Hxc||^2/N_Hxc + ||-U*r_N||^2/N_U.
    ! Only the onsite coordinates are rescaled, leaving the core's physical
    ! Hxc dr2, stopping criterion, and QE adaptive-CG input exactly unchanged.
    ! Counts are fixed throughout an irrep: every stored secant uses the same
    ! metric. No residual-amplitude-dependent reweighting is allowed here.
    hub_weight = SQRT(field_count/REAL(nactive,DP))
    pot_in = hub_weight*pot_in
    pot_out = hub_weight*pot_out
    IF (PRESENT(metric_info)) metric_info(2:3) = &
         [SQRT(hub_pot_r2/REAL(nactive,DP)),hub_weight]
    !
    ALLOCATE(flat_in(nhub),flat_out(nhub),field_in(nhxc),field_out(nhxc))
    ALLOCATE(mixin(nhxc+nhub),mixout(nhxc+nhub))
    flat_in = RESHAPE(pot_in,[nhub])
    flat_out = RESHAPE(pot_out,[nhub])
    field_in = RESHAPE(hxc_in,[nhxc])
    field_out = RESHAPE(hxc_out,[nhxc])
    ! Reuse the PAW auxiliary-block distribution: each Hubbard element
    ! contributes on exactly one band-group rank to every Broyden dot product.
    CALL setmixout(nhxc,nhub,mixin,field_in,flat_in,nlocal,-1)
    CALL setmixout(nhxc,nhub,mixout,field_out,flat_out,nlocal,-1)
    nreal = 2*(nhxc+nlocal)
    ALLOCATE(real_in(nreal),real_out(nreal))
    real_in = TRANSFER(mixin(1:nhxc+nlocal),real_in)
    real_out = TRANSFER(mixout(1:nhxc+nlocal),real_out)
    CALL mix_potential_core(nreal,real_out,real_in,alphamix,dr2,tr2,iter, &
         n_iter,file_extension,conv,2*nhxc,auxiliary_convt,alpha_onsite)
    mixin(1:nhxc+nlocal) = TRANSFER(real_in,mixin(1:nhxc+nlocal))
    CALL setmixout(nhxc,nhub,mixin,field_in,flat_in,nlocal,1)
    hxc_in = RESHAPE(field_in,SHAPE(hxc_in))
    pot_in = RESHAPE(flat_in,SHAPE(hub_in))
    hub_in = hub_out
    DO na = 1, SIZE(ueff)
       ld = ldim(na)
       IF (ld == 0 .OR. ueff(na) == 0.0_DP) CYCLE
       hub_in(1:ld,1:ld,:,na,:) = -pot_in(1:ld,1:ld,:,na,:)/(hub_weight*ueff(na))
    ENDDO
    ! hub_out is deliberately untouched: dynamical assembly and response
    ! checkpoints must consume the current wavefunction-derived response.
    DEALLOCATE(pot_in,pot_out,flat_in,flat_out,field_in,field_out,mixin,mixout,real_in,real_out)
    CALL stop_clock('mix_hubbard')
  END SUBROUTINE mix_hubbard_response
END MODULE hubbard_response_mixing

! Legacy ABI and convergence norm are unchanged for all existing callers.
SUBROUTINE mix_potential(ndim,vout,vin,alphamix,dr2,tr2,iter,n_iter,file_extension,conv)
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: ndim, iter, n_iter
  REAL(DP) :: vout(ndim), vin(ndim), alphamix, dr2, tr2
  CHARACTER(LEN=256) :: file_extension
  LOGICAL :: conv
  CALL mix_potential_core(ndim,vout,vin,alphamix,dr2,tr2,iter,n_iter, &
       file_extension,conv,ndim,.TRUE.,alphamix)
END SUBROUTINE mix_potential

!-----------------------------------------------------------------------
subroutine mix_potential_core (ndim, vout, vin, alphamix, dr2, tr2, &
     iter, n_iter, file_extension, conv, nconv, auxiliary_convt, alpha_aux)
  !-----------------------------------------------------------------------
  !
  ! Modified Broyden's method for potential/charge density mixing
  !             D.D.Johnson, PRB 38, 12807 (1988)
  ! On input :
  !    ndim      dimension of arrays vout, vin
  !    vout      output potential/rho at current iteration
  !    vin       potential/rho at previous iteration
  !    alphamix  mixing factor (0 < alphamix <= 1)
  !    alpha_aux mixing factor of the components after nconv (block-diagonal
  !              initial inverse Jacobian; alpha_aux = alphamix: standard)
  !    tr2       threshold for selfconsistency
  !    iter      current iteration number
  !    n_iter    number of iterations used in the mixing
  !    file_extension  if present save previous iterations on
  !                    file 'prefix'.'file_extension'
  !                    otherwise keep everything in memory
  !
  ! On output:
  !    dr2       [(vout-vin)/ndim]^2
  !    vin       mixed potential
  !    vout      vout-vin
  !    conv      true if dr2.le.tr2
  !
  USE kinds,           ONLY : DP
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum
  USE io_files,        ONLY : diropn
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  !
  implicit none
  !
  !   First the dummy variables
  !
  character (len=256) :: file_extension
  integer :: ndim, iter, n_iter, nconv
  real(DP) :: vout (ndim), vin (ndim), alphamix, dr2, tr2, alpha_aux
  logical :: conv, auxiliary_convt
  !
  !   Here the local variables
  !
  logical :: saveonfile, exst
  integer :: iunmix, n, i, j, info, iter_used, ipos, inext, ndimtot
  real(DP) :: gamma, norm
  real(DP), allocatable, save :: df (:,:), dv (:,:)
  !! work space containing info from previous iterations:
  !! must be kept in memory and saved between calls if file_extension=' '
  !
  integer, allocatable :: iwork(:)
  real(DP), allocatable :: vinsave(:), beta(:, :), work(:), w(:)
  !
  real(DP) :: w0 = 0.01_dp
  real(DP) :: w1 = 1.0_dp
  !! adjustable parameters as suggested in the original paper. w(:) is set to w1.
  !
  INTEGER, EXTERNAL :: find_free_unit
  REAL(DP), EXTERNAL :: ddot, dnrm2
  !
  call start_clock ('mix_pot')
  !
  if (iter < 1) call errore ('mix_potential', 'iter must be positive', 1)
  if (ndim < 1) call errore ('mix_potential', 'ndim must be positive', 3)
  if (nconv < 1 .OR. nconv > ndim) &
       call errore('mix_potential', 'invalid convergence-prefix length', 1)
  if (n_iter < 1) call errore('mix_potential', 'empty mixing history', 1)
  !
  saveonfile = file_extension /= ' '
  !
  do n = 1, ndim
     vout (n) = vout (n) - vin (n)
  enddo
  !
  ! The Broyden metric uses the full vector, but retain the original Hxc
  ! convergence norm/normalization when an auxiliary response is appended.
  dr2 = dnrm2 (nconv, vout, 1) **2
  ndimtot = nconv
  !
  call mp_sum (dr2, intra_bgrp_comm)
  call mp_sum (ndimtot, intra_bgrp_comm)
  !
  dr2 = (sqrt (dr2) / ndimtot) **2
  if (.NOT.ieee_is_finite(dr2)) &
       call errore('mix_potential', 'nonfinite potential residual', 1)
  !
  conv = dr2 < tr2 .AND. auxiliary_convt
  !
  if (saveonfile) then
     iunmix = find_free_unit()
     if (conv) then
        ! remove temporary file (open and close it)
        call diropn (iunmix, file_extension, ndim, exst)
        close (unit=iunmix, status='delete')
        call stop_clock ('mix_pot')
        return
     endif
     call diropn (iunmix, file_extension, ndim, exst)
     if (iter > 1 .AND. .NOT. exst) then
        call infomsg ('mix_potential', 'file not found, restarting')
        iter = 1
     endif
     allocate (df( ndim , n_iter))
     allocate (dv( ndim , n_iter))
  else
     if (iter == 1) then
        ! A recovered joint input cold-starts history; it must never reuse
        ! a previous irrep's or differently packed input's secants.
        if (allocated(df)) deallocate(df)
        if (allocated(dv)) deallocate(dv)
        allocate (df( ndim , n_iter))
        allocate (dv( ndim , n_iter))
     endif
     if (conv) then
        deallocate (dv)
        deallocate (df)
        call stop_clock ('mix_pot')
        return
     endif
     allocate (vinsave( ndim))
  endif
  !
  ! iter_used = iter-1  if iter <= n_iter
  ! iter_used = n_iter  if iter >  n_iter
  !
  iter_used = min (iter - 1, n_iter)
  !
  ! ipos is the position in which results from the present iteraction
  ! are stored. ipos=iter-1 until ipos=n_iter, then back to 1,2,...
  !
  ipos = iter - 1 - ( (iter - 2) / n_iter) * n_iter
  !
  if (iter > 1) then
     if (saveonfile) then
        call davcio (df (1, ipos), ndim, iunmix, 1, - 1)
        call davcio (dv (1, ipos), ndim, iunmix, 2, - 1)
     endif
     do n = 1, ndim
        df (n, ipos) = vout (n) - df (n, ipos)
        dv (n, ipos) = vin (n) - dv (n, ipos)
     enddo
     norm = (dnrm2 (ndim, df (1, ipos), 1) ) **2
     call mp_sum (norm, intra_bgrp_comm)
     if (.NOT.ieee_is_finite(norm)) &
          call errore('mix_potential', 'nonfinite Broyden secant', 1)
     if (norm <= TINY(1.0_DP)) then
        ! A repeated residual supplies no usable secant. Start a valid
        ! linear step/history cycle instead of normalizing a zero vector.
        iter = 1
        iter_used = 0
        ipos = 0
     else
        norm = sqrt (norm)
        call dscal (ndim, 1.d0 / norm, df (1, ipos), 1)
        call dscal (ndim, 1.d0 / norm, dv (1, ipos), 1)
     endif
  endif
  !
  if (saveonfile) then
     do i = 1, iter_used
        if (i /= ipos) then
           call davcio (df (1, i), ndim, iunmix, 2 * i + 1, - 1)
           call davcio (dv (1, i), ndim, iunmix, 2 * i + 2, - 1)
        endif
     enddo
     call davcio (vout, ndim, iunmix, 1, 1)
     call davcio (vin, ndim, iunmix, 2, 1)
     if (iter > 1) then
        call davcio (df (1, ipos), ndim, iunmix, 2 * ipos + 1, 1)
        call davcio (dv (1, ipos), ndim, iunmix, 2 * ipos + 2, 1)
     endif
  else
     call DCOPY (ndim, vin, 1, vinsave, 1)
  endif
  !
  if ( iter_used > 0 ) then
     !
     ALLOCATE(beta(iter_used, iter_used))
     ALLOCATE(w(iter_used))
     ALLOCATE(work(iter_used))
     ALLOCATE(iwork(iter_used))
     w(:) = w1
     !
     beta = 0.0_dp
     do i = 1, iter_used
        do j = i + 1, iter_used
           beta (i, j) = w (i) * w (j) * ddot (ndim, df (1, j), 1, df (1, i), 1)
           call mp_sum ( beta (i, j), intra_bgrp_comm )
        enddo
        beta (i, i) = w0**2 + w (i) **2
     enddo
     !
     call DSYTRF ('U', iter_used, beta, iter_used, iwork, work, iter_used, info)
     call errore ('broyden', 'factorization', info)
     call DSYTRI ('U', iter_used, beta, iter_used, iwork, work, info)
     call errore ('broyden', 'DSYTRI', info)
     deallocate ( iwork )
     !
     do i = 1, iter_used
        do j = i + 1, iter_used
           beta (j, i) = beta (i, j)
        enddo
     enddo
     !
     do i = 1, iter_used
        work (i) = ddot (ndim, df (1, i), 1, vout, 1)
     enddo
     call mp_sum ( work, intra_bgrp_comm )
     !
  end if
  !
  do n = 1, nconv
     vin (n) = vin (n) + alphamix * vout (n)
  enddo
  do n = nconv + 1, ndim
     vin (n) = vin (n) + alpha_aux * vout (n)
  enddo
  !
  do i = 1, iter_used
     gamma = 0.d0
     do j = 1, iter_used
        gamma = gamma + beta (j, i) * w (j) * work (j)
     enddo
     !
     do n = 1, nconv
        vin (n) = vin (n) - w (i) * gamma * (alphamix * df (n, i) + dv (n, i) )
     enddo
     do n = nconv + 1, ndim
        vin (n) = vin (n) - w (i) * gamma * (alpha_aux * df (n, i) + dv (n, i) )
     enddo
  enddo
  !
  IF ( iter_used > 0 ) THEN
     deallocate(beta)
     deallocate(w)
     deallocate(work)
  ENDIF
  !
  if (saveonfile) then
     close (iunmix, status='keep')
     deallocate(dv)
     deallocate(df)
  else
     inext = iter - ( (iter - 1) / n_iter) * n_iter
     call DCOPY (ndim, vout, 1, df (1, inext), 1)
     call DCOPY (ndim, vinsave, 1, dv (1, inext), 1)
     deallocate(vinsave)
  endif
  !
  call stop_clock ('mix_pot')
  !
  return
  !
end subroutine mix_potential_core

SUBROUTINE setmixout(in1, in2, mix, dvscfout, dbecsum, ndim, flag )
 !
 USE kinds,    ONLY : DP
 USE mp_bands, ONLY : intra_bgrp_comm
 USE mp,       ONLY : mp_sum
 !
 IMPLICIT NONE
 INTEGER :: in1, in2, flag, ndim, startb, lastb
 COMPLEX(DP) :: mix(in1+in2), dvscfout(in1), dbecsum(in2)
 !
 CALL divide (intra_bgrp_comm, in2, startb, lastb)
 ndim=lastb-startb+1
 !
 IF (flag==-1) THEN
    mix(1:in1)=dvscfout(1:in1)
    mix(in1+1:in1+ndim)=dbecsum(startb:lastb)
 ELSE
    dvscfout(1:in1)=mix(1:in1)
    dbecsum=(0.0_DP,0.0_DP)
    dbecsum(startb:lastb)=mix(in1+1:in1+ndim)
    CALL mp_sum(dbecsum, intra_bgrp_comm)
 ENDIF
 !
 RETURN
 !
END SUBROUTINE setmixout
