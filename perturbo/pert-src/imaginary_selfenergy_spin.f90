!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  calculate the imaginary part of the electronic self-energy
!
!  NOTE: We use openmp parallel for qset.
!  NOTE: Segment fault may occur if openmp private space runs out.
!    we allocated memory for large arrays in each threads. however, the 
!    private stack space each thread has is limited (defined by OMP_STACKSIZE)
!    default value is implememention dependent and is usually quite samll. 
!    Segmentation fault may appear if run of of space. We can specify the size 
!    of the private space by 'export OMP_STACKSIZE=xxK (or xxM, xxG)'
!    or use openmp runtime subroutines: kmp_{set, get}_stacksize_s().
!
! Maintenance:
!  03/16/2020  Add spin-flip routines
!===============================================================================

module imaginary_selfenergy_spin
   use pert_const, only: dp, pi, czero, cone
   use qe_mpi_mod, only: ionode, stdout, mp_sum, inter_pool_comm, npool, &
      mp_split_pools, distribute_points
   use pert_param, only: delta_smear, phfreq_cutoff, prefix
   use pert_utils, only: abs2, bose, fermi, gauss, match_table
   use pert_output,only: progressbar_init, progressbar, set_progress_step, stopwatch
   use polar_correction, only: eph_wan_longrange
   use band_structure, only: electron_wann, solve_eigenvalue_vector, diagonalize_spinmatrix
   use phonon_dispersion, only: lattice_ifc, solve_phonon_modes
   use elphon_coupling_matrix, only: elph_mat_wann, eph_fourier_el, eph_fourier_elph, &
      eph_transform_fast, eph_transform_polar
   implicit none

   public :: calc_selfel_spin
contains

subroutine calc_selfel_spin &
      (elph, el, ph, kpts, qset, qwt, tmpr, ef, bmin, enk, imsgm, snk, imsgm_spin, rm_pol)
   implicit none
   type(electron_wann), intent(in) :: el
   type(lattice_ifc),   intent(in) :: ph
   type(elph_mat_wann), intent(in) :: elph
   !tmpr: temperature array; ef: chemical potential array; their size is the same.
   !qwt: q-points (qset) weight. kpts: k-points in the current pool. kpts(3,:)
   real(dp), intent(in) :: kpts(:,:), qset(:,:), qwt(:), tmpr(:), ef(:)
   !compute imsigma for bands start from bmin.
   integer, intent(in) :: bmin 
   !imsgm(nband, nmod, ntmpr, nkpt), enk(nband, nkpt)
   real(dp), intent(out) :: enk(:,:), imsgm(:,:,:,:) 
   !imsgm_spin(nband, nmod, ntmpr, nkpt), wnk(nband, nkpt)
   real(dp), intent(out) :: snk(:,:), imsgm_spin(:,:,:,:,:) 
   logical, intent(in), optional :: rm_pol !if true, (g^2-g_pol^2) is used.
   ! local variables
   integer :: mem_max, step_max, nk_pool, qst_pool, qend_pool, nq_pool, nstep_k, ikc
   integer :: step_q, nstep_q, istep_k, kst, kend, step_k, n, icount, iter
   integer :: numk, numq, ntmpr, numb, ik, iq, it, im, ib, jb, i, bmax, ndim(4), neph
   real(dp) :: xk(3), xq(3), xkq(3), wcut, e_thr, ph_n, el_fkq, dt1, dt2, scat
   
   real(dp) :: ek(el%nb), wq(ph%nm), ekq(el%nb), sk(el%nb), skq(el%nb)
   complex(dp) :: gpol(ph%nm), ukq(el%nb, el%nb), uq(ph%nm, ph%nm)
   logical,  allocatable :: ltable(:,:,:)
   real(dp), allocatable :: g2(:,:,:)
   complex(dp), allocatable :: uk(:,:,:), gkq(:,:,:), g_kerp(:,:,:,:,:,:)
   integer, pointer :: ik_loc(:)

   call start_clock('calc_selfel_spin')
   enk = 0.0_dp;  imsgm = 0.0_dp;  ndim = shape(imsgm)
   snk = 0.0_dp;  imsgm_spin = 0.0_dp
   !#. of selected bands, #.temperature, #. of kpts, number of q-points
   numb = ndim(1);  ntmpr = ndim(3);  numk = ndim(4);  numq = size(qset,2)
   bmax = bmin + numb -1
   !array sanity check
   if(size(tmpr).ne.ntmpr .or. size(ef).ne.ntmpr .or. size(kpts,2).ne.numk &
      .or. size(enk,1).ne.numb .or. size(enk,2) .ne. numk) &
      call errore('calc_selfel_spin','arguments dimension mismatch.',1)
   if(bmax > elph%nb) &
      call errore('calc_selfel_spin','bmin and imsgm mismatch',1)
   if(elph%nb .ne. el%nb .or. elph%na .ne. ph%na) &
      call errore('calc_selfel_spin','array dimension mismatch',1)
   if(present(rm_pol))then
      if (rm_pol .and. (.not. elph%lpol)) &
         call errore('calc_selfel_spin',"polar_split='rmpol' is for polar material",1)
   endif

   e_thr = delta_smear*3.0_dp  !exp(-3*3) ~ 10^-4
   wcut  = phfreq_cutoff
   neph = el%nb * el%nb * ph%na * 3
   !max memory per process for g_kerp, hard-coded here, not sure what is the optimal value
   mem_max  = 1024*1024*1024/2 !8GB
   step_max = mem_max / ( neph * elph%max_nrp )
   
   ik_loc => null()
   if(numk > 2*step_max*npool) then
      !interleave distribution of kpoints over pools, for better load balance
      call distribute_points(numk, nk_pool, ik_loc, .true.)
      !call mp_split_pools(numk, kst_pool, kend_pool, nk_pool)
      qst_pool = 1;     qend_pool = numq;    nq_pool = numq
   else
      call distribute_points(numk, nk_pool, ik_loc, .false.)
      !kst_pool = 1;     kend_pool = numk;    nk_pool = numk
      !distrbite q among pools
      call mp_split_pools(numq, qst_pool, qend_pool, nq_pool)
   endif
   !
   call set_progress_step(nk_pool, step_k, nstep_k, step_max)
   if(nstep_k >= 8) then
      !tracking progess on k only
      step_q = nq_pool;    nstep_q = 1
   else
      !track step on k and q combined
      call set_progress_step(nq_pool, step_q, nstep_q)
   endif
   !output debug info
   !write(stdout,'(5x, 3(a,i8,1x))') &
   !   'step_k:', step_k, 'nstep_k:', nstep_k, 'step_max:', step_max
   
   !allocate work sapce
   allocate( uk(el%nb, el%nb, step_k), &
      g_kerp(3, elph%max_nrp, elph%nb, elph%nb, elph%na, step_k) )

   if(ionode) call progressbar_init('Im[Sigma]:')
   icount = 0
   do istep_k = 1, nstep_k
      kst = (istep_k-1)*step_k + 1
      kend = min(kst+step_k-1, nk_pool)
      
      uk = czero;  g_kerp = czero
!$omp parallel do schedule(guided) default(shared) private(n, ikc, ik, xk, ek, sk)
      do n = kst, kend
         ikc = n - kst + 1
         !electronic wavefunction at ik
         ik = ik_loc(n)
         xk = kpts(:, ik)
         call solve_eigenvalue_vector(el, xk, ek, uk(:,:,ikc))
         enk(:,ik) = ek(bmin:bmax)
         call diagonalize_spinmatrix(el, xk, ek, uk(:,:,ikc), sk)
         snk(:,ik) = sk(bmin:bmax)
         call eph_fourier_el(elph, xk, g_kerp(:,:,:,:,:,ikc))
      enddo
!$omp end parallel do

      ! 
!$omp parallel default(shared) private(iq, xq, wq, uq, gpol, i, n, ik, ikc, xkq, &
!$omp& ekq, skq, ukq, ltable, gkq, g2, im, ib, jb, dt1, dt2, it, ph_n, el_fkq, scat, iter)
      allocate( ltable(el%nb, numb, ph%nm) )
      allocate( gkq(el%nb, el%nb, ph%nm), g2(el%nb, el%nb, ph%nm) )
!$omp do schedule(guided)
      do iq = qst_pool, qend_pool
         xq = qset(1:3, iq)
         !get phonon frequcies and eigen-displacement at iq
         call solve_phonon_modes(ph, xq, wq, uq)
         !get e-ph matrix elements in wannier gauge and cartesian coord.
         gpol = czero
         if(elph%lpol) call eph_wan_longrange(elph%pol, xq, gpol, uq)
         
         do n = kst, kend
            ikc = n - kst + 1
            !electronic wavefunction at ik
            ik = ik_loc(n)
            xkq = kpts(:,ik) + xq
            !get electronic wavefunction at ikq
            call solve_eigenvalue_vector(el, xkq, ekq, ukq)
            !check if there is any scattering channel perserve energy conservation
            call match_table(ekq, enk(:,ik), wq, wcut, e_thr, ltable)
            if(.not. any(ltable)) cycle
            call diagonalize_spinmatrix(el, xkq, ekq, ukq, skq)
            !
            call eph_fourier_elph(elph, xq, g_kerp(:,:,:,:,:,ikc), gkq)
            !transfor to phonon mode and bloch gauge.
            call eph_transform_fast(elph, uq, uk(:,:,ikc), ukq, gkq, gpol)
            !compute |g|^2
            g2 = abs2( gkq )
            if(present(rm_pol)) then
               if(rm_pol .and. elph%lpol) then
                  call eph_transform_polar(elph, gpol, uk(:,:,ikc), ukq, gkq)
                  g2 = g2 - abs2(gkq)
              endif
            endif

            do im = 1, ph%nm
            ! compute imsgm for bands between bmin:(bmin+numb-1)
            do  i = 1, numb
            do jb = 1, el%nb
               if(.not. ltable(jb,i,im)) cycle
               dt1 = gauss( delta_smear, enk(i,ik) - ekq(jb) + wq(im) )
               dt2 = gauss( delta_smear, enk(i,ik) - ekq(jb) - wq(im) )
               ib = i + bmin - 1
            
               do it = 1, ntmpr
                  ph_n = bose(tmpr(it), wq(im)) ! N_q\mu
                  el_fkq = fermi( tmpr(it), ekq(jb)-ef(it) )
                  scat = (el_fkq + ph_n)*dt1 + (1.0_dp - el_fkq + ph_n)*dt2
                  !NOTE: the order of jb and ib are important! Tr[g^\dagger g]
                  scat = scat * pi * qwt(iq) * g2(jb, ib, im) / (2.0_dp*wq(im))
!!NB: only one thread can update this shared variable at a time. negligible overhead.
!$omp atomic update
                  imsgm(i, im, it, ik) = imsgm(i, im, it, ik) + scat
                  !
!                  if(snk(i,ik)*skq(jb).lt.0.0_dp) then
!$omp atomic update
                     imsgm_spin(i, im, it, ik, 1) = imsgm_spin(i, im, it, ik, 1) + scat*(snk(i,ik)-skq(jb))/snk(i,ik)
               enddo
            enddo; enddo; enddo
         enddo

         !track progress
!$omp atomic update
         icount = icount + 1
         iter = icount
         if( mod(iter, step_q*nstep_k).eq.0 .or. iter.eq.(nq_pool*nstep_k) ) then
!$omp critical (calc_selfel_spin_progress)
            write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/(nq_pool*nstep_k), '%'
!$omp end critical (calc_selfel_spin_progress)
         endif

      enddo
!$omp end do
      deallocate(gkq, g2, ltable)
!$omp end parallel
   enddo
   deallocate(uk, g_kerp)
   call mp_sum(imsgm, inter_pool_comm)
   if(numk .ne. nk_pool) call mp_sum(enk, inter_pool_comm)
   !
   call mp_sum(imsgm_spin, inter_pool_comm)
   if(numk .ne. nk_pool) call mp_sum(snk, inter_pool_comm)
   !release space
   if( associated(ik_loc) ) deallocate( ik_loc )

   call stop_clock('calc_selfel_spin')
end subroutine calc_selfel_spin

end module imaginary_selfenergy_spin
