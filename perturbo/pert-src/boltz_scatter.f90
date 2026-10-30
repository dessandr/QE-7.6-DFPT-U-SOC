!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
! Copyright (C) 2021-2023 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Ivan Maliyov
!                         Dhruv Desai, Sergio Pineda Flores, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!
! Maintenance:
!===============================================================================

#include "pert_config.h"

module boltz_scatter
   use kinds, only: i8b
   use boltz_grid, only: grid
   use boltz_utils, only: idx2band
   use yaml_utils, only: ymlout, output_scatter_info_yaml
   use pert_const, only: dp, czero
   use pert_data,  only: epr_fid, qc_dim, kc_dim
   use pert_param, only: phfreq_cutoff, delta_smear, use_mem, load_scatter_eph, ntemper
   use band_structure, only: electron_wann, solve_eigenvalue_vector
   use phonon_dispersion, only: lattice_ifc, solve_phonon_modes
   use qe_mpi_mod, only: stdout, ionode, mp_barrier, inter_pool_comm, npool, my_pool_id, mp_sum
   use elphon_coupling_matrix, only: elph_mat_wann, release_elph_mat_wann, &
      init_elph_mat_wann, eph_fourier_el_para, eph_fourier_elph, eph_transform_fast
   use boltz_scatter_sizest
#if defined(STORE_G2DT)
   use pert_const, only: pi
   use pert_utils, only: gauss
#endif

   implicit none
   public

   !=========================================================================
   ! Types for kq_grid data structure

   type, private :: kq_channel
      integer :: ik
      integer :: nkq
      !
      integer, pointer :: iq(:) => null()
      integer, pointer :: ikq(:) => null()
      integer, pointer :: nchl(:) => null()
   end type

   type :: phonon_grid
      integer :: ndim(3)  !size of the grid
      !
      integer :: npts  !number of points selected
      integer, pointer :: list(:) => null()
      real(dp), pointer :: freq(:,:) => null()
      !
      real(dp):: weight  ! 1/(ndim(1) * ndim(2) * ndim(3))
      real(dp), pointer :: disp(:,:,:) => null() !> calculate the bose distribution for saving runtime
   end type

   !=========================================================================
   ! Types for hybrid scattering data structures

   !user-defined data structure for scattering channel. based on kgrid
   type :: t_scatter
      !location of k and kq in kgrid%kpt: xk(1:3) = num2kpt(kgrid%kpt(ik), ndim)
      !NOTE: (k, kq) and (kq, k) are equivalent, only one of them are stored.
      integer :: ik
      integer :: ikq
      !location of q = (k+q) - k in qgrid%list: xq(1:3) = num2kpt(qgrid%list(iq), ndim)
      integer :: iq

#if defined(SCAT_FWD)
      ! It's possible to eliminate one of these two integer values, and just
      ! store the end-index of each kq-pair's sequence of channels.  But,
      ! this approach is just simpler to parallelize with OpenMP/OpenACC.

      integer(kind=i8b) :: chl_start_index

      !number of (mkq, nk, mu) terms, mu is phonon modes index
      !only g_m,n,mu of (k, kq) satisfy the following condition are stored.
      !abs( abs(e_nk - e_mkq) - w_qmu ) < delta_E. (delta_E = e_thr*broadening)
      integer :: nchl
#endif
   end type t_scatter

   ! A "flattened" and denormalized scattering-channel array for rapid,
   ! (hopefully) OpenMP/OpenACC-friendly traversal.  The main difference
   ! between t_scatt and t_flatscat is that t_flatscat doesn't have any
   ! sub-arrays to access.  The downside is that the (ik, ikq, iq) values
   ! are repeated across multiple rows.
   type :: t_scatter_channel
#if defined(SCAT_REV)
      !location of k and kq in kgrid%kpt: xk(1:3) = num2kpt(kgrid%kpt(ik), ndim)
      !NOTE: (k, kq) and (kq, k) are equivalent, only one of them are stored.
      integer :: scat_index
#endif

      integer :: bands_index   ! (m, n, im) packed into a single integer

#if !defined(STORE_G2DT)
      real(dp) :: eph_g2       ! g2
#else
      real(dp) :: g2_dt1       ! g2*dt1
      real(dp) :: g2_dt2       ! g2*dt2
#endif
   end type t_scatter_channel

   !=========================================================================
   ! Data structures for dynamics simulation

   integer, private, save :: nk_loc
   type(kq_channel), allocatable, private, save :: kq_pair(:)
   !
   type(phonon_grid), save :: qgrid
   !
   integer, save :: num_kq_pairs         ! replaces num_scat
   integer(kind=i8b), save :: num_scatter_channels ! replaces num_flatscat
   !
   type(t_scatter), allocatable, save :: scatter(:)
   type(t_scatter_channel), allocatable, save :: scatter_channels(:)
contains


subroutine boltz_scatter_setup(kg, el, ph, qdim)
   use pert_param, only: ntemper, temper, efermi
   use pert_utils, only: bose
   implicit none
   type(grid), intent(in) :: kg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   !dimension for the q-grid. q-should be commensurate with the k-grid
   integer, intent(in) :: qdim(3)
   !local
   integer :: nkpt, i, ik, collect_pools(npool), nq_col(npool)
   integer :: numq, numb, itemp, iq
   integer :: prev_kqp
   integer(kind=i8b) :: prev_sch
   logical :: overflow
   real(dp) :: wcut, e_thr
   type(elph_mat_wann) :: elph

   !init
   nkpt = kg%nk
   wcut  = phfreq_cutoff
   e_thr = delta_smear*3.0_dp  !exp(-3*3) ~ 10^-4

   !interleave distribution of kpoints over pools, for better load balance
   nk_loc = nkpt/npool + merge(1, 0, my_pool_id < mod(nkpt, npool))
   allocate( kq_pair(nk_loc) )
   !store the allocated kpts (location in kg%kpt) of the current pool
   i = 0
   do ik = my_pool_id+1, nkpt, npool
      i = i + 1
      kq_pair(i)%ik = ik
   enddo
   !debug
   if(i .ne. nk_loc) call errore('boltz_scatter_setup', 'i .ne. nk_loc', 1)
   !
   qgrid%ndim(1:3) = qdim(1:3)
   qgrid%weight = 1.0_dp / dble(qdim(1)) / dble(qdim(2)) / dble(qdim(3))
   !generate (k,q) pair (initialized kq_pair, qgrid)
   call setup_kq_pair(kg, el, ph, wcut, e_thr, nk_loc, kq_pair, qgrid)

   !initialize the phonon distribution
   numb = size(qgrid%freq,1)
   numq = size(qgrid%freq,2)
#if defined(CAL_FERMI_BOSE)
   allocate(qgrid%disp(numb,numq,ntemper))
   do iq = 1, numq
      do itemp = 1, ntemper
         qgrid%disp(:,iq,itemp) = bose(temper(itemp), qgrid%freq(:,iq)) 
      enddo
   enddo
#endif

   !
   ! compute and save g2 to disk
   if(.not. load_scatter_eph) then
      call init_elph_mat_wann(epr_fid, kc_dim, qc_dim, elph)
      !compute g2 and sove them to temperary files
      call compute_scatter_eph_g2 &
         (kg, el, ph, elph, qgrid, wcut, e_thr, use_mem)
      !release space
      call release_elph_mat_wann(elph)
   endif

   !!! END KQ-GRID SETUP !!!

   !
   ! Compute the total number of (k,q) pairs and the total number of scatter
   ! channels on this MPI node.
   num_kq_pairs = 0
   num_scatter_channels = 0
   overflow = .false.
   do i = 1, nk_loc
      if (kq_pair(i)%nkq .eq. 0) cycle

      ! For detecting overflows
      prev_kqp = num_kq_pairs
      prev_sch = num_scatter_channels

      num_kq_pairs = num_kq_pairs + kq_pair(i)%nkq
      num_scatter_channels = num_scatter_channels + sum( kq_pair(i)%nchl(:) )

      if (prev_kqp >= 0 .and. num_kq_pairs < 0) then
         write(stdout, *) 'ERROR:  Detected overflow on num_kq_pairs!'
         flush(stdout)
         overflow = .true.
      endif

      if (prev_sch >= 0 .and. num_scatter_channels < 0) then
         write(stdout, *) 'ERROR:  Detected overflow on num_scatter_channels!'
         flush(stdout)
         overflow = .true.
      endif
   enddo

   if (overflow) call errore('setup_scatter', 'Overflows detected!', 1)

   !
   ! Collect the number of q-points and (k,q)-pairs in each pool so we can
   ! print it and save it to the YAML file.
   nq_col = 0
   nq_col( my_pool_id+1 ) = qgrid%npts
   call mp_sum( nq_col, inter_pool_comm )
   !
   collect_pools = 0
   collect_pools( my_pool_id+1 ) = num_kq_pairs
   call mp_sum( collect_pools, inter_pool_comm )
   !output info
   write(stdout, '(5x, a)') "No. of q-points and (k,q) pairs on each processor:"
   do i = 1, npool
      write(stdout, '(5x, a, i4, a, i9, a, i12)')  &
         "Proc.", i, ':  #.q ', nq_col(i), ' |  #.(k,q)', collect_pools(i)
   enddo

   write(stdout, '(/, 5x, a, i23)') "Total No. of q-points:", sum(nq_col)
   write(stdout, '(5x, a, i20, /)') "Total No. of (k,q) pairs:", sum(INT8(collect_pools))
   flush(stdout)

   ! output to the YAML file
   call output_scatter_info_yaml(npool, nq_col, collect_pools)

   call setup_scatter(kg, el, ph%nm)

   if (ionode) then
      ! Log data-structure allocation details into the YAML output file.
      write(ymlout, '(/,a)') 'boltz_scatter memory:'
      call output_kq_pair_alloc_stats_to_yaml()
      call output_scatter_alloc_stats_to_yaml()
      call output_scatter_channels_alloc_stats_to_yaml()
   endif

   do i = 1, nk_loc
      if( kq_pair(i)%nkq == 0 ) cycle
      deallocate( kq_pair(i)%iq, kq_pair(i)%ikq, kq_pair(i)%nchl )
   enddo
   deallocate( kq_pair )

end subroutine boltz_scatter_setup


! Generate the initial kq-pair data structure that drives the rest of the
! dynamics simulation.
subroutine setup_kq_pair(kg, el, ph, wcut, e_thr, nk_loc, kq_pair, qg)
   use pert_utils, only: match_table
   use boltz_utils, only: kpts_minus, num2kpt, band2idx, kpts_plus
   implicit none
   type(grid), intent(in) :: kg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   real(dp), intent(in) :: wcut, e_thr
   !
   integer, intent(in) :: nk_loc
   type(kq_channel), intent(inout) :: kq_pair(:)
   !
   type(phonon_grid), intent(inout) :: qg
   !
   !local variables
   integer :: nkpt, nmod, maxq, nkq, nq_tmp, i, ik, ikq, j, qidx, iq, k
   logical :: ltable(kg%numb, kg%numb, ph%nm)
   integer, allocatable :: qcollect(:)
   real(dp) :: xq(3)
   !
   integer,  pointer :: ql_tmp(:), q_check(:)
   real(dp), pointer :: qw_tmp(:,:)
   type(kq_channel), pointer :: kq_tmp(:)

   call start_clock('set_kq_pair')
   !init
   nkpt = kg%nk
   nmod = ph%nm

   allocate( kq_tmp(nk_loc) )
   !find iq
   !maxq = qg%ndim(1)*qg%ndim(2)*qg%ndim(3) - 1
   maxq = qg%ndim(3) * (1 +  qg%ndim(2) * (1 + qg%ndim(1)) ) / 2
   allocate(qcollect(0:maxq))
   qcollect = 0
   !(k, k') and (k', k) are connected, only one of them is stored.
   !to balance load, each kpoint is associated with equal number of k+q.
   !so for each ik in grid%kpt, the (k+q) are choose as ik+1, ik+2, ...., ik+(nkpt-1)/2
   !  if (ik + i) is larger than nkpt, restart from 1. for example, nkpt + 5 -> 5
   !
   !initialize and pre-allocate space
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      kq_tmp(i)%ik = ik
      !assign (k+q) evenly for each k, for load balance
      nkq = (nkpt-1)/2 + merge(1, 0, ((mod(nkpt-1,2)>0) .and. ik<=nkpt/2))
      kq_tmp(i)%nkq = nkq

      allocate( kq_tmp(i)%iq(nkq) )
      kq_tmp(i)%iq(:) = -1
   enddo
   !
   !collect q-points and store them in qg%list.
!$omp parallel do schedule(guided) default(shared) private(i, ik, nkq, j, ikq, qidx)
   do i = 1, nk_loc
      ik = kq_tmp(i)%ik
      nkq = kq_tmp(i)%nkq

      !loop over (k+q) points assigned to this current k-points
      do j = 1, nkq
         ikq = mod(ik+j-1, nkpt) + 1
         !store (k+q) - k in kq_pair
         !(since q and -q has the same freqencies, we only store one of them)
         ! if xq = (ikq - iq) is not on the q-grid (qg), then qidx = -1.
         qidx = min( kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim), &
               kpts_minus(kg%kpt(ik), kg%kpt(ikq), kg%ndim, qg%ndim) )
         !debug
         !if(qidx > maxq) call errore('setup_kq_pair', 'qidx > maxq', 1)
         !
         !if this xq is on the q-grid qg.
         if(qidx > -1) then
            kq_tmp(i)%iq(j)  = qidx
!$omp atomic
            qcollect(qidx) = qcollect(qidx) + 1
         endif
      enddo
   enddo
!$omp end parallel do

   !collect q-points
   iq = 0
   do qidx = 0, maxq
      if(qcollect(qidx) > 0) then
         iq = iq + 1
         qcollect(qidx) = iq
      endif
   enddo
   nq_tmp = iq

   allocate( ql_tmp(nq_tmp) )
   do qidx = 0, maxq
      iq = qcollect(qidx)
      if( iq > 0) ql_tmp( iq ) = qidx
   enddo

   !update kq_pair%iq
   do i = 1, nk_loc
      do j = 1, kq_tmp(i)%nkq
         qidx = kq_tmp(i)%iq(j)
         !qidx = -1 if xq is not on the q-grid.
         if(qidx > -1) then
            iq = qcollect(qidx)
            kq_tmp(i)%iq(j) = iq
         endif
      enddo
   enddo
   !release space
   deallocate(qcollect)

   !allocate space
   allocate(qw_tmp(nmod, nq_tmp))
   do i = 1, nk_loc
      nkq = kq_tmp(i)%nkq
      allocate( kq_tmp(i)%ikq(nkq)  )
      allocate( kq_tmp(i)%nchl(nkq) )
      !initialize
      kq_tmp(i)%nchl(:) = 0
   enddo

   !compute phonon freq for qpts
!$omp parallel default(shared) private(iq, xq, i, ik, j, ikq, ltable)
!$omp do schedule(guided)
   do iq = 1, nq_tmp
      xq  = num2kpt( ql_tmp(iq),  qg%ndim )
      call solve_phonon_modes(ph, xq, qw_tmp(:,iq))
   enddo
!$omp end do

!$omp do schedule(guided)
   do i = 1, nk_loc
      ik = kq_tmp(i)%ik
      !nkq = (nkpt-1)/2 + merge(1, 0, ((mod(nkpt-1,2)>0) .and. ik<=nkpt/2))
      do j = 1, kq_tmp(i)%nkq
         ikq = mod(ik+j-1, nkpt) + 1
         kq_tmp(i)%ikq(j) = ikq
         !
         iq = kq_tmp(i)%iq(j)
         if(iq > 0) then
            call match_table(kg%enk(:,ikq), kg%enk(:,ik), qw_tmp(:,iq), wcut, e_thr, ltable)
            kq_tmp(i)%nchl(j) = count(ltable)
         endif
      enddo
   enddo
!$omp end do
!$omp end parallel

   !update qg%list
   allocate( q_check(nq_tmp) )
   q_check = 0
   !collect (k,q) pair that has at least one valid scattering channel (nchl > 0).
   do i= 1, nk_loc
      nkq = count( kq_tmp(i)%nchl(:) > 0 )
      kq_pair(i)%nkq = nkq
      ! only needed when nkq > 0
      if( nkq > 0) then
         allocate( kq_pair(i)%iq(nkq), kq_pair(i)%ikq(nkq), kq_pair(i)%nchl(nkq) )
         !
         k = 0
         do j = 1, kq_tmp(i)%nkq
            if( kq_tmp(i)%nchl(j) > 0 ) then
               k = k + 1
               !
               kq_pair(i)%ikq(k) = kq_tmp(i)%ikq(j)
               kq_pair(i)%nchl(k) = kq_tmp(i)%nchl(j)
               iq = kq_tmp(i)%iq(j)
               kq_pair(i)%iq(k) = iq
               !
               q_check(iq) = q_check(iq) + 1
            endif
         enddo
      endif
      !
      deallocate(kq_tmp(i)%iq, kq_tmp(i)%ikq, kq_tmp(i)%nchl)
   enddo
   deallocate( kq_tmp )

   !if number of q-points is reduced.
   if( count(q_check > 0) .ne. nq_tmp ) then
      k = 0
      do i = 1, nq_tmp
         if(q_check(i) > 0) then
            k = k + 1
            q_check(i) = k
         endif
      enddo
      qg%npts = k

      !update qg%list
      allocate( qg%list( qg%npts ), qg%freq(nmod, qg%npts) )
      do i = 1, nq_tmp
         if(q_check(i) > 0) then
            qg%list( q_check(i) ) = ql_tmp(i)
            qg%freq(:, q_check(i)) = qw_tmp(:, i)
         endif
      enddo
      deallocate( ql_tmp, qw_tmp )

      !update kq_pair
      do i = 1, nk_loc
         do j = 1, kq_pair(i)%nkq
            iq = kq_pair(i)%iq(j)
            kq_pair(i)%iq(j) = q_check(iq)
         enddo
      enddo
   else
      qg%npts = nq_tmp
      qg%list => ql_tmp
      qg%freq => qw_tmp
   endif
   deallocate( q_check )

   call stop_clock('set_kq_pair')
end subroutine setup_kq_pair


! Generate an HDF5 file for the current MPI pool (i.e. MPI node), and generate
! the eph and g2 values into the HDF5 file.
!
! Each MPI node will generate its own HDF5 file, so the working directory will
! end up with N files for N nodes.  These files can then be reloaded on
! subsequent perturbo runs, as long as the number of MPI nodes doesn't change.
subroutine compute_scatter_eph_g2(kg, el, ph, elph, qg, wcut, e_thr, use_mem)
   use pert_param, only: prefix, tmp_dir
   use pert_utils, only: abs2, match_table
   use polar_correction,  only: eph_wan_longrange
   use boltz_utils, only: kpts_minus, num2kpt, band2idx, kpts_plus
   use pert_output, only: progressbar_init, progressbar
   use hdf5_utils
   implicit none
   type(grid), intent(in) :: kg
   type(phonon_grid), intent(in) :: qg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   type(elph_mat_wann), intent(in) :: elph
   !
   logical, intent(in) :: use_mem
   real(dp), intent(in) :: wcut, e_thr
   !
   integer :: nmod, qidx, numq
   integer :: i, j, ik, iq, ikq, it, im, n, m, nkq, tot_chnl, bands(3)
   real(dp) :: xq(3), xk(3), xkq(3), wq(ph%nm), ek(el%nb), ekq(el%nb)
   complex(dp) :: uqt(ph%nm, ph%nm), uk(el%nb, el%nb), ukq(el%nb, el%nb), gpol(ph%nm)
   logical :: ltable(kg%numb, kg%numb, ph%nm)
   !
   integer, allocatable :: ist(:), bnd_idx(:)
   real(dp), allocatable :: eph_g2(:), g2(:,:,:)
   complex(dp), allocatable :: uq(:,:,:), gp(:,:), g_kerp(:,:,:,:,:), gkq(:,:,:)
   !
   integer(HID_T) :: file_id
   character(len=120) :: fname, dset_name, msg
   character(len=6), external :: int_to_char

   call start_clock('compute_g2')
   !sanity check
   if(elph%nb .ne. el%nb .or. elph%na .ne. ph%na) &
      call errore('compute_scatter_eph_g2','array dimension mismatch',1)
   nmod = ph%nm
   numq = qg%npts

   if( use_mem ) allocate( uq(nmod, nmod, numq) )
   if(elph%lpol) allocate( gp(nmod, numq) )
   if(elph%lpol .or. use_mem) then
!$omp parallel do schedule(guided) default(shared) private(iq, xq, wq, uqt)
      do iq = 1, numq
         xq = num2kpt(qg%list(iq), qg%ndim)
         call solve_phonon_modes(ph, xq, wq, uqt)
         if(elph%lpol) call eph_wan_longrange(elph%pol, xq, gp(:, iq), uqt)
         if(use_mem)  uq(:,:, iq) = uqt(:,:)
      enddo
!$omp end parallel do
   endif

   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_eph_g2_p"
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   call hdf_open_file(file_id, trim(fname), status='NEW')

   ! allocate work space
   allocate( g_kerp(3, elph%max_nrp, elph%nb, elph%nb, elph%na) )
   ! print out progress info
   call progressbar_init('Computing EPhMatrix:')
   ! main loop
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      nkq = kq_pair(i)%nkq
      ! skip if nkq = 0
      if( nkq < 1 ) cycle

      !kpts in crystal coordinate
      xk = num2kpt( kg%kpt(ik), kg%ndim )
      call solve_eigenvalue_vector(el, xk, ek, uk)
      call eph_fourier_el_para(elph, xk, g_kerp)
      !
      allocate( ist(nkq) )
      ist(1) = 0
      do j = 1, nkq-1
         ist(j+1) = ist(j) + kq_pair(i)%nchl(j)
      enddo
      !allocate space
      tot_chnl = sum(kq_pair(i)%nchl(:))
      allocate( bnd_idx(tot_chnl), eph_g2(tot_chnl) )
      bnd_idx = 0;   eph_g2 = 0.0_dp

!$omp parallel default(shared) private(j, iq, ikq, xkq, ekq, ukq, xq, wq, uqt, &
!$omp& qidx, gpol, gkq, g2, ltable, it, im, n, m, bands)
      allocate( gkq(el%nb, el%nb, nmod), g2(el%nb, el%nb, nmod) )
!$omp do schedule(guided)
      do j = 1, nkq
         iq = kq_pair(i)%iq(j)
         ikq = kq_pair(i)%ikq(j)
         !
         xkq = num2kpt(kg%kpt(ikq), kg%ndim)
         call solve_eigenvalue_vector(el, xkq, ekq, ukq)
         !
         xq = num2kpt( qg%list(iq), qg%ndim )
         if(use_mem) then
            uqt(:,:) = uq(:,:, iq)
         else
            call solve_phonon_modes(ph, xq, wq, uqt)
         endif
         !
         gpol = (0.0_dp,0.0_dp)
         !NOTE: qg%list(iq) could be (k+q)-k or k-(k+q)
         !  Here we need to check which case it is
         qidx= kpts_minus( kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim)
         ! iq -> (k+q) - k
         if( qidx .eq. qg%list(iq) ) then
            !uqt(:,:) = uq(:,:, iq)
            if(elph%lpol) gpol(:) = gp(:,iq)
         ! iq -> k - (k+q) .e.g -q
         else if( kpts_plus(qidx, qg%list(iq), qg%ndim, qg%ndim) .eq. 0 ) then
            !NOTE: it's important to use xq = -xq; alternatively we can
            ! use xq'=num2kpt(qidx), but xq' might be equal to xq + G
            ! (reciprocal lattice vectors). a small numerical noise due
            ! to G might break the relation u(q)* = u(-q);
            xq  = -xq
            !u(q)* = u(-q);  g^L(q)* = g^L(-q)
            uqt(:,:) = dconjg( uqt(:,:) )
            if(elph%lpol) gpol(:) = dconjg( gp(:,iq) )
         else
            !write(stdout,*) kpts_plus(qidx, qg%list(iq), qg%ndim, qg%ndim), qg%list(iq), qidx
            write(msg,'(2(3f12.6,2x))') num2kpt(qidx, qg%ndim), xq
            call errore('compute_scatter_eph_g2','failed: ' // trim(msg), 1)
         endif

         !get e-ph matrix elements in wannier gauge and cartesian coord.
         call eph_fourier_elph(elph, xq, g_kerp, gkq)
         call eph_transform_fast(elph, uqt, uk, ukq, gkq, gpol)
         !compute |g|^2
         g2 = abs2( gkq )
         !
         call match_table(kg%enk(:,ikq), kg%enk(:,ik), qg%freq(:,iq), wcut, e_thr, ltable)
         !sanity check
         if( kq_pair(i)%nchl(j) .ne. count(ltable) ) &
            call errore('compute_scatter_eph_g2','kq_pair(i)%nchl(j) .ne. count(ltable)',1)

         it = ist(j)
         do im = 1, nmod
         do n = 1, kg%numb
         do m = 1, kg%numb
            if( ltable(m,n,im) ) then
               it = it + 1
               bands = (/m, n, im/)
               bnd_idx(it) = band2idx(bands, kg%numb)
               eph_g2(it) = g2(m+kg%bmin-1, n+kg%bmin-1, im) / (2.0_dp*qg%freq(im,iq))
            endif
         enddo; enddo; enddo
      enddo
!$omp end do
      deallocate(gkq, g2)
!$omp end parallel
      !
      ! output to hdf5 files
      dset_name = "bands_index_" // trim( int_to_char(i) )
      call hdf_write_dataset(file_id, trim(dset_name), bnd_idx)
      dset_name = "eph_g2_" // trim( int_to_char(i) )
      call hdf_write_dataset(file_id, trim(dset_name), eph_g2)
      ! show progress
      call progressbar(i, nk_loc)
      !
      deallocate(bnd_idx, eph_g2, ist)
   enddo
   call hdf_close_file(file_id)

   deallocate(g_kerp)
   if( use_mem ) deallocate( uq )
   if(elph%lpol) deallocate( gp )

   call mp_barrier(inter_pool_comm)
   call stop_clock('compute_g2')
end subroutine compute_scatter_eph_g2


! Open the HDF5 file generated for this MPI pool (i.e. MPI node) that
! contains the eph_g2 and bands_index values.  The file-handle is
! returned back in the file_id out-parameter.
!
! The subroutine doesn't report any errors because if the file can't
! be opened then the subroutine halts the program with an error.
subroutine open_pool_eph_g2_file(file_id)
   use hdf5_utils
   use pert_param, only: prefix, tmp_dir
   implicit none
   !
   ! For reading the HDF5 file data
   integer(HID_T), intent(out) :: file_id
   logical :: has_file
   character(len=120) :: fname
   character(len=6), external :: int_to_char

   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_eph_g2_p"
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   ! check if file exist
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('load_scatter_eph_g2', "Missing file: "//trim(fname), 1)
   !
   ! open file
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
end subroutine open_pool_eph_g2_file


! Read the eph_g2 array and bands_index array for a given k-grid location
! from the HDF5 file generated for this MPI pool (i.e. MPI node).  The
! file-handle and k-grid location are passed in as input-parameters.  The
! bnd_idx and eph_g2 pointers are expected to point to arrays that are
! large enough to receive the values stored in the file.
subroutine read_pool_eph_g2_data(file_id, i_kloc, bnd_idx, eph_g2)
   use hdf5_utils
   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: i_kloc
   integer, allocatable, intent(inout) :: bnd_idx(:)
   real(dp), allocatable, intent(inout) :: eph_g2(:)
   !
   character(len=120) :: dset_name
   character(len=6), external :: int_to_char

   dset_name = "bands_index_" // trim( int_to_char(i_kloc) )
   call hdf_read_dataset(file_id, trim(dset_name), bnd_idx)
   dset_name = "eph_g2_" // trim( int_to_char(i_kloc) )
   call hdf_read_dataset(file_id, trim(dset_name), eph_g2)
end subroutine


subroutine setup_scatter(kg, el, nmodes)
   use hdf5_utils
   use boltz_utils, only: idx2band
   use qe_mpi_mod, only: ionode, stdout
   use pert_param, only: prefix, tmp_dir, ph_mode_exclude_ranges
   implicit none

   type(grid), intent(in) :: kg
   type(electron_wann), intent(in) :: el
   integer, intent(in) :: nmodes

   !local
   integer :: i_scat, tot_chl
   integer(kind=i8b) :: i_scat_chl
   integer :: n, i, ik, nkq, ikq, iq, j, nchl, it, ib
#if defined(STORE_G2DT)
   ! Variables needed to compute -g2*dt1 and -g2*dt2
   integer :: m
   real(dp) :: wq, enk, emkq, g2, dt1, dt2
#endif
   !
   ! For reading the HDF5 file data
   integer(HID_T) :: file_id
   integer, allocatable :: bnd_idx(:) ! For reading band-indexes from HDF5 file
   real(dp), allocatable :: eph_g2(:) ! For reading eph-g2 values from HDF5 file
   !
   ! Exclude ph mode variables
   integer :: im, jm, im1, im2
   integer :: bands(3)
   logical :: include_phonon_modes(nmodes)

   call start_clock('scat_setup')
   write(stdout,'(5x,a20)') '>start scatter_setup'

   allocate( scatter(num_kq_pairs) )
   allocate( scatter_channels(num_scatter_channels) )

   ! TODO(donnie):  Remove this output?
#if defined(SCAT_FWD)
   write (*, '(5x, a)') ' * scatter will support forward-traversal'
#endif

#if defined(SCAT_REV)
   write (*, '(5x, a)') ' * scatter will support reverse-traversal'
#endif

   write (*, '(5x, a, i12)') ' * Number of kq-pairs:         ', num_kq_pairs
   write (*, '(5x, a, i12)') ' * Number of scatter-channels: ', num_scatter_channels

   !read g2 from files

   call start_clock("load_eph_g2")

   call open_pool_eph_g2_file(file_id)

   i_scat = 0
   i_scat_chl = 0
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      nkq = kq_pair(i)%nkq
      !skip if nkq = 0
      if(nkq < 1) cycle

      ! Arrays for reading values from HDF5 file
      tot_chl = sum( kq_pair(i)%nchl(:) )
      allocate( bnd_idx(tot_chl), eph_g2(tot_chl) )
      bnd_idx = 0
      eph_g2 = 0.0_dp
      it = 0

      ! Read bnd_idx and eph_g2 from HDF5 file
      call read_pool_eph_g2_data(file_id, i, bnd_idx, eph_g2)

      do j = 1, nkq
         i_scat = i_scat + 1
         iq = kq_pair(i)%iq(j)
         ikq = kq_pair(i)%ikq(j)
         !
         scatter(i_scat)%ik = ik
         scatter(i_scat)%iq = iq
         scatter(i_scat)%ikq = ikq

         nchl = kq_pair(i)%nchl(j)

#if defined(SCAT_FWD)
         ! must add 1 since i_scat_chl is incremented in the innermost loop
         scatter(i_scat)%chl_start_index = i_scat_chl + 1
         scatter(i_scat)%nchl = nchl
#endif

         do ib = 1, nchl
            i_scat_chl = i_scat_chl + 1
            it = it + 1
#if defined(SCAT_REV)
            scatter_channels(i_scat_chl)%scat_index = i_scat
#endif
            scatter_channels(i_scat_chl)%bands_index = bnd_idx(it)

#if defined(STORE_G2DT)
            ! Precompute and store the -g2*dt1 and -g2*dt2 terms, since they
            ! don't change from iteration to iteration.

            bands = idx2band( bnd_idx(it), kg%numb )
            m = bands(1);  n = bands(2);  im = bands(3)

            ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
            g2 = eph_g2(it)

            ! eigenvalue
            wq   = qgrid%freq(im, iq)
            enk  = kg%enk(n, ik)
            emkq = kg%enk(m, ikq)

            ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
            ! delta(Enk - Emkq + wq)
            dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
            ! delta(Enk - Emkq -wq)
            dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)

            scatter_channels(i_scat_chl)%g2_dt1 = g2 * dt1
            scatter_channels(i_scat_chl)%g2_dt2 = g2 * dt2
#else
            scatter_channels(i_scat_chl)%eph_g2 = eph_g2(it)
#endif
         enddo
      enddo
      deallocate(bnd_idx, eph_g2)
   enddo
   !

   call hdf_close_file(file_id)
   call mp_barrier( inter_pool_comm )
   call stop_clock("load_eph_g2")

   ! Set eph_g2 to zero for the phonon modes to exclude
   ! This is done in separate nested do loops to not to effect
   ! the efficiency for the case where no phonon modes have to be excluded.
   if( any(ph_mode_exclude_ranges(:, :) > 0) ) then

      call start_clock("exclude_ph_modes")

      write(stdout, '(/, 5x, a)') '>excluding phonon modes:'

      call interpret_ph_exclude_ranges(nmodes, include_phonon_modes)

      ! Iterate over all scatter-channels to zero out eph_g2 values for
      ! phonon modes to be excluded.
!$omp parallel do schedule(guided) default(shared) private(i_scat_chl, bands, im)
      do i_scat_chl = 1, num_scatter_channels
         bands = idx2band( scatter_channels(i_scat_chl)%bands_index, kg%numb )
         im = bands(3)

         ! Set the g2 elements to zero for the phonon modes to exclude
         if (.not. include_phonon_modes(im)) then
#if !defined(STORE_G2DT)
            scatter_channels(i_scat_chl)%eph_g2 = 0.0_dp
#else
            scatter_channels(i_scat_chl)%g2_dt1 = 0.0_dp
            scatter_channels(i_scat_chl)%g2_dt2 = 0.0_dp
#endif
         endif
      enddo
!$omp end parallel do

      call stop_clock("exclude_ph_modes")
   end if

   !
   write(stdout,'(5x,a20)') '>scatter_setup done.'
   call stop_clock('scat_setup')
   flush(stdout)

end subroutine setup_scatter


!> Based on the array from the input file ph_mode_exclude_ranges
!! build the include_phonon_modes bollean
!! array. If an element of this array is false, this phonon mode
!! will be omitted in the scat: g2 will be set to zero.
subroutine interpret_ph_exclude_ranges(nmodes, include_phonon_modes)
   use yaml_utils, only: ymlout
   use qe_mpi_mod, only: ionode, stdout
   use pert_param, only: ph_mode_exclude_ranges
   implicit none
   integer, intent(in) :: nmodes
   logical, intent(inout) :: include_phonon_modes(nmodes)
   !
   integer :: im, jm, im1, im2
   character(len=6), external :: int_to_char
   character(len=120) :: msg

   include_phonon_modes(:) = .true.

   do im = 1, size(ph_mode_exclude_ranges, 1)

      im1 = ph_mode_exclude_ranges(im, 1)
      if( im1 > 0 ) then

         im2 = ph_mode_exclude_ranges(im, 2)
         if ( im1 > nmodes .or. im2 > nmodes ) then
            write(msg, '(a,I4)') &
            'exclude mode indices must be smaller than the number of phonon modes:', nmodes
            call errore('load_scatter_eph_g2', trim(msg), 1)

         else if (im2 < im1) then
            write(msg, '(a,I2,a,I2,a)') &
            'ph_mode_exclude_ranges(',im, &
            ', 2) is smaller than ph_mode_exclude_ranges(',im, ', 1)'
            call errore('load_scatter_eph_g2', trim(msg), 1)
         end if

         do jm = im1, im2
            include_phonon_modes(jm) = .false.
         end do

      end if
   end do

   ! output which phonon modes will be included or excluded
   if (ionode) then
      write(stdout, '(/,8x,a)') 'mode number | included/excluded (+/-)'
      write(stdout, '(8x,a)')   '-------------------------------------'
      do im = 1, nmodes
         if(include_phonon_modes(im)) then
            write(stdout, '(10x,I8,2x,a,12x,a)') im, '|', '+'
         else
            write(stdout, '(10x,I8,2x,a,12x,a)') im, '|', '-'
         end if
      end do
      write(stdout, '(8x,a,/)')   '-------------------------------------'

      ! output phonon mode include into YAML
      write(ymlout, '(/,a)') 'include phonon modes:'
      do im = 1, nmodes
         if(include_phonon_modes(im)) then
            write(ymlout, '(6x,a,a,2x,a)') trim(int_to_char(im)), ':', 'True' 
         else
            write(ymlout, '(6x,a,a,2x,a)') trim(int_to_char(im)), ':', 'False' 
         end if
      end do
   end if
end subroutine

!------
! Functions for reporting memory usage of data structures

! Output details regarding the memory allocation of the kq-pair
! data structure.
subroutine output_kq_pair_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes

   ! kq_pair

   bytes_per_entry = STORAGE_SIZE(kq_pair(1)) / 8
   total_bytes = bytes_per_entry * SIZE(kq_pair, kind=8)

   ! subarrays of kq_pair

   total_allocs = 1
   do i = 1, nk_loc
      if (kq_pair(i)%nkq == 0) cycle

      ! Each element in kq_pair has three subarrays
      total_allocs = total_allocs + 3
      total_bytes = total_bytes + &
         SIZE(kq_pair(i)%iq  , kind=8) * STORAGE_SIZE(kq_pair(i)%iq(1)  ) / 8 + &
         SIZE(kq_pair(i)%ikq , kind=8) * STORAGE_SIZE(kq_pair(i)%ikq(1) ) / 8 + &
         SIZE(kq_pair(i)%nchl, kind=8) * STORAGE_SIZE(kq_pair(i)%nchl(1)) / 8
   enddo

   write(ymlout,'(3x, a)') 'kq_pair [intermediate]:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', nk_loc
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i12/)') 'allocations: ', total_allocs

   call boltzscat_accum_cpu_bytes(total_bytes)
end subroutine output_kq_pair_alloc_stats_to_yaml

! Output details regarding the memory allocation of the scatter
! data structure.
subroutine output_scatter_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes, subelems

   ! scatter

   bytes_per_entry = STORAGE_SIZE(scatter(1)) / 8
   total_bytes = bytes_per_entry * SIZE(scatter, kind=8)

   write(ymlout,'(3x, a)') 'scatter:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_kq_pairs
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i5/)') 'allocations: ', 1

   call boltzscat_accum_cpugpu_bytes(total_bytes)
end subroutine output_scatter_alloc_stats_to_yaml

! Output details regarding the memory allocation of the scatter_channels
! data structure.
subroutine output_scatter_channels_alloc_stats_to_yaml()
   implicit none
   integer :: i, total_allocs
   integer(8) :: bytes_per_entry, total_bytes, subelems

   ! scatter_channels

   bytes_per_entry = STORAGE_SIZE(scatter_channels(1)) / 8
   total_bytes = bytes_per_entry * SIZE(scatter_channels, kind=8)

   write(ymlout,'(3x, a)') 'scatter_channels:'
   write(ymlout,'(6x, a, i5)') 'bytes per entry: ', bytes_per_entry
   write(ymlout,'(6x, a, i12)') 'entries: ', num_scatter_channels
   write(ymlout,'(6x, a, i12)') 'total bytes allocated: ', total_bytes
   write(ymlout,'(6x, a, i5/)') 'allocations: ', 1

   call boltzscat_accum_cpugpu_bytes(total_bytes)
end subroutine output_scatter_channels_alloc_stats_to_yaml

end module boltz_scatter
