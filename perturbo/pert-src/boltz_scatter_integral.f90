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
!  solve iterative linearized Boltzmann equation, and real-time dynamics
!
! Maintenance:
!===============================================================================

#include "pert_config.h"

#undef _OPENACC

module boltz_scatter_integral
   use iso_c_binding
   use kinds, only: i8b
   use pert_const, only: dp, pi
   use boltz_grid, only: grid
   use boltz_utils, only: idx2band
   use pert_utils, only: bose, fermi, gauss
   use pert_param, only: delta_smear
   use pert_output, only: stopwatch
   use boltz_scatter, only: num_kq_pairs, num_scatter_channels, scatter, scatter_channels, qgrid
   use qe_mpi_mod, only: stdout, ionode, mp_barrier, mp_sum, inter_pool_comm
   implicit none
   private

   public :: rates_scat_int, trans_scat_int, cdyna_scat_int
contains

! compute scattering rate
subroutine rates_scat_int(kg, eptemp, efermi, rates)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: eptemp, efermi
   real(dp), intent(out) :: rates(:,:) !rates(kg%numb, kg%numk)
   !local variables
   integer :: i, iq, ik, ikq, im, n, m, bands(3)
   integer(kind=i8b) :: ns
   real(dp) :: wq, enk, emkq, ph_n, fkq, nk2mkq, mkq2nk
#if !defined(STORE_G2DT)
   real(dp) :: g2, dt1, dt2
#else
   real(dp) :: g2_dt1, g2_dt2
#endif

  !call stopwatch('rates_scat_int', 'start')
   rates(:,:) = 0.0_dp

#if defined(_OPENACC)
! OpenACC acceleration
!$acc kernels
#else
! OpenMP acceleration
!$omp parallel do schedule(guided) &

#if defined(CDYN_USE_ARRAY_REDUCE)
!$omp&  reduction(+:rates) &
#endif

#if !defined(STORE_G2DT)
!$omp& default(shared) private(i, iq, ik, ikq, ns, &
!$omp& im, n, m, wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, mkq2nk, nk2mkq, bands)
#else
!$omp& default(shared) private(i, iq, ik, ikq, ns, &
!$omp& im, n, m, wq, enk, emkq, g2_dt1, g2_dt2, ph_n, fkq, mkq2nk, nk2mkq, bands)
#endif

#endif

#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif

         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! eigenvalue
         wq   = qgrid%freq(im, iq)
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)

         ! phonon occupation
         ! the constrain of 'wq > eps_acustic' is aready done in scat_setup.
         ph_n = bose(eptemp, wq) ! N_q\mu

#if !defined(STORE_G2DT)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scatter_channels(ns)%eph_g2

         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq)
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
#else
         g2_dt1 = scatter_channels(ns)%g2_dt1
         g2_dt2 = scatter_channels(ns)%g2_dt2
#endif

         ! compute scattering rates
         fkq = fermi(eptemp, emkq-efermi)
#if !defined(STORE_G2DT)
         mkq2nk = g2 * (dt1*(ph_n + fkq) + dt2*(ph_n + 1.0_dp - fkq))
#else
         mkq2nk = (g2_dt1 * (ph_n + fkq) + g2_dt2 * (ph_n + 1.0_dp - fkq))
#endif

! lock memory when updating this variable
#if defined(_OPENACC)
!$acc atomic
#else
#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
#endif
#endif
         rates(n, ik)  = rates(n, ik)  + mkq2nk

         !update mk+q with contribution from k
         fkq = fermi(eptemp, enk-efermi)
         !delta(Ekq-Ek+wq)=delta(Ek-Ekq-wq); delta(Ekq-Ek-wq)=delta(Ek-Ekq+wq)
#if !defined(STORE_G2DT)
         nk2mkq = g2 * (dt2*(ph_n + fkq) + dt1*(ph_n + 1.0_dp - fkq))
#else
         nk2mkq = (g2_dt2 * (ph_n + fkq) + g2_dt1 * (ph_n + 1.0_dp - fkq))
#endif

#if defined(_OPENACC)
!$acc atomic
#else
#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
#endif
#endif
         rates(m, ikq) = rates(m, ikq) + nk2mkq
      enddo  !; enddo; enddo
#if defined(SCAT_FWD)
   enddo ! i loop
#endif

#if defined(_OPENACC)
!$acc end kernels
#else
!$omp end parallel do
#endif
   ! combine contributions from all the pools/processors.
   call mp_sum(rates, inter_pool_comm)
  !call stopwatch('rates_scat_int', 'stop')
end subroutine rates_scat_int

! iterative process for tranport calculation
subroutine trans_scat_int(kg, eptemp, efermi, mfd, epint)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: efermi, eptemp, mfd(:,:,:) !mfd(3 or 6, kg%numb, kg%nk)
   real(dp), intent(out) :: epint(:,:,:)  !epint(3 or 6, kg%numb, kg%nk)
   !local variables
   integer :: i, iq, ik, ikq, im, n, m, ic, bands(3), ncomp
   integer(kind=i8b) :: ns
   real(dp) :: wq, enk, emkq, ph_n, nk2mkq, mkq2nk, fkq
#if !defined(STORE_G2DT)
   real(dp) :: g2, dt1, dt2
#else
   real(dp) :: g2_dt1, g2_dt2
#endif

  !call stopwatch('trans_scat_int','start')
   !sanity check
   ncomp = size(mfd, 1)
   if( ncomp .ne. size(epint, 1) ) &
      call errore('trans_scat_int',"mismatch dimension in input arguments!",1)

   epint(:,:,:) = 0.0_dp
#if defined(_OPENACC)
! OpenACC acceleration
!$acc kernels
#else
! OpenMP acceleration
!$omp parallel do schedule(guided) &

#if defined(CDYN_USE_ARRAY_REDUCE)
!$omp&  reduction(+:epint) &
#endif

#if !defined(STORE_G2DT)
!$omp& default(shared) private(i, iq, ik, ikq, ns,&
!$omp&  im, n, m, wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, mkq2nk, nk2mkq, ic, bands)
#else
!$omp& default(shared) private(i, iq, ik, ikq, ns,&
!$omp&  im, n, m, wq, enk, emkq, g2_dt1, g2_dt2, ph_n, fkq, mkq2nk, nk2mkq, ic, bands)
#endif

#endif

#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif

         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! eigenvalue
         wq   = qgrid%freq(im, iq)
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)

         ! phonon occupation
         ! the constrain of 'wq > eps_acustic' is aready done in scat_setup.
         ph_n = bose(eptemp, wq) ! N_q\mu

#if !defined(STORE_G2DT)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scatter_channels(ns)%eph_g2

         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq)
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
#else
         g2_dt1 = scatter_channels(ns)%g2_dt1
         g2_dt2 = scatter_channels(ns)%g2_dt2
#endif

            !update nk with contribution from k+q
         fkq = fermi(eptemp, emkq - efermi)
#if !defined(STORE_G2DT)
         mkq2nk = g2  * (dt1*(ph_n + fkq) + dt2*(ph_n + 1.0_dp - fkq))
#else
         mkq2nk = g2_dt1 * (ph_n + fkq) + g2_dt2 * (ph_n + 1.0_dp - fkq)
#endif

#if defined(_OPENACC)
!$acc atomic
#else
#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
#endif
#endif
!           epint(ic, n, ik)  = epint(ic, n, ik)  + mkq2nk

            ! update mk+q with contribution from k
            fkq = fermi(eptemp, enk-efermi);
#if !defined(STORE_G2DT)
            nk2mkq = g2 * (dt2*(ph_n + fkq) + dt1*(ph_n + 1.0_dp - fkq))
#else
            nk2mkq = g2_dt2 * (ph_n + fkq) + g2_dt1 * (ph_n + 1.0_dp - fkq)
#endif
         do ic = 1, ncomp
            !update nk with contribution from k+q
!            fkq = fermi(eptemp, emkq - efermi)
!#if !defined(STORE_G2DT)
!            mkq2nk = g2 * mfd(ic, m, ikq) * (dt1*(ph_n + fkq) + dt2*(ph_n + 1.0_dp - fkq))
!#else
!            mkq2nk = mfd(ic, m, ikq) * (g2_dt1 * (ph_n + fkq) + g2_dt2 * (ph_n + 1.0_dp - fkq))
!#endif

#if defined(_OPENACC)
!$acc atomic
#else
#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
#endif
#endif
            epint(ic, n, ik)  = epint(ic, n, ik)  + mkq2nk * mfd(ic, m, ikq)

            ! update mk+q with contribution from k
!            fkq = fermi(eptemp, enk-efermi);
!#if !defined(STORE_G2DT)
!            nk2mkq = g2 * mfd(ic, n, ik) * (dt2*(ph_n + fkq) + dt1*(ph_n + 1.0_dp - fkq))
!#else
!            nk2mkq = mfd(ic, n, ik) * (g2_dt2 * (ph_n + fkq) + g2_dt1 * (ph_n + 1.0_dp - fkq))
!#endif

! lock memory when updating this variable
#if defined(_OPENACC)
!$acc atomic
#else
#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
#endif
#endif
            epint(ic, m, ikq) = epint(ic, m, ikq) + nk2mkq * mfd(ic, n, ik)
         enddo
      enddo
#if defined(SCAT_FWD)
   enddo
#endif

#if defined(_OPENACC)
!$acc end kernels
#else
!$omp end parallel do
#endif
   !combine contributions from all the pools/processors.
   call mp_sum(epint, inter_pool_comm)
  !call stopwatch('trans_scat_int','stop')
end subroutine trans_scat_int

! Real time dynamics, Refer to eq.17. 18 in Eur. Phys. J. B (2016) 89: 239
! Bernardi, First-principles dynamics of electrons and phonons
subroutine cdyna_scat_int(kg, eptemp, dist, epcol)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: eptemp, dist(kg%numb, kg%nk)
   ! e-ph collision term computed from given distribution function dist(:,:)
   real(dp), intent(out) :: epcol(kg%numb, kg%nk)
   !local variables
   integer :: i, iq, ik, ikq, im, n, m, bands(3)
   integer(kind=i8b) :: ns
   real(dp):: wq, enk, emkq, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq
#if !defined(STORE_G2DT)
   real(dp) :: g2, dt1, dt2
#else
   real(dp) :: g2_dt1, g2_dt2
#endif

   epcol(:,:) = 0.0_dp

#if defined(_OPENACC)
! OpenACC acceleration
!$acc kernels
#else
! OpenMP acceleration
!$omp parallel do schedule(guided) &

#if defined(CDYN_USE_ARRAY_REDUCE)
!$omp&  reduction(+:epcol) &
#endif

!$omp&  default(shared) private(i, iq, ik, ikq, ns, bands, im, n, m, wq, enk, emkq, &

#if !defined(STORE_G2DT)
!$omp&  g2, dt1, dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq)
#else
!$omp&  g2_dt1, g2_dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq)
#endif

#endif
#if defined(SCAT_FWD)
   do i = 1, num_kq_pairs
      iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq

      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = scatter(i)%chl_start_index, scatter(i)%chl_start_index + scatter(i)%nchl - 1
#elif defined(SCAT_REV)
      ! loop over all scatter-channels, using the scat_index to lookup
      ! the corresponding (iq, ik, ikq) values
      do ns = 1, num_scatter_channels
         i = scatter_channels(ns)%scat_index
         iq = scatter(i)%iq;   ik = scatter(i)%ik;    ikq = scatter(i)%ikq
#else
#error At least one of SCAT_FWD or SCAT_REV must be defined.
#endif

         ! map index to (mkq, nk, mu)
         bands = idx2band( scatter_channels(ns)%bands_index, kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)

         ! phonon occupation
         wq   = qgrid%freq(im, iq)
         ph_n = bose(eptemp, wq) ! N_q\mu

         ! electron distributaion
         fnk  = dist(n, ik)
         fmkq = dist(m, ikq)
         ! scattering process:  formula 1
         !fabs = fnk*(one - fmkq)*wgq - fmkq*(one - fnk)*(one + wgq)
         !fem  = fnk*(one - fmkq)*(one + wgq) - fmkq*(one - fnk)*wgq
         ! formula 2
         fabs = fnk*(ph_n + fmkq) - fmkq*(ph_n + 1.0_dp)
         fem  = fnk*(ph_n + 1.0_dp ) - fmkq*(ph_n + fnk)

         ! update nk with contribution from k+q
#if defined(STORE_G2DT)
         ! Use the precomputed -g2*dt1 and -g2*dt2 terms
         mkq2nk = -(scatter_channels(ns)%g2_dt1 * fabs + &
                    scatter_channels(ns)%g2_dt2 * fem)
#else
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scatter_channels(ns)%eph_g2

         ! eigenvalue
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)
         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq)
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)

         mkq2nk = - g2 * (dt1*fabs + dt2*fem)
#endif
         !update nk with contribution from k+q
#if defined(_OPENACC)
!$acc atomic
#else
#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
#endif
#endif
         epcol(n, ik) = epcol(n, ik) + mkq2nk
         ! the inverse: update mk+q with contribution from k
         nk2mkq = - mkq2nk
#if defined(_OPENACC)
!$acc atomic
#else
#if !defined(CDYN_USE_ARRAY_REDUCE)
!$omp atomic
#endif
#endif
         epcol(m, ikq) = epcol(m, ikq) + nk2mkq
      enddo   !; enddo; enddo ! ns loop
#if defined(SCAT_FWD)
   enddo ! i loop
#endif
#if defined(_OPENACC)
!$acc end kernels
#else
!$omp end parallel do
#endif
   ! combine contributions from all the pools/processes.
   call mp_sum(epcol, inter_pool_comm)
end subroutine cdyna_scat_int

end module boltz_scatter_integral
