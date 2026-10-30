!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jspark <jpark7@caltech.edu>
! Comment:
!
! Maintenance:
!===============================================================================
subroutine calc_average_spin_lifetime()
   use pert_const, only: dp, ryd2ev
   use pert_output,only: progressbar_init, progressbar, load_imsigma_spin
   use pert_data,  only: volume, epr_fid, kc_dim
   use boltz_grid, only: grid, init_boltz_grid, boltz_grid_load
   use qe_mpi_mod, only: ionode, stdout, mp_barrier, inter_pool_comm
   use pert_param, only: ntemper, temper, efermi, band_min, band_max, boltz_de,& 
                          boltz_emax, boltz_emin, spinor, boltz_kdim, prefix
   !use boltz_trans_mod, only: trans_rateavg_calc
   use boltz_trans_output, only: output_spin_lifetime
   use HDF5_utils
   !
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
 
   implicit none
   type(grid) :: kg
   logical :: read_rate
   integer :: nestep, ntmpr, i, it, ik, ib, iv
   real(dp) :: emin, emax, ef, enk, mevhbar
   real(dp), allocatable :: ene(:), rates(:,:,:,:), rateavg(:,:), imsgm(:,:,:), dor(:,:,:)
   !
   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec

   mevhbar = 658.212E-6_dp !hbar in meV * ns
   if(ionode) write(stdout,'(1x,a)') 'Start calc_average_spin_lifetime...'

   if(ntemper .eq. 0) call errore('calc_average_spin_lifetime', &
      'ftemper is not specified, missing temperatures and chemical potential',1)  

   call init_electron_wann(epr_fid, kc_dim, elec)
   call init_boltz_grid(kg, elec, band_min, band_max)

   ntmpr = size(temper)
   call boltz_grid_load(kg, 'k')
   !call init_boltz_grid(kg, elec, band_min, band_max)
   
   ! setup up energy windows for calc_average_spin_lifetime calculations
   emin = max(boltz_emin, minval(kg%enk)) !merge(boltz_emin, max(boltz_emin, minval(kg%enk)), .false.)
   emax = min(boltz_emax, maxval(kg%enk))  !)merge(boltz_emax, min(boltz_emax, maxval(kg%enk)), .false.)
   if(ionode) write(stdout,'(2x, a, 2(1x,f12.6))') &
      'Energy window (eV):', emin*ryd2ev, emax*ryd2ev
   if(emin > emax) call errore('calc_average_spin_lifetime','illegal energy window',1)
   !
   call mp_barrier(inter_pool_comm)
   write(stdout,'(5x,a)') '>finish init grid.'
 
   !setup energy grid
   nestep = int( floor((emax -emin)/boltz_de) ) + 3
   allocate(ene(nestep))
   do i = 1, nestep
      ene(i) = emin + (i-2)*boltz_de
   enddo
   
   allocate(rates(kg%numb, kg%nk, ntmpr,4), imsgm(kg%numb, kg%nk_irr, ntmpr))
   !read imsigma files: imaginary part of the selfenergy in meV.
   do iv = 1, 1
      call load_imsigma_spin(imsgm, read_rate, iv)  
      !call load_imsigma_tot(imsgm, read_rate, iv)  
      if(read_rate) then
         do it = 1, ntmpr
         do ik = 1, kg%nk
         do ib = 1, kg%numb
            !Scattering rates in (ns^-1) unit, (Rate = 2/hbar * ImSigma)
            rates(ib,ik,it,iv) = 2.0_dp/mevhbar*imsgm(ib, kg%kpt2ir(ik), it) !for spin
         enddo; enddo; enddo
      elseif (.not. read_rate) then
         call errore('calc_average_spin_lifetime',"Can not find "//trim(prefix)//".imsigma_flip",1)
      endif
   enddo
   deallocate( imsgm )
   
   write(stdout,'(5x,a)') '>finish init rates'

   !! now all data we needed are ready. start to do the actual work.
   ! allocate work space first.
   allocate( rateavg(ntmpr,4), dor(nestep,ntmpr,4) )
   
   rateavg = 0.0_dp
   if(ionode) call progressbar_init('calc_average_spin_lifetime:')
   
   do it = 1, ntmpr
      do iv = 1, 1
      !compute thermal average of scattering rates
      
      !call trans_rateavg_calc(kg, ene, temper(it), efermi(it), spinor, &
      !  volume, rates(:,:,it,iv), rateavg(it,iv), dor(:,it,iv))
      !rateavg(it,iv) = 1.0_dp/(sum(1.0_dp/rates(:,:,it,iv)*exp(-(kg%enk(:,:)-efermi(it))/temper(it)))/sum(exp(-(kg%enk(:,:)-efermi(it))/temper(it))))
      rateavg(it,iv) = (sum(rates(:,:,it,iv)*exp(-(kg%enk(:,:)-efermi(it))/temper(it)))/sum(exp(-(kg%enk(:,:)-efermi(it))/temper(it))))
      !output progress
      enddo
      if(ionode) call progressbar(it, ntmpr)
   enddo
   
   if(ionode) write(stdout,'(5x, a)')   '>finish calc_average_spin_lifetime '
   
   !output spin lifetime
   call output_spin_lifetime(temper, efermi, rateavg)
   !kcall output_density_of_rates(ene, dor, temper, efermi)
   deallocate(ene, rates, rateavg, dor)
end subroutine calc_average_spin_lifetime
