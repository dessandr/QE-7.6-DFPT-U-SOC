!===============================================================================
! Author: jspark <jinsoop412@gmail.com>
! 
! File: calc_spin_wann.f90
! Create Date: 2019-10-14 17:16:41
!
! Comment:
!
!===============================================================================

subroutine calc_spin_wann(el, kdim, nwan, at, wan_center, nband, numk, xk, rot_wan)
   use kinds, only: dp
   use constants, only: tpi
   use qe_mpi_mod, only: mp_sum, inter_pool_comm, mp_split_pools, mp_bcast, mp_sum
   use electron_wannier, only: electron_wann, set_ws_cell_el_spin
   use io_global, only : ionode, ionode_id, stdout
   use io_files, only: prefix
   use input_param, only: spn_formatted
   implicit none
   type(electron_wann), intent(inout) :: el
   ! number of k-point and bands (may differ from number of wannier functons)
   integer, intent(in) :: kdim(3), nwan, nband, numk
   ! wannier center in crystal coordinate
   real(dp), intent(in) :: at(3,3), wan_center(3, nwan)
   real(dp), intent(in) :: xk(3,numk) ! k-points in crystal coordinate
   complex(dp), intent(in) :: rot_wan(nband, nwan, numk) ! rot_wan(nband, el%nb, numk)
   !
   integer :: m, ib, jb, ik, n, ir, ire, nk_tmp, ik_st, ik_end, nk_loc, iib, jjb, ios
   real(dp) :: fac, rdotk, rvec(3)
   complex(dp) :: cfac, spin_k_1, spin_k_2, spin_k_3
   character(len=256) :: fname

   !variables imported from wannier90
   complex(kind=dp), allocatable :: spn_o(:, :, :, :), spn_temp(:, :)
   real(kind=dp)                 :: s_real, s_img
   integer                       :: spn_in, nb_tmp, nkp_tmp, ierr, s, counter
   character(len=60)             :: header

   integer, external :: find_free_unit

   nk_tmp = kdim(1) * kdim(2) * kdim(3)
   !sanity check
   if(nk_tmp .ne. numk) &
      call errore('calc_spin_wann','inconsistent kdim and numk',1)
   ! setup wigner seitz cell
   call set_ws_cell_el_spin(el, kdim, nwan, at, wan_center)

   ! 1. Read from .spn file the original spin matrices <psi_nk|sigma_i|psi_mk>
   !    (sigma_i = Pauli matrix) between ab initio eigenstates 
   if (ionode) then
      allocate (spn_o(nband, nband, numk, 3))

      fname = trim(prefix)//'.spn'
      spn_in = find_free_unit()
      
      if (spn_formatted) then
         open (unit=spn_in, file=trim(prefix)//'.spn', form='formatted', &
               status='old', err=100, iostat=ios)
         !write (stdout, '(/a)', advance='no') &
         !  ' Reading spin matrices from '//trim(prefix)//'.spn : '
         read (spn_in, *, err=110, end=110) header
         !write (stdout, '(a)') trim(header)
         read (spn_in, *, err=110, end=110) nb_tmp, nkp_tmp
      else
         open (unit=spn_in, file=trim(prefix)//'.spn', form='unformatted', &
               status='old', err=100, iostat=ios)
         !write (stdout, '(/a)', advance='no') &
         !  ' Reading spin matrices from '//trim(prefix)//'.spn : '
         read (spn_in, err=110, end=110) header
         !write (stdout, '(a)') trim(header)
         read (spn_in, err=110, end=110) nb_tmp, nkp_tmp
      endif
      if (nb_tmp .ne. nband) &
          call errore('calc_spin_wann', trim(prefix)//'.spn has wrong number of bands', 1)
      if (nkp_tmp .ne. numk) &
          call errore('calc_spin_wann', trim(prefix)//'.spn has wrong number of k-points', 1)
      if (spn_formatted) then
         do ik = 1, numk
            do m = 1, nband
               do n = 1, m
                  read (spn_in, *, err=110, end=110) s_real, s_img
                  spn_o(n, m, ik, 1) = cmplx(s_real, s_img, dp)
                  read (spn_in, *, err=110, end=110) s_real, s_img
                  spn_o(n, m, ik, 2) = cmplx(s_real, s_img, dp)
                  read (spn_in, *, err=110, end=110) s_real, s_img
                  spn_o(n, m, ik, 3) = cmplx(s_real, s_img, dp)
                  ! Read upper-triangular part, now build the rest
                  spn_o(m, n, ik, 1) = conjg(spn_o(n, m, ik, 1))
                  spn_o(m, n, ik, 2) = conjg(spn_o(n, m, ik, 2))
                  spn_o(m, n, ik, 3) = conjg(spn_o(n, m, ik, 3))
               end do
            end do
         enddo
      else
         allocate (spn_temp(3, (nband*(nband + 1))/2), stat=ierr)
         if (ierr /= 0) call errore('calc_spin_wann','Error in allocating spn_temp',1)
         do ik = 1, numk
            read (spn_in) ((spn_temp(s, m), s=1, 3), m=1, (nband*(nband + 1))/2)
            counter = 0
            do m = 1, nband
               do n = 1, m
                  counter = counter + 1
                  spn_o(n, m, ik, 1) = spn_temp(1, counter)
                  spn_o(m, n, ik, 1) = conjg(spn_temp(1, counter))
                  spn_o(n, m, ik, 2) = spn_temp(2, counter)
                  spn_o(m, n, ik, 2) = conjg(spn_temp(2, counter))
                  spn_o(n, m, ik, 3) = spn_temp(3, counter)
                  spn_o(m, n, ik, 3) = conjg(spn_temp(3, counter))
               end do
            end do
         end do
         deallocate (spn_temp, stat=ierr)
         if (ierr /= 0) call errore('calc_spin_wann', 'Error in deallocating spm_temp',1)
      endif
      close (spn_in)
   endif
   ! broacdast : spinmatrix in Bloch basis, <psink|sigma_i|psi_mk>
   if(.not. allocated(spn_o) )  allocate (spn_o(nband, nband, numk, 3))
   call mp_bcast(spn_o, ionode_id, inter_pool_comm)
   
   ! for mpi parallel
   call mp_split_pools(numk, ik_st, ik_end, nk_loc)
   

   ! 2. Go to Wannier represenatation and do Fourier transform
   !
   m = 0
   do jb = 1, nwan
   do ib = 1, jb
      !m = ( jb * (jb-1) ) / 2 + ib
      m = m + 1
      !init
      el%spin_r(m)%hop(:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

      do ik = ik_st, ik_end
         !Preparing S(k)_{ij} = [U^{\dagger} S U]_{ij}
         spin_k_1 = cmplx(0.0_dp, 0.0_dp, kind=dp)
         spin_k_2 = cmplx(0.0_dp, 0.0_dp, kind=dp)
         spin_k_3 = cmplx(0.0_dp, 0.0_dp, kind=dp)
         !
         do iib = 1, nband
            do jjb = 1, nband
               spin_k_1 = spin_k_1 + conjg(rot_wan(jjb,ib,ik)) * spn_o(jjb,iib,ik,1) * rot_wan(iib,jb,ik)
               spin_k_2 = spin_k_2 + conjg(rot_wan(jjb,ib,ik)) * spn_o(jjb,iib,ik,2) * rot_wan(iib,jb,ik)
               spin_k_3 = spin_k_3 + conjg(rot_wan(jjb,ib,ik)) * spn_o(jjb,iib,ik,3) * rot_wan(iib,jb,ik)
            enddo
         enddo
         !N.B. S(k)_{ji} = conjg( S(k)_{ij} ) since S(k) is Hermitian

!$omp parallel do schedule(guided) default(shared) private(ir, ire, rvec, rdotk, cfac)
         do ir = 1, el%spin_r(m)%ws_el%nr
            ire = el%spin_r(m)%ws_el%rvec(ir)
            rvec(:) = el%rvec_set(:, ire)
            !
            rdotk = tpi * dot_product(xk(:,ik), rvec)
            ! exp(-i k*R) / N_k
            cfac = cmplx(cos(rdotk), -sin(rdotk), kind=dp)
            ! spin_r(ir,1) = spin_r(ir,1) + spin_k_1 * factor
            el%spin_r(m)%hop(ir,1) = el%spin_r(m)%hop(ir,1) + spin_k_1 * cfac
            el%spin_r(m)%hop(ir,2) = el%spin_r(m)%hop(ir,2) + spin_k_2 * cfac
            el%spin_r(m)%hop(ir,3) = el%spin_r(m)%hop(ir,3) + spin_k_3 * cfac
         enddo
!$omp end parallel do

      enddo
   enddo; enddo

   deallocate(spn_o)
   
   !collect results from different pools
   do m = 1, (nwan*(nwan+1))/2
      call mp_sum(el%spin_r(m)%hop(:,:), inter_pool_comm)
      ! include the ndeg
      do ir = 1, el%spin_r(m)%ws_el%nr
         fac = 1.0_dp / real( el%spin_r(m)%ws_el%ndeg(ir) * numk, dp)
         el%spin_r(m)%hop(ir,1) = el%spin_r(m)%hop(ir,1) * fac
         el%spin_r(m)%hop(ir,2) = el%spin_r(m)%hop(ir,2) * fac
         el%spin_r(m)%hop(ir,3) = el%spin_r(m)%hop(ir,3) * fac
      enddo
   enddo
   return
100 call errore('calc_spin_wann', 'Problem opening '//trim(fname), abs(ios))   
110 call errore('calc_spin_wann', 'Problem reading '//trim(fname), abs(ios))   
end subroutine calc_spin_wann
