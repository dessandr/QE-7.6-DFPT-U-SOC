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

module band_structure
   use pert_const, only: dp, ci
   use wigner_seitz_cell, only: ws_cell
   use electron_wannier, only: electron_wann, set_ws_cell_el, set_ws_cell_el_spin
   use qe_mpi_mod, only: ionode, mp_bcast, inter_pool_comm, ionode_id, ionode, stdout
   implicit none
   public
   
   type(electron_wann), save :: elec

   public :: init_electron_wann, solve_eigenvalue_vector, solve_band_velocity, &
             init_spin_wann, diagonalize_spinmatrix
contains

subroutine init_electron_wann(file_id, kdim, el)
   use hdf5_utils
   use pert_data, only: num_wann, at, wannier_center_cryst
   use epr_hdf5_io, only: read_electron_wannier
   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: kdim(3)
   type(electron_wann), intent(out) :: el
   !
   integer :: i, nelem

   ! wannier_center_cryst(:, i), the center of the i-th wannier function in crystal coordinate
   call set_ws_cell_el(el, kdim, num_wann, at, wannier_center_cryst)
   if(ionode) call read_electron_wannier(file_id, el)
   !bcast
   nelem = el%nb * (el%nb+1) / 2
   do i = 1, nelem
      call mp_bcast(el%ham_r(i)%hop(:), ionode_id, inter_pool_comm)
   enddo
   
   ! check if rset_cart is initialized
   if( .not. associated(el%rset_cart) ) then
      allocate( el%rset_cart(3, el%nrvec) )
      el%rset_cart(:,:) = el%rvec_set(:,:)
      ! compute R in cartesian coordinate, needed in band velocity calculations
      call cryst_to_cart(el%nrvec, el%rset_cart, at, 1)
   endif
end subroutine init_electron_wann
   

subroutine solve_eigenvalue_vector(el, kpt, eig, evec)
   use pert_utils, only: hmat_upper_diag, get_exp_ikr
   implicit none
   type(electron_wann), intent(in) :: el
   real(dp), intent(in) :: kpt(3) ! in crysetal coordinate
   real(dp), intent(out) :: eig(el%nb)
   complex(dp), intent(out), optional :: evec(el%nb, el%nb)
   !local
   integer :: m, ir, ire, nelem
   complex(dp), allocatable :: hamk_upper(:), exp_ikr(:)
   type(ws_cell), pointer :: ws
   
   nelem = el%nb * (el%nb+1) / 2
   allocate( hamk_upper(nelem), exp_ikr(el%nrvec) )
   hamk_upper = cmplx(0.0_dp, 0.0_dp, kind=dp)
   exp_ikr = cmplx(0.0_dp, 0.0_dp, kind=dp)

   call get_exp_ikr(kpt, el%rvec_set, exp_ikr)
   ! prepare hamk_upper(k)
   do m = 1, nelem
      ws => el%ham_r(m)%ws_el

      do ir = 1, ws%nr
         ire = ws%rvec(ir)
         ! the 1/ws%ndeg(ir) is already included in hop(:) ! / real(ws%ndeg(ir), dp)
         hamk_upper(m) = hamk_upper(m) + exp_ikr(ire) * el%ham_r(m)%hop(ir) 
      enddo
   enddo
   
   !diagonialzed 
   if( present(evec) ) then
      call hmat_upper_diag(hamk_upper, el%nb, eig, evec)
   else
      call hmat_upper_diag(hamk_upper, el%nb, eig)
   endif

   deallocate( hamk_upper, exp_ikr)
end subroutine solve_eigenvalue_vector


subroutine solve_band_velocity(el, kpt, velocity, eigvalue, eigvector)
   use pert_utils, only: hmat_upper_diag, get_exp_ikr
   implicit none
   type(electron_wann), intent(in) :: el
   real(dp), intent(in) :: kpt(3) ! in crysetal coordinate
   real(dp), intent(out):: velocity(3, el%nb)
   real(dp), intent(out), optional :: eigvalue(el%nb)
   complex(dp), intent(out), optional  :: eigvector(el%nb, el%nb)
   ! local array
   integer :: i, j, m, ir, ire, ia, nelem
   real(dp) :: del_eig(el%nb), ene(el%nb)
   complex(dp) :: evec(el%nb, el%nb), hamk(el%nb, el%nb, 3)
   complex(dp), allocatable :: hamk_upper(:), exp_ikr(:)
   type(ws_cell), pointer :: ws
   
   nelem = el%nb * (el%nb+1) / 2
   allocate( hamk_upper(nelem), exp_ikr(el%nrvec) )
   hamk_upper = cmplx(0.0_dp, 0.0_dp, kind=dp)
   exp_ikr = cmplx(0.0_dp, 0.0_dp, kind=dp)
   
   call get_exp_ikr(kpt, el%rvec_set, exp_ikr)
   ! prepare hamk_upper(k)
   do m = 1, nelem
      ws => el%ham_r(m)%ws_el

      do ir = 1, ws%nr
         ire = ws%rvec(ir)
         ! the 1/ws%ndeg(ir) is already included in hop(:) ! / real(ws%ndeg(ir), dp)
         hamk_upper(m) = hamk_upper(m) + exp_ikr(ire) * el%ham_r(m)%hop(ir) 
      enddo
   enddo
   call hmat_upper_diag(hamk_upper, el%nb, ene, evec)
   !
   if( present(eigvalue)  ) eigvalue(:) = ene(:)
   if( present(eigvector) ) eigvector(:,:) = evec(:,:)

   !compute band velocity dH/dk_i: delHH(:,:,i), i=1, 2, 3, or x, y, z
   !corresponding to 'fourier_R_to_k(kpt,HH_R,delHH(:,:,1),1)' in postw90
   ! H(k)^a = dH(k)/dk_a  = sum_R 1/ndegen(R) e^{ikR} H_k(R) *{i*R_a}
   ! NOTE: since H(k)_{ij} = conjg(H(k)_ji), H(k)^a_{ij} = conjg(H(k)^{a}_{ji})
   !
   hamk = cmplx(0.0_dp, 0.0_dp, kind=dp)
   m = 0
   do j = 1, el%nb
   do i = 1, j
      m = m + 1
      ws => el%ham_r(m)%ws_el
      
      do ir = 1, ws%nr
         ire = ws%rvec(ir)
         !e^{ikR} * {i*R_a}, R_a in cart coordinate, in the unit of alat. 
         ! corresponding to k in the unit of 1/alat.
         ! now hamk corresponding to delHH(:,:,ia) in postw90.
         ! sum_R 1/ndegen(R) e^{ikR} H_k(R) *{i*R_a}, 1/ndegen(R) is included in hop(ir)
         hamk(i,j,:) = hamk(i,j,:) + &
            (el%ham_r(m)%hop(ir) * exp_ikr(ire) * ci) * el%rset_cart(:,ire)
      enddo
      
      ! H(k)^a_{ji} = conjg( H(k)^{a}_{ij} )
      if(i .ne. j) then
         hamk(j,i,:) = conjg( hamk(i,j,:) )
      else
         !hamk is real when i .eq. j
         hamk(j,i,:) = cmplx(real(hamk(i,j,:)), 0.0_dp, kind=dp)
      endif
   enddo; enddo
   
   ! similar to get_deleig_a(del_eig(:,i),eig,delHH(:,:,i),UU) in postw90
   do ia = 1, 3
      call calc_deleig_a(hamk(:,:,ia), el%nb, ene, evec, del_eig)
      ! the unit of velocity here is Ryd*Bohr*alat. 
      ! checked by compare with postw90.
      velocity(ia,:) = del_eig
   enddo

   deallocate( hamk_upper, exp_ikr )
end subroutine solve_band_velocity


  !==========================!
  ! Band derivatives dE/dk_a , adapted from postw90
  !==========================!
subroutine calc_deleig_a(delHH_a, nbnd, eig, UU, deleig_a)
  use pert_utils, only: hmat_diag
  use pert_const, only: dp, ci, czero, ryd2mev, cone
  !use w90_parameters, only : num_wann,use_degen_pert,degen_thr
  implicit none 
  logical, parameter :: use_degen_pert = .true.  ! turn on by default
  real(dp), parameter :: degen_thr = 0.1/ryd2mev ! 0.1 mev

  integer, intent(in) :: nbnd
  real(dp),    intent(in) :: eig(nbnd)
  complex(dp), intent(in) :: delHH_a(nbnd, nbnd)
  complex(dp), intent(in) :: UU(nbnd, nbnd)
  real(dp),    intent(out):: deleig_a(nbnd)

  ! Misc/Dummy
  integer :: i, degen_min, degen_max, ndim
  real(kind=dp) :: diff
  complex(dp) :: delHH_bar_a(nbnd,nbnd), U_deg(nbnd,nbnd)
  
  if(use_degen_pert) then
     !delHH_bar_a=utility_rotate(delHH_a,UU,num_wann)
     ! (UU)^dagger.mat.UU, where 'UU' is a unitary matrix     
     ! delHH_bar_a=matmul(matmul(transpose(dconjg(UU)), delHH_a), UU)
     call zgemm ('c', 'n', nbnd, nbnd, nbnd, cone, UU, &
         nbnd, delHH_a, nbnd, czero, U_deg, nbnd)
     call zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, U_deg, &
         nbnd, UU, nbnd, czero, delHH_bar_a, nbnd)

     ! Assuming that the energy eigenvalues are stored in eig(:) in
     ! increasing order (diff >= 0)
     i=0
     do 
        i=i+1
        if(i>nbnd) exit
        if(i+1 <= nbnd) then
           diff=eig(i+1)-eig(i)
        else
           ! i-th is the highest band, and it is non-degenerate
           diff =degen_thr+1.0_dp
        end if
        if(diff<degen_thr) then
           ! Bands i and i+1 are degenerate 
           degen_min=i
           degen_max=degen_min+1
           ! See if any higher bands are in the same degenerate group
           do
              if(degen_max+1>nbnd) exit
              diff=eig(degen_max+1)-eig(degen_max)
              if(diff<degen_thr) then
                 degen_max=degen_max+1
              else
                 exit
              end if
           end do
           ! Bands from degen_min to degen_max are degenerate. Diagonalize 
           ! the submatrix in Eq.(31) YWVS07 over this degenerate subspace.
           ! The eigenvalues are the band gradients
           !
           ndim=degen_max-degen_min+1
           call hmat_diag(delHH_bar_a(degen_min:degen_max,degen_min:degen_max),&
              ndim, deleig_a(degen_min:degen_max), U_deg(1:ndim,1:ndim))
           ! Scanned bands up to degen_max
           i=degen_max
        else
           ! Use non-degenerate form [Eq.(27) YWVS07] for current (i-th) band
           deleig_a(i) = aimag(ci*delHH_bar_a(i,i))
        end if
     end do
  else
     ! Use non-degenerate form for all bands
     !deleig_a(:)=aimag(ci*utility_rotate_diag(delHH_a(:,:),UU,nbnd))
     call zgemm ('c', 'n', nbnd, nbnd, nbnd, cone, UU, &
         nbnd, delHH_a, nbnd, czero, U_deg, nbnd)
     call zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, U_deg, &
         nbnd, UU, nbnd, czero, delHH_bar_a, nbnd)
     do i = 1, nbnd
        deleig_a(i) = aimag(ci*delHH_bar_a(i,i))
     enddo
  end if
end subroutine calc_deleig_a


subroutine init_spin_wann(file_id, kdim, el)
   use hdf5_utils
   use pert_data, only: num_wann, at, wannier_center_cryst
   use epr_hdf5_io, only: read_spin_wannier
   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: kdim(3)
   type(electron_wann), intent(inout) :: el
   !
   integer :: i, nelem

   ! wannier_center_cryst(:, i), the center of the i-th wannier function in crystal coordinate
   call set_ws_cell_el_spin(el, kdim, num_wann, at, wannier_center_cryst)
   if(ionode) call read_spin_wannier(file_id, el, 0, 0, 1)
   !bcast
   nelem = el%nb * (el%nb+1) / 2
   do i = 1, nelem
      call mp_bcast(el%spin_r(i)%hop(:,:), ionode_id, inter_pool_comm)
   enddo
end subroutine init_spin_wann

subroutine diagonalize_spinmatrix(el, kpt, eig, evec, spin_diag)
! 1. Diagonalize the spin matrix along the each principal axis and update evec
! 2. Stores the spin expectation value of each bands to spin_diag
   use pert_const, only: czero, cone, ryd2mev
   use pert_data, only: linvsym
   use pert_utils, only: get_exp_ikr, hmat_diag
   implicit none
   type(electron_wann), intent(in) :: el
   real(dp), intent(in) :: kpt(3) ! in crystal coordinate
   real(dp), intent(in) :: eig(el%nb)
   complex(dp), intent(inout) :: evec(el%nb, el%nb)
   !integer :: quantization_axis ! 1:x 2:y 3:z 4:alpha
   real(dp), intent(out) :: spin_diag(el%nb) ! spin expectation value
   complex(dp) :: evec_spin(el%nb, el%nb)
   !blas/zdotu: product of two complex vector: sum(a(:) * b(:))
   complex(dp), external :: zdotu
   integer :: i, j, m, nelem, ir, ire, ndeg, is
   complex(dp) :: spink(el%nb, el%nb)
   ! temporary matrices to store 2x2 spin matrices and gauge transformation matrices
   complex(dp) :: spink_sub(2,2), evec_spin_sub(2,2)
   ! temporary matrices to store 4x4 spin matrices and gauge transformation matrices
   complex(dp) :: spink_sub4(4,4), evec_spin_sub4(4,4)

   logical :: lbisect

   ! temporary variables
   complex(dp) :: spinktmp(el%nb, el%nb)
   complex(dp) :: evec_new(el%nb, el%nb)
   real(dp), parameter :: degen_thr = 0.01_dp/ryd2mev ! 0.01 mev

   complex(dp), allocatable :: exp_ikr(:)
   type(ws_cell), pointer :: ws
 
   nelem = el%nb * (el%nb+1)/2
   allocate(exp_ikr(el%nrvec) )
   exp_ikr = cmplx(0.0_dp, 0.0_dp, kind=dp)
   spink = cmplx(0.0_dp, 0.0_dp, kind=dp)

   ! Inverse fourier transform the spin matrix to fine k mesh (Real Wann -> Recip Wann)
   ! S~_k'   = sum_R 1/n_degen(R) e^{-ik'R    } S_k(R)
   call get_exp_ikr(kpt, el%rvec_set, exp_ikr)
   m = 0
   do j = 1, el%nb
   do i = 1, j
      m = m + 1
      ws => el%spin_r(m)%ws_el

      do ir = 1, ws%nr
         ire = ws%rvec(ir)
         ! the 1/ws%ndeg(ir) is already included in hop(:) ! / real(ws%ndeg(ir), dp)
         spink(i,j) = spink(i,j) + exp_ikr(ire) * el%spin_r(m)%hop(ir,4)
      enddo

      ! S(k)^a_{ji} = conjg( S(k)^{a}_{ij} )
      if(i .ne. j) then
         spink(j,i) = conjg( spink(i,j) )
      else
         !spink is real when i .eq. j
         spink(j,i) = cmplx(real(spink(i,j)), 0.0_dp, kind=dp)
      endif
   enddo
   enddo
   
   ! Un-rotate the spin matrix spink to Bloch space, fine k mesh (Recip Wann -> Recip Bloch)
   ! spink  = [evec^\dagger * spink] * evec
   call zgemm ('c', 'n', el%nb, el%nb, el%nb, cone, evec, &
               el%nb, spink, el%nb, czero, spinktmp, el%nb)
   call zgemm ('n', 'n', el%nb, el%nb, el%nb, cone, spinktmp, &
               el%nb, evec, el%nb, czero, spink, el%nb)

   if(.true.) then
      ! Diagonalize spin matrix in all Kramers degenerate space
      ! The degeneracy is assumed to be 2, except for the high symmetry points bisecting the b's
      ! evec_spin is the eivengector that diagonalizes the Kramer's subspaces in the spin matrix
      evec_spin(:,:) = czero
      i = 1 ! starting band index of the current degenerate subspace
      ndeg = 2 ! degeneracy of the current degenerate subspace
      ! check if this kpt is a high-symmetry point which bisects the crystal b vector
      lbisect = all(abs(2.d0*kpt-nint(2.d0*(kpt)))/2.d0.le.1.0E-12_dp)
      do while (i .le. el%nb)
         j = i+1 ! ending band index of the current degenerate subspace
         ! if the current k-point bisects a reciprocal lattice vector then check the degeneracy of
         ! band i
         if(lbisect) then
           !do while((j+1).le.el%nb)
           if( (j+2).le.el%nb) then
           if(abs(eig(j+1)-eig(i)).le.degen_thr) then
              ! increase the degeneacy of the current subspace by 2 : degeneracy of 6 is unlikely
              j = j+2
           endif
           endif
           !enddo
         endif
         ndeg = j-i+1
         ! diagonalize spink(i:j,i:j) in this degenerace subspace
         if(ndeg .eq. 2) then ! to suppress runtime warning, pass the copy of the matrix
            spink_sub = spink(i:j,i:j)
            evec_spin_sub = evec_spin(i:j,i:j)
            call hmat_diag(spink_sub, ndeg, spin_diag(i:j), evec_spin_sub)
            evec_spin(i:j,i:j) = evec_spin_sub
         else if(ndeg.eq.4) then
            spink_sub4 = spink(i:j,i:j)
            evec_spin_sub4 = evec_spin(i:j,i:j)
            call hmat_diag(spink_sub4, ndeg, spin_diag(i:j), evec_spin_sub4)
            evec_spin(i:j,i:j) = evec_spin_sub4
         else
            call hmat_diag(spink(i:j,i:j), ndeg, spin_diag(i:j), evec_spin(i:j,i:j))
         endif
         i = j+1
      enddo
      ! Rotate evec with evec_spin
      ! U_k_new (k') = U_k (k') * U_k_spin (k')
      call zgemm ('n','n',el%nb,el%nb,el%nb,cone,evec,el%nb,evec_spin,el%nb,czero,evec_new,el%nb)
      evec=evec_new
   else
      !set spin_diag to the diagonal part of spink
      forall (i=1:el%nb) spin_diag(i)=real(spink(i,i))
   endif
   deallocate(exp_ikr)
end subroutine diagonalize_spinmatrix

end module band_structure
