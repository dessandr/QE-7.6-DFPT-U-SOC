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
!   basic functions used in boltzmann subroutines
!
! Maintenance:
!===============================================================================

!> Various functions used in boltzmann subroutines
module boltz_utils
   use pert_const, only: dp
   private
   public :: num2kpt
   public :: kpt2num
   public :: idx2band
   public :: band2idx
   public :: kpts_plus
   public :: kpts_minus
   public :: inside_win
   public :: velprod
contains

!> Take a non-negative integer; `num`, and the number of k points
!! in 3D; `nk`= (nk1, nk2, nk3). `num2kpt` finds the corresponding
!! k point; `kpt`, based on splitting up `num` as follows
!! `num` = k + j * nk3 + i * nk2 * nk3
!! where \f$ i \in [0, nk1 - 1] \f$, \f$ j \in [0, nk2 - 1] \f$
!! and \f$ k \in [0, nk3 - 1] \f$. 
!!
!! The k point; `kpt`, is given as a fraction [i/nk1, j/nk2, k/nk3]
!! and folded to the Gamma-centered FBZ in crystal coordinates such
!! that each `kpt` value is between -0.5 and 0.5
pure function num2kpt(num, nk) result(kpt)
   implicit none
   integer, intent(in) :: num   !< non-negative integer used to determine `kpt`
   integer, intent(in) :: nk(3) !< number of k points in each direction 
   integer :: ipt(3)            !< used to store [i, j, k]
   integer :: i                 !< used in loop
   real(dp) :: kpt(3)           !< k point in Gamma-centered FBZ crystal coordinates

   ! k = mod(num, nk3) 
   ipt(3) = mod(num, nk(3))
   ! j = mod(num / nk3, nk2)
   ipt(2) = mod(num / nk(3), nk(2))
   ! i = num / (nk2*nk3)
   ipt(1) = num / (nk(2) * nk(3))

   ! kpt as a fraction of total k points
   kpt(:) = real(ipt(:), dp) / real(nk(:), dp)

   ! fold to the Gamma-centered FBZ, in crystal coordinate
   ! if kpt(i) rounds to 1, then subtract nk(i) so it is 
   ! between -nk(i)/2 and 0
   do i = 1, 3
      ipt(i) = merge(ipt(i)-nk(i), ipt(i), nint(kpt(i)) > 0)
   enddo

   ! get the folded kpt as a fraction
   kpt(:) = real(ipt(:), dp) / real(nk(:), dp)

end function num2kpt

! map (kx, ky, kz) to num, if num < 0, then (kx, ky, kz) is not in the list
pure function kpt2num(kpt, nk) result(num)
   implicit none
   real(dp), intent(in) :: kpt(3)
   integer, intent(in) :: nk(3)
   integer :: num
   real(dp), parameter :: eps = 1.0E-5_dp
   ! local variables
   integer :: r(3)  !, i 
   real(dp) :: xkr(3), dis(3)
   ! init num to the default value 
   num = -1
   ! fold to the Gamma-centered FBZ, in crystal coordinate
   ! and check if kpt(i)*nk(i) is a integer or not.
   xkr(:) = (kpt(:) - nint(kpt(:))) * nk(:)
   dis(:) =  xkr(:) - nint(xkr(:))
   ! return -1 if (kx, ky, kz) is not in the k-mesh.
   if( sqrt(dot_product(dis, dis)) > eps ) return
   ! ri = 0...nki-1; r(1)->i; r(2)->j; r(3)->k
   r(:) = mod( nint(xkr(:)+2*nk(:)), nk(:) )
   ! num = k + j*nk3 + i*nk2*nk3 + 1
   num = r(3) + r(2)*nk(3) + r(1)*nk(2)*nk(3)
end function kpt2num

! given index of ikq and ik, compute the index of xq = ikq-ik on q-grid.
pure function kpts_minus(ikq, ik, nk, nq) result(iq)
   implicit none
   integer, intent(in) :: ikq, ik, nk(3), nq(3)
   integer :: iq, q(3) !, i
   real(dp) :: xq(3)
   ! compute difference in i, j, k
   q(1) =     ikq/(nk(2)*nk(3)) -     ik/(nk(2)*nk(3))
   q(2) = mod(ikq/nk(3), nk(2)) - mod(ik/nk(3), nk(2))
   q(3) = mod(ikq, nk(3))       - mod(ik, nk(3))
   !
   !coordinates of q-points
   xq(:) = real(q(:), dp) / real(nk(:), dp)
   !get the index, if xq is not on the q-grid, then iq = -1
   iq = kpt2num(xq, nq)
end function kpts_minus

! given index of ikq and ik, compute the index of xq = ikq+ik on q-grid.
pure function kpts_plus(ikq, ik, nk, nq) result(iq)
   implicit none
   integer, intent(in) :: ikq, ik, nk(3), nq(3)
   integer :: iq, q(3) !, i
   real(dp) :: xq(3)
   ! compute difference in i, j, k
   q(1) =     ikq/(nk(2)*nk(3)) +     ik/(nk(2)*nk(3))
   q(2) = mod(ikq/nk(3), nk(2)) + mod(ik/nk(3), nk(2))
   q(3) = mod(ikq, nk(3))       + mod(ik, nk(3))
   ! 
   !coordinates of q-points
   xq(:) = real(q(:), dp) / real(nk(:), dp)
   !get the index, if xq is not on the q-grid, then iq = -1
   iq = kpt2num(xq, nq)
end function kpts_plus

! map index to (mkq, nk, mu)
!
! Given the number of bands and phonon modes, these values are packed into a
! 32-bit integer as follows:
! * Bits 0..7 = mkq
! * Bits 8..15 = nk
! * Bits 16-31 = mu
pure function idx2band(idx, nbnd) result(band)
#if defined(_OPENACC)
!$acc routine seq
#endif
   implicit none
   integer, intent(in), value :: idx
   integer, intent(in), value :: nbnd
   integer :: band(3)
   ! mu <= nmodes, nk <= nbnd, mkq <= nbnd
   band(1) = ibits(idx,  0,  8) + 1
   band(2) = ibits(idx,  8,  8) + 1
   band(3) = ibits(idx, 16, 15) + 1
end function idx2band

! map (mkq, nk, mu) to index
!
! Given the number of bands and phonon modes, these values are packed into a
! 32-bit integer as follows:
! * Bits 0..7 = mkq
! * Bits 8..15 = nk
! * Bits 16-31 = mu
pure function band2idx(band, nbnd) result(idx)
   implicit none
   integer, intent(in) :: band(3) ! mkq, nk, mu
   integer, intent(in) :: nbnd ! numb, numb, nmodes
   integer :: idx
   idx = 0
   call mvbits(band(1) - 1, 0,  8, idx,  0)
   call mvbits(band(2) - 1, 0,  8, idx,  8)
   call mvbits(band(3) - 1, 0, 15, idx, 16)
end function band2idx

! check if xkg has energy levels inside [emin, emax]
! only bands in [bmin, bmax] are considered.
function inside_win(el, xkg, emin, emax, bmin, bmax)
   use band_structure, only: electron_wann, solve_eigenvalue_vector
   implicit none
   type(electron_wann), intent(in) :: el
   integer,  intent(in) :: bmin, bmax
   real(dp), intent(in) :: emin, emax, xkg(3)
   integer :: inside_win
   ! local variables
   real(dp) ::  eval(el%nb)
   integer :: lower, upper, ib 

   call solve_eigenvalue_vector(el, xkg, eval)
   ! only bands in [bmin, bmax] are considered.
   lower = bmin - 1
   upper = bmax + 1
   do ib = bmin, bmax
      if( eval(ib) < emin ) lower = lower + 1
      if( eval(ib) > emax ) upper = upper - 1
   enddo
   ! no level inside [emin, emax]
   ! e.g. all levels are lower than emin .or larger than emax
   ! or have levels lower than emin and levels larger than emax.
   if( upper .eq. (lower+1) ) then
      inside_win = lower
   ! have levels inside [emin, emax]
   elseif(upper > (lower+1) ) then
      inside_win = -2
   else
      ! error appears, not a logical result
      inside_win = -4
   endif
   return
end function inside_win

pure function velprod(v1, v2) result(vv)
   implicit none
   real(dp), intent(in) :: v1(3)
   real(dp), intent(in) :: v2(3)
   !
   real(dp) :: vv(6)

   vv(1) = v1(1) * v2(1)  ! XX: v_x * v_x
   vv(2) = v1(1) * v2(2)  ! XY: v_x * v_y
   vv(3) = v1(2) * v2(2)  ! YY: v_y * v_y
   vv(4) = v1(1) * v2(3)  ! XZ: v_x * v_z
   vv(5) = v1(2) * v2(3)  ! YZ: v_y * v_z
   vv(6) = v1(3) * v2(3)  ! ZZ: v_z * v_z
end function velprod

end module boltz_utils
