!===============================================================================
! Author: jyyao <jyyao27@gmail.com>
! 
! File: tdep_brents_method.f90
! Create Date: 2024-12-05
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Comment:
!      Brent's method for finding the zero of a function
!      Adapted from TDEP/src/libolle/lo_brents_method.f90
!  Refer to the following literature:
!      F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)
!
! Reverse communication interface to Brent's method of finding the zero of a function.
!
! ! pick x0 and x1 so that they safely bracket 0
! y0=function(x0)
! y1=function(x1)
! ! initialize the solver
! call bt%init( x0,x1,y0,y1,tol )
! ! iterate to find zero
! do iter=1,maxiter
!     ! check for convergence
!     if ( abs(bt%current_y()) .lt. bt%tolerance ) then
!         ! found a zero
!         x=bt%current_x()
!         exit
!     endif
!     ! pick next x-value
!     x=bt%next_x()
!     ! evaluate function
!     y=function(x)
!     ! tell the solver what happened
!     call bt%update(x,y)
! enddo
!
!===============================================================================
!> Module for finding the zero of a function adapted from TDEP utils
module brents_method
use kinds, only : dp

implicit none
private

integer, parameter :: lo_hugeint = huge(1)
real(dp), parameter :: lo_huge = huge(1.0_dp)

public :: brent_helper

!> helper for Brent's method of finding zeros
type brent_helper
    !> x-value, lower bound
    real(dp), private :: xlo=-lo_huge
    !> x-value, upper bound
    real(dp), private :: xhi=-lo_huge
    !> y-value, lower bound
    real(dp), private :: ylo=-lo_huge
    !> y-value, upper bound
    real(dp), private :: yhi=-lo_huge
    !> tolerance
    real(dp) :: tolerance=-lo_huge
    !> iteration counter
    integer, private :: niter=-lo_hugeint
    !> helper variables
    logical, private :: mflag=.false.

    real(dp), private :: xc=-lo_huge
    real(dp), private :: fc=-lo_huge
    real(dp), private :: fd=-lo_huge
    real(dp), private :: s =-lo_huge
    real(dp), private :: xd=-lo_huge
    real(dp), private :: fs=-lo_huge
    real(dp), private :: sd=-lo_huge
    contains
        !> initialize the helper
        procedure :: init
        !> get the current function value/argument
        procedure :: current_y
        procedure :: current_x
        !> next argument to evaluate
        procedure :: next_x
        !> update solver with results
        procedure :: update
end type

contains

!> initialize Brent's method
subroutine init(bt,xlo,xhi,ylo,yhi,tolerance)
    !> Brent's method handle
    class(brent_helper), intent(out) :: bt
    !> safe lower bound
    real(dp), intent(in) :: xlo
    !> safe upper bound
    real(dp), intent(in) :: xhi
    !> safe lower bound
    real(dp), intent(in) :: ylo
    !> safe upper bound
    real(dp), intent(in) :: yhi
    !> tolerance when to quit
    real(dp), intent(in) :: tolerance

    ! Get them in the right order
    if ( abs(ylo) .lt. abs(yhi) ) then
        bt%xlo=xlo
        bt%xhi=xhi
        bt%ylo=ylo
        bt%yhi=yhi
    else
        bt%xlo=xhi
        bt%xhi=xlo
        bt%ylo=yhi
        bt%yhi=ylo
    endif

    ! Set the helper variables
    bt%tolerance=tolerance
    bt%niter=0
    bt%xc=bt%xlo
    bt%fc=bt%ylo
    bt%mflag=.true.
    bt%s=0.0_dp
    bt%xd=0.0_dp
    bt%fs=-lo_huge
    bt%fd=-lo_huge
end subroutine

!> fetch the next argument to evaluate function at
function next_x(bt) result(x)
    !> Brent's method handle
    class(brent_helper), intent(inout) :: bt
    !> next x-value to evaluate the function at
    real(dp) :: x

    real(dp) :: f0,f1
    logical :: sflag

    ! Figure out where to put the next point
    if ( abs(bt%ylo-bt%fc) .gt. bt%tolerance .and. abs(bt%yhi-bt%fc) .gt. bt%tolerance ) then
        x=bt%xlo*bt%yhi*bt%fc/( (bt%ylo-bt%yhi)*(bt%ylo-bt%fc) ) + &
        bt%xhi*bt%ylo*bt%fc/( (bt%yhi-bt%ylo)*(bt%yhi-bt%fc) )+&
        bt%xc*bt%ylo*bt%yhi/( (bt%fc-bt%ylo)*(bt%fc-bt%yhi) )
    else
        x=bt%xhi-bt%yhi*(bt%xhi-bt%xlo)/(bt%yhi-bt%ylo)
    endif

    ! Also try to figure out is that was a bad idea, and revert to just slicing the interval in half.
    sflag=.false.
    f0=min((3*bt%xlo+bt%xhi)*0.25_dp,bt%xhi)
    f1=max((3*bt%xlo+bt%xhi)*0.25_dp,bt%xhi)
    if ( x .lt. f0      .or.  x                 .gt. f1                       ) sflag=.true.
    if ( bt%mflag       .and. abs(x-bt%xhi)     .ge. abs(bt%xhi-bt%xc)*0.5_dp ) sflag=.true.
    if ( .not.bt%mflag  .and. abs(x-bt%xhi)     .ge. abs(bt%xc-bt%xd)*0.5_dp  ) sflag=.true.
    if ( bt%mflag       .and. abs(bt%xhi-bt%xc) .lt. bt%tolerance             ) sflag=.true.
    if ( .not.bt%mflag  .and. abs(bt%xc-bt%xd)  .lt. bt%tolerance             ) sflag=.true.

    sflag=.true.
    if ( sflag ) then
        x=(bt%xlo+bt%xhi)*0.5_dp
        bt%mflag=.true.
    else
        bt%mflag=.false.
    endif
end function

!> update solver based on function value
subroutine update(bt,x,y)
    !> Brent's method handle
    class(brent_helper), intent(inout) :: bt
    !> function argument
    real(dp), intent(in) :: x
    !> function value
    real(dp), intent(in) :: y

    real(dp) :: f0

    ! Update the solver settings
    bt%xd=bt%xc
    bt%fd=bt%fc
    bt%xc=bt%xhi
    bt%fc=bt%yhi
    if ( bt%ylo*y .lt. 0.0_dp ) then
        bt%xhi=x
        bt%yhi=y
    else
        bt%xlo=x
        bt%ylo=y
    endif

    ! Maybe flip?
    if ( abs(bt%ylo) .lt. abs(bt%ylo) ) then
        f0=bt%ylo
        bt%ylo=bt%yhi
        bt%yhi=f0
        f0=bt%xlo
        bt%xlo=bt%xhi
        bt%xhi=f0
    endif

    ! Count iterations?
    bt%niter=bt%niter+1
end subroutine

!> fetch current argument
function current_x(bt) result(x)
    !> Brent's method handle
    class(brent_helper), intent(in) :: bt
    !> function argument
    real(dp) :: x

    ! check for convergence
    if ( abs(bt%xhi) .lt. abs(bt%fs) ) then
        x=bt%xhi
    else
        x=bt%s
    endif
end function

!> fetch current function value
function current_y(bt) result(y)
    !> Brent's method handle
    class(brent_helper), intent(in) :: bt
    !> function value
    real(dp) :: y

    ! check for convergence
    if ( abs(bt%xhi) .lt. abs(bt%fs) ) then
        y=bt%yhi
    else
        y=bt%fs
    endif
end function

end module
