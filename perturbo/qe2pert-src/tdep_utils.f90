!===============================================================================
! Author: jyyao <jyyao27@gmail.com>
! 
! File: tdep_utils.f90
! Create Date: 2024-12-05
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Comment:
!  utility module to the TDEP interface.
!  adapted from TDEP/src/libolle
!  Refer to the following literature:
!      O. Hellman et al., Phys Rev B 87, 104111 (2013)
!      O. Hellman and I. A. Abrikosov, Phys Rev B 88, 144301 (2013)
!      F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)
!
! Maintenance:
!================================================================================
module tdep_utils
   use kinds, only: dp
   use constants, only : pi
   implicit none
   public
   real(dp), parameter, private :: lo_tol = 1.0E-5_dp
   real(dp), parameter, private :: sq_tol = 1.0E-10_dp
   real(dp), parameter, private :: lo_huge = huge(1.0_dp)
   integer, parameter, private :: lo_hugeint = huge(1)
   integer, parameter, private :: lo_status = 0
   complex(dp), parameter :: lo_imag = (0.0_dp,1.0_dp)

   !> helper to keep track of nonlinear minimization. Adapted directly from tdep
   type nonlinear_minimization_helper
      !> which iteration are we on?
      integer :: iter = -lo_hugeint
      !> how many variables
      integer :: nx = -lo_hugeint
      !> how many directions to provide
      integer :: n_direction = -lo_hugeint
      !> which algo did I use this step?
      integer :: algo = -lo_hugeint
      !> Errors (the thing we want to minimize)
      real(dp) :: e_curr, e_prev
      !> Current x
      real(dp), dimension(:), allocatable :: x_curr
      !> Previous x
      real(dp), dimension(:), allocatable :: x_prev
      !> Current gradient
      real(dp), dimension(:), allocatable :: g_curr
      !> Previous gradient
      real(dp), dimension(:), allocatable :: g_prev
   
      !> Search direction things
      real(dp), dimension(:, :), allocatable :: sdir
      real(dp), dimension(:, :), allocatable :: sx
      real(dp), dimension(:), allocatable :: serr
      real(dp), dimension(:), allocatable :: sstep
   
      !> BFGS helpers
      real(dp), dimension(:), allocatable :: bfgs_oldx
      real(dp), dimension(:), allocatable :: bfgs_oldg
      real(dp), dimension(:), allocatable :: bfgs_s
      real(dp), dimension(:), allocatable :: bfgs_y
      real(dp), dimension(:, :), allocatable :: bfgs_B
      integer :: bfgs_n_update = -lo_hugeint
   
      !> CG helpers
      real(dp), dimension(:), allocatable :: cg_oldx
      real(dp), dimension(:), allocatable :: cg_oldg
   end type

contains
   !> Returns the atomic symbol from atomic number
   function z_to_symbol(z_nucleus) result(symbol)
      integer, intent(in) :: z_nucleus  !< Nuclear charge, e.g. 14
      character(len=2) :: symbol        !< Atomic symbol, e.g. "Si"

      character(len=2), parameter :: symbols_of_z(103) = [&
         &'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
         &'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
         &'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
         &'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
         &'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
         &'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
         &'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
         &'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
         &'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
         &'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
         &'Md', 'No', 'Lr']
      symbol = symbols_of_z(z_nucleus)
   end function

   !> Translates an atomic symbol into atomic number
   function symbol_to_z(symbol) result(z_nucleus)
      character(len=*), intent(in) :: symbol !< Atomic symbol, e.g. "Si"
      integer :: z_nucleus !< Nuclear charge, e.g. 14
      select case(trim(adjustl(symbol)))
         case('H' )  ;  z_nucleus=1
         case('He')  ;  z_nucleus=2
         case('Li')  ;  z_nucleus=3
         case('Be')  ;  z_nucleus=4
         case('B' )  ;  z_nucleus=5
         case('C' )  ;  z_nucleus=6
         case('N' )  ;  z_nucleus=7
         case('O' )  ;  z_nucleus=8
         case('F' )  ;  z_nucleus=9

         case('Ne')  ;  z_nucleus=10
         case('Na')  ;  z_nucleus=11
         case('Mg')  ;  z_nucleus=12
         case('Al')  ;  z_nucleus=13
         case('Si')  ;  z_nucleus=14
         case('P' )  ;  z_nucleus=15
         case('S' )  ;  z_nucleus=16
         case('Cl')  ;  z_nucleus=17
         case('Ar')  ;  z_nucleus=18
         case('K' )  ;  z_nucleus=19

         case('Ca')  ;  z_nucleus=20
         case('Sc')  ;  z_nucleus=21
         case('Ti')  ;  z_nucleus=22
         case('V' )  ;  z_nucleus=23
         case('Cr')  ;  z_nucleus=24
         case('Mn')  ;  z_nucleus=25
         case('Fe')  ;  z_nucleus=26
         case('Co')  ;  z_nucleus=27
         case('Ni')  ;  z_nucleus=28
         case('Cu')  ;  z_nucleus=29

         case('Zn')  ;  z_nucleus=30
         case('Ga')  ;  z_nucleus=31
         case('Ge')  ;  z_nucleus=32
         case('As')  ;  z_nucleus=33
         case('Se')  ;  z_nucleus=34
         case('Br')  ;  z_nucleus=35
         case('Kr')  ;  z_nucleus=36
         case('Rb')  ;  z_nucleus=37
         case('Sr')  ;  z_nucleus=38
         case('Y' )  ;  z_nucleus=39

         case('Zr')  ;  z_nucleus=40
         case('Nb')  ;  z_nucleus=41
         case('Mo')  ;  z_nucleus=42
         case('Tc')  ;  z_nucleus=43
         case('Ru')  ;  z_nucleus=44
         case('Rh')  ;  z_nucleus=45
         case('Pd')  ;  z_nucleus=46
         case('Ag')  ;  z_nucleus=47
         case('Cd')  ;  z_nucleus=48
         case('In')  ;  z_nucleus=49

         case('Sn')  ;  z_nucleus=50
         case('Sb')  ;  z_nucleus=51
         case('Te')  ;  z_nucleus=52
         case('I' )  ;  z_nucleus=53
         case('Xe')  ;  z_nucleus=54
         case('Cs')  ;  z_nucleus=55
         case('Ba')  ;  z_nucleus=56
         case('La')  ;  z_nucleus=57
         case('Ce')  ;  z_nucleus=58
         case('Pr')  ;  z_nucleus=59

         case('Nd')  ;  z_nucleus=60
         case('Pm')  ;  z_nucleus=61
         case('Sm')  ;  z_nucleus=62
         case('Eu')  ;  z_nucleus=63
         case('Gd')  ;  z_nucleus=64
         case('Tb')  ;  z_nucleus=65
         case('Dy')  ;  z_nucleus=66
         case('Ho')  ;  z_nucleus=67
         case('Er')  ;  z_nucleus=68
         case('Tm')  ;  z_nucleus=69

         case('Yb')  ;  z_nucleus=70
         case('Lu')  ;  z_nucleus=71
         case('Hf')  ;  z_nucleus=72
         case('Ta')  ;  z_nucleus=73
         case('W' )  ;  z_nucleus=74
         case('Re')  ;  z_nucleus=75
         case('Os')  ;  z_nucleus=76
         case('Ir')  ;  z_nucleus=77
         case('Pt')  ;  z_nucleus=78
         case('Au')  ;  z_nucleus=79

         case('Hg')  ;  z_nucleus=80
         case('Tl')  ;  z_nucleus=81
         case('Pb')  ;  z_nucleus=82
         case('Bi')  ;  z_nucleus=83
         case('Po')  ;  z_nucleus=84
         case('At')  ;  z_nucleus=85
         case('Rn')  ;  z_nucleus=86
         case('Fr')  ;  z_nucleus=87
         case('Ra')  ;  z_nucleus=88
         case('Ac')  ;  z_nucleus=89

         case('Th')  ;  z_nucleus=90
         case('Pa')  ;  z_nucleus=91
         case('U' )  ;  z_nucleus=92
         case('Np')  ;  z_nucleus=93
         case('Pu')  ;  z_nucleus=94
         case('Am')  ;  z_nucleus=95
         case('Cm')  ;  z_nucleus=96
         case('Bk')  ;  z_nucleus=97
         case('Cf')  ;  z_nucleus=98
         case('Es')  ;  z_nucleus=89

         case('Fm')  ;  z_nucleus=100
         case('Md')  ;  z_nucleus=101
         case('No')  ;  z_nucleus=102
         case('Lr')  ;  z_nucleus=103
         case default
            write(*,*) 'TROUBLE!!!, no atomic number available for '//trim(adjustl(symbol))
            write(*,*) 'Are you really sure that is an element?'
            stop
      end select
   end function symbol_to_z

   !> Determinant of a 3x3 matrix, doubles
   pure function determ_real(a) result(det)
      real(dp), intent(in) :: a(3,3) !< 3x3 matrix
      real(dp) :: det                !< the determinant
      det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
         + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) &
         + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
   end function

   !> Invert 3x3 matrix, without lapack
   pure function invert3x3matrix(m) result(n) 
      real(dp), intent(in) :: m(3,3) !< 3x3 matrix to invert
      real(dp) :: n(3,3)             !< the inverse
      real(dp) :: det                !< determinant

      det = 1.0_dp / determ_real(m)
      ! Calculate the inverse of the matrix
      n(1,1)=+det*(m(2,2)*m(3,3)-m(2,3)*m(3,2))
      n(2,1)=-det*(m(2,1)*m(3,3)-m(2,3)*m(3,1))
      n(3,1)=+det*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
      n(1,2)=-det*(m(1,2)*m(3,3)-m(1,3)*m(3,2))
      n(2,2)=+det*(m(1,1)*m(3,3)-m(1,3)*m(3,1))
      n(3,2)=-det*(m(1,1)*m(3,2)-m(1,2)*m(3,1))
      n(1,3)=+det*(m(1,2)*m(2,3)-m(1,3)*m(2,2))
      n(2,3)=-det*(m(1,1)*m(2,3)-m(1,3)*m(2,1))
      n(3,3)=+det*(m(1,1)*m(2,2)-m(1,2)*m(2,1))
   end function

   !> Clean up fractional coordinates in a way that I like, so that they become 0-1, inclusive 0 not inclusive 1. Some quick testing suggested this was faster than mod if I am close to 0-1, which almost always is the case. Also means that you should never use this for large numbers.
   elemental function clean_fractional_coordinates(x,tol) result(y)

      real(dp), intent(in) :: x !< unclean coordinate
      real(dp), intent(in) :: tol !< tolerance
      real(dp) :: y !< the clean coordinate
      y=x
      do
         if ( abs(y) .lt. tol ) y=0.0_dp
         if ( abs(y-1.0_dp) .lt. tol ) y=0.0_dp
         if ( y .gt. 1.0_dp ) y=y-1.0_dp
         if ( y .lt. 0.0_dp ) y=y+1.0_dp
         if ( y .gt. -tol .and. y .lt. 1.0_dp ) exit
      enddo
   end function

   !> find the reciprocal lattice vectors
   pure function reciprocal_basis(a) result(b)
      real(dp), dimension(3,3), intent(in) :: a !< the basis
      real(dp), dimension(3,3) :: b !< reciprocal basis
      real(dp), dimension(3) :: a1,a2,a3,b1,b2,b3
      real(dp) :: vol
      a1=a(:,1)
      a2=a(:,2)
      a3=a(:,3)
      vol=abs(determ_real(a))
      b1=cross(a2,a3)
      b2=cross(a3,a1)
      b3=cross(a1,a2)
      b(:,1)=b1/vol
      b(:,2)=b2/vol
      b(:,3)=b3/vol
   end function

   !> compute cross product 
   pure function cross(b,c) result(a)
      real(dp), dimension(3), intent(in) :: b !< first vector
      real(dp), dimension(3), intent(in) :: c !< second vector
      real(dp), dimension(3) :: a !< cross product

      a(1)=b(2)*c(3)-b(3)*c(2)
      a(2)=b(3)*c(1)-b(1)*c(3)
      a(3)=b(1)*c(2)-b(2)*c(1)
   end function


   !> From an integer array, return only the unique values
   subroutine return_unique_integers(a,u)
      integer, intent(in) :: a(:) !< list of integers
      integer, allocatable, intent(out) :: u(:) !< unique integers
      integer, allocatable :: di(:), dj(:)
      integer :: i,j,n,ctr

      n=size(a)

      ! This version is faster if we have a constant small number of unique, then it's strictly O(N).
      allocate(di(n))
      di=-lo_hugeint
      ctr=0
      ! Start checking list
      vl: do i=1,n
         ! Check if this lines is in the list of unique lines
         do j=1,ctr
            if ( a(i)-di(j) .eq. 0 ) cycle vl
         enddo
         ! If we made it here it's a new entry
         ctr=ctr+1
         di(ctr)=a(i)
      enddo vl
      ! And return the list of unique
      allocate(u(ctr))
      u=di(1:ctr)
      deallocate(di)
   end subroutine


   !> Chop real numbers during reading grid points 
   elemental function chop_real(x,tol) result(y)
      real(dp), intent(in) :: x !< number to be chopped
      real(dp), intent(in) :: tol !< tolerance
      real(dp) :: y !< the chopped number
      integer :: i
      real(dp), dimension(69), parameter :: welldefined=[&
         0.0000000000000000_dp,&
         0.1000000000000000_dp,&
        -0.1000000000000000_dp,&
         0.1111111111111111_dp,&
        -0.1111111111111111_dp,&
         0.1250000000000000_dp,&
        -0.1250000000000000_dp,&
         0.1428571428571428_dp,&
        -0.1428571428571428_dp,&
         0.1666666666666667_dp,&
        -0.1666666666666667_dp,&
         0.2000000000000000_dp,&
        -0.2000000000000000_dp,&
         0.2222222222222222_dp,&
        -0.2222222222222222_dp,&
         0.2500000000000000_dp,&
        -0.2500000000000000_dp,&
         0.2857142857142857_dp,&
        -0.2857142857142857_dp,&
         0.3000000000000000_dp,&
        -0.3000000000000000_dp,&
         0.3333333333333333_dp,&
        -0.3333333333333333_dp,&
         0.3750000000000000_dp,&
        -0.3750000000000000_dp,&
         0.4000000000000000_dp,&
        -0.4000000000000000_dp,&
         0.4285714285714285_dp,&
        -0.4285714285714285_dp,&
         0.4444444444444444_dp,&
        -0.4444444444444444_dp,&
         0.5000000000000000_dp,&
        -0.5000000000000000_dp,&
         0.5555555555555556_dp,&
        -0.5555555555555556_dp,&
         0.5714285714285714_dp,&
        -0.5714285714285714_dp,&
         0.6000000000000000_dp,&
        -0.6000000000000000_dp,&
         0.6250000000000000_dp,&
        -0.6250000000000000_dp,&
         0.6666666666666667_dp,&
        -0.6666666666666667_dp,&
         0.7000000000000000_dp,&
        -0.7000000000000000_dp,&
         0.7071067811865475_dp,&
        -0.7071067811865475_dp,&
         0.7142857142857143_dp,&
        -0.7142857142857143_dp,&
         0.7500000000000000_dp,&
        -0.7500000000000000_dp,&
         0.7777777777777778_dp,&
        -0.7777777777777778_dp,&
         0.8000000000000000_dp,&
        -0.8000000000000000_dp,&
         0.8333333333333334_dp,&
        -0.8333333333333334_dp,&
         0.8571428571428571_dp,&
        -0.8571428571428571_dp,&
         0.8660254037844386_dp,&
        -0.8660254037844386_dp,&
         0.8750000000000000_dp,&
        -0.8750000000000000_dp,&
         0.8888888888888888_dp,&
        -0.8888888888888888_dp,&
         0.9000000000000000_dp,&
        -0.9000000000000000_dp,&
         1.0000000000000000_dp,&
        -1.0000000000000000_dp]
      ! my list of well defined values. Could very well grow a little bit.
      ! starts with numbers that are common in symmetry operations
      y=x
      do i=1,69
         if ( abs(y-welldefined(i)) .lt. tol ) then
            y=welldefined(i)
            return
         endif
      enddo
   end function

   !> cutoff at small numbers for complex numbers
   elemental function chop_complex(x,tol) result(y)
      complex(dp), intent(in) :: x !< number to be chopped
      real(dp), intent(in) :: tol !< tolerance
      complex(dp) :: y !< chopped number
      real(dp) :: re,im

      re=real(x)
      im=aimag(x)
      if ( abs(re) .lt. tol ) re=0.0_dp
      if ( abs(im) .lt. tol ) im=0.0_dp
      y=cmplx(re,im,dp)
   end function

   !> find points on the sphere 
   pure subroutine points_on_sphere(r,randomseed)
      real(dp), dimension(:,:), intent(out) :: r    !< the points
      real(dp), intent(in), optional :: randomseed  !< a random 0-1 seed to shift the points around
      real(dp) :: offset,increment,x,y,z,rad,phi,rnd
      integer :: i,n
   
      if ( present(randomseed) ) then
         rnd=randomseed
      else
         rnd=0.0_dp
      endif
   
      n=size(r,2)
      r=0.0_dp
      offset=2.0_dp/(1.0_dp*n)
      increment=pi*(3.0_dp-sqrt(5.0_dp))
      do i=1,n
         ! Negative 2 here compared to what you might find, since
         ! Fortran is 1-indexed.
         y=((i-1)*offset-1)+offset*0.5_dp;
         rad=sqrt( max(1.0_dp-y**2,0.0_dp) );
         phi=mod(real(i+rnd,dp),real(n,dp))*increment;
         x=cos(phi)*rad;
         z=sin(phi)*rad;
         r(:,i)=[x,y,z]/norm2([x,y,z]);
      enddo
   end subroutine

   
   !> given a parallelepiped, what is the radius of the largest sphere that fits inside?
   pure function inscribed_sphere_in_box(box) result(r)
      real(dp), dimension(3,3), intent(in) :: box !< vectors that define the box
      real(dp) :: r !< largest sphere that fit inside box
   
      real(dp), dimension(3) :: na,nb,nc
      ! the normals of the faces of the box
      na=cross(box(:,2),box(:,3))
      nb=cross(box(:,3),box(:,1))
      nc=cross(box(:,1),box(:,2))
      na=na/norm2(na)
      nb=nb/norm2(nb)
      nc=nc/norm2(nc)
      ! distances between opposing planes
      r=lo_huge
      r=min(r,abs(dot_product(na,box(:,1))))
      r=min(r,abs(dot_product(nb,box(:,2))))
      r=min(r,abs(dot_product(nc,box(:,3))))
      r=r*0.5_dp
   end function

   !> find the minimal radius of a sphere that constains the box
   pure function bounding_sphere_of_box(box) result(r)
      real(dp), dimension(3,3), intent(in) :: box !< vectors that define the box
      real(dp) :: r !< radius of smallest sphere that contains the box
      real(dp), dimension(3) :: a,b,c,v
      real(dp), dimension(4) :: diags
   
      a=box(:,1)
      b=box(:,2)
      c=box(:,3)
      ! had to do this in an awkward way to avoid segfaults on intel compilers. Odd.
      v=a-b-c
      diags(1)=dot_product(v,v)
      v=a+b-c
      diags(2)=dot_product(v,v)
      v=a-b+c
      diags(3)=dot_product(v,v)
      v=a+b+c
      diags(4)=dot_product(v,v)
      r=maxval(diags)
      r=sqrt(r)*0.5_dp
   end function

   !> Given a certain cell, and certain repetitions of it, make the repetitions 
   !! larger in the direction that gets a bigger inscribed sphere in there.
   function increment_dimensions(dimin, box) result(dimut)
      integer, dimension(3), intent(in) :: dimin !< reptitions (half of the length of the box)
      real(dp), dimension(3, 3), intent(in) :: box !< lattice vectors of the unit cell
      integer, dimension(3) :: dimut
      real(dp), dimension(3, 3) :: m0
      integer, dimension(3) :: di
      integer :: i, j, k
      real(dp) :: f0, f1, ff0
   
      ! try with the sphere thing. First get a baseline
      do j = 1, 3
          m0(:, j) = box(:, j)*(2*dimin(j) + 1)
      enddo
      ff0 = inscribed_sphere_in_box(m0)
      ! Increment the dimension that gives the biggest increase in radii
      f0 = 0.0_dp
      dimut = 0
      do i = 1, 3
         di = dimin
         di(i) = di(i) + 1
         do j = 1, 3
            m0(:, j) = box(:, j)*(2*di(j) + 1)
         enddo
         f1 = inscribed_sphere_in_box(m0)
         if (f1 .gt. f0 .and. abs(f1 - ff0) .gt. lo_tol) then
            dimut = di
            f0 = f1
         endif
      enddo
   
      ! if nothing helped, increment the lowest number
      if (dimut(1) .eq. 0) then
         j = lo_hugeint
         do i = 1, 3
            if (di(i) .lt. j) then
               j = di(i)
               k = i
            endif
         enddo
         dimut = dimin
         dimut(k) = dimut(k) + 1
      endif
   end function


   !> Computes the dipole-dipole interaction H for an anistropic medium
   pure function anisotropic_dipole_int(x, y, inveps) result(H)
      real(dp), dimension(3), intent(in) :: x !< lambda*r of dipoles
      real(dp), intent(in) :: y               !< lambda*sqrt(r*inveps*r)
      real(dp), dimension(3, 3), intent(in) :: inveps !< inverse dielectric tensor
      real(dp), dimension(3, 3) :: H          !< dipole-dipole interaction matrix 
      real(dp) :: f0, f1
      real(dp) :: erfcy, invy2, invy3, expminvy2
      real(dp), parameter :: twooversqrtpi = 1.128379167095513_dp
      ! Keep everything at 0 if I'm at a self term
      if (abs(y) .lt. 1.0E-10_dp) then
         H = 0.0_dp
         return
      endif
      ! Pre-calculate f0 and f1 factors
      erfcy = erfc(y)
      invy2 = 1.0_dp/(y*y)
      invy3 = 1.0_dp/(y*y*y)
      expminvy2 = exp(-y*y)
      f0 = invy2*(3*erfcy*invy3 + twooversqrtpi*expminvy2*(3*invy2 + 2))
      f1 = erfcy*invy3 + twooversqrtpi*expminvy2*invy2
      ! And H_ij = f0 * ri * rj - f1 * inveps_ij
      H(1, 1) = x(1)*x(1)*f0 - inveps(1, 1)*f1
      H(2, 1) = x(2)*x(1)*f0 - inveps(2, 1)*f1
      H(3, 1) = x(3)*x(1)*f0 - inveps(3, 1)*f1
      H(1, 2) = x(1)*x(2)*f0 - inveps(1, 2)*f1
      H(2, 2) = x(2)*x(2)*f0 - inveps(2, 2)*f1
      H(3, 2) = x(3)*x(2)*f0 - inveps(3, 2)*f1
      H(1, 3) = x(1)*x(3)*f0 - inveps(1, 3)*f1
      H(2, 3) = x(2)*x(3)*f0 - inveps(2, 3)*f1
      H(3, 3) = x(3)*x(3)*f0 - inveps(3, 3)*f1
   end function

   !> find the square of a norm
   pure function sqnorm(a) result(nrm)
      real(dp), dimension(3), intent(in) :: a !< vector
      real(dp) :: nrm !< squared norm

      nrm=a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
   end function
   
   !> quick integer sorting
   pure subroutine quick_sort_int(list,order)
      integer, dimension(:), intent(inout) :: list !< list of integer
      integer, dimension(:), intent(inout), optional :: order !< sorted order
      integer, parameter  :: max_simple_sort_size=6
      integer :: k
      ! initialize the order thingy
      if ( present(order) ) then
         do k=1,size(list)
            order(k) = k
         enddo
      endif
   
      if ( present(order) ) then
         call qsint(list,1,size(list),order)
      else
         call qsint(list,1,size(list))
      endif
   
      contains
      pure recursive subroutine qsint(buf,left_end,right_end,ind)
         integer, dimension(:), intent(inout) :: buf
         integer, dimension(:), intent(inout), optional :: ind
         integer, intent(in) :: left_end, right_end
   
         integer :: i,j,reference,tmp
   
         if (right_end < left_end + max_simple_sort_size) then
            ! Use interchange sort for small lists
            if ( present(ind) ) then
               call interchange_sort_int(buf,left_end, right_end,ind)
            else
               call interchange_sort_int(buf,left_end, right_end)
            endif
         else
            ! Use partition quicksort
            reference = buf((left_end+right_end)/2)
            i=left_end-1
            j=right_end+1
            do
               ! Scan list from left end until element >= reference is found
               do
                  i = i + 1
                  if ( buf(i) >= reference) exit
               enddo
               ! Scan list from right end until element <= reference is found
               do
                  j = j - 1
                  if ( buf(j) <= reference) exit
               enddo
               if (i < j) then
                  ! Swap two out-of-order elements
                  tmp=buf(i)
                  buf(i)=buf(j)
                  buf(j)=tmp
                  if ( present(order) ) then
                     tmp=ind(i)
                     ind(i)=ind(j)
                     ind(j)=tmp
                  endif
               elseif (i == j) then
                  i = i + 1
                  exit
               else
                  exit
               endif
            enddo
            if ( present(ind) ) then
               if (left_end < j) call qsint(buf,left_end,j,ind)
               if (i < right_end) call qsint(buf,i,right_end,ind)
            else
               if (left_end < j) call qsint(buf,left_end,j)
               if (i < right_end) call qsint(buf,i,right_end)
            endif
         endif
      end subroutine
   end subroutine

   
   !> Interchange sort for small integer arrays
   pure subroutine interchange_sort_int(list,left_end,right_end,order)
      integer, dimension(:), intent(inout) :: list !< list to be sorted
      integer, intent(in) :: left_end, right_end !< left and right ends
      integer, dimension(:), intent(inout), optional :: order !< sorted order
   
      integer :: i,j,tmp
   
      if ( present(order) ) then
         do i=left_end,right_end-1
         do j=i+1,right_end
            if ( list(i) > list(j) ) then
               tmp=list(i)
               list(i)=list(j)
               list(j)=tmp
               tmp=order(i)
               order(i)=order(j)
               order(j)=tmp
            endif
         enddo
         enddo
      else
         do i=left_end,right_end-1
         do j=i+1,right_end
            if ( list(i) > list(j) ) then
               tmp=list(i)
               list(i)=list(j)
               list(j)=tmp
            endif
         enddo
         enddo
      endif
   end subroutine

   !> quicksort, but with vectors instead.
   pure subroutine quick_sortv(A,ind,tol)
       real(dp), dimension(:,:), intent(inout) :: A !< list of vectors
       integer, dimension(:), intent(out), optional :: ind !< index that sorts A
       real(dp), intent(in) :: tol !< tolerance
       real(dp) :: diff
       real(dp), dimension(:), allocatable :: dr
       integer, dimension(:), allocatable :: di1
       integer :: i,n,m,col,ctr,istart
   
       ! size of the input
       n=size(A,2)
       m=size(A,1)
       if ( n .eq. 1 ) then
          ! might get weird, and also, we have to do nothing
          if ( present(ind) ) ind=1
          return
       endif
       ! some temporary things
       allocate(dr(n))
       allocate(di1(n))
   
       dr=A(1,:) ! all the x-components
       call quick_sort_real(dr,di1)
       ! rearrange A so that the first element is in descending order
       A=A(:,di1)
       if ( present(ind) ) then
          do i=1,n
             ind(i)=i
          enddo
          ind=ind(di1)
       endif
   
       ! Now rearrange stuff per each column, in smaller and smaller groups
       dr=0.0_dp
       di1=0
       cloop: do col=2,m
          istart=0
          ctr=1
          dr=0.0_dp
          dr(ctr)=A(col,1)
          do i=2,n+1
             if ( i .le. n ) then
                diff=sum(abs(A(1:col-1,i)-A(1:col-1,i-1)))
             else
                diff=lo_huge
             endif
             if ( diff .lt. tol .and. i .le. n ) then
                ! one more of the same element
                ctr=ctr+1
                dr(ctr)=A(col,i)
             else
                ! now we are finished with this group, sort it
                ! we only have to do that if there is more than one row in that group.
                if ( ctr .gt. 1 ) then
                   call quick_sort_real(dr(1:ctr),di1(1:ctr))
                   ! Rearrange the array
                   A(:,istart+1:istart+ctr)=A(:,di1(1:ctr)+istart)
                   if ( present(ind) ) ind(istart+1:istart+ctr)=ind(istart+di1(1:ctr))
                endif
                ! reset counters and stuff
                if ( i .le. n ) then
                   ctr=1
                   istart=i-1
                   di1=0
                   dr=0.0_dp
                   dr(1)=A(col,i)
                endif
             endif
          enddo
       enddo cloop
       deallocate(dr)
       deallocate(di1)
   end subroutine

   
   !> qsort for double precision
   pure subroutine quick_sort_real(list,order)
      real(dp), dimension(:), intent(inout)  :: list !< array to sort
      integer, dimension(:), intent(inout), optional :: order !< optionally, an index array
      integer :: k
   
      if ( present(order) ) then
         do k=1,size(list)
             order(k) = k
         enddo
      endif
   
      if ( present(order) ) then
         call qsreal(1,size(list),list,order)
      else
         call qsreal(1,size(list),list)
      endif
      contains
   
      pure recursive subroutine qsreal(left_end,right_end,buf,ind)
         integer, intent(in) :: left_end,right_end
         real(dp), dimension(:), intent(inout) :: buf
         integer, dimension(:), intent(inout), optional :: ind
   
         integer, parameter :: max_simple_sort_size=6
         real(dp) :: reference, temp
         integer :: i,j,itemp
   
         if (right_end < left_end + max_simple_sort_size) then
            ! Use interchange sort for small lists
            if ( present(order) ) then
               call interchange_sort_real(buf, left_end, right_end, ind)
            else
               call interchange_sort_real(buf, left_end, right_end)
            endif
         else
            ! Use partition ("quick") sort
            reference = buf((left_end + right_end)/2)
            i = left_end - 1; j = right_end + 1
            do
               ! Scan list from left end until element >= reference is
               ! found
               do
                  i = i + 1
                  if (buf(i) >= reference) exit
               enddo
               ! Scan list from right end until element <= reference is found
               do
                  j = j - 1
                  if (buf(j) <= reference) exit
               enddo
   
               if (i < j) then
                  ! Swap two out-of-order elements
                  temp = buf(i); buf(i) = buf(j); buf(j) = temp
                  if ( present(order) ) then
                     itemp = ind(i); ind(i) = ind(j); ind(j) = itemp
                  endif
               elseif (i == j) then
                  i = i + 1
                  exit
               else
                  exit
               endif
            enddo
   
            if ( present(ind) ) then
               if (left_end < j) call qsreal(left_end, j,buf,ind)
               if (i < right_end) call qsreal(i, right_end,buf,ind)
            else
               if (left_end < j) call qsreal(left_end, j,buf)
               if (i < right_end) call qsreal(i, right_end,buf)
            endif
         endif
      end subroutine
   end subroutine

   
   !> Interchange sort to be used for small real arrays
   pure subroutine interchange_sort_real(list, left_end, right_end, order)
      real(dp), dimension(:), intent(inout) :: list !< array to sort
      integer, intent(in) :: left_end, right_end    !< left and right ends
      integer, dimension(:), intent(inout), optional :: order !< sorted order
   
      real(dp) :: temp
      integer :: i,j,itemp
   
      if ( present(order) ) then
         do i = left_end, right_end-1
         do j = i+1, right_end
            if (list(i) > list(j)) then
               temp = list(i); list(i) = list(j); list(j) = temp
               itemp = order(i); order(i) = order(j); order(j) = itemp
            endif
         enddo
         enddo
      else
         do i = left_end, right_end-1
         do j = i+1, right_end
            if (list(i) > list(j)) then
               temp = list(i)
               list(i) = list(j)
               list(j) = temp
            endif
         enddo
         enddo
      endif
   end subroutine


   !> gradient for fixing the Born charges
   subroutine Z_error_gradient(coeff, x, delta, grad, wf, wZ, wD, dm)
      real(dp), dimension(:, :), intent(in) :: coeff !< coefficient matrix to convert from irreducible
      real(dp), dimension(:), intent(in) :: x !< irreducible
      real(dp), intent(in) :: delta !< delta in stencil
      real(dp), dimension(:), intent(out) :: grad !< gradient
      real(dp), dimension(:), intent(inout) :: wf !< workspace
      real(dp), dimension(:, :, :), intent(inout) :: wZ
      real(dp), dimension(:, :, :, :), intent(inout) :: wD
      real(dp), dimension(:, :, :, :), intent(in) :: dm
   
      integer, parameter :: order=3   ! fixed order
      real(dp), dimension(2*order + 1, 4) :: sc
      real(dp), dimension(2*order + 1) :: vl
      real(dp), dimension(size(x)) :: wx
      integer :: i, j
   
      do i = 1, size(x)
         call centraldifference(x(i), delta, sc)
         do j = 1, 2*order + 1
            wx = x
            wx(i) = sc(j, 1)
            vl(j) = Z_error_function(coeff, wx, wf, wZ, wD, dm)
         enddo
         grad(i) = sum(sc(:, 2)*vl)
      enddo
   end subroutine


   !> Calculate the squared violation of hermiticity
   function Z_error_function(coeff, x, wf, wZ, wD, dm) result(hermerr)
      real(dp), dimension(:, :), intent(in) :: coeff !< coefficient matrix to convert from irreducible
      real(dp), dimension(:), intent(in) :: x !< irreducible
      real(dp), dimension(:), intent(inout) :: wf !<work arrays
      real(dp), dimension(:, :, :), intent(inout) :: wZ
      real(dp), dimension(:, :, :, :), intent(inout) :: wD
      real(dp), dimension(:, :, :, :), intent(in) :: dm
      real(dp) :: hermerr !< squared error
   
      real(dp), dimension(3, 3) :: m0
      integer :: a1, a2, i, j, ii, jj
   
      ! Unflatten the Born charges
      wf = matmul(coeff, x)
      do a1 = 1, size(wZ, 3)
         wZ(:, :, a1) = unflatten_2tensor(wf((a1 - 1)*9 + 1:a1*9))
      enddo
   
      ! Multiply in the charges. This can probably be done smarter.
      wD = 0.0_dp
      do a2 = 1, size(wZ, 3)
      do a1 = 1, size(wZ, 3)
         do j = 1, 3
         do i = 1, 3
            do jj = 1, 3
            do ii = 1, 3
               wD(i, j, a1, a2) = wD(i, j, a1, a2) + wZ(i, ii, a1)*wZ(j, jj, a2)*dm(ii, jj, a1, a2)
            enddo
            enddo
         enddo
         enddo
      enddo
      enddo
   
      hermerr = 0.0_dp
      do a1 = 1, size(wZ, 3)
         m0 = 0.0_dp
         do a2 = 1, size(wZ, 3)
            m0 = m0 + wD(:, :, a1, a2)
         enddo
         hermerr = hermerr + sum(abs(m0 - transpose(m0)))**2
      enddo
      hermerr = hermerr/((size(wZ, 3)*3)**2)
   end function


   !> initialize the nonlinear helper
   subroutine init_nonlinear_helper(mh, nx, x0, g0, initial_hessian)
       type(nonlinear_minimization_helper), intent(out) :: mh !< helper
       integer, intent(in) :: nx !< number of variables
       real(dp), dimension(:), intent(in) :: x0 !< starting x
       real(dp), dimension(:), intent(in) :: g0 !< starting gradient
       real(dp), intent(in) :: initial_hessian !< starting Hessian
       integer :: i
   
       mh%iter = 0
       mh%nx = nx
       mh%n_direction = 3 ! how many different search directions to provide
       allocate(mh%sdir(nx, mh%n_direction))
       allocate(mh%sx(nx, mh%n_direction))
       allocate(mh%serr(mh%n_direction))
       mh%sdir = 0.0_dp
       mh%sx = 0.0_dp
       mh%serr = 0.0_dp
   
       allocate(mh%x_curr(mh%nx))
       allocate(mh%x_prev(mh%nx))
       allocate(mh%g_curr(mh%nx))
       allocate(mh%g_prev(mh%nx))
       mh%x_curr = x0
       mh%x_prev = x0
       mh%g_prev = g0
       mh%g_prev = g0
   
       ! BFGS helper things
       allocate(mh%bfgs_oldg(nx))
       allocate(mh%bfgs_oldx(nx))
       allocate(mh%bfgs_s(nx))
       allocate(mh%bfgs_y(nx))
       allocate(mh%bfgs_B(nx, nx))
       mh%bfgs_oldg = g0
       mh%bfgs_oldx = x0
       mh%bfgs_s = 0.0_dp
       mh%bfgs_y = 0.0_dp
       mh%bfgs_B = 0.0_dp
       do i = 1, nx
          mh%bfgs_B(i, i) = initial_hessian
       enddo
       mh%bfgs_n_update = 0
   
       ! CG helper things
       allocate(mh%cg_oldg(nx))
       allocate(mh%cg_oldx(nx))
       mh%cg_oldg = g0
       mh%cg_oldx = g0
   end subroutine


   !> get a sensible step size to seed the linear search?
   function get_step(mh, searchdir, gradient, y) result(step)
      type(nonlinear_minimization_helper), intent(in) :: mh !< helper
      real(dp), dimension(:), intent(in) :: searchdir !< direction to look
      real(dp), dimension(:), intent(in) :: gradient !< gradient at this point
      real(dp), intent(in) :: y !< function value at this point
      real(dp) :: step !< step size
      real(dp) :: A, B, f0
   
      ! Try a few different things?
      A = dot_product(searchdir, matmul(mh%bfgs_B, searchdir))*0.5_dp
      B = dot_product(searchdir, gradient)
   
      f0 = B**2 - 4*A*y
      ! try to minimize the quadratic equation?
      if (f0 .gt. 1.0E-30_dp) then
          ! maybe use quadratic equation?
          if (sqrt(f0) .gt. b) then
              ! yup, makes sense
              step = (sqrt(f0) - B)/(2*A)
              return
          endif
      endif
   
      ! quadratic no good. Linear?
      if (B .lt. -1.0E-15_dp) then
          step = -y/B
          return
      endif
   
      ! nothing worked. Hmm. Pick 1.
      step = 1.0_dp
   end function


   !> get a sensible derivative delta?
   function get_delta(mh, y, prefactor) result(delta)
      type(nonlinear_minimization_helper), intent(in) :: mh !< helper
      real(dp), intent(in) :: y !< current function value
      real(dp), intent(in) :: prefactor !< scaling factor
      real(dp) :: delta !< stencil delta
   
      real(dp) :: f0
   
      ! Norm of the Hessian or something like that
      f0 = norm2(mh%bfgs_B)/(mh%nx**2)
      ! Probably need some tolerance here?
      if (f0 .gt. 1.0E-15_dp) then
         delta = prefactor*sqrt(y/f0)
      else
         delta = prefactor*sqrt(y/1.0E-15_dp)
      endif
   end function
   
   !> pick the next x value in dipole calculation 
   subroutine next_x(mh, x1, e1)
      type(nonlinear_minimization_helper), intent(inout) :: mh !< helper
      real(dp), dimension(:), allocatable, intent(inout) :: x1 !< next position
      real(dp), intent(out) :: e1 !< lowest error
      real(dp) :: f0
      integer :: i
   
      f0 = lo_huge
      do i = 1, mh%n_direction
         if (mh%serr(i) .lt. f0) then
            f0 = mh%serr(i)
            x1 = mh%sx(:, i)
            e1 = f0
            mh%algo = i
         endif
      enddo
   end subroutine

   !> fetch search directions (yes several candidates)
   subroutine get_searchdir(mh, x1, g1, e1)
      type(nonlinear_minimization_helper), intent(inout) :: mh !< helper
      real(dp), dimension(:), intent(in) :: x1 !< updated positions
      real(dp), dimension(:), intent(in) :: g1 !< updated gradient
      real(dp), intent(in) :: e1 !< updated error
      real(dp), parameter :: betatol = 1.0E-20_dp
      real(dp), parameter :: eigvaltol = 1.0E-10_dp
      real(dp) :: beta1, beta2
      integer :: i, j
   
      ! Update where we are?
      mh%x_prev = mh%x_curr
      mh%g_prev = mh%g_curr
      mh%e_prev = mh%e_curr
      mh%x_curr = x1
      mh%g_curr = g1
      mh%e_curr = e1
   
      if (mh%iter .eq. 0) then
         ! If it's the first step, all the directions will be steepest descent
         do i = 1, mh%n_direction
             mh%sdir(:, i) = -g1
         enddo
      else
         ! This is a normal step, after the first.
         ! Try to update the search directions in fancy ways. First it's just
         ! steepest descent:
         steepdesc: block
             mh%sdir(:, 1) = -g1
         end block steepdesc
         ! Second one is BFGS
         bfgs: block
            real(dp), dimension(:, :), allocatable :: wB, eigvec
            real(dp), dimension(:), allocatable :: eigval
   
            allocate(eigvec(mh%nx, mh%nx))
            allocate(eigval(mh%nx))
            allocate(wB(mh%nx, mh%nx))
            eigvec = 0.0_dp
            eigval = 0.0_dp
            wB = 0.0_dp
   
            mh%bfgs_y = mh%bfgs_oldg - g1
            mh%bfgs_s = mh%bfgs_oldx - x1
            ! Start with a straight update
            beta1 = dot_product(mh%bfgs_y, mh%bfgs_s)
            beta2 = dot_product(mh%bfgs_s, matmul(mh%bfgs_B, mh%bfgs_s))
            if (abs(beta1) .gt. betatol .and. beta2 .gt. betatol) then
               ! Ok numerics seem fine, try BFGS update
               mh%bfgs_s = matmul(mh%bfgs_B, mh%bfgs_s)
               wB = 0.0_dp
               do i = 1, mh%nx
               do j = 1, mh%nx
                  wB(i, j) = mh%bfgs_B(j, i) + mh%bfgs_y(i)*mh%bfgs_y(j)/beta1 - mh%bfgs_s(i)*mh%bfgs_s(j)/beta2
               enddo
               enddo
               ! Check the eigenvalues
               call real_symmetric_eigenvalues_eigenvectors(wB, eigval, eigvec)
               if (minval(eigval) .gt. eigvaltol) then
                  ! This was a proper update, update everything
                  mh%bfgs_B = wB
                  mh%bfgs_oldg = g1
                  mh%bfgs_oldx = x1
                  mh%bfgs_n_update = mh%bfgs_n_update + 1
               endif
            endif
   
            ! Anyhow, solve for the search direction
            call linear_least_squares(mh%bfgs_B, -g1, mh%sdir(:, 2))
   
            deallocate(eigvec)
            deallocate(eigval)
            deallocate(wB)
         end block bfgs
   
          ! Third one is conjugate gradient
         cg: block
            real(dp) :: f0, f1
   
            ! Start with Fletcher-Reeves?
            f0 = dot_product(g1, g1)
            f1 = dot_product(mh%cg_oldg, mh%cg_oldg)
            if (f1 .gt. betatol) then
               mh%sdir(:, 3) = -g1 - (f0/f1)*mh%cg_oldg
               mh%cg_oldg = g1
               mh%cg_oldx = x1
            else
               mh%sdir(:, 3) = -g1
            endif
            ! Maybe something else when things are annoying.
         end block cg
      endif
   
      ! And update the number of iterations
      mh%iter = mh%iter + 1
   end subroutine

   !> flatten a 2-tensor
   pure function flatten_2tensor(m) result(fm)
      real(dp), dimension(3,3), intent(in) :: m !< 2-tensor 
      real(dp), dimension(9) :: fm              !< flattened array

      integer :: i,j,l
      l=0
      do i=1,3
      do j=1,3
          l=l+1
          fm(l)=m(i,j)
      enddo
      enddo
   end function


   !> put an array in 1D into a 2-tensor
   pure function unflatten_2tensor(fm) result(m)
      real(dp), dimension(9), intent(in) :: fm !< 1D array
      real(dp), dimension(3,3) :: m            !< 2-tensor

      integer :: i,j,l
      l=0
      do i=1,3
      do j=1,3
          l=l+1
          m(i,j)=fm(l)
      enddo
      enddo
   end function

   !> compute abscissas and weights of central difference
   subroutine centraldifference(x0,h,sc)
      real(dp), intent(in) :: x0 !< point to differentiate around
      real(dp), intent(in) :: h !< step
      real(dp), dimension(:,:), intent(out) :: sc !< abscissas and weights. 
      !! Should be dimension (2*n+1,4), where the columns are abscissa, 
      !! and weights for 1st,2nd and 3rd derivatives
   
      sc=0.0_dp
      sc(1,:)=[-3.0_dp, -1.666666666666667e-2_dp, 1.111111111111111e-2_dp, 1.250000000000000e-1_dp];
      sc(2,:)=[-2.0_dp, 1.500000000000000e-1_dp, -1.500000000000000e-1_dp, -1.000000000000000_dp];
      sc(3,:)=[-1.0_dp, -7.500000000000000e-1_dp, 1.500000000000000_dp, 1.625000000000000_dp];
      sc(4,:)=[0.0_dp, 0.0_dp, -2.722222222222222_dp, 0.0_dp];
      sc(5,:)=[1.0_dp, 7.500000000000000e-1_dp, 1.500000000000000_dp, -1.625000000000000_dp];
      sc(6,:)=[2.0_dp, -1.500000000000000e-1_dp, -1.500000000000000e-1_dp, 1.000000000000000_dp];
      sc(7,:)=[3.0_dp, 1.666666666666667e-2_dp, 1.111111111111111e-2_dp, -1.250000000000000e-1_dp];
      ! Now adjust it to the interval we care about
      sc(:,1)=sc(:,1)*h+x0
      sc(:,2)=sc(:,2)/h
      sc(:,3)=sc(:,3)/h**2
      sc(:,4)=sc(:,4)/h**3
   end subroutine


   !> Real symmetric eigenvalue problem, with careful checks and some neat options
   subroutine real_symmetric_eigenvalues_eigenvectors(A,eigenvalues,eigenvectors)
      real(dp), dimension(:,:), intent(in) :: A !< matrix to be decomposed
      real(dp), dimension(:), intent(out) :: eigenvalues !< eigenvalues
      real(dp), dimension(:,:), intent(out), optional :: eigenvectors !< eigenvectors
   
      real(dp), dimension(:,:), allocatable :: wA
      real(dp) :: tol
      integer :: n
   
      ! Size of problem
      n=size(A,1)
      tol = sq_tol
   
      ! Zero out the things to be returned
      eigenvalues=0.0_dp
      if ( present(eigenvectors) ) then
         eigenvectors=0.0_dp
      endif
      ! And a copy of the matrix
      allocate(wA(n,n))
      wA=A
   
      ! Actually diagonalize
      diag: block
         ! Dummy things for lapack
         character(len=1) :: o_uplo,jobz
         real(dp), allocatable :: work(:),rwork(:)
         real(dp) :: s_work(1),o_z(1),s_rwork(1),o_vl,o_vu
         integer, allocatable :: iwork(:),o_isuppz(:)
         integer :: s_iwork(1),o_il,o_iu,o_m,o_info,lwork,lrwork,liwork,l_stat_alloc
   
         o_info=0
         o_uplo = 'U'
         o_vl = -lo_huge
         o_vu = lo_huge
         if( present(eigenvectors) ) then
            jobz = 'V'
         else
            jobz = 'N'
         endif
         o_iu = n
         l_stat_alloc = 0
         if( .not.present(eigenvectors) ) then
            allocate(o_isuppz(1))
            o_isuppz=0
         else
            allocate(o_isuppz(2*n), stat=l_stat_alloc)
         endif
   
         liwork = -1
         lrwork = -1
         lwork = -1
         if ( present(eigenvectors) ) then
            call dsyevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,n,0.0_dp,o_m,eigenvalues,&
                 eigenvectors,n,o_isuppz,s_work,lwork,s_iwork,liwork,o_info)
         else
            call dsyevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,n,0.0_dp,o_m,eigenvalues,&
                 o_z,1,o_isuppz,s_work,lwork,s_iwork,liwork,o_info)
         endif
         if (o_info == 0 ) then
            liwork = int(s_iwork(1))
            lrwork = int(s_rwork(1))
            lwork = int(s_work(1))
            allocate(iwork(liwork))
            allocate(rwork(lrwork))
            allocate(work(lwork))
            if ( present(eigenvectors) ) then
                call dsyevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,o_iu,0.0_dp,o_m,&
                     eigenvalues,eigenvectors,n,o_isuppz,work,lwork,iwork,liwork,o_info)
            else
                call dsyevr(jobz,'A','U',n,wA,n,o_vl,o_vu,o_il,o_iu,0.0_dp,o_m,&
                     eigenvalues,o_z,1,o_isuppz,work,lwork,iwork,liwork,o_info)
            endif
            deallocate(iwork)
            deallocate(rwork)
            deallocate(work)
         endif
         if ( allocated(o_isuppz) ) deallocate(o_isuppz)
         if( o_info <= -1000 ) then
            call xerbla('HEEVR',-o_info)
            call errore('real_symmetric_eigenvalues_eigenvectors', 'zheevr exit code', 1)
         endif
         ! Ok, now w should hold the eigenvalues, eigenvectors in eigenvectors.
      end block diag
   
      ! And chop off annoying things
      eigenvalues=chop_real(eigenvalues,tol)
      if ( present(eigenvectors) ) eigenvectors=chop_real(eigenvectors,tol)
   end subroutine


   !> Gram-Schmidt orthoginalization of a real matrix
   pure subroutine real_gram_schmidt(X)
      real(dp), dimension(:,:), intent(inout) :: X !< matrix to orthoginalize. Result should be X^TX=I
      real(dp), dimension(:,:), allocatable :: Q,R
      integer :: i,k,nr,nc
   
      nr=size(X,1)
      nc=size(X,2)
      allocate(Q(nr,nc))
      allocate(R(nc,nc))
      Q=0.0_dp
      R=0.0_dp
      do k = 1,nc
         Q(:,k) = X(:,k)
         do i = 1,k-1
            R(i,k)=dot_product(Q(:,i),Q(:,k))
            Q(:,k)=Q(:,k)-R(i,k)*Q(:,i)
         enddo
         R(k,k) = norm2(Q(:,k))
         if ( abs(R(k,k)) .gt. sq_tol ) then
            Q(:,k)=Q(:,k)/R(k,k)
         else
            Q(:,k)=0.0_dp
         endif
      enddo
      X=chop_real(Q,1.0E-13_dp)
      deallocate(Q,R)
   end subroutine


   !> wrapper around LAPACK so that it is not destructive. Solves Ax=B, perhaps with constraints such that Cx=0, with a bunch of fancy options.
   subroutine linear_least_squares(A,B,x)
       real(dp), dimension(:,:), intent(in) :: A !< matrix A
       real(dp), dimension(:), intent(in) :: B !< vector B
       real(dp), dimension(:), intent(out) :: x !< solution
   
       real(dp), dimension(:,:), allocatable :: wA,wC
       real(dp), dimension(:), allocatable :: sol,vB
       real(dp) :: tol
       integer, dimension(:), allocatable :: xind
       integer :: nx
   
       init: block
          real(dp), dimension(:,:), allocatable :: vC
          integer :: rankA,i,j,l,n
   
          ! First thing to do is figure out a whole bunch of heuristics. Like how do I want to solve it, what variables to solve for,
   
          ! Initial sanity test for the constraints
   
          ! Sanity test, double-check the dimensions?
          if ( size(A,2) .ne. size(x) ) then
             call errore('linear_least_squares', 'Dimension mismatch for linear least squares', 1)
          endif
          if ( size(A,1) .ne. size(B) ) then
             call errore('linear_least_squares', 'Dimension mismatch for linear least squares', 1)
          endif
   
          ! Now I know the number of variables
          nx=size(A,2)
   
          ! Just a straight copy
          allocate(wA(size(A,1),size(A,2)))
          wA=A
          ! And a copy of the prediction vector
          allocate(vB(size(B)))
          vB=B
   
          ! Some sanity checks
          if ( size(wA,1) .ne. size(vB) ) then
             call errore('linear_least_squares', 'Incompatible dimensions for linear least squares', 1)
          endif
   
          tol=sq_tol
   
          ! And some space for the solution
          allocate(sol(nx))
          sol=0.0_dp
          X=0.0_dp
       end block init
   
       ! Now actually solve it. A few variants.
       ! This is the normal, SVD-based least squares solution with no constraints.
       sl1: block
          real(dp), dimension(:,:), allocatable :: wB
          real(dp), dimension(:), allocatable :: work,o_s
          real(dp) :: s_work(1),o_rcond
          integer :: nrhs,lwork,o_rank,lda,ldb,m,n
   
          nrhs = 1
          lwork = -1
          m = size(wA,1)
          n = size(wA,2)
          o_rcond=100*epsilon(1.0_dp)
          allocate(wB(max(m,n),1))
          allocate(o_s(min(m,n)))
          o_s=0.0_dp
          wB=0.0_dp
          wB(1:m,1)=vB
          lda = max(1,size(wA,1))
          ldb = max(1,size(wB,1))
          call dgelss(m,n,nrhs,wA,lda,wB,ldb,o_s,o_rcond,o_rank,s_work,lwork,lo_status)
          if ( lo_status .ne. 0 ) then
             call errore('linear_least_squares', 'dgelss exit status', 1)
          endif
          lwork = int(anint( s_work(1)))
          allocate(work(lwork))
          work=0.0_dp
          call dgelss(m,n,nrhs,wA,lda,wB,ldb,o_s,o_rcond,o_rank,work,lwork,lo_status)
          if ( lo_status .ne. 0 ) then
             call errore('linear_least_squares', 'dgelss exit status', 1)
          endif
          sol=wB(1:nx,1)
          x=sol
          deallocate(wA)
          deallocate(wB)
          deallocate(o_s)
          deallocate(work)
       end block sl1
   
       X=chop_real(X,tol*1.0E-2_dp)
   
   end subroutine

   !> vector outer (dyadic) product, three-dimensional
   pure function real_outerproduct(a,b) result(m)
      real(dp), dimension(3), intent(in) :: a !< first vector
      real(dp), dimension(3), intent(in) :: b !< second vector
      real(dp), dimension(3,3) :: m !< cross product
   
      m(1,1)=a(1)*b(1)
      m(2,1)=a(1)*b(2)
      m(3,1)=a(1)*b(3)
      m(1,2)=a(2)*b(1)
      m(2,2)=a(2)*b(2)
      m(3,2)=a(2)*b(3)
      m(1,3)=a(3)*b(1)
      m(2,3)=a(3)*b(2)
      m(3,3)=a(3)*b(3)
   end function

end module tdep_utils
