!===============================================================================
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Comment:
!  interface to the TDEP code.
!  mainly interact with the following code within TDEP/src/libolle
!  type_crystalstructure.f90
!  type_forceconstant_secondorder.f90
!  type_forceconstant_secondorder_io.f90
!  type_forceconstant_secondorder_loto.f90
!  type_forceconstant_secondorder_dynamicalmatrix.f90
!  Refer to the following literature:
!      O. Hellman et al., Phys Rev B 87, 104111 (2013)
!      O. Hellman and I. A. Abrikosov, Phys Rev B 88, 144301 (2013)
!      F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)
!
! Maintenance:
!===============================================================================
!> Interface to TDEP. 
!! Reads tdep unit cell and force constants into uc_tdep and fc_tdep structure
!! check agreement with DFPT grid
module lattice_tdep
   use kinds, only: dp
   use tdep_utils
   use qe_mpi_mod, only: stdout
   use cell_base, only: alat, tpiba
   use force_constant, only: lattice_ifc, force_const
   use brents_method, only: brent_helper
   use constants, only : pi, tpi
   use qe_mpi_mod, only: stdout
  
   implicit none
   public
   real(dp), parameter, private :: bohr2ang = 0.52917721067_dp
   real(dp), parameter, private :: hartree2ev = 27.21138602_dp
   real(dp), parameter, private :: sq_tol = 1.0E-10_dp
   
   type :: uc_tdep          !< lattice information for reading tdep inputs
      integer :: nelements                                       !< number of elements per unit cell
      integer :: na                                              !< number of atoms per unit cell
      integer, allocatable :: element_counter(:)                 !< counter for each element type
      character(len=8), allocatable :: atomic_symbol(:)          !< atomic symbol of each element
      real(dp), allocatable :: r(:,:)                            !< positions in fractional coordinates
      real(dp) :: latticevectors(3,3)                            !< lattice vectors
      real(dp) :: inv_latticevectors(3,3)                        !< inverse_lattice vectors
      real(dp) :: reciprocal_latticevectors(3,3)                 !< reciprocal lattice vectors
   end type

   type fc2_pair            !< single pair of atom
      integer :: i1                                !< index in the unit cell to atom 1
      integer :: i2                                !< index in the unit cell to atom 2
      real(dp) :: lv1(3)                           !< lattice vectors positioning the unit cell of atom1
      real(dp) :: lv2(3)                           !< lattice vectors positioning the unit cell of atom2
      real(dp) :: r(3)                             !< vector between atom 1 and 2 
      real(dp) :: m(3,3)                           !< force constants
   end type

   type fc2_atom           !< list of paris per atom
      integer :: n                                               !< number of pairs
      type(fc2_pair), allocatable :: pair(:)                     !< info about each pair
   end type
   
   type ewald_info         !< ewald information for polar materials
      real(dp) :: lambda                                         !< Ewald lambda parameter
      real(dp) :: eps(3,3)                                       !< dielectric tensor
      real(dp) :: inveps(3,3)                                    !< inverse dielectric tensor
      integer :: n_Rvector                                       !< number of R vectors to sum over
      integer :: n_Gvector                                       !< number of G vectors to sum over
      real(dp), allocatable :: Rvec(:,:)                         !< R vectors
      real(dp), allocatable :: Gvec(:,:)                         !< G vectors
   end type

   type fc_polarinfo       !< born charges and other information for polar materials
      integer :: correctiontype                                  !< correction type
      real(dp) :: eps(3,3)                                       !< dielectric tensor
      real(dp), allocatable :: born_effective_charges(:,:,:)     !< Born effective charges
      real(dp), allocatable :: born_onsite_correction(:,:,:)     !< Born onsite correction
      real(dp), allocatable :: Z_buff(:,:,:)                     !< Born onsite correction buffer
      integer :: nx_Z                                            !< number of irreducible Born effective charges
      real(dp), allocatable :: x_Z(:)                            !< irreducible representation of Born effective charges
      real(dp), allocatable :: coeff_Z(:,:)                      !< coefficient matrix to generate irreducible effective charge
      real(dp) :: lambda                                         !< Ewald lambda parameter
      logical :: symmborn                                      !< if true, do not enforce born charge hermiticity
   end type

   type :: fc_tdep         !< 2nd order forceconstants from tdep
      integer :: na                                              !< number of atoms per unit cell
      logical :: polar = .false.                                 !< is this a polar material
      type(fc_polarinfo) :: loto                                 !< born charges, dielectric tensor, etc.
      type(ewald_info) :: ew                                     !< polar dipole ewald information
      type(fc2_atom), allocatable :: atom(:)                     !< atoms
      real(dp) :: cutoff                                         !< length of longest pair + a tiny margin
   end type
      
   type(uc_tdep), save :: uc
   type(fc_tdep), save :: fc

   public :: init_ifc_tdep
   public :: uc_readfromfile
   public :: fc_readfromfile
contains
   !> Initialize lattice_ifc with data from tdep/forceconstant fc.
   subroutine init_ifc_tdep(uc, fc, ph, at, cryst_tau)
      implicit none
      type(uc_tdep), intent(in) :: uc        !< tdep unit cell structure
      type(fc_tdep), intent(in) :: fc        !< tdep force constant structure
      type(lattice_ifc), intent(inout) :: ph !< perturbo force constant structure 
      real(dp), intent(in) :: at(3,3) !< lattice vector in the unit of alat
      real(dp), intent(in) :: cryst_tau(:,:) !< atomic positions in crystal coordinates 
      integer :: idx, a1, a2, ir, i
      integer, allocatable :: ph2fc(:,:)
      type(force_const), pointer :: ptr
      !     
      !unit of ifc should be rydberg/bohr^2, to be consistent with QE.
      real(dp) :: fc_m(3, 3)
      !
      !the latest version of TDEP using Hatree atomic unit.
      !check "type_forceconstant_secondorder_io.f90/readfromfile"
      ! > fc%atom(a1)%pair(i)%m(j,:)=v1*lo_forceconstant_2nd_eVA_to_HartreeBohr
      ! 
      !conversion factor from Hartree/Borh^2 to Ryd/Bohr^2
      real(dp), parameter :: conv = 2.0_dp
      
      !make sure the crystal structure in TDEP and QE are consistent.
      call check_consistency_tdep(uc, fc, ph, at, cryst_tau)
      !make sure any pairs in infile.forceconstant cann be mapped to wscell
      allocate( ph2fc(ph%max_nr, ph%na*(ph%na+1)/2) )
      call check_wscell_coverage(uc, fc, ph, ph2fc)
      idx = 0
      do a2 = 1, ph%na
      do a1 = 1, a2
         ! idx = (a2 * (a2 - 1)) / 2 + a1;  a2 >= a1
         idx = idx + 1
         ptr => ph%phi(idx)

         do ir = 1, ptr%ws_ph%nr
            fc_m(:,:) = 0.0_dp
            ! map (ir, idx) in ph to (a1, i) in fc
            i = ph2fc(ir, idx)
            !the unit of ifc in tdep is Hartree/Borh^2. Converted to rydberg/bohr^2.
            if(i > 0) fc_m(:,:) = real(fc%atom(a1)%pair(i)%m(:,:), dp) * conv
            !assign fc_m to perturbo/lattice_ifc
            ptr%ifc(:,:,ir) = cmplx(fc_m, 0.0_dp, kind=dp)
         enddo
      enddo; enddo

      !now some polar parameter
      !check type_forceconstant_secondorder_io.f90/readfromfile
      !and type definition in type_forceconstant_secondorder.f90
      ph%lpol = fc%polar
      allocate( ph%pol%bcharge(3,3,ph%na) )
      !only polarcorrectiontype=3 is supported
      if(fc%polar .and. (fc%loto%correctiontype .ne. 3)) &
         call errore('init_ifc_tdep','Only support polarcorrectiontype=3 in TDEP', 1)
      
      if( ph%lpol ) then
         ph%pol%epsil(:,:) = real(fc%loto%eps(:,:), dp)
         ! i need to read lambda from file
         ph%pol%alpha = ( real(fc%loto%lambda, dp) / tpiba )**2
         !
         do a1 = 1, ph%na
            ph%pol%bcharge(:,:,a1) = real(fc%loto%born_effective_charges(:,:,a1), dp)
         enddo
      else
         ph%pol%epsil(:,:) = 1.0_dp
         ph%pol%alpha = 1.0_dp
         !
         ph%pol%bcharge(:,:,:) = 0.0_dp
      endif
      
      deallocate( ph2fc )
   end subroutine 


   !> Check consistency with phonon information from DFPT 
   subroutine check_consistency_tdep(uc, fc, ph, at, cryst_tau)
      implicit none
      type(uc_tdep), intent(in) :: uc        !< tdep unit cell structure
      type(fc_tdep), intent(in) :: fc        !< tdep force constant structure
      type(lattice_ifc), intent(inout) :: ph !< perturbo force constant structure 
      real(dp), intent(in) :: at(3,3) !< lattice vector in the unit of alat
      real(dp), intent(in) :: cryst_tau(:,:) !< atomic positions in crystal coordinates 
      integer :: i, a1
      real(dp) :: diff(3)
      character(len=120) :: trams, ctmp

      !print out basic structure info
      !Adapted from TDEP/type_crystalstructure.f90
      write(stdout,'(4x,a)') ' lattice vectors in TDEP calculations (in Bohr):'
      write(stdout,"(4X,A,3(2X,F18.12))") ' a1:',uc%latticevectors(:,1)
      write(stdout,"(4X,A,3(2X,F18.12))") ' a2:',uc%latticevectors(:,2)
      write(stdout,"(4X,A,3(2X,F18.12))") ' a3:',uc%latticevectors(:,3)
      write(stdout,'(4x,a)') ' ... reciprocal lattice vectors:'
      write(stdout,"(4X,A,3(2X,F18.12))") ' b1:',uc%reciprocal_latticevectors(:,1)
      write(stdout,"(4X,A,3(2X,F18.12))") ' b2:',uc%reciprocal_latticevectors(:,2)
      write(stdout,"(4X,A,3(2X,F18.12))") ' b3:',uc%reciprocal_latticevectors(:,3)
      trams=''
      do i = 1, uc%nelements
         write(ctmp,'(a,i2)') trim(uc%atomic_symbol(i))//':',uc%element_counter(i)
         trams=trim(trams)//'   '//trim(ctmp)
      enddo
      write(stdout,'(5x,a )') '... with composition   '//trim(adjustl(trams))
      write(stdout,'(5x,a/)') "----------------------------------------------"
      !
      if(fc%na.ne.ph%na) call errore("lattice_tdep","different number of atoms in TDEP and QE",1)
      !atomic position
      do a1 = 1, fc%na
         ! uc%r(:,a1): atomic position in crystal coordinate
         diff = real(uc%r(:,a1), dp) - cryst_tau(:,a1)
         if(norm2(diff) > 1.0E-6_dp) &
            call errore('lattice_tdep','different atomic position in TDEP and QE',1)
      enddo

      !lattice vector
      do i = 1, 3
         diff = real(uc%latticevectors(:,i), dp) - at(:,i)*alat
         if(norm2(diff) > 1.0E-6_dp*alat) &
            call errore('lattice_tdep','different lattice vectors in TDEP and QE',1)
      enddo
   end subroutine

   !> Check converage of Wigner Seitz Cell 
   subroutine check_wscell_coverage(uc, fc, ph, map_ph2fc)
      implicit none
      type(uc_tdep), intent(in) :: uc        !< tdep unit cell structure
      type(fc_tdep), intent(in) :: fc        !< tdep force constant structure
      type(lattice_ifc), intent(in) :: ph !< perturbo force constant structure 
      integer, intent(out) :: map_ph2fc(:,:) !< mapping map_ph2fc( ph%max_nr, na((na+1)/2) )
      !
      logical :: found
      integer :: a1, a2, i, idx, ir
      real(dp) :: rvec(3), lv2(3), diff(3)
      type(force_const), pointer :: ptr
      !
      real(dp), parameter :: r_tol = 1.0E-6_dp
      !init
      map_ph2fc(:,:) = 0
      do a1 = 1, fc%na
      do i = 1, fc%atom(a1)%n
         !this pair: a1, a2 + Rvec
         a2 = fc%atom(a1)%pair(i)%i2
         !Note that lv2 is in fractional in infile.forceconstant, but is converted to 
         ! Cartesian after readin via readfromfile 
         ! (see type_forceconstant_secondorder_io.f90/readfromfile)
         !here we convert it back to fractional for convenienice.
         lv2 = fc%atom(a1)%pair(i)%lv2
         call cryst_to_cart(1, lv2, uc%inv_latticevectors, 1)
         rvec = real(lv2, dp)
         !since (a1, a2+Rvec) and (a2, a1-Rvec) are equivalent, and only one of them
         ! are stored in Perturbo/lattice_ifc (a2 >= a1), we only check one case.
         if ( a1 > a2 ) cycle
         ! map (a1, a2) to the index used in lattice_ifc (see set_ws_cell_ph)
         ! e.g. (1, 1)->1, (1,2)->2, (2,2)->3, (1,3)->4, (2,3)->5, ...
         idx = a2 * (a2 - 1) / 2 + a1
         ptr => ph%phi(idx)
         
         found = .false.
         do ir = 1, ptr%ws_ph%nr
            !skip these already mapped (ir, idx)
            if( map_ph2fc(ir, idx) > 0 ) cycle
            diff(:) = ph%rvec_set(:, ptr%ws_ph%rvec(ir)) - rvec(:)
            ! if matches, label the matched (ir, idx), and exit the loop
            if( norm2(diff) < r_tol ) then
               found = .true.
               ! map (ir, idx) to (a1, i), and visa versa
               map_ph2fc(ir, idx) = i
               exit
            endif
         enddo
         !if not found, throw an error
         if(.not. found) call errore('check_wscell_coverage', &
            'Map tdep-ifc to perturbo failed. Try smaller rc2 in TDEP or larger q-grid in Perturbo',1)
      enddo
      enddo
   end subroutine

   !> Read unit cell information from TDEP outputs
   subroutine uc_readfromfile(filename, uc)
      character(len=*), intent(in) :: filename    !< infile.ucposcar file name
      type(uc_tdep), intent(out) :: uc            !< output uc_tdep
      !local variables
      logical :: seldyn                           !< if file is selective dynamics 
      logical :: alloy                            !< if material is an alloy (does not support)
      logical :: collmag                          !< if material is collinear magnetic
      logical :: noncollmag                       !< if material is noncollinear magnetic
      integer :: nelem                            !< element counter
      integer :: ios                              !< file unit
      integer :: u                                !< unit number for file
      integer, external :: find_free_unit         !< finding a free unit to open file

      real(dp) :: m(3,3)                          !< lattice vector m(:,i) is the ith lattice vector
      character(len=1000), allocatable :: symbols(:) !< element symbols
      real(dp), allocatable :: r(:,:)             !< position vectors      
      real(dp), allocatable :: noncollmagmom(:,:) !< non collinear magnetic moments
      real(dp), allocatable :: collmagmom(:)      !< collinear magnetic moments
      integer, allocatable :: atomic_number(:)    !< atomic number
      logical, allocatable :: cmatom(:)           !< collinear magnetism
      real(dp) :: im(3,3)                         !< inverse of m

      integer :: n_elem                           !< number of element counter
      integer :: na                               !< number of atoms per unit cell
      integer :: ia                               !< loop variable
      integer :: ja                               !< loop variable

      dummyopen: block
         integer, parameter :: max_n_elements=20000
         integer :: i,j,k
         character(len=10*max_n_elements) :: trams
         u = find_free_unit()
         open(u, file=trim(filename), status='old', action='read', err=101, iostat=ios)
         read(u,*) trams
         read(u,*) trams
         read(u,*) trams
         read(u,*) trams
         read(u,*) trams
         read(u,*) trams
         read(u,'(A)') trams
         j=0
         i=0
         ! figure out how many specified elements there are, this is surprisingly robust.
         do
            if ( j > max_n_elements ) exit
            j=j+1
            if ( trams(j:j) .ne. ' ' ) then
               i=i+1
               do k=1,max_n_elements
                  if ( trams(j+1:j+1) .eq. ' ' ) exit
                  if ( trams(j+1:j+1) .ne. ' ' ) j=j+1
                  if ( k .eq. max_n_elements ) then
                     call errore('lattice_tdep', ' different elements in the structure? really?', 1)
                     stop
                  endif
               enddo
            endif
         enddo
         n_elem=i
         ! perhaps it is a POSCAR with selective dynamics
         ! I will ignore that, but it's nice if it does not crash
         read(u,*) trams
         seldyn=.false.
         if ( trams(1:1) .eq. 's' .or. trams(1:1) .eq. 'S' ) seldyn=.true.
         close(u)
      end block dummyopen

      readstuff: block
         real(dp) :: latpar_or_volume         !< factor for lattice vector
         integer, allocatable :: elemcount(:) !< element counter   
         character(len=1) :: trams            !< reading buffer
         character(len=10) :: dumsp           !< DLM atom
         character(len=1) :: coord_type       !< coordinate type (Cartesian or Direct)
         logical :: cartesian                 !< if coordinate = Cartesian
         integer :: i                         !< loop index
         integer :: j                         !< loop index
         integer :: l                         !< loop index
         real(dp) :: f0                       !< volumn inverse
         real(dp) :: v(3)                     !< temporary vector
         allocate(symbols(n_elem))
         allocate(elemcount(n_elem))
         
         u = find_free_unit()
         open(u, file=trim(filename), status='old', action='read', err=101, iostat=ios)
         read(u, *)    ! header 
         read(u, *) latpar_or_volume   !< read latpar to volume
         read(u, *) m(:,1)
         read(u, *) m(:,2)
         read(u, *) m(:,3)
         read(u, *) symbols
         read(u, *) elemcount
         read(u, *) trams
         ! figure out cartesian or fractional coordiates
         if ((trams .eq. 'C') .or. (trams .eq. 'c')) then
            cartesian = .true.
         elseif ((trams .eq. 'D') .or. (trams .eq. 'd')) then
            cartesian = .false.
         else
            call errore('lattice_tdep', 'specify either cartesian or direct coordinates in infile.ucposcar', 1)
         endif
         ! skip a line for self dynamics
         if (seldyn) read(u,*) trams 
         alloy = .false.
         do i = 1, n_elem
            if ( symbols(i)(1:5) .eq. 'alloy' .or. symbols(i)(1:5) .eq. 'Alloy' .or. symbols(i)(1:5) .eq. 'ALLOY' ) then
               symbols(i) = 'ALLOY'
               alloy = .true.
            endif
         enddo
         collmag = .false.
         do i = 1, n_elem
            if ( symbols(i)(1:2) .eq. 'CM')  then
               symbols(i) = 'CM'
               collmag =.true.
            endif
         enddo
         noncollmag = .false.
         do i = 1, n_elem
            if ( symbols(i)(1:3) .eq. 'NCM')  then
               symbols(i) = 'NCM'
               noncollmag =.true.
            endif
         enddo

         ! do not treat alloy in qe2pert
         if (alloy) then
            call errore('lattice_tdep', 'perturbo does not support alloy right now', 1)
         endif
         ! no magnetism in qe2pert
         if (collmag .or. noncollmag) then
            write(*,*) 'ignore magnetism for forceconstant reading'
         endif
         
         if (latpar_or_volume .gt. 0.0_dp) then
            m = m * latpar_or_volume
         else
            f0 = (abs(latpar_or_volume) / abs(determ_real(m))) ** (1.0_dp/3.0_dp)
            m = m * f0
         endif
         ! find inverse of m which is transpose of reciprocal lattice vector
         im = invert3x3matrix(m)

         ! allocate temporary space to save variables
         na = sum(elemcount)
         allocate(r(3, na))
         allocate(atomic_number(na))
         allocate(cmatom(na))
         allocate(collmagmom(na))
         allocate(noncollmagmom(3,na))
         r = 0.0_dp
         atomic_number = 0
         cmatom = .false.
         collmagmom = 0.0_dp
         noncollmagmom = 0.0_dp
         
         if (collmag .or. noncollmag) then
            cmatom = .false.
            collmagmom = 0.0_dp
            noncollmagmom = 0.0_dp
            l = 0
            do i = 1, n_elem
               do j = 1, elemcount(i)
                  l = l + 1
                  if ((collmag) .and. (trim(symbols(i)) .eq. 'CM')) then
                     ! read DLM atom
                     read(u, *) v, dumsp, collmagmom(l)
                     cmatom(l) = .true.
                     atomic_number(l) = symbol_to_z(trim(dumsp))
                  else if ((noncollmag) .and. (trim(symbols(i)) .eq. 'NCM')) then
                     ! read DLM atom
                     read(u, *) v, dumsp, noncollmagmom(:,l)
                     cmatom(l) = .true.
                     atomic_number(l) = symbol_to_z(trim(dumsp))
                  else
                     read(u, *) v
                     atomic_number(l) = symbol_to_z(trim( adjustl(symbols(i))))
                  endif
                  if (cartesian) then
                     ! transform v into crystal coodinates
                     ! find bg which is m^{-T}
                     call cryst_to_cart(1, v, m, -1)
                  endif
                  r(:,l) = clean_fractional_coordinates(v, sq_tol)
               enddo
            enddo
         else ! non magnetic
            l = 0
            do i = 1, n_elem
               do j = 1, elemcount(i)
                  l = l + 1
                  atomic_number(l) = symbol_to_z(trim( adjustl(symbols(i))))
               enddo
            enddo
            do i = 1, na
               read(u,*) v
               if (cartesian) then
                  call cryst_to_cart(1, v, m, -1)
               endif
               r(:,i) = clean_fractional_coordinates(v, sq_tol)
            enddo
         endif
         close(u)
         deallocate(elemcount)
      end block readstuff
 
      ! put all parameters into uc_tdep 
      uc%na = na
      ! convert position from Angstrom to bohr radius
      uc%latticevectors = m / bohr2ang
      uc%inv_latticevectors = chop_real(invert3x3matrix(uc%latticevectors), sq_tol)
      uc%latticevectors = chop_real(uc%latticevectors, sq_tol)
      uc%reciprocal_latticevectors = chop_real(reciprocal_basis(uc%latticevectors), sq_tol)
      allocate(uc%r(3, uc%na))
      uc%r = 0.0_dp
      do ia = 1, uc%na
         uc%r(:,ia) = chop_real(clean_fractional_coordinates(r(:,ia), sq_tol), sq_tol)
      enddo
      ! update nelements, element counter and atomic symbol
      fixatoms: block
         integer, allocatable :: dum(:)
         call return_unique_integers(atomic_number, dum)
         uc%nelements = size(dum, 1)
         allocate(uc%atomic_symbol(uc%nelements)) 
         allocate(uc%element_counter(uc%nelements)) 
         uc%element_counter = 0
         do ia = 1, uc%na
            do ja = 1, uc%nelements
               if ( dum(ja) .eq. atomic_number(ia)) then
                  uc%element_counter(ja) = uc%element_counter(ja) + 1
               endif
            enddo
         enddo
         do ja = 1, uc%nelements
            uc%atomic_symbol(ja) = trim(adjustl(z_to_symbol( dum(ja))))
         enddo
      end block fixatoms
   
      deallocate(symbols)
      deallocate(r)
      deallocate(atomic_number)
      deallocate(cmatom)
      deallocate(collmagmom)
      deallocate(noncollmagmom)
         
101 call errore('uc_readfromfile', 'open '//trim(filename), abs(ios))
   end subroutine uc_readfromfile
   
   !> Read 2nd order force constant information from TDEP outputs
   subroutine fc_readfromfile(filename, uc, fc)
      use input_param, only : symmborn
      character(len=*), intent(in) :: filename    !< infile.ucposcar file name
      type(uc_tdep), intent(inout) :: uc          !< in and output unit cell info uc_tdep
      type(fc_tdep), intent(out) :: fc            !< output forceconstant fc_tdep
     
      integer :: correction 
      integer :: i, j, a1, a2, u
      integer :: ios
      real(dp) :: v1(3), v2(3), lv1(3), lv2(3), ucv1(3), ucv2(3), r(3)
      real(dp) :: lambda
      real(dp), parameter :: ewaldtol = 1.0E-20_dp
      integer, external :: find_free_unit

      real(dp), dimension(:), allocatable :: Zflat
      complex(dp), dimension(:, :, :, :), allocatable :: Dc
      real(dp), dimension(3, 3) :: acc_chg
      ! set everything to nothing
      fc%na = -1
      fc%cutoff = -1
      fc%polar = .false.
      fc%loto%correctiontype = -1
      fc%loto%eps = -1.0_dp
      fc%loto%nx_Z = -1

      ! read file
      u = find_free_unit()
      open(u, file=trim(filename), status='old', action='read', err=102, iostat=ios)
      read(u, *) fc%na
      if (fc%na .ne. uc%na) then
         call errore('fc_readfromfile', 'number of atoms mismatch between two files', 1)
      endif
      read(u, *) fc%cutoff
      allocate(fc%atom(fc%na))
      do a1 = 1, fc%na
         read(u, *) fc%atom(a1)%n
         allocate( fc%atom(a1)%pair(fc%atom(a1)%n) )
         do i = 1, fc%atom(a1)%n
            fc%atom(a1)%pair(i)%i1 = -1 
            fc%atom(a1)%pair(i)%i2 = -1 
            fc%atom(a1)%pair(i)%lv1 = -1
            fc%atom(a1)%pair(i)%lv2 = -1
            fc%atom(a1)%pair(i)%r = -1
            fc%atom(a1)%pair(i)%m = -1
         enddo

         do i = 1, fc%atom(a1)%n
            fc%atom(a1)%pair(i)%i1 = a1
            read(u, *) fc%atom(a1)%pair(i)%i2
            read(u, *) fc%atom(a1)%pair(i)%lv2
            do j = 1, 3
               read(u, *) v1
               fc%atom(a1)%pair(i)%m(j, :) = v1 * bohr2ang**2 /hartree2ev
            enddo
         enddo
      enddo
      ! check if polar
      read(u, *) correction
      if (correction .eq. 0) then
         ! not polar
         fc%polar = .false.
         fc%loto%correctiontype = 0
      else
         ! polar
         fc%polar = .true.
         fc%loto%correctiontype = correction
         do i = 1, 3
           read(u, *) fc%loto%eps(:, i)
         enddo
         read(u, *) lambda
         read(u, *) fc%loto%nx_Z
         allocate( fc%loto%x_Z(fc%loto%nx_Z))  
         allocate( fc%loto%coeff_Z(fc%na*9, fc%loto%nx_Z))  
         fc%loto%x_Z = 0.0_dp
         fc%loto%coeff_Z = 0.0_dp
         do i = 1, fc%loto%nx_Z
            read(u, *) fc%loto%x_Z(i)
         enddo
         do i = 1, fc%na*9
            read(u, *) fc%loto%coeff_Z(i, :)
         enddo

         allocate( fc%loto%born_effective_charges(3, 3, fc%na))
         allocate( fc%loto%born_onsite_correction(3, 3, fc%na))
         fc%loto%born_effective_charges = 0.0_dp
         fc%loto%born_onsite_correction = 0.0_dp

         ! read raw Born effective charges
         allocate( fc%loto%Z_buff(3, 3, fc%na))
         do a1 = 1, fc%na
            do i = 1, 3
               read(u, *) fc%loto%Z_buff(1:3, i, a1)
            enddo
         enddo
      endif 
      close(u)

      ! convert from fractional to cartesian and calculate actual cutoff in cartesian coordinates
      fc%cutoff = 0.0_dp
      do a1 = 1, fc%na
         do i = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%pair(i)%i2
            ! get all the vectors right 
            ucv1 = uc%r(:, a1)
            ucv2 = uc%r(:, a2)
            lv1 = 0.0_dp
            lv2 = fc%atom(a1)%pair(i)%lv2
            v1 = lv1 + ucv1
            v2 = lv2 + ucv2
            r = v2 - v1
            ! convert to cartesian coordinates
            call cryst_to_cart(1, lv1, uc%latticevectors, 1)
            call cryst_to_cart(1, lv2, uc%latticevectors, 1)
            call cryst_to_cart(1, r, uc%latticevectors, 1)
            fc%atom(a1)%pair(i)%lv1 = chop_real(lv1, 1.0E-13_dp)
            fc%atom(a1)%pair(i)%lv2 = chop_real(lv2, 1.0E-13_dp)
            fc%atom(a1)%pair(i)%r = chop_real(r, 1.0E-13_dp)
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%pair(i)%r))
         enddo
      enddo
      fc%cutoff = fc%cutoff + sq_tol
      ! set ewald parameters
      if (fc%polar) then
         ! check correction type
         select case (fc%loto%correctiontype)
         case (0)
            ! fine
         case (1:2)
            ! not fine
            call errore('fc_readfromfile', 'This polar correction type is deprecated', 1)
         case (3)
            ! also fine
         case default
            ! This should never happen
            call errore('fc_readfromfile', 'Undefined polar correction type', 1)
         end select 

         fc%loto%lambda = lambda 
         ! set ewald parameters
         call set_dipole_ewald_parameters(uc, fc%ew, fc%loto%eps, ewaldtol, lambda)

         ! set forces to be hermitian
         ! if symmborn enforce hermiticity - jjzhou
         if ( symmborn ) then
            ! enforce hermiticity
            write(stdout, *) 'Enforcing Hermiticity on Born effective charges'
            fc%loto%symmborn = .true.
         else
            write(stdout, *) 'Skipping Hermiticity enforcement on Born effective charges'
            ! do not enforce hermiticity
            fc%loto%symmborn = .false. 
         endif

         if (fc%loto%symmborn) then
   
            ! enforce hermiticity on irreducible charges
            call force_borncharges_Hermitian(fc%ew, fc%loto%x_Z, fc%loto%coeff_Z, fc%loto%eps, uc)

            ! Now set the charges and so on.

            ! Convert from flat to sensible charges
            allocate(Zflat(9*fc%na))
            Zflat = 0.0_dp
            call dgemv('N', 9*uc%na, fc%loto%nx_Z, 1.0_dp, fc%loto%coeff_Z, 9*uc%na, fc%loto%x_z, 1, 0.0_dp, Zflat, 1)
            acc_chg = 0.0_dp
            do a1 = 1, uc%na
               fc%loto%born_effective_charges(:, :, a1) = unflatten_2tensor(Zflat((a1 - 1)*9 + 1:a1*9))
               acc_chg = acc_chg + fc%loto%born_effective_charges(:, :, a1)
            enddo
            if (sum(abs(acc_chg)) .gt. sq_tol) then
               call errore('fc_readfromfile', 'Born charges do not add up to zero', 1)
            endif
            deallocate(Zflat)

         else
            fc%loto%born_effective_charges(:,:,:) = fc%loto%Z_buff(:,:,:)
         endif

         ! now add the charges for both cases
         allocate(Dc(3, 3, uc%na, uc%na))
         Dc = 0.0_dp
         select case (fc%loto%correctiontype)
         case (3)
            call longrange_dynamical_matrix(fc%ew, uc, [0.0_dp, 0.0_dp, 0.0_dp], fc%loto%born_effective_charges, fc%loto%born_onsite_correction, fc%loto%eps, Dc, reconly=.true.)
         case default
            call errore('fc_readfromfile', 'Undefined polar correction type', 1)
         end select
         ! Slight sanity check
         if (sum(abs(aimag(Dc))) .gt. sq_tol*fc%na**2) then
            call errore('fc_readfromfile', 'Dynamical matrix at gamma point is not real', 1)
         endif

         ! Now I can determine the onsite correction term!
         do a1 = 1, fc%na
            acc_chg = 0.0_dp
            do a2 = 1, fc%na
               acc_chg = acc_chg + real(Dc(:, :, a1, a2), dp)
            enddo
            fc%loto%born_onsite_correction(:, :, a1) = -acc_chg
         enddo
         fc%loto%born_onsite_correction = chop_real(fc%loto%born_onsite_correction, 1.0E-13_dp)
         deallocate(Dc)

      endif ! if polar

102 call errore('fc_readfromfile', 'open '//trim(filename), abs(ios))
   end subroutine fc_readfromfile
   
   !> pick optimum Ewald parameters for dipole terms
   subroutine set_dipole_ewald_parameters(uc, ew, eps, tol, forced_lambda)
      type(uc_tdep), intent(in) :: uc                 !< crystal structure
      type(ewald_info), intent(out) :: ew             !< ewald parameters
      real(dp), dimension(3, 3), intent(in) :: eps    !< dielectric tensor
      real(dp), intent(in) :: tol                     !< tolerance
      real(dp), intent(in), optional :: forced_lambda !< force a certain lambda parameter
      integer, parameter :: npts = 40
      real(dp), dimension(3, npts) :: pts
       
      call points_on_sphere(pts)
      ! store the dielectric tensor, and the inverse
      ew%eps = eps
      ew%inveps = invert3x3matrix(eps)
      ew%lambda = forced_lambda
      ! strategy = 3 enforced

      ! now build the vectors
      buildvecs: block
         integer, parameter :: niter = 1000
         integer, dimension(3) :: bd
         real(dp), dimension(:, :), allocatable :: dr
         real(dp), dimension(:), allocatable :: vn
         real(dp), dimension(3, 3) :: m0
         real(dp), dimension(3) :: v0
         real(dp) :: krad, rrad, f0, f1
         integer, dimension(:), allocatable :: di
         integer :: i, j, k, l
   
         ! Get the radii in real and reciprocal space
         call ewald_dipole_k_r(ew%lambda, pts, ew%eps, ew%inveps, uc, tol, rrad, krad)

         ! Start with realspace, get a supercell large enough to fit the cutoff
         bd = 1
         do
            do i = 1, 3
               m0(:, i) = uc%latticevectors(:, i)*(2*bd(i) + 1)
            enddo
            f0 = inscribed_sphere_in_box(m0)
            if (f0 .gt. rrad) exit
            bd = increment_dimensions(bd, uc%latticevectors)
         enddo
         ! count realspace vectors
         f0 = rrad**2
         l = 0
         do i = -bd(1), bd(1)
         do j = -bd(2), bd(2)
         do k = -bd(3), bd(3)
            v0 = real([i, j, k], dp)
            call cryst_to_cart(1, v0, uc%latticevectors, 1)
            if (sqnorm(v0) .lt. f0) then
               l = l + 1
            endif
         enddo
         enddo
         enddo
         ew%n_Rvector = l
         allocate(ew%Rvec(3, l))
         allocate(dr(3, l))
         allocate(vn(l))
         allocate(di(l))
         ew%rvec = 0.0_dp
         dr = 0.0_dp
         vn = 0.0_dp
         di = 0
         l = 0
         do i = -bd(1), bd(1)
         do j = -bd(2), bd(2)
         do k = -bd(3), bd(3)
            v0 = real([i, j, k], dp)
            call cryst_to_cart(1, v0, uc%latticevectors, 1)
            f1 = sqnorm(v0)
            if (f1 .lt. f0) then
               l = l + 1
               vn(l) = -f1
               dr(:, l) = v0
            endif
         enddo
         enddo
         enddo

         ! Sort backwards by length
         call quick_sort_real(vn, di)
         do i = 1, ew%n_Rvector
             ew%Rvec(:, i) = dr(:, di(i))
         enddo
         deallocate(dr)
         deallocate(vn)
         deallocate(di)

         ! Then exactly the same thing again for G-vectors
         bd = 1
         do
            do i = 1, 3
               m0(:, i) = uc%reciprocal_latticevectors(:, i)*(2*bd(i) + 1)
            enddo
            f0 = inscribed_sphere_in_box(m0)
            if (f0 .gt. krad) exit
            bd = increment_dimensions(bd, uc%reciprocal_latticevectors)
         enddo
         ! count realspace vectors
         f0 = krad**2
         l = 0
         do i = -bd(1), bd(1)
         do j = -bd(2), bd(2)
         do k = -bd(3), bd(3)
            v0 = real([i, j, k], dp)
            call cryst_to_cart(1, v0, uc%reciprocal_latticevectors, 1)
            if (sqnorm(v0) .lt. f0) then
               l = l + 1
            endif
         enddo
         enddo
         enddo

         ew%n_Gvector = l
         allocate(ew%Gvec(3, l))
         allocate(dr(3, l))
         allocate(vn(l))
         allocate(di(l))
         ew%Gvec = 0.0_dp
         dr = 0.0_dp
         vn = 0.0_dp
         di = 0
         l = 0
         do i = -bd(1), bd(1)
         do j = -bd(2), bd(2)
         do k = -bd(3), bd(3)
            v0 = real([i, j, k], dp)
            call cryst_to_cart(1, v0, uc%reciprocal_latticevectors, 1)
            f1 = sqnorm(v0)
            if (f1 .lt. f0) then
               l = l + 1
               vn(l) = -f1
               dr(:, l) = v0
            endif
         enddo
         enddo
         enddo
         ! Sort backwards by length
         call quick_sort_real(vn, di)
         do i = 1, ew%n_Gvector
             ew%Gvec(:, i) = dr(:, di(i))
         enddo
         deallocate(dr)
         deallocate(vn)
         deallocate(di)
      end block buildvecs
   
   end subroutine
   
   !> turn the ratio of k-space volume to r-space volume into a function to be minimized
   subroutine ewald_dipole_k_r(lambda, pts, eps, inveps, uc, tol, rrad, krad)
      real(dp), intent(in) :: lambda !< ewald parameter
      real(dp), dimension(:, :), intent(in) :: pts !< points to test with
      real(dp), dimension(3, 3), intent(in) :: eps !< dielectric tensor
      real(dp), dimension(3, 3), intent(in) :: inveps !< inverse dielectric tensor
      type(uc_tdep), intent(in) :: uc !< crystal structure
      real(dp), intent(in) :: tol !< tolerance to target
      real(dp), intent(out) :: rrad !< radius in r-space
      real(dp), intent(out) :: krad !< radius in k-space
   
      integer, parameter :: niter = 1000
      type(brent_helper) :: bt
   
      ! Init to something
      krad = -1.0_dp
      rrad = -1.0_dp
   
      ! find radius in K-space
      findkrad: block
         real(dp), dimension(3) :: v0, v1
         real(dp) :: r0, r1, y0, y1, f0, f1, x, y, inv4lambdasq
         integer :: iter, i, j, k
   
         f0 = inscribed_sphere_in_box(uc%reciprocal_latticevectors)
         inv4lambdasq = 1.0_dp/(4.0_dp*lambda*lambda)
         ! starting guesses
         r0 = f0*1.0E-9_dp
         r1 = f0*1.0E9_dp
         y0 = 0.0_dp
         y1 = 0.0_dp
         do i = 1, size(pts, 2)
            v0 = pts(:, i)*r0*tpi
            v1 = pts(:, i)*r1*tpi
            f0 = 0.0_dp
            f1 = 0.0_dp
            do j = 1, 3
            do k = 1, 3
               f0 = f0 + v0(k)*eps(k, j)*v0(j)
               f1 = f1 + v1(k)*eps(k, j)*v1(j)
            enddo
            enddo
            y0 = y0 + exp(-f0*inv4lambdasq)*f0
            y1 = y1 + exp(-f1*inv4lambdasq)*f1
         enddo
         y0 = y0 - tol
         y1 = y1 - tol
         ! Make sure we bracket the target
         if (y0*y1 .gt. 0.0_dp) then
            call errore('ewald_dipole_k_r', 'Could not bracket values', 1)
         endif
   
         call bt%init(r0, r1, y0, y1, tol*1.0E-3_dp)
         do iter = 1, niter
            ! check for convergence
            if (abs(bt%current_y()) .lt. tol*1.0E-3_dp) then
               ! found a zero
               krad = bt%current_x()
               exit
            endif
            ! pick next x-value
            x = bt%next_x()
            ! evaluate function
            y = 0.0_dp
            do i = 1, size(pts, 2)
               v0 = pts(:, i)*x*tpi
               f0 = 0.0_dp
               do j = 1, 3
               do k = 1, 3
                  f0 = f0 + v0(k)*eps(k, j)*v0(j)
               enddo
               enddo
               y = y + exp(-f0*inv4lambdasq)*f0
            enddo
            y = y - tol
            ! tell the solver what happened
            call bt%update(x, y)
            if (iter .eq. niter) then
               call errore('ewald_dipole_k_r', 'Could not converge Brents method', 1)
            endif
         enddo
         ! then add the bounding sphere to be on the safe side, the longest possible q-vector added later.
         krad = krad + bounding_sphere_of_box(uc%reciprocal_latticevectors)
      end block findkrad
   
      ! Get the radius in R-space, same thing
      findrrad: block
         real(dp), dimension(3) :: v0, v1
         real(dp) :: r0, r1, y0, y1, x, y, f0, f1, invdete
         integer :: iter, i
   
         invdete = 1.0_dp/sqrt(determ_real(eps))
         f0 = inscribed_sphere_in_box(uc%latticevectors)
         ! This should bracket any sensible choice
         r0 = f0*1.0E-4_dp
         r1 = f0*1.0E4_dp
         ! Get the function values
         y0 = 0.0_dp
         y1 = 0.0_dp
         do i = 1, size(pts, 2)
            v0 = matmul(inveps, pts(:, i)*r0)
            v1 = matmul(inveps, pts(:, i)*r1)
            f0 = sqrt(dot_product(v0, pts(:, i)*r0))
            f1 = sqrt(dot_product(v1, pts(:, i)*r1))
            y0 = y0 + (lambda**3)*sum(abs(anisotropic_dipole_int(v0*lambda, f0*lambda, inveps)))*invdete
            y1 = y1 + (lambda**3)*sum(abs(anisotropic_dipole_int(v1*lambda, f1*lambda, inveps)))*invdete
         enddo
         y0 = y0 - tol
         y1 = y1 - tol
         ! Be sure we bracket the target
         if (y0*y1 .gt. 0.0_dp) then
             call errore('ewald_dipole_k_r', 'Could not converge braket values', 1)
         endif
         call bt%init(r0, r1, y0, y1, tol*1.0E-3_dp)
         do iter = 1, niter
             ! check for convergence
             if (abs(bt%current_y()) .lt. tol*1.0E-3_dp) then
                 ! found a zero
                 rrad = bt%current_x()
                 exit
             endif
             ! pick next x-value
             x = bt%next_x()
             ! evaluate function
             y = 0.0_dp
             do i = 1, size(pts, 2)
                 v0 = matmul(inveps, pts(:, i)*x)
                 f0 = sqrt(dot_product(v0, pts(:, i)*x))
                 y = y + (lambda**3)*sum(abs(anisotropic_dipole_int(v0*lambda, f0*lambda, inveps)))*invdete
             enddo
             y = y - tol
             ! tell the solver what happened
             call bt%update(x, y)
             if (iter .eq. niter) then
                 call errore('ewald_dipole_k_r', 'Could not converge Brents method', 1)
             endif
         enddo
         ! then add the bounding sphere, the longest possible rij-vector added later.
         rrad = rrad + bounding_sphere_of_box(uc%latticevectors)
      end block findrrad
   end subroutine


   !> force Born charges to become Hermitian. Serial for now.
   subroutine force_borncharges_Hermitian(ew, x_Z, coeff_Z, eps, uc)
      type(ewald_info), intent(in) :: ew  !< ewald information 
      real(dp), dimension(:), intent(inout) :: x_Z !< irreducible representation of Born charges
      real(dp), dimension(:, :), intent(in) :: coeff_Z !< coefficient matrix
      !! that can transform the irreducible representation to actual Born charges.
      real(dp), dimension(3, 3), intent(in) :: eps !< dielectric constant
      type(uc_tdep), intent(in) :: uc !< crystal structure
   
      ! Define some magical parameters for the nonlinear solver
      real(dp), parameter :: errtol = 1.0E-28_dp
      integer, parameter :: derivorder = 3
   
      ! Normal temporary things
      real(dp), dimension(:, :, :, :), allocatable :: wD, dynmat_reciprocal
      real(dp), dimension(:, :, :), allocatable :: wZ
      real(dp), dimension(:), allocatable :: wf
      real(dp) :: derivdelta
      integer :: nx
   
      ! Set some sensible things
      init: block
         complex(dp), dimension(:, :, :, :), allocatable :: D0
         real(dp), dimension(:, :, :), allocatable :: w0
         integer :: a1
         ! Number of variables to optimize for
         nx = size(x_Z)
         ! Workspace for born charges, tensor form
         allocate(wZ(3, 3, uc%na))
         wZ = 0.0_dp
         ! Workspace for born charges, flattened
         allocate(wf(9*uc%na))
         wf = 0.0_dp
         ! Unflatten the initial born charges
         !@TODO think about how annoying this is to change such that it is automagical?
         call dgemv('N', 9*uc%na, nx, 1.0_dp, coeff_Z, 9*uc%na, x_z, 1, 0.0_dp, wf, 1)
         do a1 = 1, uc%na
             wZ(:, :, a1) = unflatten_2tensor(wf((a1 - 1)*9 + 1:a1*9))
         enddo
         ! Get the dynamical matrix at Gamma, but without any charges multiplied in.
         ! Also create some workspace for the error and gradient functions.
         allocate(D0(3, 3, uc%na, uc%na))
         allocate(wD(3, 3, uc%na, uc%na))
         allocate(dynmat_reciprocal(3, 3, uc%na, uc%na))
         allocate(w0(3, 3, uc%na))
         D0 = 0.0_dp
         wD = 0.0_dp
         w0 = 0.0_dp
         dynmat_reciprocal = 0.0_dp
         call longrange_dynamical_matrix(ew, uc, [0.0_dp, 0.0_dp, 0.0_dp], wZ, w0, eps, D0, reconly=.true., chgmult=.false.)
         ! Make sure it's real. Has to be.
         if (sum(abs(aimag(D0))) .gt. sq_tol) then
            call errore('force_borncharges_hermitian', 'Dipole dynamical matrix at Gamma is not real', 1) 
         else
             ! keep only the real part to work with, slightly faster.
             dynmat_reciprocal = real(D0, dp)
         endif
         ! cleanup
         deallocate(D0)
         deallocate(w0)
      end block init
   
      ! Figure out settings for the nonlinear solver. Also possibly kill early, in
      ! case things are already fine.
      preopt: block
         real(dp) :: err0
   
         ! Calculate the current error?
         err0 = Z_error_function(coeff_Z, x_Z, wf, wZ, wD, dynmat_reciprocal)
         ! We could be fine?
         if (err0 .lt. errtol) then
            write(stdout, *) '... Born charges already Hermitian'
            return
         endif
   
         ! We were not fine. Too bad.
         write(stdout, *) '... Born charges not Hermitian, error:', err0
         write(stdout, *) '... DISPERSIONS MIGHT BE IRREGULAR!'
         write(stdout, *) '... FIXME!'
   
         ! fkdev: we currently don't enfore anything until we now how to do this properly
         return
      end block preopt
   
   
      ! Apparently we have to optimize for reals. I hate this. 50% of the time, it works 100% of the time.
      opt: block
          ! First some parameters for the minimization
          integer, parameter :: nouter = 1000
          integer, parameter :: ninner = 5000
          integer, parameter :: maxctr = 25
          real(dp), parameter :: gradtol = 1.0E-20_dp
          real(dp), parameter :: longtime = 1.0_dp ! if I waited this long, print how it's going.
          ! Things needed to set stuff
          type(nonlinear_minimization_helper) :: mh
          real(dp), dimension(:), allocatable :: x, x1, x2, grad, searchdir
          real(dp) :: err0, err1, err2, step
          real(dp) :: t0, t1
          integer :: outiter, initer, ctr, ictr
          integer :: idirection
   
          ! Some workspace
          allocate(x(nx))
          allocate(x1(nx))
          allocate(x2(nx))
          allocate(grad(nx))
          allocate(searchdir(nx))
          x = 0.0_dp
          x1 = 0.0_dp
          x2 = 0.0_dp
          grad = 0.0_dp
          searchdir = 0.0_dp
          ! Start guessing things for the solver
          derivdelta = 1.0E-6_dp  ! or some other random number
          x = x_Z               ! starting point
          ! Calculate error and gradient
          err0 = Z_error_function(coeff_Z, x, wf, wZ, wD, dynmat_reciprocal)
          !derivdelta=mh%get_delta(err0,1.0_dp)
          ! Just something to guess with.
          derivdelta = 1.0E-2_dp
          call Z_error_gradient(coeff_Z, x, derivdelta, grad, wf, wZ, wD, dynmat_reciprocal)
          ! Initialize solver
          call init_nonlinear_helper(mh, nx, x, grad, initial_hessian=10.0_dp)
   
          write(stdout, *) '... optimizing Born charges to remove error'
          write(stdout, "(5X,A10,3X,A14,3X,A14)") 'iter', 'err', 'grad'
          ictr = 0
          outloop: do outiter = 1, nouter
             ! Get the current error
             err0 = Z_error_function(coeff_Z, x, wf, wZ, wD, dynmat_reciprocal)
             ! Update the derivative, if we have enough updates of the Hessian
             derivdelta = min(derivdelta, get_delta(mh, err0, 1.0E4_dp))
             call Z_error_gradient(coeff_Z, x, derivdelta, grad, wf, wZ, wD, dynmat_reciprocal)
   
             ! Check if we are converged
             if (err0 .lt. errtol) then
                ! Converged
                exit outloop
             endif
             if (norm2(grad) .lt. gradtol) then
                 ! Converged, I think
                 exit outloop
             endif
             ! Get the search direction
             call get_searchdir(mh, x, grad, err0)
   
             ! Do a line search in this direction?
             do idirection = 1, mh%n_direction
                 x1 = mh%x_curr
                 err0 = mh%e_curr
                 err2 = huge(1.0_dp)
                 ctr = 0
                 !step=1.0_dp
                 searchdir = mh%sdir(:, idirection)
                 searchdir = searchdir/norm2(searchdir)
                 ! pick a sensible step size?
                 step = get_step(mh, searchdir, grad, err0) !/1.5_dp**3
                 inloop: do initer = 1, ninner
                    ictr = ictr + 1
                    ! take a step
                    x1 = x1 + step*searchdir
                    ! evaluate error
                    err1 = Z_error_function(coeff_Z, x1, wf, wZ, wD, dynmat_reciprocal)
                    if (err1 .lt. errtol) then
                        ! Could happen that we are really converged here!
                        x = x1
                        err0 = err1
                        exit outloop
                    elseif (err1 .lt. err0) then
                        ! Good step, keep going
                        step = step*1.5_dp
                        ! If it was better than the previous best, keep it!
                        if (err1 .lt. err2) then
                            err2 = err1
                            x2 = x1
                        endif
                    else
                        ! Bad step, backstep and decrease the step size
                        x1 = x1 - step*searchdir
                        step = step*0.25_dp
                        ctr = ctr + 1
                    endif
                    ! Start over with new direction if we tried too many times
                    if (ctr .ge. maxctr) exit inloop
   
                 enddo inloop
                 mh%sx(:, idirection) = x2
                 mh%serr(idirection) = err2
             enddo
   
             ! Pick the next x
             call next_x(mh, x, err1)
   
          enddo outloop
   
          ! Seems we have fixed the Born charges! Return the neat ones
          x_Z = x
          ! And cleanup
          deallocate(x)
          deallocate(x1)
          deallocate(x2)
          deallocate(grad)
          deallocate(searchdir)
          deallocate(wD)
          deallocate(dynmat_reciprocal)
          deallocate(wZ)
          deallocate(wf)
      end block opt
   end subroutine
   
   
   !> longrange dipole-dipole dynamical matrix
   subroutine longrange_dynamical_matrix( &
     ew, uc, q, born_effective_charges, born_onsite_correction, eps, D, reconly, chgmult)
     type(ewald_info), intent(in) :: ew !< Ewald information 
     type(uc_tdep), intent(in) :: uc !< crystal structure
     real(dp), dimension(3), intent(in) :: q !< q-vector
     real(dp), dimension(:, :, :), intent(in) :: born_effective_charges !< Born effective charges
     real(dp), dimension(:, :, :), intent(in) :: born_onsite_correction !< On-site correction for Born charges
     real(dp), dimension(3, 3), intent(in) :: eps !< dielectric constant
     complex(dp), dimension(:, :, :, :), intent(out) :: D !< dynamical matrix
     logical, intent(in), optional :: reconly !< only calculate the reciprocal sum
     logical, intent(in), optional :: chgmult !< multiply charges
   
     complex(dp), dimension(:, :, :, :), allocatable :: D0, Dx0, Dy0, Dz0
     complex(dp), dimension(:, :, :), allocatable :: DQLx1, DQLx2, DQLy1, DQLy2, DQLz1, DQLz2
     complex(dp), dimension(:, :, :), allocatable :: DL1, DL2, DL3
     real(dp), dimension(:, :), allocatable :: ucvl
     real(dp), dimension(3, 3) :: inveps
     real(dp) :: dete
     real(dp) :: rcart1(3), rcart2(3) 
     integer :: npair
     logical :: longrange, mult
   
     ! Set some things and precalculate
     init: block
        integer :: a1, a2
        ! First figure out what we are going to calculate
        if (present(reconly)) then
           longrange = reconly
        else
           longrange = .false.
        endif
        if (present(chgmult)) then
           mult = chgmult
        else
           mult = .true.
        endif
   
        ! inverse dielectric tensor
        inveps = invert3x3matrix(eps)
        inveps = chop_real(inveps, sq_tol)
        ! inverse determinant
        dete = 1.0_dp/sqrt(determ_real(eps))
   
        ! Make some space for temporary things
        allocate(D0(3, 3, uc%na, uc%na))
        D0 = 0.0_dp
   
        ! Count the ucvec thing
        npair = 0
        do a1 = 1, uc%na
        do a2 = a1, uc%na
           npair = npair + 1
        enddo
        enddo
        allocate(ucvl(3, npair))
        ucvl = 0.0_dp
        npair = 0
        do a1 = 1, uc%na
        do a2 = a1, uc%na
           npair = npair + 1
           rcart1 = uc%r(:, a1)
           rcart2 = uc%r(:, a2)
           call cryst_to_cart(1, rcart1, uc%latticevectors, 1)
           call cryst_to_cart(1, rcart2, uc%latticevectors, 1)
           rcart1 = chop_real(rcart1, sq_tol)
           rcart2 = chop_real(rcart2, sq_tol)
           ucvl(:, npair) = chop_real(rcart2 - rcart1, sq_tol)
           !ucvl(:, npair) = chop_real(p%rcart(:, a2) - p%rcart(:, a1), sq_tol)
        enddo
        enddo
        allocate(DL1(3, 3, npair))
        allocate(DL2(3, 3, npair))
        allocate(DL3(3, 3, uc%na))
        DL1 = 0.0_dp
        DL2 = 0.0_dp
        DL3 = 0.0_dp
     end block init
   
     ! First the reciprocal part
     reciprocalspace: block
        complex(dp), dimension(3) :: v0
        complex(dp) :: expikr, Chi
        real(dp), dimension(3, 3) :: kkx, kky, kkz, kk
        real(dp), dimension(3) :: Gvec, Kvec, tauvec, Keps
        real(dp) :: inv4lambda2, knorm, expLambdaKnorm, ikr, invknorm, Kx, Ky, Kz, partialChi, f0
        integer :: ig, ip
   
        inv4lambda2 = 1.0_dp/(4.0_dp*(ew%lambda**2))
        do ig = 1, ew%n_Gvector
           Gvec = ew%Gvec(:, ig)
           Kvec = (Gvec + q)*tpi
           Kvec = chop_real((Gvec + q)*tpi, sq_tol)
           if (sqnorm(Kvec) .lt. sq_tol) cycle
           ! K-norm with dielectric metric
           knorm = dot_product(Kvec, matmul(eps, Kvec))
           invknorm = 1.0_dp/knorm
           expLambdaKnorm = exp(-knorm*inv4lambda2)
           ! Half of the Chi-function
           partialChi = expLambdaKnorm*invknorm
           ! do only half of all pairs
           ! Outer product of the vectors
           kk = real_outerproduct(Kvec, Kvec)
           do ip = 1, npair
              ikr = -dot_product(Kvec, ucvl(:, ip))
              expikr = cmplx(cos(ikr), sin(ikr), dp)
              Chi = partialChi*expikr
              DL1(:, :, ip) = DL1(:, :, ip) + kk*Chi
           enddo
        enddo
        ! and multiply in the prefactor
        !f0 = 4.0_dp*pi/p%volume
        f0 = 4.0_dp*pi/abs(determ_real(uc%latticevectors))
        DL1 = DL1*f0
     end block reciprocalspace
   
     if (longrange .eqv. .false.) then
        ! Realspace part
        realspace: block
           complex(dp) :: expikr
           complex(dp), dimension(3, 3) :: Hp
           real(dp), dimension(3, 3) :: H
           real(dp), dimension(3) :: Rvec, Dvec, delta
           real(dp) :: ikr, bigD, lambdacub
           integer :: ir, ip
   
           lambdacub = ew%lambda**3
           DL2 = 0.0_dp
           do ir = 1, ew%n_Rvector
              Rvec = ew%Rvec(:, ir)
              do ip = 1, npair
                 Dvec = chop_real(Rvec + ucvl(:, ip), sq_tol)
                 ! exp(iq.r)
                 ikr = dot_product(q*tpi, Rvec)
                 expikr = cmplx(cos(ikr), sin(ikr), dp)
                 ! get the delta, something with inverse dielectric tensor as the metric:
                 delta = matmul(inveps, Dvec)
                 ! get D, distance 
                 bigD = sqrt(dot_product(delta, Dvec))
                 ! compute anisotropic dipole interaction 
                 H = anisotropic_dipole_int(ew%lambda*delta, ew%lambda*bigD, inveps)
                 ! Multiply in all the prefactors
                 Hp = H*lambdacub*expikr*dete
                 ! Add it up
                 DL2(:, :, ip) = DL2(:, :, ip) + Hp
              enddo
           enddo
        end block realspace
     endif
   
     if (longrange .eqv. .false.) then
        ! Connecting part
        connect: block
           integer :: a1
           real(dp), parameter :: fouroverthreesqrtpi = 0.7522527780636750_dp
           DL3 = 0.0_dp
           do a1 = 1, uc%na
              DL3(:, :, a1) = (ew%lambda**3)*inveps*dete*fouroverthreesqrtpi
           enddo
        end block connect
     endif
   
     ! And the final touches
     finalize: block
        integer :: a1, a2, i, j, ii, jj, ip
        ! I only calculated half of the things, so fix that. Also add things
        ! together.
        D0 = 0.0_dp
        ip = 0
        do a1 = 1, uc%na
        do a2 = a1, uc%na
           ip = ip + 1
           if (a1 .eq. a2) then
               D0(:, :, a1, a2) = DL1(:, :, ip) - DL2(:, :, ip) - DL3(:, :, a1)
           else
               D0(:, :, a1, a2) = DL1(:, :, ip) - DL2(:, :, ip)
               do i = 1, 3
               do j = 1, 3
                   D0(j, i, a2, a1) = conjg(D0(i, j, a1, a2))
               enddo
               enddo
           endif
        enddo
        enddo
        if (mult) then
           ! Multiply in the charges
           D = 0.0_dp
           do a2 = 1, uc%na
           do a1 = 1, uc%na
               do j = 1, 3
               do i = 1, 3
                   do jj = 1, 3
                   do ii = 1, 3
                       ! Reasonably sure this is the correct way of doing it. Or it could be transposed, who knows.
                       !D(i,j,a1,a2)=D(i,j,a1,a2)+fc%atom(a1)%Z(i,ii)*fc%atom(a2)%Z(j,jj)*D0(ii,jj,a1,a2)
                       D(i, j, a1, a2) &
                           = D(i, j, a1, a2) &
                             + born_effective_charges(i, ii, a1)*born_effective_charges(j, jj, a2)*D0(ii, jj, a1, a2)
                   enddo
                   enddo
               enddo
               enddo
           enddo
           enddo
           ! Add the on-site correction
           do a1 = 1, uc%na
               D(:, :, a1, a1) = D(:, :, a1, a1) + born_onsite_correction(:, :, a1)
           enddo
        else
           ! Skip multiplying in the charges
           D = D0
        endif
        ! Remove tiny ugly numbers
        D = chop_complex(D, sq_tol)
     end block finalize
      ! Cleanup
      deallocate(D0)
      deallocate(Dl1)
      deallocate(Dl2)
      deallocate(Dl3)
      deallocate(ucvl)
   end subroutine

end module lattice_tdep
