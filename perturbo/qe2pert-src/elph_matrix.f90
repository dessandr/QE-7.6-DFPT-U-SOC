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
!  based on the workflow in Phonon/PH/elphon.f90
!  Calculation of the electron-phonon matrix elements el_ph_mat
!     <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
!  
! Note: the noncollinear magnetic local-potential action supports
! nspin_mag=4.  Magnetic q-star reconstruction remains unsupported, so
! SOC+U currently requires one explicitly calculated response per q point.
!
! NOTE: Segment fault may occur if OpenMP private space runs out. 
!   We allocated memory for large arrays in each thread. however, 
!   the private stack space each thread has is limited (defined by OMP_STACKSIZE).
!   The default value is implementation dependent and is usually quite small. 
!   A segmentation fault may appear if running out of space. We can specify 
!   the size of the private space by 'export OMP_STACKSIZE=xxK (or xxM, xxG)' 
!   or use OpenMP runtime subroutines: kmp_{set, get}_stacksize_s().
!
! Maintenance:
!===============================================================================

module elph_matrix
   use kinds, only: dp
   use constants, only: tpi
   use gvect, only: ngm, g
   use fft_base, only: dfftp, dffts
   use fft_interfaces, only: fwfft, invfft
   use fft_wave, only: fwfft_wave, invfft_wave  ! QE 7.6: module, explicit npwx
   use cell_base, only: tpiba2, omega, at, bg
   use atom,   ONLY : msh, rgrid
   use pwcom, only: npwx
   use gvecs, only: doublegrid
   use uspp_param, only: upf, nhm, nh
   use uspp,  only: nlcc_any, nkb, okvan
   use lsda_mod, only : nspin, lsda, current_spin
   use ldaU, only: lda_plus_u, nwfcU
   use noncollin_module, only: noncolin, npol, nspin_mag, lspinorb
   use ions_base, only: nat, ityp, tau, ntyp => nsp
   use uspp_init, only: init_us_2
   use eqv,       only: vlocq
   implicit none
   public

   complex(dp), parameter :: cone = (1.0_dp, 0.0_dp)
   complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)

   complex(dp), allocatable, save :: int1(:,:,:,:,:),    int2(:,:,:,:,:)
   complex(dp), allocatable, save :: int1_nc(:,:,:,:,:), int2_so(:,:,:,:,:,:)
   complex(dp), allocatable, save :: int3(:,:,:,:,:),    int3_nc(:,:,:,:,:)
   complex(dp), allocatable, save :: drc(:,:)
   
   public :: compute_elph_matrix, validate_soc_u_load_only, &
      validate_soc_u_production
contains
! 
subroutine compute_elph_matrix()
   use buffers, only: open_buffer, save_buffer, close_buffer
   use fft_interfaces, only: fft_interpolate
   use io_global, only: stdout
   use symm_base, only: nsym, s, sr, invs, irt
   !
   use eqv, only: dmuxc
   ! xq: q in cartesian coordinate
   use qpoint, only: xq, eigqts
   !
   use input_param,  only: lwannier, num_band
   use qe_mpi_mod,   only: mp_split_pools
   use lattice_data, only: numq, xq_tot, iq2irrq
   use electronic_data, only: num_kpts, spin_updn, evc_sub, xk_tot, ktable, ngk_tot, &
      igk_k_tot, kdim, rot_wan, kpt2num, becp1_tot, nb_sub
   !
   use dvscf_mod,    only: obtain_dvscf
   !
   !for LDA+U 
   use dvhub_mod, only: ldim, setup_hubbard, obtain_dnshub

   implicit none
   complex(dp), pointer :: dvscfins(:,:,:)
   complex(dp), allocatable, target :: dvscf(:,:,:)
   complex(dp), allocatable :: uact(:), vkb_kq(:,:), ephmat(:,:,:,:), dvscf_irrq(:,:,:)
   complex(dp), allocatable :: dvpsi_noloc(:,:), dvpsi(:,:), aux1(:,:), aux2(:), aux(:)
   complex(dp), allocatable :: vkb_kk(:,:), wfck(:,:), wfckq(:,:), swfck(:,:), swfckq(:,:)
   integer, external :: find_free_unit
   complex(dp), external :: zdotc
   !for DFPT+U
   complex(dp), allocatable :: dnsbare_irrq(:,:,:,:,:,:), dnsbare(:,:,:,:,:,:), &
      dnsscf_irrq(:,:,:,:,:,:), dnsscf(:,:,:,:,:,:)
   
   logical :: exst, lnew_irrq
   integer :: nmodes, iun_eph, eph_nword, i, iq, iq_st, iq_end, nq_loc, irec, icart
   integer :: na, ib, jb, ik, ikq, ip, npwq, nt, current_irrq, is, icount, iter, irrq
   real(dp) :: arg, xkk(3), xkq(3), xkqg(3)
   integer,  allocatable :: igk_kqg(:)

   write(stdout,'(a)') "start elph_matrix:"

   !check if all xk+xq are on the grid
   do ik = 1, num_kpts
   do iq = 1, numq
      xkq = xk_tot(:,ik) + xq_tot(:,iq)
      ikq = kpt2num(xkq, kdim)
      if( ikq < 0 ) call errore('elph_matrix',"k- and q-grid are incommensurate",1)
   enddo; enddo
   
   !setup spin channel, this is important. Although all the subroutines in qe2pert-src
   ! use spin_updn, there is subroutines from PW/PH called here use lsda_mod/current_spin
   ! to determine the spin chnnel, such as 'compute_deff' (from PW/src/compute_deff.f90) 
   ! called in 'ph_dvqpsi_us_noloc'. so we set current_spin to the value of spin_updn
   current_spin = spin_updn
   !
   !setup basic parameters
   nmodes = 3*nat
   iun_eph = find_free_unit()
   eph_nword = nb_sub * nb_sub * num_kpts
   !io_level = 1, save eph matrix to disk
   call open_buffer(iun_eph, 'ephmat', eph_nword, 1, exst)
   
   ! Set non linear core correction stuff
   nlcc_any = any( upf(1:ntyp)%nlcc )
   !
   call allocate_work_arrays( nlcc_any )
   allocate( ephmat(nb_sub, nb_sub, num_kpts, nmodes), vlocq(ngm, ntyp) )
   !
   if(.not. allocated(eigqts)) allocate(eigqts(nat))
   if(.not. allocated(dmuxc) ) allocate(dmuxc(dfftp%nnr, nspin_mag, nspin_mag))
   !work arrays
   allocate( dvscf_irrq(dfftp%nnr, nspin_mag, nmodes) )
   allocate( dvscf(dfftp%nnr, nspin_mag, nmodes) )
   if(doublegrid) allocate( dvscfins(dffts%nnr, nspin_mag, nmodes) )
   ! pre-run fwfft to setup FFT tables for dffts, to avoid race condition 
   !   when calling fwfft and invffw in openMP parallel region.
   allocate( aux(dffts%nnr) );  aux = (0.0_dp,  0.0_dp)
   call fwfft('Wave', aux, dffts)
   deallocate( aux )
   
   ! setup, see PH/phq_setup.f90
   !3) Computes the derivative of the XC potential
   ! explicit output to eqv: dmuxc
   call setup_dmuxc()
   ! setup all gradient correction stuff
   call setup_dgc()
   if(lda_plus_u) then
      ! DFPT+U, see Phonon/PH/phq_setup.f90, set ldim, effU, ofsbeta, and d1,d2,d3
      call setup_hubbard()
      
      !allocate space
      allocate(dnsbare_irrq (ldim, ldim, nspin, nat, 3, nat), &
                dnsscf_irrq (ldim, ldim, nspin, nat, 3, nat), &
                    dnsbare (ldim, ldim, nspin, nat, 3, nat), &
                     dnsscf (ldim, ldim, nspin, nat, 3, nat)  )
      !init
      dnsbare_irrq = cmplx(0._dp, 0._dp, kind=dp)
      dnsscf_irrq  = cmplx(0._dp, 0._dp, kind=dp)
   endif
   
   current_irrq = -1
   dvscf_irrq = cmplx(0._dp, 0._dp, kind=dp)
   !
   !mpi parallelization
   call mp_split_pools(numq, iq_st, iq_end, nq_loc)
   ! for progress tracking
   icount = 0
   !
   do i = 1, nq_loc
      iq = iq_st + i - 1
      
      !Note: xq and eigqts are global variables defined in module qpoint
      ! cutoff_fact_qg() and cutoff_lr_Vlocq() requires xq in qpoint (in cartesian coord.)
      ! xq_tot(:,iq): the iq-th point in crystal coordinate; TC not sure this comment 
      ! applies anymore with the new code (all inside init_vlocq)
      xq = xq_tot(:,iq)
      !transfer to cartesian coordinate
      call cryst_to_cart(1, xq, bg, 1)

      !(following QE/PHonon/PH/phq_init)
      !compute  exp(-i xq times tau)
      do na = 1, nat
         arg = dot_product(xq, tau(:,na)) * tpi
         eigqts(na) = cmplx(cos(arg), -sin(arg), kind=dp)
      enddo
      
      !initializatiion (q-dependent) (following QE/PHonon/PH/phq_init)
      ! a0) compute rhocore for each atomic-type if needed for nlcc
      if(nlcc_any)  call set_drhoc(xq, drc)
      !
      ! b) the fourier components of the local potential at q+G
      call init_vlocq ( xq )
      ! index of the corresponding irreducible-q of xq
      irrq = iq2irrq(iq)
      !if the current_irrq is not the irreducible-q of iq, 
      lnew_irrq = (irrq .ne. current_irrq)
      !  then update current_irrq to the new one
      if(lnew_irrq)  current_irrq = irrq

      ! self-consistent part:Eq. B30 of PRB 64 235118
      call obtain_dvscf(current_irrq, dvscf_irrq, iq, dvscf, lnew_irrq)

      ! DFPT+U, obtain dnsbare and dnsscf (include dnsorth) for the current q
      if(lda_plus_u) then
         !get dnsbare, i.e. the bare variation of ns for all cartesian coordinates
         call obtain_dnshub("bare_pattern", current_irrq, dnsbare_irrq, iq, dnsbare, lnew_irrq)
         !get dnsscf for all cartesian coordinates, dnsorrth in also included in USPP case.
         call obtain_dnshub("scf",  current_irrq, dnsscf_irrq,  iq, dnsscf, lnew_irrq)
      endif
      
      ! in the case of us-pp
      if(okvan) then
         ! compute int1, int2, int1_nc, int2_so for xq
         call ph_dvanqq(xq, eigqts, vlocq)
         ! preparing int3 for the 2nd term in Eq.B30 of PRB 64 235118
         ! compute int3 and int3_nc
         call ph_newdq(xq, eigqts, dvscf, nmodes)
      endif
      !!end of initialization for the current q
      
      !compute 1st term in Eq.B30 of PRB 64 235118
      IF (doublegrid) THEN
         do ip = 1, nmodes
         do is = 1, nspin_mag
            call fft_interpolate(dfftp, dvscf(:,is,ip), dffts, dvscfins(:,is,ip))
         enddo; enddo
      ELSE
         dvscfins => dvscf
      ENDIF
      !
      ephmat = cmplx(0.E0_dp, 0.E0_dp, kind=dp)
      !
!$omp parallel default(shared) private(ik, xkk, xkqg, ikq, npwq, xkq, ib, jb, &
!$omp&   dvpsi, dvpsi_noloc, vkb_kq, aux1, aux2, igk_kqg, iter, ip, uact, &
!$omp&   vkb_kk, wfck, wfckq, swfck, swfckq, icart, na)
      allocate( dvpsi_noloc(npwx*npol, num_band), dvpsi(npwx*npol, nb_sub) )
      allocate( aux1(dffts%nnr, npol), aux2(npwx*npol), vkb_kq(npwx, nkb), uact(nmodes))
      allocate( igk_kqg(npwx) )
      !
      if(lda_plus_u) then
         allocate( vkb_kk(npwx, nkb) )
         allocate( wfck(npwx*npol, nwfcU), wfckq(npwx*npol, nwfcU) )
         allocate( swfck(npwx*npol, nwfcU), swfckq(npwx*npol, nwfcU) )
      endif
      !
!$omp do schedule(guided)
      do ik = 1, num_kpts
         !xk+xq could be xk_tot(:,ikq)+G, the G need to be taken into account
         xkqg = xk_tot(:,ik) + xq_tot(:,iq)
         ikq = ktable( kpt2num(xkqg, kdim) )
         !transfer from crystal to cartesian coordinate, xk+xq with possible G
         call cryst_to_cart(1, xkqg, bg, 1)
         npwq = ngk_tot(ikq)
         ! xk+xq removed possible G, xkq is on k-grid.
         xkq = xk_tot(:,ikq)
         ! to cartesian coordiante
         call cryst_to_cart(1, xkq, bg, 1)
         !generate igk_k for xk+xq
         igk_kqg(:) = 0
         call get_igk_kqg(npwq, igk_k_tot(1,ikq), xkq, xkqg, igk_kqg(1))
         !
         !xk transfer from crystal to cartesian coordinate
         xkk  = xk_tot(:,ik)
         call cryst_to_cart(1, xkk, bg, 1)
            
         !compute vkb_kq(k+q+G)
         call init_us_2(npwq, igk_kqg(1), xkqg, vkb_kq)
         !
         if(lda_plus_u) then
            call init_us_2(ngk_tot(ik), igk_k_tot(1,ik), xkk, vkb_kk)
            !compute the atomic orbitals at k and k+q: 
            !                      \phi(k), \phi(k+q), S\phi(k), S\phi(k+q)
            call ph_lr_orthoUwfc(xkk, ngk_tot(ik), igk_k_tot(1,ik), xkqg, npwq, &
               igk_kqg(1), vkb_kk, vkb_kq, wfck, swfck, wfckq, swfckq)
         endif
      
         do na = 1, nat
         do icart = 1, 3
            ip = (na-1)*3 + icart
            !
            uact     = cmplx(0.0_dp, 0.0_dp, kind=dp)
            uact(ip) = cmplx(1.0_dp, 0.0_dp, kind=dp)
            !
            !compute dvscf_q*psi_{k,i}>
            !
            !non-self-consistent part: Eq. B29 of PRB 64 235118
            call ph_dvqpsi_us_local &
               (ik, xq, eigqts, npwq, igk_kqg(1), uact, dmuxc, vlocq, dvpsi, .false.)
            
            if(lda_plus_u) then
               if(noncolin) then
                  call ph_dvqhub_barepsi_nc(ik, xkk, xkqg, npwq, &
                     igk_kqg(1), wfck, swfck, wfckq, swfckq, &
                     dnsbare(1,1,1,1,icart,na), icart, na, &
                     dvpsi_noloc(1,1))
                  dvpsi(:,1:nb_sub) = dvpsi(:,1:nb_sub) + &
                     dvpsi_noloc(:,1:nb_sub)

                  call ph_adddvhubscf_nc(ik, swfck, swfckq, &
                     dnsscf(1,1,1,1,icart,na), dvpsi_noloc(1,1))
                  dvpsi(:,1:nb_sub) = dvpsi(:,1:nb_sub) + &
                     dvpsi_noloc(:,1:nb_sub)
               else
                  ! bare part of dV_hub (temporarily to dvpsi_noloc)
                  call ph_dvqhub_barepsi_us(ik, xkk, xkqg, npwq, &
                     igk_kqg(1), vkb_kk, vkb_kq, wfck, swfck, wfckq, &
                     swfckq, dnsbare(1,1,1,1,icart,na), icart, na, &
                     dvpsi_noloc(1,1))
                  dvpsi(:,1:nb_sub) = dvpsi(:,1:nb_sub) + &
                     dvpsi_noloc(:,1:nb_sub)

                  ! scf part of dV_hub (including dnsorth for USPP)
                  call ph_adddvhubscf(ik, npwq, swfck, swfckq, &
                     dnsscf(1,1,1,1,icart,na), dvpsi_noloc(1,1))
                  dvpsi(:,1:nb_sub) = dvpsi(:,1:nb_sub) + &
                     dvpsi_noloc(:,1:nb_sub)
               endif
            endif
            
            !non-local part
            call ph_dvqpsi_us_noloc &
               (ik, xkqg, npwq, igk_kqg(1), uact, vkb_kq, dvpsi_noloc)
            
            !apply the rotation if compute in wannier gauge
            if(lwannier) then
               ! if lwannier = .true. then nb_sub is the same as num_wann
               ! dvpsi = dvpsi_noloc * rot_wan(1:num_band,1:num_wann,ik) + dvpsi
               call zgemm('N','N', npwx*npol, nb_sub, num_band, cone, dvpsi_noloc, &
                        npwx*npol, rot_wan(1,1,ik), num_band, cone, dvpsi, npwx*npol)
            else
               ! if lwannier = .false. then nb_sub is the same as num_band
               dvpsi(:,1:nb_sub) = dvpsi(:,1:nb_sub) + dvpsi_noloc(:,1:nb_sub)
            endif
               
            if(okvan) then
               !add the 2nd term in Eq.B30 of PRB 64 235118
               call ph_adddvscf(npwq, ip, becp1_tot(ik), vkb_kq, dvpsi_noloc)
               !apply the rotation
               if(lwannier) then
                  ! dvpsi = dvpsi_noloc * rot_wan(1:num_band,1:num_wann,ik) + dvpsi
                  call zgemm('N','N', npwx*npol, nb_sub, num_band, cone, dvpsi_noloc, &
                        npwx*npol, rot_wan(1,1,ik), num_band, cone, dvpsi, npwx*npol)
               else
                  dvpsi(:,1:nb_sub) = dvpsi(:,1:nb_sub) + dvpsi_noloc(:,1:nb_sub)
               endif
            endif

            !NOTE: we enforce dffts%has_task_groups =.false. (only support 1 processor/pool )
            do ib = 1, nb_sub
               aux2(:) = czero
               !self-consistent part: first term of Eq. B30
               !call cft_wave (ik, evc_sub(1, ib, ik), aux1, +1)
               ! Inverse Fourier transform of a wavefunction: G-space to r-space
               call invfft_wave(npwx, ngk_tot(ik), igk_k_tot(1,ik), evc_sub(1,ib,ik), aux1)
               !
               !! NOTE: current_spin and isk only matter in lsda calc., 
               ! and the number of k-points is doubled:  the first half for spin up (current_spin=1),
               !  the second half for spin down (current_spin=2).
               ! for non-spin and non-collinear calc. current_spin (=1) does not matter.
               call apply_dpot(dffts%nnr, aux1, dvscfins(1,1,ip), spin_updn)
               !
               ! Fourier transform of a wavefunction: r-space to G-space 
               !call cft_wave (ik, aux2, aux1, -1)
               call fwfft_wave(npwx, npwq, igk_kqg(1), aux2, aux1)
               !
               dvpsi(:,ib) = dvpsi(:,ib) + aux2(:)
            enddo
            
            ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
            do ib = 1, nb_sub
            do jb = 1, nb_sub
               ephmat(jb,ib,ik,ip) = zdotc(npwq, evc_sub(1,jb,ikq), 1, dvpsi(1,ib), 1)
               if(noncolin) then
                  ephmat(jb,ib,ik,ip) = ephmat(jb,ib,ik,ip) + &
                    zdotc(npwq, evc_sub(npwx+1,jb, ikq), 1, dvpsi(npwx+1,ib), 1)
               endif
            enddo
            enddo
         enddo ! icart
         enddo ! na
   
         !tracking progress
!$omp atomic update
         icount = icount + 1
         iter = icount
         if( mod(iter, merge(num_kpts/nat,1,num_kpts>nat)).eq.0 .or. iter.eq.(nq_loc*num_kpts) ) then
!$omp critical (elph_matrix_progress)
            write(stdout, '(2x, a, f7.2, a1)') &
               "progress: ", (100.0_dp*iter)/real(nq_loc*num_kpts, dp), "%"
!$omp end critical (elph_matrix_progress)
         endif
         !
      enddo !ik
!$omp end do
      deallocate(vkb_kq, dvpsi, dvpsi_noloc, aux1, aux2, igk_kqg, uact)
      if(lda_plus_u) deallocate(vkb_kk, wfck, wfckq, swfck, swfckq)
!$omp end parallel
      
      !save result
      do ip = 1, nmodes
         irec = (i - 1) * nmodes + ip
         call save_buffer(ephmat(:,:,:,ip), eph_nword, iun_eph, irec)
      enddo
   enddo ! iq
   call close_buffer(iun_eph, 'KEEP')
   
   !release memory
   call deallocate_work_arrays()
   if(doublegrid) deallocate(dvscfins)
   if(lda_plus_u) deallocate(dnsbare_irrq, dnsbare, dnsscf_irrq, dnsscf)
   deallocate(dvscf_irrq, dvscf)
   deallocate(dmuxc, eigqts, ephmat, vlocq)
end subroutine compute_elph_matrix

subroutine validate_soc_u_capability()
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use paw_variables, only: okpaw
   use klist, only: lgauss, ltetra, two_fermi_energies
   use ldaU, only: lda_plus_u_kind, Hubbard_projectors, Hubbard_lmax, &
      Hubbard_U2, Hubbard_J0, Hubbard_J, Hubbard_alpha, Hubbard_alpha_back, &
      Hubbard_beta, is_hubbard, is_hubbard_back, hub_pot_fix
   use scf, only: rho, v
   use lattice_data, only: numq, nq_irr, iq2irrq, identity_op, &
      q2minus_q, q2irq_symop, pattern, magnetic_qstar, xq_tot, xq_irr
   use symm_base, only: nsym, sr, invs, t_rev
   implicit none

   integer :: expected_ldim, irrq, iq, isym
   real(dp) :: qtarget(3), qsource(3)

   if(.not. lda_plus_u) &
      call errore('validate_soc_u_capability', 'Hubbard U is not active', 1)
   if(.not. noncolin .or. .not. lspinorb) &
      call errore('validate_soc_u_capability', 'noncollinear SOC is required', 1)
   if(npol /= 2 .or. nspin /= 4 .or. nspin_mag /= 4) &
      call errore('validate_soc_u_capability', &
         'expected npol=2, nspin=4, nspin_mag=4', 1)
   if(okvan .or. okpaw) &
      call errore('validate_soc_u_capability', &
         'only norm-conserving pseudopotentials are supported', 1)
   if(any(.not. upf(1:ntyp)%has_so)) &
      call errore('validate_soc_u_capability', &
         'fully relativistic pseudopotentials are required', 1)
   if(lda_plus_u_kind /= 0) &
      call errore('validate_soc_u_capability', &
         'only simplified Dudarev DFT+U is supported', 1)
   if(trim(Hubbard_projectors) /= 'atomic') &
      call errore('validate_soc_u_capability', &
         'Hubbard_projectors must be atomic', 1)
   if(.not. any(is_hubbard(1:ntyp))) &
      call errore('validate_soc_u_capability', &
         'no Hubbard manifold is active', 1)
   if(any(Hubbard_U2(1:ntyp) /= 0.0_dp) .or. &
      any(is_hubbard_back(1:ntyp))) &
      call errore('validate_soc_u_capability', &
         'multiple Hubbard manifolds are not supported', 1)
   if(any(Hubbard_J0(1:ntyp) /= 0.0_dp) .or. &
      any(Hubbard_J(:,1:ntyp) /= 0.0_dp) .or. &
      any(Hubbard_alpha(1:ntyp) /= 0.0_dp) .or. &
      any(Hubbard_alpha_back(1:ntyp) /= 0.0_dp) .or. &
      any(Hubbard_beta(1:ntyp) /= 0.0_dp) .or. hub_pot_fix) &
      call errore('validate_soc_u_capability', &
         'only fixed Dudarev Ueff is supported', 1)
   if(lgauss .or. ltetra .or. two_fermi_energies) &
      call errore('validate_soc_u_capability', &
         'only insulating occupations are supported', 1)

   if(.not. allocated(iq2irrq) .or. .not. allocated(q2irq_symop) .or. &
      .not. allocated(q2minus_q) .or. .not. allocated(pattern)) &
      call errore('validate_soc_u_capability', &
         'q-grid metadata are not allocated', 1)
   if(numq < 1 .or. nq_irr < 1 .or. nq_irr > numq) &
      call errore('validate_soc_u_capability', &
         'invalid full-grid/irreducible q counts', 1)
   if(size(iq2irrq) /= numq .or. size(q2irq_symop) /= numq .or. &
      size(q2minus_q) /= numq) &
      call errore('validate_soc_u_capability', &
         'inconsistent q-grid metadata dimensions', 1)
   if(size(pattern,1) /= 3*nat .or. size(pattern,2) /= 3*nat .or. &
      size(pattern,3) /= nq_irr) &
      call errore('validate_soc_u_capability', &
         'inconsistent displacement-pattern dimensions', 1)
   if(.not. magnetic_qstar) call errore('validate_soc_u_capability', &
      'magnetic response symmetry is not initialized', 1)
   if(any(q2irq_symop < 1) .or. any(q2irq_symop > nsym)) &
      call errore('validate_soc_u_capability', 'invalid magnetic symmetry index', 1)
   if(any(iq2irrq < 1) .or. any(iq2irrq > nq_irr)) &
      call errore('validate_soc_u_capability', &
         'invalid irreducible-q mapping', 1)
   do irrq = 1, nq_irr
      if(count(iq2irrq == irrq) < 1) &
         call errore('validate_soc_u_capability', &
            'each irreducible q must own a nonempty star', 1)
   enddo

   do iq = 1, numq
      isym = q2irq_symop(iq)
      if(invs(isym)<1 .or. invs(isym)>nsym) &
         call errore('validate_soc_u_capability', 'invalid inverse symmetry', 1)
      if(q2minus_q(iq) .neqv. (t_rev(isym)==1)) &
         call errore('validate_soc_u_capability', 'q mapping requires a non-magnetic Fourier reversal', 1)
      if(t_rev(isym)/=t_rev(invs(isym))) &
         call errore('validate_soc_u_capability', 'inconsistent inverse magnetic parity', 1)
      qtarget = matmul(bg, xq_tot(:,iq))
      qsource = matmul(bg, xq_irr(:,iq2irrq(iq)))
      if(q2minus_q(iq)) qtarget = -qtarget
      if(norm2(matmul(sr(:,:,isym),qtarget)-qsource)>1.e-5_dp) &
         call errore('validate_soc_u_capability', 'inconsistent magnetic q mapping', 1)
   enddo

   expected_ldim = 2*Hubbard_lmax + 1
   if(.not. allocated(rho%ns_nc) .or. .not. allocated(v%ns_nc)) &
      call errore('validate_soc_u_capability', &
         'ground-state Hubbard matrices are not allocated', 1)
   if(size(rho%ns_nc,1) /= expected_ldim .or. &
      size(rho%ns_nc,2) /= expected_ldim .or. &
      size(rho%ns_nc,3) /= nspin .or. size(rho%ns_nc,4) /= nat .or. &
      any(shape(rho%ns_nc) /= shape(v%ns_nc))) &
      call errore('validate_soc_u_capability', &
         'invalid ground-state Hubbard matrix dimensions', 1)
   if(any(.not. ieee_is_finite(real(rho%ns_nc,dp))) .or. &
      any(.not. ieee_is_finite(aimag(rho%ns_nc))) .or. &
      any(.not. ieee_is_finite(real(v%ns_nc,dp))) .or. &
      any(.not. ieee_is_finite(aimag(v%ns_nc)))) &
      call errore('validate_soc_u_capability', &
         'non-finite ground-state Hubbard data', 1)
end subroutine validate_soc_u_capability

subroutine validate_soc_u_production()
   use io_global, only: stdout, ionode
   implicit none

   call validate_soc_u_capability()
   call validate_soc_u_response_files()
   if(ionode) write(stdout,'(5x,a)') &
      'SOC_U_QE2PERT_PRODUCTION_CAPABILITY PASS'
end subroutine validate_soc_u_production

subroutine validate_soc_u_response_files()
   use, intrinsic :: iso_fortran_env, only: int64
   use io_global, only: ionode
   use mp, only: mp_barrier
   use mp_world, only: world_comm
   use input_param, only: prefix, phdir
   use lattice_data, only: nq_irr
   implicit none

   integer :: irrq
   integer(int64) :: file_size
   logical :: file_exists
   character(len=512) :: filename
   character(len=6), external :: int_to_char

   if(ionode) then
      do irrq = 1, nq_irr
         filename = trim(phdir) // trim(prefix) // '.dvscf_q' // &
            trim(int_to_char(irrq))
         file_size = -1_int64
         inquire(file=trim(filename), exist=file_exists, size=file_size)
         if(.not. file_exists) call errore( &
            'validate_soc_u_response_files', &
            trim(filename)//' does not exist', 1)
         if(file_exists .and. file_size <= 0_int64) call errore( &
            'validate_soc_u_response_files', &
            trim(filename)//' is empty', 1)

         filename = trim(phdir) // trim(prefix) // &
            '.dnsbare_pattern_q' // trim(int_to_char(irrq))
         file_size = -1_int64
         inquire(file=trim(filename), exist=file_exists, size=file_size)
         if(.not. file_exists) call errore( &
            'validate_soc_u_response_files', &
            trim(filename)//' does not exist', 1)
         if(file_exists .and. file_size <= 0_int64) call errore( &
            'validate_soc_u_response_files', &
            trim(filename)//' is empty', 1)

         filename = trim(phdir) // trim(prefix) // '.dnsscf_q' // &
            trim(int_to_char(irrq))
         file_size = -1_int64
         inquire(file=trim(filename), exist=file_exists, size=file_size)
         if(.not. file_exists) call errore( &
            'validate_soc_u_response_files', &
            trim(filename)//' does not exist', 1)
         if(file_exists .and. file_size <= 0_int64) call errore( &
            'validate_soc_u_response_files', &
            trim(filename)//' is empty', 1)
      enddo
   endif
   call mp_barrier(world_comm)
end subroutine validate_soc_u_response_files

subroutine validate_soc_u_load_only()
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use io_global, only: stdout, ionode
   use mp, only: mp_barrier
   use mp_world, only: world_comm
   use ldaU, only: Hubbard_lmax
   use scf, only: rho, v
   use lattice_data, only: numq, nq_irr, xq_tot, xq_irr, iq2irrq, pattern, q2irq_symop, q2minus_q
   use symm_base, only: invs
   use dvscf_mod, only: load_dvscf, get_dvscf_star
   use dvhub_mod, only: ldim, setup_hubbard, load_dns, get_dns_star
   implicit none

   complex(dp), allocatable :: dvscf(:,:,:), dvrot(:,:,:), dvback(:,:,:)
   complex(dp), allocatable :: dns(:,:,:,:,:,:), dnrot(:,:,:,:,:,:), dnback(:,:,:,:,:,:)
   real(dp) :: qtarget(3), qsource(3), residual, norm_residual
   integer :: isym
   integer :: expected_ldim, iq, irrq, is
   real(dp) :: checksum_real, checksum_imag, checksum_norm2, checksum_max

   call validate_soc_u_capability()
   call validate_soc_u_response_files()
   expected_ldim = 2*Hubbard_lmax + 1

   if(.not. ionode) then
      call mp_barrier(world_comm)
      return
   endif

   write(stdout,'(/5x,a)') 'SOC_U_QE2PERT_LOAD_ONLY_BEGIN'
   write(stdout,'(5x,a,6(1x,i0))') 'SOC_U_QE2PERT_CAPABILITY_CHECK', &
      npol, nspin, nspin_mag, nat, numq, nq_irr
   write(stdout,'(5x,a,4(1x,i0))') 'SOC_U_QE2PERT_GROUND_DIM_CHECK', &
      shape(rho%ns_nc)

   do is = 1, nspin
      checksum_real = sum(real(rho%ns_nc(:,:,is,:),dp))
      checksum_imag = sum(aimag(rho%ns_nc(:,:,is,:)))
      checksum_norm2 = sum(abs(rho%ns_nc(:,:,is,:))**2)
      checksum_max = maxval(abs(rho%ns_nc(:,:,is,:)))
      write(stdout,'(5x,a,1x,i0,4(1x,es24.16))') &
         'SOC_U_QE2PERT_RHO_NS_NC_CHECK', is, checksum_real, checksum_imag, &
         checksum_norm2, checksum_max

      checksum_real = sum(real(v%ns_nc(:,:,is,:),dp))
      checksum_imag = sum(aimag(v%ns_nc(:,:,is,:)))
      checksum_norm2 = sum(abs(v%ns_nc(:,:,is,:))**2)
      checksum_max = maxval(abs(v%ns_nc(:,:,is,:)))
      write(stdout,'(5x,a,1x,i0,4(1x,es24.16))') &
         'SOC_U_QE2PERT_V_NS_NC_CHECK', is, checksum_real, checksum_imag, &
         checksum_norm2, checksum_max
   enddo

   call setup_hubbard()
   if(ldim /= expected_ldim) &
      call errore('validate_soc_u_load_only', 'inconsistent Hubbard leading dimension', 1)

   allocate(dvscf(dfftp%nnr,nspin_mag,3*nat))
   allocate(dns(ldim,ldim,nspin,nat,3,nat))
   allocate(dvrot(dfftp%nnr,nspin_mag,3*nat), dvback(dfftp%nnr,nspin_mag,3*nat))
   allocate(dnrot(ldim,ldim,nspin,nat,3,nat), dnback(ldim,ldim,nspin,nat,3,nat))
   write(stdout,'(5x,a,3(1x,i0))') 'SOC_U_QE2PERT_DVSCF_DIM_CHECK', &
      shape(dvscf)
   write(stdout,'(5x,a,6(1x,i0))') 'SOC_U_QE2PERT_DNS_DIM_CHECK', &
      shape(dns)

   do iq = 1, numq
      irrq = iq2irrq(iq)
      isym = q2irq_symop(iq)
      qtarget = matmul(bg,xq_tot(:,iq))
      qsource = matmul(bg,xq_irr(:,irrq))
      write(stdout,'(5x,a,2(1x,i0),3(1x,f16.10))') &
         'SOC_U_QE2PERT_Q_CHECK', iq, irrq, xq_tot(:,iq)

      call load_dvscf(irrq, pattern(:,:,irrq), dvscf)
      if(any(.not. ieee_is_finite(real(dvscf,dp))) .or. &
         any(.not. ieee_is_finite(aimag(dvscf)))) &
         call errore('validate_soc_u_load_only', 'non-finite dvscf data', 1)
      do is = 1, nspin_mag
         checksum_real = sum(real(dvscf(:,is,:),dp))
         checksum_imag = sum(aimag(dvscf(:,is,:)))
         checksum_norm2 = sum(abs(dvscf(:,is,:))**2)
         checksum_max = maxval(abs(dvscf(:,is,:)))
         write(stdout,'(5x,a,3(1x,i0),4(1x,es24.16))') &
            'SOC_U_QE2PERT_DVSCF_CHECK', iq, irrq, is, checksum_real, &
            checksum_imag, checksum_norm2, checksum_max
      enddo

      call get_dvscf_star(qtarget, qsource, isym, q2minus_q(iq), dvscf, dvrot)
      call get_dvscf_star(qsource, qtarget, invs(isym), q2minus_q(iq), dvrot, dvback)
      residual = maxval(abs(dvback-dvscf))/max(1.e-30_dp,maxval(abs(dvscf)))
      norm_residual = abs(sum(abs(dvrot)**2)-sum(abs(dvscf)**2)) / &
         max(1.e-30_dp,sum(abs(dvscf)**2))
      write(stdout,'(5x,a,1x,i0,2(1x,es24.16))') 'SOC_U_QSTAR_DVSCF_ROUNDTRIP', iq, residual, norm_residual
      if(.not. ieee_is_finite(residual) .or. .not. ieee_is_finite(norm_residual) .or. &
         residual>1.e-8_dp .or. norm_residual>1.e-8_dp) &
         call errore('validate_soc_u_load_only', 'dvscf magnetic rotation roundtrip failed', 1)

      call load_dns('bare_pattern', irrq, pattern(:,:,irrq), dns)
      if(any(.not. ieee_is_finite(real(dns,dp))) .or. &
         any(.not. ieee_is_finite(aimag(dns)))) &
         call errore('validate_soc_u_load_only', 'non-finite dnsbare data', 1)
      call get_dns_star(qtarget, qsource, isym, q2minus_q(iq), dns, dnrot)
      call get_dns_star(qsource, qtarget, invs(isym), q2minus_q(iq), dnrot, dnback)
      residual = maxval(abs(dnback-dns))/max(1.e-30_dp,maxval(abs(dns)))
      norm_residual = abs(sum(abs(dnrot)**2)-sum(abs(dns)**2))/max(1.e-30_dp,sum(abs(dns)**2))
      write(stdout,'(5x,a,1x,i0,2(1x,es24.16))') 'SOC_U_QSTAR_DNSBARE_ROUNDTRIP', iq, residual, norm_residual
      if(.not. ieee_is_finite(residual) .or. .not. ieee_is_finite(norm_residual) .or. &
         residual>1.e-8_dp .or. norm_residual>1.e-8_dp) &
         call errore('validate_soc_u_load_only', 'dnsbare magnetic rotation roundtrip failed', 1)
      do is = 1, nspin
         checksum_real = sum(real(dns(:,:,is,:,:,:),dp))
         checksum_imag = sum(aimag(dns(:,:,is,:,:,:)))
         checksum_norm2 = sum(abs(dns(:,:,is,:,:,:))**2)
         checksum_max = maxval(abs(dns(:,:,is,:,:,:)))
         write(stdout,'(5x,a,3(1x,i0),4(1x,es24.16))') &
            'SOC_U_QE2PERT_DNSBARE_CHECK', iq, irrq, is, checksum_real, &
            checksum_imag, checksum_norm2, checksum_max
      enddo

      call load_dns('scf', irrq, pattern(:,:,irrq), dns)
      if(any(.not. ieee_is_finite(real(dns,dp))) .or. &
         any(.not. ieee_is_finite(aimag(dns)))) &
         call errore('validate_soc_u_load_only', 'non-finite dnsscf data', 1)
      call get_dns_star(qtarget, qsource, isym, q2minus_q(iq), dns, dnrot)
      call get_dns_star(qsource, qtarget, invs(isym), q2minus_q(iq), dnrot, dnback)
      residual = maxval(abs(dnback-dns))/max(1.e-30_dp,maxval(abs(dns)))
      norm_residual = abs(sum(abs(dnrot)**2)-sum(abs(dns)**2))/max(1.e-30_dp,sum(abs(dns)**2))
      write(stdout,'(5x,a,1x,i0,2(1x,es24.16))') 'SOC_U_QSTAR_DNSSCF_ROUNDTRIP', iq, residual, norm_residual
      if(.not. ieee_is_finite(residual) .or. .not. ieee_is_finite(norm_residual) .or. &
         residual>1.e-8_dp .or. norm_residual>1.e-8_dp) &
         call errore('validate_soc_u_load_only', 'dnsscf magnetic rotation roundtrip failed', 1)
      do is = 1, nspin
         checksum_real = sum(real(dns(:,:,is,:,:,:),dp))
         checksum_imag = sum(aimag(dns(:,:,is,:,:,:)))
         checksum_norm2 = sum(abs(dns(:,:,is,:,:,:))**2)
         checksum_max = maxval(abs(dns(:,:,is,:,:,:)))
         write(stdout,'(5x,a,3(1x,i0),4(1x,es24.16))') &
            'SOC_U_QE2PERT_DNSSCF_CHECK', iq, irrq, is, checksum_real, &
            checksum_imag, checksum_norm2, checksum_max
      enddo
   enddo

   deallocate(dvscf, dns, dvrot, dvback, dnrot, dnback)
   write(stdout,'(5x,a,1x,i0)') 'SOC_U_QE2PERT_LOAD_ONLY_RESULT PASS', numq
   write(stdout,'(5x,a/)') 'SOC_U_QE2PERT_LOAD_ONLY_END'
   call mp_barrier(world_comm)
end subroutine validate_soc_u_load_only

subroutine allocate_work_arrays( nlcc_any_ )
   implicit none
   logical, intent(in) :: nlcc_any_
   !work arrays, see PHonon/PH/allocate_phq.f90
   !if .TRUE. at least one pseudo is Vanderbilt
   if( okvan ) then
      allocate( int1(nhm, nhm, 3, nat, nspin_mag) )
      allocate( int2(nhm, nhm, 3, nat, nat) )
      allocate( int3(nhm, nhm, nat, nspin_mag, nat*3) )
      
      if( noncolin ) then
         allocate( int1_nc(nhm, nhm, 3, nat, nspin) )
         allocate( int3_nc(nhm, nhm, nat, nspin, nat*3) )
         if( lspinorb ) then
            allocate( int2_so(nhm, nhm, 3, nat, nat, nspin) )
         endif
      endif
   endif
   if( nlcc_any_ ) allocate ( drc(ngm, ntyp) )
end subroutine allocate_work_arrays

subroutine deallocate_work_arrays()
   implicit none
   !work arrays
   !if .TRUE. at least one pseudo is Vanderbilt
   if (okvan) then
      deallocate( int1, int2, int3)
      if( noncolin ) then
         deallocate( int1_nc, int3_nc )
         if( lspinorb ) deallocate( int2_so )
      endif
   endif
   if(allocated(drc)) deallocate(drc)
end subroutine deallocate_work_arrays

end module elph_matrix
