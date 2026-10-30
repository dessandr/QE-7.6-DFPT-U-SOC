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
!  compute inter-atomic force constants from the dynamical matrice
!
! Maintenance:
!   jpark: add acoustic sum rule adapted from PHonon/PH/matdyn.f90
!===============================================================================

subroutine calc_lattice_ifc(ph, qdim, nat, at, cryst_tau, numq, xq_tot, dynmat, dynmat_lr)
   use kinds, only: dp
   use constants, only: tpi
   use qe_mpi_mod, only: mp_sum, stdout, inter_pool_comm, mp_split_pools
   use force_constant, only: lattice_ifc, force_const, set_ws_cell_ph
   use lattice_data, only: ph_zstar, count_q2m1, count_q2m2, count_q2m3
   use input_param, only : asr
   use fft_scalar, only : cfft3d   
   use cell_base, only : ibrav
   implicit none
   type(lattice_ifc), intent(out) :: ph
   integer, intent(in) :: qdim(3), nat, numq
   ! atomic position in crystal coord.,  q-point in crystal coord.
   real(dp), intent(in) :: at(3,3), cryst_tau(3, nat), xq_tot(3, numq)
   complex(dp), intent(in) :: dynmat(numq, 3, 3, nat, nat), dynmat_lr(numq, 3, 3, nat, nat)
   !
   integer :: nq_tmp, ia, ja, iq, iq_st, iq_end, m, ir, nelem, nq_loc, &
              j1, j2, na1, na2
   real(dp) :: fac, rfac, rdotk, rvec(3), tau(3, nat)
   complex(dp) :: cfac
   complex(dp), allocatable :: dyn(:,:,:,:), dynmat_asr(:,:,:,:,:), &
                               phid(:,:,:,:,:), phid_lr(:,:,:,:,:)
   real(dp), allocatable :: frc(:,:,:,:,:), frc_lr(:,:,:,:,:)
   type(force_const), pointer :: ptr
   !set huang to true by default
   logical, parameter :: huang = .true.
   
   nelem  = ( nat * (nat + 1) ) / 2
   nq_tmp = qdim(1) * qdim(2) * qdim(3)
   !
   !tau: atomic position in cartesian coordinate
   tau(:,:) = cryst_tau(:,:)
   call cryst_to_cart(nat, tau, at, 1)

   if( nq_tmp .ne. numq ) &
      call errore('calc_lattice_ifc','inconsistent qdim and numq',1)
   !sanity check
   if(norm2(xq_tot(:,1)) > 1.0E-6_dp) &
      call errore('get_lattice_ifc','the first q is not Gamma',1)
   
   !set up wigner seitz cell
   call set_ws_cell_ph(ph, qdim, nat, at, cryst_tau)
   !init
   do m = 1, nelem
      ph%phi(m)%ifc(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)
   enddo
   
   allocate( dyn(3, 3, nat, nat) )
   allocate( dynmat_asr(numq, 3, 3, nat, nat) )
   dynmat_asr=dynmat

   if (asr /= 'no' ) then 
      !
      ! ASR in real space
      !
      allocate( phid(numq, 3, 3, nat, nat), phid_lr(numq, 3, 3, nat, nat))
      allocate( frc(numq, 3, 3, nat, nat), frc_lr(numq, 3,3, nat,nat) )
      phid(:,:,:,:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)
      phid_lr(:,:,:,:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)
      frc(:,:,:,:,:) = 0.0_dp
      frc_lr(:,:,:,:,:) = 0.0_dp
      !
      ! 1. dynmat (Perturbo format) -> phid (format used in q2r.f90)
      !       The q-grid index is rearranged here.
      !       Note that phid(numq,:,:,:,:) is read as phid(qdim(1),qdim(2),qdim(3),:,:,:,:)
      !       in dynmat2phid
      do iq = 1, numq
         call dynmat2phid(phid, dynmat, iq, numq, qdim(1), qdim(2), qdim(3), nat, &
                          count_q2m1(iq),count_q2m2(iq),count_q2m3(iq))
         call dynmat2phid(phid_lr, dynmat_lr, iq, numq, qdim(1), qdim(2), qdim(3), nat, &
                          count_q2m1(iq),count_q2m2(iq),count_q2m3(iq))
      enddo
      !
      ! 2. phid: reciprocal space -> real space
      !          required because set_asr is done in real space
      DO j1=1,3
         DO j2=1,3
            DO na1=1,nat
               DO na2=1,nat
                  CALL cfft3d ( phid (:,j1,j2,na1,na2), &
                       qdim(1),qdim(2),qdim(3), qdim(1),qdim(2),qdim(3), 1, 1 )
                  phid(:,j1,j2,na1,na2) = &
                       phid(:,j1,j2,na1,na2) / DBLE(qdim(1)*qdim(2)*qdim(3))

                  !for long-ranged part
                  CALL cfft3d ( phid_lr(:,j1,j2,na1,na2), &
                       qdim(1),qdim(2),qdim(3), qdim(1),qdim(2),qdim(3), 1, 1 )
                  phid_lr(:,j1,j2,na1,na2) = &
                       phid_lr(:,j1,j2,na1,na2) / DBLE(qdim(1)*qdim(2)*qdim(3))
               END DO
            END DO
         END DO
      END DO
      !
      ! 3. Apply ASR in real space: adapted from PHonon/PH/matdyn.f90
      !
      ! copy real part of phid to frc
      frc=real(phid)
      frc_lr = real(phid_lr)
      !
      CALL set_asr (asr, qdim(1), qdim(2), qdim(3), frc, frc_lr, ph_zstar, nat, ibrav, tau, at, huang)
      !
      !copy frc to phid, acoutic sum rule applied to the short-ranged part frc
      phid=frc
      !
      ! 4. phid: real space -> reciprocal space
      DO j1=1,3
         DO j2=1,3
            DO na1=1,nat
               DO na2=1,nat
                  CALL cfft3d ( phid (:,j1,j2,na1,na2), &
                       qdim(1),qdim(2),qdim(3), qdim(1),qdim(2),qdim(3), 1, -1 )
                  phid(:,j1,j2,na1,na2) = &
                       phid(:,j1,j2,na1,na2) * DBLE(qdim(1)*qdim(2)*qdim(3))
               END DO
            END DO
         END DO
      END DO
      !
      ! 5. phid (q2r format) -> dynmat(Perturbo format)
      do iq = 1, numq
          call phid2dynmat(phid, dynmat_asr, iq, numq, qdim(1), qdim(2), qdim(3), nat, &
                          count_q2m1(iq),count_q2m2(iq),count_q2m3(iq))
      enddo
      !
      deallocate(phid,frc, phid_lr, frc_lr)
      !
      ! ASR in real space finished
      !
   endif
 
   call mp_split_pools(numq, iq_st, iq_end, nq_loc)
   do iq = iq_st, iq_end
      !preparing D(q) and apply acoutic sum rule
      do ja = 1, nat
         do ia = 1, nat
            ! The hermiticity of dynmat is enforced in lattice_data/read_dynmat
            dyn(:, :, ia, ja) = dynmat_asr(iq, 1:3, 1:3, ia, ja)
         enddo
      enddo
      
      m = 0
      do ja = 1, nat
      do ia = 1, ja
         m = m + 1
         ptr => ph%phi(m)
         !
!$omp parallel do schedule(guided) default(shared) private(ir, rvec, rdotk, cfac)
         do ir = 1, ptr%ws_ph%nr
            rvec = ph%rvec_set(:, ptr%ws_ph%rvec(ir))
            rdotk = tpi * dot_product(xq_tot(:,iq), rvec)
            ! exp(-i k*R) / N_k
            cfac = cmplx( cos(rdotk), -sin(rdotk), kind=dp )
            ! ifc(:,:,ir) = ifc(:,:,ir) + dyn(:.:) * factor
            call zaxpy(9, cfac, dyn(1,1,ia,ja), 1, ptr%ifc(1,1,ir), 1)
         enddo
!$omp end parallel do
      !
      enddo; enddo
   enddo

   fac = 1.0_dp / real(numq, dp)
   !collect results
   do m = 1, nelem
      call mp_sum( ph%phi(m)%ifc, inter_pool_comm )
      ! including the ndeg and 1/numq
      do ir = 1, ph%phi(m)%ws_ph%nr
         rfac = fac / real( ph%phi(m)%ws_ph%ndeg(ir), dp)
         ph%phi(m)%ifc(:,:,ir) = ph%phi(m)%ifc(:,:,ir) * rfac
      enddo
      ! D(-q) = conjg(D(q)) => fc_tmp is real, just to double check.
      if( maxval(abs(aimag(ph%phi(m)%ifc(:,:,:)))) > 1.0E-8_dp ) then
         write(stdout,'(5x,a)') &
         "Warning (calc_lattice_ifc): non-negligible imaginary part in force constant."
      endif
   enddo
   
   deallocate( dyn, dynmat_asr )
end subroutine calc_lattice_ifc

subroutine dynmat2phid(phid, dynmat, iq, numq, nr1, nr2, nr3, nat, m1, m2, m3)
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, intent(in) ::  nr1, nr2, nr3, m1, m2, m3, nat, numq, iq
  COMPLEX(DP), intent(in) :: dynmat(numq,3,3,nat,nat)
  COMPLEX(DP), intent(inout) :: phid(nr1,nr2,nr3,3,3,nat,nat)
  !
  phid(m1,m2,m3,:,:,:,:)=dynmat(iq,:,:,:,:)
  !
  RETURN
end subroutine dynmat2phid

subroutine phid2dynmat(phid, dynmat, iq, numq, nr1, nr2, nr3, nat, m1, m2, m3)
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, intent(in) ::  nr1, nr2, nr3, m1, m2, m3, nat, numq, iq
  COMPLEX(DP), intent(inout) :: dynmat(numq,3,3,nat,nat)
  COMPLEX(DP), intent(in) :: phid(nr1,nr2,nr3,3,3,nat,nat)
  !
  dynmat(iq,:,:,:,:)=phid(m1,m2,m3,:,:,:,:)
  !
  RETURN
end subroutine phid2dynmat

! adapted from PHonon/PH/matdyn.f90
!----------------------------------------------------------------------
SUBROUTINE set_asr (asr, nr1, nr2, nr3, frc, frc_lr, zeu, nat, ibrav, tau_blk, at_blk, huang)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  CHARACTER (LEN=10), intent(in) :: asr
  LOGICAL, intent(in) :: huang
  INTEGER, intent(in) :: nr1, nr2, nr3, nat, ibrav
  REAL(DP), intent(in) :: tau_blk(3,nat), at_blk(3,3), frc_lr(nr1,nr2,nr3,3,3,nat,nat)
  REAL(DP), intent(inout) :: frc(nr1,nr2,nr3,3,3,nat,nat), zeu(3,3,nat)
  !
  INTEGER :: axis, n, i, j, na, nb, n1,n2,n3, m,p,k,l,q,r, i1,j1,na1, ip,ieq,neq
  INTEGER :: huang_set(4,15)
  REAL(DP) :: zeu_new(3,3,nat), tau(3,48), rcell(3), r_ws(3), eps
  REAL(DP), ALLOCATABLE :: frc_new(:,:,:,:,:,:,:)
  parameter (eps=1.0d-8)
  type vector
     real(DP),pointer :: vec(:,:,:,:,:,:,:)
  end type vector
  !
  type (vector) u(6*3*nat+15)
  ! These are the "vectors" associated with the sum rules on force-constants
  !
  integer :: u_less(6*3*nat),n_less,i_less
  ! indices of the vectors u that are not independent to the preceding ones,
  ! n_less = number of such vectors, i_less = temporary parameter
  !
  integer, allocatable :: ind_v(:,:,:)
  real(DP), allocatable :: v(:,:)
  ! These are the "vectors" associated with symmetry conditions, coded by
  ! indicating the positions (i.e. the seven indices) of the non-zero elements (there
  ! should be only 2 of them) and the value of that element. We do so in order
  ! to limit the amount of memory used.
  !
  real(DP), allocatable :: w(:,:,:,:,:,:,:), x(:,:,:,:,:,:,:)
  ! temporary vectors and parameters
  real(DP) :: scal,norm2, sum
  !
  real(DP) :: zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer :: zeu_less(6*3),nzeu_less,izeu_less
  ! indices of the vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP) :: zeu_w(3,3,nat), zeu_x(3,3,nat)
  ! temporary vectors

  ! Initialization. n is the number of sum rules to be considered (if asr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if axis='3')
  !
  if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'one-dim') &
                      .and.(asr.ne.'zero-dim').and.(asr.ne.'all')) then
     call errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
  endif
  !
  if(asr.eq.'simple') then
     !
     ! Simple Acoustic Sum Rule on effective charges
     !
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
              sum = sum + zeu(i,j,na)
           end do
           do na=1,nat
              zeu(i,j,na) = zeu(i,j,na) - sum/nat
           end do
        end do
     end do
     !
     ! Simple Acoustic Sum Rule on force constants in real space
     !
     do i=1,3
        do j=1,3
           do na=1,nat
              sum=0.0d0
               do nb=1,nat
                  do n1=1,nr1
                     do n2=1,nr2
                        do n3=1,nr3
                           sum=sum+frc(n1,n2,n3,i,j,na,nb)
                        end do
                     end do
                  end do
               end do
               frc(1,1,1,i,j,na,na) = frc(1,1,1,i,j,na,na) - sum
               !               write(6,*) ' na, i, j, sum = ',na,i,j,sum
            end do
         end do
      end do
      !
      return
      !
  end if
  if(asr.eq.'crystal') n=3
  if(asr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
        call errore('set_asr','too many directions of &
             & periodicity in 1D system',axis)
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(stdout,*) 'asr: rotational axis may be wrong'
     endif
     write(stdout,'("asr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(asr.eq.'zero-dim') n=6
  if(asr.eq.'all') n=21
  !
  ! Acoustic Sum Rule on effective charges
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the effective charges matrix on
  !
  zeu_u(:,:,:,:)=0.0d0
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu_new(i,j,na)=zeu(i,j,na)
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        ! These are the 3*3 vectors associated with the
        ! translational acoustic sum rules
        p=p+1
        zeu_u(p,i,j,:)=1.0d0
        !
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p=p+1
        do na=1,nat
           zeu_u(p,i,MOD(axis,3)+1,na)=-tau_blk(MOD(axis+1,3)+1,na)
           zeu_u(p,i,MOD(axis+1,3)+1,na)=tau_blk(MOD(axis,3)+1,na)
        enddo
        !
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           ! These are the 3*3 vectors associated with the
           ! three rotational sum rules (0D system - typ. molecule)
           p=p+1
           do na=1,nat
              zeu_u(p,i,MOD(j,3)+1,na)=-tau_blk(MOD(j+1,3)+1,na)
              zeu_u(p,i,MOD(j+1,3)+1,na)=tau_blk(MOD(j,3)+1,na)
           enddo
           !
        enddo
     enddo
  endif
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  !
  nzeu_less=0
  do k=1,p
     zeu_w(:,:,:)=zeu_u(k,:,:,:)
     zeu_x(:,:,:)=zeu_u(k,:,:,:)
     do q=1,k-1
        r=1
        do izeu_less=1,nzeu_less
           if (zeu_less(izeu_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
           zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
        endif
     enddo
     call sp_zeu(zeu_w,zeu_w,nat,norm2)
     if (norm2.gt.1.0d-16) then
        zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
     else
        nzeu_less=nzeu_less+1
        zeu_less(nzeu_less)=k
     endif
  enddo
  !
  ! Projection of the effective charge "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules
  !
  zeu_w(:,:,:)=0.0d0
  do k=1,p
     r=1
     do izeu_less=1,nzeu_less
        if (zeu_less(izeu_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        zeu_x(:,:,:)=zeu_u(k,:,:,:)
        call sp_zeu(zeu_x,zeu_new,nat,scal)
        zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
     endif
  enddo
  !
  ! Final substraction of the former projection to the initial zeu, to get
  ! the new "projected" zeu
  !
  zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
  call sp_zeu(zeu_w,zeu_w,nat,norm2)
  write(stdout,'(" Norm of the difference between old and new effective ", &
       & "charges: ",F25.20/)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection of zeu")')
  !do k=1,p
  !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
  !  call sp_zeu(zeu_x,zeu_new,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu(i,j,na)=zeu_new(i,j,na)
        enddo
     enddo
  enddo
  !
  ! Acoustic Sum Rule on force constants
  !
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the force-constants matrix on
  !
  do k=1,18*nat+15
     allocate(u(k) % vec(nr1,nr2,nr3,3,3,nat,nat))
     u(k) % vec (:,:,:,:,:,:,:)=0.0d0
  enddo
  ALLOCATE (frc_new(nr1,nr2,nr3,3,3,nat,nat))

  if (asr.eq.'all') then
     frc_new(:,:,:,:,:,:,:)=frc(:,:,:,:,:,:,:)+frc_lr(:,:,:,:,:,:,:)
  else
     frc_new(:,:,:,:,:,:,:)=frc(:,:,:,:,:,:,:)
  endif
  !
  p=0
  do i=1,3
     do j=1,3
        do na=1,nat
           ! These are the 3*3*nat vectors associated with the
           ! translational acoustic sum rules
           p=p+1
           u(p) % vec (:,:,:,i,j,na,:)=1.0d0
           !
        enddo
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        do na=1,nat
           ! These are the 3*nat vectors associated with the
           ! single rotational sum rule (1D system)
           p=p+1
           do nb=1,nat
              u(p) % vec (:,:,:,i,MOD(axis,3)+1,na,nb)=-tau_blk(MOD(axis+1,3)+1,nb)
              u(p) % vec (:,:,:,i,MOD(axis+1,3)+1,na,nb)=tau_blk(MOD(axis,3)+1,nb)
           enddo
           !
        enddo
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           do na=1,nat
              ! These are the 3*3*nat vectors associated with the
              ! three rotational sum rules (0D system - typ. molecule)
              p=p+1
              do nb=1,nat
                 u(p) % vec (:,:,:,i,MOD(j,3)+1,na,nb)=-tau_blk(MOD(j+1,3)+1,nb)
                 u(p) % vec (:,:,:,i,MOD(j+1,3)+1,na,nb)=tau_blk(MOD(j,3)+1,nb)
              enddo
              !
           enddo
        enddo
     enddo
  endif
  !
  if (n.eq.21) then
     !
     !! Please consider citing C. Lin, S. Ponc\'e and N. Marzari, npj Comput Mater 8, 236 (2022)
     !! if asr='all' is used.
     !
     ! Born-Huang invariance conditions
     do i=1,3
        do j=1,3
           do na=1,nat
              ! These are 3*3*nat vectors associated with the three
              ! rotational sum rules (valid for system of any dimension).
              ! These differ from the case of n=6 (zero-dim), where atom b is
              ! in supercell and its position is wrapped to its nearest periodic
              ! image with the correct weight in the case of degeneracy.
              p=p+1
              do nb=1,nat
                 do n1=1,nr1
                    do n2=1,nr2
                       do n3=1,nr3
                          rcell=matmul(at_blk,(/n1,n2,n3/)-(/1,1,1/))
                          r_ws=rcell-tau_blk(:,nb)+tau_blk(:,na)
                          call ws_all(tau,neq,nr1,nr2,nr3,r_ws,at_blk)
                          do ieq=1,neq
                             u(p) % vec (n1,n2,n3,i,MOD(j,3)+1,na,nb)=&
                             u(p) % vec (n1,n2,n3,i,MOD(j,3)+1,na,nb)-(tau(MOD(j+1,3)+1,ieq)-tau_blk(MOD(j+1,3)+1,na))/DBLE(neq)
                             u(p) % vec (n1,n2,n3,i,MOD(j+1,3)+1,na,nb)=&
                             u(p) % vec (n1,n2,n3,i,MOD(j+1,3)+1,na,nb)+(tau(MOD(j,3)+1,ieq)-tau_blk(MOD(j,3)+1,na))/DBLE(neq)
                          enddo  
                       enddo
                    enddo
                 enddo
              enddo
              !call sp1(frc_new, u(p) % vec (:,:,:,:,:,:,:), nr1,nr2,nr3,nat,scal)
              !write(*,*) scal
              !
              where(abs(u(p) % vec)<eps) u(p) % vec=0.d0
           enddo
        enddo
     enddo
     !
     ! Huang conditions
     if (huang) then
     do na=1,nat
        do nb=1,nat
           do n1=1,nr1
              do n2=1,nr2
                 do n3=1,nr3
                    ! These are 15 vectors for vanishing stress tensor
                    rcell=matmul(at_blk,(/n1,n2,n3/)-(/1,1,1/))
                    r_ws=rcell-tau_blk(:,nb)+tau_blk(:,na)
                    call ws_all(tau,neq,nr1,nr2,nr3,r_ws,at_blk)
                    !
                    do ieq=1,neq
                       ! 1 1 1 2, yx
                       huang_set(:,1)=(/1,1,1,2/)
                       u(p+1) % vec (n1,n2,n3,1,1,na,nb)=&
                       u(p+1) % vec (n1,n2,n3,1,1,na,nb)-tau(1,ieq)*tau(2,ieq)/DBLE(neq)
                       u(p+1) % vec (n1,n2,n3,1,2,na,nb)=&
                       u(p+1) % vec (n1,n2,n3,1,2,na,nb)+tau(1,ieq)*tau(1,ieq)/DBLE(neq)
                       ! 1 1 1 3, zx
                       huang_set(:,2)=(/1,1,1,3/)
                       u(p+2) % vec (n1,n2,n3,1,1,na,nb)=&
                       u(p+2) % vec (n1,n2,n3,1,1,na,nb)-tau(1,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+2) % vec (n1,n2,n3,1,3,na,nb)=&
                       u(p+2) % vec (n1,n2,n3,1,3,na,nb)+tau(1,ieq)*tau(1,ieq)/DBLE(neq)
                       ! 1 1 2 2, xx-yy
                       huang_set(:,3)=(/1,1,2,2/)
                       u(p+3) % vec (n1,n2,n3,1,1,na,nb)=&
                       u(p+3) % vec (n1,n2,n3,1,1,na,nb)-tau(2,ieq)*tau(2,ieq)/DBLE(neq)
                       u(p+3) % vec (n1,n2,n3,2,2,na,nb)=&
                       u(p+3) % vec (n1,n2,n3,2,2,na,nb)+tau(1,ieq)*tau(1,ieq)/DBLE(neq)
                       ! 1 1 2 3
                       huang_set(:,4)=(/1,1,2,3/)
                       u(p+4) % vec (n1,n2,n3,1,1,na,nb)=&
                       u(p+4) % vec (n1,n2,n3,1,1,na,nb)-tau(2,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+4) % vec (n1,n2,n3,2,3,na,nb)=&
                       u(p+4) % vec (n1,n2,n3,2,3,na,nb)+tau(1,ieq)*tau(1,ieq)/DBLE(neq)
                       ! 1 1 3 3, xx-zz
                       huang_set(:,5)=(/1,1,3,3/)
                       u(p+5) % vec (n1,n2,n3,1,1,na,nb)=&
                       u(p+5) % vec (n1,n2,n3,1,1,na,nb)-tau(3,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+5) % vec (n1,n2,n3,3,3,na,nb)=&
                       u(p+5) % vec (n1,n2,n3,3,3,na,nb)+tau(1,ieq)*tau(1,ieq)/DBLE(neq)
                       ! 1 2 1 3
                       huang_set(:,6)=(/1,2,1,3/)
                       u(p+6) % vec (n1,n2,n3,1,2,na,nb)=&
                       u(p+6) % vec (n1,n2,n3,1,2,na,nb)-tau(1,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+6) % vec (n1,n2,n3,1,3,na,nb)=&
                       u(p+6) % vec (n1,n2,n3,1,3,na,nb)+tau(1,ieq)*tau(2,ieq)/DBLE(neq)
                       ! 1 2 2 2, xy
                       huang_set(:,7)=(/1,2,2,2/)
                       u(p+7) % vec (n1,n2,n3,1,2,na,nb)=&
                       u(p+7) % vec (n1,n2,n3,1,2,na,nb)-tau(2,ieq)*tau(2,ieq)/DBLE(neq)
                       u(p+7) % vec (n1,n2,n3,2,2,na,nb)=&
                       u(p+7) % vec (n1,n2,n3,2,2,na,nb)+tau(1,ieq)*tau(2,ieq)/DBLE(neq)
                       ! 1 2 2 3
                       huang_set(:,8)=(/1,2,2,3/)
                       u(p+8) % vec (n1,n2,n3,1,2,na,nb)=&
                       u(p+8) % vec (n1,n2,n3,1,2,na,nb)-tau(2,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+8) % vec (n1,n2,n3,2,3,na,nb)=&
                       u(p+8) % vec (n1,n2,n3,2,3,na,nb)+tau(1,ieq)*tau(2,ieq)/DBLE(neq)
                       ! 1 2 3 3
                       huang_set(:,9)=(/1,2,3,3/)
                       u(p+9) % vec (n1,n2,n3,1,2,na,nb)=&
                       u(p+9) % vec (n1,n2,n3,1,2,na,nb)-tau(3,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+9) % vec (n1,n2,n3,3,3,na,nb)=&
                       u(p+9) % vec (n1,n2,n3,3,3,na,nb)+tau(1,ieq)*tau(2,ieq)/DBLE(neq)
                       ! 1 3 2 2
                       huang_set(:,10)=(/1,3,2,2/)
                       u(p+10) % vec (n1,n2,n3,1,3,na,nb)=&
                       u(p+10) % vec (n1,n2,n3,1,3,na,nb)-tau(2,ieq)*tau(2,ieq)/DBLE(neq)
                       u(p+10) % vec (n1,n2,n3,2,2,na,nb)=&
                       u(p+10) % vec (n1,n2,n3,2,2,na,nb)+tau(1,ieq)*tau(3,ieq)/DBLE(neq)
                       ! 1 3 2 3
                       huang_set(:,11)=(/1,3,2,3/)
                       u(p+11) % vec (n1,n2,n3,1,3,na,nb)=&
                       u(p+11) % vec (n1,n2,n3,1,3,na,nb)-tau(2,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+11) % vec (n1,n2,n3,2,3,na,nb)=&
                       u(p+11) % vec (n1,n2,n3,2,3,na,nb)+tau(1,ieq)*tau(3,ieq)/DBLE(neq)
                       ! 1 3 3 3, xz
                       huang_set(:,12)=(/1,3,3,3/)
                       u(p+12) % vec (n1,n2,n3,1,3,na,nb)=&
                       u(p+12) % vec (n1,n2,n3,1,3,na,nb)-tau(3,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+12) % vec (n1,n2,n3,3,3,na,nb)=&
                       u(p+12) % vec (n1,n2,n3,3,3,na,nb)+tau(1,ieq)*tau(3,ieq)/DBLE(neq)
                       ! 2 2 2 3, zy
                       huang_set(:,13)=(/2,2,2,3/)
                       u(p+13) % vec (n1,n2,n3,2,2,na,nb)=&
                       u(p+13) % vec (n1,n2,n3,2,2,na,nb)-tau(2,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+13) % vec (n1,n2,n3,2,3,na,nb)=&
                       u(p+13) % vec (n1,n2,n3,2,3,na,nb)+tau(2,ieq)*tau(2,ieq)/DBLE(neq)
                       ! 2 2 3 3, yy-zz
                       huang_set(:,14)=(/2,2,3,3/)
                       u(p+14) % vec (n1,n2,n3,2,2,na,nb)=&
                       u(p+14) % vec (n1,n2,n3,2,2,na,nb)-tau(3,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+14) % vec (n1,n2,n3,3,3,na,nb)=&
                       u(p+14) % vec (n1,n2,n3,3,3,na,nb)+tau(2,ieq)*tau(2,ieq)/DBLE(neq)
                       ! 2 3 3 3, yz
                       huang_set(:,15)=(/2,3,3,3/)
                       u(p+15) % vec (n1,n2,n3,2,3,na,nb)=&
                       u(p+15) % vec (n1,n2,n3,2,3,na,nb)-tau(3,ieq)*tau(3,ieq)/DBLE(neq)
                       u(p+15) % vec (n1,n2,n3,3,3,na,nb)=&
                       u(p+15) % vec (n1,n2,n3,3,3,na,nb)+tau(2,ieq)*tau(3,ieq)/DBLE(neq)
                    enddo
                    !
                 enddo
              enddo
           enddo
        enddo
     enddo
     !
     do ip=1,15
        where(abs(u(p+ip) % vec)<eps) u(p+ip) % vec=0.d0
        call sp1(u(p+ip) % vec (:,:,:,:,:,:,:), frc_new, nr1,nr2,nr3,nat,scal)
        write(stdout,'(A, 4I4, A, F20.16)') " Huang before: ", huang_set(:,ip), ",   residual stress: ", scal
     enddo
     write(stdout,*)
     p=p+15
     endif
     !
  endif
  !
  allocate (ind_v(9*nat*nat*nr1*nr2*nr3,2,7), v(9*nat*nat*nr1*nr2*nr3,2) )
  m=0
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       ! These are the vectors associated with the symmetry constraints
                       q=1
                       l=1
                       do while((l.le.m).and.(q.ne.0))
                          if ((ind_v(l,1,1).eq.n1).and.(ind_v(l,1,2).eq.n2).and. &
                               (ind_v(l,1,3).eq.n3).and.(ind_v(l,1,4).eq.i).and. &
                               (ind_v(l,1,5).eq.j).and.(ind_v(l,1,6).eq.na).and. &
                               (ind_v(l,1,7).eq.nb)) q=0
                          if ((ind_v(l,2,1).eq.n1).and.(ind_v(l,2,2).eq.n2).and. &
                               (ind_v(l,2,3).eq.n3).and.(ind_v(l,2,4).eq.i).and. &
                               (ind_v(l,2,5).eq.j).and.(ind_v(l,2,6).eq.na).and. &
                               (ind_v(l,2,7).eq.nb)) q=0
                          l=l+1
                       enddo
                       if ((n1.eq.MOD(nr1+1-n1,nr1)+1).and.(n2.eq.MOD(nr2+1-n2,nr2)+1) &
                            .and.(n3.eq.MOD(nr3+1-n3,nr3)+1).and.(i.eq.j).and.(na.eq.nb)) q=0
                       if (q.ne.0) then
                          m=m+1
                          ind_v(m,1,1)=n1
                          ind_v(m,1,2)=n2
                          ind_v(m,1,3)=n3
                          ind_v(m,1,4)=i
                          ind_v(m,1,5)=j
                          ind_v(m,1,6)=na
                          ind_v(m,1,7)=nb
                          v(m,1)=1.0d0/DSQRT(2.0d0)
                          ind_v(m,2,1)=MOD(nr1+1-n1,nr1)+1
                          ind_v(m,2,2)=MOD(nr2+1-n2,nr2)+1
                          ind_v(m,2,3)=MOD(nr3+1-n3,nr3)+1
                          ind_v(m,2,4)=j
                          ind_v(m,2,5)=i
                          ind_v(m,2,6)=nb
                          ind_v(m,2,7)=na
                          v(m,2)=-1.0d0/DSQRT(2.0d0)
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  ! Note that the vectors corresponding to symmetry constraints are already
  ! orthonormalized by construction.
  !
  n_less=0
  allocate (w(nr1,nr2,nr3,3,3,nat,nat), x(nr1,nr2,nr3,3,3,nat,nat))
  do k=1,p
     w(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     do l=1,m
        !
        call sp2(x,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
        do r=1,2
           n1=ind_v(l,r,1)
           n2=ind_v(l,r,2)
           n3=ind_v(l,r,3)
           i=ind_v(l,r,4)
           j=ind_v(l,r,5)
           na=ind_v(l,r,6)
           nb=ind_v(l,r,7)
           w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)-scal*v(l,r)
        enddo
     enddo
     if (k.le.(9*nat)) then
        na1=MOD(k,nat)
        if (na1.eq.0) na1=nat
        j1=MOD((k-na1)/nat,3)+1
        i1=MOD((((k-na1)/nat)-j1+1)/3,3)+1
     else
        q=k-9*nat
        if (k.le.(18*nat)) then
           if (n.eq.4) then
              na1=MOD(q,nat)
              if (na1.eq.0) na1=nat
              i1=MOD((q-na1)/nat,3)+1
           else
              na1=MOD(q,nat)
              if (na1.eq.0) na1=nat
              j1=MOD((q-na1)/nat,3)+1
              i1=MOD((((q-na1)/nat)-j1+1)/3,3)+1
           endif
        endif
     endif
     do q=1,k-1
        r=1
        do i_less=1,n_less
           if (u_less(i_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           if (k.le.(18*nat)) then
              call sp3(x,u(q) % vec (:,:,:,:,:,:,:), i1,na1,nr1,nr2,nr3,nat,scal)
           else
              call sp1(x,u(q) % vec (:,:,:,:,:,:,:), nr1,nr2,nr3,nat,scal)
           endif
           w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) - scal * u(q) % vec (:,:,:,:,:,:,:)
        endif
     enddo
     call sp1(w,w,nr1,nr2,nr3,nat,norm2)
     if (norm2.gt.eps) then
        u(k) % vec (:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) / DSQRT(norm2)
     else
        n_less=n_less+1
        u_less(n_less)=k
     endif
  enddo
  !
  ! Projection of the force-constants "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules and symmetry contraints
  !
  w(:,:,:,:,:,:,:)=0.0d0
  do l=1,m
     call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
     do r=1,2
        n1=ind_v(l,r,1)
        n2=ind_v(l,r,2)
        n3=ind_v(l,r,3)
        i=ind_v(l,r,4)
        j=ind_v(l,r,5)
        na=ind_v(l,r,6)
        nb=ind_v(l,r,7)
        w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)+scal*v(l,r)
     enddo
  enddo
  do k=1,p
     r=1
     do i_less=1,n_less
        if (u_less(i_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
        call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
        w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) + scal*u(k)%vec(:,:,:,:,:,:,:)
     endif
     deallocate(u(k) % vec)
  enddo
  !
  ! Final substraction of the former projection to the initial frc, to get
  ! the new "projected" frc
  !
  frc_new(:,:,:,:,:,:,:)=frc_new(:,:,:,:,:,:,:) - w(:,:,:,:,:,:,:)
  if (asr.eq.'all') frc_new(:,:,:,:,:,:,:)=frc_new(:,:,:,:,:,:,:) - frc_lr(:,:,:,:,:,:,:)
  call sp1(w,w,nr1,nr2,nr3,nat,norm2)
  write(stdout,'(" Norm of the difference between old and new force-constants:",&
       &     F25.20/)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection IFC")')
  !do l=1,m
  !  call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("l= ",I8," frc_new|v(l)= ",F15.10)') l,scal
  !enddo
  !do k=1,p
  !  x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
  !  call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," frc_new|u(k)= ",F15.10)') k,scal
  !  deallocate(u(k) % vec)
  !enddo
  !
  frc(:,:,:,:,:,:,:)=frc_new(:,:,:,:,:,:,:)
  !
  deallocate (x, w)
  deallocate (v, ind_v)
  deallocate (frc_new)
  !
  return
end subroutine set_asr
!
!
!----------------------------------------------------------------------
! does the scalar product of two effective charges matrices zeu_u and zeu_v
! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  USE kinds, ONLY: DP
  implicit none
  integer i,j,na,nat
  real(DP) zeu_u(3,3,nat)
  real(DP) zeu_v(3,3,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp_zeu
!
!
!----------------------------------------------------------------------
! does the scalar product of two force-constants matrices u and v (considered as
! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual way)
subroutine sp1(u,v,nr1,nr2,nr3,nat,scal)
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        do nb=1,nat
          do n1=1,nr1
            do n2=1,nr2
              do n3=1,nr3
                scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp1
!
!
!-----------------------------------------------------------------------
! does the scalar product of two force-constants matrices u and v (considered as
! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual way
! but v is coded as explained when defining the vectors corresponding to the
! symmetry constraints
subroutine sp2(u,v,ind_v,nr1,nr2,nr3,nat,scal)
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  integer ind_v(2,7)
  real(DP) v(2)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,2
    scal=scal+u(ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4),ind_v(i,5),ind_v(i,6), &
         ind_v(i,7))*v(i)
  enddo
  !
  return
  !
end subroutine sp2
!
!
!
!-----------------------------------------------------------------------
! like sp1, but in the particular case when u is one of the u(k)%vec
! defined in set_asr (before orthonormalization). In this case most of the
! terms are zero (the ones that are not are characterized by i and na), so
! that a lot of computer time can be saved (during Gram-Schmidt).
subroutine sp3(u,v,i,na,nr1,nr2,nr3,nat,scal)
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do j=1,3
    do nb=1,nat
      do n1=1,nr1
        do n2=1,nr2
          do n3=1,nr3
            scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp3


!adapted from PHonon/PH/matdyn.f90
!-----------------------------------------------------------------------
SUBROUTINE ws_all(tau,neq,nr1,nr2,nr3,r_ws,at_blk)
  !-----------------------------------------------------------------------
  !! Used in set_asr
  !! Determine the nearest periodic image for an atom pair
  !! Calculte the number of degeneracies as weight when
  !! on the surface of the optimal Wigner-Seitz cell
  !
  USE kinds, ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nr1,nr2,nr3
  !! Supercell size
  INTEGER, INTENT(out) :: neq
  !! Weight of the degeneracy
  REAL(DP), INTENT(in) :: r_ws(3)
  !! Original position in the supercell
  REAL(DP), INTENT(in) :: at_blk(3, 3)
  !! Lattice vector of the primitive cell
  REAL(DP), INTENT(out) :: tau(3, 48)
  !! Nearest periodic atomic position 
  !
  ! Local variable
  INTEGER :: n1, n2, n3
  !! Search around the cell to a maximum of 2 around
  REAL(DP) :: rnorm
  !! Distance betweem origin cell and current cell
  REAL(DP) :: eps
  !! Tolerence
  REAL(DP) :: tau_tmp(3)
  !! Position of the current cell 
  REAL(DP) :: at(3,3)
  !! Lattice vector of the supercell
  REAL(DP) :: rmin
  !! Minimal distance
  ! 
  ! Large tolerence to find degeneracy position in case of bad relaxation
  eps = 1.0d-5
  !
  at(:, 1) = at_blk(:, 1) * DBLE(nr1)
  at(:, 2) = at_blk(:, 2) * DBLE(nr2)
  at(:, 3) = at_blk(:, 3) * DBLE(nr3)
  !
  rmin = HUGE(rmin)
  DO n1 = -2,2
    DO n2 = -2,2
      DO n3 = -2,2
        tau_tmp = r_ws + n1 * at(:, 1) + n2 * at(:, 2) + n3 * at(:, 3)
        rnorm = NORM2(tau_tmp)
        IF (ABS(rnorm - rmin) .gt. eps) THEN
          IF (rnorm .lt. rmin) THEN
            neq = 1
            rmin = rnorm
            tau(:, neq) = tau_tmp
          ENDIF
        ELSE
          neq = neq + 1
          tau(:, neq) = tau_tmp
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !
  RETURN
  !
!-----------------------------------------------------------------------
END SUBROUTINE ws_all
!-----------------------------------------------------------------------
