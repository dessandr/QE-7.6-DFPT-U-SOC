!
! Finite-difference check of the DFPT xc kernel in the noncollinear magnetic
! case (SOC+U port, E-field / GGA investigation). Test use only.
!
! Reads a pw.x save directory, sets up the xc kernel as ph.x does (phq_setup:
! compute_ux, setup_dmuxc, setup_dgc) and compares, for a few q = 0 density
! perturbations drho, the kernel response dv_of_drho_xc with the central
! difference [v_xc(rho+h drho) - v_xc(rho-h drho)]/2h of the pw.x potential
! (Richardson h, 2h). Kernel variants:
!   lib  : LR_Modules dv_of_drho_xc as is (gc_lr%gmag is never filled)
!   gmag : the same with gmag = grad m filled in before the call
!   cons : grad of delta(rho +- s|m|)/2 by FFT (dgradcorr_cons.f90), the
!          discretisation used by PW gradcorr
! It also compares v_xc with lsign = .false. (what ph.x keeps in v%of_r and
! vrs: post_xml_init runs v_of_rho before phq_setup calls compute_ux) and
! lsign = .true. (pw.x SCF and the DFPT kernel).
!
! Input namelist &xckfd: prefix, outdir, h (default 1d-3).
!
PROGRAM xckfd_main
  USE kinds,            ONLY : DP
  USE mp_global,        ONLY : mp_startup, mp_global_end
  USE mp_world,         ONLY : world_comm
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_bcast, mp_sum, mp_max
  USE environment,      ONLY : environment_start, environment_end
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE io_files,         ONLY : prefix, tmp_dir
  USE fft_base,         ONLY : dfftp
  USE fft_rho,          ONLY : rho_r2g
  USE fft_types,        ONLY : fft_index_to_3d
  USE gvect,            ONLY : g
  USE cell_base,        ONLY : omega
  USE scf,              ONLY : rho, rho_core, rhog_core, v, scf_type, &
                               create_scf_type, destroy_scf_type
  USE lsda_mod,         ONLY : nspin, starting_magnetization
  USE noncollin_module, ONLY : noncolin, domag, nspin_mag, nspin_gga, &
                               m_loc, angle1, angle2, ux, lsign
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp
  USE uspp,             ONLY : nlcc_any
  USE uspp_param,       ONLY : upf
  USE xc_lib,           ONLY : xclib_dft_is
  USE funct,            ONLY : get_dft_name
  USE qpoint,           ONLY : xq
  USE gc_lr,            ONLY : grho, gmag, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE eqv,              ONLY : dmuxc
  USE dv_of_drho_lr,    ONLY : dv_of_drho_xc
  !
  IMPLICIT NONE
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(LEN=256) :: outdir
  REAL(DP) :: h
  NAMELIST / xckfd / prefix, outdir, h
  !
  INTEGER, PARAMETER :: npert = 4, nkern = 3
  CHARACTER(LEN=8), PARAMETER :: pname(npert) = &
       (/ 'charge  ', 'long    ', 'rot     ', 'rotmod  ' /)
  CHARACTER(LEN=4), PARAMETER :: kname(nkern) = (/ 'lib ', 'gmag', 'cons' /)
  INTEGER :: ios, na, ir, is, ip, ik, i1, i2, i3, nrtot
  LOGICAL :: needwf, offrange, lsign_t
  REAL(DP) :: axis(3), w(3), etxc, vtxc, etxc_f, amag, f
  REAL(DP), ALLOCATABLE :: v0(:,:), v0f(:,:), vp(:,:), vm(:,:), dfd(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: drho(:,:), dv(:,:,:), fd(:,:), ref(:,:)
  TYPE(scf_type) :: rhop
  !
  CALL mp_startup()
  CALL environment_start('XCKFD')
  prefix = 'pwscf'
  outdir = './'
  h = 1.d-3
  IF (ionode) THEN
     READ(5, xckfd, IOSTAT=ios)
     IF (ios /= 0) CALL errore('xckfd', 'reading namelist xckfd', ABS(ios))
  ENDIF
  CALL mp_bcast(prefix, ionode_id, world_comm)
  CALL mp_bcast(outdir, ionode_id, world_comm)
  CALL mp_bcast(h, ionode_id, world_comm)
  tmp_dir = trimcheck(outdir)
  needwf = .false.
  CALL read_file_new(needwf)
  IF (.NOT. (noncolin .AND. domag)) CALL errore('xckfd', 'needs noncolin and domag', 1)
  nrtot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
  WRITE(stdout, '(/5x,"dft = ",a,"  gradient = ",l1,"  nlcc = ",l1,"  nr = ",3i5,"  h = ",es9.2)') &
       TRIM(get_dft_name()), xclib_dft_is('gradient'), ANY(upf(1:ntyp)%nlcc), &
       dfftp%nr1, dfftp%nr2, dfftp%nr3, h
  !
  ! ph.x (phq_setup): m_loc, ux, lsign; then the xc kernel
  nlcc_any = ANY(upf(1:ntyp)%nlcc)
  IF (.NOT. ALLOCATED(m_loc)) ALLOCATE(m_loc(3, nat))
  DO na = 1, nat
     m_loc(1,na) = starting_magnetization(ityp(na)) * SIN(angle1(ityp(na))) * COS(angle2(ityp(na)))
     m_loc(2,na) = starting_magnetization(ityp(na)) * SIN(angle1(ityp(na))) * SIN(angle2(ityp(na)))
     m_loc(3,na) = starting_magnetization(ityp(na)) * COS(angle1(ityp(na)))
  ENDDO
  ux = 0.0_DP
  IF (xclib_dft_is('gradient')) CALL compute_ux(m_loc, ux, nat)
  lsign_t = lsign
  WRITE(stdout, '(5x,"lsign = ",l1,"  ux = ",3f10.6)') lsign, ux
  xq = 0.0_DP
  ALLOCATE(dmuxc(dfftp%nnr, nspin_mag, nspin_mag))
  CALL setup_dmuxc()
  CALL setup_dgc()
  !
  ! ground-state v_xc with lsign as in pw.x / the kernel, and with lsign = .false.
  ALLOCATE(v0(dfftp%nnr, nspin), v0f(dfftp%nnr, nspin), vp(dfftp%nnr, nspin), vm(dfftp%nnr, nspin))
  CALL v_xc(rho, rho_core, rhog_core, etxc, vtxc, v0)
  lsign = .false.
  CALL v_xc(rho, rho_core, rhog_core, etxc_f, vtxc, v0f)
  lsign = lsign_t
  WRITE(stdout, '(/5x,"H0 check: v_xc(lsign=F) - v_xc(lsign=",l1,")")') lsign_t
  WRITE(stdout, '(7x,"E_xc(T) = ",f18.10," Ry  E_xc(F) - E_xc(T) = ",es12.4," Ry")') etxc, etxc_f - etxc
  CALL report_real('   v (charge)', v0f(:,1:1) - v0(:,1:1), v0(:,1:1))
  CALL report_real('   B_xc      ', v0f(:,2:4) - v0(:,2:4), v0(:,2:4))
  WRITE(stdout, '(5x,"B in v%of_r from read_file (ph.x vrs) - B_xc(lsign=",l1,")")') lsign_t
  CALL report_real('   B_xc      ', v%of_r(:,2:4) - v0(:,2:4), v0(:,2:4))
  !
  ! transverse xc kernel K = (B_xc.m^)/|m| of the local-frame functional
  ! (dB_perp = K dm_perp): finite for LDA, ~ vsgga/|m| at the nodes for GGA
  CALL kperp_stats(v0)
  ! nodal surfaces: with lsign, s|m| (s = sign(m.ux)) jumps by 2|m_perp| where
  ! m.ux changes sign while the canting m_perp = m - (m.ux)ux does not vanish
  IF (lsign) CALL node_stats()
  ! perturbation axis: perpendicular to the quantization axis (or to the
  ! dominant moment direction when no axis is fixed, e.g. LDA)
  IF (.NOT. lsign) THEN
     w = (/ SUM(rho%of_r(:,2)**2), SUM(rho%of_r(:,3)**2), SUM(rho%of_r(:,4)**2) /)
     CALL mp_sum(w, intra_bgrp_comm)
     ux = 0.0_DP; ux(MAXLOC(w, 1)) = 1.0_DP
  ENDIF
  IF (ABS(ux(3)) > 0.9_DP) THEN
     axis = (/ 1.0_DP, 0.0_DP, 0.0_DP /)
  ELSE
     axis = (/ 0.0_DP, 0.0_DP, 1.0_DP /)
  ENDIF
  WRITE(stdout, '(5x,"rotation axis a = ",3f6.2)') axis
  !
  ALLOCATE(drho(dfftp%nnr, nspin_mag), dv(dfftp%nnr, nspin_mag, nkern), fd(dfftp%nnr, nspin_mag), &
           ref(dfftp%nnr, nspin_mag), dfd(dfftp%nnr, nspin_mag, 2))
  CALL create_scf_type(rhop)
  !
  DO ip = 1, npert
     drho = (0.0_DP, 0.0_DP)
     SELECT CASE (ip)
     CASE (1)
        drho(:,1) = rho%of_r(:,1)
     CASE (2)
        drho(:,2:4) = rho%of_r(:,2:4)
     CASE (3, 4)
        DO ir = 1, dfftp%nnr
           f = 1.0_DP
           IF (ip == 4) THEN
              CALL fft_index_to_3d(ir, dfftp, i1, i2, i3, offrange)
              IF (offrange) CYCLE
              f = COS(8.0_DP * ATAN(1.0_DP) * DBLE(i1) / DBLE(dfftp%nr1))
           ENDIF
           drho(ir,2) = f * (axis(2)*rho%of_r(ir,4) - axis(3)*rho%of_r(ir,3))
           drho(ir,3) = f * (axis(3)*rho%of_r(ir,2) - axis(1)*rho%of_r(ir,4))
           drho(ir,4) = f * (axis(1)*rho%of_r(ir,3) - axis(2)*rho%of_r(ir,2))
        ENDDO
     END SELECT
     !
     ! finite differences of the pw.x v_xc, step h and 2h
     DO is = 1, 2
        CALL vxc_at(DBLE(is) * h, vp)
        CALL vxc_at(-DBLE(is) * h, vm)
        dfd(:,:,is) = (vp - vm) / (2.0_DP * DBLE(is) * h)
     ENDDO
     fd = CMPLX((4.0_DP * dfd(:,:,1) - dfd(:,:,2)) / 3.0_DP, 0.0_DP, KIND=DP)
     !
     DO ik = 1, nkern
        CALL kernel(ik, drho, dv(:,:,ik))
     ENDDO
     !
     WRITE(stdout, '(/5x,"== perturbation ",a," ==   rel. L2 error over (n, mx, my, mz) and per component;", &
          &" q = int drho.dv")') TRIM(pname(ip))
     CALL report_cmplx('FD h vs 2h', CMPLX(dfd(:,:,1) - dfd(:,:,2), 0.0_DP, KIND=DP), fd)
     CALL report_cmplx('FD (Rich.)', fd - fd, fd)
     DO ik = 1, nkern
        CALL report_cmplx('DFPT '//kname(ik), dv(:,:,ik) - fd, fd)
     ENDDO
     IF (ip == 3) THEN
        ! global rotation: exact answer a x B_xc (rotation invariance)
        ref = (0.0_DP, 0.0_DP)
        ref(:,2) = axis(2)*v0(:,4) - axis(3)*v0(:,3)
        ref(:,3) = axis(3)*v0(:,2) - axis(1)*v0(:,4)
        ref(:,4) = axis(1)*v0(:,3) - axis(2)*v0(:,2)
        WRITE(stdout, '(7x,"vs exact a x B_xc(lsign=",l1,"):")') lsign_t
        CALL report_cmplx('FD (Rich.)', fd - ref, ref)
        DO ik = 1, nkern
           CALL report_cmplx('DFPT '//kname(ik), dv(:,:,ik) - ref, ref)
        ENDDO
        ref(:,2) = axis(2)*v0f(:,4) - axis(3)*v0f(:,3)
        ref(:,3) = axis(3)*v0f(:,2) - axis(1)*v0f(:,4)
        ref(:,4) = axis(1)*v0f(:,3) - axis(2)*v0f(:,2)
        WRITE(stdout, '(7x,"vs a x B_xc(lsign=F) (the B_xc in ph.x vrs):")')
        DO ik = 1, nkern
           CALL report_cmplx('DFPT '//kname(ik), dv(:,:,ik) - ref, ref)
        ENDDO
     ENDIF
  ENDDO
  !
  CALL destroy_scf_type(rhop)
  CALL environment_end('XCKFD')
  CALL mp_global_end()
  !
CONTAINS
  !
  SUBROUTINE kperp_stats(vv)
    ! K = (B.m)/|m|^2 on the grid: range, B.m > 0 (anti-Stoner, QE sign:
    ! B.m < 0 for exchange along m) and the rho-weighted share of |K| > 10
    ! max|K_LDA-like| where |m| < 1e-3
    REAL(DP), INTENT(IN) :: vv(dfftp%nnr, nspin)
    REAL(DP) :: am2, k, kmin, kmax, kmax_big, wbig, wtot, wpos, kref
    INTEGER :: ir, npos, nsmall, nbig
    kmin = HUGE(1.0_DP); kmax = -HUGE(1.0_DP); npos = 0; nsmall = 0
    wpos = 0.0_DP; wtot = 0.0_DP
    ! reference scale: K over points with |m| > 0.05 (inside the moments)
    kref = 0.0_DP
    DO ir = 1, dfftp%nnr
       am2 = SUM(rho%of_r(ir,2:4)**2)
       IF (am2 > 0.05_DP**2) kref = MAX(kref, ABS(SUM(vv(ir,2:4) * rho%of_r(ir,2:4))) / am2)
    ENDDO
    nbig = 0; wbig = 0.0_DP; kmax_big = 0.0_DP
    DO ir = 1, dfftp%nnr
       am2 = SUM(rho%of_r(ir,2:4)**2)
       IF (am2 < 1.d-24) CYCLE
       k = SUM(vv(ir,2:4) * rho%of_r(ir,2:4)) / am2
       kmin = MIN(kmin, k); kmax = MAX(kmax, k)
       wtot = wtot + rho%of_r(ir,1)
       IF (k > 0.0_DP) THEN
          npos = npos + 1; wpos = wpos + rho%of_r(ir,1)
       ENDIF
       IF (ABS(k) > 10.0_DP * kref) THEN
          nbig = nbig + 1; wbig = wbig + rho%of_r(ir,1)
       ENDIF
    ENDDO
    WRITE(stdout, '(/5x,"K_perp = (B.m)/|m|^2 [Ry]: min ",es10.3," max ",es10.3,"; |K| scale inside moments ",es10.3)') &
         kmin, kmax, kref
    WRITE(stdout, '(7x,"anti-Stoner points (B.m > 0): ",i8," of ",i8," (rho share ",f7.4,");", &
         &" |K| > 10x scale: ",i7," (rho share ",f7.4,")")') npos, nrtot, wpos / wtot, nbig, wbig / wtot
  END SUBROUTINE kperp_stats
  !
  SUBROUTINE node_stats()
    ! grid points next to a sign change of m.ux (along any grid axis, serial
    ! layout) and the size of the jump of s|m| there
    INTEGER :: ir, jr, i1, i2, i3, j(3), nr(3), d, n_node, n_cant
    REAL(DP) :: mpar, mpar2, mperp, jmax, jrho
    LOGICAL :: off
    n_node = 0; n_cant = 0; jmax = 0.0_DP; jrho = 0.0_DP
    nr = (/ dfftp%nr1, dfftp%nr2, dfftp%nr3 /)
    DO ir = 1, dfftp%nnr
       CALL fft_index_to_3d(ir, dfftp, i1, i2, i3, off)
       IF (off) CYCLE
       mpar = SUM(rho%of_r(ir,2:4) * ux)
       DO d = 1, 3
          j = (/ i1, i2, i3 /); j(d) = MOD(j(d) + 1, nr(d))
          jr = 1 + j(1) + j(2) * dfftp%nr1x + j(3) * dfftp%nr1x * dfftp%my_nr2p
          mpar2 = SUM(rho%of_r(jr,2:4) * ux)
          IF (mpar * mpar2 < 0.0_DP) THEN
             n_node = n_node + 1
             mperp = SQRT(MAX(SUM(rho%of_r(ir,2:4)**2) - mpar**2, 0.0_DP))
             IF (mperp > ABS(mpar)) n_cant = n_cant + 1
             IF (mperp > jmax) THEN
                jmax = mperp; jrho = rho%of_r(ir,1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    WRITE(stdout, '(5x,"nodal links (m.ux sign change) = ",i8," of ",i9,", with |m_perp|>|m.ux|: ",i7, &
         &"; max |m_perp| at a node = ",es10.3," (rho there ",es10.3,", max |m| ",es10.3,")")') &
         n_node, 3 * nrtot, n_cant, jmax, jrho, SQRT(MAXVAL(SUM(rho%of_r(:,2:4)**2, 2)))
  END SUBROUTINE node_stats
  !
  SUBROUTINE vxc_at(t, v)
    ! v_xc(rho + t drho) with the lsign of the kernel
    REAL(DP), INTENT(IN) :: t
    REAL(DP), INTENT(OUT) :: v(dfftp%nnr, nspin)
    REAL(DP) :: e, vt
    rhop%of_r(:,:) = rho%of_r(:,:) + t * DBLE(drho(:,:))
    CALL rho_r2g(dfftp, rhop%of_r, rhop%of_g)
    CALL v_xc(rhop, rho_core, rhog_core, e, vt, v)
  END SUBROUTINE vxc_at
  !
  SUBROUTINE kernel(ik, drho, dv)
    INTEGER, INTENT(IN) :: ik
    COMPLEX(DP), INTENT(IN) :: drho(dfftp%nnr, nspin_mag)
    COMPLEX(DP), INTENT(OUT) :: dv(dfftp%nnr, nspin_mag)
    INTEGER :: is, js
    dv = (0.0_DP, 0.0_DP)
    SELECT CASE (ik)
    CASE (1)
       CALL dv_of_drho_xc(dv, drho=drho)
    CASE (2)
       IF (ALLOCATED(gmag)) THEN
          DO is = 2, nspin_mag
             CALL fft_gradient_g2r(dfftp, rho%of_g(1,is), g, gmag(1,1,is))
          ENDDO
       ENDIF
       CALL dv_of_drho_xc(dv, drho=drho)
       IF (ALLOCATED(gmag)) gmag = 0.0_DP
    CASE (3)
       DO js = 1, nspin_mag
          DO is = 1, nspin_mag
             dv(:,is) = dv(:,is) + dmuxc(:,is,js) * drho(:,js)
          ENDDO
       ENDDO
       IF (xclib_dft_is('gradient')) THEN
          IF (nlcc_any) rho%of_r(:,1) = rho%of_r(:,1) + rho_core(:)
          CALL dgradcorr_cons(dfftp, rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, &
                              xq, drho, nspin_mag, nspin_gga, g, dv)
          IF (nlcc_any) rho%of_r(:,1) = rho%of_r(:,1) - rho_core(:)
       ENDIF
    END SELECT
  END SUBROUTINE kernel
  !
  SUBROUTINE report_cmplx(label, d, r)
    ! ||d|| / ||r|| over all components, per component, and the quadratic
    ! form q = int drho.(r + d) dr (the second derivative of E_xc along drho)
    CHARACTER(LEN=*), INTENT(IN) :: label
    COMPLEX(DP), INTENT(IN) :: d(:,:), r(:,:)
    REAL(DP) :: nd(SIZE(d,2)), nr(SIZE(d,2)), q, dmax, rmax
    INTEGER :: is
    DO is = 1, SIZE(d,2)
       nd(is) = SUM(ABS(d(:,is))**2)
       nr(is) = SUM(ABS(r(:,is))**2)
    ENDDO
    q = SUM(DBLE(CONJG(drho(:,1:SIZE(d,2))) * (r + d))) * omega / DBLE(nrtot)
    dmax = MAXVAL(ABS(d)); rmax = MAXVAL(ABS(r))
    CALL mp_sum(nd, intra_bgrp_comm); CALL mp_sum(nr, intra_bgrp_comm)
    CALL mp_sum(q, intra_bgrp_comm)
    CALL mp_max(dmax, intra_bgrp_comm); CALL mp_max(rmax, intra_bgrp_comm)
    WRITE(stdout, '(7x,a12," rel=",es10.3,"  comp=",4es10.3,"  max|d|/max|r|=",es10.3,"  q=",es16.8)') &
         label, SQRT(SUM(nd) / MAX(SUM(nr), 1.d-300)), SQRT(nd / MAX(SUM(nr), 1.d-300)), &
         dmax / MAX(rmax, 1.d-300), q
  END SUBROUTINE report_cmplx
  !
  SUBROUTINE report_real(label, d, r)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(DP), INTENT(IN) :: d(:,:), r(:,:)
    REAL(DP) :: nd, nr, dmax, rmax
    INTEGER :: imax(2)
    nd = SUM(d**2); nr = SUM(r**2)
    dmax = MAXVAL(ABS(d)); rmax = MAXVAL(ABS(r))
    imax = MAXLOC(ABS(d))
    CALL mp_sum(nd, intra_bgrp_comm); CALL mp_sum(nr, intra_bgrp_comm)
    WRITE(stdout, '(7x,a," rel L2=",es10.3,"  max|d|=",es10.3," Ry (max|r|=",es10.3,")", &
         &"  at rho=",es10.3," |m|=",es10.3)') label, SQRT(nd / MAX(nr, 1.d-300)), dmax, rmax, &
         rho%of_r(imax(1),1), SQRT(SUM(rho%of_r(imax(1),2:4)**2))
  END SUBROUTINE report_real
  !
END PROGRAM xckfd_main
