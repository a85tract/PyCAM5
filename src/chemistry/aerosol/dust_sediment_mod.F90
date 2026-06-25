module dust_sediment_mod

!---------------------------------------------------------------------------------
! Purpose:
!
! Contains routines to compute tendencies from sedimentation of dust
!
! Author: Phil Rasch
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,      only: r8=>shr_kind_r8
  use ppgrid,            only: pcols, pver, pverp
  use physconst,         only: gravit, rair
  use cam_logfile,       only: iulog
  use cam_abortutils,    only: endrun
  use spmd_utils,        only: masterproc

  private
  public :: dust_sediment_vel, dust_sediment_tend


  real (r8), parameter :: vland  = 2.8_r8            ! dust fall velocity over land  (cm/s)
  real (r8), parameter :: vocean = 1.5_r8            ! dust fall velocity over ocean (cm/s)
  real (r8), parameter :: mxsedfac   = 0.99_r8       ! maximum sedimentation flux factor
  logical, save :: dust_sediment_tend_use_native_impl = .false.
  logical, save :: dust_sediment_tend_impl_selected = .false.
  logical, save :: getflx_use_native_impl = .false.
  logical, save :: getflx_impl_selected = .false.
  logical, save :: cfint2_use_native_impl = .false.
  logical, save :: cfint2_impl_selected = .false.
  logical, save :: cfdotmc_pro_use_native_impl = .false.
  logical, save :: cfdotmc_pro_impl_selected = .false.

contains

!===============================================================================
  subroutine dust_sediment_vel (ncol,                               &
       icefrac , landfrac, ocnfrac , pmid    , pdel    , t       , &
       dustmr  , pvdust   )

!----------------------------------------------------------------------

! Compute gravitational sedimentation velocities for dust

    implicit none

! Arguments
    integer, intent(in) :: ncol                     ! number of colums to process

    real(r8), intent(in)  :: icefrac (pcols)        ! sea ice fraction (fraction)
    real(r8), intent(in)  :: landfrac(pcols)        ! land fraction (fraction)
    real(r8), intent(in)  :: ocnfrac (pcols)        ! ocean fraction (fraction)
    real(r8), intent(in)  :: pmid  (pcols,pver)     ! pressure of midpoint levels (Pa)
    real(r8), intent(in)  :: pdel  (pcols,pver)     ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: t     (pcols,pver)     ! temperature (K)
    real(r8), intent(in)  :: dustmr(pcols,pver)     ! dust (kg/kg)

    real(r8), intent(out) :: pvdust (pcols,pverp)    ! vertical velocity of dust (Pa/s)
! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

! Local variables
    real (r8) :: rho(pcols,pver)                    ! air density in kg/m3
    real (r8) :: vfall(pcols)                       ! settling velocity of dust particles (m/s)

    integer i,k

    real (r8) :: lbound, ac, bc, cc

!-----------------------------------------------------------------------
!--------------------- dust fall velocity ----------------------------
!-----------------------------------------------------------------------

    do k = 1,pver
       do i = 1,ncol

          ! merge the dust fall velocities for land and ocean (cm/s)
          ! SHOULD ALSO ACCOUNT FOR ICEFRAC
          vfall(i) = vland*landfrac(i) + vocean*(1._r8-landfrac(i))
          !!         vfall(i) = vland*landfrac(i) + vocean*ocnfrac(i) + vseaice*icefrac(i)

          ! fall velocity (assume positive downward)
          pvdust(i,k+1) = vfall(i)     
       end do
    end do

    return
  end subroutine dust_sediment_vel


!===============================================================================
  subroutine dust_sediment_tend ( &
       ncol,   dtime,  pint,     pmid,    pdel,  t,   &
       dustmr ,pvdust, dusttend, sfdust )

!----------------------------------------------------------------------
!     Apply Particle Gravitational Sedimentation 
!----------------------------------------------------------------------

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

! Arguments
    integer,  intent(in)  :: ncol                      ! number of colums to process

    real(r8), intent(in)  :: dtime                     ! time step
    real(r8), target, intent(in)  :: pint  (pcols,pverp)       ! interfaces pressure (Pa)
    real(r8), target, intent(in)  :: pmid  (pcols,pver)        ! midpoint pressures (Pa)
    real(r8), target, intent(in)  :: pdel  (pcols,pver)        ! pressure diff across layer (Pa)
    real(r8), target, intent(in)  :: t     (pcols,pver)        ! temperature (K)
    real(r8), target, intent(in)  :: dustmr(pcols,pver)        ! dust (kg/kg)
    real(r8), target, intent(in)  :: pvdust (pcols,pverp)      ! vertical velocity of dust drops  (Pa/s)
! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

    real(r8), target, intent(out) :: dusttend(pcols,pver)      ! dust tend
    real(r8), target, intent(out) :: sfdust  (pcols)           ! surface flux of dust (rain, kg/m/s)

! Local variables
    real(r8), target :: fxdust(pcols,pverp)             ! fluxes at the interfaces, dust (positive = down)
    real(r8), target :: psi(pcols,pverp)
    real(r8), target :: fdot(pcols,pverp)
    real(r8), target :: xxk(pcols,pver)
    real(r8), target :: fxdot(pcols)
    real(r8), target :: fxdd(pcols)
    real(r8), target :: psistar(pcols)
    real(r8), target :: xins(pcols)
    real(r8), target :: s(pcols,pverp)
    real(r8), target :: sh(pcols,pverp)
    real(r8), target :: d(pcols,pverp)
    real(r8), target :: dh(pcols,pverp)
    real(r8), target :: e(pcols,pverp)
    real(r8), target :: eh(pcols,pverp)
    real(r8), target :: ppl(pcols,pverp)
    real(r8), target :: ppr(pcols,pverp)
    real(r8), target :: delxh(pcols,pverp)
    integer(c_int64_t), target :: intz(pcols)
    integer(c_int64_t), target :: status_code
    integer(c_int64_t), target :: fail_i
    integer(c_int64_t), target :: fail_k

    interface
       subroutine dust_sediment_tend_codon(ncol_c, pcols_c, pver_c, pverp_c, dtime_c, mxsedfac_c, gravit_c, &
            pint_p, pdel_p, dustmr_p, pvdust_p, dusttend_p, sfdust_p, fxdust_p, psi_p, fdot_p, xxk_p, fxdot_p, &
            fxdd_p, psistar_p, s_p, sh_p, d_p, dh_p, e_p, eh_p, ppl_p, ppr_p, delxh_p, xins_p, intz_p, &
            status_p, fail_i_p, fail_k_p) bind(c, name="dust_sediment_tend_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c
         real(c_double), value :: dtime_c, mxsedfac_c, gravit_c
         type(c_ptr), value :: pint_p, pdel_p, dustmr_p, pvdust_p, dusttend_p, sfdust_p, fxdust_p, psi_p
         type(c_ptr), value :: fdot_p, xxk_p, fxdot_p, fxdd_p, psistar_p, s_p, sh_p, d_p, dh_p, e_p, eh_p
         type(c_ptr), value :: ppl_p, ppr_p, delxh_p, xins_p, intz_p, status_p, fail_i_p, fail_k_p
       end subroutine dust_sediment_tend_codon
    end interface

!----------------------------------------------------------------------

    call dust_sediment_tend_select_impl()

    if (dust_sediment_tend_use_native_impl) then
       call dust_sediment_tend_native(ncol, dtime, pint, pmid, pdel, t, dustmr, pvdust, dusttend, sfdust)
       return
    end if

    status_code = 0_c_int64_t
    fail_i = 0_c_int64_t
    fail_k = 0_c_int64_t

    call dust_sediment_tend_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), &
         real(dtime, c_double), real(mxsedfac, c_double), real(gravit, c_double), &
         c_loc(pint), c_loc(pdel), c_loc(dustmr), c_loc(pvdust), c_loc(dusttend), c_loc(sfdust), c_loc(fxdust), &
         c_loc(psi), c_loc(fdot), c_loc(xxk), c_loc(fxdot), c_loc(fxdd), c_loc(psistar), c_loc(s), c_loc(sh), &
         c_loc(d), c_loc(dh), c_loc(e), c_loc(eh), c_loc(ppl), c_loc(ppr), c_loc(delxh), c_loc(xins), c_loc(intz), &
         c_loc(status_code), c_loc(fail_i), c_loc(fail_k) &
    )

    if (status_code /= 0_c_int64_t) then
       write(iulog,*) 'DUST_SEDIMENT_MOD:dust_sediment_tend -- interval was not found ', int(fail_i), int(fail_k)
       call endrun('DUST_SEDIMENT_MOD:dust_sediment_tend -- interval was not found ')
    end if

    return
  end subroutine dust_sediment_tend

!===============================================================================
  subroutine dust_sediment_tend_select_impl()

    implicit none

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (dust_sediment_tend_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('DUST_SEDIMENT_TEND_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       dust_sediment_tend_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       dust_sediment_tend_use_native_impl = .false.
    end if

    dust_sediment_tend_impl_selected = .true.

    if (masterproc) then
       if (dust_sediment_tend_use_native_impl) then
          write(iulog,*) 'dust_sediment_tend implementation = native'
       else
          write(iulog,*) 'dust_sediment_tend implementation = codon'
       end if
    end if

  end subroutine dust_sediment_tend_select_impl

!===============================================================================
  subroutine getflx_select_impl()

    implicit none

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (getflx_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('DUST_GETFLX_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       getflx_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       getflx_use_native_impl = .false.
    end if

    getflx_impl_selected = .true.

    if (masterproc) then
       if (getflx_use_native_impl) then
          write(iulog,*) 'getflx implementation = native'
       else
          write(iulog,*) 'getflx implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine getflx_select_impl

!===============================================================================
  subroutine cfint2_select_impl()

    implicit none

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (cfint2_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('DUST_CFINT2_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       cfint2_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       cfint2_use_native_impl = .false.
    end if

    cfint2_impl_selected = .true.

    if (masterproc) then
       if (cfint2_use_native_impl) then
          write(iulog,*) 'cfint2 implementation = native'
       else
          write(iulog,*) 'cfint2 implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine cfint2_select_impl

!===============================================================================
  subroutine cfdotmc_pro_select_impl()

    implicit none

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (cfdotmc_pro_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('DUST_CFDOTMC_PRO_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       cfdotmc_pro_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       cfdotmc_pro_use_native_impl = .false.
    end if

    cfdotmc_pro_impl_selected = .true.

    if (masterproc) then
       if (cfdotmc_pro_use_native_impl) then
          write(iulog,*) 'cfdotmc_pro implementation = native'
       else
          write(iulog,*) 'cfdotmc_pro implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine cfdotmc_pro_select_impl

!===============================================================================
  subroutine dust_sediment_tend_native ( &
       ncol,   dtime,  pint,     pmid,    pdel,  t,   &
       dustmr ,pvdust, dusttend, sfdust )

!----------------------------------------------------------------------
!     Apply Particle Gravitational Sedimentation
!----------------------------------------------------------------------

    implicit none

! Arguments
    integer,  intent(in)  :: ncol                      ! number of colums to process

    real(r8), intent(in)  :: dtime                     ! time step
    real(r8), intent(in)  :: pint  (pcols,pverp)       ! interfaces pressure (Pa)
    real(r8), intent(in)  :: pmid  (pcols,pver)        ! midpoint pressures (Pa)
    real(r8), intent(in)  :: pdel  (pcols,pver)        ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: t     (pcols,pver)        ! temperature (K)
    real(r8), intent(in)  :: dustmr(pcols,pver)        ! dust (kg/kg)
    real(r8), intent(in)  :: pvdust (pcols,pverp)      ! vertical velocity of dust drops  (Pa/s)
! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

    real(r8), intent(out) :: dusttend(pcols,pver)      ! dust tend
    real(r8), intent(out) :: sfdust  (pcols)           ! surface flux of dust (rain, kg/m/s)

! Local variables
    real(r8) :: fxdust(pcols,pverp)                     ! fluxes at the interfaces, dust (positive = down)

    integer :: i,k
!----------------------------------------------------------------------

! initialize variables
    fxdust  (:ncol,:) = 0._r8 ! flux at interfaces (dust)
    dusttend(:ncol,:) = 0._r8 ! tend (dust)
    sfdust(:ncol)     = 0._r8 ! sedimentation flux out bot of column (dust)

! fluxes at interior points
    call getflx(ncol, pint, dustmr, pvdust, dtime, fxdust)

! calculate fluxes at boundaries
    do i = 1,ncol
       fxdust(i,1) = 0
! surface flux by upstream scheme
       fxdust(i,pverp) = dustmr(i,pver) * pvdust(i,pverp) * dtime
    end do

! filter out any negative fluxes from the getflx routine
    do k = 2,pver
       fxdust(:ncol,k) = max(0._r8, fxdust(:ncol,k))
    end do

! Limit the flux out of the bottom of each cell to the water content in each phase.
! Apply mxsedfac to prevent generating very small negative cloud water/ice
! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
! ***Should we include the flux in the top, to allow for thin surface layers?
! ***Requires simple treatment of cloud overlap, already included below.
    do k = 1,pver
       do i = 1,ncol
          fxdust(i,k+1) = min( fxdust(i,k+1), mxsedfac * dustmr(i,k) * pdel(i,k) )
!!$        fxdust(i,k+1) = min( fxdust(i,k+1), dustmr(i,k) * pdel(i,k) + fxdust(i,k))
       end do
    end do

! Now calculate the tendencies 
    do k = 1,pver
       do i = 1,ncol
! net flux into cloud changes cloud dust/ice (all flux is out of cloud)
          dusttend(i,k)  = (fxdust(i,k) - fxdust(i,k+1)) / (dtime * pdel(i,k))
       end do
    end do

! convert flux out the bottom to mass units Pa -> kg/m2/s
    sfdust(:ncol) = fxdust(:ncol,pverp) / (dtime*gravit)

    return
  end subroutine dust_sediment_tend_native

!===============================================================================
  subroutine getflx(ncol, xw, phi, vel, deltat, flux)
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

!.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
!....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
!....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
!.........phi1......phi2.......phi3.....phi4.......phi5.......


    implicit none

    integer ncol                      ! number of colums to process

    integer i
    integer k

    real (r8), target :: vel(pcols,pverp)
    real (r8), target :: flux(pcols,pverp)
    real (r8), target :: xw(pcols,pverp)
    real (r8), target :: psi(pcols,pverp)
    real (r8), target :: phi(pcols,pverp-1)
    real (r8), target :: fdot(pcols,pverp)
    real (r8) :: xx(pcols)
    real (r8), target :: fxdot(pcols)
    real (r8), target :: fxdd(pcols)

    real (r8), target :: psistar(pcols)
    real (r8) deltat

    real (r8), target :: xxk(pcols,pver)
    real (r8), target :: xins(pcols)
    real (r8), target :: s(pcols,pverp)
    real (r8), target :: sh(pcols,pverp)
    real (r8), target :: d(pcols,pverp)
    real (r8), target :: dh(pcols,pverp)
    real (r8), target :: e(pcols,pverp)
    real (r8), target :: eh(pcols,pverp)
    real (r8), target :: ppl(pcols,pverp)
    real (r8), target :: ppr(pcols,pverp)
    real (r8), target :: delxh(pcols,pverp)
    integer(c_int64_t), target :: intz(pcols)
    integer(c_int64_t), target :: status_code
    integer(c_int64_t), target :: fail_i
    integer(c_int64_t), target :: fail_k

    interface
       subroutine getflx_codon(ncol_c, pcols_c, pver_c, pverp_c, deltat_c, xw_p, phi_p, vel_p, flux_p, &
            psi_p, fdot_p, xxk_p, fxdot_p, fxdd_p, psistar_p, xins_p, intz_p, status_p, fail_i_p, fail_k_p, &
            s_p, sh_p, d_p, dh_p, e_p, eh_p, ppl_p, ppr_p, delxh_p) bind(c, name="getflx_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c
         real(c_double), value :: deltat_c
         type(c_ptr), value :: xw_p, phi_p, vel_p, flux_p, psi_p, fdot_p, xxk_p, fxdot_p, fxdd_p, psistar_p
         type(c_ptr), value :: xins_p, intz_p, status_p, fail_i_p, fail_k_p
         type(c_ptr), value :: s_p, sh_p, d_p, dh_p, e_p, eh_p, ppl_p, ppr_p, delxh_p
       end subroutine getflx_codon
    end interface

    call getflx_select_impl()
    if (.not. getflx_use_native_impl) then
       status_code = 0_c_int64_t
       fail_i = 0_c_int64_t
       fail_k = 0_c_int64_t
       call getflx_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), &
            real(deltat, c_double), c_loc(xw), c_loc(phi), c_loc(vel), c_loc(flux), c_loc(psi), c_loc(fdot), &
            c_loc(xxk), c_loc(fxdot), c_loc(fxdd), c_loc(psistar), c_loc(xins), c_loc(intz), c_loc(status_code), &
            c_loc(fail_i), c_loc(fail_k), c_loc(s), c_loc(sh), c_loc(d), c_loc(dh), c_loc(e), c_loc(eh), &
            c_loc(ppl), c_loc(ppr), c_loc(delxh) &
       )
       if (status_code /= 0_c_int64_t) then
          write(iulog,*) 'DUST_SEDIMENT_MOD:getflx -- interval was not found ', int(fail_i), int(fail_k)
          call endrun('DUST_SEDIMENT_MOD:getflx -- interval was not found ')
       end if
       return
    end if

    do i = 1,ncol
!        integral of phi
       psi(i,1) = 0._r8
!        fluxes at boundaries
       flux(i,1) = 0
       flux(i,pverp) = 0._r8
    end do

!     integral function
    do k = 2,pverp
       do i = 1,ncol
          psi(i,k) = phi(i,k-1)*(xw(i,k)-xw(i,k-1)) + psi(i,k-1)
       end do
    end do


!     calculate the derivatives for the interpolating polynomial
    call cfdotmc_pro (ncol, xw, psi, fdot)

!  NEW WAY
!     calculate fluxes at interior pts
    do k = 2,pver
       do i = 1,ncol
          xxk(i,k) = xw(i,k)-vel(i,k)*deltat
       end do
    end do
    do k = 2,pver
       call cfint2(ncol, xw, psi, fdot, xxk(1,k), fxdot, fxdd, psistar)
       do i = 1,ncol
          flux(i,k) = (psi(i,k)-psistar(i))
       end do
    end do


    return
  end subroutine getflx



!##############################################################################

  subroutine cfint2 (ncol, x, f, fdot, xin, fxdot, fxdd, psistar)
    use iso_c_binding, only: c_int64_t, c_loc, c_ptr


    implicit none

! input
    integer ncol                      ! number of colums to process

    real (r8), target :: x(pcols, pverp)
    real (r8), target :: f(pcols, pverp)
    real (r8), target :: fdot(pcols, pverp)
    real (r8), target :: xin(pcols)

! output
    real (r8), target :: fxdot(pcols)
    real (r8), target :: fxdd(pcols)
    real (r8), target :: psistar(pcols)

    integer i
    integer k
    integer(c_int64_t), target :: intz(pcols)
    integer(c_int64_t), target :: status_code
    integer(c_int64_t), target :: fail_i
    real (r8) dx
    real (r8) s
    real (r8) c2
    real (r8) c3
    real (r8) xx
    real (r8) xinf
    real (r8) psi1, psi2, psi3, psim
    real (r8) cfint
    real (r8) cfnew
    real (r8), target :: xins(pcols)

    interface
       subroutine cfint2_codon(ncol_c, pcols_c, pverp_c, x_p, f_p, fdot_p, xin_p, fxdot_p, fxdd_p, &
            psistar_p, xins_p, intz_p, status_p, fail_i_p) bind(c, name="cfint2_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pverp_c
         type(c_ptr), value :: x_p, f_p, fdot_p, xin_p, fxdot_p, fxdd_p, psistar_p, xins_p, intz_p
         type(c_ptr), value :: status_p, fail_i_p
       end subroutine cfint2_codon
    end interface

!     the minmod function 
    real (r8) a, b, c
    real (r8) minmod
    real (r8) medan
    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    call cfint2_select_impl()
    if (.not. cfint2_use_native_impl) then
       status_code = 0_c_int64_t
       fail_i = 0_c_int64_t
       call cfint2_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pverp, c_int64_t), &
            c_loc(x), c_loc(f), c_loc(fdot), c_loc(xin), c_loc(fxdot), c_loc(fxdd), c_loc(psistar), &
            c_loc(xins), c_loc(intz), c_loc(status_code), c_loc(fail_i) &
       )
       if (status_code /= 0_c_int64_t) then
          write(iulog,*) 'DUST_SEDIMENT_MOD:cfint2 -- interval was not found ', int(fail_i)
          call endrun('DUST_SEDIMENT_MOD:cfint2 -- interval was not found ')
       end if
       return
    end if

    do i = 1,ncol
       xins(i) = medan(x(i,1), xin(i), x(i,pverp))
       intz(i) = 0
    end do

! first find the interval 
    do k =  1,pverp-1
       do i = 1,ncol
          if ((xins(i)-x(i,k))*(x(i,k+1)-xins(i)).ge.0._r8) then
             intz(i) = k
          endif
       end do
    end do

    do i = 1,ncol
       if (intz(i).eq.0) then
          write(iulog,*) ' interval was not found for col i ', i
          call endrun('DUST_SEDIMENT_MOD:cfint2 -- interval was not found ')
       endif
    end do

! now interpolate
    do i = 1,ncol
       k = intz(i)
       dx = (x(i,k+1)-x(i,k))
       s = (f(i,k+1)-f(i,k))/dx
       c2 = (3*s-2*fdot(i,k)-fdot(i,k+1))/dx
       c3 = (fdot(i,k)+fdot(i,k+1)-2*s)/dx**2
       xx = (xins(i)-x(i,k))
       fxdot(i) =  (3*c3*xx + 2*c2)*xx + fdot(i,k)
       fxdd(i) = 6*c3*xx + 2*c2
       cfint = ((c3*xx + c2)*xx + fdot(i,k))*xx + f(i,k)

!        limit the interpolant
       psi1 = f(i,k)+(f(i,k+1)-f(i,k))*xx/dx
       if (k.eq.1) then
          psi2 = f(i,1)
       else
          psi2 = f(i,k) + (f(i,k)-f(i,k-1))*xx/(x(i,k)-x(i,k-1))
       endif
       if (k+1.eq.pverp) then
          psi3 = f(i,pverp)
       else
          psi3 = f(i,k+1) - (f(i,k+2)-f(i,k+1))*(dx-xx)/(x(i,k+2)-x(i,k+1))
       endif
       psim = medan(psi1, psi2, psi3)
       cfnew = medan(cfint, psi1, psim)
       if (abs(cfnew-cfint)/(abs(cfnew)+abs(cfint)+1.e-36_r8)  .gt..03_r8) then
!     CHANGE THIS BACK LATER!!!
!     $        .gt..1) then


!     UNCOMMENT THIS LATER!!!
!            write(iulog,*) ' cfint2 limiting important ', cfint, cfnew


       endif
       psistar(i) = cfnew
    end do

    return
  end subroutine cfint2



!##############################################################################

  subroutine cfdotmc_pro (ncol, x, f, fdot)
    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

!     prototype version; eventually replace with final SPITFIRE scheme

!     calculate the derivative for the interpolating polynomial
!     multi column version


    implicit none

! input
    integer ncol                      ! number of colums to process

    real (r8), target :: x(pcols, pverp)
    real (r8), target :: f(pcols, pverp)
! output
    real (r8), target :: fdot(pcols, pverp)          ! derivative at nodes

! assumed variable distribution
!     x1.......x2.......x3.......x4.......x5.......x6     1,pverp points
!     f1.......f2.......f3.......f4.......f5.......f6     1,pverp points
!     ...sh1.......sh2......sh3......sh4......sh5....     1,pver points
!     .........d2.......d3.......d4.......d5.........     2,pver points
!     .........s2.......s3.......s4.......s5.........     2,pver points
!     .............dh2......dh3......dh4.............     2,pver-1 points
!     .............eh2......eh3......eh4.............     2,pver-1 points
!     ..................e3.......e4..................     3,pver-1 points
!     .................ppl3......ppl4................     3,pver-1 points
!     .................ppr3......ppr4................     3,pver-1 points
!     .................t3........t4..................     3,pver-1 points
!     ................fdot3.....fdot4................     3,pver-1 points


! work variables


    integer i
    integer k

    real (r8) a                    ! work var
    real (r8) b                    ! work var
    real (r8) c                    ! work var
    real (r8), target :: s(pcols,pverp)             ! first divided differences at nodes
    real (r8), target :: sh(pcols,pverp)            ! first divided differences between nodes
    real (r8), target :: d(pcols,pverp)             ! second divided differences at nodes
    real (r8), target :: dh(pcols,pverp)            ! second divided differences between nodes
    real (r8), target :: e(pcols,pverp)             ! third divided differences at nodes
    real (r8), target :: eh(pcols,pverp)            ! third divided differences between nodes
    real (r8) pp                   ! p prime
    real (r8), target :: ppl(pcols,pverp)           ! p prime on left
    real (r8), target :: ppr(pcols,pverp)           ! p prime on right
    real (r8) qpl
    real (r8) qpr
    real (r8) ttt
    real (r8) t
    real (r8) tmin
    real (r8) tmax
    real (r8), target :: delxh(pcols,pverp)

    interface
       subroutine cfdotmc_pro_codon(ncol_c, pcols_c, pver_c, pverp_c, x_p, f_p, fdot_p, &
            s_p, sh_p, d_p, dh_p, e_p, eh_p, ppl_p, ppr_p, delxh_p) bind(c, name="cfdotmc_pro_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c
         type(c_ptr), value :: x_p, f_p, fdot_p, s_p, sh_p, d_p, dh_p, e_p, eh_p, ppl_p, ppr_p, delxh_p
       end subroutine cfdotmc_pro_codon
    end interface


!     the minmod function 
    real (r8) minmod
    real (r8) medan
    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    call cfdotmc_pro_select_impl()
    if (.not. cfdotmc_pro_use_native_impl) then
       call cfdotmc_pro_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), &
            c_loc(x), c_loc(f), c_loc(fdot), c_loc(s), c_loc(sh), c_loc(d), c_loc(dh), c_loc(e), c_loc(eh), &
            c_loc(ppl), c_loc(ppr), c_loc(delxh) &
       )
       return
    end if

    do k = 1,pver


!        first divided differences between nodes
       do i = 1, ncol
          delxh(i,k) = (x(i,k+1)-x(i,k))
          sh(i,k) = (f(i,k+1)-f(i,k))/delxh(i,k)
       end do

!        first and second divided differences at nodes
       if (k.ge.2) then
          do i = 1,ncol
             d(i,k) = (sh(i,k)-sh(i,k-1))/(x(i,k+1)-x(i,k-1))
             s(i,k) = minmod(sh(i,k),sh(i,k-1))
          end do
       endif
    end do

!     second and third divided diffs between nodes
    do k = 2,pver-1
       do i = 1, ncol
          eh(i,k) = (d(i,k+1)-d(i,k))/(x(i,k+2)-x(i,k-1))
          dh(i,k) = minmod(d(i,k),d(i,k+1))
       end do
    end do

!     treat the boundaries
    do i = 1,ncol
       e(i,2) = eh(i,2)
       e(i,pver) = eh(i,pver-1)
!        outside level
       fdot(i,1) = sh(i,1) - d(i,2)*delxh(i,1)  &
            - eh(i,2)*delxh(i,1)*(x(i,1)-x(i,3))
       fdot(i,1) = minmod(fdot(i,1),3*sh(i,1))
       fdot(i,pverp) = sh(i,pver) + d(i,pver)*delxh(i,pver)  &
            + eh(i,pver-1)*delxh(i,pver)*(x(i,pverp)-x(i,pver-1))
       fdot(i,pverp) = minmod(fdot(i,pverp),3*sh(i,pver))
!        one in from boundary
       fdot(i,2) = sh(i,1) + d(i,2)*delxh(i,1) - eh(i,2)*delxh(i,1)*delxh(i,2)
       fdot(i,2) = minmod(fdot(i,2),3*s(i,2))
       fdot(i,pver) = sh(i,pver) - d(i,pver)*delxh(i,pver)   &
            - eh(i,pver-1)*delxh(i,pver)*delxh(i,pver-1)
       fdot(i,pver) = minmod(fdot(i,pver),3*s(i,pver))
    end do


    do k = 3,pver-1
       do i = 1,ncol
          e(i,k) = minmod(eh(i,k),eh(i,k-1))
       end do
    end do



    do k = 3,pver-1

       do i = 1,ncol

!           p prime at k-0.5
          ppl(i,k)=sh(i,k-1) + dh(i,k-1)*delxh(i,k-1)  
!           p prime at k+0.5
          ppr(i,k)=sh(i,k)   - dh(i,k)  *delxh(i,k)

          t = minmod(ppl(i,k),ppr(i,k))

!           derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
          pp = sh(i,k-1) + d(i,k)*delxh(i,k-1) 

!           quartic estimate of fdot
          fdot(i,k) = pp                            &
               - delxh(i,k-1)*delxh(i,k)            &
               *(  eh(i,k-1)*(x(i,k+2)-x(i,k  ))    &
               + eh(i,k  )*(x(i,k  )-x(i,k-2))      &
               )/(x(i,k+2)-x(i,k-2))

!           now limit it
          qpl = sh(i,k-1)       &
               + delxh(i,k-1)*minmod(d(i,k-1)+e(i,k-1)*(x(i,k)-x(i,k-2)), &
               d(i,k)  -e(i,k)*delxh(i,k))
          qpr = sh(i,k)         &
               + delxh(i,k  )*minmod(d(i,k)  +e(i,k)*delxh(i,k-1),        &
               d(i,k+1)+e(i,k+1)*(x(i,k)-x(i,k+2)))

          fdot(i,k) = medan(fdot(i,k), qpl, qpr)

          ttt = minmod(qpl, qpr)
          tmin = min(0._r8,3*s(i,k),1.5_r8*t,ttt)
          tmax = max(0._r8,3*s(i,k),1.5_r8*t,ttt)

          fdot(i,k) = fdot(i,k) + minmod(tmin-fdot(i,k), tmax-fdot(i,k))

       end do

    end do

    return
  end subroutine cfdotmc_pro
end module dust_sediment_mod
