module ap_dust_sediment_tend_scheme

!---------------------------------------------------------------------------------
! Purpose:
!
! Contains routines to compute tendencies from sedimentation of dust
!
! Author: Phil Rasch
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,      only: r8=>shr_kind_r8

  private
  public :: dust_sediment_vel, dust_sediment_tend_run


  real (r8), parameter :: vland  = 2.8_r8            ! dust fall velocity over land  (cm/s)
  real (r8), parameter :: vocean = 1.5_r8            ! dust fall velocity over ocean (cm/s)
  real (r8), parameter :: mxsedfac   = 0.99_r8       ! maximum sedimentation flux factor
  integer :: pcols_cfg = 0
  integer :: pver_cfg = 0
  integer :: pverp_cfg = 0
  character(len=512) :: core_errmsg = ''
  integer :: core_errflg = 0

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

    real(r8), intent(in)  :: icefrac (:)        ! sea ice fraction (fraction)
    real(r8), intent(in)  :: landfrac(:)        ! land fraction (fraction)
    real(r8), intent(in)  :: ocnfrac (:)        ! ocean fraction (fraction)
    real(r8), intent(in)  :: pmid  (:,:)        ! pressure of midpoint levels (Pa)
    real(r8), intent(in)  :: pdel  (:,:)        ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: t     (:,:)        ! temperature (K)
    real(r8), intent(in)  :: dustmr(:,:)        ! dust (kg/kg)

    real(r8), intent(out) :: pvdust (:,:)        ! vertical velocity of dust (Pa/s)
! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

! Local variables
    real (r8) :: vfall(size(landfrac))              ! settling velocity of dust particles (m/s)

    integer i,k

    real (r8) :: lbound, ac, bc, cc

!-----------------------------------------------------------------------
!--------------------- dust fall velocity ----------------------------
!-----------------------------------------------------------------------

    pcols_cfg = size(pmid, 1)
    pver_cfg = size(pmid, 2)
    pverp_cfg = size(pvdust, 2)

    ! Interface 1 and padded columns are not used by the sedimentation
    ! tendency, but intent(out) requires them to be defined on return.
    pvdust(:,:) = 0._r8

    do k = 1,pver_cfg
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
  !> \section arg_table_dust_sediment_tend_run Argument Table
  !! \htmlinclude dust_sediment_tend_run.html
  subroutine dust_sediment_tend_run ( &
       ncol,   dtime,  pint,     pmid,    pdel,  t,   &
       dustmr ,pvdust, gravit,   dusttend, sfdust, errmsg, errflg )

!----------------------------------------------------------------------
!     Apply Particle Gravitational Sedimentation
!----------------------------------------------------------------------

    implicit none

! Arguments
    integer,  intent(in)  :: ncol                      ! number of colums to process

    real(r8), intent(in)  :: dtime                     ! time step
    real(r8), intent(in)  :: pint  (:,:)       ! interfaces pressure (Pa)
    real(r8), intent(in)  :: pmid  (:,:)        ! midpoint pressures (Pa)
    real(r8), intent(in)  :: pdel  (:,:)        ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: t     (:,:)        ! temperature (K)
    real(r8), intent(in)  :: dustmr(:,:)        ! dust (kg/kg)
    real(r8), intent(in)  :: pvdust (:,:)      ! vertical velocity of dust drops  (Pa/s)
    real(r8), intent(in)  :: gravit            ! gravitational acceleration
! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

    real(r8), intent(out) :: dusttend(:,:)      ! dust tend
    real(r8), intent(out) :: sfdust  (:)           ! surface flux of dust (rain, kg/m/s)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

! Local variables
    real(r8) :: fxdust(size(pint,1),size(pint,2)) ! fluxes at interfaces, positive down

    integer :: i,k
!----------------------------------------------------------------------

    pcols_cfg = size(pint, 1)
    pverp_cfg = size(pint, 2)
    pver_cfg = pverp_cfg - 1
    core_errmsg = ''
    core_errflg = 0
    errmsg = ''
    errflg = 0

! initialize variables
    fxdust  (:,:) = 0._r8 ! flux at interfaces (dust)
    dusttend(:,:) = 0._r8 ! tend (dust)
    sfdust(:)     = 0._r8 ! sedimentation flux out bot of column (dust)

! fluxes at interior points
    call getflx(ncol, pint, dustmr, pvdust, dtime, fxdust)
    if (core_errflg /= 0) then
       errmsg = core_errmsg
       errflg = core_errflg
       return
    end if

! calculate fluxes at boundaries
    do i = 1,ncol
       fxdust(i,1) = 0
! surface flux by upstream scheme
       fxdust(i,pverp_cfg) = dustmr(i,pver_cfg) * pvdust(i,pverp_cfg) * dtime
    end do

! filter out any negative fluxes from the getflx routine
    do k = 2,pver_cfg
       fxdust(:ncol,k) = max(0._r8, fxdust(:ncol,k))
    end do

! Limit the flux out of the bottom of each cell to the water content in each phase.
! Apply mxsedfac to prevent generating very small negative cloud water/ice
! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
! ***Should we include the flux in the top, to allow for thin surface layers?
! ***Requires simple treatment of cloud overlap, already included below.
    do k = 1,pver_cfg
       do i = 1,ncol
          fxdust(i,k+1) = min( fxdust(i,k+1), mxsedfac * dustmr(i,k) * pdel(i,k) )
!!$        fxdust(i,k+1) = min( fxdust(i,k+1), dustmr(i,k) * pdel(i,k) + fxdust(i,k))
       end do
    end do

! Now calculate the tendencies
    do k = 1,pver_cfg
       do i = 1,ncol
! net flux into cloud changes cloud dust/ice (all flux is out of cloud)
          dusttend(i,k)  = (fxdust(i,k) - fxdust(i,k+1)) / (dtime * pdel(i,k))
       end do
    end do

! convert flux out the bottom to mass units Pa -> kg/m2/s
    sfdust(:ncol) = fxdust(:ncol,pverp_cfg) / (dtime*gravit)

    return
  end subroutine dust_sediment_tend_run

!===============================================================================
  subroutine getflx(ncol, xw, phi, vel, deltat, flux)

!.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
!....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
!....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
!.........phi1......phi2.......phi3.....phi4.......phi5.......


    implicit none

    integer ncol                      ! number of colums to process

    integer i
    integer k

    real (r8) vel(pcols_cfg,pverp_cfg)
    real (r8) flux(pcols_cfg,pverp_cfg)
    real (r8) xw(pcols_cfg,pverp_cfg)
    real (r8) psi(pcols_cfg,pverp_cfg)
    real (r8) phi(pcols_cfg,pverp_cfg-1)
    real (r8) fdot(pcols_cfg,pverp_cfg)
    real (r8) xx(pcols_cfg)
    real (r8) fxdot(pcols_cfg)
    real (r8) fxdd(pcols_cfg)

    real (r8) psistar(pcols_cfg)
    real (r8) deltat

    real (r8) xxk(pcols_cfg,pver_cfg)

    do i = 1,ncol
!        integral of phi
       psi(i,1) = 0._r8
!        fluxes at boundaries
       flux(i,1) = 0
       flux(i,pverp_cfg) = 0._r8
    end do

!     integral function
    do k = 2,pverp_cfg
       do i = 1,ncol
          psi(i,k) = phi(i,k-1)*(xw(i,k)-xw(i,k-1)) + psi(i,k-1)
       end do
    end do


!     calculate the derivatives for the interpolating polynomial
    call cfdotmc_pro (ncol, xw, psi, fdot)

!  NEW WAY
!     calculate fluxes at interior pts
    do k = 2,pver_cfg
       do i = 1,ncol
          xxk(i,k) = xw(i,k)-vel(i,k)*deltat
       end do
    end do
    do k = 2,pver_cfg
       call cfint2(ncol, xw, psi, fdot, xxk(1,k), fxdot, fxdd, psistar)
       do i = 1,ncol
          flux(i,k) = (psi(i,k)-psistar(i))
       end do
    end do


    return
  end subroutine getflx



!##############################################################################

  subroutine cfint2 (ncol, x, f, fdot, xin, fxdot, fxdd, psistar)


    implicit none

! input
    integer ncol                      ! number of colums to process

    real (r8) x(pcols_cfg, pverp_cfg)
    real (r8) f(pcols_cfg, pverp_cfg)
    real (r8) fdot(pcols_cfg, pverp_cfg)
    real (r8) xin(pcols_cfg)

! output
    real (r8) fxdot(pcols_cfg)
    real (r8) fxdd(pcols_cfg)
    real (r8) psistar(pcols_cfg)

    integer i
    integer k
    integer intz(pcols_cfg)
    real (r8) dx
    real (r8) s
    real (r8) c2
    real (r8) c3
    real (r8) xx
    real (r8) xinf
    real (r8) psi1, psi2, psi3, psim
    real (r8) cfint
    real (r8) cfnew
    real (r8) xins(pcols_cfg)

!     the minmod function
    real (r8) a, b, c
    real (r8) minmod
    real (r8) medan
    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do i = 1,ncol
       xins(i) = medan(x(i,1), xin(i), x(i,pverp_cfg))
       intz(i) = 0
    end do

! first find the interval
    do k =  1,pverp_cfg-1
       do i = 1,ncol
          if ((xins(i)-x(i,k))*(x(i,k+1)-xins(i)).ge.0._r8) then
             intz(i) = k
          endif
       end do
    end do

    do i = 1,ncol
       if (intz(i).eq.0) then
          write(core_errmsg,'(a,i0)') &
               'dust_sediment_tend_run: interpolation interval missing for column ', i
          core_errflg = 1
          return
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
       if (k+1.eq.pverp_cfg) then
          psi3 = f(i,pverp_cfg)
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

!     prototype version; eventually replace with final SPITFIRE scheme

!     calculate the derivative for the interpolating polynomial
!     multi column version


    implicit none

! input
    integer ncol                      ! number of colums to process

    real (r8) x(pcols_cfg, pverp_cfg)
    real (r8) f(pcols_cfg, pverp_cfg)
! output
    real (r8) fdot(pcols_cfg, pverp_cfg)          ! derivative at nodes

! assumed variable distribution
!     x1.......x2.......x3.......x4.......x5.......x6     1,pverp_cfg points
!     f1.......f2.......f3.......f4.......f5.......f6     1,pverp_cfg points
!     ...sh1.......sh2......sh3......sh4......sh5....     1,pver_cfg points
!     .........d2.......d3.......d4.......d5.........     2,pver_cfg points
!     .........s2.......s3.......s4.......s5.........     2,pver_cfg points
!     .............dh2......dh3......dh4.............     2,pver_cfg-1 points
!     .............eh2......eh3......eh4.............     2,pver_cfg-1 points
!     ..................e3.......e4..................     3,pver_cfg-1 points
!     .................ppl3......ppl4................     3,pver_cfg-1 points
!     .................ppr3......ppr4................     3,pver_cfg-1 points
!     .................t3........t4..................     3,pver_cfg-1 points
!     ................fdot3.....fdot4................     3,pver_cfg-1 points


! work variables


    integer i
    integer k

    real (r8) a                    ! work var
    real (r8) b                    ! work var
    real (r8) c                    ! work var
    real (r8) s(pcols_cfg,pverp_cfg)             ! first divided differences at nodes
    real (r8) sh(pcols_cfg,pverp_cfg)            ! first divided differences between nodes
    real (r8) d(pcols_cfg,pverp_cfg)             ! second divided differences at nodes
    real (r8) dh(pcols_cfg,pverp_cfg)            ! second divided differences between nodes
    real (r8) e(pcols_cfg,pverp_cfg)             ! third divided differences at nodes
    real (r8) eh(pcols_cfg,pverp_cfg)            ! third divided differences between nodes
    real (r8) pp                   ! p prime
    real (r8) ppl(pcols_cfg,pverp_cfg)           ! p prime on left
    real (r8) ppr(pcols_cfg,pverp_cfg)           ! p prime on right
    real (r8) qpl
    real (r8) qpr
    real (r8) ttt
    real (r8) t
    real (r8) tmin
    real (r8) tmax
    real (r8) delxh(pcols_cfg,pverp_cfg)


!     the minmod function
    real (r8) minmod
    real (r8) medan
    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do k = 1,pver_cfg


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
    do k = 2,pver_cfg-1
       do i = 1, ncol
          eh(i,k) = (d(i,k+1)-d(i,k))/(x(i,k+2)-x(i,k-1))
          dh(i,k) = minmod(d(i,k),d(i,k+1))
       end do
    end do

!     treat the boundaries
    do i = 1,ncol
       e(i,2) = eh(i,2)
       e(i,pver_cfg) = eh(i,pver_cfg-1)
!        outside level
       fdot(i,1) = sh(i,1) - d(i,2)*delxh(i,1)  &
            - eh(i,2)*delxh(i,1)*(x(i,1)-x(i,3))
       fdot(i,1) = minmod(fdot(i,1),3*sh(i,1))
       fdot(i,pverp_cfg) = sh(i,pver_cfg) + d(i,pver_cfg)*delxh(i,pver_cfg)  &
            + eh(i,pver_cfg-1)*delxh(i,pver_cfg)*(x(i,pverp_cfg)-x(i,pver_cfg-1))
       fdot(i,pverp_cfg) = minmod(fdot(i,pverp_cfg),3*sh(i,pver_cfg))
!        one in from boundary
       fdot(i,2) = sh(i,1) + d(i,2)*delxh(i,1) - eh(i,2)*delxh(i,1)*delxh(i,2)
       fdot(i,2) = minmod(fdot(i,2),3*s(i,2))
       fdot(i,pver_cfg) = sh(i,pver_cfg) - d(i,pver_cfg)*delxh(i,pver_cfg)   &
            - eh(i,pver_cfg-1)*delxh(i,pver_cfg)*delxh(i,pver_cfg-1)
       fdot(i,pver_cfg) = minmod(fdot(i,pver_cfg),3*s(i,pver_cfg))
    end do


    do k = 3,pver_cfg-1
       do i = 1,ncol
          e(i,k) = minmod(eh(i,k),eh(i,k-1))
       end do
    end do



    do k = 3,pver_cfg-1

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
end module ap_dust_sediment_tend_scheme
