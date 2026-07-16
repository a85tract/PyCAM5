module drydep_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid

      ! Shared Data for dry deposition calculation.

      real(r8) rair                ! Gas constant for dry air (J/K/kg)
      real(r8) gravit              ! Gravitational acceleration
!      real(r8), allocatable :: phi(:)           ! grid latitudes (radians)11

contains

!##############################################################################

! $Id$

      subroutine inidrydep( xrair, xgravit) !, xphi )

! Initialize dry deposition parameterization.

      implicit none

! Input arguments:
      real(r8), intent(in) :: xrair                ! Gas constant for dry air
      real(r8), intent(in) :: xgravit              ! Gravitational acceleration
!      real(r8), intent(in) :: xphi(:)           ! grid latitudes (radians)

! Local variables:
      integer i, j, ncid, vid, ns
!-----------------------------------------------------------------------
!      ns = size(xphi)
!      allocate(phi(ns))
      rair = xrair
      gravit = xgravit
!      do j = 1, ns
!         phi(j) = xphi(j)
!      end do

      return
      end subroutine inidrydep

!##############################################################################

      subroutine setdvel( ncol, landfrac, icefrac, ocnfrac, vgl, vgo, vgsi, vg )

! Set the deposition velocity depending on whether we are over
! land, ocean, and snow/ice


      implicit none

! Input arguments:

      integer, intent(in) :: ncol
      real (r8), intent(in) :: landfrac(pcols)       ! land fraction
      real (r8), intent(in) :: icefrac(pcols)       ! ice fraction
      real (r8), intent(in) :: ocnfrac(pcols)       ! ocean fraction

      real(r8), intent(in) :: vgl                  ! dry deposition velocity in m/s (land)
      real(r8), intent(in) :: vgo                  ! dry deposition velocity in m/s (ocean)
      real(r8), intent(in) :: vgsi                 ! dry deposition velocity in m/s (snow/ice)

! Output arguments:
      real(r8), intent(out) ::  vg(pcols) ! dry deposition velocity in m/s

! Local variables:

      integer i
      real(r8) a


      do i = 1, ncol
         vg(i) = landfrac(i)*vgl + ocnfrac(i)*vgo + icefrac(i)*vgsi
!         if (ioro(i).eq.0) then
!            vg(i) = vgo
!         else if (ioro(i).eq.1) then
!            vg(i) = vgl
!         else
!            vg(i) = vgsi
!         endif
      end do

      return
      end subroutine setdvel

!##############################################################################

      subroutine ddflux( ncol, vg, q, p, tv, flux )

! Compute surface flux due to dry deposition processes.


      implicit none

! Input arguments:
      integer , intent(in) :: ncol
      real(r8), intent(in) ::    vg(pcols)  ! dry deposition velocity in m/s
      real(r8), intent(in) ::    q(pcols)   ! tracer conc. in surface layer (kg tracer/kg moist air)
      real(r8), intent(in) ::    p(pcols)   ! midpoint pressure in surface layer (Pa)
      real(r8), intent(in) ::    tv(pcols)  ! midpoint virtual temperature in surface layer (K)

! Output arguments:

      real(r8), intent(out) ::    flux(pcols) ! flux due to dry deposition in kg/m^s/sec

! Local variables:

      integer i

      do i = 1, ncol
         flux(i) = -vg(i) * q(i) * p(i) /(tv(i) * rair)
      end do

      return
      end subroutine ddflux

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine d3ddflux
!
! !INTERFACE:
!
   subroutine  d3ddflux ( ncol, vlc_dry, q,pmid,pdel, tv, dep_dry,dep_dry_tend,dt)
! Description:
!Do 3d- settling deposition calculations following Zender's dust codes, Dec 02.
!
! Author: Natalie Mahowald
!
      use ap_d3ddflux_scheme, only: d3ddflux_run
      use perf_mod, only: t_startf, t_stopf
      implicit none

! Input arguments:
      integer , intent(in) :: ncol
      real(r8), intent(in) ::    vlc_dry(pcols,pver)  ! dry deposition velocity in m/s
      real(r8), intent(in) ::    q(pcols,pver)   ! tracer conc. in surface layer (kg tracer/kg moist air)
      real(r8), intent(in) ::    pmid(pcols,pver)   ! midpoint pressure in surface layer (Pa)
      real(r8), intent(in) ::    pdel(pcols,pver)   ! delta pressure across level (Pa)
      real(r8), intent(in) ::    tv(pcols,pver)  ! midpoint virtual temperature in surface layer (K)
    real(r8),            intent(in)  :: dt             ! time step

! Output arguments:

      real(r8), intent(out) ::    dep_dry(pcols) ! flux due to dry deposition in kg /m^s/sec
      real(r8), intent(out) ::    dep_dry_tend(pcols,pver) ! flux due to dry deposition in kg /m^s/sec

      call t_startf('ap_d3ddflux_run')
      call d3ddflux_run(pcols, pver, ncol, vlc_dry, q, pmid, pdel, tv, &
           dep_dry, dep_dry_tend, dt, rair, gravit)
      call t_stopf('ap_d3ddflux_run')
      return
      end subroutine d3ddflux



!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subroutine Calcram
!
! !INTERFACE:
!

      subroutine  calcram(ncol,landfrac,icefrac,ocnfrac,obklen,&
           ustar,ram1in,ram1,t,pmid,&
           pdel,fvin,fv)
        use ap_calcram_scheme, only: calcram_run
        use perf_mod, only: t_startf, t_stopf
        !
        ! !DESCRIPTION: 
        !  
        ! Calc aerodynamic resistance over oceans and sea ice (comes in from land model)
        ! from Seinfeld and Pandis, p.963.
        !  
        ! Author: Natalie Mahowald
        !
        implicit none
        integer, intent(in) :: ncol
        real(r8),intent(in) :: ram1in(pcols)         !aerodynamical resistance (s/m)
        real(r8),intent(in) :: fvin(pcols)                 ! sfc frc vel from land
        real(r8),intent(out) :: ram1(pcols)         !aerodynamical resistance (s/m)
        real(r8),intent(out) :: fv(pcols)                 ! sfc frc vel from land
        real(r8), intent(in) :: obklen(pcols)                 ! obklen
        real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
        real(r8), intent(in) :: landfrac(pcols)               ! land fraction
        real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
        real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
        real(r8), intent(in) :: t(pcols)       !atm temperature (K)
        real(r8), intent(in) :: pmid(pcols)    !atm pressure (Pa)
        real(r8), intent(in) :: pdel(pcols)    !atm pressure (Pa)
        call t_startf('ap_calcram_run')
        call calcram_run(pcols, rair, gravit, ncol, landfrac, icefrac, ocnfrac, obklen, &
             ustar, ram1in, ram1, t, pmid, pdel, fvin, fv)
        call t_stopf('ap_calcram_run')

        return
      end subroutine calcram


!##############################################################################
end module drydep_mod
