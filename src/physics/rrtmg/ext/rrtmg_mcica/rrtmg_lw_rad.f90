!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.6 $
!     created:   $Date: 2008/04/24 16:17:27 $
!

       module rrtmg_lw_rad

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                              RRTMG_LW                                    *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                   a rapid radiative transfer model                       *
! *                       for the longwave region                            * 
! *             for application to general circulation models                *
! *                                                                          *
! *                                                                          *
! *            Atmospheric and Environmental Research, Inc.                  *
! *                        131 Hartwell Avenue                               *
! *                        Lexington, MA 02421                               *
! *                                                                          *
! *                                                                          *
! *                           Eli J. Mlawer                                  *
! *                        Jennifer S. Delamere                              *
! *                         Michael J. Iacono                                *
! *                         Shepard A. Clough                                *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                       email:  miacono@aer.com                            *
! *                       email:  emlawer@aer.com                            *
! *                       email:  jdelamer@aer.com                           *
! *                                                                          *
! *        The authors wish to acknowledge the contributions of the          *
! *        following people:  Steven J. Taubman, Karen Cady-Pereira,         *
! *        Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert Bergstrom.  *
! *                                                                          *
! ****************************************************************************

! -------- Modules --------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use ppgrid,       only: pcols, begchunk, endchunk

!      use parkind, only : jpim, jprb 
      use rrlw_vsn
      use mcica_subcol_gen_lw, only: mcica_subcol_lw
      use rrtmg_lw_cldprmc, only: cldprmc
! Move call to rrtmg_lw_ini and following use association to 
! GCM initialization area
!      use rrtmg_lw_init, only: rrtmg_lw_ini
      use rrtmg_lw_rtrnmc, only: rtrnmc
      use rrtmg_lw_setcoef, only: setcoef
      use rrtmg_lw_taumol, only: taumol

      implicit none

      logical :: use_native_rrtmg_lw_rad_pack_impl = .false.
      logical :: rrtmg_lw_rad_pack_impl_selected = .false.
      logical :: rrtmg_lw_rad_pack_entered_logged = .false.
      logical :: use_native_rrtmg_lw_inatm_impl = .false.
      logical :: rrtmg_lw_inatm_impl_selected = .false.
      logical :: rrtmg_lw_inatm_entered_logged = .false.

! public interfaces/functions/subroutines
      public :: rrtmg_lw, inatm

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_lw &
            (lchnk   ,ncol    ,nlay    ,icld    ,                   &
             play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr  , &
             o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  ,&
             cfc11vmr,cfc12vmr, &
             cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
             cldfmcl ,taucmcl ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, uflxs, dflxs )

! -------- Description --------

! This program is the driver subroutine for RRTMG_LW, the AER LW radiation 
! model for application to GCMs, that has been adapted from RRTM_LW for
! improved efficiency.
!
! NOTE: The call to RRTMG_LW_INI should be moved to the GCM initialization
!  area, since this has to be called only once. 
!
! This routine:
!    a) calls INATM to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa. 
!    b) calls CLDPRMC to set cloud optical depth for McICA based 
!       on input cloud properties 
!    c) calls SETCOEF to calculate various quantities needed for 
!       the radiative transfer algorithm
!    d) calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands
!    e) calls RTRNMC (for both clear and cloudy profiles) to perform the
!       radiative transfer calculation using McICA, the Monte-Carlo 
!       Independent Column Approximation, to represent sub-grid scale 
!       cloud variability
!    f) passes the necessary fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_lw.nomcica.f90 (to not use
!     McICA) or rrtmg_lw.f90 (to use McICA) to interface with a GCM. 
!
!    1) Standard, single forward model calculation (imca = 0)
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!
! This call to RRTMG_LW must be preceeded by a call to the module
!     mcica_subcol_gen_lw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngpt) dimension.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflglw, iceflglw, and liqflglw; see text file rrtmg_lw_instructions
!     and subroutine rrtmg_lw_cldprop.f90 for further details):
!
!    1) Input cloud fraction and cloud optical depth directly (inflglw = 0)
!    2) Input cloud fraction and cloud physical properties (inflglw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflglw and liqflglw
!
! One method of aerosol property input is possible:
!     Aerosol properties can be input in only one way (controlled by input 
!     flag iaer, see text file rrtmg_lw_instructions for further details):
!
!    1) Input aerosol optical depth directly by layer and spectral band (iaer=10);
!       band average optical depth at the mid-point of each spectral band.
!       RRTMG_LW currently treats only aerosol absorption;
!       scattering capability is not presently available. 
!
!
! ------- Modifications -------
!
! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced 
! set of g-points for application to GCMs.  
!
!-- Original version (derived from RRTM_LW), reduction of g-points, other
!   revisions for use with GCMs.  
!     1999: M. J. Iacono, AER, Inc.
!-- Adapted for use with NCAR/CAM.
!     May 2004: M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability. 
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Conversion to F90 formatting for consistency with rrtmg_sw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays.
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to add longwave aerosol absorption.
!     Apr 2008: M. J. Iacono, AER, Inc.

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol
      use iso_c_binding, only: c_int64_t, c_loc, c_ptr
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi
      use rrlw_wvn, only: ng, ngb, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: lchnk                      ! chunk identifier
      integer, intent(in) :: ncol                       ! Number of horizontal columns
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(inout) :: icld                    ! Cloud overlap method
                                                        !    0: Clear only
                                                        !    1: Random
                                                        !    2: Maximum/random
                                                        !    3: Maximum
      real(kind=r8), intent(in) :: play(:,:)            ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(:,:)            ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(:,:)            ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(:)              ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc11vmr(:,:)        ! CFC11 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc12vmr(:,:)        ! CFC12 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc22vmr(:,:)        ! CFC22 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ccl4vmr(:,:)         ! CCL4 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: emis(:,:)            ! Surface emissivity
                                                        !    Dimensions: (ncol,nbndlw)

      integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfmcl(:,:,:)       ! Cloud fraction
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: ciwpmcl(:,:,:)       ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: clwpmcl(:,:,:)       ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: reicmcl(:,:)         ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: relqmcl(:,:)         ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: taucmcl(:,:,:)       ! Cloud optical depth
                                                        !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=r8), intent(in) :: ssacmcl(:,:,:)      ! Cloud single scattering albedo
                                                        !    Dimensions: (ngptlw,ncol,nlay)
                                                        !   for future expansion
                                                        !   lw scattering not yet available
!      real(kind=r8), intent(in) :: asmcmcl(:,:,:)      ! Cloud asymmetry parameter
                                                        !    Dimensions: (ngptlw,ncol,nlay)
                                                        !   for future expansion
                                                        !   lw scattering not yet available
      real(kind=r8), intent(in) :: tauaer(:,:,:)        ! aerosol optical depth
                                                        !   at mid-point of LW spectral bands
                                                        !    Dimensions: (ncol,nlay,nbndlw)
!      real(kind=r8), intent(in) :: ssaaer(:,:,:)       ! aerosol single scattering albedo
                                                        !    Dimensions: (ncol,nlay,nbndlw)
                                                        !   for future expansion 
                                                        !   (lw aerosols/scattering not yet available)
!      real(kind=r8), intent(in) :: asmaer(:,:,:)       ! aerosol asymmetry parameter
                                                        !    Dimensions: (ncol,nlay,nbndlw)
                                                        !   for future expansion 
                                                        !   (lw aerosols/scattering not yet available)

! ----- Output -----

      real(kind=r8), target, intent(out) :: uflx(:,:)   ! Total sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), target, intent(out) :: dflx(:,:)   ! Total sky longwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), target, intent(out) :: hr(:,:)     ! Total sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), target, intent(out) :: uflxc(:,:)  ! Clear sky longwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), target, intent(out) :: dflxc(:,:)  ! Clear sky longwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), target, intent(out) :: hrc(:,:)    ! Clear sky longwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), target, intent(out) :: uflxs(:,:,:)! Total sky longwave upward flux spectral (W/m2)
                                                        !    Dimensions: (nbndlw,ncol,nlay+1)
      real(kind=r8), target, intent(out) :: dflxs(:,:,:)! Total sky longwave downward flux spectral (W/m2)
                                                        !    Dimensions: (nbndlw,ncol,nlay+1)

! ----- Local -----

! Control
      integer :: istart                         ! beginning band of calculation
      integer :: iend                           ! ending band of calculation
      integer :: iout                           ! output option flag (inactive)
      integer :: iaer                           ! aerosol option flag
      integer :: iplon                          ! column loop index
      integer :: imca                           ! flag for mcica [0=off, 1=on]
      integer :: ims                            ! value for changing mcica permute seed
      integer :: k                              ! layer loop index
      integer :: ig                             ! g-point loop index

! Atmosphere
      real(kind=r8) :: pavel(nlay)              ! layer pressures (mb) 
      real(kind=r8) :: tavel(nlay)              ! layer temperatures (K)
      real(kind=r8) :: pz(0:nlay)               ! level (interface) pressures (hPa, mb)
      real(kind=r8) :: tz(0:nlay)               ! level (interface) temperatures (K)
      real(kind=r8) :: tbound                   ! surface temperature (K)
      real(kind=r8) :: coldry(nlay)             ! dry air column density (mol/cm2)
      real(kind=r8) :: wbrodl(nlay)             ! broadening gas column density (mol/cm2)
      real(kind=r8) :: wkl(mxmol,nlay)          ! molecular amounts (mol/cm-2)
      real(kind=r8) :: wx(maxxsec,nlay)         ! cross-section amounts (mol/cm-2)
      real(kind=r8) :: pwvcm                    ! precipitable water vapor (cm)
      real(kind=r8) :: semiss(nbndlw)           ! lw surface emissivity
      real(kind=r8) :: fracs(nlay,ngptlw)       ! 
      real(kind=r8), target :: taug(nlay,ngptlw)! gaseous optical depths
      real(kind=r8), target :: taut(nlay,ngptlw)! gaseous + aerosol optical depths

      real(kind=r8), target :: taua(nlay,nbndlw)! aerosol optical depth
!      real(kind=r8) :: ssaa(nlay,nbndlw)        ! aerosol single scattering albedo
                                                 !   for future expansion 
                                                 !   (lw aerosols/scattering not yet available)
!      real(kind=r8) :: asma(nlay+1,nbndlw)      ! aerosol asymmetry parameter
                                                 !   for future expansion 
                                                 !   (lw aerosols/scattering not yet available)

! Atmosphere - setcoef
      integer :: laytrop                          ! tropopause layer index
      integer :: jp(nlay)                         ! lookup table index 
      integer :: jt(nlay)                         ! lookup table index 
      integer :: jt1(nlay)                        ! lookup table index 
      real(kind=r8) :: planklay(nlay,nbndlw)      ! 
      real(kind=r8) :: planklev(0:nlay,nbndlw)    ! 
      real(kind=r8) :: plankbnd(nbndlw)           ! 

      real(kind=r8) :: colh2o(nlay)               ! column amount (h2o)
      real(kind=r8) :: colco2(nlay)               ! column amount (co2)
      real(kind=r8) :: colo3(nlay)                ! column amount (o3)
      real(kind=r8) :: coln2o(nlay)               ! column amount (n2o)
      real(kind=r8) :: colco(nlay)                ! column amount (co)
      real(kind=r8) :: colch4(nlay)               ! column amount (ch4)
      real(kind=r8) :: colo2(nlay)                ! column amount (o2)
      real(kind=r8) :: colbrd(nlay)               ! column amount (broadening gases)

      integer :: indself(nlay)
      integer :: indfor(nlay)
      real(kind=r8) :: selffac(nlay)
      real(kind=r8) :: selffrac(nlay)
      real(kind=r8) :: forfac(nlay)
      real(kind=r8) :: forfrac(nlay)

      integer :: indminor(nlay)
      real(kind=r8) :: minorfrac(nlay)
      real(kind=r8) :: scaleminor(nlay)
      real(kind=r8) :: scaleminorn2(nlay)

      real(kind=r8) :: &                          !
                         fac00(nlay), fac01(nlay), &
                         fac10(nlay), fac11(nlay) 
      real(kind=r8) :: &                          !
                         rat_h2oco2(nlay),rat_h2oco2_1(nlay), &
                         rat_h2oo3(nlay),rat_h2oo3_1(nlay), &
                         rat_h2on2o(nlay),rat_h2on2o_1(nlay), &
                         rat_h2och4(nlay),rat_h2och4_1(nlay), &
                         rat_n2oco2(nlay),rat_n2oco2_1(nlay), &
                         rat_o3co2(nlay),rat_o3co2_1(nlay)

! Atmosphere/clouds - cldprop
      integer :: ncbands                          ! number of cloud spectral bands
      integer :: inflag                           ! flag for cloud property method
      integer :: iceflag                          ! flag for ice cloud properties
      integer :: liqflag                          ! flag for liquid cloud properties

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=r8) :: cldfmc(ngptlw,nlay)      ! cloud fraction [mcica]
      real(kind=r8) :: ciwpmc(ngptlw,nlay)      ! cloud ice water path [mcica]
      real(kind=r8) :: clwpmc(ngptlw,nlay)      ! cloud liquid water path [mcica]
      real(kind=r8) :: relqmc(nlay)             ! liquid particle size (microns)
      real(kind=r8) :: reicmc(nlay)             ! ice particle effective radius (microns)
      real(kind=r8) :: dgesmc(nlay)             ! ice particle generalized effective size (microns)
      real(kind=r8) :: taucmc(ngptlw,nlay)      ! cloud optical depth [mcica]
!      real(kind=r8) :: ssacmc(ngptlw,nlay)     ! cloud single scattering albedo [mcica]
                                                !   for future expansion 
                                                !   (lw scattering not yet available)
!      real(kind=r8) :: asmcmc(ngptlw,nlay)     ! cloud asymmetry parameter [mcica]
                                                !   for future expansion 
                                                !   (lw scattering not yet available)

! Output
      real(kind=r8), target :: totuflux(0:nlay) ! upward longwave flux (w/m2)
      real(kind=r8), target :: totdflux(0:nlay) ! downward longwave flux (w/m2)
      real(kind=r8), target :: totufluxs(nbndlw,0:nlay)! upward longwave flux spectral (w/m2)
      real(kind=r8), target :: totdfluxs(nbndlw,0:nlay)! downward longwave flux spectral (w/m2)
      real(kind=r8) :: fnet(0:nlay)             ! net longwave flux (w/m2)
      real(kind=r8), target :: htr(0:nlay)      ! longwave heating rate (k/day)
      real(kind=r8), target :: totuclfl(0:nlay) ! clear sky upward longwave flux (w/m2)
      real(kind=r8), target :: totdclfl(0:nlay) ! clear sky downward longwave flux (w/m2)
      real(kind=r8) :: fnetc(0:nlay)            ! clear sky net longwave flux (w/m2)
      real(kind=r8), target :: htrc(0:nlay)     ! clear sky longwave heating rate (k/day)
      integer(c_int64_t), target :: ngb64(ngptlw)

      interface
         subroutine rrtmg_lw_rad_taut_codon(nlay_c, ngptlw_c, iaer_c, taug_p, taua_p, ngb_p, &
              taut_p) bind(c, name="rrtmg_lw_rad_taut_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlay_c, ngptlw_c, iaer_c
            type(c_ptr), value :: taug_p, taua_p, ngb_p, taut_p
         end subroutine rrtmg_lw_rad_taut_codon

         subroutine rrtmg_lw_rad_store_flux_codon(nlay_c, nbndlw_c, ncol_c, iplon_c, &
              totuflux_p, totdflux_p, totuclfl_p, totdclfl_p, totufluxs_p, totdfluxs_p, &
              htr_p, htrc_p, uflx_p, dflx_p, uflxc_p, dflxc_p, uflxs_p, dflxs_p, hr_p, &
              hrc_p) bind(c, name="rrtmg_lw_rad_store_flux_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlay_c, nbndlw_c, ncol_c, iplon_c
            type(c_ptr), value :: totuflux_p, totdflux_p, totuclfl_p, totdclfl_p
            type(c_ptr), value :: totufluxs_p, totdfluxs_p, htr_p, htrc_p
            type(c_ptr), value :: uflx_p, dflx_p, uflxc_p, dflxc_p, uflxs_p, dflxs_p
            type(c_ptr), value :: hr_p, hrc_p
         end subroutine rrtmg_lw_rad_store_flux_codon
      end interface

! Initializations

      oneminus = 1._r8 - 1.e-6_r8
      pi = 2._r8 * asin(1._r8)
      fluxfac = pi * 2.e4_r8                    ! orig:   fluxfac = pi * 2.d4  
      istart = 1
      iend = 16
      iout = 0
      ims = 1

! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability

! *** This version uses McICA (imca = 1) ***

! Set icld to select of clear or cloud calculation and cloud overlap method  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap
! icld = 2, with clouds using maximum/random cloud overlap
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (icld.lt.0.or.icld.gt.3) icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 10, input total aerosol optical depth (tauaer) directly 
      iaer = 10

      call rrtmg_lw_rad_pack_select_impl()
      do ig = 1, ngptlw
         ngb64(ig) = int(ngb(ig), c_int64_t)
      enddo

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 256 to 140 for input absorption coefficient 
! data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_lw_ini

!  This is the main longitude/column loop within RRTMG.
      do iplon = 1, ncol

!  Prepare atmospheric profile from GCM for use in RRTMG, and define
!  other input parameters.  

         call inatm (iplon, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, h2ovmr, &
              o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, cfc11vmr, cfc12vmr, &
              cfc22vmr, ccl4vmr, emis, inflglw, iceflglw, liqflglw, &
              cldfmcl, taucmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, tauaer, &
              pavel, pz, tavel, tz, tbound, semiss, coldry, &
              wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, &
              cldfmc, taucmc, ciwpmc, clwpmc, reicmc, dgesmc, relqmc, taua)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed into cldprop.  Cloud fraction and cloud
!  optical depth are transferred to rrtmg_lw arrays in cldprop.  

         call cldprmc(nlay, inflag, iceflag, liqflag, cldfmc, ciwpmc, &
                      clwpmc, reicmc, dgesmc, relqmc, ncbands, taucmc)

! Calculate information needed by the radiative transfer routine
! that is specific to this atmosphere, especially some of the 
! coefficients and indices needed to compute the optical depths
! by interpolating data from stored reference atmospheres. 

         call setcoef(nlay, istart, pavel, tavel, tz, tbound, semiss, &
                      coldry, wkl, wbrodl, &
                      laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                      colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                      colbrd, fac00, fac01, fac10, fac11, &
                      rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                      rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                      rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                      selffac, selffrac, indself, forfac, forfrac, indfor, &
                      minorfrac, scaleminor, scaleminorn2, indminor)

!  Calculate the gaseous optical depths and Planck fractions for 
!  each longwave spectral band.

         call taumol(nlay, pavel, wx, coldry, &
                     laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                     colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                     colbrd, fac00, fac01, fac10, fac11, &
                     rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                     rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                     rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     minorfrac, scaleminor, scaleminorn2, indminor, &
                     fracs, taug)



! Combine gaseous and aerosol optical depths, if aerosol active
         if (use_native_rrtmg_lw_rad_pack_impl) then
            if (iaer .eq. 0) then
               do k = 1, nlay
                  do ig = 1, ngptlw
                     taut(k,ig) = taug(k,ig)
                  enddo
               enddo
            elseif (iaer .eq. 10) then
               do k = 1, nlay
                  do ig = 1, ngptlw
                     taut(k,ig) = taug(k,ig) + taua(k,ngb(ig))
                  enddo
               enddo
            endif
         else
            call rrtmg_lw_rad_pack_log_entered()
            call rrtmg_lw_rad_taut_codon( &
                 int(nlay, c_int64_t), int(ngptlw, c_int64_t), int(iaer, c_int64_t), &
                 c_loc(taug(1,1)), c_loc(taua(1,1)), c_loc(ngb64(1)), c_loc(taut(1,1)) &
            )
         endif

! Call the radiative transfer routine.
! Either routine can be called to do clear sky calculation.  If clouds
! are present, then select routine based on cloud overlap assumption
! to be used.  Clear sky calculation is done simultaneously.
! For McICA, RTRNMC is called for clear and cloudy calculations.

         call rtrnmc(nlay, istart, iend, iout, pz, semiss, ncbands, &
                     cldfmc, taucmc, planklay, planklev, plankbnd, &
                     pwvcm, fracs, taut, &
                     totuflux, totdflux, fnet, htr, &
                     totuclfl, totdclfl, fnetc, htrc, totufluxs, totdfluxs )

!  Transfer up and down fluxes and heating rate to output arrays.
!  Vertical indexing goes from bottom to top

         if (use_native_rrtmg_lw_rad_pack_impl) then
            do k = 0, nlay
               uflx(iplon,k+1) = totuflux(k)
               dflx(iplon,k+1) = totdflux(k)
               uflxc(iplon,k+1) = totuclfl(k)
               dflxc(iplon,k+1) = totdclfl(k)
               uflxs(:,iplon,k+1) = totufluxs(:,k)
               dflxs(:,iplon,k+1) = totdfluxs(:,k)
            enddo
            do k = 0, nlay-1
               hr(iplon,k+1) = htr(k)
               hrc(iplon,k+1) = htrc(k)
            enddo
         else
            call rrtmg_lw_rad_store_flux_codon( &
                 int(nlay, c_int64_t), int(nbndlw, c_int64_t), int(size(uflx,1), c_int64_t), &
                 int(iplon, c_int64_t), c_loc(totuflux(0)), c_loc(totdflux(0)), &
                 c_loc(totuclfl(0)), c_loc(totdclfl(0)), c_loc(totufluxs(1,0)), &
                 c_loc(totdfluxs(1,0)), c_loc(htr(0)), c_loc(htrc(0)), c_loc(uflx(1,1)), &
                 c_loc(dflx(1,1)), c_loc(uflxc(1,1)), c_loc(dflxc(1,1)), &
                 c_loc(uflxs(1,1,1)), c_loc(dflxs(1,1,1)), c_loc(hr(1,1)), &
                 c_loc(hrc(1,1)) &
            )
         endif

      enddo

      end subroutine rrtmg_lw

! --------------------------------------------------------------------------
      subroutine rrtmg_lw_rad_pack_select_impl()

      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (rrtmg_lw_rad_pack_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('RRTMG_LW_RAD_PACK_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_rrtmg_lw_rad_pack_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_rrtmg_lw_rad_pack_impl = .false.
      end if

      rrtmg_lw_rad_pack_impl_selected = .true.

      if (masterproc) then
         if (use_native_rrtmg_lw_rad_pack_impl) then
            write(iulog,*) 'rrtmg_lw_rad_pack implementation = native'
         else
            write(iulog,*) 'rrtmg_lw_rad_pack implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine rrtmg_lw_rad_pack_select_impl

! --------------------------------------------------------------------------
      subroutine rrtmg_lw_rad_pack_log_entered()

      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

      if (rrtmg_lw_rad_pack_entered_logged) return
      rrtmg_lw_rad_pack_entered_logged = .true.

      if (masterproc) then
         write(iulog,*) 'rrtmg_lw_rad_pack entered (mcica lw taut/flux transfer = codon)'
         call flush(iulog)
      end if

      end subroutine rrtmg_lw_rad_pack_log_entered

! --------------------------------------------------------------------------
      subroutine rrtmg_lw_inatm_select_impl()

      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (rrtmg_lw_inatm_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('RRTMG_LW_INATM_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_rrtmg_lw_inatm_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_rrtmg_lw_inatm_impl = .false.
      end if

      rrtmg_lw_inatm_impl_selected = .true.

      if (masterproc) then
         if (use_native_rrtmg_lw_inatm_impl) then
            write(iulog,*) 'rrtmg_lw_inatm implementation = native'
         else
            write(iulog,*) 'rrtmg_lw_inatm implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine rrtmg_lw_inatm_select_impl

! --------------------------------------------------------------------------
      subroutine rrtmg_lw_inatm_log_entered()

      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

      if (rrtmg_lw_inatm_entered_logged) return
      rrtmg_lw_inatm_entered_logged = .true.

      if (masterproc) then
         write(iulog,*) 'rrtmg_lw_inatm entered (mcica lw atmosphere packing = codon)'
         call flush(iulog)
      end if

      end subroutine rrtmg_lw_inatm_log_entered

!***************************************************************************
      subroutine inatm (iplon, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, h2ovmr, &
              o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, cfc11vmr, cfc12vmr, &
              cfc22vmr, ccl4vmr, emis, inflglw, iceflglw, liqflglw, &
              cldfmcl, taucmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, tauaer, &
              pavel, pz, tavel, tz, tbound, semiss, coldry, &
              wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, &
              cldfmc, taucmc, ciwpmc, clwpmc, reicmc, dgesmc, relqmc, taua)
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_LW.
!  Set other RRTMG_LW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, nmol, maxxsec, mxmol
      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi, grav, avogad
      use rrlw_wvn, only: ng, nspa, nspb, wavenum1, wavenum2, delwave, ixindx

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: iplon                      ! column loop index
      integer, intent(in) :: nlay                       ! Number of model layers
      integer, intent(in) :: icld                       ! clear/cloud and cloud overlap flag
      integer, intent(in) :: iaer                       ! aerosol option flag

      real(kind=r8), intent(in) :: play(:,:)            ! Layer pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(:,:)            ! Layer temperatures (K)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(:,:)            ! Interface temperatures (K)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(:)              ! Surface temperature (K)
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc11vmr(:,:)        ! CFC11 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc12vmr(:,:)        ! CFC12 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: cfc22vmr(:,:)        ! CFC22 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ccl4vmr(:,:)         ! CCL4 volume mixing ratio
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: emis(:,:)            ! Surface emissivity
                                                        !    Dimensions: (ncol,nbndlw)

      integer, intent(in) :: inflglw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflglw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflglw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfmcl(:,:,:)       ! Cloud fraction
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: ciwpmcl(:,:,:)       ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: clwpmcl(:,:,:)       ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: reicmcl(:,:)         ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: relqmcl(:,:)         ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: taucmcl(:,:,:)       ! Cloud optical depth
                                                        !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=r8), intent(in) :: tauaer(:,:,:)        ! Aerosol optical depth
                                                        !    Dimensions: (ncol,nlay,nbndlw)

! ----- Output -----
! Atmosphere
      real(kind=r8), intent(out) :: pavel(:)            ! layer pressures (mb) 
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: tavel(:)            ! layer temperatures (K)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: pz(0:)              ! level (interface) pressures (hPa, mb)
                                                        !    Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tz(0:)              ! level (interface) temperatures (K)
                                                        !    Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tbound              ! surface temperature (K)
      real(kind=r8), intent(out) :: coldry(:)           ! dry air column density (mol/cm2)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: wbrodl(:)           ! broadening gas column density (mol/cm2)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: wkl(:,:)            ! molecular amounts (mol/cm-2)
                                                        !    Dimensions: (mxmol,nlay)
      real(kind=r8), intent(out) :: wx(:,:)             ! cross-section amounts (mol/cm-2)
                                                        !    Dimensions: (maxxsec,nlay)
      real(kind=r8), intent(out) :: pwvcm               ! precipitable water vapor (cm)
      real(kind=r8), intent(out) :: semiss(:)           ! lw surface emissivity
                                                        !    Dimensions: (nbndlw)

! Atmosphere/clouds - cldprop
      integer, intent(out) :: inflag                    ! flag for cloud property method
      integer, intent(out) :: iceflag                   ! flag for ice cloud properties
      integer, intent(out) :: liqflag                   ! flag for liquid cloud properties

      real(kind=r8), intent(out) :: cldfmc(:,:)         ! cloud fraction [mcica]
                                                        !    Dimensions: (ngptlw,nlay)
      real(kind=r8), intent(out) :: ciwpmc(:,:)         ! cloud ice water path [mcica]
                                                        !    Dimensions: (ngptlw,nlay)
      real(kind=r8), intent(out) :: clwpmc(:,:)         ! cloud liquid water path [mcica]
                                                        !    Dimensions: (ngptlw,nlay)
      real(kind=r8), intent(out) :: relqmc(:)           ! liquid particle effective radius (microns)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: reicmc(:)           ! ice particle effective radius (microns)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: dgesmc(:)           ! ice particle generalized effective size (microns)
                                                        !    Dimensions: (nlay)
      real(kind=r8), intent(out) :: taucmc(:,:)         ! cloud optical depth [mcica]
                                                        !    Dimensions: (ngptlw,nlay)
      real(kind=r8), intent(out) :: taua(:,:)           ! Aerosol optical depth
                                                        ! Dimensions: (nlay,nbndlw)


! ----- Local -----
      real(kind=r8), parameter :: amd = 28.9660_r8      ! Effective molecular weight of dry air (g/mol)
      real(kind=r8), parameter :: amw = 18.0160_r8      ! Molecular weight of water vapor (g/mol)
!      real(kind=r8), parameter :: amc = 44.0098_r8      ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=r8), parameter :: amo = 47.9998_r8      ! Molecular weight of ozone (g/mol)
!      real(kind=r8), parameter :: amo2 = 31.9999_r8     ! Molecular weight of oxygen (g/mol)
!      real(kind=r8), parameter :: amch4 = 16.0430_r8    ! Molecular weight of methane (g/mol)
!      real(kind=r8), parameter :: amn2o = 44.0128_r8    ! Molecular weight of nitrous oxide (g/mol)
!      real(kind=r8), parameter :: amc11 = 137.3684_r8   ! Molecular weight of CFC11 (g/mol) - CCL3F
!      real(kind=r8), parameter :: amc12 = 120.9138_r8   ! Molecular weight of CFC12 (g/mol) - CCL2F2
!      real(kind=r8), parameter :: amc22 = 86.4688_r8    ! Molecular weight of CFC22 (g/mol) - CHCLF2
!      real(kind=r8), parameter :: amcl4 = 153.823_r8    ! Molecular weight of CCL4 (g/mol) - CCL4

! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
      real(kind=r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
      real(kind=r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
      real(kind=r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
      real(kind=r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
      real(kind=r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
      real(kind=r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

      real(kind=r8), parameter :: sbc = 5.67e-08_r8     ! Stefan-Boltzmann constant (W/m2K4)

      integer :: isp, l, ix, n, imol, ib, ig            ! Loop indices
      real(kind=r8) :: amm, amttl, wvttl, wvsh, summol  
      integer(c_int64_t), target :: inflag64, iceflag64, liqflag64
      integer(c_int64_t), target :: ixindx64(maxxsec)

      interface
         subroutine rrtmg_lw_inatm_codon(iplon_c, nlay_c, ldcol_c, icld_c, iaer_c, &
              nbndlw_c, ngptlw_c, nmol_c, maxxsec_c, mxmol_c, grav_c, avogad_c, &
              play_p, plev_p, tlay_p, tlev_p, tsfc_p, h2ovmr_p, o3vmr_p, co2vmr_p, &
              ch4vmr_p, o2vmr_p, n2ovmr_p, cfc11vmr_p, cfc12vmr_p, cfc22vmr_p, &
              ccl4vmr_p, emis_p, inflglw_c, iceflglw_c, liqflglw_c, cldfmcl_p, &
              taucmcl_p, ciwpmcl_p, clwpmcl_p, reicmcl_p, relqmcl_p, tauaer_p, &
              pavel_p, pz_p, tavel_p, tz_p, tbound_p, semiss_p, coldry_p, wbrodl_p, &
              wkl_p, wx_p, pwvcm_p, inflag_p, iceflag_p, liqflag_p, cldfmc_p, &
              taucmc_p, ciwpmc_p, clwpmc_p, reicmc_p, dgesmc_p, relqmc_p, taua_p, &
              ixindx_p) bind(c, name="rrtmg_lw_inatm_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: iplon_c, nlay_c, ldcol_c, icld_c, iaer_c
            integer(c_int64_t), value :: nbndlw_c, ngptlw_c, nmol_c, maxxsec_c, mxmol_c
            integer(c_int64_t), value :: inflglw_c, iceflglw_c, liqflglw_c
            real(c_double), value :: grav_c, avogad_c
            type(c_ptr), value :: play_p, plev_p, tlay_p, tlev_p, tsfc_p, h2ovmr_p
            type(c_ptr), value :: o3vmr_p, co2vmr_p, ch4vmr_p, o2vmr_p, n2ovmr_p
            type(c_ptr), value :: cfc11vmr_p, cfc12vmr_p, cfc22vmr_p, ccl4vmr_p, emis_p
            type(c_ptr), value :: cldfmcl_p, taucmcl_p, ciwpmcl_p, clwpmcl_p
            type(c_ptr), value :: reicmcl_p, relqmcl_p, tauaer_p, pavel_p, pz_p, tavel_p
            type(c_ptr), value :: tz_p, tbound_p, semiss_p, coldry_p, wbrodl_p, wkl_p
            type(c_ptr), value :: wx_p, pwvcm_p, inflag_p, iceflag_p, liqflag_p, cldfmc_p
            type(c_ptr), value :: taucmc_p, ciwpmc_p, clwpmc_p, reicmc_p, dgesmc_p
            type(c_ptr), value :: relqmc_p, taua_p, ixindx_p
         end subroutine rrtmg_lw_inatm_codon
      end interface

      call rrtmg_lw_inatm_select_impl()
      if (.not. use_native_rrtmg_lw_inatm_impl) then
         call rrtmg_lw_inatm_log_entered()
         inflag64 = 0_c_int64_t
         iceflag64 = 0_c_int64_t
         liqflag64 = 0_c_int64_t
         do ix = 1, maxxsec
            ixindx64(ix) = int(ixindx(ix), c_int64_t)
         enddo
         call rrtmg_lw_inatm_codon( &
              int(iplon, c_int64_t), int(nlay, c_int64_t), int(size(play,1), c_int64_t), &
              int(icld, c_int64_t), int(iaer, c_int64_t), int(nbndlw, c_int64_t), &
              int(ngptlw, c_int64_t), int(nmol, c_int64_t), int(maxxsec, c_int64_t), &
              int(mxmol, c_int64_t), real(grav, c_double), real(avogad, c_double), &
              transfer(loc(play(1,1)), c_null_ptr), transfer(loc(plev(1,1)), c_null_ptr), &
              transfer(loc(tlay(1,1)), c_null_ptr), transfer(loc(tlev(1,1)), c_null_ptr), &
              transfer(loc(tsfc(1)), c_null_ptr), transfer(loc(h2ovmr(1,1)), c_null_ptr), &
              transfer(loc(o3vmr(1,1)), c_null_ptr), transfer(loc(co2vmr(1,1)), c_null_ptr), &
              transfer(loc(ch4vmr(1,1)), c_null_ptr), transfer(loc(o2vmr(1,1)), c_null_ptr), &
              transfer(loc(n2ovmr(1,1)), c_null_ptr), transfer(loc(cfc11vmr(1,1)), c_null_ptr), &
              transfer(loc(cfc12vmr(1,1)), c_null_ptr), transfer(loc(cfc22vmr(1,1)), c_null_ptr), &
              transfer(loc(ccl4vmr(1,1)), c_null_ptr), transfer(loc(emis(1,1)), c_null_ptr), &
              int(inflglw, c_int64_t), int(iceflglw, c_int64_t), int(liqflglw, c_int64_t), &
              transfer(loc(cldfmcl(1,1,1)), c_null_ptr), transfer(loc(taucmcl(1,1,1)), c_null_ptr), &
              transfer(loc(ciwpmcl(1,1,1)), c_null_ptr), transfer(loc(clwpmcl(1,1,1)), c_null_ptr), &
              transfer(loc(reicmcl(1,1)), c_null_ptr), transfer(loc(relqmcl(1,1)), c_null_ptr), &
              transfer(loc(tauaer(1,1,1)), c_null_ptr), transfer(loc(pavel(1)), c_null_ptr), &
              transfer(loc(pz(0)), c_null_ptr), transfer(loc(tavel(1)), c_null_ptr), &
              transfer(loc(tz(0)), c_null_ptr), transfer(loc(tbound), c_null_ptr), &
              transfer(loc(semiss(1)), c_null_ptr), transfer(loc(coldry(1)), c_null_ptr), &
              transfer(loc(wbrodl(1)), c_null_ptr), transfer(loc(wkl(1,1)), c_null_ptr), &
              transfer(loc(wx(1,1)), c_null_ptr), transfer(loc(pwvcm), c_null_ptr), &
              c_loc(inflag64), c_loc(iceflag64), c_loc(liqflag64), &
              transfer(loc(cldfmc(1,1)), c_null_ptr), transfer(loc(taucmc(1,1)), c_null_ptr), &
              transfer(loc(ciwpmc(1,1)), c_null_ptr), transfer(loc(clwpmc(1,1)), c_null_ptr), &
              transfer(loc(reicmc(1)), c_null_ptr), transfer(loc(dgesmc(1)), c_null_ptr), &
              transfer(loc(relqmc(1)), c_null_ptr), transfer(loc(taua(1,1)), c_null_ptr), &
              c_loc(ixindx64(1)) &
         )
         if (icld .ge. 1) then
            inflag = int(inflag64)
            iceflag = int(iceflag64)
            liqflag = int(liqflag64)
         endif
         return
      endif

!  Initialize all molecular amounts and cloud properties to zero here, then pass input amounts
!  into RRTM arrays below.

      wkl(:,:) = 0.0_r8
      wx(:,:) = 0.0_r8
      cldfmc(:,:) = 0.0_r8
      taucmc(:,:) = 0.0_r8
      ciwpmc(:,:) = 0.0_r8
      clwpmc(:,:) = 0.0_r8
      reicmc(:) = 0.0_r8
      dgesmc(:) = 0.0_r8
      relqmc(:) = 0.0_r8
      taua(:,:) = 0.0_r8
      amttl = 0.0_r8
      wvttl = 0.0_r8
 
!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_LW arrays for pressure, temperature,
!  and molecular amounts.  
!  Pressures are input in mb, or are converted to mb here.
!  Molecular amounts are input in volume mixing ratio, or are converted from 
!  mass mixing ratio (or specific humidity for h2o) to volume mixing ratio
!  here. These are then converted to molecular amount (molec/cm2) below.  
!  The dry air column COLDRY (in molec/cm2) is calculated from the level 
!  pressures, pz (in mb), based on the hydrostatic equation and includes a 
!  correction to account for h2o in the layer.  The molecular weight of moist 
!  air (amm) is calculated for each layer.  
!  Note: In RRTMG, layer indexing goes from bottom to top, and coding below
!  assumes GCM input fields are also bottom to top. Input layer indexing
!  from GCM fields should be reversed here if necessary.

      pz(0) = plev(iplon,nlay+1)
      tz(0) = tlev(iplon,nlay+1)
      do l = 1, nlay
         pavel(l) = play(iplon,nlay-l+1)
         tavel(l) = tlay(iplon,nlay-l+1)
         pz(l) = plev(iplon,nlay-l+1)
         tz(l) = tlev(iplon,nlay-l+1)
! For h2o input in vmr:
         wkl(1,l) = h2ovmr(iplon,nlay-l+1)
! For h2o input in mmr:
!         wkl(1,l) = h2o(iplon,nlay-l)*amdw
! For h2o input in specific humidity;
!         wkl(1,l) = (h2o(iplon,nlay-l)/(1._r8 - h2o(iplon,nlay-l)))*amdw
         wkl(2,l) = co2vmr(iplon,nlay-l+1)
         wkl(3,l) = o3vmr(iplon,nlay-l+1)
         wkl(4,l) = n2ovmr(iplon,nlay-l+1)
         wkl(6,l) = ch4vmr(iplon,nlay-l+1)
         wkl(7,l) = o2vmr(iplon,nlay-l+1)

         amm = (1._r8 - wkl(1,l)) * amd + wkl(1,l) * amw            

         coldry(l) = (pz(l-1)-pz(l)) * 1.e3_r8 * avogad / &
                     (1.e2_r8 * grav * amm * (1._r8 + wkl(1,l)))

! Set cross section molecule amounts from input; convert to vmr if necessary
         wx(1,l) = ccl4vmr(iplon,nlay-l+1)
         wx(2,l) = cfc11vmr(iplon,nlay-l+1)
         wx(3,l) = cfc12vmr(iplon,nlay-l+1)
         wx(4,l) = cfc22vmr(iplon,nlay-l+1)

      enddo

      coldry(nlay) = (pz(nlay-1)) * 1.e3_r8 * avogad / &
                        (1.e2_r8 * grav * amm * (1._r8 + wkl(1,nlay-1)))

! At this point all molecular amounts in wkl and wx are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  also, compute precipitable
! water vapor for diffusivity angle adjustments in rtrn and rtrnmr.

      do l = 1, nlay
         summol = 0.0_r8
         do imol = 2, nmol
            summol = summol + wkl(imol,l)
         enddo
         wbrodl(l) = coldry(l) * (1._r8 - summol)
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
         amttl = amttl + coldry(l)+wkl(1,l)
         wvttl = wvttl + wkl(1,l)
         do ix = 1,maxxsec
            if (ixindx(ix) .ne. 0) then
               wx(ixindx(ix),l) = coldry(l) * wx(ix,l) * 1.e-20_r8
            endif
         enddo
      enddo

      wvsh = (amw * wvttl) / (amd * amttl)
      pwvcm = wvsh * (1.e3_r8 * pz(0)) / (1.e2_r8 * grav)

! Set spectral surface emissivity for each longwave band.  

      do n=1,nbndlw
         semiss(n) = emis(iplon,n)
!          semiss(n) = 1.0_r8
      enddo

! Transfer aerosol optical properties to RRTM variable;
! modify to reverse layer indexing here if necessary.

      if (iaer .ge. 1) then 
         do l = 1, nlay-1
            do ib = 1, nbndlw
               taua(l,ib) = tauaer(iplon,nlay-l,ib)
            enddo
         enddo
      endif

! Transfer cloud fraction and cloud optical properties to RRTM variables,
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflglw
         iceflag = iceflglw
         liqflag = liqflglw

! Move incoming GCM cloud arrays to RRTMG cloud arrays.
! For GCM input, incoming reice is in effective radius; for Fu parameterization (iceflag = 3)
! convert effective radius to generalized effective size using method of Mitchell, JAS, 2002:

         do l = 1, nlay-1
            do ig = 1, ngptlw
               cldfmc(ig,l) = cldfmcl(ig,iplon,nlay-l)
               taucmc(ig,l) = taucmcl(ig,iplon,nlay-l)
               ciwpmc(ig,l) = ciwpmcl(ig,iplon,nlay-l)
               clwpmc(ig,l) = clwpmcl(ig,iplon,nlay-l)
            enddo
            reicmc(l) = reicmcl(iplon,nlay-l)
            if (iceflag .eq. 3) then
               dgesmc(l) = 1.5396_r8 * reicmcl(iplon,nlay-l)
            endif
            relqmc(l) = relqmcl(iplon,nlay-l)
         enddo

! If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.

         cldfmc(:,nlay) = 0.0_r8
         taucmc(:,nlay) = 0.0_r8
         ciwpmc(:,nlay) = 0.0_r8
         clwpmc(:,nlay) = 0.0_r8
         reicmc(nlay) = 0.0_r8
         dgesmc(nlay) = 0.0_r8
         relqmc(nlay) = 0.0_r8
         taua(nlay,:) = 0.0_r8

      endif
      
      end subroutine inatm

      end module rrtmg_lw_rad
