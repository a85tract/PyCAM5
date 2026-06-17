!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.6 $
!     created:   $Date: 2008/01/03 21:35:35 $
!

       module rrtmg_sw_rad

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
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  miacono@aer.com                             *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people:  Steven J. Taubman, Patrick D. Brown,            *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

! --------- Modules ---------

      use shr_kind_mod, only: r8 => shr_kind_r8

!      use parkind, only : jpim, jprb
      use rrsw_vsn
      use mcica_subcol_gen_sw, only: mcica_subcol_sw
      use rrtmg_sw_cldprop, only: cldprop_sw
      use rrtmg_sw_cldprmc, only: cldprmc_sw
! Move call to rrtmg_sw_ini and following use association to 
! GCM initialization area
!      use rrtmg_sw_init, only: rrtmg_sw_ini
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvmc, only: spcvmc_sw

      use perf_mod

      implicit none

      logical :: use_native_rrtmg_sw_rad_pack_impl = .false.
      logical :: rrtmg_sw_rad_pack_impl_selected = .false.
      logical :: rrtmg_sw_rad_pack_entered_logged = .false.
      logical :: use_native_rrtmg_sw_inatm_impl = .false.
      logical :: rrtmg_sw_inatm_impl_selected = .false.
      logical :: rrtmg_sw_inatm_entered_logged = .false.

! public interfaces/functions/subroutines
!      public :: rrtmg_sw, inatm_sw, earth_sun
      public :: rrtmg_sw

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_sw &
            (lchnk   ,ncol    ,nlay    ,icld    ,          &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  , &
             asdir   ,asdif   ,aldir   ,aldif   , &
             coszen  ,adjes   ,dyofyr  ,solvar, &
             inflgsw ,iceflgsw,liqflgsw, &
             cldfmcl ,taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl, &
             ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  ,ssaaer  ,asmaer  , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
             dirdnuv, dirdnir, difdnuv, difdnir, ninflx, ninflxc, &
             swuflxs, swdflxs)


! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for 
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.
!
! Note: The call to RRTMG_SW_INI should be moved to the GCM initialization 
!  area, since this has to be called only once. 
!
! This routine
!    b) calls INATM_SW to read in the atmospheric profile;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPRMC_SW to set cloud optical depth for McICA based
!       on input cloud properties
!    d) calls SETCOEF_SW to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls SPCVMC to call the two-stream model that in turn 
!       calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands and to perform the radiative transfer
!       using McICA, the Monte-Carlo Independent Column Approximation,
!       to represent sub-grid scale cloud variability
!    f) passes the calculated fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_sw.nomcica.f90 (to not use
!     McICA) or rrtmg_sw.f90 (to use McICA) to interface with a GCM.
!
!    1) Standard, single forward model calculation (imca = 0); this is 
!       valid only for clear sky or fully overcast clouds
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!       This method is valid for clear sky or partial cloud conditions.
!
! This call to RRTMG_SW must be preceeded by a call to the module
!     mcica_subcol_gen_sw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngptsw) dimension.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflag, iceflag and liqflag; see text file rrtmg_sw_instructions
!     and subroutine rrtmg_sw_cldprop.f90 for further details):
!
!    1) Input cloud fraction, cloud optical depth, single scattering albedo 
!       and asymmetry parameter directly (inflgsw = 0)
!    2) Input cloud fraction and cloud physical properties: ice fracion,
!       ice and liquid particle sizes (inflgsw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflgsw and liqflgsw
!
! Two methods of aerosol property input are possible:
!     Aerosol properties can be input in one of two ways (controlled by input 
!     flag iaer, see text file rrtmg_sw_instructions for further details):
!
!    1) Input aerosol optical depth, single scattering albedo and asymmetry
!       parameter directly by layer and spectral band (iaer=10)
!    2) Input aerosol optical depth and 0.55 micron directly by layer and use
!       one or more of six ECMWF aerosol types (iaer=6)
!
!
! ------- Modifications -------
!
! This version of RRTMG_SW has been modified from RRTM_SW to use a reduced
! set of g-point intervals and a two-stream model for application to GCMs. 
!
!-- Original version (derived from RRTM_SW)
!     2002: AER. Inc.
!-- Conversion to F90 formatting; addition of 2-stream radiative transfer
!     Feb 2003: J.-J. Morcrette, ECMWF
!-- Additional modifications for GCM application
!     Aug 2003: M. J. Iacono, AER Inc.
!-- Total number of g-points reduced from 224 to 112.  Original
!   set of 224 can be restored by exchanging code in module parrrsw.f90 
!   and in file rrtmg_sw_init.f90.
!     Apr 2004: M. J. Iacono, AER, Inc.
!-- Modifications to include output for direct and diffuse 
!   downward fluxes.  There are output as "true" fluxes without
!   any delta scaling applied.  Code can be commented to exclude
!   this calculation in source file rrtmg_sw_spcvrt.f90.
!     Jan 2005: E. J. Mlawer, M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability.
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays. 
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to output direct and diffuse fluxes either with or without
!   delta scaling based on setting of idelm flag
!     Dec 2008: M. J. Iacono, AER, Inc.

! --------- Modules ---------

      use parrrsw, only : nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2
      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi
      use rrsw_wvn, only : wavenum1, wavenum2

! ------- Declarations

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
      real(kind=r8), intent(in) :: asdir(:)             ! UV/vis surface albedo direct rad
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: aldir(:)             ! Near-IR surface albedo direct rad
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: asdif(:)             ! UV/vis surface albedo: diffuse rad
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: aldif(:)             ! Near-IR surface albedo: diffuse rad
                                                        !    Dimensions: (ncol)

      integer, intent(in) :: dyofyr                     ! Day of the year (used to get Earth/Sun
                                                        !  distance if adjflx not provided)
      real(kind=r8), intent(in) :: adjes                ! Flux adjustment for Earth/Sun distance
      real(kind=r8), intent(in) :: coszen(:)            ! Cosine of solar zenith angle
                                                        !    Dimensions: (ncol)
      real(kind=r8), intent(in) :: solvar(1:nbndsw)     ! Solar constant (Wm-2) scaling per band

      integer, intent(in) :: inflgsw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflgsw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflgsw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfmcl(:,:,:)       ! Cloud fraction
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: taucmcl(:,:,:)       ! Cloud optical depth
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: ssacmcl(:,:,:)       ! Cloud single scattering albedo
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: asmcmcl(:,:,:)       ! Cloud asymmetry parameter
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: fsfcmcl(:,:,:)       ! Cloud forward scattering parameter
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: ciwpmcl(:,:,:)       ! Cloud ice water path (g/m2)
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: clwpmcl(:,:,:)       ! Cloud liquid water path (g/m2)
                                                        !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: reicmcl(:,:)         ! Cloud ice effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: relqmcl(:,:)         ! Cloud water drop effective radius (microns)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tauaer(:,:,:)        ! Aerosol optical depth (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=r8), intent(in) :: ssaaer(:,:,:)        ! Aerosol single scattering albedo (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
      real(kind=r8), intent(in) :: asmaer(:,:,:)        ! Aerosol asymmetry parameter (iaer=10 only)
                                                        !    Dimensions: (ncol,nlay,nbndsw)
                                                        ! (non-delta scaled)      
!      real(kind=r8), intent(in) :: ecaer(:,:,:)         ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                                        !    Dimensions: (ncol,nlay,naerec)
                                                        ! (non-delta scaled)      

! ----- Output -----

      real(kind=r8), target, intent(out) :: swuflx(:,:) ! Total sky shortwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), target, intent(out) :: swdflx(:,:) ! Total sky shortwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), target, intent(out) :: swhr(:,:)   ! Total sky shortwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)
      real(kind=r8), target, intent(out) :: swuflxc(:,:)! Clear sky shortwave upward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), target, intent(out) :: swdflxc(:,:)! Clear sky shortwave downward flux (W/m2)
                                                        !    Dimensions: (ncol,nlay+1)
      real(kind=r8), target, intent(out) :: swhrc(:,:)  ! Clear sky shortwave radiative heating rate (K/d)
                                                        !    Dimensions: (ncol,nlay)

      real(kind=r8), target, intent(out) :: dirdnuv(:,:)! Direct downward shortwave flux, UV/vis
      real(kind=r8), target, intent(out) :: difdnuv(:,:)! Diffuse downward shortwave flux, UV/vis
      real(kind=r8), target, intent(out) :: dirdnir(:,:)! Direct downward shortwave flux, near-IR
      real(kind=r8), target, intent(out) :: difdnir(:,:)! Diffuse downward shortwave flux, near-IR

      real(kind=r8), target, intent(out) :: ninflx(:,:) ! Net shortwave flux, near-IR
      real(kind=r8), target, intent(out) :: ninflxc(:,:)! Net clear sky shortwave flux, near-IR

      real(kind=r8), target, intent(out) :: swuflxs(:,:,:)! shortwave spectral flux up
      real(kind=r8), target, intent(out) :: swdflxs(:,:,:)! shortwave spectral flux down

! ----- Local -----

! Control
      integer :: istart                         ! beginning band of calculation
      integer :: iend                           ! ending band of calculation
      integer :: icpr                           ! cldprop/cldprmc use flag
      integer :: iout = 0                       ! output option flag (inactive)
      integer :: iaer                           ! aerosol option flag
      integer :: idelm                          ! delta-m scaling flag
                                                ! [0 = direct and diffuse fluxes are unscaled]
                                                ! [1 = direct and diffuse fluxes are scaled]
                                                ! (total downward fluxes are always delta scaled)
      integer :: isccos                         ! instrumental cosine response flag (inactive)
      integer :: iplon                          ! column loop index
      integer :: i                              ! layer loop index                       ! jk
      integer :: ib                             ! band loop index                        ! jsw
      integer :: ia, ig                         ! indices
      integer :: k                              ! layer loop index
      integer :: ims                            ! value for changing mcica permute seed
      integer :: imca                           ! flag for mcica [0=off, 1=on]

      real(kind=r8) :: zepsec, zepzen           ! epsilon
      real(kind=r8) :: zdpgcp                   ! flux to heating conversion ratio

! Atmosphere
      real(kind=r8) :: pavel(nlay)            ! layer pressures (mb) 
      real(kind=r8) :: tavel(nlay)            ! layer temperatures (K)
      real(kind=r8) :: pz(0:nlay)             ! level (interface) pressures (hPa, mb)
      real(kind=r8) :: tz(0:nlay)             ! level (interface) temperatures (K)
      real(kind=r8) :: tbound                   ! surface temperature (K)
      real(kind=r8), target :: pdp(nlay)      ! layer pressure thickness (hPa, mb)
      real(kind=r8) :: coldry(nlay)           ! dry air column amount
      real(kind=r8) :: wkl(mxmol,nlay)        ! molecular amounts (mol/cm-2)

!      real(kind=r8) :: earth_sun               ! function for Earth/Sun distance factor
      real(kind=r8) :: cossza                   ! Cosine of solar zenith angle
      real(kind=r8) :: adjflux(jpband)          ! adjustment for current Earth/Sun distance
!      real(kind=r8) :: solvar(jpband)           ! solar constant scaling factor from rrtmg_sw
                                                !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=r8), target :: albdir(nbndsw)   ! surface albedo, direct          ! zalbp
      real(kind=r8), target :: albdif(nbndsw)   ! surface albedo, diffuse         ! zalbd

      real(kind=r8), target :: taua(nlay,nbndsw)! Aerosol optical depth
      real(kind=r8), target :: ssaa(nlay,nbndsw)! Aerosol single scattering albedo
      real(kind=r8), target :: asma(nlay,nbndsw)! Aerosol asymmetry parameter

! Atmosphere - setcoef
      integer :: laytrop                        ! tropopause layer index
      integer :: layswtch                       ! 
      integer :: laylow                         ! 
      integer :: jp(nlay)                     ! 
      integer :: jt(nlay)                     !
      integer :: jt1(nlay)                    !

      real(kind=r8) :: colh2o(nlay)           ! column amount (h2o)
      real(kind=r8) :: colco2(nlay)           ! column amount (co2)
      real(kind=r8) :: colo3(nlay)            ! column amount (o3)
      real(kind=r8) :: coln2o(nlay)           ! column amount (n2o)
      real(kind=r8) :: colch4(nlay)           ! column amount (ch4)
      real(kind=r8) :: colo2(nlay)            ! column amount (o2)
      real(kind=r8) :: colmol(nlay)           ! column amount
      real(kind=r8) :: co2mult(nlay)          ! column amount 

      integer :: indself(nlay)
      integer :: indfor(nlay)
      real(kind=r8) :: selffac(nlay)
      real(kind=r8) :: selffrac(nlay)
      real(kind=r8) :: forfac(nlay)
      real(kind=r8) :: forfrac(nlay)

      real(kind=r8) :: &                        !
                         fac00(nlay), fac01(nlay), &
                         fac10(nlay), fac11(nlay) 

! Atmosphere/clouds - cldprop
      integer :: ncbands                        ! number of cloud spectral bands
      integer :: inflag                         ! flag for cloud property method
      integer :: iceflag                        ! flag for ice cloud properties
      integer :: liqflag                        ! flag for liquid cloud properties

!      real(kind=r8) :: cldfrac(nlay)            ! layer cloud fraction
!      real(kind=r8) :: tauc(nlay)               ! cloud optical depth (non-delta scaled)
!      real(kind=r8) :: ssac(nlay)               ! cloud single scattering albedo (non-delta scaled)
!      real(kind=r8) :: asmc(nlay)               ! cloud asymmetry parameter (non-delta scaled)
!      real(kind=r8) :: ciwp(nlay)               ! cloud ice water path
!      real(kind=r8) :: clwp(nlay)               ! cloud liquid water path
!      real(kind=r8) :: rei(nlay)                ! cloud ice particle size
!      real(kind=r8) :: rel(nlay)                ! cloud liquid particle size

!      real(kind=r8) :: taucloud(nlay,jpband)    ! cloud optical depth
!      real(kind=r8) :: taucldorig(nlay,jpband)  ! cloud optical depth (non-delta scaled)
!      real(kind=r8) :: ssacloud(nlay,jpband)    ! cloud single scattering albedo
!      real(kind=r8) :: asmcloud(nlay,jpband)    ! cloud asymmetry parameter

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=r8), target :: cldfmc(ngptsw,nlay)! cloud fraction [mcica]
      real(kind=r8) :: ciwpmc(ngptsw,nlay)    ! cloud ice water path [mcica]
      real(kind=r8) :: clwpmc(ngptsw,nlay)    ! cloud liquid water path [mcica]
      real(kind=r8) :: relqmc(nlay)           ! liquid particle size (microns)
      real(kind=r8) :: reicmc(nlay)           ! ice particle effective radius (microns)
      real(kind=r8) :: dgesmc(nlay)           ! ice particle generalized effective size (microns)
      real(kind=r8), target :: taucmc(ngptsw,nlay)! cloud optical depth [mcica]
      real(kind=r8), target :: taormc(ngptsw,nlay)! unscaled cloud optical depth [mcica]
      real(kind=r8), target :: ssacmc(ngptsw,nlay)! cloud single scattering albedo [mcica]
      real(kind=r8), target :: asmcmc(ngptsw,nlay)! cloud asymmetry parameter [mcica]
      real(kind=r8) :: fsfcmc(ngptsw,nlay)    ! cloud forward scattering fraction [mcica]

! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real(kind=r8) :: ztauc(nlay,nbndsw)     ! cloud optical depth
      real(kind=r8) :: ztaucorig(nlay,nbndsw) ! unscaled cloud optical depth
      real(kind=r8) :: zasyc(nlay,nbndsw)     ! cloud asymmetry parameter 
                                                !  (first moment of phase function)
      real(kind=r8) :: zomgc(nlay,nbndsw)     ! cloud single scattering albedo
      real(kind=r8), target :: ztaua(nlay,nbndsw)! total aerosol optical depth
      real(kind=r8), target :: zasya(nlay,nbndsw)! total aerosol asymmetry parameter
      real(kind=r8), target :: zomga(nlay,nbndsw)! total aerosol single scattering albedo

      real(kind=r8), target :: zcldfmc(nlay,ngptsw)! cloud fraction [mcica]
      real(kind=r8), target :: ztaucmc(nlay,ngptsw)! cloud optical depth [mcica]
      real(kind=r8), target :: ztaormc(nlay,ngptsw)! unscaled cloud optical depth [mcica]
      real(kind=r8), target :: zasycmc(nlay,ngptsw)! cloud asymmetry parameter [mcica]
      real(kind=r8), target :: zomgcmc(nlay,ngptsw)! cloud single scattering albedo [mcica]

      real(kind=r8), target :: zbbfu(nlay+2)  ! temporary upward shortwave flux (w/m2)
      real(kind=r8), target :: zbbfd(nlay+2)  ! temporary downward shortwave flux (w/m2)
      real(kind=r8), target :: zbbcu(nlay+2)  ! temporary clear sky upward shortwave flux (w/m2)
      real(kind=r8), target :: zbbcd(nlay+2)  ! temporary clear sky downward shortwave flux (w/m2)
      real(kind=r8), target :: zbbfddir(nlay+2)! temporary downward direct shortwave flux (w/m2)
      real(kind=r8), target :: zbbcddir(nlay+2)! temporary clear sky downward direct shortwave flux (w/m2)
      real(kind=r8), target :: zuvfd(nlay+2)  ! temporary UV downward shortwave flux (w/m2)
      real(kind=r8), target :: zuvcd(nlay+2)  ! temporary clear sky UV downward shortwave flux (w/m2)
      real(kind=r8), target :: zuvfddir(nlay+2)! temporary UV downward direct shortwave flux (w/m2)
      real(kind=r8), target :: zuvcddir(nlay+2)! temporary clear sky UV downward direct shortwave flux (w/m2)
      real(kind=r8), target :: znifd(nlay+2)  ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=r8), target :: znicd(nlay+2)  ! temporary clear sky near-IR downward shortwave flux (w/m2)
      real(kind=r8), target :: znifddir(nlay+2)! temporary near-IR downward direct shortwave flux (w/m2)
      real(kind=r8), target :: znicddir(nlay+2)! temporary clear sky near-IR downward direct shortwave flux (w/m2)
! Added for near-IR flux diagnostic
      real(kind=r8), target :: znifu(nlay+2)  ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=r8), target :: znicu(nlay+2)  ! temporary clear sky near-IR downward shortwave flux (w/m2)

! Optional output fields 
      real(kind=r8), target :: swnflx(nlay+2) ! Total sky shortwave net flux (W/m2)
      real(kind=r8), target :: swnflxc(nlay+2)! Clear sky shortwave net flux (W/m2)
      real(kind=r8), target :: dirdflux(nlay+2)! Direct downward shortwave surface flux
      real(kind=r8), target :: difdflux(nlay+2)! Diffuse downward shortwave surface flux
      real(kind=r8), target :: uvdflx(nlay+2)! Total sky downward shortwave flux, UV/vis
      real(kind=r8), target :: nidflx(nlay+2)! Total sky downward shortwave flux, near-IR
      real(kind=r8), target :: zbbfsu(nbndsw,nlay+2)! temporary upward shortwave flux spectral (w/m2)
      real(kind=r8), target :: zbbfsd(nbndsw,nlay+2)! temporary downward shortwave flux spectral (w/m2)

! Output - inactive
!      real(kind=r8) :: zuvfu(nlay+2)         ! temporary upward UV shortwave flux (w/m2)
!      real(kind=r8) :: zuvfd(nlay+2)         ! temporary downward UV shortwave flux (w/m2)
!      real(kind=r8) :: zuvcu(nlay+2)         ! temporary clear sky upward UV shortwave flux (w/m2)
!      real(kind=r8) :: zuvcd(nlay+2)         ! temporary clear sky downward UV shortwave flux (w/m2)
!      real(kind=r8) :: zvsfu(nlay+2)         ! temporary upward visible shortwave flux (w/m2)
!      real(kind=r8) :: zvsfd(nlay+2)         ! temporary downward visible shortwave flux (w/m2)
!      real(kind=r8) :: zvscu(nlay+2)         ! temporary clear sky upward visible shortwave flux (w/m2)
!      real(kind=r8) :: zvscd(nlay+2)         ! temporary clear sky downward visible shortwave flux (w/m2)
!      real(kind=r8) :: znifu(nlay+2)         ! temporary upward near-IR shortwave flux (w/m2)
!      real(kind=r8) :: znifd(nlay+2)         ! temporary downward near-IR shortwave flux (w/m2)
!      real(kind=r8) :: znicu(nlay+2)         ! temporary clear sky upward near-IR shortwave flux (w/m2)
!      real(kind=r8) :: znicd(nlay+2)         ! temporary clear sky downward near-IR shortwave flux (w/m2)

      interface
         subroutine rrtmg_sw_rad_setup_codon(nlay_c, ngptsw_c, nbndsw_c, icld_c, iaer_c, &
              aldir_i_c, aldif_i_c, asdir_i_c, asdif_i_c, albdir_p, albdif_p, cldfmc_p, &
              taucmc_p, taormc_p, asmcmc_p, ssacmc_p, zcldfmc_p, ztaucmc_p, ztaormc_p, &
              zasycmc_p, zomgcmc_p, taua_p, ssaa_p, asma_p, ztaua_p, zasya_p, zomga_p) &
              bind(c, name="rrtmg_sw_rad_setup_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlay_c, ngptsw_c, nbndsw_c, icld_c, iaer_c
            real(c_double), value :: aldir_i_c, aldif_i_c, asdir_i_c, asdif_i_c
            type(c_ptr), value :: albdir_p, albdif_p, cldfmc_p, taucmc_p, taormc_p
            type(c_ptr), value :: asmcmc_p, ssacmc_p, zcldfmc_p, ztaucmc_p, ztaormc_p
            type(c_ptr), value :: zasycmc_p, zomgcmc_p, taua_p, ssaa_p, asma_p
            type(c_ptr), value :: ztaua_p, zasya_p, zomga_p
         end subroutine rrtmg_sw_rad_setup_codon

         subroutine rrtmg_sw_rad_zero_flux_codon(nlay_c, nbndsw_c, zbbcu_p, zbbcd_p, &
              zbbfu_p, zbbfd_p, zbbcddir_p, zbbfddir_p, zuvcd_p, zuvfd_p, zuvcddir_p, &
              zuvfddir_p, znicd_p, znifd_p, znicddir_p, znifddir_p, znicu_p, znifu_p, &
              zbbfsu_p, zbbfsd_p) bind(c, name="rrtmg_sw_rad_zero_flux_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: nlay_c, nbndsw_c
            type(c_ptr), value :: zbbcu_p, zbbcd_p, zbbfu_p, zbbfd_p, zbbcddir_p
            type(c_ptr), value :: zbbfddir_p, zuvcd_p, zuvfd_p, zuvcddir_p, zuvfddir_p
            type(c_ptr), value :: znicd_p, znifd_p, znicddir_p, znifddir_p, znicu_p
            type(c_ptr), value :: znifu_p, zbbfsu_p, zbbfsd_p
         end subroutine rrtmg_sw_rad_zero_flux_codon

         subroutine rrtmg_sw_rad_store_flux_codon(nlay_c, nbndsw_c, ncol_c, iplon_c, heatfac_c, &
              zbbcu_p, zbbcd_p, zbbfu_p, zbbfd_p, zbbcddir_p, zbbfddir_p, zuvfd_p, &
              zuvfddir_p, znicd_p, znifd_p, znifddir_p, znicu_p, znifu_p, zbbfsu_p, &
              zbbfsd_p, pdp_p, swuflxc_p, swdflxc_p, swuflx_p, swdflx_p, swuflxs_p, &
              swdflxs_p, uvdflx_p, nidflx_p, dirdflux_p, difdflux_p, dirdnuv_p, &
              difdnuv_p, dirdnir_p, difdnir_p, ninflx_p, ninflxc_p, swnflxc_p, &
              swnflx_p, swhrc_p, swhr_p) bind(c, name="rrtmg_sw_rad_store_flux_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: nlay_c, nbndsw_c, ncol_c, iplon_c
            real(c_double), value :: heatfac_c
            type(c_ptr), value :: zbbcu_p, zbbcd_p, zbbfu_p, zbbfd_p, zbbcddir_p
            type(c_ptr), value :: zbbfddir_p, zuvfd_p, zuvfddir_p, znicd_p, znifd_p
            type(c_ptr), value :: znifddir_p, znicu_p, znifu_p, zbbfsu_p, zbbfsd_p
            type(c_ptr), value :: pdp_p, swuflxc_p, swdflxc_p, swuflx_p, swdflx_p
            type(c_ptr), value :: swuflxs_p, swdflxs_p, uvdflx_p, nidflx_p, dirdflux_p
            type(c_ptr), value :: difdflux_p, dirdnuv_p, difdnuv_p, dirdnir_p, difdnir_p
            type(c_ptr), value :: ninflx_p, ninflxc_p, swnflxc_p, swnflx_p, swhrc_p
            type(c_ptr), value :: swhr_p
         end subroutine rrtmg_sw_rad_store_flux_codon
      end interface

! Initializations

      zepsec = 1.e-06_r8
      zepzen = 1.e-10_r8
      oneminus = 1.0_r8 - zepsec
      pi = 2._r8 * asin(1._r8)

      istart = jpb1
      iend = jpb2
      icpr = 0
      ims = 2

! In a GCM with or without McICA, set nlon to the longitude dimension
!
! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation (clear and overcast only)
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
!            (clear, overcast or partial cloud conditions)

! *** This version uses McICA (imca = 1) ***

! Set icld to select of clear or cloud calculation and cloud 
! overlap method (read by subroutine readprof from input file INPUT_RRTM):  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap (McICA only)
! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (icld.lt.0.or.icld.gt.3) icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 6, use six ECMWF aerosol types
!           input aerosol optical depth at 0.55 microns for each aerosol type (ecaer)
! iaer = 10, input total aerosol optical depth, single scattering albedo 
!            and asymmetry parameter (tauaer, ssaaer, asmaer) directly
      iaer = 10

      call rrtmg_sw_rad_pack_select_impl()

! Set idelm to select between delta-M scaled or unscaled output direct and diffuse fluxes
! NOTE: total downward fluxes are always delta scaled
! idelm = 0, output direct and diffuse flux components are not delta scaled
!            (direct flux does not include forward scattering peak)
! idelm = 1, output direct and diffuse flux components are delta scaled (default)
!            (direct flux includes part or most of forward scattering peak)
      idelm = 1

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 224 to 112 for input absorption
! coefficient data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_sw_ini

! This is the main longitude/column loop in RRTMG.
! Modify to loop over all columns (nlon) or over daylight columns

      do iplon = 1, ncol

! Prepare atmosphere profile from GCM for use in RRTMG, and define
! other input parameters

         call inatm_sw (iplon, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, &
              h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, adjes, dyofyr, solvar, &
              inflgsw, iceflgsw, liqflgsw, &
              cldfmcl, taucmcl, ssacmcl, asmcmcl, fsfcmcl, ciwpmcl, clwpmcl, &
              reicmcl, relqmcl, tauaer, ssaaer, asmaer, &
              pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
              adjflux, inflag, iceflag, liqflag, cldfmc, taucmc, &
              ssacmc, asmcmc, fsfcmc, ciwpmc, clwpmc, reicmc, dgesmc, relqmc, &
              taua, ssaa, asma)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprop.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprop.  

         call cldprmc_sw(nlay, inflag, iceflag, liqflag, cldfmc, &
                         ciwpmc, clwpmc, reicmc, dgesmc, relqmc, &
                         taormc, taucmc, ssacmc, asmcmc, fsfcmc)
         icpr = 1

! Calculate coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients by interpolating data from stored
! reference atmospheres.

         call setcoef_sw(nlay, pavel, tavel, pz, tz, tbound, coldry, wkl, &
                         laytrop, layswtch, laylow, jp, jt, jt1, &
                         co2mult, colch4, colco2, colh2o, colmol, coln2o, &
                         colo2, colo3, fac00, fac01, fac10, fac11, &
                         selffac, selffrac, indself, forfac, forfrac, indfor)

! Cosine of the solar zenith angle 
!  Prevent using value of zero; ideally, SW model is not called from host model when sun 
!  is below horizon

         cossza = coszen(iplon)
         if (cossza .lt. zepzen) cossza = zepzen

         if (use_native_rrtmg_sw_rad_pack_impl) then

! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer 

! Surface albedo
!  Near-IR bands 16-24 and 29 (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
!         do ib=1,9
         do ib=1,8
            albdir(ib) = aldir(iplon)
            albdif(ib) = aldif(iplon)
         enddo
         albdir(nbndsw) = aldir(iplon)
         albdif(nbndsw) = aldif(iplon)
!  Set band 24 (or, band 9 counting from 1) to use linear average of UV/visible
!  and near-IR values, since this band straddles 0.7 microns: 
         albdir(9) = 0.5*(aldir(iplon) + asdir(iplon))
         albdif(9) = 0.5*(aldif(iplon) + asdif(iplon))
!  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
         do ib=10,13
            albdir(ib) = asdir(iplon)
            albdif(ib) = asdif(iplon)
         enddo


! Clouds
         if (icld.eq.0) then

            zcldfmc(:,:) = 0._r8
            ztaucmc(:,:) = 0._r8
            ztaormc(:,:) = 0._r8
            zasycmc(:,:) = 0._r8
            zomgcmc(:,:) = 1._r8

         elseif (icld.ge.1) then
            do i=1,nlay
               do ig=1,ngptsw
                  zcldfmc(i,ig) = cldfmc(ig,i)
                  ztaucmc(i,ig) = taucmc(ig,i)
                  ztaormc(i,ig) = taormc(ig,i)
                  zasycmc(i,ig) = asmcmc(ig,i)
                  zomgcmc(i,ig) = ssacmc(ig,i)
               enddo
            enddo

         endif   

! Aerosol
! IAER = 0: no aerosols
         if (iaer.eq.0) then

            ztaua(:,:) = 0._r8
            zasya(:,:) = 0._r8
            zomga(:,:) = 1._r8

! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer), 
! or set manually here for each aerosol and layer.
         elseif (iaer.eq.6) then

!            do i = 1, nlay
!               do ia = 1, naerec
!                  ecaer(iplon,i,ia) = 1.0e-15_r8
!               enddo
!            enddo

!            do i = 1, nlay
!               do ib = 1, nbndsw
!                  ztaua(i,ib) = 0._r8
!                  zasya(i,ib) = 0._r8
!                  zomga(i,ib) = 1._r8
!                  do ia = 1, naerec
!                     ztaua(i,ib) = ztaua(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia)
!                     zomga(i,ib) = zomga(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
!                                   rsrpiza(ib,ia)
!                     zasya(i,ib) = zasya(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
!                                   rsrpiza(ib,ia) * rsrasya(ib,ia)
!                  enddo
!                  if (zomga(i,ib) /= 0._r8) then
!                     zasya(i,ib) = zasya(i,ib) / zomga(i,ib)
!                  endif
!                  if (ztaua(i,ib) /= 0._r8) then
!                     zomga(i,ib) = zomga(i,ib) / ztaua(i,ib)
!                  endif
!               enddo
!            enddo

! IAER=10: Direct specification of aerosol optical properties from GCM
         elseif (iaer.eq.10) then

            do i = 1 ,nlay
               do ib = 1 ,nbndsw
                  ztaua(i,ib) = taua(i,ib)
                  zasya(i,ib) = asma(i,ib)
                  zomga(i,ib) = ssaa(i,ib)
               enddo
            enddo

         endif

         else
            call rrtmg_sw_rad_pack_log_entered()
            call rrtmg_sw_rad_setup_codon( &
                 int(nlay, c_int64_t), int(ngptsw, c_int64_t), int(nbndsw, c_int64_t), &
                 int(icld, c_int64_t), int(iaer, c_int64_t), real(aldir(iplon), c_double), &
                 real(aldif(iplon), c_double), real(asdir(iplon), c_double), &
                 real(asdif(iplon), c_double), c_loc(albdir(1)), c_loc(albdif(1)), &
                 c_loc(cldfmc(1,1)), c_loc(taucmc(1,1)), c_loc(taormc(1,1)), &
                 c_loc(asmcmc(1,1)), c_loc(ssacmc(1,1)), c_loc(zcldfmc(1,1)), &
                 c_loc(ztaucmc(1,1)), c_loc(ztaormc(1,1)), c_loc(zasycmc(1,1)), &
                 c_loc(zomgcmc(1,1)), c_loc(taua(1,1)), c_loc(ssaa(1,1)), &
                 c_loc(asma(1,1)), c_loc(ztaua(1,1)), c_loc(zasya(1,1)), &
                 c_loc(zomga(1,1)) &
            )
         endif


! Call the 2-stream radiation transfer model

         if (use_native_rrtmg_sw_rad_pack_impl) then
            do i=1,nlay+1
               zbbcu(i) = 0._r8
               zbbcd(i) = 0._r8
               zbbfu(i) = 0._r8
               zbbfd(i) = 0._r8
               zbbcddir(i) = 0._r8
               zbbfddir(i) = 0._r8
               zuvcd(i) = 0._r8
               zuvfd(i) = 0._r8
               zuvcddir(i) = 0._r8
               zuvfddir(i) = 0._r8
               znicd(i) = 0._r8
               znifd(i) = 0._r8
               znicddir(i) = 0._r8
               znifddir(i) = 0._r8
               znicu(i) = 0._r8
               znifu(i) = 0._r8
               zbbfsu(:,i) = 0._r8
               zbbfsd(:,i) = 0._r8
            enddo
         else
            call rrtmg_sw_rad_zero_flux_codon( &
                 int(nlay, c_int64_t), int(nbndsw, c_int64_t), c_loc(zbbcu(1)), &
                 c_loc(zbbcd(1)), c_loc(zbbfu(1)), c_loc(zbbfd(1)), c_loc(zbbcddir(1)), &
                 c_loc(zbbfddir(1)), c_loc(zuvcd(1)), c_loc(zuvfd(1)), c_loc(zuvcddir(1)), &
                 c_loc(zuvfddir(1)), c_loc(znicd(1)), c_loc(znifd(1)), c_loc(znicddir(1)), &
                 c_loc(znifddir(1)), c_loc(znicu(1)), c_loc(znifu(1)), c_loc(zbbfsu(1,1)), &
                 c_loc(zbbfsd(1,1)) &
            )
         endif

         call spcvmc_sw &
             (lchnk, iplon, nlay, istart, iend, icpr, idelm, iout, &
              pavel, tavel, pz, tz, tbound, albdif, albdir, &
              zcldfmc, ztaucmc, zasycmc, zomgcmc, ztaormc, &
              ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &	 
              laytrop, layswtch, laylow, jp, jt, jt1, &
              co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
              fac00, fac01, fac10, fac11, &
              selffac, selffrac, indself, forfac, forfrac, indfor, &
              zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd, znifu, znicu, &
              zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir, zbbfsu, zbbfsd)

! Transfer up and down, clear and total sky fluxes to output arrays.
! Vertical indexing goes from bottom to top

         if (use_native_rrtmg_sw_rad_pack_impl) then
            do i = 1, nlay+1
               swuflxc(iplon,i) = zbbcu(i)
               swdflxc(iplon,i) = zbbcd(i)
               swuflx(iplon,i) = zbbfu(i)
               swdflx(iplon,i) = zbbfd(i)
               swuflxs(:,iplon,i) = zbbfsu(:,i)
               swdflxs(:,iplon,i) = zbbfsd(:,i)
               uvdflx(i) = zuvfd(i)
               nidflx(i) = znifd(i)
!  Direct/diffuse fluxes
               dirdflux(i) = zbbfddir(i)
               difdflux(i) = swdflx(iplon,i) - dirdflux(i)
!  UV/visible direct/diffuse fluxes
               dirdnuv(iplon,i) = zuvfddir(i)
               difdnuv(iplon,i) = zuvfd(i) - dirdnuv(iplon,i)
!  Near-IR direct/diffuse fluxes
               dirdnir(iplon,i) = znifddir(i)
               difdnir(iplon,i) = znifd(i) - dirdnir(iplon,i)
!  Added for net near-IR diagnostic
               ninflx(iplon,i) = znifd(i) - znifu(i)
               ninflxc(iplon,i) = znicd(i) - znicu(i)
            enddo

!  Total and clear sky net fluxes
            do i = 1, nlay+1
               swnflxc(i) = swdflxc(iplon,i) - swuflxc(iplon,i)
               swnflx(i) = swdflx(iplon,i) - swuflx(iplon,i)
            enddo

!  Total and clear sky heating rates
!  Heating units are in K/d. Flux units are in W/m2.
            do i = 1, nlay
               zdpgcp = heatfac / pdp(i)
               swhrc(iplon,i) = (swnflxc(i+1) - swnflxc(i)) * zdpgcp
               swhr(iplon,i) = (swnflx(i+1) - swnflx(i)) * zdpgcp
            enddo
            swhrc(iplon,nlay) = 0._r8
            swhr(iplon,nlay) = 0._r8
         else
            call rrtmg_sw_rad_store_flux_codon( &
                 int(nlay, c_int64_t), int(nbndsw, c_int64_t), int(size(swuflx,1), c_int64_t), &
                 int(iplon, c_int64_t), real(heatfac, c_double), c_loc(zbbcu(1)), &
                 c_loc(zbbcd(1)), c_loc(zbbfu(1)), c_loc(zbbfd(1)), c_loc(zbbcddir(1)), &
                 c_loc(zbbfddir(1)), c_loc(zuvfd(1)), c_loc(zuvfddir(1)), c_loc(znicd(1)), &
                 c_loc(znifd(1)), c_loc(znifddir(1)), c_loc(znicu(1)), c_loc(znifu(1)), &
                 c_loc(zbbfsu(1,1)), c_loc(zbbfsd(1,1)), c_loc(pdp(1)), &
                 c_loc(swuflxc(1,1)), c_loc(swdflxc(1,1)), c_loc(swuflx(1,1)), &
                 c_loc(swdflx(1,1)), c_loc(swuflxs(1,1,1)), c_loc(swdflxs(1,1,1)), &
                 c_loc(uvdflx(1)), c_loc(nidflx(1)), c_loc(dirdflux(1)), c_loc(difdflux(1)), &
                 c_loc(dirdnuv(1,1)), c_loc(difdnuv(1,1)), c_loc(dirdnir(1,1)), &
                 c_loc(difdnir(1,1)), c_loc(ninflx(1,1)), c_loc(ninflxc(1,1)), &
                 c_loc(swnflxc(1)), c_loc(swnflx(1)), c_loc(swhrc(1,1)), c_loc(swhr(1,1)) &
            )
         endif

! End longitude loop
      enddo

      end subroutine rrtmg_sw

! --------------------------------------------------------------------------
      subroutine rrtmg_sw_rad_pack_select_impl()

      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (rrtmg_sw_rad_pack_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('RRTMG_SW_RAD_PACK_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_rrtmg_sw_rad_pack_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_rrtmg_sw_rad_pack_impl = .false.
      end if

      rrtmg_sw_rad_pack_impl_selected = .true.

      if (masterproc) then
         if (use_native_rrtmg_sw_rad_pack_impl) then
            write(iulog,*) 'rrtmg_sw_rad_pack implementation = native'
         else
            write(iulog,*) 'rrtmg_sw_rad_pack implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine rrtmg_sw_rad_pack_select_impl

! --------------------------------------------------------------------------
      subroutine rrtmg_sw_rad_pack_log_entered()

      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

      if (rrtmg_sw_rad_pack_entered_logged) return
      rrtmg_sw_rad_pack_entered_logged = .true.

      if (masterproc) then
         write(iulog,*) 'rrtmg_sw_rad_pack entered (mcica sw setup/flux transfer = codon)'
         call flush(iulog)
      end if

      end subroutine rrtmg_sw_rad_pack_log_entered

! --------------------------------------------------------------------------
      subroutine rrtmg_sw_inatm_select_impl()

      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (rrtmg_sw_inatm_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('RRTMG_SW_INATM_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_rrtmg_sw_inatm_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_rrtmg_sw_inatm_impl = .false.
      end if

      rrtmg_sw_inatm_impl_selected = .true.

      if (masterproc) then
         if (use_native_rrtmg_sw_inatm_impl) then
            write(iulog,*) 'rrtmg_sw_inatm implementation = native'
         else
            write(iulog,*) 'rrtmg_sw_inatm implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine rrtmg_sw_inatm_select_impl

! --------------------------------------------------------------------------
      subroutine rrtmg_sw_inatm_log_entered()

      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

      if (rrtmg_sw_inatm_entered_logged) return
      rrtmg_sw_inatm_entered_logged = .true.

      if (masterproc) then
         write(iulog,*) 'rrtmg_sw_inatm entered (mcica sw atmosphere packing = codon)'
         call flush(iulog)
      end if

      end subroutine rrtmg_sw_inatm_log_entered

!*************************************************************************
      real(kind=r8) function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

! ------- Modules -------

      use rrsw_con, only : pi

      integer, intent(in) :: idn

      real(kind=r8) :: gamma

      gamma = 2._r8*pi*(idn-1)/365._r8

! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110_r8 + .034221_r8 * cos(gamma) + .001289_r8 * sin(gamma) + &
                   .000719_r8 * cos(2._r8*gamma) + .000077_r8 * sin(2._r8*gamma)

      end function earth_sun

!***************************************************************************
      subroutine inatm_sw (iplon, nlay, icld, iaer, &
            play, plev, tlay, tlev, tsfc, &
            h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, adjes, dyofyr, solvar, &
            inflgsw, iceflgsw, liqflgsw, &
            cldfmcl, taucmcl, ssacmcl, asmcmcl, fsfcmcl, ciwpmcl, clwpmcl, &
            reicmcl, relqmcl, tauaer, ssaaer, asmaer, &
            pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
            adjflux, inflag, iceflag, liqflag, cldfmc, taucmc, &
            ssacmc, asmcmc, fsfcmc, ciwpmc, clwpmc, reicmc, dgesmc, relqmc, &
            taua, ssaa, asma)
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_SW.
!  Set other RRTMG_SW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parrrsw, only : nbndsw, ngptsw, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2, rrsw_scon
      use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
      use rrsw_con, only : heatfac, oneminus, pi, grav, avogad
      use rrsw_wvn, only : ng, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: iplon                      ! column loop index
      integer, intent(in) :: nlay                       ! number of model layers
      integer, intent(in) :: icld                       ! clear/cloud and cloud overlap flag
      integer, intent(in) :: iaer                       ! aerosol option flag

      real(kind=r8), intent(in) :: play(:,:)            ! Layer pressures (hPa, mb)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: plev(:,:)            ! Interface pressures (hPa, mb)
                                                        ! Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tlay(:,:)            ! Layer temperatures (K)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: tlev(:,:)            ! Interface temperatures (K)
                                                        ! Dimensions: (ncol,nlay+1)
      real(kind=r8), intent(in) :: tsfc(:)              ! Surface temperature (K)
                                                        ! Dimensions: (ncol)
      real(kind=r8), intent(in) :: h2ovmr(:,:)          ! H2O volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o3vmr(:,:)           ! O3 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: co2vmr(:,:)          ! CO2 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: ch4vmr(:,:)          ! Methane volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: o2vmr(:,:)           ! O2 volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: n2ovmr(:,:)          ! Nitrous oxide volume mixing ratio
                                                        ! Dimensions: (ncol,nlay)

      integer, intent(in) :: dyofyr                     ! Day of the year (used to get Earth/Sun
                                                        !  distance if adjflx not provided)
      real(kind=r8), intent(in) :: adjes                ! Flux adjustment for Earth/Sun distance
      real(kind=r8), intent(in) :: solvar(jpb1:jpb2)    ! Solar constant (Wm-2) scaling per band

      integer, intent(in) :: inflgsw                    ! Flag for cloud optical properties
      integer, intent(in) :: iceflgsw                   ! Flag for ice particle specification
      integer, intent(in) :: liqflgsw                   ! Flag for liquid droplet specification

      real(kind=r8), intent(in) :: cldfmcl(:,:,:)       ! Cloud fraction
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: taucmcl(:,:,:)       ! Cloud optical depth (optional)
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: ssacmcl(:,:,:)       ! Cloud single scattering albedo
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: asmcmcl(:,:,:)       ! Cloud asymmetry parameter
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: fsfcmcl(:,:,:)       ! Cloud forward scattering fraction
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: ciwpmcl(:,:,:)       ! Cloud ice water path (g/m2)
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: clwpmcl(:,:,:)       ! Cloud liquid water path (g/m2)
                                                        ! Dimensions: (ngptsw,ncol,nlay)
      real(kind=r8), intent(in) :: reicmcl(:,:)         ! Cloud ice effective radius (microns)
                                                        ! Dimensions: (ncol,nlay)
      real(kind=r8), intent(in) :: relqmcl(:,:)         ! Cloud water drop effective radius (microns)
                                                        ! Dimensions: (ncol,nlay)

      real(kind=r8), intent(in) :: tauaer(:,:,:)        ! Aerosol optical depth
                                                        ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=r8), intent(in) :: ssaaer(:,:,:)        ! Aerosol single scattering albedo
                                                        ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=r8), intent(in) :: asmaer(:,:,:)        ! Aerosol asymmetry parameter
                                                        ! Dimensions: (ncol,nlay,nbndsw)

! Atmosphere

      real(kind=r8), intent(out) :: pavel(:)            ! layer pressures (mb) 
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: tavel(:)            ! layer temperatures (K)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: pz(0:)              ! level (interface) pressures (hPa, mb)
                                                        ! Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tz(0:)              ! level (interface) temperatures (K)
                                                        ! Dimensions: (0:nlay)
      real(kind=r8), intent(out) :: tbound              ! surface temperature (K)
      real(kind=r8), intent(out) :: pdp(:)              ! layer pressure thickness (hPa, mb)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: coldry(:)           ! dry air column density (mol/cm2)
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: wkl(:,:)            ! molecular amounts (mol/cm-2)
                                                        ! Dimensions: (mxmol,nlay)

      real(kind=r8), intent(out) :: adjflux(:)          ! adjustment for current Earth/Sun distance
                                                        ! Dimensions: (jpband)
!      real(kind=r8), intent(out) :: solvar(:)           ! solar constant scaling factor from rrtmg_sw
                                                        ! Dimensions: (jpband)
                                                        !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=r8), intent(out) :: taua(:,:)           ! Aerosol optical depth
                                                        ! Dimensions: (nlay,nbndsw)
      real(kind=r8), intent(out) :: ssaa(:,:)           ! Aerosol single scattering albedo
                                                        ! Dimensions: (nlay,nbndsw)
      real(kind=r8), intent(out) :: asma(:,:)           ! Aerosol asymmetry parameter
                                                        ! Dimensions: (nlay,nbndsw)

! Atmosphere/clouds - cldprop
      integer, intent(out) :: inflag                    ! flag for cloud property method
      integer, intent(out) :: iceflag                   ! flag for ice cloud properties
      integer, intent(out) :: liqflag                   ! flag for liquid cloud properties

      real(kind=r8), intent(out) :: cldfmc(:,:)         ! layer cloud fraction
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: taucmc(:,:)         ! cloud optical depth (non-delta scaled)
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: ssacmc(:,:)         ! cloud single scattering albedo (non-delta-scaled)
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: asmcmc(:,:)         ! cloud asymmetry parameter (non-delta scaled)
      real(kind=r8), intent(out) :: fsfcmc(:,:)         ! cloud forward scattering fraction (non-delta scaled)
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: ciwpmc(:,:)         ! cloud ice water path
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: clwpmc(:,:)         ! cloud liquid water path
                                                        ! Dimensions: (ngptsw,nlay)
      real(kind=r8), intent(out) :: reicmc(:)           ! cloud ice particle effective radius
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: dgesmc(:)           ! cloud ice particle effective radius
                                                        ! Dimensions: (nlay)
      real(kind=r8), intent(out) :: relqmc(:)           ! cloud liquid particle size
                                                        ! Dimensions: (nlay)

! ----- Local -----
      real(kind=r8), parameter :: amd = 28.9660_r8      ! Effective molecular weight of dry air (g/mol)
      real(kind=r8), parameter :: amw = 18.0160_r8      ! Molecular weight of water vapor (g/mol)
!      real(kind=r8), parameter :: amc = 44.0098_r8      ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=r8), parameter :: amo = 47.9998_r8      ! Molecular weight of ozone (g/mol)
!      real(kind=r8), parameter :: amo2 = 31.9999_r8     ! Molecular weight of oxygen (g/mol)
!      real(kind=r8), parameter :: amch4 = 16.0430_r8    ! Molecular weight of methane (g/mol)
!      real(kind=r8), parameter :: amn2o = 44.0128_r8    ! Molecular weight of nitrous oxide (g/mol)

! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
      real(kind=r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
      real(kind=r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
      real(kind=r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
      real(kind=r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide

      real(kind=r8), parameter :: sbc = 5.67e-08_r8     ! Stefan-Boltzmann constant (W/m2K4)

      integer :: isp, l, ix, n, imol, ib, ig   ! Loop indices
      real(kind=r8) :: amm, summol                      ! 
      real(kind=r8) :: adjflx                           ! flux adjustment for Earth/Sun distance
      integer(c_int64_t), target :: inflag64, iceflag64, liqflag64

      interface
         subroutine rrtmg_sw_inatm_codon(iplon_c, nlay_c, ldcol_c, icld_c, iaer_c, &
              nbndsw_c, ngptsw_c, nmol_c, mxmol_c, jpband_c, jpb1_c, jpb2_c, &
              grav_c, avogad_c, adjflx_c, play_p, plev_p, tlay_p, tlev_p, &
              tsfc_p, h2ovmr_p, o3vmr_p, co2vmr_p, ch4vmr_p, o2vmr_p, n2ovmr_p, &
              solvar_p, inflgsw_c, iceflgsw_c, liqflgsw_c, cldfmcl_p, taucmcl_p, &
              ssacmcl_p, asmcmcl_p, fsfcmcl_p, ciwpmcl_p, clwpmcl_p, reicmcl_p, &
              relqmcl_p, tauaer_p, ssaaer_p, asmaer_p, pavel_p, pz_p, pdp_p, &
              tavel_p, tz_p, tbound_p, coldry_p, wkl_p, adjflux_p, inflag_p, &
              iceflag_p, liqflag_p, cldfmc_p, taucmc_p, ssacmc_p, asmcmc_p, &
              fsfcmc_p, ciwpmc_p, clwpmc_p, reicmc_p, dgesmc_p, relqmc_p, &
              taua_p, ssaa_p, asma_p) bind(c, name="rrtmg_sw_inatm_codon")
            use iso_c_binding, only: c_double, c_int64_t, c_ptr
            integer(c_int64_t), value :: iplon_c, nlay_c, ldcol_c, icld_c, iaer_c
            integer(c_int64_t), value :: nbndsw_c, ngptsw_c, nmol_c, mxmol_c
            integer(c_int64_t), value :: jpband_c, jpb1_c, jpb2_c
            integer(c_int64_t), value :: inflgsw_c, iceflgsw_c, liqflgsw_c
            real(c_double), value :: grav_c, avogad_c, adjflx_c
            type(c_ptr), value :: play_p, plev_p, tlay_p, tlev_p, tsfc_p
            type(c_ptr), value :: h2ovmr_p, o3vmr_p, co2vmr_p, ch4vmr_p, o2vmr_p
            type(c_ptr), value :: n2ovmr_p, solvar_p, cldfmcl_p, taucmcl_p, ssacmcl_p
            type(c_ptr), value :: asmcmcl_p, fsfcmcl_p, ciwpmcl_p, clwpmcl_p
            type(c_ptr), value :: reicmcl_p, relqmcl_p, tauaer_p, ssaaer_p, asmaer_p
            type(c_ptr), value :: pavel_p, pz_p, pdp_p, tavel_p, tz_p, tbound_p
            type(c_ptr), value :: coldry_p, wkl_p, adjflux_p, inflag_p, iceflag_p
            type(c_ptr), value :: liqflag_p, cldfmc_p, taucmc_p, ssacmc_p, asmcmc_p
            type(c_ptr), value :: fsfcmc_p, ciwpmc_p, clwpmc_p, reicmc_p, dgesmc_p
            type(c_ptr), value :: relqmc_p, taua_p, ssaa_p, asma_p
         end subroutine rrtmg_sw_inatm_codon
      end interface
!      real(kind=r8) :: earth_sun                        ! function for Earth/Sun distance adjustment
!      real(kind=r8) :: solar_band_irrad(jpb1:jpb2) ! rrtmg assumed-solar irradiance in each sw band

      call rrtmg_sw_inatm_select_impl()
      if (.not. use_native_rrtmg_sw_inatm_impl) then
         adjflx = adjes
         if (dyofyr .gt. 0) then
            adjflx = earth_sun(dyofyr)
         endif
         call rrtmg_sw_inatm_log_entered()
         inflag64 = 0_c_int64_t
         iceflag64 = 0_c_int64_t
         liqflag64 = 0_c_int64_t
         call rrtmg_sw_inatm_codon( &
              int(iplon, c_int64_t), int(nlay, c_int64_t), int(size(play,1), c_int64_t), &
              int(icld, c_int64_t), int(iaer, c_int64_t), int(nbndsw, c_int64_t), &
              int(ngptsw, c_int64_t), int(nmol, c_int64_t), int(mxmol, c_int64_t), &
              int(jpband, c_int64_t), int(jpb1, c_int64_t), int(jpb2, c_int64_t), &
              real(grav, c_double), real(avogad, c_double), real(adjflx, c_double), &
              transfer(loc(play(1,1)), c_null_ptr), transfer(loc(plev(1,1)), c_null_ptr), &
              transfer(loc(tlay(1,1)), c_null_ptr), transfer(loc(tlev(1,1)), c_null_ptr), &
              transfer(loc(tsfc(1)), c_null_ptr), transfer(loc(h2ovmr(1,1)), c_null_ptr), &
              transfer(loc(o3vmr(1,1)), c_null_ptr), transfer(loc(co2vmr(1,1)), c_null_ptr), &
              transfer(loc(ch4vmr(1,1)), c_null_ptr), transfer(loc(o2vmr(1,1)), c_null_ptr), &
              transfer(loc(n2ovmr(1,1)), c_null_ptr), transfer(loc(solvar(jpb1)), c_null_ptr), &
              int(inflgsw, c_int64_t), int(iceflgsw, c_int64_t), int(liqflgsw, c_int64_t), &
              transfer(loc(cldfmcl(1,1,1)), c_null_ptr), transfer(loc(taucmcl(1,1,1)), c_null_ptr), &
              transfer(loc(ssacmcl(1,1,1)), c_null_ptr), transfer(loc(asmcmcl(1,1,1)), c_null_ptr), &
              transfer(loc(fsfcmcl(1,1,1)), c_null_ptr), transfer(loc(ciwpmcl(1,1,1)), c_null_ptr), &
              transfer(loc(clwpmcl(1,1,1)), c_null_ptr), transfer(loc(reicmcl(1,1)), c_null_ptr), &
              transfer(loc(relqmcl(1,1)), c_null_ptr), transfer(loc(tauaer(1,1,1)), c_null_ptr), &
              transfer(loc(ssaaer(1,1,1)), c_null_ptr), transfer(loc(asmaer(1,1,1)), c_null_ptr), &
              transfer(loc(pavel(1)), c_null_ptr), transfer(loc(pz(0)), c_null_ptr), &
              transfer(loc(pdp(1)), c_null_ptr), transfer(loc(tavel(1)), c_null_ptr), &
              transfer(loc(tz(0)), c_null_ptr), transfer(loc(tbound), c_null_ptr), &
              transfer(loc(coldry(1)), c_null_ptr), transfer(loc(wkl(1,1)), c_null_ptr), &
              transfer(loc(adjflux(1)), c_null_ptr), c_loc(inflag64), c_loc(iceflag64), &
              c_loc(liqflag64), transfer(loc(cldfmc(1,1)), c_null_ptr), &
              transfer(loc(taucmc(1,1)), c_null_ptr), transfer(loc(ssacmc(1,1)), c_null_ptr), &
              transfer(loc(asmcmc(1,1)), c_null_ptr), transfer(loc(fsfcmc(1,1)), c_null_ptr), &
              transfer(loc(ciwpmc(1,1)), c_null_ptr), transfer(loc(clwpmc(1,1)), c_null_ptr), &
              transfer(loc(reicmc(1)), c_null_ptr), transfer(loc(dgesmc(1)), c_null_ptr), &
              transfer(loc(relqmc(1)), c_null_ptr), transfer(loc(taua(1,1)), c_null_ptr), &
              transfer(loc(ssaa(1,1)), c_null_ptr), transfer(loc(asma(1,1)), c_null_ptr) &
         )
         if (icld .ge. 1) then
            inflag = int(inflag64)
            iceflag = int(iceflag64)
            liqflag = int(liqflag64)
         endif
         return
      endif

!  Initialize all molecular amounts to zero here, then pass input amounts
!  into RRTM array WKL below.

       wkl(:,:) = 0.0_r8
       cldfmc(:,:) = 0.0_r8
       taucmc(:,:) = 0.0_r8
       ssacmc(:,:) = 1.0_r8
       asmcmc(:,:) = 0.0_r8
       fsfcmc(:,:) = 0.0_r8
       ciwpmc(:,:) = 0.0_r8
       clwpmc(:,:) = 0.0_r8
       reicmc(:) = 0.0_r8
       dgesmc(:) = 0.0_r8
       relqmc(:) = 0.0_r8
       taua(:,:) = 0.0_r8
       ssaa(:,:) = 1.0_r8
       asma(:,:) = 0.0_r8
 
! Set flux adjustment for current Earth/Sun distance (two options).
! 1) Use Earth/Sun distance flux adjustment provided by GCM (input as adjes);
      adjflx = adjes
!
! 2) Calculate Earth/Sun distance from DYOFYR, the cumulative day of the year.
!    (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU). 
      if (dyofyr .gt. 0) then
         adjflx = earth_sun(dyofyr)
      endif

! Set incoming solar flux adjustment to include adjustment for
! current Earth/Sun distance (ADJFLX) and scaling of default internal
! solar constant (rrsw_scon = 1368.22 Wm-2) by band (SOLVAR).  SOLVAR can be set 
! to a single scaling factor as needed, or to a different value in each 
! band, which may be necessary for paleoclimate simulations. 
! 

      adjflux(:) = 0._r8
      do ib = jpb1,jpb2
         adjflux(ib) = adjflx * solvar(ib)
      enddo

!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_SW arrays for pressure, temperature,
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
         pdp(l) = pz(l-1) - pz(l)
! For h2o input in vmr:
         wkl(1,l) = h2ovmr(iplon,nlay-l+1)
! For h2o input in mmr:
!         wkl(1,l) = h2o(iplon,nlayers-l)*amdw
! For h2o input in specific humidity;
!         wkl(1,l) = (h2o(iplon,nlayers-l)/(1._r8 - h2o(iplon,nlayers-l)))*amdw
         wkl(2,l) = co2vmr(iplon,nlay-l+1)
         wkl(3,l) = o3vmr(iplon,nlay-l+1)
         wkl(4,l) = n2ovmr(iplon,nlay-l+1)
         wkl(6,l) = ch4vmr(iplon,nlay-l+1)
         wkl(7,l) = o2vmr(iplon,nlay-l+1) 
         amm = (1._r8 - wkl(1,l)) * amd + wkl(1,l) * amw            
         coldry(l) = (pz(l-1)-pz(l)) * 1.e3_r8 * avogad / &
                     (1.e2_r8 * grav * amm * (1._r8 + wkl(1,l)))
      enddo

      coldry(nlay) = (pz(nlay-1)) * 1.e3_r8 * avogad / &
                        (1.e2_r8 * grav * amm * (1._r8 + wkl(1,nlay-1)))

! At this point all molecular amounts in wkl are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  

      do l = 1, nlay
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
      enddo

! Transfer aerosol optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (iaer .ge. 1) then 
         do l = 1, nlay-1
            do ib = 1, nbndsw
               taua(l,ib) = tauaer(iplon,nlay-l,ib)
               ssaa(l,ib) = ssaaer(iplon,nlay-l,ib)
               asma(l,ib) = asmaer(iplon,nlay-l,ib)
            enddo
         enddo
      endif

! Transfer cloud fraction and cloud optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflgsw
         iceflag = iceflgsw
         liqflag = liqflgsw

! Move incoming GCM cloud arrays to RRTMG cloud arrays.
! For GCM input, incoming reice is in effective radius; for Fu parameterization (iceflag = 3)
! convert effective radius to generalized effective size using method of Mitchell, JAS, 2002:

         do l = 1, nlay-1
            do ig = 1, ngptsw
               cldfmc(ig,l) = cldfmcl(ig,iplon,nlay-l)
               taucmc(ig,l) = taucmcl(ig,iplon,nlay-l)
               ssacmc(ig,l) = ssacmcl(ig,iplon,nlay-l)
               asmcmc(ig,l) = asmcmcl(ig,iplon,nlay-l)
               fsfcmc(ig,l) = fsfcmcl(ig,iplon,nlay-l)
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
         ssacmc(:,nlay) = 1.0_r8
         asmcmc(:,nlay) = 0.0_r8
         fsfcmc(:,nlay) = 0.0_r8
         ciwpmc(:,nlay) = 0.0_r8
         clwpmc(:,nlay) = 0.0_r8
         reicmc(nlay) = 0.0_r8
         dgesmc(nlay) = 0.0_r8
         relqmc(nlay) = 0.0_r8
         taua(nlay,:) = 0.0_r8
         ssaa(nlay,:) = 1.0_r8
         asma(nlay,:) = 0.0_r8
     
      endif

      end subroutine inatm_sw

      end module rrtmg_sw_rad
