
module zm_conv

!---------------------------------------------------------------------------------
! Purpose:
!
! Interface from Zhang-McFarlane convection scheme, includes evaporation of convective 
! precip from the ZM scheme
!
! Apr 2006: RBN: Code added to perform a dilute ascent for closure of the CM mass flux
!                based on an entraining plume a la Raymond and Blythe (1992)
!
! Author: Byron Boville, from code in tphysbc
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use spmd_utils,      only: masterproc
  use ppgrid,          only: pcols, pver, pverp
  use cloud_fraction,  only: cldfrc_fice
  use physconst,       only: cpair, epsilo, gravit, latice, latvap, tmelt, rair, &
                             cpwv, cpliq, rh2o
  use cam_abortutils,  only: endrun
  use cam_logfile,     only: iulog
  use perf_mod,        only: t_startf, t_stopf
  use ap_zm_convr_scheme, only: zm_convr_init, zm_convr_run
  use ap_zm_conv_evap_scheme, only: zm_conv_evap_run
  use ap_zm_conv_convtran_scheme, only: zm_conv_convtran_run
  use ap_zm_conv_momtran_scheme, only: zm_conv_momtran_run

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zm_convi                 ! ZM schemea
  public zm_convr                 ! ZM schemea
  public zm_conv_evap             ! evaporation of precip from ZM schemea
  public convtran                 ! convective transport
  public momtran                  ! convective momentum transport

!
! Private data
!
   real(r8) rl         ! wg latent heat of vaporization.
   real(r8) cpres      ! specific heat at constant pressure in j/kg-degk.
   real(r8), parameter :: capelmt = 70._r8  ! threshold value for cape for deep convection.
   real(r8) :: ke           ! Tunable evaporation efficiency set from namelist input zmconv_ke
   real(r8) :: ke_lnd
   real(r8) :: c0_lnd       ! set from namelist input zmconv_c0_lnd
   real(r8) :: c0_ocn       ! set from namelist input zmconv_c0_ocn
   logical  :: zm_org
   real(r8) tau   ! convective time scale
   real(r8),parameter :: c1 = 6.112_r8
   real(r8),parameter :: c2 = 17.67_r8
   real(r8),parameter :: c3 = 243.5_r8
   real(r8) :: tfreez
   real(r8) :: eps1
      

   logical :: no_deep_pbl ! default = .false.
                          ! no_deep_pbl = .true. eliminates deep convection entirely within PBL 
   

!moved from moistconvection.F90
   real(r8) :: rgrav       ! reciprocal of grav
   real(r8) :: rgas        ! gas constant for dry air
   real(r8) :: grav        ! = gravit
   real(r8) :: cp          ! = cpres = cpair
   
   integer  limcnv       ! top interface level limit for convection

   real(r8),parameter ::  tiedke_add = 0.5_r8   

contains


subroutine zm_convi(limcnv_in, zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_ke_lnd, &
                    zmconv_org, no_deep_pbl_in)

   use dycore,       only: dycore_is, get_resolution

   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection
   real(r8),intent(in)           :: zmconv_c0_lnd
   real(r8),intent(in)           :: zmconv_c0_ocn
   real(r8),intent(in)           :: zmconv_ke
   real(r8),intent(in)           :: zmconv_ke_lnd
   logical                       :: zmconv_org
   logical, intent(in), optional :: no_deep_pbl_in  ! no_deep_pbl = .true. eliminates ZM convection entirely within PBL 

   ! local variables
   character(len=32)   :: hgrid           ! horizontal grid specifier

   ! Initialization of ZM constants
   limcnv = limcnv_in
   tfreez = tmelt
   eps1   = epsilo
   rl     = latvap
   cpres  = cpair
   rgrav  = 1.0_r8/gravit
   rgas   = rair
   grav   = gravit
   cp     = cpres

   c0_lnd = zmconv_c0_lnd 
   c0_ocn = zmconv_c0_ocn
   ke     = zmconv_ke
   ke_lnd = zmconv_ke_lnd
   zm_org = zmconv_org

   if ( present(no_deep_pbl_in) )  then
      no_deep_pbl = no_deep_pbl_in
   else
      no_deep_pbl = .false.
   endif

   ! tau=4800. were used in canadian climate center. however, in echam3 t42, 
   ! convection is too weak, thus adjusted to 2400.

   hgrid = get_resolution()
   tau = 3600._r8

   call zm_convr_init(rl, cpres, c0_lnd, c0_ocn, zm_org, tau, &
                      tfreez, eps1, no_deep_pbl, rgrav, rgas, grav, &
                      cp, limcnv)

   if ( masterproc ) then
      write(iulog,*) 'tuning parameters zm_convi: tau',tau
      write(iulog,*) 'tuning parameters zm_convi: c0_lnd',c0_lnd, ', c0_ocn', c0_ocn 
      write(iulog,*) 'tuning parameters zm_convi: ke',ke
      write(iulog,*) 'tuning parameters zm_convi: no_deep_pbl',no_deep_pbl
   endif

   if (masterproc) write(iulog,*)'**** ZM: DILUTE Buoyancy Calculation ****'

end subroutine zm_convi



subroutine zm_convr(lchnk   ,ncol    , &
                    t       ,qh      ,prec    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    ,dpp     , &
                    delt    ,mcon    ,cme     ,cape    , &
                    tpert   ,dlf     ,pflx    ,zdu     ,rprd    , &
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    lengath ,ql      ,rliq    ,landfrac,          &
                    org     ,orgt    , org2d,                     &
                    qu      , &
                    qd      ,dz      ,rppe    ,eps0    , &     
                    cu      ,evp     ,tu     ,td       ,jd      , &
                    done    ,lel     ,jlcl   ,qs       ,qsthat  , &
                    hmn     ,hsat   ,hsthat  ,wteu    , &
                    wted    ,wtdu   ,wtmu    ,wtmd    ,wtcu     , &
                    c0mask  ,wtrpd  ,qds     ,wtevp  )
                
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Main driver for zhang-mcfarlane convection scheme 
! 
! Method: 
! performs deep convective adjustment based on mass-flux closure
! algorithm.
! 
! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
!
! This is contributed code not fully standardized by the CAM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version
! of the CAM where it will include a more straightforward formulation
! and will make use of the standard CAM nomenclature
! 
!-----------------------------------------------------------------------
!
! Compatibility adapter arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns

   real(r8), intent(in) :: t(pcols,pver)          ! grid slice of temperature at mid-layer.
   real(r8), intent(in) :: qh(pcols,pver)   ! grid slice of specific humidity.
   real(r8), intent(in) :: pap(pcols,pver)     
   real(r8), intent(in) :: paph(pcols,pver+1)
   real(r8), intent(in) :: dpp(pcols,pver)        ! local sigma half-level thickness (i.e. dshj).
   real(r8), intent(in) :: zm(pcols,pver)
   real(r8), intent(in) :: geos(pcols)
   real(r8), intent(in) :: zi(pcols,pver+1)
   real(r8), intent(in) :: pblh(pcols)
   real(r8), intent(in) :: tpert(pcols)
   real(r8), intent(in) :: landfrac(pcols) ! RBN Landfrac
   real(r8), intent(in) :: delt            ! length of model time-step in seconds

! output arguments
!
   real(r8), intent(out) :: qtnd(pcols,pver)           ! specific humidity tendency (kg/kg/s)
   real(r8), intent(out) :: heat(pcols,pver)           ! heating rate (dry static energy tendency, W/kg)
   real(r8), intent(out) :: mcon(pcols,pverp)
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)
   real(r8), intent(out) :: cape(pcols)        ! w  convective available potential energy.
   real(r8), intent(out) :: zdu(pcols,pver)
   real(r8), intent(out) :: rprd(pcols,pver)     ! rain production rate
! move these vars from local storage to output so that convective
! transports can be done in outside of conv_cam.
   real(r8), intent(out) :: mu(pcols,pver)
   real(r8), intent(out) :: eu(pcols,pver)
   real(r8), intent(out) :: du(pcols,pver)
   real(r8), intent(out) :: md(pcols,pver)
   real(r8), intent(out) :: ed(pcols,pver)
   real(r8), intent(out) :: dp(pcols,pver)       ! wg layer thickness in mbs (between upper/lower interface).
   real(r8), intent(out) :: dsubcld(pcols)       ! wg layer thickness in mbs between lcl and maxi.
   real(r8), intent(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8), intent(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
   real(r8), intent(out) :: prec(pcols)
   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals

   real(r8), pointer, intent(in) :: org(:,:)       ! Only used if zm_org is true
   real(r8), pointer, intent(inout) :: orgt(:,:)   ! Only used if zm_org is true
   real(r8), pointer, intent(inout) :: org2d(:,:)  ! Only used if zm_org is true

   !Variables needed for Water tracers:
   real(r8), intent(out) :: tu(pcols,pver)     !scattered updraft temperature
   real(r8), intent(out) :: td(pcols,pver)     !scattered downdraft temperature
   real(r8), intent(out) :: cu(pcols,pver)     !scattered version of cug 
   real(r8), intent(out) :: evp(pcols,pver)    !scattered version of evp 
   real(r8), intent(out) :: rppe(pcols,pver)   !g: rain production pre-evap
   real(r8), intent(out) :: qsthat(pcols,pver) !g: qst at interfaces
   real(r8), intent(out) :: hmn(pcols,pver)    !g: environmental moist static energy
   real(r8), intent(out) :: hsat(pcols,pver)   !g: environmental saturation moist static energy
   real(r8), intent(out) :: hsthat(pcols,pver) !g: hsat at interfaces
   real(r8), intent(out) :: dz(pcols,pver)     !g: layer thickness [m]
   real(r8), intent(out) :: eps0(pcols)        !g: not sure, comment later...
   real(r8), intent(out) :: wteu(pcols,pver)   !g: eu pre-unit conversion
   real(r8), intent(out) :: wted(pcols,pver)   !g: ed pre-unit conversion
   real(r8), intent(out) :: wtdu(pcols,pver)   !g: du pre-unit conversion
   real(r8), intent(out) :: wtmu(pcols,pver)   !g: mu pre-unit conversion
   real(r8), intent(out) :: wtmd(pcols,pver)   !g: md pre-unit conversion
   real(r8), intent(out) :: wtcu(pcols,pver)   !g: cu pre-unit conversion
   real(r8), intent(out) :: wtevp(pcols,pver)  !g: evp pre-unit conversion
   real(r8), intent(out) :: wtrpd(pcols,pver)  !g: rprd pre-unit conversion
   real(r8), intent(out) :: c0mask(pcols)      !g: auto-conversion rates
   logical, intent(out)  :: done(pcols,pver)   !g: end of updraft calculation loop

   integer, intent(out) :: jt(pcols)
   integer, intent(out) :: maxg(pcols)
   integer, intent(out) :: ideep(pcols)
   integer, intent(out) :: lengath
   real(r8), intent(out) :: ql(pcols,pver)
   integer, intent(out) :: lel(pcols)
   integer, intent(out) :: jlcl(pcols)
   integer, intent(out) :: jd(pcols)
   real(r8), intent(out) :: qd(pcols,pver)
   real(r8), intent(out) :: qu(pcols,pver)
   real(r8), intent(out) :: qs(pcols,pver)
   real(r8), intent(out) :: qds(pcols,pver)
   call t_startf('ap_zm_convr_run')
   call zm_convr_run(lchnk, ncol, pcols, pver, pverp, &
                     t, qh, prec, jctop, jcbot, &
                     pblh, zm, geos, zi, qtnd, &
                     heat, pap, paph, dpp, &
                     delt, mcon, cme, cape, &
                     tpert, dlf, pflx, zdu, rprd, &
                     mu, md, du, eu, ed, &
                     dp, dsubcld, jt, maxg, ideep, &
                     lengath, ql, rliq, landfrac, &
                     org, orgt, org2d, qu, &
                     qd, dz, rppe, eps0, &
                     cu, evp, tu, td, jd, &
                     done, lel, jlcl, qs, qsthat, &
                     hmn, hsat, hsthat, wteu, &
                     wted, wtdu, wtmu, wtmd, wtcu, &
                     c0mask, wtrpd, qds, wtevp)
   call t_stopf('ap_zm_convr_run')

end subroutine zm_convr

!===============================================================================
subroutine zm_conv_evap(ncol,lchnk, &
     t,pmid,pdel,q, &
     landfrac, &
     tend_s, tend_s_snwprd, tend_s_snwevmlt, tend_q, &
     prdprec, cldfrc, deltat, prec, snow, &
     evpstore, substore, ntprprd, ntsnprd, flxprec, flxsnow )

!-----------------------------------------------------------------------
! Compute tendencies due to evaporation of rain from ZM scheme
!--
! Compute the total precipitation and snow fluxes at the surface.
! Add in the latent heat of fusion for snow formation and melt, since it not dealt with
! in the Zhang-MacFarlane parameterization.
! Evaporate some of the precip directly into the environment using a Sundqvist type algorithm
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    integer,intent(in) :: ncol, lchnk                        ! number of columns and chunk index
    real(r8),intent(in), dimension(pcols,pver) :: t          ! temperature (K)
    real(r8),intent(in), dimension(pcols,pver) :: pmid       ! midpoint pressure (Pa) 
    real(r8),intent(in), dimension(pcols,pver) :: pdel       ! layer thickness (Pa)
    real(r8),intent(in), dimension(pcols,pver) :: q          ! water vapor (kg/kg)
    real(r8),intent(in), dimension(pcols) :: landfrac
    real(r8),intent(inout), dimension(pcols,pver) :: tend_s     ! heating rate (J/kg/s)
    real(r8),intent(inout), dimension(pcols,pver) :: tend_q     ! water vapor tendency (kg/kg/s)
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwprd ! Heating rate of snow production
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow
    


    real(r8), intent(in   ) :: prdprec(pcols,pver)! precipitation production (kg/ks/s)
    real(r8), intent(in   ) :: cldfrc(pcols,pver) ! cloud fraction
    real(r8), intent(in   ) :: deltat             ! time step

    real(r8), intent(inout) :: prec(pcols)        ! Convective-scale preciptn rate
    real(r8), intent(out)   :: snow(pcols)        ! Convective-scale snowfall rate
!
! Output arguments
    real(r8),intent(out) :: flxprec(pcols,pverp)   ! Convective-scale flux of precip at interfaces (kg/m2/s)
    real(r8),intent(out) :: flxsnow(pcols,pverp)   ! Convective-scale flux of snow   at interfaces (kg/m2/s)
    real(r8),intent(out) :: ntprprd(pcols,pver)    ! net precip production in layer
    real(r8),intent(out) :: ntsnprd(pcols,pver)    ! net snow production in layer
 
    !Needed for water tracers:   
    real(r8), intent(out) :: evpstore(pcols,pver) !preciptation evaporation
    real(r8), intent(out) :: substore(pcols,pver) !snow sublimation

    call t_startf('ap_zm_conv_evap_run')
    call zm_conv_evap_run(ncol, lchnk, pcols, pver, pverp, &
                          ke, ke_lnd, zm_org, &
                          t, pmid, pdel, q, landfrac, &
                          tend_s, tend_s_snwprd, tend_s_snwevmlt, tend_q, &
                          prdprec, cldfrc, deltat, prec, snow, &
                          evpstore, substore, ntprprd, ntsnprd, &
                          flxprec, flxsnow)
    call t_stopf('ap_zm_conv_evap_run')

  end subroutine zm_conv_evap



subroutine convtran(lchnk   , &
                    doconvtran,q       ,ncnst   ,mu      ,md      , &
                    du      ,eu      ,ed      ,dp      ,dsubcld , &
                    jt      ,mx      ,ideep   ,il1g    ,il2g    , &
                    nstep   ,fracis  ,dqdt    ,dpdry, Rwt )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective transport of trace species
!
! Mixing ratios may be with respect to either dry or moist air
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: P. Rasch
! 
! Rwt added by J. Nusbaumer to fix mystery variable passing error 
!
!-----------------------------------------------------------------------
   use water_tracer_vars, only: wtrc_ntype
   use water_types, only: iwtice

   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   logical, intent(in) :: doconvtran(ncnst)     ! flag for doing convective transport
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Tracer array including moisture
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: dsubcld(pcols)       ! Delta pressure from cloud base to sfc
   real(r8), intent(in) :: fracis(pcols,pver,ncnst) ! fraction of tracer that is insoluble

   integer, intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: mx(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index

   real(r8), intent(in) :: dpdry(pcols,pver)       ! Delta pressure between interfaces


! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array
  
!water tracers:
   !NOTE:  the 2 at the end is for liquid and ice - JN
   real(r8), intent(out) :: Rwt(pcols,pver,wtrc_ntype(iwtice),2) !water tracer ratio

   call t_startf('ap_zm_conv_convtran_run')
   call zm_conv_convtran_run(lchnk, pcols, pver, size(Rwt, 3), size(Rwt, 4), &
                             doconvtran, q, ncnst, mu, md, &
                             du, eu, ed, dp, dsubcld, &
                             jt, mx, ideep, il1g, il2g, &
                             nstep, fracis, dqdt, dpdry, Rwt)
   call t_stopf('ap_zm_conv_convtran_run')

end subroutine convtran

!=========================================================================================

subroutine momtran(lchnk, ncol, &
                    domomtran,q       ,ncnst   ,mu      ,md    , &
                    du      ,eu      ,ed      ,dp      ,dsubcld , &
                    jt      ,mx      ,ideep   ,il1g    ,il2g    , &
                    nstep   ,dqdt    ,pguall     ,pgdall, icwu, icwd, dt, seten    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective transport of momentum
!
! Mixing ratios may be with respect to either dry or moist air
! 
! Method: 
! Based on the convtran subroutine by P. Rasch
! <Also include any applicable external references.> 
! 
! Author: J. Richter and P. Rasch
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   logical, intent(in) :: domomtran(ncnst)      ! flag for doing convective transport
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Wind array
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: dsubcld(pcols)       ! Delta pressure from cloud base to sfc
   real(r8), intent(in)    :: dt                       !  time step in seconds : 2*delta_t

   integer, intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: mx(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index



! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

   real(r8),intent(out) ::  pguall(pcols,pver,ncnst)      ! Apparent force from  updraft PG
   real(r8),intent(out) ::  pgdall(pcols,pver,ncnst)      ! Apparent force from  downdraft PG

   real(r8),intent(out) ::  icwu(pcols,pver,ncnst)      ! In-cloud winds in updraft
   real(r8),intent(out) ::  icwd(pcols,pver,ncnst)      ! In-cloud winds in downdraft

   real(r8),intent(out) ::  seten(pcols,pver) ! Dry static energy tendency
   
   call t_startf('ap_zm_conv_momtran_run')
   call zm_conv_momtran_run(lchnk, ncol, pcols, pver, pverp, &
                            domomtran, q, ncnst, mu, md, &
                            du, eu, ed, dp, dsubcld, &
                            jt, mx, ideep, il1g, il2g, &
                            nstep, dqdt, pguall, pgdall, &
                            icwu, icwd, dt, seten)
   call t_stopf('ap_zm_conv_momtran_run')

end subroutine momtran

!=========================================================================================

end module zm_conv
