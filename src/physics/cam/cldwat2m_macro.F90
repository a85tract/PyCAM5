
  module cldwat2m_macro

  !--------------------------------------------------- !
  ! Purpose     : CAM Interface for Cloud Macrophysics !
  ! Author      : Sungsu Park                          !
  ! Description : Park et al. 2010.                    !
  ! For questions, contact Sungsu Park                 !
  !                        e-mail : sungsup@ucar.edu   !
  !                        phone  : 303-497-1375       !
  !--------------------------------------------------- !

   use shr_kind_mod,     only: r8=>shr_kind_r8
   use spmd_utils,       only: masterproc
   use ppgrid,           only: pcols, pver, pverp
   use cam_abortutils,   only: endrun
   use physconst,        only: cpair, latvap, latice, rh2o, gravit, rair
   use wv_saturation,    only: qsat_water, svp_water, svp_ice, qsat_ice
   use cam_history,      only: addfld, phys_decomp, outfld, hist_fld_active
   use cam_logfile,      only: iulog
   use ref_pres,         only: top_lev=>trop_cloud_top_lev
   use cldfrc2m,         only: astG_PDF_single, astG_PDF, astG_RHU_single, &
                               astG_RHU, aist_single, aist_vector,         &
                               rhmini_const, rhmaxi=>rhmaxi_const
   use iso_c_binding,    only: c_int64_t, c_loc, c_ptr

   implicit none
   private
   save

   public ::           &
      ini_macro,       &
      mmacro_pcond

   ! -------------- !
   ! Set Parameters !
   ! -------------- !

   ! ------------------------------------------------------------------------------- !
   ! Parameter used for selecting generalized critical RH for liquid and ice stratus !
   ! ------------------------------------------------------------------------------- !

   integer :: i_rhminl ! This is for liquid stratus fraction.
                       ! If 0 : Original fixed critical RH from the namelist.
                       ! If 1 : Add convective detrainment effect on the above '0' option. 
                       !        In this case, 'tau_detw' [s] should be specified below.
                       ! If 2 : Use fully scale-adaptive method.
                       !        In this case, 'tau_detw' [s] and 'c_aniso' [no unit] should
                       !        be specified below. 

   integer :: i_rhmini ! This is for ice stratus fraction.
                       ! If 0 : Original fixed critical RH from the namelist.
                       ! If 1 : Add convective detrainment effect on the above '0' option. 
                       !        In this case, 'tau_deti' [s] should be specified below.
                       ! If 2 : Use fully scale-adaptive method.
                       !        In this case, 'tau_deti' [s] and 'c_aniso' [no unit] should
                       !        be specified below. 
                       ! Note that 'micro_mg_cam' is using below 'rhmini_const', regardless
                       ! of 'i_rhmini'.  This connection should be built in future.

   real(r8), parameter :: tau_detw =100._r8   ! Dissipation time scale of convective liquid condensate detrained
                                              !  into the clear portion. [hr]. 0.5-3 hr is possible.
   real(r8), parameter :: tau_deti =  1._r8   ! Dissipation time scale of convective ice    condensate detrained
                                              !  into the clear portion. [hr]. 0.5-3 hr is possible.
   real(r8), parameter :: c_aniso  =  1._r8   ! Inverse of anisotropic factor of PBL turbulence

   ! ----------------------------- !
   ! Parameters for Liquid Stratus !
   ! ----------------------------- !

   logical,  parameter  :: CAMstfrac    = .false.    ! If .true. (.false.),
                                                     ! use Slingo (triangular PDF-based) liquid stratus fraction
   real(r8), parameter  :: qlst_min     = 2.e-5_r8   ! Minimum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter  :: qlst_max     = 3.e-3_r8   ! Maximum in-stratus LWC constraint [ kg/kg ]
   real(r8), parameter  :: cc           = 0.1_r8     ! For newly formed/dissipated in-stratus CWC ( 0 <= cc <= 1 )
   integer,  parameter  :: niter        = 2          ! For iterative computation of QQ with 'ramda' below.
   real(r8), parameter  :: ramda        = 0.5_r8     ! Explicit : ramda = 0, Implicit : ramda = 1 ( 0<= ramda <= 1 )
   real(r8), private    :: rhminl_const              ! Critical RH for low-level  liquid stratus clouds
   real(r8), private    :: rhminl_adj_land_const     ! rhminl adjustment for snowfree land
   real(r8), private    :: rhminh_const              ! Critical RH for high-level liquid stratus clouds
   real(r8), private    :: premit                    ! Top    height for mid-level liquid stratus fraction
   real(r8), private    :: premib                    ! Bottom height for mid-level liquid stratus fraction

   real(r8), parameter :: qsmall = 1.e-18_r8         ! Smallest mixing ratio considered in the macrophysics
   logical :: use_native_positive_moisture_impl = .false.
   logical :: positive_moisture_impl_selected = .false.
   logical :: positive_moisture_entered_logged = .false.
   logical :: use_native_rhcrit_const_impl = .false.
   logical :: rhcrit_const_impl_selected = .false.
   logical :: rhcrit_const_entered_logged = .false.
   logical :: use_native_instratus_tendency_impl = .false.
   logical :: instratus_tendency_impl_selected = .false.
   logical :: instratus_tendency_entered_logged = .false.
   logical :: use_native_dropnum_limit_impl = .false.
   logical :: dropnum_limit_impl_selected = .false.
   logical :: dropnum_limit_entered_logged = .false.
   logical :: use_native_final_tendency_impl = .false.
   logical :: final_tendency_impl_selected = .false.
   logical :: final_tendency_entered_logged = .false.
   logical :: use_native_iter_state_impl = .false.
   logical :: iter_state_impl_selected = .false.
   logical :: iter_state_entered_logged = .false.
   logical :: use_native_advective_state_impl = .false.
   logical :: advective_state_impl_selected = .false.
   logical :: advective_state_entered_logged = .false.
   logical :: use_native_ref_state_impl = .false.
   logical :: ref_state_impl_selected = .false.
   logical :: ref_state_entered_logged = .false.
   logical :: use_native_linear_state_impl = .false.
   logical :: linear_state_impl_selected = .false.
   logical :: linear_state_entered_logged = .false.
   logical :: use_native_qq_limiter_impl = .false.
   logical :: qq_limiter_impl_selected = .false.
   logical :: qq_limiter_entered_logged = .false.

   contains

   ! -------------- !
   ! Initialization !
   ! -------------- !

   subroutine ini_macro(rhminl_opt_in, rhmini_opt_in)

   !--------------------------------------------------------------------- !
   !                                                                      ! 
   ! Purpose: Initialize constants for the liquid stratiform macrophysics !
   !                                                                      !
   ! Author:  Sungsu Park, Dec.01.2009.                                   !
   !                                                                      !
   !--------------------------------------------------------------------- !

   use cloud_fraction, only: cldfrc_getparams
   use cam_history,    only: addfld, phys_decomp

   integer,  intent(in) :: rhminl_opt_in
   integer,  intent(in) :: rhmini_opt_in

   i_rhminl   = rhminl_opt_in
   i_rhmini   = rhmini_opt_in

   call cldfrc_getparams(rhminl_out=rhminl_const, rhminl_adj_land_out=rhminl_adj_land_const,  &
                         rhminh_out=rhminh_const, premit_out=premit, premib_out=premib)

   if( masterproc ) then
       write(iulog,*) 'Park Macrophysics Parameters'
       write(iulog,*) '  rhminl          = ', rhminl_const
       write(iulog,*) '  rhminl_adj_land = ', rhminl_adj_land_const
       write(iulog,*) '  rhminh          = ', rhminh_const
       write(iulog,*) '  premit          = ', premit
       write(iulog,*) '  premib          = ', premib
       write(iulog,*) '  i_rhminl        = ', i_rhminl
       write(iulog,*) '  i_rhmini        = ', i_rhmini
   end if


   call addfld ('RHMIN_LIQ',     'fraction', pver, 'A', 'Default critical RH for liquid-stratus', phys_decomp)
   call addfld ('RHMIN_ICE',     'fraction', pver, 'A', 'Default critical RH for    ice-stratus', phys_decomp)
   call addfld ('DRHMINPBL_LIQ', 'fraction', pver, 'A', 'Drop of liquid-stratus critical RH by PBL turbulence', phys_decomp)
   call addfld ('DRHMINPBL_ICE', 'fraction', pver, 'A', 'Drop of    ice-stratus critical RH by PBL turbulence', phys_decomp)
   call addfld ('DRHMINDET_LIQ', 'fraction', pver, 'A', 'Drop of liquid-stratus critical RH by convective detrainment', phys_decomp)
   call addfld ('DRHMINDET_ICE', 'fraction', pver, 'A', 'Drop of    ice-stratus critical RH by convective detrainment', phys_decomp)

   end subroutine ini_macro

   ! ------------------------------ !
   ! Stratiform Liquid Macrophysics !
   ! ------------------------------ !

   ! In the version, 'macro --> micro --> advective forcing --> macro...'
   ! A_...: only 'advective forcing' without 'microphysical tendency'
   ! C_...: only 'microphysical tendency'
   ! D_...: only 'detrainment of cumulus condensate'  
   ! So, 'A' and 'C' are exclusive. 

   subroutine mmacro_pcond( lchnk      , ncol       , dt         , p            , dp         ,              &
                            T0         , qv0        , ql0        , qi0          , nl0        , ni0        , &
                            A_T        , A_qv       , A_ql       , A_qi         , A_nl       , A_ni       , &
                            C_T        , C_qv       , C_ql       , C_qi         , C_nl       , C_ni       , C_qlst, &
                            D_T        , D_qv       , D_ql       , D_qi         , D_nl       , D_ni       , &
                            a_cud      , a_cu0      , clrw_old   , clri_old     , landfrac   , snowh      , & 
                            tke        , qtl_flx    , qti_flx    , cmfr_det     , qlr_det    , qir_det    , &
                            s_tendout  , qv_tendout , ql_tendout , qi_tendout   , nl_tendout , ni_tendout , &
                            qme        , qvadj      , qladj      , qiadj        , qllim      , qilim      , &
                            cld        , al_st_star , ai_st_star , ql_st_star   , qi_st_star , do_cldice  )

   use constituents,     only : qmin, cnst_get_ind
   use wv_saturation,    only : findsp_vc

   integer   icol
   integer,  intent(in)    :: lchnk                        ! Chunk number
   integer,  intent(in)    :: ncol                         ! Number of active columns

   ! Input-Output variables

   real(r8), intent(inout) :: T0(pcols,pver)               ! Temperature [K]
   real(r8), intent(inout) :: qv0(pcols,pver)              ! Grid-mean water vapor specific humidity [kg/kg]
   real(r8), intent(inout) :: ql0(pcols,pver)              ! Grid-mean liquid water content [kg/kg]
   real(r8), intent(inout) :: qi0(pcols,pver)              ! Grid-mean ice water content [kg/kg]
   real(r8), intent(inout) :: nl0(pcols,pver)              ! Grid-mean number concentration of cloud liquid droplet [#/kg]
   real(r8), intent(inout) :: ni0(pcols,pver)              ! Grid-mean number concentration of cloud ice    droplet [#/kg]

   ! Input variables

   real(r8), intent(in)    :: dt                           ! Model integration time step [s]
   real(r8), intent(in)    :: p(pcols,pver)                ! Pressure at the layer mid-point [Pa]
   real(r8), intent(in)    :: dp(pcols,pver)               ! Pressure thickness [Pa] > 0

   real(r8), intent(in)    :: A_T(pcols,pver)              ! Non-microphysical advective external forcing of T  [K/s]
   real(r8), intent(in)    :: A_qv(pcols,pver)             ! Non-microphysical advective external forcing of qv [kg/kg/s]
   real(r8), intent(in)    :: A_ql(pcols,pver)             ! Non-microphysical advective external forcing of ql [kg/kg/s]
   real(r8), intent(in)    :: A_qi(pcols,pver)             ! Non-microphysical advective external forcing of qi [kg/kg/s]
   real(r8), intent(in)    :: A_nl(pcols,pver)             ! Non-microphysical advective external forcing of nl [#/kg/s]
   real(r8), intent(in)    :: A_ni(pcols,pver)             ! Non-microphysical advective external forcing of ni [#/kg/s] 

   real(r8), intent(in)    :: C_T(pcols,pver)              ! Microphysical advective external forcing of T  [K/s]
   real(r8), intent(in)    :: C_qv(pcols,pver)             ! Microphysical advective external forcing of qv [kg/kg/s]
   real(r8), intent(in)    :: C_ql(pcols,pver)             ! Microphysical advective external forcing of ql [kg/kg/s]
   real(r8), intent(in)    :: C_qi(pcols,pver)             ! Microphysical advective external forcing of qi [kg/kg/s]
   real(r8), intent(in)    :: C_nl(pcols,pver)             ! Microphysical advective external forcing of nl [#/kg/s]
   real(r8), intent(in)    :: C_ni(pcols,pver)             ! Microphysical advective external forcing of ni [#/kg/s] 
   real(r8), intent(in)    :: C_qlst(pcols,pver)           ! Microphysical advective external forcing of ql
                                                           ! within liquid stratus [kg/kg/s]

   real(r8), intent(in)    :: D_T(pcols,pver)              ! Cumulus detrainment external forcing of T  [K/s]
   real(r8), intent(in)    :: D_qv(pcols,pver)             ! Cumulus detrainment external forcing of qv [kg/kg/s]
   real(r8), intent(in)    :: D_ql(pcols,pver)             ! Cumulus detrainment external forcing of ql [kg/kg/s]
   real(r8), intent(in)    :: D_qi(pcols,pver)             ! Cumulus detrainment external forcing of qi [kg/kg/s]
   real(r8), intent(in)    :: D_nl(pcols,pver)             ! Cumulus detrainment external forcing of nl [#/kg/s]
   real(r8), intent(in)    :: D_ni(pcols,pver)             ! Cumulus detrainment external forcing of qi [#/kg/s] 

   real(r8), intent(in)    :: a_cud(pcols,pver)            ! Old cumulus fraction before update
   real(r8), intent(in)    :: a_cu0(pcols,pver)            ! New cumulus fraction after update

   real(r8), intent(in)    :: clrw_old(pcols,pver)         ! Clear sky fraction at the previous time step for liquid stratus process
   real(r8), intent(in)    :: clri_old(pcols,pver)         ! Clear sky fraction at the previous time step for    ice stratus process
   real(r8), pointer       :: tke(:,:)                     ! (pcols,pverp) TKE from the PBL scheme
   real(r8), pointer       :: qtl_flx(:,:)                 ! (pcols,pverp) overbar(w'qtl') from PBL scheme where qtl = qv + ql
   real(r8), pointer       :: qti_flx(:,:)                 ! (pcols,pverp) overbar(w'qti') from PBL scheme where qti = qv + qi
   real(r8), pointer       :: cmfr_det(:,:)                ! (pcols,pver)  Detrained mass flux from the convection scheme
   real(r8), pointer       :: qlr_det(:,:)                 ! (pcols,pver)  Detrained        ql from the convection scheme
   real(r8), pointer       :: qir_det(:,:)                 ! (pcols,pver)  Detrained        qi from the convection scheme

   real(r8), intent(in)    :: landfrac(pcols)              ! Land fraction
   real(r8), intent(in)    :: snowh(pcols)                 ! Snow depth (liquid water equivalent)
   logical,  intent(in)    :: do_cldice                    ! Whether or not cldice should be prognosed

   ! Output variables

   real(r8), intent(out)   :: s_tendout(pcols,pver)        ! Net tendency of grid-mean s  from 'Micro+Macro' processes [J/kg/s]
   real(r8), intent(out)   :: qv_tendout(pcols,pver)       ! Net tendency of grid-mean qv from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: ql_tendout(pcols,pver)       ! Net tendency of grid-mean ql from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: qi_tendout(pcols,pver)       ! Net tendency of grid-mean qi from 'Micro+Macro' processes [kg/kg/s]
   real(r8), intent(out)   :: nl_tendout(pcols,pver)       ! Net tendency of grid-mean nl from 'Micro+Macro' processes [#/kg/s]
   real(r8), intent(out)   :: ni_tendout(pcols,pver)       ! Net tendency of grid-mean ni from 'Micro+Macro' processes [#/kg/s]

   real(r8), intent(out)   :: qme  (pcols,pver)            ! Net condensation rate [kg/kg/s]
   real(r8), intent(out)   :: qvadj(pcols,pver)            ! adjustment tendency from "positive_moisture" call (vapor)
   real(r8), intent(out)   :: qladj(pcols,pver)            ! adjustment tendency from "positive_moisture" call (liquid)
   real(r8), intent(out)   :: qiadj(pcols,pver)            ! adjustment tendency from "positive_moisture" call (ice)
   real(r8), intent(out)   :: qllim(pcols,pver)            ! tendency from "instratus_condensate" call (liquid)
   real(r8), intent(out)   :: qilim(pcols,pver)            ! tendency from "instratus_condensate" call (ice)

   real(r8), intent(out)   :: cld(pcols,pver)              ! Net cloud fraction ( 0 <= cld <= 1 )
   real(r8), intent(out)   :: al_st_star(pcols,pver)       ! Physical liquid stratus fraction
   real(r8), intent(out)   :: ai_st_star(pcols,pver)       ! Physical ice stratus fraction
   real(r8), intent(out)   :: ql_st_star(pcols,pver)       ! In-stratus LWC [kg/kg] 
   real(r8), intent(out)   :: qi_st_star(pcols,pver)       ! In-stratus IWC [kg/kg] 

   ! --------------- !
   ! Local variables !
   ! --------------- !
   integer :: ixcldliq, ixcldice
 
   integer :: i, j, k, iter, ii, jj                        ! Loop indexes

   ! Thermodynamic state variables

   real(r8) T(pcols,pver)                                  ! Temperature of equilibrium reference state
                                                           ! from which 'Micro & Macro' are computed [K]
   real(r8) T1(pcols,pver)                                 ! Temperature after 'fice_force' on T01  
   real(r8) T_0(pcols,pver)                                ! Temperature after 'instratus_condensate' on T1
   real(r8) T_05(pcols,pver)                               ! Temperature after 'advection' on T_0 
   real(r8) T_prime0(pcols,pver)                           ! Temperature after 'Macrophysics (QQ)' on T_05star
   real(r8) T_dprime(pcols,pver)                           ! Temperature after 'fice_force' on T_prime
   real(r8) T_star(pcols,pver)                             ! Temperature after 'instratus_condensate' on T_dprime

   real(r8) qv(pcols,pver)                                 ! Grid-mean qv of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) qv1(pcols,pver)                                ! Grid-mean qv after 'fice_force' on qv01  
   real(r8) qv_0(pcols,pver)                               ! Grid-mean qv after 'instratus_condensate' on qv1
   real(r8) qv_05(pcols,pver)                              ! Grid-mean qv after 'advection' on qv_0 
   real(r8) qv_prime0(pcols,pver)                          ! Grid-mean qv after 'Macrophysics (QQ)' on qv_05star
   real(r8) qv_dprime(pcols,pver)                          ! Grid-mean qv after 'fice_force' on qv_prime
   real(r8) qv_star(pcols,pver)                            ! Grid-mean qv after 'instratus_condensate' on qv_dprime

   real(r8) ql(pcols,pver)                                 ! Grid-mean ql of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) ql1(pcols,pver)                                ! Grid-mean ql after 'fice_force' on ql01  
   real(r8) ql_0(pcols,pver)                               ! Grid-mean ql after 'instratus_condensate' on ql1
   real(r8) ql_05(pcols,pver)                              ! Grid-mean ql after 'advection' on ql_0 
   real(r8) ql_prime0(pcols,pver)                          ! Grid-mean ql after 'Macrophysics (QQ)' on ql_05star
   real(r8) ql_dprime(pcols,pver)                          ! Grid-mean ql after 'fice_force' on ql_prime
   real(r8) ql_star(pcols,pver)                            ! Grid-mean ql after 'instratus_condensate' on ql_dprime

   real(r8) qi(pcols,pver)                                 ! Grid-mean qi of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) qi1(pcols,pver)                                ! Grid-mean qi after 'fice_force' on qi01  
   real(r8) qi_0(pcols,pver)                               ! Grid-mean qi after 'instratus_condensate' on qi1
   real(r8) qi_05(pcols,pver)                              ! Grid-mean qi after 'advection' on qi_0 
   real(r8) qi_prime0(pcols,pver)                          ! Grid-mean qi after 'Macrophysics (QQ)' on qi_05star
   real(r8) qi_dprime(pcols,pver)                          ! Grid-mean qi after 'fice_force' on qi_prime
   real(r8) qi_star(pcols,pver)                            ! Grid-mean qi after 'instratus_condensate' on qi_dprime

   real(r8) nl(pcols,pver)                                 ! Grid-mean nl of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) nl1(pcols,pver)                                ! Grid-mean nl after 'fice_force' on nl01  
   real(r8) nl_0(pcols,pver)                               ! Grid-mean nl after 'instratus_condensate' on nl1
   real(r8) nl_05(pcols,pver)                              ! Grid-mean nl after 'advection' on nl_0 
   real(r8) nl_prime0(pcols,pver)                          ! Grid-mean nl after 'Macrophysics (QQ)' on nl_05star
   real(r8) nl_dprime(pcols,pver)                          ! Grid-mean nl after 'fice_force' on nl_prime
   real(r8) nl_star(pcols,pver)                            ! Grid-mean nl after 'instratus_condensate' on nl_dprime

   real(r8) ni(pcols,pver)                                 ! Grid-mean ni of equilibrium reference state from which
                                                           ! 'Micro & Macro' are computed [kg/kg]
   real(r8) ni1(pcols,pver)                                ! Grid-mean ni after 'fice_force' on ni01  
   real(r8) ni_0(pcols,pver)                               ! Grid-mean ni after 'instratus_condensate' on ni1
   real(r8) ni_05(pcols,pver)                              ! Grid-mean ni after 'advection' on ni_0 
   real(r8) ni_prime0(pcols,pver)                          ! Grid-mean ni after 'Macrophysics (QQ)' on ni_05star
   real(r8) ni_dprime(pcols,pver)                          ! Grid-mean ni after 'fice_force' on ni_prime
   real(r8) ni_star(pcols,pver)                            ! Grid-mean ni after 'instratus_condensate' on ni_dprime

   real(r8) a_st(pcols,pver)                               ! Stratus fraction of equilibrium reference state 
   real(r8) a_st_0(pcols,pver)                             ! Stratus fraction at '_0' state
   real(r8) a_st_star(pcols,pver)                          ! Stratus fraction at '_star' state

   real(r8) al_st(pcols,pver)                              ! Liquid stratus fraction of equilibrium reference state 
   real(r8) al_st_0(pcols,pver)                            ! Liquid stratus fraction at '_0' state
   real(r8) al_st_nc(pcols,pver)                           ! Non-physical liquid stratus fraction in the non-cumulus pixels

   real(r8) ai_st(pcols,pver)                              ! Ice stratus fraction of equilibrium reference state 
   real(r8) ai_st_0(pcols,pver)                            ! Ice stratus fraction at '_0' state
   real(r8) ai_st_nc(pcols,pver)                           ! Non-physical ice stratus fraction in the non-cumulus pixels

   real(r8) ql_st(pcols,pver)                              ! In-stratus LWC of equilibrium reference state [kg/kg] 
   real(r8) ql_st_0(pcols,pver)                            ! In-stratus LWC at '_0' state

   real(r8) qi_st(pcols,pver)                              ! In-stratus IWC of equilibrium reference state [kg/kg] 
   real(r8) qi_st_0(pcols,pver)                            ! In-stratus IWC at '_0' state

 ! Cumulus properties 

   real(r8) dacudt(pcols,pver)
   real(r8) a_cu(pcols,pver)

 ! Adjustment tendency in association with 'positive_moisture'

   real(r8) Tten_pwi1(pcols,pver)                          ! Pre-process T  tendency of input equilibrium state [K/s] 
   real(r8) qvten_pwi1(pcols,pver)                         ! Pre-process qv tendency of input equilibrium state [kg/kg/s]
   real(r8) qlten_pwi1(pcols,pver)                         ! Pre-process ql tendency of input equilibrium state [kg/kg/s]
   real(r8) qiten_pwi1(pcols,pver)                         ! Pre-process qi tendency of input equilibrium state [kg/kg/s]
   real(r8) nlten_pwi1(pcols,pver)                         ! Pre-process nl tendency of input equilibrium state [#/kg/s]
   real(r8) niten_pwi1(pcols,pver)                         ! Pre-process ni tendency of input equilibrium state [#/kg/s] 

   real(r8) Tten_pwi2(pcols,pver)                          ! Post-process T  tendency of provisional equilibrium state [K/s] 
   real(r8) qvten_pwi2(pcols,pver)                         ! Post-process qv tendency of provisional equilibrium state [kg/kg/s]
   real(r8) qlten_pwi2(pcols,pver)                         ! Post-process ql tendency of provisional equilibrium state [kg/kg/s]
   real(r8) qiten_pwi2(pcols,pver)                         ! Post-process qi tendency of provisional equilibrium state [kg/kg/s]
   real(r8) nlten_pwi2(pcols,pver)                         ! Post-process nl tendency of provisoonal equilibrium state [#/kg/s]
   real(r8) niten_pwi2(pcols,pver)                         ! Post-process ni tendency of provisional equilibrium state [#/kg/s] 

   real(r8) A_T_adj(pcols,pver)                            ! After applying external advective forcing [K/s]
   real(r8) A_qv_adj(pcols,pver)                           ! After applying external advective forcing [kg/kg/s]
   real(r8) A_ql_adj(pcols,pver)                           ! After applying external advective forcing [kg/kg/s]
   real(r8) A_qi_adj(pcols,pver)                           ! After applying external advective forcing [kg/kg/s]
   real(r8) A_nl_adj(pcols,pver)                           ! After applying external advective forcing [#/kg/s]
   real(r8) A_ni_adj(pcols,pver)                           ! After applying external advective forcing [#/kg/s]

 ! Adjustment tendency in association with 'instratus_condensate'

   real(r8) QQw1(pcols,pver)           ! Effective adjustive condensation into water due to 'instratus_condensate' [kg/kg/s]
   real(r8) QQi1(pcols,pver)           ! Effective adjustive condensation into ice   due to 'instratus_condensate' [kg/kg/s]
   real(r8) QQw2(pcols,pver)           ! Effective adjustive condensation into water due to 'instratus_condensate' [kg/kg/s]
   real(r8) QQi2(pcols,pver)           ! Effective adjustive condensation into ice   due to 'instratus_condensate' [kg/kg/s]

   real(r8) QQnl1(pcols,pver)          ! Tendency of nl associated with QQw1 only when QQw1<0 (net evaporation) [#/kg/s]
   real(r8) QQni1(pcols,pver)          ! Tendency of ni associated with QQi1 only when QQw1<0 (net evaporation) [#/kg/s]
   real(r8) QQnl2(pcols,pver)          ! Tendency of nl associated with QQw2 only when QQw2<0 (net evaporation) [#/kg/s]
   real(r8) QQni2(pcols,pver)          ! Tendency of ni associated with QQi2 only when QQw2<0 (net evaporation) [#/kg/s]

 ! Macrophysical process tendency variables

   real(r8) QQ(pcols,pver)             ! Net condensation rate into water+ice           [kg/kg/s] 
   real(r8) QQw(pcols,pver)            ! Net condensation rate into water               [kg/kg/s] 
   real(r8) QQi(pcols,pver)            ! Net condensation rate into ice                 [kg/kg/s]
   real(r8) QQnl(pcols,pver)           ! Tendency of nl associated with QQw both for condensation and evaporation [#/kg/s]
   real(r8) QQni(pcols,pver)           ! Tendency of ni associated with QQi both for condensation and evaporation [#/kg/s]
   real(r8) ACnl(pcols,pver)           ! Cloud liquid droplet (nl) activation tendency [#/kg/s]
   real(r8) ACni(pcols,pver)           ! Cloud ice    droplet (ni) activation tendency [#/kg/s]

   real(r8) QQw_prev(pcols,pver)   
   real(r8) QQi_prev(pcols,pver)   
   real(r8) QQnl_prev(pcols,pver)  
   real(r8) QQni_prev(pcols,pver)  

   real(r8) QQw_prog(pcols,pver)   
   real(r8) QQi_prog(pcols,pver)   
   real(r8) QQnl_prog(pcols,pver)  
   real(r8) QQni_prog(pcols,pver)  

   real(r8) QQ_final(pcols,pver)                           
   real(r8) QQw_final(pcols,pver)                           
   real(r8) QQi_final(pcols,pver)                           
   real(r8) QQn_final(pcols,pver)                           
   real(r8) QQnl_final(pcols,pver)                          
   real(r8) QQni_final(pcols,pver)                          

   real(r8) QQ_all(pcols,pver)         ! QQw_all    + QQi_all
   real(r8) QQw_all(pcols,pver)        ! QQw_final  + QQw1  + QQw2  + qlten_pwi1 + qlten_pwi2 + A_ql_adj [kg/kg/s]
   real(r8) QQi_all(pcols,pver)        ! QQi_final  + QQi1  + QQi2  + qiten_pwi1 + qiten_pwi2 + A_qi_adj [kg/kg/s]
   real(r8) QQn_all(pcols,pver)        ! QQnl_all   + QQni_all
   real(r8) QQnl_all(pcols,pver)       ! QQnl_final + QQnl1 + QQnl2 + nlten_pwi1 + nlten_pwi2 + ACnl [#/kg/s]
   real(r8) QQni_all(pcols,pver)       ! QQni_final + QQni1 + QQni2 + niten_pwi1 + niten_pwi2 + ACni [#/kg/s]

 ! Coefficient for computing QQ and related processes

   real(r8) U(pcols,pver)                                  ! Grid-mean RH
   real(r8) U_nc(pcols,pver)                               ! Mean RH of non-cumulus pixels
   real(r8) G_nc(pcols,pver)                               ! d(U_nc)/d(a_st_nc)
   real(r8) F_nc(pcols,pver)                               ! A function of second parameter for a_st_nc
   real(r8) alpha                                          ! = 1/qs
   real(r8) beta                                           ! = (qv/qs**2)*dqsdT
   real(r8) betast                                         ! = alpha*dqsdT
   real(r8) gammal                                         ! = alpha + (latvap/cpair)*beta
   real(r8) gammai                                         ! = alpha + ((latvap+latice)/cpair)*beta
   real(r8) gammaQ                                         ! = alpha + (latvap/cpair)*beta
   real(r8) deltal                                         ! = 1 + a_st*(latvap/cpair)*(betast/alpha)
   real(r8) deltai                                         ! = 1 + a_st*((latvap+latice)/cpair)*(betast/alpha)
   real(r8) A_Tc                                           ! Advective external forcing of Tc [K/s]
   real(r8) A_qt                                           ! Advective external forcing of qt [kg/kg/s]
   real(r8) C_Tc                                           ! Microphysical forcing of Tc [K/s]
   real(r8) C_qt                                           ! Microphysical forcing of qt [kg/kg/s]
   real(r8) dTcdt                                          ! d(Tc)/dt      [K/s]
   real(r8) dqtdt                                          ! d(qt)/dt      [kg/kg/s]
   real(r8) dqtstldt                                       ! d(qt_alst)/dt [kg/kg/s]
   real(r8) dqidt                                          ! d(qi)/dt      [kg/kg/s]

   real(r8) dqlstdt                                        ! d(ql_st)/dt [kg/kg/s]
   real(r8) dalstdt                                        ! d(al_st)/dt  [1/s]
   real(r8) dastdt                                         ! d(a_st)/dt  [1/s]

   real(r8) anic                                           ! Fractional area of non-cumulus and non-ice stratus fraction
   real(r8) GG                                             ! G_nc(i,k)/anic

   real(r8) aa(2,2)
   real(r8) bb(2,1)

   real(r8) zeros(pcols,pver)

   real(r8) qmin1(pcols,pver)
   real(r8) qmin2(pcols,pver)
   real(r8) qmin3(pcols,pver)

   real(r8) esat_a(pcols)                             ! Saturation water vapor pressure [Pa]
   real(r8) qsat_a(pcols,pver)                        ! Saturation water vapor specific humidity [kg/kg]
   real(r8) Twb_aw(pcols)                             ! Wet-bulb temperature [K]
   real(r8) qvwb_aw(pcols,pver)                       ! Wet-bulb water vapor specific humidity [kg/kg]

   real(r8) esat_b(pcols)                                 
   real(r8) qsat_b(pcols)                                 
   real(r8) dqsdT_b(pcols)                                 

   logical  land
   real(r8) tmp

   real(r8) d_rhmin_liq_PBL(pcols,pver)
   real(r8) d_rhmin_ice_PBL(pcols,pver)
   real(r8) d_rhmin_liq_det(pcols,pver)
   real(r8) d_rhmin_ice_det(pcols,pver)
   real(r8) rhmini_arr(pcols,pver)
   real(r8) rhminl_arr(pcols,pver)
   real(r8) rhminl_adj_land_arr(pcols,pver)
   real(r8) rhminh_arr(pcols,pver) 
   real(r8) rhmin_liq_diag(pcols,pver)
   real(r8) rhmin_ice_diag(pcols,pver)

   real(r8) QQmax,QQmin,QQwmin,QQimin                      ! For limiting QQ
   real(r8) cone                                           ! Number close to but smaller than 1

   cone            = 0.999_r8
   zeros(:ncol,:)  = 0._r8

   ! ------------------------------------ !
   ! Global initialization of main output !
   ! ------------------------------------ !

     s_tendout(:ncol,:)     = 0._r8
     qv_tendout(:ncol,:)    = 0._r8
     ql_tendout(:ncol,:)    = 0._r8
     qi_tendout(:ncol,:)    = 0._r8
     nl_tendout(:ncol,:)    = 0._r8
     ni_tendout(:ncol,:)    = 0._r8

     qme(:ncol,:)           = 0._r8

     cld(:ncol,:)           = 0._r8
     al_st_star(:ncol,:)    = 0._r8
     ai_st_star(:ncol,:)    = 0._r8
     ql_st_star(:ncol,:)    = 0._r8
     qi_st_star(:ncol,:)    = 0._r8

   ! --------------------------------------- !
   ! Initialization of internal 2D variables !
   ! --------------------------------------- !

     T(:ncol,:)             = 0._r8
     T1(:ncol,:)            = 0._r8
     T_0(:ncol,:)           = 0._r8
     T_05(:ncol,:)          = 0._r8
     T_prime0(:ncol,:)      = 0._r8
     T_dprime(:ncol,:)      = 0._r8
     T_star(:ncol,:)        = 0._r8

     qv(:ncol,:)            = 0._r8
     qv1(:ncol,:)           = 0._r8
     qv_0(:ncol,:)          = 0._r8
     qv_05(:ncol,:)         = 0._r8
     qv_prime0(:ncol,:)     = 0._r8
     qv_dprime(:ncol,:)     = 0._r8
     qv_star(:ncol,:)       = 0._r8

     ql(:ncol,:)            = 0._r8
     ql1(:ncol,:)           = 0._r8
     ql_0(:ncol,:)          = 0._r8
     ql_05(:ncol,:)         = 0._r8
     ql_prime0(:ncol,:)     = 0._r8
     ql_dprime(:ncol,:)     = 0._r8
     ql_star(:ncol,:)       = 0._r8

     qi(:ncol,:)            = 0._r8
     qi1(:ncol,:)           = 0._r8
     qi_0(:ncol,:)          = 0._r8
     qi_05(:ncol,:)         = 0._r8
     qi_prime0(:ncol,:)     = 0._r8
     qi_dprime(:ncol,:)     = 0._r8
     qi_star(:ncol,:)       = 0._r8

     nl(:ncol,:)            = 0._r8
     nl1(:ncol,:)           = 0._r8
     nl_0(:ncol,:)          = 0._r8
     nl_05(:ncol,:)         = 0._r8
     nl_prime0(:ncol,:)     = 0._r8
     nl_dprime(:ncol,:)     = 0._r8
     nl_star(:ncol,:)       = 0._r8

     ni(:ncol,:)            = 0._r8
     ni1(:ncol,:)           = 0._r8
     ni_0(:ncol,:)          = 0._r8
     ni_05(:ncol,:)         = 0._r8
     ni_prime0(:ncol,:)     = 0._r8
     ni_dprime(:ncol,:)     = 0._r8
     ni_star(:ncol,:)       = 0._r8

     a_st(:ncol,:)          = 0._r8
     a_st_0(:ncol,:)        = 0._r8
     a_st_star(:ncol,:)     = 0._r8

     al_st(:ncol,:)         = 0._r8
     al_st_0(:ncol,:)       = 0._r8
     al_st_nc(:ncol,:)      = 0._r8

     ai_st(:ncol,:)         = 0._r8
     ai_st_0(:ncol,:)       = 0._r8
     ai_st_nc(:ncol,:)      = 0._r8

     ql_st(:ncol,:)         = 0._r8
     ql_st_0(:ncol,:)       = 0._r8

     qi_st(:ncol,:)         = 0._r8
     qi_st_0(:ncol,:)       = 0._r8

 ! Cumulus properties 

     dacudt(:ncol,:)        = 0._r8
     a_cu(:ncol,:)          = 0._r8

 ! Adjustment tendency in association with 'positive_moisture'

     Tten_pwi1(:ncol,:)     = 0._r8
     qvten_pwi1(:ncol,:)    = 0._r8
     qlten_pwi1(:ncol,:)    = 0._r8
     qiten_pwi1(:ncol,:)    = 0._r8
     nlten_pwi1(:ncol,:)    = 0._r8
     niten_pwi1(:ncol,:)    = 0._r8

     Tten_pwi2(:ncol,:)     = 0._r8
     qvten_pwi2(:ncol,:)    = 0._r8
     qlten_pwi2(:ncol,:)    = 0._r8
     qiten_pwi2(:ncol,:)    = 0._r8
     nlten_pwi2(:ncol,:)    = 0._r8
     niten_pwi2(:ncol,:)    = 0._r8

     A_T_adj(:ncol,:)       = 0._r8
     A_qv_adj(:ncol,:)      = 0._r8
     A_ql_adj(:ncol,:)      = 0._r8
     A_qi_adj(:ncol,:)      = 0._r8
     A_nl_adj(:ncol,:)      = 0._r8
     A_ni_adj(:ncol,:)      = 0._r8

     qvadj   (:ncol,:)      = 0._r8
     qladj   (:ncol,:)      = 0._r8
     qiadj   (:ncol,:)      = 0._r8

 ! Adjustment tendency in association with 'instratus_condensate'

     QQw1(:ncol,:)          = 0._r8
     QQi1(:ncol,:)          = 0._r8
     QQw2(:ncol,:)          = 0._r8
     QQi2(:ncol,:)          = 0._r8

     QQnl1(:ncol,:)         = 0._r8
     QQni1(:ncol,:)         = 0._r8
     QQnl2(:ncol,:)         = 0._r8
     QQni2(:ncol,:)         = 0._r8

     QQnl(:ncol,:)          = 0._r8
     QQni(:ncol,:)          = 0._r8

 ! Macrophysical process tendency variables

     QQ(:ncol,:)            = 0._r8
     QQw(:ncol,:)           = 0._r8
     QQi(:ncol,:)           = 0._r8
     QQnl(:ncol,:)          = 0._r8
     QQni(:ncol,:)          = 0._r8
     ACnl(:ncol,:)          = 0._r8
     ACni(:ncol,:)          = 0._r8

     QQw_prev(:ncol,:)      = 0._r8
     QQi_prev(:ncol,:)      = 0._r8
     QQnl_prev(:ncol,:)     = 0._r8
     QQni_prev(:ncol,:)     = 0._r8

     QQw_prog(:ncol,:)      = 0._r8
     QQi_prog(:ncol,:)      = 0._r8
     QQnl_prog(:ncol,:)     = 0._r8
     QQni_prog(:ncol,:)     = 0._r8

     QQ_final(:ncol,:)      = 0._r8                        
     QQw_final(:ncol,:)     = 0._r8                  
     QQi_final(:ncol,:)     = 0._r8           
     QQn_final(:ncol,:)     = 0._r8    
     QQnl_final(:ncol,:)    = 0._r8
     QQni_final(:ncol,:)    = 0._r8

     QQ_all(:ncol,:)        = 0._r8
     QQw_all(:ncol,:)       = 0._r8
     QQi_all(:ncol,:)       = 0._r8
     QQn_all(:ncol,:)       = 0._r8
     QQnl_all(:ncol,:)      = 0._r8
     QQni_all(:ncol,:)      = 0._r8

 ! Coefficient for computing QQ and related processes

     U(:ncol,:)             = 0._r8
     U_nc(:ncol,:)          = 0._r8
     G_nc(:ncol,:)          = 0._r8
     F_nc(:ncol,:)          = 0._r8

 ! Other

     qmin1(:ncol,:)         = 0._r8
     qmin2(:ncol,:)         = 0._r8
     qmin3(:ncol,:)         = 0._r8

   ! ---------------- !
   ! Main computation ! 
   ! ---------------- !

   ! Compute critical RH for stratus
   call rhcrit_calc( &
      ncol, dp, T0, p, &
      clrw_old, clri_old, tke, qtl_flx, &
      qti_flx, cmfr_det, qlr_det, qir_det, &
      rhmini_arr, rhminl_arr, rhminl_adj_land_arr, rhminh_arr, &
      d_rhmin_liq_PBL, d_rhmin_ice_PBL, d_rhmin_liq_det, d_rhmin_ice_det)

   ! -------------------------------------------------------------------- !
   ! Compute cumulus-related properties and prepare input reference state. !
   ! -------------------------------------------------------------------- !

   call input_state_codon_wrap(ncol, dt, a_cu0, a_cud, T0, qv0, ql0, qi0, nl0, ni0, &
        dacudt, T1, qv1, ql1, qi1, nl1, ni1)
   
   ! ---------------------------------------------------------------------- !
   ! Check if input non-cumulus pixels satisfie a non-negative constraint.  !
   ! If not, force all water vapor substances to be positive in all layers. !
   ! We should use 'old' cumulus properties for this routine.               !                
   ! ---------------------------------------------------------------------- !
   
   call cnst_get_ind( 'CLDLIQ', ixcldliq )
   call cnst_get_ind( 'CLDICE', ixcldice )


   call qmin_fill_codon_wrap(ncol, qmin(1), qmin(ixcldliq), qmin(ixcldice), qmin1, qmin2, qmin3)

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp, & 
                           qv1, ql1, qi1, T1, qvten_pwi1, qlten_pwi1, &
                           qiten_pwi1, Tten_pwi1, do_cldice)

   call dropnum_limit_codon_wrap(1, ncol, dt, ql1, qi1, nl1, ni1, nlten_pwi1, niten_pwi1)

   ! ------------------------------------------------------------- !
   ! Impose 'in-stratus condensate amount constraint'              !
   ! such that it is bounded by two limiting values.               !      
   ! This should also use 'old' cumulus properties since it is     !
   ! before applying external forcings.                            ! 
   ! Below 'QQw1,QQi1' are effective adjustive condensation        ! 
   ! Although this process also involves freezing of cloud         !
   ! liquid into ice, they can be and only can be expressed        !
   ! in terms of effective condensation.                           !
   ! ------------------------------------------------------------- !

   do k = top_lev, pver
      call instratus_condensate( lchnk, ncol, k,                                   &
                                 p(:,k), T1(:,k), qv1(:,k), ql1(:,k), qi1(:,k),    &
                                 ni1(:,k),                                         &
                                 a_cud(:,k), zeros(:,k), zeros(:,k),               &
                                 zeros(:,k), zeros(:,k), zeros(:,k),               &
                                 landfrac, snowh,                                  &
                                 rhmini_arr(:,k), rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k), &
                                 T_0(:,k), qv_0(:,k), ql_0(:,k), qi_0(:,k),        & 
                                 al_st_0(:,k), ai_st_0(:,k), ql_st_0(:,k), qi_st_0(:,k) )
      call instratus_tendency_codon_wrap(1, ncol, dt, cone, ql_0(:,k), qi_0(:,k), ql1(:,k), qi1(:,k), &
                                         nl1(:,k), ni1(:,k), al_st_0(:,k), ai_st_0(:,k), a_st_0(:,k), &
                                         QQw1(:,k), QQi1(:,k), QQnl1(:,k), QQni1(:,k), nl_0(:,k), ni_0(:,k))
   enddo

   ! ----------------------------------------------------------------------------- !
   ! Check if non-cumulus pixels of '_05' state satisfies non-negative constraint. !
   ! If not, force all water substances of '_05' state to be positive by imposing  !
   ! adjustive advection. We should use 'new' cumulus properties for this routine. !                
   ! ----------------------------------------------------------------------------- !

   call advective_state_codon_wrap(ncol, dt, T_0, qv_0, ql_0, qi_0, nl_0, ni_0, &
        A_T, C_T, A_qv, C_qv, A_ql, C_ql, A_qi, C_qi, A_nl, C_nl, A_ni, C_ni, &
        T_05, qv_05, ql_05, qi_05, nl_05, ni_05)

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp, & 
                           qv_05, ql_05, qi_05, T_05, A_qv_adj, &
                           A_ql_adj, A_qi_adj, A_T_adj, do_cldice)

   ! -------------------------------------------------------------- !
   ! Define reference state at the first iteration. This will be    !
   ! continuously updated within the iteration loop below.          !
   ! While equlibrium state properties are already output from the  !
   ! 'instratus_condensate', they will be re-computed within the    !
   ! each iteration process. At the first iteration, they will      !
   ! produce exactly identical results. Note that except at the     !
   ! very first iteration iter = 1, we must use updated cumulus     !
   ! properties at all the other iteration processes. Even at the   !
   ! first iteration, we should use updated cumulus properties      !
   ! when computing limiters for (Q,P,E).                           !
   ! -------------------------------------------------------------- !

   ! -------------------------------------------------------------- !
   ! Define variables at the reference state of the first iteration !
   ! -------------------------------------------------------------- !

   call ref_state_codon_wrap(ncol, T_0, qv_0, ql_0, qi_0, al_st_0, ai_st_0, a_st_0, &
        ql_st_0, qi_st_0, nl_0, ni_0, T, qv, ql, qi, al_st, ai_st, a_st, ql_st, qi_st, nl, ni)

   ! -------------------------- !
   ! Main iterative computation !
   ! -------------------------- !

   do k = top_lev, pver
      call findsp_vc(qv_05(:ncol,k), T_05(:ncol,k), p(:ncol,k), .false., &
                     Twb_aw(:ncol), qvwb_aw(:ncol,k))
      call qsat_water(T_05(1:ncol,k), p(1:ncol,k), &
                      esat_a(1:ncol), qsat_a(1:ncol,k))
   enddo

   do iter = 1, niter

      ! ------------------------------------------ !
      ! Initialize array within the iteration loop !
      ! ------------------------------------------ !

      QQ(:,:)         = 0._r8
      QQw(:,:)        = 0._r8
      QQi(:,:)        = 0._r8
      QQnl(:,:)       = 0._r8
      QQni(:,:)       = 0._r8 
      QQw2(:,:)       = 0._r8
      QQi2(:,:)       = 0._r8
      QQnl2(:,:)      = 0._r8
      QQni2(:,:)      = 0._r8
      nlten_pwi2(:,:) = 0._r8
      niten_pwi2(:,:) = 0._r8
      ACnl(:,:)       = 0._r8
      ACni(:,:)       = 0._r8 
      aa(:,:)         = 0._r8
      bb(:,:)         = 0._r8

      do k = top_lev, pver

      call qsat_water(T(1:ncol,k), p(1:ncol,k), &
                      esat_b(1:ncol), qsat_b(1:ncol), dqsdt=dqsdT_b(1:ncol))

      if( iter .eq. 1 ) then
          a_cu(:ncol,k) = a_cud(:ncol,k)
      else
          a_cu(:ncol,k) = a_cu0(:ncol,k)
      endif
      do i = 1, ncol
         U(i,k)    =  qv(i,k)/qsat_b(i)
         U_nc(i,k) =  U(i,k)
      enddo
      if( CAMstfrac ) then
          call astG_RHU(U_nc(:,k),p(:,k),qv(:,k),landfrac(:),snowh(:),al_st_nc(:,k),G_nc(:,k),ncol,&
                        rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k))                          
      else
          call astG_PDF(U_nc(:,k),p(:,k),qv(:,k),landfrac(:),snowh(:),al_st_nc(:,k),G_nc(:,k),ncol,&
                        rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k))
      endif
      call aist_vector(qv(:,k),T(:,k),p(:,k),qi(:,k),ni(:,k),landfrac(:),snowh(:),ai_st_nc(:,k),ncol,&
                       rhmaxi, rhmini_arr(:,k), rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k))

      ai_st(:ncol,k)  =  (1._r8-a_cu(:ncol,k))*ai_st_nc(:ncol,k)
      al_st(:ncol,k)  =  (1._r8-a_cu(:ncol,k))*al_st_nc(:ncol,k)
      a_st(:ncol,k)   =  max(al_st(:ncol,k),ai_st(:ncol,k))  

      do i = 1, ncol

         ! -------------------------------------------------------- !
         ! Compute basic thermodynamic coefficients for computing Q !
         ! -------------------------------------------------------- !

         alpha  =  1._r8/qsat_b(i)
         beta   =  dqsdT_b(i)*(qv(i,k)/qsat_b(i)**2)
         betast =  alpha*dqsdT_b(i) 
         gammal =  alpha + (latvap/cpair)*beta
         gammai =  alpha + ((latvap+latice)/cpair)*beta
         gammaQ =  alpha + (latvap/cpair)*beta
         deltal =  1._r8 + a_st(i,k)*(latvap/cpair)*(betast/alpha)
         deltai =  1._r8 + a_st(i,k)*((latvap+latice)/cpair)*(betast/alpha)
         A_Tc   =  A_T(i,k)+A_T_adj(i,k)-(latvap/cpair)*(A_ql(i,k)+A_ql_adj(i,k))-((latvap+latice)/cpair)*(A_qi(i,k)+A_qi_adj(i,k))
         A_qt   =  A_qv(i,k) + A_qv_adj(i,k) + A_ql(i,k) + A_ql_adj(i,k) + A_qi(i,k) + A_qi_adj(i,k)
         C_Tc   =  C_T(i,k) - (latvap/cpair)*C_ql(i,k) - ((latvap+latice)/cpair)*C_qi(i,k)
         C_qt   =  C_qv(i,k) + C_ql(i,k) + C_qi(i,k)
         dTcdt  =  A_Tc + C_Tc
         dqtdt  =  A_qt + C_qt
       ! dqtstldt = A_qt + C_ql(i,k)/max(1.e-2_r8,al_st(i,k))                             ! Original  
       ! dqtstldt = A_qt - A_qi(i,k) - A_qi_adj(i,k) + C_ql(i,k)/max(1.e-2_r8,al_st(i,k)) ! New 1 on Dec.30.2009.
         dqtstldt = A_qt - A_qi(i,k) - A_qi_adj(i,k) + C_qlst(i,k)                        ! New 2 on Dec.30.2009.
       ! dqtstldt = A_qt + C_qt                                                           ! Original Conservative treatment
       ! dqtstldt = A_qt - A_qi(i,k) - A_qi_adj(i,k) + C_qt - C_qi(i,k)            ! New Conservative treatment on Dec.30.2009
         dqidt = A_qi(i,k) + A_qi_adj(i,k) + C_qi(i,k) 

         anic    = max(1.e-8_r8,(1._r8-a_cu(i,k)))
         GG      = G_nc(i,k)/anic
         aa(1,1) = gammal*al_st(i,k)
         aa(1,2) = GG + gammal*cc*ql_st(i,k)          
         aa(2,1) = alpha + (latvap/cpair)*betast*al_st(i,k)
         aa(2,2) = (latvap/cpair)*betast*cc*ql_st(i,k) 
         bb(1,1) = alpha*dqtdt - beta*dTcdt - gammai*dqidt - GG*al_st_nc(i,k)*dacudt(i,k) + F_nc(i,k) 
         bb(2,1) = alpha*dqtstldt - betast*(dTcdt + ((latvap+latice)/cpair)*dqidt) 
         call gaussj(aa(1:2,1:2),2,2,bb(1:2,1),1,1)
         dqlstdt = bb(1,1)
         dalstdt = bb(2,1)
         QQ(i,k) = al_st(i,k)*dqlstdt + cc*ql_st(i,k)*dalstdt - ( A_ql(i,k) + A_ql_adj(i,k) + C_ql(i,k) )

	      enddo
	      enddo

      call qq_limiter_codon_wrap(ncol, dt, cone, qmin(1), qv_05, ql_05, qi_05, nl_05, ni_05, &
           qsat_a, qvwb_aw, QQ, QQw, QQi, QQnl, QQni)

	    ! -------------------------------------------------------------------- !
    ! Until now, we have finished computing all necessary tendencies       ! 
    ! from the equilibrium input state (T_0).                              !
    ! If ramda = 0 : fully explicit scheme                                 !
    !    ramda = 1 : fully implicit scheme                                 !
    ! Note that 'ramda = 0.5 with niter = 2' can mimic                     !
    ! -------------------------------------------------------------------- !

      call iter_state_codon_wrap(iter, ncol, dt, &
           QQw, QQi, QQnl, QQni, QQw_prev, QQi_prev, QQnl_prev, QQni_prev, &
           QQw_prog, QQi_prog, QQnl_prog, QQni_prog, &
           T_0, qv_0, ql_0, qi_0, nl_0, ni_0, &
           A_T, A_T_adj, C_T, A_qv, A_qv_adj, C_qv, A_ql, A_ql_adj, C_ql, &
           A_qi, A_qi_adj, C_qi, A_nl, C_nl, A_ni, C_ni, &
           T_prime0, qv_prime0, ql_prime0, qi_prime0, nl_prime0, ni_prime0)

   ! -------------------------------------------------- !
   ! Perform diagnostic 'positive_moisture' constraint. !
   ! -------------------------------------------------- !

   call detrain_state_codon_wrap(ncol, dt, T_prime0, qv_prime0, ql_prime0, qi_prime0, nl_prime0, ni_prime0, &
        D_T, D_qv, D_ql, D_qi, D_nl, D_ni, T_dprime, qv_dprime, ql_dprime, qi_dprime, nl_dprime, ni_dprime)

   call positive_moisture( ncol, dt, qmin1, qmin2, qmin3, dp,          & 
                           qv_dprime, ql_dprime, qi_dprime, T_dprime,  &
                           qvten_pwi2, qlten_pwi2, qiten_pwi2, Tten_pwi2, do_cldice)

   call dropnum_limit_codon_wrap(2, ncol, dt, ql_dprime, qi_dprime, nl_dprime, ni_dprime, &
                                 nlten_pwi2, niten_pwi2)

   ! -------------------------------------------------------------- !
   ! Add tendency associated with detrainment of cumulus condensate !
   ! This tendency is not used in computing Q                       !
   ! Since D_ql,D_qi,D_nl,D_ni > 0, don't need to worry about       !
   ! negative scalar.                                               !
   ! This tendency is not reflected into Fzs2, which is OK.         !
   ! -------------------------------------------------------------- !

   ! ---------------------------------------------------------- !
   ! Impose diagnostic upper and lower limits on the in-stratus !
   ! condensate amount. This produces a final equilibrium state !
   ! at the end of each iterative process.                      !
   ! ---------------------------------------------------------- !

   do k = top_lev, pver
      call instratus_condensate( lchnk          , ncol           , k              , p(:,k)        , &
                                 T_dprime(:,k)  , qv_dprime(:,k) , ql_dprime(:,k) , qi_dprime(:,k), &
                                 ni_dprime(:,k) ,                                                   &
                                 a_cu0(:,k)     , zeros(:,k)     , zeros(:,k)     ,                 & 
                                 zeros(:,k)     , zeros(:,k)     , zeros(:,k)     ,                 &
                                 landfrac       , snowh          ,                                  &
                                 rhmini_arr(:,k), rhminl_arr(:,k), rhminl_adj_land_arr(:,k), rhminh_arr(:,k), &
                                 T_star(:,k)    , qv_star(:,k)   , ql_star(:,k)   , qi_star(:,k)  , & 
                                 al_st_star(:,k), ai_st_star(:,k), ql_st_star(:,k), qi_st_star(:,k) )
      call instratus_tendency_codon_wrap(2, ncol, dt, cone, ql_star(:,k), qi_star(:,k), ql_dprime(:,k), &
                                         qi_dprime(:,k), nl_dprime(:,k), ni_dprime(:,k), al_st_star(:,k), &
                                         ai_st_star(:,k), a_st_star(:,k), QQw2(:,k), QQi2(:,k), QQnl2(:,k), &
                                         QQni2(:,k), nl_star(:,k), ni_star(:,k))
   enddo

   ! ------------------------------------------ !
   ! Final adjustment of droplet concentration. !
   ! Set # to zero if there is no cloud.        !
   ! ------------------------------------------ !

   call dropnum_limit_codon_wrap(3, ncol, dt, ql_star, qi_star, nl_star, ni_star, ACnl, ACni)

   ! ----------------------------------------------------- !
   ! Define equilibrium reference state for next iteration !
   ! ----------------------------------------------------- !

   call ref_state_codon_wrap(ncol, T_star, qv_star, ql_star, qi_star, al_st_star, ai_st_star, a_st_star, &
        ql_st_star, qi_st_star, nl_star, ni_star, T, qv, ql, qi, al_st, ai_st, a_st, ql_st, qi_st, nl, ni)

   enddo ! End of 'iter' prognostic iterative computation

   ! ------------------------------------------------------------------------ !
   ! Compute final tendencies of main output variables and diagnostic outputs !
   ! Note that the very input state [T0,qv0,ql0,qi0] are                      !
   ! marched to [T_star,qv_star,ql_star,qi_star] with equilibrium             !
   ! stratus informations of [a_st_star,ql_st_star,qi_st_star] by             !
   ! below final tendencies and [A_T,A_qv,A_ql,A_qi]                          !
   ! ------------------------------------------------------------------------ !

   call final_tendency_codon_wrap(ncol, dt, do_cldice, &
        QQw_prog, QQi_prog, QQnl_prog, QQni_prog, QQw1, QQi1, QQw2, QQi2, &
        qlten_pwi1, qlten_pwi2, qiten_pwi1, qiten_pwi2, A_ql_adj, A_qi_adj, &
        QQnl1, QQni1, QQnl2, QQni2, nlten_pwi1, nlten_pwi2, niten_pwi1, niten_pwi2, &
        ACnl, ACni, A_nl_adj, A_ni_adj, qvten_pwi1, qvten_pwi2, A_qv_adj, &
        T_star, qv_star, ql_star, qi_star, nl_star, ni_star, &
        A_T, C_T, A_qv, C_qv, A_ql, C_ql, A_qi, C_qi, A_nl, C_nl, A_ni, C_ni, &
        a_st_star, a_cu0, QQw_final, QQi_final, QQ_final, QQw_all, QQi_all, QQ_all, &
        QQnl_final, QQni_final, QQn_final, QQnl_all, QQni_all, QQn_all, &
        qme, qvadj, qladj, qiadj, qllim, qilim, &
        s_tendout, qv_tendout, ql_tendout, qi_tendout, nl_tendout, ni_tendout, cld, &
        T0, qv0, ql0, qi0, nl0, ni0)

   if (hist_fld_active('RHMIN_LIQ')) then
      ! Compute default critical RH as a function of height and surface type as in the current code.
      rhmin_liq_diag(:,:) = 0._r8
      do k = top_lev, pver
         do i = 1, ncol
            land = nint(landfrac(i)) == 1
            if( p(i,k) .ge. premib ) then
               if( land .and. (snowh(i).le.0.000001_r8) ) then
                  rhmin_liq_diag(i,k) = rhminl_const - rhminl_adj_land_const
               else
                  rhmin_liq_diag(i,k) = rhminl_const
               endif
            elseif( p(i,k) .lt. premit ) then
               rhmin_liq_diag(i,k) = rhminh_const
            else
               tmp = (premib-(max(p(i,k),premit)))/(premib-premit)
               rhmin_liq_diag(i,k) = rhminh_const*tmp + rhminl_const*(1.0_r8-tmp)
            endif
         end do
      end do
      call outfld( 'RHMIN_LIQ',      rhmin_liq_diag,  pcols, lchnk )
   end if

   rhmin_ice_diag(:,:) = rhminh_const
   call outfld( 'RHMIN_ICE',      rhmin_ice_diag,  pcols, lchnk )

   call outfld( 'DRHMINPBL_LIQ', d_rhmin_liq_PBL,  pcols, lchnk )
   call outfld( 'DRHMINPBL_ICE', d_rhmin_ice_PBL,  pcols, lchnk )
   call outfld( 'DRHMINDET_LIQ', d_rhmin_liq_det,  pcols, lchnk )
   call outfld( 'DRHMINDET_ICE', d_rhmin_ice_det,  pcols, lchnk )

   end subroutine mmacro_pcond


!=======================================================================================================

subroutine rhcrit_calc( &
   ncol, dp, T0, p, &
   clrw_old, clri_old, tke, qtl_flx, &
   qti_flx, cmfr_det, qlr_det, qir_det, &
   rhmini_arr, rhminl_arr, rhminl_adj_land_arr, rhminh_arr, &
   d_rhmin_liq_PBL, d_rhmin_ice_PBL, d_rhmin_liq_det, d_rhmin_ice_det)

   ! ------------------------------------------------- !
   ! Compute a drop of critical RH for stratus by      !
   ! (1) PBL turbulence, and                           !
   ! (2) convective detrainment.                       !
   ! Note that all of 'd_rhmin...' terms are positive. !
   ! ------------------------------------------------- !

   integer,  intent(in) :: ncol                         ! Number of active columns
   real(r8), intent(in) :: dp(pcols,pver)               ! Pressure thickness [Pa] > 0
   real(r8), intent(in) :: T0(pcols,pver)               ! Temperature [K]
   real(r8), intent(in) :: p(pcols,pver)                ! Pressure at the layer mid-point [Pa]
   real(r8), intent(in) :: clrw_old(pcols,pver)         ! Clear sky fraction at the previous time step for liquid stratus process
   real(r8), intent(in) :: clri_old(pcols,pver)         ! Clear sky fraction at the previous time step for    ice stratus process
   real(r8), pointer    :: tke(:,:)                     ! (pcols,pverp) TKE from the PBL scheme
   real(r8), pointer    :: qtl_flx(:,:)                 ! (pcols,pverp) overbar(w'qtl') from PBL scheme where qtl = qv + ql
   real(r8), pointer    :: qti_flx(:,:)                 ! (pcols,pverp) overbar(w'qti') from PBL scheme where qti = qv + qi
   real(r8), pointer    :: cmfr_det(:,:)                ! (pcols,pver)  Detrained mass flux from the convection scheme
   real(r8), pointer    :: qlr_det(:,:)                 ! (pcols,pver)  Detrained        ql from the convection scheme
   real(r8), pointer    :: qir_det(:,:)                 ! (pcols,pver)  Detrained        qi from the convection scheme

   real(r8), target, intent(out) :: rhmini_arr(pcols,pver)
   real(r8), target, intent(out) :: rhminl_arr(pcols,pver)
   real(r8), target, intent(out) :: rhminl_adj_land_arr(pcols,pver)
   real(r8), target, intent(out) :: rhminh_arr(pcols,pver)
   real(r8), intent(out) :: d_rhmin_liq_PBL(pcols,pver)
   real(r8), intent(out) :: d_rhmin_ice_PBL(pcols,pver)
   real(r8), intent(out) :: d_rhmin_liq_det(pcols,pver)
   real(r8), intent(out) :: d_rhmin_ice_det(pcols,pver)

   ! local variables

   integer :: i, k

   real(r8) :: esat_tmp(pcols)          ! Dummy for saturation vapor pressure calc.
   real(r8) :: qsat_tmp(pcols)          ! Saturation water vapor specific humidity [kg/kg]
   real(r8) :: sig_tmp

   interface
      subroutine cldwat2m_rhcrit_const_codon(pcols_c, pver_c, rhmini_const_c, rhminl_const_c, &
           rhminl_adj_land_const_c, rhminh_const_c, rhmini_p, rhminl_p, rhminl_adj_land_p, rhminh_p) &
           bind(c, name="cldwat2m_rhcrit_const_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: pcols_c, pver_c
         real(c_double), value :: rhmini_const_c, rhminl_const_c, rhminl_adj_land_const_c, rhminh_const_c
         type(c_ptr), value :: rhmini_p, rhminl_p, rhminl_adj_land_p, rhminh_p
      end subroutine cldwat2m_rhcrit_const_codon
   end interface
   !---------------------------------------------------------------------------------------------------

   call rhcrit_const_select_impl()
   if (.not. use_native_rhcrit_const_impl) then
      call rhcrit_const_log_entered()
      call cldwat2m_rhcrit_const_codon(int(pcols, c_int64_t), int(pver, c_int64_t), &
           rhmini_const, rhminl_const, rhminl_adj_land_const, rhminh_const, &
           c_loc(rhmini_arr(1,1)), c_loc(rhminl_arr(1,1)), c_loc(rhminl_adj_land_arr(1,1)), &
           c_loc(rhminh_arr(1,1)))
      return
   end if


   ! ---------------------------------- !
   ! Calc critical RH for ice stratus   !
   ! ---------------------------------- !

   rhmini_arr(:,:) = rhmini_const

   if (i_rhmini > 0) then

      ! Compute the drop of critical RH by convective detrainment of cloud condensate

      do k = top_lev, pver
         do i = 1, ncol
            d_rhmin_ice_det(i,k) = tau_deti*(gravit/dp(i,k))*cmfr_det(i,k)*clri_old(i,k)*qir_det(i,k)*3.6e6_r8 
            d_rhmin_ice_det(i,k) = max(0._r8,min(0.5_r8,d_rhmin_ice_det(i,k)))
         end do
      end do

      if (i_rhmini == 1) then
         rhmini_arr(:ncol,:) = rhmini_const - d_rhmin_ice_det(:ncol,:)
      end if

   end if

   if (i_rhmini == 2) then

      ! Compute the drop of critical RH by the variability induced by PBL turbulence

      do k = top_lev, pver
         call qsat_ice(T0(1:ncol,k), p(1:ncol,k), esat_tmp(1:ncol), qsat_tmp(1:ncol))

         do i = 1, ncol
            sig_tmp = 0.5_r8 * ( qti_flx(i,k)   / sqrt(max(qsmall,tke(i,k))) + & 
                                 qti_flx(i,k+1) / sqrt(max(qsmall,tke(i,k+1))) )
            d_rhmin_ice_PBL(i,k) = c_aniso*sig_tmp/max(qsmall,qsat_tmp(i)) 
            d_rhmin_ice_PBL(i,k) = max(0._r8,min(0.5_r8,d_rhmin_ice_PBL(i,k)))

            rhmini_arr(i,k) = 1._r8 - d_rhmin_ice_PBL(i,k) - d_rhmin_ice_det(i,k)
         end do
      end do
   end if

   if (i_rhmini > 0) then
      do k = top_lev, pver
         do i = 1, ncol
            rhmini_arr(i,k) = max(0._r8,min(rhmaxi,rhmini_arr(i,k))) 
         end do
      end do
   end if

   ! ------------------------------------- !
   ! Choose critical RH for liquid stratus !
   ! ------------------------------------- !

   rhminl_arr(:,:)          = rhminl_const
   rhminl_adj_land_arr(:,:) = rhminl_adj_land_const
   rhminh_arr(:,:)          = rhminh_const

   if (i_rhminl > 0) then

      ! Compute the drop of critical RH by convective detrainment of cloud condensate

      do k = top_lev, pver
         do i = 1, ncol
            d_rhmin_liq_det(i,k) = tau_detw*(gravit/dp(i,k))*cmfr_det(i,k)*clrw_old(i,k)*qlr_det(i,k)*3.6e6_r8 
            d_rhmin_liq_det(i,k) = max(0._r8,min(0.5_r8,d_rhmin_liq_det(i,k)))
         end do
      end do

      if (i_rhminl == 1) then
         rhminl_arr(:ncol,top_lev:) = rhminl_const - d_rhmin_liq_det(:ncol,top_lev:)
         rhminh_arr(:ncol,top_lev:) = rhminh_const - d_rhmin_liq_det(:ncol,top_lev:)
      end if

   end if

   if (i_rhminl == 2) then

      ! Compute the drop of critical RH by the variability induced by PBL turbulence

      do k = top_lev, pver
         call qsat_water(T0(1:ncol,k), p(1:ncol,k), esat_tmp(1:ncol), qsat_tmp(1:ncol))

         do i = 1, ncol
            sig_tmp = 0.5_r8 * ( qtl_flx(i,k)   / sqrt(max(qsmall,tke(i,k))) + & 
                                 qtl_flx(i,k+1) / sqrt(max(qsmall,tke(i,k+1))) )
            d_rhmin_liq_PBL(i,k) = c_aniso*sig_tmp/max(qsmall,qsat_tmp(i)) 
            d_rhmin_liq_PBL(i,k) = max(0._r8,min(0.5_r8,d_rhmin_liq_PBL(i,k)))

            rhminl_arr(i,k) = 1._r8 - d_rhmin_liq_PBL(i,k) - d_rhmin_liq_det(i,k)
            rhminl_adj_land_arr(i,k) = 0._r8
            rhminh_arr(i,k) = rhminl_arr(i,k)
         end do
      end do
   end if

   if (i_rhminl > 0) then
      do k = top_lev, pver
         do i = 1, ncol
            rhminl_arr(i,k) = max(rhminl_adj_land_arr(i,k),min(1._r8,rhminl_arr(i,k))) 
            rhminh_arr(i,k) = max(0._r8,min(1._r8,rhminh_arr(i,k))) 
         end do
      end do
   end if

end subroutine rhcrit_calc

!=======================================================================================================

   subroutine rhcrit_const_select_impl()
   character(len=32) :: impl_name
   integer :: n, status

   if (rhcrit_const_impl_selected) return
   call get_environment_variable('CLDWAT2M_RHCRIT_CONST_IMPL', value=impl_name, length=n, status=status)
   use_native_rhcrit_const_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_rhcrit_const_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_rhcrit_const_impl = .false.
      case default
         use_native_rhcrit_const_impl = .false.
      end select
   end if
   if (i_rhmini /= 0 .or. i_rhminl /= 0) use_native_rhcrit_const_impl = .true.
   rhcrit_const_impl_selected = .true.
   if (masterproc) then
      if (use_native_rhcrit_const_impl) then
         write(iulog,*) 'cldwat2m_rhcrit_const implementation = native'
      else
         write(iulog,*) 'cldwat2m_rhcrit_const implementation = codon'
      end if
   end if
   end subroutine rhcrit_const_select_impl

   subroutine rhcrit_const_log_entered()
   if (rhcrit_const_entered_logged) return
   rhcrit_const_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_rhcrit_const entered (macrophysics constant rhcrit fields = codon)'
   end if
   end subroutine rhcrit_const_log_entered

!=======================================================================================================

   subroutine instratus_tendency_select_impl()

   character(len=32) :: impl_name
   integer :: n, status

   if (instratus_tendency_impl_selected) return
   call get_environment_variable('CLDWAT2M_INSTRATUS_TENDENCY_IMPL', value=impl_name, length=n, status=status)
   use_native_instratus_tendency_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_instratus_tendency_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_instratus_tendency_impl = .false.
      case default
         use_native_instratus_tendency_impl = .false.
      end select
   end if
   instratus_tendency_impl_selected = .true.
   if (masterproc) then
      if (use_native_instratus_tendency_impl) then
         write(iulog,*) 'cldwat2m_instratus_tendency implementation = native'
      else
         write(iulog,*) 'cldwat2m_instratus_tendency implementation = codon'
      end if
   end if
   end subroutine instratus_tendency_select_impl

   subroutine instratus_tendency_log_entered()
   if (instratus_tendency_entered_logged) return
   instratus_tendency_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_instratus_tendency entered (macrophysics instratus condensate tendency = codon)'
   end if
   end subroutine instratus_tendency_log_entered

   subroutine instratus_tendency_codon_wrap(stage, ncol, dt, cone, ql_new, qi_new, ql_old, qi_old, &
                                            nl_old, ni_old, al_st_new, ai_st_new, a_st_new, &
                                            QQw, QQi, QQnl, QQni, nl_new, ni_new)

   integer, intent(in) :: stage
   integer, intent(in) :: ncol
   real(r8), intent(in) :: dt
   real(r8), intent(in) :: cone
   real(r8), target, intent(in) :: ql_new(pcols), qi_new(pcols), ql_old(pcols), qi_old(pcols)
   real(r8), target, intent(in) :: nl_old(pcols), ni_old(pcols), al_st_new(pcols), ai_st_new(pcols)
   real(r8), target, intent(out) :: a_st_new(pcols), QQw(pcols), QQi(pcols), QQnl(pcols), QQni(pcols)
   real(r8), target, intent(out) :: nl_new(pcols), ni_new(pcols)

   integer :: i

   interface
      subroutine cldwat2m_instratus_tendency_codon(stage_c, ncol_c, dt_c, qsmall_c, cone_c, &
           ql_new_p, qi_new_p, ql_old_p, qi_old_p, nl_old_p, ni_old_p, al_st_new_p, ai_st_new_p, &
           a_st_new_p, QQw_p, QQi_p, QQnl_p, QQni_p, nl_new_p, ni_new_p) &
           bind(c, name="cldwat2m_instratus_tendency_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c
         real(c_double), value :: dt_c, qsmall_c, cone_c
         type(c_ptr), value :: ql_new_p, qi_new_p, ql_old_p, qi_old_p, nl_old_p, ni_old_p
         type(c_ptr), value :: al_st_new_p, ai_st_new_p, a_st_new_p, QQw_p, QQi_p, QQnl_p, QQni_p
         type(c_ptr), value :: nl_new_p, ni_new_p
      end subroutine cldwat2m_instratus_tendency_codon
   end interface

   call instratus_tendency_select_impl()
   if (.not. use_native_instratus_tendency_impl) then
      call instratus_tendency_log_entered()
      call cldwat2m_instratus_tendency_codon(int(stage, c_int64_t), int(ncol, c_int64_t), dt, qsmall, cone, &
           c_loc(ql_new(1)), c_loc(qi_new(1)), c_loc(ql_old(1)), c_loc(qi_old(1)), &
           c_loc(nl_old(1)), c_loc(ni_old(1)), c_loc(al_st_new(1)), c_loc(ai_st_new(1)), &
           c_loc(a_st_new(1)), c_loc(QQw(1)), c_loc(QQi(1)), c_loc(QQnl(1)), c_loc(QQni(1)), &
           c_loc(nl_new(1)), c_loc(ni_new(1)))
      return
   endif

   do i = 1, ncol
      a_st_new(i) = max(al_st_new(i), ai_st_new(i))
      QQw(i) = (ql_new(i) - ql_old(i))/dt
      QQi(i) = (qi_new(i) - qi_old(i))/dt
      QQnl(i) = 0._r8
      QQni(i) = 0._r8
      if( QQw(i) .le. 0._r8 ) then
         if( (stage .eq. 2 .and. ql_old(i) .ge. qsmall) .or. &
             (stage .ne. 2 .and. ql_old(i) .gt. qsmall) ) then
            QQnl(i) = QQw(i)*nl_old(i)/ql_old(i)
            QQnl(i) = min(0._r8,cone*max(QQnl(i),-nl_old(i)/dt))
         endif
      endif
      if( QQi(i) .le. 0._r8 ) then
         if( qi_old(i) .gt. qsmall ) then
            QQni(i) = QQi(i)*ni_old(i)/qi_old(i)
            QQni(i) = min(0._r8,cone*max(QQni(i),-ni_old(i)/dt))
         endif
      endif
      nl_new(i) = max(0._r8,nl_old(i)+QQnl(i)*dt)
      ni_new(i) = max(0._r8,ni_old(i)+QQni(i)*dt)
   enddo

   end subroutine instratus_tendency_codon_wrap

!=======================================================================================================

   subroutine dropnum_limit_select_impl()

   character(len=32) :: impl_name
   integer :: n, status

   if (dropnum_limit_impl_selected) return
   call get_environment_variable('CLDWAT2M_DROPNUM_LIMIT_IMPL', value=impl_name, length=n, status=status)
   use_native_dropnum_limit_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_dropnum_limit_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_dropnum_limit_impl = .false.
      case default
         use_native_dropnum_limit_impl = .false.
      end select
   end if
   dropnum_limit_impl_selected = .true.
   if (masterproc) then
      if (use_native_dropnum_limit_impl) then
         write(iulog,*) 'cldwat2m_dropnum_limit implementation = native'
      else
         write(iulog,*) 'cldwat2m_dropnum_limit implementation = codon'
      end if
   end if
   end subroutine dropnum_limit_select_impl

   subroutine dropnum_limit_log_entered()
   if (dropnum_limit_entered_logged) return
   dropnum_limit_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_dropnum_limit entered (macrophysics cloud number limiter = codon)'
   end if
   end subroutine dropnum_limit_log_entered

   subroutine dropnum_limit_codon_wrap(stage, ncol, dt, ql, qi, nl, ni, nlten, niten)

   integer, intent(in) :: stage
   integer, intent(in) :: ncol
   real(r8), intent(in) :: dt
   real(r8), target, intent(in) :: ql(pcols,pver), qi(pcols,pver)
   real(r8), target, intent(inout) :: nl(pcols,pver), ni(pcols,pver)
   real(r8), target, intent(out) :: nlten(pcols,pver), niten(pcols,pver)

   integer :: i, k

   interface
      subroutine cldwat2m_dropnum_limit_codon(stage_c, ncol_c, pcols_c, pver_c, top_lev_c, dt_c, qsmall_c, &
           ql_p, qi_p, nl_p, ni_p, nlten_p, niten_p) bind(c, name="cldwat2m_dropnum_limit_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: dt_c, qsmall_c
         type(c_ptr), value :: ql_p, qi_p, nl_p, ni_p, nlten_p, niten_p
      end subroutine cldwat2m_dropnum_limit_codon
   end interface

   call dropnum_limit_select_impl()
   if (.not. use_native_dropnum_limit_impl) then
      call dropnum_limit_log_entered()
      call cldwat2m_dropnum_limit_codon(int(stage, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(top_lev, c_int64_t), dt, qsmall, c_loc(ql(1,1)), c_loc(qi(1,1)), &
           c_loc(nl(1,1)), c_loc(ni(1,1)), c_loc(nlten(1,1)), c_loc(niten(1,1)))
      return
   endif

   do k = top_lev, pver
      do i = 1, ncol
         nlten(i,k) = 0._r8
         niten(i,k) = 0._r8
         if( ql(i,k) .lt. qsmall ) then
            nlten(i,k) = -nl(i,k)/dt
            nl(i,k) = 0._r8
         endif
         if( qi(i,k) .lt. qsmall ) then
            niten(i,k) = -ni(i,k)/dt
            ni(i,k) = 0._r8
         endif
      enddo
   enddo

   end subroutine dropnum_limit_codon_wrap

!=======================================================================================================

   subroutine ref_state_select_impl()

   character(len=32) :: impl_name
   integer :: n, status

   if (ref_state_impl_selected) return
   call get_environment_variable('CLDWAT2M_REF_STATE_IMPL', value=impl_name, length=n, status=status)
   use_native_ref_state_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_ref_state_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_ref_state_impl = .false.
      case default
         use_native_ref_state_impl = .false.
      end select
   end if
   ref_state_impl_selected = .true.
   if (masterproc) then
      if (use_native_ref_state_impl) then
         write(iulog,*) 'cldwat2m_ref_state implementation = native'
      else
         write(iulog,*) 'cldwat2m_ref_state implementation = codon'
      end if
   end if
   end subroutine ref_state_select_impl

   subroutine ref_state_log_entered()
   if (ref_state_entered_logged) return
   ref_state_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_ref_state entered (macrophysics reference state copy = codon)'
   end if
   end subroutine ref_state_log_entered

   subroutine ref_state_codon_wrap(ncol, T_src, qv_src, ql_src, qi_src, al_st_src, ai_st_src, a_st_src, &
        ql_st_src, qi_st_src, nl_src, ni_src, T_dst, qv_dst, ql_dst, qi_dst, al_st_dst, ai_st_dst, &
        a_st_dst, ql_st_dst, qi_st_dst, nl_dst, ni_dst)

   integer, intent(in) :: ncol
   real(r8), target, intent(in) :: T_src(pcols,pver), qv_src(pcols,pver), ql_src(pcols,pver)
   real(r8), target, intent(in) :: qi_src(pcols,pver), al_st_src(pcols,pver), ai_st_src(pcols,pver)
   real(r8), target, intent(in) :: a_st_src(pcols,pver), ql_st_src(pcols,pver), qi_st_src(pcols,pver)
   real(r8), target, intent(in) :: nl_src(pcols,pver), ni_src(pcols,pver)
   real(r8), target, intent(out) :: T_dst(pcols,pver), qv_dst(pcols,pver), ql_dst(pcols,pver)
   real(r8), target, intent(out) :: qi_dst(pcols,pver), al_st_dst(pcols,pver), ai_st_dst(pcols,pver)
   real(r8), target, intent(out) :: a_st_dst(pcols,pver), ql_st_dst(pcols,pver), qi_st_dst(pcols,pver)
   real(r8), target, intent(out) :: nl_dst(pcols,pver), ni_dst(pcols,pver)

   integer :: i, k

   interface
      subroutine cldwat2m_ref_state_codon(ncol_c, pcols_c, pver_c, top_lev_c, &
           T_src_p, qv_src_p, ql_src_p, qi_src_p, al_st_src_p, ai_st_src_p, a_st_src_p, &
           ql_st_src_p, qi_st_src_p, nl_src_p, ni_src_p, T_dst_p, qv_dst_p, ql_dst_p, qi_dst_p, &
           al_st_dst_p, ai_st_dst_p, a_st_dst_p, ql_st_dst_p, qi_st_dst_p, nl_dst_p, ni_dst_p) &
           bind(c, name="cldwat2m_ref_state_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         type(c_ptr), value :: T_src_p, qv_src_p, ql_src_p, qi_src_p, al_st_src_p, ai_st_src_p, a_st_src_p
         type(c_ptr), value :: ql_st_src_p, qi_st_src_p, nl_src_p, ni_src_p
         type(c_ptr), value :: T_dst_p, qv_dst_p, ql_dst_p, qi_dst_p, al_st_dst_p, ai_st_dst_p
         type(c_ptr), value :: a_st_dst_p, ql_st_dst_p, qi_st_dst_p, nl_dst_p, ni_dst_p
      end subroutine cldwat2m_ref_state_codon
   end interface

   call ref_state_select_impl()
   if (.not. use_native_ref_state_impl) then
      call ref_state_log_entered()
      call cldwat2m_ref_state_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(top_lev, c_int64_t), c_loc(T_src(1,1)), c_loc(qv_src(1,1)), &
           c_loc(ql_src(1,1)), c_loc(qi_src(1,1)), c_loc(al_st_src(1,1)), c_loc(ai_st_src(1,1)), &
           c_loc(a_st_src(1,1)), c_loc(ql_st_src(1,1)), c_loc(qi_st_src(1,1)), c_loc(nl_src(1,1)), &
           c_loc(ni_src(1,1)), c_loc(T_dst(1,1)), c_loc(qv_dst(1,1)), c_loc(ql_dst(1,1)), &
           c_loc(qi_dst(1,1)), c_loc(al_st_dst(1,1)), c_loc(ai_st_dst(1,1)), c_loc(a_st_dst(1,1)), &
           c_loc(ql_st_dst(1,1)), c_loc(qi_st_dst(1,1)), c_loc(nl_dst(1,1)), c_loc(ni_dst(1,1)))
      return
   endif

   do k = top_lev, pver
      do i = 1, ncol
         T_dst(i,k) = T_src(i,k)
         qv_dst(i,k) = qv_src(i,k)
         ql_dst(i,k) = ql_src(i,k)
         qi_dst(i,k) = qi_src(i,k)
         al_st_dst(i,k) = al_st_src(i,k)
         ai_st_dst(i,k) = ai_st_src(i,k)
         a_st_dst(i,k) = a_st_src(i,k)
         ql_st_dst(i,k) = ql_st_src(i,k)
         qi_st_dst(i,k) = qi_st_src(i,k)
         nl_dst(i,k) = nl_src(i,k)
         ni_dst(i,k) = ni_src(i,k)
      enddo
   enddo

   end subroutine ref_state_codon_wrap

!=======================================================================================================

   subroutine linear_state_select_impl()

   character(len=32) :: impl_name
   integer :: n, status

   if (linear_state_impl_selected) return
   call get_environment_variable('CLDWAT2M_LINEAR_STATE_IMPL', value=impl_name, length=n, status=status)
   use_native_linear_state_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_linear_state_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_linear_state_impl = .false.
      case default
         use_native_linear_state_impl = .false.
      end select
   end if
   linear_state_impl_selected = .true.
   if (masterproc) then
      if (use_native_linear_state_impl) then
         write(iulog,*) 'cldwat2m_linear_state implementation = native'
      else
         write(iulog,*) 'cldwat2m_linear_state implementation = codon'
      end if
   end if
   end subroutine linear_state_select_impl

   subroutine linear_state_log_entered()
   if (linear_state_entered_logged) return
   linear_state_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_linear_state entered (macrophysics input/qmin/detrainment state = codon)'
   end if
   end subroutine linear_state_log_entered

   subroutine input_state_codon_wrap(ncol, dt, a_cu0, a_cud, T0, qv0, ql0, qi0, nl0, ni0, &
        dacudt, T1, qv1, ql1, qi1, nl1, ni1)

   integer, intent(in) :: ncol
   real(r8), intent(in) :: dt
   real(r8), target, intent(in) :: a_cu0(pcols,pver), a_cud(pcols,pver), T0(pcols,pver), qv0(pcols,pver)
   real(r8), target, intent(inout) :: ql0(pcols,pver), qi0(pcols,pver), nl0(pcols,pver), ni0(pcols,pver)
   real(r8), target, intent(inout) :: dacudt(pcols,pver)
   real(r8), target, intent(out) :: T1(pcols,pver), qv1(pcols,pver), ql1(pcols,pver)
   real(r8), target, intent(out) :: qi1(pcols,pver), nl1(pcols,pver), ni1(pcols,pver)

   integer :: i, k

   interface
      subroutine cldwat2m_input_state_codon(ncol_c, pcols_c, pver_c, top_lev_c, dt_c, &
           a_cu0_p, a_cud_p, T0_p, qv0_p, ql0_p, qi0_p, nl0_p, ni0_p, dacudt_p, &
           T1_p, qv1_p, ql1_p, qi1_p, nl1_p, ni1_p) bind(c, name="cldwat2m_input_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: dt_c
         type(c_ptr), value :: a_cu0_p, a_cud_p, T0_p, qv0_p, ql0_p, qi0_p, nl0_p, ni0_p
         type(c_ptr), value :: dacudt_p, T1_p, qv1_p, ql1_p, qi1_p, nl1_p, ni1_p
      end subroutine cldwat2m_input_state_codon
   end interface

   call linear_state_select_impl()
   if (.not. use_native_linear_state_impl) then
      call linear_state_log_entered()
      call cldwat2m_input_state_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(top_lev, c_int64_t), dt, c_loc(a_cu0(1,1)), c_loc(a_cud(1,1)), &
           c_loc(T0(1,1)), c_loc(qv0(1,1)), c_loc(ql0(1,1)), c_loc(qi0(1,1)), c_loc(nl0(1,1)), &
           c_loc(ni0(1,1)), c_loc(dacudt(1,1)), c_loc(T1(1,1)), c_loc(qv1(1,1)), c_loc(ql1(1,1)), &
           c_loc(qi1(1,1)), c_loc(nl1(1,1)), c_loc(ni1(1,1)))
      return
   endif

   dacudt(:ncol,top_lev:pver) = &
        (a_cu0(:ncol,top_lev:pver) - a_cud(:ncol,top_lev:pver))/dt

   ql0(:ncol,:top_lev-1) = 0._r8
   qi0(:ncol,:top_lev-1) = 0._r8
   nl0(:ncol,:top_lev-1) = 0._r8
   ni0(:ncol,:top_lev-1) = 0._r8

   T1(:ncol,:)    =  T0(:ncol,:)
   qv1(:ncol,:)   = qv0(:ncol,:)
   ql1(:ncol,:)   = ql0(:ncol,:)
   qi1(:ncol,:)   = qi0(:ncol,:)
   nl1(:ncol,:)   = nl0(:ncol,:)
   ni1(:ncol,:)   = ni0(:ncol,:)

   end subroutine input_state_codon_wrap

   subroutine qmin_fill_codon_wrap(ncol, qvmin, qlmin, qimin, qmin1, qmin2, qmin3)

   integer, intent(in) :: ncol
   real(r8), intent(in) :: qvmin, qlmin, qimin
   real(r8), target, intent(out) :: qmin1(pcols,pver), qmin2(pcols,pver), qmin3(pcols,pver)

   interface
      subroutine cldwat2m_qmin_fill_codon(ncol_c, pcols_c, pver_c, qvmin_c, qlmin_c, qimin_c, &
           qmin1_p, qmin2_p, qmin3_p) bind(c, name="cldwat2m_qmin_fill_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: qvmin_c, qlmin_c, qimin_c
         type(c_ptr), value :: qmin1_p, qmin2_p, qmin3_p
      end subroutine cldwat2m_qmin_fill_codon
   end interface

   call linear_state_select_impl()
   if (.not. use_native_linear_state_impl) then
      call linear_state_log_entered()
      call cldwat2m_qmin_fill_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           qvmin, qlmin, qimin, c_loc(qmin1(1,1)), c_loc(qmin2(1,1)), c_loc(qmin3(1,1)))
      return
   endif

   qmin1(:ncol,:) = qvmin
   qmin2(:ncol,:) = qlmin
   qmin3(:ncol,:) = qimin

   end subroutine qmin_fill_codon_wrap

   subroutine detrain_state_codon_wrap(ncol, dt, T_prime0, qv_prime0, ql_prime0, qi_prime0, &
        nl_prime0, ni_prime0, D_T, D_qv, D_ql, D_qi, D_nl, D_ni, T_dprime, qv_dprime, &
        ql_dprime, qi_dprime, nl_dprime, ni_dprime)

   integer, intent(in) :: ncol
   real(r8), intent(in) :: dt
   real(r8), target, intent(in) :: T_prime0(pcols,pver), qv_prime0(pcols,pver), ql_prime0(pcols,pver)
   real(r8), target, intent(in) :: qi_prime0(pcols,pver), nl_prime0(pcols,pver), ni_prime0(pcols,pver)
   real(r8), target, intent(in) :: D_T(pcols,pver), D_qv(pcols,pver), D_ql(pcols,pver)
   real(r8), target, intent(in) :: D_qi(pcols,pver), D_nl(pcols,pver), D_ni(pcols,pver)
   real(r8), target, intent(out) :: T_dprime(pcols,pver), qv_dprime(pcols,pver), ql_dprime(pcols,pver)
   real(r8), target, intent(out) :: qi_dprime(pcols,pver), nl_dprime(pcols,pver), ni_dprime(pcols,pver)

   interface
      subroutine cldwat2m_detrain_state_codon(ncol_c, pcols_c, pver_c, top_lev_c, dt_c, &
           T_prime0_p, qv_prime0_p, ql_prime0_p, qi_prime0_p, nl_prime0_p, ni_prime0_p, &
           D_T_p, D_qv_p, D_ql_p, D_qi_p, D_nl_p, D_ni_p, &
           T_dprime_p, qv_dprime_p, ql_dprime_p, qi_dprime_p, nl_dprime_p, ni_dprime_p) &
           bind(c, name="cldwat2m_detrain_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: dt_c
         type(c_ptr), value :: T_prime0_p, qv_prime0_p, ql_prime0_p, qi_prime0_p, nl_prime0_p, ni_prime0_p
         type(c_ptr), value :: D_T_p, D_qv_p, D_ql_p, D_qi_p, D_nl_p, D_ni_p
         type(c_ptr), value :: T_dprime_p, qv_dprime_p, ql_dprime_p, qi_dprime_p, nl_dprime_p, ni_dprime_p
      end subroutine cldwat2m_detrain_state_codon
   end interface

   call linear_state_select_impl()
   if (.not. use_native_linear_state_impl) then
      call linear_state_log_entered()
      call cldwat2m_detrain_state_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(top_lev, c_int64_t), dt, c_loc(T_prime0(1,1)), &
           c_loc(qv_prime0(1,1)), c_loc(ql_prime0(1,1)), c_loc(qi_prime0(1,1)), &
           c_loc(nl_prime0(1,1)), c_loc(ni_prime0(1,1)), c_loc(D_T(1,1)), c_loc(D_qv(1,1)), &
           c_loc(D_ql(1,1)), c_loc(D_qi(1,1)), c_loc(D_nl(1,1)), c_loc(D_ni(1,1)), &
           c_loc(T_dprime(1,1)), c_loc(qv_dprime(1,1)), c_loc(ql_dprime(1,1)), &
           c_loc(qi_dprime(1,1)), c_loc(nl_dprime(1,1)), c_loc(ni_dprime(1,1)))
      return
   endif

   T_dprime(:ncol,top_lev:)  =  T_prime0(:ncol,top_lev:)
   qv_dprime(:ncol,top_lev:) = qv_prime0(:ncol,top_lev:)
   ql_dprime(:ncol,top_lev:) = ql_prime0(:ncol,top_lev:)
   qi_dprime(:ncol,top_lev:) = qi_prime0(:ncol,top_lev:)
   nl_dprime(:ncol,top_lev:) = nl_prime0(:ncol,top_lev:)
   ni_dprime(:ncol,top_lev:) = ni_prime0(:ncol,top_lev:)

   T_dprime(:ncol,top_lev:)   =  T_dprime(:ncol,top_lev:)  + D_T(:ncol,top_lev:) * dt
   qv_dprime(:ncol,top_lev:)  = qv_dprime(:ncol,top_lev:) + D_qv(:ncol,top_lev:) * dt
   ql_dprime(:ncol,top_lev:)  = ql_dprime(:ncol,top_lev:) + D_ql(:ncol,top_lev:) * dt
   qi_dprime(:ncol,top_lev:)  = qi_dprime(:ncol,top_lev:) + D_qi(:ncol,top_lev:) * dt
   nl_dprime(:ncol,top_lev:)  = nl_dprime(:ncol,top_lev:) + D_nl(:ncol,top_lev:) * dt
   ni_dprime(:ncol,top_lev:)  = ni_dprime(:ncol,top_lev:) + D_ni(:ncol,top_lev:) * dt

   end subroutine detrain_state_codon_wrap

!=======================================================================================================

   subroutine qq_limiter_select_impl()

   character(len=32) :: impl_name
   integer :: n, status

   if (qq_limiter_impl_selected) return
   call get_environment_variable('CLDWAT2M_QQ_LIMITER_IMPL', value=impl_name, length=n, status=status)
   use_native_qq_limiter_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_qq_limiter_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_qq_limiter_impl = .false.
      case default
         use_native_qq_limiter_impl = .false.
      end select
   end if
   qq_limiter_impl_selected = .true.
   if (masterproc) then
      if (use_native_qq_limiter_impl) then
         write(iulog,*) 'cldwat2m_qq_limiter implementation = native'
      else
         write(iulog,*) 'cldwat2m_qq_limiter implementation = codon'
      end if
   end if
   end subroutine qq_limiter_select_impl

   subroutine qq_limiter_log_entered()
   if (qq_limiter_entered_logged) return
   qq_limiter_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_qq_limiter entered (macrophysics QQ limiter = codon)'
   end if
   end subroutine qq_limiter_log_entered

   subroutine qq_limiter_codon_wrap(ncol, dt, cone, qvmin, qv_05, ql_05, qi_05, nl_05, ni_05, &
        qsat_a, qvwb_aw, QQ, QQw, QQi, QQnl, QQni)

   integer, intent(in) :: ncol
   real(r8), intent(in) :: dt, cone, qvmin
   real(r8), target, intent(in) :: qv_05(pcols,pver), ql_05(pcols,pver), qi_05(pcols,pver)
   real(r8), target, intent(in) :: nl_05(pcols,pver), ni_05(pcols,pver), qsat_a(pcols,pver)
   real(r8), target, intent(in) :: qvwb_aw(pcols,pver)
   real(r8), target, intent(inout) :: QQ(pcols,pver), QQw(pcols,pver), QQi(pcols,pver)
   real(r8), target, intent(inout) :: QQnl(pcols,pver), QQni(pcols,pver)

   integer :: i, k
   real(r8) :: QQmax, QQmin, QQwmin, QQimin

   interface
      subroutine cldwat2m_qq_limiter_codon(ncol_c, pcols_c, pver_c, top_lev_c, dt_c, qsmall_c, &
           cone_c, qvmin_c, qv_05_p, ql_05_p, qi_05_p, nl_05_p, ni_05_p, qsat_a_p, qvwb_aw_p, &
           QQ_p, QQw_p, QQi_p, QQnl_p, QQni_p) bind(c, name="cldwat2m_qq_limiter_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: dt_c, qsmall_c, cone_c, qvmin_c
         type(c_ptr), value :: qv_05_p, ql_05_p, qi_05_p, nl_05_p, ni_05_p, qsat_a_p, qvwb_aw_p
         type(c_ptr), value :: QQ_p, QQw_p, QQi_p, QQnl_p, QQni_p
      end subroutine cldwat2m_qq_limiter_codon
   end interface

   call qq_limiter_select_impl()
   if (.not. use_native_qq_limiter_impl) then
      call qq_limiter_log_entered()
      call cldwat2m_qq_limiter_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(top_lev, c_int64_t), dt, qsmall, cone, qvmin, &
           c_loc(qv_05(1,1)), c_loc(ql_05(1,1)), c_loc(qi_05(1,1)), c_loc(nl_05(1,1)), &
           c_loc(ni_05(1,1)), c_loc(qsat_a(1,1)), c_loc(qvwb_aw(1,1)), c_loc(QQ(1,1)), &
           c_loc(QQw(1,1)), c_loc(QQi(1,1)), c_loc(QQnl(1,1)), c_loc(QQni(1,1)))
      return
   endif

   do k = top_lev, pver
      do i = 1, ncol
         QQnl(i,k) = 0._r8
         QQni(i,k) = 0._r8
         if( QQ(i,k) .ge. 0._r8 ) then
             QQmax    = (qv_05(i,k) - qvmin)/dt
             QQmax    = max(0._r8,QQmax)
             QQ(i,k)  = min(QQ(i,k),QQmax)
             QQw(i,k) = QQ(i,k)
             QQi(i,k) = 0._r8
         else
             QQmin  = 0._r8
             if( qv_05(i,k) .lt. qsat_a(i,k) ) QQmin = min(0._r8,cone*(qv_05(i,k)-qvwb_aw(i,k))/dt)
             QQ(i,k)  = max(QQ(i,k),QQmin)
             QQw(i,k) = QQ(i,k)
             QQi(i,k) = 0._r8
             QQwmin   = min(0._r8,-cone*ql_05(i,k)/dt)
             QQimin   = min(0._r8,-cone*qi_05(i,k)/dt)
             QQw(i,k) = min(0._r8,max(QQw(i,k),QQwmin))
             QQi(i,k) = min(0._r8,max(QQi(i,k),QQimin))
         endif

         if( QQw(i,k) .lt. 0._r8 ) then
             if( ql_05(i,k) .gt. qsmall ) then
                 QQnl(i,k) = QQw(i,k)*nl_05(i,k)/ql_05(i,k)
                 QQnl(i,k) = min(0._r8,cone*max(QQnl(i,k),-nl_05(i,k)/dt))
             else
                 QQnl(i,k) = 0._r8
             endif
         endif

         if( QQi(i,k) .lt. 0._r8 ) then
             if( qi_05(i,k) .gt. qsmall ) then
                 QQni(i,k) = QQi(i,k)*ni_05(i,k)/qi_05(i,k)
                 QQni(i,k) = min(0._r8,cone*max(QQni(i,k),-ni_05(i,k)/dt))
             else
                 QQni(i,k) = 0._r8
             endif
         endif
      enddo
   enddo

   end subroutine qq_limiter_codon_wrap

!=======================================================================================================

   subroutine advective_state_select_impl()

   character(len=32) :: impl_name
   integer :: n, status

   if (advective_state_impl_selected) return
   call get_environment_variable('CLDWAT2M_ADVECTIVE_STATE_IMPL', value=impl_name, length=n, status=status)
   use_native_advective_state_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_advective_state_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_advective_state_impl = .false.
      case default
         use_native_advective_state_impl = .false.
      end select
   end if
   advective_state_impl_selected = .true.
   if (masterproc) then
      if (use_native_advective_state_impl) then
         write(iulog,*) 'cldwat2m_advective_state implementation = native'
      else
         write(iulog,*) 'cldwat2m_advective_state implementation = codon'
      end if
   end if
   end subroutine advective_state_select_impl

   subroutine advective_state_log_entered()
   if (advective_state_entered_logged) return
   advective_state_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_advective_state entered (macrophysics advective state update = codon)'
   end if
   end subroutine advective_state_log_entered

   subroutine advective_state_codon_wrap(ncol, dt, T_0, qv_0, ql_0, qi_0, nl_0, ni_0, &
        A_T, C_T, A_qv, C_qv, A_ql, C_ql, A_qi, C_qi, A_nl, C_nl, A_ni, C_ni, &
        T_05, qv_05, ql_05, qi_05, nl_05, ni_05)

   integer, intent(in) :: ncol
   real(r8), intent(in) :: dt
   real(r8), target, intent(in) :: T_0(pcols,pver), qv_0(pcols,pver), ql_0(pcols,pver)
   real(r8), target, intent(in) :: qi_0(pcols,pver), nl_0(pcols,pver), ni_0(pcols,pver)
   real(r8), target, intent(in) :: A_T(pcols,pver), C_T(pcols,pver), A_qv(pcols,pver), C_qv(pcols,pver)
   real(r8), target, intent(in) :: A_ql(pcols,pver), C_ql(pcols,pver), A_qi(pcols,pver), C_qi(pcols,pver)
   real(r8), target, intent(in) :: A_nl(pcols,pver), C_nl(pcols,pver), A_ni(pcols,pver), C_ni(pcols,pver)
   real(r8), target, intent(out) :: T_05(pcols,pver), qv_05(pcols,pver), ql_05(pcols,pver)
   real(r8), target, intent(out) :: qi_05(pcols,pver), nl_05(pcols,pver), ni_05(pcols,pver)

   integer :: i, k

   interface
      subroutine cldwat2m_advective_state_codon(ncol_c, pcols_c, pver_c, top_lev_c, dt_c, &
           T_0_p, qv_0_p, ql_0_p, qi_0_p, nl_0_p, ni_0_p, &
           A_T_p, C_T_p, A_qv_p, C_qv_p, A_ql_p, C_ql_p, A_qi_p, C_qi_p, A_nl_p, C_nl_p, A_ni_p, C_ni_p, &
           T_05_p, qv_05_p, ql_05_p, qi_05_p, nl_05_p, ni_05_p) bind(c, name="cldwat2m_advective_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: dt_c
         type(c_ptr), value :: T_0_p, qv_0_p, ql_0_p, qi_0_p, nl_0_p, ni_0_p
         type(c_ptr), value :: A_T_p, C_T_p, A_qv_p, C_qv_p, A_ql_p, C_ql_p, A_qi_p, C_qi_p
         type(c_ptr), value :: A_nl_p, C_nl_p, A_ni_p, C_ni_p, T_05_p, qv_05_p, ql_05_p
         type(c_ptr), value :: qi_05_p, nl_05_p, ni_05_p
      end subroutine cldwat2m_advective_state_codon
   end interface

   call advective_state_select_impl()
   if (.not. use_native_advective_state_impl) then
      call advective_state_log_entered()
      call cldwat2m_advective_state_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(top_lev, c_int64_t), dt, &
           c_loc(T_0(1,1)), c_loc(qv_0(1,1)), c_loc(ql_0(1,1)), c_loc(qi_0(1,1)), &
           c_loc(nl_0(1,1)), c_loc(ni_0(1,1)), c_loc(A_T(1,1)), c_loc(C_T(1,1)), &
           c_loc(A_qv(1,1)), c_loc(C_qv(1,1)), c_loc(A_ql(1,1)), c_loc(C_ql(1,1)), &
           c_loc(A_qi(1,1)), c_loc(C_qi(1,1)), c_loc(A_nl(1,1)), c_loc(C_nl(1,1)), &
           c_loc(A_ni(1,1)), c_loc(C_ni(1,1)), c_loc(T_05(1,1)), c_loc(qv_05(1,1)), &
           c_loc(ql_05(1,1)), c_loc(qi_05(1,1)), c_loc(nl_05(1,1)), c_loc(ni_05(1,1)))
      return
   endif

   do k = top_lev, pver
      do i = 1, ncol
         T_05(i,k) = T_0(i,k) + ( A_T(i,k) + C_T(i,k) ) * dt
         qv_05(i,k) = qv_0(i,k) + ( A_qv(i,k) + C_qv(i,k) ) * dt
         ql_05(i,k) = ql_0(i,k) + ( A_ql(i,k) + C_ql(i,k) ) * dt
         qi_05(i,k) = qi_0(i,k) + ( A_qi(i,k) + C_qi(i,k) ) * dt
         nl_05(i,k) = max(0._r8, nl_0(i,k) + ( A_nl(i,k) + C_nl(i,k) ) * dt )
         ni_05(i,k) = max(0._r8, ni_0(i,k) + ( A_ni(i,k) + C_ni(i,k) ) * dt )
      enddo
   enddo

   end subroutine advective_state_codon_wrap

!=======================================================================================================

   subroutine iter_state_select_impl()

   character(len=32) :: impl_name
   integer :: n, status

   if (iter_state_impl_selected) return
   call get_environment_variable('CLDWAT2M_ITER_STATE_IMPL', value=impl_name, length=n, status=status)
   use_native_iter_state_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_iter_state_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_iter_state_impl = .false.
      case default
         use_native_iter_state_impl = .false.
      end select
   end if
   iter_state_impl_selected = .true.
   if (masterproc) then
      if (use_native_iter_state_impl) then
         write(iulog,*) 'cldwat2m_iter_state implementation = native'
      else
         write(iulog,*) 'cldwat2m_iter_state implementation = codon'
      end if
   end if
   end subroutine iter_state_select_impl

   subroutine iter_state_log_entered()
   if (iter_state_entered_logged) return
   iter_state_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_iter_state entered (macrophysics iterative prognostic state = codon)'
   end if
   end subroutine iter_state_log_entered

   subroutine iter_state_codon_wrap(iter, ncol, dt, &
        QQw, QQi, QQnl, QQni, QQw_prev, QQi_prev, QQnl_prev, QQni_prev, &
        QQw_prog, QQi_prog, QQnl_prog, QQni_prog, &
        T_0, qv_0, ql_0, qi_0, nl_0, ni_0, &
        A_T, A_T_adj, C_T, A_qv, A_qv_adj, C_qv, A_ql, A_ql_adj, C_ql, &
        A_qi, A_qi_adj, C_qi, A_nl, C_nl, A_ni, C_ni, &
        T_prime0, qv_prime0, ql_prime0, qi_prime0, nl_prime0, ni_prime0)

   integer, intent(in) :: iter
   integer, intent(in) :: ncol
   real(r8), intent(in) :: dt
   real(r8), target, intent(in) :: QQw(pcols,pver), QQi(pcols,pver), QQnl(pcols,pver), QQni(pcols,pver)
   real(r8), target, intent(inout) :: QQw_prev(pcols,pver), QQi_prev(pcols,pver)
   real(r8), target, intent(inout) :: QQnl_prev(pcols,pver), QQni_prev(pcols,pver)
   real(r8), target, intent(out) :: QQw_prog(pcols,pver), QQi_prog(pcols,pver)
   real(r8), target, intent(out) :: QQnl_prog(pcols,pver), QQni_prog(pcols,pver)
   real(r8), target, intent(in) :: T_0(pcols,pver), qv_0(pcols,pver), ql_0(pcols,pver)
   real(r8), target, intent(in) :: qi_0(pcols,pver), nl_0(pcols,pver), ni_0(pcols,pver)
   real(r8), target, intent(in) :: A_T(pcols,pver), A_T_adj(pcols,pver), C_T(pcols,pver)
   real(r8), target, intent(in) :: A_qv(pcols,pver), A_qv_adj(pcols,pver), C_qv(pcols,pver)
   real(r8), target, intent(in) :: A_ql(pcols,pver), A_ql_adj(pcols,pver), C_ql(pcols,pver)
   real(r8), target, intent(in) :: A_qi(pcols,pver), A_qi_adj(pcols,pver), C_qi(pcols,pver)
   real(r8), target, intent(in) :: A_nl(pcols,pver), C_nl(pcols,pver), A_ni(pcols,pver), C_ni(pcols,pver)
   real(r8), target, intent(out) :: T_prime0(pcols,pver), qv_prime0(pcols,pver), ql_prime0(pcols,pver)
   real(r8), target, intent(out) :: qi_prime0(pcols,pver), nl_prime0(pcols,pver), ni_prime0(pcols,pver)

   integer :: i, k

   interface
      subroutine cldwat2m_iter_state_codon(iter_c, ncol_c, pcols_c, pver_c, top_lev_c, &
           dt_c, ramda_c, qsmall_c, latvap_c, latice_c, cpair_c, &
           QQw_p, QQi_p, QQnl_p, QQni_p, QQw_prev_p, QQi_prev_p, QQnl_prev_p, QQni_prev_p, &
           QQw_prog_p, QQi_prog_p, QQnl_prog_p, QQni_prog_p, &
           T_0_p, qv_0_p, ql_0_p, qi_0_p, nl_0_p, ni_0_p, &
           A_T_p, A_T_adj_p, C_T_p, A_qv_p, A_qv_adj_p, C_qv_p, A_ql_p, A_ql_adj_p, C_ql_p, &
           A_qi_p, A_qi_adj_p, C_qi_p, A_nl_p, C_nl_p, A_ni_p, C_ni_p, &
           T_prime0_p, qv_prime0_p, ql_prime0_p, qi_prime0_p, nl_prime0_p, ni_prime0_p) &
           bind(c, name="cldwat2m_iter_state_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: iter_c, ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: dt_c, ramda_c, qsmall_c, latvap_c, latice_c, cpair_c
         type(c_ptr), value :: QQw_p, QQi_p, QQnl_p, QQni_p, QQw_prev_p, QQi_prev_p, QQnl_prev_p, QQni_prev_p
         type(c_ptr), value :: QQw_prog_p, QQi_prog_p, QQnl_prog_p, QQni_prog_p
         type(c_ptr), value :: T_0_p, qv_0_p, ql_0_p, qi_0_p, nl_0_p, ni_0_p
         type(c_ptr), value :: A_T_p, A_T_adj_p, C_T_p, A_qv_p, A_qv_adj_p, C_qv_p
         type(c_ptr), value :: A_ql_p, A_ql_adj_p, C_ql_p, A_qi_p, A_qi_adj_p, C_qi_p
         type(c_ptr), value :: A_nl_p, C_nl_p, A_ni_p, C_ni_p
         type(c_ptr), value :: T_prime0_p, qv_prime0_p, ql_prime0_p, qi_prime0_p, nl_prime0_p, ni_prime0_p
      end subroutine cldwat2m_iter_state_codon
   end interface

   call iter_state_select_impl()
   if (.not. use_native_iter_state_impl) then
      call iter_state_log_entered()
      call cldwat2m_iter_state_codon(int(iter, c_int64_t), int(ncol, c_int64_t), &
           int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
           dt, ramda, qsmall, latvap, latice, cpair, &
           c_loc(QQw(1,1)), c_loc(QQi(1,1)), c_loc(QQnl(1,1)), c_loc(QQni(1,1)), &
           c_loc(QQw_prev(1,1)), c_loc(QQi_prev(1,1)), c_loc(QQnl_prev(1,1)), c_loc(QQni_prev(1,1)), &
           c_loc(QQw_prog(1,1)), c_loc(QQi_prog(1,1)), c_loc(QQnl_prog(1,1)), c_loc(QQni_prog(1,1)), &
           c_loc(T_0(1,1)), c_loc(qv_0(1,1)), c_loc(ql_0(1,1)), c_loc(qi_0(1,1)), &
           c_loc(nl_0(1,1)), c_loc(ni_0(1,1)), c_loc(A_T(1,1)), c_loc(A_T_adj(1,1)), c_loc(C_T(1,1)), &
           c_loc(A_qv(1,1)), c_loc(A_qv_adj(1,1)), c_loc(C_qv(1,1)), &
           c_loc(A_ql(1,1)), c_loc(A_ql_adj(1,1)), c_loc(C_ql(1,1)), &
           c_loc(A_qi(1,1)), c_loc(A_qi_adj(1,1)), c_loc(C_qi(1,1)), &
           c_loc(A_nl(1,1)), c_loc(C_nl(1,1)), c_loc(A_ni(1,1)), c_loc(C_ni(1,1)), &
           c_loc(T_prime0(1,1)), c_loc(qv_prime0(1,1)), c_loc(ql_prime0(1,1)), c_loc(qi_prime0(1,1)), &
           c_loc(nl_prime0(1,1)), c_loc(ni_prime0(1,1)))
      return
   endif

   if( iter .eq. 1 ) then
      do k = top_lev, pver
         do i = 1, ncol
            QQw_prev(i,k) = QQw(i,k)
            QQi_prev(i,k) = QQi(i,k)
            QQnl_prev(i,k) = QQnl(i,k)
            QQni_prev(i,k) = QQni(i,k)
         enddo
      enddo
   endif

   do k = top_lev, pver
      do i = 1, ncol
         QQw_prog(i,k) = ramda*QQw(i,k) + (1._r8-ramda)*QQw_prev(i,k)
         QQi_prog(i,k) = ramda*QQi(i,k) + (1._r8-ramda)*QQi_prev(i,k)
         QQnl_prog(i,k) = ramda*QQnl(i,k) + (1._r8-ramda)*QQnl_prev(i,k)
         QQni_prog(i,k) = ramda*QQni(i,k) + (1._r8-ramda)*QQni_prev(i,k)
         QQw_prev(i,k) = QQw_prog(i,k)
         QQi_prev(i,k) = QQi_prog(i,k)
         QQnl_prev(i,k) = QQnl_prog(i,k)
         QQni_prev(i,k) = QQni_prog(i,k)
         T_prime0(i,k) = T_0(i,k) + dt*( A_T(i,k) + A_T_adj(i,k) + C_T(i,k) + &
              (latvap*QQw_prog(i,k)+(latvap+latice)*QQi_prog(i,k))/cpair )
         qv_prime0(i,k) = qv_0(i,k) + dt*( A_qv(i,k) + A_qv_adj(i,k) + C_qv(i,k) - QQw_prog(i,k) - QQi_prog(i,k) )
         ql_prime0(i,k) = ql_0(i,k) + dt*( A_ql(i,k) + A_ql_adj(i,k) + C_ql(i,k) + QQw_prog(i,k) )
         qi_prime0(i,k) = qi_0(i,k) + dt*( A_qi(i,k) + A_qi_adj(i,k) + C_qi(i,k) + QQi_prog(i,k) )
         nl_prime0(i,k) = max(0._r8,nl_0(i,k) + dt*( A_nl(i,k) + C_nl(i,k) + QQnl_prog(i,k) ))
         ni_prime0(i,k) = max(0._r8,ni_0(i,k) + dt*( A_ni(i,k) + C_ni(i,k) + QQni_prog(i,k) ))
         if( ql_prime0(i,k) .lt. qsmall ) nl_prime0(i,k) = 0._r8
         if( qi_prime0(i,k) .lt. qsmall ) ni_prime0(i,k) = 0._r8
      enddo
   enddo

   end subroutine iter_state_codon_wrap

!=======================================================================================================

   subroutine final_tendency_select_impl()

   character(len=32) :: impl_name
   integer :: n, status

   if (final_tendency_impl_selected) return
   call get_environment_variable('CLDWAT2M_FINAL_TENDENCY_IMPL', value=impl_name, length=n, status=status)
   use_native_final_tendency_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_final_tendency_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_final_tendency_impl = .false.
      case default
         use_native_final_tendency_impl = .false.
      end select
   end if
   final_tendency_impl_selected = .true.
   if (masterproc) then
      if (use_native_final_tendency_impl) then
         write(iulog,*) 'cldwat2m_final_tendency implementation = native'
      else
         write(iulog,*) 'cldwat2m_final_tendency implementation = codon'
      end if
   end if
   end subroutine final_tendency_select_impl

   subroutine final_tendency_log_entered()
   if (final_tendency_entered_logged) return
   final_tendency_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_final_tendency entered (macrophysics final tendency/update = codon)'
   end if
   end subroutine final_tendency_log_entered

   subroutine final_tendency_codon_wrap(ncol, dt, do_cldice, &
        QQw_prog, QQi_prog, QQnl_prog, QQni_prog, QQw1, QQi1, QQw2, QQi2, &
        qlten_pwi1, qlten_pwi2, qiten_pwi1, qiten_pwi2, A_ql_adj, A_qi_adj, &
        QQnl1, QQni1, QQnl2, QQni2, nlten_pwi1, nlten_pwi2, niten_pwi1, niten_pwi2, &
        ACnl, ACni, A_nl_adj, A_ni_adj, qvten_pwi1, qvten_pwi2, A_qv_adj, &
        T_star, qv_star, ql_star, qi_star, nl_star, ni_star, &
        A_T, C_T, A_qv, C_qv, A_ql, C_ql, A_qi, C_qi, A_nl, C_nl, A_ni, C_ni, &
        a_st_star, a_cu0, QQw_final, QQi_final, QQ_final, QQw_all, QQi_all, QQ_all, &
        QQnl_final, QQni_final, QQn_final, QQnl_all, QQni_all, QQn_all, &
        qme, qvadj, qladj, qiadj, qllim, qilim, &
        s_tendout, qv_tendout, ql_tendout, qi_tendout, nl_tendout, ni_tendout, cld, &
        T0, qv0, ql0, qi0, nl0, ni0)

   integer, intent(in) :: ncol
   real(r8), intent(in) :: dt
   logical, intent(in) :: do_cldice
   real(r8), target, intent(in) :: QQw_prog(pcols,pver), QQi_prog(pcols,pver)
   real(r8), target, intent(in) :: QQnl_prog(pcols,pver), QQni_prog(pcols,pver)
   real(r8), target, intent(in) :: QQw1(pcols,pver), QQi1(pcols,pver), QQw2(pcols,pver), QQi2(pcols,pver)
   real(r8), target, intent(in) :: qlten_pwi1(pcols,pver), qlten_pwi2(pcols,pver)
   real(r8), target, intent(in) :: qiten_pwi1(pcols,pver), qiten_pwi2(pcols,pver)
   real(r8), target, intent(in) :: A_ql_adj(pcols,pver), A_qi_adj(pcols,pver)
   real(r8), target, intent(in) :: QQnl1(pcols,pver), QQni1(pcols,pver), QQnl2(pcols,pver), QQni2(pcols,pver)
   real(r8), target, intent(in) :: nlten_pwi1(pcols,pver), nlten_pwi2(pcols,pver)
   real(r8), target, intent(in) :: niten_pwi1(pcols,pver), niten_pwi2(pcols,pver)
   real(r8), target, intent(in) :: ACnl(pcols,pver), ACni(pcols,pver), A_nl_adj(pcols,pver), A_ni_adj(pcols,pver)
   real(r8), target, intent(in) :: qvten_pwi1(pcols,pver), qvten_pwi2(pcols,pver), A_qv_adj(pcols,pver)
   real(r8), target, intent(in) :: T_star(pcols,pver), qv_star(pcols,pver), ql_star(pcols,pver)
   real(r8), target, intent(in) :: qi_star(pcols,pver), nl_star(pcols,pver), ni_star(pcols,pver)
   real(r8), target, intent(in) :: A_T(pcols,pver), C_T(pcols,pver), A_qv(pcols,pver), C_qv(pcols,pver)
   real(r8), target, intent(in) :: A_ql(pcols,pver), C_ql(pcols,pver), A_qi(pcols,pver), C_qi(pcols,pver)
   real(r8), target, intent(in) :: A_nl(pcols,pver), C_nl(pcols,pver), A_ni(pcols,pver), C_ni(pcols,pver)
   real(r8), target, intent(in) :: a_st_star(pcols,pver), a_cu0(pcols,pver)
   real(r8), target, intent(out) :: QQw_final(pcols,pver), QQi_final(pcols,pver), QQ_final(pcols,pver)
   real(r8), target, intent(out) :: QQw_all(pcols,pver), QQi_all(pcols,pver), QQ_all(pcols,pver)
   real(r8), target, intent(out) :: QQnl_final(pcols,pver), QQni_final(pcols,pver), QQn_final(pcols,pver)
   real(r8), target, intent(out) :: QQnl_all(pcols,pver), QQni_all(pcols,pver), QQn_all(pcols,pver)
   real(r8), target, intent(out) :: qme(pcols,pver), qvadj(pcols,pver), qladj(pcols,pver)
   real(r8), target, intent(out) :: qiadj(pcols,pver), qllim(pcols,pver), qilim(pcols,pver)
   real(r8), target, intent(out) :: s_tendout(pcols,pver), qv_tendout(pcols,pver), ql_tendout(pcols,pver)
   real(r8), target, intent(out) :: qi_tendout(pcols,pver), nl_tendout(pcols,pver), ni_tendout(pcols,pver)
   real(r8), target, intent(out) :: cld(pcols,pver)
   real(r8), target, intent(inout) :: T0(pcols,pver), qv0(pcols,pver), ql0(pcols,pver)
   real(r8), target, intent(inout) :: qi0(pcols,pver), nl0(pcols,pver), ni0(pcols,pver)

   integer :: i, k
   integer(c_int64_t) :: do_cldice_c

   interface
      subroutine cldwat2m_final_tendency_codon(ncol_c, pcols_c, pver_c, top_lev_c, do_cldice_c, dt_c, cpair_c, &
           QQw_prog_p, QQi_prog_p, QQnl_prog_p, QQni_prog_p, QQw1_p, QQi1_p, QQw2_p, QQi2_p, &
           qlten_pwi1_p, qlten_pwi2_p, qiten_pwi1_p, qiten_pwi2_p, A_ql_adj_p, A_qi_adj_p, &
           QQnl1_p, QQni1_p, QQnl2_p, QQni2_p, nlten_pwi1_p, nlten_pwi2_p, niten_pwi1_p, niten_pwi2_p, &
           ACnl_p, ACni_p, A_nl_adj_p, A_ni_adj_p, qvten_pwi1_p, qvten_pwi2_p, A_qv_adj_p, &
           T_star_p, qv_star_p, ql_star_p, qi_star_p, nl_star_p, ni_star_p, &
           A_T_p, C_T_p, A_qv_p, C_qv_p, A_ql_p, C_ql_p, A_qi_p, C_qi_p, A_nl_p, C_nl_p, A_ni_p, C_ni_p, &
           a_st_star_p, a_cu0_p, QQw_final_p, QQi_final_p, QQ_final_p, QQw_all_p, QQi_all_p, QQ_all_p, &
           QQnl_final_p, QQni_final_p, QQn_final_p, QQnl_all_p, QQni_all_p, QQn_all_p, &
           qme_p, qvadj_p, qladj_p, qiadj_p, qllim_p, qilim_p, &
           s_tendout_p, qv_tendout_p, ql_tendout_p, qi_tendout_p, nl_tendout_p, ni_tendout_p, cld_p, &
           T0_p, qv0_p, ql0_p, qi0_p, nl0_p, ni0_p) bind(c, name="cldwat2m_final_tendency_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, do_cldice_c
         real(c_double), value :: dt_c, cpair_c
         type(c_ptr), value :: QQw_prog_p, QQi_prog_p, QQnl_prog_p, QQni_prog_p, QQw1_p, QQi1_p, QQw2_p, QQi2_p
         type(c_ptr), value :: qlten_pwi1_p, qlten_pwi2_p, qiten_pwi1_p, qiten_pwi2_p, A_ql_adj_p, A_qi_adj_p
         type(c_ptr), value :: QQnl1_p, QQni1_p, QQnl2_p, QQni2_p, nlten_pwi1_p, nlten_pwi2_p
         type(c_ptr), value :: niten_pwi1_p, niten_pwi2_p, ACnl_p, ACni_p, A_nl_adj_p, A_ni_adj_p
         type(c_ptr), value :: qvten_pwi1_p, qvten_pwi2_p, A_qv_adj_p
         type(c_ptr), value :: T_star_p, qv_star_p, ql_star_p, qi_star_p, nl_star_p, ni_star_p
         type(c_ptr), value :: A_T_p, C_T_p, A_qv_p, C_qv_p, A_ql_p, C_ql_p, A_qi_p, C_qi_p
         type(c_ptr), value :: A_nl_p, C_nl_p, A_ni_p, C_ni_p, a_st_star_p, a_cu0_p
         type(c_ptr), value :: QQw_final_p, QQi_final_p, QQ_final_p, QQw_all_p, QQi_all_p, QQ_all_p
         type(c_ptr), value :: QQnl_final_p, QQni_final_p, QQn_final_p, QQnl_all_p, QQni_all_p, QQn_all_p
         type(c_ptr), value :: qme_p, qvadj_p, qladj_p, qiadj_p, qllim_p, qilim_p
         type(c_ptr), value :: s_tendout_p, qv_tendout_p, ql_tendout_p, qi_tendout_p
         type(c_ptr), value :: nl_tendout_p, ni_tendout_p, cld_p
         type(c_ptr), value :: T0_p, qv0_p, ql0_p, qi0_p, nl0_p, ni0_p
      end subroutine cldwat2m_final_tendency_codon
   end interface

   call final_tendency_select_impl()
   if (.not. use_native_final_tendency_impl) then
      do_cldice_c = 0_c_int64_t
      if (do_cldice) do_cldice_c = 1_c_int64_t
      call final_tendency_log_entered()
      call cldwat2m_final_tendency_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(top_lev, c_int64_t), do_cldice_c, dt, cpair, &
           c_loc(QQw_prog(1,1)), c_loc(QQi_prog(1,1)), c_loc(QQnl_prog(1,1)), c_loc(QQni_prog(1,1)), &
           c_loc(QQw1(1,1)), c_loc(QQi1(1,1)), c_loc(QQw2(1,1)), c_loc(QQi2(1,1)), &
           c_loc(qlten_pwi1(1,1)), c_loc(qlten_pwi2(1,1)), c_loc(qiten_pwi1(1,1)), c_loc(qiten_pwi2(1,1)), &
           c_loc(A_ql_adj(1,1)), c_loc(A_qi_adj(1,1)), c_loc(QQnl1(1,1)), c_loc(QQni1(1,1)), &
           c_loc(QQnl2(1,1)), c_loc(QQni2(1,1)), c_loc(nlten_pwi1(1,1)), c_loc(nlten_pwi2(1,1)), &
           c_loc(niten_pwi1(1,1)), c_loc(niten_pwi2(1,1)), c_loc(ACnl(1,1)), c_loc(ACni(1,1)), &
           c_loc(A_nl_adj(1,1)), c_loc(A_ni_adj(1,1)), c_loc(qvten_pwi1(1,1)), c_loc(qvten_pwi2(1,1)), &
           c_loc(A_qv_adj(1,1)), c_loc(T_star(1,1)), c_loc(qv_star(1,1)), c_loc(ql_star(1,1)), &
           c_loc(qi_star(1,1)), c_loc(nl_star(1,1)), c_loc(ni_star(1,1)), c_loc(A_T(1,1)), c_loc(C_T(1,1)), &
           c_loc(A_qv(1,1)), c_loc(C_qv(1,1)), c_loc(A_ql(1,1)), c_loc(C_ql(1,1)), &
           c_loc(A_qi(1,1)), c_loc(C_qi(1,1)), c_loc(A_nl(1,1)), c_loc(C_nl(1,1)), &
           c_loc(A_ni(1,1)), c_loc(C_ni(1,1)), c_loc(a_st_star(1,1)), c_loc(a_cu0(1,1)), &
           c_loc(QQw_final(1,1)), c_loc(QQi_final(1,1)), c_loc(QQ_final(1,1)), &
           c_loc(QQw_all(1,1)), c_loc(QQi_all(1,1)), c_loc(QQ_all(1,1)), &
           c_loc(QQnl_final(1,1)), c_loc(QQni_final(1,1)), c_loc(QQn_final(1,1)), &
           c_loc(QQnl_all(1,1)), c_loc(QQni_all(1,1)), c_loc(QQn_all(1,1)), &
           c_loc(qme(1,1)), c_loc(qvadj(1,1)), c_loc(qladj(1,1)), c_loc(qiadj(1,1)), &
           c_loc(qllim(1,1)), c_loc(qilim(1,1)), c_loc(s_tendout(1,1)), c_loc(qv_tendout(1,1)), &
           c_loc(ql_tendout(1,1)), c_loc(qi_tendout(1,1)), c_loc(nl_tendout(1,1)), c_loc(ni_tendout(1,1)), &
           c_loc(cld(1,1)), c_loc(T0(1,1)), c_loc(qv0(1,1)), c_loc(ql0(1,1)), &
           c_loc(qi0(1,1)), c_loc(nl0(1,1)), c_loc(ni0(1,1)))
      return
   endif

   do k = top_lev, pver
      do i = 1, ncol
         QQw_final(i,k) = QQw_prog(i,k)
         QQi_final(i,k) = QQi_prog(i,k)
         QQ_final(i,k) = QQw_final(i,k) + QQi_final(i,k)
         QQw_all(i,k) = QQw_prog(i,k) + QQw1(i,k) + QQw2(i,k) + qlten_pwi1(i,k) + &
              qlten_pwi2(i,k) + A_ql_adj(i,k)
         QQi_all(i,k) = QQi_prog(i,k) + QQi1(i,k) + QQi2(i,k) + qiten_pwi1(i,k) + &
              qiten_pwi2(i,k) + A_qi_adj(i,k)
         QQ_all(i,k) = QQw_all(i,k) + QQi_all(i,k)
         QQnl_final(i,k) = QQnl_prog(i,k)
         QQni_final(i,k) = QQni_prog(i,k)
         QQn_final(i,k) = QQnl_final(i,k) + QQni_final(i,k)
         QQnl_all(i,k) = QQnl_prog(i,k) + QQnl1(i,k) + QQnl2(i,k) + nlten_pwi1(i,k) + &
              nlten_pwi2(i,k) + ACnl(i,k) + A_nl_adj(i,k)
         QQni_all(i,k) = QQni_prog(i,k) + QQni1(i,k) + QQni2(i,k) + niten_pwi1(i,k) + &
              niten_pwi2(i,k) + ACni(i,k) + A_ni_adj(i,k)
         QQn_all(i,k) = QQnl_all(i,k) + QQni_all(i,k)
         qme(i,k) = QQ_final(i,k)
         qvadj(i,k) = qvten_pwi1(i,k) + qvten_pwi2(i,k) + A_qv_adj(i,k)
         qladj(i,k) = qlten_pwi1(i,k) + qlten_pwi2(i,k) + A_ql_adj(i,k)
         qiadj(i,k) = qiten_pwi1(i,k) + qiten_pwi2(i,k) + A_qi_adj(i,k)
         qllim(i,k) = QQw1(i,k) + QQw2(i,k)
         qilim(i,k) = QQi1(i,k) + QQi2(i,k)
         s_tendout(i,k) = cpair*(T_star(i,k) - T0(i,k))/dt - cpair*(A_T(i,k)+C_T(i,k))
         qv_tendout(i,k) = (qv_star(i,k) - qv0(i,k))/dt - (A_qv(i,k)+C_qv(i,k))
         ql_tendout(i,k) = (ql_star(i,k) - ql0(i,k))/dt - (A_ql(i,k)+C_ql(i,k))
         qi_tendout(i,k) = (qi_star(i,k) - qi0(i,k))/dt - (A_qi(i,k)+C_qi(i,k))
         nl_tendout(i,k) = (nl_star(i,k) - nl0(i,k))/dt - (A_nl(i,k)+C_nl(i,k))
         ni_tendout(i,k) = (ni_star(i,k) - ni0(i,k))/dt - (A_ni(i,k)+C_ni(i,k))
         if (.not. do_cldice) then
            qi_tendout(i,k) = 0._r8
            ni_tendout(i,k) = 0._r8
         end if
         cld(i,k) = a_st_star(i,k) + a_cu0(i,k)
         T0(i,k) = T_star(i,k)
         qv0(i,k) = qv_star(i,k)
         ql0(i,k) = ql_star(i,k)
         qi0(i,k) = qi_star(i,k)
         nl0(i,k) = nl_star(i,k)
         ni0(i,k) = ni_star(i,k)
      end do
   end do

   end subroutine final_tendency_codon_wrap

!=======================================================================================================

   subroutine instratus_condensate( lchnk, ncol, k,                      &  
                                    p_in, T0_in, qv0_in, ql0_in, qi0_in, & 
                                    ni0_in,                              &
                                    a_dc_in, ql_dc_in, qi_dc_in,         &
                                    a_sc_in, ql_sc_in, qi_sc_in,         & 
                                    landfrac, snowh,                     &
                                    rhmini_in, rhminl_in, rhminl_adj_land_in, rhminh_in, &
                                    T_out, qv_out, ql_out, qi_out,       &
                                    al_st_out, ai_st_out, ql_st_out, qi_st_out )

   ! ------------------------------------------------------- !
   ! Diagnostically force in-stratus condensate to be        ! 
   ! in the range of 'qlst_min < qc_st < qlst_max'           !
   ! whenever stratus exists in the equilibrium state        !
   ! ------------------------------------------------------- !

   integer,  intent(in)  :: lchnk                ! Chunk identifier
   integer,  intent(in)  :: ncol                 ! Number of atmospheric columns
   integer,  intent(in)  :: k                    ! Layer index

   real(r8), intent(in)  :: p_in(pcols)          ! Pressure [Pa]
   real(r8), intent(in)  :: T0_in(pcols)         ! Temperature [K]
   real(r8), intent(in)  :: qv0_in(pcols)        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(in)  :: ql0_in(pcols)        ! Grid-mean LWC [kg/kg]
   real(r8), intent(in)  :: qi0_in(pcols)        ! Grid-mean IWC [kg/kg]
   real(r8), intent(in)  :: ni0_in(pcols)

   real(r8), intent(in)  :: a_dc_in(pcols)       ! Deep cumulus cloud fraction
   real(r8), intent(in)  :: ql_dc_in(pcols)      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_dc_in(pcols)      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)  :: a_sc_in(pcols)       ! Shallow cumulus cloud fraction
   real(r8), intent(in)  :: ql_sc_in(pcols)      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_sc_in(pcols)      ! In-shallow cumulus IWC [kg/kg]

   real(r8), intent(in)  :: landfrac(pcols)      ! Land fraction
   real(r8), intent(in)  :: snowh(pcols)         ! Snow depth (liquid water equivalent)

   real(r8), intent(in)  :: rhmini_in(pcols) 
   real(r8), intent(in)  :: rhminl_in(pcols)
   real(r8), intent(in)  :: rhminl_adj_land_in(pcols)
   real(r8), intent(in)  :: rhminh_in(pcols)     

   real(r8), intent(out) :: T_out(pcols)         ! Temperature [K]
   real(r8), intent(out) :: qv_out(pcols)        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(out) :: ql_out(pcols)        ! Grid-mean LWC [kg/kg]
   real(r8), intent(out) :: qi_out(pcols)        ! Grid-mean IWC [kg/kg]

   real(r8), intent(out) :: al_st_out(pcols)     ! Liquid stratus fraction
   real(r8), intent(out) :: ai_st_out(pcols)     ! Ice stratus fraction
   real(r8), intent(out) :: ql_st_out(pcols)     ! In-stratus LWC [kg/kg]
   real(r8), intent(out) :: qi_st_out(pcols)     ! In-stratus IWC [kg/kg]

   ! Local variables

   integer i                                     ! Column    index

   real(r8) p    
   real(r8) T0   
   real(r8) qv0    
   real(r8) ql0    
   real(r8) qi0    
   real(r8) a_dc   
   real(r8) ql_dc  
   real(r8) qi_dc  
   real(r8) a_sc   
   real(r8) ql_sc  
   real(r8) qi_sc  
   real(r8) esat0  
   real(r8) qsat0  
   real(r8) U0     
   real(r8) U0_nc  
   real(r8) G0_nc
   real(r8) al0_st_nc            
   real(r8) al0_st
   real(r8) ai0_st_nc            
   real(r8) ai0_st               
   real(r8) a0_st               
   real(r8) ql0_nc
   real(r8) qi0_nc
   real(r8) qc0_nc
   real(r8) ql0_st
   real(r8) qi0_st
   real(r8) qc0_st
   real(r8) T   
   real(r8) qv    
   real(r8) ql    
   real(r8) qi
   real(r8) ql_st
   real(r8) qi_st
   real(r8) es  
   real(r8) qs  
   real(r8) esat_in(pcols)  
   real(r8) qsat_in(pcols)  
   real(r8) U0_in(pcols)  
   real(r8) al0_st_nc_in(pcols)
   real(r8) ai0_st_nc_in(pcols)
   real(r8) G0_nc_in(pcols)
   integer  idxmod 
   real(r8) U
   real(r8) U_nc
   real(r8) al_st_nc
   real(r8) ai_st_nc
   real(r8) G_nc
   real(r8) a_st
   real(r8) al_st
   real(r8) ai_st
   real(r8) Tmin0
   real(r8) Tmax0
   real(r8) Tmin
   real(r8) Tmax
   integer caseid

   real(r8) rhmini
   real(r8) rhminl
   real(r8) rhminl_adj_land
   real(r8) rhminh

   ! ---------------- !
   ! Main Computation ! 
   ! ---------------- !

   call qsat_water(T0_in(1:ncol), p_in(1:ncol), &
        esat_in(1:ncol), qsat_in(1:ncol))
   U0_in(:ncol) = qv0_in(:ncol)/qsat_in(:ncol)
   if( CAMstfrac ) then
       call astG_RHU(U0_in(:),p_in(:),qv0_in(:),landfrac(:),snowh(:),al0_st_nc_in(:),G0_nc_in(:),ncol,&
                     rhminl_in(:), rhminl_adj_land_in(:), rhminh_in(:))
   else
       call astG_PDF(U0_in(:),p_in(:),qv0_in(:),landfrac(:),snowh(:),al0_st_nc_in(:),G0_nc_in(:),ncol,&
                     rhminl_in(:), rhminl_adj_land_in(:), rhminh_in(:))
   endif
   call aist_vector(qv0_in(:),T0_in(:),p_in(:),qi0_in(:),ni0_in(:),landfrac(:),snowh(:),ai0_st_nc_in(:),ncol,&
                    rhmaxi, rhmini_in(:), rhminl_in(:), rhminl_adj_land_in(:), rhminh_in(:))

   do i = 1, ncol

      ! ---------------------- !
      ! Define local variables !
      ! ---------------------- !

      p   = p_in(i)

      T0  = T0_in(i)
      qv0 = qv0_in(i)
      ql0 = ql0_in(i)
      qi0 = qi0_in(i)

      a_dc  = a_dc_in(i)
      ql_dc = ql_dc_in(i)
      qi_dc = qi_dc_in(i)

      a_sc  = a_sc_in(i)
      ql_sc = ql_sc_in(i)
      qi_sc = qi_sc_in(i)

      ql_dc = 0._r8
      qi_dc = 0._r8
      ql_sc = 0._r8
      qi_sc = 0._r8

      es  = esat_in(i) 
      qs  = qsat_in(i) 
 
      rhmini = rhmini_in(i)     
      rhminl = rhminl_in(i)     
      rhminl_adj_land = rhminl_adj_land_in(i)     
      rhminh = rhminh_in(i)     

      idxmod = 0
      caseid = -1

      ! ------------------------------------------------------------ !
      ! Force the grid-mean RH to be smaller than 1 if oversaturated !
      ! In order to be compatible with reduced 3x3 QQ, condensation  !
      ! should occur only into the liquid in gridmean_RH.            !
      ! ------------------------------------------------------------ !

      if( qv0 .gt. qs ) then
          call gridmean_RH( lchnk, i, k, p, T0, qv0, ql0, qi0,      &
                            a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, &
                            landfrac(i), snowh(i) )
          call qsat_water(T0, p, esat0, qsat0)
          U0      = (qv0/qsat0)
          U0_nc   =  U0 
          if( CAMstfrac ) then
              call astG_RHU_single(U0_nc, p, qv0, landfrac(i), snowh(i), al0_st_nc, G0_nc, &
                 rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
          else
              call astG_PDF_single(U0_nc, p, qv0, landfrac(i), snowh(i), al0_st_nc, G0_nc, &
                 rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
          endif
          call aist_single(qv0,T0,p,qi0,landfrac(i),snowh(i),ai0_st_nc,&
                           rhmaxi, rhmini, rhminl, rhminl_adj_land, rhminh)
          ai0_st  = (1._r8-a_dc-a_sc)*ai0_st_nc
          al0_st  = (1._r8-a_dc-a_sc)*al0_st_nc
          a0_st   = max(ai0_st,al0_st)         
          idxmod  = 1 
      else
          ai0_st  = (1._r8-a_dc-a_sc)*ai0_st_nc_in(i)
          al0_st  = (1._r8-a_dc-a_sc)*al0_st_nc_in(i)
      endif    
      a0_st   = max(ai0_st,al0_st)         

      ! ----------------------- ! 
      ! Handling of input state !
      ! ----------------------- !

      ql0_nc  = max(0._r8,ql0-a_dc*ql_dc-a_sc*ql_sc)
      qi0_nc  = max(0._r8,qi0-a_dc*qi_dc-a_sc*qi_sc)
      qc0_nc  = ql0_nc + qi0_nc 

      Tmin0 = T0 - (latvap/cpair)*ql0
      Tmax0 = T0 + ((latvap+latice)/cpair)*qv0

      ! ------------------------------------------------------------- !
      ! Do nothing and just exit if generalized in-stratus condensate !
      ! condition is satisfied. This includes the case I.             !
      ! For 4x4 liquid stratus, a0_st --> al0_st.                     ! 
      ! ------------------------------------------------------------- !
      if( ( ql0_nc .ge. qlst_min*al0_st ) .and. ( ql0_nc .le. qlst_max*al0_st ) ) then

          ! ------------------ !
          ! This is the case I !
          ! ------------------ ! 
             T = T0
             qv = qv0
             ql = ql0
             qi = qi0
             caseid = 0
             goto 10
      else
         ! ----------------------------- !
         ! This is case II : Dense Cloud !
         ! ----------------------------- !   
         if( al0_st .eq. 0._r8 .and. ql0_nc .gt. 0._r8 ) then
             ! ------------------------------------- !
             ! Compute hypothetical full evaporation !
             ! ------------------------------------- !
             T  = Tmin0
             qv = qv0 + ql0 
             call qsat_water(T, p, es, qs)
             U  = qv/qs
             U_nc = U  
             if( CAMstfrac ) then
                 call astG_RHU_single(U_nc, p, qv, landfrac(i), snowh(i), al_st_nc, G_nc, &
                    rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
             else
                 call astG_PDF_single(U_nc, p, qv, landfrac(i), snowh(i), al_st_nc, G_nc, &
                    rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
             endif
             al_st = (1._r8-a_dc-a_sc)*al_st_nc  
             caseid = 0

             if( al_st .eq. 0._r8 ) then
                 ql = 0._r8
                 qi = qi0
                 idxmod = 1
                 caseid = 1
                 goto 10
             else
                 ! ------------------------------------------- !
                 ! Evaporate until qc_st decreases to qlst_max !
                 ! ------------------------------------------- !
                 Tmin = Tmin0
                 Tmax = T0 
                 call instratus_core( lchnk, i, k, p,                              &
                                      T0, qv0, ql0, 0._r8,                         &
                                      a_dc, ql_dc, qi_dc,                          &
                                      a_sc, ql_sc, qi_sc, ai0_st,                  &
                                      qlst_max, Tmin, Tmax, landfrac(i), snowh(i), &
                                      rhminl, rhminl_adj_land, rhminh,             &
                                      T, qv, ql, qi )   
                 idxmod = 1
                 caseid = 2
                 goto 10
             endif
         ! ------------------------------ !
         ! This is case III : Empty Cloud !
         ! ------------------------------ !  
         elseif( al0_st .gt. 0._r8 .and. ql0_nc .eq. 0._r8 ) then
              ! ------------------------------------------ ! 
              ! Condense until qc_st increases to qlst_min !
              ! ------------------------------------------ !
              Tmin = Tmin0
              Tmax = Tmax0  
              call instratus_core( lchnk, i, k, p,                              &
                                   T0, qv0, ql0, 0._r8,                         &
                                   a_dc, ql_dc, qi_dc,                          &
                                   a_sc, ql_sc, qi_sc, ai0_st,                  &
                                   qlst_min, Tmin, Tmax, landfrac(i), snowh(i), &
                                   rhminl, rhminl_adj_land, rhminh,             &
                                   T, qv, ql, qi )   
              idxmod = 1 
              caseid = 3
              goto 10
         ! --------------- !
         ! This is case IV !
         ! --------------- !   
         elseif( al0_st .gt. 0._r8 .and. ql0_nc .gt. 0._r8 ) then

             if( ql0_nc .gt. qlst_max*al0_st ) then
                 ! --------------------------------------- !
                 ! Evaporate until qc_st drops to qlst_max !
                 ! --------------------------------------- !
                 Tmin = Tmin0
                 Tmax = Tmax0
                 call instratus_core( lchnk, i, k, p,                              &
                                      T0, qv0, ql0, 0._r8,                         &
                                      a_dc, ql_dc, qi_dc,                          &
                                      a_sc, ql_sc, qi_sc, ai0_st,                  &
                                      qlst_max, Tmin, Tmax, landfrac(i), snowh(i), &
                                      rhminl, rhminl_adj_land, rhminh,             &
                                      T, qv, ql, qi )   
                 idxmod = 1
                 caseid = 4
                 goto 10
             elseif( ql0_nc .lt. qlst_min*al0_st ) then
                 ! -------------------------------------------- !
                 ! Condensate until qc_st increases to qlst_min !
                 ! -------------------------------------------- !
                 Tmin = Tmin0
                 Tmax = Tmax0 
                 call instratus_core( lchnk, i, k, p,                              &
                                      T0, qv0, ql0, 0._r8,                         &
                                      a_dc, ql_dc, qi_dc,                          &
                                      a_sc, ql_sc, qi_sc, ai0_st,                  &
                                      qlst_min, Tmin, Tmax, landfrac(i), snowh(i), & 
                                      rhminl, rhminl_adj_land, rhminh,             &
                                      T, qv, ql, qi )   
                 idxmod = 1
                 caseid = 5
                 goto 10
             else
                 ! ------------------------------------------------ !
                 ! This case should not happen. Issue error message !
                 ! ------------------------------------------------ !
                 write(iulog,*) 'Impossible case1 in instratus_condensate' 
                 call endrun
             endif
         ! ------------------------------------------------ !                   
         ! This case should not happen. Issue error message !
         ! ------------------------------------------------ !    
         else
             write(iulog,*) 'Impossible case2 in instratus_condensate' 
             write(iulog,*)  al0_st, a_sc, a_dc
             write(iulog,*)  1000*ql0_nc, 1000*(ql0+qi0)
             call endrun
         endif
      endif

10 continue   

   ! -------------------------------------------------- !
   ! Force final energy-moisture conserving consistency !
   ! -------------------------------------------------- !

     qi = qi0

     if( idxmod .eq. 1 ) then
         call aist_single(qv,T,p,qi,landfrac(i),snowh(i),ai_st_nc,&
                          rhmaxi, rhmini, rhminl, rhminl_adj_land, rhminh)
         ai_st = (1._r8-a_dc-a_sc)*ai_st_nc
         call qsat_water(T, p, es, qs)
         U     = (qv/qs)
         U_nc  =  U
         if( CAMstfrac ) then
             call astG_RHU_single(U_nc, p, qv, landfrac(i), snowh(i), al_st_nc, G_nc, &
                rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
         else
             call astG_PDF_single(U_nc, p, qv, landfrac(i), snowh(i), al_st_nc, G_nc, &
                rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
         endif
         al_st = (1._r8-a_dc-a_sc)*al_st_nc
     else
         ai_st  = (1._r8-a_dc-a_sc)*ai0_st_nc_in(i)
         al_st  = (1._r8-a_dc-a_sc)*al0_st_nc_in(i)
     endif

     a_st  = max(ai_st,al_st)

     if( al_st .eq. 0._r8 ) then
         ql_st = 0._r8
     else
         ql_st = ql/al_st
         ql_st = min(qlst_max,max(qlst_min,ql_st)) ! PJR
     endif
     if( ai_st .eq. 0._r8 ) then
         qi_st = 0._r8
     else
         qi_st = qi/ai_st
     endif

     qi    = ai_st*qi_st
     ql    = al_st*ql_st

     T     = T0 - (latvap/cpair)*(ql0-ql) - ((latvap+latice)/cpair)*(qi0-qi)
     qv    = qv0 + ql0 - ql + qi0 - qi

   ! -------------- !
   ! Send to output !
   ! -------------- !

   T_out(i)  = T
   qv_out(i) = qv
   ql_out(i) = ql
   qi_out(i) = qi
   al_st_out(i) = al_st
   ai_st_out(i) = ai_st
   ql_st_out(i) = ql_st
   qi_st_out(i) = qi_st

   enddo 

   return
   end subroutine instratus_condensate

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine instratus_core( lchnk, icol, k, p,                      &
                              T0, qv0, ql0, qi0,                      &
                              a_dc, ql_dc, qi_dc,                     & 
                              a_sc, ql_sc, qi_sc, ai_st,              &
                              qcst_crit, Tmin, Tmax, landfrac, snowh, &
                              rhminl, rhminl_adj_land, rhminh,        &
                              T, qv, ql, qi )

   ! ------------------------------------------------------ !
   ! Subroutine to find saturation equilibrium state using  ! 
   ! a Newton iteration method, so that 'qc_st = qcst_crit' !
   ! is satisfied.                                          !
   ! ------------------------------------------------------ !

   integer,  intent(in)  :: lchnk      ! Chunk identifier
   integer,  intent(in)  :: icol       ! Number of atmospheric columns
   integer,  intent(in)  :: k          ! Layer index

   real(r8), intent(in)  :: p          ! Pressure [Pa]
   real(r8), intent(in)  :: T0         ! Temperature [K]
   real(r8), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
   real(r8), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]

   real(r8), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
   real(r8), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
   real(r8), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(r8), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)

   real(r8), intent(in)  :: Tmin       ! Minimum temperature system can have [K]
   real(r8), intent(in)  :: Tmax       ! Maximum temperature system can have [K]
   real(r8), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
   real(r8), intent(in)  :: landfrac   ! Land fraction
   real(r8), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)

   real(r8), intent(in)  :: rhminl
   real(r8), intent(in)  :: rhminl_adj_land
   real(r8), intent(in)  :: rhminh

   real(r8), intent(out) :: T          ! Temperature [K]
   real(r8), intent(out) :: qv         ! Grid-mean water vapor [kg/kg]
   real(r8), intent(out) :: ql         ! Grid-mean LWC [kg/kg]
   real(r8), intent(out) :: qi         ! Grid-mean IWC [kg/kg]

   ! Local variables

   integer i                           ! Iteration index

   real(r8) muQ0, muQ
   real(r8) ql_nc0, qi_nc0, qc_nc0, qc_nc    
   real(r8) fice0, fice    
   real(r8) ficeg0, ficeg   
   real(r8) esat0
   real(r8) qsat0
   real(r8) dqcncdt, dastdt, dUdt
   real(r8) alpha, beta
   real(r8) U, U_nc
   real(r8) al_st_nc, G_nc
   real(r8) al_st

   ! Variables for root-finding algorithm

   integer j                          
   real(r8)  x1, x2
   real(r8)  rtsafe
   real(r8)  df, dx, dxold, f, fh, fl, temp, xh, xl
   real(r8), parameter :: xacc = 1.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   ql_nc0 = max(0._r8,ql0-a_dc*ql_dc-a_sc*ql_sc)
   qi_nc0 = max(0._r8,qi0-a_dc*qi_dc-a_sc*qi_sc)
   qc_nc0 = max(0._r8,ql0+qi0-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
   fice0  = 0._r8
   ficeg0 = 0._r8
   muQ0   = 1._r8

   ! ------------ !
   ! Root finding !
   ! ------------ !

   x1 = Tmin
   x2 = Tmax
   call funcd_instratus( x1, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                         qcst_crit, landfrac, snowh,                    &
                         rhminl, rhminl_adj_land, rhminh,               &
                         fl, df, qc_nc, fice, al_st )
   call funcd_instratus( x2, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                         qcst_crit, landfrac, snowh,                    &
                         rhminl, rhminl_adj_land, rhminh,               &
                         fh, df, qc_nc, fice, al_st )
   if((fl > 0._r8 .and. fh > 0._r8) .or. (fl < 0._r8 .and. fh < 0._r8)) then
       call funcd_instratus( T0, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                             a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st, &
                             qcst_crit, landfrac, snowh,                    &
                             rhminl, rhminl_adj_land, rhminh,               &
                             fl, df, qc_nc, fice, al_st )
       rtsafe = T0 
       goto 10       
   endif
   if( fl == 0._r8) then
           rtsafe = x1
           goto 10
   elseif( fh == 0._r8) then
           rtsafe = x2
           goto 10
   elseif( fl < 0._r8) then
           xl = x1
           xh = x2
   else
           xh = x1
           xl = x2
   end if
   rtsafe = 0.5_r8*(x1+x2)
   dxold = abs(x2-x1)
   dx = dxold
   call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                         a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
                         qcst_crit, landfrac, snowh,                        &
                         rhminl, rhminl_adj_land, rhminh,                   &
                         f, df, qc_nc, fice, al_st )
   do j = 1, 20
      if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0._r8 .or. abs(2.0_r8*f) > abs(dxold*df) ) then
           dxold = dx
           dx = 0.5_r8*(xh-xl)
           rtsafe = xl + dx
           if(xl == rtsafe) goto 10
      else
           dxold = dx
           dx = f/df
           temp = rtsafe
           rtsafe = rtsafe - dx
           if (temp == rtsafe) goto 10
      end if
    ! if(abs(dx) < xacc) goto 10
      call funcd_instratus( rtsafe, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0, &
                            a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,     &
                            qcst_crit, landfrac, snowh,                        &
                            rhminl, rhminl_adj_land, rhminh,                   &
                            f, df, qc_nc, fice, al_st )
    ! Sep.21.2010. Sungsu modified to enhance convergence and guarantee 'qlst_min <  qlst < qlst_max'.
      if( qcst_crit < 0.5_r8 * ( qlst_min + qlst_max ) ) then
          if( ( qc_nc*(1._r8-fice) .gt.          qlst_min*al_st .and. &
                qc_nc*(1._r8-fice) .lt. 1.1_r8 * qlst_min*al_st ) ) goto 10
      else
          if( ( qc_nc*(1._r8-fice) .gt. 0.9_r8 * qlst_max*al_st .and. &
                qc_nc*(1._r8-fice) .lt.          qlst_max*al_st ) ) goto 10
      endif
      if(f < 0._r8) then
          xl = rtsafe
      else
          xh = rtsafe
      endif

   enddo

10 continue

   ! ------------------------------------------- !
   ! Final safety check before sending to output !
   ! ------------------------------------------- !

   qc_nc = max(0._r8,qc_nc)

   T  = rtsafe
   ql = qc_nc*(1._r8-fice) + a_dc*ql_dc + a_sc*ql_sc
   qi = qc_nc*fice + a_dc*qi_dc + a_sc*qi_sc
   qv = qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc))
   qv = max(qv,1.e-12_r8) 

   return
   end subroutine instratus_core

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine funcd_instratus( T, p, T0, qv0, ql0, qi0, fice0, muQ0, qc_nc0,   &
                               a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, ai_st,  &
                               qcst_crit, landfrac, snowh,                     &
                               rhminl, rhminl_adj_land, rhminh,                &
                               f, fg, qc_nc, fice, al_st ) 

   ! --------------------------------------------------- !
   ! Subroutine to find function value and gradient at T !
   ! --------------------------------------------------- !

   implicit none

   real(r8), intent(in)  :: T          ! Iteration temperature [K]

   real(r8), intent(in)  :: p          ! Pressure [Pa]
   real(r8), intent(in)  :: T0         ! Initial temperature [K]
   real(r8), intent(in)  :: qv0        ! Grid-mean water vapor [kg/kg]
   real(r8), intent(in)  :: ql0        ! Grid-mean LWC [kg/kg]
   real(r8), intent(in)  :: qi0        ! Grid-mean IWC [kg/kg]
   real(r8), intent(in)  :: fice0      ! 
   real(r8), intent(in)  :: muQ0       ! 
   real(r8), intent(in)  :: qc_nc0     ! 

   real(r8), intent(in)  :: a_dc       ! Deep cumulus cloud fraction
   real(r8), intent(in)  :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)  :: a_sc       ! Shallow cumulus cloud fraction
   real(r8), intent(in)  :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)  :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(r8), intent(in)  :: ai_st      ! Ice stratus fraction (fixed)

   real(r8), intent(in)  :: qcst_crit  ! Critical in-stratus condensate [kg/kg]
   real(r8), intent(in)  :: landfrac   ! Land fraction
   real(r8), intent(in)  :: snowh      ! Snow depth (liquid water equivalent)

   real(r8), intent(in)  :: rhminl
   real(r8), intent(in)  :: rhminl_adj_land
   real(r8), intent(in)  :: rhminh

   real(r8), intent(out) :: f          ! Value of minimization function at T
   real(r8), intent(out) :: fg         ! Gradient of minimization function 
   real(r8), intent(out) :: qc_nc      !
   real(r8), intent(out) :: al_st      !
   real(r8), intent(out) :: fice       !

   ! Local variables

   real(r8) es
   real(r8) qs
   real(r8) dqsdT
   real(r8) dqcncdt
   real(r8) alpha
   real(r8) beta
   real(r8) U
   real(r8) U_nc
   real(r8) al_st_nc
   real(r8) G_nc
   real(r8) dUdt
   real(r8) dalstdt
   real(r8) qv

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   call qsat_water(T, p, es, qs, dqsdt=dqsdT)

   fice    = fice0 
   qc_nc   = (cpair/latvap)*(T-T0)+muQ0*qc_nc0       
   dqcncdt = (cpair/latvap) 
   qv      = (qv0 + ql0 + qi0 - (qc_nc + a_dc*(ql_dc+qi_dc) + a_sc*(ql_sc+qi_sc)))
   alpha   = (1._r8/qs)
   beta    = (qv/qs**2._r8)*dqsdT 

   U      =  (qv/qs)
   U_nc   =   U
   if( CAMstfrac ) then
       call astG_RHU_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
          rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
   else
       call astG_PDF_single(U_nc, p, qv, landfrac, snowh, al_st_nc, G_nc, &
          rhminl_in=rhminl, rhminl_adj_land_in=rhminl_adj_land, rhminh_in=rhminh)
   endif
   al_st   =  (1._r8-a_dc-a_sc)*al_st_nc 
   dUdt    = -(alpha*dqcncdt+beta)
   dalstdt =  (1._r8/G_nc)*dUdt
   if( U_nc .eq. 1._r8 ) dalstdt = 0._r8

   f  = qc_nc   - qcst_crit*al_st
   fg = dqcncdt - qcst_crit*dalstdt

   return
   end subroutine funcd_instratus

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine gridmean_RH( lchnk, icol, k, p, T, qv, ql, qi,       &
                           a_dc, ql_dc, qi_dc, a_sc, ql_sc, qi_sc, &
                           landfrac, snowh )

   ! ------------------------------------------------------------- !
   ! Subroutine to force grid-mean RH = 1 when RH > 1              !
   ! This is condensation process similar to instratus_condensate. !
   ! During condensation, we assume 'fice' is maintained in this   !
   ! verison for MG not for RK.                                    !
   ! ------------------------------------------------------------- !

   integer,  intent(in)    :: lchnk      ! Chunk identifier
   integer,  intent(in)    :: icol       ! Number of atmospheric columns
   integer,  intent(in)    :: k          ! Layer index

   real(r8), intent(in)    :: p          ! Pressure [Pa]
   real(r8), intent(inout) :: T          ! Temperature [K]
   real(r8), intent(inout) :: qv         ! Grid-mean water vapor [kg/kg]
   real(r8), intent(inout) :: ql         ! Grid-mean LWC [kg/kg]
   real(r8), intent(inout) :: qi         ! Grid-mean IWC [kg/kg]

   real(r8), intent(in)    :: a_dc       ! Deep cumulus cloud fraction
   real(r8), intent(in)    :: ql_dc      ! In-deep cumulus LWC [kg/kg]
   real(r8), intent(in)    :: qi_dc      ! In-deep cumulus IWC [kg/kg]
   real(r8), intent(in)    :: a_sc       ! Shallow cumulus cloud fraction
   real(r8), intent(in)    :: ql_sc      ! In-shallow cumulus LWC [kg/kg]
   real(r8), intent(in)    :: qi_sc      ! In-shallow cumulus IWC [kg/kg]

   real(r8), intent(in)    :: landfrac   ! Land fraction
   real(r8), intent(in)    :: snowh      ! Snow depth (liquid water equivalent)

   ! Local variables

   integer m                             ! Iteration index

   real(r8)  ql_nc0, qi_nc0, qc_nc0
   real(r8)  Tscale
   real(r8)  Tc, qt, qc, dqcdt, qc_nc    
   real(r8)  es, qs, dqsdT
   real(r8)  al_st_nc, G_nc
   real(r8)  f, fg
   real(r8), parameter :: xacc = 1.e-3_r8

   ! ---------------- !
   ! Main computation !
   ! ---------------- !

   ql_nc0 = max(0._r8,ql-a_dc*ql_dc-a_sc*ql_sc)
   qi_nc0 = max(0._r8,qi-a_dc*qi_dc-a_sc*qi_sc)
   qc_nc0 = max(0._r8,ql+qi-a_dc*(ql_dc+qi_dc)-a_sc*(ql_sc+qi_sc))
   Tc    = T - (latvap/cpair)*ql
   qt    = qv + ql

   do m = 1, 20
      call qsat_water(T, p, es, qs, dqsdt=dqsdT)
      Tscale = latvap/cpair
      qc     = (T-Tc)/Tscale
      dqcdt  = 1._r8/Tscale
      f      = qs + qc - qt 
      fg     = dqsdT + dqcdt
      fg     = sign(1._r8,fg)*max(1.e-10_r8,abs(fg))
    ! Sungsu modified convergence criteria to speed up convergence and guarantee RH <= 1.
      if( qc .ge. 0._r8 .and. ( qt - qc ) .ge. 0.999_r8*qs .and. ( qt - qc ) .le. 1._r8*qs ) then
          goto 10
      endif
      T = T - f/fg
   enddo
 ! write(iulog,*) 'Convergence in gridmean_RH is not reached. RH = ', ( qt - qc ) / qs
10 continue

   call qsat_water(T, p, es, qs)
 ! Sungsu modified 'qv = qs' in consistent with the modified convergence criteria above.
   qv = min(qt,qs) ! Modified
   ql = qt - qv
   T  = Tc + (latvap/cpair)*ql

   return
   end subroutine gridmean_RH

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

   subroutine positive_moisture_select_impl()
   character(len=32) :: impl_name
   integer :: n, status

   if (positive_moisture_impl_selected) return
   call get_environment_variable('CLDWAT2M_POSITIVE_MOISTURE_IMPL', value=impl_name, length=n, status=status)
   use_native_positive_moisture_impl = .false.
   if (status == 0 .and. n > 0) then
      select case (adjustl(impl_name(:n)))
      case ('native', 'Native', 'NATIVE')
         use_native_positive_moisture_impl = .true.
      case ('codon', 'Codon', 'CODON')
         use_native_positive_moisture_impl = .false.
      case default
         use_native_positive_moisture_impl = .false.
      end select
   end if
   positive_moisture_impl_selected = .true.
   if (masterproc) then
      if (use_native_positive_moisture_impl) then
         write(iulog,*) 'cldwat2m_positive_moisture implementation = native'
      else
         write(iulog,*) 'cldwat2m_positive_moisture implementation = codon'
      end if
   end if
   end subroutine positive_moisture_select_impl

   subroutine positive_moisture_log_entered()
   if (positive_moisture_entered_logged) return
   positive_moisture_entered_logged = .true.
   if (masterproc) then
      write(iulog,*) 'cldwat2m_positive_moisture entered (macrophysics positive moisture limiter = codon)'
   end if
   end subroutine positive_moisture_log_entered

   subroutine positive_moisture( ncol, dt, qvmin, qlmin, qimin, dp, &
                                 qv,   ql, qi,    t,     qvten, &
                                 qlten,    qiten, tten,  do_cldice)

   ! ------------------------------------------------------------------------------- !
   ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
   ! force them to be larger than minimum value by (1) condensating water vapor      !
   ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
   ! layer. '2._r8' is multiplied to the minimum values for safety.                  !
   ! Update final state variables and tendencies associated with this correction.    !
   ! If any condensation happens, update (s,t) too.                                  !
   ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
   ! input tendencies.                                                               !
   ! Be careful the order of k : '1': top layer, 'pver' : near-surface layer         ! 
   ! ------------------------------------------------------------------------------- !

   implicit none
   integer,  intent(in)     :: ncol
   real(r8), intent(in)     :: dt
   real(r8), target, intent(in)     :: dp(pcols,pver), qvmin(pcols,pver), qlmin(pcols,pver), qimin(pcols,pver)
   real(r8), target, intent(inout)  :: qv(pcols,pver), ql(pcols,pver), qi(pcols,pver), t(pcols,pver)
   real(r8), target, intent(out)    :: qvten(pcols,pver), qlten(pcols,pver), qiten(pcols,pver), tten(pcols,pver)
   logical, intent(in)      :: do_cldice
   integer   i, k
   integer(c_int64_t) :: do_cldice_c
   real(r8)  dql, dqi, dqv, sum, aa, dum 

   interface
      subroutine cldwat2m_positive_moisture_codon(ncol_c, pcols_c, pver_c, top_lev_c, do_cldice_c, &
           dt_c, latvap_c, latice_c, cpair_c, dp_p, qvmin_p, qlmin_p, qimin_p, qv_p, ql_p, qi_p, &
           t_p, qvten_p, qlten_p, qiten_p, tten_p) bind(c, name="cldwat2m_positive_moisture_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, do_cldice_c
         real(c_double), value :: dt_c, latvap_c, latice_c, cpair_c
         type(c_ptr), value :: dp_p, qvmin_p, qlmin_p, qimin_p, qv_p, ql_p, qi_p, t_p
         type(c_ptr), value :: qvten_p, qlten_p, qiten_p, tten_p
      end subroutine cldwat2m_positive_moisture_codon
   end interface

   call positive_moisture_select_impl()
   if (.not. use_native_positive_moisture_impl) then
      call positive_moisture_log_entered()
      do_cldice_c = 0_c_int64_t
      if (do_cldice) do_cldice_c = 1_c_int64_t
      call cldwat2m_positive_moisture_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(top_lev, c_int64_t), do_cldice_c, dt, latvap, latice, cpair, &
           c_loc(dp(1,1)), c_loc(qvmin(1,1)), c_loc(qlmin(1,1)), c_loc(qimin(1,1)), &
           c_loc(qv(1,1)), c_loc(ql(1,1)), c_loc(qi(1,1)), c_loc(t(1,1)), &
           c_loc(qvten(1,1)), c_loc(qlten(1,1)), c_loc(qiten(1,1)), c_loc(tten(1,1)))
      return
   end if

   tten(:ncol,:pver)  = 0._r8
   qvten(:ncol,:pver) = 0._r8
   qlten(:ncol,:pver) = 0._r8
   qiten(:ncol,:pver) = 0._r8

   do i = 1, ncol
      do k = top_lev, pver
         if( qv(i,k) .lt. qvmin(i,k) .or. ql(i,k) .lt. qlmin(i,k) .or. qi(i,k) .lt. qimin(i,k) ) then
             goto 10
         endif
      enddo
      goto 11
   10 continue
      do k = top_lev, pver    ! From the top to the 1st (lowest) layer from the surface
         dql = max(0._r8,1._r8*qlmin(i,k)-ql(i,k))

         if (do_cldice) then
         dqi = max(0._r8,1._r8*qimin(i,k)-qi(i,k))
         else
           dqi = 0._r8
         end if

         qlten(i,k) = qlten(i,k) +  dql/dt
         qiten(i,k) = qiten(i,k) +  dqi/dt
         qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
         tten(i,k)  = tten(i,k)  + (latvap/cpair)*(dql/dt) + ((latvap+latice)/cpair)*(dqi/dt)
         ql(i,k)    = ql(i,k) + dql
         qi(i,k)    = qi(i,k) + dqi
         qv(i,k)    = qv(i,k) - dql - dqi
         t(i,k)     = t(i,k)  + (latvap * dql + (latvap+latice) * dqi)/cpair
         dqv        = max(0._r8,1._r8*qvmin(i,k)-qv(i,k))
         qvten(i,k) = qvten(i,k) + dqv/dt
         qv(i,k)    = qv(i,k)    + dqv
         if( k .ne. pver ) then 
             qv(i,k+1)    = qv(i,k+1)    - dqv*dp(i,k)/dp(i,k+1)
             qvten(i,k+1) = qvten(i,k+1) - dqv*dp(i,k)/dp(i,k+1)/dt
         endif
         qv(i,k) = max(qv(i,k),qvmin(i,k))
         ql(i,k) = max(ql(i,k),qlmin(i,k))
         qi(i,k) = max(qi(i,k),qimin(i,k))
      end do
      ! Extra moisture used to satisfy 'qv(i,pver)=qvmin' is proportionally 
      ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
      ! preserves column moisture. 
      if( dqv .gt. 1.e-20_r8 ) then
          sum = 0._r8
          do k = top_lev, pver
             if( qv(i,k) .gt. 2._r8*qvmin(i,k) ) sum = sum + qv(i,k)*dp(i,k)
          enddo
          aa = dqv*dp(i,pver)/max(1.e-20_r8,sum)
          if( aa .lt. 0.5_r8 ) then
              do k = top_lev, pver
                 if( qv(i,k) .gt. 2._r8*qvmin(i,k) ) then
                     dum        = aa*qv(i,k)
                     qv(i,k)    = qv(i,k) - dum
                     qvten(i,k) = qvten(i,k) - dum/dt
                 endif
              enddo 
          else 
              write(iulog,*) 'Full positive_moisture is impossible in Park Macro'
          endif
      endif 
11 continue
   enddo
   return

   end subroutine positive_moisture

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      real(r8) a(np,np),b(np,mp)
      real(r8) aa(np,np),bb(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,ii,jj,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      real(r8) big,dum,pivinv

      aa(:,:) = a(:,:)
      bb(:,:) = b(:,:)

      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0._r8
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(iulog,*) 'singular matrix in gaussj 1'
                do ii = 1, np
                do jj = 1, np
                   write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
                end do
                end do   
                call endrun
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0._r8) then
            write(iulog,*) 'singular matrix in gaussj 2'
            do ii = 1, np
            do jj = 1, np
               write(iulog,*) ii, jj, aa(ii,jj), bb(ii,1)
            end do
            end do   
            call endrun
        endif 
        pivinv=1._r8/a(icol,icol)
        a(icol,icol)=1._r8
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0._r8
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue

      return
      end subroutine gaussj

   ! ----------------- !
   ! End of subroutine !
   ! ----------------- !

end module cldwat2m_macro
