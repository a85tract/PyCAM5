   module convect_shallow

   !----------------------------------------------- !
   ! Purpose:                                       !
   !                                                !
   ! CAM interface to the shallow convection scheme !
   !                                                !
   ! Author: D.B. Coleman                           !
   !         Sungsu Park. Jan. 2010.                !
   !                                                !
   !----------------------------------------------- !

   use shr_kind_mod,      only : r8=>shr_kind_r8
   use physconst,         only : cpair, zvir
   use ppgrid,            only : pver, pcols, pverp
   use zm_conv,           only : zm_conv_evap
   use cam_history,       only : outfld, addfld, phys_decomp
   use cam_logfile,       only : iulog
   use phys_control,      only : phys_getopts
   use iso_c_binding,     only : c_int64_t

   implicit none
   private                 
   save

   public :: &
             convect_shallow_register,       & ! Register fields in physics buffer
             convect_shallow_init,           & ! Initialize shallow module
             convect_shallow_init_cnst,	     & ! 
             convect_shallow_implements_cnst,&
             convect_shallow_tend,           & ! Return tendencies
             convect_shallow_use_shfrc	       ! 

   ! The following namelist variable controls which shallow convection package is used.
   !        'Hack'   = Hack shallow convection (default)
   !        'UW'     = UW shallow convection by Sungsu Park and Christopher S. Bretherton
   !        'UNICON' = General Convection Model by Sungsu Park  
   !        'off'    = No shallow convection

   character(len=16) :: shallow_scheme      ! Default set in phys_control.F90, use namelist to change
   character(len=16) :: microp_scheme       ! Microphysics scheme
   logical           :: history_amwg        ! output the variables used by the AMWG diag package
   logical           :: history_budget      ! Output tendencies and state variables for CAM4 T, qv, ql, qi
   integer           :: history_budget_histfile_num ! output history file number for budget fields

   ! Physics buffer indices 
   integer    ::     icwmrsh_idx    = 0  
   integer    ::      rprdsh_idx    = 0 
   integer    ::     rprdtot_idx    = 0 
   integer    ::      cldtop_idx    = 0 
   integer    ::      cldbot_idx    = 0 
   integer    ::        cush_idx    = 0 
   integer    :: nevapr_shcu_idx    = 0
   integer    ::       shfrc_idx    = 0 
   integer    ::         cld_idx    = 0 
   integer    ::      concld_idx    = 0
   integer    ::      rprddp_idx    = 0
   integer    ::         tke_idx    = 0

   integer    ::       qpert_idx    = 0
   integer    ::       pblh_idx     = 0
   integer    ::    prec_sh_idx     = 0
   integer    ::    snow_sh_idx     = 0

   integer    ::  ttend_sh_idx      = 0
   logical    :: use_native_impl    = .false.
   logical    :: impl_selected      = .false.
   integer    :: codon_scheme_code  = 0
   logical    :: codon_scheme_selected = .false.
   logical    :: convect_shallow_diag_shell_logged = .false.
   logical    :: convect_shallow_init_direct_logged = .false.
   logical    :: convect_shallow_init_shell_logged = .false.
   logical    :: convect_shallow_ptend_lq_mask_shell_logged = .false.
   logical    :: convect_shallow_uw_post_shell_logged = .false.
   logical    :: convect_shallow_wtrc_precip_shell_logged = .false.
   logical    :: convect_shallow_tend_logged = .false.
   logical    :: use_native_init_impl = .false.
   logical    :: init_impl_selected = .false.
   logical    :: init_mw_ratio_logged = .false.
   logical    :: convect_shallow_use_shfrc_logged = .false.
   logical    :: convect_shallow_register_logged = .false.

   integer :: & ! field index in physics buffer
      sh_flxprc_idx, &
      sh_flxsnw_idx, &
      sh_cldliq_idx, &
      sh_cldice_idx

   interface
	     function convect_shallow_use_shfrc_codon(scheme_len_c, scheme_ascii_p) result(out_c) &
	          bind(c, name="convect_shallow_use_shfrc_codon")
	        use iso_c_binding, only: c_int64_t, c_ptr
	        integer(c_int64_t), value :: scheme_len_c
	        type(c_ptr), value :: scheme_ascii_p
	        integer(c_int64_t) :: out_c
	     end function convect_shallow_use_shfrc_codon
	     function convect_shallow_register_decision_codon(scheme_len_c, scheme_ascii_p, use_gw_convect_sh_c) &
	          result(mask_c) bind(c, name="convect_shallow_register_decision_codon")
	        use iso_c_binding, only: c_int64_t, c_ptr
	        integer(c_int64_t), value :: scheme_len_c, use_gw_convect_sh_c
	        type(c_ptr), value :: scheme_ascii_p
	        integer(c_int64_t) :: mask_c
	     end function convect_shallow_register_decision_codon
	     subroutine convect_shallow_select_scheme_codon(scheme_len_c, scheme_ascii_p, scheme_code_p, status_p) &
	          bind(c, name="convect_shallow_select_scheme_codon")
	        use iso_c_binding, only: c_int64_t, c_ptr
	        integer(c_int64_t), value :: scheme_len_c
	        type(c_ptr), value :: scheme_ascii_p, scheme_code_p, status_p
	     end subroutine convect_shallow_select_scheme_codon
	  end interface

   contains

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine convect_shallow_register

  !-------------------------------------------------- !
  ! Purpose : Register fields with the physics buffer !
  !-------------------------------------------------- !

	  use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls
	  use phys_control, only: use_gw_convect_sh
	  use unicon_cam,     only: unicon_cam_register
	  use spmd_utils,     only: masterproc
	  use iso_c_binding,  only: c_int64_t, c_loc

	  integer(c_int64_t), target :: scheme_ascii(len(shallow_scheme))
	  integer(c_int64_t) :: register_mask
	  integer :: i

	  call phys_getopts( shallow_scheme_out = shallow_scheme, microp_scheme_out = microp_scheme)

	  do i = 1, len(shallow_scheme)
	     scheme_ascii(i) = int(iachar(shallow_scheme(i:i)), c_int64_t)
	  end do
	  register_mask = convect_shallow_register_decision_codon( &
	       int(len(shallow_scheme), c_int64_t), c_loc(scheme_ascii(1)), &
	       merge(1_c_int64_t, 0_c_int64_t, use_gw_convect_sh))

	  call pbuf_add_field('ICWMRSH',    'physpkg' ,dtype_r8,(/pcols,pver/),       icwmrsh_idx )
  call pbuf_add_field('RPRDSH',     'physpkg' ,dtype_r8,(/pcols,pver/),       rprdsh_idx )
  call pbuf_add_field('RPRDTOT',    'physpkg' ,dtype_r8,(/pcols,pver/),       rprdtot_idx )
  call pbuf_add_field('CLDTOP',     'physpkg' ,dtype_r8,(/pcols,1/),          cldtop_idx )
  call pbuf_add_field('CLDBOT',     'physpkg' ,dtype_r8,(/pcols,1/),          cldbot_idx )
  call pbuf_add_field('cush',       'global'  ,dtype_r8,(/pcols,dyn_time_lvls/), cush_idx ) 	
  call pbuf_add_field('NEVAPR_SHCU','physpkg' ,dtype_r8,(/pcols,pver/),       nevapr_shcu_idx )
  call pbuf_add_field('PREC_SH',    'physpkg' ,dtype_r8,(/pcols/),            prec_sh_idx )
  call pbuf_add_field('SNOW_SH',    'physpkg' ,dtype_r8,(/pcols/),            snow_sh_idx )

	  if (iand(int(register_mask), 1) /= 0) then
	     call pbuf_add_field('shfrc', 'physpkg', dtype_r8, (/pcols,pver/), shfrc_idx)
	  end if

! shallow interface gbm flux_convective_cloud_rain+snow (kg/m2/s)
  call pbuf_add_field('SH_FLXPRC','physpkg',dtype_r8,(/pcols,pverp/),sh_flxprc_idx)  

! shallow interface gbm flux_convective_cloud_snow (kg/m2/s)
  call pbuf_add_field('SH_FLXSNW','physpkg',dtype_r8,(/pcols,pverp/),sh_flxsnw_idx)  

! shallow gbm cloud liquid water (kg/kg)
  call pbuf_add_field('SH_CLDLIQ','physpkg',dtype_r8,(/pcols,pver/),sh_cldliq_idx)  

! shallow gbm cloud ice water (kg/kg)
  call pbuf_add_field('SH_CLDICE','physpkg',dtype_r8,(/pcols,pver/),sh_cldice_idx)  

  ! If gravity waves from shallow convection are on, output this field.
	  if (iand(int(register_mask), 2) /= 0) then
	     call pbuf_add_field('TTEND_SH','physpkg',dtype_r8,(/pcols,pver/),ttend_sh_idx)
	  end if

	  if (iand(int(register_mask), 4) /= 0) then
	     call unicon_cam_register()
	  end if

	  if (masterproc .and. .not. convect_shallow_register_logged) then
	     convect_shallow_register_logged = .true.
	     write(iulog,'(A)') &
	          'convect_shallow_register direct = codon; pbuf_add_field/unicon_cam_register native CAM API islands'
	     call convect_shallow_append_proof( &
	          'convect_shallow_register direct = codon; pbuf_add_field/unicon_cam_register native CAM API islands')
	     call flush(iulog)
	  end if

	  end subroutine convect_shallow_register

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !


  subroutine convect_shallow_init(pref_edge, pbuf2d)

  !------------------------------------------------------------------------------- !
  ! Purpose : Declare output fields, and initialize variables needed by convection !
  !------------------------------------------------------------------------------- !

  use iso_c_binding,    only : c_double, c_int64_t, c_loc
  use cam_history,       only : addfld, add_default, phys_decomp
  use ppgrid,            only : pcols, pver
  use hk_conv,           only : mfinti
  use uwshcu,            only : init_uwshcu
  use unicon_cam,        only : unicon_cam_init
  use physconst,         only : rair, gravit, latvap, rhoh2o, zvir, &
                                cappa, latice, mwdry, mwh2o
  use pmgrid,            only : plev, plevp
  use spmd_utils,        only : masterproc
  use cam_abortutils,    only : endrun
  use phys_control,      only : cam_physpkg_is
  
  use physics_buffer,    only : pbuf_get_index, physics_buffer_desc, pbuf_set_field
  use time_manager,      only : is_first_step

  real(r8),                  intent(in) :: pref_edge(plevp)  ! Reference pressures at interfaces
  type(physics_buffer_desc), pointer    :: pbuf2d(:,:)

  integer limcnv                                   ! Top interface level limit for convection
  integer i, k
  character(len=16)          :: eddy_scheme
  real(r8)                   :: mwh2o_mwdry_ratio
  integer                    :: scheme_action
  integer(c_int64_t)         :: init_action_c
  integer(c_int64_t), target :: scheme_ascii(len(shallow_scheme))

  interface
     function convect_shallow_init_action_codon(scheme_len_c, scheme_ascii_p) result(action_c) &
          bind(c, name="convect_shallow_init_action_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: scheme_len_c
       type(c_ptr), value :: scheme_ascii_p
       integer(c_int64_t) :: action_c
     end function convect_shallow_init_action_codon
     function convect_shallow_init_mw_ratio_codon(mwh2o_c, mwdry_c) result(ratio_c) &
          bind(c, name="convect_shallow_init_mw_ratio_codon")
        use iso_c_binding, only: c_double
        real(c_double), value :: mwh2o_c, mwdry_c
        real(c_double) :: ratio_c
     end function convect_shallow_init_mw_ratio_codon
  end interface

  call convect_shallow_init_select_impl()
  if (use_native_init_impl) then
     select case (shallow_scheme)
     case('off')
        scheme_action = 1
     case('Hack')
        scheme_action = 2
     case('UW')
        scheme_action = 3
     case('UNICON')
        scheme_action = 4
     case default
        scheme_action = 0
     end select
  else
     do i = 1, len(shallow_scheme)
        scheme_ascii(i) = int(iachar(shallow_scheme(i:i)), c_int64_t)
     end do
     init_action_c = convect_shallow_init_action_codon(int(len(shallow_scheme), c_int64_t), c_loc(scheme_ascii(1)))
     scheme_action = int(init_action_c)
     call convect_shallow_log_init_direct()
  end if
    
  ! ------------------------------------------------- !
  ! Variables for detailed abalysis of UW-ShCu scheme !
  ! ------------------------------------------------- !

  call addfld( 'qt_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qt_preCU'                                         ,  phys_decomp )
  call addfld( 'sl_pre_Cu    ', 'J/kg'    ,  pver ,  'I' , 'sl_preCU'                                         ,  phys_decomp )
  call addfld( 'slv_pre_Cu   ', 'J/kg'    ,  pver ,  'I' , 'slv_preCU'                                        ,  phys_decomp )
  call addfld( 'u_pre_Cu     ', 'm/s'     ,  pver ,  'I' , 'u_preCU'                                          ,  phys_decomp )
  call addfld( 'v_pre_Cu     ', 'm/s'     ,  pver ,  'I' , 'v_preCU'                                          ,  phys_decomp )
  call addfld( 'qv_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qv_preCU'                                         ,  phys_decomp )
  call addfld( 'ql_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'ql_preCU'                                         ,  phys_decomp )
  call addfld( 'qi_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qi_preCU'                                         ,  phys_decomp )
  call addfld( 't_pre_Cu     ', 'K'       ,  pver ,  'I' , 't_preCU'                                          ,  phys_decomp )
  call addfld( 'rh_pre_Cu    ', '%'       ,  pver ,  'I' , 'rh_preCU'                                         ,  phys_decomp )

  call addfld( 'qt_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qt_afterCU'                                       ,  phys_decomp )
  call addfld( 'sl_aft_Cu    ', 'J/kg'    ,  pver ,  'I' , 'sl_afterCU'                                       ,  phys_decomp )
  call addfld( 'slv_aft_Cu   ', 'J/kg'    ,  pver ,  'I' , 'slv_afterCU'                                      ,  phys_decomp )
  call addfld( 'u_aft_Cu     ', 'm/s'     ,  pver ,  'I' , 'u_afterCU'                                        ,  phys_decomp )
  call addfld( 'v_aft_Cu     ', 'm/s'     ,  pver ,  'I' , 'v_afterCU'                                        ,  phys_decomp )
  call addfld( 'qv_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qv_afterCU'                                       ,  phys_decomp )
  call addfld( 'ql_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'ql_afterCU'                                       ,  phys_decomp )
  call addfld( 'qi_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qi_afterCU'                                       ,  phys_decomp )
  call addfld( 't_aft_Cu     ', 'K'       ,  pver ,  'I' , 't_afterCU'                                        ,  phys_decomp )
  call addfld( 'rh_aft_Cu    ', '%'       ,  pver ,  'I' , 'rh_afterCU'                                       ,  phys_decomp )

  call addfld( 'tten_Cu      ', 'K/s'     ,  pver ,  'I' , 'Temperature tendency by cumulus convection'       ,  phys_decomp )
  call addfld( 'rhten_Cu     ', '%/s'     ,  pver ,  'I' , 'RH tendency by cumumus convection'                ,  phys_decomp )

  ! ------------------------------------------- !
  ! Common Output for Shallow Convection Scheme !
  ! ------------------------------------------- !

  call addfld( 'CMFDT   '     , 'K/s     ',  pver ,  'A' , &
       'T tendency - shallow convection'                           ,  phys_decomp )
  call addfld( 'CMFDQ   '     , 'kg/kg/s ',  pver ,  'A' , &
       'QV tendency - shallow convection'                          ,  phys_decomp )
  call addfld( 'CMFDLIQ '     , 'kg/kg/s ',  pver ,  'A' , &
       'Cloud liq tendency - shallow convection'                   ,  phys_decomp )
  call addfld( 'CMFDICE '     , 'kg/kg/s ',  pver ,  'A' , &
       'Cloud ice tendency - shallow convection'                   ,  phys_decomp )
  call addfld( 'CMFDQR  '     , 'kg/kg/s ',  pver ,  'A' , &
       'Q tendency - shallow convection rainout'                   ,  phys_decomp )
  call addfld( 'EVAPTCM '     , 'K/s     ',  pver ,  'A' , &
       'T tendency - Evaporation/snow prod from Hack convection'   ,  phys_decomp )
  call addfld( 'FZSNTCM '     , 'K/s     ',  pver ,  'A' , &
       'T tendency - Rain to snow conversion from Hack convection' ,  phys_decomp )
  call addfld( 'EVSNTCM '     , 'K/s     ',  pver ,  'A' , &
       'T tendency - Snow to rain prod from Hack convection'       ,  phys_decomp )
  call addfld( 'EVAPQCM '     , 'kg/kg/s ',  pver ,  'A' , &
       'Q tendency - Evaporation from Hack convection'             ,  phys_decomp )
  call addfld( 'QC      '     , 'kg/kg/s ',  pver ,  'A' , &
       'Q tendency - shallow convection LW export'                 ,  phys_decomp )
  call addfld( 'PRECSH  '     , 'm/s     ',  1,      'A' , &
       'Shallow Convection precipitation rate'                     ,  phys_decomp )
  call addfld( 'CMFMC   '     , 'kg/m2/s ',  pverp,  'A' , &
       'Moist convection (deep+shallow) mass flux'                 ,  phys_decomp )
  call addfld( 'CMFSL   '     , 'W/m2    ',  pverp,  'A' , &
       'Moist shallow convection liquid water static energy flux'  ,  phys_decomp )
  call addfld( 'CMFLQ   '     , 'W/m2    ',  pverp,  'A' , &
       'Moist shallow convection total water flux'                 ,  phys_decomp )
  call addfld( 'CIN     '     , 'J/kg    ',  1    ,  'A' , &
       'Convective inhibition'                                     ,  phys_decomp )
  call addfld( 'CBMF    '     , 'kg/m2/s ',  1    ,  'A' , &
       'Cloud base mass flux'                                      ,  phys_decomp )
  call addfld( 'CLDTOP  '     , '1       ',  1    ,  'I' , &
       'Vertical index of cloud top'                               ,  phys_decomp )
  call addfld( 'CLDBOT  '     , '1       ',  1    ,  'I' , &
       'Vertical index of cloud base'                              ,  phys_decomp )
  call addfld( 'PCLDTOP '     , '1       ',  1    ,  'A' , &
       'Pressure of cloud top'                                     ,  phys_decomp )
  call addfld( 'PCLDBOT '     , '1       ',  1    ,  'A' , &
       'Pressure of cloud base'                                    ,  phys_decomp )

  call addfld( 'FREQSH '      , 'fraction',  1    ,  'A' , &
       'Fractional occurance of shallow convection'                ,  phys_decomp )
                                                                                                                    
  call addfld( 'HKFLXPRC'     , 'kg/m2/s ',  pverp,  'A' , &
       'Flux of precipitation from HK convection'                  ,  phys_decomp )
  call addfld( 'HKFLXSNW'     , 'kg/m2/s ',  pverp,  'A' , &
       'Flux of snow from HK convection'                           ,  phys_decomp )
  call addfld( 'HKNTPRPD'     , 'kg/kg/s ',  pver ,  'A' , &
       'Net precipitation production from HK convection'           ,  phys_decomp )
  call addfld( 'HKNTSNPD'     , 'kg/kg/s ',  pver ,  'A' , &
       'Net snow production from HK convection'                    ,  phys_decomp )
  call addfld( 'HKEIHEAT'     , 'W/kg'    ,  pver ,  'A' , &
       'Heating by ice and evaporation in HK convection'           ,  phys_decomp )

  call addfld ('ICWMRSH  '    , 'kg/kg   ',  pver,   'A' , &
       'Shallow Convection in-cloud water mixing ratio '           ,  phys_decomp )

  if( scheme_action == 3 ) then
     call addfld( 'UWFLXPRC'     , 'kg/m2/s ',  pverp,  'A' , &
          'Flux of precipitation from UW shallow convection'          ,  phys_decomp )
     call addfld( 'UWFLXSNW'     , 'kg/m2/s ',  pverp,  'A' , &
          'Flux of snow from UW shallow convection'                   ,  phys_decomp )
  end if



  call phys_getopts( eddy_scheme_out = eddy_scheme      , &
                     history_amwg_out = history_amwg    , &
                     history_budget_out = history_budget, &
                     history_budget_histfile_num_out = history_budget_histfile_num)


  if( history_budget ) then
      call add_default( 'CMFDLIQ  ', history_budget_histfile_num, ' ' )
      call add_default( 'CMFDICE  ', history_budget_histfile_num, ' ' )
      call add_default( 'CMFDT   ', history_budget_histfile_num, ' ' )
      call add_default( 'CMFDQ   ', history_budget_histfile_num, ' ' )
      if( cam_physpkg_is('cam3') .or. cam_physpkg_is('cam4') ) then
         call add_default( 'EVAPQCM  ', history_budget_histfile_num, ' ' )
         call add_default( 'EVAPTCM  ', history_budget_histfile_num, ' ' )
      end if
  end if
  pblh_idx  = pbuf_get_index('pblh')


  select case (scheme_action)

  case(1)  ! None

     if( masterproc ) write(iulog,*) 'convect_shallow_init: shallow convection OFF'
     continue

  case(2) ! Hack scheme

     qpert_idx = pbuf_get_index('qpert')

     if( masterproc ) write(iulog,*) 'convect_shallow_init: Hack shallow convection'
   ! Limit shallow convection to regions below 40 mb
   ! Note this calculation is repeated in the deep convection interface
     if( pref_edge(1) >= 4.e3_r8 ) then
         limcnv = 1
     else
         do k = 1, plev
            if( pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8 ) then
                limcnv = k
                goto 10
            end if
         end do
         limcnv = plevp
     end if
10   continue

     if( masterproc ) then
         write(iulog,*) 'MFINTI: Convection will be capped at intfc ', limcnv, ' which is ', pref_edge(limcnv), ' pascals'
     end if
     
     call mfinti( rair, cpair, gravit, latvap, rhoh2o, limcnv) ! Get args from inti.F90

  case(3) ! Park and Bretherton shallow convection scheme

     if( masterproc ) write(iulog,*) 'convect_shallow_init: UW shallow convection scheme (McCaa)'
     if( eddy_scheme .ne. 'diag_TKE' ) then
         write(iulog,*) 'ERROR: shallow convection scheme ', shallow_scheme, ' is incompatible with eddy scheme ', eddy_scheme
         call endrun( 'convect_shallow_init: shallow_scheme and eddy_scheme are incompatible' )
     endif
     if (use_native_init_impl) then
        mwh2o_mwdry_ratio = mwh2o/mwdry
     else
        call convect_shallow_log_init_mw_ratio_entered()
        mwh2o_mwdry_ratio = real(convect_shallow_init_mw_ratio_codon(real(mwh2o, c_double), &
             real(mwdry, c_double)), r8)
     end if
     call init_uwshcu( r8, latvap, cpair, latice, zvir, rair, gravit, mwh2o_mwdry_ratio )

     tke_idx = pbuf_get_index('tke')

  case(4) ! Sungsu Park's General Convection Model

     if ( masterproc ) write(iulog,*) 'convect_shallow_init: General Convection Model by Sungsu Park'
     if ( eddy_scheme .ne. 'diag_TKE' ) then
          write(iulog,*)  eddy_scheme
          write(iulog,*) 'ERROR: shallow convection scheme ',shallow_scheme,' is incompatible with eddy scheme ', eddy_scheme
          call endrun( 'convect_shallow_init: shallow_scheme and eddy_scheme are incompatible' )
     endif
     call unicon_cam_init(pbuf2d)

  end select

  cld_idx      = pbuf_get_index('CLD')
  concld_idx   = pbuf_get_index('CONCLD')
  rprddp_idx   = pbuf_get_index('RPRDDP')

  end subroutine convect_shallow_init

!==================================================================================================

function convect_shallow_implements_cnst(name)

   ! Return true if specified constituent is implemented by a shallow convetion package

   use unicon_cam, only: unicon_implements_cnst

   character(len=*), intent(in) :: name          ! constituent name
   logical :: convect_shallow_implements_cnst    ! return value
   
   integer :: m
   !-----------------------------------------------------------------------

   select case (shallow_scheme)

   case('UNICON')
      convect_shallow_implements_cnst = unicon_implements_cnst(name)

   case default
      convect_shallow_implements_cnst = .false.

   end select

end function convect_shallow_implements_cnst

!==================================================================================================

subroutine convect_shallow_init_cnst(name, q, gcid)

  ! Initialize constituents if they are not read from the initial file

   use unicon_cam, only: unicon_init_cnst

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------

   select case (shallow_scheme)

   case('UNICON')
      call unicon_init_cnst(name, q, gcid)

   case default

   end select

end subroutine convect_shallow_init_cnst

!==================================================================================================

  function convect_shallow_use_shfrc()
  !-------------------------------------------------------------- !
  ! Return true if cloud fraction should use shallow convection   !
  !          calculated convective clouds.                        !
  !-------------------------------------------------------------- !
     use iso_c_binding, only : c_int64_t, c_loc
     use spmd_utils,    only : masterproc

     implicit none
     logical :: convect_shallow_use_shfrc     ! Return value
     integer(c_int64_t) :: out_c
     integer(c_int64_t), target :: scheme_ascii(len(shallow_scheme))
     integer :: i

     do i = 1, len(shallow_scheme)
        scheme_ascii(i) = int(iachar(shallow_scheme(i:i)), c_int64_t)
     end do

     out_c = convect_shallow_use_shfrc_codon(int(len(shallow_scheme), c_int64_t), c_loc(scheme_ascii(1)))
     convect_shallow_use_shfrc = out_c /= 0_c_int64_t
     if (.not. convect_shallow_use_shfrc_logged) then
        convect_shallow_use_shfrc_logged = .true.
        if (masterproc) then
           write(iulog,'(A)') 'convect_shallow_use_shfrc direct = codon'
           call convect_shallow_append_proof('convect_shallow_use_shfrc direct = codon')
           call flush(iulog)
        end if
     end if

     return

  end function convect_shallow_use_shfrc

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine convect_shallow_tend( ztodt  , cmfmc   , cmfmc2   , &
                                   qc     , qc2     , rliq     , rliq2    , & 
                                   state  , ptend_all, pbuf, sgh30, cam_in, wtdlf)

   use physics_buffer,  only : physics_buffer_desc, pbuf_get_field, pbuf_set_field, &
                               pbuf_old_tim_idx, pbuf_get_index
   use cam_history,     only : outfld
   use physics_types,   only : physics_state, physics_ptend
   use physics_types,   only : physics_ptend_init, physics_update
   use physics_types,   only : physics_state_copy, physics_state_dealloc
   use physics_types,   only : physics_ptend_dealloc
   use physics_types,   only : physics_ptend_sum
   use camsrfexch,      only : cam_in_t
   
   use constituents,    only : pcnst, cnst_get_ind, cnst_get_type_byind, cnst_name
   use hk_conv,         only : cmfmca
   use uwshcu,          only : compute_uwshcu_inv
   use unicon_cam,      only : unicon_out_t, unicon_cam_tend

   use time_manager,    only : get_nstep, is_first_step
   use wv_saturation,   only : qsat
   use physconst,       only : latice, latvap, rhoh2o
   use spmd_utils, only : iam

  !water tracers:
   use water_tracer_vars,only: trace_water, wtrc_ntype, wtrc_srfpcp_indices,&
                                wtrc_iatype, wtrc_nwset, wtrc_bulk_indices
   use water_tracers,   only : wtrc_check_h2o
   use water_types,     only : iwtvap, iwtliq, iwtice, iwtcvrain, iwtcvsnow


   implicit none

   ! ---------------------- !
   ! Input-Output Arguments !
   ! ---------------------- !
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state), intent(in)    :: state                           ! Physics state variables
   real(r8),            intent(in)    :: ztodt                           ! 2 delta-t  [ s ]

   type(physics_ptend), intent(out)   :: ptend_all                       ! Indivdual parameterization tendencies
   real(r8),            intent(out)   :: cmfmc2(pcols,pverp)             ! Updraft mass flux by shallow convection [ kg/s/m2 ]
   real(r8),            intent(out)   :: rliq2(pcols)                    ! Vertically-integrated reserved cloud condensate [ m/s ]
   real(r8),            intent(out)   :: qc2(pcols,pver)                 ! Same as qc but only from shallow convection scheme

   

   real(r8),            intent(inout) :: cmfmc(pcols,pverp)    ! Moist deep + shallow convection cloud mass flux [ kg/s/m2 ]
   real(r8),            intent(inout) :: qc(pcols,pver)        ! dq/dt due to export of cloud water into environment by shallow
                                                               ! and deep convection [ kg/kg/s ]
   real(r8),            intent(inout) :: rliq(pcols)           ! Vertical integral of qc [ m/s ]

   real(r8),            intent(in) :: sgh30(pcols)             ! Std. deviation of 30 s orography for tms
   type(cam_in_t),      intent(in) :: cam_in


   !Water tracers:
   real(r8),            intent(inout) :: wtdlf(pcols,pver,wtrc_nwset)    ! dqdt for water tracers due to export of cloud water from conv.

   ! --------------- !
   ! Local Variables ! 
   ! --------------- !
   integer  :: i, k, m
   integer  :: n, x
   integer  :: ilon                                                      ! Global longitude index of a column
   integer  :: ilat                                                      ! Global latitude  index of a column
   integer  :: lchnk                                                     ! Chunk identifier
   integer  :: ncol                                                      ! Number of atmospheric columns
   integer  :: nstep                                                     ! Current time step index
   integer  :: ixcldice, ixcldliq                                        ! Constituent indices for cloud liquid and ice water.
   integer  :: ixnumice, ixnumliq                                        ! Constituent indices for cloud liquid and ice number concentration

   real(r8),  pointer   :: precc(:)                                      ! Shallow convective precipitation (rain+snow) rate at surface [ m/s ]
   real(r8),  pointer   :: snow(:)                                       ! Shallow convective snow rate at surface [ m/s ]

   real(r8), target :: ftem(pcols,pver)                                  ! Temporary workspace for outfld variables
   real(r8) :: cnt2(pcols)                                               ! Top level of shallow convective activity
   real(r8) :: cnb2(pcols)                                               ! Bottom level of convective activity
   real(r8) :: tpert(pcols)                                              ! PBL perturbation theta

   real(r8), pointer   :: pblh(:)                                        ! PBL height [ m ]
   real(r8), pointer   :: qpert(:,:)                                     ! PBL perturbation specific humidity

   ! Temperature tendency from shallow convection (pbuf pointer).
   real(r8), pointer, dimension(:,:) :: ttend_sh

   real(r8) :: ntprprd(pcols,pver)                                       ! Net precip production in layer
   real(r8) :: ntsnprd(pcols,pver)                                       ! Net snow   production in layer
   real(r8) :: tend_s_snwprd(pcols,pver)                                 ! Heating rate of snow production
   real(r8) :: tend_s_snwevmlt(pcols,pver)                               ! Heating rate of evap/melting of snow
   real(r8) :: slflx(pcols,pverp)                                        ! Shallow convective liquid water static energy flux
   real(r8) :: qtflx(pcols,pverp)                                        ! Shallow convective total water flux
   real(r8) :: cmfdqs(pcols, pver)                                       ! Shallow convective snow production
   real(r8) :: zero(pcols)                                               ! Array of zeros
   real(r8) :: cbmf(pcols)                                               ! Shallow cloud base mass flux [ kg/s/m2 ]
   real(r8) :: freqsh(pcols)                                             ! Frequency of shallow convection occurence
   real(r8) :: pcnt(pcols)                                               ! Top    pressure level of shallow + deep convective activity
   real(r8) :: pcnb(pcols)                                               ! Bottom pressure level of shallow + deep convective activity
   real(r8) :: cmfsl(pcols,pverp )                                       ! Convective flux of liquid water static energy
   real(r8) :: cmflq(pcols,pverp )                                       ! Convective flux of total water in energy unit
   real(r8) :: rprdtot(pcols,pver)                                       ! Total shallow+deep rain production tendency
   
   real(r8), target :: ftem_preCu(pcols,pver)                            ! Saturation vapor pressure after shallow Cu convection
   real(r8), target :: tem2(pcols,pver)                                  ! Saturation specific humidity and RH
   real(r8), target :: t_preCu(pcols,pver)                               ! Temperature after shallow Cu convection
   real(r8), target :: tten(pcols,pver)                                  ! Temperature tendency after shallow Cu convection
   real(r8), target :: rhten(pcols,pver)                                 ! RH tendency after shallow Cu convection
   real(r8) :: iccmr_UW(pcols,pver)                                      ! In-cloud Cumulus LWC+IWC [ kg/m2 ]
   real(r8) :: icwmr_UW(pcols,pver)                                      ! In-cloud Cumulus LWC     [ kg/m2 ]
   real(r8) :: icimr_UW(pcols,pver)                                      ! In-cloud Cumulus IWC     [ kg/m2 ]
   real(r8) :: ptend_tracer(pcols,pver,pcnst)                            ! Tendencies of tracers
   real(r8) :: sum1, sum2, sum3, pdelx 
   real(r8) :: landfracdum(pcols)

   !water tracer variables:
   !**********************
   integer  :: wtpcidx                       !Physics Buffer index
   integer  :: wtsnidx                       !Physics Buffer index
   real(r8), pointer,dimension(:) :: wtprec  !tracer total convective precipitation
   real(r8), pointer,dimension(:) :: wtsnow  !tracer total convective snow
   real(r8) :: wtprect(pcols,pcnst)          !Water tracer surface precipitation
   real(r8) :: wtsnowt(pcols,pcnst)          !Water tracer surface snow
   real(r8) :: evpstore(pcols,pver)          !Precipitation Evaporation
   real(r8) :: substore(pcols,pver)          !Snow Sublimation
   real(r8) :: wtqc(pcols,pver,pcnst)        !tendency of detrained cloud condensate
   logical  :: isOk                          !Used to check mass balance
   !**********************

   real(r8), target, dimension(pcols,pver) :: sl, qt, slv
   real(r8), target, dimension(pcols,pver) :: sl_preCu, qt_preCu, slv_preCu

   type(physics_state), target :: state1                                 ! Locally modify for evaporation to use, not returned
   type(physics_ptend) :: ptend_loc                                      ! Local tendency from processes, added up to return as ptend_all

   integer itim_old, ifld
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: concld
   real(r8), pointer, dimension(:,:) :: icwmr                            ! In cloud water + ice mixing ratio
   real(r8), pointer, dimension(:,:) :: rprddp                           ! dq/dt due to deep convective rainout
   real(r8), pointer, dimension(:,:) :: rprdsh                           ! dq/dt due to deep and shallow convective rainout
   real(r8), pointer, dimension(:,:) :: evapcsh                          ! Evaporation of shallow convective precipitation >= 0.
   real(r8), pointer, dimension(:)   :: cnt
   real(r8), pointer, dimension(:)   :: cnb
   real(r8), pointer, dimension(:)   :: cush
   real(r8), pointer, dimension(:,:) :: tke
   real(r8), pointer, dimension(:,:) :: shfrc
   real(r8), pointer, dimension(:,:) :: flxprec                          ! Shallow convective-scale flux of precip (rain+snow) at interfaces [ kg/m2/s ]
   real(r8), pointer, dimension(:,:) :: flxsnow                          ! Shallow convective-scale flux of snow at interfaces [ kg/m2/s ]
   real(r8), pointer, dimension(:,:) :: sh_cldliq
   real(r8), pointer, dimension(:,:) :: sh_cldice

   logical                           :: lq(pcnst)

   type(unicon_out_t) :: unicon_out
   integer :: scheme_code

   ! ----------------------- !
   ! Main Computation Begins ! 
   ! ----------------------- !

   zero  = 0._r8
   nstep = get_nstep()
   lchnk = state%lchnk
   ncol  = state%ncol

   call convect_shallow_select_impl()
   if (.not. use_native_impl) then
      call convect_shallow_select_codon_scheme()
      call convect_shallow_log_tend_direct()
   end if
  
   call physics_state_copy( state, state1 )          ! Copy state to local state1.

   ! Associate pointers with physics buffer fields


   itim_old   =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, concld_idx,      concld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, icwmrsh_idx,     icwmr)

   call pbuf_get_field(pbuf, rprddp_idx,      rprddp )

   call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh )

   call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh  )

   call pbuf_get_field(pbuf, cldtop_idx,      cnt )

   call pbuf_get_field(pbuf, cldbot_idx,      cnb )

   call pbuf_get_field(pbuf, prec_sh_idx,   precc )

   call pbuf_get_field(pbuf, snow_sh_idx,    snow )

   if( convect_shallow_use_shfrc() ) then
       call pbuf_get_field(pbuf, shfrc_idx,  shfrc  )
   endif

   ! Initialization


   call cnst_get_ind( 'CLDLIQ', ixcldliq )
   call cnst_get_ind( 'CLDICE', ixcldice )

   call pbuf_get_field(pbuf, pblh_idx, pblh)

   !  This field probably should reference the pbuf tpert field but it doesnt
   call convect_shallow_init_shell(ncol, tpert, landfracdum)

   if (use_native_impl) then
      select case (shallow_scheme)
      case('off', 'CLUBB_SGS')
         scheme_code = 1
      case('Hack')
         scheme_code = 2
      case('UW')
         scheme_code = 3
      case('UNICON')
         scheme_code = 4
      case default
         scheme_code = -1
      end select
   else
      scheme_code = codon_scheme_code
   end if

   select case (scheme_code)

   case(1) ! off, CLUBB_SGS

      call convect_shallow_ptend_lq_mask_shell(pcnst, lq)
      call physics_ptend_init( ptend_loc, state%psetcols, 'convect_shallow (off)', ls=.true., lq=lq ) ! Initialize local ptend type

      cmfmc2      = 0._r8
      ptend_loc%q = 0._r8
      ptend_loc%s = 0._r8
      rprdsh      = 0._r8
      cmfdqs      = 0._r8
      precc       = 0._r8
      slflx       = 0._r8
      qtflx       = 0._r8
      icwmr       = 0._r8
      rliq2       = 0._r8
      qc2         = 0._r8
      cmfsl       = 0._r8
      cmflq       = 0._r8
      cnt2        = pver
      cnb2        = 1._r8
      evapcsh     = 0._r8
      snow        = 0._r8

     !water tracers:
      wtqc(:,:,:) = 0._r8
      wtprect(:,:)= 0._r8 
      wtsnowt(:,:)= 0._r8

   case(2) ! Hack scheme
                                   
      call convect_shallow_ptend_lq_mask_shell(pcnst, lq)
      call physics_ptend_init( ptend_loc, state%psetcols, 'cmfmca', ls=.true., lq=lq  ) ! Initialize local ptend type

      call pbuf_get_field(pbuf, qpert_idx, qpert)
      qpert(:ncol,2:pcnst) = 0._r8

      call cmfmca( lchnk        ,  ncol         ,                                               &
                   nstep        ,  ztodt        ,  state%pmid ,  state%pdel  ,                  &
                   state%rpdel  ,  state%zm     ,  tpert      ,  qpert       ,  state%phis  ,   &
                   pblh         ,  state%t      ,  state%q    ,  ptend_loc%s ,  ptend_loc%q ,   &
                   cmfmc2       ,  rprdsh       ,  cmfsl      ,  cmflq       ,  precc       ,   &
                   qc2          ,  cnt2         ,  cnb2       ,  icwmr       ,  rliq2       ,   & 
                   state%pmiddry,  state%pdeldry,  state%rpdeldry )

   case(3)   ! UW shallow convection scheme

      ! -------------------------------------- !
      ! uwshcu does momentum transport as well !
      ! -------------------------------------- !

      ! Initialize local ptend type
      call convect_shallow_ptend_lq_mask_shell(pcnst, lq)
      call physics_ptend_init( ptend_loc, state%psetcols, 'UWSHCU', ls=.true., lu=.true., lv=.true., lq=lq  ) 

      call pbuf_get_field(pbuf, cush_idx, cush  ,(/1,itim_old/),  (/pcols,1/))
      call pbuf_get_field(pbuf, tke_idx,  tke)


      call pbuf_get_field(pbuf, sh_flxprc_idx, flxprec)
      call pbuf_get_field(pbuf, sh_flxsnw_idx, flxsnow)

      call compute_uwshcu_inv( pcols     , pver    , ncol           , pcnst         , ztodt         ,                   &
                               state%pint, state%zi, state%pmid     , state%zm      , state%pdel    ,                   & 
                               state%u   , state%v , state%q(:,:,1) , state%q(:,:,ixcldliq), state%q(:,:,ixcldice),     &
                               state%t   , state%s , state%q(:,:,:) ,                                                   &
                               tke       , cld     , concld         , pblh          , cush          ,                   &
                               cmfmc2    , slflx   , qtflx          , 							&
			       flxprec, flxsnow, 			         					&
                               ptend_loc%q(:,:,1)  , ptend_loc%q(:,:,ixcldliq), ptend_loc%q(:,:,ixcldice),              &
                               ptend_loc%s         , ptend_loc%u    , ptend_loc%v   , ptend_tracer  ,                   &
                               rprdsh              , cmfdqs         , precc         , snow          ,                   &
                               evapcsh             , shfrc          , iccmr_UW      , icwmr_UW      ,                   &
                               icimr_UW            , cbmf           , qc2           , rliq2         ,                   &
                               cnt2                , cnb2           , lchnk         , state%pdeldry ,                   &
                               wtprect             , wtsnowt        , wtqc )

      ! --------------------------------------------------------------------- !
      ! Here, 'rprdsh = qrten', 'cmfdqs = qsten' both in unit of [ kg/kg/s ]  !
      ! In addition, define 'icwmr' which includes both liquid and ice.       !
      ! --------------------------------------------------------------------- !

      if (use_native_impl) then
         icwmr(:ncol,:)  = iccmr_UW(:ncol,:)
         rprdsh(:ncol,:) = rprdsh(:ncol,:) + cmfdqs(:ncol,:)
         do m = 4, pcnst
            ptend_loc%q(:ncol,:pver,m) = ptend_tracer(:ncol,:pver,m)
         end do
      end if

      ! Conservation check
      
      !  do i = 1, ncol
      !  do m = 1, pcnst
      !     sum1 = 0._r8
      !     sum2 = 0._r8
      !     sum3 = 0._r8
      !  do k = 1, pver
      !       if(cnst_get_type_byind(m).eq.'wet') then
      !          pdelx = state%pdel(i,k)
      !       else
      !          pdelx = state%pdeldry(i,k)
      !       endif
      !       sum1 = sum1 + state%q(i,k,m)*pdelx
      !       sum2 = sum2 +(state%q(i,k,m)+ptend_loc%q(i,k,m)*ztodt)*pdelx  
      !       sum3 = sum3 + ptend_loc%q(i,k,m)*pdelx 
      !  enddo
      !  if( m .gt. 3 .and. abs(sum1) .gt. 1.e-13_r8 .and. abs(sum2-sum1)/sum1 .gt. 1.e-12_r8 ) then
      !! if( m .gt. 3 .and. abs(sum3) .gt. 1.e-13_r8 ) then
      !      write(iulog,*) 'Sungsu : convect_shallow.F90 does not conserve tracers : ', m, sum1, sum2, abs(sum2-sum1)/sum1
      !!     write(iulog,*) 'Sungsu : convect_shallow.F90 does not conserve tracers : ', m, sum3
      !  endif
      !  enddo
      !  enddo

      ! ------------------------------------------------- !
      ! Convective fluxes of 'sl' and 'qt' in energy unit !
      ! ------------------------------------------------- !

      if (use_native_impl) then
         cmfsl(:ncol,:) = slflx(:ncol,:)
         cmflq(:ncol,:) = qtflx(:ncol,:) * latvap
      end if

      call outfld( 'PRECSH' , precc  , pcols, lchnk )


   case(4)

      icwmr = 0.0_r8

      call unicon_cam_tend(ztodt, state, cam_in, sgh30, &
                           pbuf, ptend_loc, unicon_out)

      cmfmc2(:ncol,:) = unicon_out%cmfmc(:ncol,:)
      qc2(:ncol,:)    = unicon_out%rqc(:ncol,:)
      rliq2(:ncol)    = unicon_out%rliq(:ncol)
      cnt2(:ncol)     = unicon_out%cnt(:ncol)
      cnb2(:ncol)     = unicon_out%cnb(:ncol)

      ! ------------------------------------------------- !
      ! Convective fluxes of 'sl' and 'qt' in energy unit !
      ! ------------------------------------------------- !

      cmfsl(:ncol,:) = unicon_out%slflx(:ncol,:)
      cmflq(:ncol,:) = unicon_out%qtflx(:ncol,:) * latvap

      call outfld( 'PRECSH' , precc  , pcols, lchnk )

   end select

   ! --------------------------------------------------------!     
   ! Calculate fractional occurance of shallow convection    !
   ! --------------------------------------------------------!

 ! Modification : I should check whether below computation of freqsh is correct.

   if (.not. use_native_impl .and. scheme_code == 3) then
      call convect_shallow_uw_post_shell(ncol, state%pmid, cmfmc, cmfmc2, cnt, cnt2, cnb, cnb2, pcnt, pcnb, qc, qc2, rliq, &
           rliq2, wtqc, wtdlf, freqsh, icwmr, iccmr_UW, rprdsh, cmfdqs, ptend_loc%q, ptend_tracer, cmfsl, cmflq, slflx, &
           qtflx, rprddp, rprdtot, ptend_loc%s, ftem)
   else
      if (use_native_impl) then
         freqsh(:) = 0._r8
         do i = 1, ncol
            if( maxval(cmfmc2(i,:pver)) <= 0._r8 ) then
                freqsh(i) = 1._r8
            end if
         end do

         cmfmc(:ncol,:) = cmfmc(:ncol,:) + cmfmc2(:ncol,:)

         do i = 1, ncol
            if( cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
            if( cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
            pcnt(i) = state%pmid(i,int(cnt(i)))
            pcnb(i) = state%pmid(i,int(cnb(i)))
         end do

         qc(:ncol,:pver) = qc(:ncol,:pver) + qc2(:ncol,:pver)
         rliq(:ncol)     = rliq(:ncol) + rliq2(:ncol)

         do m=1,wtrc_nwset
           wtdlf(:ncol,:pver,m) = wtdlf(:ncol,:pver,m) + &
                                 (wtqc(:ncol,:pver,wtrc_iatype(m,iwtliq)) + wtqc(:ncol,:pver,wtrc_iatype(m,iwtice)))
         end do
      else
         call convect_shallow_postmerge(ncol, state%pmid, cmfmc, cmfmc2, cnt, cnt2, cnb, cnb2, pcnt, pcnb, qc, qc2, rliq, &
              rliq2, wtqc, wtdlf, freqsh)
      end if

      rprdtot(:ncol,:pver) = rprdsh(:ncol,:pver) + rprddp(:ncol,:pver)
      ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair
   end if
   
   ! ----------------------------------------------- !
   ! This quantity was previously known as CMFDQR.   !
   ! Now CMFDQR is the shallow rain production only. !
   ! ----------------------------------------------- !

   
   call pbuf_set_field(pbuf, rprdtot_idx, rprdtot(:ncol,:pver), start=(/1,1/), kount=(/ncol,pver/))
 
   ! ----------------------------------------------------------------------- ! 
   ! Add shallow reserved cloud condensate to deep reserved cloud condensate !
   !     qc [ kg/kg/s] , rliq [ m/s ]                                        !
   ! ----------------------------------------------------------------------- !

   ! ---------------------------------------------------------------------------- !
   ! Output new partition of cloud condensate variables, as well as precipitation !
   ! ---------------------------------------------------------------------------- ! 

   if( microp_scheme == 'MG' ) then
       call cnst_get_ind( 'NUMLIQ', ixnumliq )
       call cnst_get_ind( 'NUMICE', ixnumice )
   endif

   call outfld( 'ICWMRSH ', icwmr                    , pcols   , lchnk )

   call outfld( 'CMFDT  ', ftem                      , pcols   , lchnk )
   call outfld( 'CMFDQ  ', ptend_loc%q(1,1,1)        , pcols   , lchnk )
   call outfld( 'CMFDICE', ptend_loc%q(1,1,ixcldice) , pcols   , lchnk )
   call outfld( 'CMFDLIQ', ptend_loc%q(1,1,ixcldliq) , pcols   , lchnk )
   call outfld( 'CMFMC'  , cmfmc                     , pcols   , lchnk )
   call outfld( 'QC'     , qc2                       , pcols   , lchnk )
   call outfld( 'CMFDQR' , rprdsh                    , pcols   , lchnk )
   call outfld( 'CMFSL'  , cmfsl                     , pcols   , lchnk )
   call outfld( 'CMFLQ'  , cmflq                     , pcols   , lchnk )
   call outfld( 'DQP'    , qc2                       , pcols   , lchnk )
   call outfld( 'CLDTOP' , cnt                       , pcols   , lchnk )
   call outfld( 'CLDBOT' , cnb                       , pcols   , lchnk )
   call outfld( 'PCLDTOP', pcnt                      , pcols   , lchnk )
   call outfld( 'PCLDBOT', pcnb                      , pcols   , lchnk )  
   call outfld( 'FREQSH' , freqsh                    , pcols   , lchnk )

   if (scheme_code == 3) then
      call outfld( 'CBMF'   , cbmf                      , pcols   , lchnk )
      call outfld( 'UWFLXPRC', flxprec                  , pcols   , lchnk )  
      call outfld( 'UWFLXSNW' , flxsnow                 , pcols   , lchnk )
   endif

   ! ---------------------------------------------------------------- !
   ! Add tendency from this process to tend from other processes here !
   ! ---------------------------------------------------------------- !

   call physics_ptend_init(ptend_all, state1%psetcols, 'convect_shallow')
   call physics_ptend_sum( ptend_loc, ptend_all, ncol )

   ! ----------------------------------------------------------------------------- !
   ! For diagnostic purpose, print out 'QT,SL,SLV,T,RH' just before cumulus scheme !
   ! ----------------------------------------------------------------------------- !

   if (use_native_impl) then
      sl_preCu(:ncol,:pver)  = state1%s(:ncol,:pver) -   latvap           * state1%q(:ncol,:pver,ixcldliq) &
                                                     - ( latvap + latice) * state1%q(:ncol,:pver,ixcldice)
      qt_preCu(:ncol,:pver)  = state1%q(:ncol,:pver,1) + state1%q(:ncol,:pver,ixcldliq) &
                                                       + state1%q(:ncol,:pver,ixcldice)
      slv_preCu(:ncol,:pver) = sl_preCu(:ncol,:pver) * ( 1._r8 + zvir * qt_preCu(:ncol,:pver) )

      t_preCu(:ncol,:)       = state1%t(:ncol,:pver)
   else
      call convect_shallow_diag_shell(1, ncol, ixcldliq, ixcldice, ztodt, state1%s, state1%t, state1%q, ftem, &
           sl_preCu, qt_preCu, slv_preCu, t_preCu, ftem_preCu, tten, rhten)
   end if
   call qsat(state1%t(:ncol,:), state1%pmid(:ncol,:), &
        tem2(:ncol,:), ftem(:ncol,:))
   if (use_native_impl) then
      ftem_preCu(:ncol,:)    = state1%q(:ncol,:,1) / ftem(:ncol,:) * 100._r8
   else
      call convect_shallow_diag_shell(2, ncol, ixcldliq, ixcldice, ztodt, state1%s, state1%t, state1%q, ftem, &
           sl_preCu, qt_preCu, slv_preCu, t_preCu, ftem_preCu, tten, rhten)
   end if

   call outfld( 'qt_pre_Cu      ', qt_preCu               , pcols, lchnk )
   call outfld( 'sl_pre_Cu      ', sl_preCu               , pcols, lchnk )
   call outfld( 'slv_pre_Cu     ', slv_preCu              , pcols, lchnk )
   call outfld( 'u_pre_Cu       ', state1%u               , pcols, lchnk )
   call outfld( 'v_pre_Cu       ', state1%v               , pcols, lchnk )
   call outfld( 'qv_pre_Cu      ', state1%q(:,:,1)        , pcols, lchnk )
   call outfld( 'ql_pre_Cu      ', state1%q(:,:,ixcldliq) , pcols, lchnk )
   call outfld( 'qi_pre_Cu      ', state1%q(:,:,ixcldice) , pcols, lchnk )
   call outfld( 't_pre_Cu       ', state1%t               , pcols, lchnk )
   call outfld( 'rh_pre_Cu      ', ftem_preCu             , pcols, lchnk )

   ! ----------------------------------------------- ! 
   ! Update physics state type state1 with ptend_loc ! 
   ! ----------------------------------------------- !

   call physics_update( state1, ptend_loc, ztodt )

  !***********************
  !Check water tracer mass
  !***********************
  if (trace_water) then
    isOk = wtrc_check_h2o("after-shallow UW", state1, state1%q, ztodt)

    !Check precipitation:
    !write(iulog,*) 'wtprect',sum(precc(:)),sum(wtprect(:,wtrc_iatype(1,iwtvap))) !,sum(wtprect(:,wtrc_iatype(2,iwtvap)))
    !write(iulog,*) 'wtsnowt',sum(snow(:)),sum(wtsnowt(:,wtrc_iatype(1,iwtvap)))
  end if
  !***********************
  !assign values to physics buffer variables
  !*********************** 
    do m=1,wtrc_ntype(iwtcvrain)
       call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvrain,m), wtprec)
       call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvsnow,m), wtsnow)
       if (use_native_impl) then
          wtprec(:) = wtprec(:) + (wtprect(:,wtrc_iatype(m,iwtvap)) - wtsnowt(:,wtrc_iatype(m,iwtvap))) !assign values (should be rain only)
          wtsnow(:) = wtsnow(:) + wtsnowt(:,wtrc_iatype(m,iwtvap))                                      !(snow only)
       else
          call convect_shallow_wtrc_precip_shell(wtrc_iatype(m,iwtvap), wtprect, wtsnowt, wtprec, wtsnow)
       end if
    end do
  !**********************

   ! ----------------------------------------------------------------------------- !
   ! For diagnostic purpose, print out 'QT,SL,SLV,t,RH' just after cumulus scheme  !
   ! ----------------------------------------------------------------------------- !

   if (use_native_impl) then
      sl(:ncol,:pver)  = state1%s(:ncol,:pver) -   latvap           * state1%q(:ncol,:pver,ixcldliq) &
                                               - ( latvap + latice) * state1%q(:ncol,:pver,ixcldice)
      qt(:ncol,:pver)  = state1%q(:ncol,:pver,1) + state1%q(:ncol,:pver,ixcldliq) &
                                                 + state1%q(:ncol,:pver,ixcldice)
      slv(:ncol,:pver) = sl(:ncol,:pver) * ( 1._r8 + zvir * qt(:ncol,:pver) )
   else
      call convect_shallow_diag_shell(3, ncol, ixcldliq, ixcldice, ztodt, state1%s, state1%t, state1%q, ftem, &
           sl, qt, slv, t_preCu, ftem_preCu, tten, rhten)
   end if

   call qsat(state1%t(:ncol,:), state1%pmid(:ncol,:), &
        tem2(:ncol,:), ftem(:ncol,:))
   if (use_native_impl) then
      ftem(:ncol,:)    = state1%q(:ncol,:,1) / ftem(:ncol,:) * 100._r8
   else
      call convect_shallow_diag_shell(4, ncol, ixcldliq, ixcldice, ztodt, state1%s, state1%t, state1%q, ftem, &
           sl, qt, slv, t_preCu, ftem_preCu, tten, rhten)
   end if

   call outfld( 'qt_aft_Cu      ', qt                     , pcols, lchnk )
   call outfld( 'sl_aft_Cu      ', sl                     , pcols, lchnk )
   call outfld( 'slv_aft_Cu     ', slv                    , pcols, lchnk )
   call outfld( 'u_aft_Cu       ', state1%u               , pcols, lchnk )
   call outfld( 'v_aft_Cu       ', state1%v               , pcols, lchnk )
   call outfld( 'qv_aft_Cu      ', state1%q(:,:,1)        , pcols, lchnk )
   call outfld( 'ql_aft_Cu      ', state1%q(:,:,ixcldliq) , pcols, lchnk )
   call outfld( 'qi_aft_Cu      ', state1%q(:,:,ixcldice) , pcols, lchnk )
   call outfld( 't_aft_Cu       ', state1%t               , pcols, lchnk )
   call outfld( 'rh_aft_Cu      ', ftem                   , pcols, lchnk )

   if (use_native_impl) then
      tten(:ncol,:)  = ( state1%t(:ncol,:pver) - t_preCu(:ncol,:) ) / ztodt
      rhten(:ncol,:) = ( ftem(:ncol,:) - ftem_preCu(:ncol,:) ) / ztodt
   end if

   call outfld( 'tten_Cu        ', tten                           , pcols, lchnk )
   call outfld( 'rhten_Cu       ', rhten                          , pcols, lchnk )


   ! ------------------------------------------------------------------------ !
   ! UW-Shallow Cumulus scheme includes                                       !
   ! evaporation physics inside in it. So when 'shallow_scheme = UW', we must !
   ! NOT perform below 'zm_conv_evap'.                                        !
   ! ------------------------------------------------------------------------ !

   if (scheme_code == 2) then

   ! ------------------------------------------------------------------------------- !
   ! Determine the phase of the precipitation produced and add latent heat of fusion !
   ! Evaporate some of the precip directly into the environment (Sundqvist)          !
   ! Allow this to use the updated state1 and a fresh ptend_loc type                 !
   ! Heating and specific humidity tendencies produced                               !
   ! ------------------------------------------------------------------------------- !

   ! --------------------------------- !
   ! initialize ptend for next process !
   ! --------------------------------- !

    lq(1) = .TRUE.
    lq(2:) = .FALSE.
    call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

    call pbuf_get_field(pbuf, sh_flxprc_idx, flxprec    )
    call pbuf_get_field(pbuf, sh_flxsnw_idx, flxsnow    )
    call pbuf_get_field(pbuf, sh_cldliq_idx, sh_cldliq  )
    call pbuf_get_field(pbuf, sh_cldice_idx, sh_cldice  )

    !! clouds have no water... :)
    sh_cldliq(:ncol,:) = 0._r8
    sh_cldice(:ncol,:) = 0._r8

    call zm_conv_evap( state1%ncol, state1%lchnk,                                    &
                       state1%t, state1%pmid, state1%pdel, state1%q(:pcols,:pver,1), &
		       landfracdum, &
                       ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt,                  & 
                       ptend_loc%q(:pcols,:pver,1),                                  &
                       rprdsh, cld, ztodt, precc, snow,                         &
                       evpstore, substore, ntprprd, ntsnprd, flxprec, flxsnow )

   ! ------------------------------------------ !
   ! record history variables from zm_conv_evap !
   ! ------------------------------------------ !

   evapcsh(:ncol,:pver) = ptend_loc%q(:ncol,:pver,1)

   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver) / cpair
   call outfld( 'EVAPTCM '       , ftem                           , pcols, lchnk )
   ftem(:ncol,:pver) = tend_s_snwprd(:ncol,:pver) / cpair
   call outfld( 'FZSNTCM '       , ftem                           , pcols, lchnk )
   ftem(:ncol,:pver) = tend_s_snwevmlt(:ncol,:pver) / cpair
   call outfld( 'EVSNTCM '       , ftem                           , pcols, lchnk )
   call outfld( 'EVAPQCM '       , ptend_loc%q(1,1,1)             , pcols, lchnk )
   call outfld( 'PRECSH  '       , precc                          , pcols, lchnk )
   call outfld( 'HKFLXPRC'       , flxprec                        , pcols, lchnk )
   call outfld( 'HKFLXSNW'       , flxsnow                        , pcols, lchnk )
   call outfld( 'HKNTPRPD'       , ntprprd                        , pcols, lchnk )
   call outfld( 'HKNTSNPD'       , ntsnprd                        , pcols, lchnk )
   call outfld( 'HKEIHEAT'       , ptend_loc%s                    , pcols, lchnk )

   ! ---------------------------------------------------------------- !      
   ! Add tendency from this process to tend from other processes here !
   ! ---------------------------------------------------------------- !

   call physics_ptend_sum( ptend_loc, ptend_all, ncol )
   call physics_ptend_dealloc(ptend_loc)

   ! -------------------------------------------- !
   ! Do not perform evaporation process for UW-Cu !
   ! -------------------------------------------- !

   end if

   ! ------------------------------------------------------------- !
   ! Update name of parameterization tendencies to send to tphysbc !
   ! ------------------------------------------------------------- !

   call physics_state_dealloc(state1)

   ! If we added temperature tendency to pbuf, set it now.
   if (ttend_sh_idx > 0) then
      call pbuf_get_field(pbuf, ttend_sh_idx, ttend_sh)
      ttend_sh(:ncol,:pver) = ptend_all%s(:ncol,:pver)/cpair
   end if

  end subroutine convect_shallow_tend

  !=============================================================================== !

subroutine convect_shallow_select_codon_scheme()

   use iso_c_binding, only: c_int64_t, c_loc

   integer :: i
   integer(c_int64_t), target :: scheme_ascii(len(shallow_scheme))
   integer(c_int64_t), target :: scheme_code_c, status_c

   if (codon_scheme_selected) return

   do i = 1, len(shallow_scheme)
      scheme_ascii(i) = int(iachar(shallow_scheme(i:i)), c_int64_t)
   end do

   call convect_shallow_select_scheme_codon(int(len(shallow_scheme), c_int64_t), c_loc(scheme_ascii(1)), &
        c_loc(scheme_code_c), c_loc(status_c))

   if (status_c == 0_c_int64_t) then
      codon_scheme_code = int(scheme_code_c)
   else
      codon_scheme_code = -1
   end if

   codon_scheme_selected = .true.

end subroutine convect_shallow_select_codon_scheme

subroutine convect_shallow_log_tend_direct()

   use spmd_utils, only: masterproc

   if (convect_shallow_tend_logged) return
   convect_shallow_tend_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') &
           'convect_shallow_tend direct = codon; scheme/action dispatch and outer stage shells direct = codon; ' // &
           'cmfmca/compute_uwshcu_inv/unicon/zm_conv_evap/pbuf/outfld/physics_update native CAM API islands'
      call convect_shallow_append_proof( &
           'convect_shallow_tend direct = codon; scheme/action dispatch and outer stage shells direct = codon; ' // &
           'cmfmca/compute_uwshcu_inv/unicon/zm_conv_evap/pbuf/outfld/physics_update native CAM API islands')
      call flush(iulog)
   end if

end subroutine convect_shallow_log_tend_direct

  !=============================================================================== !

subroutine convect_shallow_init_select_impl()

   use spmd_utils, only: masterproc

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (init_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CONVECT_SHALLOW_INIT_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_init_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_init_impl = .false.
   end if

   init_impl_selected = .true.

   if (masterproc) then
      if (use_native_init_impl) then
         write(iulog,*) 'convect_shallow_init implementation = native'
         call convect_shallow_append_proof('convect_shallow_init implementation = native')
      else
         write(iulog,*) 'convect_shallow_init implementation = codon'
         call convect_shallow_append_proof('convect_shallow_init implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine convect_shallow_init_select_impl

subroutine convect_shallow_log_init_mw_ratio_entered()

   use spmd_utils, only: masterproc

   if (init_mw_ratio_logged) return
   init_mw_ratio_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'convect_shallow_init mw ratio entered (UW molecular-weight ratio = codon)'
      call convect_shallow_append_proof('convect_shallow_init mw ratio entered (UW molecular-weight ratio = codon)')
      call flush(iulog)
   end if

end subroutine convect_shallow_log_init_mw_ratio_entered

subroutine convect_shallow_log_init_direct()

   use spmd_utils, only: masterproc

   if (convect_shallow_init_direct_logged) return
   convect_shallow_init_direct_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') &
           'convect_shallow_init direct = codon; scheme/action dispatch direct = codon; addfld/pbuf/init_uwshcu/unicon native CAM API islands'
      call convect_shallow_append_proof( &
           'convect_shallow_init direct = codon; scheme/action dispatch direct = codon; addfld/pbuf/init_uwshcu/unicon native CAM API islands')
      call flush(iulog)
   end if

end subroutine convect_shallow_log_init_direct

  !=============================================================================== !

subroutine convect_shallow_select_impl()

   use spmd_utils, only: masterproc

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CONVECT_SHALLOW_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_impl = .false.
   end if

   impl_selected = .true.

   if (masterproc) then
      if (use_native_impl) then
         write(iulog,*) 'convect_shallow_tend implementation = native'
         call convect_shallow_append_proof('convect_shallow selector entered implementation = native')
      else
         write(iulog,*) 'convect_shallow_tend implementation = codon'
         call convect_shallow_append_proof('convect_shallow selector entered implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine convect_shallow_select_impl

subroutine convect_shallow_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line
   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('CONVECT_SHALLOW_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine convect_shallow_append_proof

subroutine convect_shallow_outer_stage_dispatch_call(stage_c, mode_c, ncol_c, pcols_c, pver_c, pverp_c, pcnst_c, &
     wtrc_nwset_c, ixcldliq_c, ixcldice_c, vap_idx_c, ztodt_c, latvap_c, latice_c, zvir_c, cpair_c, &
     p1_p, p2_p, p3_p, p4_p, p5_p, p6_p, p7_p, p8_p, p9_p, p10_p, p11_p, p12_p, p13_p, p14_p, p15_p, p16_p, &
     p17_p, p18_p, p19_p, p20_p, p21_p, p22_p, p23_p, p24_p, p25_p, p26_p, p27_p, p28_p, p29_p, p30_p, p31_p, p32_p)

   use iso_c_binding, only: c_double, c_int64_t, c_null_ptr, c_ptr

   integer(c_int64_t), intent(in) :: stage_c, mode_c, ncol_c, pcols_c, pver_c, pverp_c, pcnst_c
   integer(c_int64_t), intent(in) :: wtrc_nwset_c, ixcldliq_c, ixcldice_c, vap_idx_c
   real(c_double), intent(in) :: ztodt_c, latvap_c, latice_c, zvir_c, cpair_c
   type(c_ptr), intent(in), optional :: p1_p, p2_p, p3_p, p4_p, p5_p, p6_p, p7_p, p8_p
   type(c_ptr), intent(in), optional :: p9_p, p10_p, p11_p, p12_p, p13_p, p14_p, p15_p, p16_p
   type(c_ptr), intent(in), optional :: p17_p, p18_p, p19_p, p20_p, p21_p, p22_p, p23_p, p24_p
   type(c_ptr), intent(in), optional :: p25_p, p26_p, p27_p, p28_p, p29_p, p30_p, p31_p, p32_p
   type(c_ptr) :: q(32)

   interface
      subroutine convect_shallow_outer_stage_dispatch_codon(stage_c, mode_c, ncol_c, pcols_c, pver_c, pverp_c, &
           pcnst_c, wtrc_nwset_c, ixcldliq_c, ixcldice_c, vap_idx_c, ztodt_c, latvap_c, latice_c, zvir_c, cpair_c, &
           p1_p, p2_p, p3_p, p4_p, p5_p, p6_p, p7_p, p8_p, p9_p, p10_p, p11_p, p12_p, p13_p, p14_p, p15_p, p16_p, &
           p17_p, p18_p, p19_p, p20_p, p21_p, p22_p, p23_p, p24_p, p25_p, p26_p, p27_p, p28_p, p29_p, p30_p, &
           p31_p, p32_p) bind(c, name="convect_shallow_outer_stage_dispatch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, mode_c, ncol_c, pcols_c, pver_c, pverp_c, pcnst_c
         integer(c_int64_t), value :: wtrc_nwset_c, ixcldliq_c, ixcldice_c, vap_idx_c
         real(c_double), value :: ztodt_c, latvap_c, latice_c, zvir_c, cpair_c
         type(c_ptr), value :: p1_p, p2_p, p3_p, p4_p, p5_p, p6_p, p7_p, p8_p
         type(c_ptr), value :: p9_p, p10_p, p11_p, p12_p, p13_p, p14_p, p15_p, p16_p
         type(c_ptr), value :: p17_p, p18_p, p19_p, p20_p, p21_p, p22_p, p23_p, p24_p
         type(c_ptr), value :: p25_p, p26_p, p27_p, p28_p, p29_p, p30_p, p31_p, p32_p
      end subroutine convect_shallow_outer_stage_dispatch_codon
   end interface

   q(:) = c_null_ptr
   if (present(p1_p)) q(1) = p1_p
   if (present(p2_p)) q(2) = p2_p
   if (present(p3_p)) q(3) = p3_p
   if (present(p4_p)) q(4) = p4_p
   if (present(p5_p)) q(5) = p5_p
   if (present(p6_p)) q(6) = p6_p
   if (present(p7_p)) q(7) = p7_p
   if (present(p8_p)) q(8) = p8_p
   if (present(p9_p)) q(9) = p9_p
   if (present(p10_p)) q(10) = p10_p
   if (present(p11_p)) q(11) = p11_p
   if (present(p12_p)) q(12) = p12_p
   if (present(p13_p)) q(13) = p13_p
   if (present(p14_p)) q(14) = p14_p
   if (present(p15_p)) q(15) = p15_p
   if (present(p16_p)) q(16) = p16_p
   if (present(p17_p)) q(17) = p17_p
   if (present(p18_p)) q(18) = p18_p
   if (present(p19_p)) q(19) = p19_p
   if (present(p20_p)) q(20) = p20_p
   if (present(p21_p)) q(21) = p21_p
   if (present(p22_p)) q(22) = p22_p
   if (present(p23_p)) q(23) = p23_p
   if (present(p24_p)) q(24) = p24_p
   if (present(p25_p)) q(25) = p25_p
   if (present(p26_p)) q(26) = p26_p
   if (present(p27_p)) q(27) = p27_p
   if (present(p28_p)) q(28) = p28_p
   if (present(p29_p)) q(29) = p29_p
   if (present(p30_p)) q(30) = p30_p
   if (present(p31_p)) q(31) = p31_p
   if (present(p32_p)) q(32) = p32_p

   call convect_shallow_outer_stage_dispatch_codon(stage_c, mode_c, ncol_c, pcols_c, pver_c, pverp_c, &
        pcnst_c, wtrc_nwset_c, ixcldliq_c, ixcldice_c, vap_idx_c, ztodt_c, latvap_c, latice_c, zvir_c, cpair_c, &
        q(1), q(2), q(3), q(4), q(5), q(6), q(7), q(8), q(9), q(10), q(11), q(12), q(13), q(14), q(15), q(16), &
        q(17), q(18), q(19), q(20), q(21), q(22), q(23), q(24), q(25), q(26), q(27), q(28), q(29), q(30), q(31), q(32))

end subroutine convect_shallow_outer_stage_dispatch_call

subroutine convect_shallow_log_init_shell_entered()

   use spmd_utils, only: masterproc

   if (convect_shallow_init_shell_logged) return
   convect_shallow_init_shell_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'convect_shallow init shell entered (unified shallow-stage dispatch = codon)'
      call convect_shallow_append_proof('convect_shallow init shell entered (unified shallow-stage dispatch = codon)')
      call flush(iulog)
   end if

end subroutine convect_shallow_log_init_shell_entered

subroutine convect_shallow_init_shell(ncol_local, tpert_local, landfracdum_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   integer, intent(in) :: ncol_local
   real(r8), target, intent(inout) :: tpert_local(pcols), landfracdum_local(pcols)

   interface
      subroutine convect_shallow_init_shell_codon(ncol_c, pcols_c, tpert_p, landfracdum_p) &
           bind(c, name="convect_shallow_init_shell_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c
         type(c_ptr), value :: tpert_p, landfracdum_p
      end subroutine convect_shallow_init_shell_codon
   end interface

   if (use_native_impl) then
      call convect_shallow_init_shell_native(ncol_local, tpert_local, landfracdum_local)
      return
   end if

   call convect_shallow_log_init_shell_entered()
   call convect_shallow_outer_stage_dispatch_call(1_c_int64_t, 0_c_int64_t, int(ncol_local, c_int64_t), &
        int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), 0_c_int64_t, &
        0_c_int64_t, 0_c_int64_t, 0_c_int64_t, 0_c_int64_t, 0._r8, 0._r8, 0._r8, 0._r8, 0._r8, &
        c_loc(tpert_local), c_loc(landfracdum_local))

end subroutine convect_shallow_init_shell

subroutine convect_shallow_init_shell_native(ncol_local, tpert_local, landfracdum_local)

   integer, intent(in) :: ncol_local
   real(r8), intent(inout) :: tpert_local(pcols), landfracdum_local(pcols)

   tpert_local(:ncol_local) = 0._r8
   landfracdum_local(:ncol_local) = 0._r8

end subroutine convect_shallow_init_shell_native

subroutine convect_shallow_log_ptend_lq_mask_shell_entered()

   use spmd_utils, only: masterproc

   if (convect_shallow_ptend_lq_mask_shell_logged) return
   convect_shallow_ptend_lq_mask_shell_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'convect_shallow ptend lq mask shell entered (unified shallow-stage dispatch = codon)'
      call convect_shallow_append_proof( &
           'convect_shallow ptend lq mask shell entered (unified shallow-stage dispatch = codon)')
      call flush(iulog)
   end if

end subroutine convect_shallow_log_ptend_lq_mask_shell_entered

subroutine convect_shallow_ptend_lq_mask_shell(pcnst_local, lq_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   integer, intent(in) :: pcnst_local
   logical, intent(out) :: lq_local(pcnst_local)

   integer :: m
   integer(c_int64_t), target :: lq_mask_c(pcnst_local)

   interface
      subroutine convect_shallow_ptend_lq_mask_shell_codon(pcnst_c, lq_mask_p) &
           bind(c, name="convect_shallow_ptend_lq_mask_shell_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: pcnst_c
         type(c_ptr), value :: lq_mask_p
      end subroutine convect_shallow_ptend_lq_mask_shell_codon
   end interface

   if (use_native_impl) then
      lq_local(:) = .TRUE.
      return
   end if

   call convect_shallow_log_ptend_lq_mask_shell_entered()
   call convect_shallow_outer_stage_dispatch_call(2_c_int64_t, 0_c_int64_t, 0_c_int64_t, int(pcols, c_int64_t), &
        int(pver, c_int64_t), int(pverp, c_int64_t), int(pcnst_local, c_int64_t), 0_c_int64_t, &
        0_c_int64_t, 0_c_int64_t, 0_c_int64_t, 0._r8, 0._r8, 0._r8, 0._r8, 0._r8, c_loc(lq_mask_c))

   do m = 1, pcnst_local
      lq_local(m) = lq_mask_c(m) /= 0_c_int64_t
   end do

end subroutine convect_shallow_ptend_lq_mask_shell

subroutine convect_shallow_log_wtrc_precip_shell_entered()

   use spmd_utils, only: masterproc

   if (convect_shallow_wtrc_precip_shell_logged) return
   convect_shallow_wtrc_precip_shell_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'convect_shallow wtrc precip shell entered (unified shallow-stage dispatch = codon)'
      call convect_shallow_append_proof( &
           'convect_shallow wtrc precip shell entered (unified shallow-stage dispatch = codon)')
      call flush(iulog)
   end if

end subroutine convect_shallow_log_wtrc_precip_shell_entered

subroutine convect_shallow_wtrc_precip_shell(vap_idx, wtprect_local, wtsnowt_local, wtprec_local, wtsnow_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr
   use constituents, only: pcnst

   integer, intent(in) :: vap_idx
   real(r8), target, intent(in) :: wtprect_local(pcols,pcnst), wtsnowt_local(pcols,pcnst)
   real(r8), pointer, intent(inout) :: wtprec_local(:), wtsnow_local(:)

   interface
      subroutine convect_shallow_wtrc_precip_shell_codon(pcols_c, vap_idx_c, wtprect_p, wtsnowt_p, &
           wtprec_p, wtsnow_p) bind(c, name="convect_shallow_wtrc_precip_shell_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: pcols_c, vap_idx_c
         type(c_ptr), value :: wtprect_p, wtsnowt_p, wtprec_p, wtsnow_p
      end subroutine convect_shallow_wtrc_precip_shell_codon
   end interface

   call convect_shallow_log_wtrc_precip_shell_entered()
   call convect_shallow_outer_stage_dispatch_call(3_c_int64_t, 0_c_int64_t, 0_c_int64_t, int(pcols, c_int64_t), &
        int(pver, c_int64_t), int(pverp, c_int64_t), 0_c_int64_t, 0_c_int64_t, 0_c_int64_t, &
        0_c_int64_t, int(vap_idx, c_int64_t), 0._r8, 0._r8, 0._r8, 0._r8, 0._r8, &
        c_loc(wtprect_local), c_loc(wtsnowt_local), c_loc(wtprec_local), c_loc(wtsnow_local))

end subroutine convect_shallow_wtrc_precip_shell

subroutine convect_shallow_log_diag_shell_entered()

   use spmd_utils, only: masterproc

   if (convect_shallow_diag_shell_logged) return
   convect_shallow_diag_shell_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'convect_shallow diag shell entered (unified shallow-stage dispatch = codon)'
      call convect_shallow_append_proof('convect_shallow diag shell entered (unified shallow-stage dispatch = codon)')
      call flush(iulog)
   end if

end subroutine convect_shallow_log_diag_shell_entered

subroutine convect_shallow_diag_shell(mode, ncol_local, ixcldliq_local, ixcldice_local, ztodt_local, &
     state_s_local, state_t_local, state_q_local, sat_rh_local, sl_local, qt_local, slv_local, &
     t_precu_local, rh_precu_local, tten_local, rhten_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
   use physconst, only: latvap, latice, zvir
   use constituents, only: pcnst

   integer, intent(in) :: mode, ncol_local, ixcldliq_local, ixcldice_local
   real(r8), intent(in) :: ztodt_local
   real(r8), target, intent(in) :: state_s_local(pcols,pver), state_t_local(pcols,pver)
   real(r8), target, intent(in) :: state_q_local(pcols,pver,pcnst)
   real(r8), target, intent(inout) :: sat_rh_local(pcols,pver)
   real(r8), target, intent(inout) :: sl_local(pcols,pver), qt_local(pcols,pver), slv_local(pcols,pver)
   real(r8), target, intent(inout) :: t_precu_local(pcols,pver), rh_precu_local(pcols,pver)
   real(r8), target, intent(inout) :: tten_local(pcols,pver), rhten_local(pcols,pver)

   interface
      subroutine convect_shallow_diag_shell_codon(mode_c, ncol_c, pcols_c, pver_c, ixcldliq_c, ixcldice_c, &
           ztodt_c, latvap_c, latice_c, zvir_c, state_s_p, state_t_p, state_q_p, sat_rh_p, sl_p, qt_p, slv_p, &
           t_precu_p, rh_precu_p, tten_p, rhten_p) bind(c, name="convect_shallow_diag_shell_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c, ixcldliq_c, ixcldice_c
         real(c_double), value :: ztodt_c, latvap_c, latice_c, zvir_c
         type(c_ptr), value :: state_s_p, state_t_p, state_q_p, sat_rh_p, sl_p, qt_p, slv_p
         type(c_ptr), value :: t_precu_p, rh_precu_p, tten_p, rhten_p
      end subroutine convect_shallow_diag_shell_codon
   end interface

   call convect_shallow_log_diag_shell_entered()

   call convect_shallow_outer_stage_dispatch_call(4_c_int64_t, int(mode, c_int64_t), int(ncol_local, c_int64_t), &
        int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), int(pcnst, c_int64_t), &
        0_c_int64_t, int(ixcldliq_local, c_int64_t), int(ixcldice_local, c_int64_t), 0_c_int64_t, &
        real(ztodt_local, c_double), real(latvap, c_double), real(latice, c_double), real(zvir, c_double), &
        0._r8, c_loc(state_s_local), c_loc(state_t_local), c_loc(state_q_local), c_loc(sat_rh_local), &
        c_loc(sl_local), c_loc(qt_local), c_loc(slv_local), c_loc(t_precu_local), c_loc(rh_precu_local), &
        c_loc(tten_local), c_loc(rhten_local))

end subroutine convect_shallow_diag_shell

subroutine convect_shallow_log_uw_post_shell_entered()

   use spmd_utils, only: masterproc

   if (convect_shallow_uw_post_shell_logged) return
   convect_shallow_uw_post_shell_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'convect_shallow uw post shell entered (unified shallow-stage dispatch = codon)'
      call convect_shallow_append_proof('convect_shallow uw post shell entered (unified shallow-stage dispatch = codon)')
      call flush(iulog)
   end if

end subroutine convect_shallow_log_uw_post_shell_entered

subroutine convect_shallow_uw_post_shell(ncol_local, state_pmid_local, cmfmc_local, cmfmc2_local, cnt_local, cnt2_local, cnb_local, &
     cnb2_local, pcnt_local, pcnb_local, qc_local, qc2_local, rliq_local, rliq2_local, wtqc_local, wtdlf_local, freqsh_local, &
     icwmr_local, iccmr_uw_local, rprdsh_local, cmfdqs_local, ptend_q_local, ptend_tracer_local, cmfsl_local, cmflq_local, &
     slflx_local, qtflx_local, rprddp_local, rprdtot_local, ptend_s_local, ftem_local)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
   use physconst, only: latvap, cpair
   use constituents, only: pcnst
   use water_tracer_vars, only: wtrc_nwset, wtrc_iatype
   use water_types, only: iwtliq, iwtice

   integer, intent(in) :: ncol_local
   real(r8), target, intent(in) :: state_pmid_local(pcols,pver), cmfmc2_local(pcols,pverp), cnt2_local(pcols), cnb2_local(pcols)
   real(r8), target, intent(in) :: qc2_local(pcols,pver), rliq2_local(pcols), wtqc_local(pcols,pver,pcnst)
   real(r8), target, intent(in) :: iccmr_uw_local(pcols,pver), cmfdqs_local(pcols,pver), ptend_tracer_local(pcols,pver,pcnst)
   real(r8), target, intent(in) :: slflx_local(pcols,pverp), qtflx_local(pcols,pverp), rprddp_local(pcols,pver), ptend_s_local(pcols,pver)
   real(r8), target, intent(inout) :: cmfmc_local(pcols,pverp), cnt_local(pcols), cnb_local(pcols), pcnt_local(pcols), pcnb_local(pcols)
   real(r8), target, intent(inout) :: qc_local(pcols,pver), rliq_local(pcols), wtdlf_local(pcols,pver,wtrc_nwset), freqsh_local(pcols)
   real(r8), target, intent(inout) :: icwmr_local(pcols,pver), rprdsh_local(pcols,pver), ptend_q_local(pcols,pver,pcnst)
   real(r8), target, intent(inout) :: cmfsl_local(pcols,pverp), cmflq_local(pcols,pverp), rprdtot_local(pcols,pver), ftem_local(pcols,pver)

   integer(c_int64_t), target :: liq_type_c(wtrc_nwset), ice_type_c(wtrc_nwset)
   integer :: m

   interface
      subroutine convect_shallow_uw_post_shell_codon(ncol_c, pcols_c, pver_c, pverp_c, pcnst_c, wtrc_nwset_c, latvap_c, cpair_c, &
           state_pmid_p, cmfmc_p, cmfmc2_p, cnt_p, cnt2_p, cnb_p, cnb2_p, pcnt_p, pcnb_p, qc_p, qc2_p, rliq_p, rliq2_p, &
           wtqc_p, wtdlf_p, freqsh_p, icwmr_p, iccmr_uw_p, rprdsh_p, cmfdqs_p, ptend_q_p, ptend_tracer_p, cmfsl_p, cmflq_p, &
           slflx_p, qtflx_p, rprddp_p, rprdtot_p, ptend_s_p, ftem_p, liq_type_p, ice_type_p) &
           bind(c, name="convect_shallow_uw_post_shell_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c, pcnst_c, wtrc_nwset_c
         real(c_double), value :: latvap_c, cpair_c
         type(c_ptr), value :: state_pmid_p, cmfmc_p, cmfmc2_p, cnt_p, cnt2_p, cnb_p, cnb2_p, pcnt_p, pcnb_p, qc_p, qc2_p
         type(c_ptr), value :: rliq_p, rliq2_p, wtqc_p, wtdlf_p, freqsh_p, icwmr_p, iccmr_uw_p, rprdsh_p, cmfdqs_p
         type(c_ptr), value :: ptend_q_p, ptend_tracer_p, cmfsl_p, cmflq_p, slflx_p, qtflx_p, rprddp_p, rprdtot_p
         type(c_ptr), value :: ptend_s_p, ftem_p, liq_type_p, ice_type_p
      end subroutine convect_shallow_uw_post_shell_codon
   end interface

   do m = 1, wtrc_nwset
      liq_type_c(m) = int(wtrc_iatype(m,iwtliq), c_int64_t)
      ice_type_c(m) = int(wtrc_iatype(m,iwtice), c_int64_t)
   end do

   call convect_shallow_log_uw_post_shell_entered()
   call convect_shallow_outer_stage_dispatch_call(5_c_int64_t, 0_c_int64_t, int(ncol_local, c_int64_t), &
        int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), int(pcnst, c_int64_t), &
        int(wtrc_nwset, c_int64_t), 0_c_int64_t, 0_c_int64_t, 0_c_int64_t, 0._r8, real(latvap, c_double), &
        0._r8, 0._r8, real(cpair, c_double), c_loc(state_pmid_local), c_loc(cmfmc_local), c_loc(cmfmc2_local), &
        c_loc(cnt_local), c_loc(cnt2_local), c_loc(cnb_local), c_loc(cnb2_local), c_loc(pcnt_local), c_loc(pcnb_local), &
        c_loc(qc_local), c_loc(qc2_local), c_loc(rliq_local), c_loc(rliq2_local), c_loc(wtqc_local), c_loc(wtdlf_local), &
        c_loc(freqsh_local), c_loc(icwmr_local), c_loc(iccmr_uw_local), c_loc(rprdsh_local), c_loc(cmfdqs_local), &
        c_loc(ptend_q_local), c_loc(ptend_tracer_local), c_loc(cmfsl_local), c_loc(cmflq_local), c_loc(slflx_local), &
        c_loc(qtflx_local), c_loc(rprddp_local), c_loc(rprdtot_local), c_loc(ptend_s_local), c_loc(ftem_local), &
        c_loc(liq_type_c), c_loc(ice_type_c))

end subroutine convect_shallow_uw_post_shell

subroutine convect_shallow_postmerge(ncol_local, state_pmid_local, cmfmc_local, cmfmc2_local, cnt_local, cnt2_local, cnb_local, &
     cnb2_local, pcnt_local, pcnb_local, qc_local, qc2_local, rliq_local, rliq2_local, wtqc_local, wtdlf_local, freqsh_local)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr
   use constituents, only: pcnst
   use water_tracer_vars, only: wtrc_nwset, wtrc_iatype
   use water_types, only: iwtliq, iwtice

   integer, intent(in) :: ncol_local
   real(r8), target, intent(in) :: state_pmid_local(pcols,pver), cmfmc2_local(pcols,pverp), cnt2_local(pcols), cnb2_local(pcols)
   real(r8), target, intent(in) :: qc2_local(pcols,pver), rliq2_local(pcols), wtqc_local(pcols,pver,pcnst)
   real(r8), target, intent(inout) :: cmfmc_local(pcols,pverp), cnt_local(pcols), cnb_local(pcols), pcnt_local(pcols), pcnb_local(pcols)
   real(r8), target, intent(inout) :: qc_local(pcols,pver), rliq_local(pcols), wtdlf_local(pcols,pver,wtrc_nwset), freqsh_local(pcols)

   integer(c_int64_t), target :: liq_type_c(wtrc_nwset), ice_type_c(wtrc_nwset)
   integer :: m

   interface
      subroutine convect_shallow_postmerge_codon(ncol_c, pcols_c, pver_c, pverp_c, pcnst_c, wtrc_nwset_c, state_pmid_p, cmfmc_p, &
           cmfmc2_p, cnt_p, cnt2_p, cnb_p, cnb2_p, pcnt_p, pcnb_p, qc_p, qc2_p, rliq_p, rliq2_p, wtqc_p, wtdlf_p, freqsh_p, &
           liq_type_p, ice_type_p) bind(c, name="convect_shallow_postmerge_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c, pcnst_c, wtrc_nwset_c
         type(c_ptr), value :: state_pmid_p, cmfmc_p, cmfmc2_p, cnt_p, cnt2_p, cnb_p, cnb2_p, pcnt_p, pcnb_p, qc_p, qc2_p
         type(c_ptr), value :: rliq_p, rliq2_p, wtqc_p, wtdlf_p, freqsh_p, liq_type_p, ice_type_p
      end subroutine convect_shallow_postmerge_codon
   end interface

   do m = 1, wtrc_nwset
      liq_type_c(m) = int(wtrc_iatype(m,iwtliq), c_int64_t)
      ice_type_c(m) = int(wtrc_iatype(m,iwtice), c_int64_t)
   end do

   call convect_shallow_postmerge_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(pverp, c_int64_t), int(pcnst, c_int64_t), int(wtrc_nwset, c_int64_t), c_loc(state_pmid_local), c_loc(cmfmc_local), &
        c_loc(cmfmc2_local), c_loc(cnt_local), c_loc(cnt2_local), c_loc(cnb_local), c_loc(cnb2_local), c_loc(pcnt_local), &
        c_loc(pcnb_local), c_loc(qc_local), c_loc(qc2_local), c_loc(rliq_local), c_loc(rliq2_local), c_loc(wtqc_local), &
        c_loc(wtdlf_local), c_loc(freqsh_local), c_loc(liq_type_c), c_loc(ice_type_c))

end subroutine convect_shallow_postmerge

  end module convect_shallow
