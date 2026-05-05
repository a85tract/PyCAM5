module micro_mg_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for MG microphysics
!
!---------------------------------------------------------------------------------
!
! How to add new packed MG inputs to micro_mg_cam_tend:
!
! If you have an input with first dimension [psetcols, pver], the procedure
! for adding inputs is as follows:
!
! 1) In addition to any variables you need to declare for the "unpacked"
!    (CAM format) version, you must declare an allocatable or pointer array
!    for the "packed" (MG format) version.
!
! 2) After micro_mg_get_cols is called, allocate the "packed" array with
!    size [mgncol, nlev].
!
! 3) Add a call similar to the following line (look before the
!    micro_mg_tend calls to see similar lines):
!
!      packed_array = packer%pack(original_array)
!
!    The packed array can then be passed into any of the MG schemes.
!
! This same procedure will also work for 1D arrays of size psetcols, 3-D
! arrays with psetcols and pver as the first dimensions, and for arrays of
! dimension [psetcols, pverp]. You only have to modify the allocation of
! the packed array before the "pack" call.
!
!---------------------------------------------------------------------------------
!
! How to add new packed MG outputs to micro_mg_cam_tend:
!
! 1) As with inputs, in addition to the unpacked outputs you must declare
!    an allocatable or pointer array for packed data. The unpacked and
!    packed arrays must *also* be targets or pointers (but cannot be both).
!
! 2) Again as for inputs, allocate the packed array using mgncol and nlev,
!    which are set in micro_mg_get_cols.
!
! 3) Add the field to post-processing as in the following line (again,
!    there are many examples before the micro_mg_tend calls):
!
!      call post_proc%add_field(p(final_array),p(packed_array))
!
!    This registers the field for post-MG averaging, and to scatter to the
!    final, unpacked version of the array.
!
!    By default, any columns/levels that are not operated on by MG will be
!    set to 0 on output; this value can be adjusted using the "fillvalue"
!    optional argument to post_proc%add_field.
!
!    Also by default, outputs from multiple substeps will be averaged after
!    MG's substepping is complete. Passing the optional argument
!    "accum_method=accum_null" will change this behavior so that the last
!    substep is always output.
!
! This procedure works on 1-D and 2-D outputs. Note that the final,
! unpacked arrays are not set until the call to
! "post_proc%process_and_unpack", which sets every single field that was
! added with post_proc%add_field.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver, pverp, psubcols
use physconst,      only: gravit, rair, tmelt, cpair, rh2o, rhoh2o, &
                          latvap, latice, mwh2o
use phys_control,   only: phys_getopts, use_hetfrz_classnuc


use physics_types,  only: physics_state, physics_ptend, &
                          physics_ptend_init, physics_state_copy, &
                          physics_update, physics_state_dealloc, &
                          physics_ptend_sum, physics_ptend_scale

use physics_buffer, only: physics_buffer_desc, pbuf_add_field, dyn_time_lvls, &
                          pbuf_old_tim_idx, pbuf_get_index, dtype_r8, dtype_i4, &
                          pbuf_get_field, pbuf_set_field, col_type_subcol, &
                          pbuf_register_subcol
use constituents,   only: cnst_add, cnst_get_ind, &
                          cnst_name, cnst_longname, sflxnam, apcnst, bpcnst, pcnst

use cldfrc2m,       only: rhmini=>rhmini_const

use cam_history,    only: addfld, add_default, phys_decomp, outfld

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use error_messages, only: handle_errmsg
use ref_pres,       only: top_lev=>trop_cloud_top_lev

use subcol_utils,   only: subcol_get_scheme

implicit none
private
save

public :: &
   micro_mg_cam_readnl,          &
   micro_mg_cam_register,        &
   micro_mg_cam_init_cnst,       &
   micro_mg_cam_implements_cnst, &
   micro_mg_cam_init,            &
   micro_mg_cam_tend,            &
   micro_mg_version

integer :: micro_mg_version     = 1      ! Version number for MG.
integer :: micro_mg_sub_version = 0      ! Second part of version number.

real(r8) :: micro_mg_dcs = -1._r8

logical :: microp_uniform

character(len=16) :: micro_mg_precip_frac_method = 'max_overlap' ! type of precipitation fraction method

real(r8)          :: micro_mg_berg_eff_factor    = 1.0_r8        ! berg efficiency factor

logical, public :: do_cldliq ! Prognose cldliq flag
logical, public :: do_cldice ! Prognose cldice flag

logical :: use_native_postmg_diag_impl = .false.
logical :: postmg_diag_impl_selected = .false.
logical :: use_native_grid_diag_impl = .false.
logical :: grid_diag_impl_selected = .false.
logical :: use_native_tail_shell_impl = .false.
logical :: tail_shell_impl_selected = .false.
logical :: use_native_wtrc_shell_impl = .false.
logical :: wtrc_shell_impl_selected = .false.
logical :: use_native_wtrc_prep_impl = .false.
logical :: wtrc_prep_impl_selected = .false.
logical :: use_native_budget_diag_impl = .false.
logical :: budget_diag_impl_selected = .false.
logical :: use_native_reff_calc_impl = .false.
logical :: reff_calc_impl_selected = .false.
logical :: use_native_diag_shell_impl = .false.
logical :: diag_shell_impl_selected = .false.
logical :: use_native_pbuf_copy_impl = .false.
logical :: pbuf_copy_impl_selected = .false.
logical :: pbuf_copy_entered_logged = .false.
logical :: use_reff_calc_compare = .false.
logical :: reff_calc_compare_selected = .false.
logical :: reff_calc_compare_done = .false.

integer :: num_steps ! Number of MG substeps

integer :: ncnst = 4       ! Number of constituents

character(len=8), parameter :: &      ! Constituent names
   cnst_names(8) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE', &
                     'RAINQM', 'SNOWQM','NUMRAI','NUMSNO'/)

integer :: &
   ixcldliq = -1,      &! cloud liquid amount index
   ixcldice = -1,      &! cloud ice amount index
   ixnumliq = -1,      &! cloud liquid number index
   ixnumice = -1,      &! cloud ice water index
   ixrain = -1,        &! rain index
   ixsnow = -1,        &! snow index
   ixnumrain = -1,     &! rain number index
   ixnumsnow = -1       ! snow number index

! Physics buffer indices for fields registered by this module
integer :: &
   cldo_idx,           &
   qme_idx,            &
   prain_idx,          &
   nevapr_idx,         &
   wsedl_idx,          &
   rei_idx,            &
   rel_idx,            &
   dei_idx,            &
   mu_idx,             &
   prer_evap_idx,            &
   lambdac_idx,        &
   iciwpst_idx,        &
   iclwpst_idx,        &
   des_idx,            &
   icswp_idx,          &
   cldfsnow_idx,       &
   rate1_cw2pr_st_idx = -1, &
   ls_flxprc_idx,      &
   ls_flxsnw_idx,      &
   relvar_idx,         &
   cmeliq_idx,         &
   accre_enhan_idx

! Fields for UNICON
integer :: &
     am_evp_st_idx,      &! Evaporation area of stratiform precipitation
     evprain_st_idx,     &! Evaporation rate of stratiform rain [kg/kg/s]. >= 0.
     evpsnow_st_idx       ! Evaporation rate of stratiform snow [kg/kg/s]. >= 0.

! Fields needed as inputs to COSP
integer :: &
     ls_mrprc_idx,    ls_mrsnw_idx,    &
     ls_reffrain_idx, ls_reffsnow_idx, &
     cv_reffliq_idx,  cv_reffice_idx

! Fields needed by Park macrophysics
integer :: &
     cc_t_idx,  cc_qv_idx, &
     cc_ql_idx, cc_qi_idx, &
     cc_nl_idx, cc_ni_idx, &
     cc_qlst_idx

! Used to replace aspects of MG microphysics
! (e.g. by CARMA)
integer :: &
     tnd_qsnow_idx = -1, &
     tnd_nsnow_idx = -1, &
     re_ice_idx = -1

! Index fields for precipitation efficiency.
integer :: &
     acpr_idx = -1, &
     acgcme_idx = -1, &
     acnum_idx = -1

! Physics buffer indices for fields registered by other modules
integer :: &
   ast_idx = -1,            &
   cld_idx = -1,            &
   concld_idx = -1

! Pbuf fields needed for subcol_SILHS
integer :: &
     qrain_idx=-1, qsnow_idx=-1,    &
     nrain_idx=-1, nsnow_idx=-1

integer :: &
   naai_idx = -1,           &
   naai_hom_idx = -1,       &
   npccn_idx = -1,          &
   rndst_idx = -1,          &
   nacon_idx = -1,          &
   prec_str_idx = -1,       &
   snow_str_idx = -1,       &
   prec_pcw_idx = -1,       &
   snow_pcw_idx = -1,       &
   prec_sed_idx = -1,       &
   snow_sed_idx = -1

! pbuf fields for heterogeneous freezing
integer :: &
   frzimm_idx = -1, &
   frzcnt_idx = -1, &
   frzdep_idx = -1

interface p
   module procedure p1
   module procedure p2
end interface p


!===============================================================================
contains
!===============================================================================

subroutine micro_mg_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Namelist variables
  logical :: micro_mg_do_cldice = .true. ! do_cldice = .true., MG microphysics is prognosing cldice
  logical :: micro_mg_do_cldliq = .true. ! do_cldliq = .true., MG microphysics is prognosing cldliq
  integer :: micro_mg_num_steps = 1      ! Number of substepping iterations done by MG (1.5 only for now).


  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'micro_mg_cam_readnl'

  namelist /micro_mg_nl/ micro_mg_version, micro_mg_sub_version, &
       micro_mg_do_cldice, micro_mg_do_cldliq, micro_mg_num_steps, &
       microp_uniform, micro_mg_dcs, micro_mg_precip_frac_method, micro_mg_berg_eff_factor

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'micro_mg_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, micro_mg_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

     ! set local variables
     do_cldice = micro_mg_do_cldice
     do_cldliq = micro_mg_do_cldliq
     num_steps = micro_mg_num_steps

     ! Verify that version numbers are valid.
     select case (micro_mg_version)
     case (1)
        select case (micro_mg_sub_version)
        case(0)
           ! MG version 1.0
        case(5)
           ! MG version 1.5 - MG2 development
        case default
           call bad_version_endrun()
        end select
     case (2)
        select case (micro_mg_sub_version)
        case(0)
           ! MG version 2.0
        case default
           call bad_version_endrun()
        end select
     case default
        call bad_version_endrun()
     end select

     if (micro_mg_dcs < 0._r8) call endrun( "micro_mg_cam_readnl: &
              &micro_mg_dcs has not been set to a valid value.")
  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(micro_mg_version,            1, mpiint, 0, mpicom)
  call mpibcast(micro_mg_sub_version,        1, mpiint, 0, mpicom)
  call mpibcast(do_cldice,                   1, mpilog, 0, mpicom)
  call mpibcast(do_cldliq,                   1, mpilog, 0, mpicom)
  call mpibcast(num_steps,                   1, mpiint, 0, mpicom)
  call mpibcast(microp_uniform,              1, mpilog, 0, mpicom)
  call mpibcast(micro_mg_dcs,                1, mpir8,  0, mpicom)
  call mpibcast(micro_mg_berg_eff_factor,    1, mpir8,  0, mpicom)
  call mpibcast(micro_mg_precip_frac_method, 16, mpichar,0, mpicom)

#endif

contains

  subroutine bad_version_endrun
    ! Endrun wrapper with a more useful error message.
    character(len=128) :: errstring
    write(errstring,*) "Invalid version number specified for MG microphysics: ", &
         micro_mg_version,".",micro_mg_sub_version
    call endrun(errstring)
  end subroutine bad_version_endrun

end subroutine micro_mg_cam_readnl

!================================================================================================

subroutine micro_mg_cam_register

   ! Register microphysics constituents and fields in the physics buffer.
   !-----------------------------------------------------------------------

   logical :: prog_modal_aero
   logical :: use_subcol_microp  ! If true, then are using subcolumns in microphysics

   call phys_getopts(use_subcol_microp_out    = use_subcol_microp, &
                     prog_modal_aero_out      = prog_modal_aero)

   ! Register microphysics constituents and save indices.

   call cnst_add(cnst_names(1), mwh2o, cpair, 0._r8, ixcldliq, &
      longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
   call cnst_add(cnst_names(2), mwh2o, cpair, 0._r8, ixcldice, &
      longname='Grid box averaged cloud ice amount', is_convtran1=.true.)

   ! The next statements should have "is_convtran1=.true.", but this would change
   ! answers for MG 1.0. Thus make an exception for that version only.
   if (micro_mg_version == 1 .and. micro_mg_sub_version == 0) then
      call cnst_add(cnst_names(3), mwh2o, cpair, 0._r8, ixnumliq, &
           longname='Grid box averaged cloud liquid number', is_convtran1=.false.)
      call cnst_add(cnst_names(4), mwh2o, cpair, 0._r8, ixnumice, &
           longname='Grid box averaged cloud ice number', is_convtran1=.false.)
   else
      call cnst_add(cnst_names(3), mwh2o, cpair, 0._r8, ixnumliq, &
           longname='Grid box averaged cloud liquid number', is_convtran1=.true.)
      call cnst_add(cnst_names(4), mwh2o, cpair, 0._r8, ixnumice, &
           longname='Grid box averaged cloud ice number', is_convtran1=.true.)
   end if

   ! Note is_convtran1 is set to .true.
   if (micro_mg_version > 1) then
      call cnst_add(cnst_names(5), mwh2o, cpair, 0._r8, ixrain, &
           longname='Grid box averaged rain amount', is_convtran1=.true.)
      call cnst_add(cnst_names(6), mwh2o, cpair, 0._r8, ixsnow, &
           longname='Grid box averaged snow amount', is_convtran1=.true.)
      call cnst_add(cnst_names(7), mwh2o, cpair, 0._r8, ixnumrain, &
           longname='Grid box averaged rain number', is_convtran1=.true.)
      call cnst_add(cnst_names(8), mwh2o, cpair, 0._r8, ixnumsnow, &
           longname='Grid box averaged snow number', is_convtran1=.true.)
   end if

   ! Request physics buffer space for fields that persist across timesteps.

   call pbuf_add_field('CLDO','global',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldo_idx)

   ! Physics buffer variables for convective cloud properties.

   call pbuf_add_field('QME',        'physpkg',dtype_r8,(/pcols,pver/), qme_idx)
   call pbuf_add_field('PRAIN',      'physpkg',dtype_r8,(/pcols,pver/), prain_idx)
   call pbuf_add_field('NEVAPR',     'physpkg',dtype_r8,(/pcols,pver/), nevapr_idx)
   call pbuf_add_field('PRER_EVAP',  'global', dtype_r8,(/pcols,pver/), prer_evap_idx)

   call pbuf_add_field('WSEDL',      'physpkg',dtype_r8,(/pcols,pver/), wsedl_idx)

   call pbuf_add_field('REI',        'physpkg',dtype_r8,(/pcols,pver/), rei_idx)
   call pbuf_add_field('REL',        'physpkg',dtype_r8,(/pcols,pver/), rel_idx)

   ! Mitchell ice effective diameter for radiation
   call pbuf_add_field('DEI',        'physpkg',dtype_r8,(/pcols,pver/), dei_idx)
   ! Size distribution shape parameter for radiation
   call pbuf_add_field('MU',         'physpkg',dtype_r8,(/pcols,pver/), mu_idx)
   ! Size distribution shape parameter for radiation
   call pbuf_add_field('LAMBDAC',    'physpkg',dtype_r8,(/pcols,pver/), lambdac_idx)

   ! Stratiform only in cloud ice water path for radiation
   call pbuf_add_field('ICIWPST',    'physpkg',dtype_r8,(/pcols,pver/), iciwpst_idx)
   ! Stratiform in cloud liquid water path for radiation
   call pbuf_add_field('ICLWPST',    'physpkg',dtype_r8,(/pcols,pver/), iclwpst_idx)

   ! Snow effective diameter for radiation
   call pbuf_add_field('DES',        'physpkg',dtype_r8,(/pcols,pver/), des_idx)
   ! In cloud snow water path for radiation
   call pbuf_add_field('ICSWP',      'physpkg',dtype_r8,(/pcols,pver/), icswp_idx)
   ! Cloud fraction for liquid drops + snow
   call pbuf_add_field('CLDFSNOW ',  'physpkg',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldfsnow_idx)

   if (prog_modal_aero) then
      call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,(/pcols,pver/), rate1_cw2pr_st_idx)
   endif

   call pbuf_add_field('LS_FLXPRC',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxprc_idx)
   call pbuf_add_field('LS_FLXSNW',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxsnw_idx)


   ! Fields needed as inputs to COSP
   call pbuf_add_field('LS_MRPRC',   'physpkg',dtype_r8,(/pcols,pver/), ls_mrprc_idx)
   call pbuf_add_field('LS_MRSNW',   'physpkg',dtype_r8,(/pcols,pver/), ls_mrsnw_idx)
   call pbuf_add_field('LS_REFFRAIN','physpkg',dtype_r8,(/pcols,pver/), ls_reffrain_idx)
   call pbuf_add_field('LS_REFFSNOW','physpkg',dtype_r8,(/pcols,pver/), ls_reffsnow_idx)
   call pbuf_add_field('CV_REFFLIQ', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffliq_idx)
   call pbuf_add_field('CV_REFFICE', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffice_idx)

   ! CC_* Fields needed by Park macrophysics
   call pbuf_add_field('CC_T',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_t_idx)
   call pbuf_add_field('CC_qv',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qv_idx)
   call pbuf_add_field('CC_ql',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ql_idx)
   call pbuf_add_field('CC_qi',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qi_idx)
   call pbuf_add_field('CC_nl',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_nl_idx)
   call pbuf_add_field('CC_ni',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ni_idx)
   call pbuf_add_field('CC_qlst',  'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qlst_idx)

   ! Fields for UNICON
   call pbuf_add_field('am_evp_st',  'global', dtype_r8, (/pcols,pver/), am_evp_st_idx)
   call pbuf_add_field('evprain_st', 'global', dtype_r8, (/pcols,pver/), evprain_st_idx)
   call pbuf_add_field('evpsnow_st', 'global', dtype_r8, (/pcols,pver/), evpsnow_st_idx)

   ! Register subcolumn pbuf fields
   if (use_subcol_microp) then
      ! Global pbuf fields
      call pbuf_register_subcol('CLDO',        'micro_mg_cam_register', cldo_idx)

      ! CC_* Fields needed by Park macrophysics
      call pbuf_register_subcol('CC_T',        'micro_mg_cam_register', cc_t_idx)
      call pbuf_register_subcol('CC_qv',       'micro_mg_cam_register', cc_qv_idx)
      call pbuf_register_subcol('CC_ql',       'micro_mg_cam_register', cc_ql_idx)
      call pbuf_register_subcol('CC_qi',       'micro_mg_cam_register', cc_qi_idx)
      call pbuf_register_subcol('CC_nl',       'micro_mg_cam_register', cc_nl_idx)
      call pbuf_register_subcol('CC_ni',       'micro_mg_cam_register', cc_ni_idx)
      call pbuf_register_subcol('CC_qlst',     'micro_mg_cam_register', cc_qlst_idx)

      ! Physpkg pbuf fields
      ! Physics buffer variables for convective cloud properties.

      call pbuf_register_subcol('QME',         'micro_mg_cam_register', qme_idx)
      call pbuf_register_subcol('PRAIN',       'micro_mg_cam_register', prain_idx)
      call pbuf_register_subcol('NEVAPR',      'micro_mg_cam_register', nevapr_idx)
      call pbuf_register_subcol('PRER_EVAP',   'micro_mg_cam_register', prer_evap_idx)

      call pbuf_register_subcol('WSEDL',       'micro_mg_cam_register', wsedl_idx)

      call pbuf_register_subcol('REI',         'micro_mg_cam_register', rei_idx)
      call pbuf_register_subcol('REL',         'micro_mg_cam_register', rel_idx)

      ! Mitchell ice effective diameter for radiation
      call pbuf_register_subcol('DEI',         'micro_mg_cam_register', dei_idx)
      ! Size distribution shape parameter for radiation
      call pbuf_register_subcol('MU',          'micro_mg_cam_register', mu_idx)
      ! Size distribution shape parameter for radiation
      call pbuf_register_subcol('LAMBDAC',     'micro_mg_cam_register', lambdac_idx)

      ! Stratiform only in cloud ice water path for radiation
      call pbuf_register_subcol('ICIWPST',     'micro_mg_cam_register', iciwpst_idx)
      ! Stratiform in cloud liquid water path for radiation
      call pbuf_register_subcol('ICLWPST',     'micro_mg_cam_register', iclwpst_idx)

      ! Snow effective diameter for radiation
      call pbuf_register_subcol('DES',         'micro_mg_cam_register', des_idx)
      ! In cloud snow water path for radiation
      call pbuf_register_subcol('ICSWP',       'micro_mg_cam_register', icswp_idx)
      ! Cloud fraction for liquid drops + snow
      call pbuf_register_subcol('CLDFSNOW ',   'micro_mg_cam_register', cldfsnow_idx)

      if (prog_modal_aero) then
         call pbuf_register_subcol('RATE1_CW2PR_ST', 'micro_mg_cam_register', rate1_cw2pr_st_idx)
      end if

      call pbuf_register_subcol('LS_FLXPRC',   'micro_mg_cam_register', ls_flxprc_idx)
      call pbuf_register_subcol('LS_FLXSNW',   'micro_mg_cam_register', ls_flxsnw_idx)

      ! Fields needed as inputs to COSP
      call pbuf_register_subcol('LS_MRPRC',    'micro_mg_cam_register', ls_mrprc_idx)
      call pbuf_register_subcol('LS_MRSNW',    'micro_mg_cam_register', ls_mrsnw_idx)
      call pbuf_register_subcol('LS_REFFRAIN', 'micro_mg_cam_register', ls_reffrain_idx)
      call pbuf_register_subcol('LS_REFFSNOW', 'micro_mg_cam_register', ls_reffsnow_idx)
      call pbuf_register_subcol('CV_REFFLIQ',  'micro_mg_cam_register', cv_reffliq_idx)
      call pbuf_register_subcol('CV_REFFICE',  'micro_mg_cam_register', cv_reffice_idx)
   end if

   ! Additional pbuf for CARMA interface
   if (.not. do_cldice) then
      call pbuf_add_field('TND_QSNOW',  'physpkg',dtype_r8,(/pcols,pver/), tnd_qsnow_idx)
      call pbuf_add_field('TND_NSNOW',  'physpkg',dtype_r8,(/pcols,pver/), tnd_nsnow_idx)
      call pbuf_add_field('RE_ICE',     'physpkg',dtype_r8,(/pcols,pver/), re_ice_idx)
   end if

   ! Precipitation efficiency fields across timesteps.
   call pbuf_add_field('ACPRECL',    'global',dtype_r8,(/pcols/), acpr_idx)   ! accumulated precip
   call pbuf_add_field('ACGCME',     'global',dtype_r8,(/pcols/), acgcme_idx) ! accumulated condensation
   call pbuf_add_field('ACNUM',      'global',dtype_i4,(/pcols/), acnum_idx)  ! counter for accumulated # timesteps

   ! SGS variability  -- These could be reset by CLUBB so they need to be grid only
   call pbuf_add_field('RELVAR',     'global',dtype_r8,(/pcols,pver/), relvar_idx)
   call pbuf_add_field('ACCRE_ENHAN','global',dtype_r8,(/pcols,pver/), accre_enhan_idx)

   ! Diagnostic fields needed for subcol_SILHS, need to be grid-only
   if (subcol_get_scheme() == 'SILHS') then
      call pbuf_add_field('QRAIN',   'global',dtype_r8,(/pcols,pver/), qrain_idx)
      call pbuf_add_field('QSNOW',   'global',dtype_r8,(/pcols,pver/), qsnow_idx)
      call pbuf_add_field('NRAIN',   'global',dtype_r8,(/pcols,pver/), nrain_idx)
      call pbuf_add_field('NSNOW',   'global',dtype_r8,(/pcols,pver/), nsnow_idx)
   end if

end subroutine micro_mg_cam_register

!===============================================================================

function micro_mg_cam_implements_cnst(name)

   ! Return true if specified constituent is implemented by the
   ! microphysics package

   character(len=*), intent(in) :: name        ! constituent name
   logical :: micro_mg_cam_implements_cnst    ! return value

   !-----------------------------------------------------------------------

   micro_mg_cam_implements_cnst = any(name == cnst_names)

end function micro_mg_cam_implements_cnst

!===============================================================================

subroutine micro_mg_cam_init_cnst(name, q, gcid)

   ! Initialize the microphysics constituents, if they are
   ! not read from the initial file.

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------

   if (micro_mg_cam_implements_cnst(name)) q = 0.0_r8

end subroutine micro_mg_cam_init_cnst

!===============================================================================

subroutine micro_mg_cam_init(pbuf2d)
   use time_manager,   only: is_first_step
   use micro_mg_utils, only: micro_mg_utils_init
   use micro_mg1_0, only: micro_mg_init1_0 => micro_mg_init
   use micro_mg1_5, only: micro_mg_init1_5 => micro_mg_init
   use micro_mg2_0, only: micro_mg_init2_0 => micro_mg_init

   !-----------------------------------------------------------------------
   !
   ! Initialization for MG microphysics
   !
   !-----------------------------------------------------------------------

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   integer :: m, mm
   logical :: history_amwg         ! output the variables used by the AMWG diag package
   logical :: history_budget       ! Output tendencies and state variables for CAM4
                                   ! temperature, water vapor, cloud ice and cloud
                                   ! liquid budgets.
   logical :: use_subcol_microp
   integer :: budget_histfile      ! output history file number for budget fields
   integer :: ierr
   character(128) :: errstring     ! return status (non-blank for error return)

   !-----------------------------------------------------------------------

   call phys_getopts(use_subcol_microp_out=use_subcol_microp)

   if (masterproc) then
      write(iulog,"(A,I2,A,I2)") "Initializing MG version ",micro_mg_version,".",micro_mg_sub_version
      if (.not. do_cldliq) &
           write(iulog,*) "MG prognostic cloud liquid has been turned off via namelist."
      if (.not. do_cldice) &
           write(iulog,*) "MG prognostic cloud ice has been turned off via namelist."
      write(iulog,*) "Number of microphysics substeps is: ",num_steps
   end if

   select case (micro_mg_version)
   case (1)
      ! Set constituent number for later loops.
      ncnst = 4

      select case (micro_mg_sub_version)
      case (0)
         ! MG 1 does not initialize micro_mg_utils, so have to do it here.
         call micro_mg_utils_init(r8, rh2o, cpair, tmelt, latvap, latice, &
              micro_mg_dcs, errstring)
         call handle_errmsg(errstring, subname="micro_mg_utils_init")

         call micro_mg_init1_0( &
              r8, gravit, rair, rh2o, cpair, &
              rhoh2o, tmelt, latvap, latice, &
              rhmini, micro_mg_dcs, use_hetfrz_classnuc, &
              micro_mg_precip_frac_method, micro_mg_berg_eff_factor, errstring)
      case (5)
         ! MG 1 does not initialize micro_mg_utils, so have to do it here.
         call micro_mg_utils_init(r8, rh2o, cpair, tmelt, latvap, latice, &
              micro_mg_dcs, errstring)
         call handle_errmsg(errstring, subname="micro_mg_utils_init")

         call micro_mg_init1_5( &
              r8, gravit, rair, rh2o, cpair, &
              tmelt, latvap, latice, rhmini, &
              micro_mg_dcs,                  &
              microp_uniform, do_cldice, use_hetfrz_classnuc, &
              micro_mg_precip_frac_method, micro_mg_berg_eff_factor, errstring)
      end select
   case (2)
      ! Set constituent number for later loops.
      ncnst = 8

      select case (micro_mg_sub_version)
      case (0)
         call micro_mg_init2_0( &
              r8, gravit, rair, rh2o, cpair, &
              tmelt, latvap, latice, rhmini, &
              micro_mg_dcs,                  &
              microp_uniform, do_cldice, use_hetfrz_classnuc, &
              micro_mg_precip_frac_method, micro_mg_berg_eff_factor, errstring)
      end select
   end select

   call handle_errmsg(errstring, subname="micro_mg_init")

   ! Register history variables
   do m = 1, ncnst
      call cnst_get_ind(cnst_names(m), mm)
      if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixsnow /)) ) then
         ! mass mixing ratios
         call addfld(cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm)                   , phys_decomp)
         call addfld(sflxnam(mm),   'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
      else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain, ixnumsnow /)) ) then
         ! number concentrations
         call addfld(cnst_name(mm), '1/kg    ', pver, 'A', cnst_longname(mm)                   , phys_decomp)
         call addfld(sflxnam(mm),   '1/m2/s  ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
      else
         call endrun( "micro_mg_cam_init: &
              &Could not call addfld for constituent with unknown units.")
      endif
   end do

   call addfld(apcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' after physics'  , phys_decomp)
   call addfld(apcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' after physics'  , phys_decomp)
   call addfld(bpcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' before physics' , phys_decomp)
   call addfld(bpcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' before physics' , phys_decomp)

   if (micro_mg_version > 1) then
      call addfld(apcnst(ixrain), 'kg/kg   ', pver, 'A', trim(cnst_name(ixrain))//' after physics'  , phys_decomp)
      call addfld(apcnst(ixsnow), 'kg/kg   ', pver, 'A', trim(cnst_name(ixsnow))//' after physics'  , phys_decomp)
      call addfld(bpcnst(ixrain), 'kg/kg   ', pver, 'A', trim(cnst_name(ixrain))//' before physics' , phys_decomp)
      call addfld(bpcnst(ixsnow), 'kg/kg   ', pver, 'A', trim(cnst_name(ixsnow))//' before physics' , phys_decomp)
   end if

   call addfld ('CME      ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap within the cloud'                      ,phys_decomp)
   call addfld ('PRODPREC ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of condensate to precip'              ,phys_decomp)
   call addfld ('EVAPPREC ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling precip'                   ,phys_decomp)
   call addfld ('EVAPSNOW ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling snow'                     ,phys_decomp)
   call addfld ('HPROGCLD ', 'W/kg'    , pver, 'A', 'Heating from prognostic clouds'                          ,phys_decomp)
   call addfld ('FICE     ', 'fraction', pver, 'A', 'Fractional ice content within cloud'                     ,phys_decomp)
   call addfld ('ICWMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus water mixing ratio'                ,phys_decomp)
   call addfld ('ICIMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus ice mixing ratio'                  ,phys_decomp)

   ! MG microphysics diagnostics
   call addfld ('QCSEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling cloud water'              ,phys_decomp)
   call addfld ('QISEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of sublimation of falling cloud ice'                ,phys_decomp)
   call addfld ('QVRES    ', 'kg/kg/s ', pver, 'A', 'Rate of residual condensation term'                      ,phys_decomp)
   call addfld ('CMEIOUT  ', 'kg/kg/s ', pver, 'A', 'Rate of deposition/sublimation of cloud ice'             ,phys_decomp)
   call addfld ('VTRMC    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud water fallspeed'                     ,phys_decomp)
   call addfld ('VTRMI    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud ice fallspeed'                       ,phys_decomp)
   call addfld ('QCSEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud water mixing ratio tendency from sedimentation'    ,phys_decomp)
   call addfld ('QISEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud ice mixing ratio tendency from sedimentation'      ,phys_decomp)
   call addfld ('PRAO     ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud water by rain'                        ,phys_decomp)
   call addfld ('PRCO     ', 'kg/kg/s ', pver, 'A', 'Autoconversion of cloud water'                           ,phys_decomp)
   call addfld ('MNUCCCO  ', 'kg/kg/s ', pver, 'A', 'Immersion freezing of cloud water'                       ,phys_decomp)
   call addfld ('MNUCCTO  ', 'kg/kg/s ', pver, 'A', 'Contact freezing of cloud water'                         ,phys_decomp)
   call addfld ('MNUCCDO  ', 'kg/kg/s ', pver, 'A', 'Homogeneous and heterogeneous nucleation from vapor'     ,phys_decomp)
   call addfld ('MNUCCDOhet','kg/kg/s ', pver, 'A', 'Heterogeneous nucleation from vapor'                     ,phys_decomp)
   call addfld ('MSACWIO  ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water from rime-splintering'         ,phys_decomp)
   call addfld ('PSACWSO  ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud water by snow'                        ,phys_decomp)
   call addfld ('BERGSO   ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water to snow from bergeron'         ,phys_decomp)
   call addfld ('BERGO    ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water to cloud ice from bergeron'    ,phys_decomp)
   call addfld ('MELTO    ', 'kg/kg/s ', pver, 'A', 'Melting of cloud ice'                                    ,phys_decomp)
   call addfld ('HOMOO    ', 'kg/kg/s ', pver, 'A', 'Homogeneous freezing of cloud water'                     ,phys_decomp)
   call addfld ('QCRESO   ', 'kg/kg/s ', pver, 'A', 'Residual condensation term for cloud water'              ,phys_decomp)
   call addfld ('PRCIO    ', 'kg/kg/s ', pver, 'A', 'Autoconversion of cloud ice'                             ,phys_decomp)
   call addfld ('PRAIO    ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud ice by rain'                          ,phys_decomp)
   call addfld ('QIRESO   ', 'kg/kg/s ', pver, 'A', 'Residual deposition term for cloud ice'                  ,phys_decomp)
   call addfld ('MNUCCRO  ', 'kg/kg/s ', pver, 'A', 'Heterogeneous freezing of rain to snow'                  ,phys_decomp)
   call addfld ('PRACSO   ', 'kg/kg/s ', pver, 'A', 'Accretion of rain by snow'                               ,phys_decomp)
   call addfld ('MELTSDT  ', 'W/kg    ', pver, 'A', 'Latent heating rate due to melting of snow'              ,phys_decomp)
   call addfld ('FRZRDT   ', 'W/kg    ', pver, 'A', 'Latent heating rate due to homogeneous freezing of rain' ,phys_decomp)
   if (micro_mg_version > 1) then
      call addfld ('QRSEDTEN ', 'kg/kg/s ', pver, 'A', 'Rain mixing ratio tendency from sedimentation'           ,phys_decomp)
      call addfld ('QSSEDTEN ', 'kg/kg/s ', pver, 'A', 'Snow mixing ratio tendency from sedimentation'           ,phys_decomp)
   end if

   ! History variables for CAM5 microphysics
   call addfld ('MPDT     ', 'W/kg    ', pver, 'A', 'Heating tendency - Morrison microphysics'                ,phys_decomp)
   call addfld ('MPDQ     ', 'kg/kg/s ', pver, 'A', 'Q tendency - Morrison microphysics'                      ,phys_decomp)
   call addfld ('MPDLIQ   ', 'kg/kg/s ', pver, 'A', 'CLDLIQ tendency - Morrison microphysics'                 ,phys_decomp)
   call addfld ('MPDICE   ', 'kg/kg/s ', pver, 'A', 'CLDICE tendency - Morrison microphysics'                 ,phys_decomp)
   call addfld ('MPDW2V   ', 'kg/kg/s ', pver, 'A', 'Water <--> Vapor tendency - Morrison microphysics'       ,phys_decomp)
   call addfld ('MPDW2I   ', 'kg/kg/s ', pver, 'A', 'Water <--> Ice tendency - Morrison microphysics'         ,phys_decomp)
   call addfld ('MPDW2P   ', 'kg/kg/s ', pver, 'A', 'Water <--> Precip tendency - Morrison microphysics'      ,phys_decomp)
   call addfld ('MPDI2V   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Vapor tendency - Morrison microphysics'         ,phys_decomp)
   call addfld ('MPDI2W   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Water tendency - Morrison microphysics'         ,phys_decomp)
   call addfld ('MPDI2P   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Precip tendency - Morrison microphysics'        ,phys_decomp)
   call addfld ('ICWNC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud water number conc'                   ,phys_decomp)
   call addfld ('ICINC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud ice number conc'                     ,phys_decomp)
   call addfld ('EFFLIQ_IND','Micron  ', pver, 'A', 'Prognostic droplet effective radius (indirect effect)'   ,phys_decomp)
   call addfld ('CDNUMC   ', '1/m2    ', 1,    'A', 'Vertically-integrated droplet concentration'             ,phys_decomp)
   call addfld ('MPICLWPI ', 'kg/m2   ', 1,    'A', 'Vertically-integrated &
        &in-cloud Initial Liquid WP (Before Micro)' ,phys_decomp)
   call addfld ('MPICIWPI ', 'kg/m2   ', 1,    'A', 'Vertically-integrated &
        &in-cloud Initial Ice WP (Before Micro)'    ,phys_decomp)

   ! This is provided as an example on how to write out subcolumn output
   ! NOTE -- only 'I' should be used for sub-column fields as subc-columns could shift from time-step to time-step
   if (use_subcol_microp) then
      call addfld('FICE_SCOL', 'fraction', psubcols*pver, 'I', &
           'Sub-column fractional ice content within cloud', phys_decomp, &
           mdimnames=(/'psubcols','lev     '/), flag_xyfill=.true., fill_value=1.e30_r8)
   end if

   ! Averaging for cloud particle number and size
   call addfld ('AWNC     ', 'm-3     ', pver, 'A', 'Average cloud water number conc'                         ,phys_decomp)
   call addfld ('AWNI     ', 'm-3     ', pver, 'A', 'Average cloud ice number conc'                           ,phys_decomp)
   call addfld ('AREL     ', 'Micron  ', pver, 'A', 'Average droplet effective radius'                        ,phys_decomp)
   call addfld ('AREI     ', 'Micron  ', pver, 'A', 'Average ice effective radius'                            ,phys_decomp)
   ! Frequency arrays for above
   call addfld ('FREQL    ', 'fraction', pver, 'A', 'Fractional occurrence of liquid'                          ,phys_decomp)
   call addfld ('FREQI    ', 'fraction', pver, 'A', 'Fractional occurrence of ice'                             ,phys_decomp)

   ! Average cloud top particle size and number (liq, ice) and frequency
   call addfld ('ACTREL   ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet effective radius'              ,phys_decomp)
   call addfld ('ACTREI   ', 'Micron  ', 1,    'A', 'Average Cloud Top ice effective radius'                  ,phys_decomp)
   call addfld ('ACTNL    ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet number'                        ,phys_decomp)
   call addfld ('ACTNI    ', 'Micron  ', 1,    'A', 'Average Cloud Top ice number'                            ,phys_decomp)

   call addfld ('FCTL     ', 'fraction', 1,    'A', 'Fractional occurrence of cloud top liquid'                ,phys_decomp)
   call addfld ('FCTI     ', 'fraction', 1,    'A', 'Fractional occurrence of cloud top ice'                   ,phys_decomp)

   call addfld ('LS_FLXPRC', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface rain+snow flux', phys_decomp)
   call addfld ('LS_FLXSNW', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface snow flux', phys_decomp)

   call addfld ('REL', 'micron', pver, 'A', 'MG REL stratiform cloud effective radius liquid', phys_decomp)
   call addfld ('REI', 'micron', pver, 'A', 'MG REI stratiform cloud effective radius ice', phys_decomp)
   call addfld ('LS_REFFRAIN', 'micron', pver, 'A', 'ls stratiform rain effective radius', phys_decomp)
   call addfld ('LS_REFFSNOW', 'micron', pver, 'A', 'ls stratiform snow effective radius', phys_decomp)
   call addfld ('CV_REFFLIQ', 'micron', pver, 'A', 'convective cloud liq effective radius', phys_decomp)
   call addfld ('CV_REFFICE', 'micron', pver, 'A', 'convective cloud ice effective radius', phys_decomp)

   ! diagnostic precip
   call addfld ('QRAIN   ','kg/kg   ',pver, 'A','Diagnostic grid-mean rain mixing ratio'         ,phys_decomp)
   call addfld ('QSNOW   ','kg/kg   ',pver, 'A','Diagnostic grid-mean snow mixing ratio'         ,phys_decomp)
   call addfld ('NRAIN   ','m-3     ',pver, 'A','Diagnostic grid-mean rain number conc'         ,phys_decomp)
   call addfld ('NSNOW   ','m-3     ',pver, 'A','Diagnostic grid-mean snow number conc'         ,phys_decomp)

   ! size of precip
   call addfld ('RERCLD   ','m      ',pver, 'A','Diagnostic effective radius of Liquid Cloud and Rain' ,phys_decomp)
   call addfld ('DSNOW   ','m       ',pver, 'A','Diagnostic grid-mean snow diameter'         ,phys_decomp)

   ! diagnostic radar reflectivity, cloud-averaged
   call addfld ('REFL  ','DBz  ',pver, 'A','94 GHz radar reflectivity'       ,phys_decomp)
   call addfld ('AREFL  ','DBz  ',pver, 'A','Average 94 GHz radar reflectivity'       ,phys_decomp)
   call addfld ('FREFL  ','fraction  ',pver, 'A','Fractional occurrence of radar reflectivity'       ,phys_decomp)

   call addfld ('CSRFL  ','DBz  ',pver, 'A','94 GHz radar reflectivity (CloudSat thresholds)'       ,phys_decomp)
   call addfld ('ACSRFL  ','DBz  ',pver, 'A','Average 94 GHz radar reflectivity (CloudSat thresholds)'       ,phys_decomp)
   call addfld ('FCSRFL  ','fraction  ',pver, 'A','Fractional occurrence of radar reflectivity (CloudSat thresholds)' &
        ,phys_decomp)

   call addfld ('AREFLZ ','mm^6/m^3 ',pver, 'A','Average 94 GHz radar reflectivity'       ,phys_decomp)

   ! Aerosol information
   call addfld ('NCAL    ','1/m3   ',pver, 'A','Number Concentation Activated for Liquid',phys_decomp)
   call addfld ('NCAI    ','1/m3   ',pver, 'A','Number Concentation Activated for Ice',phys_decomp)

   ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
   call addfld ('AQRAIN   ','kg/kg   ',pver, 'A','Average rain mixing ratio'         ,phys_decomp)
   call addfld ('AQSNOW   ','kg/kg   ',pver, 'A','Average snow mixing ratio'         ,phys_decomp)
   call addfld ('ANRAIN   ','m-3     ',pver, 'A','Average rain number conc'         ,phys_decomp)
   call addfld ('ANSNOW   ','m-3     ',pver, 'A','Average snow number conc'         ,phys_decomp)
   call addfld ('ADRAIN   ','Micron  ',pver, 'A','Average rain effective Diameter'         ,phys_decomp)
   call addfld ('ADSNOW   ','Micron  ',pver, 'A','Average snow effective Diameter'         ,phys_decomp)
   call addfld ('FREQR  ','fraction  ',pver, 'A','Fractional occurrence of rain'       ,phys_decomp)
   call addfld ('FREQS  ','fraction  ',pver, 'A','Fractional occurrence of snow'       ,phys_decomp)

   ! precipitation efficiency & other diagnostic fields
   call addfld('PE'    , '1',       1, 'A', 'Stratiform Precipitation Efficiency  (precip/cmeliq)',       phys_decomp )
   call addfld('APRL'  , 'm/s',     1, 'A', 'Average Stratiform Precip Rate over efficiency calculation', phys_decomp )
   call addfld('PEFRAC', '1',       1, 'A', 'Fraction of timesteps precip efficiency reported',           phys_decomp )
   call addfld('VPRCO' , 'kg/kg/s', 1, 'A', 'Vertical average of autoconversion rate',                    phys_decomp )
   call addfld('VPRAO' , 'kg/kg/s', 1, 'A', 'Vertical average of accretion rate',                         phys_decomp )
   call addfld('RACAU' , 'kg/kg/s', 1, 'A', 'Accretion/autoconversion ratio from vertical average',       phys_decomp )

   if (micro_mg_version > 1) then
      call addfld('UMR',   'm/s     ', pver, 'A', 'Mass-weighted rain  fallspeed'              , phys_decomp)
      call addfld('UMS',   'm/s     ', pver, 'A', 'Mass-weighted snow fallspeed'               , phys_decomp)
   end if

   ! qc limiter (only output in versions 1.5 and later)
   if (.not. (micro_mg_version == 1 .and. micro_mg_sub_version == 0)) then
      call addfld('QCRAT', 'fraction', pver, 'A', 'Qc Limiter: Fraction of qc tendency applied', phys_decomp)
   end if

   ! determine the add_default fields
   call phys_getopts(history_amwg_out           = history_amwg         , &
                     history_budget_out         = history_budget       , &
                     history_budget_histfile_num_out = budget_histfile)

   if (history_amwg) then
      call add_default ('FICE    ', 1, ' ')
      call add_default ('AQRAIN   ', 1, ' ')
      call add_default ('AQSNOW   ', 1, ' ')
      call add_default ('ANRAIN   ', 1, ' ')
      call add_default ('ANSNOW   ', 1, ' ')
      call add_default ('ADRAIN   ', 1, ' ')
      call add_default ('ADSNOW   ', 1, ' ')
      call add_default ('AREI     ', 1, ' ')
      call add_default ('AREL     ', 1, ' ')
      call add_default ('AWNC     ', 1, ' ')
      call add_default ('AWNI     ', 1, ' ')
      call add_default ('CDNUMC   ', 1, ' ')
      call add_default ('FREQR    ', 1, ' ')
      call add_default ('FREQS    ', 1, ' ')
      call add_default ('FREQL    ', 1, ' ')
      call add_default ('FREQI    ', 1, ' ')
      do m = 1, ncnst
         call cnst_get_ind(cnst_names(m), mm)
         call add_default(cnst_name(mm), 1, ' ')
         ! call add_default(sflxnam(mm),   1, ' ')
      end do
   end if

   if ( history_budget ) then
      call add_default ('EVAPSNOW ', budget_histfile, ' ')
      call add_default ('EVAPPREC ', budget_histfile, ' ')
      call add_default ('QVRES    ', budget_histfile, ' ')
      call add_default ('QISEVAP  ', budget_histfile, ' ')
      call add_default ('QCSEVAP  ', budget_histfile, ' ')
      call add_default ('QISEDTEN ', budget_histfile, ' ')
      call add_default ('QCSEDTEN ', budget_histfile, ' ')
      call add_default ('QIRESO   ', budget_histfile, ' ')
      call add_default ('QCRESO   ', budget_histfile, ' ')
      if (micro_mg_version > 1) then
         call add_default ('QRSEDTEN ', budget_histfile, ' ')
         call add_default ('QSSEDTEN ', budget_histfile, ' ')
      end if
      call add_default ('PSACWSO  ', budget_histfile, ' ')
      call add_default ('PRCO     ', budget_histfile, ' ')
      call add_default ('PRCIO    ', budget_histfile, ' ')
      call add_default ('PRAO     ', budget_histfile, ' ')
      call add_default ('PRAIO    ', budget_histfile, ' ')
      call add_default ('PRACSO   ', budget_histfile, ' ')
      call add_default ('MSACWIO  ', budget_histfile, ' ')
      call add_default ('MPDW2V   ', budget_histfile, ' ')
      call add_default ('MPDW2P   ', budget_histfile, ' ')
      call add_default ('MPDW2I   ', budget_histfile, ' ')
      call add_default ('MPDT     ', budget_histfile, ' ')
      call add_default ('MPDQ     ', budget_histfile, ' ')
      call add_default ('MPDLIQ   ', budget_histfile, ' ')
      call add_default ('MPDICE   ', budget_histfile, ' ')
      call add_default ('MPDI2W   ', budget_histfile, ' ')
      call add_default ('MPDI2V   ', budget_histfile, ' ')
      call add_default ('MPDI2P   ', budget_histfile, ' ')
      call add_default ('MNUCCTO  ', budget_histfile, ' ')
      call add_default ('MNUCCRO  ', budget_histfile, ' ')
      call add_default ('MNUCCCO  ', budget_histfile, ' ')
      call add_default ('MELTSDT  ', budget_histfile, ' ')
      call add_default ('MELTO    ', budget_histfile, ' ')
      call add_default ('HOMOO    ', budget_histfile, ' ')
      call add_default ('FRZRDT   ', budget_histfile, ' ')
      call add_default ('CMEIOUT  ', budget_histfile, ' ')
      call add_default ('BERGSO   ', budget_histfile, ' ')
      call add_default ('BERGO    ', budget_histfile, ' ')

      call add_default(cnst_name(ixcldliq), budget_histfile, ' ')
      call add_default(cnst_name(ixcldice), budget_histfile, ' ')
      call add_default(apcnst   (ixcldliq), budget_histfile, ' ')
      call add_default(apcnst   (ixcldice), budget_histfile, ' ')
      call add_default(bpcnst   (ixcldliq), budget_histfile, ' ')
      call add_default(bpcnst   (ixcldice), budget_histfile, ' ')
      if (micro_mg_version > 1) then
         call add_default(cnst_name(ixrain), budget_histfile, ' ')
         call add_default(cnst_name(ixsnow), budget_histfile, ' ')
         call add_default(apcnst   (ixrain), budget_histfile, ' ')
         call add_default(apcnst   (ixsnow), budget_histfile, ' ')
         call add_default(bpcnst   (ixrain), budget_histfile, ' ')
         call add_default(bpcnst   (ixsnow), budget_histfile, ' ')
      end if

   end if

   ! physics buffer indices
   ast_idx      = pbuf_get_index('AST')
   cld_idx      = pbuf_get_index('CLD')
   concld_idx   = pbuf_get_index('CONCLD')

   naai_idx     = pbuf_get_index('NAAI')
   naai_hom_idx = pbuf_get_index('NAAI_HOM')
   npccn_idx    = pbuf_get_index('NPCCN')
   rndst_idx    = pbuf_get_index('RNDST')
   nacon_idx    = pbuf_get_index('NACON')

   prec_str_idx = pbuf_get_index('PREC_STR')
   snow_str_idx = pbuf_get_index('SNOW_STR')
   prec_sed_idx = pbuf_get_index('PREC_SED')
   snow_sed_idx = pbuf_get_index('SNOW_SED')
   prec_pcw_idx = pbuf_get_index('PREC_PCW')
   snow_pcw_idx = pbuf_get_index('SNOW_PCW')

   cmeliq_idx = pbuf_get_index('CMELIQ')

   ! These fields may have been added, so don't abort if they have not been
   qrain_idx    = pbuf_get_index('QRAIN', ierr)
   qsnow_idx    = pbuf_get_index('QSNOW', ierr)
   nrain_idx    = pbuf_get_index('NRAIN', ierr)
   nsnow_idx    = pbuf_get_index('NSNOW', ierr)

  ! fields for heterogeneous freezing
  frzimm_idx = pbuf_get_index('FRZIMM', ierr)
  frzcnt_idx = pbuf_get_index('FRZCNT', ierr)
  frzdep_idx = pbuf_get_index('FRZDEP', ierr)

  ! Initialize physics buffer grid fields for accumulating precip and condensation
   if (is_first_step()) then
      call pbuf_set_field(pbuf2d, cldo_idx,   0._r8)
      call pbuf_set_field(pbuf2d, cc_t_idx,   0._r8)
      call pbuf_set_field(pbuf2d, cc_qv_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_ql_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_qi_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_nl_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_ni_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_qlst_idx,0._r8)
      call pbuf_set_field(pbuf2d, acpr_idx,   0._r8)
      call pbuf_set_field(pbuf2d, acgcme_idx, 0._r8)
      call pbuf_set_field(pbuf2d, acnum_idx,  0)
      call pbuf_set_field(pbuf2d, relvar_idx, 2._r8)
      call pbuf_set_field(pbuf2d, accre_enhan_idx, 1._r8)
      call pbuf_set_field(pbuf2d, am_evp_st_idx,  0._r8)
      call pbuf_set_field(pbuf2d, evprain_st_idx, 0._r8)
      call pbuf_set_field(pbuf2d, evpsnow_st_idx, 0._r8)
      call pbuf_set_field(pbuf2d, prer_evap_idx,  0._r8)

      if (qrain_idx > 0)   call pbuf_set_field(pbuf2d, qrain_idx, 0._r8)
      if (qsnow_idx > 0)   call pbuf_set_field(pbuf2d, qsnow_idx, 0._r8)
      if (nrain_idx > 0)   call pbuf_set_field(pbuf2d, nrain_idx, 0._r8)
      if (nsnow_idx > 0)   call pbuf_set_field(pbuf2d, nsnow_idx, 0._r8)

      ! If sub-columns turned on, need to set the sub-column fields as well
      if (use_subcol_microp) then
         call pbuf_set_field(pbuf2d, cldo_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_t_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_qv_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_ql_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_qi_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_nl_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_ni_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_qlst_idx,0._r8, col_type=col_type_subcol)
      end if

   end if

end subroutine micro_mg_cam_init

!===============================================================================

subroutine micro_mg_cam_tend(state, ptend, dtime, pbuf)

   use micro_mg_utils, only: size_dist_param_basic, size_dist_param_liq, &
        mg_liq_props, mg_ice_props, avg_diameter, rhoi, rhosn, rhow, rhows, &
        qsmall, mincld

   use micro_mg_data, only: MGPacker, MGPostProc, accum_null, accum_mean

   use micro_mg1_0, only: micro_mg_tend1_0 => micro_mg_tend, &
        micro_mg_get_cols1_0 => micro_mg_get_cols
   use micro_mg1_5, only: micro_mg_tend1_5 => micro_mg_tend, &
        micro_mg_get_cols1_5 => micro_mg_get_cols
   use micro_mg2_0, only: micro_mg_tend2_0 => micro_mg_tend, &
        micro_mg_get_cols2_0 => micro_mg_get_cols

   use physics_buffer,  only: pbuf_col_type_index
   use subcol,          only: subcol_field_avg

   use water_tracer_vars, only: trace_water, wtrc_add_stprecip, wtrc_bulk_indices, &
                                wtrc_iatype, wtrc_indices, wtrc_ncnst
   use water_tracers,     only: wtrc_apply_rates, wtrc_init_rates, wtrc_add_rate, &
                                wtrc_output_precip
   use water_types,       only: pwtype, iwtvap, iwtliq, iwtice, iwtstrain, iwtstsnow

   type(physics_state),         intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: dtime
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! Local variables
   integer :: lchnk, ncol, psetcols, ngrdcol

   integer :: i, k, m, itim_old, it

   real(r8), pointer :: naai(:,:)      ! ice nucleation number
   real(r8), pointer :: naai_hom(:,:)  ! ice nucleation number (homogeneous)
   real(r8), pointer :: npccn(:,:)     ! liquid activation number tendency
   real(r8), pointer :: rndst(:,:,:)
   real(r8), pointer :: nacon(:,:,:)
   real(r8), pointer :: am_evp_st_grid(:,:)    ! Evaporation area of stratiform precipitation. 0<= am_evp_st <=1.
   real(r8), pointer :: evprain_st_grid(:,:)   ! Evaporation rate of stratiform rain [kg/kg/s]
   real(r8), pointer :: evpsnow_st_grid(:,:)   ! Evaporation rate of stratiform snow [kg/kg/s]

   real(r8), pointer :: prec_str(:)          ! [Total] Sfc flux of precip from stratiform [ m/s ]
   real(r8), pointer :: snow_str(:)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
   real(r8), pointer :: prec_sed(:)          ! Surface flux of total cloud water from sedimentation
   real(r8), pointer :: snow_sed(:)          ! Surface flux of cloud ice from sedimentation
   real(r8), pointer :: prec_pcw(:)          ! Sfc flux of precip from microphysics [ m/s ]
   real(r8), pointer :: snow_pcw(:)          ! Sfc flux of snow from microphysics [ m/s ]

   real(r8), pointer :: ast(:,:)          ! Relative humidity cloud fraction
   real(r8), pointer :: alst_mic(:,:)
   real(r8), pointer :: aist_mic(:,:)
   real(r8), pointer :: cldo(:,:)         ! Old cloud fraction
   real(r8), pointer :: nevapr(:,:)       ! Evaporation of total precipitation (rain + snow)
   real(r8), pointer :: prer_evap(:,:)    ! precipitation evaporation rate
   real(r8), pointer :: relvar(:,:)       ! relative variance of cloud water
   real(r8), pointer :: accre_enhan(:,:)  ! optional accretion enhancement for experimentation
   real(r8), pointer :: prain(:,:)        ! Total precipitation (rain + snow)
   real(r8), pointer :: dei(:,:)          ! Ice effective diameter (meters) (AG: microns?)
   real(r8), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
   real(r8), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation
   real(r8), pointer :: des(:,:)          ! Snow effective diameter (m)

   real(r8) :: rho(state%psetcols,pver)
   real(r8) :: cldmax(state%psetcols,pver)

   real(r8), target :: rate1cld(state%psetcols,pver) ! array to hold rate1ord_cw2pr_st from microphysics

   real(r8), target :: tlat(state%psetcols,pver)
   real(r8), target :: qvlat(state%psetcols,pver)
   real(r8), target :: qcten(state%psetcols,pver)
   real(r8), target :: qiten(state%psetcols,pver)
   real(r8), target :: ncten(state%psetcols,pver)
   real(r8), target :: niten(state%psetcols,pver)

   real(r8), target :: qrten(state%psetcols,pver)
   real(r8), target :: qsten(state%psetcols,pver)
   real(r8), target :: nrten(state%psetcols,pver)
   real(r8), target :: nsten(state%psetcols,pver)

   real(r8), target :: prect(state%psetcols)
   real(r8), target :: preci(state%psetcols)
   real(r8), target :: am_evp_st(state%psetcols,pver)  ! Area over which precip evaporates
   real(r8), target :: evapsnow(state%psetcols,pver)   ! Local evaporation of snow
   real(r8), target :: prodsnow(state%psetcols,pver)   ! Local production of snow
   real(r8), target :: cmeice(state%psetcols,pver)     ! Rate of cond-evap of ice within the cloud
   real(r8), target :: qsout(state%psetcols,pver)      ! Snow mixing ratio
   real(r8), target :: rflx(state%psetcols,pverp)      ! grid-box average rain flux (kg m^-2 s^-1)
   real(r8), target :: sflx(state%psetcols,pverp)      ! grid-box average snow flux (kg m^-2 s^-1)
   real(r8), target :: qrout(state%psetcols,pver)      ! Rain mixing ratio
   real(r8), target :: qcsevap(state%psetcols,pver)    ! Evaporation of falling cloud water
   real(r8), target :: qisevap(state%psetcols,pver)    ! Sublimation of falling cloud ice
   real(r8), target :: qvres(state%psetcols,pver)      ! Residual condensation term to remove excess saturation
   real(r8), target :: cmeiout(state%psetcols,pver)    ! Deposition/sublimation rate of cloud ice
   real(r8), target :: vtrmc(state%psetcols,pver)      ! Mass-weighted cloud water fallspeed
   real(r8), target :: vtrmi(state%psetcols,pver)      ! Mass-weighted cloud ice fallspeed
   real(r8), target :: umr(state%psetcols,pver)        ! Mass-weighted rain fallspeed
   real(r8), target :: ums(state%psetcols,pver)        ! Mass-weighted snow fallspeed
   real(r8), target :: qcsedten(state%psetcols,pver)   ! Cloud water mixing ratio tendency from sedimentation
   real(r8), target :: qisedten(state%psetcols,pver)   ! Cloud ice mixing ratio tendency from sedimentation
   real(r8), target :: qrsedten(state%psetcols,pver)   ! Rain mixing ratio tendency from sedimentation
   real(r8), target :: qssedten(state%psetcols,pver)   ! Snow mixing ratio tendency from sedimentation

   real(r8), target :: prao(state%psetcols,pver)
   real(r8), target :: prco(state%psetcols,pver)
   real(r8), target :: mnuccco(state%psetcols,pver)
   real(r8), target :: mnuccto(state%psetcols,pver)
   real(r8), target :: msacwio(state%psetcols,pver)
   real(r8), target :: psacwso(state%psetcols,pver)
   real(r8), target :: bergso(state%psetcols,pver)
   real(r8), target :: bergo(state%psetcols,pver)
   real(r8), target :: melto(state%psetcols,pver)
   real(r8), target :: homoo(state%psetcols,pver)
   real(r8), target :: qcreso(state%psetcols,pver)
   real(r8), target :: prcio(state%psetcols,pver)
   real(r8), target :: praio(state%psetcols,pver)
   real(r8), target :: qireso(state%psetcols,pver)
   real(r8), target :: mnuccro(state%psetcols,pver)
   real(r8), target :: pracso (state%psetcols,pver)
   real(r8), target :: meltsdt(state%psetcols,pver)
   real(r8), target :: frzrdt (state%psetcols,pver)
   real(r8), target :: mnuccdo(state%psetcols,pver)
   real(r8), target :: nrout(state%psetcols,pver)
   real(r8), target :: nsout(state%psetcols,pver)
   real(r8), target :: refl(state%psetcols,pver)    ! analytic radar reflectivity
   real(r8), target :: arefl(state%psetcols,pver)   ! average reflectivity will zero points outside valid range
   real(r8), target :: areflz(state%psetcols,pver)  ! average reflectivity in z.
   real(r8), target :: frefl(state%psetcols,pver)
   real(r8), target :: csrfl(state%psetcols,pver)   ! cloudsat reflectivity
   real(r8), target :: acsrfl(state%psetcols,pver)  ! cloudsat average
   real(r8), target :: fcsrfl(state%psetcols,pver)
   real(r8), target :: rercld(state%psetcols,pver)  ! effective radius calculation for rain + cloud
   real(r8), target :: ncai(state%psetcols,pver)    ! output number conc of ice nuclei available (1/m3)
   real(r8), target :: ncal(state%psetcols,pver)    ! output number conc of CCN (1/m3)
   real(r8), target :: qrout2(state%psetcols,pver)
   real(r8), target :: qsout2(state%psetcols,pver)
   real(r8), target :: nrout2(state%psetcols,pver)
   real(r8), target :: nsout2(state%psetcols,pver)
   real(r8), target :: freqs(state%psetcols,pver)
   real(r8), target :: freqr(state%psetcols,pver)
   real(r8), target :: nfice(state%psetcols,pver)
   real(r8), target :: qcrat(state%psetcols,pver)   ! qc limiter ratio (1=no limit)

   ! Object that packs columns with clouds/precip.
   type(MGPacker) :: packer

   ! Packed versions of inputs.
   real(r8), allocatable :: packed_t(:,:)
   real(r8), allocatable :: packed_q(:,:)
   real(r8), allocatable :: packed_qc(:,:)
   real(r8), allocatable :: packed_nc(:,:)
   real(r8), allocatable :: packed_qi(:,:)
   real(r8), allocatable :: packed_ni(:,:)
   real(r8), allocatable :: packed_qr(:,:)
   real(r8), allocatable :: packed_nr(:,:)
   real(r8), allocatable :: packed_qs(:,:)
   real(r8), allocatable :: packed_ns(:,:)

   real(r8), allocatable :: packed_relvar(:,:)
   real(r8), allocatable :: packed_accre_enhan(:,:)

   real(r8), allocatable :: packed_p(:,:)
   real(r8), allocatable :: packed_pdel(:,:)

   ! This is only needed for MG1.5, and can be removed when support for
   ! that version is dropped.
   real(r8), allocatable :: packed_pint(:,:)

   real(r8), allocatable :: packed_cldn(:,:)
   real(r8), allocatable :: packed_liqcldf(:,:)
   real(r8), allocatable :: packed_icecldf(:,:)

   real(r8), allocatable :: packed_naai(:,:)
   real(r8), allocatable :: packed_npccn(:,:)

   real(r8), allocatable :: packed_rndst(:,:,:)
   real(r8), allocatable :: packed_nacon(:,:,:)

   ! Optional outputs.
   real(r8), pointer :: packed_tnd_qsnow(:,:)
   real(r8), pointer :: packed_tnd_nsnow(:,:)
   real(r8), pointer :: packed_re_ice(:,:)

   real(r8), pointer :: packed_frzimm(:,:)
   real(r8), pointer :: packed_frzcnt(:,:)
   real(r8), pointer :: packed_frzdep(:,:)

   ! Output field post-processing.
   type(MGPostProc) :: post_proc

   ! Packed versions of outputs.
   real(r8), allocatable, target :: packed_rate1ord_cw2pr_st(:,:)
   real(r8), allocatable, target :: packed_tlat(:,:)
   real(r8), allocatable, target :: packed_qvlat(:,:)
   real(r8), allocatable, target :: packed_qctend(:,:)
   real(r8), allocatable, target :: packed_qitend(:,:)
   real(r8), allocatable, target :: packed_nctend(:,:)
   real(r8), allocatable, target :: packed_nitend(:,:)

   real(r8), allocatable, target :: packed_qrtend(:,:)
   real(r8), allocatable, target :: packed_qstend(:,:)
   real(r8), allocatable, target :: packed_nrtend(:,:)
   real(r8), allocatable, target :: packed_nstend(:,:)

   real(r8), allocatable, target :: packed_prect(:)
   real(r8), allocatable, target :: packed_preci(:)
   real(r8), allocatable, target :: packed_nevapr(:,:)
   real(r8), allocatable, target :: packed_am_evp_st(:,:)
   real(r8), allocatable, target :: packed_evapsnow(:,:)
   real(r8), allocatable, target :: packed_prain(:,:)
   real(r8), allocatable, target :: packed_prodsnow(:,:)
   real(r8), allocatable, target :: packed_cmeout(:,:)
   real(r8), allocatable, target :: packed_qsout(:,:)
   real(r8), allocatable, target :: packed_rflx(:,:)
   real(r8), allocatable, target :: packed_sflx(:,:)
   real(r8), allocatable, target :: packed_qrout(:,:)
   real(r8), allocatable, target :: packed_qcsevap(:,:)
   real(r8), allocatable, target :: packed_qisevap(:,:)
   real(r8), allocatable, target :: packed_qvres(:,:)
   real(r8), allocatable, target :: packed_cmei(:,:)
   real(r8), allocatable, target :: packed_vtrmc(:,:)
   real(r8), allocatable, target :: packed_vtrmi(:,:)
   real(r8), allocatable, target :: packed_qcsedten(:,:)
   real(r8), allocatable, target :: packed_qisedten(:,:)
   real(r8), allocatable, target :: packed_qrsedten(:,:)
   real(r8), allocatable, target :: packed_qssedten(:,:)
   real(r8), allocatable, target :: packed_umr(:,:)
   real(r8), allocatable, target :: packed_ums(:,:)
   real(r8), allocatable, target :: packed_pra(:,:)
   real(r8), allocatable, target :: packed_prc(:,:)
   real(r8), allocatable, target :: packed_mnuccc(:,:)
   real(r8), allocatable, target :: packed_mnucct(:,:)
   real(r8), allocatable, target :: packed_msacwi(:,:)
   real(r8), allocatable, target :: packed_psacws(:,:)
   real(r8), allocatable, target :: packed_bergs(:,:)
   real(r8), allocatable, target :: packed_berg(:,:)
   real(r8), allocatable, target :: packed_melt(:,:)
   real(r8), allocatable, target :: packed_homo(:,:)
   real(r8), allocatable, target :: packed_qcres(:,:)
   real(r8), allocatable, target :: packed_prci(:,:)
   real(r8), allocatable, target :: packed_prai(:,:)
   real(r8), allocatable, target :: packed_qires(:,:)
   real(r8), allocatable, target :: packed_mnuccr(:,:)
   real(r8), allocatable, target :: packed_pracs(:,:)
   real(r8), allocatable, target :: packed_meltsdt(:,:)
   real(r8), allocatable, target :: packed_frzrdt(:,:)
   real(r8), allocatable, target :: packed_mnuccd(:,:)
   real(r8), allocatable, target :: packed_nrout(:,:)
   real(r8), allocatable, target :: packed_nsout(:,:)
   real(r8), allocatable, target :: packed_refl(:,:)
   real(r8), allocatable, target :: packed_arefl(:,:)
   real(r8), allocatable, target :: packed_areflz(:,:)
   real(r8), allocatable, target :: packed_frefl(:,:)
   real(r8), allocatable, target :: packed_csrfl(:,:)
   real(r8), allocatable, target :: packed_acsrfl(:,:)
   real(r8), allocatable, target :: packed_fcsrfl(:,:)
   real(r8), allocatable, target :: packed_rercld(:,:)
   real(r8), allocatable, target :: packed_ncai(:,:)
   real(r8), allocatable, target :: packed_ncal(:,:)
   real(r8), allocatable, target :: packed_qrout2(:,:)
   real(r8), allocatable, target :: packed_qsout2(:,:)
   real(r8), allocatable, target :: packed_nrout2(:,:)
   real(r8), allocatable, target :: packed_nsout2(:,:)
   real(r8), allocatable, target :: packed_freqs(:,:)
   real(r8), allocatable, target :: packed_freqr(:,:)
   real(r8), allocatable, target :: packed_nfice(:,:)
   real(r8), allocatable, target :: packed_prer_evap(:,:)
   real(r8), allocatable, target :: packed_qcrat(:,:)

   real(r8), allocatable, target :: packed_rel(:,:)
   real(r8), allocatable, target :: packed_rei(:,:)
   real(r8), allocatable, target :: packed_lambdac(:,:)
   real(r8), allocatable, target :: packed_mu(:,:)
   real(r8), allocatable, target :: packed_des(:,:)
   real(r8), allocatable, target :: packed_dei(:,:)

   ! Dummy arrays for cases where we throw away the MG version and
   ! recalculate sizes on the CAM grid to avoid time/subcolumn averaging
   ! issues.
   real(r8), allocatable :: rel_fn_dum(:,:)
   real(r8), allocatable :: dsout2_dum(:,:)
   real(r8), allocatable :: drout_dum(:,:)
   real(r8), allocatable :: reff_rain_dum(:,:)
   real(r8), allocatable :: reff_snow_dum(:,:)

   ! Heterogeneous-only version of mnuccdo.
   real(r8) :: mnuccdohet(state%psetcols,pver)

   ! physics buffer fields for COSP simulator
   real(r8), pointer :: mgflxprc(:,:)     ! MG grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
   real(r8), pointer :: mgflxsnw(:,:)     ! MG grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
   real(r8), pointer :: mgmrprc(:,:)      ! MG grid-box mean mixingratio_large_scale_cloud_rain+snow at interfaces (kg/kg)
   real(r8), pointer :: mgmrsnw(:,:)      ! MG grid-box mean mixingratio_large_scale_cloud_snow at interfaces (kg/kg)
   real(r8), pointer :: mgreffrain_grid(:,:)   ! MG diagnostic rain effective radius (um)
   real(r8), pointer :: mgreffsnow_grid(:,:)   ! MG diagnostic snow effective radius (um)
   real(r8), pointer :: cvreffliq(:,:)    ! convective cloud liquid effective radius (um)
   real(r8), pointer :: cvreffice(:,:)    ! convective cloud ice effective radius (um)

   ! physics buffer fields used with CARMA
   real(r8), pointer, dimension(:,:) :: tnd_qsnow    ! external tendency on snow mass (kg/kg/s)
   real(r8), pointer, dimension(:,:) :: tnd_nsnow    ! external tendency on snow number(#/kg/s)
   real(r8), pointer, dimension(:,:) :: re_ice       ! ice effective radius (m)

   real(r8), pointer :: rate1ord_cw2pr_st(:,:) ! 1st order rate for direct conversion of
                                               ! strat. cloud water to precip (1/s)    ! rce 2010/05/01
   real(r8), pointer :: wsedl(:,:)        ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]


   real(r8), pointer :: CC_T(:,:)         ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_qv(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_ql(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_qi(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_nl(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_ni(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_qlst(:,:)      ! In-liquid stratus microphysical tendency

  ! variables for heterogeneous freezing
  real(r8), pointer :: frzimm(:,:)
  real(r8), pointer :: frzcnt(:,:)
  real(r8), pointer :: frzdep(:,:)

   real(r8), pointer :: qme(:,:)

   ! A local copy of state is used for diagnostic calculations
   type(physics_state) :: state_loc
   type(physics_ptend) :: ptend_loc

   real(r8) :: icecldf(state%psetcols,pver) ! Ice cloud fraction
   real(r8) :: liqcldf(state%psetcols,pver) ! Liquid cloud fraction (combined into cloud)

   real(r8), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
   real(r8), pointer :: rei(:,:)          ! Ice effective drop size (microns)

   real(r8), pointer :: cmeliq(:,:)

   real(r8), pointer :: cld(:,:)          ! Total cloud fraction
   real(r8), pointer :: concld(:,:)       ! Convective cloud fraction
   real(r8), pointer :: iciwpst(:,:)      ! Stratiform in-cloud ice water path for radiation
   real(r8), pointer :: iclwpst(:,:)      ! Stratiform in-cloud liquid water path for radiation
   real(r8), pointer :: cldfsnow(:,:)     ! Cloud fraction for liquid+snow
   real(r8), pointer :: icswp(:,:)        ! In-cloud snow water path

   real(r8) :: icimrst(state%psetcols,pver) ! In stratus ice mixing ratio
   real(r8) :: icwmrst(state%psetcols,pver) ! In stratus water mixing ratio
   real(r8) :: icinc(state%psetcols,pver)   ! In cloud ice number conc
   real(r8) :: icwnc(state%psetcols,pver)   ! In cloud water number conc

   real(r8) :: iclwpi(state%psetcols)       ! Vertically-integrated in-cloud Liquid WP before microphysics
   real(r8) :: iciwpi(state%psetcols)       ! Vertically-integrated in-cloud Ice WP before microphysics

   ! Averaging arrays for effective radius and number....
   real(r8) :: efiout_grid(pcols,pver)
   real(r8) :: efcout_grid(pcols,pver)
   real(r8) :: ncout_grid(pcols,pver)
   real(r8) :: niout_grid(pcols,pver)
   real(r8) :: freqi_grid(pcols,pver)
   real(r8) :: freql_grid(pcols,pver)

   real(r8) :: cdnumc_grid(pcols)           ! Vertically-integrated droplet concentration
   real(r8) :: icimrst_grid_out(pcols,pver) ! In stratus ice mixing ratio
   real(r8) :: icwmrst_grid_out(pcols,pver) ! In stratus water mixing ratio

   ! Cloud fraction used for precipitation.
   real(r8) :: cldmax_grid(pcols,pver)

   ! Average cloud top radius & number
   real(r8) :: ctrel_grid(pcols)
   real(r8) :: ctrei_grid(pcols)
   real(r8) :: ctnl_grid(pcols)
   real(r8) :: ctni_grid(pcols)
   real(r8) :: fcti_grid(pcols)
   real(r8) :: fctl_grid(pcols)

   real(r8) :: budget_ftem_grid(pcols,pver,6)

   ! Variables for precip efficiency calculation
   real(r8) :: minlwp        ! LWP threshold

   real(r8), pointer, dimension(:) :: acprecl_grid ! accumulated precip across timesteps
   real(r8), pointer, dimension(:) :: acgcme_grid  ! accumulated condensation across timesteps
   integer,  pointer, dimension(:) :: acnum_grid   ! counter for # timesteps accumulated

   ! Variables for liquid water path and column condensation
   real(r8) :: tgliqwp_grid(pcols)   ! column liquid
   real(r8) :: tgcmeliq_grid(pcols)  ! column condensation rate (units)

   real(r8) :: pe_grid(pcols)        ! precip efficiency for output
   real(r8) :: pefrac_grid(pcols)    ! fraction of time precip efficiency is written out
   real(r8) :: tpr_grid(pcols)       ! average accumulated precipitation rate in pe calculation

   ! variables for autoconversion and accretion vertical averages
   real(r8) :: vprco_grid(pcols)     ! vertical average autoconversion
   real(r8) :: vprao_grid(pcols)     ! vertical average accretion
   real(r8) :: racau_grid(pcols)     ! ratio of vertical averages
   integer  :: cnt_grid(pcols)       ! counters

   logical  :: lq(pcnst)

   real(r8) :: icimrst_grid(pcols,pver) ! stratus ice mixing ratio - on grid
   real(r8) :: icwmrst_grid(pcols,pver) ! stratus water mixing ratio - on grid

   real(r8), pointer :: lambdac_grid(:,:)
   real(r8), pointer :: mu_grid(:,:)
   real(r8), pointer :: rel_grid(:,:)
   real(r8), pointer :: rei_grid(:,:)
   real(r8), pointer :: dei_grid(:,:)
   real(r8), pointer :: des_grid(:,:)
   real(r8), pointer :: iclwpst_grid(:,:)

   real(r8) :: rho_grid(pcols,pver)
   real(r8) :: liqcldf_grid(pcols,pver)
   real(r8) :: qsout_grid(pcols,pver)
   real(r8) :: ncic_grid(pcols,pver)
   real(r8) :: niic_grid(pcols,pver)
   real(r8) :: rel_fn_grid(pcols,pver)    ! Ice effective drop size at fixed number (indirect effect) (microns) - on grid
   real(r8) :: qrout_grid(pcols,pver)
   real(r8) :: drout2_grid(pcols,pver)
   real(r8) :: dsout2_grid(pcols,pver)
   real(r8) :: nsout_grid(pcols,pver)
   real(r8) :: nrout_grid(pcols,pver)
   real(r8) :: reff_rain_grid(pcols,pver)
   real(r8) :: reff_snow_grid(pcols,pver)
   real(r8) :: cld_grid(pcols,pver)
   real(r8) :: pdel_grid(pcols,pver)
   real(r8) :: prco_grid(pcols,pver)
   real(r8) :: prao_grid(pcols,pver)
   real(r8) :: icecldf_grid(pcols,pver)
   real(r8) :: icwnc_grid(pcols,pver)
   real(r8) :: icinc_grid(pcols,pver)
   real(r8) :: qcreso_grid(pcols,pver)
   real(r8) :: melto_grid(pcols,pver)
   real(r8) :: mnuccco_grid(pcols,pver)
   real(r8) :: mnuccto_grid(pcols,pver)
   real(r8) :: bergo_grid(pcols,pver)
   real(r8) :: homoo_grid(pcols,pver)
   real(r8) :: msacwio_grid(pcols,pver)
   real(r8) :: psacwso_grid(pcols,pver)
   real(r8) :: bergso_grid(pcols,pver)
   real(r8) :: cmeiout_grid(pcols,pver)
   real(r8) :: qireso_grid(pcols,pver)
   real(r8) :: prcio_grid(pcols,pver)
   real(r8) :: praio_grid(pcols,pver)

   real(r8) :: nc_grid(pcols,pver)
   real(r8) :: ni_grid(pcols,pver)
   real(r8) :: qr_grid(pcols,pver)
   real(r8) :: nr_grid(pcols,pver)
   real(r8) :: qs_grid(pcols,pver)
   real(r8) :: ns_grid(pcols,pver)

   real(r8), pointer :: cmeliq_grid(:,:)

   real(r8), pointer :: prec_str_grid(:)
   real(r8), pointer :: snow_str_grid(:)
   real(r8), pointer :: prec_pcw_grid(:)
   real(r8), pointer :: snow_pcw_grid(:)
   real(r8), pointer :: prec_sed_grid(:)
   real(r8), pointer :: snow_sed_grid(:)
   real(r8), pointer :: cldo_grid(:,:)
   real(r8), pointer :: nevapr_grid(:,:)
   real(r8), pointer :: prain_grid(:,:)
   real(r8), pointer :: mgflxprc_grid(:,:)
   real(r8), pointer :: mgflxsnw_grid(:,:)
   real(r8), pointer :: mgmrprc_grid(:,:)
   real(r8), pointer :: mgmrsnw_grid(:,:)
   real(r8), pointer :: cvreffliq_grid(:,:)
   real(r8), pointer :: cvreffice_grid(:,:)
   real(r8), pointer :: rate1ord_cw2pr_st_grid(:,:)
   real(r8), pointer :: wsedl_grid(:,:)
   real(r8), pointer :: CC_t_grid(:,:)
   real(r8), pointer :: CC_qv_grid(:,:)
   real(r8), pointer :: CC_ql_grid(:,:)
   real(r8), pointer :: CC_qi_grid(:,:)
   real(r8), pointer :: CC_nl_grid(:,:)
   real(r8), pointer :: CC_ni_grid(:,:)
   real(r8), pointer :: CC_qlst_grid(:,:)
   real(r8), pointer :: qme_grid(:,:)
   real(r8), pointer :: iciwpst_grid(:,:)
   real(r8), pointer :: icswp_grid(:,:)
   real(r8), pointer :: ast_grid(:,:)
   real(r8), pointer :: cldfsnow_grid(:,:)

   real(r8), pointer :: qrout_grid_ptr(:,:)
   real(r8), pointer :: qsout_grid_ptr(:,:)
   real(r8), pointer :: nrout_grid_ptr(:,:)
   real(r8), pointer :: nsout_grid_ptr(:,:)

   integer :: nlev   ! number of levels where cloud physics is done
   integer :: mgncol ! size of mgcols
   integer, allocatable :: mgcols(:) ! Columns with microphysics performed

   logical :: use_subcol_microp
   integer :: col_type ! Flag to store whether accessing grid or sub-columns in pbuf_get_field

   character(128) :: errstring   ! return status (non-blank for error return)
   real(r8) :: rate_local

   ! For rrtmg optics. specified distribution.
   real(r8), parameter :: dcon   = 25.e-6_r8         ! Convective size distribution effective radius (meters)
   real(r8), parameter :: mucon  = 5.3_r8            ! Convective size distribution shape parameter
   real(r8), parameter :: deicon = 50._r8            ! Convective ice effective diameter (meters)

   real(r8), pointer :: pckdptr(:,:)

   ! Local variables for water tracers/isotopes

   real(r8), target :: preo(state%psetcols,pver)                              ! rain re-evaporation (kg/kg/sec)
   real(r8), target :: prdso(state%psetcols,pver)                             ! snow sublimation (kg/kg/sec)
   real(r8), target :: frzro(state%psetcols,pver)                             ! rain freezing (kg/kg/sec)
   real(r8), target :: meltso(state%psetcols,pver)                            ! snow melting  (kg/kg/sec)
   real(r8), target :: wtfc(state%psetcols,pver)                              ! Initial cloud liquid fall velocity
   real(r8), target :: wtfi(state%psetcols,pver)                              ! Initial cloud ice fall velocity
   real(r8), target :: wtprelat(state%psetcols,pver)                          ! Latent heat change due to pre_rates
   real(r8), target :: wtpostlat(state%psetcols,pver)                         ! Latent heat change due to post_rates 

   ! Water tracers/isotopes on the grid level
   real(r8) :: pre_rates_grid(pcols,pver,pwtype,pwtype,pwtype)    ! Process rates (kg/kg/sec)
   real(r8) :: sed_rates_grid(pcols,pver,pwtype)                  ! Sedimentation rates (kg/kg/sec)
   real(r8) :: post_rates_grid(pcols,pver,pwtype,pwtype,pwtype)   ! Process rates (kg/kg/sec)
   real(r8) :: pcmei_grid(pcols,pver)                             ! Positive cmeiout - deposition
   real(r8) :: ncmei_grid(pcols,pver)                             ! Negative cmeiout - sublimation
   real(r8) :: pmelts_grid(pcols,pver)                            ! Positive melts - melting
   real(r8) :: nmelts_grid(pcols,pver)                            ! Negative melts (freezing?)
   logical  :: isOk                                               ! Flag indicating test success

   ! Local packed arrays for water tracers/isotopes

   real(r8), allocatable, target :: packed_preo(:,:)              ! rain re-evaporation (kg/kg/sec)
   real(r8), allocatable, target :: packed_prdso(:,:)             ! snow sublimation (kg/kg/sec)
   real(r8), allocatable, target :: packed_frzro(:,:)             ! rain freezing (kg/kg/sec)
   real(r8), allocatable, target :: packed_meltso(:,:)            ! snow melting  (kg/kg/sec)
   real(r8), allocatable, target :: packed_wtfc(:,:)              ! Initial cloud liquid fall velocity
   real(r8), allocatable, target :: packed_wtfi(:,:)              ! Initial cloud ice fall velocity
   real(r8), allocatable, target :: packed_wtprelat(:,:)          ! Latent heat change due to pre_rates
   real(r8), allocatable, target :: packed_wtpostlat(:,:)         ! Latent heat change due to post_rates 

   ! above water tracers/isotopes arrays on the grid level

   real(r8), pointer :: preo_grid(:,:)                            ! rain re-evaporation (kg/kg/sec)
   real(r8), pointer :: prdso_grid(:,:)                           ! snow sublimation (kg/kg/sec)
   real(r8), pointer :: frzro_grid(:,:)                           ! rain freezing (kg/kg/sec)
   real(r8), pointer :: meltso_grid(:,:)                          ! snow melting  (kg/kg/sec)
   real(r8), pointer :: wtfc_grid(:,:)                            ! Initial cloud liquid fall velocity
   real(r8), pointer :: wtfi_grid(:,:)                            ! Initial cloud ice fall velocity
   real(r8), pointer :: wtprelat_grid(:,:)                        ! Latent heat change due to pre_rates
   real(r8), pointer :: wtpostlat_grid(:,:)                       ! Latent heat change due to post_rates 

   real(r8), pointer :: mnuccro_grid(:,:)
   real(r8), pointer :: pracso_grid(:,:)
   real(r8), pointer :: qcsedten_grid(:,:)
   real(r8), pointer :: qisedten_grid(:,:)
   real(r8), pointer :: alst_mic_grid(:,:)
   real(r8), pointer :: aist_mic_grid(:,:)
   character(len=*), parameter :: subname = 'micro_mg_cam_tend'
   !-------------------------------------------------------------------------------

   ! Find the number of levels used in the microphysics.
   nlev  = pver - top_lev + 1

   lchnk = state%lchnk
   ncol  = state%ncol
   psetcols = state%psetcols
   ngrdcol  = state%ngrdcol

   itim_old = pbuf_old_tim_idx()

   call phys_getopts(use_subcol_microp_out=use_subcol_microp)

   ! Set the col_type flag to grid or subcolumn dependent on the value of use_subcol_microp
   call pbuf_col_type_index(use_subcol_microp, col_type=col_type)

   !-----------------------
   ! These physics buffer fields are read only and not set in this parameterization
   ! If these fields do not have subcolumn data, copy the grid to the subcolumn if subcolumns is turned on
   ! If subcolumns is not turned on, then these fields will be grid data

   call pbuf_get_field(pbuf, naai_idx,        naai,        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, naai_hom_idx,    naai_hom,    col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, npccn_idx,       npccn,       col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, rndst_idx,       rndst,       col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, nacon_idx,       nacon,       col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, relvar_idx,      relvar,      col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan, col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq,      col_type=col_type, copy_if_needed=use_subcol_microp)

   call pbuf_get_field(pbuf, cld_idx,         cld,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, concld_idx,      concld,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, ast_idx,         ast,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)

   if (.not. do_cldice) then
      call pbuf_get_field(pbuf, tnd_qsnow_idx,   tnd_qsnow,   col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, tnd_nsnow_idx,   tnd_nsnow,   col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, re_ice_idx,      re_ice,      col_type=col_type, copy_if_needed=use_subcol_microp)
   end if

   if (use_hetfrz_classnuc) then
      call pbuf_get_field(pbuf, frzimm_idx, frzimm, col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, frzcnt_idx, frzcnt, col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, frzdep_idx, frzdep, col_type=col_type, copy_if_needed=use_subcol_microp)
   end if

   !-----------------------
   ! These physics buffer fields are calculated and set in this parameterization
   ! If subcolumns is turned on, then these fields will be calculated on a subcolumn grid, otherwise they will be a normal grid

   call pbuf_get_field(pbuf, prec_str_idx,    prec_str,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_str_idx,    snow_str,    col_type=col_type)
   call pbuf_get_field(pbuf, prec_pcw_idx,    prec_pcw,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_pcw_idx,    snow_pcw,    col_type=col_type)
   call pbuf_get_field(pbuf, prec_sed_idx,    prec_sed,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_sed_idx,    snow_sed,    col_type=col_type)
   call pbuf_get_field(pbuf, nevapr_idx,      nevapr,      col_type=col_type)
   call pbuf_get_field(pbuf, prer_evap_idx,   prer_evap,   col_type=col_type)
   call pbuf_get_field(pbuf, prain_idx,       prain,       col_type=col_type)
   call pbuf_get_field(pbuf, dei_idx,         dei,         col_type=col_type)
   call pbuf_get_field(pbuf, mu_idx,          mu,          col_type=col_type)
   call pbuf_get_field(pbuf, lambdac_idx,     lambdac,     col_type=col_type)
   call pbuf_get_field(pbuf, des_idx,         des,         col_type=col_type)
   call pbuf_get_field(pbuf, ls_flxprc_idx,   mgflxprc,    col_type=col_type)
   call pbuf_get_field(pbuf, ls_flxsnw_idx,   mgflxsnw,    col_type=col_type)
   call pbuf_get_field(pbuf, ls_mrprc_idx,    mgmrprc,     col_type=col_type)
   call pbuf_get_field(pbuf, ls_mrsnw_idx,    mgmrsnw,     col_type=col_type)
   call pbuf_get_field(pbuf, cv_reffliq_idx,  cvreffliq,   col_type=col_type)
   call pbuf_get_field(pbuf, cv_reffice_idx,  cvreffice,   col_type=col_type)
   call pbuf_get_field(pbuf, iciwpst_idx,     iciwpst,     col_type=col_type)
   call pbuf_get_field(pbuf, iclwpst_idx,     iclwpst,     col_type=col_type)
   call pbuf_get_field(pbuf, icswp_idx,       icswp,       col_type=col_type)
   call pbuf_get_field(pbuf, rel_idx,         rel,         col_type=col_type)
   call pbuf_get_field(pbuf, rei_idx,         rei,         col_type=col_type)
   call pbuf_get_field(pbuf, wsedl_idx,       wsedl,       col_type=col_type)
   call pbuf_get_field(pbuf, qme_idx,         qme,         col_type=col_type)

   call pbuf_get_field(pbuf, cldo_idx,        cldo,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cldfsnow_idx,    cldfsnow, start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_t_idx,        CC_t,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qv_idx,       CC_qv,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_ql_idx,       CC_ql,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qi_idx,       CC_qi,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_nl_idx,       CC_nl,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_ni_idx,       CC_ni,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qlst_idx,     CC_qlst,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)

   if (rate1_cw2pr_st_idx > 0) then
      call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st, col_type=col_type)
   end if

   if (qrain_idx > 0) call pbuf_get_field(pbuf, qrain_idx, qrout_grid_ptr)
   if (qsnow_idx > 0) call pbuf_get_field(pbuf, qsnow_idx, qsout_grid_ptr)
   if (nrain_idx > 0) call pbuf_get_field(pbuf, nrain_idx, nrout_grid_ptr)
   if (nsnow_idx > 0) call pbuf_get_field(pbuf, nsnow_idx, nsout_grid_ptr)

   !-----------------------
   ! If subcolumns is turned on, all calculated fields which are on subcolumns
   ! need to be retrieved on the grid as well for storing averaged values

   if (use_subcol_microp) then
      call pbuf_get_field(pbuf, prec_str_idx,    prec_str_grid)
      call pbuf_get_field(pbuf, snow_str_idx,    snow_str_grid)
      call pbuf_get_field(pbuf, prec_pcw_idx,    prec_pcw_grid)
      call pbuf_get_field(pbuf, snow_pcw_idx,    snow_pcw_grid)
      call pbuf_get_field(pbuf, prec_sed_idx,    prec_sed_grid)
      call pbuf_get_field(pbuf, snow_sed_idx,    snow_sed_grid)
      call pbuf_get_field(pbuf, nevapr_idx,      nevapr_grid)
      call pbuf_get_field(pbuf, prain_idx,       prain_grid)
      call pbuf_get_field(pbuf, dei_idx,         dei_grid)
      call pbuf_get_field(pbuf, mu_idx,          mu_grid)
      call pbuf_get_field(pbuf, lambdac_idx,     lambdac_grid)
      call pbuf_get_field(pbuf, des_idx,         des_grid)
      call pbuf_get_field(pbuf, ls_flxprc_idx,   mgflxprc_grid)
      call pbuf_get_field(pbuf, ls_flxsnw_idx,   mgflxsnw_grid)
      call pbuf_get_field(pbuf, ls_mrprc_idx,    mgmrprc_grid)
      call pbuf_get_field(pbuf, ls_mrsnw_idx,    mgmrsnw_grid)
      call pbuf_get_field(pbuf, cv_reffliq_idx,  cvreffliq_grid)
      call pbuf_get_field(pbuf, cv_reffice_idx,  cvreffice_grid)
      call pbuf_get_field(pbuf, iciwpst_idx,     iciwpst_grid)
      call pbuf_get_field(pbuf, iclwpst_idx,     iclwpst_grid)
      call pbuf_get_field(pbuf, icswp_idx,       icswp_grid)
      call pbuf_get_field(pbuf, rel_idx,         rel_grid)
      call pbuf_get_field(pbuf, rei_idx,         rei_grid)
      call pbuf_get_field(pbuf, wsedl_idx,       wsedl_grid)
      call pbuf_get_field(pbuf, qme_idx,         qme_grid)

      call pbuf_get_field(pbuf, cldo_idx,     cldo_grid,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cldfsnow_idx, cldfsnow_grid, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_t_idx,     CC_t_grid,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_qv_idx,    CC_qv_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_ql_idx,    CC_ql_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_qi_idx,    CC_qi_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_nl_idx,    CC_nl_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_ni_idx,    CC_ni_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_qlst_idx,  CC_qlst_grid,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

      if (rate1_cw2pr_st_idx > 0) then
         call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st_grid)
      end if

   end if

   !-----------------------
   ! These are only on the grid regardless of whether subcolumns are turned on or not
   call pbuf_get_field(pbuf, ls_reffrain_idx, mgreffrain_grid)
   call pbuf_get_field(pbuf, ls_reffsnow_idx, mgreffsnow_grid)
   call pbuf_get_field(pbuf, acpr_idx,        acprecl_grid)
   call pbuf_get_field(pbuf, acgcme_idx,      acgcme_grid)
   call pbuf_get_field(pbuf, acnum_idx,       acnum_grid)
   call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq_grid)
   call pbuf_get_field(pbuf, ast_idx,         ast_grid, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, evprain_st_idx,  evprain_st_grid)
   call pbuf_get_field(pbuf, evpsnow_st_idx,  evpsnow_st_grid)

   ! Only MG 1 defines this field so far.
   if (micro_mg_version == 1 .and. micro_mg_sub_version == 0) then
      call pbuf_get_field(pbuf, am_evp_st_idx,   am_evp_st_grid)
   end if

   !-------------------------------------------------------------------------------------
   ! Microphysics assumes 'liquid stratus frac = ice stratus frac
   !                      = max( liquid stratus frac, ice stratus frac )'.
   alst_mic => ast
   aist_mic => ast

   ! Output initial in-cloud LWP (before microphysics)

   iclwpi = 0._r8
   iciwpi = 0._r8

   do i = 1, ncol
      do k = top_lev, pver
         iclwpi(i) = iclwpi(i) + &
              min(state%q(i,k,ixcldliq) / max(mincld,ast(i,k)),0.005_r8) &
              * state%pdel(i,k) / gravit
         iciwpi(i) = iciwpi(i) + &
              min(state%q(i,k,ixcldice) / max(mincld,ast(i,k)),0.005_r8) &
              * state%pdel(i,k) / gravit
      end do
   end do

   cldo(:ncol,top_lev:pver)=ast(:ncol,top_lev:pver)

   ! Initialize local state from input.
   call physics_state_copy(state, state_loc)

   ! Initialize ptend for output.
   lq = .false.
   lq(1) = .true.
   lq(ixcldliq) = .true.
   lq(ixcldice) = .true.
   lq(ixnumliq) = .true.
   lq(ixnumice) = .true.
   if (micro_mg_version > 1) then
      lq(ixrain) = .true.
      lq(ixsnow) = .true.
      lq(ixnumrain) = .true.
      lq(ixnumsnow) = .true.
   end if

   !Water tracers:
   if ( trace_water ) then
      do i=1,wtrc_ncnst
        lq(wtrc_indices(i)) = .true.
      end do
   end if

   ! the name 'cldwat' triggers special tests on cldliq
   ! and cldice in physics_update
   call physics_ptend_init(ptend, psetcols, "cldwat", ls=.true., lq=lq)

   select case (micro_mg_version)
   case (1)
      select case (micro_mg_sub_version)
      case (0)
         call micro_mg_get_cols1_0(ncol, nlev, top_lev, state%q(:,:,ixcldliq), &
              state%q(:,:,ixcldice), mgncol, mgcols)
      case (5)
         call micro_mg_get_cols1_5(ncol, nlev, top_lev, state%q(:,:,ixcldliq), &
              state%q(:,:,ixcldice), mgncol, mgcols)
      end select
   case (2)
      call micro_mg_get_cols2_0(ncol, nlev, top_lev, state%q(:,:,ixcldliq), &
           state%q(:,:,ixcldice), state%q(:,:,ixrain), state%q(:,:,ixsnow), &
           mgncol, mgcols)
   end select

   packer = MGPacker(psetcols, pver, mgcols, top_lev)
   post_proc = MGPostProc(packer)

   allocate(packed_rate1ord_cw2pr_st(mgncol,nlev))
   pckdptr => packed_rate1ord_cw2pr_st ! workaround an apparent pgi compiler bug on goldbach
   call post_proc%add_field(p(rate1cld), pckdptr)
   allocate(packed_tlat(mgncol,nlev))
   call post_proc%add_field(p(tlat), p(packed_tlat))
   allocate(packed_qvlat(mgncol,nlev))
   call post_proc%add_field(p(qvlat), p(packed_qvlat))
   allocate(packed_qctend(mgncol,nlev))
   call post_proc%add_field(p(qcten), p(packed_qctend))
   allocate(packed_qitend(mgncol,nlev))
   call post_proc%add_field(p(qiten), p(packed_qitend))
   allocate(packed_nctend(mgncol,nlev))
   call post_proc%add_field(p(ncten), p(packed_nctend))
   allocate(packed_nitend(mgncol,nlev))
   call post_proc%add_field(p(niten), p(packed_nitend))

   if (micro_mg_version > 1) then
      allocate(packed_qrtend(mgncol,nlev))
      call post_proc%add_field(p(qrten), p(packed_qrtend))
      allocate(packed_qstend(mgncol,nlev))
      call post_proc%add_field(p(qsten), p(packed_qstend))
      allocate(packed_nrtend(mgncol,nlev))
      call post_proc%add_field(p(nrten), p(packed_nrtend))
      allocate(packed_nstend(mgncol,nlev))
      call post_proc%add_field(p(nsten), p(packed_nstend))
      allocate(packed_umr(mgncol,nlev))
      call post_proc%add_field(p(umr), p(packed_umr))
      allocate(packed_ums(mgncol,nlev))
      call post_proc%add_field(p(ums), p(packed_ums))
   else if (micro_mg_sub_version == 0) then
      allocate(packed_am_evp_st(mgncol,nlev))
      call post_proc%add_field(p(am_evp_st), p(packed_am_evp_st))
   end if

   allocate(packed_prect(mgncol))
   call post_proc%add_field(p(prect), p(packed_prect))
   allocate(packed_preci(mgncol))
   call post_proc%add_field(p(preci), p(packed_preci))
   allocate(packed_nevapr(mgncol,nlev))
   call post_proc%add_field(p(nevapr), p(packed_nevapr))
   allocate(packed_evapsnow(mgncol,nlev))
   call post_proc%add_field(p(evapsnow), p(packed_evapsnow))
   allocate(packed_prain(mgncol,nlev))
   call post_proc%add_field(p(prain), p(packed_prain))
   allocate(packed_prodsnow(mgncol,nlev))
   call post_proc%add_field(p(prodsnow), p(packed_prodsnow))
   allocate(packed_cmeout(mgncol,nlev))
   call post_proc%add_field(p(cmeice), p(packed_cmeout))
   allocate(packed_qsout(mgncol,nlev))
   call post_proc%add_field(p(qsout), p(packed_qsout))
   allocate(packed_rflx(mgncol,nlev+1))
   call post_proc%add_field(p(rflx), p(packed_rflx))
   allocate(packed_sflx(mgncol,nlev+1))
   call post_proc%add_field(p(sflx), p(packed_sflx))
   allocate(packed_qrout(mgncol,nlev))
   call post_proc%add_field(p(qrout), p(packed_qrout))
   allocate(packed_qcsevap(mgncol,nlev))
   call post_proc%add_field(p(qcsevap), p(packed_qcsevap))
   allocate(packed_qisevap(mgncol,nlev))
   call post_proc%add_field(p(qisevap), p(packed_qisevap))
   allocate(packed_qvres(mgncol,nlev))
   call post_proc%add_field(p(qvres), p(packed_qvres))
   allocate(packed_cmei(mgncol,nlev))
   call post_proc%add_field(p(cmeiout), p(packed_cmei))
   allocate(packed_vtrmc(mgncol,nlev))
   call post_proc%add_field(p(vtrmc), p(packed_vtrmc))
   allocate(packed_vtrmi(mgncol,nlev))
   call post_proc%add_field(p(vtrmi), p(packed_vtrmi))
   allocate(packed_qcsedten(mgncol,nlev))
   call post_proc%add_field(p(qcsedten), p(packed_qcsedten))
   allocate(packed_qisedten(mgncol,nlev))
   call post_proc%add_field(p(qisedten), p(packed_qisedten))
   if (micro_mg_version > 1) then
      allocate(packed_qrsedten(mgncol,nlev))
      call post_proc%add_field(p(qrsedten), p(packed_qrsedten))
      allocate(packed_qssedten(mgncol,nlev))
      call post_proc%add_field(p(qssedten), p(packed_qssedten))
   end if

   allocate(packed_pra(mgncol,nlev))
   call post_proc%add_field(p(prao), p(packed_pra))
   allocate(packed_prc(mgncol,nlev))
   call post_proc%add_field(p(prco), p(packed_prc))
   allocate(packed_mnuccc(mgncol,nlev))
   call post_proc%add_field(p(mnuccco), p(packed_mnuccc))
   allocate(packed_mnucct(mgncol,nlev))
   call post_proc%add_field(p(mnuccto), p(packed_mnucct))
   allocate(packed_msacwi(mgncol,nlev))
   call post_proc%add_field(p(msacwio), p(packed_msacwi))
   allocate(packed_psacws(mgncol,nlev))
   call post_proc%add_field(p(psacwso), p(packed_psacws))
   allocate(packed_bergs(mgncol,nlev))
   call post_proc%add_field(p(bergso), p(packed_bergs))
   allocate(packed_berg(mgncol,nlev))
   call post_proc%add_field(p(bergo), p(packed_berg))
   allocate(packed_melt(mgncol,nlev))
   call post_proc%add_field(p(melto), p(packed_melt))
   allocate(packed_homo(mgncol,nlev))
   call post_proc%add_field(p(homoo), p(packed_homo))
   allocate(packed_qcres(mgncol,nlev))
   call post_proc%add_field(p(qcreso), p(packed_qcres))
   allocate(packed_prci(mgncol,nlev))
   call post_proc%add_field(p(prcio), p(packed_prci))
   allocate(packed_prai(mgncol,nlev))
   call post_proc%add_field(p(praio), p(packed_prai))
   allocate(packed_qires(mgncol,nlev))
   call post_proc%add_field(p(qireso), p(packed_qires))
   allocate(packed_mnuccr(mgncol,nlev))
   call post_proc%add_field(p(mnuccro), p(packed_mnuccr))
   allocate(packed_pracs(mgncol,nlev))
   call post_proc%add_field(p(pracso), p(packed_pracs))
   allocate(packed_meltsdt(mgncol,nlev))
   call post_proc%add_field(p(meltsdt), p(packed_meltsdt))
   allocate(packed_frzrdt(mgncol,nlev))
   call post_proc%add_field(p(frzrdt), p(packed_frzrdt))
   allocate(packed_mnuccd(mgncol,nlev))
   call post_proc%add_field(p(mnuccdo), p(packed_mnuccd))
   allocate(packed_nrout(mgncol,nlev))
   call post_proc%add_field(p(nrout), p(packed_nrout))
   allocate(packed_nsout(mgncol,nlev))
   call post_proc%add_field(p(nsout), p(packed_nsout))

   allocate(packed_refl(mgncol,nlev))
   call post_proc%add_field(p(refl), p(packed_refl), fillvalue=-9999._r8)
   allocate(packed_arefl(mgncol,nlev))
   call post_proc%add_field(p(arefl), p(packed_arefl))
   allocate(packed_areflz(mgncol,nlev))
   call post_proc%add_field(p(areflz), p(packed_areflz))
   allocate(packed_frefl(mgncol,nlev))
   call post_proc%add_field(p(frefl), p(packed_frefl))
   allocate(packed_csrfl(mgncol,nlev))
   call post_proc%add_field(p(csrfl), p(packed_csrfl), fillvalue=-9999._r8)
   allocate(packed_acsrfl(mgncol,nlev))
   call post_proc%add_field(p(acsrfl), p(packed_acsrfl))
   allocate(packed_fcsrfl(mgncol,nlev))
   call post_proc%add_field(p(fcsrfl), p(packed_fcsrfl))

   allocate(packed_rercld(mgncol,nlev))
   call post_proc%add_field(p(rercld), p(packed_rercld))
   allocate(packed_ncai(mgncol,nlev))
   call post_proc%add_field(p(ncai), p(packed_ncai))
   allocate(packed_ncal(mgncol,nlev))
   call post_proc%add_field(p(ncal), p(packed_ncal))
   allocate(packed_qrout2(mgncol,nlev))
   call post_proc%add_field(p(qrout2), p(packed_qrout2))
   allocate(packed_qsout2(mgncol,nlev))
   call post_proc%add_field(p(qsout2), p(packed_qsout2))
   allocate(packed_nrout2(mgncol,nlev))
   call post_proc%add_field(p(nrout2), p(packed_nrout2))
   allocate(packed_nsout2(mgncol,nlev))
   call post_proc%add_field(p(nsout2), p(packed_nsout2))
   allocate(packed_freqs(mgncol,nlev))
   call post_proc%add_field(p(freqs), p(packed_freqs))
   allocate(packed_freqr(mgncol,nlev))
   call post_proc%add_field(p(freqr), p(packed_freqr))
   allocate(packed_nfice(mgncol,nlev))
   call post_proc%add_field(p(nfice), p(packed_nfice))
   if (micro_mg_version /= 1 .or. micro_mg_sub_version /= 0) then
      allocate(packed_qcrat(mgncol,nlev))
      call post_proc%add_field(p(qcrat), p(packed_qcrat), fillvalue=1._r8)
   end if

   ! The following are all variables related to sizes, where it does not
   ! necessarily make sense to average over time steps. Instead, we keep
   ! the value from the last substep, which is what "accum_null" does.
   allocate(packed_rel(mgncol,nlev))
   call post_proc%add_field(p(rel), p(packed_rel), &
        fillvalue=10._r8, accum_method=accum_null)
   allocate(packed_rei(mgncol,nlev))
   call post_proc%add_field(p(rei), p(packed_rei), &
        fillvalue=25._r8, accum_method=accum_null)
   allocate(packed_lambdac(mgncol,nlev))
   call post_proc%add_field(p(lambdac), p(packed_lambdac), &
        accum_method=accum_null)
   allocate(packed_mu(mgncol,nlev))
   call post_proc%add_field(p(mu), p(packed_mu), &
        accum_method=accum_null)
   allocate(packed_des(mgncol,nlev))
   call post_proc%add_field(p(des), p(packed_des), &
        accum_method=accum_null)
   allocate(packed_dei(mgncol,nlev))
   call post_proc%add_field(p(dei), p(packed_dei), &
        accum_method=accum_null)
   allocate(packed_prer_evap(mgncol,nlev))
   call post_proc%add_field(p(prer_evap), p(packed_prer_evap), &
        accum_method=accum_null)

   if (micro_mg_version == 1 .and. micro_mg_sub_version == 0 ) then
      allocate(packed_preo(mgncol,nlev))
      call post_proc%add_field(p(preo), p(packed_preo))
      allocate(packed_prdso(mgncol,nlev))
      call post_proc%add_field(p(prdso), p(packed_prdso))
      allocate(packed_meltso(mgncol,nlev))
      call post_proc%add_field(p(meltso), p(packed_meltso))
      allocate(packed_wtfc(mgncol,nlev))
      call post_proc%add_field(p(wtfc), p(packed_wtfc))
      allocate(packed_wtfi(mgncol,nlev))
      call post_proc%add_field(p(wtfi), p(packed_wtfi))
      allocate(packed_wtprelat(mgncol,nlev))
      call post_proc%add_field(p(wtprelat), p(packed_wtprelat))
      allocate(packed_wtpostlat(mgncol,nlev))
      call post_proc%add_field(p(wtpostlat), p(packed_wtpostlat))
      allocate(packed_frzro(mgncol,nlev))
      call post_proc%add_field(p(frzro), p(packed_frzro))
   end if

   ! Allocate all the dummies with MG sizes.
   allocate(rel_fn_dum(mgncol,nlev))
   allocate(dsout2_dum(mgncol,nlev))
   allocate(drout_dum(mgncol,nlev))
   allocate(reff_rain_dum(mgncol,nlev))
   allocate(reff_snow_dum(mgncol,nlev))

   ! Pack input variables that are not updated during substeps.
   allocate(packed_relvar(mgncol,nlev))
   packed_relvar = packer%pack(relvar)
   allocate(packed_accre_enhan(mgncol,nlev))
   packed_accre_enhan = packer%pack(accre_enhan)

   allocate(packed_p(mgncol,nlev))
   packed_p = packer%pack(state_loc%pmid)
   allocate(packed_pdel(mgncol,nlev))
   packed_pdel = packer%pack(state_loc%pdel)

   allocate(packed_pint(mgncol,nlev+1))
   packed_pint = packer%pack_interface(state_loc%pint)

   allocate(packed_cldn(mgncol,nlev))
   packed_cldn = packer%pack(ast)
   allocate(packed_liqcldf(mgncol,nlev))
   packed_liqcldf = packer%pack(alst_mic)
   allocate(packed_icecldf(mgncol,nlev))
   packed_icecldf = packer%pack(aist_mic)

   allocate(packed_naai(mgncol,nlev))
   packed_naai = packer%pack(naai)
   allocate(packed_npccn(mgncol,nlev))
   packed_npccn = packer%pack(npccn)

   allocate(packed_rndst(mgncol,nlev,size(rndst, 3)))
   packed_rndst = packer%pack(rndst)
   allocate(packed_nacon(mgncol,nlev,size(nacon, 3)))
   packed_nacon = packer%pack(nacon)

   if (.not. do_cldice) then
      allocate(packed_tnd_qsnow(mgncol,nlev))
      packed_tnd_qsnow = packer%pack(tnd_qsnow)
      allocate(packed_tnd_nsnow(mgncol,nlev))
      packed_tnd_nsnow = packer%pack(tnd_nsnow)
      allocate(packed_re_ice(mgncol,nlev))
      packed_re_ice = packer%pack(re_ice)
   else
      nullify(packed_tnd_qsnow)
      nullify(packed_tnd_nsnow)
      nullify(packed_re_ice)
   end if

   if (use_hetfrz_classnuc) then
      allocate(packed_frzimm(mgncol,nlev))
      packed_frzimm = packer%pack(frzimm)
      allocate(packed_frzcnt(mgncol,nlev))
      packed_frzcnt = packer%pack(frzcnt)
      allocate(packed_frzdep(mgncol,nlev))
      packed_frzdep = packer%pack(frzdep)
   else
      nullify(packed_frzimm)
      nullify(packed_frzcnt)
      nullify(packed_frzdep)
   end if

   ! Allocate input variables that are updated during substeps.
   allocate(packed_t(mgncol,nlev))
   allocate(packed_q(mgncol,nlev))
   allocate(packed_qc(mgncol,nlev))
   allocate(packed_nc(mgncol,nlev))
   allocate(packed_qi(mgncol,nlev))
   allocate(packed_ni(mgncol,nlev))
   if (micro_mg_version > 1) then
      allocate(packed_qr(mgncol,nlev))
      allocate(packed_nr(mgncol,nlev))
      allocate(packed_qs(mgncol,nlev))
      allocate(packed_ns(mgncol,nlev))
   end if

   do it = 1, num_steps

      ! Pack input variables that are updated during substeps.
      packed_t = packer%pack(state_loc%t)
      packed_q = packer%pack(state_loc%q(:,:,1))
      packed_qc = packer%pack(state_loc%q(:,:,ixcldliq))
      packed_nc = packer%pack(state_loc%q(:,:,ixnumliq))
      packed_qi = packer%pack(state_loc%q(:,:,ixcldice))
      packed_ni = packer%pack(state_loc%q(:,:,ixnumice))
      if (micro_mg_version > 1) then
         packed_qr = packer%pack(state_loc%q(:,:,ixrain))
         packed_nr = packer%pack(state_loc%q(:,:,ixnumrain))
         packed_qs = packer%pack(state_loc%q(:,:,ixsnow))
         packed_ns = packer%pack(state_loc%q(:,:,ixnumsnow))
      end if

      select case (micro_mg_version)
      case (1)
         select case (micro_mg_sub_version)
         case (0)

            call micro_mg_tend1_0( &
                 microp_uniform, mgncol, nlev, mgncol, 1, dtime/num_steps, &
                 packed_t, packed_q, packed_qc, packed_qi, packed_nc,     &
                 packed_ni, packed_p, packed_pdel, packed_cldn, packed_liqcldf,&
                 packed_relvar, packed_accre_enhan,                             &
                 packed_icecldf, packed_rate1ord_cw2pr_st, packed_naai, packed_npccn,                 &
                 packed_rndst, packed_nacon, packed_tlat, packed_qvlat, packed_qctend,                &
                 packed_qitend, packed_nctend, packed_nitend, packed_rel, rel_fn_dum,      &
                 packed_rei, packed_prect, packed_preci, packed_nevapr, packed_evapsnow, packed_am_evp_st, &
                 packed_prain, packed_prodsnow, packed_cmeout, packed_dei, packed_mu,                &
                 packed_lambdac, packed_qsout, packed_des, packed_rflx, packed_sflx,                 &
                 packed_qrout, reff_rain_dum, reff_snow_dum, packed_qcsevap, packed_qisevap,   &
                 packed_qvres, packed_cmei, packed_vtrmc, packed_vtrmi, packed_qcsedten,          &
                 packed_qisedten, packed_pra, packed_prc, packed_mnuccc, packed_mnucct,          &
                 packed_msacwi, packed_psacws, packed_bergs, packed_berg, packed_melt,          &
                 packed_homo, packed_qcres, packed_prci, packed_prai, packed_qires,             &
                 packed_mnuccr, packed_pracs, packed_meltsdt, packed_frzrdt, packed_mnuccd,       &
                 packed_nrout, packed_nsout, packed_refl, packed_arefl, packed_areflz,               &
                 packed_frefl, packed_csrfl, packed_acsrfl, packed_fcsrfl, packed_rercld,            &
                 packed_ncai, packed_ncal, packed_qrout2, packed_qsout2, packed_nrout2,              &
                 packed_nsout2, drout_dum, dsout2_dum, packed_freqs,packed_freqr,            &
                 packed_nfice, packed_prer_evap, do_cldice, errstring, &
                 packed_tnd_qsnow, packed_tnd_nsnow, packed_re_ice, &
                 packed_frzimm, packed_frzcnt, packed_frzdep, packed_preo, packed_prdso,     &
                 packed_frzro, packed_meltso, packed_wtfc, packed_wtfi, packed_wtprelat,     &
                 packed_wtpostlat)

         case (5)

            call micro_mg_tend1_5( &
                 mgncol,   nlev,     dtime/num_steps,    &
                 packed_t,       packed_q,                     &
                 packed_qc,      packed_qi,    &
                 packed_nc,      packed_ni,    &
                 packed_relvar,             packed_accre_enhan,                            &
                 packed_p,     packed_pdel,     packed_pint,     &
                 packed_cldn,                packed_liqcldf,           packed_icecldf,           &
                 packed_rate1ord_cw2pr_st,           packed_naai,     packed_npccn,    packed_rndst,    packed_nacon,    &
                 packed_tlat,     packed_qvlat,    packed_qctend,    packed_qitend,    packed_nctend,    packed_nitend,    &
                 packed_rel,      rel_fn_dum,   packed_rei,                packed_prect,    packed_preci,    &
                 packed_nevapr,   packed_evapsnow, packed_prain,    packed_prodsnow, packed_cmeout,   packed_dei,      &
                 packed_mu,       packed_lambdac,  packed_qsout,    packed_des,      packed_rflx,     packed_sflx,     &
                 packed_qrout,              reff_rain_dum,          reff_snow_dum,          &
                 packed_qcsevap,  packed_qisevap,  packed_qvres,    packed_cmei,  packed_vtrmc,   packed_vtrmi,    &
                 packed_qcsedten, packed_qisedten, packed_pra,     packed_prc,     packed_mnuccc,  packed_mnucct,  &
                 packed_msacwi,  packed_psacws,  packed_bergs,   packed_berg,    packed_melt,    packed_homo,    &
                 packed_qcres,             packed_prci,    packed_prai,    packed_qires,             &
                 packed_mnuccr,  packed_pracs,   packed_meltsdt,  packed_frzrdt,   packed_mnuccd,            &
                 packed_nrout,   packed_nsout,    packed_refl,     packed_arefl,    packed_areflz,   packed_frefl,    &
                 packed_csrfl,    packed_acsrfl,   packed_fcsrfl,             packed_rercld,             &
                 packed_ncai,     packed_ncal,     packed_qrout2,   packed_qsout2,   packed_nrout2,   packed_nsout2,   &
                 drout_dum,   dsout2_dum,   packed_freqs,    packed_freqr,    packed_nfice,    packed_qcrat,    &
                 errstring, &
                 packed_tnd_qsnow,          packed_tnd_nsnow,          packed_re_ice, packed_prer_evap,             &
                 packed_frzimm, packed_frzcnt, packed_frzdep)

         end select
      case(2)
         select case (micro_mg_sub_version)
         case (0)

            call micro_mg_tend2_0( &
                 mgncol,         nlev,           dtime/num_steps,&
                 packed_t,               packed_q,               &
                 packed_qc,              packed_qi,              &
                 packed_nc,              packed_ni,              &
                 packed_qr,              packed_qs,              &
                 packed_nr,              packed_ns,              &
                 packed_relvar,          packed_accre_enhan,     &
                 packed_p,               packed_pdel,            &
                 packed_cldn,    packed_liqcldf, packed_icecldf, &
                 packed_rate1ord_cw2pr_st,                       &
                 packed_naai,            packed_npccn,           &
                 packed_rndst,           packed_nacon,           &
                 packed_tlat,            packed_qvlat,           &
                 packed_qctend,          packed_qitend,          &
                 packed_nctend,          packed_nitend,          &
                 packed_qrtend,          packed_qstend,          &
                 packed_nrtend,          packed_nstend,          &
                 packed_rel,     rel_fn_dum,     packed_rei,     &
                 packed_prect,           packed_preci,           &
                 packed_nevapr,          packed_evapsnow,        &
                 packed_prain,           packed_prodsnow,        &
                 packed_cmeout,          packed_dei,             &
                 packed_mu,              packed_lambdac,         &
                 packed_qsout,           packed_des,             &
                 packed_rflx,    packed_sflx,    packed_qrout,   &
                 reff_rain_dum,          reff_snow_dum,          &
                 packed_qcsevap, packed_qisevap, packed_qvres,   &
                 packed_cmei,    packed_vtrmc,   packed_vtrmi,   &
                 packed_umr,             packed_ums,             &
                 packed_qcsedten,        packed_qisedten,        &
                 packed_pra,             packed_prc,             &
                 packed_mnuccc,  packed_mnucct,  packed_msacwi,  &
                 packed_psacws,  packed_bergs,   packed_berg,    &
                 packed_melt,            packed_homo,            &
                 packed_qcres,   packed_prci,    packed_prai,    &
                 packed_qires,   packed_mnuccr,  packed_pracs,   &
                 packed_meltsdt, packed_frzrdt,  packed_mnuccd,  &
                 packed_nrout,           packed_nsout,           &
                 packed_refl,    packed_arefl,   packed_areflz,  &
                 packed_frefl,   packed_csrfl,   packed_acsrfl,  &
                 packed_fcsrfl,          packed_rercld,          &
                 packed_ncai,            packed_ncal,            &
                 packed_qrout2,          packed_qsout2,          &
                 packed_nrout2,          packed_nsout2,          &
                 drout_dum,              dsout2_dum,             &
                 packed_freqs,           packed_freqr,           &
                 packed_nfice,           packed_qcrat,           &
                 errstring, &
                 packed_tnd_qsnow,packed_tnd_nsnow,packed_re_ice,&
		 packed_prer_evap,                                     &
                 packed_frzimm,  packed_frzcnt,  packed_frzdep   )
         end select
      end select

      call handle_errmsg(errstring, subname="micro_mg_tend")

      call physics_ptend_init(ptend_loc, psetcols, "micro_mg", &
                              ls=.true., lq=lq)

      ! Set local tendency.
      ptend_loc%s               = packer%unpack(packed_tlat, 0._r8)
      ptend_loc%q(:,:,1)        = packer%unpack(packed_qvlat, 0._r8)
      ptend_loc%q(:,:,ixcldliq) = packer%unpack(packed_qctend, 0._r8)
      ptend_loc%q(:,:,ixcldice) = packer%unpack(packed_qitend, 0._r8)
      ptend_loc%q(:,:,ixnumliq) = packer%unpack(packed_nctend, &
           -state_loc%q(:,:,ixnumliq)/(dtime/num_steps))
      if (do_cldice) then
         ptend_loc%q(:,:,ixnumice) = packer%unpack(packed_nitend, &
              -state_loc%q(:,:,ixnumice)/(dtime/num_steps))
      else
         ! In this case, the tendency should be all 0.
         if (any(packed_nitend /= 0._r8)) &
              call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
              " but micro_mg_tend has ice number tendencies.")
         ptend_loc%q(:,:,ixnumice) = 0._r8
      end if

      if (micro_mg_version > 1) then
         ptend_loc%q(:,:,ixrain)    = packer%unpack(packed_qrtend, 0._r8)
         ptend_loc%q(:,:,ixsnow)    = packer%unpack(packed_qstend, 0._r8)
         ptend_loc%q(:,:,ixnumrain) = packer%unpack(packed_nrtend, &
              -state_loc%q(:,:,ixnumrain)/(dtime/num_steps))
         ptend_loc%q(:,:,ixnumsnow) = packer%unpack(packed_nstend, &
              -state_loc%q(:,:,ixnumsnow)/(dtime/num_steps))
      end if

      ! Sum into overall ptend
      call physics_ptend_sum(ptend_loc, ptend, ncol)

      ! Update local state
      call physics_update(state_loc, ptend_loc, dtime/num_steps)

      ! Sum all outputs for averaging.
      call post_proc%accumulate()

   end do

   ! Divide ptend by substeps.
   call physics_ptend_scale(ptend, 1._r8/num_steps, ncol)

   ! Use summed outputs to produce averages
   call post_proc%process_and_unpack()

   call post_proc%finalize()

   if (associated(packed_tnd_qsnow)) deallocate(packed_tnd_qsnow)
   if (associated(packed_tnd_nsnow)) deallocate(packed_tnd_nsnow)
   if (associated(packed_re_ice)) deallocate(packed_re_ice)
   if (associated(packed_frzimm)) deallocate(packed_frzimm)
   if (associated(packed_frzcnt)) deallocate(packed_frzcnt)
   if (associated(packed_frzdep)) deallocate(packed_frzdep)

   ! Check to make sure that the microphysics code is respecting the flags that control
   ! whether MG should be prognosing cloud ice and cloud liquid or not.
   if (.not. do_cldice) then
      if (any(ptend%q(:ncol,top_lev:pver,ixcldice) /= 0.0_r8)) &
           call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
           " but micro_mg_tend has ice mass tendencies.")
      if (any(ptend%q(:ncol,top_lev:pver,ixnumice) /= 0.0_r8)) &
           call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
           " but micro_mg_tend has ice number tendencies.")
   end if
   if (.not. do_cldliq) then
      if (any(ptend%q(:ncol,top_lev:pver,ixcldliq) /= 0.0_r8)) &
           call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
           " but micro_mg_tend has liquid mass tendencies.")
      if (any(ptend%q(:ncol,top_lev:pver,ixnumliq) /= 0.0_r8)) &
           call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
           " but micro_mg_tend has liquid number tendencies.")
   end if

   call micro_mg_cam_postmg_diag(ncol, psetcols, micro_mg_version, rate1_cw2pr_st_idx, ixcldliq, ixcldice, ixnumliq, &
        ixnumice, ixrain, ixsnow, state_loc%q, state_loc%t, state_loc%pmid, state_loc%pdel, naai, naai_hom, mnuccdo, &
        rflx, sflx, qrout, qsout, prect, preci, rate1cld, vtrmc, tlat, qvlat, qcten, qiten, ncten, niten, alst_mic, &
        cmeliq, cmeiout, ast, cld, concld, mnuccdohet, mgflxprc, mgflxsnw, mgmrprc, mgmrsnw, cvreffliq, cvreffice, &
        rate1ord_cw2pr_st, wsedl, CC_T, CC_qv, CC_ql, CC_qi, CC_nl, CC_ni, CC_qlst, qme, prec_pcw, snow_pcw, prec_sed, &
        snow_sed, prec_str, snow_str, icecldf, liqcldf, icinc, icwnc, iciwpst, iclwpst, icswp, cldfsnow, icimrst, &
        icwmrst, cldmax)

   ! ------------------------------------------------------ !
   ! ------------------------------------------------------ !
   ! All code from here to the end is on grid columns only  !
   ! ------------------------------------------------------ !
   ! ------------------------------------------------------ !

   ! Average the fields which are needed later in this paramterization to be on the grid
   if (use_subcol_microp) then
      call subcol_field_avg(prec_str,  ngrdcol, lchnk, prec_str_grid)
      call subcol_field_avg(iclwpst,   ngrdcol, lchnk, iclwpst_grid)
      call subcol_field_avg(cvreffliq, ngrdcol, lchnk, cvreffliq_grid)
      call subcol_field_avg(cvreffice, ngrdcol, lchnk, cvreffice_grid)
      call subcol_field_avg(mgflxprc,  ngrdcol, lchnk, mgflxprc_grid)
      call subcol_field_avg(mgflxsnw,  ngrdcol, lchnk, mgflxsnw_grid)
      call subcol_field_avg(qme,       ngrdcol, lchnk, qme_grid)
      call subcol_field_avg(nevapr,    ngrdcol, lchnk, nevapr_grid)
      call subcol_field_avg(prain,     ngrdcol, lchnk, prain_grid)
      call subcol_field_avg(evapsnow,  ngrdcol, lchnk, evpsnow_st_grid)

      if (micro_mg_version == 1 .and. micro_mg_sub_version == 0) then
         call subcol_field_avg(am_evp_st, ngrdcol, lchnk, am_evp_st_grid)
      end if

      ! Average fields which are not in pbuf
      call subcol_field_avg(qrout,     ngrdcol, lchnk, qrout_grid)
      call subcol_field_avg(qsout,     ngrdcol, lchnk, qsout_grid)
      call subcol_field_avg(nsout,     ngrdcol, lchnk, nsout_grid)
      call subcol_field_avg(nrout,     ngrdcol, lchnk, nrout_grid)
      call subcol_field_avg(cld,       ngrdcol, lchnk, cld_grid)
      call subcol_field_avg(qcreso,    ngrdcol, lchnk, qcreso_grid)
      call subcol_field_avg(melto,     ngrdcol, lchnk, melto_grid)
      call subcol_field_avg(mnuccco,   ngrdcol, lchnk, mnuccco_grid)
      call subcol_field_avg(mnuccto,   ngrdcol, lchnk, mnuccto_grid)
      call subcol_field_avg(bergo,     ngrdcol, lchnk, bergo_grid)
      call subcol_field_avg(homoo,     ngrdcol, lchnk, homoo_grid)
      call subcol_field_avg(msacwio,   ngrdcol, lchnk, msacwio_grid)
      call subcol_field_avg(psacwso,   ngrdcol, lchnk, psacwso_grid)
      call subcol_field_avg(bergso,    ngrdcol, lchnk, bergso_grid)
      call subcol_field_avg(cmeiout,   ngrdcol, lchnk, cmeiout_grid)
      call subcol_field_avg(qireso,    ngrdcol, lchnk, qireso_grid)
      call subcol_field_avg(prcio,     ngrdcol, lchnk, prcio_grid)
      call subcol_field_avg(praio,     ngrdcol, lchnk, praio_grid)
      call subcol_field_avg(icwmrst,   ngrdcol, lchnk, icwmrst_grid)
      call subcol_field_avg(icimrst,   ngrdcol, lchnk, icimrst_grid)
      call subcol_field_avg(liqcldf,   ngrdcol, lchnk, liqcldf_grid)
      call subcol_field_avg(icecldf,   ngrdcol, lchnk, icecldf_grid)
      call subcol_field_avg(icwnc,     ngrdcol, lchnk, icwnc_grid)
      call subcol_field_avg(icinc,     ngrdcol, lchnk, icinc_grid)
      call subcol_field_avg(state_loc%pdel,            ngrdcol, lchnk, pdel_grid)
      call subcol_field_avg(prao,      ngrdcol, lchnk, prao_grid)
      call subcol_field_avg(prco,      ngrdcol, lchnk, prco_grid)

      call subcol_field_avg(state_loc%q(:,:,ixnumliq), ngrdcol, lchnk, nc_grid)
      call subcol_field_avg(state_loc%q(:,:,ixnumice), ngrdcol, lchnk, ni_grid)

      if (micro_mg_version > 1) then
         call subcol_field_avg(cldmax,    ngrdcol, lchnk, cldmax_grid)

         call subcol_field_avg(state_loc%q(:,:,ixrain),    ngrdcol, lchnk, qr_grid)
         call subcol_field_avg(state_loc%q(:,:,ixnumrain), ngrdcol, lchnk, nr_grid)
         call subcol_field_avg(state_loc%q(:,:,ixsnow),    ngrdcol, lchnk, qs_grid)
         call subcol_field_avg(state_loc%q(:,:,ixnumsnow), ngrdcol, lchnk, ns_grid)
      end if

   else
      ! These pbuf fields need to be assigned.  There is no corresponding subcol_field_avg
      ! as they are reset before being used, so it would be a needless calculation
      lambdac_grid    => lambdac
      mu_grid         => mu
      rel_grid        => rel
      rei_grid        => rei
      dei_grid        => dei
      des_grid        => des

      ! fields already on grids, so just assign
      prec_str_grid   => prec_str
      iclwpst_grid    => iclwpst
      cvreffliq_grid  => cvreffliq
      cvreffice_grid  => cvreffice
      mgflxprc_grid   => mgflxprc
      mgflxsnw_grid   => mgflxsnw
      qme_grid        => qme
      nevapr_grid     => nevapr
      prain_grid      => prain

      if (micro_mg_version == 1 .and. micro_mg_sub_version == 0) then
         am_evp_st_grid  = am_evp_st
      end if

      evpsnow_st_grid = evapsnow
      qrout_grid      = qrout
      qsout_grid      = qsout
      nsout_grid      = nsout
      nrout_grid      = nrout
      cld_grid        = cld
      qcreso_grid     = qcreso
      melto_grid      = melto
      mnuccco_grid    = mnuccco
      mnuccto_grid    = mnuccto
      bergo_grid      = bergo
      homoo_grid      = homoo
      msacwio_grid    = msacwio
      psacwso_grid    = psacwso
      bergso_grid     = bergso
      cmeiout_grid    = cmeiout
      qireso_grid     = qireso
      prcio_grid      = prcio
      praio_grid      = praio
      icwmrst_grid    = icwmrst
      icimrst_grid    = icimrst
      liqcldf_grid    = liqcldf
      icecldf_grid    = icecldf
      icwnc_grid      = icwnc
      icinc_grid      = icinc
      pdel_grid       = state_loc%pdel
      prao_grid       = prao
      prco_grid       = prco

      nc_grid = state_loc%q(:,:,ixnumliq)
      ni_grid = state_loc%q(:,:,ixnumice)

      if (micro_mg_version > 1) then
         cldmax_grid = cldmax

         qr_grid = state_loc%q(:,:,ixrain)
         nr_grid = state_loc%q(:,:,ixnumrain)
         qs_grid = state_loc%q(:,:,ixsnow)
         ns_grid = state_loc%q(:,:,ixnumsnow)
      end if

   end if

   ! If on subcolumns, average the rest of the pbuf fields which were modified on subcolumns but are not used further in
   ! this parameterization  (no need to assign in the non-subcolumn case -- the else step)
   if (use_subcol_microp) then
      call subcol_field_avg(snow_str,    ngrdcol, lchnk, snow_str_grid)
      call subcol_field_avg(prec_pcw,    ngrdcol, lchnk, prec_pcw_grid)
      call subcol_field_avg(snow_pcw,    ngrdcol, lchnk, snow_pcw_grid)
      call subcol_field_avg(prec_sed,    ngrdcol, lchnk, prec_sed_grid)
      call subcol_field_avg(snow_sed,    ngrdcol, lchnk, snow_sed_grid)
      call subcol_field_avg(cldo,        ngrdcol, lchnk, cldo_grid)
      call subcol_field_avg(mgmrprc,     ngrdcol, lchnk, mgmrprc_grid)
      call subcol_field_avg(mgmrsnw,     ngrdcol, lchnk, mgmrsnw_grid)
      call subcol_field_avg(wsedl,       ngrdcol, lchnk, wsedl_grid)
      call subcol_field_avg(cc_t,        ngrdcol, lchnk, cc_t_grid)
      call subcol_field_avg(cc_qv,       ngrdcol, lchnk, cc_qv_grid)
      call subcol_field_avg(cc_ql,       ngrdcol, lchnk, cc_ql_grid)
      call subcol_field_avg(cc_qi,       ngrdcol, lchnk, cc_qi_grid)
      call subcol_field_avg(cc_nl,       ngrdcol, lchnk, cc_nl_grid)
      call subcol_field_avg(cc_ni,       ngrdcol, lchnk, cc_ni_grid)
      call subcol_field_avg(cc_qlst,     ngrdcol, lchnk, cc_qlst_grid)
      call subcol_field_avg(iciwpst,     ngrdcol, lchnk, iciwpst_grid)
      call subcol_field_avg(icswp,       ngrdcol, lchnk, icswp_grid)
      call subcol_field_avg(cldfsnow,    ngrdcol, lchnk, cldfsnow_grid)

      if (rate1_cw2pr_st_idx > 0) then
         call subcol_field_avg(rate1ord_cw2pr_st,    ngrdcol, lchnk, rate1ord_cw2pr_st_grid)
      end if

   end if

   call micro_mg_cam_select_tail_shell_impl()

   if (use_native_tail_shell_impl) then

      !----------------------------------------
      !water tracers/isotopes   (on gridlevel)
      !----------------------------------------

      ! Convert fields to grid level early, that are needed by water tracers

      if (trace_water) then

         ! Average isotope fields to the grid level , so they can be operated on
         if (use_subcol_microp) then
            !
            ! EBK Apr/21/2015
            ! In order to run on sub-columns all fields would need to
            ! be averaged to the grid level like so...
            !call subcol_field_avg(preo,      ngrdcol, lchnk, preo_grid)
            ! For the list of fields, see the "else" statement
            ! Also "state", "ptend" and "pbuf" would all be on
            ! sub-columns and would need to be averaged to grid level
            ! and then copied to the sub-column level. Or all the operations
            ! below would need to be on the sub-column level rather than
            ! grid level. Some of this would also require changes in water_tracers.F90.
            !
            ! (For now just terminate early)
            call endrun(subname // ':: ERROR water tracers are NOT configured to work with subcolumns')
         else
            preo_grid      => preo
            prdso_grid     => prdso
            frzro_grid     => frzro
            meltso_grid    => meltso
            wtfc_grid      => wtfc
            wtfi_grid      => wtfi
            wtprelat_grid  => wtprelat
            wtpostlat_grid => wtpostlat

            mnuccro_grid  => mnuccro
            pracso_grid   => pracso
            qcsedten_grid => qcsedten
            qisedten_grid => qisedten
            alst_mic_grid => alst_mic
            aist_mic_grid => aist_mic
         end if

         call micro_mg_cam_select_wtrc_shell_impl()

         if (use_native_wtrc_shell_impl) then
            call wtrc_init_rates(top_lev, pre_rates_grid)
            call wtrc_init_rates(top_lev, post_rates_grid)

            do k = top_lev, pver
               do i = 1, pcols
                  pcmei_grid(i,k) = 0._r8
                  ncmei_grid(i,k) = 0._r8
                  pmelts_grid(i,k) = 0._r8
                  nmelts_grid(i,k) = 0._r8
                  do m = 1, pwtype
                     sed_rates_grid(i,k,m) = 0._r8
                  end do
               end do
            end do

            do k = top_lev, pver
               do i = 1, ncol
                  if (cmeiout_grid(i,k) < 0._r8) then
                     ncmei_grid(i,k) = cmeiout_grid(i,k)
                  else
                     pcmei_grid(i,k) = cmeiout_grid(i,k)
                  end if
                  if (meltso_grid(i,k) < 0._r8) then
                     nmelts_grid(i,k) = meltso_grid(i,k)
                  else
                     pmelts_grid(i,k) = meltso_grid(i,k)
                  end if
                  sed_rates_grid(i,k,iwtliq) = qcsedten_grid(i,k)
                  sed_rates_grid(i,k,iwtice) = qisedten_grid(i,k)
               end do
            end do

            do k = top_lev, pver
               do i = 1, ncol
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtvap, iwtice, iwtvap, pcmei_grid(i,k))
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtvap, iwtice, iwtice, ncmei_grid(i,k))
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtvap, iwtstrain, iwtstrain, preo_grid(i,k))
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtvap, iwtstsnow, iwtstsnow, prdso_grid(i,k))

                  rate_local = mnuccco_grid(i,k) + mnuccto_grid(i,k)
                  rate_local = rate_local + msacwio_grid(i,k)
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtliq, iwtice, iwtliq, rate_local)

                  rate_local = prao_grid(i,k) + prco_grid(i,k)
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtliq, iwtstrain, iwtliq, rate_local)

                  call wtrc_add_rate(pre_rates_grid, i, k, iwtliq, iwtstsnow, iwtliq, psacwso_grid(i,k))
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtliq, iwtliq, iwtliq, bergo_grid(i,k))
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtice, iwtice, iwtice, bergso_grid(i,k))

                  rate_local = praio_grid(i,k) + prcio_grid(i,k)
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtice, iwtstsnow, iwtice, rate_local)

                  rate_local = pracso_grid(i,k) + mnuccro_grid(i,k)
                  call wtrc_add_rate(pre_rates_grid, i, k, iwtstrain, iwtstsnow, iwtstrain, rate_local)

                  call wtrc_add_rate(post_rates_grid, i, k, iwtvap, iwtliq, iwtvap, qcreso_grid(i,k))
                  call wtrc_add_rate(post_rates_grid, i, k, iwtvap, iwtice, iwtvap, qireso_grid(i,k))
                  call wtrc_add_rate(post_rates_grid, i, k, iwtliq, iwtice, iwtliq, homoo_grid(i,k))
                  call wtrc_add_rate(post_rates_grid, i, k, iwtice, iwtliq, iwtice, melto_grid(i,k))
               end do
            end do
         else
            call micro_mg_cam_wtrc_shell_codon_wrap(ncol, preo_grid, prdso_grid, cmeiout_grid, meltso_grid, qcsedten_grid, &
                 qisedten_grid, mnuccco_grid, mnuccto_grid, msacwio_grid, prao_grid, prco_grid, psacwso_grid, bergo_grid, &
                 bergso_grid, praio_grid, prcio_grid, pracso_grid, mnuccro_grid, qcreso_grid, qireso_grid, homoo_grid, &
                 melto_grid, pre_rates_grid, sed_rates_grid, post_rates_grid, pcmei_grid, ncmei_grid, pmelts_grid, nmelts_grid)
         end if

         ! Apply the microphysical process to the isotopes. rates.
         call wtrc_apply_rates(state, ptend, pbuf, top_lev, dtime, .true., pre_rates=pre_rates_grid, sed_rates=sed_rates_grid, &
              post_rates=post_rates_grid, do_stprecip=.true., liqcldf=alst_mic, icecldf=aist_mic, fc=wtfc_grid, fi=wtfi_grid, &
              prelat=wtprelat_grid, postlat=wtpostlat_grid, frzro=frzro_grid, meltso=meltso_grid)

      end if !water tracers

      !-------------------------------------
      ! ------------------------------------- !
      ! Size distribution calculation         !
      ! ------------------------------------- !

      ! Calculate rho (on subcolumns if turned on) for size distribution
      ! parameter calculations and average it if needed
      !
      ! State instead of state_loc to preserve answers for MG1 (and in any
      ! case, it is unlikely to make much difference).
      rho(:ncol,top_lev:) = state%pmid(:ncol,top_lev:) / &
           (rair*state%t(:ncol,top_lev:))
      if (use_subcol_microp) then
         call subcol_field_avg(rho, ngrdcol, lchnk, rho_grid)
      else
         rho_grid = rho
      end if

      call micro_mg_cam_select_diag_shell_impl()

      if (use_native_diag_shell_impl) then
         call micro_mg_cam_reff_calc(ngrdcol, micro_mg_version, rho_grid, icwmrst_grid, liqcldf_grid, nc_grid, qr_grid, nr_grid, &
              qs_grid, ns_grid, qrout_grid, nrout_grid, qsout_grid, nsout_grid, ni_grid, icecldf_grid, icimrst_grid, ast_grid, &
              mu_grid, lambdac_grid, rel_fn_grid, ncic_grid, rel_grid, drout2_grid, reff_rain_grid, des_grid, dsout2_grid, &
              reff_snow_grid, rei_grid, niic_grid, dei_grid, mgreffrain_grid, mgreffsnow_grid)

         ! ------------------------------------- !
         ! Precipitation efficiency Calculation  !
         ! ------------------------------------- !

         !-----------------------------------------------------------------------
         ! Liquid water path

         ! Compute liquid water paths, and column condensation
         minlwp = 0.01_r8        !minimum lwp threshold (kg/m3)

         call micro_mg_cam_grid_diag(ngrdcol, minlwp, iclwpst_grid, cld_grid, cmeliq_grid, pdel_grid, prec_str_grid, &
              acgcme_grid, acprecl_grid, acnum_grid, prao_grid, prco_grid, nc_grid, liqcldf_grid, icwmrst_grid, rel_grid, &
              icwnc_grid, icecldf_grid, icimrst_grid, rei_grid, icinc_grid, nevapr_grid, evpsnow_st_grid, tgliqwp_grid, &
              tgcmeliq_grid, pe_grid, tpr_grid, pefrac_grid, vprao_grid, vprco_grid, racau_grid, cnt_grid, cdnumc_grid, &
              efcout_grid, efiout_grid, ncout_grid, niout_grid, freql_grid, freqi_grid, icwmrst_grid_out, icimrst_grid_out, &
              fcti_grid, fctl_grid, ctrel_grid, ctrei_grid, ctnl_grid, ctni_grid, evprain_st_grid)
      else
         minlwp = 0.01_r8        !minimum lwp threshold (kg/m3)

         ! Keep the liquid effective-radius branch native; the remaining active
         ! diagnostic tail is combined into one Codon shell.
         call micro_mg_cam_reff_liq_native(ngrdcol, rho_grid, icwmrst_grid, liqcldf_grid, nc_grid, mu_grid, lambdac_grid, &
              rel_fn_grid, ncic_grid, rel_grid)

         call micro_mg_cam_diag_shell_codon_wrap(ngrdcol, micro_mg_version, minlwp, rho_grid, icwmrst_grid, liqcldf_grid, &
              nc_grid, qr_grid, nr_grid, qs_grid, ns_grid, qrout_grid, nrout_grid, qsout_grid, nsout_grid, ni_grid, &
              icecldf_grid, icimrst_grid, ast_grid, mu_grid, lambdac_grid, rel_fn_grid, ncic_grid, rel_grid, drout2_grid, &
              reff_rain_grid, des_grid, dsout2_grid, reff_snow_grid, rei_grid, niic_grid, dei_grid, mgreffrain_grid, &
              mgreffsnow_grid, iclwpst_grid, cld_grid, cmeliq_grid, pdel_grid, prec_str_grid, acgcme_grid, acprecl_grid, &
              acnum_grid, prao_grid, prco_grid, icwnc_grid, icinc_grid, nevapr_grid, evpsnow_st_grid, tgliqwp_grid, &
              tgcmeliq_grid, pe_grid, tpr_grid, pefrac_grid, vprao_grid, vprco_grid, racau_grid, cnt_grid, cdnumc_grid, &
              efcout_grid, efiout_grid, ncout_grid, niout_grid, freql_grid, freqi_grid, icwmrst_grid_out, icimrst_grid_out, &
              fcti_grid, fctl_grid, ctrel_grid, ctrei_grid, ctnl_grid, ctni_grid, evprain_st_grid, qcreso_grid, melto_grid, &
              mnuccco_grid, mnuccto_grid, bergo_grid, homoo_grid, msacwio_grid, psacwso_grid, bergso_grid, cmeiout_grid, &
              qireso_grid, prcio_grid, praio_grid, budget_ftem_grid)
      end if

   else

      if (trace_water) then

         if (use_subcol_microp) then
            call endrun(subname // ':: ERROR water tracers are NOT configured to work with subcolumns')
         else
            preo_grid      => preo
            prdso_grid     => prdso
            frzro_grid     => frzro
            meltso_grid    => meltso
            wtfc_grid      => wtfc
            wtfi_grid      => wtfi
            wtprelat_grid  => wtprelat
            wtpostlat_grid => wtpostlat

            mnuccro_grid  => mnuccro
            pracso_grid   => pracso
            qcsedten_grid => qcsedten
            qisedten_grid => qisedten
            alst_mic_grid => alst_mic
            aist_mic_grid => aist_mic
         end if

         call micro_mg_cam_wtrc_shell_codon_wrap(ncol, preo_grid, prdso_grid, cmeiout_grid, meltso_grid, qcsedten_grid, &
              qisedten_grid, mnuccco_grid, mnuccto_grid, msacwio_grid, prao_grid, prco_grid, psacwso_grid, bergo_grid, &
              bergso_grid, praio_grid, prcio_grid, pracso_grid, mnuccro_grid, qcreso_grid, qireso_grid, homoo_grid, melto_grid, &
              pre_rates_grid, sed_rates_grid, post_rates_grid, pcmei_grid, ncmei_grid, pmelts_grid, nmelts_grid)

         call wtrc_apply_rates(state, ptend, pbuf, top_lev, dtime, .true., pre_rates=pre_rates_grid, sed_rates=sed_rates_grid, &
              post_rates=post_rates_grid, do_stprecip=.true., liqcldf=alst_mic, icecldf=aist_mic, fc=wtfc_grid, fi=wtfi_grid, &
              prelat=wtprelat_grid, postlat=wtpostlat_grid, frzro=frzro_grid, meltso=meltso_grid)

      end if

      rho(:ncol,top_lev:) = state%pmid(:ncol,top_lev:) / &
           (rair*state%t(:ncol,top_lev:))
      if (use_subcol_microp) then
         call subcol_field_avg(rho, ngrdcol, lchnk, rho_grid)
      else
         rho_grid = rho
      end if

      minlwp = 0.01_r8
      call micro_mg_cam_reff_liq_native(ngrdcol, rho_grid, icwmrst_grid, liqcldf_grid, nc_grid, mu_grid, lambdac_grid, &
           rel_fn_grid, ncic_grid, rel_grid)

      call micro_mg_cam_diag_shell_codon_wrap(ngrdcol, micro_mg_version, minlwp, rho_grid, icwmrst_grid, liqcldf_grid, &
           nc_grid, qr_grid, nr_grid, qs_grid, ns_grid, qrout_grid, nrout_grid, qsout_grid, nsout_grid, ni_grid, &
           icecldf_grid, icimrst_grid, ast_grid, mu_grid, lambdac_grid, rel_fn_grid, ncic_grid, rel_grid, drout2_grid, &
           reff_rain_grid, des_grid, dsout2_grid, reff_snow_grid, rei_grid, niic_grid, dei_grid, mgreffrain_grid, &
           mgreffsnow_grid, iclwpst_grid, cld_grid, cmeliq_grid, pdel_grid, prec_str_grid, acgcme_grid, acprecl_grid, &
           acnum_grid, prao_grid, prco_grid, icwnc_grid, icinc_grid, nevapr_grid, evpsnow_st_grid, tgliqwp_grid, &
           tgcmeliq_grid, pe_grid, tpr_grid, pefrac_grid, vprao_grid, vprco_grid, racau_grid, cnt_grid, cdnumc_grid, &
           efcout_grid, efiout_grid, ncout_grid, niout_grid, freql_grid, freqi_grid, icwmrst_grid_out, icimrst_grid_out, &
           fcti_grid, fctl_grid, ctrel_grid, ctrei_grid, ctnl_grid, ctni_grid, evprain_st_grid, qcreso_grid, melto_grid, &
           mnuccco_grid, mnuccto_grid, bergo_grid, homoo_grid, msacwio_grid, psacwso_grid, bergso_grid, cmeiout_grid, &
           qireso_grid, prcio_grid, praio_grid, budget_ftem_grid)

   end if

   ! Assign the values to the pbuf pointers if they exist in pbuf
   call micro_mg_cam_pbuf_copy(ncol, qrain_idx > 0, qsnow_idx > 0, nrain_idx > 0, nsnow_idx > 0, &
        qrout_grid, qsout_grid, nrout_grid, nsout_grid, qrout_grid_ptr, qsout_grid_ptr, nrout_grid_ptr, nsout_grid_ptr)

   ! --------------------------------------------- !
   ! General outfield calls for microphysics       !
   ! --------------------------------------------- !

   if (use_native_diag_shell_impl) then
      call micro_mg_cam_budget_diag(1, ngrdcol, qcreso_grid, melto_grid, mnuccco_grid, mnuccto_grid, bergo_grid, homoo_grid, &
           msacwio_grid, prao_grid, prco_grid, psacwso_grid, bergso_grid, cmeiout_grid, qireso_grid, prcio_grid, praio_grid, &
           budget_ftem_grid(:,:,1))

      call micro_mg_cam_budget_diag(2, ngrdcol, qcreso_grid, melto_grid, mnuccco_grid, mnuccto_grid, bergo_grid, homoo_grid, &
           msacwio_grid, prao_grid, prco_grid, psacwso_grid, bergso_grid, cmeiout_grid, qireso_grid, prcio_grid, praio_grid, &
           budget_ftem_grid(:,:,2))

      call micro_mg_cam_budget_diag(3, ngrdcol, qcreso_grid, melto_grid, mnuccco_grid, mnuccto_grid, bergo_grid, homoo_grid, &
           msacwio_grid, prao_grid, prco_grid, psacwso_grid, bergso_grid, cmeiout_grid, qireso_grid, prcio_grid, praio_grid, &
           budget_ftem_grid(:,:,3))

      call micro_mg_cam_budget_diag(4, ngrdcol, qcreso_grid, melto_grid, mnuccco_grid, mnuccto_grid, bergo_grid, homoo_grid, &
           msacwio_grid, prao_grid, prco_grid, psacwso_grid, bergso_grid, cmeiout_grid, qireso_grid, prcio_grid, praio_grid, &
           budget_ftem_grid(:,:,4))

      call micro_mg_cam_budget_diag(5, ngrdcol, qcreso_grid, melto_grid, mnuccco_grid, mnuccto_grid, bergo_grid, homoo_grid, &
           msacwio_grid, prao_grid, prco_grid, psacwso_grid, bergso_grid, cmeiout_grid, qireso_grid, prcio_grid, praio_grid, &
           budget_ftem_grid(:,:,5))

      call micro_mg_cam_budget_diag(6, ngrdcol, qcreso_grid, melto_grid, mnuccco_grid, mnuccto_grid, bergo_grid, homoo_grid, &
           msacwio_grid, prao_grid, prco_grid, psacwso_grid, bergso_grid, cmeiout_grid, qireso_grid, prcio_grid, praio_grid, &
           budget_ftem_grid(:,:,6))
   end if
   call outfld( 'MPDW2V', budget_ftem_grid(:,:,1), pcols, lchnk)
   call outfld( 'MPDW2I', budget_ftem_grid(:,:,2), pcols, lchnk)
   call outfld( 'MPDW2P', budget_ftem_grid(:,:,3), pcols, lchnk)
   call outfld( 'MPDI2V', budget_ftem_grid(:,:,4), pcols, lchnk)
   call outfld( 'MPDI2W', budget_ftem_grid(:,:,5), pcols, lchnk)
   call outfld( 'MPDI2P', budget_ftem_grid(:,:,6), pcols, lchnk)

   ! Output fields which have not been averaged already, averaging if use_subcol_microp is true
   call outfld('MPICLWPI',    iclwpi,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPICIWPI',    iciwpi,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('REFL',        refl,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('AREFL',       arefl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('AREFLZ',      areflz,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FREFL',       frefl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('CSRFL',       csrfl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('ACSRFL',      acsrfl,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FCSRFL',      fcsrfl,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('RERCLD',      rercld,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('NCAL',        ncal,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('NCAI',        ncai,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('AQRAIN',      qrout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('AQSNOW',      qsout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('ANRAIN',      nrout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('ANSNOW',      nsout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FREQR',       freqr,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FREQS',       freqs,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDT',        tlat,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDQ',        qvlat,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDLIQ',      qcten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDICE',      qiten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('EVAPSNOW',    evapsnow,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QCSEVAP',     qcsevap,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QISEVAP',     qisevap,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QVRES',       qvres,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('VTRMC',       vtrmc,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('VTRMI',       vtrmi,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QCSEDTEN',    qcsedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QISEDTEN',    qisedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   if (micro_mg_version > 1) then
      call outfld('QRSEDTEN',    qrsedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('QSSEDTEN',    qssedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   end if
   call outfld('MNUCCDO',     mnuccdo,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MNUCCDOhet',  mnuccdohet,  psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MNUCCRO',     mnuccro,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('PRACSO',      pracso ,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MELTSDT',     meltsdt,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FRZRDT',      frzrdt ,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FICE',        nfice,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)

   if (micro_mg_version > 1) then
      call outfld('UMR',      umr,         psetcols, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('UMS',      ums,         psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   end if

   if (.not. (micro_mg_version == 1 .and. micro_mg_sub_version == 0)) then
      call outfld('QCRAT',    qcrat,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   end if

   ! Example subcolumn outfld call
   if (use_subcol_microp) then
      call outfld('FICE_SCOL',   nfice,       psubcols*pcols, lchnk)
   end if

   ! Output fields which are already on the grid
   call outfld('QRAIN',       qrout_grid,       pcols, lchnk)
   call outfld('QSNOW',       qsout_grid,       pcols, lchnk)
   call outfld('NRAIN',       nrout_grid,       pcols, lchnk)
   call outfld('NSNOW',       nsout_grid,       pcols, lchnk)
   call outfld('CV_REFFLIQ',  cvreffliq_grid,   pcols, lchnk)
   call outfld('CV_REFFICE',  cvreffice_grid,   pcols, lchnk)
   call outfld('LS_FLXPRC',   mgflxprc_grid,    pcols, lchnk)
   call outfld('LS_FLXSNW',   mgflxsnw_grid,    pcols, lchnk)
   call outfld('CME',         qme_grid,         pcols, lchnk)
   call outfld('PRODPREC',    prain_grid,       pcols, lchnk)
   call outfld('EVAPPREC',    nevapr_grid,      pcols, lchnk)
   call outfld('QCRESO',      qcreso_grid,      pcols, lchnk)
   call outfld('LS_REFFRAIN', mgreffrain_grid,  pcols, lchnk)
   call outfld('LS_REFFSNOW', mgreffsnow_grid,  pcols, lchnk)
   call outfld('DSNOW',       des_grid,         pcols, lchnk)
   call outfld('ADRAIN',      drout2_grid,      pcols, lchnk)
   call outfld('ADSNOW',      dsout2_grid,      pcols, lchnk)
   call outfld('PE',          pe_grid,          pcols, lchnk)
   call outfld('PEFRAC',      pefrac_grid,      pcols, lchnk)
   call outfld('APRL',        tpr_grid,         pcols, lchnk)
   call outfld('VPRAO',       vprao_grid,       pcols, lchnk)
   call outfld('VPRCO',       vprco_grid,       pcols, lchnk)
   call outfld('RACAU',       racau_grid,       pcols, lchnk)
   call outfld('AREL',        efcout_grid,      pcols, lchnk)
   call outfld('AREI',        efiout_grid,      pcols, lchnk)
   call outfld('AWNC' ,       ncout_grid,       pcols, lchnk)
   call outfld('AWNI' ,       niout_grid,       pcols, lchnk)
   call outfld('FREQL',       freql_grid,       pcols, lchnk)
   call outfld('FREQI',       freqi_grid,       pcols, lchnk)
   call outfld('ACTREL',      ctrel_grid,       pcols, lchnk)
   call outfld('ACTREI',      ctrei_grid,       pcols, lchnk)
   call outfld('ACTNL',       ctnl_grid,        pcols, lchnk)
   call outfld('ACTNI',       ctni_grid,        pcols, lchnk)
   call outfld('FCTL',        fctl_grid,        pcols, lchnk)
   call outfld('FCTI',        fcti_grid,        pcols, lchnk)
   call outfld('ICINC',       icinc_grid,       pcols, lchnk)
   call outfld('ICWNC',       icwnc_grid,       pcols, lchnk)
   call outfld('EFFLIQ_IND',  rel_fn_grid,      pcols, lchnk)
   call outfld('CDNUMC',      cdnumc_grid,      pcols, lchnk)
   call outfld('REL',         rel_grid,         pcols, lchnk)
   call outfld('REI',         rei_grid,         pcols, lchnk)
   call outfld('ICIMRST',     icimrst_grid_out, pcols, lchnk)
   call outfld('ICWMRST',     icwmrst_grid_out, pcols, lchnk)
   call outfld('CMEIOUT',     cmeiout_grid,     pcols, lchnk)
   call outfld('PRAO',        prao_grid,        pcols, lchnk)
   call outfld('PRCO',        prco_grid,        pcols, lchnk)
   call outfld('MNUCCCO',     mnuccco_grid,     pcols, lchnk)
   call outfld('MNUCCTO',     mnuccto_grid,     pcols, lchnk)
   call outfld('MSACWIO',     msacwio_grid,     pcols, lchnk)
   call outfld('PSACWSO',     psacwso_grid,     pcols, lchnk)
   call outfld('BERGSO',      bergso_grid,      pcols, lchnk)
   call outfld('BERGO',       bergo_grid,       pcols, lchnk)
   call outfld('MELTO',       melto_grid,       pcols, lchnk)
   call outfld('HOMOO',       homoo_grid,       pcols, lchnk)
   call outfld('PRCIO',       prcio_grid,       pcols, lchnk)
   call outfld('PRAIO',       praio_grid,       pcols, lchnk)
   call outfld('QIRESO',      qireso_grid,      pcols, lchnk)

   ! Output fields for the water traacers.
   call wtrc_output_precip(state_loc, pbuf)

   ! ptend_loc is deallocated in physics_update above
   call physics_state_dealloc(state_loc)

end subroutine micro_mg_cam_tend

subroutine micro_mg_cam_select_tail_shell_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tail_shell_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MICRO_MG_CAM_TAIL_SHELL_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tail_shell_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tail_shell_impl = .false.
  end if

  tail_shell_impl_selected = .true.

  if (use_native_tail_shell_impl) then
     write(iulog,*) 'micro_mg_cam_tail_shell implementation = native'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_TAIL_SHELL_PROOF_FILE', &
          'micro_mg_cam_tail_shell implementation = native')
  else
     write(iulog,*) 'micro_mg_cam_tail_shell implementation = codon'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_TAIL_SHELL_PROOF_FILE', &
          'micro_mg_cam_tail_shell implementation = codon')
  end if
  call flush(iulog)

end subroutine micro_mg_cam_select_tail_shell_impl

subroutine micro_mg_cam_select_diag_shell_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (diag_shell_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MICRO_MG_CAM_DIAG_SHELL_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_diag_shell_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_diag_shell_impl = .false.
  end if

  diag_shell_impl_selected = .true.

  if (use_native_diag_shell_impl) then
     write(iulog,*) 'micro_mg_cam_diag_shell implementation = native'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_DIAG_SHELL_PROOF_FILE', &
          'micro_mg_cam_diag_shell implementation = native')
  else
     write(iulog,*) 'micro_mg_cam_diag_shell implementation = codon'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_DIAG_SHELL_PROOF_FILE', &
          'micro_mg_cam_diag_shell implementation = codon')
  end if
  call flush(iulog)

end subroutine micro_mg_cam_select_diag_shell_impl

subroutine micro_mg_cam_diag_shell_codon_wrap(ngrdcol_local, micro_mg_version_local, minlwp_local, rho_grid_local, &
     icwmrst_grid_local, liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, &
     qrout_grid_local, nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, &
     icimrst_grid_local, ast_grid_local, mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local, &
     drout2_grid_local, reff_rain_grid_local, des_grid_local, dsout2_grid_local, reff_snow_grid_local, rei_grid_local, &
     niic_grid_local, dei_grid_local, mgreffrain_grid_local, mgreffsnow_grid_local, iclwpst_grid_local, cld_grid_local, &
     cmeliq_grid_local, pdel_grid_local, prec_str_grid_local, acgcme_grid_local, acprecl_grid_local, acnum_grid_local, &
     prao_grid_local, prco_grid_local, icwnc_grid_local, icinc_grid_local, nevapr_grid_local, evpsnow_st_grid_local, &
     tgliqwp_grid_local, tgcmeliq_grid_local, pe_grid_local, tpr_grid_local, pefrac_grid_local, vprao_grid_local, &
     vprco_grid_local, racau_grid_local, cnt_grid_local, cdnumc_grid_local, efcout_grid_local, efiout_grid_local, &
     ncout_grid_local, niout_grid_local, freql_grid_local, freqi_grid_local, icwmrst_grid_out_local, icimrst_grid_out_local, &
     fcti_grid_local, fctl_grid_local, ctrel_grid_local, ctrei_grid_local, ctnl_grid_local, ctni_grid_local, &
     evprain_st_grid_local, qcreso_grid_local, melto_grid_local, mnuccco_grid_local, mnuccto_grid_local, bergo_grid_local, &
     homoo_grid_local, msacwio_grid_local, psacwso_grid_local, bergso_grid_local, cmeiout_grid_local, qireso_grid_local, &
     prcio_grid_local, praio_grid_local, budget_ftem_grid_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use micro_mg_utils, only: mg_liq_props, mg_ice_props, qsmall, mincld, rhosn, rhoi, rhow, rhows
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ngrdcol_local, micro_mg_version_local
  real(r8), intent(in) :: minlwp_local
  real(r8), parameter :: dcon_local = 25.e-6_r8
  real(r8), parameter :: mucon_local = 5.3_r8
  real(r8), parameter :: deicon_local = 50._r8
  real(r8), target, intent(in) :: rho_grid_local(pcols,pver), icwmrst_grid_local(pcols,pver)
  real(r8), target, intent(in) :: liqcldf_grid_local(pcols,pver), nc_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qr_grid_local(pcols,pver), nr_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qs_grid_local(pcols,pver), ns_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qrout_grid_local(pcols,pver), nrout_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qsout_grid_local(pcols,pver), nsout_grid_local(pcols,pver)
  real(r8), target, intent(in) :: ni_grid_local(pcols,pver), icecldf_grid_local(pcols,pver)
  real(r8), target, intent(in) :: icimrst_grid_local(pcols,pver), ast_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: mu_grid_local(pcols,pver), lambdac_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: rel_fn_grid_local(pcols,pver), ncic_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: rel_grid_local(pcols,pver), drout2_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: reff_rain_grid_local(pcols,pver), des_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: dsout2_grid_local(pcols,pver), reff_snow_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: rei_grid_local(pcols,pver), niic_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: dei_grid_local(pcols,pver), mgreffrain_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: mgreffsnow_grid_local(pcols,pver)
  real(r8), target, intent(in) :: iclwpst_grid_local(pcols,pver), cld_grid_local(pcols,pver)
  real(r8), target, intent(in) :: cmeliq_grid_local(pcols,pver), pdel_grid_local(pcols,pver)
  real(r8), target, intent(in) :: prec_str_grid_local(pcols), prao_grid_local(pcols,pver), prco_grid_local(pcols,pver)
  real(r8), target, intent(in) :: icwnc_grid_local(pcols,pver), icinc_grid_local(pcols,pver)
  real(r8), target, intent(in) :: nevapr_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: evpsnow_st_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: acgcme_grid_local(:), acprecl_grid_local(:)
  integer, target, intent(inout) :: acnum_grid_local(:)
  real(r8), target, intent(inout) :: tgliqwp_grid_local(pcols), tgcmeliq_grid_local(pcols), pe_grid_local(pcols)
  real(r8), target, intent(inout) :: tpr_grid_local(pcols), pefrac_grid_local(pcols), vprao_grid_local(pcols)
  real(r8), target, intent(inout) :: vprco_grid_local(pcols), racau_grid_local(pcols), cdnumc_grid_local(pcols)
  integer, target, intent(inout) :: cnt_grid_local(pcols)
  real(r8), target, intent(inout) :: efcout_grid_local(pcols,pver), efiout_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: ncout_grid_local(pcols,pver), niout_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: freql_grid_local(pcols,pver), freqi_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: icwmrst_grid_out_local(pcols,pver), icimrst_grid_out_local(pcols,pver)
  real(r8), target, intent(inout) :: fcti_grid_local(pcols), fctl_grid_local(pcols), ctrel_grid_local(pcols)
  real(r8), target, intent(inout) :: ctrei_grid_local(pcols), ctnl_grid_local(pcols), ctni_grid_local(pcols)
  real(r8), target, intent(inout) :: evprain_st_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qcreso_grid_local(pcols,pver), melto_grid_local(pcols,pver)
  real(r8), target, intent(in) :: mnuccco_grid_local(pcols,pver), mnuccto_grid_local(pcols,pver)
  real(r8), target, intent(in) :: bergo_grid_local(pcols,pver), homoo_grid_local(pcols,pver)
  real(r8), target, intent(in) :: msacwio_grid_local(pcols,pver), psacwso_grid_local(pcols,pver)
  real(r8), target, intent(in) :: bergso_grid_local(pcols,pver), cmeiout_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qireso_grid_local(pcols,pver), prcio_grid_local(pcols,pver), praio_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: budget_ftem_grid_local(pcols,pver,6)

  interface
     subroutine micro_mg_cam_diag_shell_codon(ngrdcol_c, pcols_c, pver_c, top_lev_c, micro_mg_version_c, qsmall_c, &
          mincld_c, liq_rho_c, liq_eff_dim_c, liq_min_mean_mass_c, ice_eff_dim_c, ice_shape_coef_c, ice_lambda_lo_c, &
          ice_lambda_hi_c, ice_min_mean_mass_c, rhosn_c, rhoi_c, rhow_c, rhows_c, mucon_c, dcon_c, deicon_c, minlwp_c, &
          gravit_c, rhoh2o_c, rho_grid_p, icwmrst_grid_p, liqcldf_grid_p, nc_grid_p, qr_grid_p, nr_grid_p, qs_grid_p, &
          ns_grid_p, qrout_grid_p, nrout_grid_p, qsout_grid_p, nsout_grid_p, ni_grid_p, icecldf_grid_p, icimrst_grid_p, &
          ast_grid_p, mu_grid_p, lambdac_grid_p, rel_fn_grid_p, ncic_grid_p, rel_grid_p, drout2_grid_p, &
          reff_rain_grid_p, des_grid_p, dsout2_grid_p, reff_snow_grid_p, rei_grid_p, niic_grid_p, dei_grid_p, &
          mgreffrain_grid_p, mgreffsnow_grid_p, iclwpst_grid_p, cld_grid_p, cmeliq_grid_p, pdel_grid_p, prec_str_grid_p, &
          acgcme_grid_p, acprecl_grid_p, acnum_grid_p, prao_grid_p, prco_grid_p, icwnc_grid_p, icinc_grid_p, nevapr_grid_p, &
          evpsnow_st_grid_p, tgliqwp_grid_p, tgcmeliq_grid_p, pe_grid_p, tpr_grid_p, pefrac_grid_p, vprao_grid_p, &
          vprco_grid_p, racau_grid_p, cnt_grid_p, cdnumc_grid_p, efcout_grid_p, efiout_grid_p, ncout_grid_p, niout_grid_p, &
          freql_grid_p, freqi_grid_p, icwmrst_grid_out_p, icimrst_grid_out_p, fcti_grid_p, fctl_grid_p, ctrel_grid_p, &
          ctrei_grid_p, ctnl_grid_p, ctni_grid_p, evprain_st_grid_p, qcreso_grid_p, melto_grid_p, mnuccco_grid_p, &
          mnuccto_grid_p, bergo_grid_p, homoo_grid_p, msacwio_grid_p, psacwso_grid_p, bergso_grid_p, cmeiout_grid_p, &
          qireso_grid_p, prcio_grid_p, praio_grid_p, budget_ftem_grid_p) bind(c, name="micro_mg_cam_diag_shell_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ngrdcol_c, pcols_c, pver_c, top_lev_c, micro_mg_version_c
       real(c_double), value :: qsmall_c, mincld_c, liq_rho_c, liq_eff_dim_c, liq_min_mean_mass_c
       real(c_double), value :: ice_eff_dim_c, ice_shape_coef_c, ice_lambda_lo_c, ice_lambda_hi_c
       real(c_double), value :: ice_min_mean_mass_c, rhosn_c, rhoi_c, rhow_c, rhows_c, mucon_c, dcon_c, deicon_c
       real(c_double), value :: minlwp_c, gravit_c, rhoh2o_c
       type(c_ptr), value :: rho_grid_p, icwmrst_grid_p, liqcldf_grid_p, nc_grid_p, qr_grid_p, nr_grid_p
       type(c_ptr), value :: qs_grid_p, ns_grid_p, qrout_grid_p, nrout_grid_p, qsout_grid_p, nsout_grid_p
       type(c_ptr), value :: ni_grid_p, icecldf_grid_p, icimrst_grid_p, ast_grid_p, mu_grid_p, lambdac_grid_p
       type(c_ptr), value :: rel_fn_grid_p, ncic_grid_p, rel_grid_p, drout2_grid_p, reff_rain_grid_p, des_grid_p
       type(c_ptr), value :: dsout2_grid_p, reff_snow_grid_p, rei_grid_p, niic_grid_p, dei_grid_p
       type(c_ptr), value :: mgreffrain_grid_p, mgreffsnow_grid_p, iclwpst_grid_p, cld_grid_p, cmeliq_grid_p
       type(c_ptr), value :: pdel_grid_p, prec_str_grid_p, acgcme_grid_p, acprecl_grid_p, acnum_grid_p, prao_grid_p
       type(c_ptr), value :: prco_grid_p, icwnc_grid_p, icinc_grid_p, nevapr_grid_p, evpsnow_st_grid_p
       type(c_ptr), value :: tgliqwp_grid_p, tgcmeliq_grid_p, pe_grid_p, tpr_grid_p, pefrac_grid_p, vprao_grid_p
       type(c_ptr), value :: vprco_grid_p, racau_grid_p, cnt_grid_p, cdnumc_grid_p, efcout_grid_p, efiout_grid_p
       type(c_ptr), value :: ncout_grid_p, niout_grid_p, freql_grid_p, freqi_grid_p, icwmrst_grid_out_p
       type(c_ptr), value :: icimrst_grid_out_p, fcti_grid_p, fctl_grid_p, ctrel_grid_p, ctrei_grid_p, ctnl_grid_p
       type(c_ptr), value :: ctni_grid_p, evprain_st_grid_p, qcreso_grid_p, melto_grid_p, mnuccco_grid_p
       type(c_ptr), value :: mnuccto_grid_p, bergo_grid_p, homoo_grid_p, msacwio_grid_p, psacwso_grid_p, bergso_grid_p
       type(c_ptr), value :: cmeiout_grid_p, qireso_grid_p, prcio_grid_p, praio_grid_p, budget_ftem_grid_p
     end subroutine micro_mg_cam_diag_shell_codon
  end interface

  call micro_mg_cam_diag_shell_codon(int(ngrdcol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), int(micro_mg_version_local, c_int64_t), qsmall, mincld, mg_liq_props%rho, &
       mg_liq_props%eff_dim, mg_liq_props%min_mean_mass, mg_ice_props%eff_dim, mg_ice_props%shape_coef, &
       mg_ice_props%lambda_bounds(1), mg_ice_props%lambda_bounds(2), mg_ice_props%min_mean_mass, rhosn, rhoi, rhow, rhows, &
       mucon_local, dcon_local, deicon_local, minlwp_local, gravit, rhoh2o, c_loc(rho_grid_local), c_loc(icwmrst_grid_local), &
       c_loc(liqcldf_grid_local), c_loc(nc_grid_local), c_loc(qr_grid_local), c_loc(nr_grid_local), c_loc(qs_grid_local), &
       c_loc(ns_grid_local), c_loc(qrout_grid_local), c_loc(nrout_grid_local), c_loc(qsout_grid_local), c_loc(nsout_grid_local), &
       c_loc(ni_grid_local), c_loc(icecldf_grid_local), c_loc(icimrst_grid_local), c_loc(ast_grid_local), c_loc(mu_grid_local), &
       c_loc(lambdac_grid_local), c_loc(rel_fn_grid_local), c_loc(ncic_grid_local), c_loc(rel_grid_local), &
       c_loc(drout2_grid_local), c_loc(reff_rain_grid_local), c_loc(des_grid_local), c_loc(dsout2_grid_local), &
       c_loc(reff_snow_grid_local), c_loc(rei_grid_local), c_loc(niic_grid_local), c_loc(dei_grid_local), &
       c_loc(mgreffrain_grid_local), c_loc(mgreffsnow_grid_local), c_loc(iclwpst_grid_local), c_loc(cld_grid_local), &
       c_loc(cmeliq_grid_local), c_loc(pdel_grid_local), c_loc(prec_str_grid_local), c_loc(acgcme_grid_local), &
       c_loc(acprecl_grid_local), c_loc(acnum_grid_local(1)), c_loc(prao_grid_local), c_loc(prco_grid_local), &
       c_loc(icwnc_grid_local), c_loc(icinc_grid_local), c_loc(nevapr_grid_local), c_loc(evpsnow_st_grid_local), &
       c_loc(tgliqwp_grid_local), c_loc(tgcmeliq_grid_local), c_loc(pe_grid_local), c_loc(tpr_grid_local), &
       c_loc(pefrac_grid_local), c_loc(vprao_grid_local), c_loc(vprco_grid_local), c_loc(racau_grid_local), &
       c_loc(cnt_grid_local(1)), c_loc(cdnumc_grid_local), c_loc(efcout_grid_local), c_loc(efiout_grid_local), &
       c_loc(ncout_grid_local), c_loc(niout_grid_local), c_loc(freql_grid_local), c_loc(freqi_grid_local), &
       c_loc(icwmrst_grid_out_local), c_loc(icimrst_grid_out_local), c_loc(fcti_grid_local), c_loc(fctl_grid_local), &
       c_loc(ctrel_grid_local), c_loc(ctrei_grid_local), c_loc(ctnl_grid_local), c_loc(ctni_grid_local), &
       c_loc(evprain_st_grid_local), c_loc(qcreso_grid_local), c_loc(melto_grid_local), c_loc(mnuccco_grid_local), &
       c_loc(mnuccto_grid_local), c_loc(bergo_grid_local), c_loc(homoo_grid_local), c_loc(msacwio_grid_local), &
       c_loc(psacwso_grid_local), c_loc(bergso_grid_local), c_loc(cmeiout_grid_local), c_loc(qireso_grid_local), &
       c_loc(prcio_grid_local), c_loc(praio_grid_local), c_loc(budget_ftem_grid_local))

end subroutine micro_mg_cam_diag_shell_codon_wrap

subroutine micro_mg_cam_select_pbuf_copy_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (pbuf_copy_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MICRO_MG_CAM_PBUF_COPY_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_pbuf_copy_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_pbuf_copy_impl = .false.
  end if

  pbuf_copy_impl_selected = .true.

  if (masterproc) then
     if (use_native_pbuf_copy_impl) then
        write(iulog,*) 'micro_mg_cam_pbuf_copy implementation = native'
        call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_PBUF_COPY_PROOF_FILE', &
             'micro_mg_cam_pbuf_copy implementation = native')
     else
        write(iulog,*) 'micro_mg_cam_pbuf_copy implementation = codon'
        call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_PBUF_COPY_PROOF_FILE', &
             'micro_mg_cam_pbuf_copy implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine micro_mg_cam_select_pbuf_copy_impl

subroutine micro_mg_cam_pbuf_copy_log_entered()

  if (pbuf_copy_entered_logged) return
  pbuf_copy_entered_logged = .true.

  if (masterproc) then
     write(iulog,*) 'micro_mg_cam_pbuf_copy entered (QRAIN/QSNOW/NRAIN/NSNOW pbuf copies direct = codon)'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_PBUF_COPY_PROOF_FILE', &
          'micro_mg_cam_pbuf_copy entered (QRAIN/QSNOW/NRAIN/NSNOW pbuf copies direct = codon)')
     call flush(iulog)
  end if

end subroutine micro_mg_cam_pbuf_copy_log_entered

subroutine micro_mg_cam_pbuf_copy(ncol_local, copy_qrain_local, copy_qsnow_local, copy_nrain_local, copy_nsnow_local, &
     qrout_grid_local, qsout_grid_local, nrout_grid_local, nsout_grid_local, qrout_grid_ptr_local, qsout_grid_ptr_local, &
     nrout_grid_ptr_local, nsout_grid_ptr_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr

  integer, intent(in) :: ncol_local
  logical, intent(in) :: copy_qrain_local, copy_qsnow_local, copy_nrain_local, copy_nsnow_local
  real(r8), target, intent(in) :: qrout_grid_local(pcols,pver), qsout_grid_local(pcols,pver)
  real(r8), target, intent(in) :: nrout_grid_local(pcols,pver), nsout_grid_local(pcols,pver)
  real(r8), pointer, intent(inout) :: qrout_grid_ptr_local(:,:), qsout_grid_ptr_local(:,:)
  real(r8), pointer, intent(inout) :: nrout_grid_ptr_local(:,:), nsout_grid_ptr_local(:,:)
  type(c_ptr) :: qrout_dst_p, qsout_dst_p, nrout_dst_p, nsout_dst_p

  interface
     subroutine micro_mg_cam_pbuf_copy_codon(ncol_c, pcols_c, pver_c, copy_qrain_c, copy_qsnow_c, copy_nrain_c, &
          copy_nsnow_c, qrout_grid_p, qsout_grid_p, nrout_grid_p, nsout_grid_p, qrout_grid_ptr_p, qsout_grid_ptr_p, &
          nrout_grid_ptr_p, nsout_grid_ptr_p) bind(c, name="micro_mg_cam_pbuf_copy_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, copy_qrain_c, copy_qsnow_c, copy_nrain_c, copy_nsnow_c
       type(c_ptr), value :: qrout_grid_p, qsout_grid_p, nrout_grid_p, nsout_grid_p
       type(c_ptr), value :: qrout_grid_ptr_p, qsout_grid_ptr_p, nrout_grid_ptr_p, nsout_grid_ptr_p
     end subroutine micro_mg_cam_pbuf_copy_codon
  end interface

  call micro_mg_cam_select_pbuf_copy_impl()

  if (use_native_pbuf_copy_impl) then
     if (copy_qrain_local) qrout_grid_ptr_local = qrout_grid_local
     if (copy_qsnow_local) qsout_grid_ptr_local = qsout_grid_local
     if (copy_nrain_local) nrout_grid_ptr_local = nrout_grid_local
     if (copy_nsnow_local) nsout_grid_ptr_local = nsout_grid_local
     return
  end if

  call micro_mg_cam_pbuf_copy_log_entered()

  qrout_dst_p = c_loc(qrout_grid_local)
  qsout_dst_p = c_loc(qsout_grid_local)
  nrout_dst_p = c_loc(nrout_grid_local)
  nsout_dst_p = c_loc(nsout_grid_local)
  if (copy_qrain_local) qrout_dst_p = c_loc(qrout_grid_ptr_local)
  if (copy_qsnow_local) qsout_dst_p = c_loc(qsout_grid_ptr_local)
  if (copy_nrain_local) nrout_dst_p = c_loc(nrout_grid_ptr_local)
  if (copy_nsnow_local) nsout_dst_p = c_loc(nsout_grid_ptr_local)

  call micro_mg_cam_pbuf_copy_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(merge(1, 0, copy_qrain_local), c_int64_t), int(merge(1, 0, copy_qsnow_local), c_int64_t), &
       int(merge(1, 0, copy_nrain_local), c_int64_t), int(merge(1, 0, copy_nsnow_local), c_int64_t), &
       c_loc(qrout_grid_local), c_loc(qsout_grid_local), c_loc(nrout_grid_local), c_loc(nsout_grid_local), &
       qrout_dst_p, qsout_dst_p, nrout_dst_p, nsout_dst_p)

end subroutine micro_mg_cam_pbuf_copy

subroutine micro_mg_cam_select_wtrc_shell_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (wtrc_shell_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MICRO_MG_CAM_WTRC_SHELL_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_wtrc_shell_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_wtrc_shell_impl = .false.
  end if

  wtrc_shell_impl_selected = .true.

  if (use_native_wtrc_shell_impl) then
     write(iulog,*) 'micro_mg_cam_wtrc_shell implementation = native'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_WTRC_SHELL_PROOF_FILE', &
          'micro_mg_cam_wtrc_shell implementation = native')
  else
     write(iulog,*) 'micro_mg_cam_wtrc_shell implementation = codon'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_WTRC_SHELL_PROOF_FILE', &
          'micro_mg_cam_wtrc_shell implementation = codon')
  end if
  call flush(iulog)

end subroutine micro_mg_cam_select_wtrc_shell_impl

subroutine micro_mg_cam_wtrc_shell_codon_wrap(ncol_local, preo_grid_local, prdso_grid_local, cmeiout_grid_local, meltso_grid_local, &
     qcsedten_grid_local, qisedten_grid_local, mnuccco_grid_local, mnuccto_grid_local, msacwio_grid_local, prao_grid_local, &
     prco_grid_local, psacwso_grid_local, bergo_grid_local, bergso_grid_local, praio_grid_local, prcio_grid_local, &
     pracso_grid_local, mnuccro_grid_local, qcreso_grid_local, qireso_grid_local, homoo_grid_local, melto_grid_local, &
     pre_rates_grid_local, sed_rates_grid_local, post_rates_grid_local, pcmei_grid_local, ncmei_grid_local, pmelts_grid_local, &
     nmelts_grid_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use water_types, only: pwtype, iwtvap, iwtliq, iwtice, iwtstrain, iwtstsnow
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: preo_grid_local(pcols,pver), prdso_grid_local(pcols,pver)
  real(r8), target, intent(in) :: cmeiout_grid_local(pcols,pver), meltso_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qcsedten_grid_local(pcols,pver), qisedten_grid_local(pcols,pver)
  real(r8), target, intent(in) :: mnuccco_grid_local(pcols,pver), mnuccto_grid_local(pcols,pver), msacwio_grid_local(pcols,pver)
  real(r8), target, intent(in) :: prao_grid_local(pcols,pver), prco_grid_local(pcols,pver), psacwso_grid_local(pcols,pver)
  real(r8), target, intent(in) :: bergo_grid_local(pcols,pver), bergso_grid_local(pcols,pver)
  real(r8), target, intent(in) :: praio_grid_local(pcols,pver), prcio_grid_local(pcols,pver)
  real(r8), target, intent(in) :: pracso_grid_local(pcols,pver), mnuccro_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qcreso_grid_local(pcols,pver), qireso_grid_local(pcols,pver), homoo_grid_local(pcols,pver)
  real(r8), target, intent(in) :: melto_grid_local(pcols,pver)
  real(r8), target, intent(out) :: pre_rates_grid_local(pcols,pver,pwtype,pwtype,pwtype)
  real(r8), target, intent(out) :: sed_rates_grid_local(pcols,pver,pwtype)
  real(r8), target, intent(out) :: post_rates_grid_local(pcols,pver,pwtype,pwtype,pwtype)
  real(r8), target, intent(inout) :: pcmei_grid_local(pcols,pver), ncmei_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: pmelts_grid_local(pcols,pver), nmelts_grid_local(pcols,pver)

  interface
     subroutine micro_mg_cam_wtrc_shell_codon(ncol_c, pcols_c, pver_c, top_lev_c, pwtype_c, iwtvap_c, iwtliq_c, iwtice_c, &
          iwtstrain_c, iwtstsnow_c, preo_grid_p, prdso_grid_p, cmeiout_grid_p, meltso_grid_p, qcsedten_grid_p, &
          qisedten_grid_p, mnuccco_grid_p, mnuccto_grid_p, msacwio_grid_p, prao_grid_p, prco_grid_p, psacwso_grid_p, &
          bergo_grid_p, bergso_grid_p, praio_grid_p, prcio_grid_p, pracso_grid_p, mnuccro_grid_p, qcreso_grid_p, &
          qireso_grid_p, homoo_grid_p, melto_grid_p, pre_rates_grid_p, sed_rates_grid_p, post_rates_grid_p, pcmei_grid_p, &
          ncmei_grid_p, pmelts_grid_p, nmelts_grid_p) bind(c, name="micro_mg_cam_wtrc_shell_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, pwtype_c, iwtvap_c, iwtliq_c, iwtice_c
       integer(c_int64_t), value :: iwtstrain_c, iwtstsnow_c
       type(c_ptr), value :: preo_grid_p, prdso_grid_p, cmeiout_grid_p, meltso_grid_p, qcsedten_grid_p, qisedten_grid_p
       type(c_ptr), value :: mnuccco_grid_p, mnuccto_grid_p, msacwio_grid_p, prao_grid_p, prco_grid_p, psacwso_grid_p
       type(c_ptr), value :: bergo_grid_p, bergso_grid_p, praio_grid_p, prcio_grid_p, pracso_grid_p, mnuccro_grid_p
       type(c_ptr), value :: qcreso_grid_p, qireso_grid_p, homoo_grid_p, melto_grid_p, pre_rates_grid_p, sed_rates_grid_p
       type(c_ptr), value :: post_rates_grid_p, pcmei_grid_p, ncmei_grid_p, pmelts_grid_p, nmelts_grid_p
     end subroutine micro_mg_cam_wtrc_shell_codon
  end interface

  call micro_mg_cam_wtrc_shell_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), int(pwtype, c_int64_t), int(iwtvap, c_int64_t), int(iwtliq, c_int64_t), int(iwtice, c_int64_t), &
       int(iwtstrain, c_int64_t), int(iwtstsnow, c_int64_t), c_loc(preo_grid_local), c_loc(prdso_grid_local), &
       c_loc(cmeiout_grid_local), c_loc(meltso_grid_local), c_loc(qcsedten_grid_local), c_loc(qisedten_grid_local), &
       c_loc(mnuccco_grid_local), c_loc(mnuccto_grid_local), c_loc(msacwio_grid_local), c_loc(prao_grid_local), &
       c_loc(prco_grid_local), c_loc(psacwso_grid_local), c_loc(bergo_grid_local), c_loc(bergso_grid_local), &
       c_loc(praio_grid_local), c_loc(prcio_grid_local), c_loc(pracso_grid_local), c_loc(mnuccro_grid_local), &
       c_loc(qcreso_grid_local), c_loc(qireso_grid_local), c_loc(homoo_grid_local), c_loc(melto_grid_local), &
       c_loc(pre_rates_grid_local), c_loc(sed_rates_grid_local), c_loc(post_rates_grid_local), c_loc(pcmei_grid_local), &
       c_loc(ncmei_grid_local), c_loc(pmelts_grid_local), c_loc(nmelts_grid_local))

end subroutine micro_mg_cam_wtrc_shell_codon_wrap

subroutine micro_mg_cam_select_wtrc_prep_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (wtrc_prep_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MICRO_MG_CAM_WTRC_PREP_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_wtrc_prep_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_wtrc_prep_impl = .false.
  end if

  wtrc_prep_impl_selected = .true.

  if (use_native_wtrc_prep_impl) then
     write(iulog,*) 'micro_mg_cam_wtrc_prep implementation = native'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_WTRC_PREP_PROOF_FILE', &
          'micro_mg_cam_wtrc_prep implementation = native')
  else
     write(iulog,*) 'micro_mg_cam_wtrc_prep implementation = codon'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_WTRC_PREP_PROOF_FILE', &
          'micro_mg_cam_wtrc_prep implementation = codon')
  end if
  call flush(iulog)

end subroutine micro_mg_cam_select_wtrc_prep_impl

subroutine micro_mg_cam_wtrc_prep(ncol_local, cmeiout_grid_local, meltso_grid_local, qcsedten_grid_local, qisedten_grid_local, &
     pcmei_grid_local, ncmei_grid_local, pmelts_grid_local, nmelts_grid_local, sed_rates_grid_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use water_types, only: pwtype, iwtliq, iwtice
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), target, intent(in) :: cmeiout_grid_local(pcols,pver), meltso_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qcsedten_grid_local(pcols,pver), qisedten_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: pcmei_grid_local(pcols,pver), ncmei_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: pmelts_grid_local(pcols,pver), nmelts_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: sed_rates_grid_local(pcols,pver,pwtype)

  interface
     subroutine micro_mg_cam_wtrc_prep_codon(ncol_c, pcols_c, pver_c, top_lev_c, pwtype_c, iwtliq_c, iwtice_c, &
          cmeiout_grid_p, meltso_grid_p, qcsedten_grid_p, qisedten_grid_p, pcmei_grid_p, ncmei_grid_p, pmelts_grid_p, &
          nmelts_grid_p, sed_rates_grid_p) bind(c, name="micro_mg_cam_wtrc_prep_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, pwtype_c, iwtliq_c, iwtice_c
       type(c_ptr), value :: cmeiout_grid_p, meltso_grid_p, qcsedten_grid_p, qisedten_grid_p
       type(c_ptr), value :: pcmei_grid_p, ncmei_grid_p, pmelts_grid_p, nmelts_grid_p, sed_rates_grid_p
     end subroutine micro_mg_cam_wtrc_prep_codon
  end interface

  call micro_mg_cam_select_wtrc_prep_impl()

  if (use_native_wtrc_prep_impl) then
     call micro_mg_cam_wtrc_prep_native(ncol_local, cmeiout_grid_local, meltso_grid_local, qcsedten_grid_local, &
          qisedten_grid_local, pcmei_grid_local, ncmei_grid_local, pmelts_grid_local, nmelts_grid_local, sed_rates_grid_local)
     return
  end if

  call micro_mg_cam_wtrc_prep_codon(int(ncol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), int(pwtype, c_int64_t), int(iwtliq, c_int64_t), int(iwtice, c_int64_t), &
       c_loc(cmeiout_grid_local), c_loc(meltso_grid_local), c_loc(qcsedten_grid_local), c_loc(qisedten_grid_local), &
       c_loc(pcmei_grid_local), c_loc(ncmei_grid_local), c_loc(pmelts_grid_local), c_loc(nmelts_grid_local), &
       c_loc(sed_rates_grid_local))

end subroutine micro_mg_cam_wtrc_prep

subroutine micro_mg_cam_wtrc_prep_native(ncol_local, cmeiout_grid_local, meltso_grid_local, qcsedten_grid_local, &
     qisedten_grid_local, pcmei_grid_local, ncmei_grid_local, pmelts_grid_local, nmelts_grid_local, sed_rates_grid_local)

  use water_types, only: pwtype, iwtliq, iwtice
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local
  real(r8), intent(in) :: cmeiout_grid_local(pcols,pver), meltso_grid_local(pcols,pver)
  real(r8), intent(in) :: qcsedten_grid_local(pcols,pver), qisedten_grid_local(pcols,pver)
  real(r8), intent(inout) :: pcmei_grid_local(pcols,pver), ncmei_grid_local(pcols,pver)
  real(r8), intent(inout) :: pmelts_grid_local(pcols,pver), nmelts_grid_local(pcols,pver)
  real(r8), intent(inout) :: sed_rates_grid_local(pcols,pver,pwtype)
  integer :: i, k, m

  do m = 1, pwtype
     do k = top_lev, pver
        do i = 1, pcols
           sed_rates_grid_local(i,k,m) = 0._r8
        end do
     end do
  end do

  do k = top_lev, pver
     do i = 1, pcols
        pcmei_grid_local(i,k) = 0._r8
        ncmei_grid_local(i,k) = 0._r8
        pmelts_grid_local(i,k) = 0._r8
        nmelts_grid_local(i,k) = 0._r8
     end do
  end do

  do k = top_lev, pver
     do i = 1, ncol_local
        if (cmeiout_grid_local(i,k) < 0._r8) then
           ncmei_grid_local(i,k) = cmeiout_grid_local(i,k)
        else
           pcmei_grid_local(i,k) = cmeiout_grid_local(i,k)
        end if
        if (meltso_grid_local(i,k) < 0._r8) then
           nmelts_grid_local(i,k) = meltso_grid_local(i,k)
        else
           pmelts_grid_local(i,k) = meltso_grid_local(i,k)
        end if
        sed_rates_grid_local(i,k,iwtliq) = qcsedten_grid_local(i,k)
        sed_rates_grid_local(i,k,iwtice) = qisedten_grid_local(i,k)
     end do
  end do

end subroutine micro_mg_cam_wtrc_prep_native

subroutine micro_mg_cam_select_budget_diag_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (budget_diag_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MICRO_MG_CAM_BUDGET_DIAG_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_budget_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_budget_diag_impl = .false.
  end if

  budget_diag_impl_selected = .true.

  if (use_native_budget_diag_impl) then
     write(iulog,*) 'micro_mg_cam_budget_diag implementation = native'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_BUDGET_DIAG_PROOF_FILE', &
          'micro_mg_cam_budget_diag implementation = native')
  else
     write(iulog,*) 'micro_mg_cam_budget_diag implementation = codon'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_BUDGET_DIAG_PROOF_FILE', &
          'micro_mg_cam_budget_diag implementation = codon')
  end if
  call flush(iulog)

end subroutine micro_mg_cam_select_budget_diag_impl

subroutine micro_mg_cam_budget_diag(mode_local, ncol_local, qcreso_grid_local, melto_grid_local, mnuccco_grid_local, &
     mnuccto_grid_local, bergo_grid_local, homoo_grid_local, msacwio_grid_local, prao_grid_local, prco_grid_local, &
     psacwso_grid_local, bergso_grid_local, cmeiout_grid_local, qireso_grid_local, prcio_grid_local, praio_grid_local, &
     ftem_grid_local)

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: mode_local, ncol_local
  real(r8), target, intent(in) :: qcreso_grid_local(pcols,pver), melto_grid_local(pcols,pver)
  real(r8), target, intent(in) :: mnuccco_grid_local(pcols,pver), mnuccto_grid_local(pcols,pver)
  real(r8), target, intent(in) :: bergo_grid_local(pcols,pver), homoo_grid_local(pcols,pver)
  real(r8), target, intent(in) :: msacwio_grid_local(pcols,pver), prao_grid_local(pcols,pver)
  real(r8), target, intent(in) :: prco_grid_local(pcols,pver), psacwso_grid_local(pcols,pver)
  real(r8), target, intent(in) :: bergso_grid_local(pcols,pver), cmeiout_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qireso_grid_local(pcols,pver), prcio_grid_local(pcols,pver), praio_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: ftem_grid_local(pcols,pver)

  interface
     subroutine micro_mg_cam_budget_diag_codon(mode_c, ncol_c, pcols_c, pver_c, top_lev_c, qcreso_grid_p, melto_grid_p, &
          mnuccco_grid_p, mnuccto_grid_p, bergo_grid_p, homoo_grid_p, msacwio_grid_p, prao_grid_p, prco_grid_p, &
          psacwso_grid_p, bergso_grid_p, cmeiout_grid_p, qireso_grid_p, prcio_grid_p, praio_grid_p, ftem_grid_p) &
          bind(c, name="micro_mg_cam_budget_diag_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: qcreso_grid_p, melto_grid_p, mnuccco_grid_p, mnuccto_grid_p, bergo_grid_p, homoo_grid_p
       type(c_ptr), value :: msacwio_grid_p, prao_grid_p, prco_grid_p, psacwso_grid_p, bergso_grid_p, cmeiout_grid_p
       type(c_ptr), value :: qireso_grid_p, prcio_grid_p, praio_grid_p, ftem_grid_p
     end subroutine micro_mg_cam_budget_diag_codon
  end interface

  call micro_mg_cam_select_budget_diag_impl()

  if (use_native_budget_diag_impl) then
     call micro_mg_cam_budget_diag_native(mode_local, ncol_local, qcreso_grid_local, melto_grid_local, mnuccco_grid_local, &
          mnuccto_grid_local, bergo_grid_local, homoo_grid_local, msacwio_grid_local, prao_grid_local, prco_grid_local, &
          psacwso_grid_local, bergso_grid_local, cmeiout_grid_local, qireso_grid_local, prcio_grid_local, praio_grid_local, &
          ftem_grid_local)
     return
  end if

  call micro_mg_cam_budget_diag_codon(int(mode_local, c_int64_t), int(ncol_local, c_int64_t), int(pcols, c_int64_t), &
       int(pver, c_int64_t), int(top_lev, c_int64_t), c_loc(qcreso_grid_local), c_loc(melto_grid_local), &
       c_loc(mnuccco_grid_local), c_loc(mnuccto_grid_local), c_loc(bergo_grid_local), c_loc(homoo_grid_local), &
       c_loc(msacwio_grid_local), c_loc(prao_grid_local), c_loc(prco_grid_local), c_loc(psacwso_grid_local), &
       c_loc(bergso_grid_local), c_loc(cmeiout_grid_local), c_loc(qireso_grid_local), c_loc(prcio_grid_local), &
       c_loc(praio_grid_local), c_loc(ftem_grid_local))

end subroutine micro_mg_cam_budget_diag

subroutine micro_mg_cam_budget_diag_native(mode_local, ncol_local, qcreso_grid_local, melto_grid_local, mnuccco_grid_local, &
     mnuccto_grid_local, bergo_grid_local, homoo_grid_local, msacwio_grid_local, prao_grid_local, prco_grid_local, &
     psacwso_grid_local, bergso_grid_local, cmeiout_grid_local, qireso_grid_local, prcio_grid_local, praio_grid_local, &
     ftem_grid_local)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: mode_local, ncol_local
  real(r8), intent(in) :: qcreso_grid_local(pcols,pver), melto_grid_local(pcols,pver)
  real(r8), intent(in) :: mnuccco_grid_local(pcols,pver), mnuccto_grid_local(pcols,pver)
  real(r8), intent(in) :: bergo_grid_local(pcols,pver), homoo_grid_local(pcols,pver)
  real(r8), intent(in) :: msacwio_grid_local(pcols,pver), prao_grid_local(pcols,pver)
  real(r8), intent(in) :: prco_grid_local(pcols,pver), psacwso_grid_local(pcols,pver)
  real(r8), intent(in) :: bergso_grid_local(pcols,pver), cmeiout_grid_local(pcols,pver)
  real(r8), intent(in) :: qireso_grid_local(pcols,pver), prcio_grid_local(pcols,pver), praio_grid_local(pcols,pver)
  real(r8), intent(inout) :: ftem_grid_local(pcols,pver)

  ftem_grid_local = 0._r8

  select case (mode_local)
  case (1)
     ftem_grid_local(:ncol_local,top_lev:pver) = qcreso_grid_local(:ncol_local,top_lev:pver)
  case (2)
     ftem_grid_local(:ncol_local,top_lev:pver) = melto_grid_local(:ncol_local,top_lev:pver) - &
          mnuccco_grid_local(:ncol_local,top_lev:pver) - mnuccto_grid_local(:ncol_local,top_lev:pver) - &
          bergo_grid_local(:ncol_local,top_lev:pver) - homoo_grid_local(:ncol_local,top_lev:pver) - &
          msacwio_grid_local(:ncol_local,top_lev:pver)
  case (3)
     ftem_grid_local(:ncol_local,top_lev:pver) = -prao_grid_local(:ncol_local,top_lev:pver) - &
          prco_grid_local(:ncol_local,top_lev:pver) - psacwso_grid_local(:ncol_local,top_lev:pver) - &
          bergso_grid_local(:ncol_local,top_lev:pver)
  case (4)
     ftem_grid_local(:ncol_local,top_lev:pver) = cmeiout_grid_local(:ncol_local,top_lev:pver) + &
          qireso_grid_local(:ncol_local,top_lev:pver)
  case (5)
     ftem_grid_local(:ncol_local,top_lev:pver) = -melto_grid_local(:ncol_local,top_lev:pver) + &
          mnuccco_grid_local(:ncol_local,top_lev:pver) + mnuccto_grid_local(:ncol_local,top_lev:pver) + &
          bergo_grid_local(:ncol_local,top_lev:pver) + homoo_grid_local(:ncol_local,top_lev:pver) + &
          msacwio_grid_local(:ncol_local,top_lev:pver)
  case (6)
     ftem_grid_local(:ncol_local,top_lev:pver) = -prcio_grid_local(:ncol_local,top_lev:pver) - &
          praio_grid_local(:ncol_local,top_lev:pver)
  case default
     call endrun('micro_mg_cam_budget_diag: invalid mode')
  end select

end subroutine micro_mg_cam_budget_diag_native

subroutine micro_mg_cam_append_impl_proof(env_name, proof_line)

  character(len=*), intent(in) :: env_name, proof_line
  character(len=512) :: proof_path
  integer :: status, n, unit_id

  proof_path = ''
  call get_environment_variable(env_name, value=proof_path, length=n, status=status)
  if (status /= 0 .or. n <= 0) return

  open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
       position='append', iostat=status)
  if (status /= 0) return

  write(unit_id,'(A)') trim(proof_line)
  close(unit_id)

end subroutine micro_mg_cam_append_impl_proof

subroutine micro_mg_cam_select_reff_calc_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (reff_calc_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MICRO_MG_CAM_REFF_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_reff_calc_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_reff_calc_impl = .false.
  end if

  reff_calc_impl_selected = .true.

  if (use_native_reff_calc_impl) then
     write(iulog,*) 'micro_mg_cam_reff_calc implementation = native'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_REFF_PROOF_FILE', &
          'micro_mg_cam_reff_calc implementation = native')
  else
     write(iulog,*) 'micro_mg_cam_reff_calc implementation = codon'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_REFF_PROOF_FILE', &
          'micro_mg_cam_reff_calc implementation = codon')
  end if
  call flush(iulog)

end subroutine micro_mg_cam_select_reff_calc_impl

subroutine micro_mg_cam_select_reff_calc_compare()

  integer :: status, n

  if (reff_calc_compare_selected) return

  call get_environment_variable('MICRO_MG_CAM_REFF_COMPARE_FILE', length=n, status=status)
  use_reff_calc_compare = status == 0 .and. n > 0
  reff_calc_compare_selected = .true.

end subroutine micro_mg_cam_select_reff_calc_compare

subroutine micro_mg_cam_reff_calc(ngrdcol_local, micro_mg_version_local, rho_grid_local, icwmrst_grid_local, &
     liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, qrout_grid_local, &
     nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, icimrst_grid_local, &
     ast_grid_local, mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local, &
     drout2_grid_local, reff_rain_grid_local, des_grid_local, dsout2_grid_local, reff_snow_grid_local, rei_grid_local, &
     niic_grid_local, dei_grid_local, mgreffrain_grid_local, mgreffsnow_grid_local)

  integer, intent(in) :: ngrdcol_local, micro_mg_version_local
  real(r8), intent(in) :: rho_grid_local(pcols,pver), icwmrst_grid_local(pcols,pver)
  real(r8), intent(in) :: liqcldf_grid_local(pcols,pver), nc_grid_local(pcols,pver)
  real(r8), intent(in) :: qr_grid_local(pcols,pver), nr_grid_local(pcols,pver)
  real(r8), intent(in) :: qs_grid_local(pcols,pver), ns_grid_local(pcols,pver)
  real(r8), intent(in) :: qrout_grid_local(pcols,pver), nrout_grid_local(pcols,pver)
  real(r8), intent(in) :: qsout_grid_local(pcols,pver), nsout_grid_local(pcols,pver)
  real(r8), intent(in) :: ni_grid_local(pcols,pver), icecldf_grid_local(pcols,pver)
  real(r8), intent(in) :: icimrst_grid_local(pcols,pver), ast_grid_local(pcols,pver)
  real(r8), intent(inout) :: mu_grid_local(pcols,pver), lambdac_grid_local(pcols,pver)
  real(r8), intent(inout) :: rel_fn_grid_local(pcols,pver), ncic_grid_local(pcols,pver)
  real(r8), intent(inout) :: rel_grid_local(pcols,pver), drout2_grid_local(pcols,pver)
  real(r8), intent(inout) :: reff_rain_grid_local(pcols,pver), des_grid_local(pcols,pver)
  real(r8), intent(inout) :: dsout2_grid_local(pcols,pver), reff_snow_grid_local(pcols,pver)
  real(r8), intent(inout) :: rei_grid_local(pcols,pver), niic_grid_local(pcols,pver)
  real(r8), intent(inout) :: dei_grid_local(pcols,pver), mgreffrain_grid_local(pcols,pver)
  real(r8), intent(inout) :: mgreffsnow_grid_local(pcols,pver)

  call micro_mg_cam_select_reff_calc_impl()
  call micro_mg_cam_select_reff_calc_compare()

  if (use_reff_calc_compare .and. .not. reff_calc_compare_done) then
     call micro_mg_cam_reff_calc_compare(ngrdcol_local, micro_mg_version_local, rho_grid_local, icwmrst_grid_local, &
          liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, qrout_grid_local, &
          nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, icimrst_grid_local, &
          ast_grid_local, mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local, &
          drout2_grid_local, reff_rain_grid_local, des_grid_local, dsout2_grid_local, reff_snow_grid_local, rei_grid_local, &
          niic_grid_local, dei_grid_local, mgreffrain_grid_local, mgreffsnow_grid_local)
     reff_calc_compare_done = .true.
     return
  end if

  if (use_native_reff_calc_impl) then
     call micro_mg_cam_reff_calc_native(ngrdcol_local, micro_mg_version_local, rho_grid_local, icwmrst_grid_local, &
          liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, qrout_grid_local, &
          nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, icimrst_grid_local, &
          ast_grid_local, mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local, &
          drout2_grid_local, reff_rain_grid_local, des_grid_local, dsout2_grid_local, reff_snow_grid_local, rei_grid_local, &
          niic_grid_local, dei_grid_local, mgreffrain_grid_local, mgreffsnow_grid_local)
     return
  end if

  call micro_mg_cam_reff_calc_codon_invoke(ngrdcol_local, micro_mg_version_local, rho_grid_local, icwmrst_grid_local, &
       liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, qrout_grid_local, &
       nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, icimrst_grid_local, &
       ast_grid_local, mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local, &
       drout2_grid_local, reff_rain_grid_local, des_grid_local, dsout2_grid_local, reff_snow_grid_local, rei_grid_local, &
       niic_grid_local, dei_grid_local, mgreffrain_grid_local, mgreffsnow_grid_local)

end subroutine micro_mg_cam_reff_calc

subroutine micro_mg_cam_reff_calc_codon_invoke(ngrdcol_local, micro_mg_version_local, rho_grid_local, icwmrst_grid_local, &
     liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, qrout_grid_local, &
     nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, icimrst_grid_local, &
     ast_grid_local, mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local, &
     drout2_grid_local, reff_rain_grid_local, des_grid_local, dsout2_grid_local, reff_snow_grid_local, rei_grid_local, &
     niic_grid_local, dei_grid_local, mgreffrain_grid_local, mgreffsnow_grid_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use micro_mg_utils, only: mg_liq_props, mg_ice_props, qsmall, mincld, rhosn, rhoi, rhow, rhows
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ngrdcol_local, micro_mg_version_local
  real(r8), parameter :: dcon_local = 25.e-6_r8
  real(r8), parameter :: mucon_local = 5.3_r8
  real(r8), parameter :: deicon_local = 50._r8
  real(r8), target, intent(in) :: rho_grid_local(pcols,pver), icwmrst_grid_local(pcols,pver)
  real(r8), target, intent(in) :: liqcldf_grid_local(pcols,pver), nc_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qr_grid_local(pcols,pver), nr_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qs_grid_local(pcols,pver), ns_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qrout_grid_local(pcols,pver), nrout_grid_local(pcols,pver)
  real(r8), target, intent(in) :: qsout_grid_local(pcols,pver), nsout_grid_local(pcols,pver)
  real(r8), target, intent(in) :: ni_grid_local(pcols,pver), icecldf_grid_local(pcols,pver)
  real(r8), target, intent(in) :: icimrst_grid_local(pcols,pver), ast_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: mu_grid_local(pcols,pver), lambdac_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: rel_fn_grid_local(pcols,pver), ncic_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: rel_grid_local(pcols,pver), drout2_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: reff_rain_grid_local(pcols,pver), des_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: dsout2_grid_local(pcols,pver), reff_snow_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: rei_grid_local(pcols,pver), niic_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: dei_grid_local(pcols,pver), mgreffrain_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: mgreffsnow_grid_local(pcols,pver)

  interface
     subroutine micro_mg_cam_reff_calc_codon(ngrdcol_c, pcols_c, pver_c, top_lev_c, micro_mg_version_c, qsmall_c, &
          mincld_c, liq_rho_c, liq_eff_dim_c, liq_min_mean_mass_c, ice_eff_dim_c, ice_shape_coef_c, ice_lambda_lo_c, &
          ice_lambda_hi_c, ice_min_mean_mass_c, rhosn_c, rhoi_c, rhow_c, rhows_c, mucon_c, dcon_c, deicon_c, rho_grid_p, &
          icwmrst_grid_p, liqcldf_grid_p, nc_grid_p, qr_grid_p, nr_grid_p, qs_grid_p, ns_grid_p, qrout_grid_p, &
          nrout_grid_p, qsout_grid_p, nsout_grid_p, ni_grid_p, icecldf_grid_p, icimrst_grid_p, ast_grid_p, mu_grid_p, &
          lambdac_grid_p, rel_fn_grid_p, ncic_grid_p, rel_grid_p, drout2_grid_p, reff_rain_grid_p, des_grid_p, &
          dsout2_grid_p, reff_snow_grid_p, rei_grid_p, niic_grid_p, dei_grid_p, mgreffrain_grid_p, &
          mgreffsnow_grid_p) bind(c, name="micro_mg_cam_reff_calc_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ngrdcol_c, pcols_c, pver_c, top_lev_c, micro_mg_version_c
       real(c_double), value :: qsmall_c, mincld_c, liq_rho_c, liq_eff_dim_c, liq_min_mean_mass_c
       real(c_double), value :: ice_eff_dim_c, ice_shape_coef_c, ice_lambda_lo_c, ice_lambda_hi_c
       real(c_double), value :: ice_min_mean_mass_c, rhosn_c, rhoi_c, rhow_c, rhows_c
       real(c_double), value :: mucon_c, dcon_c, deicon_c
       type(c_ptr), value :: rho_grid_p, icwmrst_grid_p, liqcldf_grid_p, nc_grid_p, qr_grid_p, nr_grid_p
       type(c_ptr), value :: qs_grid_p, ns_grid_p, qrout_grid_p, nrout_grid_p, qsout_grid_p, nsout_grid_p
       type(c_ptr), value :: ni_grid_p, icecldf_grid_p, icimrst_grid_p, ast_grid_p, mu_grid_p, lambdac_grid_p
       type(c_ptr), value :: rel_fn_grid_p, ncic_grid_p, rel_grid_p, drout2_grid_p, reff_rain_grid_p, des_grid_p
       type(c_ptr), value :: dsout2_grid_p, reff_snow_grid_p, rei_grid_p, niic_grid_p, dei_grid_p
       type(c_ptr), value :: mgreffrain_grid_p, mgreffsnow_grid_p
     end subroutine micro_mg_cam_reff_calc_codon
  end interface

  call micro_mg_cam_reff_liq_native(ngrdcol_local, rho_grid_local, icwmrst_grid_local, liqcldf_grid_local, nc_grid_local, &
       mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local)

  call micro_mg_cam_reff_calc_codon(int(ngrdcol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), int(micro_mg_version_local, c_int64_t), qsmall, mincld, mg_liq_props%rho, &
       mg_liq_props%eff_dim, mg_liq_props%min_mean_mass, mg_ice_props%eff_dim, mg_ice_props%shape_coef, &
       mg_ice_props%lambda_bounds(1), mg_ice_props%lambda_bounds(2), mg_ice_props%min_mean_mass, rhosn, rhoi, rhow, rhows, &
       mucon_local, dcon_local, deicon_local, c_loc(rho_grid_local), c_loc(icwmrst_grid_local), c_loc(liqcldf_grid_local), &
       c_loc(nc_grid_local), c_loc(qr_grid_local), c_loc(nr_grid_local), c_loc(qs_grid_local), c_loc(ns_grid_local), &
       c_loc(qrout_grid_local), c_loc(nrout_grid_local), c_loc(qsout_grid_local), c_loc(nsout_grid_local), c_loc(ni_grid_local), &
       c_loc(icecldf_grid_local), c_loc(icimrst_grid_local), c_loc(ast_grid_local), c_loc(mu_grid_local), c_loc(lambdac_grid_local), &
       c_loc(rel_fn_grid_local), c_loc(ncic_grid_local), c_loc(rel_grid_local), c_loc(drout2_grid_local), c_loc(reff_rain_grid_local), &
       c_loc(des_grid_local), c_loc(dsout2_grid_local), c_loc(reff_snow_grid_local), c_loc(rei_grid_local), c_loc(niic_grid_local), &
       c_loc(dei_grid_local), c_loc(mgreffrain_grid_local), c_loc(mgreffsnow_grid_local))

end subroutine micro_mg_cam_reff_calc_codon_invoke

subroutine micro_mg_cam_reff_calc_compare(ngrdcol_local, micro_mg_version_local, rho_grid_local, icwmrst_grid_local, &
     liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, qrout_grid_local, &
     nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, icimrst_grid_local, &
     ast_grid_local, mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local, &
     drout2_grid_local, reff_rain_grid_local, des_grid_local, dsout2_grid_local, reff_snow_grid_local, rei_grid_local, &
     niic_grid_local, dei_grid_local, mgreffrain_grid_local, mgreffsnow_grid_local)

  integer, intent(in) :: ngrdcol_local, micro_mg_version_local
  real(r8), intent(in) :: rho_grid_local(pcols,pver), icwmrst_grid_local(pcols,pver)
  real(r8), intent(in) :: liqcldf_grid_local(pcols,pver), nc_grid_local(pcols,pver)
  real(r8), intent(in) :: qr_grid_local(pcols,pver), nr_grid_local(pcols,pver)
  real(r8), intent(in) :: qs_grid_local(pcols,pver), ns_grid_local(pcols,pver)
  real(r8), intent(in) :: qrout_grid_local(pcols,pver), nrout_grid_local(pcols,pver)
  real(r8), intent(in) :: qsout_grid_local(pcols,pver), nsout_grid_local(pcols,pver)
  real(r8), intent(in) :: ni_grid_local(pcols,pver), icecldf_grid_local(pcols,pver)
  real(r8), intent(in) :: icimrst_grid_local(pcols,pver), ast_grid_local(pcols,pver)
  real(r8), intent(inout) :: mu_grid_local(pcols,pver), lambdac_grid_local(pcols,pver)
  real(r8), intent(inout) :: rel_fn_grid_local(pcols,pver), ncic_grid_local(pcols,pver)
  real(r8), intent(inout) :: rel_grid_local(pcols,pver), drout2_grid_local(pcols,pver)
  real(r8), intent(inout) :: reff_rain_grid_local(pcols,pver), des_grid_local(pcols,pver)
  real(r8), intent(inout) :: dsout2_grid_local(pcols,pver), reff_snow_grid_local(pcols,pver)
  real(r8), intent(inout) :: rei_grid_local(pcols,pver), niic_grid_local(pcols,pver)
  real(r8), intent(inout) :: dei_grid_local(pcols,pver), mgreffrain_grid_local(pcols,pver)
  real(r8), intent(inout) :: mgreffsnow_grid_local(pcols,pver)
  real(r8), allocatable, target :: native_mu_grid(:,:), codon_mu_grid(:,:)
  real(r8), allocatable, target :: native_lambdac_grid(:,:), codon_lambdac_grid(:,:)
  real(r8), allocatable, target :: native_rel_fn_grid(:,:), codon_rel_fn_grid(:,:)
  real(r8), allocatable, target :: native_ncic_grid(:,:), codon_ncic_grid(:,:)
  real(r8), allocatable, target :: native_rel_grid(:,:), codon_rel_grid(:,:)
  real(r8), allocatable, target :: native_drout2_grid(:,:), codon_drout2_grid(:,:)
  real(r8), allocatable, target :: native_reff_rain_grid(:,:), codon_reff_rain_grid(:,:)
  real(r8), allocatable, target :: native_des_grid(:,:), codon_des_grid(:,:)
  real(r8), allocatable, target :: native_dsout2_grid(:,:), codon_dsout2_grid(:,:)
  real(r8), allocatable, target :: native_reff_snow_grid(:,:), codon_reff_snow_grid(:,:)
  real(r8), allocatable, target :: native_rei_grid(:,:), codon_rei_grid(:,:)
  real(r8), allocatable, target :: native_niic_grid(:,:), codon_niic_grid(:,:)
  real(r8), allocatable, target :: native_dei_grid(:,:), codon_dei_grid(:,:)
  real(r8), allocatable, target :: native_mgreffrain_grid(:,:), codon_mgreffrain_grid(:,:)
  real(r8), allocatable, target :: native_mgreffsnow_grid(:,:), codon_mgreffsnow_grid(:,:)

  allocate(native_mu_grid(pcols,pver), codon_mu_grid(pcols,pver))
  allocate(native_lambdac_grid(pcols,pver), codon_lambdac_grid(pcols,pver))
  allocate(native_rel_fn_grid(pcols,pver), codon_rel_fn_grid(pcols,pver))
  allocate(native_ncic_grid(pcols,pver), codon_ncic_grid(pcols,pver))
  allocate(native_rel_grid(pcols,pver), codon_rel_grid(pcols,pver))
  allocate(native_drout2_grid(pcols,pver), codon_drout2_grid(pcols,pver))
  allocate(native_reff_rain_grid(pcols,pver), codon_reff_rain_grid(pcols,pver))
  allocate(native_des_grid(pcols,pver), codon_des_grid(pcols,pver))
  allocate(native_dsout2_grid(pcols,pver), codon_dsout2_grid(pcols,pver))
  allocate(native_reff_snow_grid(pcols,pver), codon_reff_snow_grid(pcols,pver))
  allocate(native_rei_grid(pcols,pver), codon_rei_grid(pcols,pver))
  allocate(native_niic_grid(pcols,pver), codon_niic_grid(pcols,pver))
  allocate(native_dei_grid(pcols,pver), codon_dei_grid(pcols,pver))
  allocate(native_mgreffrain_grid(pcols,pver), codon_mgreffrain_grid(pcols,pver))
  allocate(native_mgreffsnow_grid(pcols,pver), codon_mgreffsnow_grid(pcols,pver))

  native_mu_grid = mu_grid_local
  codon_mu_grid = mu_grid_local
  native_lambdac_grid = lambdac_grid_local
  codon_lambdac_grid = lambdac_grid_local
  native_rel_fn_grid = rel_fn_grid_local
  codon_rel_fn_grid = rel_fn_grid_local
  native_ncic_grid = ncic_grid_local
  codon_ncic_grid = ncic_grid_local
  native_rel_grid = rel_grid_local
  codon_rel_grid = rel_grid_local
  native_drout2_grid = drout2_grid_local
  codon_drout2_grid = drout2_grid_local
  native_reff_rain_grid = reff_rain_grid_local
  codon_reff_rain_grid = reff_rain_grid_local
  native_des_grid = des_grid_local
  codon_des_grid = des_grid_local
  native_dsout2_grid = dsout2_grid_local
  codon_dsout2_grid = dsout2_grid_local
  native_reff_snow_grid = reff_snow_grid_local
  codon_reff_snow_grid = reff_snow_grid_local
  native_rei_grid = rei_grid_local
  codon_rei_grid = rei_grid_local
  native_niic_grid = niic_grid_local
  codon_niic_grid = niic_grid_local
  native_dei_grid = dei_grid_local
  codon_dei_grid = dei_grid_local
  native_mgreffrain_grid = mgreffrain_grid_local
  codon_mgreffrain_grid = mgreffrain_grid_local
  native_mgreffsnow_grid = mgreffsnow_grid_local
  codon_mgreffsnow_grid = mgreffsnow_grid_local

  call micro_mg_cam_reff_calc_native(ngrdcol_local, micro_mg_version_local, rho_grid_local, icwmrst_grid_local, &
       liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, qrout_grid_local, &
       nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, icimrst_grid_local, &
       ast_grid_local, native_mu_grid, native_lambdac_grid, native_rel_fn_grid, native_ncic_grid, native_rel_grid, &
       native_drout2_grid, native_reff_rain_grid, native_des_grid, native_dsout2_grid, native_reff_snow_grid, native_rei_grid, &
       native_niic_grid, native_dei_grid, native_mgreffrain_grid, native_mgreffsnow_grid)

  call micro_mg_cam_reff_calc_codon_invoke(ngrdcol_local, micro_mg_version_local, rho_grid_local, icwmrst_grid_local, &
       liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, qrout_grid_local, &
       nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, icimrst_grid_local, &
       ast_grid_local, codon_mu_grid, codon_lambdac_grid, codon_rel_fn_grid, codon_ncic_grid, codon_rel_grid, &
       codon_drout2_grid, codon_reff_rain_grid, codon_des_grid, codon_dsout2_grid, codon_reff_snow_grid, codon_rei_grid, &
       codon_niic_grid, codon_dei_grid, codon_mgreffrain_grid, codon_mgreffsnow_grid)

  call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_REFF_COMPARE_FILE', &
       'micro_mg_cam_reff_calc direct compare begin')
  call report_field('mu_grid', native_mu_grid, codon_mu_grid)
  call report_field('lambdac_grid', native_lambdac_grid, codon_lambdac_grid)
  call report_field('rel_fn_grid', native_rel_fn_grid, codon_rel_fn_grid)
  call report_field('ncic_grid', native_ncic_grid, codon_ncic_grid)
  call report_field('rel_grid', native_rel_grid, codon_rel_grid)
  call report_field('drout2_grid', native_drout2_grid, codon_drout2_grid)
  call report_field('reff_rain_grid', native_reff_rain_grid, codon_reff_rain_grid)
  call report_field('des_grid', native_des_grid, codon_des_grid)
  call report_field('dsout2_grid', native_dsout2_grid, codon_dsout2_grid)
  call report_field('reff_snow_grid', native_reff_snow_grid, codon_reff_snow_grid)
  call report_field('rei_grid', native_rei_grid, codon_rei_grid)
  call report_field('niic_grid', native_niic_grid, codon_niic_grid)
  call report_field('dei_grid', native_dei_grid, codon_dei_grid)
  call report_field('mgreffrain_grid', native_mgreffrain_grid, codon_mgreffrain_grid)
  call report_field('mgreffsnow_grid', native_mgreffsnow_grid, codon_mgreffsnow_grid)
  call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_REFF_COMPARE_FILE', &
       'micro_mg_cam_reff_calc direct compare end')

  if (use_native_reff_calc_impl) then
     mu_grid_local = native_mu_grid
     lambdac_grid_local = native_lambdac_grid
     rel_fn_grid_local = native_rel_fn_grid
     ncic_grid_local = native_ncic_grid
     rel_grid_local = native_rel_grid
     drout2_grid_local = native_drout2_grid
     reff_rain_grid_local = native_reff_rain_grid
     des_grid_local = native_des_grid
     dsout2_grid_local = native_dsout2_grid
     reff_snow_grid_local = native_reff_snow_grid
     rei_grid_local = native_rei_grid
     niic_grid_local = native_niic_grid
     dei_grid_local = native_dei_grid
     mgreffrain_grid_local = native_mgreffrain_grid
     mgreffsnow_grid_local = native_mgreffsnow_grid
  else
     mu_grid_local = codon_mu_grid
     lambdac_grid_local = codon_lambdac_grid
     rel_fn_grid_local = codon_rel_fn_grid
     ncic_grid_local = codon_ncic_grid
     rel_grid_local = codon_rel_grid
     drout2_grid_local = codon_drout2_grid
     reff_rain_grid_local = codon_reff_rain_grid
     des_grid_local = codon_des_grid
     dsout2_grid_local = codon_dsout2_grid
     reff_snow_grid_local = codon_reff_snow_grid
     rei_grid_local = codon_rei_grid
     niic_grid_local = codon_niic_grid
     dei_grid_local = codon_dei_grid
     mgreffrain_grid_local = codon_mgreffrain_grid
     mgreffsnow_grid_local = codon_mgreffsnow_grid
  end if

  deallocate(native_mu_grid, codon_mu_grid)
  deallocate(native_lambdac_grid, codon_lambdac_grid)
  deallocate(native_rel_fn_grid, codon_rel_fn_grid)
  deallocate(native_ncic_grid, codon_ncic_grid)
  deallocate(native_rel_grid, codon_rel_grid)
  deallocate(native_drout2_grid, codon_drout2_grid)
  deallocate(native_reff_rain_grid, codon_reff_rain_grid)
  deallocate(native_des_grid, codon_des_grid)
  deallocate(native_dsout2_grid, codon_dsout2_grid)
  deallocate(native_reff_snow_grid, codon_reff_snow_grid)
  deallocate(native_rei_grid, codon_rei_grid)
  deallocate(native_niic_grid, codon_niic_grid)
  deallocate(native_dei_grid, codon_dei_grid)
  deallocate(native_mgreffrain_grid, codon_mgreffrain_grid)
  deallocate(native_mgreffsnow_grid, codon_mgreffsnow_grid)

contains

  subroutine report_field(name, native_field, codon_field)

    character(len=*), intent(in) :: name
    real(r8), intent(in) :: native_field(pcols,pver), codon_field(pcols,pver)
    integer :: i, k, diff_count, first_i, first_k
    real(r8) :: max_abs_diff, abs_diff, first_native, first_codon
    character(len=512) :: line

    diff_count = 0
    first_i = 0
    first_k = 0
    max_abs_diff = 0._r8
    first_native = 0._r8
    first_codon = 0._r8

    do k = 1, pver
       do i = 1, pcols
          if (native_field(i,k) /= codon_field(i,k)) then
             diff_count = diff_count + 1
             if (first_i == 0) then
                first_i = i
                first_k = k
                first_native = native_field(i,k)
                first_codon = codon_field(i,k)
             end if
             abs_diff = abs(native_field(i,k) - codon_field(i,k))
             if (abs_diff == abs_diff) then
                max_abs_diff = max(max_abs_diff, abs_diff)
             end if
          end if
       end do
    end do

    if (diff_count == 0) then
      write(line,'(A,": diff_count=0")') trim(name)
    else
      write(line,'(A,": diff_count=",I0,", max_abs_diff=",ES24.16E3,", first_diff=(",I0,",",I0,")",", native=",ES24.16E3,", codon=",ES24.16E3)') &
           trim(name), diff_count, max_abs_diff, first_i, first_k, first_native, first_codon
    end if
    call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_REFF_COMPARE_FILE', trim(line))

  end subroutine report_field

end subroutine micro_mg_cam_reff_calc_compare

subroutine micro_mg_cam_reff_liq_native(ngrdcol_local, rho_grid_local, icwmrst_grid_local, liqcldf_grid_local, nc_grid_local, &
     mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local)

  use micro_mg_utils, only: mg_liq_props, size_dist_param_liq, qsmall, mincld
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ngrdcol_local
  real(r8), intent(in) :: rho_grid_local(pcols,pver), icwmrst_grid_local(pcols,pver)
  real(r8), intent(in) :: liqcldf_grid_local(pcols,pver), nc_grid_local(pcols,pver)
  real(r8), intent(inout) :: mu_grid_local(pcols,pver), lambdac_grid_local(pcols,pver)
  real(r8), intent(inout) :: rel_fn_grid_local(pcols,pver), ncic_grid_local(pcols,pver)
  real(r8), intent(inout) :: rel_grid_local(pcols,pver)

  mu_grid_local = 0._r8
  lambdac_grid_local = 0._r8
  rel_fn_grid_local = 10._r8

  ncic_grid_local = 1.e8_r8
  call size_dist_param_liq(mg_liq_props, icwmrst_grid_local(:ngrdcol_local,top_lev:), ncic_grid_local(:ngrdcol_local,top_lev:), &
       rho_grid_local(:ngrdcol_local,top_lev:), mu_grid_local(:ngrdcol_local,top_lev:), lambdac_grid_local(:ngrdcol_local,top_lev:))
  where (icwmrst_grid_local(:ngrdcol_local,top_lev:) > qsmall)
     rel_fn_grid_local(:ngrdcol_local,top_lev:) = (mu_grid_local(:ngrdcol_local,top_lev:) + 3._r8) / &
          lambdac_grid_local(:ngrdcol_local,top_lev:) / 2._r8 * 1.e6_r8
  end where

  mu_grid_local = 0._r8
  lambdac_grid_local = 0._r8
  rel_grid_local = 10._r8

  ncic_grid_local(:ngrdcol_local,top_lev:) = nc_grid_local(:ngrdcol_local,top_lev:) / max(mincld, liqcldf_grid_local(:ngrdcol_local,top_lev:))
  call size_dist_param_liq(mg_liq_props, icwmrst_grid_local(:ngrdcol_local,top_lev:), ncic_grid_local(:ngrdcol_local,top_lev:), &
       rho_grid_local(:ngrdcol_local,top_lev:), mu_grid_local(:ngrdcol_local,top_lev:), lambdac_grid_local(:ngrdcol_local,top_lev:))
  where (icwmrst_grid_local(:ngrdcol_local,top_lev:) >= qsmall)
     rel_grid_local(:ngrdcol_local,top_lev:) = (mu_grid_local(:ngrdcol_local,top_lev:) + 3._r8) / &
          lambdac_grid_local(:ngrdcol_local,top_lev:) / 2._r8 * 1.e6_r8
  elsewhere
     mu_grid_local(:ngrdcol_local,top_lev:) = 0._r8
  end where

end subroutine micro_mg_cam_reff_liq_native

subroutine micro_mg_cam_reff_calc_native(ngrdcol_local, micro_mg_version_local, rho_grid_local, icwmrst_grid_local, &
     liqcldf_grid_local, nc_grid_local, qr_grid_local, nr_grid_local, qs_grid_local, ns_grid_local, qrout_grid_local, &
     nrout_grid_local, qsout_grid_local, nsout_grid_local, ni_grid_local, icecldf_grid_local, icimrst_grid_local, &
     ast_grid_local, mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local, &
     drout2_grid_local, reff_rain_grid_local, des_grid_local, dsout2_grid_local, reff_snow_grid_local, rei_grid_local, &
     niic_grid_local, dei_grid_local, mgreffrain_grid_local, mgreffsnow_grid_local)

  use micro_mg_utils, only: mg_liq_props, mg_ice_props, size_dist_param_liq, size_dist_param_basic, avg_diameter, qsmall, &
       mincld, rhosn, rhoi, rhow, rhows
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ngrdcol_local, micro_mg_version_local
  real(r8), parameter :: dcon_local = 25.e-6_r8
  real(r8), parameter :: mucon_local = 5.3_r8
  real(r8), parameter :: deicon_local = 50._r8
  real(r8), intent(in) :: rho_grid_local(pcols,pver), icwmrst_grid_local(pcols,pver)
  real(r8), intent(in) :: liqcldf_grid_local(pcols,pver), nc_grid_local(pcols,pver)
  real(r8), intent(in) :: qr_grid_local(pcols,pver), nr_grid_local(pcols,pver)
  real(r8), intent(in) :: qs_grid_local(pcols,pver), ns_grid_local(pcols,pver)
  real(r8), intent(in) :: qrout_grid_local(pcols,pver), nrout_grid_local(pcols,pver)
  real(r8), intent(in) :: qsout_grid_local(pcols,pver), nsout_grid_local(pcols,pver)
  real(r8), intent(in) :: ni_grid_local(pcols,pver), icecldf_grid_local(pcols,pver)
  real(r8), intent(in) :: icimrst_grid_local(pcols,pver), ast_grid_local(pcols,pver)
  real(r8), intent(inout) :: mu_grid_local(pcols,pver), lambdac_grid_local(pcols,pver)
  real(r8), intent(inout) :: rel_fn_grid_local(pcols,pver), ncic_grid_local(pcols,pver)
  real(r8), intent(inout) :: rel_grid_local(pcols,pver), drout2_grid_local(pcols,pver)
  real(r8), intent(inout) :: reff_rain_grid_local(pcols,pver), des_grid_local(pcols,pver)
  real(r8), intent(inout) :: dsout2_grid_local(pcols,pver), reff_snow_grid_local(pcols,pver)
  real(r8), intent(inout) :: rei_grid_local(pcols,pver), niic_grid_local(pcols,pver)
  real(r8), intent(inout) :: dei_grid_local(pcols,pver), mgreffrain_grid_local(pcols,pver)
  real(r8), intent(inout) :: mgreffsnow_grid_local(pcols,pver)
  integer :: i, k

  call micro_mg_cam_reff_liq_native(ngrdcol_local, rho_grid_local, icwmrst_grid_local, liqcldf_grid_local, nc_grid_local, &
       mu_grid_local, lambdac_grid_local, rel_fn_grid_local, ncic_grid_local, rel_grid_local)

  drout2_grid_local = 0._r8
  reff_rain_grid_local = 0._r8
  des_grid_local = 0._r8
  dsout2_grid_local = 0._r8
  reff_snow_grid_local = 0._r8

  if (micro_mg_version_local > 1) then
     where (qr_grid_local(:ngrdcol_local,top_lev:) >= 1.e-7_r8)
        drout2_grid_local(:ngrdcol_local,top_lev:) = avg_diameter(qr_grid_local(:ngrdcol_local,top_lev:), &
             nr_grid_local(:ngrdcol_local,top_lev:) * rho_grid_local(:ngrdcol_local,top_lev:), &
             rho_grid_local(:ngrdcol_local,top_lev:), rhow)
        reff_rain_grid_local(:ngrdcol_local,top_lev:) = drout2_grid_local(:ngrdcol_local,top_lev:) * 1.5_r8 * 1.e6_r8
     end where

     where (qs_grid_local(:ngrdcol_local,top_lev:) >= 1.e-7_r8)
        dsout2_grid_local(:ngrdcol_local,top_lev:) = avg_diameter(qs_grid_local(:ngrdcol_local,top_lev:), &
             ns_grid_local(:ngrdcol_local,top_lev:) * rho_grid_local(:ngrdcol_local,top_lev:), &
             rho_grid_local(:ngrdcol_local,top_lev:), rhosn)
        des_grid_local(:ngrdcol_local,top_lev:) = dsout2_grid_local(:ngrdcol_local,top_lev:) * 3._r8 * rhosn/rhows
        reff_snow_grid_local(:ngrdcol_local,top_lev:) = dsout2_grid_local(:ngrdcol_local,top_lev:) * 1.5_r8 * 1.e6_r8
     end where
  else
     where (qrout_grid_local(:ngrdcol_local,top_lev:) >= 1.e-7_r8)
        drout2_grid_local(:ngrdcol_local,top_lev:) = avg_diameter(qrout_grid_local(:ngrdcol_local,top_lev:), &
             nrout_grid_local(:ngrdcol_local,top_lev:) * rho_grid_local(:ngrdcol_local,top_lev:), &
             rho_grid_local(:ngrdcol_local,top_lev:), rhow)
        reff_rain_grid_local(:ngrdcol_local,top_lev:) = drout2_grid_local(:ngrdcol_local,top_lev:) * 1.5_r8 * 1.e6_r8
     end where

     where (qsout_grid_local(:ngrdcol_local,top_lev:) >= 1.e-7_r8)
        dsout2_grid_local(:ngrdcol_local,top_lev:) = avg_diameter(qsout_grid_local(:ngrdcol_local,top_lev:), &
             nsout_grid_local(:ngrdcol_local,top_lev:) * rho_grid_local(:ngrdcol_local,top_lev:), &
             rho_grid_local(:ngrdcol_local,top_lev:), rhosn)
        des_grid_local(:ngrdcol_local,top_lev:) = dsout2_grid_local(:ngrdcol_local,top_lev:) * 3._r8 * rhosn/rhows
        reff_snow_grid_local(:ngrdcol_local,top_lev:) = dsout2_grid_local(:ngrdcol_local,top_lev:) * 1.5_r8 * 1.e6_r8
     end where
  end if

  rei_grid_local = 25._r8
  niic_grid_local(:ngrdcol_local,top_lev:) = ni_grid_local(:ngrdcol_local,top_lev:) / max(mincld, icecldf_grid_local(:ngrdcol_local,top_lev:))
  call size_dist_param_basic(mg_ice_props, icimrst_grid_local(:ngrdcol_local,top_lev:), niic_grid_local(:ngrdcol_local,top_lev:), &
       rei_grid_local(:ngrdcol_local,top_lev:))
  where (icimrst_grid_local(:ngrdcol_local,top_lev:) >= qsmall)
     rei_grid_local(:ngrdcol_local,top_lev:) = 1.5_r8 / rei_grid_local(:ngrdcol_local,top_lev:) * 1.e6_r8
  elsewhere
     rei_grid_local(:ngrdcol_local,top_lev:) = 25._r8
  end where

  dei_grid_local = rei_grid_local * rhoi/rhows * 2._r8
  do k = top_lev, pver
     do i = 1, ngrdcol_local
        des_grid_local(i,k) = des_grid_local(i,k) * 1.e6_r8
        if (ast_grid_local(i,k) < 1.e-4_r8) then
           mu_grid_local(i,k) = mucon_local
           lambdac_grid_local(i,k) = (mucon_local + 1._r8) / dcon_local
           dei_grid_local(i,k) = deicon_local
        end if
     end do
  end do

  mgreffrain_grid_local(:ngrdcol_local,top_lev:pver) = reff_rain_grid_local(:ngrdcol_local,top_lev:pver)
  mgreffsnow_grid_local(:ngrdcol_local,top_lev:pver) = reff_snow_grid_local(:ngrdcol_local,top_lev:pver)

end subroutine micro_mg_cam_reff_calc_native

subroutine micro_mg_cam_select_postmg_diag_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (postmg_diag_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MICRO_MG_CAM_POSTMG_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_postmg_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_postmg_diag_impl = .false.
  end if

  postmg_diag_impl_selected = .true.

  if (use_native_postmg_diag_impl) then
     write(iulog,*) 'micro_mg_cam_postmg_diag implementation = native'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_POSTMG_PROOF_FILE', &
          'micro_mg_cam_postmg_diag implementation = native')
  else
     write(iulog,*) 'micro_mg_cam_postmg_diag implementation = codon'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_POSTMG_PROOF_FILE', &
          'micro_mg_cam_postmg_diag implementation = codon')
  end if
  call flush(iulog)

end subroutine micro_mg_cam_select_postmg_diag_impl

subroutine micro_mg_cam_postmg_diag(ncol_local, psetcols_local, micro_mg_version_local, rate1_cw2pr_st_idx_local, &
     ixcldliq_local, ixcldice_local, ixnumliq_local, ixnumice_local, ixrain_local, ixsnow_local, state_q_local, &
     state_t_local, state_pmid_local, state_pdel_local, naai_local, naai_hom_local, mnuccdo_local, rflx_local, sflx_local, &
     qrout_local, qsout_local, prect_local, preci_local, rate1cld_local, vtrmc_local, tlat_local, qvlat_local, qcten_local, &
     qiten_local, ncten_local, niten_local, alst_mic_local, cmeliq_local, cmeiout_local, ast_local, cld_local, &
     concld_local, mnuccdohet_local, mgflxprc_local, mgflxsnw_local, mgmrprc_local, mgmrsnw_local, cvreffliq_local, &
     cvreffice_local, rate1ord_cw2pr_st_local, wsedl_local, cc_t_local, cc_qv_local, cc_ql_local, cc_qi_local, cc_nl_local, &
     cc_ni_local, cc_qlst_local, qme_local, prec_pcw_local, snow_pcw_local, prec_sed_local, snow_sed_local, prec_str_local, &
     snow_str_local, icecldf_local, liqcldf_local, icinc_local, icwnc_local, iciwpst_local, iclwpst_local, icswp_local, &
     cldfsnow_local, icimrst_local, icwmrst_local, cldmax_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use micro_mg_utils, only: qsmall, mincld
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local, psetcols_local, micro_mg_version_local, rate1_cw2pr_st_idx_local
  integer, intent(in) :: ixcldliq_local, ixcldice_local, ixnumliq_local, ixnumice_local, ixrain_local, ixsnow_local
  real(r8), target, intent(in) :: state_q_local(psetcols_local,pver,pcnst)
  real(r8), target, intent(in) :: state_t_local(psetcols_local,pver), state_pmid_local(psetcols_local,pver)
  real(r8), target, intent(in) :: state_pdel_local(psetcols_local,pver)
  real(r8), target, intent(in) :: naai_local(psetcols_local,pver), naai_hom_local(psetcols_local,pver)
  real(r8), target, intent(in) :: mnuccdo_local(psetcols_local,pver), rflx_local(psetcols_local,pverp)
  real(r8), target, intent(in) :: sflx_local(psetcols_local,pverp), qrout_local(psetcols_local,pver)
  real(r8), target, intent(in) :: qsout_local(psetcols_local,pver), prect_local(psetcols_local), preci_local(psetcols_local)
  real(r8), target, intent(in) :: rate1cld_local(psetcols_local,pver), vtrmc_local(psetcols_local,pver)
  real(r8), target, intent(in) :: tlat_local(psetcols_local,pver), qvlat_local(psetcols_local,pver)
  real(r8), target, intent(in) :: qcten_local(psetcols_local,pver), qiten_local(psetcols_local,pver)
  real(r8), target, intent(in) :: ncten_local(psetcols_local,pver), niten_local(psetcols_local,pver)
  real(r8), target, intent(in) :: alst_mic_local(psetcols_local,pver), cmeliq_local(psetcols_local,pver)
  real(r8), target, intent(in) :: cmeiout_local(psetcols_local,pver), ast_local(psetcols_local,pver)
  real(r8), target, intent(in) :: cld_local(psetcols_local,pver), concld_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: mnuccdohet_local(psetcols_local,pver), mgflxprc_local(psetcols_local,pverp)
  real(r8), target, intent(inout) :: mgflxsnw_local(psetcols_local,pverp), mgmrprc_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: mgmrsnw_local(psetcols_local,pver), cvreffliq_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: cvreffice_local(psetcols_local,pver), rate1ord_cw2pr_st_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: wsedl_local(psetcols_local,pver), cc_t_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: cc_qv_local(psetcols_local,pver), cc_ql_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: cc_qi_local(psetcols_local,pver), cc_nl_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: cc_ni_local(psetcols_local,pver), cc_qlst_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: qme_local(psetcols_local,pver), prec_pcw_local(psetcols_local)
  real(r8), target, intent(inout) :: snow_pcw_local(psetcols_local), prec_sed_local(psetcols_local)
  real(r8), target, intent(inout) :: snow_sed_local(psetcols_local), prec_str_local(psetcols_local)
  real(r8), target, intent(inout) :: snow_str_local(psetcols_local), icecldf_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: liqcldf_local(psetcols_local,pver), icinc_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: icwnc_local(psetcols_local,pver), iciwpst_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: iclwpst_local(psetcols_local,pver), icswp_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: cldfsnow_local(psetcols_local,pver), icimrst_local(psetcols_local,pver)
  real(r8), target, intent(inout) :: icwmrst_local(psetcols_local,pver), cldmax_local(psetcols_local,pver)

  interface
     subroutine micro_mg_cam_postmg_diag_codon(ncol_c, psetcols_c, pver_c, pverp_c, top_lev_c, micro_mg_version_c, &
          rate1_cw2pr_st_idx_c, ixcldliq_c, ixcldice_c, ixnumliq_c, ixnumice_c, ixrain_c, ixsnow_c, cpair_c, gravit_c, &
          mincld_c, qsmall_c, state_q_p, state_t_p, state_pmid_p, state_pdel_p, naai_p, naai_hom_p, mnuccdo_p, rflx_p, &
          sflx_p, qrout_p, qsout_p, prect_p, preci_p, rate1cld_p, vtrmc_p, tlat_p, qvlat_p, qcten_p, qiten_p, ncten_p, &
          niten_p, alst_mic_p, cmeliq_p, cmeiout_p, ast_p, cld_p, concld_p, mnuccdohet_p, mgflxprc_p, mgflxsnw_p, &
          mgmrprc_p, mgmrsnw_p, cvreffliq_p, cvreffice_p, rate1ord_cw2pr_st_p, wsedl_p, cc_t_p, cc_qv_p, cc_ql_p, cc_qi_p, &
          cc_nl_p, cc_ni_p, cc_qlst_p, qme_p, prec_pcw_p, snow_pcw_p, prec_sed_p, snow_sed_p, prec_str_p, snow_str_p, &
          icecldf_p, liqcldf_p, icinc_p, icwnc_p, iciwpst_p, iclwpst_p, icswp_p, cldfsnow_p, icimrst_p, icwmrst_p, &
          cldmax_p) bind(c, name="micro_mg_cam_postmg_diag_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, pverp_c, top_lev_c, micro_mg_version_c
       integer(c_int64_t), value :: rate1_cw2pr_st_idx_c, ixcldliq_c, ixcldice_c, ixnumliq_c, ixnumice_c, ixrain_c
       integer(c_int64_t), value :: ixsnow_c
       real(c_double), value :: cpair_c, gravit_c, mincld_c, qsmall_c
       type(c_ptr), value :: state_q_p, state_t_p, state_pmid_p, state_pdel_p, naai_p, naai_hom_p, mnuccdo_p, rflx_p
       type(c_ptr), value :: sflx_p, qrout_p, qsout_p, prect_p, preci_p, rate1cld_p, vtrmc_p, tlat_p, qvlat_p, qcten_p
       type(c_ptr), value :: qiten_p, ncten_p, niten_p, alst_mic_p, cmeliq_p, cmeiout_p, ast_p, cld_p, concld_p
       type(c_ptr), value :: mnuccdohet_p, mgflxprc_p, mgflxsnw_p, mgmrprc_p, mgmrsnw_p, cvreffliq_p, cvreffice_p
       type(c_ptr), value :: rate1ord_cw2pr_st_p, wsedl_p, cc_t_p, cc_qv_p, cc_ql_p, cc_qi_p, cc_nl_p, cc_ni_p
       type(c_ptr), value :: cc_qlst_p, qme_p, prec_pcw_p, snow_pcw_p, prec_sed_p, snow_sed_p, prec_str_p, snow_str_p
       type(c_ptr), value :: icecldf_p, liqcldf_p, icinc_p, icwnc_p, iciwpst_p, iclwpst_p, icswp_p, cldfsnow_p
       type(c_ptr), value :: icimrst_p, icwmrst_p, cldmax_p
     end subroutine micro_mg_cam_postmg_diag_codon
  end interface

  call micro_mg_cam_select_postmg_diag_impl()

  if (use_native_postmg_diag_impl) then
     call micro_mg_cam_postmg_diag_native(ncol_local, psetcols_local, micro_mg_version_local, rate1_cw2pr_st_idx_local, &
          ixcldliq_local, ixcldice_local, ixnumliq_local, ixnumice_local, ixrain_local, ixsnow_local, state_q_local, &
          state_t_local, state_pmid_local, state_pdel_local, naai_local, naai_hom_local, mnuccdo_local, rflx_local, &
          sflx_local, qrout_local, qsout_local, prect_local, preci_local, rate1cld_local, vtrmc_local, tlat_local, &
          qvlat_local, qcten_local, qiten_local, ncten_local, niten_local, alst_mic_local, cmeliq_local, cmeiout_local, &
          ast_local, cld_local, concld_local, mnuccdohet_local, mgflxprc_local, mgflxsnw_local, mgmrprc_local, &
          mgmrsnw_local, cvreffliq_local, cvreffice_local, rate1ord_cw2pr_st_local, wsedl_local, cc_t_local, cc_qv_local, &
          cc_ql_local, cc_qi_local, cc_nl_local, cc_ni_local, cc_qlst_local, qme_local, prec_pcw_local, snow_pcw_local, &
          prec_sed_local, snow_sed_local, prec_str_local, snow_str_local, icecldf_local, liqcldf_local, icinc_local, &
          icwnc_local, iciwpst_local, iclwpst_local, icswp_local, cldfsnow_local, icimrst_local, icwmrst_local, cldmax_local)
     return
  end if

  call micro_mg_cam_postmg_diag_codon(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
       int(pverp, c_int64_t), int(top_lev, c_int64_t), int(micro_mg_version_local, c_int64_t), &
       int(rate1_cw2pr_st_idx_local, c_int64_t), int(ixcldliq_local, c_int64_t), int(ixcldice_local, c_int64_t), &
       int(ixnumliq_local, c_int64_t), int(ixnumice_local, c_int64_t), int(ixrain_local, c_int64_t), &
       int(ixsnow_local, c_int64_t), cpair, gravit, mincld, qsmall, c_loc(state_q_local), c_loc(state_t_local), &
       c_loc(state_pmid_local), c_loc(state_pdel_local), c_loc(naai_local), c_loc(naai_hom_local), c_loc(mnuccdo_local), &
       c_loc(rflx_local), c_loc(sflx_local), c_loc(qrout_local), c_loc(qsout_local), c_loc(prect_local), c_loc(preci_local), &
       c_loc(rate1cld_local), c_loc(vtrmc_local), c_loc(tlat_local), c_loc(qvlat_local), c_loc(qcten_local), &
       c_loc(qiten_local), c_loc(ncten_local), c_loc(niten_local), c_loc(alst_mic_local), c_loc(cmeliq_local), &
       c_loc(cmeiout_local), c_loc(ast_local), c_loc(cld_local), c_loc(concld_local), c_loc(mnuccdohet_local), &
       c_loc(mgflxprc_local), c_loc(mgflxsnw_local), c_loc(mgmrprc_local), c_loc(mgmrsnw_local), c_loc(cvreffliq_local), &
       c_loc(cvreffice_local), c_loc(rate1ord_cw2pr_st_local), c_loc(wsedl_local), c_loc(cc_t_local), c_loc(cc_qv_local), &
       c_loc(cc_ql_local), c_loc(cc_qi_local), c_loc(cc_nl_local), c_loc(cc_ni_local), c_loc(cc_qlst_local), c_loc(qme_local), &
       c_loc(prec_pcw_local), c_loc(snow_pcw_local), c_loc(prec_sed_local), c_loc(snow_sed_local), c_loc(prec_str_local), &
       c_loc(snow_str_local), c_loc(icecldf_local), c_loc(liqcldf_local), c_loc(icinc_local), c_loc(icwnc_local), &
       c_loc(iciwpst_local), c_loc(iclwpst_local), c_loc(icswp_local), c_loc(cldfsnow_local), c_loc(icimrst_local), &
       c_loc(icwmrst_local), c_loc(cldmax_local))

end subroutine micro_mg_cam_postmg_diag

subroutine micro_mg_cam_postmg_diag_native(ncol_local, psetcols_local, micro_mg_version_local, rate1_cw2pr_st_idx_local, &
     ixcldliq_local, ixcldice_local, ixnumliq_local, ixnumice_local, ixrain_local, ixsnow_local, state_q_local, &
     state_t_local, state_pmid_local, state_pdel_local, naai_local, naai_hom_local, mnuccdo_local, rflx_local, sflx_local, &
     qrout_local, qsout_local, prect_local, preci_local, rate1cld_local, vtrmc_local, tlat_local, qvlat_local, qcten_local, &
     qiten_local, ncten_local, niten_local, alst_mic_local, cmeliq_local, cmeiout_local, ast_local, cld_local, &
     concld_local, mnuccdohet_local, mgflxprc_local, mgflxsnw_local, mgmrprc_local, mgmrsnw_local, cvreffliq_local, &
     cvreffice_local, rate1ord_cw2pr_st_local, wsedl_local, cc_t_local, cc_qv_local, cc_ql_local, cc_qi_local, cc_nl_local, &
     cc_ni_local, cc_qlst_local, qme_local, prec_pcw_local, snow_pcw_local, prec_sed_local, snow_sed_local, prec_str_local, &
     snow_str_local, icecldf_local, liqcldf_local, icinc_local, icwnc_local, iciwpst_local, iclwpst_local, icswp_local, &
     cldfsnow_local, icimrst_local, icwmrst_local, cldmax_local)

  use micro_mg_utils, only: qsmall, mincld
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ncol_local, psetcols_local, micro_mg_version_local, rate1_cw2pr_st_idx_local
  integer, intent(in) :: ixcldliq_local, ixcldice_local, ixnumliq_local, ixnumice_local, ixrain_local, ixsnow_local
  real(r8), intent(in) :: state_q_local(psetcols_local,pver,pcnst)
  real(r8), intent(in) :: state_t_local(psetcols_local,pver), state_pmid_local(psetcols_local,pver)
  real(r8), intent(in) :: state_pdel_local(psetcols_local,pver)
  real(r8), intent(in) :: naai_local(psetcols_local,pver), naai_hom_local(psetcols_local,pver)
  real(r8), intent(in) :: mnuccdo_local(psetcols_local,pver), rflx_local(psetcols_local,pverp)
  real(r8), intent(in) :: sflx_local(psetcols_local,pverp), qrout_local(psetcols_local,pver)
  real(r8), intent(in) :: qsout_local(psetcols_local,pver), prect_local(psetcols_local), preci_local(psetcols_local)
  real(r8), intent(in) :: rate1cld_local(psetcols_local,pver), vtrmc_local(psetcols_local,pver)
  real(r8), intent(in) :: tlat_local(psetcols_local,pver), qvlat_local(psetcols_local,pver)
  real(r8), intent(in) :: qcten_local(psetcols_local,pver), qiten_local(psetcols_local,pver)
  real(r8), intent(in) :: ncten_local(psetcols_local,pver), niten_local(psetcols_local,pver)
  real(r8), intent(in) :: alst_mic_local(psetcols_local,pver), cmeliq_local(psetcols_local,pver)
  real(r8), intent(in) :: cmeiout_local(psetcols_local,pver), ast_local(psetcols_local,pver)
  real(r8), intent(in) :: cld_local(psetcols_local,pver), concld_local(psetcols_local,pver)
  real(r8), intent(inout) :: mnuccdohet_local(psetcols_local,pver), mgflxprc_local(psetcols_local,pverp)
  real(r8), intent(inout) :: mgflxsnw_local(psetcols_local,pverp), mgmrprc_local(psetcols_local,pver)
  real(r8), intent(inout) :: mgmrsnw_local(psetcols_local,pver), cvreffliq_local(psetcols_local,pver)
  real(r8), intent(inout) :: cvreffice_local(psetcols_local,pver), rate1ord_cw2pr_st_local(psetcols_local,pver)
  real(r8), intent(inout) :: wsedl_local(psetcols_local,pver), cc_t_local(psetcols_local,pver)
  real(r8), intent(inout) :: cc_qv_local(psetcols_local,pver), cc_ql_local(psetcols_local,pver)
  real(r8), intent(inout) :: cc_qi_local(psetcols_local,pver), cc_nl_local(psetcols_local,pver)
  real(r8), intent(inout) :: cc_ni_local(psetcols_local,pver), cc_qlst_local(psetcols_local,pver)
  real(r8), intent(inout) :: qme_local(psetcols_local,pver), prec_pcw_local(psetcols_local)
  real(r8), intent(inout) :: snow_pcw_local(psetcols_local), prec_sed_local(psetcols_local)
  real(r8), intent(inout) :: snow_sed_local(psetcols_local), prec_str_local(psetcols_local)
  real(r8), intent(inout) :: snow_str_local(psetcols_local), icecldf_local(psetcols_local,pver)
  real(r8), intent(inout) :: liqcldf_local(psetcols_local,pver), icinc_local(psetcols_local,pver)
  real(r8), intent(inout) :: icwnc_local(psetcols_local,pver), iciwpst_local(psetcols_local,pver)
  real(r8), intent(inout) :: iclwpst_local(psetcols_local,pver), icswp_local(psetcols_local,pver)
  real(r8), intent(inout) :: cldfsnow_local(psetcols_local,pver), icimrst_local(psetcols_local,pver)
  real(r8), intent(inout) :: icwmrst_local(psetcols_local,pver), cldmax_local(psetcols_local,pver)
  integer :: i, k

  mnuccdohet_local = 0._r8
  do k = top_lev, pver
     do i = 1, ncol_local
        if (naai_local(i,k) > 0._r8) then
           mnuccdohet_local(i,k) = mnuccdo_local(i,k) - (naai_hom_local(i,k)/naai_local(i,k))*mnuccdo_local(i,k)
        end if
     end do
  end do

  mgflxprc_local(:ncol_local,top_lev:pverp) = rflx_local(:ncol_local,top_lev:pverp) + sflx_local(:ncol_local,top_lev:pverp)
  mgflxsnw_local(:ncol_local,top_lev:pverp) = sflx_local(:ncol_local,top_lev:pverp)

  mgmrprc_local(:ncol_local,top_lev:pver) = qrout_local(:ncol_local,top_lev:pver) + qsout_local(:ncol_local,top_lev:pver)
  mgmrsnw_local(:ncol_local,top_lev:pver) = qsout_local(:ncol_local,top_lev:pver)

  cvreffliq_local(:ncol_local,top_lev:pver) = 9.0_r8
  cvreffice_local(:ncol_local,top_lev:pver) = 37.0_r8

  if (rate1_cw2pr_st_idx_local > 0) then
     rate1ord_cw2pr_st_local(:ncol_local,top_lev:pver) = rate1cld_local(:ncol_local,top_lev:pver)
  end if

  wsedl_local(:ncol_local,top_lev:pver) = vtrmc_local(:ncol_local,top_lev:pver)

  cc_t_local(:ncol_local,top_lev:pver)    = tlat_local(:ncol_local,top_lev:pver)/cpair
  cc_qv_local(:ncol_local,top_lev:pver)   = qvlat_local(:ncol_local,top_lev:pver)
  cc_ql_local(:ncol_local,top_lev:pver)   = qcten_local(:ncol_local,top_lev:pver)
  cc_qi_local(:ncol_local,top_lev:pver)   = qiten_local(:ncol_local,top_lev:pver)
  cc_nl_local(:ncol_local,top_lev:pver)   = ncten_local(:ncol_local,top_lev:pver)
  cc_ni_local(:ncol_local,top_lev:pver)   = niten_local(:ncol_local,top_lev:pver)
  cc_qlst_local(:ncol_local,top_lev:pver) = qcten_local(:ncol_local,top_lev:pver) / &
       max(0.01_r8, alst_mic_local(:ncol_local,top_lev:pver))

  qme_local(:ncol_local,top_lev:pver) = cmeliq_local(:ncol_local,top_lev:pver) + cmeiout_local(:ncol_local,top_lev:pver)

  prec_pcw_local = prect_local
  snow_pcw_local = preci_local
  prec_sed_local = 0._r8
  snow_sed_local = 0._r8
  prec_str_local = prec_pcw_local + prec_sed_local
  snow_str_local = snow_pcw_local + snow_sed_local

  icecldf_local(:ncol_local,top_lev:pver) = ast_local(:ncol_local,top_lev:pver)
  liqcldf_local(:ncol_local,top_lev:pver) = ast_local(:ncol_local,top_lev:pver)

  icinc_local = 0._r8
  icwnc_local = 0._r8
  iciwpst_local = 0._r8
  iclwpst_local = 0._r8
  icswp_local = 0._r8
  cldfsnow_local = 0._r8

  do k = top_lev, pver
     do i = 1, ncol_local
        icimrst_local(i,k) = min(state_q_local(i,k,ixcldice_local) / max(mincld,icecldf_local(i,k)), 0.005_r8)
        icwmrst_local(i,k) = min(state_q_local(i,k,ixcldliq_local) / max(mincld,liqcldf_local(i,k)), 0.005_r8)
        icinc_local(i,k) = state_q_local(i,k,ixnumice_local) / max(mincld,icecldf_local(i,k)) * state_pmid_local(i,k) / &
             (287.15_r8*state_t_local(i,k))
        icwnc_local(i,k) = state_q_local(i,k,ixnumliq_local) / max(mincld,liqcldf_local(i,k)) * state_pmid_local(i,k) / &
             (287.15_r8*state_t_local(i,k))
        iciwpst_local(i,k) = min(state_q_local(i,k,ixcldice_local) / max(mincld,ast_local(i,k)), 0.005_r8) * &
             state_pdel_local(i,k) / gravit
        iclwpst_local(i,k) = min(state_q_local(i,k,ixcldliq_local) / max(mincld,ast_local(i,k)), 0.005_r8) * &
             state_pdel_local(i,k) / gravit

        cldfsnow_local(i,k) = cld_local(i,k)
        if ((cldfsnow_local(i,k) .gt. 1.e-4_r8) .and. (concld_local(i,k) .lt. 1.e-4_r8) .and. &
             (state_q_local(i,k,ixcldliq_local) .lt. 1.e-10_r8)) then
           cldfsnow_local(i,k) = 0._r8
        end if
        if ((cldfsnow_local(i,k) .le. 1.e-4_r8) .and. (qsout_local(i,k) .gt. 1.e-6_r8)) then
           cldfsnow_local(i,k) = 0.25_r8
        end if
        icswp_local(i,k) = qsout_local(i,k) / max(mincld, cldfsnow_local(i,k)) * state_pdel_local(i,k) / gravit
     end do
  end do

  if (micro_mg_version_local > 1) then
     cldmax_local = max(mincld, ast_local)
     do k = top_lev+1, pver
        where (state_q_local(:ncol_local,k-1,ixrain_local) >= qsmall .or. &
             state_q_local(:ncol_local,k-1,ixsnow_local) >= qsmall)
           cldmax_local(:ncol_local,k) = max(cldmax_local(:ncol_local,k-1), cldmax_local(:ncol_local,k))
        end where
     end do
  end if

end subroutine micro_mg_cam_postmg_diag_native

subroutine micro_mg_cam_select_grid_diag_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (grid_diag_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('MICRO_MG_CAM_GRID_DIAG_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_grid_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_grid_diag_impl = .false.
  end if

  grid_diag_impl_selected = .true.

  if (use_native_grid_diag_impl) then
     write(iulog,*) 'micro_mg_cam_grid_diag implementation = native'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_GRID_DIAG_PROOF_FILE', &
          'micro_mg_cam_grid_diag implementation = native')
  else
     write(iulog,*) 'micro_mg_cam_grid_diag implementation = codon'
     call micro_mg_cam_append_impl_proof('MICRO_MG_CAM_GRID_DIAG_PROOF_FILE', &
          'micro_mg_cam_grid_diag implementation = codon')
  end if
  call flush(iulog)

end subroutine micro_mg_cam_select_grid_diag_impl

subroutine micro_mg_cam_grid_diag(ngrdcol_local, minlwp_local, iclwpst_grid_local, cld_grid_local, cmeliq_grid_local, &
     pdel_grid_local, prec_str_grid_local, acgcme_grid_local, acprecl_grid_local, acnum_grid_local, prao_grid_local, &
     prco_grid_local, nc_grid_local, liqcldf_grid_local, icwmrst_grid_local, rel_grid_local, icwnc_grid_local, &
     icecldf_grid_local, icimrst_grid_local, rei_grid_local, icinc_grid_local, nevapr_grid_local, evpsnow_st_grid_local, &
     tgliqwp_grid_local, tgcmeliq_grid_local, pe_grid_local, tpr_grid_local, pefrac_grid_local, vprao_grid_local, &
     vprco_grid_local, racau_grid_local, cnt_grid_local, cdnumc_grid_local, efcout_grid_local, efiout_grid_local, &
     ncout_grid_local, niout_grid_local, freql_grid_local, freqi_grid_local, icwmrst_grid_out_local, &
     icimrst_grid_out_local, fcti_grid_local, fctl_grid_local, ctrel_grid_local, ctrei_grid_local, ctnl_grid_local, &
     ctni_grid_local, evprain_st_grid_local)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ngrdcol_local
  real(r8), intent(in) :: minlwp_local
  real(r8), target, intent(in) :: iclwpst_grid_local(pcols,pver), cld_grid_local(pcols,pver)
  real(r8), target, intent(in) :: cmeliq_grid_local(pcols,pver), pdel_grid_local(pcols,pver)
  real(r8), target, intent(in) :: prec_str_grid_local(pcols), prao_grid_local(pcols,pver), prco_grid_local(pcols,pver)
  real(r8), target, intent(in) :: nc_grid_local(pcols,pver), liqcldf_grid_local(pcols,pver), icwmrst_grid_local(pcols,pver)
  real(r8), target, intent(in) :: rel_grid_local(pcols,pver), icwnc_grid_local(pcols,pver), icecldf_grid_local(pcols,pver)
  real(r8), target, intent(in) :: icimrst_grid_local(pcols,pver), rei_grid_local(pcols,pver), icinc_grid_local(pcols,pver)
  real(r8), target, intent(in) :: nevapr_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: evpsnow_st_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: acgcme_grid_local(:), acprecl_grid_local(:)
  integer, target, intent(inout) :: acnum_grid_local(:)
  real(r8), target, intent(inout) :: tgliqwp_grid_local(pcols), tgcmeliq_grid_local(pcols), pe_grid_local(pcols)
  real(r8), target, intent(inout) :: tpr_grid_local(pcols), pefrac_grid_local(pcols), vprao_grid_local(pcols)
  real(r8), target, intent(inout) :: vprco_grid_local(pcols), racau_grid_local(pcols), cdnumc_grid_local(pcols)
  integer, target, intent(inout) :: cnt_grid_local(pcols)
  real(r8), target, intent(inout) :: efcout_grid_local(pcols,pver), efiout_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: ncout_grid_local(pcols,pver), niout_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: freql_grid_local(pcols,pver), freqi_grid_local(pcols,pver)
  real(r8), target, intent(inout) :: icwmrst_grid_out_local(pcols,pver), icimrst_grid_out_local(pcols,pver)
  real(r8), target, intent(inout) :: fcti_grid_local(pcols), fctl_grid_local(pcols), ctrel_grid_local(pcols)
  real(r8), target, intent(inout) :: ctrei_grid_local(pcols), ctnl_grid_local(pcols), ctni_grid_local(pcols)
  real(r8), target, intent(inout) :: evprain_st_grid_local(pcols,pver)

  interface
     subroutine micro_mg_cam_grid_diag_codon(ngrdcol_c, pcols_c, pver_c, top_lev_c, minlwp_c, gravit_c, rhoh2o_c, &
          iclwpst_grid_p, cld_grid_p, cmeliq_grid_p, pdel_grid_p, prec_str_grid_p, acgcme_grid_p, acprecl_grid_p, &
          acnum_grid_p, prao_grid_p, prco_grid_p, nc_grid_p, liqcldf_grid_p, icwmrst_grid_p, rel_grid_p, icwnc_grid_p, &
          icecldf_grid_p, icimrst_grid_p, rei_grid_p, icinc_grid_p, nevapr_grid_p, evpsnow_st_grid_p, tgliqwp_grid_p, &
          tgcmeliq_grid_p, pe_grid_p, tpr_grid_p, pefrac_grid_p, vprao_grid_p, vprco_grid_p, racau_grid_p, cnt_grid_p, &
          cdnumc_grid_p, efcout_grid_p, efiout_grid_p, ncout_grid_p, niout_grid_p, freql_grid_p, freqi_grid_p, &
          icwmrst_grid_out_p, icimrst_grid_out_p, fcti_grid_p, fctl_grid_p, ctrel_grid_p, ctrei_grid_p, ctnl_grid_p, &
          ctni_grid_p, evprain_st_grid_p) bind(c, name="micro_mg_cam_grid_diag_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ngrdcol_c, pcols_c, pver_c, top_lev_c
       real(c_double), value :: minlwp_c, gravit_c, rhoh2o_c
       type(c_ptr), value :: iclwpst_grid_p, cld_grid_p, cmeliq_grid_p, pdel_grid_p, prec_str_grid_p, acgcme_grid_p
       type(c_ptr), value :: acprecl_grid_p, acnum_grid_p, prao_grid_p, prco_grid_p, nc_grid_p, liqcldf_grid_p
       type(c_ptr), value :: icwmrst_grid_p, rel_grid_p, icwnc_grid_p, icecldf_grid_p, icimrst_grid_p, rei_grid_p
       type(c_ptr), value :: icinc_grid_p, nevapr_grid_p, evpsnow_st_grid_p, tgliqwp_grid_p, tgcmeliq_grid_p, pe_grid_p
       type(c_ptr), value :: tpr_grid_p, pefrac_grid_p, vprao_grid_p, vprco_grid_p, racau_grid_p, cnt_grid_p
       type(c_ptr), value :: cdnumc_grid_p, efcout_grid_p, efiout_grid_p, ncout_grid_p, niout_grid_p, freql_grid_p
       type(c_ptr), value :: freqi_grid_p, icwmrst_grid_out_p, icimrst_grid_out_p, fcti_grid_p, fctl_grid_p
       type(c_ptr), value :: ctrel_grid_p, ctrei_grid_p, ctnl_grid_p, ctni_grid_p, evprain_st_grid_p
     end subroutine micro_mg_cam_grid_diag_codon
  end interface

  call micro_mg_cam_select_grid_diag_impl()

  if (use_native_grid_diag_impl) then
     call micro_mg_cam_grid_diag_native(ngrdcol_local, minlwp_local, iclwpst_grid_local, cld_grid_local, cmeliq_grid_local, &
          pdel_grid_local, prec_str_grid_local, acgcme_grid_local, acprecl_grid_local, acnum_grid_local, prao_grid_local, &
          prco_grid_local, nc_grid_local, liqcldf_grid_local, icwmrst_grid_local, rel_grid_local, icwnc_grid_local, &
          icecldf_grid_local, icimrst_grid_local, rei_grid_local, icinc_grid_local, nevapr_grid_local, evpsnow_st_grid_local, &
          tgliqwp_grid_local, tgcmeliq_grid_local, pe_grid_local, tpr_grid_local, pefrac_grid_local, vprao_grid_local, &
          vprco_grid_local, racau_grid_local, cnt_grid_local, cdnumc_grid_local, efcout_grid_local, efiout_grid_local, &
          ncout_grid_local, niout_grid_local, freql_grid_local, freqi_grid_local, icwmrst_grid_out_local, &
          icimrst_grid_out_local, fcti_grid_local, fctl_grid_local, ctrel_grid_local, ctrei_grid_local, ctnl_grid_local, &
          ctni_grid_local, evprain_st_grid_local)
     return
  end if

  call micro_mg_cam_grid_diag_codon(int(ngrdcol_local, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(top_lev, c_int64_t), minlwp_local, gravit, rhoh2o, c_loc(iclwpst_grid_local), c_loc(cld_grid_local), &
       c_loc(cmeliq_grid_local), c_loc(pdel_grid_local), c_loc(prec_str_grid_local), c_loc(acgcme_grid_local), &
       c_loc(acprecl_grid_local), c_loc(acnum_grid_local(1)), c_loc(prao_grid_local), c_loc(prco_grid_local), &
       c_loc(nc_grid_local), c_loc(liqcldf_grid_local), c_loc(icwmrst_grid_local), c_loc(rel_grid_local), c_loc(icwnc_grid_local), &
       c_loc(icecldf_grid_local), c_loc(icimrst_grid_local), c_loc(rei_grid_local), c_loc(icinc_grid_local), &
       c_loc(nevapr_grid_local), c_loc(evpsnow_st_grid_local), c_loc(tgliqwp_grid_local), c_loc(tgcmeliq_grid_local), &
       c_loc(pe_grid_local), c_loc(tpr_grid_local), c_loc(pefrac_grid_local), c_loc(vprao_grid_local), c_loc(vprco_grid_local), &
       c_loc(racau_grid_local), c_loc(cnt_grid_local(1)), c_loc(cdnumc_grid_local), c_loc(efcout_grid_local), &
       c_loc(efiout_grid_local), c_loc(ncout_grid_local), c_loc(niout_grid_local), c_loc(freql_grid_local), c_loc(freqi_grid_local), &
       c_loc(icwmrst_grid_out_local), c_loc(icimrst_grid_out_local), c_loc(fcti_grid_local), c_loc(fctl_grid_local), &
       c_loc(ctrel_grid_local), c_loc(ctrei_grid_local), c_loc(ctnl_grid_local), c_loc(ctni_grid_local), &
       c_loc(evprain_st_grid_local))

end subroutine micro_mg_cam_grid_diag

subroutine micro_mg_cam_grid_diag_native(ngrdcol_local, minlwp_local, iclwpst_grid_local, cld_grid_local, cmeliq_grid_local, &
     pdel_grid_local, prec_str_grid_local, acgcme_grid_local, acprecl_grid_local, acnum_grid_local, prao_grid_local, &
     prco_grid_local, nc_grid_local, liqcldf_grid_local, icwmrst_grid_local, rel_grid_local, icwnc_grid_local, &
     icecldf_grid_local, icimrst_grid_local, rei_grid_local, icinc_grid_local, nevapr_grid_local, evpsnow_st_grid_local, &
     tgliqwp_grid_local, tgcmeliq_grid_local, pe_grid_local, tpr_grid_local, pefrac_grid_local, vprao_grid_local, &
     vprco_grid_local, racau_grid_local, cnt_grid_local, cdnumc_grid_local, efcout_grid_local, efiout_grid_local, &
     ncout_grid_local, niout_grid_local, freql_grid_local, freqi_grid_local, icwmrst_grid_out_local, &
     icimrst_grid_out_local, fcti_grid_local, fctl_grid_local, ctrel_grid_local, ctrei_grid_local, ctnl_grid_local, &
     ctni_grid_local, evprain_st_grid_local)

  use ref_pres, only: top_lev => trop_cloud_top_lev

  integer, intent(in) :: ngrdcol_local
  real(r8), intent(in) :: minlwp_local
  real(r8), intent(in) :: iclwpst_grid_local(pcols,pver), cld_grid_local(pcols,pver)
  real(r8), intent(in) :: cmeliq_grid_local(pcols,pver), pdel_grid_local(pcols,pver)
  real(r8), intent(in) :: prec_str_grid_local(pcols), prao_grid_local(pcols,pver), prco_grid_local(pcols,pver)
  real(r8), intent(in) :: nc_grid_local(pcols,pver), liqcldf_grid_local(pcols,pver), icwmrst_grid_local(pcols,pver)
  real(r8), intent(in) :: rel_grid_local(pcols,pver), icwnc_grid_local(pcols,pver), icecldf_grid_local(pcols,pver)
  real(r8), intent(in) :: icimrst_grid_local(pcols,pver), rei_grid_local(pcols,pver), icinc_grid_local(pcols,pver)
  real(r8), intent(in) :: nevapr_grid_local(pcols,pver)
  real(r8), intent(inout) :: evpsnow_st_grid_local(pcols,pver)
  real(r8), intent(inout) :: acgcme_grid_local(:), acprecl_grid_local(:)
  integer, intent(inout) :: acnum_grid_local(:)
  real(r8), intent(inout) :: tgliqwp_grid_local(pcols), tgcmeliq_grid_local(pcols), pe_grid_local(pcols)
  real(r8), intent(inout) :: tpr_grid_local(pcols), pefrac_grid_local(pcols), vprao_grid_local(pcols)
  real(r8), intent(inout) :: vprco_grid_local(pcols), racau_grid_local(pcols), cdnumc_grid_local(pcols)
  integer, intent(inout) :: cnt_grid_local(pcols)
  real(r8), intent(inout) :: efcout_grid_local(pcols,pver), efiout_grid_local(pcols,pver)
  real(r8), intent(inout) :: ncout_grid_local(pcols,pver), niout_grid_local(pcols,pver)
  real(r8), intent(inout) :: freql_grid_local(pcols,pver), freqi_grid_local(pcols,pver)
  real(r8), intent(inout) :: icwmrst_grid_out_local(pcols,pver), icimrst_grid_out_local(pcols,pver)
  real(r8), intent(inout) :: fcti_grid_local(pcols), fctl_grid_local(pcols), ctrel_grid_local(pcols)
  real(r8), intent(inout) :: ctrei_grid_local(pcols), ctnl_grid_local(pcols), ctni_grid_local(pcols)
  real(r8), intent(inout) :: evprain_st_grid_local(pcols,pver)
  integer :: i, k

  tgliqwp_grid_local(:ngrdcol_local) = 0._r8
  tgcmeliq_grid_local(:ngrdcol_local) = 0._r8
  do k = top_lev, pver
     do i = 1, ngrdcol_local
        tgliqwp_grid_local(i) = tgliqwp_grid_local(i) + iclwpst_grid_local(i,k)*cld_grid_local(i,k)
        if (cmeliq_grid_local(i,k) > 1.e-12_r8) then
           tgcmeliq_grid_local(i) = tgcmeliq_grid_local(i) + cmeliq_grid_local(i,k) * &
                (pdel_grid_local(i,k) / gravit) / rhoh2o
        end if
     end do
  end do

  pe_grid_local(:ngrdcol_local) = 0._r8
  tpr_grid_local(:ngrdcol_local) = 0._r8
  pefrac_grid_local(:ngrdcol_local) = 0._r8

  do i = 1, ngrdcol_local
     acgcme_grid_local(i) = acgcme_grid_local(i) + tgcmeliq_grid_local(i)
     acprecl_grid_local(i) = acprecl_grid_local(i) + prec_str_grid_local(i)
     acnum_grid_local(i) = acnum_grid_local(i) + 1

     if (tgliqwp_grid_local(i) < minlwp_local) then
        if (acprecl_grid_local(i) > 5.e-8_r8) then
           tpr_grid_local(i) = max(acprecl_grid_local(i)/acnum_grid_local(i), 1.e-15_r8)
           if (acgcme_grid_local(i) > 1.e-10_r8) then
              pe_grid_local(i) = min(max(acprecl_grid_local(i)/acgcme_grid_local(i), 1.e-15_r8), 1.e5_r8)
              pefrac_grid_local(i) = 1._r8
           end if
        end if

        acprecl_grid_local(i) = 0._r8
        acgcme_grid_local(i)  = 0._r8
        acnum_grid_local(i)   = 0
     end if

     if (acnum_grid_local(i) > 1000) then
        acnum_grid_local(i)   = 0
        acprecl_grid_local(i) = 0._r8
        acgcme_grid_local(i)  = 0._r8
     end if
  end do

  vprao_grid_local = 0._r8
  cnt_grid_local = 0
  do k = top_lev, pver
     vprao_grid_local(:ngrdcol_local) = vprao_grid_local(:ngrdcol_local) + prao_grid_local(:ngrdcol_local,k)
     where (prao_grid_local(:ngrdcol_local,k) /= 0._r8) cnt_grid_local(:ngrdcol_local) = cnt_grid_local(:ngrdcol_local) + 1
  end do
  where (cnt_grid_local > 0) vprao_grid_local = vprao_grid_local/cnt_grid_local

  vprco_grid_local = 0._r8
  cnt_grid_local = 0
  do k = top_lev, pver
     vprco_grid_local(:ngrdcol_local) = vprco_grid_local(:ngrdcol_local) + prco_grid_local(:ngrdcol_local,k)
     where (prco_grid_local(:ngrdcol_local,k) /= 0._r8) cnt_grid_local(:ngrdcol_local) = cnt_grid_local(:ngrdcol_local) + 1
  end do
  where (cnt_grid_local > 0)
     vprco_grid_local = vprco_grid_local/cnt_grid_local
     racau_grid_local = vprao_grid_local/vprco_grid_local
  elsewhere
     racau_grid_local = 0._r8
  end where
  racau_grid_local = min(racau_grid_local, 1.e10_r8)

  cdnumc_grid_local(:ngrdcol_local) = 0._r8
  do k = top_lev, pver
     do i = 1, ngrdcol_local
        cdnumc_grid_local(i) = cdnumc_grid_local(i) + nc_grid_local(i,k) * pdel_grid_local(i,k) / gravit
     end do
  end do

  efcout_grid_local = 0._r8
  efiout_grid_local = 0._r8
  ncout_grid_local = 0._r8
  niout_grid_local = 0._r8
  freql_grid_local = 0._r8
  freqi_grid_local = 0._r8
  icwmrst_grid_out_local = 0._r8
  icimrst_grid_out_local = 0._r8

  do k = top_lev, pver
     do i = 1, ngrdcol_local
        if (liqcldf_grid_local(i,k) > 0.01_r8 .and. icwmrst_grid_local(i,k) > 5.e-5_r8) then
           efcout_grid_local(i,k) = rel_grid_local(i,k) * liqcldf_grid_local(i,k)
           ncout_grid_local(i,k) = icwnc_grid_local(i,k) * liqcldf_grid_local(i,k)
           freql_grid_local(i,k) = liqcldf_grid_local(i,k)
           icwmrst_grid_out_local(i,k) = icwmrst_grid_local(i,k)
        end if
        if (icecldf_grid_local(i,k) > 0.01_r8 .and. icimrst_grid_local(i,k) > 1.e-6_r8) then
           efiout_grid_local(i,k) = rei_grid_local(i,k) * icecldf_grid_local(i,k)
           niout_grid_local(i,k) = icinc_grid_local(i,k) * icecldf_grid_local(i,k)
           freqi_grid_local(i,k) = icecldf_grid_local(i,k)
           icimrst_grid_out_local(i,k) = icimrst_grid_local(i,k)
        end if
     end do
  end do

  fcti_grid_local = 0._r8
  fctl_grid_local = 0._r8
  ctrel_grid_local = 0._r8
  ctrei_grid_local = 0._r8
  ctnl_grid_local = 0._r8
  ctni_grid_local = 0._r8
  do i = 1, ngrdcol_local
     do k = top_lev, pver
        if (liqcldf_grid_local(i,k) > 0.01_r8 .and. icwmrst_grid_local(i,k) > 1.e-7_r8) then
           ctrel_grid_local(i) = rel_grid_local(i,k) * liqcldf_grid_local(i,k)
           ctnl_grid_local(i) = icwnc_grid_local(i,k) * liqcldf_grid_local(i,k)
           fctl_grid_local(i) = liqcldf_grid_local(i,k)
           exit
        end if
        if (icecldf_grid_local(i,k) > 0.01_r8 .and. icimrst_grid_local(i,k) > 1.e-7_r8) then
           ctrei_grid_local(i) = rei_grid_local(i,k) * icecldf_grid_local(i,k)
           ctni_grid_local(i) = icinc_grid_local(i,k) * icecldf_grid_local(i,k)
           fcti_grid_local(i) = icecldf_grid_local(i,k)
           exit
        end if
     end do
  end do

  evprain_st_grid_local(:ngrdcol_local,:pver) = nevapr_grid_local(:ngrdcol_local,:pver) - evpsnow_st_grid_local(:ngrdcol_local,:pver)
  do k = top_lev, pver
     do i = 1, ngrdcol_local
        evprain_st_grid_local(i,k) = max(evprain_st_grid_local(i,k), 0._r8)
        evpsnow_st_grid_local(i,k) = max(evpsnow_st_grid_local(i,k), 0._r8)
     end do
  end do

end subroutine micro_mg_cam_grid_diag_native

function p1(tin) result(pout)
  real(r8), target, intent(in) :: tin(:)
  real(r8), pointer :: pout(:)
  pout => tin
end function p1

function p2(tin) result(pout)
  real(r8), target, intent(in) :: tin(:,:)
  real(r8), pointer :: pout(:,:)
  pout => tin
end function p2

end module micro_mg_cam
