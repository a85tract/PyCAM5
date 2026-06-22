!-------------------------------------------------------------------------------
!physics data types module
!-------------------------------------------------------------------------------
module physics_types

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use ppgrid,           only: pcols, pver, psubcols
  use constituents,     only: pcnst, qmin, cnst_name
  use geopotential,     only: geopotential_dse
  use physconst,        only: zvir, gravit, cpair, rair, cpairv, rairv
  use dycore,           only: dycore_is
  use phys_grid,        only: get_ncols_p, get_rlon_all_p, get_rlat_all_p, get_gcol_all_p
  use cam_logfile,      only: iulog
  use cam_abortutils,   only: endrun
  use spmd_utils,       only: masterproc
  use phys_control,     only: waccmx_is
  use shr_const_mod,    only: shr_const_rwv

  implicit none
  private          ! Make default type private to the module

  logical, parameter :: adjust_te = .FALSE.

! Public types:

  public physics_state
  public physics_tend
  public physics_ptend

! Public interfaces

  public physics_update
  public physics_state_check ! Check state object for invalid data.
  public physics_ptend_reset
  public physics_ptend_init
  public physics_state_set_grid
  public physics_dme_adjust  ! adjust dry mass and energy for change in water
                             ! cannot be applied to eul or sld dycores
  public physics_state_copy  ! copy a physics_state object
  public physics_ptend_copy  ! copy a physics_ptend object
  public physics_ptend_sum   ! accumulate physics_ptend objects
  public physics_ptend_scale ! Multiply physics_ptend objects by a constant factor.
  public physics_tend_init   ! initialize a physics_tend object

  public set_state_pdry      ! calculate dry air masses in state variable
  public set_wet_to_dry
  public set_dry_to_wet
  public physics_type_alloc

  public physics_state_alloc   ! allocate individual components within state
  public physics_state_dealloc ! deallocate individual components within state
  public physics_tend_alloc    ! allocate individual components within tend
  public physics_tend_dealloc  ! deallocate individual components within tend
  public physics_ptend_alloc   ! allocate individual components within tend
  public physics_ptend_dealloc ! deallocate individual components within tend

!-------------------------------------------------------------------------------
  type physics_state
     integer                                     :: &
          lchnk,                &! chunk index
          ngrdcol,              &! -- Grid        -- number of active columns (on the grid)
          psetcols=0,           &! --             -- max number of columns set - if subcols = pcols*psubcols, else = pcols
          ncol=0                 ! --             -- sum of nsubcol for all ngrdcols - number of active columns
     real(r8), dimension(:), allocatable         :: &
          lat,     &! latitude (radians)
          lon,     &! longitude (radians)
          ps,      &! surface pressure
          psdry,   &! dry surface pressure
          phis,    &! surface geopotential
          ulat,    &! unique latitudes  (radians)
          ulon      ! unique longitudes (radians)
     real(r8), dimension(:,:),allocatable        :: &
          t,       &! temperature (K)
          u,       &! zonal wind (m/s)
          v,       &! meridional wind (m/s)
          s,       &! dry static energy
          omega,   &! vertical pressure velocity (Pa/s)
          pmid,    &! midpoint pressure (Pa)
          pmiddry, &! midpoint pressure dry (Pa)
          pdel,    &! layer thickness (Pa)
          pdeldry, &! layer thickness dry (Pa)
          rpdel,   &! reciprocal of layer thickness (Pa)
          rpdeldry,&! recipricol layer thickness dry (Pa)
          lnpmid,  &! ln(pmid)
          lnpmiddry,&! log midpoint pressure dry (Pa)
          exner,   &! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
          zm        ! geopotential height above surface at midpoints (m)

     real(r8), dimension(:,:,:),allocatable      :: &
          q         ! constituent mixing ratio (kg/kg moist or dry air depending on type)

     real(r8), dimension(:,:),allocatable        :: &
          pint,    &! interface pressure (Pa)
          pintdry, &! interface pressure dry (Pa)
          lnpint,  &! ln(pint)
          lnpintdry,&! log interface pressure dry (Pa)
          zi        ! geopotential height above surface at interfaces (m)

     real(r8), dimension(:),allocatable          :: &
          te_ini,  &! vertically integrated total (kinetic + static) energy of initial state
          te_cur,  &! vertically integrated total (kinetic + static) energy of current state
          tw_ini,  &! vertically integrated total water of initial state
          tw_cur    ! vertically integrated total water of new state
     integer :: count ! count of values with significant energy or water imbalances
     integer, dimension(:),allocatable           :: &
          latmapback, &! map from column to unique lat for that column
          lonmapback, &! map from column to unique lon for that column
          cid        ! unique column id
     integer :: ulatcnt, &! number of unique lats in chunk
                uloncnt   ! number of unique lons in chunk

  end type physics_state

!-------------------------------------------------------------------------------
  type physics_tend

     integer   ::   psetcols=0 ! max number of columns set- if subcols = pcols*psubcols, else = pcols

     real(r8), dimension(:,:),allocatable        :: dtdt, dudt, dvdt
     real(r8), dimension(:),  allocatable        :: flx_net
     real(r8), dimension(:),  allocatable        :: &
          te_tnd,  &! cumulative boundary flux of total energy
          tw_tnd    ! cumulative boundary flux of total water
  end type physics_tend

!-------------------------------------------------------------------------------
! This is for tendencies returned from individual parameterizations
  type physics_ptend

     integer   ::   psetcols=0 ! max number of columns set- if subcols = pcols*psubcols, else = pcols

     character*24 :: name    ! name of parameterization which produced tendencies.

     logical ::             &
          ls = .false.,               &! true if dsdt is returned
          lu = .false.,               &! true if dudt is returned
          lv = .false.                 ! true if dvdt is returned

     logical,dimension(pcnst) ::  lq = .false.  ! true if dqdt() is returned

     integer ::             &
          top_level,        &! top level index for which nonzero tendencies have been set
          bot_level          ! bottom level index for which nonzero tendencies have been set

     real(r8), dimension(:,:),allocatable   :: &
          s,                &! heating rate (J/kg/s)
          u,                &! u momentum tendency (m/s/s)
          v                  ! v momentum tendency (m/s/s)
     real(r8), dimension(:,:,:),allocatable :: &
          q                  ! consituent tendencies (kg/kg/s)

! boundary fluxes
     real(r8), dimension(:),allocatable     ::&
          hflux_srf,     &! net heat flux at surface (W/m2)
          hflux_top,     &! net heat flux at top of model (W/m2)
          taux_srf,      &! net zonal stress at surface (Pa)
          taux_top,      &! net zonal stress at top of model (Pa)
          tauy_srf,      &! net meridional stress at surface (Pa)
          tauy_top        ! net meridional stress at top of model (Pa)
     real(r8), dimension(:,:),allocatable   ::&
          cflx_srf,      &! constituent flux at surface (kg/m2/s)
          cflx_top        ! constituent flux top of model (kg/m2/s)

  end type physics_ptend

  logical :: use_native_zero_impl = .false.
  logical :: zero_impl_selected = .false.
  logical :: zero_proof_written = .false.
  logical :: physics_ptend_reset_logged = .false.
  logical :: physics_ptend_scale_logged = .false.
  logical :: physics_ptend_sum_logged = .false.
  logical :: set_state_pdry_logged = .false.
  logical :: set_wet_to_dry_logged = .false.
  logical :: set_dry_to_wet_logged = .false.
  logical :: init_geo_unique_logged = .false.
  logical :: physics_state_set_grid_logged = .false.
  logical :: physics_state_copy_logged = .false.
  logical :: physics_state_alloc_logged = .false.
  logical :: physics_tend_alloc_logged = .false.
  logical :: state_cnst_min_nz_logged = .false.
  logical :: physics_dme_adjust_logged = .false.
  logical :: physics_update_logged = .false.
  logical :: physics_update_selector_selected = .false.
  logical :: physics_update_native_field = .false.
  logical :: physics_update_native_q = .false.
  logical :: physics_update_native_s = .false.
  logical :: physics_type_alloc_logged = .false.
  logical :: physics_ptend_init_logged = .false.
  logical :: physics_ptend_alloc_logged = .false.
  logical :: physics_ptend_dealloc_logged = .false.
  logical :: physics_state_dealloc_logged = .false.
  logical :: multiply_tendency_logged = .false.

  interface
     function physics_dme_adjust_codon(is_lr_c) result(active_c) &
          bind(c, name="physics_dme_adjust_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: is_lr_c
       integer(c_int64_t) :: active_c
     end function physics_dme_adjust_codon

     function physics_types_touch_codon(stage_c) result(stage_out) bind(c, name="physics_types_touch_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: stage_out
     end function physics_types_touch_codon

     function physics_type_alloc_codon_raw(stage_c) result(stage_out) bind(c, name="physics_type_alloc_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: stage_out
     end function physics_type_alloc_codon_raw

     function physics_state_dealloc_codon_raw(stage_c) result(stage_out) bind(c, name="physics_state_dealloc_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: stage_out
     end function physics_state_dealloc_codon_raw

     function physics_ptend_alloc_codon_raw(stage_c) result(stage_out) bind(c, name="physics_ptend_alloc_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: stage_out
     end function physics_ptend_alloc_codon_raw

     function physics_ptend_dealloc_codon_raw(stage_c) result(stage_out) bind(c, name="physics_ptend_dealloc_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: stage_c
       integer(c_int64_t) :: stage_out
     end function physics_ptend_dealloc_codon_raw

     subroutine physics_tend_init_codon_raw(psetcols_c, pver_c, dtdt_p, dudt_p, dvdt_p, flx_net_p, te_tnd_p, tw_tnd_p) &
          bind(c, name="physics_tend_init_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c
       type(c_ptr), value :: dtdt_p, dudt_p, dvdt_p, flx_net_p, te_tnd_p, tw_tnd_p
     end subroutine physics_tend_init_codon_raw

     subroutine physics_state_alloc_codon_raw(psetcols_c, pver_c, pcnst_c, value_c, &
          lat_p, lon_p, ulat_p, ulon_p, ps_p, psdry_p, phis_p, &
          t_p, u_p, v_p, s_p, omega_p, pmid_p, pmiddry_p, pdel_p, pdeldry_p, rpdel_p, rpdeldry_p, &
          lnpmid_p, lnpmiddry_p, exner_p, zm_p, q_p, &
          pint_p, pintdry_p, lnpint_p, lnpintdry_p, zi_p, &
          te_ini_p, te_cur_p, tw_ini_p, tw_cur_p) bind(c, name="physics_state_alloc_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c, pcnst_c
       real(c_double), value :: value_c
       type(c_ptr), value :: lat_p, lon_p, ulat_p, ulon_p, ps_p, psdry_p, phis_p
       type(c_ptr), value :: t_p, u_p, v_p, s_p, omega_p, pmid_p, pmiddry_p, pdel_p, pdeldry_p, rpdel_p, rpdeldry_p
       type(c_ptr), value :: lnpmid_p, lnpmiddry_p, exner_p, zm_p, q_p
       type(c_ptr), value :: pint_p, pintdry_p, lnpint_p, lnpintdry_p, zi_p
       type(c_ptr), value :: te_ini_p, te_cur_p, tw_ini_p, tw_cur_p
     end subroutine physics_state_alloc_codon_raw

     subroutine physics_tend_alloc_codon_raw(psetcols_c, pver_c, value_c, dtdt_p, dudt_p, dvdt_p, &
          flx_net_p, te_tnd_p, tw_tnd_p) bind(c, name="physics_tend_alloc_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c
       real(c_double), value :: value_c
       type(c_ptr), value :: dtdt_p, dudt_p, dvdt_p, flx_net_p, te_tnd_p, tw_tnd_p
     end subroutine physics_tend_alloc_codon_raw

     function physics_ptend_init_codon_raw(ls_present_c, ls_value_c, lu_present_c, lu_value_c, &
          lv_present_c, lv_value_c, lq_present_c) result(policy_c) bind(c, name="physics_ptend_init_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: ls_present_c, ls_value_c, lu_present_c, lu_value_c
       integer(c_int64_t), value :: lv_present_c, lv_value_c, lq_present_c
       integer(c_int64_t) :: policy_c
     end function physics_ptend_init_codon_raw

     subroutine physics_fill_real_1d_codon_raw(n_c, value_c, arr_p) bind(c, name="physics_fill_real_1d_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: n_c
       real(c_double), value :: value_c
       type(c_ptr), value :: arr_p
     end subroutine physics_fill_real_1d_codon_raw

     subroutine physics_fill_real_2d_codon_raw(ld1_c, nlev_c, value_c, arr_p) bind(c, name="physics_fill_real_2d_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ld1_c, nlev_c
       real(c_double), value :: value_c
       type(c_ptr), value :: arr_p
     end subroutine physics_fill_real_2d_codon_raw

     subroutine physics_fill_real_3d_codon_raw(ld1_c, nlev_c, n3_c, value_c, arr_p) bind(c, name="physics_fill_real_3d_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ld1_c, nlev_c, n3_c
       real(c_double), value :: value_c
       type(c_ptr), value :: arr_p
     end subroutine physics_fill_real_3d_codon_raw

     subroutine physics_ptend_reset_s_codon_raw(psetcols_c, pver_c, s_p, hflux_srf_p, hflux_top_p) &
          bind(c, name="physics_ptend_reset_s_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c
       type(c_ptr), value :: s_p, hflux_srf_p, hflux_top_p
     end subroutine physics_ptend_reset_s_codon_raw

     subroutine physics_ptend_reset_u_codon_raw(psetcols_c, pver_c, u_p, taux_srf_p, taux_top_p) &
          bind(c, name="physics_ptend_reset_u_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c
       type(c_ptr), value :: u_p, taux_srf_p, taux_top_p
     end subroutine physics_ptend_reset_u_codon_raw

     subroutine physics_ptend_reset_v_codon_raw(psetcols_c, pver_c, v_p, tauy_srf_p, tauy_top_p) &
          bind(c, name="physics_ptend_reset_v_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c
       type(c_ptr), value :: v_p, tauy_srf_p, tauy_top_p
     end subroutine physics_ptend_reset_v_codon_raw

     subroutine physics_ptend_reset_q_codon_raw(psetcols_c, pver_c, pcnst_c, q_p, cflx_srf_p, cflx_top_p) &
          bind(c, name="physics_ptend_reset_q_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c, pcnst_c
       type(c_ptr), value :: q_p, cflx_srf_p, cflx_top_p
     end subroutine physics_ptend_reset_q_codon_raw

     subroutine physics_ptend_copy_s_codon_raw(psetcols_c, pver_c, src_s_p, src_hflux_srf_p, src_hflux_top_p, &
          dst_s_p, dst_hflux_srf_p, dst_hflux_top_p) bind(c, name="physics_ptend_copy_s_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c
       type(c_ptr), value :: src_s_p, src_hflux_srf_p, src_hflux_top_p, dst_s_p, dst_hflux_srf_p, dst_hflux_top_p
     end subroutine physics_ptend_copy_s_codon_raw

     subroutine physics_ptend_copy_u_codon_raw(psetcols_c, pver_c, src_u_p, src_taux_srf_p, src_taux_top_p, &
          dst_u_p, dst_taux_srf_p, dst_taux_top_p) bind(c, name="physics_ptend_copy_u_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c
       type(c_ptr), value :: src_u_p, src_taux_srf_p, src_taux_top_p, dst_u_p, dst_taux_srf_p, dst_taux_top_p
     end subroutine physics_ptend_copy_u_codon_raw

     subroutine physics_ptend_copy_v_codon_raw(psetcols_c, pver_c, src_v_p, src_tauy_srf_p, src_tauy_top_p, &
          dst_v_p, dst_tauy_srf_p, dst_tauy_top_p) bind(c, name="physics_ptend_copy_v_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c
       type(c_ptr), value :: src_v_p, src_tauy_srf_p, src_tauy_top_p, dst_v_p, dst_tauy_srf_p, dst_tauy_top_p
     end subroutine physics_ptend_copy_v_codon_raw

     subroutine physics_ptend_copy_q_codon_raw(psetcols_c, pver_c, pcnst_c, src_q_p, src_cflx_srf_p, src_cflx_top_p, &
          dst_q_p, dst_cflx_srf_p, dst_cflx_top_p) bind(c, name="physics_ptend_copy_q_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c, pcnst_c
       type(c_ptr), value :: src_q_p, src_cflx_srf_p, src_cflx_top_p, dst_q_p, dst_cflx_srf_p, dst_cflx_top_p
     end subroutine physics_ptend_copy_q_codon_raw

     subroutine physics_ptend_scale_field_codon_raw(ncol_c, psetcols_c, top_level_c, bot_level_c, fac_c, field_p, &
          flx_srf_p, flx_top_p) bind(c, name="physics_ptend_scale_field_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, top_level_c, bot_level_c
       real(c_double), value :: fac_c
       type(c_ptr), value :: field_p, flx_srf_p, flx_top_p
     end subroutine physics_ptend_scale_field_codon_raw

     subroutine physics_ptend_sum_field_codon_raw(ncol_c, psetcols_c, top_level_c, bot_level_c, src_field_p, &
          dst_field_p, src_flx_srf_p, dst_flx_srf_p, src_flx_top_p, dst_flx_top_p) &
          bind(c, name="physics_ptend_sum_field_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, top_level_c, bot_level_c
       type(c_ptr), value :: src_field_p, dst_field_p, src_flx_srf_p, dst_flx_srf_p
       type(c_ptr), value :: src_flx_top_p, dst_flx_top_p
     end subroutine physics_ptend_sum_field_codon_raw

     subroutine physics_ptend_reset_codon_raw(psetcols_c, pver_c, pcnst_c, ls_c, lu_c, lv_c, lq_p, &
          s_p, u_p, v_p, q_p, hflux_srf_p, hflux_top_p, taux_srf_p, taux_top_p, tauy_srf_p, tauy_top_p, &
          cflx_srf_p, cflx_top_p) bind(c, name="physics_ptend_reset_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: psetcols_c, pver_c, pcnst_c, ls_c, lu_c, lv_c
       type(c_ptr), value :: lq_p, s_p, u_p, v_p, q_p
       type(c_ptr), value :: hflux_srf_p, hflux_top_p, taux_srf_p, taux_top_p
       type(c_ptr), value :: tauy_srf_p, tauy_top_p, cflx_srf_p, cflx_top_p
     end subroutine physics_ptend_reset_codon_raw

     subroutine physics_ptend_scale_codon_raw(ncol_c, psetcols_c, pver_c, pcnst_c, top_level_c, &
          bot_level_c, fac_c, ls_c, lu_c, lv_c, lq_p, s_p, u_p, v_p, q_p, hflux_srf_p, hflux_top_p, &
          taux_srf_p, taux_top_p, tauy_srf_p, tauy_top_p, cflx_srf_p, cflx_top_p) &
          bind(c, name="physics_ptend_scale_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, pcnst_c, top_level_c, bot_level_c
       integer(c_int64_t), value :: ls_c, lu_c, lv_c
       real(c_double), value :: fac_c
       type(c_ptr), value :: lq_p, s_p, u_p, v_p, q_p
       type(c_ptr), value :: hflux_srf_p, hflux_top_p, taux_srf_p, taux_top_p
       type(c_ptr), value :: tauy_srf_p, tauy_top_p, cflx_srf_p, cflx_top_p
     end subroutine physics_ptend_scale_codon_raw

     subroutine physics_ptend_sum_codon_raw(ncol_c, psetcols_c, pver_c, pcnst_c, top_level_c, bot_level_c, &
          ls_c, lu_c, lv_c, lq_p, src_s_p, dst_s_p, src_u_p, dst_u_p, src_v_p, dst_v_p, src_q_p, dst_q_p, &
          src_hflux_srf_p, dst_hflux_srf_p, src_hflux_top_p, dst_hflux_top_p, src_taux_srf_p, &
          dst_taux_srf_p, src_taux_top_p, dst_taux_top_p, src_tauy_srf_p, dst_tauy_srf_p, &
          src_tauy_top_p, dst_tauy_top_p, src_cflx_srf_p, dst_cflx_srf_p, src_cflx_top_p, &
          dst_cflx_top_p) bind(c, name="physics_ptend_sum_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, pcnst_c, top_level_c, bot_level_c
       integer(c_int64_t), value :: ls_c, lu_c, lv_c
       type(c_ptr), value :: lq_p, src_s_p, dst_s_p, src_u_p, dst_u_p, src_v_p, dst_v_p, src_q_p, dst_q_p
       type(c_ptr), value :: src_hflux_srf_p, dst_hflux_srf_p, src_hflux_top_p, dst_hflux_top_p
       type(c_ptr), value :: src_taux_srf_p, dst_taux_srf_p, src_taux_top_p, dst_taux_top_p
       type(c_ptr), value :: src_tauy_srf_p, dst_tauy_srf_p, src_tauy_top_p, dst_tauy_top_p
       type(c_ptr), value :: src_cflx_srf_p, dst_cflx_srf_p, src_cflx_top_p, dst_cflx_top_p
     end subroutine physics_ptend_sum_codon_raw

     subroutine set_state_pdry_codon_raw(ncol_c, psetcols_c, pver_c, do_pdeld_calc_c, psdry_p, pint_p, &
          pdel_p, q_p, pdeldry_p, pintdry_p, pmiddry_p, rpdeldry_p, lnpmiddry_p, lnpintdry_p) &
          bind(c, name="set_state_pdry_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, do_pdeld_calc_c
       type(c_ptr), value :: psdry_p, pint_p, pdel_p, q_p, pdeldry_p, pintdry_p
       type(c_ptr), value :: pmiddry_p, rpdeldry_p, lnpmiddry_p, lnpintdry_p
     end subroutine set_state_pdry_codon_raw

     subroutine physics_set_wet_to_dry_constituent_codon_raw(ncol_c, psetcols_c, pver_c, q_p, pdel_p, &
          pdeldry_p) bind(c, name="physics_set_wet_to_dry_constituent_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c
       type(c_ptr), value :: q_p, pdel_p, pdeldry_p
     end subroutine physics_set_wet_to_dry_constituent_codon_raw

     subroutine physics_set_dry_to_wet_constituent_codon_raw(ncol_c, psetcols_c, pver_c, q_p, pdeldry_p, &
          pdel_p) bind(c, name="physics_set_dry_to_wet_constituent_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c
       type(c_ptr), value :: q_p, pdeldry_p, pdel_p
     end subroutine physics_set_dry_to_wet_constituent_codon_raw

     subroutine set_wet_to_dry_codon_raw(ncol_c, psetcols_c, pver_c, pcnst_c, dry_mask_p, q_p, pdel_p, &
          pdeldry_p) bind(c, name="set_wet_to_dry_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, pcnst_c
       type(c_ptr), value :: dry_mask_p, q_p, pdel_p, pdeldry_p
     end subroutine set_wet_to_dry_codon_raw

     subroutine set_dry_to_wet_codon_raw(ncol_c, psetcols_c, pver_c, pcnst_c, dry_mask_p, q_p, pdeldry_p, &
          pdel_p) bind(c, name="set_dry_to_wet_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, pcnst_c
       type(c_ptr), value :: dry_mask_p, q_p, pdeldry_p, pdel_p
     end subroutine set_dry_to_wet_codon_raw

     subroutine init_geo_unique_codon_raw(ncol_c, psetcols_c, lat_p, lon_p, ulat_p, ulon_p, &
          latmapback_p, lonmapback_p, ulatcnt_p, uloncnt_p) bind(c, name="init_geo_unique_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c
       type(c_ptr), value :: lat_p, lon_p, ulat_p, ulon_p, latmapback_p, lonmapback_p, ulatcnt_p, uloncnt_p
     end subroutine init_geo_unique_codon_raw

     subroutine physics_state_set_grid_codon_raw(ncol_c, psetcols_c, rlat_p, rlon_p, lat_p, lon_p, ulat_p, ulon_p, &
          latmapback_p, lonmapback_p, ulatcnt_p, uloncnt_p) bind(c, name="physics_state_set_grid_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c
       type(c_ptr), value :: rlat_p, rlon_p, lat_p, lon_p, ulat_p, ulon_p
       type(c_ptr), value :: latmapback_p, lonmapback_p, ulatcnt_p, uloncnt_p
     end subroutine physics_state_set_grid_codon_raw

     subroutine physics_copy_real_1d_codon_raw(n_c, src_p, dst_p) bind(c, name="physics_copy_real_1d_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: n_c
       type(c_ptr), value :: src_p, dst_p
     end subroutine physics_copy_real_1d_codon_raw

     subroutine physics_copy_real_2d_codon_raw(ncol_c, ld1_c, nlev_c, src_p, dst_p) &
          bind(c, name="physics_copy_real_2d_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, ld1_c, nlev_c
       type(c_ptr), value :: src_p, dst_p
     end subroutine physics_copy_real_2d_codon_raw

     subroutine physics_copy_real_3d_codon_raw(ncol_c, ld1_c, nlev_c, n3_c, src_p, dst_p) &
          bind(c, name="physics_copy_real_3d_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, ld1_c, nlev_c, n3_c
       type(c_ptr), value :: src_p, dst_p
     end subroutine physics_copy_real_3d_codon_raw

     subroutine physics_state_copy_codon_raw(ncol_c, psetcols_c, pver_c, pcnst_c, &
          src_lat_p, dst_lat_p, src_lon_p, dst_lon_p, src_ps_p, dst_ps_p, src_phis_p, dst_phis_p, &
          src_te_ini_p, dst_te_ini_p, src_te_cur_p, dst_te_cur_p, src_tw_ini_p, dst_tw_ini_p, &
          src_tw_cur_p, dst_tw_cur_p, src_psdry_p, dst_psdry_p, &
          src_t_p, dst_t_p, src_u_p, dst_u_p, src_v_p, dst_v_p, src_s_p, dst_s_p, &
          src_omega_p, dst_omega_p, src_pmid_p, dst_pmid_p, src_pdel_p, dst_pdel_p, &
          src_rpdel_p, dst_rpdel_p, src_lnpmid_p, dst_lnpmid_p, src_exner_p, dst_exner_p, &
          src_zm_p, dst_zm_p, src_lnpmiddry_p, dst_lnpmiddry_p, src_pmiddry_p, dst_pmiddry_p, &
          src_pdeldry_p, dst_pdeldry_p, src_rpdeldry_p, dst_rpdeldry_p, &
          src_pint_p, dst_pint_p, src_lnpint_p, dst_lnpint_p, src_zi_p, dst_zi_p, &
          src_pintdry_p, dst_pintdry_p, src_lnpintdry_p, dst_lnpintdry_p, &
          src_q_p, dst_q_p) bind(c, name="physics_state_copy_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, pcnst_c
       type(c_ptr), value :: src_lat_p, dst_lat_p, src_lon_p, dst_lon_p, src_ps_p, dst_ps_p
       type(c_ptr), value :: src_phis_p, dst_phis_p, src_te_ini_p, dst_te_ini_p, src_te_cur_p, dst_te_cur_p
       type(c_ptr), value :: src_tw_ini_p, dst_tw_ini_p, src_tw_cur_p, dst_tw_cur_p, src_psdry_p, dst_psdry_p
       type(c_ptr), value :: src_t_p, dst_t_p, src_u_p, dst_u_p, src_v_p, dst_v_p, src_s_p, dst_s_p
       type(c_ptr), value :: src_omega_p, dst_omega_p, src_pmid_p, dst_pmid_p, src_pdel_p, dst_pdel_p
       type(c_ptr), value :: src_rpdel_p, dst_rpdel_p, src_lnpmid_p, dst_lnpmid_p, src_exner_p, dst_exner_p
       type(c_ptr), value :: src_zm_p, dst_zm_p, src_lnpmiddry_p, dst_lnpmiddry_p, src_pmiddry_p, dst_pmiddry_p
       type(c_ptr), value :: src_pdeldry_p, dst_pdeldry_p, src_rpdeldry_p, dst_rpdeldry_p
       type(c_ptr), value :: src_pint_p, dst_pint_p, src_lnpint_p, dst_lnpint_p, src_zi_p, dst_zi_p
       type(c_ptr), value :: src_pintdry_p, dst_pintdry_p, src_lnpintdry_p, dst_lnpintdry_p
       type(c_ptr), value :: src_q_p, dst_q_p
     end subroutine physics_state_copy_codon_raw

     subroutine state_cnst_min_nz_codon_raw(ncol_c, psetcols_c, pver_c, pcnst_c, q_p, lim_c, qix_c, numix_c) &
          bind(c, name="state_cnst_min_nz_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, pcnst_c, qix_c, numix_c
       real(c_double), value :: lim_c
       type(c_ptr), value :: q_p
     end subroutine state_cnst_min_nz_codon_raw

     subroutine physics_update_field_codon_raw(ncol_c, psetcols_c, pver_c, top_level_c, bot_level_c, dt_c, &
          update_tend_c, state_field_p, ptend_field_p, tend_field_p) bind(c, name="physics_update_field_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, top_level_c, bot_level_c, update_tend_c
       real(c_double), value :: dt_c
       type(c_ptr), value :: state_field_p, ptend_field_p, tend_field_p
     end subroutine physics_update_field_codon_raw

     subroutine physics_update_q_codon_raw(ncol_c, psetcols_c, pver_c, pcnst_c, top_level_c, bot_level_c, &
          dt_c, m_c, is_number_c, state_q_p, ptend_q_p) bind(c, name="physics_update_q_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, pcnst_c, top_level_c, bot_level_c
       integer(c_int64_t), value :: m_c, is_number_c
       real(c_double), value :: dt_c
       type(c_ptr), value :: state_q_p, ptend_q_p
     end subroutine physics_update_q_codon_raw

     subroutine physics_update_s_codon_raw(ncol_c, psetcols_c, pver_c, top_level_c, bot_level_c, dt_c, &
          update_tend_c, state_s_p, ptend_s_p, tend_dtdt_p, cpairv_loc_p) bind(c, name="physics_update_s_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, top_level_c, bot_level_c, update_tend_c
       real(c_double), value :: dt_c
       type(c_ptr), value :: state_s_p, ptend_s_p, tend_dtdt_p, cpairv_loc_p
     end subroutine physics_update_s_codon_raw
  end interface


!===============================================================================
contains
!===============================================================================
  subroutine physics_types_zero_select_impl()
    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (zero_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('PHYSICS_TYPES_ZERO_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_zero_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_zero_impl = .false.
    end if

    zero_impl_selected = .true.

    if (masterproc) then
       if (use_native_zero_impl) then
          write(iulog,*) 'physics_types_zero implementation = native'
       else
          write(iulog,*) 'physics_types_zero implementation = codon'
       end if
    end if
  end subroutine physics_types_zero_select_impl

  subroutine physics_update_select_impl()
    character(len=96) :: mask, normalized
    integer :: status, n, i, code

    if (physics_update_selector_selected) return

    mask = 'q'
    n = len_trim(mask)
    call get_environment_variable('PHYSICS_UPDATE_NATIVE', value=mask, length=n, status=status)

    if (status /= 0 .or. n <= 0) then
       mask = 'q'
       n = len_trim(mask)
    end if

    normalized = ',' // trim(adjustl(mask(:n))) // ','
    do i = 1, len_trim(normalized)
       code = iachar(normalized(i:i))
       if (code >= iachar('A') .and. code <= iachar('Z')) then
          normalized(i:i) = achar(code + iachar('a') - iachar('A'))
       else if (normalized(i:i) == ' ' .or. normalized(i:i) == ';' .or. &
                normalized(i:i) == ':' .or. normalized(i:i) == '/') then
          normalized(i:i) = ','
       end if
    end do

    if (index(normalized, ',all,') > 0) then
       physics_update_native_field = .true.
       physics_update_native_q = .true.
       physics_update_native_s = .true.
    else
       physics_update_native_field = index(normalized, ',field,') > 0 .or. index(normalized, ',uv,') > 0
       physics_update_native_q = index(normalized, ',q,') > 0
       physics_update_native_s = index(normalized, ',s,') > 0
    end if

    physics_update_selector_selected = .true.

    if (masterproc .and. (physics_update_native_field .or. physics_update_native_q .or. physics_update_native_s)) then
       write(iulog,*) 'physics_update native selector = ', trim(adjustl(mask(:n)))
    end if
  end subroutine physics_update_select_impl

  subroutine physics_types_zero_proof_once()
    if (zero_proof_written) return
    zero_proof_written = .true.
    if (masterproc) then
       write(iulog,'(A)') 'physics_types_zero helpers entered (tend init/ptend reset/copy/scale = codon)'
    end if
  end subroutine physics_types_zero_proof_once

  subroutine physics_types_log_direct(logged, proof_line)
    logical, intent(inout) :: logged
    character(len=*), intent(in) :: proof_line

    if (logged) return
    logged = .true.

    if (masterproc) then
       write(iulog,'(A)') trim(proof_line)
       call flush(iulog)
    end if
  end subroutine physics_types_log_direct

  subroutine physics_tend_init_codon(psetcols_local, dtdt, dudt, dvdt, flx_net, te_tnd, tw_tnd)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(inout) :: dtdt(psetcols_local,pver)
    real(r8), target, intent(inout) :: dudt(psetcols_local,pver)
    real(r8), target, intent(inout) :: dvdt(psetcols_local,pver)
    real(r8), target, intent(inout) :: flx_net(psetcols_local)
    real(r8), target, intent(inout) :: te_tnd(psetcols_local)
    real(r8), target, intent(inout) :: tw_tnd(psetcols_local)

    call physics_tend_init_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
         c_loc(dtdt), c_loc(dudt), c_loc(dvdt), c_loc(flx_net), c_loc(te_tnd), c_loc(tw_tnd))
  end subroutine physics_tend_init_codon

  subroutine physics_state_alloc_codon(state, psetcols_local, value_local)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    type(physics_state), intent(inout), target :: state
    integer, intent(in) :: psetcols_local
    real(r8), intent(in) :: value_local

    call physics_state_alloc_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), int(pcnst, c_int64_t), &
         real(value_local, c_double), &
         c_loc(state%lat), c_loc(state%lon), c_loc(state%ulat), c_loc(state%ulon), c_loc(state%ps), c_loc(state%psdry), &
         c_loc(state%phis), c_loc(state%t), c_loc(state%u), c_loc(state%v), c_loc(state%s), c_loc(state%omega), &
         c_loc(state%pmid), c_loc(state%pmiddry), c_loc(state%pdel), c_loc(state%pdeldry), &
         c_loc(state%rpdel), c_loc(state%rpdeldry), c_loc(state%lnpmid), c_loc(state%lnpmiddry), &
         c_loc(state%exner), c_loc(state%zm), c_loc(state%q), c_loc(state%pint), c_loc(state%pintdry), &
         c_loc(state%lnpint), c_loc(state%lnpintdry), c_loc(state%zi), &
         c_loc(state%te_ini), c_loc(state%te_cur), c_loc(state%tw_ini), c_loc(state%tw_cur))
  end subroutine physics_state_alloc_codon

  subroutine physics_tend_alloc_codon(tend, psetcols_local, value_local)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    type(physics_tend), intent(inout), target :: tend
    integer, intent(in) :: psetcols_local
    real(r8), intent(in) :: value_local

    call physics_tend_alloc_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
         real(value_local, c_double), c_loc(tend%dtdt), c_loc(tend%dudt), c_loc(tend%dvdt), &
         c_loc(tend%flx_net), c_loc(tend%te_tnd), c_loc(tend%tw_tnd))
  end subroutine physics_tend_alloc_codon

  function physics_ptend_init_codon(ls_present, ls_value, lu_present, lu_value, lv_present, lv_value, lq_present) result(policy)
    use iso_c_binding, only: c_int64_t
    logical, intent(in) :: ls_present, ls_value, lu_present, lu_value, lv_present, lv_value, lq_present
    integer :: policy

    policy = int(physics_ptend_init_codon_raw(merge(1_c_int64_t, 0_c_int64_t, ls_present), &
         merge(1_c_int64_t, 0_c_int64_t, ls_value), merge(1_c_int64_t, 0_c_int64_t, lu_present), &
         merge(1_c_int64_t, 0_c_int64_t, lu_value), merge(1_c_int64_t, 0_c_int64_t, lv_present), &
         merge(1_c_int64_t, 0_c_int64_t, lv_value), merge(1_c_int64_t, 0_c_int64_t, lq_present)))
  end function physics_ptend_init_codon

  subroutine physics_fill_real_1d_codon(n_local, value_local, arr)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: n_local
    real(r8), intent(in) :: value_local
    real(r8), target, intent(inout) :: arr(:)

    call physics_fill_real_1d_codon_raw(int(n_local, c_int64_t), real(value_local, c_double), c_loc(arr))
  end subroutine physics_fill_real_1d_codon

  subroutine physics_fill_real_2d_codon(ld1_local, nlev_local, value_local, arr)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: ld1_local, nlev_local
    real(r8), intent(in) :: value_local
    real(r8), target, intent(inout) :: arr(:,:)

    call physics_fill_real_2d_codon_raw(int(ld1_local, c_int64_t), int(nlev_local, c_int64_t), &
         real(value_local, c_double), c_loc(arr))
  end subroutine physics_fill_real_2d_codon

  subroutine physics_fill_real_3d_codon(ld1_local, nlev_local, n3_local, value_local, arr)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: ld1_local, nlev_local, n3_local
    real(r8), intent(in) :: value_local
    real(r8), target, intent(inout) :: arr(:,:,:)

    call physics_fill_real_3d_codon_raw(int(ld1_local, c_int64_t), int(nlev_local, c_int64_t), int(n3_local, c_int64_t), &
         real(value_local, c_double), c_loc(arr))
  end subroutine physics_fill_real_3d_codon

  subroutine physics_ptend_reset_s_codon(psetcols_local, s, hflux_srf, hflux_top)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(inout) :: s(psetcols_local,pver)
    real(r8), target, intent(inout) :: hflux_srf(psetcols_local)
    real(r8), target, intent(inout) :: hflux_top(psetcols_local)

    call physics_ptend_reset_s_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
         c_loc(s), c_loc(hflux_srf), c_loc(hflux_top))
  end subroutine physics_ptend_reset_s_codon

  subroutine physics_ptend_reset_u_codon(psetcols_local, u, taux_srf, taux_top)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(inout) :: u(psetcols_local,pver)
    real(r8), target, intent(inout) :: taux_srf(psetcols_local)
    real(r8), target, intent(inout) :: taux_top(psetcols_local)

    call physics_ptend_reset_u_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
         c_loc(u), c_loc(taux_srf), c_loc(taux_top))
  end subroutine physics_ptend_reset_u_codon

  subroutine physics_ptend_reset_v_codon(psetcols_local, v, tauy_srf, tauy_top)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(inout) :: v(psetcols_local,pver)
    real(r8), target, intent(inout) :: tauy_srf(psetcols_local)
    real(r8), target, intent(inout) :: tauy_top(psetcols_local)

    call physics_ptend_reset_v_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
         c_loc(v), c_loc(tauy_srf), c_loc(tauy_top))
  end subroutine physics_ptend_reset_v_codon

  subroutine physics_ptend_reset_q_codon(psetcols_local, q, cflx_srf, cflx_top)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(inout) :: q(psetcols_local,pver,pcnst)
    real(r8), target, intent(inout) :: cflx_srf(psetcols_local,pcnst)
    real(r8), target, intent(inout) :: cflx_top(psetcols_local,pcnst)

    call physics_ptend_reset_q_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
         int(pcnst, c_int64_t), c_loc(q), c_loc(cflx_srf), c_loc(cflx_top))
  end subroutine physics_ptend_reset_q_codon

  subroutine physics_ptend_copy_s_codon(psetcols_local, src_s, src_hflux_srf, src_hflux_top, &
       dst_s, dst_hflux_srf, dst_hflux_top)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(in) :: src_s(psetcols_local,pver)
    real(r8), target, intent(in) :: src_hflux_srf(psetcols_local)
    real(r8), target, intent(in) :: src_hflux_top(psetcols_local)
    real(r8), target, intent(inout) :: dst_s(psetcols_local,pver)
    real(r8), target, intent(inout) :: dst_hflux_srf(psetcols_local)
    real(r8), target, intent(inout) :: dst_hflux_top(psetcols_local)

    call physics_ptend_copy_s_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
         c_loc(src_s), c_loc(src_hflux_srf), c_loc(src_hflux_top), &
         c_loc(dst_s), c_loc(dst_hflux_srf), c_loc(dst_hflux_top))
  end subroutine physics_ptend_copy_s_codon

  subroutine physics_ptend_copy_u_codon(psetcols_local, src_u, src_taux_srf, src_taux_top, &
       dst_u, dst_taux_srf, dst_taux_top)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(in) :: src_u(psetcols_local,pver)
    real(r8), target, intent(in) :: src_taux_srf(psetcols_local)
    real(r8), target, intent(in) :: src_taux_top(psetcols_local)
    real(r8), target, intent(inout) :: dst_u(psetcols_local,pver)
    real(r8), target, intent(inout) :: dst_taux_srf(psetcols_local)
    real(r8), target, intent(inout) :: dst_taux_top(psetcols_local)

    call physics_ptend_copy_u_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
         c_loc(src_u), c_loc(src_taux_srf), c_loc(src_taux_top), &
         c_loc(dst_u), c_loc(dst_taux_srf), c_loc(dst_taux_top))
  end subroutine physics_ptend_copy_u_codon

  subroutine physics_ptend_copy_v_codon(psetcols_local, src_v, src_tauy_srf, src_tauy_top, &
       dst_v, dst_tauy_srf, dst_tauy_top)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(in) :: src_v(psetcols_local,pver)
    real(r8), target, intent(in) :: src_tauy_srf(psetcols_local)
    real(r8), target, intent(in) :: src_tauy_top(psetcols_local)
    real(r8), target, intent(inout) :: dst_v(psetcols_local,pver)
    real(r8), target, intent(inout) :: dst_tauy_srf(psetcols_local)
    real(r8), target, intent(inout) :: dst_tauy_top(psetcols_local)

    call physics_ptend_copy_v_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), &
         c_loc(src_v), c_loc(src_tauy_srf), c_loc(src_tauy_top), &
         c_loc(dst_v), c_loc(dst_tauy_srf), c_loc(dst_tauy_top))
  end subroutine physics_ptend_copy_v_codon

  subroutine physics_ptend_copy_q_codon(psetcols_local, src_q, src_cflx_srf, src_cflx_top, &
       dst_q, dst_cflx_srf, dst_cflx_top)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: psetcols_local
    real(r8), target, intent(in) :: src_q(psetcols_local,pver,pcnst)
    real(r8), target, intent(in) :: src_cflx_srf(psetcols_local,pcnst)
    real(r8), target, intent(in) :: src_cflx_top(psetcols_local,pcnst)
    real(r8), target, intent(inout) :: dst_q(psetcols_local,pver,pcnst)
    real(r8), target, intent(inout) :: dst_cflx_srf(psetcols_local,pcnst)
    real(r8), target, intent(inout) :: dst_cflx_top(psetcols_local,pcnst)

    call physics_ptend_copy_q_codon_raw(int(psetcols_local, c_int64_t), int(pver, c_int64_t), int(pcnst, c_int64_t), &
         c_loc(src_q), c_loc(src_cflx_srf), c_loc(src_cflx_top), &
         c_loc(dst_q), c_loc(dst_cflx_srf), c_loc(dst_cflx_top))
  end subroutine physics_ptend_copy_q_codon

  subroutine physics_ptend_scale_field_codon(ncol_local, psetcols_local, top_level_local, bot_level_local, fac_local, &
                                                  field, flx_srf, flx_top)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local, top_level_local, bot_level_local
    real(r8), intent(in) :: fac_local
    real(r8), target, intent(inout) :: field(psetcols_local,pver)
    real(r8), target, intent(inout) :: flx_srf(psetcols_local)
    real(r8), target, intent(inout) :: flx_top(psetcols_local)

    call physics_ptend_scale_field_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(top_level_local, c_int64_t), int(bot_level_local, c_int64_t), real(fac_local, c_double), &
         c_loc(field), c_loc(flx_srf), c_loc(flx_top))
  end subroutine physics_ptend_scale_field_codon

  subroutine physics_ptend_sum_field_codon(ncol_local, psetcols_local, top_level_local, bot_level_local, &
                                                src_field, dst_field, src_flx_srf, dst_flx_srf, &
                                                src_flx_top, dst_flx_top)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local, top_level_local, bot_level_local
    real(r8), target, intent(in) :: src_field(:,:), src_flx_srf(:), src_flx_top(:)
    real(r8), target, intent(inout) :: dst_field(:,:), dst_flx_srf(:), dst_flx_top(:)

    call physics_ptend_sum_field_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(top_level_local, c_int64_t), int(bot_level_local, c_int64_t), c_loc(src_field), &
         c_loc(dst_field), c_loc(src_flx_srf), c_loc(dst_flx_srf), c_loc(src_flx_top), c_loc(dst_flx_top))
  end subroutine physics_ptend_sum_field_codon

  subroutine set_state_pdry_codon(ncol_local, psetcols_local, do_pdeld_calc, psdry, pint, pdel, q, &
                                               pdeldry, pintdry, pmiddry, rpdeldry, lnpmiddry, lnpintdry)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local
    logical, intent(in) :: do_pdeld_calc
    real(r8), target, intent(inout) :: psdry(:), pdeldry(:,:), pintdry(:,:), pmiddry(:,:)
    real(r8), target, intent(inout) :: rpdeldry(:,:), lnpmiddry(:,:), lnpintdry(:,:)
    real(r8), target, intent(in) :: pint(:,:), pdel(:,:), q(:,:,:)

    call set_state_pdry_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(pver, c_int64_t), merge(1_c_int64_t, 0_c_int64_t, do_pdeld_calc), c_loc(psdry), &
         c_loc(pint), c_loc(pdel), c_loc(q), c_loc(pdeldry), c_loc(pintdry), c_loc(pmiddry), &
         c_loc(rpdeldry), c_loc(lnpmiddry), c_loc(lnpintdry))
  end subroutine set_state_pdry_codon

  subroutine physics_set_wet_to_dry_constituent_codon(ncol_local, psetcols_local, q_m, pdel, pdeldry)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local
    real(r8), target, intent(inout) :: q_m(:,:)
    real(r8), target, intent(in) :: pdel(:,:), pdeldry(:,:)

    call physics_set_wet_to_dry_constituent_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(pver, c_int64_t), c_loc(q_m), c_loc(pdel), c_loc(pdeldry))
  end subroutine physics_set_wet_to_dry_constituent_codon

  subroutine physics_set_dry_to_wet_constituent_codon(ncol_local, psetcols_local, q_m, pdeldry, pdel)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local
    real(r8), target, intent(inout) :: q_m(:,:)
    real(r8), target, intent(in) :: pdeldry(:,:), pdel(:,:)

    call physics_set_dry_to_wet_constituent_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(pver, c_int64_t), c_loc(q_m), c_loc(pdeldry), c_loc(pdel))
  end subroutine physics_set_dry_to_wet_constituent_codon

  subroutine set_wet_to_dry_codon(ncol_local, psetcols_local, dry_mask, q, pdel, pdeldry)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local
    integer(c_int64_t), target, intent(in) :: dry_mask(pcnst)
    real(r8), target, intent(inout) :: q(psetcols_local,pver,pcnst)
    real(r8), target, intent(in) :: pdel(psetcols_local,pver), pdeldry(psetcols_local,pver)

    call set_wet_to_dry_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(pver, c_int64_t), int(pcnst, c_int64_t), c_loc(dry_mask), c_loc(q), c_loc(pdel), c_loc(pdeldry))
  end subroutine set_wet_to_dry_codon

  subroutine set_dry_to_wet_codon(ncol_local, psetcols_local, dry_mask, q, pdeldry, pdel)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local
    integer(c_int64_t), target, intent(in) :: dry_mask(pcnst)
    real(r8), target, intent(inout) :: q(psetcols_local,pver,pcnst)
    real(r8), target, intent(in) :: pdeldry(psetcols_local,pver), pdel(psetcols_local,pver)

    call set_dry_to_wet_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(pver, c_int64_t), int(pcnst, c_int64_t), c_loc(dry_mask), c_loc(q), c_loc(pdeldry), c_loc(pdel))
  end subroutine set_dry_to_wet_codon

  subroutine init_geo_unique_codon(ncol_local, psetcols_local, lat, lon, ulat, ulon, &
                                                     latmapback, lonmapback, ulatcnt, uloncnt)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local
    real(r8), target, intent(in) :: lat(:), lon(:)
    real(r8), target, intent(inout) :: ulat(:), ulon(:)
    integer, target, intent(inout) :: latmapback(:), lonmapback(:)
    integer, intent(out) :: ulatcnt, uloncnt
    integer(c_int64_t), target :: ulatcnt_out, uloncnt_out

    call init_geo_unique_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         c_loc(lat), c_loc(lon), c_loc(ulat), c_loc(ulon), c_loc(latmapback), c_loc(lonmapback), &
         c_loc(ulatcnt_out), c_loc(uloncnt_out))
    ulatcnt = int(ulatcnt_out)
    uloncnt = int(uloncnt_out)
  end subroutine init_geo_unique_codon

  subroutine physics_state_set_grid_codon(ncol_local, psetcols_local, rlat, rlon, lat, lon, ulat, ulon, &
                                                     latmapback, lonmapback, ulatcnt, uloncnt)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local
    real(r8), target, intent(in) :: rlat(:), rlon(:)
    real(r8), target, intent(inout) :: lat(:), lon(:), ulat(:), ulon(:)
    integer, target, intent(inout) :: latmapback(:), lonmapback(:)
    integer, intent(out) :: ulatcnt, uloncnt
    integer(c_int64_t), target :: ulatcnt_out, uloncnt_out

    call physics_state_set_grid_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         c_loc(rlat), c_loc(rlon), c_loc(lat), c_loc(lon), c_loc(ulat), c_loc(ulon), &
         c_loc(latmapback), c_loc(lonmapback), c_loc(ulatcnt_out), c_loc(uloncnt_out))
    ulatcnt = int(ulatcnt_out)
    uloncnt = int(uloncnt_out)
  end subroutine physics_state_set_grid_codon

  subroutine physics_copy_real_1d_codon(n_local, src, dst)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: n_local
    real(r8), target, intent(in) :: src(:)
    real(r8), target, intent(inout) :: dst(:)

    call physics_copy_real_1d_codon_raw(int(n_local, c_int64_t), c_loc(src), c_loc(dst))
  end subroutine physics_copy_real_1d_codon

  subroutine physics_copy_real_2d_codon(ncol_local, ld1_local, nlev_local, src, dst)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, ld1_local, nlev_local
    real(r8), target, intent(in) :: src(:,:)
    real(r8), target, intent(inout) :: dst(:,:)

    call physics_copy_real_2d_codon_raw(int(ncol_local, c_int64_t), int(ld1_local, c_int64_t), &
         int(nlev_local, c_int64_t), c_loc(src), c_loc(dst))
  end subroutine physics_copy_real_2d_codon

  subroutine physics_copy_real_3d_codon(ncol_local, ld1_local, nlev_local, n3_local, src, dst)
    use iso_c_binding, only: c_int64_t, c_loc
    integer, intent(in) :: ncol_local, ld1_local, nlev_local, n3_local
    real(r8), target, intent(in) :: src(:,:,:)
    real(r8), target, intent(inout) :: dst(:,:,:)

    call physics_copy_real_3d_codon_raw(int(ncol_local, c_int64_t), int(ld1_local, c_int64_t), &
         int(nlev_local, c_int64_t), int(n3_local, c_int64_t), c_loc(src), c_loc(dst))
  end subroutine physics_copy_real_3d_codon

  subroutine physics_state_copy_codon(state_in, state_out, ncol_local)
    use iso_c_binding, only: c_int64_t, c_loc
    type(physics_state), intent(in), target :: state_in
    type(physics_state), intent(inout), target :: state_out
    integer, intent(in) :: ncol_local

    call physics_state_copy_codon_raw(int(ncol_local, c_int64_t), int(state_in%psetcols, c_int64_t), &
         int(pver, c_int64_t), int(pcnst, c_int64_t), &
         c_loc(state_in%lat), c_loc(state_out%lat), c_loc(state_in%lon), c_loc(state_out%lon), &
         c_loc(state_in%ps), c_loc(state_out%ps), c_loc(state_in%phis), c_loc(state_out%phis), &
         c_loc(state_in%te_ini), c_loc(state_out%te_ini), c_loc(state_in%te_cur), c_loc(state_out%te_cur), &
         c_loc(state_in%tw_ini), c_loc(state_out%tw_ini), c_loc(state_in%tw_cur), c_loc(state_out%tw_cur), &
         c_loc(state_in%psdry), c_loc(state_out%psdry), &
         c_loc(state_in%t), c_loc(state_out%t), c_loc(state_in%u), c_loc(state_out%u), &
         c_loc(state_in%v), c_loc(state_out%v), c_loc(state_in%s), c_loc(state_out%s), &
         c_loc(state_in%omega), c_loc(state_out%omega), c_loc(state_in%pmid), c_loc(state_out%pmid), &
         c_loc(state_in%pdel), c_loc(state_out%pdel), c_loc(state_in%rpdel), c_loc(state_out%rpdel), &
         c_loc(state_in%lnpmid), c_loc(state_out%lnpmid), c_loc(state_in%exner), c_loc(state_out%exner), &
         c_loc(state_in%zm), c_loc(state_out%zm), c_loc(state_in%lnpmiddry), c_loc(state_out%lnpmiddry), &
         c_loc(state_in%pmiddry), c_loc(state_out%pmiddry), c_loc(state_in%pdeldry), c_loc(state_out%pdeldry), &
         c_loc(state_in%rpdeldry), c_loc(state_out%rpdeldry), c_loc(state_in%pint), c_loc(state_out%pint), &
         c_loc(state_in%lnpint), c_loc(state_out%lnpint), c_loc(state_in%zi), c_loc(state_out%zi), &
         c_loc(state_in%pintdry), c_loc(state_out%pintdry), c_loc(state_in%lnpintdry), c_loc(state_out%lnpintdry), &
         c_loc(state_in%q), c_loc(state_out%q))
  end subroutine physics_state_copy_codon

  subroutine state_cnst_min_nz_codon(ncol_local, psetcols_local, q, lim, qix, numix)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local, qix, numix
    real(r8), target, intent(inout) :: q(:,:,:)
    real(r8), intent(in) :: lim

    call state_cnst_min_nz_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(pver, c_int64_t), int(pcnst, c_int64_t), c_loc(q), real(lim, c_double), &
         int(qix, c_int64_t), int(numix, c_int64_t))
  end subroutine state_cnst_min_nz_codon

  subroutine physics_update_field_codon(ncol_local, psetcols_local, top_level_local, bot_level_local, &
       dt_local, update_tend, state_field, ptend_field, tend_field)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local, top_level_local, bot_level_local
    real(r8), intent(in) :: dt_local
    logical, intent(in) :: update_tend
    real(r8), target, intent(inout) :: state_field(:,:)
    real(r8), target, intent(in) :: ptend_field(:,:)
    real(r8), target, intent(inout) :: tend_field(:,:)

    call physics_update_field_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(pver, c_int64_t), int(top_level_local, c_int64_t), int(bot_level_local, c_int64_t), &
         real(dt_local, c_double), merge(1_c_int64_t, 0_c_int64_t, update_tend), &
         c_loc(state_field), c_loc(ptend_field), c_loc(tend_field))
  end subroutine physics_update_field_codon

  subroutine physics_update_q_codon(ncol_local, psetcols_local, top_level_local, bot_level_local, dt_local, &
       m_local, is_number, state_q, ptend_q)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local, top_level_local, bot_level_local, m_local
    logical, intent(in) :: is_number
    real(r8), intent(in) :: dt_local
    real(r8), target, intent(inout) :: state_q(:,:,:)
    real(r8), target, intent(in) :: ptend_q(:,:,:)

    call physics_update_q_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(pver, c_int64_t), int(pcnst, c_int64_t), int(top_level_local, c_int64_t), &
         int(bot_level_local, c_int64_t), real(dt_local, c_double), int(m_local, c_int64_t), &
         merge(1_c_int64_t, 0_c_int64_t, is_number), c_loc(state_q), c_loc(ptend_q))
  end subroutine physics_update_q_codon

  subroutine physics_update_s_codon(ncol_local, psetcols_local, top_level_local, bot_level_local, &
       dt_local, update_tend, state_s, ptend_s, tend_dtdt, cpairv_loc_chunk)
    use iso_c_binding, only: c_double, c_int64_t, c_loc
    integer, intent(in) :: ncol_local, psetcols_local, top_level_local, bot_level_local
    real(r8), intent(in) :: dt_local
    logical, intent(in) :: update_tend
    real(r8), target, intent(inout) :: state_s(:,:)
    real(r8), target, intent(in) :: ptend_s(:,:), cpairv_loc_chunk(:,:)
    real(r8), target, intent(inout) :: tend_dtdt(:,:)

    call physics_update_s_codon_raw(int(ncol_local, c_int64_t), int(psetcols_local, c_int64_t), &
         int(pver, c_int64_t), int(top_level_local, c_int64_t), int(bot_level_local, c_int64_t), &
         real(dt_local, c_double), merge(1_c_int64_t, 0_c_int64_t, update_tend), &
         c_loc(state_s), c_loc(ptend_s), c_loc(tend_dtdt), c_loc(cpairv_loc_chunk))
  end subroutine physics_update_s_codon

!===============================================================================
  subroutine physics_type_alloc(phys_state, phys_tend, begchunk, endchunk, psetcols)
    use iso_c_binding, only: c_int64_t
    implicit none
    type(physics_state), pointer :: phys_state(:)
    type(physics_tend), pointer :: phys_tend(:)
    integer, intent(in) :: begchunk, endchunk
    integer, intent(in) :: psetcols

    integer :: ierr=0, lchnk
    type(physics_state), pointer :: state
    type(physics_tend), pointer :: tend
    integer(c_int64_t) :: touch_c

    call physics_types_zero_select_impl()
    if (.not. use_native_zero_impl) then
       touch_c = physics_type_alloc_codon_raw(1_c_int64_t)
       if (touch_c == 1_c_int64_t) then
          call physics_types_zero_proof_once()
          call physics_types_log_direct(physics_type_alloc_logged, &
               'physics_type_alloc direct = codon; allocation parent selector/touch direct = codon; ' // &
               'native allocate/CAM state container island')
       end if
    end if

    allocate(phys_state(begchunk:endchunk), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'physics_types: phys_state allocation error = ',ierr
       call endrun('physics_types: failed to allocate physics_state array')
    end if

    do lchnk=begchunk,endchunk
       call physics_state_alloc(phys_state(lchnk),lchnk,pcols)
    end do

    allocate(phys_tend(begchunk:endchunk), stat=ierr)
    if( ierr /= 0 ) then
       write(iulog,*) 'physics_types: phys_tend allocation error = ',ierr
       call endrun('physics_types: failed to allocate physics_tend array')
    end if

    do lchnk=begchunk,endchunk
       call physics_tend_alloc(phys_tend(lchnk),phys_state(lchnk)%psetcols)
    end do

  end subroutine physics_type_alloc
!===============================================================================
  subroutine physics_update(state, ptend, dt, tend)
!-----------------------------------------------------------------------
! Update the state and or tendency structure with the parameterization tendencies
!-----------------------------------------------------------------------
    use iso_c_binding, only: c_int64_t
    use shr_sys_mod,  only: shr_sys_flush
    use geopotential, only: geopotential_dse
    use constituents, only: cnst_get_ind, cnst_mw
    use scamMod,      only: scm_crm_mode, single_column
    use phys_control, only: phys_getopts
    use physconst,    only: physconst_update ! Routine which updates physconst variables (WACCM-X)
    use ppgrid,       only: begchunk, endchunk
    use water_tracer_vars, only: trace_water, wtrc_nwset, wtrc_iatype
    use water_types,  only: iwtliq, iwtice

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(inout)  :: ptend   ! Parameterization tendencies

    type(physics_state), intent(inout)  :: state   ! Physics state variables

    real(r8), intent(in) :: dt                     ! time step

    type(physics_tend ), intent(inout), optional  :: tend  ! Physics tendencies over timestep
                    ! This is usually only needed by calls from physpkg.
!
!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: ixcldice, ixcldliq                  ! indices for CLDICE and CLDLIQ
    integer :: ixnumice, ixnumliq
    integer :: ixnumsnow, ixnumrain
    integer :: ncol                                ! number of columns
    character*40 :: name    ! param and tracer name for qneg3

    integer :: ixo, ixo2, ixh, ixh2, ixn    ! indices for O, O2, H2, and N

    real(r8) :: zvirv(state%psetcols,pver)  ! Local zvir array pointer

    real(r8),allocatable :: cpairv_loc(:,:,:)
    real(r8),allocatable :: rairv_loc(:,:,:)
    real(r8), target :: dummy_tend_field(1,1)

    ! PERGRO limits cldliq/ice for macro/microphysics:
    character(len=24), parameter :: pergro_cldlim_names(4) = &
         (/ "stratiform", "cldwat    ", "micro_mg  ", "macro_park" /)

    ! cldliq/ice limits that are always on.
    character(len=24), parameter :: cldlim_names(2) = &
         (/ "convect_deep", "zm_conv_tend" /)

    ! Whether to do validation of state on each call.
    logical :: state_debug_checks

    !-----------------------------------------------------------------------

    call physics_types_zero_select_impl()
    call physics_update_select_impl()
    if (.not. use_native_zero_impl) then
       if (physics_types_touch_codon(2_c_int64_t) == 2_c_int64_t) then
          call physics_types_zero_proof_once()
          call physics_types_log_direct(physics_update_logged, &
               'physics_update direct = codon; active state/tendency update loops direct = codon; ' // &
               'qneg3/geopotential/WACCM native boundaries')
       end if
    end if

    ! The column radiation model does not update the state
    if(single_column.and.scm_crm_mode) return


    !-----------------------------------------------------------------------
    ! If no fields are set, then return
    if (.not. (any(ptend%lq(:)) .or. ptend%ls .or. ptend%lu .or. ptend%lv)) then
       ptend%name  = "none"
       ptend%psetcols = 0
       return
    end if

    !-----------------------------------------------------------------------
    ! Check that the state/tend/ptend are all dimensioned with the same number of columns
    if (state%psetcols /= ptend%psetcols) then
       call endrun('ERROR in physics_update with ptend%name='//trim(ptend%name) &
            //': state and ptend must have the same number of psetcols.')
    end if

    if (present(tend)) then
       if (state%psetcols /= tend%psetcols) then
          call endrun('ERROR in physics_update with ptend%name='//trim(ptend%name) &
               //': state and tend must have the same number of psetcols.')
       end if
    end if

    !-----------------------------------------------------------------------
    ! cpairv_loc and rairv_loc need to be allocated to a size which matches state and ptend
    ! If psetcols == pcols, the cpairv is the correct size and just copy
    ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair
    if (state%psetcols == pcols) then
       allocate (cpairv_loc(state%psetcols,pver,begchunk:endchunk))
       cpairv_loc(:,:,:) = cpairv(:,:,:)
    else if (state%psetcols > pcols .and. all(cpairv(:,:,:) == cpair)) then
       allocate(cpairv_loc(state%psetcols,pver,begchunk:endchunk))
       cpairv_loc(:,:,:) = cpair
    else
       call endrun('physics_update: cpairv is not allowed to vary when subcolumns are turned on')
    end if
    if (state%psetcols == pcols) then
       allocate (rairv_loc(state%psetcols,pver,begchunk:endchunk))
       rairv_loc(:,:,:) = rairv(:,:,:)
    else if (state%psetcols > pcols .and. all(rairv(:,:,:) == rair)) then
       allocate(rairv_loc(state%psetcols,pver,begchunk:endchunk))
       rairv_loc(:,:,:) = rair
    else
       call endrun('physics_update: rairv_loc is not allowed to vary when subcolumns are turned on')
    end if

    !-----------------------------------------------------------------------
    call phys_getopts(state_debug_checks_out=state_debug_checks)

    ncol = state%ncol

    ! Update u,v fields
    if(ptend%lu) then
       if (use_native_zero_impl .or. physics_update_native_field) then
          do k = ptend%top_level, ptend%bot_level
             state%u  (:ncol,k) = state%u  (:ncol,k) + ptend%u(:ncol,k) * dt
             if (present(tend)) &
                  tend%dudt(:ncol,k) = tend%dudt(:ncol,k) + ptend%u(:ncol,k)
          end do
       else if (present(tend)) then
          call physics_update_field_codon(ncol, state%psetcols, ptend%top_level, ptend%bot_level, &
               dt, .true., state%u, ptend%u, tend%dudt)
       else
          call physics_update_field_codon(ncol, state%psetcols, ptend%top_level, ptend%bot_level, &
               dt, .false., state%u, ptend%u, dummy_tend_field)
       end if
    end if

    if(ptend%lv) then
       if (use_native_zero_impl .or. physics_update_native_field) then
          do k = ptend%top_level, ptend%bot_level
             state%v  (:ncol,k) = state%v  (:ncol,k) + ptend%v(:ncol,k) * dt
             if (present(tend)) &
                  tend%dvdt(:ncol,k) = tend%dvdt(:ncol,k) + ptend%v(:ncol,k)
          end do
       else if (present(tend)) then
          call physics_update_field_codon(ncol, state%psetcols, ptend%top_level, ptend%bot_level, &
               dt, .true., state%v, ptend%v, tend%dvdt)
       else
          call physics_update_field_codon(ncol, state%psetcols, ptend%top_level, ptend%bot_level, &
               dt, .false., state%v, ptend%v, dummy_tend_field)
       end if
    end if

   ! Update constituents, all schemes use time split q: no tendency kept
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    ! Check for number concentration of cloud liquid and cloud ice (if not present
    ! the indices will be set to -1)
    call cnst_get_ind('NUMICE', ixnumice, abort=.false.)
    call cnst_get_ind('NUMLIQ', ixnumliq, abort=.false.)
    call cnst_get_ind('NUMRAI', ixnumrain, abort=.false.)
    call cnst_get_ind('NUMSNO', ixnumsnow, abort=.false.)

    do m = 1, pcnst
       if(ptend%lq(m)) then
          ! now test for mixing ratios which are too small
          ! don't call qneg3 for number concentration variables
          if (m /= ixnumice  .and.  m /= ixnumliq .and. &
              m /= ixnumrain .and.  m /= ixnumsnow ) then
             if (use_native_zero_impl .or. physics_update_native_q) then
                do k = ptend%top_level, ptend%bot_level
                   state%q(:ncol,k,m) = state%q(:ncol,k,m) + ptend%q(:ncol,k,m) * dt
                end do
             else
                call physics_update_q_codon(ncol, state%psetcols, ptend%top_level, ptend%bot_level, &
                     dt, m, .false., state%q, ptend%q)
             end if
             name = trim(ptend%name) // '/' // trim(cnst_name(m))
             call qneg3(trim(name), state%lchnk, ncol, state%psetcols, pver, m, m, qmin(m), state%q(1,1,m))
          else
             if (use_native_zero_impl .or. physics_update_native_q) then
                do k = ptend%top_level, ptend%bot_level
                   state%q(:ncol,k,m) = state%q(:ncol,k,m) + ptend%q(:ncol,k,m) * dt
                   ! checks for number concentration
                   state%q(:ncol,k,m) = max(1.e-12_r8,state%q(:ncol,k,m))
                   state%q(:ncol,k,m) = min(1.e10_r8,state%q(:ncol,k,m))
                end do
             else
                call physics_update_q_codon(ncol, state%psetcols, ptend%top_level, ptend%bot_level, &
                     dt, m, .true., state%q, ptend%q)
             end if
          end if

       end if

    end do

    !------------------------------------------------------------------------
    ! This is a temporary fix for the large H, H2 in WACCM-X
    ! Well, it was supposed to be temporary, but it has been here
    ! for a while now.
    !------------------------------------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       call cnst_get_ind('H', ixh)
       do k = ptend%top_level, ptend%bot_level
          state%q(:ncol,k,ixh) = min(state%q(:ncol,k,ixh), 0.01_r8)
       end do

       call cnst_get_ind('H2', ixh2)
       do k = ptend%top_level, ptend%bot_level
          state%q(:ncol,k,ixh2) = min(state%q(:ncol,k,ixh2), 6.e-5_r8)
       end do
    endif

    ! Special tests for cloud liquid and ice:
    ! Enforce a minimum non-zero value.
    if (ixcldliq > 1) then
       if(ptend%lq(ixcldliq)) then
#ifdef PERGRO
          if ( any(ptend%name == pergro_cldlim_names) ) &
               call state_cnst_min_nz(1.e-12_r8, ixcldliq, ixnumliq)
#endif
          if ( any(ptend%name == cldlim_names) ) &
               call state_cnst_min_nz(1.e-36_r8, ixcldliq, ixnumliq)
       end if
    end if

    if (ixcldice > 1) then
       if(ptend%lq(ixcldice)) then
#ifdef PERGRO
          if ( any(ptend%name == pergro_cldlim_names) ) &
               call state_cnst_min_nz(1.e-12_r8, ixcldice, ixnumice)
#endif
          if ( any(ptend%name == cldlim_names) ) &
               call state_cnst_min_nz(1.e-36_r8, ixcldice, ixnumice)
       end if
    end if

    !**************************************************************
    !special tests for water tracers (to match cloud water and ice)
    !**************************************************************
    if(trace_water) then !are water tracers on?

      !NOTE:  This code might not work for water isotopes, as the values
      !will be significantly smaller.  If a problem occurs when doing
      !ratio or similar tests, this could likely be the culprit. - JN

      !NOTE:  Make sure to zero out the water tracer only where the actual
      !bulk water satisfies the logical condition, not where the water tracer
      !itself passes. -JN

      do m=1,wtrc_nwset !loop over water tracers
        !-------------
        !Cloud liquid:
        !-------------
        if(ptend%lq(wtrc_iatype(m,iwtliq))) then
#ifdef PERGRO
          if ( any(ptend%name == pergro_cldlim_names) ) &
            call state_cnst_min_nz(1.e-12_r8, ixcldliq, wtrc_iatype(m,iwtliq))
#endif
          if ( any(ptend%name == cldlim_names) ) &
            call state_cnst_min_nz(1.e-36_r8, ixcldliq, wtrc_iatype(m,iwtliq))
        end if
        !---------
        !Cloud Ice:
        !---------
        if(ptend%lq(wtrc_iatype(m,iwtice))) then
#ifdef PERGRO
          if ( any(ptend%name == pergro_cldlim_names) ) &
            call state_cnst_min_nz(1.e-12_r8, ixcldice, wtrc_iatype(m,iwtice))
#endif
          if ( any(ptend%name == cldlim_names) ) &
            call state_cnst_min_nz(1.e-36_r8, ixcldice, wtrc_iatype(m,iwtice))
        end if
        !---------
      end do !water tracers
    end if
    !**************************************************************

    !------------------------------------------------------------------------
    ! Get indices for molecular weights and call WACCM-X physconst_update
    !------------------------------------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
      call cnst_get_ind('O', ixo)
      call cnst_get_ind('O2', ixo2)
      call cnst_get_ind('N', ixn)

      call physconst_update(state%q, state%t, &
                            cnst_mw(ixo), cnst_mw(ixo2), cnst_mw(ixh), cnst_mw(ixn), &
                            ixo, ixo2, ixh, pcnst, state%lchnk, ncol)
    endif

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
      zvirv(:,:) = shr_const_rwv / rairv_loc(:,:,state%lchnk) - 1._r8
    else
      zvirv(:,:) = zvir
    endif

    !-------------------------------------------------------------------------------------------
    ! Update dry static energy(moved from above for WACCM-X so updating after cpairv_loc update)
    !-------------------------------------------------------------------------------------------
    if(ptend%ls) then
       if (use_native_zero_impl .or. physics_update_native_s) then
          do k = ptend%top_level, ptend%bot_level
             state%s(:ncol,k)   = state%s(:ncol,k)   + ptend%s(:ncol,k) * dt
             if (present(tend)) &
                  tend%dtdt(:ncol,k) = tend%dtdt(:ncol,k) + ptend%s(:ncol,k)/cpairv_loc(:ncol,k,state%lchnk)
          end do
       else if (present(tend)) then
          call physics_update_s_codon(ncol, state%psetcols, ptend%top_level, ptend%bot_level, &
               dt, .true., state%s, ptend%s, tend%dtdt, cpairv_loc(:,:,state%lchnk))
       else
          call physics_update_s_codon(ncol, state%psetcols, ptend%top_level, ptend%bot_level, &
               dt, .false., state%s, ptend%s, dummy_tend_field, cpairv_loc(:,:,state%lchnk))
       end if
    end if

    ! Derive new temperature and geopotential fields if heating or water tendency not 0.
    if (ptend%ls .or. ptend%lq(1)) then
       call geopotential_dse(  &
            state%lnpint, state%lnpmid, state%pint  , state%pmid  , state%pdel  , state%rpdel  , &
            state%s     , state%q(:,:,1),state%phis , rairv_loc(:,:,state%lchnk), gravit  , cpairv_loc(:,:,state%lchnk), &
            zvirv    , state%t     , state%zi    , state%zm    , ncol         )
    end if

    ! Good idea to do this regularly.
    call shr_sys_flush(iulog)

    if (state_debug_checks) call physics_state_check(state, ptend%name)

    deallocate(cpairv_loc, rairv_loc)

    ! Deallocate ptend
    call physics_ptend_dealloc(ptend)

    ptend%name  = "none"
    ptend%lq(:) = .false.
    ptend%ls    = .false.
    ptend%lu    = .false.
    ptend%lv    = .false.
    ptend%psetcols = 0

  contains

    subroutine state_cnst_min_nz(lim, qix, numix)
      ! Small utility function for setting minimum nonzero
      ! constituent concentrations.

      ! Lower limit and constituent index
      real(r8), intent(in) :: lim
      integer,  intent(in) :: qix
      ! Number concentration that goes with qix.
      ! Ignored if <= 0 (and therefore constituent is not present).
      integer,  intent(in) :: numix
      integer :: i, k

      if (use_native_zero_impl) then
         do k = 1, pver
            do i = 1, ncol
               if (state%q(i,k,qix) < lim) then
                  state%q(i,k,qix) = 0._r8
                  if (numix > 0) then
                     state%q(i,k,numix) = 0._r8
                  end if
               end if
            end do
         end do
         if (.not. state_cnst_min_nz_logged) then
            call physics_types_log_direct(state_cnst_min_nz_logged, 'state_cnst_min_nz direct = native')
         end if
      else
         call state_cnst_min_nz_codon(ncol, state%psetcols, state%q, lim, qix, numix)
         if (.not. state_cnst_min_nz_logged) then
            call physics_types_zero_proof_once()
            call physics_types_log_direct(state_cnst_min_nz_logged, 'state_cnst_min_nz direct = codon')
         end if
      end if

    end subroutine state_cnst_min_nz


  end subroutine physics_update

!===============================================================================

  subroutine physics_state_check(state, name)
!-----------------------------------------------------------------------
! Check a physics_state object for invalid data (e.g NaNs, negative
! temperatures).
!-----------------------------------------------------------------------
    use shr_infnan_mod, only: shr_infnan_inf_type, assignment(=), &
                              shr_infnan_posinf, shr_infnan_neginf
    use shr_assert_mod, only: shr_assert, shr_assert_in_domain
    use physconst,      only: pi
    use constituents,   only: pcnst, qmin

!------------------------------Arguments--------------------------------
    ! State to check.
    type(physics_state), intent(in)           :: state
    ! Name of the package responsible for this state.
    character(len=*),    intent(in), optional :: name

!---------------------------Local storage-------------------------------
    ! Shortened name for ncol.
    integer :: ncol
    ! Double precision positive/negative infinity.
    real(r8) :: posinf_r8, neginf_r8
    ! Canned message.
    character(len=64) :: msg
    ! Constituent index
    integer :: m

!-----------------------------------------------------------------------

    ncol = state%ncol

    posinf_r8 = shr_infnan_posinf
    neginf_r8 = shr_infnan_neginf

    ! It may be reasonable to check some of the integer components of the
    ! state as well, but this is not yet implemented.

    ! Check for NaN first to avoid any IEEE exceptions.

    if (present(name)) then
       msg = "NaN produced in physics_state by package "// &
            trim(name)//"."
    else
       msg = "NaN found in physics_state."
    end if

    ! 1-D variables
    call shr_assert_in_domain(state%ps(:ncol),          is_nan=.false., &
         varname="state%ps",        msg=msg)
    call shr_assert_in_domain(state%psdry(:ncol),       is_nan=.false., &
         varname="state%psdry",     msg=msg)
    call shr_assert_in_domain(state%phis(:ncol),        is_nan=.false., &
         varname="state%phis",      msg=msg)
    call shr_assert_in_domain(state%te_ini(:ncol),      is_nan=.false., &
         varname="state%te_ini",    msg=msg)
    call shr_assert_in_domain(state%te_cur(:ncol),      is_nan=.false., &
         varname="state%te_cur",    msg=msg)
    call shr_assert_in_domain(state%tw_ini(:ncol),      is_nan=.false., &
         varname="state%tw_ini",    msg=msg)
    call shr_assert_in_domain(state%tw_cur(:ncol),      is_nan=.false., &
         varname="state%tw_cur",    msg=msg)

    ! 2-D variables (at midpoints)
    call shr_assert_in_domain(state%t(:ncol,:),         is_nan=.false., &
         varname="state%t",         msg=msg)
    call shr_assert_in_domain(state%u(:ncol,:),         is_nan=.false., &
         varname="state%u",         msg=msg)
    call shr_assert_in_domain(state%v(:ncol,:),         is_nan=.false., &
         varname="state%v",         msg=msg)
    call shr_assert_in_domain(state%s(:ncol,:),         is_nan=.false., &
         varname="state%s",         msg=msg)
    call shr_assert_in_domain(state%omega(:ncol,:),     is_nan=.false., &
         varname="state%omega",     msg=msg)
    call shr_assert_in_domain(state%pmid(:ncol,:),      is_nan=.false., &
         varname="state%pmid",      msg=msg)
    call shr_assert_in_domain(state%pmiddry(:ncol,:),   is_nan=.false., &
         varname="state%pmiddry",   msg=msg)
    call shr_assert_in_domain(state%pdel(:ncol,:),      is_nan=.false., &
         varname="state%pdel",      msg=msg)
    call shr_assert_in_domain(state%pdeldry(:ncol,:),   is_nan=.false., &
         varname="state%pdeldry",   msg=msg)
    call shr_assert_in_domain(state%rpdel(:ncol,:),     is_nan=.false., &
         varname="state%rpdel",     msg=msg)
    call shr_assert_in_domain(state%rpdeldry(:ncol,:),  is_nan=.false., &
         varname="state%rpdeldry",  msg=msg)
    call shr_assert_in_domain(state%lnpmid(:ncol,:),    is_nan=.false., &
         varname="state%lnpmid",    msg=msg)
    call shr_assert_in_domain(state%lnpmiddry(:ncol,:), is_nan=.false., &
         varname="state%lnpmiddry", msg=msg)
    call shr_assert_in_domain(state%exner(:ncol,:),     is_nan=.false., &
         varname="state%exner",     msg=msg)
    call shr_assert_in_domain(state%zm(:ncol,:),        is_nan=.false., &
         varname="state%zm",        msg=msg)

    ! 2-D variables (at interfaces)
    call shr_assert_in_domain(state%pint(:ncol,:),      is_nan=.false., &
         varname="state%pint",      msg=msg)
    call shr_assert_in_domain(state%pintdry(:ncol,:),   is_nan=.false., &
         varname="state%pintdry",   msg=msg)
    call shr_assert_in_domain(state%lnpint(:ncol,:),    is_nan=.false., &
         varname="state%lnpint",    msg=msg)
    call shr_assert_in_domain(state%lnpintdry(:ncol,:), is_nan=.false., &
         varname="state%lnpintdry", msg=msg)
    call shr_assert_in_domain(state%zi(:ncol,:),        is_nan=.false., &
         varname="state%zi",        msg=msg)

    ! 3-D variables
    call shr_assert_in_domain(state%q(:ncol,:,:),       is_nan=.false., &
         varname="state%q",         msg=msg)

    ! Now run other checks (i.e. values are finite and within a range that
    ! is physically meaningful).

    if (present(name)) then
       msg = "Invalid value produced in physics_state by package "// &
            trim(name)//"."
    else
       msg = "Invalid value found in physics_state."
    end if

    ! 1-D variables
    call shr_assert_in_domain(state%ps(:ncol),          lt=posinf_r8, gt=0._r8, &
         varname="state%ps",        msg=msg)
    call shr_assert_in_domain(state%psdry(:ncol),       lt=posinf_r8, gt=0._r8, &
         varname="state%psdry",     msg=msg)
    call shr_assert_in_domain(state%phis(:ncol),        lt=posinf_r8, gt=neginf_r8, &
         varname="state%phis",      msg=msg)
    call shr_assert_in_domain(state%te_ini(:ncol),      lt=posinf_r8, gt=neginf_r8, &
         varname="state%te_ini",    msg=msg)
    call shr_assert_in_domain(state%te_cur(:ncol),      lt=posinf_r8, gt=neginf_r8, &
         varname="state%te_cur",    msg=msg)
    call shr_assert_in_domain(state%tw_ini(:ncol),      lt=posinf_r8, gt=neginf_r8, &
         varname="state%tw_ini",    msg=msg)
    call shr_assert_in_domain(state%tw_cur(:ncol),      lt=posinf_r8, gt=neginf_r8, &
         varname="state%tw_cur",    msg=msg)

    ! 2-D variables (at midpoints)
    call shr_assert_in_domain(state%t(:ncol,:),         lt=posinf_r8, gt=0._r8, &
         varname="state%t",         msg=msg)
    call shr_assert_in_domain(state%u(:ncol,:),         lt=posinf_r8, gt=neginf_r8, &
         varname="state%u",         msg=msg)
    call shr_assert_in_domain(state%v(:ncol,:),         lt=posinf_r8, gt=neginf_r8, &
         varname="state%v",         msg=msg)
    call shr_assert_in_domain(state%s(:ncol,:),         lt=posinf_r8, gt=neginf_r8, &
         varname="state%s",         msg=msg)
    call shr_assert_in_domain(state%omega(:ncol,:),     lt=posinf_r8, gt=neginf_r8, &
         varname="state%omega",     msg=msg)
    call shr_assert_in_domain(state%pmid(:ncol,:),      lt=posinf_r8, gt=0._r8, &
         varname="state%pmid",      msg=msg)
    call shr_assert_in_domain(state%pmiddry(:ncol,:),   lt=posinf_r8, gt=0._r8, &
         varname="state%pmiddry",   msg=msg)
    call shr_assert_in_domain(state%pdel(:ncol,:),      lt=posinf_r8, gt=neginf_r8, &
         varname="state%pdel",      msg=msg)
    call shr_assert_in_domain(state%pdeldry(:ncol,:),   lt=posinf_r8, gt=neginf_r8, &
         varname="state%pdeldry",   msg=msg)
    call shr_assert_in_domain(state%rpdel(:ncol,:),     lt=posinf_r8, gt=neginf_r8, &
         varname="state%rpdel",     msg=msg)
    call shr_assert_in_domain(state%rpdeldry(:ncol,:),  lt=posinf_r8, gt=neginf_r8, &
         varname="state%rpdeldry",  msg=msg)
    call shr_assert_in_domain(state%lnpmid(:ncol,:),    lt=posinf_r8, gt=neginf_r8, &
         varname="state%lnpmid",    msg=msg)
    call shr_assert_in_domain(state%lnpmiddry(:ncol,:), lt=posinf_r8, gt=neginf_r8, &
         varname="state%lnpmiddry", msg=msg)
    call shr_assert_in_domain(state%exner(:ncol,:),     lt=posinf_r8, gt=0._r8, &
         varname="state%exner",     msg=msg)
    call shr_assert_in_domain(state%zm(:ncol,:),        lt=posinf_r8, gt=neginf_r8, &
         varname="state%zm",        msg=msg)

    ! 2-D variables (at interfaces)
    call shr_assert_in_domain(state%pint(:ncol,:),      lt=posinf_r8, gt=0._r8, &
         varname="state%pint",      msg=msg)
    call shr_assert_in_domain(state%pintdry(:ncol,:),   lt=posinf_r8, gt=0._r8, &
         varname="state%pintdry",   msg=msg)
    call shr_assert_in_domain(state%lnpint(:ncol,:),    lt=posinf_r8, gt=neginf_r8, &
         varname="state%lnpint",    msg=msg)
    call shr_assert_in_domain(state%lnpintdry(:ncol,:), lt=posinf_r8, gt=neginf_r8, &
         varname="state%lnpintdry", msg=msg)
    call shr_assert_in_domain(state%zi(:ncol,:),        lt=posinf_r8, gt=neginf_r8, &
         varname="state%zi",        msg=msg)

    ! 3-D variables
    do m = 1,pcnst
       call shr_assert_in_domain(state%q(:ncol,:,m),    lt=posinf_r8, ge=qmin(m), &
            varname="state%q ("//trim(cnst_name(m))//")", msg=msg)
    end do

  end subroutine physics_state_check

!===============================================================================

  subroutine physics_ptend_sum(ptend, ptend_sum, ncol)
    use iso_c_binding, only: c_int64_t, c_loc, c_ptr, c_null_ptr
!-----------------------------------------------------------------------
! Add ptend fields to ptend_sum for ptend logical flags = .true.
! Where ptend logical flags = .false, don't change ptend_sum
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), target, intent(in)     :: ptend   ! New parameterization tendencies
    type(physics_ptend), target, intent(inout)  :: ptend_sum   ! Sum of incoming ptend_sum and ptend
    integer, intent(in)                 :: ncol    ! number of columns

!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices
    integer :: psetcols                            ! maximum number of columns
    integer :: ierr = 0
    logical :: used_codon
    integer(c_int64_t), target :: lq_mask(pcnst)
    type(c_ptr) :: src_s_p, dst_s_p, src_u_p, dst_u_p, src_v_p, dst_v_p, src_q_p, dst_q_p
    type(c_ptr) :: src_hflux_srf_p, dst_hflux_srf_p, src_hflux_top_p, dst_hflux_top_p
    type(c_ptr) :: src_taux_srf_p, dst_taux_srf_p, src_taux_top_p, dst_taux_top_p
    type(c_ptr) :: src_tauy_srf_p, dst_tauy_srf_p, src_tauy_top_p, dst_tauy_top_p
    type(c_ptr) :: src_cflx_srf_p, dst_cflx_srf_p, src_cflx_top_p, dst_cflx_top_p

!-----------------------------------------------------------------------
    if (ptend%psetcols /= ptend_sum%psetcols) then
       call endrun('physics_ptend_sum error: ptend and ptend_sum must have the same value for psetcols')
    end if

    if (ncol > ptend_sum%psetcols) then
       call endrun('physics_ptend_sum error: ncol must be less than or equal to psetcols')
    end if

    psetcols = ptend_sum%psetcols
    used_codon = .false.
    lq_mask = 0_c_int64_t
    call physics_types_zero_select_impl()

    ptend_sum%top_level = ptend%top_level
    ptend_sum%bot_level = ptend%bot_level

! Update u,v fields
    if(ptend%lu) then
       if (.not. allocated(ptend_sum%u)) then
          allocate(ptend_sum%u(psetcols,pver), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%u')
          ptend_sum%u=0.0_r8

          allocate(ptend_sum%taux_srf(psetcols), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%taux_srf')
          ptend_sum%taux_srf=0.0_r8

          allocate(ptend_sum%taux_top(psetcols), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%taux_top')
          ptend_sum%taux_top=0.0_r8
       end if
       ptend_sum%lu = .true.

       if (use_native_zero_impl) then
          do k = ptend%top_level, ptend%bot_level
             do i = 1, ncol
                ptend_sum%u(i,k) = ptend_sum%u(i,k) + ptend%u(i,k)
             end do
          end do
          do i = 1, ncol
             ptend_sum%taux_srf(i) = ptend_sum%taux_srf(i) + ptend%taux_srf(i)
             ptend_sum%taux_top(i) = ptend_sum%taux_top(i) + ptend%taux_top(i)
          end do
       else
          used_codon = .true.
       end if
    end if

    if(ptend%lv) then
       if (.not. allocated(ptend_sum%v)) then
          allocate(ptend_sum%v(psetcols,pver), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%v')
          ptend_sum%v=0.0_r8

          allocate(ptend_sum%tauy_srf(psetcols), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%tauy_srf')
          ptend_sum%tauy_srf=0.0_r8

          allocate(ptend_sum%tauy_top(psetcols), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%tauy_top')
          ptend_sum%tauy_top=0.0_r8
       end if
       ptend_sum%lv = .true.

       if (use_native_zero_impl) then
          do k = ptend%top_level, ptend%bot_level
             do i = 1, ncol
                ptend_sum%v(i,k) = ptend_sum%v(i,k) + ptend%v(i,k)
             end do
          end do
          do i = 1, ncol
             ptend_sum%tauy_srf(i) = ptend_sum%tauy_srf(i) + ptend%tauy_srf(i)
             ptend_sum%tauy_top(i) = ptend_sum%tauy_top(i) + ptend%tauy_top(i)
          end do
       else
          used_codon = .true.
       end if
    end if


    if(ptend%ls) then
       if (.not. allocated(ptend_sum%s)) then
          allocate(ptend_sum%s(psetcols,pver), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%s')
          ptend_sum%s=0.0_r8

          allocate(ptend_sum%hflux_srf(psetcols), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%hflux_srf')
          ptend_sum%hflux_srf=0.0_r8

          allocate(ptend_sum%hflux_top(psetcols), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%hflux_top')
          ptend_sum%hflux_top=0.0_r8
       end if
       ptend_sum%ls = .true.

       if (use_native_zero_impl) then
          do k = ptend%top_level, ptend%bot_level
             do i = 1, ncol
                ptend_sum%s(i,k) = ptend_sum%s(i,k) + ptend%s(i,k)
             end do
          end do
          do i = 1, ncol
             ptend_sum%hflux_srf(i) = ptend_sum%hflux_srf(i) + ptend%hflux_srf(i)
             ptend_sum%hflux_top(i) = ptend_sum%hflux_top(i) + ptend%hflux_top(i)
          end do
       else
          used_codon = .true.
       end if
    end if

    if (any(ptend%lq(:))) then

       if (.not. allocated(ptend_sum%q)) then
          allocate(ptend_sum%q(psetcols,pver,pcnst), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%q')
          ptend_sum%q=0.0_r8

          allocate(ptend_sum%cflx_srf(psetcols,pcnst), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%cflx_srf')
          ptend_sum%cflx_srf=0.0_r8

          allocate(ptend_sum%cflx_top(psetcols,pcnst), stat=ierr)
          if ( ierr /= 0 ) call endrun('physics_ptend_sum error: allocation error for ptend_sum%cflx_top')
          ptend_sum%cflx_top=0.0_r8
       end if

       do m = 1, pcnst
          if(ptend%lq(m)) then
             ptend_sum%lq(m) = .true.
             lq_mask(m) = 1_c_int64_t
             if (use_native_zero_impl) then
                do k = ptend%top_level, ptend%bot_level
                   do i = 1,ncol
                      ptend_sum%q(i,k,m) = ptend_sum%q(i,k,m) + ptend%q(i,k,m)
                   end do
                end do
                do i = 1,ncol
                   ptend_sum%cflx_srf(i,m) = ptend_sum%cflx_srf(i,m) + ptend%cflx_srf(i,m)
                   ptend_sum%cflx_top(i,m) = ptend_sum%cflx_top(i,m) + ptend%cflx_top(i,m)
                end do
             else
                used_codon = .true.
             end if
          end if
       end do

    end if

    if (used_codon) then
       src_s_p = c_null_ptr; dst_s_p = c_null_ptr
       src_u_p = c_null_ptr; dst_u_p = c_null_ptr
       src_v_p = c_null_ptr; dst_v_p = c_null_ptr
       src_q_p = c_null_ptr; dst_q_p = c_null_ptr
       src_hflux_srf_p = c_null_ptr; dst_hflux_srf_p = c_null_ptr
       src_hflux_top_p = c_null_ptr; dst_hflux_top_p = c_null_ptr
       src_taux_srf_p = c_null_ptr; dst_taux_srf_p = c_null_ptr
       src_taux_top_p = c_null_ptr; dst_taux_top_p = c_null_ptr
       src_tauy_srf_p = c_null_ptr; dst_tauy_srf_p = c_null_ptr
       src_tauy_top_p = c_null_ptr; dst_tauy_top_p = c_null_ptr
       src_cflx_srf_p = c_null_ptr; dst_cflx_srf_p = c_null_ptr
       src_cflx_top_p = c_null_ptr; dst_cflx_top_p = c_null_ptr

       if (ptend%ls) then
          src_s_p = c_loc(ptend%s); dst_s_p = c_loc(ptend_sum%s)
          src_hflux_srf_p = c_loc(ptend%hflux_srf); dst_hflux_srf_p = c_loc(ptend_sum%hflux_srf)
          src_hflux_top_p = c_loc(ptend%hflux_top); dst_hflux_top_p = c_loc(ptend_sum%hflux_top)
       end if
       if (ptend%lu) then
          src_u_p = c_loc(ptend%u); dst_u_p = c_loc(ptend_sum%u)
          src_taux_srf_p = c_loc(ptend%taux_srf); dst_taux_srf_p = c_loc(ptend_sum%taux_srf)
          src_taux_top_p = c_loc(ptend%taux_top); dst_taux_top_p = c_loc(ptend_sum%taux_top)
       end if
       if (ptend%lv) then
          src_v_p = c_loc(ptend%v); dst_v_p = c_loc(ptend_sum%v)
          src_tauy_srf_p = c_loc(ptend%tauy_srf); dst_tauy_srf_p = c_loc(ptend_sum%tauy_srf)
          src_tauy_top_p = c_loc(ptend%tauy_top); dst_tauy_top_p = c_loc(ptend_sum%tauy_top)
       end if
       if (any(ptend%lq(:))) then
          src_q_p = c_loc(ptend%q); dst_q_p = c_loc(ptend_sum%q)
          src_cflx_srf_p = c_loc(ptend%cflx_srf); dst_cflx_srf_p = c_loc(ptend_sum%cflx_srf)
          src_cflx_top_p = c_loc(ptend%cflx_top); dst_cflx_top_p = c_loc(ptend_sum%cflx_top)
       end if

       call physics_ptend_sum_codon_raw(int(ncol, c_int64_t), int(psetcols, c_int64_t), &
            int(pver, c_int64_t), int(pcnst, c_int64_t), int(ptend%top_level, c_int64_t), &
            int(ptend%bot_level, c_int64_t), merge(1_c_int64_t, 0_c_int64_t, ptend%ls), &
            merge(1_c_int64_t, 0_c_int64_t, ptend%lu), merge(1_c_int64_t, 0_c_int64_t, ptend%lv), &
            c_loc(lq_mask), src_s_p, dst_s_p, src_u_p, dst_u_p, src_v_p, dst_v_p, src_q_p, dst_q_p, &
            src_hflux_srf_p, dst_hflux_srf_p, src_hflux_top_p, dst_hflux_top_p, src_taux_srf_p, &
            dst_taux_srf_p, src_taux_top_p, dst_taux_top_p, src_tauy_srf_p, dst_tauy_srf_p, &
            src_tauy_top_p, dst_tauy_top_p, src_cflx_srf_p, dst_cflx_srf_p, src_cflx_top_p, &
            dst_cflx_top_p)
       call physics_types_zero_proof_once()
       call physics_types_log_direct(physics_ptend_sum_logged, &
            'physics_ptend_sum direct = codon; allocation/flag setup native derived-type boundary')
    end if

  end subroutine physics_ptend_sum

!===============================================================================

  subroutine physics_ptend_scale(ptend, fac, ncol)
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr, c_null_ptr
!-----------------------------------------------------------------------
! Scale ptend fields for ptend logical flags = .true.
! Where ptend logical flags = .false, don't change ptend.
!
! Assumes that input ptend is valid (e.g. that
! ptend%lu .eqv. allocated(ptend%u)), and therefore
! does not check allocation status of individual arrays.
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), target, intent(inout)  :: ptend   ! Incoming ptend
    real(r8), intent(in) :: fac                    ! Factor to multiply ptend by.
    integer, intent(in)                 :: ncol    ! number of columns

!---------------------------Local storage-------------------------------
    integer :: m                                   ! constituent index
    integer(c_int64_t), target :: lq_mask(pcnst)
    type(c_ptr) :: s_p, u_p, v_p, q_p, hflux_srf_p, hflux_top_p
    type(c_ptr) :: taux_srf_p, taux_top_p, tauy_srf_p, tauy_top_p
    type(c_ptr) :: cflx_srf_p, cflx_top_p

!-----------------------------------------------------------------------

    call physics_types_zero_select_impl()
    lq_mask = 0_c_int64_t

    if (use_native_zero_impl) then

! Update u,v fields
       if (ptend%lu) &
            call multiply_tendency(ptend%u, &
            ptend%taux_srf, ptend%taux_top)

       if (ptend%lv) &
            call multiply_tendency(ptend%v, &
            ptend%tauy_srf, ptend%tauy_top)

! Heat
       if (ptend%ls) &
            call multiply_tendency(ptend%s, &
            ptend%hflux_srf, ptend%hflux_top)

! Update constituents
       do m = 1, pcnst
          if (ptend%lq(m)) &
               call multiply_tendency(ptend%q(:,:,m), &
               ptend%cflx_srf(:,m), ptend%cflx_top(:,m))
       end do

    else

       do m = 1, pcnst
          if (ptend%lq(m)) then
             lq_mask(m) = 1_c_int64_t
          end if
       end do

       s_p = c_null_ptr; u_p = c_null_ptr; v_p = c_null_ptr; q_p = c_null_ptr
       hflux_srf_p = c_null_ptr; hflux_top_p = c_null_ptr
       taux_srf_p = c_null_ptr; taux_top_p = c_null_ptr
       tauy_srf_p = c_null_ptr; tauy_top_p = c_null_ptr
       cflx_srf_p = c_null_ptr; cflx_top_p = c_null_ptr

       if (ptend%ls) then
          s_p = c_loc(ptend%s)
          hflux_srf_p = c_loc(ptend%hflux_srf)
          hflux_top_p = c_loc(ptend%hflux_top)
       end if
       if (ptend%lu) then
          u_p = c_loc(ptend%u)
          taux_srf_p = c_loc(ptend%taux_srf)
          taux_top_p = c_loc(ptend%taux_top)
       end if
       if (ptend%lv) then
          v_p = c_loc(ptend%v)
          tauy_srf_p = c_loc(ptend%tauy_srf)
          tauy_top_p = c_loc(ptend%tauy_top)
       end if
       if (any(ptend%lq(:))) then
          q_p = c_loc(ptend%q)
          cflx_srf_p = c_loc(ptend%cflx_srf)
          cflx_top_p = c_loc(ptend%cflx_top)
       end if

       call physics_ptend_scale_codon_raw(int(ncol, c_int64_t), int(ptend%psetcols, c_int64_t), &
            int(pver, c_int64_t), int(pcnst, c_int64_t), int(ptend%top_level, c_int64_t), &
            int(ptend%bot_level, c_int64_t), real(fac, c_double), &
            merge(1_c_int64_t, 0_c_int64_t, ptend%ls), merge(1_c_int64_t, 0_c_int64_t, ptend%lu), &
            merge(1_c_int64_t, 0_c_int64_t, ptend%lv), c_loc(lq_mask), s_p, u_p, v_p, q_p, &
            hflux_srf_p, hflux_top_p, taux_srf_p, taux_top_p, tauy_srf_p, tauy_top_p, &
            cflx_srf_p, cflx_top_p)

       call physics_types_zero_proof_once()
       call physics_types_log_direct(physics_ptend_scale_logged, 'physics_ptend_scale direct = codon')
       call physics_types_log_direct(multiply_tendency_logged, &
            'multiply_tendency direct = codon via physics_ptend_scale_shell_codon; ' // &
            'native internal subroutine used only by native selector')
    end if


  contains

    subroutine multiply_tendency(tend_arr, flx_srf, flx_top)
      real(r8), intent(inout) :: tend_arr(:,:) ! Tendency array (pcols, plev)
      real(r8), intent(inout) :: flx_srf(:)    ! Surface flux (or stress)
      real(r8), intent(inout) :: flx_top(:)    ! Top-of-model flux (or stress)

      integer :: k

      do k = ptend%top_level, ptend%bot_level
         tend_arr(:ncol,k) = tend_arr(:ncol,k) * fac
      end do
      flx_srf(:ncol) = flx_srf(:ncol) * fac
      flx_top(:ncol) = flx_top(:ncol) * fac

    end subroutine multiply_tendency

  end subroutine physics_ptend_scale

!===============================================================================

subroutine physics_ptend_copy(ptend, ptend_cp)
   use iso_c_binding, only: c_int64_t, c_loc

   !-----------------------------------------------------------------------
   ! Copy a physics_ptend object.  Allocate ptend_cp internally before copy.
   !-----------------------------------------------------------------------

   type(physics_ptend), intent(in)    :: ptend    ! ptend source
   type(physics_ptend), intent(out)   :: ptend_cp ! copy of ptend

   !-----------------------------------------------------------------------

   ptend_cp%name      = ptend%name

   ptend_cp%ls = ptend%ls
   ptend_cp%lu = ptend%lu
   ptend_cp%lv = ptend%lv
   ptend_cp%lq = ptend%lq

   call physics_ptend_alloc(ptend_cp, ptend%psetcols)

   ptend_cp%top_level = ptend%top_level
   ptend_cp%bot_level = ptend%bot_level

   call physics_types_zero_select_impl()

   if (use_native_zero_impl) then

      if (ptend_cp%ls) then
         ptend_cp%s = ptend%s
         ptend_cp%hflux_srf = ptend%hflux_srf
         ptend_cp%hflux_top = ptend%hflux_top
      end if

      if (ptend_cp%lu) then
         ptend_cp%u = ptend%u
         ptend_cp%taux_srf  = ptend%taux_srf
         ptend_cp%taux_top  = ptend%taux_top
      end if

      if (ptend_cp%lv) then
         ptend_cp%v = ptend%v
         ptend_cp%tauy_srf  = ptend%tauy_srf
         ptend_cp%tauy_top  = ptend%tauy_top
      end if

      if (any(ptend_cp%lq(:))) then
         ptend_cp%q = ptend%q
         ptend_cp%cflx_srf  = ptend%cflx_srf
         ptend_cp%cflx_top  = ptend%cflx_top
      end if

   else

      if (ptend_cp%ls) then
         call physics_ptend_copy_s_codon(ptend%psetcols, ptend%s, ptend%hflux_srf, ptend%hflux_top, &
              ptend_cp%s, ptend_cp%hflux_srf, ptend_cp%hflux_top)
      end if

      if (ptend_cp%lu) then
         call physics_ptend_copy_u_codon(ptend%psetcols, ptend%u, ptend%taux_srf, ptend%taux_top, &
              ptend_cp%u, ptend_cp%taux_srf, ptend_cp%taux_top)
      end if

      if (ptend_cp%lv) then
         call physics_ptend_copy_v_codon(ptend%psetcols, ptend%v, ptend%tauy_srf, ptend%tauy_top, &
              ptend_cp%v, ptend_cp%tauy_srf, ptend_cp%tauy_top)
      end if

      if (any(ptend_cp%lq(:))) then
         call physics_ptend_copy_q_codon(ptend%psetcols, ptend%q, ptend%cflx_srf, ptend%cflx_top, &
              ptend_cp%q, ptend_cp%cflx_srf, ptend_cp%cflx_top)
      end if

      call physics_types_zero_proof_once()
   end if

end subroutine physics_ptend_copy

!===============================================================================

  subroutine physics_ptend_reset(ptend)
    use iso_c_binding, only: c_int64_t, c_loc, c_ptr, c_null_ptr
!-----------------------------------------------------------------------
! Reset the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), target, intent(inout)  :: ptend   ! Parameterization tendencies
!-----------------------------------------------------------------------
    integer :: m             ! Index for constiuent
    integer(c_int64_t), target :: lq_mask(pcnst)
    type(c_ptr) :: s_p, u_p, v_p, q_p, hflux_srf_p, hflux_top_p
    type(c_ptr) :: taux_srf_p, taux_top_p, tauy_srf_p, tauy_top_p
    type(c_ptr) :: cflx_srf_p, cflx_top_p
!-----------------------------------------------------------------------

    call physics_types_zero_select_impl()
    lq_mask = 0_c_int64_t

    if (use_native_zero_impl) then
       if(ptend%ls) then
          ptend%s = 0._r8
          ptend%hflux_srf = 0._r8
          ptend%hflux_top = 0._r8
       endif
       if(ptend%lu) then
          ptend%u = 0._r8
          ptend%taux_srf = 0._r8
          ptend%taux_top = 0._r8
       endif
       if(ptend%lv) then
          ptend%v = 0._r8
          ptend%tauy_srf = 0._r8
          ptend%tauy_top = 0._r8
       endif
       if(any (ptend%lq(:))) then
          ptend%q = 0._r8
          ptend%cflx_srf = 0._r8
          ptend%cflx_top = 0._r8
       end if
    else
       do m = 1, pcnst
          if (ptend%lq(m)) lq_mask(m) = 1_c_int64_t
       end do

       s_p = c_null_ptr; u_p = c_null_ptr; v_p = c_null_ptr; q_p = c_null_ptr
       hflux_srf_p = c_null_ptr; hflux_top_p = c_null_ptr
       taux_srf_p = c_null_ptr; taux_top_p = c_null_ptr
       tauy_srf_p = c_null_ptr; tauy_top_p = c_null_ptr
       cflx_srf_p = c_null_ptr; cflx_top_p = c_null_ptr

       if (ptend%ls) then
          s_p = c_loc(ptend%s)
          hflux_srf_p = c_loc(ptend%hflux_srf)
          hflux_top_p = c_loc(ptend%hflux_top)
       end if
       if (ptend%lu) then
          u_p = c_loc(ptend%u)
          taux_srf_p = c_loc(ptend%taux_srf)
          taux_top_p = c_loc(ptend%taux_top)
       end if
       if (ptend%lv) then
          v_p = c_loc(ptend%v)
          tauy_srf_p = c_loc(ptend%tauy_srf)
          tauy_top_p = c_loc(ptend%tauy_top)
       end if
       if (any(ptend%lq(:))) then
          q_p = c_loc(ptend%q)
          cflx_srf_p = c_loc(ptend%cflx_srf)
          cflx_top_p = c_loc(ptend%cflx_top)
       end if

       call physics_ptend_reset_codon_raw(int(ptend%psetcols, c_int64_t), int(pver, c_int64_t), &
            int(pcnst, c_int64_t), merge(1_c_int64_t, 0_c_int64_t, ptend%ls), &
            merge(1_c_int64_t, 0_c_int64_t, ptend%lu), merge(1_c_int64_t, 0_c_int64_t, ptend%lv), &
            c_loc(lq_mask), s_p, u_p, v_p, q_p, hflux_srf_p, hflux_top_p, taux_srf_p, taux_top_p, &
            tauy_srf_p, tauy_top_p, cflx_srf_p, cflx_top_p)
       call physics_types_zero_proof_once()
       call physics_types_log_direct(physics_ptend_reset_logged, 'physics_ptend_reset direct = codon')
    end if

    ptend%top_level = 1
    ptend%bot_level = pver

    return
  end subroutine physics_ptend_reset

!===============================================================================
  subroutine physics_ptend_init(ptend, psetcols, name, ls, lu, lv, lq)
    use iso_c_binding, only: c_int64_t
!-----------------------------------------------------------------------
! Allocate the fields in the structure which are specified.
! Initialize the parameterization tendency structure to "empty"
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    type(physics_ptend), intent(out)    :: ptend    ! Parameterization tendencies
    integer, intent(in)                 :: psetcols ! maximum number of columns
    character(len=*)                    :: name     ! optional name of parameterization which produced tendencies.
    logical, optional                   :: ls       ! if true, then fields to support dsdt are allocated
    logical, optional                   :: lu       ! if true, then fields to support dudt are allocated
    logical, optional                   :: lv       ! if true, then fields to support dvdt are allocated
    logical, dimension(pcnst),optional  :: lq       ! if true, then fields to support dqdt are allocated

    integer :: init_policy
    logical :: ls_flag, lu_flag, lv_flag

!-----------------------------------------------------------------------

    call physics_types_zero_select_impl()

    if (allocated(ptend%s)) then
       call endrun(' physics_ptend_init: ptend should not be allocated before calling this routine')
    end if

    ptend%name     = name
    ptend%psetcols =  psetcols
    ls_flag = .false.
    lu_flag = .false.
    lv_flag = .false.
    if (present(ls)) ls_flag = ls
    if (present(lu)) lu_flag = lu
    if (present(lv)) lv_flag = lv

    ! If no fields being stored, initialize all values to appropriate nulls and return
    if (use_native_zero_impl) then
       if (.not. present(ls) .and. .not. present(lu) .and. .not. present(lv) .and. .not. present(lq) ) then
          ptend%ls       = .false.
          ptend%lu       = .false.
          ptend%lv       = .false.
          ptend%lq(:)    = .false.
          ptend%top_level = 1
          ptend%bot_level = pver
          return
       end if

       if (present(ls)) then
          ptend%ls = ls
       else
          ptend%ls = .false.
       end if

       if (present(lu)) then
          ptend%lu = lu
       else
          ptend%lu = .false.
       end if

       if (present(lv)) then
          ptend%lv = lv
       else
          ptend%lv = .false.
       end if
    else
       init_policy = physics_ptend_init_codon(present(ls), ls_flag, present(lu), lu_flag, &
            present(lv), lv_flag, present(lq))
       call physics_types_zero_proof_once()
       call physics_types_log_direct(physics_ptend_init_logged, &
            'physics_ptend_init direct = codon; optional flag/no-field policy direct; allocation/name/lq copy native Fortran boundary')

       if (btest(init_policy, 3)) then
          ptend%ls       = .false.
          ptend%lu       = .false.
          ptend%lv       = .false.
          ptend%lq(:)    = .false.
          ptend%top_level = 1
          ptend%bot_level = pver
          return
       end if

       ptend%ls = btest(init_policy, 0)
       ptend%lu = btest(init_policy, 1)
       ptend%lv = btest(init_policy, 2)
    end if

    if (present(lq)) then
       ptend%lq(:) = lq(:)
    else
       ptend%lq(:) = .false.
    end if

    call physics_ptend_alloc(ptend, psetcols)

    call physics_ptend_reset(ptend)

    return
  end subroutine physics_ptend_init

!===============================================================================

  subroutine physics_state_set_grid(lchnk, phys_state)
    use iso_c_binding, only: c_int64_t
!-----------------------------------------------------------------------
! Set the grid components of the physics_state object
!-----------------------------------------------------------------------

    integer,             intent(in)    :: lchnk
    type(physics_state), intent(inout) :: phys_state

    ! local variables
    integer  :: i, ncol
    real(r8), target :: rlon(pcols)
    real(r8), target :: rlat(pcols)

    !-----------------------------------------------------------------------
    ! get_ncols_p requires a state which does not have sub-columns
    if (phys_state%psetcols .ne. pcols) then
       call endrun('physics_state_set_grid: cannot pass in a state which has sub-columns')
    end if

    ncol = get_ncols_p(lchnk)

    if(ncol<=0) then
       write(iulog,*) lchnk, ncol
       call endrun('physics_state_set_grid')
    end if

    call get_rlon_all_p(lchnk, ncol, rlon)
    call get_rlat_all_p(lchnk, ncol, rlat)
    phys_state%ncol  = ncol
    phys_state%lchnk = lchnk
    call physics_types_zero_select_impl()
    if (use_native_zero_impl) then
       do i=1,ncol
          phys_state%lat(i) = rlat(i)
          phys_state%lon(i) = rlon(i)
       end do
       call init_geo_unique(phys_state,ncol)
    else
       call physics_state_set_grid_codon(ncol, phys_state%psetcols, rlat, rlon, &
            phys_state%lat, phys_state%lon, phys_state%ulat, phys_state%ulon, &
            phys_state%latmapback, phys_state%lonmapback, phys_state%ulatcnt, phys_state%uloncnt)
       call physics_types_zero_proof_once()
       call physics_types_log_direct(physics_state_set_grid_logged, &
            'physics_state_set_grid direct = codon; init_geo_unique direct; get_gcol_all_p native CAM API island')
    end if

  end subroutine physics_state_set_grid


  subroutine init_geo_unique(phys_state,ncol)
    integer,             intent(in)    :: ncol
    type(physics_state), intent(inout) :: phys_state
    logical :: match
    integer :: i, j, ulatcnt, uloncnt

    call physics_types_zero_select_impl()

    if (use_native_zero_impl) then
       phys_state%ulat=-999._r8
       phys_state%ulon=-999._r8
       phys_state%latmapback=0
       phys_state%lonmapback=0
       match=.false.
       ulatcnt=0
       uloncnt=0
       match=.false.
       do i=1,ncol
          do j=1,ulatcnt
             if(phys_state%lat(i) .eq. phys_state%ulat(j)) then
                match=.true.
                phys_state%latmapback(i)=j
             end if
          end do
          if(.not. match) then
             ulatcnt=ulatcnt+1
             phys_state%ulat(ulatcnt)=phys_state%lat(i)
             phys_state%latmapback(i)=ulatcnt
          end if

          match=.false.
          do j=1,uloncnt
             if(phys_state%lon(i) .eq. phys_state%ulon(j)) then
                match=.true.
                phys_state%lonmapback(i)=j
             end if
          end do
          if(.not. match) then
             uloncnt=uloncnt+1
             phys_state%ulon(uloncnt)=phys_state%lon(i)
             phys_state%lonmapback(i)=uloncnt
          end if
          match=.false.

       end do
       phys_state%uloncnt=uloncnt
       phys_state%ulatcnt=ulatcnt
    else
       call init_geo_unique_codon(ncol, phys_state%psetcols, phys_state%lat, phys_state%lon, &
            phys_state%ulat, phys_state%ulon, phys_state%latmapback, phys_state%lonmapback, ulatcnt, uloncnt)
       phys_state%uloncnt=uloncnt
       phys_state%ulatcnt=ulatcnt
       call physics_types_zero_proof_once()
       call physics_types_log_direct(init_geo_unique_logged, 'init_geo_unique direct = codon')
    end if

    call get_gcol_all_p(phys_state%lchnk,pcols,phys_state%cid)


  end subroutine init_geo_unique

!===============================================================================
  subroutine physics_dme_adjust(state, tend, qini, dt)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Adjust the dry mass in each layer back to the value of physics input state
    !
    ! Method: Conserve the integrated mass, momentum and total energy in each layer
    !         by scaling the specific mass of consituents, specific momentum (velocity)
    !         and specific total energy by the relative change in layer mass. Solve for
    !         the new temperature by subtracting the new kinetic energy from total energy
    !         and inverting the hydrostatic equation
    !
    !         The mass in each layer is modified, changing the relationship of the layer
    !         interfaces and midpoints to the surface pressure. The result is no longer in
    !         the original hybrid coordinate.
    !
    !         This procedure cannot be applied to the "eul" or "sld" dycores because they
    !         require the hybrid coordinate.
    !
    ! Author: Byron Boville

    ! !REVISION HISTORY:
    !   03.03.28  Boville    Created, partly from code by Lin in p_d_adjust
    !
    !-----------------------------------------------------------------------

    use constituents, only : cnst_get_type_byind
    use iso_c_binding, only: c_int64_t

    implicit none
    !
    ! Arguments
    !
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    real(r8),            intent(in   ) :: qini(pcols,pver)    ! initial specific humidity
    real(r8),            intent(in   ) :: dt                  ! model physics timestep
    !
    !---------------------------Local workspace-----------------------------
    !
    integer(c_int64_t) :: active_c
    integer  :: lchnk         ! chunk identifier
    integer  :: ncol          ! number of atmospheric columns
    integer  :: i,k,m         ! Longitude, level indices
    real(r8) :: fdq(pcols)    ! mass adjustment factor
    real(r8) :: te(pcols)     ! total energy in a layer
    real(r8) :: utmp(pcols)   ! temp variable for recalculating the initial u values
    real(r8) :: vtmp(pcols)   ! temp variable for recalculating the initial v values

    real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer
    !
    !-----------------------------------------------------------------------
    ! verify that the dycore is FV
    call physics_types_zero_select_impl()
    if (use_native_zero_impl) then
       if (.not. dycore_is('LR') ) return
    else
       active_c = physics_dme_adjust_codon(merge(1_c_int64_t, 0_c_int64_t, dycore_is('LR')))
       call physics_types_zero_proof_once()
       call physics_types_log_direct(physics_dme_adjust_logged, 'physics_dme_adjust direct = codon')
       if (active_c == 0_c_int64_t) return
    end if

    if (state%psetcols .ne. pcols) then
       call endrun('physics_dme_adjust: cannot pass in a state which has sub-columns')
    end if

    lchnk = state%lchnk
    ncol  = state%ncol

    ! adjust dry mass in each layer back to input value, while conserving
    ! constituents, momentum, and total energy
    do k = 1, pver

       ! adjusment factor is just change in water vapor
       fdq(:ncol) = 1._r8 + state%q(:ncol,k,1) - qini(:ncol,k)

       ! adjust constituents to conserve mass in each layer
       do m = 1, pcnst
          state%q(:ncol,k,m) = state%q(:ncol,k,m) / fdq(:ncol)
       end do

       if (adjust_te) then
          ! compute specific total energy of unadjusted state (J/kg)
          te(:ncol) = state%s(:ncol,k) + 0.5_r8*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2)

          ! recompute initial u,v from the new values and the tendencies
          utmp(:ncol) = state%u(:ncol,k) - dt * tend%dudt(:ncol,k)
          vtmp(:ncol) = state%v(:ncol,k) - dt * tend%dvdt(:ncol,k)
          ! adjust specific total energy and specific momentum (velocity) to conserve each
          te     (:ncol)   = te     (:ncol)     / fdq(:ncol)
          state%u(:ncol,k) = state%u(:ncol,k  ) / fdq(:ncol)
          state%v(:ncol,k) = state%v(:ncol,k  ) / fdq(:ncol)
          ! compute adjusted u,v tendencies
          tend%dudt(:ncol,k) = (state%u(:ncol,k) - utmp(:ncol)) / dt
          tend%dvdt(:ncol,k) = (state%v(:ncol,k) - vtmp(:ncol)) / dt

          ! compute adjusted static energy
          state%s(:ncol,k) = te(:ncol) - 0.5_r8*(state%u(:ncol,k)**2 + state%v(:ncol,k)**2)
       end if

! compute new total pressure variables
       state%pdel  (:ncol,k  ) = state%pdel(:ncol,k  ) * fdq(:ncol)
       state%pint  (:ncol,k+1) = state%pint(:ncol,k  ) + state%pdel(:ncol,k)
       state%lnpint(:ncol,k+1) = log(state%pint(:ncol,k+1))
       state%rpdel (:ncol,k  ) = 1._r8/ state%pdel(:ncol,k  )
    end do

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
      zvirv(:,:) = shr_const_rwv / rairv(:,:,state%lchnk) - 1._r8
    else
      zvirv(:,:) = zvir
    endif

! compute new T,z from new s,q,dp
    if (adjust_te) then
       call geopotential_dse(state%lnpint, state%lnpmid, state%pint,  &
            state%pmid  , state%pdel    , state%rpdel,  &
            state%s     , state%q(:,:,1), state%phis , rairv(:,:,state%lchnk), &
            gravit, cpairv(:,:,state%lchnk), zvirv, &
            state%t     , state%zi      , state%zm   , ncol)
    end if

  end subroutine physics_dme_adjust
!-----------------------------------------------------------------------

!===============================================================================
  subroutine physics_state_copy(state_in, state_out)

    use ppgrid,       only: pver, pverp
    use constituents, only: pcnst

    implicit none

    !
    ! Arguments
    !
    type(physics_state), intent(in)    :: state_in
    type(physics_state), intent(out)   :: state_out

    !
    ! Local variables
    !
    integer i, k, m, ncol

    ! Allocate state_out with same subcol dimension as state_in
    call physics_state_alloc ( state_out, state_in%lchnk, state_in%psetcols)

    ncol = state_in%ncol

    state_out%psetcols = state_in%psetcols
    state_out%ngrdcol  = state_in%ngrdcol
    state_out%lchnk    = state_in%lchnk
    state_out%ncol     = state_in%ncol
    state_out%count    = state_in%count

    call physics_types_zero_select_impl()

    if (use_native_zero_impl) then
       do i = 1, ncol
          state_out%lat(i)    = state_in%lat(i)
          state_out%lon(i)    = state_in%lon(i)
          state_out%ps(i)     = state_in%ps(i)
          state_out%phis(i)   = state_in%phis(i)
          state_out%te_ini(i) = state_in%te_ini(i)
          state_out%te_cur(i) = state_in%te_cur(i)
          state_out%tw_ini(i) = state_in%tw_ini(i)
          state_out%tw_cur(i) = state_in%tw_cur(i)
       end do

       do k = 1, pver
          do i = 1, ncol
             state_out%t(i,k)         = state_in%t(i,k)
             state_out%u(i,k)         = state_in%u(i,k)
             state_out%v(i,k)         = state_in%v(i,k)
             state_out%s(i,k)         = state_in%s(i,k)
             state_out%omega(i,k)     = state_in%omega(i,k)
             state_out%pmid(i,k)      = state_in%pmid(i,k)
             state_out%pdel(i,k)      = state_in%pdel(i,k)
             state_out%rpdel(i,k)     = state_in%rpdel(i,k)
             state_out%lnpmid(i,k)    = state_in%lnpmid(i,k)
             state_out%exner(i,k)     = state_in%exner(i,k)
             state_out%zm(i,k)        = state_in%zm(i,k)
          end do
       end do

       do k = 1, pverp
          do i = 1, ncol
             state_out%pint(i,k)      = state_in%pint(i,k)
             state_out%lnpint(i,k)    = state_in%lnpint(i,k)
             state_out%zi(i,k)        = state_in% zi(i,k)
          end do
       end do


          do i = 1, ncol
             state_out%psdry(i)  = state_in%psdry(i)
          end do
          do k = 1, pver
             do i = 1, ncol
                state_out%lnpmiddry(i,k) = state_in%lnpmiddry(i,k)
                state_out%pmiddry(i,k)   = state_in%pmiddry(i,k)
                state_out%pdeldry(i,k)   = state_in%pdeldry(i,k)
                state_out%rpdeldry(i,k)  = state_in%rpdeldry(i,k)
             end do
          end do
          do k = 1, pverp
             do i = 1, ncol
                state_out%pintdry(i,k)   = state_in%pintdry(i,k)
                state_out%lnpintdry(i,k) = state_in%lnpintdry(i,k)
             end do
          end do

       do m = 1, pcnst
          do k = 1, pver
             do i = 1, ncol
                state_out%q(i,k,m) = state_in%q(i,k,m)
             end do
          end do
       end do
    else
       call physics_state_copy_codon(state_in, state_out, ncol)

       call physics_types_zero_proof_once()
       call physics_types_log_direct(physics_state_copy_logged, &
            'physics_state_copy direct = codon; state field copy direct; allocation/setup native Fortran boundary')
    end if

  end subroutine physics_state_copy
!===============================================================================

  subroutine physics_tend_init(tend)
    use iso_c_binding, only: c_int64_t, c_loc

    implicit none

    !
    ! Arguments
    !
    type(physics_tend), intent(inout) :: tend

    !
    ! Local variables
    !

    if (.not. allocated(tend%dtdt)) then
       call endrun('physics_tend_init: tend must be allocated before it can be initialized')
    end if

    call physics_types_zero_select_impl()

    if (use_native_zero_impl) then
       tend%dtdt    = 0._r8
       tend%dudt    = 0._r8
       tend%dvdt    = 0._r8
       tend%flx_net = 0._r8
       tend%te_tnd  = 0._r8
       tend%tw_tnd  = 0._r8
    else
       call physics_tend_init_codon(tend%psetcols, tend%dtdt, tend%dudt, tend%dvdt, tend%flx_net, &
            tend%te_tnd, tend%tw_tnd)
       call physics_types_zero_proof_once()
    end if

end subroutine physics_tend_init

!===============================================================================

subroutine set_state_pdry (state,pdeld_calc)

  use ppgrid,  only: pver
  use pmgrid,  only: plev, plevp
  implicit none

  type(physics_state), intent(inout) :: state
  logical, optional, intent(in) :: pdeld_calc    !  .true. do calculate pdeld [default]
                                                 !  .false. don't calculate pdeld
  integer ncol
  integer i, k
  logical do_pdeld_calc

  if ( present(pdeld_calc) ) then
     do_pdeld_calc = pdeld_calc
  else
     do_pdeld_calc = .true.
  endif

  ncol = state%ncol

  call physics_types_zero_select_impl()

  if (use_native_zero_impl) then
     state%psdry(:ncol) = state%pint(:ncol,1)
     state%pintdry(:ncol,1) = state%pint(:ncol,1)

     if (do_pdeld_calc)  then
        do k = 1, pver
           state%pdeldry(:ncol,k) = state%pdel(:ncol,k)*(1._r8-state%q(:ncol,k,1))
        end do
     endif
     do k = 1, pver
        state%pintdry(:ncol,k+1) = state%pintdry(:ncol,k)+state%pdeldry(:ncol,k)
        state%pmiddry(:ncol,k) = (state%pintdry(:ncol,k+1)+state%pintdry(:ncol,k))/2._r8
        state%psdry(:ncol) = state%psdry(:ncol) + state%pdeldry(:ncol,k)
     end do

     state%rpdeldry(:ncol,:) = 1._r8/state%pdeldry(:ncol,:)
     state%lnpmiddry(:ncol,:) = log(state%pmiddry(:ncol,:))
     state%lnpintdry(:ncol,:) = log(state%pintdry(:ncol,:))
  else
     call set_state_pdry_codon(ncol, state%psetcols, do_pdeld_calc, state%psdry, state%pint, &
          state%pdel, state%q, state%pdeldry, state%pintdry, state%pmiddry, state%rpdeldry, &
          state%lnpmiddry, state%lnpintdry)
     call physics_types_zero_proof_once()
     call physics_types_log_direct(set_state_pdry_logged, 'set_state_pdry direct = codon')
  end if

end subroutine set_state_pdry

!===============================================================================

subroutine set_wet_to_dry (state)

  use constituents,  only: pcnst, cnst_type
  use iso_c_binding, only: c_int64_t

  type(physics_state), intent(inout) :: state

  integer m, ncol
  integer(c_int64_t), target :: dry_mask(pcnst)

  ncol = state%ncol
  dry_mask = 0_c_int64_t
  call physics_types_zero_select_impl()

  if (use_native_zero_impl) then
     do m = 1,pcnst
        if (cnst_type(m).eq.'dry') then
           state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdel(:ncol,:)/state%pdeldry(:ncol,:)
        endif
     end do
  else
     do m = 1,pcnst
        if (cnst_type(m).eq.'dry') then
           dry_mask(m) = 1_c_int64_t
        endif
     end do
     if (any(dry_mask /= 0_c_int64_t)) then
        call set_wet_to_dry_codon(ncol, state%psetcols, dry_mask, state%q, state%pdel, state%pdeldry)
        call physics_types_zero_proof_once()
        call physics_types_log_direct(set_wet_to_dry_logged, 'set_wet_to_dry direct = codon')
     end if
  end if

end subroutine set_wet_to_dry

!===============================================================================

subroutine set_dry_to_wet (state)

  use constituents,  only: pcnst, cnst_type
  use iso_c_binding, only: c_int64_t

  type(physics_state), intent(inout) :: state

  integer m, ncol
  integer(c_int64_t), target :: dry_mask(pcnst)

  ncol = state%ncol
  dry_mask = 0_c_int64_t
  call physics_types_zero_select_impl()

  if (use_native_zero_impl) then
     do m = 1,pcnst
        if (cnst_type(m).eq.'dry') then
           state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdeldry(:ncol,:)/state%pdel(:ncol,:)
        endif
     end do
  else
     do m = 1,pcnst
        if (cnst_type(m).eq.'dry') then
           dry_mask(m) = 1_c_int64_t
        endif
     end do
     if (any(dry_mask /= 0_c_int64_t)) then
        call set_dry_to_wet_codon(ncol, state%psetcols, dry_mask, state%q, state%pdeldry, state%pdel)
        call physics_types_zero_proof_once()
        call physics_types_log_direct(set_dry_to_wet_logged, 'set_dry_to_wet direct = codon')
     end if
  end if

end subroutine set_dry_to_wet

!===============================================================================

subroutine physics_state_alloc(state,lchnk,psetcols)

  use infnan, only : inf, assignment(=)

! allocate the individual state components

  type(physics_state), intent(inout) :: state
  integer,intent(in)                 :: lchnk

  integer, intent(in)                :: psetcols

  integer :: ierr=0, i
  real(r8) :: inf_local

  state%lchnk    = lchnk
  state%psetcols = psetcols
  state%ngrdcol  = get_ncols_p(lchnk)  ! Number of grid columns

  !----------------------------------
  ! Following variables will be overwritten by sub-column generator, if sub-columns are being used

  !  state%ncol - is initialized in physics_state_set_grid,  if not using sub-columns

  !----------------------------------

  allocate(state%lat(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lat')

  allocate(state%lon(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lon')

  allocate(state%ps(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%ps')

  allocate(state%psdry(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%psdry')

  allocate(state%phis(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%phis')

  allocate(state%ulat(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%ulat')

  allocate(state%ulon(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%ulon')

  allocate(state%t(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%t')

  allocate(state%u(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%u')

  allocate(state%v(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%v')

  allocate(state%s(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%s')

  allocate(state%omega(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%omega')

  allocate(state%pmid(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pmid')

  allocate(state%pmiddry(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pmiddry')

  allocate(state%pdel(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pdel')

  allocate(state%pdeldry(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pdeldry')

  allocate(state%rpdel(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%rpdel')

  allocate(state%rpdeldry(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%rpdeldry')

  allocate(state%lnpmid(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lnpmid')

  allocate(state%lnpmiddry(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lnpmiddry')

  allocate(state%exner(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%exner')

  allocate(state%zm(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%zm')

  allocate(state%q(psetcols,pver,pcnst), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%q')

  allocate(state%pint(psetcols,pver+1), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pint')

  allocate(state%pintdry(psetcols,pver+1), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pintdry')

  allocate(state%lnpint(psetcols,pver+1), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lnpint')

  allocate(state%lnpintdry(psetcols,pver+1), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lnpintdry')

  allocate(state%zi(psetcols,pver+1), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%zi')

  allocate(state%te_ini(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%te_ini')

  allocate(state%te_cur(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%te_cur')

  allocate(state%tw_ini(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%tw_ini')

  allocate(state%tw_cur(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%tw_cur')

  allocate(state%latmapback(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%latmapback')

  allocate(state%lonmapback(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lonmapback')

  allocate(state%cid(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%cid')

  inf_local = inf
  call physics_types_zero_select_impl()

  if (use_native_zero_impl) then
     state%lat(:) = inf
     state%lon(:) = inf
     state%ulat(:) = inf
     state%ulon(:) = inf
     state%ps(:) = inf
     state%psdry(:) = inf
     state%phis(:) = inf
     state%t(:,:) = inf
     state%u(:,:) = inf
     state%v(:,:) = inf
     state%s(:,:) = inf
     state%omega(:,:) = inf
     state%pmid(:,:) = inf
     state%pmiddry(:,:) = inf
     state%pdel(:,:) = inf
     state%pdeldry(:,:) = inf
     state%rpdel(:,:) = inf
     state%rpdeldry(:,:) = inf
     state%lnpmid(:,:) = inf
     state%lnpmiddry(:,:) = inf
     state%exner(:,:) = inf
     state%zm(:,:) = inf
     state%q(:,:,:) = inf

     state%pint(:,:) = inf
     state%pintdry(:,:) = inf
     state%lnpint(:,:) = inf
     state%lnpintdry(:,:) = inf
     state%zi(:,:) = inf

     state%te_ini(:) = inf
     state%te_cur(:) = inf
     state%tw_ini(:) = inf
     state%tw_cur(:) = inf
  else
     call physics_state_alloc_codon(state, psetcols, inf_local)
     call physics_types_zero_proof_once()
     call physics_types_log_direct(physics_state_alloc_logged, &
          'physics_state_alloc direct = codon; state inf-fill direct; allocation/error/derived-type setup native Fortran boundary')
  end if

end subroutine physics_state_alloc

!===============================================================================

subroutine physics_state_dealloc(state)

  use iso_c_binding, only: c_int64_t

! deallocate the individual state components

  type(physics_state), intent(inout) :: state
  integer                            :: ierr = 0

  call physics_types_zero_select_impl()
  if (.not. use_native_zero_impl) then
     if (physics_state_dealloc_codon_raw(4_c_int64_t) == 4_c_int64_t) then
        call physics_types_zero_proof_once()
        call physics_types_log_direct(physics_state_dealloc_logged, &
             'physics_state_dealloc direct = codon; deallocation selector/touch direct = codon; native deallocate/error island')
     end if
  end if

  deallocate(state%lat, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lat')

  deallocate(state%lon, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lon')

  deallocate(state%ps, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%ps')

  deallocate(state%psdry, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%psdry')

  deallocate(state%phis, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%phis')

  deallocate(state%ulat, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%ulat')

  deallocate(state%ulon, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%ulon')

  deallocate(state%t, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%t')

  deallocate(state%u, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%u')

  deallocate(state%v, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%v')

  deallocate(state%s, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%s')

  deallocate(state%omega, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%omega')

  deallocate(state%pmid, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pmid')

  deallocate(state%pmiddry, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pmiddry')

  deallocate(state%pdel, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pdel')

  deallocate(state%pdeldry, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pdeldry')

  deallocate(state%rpdel, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%rpdel')

  deallocate(state%rpdeldry, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%rpdeldry')

  deallocate(state%lnpmid, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lnpmid')

  deallocate(state%lnpmiddry, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lnpmiddry')

  deallocate(state%exner, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%exner')

  deallocate(state%zm, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%zm')

  deallocate(state%q, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%q')

  deallocate(state%pint, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pint')

  deallocate(state%pintdry, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pintdry')

  deallocate(state%lnpint, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lnpint')

  deallocate(state%lnpintdry, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lnpintdry')

  deallocate(state%zi, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%zi')

  deallocate(state%te_ini, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%te_ini')

  deallocate(state%te_cur, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%te_cur')

  deallocate(state%tw_ini, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%tw_ini')

  deallocate(state%tw_cur, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%tw_cur')

  deallocate(state%latmapback, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%latmapback')

  deallocate(state%lonmapback, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lonmapback')

  deallocate(state%cid, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%cid')

end subroutine physics_state_dealloc

!===============================================================================

subroutine physics_tend_alloc(tend,psetcols)

  use infnan, only : inf, assignment(=)
! allocate the individual tend components

  type(physics_tend), intent(inout)  :: tend

  integer, intent(in)                :: psetcols

  integer :: ierr = 0
  real(r8) :: inf_local

  tend%psetcols = psetcols

  allocate(tend%dtdt(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%dtdt')

  allocate(tend%dudt(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%dudt')

  allocate(tend%dvdt(psetcols,pver), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%dvdt')

  allocate(tend%flx_net(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%flx_net')

  allocate(tend%te_tnd(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%te_tnd')

  allocate(tend%tw_tnd(psetcols), stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%tw_tnd')

  inf_local = inf
  call physics_types_zero_select_impl()

  if (use_native_zero_impl) then
     tend%dtdt(:,:) = inf
     tend%dudt(:,:) = inf
     tend%dvdt(:,:) = inf
     tend%flx_net(:) = inf
     tend%te_tnd(:) = inf
     tend%tw_tnd(:) = inf
  else
     call physics_tend_alloc_codon(tend, psetcols, inf_local)
     call physics_types_zero_proof_once()
     call physics_types_log_direct(physics_tend_alloc_logged, &
          'physics_tend_alloc direct = codon; tend inf-fill direct; allocation/error/derived-type setup native Fortran boundary')
  end if

end subroutine physics_tend_alloc

!===============================================================================

subroutine physics_tend_dealloc(tend)

! deallocate the individual tend components

  type(physics_tend), intent(inout)  :: tend
  integer :: psetcols
  integer :: ierr = 0

  deallocate(tend%dtdt, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%dtdt')

  deallocate(tend%dudt, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%dudt')

  deallocate(tend%dvdt, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%dvdt')

  deallocate(tend%flx_net, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%flx_net')

  deallocate(tend%te_tnd, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%te_tnd')

  deallocate(tend%tw_tnd, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%tw_tnd')
end subroutine physics_tend_dealloc

!===============================================================================

subroutine physics_ptend_alloc(ptend,psetcols)

  use iso_c_binding, only: c_int64_t

! allocate the individual ptend components

  type(physics_ptend), intent(inout) :: ptend

  integer, intent(in)                :: psetcols

  integer :: ierr = 0

  call physics_types_zero_select_impl()
  if (.not. use_native_zero_impl) then
     if (physics_ptend_alloc_codon_raw(5_c_int64_t) == 5_c_int64_t) then
        call physics_types_zero_proof_once()
        call physics_types_log_direct(physics_ptend_alloc_logged, &
             'physics_ptend_alloc direct = codon; allocation selector/touch direct = codon; native allocate/error island')
     end if
  end if

  ptend%psetcols = psetcols

  if (ptend%ls) then
     allocate(ptend%s(psetcols,pver), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%s')

     allocate(ptend%hflux_srf(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%hflux_srf')

     allocate(ptend%hflux_top(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%hflux_top')
  end if

  if (ptend%lu) then
     allocate(ptend%u(psetcols,pver), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%u')

     allocate(ptend%taux_srf(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%taux_srf')

     allocate(ptend%taux_top(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%taux_top')
  end if

  if (ptend%lv) then
     allocate(ptend%v(psetcols,pver), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%v')

     allocate(ptend%tauy_srf(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%tauy_srf')

     allocate(ptend%tauy_top(psetcols), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%tauy_top')
  end if

  if (any(ptend%lq)) then
     allocate(ptend%q(psetcols,pver,pcnst), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%q')

     allocate(ptend%cflx_srf(psetcols,pcnst), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%cflx_srf')

     allocate(ptend%cflx_top(psetcols,pcnst), stat=ierr)
     if ( ierr /= 0 ) call endrun('physics_ptend_alloc error: allocation error for ptend%cflx_top')
  end if

end subroutine physics_ptend_alloc

!===============================================================================

subroutine physics_ptend_dealloc(ptend)

  use iso_c_binding, only: c_int64_t

! deallocate the individual ptend components

  type(physics_ptend), intent(inout) :: ptend
  integer :: ierr = 0

  call physics_types_zero_select_impl()
  if (.not. use_native_zero_impl) then
     if (physics_ptend_dealloc_codon_raw(6_c_int64_t) == 6_c_int64_t) then
        call physics_types_zero_proof_once()
        call physics_types_log_direct(physics_ptend_dealloc_logged, &
             'physics_ptend_dealloc direct = codon; deallocation selector/touch direct = codon; native deallocate/error island')
     end if
  end if

  ptend%psetcols = 0

  if (allocated(ptend%s)) deallocate(ptend%s, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%s')

  if (allocated(ptend%hflux_srf))   deallocate(ptend%hflux_srf, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%hflux_srf')

  if (allocated(ptend%hflux_top))  deallocate(ptend%hflux_top, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%hflux_top')

  if (allocated(ptend%u))   deallocate(ptend%u, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%u')

  if (allocated(ptend%taux_srf)) deallocate(ptend%taux_srf, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%taux_srf')

  if (allocated(ptend%taux_top))   deallocate(ptend%taux_top, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%taux_top')

  if (allocated(ptend%v)) deallocate(ptend%v, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%v')

  if (allocated(ptend%tauy_srf))   deallocate(ptend%tauy_srf, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%tauy_srf')

  if (allocated(ptend%tauy_top))   deallocate(ptend%tauy_top, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%tauy_top')

  if (allocated(ptend%q))  deallocate(ptend%q, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%q')

  if (allocated(ptend%cflx_srf))   deallocate(ptend%cflx_srf, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%cflx_srf')

  if(allocated(ptend%cflx_top))   deallocate(ptend%cflx_top, stat=ierr)
  if ( ierr /= 0 ) call endrun('physics_ptend_dealloc error: deallocation error for ptend%cflx_top')

end subroutine physics_ptend_dealloc

end module physics_types
