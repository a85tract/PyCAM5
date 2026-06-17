module ndrop

!---------------------------------------------------------------------------------
! Purpose:
!   CAM Interface for droplet activation by modal aerosols
!
! ***N.B.*** This module is currently hardcoded to recognize only the modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use physconst,        only: pi, rhoh2o, mwh2o, r_universal, rh2o, &
                            gravit, latvap, cpair, rair
use constituents,     only: pcnst, cnst_get_ind
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

use wv_saturation,    only: qsat
use phys_control,     only: phys_getopts
use ref_pres,         only: top_lev => trop_cloud_top_lev
use shr_spfn_mod,     only: erf => shr_spfn_erf
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr, &
                            rad_cnst_get_aer_props, rad_cnst_get_mode_props,                &
                            rad_cnst_get_mam_mmr_idx, rad_cnst_get_mode_num_idx
use cam_history,      only: addfld, add_default, phys_decomp, fieldname_len, outfld
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog

implicit none
private
save

public ndrop_init, dropmixnuc

real(r8), allocatable, target :: alogsig(:)     ! natl log of geometric standard dev of aerosol
real(r8), allocatable, target :: exp45logsig(:)
real(r8), allocatable, target :: f1(:)          ! abdul-razzak functions of width
real(r8), allocatable, target :: f2(:)          ! abdul-razzak functions of width

real(r8) :: t0            ! reference temperature
real(r8) :: aten
real(r8) :: surften       ! surface tension of water w/respect to air (N/m)
real(r8) :: alog2, alog3, alogaten
real(r8) :: third, twothird, sixth, zero
real(r8) :: sq2, sqpi

! CCN diagnostic fields
integer,  parameter :: psat=6    ! number of supersaturations to calc ccn concentration
real(r8), parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn concentration
                       (/ 0.02_r8, 0.05_r8, 0.1_r8, 0.2_r8, 0.5_r8, 1.0_r8 /)
character(len=8) :: ccn_name(psat)= &
                    (/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'/)

! indices in state and pbuf structures
integer :: numliq_idx = -1
integer :: kvh_idx    = -1

! description of modal aerosols
integer               :: ntot_amode     ! number of aerosol modes
integer,  allocatable, target :: nspec_amode(:) ! number of chemical species in each aerosol mode
real(r8), allocatable :: sigmag_amode(:)! geometric standard deviation for each aerosol mode
real(r8), allocatable :: dgnumlo_amode(:)
real(r8), allocatable :: dgnumhi_amode(:)
real(r8), allocatable, target :: voltonumblo_amode(:)
real(r8), allocatable, target :: voltonumbhi_amode(:)

logical :: history_aerosol      ! Output the MAM aerosol tendencies
character(len=fieldname_len), allocatable :: fieldname(:)    ! names for drop nuc tendency output fields
character(len=fieldname_len), allocatable :: fieldname_cw(:) ! names for drop nuc tendency output fields

! local indexing for MAM
integer, allocatable, target :: mam_idx(:,:) ! table for local indexing of modal aero number and mmr
integer :: ncnst_tot                  ! total number of mode number conc + mode species

! Indices for MAM species in the ptend%q array.  Needed for prognostic aerosol case.
integer, allocatable, target :: mam_cnst_idx(:,:)


! ptr2d_t is used to create arrays of pointers to 2D fields
type ptr2d_t
   real(r8), pointer :: fld(:,:)
end type ptr2d_t

! modal aerosols
logical :: prog_modal_aero     ! true when modal aerosols are prognostic
logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
                               ! in the ptend object

logical :: use_native_ndrop_init_impl = .false.
logical :: ndrop_init_impl_selected = .false.
logical :: ndrop_init_proof_written = .false.
logical :: use_native_ndrop_init_props_impl = .false.
logical :: ndrop_init_props_impl_selected = .false.
logical :: ndrop_init_props_proof_written = .false.
logical :: use_native_ndrop_dropmixnuc_helpers_impl = .false.
logical :: ndrop_dropmixnuc_helpers_impl_selected = .false.
logical :: ndrop_dropmixnuc_helpers_proof_written = .false.
logical :: ndrop_dropmixnuc_parent_proof_written = .false.
logical :: use_native_ndrop_loadaer_helpers_impl = .false.
logical :: ndrop_loadaer_helpers_impl_selected = .false.
logical :: ndrop_loadaer_helpers_proof_written = .false.
logical :: use_native_ndrop_ccncalc_helpers_impl = .false.
logical :: ndrop_ccncalc_helpers_impl_selected = .false.
logical :: ndrop_ccncalc_helpers_proof_written = .false.
logical :: use_native_ndrop_explmix_impl = .false.
logical :: ndrop_explmix_impl_selected = .false.
logical :: ndrop_explmix_proof_written = .false.
logical :: use_native_ndrop_activate_modal_impl = .false.
logical :: ndrop_activate_modal_impl_selected = .false.
logical :: ndrop_activate_modal_proof_written = .false.
logical :: use_native_ndrop_maxsat_impl = .false.
logical :: ndrop_maxsat_impl_selected = .false.
logical :: ndrop_maxsat_proof_written = .false.

interface
   subroutine ndrop_init_scalars_codon(mwh2o_c, r_universal_c, rhoh2o_c, pi_c, scalars_p) &
        bind(c, name="ndrop_init_scalars_codon")
      use iso_c_binding, only: c_double, c_ptr
      real(c_double), value :: mwh2o_c, r_universal_c, rhoh2o_c, pi_c
      type(c_ptr), value :: scalars_p
   end subroutine ndrop_init_scalars_codon

   subroutine ndrop_init_counts_codon(nmode_c, nspec_amode_p, nspec_max_p, ncnst_tot_p) &
        bind(c, name="ndrop_init_counts_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: nmode_c
      type(c_ptr), value :: nspec_amode_p, nspec_max_p, ncnst_tot_p
   end subroutine ndrop_init_counts_codon

   subroutine ndrop_init_mam_idx_codon(nmode_c, nspec_max_c, nspec_amode_p, mam_idx_p) &
        bind(c, name="ndrop_init_mam_idx_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: nmode_c, nspec_max_c
      type(c_ptr), value :: nspec_amode_p, mam_idx_p
   end subroutine ndrop_init_mam_idx_codon

   subroutine ndrop_mode_props_finalize_codon(nmode_c, pi_c, sigmag_p, dgnumlo_p, dgnumhi_p, &
        alogsig_p, exp45logsig_p, f1_p, f2_p, voltonumblo_p, voltonumbhi_p) &
        bind(c, name="ndrop_mode_props_finalize_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: nmode_c
      real(c_double), value :: pi_c
      type(c_ptr), value :: sigmag_p, dgnumlo_p, dgnumhi_p
      type(c_ptr), value :: alogsig_p, exp45logsig_p, f1_p, f2_p, voltonumblo_p, voltonumbhi_p
   end subroutine ndrop_mode_props_finalize_codon

   subroutine ndrop_dropmixnuc_zero_fields_codon(pcols_c, pver_c, ntot_amode_c, factnum_p, wtke_p) &
        bind(c, name="ndrop_dropmixnuc_zero_fields_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: pcols_c, pver_c, ntot_amode_c
      type(c_ptr), value :: factnum_p, wtke_p
   end subroutine ndrop_dropmixnuc_zero_fields_codon

   subroutine ndrop_dropmixnuc_column_init_codon(i_c, pcols_c, pver_c, top_lev_c, ntot_amode_c, &
        gravit_c, rair_c, zkmin_c, zkmax_c, wmixmin_c, ncldwtr_p, temp_p, pmid_p, pint_p, rpdel_p, &
        zm_p, kvh_p, wsub_p, qcld_p, qncld_p, srcn_p, cs_p, dz_p, nact_p, mact_p, zn_p, ekd_p, &
        csbot_p, csbot_cscen_p, wtke_cen_p, wtke_p, nsource_p, zs_p) &
        bind(c, name="ndrop_dropmixnuc_column_init_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c, ntot_amode_c
      real(c_double), value :: gravit_c, rair_c, zkmin_c, zkmax_c, wmixmin_c
      type(c_ptr), value :: ncldwtr_p, temp_p, pmid_p, pint_p, rpdel_p, zm_p, kvh_p, wsub_p
      type(c_ptr), value :: qcld_p, qncld_p, srcn_p, cs_p, dz_p, nact_p, mact_p, zn_p, ekd_p
      type(c_ptr), value :: csbot_p, csbot_cscen_p, wtke_cen_p, wtke_p, nsource_p, zs_p
   end subroutine ndrop_dropmixnuc_column_init_codon

   subroutine ndrop_dropmixnuc_mix_setup_codon(i_c, pcols_c, pver_c, top_lev_c, ntot_amode_c, &
        dtmicro_c, taumix_internal_pver_inv_c, cldn_p, zs_p, zn_p, csbot_p, ekd_p, nact_p, mact_p, &
        ekk0_p, ekkp_p, ekkm_p, overlapp_p, overlapm_p, count_submix_p, nsubmix_p, dtmix_p) &
        bind(c, name="ndrop_dropmixnuc_mix_setup_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c, ntot_amode_c
      real(c_double), value :: dtmicro_c, taumix_internal_pver_inv_c
      type(c_ptr), value :: cldn_p, zs_p, zn_p, csbot_p, ekd_p, nact_p, mact_p, ekk0_p, ekkp_p
      type(c_ptr), value :: ekkm_p, overlapp_p, overlapm_p, count_submix_p, nsubmix_p, dtmix_p
   end subroutine ndrop_dropmixnuc_mix_setup_codon

   subroutine ndrop_dropmixnuc_aero_column_copy_codon(i_c, pcols_c, pver_c, top_lev_c, ncnst_tot_c, &
        mm_c, slot_c, zero_all_c, raer_fld_p, qqcw_fld_p, raercol_p, raercol_cw_p) &
        bind(c, name="ndrop_dropmixnuc_aero_column_copy_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c, ncnst_tot_c, mm_c, slot_c, zero_all_c
      type(c_ptr), value :: raer_fld_p, qqcw_fld_p, raercol_p, raercol_cw_p
   end subroutine ndrop_dropmixnuc_aero_column_copy_codon

   subroutine ndrop_dropmixnuc_aero_column_copy_all_codon(i_c, pcols_c, pver_c, top_lev_c, &
        ntot_amode_c, ncnst_tot_c, slot_c, raer_ptrs_p, qqcw_ptrs_p, nspec_amode_p, &
        mam_idx_p, raercol_p, raercol_cw_p) bind(c, name="ndrop_dropmixnuc_aero_column_copy_all_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c, ntot_amode_c, ncnst_tot_c, slot_c
      type(c_ptr), value :: raer_ptrs_p, qqcw_ptrs_p, nspec_amode_p, mam_idx_p, raercol_p, raercol_cw_p
   end subroutine ndrop_dropmixnuc_aero_column_copy_all_codon

   subroutine ndrop_dropmixnuc_aero_tend_prepare_codon(i_c, pcols_c, pver_c, top_lev_c, ncnst_tot_c, &
        mm_c, slot_c, dtinv_c, raer_fld_p, qqcw_fld_p, raercol_p, raercol_cw_p, raertend_p, qqcwtend_p) &
        bind(c, name="ndrop_dropmixnuc_aero_tend_prepare_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c, ncnst_tot_c, mm_c, slot_c
      real(c_double), value :: dtinv_c
      type(c_ptr), value :: raer_fld_p, qqcw_fld_p, raercol_p, raercol_cw_p, raertend_p, qqcwtend_p
   end subroutine ndrop_dropmixnuc_aero_tend_prepare_codon

   subroutine ndrop_dropmixnuc_aero_tend_commit_qqcw_codon(i_c, pcols_c, pver_c, top_lev_c, &
        ncnst_tot_c, mm_c, slot_c, qqcw_fld_p, raercol_cw_p) &
        bind(c, name="ndrop_dropmixnuc_aero_tend_commit_qqcw_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c, ncnst_tot_c, mm_c, slot_c
      type(c_ptr), value :: qqcw_fld_p, raercol_cw_p
   end subroutine ndrop_dropmixnuc_aero_tend_commit_qqcw_codon

   subroutine ndrop_dropmixnuc_aero_tend_commit_ptend_codon(i_c, psetcols_c, pver_c, top_lev_c, &
        pcnst_c, lptr_c, raertend_p, ptend_q_p) &
        bind(c, name="ndrop_dropmixnuc_aero_tend_commit_ptend_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: i_c, psetcols_c, pver_c, top_lev_c, pcnst_c, lptr_c
      type(c_ptr), value :: raertend_p, ptend_q_p
   end subroutine ndrop_dropmixnuc_aero_tend_commit_ptend_codon

   subroutine ndrop_dropmixnuc_aero_coltend_codon(i_c, pcols_c, pver_c, mm_c, gravit_c, &
        pdel_p, raertend_p, qqcwtend_p, coltend_out_p, coltend_cw_out_p) &
        bind(c, name="ndrop_dropmixnuc_aero_coltend_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, mm_c
      real(c_double), value :: gravit_c
      type(c_ptr), value :: pdel_p, raertend_p, qqcwtend_p, coltend_out_p, coltend_cw_out_p
   end subroutine ndrop_dropmixnuc_aero_coltend_codon

   subroutine ndrop_dropmixnuc_aero_tend_all_codon(i_c, pcols_c, psetcols_c, pver_c, top_lev_c, &
        ntot_amode_c, ncnst_tot_c, slot_c, dtinv_c, gravit_c, raer_ptrs_p, qqcw_ptrs_p, &
        nspec_amode_p, mam_idx_p, mam_cnst_idx_p, pdel_p, raercol_p, raercol_cw_p, &
        coltend_p, coltend_cw_p, ptend_q_p) bind(c, name="ndrop_dropmixnuc_aero_tend_all_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, psetcols_c, pver_c, top_lev_c
      integer(c_int64_t), value :: ntot_amode_c, ncnst_tot_c, slot_c
      real(c_double), value :: dtinv_c, gravit_c
      type(c_ptr), value :: raer_ptrs_p, qqcw_ptrs_p, nspec_amode_p, mam_idx_p, mam_cnst_idx_p
      type(c_ptr), value :: pdel_p, raercol_p, raercol_cw_p, coltend_p, coltend_cw_p, ptend_q_p
   end subroutine ndrop_dropmixnuc_aero_tend_all_codon

   subroutine ndrop_dropmixnuc_finalize_column_codon(i_c, pcols_c, pver_c, top_lev_c, dtinv_c, &
        gravit_c, qcld_p, ncldwtr_p, pdel_p, nsource_p, ndropmix_p, tendnd_p, ndropcol_p) &
        bind(c, name="ndrop_dropmixnuc_finalize_column_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c
      real(c_double), value :: dtinv_c, gravit_c
      type(c_ptr), value :: qcld_p, ncldwtr_p, pdel_p, nsource_p, ndropmix_p, tendnd_p, ndropcol_p
   end subroutine ndrop_dropmixnuc_finalize_column_codon

   subroutine ndrop_dropmixnuc_clear_old_cloud_codon(i_c, k_c, pcols_c, pver_c, ntot_amode_c, &
        ncnst_tot_c, nsav_c, dtinv_c, qcld_p, nsource_p, nspec_amode_p, mam_idx_p, raercol_p, &
        raercol_cw_p) bind(c, name="ndrop_dropmixnuc_clear_old_cloud_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, ntot_amode_c, ncnst_tot_c, nsav_c
      real(c_double), value :: dtinv_c
      type(c_ptr), value :: qcld_p, nsource_p, nspec_amode_p, mam_idx_p, raercol_p, raercol_cw_p
   end subroutine ndrop_dropmixnuc_clear_old_cloud_codon

   subroutine ndrop_dropmixnuc_factnum_store_codon(i_c, k_c, pcols_c, pver_c, ntot_amode_c, &
        fn_p, factnum_p) bind(c, name="ndrop_dropmixnuc_factnum_store_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, ntot_amode_c
      type(c_ptr), value :: fn_p, factnum_p
   end subroutine ndrop_dropmixnuc_factnum_store_codon

   subroutine ndrop_dropmixnuc_shrink_cloud_codon(i_c, k_c, pcols_c, pver_c, ntot_amode_c, &
        ncnst_tot_c, nsav_c, dtinv_c, cldn_tmp_c, cldo_tmp_c, qcld_p, nsource_p, &
        nspec_amode_p, mam_idx_p, raercol_p, raercol_cw_p) &
        bind(c, name="ndrop_dropmixnuc_shrink_cloud_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, ntot_amode_c, ncnst_tot_c, nsav_c
      real(c_double), value :: dtinv_c, cldn_tmp_c, cldo_tmp_c
      type(c_ptr), value :: qcld_p, nsource_p, nspec_amode_p, mam_idx_p, raercol_p, raercol_cw_p
   end subroutine ndrop_dropmixnuc_shrink_cloud_codon

   subroutine ndrop_dropmixnuc_grow_cloud_number_update_codon(i_c, k_c, pcols_c, pver_c, &
        ncnst_tot_c, nsav_c, mm_c, dtinv_c, dumc_c, fn_m_c, raer_fld_p, qcld_p, &
        nsource_p, raercol_p, raercol_cw_p) &
        bind(c, name="ndrop_dropmixnuc_grow_cloud_number_update_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, ncnst_tot_c, nsav_c, mm_c
      real(c_double), value :: dtinv_c, dumc_c, fn_m_c
      type(c_ptr), value :: raer_fld_p, qcld_p, nsource_p, raercol_p, raercol_cw_p
   end subroutine ndrop_dropmixnuc_grow_cloud_number_update_codon

   subroutine ndrop_dropmixnuc_grow_cloud_update_all_codon(i_c, k_c, pcols_c, pver_c, &
        ntot_amode_c, ncnst_tot_c, nsav_c, dtinv_c, dumc_c, raer_ptrs_p, nspec_amode_p, &
        mam_idx_p, fn_p, fm_p, qcld_p, nsource_p, raercol_p, raercol_cw_p, factnum_p) &
        bind(c, name="ndrop_dropmixnuc_grow_cloud_update_all_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, ntot_amode_c, ncnst_tot_c, nsav_c
      real(c_double), value :: dtinv_c, dumc_c
      type(c_ptr), value :: raer_ptrs_p, nspec_amode_p, mam_idx_p, fn_p, fm_p
      type(c_ptr), value :: qcld_p, nsource_p, raercol_p, raercol_cw_p, factnum_p
   end subroutine ndrop_dropmixnuc_grow_cloud_update_all_codon

   subroutine ndrop_dropmixnuc_grow_cloud_species_update_codon(i_c, k_c, pcols_c, pver_c, &
        ncnst_tot_c, nsav_c, mm_c, dum_c, raer_fld_p, raercol_p, raercol_cw_p) &
        bind(c, name="ndrop_dropmixnuc_grow_cloud_species_update_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, ncnst_tot_c, nsav_c, mm_c
      real(c_double), value :: dum_c
      type(c_ptr), value :: raer_fld_p, raercol_p, raercol_cw_p
   end subroutine ndrop_dropmixnuc_grow_cloud_species_update_codon

   subroutine ndrop_dropmixnuc_old_cloud_activate_update_codon(i_c, k_c, kp1_c, pcols_c, &
        pver_c, ntot_amode_c, ncnst_tot_c, nsav_c, dumc_c, dum_c, cs_ik_c, dz_ik_c, &
        taumix_internal_pver_inv_c, fluxn_p, fluxm_p, nact_p, mact_p, mam_idx_p, &
        raercol_p, raercol_cw_p, srcn_p, nsource_p) &
        bind(c, name="ndrop_dropmixnuc_old_cloud_activate_update_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, k_c, kp1_c, pcols_c, pver_c, ntot_amode_c
      integer(c_int64_t), value :: ncnst_tot_c, nsav_c
      real(c_double), value :: dumc_c, dum_c, cs_ik_c, dz_ik_c, taumix_internal_pver_inv_c
      type(c_ptr), value :: fluxn_p, fluxm_p, nact_p, mact_p, mam_idx_p
      type(c_ptr), value :: raercol_p, raercol_cw_p, srcn_p, nsource_p
   end subroutine ndrop_dropmixnuc_old_cloud_activate_update_codon

   function ndrop_dropmixnuc_activation_loops_codon(i_c, pcols_c, pver_c, top_lev_c, &
        ntot_amode_c, ncnst_tot_c, nsav_c, dtmicro_c, dtinv_c, rair_c, p0_c, t0_c, &
        rhoh2o_c, latvap_c, cpair_c, rh2o_c, gravit_c, pi_c, aten_c, twothird_c, &
        sq2_c, sqpi_c, sixth_c, zero_c, cldn_p, cldo_p, cldn_regen_p, temp_p, cs_p, &
        dz_p, wtke_p, wtke_cen_p, zs_p, ekd_p, csbot_cscen_p, qs_act_p, qcld_p, &
        srcn_p, nsource_p, factnum_p, nact_p, mact_p, taumix_internal_pver_inv_p, raer_ptrs_p, &
        qqcw_ptrs_p, nspec_amode_p, mam_idx_p, species_specdens_p, species_spechygro_p, &
        voltonumblo_p, voltonumbhi_p, alogsig_p, exp45logsig_p, f1_p, f2_p, raercol_p, &
        raercol_cw_p, naermod_p, vaerosol_p, hygro_p, fn_p, fm_p, fluxn_p, fluxm_p, &
        flux_fullact_p, zeta_p, eta_p, etafactor2_p, sqrtg_p, amcube_p, smc_p, lnsm_p, &
        sumflxn_p, sumflxm_p, sumfn_p, sumfm_p, fnold_p, fmold_p) result(status_c) &
        bind(c, name="ndrop_dropmixnuc_activation_loops_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c, ntot_amode_c
      integer(c_int64_t), value :: ncnst_tot_c, nsav_c
      real(c_double), value :: dtmicro_c, dtinv_c, rair_c, p0_c, t0_c, rhoh2o_c
      real(c_double), value :: latvap_c, cpair_c, rh2o_c, gravit_c, pi_c, aten_c
      real(c_double), value :: twothird_c, sq2_c, sqpi_c, sixth_c, zero_c
      type(c_ptr), value :: cldn_p, cldo_p, cldn_regen_p, temp_p, cs_p, dz_p
      type(c_ptr), value :: wtke_p, wtke_cen_p, zs_p, ekd_p, csbot_cscen_p, qs_act_p
      type(c_ptr), value :: qcld_p, srcn_p, nsource_p, factnum_p, nact_p, mact_p
      type(c_ptr), value :: taumix_internal_pver_inv_p, raer_ptrs_p, qqcw_ptrs_p
      type(c_ptr), value :: nspec_amode_p, mam_idx_p, species_specdens_p, species_spechygro_p
      type(c_ptr), value :: voltonumblo_p, voltonumbhi_p, alogsig_p, exp45logsig_p
      type(c_ptr), value :: f1_p, f2_p, raercol_p, raercol_cw_p, naermod_p
      type(c_ptr), value :: vaerosol_p, hygro_p, fn_p, fm_p, fluxn_p, fluxm_p
      type(c_ptr), value :: flux_fullact_p, zeta_p, eta_p, etafactor2_p, sqrtg_p
      type(c_ptr), value :: amcube_p, smc_p, lnsm_p, sumflxn_p, sumflxm_p
      type(c_ptr), value :: sumfn_p, sumfm_p, fnold_p, fmold_p
      integer(c_int64_t) :: status_c
   end function ndrop_dropmixnuc_activation_loops_codon

   subroutine ndrop_dropmixnuc_srcn_from_nact_codon(pver_c, top_lev_c, ntot_amode_c, ncnst_tot_c, &
        nsav_c, taumix_internal_pver_inv_c, nact_p, mam_idx_p, raercol_p, raercol_cw_p, srcn_p) &
        bind(c, name="ndrop_dropmixnuc_srcn_from_nact_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: pver_c, top_lev_c, ntot_amode_c, ncnst_tot_c, nsav_c
      real(c_double), value :: taumix_internal_pver_inv_c
      type(c_ptr), value :: nact_p, mam_idx_p, raercol_p, raercol_cw_p, srcn_p
   end subroutine ndrop_dropmixnuc_srcn_from_nact_codon

   subroutine ndrop_dropmixnuc_source_from_act_codon(pver_c, top_lev_c, ncnst_tot_c, m_c, mm_c, &
        nsav_c, taumix_internal_pver_inv_c, act_p, raercol_p, raercol_cw_p, source_p) &
        bind(c, name="ndrop_dropmixnuc_source_from_act_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: pver_c, top_lev_c, ncnst_tot_c, m_c, mm_c, nsav_c
      real(c_double), value :: taumix_internal_pver_inv_c
      type(c_ptr), value :: act_p, raercol_p, raercol_cw_p, source_p
   end subroutine ndrop_dropmixnuc_source_from_act_codon

   subroutine ndrop_dropmixnuc_evaporate_clear_layers_codon(i_c, pcols_c, pver_c, top_lev_c, &
        ntot_amode_c, ncnst_tot_c, nnew_c, cldn_p, qcld_p, nspec_amode_p, mam_idx_p, &
        raercol_p, raercol_cw_p) bind(c, name="ndrop_dropmixnuc_evaporate_clear_layers_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c, ntot_amode_c, ncnst_tot_c, nnew_c
      type(c_ptr), value :: cldn_p, qcld_p, nspec_amode_p, mam_idx_p, raercol_p, raercol_cw_p
   end subroutine ndrop_dropmixnuc_evaporate_clear_layers_codon

   subroutine ndrop_dropmixnuc_swap_slots_codon(nsav_p, nnew_p) &
        bind(c, name="ndrop_dropmixnuc_swap_slots_codon")
      use iso_c_binding, only: c_ptr
      type(c_ptr), value :: nsav_p, nnew_p
   end subroutine ndrop_dropmixnuc_swap_slots_codon

   subroutine ndrop_dropmixnuc_submix_iter_init_codon(pver_c, qcld_p, qncld_p, srcn_p, nsav_p, nnew_p) &
        bind(c, name="ndrop_dropmixnuc_submix_iter_init_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: pver_c
      type(c_ptr), value :: qcld_p, qncld_p, srcn_p, nsav_p, nnew_p
   end subroutine ndrop_dropmixnuc_submix_iter_init_codon

   subroutine ndrop_dropmixnuc_submix_all_codon(pver_c, top_lev_c, ntot_amode_c, ncnst_tot_c, &
        dtmix_c, taumix_internal_pver_inv_c, nact_p, mact_p, mam_idx_p, nspec_amode_p, &
        ekkp_p, ekkm_p, overlapp_p, overlapm_p, qcld_p, qncld_p, srcn_p, source_p, &
        raercol_p, raercol_cw_p, nsav_p, nnew_p) bind(c, name="ndrop_dropmixnuc_submix_all_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: pver_c, top_lev_c, ntot_amode_c, ncnst_tot_c
      real(c_double), value :: dtmix_c, taumix_internal_pver_inv_c
      type(c_ptr), value :: nact_p, mact_p, mam_idx_p, nspec_amode_p
      type(c_ptr), value :: ekkp_p, ekkm_p, overlapp_p, overlapm_p
      type(c_ptr), value :: qcld_p, qncld_p, srcn_p, source_p
      type(c_ptr), value :: raercol_p, raercol_cw_p, nsav_p, nnew_p
   end subroutine ndrop_dropmixnuc_submix_all_codon

   subroutine ndrop_dropmixnuc_zero_tendencies_codon(pver_c, raertend_p, qqcwtend_p) &
        bind(c, name="ndrop_dropmixnuc_zero_tendencies_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: pver_c
      type(c_ptr), value :: raertend_p, qqcwtend_p
   end subroutine ndrop_dropmixnuc_zero_tendencies_codon

   subroutine ndrop_loadaer_zero_codon(istart_c, istop_c, vaerosol_p, hygro_p) &
        bind(c, name="ndrop_loadaer_zero_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: istart_c, istop_c
      type(c_ptr), value :: vaerosol_p, hygro_p
   end subroutine ndrop_loadaer_zero_codon

   subroutine ndrop_loadaer_species_accum_codon(istart_c, istop_c, k_c, pcols_c, phase_c, &
        specdens_c, spechygro_c, raer_p, qqcw_p, vaerosol_p, hygro_p) &
        bind(c, name="ndrop_loadaer_species_accum_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: istart_c, istop_c, k_c, pcols_c, phase_c
      real(c_double), value :: specdens_c, spechygro_c
      type(c_ptr), value :: raer_p, qqcw_p, vaerosol_p, hygro_p
   end subroutine ndrop_loadaer_species_accum_codon

   subroutine ndrop_loadaer_species_batch_codon(istart_c, istop_c, k_c, pcols_c, nspec_c, phase_c, &
        raer_ptrs_p, qqcw_ptrs_p, specdens_p, spechygro_p, vaerosol_p, hygro_p) &
        bind(c, name="ndrop_loadaer_species_batch_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: istart_c, istop_c, k_c, pcols_c, nspec_c, phase_c
      type(c_ptr), value :: raer_ptrs_p, qqcw_ptrs_p, specdens_p, spechygro_p, vaerosol_p, hygro_p
   end subroutine ndrop_loadaer_species_batch_codon

   subroutine ndrop_loadaer_finalize_volume_codon(istart_c, istop_c, k_c, pcols_c, cs_p, &
        vaerosol_p, hygro_p) bind(c, name="ndrop_loadaer_finalize_volume_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: istart_c, istop_c, k_c, pcols_c
      type(c_ptr), value :: cs_p, vaerosol_p, hygro_p
   end subroutine ndrop_loadaer_finalize_volume_codon

   subroutine ndrop_loadaer_number_codon(istart_c, istop_c, k_c, pcols_c, phase_c, &
        voltonumblo_c, voltonumbhi_c, raer_p, qqcw_p, cs_p, vaerosol_p, naerosol_p) &
        bind(c, name="ndrop_loadaer_number_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: istart_c, istop_c, k_c, pcols_c, phase_c
      real(c_double), value :: voltonumblo_c, voltonumbhi_c
      type(c_ptr), value :: raer_p, qqcw_p, cs_p, vaerosol_p, naerosol_p
   end subroutine ndrop_loadaer_number_codon

   subroutine ndrop_loadaer_direct_codon(istart_c, istop_c, k_c, pcols_c, nspec_c, phase_c, &
        voltonumblo_c, voltonumbhi_c, species_raer_ptrs_p, species_qqcw_ptrs_p, specdens_p, &
        spechygro_p, num_raer_p, num_qqcw_p, cs_p, vaerosol_p, hygro_p, naerosol_p) &
        bind(c, name="ndrop_loadaer_direct_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: istart_c, istop_c, k_c, pcols_c, nspec_c, phase_c
      real(c_double), value :: voltonumblo_c, voltonumbhi_c
      type(c_ptr), value :: species_raer_ptrs_p, species_qqcw_ptrs_p, specdens_p, spechygro_p
      type(c_ptr), value :: num_raer_p, num_qqcw_p, cs_p, vaerosol_p, hygro_p, naerosol_p
   end subroutine ndrop_loadaer_direct_codon

   subroutine ndrop_ccncalc_zero_codon(pcols_c, pver_c, psat_c, ccn_p) &
        bind(c, name="ndrop_ccncalc_zero_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: pcols_c, pver_c, psat_c
      type(c_ptr), value :: ccn_p
   end subroutine ndrop_ccncalc_zero_codon

   subroutine ndrop_ccncalc_level_coeffs_codon(ncol_c, k_c, pcols_c, surften_coef_c, &
        smcoefcoef_c, tair_p, smcoef_p) bind(c, name="ndrop_ccncalc_level_coeffs_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, k_c, pcols_c
      real(c_double), value :: surften_coef_c, smcoefcoef_c
      type(c_ptr), value :: tair_p, smcoef_p
   end subroutine ndrop_ccncalc_level_coeffs_codon

   subroutine ndrop_ccncalc_mode_accum_codon(ncol_c, k_c, pcols_c, pver_c, psat_c, &
        amcubecoef_m_c, argfactor_m_c, naerosol_p, vaerosol_p, hygro_p, smcoef_p, &
        super_p, amcube_p, sm_p, arg_p, ccn_p) bind(c, name="ndrop_ccncalc_mode_accum_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, k_c, pcols_c, pver_c, psat_c
      real(c_double), value :: amcubecoef_m_c, argfactor_m_c
      type(c_ptr), value :: naerosol_p, vaerosol_p, hygro_p, smcoef_p, super_p
      type(c_ptr), value :: amcube_p, sm_p, arg_p, ccn_p
   end subroutine ndrop_ccncalc_mode_accum_codon

   subroutine ndrop_ccncalc_scale_codon(ncol_c, pcols_c, pver_c, psat_c, ccn_p) &
        bind(c, name="ndrop_ccncalc_scale_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, psat_c
      type(c_ptr), value :: ccn_p
   end subroutine ndrop_ccncalc_scale_codon

   subroutine ndrop_ccncalc_direct_codon(ncol_c, pcols_c, pver_c, top_lev_c, psat_c, &
        ntot_amode_c, ncnst_tot_c, pi_c, surften_coef_c, smcoefcoef_c, nspec_amode_p, &
        species_raer_ptrs_p, species_qqcw_ptrs_p, specdens_p, spechygro_p, num_raer_ptrs_p, &
        num_qqcw_ptrs_p, voltonumblo_p, voltonumbhi_p, alogsig_p, exp45logsig_p, tair_p, &
        cs_p, super_p, naerosol_p, vaerosol_p, hygro_p, amcube_p, smcoef_p, sm_p, arg_p, ccn_p) &
        bind(c, name="ndrop_ccncalc_direct_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, psat_c
      integer(c_int64_t), value :: ntot_amode_c, ncnst_tot_c
      real(c_double), value :: pi_c, surften_coef_c, smcoefcoef_c
      type(c_ptr), value :: nspec_amode_p, species_raer_ptrs_p, species_qqcw_ptrs_p
      type(c_ptr), value :: specdens_p, spechygro_p, num_raer_ptrs_p, num_qqcw_ptrs_p
      type(c_ptr), value :: voltonumblo_p, voltonumbhi_p, alogsig_p, exp45logsig_p
      type(c_ptr), value :: tair_p, cs_p, super_p, naerosol_p, vaerosol_p, hygro_p
      type(c_ptr), value :: amcube_p, smcoef_p, sm_p, arg_p, ccn_p
   end subroutine ndrop_ccncalc_direct_codon

   function ndrop_activate_modal_core_codon(wbar_c, sigw_c, wdiab_c, wminf_c, wmaxf_c, &
        tair_c, rhoair_c, qs_c, nmode_c, rair_c, p0_c, t0_c, rhoh2o_c, latvap_c, &
        cpair_c, rh2o_c, gravit_c, pi_c, aten_c, twothird_c, sq2_c, sqpi_c, &
        sixth_c, zero_c, na_p, volume_p, hygro_p, alogsig_p, exp45logsig_p, &
        f1_p, f2_p, fn_p, fm_p, fluxn_p, fluxm_p, flux_fullact_p, zeta_p, eta_p, &
        etafactor2_p, sqrtg_p, amcube_p, smc_p, lnsm_p, sumflxn_p, sumflxm_p, &
        sumfn_p, sumfm_p, fnold_p, fmold_p) result(status_c) &
        bind(c, name="ndrop_activate_modal_core_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      real(c_double), value :: wbar_c, sigw_c, wdiab_c, wminf_c, wmaxf_c
      real(c_double), value :: tair_c, rhoair_c, qs_c, rair_c, p0_c, t0_c
      real(c_double), value :: rhoh2o_c, latvap_c, cpair_c, rh2o_c, gravit_c
      real(c_double), value :: pi_c, aten_c, twothird_c, sq2_c, sqpi_c, sixth_c, zero_c
      integer(c_int64_t), value :: nmode_c
      type(c_ptr), value :: na_p, volume_p, hygro_p, alogsig_p, exp45logsig_p
      type(c_ptr), value :: f1_p, f2_p, fn_p, fm_p, fluxn_p, fluxm_p, flux_fullact_p
      type(c_ptr), value :: zeta_p, eta_p, etafactor2_p, sqrtg_p, amcube_p, smc_p, lnsm_p
      type(c_ptr), value :: sumflxn_p, sumflxm_p, sumfn_p, sumfm_p, fnold_p, fmold_p
      integer(c_int64_t) :: status_c
   end function ndrop_activate_modal_core_codon

   function ndrop_maxsat_codon(nmode_c, zeta_p, eta_p, smc_p, f1_p, f2_p) result(smax_c) &
        bind(c, name="ndrop_maxsat_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: nmode_c
      type(c_ptr), value :: zeta_p, eta_p, smc_p, f1_p, f2_p
      real(c_double) :: smax_c
   end function ndrop_maxsat_codon

   subroutine ndrop_explmix_codon(pver_c, top_lev_c, surfrate_c, flxconv_c, dt_c, is_unact_c, &
        q_p, src_p, ekkp_p, ekkm_p, overlapp_p, overlapm_p, qold_p, qactold_p) &
        bind(c, name="ndrop_explmix_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: pver_c, top_lev_c, is_unact_c
      real(c_double), value :: surfrate_c, flxconv_c, dt_c
      type(c_ptr), value :: q_p, src_p, ekkp_p, ekkm_p, overlapp_p, overlapm_p, qold_p, qactold_p
   end subroutine ndrop_explmix_codon
end interface

!===============================================================================
contains
!===============================================================================

subroutine ndrop_init_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (ndrop_init_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('NDROP_INIT_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_ndrop_init_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_ndrop_init_impl = .false.
   end if

   ndrop_init_impl_selected = .true.

   if (masterproc) then
      if (use_native_ndrop_init_impl) then
         write(iulog,*) 'ndrop_init implementation = native'
      else
         write(iulog,*) 'ndrop_init implementation = codon'
      end if
   end if

end subroutine ndrop_init_select_impl

!===============================================================================

subroutine ndrop_init_proof_once()

   if (ndrop_init_proof_written) return
   ndrop_init_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'ndrop_init direct = codon scalar constants, local mode counts, and mam_idx plan; native rad_constituents, phys_getopts, history callbacks'
   end if

end subroutine ndrop_init_proof_once

!===============================================================================

subroutine ndrop_init_props_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (ndrop_init_props_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('NDROP_INIT_PROPS_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_ndrop_init_props_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_ndrop_init_props_impl = .false.
   end if

   ndrop_init_props_impl_selected = .true.

   if (masterproc) then
      if (use_native_ndrop_init_props_impl) then
         write(iulog,*) 'ndrop_init_props implementation = native'
      else
         write(iulog,*) 'ndrop_init_props implementation = codon'
      end if
   end if

end subroutine ndrop_init_props_select_impl

!===============================================================================

subroutine ndrop_init_props_proof_once()

   if (ndrop_init_props_proof_written) return
   ndrop_init_props_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'ndrop_init_props entered (modal aerosol width conversion helper = codon)'
   end if

end subroutine ndrop_init_props_proof_once

!===============================================================================

subroutine ndrop_dropmixnuc_helpers_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (ndrop_dropmixnuc_helpers_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('NDROP_DROPMIXNUC_HELPERS_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_ndrop_dropmixnuc_helpers_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_ndrop_dropmixnuc_helpers_impl = .false.
   end if

   ndrop_dropmixnuc_helpers_impl_selected = .true.

   if (masterproc) then
      if (use_native_ndrop_dropmixnuc_helpers_impl) then
         write(iulog,*) 'ndrop_dropmixnuc_helpers implementation = native'
      else
         write(iulog,*) 'ndrop_dropmixnuc_helpers implementation = codon'
      end if
   end if

end subroutine ndrop_dropmixnuc_helpers_select_impl

!===============================================================================

subroutine ndrop_dropmixnuc_helpers_proof_once()

   if (ndrop_dropmixnuc_helpers_proof_written) return
   ndrop_dropmixnuc_helpers_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'ndrop_dropmixnuc_helpers entered (array setup/activation/grow-shrink/' // &
           'oldcloud/mix/source/submix all/aero pointer-table batches/grow batch/aero tend all/clear/finalize direct = codon)'
   end if

end subroutine ndrop_dropmixnuc_helpers_proof_once

!===============================================================================

subroutine ndrop_loadaer_helpers_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (ndrop_loadaer_helpers_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('NDROP_LOADAER_HELPERS_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_ndrop_loadaer_helpers_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_ndrop_loadaer_helpers_impl = .false.
   end if

   ndrop_loadaer_helpers_impl_selected = .true.

   if (masterproc) then
      if (use_native_ndrop_loadaer_helpers_impl) then
         write(iulog,*) 'ndrop_loadaer_helpers implementation = native'
      else
         write(iulog,*) 'ndrop_loadaer_helpers implementation = codon'
      end if
   end if

end subroutine ndrop_loadaer_helpers_select_impl

!===============================================================================

subroutine ndrop_loadaer_helpers_proof_once()

   if (ndrop_loadaer_helpers_proof_written) return
   ndrop_loadaer_helpers_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'ndrop_loadaer direct = codon; rad_constituents/pbuf pointer lookup native CAM API island'
   end if

end subroutine ndrop_loadaer_helpers_proof_once

!===============================================================================

subroutine ndrop_ccncalc_helpers_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (ndrop_ccncalc_helpers_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('NDROP_CCNCALC_HELPERS_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_ndrop_ccncalc_helpers_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_ndrop_ccncalc_helpers_impl = .false.
   end if

   ndrop_ccncalc_helpers_impl_selected = .true.

   if (masterproc) then
      if (use_native_ndrop_ccncalc_helpers_impl) then
         write(iulog,*) 'ndrop_ccncalc_helpers implementation = native'
      else
         write(iulog,*) 'ndrop_ccncalc_helpers implementation = codon'
      end if
   end if

end subroutine ndrop_ccncalc_helpers_select_impl

!===============================================================================

subroutine ndrop_ccncalc_helpers_proof_once()

   if (ndrop_ccncalc_helpers_proof_written) return
   ndrop_ccncalc_helpers_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'ndrop_ccncalc direct = codon; allocation/state/rad_constituents/pbuf native CAM API islands'
   end if

end subroutine ndrop_ccncalc_helpers_proof_once

!===============================================================================

subroutine ndrop_explmix_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (ndrop_explmix_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('NDROP_EXPLMIX_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_ndrop_explmix_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_ndrop_explmix_impl = .false.
   end if

   ndrop_explmix_impl_selected = .true.

   if (masterproc) then
      if (use_native_ndrop_explmix_impl) then
         write(iulog,*) 'ndrop_explmix implementation = native'
      else
         write(iulog,*) 'ndrop_explmix implementation = codon'
      end if
   end if

end subroutine ndrop_explmix_select_impl

!===============================================================================

subroutine ndrop_explmix_proof_once()

   if (ndrop_explmix_proof_written) return
   ndrop_explmix_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'ndrop_explmix entered (vertical explicit mixing direct = codon)'
   end if

end subroutine ndrop_explmix_proof_once

!===============================================================================

subroutine ndrop_activate_modal_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (ndrop_activate_modal_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('NDROP_ACTIVATE_MODAL_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_ndrop_activate_modal_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_ndrop_activate_modal_impl = .false.
   end if

   ndrop_activate_modal_impl_selected = .true.

   if (masterproc) then
      if (use_native_ndrop_activate_modal_impl) then
         write(iulog,*) 'ndrop_activate_modal implementation = native'
      else
         write(iulog,*) 'ndrop_activate_modal implementation = codon'
      end if
   end if

end subroutine ndrop_activate_modal_select_impl

!===============================================================================

subroutine ndrop_activate_modal_proof_once()

   if (ndrop_activate_modal_proof_written) return
   ndrop_activate_modal_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'ndrop_activate_modal direct = codon; qsat native thermo island'
   end if

end subroutine ndrop_activate_modal_proof_once

!===============================================================================

subroutine ndrop_maxsat_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (ndrop_maxsat_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('NDROP_MAXSAT_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_ndrop_maxsat_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_ndrop_maxsat_impl = .false.
   end if

   ndrop_maxsat_impl_selected = .true.

   if (masterproc) then
      if (use_native_ndrop_maxsat_impl) then
         write(iulog,*) 'ndrop_maxsat implementation = native'
      else
         write(iulog,*) 'ndrop_maxsat implementation = codon'
      end if
   end if

end subroutine ndrop_maxsat_select_impl

!===============================================================================

subroutine ndrop_maxsat_proof_once()

   if (ndrop_maxsat_proof_written) return
   ndrop_maxsat_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'ndrop_maxsat entered (activation supersaturation solve = codon)'
   end if

end subroutine ndrop_maxsat_proof_once

!===============================================================================

subroutine ndrop_init

   use iso_c_binding, only: c_int64_t, c_loc

   integer  :: ii, l, lptr, m, mm
   integer  :: nspec_max            ! max number of species in a mode
   integer(c_int64_t), target :: nspec_max_c, ncnst_tot_c
   real(r8), target :: init_scalars(12)
   character(len=32)   :: tmpname
   character(len=32)   :: tmpname_cw
   character(len=128)  :: long_name
   character(len=8)    :: unit
   logical :: history_amwg         ! output the variables used by the AMWG diag package

   !-------------------------------------------------------------------------------

   call ndrop_init_select_impl()

   ! get indices into state%q and pbuf structures
   call cnst_get_ind('NUMLIQ', numliq_idx)

   kvh_idx      = pbuf_get_index('kvh')

   if (.not. use_native_ndrop_init_impl) then
      call ndrop_init_scalars_codon(mwh2o, r_universal, rhoh2o, pi, c_loc(init_scalars(1)))
      zero     = init_scalars(1)
      third    = init_scalars(2)
      twothird = init_scalars(3)
      sixth    = init_scalars(4)
      sq2      = init_scalars(5)
      sqpi     = init_scalars(6)
      t0       = init_scalars(7)
      surften  = init_scalars(8)
      aten     = init_scalars(9)
      alogaten = init_scalars(10)
      alog2    = init_scalars(11)
      alog3    = init_scalars(12)
   else
      zero     = 0._r8
      third    = 1._r8/3._r8
      twothird = 2._r8*third
      sixth    = 1._r8/6._r8
      sq2      = sqrt(2._r8)
      sqpi     = sqrt(pi)

      t0       = 273._r8
      surften  = 0.076_r8
      aten     = 2._r8*mwh2o*surften/(r_universal*t0*rhoh2o)
      alogaten = log(aten)
      alog2    = log(2._r8)
      alog3    = log(3._r8)
   end if

   ! get info about the modal aerosols
   ! get ntot_amode
   call rad_cnst_get_info(0, nmodes=ntot_amode)

   allocate( &
      nspec_amode(ntot_amode),  &
      sigmag_amode(ntot_amode), &
      dgnumlo_amode(ntot_amode), &
      dgnumhi_amode(ntot_amode), &
      alogsig(ntot_amode),      &
      exp45logsig(ntot_amode),  &
      f1(ntot_amode),           &
      f2(ntot_amode),           &
      voltonumblo_amode(ntot_amode), &
      voltonumbhi_amode(ntot_amode)  )

   do m = 1, ntot_amode
      ! use only if width of size distribution is prescribed

      ! get mode info
      call rad_cnst_get_info(0, m, nspec=nspec_amode(m))

      ! get mode properties
      call rad_cnst_get_mode_props(0, m, sigmag=sigmag_amode(m),  &
         dgnumhi=dgnumhi_amode(m), dgnumlo=dgnumlo_amode(m))
   end do

   call ndrop_mode_props_finalize(ntot_amode, sigmag_amode, dgnumlo_amode, dgnumhi_amode, &
        alogsig, exp45logsig, f1, f2, voltonumblo_amode, voltonumbhi_amode)
      
   ! Init the table for local indexing of mam number conc and mmr.
   ! This table uses species index 0 for the number conc.

   ! Find max number of species in all the modes, and the total
   ! number of mode number concentrations + mode species
   if (.not. use_native_ndrop_init_impl) then
      call ndrop_init_counts_codon(int(ntot_amode, c_int64_t), c_loc(nspec_amode(1)), &
           c_loc(nspec_max_c), c_loc(ncnst_tot_c))
      nspec_max = int(nspec_max_c)
      ncnst_tot = int(ncnst_tot_c)
   else
      nspec_max = nspec_amode(1)
      ncnst_tot = nspec_amode(1) + 1
      do m = 2, ntot_amode
         nspec_max = max(nspec_max, nspec_amode(m))
         ncnst_tot = ncnst_tot + nspec_amode(m) + 1
      end do
   end if

   allocate( &
      mam_idx(ntot_amode,0:nspec_max),      &
      mam_cnst_idx(ntot_amode,0:nspec_max), &
      fieldname(ncnst_tot),                 &
      fieldname_cw(ncnst_tot)               )

   ! Local indexing compresses the mode and number/mass indicies into one index.
   ! This indexing is used by the pointer arrays used to reference state and pbuf
   ! fields.
   if (.not. use_native_ndrop_init_impl) then
      call ndrop_init_mam_idx_codon(int(ntot_amode, c_int64_t), int(nspec_max, c_int64_t), &
           c_loc(nspec_amode(1)), c_loc(mam_idx(1,0)))
      call ndrop_init_proof_once()
   else
      ii = 0
      do m = 1, ntot_amode
         do l = 0, nspec_amode(m)
            ii = ii + 1
            mam_idx(m,l) = ii
         end do
      end do
   end if

   ! Add dropmixnuc tendencies for all modal aerosol species

   call phys_getopts(history_amwg_out = history_amwg, &
                     history_aerosol_out = history_aerosol, &
                     prog_modal_aero_out=prog_modal_aero)


   do m = 1, ntot_amode
      do l = 0, nspec_amode(m)   ! loop over number + chem constituents

         mm = mam_idx(m,l)

         unit = 'kg/m2/s'
         if (l == 0) then   ! number
            unit = '#/m2/s'
         end if

         if (l == 0) then   ! number
            call rad_cnst_get_info(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
         else
            call rad_cnst_get_info(0, m, l, spec_name=tmpname, spec_name_cw=tmpname_cw)
         end if

         fieldname(mm)    = trim(tmpname) // '_mixnuc1'
         fieldname_cw(mm) = trim(tmpname_cw) // '_mixnuc1'

         if (prog_modal_aero) then

            ! To set tendencies in the ptend object need to get the constituent indices
            ! for the prognostic species
            if (l == 0) then   ! number
               call rad_cnst_get_mode_num_idx(m, lptr)
            else
               call rad_cnst_get_mam_mmr_idx(m, l, lptr)
            end if
            mam_cnst_idx(m,l) = lptr
            lq(lptr)          = .true.

            ! Add tendency fields to the history only when prognostic MAM is enabled.
            long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency'
            call addfld(fieldname(mm), unit, 1, 'A', long_name, phys_decomp)

            long_name = trim(tmpname_cw) // ' dropmixnuc mixnuc column tendency'
            call addfld(fieldname_cw(mm), unit, 1, 'A', long_name, phys_decomp)

            if (history_aerosol) then
               call add_default(fieldname(mm), 1, ' ')
               call add_default(fieldname_cw(mm), 1, ' ')
            end if



         end if
            
      end do
   end do

   call addfld('CCN1    ','#/cm3   ',pver, 'A','CCN concentration at S=0.02%',phys_decomp)
   call addfld('CCN2    ','#/cm3   ',pver, 'A','CCN concentration at S=0.05%',phys_decomp)
   call addfld('CCN3    ','#/cm3   ',pver, 'A','CCN concentration at S=0.1%',phys_decomp)
   call addfld('CCN4    ','#/cm3   ',pver, 'A','CCN concentration at S=0.2%',phys_decomp)
   call addfld('CCN5    ','#/cm3   ',pver, 'A','CCN concentration at S=0.5%',phys_decomp)
   call addfld('CCN6    ','#/cm3   ',pver, 'A','CCN concentration at S=1.0%',phys_decomp)


   call addfld('WTKE     ', 'm/s     ', pver, 'A', 'Standard deviation of updraft velocity', phys_decomp)
   call addfld('NDROPMIX ', '#/kg/s  ', pver, 'A', 'Droplet number mixing',                  phys_decomp)
   call addfld('NDROPSRC ', '#/kg/s  ', pver, 'A', 'Droplet number source',                  phys_decomp)
   call addfld('NDROPSNK ', '#/kg/s  ', pver, 'A', 'Droplet number loss by microphysics',    phys_decomp)
   call addfld('NDROPCOL ', '#/m2    ', 1,    'A', 'Column droplet number',                  phys_decomp)

   ! set the add_default fields  
   if (history_amwg) then
      call add_default('CCN3', 1, ' ')
   endif

   if (history_aerosol .and. prog_modal_aero) then
     do m = 1, ntot_amode
        do l = 0, nspec_amode(m)   ! loop over number + chem constituents
           mm = mam_idx(m,l)
           if (l == 0) then   ! number
              call rad_cnst_get_info(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
           else
              call rad_cnst_get_info(0, m, l, spec_name=tmpname, spec_name_cw=tmpname_cw)
           end if
           fieldname(mm)    = trim(tmpname) // '_mixnuc1'
           fieldname_cw(mm) = trim(tmpname_cw) // '_mixnuc1'
        end do
     end do
   endif



end subroutine ndrop_init

!===============================================================================

subroutine ndrop_mode_props_finalize(nmode, sigmag, dgnumlo, dgnumhi, alogsig_out, &
                                     exp45logsig_out, f1_out, f2_out, voltonumblo, voltonumbhi)

   use iso_c_binding, only: c_int64_t, c_loc

   integer,  intent(in) :: nmode
   real(r8), intent(in), target, contiguous  :: sigmag(:)
   real(r8), intent(in), target, contiguous  :: dgnumlo(:)
   real(r8), intent(in), target, contiguous  :: dgnumhi(:)
   real(r8), intent(out), target, contiguous :: alogsig_out(:)
   real(r8), intent(out), target, contiguous :: exp45logsig_out(:)
   real(r8), intent(out), target, contiguous :: f1_out(:)
   real(r8), intent(out), target, contiguous :: f2_out(:)
   real(r8), intent(out), target, contiguous :: voltonumblo(:)
   real(r8), intent(out), target, contiguous :: voltonumbhi(:)

   call ndrop_init_props_select_impl()

   if (use_native_ndrop_init_props_impl) then
      call ndrop_mode_props_finalize_native(nmode, sigmag, dgnumlo, dgnumhi, alogsig_out, &
                                           exp45logsig_out, f1_out, f2_out, voltonumblo, voltonumbhi)
      return
   end if

   call ndrop_init_props_proof_once()
   call ndrop_mode_props_finalize_codon(int(nmode, c_int64_t), pi, c_loc(sigmag(1)), &
        c_loc(dgnumlo(1)), c_loc(dgnumhi(1)), c_loc(alogsig_out(1)), c_loc(exp45logsig_out(1)), &
        c_loc(f1_out(1)), c_loc(f2_out(1)), c_loc(voltonumblo(1)), c_loc(voltonumbhi(1)))

end subroutine ndrop_mode_props_finalize

!===============================================================================

subroutine ndrop_mode_props_finalize_native(nmode, sigmag, dgnumlo, dgnumhi, alogsig_out, &
                                            exp45logsig_out, f1_out, f2_out, voltonumblo, voltonumbhi)

   integer,  intent(in) :: nmode
   real(r8), intent(in)  :: sigmag(:)
   real(r8), intent(in)  :: dgnumlo(:)
   real(r8), intent(in)  :: dgnumhi(:)
   real(r8), intent(out) :: alogsig_out(:)
   real(r8), intent(out) :: exp45logsig_out(:)
   real(r8), intent(out) :: f1_out(:)
   real(r8), intent(out) :: f2_out(:)
   real(r8), intent(out) :: voltonumblo(:)
   real(r8), intent(out) :: voltonumbhi(:)

   integer :: m

   do m = 1, nmode
      alogsig_out(m)     = log(sigmag(m))
      exp45logsig_out(m) = exp(4.5_r8*alogsig_out(m)*alogsig_out(m))
      f1_out(m)          = 0.5_r8*exp(2.5_r8*alogsig_out(m)*alogsig_out(m))
      f2_out(m)          = 1._r8 + 0.25_r8*alogsig_out(m)

      voltonumblo(m) = 1._r8 / ( (pi/6._r8)*                          &
                        (dgnumlo(m)**3._r8)*exp(4.5_r8*alogsig_out(m)**2._r8) )
      voltonumbhi(m) = 1._r8 / ( (pi/6._r8)*                          &
                        (dgnumhi(m)**3._r8)*exp(4.5_r8*alogsig_out(m)**2._r8) )
   end do

end subroutine ndrop_mode_props_finalize_native

!===============================================================================

subroutine dropmixnuc( &
   state, ptend, dtmicro, pbuf, wsub, &
   cldn, cldo, tendnd, factnum)

   ! vertical diffusion and nucleation of cloud droplets
   ! assume cloud presence controlled by cloud fraction
   ! doesn't distinguish between warm, cold clouds

   use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr, c_ptr

   ! arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: dtmicro     ! time step for microphysics (s)

   type(physics_buffer_desc), pointer :: pbuf(:)

   ! arguments
   real(r8), target, intent(in) :: wsub(pcols,pver) ! subgrid vertical velocity
   real(r8), target, intent(in) :: cldn(pcols,pver) ! cloud fraction
   real(r8), target, intent(in) :: cldo(pcols,pver)    ! cloud fraction on previous time step

   ! output arguments
   real(r8), target, intent(out) :: tendnd(pcols,pver) ! change in droplet number concentration (#/kg/s)
   real(r8), target, intent(out) :: factnum(:,:,:)     ! activation fraction for aerosol number
   !--------------------Local storage-------------------------------------

   integer  :: lchnk               ! chunk identifier
   integer  :: ncol                ! number of columns

   real(r8), pointer :: ncldwtr(:,:) ! droplet number concentration (#/kg)
   real(r8), pointer :: temp(:,:)    ! temperature (K)
   real(r8), pointer :: omega(:,:)   ! vertical velocity (Pa/s)
   real(r8), pointer :: pmid(:,:)    ! mid-level pressure (Pa)
   real(r8), pointer :: pint(:,:)    ! pressure at layer interfaces (Pa)
   real(r8), pointer :: pdel(:,:)    ! pressure thickess of layer (Pa)
   real(r8), pointer :: rpdel(:,:)   ! inverse of pressure thickess of layer (/Pa)
   real(r8), pointer :: zm(:,:)      ! geopotential height of level (m)

   real(r8), pointer :: kvh(:,:)     ! vertical diffusivity (m2/s)

   type(ptr2d_t), allocatable :: raer(:)     ! aerosol mass, number mixing ratios
   type(ptr2d_t), allocatable :: qqcw(:)
   type(c_ptr), target :: raer_ptrs(ncnst_tot)
   type(c_ptr), target :: qqcw_ptrs(ncnst_tot)
   real(r8), target :: raertend(pver)  ! tendency of aerosol mass, number mixing ratios
   real(r8), target :: qqcwtend(pver)  ! tendency of cloudborne aerosol mass, number mixing ratios


   real(r8), parameter :: zkmin = 0.01_r8, zkmax = 100._r8
   real(r8), parameter :: wmixmin = 0.1_r8        ! minimum turbulence vertical velocity (m/s)
   real(r8) :: sq2pi

   integer  :: i, k, l, m, mm, n
   integer  :: km1, kp1
   integer, target :: nnew, nsav
   integer  :: ntemp
   integer  :: lptr
   integer, target :: nsubmix, nsubmix_bnd
   integer, save, target :: count_submix(100)
   integer  :: phase ! phase of aerosol

   real(r8) :: arg
   real(r8) :: dtinv
   real(r8) :: dtmin, tinv, dtt

   real(r8), target :: zs(pver) ! inverse of distance between levels (m)
   real(r8), target :: qcld(pver) ! cloud droplet number mixing ratio (#/kg)
   real(r8), target :: qncld(pver)     ! droplet number nucleated on cloud boundaries
   real(r8), target :: srcn(pver)       ! droplet source rate (/s)
   real(r8), target :: cs(pcols,pver)      ! air density (kg/m3)
   real(r8), target :: csbot(pver)       ! air density at bottom (interface) of layer (kg/m3)
   real(r8), target :: csbot_cscen(pver) ! csbot(i)/cs(i,k)
   real(r8), target :: dz(pcols,pver)      ! geometric thickness of layers (m)

   real(r8), target :: wtke(pcols,pver)     ! turbulent vertical velocity at base of layer k (m/s)
   real(r8), target :: wtke_cen(pcols,pver) ! turbulent vertical velocity at center of layer k (m/s)
   real(r8) :: wbar, wmix, wmin, wmax

   real(r8), target :: zn(pver)   ! g/pdel (m2/g) for layer
   real(r8) :: flxconv    ! convergence of flux into lowest layer

   real(r8) :: wdiab           ! diabatic vertical velocity
   real(r8), target :: ekd(pver)       ! diffusivity for droplets (m2/s)
   real(r8), target :: ekk(0:pver)     ! density*diffusivity for droplets (kg/m3 m2/s)
   real(r8), target :: ekkp(pver)      ! zn*zs*density*diffusivity
   real(r8), target :: ekkm(pver)      ! zn*zs*density*diffusivity

   real(r8) :: dum, dumc
   real(r8) :: tmpa
   real(r8) :: dact
   real(r8) :: fluxntot         ! (#/cm2/s)
   real(r8), target :: dtmix
   real(r8) :: alogarg
   real(r8), target :: overlapp(pver), overlapm(pver) ! cloud overlap

   real(r8), target :: nsource(pcols,pver)    ! droplet number source (#/kg/s)
   real(r8), target :: ndropmix(pcols,pver)   ! droplet number mixing (#/kg/s)
   real(r8), target :: ndropcol(pcols)        ! column droplet number (#/m2)
   real(r8) :: cldo_tmp, cldn_tmp
   real(r8), target :: coltend_tmp, coltend_cw_tmp
   real(r8) :: tau_cld_regenerate
   real(r8), target :: taumix_internal_pver_inv ! 1/(internal mixing time scale for k=pver) (1/s)


   real(r8), allocatable, target :: nact(:,:) ! fractional aero. number  activation rate (/s)
   real(r8), allocatable, target :: mact(:,:) ! fractional aero. mass    activation rate (/s)

   real(r8), allocatable, target :: raercol(:,:,:)    ! single column of aerosol mass, number mixing ratios
   real(r8), allocatable, target :: raercol_cw(:,:,:) ! same as raercol but for cloud-borne phase


   real(r8) :: na(pcols), va(pcols), hy(pcols)
   real(r8), allocatable, target :: naermod(:)  ! (1/m3)
   real(r8), allocatable, target :: hygro(:)    ! hygroscopicity of aerosol mode
   real(r8), allocatable, target :: vaerosol(:) ! interstit+activated aerosol volume conc (cm3/cm3)
   real(r8), allocatable, target :: dropmixnuc_species_specdens(:,:)
   real(r8), allocatable, target :: dropmixnuc_species_spechygro(:,:)

   real(r8), target :: source(pver)

   real(r8), allocatable, target :: fn(:)      ! activation fraction for aerosol number
   real(r8), allocatable, target :: fm(:)      ! activation fraction for aerosol mass

   real(r8), allocatable, target :: fluxn(:)   ! number  activation fraction flux (cm/s)
   real(r8), allocatable, target :: fluxm(:)   ! mass    activation fraction flux (cm/s)
   real(r8), target      :: flux_fullact(pver) ! 100%    activation fraction flux (cm/s)
   !     note:  activation fraction fluxes are defined as 
   !     fluxn = [flux of activated aero. number into cloud (#/cm2/s)]
   !           / [aero. number conc. in updraft, just below cloudbase (#/cm3)]


   real(r8), allocatable, target :: coltend(:,:)    ! column tendency for diagnostic output
   real(r8), allocatable, target :: coltend_cw(:,:) ! column tendency
   real(r8) :: ccn(pcols,pver,psat)    ! number conc of aerosols activated at supersat
   real(r8), target :: es_act(pver)
   real(r8), target :: qs_act(pver)
   real(r8), target :: cldn_regen(pver)
   real(r8), allocatable, target :: act_zeta(:)
   real(r8), allocatable, target :: act_eta(:)
   real(r8), allocatable, target :: act_etafactor2(:)
   real(r8), allocatable, target :: act_sqrtg(:)
   real(r8), allocatable, target :: act_amcube(:)
   real(r8), allocatable, target :: act_smc(:)
   real(r8), allocatable, target :: act_lnsm(:)
   real(r8), allocatable, target :: act_sumflxn(:)
   real(r8), allocatable, target :: act_sumflxm(:)
   real(r8), allocatable, target :: act_sumfn(:)
   real(r8), allocatable, target :: act_sumfm(:)
   real(r8), allocatable, target :: act_fnold(:)
   real(r8), allocatable, target :: act_fmold(:)
   integer(c_int64_t) :: activation_status_c

   !-------------------------------------------------------------------------------

   sq2pi = sqrt(2._r8*pi)

   lchnk = state%lchnk
   ncol  = state%ncol

   ncldwtr  => state%q(:,:,numliq_idx)
   temp     => state%t
   omega    => state%omega
   pmid     => state%pmid
   pint     => state%pint
   pdel     => state%pdel
   rpdel    => state%rpdel
   zm       => state%zm

   call pbuf_get_field(pbuf, kvh_idx, kvh)



   arg = 1.0_r8
   if (abs(0.8427_r8 - erf(arg))/0.8427_r8 > 0.001_r8) then
      write(iulog,*) 'erf(1.0) = ',ERF(arg)
      call endrun('dropmixnuc: Error function error')
   endif
   arg = 0.0_r8
   if (erf(arg) /= 0.0_r8) then
      write(iulog,*) 'erf(0.0) = ',erf(arg)
      write(iulog,*) 'dropmixnuc: Error function error'
      call endrun('dropmixnuc: Error function error')
   endif

   dtinv = 1._r8/dtmicro

   allocate( &
      nact(pver,ntot_amode),          &
      mact(pver,ntot_amode),          &
      raer(ncnst_tot),                &
      qqcw(ncnst_tot),                &
      raercol(pver,ncnst_tot,2),      &
      raercol_cw(pver,ncnst_tot,2),   &
      coltend(pcols,ncnst_tot),       &
      coltend_cw(pcols,ncnst_tot),    &
      naermod(ntot_amode),            &
      hygro(ntot_amode),              &
      vaerosol(ntot_amode),           &
      dropmixnuc_species_specdens(ntot_amode,0:size(mam_idx,2)-1), &
      dropmixnuc_species_spechygro(ntot_amode,0:size(mam_idx,2)-1), &
      fn(ntot_amode),                 &
      fm(ntot_amode),                 &
      fluxn(ntot_amode),              &
      fluxm(ntot_amode),              &
      act_zeta(ntot_amode),           &
      act_eta(ntot_amode),            &
      act_etafactor2(ntot_amode),     &
      act_sqrtg(ntot_amode),          &
      act_amcube(ntot_amode),         &
      act_smc(ntot_amode),            &
      act_lnsm(ntot_amode),           &
      act_sumflxn(ntot_amode),        &
      act_sumflxm(ntot_amode),        &
      act_sumfn(ntot_amode),          &
      act_sumfm(ntot_amode),          &
      act_fnold(ntot_amode),          &
      act_fmold(ntot_amode)           )

   dropmixnuc_species_specdens = 0._r8
   dropmixnuc_species_spechygro = 0._r8

   ! Init pointers to mode number and specie mass mixing ratios in 
   ! intersitial and cloud borne phases.
   raer_ptrs = c_null_ptr
   qqcw_ptrs = c_null_ptr
   do m = 1, ntot_amode
      mm = mam_idx(m, 0)
      call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, raer(mm)%fld)
      call rad_cnst_get_mode_num(0, m, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
      raer_ptrs(mm) = c_loc(raer(mm)%fld(1,1))
      qqcw_ptrs(mm) = c_loc(qqcw(mm)%fld(1,1))
      do l = 1, nspec_amode(m)
         mm = mam_idx(m, l)
         call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, raer(mm)%fld)
         call rad_cnst_get_aer_mmr(0, m, l, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
         call rad_cnst_get_aer_props(0, m, l, density_aer=dropmixnuc_species_specdens(m,l), &
              hygro_aer=dropmixnuc_species_spechygro(m,l))
         raer_ptrs(mm) = c_loc(raer(mm)%fld(1,1))
         qqcw_ptrs(mm) = c_loc(qqcw(mm)%fld(1,1))
      end do
   end do

   call ndrop_dropmixnuc_helpers_select_impl()
   if (use_native_ndrop_dropmixnuc_helpers_impl) then
      factnum = 0._r8
      wtke    = 0._r8
   else
      call ndrop_dropmixnuc_helpers_proof_once()
      call ndrop_dropmixnuc_zero_fields_codon(int(pcols, c_int64_t), int(pver, c_int64_t), &
           int(ntot_amode, c_int64_t), c_loc(factnum(1,1,1)), c_loc(wtke(1,1)))
   end if

   if (prog_modal_aero) then
      ! aerosol tendencies
      call physics_ptend_init(ptend, state%psetcols, 'ndrop', lq=lq)
   else
      ! no aerosol tendencies
      call physics_ptend_init(ptend, state%psetcols, 'ndrop')
   end if

   ! overall_main_i_loop
   do i = 1, ncol

      if (use_native_ndrop_dropmixnuc_helpers_impl) then
         do k = top_lev, pver-1
            zs(k) = 1._r8/(zm(i,k) - zm(i,k+1))
         end do
         zs(pver) = zs(pver-1)

         ! load number nucleated into qcld on cloud boundaries

         do k = top_lev, pver

            qcld(k)  = ncldwtr(i,k)
            qncld(k) = 0._r8
            srcn(k)  = 0._r8
            cs(i,k)  = pmid(i,k)/(rair*temp(i,k))        ! air density (kg/m3)
            dz(i,k)  = 1._r8/(cs(i,k)*gravit*rpdel(i,k)) ! layer thickness in m

            do m = 1, ntot_amode
               nact(k,m) = 0._r8
               mact(k,m) = 0._r8
            end do

            zn(k) = gravit*rpdel(i,k)

            if (k < pver) then
               ekd(k)   = kvh(i,k+1)
               ekd(k)   = max(ekd(k), zkmin)
               ekd(k)   = min(ekd(k), zkmax)
               csbot(k) = 2.0_r8*pint(i,k+1)/(rair*(temp(i,k) + temp(i,k+1)))
               csbot_cscen(k) = csbot(k)/cs(i,k)
            else
               ekd(k)   = 0._r8
               csbot(k) = cs(i,k)
               csbot_cscen(k) = 1.0_r8
            end if

            ! rce-comment - define wtke at layer centers for new-cloud activation
            !    and at layer boundaries for old-cloud activation
            !++ag
            wtke_cen(i,k) = wsub(i,k)
            wtke(i,k)     = wsub(i,k)
            !--ag
            wtke_cen(i,k) = max(wtke_cen(i,k), wmixmin)
            wtke(i,k)     = max(wtke(i,k), wmixmin)

            nsource(i,k) = 0._r8

         end do
      else
         call ndrop_dropmixnuc_helpers_proof_once()
         call ndrop_dropmixnuc_column_init_codon(int(i, c_int64_t), int(pcols, c_int64_t), &
              int(pver, c_int64_t), int(top_lev, c_int64_t), int(ntot_amode, c_int64_t), &
              gravit, rair, zkmin, zkmax, wmixmin, c_loc(ncldwtr(1,1)), c_loc(temp(1,1)), &
              c_loc(pmid(1,1)), c_loc(pint(1,1)), c_loc(rpdel(1,1)), c_loc(zm(1,1)), &
              c_loc(kvh(1,1)), c_loc(wsub(1,1)), c_loc(qcld(1)), c_loc(qncld(1)), &
              c_loc(srcn(1)), c_loc(cs(1,1)), c_loc(dz(1,1)), c_loc(nact(1,1)), &
              c_loc(mact(1,1)), c_loc(zn(1)), c_loc(ekd(1)), c_loc(csbot(1)), &
              c_loc(csbot_cscen(1)), c_loc(wtke_cen(1,1)), c_loc(wtke(1,1)), &
              c_loc(nsource(1,1)), c_loc(zs(1)))
      end if

      nsav = 1
      nnew = 2
      if (use_native_ndrop_dropmixnuc_helpers_impl) then
         do m = 1, ntot_amode
            mm = mam_idx(m,0)
            raercol_cw(:,mm,nsav) = 0.0_r8
            raercol(:,mm,nsav)    = 0.0_r8
            raercol_cw(top_lev:pver,mm,nsav) = qqcw(mm)%fld(i,top_lev:pver)
            raercol(top_lev:pver,mm,nsav)    = raer(mm)%fld(i,top_lev:pver)
            do l = 1, nspec_amode(m)
               mm = mam_idx(m,l)
               raercol_cw(top_lev:pver,mm,nsav) = qqcw(mm)%fld(i,top_lev:pver)
               raercol(top_lev:pver,mm,nsav)    = raer(mm)%fld(i,top_lev:pver)
            end do
         end do
      else
         call ndrop_dropmixnuc_helpers_proof_once()
         call ndrop_dropmixnuc_aero_column_copy_all_codon(int(i, c_int64_t), int(pcols, c_int64_t), &
              int(pver, c_int64_t), int(top_lev, c_int64_t), int(ntot_amode, c_int64_t), &
              int(ncnst_tot, c_int64_t), int(nsav, c_int64_t), c_loc(raer_ptrs(1)), &
              c_loc(qqcw_ptrs(1)), c_loc(nspec_amode(1)), c_loc(mam_idx(1,0)), &
              c_loc(raercol(1,1,1)), c_loc(raercol_cw(1,1,1)))
      end if

      ! droplet nucleation/aerosol activation

      ! tau_cld_regenerate = time scale for regeneration of cloudy air 
      !    by (horizontal) exchange with clear air
      tau_cld_regenerate = 3600.0_r8 * 3.0_r8 

      if (.not. use_native_ndrop_dropmixnuc_helpers_impl) then
         do k = top_lev, pver
            call qsat(temp(i,k), rair*cs(i,k)*temp(i,k), es_act(k), qs_act(k))
         end do
         taumix_internal_pver_inv = 0.0_r8
         call ndrop_dropmixnuc_helpers_proof_once()
         activation_status_c = ndrop_dropmixnuc_activation_loops_codon(int(i, c_int64_t), &
              int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
              int(ntot_amode, c_int64_t), int(ncnst_tot, c_int64_t), int(nsav, c_int64_t), &
              dtmicro, dtinv, rair, 1013.25e2_r8, t0, rhoh2o, latvap, cpair, rh2o, gravit, &
              pi, aten, twothird, sq2, sqpi, sixth, zero, c_loc(cldn(1,1)), c_loc(cldo(1,1)), &
              c_loc(cldn_regen(1)), c_loc(temp(1,1)), c_loc(cs(1,1)), c_loc(dz(1,1)), &
              c_loc(wtke(1,1)), c_loc(wtke_cen(1,1)), c_loc(zs(1)), c_loc(ekd(1)), &
              c_loc(csbot_cscen(1)), c_loc(qs_act(1)), c_loc(qcld(1)), c_loc(srcn(1)), &
              c_loc(nsource(1,1)), c_loc(factnum(1,1,1)), c_loc(nact(1,1)), c_loc(mact(1,1)), &
              c_loc(taumix_internal_pver_inv), c_loc(raer_ptrs(1)), c_loc(qqcw_ptrs(1)), &
              c_loc(nspec_amode(1)), c_loc(mam_idx(1,0)), c_loc(dropmixnuc_species_specdens(1,0)), &
              c_loc(dropmixnuc_species_spechygro(1,0)), c_loc(voltonumblo_amode(1)), &
              c_loc(voltonumbhi_amode(1)), c_loc(alogsig(1)), c_loc(exp45logsig(1)), &
              c_loc(f1(1)), c_loc(f2(1)), c_loc(raercol(1,1,1)), c_loc(raercol_cw(1,1,1)), &
              c_loc(naermod(1)), c_loc(vaerosol(1)), c_loc(hygro(1)), c_loc(fn(1)), &
              c_loc(fm(1)), c_loc(fluxn(1)), c_loc(fluxm(1)), c_loc(flux_fullact(1)), &
              c_loc(act_zeta(1)), c_loc(act_eta(1)), c_loc(act_etafactor2(1)), &
              c_loc(act_sqrtg(1)), c_loc(act_amcube(1)), c_loc(act_smc(1)), c_loc(act_lnsm(1)), &
              c_loc(act_sumflxn(1)), c_loc(act_sumflxm(1)), c_loc(act_sumfn(1)), c_loc(act_sumfm(1)), &
              c_loc(act_fnold(1)), c_loc(act_fmold(1)))
         if (activation_status_c == 1_c_int64_t) then
            call endrun('dropmixnuc: Codon activation integration loop did not converge')
         else if (activation_status_c == 2_c_int64_t) then
            call endrun('dropmixnuc: Codon activation fraction exceeded one')
         else if (activation_status_c /= 0_c_int64_t) then
            call endrun('dropmixnuc: Codon activation loops returned unknown status')
         end if
         go to 220
      end if

      ! k-loop for growing/shrinking cloud calcs .............................
      ! grow_shrink_main_k_loop: &
      do k = top_lev, pver

         ! shrinking cloud ......................................................
         !    treat the reduction of cloud fraction from when cldn(i,k) < cldo(i,k)
         !    and also dissipate the portion of the cloud that will be regenerated
         cldo_tmp = cldo(i,k)
         cldn_tmp = cldn(i,k) * exp( -dtmicro/tau_cld_regenerate )
         !    alternate formulation
         !    cldn_tmp = cldn(i,k) * max( 0.0_r8, (1.0_r8-dtmicro/tau_cld_regenerate) )

         if (cldn_tmp < cldo_tmp) then
            if (use_native_ndrop_dropmixnuc_helpers_impl) then
               !  droplet loss in decaying cloud
               !++ sungsup
               nsource(i,k) = nsource(i,k) + qcld(k)*(cldn_tmp - cldo_tmp)/cldo_tmp*dtinv
               qcld(k)      = qcld(k)*(1._r8 + (cldn_tmp - cldo_tmp)/cldo_tmp)
               !-- sungsup

               ! convert activated aerosol to interstitial in decaying cloud

               dumc = (cldn_tmp - cldo_tmp)/cldo_tmp
               do m = 1, ntot_amode
                  mm = mam_idx(m,0)
                  dact   = raercol_cw(k,mm,nsav)*dumc
                  raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact   ! cloud-borne aerosol
                  raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
                  do l = 1, nspec_amode(m)
                     mm = mam_idx(m,l)
                     dact    = raercol_cw(k,mm,nsav)*dumc
                     raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
                     raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
                  end do
               end do
            else
               call ndrop_dropmixnuc_helpers_proof_once()
               call ndrop_dropmixnuc_shrink_cloud_codon(int(i, c_int64_t), int(k, c_int64_t), &
                    int(pcols, c_int64_t), int(pver, c_int64_t), int(ntot_amode, c_int64_t), &
                    int(ncnst_tot, c_int64_t), int(nsav, c_int64_t), dtinv, cldn_tmp, cldo_tmp, &
                    c_loc(qcld(1)), c_loc(nsource(1,1)), c_loc(nspec_amode(1)), &
                    c_loc(mam_idx(1,0)), c_loc(raercol(1,1,1)), c_loc(raercol_cw(1,1,1)))
            end if
         end if

         ! growing cloud ......................................................
         !    treat the increase of cloud fraction from when cldn(i,k) > cldo(i,k)
         !    and also regenerate part of the cloud 
         cldo_tmp = cldn_tmp
         cldn_tmp = cldn(i,k)

         if (cldn_tmp-cldo_tmp > 0.01_r8) then

            ! rce-comment - use wtke at layer centers for new-cloud activation
            wbar  = wtke_cen(i,k)
            wmix  = 0._r8
            wmin  = 0._r8
            wmax  = 10._r8
            wdiab = 0

            ! load aerosol properties, assuming external mixtures

            phase = 1 ! interstitial
            do m = 1, ntot_amode
               call loadaer( &
                  state, pbuf, i, i, k, &
                  m, cs, phase, na, va, &
                  hy)
               naermod(m)  = na(i)
               vaerosol(m) = va(i)
               hygro(m)    = hy(i)
            end do

            call activate_modal( &
               wbar, wmix, wdiab, wmin, wmax,                       &
               temp(i,k), cs(i,k), naermod, ntot_amode, &
               vaerosol, hygro, fn, fm, fluxn,                      &
               fluxm,flux_fullact(k))

            dumc = (cldn_tmp - cldo_tmp)
            if (use_native_ndrop_dropmixnuc_helpers_impl) then
               factnum(i,k,:) = fn
               do m = 1, ntot_amode
                  mm = mam_idx(m,0)
                  dact   = dumc*fn(m)*raer(mm)%fld(i,k) ! interstitial only
                  qcld(k) = qcld(k) + dact
                  nsource(i,k) = nsource(i,k) + dact*dtinv
                  raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
                  raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
                  dum = dumc*fm(m)
                  do l = 1, nspec_amode(m)
                     mm = mam_idx(m,l)
                     dact    = dum*raer(mm)%fld(i,k) ! interstitial only
                     raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
                     raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
                  enddo
               enddo
            else
               call ndrop_dropmixnuc_helpers_proof_once()
               call ndrop_dropmixnuc_grow_cloud_update_all_codon(int(i, c_int64_t), int(k, c_int64_t), &
                    int(pcols, c_int64_t), int(pver, c_int64_t), int(ntot_amode, c_int64_t), &
                    int(ncnst_tot, c_int64_t), int(nsav, c_int64_t), dtinv, dumc, c_loc(raer_ptrs(1)), &
                    c_loc(nspec_amode(1)), c_loc(mam_idx(1,0)), c_loc(fn(1)), c_loc(fm(1)), &
                    c_loc(qcld(1)), c_loc(nsource(1,1)), c_loc(raercol(1,1,1)), &
                    c_loc(raercol_cw(1,1,1)), c_loc(factnum(1,1,1)))
            end if
         endif

      enddo  ! grow_shrink_main_k_loop
      ! end of k-loop for growing/shrinking cloud calcs ......................

      ! ......................................................................
      ! start of k-loop for calc of old cloud activation tendencies ..........
      !
      ! rce-comment
      !    changed this part of code to use current cloud fraction (cldn) exclusively
      !    consider case of cldo(:)=0, cldn(k)=1, cldn(k+1)=0
      !    previous code (which used cldo below here) would have no cloud-base activation
      !       into layer k.  however, activated particles in k mix out to k+1,
      !       so they are incorrectly depleted with no replacement

      ! old_cloud_main_k_loop
      do k = top_lev, pver
         kp1 = min0(k+1, pver)
         taumix_internal_pver_inv = 0.0_r8

         if (cldn(i,k) > 0.01_r8) then

            wdiab = 0
            wmix  = 0._r8                       ! single updraft
            wbar  = wtke(i,k)                   ! single updraft
            if (k == pver) wbar = wtke_cen(i,k) ! single updraft
            wmax  = 10._r8
            wmin  = 0._r8

            if (cldn(i,k) - cldn(i,kp1) > 0.01_r8 .or. k == pver) then

               ! cloud base

               ! ekd(k) = wtke(i,k)*dz(i,k)/sq2pi
               ! rce-comments
               !   first, should probably have 1/zs(k) here rather than dz(i,k) because
               !      the turbulent flux is proportional to ekd(k)*zs(k),
               !      while the dz(i,k) is used to get flux divergences
               !      and mixing ratio tendency/change
               !   second and more importantly, using a single updraft velocity here
               !      means having monodisperse turbulent updraft and downdrafts.
               !      The sq2pi factor assumes a normal draft spectrum.
               !      The fluxn/fluxm from activate must be consistent with the
               !      fluxes calculated in explmix.
               ekd(k) = wbar/zs(k)

               alogarg = max(1.e-20_r8, 1/cldn(i,k) - 1._r8)
               wmin    = wbar + wmix*0.25_r8*sq2pi*log(alogarg)
               phase   = 1   ! interstitial

               do m = 1, ntot_amode
                  ! rce-comment - use kp1 here as old-cloud activation involves 
                  !   aerosol from layer below
                  call loadaer( &
                     state, pbuf, i, i, kp1,  &
                     m, cs, phase, na, va,   &
                     hy)
                  naermod(m)  = na(i)
                  vaerosol(m) = va(i)
                  hygro(m)    = hy(i)
               end do

               call activate_modal( &
                  wbar, wmix, wdiab, wmin, wmax,                       &
                  temp(i,k), cs(i,k), naermod, ntot_amode, &
                  vaerosol, hygro, fn, fm, fluxn,                      &
                  fluxm, flux_fullact(k))

               if (use_native_ndrop_dropmixnuc_helpers_impl) then
                  factnum(i,k,:) = fn
               else
                  call ndrop_dropmixnuc_helpers_proof_once()
                  call ndrop_dropmixnuc_factnum_store_codon(int(i, c_int64_t), int(k, c_int64_t), &
                       int(pcols, c_int64_t), int(pver, c_int64_t), int(ntot_amode, c_int64_t), &
                       c_loc(fn(1)), c_loc(factnum(1,1,1)))
               end if

               if (k < pver) then
                  dumc = cldn(i,k) - cldn(i,kp1)
               else
                  dumc = cldn(i,k)
               endif

               fluxntot = 0

               ! rce-comment 1
               !    flux of activated mass into layer k (in kg/m2/s)
               !       = "actmassflux" = dumc*fluxm*raercol(kp1,lmass)*csbot(k)
               !    source of activated mass (in kg/kg/s) = flux divergence
               !       = actmassflux/(cs(i,k)*dz(i,k))
               !    so need factor of csbot_cscen = csbot(k)/cs(i,k)
               !                   dum=1./(dz(i,k))
               dum=csbot_cscen(k)/(dz(i,k))

               ! rce-comment 2
               !    code for k=pver was changed to use the following conceptual model
               !    in k=pver, there can be no cloud-base activation unless one considers
               !       a scenario such as the layer being partially cloudy, 
               !       with clear air at bottom and cloudy air at top
               !    assume this scenario, and that the clear/cloudy portions mix with 
               !       a timescale taumix_internal = dz(i,pver)/wtke_cen(i,pver)
               !    in the absence of other sources/sinks, qact (the activated particle 
               !       mixratio) attains a steady state value given by
               !          qact_ss = fcloud*fact*qtot
               !       where fcloud is cloud fraction, fact is activation fraction, 
               !       qtot=qact+qint, qint is interstitial particle mixratio
               !    the activation rate (from mixing within the layer) can now be
               !       written as
               !          d(qact)/dt = (qact_ss - qact)/taumix_internal
               !                     = qtot*(fcloud*fact*wtke/dz) - qact*(wtke/dz)
               !    note that (fcloud*fact*wtke/dz) is equal to the nact/mact
               !    also, d(qact)/dt can be negative.  in the code below
               !       it is forced to be >= 0
               !
               ! steve -- 
               !    you will likely want to change this.  i did not really understand 
               !       what was previously being done in k=pver
               !    in the cam3_5_3 code, wtke(i,pver) appears to be equal to the
               !       droplet deposition velocity which is quite small
               !    in the cam3_5_37 version, wtke is done differently and is much
               !       larger in k=pver, so the activation is stronger there
               !
               if (k == pver) then
                  taumix_internal_pver_inv = flux_fullact(k)/dz(i,k)
               end if

               if (use_native_ndrop_dropmixnuc_helpers_impl) then
                  do m = 1, ntot_amode
                     mm = mam_idx(m,0)
                     fluxn(m) = fluxn(m)*dumc
                     fluxm(m) = fluxm(m)*dumc
                     nact(k,m) = nact(k,m) + fluxn(m)*dum
                     mact(k,m) = mact(k,m) + fluxm(m)*dum
                     if (k < pver) then
                        ! note that kp1 is used here
                        fluxntot = fluxntot &
                           + fluxn(m)*raercol(kp1,mm,nsav)*cs(i,k)
                     else
                        tmpa = raercol(kp1,mm,nsav)*fluxn(m) &
                             + raercol_cw(kp1,mm,nsav)*(fluxn(m) &
                             - taumix_internal_pver_inv*dz(i,k))
                        fluxntot = fluxntot + max(0.0_r8, tmpa)*cs(i,k)
                     end if
                  end do
                  srcn(k)      = srcn(k) + fluxntot/(cs(i,k)*dz(i,k))
                  nsource(i,k) = nsource(i,k) + fluxntot/(cs(i,k)*dz(i,k))
               else
                  call ndrop_dropmixnuc_helpers_proof_once()
                  call ndrop_dropmixnuc_old_cloud_activate_update_codon(int(i, c_int64_t), &
                       int(k, c_int64_t), int(kp1, c_int64_t), int(pcols, c_int64_t), &
                       int(pver, c_int64_t), int(ntot_amode, c_int64_t), int(ncnst_tot, c_int64_t), &
                       int(nsav, c_int64_t), dumc, dum, cs(i,k), dz(i,k), taumix_internal_pver_inv, &
                       c_loc(fluxn(1)), c_loc(fluxm(1)), c_loc(nact(1,1)), c_loc(mact(1,1)), &
                       c_loc(mam_idx(1,0)), c_loc(raercol(1,1,1)), c_loc(raercol_cw(1,1,1)), &
                       c_loc(srcn(1)), c_loc(nsource(1,1)))
               end if

            endif  ! (cldn(i,k) - cldn(i,kp1) > 0.01 .or. k == pver)

         else

            ! no cloud

            if (use_native_ndrop_dropmixnuc_helpers_impl) then
               nsource(i,k) = nsource(i,k) - qcld(k)*dtinv
               qcld(k)      = 0

               ! convert activated aerosol to interstitial in decaying cloud

               do m = 1, ntot_amode
                  mm = mam_idx(m,0)
                  raercol(k,mm,nsav)    = raercol(k,mm,nsav) + raercol_cw(k,mm,nsav)  ! cloud-borne aerosol
                  raercol_cw(k,mm,nsav) = 0._r8

                  do l = 1, nspec_amode(m)
                     mm = mam_idx(m,l)
                     raercol(k,mm,nsav)    = raercol(k,mm,nsav) + raercol_cw(k,mm,nsav) ! cloud-borne aerosol
                     raercol_cw(k,mm,nsav) = 0._r8
                  end do
               end do
            else
               call ndrop_dropmixnuc_helpers_proof_once()
               call ndrop_dropmixnuc_clear_old_cloud_codon(int(i, c_int64_t), int(k, c_int64_t), &
                    int(pcols, c_int64_t), int(pver, c_int64_t), int(ntot_amode, c_int64_t), &
                    int(ncnst_tot, c_int64_t), int(nsav, c_int64_t), dtinv, c_loc(qcld(1)), &
                    c_loc(nsource(1,1)), c_loc(nspec_amode(1)), c_loc(mam_idx(1,0)), &
                    c_loc(raercol(1,1,1)), c_loc(raercol_cw(1,1,1)))
            end if
         end if

      end do  ! old_cloud_main_k_loop

220   continue

      ! switch nsav, nnew so that nnew is the updated aerosol
      if (use_native_ndrop_dropmixnuc_helpers_impl) then
         ntemp = nsav
         nsav  = nnew
         nnew  = ntemp
      else
         call ndrop_dropmixnuc_helpers_proof_once()
         call ndrop_dropmixnuc_swap_slots_codon(c_loc(nsav), c_loc(nnew))
      end if

      ! load new droplets in layers above, below clouds

      if (use_native_ndrop_dropmixnuc_helpers_impl) then
         dtmin     = dtmicro
         ekk(top_lev-1)    = 0.0_r8
         ekk(pver) = 0.0_r8
         do k = top_lev, pver-1
            ! rce-comment -- ekd(k) is eddy-diffusivity at k/k+1 interface
            !   want ekk(k) = ekd(k) * (density at k/k+1 interface)
            !   so use pint(i,k+1) as pint is 1:pverp
            !           ekk(k)=ekd(k)*2.*pint(i,k)/(rair*(temp(i,k)+temp(i,k+1)))
            !           ekk(k)=ekd(k)*2.*pint(i,k+1)/(rair*(temp(i,k)+temp(i,k+1)))
            ekk(k) = ekd(k)*csbot(k)
         end do

         do k = top_lev, pver
            km1     = max0(k-1, top_lev)
            ekkp(k) = zn(k)*ekk(k)*zs(k)
            ekkm(k) = zn(k)*ekk(k-1)*zs(km1)
            tinv    = ekkp(k) + ekkm(k)

            ! rce-comment -- tinv is the sum of all first-order-loss-rates
            !    for the layer.  for most layers, the activation loss rate
            !    (for interstitial particles) is accounted for by the loss by
            !    turb-transfer to the layer above.
            !    k=pver is special, and the loss rate for activation within
            !    the layer must be added to tinv.  if not, the time step
            !    can be too big, and explmix can produce negative values.
            !    the negative values are reset to zero, resulting in an
            !    artificial source.
            if (k == pver) tinv = tinv + taumix_internal_pver_inv

            if (tinv .gt. 1.e-6_r8) then
               dtt   = 1._r8/tinv
               dtmin = min(dtmin, dtt)
            end if
         end do

         dtmix   = 0.9_r8*dtmin
         nsubmix = dtmicro/dtmix + 1
         if (nsubmix > 100) then
            nsubmix_bnd = 100
         else
            nsubmix_bnd = nsubmix
         end if
         count_submix(nsubmix_bnd) = count_submix(nsubmix_bnd) + 1
         dtmix = dtmicro/nsubmix

         do k = top_lev, pver
            kp1 = min(k+1, pver)
            km1 = max(k-1, top_lev)
            ! maximum overlap assumption
            if (cldn(i,kp1) > 1.e-10_r8) then
               overlapp(k) = min(cldn(i,k)/cldn(i,kp1), 1._r8)
            else
               overlapp(k) = 1._r8
            end if
            if (cldn(i,km1) > 1.e-10_r8) then
               overlapm(k) = min(cldn(i,k)/cldn(i,km1), 1._r8)
            else
               overlapm(k) = 1._r8
            end if
         end do


         ! rce-comment
         !    the activation source(k) = mact(k,m)*raercol(kp1,lmass)
         !       should not exceed the rate of transfer of unactivated particles
         !       from kp1 to k which = ekkp(k)*raercol(kp1,lmass)
         !    however it might if things are not "just right" in subr activate
         !    the following is a safety measure to avoid negatives in explmix
         do k = top_lev, pver-1
            do m = 1, ntot_amode
               nact(k,m) = min( nact(k,m), ekkp(k) )
               mact(k,m) = min( mact(k,m), ekkp(k) )
            end do
         end do
      else
         call ndrop_dropmixnuc_helpers_proof_once()
         call ndrop_dropmixnuc_mix_setup_codon(int(i, c_int64_t), int(pcols, c_int64_t), &
              int(pver, c_int64_t), int(top_lev, c_int64_t), int(ntot_amode, c_int64_t), &
              dtmicro, taumix_internal_pver_inv, c_loc(cldn(1,1)), c_loc(zs(1)), c_loc(zn(1)), &
              c_loc(csbot(1)), c_loc(ekd(1)), c_loc(nact(1,1)), c_loc(mact(1,1)), c_loc(ekk(0)), &
              c_loc(ekkp(1)), c_loc(ekkm(1)), c_loc(overlapp(1)), c_loc(overlapm(1)), &
              c_loc(count_submix(1)), c_loc(nsubmix), c_loc(dtmix))
      end if


      ! old_cloud_nsubmix_loop
      call ndrop_explmix_select_impl()
      do n = 1, nsubmix
         if (use_native_ndrop_dropmixnuc_helpers_impl .or. use_native_ndrop_explmix_impl) then
            qncld(:) = qcld(:)
            ! switch nsav, nnew so that nsav is the updated aerosol
            ntemp   = nsav
            nsav    = nnew
            nnew    = ntemp
            srcn(:) = 0.0_r8

            do m = 1, ntot_amode
               mm = mam_idx(m,0)

               ! update droplet source
               ! rce-comment- activation source in layer k involves particles from k+1
               !	       srcn(:)=srcn(:)+nact(:,m)*(raercol(:,mm,nsav))
               srcn(top_lev:pver-1) = srcn(top_lev:pver-1) + &
                    nact(top_lev:pver-1,m)*(raercol(top_lev+1:pver,mm,nsav))

               ! rce-comment- new formulation for k=pver
               !              srcn(  pver  )=srcn(  pver  )+nact(  pver  ,m)*(raercol(  pver,mm,nsav))
               tmpa = raercol(pver,mm,nsav)*nact(pver,m) &
                    + raercol_cw(pver,mm,nsav)*(nact(pver,m) - taumix_internal_pver_inv)
               srcn(pver) = srcn(pver) + max(0.0_r8,tmpa)
            end do

            call explmix(  &
               qcld, srcn, ekkp, ekkm, overlapp,  &
               overlapm, qncld, zero, zero, pver, &
               dtmix, .false.)

            ! rce-comment
            !    the interstitial particle mixratio is different in clear/cloudy portions
            !    of a layer, and generally higher in the clear portion.  (we have/had
            !    a method for diagnosing the the clear/cloudy mixratios.)  the activation
            !    source terms involve clear air (from below) moving into cloudy air (above).
            !    in theory, the clear-portion mixratio should be used when calculating
            !    source terms
            do m = 1, ntot_amode
               mm = mam_idx(m,0)
               ! rce-comment -   activation source in layer k involves particles from k+1
               !	              source(:)= nact(:,m)*(raercol(:,mm,nsav))
               source(top_lev:pver-1) = nact(top_lev:pver-1,m)*(raercol(top_lev+1:pver,mm,nsav))
               ! rce-comment - new formulation for k=pver
               !               source(  pver  )= nact(  pver,  m)*(raercol(  pver,mm,nsav))
               tmpa = raercol(pver,mm,nsav)*nact(pver,m) &
                    + raercol_cw(pver,mm,nsav)*(nact(pver,m) - taumix_internal_pver_inv)
               source(pver) = max(0.0_r8, tmpa)

               flxconv = 0._r8

               call explmix( &
                  raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
                  overlapm, raercol_cw(:,mm,nsav), zero, zero, pver,   &
                  dtmix, .false.)

               call explmix( &
                  raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
                  overlapm, raercol(:,mm,nsav), zero, flxconv, pver, &
                  dtmix, .true., raercol_cw(:,mm,nsav))

               do l = 1, nspec_amode(m)
                  mm = mam_idx(m,l)
                  ! rce-comment -   activation source in layer k involves particles from k+1
                  !	          source(:)= mact(:,m)*(raercol(:,mm,nsav))
                  source(top_lev:pver-1) = mact(top_lev:pver-1,m)*(raercol(top_lev+1:pver,mm,nsav))
                  ! rce-comment- new formulation for k=pver
                  !                 source(  pver  )= mact(  pver  ,m)*(raercol(  pver,mm,nsav))
                  tmpa = raercol(pver,mm,nsav)*mact(pver,m) &
                       + raercol_cw(pver,mm,nsav)*(mact(pver,m) - taumix_internal_pver_inv)
                  source(pver) = max(0.0_r8, tmpa)

                  flxconv = 0._r8

                  call explmix( &
                     raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
                     overlapm, raercol_cw(:,mm,nsav), zero, zero, pver,   &
                     dtmix, .false.)

                  call explmix( &
                     raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
                     overlapm, raercol(:,mm,nsav), zero, flxconv, pver, &
                     dtmix, .true., raercol_cw(:,mm,nsav))

               end do
            end do
         else
            call ndrop_dropmixnuc_helpers_proof_once()
            call ndrop_explmix_proof_once()
            call ndrop_dropmixnuc_submix_all_codon(int(pver, c_int64_t), int(top_lev, c_int64_t), &
                 int(ntot_amode, c_int64_t), int(ncnst_tot, c_int64_t), dtmix, taumix_internal_pver_inv, &
                 c_loc(nact(1,1)), c_loc(mact(1,1)), c_loc(mam_idx(1,0)), c_loc(nspec_amode(1)), &
                 c_loc(ekkp(1)), c_loc(ekkm(1)), c_loc(overlapp(1)), c_loc(overlapm(1)), &
                 c_loc(qcld(1)), c_loc(qncld(1)), c_loc(srcn(1)), c_loc(source(1)), &
                 c_loc(raercol(1,1,1)), c_loc(raercol_cw(1,1,1)), c_loc(nsav), c_loc(nnew))
         end if

      end do ! old_cloud_nsubmix_loop

      ! evaporate particles again if no cloud

      if (use_native_ndrop_dropmixnuc_helpers_impl) then
         do k = top_lev, pver
            if (cldn(i,k) == 0._r8) then
               ! no cloud
               qcld(k)=0._r8

               ! convert activated aerosol to interstitial in decaying cloud
               do m = 1, ntot_amode
                  mm = mam_idx(m,0)
                  raercol(k,mm,nnew)    = raercol(k,mm,nnew) + raercol_cw(k,mm,nnew)
                  raercol_cw(k,mm,nnew) = 0._r8

                  do l = 1, nspec_amode(m)
                     mm = mam_idx(m,l)
                     raercol(k,mm,nnew)    = raercol(k,mm,nnew) + raercol_cw(k,mm,nnew)
                     raercol_cw(k,mm,nnew) = 0._r8
                  end do
               end do
            end if
         end do
      else
         call ndrop_dropmixnuc_helpers_proof_once()
         call ndrop_dropmixnuc_evaporate_clear_layers_codon(int(i, c_int64_t), int(pcols, c_int64_t), &
              int(pver, c_int64_t), int(top_lev, c_int64_t), int(ntot_amode, c_int64_t), &
              int(ncnst_tot, c_int64_t), int(nnew, c_int64_t), c_loc(cldn(1,1)), c_loc(qcld(1)), &
              c_loc(nspec_amode(1)), c_loc(mam_idx(1,0)), c_loc(raercol(1,1,1)), &
              c_loc(raercol_cw(1,1,1)))
      end if

      ! droplet number

      if (use_native_ndrop_dropmixnuc_helpers_impl) then
         ndropcol(i) = 0._r8
         do k = top_lev, pver
            ndropmix(i,k) = (qcld(k) - ncldwtr(i,k))*dtinv - nsource(i,k)
            tendnd(i,k)   = (max(qcld(k), 1.e-6_r8) - ncldwtr(i,k))*dtinv
            ndropcol(i)   = ndropcol(i) + ncldwtr(i,k)*pdel(i,k)
         end do
         ndropcol(i) = ndropcol(i)/gravit
      else
         call ndrop_dropmixnuc_helpers_proof_once()
         call ndrop_dropmixnuc_finalize_column_codon(int(i, c_int64_t), int(pcols, c_int64_t), &
              int(pver, c_int64_t), int(top_lev, c_int64_t), dtinv, gravit, c_loc(qcld(1)), &
              c_loc(ncldwtr(1,1)), c_loc(pdel(1,1)), c_loc(nsource(1,1)), c_loc(ndropmix(1,1)), &
              c_loc(tendnd(1,1)), c_loc(ndropcol(1)))
      end if

      if (prog_modal_aero) then

         if (use_native_ndrop_dropmixnuc_helpers_impl) then
            raertend = 0._r8
            qqcwtend = 0._r8

            do m = 1, ntot_amode
               do l = 0, nspec_amode(m)

                  mm   = mam_idx(m,l)
                  lptr = mam_cnst_idx(m,l)

                  raertend(top_lev:pver) = (raercol(top_lev:pver,mm,nnew) - raer(mm)%fld(i,top_lev:pver))*dtinv
                  qqcwtend(top_lev:pver) = (raercol_cw(top_lev:pver,mm,nnew) - qqcw(mm)%fld(i,top_lev:pver))*dtinv

                  coltend(i,mm)    = sum( pdel(i,:)*raertend )/gravit
                  coltend_cw(i,mm) = sum( pdel(i,:)*qqcwtend )/gravit

                  ptend%q(i,:,lptr) = 0.0_r8
                  ptend%q(i,top_lev:pver,lptr) = raertend(top_lev:pver) ! set tendencies for interstitial aerosol

                  qqcw(mm)%fld(i,:) = 0.0_r8
                  qqcw(mm)%fld(i,top_lev:pver) = raercol_cw(top_lev:pver,mm,nnew) ! update cloud-borne aerosol
               end do
            end do
         else
            call ndrop_dropmixnuc_helpers_proof_once()
            call ndrop_dropmixnuc_aero_tend_all_codon_wrap(i, pcols, ptend%psetcols, pver, top_lev, &
                 ntot_amode, ncnst_tot, size(mam_idx,2)-1, nnew, dtinv, gravit, raer_ptrs, qqcw_ptrs, &
                 nspec_amode, mam_idx, mam_cnst_idx, pdel, raercol, raercol_cw, coltend, &
                 coltend_cw, ptend%q)
         end if

      end if

   end do  ! overall_main_i_loop
   ! end of main loop over i/longitude ....................................

   call outfld('NDROPCOL', ndropcol, pcols, lchnk)
   call outfld('NDROPSRC', nsource,  pcols, lchnk)
   call outfld('NDROPMIX', ndropmix, pcols, lchnk)
   call outfld('WTKE    ', wtke,     pcols, lchnk)

   call ccncalc(state, pbuf, cs, ccn)
   do l = 1, psat
      call outfld(ccn_name(l), ccn(1,1,l), pcols, lchnk)
   enddo

   ! do column tendencies
   if (prog_modal_aero) then
      do m = 1, ntot_amode
         do l = 0, nspec_amode(m)
            mm = mam_idx(m,l)
            call outfld(fieldname(mm),    coltend(:,mm),    pcols, lchnk)
            call outfld(fieldname_cw(mm), coltend_cw(:,mm), pcols, lchnk)
         end do
      end do
   end if

   if (.not. use_native_ndrop_dropmixnuc_helpers_impl .and. &
       .not. use_native_ndrop_explmix_impl .and. &
       .not. use_native_ndrop_ccncalc_helpers_impl .and. &
       .not. ndrop_dropmixnuc_parent_proof_written) then
      ndrop_dropmixnuc_parent_proof_written = .true.
      if (masterproc) then
         write(iulog,'(A)') 'dropmixnuc parent active path = codon; column setup + activation/grow-shrink/' // &
              'oldcloud/submix/aero-tend/ccncalc direct = codon; activation core uses inlined loadaer/' // &
              'activate_modal math; native boundaries qsat/rad_constituents/pbuf/physics_ptend/outfld/allocation/endrun'
         call flush(iulog)
      end if
   end if

   deallocate( &
      nact,       &
      mact,       &
      raer,       &
      qqcw,       &
      raercol,    &
      raercol_cw, &
      coltend,    &
      coltend_cw, &
      naermod,    &
      hygro,      &
      vaerosol,   &
      fn,         &
      fm,         &
      fluxn,      &
      fluxm       )

end subroutine dropmixnuc

!===============================================================================

subroutine ndrop_dropmixnuc_aero_tend_all_codon_wrap(i, pcols_local, psetcols_local, pver_local, &
     top_lev_local, ntot_amode_local, ncnst_tot_local, nspec_max_local, slot, dtinv, gravit, &
     raer_ptrs, qqcw_ptrs, nspec_amode_local, mam_idx_local, mam_cnst_idx_local, pdel_local, &
     raercol_local, raercol_cw_local, coltend_local, coltend_cw_local, ptend_q)
   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   integer, intent(in) :: i, pcols_local, psetcols_local, pver_local, top_lev_local
   integer, intent(in) :: ntot_amode_local, ncnst_tot_local, nspec_max_local, slot
   real(r8), intent(in) :: dtinv, gravit
   type(c_ptr), target, intent(in) :: raer_ptrs(ncnst_tot_local)
   type(c_ptr), target, intent(in) :: qqcw_ptrs(ncnst_tot_local)
   integer, target, intent(in) :: nspec_amode_local(ntot_amode_local)
   integer, target, intent(in) :: mam_idx_local(ntot_amode_local,0:nspec_max_local)
   integer, target, intent(in) :: mam_cnst_idx_local(ntot_amode_local,0:nspec_max_local)
   real(r8), target, intent(in) :: pdel_local(pcols_local,pver_local)
   real(r8), target, intent(inout) :: raercol_local(pver_local,ncnst_tot_local,2)
   real(r8), target, intent(inout) :: raercol_cw_local(pver_local,ncnst_tot_local,2)
   real(r8), target, intent(inout) :: coltend_local(pcols_local,ncnst_tot_local)
   real(r8), target, intent(inout) :: coltend_cw_local(pcols_local,ncnst_tot_local)
   real(r8), target, intent(inout) :: ptend_q(psetcols_local,pver_local,pcnst)

   call ndrop_dropmixnuc_aero_tend_all_codon(int(i, c_int64_t), int(pcols_local, c_int64_t), &
        int(psetcols_local, c_int64_t), int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), &
        int(ntot_amode_local, c_int64_t), int(ncnst_tot_local, c_int64_t), int(slot, c_int64_t), &
        dtinv, gravit, c_loc(raer_ptrs(1)), c_loc(qqcw_ptrs(1)), c_loc(nspec_amode_local(1)), &
        c_loc(mam_idx_local(1,0)), c_loc(mam_cnst_idx_local(1,0)), c_loc(pdel_local(1,1)), &
        c_loc(raercol_local(1,1,1)), c_loc(raercol_cw_local(1,1,1)), c_loc(coltend_local(1,1)), &
        c_loc(coltend_cw_local(1,1)), c_loc(ptend_q(1,1,1)))
end subroutine ndrop_dropmixnuc_aero_tend_all_codon_wrap

!===============================================================================

subroutine ndrop_dropmixnuc_aero_tend_commit_ptend_codon_wrap(i, psetcols_local, pver_local, &
     top_lev_local, pcnst_local, lptr, raertend, ptend_q)
   use iso_c_binding, only: c_int64_t, c_loc

   integer, intent(in) :: i, psetcols_local, pver_local, top_lev_local, pcnst_local, lptr
   real(r8), target, intent(in) :: raertend(pver_local)
   real(r8), target, intent(inout) :: ptend_q(psetcols_local, pver_local, pcnst_local)

   call ndrop_dropmixnuc_aero_tend_commit_ptend_codon(int(i, c_int64_t), &
        int(psetcols_local, c_int64_t), int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), &
        int(pcnst_local, c_int64_t), int(lptr, c_int64_t), c_loc(raertend(1)), c_loc(ptend_q(1,1,1)))
end subroutine ndrop_dropmixnuc_aero_tend_commit_ptend_codon_wrap

!===============================================================================

subroutine explmix( q, src, ekkp, ekkm, overlapp, overlapm, &
   qold, surfrate, flxconv, pver, dt, is_unact, qactold )

   use iso_c_binding, only: c_int64_t, c_loc, c_null_ptr

   integer, intent(in) :: pver ! number of levels
   real(r8), target, intent(out) :: q(pver) ! mixing ratio to be updated
   real(r8), target, intent(in) :: qold(pver) ! mixing ratio from previous time step
   real(r8), target, intent(in) :: src(pver) ! source due to activation/nucleation (/s)
   real(r8), target, intent(in) :: ekkp(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s)
   ! below layer k  (k,k+1 interface)
   real(r8), target, intent(in) :: ekkm(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s)
   ! above layer k  (k,k+1 interface)
   real(r8), target, intent(in) :: overlapp(pver) ! cloud overlap below
   real(r8), target, intent(in) :: overlapm(pver) ! cloud overlap above
   real(r8), intent(in) :: surfrate ! surface exchange rate (/s)
   real(r8), intent(in) :: flxconv ! convergence of flux from surface
   real(r8), intent(in) :: dt ! time step (s)
   logical, intent(in) :: is_unact ! true if this is an unactivated species
   real(r8), target, intent(in),optional :: qactold(pver)
   ! mixing ratio of ACTIVATED species from previous step
   ! *** this should only be present
   !     if the current species is unactivated number/sfc/mass

   call ndrop_explmix_select_impl()

   if (use_native_ndrop_explmix_impl) then
      call explmix_native(q, src, ekkp, ekkm, overlapp, overlapm, &
           qold, surfrate, flxconv, pver, dt, is_unact, qactold)
   else if (is_unact) then
      if (.not. present(qactold)) then
         call endrun('explmix: qactold is required when is_unact is true')
      end if
      call ndrop_explmix_proof_once()
      call ndrop_explmix_codon(int(pver, c_int64_t), int(top_lev, c_int64_t), &
           surfrate, flxconv, dt, 1_c_int64_t, c_loc(q(1)), c_loc(src(1)), &
           c_loc(ekkp(1)), c_loc(ekkm(1)), c_loc(overlapp(1)), c_loc(overlapm(1)), &
           c_loc(qold(1)), c_loc(qactold(1)))
   else
      call ndrop_explmix_proof_once()
      call ndrop_explmix_codon(int(pver, c_int64_t), int(top_lev, c_int64_t), &
           surfrate, flxconv, dt, 0_c_int64_t, c_loc(q(1)), c_loc(src(1)), &
           c_loc(ekkp(1)), c_loc(ekkm(1)), c_loc(overlapp(1)), c_loc(overlapm(1)), &
           c_loc(qold(1)), c_null_ptr)
   end if

end subroutine explmix

!===============================================================================

subroutine explmix_native( q, src, ekkp, ekkm, overlapp, overlapm, &
   qold, surfrate, flxconv, pver, dt, is_unact, qactold )

   !  explicit integration of droplet/aerosol mixing
   !     with source due to activation/nucleation


   integer, intent(in) :: pver ! number of levels
   real(r8), intent(out) :: q(pver) ! mixing ratio to be updated
   real(r8), intent(in) :: qold(pver) ! mixing ratio from previous time step
   real(r8), intent(in) :: src(pver) ! source due to activation/nucleation (/s)
   real(r8), intent(in) :: ekkp(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
   ! below layer k  (k,k+1 interface)
   real(r8), intent(in) :: ekkm(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
   ! above layer k  (k,k+1 interface)
   real(r8), intent(in) :: overlapp(pver) ! cloud overlap below
   real(r8), intent(in) :: overlapm(pver) ! cloud overlap above
   real(r8), intent(in) :: surfrate ! surface exchange rate (/s)
   real(r8), intent(in) :: flxconv ! convergence of flux from surface
   real(r8), intent(in) :: dt ! time step (s)
   logical, intent(in) :: is_unact ! true if this is an unactivated species
   real(r8), intent(in),optional :: qactold(pver)
   ! mixing ratio of ACTIVATED species from previous step
   ! *** this should only be present
   !     if the current species is unactivated number/sfc/mass

   integer k,kp1,km1

   if ( is_unact ) then
      !     the qactold*(1-overlap) terms are resuspension of activated material
      do k=top_lev,pver
         kp1=min(k+1,pver)
         km1=max(k-1,top_lev)
         q(k) = qold(k) + dt*( - src(k) + ekkp(k)*(qold(kp1) - qold(k) +       &
            qactold(kp1)*(1.0_r8-overlapp(k)))               &
            + ekkm(k)*(qold(km1) - qold(k) +     &
            qactold(km1)*(1.0_r8-overlapm(k))) )
         !        force to non-negative
         !        if(q(k)<-1.e-30)then
         !           write(iulog,*)'q=',q(k),' in explmix'
         q(k)=max(q(k),0._r8)
         !        endif
      end do

      !     diffusion loss at base of lowest layer
      q(pver)=q(pver)-surfrate*qold(pver)*dt+flxconv*dt
      !        force to non-negative
      !        if(q(pver)<-1.e-30)then
      !           write(iulog,*)'q=',q(pver),' in explmix'
      q(pver)=max(q(pver),0._r8)
      !        endif
   else
      do k=top_lev,pver
         kp1=min(k+1,pver)
         km1=max(k-1,top_lev)
         q(k) = qold(k) + dt*(src(k) + ekkp(k)*(overlapp(k)*qold(kp1)-qold(k)) +      &
            ekkm(k)*(overlapm(k)*qold(km1)-qold(k)) )
         !        force to non-negative
         !        if(q(k)<-1.e-30)then
         !           write(iulog,*)'q=',q(k),' in explmix'
         q(k)=max(q(k),0._r8)
         !        endif
      end do
      !     diffusion loss at base of lowest layer
      q(pver)=q(pver)-surfrate*qold(pver)*dt+flxconv*dt
      !        force to non-negative
      !        if(q(pver)<-1.e-30)then
      !           write(iulog,*)'q=',q(pver),' in explmix'
      q(pver)=max(q(pver),0._r8)

   end if

end subroutine explmix_native

!===============================================================================

subroutine activate_modal(wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,  &
   na, nmode, volume, hygro, &
   fn, fm, fluxn, fluxm, flux_fullact )

   use iso_c_binding, only: c_int64_t, c_loc

   !      calculates number, surface, and mass fraction of aerosols activated as CCN
   !      calculates flux of cloud droplets, surface area, and aerosol mass into cloud
   !      assumes an internal mixture within each of up to nmode multiple aerosol modes
   !      a gaussiam spectrum of updrafts can be treated.

   !      mks units

   !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.


   !      input

   real(r8) :: wbar          ! grid cell mean vertical velocity (m/s)
   real(r8) :: sigw          ! subgrid standard deviation of vertical vel (m/s)
   real(r8) :: wdiab         ! diabatic vertical velocity (0 if adiabatic)
   real(r8) :: wminf         ! minimum updraft velocity for integration (m/s)
   real(r8) :: wmaxf         ! maximum updraft velocity for integration (m/s)
   real(r8) :: tair          ! air temperature (K)
   real(r8) :: rhoair        ! air density (kg/m3)
   real(r8), target :: na(:)      ! aerosol number concentration (/m3)
   integer  :: nmode      ! number of aerosol modes
   real(r8), target :: volume(:)  ! aerosol volume concentration (m3/m3)
   real(r8), target :: hygro(:)   ! hygroscopicity of aerosol mode

   !      output

   real(r8), target :: fn(:)      ! number fraction of aerosols activated
   real(r8), target :: fm(:)      ! mass fraction of aerosols activated
   real(r8), target :: fluxn(:)   ! flux of activated aerosol number fraction into cloud (cm/s)
   real(r8), target :: fluxm(:)   ! flux of activated aerosol mass fraction into cloud (cm/s)
   real(r8), target :: flux_fullact   ! flux of activated aerosol fraction assuming 100% activation (cm/s)
   !    rce-comment
   !    used for consistency check -- this should match (ekd(k)*zs(k))
   !    also, fluxm/flux_fullact gives fraction of aerosol mass flux
   !       that is activated

   !      local

   integer, parameter:: nx=200
   integer iquasisect_option, isectional
   real(r8) integ,integf
   real(r8), parameter :: p0 = 1013.25e2_r8    ! reference pressure (Pa)
   real(r8) xmin(nmode),xmax(nmode) ! ln(r) at section interfaces
   real(r8) volmin(nmode),volmax(nmode) ! volume at interfaces
   real(r8) tmass ! total aerosol mass concentration (g/cm3)
   real(r8) sign(nmode)    ! geometric standard deviation of size distribution
   real(r8) rm ! number mode radius of aerosol at max supersat (cm)
   real(r8) pres ! pressure (Pa)
   real(r8) path ! mean free path (m)
   real(r8) diff ! diffusivity (m2/s)
   real(r8) conduct ! thermal conductivity (Joule/m/sec/deg)
   real(r8) diff0,conduct0
   real(r8) es ! saturation vapor pressure
   real(r8) qs ! water vapor saturation mixing ratio
   real(r8) dqsdt ! change in qs with temperature
   real(r8) dqsdp ! change in qs with pressure
   real(r8) g ! thermodynamic function (m2/s)
   real(r8), target :: zeta(nmode), eta(nmode)
   real(r8) lnsmax ! ln(smax)
   real(r8) alpha
   real(r8) gamma
   real(r8) beta
   real(r8), target :: sqrtg(nmode)
   real(r8), target :: amcube(nmode) ! cube of dry mode radius (m)
   real(r8) :: smcrit(nmode) ! critical supersatuation for activation
   real(r8), target :: lnsm(nmode) ! ln(smcrit)
   real(r8), target :: smc(nmode) ! critical supersaturation for number mode radius
   real(r8) sumflx_fullact
   real(r8), target :: sumflxn(nmode)
   real(r8), target :: sumflxm(nmode)
   real(r8), target :: sumfn(nmode)
   real(r8), target :: sumfm(nmode)
   real(r8), target :: fnold(nmode)   ! number fraction activated
   real(r8), target :: fmold(nmode)   ! mass fraction activated
   real(r8) wold,gold
   real(r8) alogam
   real(r8) rlo,rhi,xint1,xint2,xint3,xint4
   real(r8) wmin,wmax,w,dw,dwmax,dwmin,wnuc,dwnew,wb
   real(r8) dfmin,dfmax,fnew,fold,fnmin,fnbar,fsbar,fmbar
   real(r8) alw,sqrtalw
   real(r8) smax
   real(r8) x,arg
   real(r8) xmincoeff,xcut,volcut,surfcut
   real(r8) z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
   real(r8) etafactor1,etafactor2max
   real(r8), target :: etafactor2(nmode)
   integer m,n
   integer(c_int64_t) :: codon_status
   !      numerical integration parameters
   real(r8), parameter :: eps=0.3_r8,fmax=0.99_r8,sds=3._r8

   real(r8), parameter :: namin=1.e6_r8   ! minimum aerosol number concentration (/m3)

   integer ndist(nx)  ! accumulates frequency distribution of integration bins required
   data ndist/nx*0/
   save ndist

   call ndrop_activate_modal_select_impl()
   if (.not. use_native_ndrop_activate_modal_impl) then
      pres=rair*rhoair*tair
      call qsat(tair, pres, es, qs)
      call ndrop_activate_modal_proof_once()
      codon_status = ndrop_activate_modal_core_codon(wbar, sigw, wdiab, wminf, wmaxf, &
           tair, rhoair, qs, int(nmode, c_int64_t), rair, p0, t0, rhoh2o, latvap, &
           cpair, rh2o, gravit, pi, aten, twothird, sq2, sqpi, sixth, zero, &
           c_loc(na(1)), c_loc(volume(1)), c_loc(hygro(1)), c_loc(alogsig(1)), &
           c_loc(exp45logsig(1)), c_loc(f1(1)), c_loc(f2(1)), c_loc(fn(1)), &
           c_loc(fm(1)), c_loc(fluxn(1)), c_loc(fluxm(1)), c_loc(flux_fullact), &
           c_loc(zeta(1)), c_loc(eta(1)), c_loc(etafactor2(1)), c_loc(sqrtg(1)), &
           c_loc(amcube(1)), c_loc(smc(1)), c_loc(lnsm(1)), c_loc(sumflxn(1)), &
           c_loc(sumflxm(1)), c_loc(sumfn(1)), c_loc(sumfm(1)), c_loc(fnold(1)), &
           c_loc(fmold(1)))
      if (codon_status == 0_c_int64_t) return
      if (codon_status == 1_c_int64_t) then
         call endrun('activate: Codon activation integration loop did not converge')
      else if (codon_status == 2_c_int64_t) then
         call endrun('activate: Codon activation fraction exceeded one')
      else
         call endrun('activate: Codon activation returned unknown status')
      end if
   end if

   fn(:)=0._r8
   fm(:)=0._r8
   fluxn(:)=0._r8
   fluxm(:)=0._r8
   flux_fullact=0._r8

   if(nmode.eq.1.and.na(1).lt.1.e-20_r8)return

   if(sigw.le.1.e-5_r8.and.wbar.le.0._r8)return

   pres=rair*rhoair*tair
   diff0=0.211e-4_r8*(p0/pres)*(tair/t0)**1.94_r8
   conduct0=(5.69_r8+0.017_r8*(tair-t0))*4.186e2_r8*1.e-5_r8 ! convert to J/m/s/deg
   call qsat(tair, pres, es, qs)
   dqsdt=latvap/(rh2o*tair*tair)*qs
   alpha=gravit*(latvap/(cpair*rh2o*tair*tair)-1._r8/(rair*tair))
   gamma=(1+latvap/cpair*dqsdt)/(rhoair*qs)
   etafactor2max=1.e10_r8/(alpha*wmaxf)**1.5_r8 ! this should make eta big if na is very small.

   do m=1,nmode
      if(volume(m).gt.1.e-39_r8.and.na(m).gt.1.e-39_r8)then
         !            number mode radius (m)
         !           write(iulog,*)'alogsig,volc,na=',alogsig(m),volc(m),na(m)
         amcube(m)=(3._r8*volume(m)/(4._r8*pi*exp45logsig(m)*na(m)))  ! only if variable size dist
         !           growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
         !           should depend on mean radius of mode to account for gas kinetic effects
         !           see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
         !           for approriate size to use for effective diffusivity.
         g=1._r8/(rhoh2o/(diff0*rhoair*qs)                                    &
            +latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair)-1._r8))
         sqrtg(m)=sqrt(g)
         beta=2._r8*pi*rhoh2o*g*gamma
         etafactor2(m)=1._r8/(na(m)*beta*sqrtg(m))
         if(hygro(m).gt.1.e-10_r8)then
            smc(m)=2._r8*aten*sqrt(aten/(27._r8*hygro(m)*amcube(m))) ! only if variable size dist
         else
            smc(m)=100._r8
         endif
         !	    write(iulog,*)'sm,hygro,amcube=',smcrit(m),hygro(m),amcube(m)
      else
         g=1._r8/(rhoh2o/(diff0*rhoair*qs)                                    &
            +latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair)-1._r8))
         sqrtg(m)=sqrt(g)
         smc(m)=1._r8
         etafactor2(m)=etafactor2max ! this should make eta big if na is very small.
      endif
      lnsm(m)=log(smc(m)) ! only if variable size dist
      !	 write(iulog,'(a,i4,4g12.2)')'m,na,amcube,hygro,sm,lnsm=', &
      !                   m,na(m),amcube(m),hygro(m),sm(m),lnsm(m)
   enddo

   if(sigw.gt.1.e-5_r8)then ! spectrum of updrafts

      wmax=min(wmaxf,wbar+sds*sigw)
      wmin=max(wminf,-wdiab)
      wmin=max(wmin,wbar-sds*sigw)
      w=wmin
      dwmax=eps*sigw
      dw=dwmax
      dfmax=0.2_r8
      dfmin=0.1_r8
      if(wmax.le.w)then
         do m=1,nmode
            fluxn(m)=0._r8
            fn(m)=0._r8
            fluxm(m)=0._r8
            fm(m)=0._r8
         enddo
         flux_fullact=0._r8
         return
      endif
      do m=1,nmode
         sumflxn(m)=0._r8
         sumfn(m)=0._r8
         fnold(m)=0._r8
         sumflxm(m)=0._r8
         sumfm(m)=0._r8
         fmold(m)=0._r8
      enddo
      sumflx_fullact=0._r8

      fold=0._r8
      wold=0._r8
      gold=0._r8

      dwmin = min( dwmax, 0.01_r8 )

      do n=1,200
100      wnuc=w+wdiab
         !           write(iulog,*)'wnuc=',wnuc
         alw=alpha*wnuc
         sqrtalw=sqrt(alw)
         etafactor1=alw*sqrtalw

         do m=1,nmode
            eta(m)=etafactor1*etafactor2(m)
            zeta(m)=twothird*sqrtalw*aten/sqrtg(m)
         enddo

         call maxsat(zeta,eta,nmode,smc,smax)
         !	      write(iulog,*)'w,smax=',w,smax

         lnsmax=log(smax)

         x=twothird*(lnsm(nmode)-lnsmax)/(sq2*alogsig(nmode))
         fnew=0.5_r8*(1._r8-erf(x))


         dwnew = dw
         if(fnew-fold.gt.dfmax.and.n.gt.1)then
            !              reduce updraft increment for greater accuracy in integration
            if (dw .gt. 1.01_r8*dwmin) then
               dw=0.7_r8*dw
               dw=max(dw,dwmin)
               w=wold+dw
               go to 100
            else
               dwnew = dwmin
            endif
         endif

         if(fnew-fold.lt.dfmin)then
            !              increase updraft increment to accelerate integration
            dwnew=min(1.5_r8*dw,dwmax)
         endif
         fold=fnew

         z=(w-wbar)/(sigw*sq2)
         g=exp(-z*z)
         fnmin=1._r8
         xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3

         do m=1,nmode
            !              modal
            x=twothird*(lnsm(m)-lnsmax)/(sq2*alogsig(m))
            fn(m)=0.5_r8*(1._r8-erf(x))
            fnmin=min(fn(m),fnmin)
            !               integration is second order accurate
            !               assumes linear variation of f*g with w
            fnbar=(fn(m)*g+fnold(m)*gold)
            arg=x-1.5_r8*sq2*alogsig(m)
            fm(m)=0.5_r8*(1._r8-erf(arg))
            fmbar=(fm(m)*g+fmold(m)*gold)
            wb=(w+wold)
            if(w.gt.0._r8)then
               sumflxn(m)=sumflxn(m)+sixth*(wb*fnbar           &
                  +(fn(m)*g*w+fnold(m)*gold*wold))*dw
               sumflxm(m)=sumflxm(m)+sixth*(wb*fmbar           &
                  +(fm(m)*g*w+fmold(m)*gold*wold))*dw
            endif
            sumfn(m)=sumfn(m)+0.5_r8*fnbar*dw
            !	       write(iulog,'(a,9g10.2)')'lnsmax,lnsm(m),x,fn(m),fnold(m),g,gold,fnbar,dw=',lnsmax,lnsm(m),x,fn(m),fnold(m),g,gold,fnbar,dw
            fnold(m)=fn(m)
            sumfm(m)=sumfm(m)+0.5_r8*fmbar*dw
            fmold(m)=fm(m)
         enddo
         !           same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact &
            + sixth*(wb*(g+gold) + (g*w+gold*wold))*dw
         !            sumg=sumg+0.5_r8*(g+gold)*dw
         gold=g
         wold=w
         dw=dwnew
         if(n.gt.1.and.(w.gt.wmax.or.fnmin.gt.fmax))go to 20
         w=w+dw
      enddo
      write(iulog,*)'do loop is too short in activate'
      write(iulog,*)'wmin=',wmin,' w=',w,' wmax=',wmax,' dw=',dw
      write(iulog,*)'wbar=',wbar,' sigw=',sigw,' wdiab=',wdiab
      write(iulog,*)'wnuc=',wnuc
      write(iulog,*)'na=',(na(m),m=1,nmode)
      write(iulog,*)'fn=',(fn(m),m=1,nmode)
      !   dump all subr parameters to allow testing with standalone code
      !   (build a driver that will read input and call activate)
      write(iulog,*)'wbar,sigw,wdiab,tair,rhoair,nmode='
      write(iulog,*) wbar,sigw,wdiab,tair,rhoair,nmode
      write(iulog,*)'na=',na
      write(iulog,*)'volume=', (volume(m),m=1,nmode)
      write(iulog,*)'hydro='
      write(iulog,*) hygro

      call endrun
20    continue
      ndist(n)=ndist(n)+1
      if(w.lt.wmaxf)then

         !            contribution from all updrafts stronger than wmax
         !            assuming constant f (close to fmax)
         wnuc=w+wdiab

         z1=(w-wbar)/(sigw*sq2)
         z2=(wmaxf-wbar)/(sigw*sq2)
         g=exp(-z1*z1)
         integ=sigw*0.5_r8*sq2*sqpi*(erf(z2)-erf(z1))
         !            consider only upward flow into cloud base when estimating flux
         wf1=max(w,zero)
         zf1=(wf1-wbar)/(sigw*sq2)
         gf1=exp(-zf1*zf1)
         wf2=max(wmaxf,zero)
         zf2=(wf2-wbar)/(sigw*sq2)
         gf2=exp(-zf2*zf2)
         gf=(gf1-gf2)
         integf=wbar*sigw*0.5_r8*sq2*sqpi*(erf(zf2)-erf(zf1))+sigw*sigw*gf

         do m=1,nmode
            sumflxn(m)=sumflxn(m)+integf*fn(m)
            sumfn(m)=sumfn(m)+fn(m)*integ
            sumflxm(m)=sumflxm(m)+integf*fm(m)
            sumfm(m)=sumfm(m)+fm(m)*integ
         enddo
         !           same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact + integf
         !            sumg=sumg+integ
      endif


      do m=1,nmode
         fn(m)=sumfn(m)/(sq2*sqpi*sigw)
         !            fn(m)=sumfn(m)/(sumg)
         if(fn(m).gt.1.01_r8)then
            write(iulog,*)'fn=',fn(m),' > 1 in activate'
            write(iulog,*)'w,m,na,amcube=',w,m,na(m),amcube(m)
            write(iulog,*)'integ,sumfn,sigw=',integ,sumfn(m),sigw
            call endrun('activate')
         endif
         fluxn(m)=sumflxn(m)/(sq2*sqpi*sigw)
         fm(m)=sumfm(m)/(sq2*sqpi*sigw)
         !            fm(m)=sumfm(m)/(sumg)
         if(fm(m).gt.1.01_r8)then
            write(iulog,*)'fm=',fm(m),' > 1 in activate'
         endif
         fluxm(m)=sumflxm(m)/(sq2*sqpi*sigw)
      enddo
      !        same form as fluxm
      flux_fullact = sumflx_fullact/(sq2*sqpi*sigw)

   else

      !        single updraft
      wnuc=wbar+wdiab

      if(wnuc.gt.0._r8)then

         w=wbar
         alw=alpha*wnuc
         sqrtalw=sqrt(alw)
         etafactor1=alw*sqrtalw

         do m=1,nmode
            eta(m)=etafactor1*etafactor2(m)
            zeta(m)=twothird*sqrtalw*aten/sqrtg(m)
         enddo

         call maxsat(zeta,eta,nmode,smc,smax)

         lnsmax=log(smax)
         xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3


         do m=1,nmode
            !                 modal
            x=twothird*(lnsm(m)-lnsmax)/(sq2*alogsig(m))
            fn(m)=0.5_r8*(1._r8-erf(x))
            arg=x-1.5_r8*sq2*alogsig(m)
            fm(m)=0.5_r8*(1._r8-erf(arg))
            if(wbar.gt.0._r8)then
               fluxn(m)=fn(m)*w
               fluxm(m)=fm(m)*w
            endif
         enddo
         flux_fullact = w
      endif

   endif

end subroutine activate_modal

!===============================================================================

subroutine maxsat(zeta,eta,nmode,smc,smax)

   use iso_c_binding, only: c_int64_t, c_loc

   integer  :: nmode
   real(r8), target :: smc(nmode)
   real(r8), target :: zeta(nmode)
   real(r8), target :: eta(nmode)
   real(r8) :: smax

   call ndrop_maxsat_select_impl()
   if (use_native_ndrop_maxsat_impl) then
      call maxsat_native(zeta, eta, nmode, smc, smax)
   else
      call ndrop_maxsat_proof_once()
      smax = ndrop_maxsat_codon(int(nmode, c_int64_t), c_loc(zeta(1)), c_loc(eta(1)), &
           c_loc(smc(1)), c_loc(f1(1)), c_loc(f2(1)))
   end if

end subroutine maxsat

!===============================================================================

subroutine maxsat_native(zeta,eta,nmode,smc,smax)

   !      calculates maximum supersaturation for multiple
   !      competing aerosol modes.

   !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

   integer  :: nmode ! number of modes
   real(r8) :: smc(nmode) ! critical supersaturation for number mode radius
   real(r8) :: zeta(nmode)
   real(r8) :: eta(nmode)
   real(r8) :: smax ! maximum supersaturation
   integer  :: m  ! mode index
   real(r8) :: sum, g1, g2, g1sqrt, g2sqrt

   do m=1,nmode
      if(zeta(m).gt.1.e5_r8*eta(m).or.smc(m)*smc(m).gt.1.e5_r8*eta(m))then
         !            weak forcing. essentially none activated
         smax=1.e-20_r8
      else
         !            significant activation of this mode. calc activation all modes.
         go to 1
      endif
   enddo

   return

1  continue

   sum=0
   do m=1,nmode
      if(eta(m).gt.1.e-20_r8)then
         g1=zeta(m)/eta(m)
         g1sqrt=sqrt(g1)
         g1=g1sqrt*g1
         g2=smc(m)/sqrt(eta(m)+3._r8*zeta(m))
         g2sqrt=sqrt(g2)
         g2=g2sqrt*g2
         sum=sum+(f1(m)*g1+f2(m)*g2)/(smc(m)*smc(m))
      else
         sum=1.e20_r8
      endif
   enddo

   smax=1._r8/sqrt(sum)

end subroutine maxsat_native

!===============================================================================

subroutine ccncalc(state, pbuf, cs, ccn)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   ! calculates number concentration of aerosols activated as CCN at
   ! supersaturation supersat.
   ! assumes an internal mixture of a multiple externally-mixed aerosol modes
   ! cgs units

   ! Ghan et al., Atmos. Res., 1993, 198-221.

   ! arguments

   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)


   real(r8), target, intent(in)  :: cs(pcols,pver)       ! air density (kg/m3)
   real(r8), target, intent(out) :: ccn(pcols,pver,psat) ! number conc of aerosols activated at supersat (#/m3)

   ! local

   integer :: lchnk ! chunk index
   integer :: ncol  ! number of columns
   real(r8), pointer :: tair(:,:)     ! air temperature (K)
   real(r8), pointer :: specmmr(:,:)

   real(r8), target :: naerosol(pcols) ! interstit+activated aerosol number conc (/m3)
   real(r8), target :: vaerosol(pcols) ! interstit+activated aerosol volume conc (m3/m3)

   real(r8), target :: amcube(pcols)
   real(r8), target :: super(psat) ! supersaturation
   real(r8), allocatable, target :: amcubecoef(:)
   real(r8), allocatable, target :: argfactor(:)
   real(r8) :: surften       ! surface tension of water w/respect to air (N/m)
   real(r8) surften_coef
   real(r8), target :: a(pcols) ! surface tension parameter
   real(r8), target :: hygro(pcols)  ! aerosol hygroscopicity
   real(r8), target :: sm(pcols)  ! critical supersaturation at mode radius
   real(r8), target :: arg(pcols)
   real(r8), pointer :: qqcw(:,:)
   !     mathematical constants
   real(r8) twothird,sq2
   integer l,m,n,i,k
   real(r8) log,cc
   real(r8) smcoefcoef
   real(r8), target :: smcoef(pcols)
   integer phase ! phase of aerosol
   type(c_ptr), target :: ccn_species_raer_ptrs(ncnst_tot,ntot_amode)
   type(c_ptr), target :: ccn_species_qqcw_ptrs(ncnst_tot,ntot_amode)
   type(c_ptr), target :: ccn_num_raer_ptrs(ntot_amode)
   type(c_ptr), target :: ccn_num_qqcw_ptrs(ntot_amode)
   real(r8), target :: ccn_species_specdens(ncnst_tot,ntot_amode)
   real(r8), target :: ccn_species_spechygro(ncnst_tot,ntot_amode)
   real(r8), target :: ccn_voltonumblo(ntot_amode)
   real(r8), target :: ccn_voltonumbhi(ntot_amode)
   !-------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   tair  => state%t

   allocate( &
      amcubecoef(ntot_amode), &
      argfactor(ntot_amode)   )

   call ndrop_ccncalc_helpers_select_impl()

   super(:)=supersat(:)*0.01_r8
   sq2=sqrt(2._r8)
   twothird=2._r8/3._r8
   surften=0.076_r8
   surften_coef=2._r8*mwh2o*surften/(r_universal*rhoh2o)
   smcoefcoef=2._r8/sqrt(27._r8)

   do m=1,ntot_amode
      amcubecoef(m)=3._r8/(4._r8*pi*exp45logsig(m))
      argfactor(m)=twothird/(sq2*alogsig(m))
   end do

   if (.not. use_native_ndrop_ccncalc_helpers_impl) then
      do m=1,ntot_amode
         do l=1,nspec_amode(m)
            call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, specmmr)
            call rad_cnst_get_aer_mmr(0, m, l, 'c', state, pbuf, qqcw)
            call rad_cnst_get_aer_props(0, m, l, density_aer=ccn_species_specdens(l,m), &
                 hygro_aer=ccn_species_spechygro(l,m))
            ccn_species_raer_ptrs(l,m) = c_loc(specmmr(1,1))
            ccn_species_qqcw_ptrs(l,m) = c_loc(qqcw(1,1))
         end do
         call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, specmmr)
         call rad_cnst_get_mode_num(0, m, 'c', state, pbuf, qqcw)
         ccn_num_raer_ptrs(m) = c_loc(specmmr(1,1))
         ccn_num_qqcw_ptrs(m) = c_loc(qqcw(1,1))
         ccn_voltonumblo(m) = voltonumblo_amode(m)
         ccn_voltonumbhi(m) = voltonumbhi_amode(m)
      end do

      call ndrop_ccncalc_helpers_proof_once()
      call ndrop_ccncalc_direct_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(top_lev, c_int64_t), int(psat, c_int64_t), &
           int(ntot_amode, c_int64_t), int(ncnst_tot, c_int64_t), pi, surften_coef, &
           smcoefcoef, c_loc(nspec_amode(1)), c_loc(ccn_species_raer_ptrs(1,1)), &
           c_loc(ccn_species_qqcw_ptrs(1,1)), c_loc(ccn_species_specdens(1,1)), &
           c_loc(ccn_species_spechygro(1,1)), c_loc(ccn_num_raer_ptrs(1)), &
           c_loc(ccn_num_qqcw_ptrs(1)), c_loc(ccn_voltonumblo(1)), c_loc(ccn_voltonumbhi(1)), &
           c_loc(alogsig(1)), c_loc(exp45logsig(1)), c_loc(tair(1,1)), c_loc(cs(1,1)), &
           c_loc(super(1)), c_loc(naerosol(1)), c_loc(vaerosol(1)), c_loc(hygro(1)), &
           c_loc(amcube(1)), c_loc(smcoef(1)), c_loc(sm(1)), c_loc(arg(1)), c_loc(ccn(1,1,1)))

      deallocate( &
         amcubecoef, &
         argfactor   )
      return
   end if

   if (use_native_ndrop_ccncalc_helpers_impl) then
      ccn = 0._r8
   else
      call ndrop_ccncalc_helpers_proof_once()
      call ndrop_ccncalc_zero_codon(int(pcols, c_int64_t), int(pver, c_int64_t), &
           int(psat, c_int64_t), c_loc(ccn(1,1,1)))
   end if

   do k=top_lev,pver

      if (use_native_ndrop_ccncalc_helpers_impl) then
         do i=1,ncol
            a(i)=surften_coef/tair(i,k)
            smcoef(i)=smcoefcoef*a(i)*sqrt(a(i))
         end do
      else
         call ndrop_ccncalc_helpers_proof_once()
         call ndrop_ccncalc_level_coeffs_codon(int(ncol, c_int64_t), int(k, c_int64_t), &
              int(pcols, c_int64_t), surften_coef, smcoefcoef, c_loc(tair(1,1)), c_loc(smcoef(1)))
      end if

      do m=1,ntot_amode

         phase=3 ! interstitial+cloudborne

         call loadaer( &
            state, pbuf, 1, ncol, k, &
            m, cs, phase, naerosol, vaerosol, &
            hygro)

         if (use_native_ndrop_ccncalc_helpers_impl) then
            where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol)) ! critical supersaturation
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
            endwhere
            do l=1,psat
               do i=1,ncol
                  arg(i)=argfactor(m)*log(sm(i)/super(l))
                  ccn(i,k,l)=ccn(i,k,l)+naerosol(i)*0.5_r8*(1._r8-erf(arg(i)))
               enddo
            enddo
         else
            call ndrop_ccncalc_helpers_proof_once()
            call ndrop_ccncalc_mode_accum_codon(int(ncol, c_int64_t), int(k, c_int64_t), &
                 int(pcols, c_int64_t), int(pver, c_int64_t), int(psat, c_int64_t), &
                 amcubecoef(m), argfactor(m), c_loc(naerosol(1)), c_loc(vaerosol(1)), &
                 c_loc(hygro(1)), c_loc(smcoef(1)), c_loc(super(1)), c_loc(amcube(1)), &
                 c_loc(sm(1)), c_loc(arg(1)), c_loc(ccn(1,1,1)))
         end if
      enddo
   enddo

   if (use_native_ndrop_ccncalc_helpers_impl) then
      ccn(:ncol,:,:)=ccn(:ncol,:,:)*1.e-6_r8 ! convert from #/m3 to #/cm3
   else
      call ndrop_ccncalc_helpers_proof_once()
      call ndrop_ccncalc_scale_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(psat, c_int64_t), c_loc(ccn(1,1,1)))
   end if

   deallocate( &
      amcubecoef, &
      argfactor   )

end subroutine ccncalc

!===============================================================================

subroutine loadaer( &
   state, pbuf, istart, istop, k, &
   m, cs, phase, naerosol, &
   vaerosol, hygro)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   ! return aerosol number, volume concentrations, and bulk hygroscopicity

   ! input arguments
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)

   integer,  intent(in) :: istart      ! start column index (1 <= istart <= istop <= pcols)
   integer,  intent(in) :: istop       ! stop column index  
   integer,  intent(in) :: m           ! mode index
   integer,  intent(in) :: k           ! level index
   real(r8), target, intent(in) :: cs(:,:)     ! air density (kg/m3)
   integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

   ! output arguments
   real(r8), target, intent(out) :: naerosol(:)  ! number conc (1/m3)
   real(r8), target, intent(out) :: vaerosol(:)  ! volume conc (m3/m3)
   real(r8), target, intent(out) :: hygro(:)     ! bulk hygroscopicity of mode

   ! internal
   integer  :: lchnk               ! chunk identifier

   real(r8), pointer :: raer(:,:) ! interstitial aerosol mass, number mixing ratios
   real(r8), pointer :: qqcw(:,:) ! cloud-borne aerosol mass, number mixing ratios
   real(r8) :: specdens, spechygro

   real(r8) :: vol(pcols) ! aerosol volume mixing ratio
   type(c_ptr), target :: species_raer_ptrs(ncnst_tot)
   type(c_ptr), target :: species_qqcw_ptrs(ncnst_tot)
   real(r8), target :: species_specdens(ncnst_tot)
   real(r8), target :: species_spechygro(ncnst_tot)
   integer  :: i, l
   !-------------------------------------------------------------------------------

   lchnk = state%lchnk

   call ndrop_loadaer_helpers_select_impl()
   if (use_native_ndrop_loadaer_helpers_impl) then
      do i = istart, istop
         vaerosol(i) = 0._r8
         hygro(i)    = 0._r8
      end do
   else
      call ndrop_loadaer_helpers_proof_once()
      if (phase < 1 .or. phase > 3) then
         write(iulog,*)'phase=',phase,' in loadaer'
         call endrun('phase error in loadaer')
      end if
   end if

   if (use_native_ndrop_loadaer_helpers_impl) then
      do l = 1, nspec_amode(m)

         call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, raer)
         call rad_cnst_get_aer_mmr(0, m, l, 'c', state, pbuf, qqcw)
         call rad_cnst_get_aer_props(0, m, l, density_aer=specdens, hygro_aer=spechygro)

         if (phase == 3) then
            do i = istart, istop
               vol(i) = max(raer(i,k) + qqcw(i,k), 0._r8)/specdens
            end do
         else if (phase == 2) then
            do i = istart, istop
               vol(i) = max(qqcw(i,k), 0._r8)/specdens
            end do
         else if (phase == 1) then
            do i = istart, istop
               vol(i) = max(raer(i,k), 0._r8)/specdens
            end do
         else
            write(iulog,*)'phase=',phase,' in loadaer'
            call endrun('phase error in loadaer')
         end if

         do i = istart, istop
            vaerosol(i) = vaerosol(i) + vol(i)
            hygro(i)    = hygro(i) + vol(i)*spechygro
         end do

      end do
   else
      do l = 1, nspec_amode(m)
         call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, raer)
         call rad_cnst_get_aer_mmr(0, m, l, 'c', state, pbuf, qqcw)
         call rad_cnst_get_aer_props(0, m, l, density_aer=species_specdens(l), &
              hygro_aer=species_spechygro(l))
         species_raer_ptrs(l) = c_loc(raer(1,1))
         species_qqcw_ptrs(l) = c_loc(qqcw(1,1))
      end do
   end if

   if (use_native_ndrop_loadaer_helpers_impl) then
      do i = istart, istop
         if (vaerosol(i) > 1.0e-30_r8) then   ! +++xl add 8/2/2007
            hygro(i)    = hygro(i)/(vaerosol(i))
            vaerosol(i) = vaerosol(i)*cs(i,k)
         else
            hygro(i)    = 0.0_r8
            vaerosol(i) = 0.0_r8
         end if
      end do
   end if

   ! aerosol number
   call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, raer)
   call rad_cnst_get_mode_num(0, m, 'c', state, pbuf, qqcw)

   if (.not. use_native_ndrop_loadaer_helpers_impl) then
      call ndrop_loadaer_helpers_proof_once()
      call ndrop_loadaer_direct_codon(int(istart, c_int64_t), int(istop, c_int64_t), &
           int(k, c_int64_t), int(pcols, c_int64_t), int(nspec_amode(m), c_int64_t), &
           int(phase, c_int64_t), voltonumblo_amode(m), voltonumbhi_amode(m), &
           c_loc(species_raer_ptrs(1)), c_loc(species_qqcw_ptrs(1)), c_loc(species_specdens(1)), &
           c_loc(species_spechygro(1)), c_loc(raer(1,1)), c_loc(qqcw(1,1)), c_loc(cs(1,1)), &
           c_loc(vaerosol(1)), c_loc(hygro(1)), c_loc(naerosol(1)))
   else if (phase == 3) then
      do i = istart, istop
         naerosol(i) = (raer(i,k) + qqcw(i,k))*cs(i,k)
      end do
   else if (phase == 2) then
      do i = istart, istop
         naerosol(i) = qqcw(i,k)*cs(i,k)
      end do
   else
      do i = istart, istop
         naerosol(i) = raer(i,k)*cs(i,k)
      end do
   end if
   ! adjust number so that dgnumlo < dgnum < dgnumhi
   if (use_native_ndrop_loadaer_helpers_impl) then
      do i = istart, istop
         naerosol(i) = max(naerosol(i), vaerosol(i)*voltonumbhi_amode(m))
         naerosol(i) = min(naerosol(i), vaerosol(i)*voltonumblo_amode(m))
      end do
   end if

end subroutine loadaer

!===============================================================================

end module ndrop
