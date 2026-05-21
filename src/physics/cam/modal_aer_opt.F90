module modal_aer_opt

! parameterizes aerosol coefficients using chebychev polynomial
! parameterize aerosol radiative properties in terms of
! surface mode wet radius and wet refractive index

! Ghan and Zaveri, JGR 2007.

! uses Wiscombe's (1979) mie scattering code


use shr_kind_mod,      only: r8 => shr_kind_r8, shr_kind_cl
use ppgrid,            only: pcols, pver, pverp
use constituents,      only: pcnst
use spmd_utils,        only: masterproc
use ref_pres,          only: top_lev => clim_modal_aero_top_lev
use physconst,         only: rhoh2o, rga, rair
use radconstants,      only: nswbands, nlwbands, idx_sw_diag, idx_uv_diag, idx_nir_diag
use rad_constituents,  only: n_diag, rad_cnst_get_call_list, rad_cnst_get_info, rad_cnst_get_aer_mmr, &
                             rad_cnst_get_aer_props, rad_cnst_get_mode_props
use physics_types,     only: physics_state

use physics_buffer, only : pbuf_get_index,physics_buffer_desc, pbuf_get_field
use pio,               only: file_desc_t, var_desc_t, pio_inq_dimlen, pio_inq_dimid, pio_inq_varid, &
                             pio_get_var, pio_nowrite, pio_closefile
use cam_pio_utils,     only: cam_pio_openfile
use cam_history,       only: phys_decomp, addfld, add_default, outfld
use cam_history_support, only: fillvalue
use cam_logfile,       only: iulog
use perf_mod,          only: t_startf, t_stopf
use cam_abortutils,    only: endrun

use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
use modal_aero_calcsize,    only: modal_aero_calcsize_diag

implicit none
private
save

public :: modal_aer_opt_readnl, modal_aer_opt_init, modal_aero_sw, modal_aero_lw


character(len=*), parameter :: unset_str = 'UNSET'

! Namelist variables:
character(shr_kind_cl)      :: modal_optics_file = unset_str   ! full pathname for modal optics dataset
character(shr_kind_cl)      :: water_refindex_file = unset_str ! full pathname for water refractive index dataset

! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
! in terms of refractive index and wet radius
integer, parameter :: ncoef=5, prefr=7, prefi=10

real(r8) :: xrmin, xrmax

! refractive index for water read in read_water_refindex
complex(r8) :: crefwsw(nswbands) ! complex refractive index for water visible
complex(r8) :: crefwlw(nlwbands) ! complex refractive index for water infrared

! physics buffer indices
integer :: dgnumwet_idx = -1
integer :: qaerwat_idx  = -1

character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                       '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

logical :: use_native_modal_aer_opt_helpers_impl = .false.
logical :: modal_aer_opt_helpers_impl_selected = .false.
logical :: modal_aer_opt_helpers_proof_written = .false.
logical :: modal_aer_opt_lw_helpers_proof_written = .false.
logical :: modal_aer_opt_sw_guard_helpers_proof_written = .false.
logical :: modal_aer_opt_sw_water_refr_proof_written = .false.
logical :: modal_aer_opt_sw_optics_tau_proof_written = .false.

interface
   subroutine modal_aer_opt_size_parameters_codon(pcols_c, pver_c, top_lev_c, ncol_c, ncoef_c, &
        sigma_logr_aer_c, xrmin_c, xrmax_c, dgnumwet_p, radsurf_p, logradsurf_p, cheb_p) &
        bind(c, name="modal_aer_opt_size_parameters_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: pcols_c, pver_c, top_lev_c, ncol_c, ncoef_c
      real(c_double), value :: sigma_logr_aer_c, xrmin_c, xrmax_c
      type(c_ptr), value :: dgnumwet_p, radsurf_p, logradsurf_p, cheb_p
   end subroutine modal_aer_opt_size_parameters_codon

   subroutine modal_aer_opt_lw_size_parameters_codon(pcols_c, pver_c, top_lev_c, ncol_c, ncoef_c, &
        sigma_logr_aer_c, xrmin_c, xrmax_c, dgnumwet_p, cheby_p) &
        bind(c, name="modal_aer_opt_lw_size_parameters_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: pcols_c, pver_c, top_lev_c, ncol_c, ncoef_c
      real(c_double), value :: sigma_logr_aer_c, xrmin_c, xrmax_c
      type(c_ptr), value :: dgnumwet_p, cheby_p
   end subroutine modal_aer_opt_lw_size_parameters_codon

   subroutine modal_aer_opt_binterp_codon(pcols_c, ncol_c, km_c, im_c, jm_c, table_p, x_p, y_p, &
        xtab_p, ytab_p, ix_p, jy_p, t_p, u_p, out_p) bind(c, name="modal_aer_opt_binterp_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: pcols_c, ncol_c, km_c, im_c, jm_c
      type(c_ptr), value :: table_p, x_p, y_p, xtab_p, ytab_p, ix_p, jy_p, t_p, u_p, out_p
   end subroutine modal_aer_opt_binterp_codon

   subroutine modal_aer_opt_sw_binterp3_codon(pcols_c, ncol_c, ncoef_c, prefr_c, prefi_c, &
        extpsw_p, abspsw_p, asmpsw_p, refr_p, refi_p, refrtabsw_p, refitabsw_p, itab_p, &
        jtab_p, ttab_p, utab_p, cext_p, cabs_p, casm_p) bind(c, name="modal_aer_opt_sw_binterp3_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: pcols_c, ncol_c, ncoef_c, prefr_c, prefi_c
      type(c_ptr), value :: extpsw_p, abspsw_p, asmpsw_p, refr_p, refi_p
      type(c_ptr), value :: refrtabsw_p, refitabsw_p, itab_p, jtab_p, ttab_p, utab_p
      type(c_ptr), value :: cext_p, cabs_p, casm_p
   end subroutine modal_aer_opt_sw_binterp3_codon

   subroutine modal_aer_opt_sw_init_state_codon(ncol_c, pcols_c, pver_c, nswbands_c, rga_c, rair_c, &
        pdeldry_p, pmid_p, state_t_p, tauxar_p, wa_p, ga_p, fa_p, mass_p, air_density_p) &
        bind(c, name="modal_aer_opt_sw_init_state_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, nswbands_c
      real(c_double), value :: rga_c, rair_c
      type(c_ptr), value :: pdeldry_p, pmid_p, state_t_p, tauxar_p, wa_p, ga_p, fa_p
      type(c_ptr), value :: mass_p, air_density_p
   end subroutine modal_aer_opt_sw_init_state_codon

   subroutine modal_aer_opt_sw_zero_diagnostics_codon(ncol_c, pcols_c, pver_c, extinct_p, absorb_p, &
        extinctuv_p, extinctnir_p, aodvis_p, aodvisst_p, aodabs_p, aodabsbc_p, ssavis_p, &
        burdendust_p, burdenso4_p, burdenpom_p, burdensoa_p, burdenbc_p, burdenseasalt_p, &
        dustaod_p, so4aod_p, pomaod_p, soaaod_p, bcaod_p, seasaltaod_p, aoduv_p, aodnir_p, &
        aoduvst_p, aodnirst_p) bind(c, name="modal_aer_opt_sw_zero_diagnostics_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
      type(c_ptr), value :: extinct_p, absorb_p, extinctuv_p, extinctnir_p
      type(c_ptr), value :: aodvis_p, aodvisst_p, aodabs_p, aodabsbc_p, ssavis_p
      type(c_ptr), value :: burdendust_p, burdenso4_p, burdenpom_p, burdensoa_p, burdenbc_p
      type(c_ptr), value :: burdenseasalt_p, dustaod_p, so4aod_p, pomaod_p, soaaod_p
      type(c_ptr), value :: bcaod_p, seasaltaod_p, aoduv_p, aodnir_p, aoduvst_p, aodnirst_p
   end subroutine modal_aer_opt_sw_zero_diagnostics_codon

   subroutine modal_aer_opt_sw_mode_diag_init_codon(ncol_c, burden_p, aodmode_p, dustaodmode_p) &
        bind(c, name="modal_aer_opt_sw_mode_diag_init_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c
      type(c_ptr), value :: burden_p, aodmode_p, dustaodmode_p
   end subroutine modal_aer_opt_sw_mode_diag_init_codon

   subroutine modal_aer_opt_sw_reset_layer_codon(ncol_c, dryvol_p, dustvol_p, scatdust_p, &
        absdust_p, hygrodust_p, scatso4_p, absso4_p, hygroso4_p, scatbc_p, absbc_p, hygrobc_p, &
        scatpom_p, abspom_p, hygropom_p, scatsoa_p, abssoa_p, hygrosoa_p, scatseasalt_p, &
        absseasalt_p, hygroseasalt_p, crefin_re_p, crefin_im_p) &
        bind(c, name="modal_aer_opt_sw_reset_layer_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c
      type(c_ptr), value :: dryvol_p, dustvol_p, scatdust_p, absdust_p, hygrodust_p
      type(c_ptr), value :: scatso4_p, absso4_p, hygroso4_p, scatbc_p, absbc_p, hygrobc_p
      type(c_ptr), value :: scatpom_p, abspom_p, hygropom_p, scatsoa_p, abssoa_p, hygrosoa_p
      type(c_ptr), value :: scatseasalt_p, absseasalt_p, hygroseasalt_p, crefin_re_p, crefin_im_p
   end subroutine modal_aer_opt_sw_reset_layer_codon

   subroutine modal_aer_opt_sw_species_volume_codon(ncol_c, pcols_c, k_c, specdens_c, &
        specrefr_c, specrefi_c, specmmr_p, vol_p, dryvol_p, crefin_re_p, crefin_im_p) &
        bind(c, name="modal_aer_opt_sw_species_volume_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, k_c
      real(c_double), value :: specdens_c, specrefr_c, specrefi_c
      type(c_ptr), value :: specmmr_p, vol_p, dryvol_p, crefin_re_p, crefin_im_p
   end subroutine modal_aer_opt_sw_species_volume_codon

   subroutine modal_aer_opt_sw_water_volume_codon(ncol_c, pcols_c, k_c, rhoh2o_c, &
        qaerwat_p, dryvol_p, watervol_p, wetvol_p) bind(c, name="modal_aer_opt_sw_water_volume_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, k_c
      real(c_double), value :: rhoh2o_c
      type(c_ptr), value :: qaerwat_p, dryvol_p, watervol_p, wetvol_p
   end subroutine modal_aer_opt_sw_water_volume_codon

   function modal_aer_opt_sw_has_negative_water_codon(ncol_c, watervol_p) result(has_c) &
        bind(c, name="modal_aer_opt_sw_has_negative_water_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c
      type(c_ptr), value :: watervol_p
      integer(c_int64_t) :: has_c
   end function modal_aer_opt_sw_has_negative_water_codon

   function modal_aer_opt_sw_water_refr_fastpath_codon(ncol_c, pcols_c, k_c, rhoh2o_c, &
        crefwsw_re_c, crefwsw_im_c, qaerwat_p, dryvol_p, watervol_p, wetvol_p, &
        crefin_re_p, crefin_im_p, refr_p, refi_p) result(has_negative_c) &
        bind(c, name="modal_aer_opt_sw_water_refr_fastpath_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, k_c
      real(c_double), value :: rhoh2o_c, crefwsw_re_c, crefwsw_im_c
      type(c_ptr), value :: qaerwat_p, dryvol_p, watervol_p, wetvol_p
      type(c_ptr), value :: crefin_re_p, crefin_im_p, refr_p, refi_p
      integer(c_int64_t) :: has_negative_c
   end function modal_aer_opt_sw_water_refr_fastpath_codon

   subroutine modal_aer_opt_sw_finalize_refr_codon(ncol_c, crefwsw_re_c, crefwsw_im_c, &
        watervol_p, wetvol_p, crefin_re_p, crefin_im_p, refr_p, refi_p) &
        bind(c, name="modal_aer_opt_sw_finalize_refr_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c
      real(c_double), value :: crefwsw_re_c, crefwsw_im_c
      type(c_ptr), value :: watervol_p, wetvol_p, crefin_re_p, crefin_im_p, refr_p, refi_p
   end subroutine modal_aer_opt_sw_finalize_refr_codon

   subroutine modal_aer_opt_sw_species_vis_diag_codon(ncol_c, pcols_c, k_c, spectype_code_c, &
        specrefr_c, specrefi_c, hygro_aer_c, specmmr_p, mass_p, vol_p, burden_p, burdendust_p, &
        burdenso4_p, burdenbc_p, burdenpom_p, burdensoa_p, burdenseasalt_p, dustvol_p, scatdust_p, &
        absdust_p, hygrodust_p, scatso4_p, absso4_p, hygroso4_p, scatbc_p, absbc_p, hygrobc_p, &
        scatpom_p, abspom_p, hygropom_p, scatsoa_p, abssoa_p, hygrosoa_p, scatseasalt_p, &
        absseasalt_p, hygroseasalt_p) bind(c, name="modal_aer_opt_sw_species_vis_diag_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, k_c, spectype_code_c
      real(c_double), value :: specrefr_c, specrefi_c, hygro_aer_c
      type(c_ptr), value :: specmmr_p, mass_p, vol_p, burden_p, burdendust_p, burdenso4_p
      type(c_ptr), value :: burdenbc_p, burdenpom_p, burdensoa_p, burdenseasalt_p, dustvol_p
      type(c_ptr), value :: scatdust_p, absdust_p, hygrodust_p, scatso4_p, absso4_p, hygroso4_p
      type(c_ptr), value :: scatbc_p, absbc_p, hygrobc_p, scatpom_p, abspom_p, hygropom_p
      type(c_ptr), value :: scatsoa_p, abssoa_p, hygrosoa_p, scatseasalt_p, absseasalt_p
      type(c_ptr), value :: hygroseasalt_p
   end subroutine modal_aer_opt_sw_species_vis_diag_codon

   subroutine modal_aer_opt_sw_optics_props_codon(ncol_c, pcols_c, k_c, ncoef_c, xrmax_c, &
        rhoh2o_c, radsurf_p, logradsurf_p, cheb_p, cext_p, cabs_p, casm_p, wetvol_p, mass_p, &
        pext_p, specpext_p, pabs_p, pasm_p, palb_p, dopaer_p) &
        bind(c, name="modal_aer_opt_sw_optics_props_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, k_c, ncoef_c
      real(c_double), value :: xrmax_c, rhoh2o_c
      type(c_ptr), value :: radsurf_p, logradsurf_p, cheb_p, cext_p, cabs_p, casm_p
      type(c_ptr), value :: wetvol_p, mass_p, pext_p, specpext_p, pabs_p, pasm_p, palb_p
      type(c_ptr), value :: dopaer_p
   end subroutine modal_aer_opt_sw_optics_props_codon

   subroutine modal_aer_opt_sw_mode_diag_night_codon(nnite_c, fillvalue_c, idxnite_p, burden_p, &
        aodmode_p, dustaodmode_p) bind(c, name="modal_aer_opt_sw_mode_diag_night_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: nnite_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: idxnite_p, burden_p, aodmode_p, dustaodmode_p
   end subroutine modal_aer_opt_sw_mode_diag_night_codon

   subroutine modal_aer_opt_sw_sum_diag_night_codon(nnite_c, pcols_c, pver_c, fillvalue_c, idxnite_p, &
        extinct_p, absorb_p, aodvis_p, aodabs_p, aodvisst_p) &
        bind(c, name="modal_aer_opt_sw_sum_diag_night_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: nnite_c, pcols_c, pver_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: idxnite_p, extinct_p, absorb_p, aodvis_p, aodabs_p, aodvisst_p
   end subroutine modal_aer_opt_sw_sum_diag_night_codon

   subroutine modal_aer_opt_sw_finalize_ssavis_codon(ncol_c, aodvis_p, ssavis_p) &
        bind(c, name="modal_aer_opt_sw_finalize_ssavis_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c
      type(c_ptr), value :: aodvis_p, ssavis_p
   end subroutine modal_aer_opt_sw_finalize_ssavis_codon

   subroutine modal_aer_opt_sw_climate_diag_night_codon(nnite_c, pcols_c, pver_c, fillvalue_c, &
        idxnite_p, ssavis_p, aoduv_p, aodnir_p, aoduvst_p, aodnirst_p, extinctuv_p, extinctnir_p, &
        burdendust_p, burdenso4_p, burdenpom_p, burdensoa_p, burdenbc_p, burdenseasalt_p, &
        aodabsbc_p, dustaod_p, so4aod_p, pomaod_p, soaaod_p, bcaod_p, seasaltaod_p) &
        bind(c, name="modal_aer_opt_sw_climate_diag_night_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: nnite_c, pcols_c, pver_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: idxnite_p, ssavis_p, aoduv_p, aodnir_p, aoduvst_p, aodnirst_p
      type(c_ptr), value :: extinctuv_p, extinctnir_p, burdendust_p, burdenso4_p, burdenpom_p
      type(c_ptr), value :: burdensoa_p, burdenbc_p, burdenseasalt_p, aodabsbc_p, dustaod_p
      type(c_ptr), value :: so4aod_p, pomaod_p, soaaod_p, bcaod_p, seasaltaod_p
   end subroutine modal_aer_opt_sw_climate_diag_night_codon

   subroutine modal_aer_opt_sw_accumulate_diagnostics_codon(ncol_c, pcols_c, k_c, do_uv_c, do_nir_c, &
        do_vis_c, crefwsw_re_c, crefwsw_im_c, troplev_p, mass_p, air_density_p, dopaer_p, pabs_p, &
        palb_p, wetvol_p, watervol_p, dustvol_p, scatdust_p, scatso4_p, scatbc_p, scatpom_p, &
        scatsoa_p, scatseasalt_p, absdust_p, absso4_p, absbc_p, abspom_p, abssoa_p, absseasalt_p, &
        hygrodust_p, hygroso4_p, hygrobc_p, hygropom_p, hygrosoa_p, hygroseasalt_p, extinctuv_p, &
        aoduv_p, aoduvst_p, extinctnir_p, aodnir_p, aodnirst_p, extinct_p, absorb_p, aodvis_p, &
        aodabs_p, aodmode_p, ssavis_p, aodvisst_p, dustaodmode_p, aodabsbc_p, dustaod_p, &
        so4aod_p, pomaod_p, soaaod_p, bcaod_p, seasaltaod_p) &
        bind(c, name="modal_aer_opt_sw_accumulate_diagnostics_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, k_c, do_uv_c, do_nir_c, do_vis_c
      real(c_double), value :: crefwsw_re_c, crefwsw_im_c
      type(c_ptr), value :: troplev_p, mass_p, air_density_p, dopaer_p, pabs_p, palb_p
      type(c_ptr), value :: wetvol_p, watervol_p, dustvol_p, scatdust_p, scatso4_p, scatbc_p
      type(c_ptr), value :: scatpom_p, scatsoa_p, scatseasalt_p, absdust_p, absso4_p, absbc_p
      type(c_ptr), value :: abspom_p, abssoa_p, absseasalt_p, hygrodust_p, hygroso4_p, hygrobc_p
      type(c_ptr), value :: hygropom_p, hygrosoa_p, hygroseasalt_p, extinctuv_p, aoduv_p, aoduvst_p
      type(c_ptr), value :: extinctnir_p, aodnir_p, aodnirst_p, extinct_p, absorb_p, aodvis_p
      type(c_ptr), value :: aodabs_p, aodmode_p, ssavis_p, aodvisst_p, dustaodmode_p, aodabsbc_p
      type(c_ptr), value :: dustaod_p, so4aod_p, pomaod_p, soaaod_p, bcaod_p, seasaltaod_p
   end subroutine modal_aer_opt_sw_accumulate_diagnostics_codon

   subroutine modal_aer_opt_sw_accumulate_tau_codon(ncol_c, pcols_c, pver_c, k_c, isw_c, &
        dopaer_p, palb_p, pasm_p, tauxar_p, wa_p, ga_p, fa_p) &
        bind(c, name="modal_aer_opt_sw_accumulate_tau_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, k_c, isw_c
      type(c_ptr), value :: dopaer_p, palb_p, pasm_p, tauxar_p, wa_p, ga_p, fa_p
   end subroutine modal_aer_opt_sw_accumulate_tau_codon

   function modal_aer_opt_sw_has_bad_dopaer_codon(ncol_c, dopaer_p) result(has_c) &
        bind(c, name="modal_aer_opt_sw_has_bad_dopaer_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c
      type(c_ptr), value :: dopaer_p
      integer(c_int64_t) :: has_c
   end function modal_aer_opt_sw_has_bad_dopaer_codon

   function modal_aer_opt_sw_optics_diag_tau_batch_codon(ncol_c, pcols_c, pver_c, k_c, isw_c, &
        ncoef_c, xrmax_c, rhoh2o_c, do_uv_c, do_nir_c, do_vis_c, crefwsw_re_c, crefwsw_im_c, &
        radsurf_p, logradsurf_p, cheb_p, cext_p, cabs_p, casm_p, wetvol_p, mass_p, pext_p, &
        specpext_p, pabs_p, pasm_p, palb_p, dopaer_p, troplev_p, air_density_p, watervol_p, &
        dustvol_p, scatdust_p, scatso4_p, scatbc_p, scatpom_p, scatsoa_p, scatseasalt_p, &
        absdust_p, absso4_p, absbc_p, abspom_p, abssoa_p, absseasalt_p, hygrodust_p, &
        hygroso4_p, hygrobc_p, hygropom_p, hygrosoa_p, hygroseasalt_p, extinctuv_p, aoduv_p, &
        aoduvst_p, extinctnir_p, aodnir_p, aodnirst_p, extinct_p, absorb_p, aodvis_p, &
        aodabs_p, aodmode_p, ssavis_p, aodvisst_p, dustaodmode_p, aodabsbc_p, dustaod_p, &
        so4aod_p, pomaod_p, soaaod_p, bcaod_p, seasaltaod_p, tauxar_p, wa_p, ga_p, fa_p) &
        result(has_c) bind(c, name="modal_aer_opt_sw_optics_diag_tau_batch_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, k_c, isw_c, ncoef_c
      integer(c_int64_t), value :: do_uv_c, do_nir_c, do_vis_c
      real(c_double), value :: xrmax_c, rhoh2o_c, crefwsw_re_c, crefwsw_im_c
      type(c_ptr), value :: radsurf_p, logradsurf_p, cheb_p, cext_p, cabs_p, casm_p
      type(c_ptr), value :: wetvol_p, mass_p, pext_p, specpext_p, pabs_p, pasm_p, palb_p
      type(c_ptr), value :: dopaer_p, troplev_p, air_density_p, watervol_p, dustvol_p
      type(c_ptr), value :: scatdust_p, scatso4_p, scatbc_p, scatpom_p, scatsoa_p, scatseasalt_p
      type(c_ptr), value :: absdust_p, absso4_p, absbc_p, abspom_p, abssoa_p, absseasalt_p
      type(c_ptr), value :: hygrodust_p, hygroso4_p, hygrobc_p, hygropom_p, hygrosoa_p
      type(c_ptr), value :: hygroseasalt_p, extinctuv_p, aoduv_p, aoduvst_p, extinctnir_p
      type(c_ptr), value :: aodnir_p, aodnirst_p, extinct_p, absorb_p, aodvis_p, aodabs_p
      type(c_ptr), value :: aodmode_p, ssavis_p, aodvisst_p, dustaodmode_p, aodabsbc_p
      type(c_ptr), value :: dustaod_p, so4aod_p, pomaod_p, soaaod_p, bcaod_p, seasaltaod_p
      type(c_ptr), value :: tauxar_p, wa_p, ga_p, fa_p
      integer(c_int64_t) :: has_c
   end function modal_aer_opt_sw_optics_diag_tau_batch_codon

   subroutine modal_aer_opt_lw_init_state_codon(ncol_c, pcols_c, pver_c, nlwbands_c, rga_c, &
        pdeldry_p, tauxar_p, mass_p) bind(c, name="modal_aer_opt_lw_init_state_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, nlwbands_c
      real(c_double), value :: rga_c
      type(c_ptr), value :: pdeldry_p, tauxar_p, mass_p
   end subroutine modal_aer_opt_lw_init_state_codon

   subroutine modal_aer_opt_lw_optics_props_codon(ncol_c, pcols_c, pver_c, k_c, ncoef_c, &
        rhoh2o_c, cheby_p, cabs_p, wetvol_p, mass_p, pabs_p, dopaer_p) &
        bind(c, name="modal_aer_opt_lw_optics_props_codon")
      use iso_c_binding, only: c_int64_t, c_double, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, k_c, ncoef_c
      real(c_double), value :: rhoh2o_c
      type(c_ptr), value :: cheby_p, cabs_p, wetvol_p, mass_p, pabs_p, dopaer_p
   end subroutine modal_aer_opt_lw_optics_props_codon

   subroutine modal_aer_opt_lw_accumulate_tau_codon(ncol_c, pcols_c, pver_c, k_c, ilw_c, &
        dopaer_p, tauxar_p) bind(c, name="modal_aer_opt_lw_accumulate_tau_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, k_c, ilw_c
      type(c_ptr), value :: dopaer_p, tauxar_p
   end subroutine modal_aer_opt_lw_accumulate_tau_codon
end interface

!===============================================================================
CONTAINS
!===============================================================================

subroutine modal_aer_opt_helpers_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (modal_aer_opt_helpers_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('MODAL_AER_OPT_HELPERS_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_modal_aer_opt_helpers_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_modal_aer_opt_helpers_impl = .false.
   end if

   modal_aer_opt_helpers_impl_selected = .true.

   if (masterproc) then
      if (use_native_modal_aer_opt_helpers_impl) then
         write(iulog,*) 'modal_aer_opt_helpers implementation = native'
      else
         write(iulog,*) 'modal_aer_opt_helpers implementation = codon'
      end if
   end if

end subroutine modal_aer_opt_helpers_select_impl

!===============================================================================

subroutine modal_aer_opt_helpers_proof_once()

   if (modal_aer_opt_helpers_proof_written) return
   modal_aer_opt_helpers_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'modal_aer_opt_helpers entered (modal aerosol size/interpolation/sw binterp/sw diagnostics/tau/' // &
           'layer reset/species volume/water refr/species vis/optics props helpers = codon)'
   end if

end subroutine modal_aer_opt_helpers_proof_once

!===============================================================================

subroutine modal_aer_opt_lw_helpers_proof_once()

   if (modal_aer_opt_lw_helpers_proof_written) return
   modal_aer_opt_lw_helpers_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'modal_aero_lw helpers entered (init state/optics props/tau = codon)'
   end if

end subroutine modal_aer_opt_lw_helpers_proof_once

!===============================================================================

subroutine modal_aer_opt_sw_guard_helpers_proof_once()

   if (modal_aer_opt_sw_guard_helpers_proof_written) return
   modal_aer_opt_sw_guard_helpers_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'modal_aero_sw guard helpers entered (layer reset/refr zero/' // &
           'water negative scan/dopaer scan = codon)'
   end if

end subroutine modal_aer_opt_sw_guard_helpers_proof_once

subroutine modal_aer_opt_sw_water_refr_proof_once()
   use spmd_utils, only: masterproc
   if (modal_aer_opt_sw_water_refr_proof_written) return
   modal_aer_opt_sw_water_refr_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'modal_aero_sw water/refr batch entered ' // &
           '(water volume/negative scan/finalize fast path = codon)'
   endif
end subroutine modal_aer_opt_sw_water_refr_proof_once

!===============================================================================

subroutine modal_aer_opt_sw_optics_tau_proof_once()
   use spmd_utils, only: masterproc
   if (modal_aer_opt_sw_optics_tau_proof_written) return
   modal_aer_opt_sw_optics_tau_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'modal_aero_sw optics/diagnostics/tau batch entered ' // &
           '(optics/diagnostics/dopaer scan/clean tau = codon)'
   endif
end subroutine modal_aer_opt_sw_optics_tau_proof_once

!===============================================================================

subroutine modal_aer_opt_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'modal_aer_opt_readnl'

   namelist /modal_aer_opt_nl/ water_refindex_file
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'modal_aer_opt_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, modal_aer_opt_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   call mpibcast(water_refindex_file, len(water_refindex_file), mpichar, 0, mpicom)
#endif


end subroutine modal_aer_opt_readnl

!===============================================================================

subroutine modal_aer_opt_init()

   use ioFileMod,        only: getfil
   use phys_control,     only: phys_getopts

   ! Local variables

   integer  :: i, m
   real(r8) :: rmmin, rmmax       ! min, max aerosol surface mode radius treated (m)
   character(len=256) :: locfile
   
   logical           :: history_amwg            ! output the variables used by the AMWG diag package
   logical           :: history_aero_optics     ! output aerosol optics diagnostics

   logical :: call_list(0:n_diag)
   integer :: ilist, nmodes, m_ncoef, m_prefr, m_prefi
   integer :: errcode

   character(len=*), parameter :: routine='modal_aer_opt_init'
   character(len=8) :: fldname
   character(len=128) :: lngname

   !----------------------------------------------------------------------------

   rmmin = 0.01e-6_r8
   rmmax = 25.e-6_r8
   xrmin = log(rmmin)
   xrmax = log(rmmax)

   ! Check that dimension sizes in the coefficient arrays used to
   ! parameterize aerosol radiative properties are consistent between this
   ! module and the mode physprop files.
   call rad_cnst_get_call_list(call_list)
   do ilist = 0, n_diag
      if (call_list(ilist)) then
         call rad_cnst_get_info(ilist, nmodes=nmodes)
         do m = 1, nmodes
            call rad_cnst_get_mode_props(ilist, m, ncoef=m_ncoef, prefr=m_prefr, prefi=m_prefi)
            if (m_ncoef /= ncoef .or. m_prefr /= prefr .or. m_prefi /= prefi) then
               write(iulog,*) routine//': ERROR - file and module values do not match:'
               write(iulog,*) '   ncoef:', ncoef, m_ncoef
               write(iulog,*) '   prefr:', prefr, m_prefr
               write(iulog,*) '   prefi:', prefi, m_prefi
               call endrun(routine//': ERROR - file and module values do not match')
            end if
         end do
      end if
   end do

   ! Initialize physics buffer indices for dgnumwet and qaerwat.  Note the implicit assumption
   ! that the loops over modes in the optics calculations will use the values for dgnumwet and qaerwat
   ! that are set in the aerosol_wet_intr code.
   dgnumwet_idx = pbuf_get_index('DGNUMWET',errcode)
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field DGNUMWET')
   end if
   qaerwat_idx  = pbuf_get_index('QAERWAT',errcode)
   if (errcode < 0) then
      call endrun(routine//' ERROR: cannot find physics buffer field QAERWAT')
   end if

   call getfil(water_refindex_file, locfile)
   call read_water_refindex(locfile)
   if (masterproc) write(iulog,*) "modal_aer_opt_init: read water refractive index file:", trim(locfile)

   call phys_getopts(history_amwg_out        = history_amwg, &
                     history_aero_optics_out = history_aero_optics )

   ! Add diagnostic fields to history output.

   call addfld ('EXTINCT','/m  ',pver,    'A','Aerosol extinction 550 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('EXTINCTUV','/m  ',pver,    'A','Aerosol extinction 350 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('EXTINCTNIR','/m  ',pver,    'A','Aerosol extinction 1020 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('ABSORB','/m  ',pver,    'A','Aerosol absorption',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODVIS','  ',1,    'A','Aerosol optical depth 550 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODVISst','  ',1,    'A','Stratospheric aerosol optical depth 550 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODUV','  ',1,    'A','Aerosol optical depth 350 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODUVst','  ',1,    'A','Stratospheric aerosol optical depth 350 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODNIR','  ',1,    'A','Aerosol optical depth 1020 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODNIRst','  ',1,    'A','Stratospheric aerosol optical depth 1020 nm',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODABS','  ',1,    'A','Aerosol absorption optical depth 550 nm',phys_decomp, flag_xyfill=.true.)

   call rad_cnst_get_info(0, nmodes=nmodes)

   do m = 1, nmodes

      write(fldname,'(a,i1)') 'BURDEN', m
      write(lngname,'(a,i1)') 'Aerosol burden mode ', m
      call addfld (fldname, 'kg/m2', 1, 'A', lngname, phys_decomp, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODMODE', m
      write(lngname,'(a,i1)') 'Aerosol optical depth 550 nm mode ', m
      call addfld (fldname, '  ', 1, 'A', lngname, phys_decomp, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

      write(fldname,'(a,i1)') 'AODDUST', m
      write(lngname,'(a,i1,a)') 'Aerosol optical depth 550 nm model ',m,' from dust'
      call addfld (fldname, '  ', 1, 'A', lngname, phys_decomp, flag_xyfill=.true.)
      if (m>3 .and. history_aero_optics) then
         call add_default (fldname, 1, ' ')
      endif

   enddo

   call addfld ('AODDUST','  ',1,    'A','Aerosol optical depth 550 nm from dust',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODSO4','  ',1,    'A','Aerosol optical depth 550 nm from SO4',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODPOM','  ',1,    'A','Aerosol optical depth 550 nm from POM',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODSOA','  ',1,    'A','Aerosol optical depth 550 nm from SOA',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODBC','  ',1,    'A','Aerosol optical depth 550 nm from BC',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODSS','  ',1,    'A','Aerosol optical depth 550 nm from seasalt',phys_decomp, flag_xyfill=.true.)
   call addfld ('AODABSBC','  ',1, 'A','Aerosol absorption optical depth 550 nm from BC',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENDUST','kg/m2'   ,1,  'A','Dust aerosol burden'        ,phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENSO4','kg/m2'    ,1,  'A','Sulfate aerosol burden'     ,phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENPOM','kg/m2'    ,1,  'A','POM aerosol burden'         ,phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENSOA','kg/m2'    ,1,  'A','SOA aerosol burden'         ,phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENBC','kg/m2'     ,1,  'A','Black carbon aerosol burden',phys_decomp, flag_xyfill=.true.)
   call addfld ('BURDENSEASALT','kg/m2',1,  'A','Seasalt aerosol burden'     ,phys_decomp, flag_xyfill=.true.)
   call addfld ('SSAVIS','  ',1,    'A','Aerosol singel-scatter albedo',phys_decomp, flag_xyfill=.true.)

 
   if (history_amwg) then 
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
      call add_default ('BURDEN1'      , 1, ' ')
      call add_default ('BURDEN2'      , 1, ' ')
      call add_default ('BURDEN3'      , 1, ' ')
      call add_default ('BURDENDUST'   , 1, ' ')
      call add_default ('BURDENSO4'    , 1, ' ')
      call add_default ('BURDENPOM'    , 1, ' ')
      call add_default ('BURDENSOA'    , 1, ' ')
      call add_default ('BURDENBC'     , 1, ' ')
      call add_default ('BURDENSEASALT', 1, ' ')
   end if

   if (history_aero_optics) then 
      call add_default ('AODDUST1'     , 1, ' ')
      call add_default ('AODDUST3'     , 1, ' ')
      call add_default ('ABSORB'       , 1, ' ')
      call add_default ('AODMODE1'     , 1, ' ')
      call add_default ('AODMODE2'     , 1, ' ')
      call add_default ('AODMODE3'     , 1, ' ')
      call add_default ('AODVIS'       , 1, ' ')
      call add_default ('AODUV'        , 1, ' ')
      call add_default ('AODNIR'       , 1, ' ')
      call add_default ('AODABS'       , 1, ' ')
      call add_default ('AODABSBC'     , 1, ' ')
      call add_default ('AODDUST'      , 1, ' ')
      call add_default ('AODSO4'       , 1, ' ')
      call add_default ('AODPOM'       , 1, ' ')
      call add_default ('AODSOA'       , 1, ' ')
      call add_default ('AODBC'        , 1, ' ')
      call add_default ('AODSS'        , 1, ' ')
      call add_default ('BURDEN1'      , 1, ' ')
      call add_default ('BURDEN2'      , 1, ' ')
      call add_default ('BURDEN3'      , 1, ' ')
      call add_default ('BURDENDUST'   , 1, ' ')
      call add_default ('BURDENSO4'    , 1, ' ')
      call add_default ('BURDENPOM'    , 1, ' ')
      call add_default ('BURDENSOA'    , 1, ' ')
      call add_default ('BURDENBC'     , 1, ' ')
      call add_default ('BURDENSEASALT', 1, ' ')
      call add_default ('SSAVIS'       , 1, ' ')
      call add_default ('EXTINCT'      , 1, ' ')
  end if

   do ilist = 1, n_diag
      if (call_list(ilist)) then
         
         call addfld ('EXTINCT'//diag(ilist),'/m  ', pver, 'A', &
              'Aerosol extinction',phys_decomp, flag_xyfill=.true.)
         call addfld ('ABSORB'//diag(ilist),'/m  ',  pver, 'A', &
              'Aerosol absorption',phys_decomp, flag_xyfill=.true.)
         call addfld ('AODVIS'//diag(ilist),'  ',       1, 'A', &
              'Aerosol optical depth 550 nm',phys_decomp, flag_xyfill=.true.)
         call addfld ('AODABS'//diag(ilist),'  ',       1, 'A', &
              'Aerosol absorption optical depth 550 nm',phys_decomp, flag_xyfill=.true.)

         if (history_aero_optics) then
            call add_default ('EXTINCT'//diag(ilist), 1, ' ')
            call add_default ('ABSORB'//diag(ilist),  1, ' ')
            call add_default ('AODVIS'//diag(ilist),  1, ' ')
            call add_default ('AODABS'//diag(ilist),  1, ' ')
         end if

      end if
   end do

end subroutine modal_aer_opt_init

!===============================================================================

subroutine modal_aero_sw(list_idx, state, pbuf, nnite, idxnite, &
                         tauxar, wa, ga, fa)

   ! calculates aerosol sw radiative properties
   
   use iso_c_binding, only: c_int64_t, c_loc
   use tropopause,     only : tropopause_find

   integer,             intent(in) :: list_idx       ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state          ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer,             intent(in) :: nnite          ! number of night columns
   integer, target,      intent(in) :: idxnite(nnite) ! local column indices of night columns

   real(r8), target, intent(out) :: tauxar(pcols,0:pver,nswbands) ! layer extinction optical depth
   real(r8), target, intent(out) :: wa(pcols,0:pver,nswbands)     ! layer single-scatter albedo
   real(r8), target, intent(out) :: ga(pcols,0:pver,nswbands)     ! asymmetry factor
   real(r8), target, intent(out) :: fa(pcols,0:pver,nswbands)     ! forward scattered fraction

   ! Local variables
   integer :: i, ifld, isw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec
   integer :: spectype_code
   integer, target :: troplev(pcols)

   real(r8), target :: mass(pcols,pver)        ! layer mass
   real(r8), target :: air_density(pcols,pver) ! (kg/m3)

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index
   character*32         :: spectype            ! species type
   real(r8)             :: hygro_aer           ! 

   real(r8), pointer :: dgnumwet(:,:)     ! number mode wet diameter
   real(r8), pointer :: qaerwat(:,:)      ! aerosol water (g/g)

   real(r8), pointer :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
   real(r8), pointer :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
   real(r8), pointer :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
   real(r8), pointer :: wetdens_m(:,:,:)  ! 

   real(r8) :: sigma_logr_aer         ! geometric standard deviation of number distribution
   real(r8), target :: radsurf(pcols,pver)    ! aerosol surface mode radius
   real(r8), target :: logradsurf(pcols,pver) ! log(aerosol surface mode radius)
   real(r8), target :: cheb(ncoef,pcols,pver)

   real(r8), target :: refr(pcols)     ! real part of refractive index
   real(r8), target :: refi(pcols)     ! imaginary part of refractive index
   complex(r8) :: crefin(pcols)   ! complex refractive index
   real(r8), target :: crefin_re(pcols), crefin_im(pcols)
   real(r8), pointer :: refrtabsw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitabsw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: extpsw(:,:,:,:) ! specific extinction
   real(r8), pointer :: abspsw(:,:,:,:) ! specific absorption
   real(r8), pointer :: asmpsw(:,:,:,:) ! asymmetry factor

   real(r8), target :: vol(pcols)    ! volume concentration of aerosol specie (m3/kg)
   real(r8), target :: dryvol(pcols) ! volume concentration of aerosol mode (m3/kg)
   real(r8), target :: watervol(pcols) ! volume concentration of water in each mode (m3/kg)
   real(r8), target :: wetvol(pcols)   ! volume concentration of wet mode (m3/kg)

   integer, target :: itab(pcols), jtab(pcols)
   real(r8), target :: ttab(pcols), utab(pcols)
   real(r8), target :: cext(pcols,ncoef), cabs(pcols,ncoef), casm(pcols,ncoef)
   real(r8), target :: pext(pcols)     ! parameterized specific extinction (m2/kg)
   real(r8), target :: specpext(pcols) ! specific extinction (m2/kg)
   real(r8), target :: dopaer(pcols)   ! aerosol optical depth in layer
   real(r8), target :: pabs(pcols)     ! parameterized specific absorption (m2/kg)
   real(r8), target :: pasm(pcols)     ! parameterized asymmetry factor
   real(r8), target :: palb(pcols)     ! parameterized single scattering albedo

   ! Diagnostics
   real(r8), target :: extinct(pcols,pver)
   real(r8), target :: extinctnir(pcols,pver)
   real(r8), target :: extinctuv(pcols,pver)
   real(r8), target :: absorb(pcols,pver)
   real(r8), target :: aodvis(pcols)        ! extinction optical depth
   real(r8), target :: aodvisst(pcols)      ! stratospheric extinction optical depth
   real(r8), target :: aodabs(pcols)        ! absorption optical depth

   real(r8), target :: aodabsbc(pcols)      ! absorption optical depth of BC

   real(r8), target :: ssavis(pcols)
   real(r8), target :: dustvol(pcols)      ! volume concentration of dust in aerosol mode (m3/kg)

   real(r8), target :: burden(pcols)
   real(r8), target :: burdendust(pcols), burdenso4(pcols), burdenbc(pcols), &
               burdenpom(pcols), burdensoa(pcols), burdenseasalt(pcols)

   real(r8), target :: aodmode(pcols)
   real(r8), target :: dustaodmode(pcols)  ! dust aod in aerosol mode

   real(r8) :: specrefr, specrefi
   real(r8), target :: scatdust(pcols), scatso4(pcols), scatbc(pcols), &
               scatpom(pcols), scatsoa(pcols), scatseasalt(pcols)
   real(r8), target :: absdust(pcols), absso4(pcols), absbc(pcols), &
               abspom(pcols), abssoa(pcols), absseasalt(pcols)
   real(r8), target :: hygrodust(pcols), hygroso4(pcols), hygrobc(pcols), &
               hygropom(pcols), hygrosoa(pcols), hygroseasalt(pcols)

   real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro
   real(r8) :: aodc                        ! aod of component

   ! total species AOD
   real(r8), target :: dustaod(pcols), so4aod(pcols), bcaod(pcols), &
               pomaod(pcols), soaaod(pcols), seasaltaod(pcols)




   logical :: savaervis ! true if visible wavelength (0.55 micron)
   logical :: savaernir ! true if near ir wavelength (~0.88 micron)
   logical :: savaeruv  ! true if uv wavelength (~0.35 micron)
   logical :: run_dopaer_diag
   logical :: run_water_fix

   real(r8), target :: aoduv(pcols)        ! extinction optical depth in uv
   real(r8), target :: aoduvst(pcols)      ! stratospheric extinction optical depth in uv
   real(r8), target :: aodnir(pcols)       ! extinction optical depth in nir
   real(r8), target :: aodnirst(pcols)     ! stratospheric extinction optical depth in nir


   character(len=32) :: outname

   ! debug output
   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf            ! volume fraction of insoluble aerosol
   character(len=*), parameter :: subname = 'modal_aero_sw'
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! initialize output variables
   call modal_aer_opt_helpers_select_impl()
   if (use_native_modal_aer_opt_helpers_impl) then
      tauxar(:ncol,:,:) = 0._r8
      wa(:ncol,:,:)     = 0._r8
      ga(:ncol,:,:)     = 0._r8
      fa(:ncol,:,:)     = 0._r8

      ! zero'th layer does not contain aerosol
      tauxar(1:ncol,0,:)  = 0._r8
      wa(1:ncol,0,:)      = 0.925_r8
      ga(1:ncol,0,:)      = 0.850_r8
      fa(1:ncol,0,:)      = 0.7225_r8

      mass(:ncol,:)        = state%pdeldry(:ncol,:)*rga
      air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))
   else
      call modal_aer_opt_helpers_proof_once()
      call modal_aer_opt_sw_init_state_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(nswbands, c_int64_t), rga, rair, c_loc(state%pdeldry(1,1)), &
           c_loc(state%pmid(1,1)), c_loc(state%t(1,1)), c_loc(tauxar(1,0,1)), c_loc(wa(1,0,1)), &
           c_loc(ga(1,0,1)), c_loc(fa(1,0,1)), c_loc(mass(1,1)), c_loc(air_density(1,1)))
   end if

   ! diagnostics for visible band summed over modes
   if (use_native_modal_aer_opt_helpers_impl) then
      extinct(1:ncol,:)     = 0.0_r8
      absorb(1:ncol,:)      = 0.0_r8
      aodvis(1:ncol)        = 0.0_r8
      aodvisst(1:ncol)      = 0.0_r8
      aodabs(1:ncol)        = 0.0_r8
      burdendust(:ncol)     = 0.0_r8
      burdenso4(:ncol)      = 0.0_r8
      burdenpom(:ncol)      = 0.0_r8
      burdensoa(:ncol)      = 0.0_r8
      burdenbc(:ncol)       = 0.0_r8
      burdenseasalt(:ncol)  = 0.0_r8
      ssavis(1:ncol)        = 0.0_r8

      aodabsbc(:ncol)       = 0.0_r8
      dustaod(:ncol)        = 0.0_r8
      so4aod(:ncol)         = 0.0_r8
      pomaod(:ncol)         = 0.0_r8
      soaaod(:ncol)         = 0.0_r8
      bcaod(:ncol)          = 0.0_r8
      seasaltaod(:ncol)     = 0.0_r8

      ! diags for other bands
      extinctuv(1:ncol,:)   = 0.0_r8
      extinctnir(1:ncol,:)  = 0.0_r8
      aoduv(:ncol)          = 0.0_r8
      aodnir(:ncol)         = 0.0_r8
      aoduvst(:ncol)        = 0.0_r8
      aodnirst(:ncol)       = 0.0_r8
   else
      call modal_aer_opt_helpers_proof_once()
      call modal_aer_opt_sw_zero_diagnostics_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), c_loc(extinct(1,1)), c_loc(absorb(1,1)), c_loc(extinctuv(1,1)), &
           c_loc(extinctnir(1,1)), c_loc(aodvis(1)), c_loc(aodvisst(1)), c_loc(aodabs(1)), &
           c_loc(aodabsbc(1)), c_loc(ssavis(1)), c_loc(burdendust(1)), c_loc(burdenso4(1)), &
           c_loc(burdenpom(1)), c_loc(burdensoa(1)), c_loc(burdenbc(1)), c_loc(burdenseasalt(1)), &
           c_loc(dustaod(1)), c_loc(so4aod(1)), c_loc(pomaod(1)), c_loc(soaaod(1)), &
           c_loc(bcaod(1)), c_loc(seasaltaod(1)), c_loc(aoduv(1)), c_loc(aodnir(1)), &
           c_loc(aoduvst(1)), c_loc(aodnirst(1)))
   end if
   call tropopause_find(state, troplev)

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx == 0) then
      ! water uptake and wet radius for the climate list has already been calculated
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_m)
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat_m)
   else
      ! If doing a diagnostic calculation then need to calculate the wet radius
      ! and water uptake for the diagnostic modes
      call modal_aero_calcsize_diag(state, pbuf, list_idx, dgnumdry_m)  
      call modal_aero_wateruptake_dr(state, pbuf, list_idx, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m)
   endif

   do m = 1, nmodes

      ! diagnostics for visible band for each mode
      if (use_native_modal_aer_opt_helpers_impl) then
         burden(:ncol)       = 0._r8
         aodmode(1:ncol)     = 0.0_r8
         dustaodmode(1:ncol) = 0.0_r8
      else
         call modal_aer_opt_helpers_proof_once()
         call modal_aer_opt_sw_mode_diag_init_codon(int(ncol, c_int64_t), c_loc(burden(1)), &
              c_loc(aodmode(1)), c_loc(dustaodmode(1)))
      end if

      dgnumwet => dgnumwet_m(:,:,m)
      qaerwat  => qaerwat_m(:,:,m)

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtabsw=refrtabsw , &
         refitabsw=refitabsw, extpsw=extpsw, abspsw=abspsw, asmpsw=asmpsw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

      do isw = 1, nswbands
         savaervis = (isw .eq. idx_sw_diag)
         savaeruv  = (isw .eq. idx_uv_diag)
         savaernir = (isw .eq. idx_nir_diag)

         do k = top_lev, pver

            if (use_native_modal_aer_opt_helpers_impl) then
               ! form bulk refractive index
               crefin(:ncol) = (0._r8, 0._r8)
               crefin_re(:ncol) = 0._r8
               crefin_im(:ncol) = 0._r8
               dryvol(:ncol) = 0._r8
               dustvol(:ncol) = 0._r8

               scatdust(:ncol)     = 0._r8
               absdust(:ncol)      = 0._r8
               hygrodust(:ncol)    = 0._r8
               scatso4(:ncol)      = 0._r8
               absso4(:ncol)       = 0._r8
               hygroso4(:ncol)     = 0._r8
               scatbc(:ncol)       = 0._r8
               absbc(:ncol)        = 0._r8
               hygrobc(:ncol)      = 0._r8
               scatpom(:ncol)      = 0._r8
               abspom(:ncol)       = 0._r8
               hygropom(:ncol)     = 0._r8
               scatsoa(:ncol)      = 0._r8
               abssoa(:ncol)       = 0._r8
               hygrosoa(:ncol)     = 0._r8
               scatseasalt(:ncol)  = 0._r8
               absseasalt(:ncol)   = 0._r8
               hygroseasalt(:ncol) = 0._r8
            else
               call modal_aer_opt_helpers_proof_once()
               call modal_aer_opt_sw_guard_helpers_proof_once()
               call modal_aer_opt_sw_reset_layer_codon(int(ncol, c_int64_t), c_loc(dryvol(1)), &
                    c_loc(dustvol(1)), c_loc(scatdust(1)), c_loc(absdust(1)), c_loc(hygrodust(1)), &
                    c_loc(scatso4(1)), c_loc(absso4(1)), c_loc(hygroso4(1)), c_loc(scatbc(1)), &
                    c_loc(absbc(1)), c_loc(hygrobc(1)), c_loc(scatpom(1)), c_loc(abspom(1)), &
                    c_loc(hygropom(1)), c_loc(scatsoa(1)), c_loc(abssoa(1)), c_loc(hygrosoa(1)), &
                    c_loc(scatseasalt(1)), c_loc(absseasalt(1)), c_loc(hygroseasalt(1)), &
                    c_loc(crefin_re(1)), c_loc(crefin_im(1)))
            end if

            ! aerosol species loop
            do l = 1, nspec
               call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
               call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                           refindex_aer_sw=specrefindex, spectype=spectype, &
                                           hygro_aer=hygro_aer)

               if (.not. use_native_modal_aer_opt_helpers_impl) then
                  specrefr = real(specrefindex(isw))
                  specrefi = aimag(specrefindex(isw))
                  call modal_aer_opt_helpers_proof_once()
                  call modal_aer_opt_sw_species_volume_codon(int(ncol, c_int64_t), &
                       int(pcols, c_int64_t), int(k, c_int64_t), specdens, specrefr, specrefi, &
                       c_loc(specmmr(1,1)), c_loc(vol(1)), c_loc(dryvol(1)), c_loc(crefin_re(1)), &
                       c_loc(crefin_im(1)))
               else
                  do i = 1, ncol
                     vol(i)      = specmmr(i,k)/specdens
                     dryvol(i)   = dryvol(i) + vol(i)
                     crefin(i)   = crefin(i) + vol(i)*specrefindex(isw)
                  end do
               end if

               ! compute some diagnostics for visible band only
               if (savaervis) then

                  specrefr = real(specrefindex(isw))
                  specrefi = aimag(specrefindex(isw))

                  if (.not. use_native_modal_aer_opt_helpers_impl) then
                     select case (trim(spectype))
                     case ('dust')
                        spectype_code = 1
                     case ('sulfate')
                        spectype_code = 2
                     case ('black-c')
                        spectype_code = 3
                     case ('p-organic')
                        spectype_code = 4
                     case ('s-organic')
                        spectype_code = 5
                     case ('seasalt')
                        spectype_code = 6
                     case default
                        spectype_code = 0
                     end select
                     call modal_aer_opt_helpers_proof_once()
                     call modal_aer_opt_sw_species_vis_diag_codon(int(ncol, c_int64_t), &
                          int(pcols, c_int64_t), int(k, c_int64_t), int(spectype_code, c_int64_t), &
                          specrefr, specrefi, hygro_aer, c_loc(specmmr(1,1)), c_loc(mass(1,1)), &
                          c_loc(vol(1)), c_loc(burden(1)), c_loc(burdendust(1)), c_loc(burdenso4(1)), &
                          c_loc(burdenbc(1)), c_loc(burdenpom(1)), c_loc(burdensoa(1)), &
                          c_loc(burdenseasalt(1)), c_loc(dustvol(1)), c_loc(scatdust(1)), &
                          c_loc(absdust(1)), c_loc(hygrodust(1)), c_loc(scatso4(1)), c_loc(absso4(1)), &
                          c_loc(hygroso4(1)), c_loc(scatbc(1)), c_loc(absbc(1)), c_loc(hygrobc(1)), &
                          c_loc(scatpom(1)), c_loc(abspom(1)), c_loc(hygropom(1)), c_loc(scatsoa(1)), &
                          c_loc(abssoa(1)), c_loc(hygrosoa(1)), c_loc(scatseasalt(1)), &
                          c_loc(absseasalt(1)), c_loc(hygroseasalt(1)))
                  else
                     do i = 1, ncol
                        burden(i) = burden(i) + specmmr(i,k)*mass(i,k)
                     end do

                     if (trim(spectype) == 'dust') then
                        do i = 1, ncol
                           burdendust(i) = burdendust(i) + specmmr(i,k)*mass(i,k)
                           dustvol(i)    = vol(i)
                           scatdust(i)   = vol(i)*specrefr
                           absdust(i)    = -vol(i)*specrefi
                           hygrodust(i)  = vol(i)*hygro_aer
                        end do
                     end if

                     if (trim(spectype) == 'sulfate') then
                        do i = 1, ncol
                           burdenso4(i) = burdenso4(i) + specmmr(i,k)*mass(i,k)
                           scatso4(i)   = vol(i)*specrefr
                           absso4(i)    = -vol(i)*specrefi
                           hygroso4(i)  = vol(i)*hygro_aer
                        end do
                     end if
                     if (trim(spectype) == 'black-c') then
                        do i = 1, ncol
                           burdenbc(i) = burdenbc(i) + specmmr(i,k)*mass(i,k)
                           scatbc(i)   = vol(i)*specrefr
                           absbc(i)    = -vol(i)*specrefi
                           hygrobc(i)  = vol(i)*hygro_aer
                      end do
                     end if
                     if (trim(spectype) == 'p-organic') then
                        do i = 1, ncol
                           burdenpom(i) = burdenpom(i) + specmmr(i,k)*mass(i,k)
                           scatpom(i)   = vol(i)*specrefr
                           abspom(i)    = -vol(i)*specrefi
                           hygropom(i)  = vol(i)*hygro_aer
                         end do
                     end if
                     if (trim(spectype) == 's-organic') then
                        do i = 1, ncol
                           burdensoa(i) = burdensoa(i) + specmmr(i,k)*mass(i,k)
                           scatsoa(i)   = vol(i)*specrefr
                           abssoa(i)    = -vol(i)*specrefi
                           hygrosoa(i)  = vol(i)*hygro_aer
                        end do
                     end if
                     if (trim(spectype) == 'seasalt') then
                        do i = 1, ncol
                           burdenseasalt(i) = burdenseasalt(i) + specmmr(i,k)*mass(i,k)
                           scatseasalt(i)   = vol(i)*specrefr
                           absseasalt(i)    = -vol(i)*specrefi
                           hygroseasalt(i)  = vol(i)*hygro_aer
                         end do
                     end if
                  end if

               end if
            end do ! species loop

            if (.not. use_native_modal_aer_opt_helpers_impl) then
               call modal_aer_opt_helpers_proof_once()
               call modal_aer_opt_sw_water_refr_proof_once()
               run_water_fix = modal_aer_opt_sw_water_refr_fastpath_codon(int(ncol, c_int64_t), &
                    int(pcols, c_int64_t), int(k, c_int64_t), rhoh2o, real(crefwsw(isw), r8), &
                    aimag(crefwsw(isw)), c_loc(qaerwat(1,1)), c_loc(dryvol(1)), c_loc(watervol(1)), &
                    c_loc(wetvol(1)), c_loc(crefin_re(1)), c_loc(crefin_im(1)), c_loc(refr(1)), &
                    c_loc(refi(1))) /= 0_c_int64_t
               if (run_water_fix) then
                  do i = 1, ncol
                     if (watervol(i) < 0._r8) then
                        if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                           write(iulog,'(a,2e10.2,a)') 'watervol,wetvol=', &
                              watervol(i), wetvol(i), ' in '//subname
                        end if
                        watervol(i) = 0._r8
                        wetvol(i) = dryvol(i)
                     end if
                  end do

                  call modal_aer_opt_helpers_proof_once()
                  call modal_aer_opt_sw_finalize_refr_codon(int(ncol, c_int64_t), real(crefwsw(isw), r8), &
                       aimag(crefwsw(isw)), c_loc(watervol(1)), c_loc(wetvol(1)), c_loc(crefin_re(1)), &
                       c_loc(crefin_im(1)), c_loc(refr(1)), c_loc(refi(1)))
               end if
            else
               do i = 1, ncol
                  watervol(i) = qaerwat(i,k)/rhoh2o
                  wetvol(i) = watervol(i) + dryvol(i)
                  if (watervol(i) < 0._r8) then
                     if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                        write(iulog,'(a,2e10.2,a)') 'watervol,wetvol=', &
                           watervol(i), wetvol(i), ' in '//subname
                     end if
                     watervol(i) = 0._r8
                     wetvol(i) = dryvol(i)
                  end if

                  ! volume mixing
                  crefin(i) = crefin(i) + watervol(i)*crefwsw(isw)
                  crefin(i) = crefin(i)/max(wetvol(i),1.e-60_r8)
                  refr(i)   = real(crefin(i))
                  refi(i)   = abs(aimag(crefin(i)))
               end do
            end if

            ! call t_startf('binterp')

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            if (.not. use_native_modal_aer_opt_helpers_impl) then
               call modal_aer_opt_helpers_proof_once()
               call modal_aer_opt_sw_binterp3_codon(int(pcols, c_int64_t), int(ncol, c_int64_t), &
                    int(ncoef, c_int64_t), int(prefr, c_int64_t), int(prefi, c_int64_t), &
                    c_loc(extpsw(1,1,1,isw)), c_loc(abspsw(1,1,1,isw)), c_loc(asmpsw(1,1,1,isw)), &
                    c_loc(refr(1)), c_loc(refi(1)), c_loc(refrtabsw(1,isw)), c_loc(refitabsw(1,isw)), &
                    c_loc(itab(1)), c_loc(jtab(1)), c_loc(ttab(1)), c_loc(utab(1)), c_loc(cext(1,1)), &
                    c_loc(cabs(1,1)), c_loc(casm(1,1)))
            else
               itab(:ncol) = 0
               call binterp(extpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                            refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                            itab, jtab, ttab, utab, cext)
               call binterp(abspsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                            refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                            itab, jtab, ttab, utab, cabs)
               call binterp(asmpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                            refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                            itab, jtab, ttab, utab, casm)
            end if

            ! call t_stopf('binterp')

            ! parameterized optical properties, diagnostics, dopaer scan, and clean tau accumulate
            if (.not. use_native_modal_aer_opt_helpers_impl) then
               call modal_aer_opt_helpers_proof_once()
               call modal_aer_opt_sw_optics_tau_proof_once()
               run_dopaer_diag = modal_aer_opt_sw_optics_diag_tau_batch_codon( &
                    int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
                    int(k, c_int64_t), int(isw, c_int64_t), int(ncoef, c_int64_t), xrmax, rhoh2o, &
                    merge(1_c_int64_t, 0_c_int64_t, savaeruv), &
                    merge(1_c_int64_t, 0_c_int64_t, savaernir), &
                    merge(1_c_int64_t, 0_c_int64_t, savaervis), real(crefwsw(isw), r8), &
                    aimag(crefwsw(isw)), c_loc(radsurf(1,1)), c_loc(logradsurf(1,1)), &
                    c_loc(cheb(1,1,1)), c_loc(cext(1,1)), c_loc(cabs(1,1)), c_loc(casm(1,1)), &
                    c_loc(wetvol(1)), c_loc(mass(1,1)), c_loc(pext(1)), c_loc(specpext(1)), &
                    c_loc(pabs(1)), c_loc(pasm(1)), c_loc(palb(1)), c_loc(dopaer(1)), &
                    c_loc(troplev(1)), c_loc(air_density(1,1)), c_loc(watervol(1)), &
                    c_loc(dustvol(1)), c_loc(scatdust(1)), c_loc(scatso4(1)), c_loc(scatbc(1)), &
                    c_loc(scatpom(1)), c_loc(scatsoa(1)), c_loc(scatseasalt(1)), &
                    c_loc(absdust(1)), c_loc(absso4(1)), c_loc(absbc(1)), c_loc(abspom(1)), &
                    c_loc(abssoa(1)), c_loc(absseasalt(1)), c_loc(hygrodust(1)), &
                    c_loc(hygroso4(1)), c_loc(hygrobc(1)), c_loc(hygropom(1)), &
                    c_loc(hygrosoa(1)), c_loc(hygroseasalt(1)), c_loc(extinctuv(1,1)), &
                    c_loc(aoduv(1)), c_loc(aoduvst(1)), c_loc(extinctnir(1,1)), &
                    c_loc(aodnir(1)), c_loc(aodnirst(1)), c_loc(extinct(1,1)), &
                    c_loc(absorb(1,1)), c_loc(aodvis(1)), c_loc(aodabs(1)), c_loc(aodmode(1)), &
                    c_loc(ssavis(1)), c_loc(aodvisst(1)), c_loc(dustaodmode(1)), &
                    c_loc(aodabsbc(1)), c_loc(dustaod(1)), c_loc(so4aod(1)), c_loc(pomaod(1)), &
                    c_loc(soaaod(1)), c_loc(bcaod(1)), c_loc(seasaltaod(1)), c_loc(tauxar(1,0,1)), &
                    c_loc(wa(1,0,1)), c_loc(ga(1,0,1)), c_loc(fa(1,0,1))) /= 0_c_int64_t
            else
               do i=1,ncol

                  if (logradsurf(i,k) .le. xrmax) then
                     pext(i) = 0.5_r8*cext(i,1)
                     do nc = 2, ncoef
                        pext(i) = pext(i) + cheb(nc,i,k)*cext(i,nc)
                     enddo
                     pext(i) = exp(pext(i))
                  else
                     pext(i) = 1.5_r8/(radsurf(i,k)*rhoh2o) ! geometric optics
                  endif

                  ! convert from m2/kg water to m2/kg aerosol
                  specpext(i) = pext(i)
                  pext(i) = pext(i)*wetvol(i)*rhoh2o
                  pabs(i) = 0.5_r8*cabs(i,1)
                  pasm(i) = 0.5_r8*casm(i,1)
                  do nc = 2, ncoef
                     pabs(i) = pabs(i) + cheb(nc,i,k)*cabs(i,nc)
                     pasm(i) = pasm(i) + cheb(nc,i,k)*casm(i,nc)
                  enddo
                  pabs(i) = pabs(i)*wetvol(i)*rhoh2o
                  pabs(i) = max(0._r8,pabs(i))
                  pabs(i) = min(pext(i),pabs(i))

                  palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
                  palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)

                  dopaer(i) = pext(i)*mass(i,k)
               end do
               if (savaeruv) then
                  do i = 1, ncol
                    extinctuv(i,k) = extinctuv(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                    aoduv(i) = aoduv(i) + dopaer(i)
                     if (k.le.troplev(i)) then
                       aoduvst(i) = aoduvst(i) + dopaer(i)
                     end if
                  end do
               end if

               if (savaernir) then
                  do i = 1, ncol
                     extinctnir(i,k) = extinctnir(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                     aodnir(i) = aodnir(i) + dopaer(i)
                     if (k.le.troplev(i)) then
                       aodnirst(i) = aodnirst(i) + dopaer(i)
                     end if
                  end do
               endif

               ! Save aerosol optical depth at longest visible wavelength
               ! sum over layers
               if (savaervis) then
                  ! aerosol extinction (/m)
                  do i = 1, ncol
                     extinct(i,k) = extinct(i,k) + dopaer(i)*air_density(i,k)/mass(i,k)
                     absorb(i,k)  = absorb(i,k) + pabs(i)*air_density(i,k)
                     aodvis(i)    = aodvis(i) + dopaer(i)
                     aodabs(i)    = aodabs(i) + pabs(i)*mass(i,k)
                     aodmode(i)   = aodmode(i) + dopaer(i)
                     ssavis(i)    = ssavis(i) + dopaer(i)*palb(i)
                     if (k.le.troplev(i)) then
                       aodvisst(i) = aodvisst(i) + dopaer(i)
                     end if

                     if (wetvol(i) > 1.e-40_r8) then

                        dustaodmode(i) = dustaodmode(i) + dopaer(i)*dustvol(i)/wetvol(i)

                        ! partition optical depth into contributions from each constituent
                        ! assume contribution is proportional to refractive index X volume

                        scath2o        = watervol(i)*real(crefwsw(isw))
                        absh2o         = -watervol(i)*aimag(crefwsw(isw))
                        sumscat        = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                                         scatdust(i) + scatseasalt(i) + scath2o
                        sumabs         = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                                         absdust(i) + absseasalt(i) + absh2o
                        sumhygro       = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                                         hygrodust(i) + hygroseasalt(i)

                        scatdust(i)    = (scatdust(i) + scath2o*hygrodust(i)/sumhygro)/sumscat
                        absdust(i)     = (absdust(i) + absh2o*hygrodust(i)/sumhygro)/sumabs

                        scatso4(i)     = (scatso4(i) + scath2o*hygroso4(i)/sumhygro)/sumscat
                        absso4(i)      = (absso4(i) + absh2o*hygroso4(i)/sumhygro)/sumabs

                        scatpom(i)     = (scatpom(i) + scath2o*hygropom(i)/sumhygro)/sumscat
                        abspom(i)      = (abspom(i) + absh2o*hygropom(i)/sumhygro)/sumabs

                        scatsoa(i)     = (scatsoa(i) + scath2o*hygrosoa(i)/sumhygro)/sumscat
                        abssoa(i)      = (abssoa(i) + absh2o*hygrosoa(i)/sumhygro)/sumabs

                        scatbc(i)      = (scatbc(i) + scath2o*hygrobc(i)/sumhygro)/sumscat
                        absbc(i)       = (absbc(i) + absh2o*hygrobc(i)/sumhygro)/sumabs

                        scatseasalt(i) = (scatseasalt(i) + scath2o*hygroseasalt(i)/sumhygro)/sumscat
                        absseasalt(i)  = (absseasalt(i) + absh2o*hygroseasalt(i)/sumhygro)/sumabs

                        aodabsbc(i)    = aodabsbc(i) + absbc(i)*dopaer(i)*(1.0_r8-palb(i))

                        aodc           = (absdust(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust(i))*dopaer(i)
                        dustaod(i)     = dustaod(i) + aodc

                        aodc           = (absso4(i)*(1.0_r8 - palb(i)) + palb(i)*scatso4(i))*dopaer(i)
                        so4aod(i)      = so4aod(i) + aodc

                        aodc           = (abspom(i)*(1.0_r8 - palb(i)) + palb(i)*scatpom(i))*dopaer(i)
                        pomaod(i)      = pomaod(i) + aodc

                        aodc           = (abssoa(i)*(1.0_r8 - palb(i)) + palb(i)*scatsoa(i))*dopaer(i)
                        soaaod(i)      = soaaod(i) + aodc

                        aodc           = (absbc(i)*(1.0_r8 - palb(i)) + palb(i)*scatbc(i))*dopaer(i)
                        bcaod(i)       = bcaod(i) + aodc

                        aodc           = (absseasalt(i)*(1.0_r8 - palb(i)) + palb(i)*scatseasalt(i))*dopaer(i)
                        seasaltaod(i)  = seasaltaod(i) + aodc

                     endif

                  end do
               endif
               run_dopaer_diag = .true.
            end if

            if (run_dopaer_diag) then
               do i = 1, ncol

                  if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 30._r8)) then

                     if (dopaer(i) <= -1.e-10_r8) then
                        write(iulog,*) "ERROR: Negative aerosol optical depth &
                             &in this layer."
                     else
                        write(iulog,*) "WARNING: Aerosol optical depth is &
                             &unreasonably high in this layer."
                     end if

                     write(iulog,*) 'dopaer(', i, ',', k, ',', m, ',', lchnk, ')=', dopaer(i)
                     ! write(iulog,*) 'itab,jtab,ttab,utab=',itab(i),jtab(i),ttab(i),utab(i)
                     write(iulog,*) 'k=', k, ' pext=', pext(i), ' specext=', specpext(i)
                     write(iulog,*) 'wetvol=', wetvol(i), ' dryvol=', dryvol(i), ' watervol=', watervol(i)
                     ! write(iulog,*) 'cext=',(cext(i,l),l=1,ncoef)
                     ! write(iulog,*) 'crefin=',crefin(i)
                     write(iulog,*) 'nspec=', nspec
                     ! write(iulog,*) 'cheb=', (cheb(nc,m,i,k),nc=2,ncoef)
                     do l = 1, nspec
                        call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
                        call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                    refindex_aer_sw=specrefindex)
                        volf = specmmr(i,k)/specdens
                        write(iulog,*) 'l=', l, 'vol(l)=', volf
                        write(iulog,*) 'isw=', isw, 'specrefindex(isw)=', specrefindex(isw)
                        write(iulog,*) 'specdens=', specdens
                     end do

                     nerr_dopaer = nerr_dopaer + 1
!                  if (nerr_dopaer >= nerrmax_dopaer) then
                     if (dopaer(i) < -1.e-10_r8) then
                        write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                        call endrun('exit from '//subname)
                     end if

                  end if
               end do
            end if

            if (.not. use_native_modal_aer_opt_helpers_impl) then
               if (run_dopaer_diag) then
                  call modal_aer_opt_helpers_proof_once()
                  call modal_aer_opt_sw_accumulate_tau_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
                       int(pver, c_int64_t), int(k, c_int64_t), int(isw, c_int64_t), c_loc(dopaer(1)), &
                       c_loc(palb(1)), c_loc(pasm(1)), c_loc(tauxar(1,0,1)), c_loc(wa(1,0,1)), &
                       c_loc(ga(1,0,1)), c_loc(fa(1,0,1)))
               end if
            else
               do i=1,ncol
                  tauxar(i,k,isw) = tauxar(i,k,isw) + dopaer(i)
                  wa(i,k,isw)     = wa(i,k,isw)     + dopaer(i)*palb(i)
                  ga(i,k,isw)     = ga(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)
                  fa(i,k,isw)     = fa(i,k,isw)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
               end do
            end if

         end do ! pver

      end do ! sw bands

      ! mode diagnostics
      ! The diagnostics are currently only output for the climate list.  Code mods will
      ! be necessary to provide output for the rad_diag lists.
      if (list_idx == 0) then
         if (use_native_modal_aer_opt_helpers_impl) then
            do i = 1, nnite
               burden(idxnite(i))  = fillvalue
               aodmode(idxnite(i)) = fillvalue
               dustaodmode(idxnite(i)) = fillvalue
            end do
         else if (nnite > 0) then
            call modal_aer_opt_helpers_proof_once()
            call modal_aer_opt_sw_mode_diag_night_codon(int(nnite, c_int64_t), fillvalue, &
                 c_loc(idxnite(1)), c_loc(burden(1)), c_loc(aodmode(1)), c_loc(dustaodmode(1)))
         end if

         write(outname,'(a,i1)') 'BURDEN', m
         call outfld(trim(outname), burden, pcols, lchnk)

         write(outname,'(a,i1)') 'AODMODE', m
         call outfld(trim(outname), aodmode, pcols, lchnk)

         write(outname,'(a,i1)') 'AODDUST', m
         call outfld(trim(outname), dustaodmode, pcols, lchnk)

      end if

   end do ! nmodes

   if (list_idx > 0) then
      deallocate(dgnumdry_m)
      deallocate(dgnumwet_m)
      deallocate(qaerwat_m)
      deallocate(wetdens_m)
   end if

   ! Output visible band diagnostics for quantities summed over the modes
   ! These fields are put out for diagnostic lists as well as the climate list.
   if (use_native_modal_aer_opt_helpers_impl) then
      do i = 1, nnite
         extinct(idxnite(i),:) = fillvalue
         absorb(idxnite(i),:)  = fillvalue
         aodvis(idxnite(i))    = fillvalue
         aodabs(idxnite(i))    = fillvalue
         aodvisst(idxnite(i))  = fillvalue
      end do
   else if (nnite > 0) then
      call modal_aer_opt_helpers_proof_once()
      call modal_aer_opt_sw_sum_diag_night_codon(int(nnite, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), fillvalue, c_loc(idxnite(1)), c_loc(extinct(1,1)), &
           c_loc(absorb(1,1)), c_loc(aodvis(1)), c_loc(aodabs(1)), c_loc(aodvisst(1)))
   end if

   call outfld('EXTINCT'//diag(list_idx),  extinct, pcols, lchnk)
   call outfld('ABSORB'//diag(list_idx),   absorb,  pcols, lchnk)
   call outfld('AODVIS'//diag(list_idx),   aodvis,  pcols, lchnk)
   call outfld('AODABS'//diag(list_idx),   aodabs,  pcols, lchnk)
   call outfld('AODVISst'//diag(list_idx), aodvisst,pcols, lchnk)

   ! These diagnostics are output only for climate list
   if (list_idx == 0) then
      if (use_native_modal_aer_opt_helpers_impl) then
         do i = 1, ncol
            if (aodvis(i) > 1.e-10_r8) then
               ssavis(i) = ssavis(i)/aodvis(i)
            else
               ssavis(i) = 0.925_r8
            endif
         end do

         do i = 1, nnite
            ssavis(idxnite(i))     = fillvalue

            aoduv(idxnite(i))      = fillvalue
            aodnir(idxnite(i))     = fillvalue
            aoduvst(idxnite(i))    = fillvalue
            aodnirst(idxnite(i))   = fillvalue
            extinctuv(idxnite(i),:)  = fillvalue
            extinctnir(idxnite(i),:) = fillvalue

            burdendust(idxnite(i)) = fillvalue
            burdenso4(idxnite(i))  = fillvalue
            burdenpom(idxnite(i))  = fillvalue
            burdensoa(idxnite(i))  = fillvalue
            burdenbc(idxnite(i))   = fillvalue
            burdenseasalt(idxnite(i)) = fillvalue

            aodabsbc(idxnite(i))   = fillvalue

            dustaod(idxnite(i))    = fillvalue
            so4aod(idxnite(i))     = fillvalue
            pomaod(idxnite(i))     = fillvalue
            soaaod(idxnite(i))     = fillvalue
            bcaod(idxnite(i))      = fillvalue
            seasaltaod(idxnite(i)) = fillvalue
          end do
      else
         call modal_aer_opt_helpers_proof_once()
         call modal_aer_opt_sw_finalize_ssavis_codon(int(ncol, c_int64_t), &
              c_loc(aodvis(1)), c_loc(ssavis(1)))
         if (nnite > 0) then
            call modal_aer_opt_sw_climate_diag_night_codon(int(nnite, c_int64_t), &
                 int(pcols, c_int64_t), int(pver, c_int64_t), fillvalue, c_loc(idxnite(1)), &
                 c_loc(ssavis(1)), c_loc(aoduv(1)), c_loc(aodnir(1)), c_loc(aoduvst(1)), &
                 c_loc(aodnirst(1)), c_loc(extinctuv(1,1)), c_loc(extinctnir(1,1)), &
                 c_loc(burdendust(1)), c_loc(burdenso4(1)), c_loc(burdenpom(1)), &
                 c_loc(burdensoa(1)), c_loc(burdenbc(1)), c_loc(burdenseasalt(1)), &
                 c_loc(aodabsbc(1)), c_loc(dustaod(1)), c_loc(so4aod(1)), c_loc(pomaod(1)), &
                 c_loc(soaaod(1)), c_loc(bcaod(1)), c_loc(seasaltaod(1)))
         end if
      end if

      call outfld('SSAVIS',        ssavis,        pcols, lchnk)

      call outfld('EXTINCTUV',     extinctuv,     pcols, lchnk)
      call outfld('EXTINCTNIR',    extinctnir,    pcols, lchnk)
      call outfld('AODUV',         aoduv,         pcols, lchnk)
      call outfld('AODNIR',        aodnir,        pcols, lchnk)
      call outfld('AODUVst',       aoduvst,       pcols, lchnk)
      call outfld('AODNIRst',      aodnirst,      pcols, lchnk)

      call outfld('BURDENDUST',    burdendust,    pcols, lchnk)
      call outfld('BURDENSO4' ,    burdenso4,     pcols, lchnk)
      call outfld('BURDENPOM' ,    burdenpom,     pcols, lchnk)
      call outfld('BURDENSOA' ,    burdensoa,     pcols, lchnk)
      call outfld('BURDENBC'  ,    burdenbc,      pcols, lchnk)
      call outfld('BURDENSEASALT', burdenseasalt, pcols, lchnk)

      call outfld('AODABSBC',      aodabsbc,      pcols, lchnk)

      call outfld('AODDUST',       dustaod,       pcols, lchnk)
      call outfld('AODSO4',        so4aod,        pcols, lchnk)
      call outfld('AODPOM',        pomaod,        pcols, lchnk)
      call outfld('AODSOA',        soaaod,        pcols, lchnk)
      call outfld('AODBC',         bcaod,         pcols, lchnk)
      call outfld('AODSS',         seasaltaod,    pcols, lchnk)
   end if

end subroutine modal_aero_sw

!===============================================================================

subroutine modal_aero_lw(list_idx, state, pbuf, tauxar)

   ! calculates aerosol lw radiative properties

   use iso_c_binding, only: c_int64_t, c_loc

   integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
   type(physics_state), intent(in), target :: state    ! state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(out), target :: tauxar(pcols,pver,nlwbands) ! layer absorption optical depth

   ! Local variables
   integer :: i, ifld, ilw, k, l, m, nc, ns
   integer :: lchnk                    ! chunk id
   integer :: ncol                     ! number of active columns in the chunk
   integer :: nmodes
   integer :: nspec

   real(r8), pointer :: dgnumwet(:,:)  ! wet number mode diameter (m)
   real(r8), pointer :: qaerwat(:,:)   ! aerosol water (g/g)

   real(r8), pointer :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
   real(r8), pointer :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
   real(r8), pointer :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
   real(r8), pointer :: wetdens_m(:,:,:)  ! 

   real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution
   real(r8) :: alnsg_amode             ! log of geometric standard deviation of number distribution
   real(r8) :: xrad(pcols)
   real(r8), target :: cheby(ncoef,pcols,pver)  ! chebychef polynomials

   real(r8), target :: mass(pcols,pver) ! layer mass

   real(r8),    pointer :: specmmr(:,:)        ! species mass mixing ratio
   real(r8)             :: specdens            ! species density (kg/m3)
   complex(r8), pointer :: specrefindex(:)     ! species refractive index

   real(r8) :: vol(pcols)       ! volume concentration of aerosol specie (m3/kg)
   real(r8) :: dryvol(pcols)    ! volume concentration of aerosol mode (m3/kg)
   real(r8), target :: wetvol(pcols)    ! volume concentration of wet mode (m3/kg)
   real(r8) :: watervol(pcols)  ! volume concentration of water in each mode (m3/kg)
   real(r8) :: refr(pcols)      ! real part of refractive index
   real(r8) :: refi(pcols)      ! imaginary part of refractive index
   complex(r8) :: crefin(pcols) ! complex refractive index
   real(r8), pointer :: refrtablw(:,:) ! table of real refractive indices for aerosols
   real(r8), pointer :: refitablw(:,:) ! table of imag refractive indices for aerosols
   real(r8), pointer :: absplw(:,:,:,:) ! specific absorption

   integer  :: itab(pcols), jtab(pcols)
   real(r8) :: ttab(pcols), utab(pcols)
   real(r8), target :: cabs(pcols,ncoef)
   real(r8), target :: pabs(pcols)      ! parameterized specific absorption (m2/kg)
   real(r8), target :: dopaer(pcols)    ! aerosol optical depth in layer

   integer, parameter :: nerrmax_dopaer=1000
   integer  :: nerr_dopaer = 0
   real(r8) :: volf             ! volume fraction of insoluble aerosol

   character(len=*), parameter :: subname = 'modal_aero_lw'
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   call modal_aer_opt_helpers_select_impl()

   if (use_native_modal_aer_opt_helpers_impl) then
      ! initialize output variables
      tauxar(:ncol,:,:) = 0._r8

      ! dry mass in each cell
      mass(:ncol,:) = state%pdeldry(:ncol,:)*rga
   else
      call modal_aer_opt_helpers_proof_once()
      call modal_aer_opt_lw_helpers_proof_once()
      call modal_aer_opt_lw_init_state_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(nlwbands, c_int64_t), rga, c_loc(state%pdeldry(1,1)), &
           c_loc(tauxar(1,1,1)), c_loc(mass(1,1)))
   end if

   ! loop over all aerosol modes
   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx == 0) then
      ! water uptake and wet radius for the climate list has already been calculated
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet_m)
      call pbuf_get_field(pbuf, qaerwat_idx,  qaerwat_m)
   else
      ! If doing a diagnostic calculation then need to calculate the wet radius
      ! and water uptake for the diagnostic modes
      call modal_aero_calcsize_diag(state, pbuf, list_idx, dgnumdry_m)  
      call modal_aero_wateruptake_dr(state, pbuf, list_idx, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m)
   endif

   do m = 1, nmodes

      dgnumwet => dgnumwet_m(:,:,m)
      qaerwat  => qaerwat_m(:,:,m)

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtablw=refrtablw , &
         refitablw=refitablw, absplw=absplw)

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      ! calc size parameter for all columns
      call modal_lw_size_parameters(ncol, sigma_logr_aer, dgnumwet, cheby)

      do ilw = 1, nlwbands

         do k = top_lev, pver

            ! form bulk refractive index. Use volume mixing for infrared
            crefin(:ncol) = (0._r8, 0._r8)
            dryvol(:ncol) = 0._r8

            ! aerosol species loop
            do l = 1, nspec
               call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
               call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                           refindex_aer_lw=specrefindex)

               do i = 1, ncol
                  vol(i)    = specmmr(i,k)/specdens
                  dryvol(i) = dryvol(i) + vol(i)
                  crefin(i) = crefin(i) + vol(i)*specrefindex(ilw)
               end do
            end do

            do i = 1, ncol
               watervol(i) = qaerwat(i,k)/rhoh2o
               wetvol(i)   = watervol(i) + dryvol(i)
               if (watervol(i) < 0.0_r8) then
                  if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                     write(iulog,*) 'watervol,wetvol,dryvol=',watervol(i),wetvol(i),dryvol(i),' in '//subname
                  end if
                  watervol(i) = 0._r8
                  wetvol(i)   = dryvol(i)
               end if

               crefin(i) = crefin(i) + watervol(i)*crefwlw(ilw)
               if (wetvol(i) > 1.e-40_r8) crefin(i) = crefin(i)/wetvol(i)
               refr(i) = real(crefin(i))
               refi(i) = aimag(crefin(i))
            end do

            ! interpolate coefficients linear in refractive index
            ! first call calcs itab,jtab,ttab,utab
            itab(:ncol) = 0
            call binterp(absplw(:,:,:,ilw), ncol, ncoef, prefr, prefi, &
                         refr, refi, refrtablw(:,ilw), refitablw(:,ilw), &
                         itab, jtab, ttab, utab, cabs)

            ! parameterized optical properties
            if (use_native_modal_aer_opt_helpers_impl) then
               do i = 1, ncol
                  pabs(i) = 0.5_r8*cabs(i,1)
                  do nc = 2, ncoef
                     pabs(i) = pabs(i) + cheby(nc,i,k)*cabs(i,nc)
                  end do
                  pabs(i)   = pabs(i)*wetvol(i)*rhoh2o
                  pabs(i)   = max(0._r8,pabs(i))
                  dopaer(i) = pabs(i)*mass(i,k)
               end do
            else
               call modal_aer_opt_helpers_proof_once()
               call modal_aer_opt_lw_helpers_proof_once()
               call modal_aer_opt_lw_optics_props_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
                    int(pver, c_int64_t), int(k, c_int64_t), int(ncoef, c_int64_t), rhoh2o, &
                    c_loc(cheby(1,1,1)), c_loc(cabs(1,1)), c_loc(wetvol(1)), c_loc(mass(1,1)), &
                    c_loc(pabs(1)), c_loc(dopaer(1)))
            end if

            do i = 1, ncol

               if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 20._r8)) then

                  if (dopaer(i) <= -1.e-10_r8) then
                     write(iulog,*) "ERROR: Negative aerosol optical depth &
                          &in this layer."
                  else
                     write(iulog,*) "WARNING: Aerosol optical depth is &
                          &unreasonably high in this layer."
                  end if

                  write(iulog,*) 'dopaer(',i,',',k,',',m,',',lchnk,')=', dopaer(i)
                  write(iulog,*) 'k=',k,' pabs=', pabs(i)
                  write(iulog,*) 'wetvol=',wetvol(i),' dryvol=',dryvol(i),     &
                     ' watervol=',watervol(i)
                  write(iulog,*) 'cabs=', (cabs(i,l),l=1,ncoef)
                  write(iulog,*) 'crefin=', crefin(i)
                  write(iulog,*) 'nspec=', nspec
                  do l = 1,nspec
                     call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, specmmr)
                     call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                 refindex_aer_lw=specrefindex)
                     volf = specmmr(i,k)/specdens
                     write(iulog,*) 'l=',l,'vol(l)=',volf
                     write(iulog,*) 'ilw=',ilw,' specrefindex(ilw)=',specrefindex(ilw)
                     write(iulog,*) 'specdens=',specdens
                  end do

                  nerr_dopaer = nerr_dopaer + 1
                  if (nerr_dopaer >= nerrmax_dopaer .or. dopaer(i) < -1.e-10_r8) then
                     write(iulog,*) '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                     call endrun()
                  end if

               end if
            end do

            if (use_native_modal_aer_opt_helpers_impl) then
               do i = 1, ncol
                  tauxar(i,k,ilw) = tauxar(i,k,ilw) + dopaer(i)
               end do
            else
               call modal_aer_opt_helpers_proof_once()
               call modal_aer_opt_lw_helpers_proof_once()
               call modal_aer_opt_lw_accumulate_tau_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
                    int(pver, c_int64_t), int(k, c_int64_t), int(ilw, c_int64_t), c_loc(dopaer(1)), &
                    c_loc(tauxar(1,1,1)))
            end if

         end do ! k = top_lev, pver

      end do  ! nlwbands

   end do ! m = 1, nmodes

   if (list_idx > 0) then
      deallocate(dgnumdry_m)
      deallocate(dgnumwet_m)
      deallocate(qaerwat_m)
      deallocate(wetdens_m)
   end if

end subroutine modal_aero_lw

!===============================================================================
! Private routines
!===============================================================================

subroutine read_water_refindex(infilename)

   ! read water refractive index file and set module data

   character*(*), intent(in) :: infilename   ! modal optics filename

   ! Local variables

   integer            :: i, ierr
   type(file_desc_t)  :: ncid              ! pio file handle
   integer            :: did               ! dimension ids
   integer            :: dimlen            ! dimension lengths
   type(var_desc_t)   :: vid               ! variable ids
   real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
   real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared
   !----------------------------------------------------------------------------

   ! open file
   call cam_pio_openfile(ncid, infilename, PIO_NOWRITE)

   ! inquire dimensions.  Check that file values match parameter values.

   ierr = pio_inq_dimid(ncid, 'lw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nlwbands) then
      write(iulog,*) 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
      call endrun('read_modal_optics: bad lw_band value')
   endif

   ierr = pio_inq_dimid(ncid, 'sw_band', did)
   ierr = pio_inq_dimlen(ncid, did, dimlen)
   if (dimlen .ne. nswbands) then
      write(iulog,*) 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
      call endrun('read_modal_optics: bad sw_band value')
   endif

   ! read variables
   ierr = pio_inq_varid(ncid, 'refindex_real_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refrwsw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_sw', vid)
   ierr = pio_get_var(ncid, vid, refiwsw)

   ierr = pio_inq_varid(ncid, 'refindex_real_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refrwlw)

   ierr = pio_inq_varid(ncid, 'refindex_im_water_lw', vid)
   ierr = pio_get_var(ncid, vid, refiwlw)

   ! set complex representation of refractive indices as module data
   do i = 1, nswbands
      crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)),kind=r8)
   end do
   do i = 1, nlwbands
      crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)),kind=r8)
   end do

   call pio_closefile(ncid)

end subroutine read_water_refindex

!===============================================================================

subroutine modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

   use iso_c_binding, only: c_int64_t, c_loc

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
   real(r8), intent(in), target, contiguous  :: dgnumwet(:,:)   ! aerosol wet number mode diameter (m)
   real(r8), intent(out), target, contiguous :: radsurf(:,:)    ! aerosol surface mode radius
   real(r8), intent(out), target, contiguous :: logradsurf(:,:) ! log(aerosol surface mode radius)
   real(r8), intent(out), target, contiguous :: cheb(:,:,:)

   call modal_aer_opt_helpers_select_impl()

   if (use_native_modal_aer_opt_helpers_impl) then
      call modal_size_parameters_native(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)
      return
   end if

   call modal_aer_opt_helpers_proof_once()
   call modal_aer_opt_size_parameters_codon(int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(top_lev, c_int64_t), int(ncol, c_int64_t), int(ncoef, c_int64_t), &
        sigma_logr_aer, xrmin, xrmax, c_loc(dgnumwet(1,1)), c_loc(radsurf(1,1)), &
        c_loc(logradsurf(1,1)), c_loc(cheb(1,1,1)))

end subroutine modal_size_parameters

!===============================================================================

subroutine modal_size_parameters_native(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
   real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter (m)
   real(r8), intent(out) :: radsurf(:,:)    ! aerosol surface mode radius
   real(r8), intent(out) :: logradsurf(:,:) ! log(aerosol surface mode radius)
   real(r8), intent(out) :: cheb(:,:,:)

   integer  :: i, k, nc
   real(r8) :: alnsg_amode
   real(r8) :: explnsigma
   real(r8) :: xrad(pcols) ! normalized aerosol radius
   !-------------------------------------------------------------------------------

   alnsg_amode = log(sigma_logr_aer)
   explnsigma = exp(2.0_r8*alnsg_amode*alnsg_amode)

   do k = top_lev, pver
      do i = 1, ncol
         ! convert from number mode diameter to surface area
         radsurf(i,k) = 0.5_r8*dgnumwet(i,k)*explnsigma
         logradsurf(i,k) = log(radsurf(i,k))
         ! normalize size parameter
         xrad(i) = max(logradsurf(i,k),xrmin)
         xrad(i) = min(xrad(i),xrmax)
         xrad(i) = (2._r8*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
         ! chebyshev polynomials
         cheb(1,i,k) = 1._r8
         cheb(2,i,k) = xrad(i)
         do nc = 3, ncoef
            cheb(nc,i,k) = 2._r8*xrad(i)*cheb(nc-1,i,k)-cheb(nc-2,i,k)
         end do
      end do
   end do

end subroutine modal_size_parameters_native

!===============================================================================

subroutine modal_lw_size_parameters(ncol, sigma_logr_aer, dgnumwet, cheby)

   use iso_c_binding, only: c_int64_t, c_loc

   integer,  intent(in) :: ncol
   real(r8), intent(in) :: sigma_logr_aer  ! geometric standard deviation of number distribution
   real(r8), intent(in), target, contiguous  :: dgnumwet(:,:) ! aerosol wet number mode diameter (m)
   real(r8), intent(out), target, contiguous :: cheby(:,:,:)

   call modal_aer_opt_helpers_select_impl()

   if (use_native_modal_aer_opt_helpers_impl) then
      call modal_lw_size_parameters_native(ncol, sigma_logr_aer, dgnumwet, cheby)
      return
   end if

   call modal_aer_opt_helpers_proof_once()
   call modal_aer_opt_lw_size_parameters_codon(int(pcols, c_int64_t), int(pver, c_int64_t), &
        int(top_lev, c_int64_t), int(ncol, c_int64_t), int(ncoef, c_int64_t), &
        sigma_logr_aer, xrmin, xrmax, c_loc(dgnumwet(1,1)), c_loc(cheby(1,1,1)))

end subroutine modal_lw_size_parameters

!===============================================================================

subroutine modal_lw_size_parameters_native(ncol, sigma_logr_aer, dgnumwet, cheby)

   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
   real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter (m)
   real(r8), intent(out) :: cheby(:,:,:)

   integer  :: i, k, nc
   real(r8) :: alnsg_amode
   real(r8) :: xrad(pcols)
   !-------------------------------------------------------------------------------

   ! This is the same calculation that's done in modal_size_parameters, but some
   ! intermediate results are saved and the chebyshev polynomials are stored in an
   ! array with different index order.
   do k = top_lev, pver
      do i = 1, ncol
         alnsg_amode = log( sigma_logr_aer )
         ! convert from number diameter to surface area
         xrad(i) = log(0.5_r8*dgnumwet(i,k)) + 2.0_r8*alnsg_amode*alnsg_amode
         ! normalize size parameter
         xrad(i) = max(xrad(i), xrmin)
         xrad(i) = min(xrad(i), xrmax)
         xrad(i) = (2*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
         ! chebyshev polynomials
         cheby(1,i,k) = 1.0_r8
         cheby(2,i,k) = xrad(i)
         do nc = 3, ncoef
            cheby(nc,i,k) = 2.0_r8*xrad(i)*cheby(nc-1,i,k)-cheby(nc-2,i,k)
         end do
      end do
   end do

end subroutine modal_lw_size_parameters_native

!===============================================================================

      subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

      use iso_c_binding, only: c_int64_t, c_loc

!     bilinear interpolation of table
!
      implicit none
      integer, intent(in) :: im,jm,km,ncol
      real(r8), intent(in), target, contiguous :: table(:,:,:),xtab(:),ytab(:)
      real(r8), intent(in), target, contiguous :: x(:),y(:)
      integer, intent(inout), target, contiguous :: ix(:),jy(:)
      real(r8), intent(inout), target, contiguous :: t(:),u(:)
      real(r8), intent(out), target, contiguous :: out(:,:)

      call modal_aer_opt_helpers_select_impl()

      if (use_native_modal_aer_opt_helpers_impl) then
         call binterp_native(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)
         return
      end if

      call modal_aer_opt_helpers_proof_once()
      call modal_aer_opt_binterp_codon(int(pcols, c_int64_t), int(ncol, c_int64_t), &
           int(km, c_int64_t), int(im, c_int64_t), int(jm, c_int64_t), c_loc(table(1,1,1)), &
           c_loc(x(1)), c_loc(y(1)), c_loc(xtab(1)), c_loc(ytab(1)), c_loc(ix(1)), c_loc(jy(1)), &
           c_loc(t(1)), c_loc(u(1)), c_loc(out(1,1)))

      end subroutine binterp

!===============================================================================

      subroutine binterp_native(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

!     bilinear interpolation of table
!
      implicit none
      integer im,jm,km,ncol
      real(r8) table(km,im,jm),xtab(im),ytab(jm),out(pcols,km)
      integer i,ix(pcols),ip1,j,jy(pcols),jp1,k,ic
      real(r8) x(pcols),dx,t(pcols),y(pcols),dy,u(pcols), &
             tu(pcols),tuc(pcols),tcu(pcols),tcuc(pcols)

      if(ix(1).gt.0)go to 30
      if(im.gt.1)then
        do ic=1,ncol
          do i=1,im
            if(x(ic).lt.xtab(i))go to 10
          enddo
   10     ix(ic)=max0(i-1,1)
          ip1=min(ix(ic)+1,im)
          dx=(xtab(ip1)-xtab(ix(ic)))
          if(abs(dx).gt.1.e-20_r8)then
             t(ic)=(x(ic)-xtab(ix(ic)))/dx
          else
             t(ic)=0._r8
          endif
	end do
      else
        ix(:ncol)=1
        t(:ncol)=0._r8
      endif
      if(jm.gt.1)then
        do ic=1,ncol
          do j=1,jm
            if(y(ic).lt.ytab(j))go to 20
          enddo
   20     jy(ic)=max0(j-1,1)
          jp1=min(jy(ic)+1,jm)
          dy=(ytab(jp1)-ytab(jy(ic)))
          if(abs(dy).gt.1.e-20_r8)then
             u(ic)=(y(ic)-ytab(jy(ic)))/dy
             if(u(ic).lt.0._r8.or.u(ic).gt.1._r8)then
                write(iulog,*) 'u,y,jy,ytab,dy=',u(ic),y(ic),jy(ic),ytab(jy(ic)),dy
             endif
          else
            u(ic)=0._r8
          endif
	end do
      else
        jy(:ncol)=1
        u(:ncol)=0._r8
      endif
   30 continue
      do ic=1,ncol
         tu(ic)=t(ic)*u(ic)
         tuc(ic)=t(ic)-tu(ic)
         tcuc(ic)=1._r8-tuc(ic)-u(ic)
         tcu(ic)=u(ic)-tu(ic)
         jp1=min(jy(ic)+1,jm)
         ip1=min(ix(ic)+1,im)
         do k=1,km
            out(ic,k)=tcuc(ic)*table(k,ix(ic),jy(ic))+tuc(ic)*table(k,ip1,jy(ic))   &
               +tu(ic)*table(k,ip1,jp1)+tcu(ic)*table(k,ix(ic),jp1)
	 end do
      enddo
      return
      end subroutine binterp_native

end module modal_aer_opt
