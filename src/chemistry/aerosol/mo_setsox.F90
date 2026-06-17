
module MO_SETSOX

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use iso_c_binding, only : c_int64_t

  private
  public :: sox_inti, setsox, setsox_shell_codon_wrap
  public :: has_sox

  save
  logical            ::  inv_o3
  integer            ::  id_msa

  integer :: id_so2, id_nh3, id_hno3, id_h2o2, id_o3, id_ho2
  integer :: id_so4, id_h2so4

  logical :: has_sox = .true.
  logical :: inv_so2, inv_nh3, inv_hno3, inv_h2o2, inv_ox, inv_nh4no3, inv_ho2

  logical :: cloud_borne = .false.
  logical :: modal_aerosols = .false.
  logical :: setsox_init_fields_use_native_impl = .false.
  logical :: setsox_init_fields_impl_selected = .false.
  logical :: setsox_init_fields_proof_written = .false.
  logical :: setsox_init_fields_wrap_proof_written = .false.
  logical :: setsox_ph_solve_use_native_impl = .false.
  logical :: setsox_ph_solve_impl_selected = .false.
  logical :: setsox_ph_solve_proof_written = .false.
  logical :: setsox_ph_solve_wrap_proof_written = .false.
  logical :: setsox_aqchem_predict_use_native_impl = .false.
  logical :: setsox_aqchem_predict_impl_selected = .false.
  logical :: setsox_aqchem_predict_proof_written = .false.
  logical :: setsox_aqchem_predict_wrap_proof_written = .false.
  logical :: setsox_xph_lwc_diag_use_native_impl = .false.
  logical :: setsox_xph_lwc_diag_impl_selected = .false.
  logical :: setsox_xph_lwc_diag_proof_written = .false.
  logical :: setsox_xph_lwc_diag_wrap_proof_written = .false.
  logical :: setsox_shell_wrap_proof_written = .false.
  logical :: setsox_shell_finalize_wrap_proof_written = .false.
  logical :: sox_inti_proof_written = .false.

  interface
     function sox_inti_active_codon(active_c) result(out_c) bind(c, name="sox_inti_active_codon")
       use iso_c_binding, only : c_int64_t
       integer(c_int64_t), value :: active_c
       integer(c_int64_t) :: out_c
     end function sox_inti_active_codon
  end interface

contains

  subroutine setsox_append_impl_proof(proof_line)

    implicit none

    character(len=*), intent(in) :: proof_line

    character(len=512) :: proof_path
    integer :: status, n, unit_id

    call get_environment_variable('SOX_CHEM_PROOF_FILE', value=proof_path, length=n, status=status)
    if (status /= 0 .or. n <= 0) return

    open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
         position='append', iostat=status)
    if (status /= 0) return

    write(unit_id,'(A)') trim(proof_line)
    close(unit_id)

  end subroutine setsox_append_impl_proof

  subroutine setsox_select_impl(env_name, helper_name, use_native_impl, impl_selected, proof_written)

    use spmd_utils, only : masterproc

    implicit none

    character(len=*), intent(in) :: env_name, helper_name
    logical, intent(inout) :: use_native_impl, impl_selected, proof_written

    character(len=48) :: impl_name
    character(len=160) :: proof_line
    integer :: status, n, i, code

    if (impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl(env_name, impl_name, n, status)

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
          proof_line = trim(helper_name)//' selector entered implementation = native'
       else
          proof_line = trim(helper_name)//' selector entered implementation = codon'
       end if
       write(iulog,'(A)') trim(proof_line)
       if (.not. proof_written) then
          call setsox_append_impl_proof(trim(proof_line))
          proof_written = .true.
       end if
       call flush(iulog)
    end if

  end subroutine setsox_select_impl

  subroutine setsox_init_fields_select_impl()

    call setsox_select_impl('SETSOX_INIT_FIELDS_IMPL', 'setsox_init_fields', &
         setsox_init_fields_use_native_impl, setsox_init_fields_impl_selected, &
         setsox_init_fields_proof_written)

  end subroutine setsox_init_fields_select_impl

  subroutine setsox_ph_solve_select_impl()

    call setsox_select_impl('SETSOX_PH_SOLVE_IMPL', 'setsox_ph_solve', &
         setsox_ph_solve_use_native_impl, setsox_ph_solve_impl_selected, &
         setsox_ph_solve_proof_written)

  end subroutine setsox_ph_solve_select_impl

  subroutine setsox_aqchem_predict_select_impl()

    call setsox_select_impl('SETSOX_AQCHEM_PREDICT_IMPL', 'setsox_aqchem_predict', &
         setsox_aqchem_predict_use_native_impl, setsox_aqchem_predict_impl_selected, &
         setsox_aqchem_predict_proof_written)

  end subroutine setsox_aqchem_predict_select_impl

  subroutine setsox_xph_lwc_diag_select_impl()

    call setsox_select_impl('SETSOX_XPH_LWC_DIAG_IMPL', 'setsox_xph_lwc_diag', &
         setsox_xph_lwc_diag_use_native_impl, setsox_xph_lwc_diag_impl_selected, &
         setsox_xph_lwc_diag_proof_written)

  end subroutine setsox_xph_lwc_diag_select_impl

  subroutine setsox_init_fields_codon_wrap(stage, ncol, ph0_in, xhnm, invariants, qin, cfact, xph, &
       xso2, xhno3, xh2o2, xnh3, xo3, xho2, xh2so4, xso4, xno3, xnh4, xmsa)

    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    use ppgrid,        only : pcols, pver
    use chem_mods,     only : gas_pcnst, nfs
    use spmd_utils,    only : masterproc

    implicit none

    integer, intent(in) :: stage, ncol
    real(r8), intent(in) :: ph0_in
    real(r8), target, intent(in) :: xhnm(ncol,pver), invariants(ncol,pver,nfs)
    real(r8), target, intent(in) :: qin(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: cfact(ncol,pver), xph(ncol,pver)
    real(r8), target, intent(inout) :: xso2(ncol,pver), xhno3(ncol,pver), xh2o2(ncol,pver)
    real(r8), target, intent(inout) :: xnh3(ncol,pver), xo3(ncol,pver), xho2(ncol,pver)
    real(r8), target, intent(inout) :: xh2so4(ncol,pver), xso4(ncol,pver), xno3(ncol,pver)
    real(r8), target, intent(inout) :: xnh4(ncol,pver), xmsa(ncol,pver)

    integer(c_int64_t) :: cloud_borne_c, inv_so2_c, inv_h2o2_c, inv_o3_c, inv_ho2_c
    character(len=96) :: proof_line

    interface
       subroutine setsox_init_fields_codon(stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, nfs_c, &
            cloud_borne_c, inv_so2_c, inv_h2o2_c, inv_o3_c, inv_ho2_c, id_so2_c, id_hno3_c, &
            id_h2o2_c, id_nh3_c, id_o3_c, id_ho2_c, id_h2so4_c, id_so4_c, id_msa_c, ph0_c, &
            xhnm_p, invariants_p, qin_p, cfact_p, xph_p, xso2_p, xhno3_p, xh2o2_p, xnh3_p, &
            xo3_p, xho2_p, xh2so4_p, xso4_p, xno3_p, xnh4_p, xmsa_p) &
            bind(c, name="setsox_init_fields_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, nfs_c
         integer(c_int64_t), value :: cloud_borne_c, inv_so2_c, inv_h2o2_c, inv_o3_c, inv_ho2_c
         integer(c_int64_t), value :: id_so2_c, id_hno3_c, id_h2o2_c, id_nh3_c, id_o3_c, id_ho2_c
         integer(c_int64_t), value :: id_h2so4_c, id_so4_c, id_msa_c
         real(c_double), value :: ph0_c
         type(c_ptr), value :: xhnm_p, invariants_p, qin_p, cfact_p, xph_p, xso2_p, xhno3_p
         type(c_ptr), value :: xh2o2_p, xnh3_p, xo3_p, xho2_p, xh2so4_p, xso4_p, xno3_p
         type(c_ptr), value :: xnh4_p, xmsa_p
       end subroutine setsox_init_fields_codon
    end interface

    if (masterproc .and. .not. setsox_init_fields_wrap_proof_written) then
       write(proof_line,'(A,I0)') 'setsox_init_fields_codon_wrap entered stage=', stage
       write(iulog,'(A)') trim(proof_line)
       call setsox_append_impl_proof(trim(proof_line))
       setsox_init_fields_wrap_proof_written = .true.
       call flush(iulog)
    end if

    cloud_borne_c = 0_c_int64_t
    inv_so2_c = 0_c_int64_t
    inv_h2o2_c = 0_c_int64_t
    inv_o3_c = 0_c_int64_t
    inv_ho2_c = 0_c_int64_t
    if (cloud_borne) cloud_borne_c = 1_c_int64_t
    if (inv_so2) inv_so2_c = 1_c_int64_t
    if (inv_h2o2) inv_h2o2_c = 1_c_int64_t
    if (inv_o3) inv_o3_c = 1_c_int64_t
    if (inv_ho2) inv_ho2_c = 1_c_int64_t

    call setsox_init_fields_codon( &
         int(stage, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(gas_pcnst, c_int64_t), int(nfs, c_int64_t), cloud_borne_c, inv_so2_c, inv_h2o2_c, &
         inv_o3_c, inv_ho2_c, int(id_so2, c_int64_t), int(id_hno3, c_int64_t), &
         int(id_h2o2, c_int64_t), int(id_nh3, c_int64_t), int(id_o3, c_int64_t), &
         int(id_ho2, c_int64_t), int(id_h2so4, c_int64_t), int(id_so4, c_int64_t), &
         int(id_msa, c_int64_t), real(ph0_in, c_double), c_loc(xhnm), c_loc(invariants), c_loc(qin), &
         c_loc(cfact), c_loc(xph), c_loc(xso2), c_loc(xhno3), c_loc(xh2o2), c_loc(xnh3), &
         c_loc(xo3), c_loc(xho2), c_loc(xh2so4), c_loc(xso4), c_loc(xno3), c_loc(xnh4), c_loc(xmsa) )

  end subroutine setsox_init_fields_codon_wrap

  subroutine setsox_ph_solve_codon_wrap(ncol, itermax_in, const0_in, ra_in, xkw_in, press, tfld, &
       cldfrc, xhnm, xlwc, xso4c, xnh4c, xno3c, xso4, xnh4, xno3, xso2, xhno3, xnh3, xph, so4_fact)

    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    use ppgrid,        only : pcols, pver
    use spmd_utils,    only : masterproc

    implicit none

    integer, intent(in) :: ncol, itermax_in
    real(r8), intent(in) :: const0_in, ra_in, xkw_in
    real(r8), target, intent(in) :: press(:,:), tfld(:,:), cldfrc(:,:), xhnm(ncol,pver)
    real(r8), target, intent(in) :: xlwc(:,:), xso4c(:,:), xnh4c(:,:), xno3c(:,:)
    real(r8), target, intent(inout) :: xso4(ncol,pver), xnh4(ncol,pver), xno3(ncol,pver)
    real(r8), target, intent(in) :: xso2(ncol,pver), xhno3(ncol,pver), xnh3(ncol,pver)
    real(r8), target, intent(inout) :: xph(ncol,pver)
    real(r8), intent(in) :: so4_fact

    integer(c_int64_t) :: cloud_borne_c
    character(len=96) :: proof_line

    interface
       subroutine setsox_ph_solve_codon(ncol_c, pcols_c, pver_c, itermax_c, cloud_borne_c, const0_c, &
            ra_c, xkw_c, so4_fact_c, press_p, tfld_p, cldfrc_p, xhnm_p, xlwc_p, xso4c_p, xnh4c_p, &
            xno3c_p, xso4_p, xnh4_p, xno3_p, xso2_p, xhno3_p, xnh3_p, xph_p) &
            bind(c, name="setsox_ph_solve_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, itermax_c, cloud_borne_c
         real(c_double), value :: const0_c, ra_c, xkw_c, so4_fact_c
         type(c_ptr), value :: press_p, tfld_p, cldfrc_p, xhnm_p, xlwc_p, xso4c_p, xnh4c_p
         type(c_ptr), value :: xno3c_p, xso4_p, xnh4_p, xno3_p, xso2_p, xhno3_p, xnh3_p, xph_p
       end subroutine setsox_ph_solve_codon
    end interface

    if (masterproc .and. .not. setsox_ph_solve_wrap_proof_written) then
       proof_line = 'setsox_ph_solve_codon_wrap entered'
       write(iulog,'(A)') trim(proof_line)
       call setsox_append_impl_proof(trim(proof_line))
       setsox_ph_solve_wrap_proof_written = .true.
       call flush(iulog)
    end if

    cloud_borne_c = 0_c_int64_t
    if (cloud_borne) cloud_borne_c = 1_c_int64_t

    call setsox_ph_solve_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(itermax_in, c_int64_t), &
         cloud_borne_c, real(const0_in, c_double), real(ra_in, c_double), real(xkw_in, c_double), &
         real(so4_fact, c_double), c_loc(press), c_loc(tfld), c_loc(cldfrc), c_loc(xhnm), c_loc(xlwc), &
         c_loc(xso4c), c_loc(xnh4c), c_loc(xno3c), c_loc(xso4), c_loc(xnh4), c_loc(xno3), &
         c_loc(xso2), c_loc(xhno3), c_loc(xnh3), c_loc(xph) )

  end subroutine setsox_ph_solve_codon_wrap

  subroutine setsox_aqchem_predict_codon_wrap(ncol, dtime, const0_in, kh0_in, kh1_in, kh2_in, kh3_in, &
       ra_in, xkw_in, press, tfld, xhnm, xlwc, xph, xho2, &
       xhno3, xno3, xh2o2, xso2, xo3, xnh3, xnh4, xso4, xso4_init, xdelso4hp, &
       hno3g, nh3g, hehno3, heh2o2, heso2, henh3, heo3)

    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    use ppgrid,        only : pcols, pver
    use spmd_utils,    only : masterproc

    implicit none

    integer, intent(in) :: ncol
    real(r8), intent(in) :: dtime, const0_in, kh0_in, kh1_in, kh2_in, kh3_in, ra_in, xkw_in
    real(r8), target, intent(in) :: press(:,:), tfld(:,:), xhnm(ncol,pver), xlwc(:,:)
    real(r8), target, intent(in) :: xph(ncol,pver), xho2(ncol,pver), xhno3(ncol,pver), xno3(ncol,pver)
    real(r8), target, intent(inout) :: xh2o2(ncol,pver), xso2(ncol,pver), xso4(ncol,pver)
    real(r8), target, intent(in) :: xo3(ncol,pver), xnh3(ncol,pver), xnh4(ncol,pver)
    real(r8), target, intent(inout) :: xso4_init(ncol,pver), xdelso4hp(ncol,pver)
    real(r8), target, intent(inout) :: hno3g(ncol,pver), nh3g(ncol,pver)
    real(r8), target, intent(inout) :: hehno3(ncol,pver), heh2o2(ncol,pver), heso2(ncol,pver)
    real(r8), target, intent(inout) :: henh3(ncol,pver), heo3(ncol,pver)

    integer(c_int64_t) :: cloud_borne_c, modal_aerosols_c
    character(len=96) :: proof_line

    interface
       subroutine setsox_aqchem_predict_codon(ncol_c, pcols_c, pver_c, cloud_borne_c, modal_aerosols_c, &
            id_nh3_c, dtime_c, const0_c, kh0_c, kh1_c, kh2_c, kh3_c, ra_c, xkw_c, press_p, tfld_p, &
            xhnm_p, xlwc_p, xph_p, xho2_p, xhno3_p, xno3_p, xh2o2_p, xso2_p, xo3_p, xnh3_p, &
            xnh4_p, xso4_p, xso4_init_p, xdelso4hp_p, hno3g_p, nh3g_p, hehno3_p, heh2o2_p, &
            heso2_p, henh3_p, heo3_p) bind(c, name="setsox_aqchem_predict_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, cloud_borne_c, modal_aerosols_c, id_nh3_c
         real(c_double), value :: dtime_c, const0_c, kh0_c, kh1_c, kh2_c, kh3_c, ra_c, xkw_c
         type(c_ptr), value :: press_p, tfld_p, xhnm_p, xlwc_p, xph_p, xho2_p, xhno3_p, xno3_p
         type(c_ptr), value :: xh2o2_p, xso2_p, xo3_p, xnh3_p, xnh4_p, xso4_p, xso4_init_p
         type(c_ptr), value :: xdelso4hp_p, hno3g_p, nh3g_p, hehno3_p, heh2o2_p, heso2_p
         type(c_ptr), value :: henh3_p, heo3_p
       end subroutine setsox_aqchem_predict_codon
    end interface

    if (masterproc .and. .not. setsox_aqchem_predict_wrap_proof_written) then
       proof_line = 'setsox_aqchem_predict_codon_wrap entered'
       write(iulog,'(A)') trim(proof_line)
       call setsox_append_impl_proof(trim(proof_line))
       setsox_aqchem_predict_wrap_proof_written = .true.
       call flush(iulog)
    end if

    cloud_borne_c = 0_c_int64_t
    modal_aerosols_c = 0_c_int64_t
    if (cloud_borne) cloud_borne_c = 1_c_int64_t
    if (modal_aerosols) modal_aerosols_c = 1_c_int64_t

    call setsox_aqchem_predict_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), cloud_borne_c, &
         modal_aerosols_c, int(id_nh3, c_int64_t), real(dtime, c_double), real(const0_in, c_double), &
         real(kh0_in, c_double), real(kh1_in, c_double), real(kh2_in, c_double), real(kh3_in, c_double), &
         real(ra_in, c_double), real(xkw_in, c_double), c_loc(press), c_loc(tfld), c_loc(xhnm), c_loc(xlwc), &
         c_loc(xph), c_loc(xho2), c_loc(xhno3), c_loc(xno3), c_loc(xh2o2), c_loc(xso2), c_loc(xo3), &
         c_loc(xnh3), c_loc(xnh4), c_loc(xso4), c_loc(xso4_init), c_loc(xdelso4hp), c_loc(hno3g), &
         c_loc(nh3g), c_loc(hehno3), c_loc(heh2o2), c_loc(heso2), c_loc(henh3), c_loc(heo3) )

  end subroutine setsox_aqchem_predict_codon_wrap

  subroutine setsox_xph_lwc_diag_codon_wrap(ncol, cldfrc, lwc, xph, xphlwc)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use ppgrid,        only : pcols, pver
    use spmd_utils,    only : masterproc

    implicit none

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: cldfrc(:,:), lwc(ncol,pver), xph(ncol,pver)
    real(r8), target, intent(inout) :: xphlwc(ncol,pver)

    character(len=96) :: proof_line

    interface
       subroutine setsox_xph_lwc_diag_codon(ncol_c, pcols_c, pver_c, cldfrc_p, lwc_p, xph_p, xphlwc_p) &
            bind(c, name="setsox_xph_lwc_diag_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: cldfrc_p, lwc_p, xph_p, xphlwc_p
       end subroutine setsox_xph_lwc_diag_codon
    end interface

    if (masterproc .and. .not. setsox_xph_lwc_diag_wrap_proof_written) then
       proof_line = 'setsox_xph_lwc_diag_codon_wrap entered'
       write(iulog,'(A)') trim(proof_line)
       call setsox_append_impl_proof(trim(proof_line))
       setsox_xph_lwc_diag_wrap_proof_written = .true.
       call flush(iulog)
    end if

    call setsox_xph_lwc_diag_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         c_loc(cldfrc), c_loc(lwc), c_loc(xph), c_loc(xphlwc) )

  end subroutine setsox_xph_lwc_diag_codon_wrap

  subroutine setsox_shell_codon_wrap(ncol, lchnk, loffset, dtime, press, pdel, tfld, mbar, lwc, cldfrc, &
       cldnum, xhnm, invariants, qcw, qin)

    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    use ppgrid,        only : pcols, pver
    use chem_mods,     only : gas_pcnst, nfs, adv_mass
    use modal_aero_data, only : ntot_amode, modeptr_accum, numptrcw_amode, lptr_so4_cw_amode, &
         lptr_msa_cw_amode, lptr_nh4_cw_amode, cnst_name_cw, specmw_so4_amode
    use cam_history,   only : outfld
    use physconst,     only : gravit
    use spmd_utils,    only : masterproc
    use sox_cldaero_mod, only : sox_cldaero_create_obj, sox_cldaero_destroy_obj
    use cldaero_mod,   only : cldaero_conc_t
    use mo_constants,  only : pi

    implicit none

    integer, intent(in) :: ncol, lchnk, loffset
    real(r8), intent(in) :: dtime
    real(r8), target, intent(in) :: press(:,:), pdel(:,:), tfld(:,:), mbar(:,:)
    real(r8), target, intent(in) :: lwc(ncol,pver), cldfrc(:,:), cldnum(:,:), xhnm(ncol,pver)
    real(r8), target, intent(in) :: invariants(ncol,pver,nfs)
    real(r8), target, intent(inout) :: qcw(ncol,pver,gas_pcnst), qin(ncol,pver,gas_pcnst)

    integer,  parameter :: itermax = 20
    real(r8), parameter :: ph0 = 5.0_r8
    real(r8), parameter :: const0 = 1.e3_r8/6.023e23_r8
    real(r8), parameter :: kh0 = 9.e3_r8
    real(r8), parameter :: kh1 = 2.05e-5_r8
    real(r8), parameter :: kh2 = 8.6e5_r8
    real(r8), parameter :: kh3 = 1.e8_r8
    real(r8), parameter :: ra = 8314._r8/101325._r8
    real(r8), parameter :: xkw = 1.e-14_r8

    real(r8), target :: xdelso4hp(ncol,pver), xphlwc(ncol,pver)
    real(r8), target :: hno3g(ncol,pver), nh3g(ncol,pver)
    real(r8), target :: xhno3(ncol,pver), xh2o2(ncol,pver), xso2(ncol,pver), xso4(ncol,pver), xno3(ncol,pver)
    real(r8), target :: xnh3(ncol,pver), xnh4(ncol,pver), xo3(ncol,pver), cfact(ncol,pver)
    real(r8), target :: xph(ncol,pver), xho2(ncol,pver), xh2so4(ncol,pver), xmsa(ncol,pver)
    real(r8), target :: xso4_init(ncol,pver), hehno3(ncol,pver), heh2o2(ncol,pver), heso2(ncol,pver)
    real(r8), target :: henh3(ncol,pver), heo3(ncol,pver)
    real(r8), target :: dqdt_aqso4(ncol,pver,gas_pcnst), dqdt_aqh2so4(ncol,pver,gas_pcnst)
    real(r8), target :: dqdt_aqhprxn(ncol,pver), dqdt_aqo3rxn(ncol,pver)
    real(r8), target :: faqgain_msa(ntot_amode), faqgain_so4(ntot_amode), qnum_c(ntot_amode)
    real(r8), target :: sflx_aqso4(ncol,ntot_amode), sflx_aqh2so4(ncol,ntot_amode)
    real(r8), target :: sflx_aqhprxn(ncol), sflx_aqo3rxn(ncol)
    real(r8), target :: adv_mass_work(gas_pcnst)
    real(r8), pointer :: xso4c(:,:), xnh4c(:,:), xno3c(:,:)
    type(cldaero_conc_t), pointer :: cldconc
    integer(c_int64_t), target :: numptrcw_amode_c(ntot_amode), lptr_so4_cw_amode_c(ntot_amode)
    integer(c_int64_t), target :: lptr_msa_cw_amode_c(ntot_amode), lptr_nh4_cw_amode_c(ntot_amode)
    integer(c_int64_t) :: cloud_borne_c, modal_aerosols_c, inv_so2_c, inv_h2o2_c, inv_o3_c, inv_ho2_c
    character(len=160) :: proof_line
    integer :: l, m, n

    interface
       subroutine setsox_shell_stage_dispatch_codon(stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, nfs_c, ntot_amode_c, &
            loffset_c, itermax_c, cloud_borne_c, modal_aerosols_c, inv_so2_c, inv_h2o2_c, inv_o3_c, &
            inv_ho2_c, id_so2_c, id_hno3_c, id_h2o2_c, id_nh3_c, id_o3_c, id_ho2_c, id_h2so4_c, &
            id_so4_c, id_msa_c, modeptr_accum_c, dtime_c, ph0_c, const0_c, kh0_c, kh1_c, kh2_c, kh3_c, ra_c, &
            xkw_c, so4_fact_c, pi_c, xhnm_p, invariants_p, qin_p, cfact_p, xph_p, xso2_p, xhno3_p, &
            xh2o2_p, xnh3_p, xo3_p, xho2_p, xh2so4_p, xso4_p, xno3_p, xnh4_p, xmsa_p, press_p, tfld_p, &
            cldfrc_p, cldnum_p, lwc_p, xlwc_p, xso4c_p, xnh4c_p, xno3c_p, xso4_init_p, xdelso4hp_p, &
            hno3g_p, nh3g_p, hehno3_p, heh2o2_p, heso2_p, henh3_p, heo3_p, xphlwc_p, qcw_p, &
            dqdt_aqso4_p, dqdt_aqh2so4_p, dqdt_aqhprxn_p, dqdt_aqo3rxn_p, faqgain_msa_p, faqgain_so4_p, &
            qnum_c_p, numptrcw_amode_p, lptr_so4_cw_amode_p, lptr_msa_cw_amode_p, lptr_nh4_cw_amode_p) &
            bind(c, name="setsox_shell_stage_dispatch_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, gas_pcnst_c, nfs_c, ntot_amode_c
         integer(c_int64_t), value :: loffset_c, itermax_c, cloud_borne_c, modal_aerosols_c
         integer(c_int64_t), value :: inv_so2_c, inv_h2o2_c, inv_o3_c, inv_ho2_c
         integer(c_int64_t), value :: id_so2_c, id_hno3_c, id_h2o2_c, id_nh3_c, id_o3_c, id_ho2_c
         integer(c_int64_t), value :: id_h2so4_c, id_so4_c, id_msa_c, modeptr_accum_c
         real(c_double), value :: dtime_c, ph0_c, const0_c, kh0_c, kh1_c, kh2_c, kh3_c, ra_c, xkw_c
         real(c_double), value :: so4_fact_c, pi_c
         type(c_ptr), value :: xhnm_p, invariants_p, qin_p, cfact_p, xph_p, xso2_p, xhno3_p
         type(c_ptr), value :: xh2o2_p, xnh3_p, xo3_p, xho2_p, xh2so4_p, xso4_p, xno3_p, xnh4_p
         type(c_ptr), value :: xmsa_p, press_p, tfld_p, cldfrc_p, cldnum_p, lwc_p, xlwc_p, xso4c_p
         type(c_ptr), value :: xnh4c_p, xno3c_p, xso4_init_p, xdelso4hp_p, hno3g_p, nh3g_p
         type(c_ptr), value :: hehno3_p, heh2o2_p, heso2_p, henh3_p, heo3_p, xphlwc_p, qcw_p
         type(c_ptr), value :: dqdt_aqso4_p, dqdt_aqh2so4_p, dqdt_aqhprxn_p, dqdt_aqo3rxn_p
         type(c_ptr), value :: faqgain_msa_p, faqgain_so4_p, qnum_c_p, numptrcw_amode_p
         type(c_ptr), value :: lptr_so4_cw_amode_p, lptr_msa_cw_amode_p, lptr_nh4_cw_amode_p
       end subroutine setsox_shell_stage_dispatch_codon

       subroutine setsox_shell_finalize_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, nfs_c, ntot_amode_c, &
            loffset_c, itermax_c, cloud_borne_c, modal_aerosols_c, inv_so2_c, inv_h2o2_c, inv_o3_c, &
            inv_ho2_c, id_so2_c, id_hno3_c, id_h2o2_c, id_nh3_c, id_o3_c, id_ho2_c, id_h2so4_c, &
            id_so4_c, id_msa_c, modeptr_accum_c, dtime_c, ph0_c, const0_c, kh0_c, kh1_c, kh2_c, kh3_c, &
            ra_c, xkw_c, so4_fact_c, pi_c, small_value_c, specmw_so4_amode_c, gravit_c, xhnm_p, &
            invariants_p, qin_p, cfact_p, xph_p, xso2_p, xhno3_p, xh2o2_p, xnh3_p, xo3_p, xho2_p, &
            xh2so4_p, xso4_p, xno3_p, xnh4_p, xmsa_p, press_p, tfld_p, cldfrc_p, cldnum_p, lwc_p, &
            xlwc_p, xso4c_p, xnh4c_p, xno3c_p, xso4_init_p, xdelso4hp_p, hno3g_p, nh3g_p, hehno3_p, &
            heh2o2_p, heso2_p, henh3_p, heo3_p, xphlwc_p, qcw_p, dqdt_aqso4_p, dqdt_aqh2so4_p, &
            dqdt_aqhprxn_p, dqdt_aqo3rxn_p, faqgain_msa_p, faqgain_so4_p, qnum_c_p, numptrcw_amode_p, &
            lptr_so4_cw_amode_p, lptr_msa_cw_amode_p, lptr_nh4_cw_amode_p, mbar_p, pdel_p, &
            sflx_aqso4_p, sflx_aqh2so4_p, sflx_aqhprxn_p, sflx_aqo3rxn_p, adv_mass_p) &
            bind(c, name="setsox_shell_finalize_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c, nfs_c, ntot_amode_c
         integer(c_int64_t), value :: loffset_c, itermax_c, cloud_borne_c, modal_aerosols_c
         integer(c_int64_t), value :: inv_so2_c, inv_h2o2_c, inv_o3_c, inv_ho2_c
         integer(c_int64_t), value :: id_so2_c, id_hno3_c, id_h2o2_c, id_nh3_c, id_o3_c, id_ho2_c
         integer(c_int64_t), value :: id_h2so4_c, id_so4_c, id_msa_c, modeptr_accum_c
         real(c_double), value :: dtime_c, ph0_c, const0_c, kh0_c, kh1_c, kh2_c, kh3_c, ra_c, xkw_c
         real(c_double), value :: so4_fact_c, pi_c, small_value_c, specmw_so4_amode_c, gravit_c
         type(c_ptr), value :: xhnm_p, invariants_p, qin_p, cfact_p, xph_p, xso2_p, xhno3_p
         type(c_ptr), value :: xh2o2_p, xnh3_p, xo3_p, xho2_p, xh2so4_p, xso4_p, xno3_p, xnh4_p
         type(c_ptr), value :: xmsa_p, press_p, tfld_p, cldfrc_p, cldnum_p, lwc_p, xlwc_p, xso4c_p
         type(c_ptr), value :: xnh4c_p, xno3c_p, xso4_init_p, xdelso4hp_p, hno3g_p, nh3g_p
         type(c_ptr), value :: hehno3_p, heh2o2_p, heso2_p, henh3_p, heo3_p, xphlwc_p, qcw_p
         type(c_ptr), value :: dqdt_aqso4_p, dqdt_aqh2so4_p, dqdt_aqhprxn_p, dqdt_aqo3rxn_p
         type(c_ptr), value :: faqgain_msa_p, faqgain_so4_p, qnum_c_p, numptrcw_amode_p
         type(c_ptr), value :: lptr_so4_cw_amode_p, lptr_msa_cw_amode_p, lptr_nh4_cw_amode_p
         type(c_ptr), value :: mbar_p, pdel_p, sflx_aqso4_p, sflx_aqh2so4_p, sflx_aqhprxn_p
         type(c_ptr), value :: sflx_aqo3rxn_p, adv_mass_p
       end subroutine setsox_shell_finalize_codon
    end interface

    if (masterproc .and. .not. setsox_shell_wrap_proof_written) then
       proof_line = 'setsox_shell_codon_wrap entered (unified setsox shell/update core stage dispatch = codon)'
       write(iulog,'(A)') trim(proof_line)
       call setsox_append_impl_proof(trim(proof_line))
       setsox_shell_wrap_proof_written = .true.
       call flush(iulog)
    end if

    cloud_borne_c = 0_c_int64_t
    modal_aerosols_c = 0_c_int64_t
    inv_so2_c = 0_c_int64_t
    inv_h2o2_c = 0_c_int64_t
    inv_o3_c = 0_c_int64_t
    inv_ho2_c = 0_c_int64_t
    if (cloud_borne) cloud_borne_c = 1_c_int64_t
    if (modal_aerosols) modal_aerosols_c = 1_c_int64_t
    if (inv_so2) inv_so2_c = 1_c_int64_t
    if (inv_h2o2) inv_h2o2_c = 1_c_int64_t
    if (inv_o3) inv_o3_c = 1_c_int64_t
    if (inv_ho2) inv_ho2_c = 1_c_int64_t

    call setsox_shell_stage_dispatch_codon( &
         1_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(gas_pcnst, c_int64_t), int(nfs, c_int64_t), int(ntot_amode, c_int64_t), int(loffset, c_int64_t), &
         int(itermax, c_int64_t), cloud_borne_c, modal_aerosols_c, inv_so2_c, inv_h2o2_c, inv_o3_c, inv_ho2_c, &
         int(id_so2, c_int64_t), int(id_hno3, c_int64_t), int(id_h2o2, c_int64_t), int(id_nh3, c_int64_t), &
         int(id_o3, c_int64_t), int(id_ho2, c_int64_t), int(id_h2so4, c_int64_t), int(id_so4, c_int64_t), &
         int(id_msa, c_int64_t), int(modeptr_accum, c_int64_t), real(dtime, c_double), real(ph0, c_double), &
         real(const0, c_double), real(kh0, c_double), real(kh1, c_double), real(kh2, c_double), real(kh3, c_double), &
         real(ra, c_double), real(xkw, c_double), real(0._r8, c_double), real(pi, c_double), c_loc(xhnm), &
         c_loc(invariants), c_loc(qin), &
         c_loc(cfact), c_loc(xph), c_loc(xso2), c_loc(xhno3), c_loc(xh2o2), c_loc(xnh3), c_loc(xo3), c_loc(xho2), &
         c_loc(xh2so4), c_loc(xso4), c_loc(xno3), c_loc(xnh4), c_loc(xmsa), c_loc(press), c_loc(tfld), c_loc(cldfrc), &
         c_loc(cldnum), c_loc(lwc), c_loc(lwc), c_loc(qcw), c_loc(qcw), c_loc(qcw), c_loc(xso4_init), c_loc(xdelso4hp), &
         c_loc(hno3g), c_loc(nh3g), c_loc(hehno3), c_loc(heh2o2), c_loc(heso2), c_loc(henh3), c_loc(heo3), c_loc(xphlwc), &
         c_loc(qcw), c_loc(dqdt_aqso4), c_loc(dqdt_aqh2so4), c_loc(dqdt_aqhprxn), c_loc(dqdt_aqo3rxn), c_loc(faqgain_msa), &
         c_loc(faqgain_so4), c_loc(qnum_c), c_loc(numptrcw_amode_c), c_loc(lptr_so4_cw_amode_c), c_loc(lptr_msa_cw_amode_c), &
         c_loc(lptr_nh4_cw_amode_c) )

    cldconc => sox_cldaero_create_obj(cldfrc, qcw, lwc, cfact, ncol, loffset)
    xso4c => cldconc%so4c
    xnh4c => cldconc%nh4c
    xno3c => cldconc%no3c

    do n = 1, ntot_amode
       numptrcw_amode_c(n) = int(numptrcw_amode(n), c_int64_t)
       lptr_so4_cw_amode_c(n) = int(lptr_so4_cw_amode(n), c_int64_t)
       lptr_msa_cw_amode_c(n) = int(lptr_msa_cw_amode(n), c_int64_t)
       lptr_nh4_cw_amode_c(n) = int(lptr_nh4_cw_amode(n), c_int64_t)
    end do
    adv_mass_work(:) = adv_mass(:)

    if (masterproc .and. .not. setsox_shell_finalize_wrap_proof_written) then
       proof_line = 'setsox shell/finalize dispatcher entered (setsox stage2/aq finalize direct = codon; outfld = native)'
       write(iulog,'(A)') trim(proof_line)
       call setsox_append_impl_proof(trim(proof_line))
       setsox_shell_finalize_wrap_proof_written = .true.
       call flush(iulog)
    end if

    call setsox_shell_finalize_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(gas_pcnst, c_int64_t), int(nfs, c_int64_t), int(ntot_amode, c_int64_t), int(loffset, c_int64_t), &
         int(itermax, c_int64_t), cloud_borne_c, modal_aerosols_c, inv_so2_c, inv_h2o2_c, inv_o3_c, inv_ho2_c, &
         int(id_so2, c_int64_t), int(id_hno3, c_int64_t), int(id_h2o2, c_int64_t), int(id_nh3, c_int64_t), &
         int(id_o3, c_int64_t), int(id_ho2, c_int64_t), int(id_h2so4, c_int64_t), int(id_so4, c_int64_t), &
         int(id_msa, c_int64_t), int(modeptr_accum, c_int64_t), real(dtime, c_double), real(ph0, c_double), &
         real(const0, c_double), real(kh0, c_double), real(kh1, c_double), real(kh2, c_double), real(kh3, c_double), &
         real(ra, c_double), real(xkw, c_double), real(cldconc%so4_fact, c_double), real(pi, c_double), &
         real(1.e-20_r8, c_double), real(specmw_so4_amode, c_double), real(gravit, c_double), &
         c_loc(xhnm), c_loc(invariants), c_loc(qin), c_loc(cfact), c_loc(xph), c_loc(xso2), c_loc(xhno3), &
         c_loc(xh2o2), c_loc(xnh3), c_loc(xo3), c_loc(xho2), c_loc(xh2so4), c_loc(xso4), c_loc(xno3), &
         c_loc(xnh4), c_loc(xmsa), c_loc(press), c_loc(tfld), c_loc(cldfrc), c_loc(cldnum), c_loc(lwc), &
         c_loc(cldconc%xlwc), c_loc(xso4c), c_loc(xnh4c), c_loc(xno3c), c_loc(xso4_init), c_loc(xdelso4hp), &
         c_loc(hno3g), c_loc(nh3g), c_loc(hehno3), c_loc(heh2o2), c_loc(heso2), &
         c_loc(henh3), c_loc(heo3), c_loc(xphlwc), c_loc(qcw), c_loc(dqdt_aqso4), c_loc(dqdt_aqh2so4), c_loc(dqdt_aqhprxn), &
         c_loc(dqdt_aqo3rxn), c_loc(faqgain_msa), c_loc(faqgain_so4), c_loc(qnum_c), c_loc(numptrcw_amode_c), &
         c_loc(lptr_so4_cw_amode_c), c_loc(lptr_msa_cw_amode_c), c_loc(lptr_nh4_cw_amode_c), &
         c_loc(mbar), c_loc(pdel), c_loc(sflx_aqso4), c_loc(sflx_aqh2so4), c_loc(sflx_aqhprxn), &
         c_loc(sflx_aqo3rxn), c_loc(adv_mass_work) )

    do n = 1, ntot_amode
       m = lptr_so4_cw_amode(n)
       l = m - loffset
       if (l > 0) then
          call outfld( trim(cnst_name_cw(m))//'AQSO4', sflx_aqso4(:ncol,n), ncol, lchnk)
          call outfld( trim(cnst_name_cw(m))//'AQH2SO4', sflx_aqh2so4(:ncol,n), ncol, lchnk)
       endif
    end do

    call outfld( 'AQSO4_H2O2', sflx_aqhprxn(:ncol), ncol, lchnk)
    call outfld( 'AQSO4_O3', sflx_aqo3rxn(:ncol), ncol, lchnk)
    call outfld('XPH_LWC', xphlwc(:ncol,:), ncol, lchnk)

    call sox_cldaero_destroy_obj(cldconc)

  end subroutine setsox_shell_codon_wrap

!-----------------------------------------------------------------------      
!-----------------------------------------------------------------------      
  subroutine sox_inti
    !-----------------------------------------------------------------------      
    !	... initialize the hetero sox routine
    !-----------------------------------------------------------------------      

    use mo_chem_utls, only : get_spc_ndx, get_inv_ndx
    use spmd_utils,   only : masterproc
    use cam_history,  only : addfld,phys_decomp
    use cam_history,  only : add_default
    use ppgrid,       only : pver
    use phys_control, only : phys_getopts
    use sox_cldaero_mod, only : sox_cldaero_init

    implicit none

    logical :: history_aerosol   ! Output aerosol diagnostics
    integer(c_int64_t) :: init_active_c

    init_active_c = sox_inti_active_codon(1_c_int64_t)
    if (init_active_c == 0_c_int64_t) return

    if (masterproc .and. .not. sox_inti_proof_written) then
       write(iulog,'(A)') 'sox_inti direct = codon; SOX active-policy direct; phys_getopts/species/history native CAM API islands'
       call setsox_append_impl_proof('sox_inti direct = codon; SOX active-policy direct; ' // &
            'phys_getopts/species/history native CAM API islands')
       sox_inti_proof_written = .true.
       call flush(iulog)
    end if

    call phys_getopts( &
         history_aerosol_out = history_aerosol, &
         prog_modal_aero_out=modal_aerosols )

    cloud_borne = modal_aerosols

    !-----------------------------------------------------------------
    !       ... get species indicies
    !-----------------------------------------------------------------
    
    if (cloud_borne) then
       id_h2so4 = get_spc_ndx( 'H2SO4' )
    else
       id_so4 = get_spc_ndx( 'SO4' )
    endif
    id_msa = get_spc_ndx( 'MSA' )

    inv_so2 = .false.
    id_so2 = get_inv_ndx( 'SO2' )
    inv_so2 = id_so2 > 0
    if ( .not. inv_so2 ) then
       id_so2 = get_spc_ndx( 'SO2' )
    endif

    inv_NH3 = .false.
    id_NH3 = get_inv_ndx( 'NH3' )
    inv_NH3 = id_NH3 > 0
    if ( .not. inv_NH3 ) then
       id_NH3 = get_spc_ndx( 'NH3' )
    endif

    inv_HNO3 = .false.
    id_HNO3 = get_inv_ndx( 'HNO3' )
    inv_HNO3 = id_hno3 > 0
    if ( .not. inv_HNO3 ) then
       id_HNO3 = get_spc_ndx( 'HNO3' )
    endif

    inv_H2O2 = .false.
    id_H2O2 = get_inv_ndx( 'H2O2' )
    inv_H2O2 = id_H2O2 > 0
    if ( .not. inv_H2O2 ) then
       id_H2O2 = get_spc_ndx( 'H2O2' )
    endif

    inv_HO2 = .false.
    id_HO2 = get_inv_ndx( 'HO2' )
    inv_HO2 = id_HO2 > 0
    if ( .not. inv_HO2 ) then
       id_HO2 = get_spc_ndx( 'HO2' )
    endif

    inv_o3 = get_inv_ndx( 'O3' ) > 0
    if (inv_o3) then
       id_o3 = get_inv_ndx( 'O3' )
    else
       id_o3 = get_spc_ndx( 'O3' )
    endif
    inv_ho2 = get_inv_ndx( 'HO2' ) > 0
    if (inv_ho2) then
       id_ho2 = get_inv_ndx( 'HO2' )
    else
       id_ho2 = get_spc_ndx( 'HO2' )
    endif

    has_sox = (id_so2>0) .and. (id_h2o2>0) .and. (id_o3>0) .and. (id_ho2>0)
    if (cloud_borne) then
       has_sox = has_sox .and. (id_h2so4>0)
    else
       has_sox = has_sox .and. (id_so4>0) .and. (id_nh3>0)
    endif

    if (masterproc) then
       write(iulog,*) 'sox_inti: has_sox = ',has_sox
    endif

    if( has_sox ) then
       if (masterproc) then
          write(iulog,*) '-----------------------------------------'
          write(iulog,*) 'mozart will do sox aerosols'
          write(iulog,*) '-----------------------------------------'
       endif
    else 
       return
    end if

    call addfld( 'XPH_LWC','kg/kg   ',pver, 'A', 'pH value multiplied by lwc', phys_decomp)
    if ( history_aerosol ) then    
       call add_default ('XPH_LWC', 1, ' ') 
    endif

    call sox_cldaero_init()

  end subroutine sox_inti
  
!-----------------------------------------------------------------------      
!-----------------------------------------------------------------------      
  subroutine SETSOX( &
       ncol,   &
       lchnk,  &
       loffset,&
       dtime,  &
       press,  &
       pdel,   &
       tfld,   &
       mbar,   &
       lwc,    &
       cldfrc, &
       cldnum, &
       xhnm,   &
       invariants, &
       qcw,    &
       qin     &
       )

    !-----------------------------------------------------------------------      
    !          ... Compute heterogeneous reactions of SOX
    !
    !       (0) using initial PH to calculate PH
    !           (a) HENRYs law constants
    !           (b) PARTIONING
    !           (c) PH values
    !
    !       (1) using new PH to repeat
    !           (a) HENRYs law constants
    !           (b) PARTIONING
    !           (c) REACTION rates
    !           (d) PREDICTION
    !-----------------------------------------------------------------------      
    !
    use ppgrid,    only : pcols, pver
    use chem_mods, only : gas_pcnst, nfs
    use chem_mods,    only : adv_mass
    use physconst,    only : mwdry, gravit
    use mo_constants, only : pi
    use cam_history,  only : outfld
    use sox_cldaero_mod, only : sox_cldaero_update, sox_cldaero_create_obj, sox_cldaero_destroy_obj
    use cldaero_mod,     only : cldaero_conc_t

    !
    implicit none
    !
    !-----------------------------------------------------------------------      
    !      ... Dummy arguments
    !-----------------------------------------------------------------------      
    integer,          intent(in)    :: ncol              ! num of columns in chunk
    integer,          intent(in)    :: lchnk             ! chunk id
    integer,          intent(in)    :: loffset           ! offset of chem tracers in the advected tracers array
    real(r8),         intent(in)    :: dtime             ! time step (sec)
    real(r8),         intent(in)    :: press(:,:)        ! midpoint pressure ( Pa )
    real(r8),         intent(in)    :: pdel(:,:)         ! pressure thickness of levels (Pa)
    real(r8),         intent(in)    :: tfld(:,:)         ! temperature
    real(r8),         intent(in)    :: mbar(:,:)         ! mean wet atmospheric mass ( amu )
    real(r8), target, intent(in)    :: lwc(:,:)          ! cloud liquid water content (kg/kg)
    real(r8), target, intent(in)    :: cldfrc(:,:)       ! cloud fraction
    real(r8),         intent(in)    :: cldnum(:,:)       ! droplet number concentration (#/kg)
    real(r8),         intent(in)    :: xhnm(:,:)         ! total atms density ( /cm**3)
    real(r8),         intent(in)    :: invariants(:,:,:)
    real(r8), target, intent(inout) :: qcw(:,:,:)        ! cloud-borne aerosol (vmr)
    real(r8),         intent(inout) :: qin(:,:,:)        ! transported species ( vmr )

    !-----------------------------------------------------------------------      
    !      ... Local variables
    !
    !           xhno3 ... in mixing ratio
    !-----------------------------------------------------------------------      
    integer,  parameter :: itermax = 20
    real(r8), parameter :: ph0 = 5.0_r8  ! INITIAL PH VALUES
    real(r8), parameter :: const0 = 1.e3_r8/6.023e23_r8
    real(r8), parameter :: xa0 = 11._r8
    real(r8), parameter :: xb0 = -.1_r8
    real(r8), parameter :: xa1 = 1.053_r8
    real(r8), parameter :: xb1 = -4.368_r8
    real(r8), parameter :: xa2 = 1.016_r8
    real(r8), parameter :: xb2 = -2.54_r8
    real(r8), parameter :: xa3 = .816e-32_r8
    real(r8), parameter :: xb3 = .259_r8

    real(r8), parameter :: kh0 = 9.e3_r8            ! HO2(g)          -> Ho2(a)
    real(r8), parameter :: kh1 = 2.05e-5_r8         ! HO2(a)          -> H+ + O2-
    real(r8), parameter :: kh2 = 8.6e5_r8           ! HO2(a) + ho2(a) -> h2o2(a) + o2
    real(r8), parameter :: kh3 = 1.e8_r8            ! HO2(a) + o2-    -> h2o2(a) + o2
    real(r8), parameter :: Ra = 8314._r8/101325._r8 ! universal constant   (atm)/(M-K)
    real(r8), parameter :: xkw = 1.e-14_r8          ! water acidity

    !
    real(r8) :: xdelso4hp(ncol,pver)
    real(r8) :: xphlwc(ncol,pver)

    integer  :: k, i, iter, file
    real(r8) :: wrk, delta
    real(r8) :: xph0, aden, xk, xe, x2
    real(r8) :: tz, xl, px, qz, pz, es, qs, patm
    real(r8) :: Eso2, Eso4, Ehno3, Eco2, Eh2o, Enh3
    real(r8) :: so2g, h2o2g, co2g, o3g
    real(r8) :: hno3a, nh3a, so2a, h2o2a, co2a, o3a
    real(r8) :: rah2o2, rao3, pso4, ccc
    real(r8) :: cnh3, chno3, com, com1, com2, xra

    real(r8) :: hno3g(ncol,pver), nh3g(ncol,pver)
    !
    !-----------------------------------------------------------------------      
    !            for Ho2(g) -> H2o2(a) formation 
    !            schwartz JGR, 1984, 11589
    !-----------------------------------------------------------------------      
    real(r8) :: kh4    ! kh2+kh3
    real(r8) :: xam    ! air density /cm3
    real(r8) :: ho2s   ! ho2s = ho2(a)+o2-
    real(r8) :: r1h2o2 ! prod(h2o2) by ho2 in mole/L(w)/s
    real(r8) :: r2h2o2 ! prod(h2o2) by ho2 in mix/s

    real(r8), dimension(ncol,pver)  ::             &
         xhno3, xh2o2, xso2, xso4, xno3, &
         xnh3, xnh4, xo3,         &
         cfact, &
         xph, xho2,         &
         xh2so4, xmsa, xso4_init, &
         hehno3, &            ! henry law const for hno3
         heh2o2, &            ! henry law const for h2o2
         heso2,  &            ! henry law const for so2
         henh3,  &            ! henry law const for nh3
         heo3              !!,   &            ! henry law const for o3

    real(r8) :: patm_x

    real(r8), dimension(ncol)  :: work1
    logical :: converged

    real(r8), pointer :: xso4c(:,:)
    real(r8), pointer :: xnh4c(:,:)
    real(r8), pointer :: xno3c(:,:)
    type(cldaero_conc_t), pointer :: cldconc

    real(r8) :: fact1_hno3, fact2_hno3, fact3_hno3
    real(r8) :: fact1_so2, fact2_so2, fact3_so2, fact4_so2
    real(r8) :: fact1_nh3, fact2_nh3, fact3_nh3
    real(r8) :: tmp_hp, tmp_hso3, tmp_hco3, tmp_nh4, tmp_no3
    real(r8) :: tmp_oh, tmp_so3, tmp_so4
    real(r8) :: tmp_neg, tmp_pos
    real(r8) :: yph, yph_lo, yph_hi
    real(r8) :: ynetpos, ynetpos_lo, ynetpos_hi

    !-----------------------------------------------------------------
    !       ... NOTE: The press array is in pascals and must be
    !                 mutiplied by 10 to yield dynes/cm**2.
    !-----------------------------------------------------------------
    !==================================================================
    !       ... First set the PH
    !==================================================================
    !      ... Initial values
    !           The values of so2, so4 are after (1) SLT, and CHEM
    !-----------------------------------------------------------------
    xph0 = 10._r8**(-ph0)                      ! initial PH value

    call setsox_init_fields_select_impl()
    if (setsox_init_fields_use_native_impl) then
       do k = 1,pver
          cfact(:,k) = xhnm(:,k)     &          ! /cm3(a)
               * 1.e6_r8             &          ! /m3(a)
               * 1.38e-23_r8/287._r8 &          ! Kg(a)/m3(a)
               * 1.e-3_r8                       ! Kg(a)/L(a)
       end do
    else
       call setsox_init_fields_codon_wrap(1, ncol, ph0, xhnm, invariants, qin, cfact, xph, &
            xso2, xhno3, xh2o2, xnh3, xo3, xho2, xh2so4, xso4, xno3, xnh4, xmsa)
    end if

    cldconc => sox_cldaero_create_obj( cldfrc,qcw,lwc, cfact, ncol, loffset )
    xso4c => cldconc%so4c
    xnh4c => cldconc%nh4c
    xno3c => cldconc%no3c

    if (setsox_init_fields_use_native_impl) then
       xso4(:,:) = 0._r8
       xno3(:,:) = 0._r8
       xnh4(:,:) = 0._r8

       do k = 1,pver
          xph(:,k) = xph0                                ! initial PH value

          if ( inv_so2 ) then
             xso2 (:,k) = invariants(:,k,id_so2)/xhnm(:,k)  ! mixing ratio
          else
             xso2 (:,k) = qin(:,k,id_so2)                   ! mixing ratio
          endif

          if (id_hno3 > 0) then
             xhno3(:,k) = qin(:,k,id_hno3)
          else
             xhno3(:,k) = 0.0_r8
          endif

          if ( inv_h2o2 ) then
             xh2o2 (:,k) = invariants(:,k,id_h2o2)/xhnm(:,k)  ! mixing ratio
          else
             xh2o2 (:,k) = qin(:,k,id_h2o2)                   ! mixing ratio
          endif

          if (id_nh3  > 0) then
             xnh3 (:,k) = qin(:,k,id_nh3)
          else
             xnh3 (:,k) = 0.0_r8
          endif

          if ( inv_o3 ) then
             xo3  (:,k) = invariants(:,k,id_o3)/xhnm(:,k) ! mixing ratio
          else
             xo3  (:,k) = qin(:,k,id_o3)                  ! mixing ratio
          endif
          if ( inv_ho2 ) then
             xho2 (:,k) = invariants(:,k,id_ho2)/xhnm(:,k)! mixing ratio
          else
             xho2 (:,k) = qin(:,k,id_ho2)                 ! mixing ratio
          endif

          if (cloud_borne) then
             xh2so4(:,k) = qin(:,k,id_h2so4)
          else
             xso4  (:,k) = qin(:,k,id_so4) ! mixing ratio
          endif
          if (id_msa > 0) xmsa (:,k) = qin(:,k,id_msa)

       end do
    else
       call setsox_init_fields_codon_wrap(2, ncol, ph0, xhnm, invariants, qin, cfact, xph, &
            xso2, xhno3, xh2o2, xnh3, xo3, xho2, xh2so4, xso4, xno3, xnh4, xmsa)
    end if
    
    !-----------------------------------------------------------------
    !       ... Temperature dependent Henry constants
    !-----------------------------------------------------------------
    call setsox_ph_solve_select_impl()
    if (setsox_ph_solve_use_native_impl) then
    ver_loop0: do k = 1,pver                               !! pver loop for STEP 0
       col_loop0: do i = 1,ncol
          
          if (cloud_borne .and. cldfrc(i,k)>0._r8) then
             xso4(i,k) = xso4c(i,k) / cldfrc(i,k)
             xnh4(i,k) = xnh4c(i,k) / cldfrc(i,k)
             xno3(i,k) = xno3c(i,k) / cldfrc(i,k)
          endif
          xl = cldconc%xlwc(i,k)

          if( xl >= 1.e-8_r8 ) then
             work1(i) = 1._r8 / tfld(i,k) - 1._r8 / 298._r8

             !-----------------------------------------------------------------
             ! 21-mar-2011 changes by rce
             ! ph calculation now uses bisection method to solve the electro-neutrality equation
             ! 3-mode aerosols (where so4 is assumed to be nh4hso4)
             !    old code set xnh4c = so4c
             !    new code sets xnh4c = 0, then uses a -1 charge (instead of -2)
             !       for so4 when solving the electro-neutrality equation
             !-----------------------------------------------------------------

             !-----------------------------------------------------------------
             !  calculations done before iterating
             !-----------------------------------------------------------------

             !-----------------------------------------------------------------
             pz = .01_r8*press(i,k)       !! pressure in mb
             tz = tfld(i,k)
             patm = pz/1013._r8
             xam  = press(i,k)/(1.38e-23_r8*tz)  !air density /M3

             !-----------------------------------------------------------------
             !        ... hno3
             !-----------------------------------------------------------------
             ! previous code
             !    hehno3(i,k)  = xk*(1._r8 + xe/xph(i,k))
             !    px = hehno3(i,k) * Ra * tz * xl
             !    hno3g = xhno3(i,k)/(1._r8 + px)
             !    Ehno3 = xk*xe*hno3g *patm
             ! equivalent new code
             !    hehno3 = xk + xk*xe/hplus
             !    hno3g = xhno3/(1 + px)
             !          = xhno3/(1 + hehno3*ra*tz*xl)
             !          = xhno3/(1 + xk*ra*tz*xl*(1 + xe/hplus)
             !    ehno3 = hno3g*xk*xe*patm
             !          = xk*xe*patm*xhno3/(1 + xk*ra*tz*xl*(1 + xe/hplus)
             !          = ( fact1_hno3    )/(1 + fact2_hno3 *(1 + fact3_hno3/hplus)
             !    [hno3-] = ehno3/hplus
             xk = 2.1e5_r8 *EXP( 8700._r8*work1(i) )
             xe = 15.4_r8
             fact1_hno3 = xk*xe*patm*xhno3(i,k)
             fact2_hno3 = xk*ra*tz*xl
             fact3_hno3 = xe

             !-----------------------------------------------------------------
             !          ... so2
             !-----------------------------------------------------------------
             ! previous code
             !    heso2(i,k)  = xk*(1._r8 + wrk*(1._r8 + x2/xph(i,k)))
             !    px = heso2(i,k) * Ra * tz * xl
             !    so2g =  xso2(i,k)/(1._r8+ px)
             !    Eso2 = xk*xe*so2g *patm
             ! equivalent new code
             !    heso2 = xk + xk*xe/hplus * xk*xe*x2/hplus**2
             !    so2g = xso2/(1 + px)
             !         = xso2/(1 + heso2*ra*tz*xl)
             !         = xso2/(1 + xk*ra*tz*xl*(1 + (xe/hplus)*(1 + x2/hplus))
             !    eso2 = so2g*xk*xe*patm
             !          = xk*xe*patm*xso2/(1 + xk*ra*tz*xl*(1 + (xe/hplus)*(1 + x2/hplus))
             !          = ( fact1_so2    )/(1 + fact2_so2 *(1 + (fact3_so2/hplus)*(1 + fact4_so2/hplus)
             !    [hso3-] + 2*[so3--] = (eso2/hplus)*(1 + 2*x2/hplus)
             xk = 1.23_r8  *EXP( 3120._r8*work1(i) )
             xe = 1.7e-2_r8*EXP( 2090._r8*work1(i) )
             x2 = 6.0e-8_r8*EXP( 1120._r8*work1(i) )
             fact1_so2 = xk*xe*patm*xso2(i,k)
             fact2_so2 = xk*ra*tz*xl
             fact3_so2 = xe
             fact4_so2 = x2

             !-----------------------------------------------------------------
             !          ... nh3
             !-----------------------------------------------------------------
             ! previous code
             !    henh3(i,k)  = xk*(1._r8 + xe*xph(i,k)/xkw)
             !    px = henh3(i,k) * Ra * tz * xl
             !    nh3g = (xnh3(i,k)+xnh4(i,k))/(1._r8+ px)
             !    Enh3 = xk*xe*nh3g/xkw *patm
             ! equivalent new code
             !    henh3 = xk + xk*xe*hplus/xkw
             !    nh3g = xnh34/(1 + px)
             !         = xnh34/(1 + henh3*ra*tz*xl)
             !         = xnh34/(1 + xk*ra*tz*xl*(1 + xe*hplus/xkw)
             !    enh3 = nh3g*xk*xe*patm/xkw
             !          = ((xk*xe*patm/xkw)*xnh34)/(1 + xk*ra*tz*xl*(1 + xe*hplus/xkw)
             !          = ( fact1_nh3            )/(1 + fact2_nh3  *(1 + fact3_nh3*hplus)
             !    [nh4+] = enh3*hplus
             xk = 58._r8   *EXP( 4085._r8*work1(i) )
             xe = 1.7e-5_r8*EXP( -4325._r8*work1(i) )

             fact1_nh3 = (xk*xe*patm/xkw)*(xnh3(i,k)+xnh4(i,k))
             fact2_nh3 = xk*ra*tz*xl
             fact3_nh3 = xe/xkw

             !-----------------------------------------------------------------
             !        ... h2o effects
             !-----------------------------------------------------------------
             Eh2o = xkw

             !-----------------------------------------------------------------
             !        ... co2 effects
             !-----------------------------------------------------------------
             co2g = 330.e-6_r8                            !330 ppm = 330.e-6 atm
             xk = 3.1e-2_r8*EXP( 2423._r8*work1(i) )
             xe = 4.3e-7_r8*EXP(-913._r8 *work1(i) )
             Eco2 = xk*xe*co2g  *patm

             !-----------------------------------------------------------------
             !         ... so4 effect
             !-----------------------------------------------------------------
             Eso4 = xso4(i,k)*xhnm(i,k)   &         ! /cm3(a)
                  *const0/xl


             !-----------------------------------------------------------------
             ! now use bisection method to solve electro-neutrality equation
             !
             ! during the iteration loop,
             !    yph_lo = lower ph value that brackets the root (i.e., correct ph)
             !    yph_hi = upper ph value that brackets the root (i.e., correct ph)
             !    yph    = current ph value
             !    yposnet_lo and yposnet_hi = net positive ions for
             !       yph_lo and yph_hi
             !-----------------------------------------------------------------
             do iter = 1,itermax

                if (iter == 1) then
                   ! 1st iteration ph = lower bound value
                   yph_lo = 2.0_r8
                   yph_hi = yph_lo
                   yph = yph_lo
                else if (iter == 2) then
                   ! 2nd iteration ph = upper bound value
                   yph_hi = 7.0_r8
                   yph = yph_hi
                else
                   ! later iteration ph = mean of the two bracketing values
                   yph = 0.5_r8*(yph_lo + yph_hi)
                end if
                ! calc current [H+] from ph
                xph(i,k) = 10.0_r8**(-yph)

                !-----------------------------------------------------------------
                !        ... hno3
                !-----------------------------------------------------------------
                Ehno3 = fact1_hno3/(1.0_r8 + fact2_hno3*(1.0_r8 + fact3_hno3/xph(i,k)))

                !-----------------------------------------------------------------
                !          ... so2
                !-----------------------------------------------------------------
                Eso2 = fact1_so2/(1.0_r8 + fact2_so2*(1.0_r8 + (fact3_so2/xph(i,k)) &
                     *(1.0_r8 +  fact4_so2/xph(i,k))))

                !-----------------------------------------------------------------
                !          ... nh3
                !-----------------------------------------------------------------
                Enh3 = fact1_nh3/(1.0_r8 + fact2_nh3*(1.0_r8 + fact3_nh3*xph(i,k)))

                tmp_nh4  = Enh3 * xph(i,k)
                tmp_hso3 = Eso2 / xph(i,k)
                tmp_so3  = tmp_hso3 * 2.0_r8*fact4_so2/xph(i,k)
                tmp_hco3 = Eco2 / xph(i,k)
                tmp_oh   = Eh2o / xph(i,k)
                tmp_no3  = Ehno3 / xph(i,k)
                tmp_so4 = cldconc%so4_fact*Eso4
                tmp_pos = xph(i,k) + tmp_nh4
                tmp_neg = tmp_oh + tmp_hco3 + tmp_no3 + tmp_hso3 + tmp_so3 + tmp_so4

                ynetpos = tmp_pos - tmp_neg


                ! yposnet = net positive ions/charge
                ! if the correct ph is bracketed by yph_lo and yph_hi (with yph_lo < yph_hi),
                !    then you will have yposnet_lo > 0 and yposnet_hi < 0
                converged = .false.
                if (iter > 2) then
                   if (ynetpos == 0.0_r8) then
                      ! the exact solution was found (very unlikely)
                      tmp_hp = xph(i,k)
                      converged = .true.
                      exit
                   else if (ynetpos >= 0.0_r8) then
                      ! net positive ions are >= 0 for both yph and yph_lo
                      !    so replace yph_lo with yph
                      yph_lo = yph
                      ynetpos_lo = ynetpos
                   else
                      ! net positive ions are <= 0 for both yph and yph_hi
                      !    so replace yph_hi with yph
                      yph_hi = yph
                      ynetpos_hi = ynetpos
                   end if

                   if (abs(yph_hi - yph_lo) .le. 0.005_r8) then
                      ! |yph_hi - yph_lo| <= convergence criterion, so set
                      !    final ph to their midpoint and exit
                      ! (.005 absolute error in pH gives .01 relative error in H+)
                      tmp_hp = xph(i,k)
                      yph = 0.5_r8*(yph_hi + yph_lo)
                      xph(i,k) = 10.0_r8**(-yph)
                      converged = .true.
                      exit
                   else 
                      ! do another iteration
                      converged = .false.
                   end if

                else if (iter == 1) then
                   if (ynetpos <= 0.0_r8) then
                      ! the lower and upper bound ph values (2.0 and 7.0) do not bracket
                      !    the correct ph, so use the lower bound
                      tmp_hp = xph(i,k)
                      converged = .true.
                      exit
                   end if
                   ynetpos_lo = ynetpos

                else ! (iter == 2)
                   if (ynetpos >= 0.0_r8) then
                      ! the lower and upper bound ph values (2.0 and 7.0) do not bracket
                      !    the correct ph, so use they upper bound
                      tmp_hp = xph(i,k)
                      converged = .true.
                      exit
                   end if
                   ynetpos_hi = ynetpos
                end if

             end do ! iter

             if( .not. converged ) then
                write(iulog,*) 'SETSOX: pH failed to converge @ (',i,',',k,'), % change=', &
                     100._r8*delta
             end if
          else
             xph(i,k) =  1.e-7_r8
          end if
       end do col_loop0
    end do ver_loop0 ! end pver loop for STEP 0
    else
       call setsox_ph_solve_codon_wrap(ncol, itermax, const0, Ra, xkw, press, tfld, cldfrc, xhnm, &
            cldconc%xlwc, xso4c, xnh4c, xno3c, xso4, xnh4, xno3, xso2, xhno3, xnh3, xph, &
            cldconc%so4_fact)
    end if

    !==============================================================
    !          ... Now use the actual PH
    !==============================================================
    call setsox_aqchem_predict_select_impl()
    if (setsox_aqchem_predict_use_native_impl) then
    ver_loop1: do k = 1,pver
       col_loop1: do i = 1,ncol
          work1(i) = 1._r8 / tfld(i,k) - 1._r8 / 298._r8
          tz = tfld(i,k)

          xl = cldconc%xlwc(i,k)

          patm = press(i,k)/101300._r8        ! press is in pascal
          xam  = press(i,k)/(1.38e-23_r8*tz)  ! air density /M3

          !-----------------------------------------------------------------------      
          !        ... hno3
          !-----------------------------------------------------------------------      
          xk = 2.1e5_r8 *EXP( 8700._r8*work1(i) )
          xe = 15.4_r8
          hehno3(i,k)  = xk*(1._r8 + xe/xph(i,k))

          !-----------------------------------------------------------------
          !        ... h2o2
          !-----------------------------------------------------------------
          xk = 7.4e4_r8   *EXP( 6621._r8*work1(i) )
          xe = 2.2e-12_r8 *EXP(-3730._r8*work1(i) )
          heh2o2(i,k)  = xk*(1._r8 + xe/xph(i,k))

          !-----------------------------------------------------------------
          !         ... so2
          !-----------------------------------------------------------------
          xk = 1.23_r8  *EXP( 3120._r8*work1(i) )
          xe = 1.7e-2_r8*EXP( 2090._r8*work1(i) )
          x2 = 6.0e-8_r8*EXP( 1120._r8*work1(i) )

          wrk = xe/xph(i,k)
          heso2(i,k)  = xk*(1._r8 + wrk*(1._r8 + x2/xph(i,k)))

          !-----------------------------------------------------------------
          !          ... nh3
          !-----------------------------------------------------------------
          xk = 58._r8   *EXP( 4085._r8*work1(i) )
          xe = 1.7e-5_r8*EXP(-4325._r8*work1(i) )
          henh3(i,k)  = xk*(1._r8 + xe*xph(i,k)/xkw)

          !-----------------------------------------------------------------
          !        ... o3
          !-----------------------------------------------------------------
          xk = 1.15e-2_r8 *EXP( 2560._r8*work1(i) )
          heo3(i,k) = xk

          !------------------------------------------------------------------------
          !       ... for Ho2(g) -> H2o2(a) formation 
          !           schwartz JGR, 1984, 11589
          !------------------------------------------------------------------------
          kh4 = (kh2 + kh3*kh1/xph(i,k)) / ((1._r8 + kh1/xph(i,k))**2)
          ho2s = kh0*xho2(i,k)*patm*(1._r8 + kh1/xph(i,k))  ! ho2s = ho2(a)+o2-
          r1h2o2 = kh4*ho2s*ho2s                         ! prod(h2o2) in mole/L(w)/s

          if ( cloud_borne ) then
             r2h2o2 = r1h2o2*xl        &    ! mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                  / const0*1.e+6_r8  &    ! correct a bug here ????
                  / xam
          else
             r2h2o2 = r1h2o2*xl  &          ! mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                  * const0     &          ! mole/fm3(a)/s * 1.e-3       = mole/cm3(a)/s
                  / xam                   ! /cm3(a)/s    / air-den     = mix-ratio/s
          endif

          if ( .not. modal_aerosols ) then
             xh2o2(i,k) = xh2o2(i,k) + r2h2o2*dtime         ! updated h2o2 by het production
          endif

          !-----------------------------------------------
          !       ... Partioning 
          !-----------------------------------------------

          !-----------------------------------------------------------------
          !        ... hno3
          !-----------------------------------------------------------------
          px = hehno3(i,k) * Ra * tz * xl
          hno3g(i,k) = (xhno3(i,k)+xno3(i,k))/(1._r8 + px)

          !------------------------------------------------------------------------
          !        ... h2o2
          !------------------------------------------------------------------------
          px = heh2o2(i,k) * Ra * tz * xl
          h2o2g =  xh2o2(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... so2
          !------------------------------------------------------------------------
          px = heso2(i,k) * Ra * tz * xl
          so2g =  xso2(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... o3
          !------------------------------------------------------------------------
          px = heo3(i,k) * Ra * tz * xl
          o3g =  xo3(i,k)/(1._r8+ px)

          !------------------------------------------------------------------------
          !         ... nh3
          !------------------------------------------------------------------------
          px = henh3(i,k) * Ra * tz * xl
          if (id_nh3>0) then
             nh3g(i,k) = (xnh3(i,k)+xnh4(i,k))/(1._r8+ px)
          else
             nh3g(i,k) = 0._r8
          endif

          !-----------------------------------------------
          !       ... Aqueous phase reaction rates
          !           SO2 + H2O2 -> SO4
          !           SO2 + O3   -> SO4
          !-----------------------------------------------

          !------------------------------------------------------------------------
          !       ... S(IV) (HSO3) + H2O2
          !------------------------------------------------------------------------
          rah2o2 = 8.e4_r8 * EXP( -3650._r8*work1(i) )  &
               / (.1_r8 + xph(i,k))

          !------------------------------------------------------------------------
          !        ... S(IV)+ O3
          !------------------------------------------------------------------------
          rao3   = 4.39e11_r8 * EXP(-4131._r8/tz)  &
               + 2.56e3_r8  * EXP(-996._r8 /tz) /xph(i,k)

          !-----------------------------------------------------------------
          !       ... Prediction after aqueous phase
          !       so4
          !       When Cloud is present 
          !   
          !       S(IV) + H2O2 = S(VI)
          !       S(IV) + O3   = S(VI)
          !
          !       reference:
          !           (1) Seinfeld
          !           (2) Benkovitz
          !-----------------------------------------------------------------
          
          !............................
          !       S(IV) + H2O2 = S(VI)
          !............................
          
          IF (XL .ge. 1.e-8_r8) THEN    !! WHEN CLOUD IS PRESENTED          

             if (cloud_borne) then
                patm_x = patm
             else
                patm_x = 1._r8
             endif

             if (modal_aerosols) then

                pso4 = rah2o2 * 7.4e4_r8*EXP(6621._r8*work1(i)) * h2o2g * patm_x &
                     * 1.23_r8 *EXP(3120._r8*work1(i)) * so2g * patm_x
             else
                pso4 = rah2o2 * heh2o2(i,k) * h2o2g * patm_x  &
                     * heso2(i,k)  * so2g  * patm_x    ! [M/s]

             endif

             pso4 = pso4 & ! [M/s] = [mole/L(w)/s]
                  * xl & ! [mole/L(a)/s]
                  / const0 & ! [/L(a)/s]
                  / xhnm(i,k)


             ccc = pso4*dtime
             ccc = max(ccc, 1.e-30_r8)

             xso4_init(i,k)=xso4(i,k)

             IF (xh2o2(i,k) .gt. xso2(i,k)) THEN
                if (ccc .gt. xso2(i,k)) then
                   xso4(i,k)=xso4(i,k)+xso2(i,k)
                   if (cloud_borne) then
                      xh2o2(i,k)=xh2o2(i,k)-xso2(i,k)
                      xso2(i,k)=1.e-20_r8
                   else       ! ???? bug ????
                      xso2(i,k)=1.e-20_r8
                      xh2o2(i,k)=xh2o2(i,k)-xso2(i,k)
                   endif
                else
                   xso4(i,k)  = xso4(i,k)  + ccc
                   xh2o2(i,k) = xh2o2(i,k) - ccc
                   xso2(i,k)  = xso2(i,k)  - ccc
                end if

             ELSE
                if (ccc  .gt. xh2o2(i,k)) then
                   xso4(i,k)=xso4(i,k)+xh2o2(i,k)
                   xso2(i,k)=xso2(i,k)-xh2o2(i,k)
                   xh2o2(i,k)=1.e-20_r8
                else
                   xso4(i,k)  = xso4(i,k)  + ccc
                   xh2o2(i,k) = xh2o2(i,k) - ccc
                   xso2(i,k)  = xso2(i,k)  - ccc
                end if
             END IF
             
             if (modal_aerosols) then
                xdelso4hp(i,k)  =  xso4(i,k) - xso4_init(i,k)
             endif
             !...........................
             !       S(IV) + O3 = S(VI)
             !...........................

             pso4 = rao3 * heo3(i,k)*o3g*patm_x * heso2(i,k)*so2g*patm_x  ! [M/s]

             pso4 = pso4        &                                ! [M/s] =  [mole/L(w)/s]
                  * xl          &                                ! [mole/L(a)/s]
                  / const0      &                                ! [/L(a)/s]
                  / xhnm(i,k)                                    ! [mixing ratio/s]
             
             ccc = pso4*dtime
             ccc = max(ccc, 1.e-30_r8)

             xso4_init(i,k)=xso4(i,k)

             if (ccc .gt. xso2(i,k)) then
                xso4(i,k) = xso4(i,k) + xso2(i,k)
                xso2(i,k) = 1.e-20_r8
             else
                xso4(i,k) = xso4(i,k) + ccc
                xso2(i,k) = xso2(i,k) - ccc
             end if

          END IF !! WHEN CLOUD IS PRESENTED

       end do col_loop1
    end do ver_loop1
    else
       call setsox_aqchem_predict_codon_wrap(ncol, dtime, const0, kh0, kh1, kh2, kh3, Ra, xkw, &
            press, tfld, xhnm, cldconc%xlwc, xph, xho2, xhno3, xno3, xh2o2, xso2, xo3, &
            xnh3, xnh4, xso4, xso4_init, xdelso4hp, hno3g, nh3g, hehno3, heh2o2, heso2, &
            henh3, heo3)
    end if

    call sox_cldaero_update( &
         ncol, lchnk, loffset, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, cldconc%xlwc, &
         xdelso4hp, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, xnh4c,  xno3c, xmsa, xso2, xh2o2, qcw, qin )
    
    call setsox_xph_lwc_diag_select_impl()
    if (setsox_xph_lwc_diag_use_native_impl) then
       xphlwc(:,:) = 0._r8
       do k = 1, pver
          do i = 1, ncol
             if (cldfrc(i,k)>=1.e-5_r8 .and. lwc(i,k)>=1.e-8_r8) then
                xphlwc(i,k) = -1._r8*log10(xph(i,k)) * lwc(i,k)
             endif
          end do
       end do
    else
       call setsox_xph_lwc_diag_codon_wrap(ncol, cldfrc, lwc, xph, xphlwc)
    end if
    call outfld( 'XPH_LWC', xphlwc(:ncol,:), ncol , lchnk )

    call sox_cldaero_destroy_obj(cldconc)

  end subroutine SETSOX

end module MO_SETSOX
