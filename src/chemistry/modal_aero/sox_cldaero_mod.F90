!----------------------------------------------------------------------------------
! Modal aerosol implementation
!----------------------------------------------------------------------------------
module sox_cldaero_mod

  use shr_kind_mod,    only : r8 => shr_kind_r8
  use cam_abortutils,  only : endrun
  use ppgrid,          only : pcols, pver
  use mo_chem_utls,    only : get_spc_ndx
  use cldaero_mod,     only : cldaero_conc_t, cldaero_allocate, cldaero_deallocate
  use modal_aero_data, only : ntot_amode, modeptr_accum, lptr_so4_cw_amode, lptr_msa_cw_amode
  use modal_aero_data, only : numptrcw_amode, lptr_nh4_cw_amode
  use modal_aero_data, only : cnst_name_cw, specmw_so4_amode
  use cam_history,     only : outfld
  use cam_history,     only : addfld, add_default, phys_decomp
  use chem_mods,       only : adv_mass
  use physconst,       only : gravit
  use phys_control,    only : phys_getopts
  use cldaero_mod,     only : cldaero_uptakerate
  use chem_mods,       only : gas_pcnst
  use mo_constants,    only : pi
  use iso_c_binding,   only : c_int64_t

  implicit none
  private

  public :: sox_cldaero_init
  public :: sox_cldaero_create_obj
  public :: sox_cldaero_update
  public :: sox_cldaero_finalize
  public :: sox_cldaero_destroy_obj

  integer :: id_msa, id_h2so4, id_so2, id_h2o2, id_nh3

  real(r8), parameter :: small_value = 1.e-20_r8
  logical :: sox_cldaero_update_core_use_native_impl = .false.
  logical :: sox_cldaero_update_core_impl_selected = .false.
  logical :: sox_cldaero_update_core_proof_written = .false.
  logical :: sox_cldaero_update_core_wrap_proof_written = .false.
  logical :: sox_cldaero_finalize_wrap_proof_written = .false.
  logical :: sox_cldaero_destroy_obj_codon_logged = .false.
  logical :: sox_cldaero_init_proof_written = .false.
  logical :: sox_cldaero_create_obj_proof_written = .false.

  interface
     function sox_cldaero_init_active_codon(active_c) result(out_c) bind(c, name="sox_cldaero_init_active_codon")
       use iso_c_binding, only : c_int64_t
       integer(c_int64_t), value :: active_c
       integer(c_int64_t) :: out_c
     end function sox_cldaero_init_active_codon
  end interface

contains

  subroutine sox_cldaero_append_impl_proof(proof_line)

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

  end subroutine sox_cldaero_append_impl_proof

  subroutine sox_cldaero_update_core_select_impl()

    use cam_logfile, only : iulog
    use spmd_utils,  only : masterproc

    implicit none

    character(len=48) :: impl_name
    character(len=160) :: proof_line
    integer :: status, n, i, code

    if (sox_cldaero_update_core_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('SOX_CLDAERO_UPDATE_CORE_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       sox_cldaero_update_core_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       sox_cldaero_update_core_use_native_impl = .false.
    end if

    sox_cldaero_update_core_impl_selected = .true.

    if (masterproc) then
       if (sox_cldaero_update_core_use_native_impl) then
          proof_line = 'sox_cldaero_update_core selector entered implementation = native'
       else
          proof_line = 'sox_cldaero_update_core selector entered implementation = codon'
       end if
       write(iulog,'(A)') trim(proof_line)
       if (.not. sox_cldaero_update_core_proof_written) then
          call sox_cldaero_append_impl_proof(trim(proof_line))
          sox_cldaero_update_core_proof_written = .true.
       end if
       call flush(iulog)
    end if

  end subroutine sox_cldaero_update_core_select_impl

  subroutine sox_cldaero_update_core_codon_wrap(ncol, loffset, dtime, press, tfld, cldnum, cldfrc, &
       cfact, xlwc, delso4_hprxn, xh2so4, xso4, xso4_init, nh3g, xnh3, xnh4c, xmsa, &
       xso2, xh2o2, qcw, qin, dqdt_aqso4, dqdt_aqh2so4, dqdt_aqhprxn, dqdt_aqo3rxn, &
       faqgain_msa, faqgain_so4, qnum_c)

    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    use cam_logfile,   only : iulog
    use spmd_utils,    only : masterproc

    implicit none

    integer, intent(in) :: ncol, loffset
    real(r8), intent(in) :: dtime
    real(r8), target, intent(in) :: press(:,:), tfld(:,:), cldnum(:,:)
    real(r8), target, intent(in) :: cldfrc(:,:), cfact(ncol,pver), xlwc(:,:)
    real(r8), target, intent(in) :: delso4_hprxn(ncol,pver), xh2so4(ncol,pver), xso4(ncol,pver)
    real(r8), target, intent(in) :: xso4_init(ncol,pver), nh3g(ncol,pver), xnh3(ncol,pver)
    real(r8), target, intent(in) :: xnh4c(:,:), xmsa(ncol,pver), xso2(ncol,pver), xh2o2(ncol,pver)
    real(r8), target, intent(inout) :: qcw(ncol,pver,gas_pcnst), qin(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: dqdt_aqso4(ncol,pver,gas_pcnst), dqdt_aqh2so4(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: dqdt_aqhprxn(ncol,pver), dqdt_aqo3rxn(ncol,pver)
    real(r8), target, intent(inout) :: faqgain_msa(ntot_amode), faqgain_so4(ntot_amode), qnum_c(ntot_amode)

    integer(c_int64_t), target :: numptrcw_amode_c(ntot_amode), lptr_so4_cw_amode_c(ntot_amode)
    integer(c_int64_t), target :: lptr_msa_cw_amode_c(ntot_amode), lptr_nh4_cw_amode_c(ntot_amode)
    character(len=96) :: proof_line

    interface
       subroutine sox_cldaero_update_core_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, ntot_amode_c, &
            loffset_c, id_msa_c, id_h2so4_c, id_so2_c, id_h2o2_c, id_nh3_c, modeptr_accum_c, &
            dtime_c, pi_c, cldfrc_p, xlwc_p, cldnum_p, cfact_p, tfld_p, press_p, delso4_hprxn_p, &
            xh2so4_p, xso4_p, xso4_init_p, nh3g_p, xnh3_p, xnh4c_p, xmsa_p, xso2_p, xh2o2_p, &
            qcw_p, qin_p, dqdt_aqso4_p, dqdt_aqh2so4_p, dqdt_aqhprxn_p, dqdt_aqo3rxn_p, &
            faqgain_msa_p, faqgain_so4_p, qnum_c_p, numptrcw_amode_p, lptr_so4_cw_amode_p, &
            lptr_msa_cw_amode_p, lptr_nh4_cw_amode_p) bind(c, name="sox_cldaero_update_core_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c, ntot_amode_c
         integer(c_int64_t), value :: loffset_c, id_msa_c, id_h2so4_c, id_so2_c, id_h2o2_c, id_nh3_c
         integer(c_int64_t), value :: modeptr_accum_c
         real(c_double), value :: dtime_c, pi_c
         type(c_ptr), value :: cldfrc_p, xlwc_p, cldnum_p, cfact_p, tfld_p, press_p
         type(c_ptr), value :: delso4_hprxn_p, xh2so4_p, xso4_p, xso4_init_p, nh3g_p, xnh3_p
         type(c_ptr), value :: xnh4c_p, xmsa_p, xso2_p, xh2o2_p, qcw_p, qin_p, dqdt_aqso4_p
         type(c_ptr), value :: dqdt_aqh2so4_p, dqdt_aqhprxn_p, dqdt_aqo3rxn_p
         type(c_ptr), value :: faqgain_msa_p, faqgain_so4_p, qnum_c_p, numptrcw_amode_p
         type(c_ptr), value :: lptr_so4_cw_amode_p, lptr_msa_cw_amode_p, lptr_nh4_cw_amode_p
       end subroutine sox_cldaero_update_core_codon
    end interface

    if (masterproc .and. .not. sox_cldaero_update_core_wrap_proof_written) then
       proof_line = 'sox_cldaero_update_core_codon_wrap entered'
       write(iulog,'(A)') trim(proof_line)
       call sox_cldaero_append_impl_proof(trim(proof_line))
       sox_cldaero_update_core_wrap_proof_written = .true.
       call flush(iulog)
    end if

    numptrcw_amode_c(:) = int(numptrcw_amode(:), c_int64_t)
    lptr_so4_cw_amode_c(:) = int(lptr_so4_cw_amode(:), c_int64_t)
    lptr_msa_cw_amode_c(:) = int(lptr_msa_cw_amode(:), c_int64_t)
    lptr_nh4_cw_amode_c(:) = int(lptr_nh4_cw_amode(:), c_int64_t)

    call sox_cldaero_update_core_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         int(ntot_amode, c_int64_t), int(loffset, c_int64_t), int(id_msa, c_int64_t), &
         int(id_h2so4, c_int64_t), int(id_so2, c_int64_t), int(id_h2o2, c_int64_t), &
         int(id_nh3, c_int64_t), int(modeptr_accum, c_int64_t), real(dtime, c_double), real(pi, c_double), &
         c_loc(cldfrc), c_loc(xlwc), c_loc(cldnum), c_loc(cfact), c_loc(tfld), c_loc(press), &
         c_loc(delso4_hprxn), c_loc(xh2so4), c_loc(xso4), c_loc(xso4_init), c_loc(nh3g), &
         c_loc(xnh3), c_loc(xnh4c), c_loc(xmsa), c_loc(xso2), c_loc(xh2o2), c_loc(qcw), c_loc(qin), &
         c_loc(dqdt_aqso4), c_loc(dqdt_aqh2so4), c_loc(dqdt_aqhprxn), c_loc(dqdt_aqo3rxn), &
         c_loc(faqgain_msa), c_loc(faqgain_so4), c_loc(qnum_c), c_loc(numptrcw_amode_c), &
         c_loc(lptr_so4_cw_amode_c), c_loc(lptr_msa_cw_amode_c), c_loc(lptr_nh4_cw_amode_c) )

  end subroutine sox_cldaero_update_core_codon_wrap

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

  subroutine sox_cldaero_init
    use cam_logfile, only : iulog
    use spmd_utils,  only : masterproc

    integer :: l, m
    logical :: history_aerosol      ! Output the MAM aerosol tendencies
    integer(c_int64_t) :: init_active_c

    init_active_c = sox_cldaero_init_active_codon(1_c_int64_t)
    if (init_active_c == 0_c_int64_t) return

    if (masterproc .and. .not. sox_cldaero_init_proof_written) then
       write(iulog,'(A)') 'sox_cldaero_init direct = codon; init active-policy direct; species/history native CAM API islands'
       call sox_cldaero_append_impl_proof('sox_cldaero_init direct = codon; init active-policy direct; ' // &
            'species/history native CAM API islands')
       sox_cldaero_init_proof_written = .true.
       call flush(iulog)
    end if

    id_msa = get_spc_ndx( 'MSA' )
    id_h2so4 = get_spc_ndx( 'H2SO4' )
    id_so2 = get_spc_ndx( 'SO2' )
    id_h2o2 = get_spc_ndx( 'H2O2' )
    id_nh3 = get_spc_ndx( 'NH3' )

    if (id_h2so4<1 .or. id_so2<1 .or. id_h2o2<1) then
      call endrun('sox_cldaero_init:MAM mech does not include necessary species' &
                  //' -- should not invoke sox_cldaero_mod ')
    endif

    call phys_getopts( history_aerosol_out        = history_aerosol   )
    !
    !   add to history
    !
    do m = 1, ntot_amode

       l = lptr_so4_cw_amode(m)
       if (l > 0) then
          call addfld (&
               trim(cnst_name_cw(l))//'AQSO4','kg/m2/s ',1,  'A', &
               trim(cnst_name_cw(l))//' aqueous phase chemistry',phys_decomp)
          call addfld (&
               trim(cnst_name_cw(l))//'AQH2SO4','kg/m2/s ',1,  'A', &
               trim(cnst_name_cw(l))//' aqueous phase chemistry',phys_decomp)
          if ( history_aerosol ) then 
             call add_default (trim(cnst_name_cw(l))//'AQSO4', 1, ' ')
             call add_default (trim(cnst_name_cw(l))//'AQH2SO4', 1, ' ')
          endif
       end if

    end do

    call addfld ('AQSO4_H2O2','kg/m2/s ',1,  'A', &
         'SO4 aqueous phase chemistry due to H2O2',phys_decomp)
    call addfld ('AQSO4_O3','kg/m2/s ',1,  'A', &
         'SO4 aqueous phase chemistry due to O3',phys_decomp)

    if ( history_aerosol ) then    
       call add_default ('AQSO4_H2O2', 1, ' ')
       call add_default ('AQSO4_O3', 1, ' ')    
    endif
  
  end subroutine sox_cldaero_init

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
  function sox_cldaero_create_obj(cldfrc, qcw, lwc, cfact, ncol, loffset) result( conc_obj )
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use cam_logfile,   only : iulog
    use spmd_utils,    only : masterproc
    
    real(r8), target, intent(in) :: cldfrc(:,:)
    real(r8), target, intent(in) :: qcw(:,:,:)
    real(r8), target, intent(in) :: lwc(:,:)
    real(r8), target, intent(in) :: cfact(:,:)
    integer,  intent(in) :: ncol
    integer,  intent(in) :: loffset

    type(cldaero_conc_t), pointer :: conc_obj

    integer(c_int64_t), target :: lptr_so4_cw_amode_c(ntot_amode)
    integer(c_int64_t), target :: lptr_nh4_cw_amode_c(ntot_amode)
    character(len=160) :: proof_line
    integer :: n, i, k
    integer :: id_so4_1a, id_so4_2a, id_so4_3a, id_so4_4a, id_so4_5a, id_so4_6a
    integer :: id_nh4_1a, id_nh4_2a, id_nh4_3a, id_nh4_4a, id_nh4_5a, id_nh4_6a
    logical :: use_native_object

    interface
       subroutine sox_cldaero_create_obj_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, ntot_amode_c, loffset_c, &
            cldfrc_p, qcw_p, lwc_p, cfact_p, so4c_p, nh4c_p, no3c_p, xlwc_p, so4_fact_p, &
            lptr_so4_cw_amode_p, lptr_nh4_cw_amode_p) bind(c, name="sox_cldaero_create_obj_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c, ntot_amode_c, loffset_c
         type(c_ptr), value :: cldfrc_p, qcw_p, lwc_p, cfact_p, so4c_p, nh4c_p, no3c_p, xlwc_p, so4_fact_p
         type(c_ptr), value :: lptr_so4_cw_amode_p, lptr_nh4_cw_amode_p
       end subroutine sox_cldaero_create_obj_codon
    end interface

    conc_obj => cldaero_allocate()

    use_native_object = sox_cldaero_object_use_native()
    if (use_native_object) then
       do k = 1, pver
          do i = 1, ncol
             if (cldfrc(i,k) > 0._r8) then
                conc_obj%xlwc(i,k) = lwc(i,k) * cfact(i,k)
                conc_obj%xlwc(i,k) = conc_obj%xlwc(i,k) / cldfrc(i,k)
             else
                conc_obj%xlwc(i,k) = 0._r8
             end if
             conc_obj%no3c(i,k) = 0._r8
          end do
       end do

       if (ntot_amode == 7) then
          id_so4_1a = lptr_so4_cw_amode(1) - loffset
          id_so4_2a = lptr_so4_cw_amode(2) - loffset
          id_so4_3a = lptr_so4_cw_amode(4) - loffset
          id_so4_4a = lptr_so4_cw_amode(5) - loffset
          id_so4_5a = lptr_so4_cw_amode(6) - loffset
          id_so4_6a = lptr_so4_cw_amode(7) - loffset
          id_nh4_1a = lptr_nh4_cw_amode(1) - loffset
          id_nh4_2a = lptr_nh4_cw_amode(2) - loffset
          id_nh4_3a = lptr_nh4_cw_amode(4) - loffset
          id_nh4_4a = lptr_nh4_cw_amode(5) - loffset
          id_nh4_5a = lptr_nh4_cw_amode(6) - loffset
          id_nh4_6a = lptr_nh4_cw_amode(7) - loffset
          do k = 1, pver
             do i = 1, ncol
                conc_obj%so4c(i,k) = qcw(i,k,id_so4_1a) + qcw(i,k,id_so4_2a) + &
                     qcw(i,k,id_so4_3a) + qcw(i,k,id_so4_4a) + qcw(i,k,id_so4_5a) + qcw(i,k,id_so4_6a)
                conc_obj%nh4c(i,k) = qcw(i,k,id_nh4_1a) + qcw(i,k,id_nh4_2a) + &
                     qcw(i,k,id_nh4_3a) + qcw(i,k,id_nh4_4a) + qcw(i,k,id_nh4_5a) + qcw(i,k,id_nh4_6a)
             end do
          end do
       else
          id_so4_1a = lptr_so4_cw_amode(1) - loffset
          id_so4_2a = lptr_so4_cw_amode(2) - loffset
          id_so4_3a = lptr_so4_cw_amode(3) - loffset
          do k = 1, pver
             do i = 1, ncol
                conc_obj%so4c(i,k) = qcw(i,k,id_so4_1a) + qcw(i,k,id_so4_2a) + qcw(i,k,id_so4_3a)
                conc_obj%nh4c(i,k) = 0._r8
             end do
          end do
          conc_obj%so4_fact = 1._r8
       end if

       if (masterproc .and. .not. sox_cldaero_create_obj_proof_written) then
          proof_line = 'sox_cldaero_create_obj implementation = native'
          write(iulog,'(A)') trim(proof_line)
          call sox_cldaero_append_impl_proof(trim(proof_line))
          sox_cldaero_create_obj_proof_written = .true.
          call flush(iulog)
       end if
       return
    end if

    do n = 1, ntot_amode
       lptr_so4_cw_amode_c(n) = int(lptr_so4_cw_amode(n), c_int64_t)
       lptr_nh4_cw_amode_c(n) = int(lptr_nh4_cw_amode(n), c_int64_t)
    end do

    call sox_cldaero_create_obj_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         int(ntot_amode, c_int64_t), int(loffset, c_int64_t), c_loc(cldfrc), c_loc(qcw), c_loc(lwc), c_loc(cfact), &
         c_loc(conc_obj%so4c(1,1)), c_loc(conc_obj%nh4c(1,1)), c_loc(conc_obj%no3c(1,1)), c_loc(conc_obj%xlwc(1,1)), &
         c_loc(conc_obj%so4_fact), c_loc(lptr_so4_cw_amode_c), c_loc(lptr_nh4_cw_amode_c) )

    if (masterproc .and. .not. sox_cldaero_create_obj_proof_written) then
       proof_line = 'sox_cldaero_create_obj implementation = codon; allocation/pointer boundary = native'
       write(iulog,'(A)') trim(proof_line)
       call sox_cldaero_append_impl_proof(trim(proof_line))
       sox_cldaero_create_obj_proof_written = .true.
       call flush(iulog)
    end if

  end function sox_cldaero_create_obj

!----------------------------------------------------------------------------------
! Update the mixing ratios
!----------------------------------------------------------------------------------
  subroutine sox_cldaero_update( &
       ncol, lchnk, loffset, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, xlwc, &
       delso4_hprxn, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, xnh4c,  xno3c, xmsa, xso2, xh2o2, qcw, qin )

    ! args 

    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk ! chunk id
    integer,  intent(in) :: loffset

    real(r8), intent(in) :: dtime ! time step (sec)

    real(r8), intent(in) :: mbar(:,:) ! mean wet atmospheric mass ( amu )
    real(r8), intent(in) :: pdel(:,:) 
    real(r8), intent(in) :: press(:,:)
    real(r8), intent(in) :: tfld(:,:)

    real(r8), intent(in) :: cldnum(:,:)
    real(r8), intent(in) :: cldfrc(:,:)
    real(r8), intent(in) :: cfact(:,:)
    real(r8), intent(in) :: xlwc(:,:)

    real(r8), intent(in) :: delso4_hprxn(:,:)
    real(r8), intent(in) :: xh2so4(:,:)
    real(r8), intent(in) :: xso4(:,:)
    real(r8), intent(in) :: xso4_init(:,:)
    real(r8), intent(in) :: nh3g(:,:)
    real(r8), intent(in) :: hno3g(:,:)
    real(r8), intent(in) :: xnh3(:,:)
    real(r8), intent(in) :: xhno3(:,:)
    real(r8), intent(in) :: xnh4c(:,:)
    real(r8), intent(in) :: xmsa(:,:)
    real(r8), intent(in) :: xso2(:,:)
    real(r8), intent(in) :: xh2o2(:,:)
    real(r8), intent(in) :: xno3c(:,:)

    real(r8), intent(inout) :: qcw(:,:,:) ! cloud-borne aerosol (vmr)
    real(r8), intent(inout) :: qin(:,:,:) ! xported species ( vmr )

    ! local vars ...

    real(r8) :: dqdt_aqso4(ncol,pver,gas_pcnst), &
         dqdt_aqh2so4(ncol,pver,gas_pcnst), &
         dqdt_aqhprxn(ncol,pver), dqdt_aqo3rxn(ncol,pver)

    real(r8) :: faqgain_msa(ntot_amode), faqgain_so4(ntot_amode), qnum_c(ntot_amode)

    real(r8) :: delso4_o3rxn, &
         dso4dt_aqrxn, dso4dt_hprxn, &
         dso4dt_gasuptk, dmsadt_gasuptk, &
         dmsadt_gasuptk_tomsa, dmsadt_gasuptk_toso4, &
         dqdt_aq, dqdt_wr, dqdt

    real(r8) :: fwetrem, sumf, uptkrate
    real(r8) :: delnh3, delnh4

    integer :: l, n, m
    integer :: ntot_msa_c

    integer :: i,k
    real(r8) :: xl

    call sox_cldaero_update_core_select_impl()
    if (sox_cldaero_update_core_use_native_impl) then
    ! make sure dqdt is zero initially, for budgets
    dqdt_aqso4(:,:,:) = 0.0_r8
    dqdt_aqh2so4(:,:,:) = 0.0_r8
    dqdt_aqhprxn(:,:) = 0.0_r8
    dqdt_aqo3rxn(:,:) = 0.0_r8

    lev_loop: do k = 1,pver
       col_loop: do i = 1,ncol
          cloud: if (cldfrc(i,k) >= 1.0e-5_r8) then
             xl = xlwc(i,k) ! / cldfrc(i,k)

             IF (XL .ge. 1.e-8_r8) THEN !! WHEN CLOUD IS PRESENTED

                delso4_o3rxn = xso4(i,k) - xso4_init(i,k)

                if (id_nh3>0) then
                   delnh3 = nh3g(i,k) - xnh3(i,k)
                   delnh4 = - delnh3
                endif

                !-------------------------------------------------------------------------
                ! compute factors for partitioning aerosol mass gains among modes
                ! the factors are proportional to the activated particle MR for each
                ! mode, which is the MR of cloud drops "associated with" the mode
                ! thus we are assuming the cloud drop size is independent of the
                ! associated aerosol mode properties (i.e., drops associated with
                ! Aitken and coarse sea-salt particles are same size)
                !
                ! qnum_c(n) = activated particle number MR for mode n (these are just
                ! used for partitioning among modes, so don't need to divide by cldfrc)

                do n = 1, ntot_amode
                   qnum_c(n) = 0.0_r8
                   l = numptrcw_amode(n) - loffset
                   if (l > 0) qnum_c(n) = max( 0.0_r8, qcw(i,k,l) )
                end do

                ! force qnum_c(n) to be positive for n=modeptr_accum or n=1
                n = modeptr_accum
                if (n <= 0) n = 1
                qnum_c(n) = max( 1.0e-10_r8, qnum_c(n) )

                ! faqgain_so4(n) = fraction of total so4_c gain going to mode n
                ! these are proportional to the activated particle MR for each mode
                sumf = 0.0_r8
                do n = 1, ntot_amode
                   faqgain_so4(n) = 0.0_r8
                   if (lptr_so4_cw_amode(n) > 0) then
                      faqgain_so4(n) = qnum_c(n)
                      sumf = sumf + faqgain_so4(n)
                   end if
                end do

                if (sumf > 0.0_r8) then
                   do n = 1, ntot_amode
                      faqgain_so4(n) = faqgain_so4(n) / sumf
                   end do
                end if
                ! at this point (sumf <= 0.0) only when all the faqgain_so4 are zero

                ! faqgain_msa(n) = fraction of total msa_c gain going to mode n
                ntot_msa_c = 0
                sumf = 0.0_r8
                do n = 1, ntot_amode
                   faqgain_msa(n) = 0.0_r8
                   if (lptr_msa_cw_amode(n) > 0) then
                      faqgain_msa(n) = qnum_c(n)
                      ntot_msa_c = ntot_msa_c + 1
                   end if
                   sumf = sumf + faqgain_msa(n)
                end do

                if (sumf > 0.0_r8) then
                   do n = 1, ntot_amode
                      faqgain_msa(n) = faqgain_msa(n) / sumf
                   end do
                end if
                ! at this point (sumf <= 0.0) only when all the faqgain_msa are zero

                uptkrate = cldaero_uptakerate( xl, cldnum(i,k), cfact(i,k), cldfrc(i,k), tfld(i,k),  press(i,k) )
                ! average uptake rate over dtime
                uptkrate = (1.0_r8 - exp(-min(100._r8,dtime*uptkrate))) / dtime

                ! dso4dt_gasuptk = so4_c tendency from h2so4 gas uptake (mol/mol/s)
                ! dmsadt_gasuptk = msa_c tendency from msa gas uptake (mol/mol/s)
                dso4dt_gasuptk = xh2so4(i,k) * uptkrate
                if (id_msa > 0) then
                   dmsadt_gasuptk = xmsa(i,k) * uptkrate
                else
                   dmsadt_gasuptk = 0.0_r8
                end if

                ! if no modes have msa aerosol, then "rename" scavenged msa gas to so4
                dmsadt_gasuptk_toso4 = 0.0_r8
                dmsadt_gasuptk_tomsa = dmsadt_gasuptk
                if (ntot_msa_c == 0) then
                   dmsadt_gasuptk_tomsa = 0.0_r8
                   dmsadt_gasuptk_toso4 = dmsadt_gasuptk
                end if

                !-----------------------------------------------------------------------
                ! now compute TMR tendencies
                ! this includes the above aqueous so2 chemistry AND
                ! the uptake of highly soluble aerosol precursor gases (h2so4, msa, ...)
                ! AND the wetremoval of dissolved, unreacted so2 and h2o2

                dso4dt_aqrxn = (delso4_o3rxn + delso4_hprxn(i,k)) / dtime
                dso4dt_hprxn = delso4_hprxn(i,k) / dtime

                ! fwetrem = fraction of in-cloud-water material that is wet removed
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*clwlrat(i,k)))) )
                fwetrem = 0.0_r8 ! don't have so4 & msa wet removal here

                ! compute TMR tendencies for so4 and msa aerosol-in-cloud-water
                do n = 1, ntot_amode
                   l = lptr_so4_cw_amode(n) - loffset
                   if (l > 0) then
                      dqdt_aqso4(i,k,l) = faqgain_so4(n)*dso4dt_aqrxn*cldfrc(i,k)
                      dqdt_aqh2so4(i,k,l) = faqgain_so4(n)* &
                           (dso4dt_gasuptk + dmsadt_gasuptk_toso4)*cldfrc(i,k)
                      dqdt_aq = dqdt_aqso4(i,k,l) + dqdt_aqh2so4(i,k,l)
                      dqdt_wr = -fwetrem*dqdt_aq
                      dqdt= dqdt_aq + dqdt_wr
                      qcw(i,k,l) = qcw(i,k,l) + dqdt*dtime
                   end if

                   l = lptr_msa_cw_amode(n) - loffset
                   if (l > 0) then
                      dqdt_aq = faqgain_msa(n)*dmsadt_gasuptk_tomsa*cldfrc(i,k)
                      dqdt_wr = -fwetrem*dqdt_aq
                      dqdt = dqdt_aq + dqdt_wr
                      qcw(i,k,l) = qcw(i,k,l) + dqdt*dtime
                   end if

                   l = lptr_nh4_cw_amode(n) - loffset
                   if (l > 0) then
                      if (delnh4 > 0.0_r8) then
                         dqdt_aq = faqgain_so4(n)*delnh4/dtime*cldfrc(i,k)
                         dqdt = dqdt_aq
                         qcw(i,k,l) = qcw(i,k,l) + dqdt*dtime
                      else
                         dqdt = (qcw(i,k,l)/max(xnh4c(i,k),1.0e-35_r8)) &
                              *delnh4/dtime*cldfrc(i,k)
                         qcw(i,k,l) = qcw(i,k,l) + dqdt*dtime
                      endif
                   end if
                end do

                ! For gas species, tendency includes
                ! reactive uptake to cloud water that essentially transforms the gas to
                ! a different species. Wet removal associated with this is applied
                ! to the "new" species (e.g., so4_c) rather than to the gas.
                ! wet removal of the unreacted gas that is dissolved in cloud water.
                ! Need to multiply both these parts by cldfrc

                ! h2so4 (g) & msa (g)
                qin(i,k,id_h2so4) = qin(i,k,id_h2so4) - dso4dt_gasuptk * dtime * cldfrc(i,k)
                if (id_msa > 0) qin(i,k,id_msa) = qin(i,k,id_msa) - dmsadt_gasuptk * dtime * cldfrc(i,k)

                ! so2 -- the first order loss rate for so2 is frso2_c*clwlrat(i,k)
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*frso2_c*clwlrat(i,k)))) )
                fwetrem = 0.0_r8 ! don't include so2 wet removal here

                dqdt_wr = -fwetrem*xso2(i,k)/dtime*cldfrc(i,k)
                dqdt_aq = -dso4dt_aqrxn*cldfrc(i,k)
                dqdt = dqdt_aq + dqdt_wr
                qin(i,k,id_so2) = qin(i,k,id_so2) + dqdt * dtime

                ! h2o2 -- the first order loss rate for h2o2 is frh2o2_c*clwlrat(i,k)
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*frh2o2_c*clwlrat(i,k)))) )
                fwetrem = 0.0_r8 ! don't include h2o2 wet removal here

                dqdt_wr = -fwetrem*xh2o2(i,k)/dtime*cldfrc(i,k)
                dqdt_aq = -dso4dt_hprxn*cldfrc(i,k)
                dqdt = dqdt_aq + dqdt_wr
                qin(i,k,id_h2o2) = qin(i,k,id_h2o2) + dqdt * dtime

                ! NH3
                if (id_nh3>0) then
                   dqdt_aq = delnh3/dtime*cldfrc(i,k)
                   dqdt = dqdt_aq
                   qin(i,k,id_nh3) = qin(i,k,id_nh3) + dqdt * dtime
                endif

                ! for SO4 from H2O2/O3 budgets
                dqdt_aqhprxn(i,k) = dso4dt_hprxn*cldfrc(i,k)
                dqdt_aqo3rxn(i,k) = (dso4dt_aqrxn - dso4dt_hprxn)*cldfrc(i,k)

             ENDIF !! WHEN CLOUD IS PRESENTED
          endif cloud
       enddo col_loop
    enddo lev_loop
    else
       call sox_cldaero_update_core_codon_wrap(ncol, loffset, dtime, press, tfld, cldnum, &
            cldfrc, cfact, xlwc, delso4_hprxn, xh2so4, xso4, xso4_init, nh3g, xnh3, &
            xnh4c, xmsa, xso2, xh2o2, qcw, qin, dqdt_aqso4, dqdt_aqh2so4, dqdt_aqhprxn, &
            dqdt_aqo3rxn, faqgain_msa, faqgain_so4, qnum_c)
    end if

    call sox_cldaero_finalize(ncol, lchnk, loffset, mbar, pdel, qcw, qin, dqdt_aqso4, &
         dqdt_aqh2so4, dqdt_aqhprxn, dqdt_aqo3rxn)

  end subroutine sox_cldaero_update

  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------
  subroutine sox_cldaero_finalize(ncol, lchnk, loffset, mbar, pdel, qcw, qin, dqdt_aqso4, &
       dqdt_aqh2so4, dqdt_aqhprxn, dqdt_aqo3rxn)

    integer, intent(in) :: ncol, lchnk, loffset
    real(r8), target, intent(in) :: mbar(:,:), pdel(:,:)
    real(r8), target, intent(inout) :: qcw(:,:,:), qin(:,:,:)
    real(r8), target, intent(in) :: dqdt_aqso4(ncol,pver,gas_pcnst), dqdt_aqh2so4(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: dqdt_aqhprxn(ncol,pver), dqdt_aqo3rxn(ncol,pver)

    real(r8) :: sflx(1:ncol)
    real(r8), target :: sflx_aqso4(ncol,ntot_amode), sflx_aqh2so4(ncol,ntot_amode)
    real(r8), target :: sflx_aqhprxn(ncol), sflx_aqo3rxn(ncol)
    integer :: i, k, l, m, n

    call sox_cldaero_update_core_select_impl()
    if (.not. sox_cldaero_update_core_use_native_impl) then
       call sox_cldaero_finalize_codon_wrap(ncol, loffset, mbar, pdel, qcw, qin, dqdt_aqso4, &
            dqdt_aqh2so4, dqdt_aqhprxn, dqdt_aqo3rxn, sflx_aqso4, sflx_aqh2so4, sflx_aqhprxn, sflx_aqo3rxn)

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
       return
    end if

    !==============================================================
    ! ... Update the mixing ratios
    !==============================================================
    do k = 1,pver

       do n = 1, ntot_amode

          l = lptr_so4_cw_amode(n) - loffset
          if (l > 0) then
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
          end if
          l = lptr_msa_cw_amode(n) - loffset
          if (l > 0) then
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
          end if
          l = lptr_nh4_cw_amode(n) - loffset
          if (l > 0) then
             qcw(:,k,l) = MAX(qcw(:,k,l), small_value )
          end if

       end do

       qin(:,k,id_so2) =  MAX( qin(:,k,id_so2),    small_value )

       if ( id_nh3 > 0 ) then
          qin(:,k,id_nh3) =  MAX( qin(:,k,id_nh3),    small_value )
       endif

    end do

    ! diagnostics

    do n = 1, ntot_amode
       m = lptr_so4_cw_amode(n)
       l = m - loffset
       if (l > 0) then
          sflx(:)=0._r8
          do k=1,pver
             do i=1,ncol
                sflx(i)=sflx(i)+dqdt_aqso4(i,k,l)*adv_mass(l)/mbar(i,k) &
                     *pdel(i,k)/gravit ! kg/m2/s
             enddo
          enddo
          call outfld( trim(cnst_name_cw(m))//'AQSO4', sflx(:ncol), ncol, lchnk)

          sflx(:)=0._r8
          do k=1,pver
             do i=1,ncol
                sflx(i)=sflx(i)+dqdt_aqh2so4(i,k,l)*adv_mass(l)/mbar(i,k) &
                     *pdel(i,k)/gravit ! kg/m2/s
             enddo
          enddo
          call outfld( trim(cnst_name_cw(m))//'AQH2SO4', sflx(:ncol), ncol, lchnk)
       endif
    end do

    sflx(:)=0._r8
    do k=1,pver
       do i=1,ncol
          sflx(i)=sflx(i)+dqdt_aqhprxn(i,k)*specmw_so4_amode/mbar(i,k) &
               *pdel(i,k)/gravit ! kg SO4 /m2/s
       enddo
    enddo
    call outfld( 'AQSO4_H2O2', sflx(:ncol), ncol, lchnk)
    sflx(:)=0._r8
    do k=1,pver
       do i=1,ncol
          sflx(i)=sflx(i)+dqdt_aqo3rxn(i,k)*specmw_so4_amode/mbar(i,k) &
               *pdel(i,k)/gravit ! kg SO4 /m2/s
       enddo
    enddo
    call outfld( 'AQSO4_O3', sflx(:ncol), ncol, lchnk)

  end subroutine sox_cldaero_finalize

  subroutine sox_cldaero_finalize_codon_wrap(ncol, loffset, mbar, pdel, qcw, qin, dqdt_aqso4, &
       dqdt_aqh2so4, dqdt_aqhprxn, dqdt_aqo3rxn, sflx_aqso4, sflx_aqh2so4, sflx_aqhprxn, sflx_aqo3rxn)

    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    use cam_logfile,   only : iulog
    use spmd_utils,    only : masterproc

    implicit none

    integer, intent(in) :: ncol, loffset
    real(r8), target, intent(in) :: mbar(:,:), pdel(:,:)
    real(r8), target, intent(inout) :: qcw(:,:,:), qin(:,:,:)
    real(r8), target, intent(in) :: dqdt_aqso4(ncol,pver,gas_pcnst), dqdt_aqh2so4(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: dqdt_aqhprxn(ncol,pver), dqdt_aqo3rxn(ncol,pver)
    real(r8), target, intent(out) :: sflx_aqso4(ncol,ntot_amode), sflx_aqh2so4(ncol,ntot_amode)
    real(r8), target, intent(out) :: sflx_aqhprxn(ncol), sflx_aqo3rxn(ncol)

    integer(c_int64_t), target :: lptr_so4_cw_amode_c(ntot_amode), lptr_msa_cw_amode_c(ntot_amode)
    integer(c_int64_t), target :: lptr_nh4_cw_amode_c(ntot_amode)
    character(len=160) :: proof_line

    interface
       subroutine sox_cldaero_finalize_codon(ncol_c, pver_c, gas_pcnst_c, ntot_amode_c, loffset_c, &
            id_so2_c, id_nh3_c, small_value_c, specmw_so4_amode_c, gravit_c, mbar_p, pdel_p, qcw_p, &
            qin_p, dqdt_aqso4_p, dqdt_aqh2so4_p, dqdt_aqhprxn_p, dqdt_aqo3rxn_p, sflx_aqso4_p, &
            sflx_aqh2so4_p, sflx_aqhprxn_p, sflx_aqo3rxn_p, adv_mass_p, lptr_so4_cw_amode_p, &
            lptr_msa_cw_amode_p, lptr_nh4_cw_amode_p) bind(c, name="sox_cldaero_finalize_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, ntot_amode_c, loffset_c
         integer(c_int64_t), value :: id_so2_c, id_nh3_c
         real(c_double), value :: small_value_c, specmw_so4_amode_c, gravit_c
         type(c_ptr), value :: mbar_p, pdel_p, qcw_p, qin_p, dqdt_aqso4_p, dqdt_aqh2so4_p
         type(c_ptr), value :: dqdt_aqhprxn_p, dqdt_aqo3rxn_p, sflx_aqso4_p, sflx_aqh2so4_p
         type(c_ptr), value :: sflx_aqhprxn_p, sflx_aqo3rxn_p, adv_mass_p, lptr_so4_cw_amode_p
         type(c_ptr), value :: lptr_msa_cw_amode_p, lptr_nh4_cw_amode_p
       end subroutine sox_cldaero_finalize_codon
    end interface

    if (masterproc .and. .not. sox_cldaero_finalize_wrap_proof_written) then
       proof_line = 'sox_cldaero_finalize_codon_wrap entered (clamp/diagnostic flux sums direct = codon; outfld = native)'
       write(iulog,'(A)') trim(proof_line)
       call sox_cldaero_append_impl_proof(trim(proof_line))
       sox_cldaero_finalize_wrap_proof_written = .true.
       call flush(iulog)
    end if

    lptr_so4_cw_amode_c(:) = int(lptr_so4_cw_amode(:), c_int64_t)
    lptr_msa_cw_amode_c(:) = int(lptr_msa_cw_amode(:), c_int64_t)
    lptr_nh4_cw_amode_c(:) = int(lptr_nh4_cw_amode(:), c_int64_t)

    call sox_cldaero_finalize_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         int(ntot_amode, c_int64_t), int(loffset, c_int64_t), int(id_so2, c_int64_t), &
         int(id_nh3, c_int64_t), real(small_value, c_double), real(specmw_so4_amode, c_double), &
         real(gravit, c_double), c_loc(mbar), c_loc(pdel), c_loc(qcw), c_loc(qin), c_loc(dqdt_aqso4), &
         c_loc(dqdt_aqh2so4), c_loc(dqdt_aqhprxn), c_loc(dqdt_aqo3rxn), c_loc(sflx_aqso4), &
         c_loc(sflx_aqh2so4), c_loc(sflx_aqhprxn), c_loc(sflx_aqo3rxn), c_loc(adv_mass), &
         c_loc(lptr_so4_cw_amode_c), c_loc(lptr_msa_cw_amode_c), c_loc(lptr_nh4_cw_amode_c) )

  end subroutine sox_cldaero_finalize_codon_wrap

  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------
  subroutine sox_cldaero_destroy_obj( conc_obj )
    use iso_c_binding, only : c_int64_t
    use cam_logfile,   only : iulog
    use spmd_utils,    only : masterproc

    type(cldaero_conc_t), pointer :: conc_obj
    integer(c_int64_t) :: active_c
    character(len=80) :: proof_line

    interface
       function sox_cldaero_destroy_obj_codon(stage_c) result(out_c) bind(c, name="sox_cldaero_destroy_obj_codon")
         import :: c_int64_t
         integer(c_int64_t), value :: stage_c
         integer(c_int64_t) :: out_c
       end function sox_cldaero_destroy_obj_codon
    end interface

    if (sox_cldaero_object_use_native()) then
       active_c = 1_c_int64_t
    else
       active_c = sox_cldaero_destroy_obj_codon(1_c_int64_t)
    end if
    if (.not. sox_cldaero_destroy_obj_codon_logged) then
       sox_cldaero_destroy_obj_codon_logged = .true.
       if (sox_cldaero_object_use_native()) then
          proof_line = 'sox_cldaero_destroy_obj implementation = native'
       else
          proof_line = 'sox_cldaero_destroy_obj implementation = codon'
       end if
       if (masterproc) then
          write(iulog,'(A)') trim(proof_line)
          call flush(iulog)
       end if
       call sox_cldaero_append_impl_proof(proof_line)
    end if
    if (active_c == 0_c_int64_t) return

    call cldaero_deallocate( conc_obj )

  end subroutine sox_cldaero_destroy_obj

  logical function sox_cldaero_object_use_native()
    character(len=32) :: impl_name
    integer :: n, status, i, code

    impl_name = 'codon'
    call cam_codon_get_impl('SOX_CLDAERO_OBJECT_IMPL', impl_name, n, status)
    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       sox_cldaero_object_use_native = trim(adjustl(impl_name(:n))) == 'native'
    else
       sox_cldaero_object_use_native = .false.
    end if
  end function sox_cldaero_object_use_native

end module sox_cldaero_mod
