! CAM compatibility adapter for the decoupled modal aerosol nucleation
! process.  The public CAM entry points are unchanged.
module modal_aero_newnuc

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : gas_pcnst
  use ppgrid, only : pcols, pver
  use perf_mod, only : t_startf, t_stopf
  use ap_modal_aero_newnuc_scheme, only : modal_aero_newnuc_sub_run

  implicit none
  private

  integer, parameter :: pcnstxx = gas_pcnst
  integer :: l_h2so4_sv = 0
  integer :: l_nh3_sv = 0
  integer :: lnumait_sv = 0
  integer :: lnh4ait_sv = 0
  integer :: lso4ait_sv = 0

  public :: modal_aero_newnuc_sub, modal_aero_newnuc_init
  ! Read-only process configuration consumed by the parent standalone driver.
  ! These values are still initialized by the CAM facade above.
  public :: l_h2so4_sv, l_nh3_sv, lnumait_sv, lnh4ait_sv, lso4ait_sv

contains

  subroutine modal_aero_newnuc_sub(                              &
       lchnk, ncol, nstep, loffset, deltat, t, pmid, pdel, zm,   &
       pblh, qv, cld, q, del_h2so4_gasprod, del_h2so4_aeruptk)

    use cam_history, only : outfld, fieldname_len
    use chem_mods, only : adv_mass
    use constituents, only : pcnst, cnst_name
    use modal_aero_data, only : ntot_amode, modeptr_aitken,      &
         lptr_nh4_a_amode, dgnumlo_amode, dgnum_amode,         &
         dgnumhi_amode, specdens_so4_amode, specmw_so4_amode,  &
         specmw_nh4_amode
    use physconst, only : gravit, mwdry
    use ref_pres, only : top_lev => clim_modal_aero_top_lev
    use wv_saturation, only : qsat

    integer, intent(in) :: lchnk, ncol, nstep, loffset
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: zm(pcols,pver)
    real(r8), intent(in) :: pblh(pcols)
    real(r8), intent(in) :: qv(pcols,pver)
    real(r8), intent(in) :: cld(ncol,pver)
    real(r8), intent(inout) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: del_h2so4_gasprod(ncol,pver)
    real(r8), intent(in) :: del_h2so4_aeruptk(ncol,pver)

    integer :: i, l, lmz, mait, lptr_nh4_aitken
    real(r8) :: dgnumlo_aitken, dgnum_aitken, dgnumhi_aitken
    real(r8) :: ev_sat(pcols,pver), qv_sat(pcols,pver)
    real(r8) :: qsrflx(pcols,pcnst)
    logical :: dotend(pcnst)
    character(len=fieldname_len+3) :: fieldname

    call t_startf('ap_modal_aero_newnuc_sub_run')

    mait = modeptr_aitken
    lptr_nh4_aitken = 0
    dgnumlo_aitken = 1.0_r8
    dgnum_aitken = 1.0_r8
    dgnumhi_aitken = 1.0_r8
    if ((mait > 0) .and. (mait <= ntot_amode)) then
       lptr_nh4_aitken = lptr_nh4_a_amode(mait)
       dgnumlo_aitken = dgnumlo_amode(mait)
       dgnum_aitken = dgnum_amode(mait)
       dgnumhi_aitken = dgnumhi_amode(mait)
    end if

    call qsat(t(1:ncol,1:pver), pmid(1:ncol,1:pver), &
         ev_sat(1:ncol,1:pver), qv_sat(1:ncol,1:pver))

    call modal_aero_newnuc_sub_run(                              &
         lchnk, ncol, nstep, pcols, pver, pcnstxx, pcnst,       &
         top_lev, loffset, deltat, l_h2so4_sv, l_nh3_sv,        &
         lnumait_sv, lnh4ait_sv, lso4ait_sv, lptr_nh4_aitken,  &
         dgnumlo_aitken, dgnum_aitken, dgnumhi_aitken,         &
         specdens_so4_amode, specmw_so4_amode,                 &
         specmw_nh4_amode, gravit, t, pmid, pdel, zm, pblh,    &
         qv, qv_sat, cld, q, del_h2so4_gasprod,                &
         del_h2so4_aeruptk, qsrflx, dotend)

    do l = loffset+1, pcnst
       lmz = l - loffset
       if (.not. dotend(lmz)) cycle
       do i = 1, ncol
          qsrflx(i,lmz) = qsrflx(i,lmz)*(adv_mass(lmz)/mwdry)
       end do
       fieldname = trim(cnst_name(l)) // '_sfnnuc1'
       call outfld(fieldname, qsrflx(:,lmz), pcols, lchnk)
    end do

    call t_stopf('ap_modal_aero_newnuc_sub_run')

  end subroutine modal_aero_newnuc_sub

  subroutine modal_aero_newnuc_init

    use cam_history, only : addfld, add_default, fieldname_len, phys_decomp
    use constituents, only : pcnst, cnst_get_ind, cnst_name
    use modal_aero_data
    use phys_control, only : phys_getopts
    use spmd_utils, only : masterproc

    integer :: l_h2so4, l_nh3
    integer :: lnumait, lnh4ait, lso4ait
    integer :: l, m, mait
    character(len=fieldname_len) :: tmpname
    character(len=fieldname_len+3) :: fieldname
    character(len=128) :: long_name
    character(len=8) :: unit
    logical :: dotend(pcnst)
    logical :: history_aerosol

    call phys_getopts(history_aerosol_out=history_aerosol)

    l_h2so4_sv = 0
    l_nh3_sv = 0
    lnumait_sv = 0
    lnh4ait_sv = 0
    lso4ait_sv = 0

    call cnst_get_ind('H2SO4', l_h2so4, .false.)
    call cnst_get_ind('NH3', l_nh3, .false.)

    mait = modeptr_aitken
    if (mait > 0) then
       lnumait = numptr_amode(mait)
       lso4ait = lptr_so4_a_amode(mait)
       lnh4ait = lptr_nh4_a_amode(mait)
    end if

    if ((l_h2so4 <= 0) .or. (l_h2so4 > pcnst)) then
       write(*,'(/a/)') '*** modal_aero_newnuc bypass -- l_h2so4 <= 0'
       return
    else if ((lso4ait <= 0) .or. (lso4ait > pcnst)) then
       write(*,'(/a/)') '*** modal_aero_newnuc bypass -- lso4ait <= 0'
       return
    else if ((lnumait <= 0) .or. (lnumait > pcnst)) then
       write(*,'(/a/)') '*** modal_aero_newnuc bypass -- lnumait <= 0'
       return
    else if ((mait <= 0) .or. (mait > ntot_amode)) then
       write(*,'(/a/)') '*** modal_aero_newnuc bypass -- modeptr_aitken <= 0'
       return
    end if

    l_h2so4_sv = l_h2so4
    l_nh3_sv = l_nh3
    lnumait_sv = lnumait
    lnh4ait_sv = lnh4ait
    lso4ait_sv = lso4ait

    dotend(:) = .false.
    dotend(lnumait) = .true.
    dotend(lso4ait) = .true.
    dotend(l_h2so4) = .true.
    if ((l_nh3 > 0) .and. (l_nh3 <= pcnst) .and. &
        (lnh4ait > 0) .and. (lnh4ait <= pcnst)) then
       dotend(lnh4ait) = .true.
       dotend(l_nh3) = .true.
    end if

    do l = 1, pcnst
       if (.not. dotend(l)) cycle
       tmpname = cnst_name(l)
       unit = 'kg/m2/s'
       do m = 1, ntot_amode
          if (l == numptr_amode(m)) unit = '#/m2/s'
       end do
       fieldname = trim(tmpname) // '_sfnnuc1'
       long_name = trim(tmpname) // &
            ' modal_aero new particle nucleation column tendency'
       call addfld(fieldname, unit, 1, 'A', long_name, phys_decomp)
       if (history_aerosol) call add_default(fieldname, 1, ' ')
       if (masterproc) write(*,'(3(a,2x))') &
            'modal_aero_newnuc_init addfld', fieldname, unit
    end do

  end subroutine modal_aero_newnuc_init

end module modal_aero_newnuc
