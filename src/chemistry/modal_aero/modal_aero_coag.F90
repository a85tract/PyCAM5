! CAM compatibility adapter for the decoupled modal aerosol coagulation
! process.  Initialization and process tables remain owned by one scheme
! module instance and are re-exported through the historical CAM module.
module modal_aero_coag

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : gas_pcnst
  use ppgrid, only : pcols, pver
  use modal_aero_data, only : ntot_amode, maxd_aspectype
  use perf_mod, only : t_startf, t_stopf
  use ap_modal_aero_coag_scheme, only : modal_aero_coag_sub_run

  implicit none
  private

  integer, parameter :: pcnstxx = gas_pcnst

#if ( defined MODAL_AERO_7MODE || defined MODAL_AERO_4MODE )
  integer, parameter :: pair_option_acoag = 3
#elif ( defined MODAL_AERO_3MODE )
  integer, parameter :: pair_option_acoag = 1
#endif
  integer, parameter :: maxpair_acoag = 10
  integer, parameter :: maxspec_acoag = maxd_aspectype

  integer :: npair_acoag = 0
  integer :: modefrm_acoag(maxpair_acoag)
  integer :: modetoo_acoag(maxpair_acoag)
  integer :: modetooeff_acoag(maxpair_acoag)
  integer :: nspecfrm_acoag(maxpair_acoag)
  integer :: lspecfrm_acoag(maxspec_acoag,maxpair_acoag)
  integer :: lspectoo_acoag(maxspec_acoag,maxpair_acoag)

  public :: modal_aero_coag_sub, modal_aero_coag_init
  public :: pair_option_acoag, maxpair_acoag, maxspec_acoag
  public :: npair_acoag, modefrm_acoag, modetoo_acoag
  public :: modetooeff_acoag, nspecfrm_acoag
  public :: lspecfrm_acoag, lspectoo_acoag

contains

  subroutine modal_aero_coag_sub(                               &
       lchnk, ncol, nstep, loffset, deltat_main, t, pmid, pdel, &
       q, dgncur_a, dgncur_awet, wetdens_a)

    use cam_abortutils, only : endrun
    use cam_history, only : outfld, fieldname_len
    use chem_mods, only : adv_mass
    use constituents, only : pcnst, cnst_name
    use modal_aero_data
    use modal_aero_gasaerexch, only : n_so4_monolayers_pcage, &
         soa_equivso4_factor
    use physconst, only : gravit, mwdry, tmelt
    use ref_pres, only : top_lev => clim_modal_aero_top_lev

    integer, intent(in) :: lchnk, ncol, nstep, loffset
    real(r8), intent(in) :: deltat_main
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(inout) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dgncur_a(pcols,pver,ntot_amode)
    real(r8), intent(in) :: dgncur_awet(pcols,pver,ntot_amode)
    real(r8), intent(in) :: wetdens_a(pcols,pver,ntot_amode)

    integer :: i, ierr, l, lmz
    real(r8) :: qsrflx(pcols,pcnst)
    logical :: dotend(pcnst)
    character(len=fieldname_len+3) :: fieldname

    call t_startf('ap_modal_aero_coag_sub_run')
    call modal_aero_coag_sub_run(                               &
         lchnk, ncol, nstep, pcols, pver, pcnstxx, ntot_amode, &
         pcnst, top_lev, maxspec_acoag, ntot_aspectype,        &
         maxpair_acoag, loffset,                               &
         deltat_main, pair_option_acoag, npair_acoag,          &
         modefrm_acoag, modetoo_acoag, nspecfrm_acoag,        &
         lspecfrm_acoag, lspectoo_acoag, modeptr_accum,       &
         modeptr_aitken, modeptr_pcarbon, numptr_amode,       &
         mprognum_amode, nspec_amode, lmassptr_amode,         &
         lspectype_amode, sigmag_amode, alnsg_amode,          &
         lptr_so4_a_amode, lptr_nh4_a_amode,                  &
         lptr_soa_a_amode, specmw_amode, specdens_amode,      &
         specmw_so4_amode, specdens_so4_amode,                &
         specmw_nh4_amode, specdens_nh4_amode,                &
         specmw_soa_amode, specdens_soa_amode,                &
         n_so4_monolayers_pcage, soa_equivso4_factor, tmelt, &
         t, pmid, pdel, q, dgncur_a, dgncur_awet, wetdens_a,  &
         qsrflx, dotend, ierr)

    if (ierr /= 0) then
       write(6,*) '*** modal_aero_coag_sub error'
       write(6,*) '    cannot do _coag_sub error pair_option_acoag =', &
            pair_option_acoag
       call endrun('modal_aero_coag_sub error')
    end if

    do l = loffset+1, pcnst
       lmz = l - loffset
       if (.not. dotend(lmz)) cycle
       qsrflx(:,lmz) = qsrflx(:,lmz)*(adv_mass(lmz)/(gravit*mwdry))
       fieldname = trim(cnst_name(l)) // '_sfcoag1'
       call outfld(fieldname, qsrflx(:,lmz), pcols, lchnk)
    end do

    call t_stopf('ap_modal_aero_coag_sub_run')

  end subroutine modal_aero_coag_sub

  subroutine modal_aero_coag_init

    use cam_abortutils, only : endrun
    use cam_history, only : addfld, add_default, fieldname_len, phys_decomp
    use constituents, only : pcnst, cnst_name
    use modal_aero_data
    use modal_aero_gasaerexch, only : modefrm_pcage,           &
         nspecfrm_pcage, lspecfrm_pcage, lspectoo_pcage
    use phys_control, only : phys_getopts
    use spmd_utils, only : masterproc

    integer :: ipair, iq, iqfrm, iqfrm_aa, iqtoo, iqtoo_aa
    integer :: l, lsfrm, lstoo, lunout
    integer :: m, mfrm, mtoo, mtef
    integer :: nsamefrm, nsametoo, nspec
    character(len=fieldname_len) :: tmpname
    character(len=fieldname_len+3) :: fieldname
    character(len=128) :: long_name
    character(len=8) :: unit
    logical :: dotend(pcnst)
    logical :: history_aerosol

    call phys_getopts(history_aerosol_out=history_aerosol)
    lunout = 6

    if (pair_option_acoag == 1) then
       npair_acoag = 1
       modefrm_acoag(1) = modeptr_aitken
       modetoo_acoag(1) = modeptr_accum
       modetooeff_acoag(1) = modeptr_accum
    else if (pair_option_acoag == 2) then
       npair_acoag = 2
       modefrm_acoag(1) = modeptr_aitken
       modetoo_acoag(1) = modeptr_accum
       modetooeff_acoag(1) = modeptr_accum
       modefrm_acoag(2) = modeptr_pcarbon
       modetoo_acoag(2) = modeptr_accum
       modetooeff_acoag(2) = modeptr_accum
    else if (pair_option_acoag == 3) then
       npair_acoag = 3
       modefrm_acoag(1) = modeptr_aitken
       modetoo_acoag(1) = modeptr_accum
       modetooeff_acoag(1) = modeptr_accum
       modefrm_acoag(2) = modeptr_pcarbon
       modetoo_acoag(2) = modeptr_accum
       modetooeff_acoag(2) = modeptr_accum
       modefrm_acoag(3) = modeptr_aitken
       modetoo_acoag(3) = modeptr_pcarbon
       modetooeff_acoag(3) = modeptr_accum
       if (modefrm_pcage <= 0) then
          write(*,*) '*** modal_aero_coag_init error'
          write(*,*) '    pair_option_acoag, modefrm_pcage mismatch'
          write(*,*) '    pair_option_acoag, modefrm_pcage =', &
               pair_option_acoag, modefrm_pcage
          call endrun('modal_aero_coag_init error')
       end if
    else
       npair_acoag = 0
       return
    end if

    do ipair = 1, npair_acoag
       mfrm = modefrm_acoag(ipair)
       mtoo = modetoo_acoag(ipair)
       mtef = modetooeff_acoag(ipair)
       if ((mfrm < 1) .or. (mfrm > ntot_amode) .or. &
           (mtoo < 1) .or. (mtoo > ntot_amode) .or. &
           (mtef < 1) .or. (mtef > ntot_amode)) then
          write(*,*) '*** modal_aero_coag_init error'
          write(*,*) '    ipair, ntot_amode =', ipair, ntot_amode
          write(*,*) '    mfrm, mtoo, mtef  =', mfrm, mtoo, mtef
          call endrun('modal_aero_coag_init error')
       end if

       mtoo = mtef
       nspec = 0
       do iqfrm = 1, nspec_amode(mfrm)
          lsfrm = lmassptr_amode(iqfrm,mfrm)
          if ((lsfrm < 1) .or. (lsfrm > pcnst)) cycle

          iqfrm_aa = 1
          iqtoo_aa = 1
          if (iqfrm > nspec_amode(mfrm)) then
             iqfrm_aa = nspec_amode(mfrm) + 1
             iqtoo_aa = nspec_amode(mtoo) + 1
          end if
          nsamefrm = 0
          do iq = iqfrm_aa, iqfrm
             if (lspectype_amode(iq,mfrm) == &
                 lspectype_amode(iqfrm,mfrm)) nsamefrm = nsamefrm + 1
          end do
          nsametoo = 0
          lstoo = 0
          do iqtoo = iqtoo_aa, nspec_amode(mtoo)
             if (lspectype_amode(iqtoo,mtoo) == &
                 lspectype_amode(iqfrm,mfrm)) then
                nsametoo = nsametoo + 1
                if (nsametoo == nsamefrm) then
                   lstoo = lmassptr_amode(iqtoo,mtoo)
                   exit
                end if
             end if
          end do

          nspec = nspec + 1
          lspecfrm_acoag(nspec,ipair) = lsfrm
          lspectoo_acoag(nspec,ipair) = lstoo
       end do
       nspecfrm_acoag(ipair) = nspec
    end do

    if (masterproc) then
       write(lunout,9310)
       do ipair = 1, npair_acoag
          mfrm = modefrm_acoag(ipair)
          mtoo = modetoo_acoag(ipair)
          mtef = modetooeff_acoag(ipair)
          write(lunout,9320) ipair, mfrm, mtoo, mtef
          do iq = 1, nspecfrm_acoag(ipair)
             lsfrm = lspecfrm_acoag(iq,ipair)
             lstoo = lspectoo_acoag(iq,ipair)
             if (lstoo > 0) then
                write(lunout,9330) lsfrm, cnst_name(lsfrm), &
                     lstoo, cnst_name(lstoo)
             else
                write(lunout,9340) lsfrm, cnst_name(lsfrm)
             end if
          end do
       end do
       write(lunout,*)
    end if

9310 format(/ 'subr. modal_aero_coag_init')
9320 format('pair', i3, 5x, 'mode', i3, ' ---> mode', i3, '   eff', i3)
9330 format(5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a)
9340 format(5x, 'spec', i3, '=', a, ' ---> LOSS')

    dotend(:) = .false.
    do ipair = 1, npair_acoag
       do iq = 1, nspecfrm_acoag(ipair)
          l = lspecfrm_acoag(iq,ipair)
          if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
          l = lspectoo_acoag(iq,ipair)
          if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
       end do
       m = modefrm_acoag(ipair)
       if ((m > 0) .and. (m <= ntot_amode)) then
          l = numptr_amode(m)
          if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
       end if
       m = modetoo_acoag(ipair)
       if ((m > 0) .and. (m <= ntot_amode)) then
          l = numptr_amode(m)
          if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
       end if
    end do

    if (pair_option_acoag == 3) then
       do iq = 1, nspecfrm_pcage
          lsfrm = lspecfrm_pcage(iq)
          lstoo = lspectoo_pcage(iq)
          if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
             dotend(lsfrm) = .true.
             if ((lstoo > 0) .and. (lstoo <= pcnst)) dotend(lstoo) = .true.
          end if
       end do
    end if

    do l = 1, pcnst
       if (.not. dotend(l)) cycle
       tmpname = cnst_name(l)
       unit = 'kg/m2/s'
       do m = 1, ntot_amode
          if (l == numptr_amode(m)) unit = '#/m2/s'
       end do
       fieldname = trim(tmpname) // '_sfcoag1'
       long_name = trim(tmpname) // ' modal_aero coagulation column tendency'
       call addfld(fieldname, unit, 1, 'A', long_name, phys_decomp)
       if (history_aerosol) call add_default(fieldname, 1, ' ')
       if (masterproc) write(*,'(3(a,2x))') &
            'modal_aero_coag_init addfld', fieldname, unit
    end do

  end subroutine modal_aero_coag_init

end module modal_aero_coag
