! CAM facade for the standalone modal aerosol gas/aerosol exchange kernel.
module modal_aero_gasaerexch

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : gas_pcnst
  use modal_aero_data, only : maxd_aspectype
  use ppgrid, only : pcols, pver
  use perf_mod, only : t_startf, t_stopf
  use ap_modal_aero_gasaerexch_scheme, only : &
       modal_aero_gasaerexch_sub_run

  implicit none
  private
  save

  integer, parameter :: pcnstxx = gas_pcnst
  integer, parameter, public :: maxspec_pcage = maxd_aspectype

  integer, public :: npair_pcage
  integer, public :: modefrm_pcage
  integer, public :: modetoo_pcage
  integer, public :: nspecfrm_pcage
  integer, public :: lspecfrm_pcage(maxspec_pcage)
  integer, public :: lspectoo_pcage(maxspec_pcage)

  real(r8), parameter, public :: n_so4_monolayers_pcage = 8.0_r8
  real(r8), parameter, public :: dr_so4_monolayers_pcage = &
       n_so4_monolayers_pcage * 4.76e-10_r8
  real(r8), public :: soa_equivso4_factor

  public :: modal_aero_gasaerexch_sub, modal_aero_gasaerexch_init

contains

  subroutine modal_aero_gasaerexch_sub(                         &
       lchnk, ncol, nstep, loffset, deltat, t, pmid, pdel, qh2o, &
       troplev, q, qqcw, dqdt_other, dqqcwdt_other, dgncur_a,    &
       dgncur_awet, sulfeq)

    use cam_abortutils, only : endrun
    use cam_history, only : outfld, fieldname_len
    use chem_mods, only : adv_mass
    use constituents, only : pcnst, cnst_get_ind, cnst_name
    use modal_aero_data
    use physconst, only : gravit, mwdry, rair
    use ref_pres, only : top_lev => clim_modal_aero_top_lev

    integer, intent(in) :: lchnk, ncol, nstep, loffset
    integer, intent(in) :: troplev(pcols)
    real(r8), intent(in) :: deltat
    real(r8), intent(inout) :: q(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: qqcw(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqdt_other(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqqcwdt_other(ncol,pver,pcnstxx)
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: qh2o(pcols,pver)
    real(r8), intent(in) :: dgncur_a(pcols,pver,ntot_amode)
    real(r8), intent(in) :: dgncur_awet(pcols,pver,ntot_amode)
    real(r8), pointer :: sulfeq(:,:,:)

    integer :: ierr, rename_error
    integer :: l, lb, l_so4g, l_nh4g, l_msag, l_soag
    real(r8) :: qsrflx(pcols,pcnstxx,2)
    real(r8) :: qqcwsrflx(pcols,pcnstxx,2)
    real(r8), target :: sulfeq_zero(pcols,pver,ntot_amode)
    real(r8), pointer :: sulfeq_arg(:,:,:)
    logical :: has_sulfeq
    logical :: dotend(pcnstxx), dotendqqcw(pcnstxx)
    logical :: dotendrn(pcnstxx), dotendqqcwrn(pcnstxx)
    character(len=fieldname_len+3) :: fieldname

    call cnst_get_ind('H2SO4', l_so4g, .false.)
    call cnst_get_ind('NH3',   l_nh4g, .false.)
    call cnst_get_ind('MSA',   l_msag, .false.)
    call cnst_get_ind('SOAG',  l_soag, .false.)
    l_so4g = l_so4g - loffset
    l_nh4g = l_nh4g - loffset
    l_msag = l_msag - loffset
    l_soag = l_soag - loffset

    has_sulfeq = associated(sulfeq)
    if (has_sulfeq) then
       sulfeq_arg => sulfeq
    else
       sulfeq_zero = 0.0_r8
       sulfeq_arg => sulfeq_zero
    end if

    call t_startf('ap_modal_aero_gasaerexch_sub_run')
    call modal_aero_gasaerexch_sub_run(                         &
         lchnk, ncol, nstep, pcols, pver, pcnstxx, ntot_amode, &
         pcnst, top_lev, maxspec_pcage, ntot_aspectype,        &
         loffset, deltat, l_so4g, l_nh4g, l_msag, l_soag,     &
         modefrm_pcage, modetoo_pcage, nspecfrm_pcage,        &
         lspecfrm_pcage, lspectoo_pcage, modeptr_pcarbon,     &
         numptr_amode, nspec_amode, lmassptr_amode,           &
         lspectype_amode, lptr_so4_a_amode,                   &
         lptr_nh4_a_amode, lptr_soa_a_amode,                  &
         lptr_pom_a_amode, sigmag_amode, alnsg_amode,        &
         specmw_amode, specdens_amode, specmw_so4_amode,     &
         specmw_nh4_amode, specmw_soa_amode,                 &
         specdens_so4_amode, specdens_nh4_amode,             &
         specdens_soa_amode, dr_so4_monolayers_pcage,        &
         soa_equivso4_factor, gravit, mwdry, rair, adv_mass, &
         t, pmid, pdel, qh2o, troplev, q, qqcw,              &
         dqdt_other, dqqcwdt_other, dgncur_a, dgncur_awet,   &
         has_sulfeq, sulfeq_arg, qsrflx, qqcwsrflx,          &
         dotend, dotendqqcw, dotendrn, dotendqqcwrn,         &
         ierr, rename_error)

    if (ierr /= 0) then
       write(*,'(/a/a,2i7)') &
            '*** modal_aero_gasaerexch_sub -- cannot find H2SO4 species', &
            '    l_so4g, loffset =', l_so4g, loffset
       call endrun('modal_aero_gasaerexch_sub error')
    end if
    if (rename_error /= 0) then
       write(*,'(/a,1x,i7)') &
            '*** modal_aero_gasaerexch_sub -- renaming error pair =', &
            rename_error
       call endrun('modal_aero_gasaerexch_sub error')
    end if

    do l = 1, pcnstxx
       lb = l + loffset
       if (dotend(l)) then
          fieldname = trim(cnst_name(lb)) // '_sfgaex1'
          call outfld(fieldname, qsrflx(:,l,1), pcols, lchnk)
       end if
       if (dotendrn(l)) then
          fieldname = trim(cnst_name(lb)) // '_sfgaex2'
          call outfld(fieldname, qsrflx(:,l,2), pcols, lchnk)
       end if
       if (dotendqqcwrn(l)) then
          fieldname = trim(cnst_name_cw(lb)) // '_sfgaex2'
          call outfld(fieldname, qqcwsrflx(:,l,2), pcols, lchnk)
       end if
    end do
    call t_stopf('ap_modal_aero_gasaerexch_sub_run')

  end subroutine modal_aero_gasaerexch_sub

  subroutine modal_aero_gasaerexch_init

    use cam_abortutils, only : endrun
    use cam_history, only : addfld, add_default, fieldname_len, phys_decomp
    use constituents, only : pcnst, cnst_get_ind, cnst_name
    use modal_aero_data
    use modal_aero_rename
    use phys_control, only : phys_getopts
    use spmd_utils, only : masterproc

    integer :: ipair, iq, iqfrm, iqtoo
    integer :: jac, l, lsfrm, lstoo, lunout
    integer :: l_so4g, l_nh4g, l_msag, l_soag
    integer :: mfrm, mtoo, nspec, n
    logical :: do_msag, do_nh4g, do_soag
    logical :: dotend(pcnst), dotendqqcw(pcnst)
    real(r8) :: tmp1, tmp2
    character(len=fieldname_len) :: tmpnamea
    character(len=fieldname_len+3) :: fieldname
    character(len=128) :: long_name
    character(len=8) :: unit
    logical :: history_aerosol

    call phys_getopts(history_aerosol_out=history_aerosol)
    lunout = 6

!   Define the primary-carbon aging pair and species mapping.
    modefrm_pcage = -999888777
    modetoo_pcage = -999888777
    if ((modeptr_pcarbon <= 0) .or. (modeptr_accum <= 0)) goto 15000
    l = lptr_so4_a_amode(modeptr_accum)
    if ((l < 1) .or. (l > pcnst)) goto 15000

    modefrm_pcage = modeptr_pcarbon
    modetoo_pcage = modeptr_accum
    mfrm = modefrm_pcage
    mtoo = modetoo_pcage
    nspec = 0
aa_iqfrm: do iqfrm = -1, nspec_amode(mfrm)
       if (iqfrm == -1) then
          lsfrm = numptr_amode(mfrm)
          lstoo = numptr_amode(mtoo)
       else if (iqfrm == 0) then
          cycle aa_iqfrm
       else
          lsfrm = lmassptr_amode(iqfrm,mfrm)
          lstoo = 0
       end if
       if ((lsfrm < 1) .or. (lsfrm > pcnst)) cycle aa_iqfrm

       if ((lsfrm > 0) .and. (iqfrm > 0)) then
          do iqtoo = 1, nspec_amode(mtoo)
             if (lspectype_amode(iqtoo,mtoo) == &
                  lspectype_amode(iqfrm,mfrm)) then
                lstoo = lmassptr_amode(iqtoo,mtoo)
                exit
             end if
          end do
       end if
       if ((lstoo < 1) .or. (lstoo > pcnst)) lstoo = 0
       nspec = nspec + 1
       lspecfrm_pcage(nspec) = lsfrm
       lspectoo_pcage(nspec) = lstoo
    end do aa_iqfrm
    nspecfrm_pcage = nspec

    if (masterproc) then
       write(lunout,9310)
       write(lunout,9320) 1, modefrm_pcage, modetoo_pcage
       do iq = 1, nspecfrm_pcage
          lsfrm = lspecfrm_pcage(iq)
          lstoo = lspectoo_pcage(iq)
          if (lstoo > 0) then
             write(lunout,9330) lsfrm, cnst_name(lsfrm), &
                  lstoo, cnst_name(lstoo)
          else
             write(lunout,9340) lsfrm, cnst_name(lsfrm)
          end if
       end do
       write(lunout,*)
    end if

9310 format(/ 'subr. modal_aero_gasaerexch_init - primary carbon aging pointers')
9320 format('pair', i3, 5x, 'mode', i3, ' ---> mode', i3)
9330 format(5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a)
9340 format(5x, 'spec', i3, '=', a, ' ---> LOSS')

15000 continue

    call cnst_get_ind('H2SO4', l_so4g, .false.)
    call cnst_get_ind('NH3',   l_nh4g, .false.)
    call cnst_get_ind('MSA',   l_msag, .false.)
    call cnst_get_ind('SOAG',  l_soag, .false.)
    if ((l_so4g <= 0) .or. (l_so4g > pcnst)) then
       write(*,'(/a/a,2i7)') &
            '*** modal_aero_gasaerexch_init -- cannot find H2SO4 species', &
            '    l_so4g=', l_so4g
       call endrun('modal_aero_gasaerexch_init error')
    end if
    do_nh4g = (l_nh4g > 0) .and. (l_nh4g <= pcnst)
    do_msag = (l_msag > 0) .and. (l_msag <= pcnst)
    do_soag = (l_soag > 0) .and. (l_soag <= pcnst)

    dotend(:) = .false.
    dotend(l_so4g) = .true.
    if (do_nh4g) dotend(l_nh4g) = .true.
    if (do_msag) dotend(l_msag) = .true.
    if (do_soag) dotend(l_soag) = .true.
    do n = 1, ntot_amode
       l = lptr_so4_a_amode(n)
       if ((l > 0) .and. (l <= pcnst)) then
          dotend(l) = .true.
          if (do_nh4g) then
             l = lptr_nh4_a_amode(n)
             if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
          end if
       end if
       l = lptr_soa_a_amode(n)
       if ((l > 0) .and. (l <= pcnst)) dotend(l) = .true.
    end do

    if (modefrm_pcage > 0) then
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
       tmpnamea = cnst_name(l)
       fieldname = trim(tmpnamea) // '_sfgaex1'
       long_name = trim(tmpnamea) // &
            ' gas-aerosol-exchange primary column tendency'
       unit = 'kg/m2/s'
       call addfld(fieldname, unit, 1, 'A', long_name, phys_decomp)
       if (history_aerosol) call add_default(fieldname, 1, ' ')
       if (masterproc) write(*,'(3(a,3x))') &
            'gasaerexch addfld', fieldname, unit
    end do

    dotend(:) = .false.
    dotendqqcw(:) = .false.
    do ipair = 1, npair_renamexf
       do iq = 1, nspecfrm_renamexf(ipair)
          lsfrm = lspecfrma_renamexf(iq,ipair)
          lstoo = lspectooa_renamexf(iq,ipair)
          if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
             dotend(lsfrm) = .true.
             if ((lstoo > 0) .and. (lstoo <= pcnst)) dotend(lstoo) = .true.
          end if
          lsfrm = lspecfrmc_renamexf(iq,ipair)
          lstoo = lspectooc_renamexf(iq,ipair)
          if ((lsfrm > 0) .and. (lsfrm <= pcnst)) then
             dotendqqcw(lsfrm) = .true.
             if ((lstoo > 0) .and. (lstoo <= pcnst)) &
                  dotendqqcw(lstoo) = .true.
          end if
       end do
    end do

    do l = 1, pcnst
       do jac = 1, 2
          if (jac == 1) then
             if (.not. dotend(l)) cycle
             tmpnamea = cnst_name(l)
          else
             if (.not. dotendqqcw(l)) cycle
             tmpnamea = cnst_name_cw(l)
          end if
          fieldname = trim(tmpnamea) // '_sfgaex2'
          long_name = trim(tmpnamea) // &
               ' gas-aerosol-exchange renaming column tendency'
          unit = 'kg/m2/s'
          if ((tmpnamea(1:3) == 'num') .or. &
              (tmpnamea(1:3) == 'NUM')) unit = '#/m2/s'
          call addfld(fieldname, unit, 1, 'A', long_name, phys_decomp)
          if (history_aerosol) call add_default(fieldname, 1, ' ')
          if (masterproc) write(*,'(3(a,3x))') &
               'gasaerexch addfld', fieldname, unit
       end do
    end do

    soa_equivso4_factor = 0.0_r8
    if (do_soag) then
       tmp1 = -1.0_r8
       tmp2 = -1.0_r8
       do l = 1, ntot_aspectype
          if (specname_amode(l) == 's-organic') tmp1 = spechygro(l)
          if (specname_amode(l) == 'sulfate') tmp2 = spechygro(l)
       end do
       if ((tmp1 > 0.0_r8) .and. (tmp2 > 0.0_r8)) then
          soa_equivso4_factor = tmp1/tmp2
       else
          write(*,'(a/a,1p,2e10.2)') &
               '*** subr modal_aero_gasaerexch_init', &
               '    cannot find hygros - tmp1/2 =', tmp1, tmp2
          call endrun()
       end if
    end if

  end subroutine modal_aero_gasaerexch_init

end module modal_aero_gasaerexch
