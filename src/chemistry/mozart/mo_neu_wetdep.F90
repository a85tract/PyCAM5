! CAM compatibility adapter for the standalone Neu wet-deposition process.
module mo_neu_wetdep

  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst
  use ppgrid, only: pcols, pver, pverp
  use perf_mod, only: t_startf, t_stopf
  use gas_wetdep_opts, only: gas_wetdep_method, gas_wetdep_list, gas_wetdep_cnt
  use neu_wetdep_scheme, only: neu_wetdep_configure, neu_wetdep_tend_run

  implicit none
  private
  save

  logical, public :: do_neu_wetdep = .false.
  logical, parameter :: debug = .false.
  logical, parameter :: do_diag = .false.
  integer, allocatable :: mapping_to_heff(:), mapping_to_mmr(:)
  real(r8), allocatable :: mol_weight(:)
  logical, allocatable :: ice_uptake(:)

  public :: neu_wetdep_init
  public :: neu_wetdep_tend

contains

  subroutine neu_wetdep_init()
    use cam_abortutils, only: endrun
    use cam_history, only: addfld, add_default, phys_decomp
    use constituents, only: cnst_get_ind, cnst_mw
    use phys_control, only: phys_getopts
    use seq_drydep_mod, only: n_species_table, species_name_table, dheff
    use spmd_utils, only: masterproc

    integer :: m, l
    integer :: index_cldice, index_cldliq, nh3_ndx, co2_ndx, hno3_ndx
    character(len=20) :: test_name
    logical :: history_chemistry

    nh3_ndx = 0
    co2_ndx = 0
    hno3_ndx = 0
    call phys_getopts(history_chemistry_out=history_chemistry)
    do_neu_wetdep = gas_wetdep_method == 'NEU' .and. gas_wetdep_cnt > 0
    if (.not. do_neu_wetdep) return

    allocate(mapping_to_heff(gas_wetdep_cnt), mapping_to_mmr(gas_wetdep_cnt))
    allocate(ice_uptake(gas_wetdep_cnt), mol_weight(gas_wetdep_cnt))

    mapping_to_heff = -99
    do m = 1, gas_wetdep_cnt
       test_name = gas_wetdep_list(m)
       select case (trim(test_name))
       case ('HYAC', 'CH3COOH', 'HCOOH', 'EOOH')
          test_name = 'CH2O'
       case ('SOGB', 'SOGI', 'SOGM', 'SOGT', 'SOGX')
          test_name = 'H2O2'
       case ('SO2t')
          test_name = 'SO2'
       case ('CLONO2', 'BRONO2', 'HCL', 'HOCL', 'HOBR', 'HBR', 'Pb', &
             'MACROOH', 'ISOPOOH', 'XOOH', 'H2SO4', 'HF', 'COF2', 'COFCL')
          test_name = 'HNO3'
       case ('NH_50W')
          test_name = 'HNO3'
       case ('ALKOOH', 'MEKOOH', 'TOLOOH', 'TERPOOH')
          test_name = 'CH3OOH'
       end select

       do l = 1, n_species_table
          if (trim(test_name) == trim(species_name_table(l))) then
             mapping_to_heff(m) = l
             exit
          end if
       end do
       if (mapping_to_heff(m) == -99 .and. masterproc) &
            print *, 'problem with mapping_to_heff of ', trim(test_name)
       if (trim(test_name) == 'NH3') nh3_ndx = m
       if (trim(test_name) == 'CO2') co2_ndx = m
       if (trim(test_name) == 'HNO3') hno3_ndx = m
    end do
    if (any(mapping_to_heff == -99)) &
         call endrun('mo_neu_wetdep: unmapped Henry-law species')

    mapping_to_mmr = -99
    do m = 1, gas_wetdep_cnt
       call cnst_get_ind(gas_wetdep_list(m), mapping_to_mmr(m), abort=.false.)
       if (mapping_to_mmr(m) <= 0) &
            call endrun('problem with mapping_to_mmr of '//trim(gas_wetdep_list(m)))
       mol_weight(m) = cnst_mw(mapping_to_mmr(m))
       ice_uptake(m) = trim(gas_wetdep_list(m)) == 'HNO3'
    end do

    call cnst_get_ind('CLDICE', index_cldice)
    call cnst_get_ind('CLDLIQ', index_cldliq)

    do m = 1, gas_wetdep_cnt
       call addfld('DTWR_'//trim(gas_wetdep_list(m)), 'mol/mol/s', pver, 'A', &
            'wet removal Neu scheme tendency', phys_decomp)
       if (history_chemistry) call add_default('DTWR_'//trim(gas_wetdep_list(m)), 1, ' ')
    end do
    if (do_diag) then
       call addfld('QT_RAIN_HNO3','kg',pver,'A','tracer in rain',phys_decomp)
       call addfld('QT_RIME_HNO3','kg',pver,'A','tracer in rimed cloud water',phys_decomp)
       call addfld('QT_WASH_HNO3','kg',pver,'A','tracer in washout',phys_decomp)
       call addfld('QT_EVAP_HNO3','kg',pver,'A','tracer release by evaporation',phys_decomp)
       call add_default('QT_RAIN_HNO3',1,' ')
       call add_default('QT_RIME_HNO3',1,' ')
       call add_default('QT_WASH_HNO3',1,' ')
       call add_default('QT_EVAP_HNO3',1,' ')
    end if

    call neu_wetdep_configure(do_neu_wetdep, mapping_to_heff, mapping_to_mmr, &
         mol_weight, ice_uptake, index_cldice, index_cldliq, nh3_ndx, &
         co2_ndx, hno3_ndx, dheff)
  end subroutine neu_wetdep_init

  subroutine neu_wetdep_tend(lchnk, ncol, mmr, pmid, pdel, zint, tfld, delt, &
       prain, nevapr, cld, cmfdqr, wd_tend)
    use cam_history, only: outfld
    use phys_grid, only: get_area_all_p, get_rlat_all_p

    integer, intent(in) :: lchnk, ncol
    real(r8), intent(in) :: mmr(pcols,pver,pcnst)
    real(r8), intent(in) :: pmid(pcols,pver), pdel(pcols,pver)
    real(r8), intent(in) :: zint(pcols,pver+1), tfld(pcols,pver)
    real(r8), intent(in) :: delt
    real(r8), intent(in) :: prain(ncol,pver), nevapr(ncol,pver)
    real(r8), intent(in) :: cld(ncol,pver), cmfdqr(ncol,pver)
    real(r8), intent(inout) :: wd_tend(pcols,pver,pcnst)

    real(r8) :: area(pcols), lats(pcols)
    real(r8) :: dtwr(pcols,pver,gas_wetdep_cnt)
    real(r8) :: qt_rain(ncol,pver), qt_rime(ncol,pver)
    real(r8) :: qt_wash(ncol,pver), qt_evap(ncol,pver)
    integer :: m

    call get_area_all_p(lchnk, ncol, area)
    call get_rlat_all_p(lchnk, pcols, lats)

    call t_startf('ap_neu_wetdep_tend_run')
    call neu_wetdep_tend_run(lchnk, ncol, pcols, pver, pverp, pcnst, &
         area, lats, mmr, pmid, pdel, zint, tfld, delt, prain, nevapr, &
         cld, cmfdqr, wd_tend, dtwr, qt_rain, qt_rime, qt_wash, qt_evap)
    call t_stopf('ap_neu_wetdep_tend_run')

    do m = 1, gas_wetdep_cnt
       call outfld('DTWR_'//trim(gas_wetdep_list(m)), dtwr(:,:,m), ncol, lchnk)
    end do
    if (do_diag) then
       call outfld('QT_RAIN_HNO3',qt_rain,ncol,lchnk)
       call outfld('QT_RIME_HNO3',qt_rime,ncol,lchnk)
       call outfld('QT_WASH_HNO3',qt_wash,ncol,lchnk)
       call outfld('QT_EVAP_HNO3',qt_evap,ncol,lchnk)
    end if
  end subroutine neu_wetdep_tend

end module mo_neu_wetdep
