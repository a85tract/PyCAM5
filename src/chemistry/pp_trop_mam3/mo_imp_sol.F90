module mo_imp_sol
  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods, only : clscnt4, gas_pcnst, clsmap
  use cam_logfile, only : iulog
  implicit none
  private
  public :: imp_slv_inti, imp_sol
  public :: itermax, epsilon, factor, small, imp_sol_scheme_entry_count
  save
  real(r8), parameter :: rel_err = 1.e-3_r8
  real(r8), parameter :: high_rel_err = 1.e-4_r8
  !-----------------------------------------------------------------------
  ! Newton-Raphson iteration limits
  !-----------------------------------------------------------------------
  integer, parameter :: itermax = 11
  integer, parameter :: cut_limit = 5
  real(r8) :: small
  real(r8) :: epsilon(clscnt4)
  logical :: factor(itermax)
  integer :: ox_ndx
  integer :: o1d_ndx = -1
  integer :: h2o_ndx = -1
  integer :: oh_ndx, ho2_ndx, ch3o2_ndx, po2_ndx, ch3co3_ndx
  integer :: c2h5o2_ndx, isopo2_ndx, macro2_ndx, mco3_ndx, c3h7o2_ndx
  integer :: ro2_ndx, xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx
  integer :: c2h4_ndx, c3h6_ndx, isop_ndx, mvk_ndx, c10h16_ndx
  integer :: ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, ox_p4_ndx, ox_p5_ndx
  integer :: ox_p6_ndx, ox_p7_ndx, ox_p8_ndx, ox_p9_ndx, ox_p10_ndx
  integer :: ox_p11_ndx
  integer :: ox_l1_ndx, ox_l2_ndx, ox_l3_ndx, ox_l4_ndx, ox_l5_ndx
  integer :: ox_l6_ndx, ox_l7_ndx, ox_l8_ndx, ox_l9_ndx, usr4_ndx
  integer :: usr16_ndx, usr17_ndx, c2o3_ndx,ole_ndx
  integer :: tolo2_ndx, terpo2_ndx, alko2_ndx, eneo2_ndx, eo2_ndx, meko2_ndx
  integer :: ox_p17_ndx,ox_p12_ndx,ox_p13_ndx,ox_p14_ndx,ox_p15_ndx,ox_p16_ndx
  integer :: lt_cnt
  logical :: full_ozone_chem = .false.
  logical :: middle_atm_chem = .false.
  logical :: reduced_ozone_chem = .false.
  ! for xnox ozone chemistry diagnostics
  integer :: o3a_ndx, xno2_ndx, no2xno3_ndx, xno2no3_ndx, xno3_ndx, o1da_ndx, xno_ndx
  integer :: usr4a_ndx, usr16a_ndx, usr16b_ndx, usr17b_ndx
  integer :: imp_sol_scheme_entry_count = 0
!$omp threadprivate(imp_sol_scheme_entry_count)
contains
  subroutine imp_slv_inti
    !-----------------------------------------------------------------------
    ! ... Initialize the implict solver
    !-----------------------------------------------------------------------
    use mo_chem_utls,   only : get_spc_ndx, get_rxt_ndx
    use cam_abortutils, only : endrun
    use cam_history,    only : addfld, add_default, phys_decomp
    use ppgrid,         only : pver
    use mo_tracname,    only : solsym
    implicit none
    !-----------------------------------------------------------------------
    ! ... Local variables
    !-----------------------------------------------------------------------
    integer :: m
    real(r8) :: eps(gas_pcnst)
    integer :: wrk(27)
    integer :: i,j
    ! small = 1.e6_r8 * tiny( small )
    small = 1.e-40_r8
    factor(:) = .true.
    eps(:) = rel_err
    ox_ndx = get_spc_ndx( 'OX' )
    h2o_ndx = get_spc_ndx( 'H2O' )
    if( ox_ndx < 1 ) then
       ox_ndx = get_spc_ndx( 'O3' )
       o1d_ndx = get_spc_ndx( 'O1D' )
       o1da_ndx = get_spc_ndx( 'O1DA' )
    end if
    if( ox_ndx > 0 ) then
       eps(ox_ndx) = high_rel_err
    end if
    m = get_spc_ndx( 'NO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'N2O5' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'OH' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    o3a_ndx = get_spc_ndx( 'O3A' )
    if( o3a_ndx > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    !do m = 1,max(1,clscnt4)
    do m = 1,clscnt4
       epsilon(m) = eps(clsmap(m,4))
    end do
    has_o3_chem: if( ox_ndx > 0 ) then
       ox_p1_ndx = get_rxt_ndx( 'ox_p1' )
       if ( ox_p1_ndx < 0 ) then
          ox_p1_ndx = get_rxt_ndx( 'NO_HO2' )
       end if
       if ( ox_p1_ndx < 0 ) then
          ox_p1_ndx = get_rxt_ndx( 'cph_NO_HO2' )
       end if
       ox_p2_ndx = get_rxt_ndx( 'ox_p2' )
       if ( ox_p2_ndx < 0 ) then
          ox_p2_ndx = get_rxt_ndx( 'CH3O2_NO' )
       end if
       ox_p3_ndx = get_rxt_ndx( 'ox_p3' )
       if ( ox_p3_ndx < 0 ) then
          ox_p3_ndx = get_rxt_ndx( 'PO2_NO' )
       end if
       ox_p4_ndx = get_rxt_ndx( 'ox_p4' )
       if ( ox_p4_ndx < 0 ) then
          ox_p4_ndx = get_rxt_ndx( 'CH3CO3_NO' )
       end if
       ox_p5_ndx = get_rxt_ndx( 'ox_p5' )
       if ( ox_p5_ndx < 0 ) then
          ox_p5_ndx = get_rxt_ndx( 'C2H5O2_NO' )
       end if
       ox_p6_ndx = get_rxt_ndx( 'ox_p6' )
       if ( ox_p6_ndx < 0 ) then
          ox_p6_ndx = get_rxt_ndx( 'ISOPO2_NO' )
       end if
       ox_p7_ndx = get_rxt_ndx( 'ox_p7' )
       if ( ox_p7_ndx < 0 ) then
          ox_p7_ndx = get_rxt_ndx( 'MACRO2_NO' )
       end if
       if ( ox_p7_ndx < 0 ) then
          ox_p7_ndx = get_rxt_ndx( 'MACRO2_NOa' )
       end if
       ox_p8_ndx = get_rxt_ndx( 'ox_p8' )
       if ( ox_p8_ndx < 0 ) then
          ox_p8_ndx = get_rxt_ndx( 'MCO3_NO' )
       end if
       ox_p9_ndx = get_rxt_ndx( 'ox_p9' )
       if ( ox_p9_ndx < 0 ) then
          ox_p9_ndx = get_rxt_ndx( 'C3H7O2_NO' )
       end if
       ox_p10_ndx = get_rxt_ndx( 'ox_p10' )
       if ( ox_p10_ndx < 0 ) then
          ox_p10_ndx = get_rxt_ndx( 'RO2_NO' )
       end if
       ox_p11_ndx = get_rxt_ndx( 'ox_p11' )
       if ( ox_p11_ndx < 0 ) then
          ox_p11_ndx = get_rxt_ndx( 'XO2_NO' )
       end if
       if ( ox_p11_ndx < 0 ) then
          ox_p11_ndx = get_rxt_ndx( 'tag_XO2_NO' )
       end if
       if ( ox_p11_ndx < 0 ) then
          ox_p11_ndx = get_rxt_ndx( 'r63' )
       end if
       ox_p12_ndx = get_rxt_ndx( 'ox_p12' )
       if ( ox_p12_ndx < 0 ) then
          ox_p12_ndx = get_rxt_ndx( 'TOLO2_NO' )
       end if
       ox_p13_ndx = get_rxt_ndx( 'ox_p13' )
       if ( ox_p13_ndx < 0 ) then
          ox_p13_ndx = get_rxt_ndx( 'TERPO2_NO' )
       end if
       ox_p14_ndx = get_rxt_ndx( 'ox_p14' )
       if ( ox_p14_ndx < 0 ) then
          ox_p14_ndx = get_rxt_ndx( 'ALKO2_NO' )
       end if
       ox_p15_ndx = get_rxt_ndx( 'ox_p15' )
       if ( ox_p15_ndx < 0 ) then
          ox_p15_ndx = get_rxt_ndx( 'ENEO2_NO' )
       end if
       ox_p16_ndx = get_rxt_ndx( 'ox_p16' )
       if ( ox_p16_ndx < 0 ) then
          ox_p16_ndx = get_rxt_ndx( 'EO2_NO' )
       end if
       ox_p17_ndx = get_rxt_ndx( 'ox_p17' )
       if ( ox_p17_ndx < 0 ) then
          ox_p17_ndx = get_rxt_ndx( 'MEKO2_NO' )
       end if
       wrk(1:17) = (/ ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, ox_p4_ndx, ox_p5_ndx, &
            ox_p6_ndx, ox_p7_ndx, ox_p8_ndx, ox_p9_ndx, ox_p10_ndx, ox_p11_ndx, &
            ox_p12_ndx, ox_p13_ndx, ox_p14_ndx, ox_p15_ndx, ox_p16_ndx, ox_p17_ndx /)
       if( all( wrk(1:17) > 0 ) ) then
          full_ozone_chem = .true.
       end if
       if ( ox_p11_ndx < 0 ) then
          ox_p11_ndx = get_rxt_ndx( 'tag_XO2_NO' )
       end if
       ox_p12_ndx = get_rxt_ndx( 'ox_p12' )
       if ( ox_p12_ndx < 0 ) then
          ox_p12_ndx = get_rxt_ndx( 'TOLO2_NO' )
       end if
       ox_p13_ndx = get_rxt_ndx( 'ox_p13' )
       if ( ox_p13_ndx < 0 ) then
          ox_p13_ndx = get_rxt_ndx( 'TERPO2_NO' )
       end if
       ox_p14_ndx = get_rxt_ndx( 'ox_p14' )
       if ( ox_p14_ndx < 0 ) then
          ox_p14_ndx = get_rxt_ndx( 'ALKO2_NO' )
       end if
       ox_p15_ndx = get_rxt_ndx( 'ox_p15' )
       if ( ox_p15_ndx < 0 ) then
          ox_p15_ndx = get_rxt_ndx( 'ENEO2_NO' )
       end if
       ox_p16_ndx = get_rxt_ndx( 'ox_p16' )
       if ( ox_p16_ndx < 0 ) then
          ox_p16_ndx = get_rxt_ndx( 'EO2_NO' )
       end if
       ox_p17_ndx = get_rxt_ndx( 'ox_p17' )
       if ( ox_p17_ndx < 0 ) then
          ox_p17_ndx = get_rxt_ndx( 'MEKO2_NO' )
       end if
       wrk(1:17) = (/ ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, ox_p4_ndx, ox_p5_ndx, &
            ox_p6_ndx, ox_p7_ndx, ox_p8_ndx, ox_p9_ndx, ox_p10_ndx, ox_p11_ndx, &
            ox_p12_ndx, ox_p13_ndx, ox_p14_ndx, ox_p15_ndx, ox_p16_ndx, ox_p17_ndx /)
       if( all( wrk(1:17) > 0 ) ) then
          full_ozone_chem = .true.
       end if
       if ( .not. full_ozone_chem ) then
          wrk(1:4) = (/ ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, ox_p11_ndx/)
          if( all( wrk(1:4) > 0 ) ) then
             reduced_ozone_chem = .true.
          end if
          if ( .not. reduced_ozone_chem ) then
             wrk(1:2) = (/ ox_p1_ndx, ox_p2_ndx/)
             if( all( wrk(1:2) > 0 ) ) then
                middle_atm_chem = .true.
             end if
          end if
       endif
       if( full_ozone_chem .or. reduced_ozone_chem .or. middle_atm_chem ) then
          ox_l1_ndx = get_rxt_ndx( 'ox_l1' )
          if ( ox_l1_ndx < 0 ) then
            ox_l1_ndx = get_rxt_ndx( 'O1D_H2O' )
          end if
          ox_l2_ndx = get_rxt_ndx( 'ox_l2' )
          if ( ox_l2_ndx < 0 ) then
            ox_l2_ndx = get_rxt_ndx( 'OH_O3' )
          end if
          if ( ox_l2_ndx < 0 ) then
            ox_l2_ndx = get_rxt_ndx( 'cph_OH_O3' )
          end if
          ox_l3_ndx = get_rxt_ndx( 'ox_l3' )
          if ( ox_l3_ndx < 0 ) then
            ox_l3_ndx = get_rxt_ndx( 'HO2_O3' )
          end if
          if ( ox_l3_ndx < 0 ) then
            ox_l3_ndx = get_rxt_ndx( 'cph_HO2_O3' )
          end if
          ox_l4_ndx = get_rxt_ndx( 'ox_l4' )
          if ( ox_l4_ndx < 0 ) then
            ox_l4_ndx = get_rxt_ndx( 'C3H6_O3' )
          end if
          ox_l5_ndx = get_rxt_ndx( 'ox_l5' )
          if ( ox_l5_ndx < 0 ) then
            ox_l5_ndx = get_rxt_ndx( 'ISOP_O3' )
          end if
          ox_l6_ndx = get_rxt_ndx( 'ox_l6' )
          if ( ox_l6_ndx < 0 ) then
            ox_l6_ndx = get_rxt_ndx( 'C2H4_O3' )
          end if
          ox_l7_ndx = get_rxt_ndx( 'ox_l7' )
          if ( ox_l7_ndx < 0 ) then
            ox_l7_ndx = get_rxt_ndx( 'MVK_O3' )
          end if
          ox_l8_ndx = get_rxt_ndx( 'ox_l8' )
          if ( ox_l8_ndx < 0 ) then
            ox_l8_ndx = get_rxt_ndx( 'MACR_O3' )
          end if
          ox_l9_ndx = get_rxt_ndx( 'ox_l9' )
          if( ox_l9_ndx < 1 ) then
             ox_l9_ndx = get_rxt_ndx( 'soa1' )
          end if
          if( ox_l9_ndx < 1 ) then
             ox_l9_ndx = get_rxt_ndx( 'C10H16_O3' )
          end if
          usr4_ndx = get_rxt_ndx( 'usr4' )
          if (usr4_ndx < 1) then
             usr4_ndx = get_rxt_ndx( 'tag_NO2_OH' )
          endif
          usr16_ndx = get_rxt_ndx( 'usr16' )
          if (usr16_ndx < 1) then
            usr16_ndx = get_rxt_ndx( 'usr_N2O5_aer' )
          endif
          usr17_ndx = get_rxt_ndx( 'usr17' )
          if (usr17_ndx < 1) then
            usr17_ndx = get_rxt_ndx( 'usr_NO3_aer' )
          endif
          usr4a_ndx = get_rxt_ndx( 'usr4a' )
          if (usr4a_ndx < 1) then
            usr4a_ndx = get_rxt_ndx( 'tag_XNO2_OH' )
          endif
          usr16b_ndx = get_rxt_ndx( 'usr16b' )
          if (usr16b_ndx < 1) then
            usr16b_ndx = get_rxt_ndx( 'usr_NO2XNO3_aer' )
          endif
          usr16a_ndx = get_rxt_ndx( 'usr16a' )
          if (usr16a_ndx < 1) then
            usr16a_ndx = get_rxt_ndx( 'usr_XNO2NO3_aer' )
          endif
          usr17b_ndx = get_rxt_ndx( 'usr17b' )
          if (usr17b_ndx < 1) then
            usr17b_ndx = get_rxt_ndx( 'usr_NO2_aer' )
          endif
          if ( full_ozone_chem ) then
             wrk(1:12) = (/ ox_l1_ndx, ox_l2_ndx, ox_l3_ndx, ox_l4_ndx, ox_l5_ndx, &
                  ox_l6_ndx, ox_l7_ndx, ox_l8_ndx, ox_l9_ndx, usr4_ndx, &
                  usr16_ndx, usr17_ndx /)
             if( any( wrk(1:12) < 1 ) ) then
                full_ozone_chem = .false.
             endif
          endif
          if ( reduced_ozone_chem ) then
             wrk(1:9) = (/ ox_l1_ndx, ox_l2_ndx, ox_l3_ndx, ox_l5_ndx, &
                  ox_l6_ndx, ox_l7_ndx, usr4_ndx, &
                  usr16_ndx, usr17_ndx /)
             if( any( wrk(1:9) < 1 ) ) then
                reduced_ozone_chem = .false.
             end if
          endif
          if ( middle_atm_chem ) then
             wrk(1:3) = (/ ox_l1_ndx, ox_l2_ndx, ox_l3_ndx /)
             if( any( wrk(1:3) < 1 ) ) then
                middle_atm_chem = .false.
             end if
          endif
          oh_ndx = get_spc_ndx( 'OH' )
          ho2_ndx = get_spc_ndx( 'HO2' )
          ch3o2_ndx = get_spc_ndx( 'CH3O2' )
          po2_ndx = get_spc_ndx( 'PO2' )
          ch3co3_ndx = get_spc_ndx( 'CH3CO3' )
          c2h5o2_ndx = get_spc_ndx( 'C2H5O2' )
          macro2_ndx = get_spc_ndx( 'MACRO2' )
          mco3_ndx = get_spc_ndx( 'MCO3' )
          c3h7o2_ndx = get_spc_ndx( 'C3H7O2' )
          ro2_ndx = get_spc_ndx( 'RO2' )
          xo2_ndx = get_spc_ndx( 'XO2' )
          no_ndx = get_spc_ndx( 'NO' )
          xno_ndx = get_spc_ndx( 'XNO' )
          no2_ndx = get_spc_ndx( 'NO2' )
          xno2_ndx = get_spc_ndx( 'XNO2' )
          no3_ndx = get_spc_ndx( 'NO3' )
          xno3_ndx = get_spc_ndx( 'XNO3' )
          n2o5_ndx = get_spc_ndx( 'N2O5' )
          xno2no3_ndx = get_spc_ndx( 'XNO2NO3' )
          no2xno3_ndx = get_spc_ndx( 'NO2XNO3' )
          c2h4_ndx = get_spc_ndx( 'C2H4' )
          c3h6_ndx = get_spc_ndx( 'C3H6' )
          isop_ndx = get_spc_ndx( 'ISOP' )
          isopo2_ndx = get_spc_ndx( 'ISOPO2' )
          mvk_ndx = get_spc_ndx( 'MVK' )
          c10h16_ndx = get_spc_ndx( 'C10H16' )
          tolo2_ndx = get_spc_ndx('TOLO2')
          terpo2_ndx =get_spc_ndx('TERPO2')
          alko2_ndx = get_spc_ndx('ALKO2')
          eneo2_ndx = get_spc_ndx('ENEO2')
          eo2_ndx = get_spc_ndx('EO2')
          meko2_ndx = get_spc_ndx('MEKO2')
          if ( full_ozone_chem ) then
             wrk(1:27) = (/ oh_ndx, ho2_ndx, ch3o2_ndx, po2_ndx, ch3co3_ndx, &
                  c2h5o2_ndx, macro2_ndx, mco3_ndx, c3h7o2_ndx, ro2_ndx, &
                  xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx, &
                  c2h4_ndx, c3h6_ndx, isop_ndx, isopo2_ndx, mvk_ndx, c10h16_ndx, &
                  tolo2_ndx, terpo2_ndx, alko2_ndx, eneo2_ndx, eo2_ndx, meko2_ndx /)
             if( any( wrk(1:27) < 1 ) ) then
                full_ozone_chem = .false.
             end if
          endif
          if ( reduced_ozone_chem ) then
             c2o3_ndx = get_spc_ndx( 'C2O3' )
             ole_ndx = get_spc_ndx( 'OLE' )
             wrk(1:12) = (/ oh_ndx, ho2_ndx, ch3o2_ndx, &
                  ole_ndx,c2o3_ndx, &
                  xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx, &
                  c2h4_ndx, isop_ndx /)
             if( any( wrk(1:12) < 1 ) ) then
                reduced_ozone_chem = .false.
             end if
          endif
       end if
    else
       reduced_ozone_chem = .false.
       full_ozone_chem = .false.
    end if has_o3_chem
    if ( reduced_ozone_chem .and. full_ozone_chem ) then
       write(iulog,*) 'can not have both full_ozone_chem and reduced_ozone_chem'
       call endrun
    endif
    do i = 1,clscnt4
       j = clsmap(i,4)
       call addfld( trim(solsym(j))//'_CHMP', '/cm3/s ', pver, 'I', 'chemical production rate', phys_decomp )
       call addfld( trim(solsym(j))//'_CHML', '/cm3/s ', pver, 'I', 'chemical loss rate',       phys_decomp )
    enddo
    call addfld('H_PEROX_CHMP', '/cm3/s ', pver, 'I', 'total ROOH production rate', phys_decomp ) !PJY changed "RO2" to "ROOH"
  end subroutine imp_slv_inti
  subroutine imp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, &
       xhnm, ncol, lchnk, ltrop, o3s_loss )
    ! CAM adapter for the pp_trop_mam3 implicit chemistry scheme.  Generated
    ! numerical work is in imp_sol_run; names, history, logging and init state
    ! remain here.
    use chem_mods, only : rxntot, extcnt, nzcnt, permute, cls_rxt_cnt
    use mo_tracname, only : solsym
    use ppgrid, only : pver
    use cam_history, only : outfld
    use perf_mod, only : t_startf, t_stopf
    use cam_abortutils, only : endrun
    use ap_imp_sol_trop_mam3_scheme, only : imp_sol_run
    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: lchnk
    real(r8), intent(in) :: delt
    real(r8), intent(in) :: reaction_rates(ncol,pver,max(1,rxntot))
    real(r8), intent(in) :: extfrc(ncol,pver,max(1,extcnt))
    real(r8), intent(in) :: het_rates(ncol,pver,max(1,gas_pcnst))
    real(r8), intent(inout) :: base_sol(ncol,pver,gas_pcnst)
    real(r8), intent(in) :: xhnm(ncol,pver)
    integer, intent(in) :: ltrop(ncol)
    real(r8), optional, intent(out) :: o3s_loss(ncol,pver)

    integer :: i, j
    integer :: errflg
    integer :: solver_failure_count
    real(r8) :: last_failure_dt
    character(len=512) :: errmsg
    real(r8) :: prod_out(ncol,pver,max(1,clscnt4))
    real(r8) :: loss_out(ncol,pver,max(1,clscnt4))
    real(r8) :: prod_hydrogen_peroxides_out(ncol,pver)

    if (present(o3s_loss)) then
       o3s_loss(:,:) = 0._r8
    end if

    call t_startf('ap_imp_sol_run')
    call imp_sol_run(ncol, pver, gas_pcnst, rxntot, extcnt, clscnt4, &
         nzcnt, itermax, base_sol, reaction_rates, het_rates, extfrc, delt, &
         xhnm, ltrop, clsmap(1:clscnt4,4), permute(1:clscnt4,4), &
         cls_rxt_cnt(1,4), epsilon, factor, small, prod_out, loss_out, &
         solver_failure_count, last_failure_dt, imp_sol_scheme_entry_count, &
         errmsg, errflg)
    call t_stopf('ap_imp_sol_run')

    if (errflg /= 0) then
       write(iulog,*) trim(errmsg)
       call endrun
    end if

    ! The normal pp_trop_mam3 path converges without this diagnostic.  Keep
    ! iulog at the CAM boundary if a future state exercises adaptive cuts.
    if (solver_failure_count > 0) then
       write(iulog,'('' imp_sol: decoupled solver recorded '',i6, &
            '' failed Newton steps; last dt = '',1p,e21.13)') &
            solver_failure_count, last_failure_dt
    end if

    prod_hydrogen_peroxides_out(:,:) = 0._r8
    do i = 1,clscnt4
       j = clsmap(i,4)
       call outfld( trim(solsym(j))//'_CHMP', prod_out(:,:,i), ncol, lchnk )
       call outfld( trim(solsym(j))//'_CHML', loss_out(:,:,i), ncol, lchnk )

       if ( trim(solsym(j)) == 'ALKOOH' &
        .or.trim(solsym(j)) == 'C2H5OOH' &
        .or.trim(solsym(j)) == 'CH3OOH' &
        .or.trim(solsym(j)) == 'CH3COOH' &
        .or.trim(solsym(j)) == 'CH3COOOH' &
        .or.trim(solsym(j)) == 'C3H7OOH' &
        .or.trim(solsym(j)) == 'EOOH' &
        .or.trim(solsym(j)) == 'ISOPOOH' &
        .or.trim(solsym(j)) == 'MACROOH' &
        .or.trim(solsym(j)) == 'MEKOOH' &
        .or.trim(solsym(j)) == 'POOH' &
        .or.trim(solsym(j)) == 'ROOH' &
        .or.trim(solsym(j)) == 'TERPOOH' &
        .or.trim(solsym(j)) == 'TOLOOH' &
        .or.trim(solsym(j)) == 'XOOH' ) then
          prod_hydrogen_peroxides_out(:,:) = &
               prod_hydrogen_peroxides_out(:,:) + prod_out(:,:,i)
       end if
    end do

    call outfld( 'H_PEROX_CHMP', prod_hydrogen_peroxides_out(:,:), ncol, lchnk )
  end subroutine imp_sol
end module mo_imp_sol
