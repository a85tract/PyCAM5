module mo_gas_phase_chemdr

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use shr_const_mod,    only : pi => shr_const_pi
  use constituents,     only : pcnst
  use cam_history,      only : fieldname_len
  use chem_mods,        only : phtcnt, rxntot, gas_pcnst
  use chem_mods,        only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map, extcnt
  use dust_model,       only : dust_names, ndust => dust_nbin
  use ppgrid,           only : pcols, pver
  use cam_logfile,      only : iulog
  use spmd_utils,       only : iam, masterproc
  use phys_control,     only : phys_getopts
  use carma_flags_mod,  only : carma_hetchem_feedback

  implicit none
  save

  private
  public :: gas_phase_chemdr, gas_phase_chemdr_inti 
  public :: map2chm

  integer :: map2chm(pcnst) = 0           ! index map to/from chemistry/constituents list

  integer :: synoz_ndx, so4_ndx, h2o_ndx, o2_ndx, o_ndx, hno3_ndx, hcl_ndx, dst_ndx, cldice_ndx
  integer :: o3_ndx, o3s_ndx
  integer :: het1_ndx
  integer :: ndx_cldfr, ndx_cmfdqr, ndx_nevapr, ndx_cldtop, ndx_prain
  integer :: ndx_h2so4
!
! CCMI
!
  integer :: st80_25_ndx
  integer :: st80_25_tau_ndx
  integer :: aoa_nh_ndx
  integer :: aoa_nh_ext_ndx
  integer :: nh_5_ndx
  integer :: nh_50_ndx
  integer :: nh_50w_ndx
  integer :: sad_pbf_ndx
  integer :: cb1_ndx,cb2_ndx,oc1_ndx,oc2_ndx,dst1_ndx,dst2_ndx,sslt1_ndx,sslt2_ndx
  integer :: soa_ndx,soai_ndx,soam_ndx,soat_ndx,soab_ndx,soax_ndx

  character(len=fieldname_len),dimension(rxntot-phtcnt) :: rxn_names
  character(len=fieldname_len),dimension(phtcnt)        :: pht_names
  character(len=fieldname_len),dimension(rxt_tag_cnt)   :: tag_names
  character(len=fieldname_len),dimension(extcnt)        :: extfrc_name

  logical :: pm25_srf_diag
  logical :: pm25_srf_diag_soa
  logical :: gas_phase_chemdr_use_native_impl = .false.
  logical :: gas_phase_chemdr_impl_selected = .false.

contains

  subroutine gas_phase_chemdr_inti()

    use mo_chem_utls,      only : get_spc_ndx, get_extfrc_ndx, get_inv_ndx, get_rxt_ndx
    use cam_history,       only : addfld,phys_decomp
    use cam_abortutils,    only : endrun
    use mo_chm_diags,      only : chm_diags_inti
    use mo_tracname,       only : solsym
    use constituents,      only : cnst_get_ind
    use physics_buffer,    only : pbuf_get_index
    use rate_diags,        only : rate_diags_init
    use rad_constituents,  only : rad_cnst_get_info

    implicit none

    character(len=3) :: string
    integer          :: n, m, err
    logical          :: history_aerosol      ! Output the MAM aerosol tendencies

    !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out = history_aerosol )
   
    ndx_h2so4 = get_spc_ndx('H2SO4')
!
! CCMI
!
    st80_25_ndx     = get_spc_ndx   ('ST80_25')
    st80_25_tau_ndx = get_rxt_ndx   ('ST80_25_tau')
    aoa_nh_ndx      = get_spc_ndx   ('AOA_NH')
    aoa_nh_ext_ndx  = get_extfrc_ndx('AOA_NH')
    nh_5_ndx        = get_spc_ndx('NH_5')
    nh_50_ndx       = get_spc_ndx('NH_50')
    nh_50w_ndx      = get_spc_ndx('NH_50W')
!
    cb1_ndx         = get_spc_ndx('CB1')
    cb2_ndx         = get_spc_ndx('CB2')
    oc1_ndx         = get_spc_ndx('OC1')
    oc2_ndx         = get_spc_ndx('OC2')
    dst1_ndx        = get_spc_ndx('DST01')
    dst2_ndx        = get_spc_ndx('DST02')
    sslt1_ndx       = get_spc_ndx('SSLT01')
    sslt2_ndx       = get_spc_ndx('SSLT02')
    soa_ndx         = get_spc_ndx('SOA')
    soam_ndx         = get_spc_ndx('SOAM')
    soai_ndx         = get_spc_ndx('SOAI')
    soat_ndx         = get_spc_ndx('SOAT')
    soab_ndx         = get_spc_ndx('SOAB')
    soax_ndx         = get_spc_ndx('SOAX')

    pm25_srf_diag = cb1_ndx>0 .and. cb2_ndx>0 .and. oc1_ndx>0 .and. oc2_ndx>0 &
              .and. dst1_ndx>0 .and. dst2_ndx>0 .and. sslt1_ndx>0 .and. sslt2_ndx>0 &
              .and. soa_ndx>0 

    pm25_srf_diag_soa = cb1_ndx>0 .and. cb2_ndx>0 .and. oc1_ndx>0 .and. oc2_ndx>0 &
              .and. dst1_ndx>0 .and. dst2_ndx>0 .and. sslt1_ndx>0 .and. sslt2_ndx>0 &
              .and. soam_ndx>0 .and. soai_ndx>0 .and. soat_ndx>0 .and. soab_ndx>0 .and. soax_ndx>0
    
    if ( pm25_srf_diag .or. pm25_srf_diag_soa) then
       call addfld('PM25_SRF','kg/kg',1,'I','bottom layer PM2.5 mixing ratio', phys_decomp )
    endif
    call addfld('U_SRF','m/s',1,'I','bottom layer wind velocity', phys_decomp )
    call addfld('V_SRF','m/s',1,'I','bottom layer wind velocity', phys_decomp )
    call addfld('Q_SRF','kg/kg',1,'I','bottom layer specific humidity', phys_decomp )
!
    call addfld('O3S_LOSS','mol/mol',pver,'I','O3S loss rate', phys_decomp )
!
    het1_ndx= get_rxt_ndx('het1')
    o3_ndx  = get_spc_ndx('O3')
    o3s_ndx = get_spc_ndx('O3S')
    o_ndx   = get_spc_ndx('O')
    o2_ndx  = get_spc_ndx('O2')
    so4_ndx = get_spc_ndx('SO4')
    h2o_ndx = get_spc_ndx('H2O')
    hno3_ndx = get_spc_ndx('HNO3')
    hcl_ndx  = get_spc_ndx('HCL')
    dst_ndx = get_spc_ndx( dust_names(1) )
    synoz_ndx = get_extfrc_ndx( 'SYNOZ' )
    call cnst_get_ind( 'CLDICE', cldice_ndx )

    do m = 1,extcnt
       WRITE(UNIT=string, FMT='(I2.2)') m
       extfrc_name(m) = 'extfrc_'// trim(string)
       call addfld( extfrc_name(m), ' ', pver, 'I', 'ext frcing', phys_decomp )
    end do

    do n = 1,rxt_tag_cnt
       tag_names(n) = trim(rxt_tag_lst(n))
       if (n<=phtcnt) then
          call addfld( tag_names(n), '/s ', pver, 'I', 'photolysis rate', phys_decomp )
       else
          call addfld( tag_names(n), '/cm3/s ', pver, 'I', 'reaction rate', phys_decomp )
       endif
    enddo

    do n = 1,phtcnt
       WRITE(UNIT=string, FMT='(I3.3)') n
       pht_names(n) = 'J_' // trim(string)
       call addfld( pht_names(n), '/s ', pver, 'I', 'photolysis rate', phys_decomp )
    enddo

    do n = 1,rxntot-phtcnt
       WRITE(UNIT=string, FMT='(I3.3)') n
       rxn_names(n) = 'R_' // trim(string)
       call addfld( rxn_names(n), '/cm3/s ', pver, 'I', 'reaction rate', phys_decomp )
    enddo

    call addfld( 'DTCBS',   ' ',  1, 'I','photolysis diagnostic black carbon OD', phys_decomp )
    call addfld( 'DTOCS',   ' ',  1, 'I','photolysis diagnostic organic carbon OD', phys_decomp )
    call addfld( 'DTSO4',   ' ',  1, 'I','photolysis diagnostic SO4 OD', phys_decomp )
    call addfld( 'DTSOA',   ' ',  1, 'I','photolysis diagnostic SOA OD', phys_decomp )
    call addfld( 'DTANT',   ' ',  1, 'I','photolysis diagnostic NH4SO4 OD', phys_decomp )
    call addfld( 'DTSAL',   ' ',  1, 'I','photolysis diagnostic salt OD', phys_decomp )
    call addfld( 'DTDUST',  ' ',  1, 'I','photolysis diagnostic dust OD', phys_decomp )
    call addfld( 'DTTOTAL', ' ',  1, 'I','photolysis diagnostic total aerosol OD', phys_decomp )   
    call addfld( 'FRACDAY', ' ',  1, 'I','photolysis diagnostic fraction of day', phys_decomp )

    call addfld( 'QDSAD', '/s ', pver, 'I', 'water vapor sad delta', phys_decomp )
    call addfld( 'SAD', 'cm2/cm3 ', pver, 'I', 'sulfate aerosol SAD', phys_decomp )
    call addfld( 'SAD_SULFC', 'cm2/cm3 ', pver, 'I', 'chemical sulfate aerosol SAD', phys_decomp )
    call addfld( 'SAD_SAGE', 'cm2/cm3 ', pver, 'I', 'SAGE sulfate aerosol SAD', phys_decomp )
    call addfld( 'SAD_LNAT', 'cm2/cm3 ', pver, 'I', 'large-mode NAT aerosol SAD', phys_decomp )
    call addfld( 'SAD_ICE', 'cm2/cm3 ', pver, 'I', 'water-ice aerosol SAD', phys_decomp )
    call addfld( 'RAD_SULFC', 'cm ', pver, 'I', 'chemical sad sulfate', phys_decomp )
    call addfld( 'RAD_LNAT', 'cm ', pver, 'I', 'large nat radius', phys_decomp )
    call addfld( 'RAD_ICE', 'cm ', pver, 'I', 'sad ice', phys_decomp )
    call addfld( 'SAD_TROP', 'cm2/cm3 ', pver, 'I', 'tropospheric aerosol SAD', phys_decomp )
    call addfld( 'QDSETT', '/s ', pver, 'I', 'water vapor settling delta', phys_decomp )
    call addfld( 'QDCHEM', '/s ', pver, 'I', 'water vapor chemistry delta', phys_decomp)
    call addfld( 'HNO3_TOTAL', 'mol/mol', pver, 'I', 'total HNO3', phys_decomp )
    call addfld( 'HNO3_STS',   'mol/mol', pver, 'I', 'STS condensed HNO3', phys_decomp )
    call addfld( 'HNO3_NAT',   'mol/mol', pver, 'I', 'NAT condensed HNO3', phys_decomp )
    call addfld( 'HNO3_GAS',   'mol/mol', pver, 'I', 'gas-phase hno3', phys_decomp )
    call addfld( 'H2O_GAS',    'mol/mol', pver, 'I', 'gas-phase h2o', phys_decomp )
    call addfld( 'HCL_TOTAL',  'mol/mol', pver, 'I', 'total hcl', phys_decomp )
    call addfld( 'HCL_GAS',    'mol/mol', pver, 'I', 'gas-phase hcl', phys_decomp )
    call addfld( 'HCL_STS',    'mol/mol', pver, 'I', 'STS condensed HCL', phys_decomp )

    if (het1_ndx>0) then
       call addfld( 'het1_total', '/s', pver, 'I', 'total N2O5 + H2O het rate constant', phys_decomp )
    endif
    call addfld( 'SZA', 'degrees', 1, 'I', 'solar zenith angle', phys_decomp )

    call chm_diags_inti()
    call rate_diags_init()

!-----------------------------------------------------------------------
! get pbuf indicies
!-----------------------------------------------------------------------
    ndx_cldfr  = pbuf_get_index('CLD')
    ndx_cmfdqr = pbuf_get_index('RPRDTOT')
    ndx_nevapr = pbuf_get_index('NEVAPR')
    ndx_prain  = pbuf_get_index('PRAIN')
    ndx_cldtop = pbuf_get_index('CLDTOP')

    sad_pbf_ndx= pbuf_get_index('VOLC_SAD',errcode=err) ! prescribed  strat aerosols (volcanic)
    if (.not.sad_pbf_ndx>0) sad_pbf_ndx = pbuf_get_index('SADSULF',errcode=err) ! CARMA's version of strat aerosols
    
  end subroutine gas_phase_chemdr_inti


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine gas_phase_chemdr(lchnk, ncol, imozart, q, &
                              phis, zm, zi, calday, &
                              tfld, pmid, pdel, pint,  &
                              cldw, troplev, &
                              ncldwtr, ufld, vfld,  &
                              delt, ps, xactive_prates, &
                              fsds, ts, asdir, ocnfrac, icefrac, &
                              precc, precl, snowhland, ghg_chem, latmapback, &
                              chem_name, drydepflx, cflx, qtend, pbuf)

    !-----------------------------------------------------------------------
    !     ... Chem_solver advances the volumetric mixing ratio
    !         forward one time step via a combination of explicit,
    !         ebi, hov, fully implicit, and/or rodas algorithms.
    !-----------------------------------------------------------------------

    use chem_mods,         only : nabscol, nfs, indexm
    use physconst,         only : rga
    use mo_photo,          only : set_ub_col, setcol, table_photo, xactive_photo
    use mo_exp_sol,        only : exp_sol
    use mo_imp_sol,        only : imp_sol
    use mo_setrxt,         only : setrxt
    use mo_adjrxt,         only : adjrxt
    use mo_phtadj,         only : phtadj
    use llnl_O1D_to_2OH_adj,only : O1D_to_2OH_adj
    use mo_usrrxt,         only : usrrxt
    use mo_setinv,         only : setinv
    use mo_negtrc,         only : negtrc
    use mo_sulf,           only : sulf_interp
    use mo_lightning,      only : prod_no
    use mo_setext,         only : setext
    use mo_sethet,         only : sethet
    use mo_drydep,         only : drydep, set_soilw
    use seq_drydep_mod,    only : DD_XLND, DD_XATM, DD_TABL, drydep_method
    use mo_fstrat,         only : set_fstrat_vals, set_fstrat_h2o
    use noy_ubc,           only : noy_ubc_set
    use mo_flbc,           only : flbc_set
    use phys_grid,         only : get_rlat_all_p, get_rlon_all_p, get_lat_all_p, get_lon_all_p
    use mo_mean_mass,      only : set_mean_mass
    use cam_history,       only : outfld
    use wv_saturation,     only : qsat
    use constituents,      only : cnst_mw
    use mo_drydep,         only : has_drydep
    use time_manager,      only : get_ref_date
    use mo_ghg_chem,       only : ghg_chem_set_rates, ghg_chem_set_flbc
    use mo_sad,            only : sad_strat_calc
    use charge_neutrality, only : charge_balance
    use mo_strato_rates,   only : ratecon_sfstrat
    use mo_aero_settling,  only : strat_aer_settling
    use shr_orb_mod,       only : shr_orb_decl
    use cam_control_mod,   only : lambm0, eccen, mvelpp, obliqr
    use mo_strato_rates,   only : has_strato_chem
    use short_lived_species,only: set_short_lived_species,get_short_lived_species
    use mo_chm_diags,      only : chm_diags, het_diags
    use perf_mod,          only : t_startf, t_stopf
    use mo_neu_wetdep,     only : do_neu_wetdep
    use physics_buffer,    only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use infnan,            only : nan, assignment(=)
    use rate_diags,        only : rate_diags_calc
    use mo_mass_xforms,    only : mmr2vmr, vmr2mmr, h2o_to_vmr, mmr2vmri
!
! LINOZ
!
    use lin_strat_chem,    only : do_lin_strat_chem, lin_strat_chem_solve
    use linoz_data,        only : has_linoz_data
!
! for aqueous chemistry and aerosol growth
!
    use aero_model,        only : aero_model_gasaerexch

    use aero_model,        only : aero_model_strat_surfarea

    implicit none

    !-----------------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------------
    integer,        intent(in)    :: lchnk                          ! chunk index
    integer,        intent(in)    :: ncol                           ! number columns in chunk
    integer,        intent(in)    :: imozart                        ! gas phase start index in q
    real(r8),       intent(in)    :: delt                           ! timestep (s)
    real(r8),       intent(in)    :: calday                         ! day of year
    real(r8),       intent(in)    :: ps(pcols)                      ! surface pressure
    real(r8),       intent(in)    :: phis(pcols)                    ! surface geopotential
    real(r8),       intent(in)    :: tfld(pcols,pver)               ! midpoint temperature (K)
    real(r8),       intent(in)    :: pmid(pcols,pver)               ! midpoint pressures (Pa)
    real(r8),       intent(in)    :: pdel(pcols,pver)               ! pressure delta about midpoints (Pa)
    real(r8),       intent(in)    :: ufld(pcols,pver)               ! zonal velocity (m/s)
    real(r8),       intent(in)    :: vfld(pcols,pver)               ! meridional velocity (m/s)
    real(r8),       intent(in)    :: cldw(pcols,pver)               ! cloud water (kg/kg)
    real(r8),       intent(in)    :: ncldwtr(pcols,pver)            ! droplet number concentration (#/kg)
    real(r8),       intent(in)    :: zm(pcols,pver)                 ! midpoint geopotential height above the surface (m)
    real(r8),       intent(in)    :: zi(pcols,pver+1)               ! interface geopotential height above the surface (m)
    real(r8),       intent(in)    :: pint(pcols,pver+1)             ! interface pressures (Pa)
    real(r8),       intent(in)    :: q(pcols,pver,pcnst)            ! species concentrations (kg/kg)
    logical,        intent(in)    :: xactive_prates
    real(r8),       intent(in)    :: fsds(pcols)                    ! longwave down at sfc
    real(r8),       intent(in)    :: icefrac(pcols)                 ! sea-ice areal fraction
    real(r8),       intent(in)    :: ocnfrac(pcols)                 ! ocean areal fraction
    real(r8),       intent(in)    :: asdir(pcols)                   ! albedo: shortwave, direct
    real(r8),       intent(in)    :: ts(pcols)                      ! sfc temp (merged w/ocean if coupled)
    real(r8),       intent(in)    :: precc(pcols)                   !
    real(r8),       intent(in)    :: precl(pcols)                   !
    real(r8),       intent(in)    :: snowhland(pcols)               !
    logical,        intent(in)    :: ghg_chem 
    integer,        intent(in)    :: latmapback(pcols)
    character(len=*), intent(in)  :: chem_name
    integer,         intent(in) ::  troplev(pcols)

    real(r8),       intent(inout) :: qtend(pcols,pver,pcnst)        ! species tendencies (kg/kg/s)
    real(r8),       intent(inout) :: cflx(pcols,pcnst)              ! constituent surface flux (kg/m^2/s)
    real(r8),       intent(out)   :: drydepflx(pcols,pcnst)         ! dry deposition flux (kg/m^2/s)

    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    real(r8), parameter :: m2km  = 1.e-3_r8
    real(r8), parameter :: Pa2mb = 1.e-2_r8

    real(r8),       pointer    :: prain(:,:)
    real(r8),       pointer    :: nevapr(:,:)
    real(r8),       pointer    :: cmfdqr(:,:)
    real(r8),       pointer    :: cldfr(:,:)
    real(r8),       pointer    :: cldtop(:)

    integer      ::  i, k, m, n
    integer      ::  tim_ndx
    real(r8)     ::  delt_inverse
    real(r8)     ::  esfact
    integer      ::  latndx(pcols)                         ! chunk lat indicies
    integer      ::  lonndx(pcols)                         ! chunk lon indicies
    real(r8)     ::  invariants(ncol,pver,nfs)
    real(r8)     ::  col_dens(ncol,pver,max(1,nabscol))    ! column densities (molecules/cm^2)
    real(r8)     ::  col_delta(ncol,0:pver,max(1,nabscol)) ! layer column densities (molecules/cm^2)
    real(r8)     ::  extfrc(ncol,pver,max(1,extcnt))
    real(r8)     ::  vmr(ncol,pver,gas_pcnst)         ! xported species (vmr)
    real(r8)     ::  reaction_rates(ncol,pver,max(1,rxntot))      ! reaction rates
    real(r8)     ::  depvel(ncol,gas_pcnst)                ! dry deposition velocity (cm/s)
    real(r8)     ::  het_rates(ncol,pver,max(1,gas_pcnst)) ! washout rate (1/s)
    real(r8), dimension(ncol,pver) :: &
         h2ovmr, &                                         ! water vapor volume mixing ratio
         mbar, &                                           ! mean wet atmospheric mass ( amu )
         zmid, &                                           ! midpoint geopotential in km
         zmidr, &                                          ! midpoint geopotential in km realitive to surf
         sulfate, &                                        ! trop sulfate aerosols
         pmb                                               ! pressure at midpoints ( hPa )
    real(r8), dimension(ncol,pver) :: &
         cwat, &                                           ! cloud water mass mixing ratio (kg/kg)
         wrk
    real(r8), dimension(ncol,pver+1) :: &
         zintr                                              ! interface geopotential in km realitive to surf
    real(r8), dimension(ncol,pver+1) :: &
         zint                                              ! interface geopotential in km
    real(r8), dimension(ncol)  :: &
         zen_angle, &                                      ! solar zenith angles
         zsurf, &                                          ! surface height (m)
         rlats, rlons                                      ! chunk latitudes and longitudes (radians)
    real(r8) :: sza(ncol)                                  ! solar zenith angles (degrees)
    real(r8), parameter :: rad2deg = 180._r8/pi                ! radians to degrees conversion factor
    real(r8) :: relhum(ncol,pver)                          ! relative humidity
    real(r8) :: satv(ncol,pver)                            ! wrk array for relative humidity
    real(r8) :: satq(ncol,pver)                            ! wrk array for relative humidity

    integer                   :: j,wd_index
    real(r8), dimension(ncol) :: noy_wk, sox_wk, nhx_wk
    integer                   ::  ltrop_sol(pcols)         ! tropopause vertical index used in chem solvers
    real(r8), pointer         ::  strato_sad(:,:)          ! stratospheric sad (1/cm)

    real(r8)                  ::  sad_total(pcols,pver)    ! total trop. sad (cm^2/cm^3)

    real(r8) :: tvs(pcols)
    integer  :: ncdate,yr,mon,day,sec
    real(r8) :: wind_speed(pcols)        ! surface wind speed (m/s)
    logical, parameter :: dyn_soilw = .false.
    logical  :: table_soilw
    real(r8) :: soilw(pcols)
    real(r8) :: prect(pcols)
    real(r8) :: sflx(pcols,gas_pcnst)
    real(r8) :: dust_vmr(ncol,pver,ndust)
    real(r8) :: noy(ncol,pver)
    real(r8) :: sox(ncol,pver)
    real(r8) :: nhx(ncol,pver)
    real(r8) :: dt_diag(pcols,8)               ! od diagnostics
    real(r8) :: fracday(pcols)                 ! fraction of day
    real(r8) :: o2mmr(ncol,pver)               ! o2 concentration (kg/kg)
    real(r8) :: ommr(ncol,pver)                ! o concentration (kg/kg)
    real(r8) :: mmr(pcols,pver,gas_pcnst)      ! chem working concentrations (kg/kg)
    real(r8) :: mmr_new(pcols,pver,gas_pcnst)      ! chem working concentrations (kg/kg)
    real(r8) :: hno3_gas(ncol,pver)            ! hno3 gas phase concentration (mol/mol)
    real(r8) :: hno3_cond(ncol,pver,2)         ! hno3 condensed phase concentration (mol/mol)
    real(r8) :: hcl_gas(ncol,pver)             ! hcl gas phase concentration (mol/mol)
    real(r8) :: hcl_cond(ncol,pver)            ! hcl condensed phase concentration (mol/mol)
    real(r8) :: h2o_gas(ncol,pver)             ! h2o gas phase concentration (mol/mol)
    real(r8) :: h2o_cond(ncol,pver)            ! h2o condensed phase concentration (mol/mol)
    real(r8) :: cldice(pcols,pver)             ! cloud water "ice" (kg/kg)
    real(r8) :: radius_strat(ncol,pver,3)      ! radius of sulfate, nat, & ice ( cm )
    real(r8) :: sad_strat(ncol,pver,3)         ! surf area density of sulfate, nat, & ice ( cm^2/cm^3 )
    real(r8) :: mmr_tend(pcols,pver,gas_pcnst) ! chemistry species tendencies (kg/kg/s)
    real(r8) :: qh2o(pcols,pver)               ! specific humidity (kg/kg)
    real(r8) :: delta

  ! for aerosol formation....  
    real(r8) :: del_h2so4_gasprod(ncol,pver)
    real(r8) :: vmr0(ncol,pver,gas_pcnst)

!
! CCMI
!
    real(r8) :: xlat
    real(r8) :: pm25(ncol)
    real(r8), dimension(ncol,pver) :: o3s_loss             ! tropospheric ozone loss for o3s
!
! jfl
!
!

    ! initialize to NaN to hopefully catch user defined rxts that go unset
    reaction_rates(:,:,:) = nan

    delt_inverse = 1._r8 / delt
    !-----------------------------------------------------------------------      
    !        ... Get chunck latitudes and longitudes
    !-----------------------------------------------------------------------      
    call get_lat_all_p( lchnk, ncol, latndx )
    call get_lon_all_p( lchnk, ncol, lonndx )
    call get_rlat_all_p( lchnk, ncol, rlats )
    call get_rlon_all_p( lchnk, ncol, rlons )
    tim_ndx = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, ndx_prain,      prain,  start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cldfr,        cldfr, start=(/1,1,tim_ndx/), kount=(/ncol,pver,1/) )
    call pbuf_get_field(pbuf, ndx_cmfdqr,     cmfdqr, start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_nevapr,     nevapr, start=(/1,1/), kount=(/ncol,pver/))
    call pbuf_get_field(pbuf, ndx_cldtop,     cldtop )

    call gas_phase_chemdr_select_impl()

    !-----------------------------------------------------------------------      
    !        ... Calculate cosine of zenith angle
    !            then cast back to angle (radians)
    !-----------------------------------------------------------------------      
    call zenith( calday, rlats, rlons, zen_angle, ncol )
    call gas_phase_chemdr_prepare_sza(ncol, rad2deg, zen_angle, sza)
    call outfld( 'SZA',   sza,    ncol, lchnk )

    !-----------------------------------------------------------------------      
    !        ... Xform geopotential height from m to km 
    !            and pressure from Pa to mb
    !-----------------------------------------------------------------------      
    call gas_phase_chemdr_prepare_state(ncol, phis, zi, zm, pmid, zsurf, zintr, zmidr, zmid, zint, pmb)

    !-----------------------------------------------------------------------      
    !        ... map incoming concentrations to working array
    !-----------------------------------------------------------------------      
    call gas_phase_chemdr_load_mmr(ncol, q, mmr)

    call get_short_lived_species( mmr, lchnk, ncol, pbuf )

    !-----------------------------------------------------------------------      
    !        ... Set atmosphere mean mass
    !-----------------------------------------------------------------------      
    call set_mean_mass( ncol, mmr, mbar )

    !-----------------------------------------------------------------------      
    !        ... Xform from mmr to vmr
    !-----------------------------------------------------------------------      
    call mmr2vmr( mmr, vmr, mbar, ncol )

!
! CCMI
!
! reset STE tracer to specific vmr of 200 ppbv
!
    if ( st80_25_ndx > 0 ) then 
       call gas_phase_chemdr_reset_ste_tracer(ncol, st80_25_ndx, 80.e+2_r8, 200.e-9_r8, pmid, vmr)
    end if
!
! reset AOA_NH, NH_5, NH_50, NH_50W surface mixing ratios between 30N and 50N
!
    if ( aoa_nh_ndx>0 .and. nh_5_ndx>0 .and. nh_50_ndx>0 .and. nh_50w_ndx>0 ) then
      do j=1,ncol
        xlat = rlats(j)*rad2deg              ! convert to degrees
        if ( xlat >= 30._r8 .and. xlat <= 50._r8 ) then
           vmr(j,pver,nh_5_ndx)   = 100.e-9_r8
           vmr(j,pver,nh_50_ndx)  = 100.e-9_r8
           vmr(j,pver,nh_50w_ndx) = 100.e-9_r8
           vmr(j,pver,aoa_nh_ndx) = 0._r8
        end if
      end do
    end if

    if (h2o_ndx>0) then
       !-----------------------------------------------------------------------      
       !        ... store water vapor in wrk variable
       !-----------------------------------------------------------------------      
       call gas_phase_chemdr_load_h2o_fields(ncol, h2o_ndx, mmr, vmr, qh2o, h2ovmr)
    else
       qh2o(:ncol,:) = q(:ncol,:,1)
       !-----------------------------------------------------------------------      
       !        ... Xform water vapor from mmr to vmr and set upper bndy values
       !-----------------------------------------------------------------------      
       call h2o_to_vmr( q(:,:,1), h2ovmr, mbar, ncol )

       call set_fstrat_h2o( h2ovmr, pmid, troplev, calday, ncol, lchnk )

    endif

    !-----------------------------------------------------------------------      
    !        ... force ion/electron balance
    !-----------------------------------------------------------------------      
    call charge_balance( ncol, vmr )

    !-----------------------------------------------------------------------      
    !        ... Set the "invariants"
    !-----------------------------------------------------------------------  
    call setinv( invariants, tfld, h2ovmr, vmr, pmid, ncol, lchnk, pbuf )

    !-----------------------------------------------------------------------      
    !        ... stratosphere aerosol surface area
    !-----------------------------------------------------------------------  
    if (sad_pbf_ndx>0) then
       call pbuf_get_field(pbuf, sad_pbf_ndx, strato_sad)
    else
       allocate(strato_sad(pcols,pver))
       strato_sad(:,:) = 0._r8

       ! Prognostic modal stratospheric sulfate: compute dry strato_sad
       call aero_model_strat_surfarea( ncol, mmr, pmid, tfld, troplev, pbuf, strato_sad )

    endif

    stratochem: if ( has_strato_chem ) then
       !-----------------------------------------------------------------------      
       !        ... initialize condensed and gas phases; all hno3 to gas
       !-----------------------------------------------------------------------    
       call gas_phase_chemdr_init_stratchem_state(ncol, hno3_ndx, hcl_ndx, cldice_ndx, vmr, h2ovmr, q, hcl_cond, &
            hcl_gas, hno3_gas, h2o_gas, wrk, cldice, hno3_cond)

       call mmr2vmri( cldice, h2o_cond, mbar, cnst_mw(cldice_ndx), ncol )

       !-----------------------------------------------------------------------      
       !        ... call SAD routine
       !-----------------------------------------------------------------------      
       call sad_strat_calc( lchnk, invariants(:ncol,:,indexm), pmb, tfld, hno3_gas, &
            hno3_cond, h2o_gas, h2o_cond, hcl_gas, hcl_cond, strato_sad(:ncol,:), radius_strat, &
            sad_strat, ncol, pbuf )

!      NOTE: output of total HNO3 is before vmr is set to gas-phase.
       call outfld( 'HNO3_TOTAL', vmr(:ncol,:,hno3_ndx), ncol ,lchnk )


       call gas_phase_chemdr_restore_strat_gases(ncol, hno3_ndx, h2o_ndx, delt_inverse, vmr, hno3_gas, h2o_gas, &
            h2ovmr, wrk)

       call outfld( 'QDSAD', wrk(:,:), ncol, lchnk )
!
       call outfld( 'SAD', strato_sad    (:ncol,:), ncol, lchnk )
       call outfld( 'SAD_SULFC',  sad_strat(:,:,1), ncol, lchnk )
       call outfld( 'SAD_LNAT',   sad_strat(:,:,2), ncol, lchnk )
       call outfld( 'SAD_ICE',    sad_strat(:,:,3), ncol, lchnk )
!
       call outfld( 'RAD_SULFC',  radius_strat(:,:,1), ncol, lchnk )
       call outfld( 'RAD_LNAT',   radius_strat(:,:,2), ncol, lchnk )
       call outfld( 'RAD_ICE',    radius_strat(:,:,3), ncol, lchnk )
!
       call outfld( 'HNO3_GAS',   vmr(:ncol,:,hno3_ndx), ncol, lchnk )
       call outfld( 'HNO3_STS',   hno3_cond(:,:,1), ncol, lchnk )
       call outfld( 'HNO3_NAT',   hno3_cond(:,:,2), ncol, lchnk )
!
       call outfld( 'HCL_TOTAL',  vmr(:ncol,:,hcl_ndx), ncol, lchnk )
       call outfld( 'HCL_GAS',    hcl_gas (:,:), ncol ,lchnk )
       call outfld( 'HCL_STS',    hcl_cond(:,:), ncol ,lchnk )

       !-----------------------------------------------------------------------      
       !        ... call aerosol reaction rates
       !-----------------------------------------------------------------------      
       call ratecon_sfstrat( invariants(:,:,indexm), pmid, tfld, &
            radius_strat(:,:,1), sad_strat(:,:,1), sad_strat(:,:,2), &
            sad_strat(:,:,3), h2ovmr, vmr, reaction_rates, ncol )

    endif stratochem

!      NOTE: For gas-phase solver only. 
!            ratecon_sfstrat needs total hcl.
    if (hcl_ndx>0) then
       call gas_phase_chemdr_restore_hcl_gas(ncol, hcl_ndx, vmr, hcl_gas)
    endif

    !-----------------------------------------------------------------------      
    !        ... Set the column densities at the upper boundary
    !-----------------------------------------------------------------------      
    call set_ub_col( col_delta, vmr, invariants, pint(:,1), pdel, ncol, lchnk)

    !-----------------------------------------------------------------------      
    !       ...  Set rates for "tabular" and user specified reactions
    !-----------------------------------------------------------------------      
    call setrxt( reaction_rates, tfld, invariants(1,1,indexm), ncol )
    
    call gas_phase_chemdr_zero_sulfate(ncol, sulfate)
    if ( .not. carma_hetchem_feedback ) then
       if( so4_ndx < 1 ) then ! get offline so4 field if not prognostic
          call sulf_interp( ncol, lchnk, sulfate )
       else
          sulfate(:,:) = vmr(:,:,so4_ndx)
       endif
    endif
    
    !-----------------------------------------------------------------
    ! ... zero out sulfate above tropopause
    !-----------------------------------------------------------------
    call gas_phase_chemdr_clip_sulfate(ncol, troplev, sulfate)

    !-----------------------------------------------------------------
    !	... compute the relative humidity
    !-----------------------------------------------------------------
    call qsat(tfld(:ncol,:), pmid(:ncol,:), satv, satq)

    call gas_phase_chemdr_compute_relhum(ncol, h2ovmr, satq, relhum)
    
    call gas_phase_chemdr_copy_cldw_to_cwat(ncol, cldw, cwat)

    call usrrxt( reaction_rates, tfld, tfld, tfld, invariants, h2ovmr, ps, &
                 pmid, invariants(:,:,indexm), sulfate, mmr, relhum, strato_sad, &
                 troplev, ncol, sad_total, cwat, mbar, pbuf )

    call outfld( 'SAD_TROP', sad_total(:ncol,:), ncol, lchnk )

    if (het1_ndx>0) then
       call outfld( 'het1_total', reaction_rates(:,:,het1_ndx), ncol, lchnk )
    endif

    if (ghg_chem) then
       call ghg_chem_set_rates( reaction_rates, latmapback, zen_angle, ncol, lchnk )
    endif

    do i = phtcnt+1,rxntot
       call outfld( rxn_names(i-phtcnt), reaction_rates(:,:,i), ncol, lchnk )
    enddo

    call adjrxt( reaction_rates, invariants, invariants(1,1,indexm), ncol )

    !-----------------------------------------------------------------------
    !        ... Compute the photolysis rates at time = t(n+1)
    !-----------------------------------------------------------------------      
    !-----------------------------------------------------------------------      
    !     	... Set the column densities
    !-----------------------------------------------------------------------      
    call setcol( col_delta, col_dens, vmr, pdel,  ncol )

    !-----------------------------------------------------------------------      
    !     	... Calculate the photodissociation rates
    !-----------------------------------------------------------------------      

    esfact = 1._r8
    call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr  , &
         delta, esfact )


    if ( xactive_prates ) then
       call gas_phase_chemdr_init_dust_vmr(ncol, dst_ndx, vmr, dust_vmr)

       !-----------------------------------------------------------------
       !	... compute the photolysis rates
       !-----------------------------------------------------------------
       call xactive_photo( reaction_rates, vmr, tfld, cwat, cldfr, &
            pmid, zmidr, col_dens, zen_angle, asdir, &
            invariants(1,1,indexm), ps, ts, &
            esfact, relhum, dust_vmr, dt_diag, fracday, ncol, lchnk )

       call outfld('DTCBS',   dt_diag(:ncol,1), ncol, lchnk )
       call outfld('DTOCS',   dt_diag(:ncol,2), ncol, lchnk )
       call outfld('DTSO4',   dt_diag(:ncol,3), ncol, lchnk )
       call outfld('DTANT',   dt_diag(:ncol,4), ncol, lchnk )
       call outfld('DTSAL',   dt_diag(:ncol,5), ncol, lchnk )
       call outfld('DTDUST',  dt_diag(:ncol,6), ncol, lchnk )
       call outfld('DTSOA',   dt_diag(:ncol,7), ncol, lchnk )
       call outfld('DTTOTAL', dt_diag(:ncol,8), ncol, lchnk )
       call outfld('FRACDAY', fracday(:ncol), ncol, lchnk )

    else
       !-----------------------------------------------------------------
       !	... lookup the photolysis rates from table
       !-----------------------------------------------------------------
       call table_photo( reaction_rates, pmid, pdel, tfld, zmid, zint, &
                         col_dens, zen_angle, asdir, cwat, cldfr, &
                         esfact, vmr, invariants, ncol, lchnk, pbuf )
    endif

    do i = 1,phtcnt
       call outfld( pht_names(i), reaction_rates(:ncol,:,i), ncol, lchnk )
       call outfld( tag_names(i), reaction_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    !-----------------------------------------------------------------------      
    !     	... Adjust the photodissociation rates
    !-----------------------------------------------------------------------  
    call O1D_to_2OH_adj( reaction_rates, invariants, invariants(:,:,indexm), ncol, tfld )
    call phtadj( reaction_rates, invariants, invariants(:,:,indexm), ncol )

    !-----------------------------------------------------------------------
    !        ... Compute the extraneous frcing at time = t(n+1)
    !-----------------------------------------------------------------------      
    if ( o2_ndx > 0 .and. o_ndx > 0 ) then
       call gas_phase_chemdr_load_oxygen_mmr(ncol, o2_ndx, o_ndx, mmr, o2mmr, ommr)
    endif
    !-----------------------------------------------------------------------
    !        ... Compute the extraneous frcing at time = t(n+1)
    !-----------------------------------------------------------------------      
    call setext( extfrc, zint, zintr, cldtop, &
                 zmid, lchnk, tfld, o2mmr, ommr, &
                 pmid, mbar, rlats, calday, ncol, rlons, pbuf )

    call gas_phase_chemdr_normalize_extfrc(ncol, extcnt, nfs, indexm, extfrc, invariants)
    do m = 1,extcnt
       call outfld( extfrc_name(m), extfrc(:ncol,:,m), ncol, lchnk )
    end do

    !-----------------------------------------------------------------------
    !        ... Form the washout rates
    !-----------------------------------------------------------------------      
    if ( do_neu_wetdep ) then
      het_rates = 0._r8
    else
      call sethet( het_rates, pmid, zmid, phis, tfld, &
                   cmfdqr, prain, nevapr, delt, invariants(:,:,indexm), &
                   vmr, ncol, lchnk )
      call het_diags( het_rates(:ncol,:,:), mmr(:ncol,:,:), pdel(:ncol,:), lchnk, ncol )
    end if
!
! CCMI
!
! set loss to below the tropopause only
!
    call gas_phase_chemdr_zero_st80_tau(ncol, rxntot, st80_25_tau_ndx, troplev, reaction_rates)

!

    do i = phtcnt+1,rxt_tag_cnt
       call outfld( tag_names(i), reaction_rates(:ncol,:,rxt_tag_map(i)), ncol, lchnk )
    enddo

    call gas_phase_chemdr_set_ltrop_sol(ncol, merge(1, 0, has_linoz_data), troplev, ltrop_sol)

    ! save h2so4 before gas phase chem (for later new particle nucleation)
    call gas_phase_chemdr_init_h2so4_gasprod(ncol, ndx_h2so4, vmr, del_h2so4_gasprod)

    call gas_phase_chemdr_store_vmr0(ncol, vmr, vmr0)

    !=======================================================================
    !        ... Call the class solution algorithms
    !=======================================================================
    !-----------------------------------------------------------------------
    !	... Solve for "Explicit" species
    !-----------------------------------------------------------------------
    call exp_sol( vmr, reaction_rates, het_rates, extfrc, delt, invariants(1,1,indexm), ncol, lchnk, ltrop_sol )

    !-----------------------------------------------------------------------
    !	... Solve for "Implicit" species
    !-----------------------------------------------------------------------
    if ( has_strato_chem ) call gas_phase_chemdr_copy_h2o_to_wrk(ncol, h2o_ndx, vmr, wrk)
    call t_startf('imp_sol')
    !
    call imp_sol( vmr, reaction_rates, het_rates, extfrc, delt, &
                  invariants(1,1,indexm), ncol, lchnk, ltrop_sol(:ncol), o3s_loss=o3s_loss )
    call t_stopf('imp_sol')

    if( h2o_ndx>0) call outfld( 'H2O_GAS',  vmr(1,1,h2o_ndx),  ncol ,lchnk )

!
! jfl : CCMI : implement O3S here because mo_fstrat is not called
!
    if ( o3_ndx > 0 .and. o3s_ndx > 0 ) then
       call gas_phase_chemdr_copy_o3_to_o3s_trop(ncol, troplev, o3_ndx, o3s_ndx, vmr)
       do i = 1,ncol
          vmr(i,troplev(i)+1:pver,o3s_ndx) = vmr(i,troplev(i)+1:pver,o3s_ndx) * exp(-delt*o3s_loss(i,troplev(i)+1:pver))
       enddo
       call outfld( 'O3S_LOSS',  o3s_loss,  ncol ,lchnk )
    end if

    ! save h2so4 change by gas phase chem (for later new particle nucleation)
    if (ndx_h2so4 > 0) then
       call gas_phase_chemdr_update_h2so4_gasprod(ncol, ndx_h2so4, vmr, del_h2so4_gasprod)
    endif

!
! Aerosol processes ...
!

    call aero_model_gasaerexch( imozart-1, ncol, lchnk, troplev, delt, reaction_rates, &
                                tfld, pmid, pdel, mbar, relhum, &
                                zm,  qh2o, cwat, cldfr, ncldwtr, &
                                invariants(:,:,indexm), invariants, del_h2so4_gasprod,  &
                                vmr0, vmr, pbuf )

    if ( has_strato_chem ) then 

       call gas_phase_chemdr_update_qdchem_wrk(ncol, h2o_ndx, delt_inverse, vmr, wrk)
       call outfld( 'QDCHEM',   wrk(:ncol,:),         ncol, lchnk )
       call outfld( 'HNO3_GAS', vmr(:ncol,:,hno3_ndx), ncol ,lchnk )

       !-----------------------------------------------------------------------      
       !         ... aerosol settling
       !             first settle hno3(2) using radius ice
       !             secnd settle hno3(3) using radius large nat
       !-----------------------------------------------------------------------      
       call gas_phase_chemdr_copy_h2o_to_wrk(ncol, h2o_ndx, vmr, wrk)
#ifdef ALT_SETTL
       where( h2o_cond(:,:) > 0._r8 )
          settl_rad(:,:) = radius_strat(:,:,3)
       elsewhere
          settl_rad(:,:) = 0._r8
       endwhere
       call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
            hno3_cond(1,1,2), settl_rad, ncol, lchnk, 1 )

       where( h2o_cond(:,:) == 0._r8 )
          settl_rad(:,:) = radius_strat(:,:,2)
       elsewhere
          settl_rad(:,:) = 0._r8
       endwhere
       call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
            hno3_cond(1,1,2), settl_rad, ncol, lchnk, 2 )
#else
       call strat_aer_settling( invariants(1,1,indexm), pmid, delt, zmid, tfld, &
            hno3_cond(1,1,2), radius_strat(1,1,2), ncol, lchnk, 2 )
#endif

       !-----------------------------------------------------------------------      
       !	... reform total hno3 and hcl = gas + all condensed
       !-----------------------------------------------------------------------      
!      NOTE: vmr for hcl and hno3 is gas-phase at this point.
!            hno3_cond(:,k,1) = STS; hno3_cond(:,k,2) = NAT
   
       call gas_phase_chemdr_reform_hno3_hcl(ncol, hno3_ndx, hcl_ndx, vmr, hno3_cond, hcl_cond)

       call gas_phase_chemdr_update_qdsett_wrk(ncol, h2o_ndx, delt_inverse, vmr, wrk)
       call outfld( 'QDSETT', wrk(:,:), ncol, lchnk )

    endif

!
! LINOZ
!
    if ( do_lin_strat_chem ) then
       call lin_strat_chem_solve( ncol, lchnk, vmr(:,:,o3_ndx), col_dens(:,:,1), tfld, zen_angle, pmid, delt, rlats, troplev )
    end if

    !-----------------------------------------------------------------------      
    !         ... Check for negative values and reset to zero
    !-----------------------------------------------------------------------      
    call negtrc( 'After chemistry ', vmr, ncol )

    !-----------------------------------------------------------------------      
    !         ... Set upper boundary mmr values
    !-----------------------------------------------------------------------      
    call set_fstrat_vals( vmr, pmid, pint, troplev, calday, ncol,lchnk )

    !-----------------------------------------------------------------------      
    !         ... Set fixed lower boundary mmr values
    !-----------------------------------------------------------------------      
    call flbc_set( vmr, ncol, lchnk, map2chm )

    !----------------------------------------------------------------------- 
    ! set NOy UBC     
    !-----------------------------------------------------------------------      
    call noy_ubc_set( lchnk, ncol, vmr )

    if ( ghg_chem ) then
       call ghg_chem_set_flbc( vmr, ncol )
    endif

    !-----------------------------------------------------------------------      
    !         ... Xform from vmr to mmr
    !-----------------------------------------------------------------------      
    call vmr2mmr( vmr, mmr_tend, mbar, ncol )

    call set_short_lived_species( mmr_tend, lchnk, ncol, pbuf )

    !-----------------------------------------------------------------------      
    !         ... Form the tendencies
    !----------------------------------------------------------------------- 
    call gas_phase_chemdr_finalize_tendencies(ncol, delt_inverse, mmr, mmr_tend, mmr_new, qtend)

    call gas_phase_chemdr_compute_tvs(ncol, tfld, qh2o, tvs)

    call gas_phase_chemdr_zero_sflx(sflx)
    call get_ref_date(yr, mon, day, sec)
    ncdate = yr*10000 + mon*100 + day
    call gas_phase_chemdr_compute_wind_speed(ncol, ufld, vfld, wind_speed)
    call gas_phase_chemdr_compute_prect(ncol, precc, precl, prect)

    if ( drydep_method == DD_XLND ) then
       soilw = -99
       call drydep( ocnfrac, icefrac, ncdate, ts, ps,  &
            wind_speed, qh2o(:,pver), tfld(:,pver), pmid(:,pver), prect, &
            snowhland, fsds, depvel, sflx, mmr, &
            tvs, soilw, relhum(:,pver:pver), ncol, lonndx, latndx, lchnk )
    else if ( drydep_method == DD_XATM ) then
       table_soilw = has_drydep( 'H2' ) .or. has_drydep( 'CO' )
       if( .not. dyn_soilw .and. table_soilw ) then
          call set_soilw( soilw, lchnk, calday )
       end if
       call drydep( ncdate, ts, ps,  &
            wind_speed, qh2o(:,pver), tfld(:,pver), pmid(:,pver), prect, &
            snowhland, fsds, depvel, sflx, mmr, &
            tvs, soilw, relhum(:,pver:pver), ncol, lonndx, latndx, lchnk )
    else if ( drydep_method == DD_TABL ) then
       call drydep( calday, ts, zen_angle, &
            depvel, sflx, mmr, pmid(:,pver), &
            tvs, ncol, icefrac, ocnfrac, lchnk )
    endif

    call gas_phase_chemdr_store_drydep(ncol, sflx, cflx, drydepflx)

    call chm_diags( lchnk, ncol, vmr(:ncol,:,:), mmr_new(:ncol,:,:), &
                    reaction_rates(:ncol,:,:), invariants(:ncol,:,:), depvel(:ncol,:),  sflx(:ncol,:), &
                    mmr_tend(:ncol,:,:), pdel(:ncol,:), pmid(:ncol,:), troplev(:ncol), pbuf  )

    call rate_diags_calc( reaction_rates(:,:,:), vmr(:,:,:), invariants(:,:,indexm), ncol, lchnk )
!
! jfl
!
! surface vmr
!
    if ( pm25_srf_diag ) then
       pm25(:ncol) = mmr_new(:ncol,pver,cb1_ndx)   &
            + mmr_new(:ncol,pver,cb2_ndx)   &
            + mmr_new(:ncol,pver,oc1_ndx)   &
            + mmr_new(:ncol,pver,oc2_ndx)   &
            + mmr_new(:ncol,pver,dst1_ndx)  &
            + mmr_new(:ncol,pver,dst2_ndx)  &
            + mmr_new(:ncol,pver,sslt1_ndx) &
            + mmr_new(:ncol,pver,sslt2_ndx) &
            + mmr_new(:ncol,pver,soa_ndx)   &
            + mmr_new(:ncol,pver,so4_ndx)
       call outfld('PM25_SRF',pm25(:ncol) , ncol, lchnk )
    endif
    if ( pm25_srf_diag_soa ) then
       pm25(:ncol) = mmr_new(:ncol,pver,cb1_ndx)   &
            + mmr_new(:ncol,pver,cb2_ndx)   &
            + mmr_new(:ncol,pver,oc1_ndx)   &
            + mmr_new(:ncol,pver,oc2_ndx)   &
            + mmr_new(:ncol,pver,dst1_ndx)  &
            + mmr_new(:ncol,pver,dst2_ndx)  &
            + mmr_new(:ncol,pver,sslt1_ndx) &
            + mmr_new(:ncol,pver,sslt2_ndx) &
            + mmr_new(:ncol,pver,soam_ndx)   &
            + mmr_new(:ncol,pver,soai_ndx)   &
            + mmr_new(:ncol,pver,soat_ndx)   &
            + mmr_new(:ncol,pver,soab_ndx)   &
            + mmr_new(:ncol,pver,soax_ndx)   &
            + mmr_new(:ncol,pver,so4_ndx)
       call outfld('PM25_SRF',pm25(:ncol) , ncol, lchnk )
    endif
!
!
    call outfld('Q_SRF',qh2o(:ncol,pver) , ncol, lchnk )
    call outfld('U_SRF',ufld(:ncol,pver) , ncol, lchnk )
    call outfld('V_SRF',vfld(:ncol,pver) , ncol, lchnk )
!
    if (.not.sad_pbf_ndx>0) then
       deallocate(strato_sad)
    endif

  end subroutine gas_phase_chemdr

  subroutine gas_phase_chemdr_finalize_tendencies(ncol, delt_inverse, mmr, mmr_tend, mmr_new, qtend)

    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), intent(in) :: delt_inverse
    real(r8), target, intent(in) :: mmr(pcols,pver,gas_pcnst)
    real(r8), target, intent(inout) :: mmr_tend(pcols,pver,gas_pcnst)
    real(r8), target, intent(out) :: mmr_new(pcols,pver,gas_pcnst)
    real(r8), target, intent(inout) :: qtend(pcols,pver,pcnst)

    integer :: m, n
    integer(c_int64_t), target :: map2chm_c(pcnst)

    interface
       subroutine gas_phase_chemdr_finalize_tendencies_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, pcnst_c, &
            delt_inverse_c, map2chm_p, mmr_p, mmr_tend_p, mmr_new_p, qtend_p) &
            bind(c, name="gas_phase_chemdr_finalize_tendencies_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c, pcnst_c
         real(c_double), value :: delt_inverse_c
         type(c_ptr), value :: map2chm_p, mmr_p, mmr_tend_p, mmr_new_p, qtend_p
       end subroutine gas_phase_chemdr_finalize_tendencies_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do m = 1, gas_pcnst
          mmr_new(:ncol,:,m) = mmr_tend(:ncol,:,m)
          mmr_tend(:ncol,:,m) = (mmr_tend(:ncol,:,m) - mmr(:ncol,:,m))*delt_inverse
       end do

       do m = 1, pcnst
          n = map2chm(m)
          if (n > 0) then
             qtend(:ncol,:,m) = qtend(:ncol,:,m) + mmr_tend(:ncol,:,n)
          end if
       end do
       return
    end if

    map2chm_c(:) = int(map2chm(:), c_int64_t)

    call gas_phase_chemdr_finalize_tendencies_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         int(pcnst, c_int64_t), real(delt_inverse, c_double), c_loc(map2chm_c), c_loc(mmr), c_loc(mmr_tend), &
         c_loc(mmr_new), c_loc(qtend) &
    )

  end subroutine gas_phase_chemdr_finalize_tendencies

  subroutine gas_phase_chemdr_prepare_sza(ncol, rad2deg_in, zen_angle, sza)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr, c_double

    integer, intent(in) :: ncol
    real(r8), intent(in) :: rad2deg_in
    real(r8), target, intent(inout) :: zen_angle(ncol)
    real(r8), target, intent(out) :: sza(ncol)

    integer :: i

    interface
       subroutine gas_phase_chemdr_prepare_sza_codon(ncol_c, rad2deg_c, zen_angle_p, sza_p) &
            bind(c, name="gas_phase_chemdr_prepare_sza_codon")
         use iso_c_binding, only : c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c
         real(c_double), value :: rad2deg_c
         type(c_ptr), value :: zen_angle_p, sza_p
       end subroutine gas_phase_chemdr_prepare_sza_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
      do i = 1, ncol
         zen_angle(i) = acos( zen_angle(i) )
         sza(i) = zen_angle(i) * rad2deg_in
      end do
      return
    end if

    call gas_phase_chemdr_prepare_sza_codon( &
         int(ncol, c_int64_t), rad2deg_in, c_loc(zen_angle), c_loc(sza) &
    )

  end subroutine gas_phase_chemdr_prepare_sza

  subroutine gas_phase_chemdr_prepare_state(ncol, phis, zi, zm, pmid, zsurf, zintr, zmidr, zmid, zint, pmb)

    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    use physconst, only : rga

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: phis(pcols)
    real(r8), target, intent(in) :: zi(pcols,pver+1)
    real(r8), target, intent(in) :: zm(pcols,pver)
    real(r8), target, intent(in) :: pmid(pcols,pver)
    real(r8), target, intent(out) :: zsurf(ncol)
    real(r8), target, intent(out) :: zintr(ncol,pver+1)
    real(r8), target, intent(out) :: zmidr(ncol,pver)
    real(r8), target, intent(out) :: zmid(ncol,pver)
    real(r8), target, intent(out) :: zint(ncol,pver+1)
    real(r8), target, intent(out) :: pmb(ncol,pver)
    real(r8), parameter :: m2km  = 1.e-3_r8
    real(r8), parameter :: Pa2mb = 1.e-2_r8

    integer :: k

    interface
       subroutine gas_phase_chemdr_prepare_state_codon(ncol_c, pcols_c, pver_c, rga_c, m2km_c, pa2mb_c, phis_p, &
            zi_p, zm_p, pmid_p, zsurf_p, zintr_p, zmidr_p, zmid_p, zint_p, pmb_p) &
            bind(c, name="gas_phase_chemdr_prepare_state_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: rga_c, m2km_c, pa2mb_c
         type(c_ptr), value :: phis_p, zi_p, zm_p, pmid_p, zsurf_p, zintr_p, zmidr_p, zmid_p, zint_p, pmb_p
       end subroutine gas_phase_chemdr_prepare_state_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       zsurf(:ncol) = rga * phis(:ncol)
       do k = 1,pver
          zintr(:ncol,k) = m2km * zi(:ncol,k)
          zmidr(:ncol,k) = m2km * zm(:ncol,k)
          zmid(:ncol,k) = m2km * (zm(:ncol,k) + zsurf(:ncol))
          zint(:ncol,k) = m2km * (zi(:ncol,k) + zsurf(:ncol))
          pmb(:ncol,k)  = Pa2mb * pmid(:ncol,k)
       end do
       zint(:ncol,pver+1) = m2km * (zi(:ncol,pver+1) + zsurf(:ncol))
       zintr(:ncol,pver+1)= m2km *  zi(:ncol,pver+1)
       return
    end if

    call gas_phase_chemdr_prepare_state_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(rga, c_double), real(m2km, c_double), &
         real(Pa2mb, c_double), c_loc(phis), c_loc(zi), c_loc(zm), c_loc(pmid), c_loc(zsurf), c_loc(zintr), &
         c_loc(zmidr), c_loc(zmid), c_loc(zint), c_loc(pmb) &
    )

  end subroutine gas_phase_chemdr_prepare_state

  subroutine gas_phase_chemdr_load_mmr(ncol, q, mmr)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: q(pcols,pver,pcnst)
    real(r8), target, intent(inout) :: mmr(pcols,pver,gas_pcnst)

    integer :: m, n
    integer(c_int64_t), target :: map2chm_c(pcnst)

    interface
       subroutine gas_phase_chemdr_load_mmr_codon(ncol_c, pcols_c, pver_c, pcnst_c, map2chm_p, q_p, mmr_p) &
            bind(c, name="gas_phase_chemdr_load_mmr_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c
         type(c_ptr), value :: map2chm_p, q_p, mmr_p
       end subroutine gas_phase_chemdr_load_mmr_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do m = 1,pcnst
          n = map2chm(m)
          if( n > 0 ) then
             mmr(:ncol,:,n) = q(:ncol,:,m)
          end if
       end do
       return
    end if

    map2chm_c(:) = int(map2chm(:), c_int64_t)

    call gas_phase_chemdr_load_mmr_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnst, c_int64_t), &
         c_loc(map2chm_c), c_loc(q), c_loc(mmr) &
    )

  end subroutine gas_phase_chemdr_load_mmr

  subroutine gas_phase_chemdr_zero_sulfate(ncol, sulfate)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(inout) :: sulfate(ncol,pver)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_zero_sulfate_codon(ncol_c, pver_c, sulfate_p) &
            bind(c, name="gas_phase_chemdr_zero_sulfate_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c
         type(c_ptr), value :: sulfate_p
       end subroutine gas_phase_chemdr_zero_sulfate_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1, pver
          do i = 1, ncol
             sulfate(i,k) = 0._r8
          end do
       end do
       return
    end if

    call gas_phase_chemdr_zero_sulfate_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), c_loc(sulfate) &
    )

  end subroutine gas_phase_chemdr_zero_sulfate

  subroutine gas_phase_chemdr_clip_sulfate(ncol, troplev, sulfate)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, target, intent(in) :: troplev(pcols)
    real(r8), target, intent(inout) :: sulfate(ncol,pver)

    integer :: i, k
    integer(c_int64_t), target :: troplev_c(pcols)

    interface
       subroutine gas_phase_chemdr_clip_sulfate_codon(ncol_c, pcols_c, pver_c, troplev_p, sulfate_p) &
            bind(c, name="gas_phase_chemdr_clip_sulfate_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: troplev_p, sulfate_p
       end subroutine gas_phase_chemdr_clip_sulfate_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1, pver
          do i = 1, ncol
             if( k < troplev(i) ) then
                sulfate(i,k) = 0.0_r8
             end if
          end do
       end do
       return
    end if

    troplev_c(:) = int(troplev(:), c_int64_t)

    call gas_phase_chemdr_clip_sulfate_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(troplev_c), c_loc(sulfate) &
    )

  end subroutine gas_phase_chemdr_clip_sulfate

  subroutine gas_phase_chemdr_load_oxygen_mmr(ncol, o2_ndx_in, o_ndx_in, mmr, o2mmr, ommr)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: o2_ndx_in, o_ndx_in
    real(r8), target, intent(in) :: mmr(pcols,pver,gas_pcnst)
    real(r8), target, intent(out) :: o2mmr(ncol,pver)
    real(r8), target, intent(out) :: ommr(ncol,pver)

    integer :: k

    interface
       subroutine gas_phase_chemdr_load_oxygen_mmr_codon(ncol_c, pcols_c, pver_c, o2_ndx_c, o_ndx_c, mmr_p, &
            o2mmr_p, ommr_p) bind(c, name="gas_phase_chemdr_load_oxygen_mmr_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, o2_ndx_c, o_ndx_c
         type(c_ptr), value :: mmr_p, o2mmr_p, ommr_p
       end subroutine gas_phase_chemdr_load_oxygen_mmr_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          o2mmr(:ncol,k) = mmr(:ncol,k,o2_ndx_in)
          ommr(:ncol,k)  = mmr(:ncol,k,o_ndx_in)
       end do
       return
    end if

    call gas_phase_chemdr_load_oxygen_mmr_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(o2_ndx_in, c_int64_t), &
         int(o_ndx_in, c_int64_t), c_loc(mmr), c_loc(o2mmr), c_loc(ommr) &
    )

  end subroutine gas_phase_chemdr_load_oxygen_mmr

  subroutine gas_phase_chemdr_set_ltrop_sol(ncol, has_linoz_data_flag_in, troplev, ltrop_sol)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: has_linoz_data_flag_in
    integer, target, intent(in) :: troplev(pcols)
    integer, target, intent(out) :: ltrop_sol(pcols)

    interface
       subroutine gas_phase_chemdr_set_ltrop_sol_codon(ncol_c, has_linoz_data_flag_c, troplev_p, ltrop_sol_p) &
            bind(c, name="gas_phase_chemdr_set_ltrop_sol_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, has_linoz_data_flag_c
         type(c_ptr), value :: troplev_p, ltrop_sol_p
       end subroutine gas_phase_chemdr_set_ltrop_sol_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       if (has_linoz_data_flag_in /= 0) then
          ltrop_sol(:ncol) = troplev(:ncol)
       else
          ltrop_sol(:ncol) = 0
       end if
       return
    end if

    call gas_phase_chemdr_set_ltrop_sol_codon( &
         int(ncol, c_int64_t), int(has_linoz_data_flag_in, c_int64_t), c_loc(troplev), c_loc(ltrop_sol) &
    )

  end subroutine gas_phase_chemdr_set_ltrop_sol

  subroutine gas_phase_chemdr_zero_st80_tau(ncol, rxntot_in, st80_25_tau_ndx_in, troplev, reaction_rates)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: rxntot_in, st80_25_tau_ndx_in
    integer, target, intent(in) :: troplev(pcols)
    real(r8), target, intent(inout) :: reaction_rates(ncol,pver,max(1,rxntot_in))

    integer :: i

    interface
       subroutine gas_phase_chemdr_zero_st80_tau_codon(ncol_c, pver_c, rxntot_c, st80_25_tau_ndx_c, troplev_p, &
            reaction_rates_p) bind(c, name="gas_phase_chemdr_zero_st80_tau_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, rxntot_c, st80_25_tau_ndx_c
         type(c_ptr), value :: troplev_p, reaction_rates_p
       end subroutine gas_phase_chemdr_zero_st80_tau_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       if (st80_25_tau_ndx_in > 0) then
          do i = 1,ncol
             reaction_rates(i,1:troplev(i),st80_25_tau_ndx_in) = 0._r8
          enddo
       end if
       return
    end if

    call gas_phase_chemdr_zero_st80_tau_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(rxntot_in, c_int64_t), int(st80_25_tau_ndx_in, c_int64_t), &
         c_loc(troplev), c_loc(reaction_rates) &
    )

  end subroutine gas_phase_chemdr_zero_st80_tau

  subroutine gas_phase_chemdr_compute_relhum(ncol, h2ovmr, satq, relhum)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: h2ovmr(ncol,pver)
    real(r8), target, intent(in) :: satq(ncol,pver)
    real(r8), target, intent(out) :: relhum(ncol,pver)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_compute_relhum_codon(ncol_c, pver_c, h2ovmr_p, satq_p, relhum_p) &
            bind(c, name="gas_phase_chemdr_compute_relhum_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c
         type(c_ptr), value :: h2ovmr_p, satq_p, relhum_p
       end subroutine gas_phase_chemdr_compute_relhum_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             relhum(i,k) = .622_r8 * h2ovmr(i,k) / satq(i,k)
             relhum(i,k) = max(0._r8, min(1._r8, relhum(i,k)))
          end do
       end do
       return
    end if

    call gas_phase_chemdr_compute_relhum_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), c_loc(h2ovmr), c_loc(satq), c_loc(relhum) &
    )

  end subroutine gas_phase_chemdr_compute_relhum

  subroutine gas_phase_chemdr_restore_strat_gases(ncol, hno3_ndx_in, h2o_ndx_in, delt_inverse_in, vmr, hno3_gas, &
       h2o_gas, h2ovmr, wrk)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: hno3_ndx_in, h2o_ndx_in
    real(r8), intent(in) :: delt_inverse_in
    real(r8), target, intent(inout) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: hno3_gas(ncol,pver), h2o_gas(ncol,pver)
    real(r8), target, intent(inout) :: h2ovmr(ncol,pver), wrk(ncol,pver)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_restore_strat_gases_codon(ncol_c, pver_c, gas_pcnst_c, hno3_ndx_c, h2o_ndx_c, &
            delt_inverse_c, vmr_p, hno3_gas_p, h2o_gas_p, h2ovmr_p, wrk_p) &
            bind(c, name="gas_phase_chemdr_restore_strat_gases_codon")
         use iso_c_binding, only : c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, hno3_ndx_c, h2o_ndx_c
         real(c_double), value :: delt_inverse_c
         type(c_ptr), value :: vmr_p, hno3_gas_p, h2o_gas_p, h2ovmr_p, wrk_p
       end subroutine gas_phase_chemdr_restore_strat_gases_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             vmr(i,k,hno3_ndx_in) = hno3_gas(i,k)
             h2ovmr(i,k) = h2o_gas(i,k)
             vmr(i,k,h2o_ndx_in) = h2o_gas(i,k)
             wrk(i,k) = (h2ovmr(i,k) - wrk(i,k))*delt_inverse_in
          end do
       end do
       return
    end if

    call gas_phase_chemdr_restore_strat_gases_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(hno3_ndx_in, c_int64_t), &
         int(h2o_ndx_in, c_int64_t), delt_inverse_in, c_loc(vmr), c_loc(hno3_gas), c_loc(h2o_gas), &
         c_loc(h2ovmr), c_loc(wrk) &
    )

  end subroutine gas_phase_chemdr_restore_strat_gases

  subroutine gas_phase_chemdr_restore_hcl_gas(ncol, hcl_ndx_in, vmr, hcl_gas)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: hcl_ndx_in
    real(r8), target, intent(inout) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: hcl_gas(ncol,pver)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_restore_hcl_gas_codon(ncol_c, pver_c, gas_pcnst_c, hcl_ndx_c, vmr_p, hcl_gas_p) &
            bind(c, name="gas_phase_chemdr_restore_hcl_gas_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, hcl_ndx_c
         type(c_ptr), value :: vmr_p, hcl_gas_p
       end subroutine gas_phase_chemdr_restore_hcl_gas_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             vmr(i,k,hcl_ndx_in) = hcl_gas(i,k)
          end do
       end do
       return
    end if

    call gas_phase_chemdr_restore_hcl_gas_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(hcl_ndx_in, c_int64_t), &
         c_loc(vmr), c_loc(hcl_gas) &
    )

  end subroutine gas_phase_chemdr_restore_hcl_gas

  subroutine gas_phase_chemdr_init_dust_vmr(ncol, dst_ndx_in, vmr, dust_vmr)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: dst_ndx_in
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(out) :: dust_vmr(ncol,pver,ndust)

    integer :: i, k, m

    interface
       subroutine gas_phase_chemdr_init_dust_vmr_codon(ncol_c, pver_c, gas_pcnst_c, ndust_c, dst_ndx_c, vmr_p, &
            dust_vmr_p) bind(c, name="gas_phase_chemdr_init_dust_vmr_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, ndust_c, dst_ndx_c
         type(c_ptr), value :: vmr_p, dust_vmr_p
       end subroutine gas_phase_chemdr_init_dust_vmr_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       if (dst_ndx_in > 0) then
          do m = 1, ndust
             do k = 1, pver
                do i = 1, ncol
                   dust_vmr(i,k,m) = vmr(i,k,dst_ndx_in + m - 1)
                end do
             end do
          end do
       else
          do m = 1, ndust
             do k = 1, pver
                do i = 1, ncol
                   dust_vmr(i,k,m) = 0.0_r8
                end do
             end do
          end do
       end if
       return
    end if

    call gas_phase_chemdr_init_dust_vmr_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(ndust, c_int64_t), &
         int(dst_ndx_in, c_int64_t), c_loc(vmr), c_loc(dust_vmr) &
    )

  end subroutine gas_phase_chemdr_init_dust_vmr

  subroutine gas_phase_chemdr_reset_ste_tracer(ncol, st80_25_ndx_in, pmid_threshold_in, st80_vmr_in, pmid, vmr)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: st80_25_ndx_in
    real(r8), intent(in) :: pmid_threshold_in
    real(r8), intent(in) :: st80_vmr_in
    real(r8), target, intent(in) :: pmid(pcols,pver)
    real(r8), target, intent(inout) :: vmr(ncol,pver,gas_pcnst)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_reset_ste_tracer_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, st80_25_ndx_c, &
            pmid_threshold_c, st80_vmr_c, pmid_p, vmr_p) bind(c, name="gas_phase_chemdr_reset_ste_tracer_codon")
         use iso_c_binding, only : c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c, st80_25_ndx_c
         real(c_double), value :: pmid_threshold_c, st80_vmr_c
         type(c_ptr), value :: pmid_p, vmr_p
       end subroutine gas_phase_chemdr_reset_ste_tracer_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             if (pmid(i,k) < pmid_threshold_in) then
                vmr(i,k,st80_25_ndx_in) = st80_vmr_in
             end if
          end do
       end do
       return
    end if

    call gas_phase_chemdr_reset_ste_tracer_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         int(st80_25_ndx_in, c_int64_t), pmid_threshold_in, st80_vmr_in, c_loc(pmid), c_loc(vmr) &
    )

  end subroutine gas_phase_chemdr_reset_ste_tracer

  subroutine gas_phase_chemdr_zero_sflx(sflx)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    real(r8), target, intent(out) :: sflx(pcols,gas_pcnst)

    integer :: i, m

    interface
       subroutine gas_phase_chemdr_zero_sflx_codon(pcols_c, gas_pcnst_c, sflx_p) &
            bind(c, name="gas_phase_chemdr_zero_sflx_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: pcols_c, gas_pcnst_c
         type(c_ptr), value :: sflx_p
       end subroutine gas_phase_chemdr_zero_sflx_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do m = 1, gas_pcnst
          do i = 1, pcols
             sflx(i,m) = 0._r8
          end do
       end do
       return
    end if

    call gas_phase_chemdr_zero_sflx_codon( &
         int(pcols, c_int64_t), int(gas_pcnst, c_int64_t), c_loc(sflx) &
    )

  end subroutine gas_phase_chemdr_zero_sflx

  subroutine gas_phase_chemdr_compute_wind_speed(ncol, ufld, vfld, wind_speed)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: ufld(pcols,pver)
    real(r8), target, intent(in) :: vfld(pcols,pver)
    real(r8), target, intent(inout) :: wind_speed(pcols)

    integer :: i

    interface
       subroutine gas_phase_chemdr_compute_wind_speed_codon(ncol_c, pcols_c, pver_c, ufld_p, vfld_p, wind_speed_p) &
            bind(c, name="gas_phase_chemdr_compute_wind_speed_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: ufld_p, vfld_p, wind_speed_p
       end subroutine gas_phase_chemdr_compute_wind_speed_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do i = 1, ncol
          wind_speed(i) = sqrt( ufld(i,pver)*ufld(i,pver) + vfld(i,pver)*vfld(i,pver) )
       end do
       return
    end if

    call gas_phase_chemdr_compute_wind_speed_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(ufld), c_loc(vfld), c_loc(wind_speed) &
    )

  end subroutine gas_phase_chemdr_compute_wind_speed

  subroutine gas_phase_chemdr_compute_prect(ncol, precc, precl, prect)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: precc(pcols)
    real(r8), target, intent(in) :: precl(pcols)
    real(r8), target, intent(out) :: prect(pcols)

    integer :: i

    interface
       subroutine gas_phase_chemdr_compute_prect_codon(ncol_c, pcols_c, precc_p, precl_p, prect_p) &
            bind(c, name="gas_phase_chemdr_compute_prect_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c
         type(c_ptr), value :: precc_p, precl_p, prect_p
       end subroutine gas_phase_chemdr_compute_prect_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do i = 1,ncol
          prect(i) = precc(i) + precl(i)
       end do
       return
    end if

    call gas_phase_chemdr_compute_prect_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), c_loc(precc), c_loc(precl), c_loc(prect) &
    )

  end subroutine gas_phase_chemdr_compute_prect

  subroutine gas_phase_chemdr_compute_tvs(ncol, tfld, qh2o, tvs)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: tfld(pcols,pver)
    real(r8), target, intent(in) :: qh2o(pcols,pver)
    real(r8), target, intent(out) :: tvs(pcols)

    integer :: i

    interface
       subroutine gas_phase_chemdr_compute_tvs_codon(ncol_c, pcols_c, pver_c, tfld_p, qh2o_p, tvs_p) &
            bind(c, name="gas_phase_chemdr_compute_tvs_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: tfld_p, qh2o_p, tvs_p
       end subroutine gas_phase_chemdr_compute_tvs_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do i = 1,ncol
          tvs(i) = tfld(i,pver) * (1._r8 + qh2o(i,pver))
       end do
       return
    end if

    call gas_phase_chemdr_compute_tvs_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(tfld), c_loc(qh2o), c_loc(tvs) &
    )

  end subroutine gas_phase_chemdr_compute_tvs

  subroutine gas_phase_chemdr_copy_cldw_to_cwat(ncol, cldw, cwat)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: cldw(pcols,pver)
    real(r8), target, intent(out) :: cwat(ncol,pver)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_copy_cldw_to_cwat_codon(ncol_c, pcols_c, pver_c, cldw_p, cwat_p) &
            bind(c, name="gas_phase_chemdr_copy_cldw_to_cwat_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         type(c_ptr), value :: cldw_p, cwat_p
       end subroutine gas_phase_chemdr_copy_cldw_to_cwat_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             cwat(i,k) = cldw(i,k)
          end do
       end do
       return
    end if

    call gas_phase_chemdr_copy_cldw_to_cwat_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), c_loc(cldw), c_loc(cwat) &
    )

  end subroutine gas_phase_chemdr_copy_cldw_to_cwat

  subroutine gas_phase_chemdr_load_h2o_fields(ncol, h2o_ndx_in, mmr, vmr, qh2o, h2ovmr)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: h2o_ndx_in
    real(r8), target, intent(in) :: mmr(pcols,pver,gas_pcnst)
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(out) :: qh2o(pcols,pver)
    real(r8), target, intent(out) :: h2ovmr(ncol,pver)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_load_h2o_fields_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, h2o_ndx_c, mmr_p, &
            vmr_p, qh2o_p, h2ovmr_p) bind(c, name="gas_phase_chemdr_load_h2o_fields_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c, h2o_ndx_c
         type(c_ptr), value :: mmr_p, vmr_p, qh2o_p, h2ovmr_p
       end subroutine gas_phase_chemdr_load_h2o_fields_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             qh2o(i,k) = mmr(i,k,h2o_ndx_in)
             h2ovmr(i,k) = vmr(i,k,h2o_ndx_in)
          end do
       end do
       return
    end if

    call gas_phase_chemdr_load_h2o_fields_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         int(h2o_ndx_in, c_int64_t), c_loc(mmr), c_loc(vmr), c_loc(qh2o), c_loc(h2ovmr) &
    )

  end subroutine gas_phase_chemdr_load_h2o_fields

  subroutine gas_phase_chemdr_copy_o3_to_o3s_trop(ncol, troplev, o3_ndx_in, o3s_ndx_in, vmr)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, target, intent(in) :: troplev(pcols)
    integer, intent(in) :: o3_ndx_in, o3s_ndx_in
    real(r8), target, intent(inout) :: vmr(ncol,pver,gas_pcnst)

    integer :: i, k
    integer(c_int64_t), target :: troplev_c(pcols)

    interface
       subroutine gas_phase_chemdr_copy_o3_to_o3s_trop_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, troplev_p, &
            o3_ndx_c, o3s_ndx_c, vmr_p) bind(c, name="gas_phase_chemdr_copy_o3_to_o3s_trop_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c, o3_ndx_c, o3s_ndx_c
         type(c_ptr), value :: troplev_p, vmr_p
       end subroutine gas_phase_chemdr_copy_o3_to_o3s_trop_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do i = 1,ncol
          do k = 1,troplev(i)
             vmr(i,k,o3s_ndx_in) = vmr(i,k,o3_ndx_in)
          end do
       end do
       return
    end if

    troplev_c(:) = int(troplev(:), c_int64_t)

    call gas_phase_chemdr_copy_o3_to_o3s_trop_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         c_loc(troplev_c), int(o3_ndx_in, c_int64_t), int(o3s_ndx_in, c_int64_t), c_loc(vmr) &
    )

  end subroutine gas_phase_chemdr_copy_o3_to_o3s_trop

  subroutine gas_phase_chemdr_copy_h2o_to_wrk(ncol, h2o_ndx_in, vmr, wrk)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: h2o_ndx_in
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(out) :: wrk(ncol,pver)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_copy_h2o_to_wrk_codon(ncol_c, pver_c, gas_pcnst_c, h2o_ndx_c, vmr_p, wrk_p) &
            bind(c, name="gas_phase_chemdr_copy_h2o_to_wrk_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, h2o_ndx_c
         type(c_ptr), value :: vmr_p, wrk_p
       end subroutine gas_phase_chemdr_copy_h2o_to_wrk_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             wrk(i,k) = vmr(i,k,h2o_ndx_in)
          end do
       end do
       return
    end if

    call gas_phase_chemdr_copy_h2o_to_wrk_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(h2o_ndx_in, c_int64_t), &
         c_loc(vmr), c_loc(wrk) &
    )

  end subroutine gas_phase_chemdr_copy_h2o_to_wrk

  subroutine gas_phase_chemdr_update_qdsett_wrk(ncol, h2o_ndx_in, delt_inverse_in, vmr, wrk)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: h2o_ndx_in
    real(r8), intent(in) :: delt_inverse_in
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: wrk(ncol,pver)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_update_qdsett_wrk_codon(ncol_c, pver_c, gas_pcnst_c, h2o_ndx_c, delt_inverse_c, &
            vmr_p, wrk_p) bind(c, name="gas_phase_chemdr_update_qdsett_wrk_codon")
         use iso_c_binding, only : c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, h2o_ndx_c
         real(c_double), value :: delt_inverse_c
         type(c_ptr), value :: vmr_p, wrk_p
       end subroutine gas_phase_chemdr_update_qdsett_wrk_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             wrk(i,k) = (vmr(i,k,h2o_ndx_in) - wrk(i,k))*delt_inverse_in
          end do
       end do
       return
    end if

    call gas_phase_chemdr_update_qdsett_wrk_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(h2o_ndx_in, c_int64_t), &
         delt_inverse_in, c_loc(vmr), c_loc(wrk) &
    )

  end subroutine gas_phase_chemdr_update_qdsett_wrk

  subroutine gas_phase_chemdr_update_qdchem_wrk(ncol, h2o_ndx_in, delt_inverse_in, vmr, wrk)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: h2o_ndx_in
    real(r8), intent(in) :: delt_inverse_in
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: wrk(ncol,pver)

    integer :: i, k

    interface
       subroutine gas_phase_chemdr_update_qdchem_wrk_codon(ncol_c, pver_c, gas_pcnst_c, h2o_ndx_c, delt_inverse_c, &
            vmr_p, wrk_p) bind(c, name="gas_phase_chemdr_update_qdchem_wrk_codon")
         use iso_c_binding, only : c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, h2o_ndx_c
         real(c_double), value :: delt_inverse_c
         type(c_ptr), value :: vmr_p, wrk_p
       end subroutine gas_phase_chemdr_update_qdchem_wrk_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             wrk(i,k) = (vmr(i,k,h2o_ndx_in) - wrk(i,k))*delt_inverse_in
          end do
       end do
       return
    end if

    call gas_phase_chemdr_update_qdchem_wrk_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(h2o_ndx_in, c_int64_t), &
         delt_inverse_in, c_loc(vmr), c_loc(wrk) &
    )

  end subroutine gas_phase_chemdr_update_qdchem_wrk

  subroutine gas_phase_chemdr_init_stratchem_state(ncol, hno3_ndx_in, hcl_ndx_in, cldice_ndx_in, vmr, h2ovmr, q, &
       hcl_cond, hcl_gas, hno3_gas, h2o_gas, wrk, cldice, hno3_cond)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: hno3_ndx_in, hcl_ndx_in, cldice_ndx_in
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst), h2ovmr(ncol,pver), q(pcols,pver,pcnst)
    real(r8), target, intent(out) :: hcl_cond(ncol,pver), hcl_gas(ncol,pver), hno3_gas(ncol,pver)
    real(r8), target, intent(out) :: h2o_gas(ncol,pver), wrk(ncol,pver), cldice(pcols,pver)
    real(r8), target, intent(out) :: hno3_cond(ncol,pver,2)

    integer :: i, k, m

    interface
       subroutine gas_phase_chemdr_init_stratchem_state_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, pcnst_c, &
            hno3_ndx_c, hcl_ndx_c, cldice_ndx_c, vmr_p, h2ovmr_p, q_p, hcl_cond_p, hcl_gas_p, hno3_gas_p, &
            h2o_gas_p, wrk_p, cldice_p, hno3_cond_p) bind(c, name="gas_phase_chemdr_init_stratchem_state_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c, pcnst_c
         integer(c_int64_t), value :: hno3_ndx_c, hcl_ndx_c, cldice_ndx_c
         type(c_ptr), value :: vmr_p, h2ovmr_p, q_p, hcl_cond_p, hcl_gas_p, hno3_gas_p, h2o_gas_p
         type(c_ptr), value :: wrk_p, cldice_p, hno3_cond_p
       end subroutine gas_phase_chemdr_init_stratchem_state_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          do i = 1,ncol
             hcl_cond(i,k) = 0.0_r8
             hno3_cond(i,k,1) = 0.0_r8
             hno3_cond(i,k,2) = 0.0_r8
             hno3_gas(i,k) = vmr(i,k,hno3_ndx_in)
             h2o_gas(i,k) = h2ovmr(i,k)
             hcl_gas(i,k) = vmr(i,k,hcl_ndx_in)
             wrk(i,k) = h2ovmr(i,k)
             cldice(i,k) = q(i,k,cldice_ndx_in)
          end do
       end do
       return
    end if

    call gas_phase_chemdr_init_stratchem_state_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         int(pcnst, c_int64_t), int(hno3_ndx_in, c_int64_t), int(hcl_ndx_in, c_int64_t), &
         int(cldice_ndx_in, c_int64_t), c_loc(vmr), c_loc(h2ovmr), c_loc(q), c_loc(hcl_cond), c_loc(hcl_gas), &
         c_loc(hno3_gas), c_loc(h2o_gas), c_loc(wrk), c_loc(cldice), c_loc(hno3_cond) &
    )

  end subroutine gas_phase_chemdr_init_stratchem_state

  subroutine gas_phase_chemdr_init_h2so4_gasprod(ncol, ndx_h2so4_in, vmr, del_h2so4_gasprod)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: ndx_h2so4_in
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(out) :: del_h2so4_gasprod(ncol,pver)

    interface
       subroutine gas_phase_chemdr_init_h2so4_gasprod_codon(ncol_c, pver_c, gas_pcnst_c, ndx_h2so4_c, vmr_p, &
            del_h2so4_gasprod_p) bind(c, name="gas_phase_chemdr_init_h2so4_gasprod_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, ndx_h2so4_c
         type(c_ptr), value :: vmr_p, del_h2so4_gasprod_p
       end subroutine gas_phase_chemdr_init_h2so4_gasprod_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       if (ndx_h2so4_in > 0) then
          del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4_in)
       else
          del_h2so4_gasprod(:,:) = 0.0_r8
       end if
       return
    end if

    call gas_phase_chemdr_init_h2so4_gasprod_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(ndx_h2so4_in, c_int64_t), &
         c_loc(vmr), c_loc(del_h2so4_gasprod) &
    )

  end subroutine gas_phase_chemdr_init_h2so4_gasprod

  subroutine gas_phase_chemdr_store_vmr0(ncol, vmr, vmr0)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(out) :: vmr0(ncol,pver,gas_pcnst)

    interface
       subroutine gas_phase_chemdr_store_vmr0_codon(ncol_c, pver_c, gas_pcnst_c, vmr_p, vmr0_p) &
            bind(c, name="gas_phase_chemdr_store_vmr0_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c
         type(c_ptr), value :: vmr_p, vmr0_p
       end subroutine gas_phase_chemdr_store_vmr0_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       vmr0(:ncol,:,:) = vmr(:ncol,:,:)
       return
    end if

    call gas_phase_chemdr_store_vmr0_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), c_loc(vmr), c_loc(vmr0) &
    )

  end subroutine gas_phase_chemdr_store_vmr0

  subroutine gas_phase_chemdr_update_h2so4_gasprod(ncol, ndx_h2so4_in, vmr, del_h2so4_gasprod)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: ndx_h2so4_in
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: del_h2so4_gasprod(ncol,pver)

    interface
       subroutine gas_phase_chemdr_update_h2so4_gasprod_codon(ncol_c, pver_c, gas_pcnst_c, ndx_h2so4_c, vmr_p, &
            del_h2so4_gasprod_p) bind(c, name="gas_phase_chemdr_update_h2so4_gasprod_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, ndx_h2so4_c
         type(c_ptr), value :: vmr_p, del_h2so4_gasprod_p
       end subroutine gas_phase_chemdr_update_h2so4_gasprod_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4_in) - del_h2so4_gasprod(1:ncol,:)
       return
    end if

    call gas_phase_chemdr_update_h2so4_gasprod_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), int(ndx_h2so4_in, c_int64_t), &
         c_loc(vmr), c_loc(del_h2so4_gasprod) &
    )

  end subroutine gas_phase_chemdr_update_h2so4_gasprod

  subroutine gas_phase_chemdr_reform_hno3_hcl(ncol, hno3_ndx_in, hcl_ndx_in, vmr, hno3_cond, hcl_cond)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: hno3_ndx_in, hcl_ndx_in
    real(r8), target, intent(inout) :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: hno3_cond(ncol,pver,2)
    real(r8), target, intent(in) :: hcl_cond(ncol,pver)

    integer :: k

    interface
       subroutine gas_phase_chemdr_reform_hno3_hcl_codon(ncol_c, pver_c, gas_pcnst_c, hno3_ndx_c, hcl_ndx_c, &
            vmr_p, hno3_cond_p, hcl_cond_p) bind(c, name="gas_phase_chemdr_reform_hno3_hcl_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, hno3_ndx_c, hcl_ndx_c
         type(c_ptr), value :: vmr_p, hno3_cond_p, hcl_cond_p
       end subroutine gas_phase_chemdr_reform_hno3_hcl_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do k = 1,pver
          vmr(:,k,hno3_ndx_in) = vmr(:,k,hno3_ndx_in) + hno3_cond(:,k,1) + hno3_cond(:,k,2)
          vmr(:,k,hcl_ndx_in) = vmr(:,k,hcl_ndx_in) + hcl_cond(:,k)
       end do
       return
    end if

    call gas_phase_chemdr_reform_hno3_hcl_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         int(hno3_ndx_in, c_int64_t), int(hcl_ndx_in, c_int64_t), &
         c_loc(vmr), c_loc(hno3_cond), c_loc(hcl_cond) &
    )

  end subroutine gas_phase_chemdr_reform_hno3_hcl

  subroutine gas_phase_chemdr_normalize_extfrc(ncol, extcnt_in, nfs_in, indexm_in, extfrc, invariants)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    integer, intent(in) :: extcnt_in, nfs_in, indexm_in
    real(r8), target, intent(inout) :: extfrc(ncol,pver,max(1,extcnt_in))
    real(r8), target, intent(in) :: invariants(ncol,pver,nfs_in)

    integer :: m, k

    interface
       subroutine gas_phase_chemdr_normalize_extfrc_codon(ncol_c, pver_c, extcnt_c, synoz_ndx_c, aoa_nh_ext_ndx_c, &
            indexm_c, extfrc_p, invariants_p) bind(c, name="gas_phase_chemdr_normalize_extfrc_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, extcnt_c, synoz_ndx_c, aoa_nh_ext_ndx_c, indexm_c
         type(c_ptr), value :: extfrc_p, invariants_p
       end subroutine gas_phase_chemdr_normalize_extfrc_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       do m = 1,extcnt_in
          if( m /= synoz_ndx .and. m /= aoa_nh_ext_ndx ) then
             do k = 1,pver
                extfrc(:ncol,k,m) = extfrc(:ncol,k,m) / invariants(:ncol,k,indexm_in)
             end do
          endif
       end do
       return
    end if

    call gas_phase_chemdr_normalize_extfrc_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(extcnt_in, c_int64_t), int(synoz_ndx, c_int64_t), &
         int(aoa_nh_ext_ndx, c_int64_t), int(indexm_in, c_int64_t), c_loc(extfrc), c_loc(invariants) &
    )

  end subroutine gas_phase_chemdr_normalize_extfrc

  subroutine gas_phase_chemdr_store_drydep(ncol, sflx, cflx, drydepflx)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: sflx(pcols,gas_pcnst)
    real(r8), target, intent(inout) :: cflx(pcols,pcnst)
    real(r8), target, intent(out) :: drydepflx(pcols,pcnst)

    integer :: m, n
    integer(c_int64_t), target :: map2chm_c(pcnst)

    interface
       subroutine gas_phase_chemdr_store_drydep_codon(ncol_c, pcols_c, gas_pcnst_c, pcnst_c, map2chm_p, sflx_p, &
            cflx_p, drydepflx_p) bind(c, name="gas_phase_chemdr_store_drydep_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, gas_pcnst_c, pcnst_c
         type(c_ptr), value :: map2chm_p, sflx_p, cflx_p, drydepflx_p
       end subroutine gas_phase_chemdr_store_drydep_codon
    end interface

    if (gas_phase_chemdr_use_native_impl) then
       drydepflx(:,:) = 0._r8
       do m = 1, pcnst
          n = map2chm(m)
          if (n > 0) then
             cflx(:ncol,m)      = cflx(:ncol,m) - sflx(:ncol,n)
             drydepflx(:ncol,m) = sflx(:ncol,n)
          end if
       end do
       return
    end if

    map2chm_c(:) = int(map2chm(:), c_int64_t)

    call gas_phase_chemdr_store_drydep_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(gas_pcnst, c_int64_t), int(pcnst, c_int64_t), &
         c_loc(map2chm_c), c_loc(sflx), c_loc(cflx), c_loc(drydepflx) &
    )

  end subroutine gas_phase_chemdr_store_drydep

  subroutine gas_phase_chemdr_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (gas_phase_chemdr_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('GAS_PHASE_CHEMDR_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       gas_phase_chemdr_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       gas_phase_chemdr_use_native_impl = .false.
    end if

    gas_phase_chemdr_impl_selected = .true.

    if (masterproc) then
       if (gas_phase_chemdr_use_native_impl) then
          write(iulog,*) 'gas_phase_chemdr helper implementation = native'
       else
          write(iulog,*) 'gas_phase_chemdr helper implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine gas_phase_chemdr_select_impl

end module mo_gas_phase_chemdr
