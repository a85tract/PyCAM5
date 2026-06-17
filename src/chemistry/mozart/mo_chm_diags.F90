module mo_chm_diags

  use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : gas_pcnst
  use mo_tracname,  only : solsym
  use chem_mods,    only : rxntot, nfs, gas_pcnst, indexm, adv_mass
  use ppgrid,       only : pver
  use mo_constants, only : pi, rgrav, rearth
  use mo_chem_utls, only : get_rxt_ndx, get_spc_ndx
  use cam_history,  only : fieldname_len
  use cam_logfile,  only : iulog
  use mo_jeuv,      only : neuv
  use mo_util,      only : chemistry_misc_codon_touch
  use spmd_utils,   only : masterproc

  private

  public :: chm_diags_inti
  public :: chm_diags
  public :: het_diags

  integer :: id_n,id_no,id_no2,id_no3,id_n2o5,id_hno3,id_ho2no2,id_clono2,id_brono2
  integer :: id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_brcl
  integer :: id_ccl4,id_cfc11,id_cfc113,id_ch3ccl3,id_cfc12,id_ch3cl,id_hcfc22,id_cf3br,id_cf2clbr
  integer :: id_cfc114,id_cfc115,id_hcfc141b,id_hcfc142b,id_h1202,id_h2402,id_ch2br2,id_chbr3
  integer :: id_hf,id_f,id_cof2,id_cofcl,id_ch3br
  integer :: id_br,id_bro,id_hbr,id_hobr,id_ch4,id_h2o,id_h2
  integer :: id_o,id_o2,id_h

  integer, parameter :: NJEUV = neuv
  integer :: rid_jeuv(NJEUV), rid_jno_i, rid_jno

  logical :: has_jeuvs, has_jno_i, has_jno

  integer :: nox_species(3),  noy_species(15)
  integer :: clox_species(6), cloy_species(9), tcly_species(21)
  integer :: brox_species(4), broy_species(6), tbry_species(13)
  integer :: foy_species(4),  tfy_species(16)
  integer :: toth_species(3)
  integer :: sox_species(3)
  integer :: nhx_species(3)
  integer :: aer_species(gas_pcnst)

  character(len=fieldname_len) :: dtchem_name(gas_pcnst)
  character(len=fieldname_len) :: depvel_name(gas_pcnst)
  character(len=fieldname_len) :: depflx_name(gas_pcnst)
  character(len=fieldname_len) :: wetdep_name(gas_pcnst)
  character(len=fieldname_len) :: wtrate_name(gas_pcnst)

  real(r8), parameter :: N_molwgt = 14.00674_r8
  real(r8), parameter :: S_molwgt = 32.066_r8

  character(len=32) :: chempkg
  logical :: chm_diags_use_native_impl = .false.
  logical :: chm_diags_impl_selected = .false.
  logical :: chm_diags_proof_written = .false.

contains

  subroutine chm_diags_inti
    !--------------------------------------------------------------------
    !	... initialize utility routine
    !--------------------------------------------------------------------

    use cam_history,  only : addfld, phys_decomp, add_default
    use constituents, only : cnst_get_ind, cnst_longname
    use phys_control, only: phys_getopts

    implicit none

    integer :: i, j, k, m, n
    character(len=16) :: jname, spc_name, attr
    character(len=2)  :: jchar
    character(len=64) :: lname
    character(len=2)  :: unit_basename  ! Units 'kg' or '1' 

    integer :: id_pan, id_onit, id_mpan, id_isopno3, id_onitr, id_nh4no3
    integer :: id_so2, id_so4, id_h2so4
    integer :: id_nh3, id_nh4
    integer :: id_dst01, id_dst02, id_dst03, id_dst04, id_sslt01, id_sslt02, id_sslt03, id_sslt04
    integer :: id_soa,  id_oc1, id_oc2, id_cb1, id_cb2
    integer :: id_soam,id_soai,id_soat,id_soab,id_soax
    integer :: id_bry, id_cly 

    logical :: history_aerosol      ! Output the MAM aerosol tendencies
    logical :: history_amwg         ! output the variables used by the AMWG diag package
    integer :: bulkaero_species(20)

    !-----------------------------------------------------------------------

    call chemistry_misc_codon_touch('mo_chm_diags', 104)

    call phys_getopts( history_aerosol_out = history_aerosol, &
                       history_amwg_out    = history_amwg,  &
                       cam_chempkg_out     = chempkg   )

    id_bry     = get_spc_ndx( 'BRY' )
    id_cly     = get_spc_ndx( 'CLY' )

    id_n       = get_spc_ndx( 'N' )
    id_no      = get_spc_ndx( 'NO' )
    id_no2     = get_spc_ndx( 'NO2' )
    id_no3     = get_spc_ndx( 'NO3' )
    id_n2o5    = get_spc_ndx( 'N2O5' )
    id_hno3    = get_spc_ndx( 'HNO3' )
    id_ho2no2  = get_spc_ndx( 'HO2NO2' )
    id_clono2  = get_spc_ndx( 'CLONO2' )
    id_brono2  = get_spc_ndx( 'BRONO2' )
    id_cl      = get_spc_ndx( 'CL' )
    id_clo     = get_spc_ndx( 'CLO' )
    id_hocl    = get_spc_ndx( 'HOCL' )
    id_cl2     = get_spc_ndx( 'CL2' )
    id_cl2o2   = get_spc_ndx( 'CL2O2' )
    id_oclo    = get_spc_ndx( 'OCLO' )
    id_hcl     = get_spc_ndx( 'HCL' )
    id_brcl    = get_spc_ndx( 'BRCL' )

    id_f       = get_spc_ndx( 'F' )
    id_hf      = get_spc_ndx( 'HF' )
    id_cofcl   = get_spc_ndx( 'COFCL' )
    id_cof2    = get_spc_ndx( 'COF2' )

    id_ccl4    = get_spc_ndx( 'CCL4' )
    id_cfc11   = get_spc_ndx( 'CFC11' )

    id_cfc113  = get_spc_ndx( 'CFC113' )
    id_cfc114  = get_spc_ndx( 'CFC114' )
    id_cfc115  = get_spc_ndx( 'CFC115' )

    id_ch3ccl3 = get_spc_ndx( 'CH3CCL3' )
    id_cfc12   = get_spc_ndx( 'CFC12' )
    id_ch3cl   = get_spc_ndx( 'CH3CL' )

    id_hcfc22  = get_spc_ndx( 'HCFC22' )
    id_hcfc141b= get_spc_ndx( 'HCFC141B' )
    id_hcfc142b= get_spc_ndx( 'HCFC142B' )

    id_cf2clbr = get_spc_ndx( 'CF2CLBR' )
    id_cf3br   = get_spc_ndx( 'CF3BR' )
    id_ch3br   = get_spc_ndx( 'CH3BR' )
    id_h1202   = get_spc_ndx( 'H1202' )
    id_h2402   = get_spc_ndx( 'H2402' )
    id_ch2br2  = get_spc_ndx( 'CH2BR2' )
    id_chbr3   = get_spc_ndx( 'CHBR3' )

    id_br      = get_spc_ndx( 'BR' )
    id_bro     = get_spc_ndx( 'BRO' )
    id_hbr     = get_spc_ndx( 'HBR' )
    id_hobr    = get_spc_ndx( 'HOBR' )
    id_ch4     = get_spc_ndx( 'CH4' )
    id_h2o     = get_spc_ndx( 'H2O' )
    id_h2      = get_spc_ndx( 'H2' )
    id_o       = get_spc_ndx( 'O' )
    id_o2      = get_spc_ndx( 'O2' )
    id_h       = get_spc_ndx( 'H' )

    id_pan     = get_spc_ndx( 'PAN' )
    id_onit    = get_spc_ndx( 'ONIT' )
    id_mpan    = get_spc_ndx( 'MPAN' )
    id_isopno3 = get_spc_ndx( 'ISOPNO3' )
    id_onitr   = get_spc_ndx( 'ONITR' )
    id_nh4no3  = get_spc_ndx( 'NH4NO3' )

    id_so2     = get_spc_ndx( 'SO2' )
    id_so4     = get_spc_ndx( 'SO4' )
    id_h2so4   = get_spc_ndx( 'H2SO4' )

    id_nh3     = get_spc_ndx( 'NH3' )
    id_nh4     = get_spc_ndx( 'NH4' )
    id_nh4no3  = get_spc_ndx( 'NH4NO3' )

    id_dst01   = get_spc_ndx( 'DST01' )
    id_dst02   = get_spc_ndx( 'DST02' )
    id_dst03   = get_spc_ndx( 'DST03' )
    id_dst04   = get_spc_ndx( 'DST04' )
    id_sslt01  = get_spc_ndx( 'SSLT01' )
    id_sslt02  = get_spc_ndx( 'SSLT02' )
    id_sslt03  = get_spc_ndx( 'SSLT03' )
    id_sslt04  = get_spc_ndx( 'SSLT04' )
    id_soa     = get_spc_ndx( 'SOA' )
    id_so4     = get_spc_ndx( 'SO4' )
    id_oc1     = get_spc_ndx( 'OC1' )
    id_oc2     = get_spc_ndx( 'OC2' )
    id_cb1     = get_spc_ndx( 'CB1' )
    id_cb2     = get_spc_ndx( 'CB2' )

    rid_jno   = get_rxt_ndx( 'jno' )
    rid_jno_i = get_rxt_ndx( 'jno_i' )

    id_soam = get_spc_ndx( 'SOAM' )
    id_soai = get_spc_ndx( 'SOAI' )
    id_soat = get_spc_ndx( 'SOAT' )
    id_soab = get_spc_ndx( 'SOAB' )
    id_soax = get_spc_ndx( 'SOAX' )


!... NOY species
    nox_species = (/ id_n, id_no, id_no2 /)
    noy_species = (/ id_n, id_no, id_no2, id_no3, id_n2o5, id_hno3, id_ho2no2, id_clono2, &
                     id_brono2, id_pan, id_onit, id_mpan, id_isopno3, id_onitr, id_nh4no3 /)
!... CLOY species
    clox_species = (/ id_cl, id_clo, id_hocl, id_cl2, id_cl2o2, id_oclo /)
    cloy_species = (/ id_cl, id_clo, id_hocl, id_cl2, id_cl2o2, id_oclo, id_hcl, id_clono2, id_brcl /)
    tcly_species = (/ id_cl, id_clo, id_hocl, id_cl2, id_cl2o2, id_oclo, id_hcl, id_clono2, id_brcl, &
                      id_ccl4, id_cfc11, id_cfc113, id_cfc114, id_cfc115, id_ch3ccl3, id_cfc12, id_ch3cl, &
                      id_hcfc22, id_hcfc141b, id_hcfc142b, id_cf2clbr /)

!... FOY species
    foy_species = (/ id_F, id_hf, id_cofcl, id_cof2 /)
    tfy_species = (/ id_f, id_hf, id_cofcl, id_cof2, id_cfc11, id_cfc12, id_cfc113, id_cfc114, id_cfc115, &
                     id_hcfc22, id_hcfc141b, id_hcfc142b, id_cf2clbr, id_cf3br, id_h1202, id_h2402 /)

!... BROY species
    brox_species = (/ id_br, id_bro, id_brcl, id_hobr /)
    broy_species = (/ id_br, id_bro, id_hbr, id_brono2, id_brcl, id_hobr /)
    tbry_species = (/ id_br, id_bro, id_hbr, id_brono2, id_brcl, id_hobr, id_cf2clbr, id_cf3br, id_ch3br, id_h1202, &
                      id_h2402, id_ch2br2, id_chbr3 /)

    sox_species = (/ id_so2, id_so4, id_h2so4 /)
    nhx_species = (/ id_nh3, id_nh4, id_nh4no3 /)
    bulkaero_species(:) = -1
    bulkaero_species(1:20) = (/ id_dst01, id_dst02, id_dst03, id_dst04, &
                                id_sslt01, id_sslt02, id_sslt03, id_sslt04, &
                                id_soa, id_so4, id_oc1, id_oc2, id_cb1, id_cb2, id_nh4no3, &
                                id_soam,id_soai,id_soat,id_soab,id_soax /)

    aer_species(:) = -1
    n = 1
    do m = 1,gas_pcnst
       k=0
       if ( any(bulkaero_species(:)==m) ) k=1
       if ( k==0 ) k = index(trim(solsym(m)), '_a')
       if ( k==0 ) k = index(trim(solsym(m)), '_c')
       if ( k>0 ) then ! must be aerosol species
          aer_species(n) = m
          n = n+1
       endif
    enddo

    toth_species = (/ id_ch4, id_h2o, id_h2 /)

    call addfld( 'NOX', 'mol/mol', pver, 'A', 'nox (N+NO+NO2)',  phys_decomp )
    call addfld( 'NOY', 'mol/mol', pver, 'A', &
                 'noy = total inorganic nitrogen (N+NO+NO2+NO3+2N2O5+HNO3+HO2NO2+ORGNOY+NH4NO3', &
                 phys_decomp )
    call addfld( 'NOY_SRF', 'mol/mol', 1, 'A', 'surface noy volume mixing ratio',  phys_decomp )

    call addfld( 'BROX','mol/mol', pver, 'A', 'brox (Br+BrO+BRCl+HOBr)', phys_decomp )
    call addfld( 'BROY','mol/mol', pver, 'A', 'total inorganic bromine (Br+BrO+HOBr+BrONO2+HBr+BrCl)', phys_decomp )
    call addfld( 'TBRY','mol/mol', pver, 'A', 'total Br (ORG+INORG) volume mixing ratio', phys_decomp )

    call addfld( 'CLOX','mol/mol', pver, 'A', 'clox (Cl+CLO+HOCl+2Cl2+2Cl2O2+OClO', phys_decomp )
    call addfld( 'CLOY','mol/mol', pver, 'A', 'total inorganic chlorine (Cl+ClO+2Cl2+2Cl2O2+OClO+HOCl+ClONO2+HCl+BrCl)', &
	phys_decomp )
    call addfld( 'TCLY','mol/mol', pver, 'A', 'total Cl (ORG+INORG) volume mixing ratio', phys_decomp )

    call addfld( 'FOY',' mol/mol', pver, 'A', 'total inorganic fluorine (F+HF+COFCL+2COF2)', phys_decomp )
    call addfld( 'TFY',' mol/mol', pver, 'A', 'total F (ORG+INORG) volume mixing ratio', phys_decomp )

    call addfld( 'TOTH','mol/mol', pver, 'A', 'total H2 volume mixing ratio', phys_decomp )

    call addfld( 'NOY_mmr', 'kg/kg', pver, 'A', 'NOy mass mixing ratio', phys_decomp )
    call addfld( 'SOX_mmr', 'kg/kg', pver, 'A', 'SOx mass mixing ratio', phys_decomp )
    call addfld( 'NHX_mmr', 'kg/kg', pver, 'A', 'NHx mass mixing ratio', phys_decomp )

    do j = 1,NJEUV
       write( jchar, '(I2)' ) j
       jname = 'jeuv_'//trim(adjustl(jchar))
       rid_jeuv(j) = get_rxt_ndx( trim(jname) )
    enddo

    has_jeuvs = all( rid_jeuv(:) > 0 )
    has_jno_i = rid_jno_i>0
    has_jno   = rid_jno>0

    if ( has_jeuvs ) then
       call addfld( 'PION_EUV','/cm^3/s', pver, 'I', 'total euv ionization rate', phys_decomp )
       call addfld( 'PEUV1',   '/cm^3/s', pver, 'I', '(j1+j2+j3)*o', phys_decomp )
       call addfld( 'PEUV1e',  '/cm^3/s', pver, 'I', '(j14+j15+j16)*o', phys_decomp )
       call addfld( 'PEUV2',   '/cm^3/s', pver, 'I', 'j4*n', phys_decomp )
       call addfld( 'PEUV3',   '/cm^3/s', pver, 'I', '(j5+j7+j8+j9)*o2', phys_decomp )
       call addfld( 'PEUV3e',  '/cm^3/s', pver, 'I', '(j17+j19+j20+j21)*o2', phys_decomp )
       call addfld( 'PEUV4',   '/cm^3/s', pver, 'I', '(j10+j11)*n2', phys_decomp )
       call addfld( 'PEUV4e',  '/cm^3/s', pver, 'I', '(j22+j23)*n2', phys_decomp )
       call addfld( 'PEUVN2D', '/cm^3/s', pver, 'I', '(j11+j13)*n2', phys_decomp )
       call addfld( 'PEUVN2De','/cm^3/s', pver, 'I', '(j23+j25)*n2', phys_decomp )
    endif
    if ( has_jno ) then
       call addfld( 'PJNO', '/cm^3/s', pver, 'I', 'jno*no', phys_decomp )
    endif
    if ( has_jno_i ) then
       call addfld( 'PJNO_I', '/cm^3/s', pver, 'I', 'jno_i*no', phys_decomp )
    endif
!
! CCMI
!
    call addfld( 'DO3CHM_TRP', 'kg/s ',   1, 'A', 'integrated net tendency from chem in troposphere', &
                 phys_decomp, flag_xyfill=.True. )
    call addfld( 'DO3CHM_LMS', 'kg/s ',   1, 'A', 'integrated net tendency from chem in lowermost stratosphere', &
                 phys_decomp, flag_xyfill=.True. )
!
    do m = 1,gas_pcnst

       spc_name = trim(solsym(m))

       call cnst_get_ind(spc_name, n, abort=.false. )
       if ( n > 0 ) then
          attr = cnst_longname(n)
       elseif ( trim(spc_name) == 'H2O' ) then
          attr = 'water vapor'
       else
          attr = spc_name
       endif

       depvel_name(m) = 'DV_'//trim(spc_name)
       depflx_name(m) = 'DF_'//trim(spc_name)
       dtchem_name(m) = 'D'//trim(spc_name)//'CHM'

       call addfld( depvel_name(m), 'cm/s ',   1,    'A', 'deposition velocity ', phys_decomp )
       call addfld( depflx_name(m), 'kg/m2/s', 1,    'A', 'dry deposition flux ', phys_decomp )
       call addfld( dtchem_name(m), 'kg/s ',   pver, 'A', 'net tendency from chem', phys_decomp )

       wetdep_name(m) = 'WD_'//trim(spc_name)
       wtrate_name(m) = 'WDR_'//trim(spc_name)

       call addfld( wetdep_name(m), 'kg/s ',   1,    'A', spc_name//' wet deposition', phys_decomp )
       call addfld( wtrate_name(m),   '/s ',   pver, 'A', spc_name//' wet deposition rate', phys_decomp )
       
       if (spc_name(1:3) == 'num') then
          unit_basename = ' 1'
       else
          unit_basename = 'kg'
       endif

       if ( any( aer_species == m ) ) then
          call addfld( spc_name, unit_basename//'/kg ',   pver, 'A', trim(attr)//' concentration', phys_decomp)
          call addfld( trim(spc_name)//'_SRF', unit_basename//'/kg', 1, 'A', trim(attr)//" in bottom layer", phys_decomp)
       else
          call addfld( spc_name, 'mol/mol ', pver, 'A', trim(attr)//' concentration', phys_decomp)
          call addfld( trim(spc_name)//'_SRF', 'mol/mol ', 1, 'A', trim(attr)//" in bottom layer", phys_decomp)
       endif

       if ((m /= id_cly) .and. (m /= id_bry)) then
          if (history_aerosol) then
             call add_default( spc_name, 1, ' ' )
             call add_default( trim(spc_name)//'_SRF', 1, ' ' )
          endif 
          if (history_amwg) then
             call add_default( trim(spc_name)//'_SRF', 1, ' ' )
          endif
       endif

    enddo

    call addfld( 'MASS', 'kg', pver, 'A', 'mass of grid box', phys_decomp )
    call addfld( 'AREA', 'm2', 1,    'A', 'area of grid box', phys_decomp )

    call addfld( 'WD_NOY', 'kg/s', 1, 'A', 'NOy wet deposition', phys_decomp )
    call addfld( 'DF_NOY', 'kg/m2/s', 1, 'I', 'NOy dry deposition flux ', phys_decomp )

    call addfld( 'WD_SOX', 'kg/s', 1, 'A', 'SOx wet deposition', phys_decomp )
    call addfld( 'DF_SOX', 'kg/m2/s', 1, 'I', 'SOx dry deposition flux ', phys_decomp )

    call addfld( 'WD_NHX', 'kg/s', 1, 'A', 'NHx wet deposition', phys_decomp )
    call addfld( 'DF_NHX', 'kg/m2/s', 1, 'I', 'NHx dry deposition flux ', phys_decomp )

  end subroutine chm_diags_inti

  subroutine chm_diags( lchnk, ncol, vmr, mmr, rxt_rates, invariants, depvel, depflx, mmr_tend, pdel, pmid, ltrop, pbuf )
    !--------------------------------------------------------------------
    !	... utility routine to output chemistry diagnostic variables
    !--------------------------------------------------------------------
    
    use cam_history,  only : outfld
    use constituents, only : pcnst
    use constituents, only : cnst_get_ind
    use phys_grid,    only : get_area_all_p, pcols
    use physics_buffer, only : physics_buffer_desc
!
! CCMI
!
    use cam_history_support, only : fillvalue
!
    
    implicit none

    !--------------------------------------------------------------------
    !	... dummy arguments
    !--------------------------------------------------------------------
    integer,  intent(in)  :: lchnk
    integer,  intent(in)  :: ncol
    real(r8), target, intent(in)  :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(in)  :: mmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(in)  :: rxt_rates(ncol,pver,rxntot)
    real(r8), target, intent(in)  :: invariants(ncol,pver,max(1,nfs))
    real(r8), intent(in)  :: depvel(ncol, gas_pcnst)
    real(r8), target, intent(in)  :: depflx(ncol, gas_pcnst)
    real(r8), target, intent(in)  :: mmr_tend(ncol,pver,gas_pcnst)
    real(r8), target, intent(in)  :: pdel(ncol,pver)
    real(r8), target, intent(in)  :: pmid(ncol,pver)
    integer,  target, intent(in)  :: ltrop(ncol)

    type(physics_buffer_desc), pointer :: pbuf(:)

    !--------------------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------------------
    integer     :: i,j,k, m, n
    integer :: plat
    real(r8), target :: wrk(ncol,pver)
    !      real(r8)    :: tmp(ncol,pver)
    !      real(r8)    :: m(ncol,pver)
    real(r8)    :: un2(ncol)
    
    real(r8), target, dimension(ncol,pver) :: vmr_nox, vmr_noy, vmr_clox, vmr_cloy, vmr_tcly, vmr_brox, vmr_broy, vmr_toth
    real(r8), target, dimension(ncol,pver) :: vmr_tbry, vmr_foy, vmr_tfy
    real(r8), target, dimension(ncol,pver) :: mmr_noy, mmr_sox, mmr_nhx, net_chem
    real(r8), target, dimension(ncol)      :: df_noy, df_sox, df_nhx, do3chm_trp, do3chm_lms

    real(r8), target :: area(ncol), mass(ncol,pver)
    real(r8) :: wgt
    character(len=16) :: spc_name

    !--------------------------------------------------------------------
    !	... "diagnostic" groups
    !--------------------------------------------------------------------
    vmr_nox(:ncol,:) = 0._r8
    vmr_noy(:ncol,:) = 0._r8
    vmr_clox(:ncol,:) = 0._r8
    vmr_cloy(:ncol,:) = 0._r8
    vmr_tcly(:ncol,:) = 0._r8
    vmr_brox(:ncol,:) = 0._r8
    vmr_broy(:ncol,:) = 0._r8
    vmr_tbry(:ncol,:) = 0._r8
    vmr_foy(:ncol,:)  = 0._r8
    vmr_tfy(:ncol,:)  = 0._r8
    vmr_toth(:ncol,:) = 0._r8
    mmr_noy(:ncol,:) = 0._r8
    mmr_sox(:ncol,:) = 0._r8
    mmr_nhx(:ncol,:) = 0._r8
    df_noy(:ncol) = 0._r8
    df_sox(:ncol) = 0._r8
    df_nhx(:ncol) = 0._r8


    call get_area_all_p(lchnk, ncol, area)
    call chm_diags_select_impl()
    if (.not. chm_diags_use_native_impl) then
       call chm_diags_write_proof_line('chm_diags implementation = codon')
       call chm_diags_mass_codon_wrap(ncol, area, pdel, mass)
    else
       area = area * rearth**2

       do k = 1,pver
          mass(:ncol,k) = pdel(:ncol,k) * area(:ncol) * rgrav
       enddo
    end if

    call outfld( 'AREA', area(:ncol),   ncol, lchnk )
    call outfld( 'MASS', mass(:ncol,:), ncol, lchnk )

    do m = 1,gas_pcnst

       if (.not. chm_diags_use_native_impl) then
          call chm_diags_species_codon_wrap(m, ncol, trim(dtchem_name(m)) == 'DO3CHM', vmr, mmr, depflx, &
               mmr_tend, mass, pmid, ltrop, vmr_nox, vmr_noy, vmr_clox, vmr_cloy, vmr_tcly, vmr_brox, &
               vmr_broy, vmr_toth, vmr_tbry, vmr_foy, vmr_tfy, mmr_noy, mmr_sox, mmr_nhx, df_noy, df_sox, &
               df_nhx, net_chem, do3chm_trp, do3chm_lms)
          if ( any( aer_species == m ) ) then
             call outfld( solsym(m), mmr(:ncol,:,m), ncol ,lchnk )
             call outfld( trim(solsym(m))//'_SRF', mmr(:ncol,pver,m), ncol ,lchnk )
          else
             call outfld( solsym(m), vmr(:ncol,:,m), ncol ,lchnk )
             call outfld( trim(solsym(m))//'_SRF', vmr(:ncol,pver,m), ncol ,lchnk )
          endif

          call outfld( depvel_name(m), depvel(:ncol,m), ncol ,lchnk )
          call outfld( depflx_name(m), depflx(:ncol,m), ncol ,lchnk )
          call outfld( dtchem_name(m), net_chem(:ncol,:), ncol, lchnk )
          if ( trim(dtchem_name(m)) == 'DO3CHM' ) then
             call outfld('DO3CHM_TRP',do3chm_trp(:ncol), ncol, lchnk )
             call outfld('DO3CHM_LMS',do3chm_lms(:ncol), ncol, lchnk )
          end if
          cycle
       end if

 !...FOY (counting Fluorines, not chlorines or bromines)
       if ( m == id_cfc12 .or. m == id_hcfc22 .or. m == id_cf2clbr .or. m == id_h1202 .or. m == id_hcfc142b &
            .or. m == id_cof2 ) then
          wgt = 2._r8
       elseif ( m == id_cfc113 .or. m == id_cf3br ) then
          wgt = 3._r8
       elseif ( m == id_cfc114 .or. m == id_h2402 ) then
          wgt = 4._r8
       elseif ( m == id_cfc115 ) then
          wgt = 5._r8
       else
          wgt = 1._r8
       endif
       if ( any( foy_species == m ) ) then
          vmr_foy(:ncol,:) = vmr_foy(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( tfy_species == m ) ) then
          vmr_tfy(:ncol,:) = vmr_tfy(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif

!... counting chlorine and bromines, etc... (and total H2 species)
       if ( m == id_ch4 .or. m == id_n2o5 .or. m == id_cfc12 .or. m == id_cl2 .or. m == id_cl2o2  ) then
          wgt = 2._r8
       elseif (m == id_cfc114 .or. m == id_hcfc141b .or. m == id_h1202 .or. m == id_h2402 .or. m == id_ch2br2 ) then
          wgt = 2._r8
       elseif ( m == id_cfc11 .or. m == id_cfc113 .or. m == id_ch3ccl3 .or. m == id_chbr3 ) then
          wgt = 3._r8
       elseif ( m == id_ccl4 ) then
          wgt = 4._r8
       else
          wgt = 1._r8
       endif
!...NOY
       if ( any( nox_species == m ) ) then
          vmr_nox(:ncol,:) = vmr_nox(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( noy_species == m ) ) then
          vmr_noy(:ncol,:) = vmr_noy(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
!...NOY, SOX, NHX
       if ( any( noy_species == m ) ) then
          mmr_noy(:ncol,:) = mmr_noy(:ncol,:) +  wgt * mmr(:ncol,:,m)
       endif
       if ( any( sox_species == m ) ) then
          mmr_sox(:ncol,:) = mmr_sox(:ncol,:) +  wgt * mmr(:ncol,:,m)
       endif
       if ( any( nhx_species == m ) ) then
          mmr_nhx(:ncol,:) = mmr_nhx(:ncol,:) +  wgt * mmr(:ncol,:,m)
       endif
!...CLOY
       if ( any( clox_species == m ) ) then
          vmr_clox(:ncol,:) = vmr_clox(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( cloy_species == m ) ) then
          vmr_cloy(:ncol,:) = vmr_cloy(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( tcly_species == m ) ) then
          vmr_tcly(:ncol,:) = vmr_tcly(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
!...BROY
       if ( any( brox_species == m ) ) then
          vmr_brox(:ncol,:) = vmr_brox(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( broy_species == m ) ) then
          vmr_broy(:ncol,:) = vmr_broy(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       if ( any( tbry_species == m ) ) then
          vmr_tbry(:ncol,:) = vmr_tbry(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
!...HOY
       if ( any ( toth_species == m ) ) then
          vmr_toth(:ncol,:) = vmr_toth(:ncol,:) +  wgt * vmr(:ncol,:,m)
       endif
       
       if ( any( aer_species == m ) ) then
          call outfld( solsym(m), mmr(:ncol,:,m), ncol ,lchnk )
          call outfld( trim(solsym(m))//'_SRF', mmr(:ncol,pver,m), ncol ,lchnk )
       else
          call outfld( solsym(m), vmr(:ncol,:,m), ncol ,lchnk )
          call outfld( trim(solsym(m))//'_SRF', vmr(:ncol,pver,m), ncol ,lchnk )
       endif

       call outfld( depvel_name(m), depvel(:ncol,m), ncol ,lchnk )
       call outfld( depflx_name(m), depflx(:ncol,m), ncol ,lchnk )

       if ( any( noy_species == m ) ) then
          df_noy(:ncol) = df_noy(:ncol) +  wgt * depflx(:ncol,m)*N_molwgt/adv_mass(m)
       endif
       if ( any( sox_species == m ) ) then
          df_sox(:ncol) = df_sox(:ncol) +  wgt * depflx(:ncol,m)*S_molwgt/adv_mass(m)
       endif
       if ( any( nhx_species == m ) ) then
          df_nhx(:ncol) = df_nhx(:ncol) +  wgt * depflx(:ncol,m)*N_molwgt/adv_mass(m)
       endif

       do k=1,pver
          do i=1,ncol
             net_chem(i,k) = mmr_tend(i,k,m) * mass(i,k) 
          end do
       end do
       call outfld( dtchem_name(m), net_chem(:ncol,:), ncol, lchnk )
!
! CCMI
!
       if ( trim(dtchem_name(m)) == 'DO3CHM' ) then
          do3chm_trp(:) = 0._r8
          do i=1,ncol
             do k=ltrop(i),pver
                do3chm_trp(i) = do3chm_trp(i) + net_chem(i,k)
             end do
          end do
          where ( do3chm_trp == 0._r8 )
             do3chm_trp = fillvalue
          end where
          call outfld('DO3CHM_TRP',do3chm_trp(:ncol), ncol, lchnk )
          do3chm_lms(:) = 0._r8
          do i=1,ncol
             do k=1,pver
                if ( pmid(i,k) > 100.e2_r8 .and. k < ltrop(i) ) then
                   do3chm_lms(i) = do3chm_lms(i) + net_chem(i,k)
                end if
             end do
          end do
          where ( do3chm_lms == 0._r8 )
             do3chm_lms = fillvalue
          end where
          call outfld('DO3CHM_LMS',do3chm_lms(:ncol), ncol, lchnk )
       end if
!
    enddo


    call outfld( 'NOX',  vmr_nox  (:ncol,:), ncol, lchnk )
    call outfld( 'NOY',  vmr_noy  (:ncol,:), ncol, lchnk )
    call outfld( 'NOY_SRF',  vmr_noy(:ncol,pver),  ncol, lchnk )
    call outfld( 'CLOX', vmr_clox (:ncol,:), ncol, lchnk )
    call outfld( 'CLOY', vmr_cloy (:ncol,:), ncol, lchnk )
    call outfld( 'BROX', vmr_brox (:ncol,:), ncol, lchnk )
    call outfld( 'BROY', vmr_broy (:ncol,:), ncol, lchnk )
    call outfld( 'TCLY', vmr_tcly (:ncol,:), ncol, lchnk )
    call outfld( 'TBRY', vmr_tbry (:ncol,:), ncol, lchnk )
    call outfld( 'FOY',  vmr_foy  (:ncol,:), ncol, lchnk )
    call outfld( 'TFY',  vmr_tfy  (:ncol,:), ncol, lchnk )
    call outfld( 'TOTH', vmr_toth (:ncol,:), ncol, lchnk )

    call outfld( 'NOY_mmr', mmr_noy(:ncol,:), ncol ,lchnk )
    call outfld( 'SOX_mmr', mmr_sox(:ncol,:), ncol ,lchnk )
    call outfld( 'NHX_mmr', mmr_nhx(:ncol,:), ncol ,lchnk )
    call outfld( 'DF_NOY', df_noy(:ncol), ncol ,lchnk )
    call outfld( 'DF_SOX', df_sox(:ncol), ncol ,lchnk )
    call outfld( 'DF_NHX', df_nhx(:ncol), ncol ,lchnk )

    !--------------------------------------------------------------------
    !	... euv ion production
    !--------------------------------------------------------------------

    jeuvs: if ( has_jeuvs ) then
       if (.not. chm_diags_use_native_impl) then
          call chm_diags_euv_codon_wrap(1, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PION_EUV', wrk, ncol, lchnk )
          call chm_diags_euv_codon_wrap(2, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PEUV1', wrk, ncol, lchnk )
          call chm_diags_euv_codon_wrap(3, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PEUV1e', wrk, ncol, lchnk )
          call chm_diags_euv_codon_wrap(4, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PEUV2', wrk, ncol, lchnk )
          call chm_diags_euv_codon_wrap(5, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PEUV3', wrk, ncol, lchnk )
          call chm_diags_euv_codon_wrap(6, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PEUV3e', wrk, ncol, lchnk )
          call chm_diags_euv_codon_wrap(7, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PEUV4', wrk, ncol, lchnk )
          call chm_diags_euv_codon_wrap(8, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PEUV4e', wrk, ncol, lchnk )
          call chm_diags_euv_codon_wrap(9, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PEUVN2D', wrk, ncol, lchnk )
          call chm_diags_euv_codon_wrap(10, ncol, 0, vmr, rxt_rates, invariants, wrk)
          call outfld( 'PEUVN2De', wrk, ncol, lchnk )
       else
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = vmr(:,k,id_o)*(rxt_rates(:,k,rid_jeuv(1)) + rxt_rates(:,k,rid_jeuv(2)) &
               + rxt_rates(:,k,rid_jeuv(3)) + rxt_rates(:,k,rid_jeuv(14)) &
               + rxt_rates(:,k,rid_jeuv(15)) + rxt_rates(:,k,rid_jeuv(16))) &
               + vmr(:,k,id_n)*rxt_rates(:,k,rid_jeuv(4)) &
               + vmr(:,k,id_o2)*(rxt_rates(:,k,rid_jeuv(5)) + rxt_rates(:,k,rid_jeuv(7)) &
               + rxt_rates(:,k,rid_jeuv(8)) + rxt_rates(:,k,rid_jeuv(9)) &
               + rxt_rates(:,k,rid_jeuv(17)) + rxt_rates(:,k,rid_jeuv(19)) &
               + rxt_rates(:,k,rid_jeuv(20)) + rxt_rates(:,k,rid_jeuv(21))) &
               + un2(:)*(rxt_rates(:,k,rid_jeuv(6)) + rxt_rates(:,k,rid_jeuv(10)) &
               + rxt_rates(:,k,rid_jeuv(11)) + rxt_rates(:,k,rid_jeuv(18)) &
               + rxt_rates(:,k,rid_jeuv(22)) + rxt_rates(:,k,rid_jeuv(23)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PION_EUV', wrk, ncol, lchnk )

       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_o)*(rxt_rates(:,k,rid_jeuv(1)) + rxt_rates(:,k,rid_jeuv(2)) &
               + rxt_rates(:,k,rid_jeuv(3)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV1', wrk, ncol, lchnk )
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_o)*(rxt_rates(:,k,rid_jeuv(14)) + rxt_rates(:,k,rid_jeuv(15)) &
               + rxt_rates(:,k,rid_jeuv(16)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV1e', wrk, ncol, lchnk )
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_n)*rxt_rates(:,k,rid_jeuv(4))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV2', wrk, ncol, lchnk )
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_o2)*(rxt_rates(:,k,rid_jeuv(5)) + rxt_rates(:,k,rid_jeuv(7)) &
               + rxt_rates(:,k,rid_jeuv(8)) + rxt_rates(:,k,rid_jeuv(9)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV3', wrk, ncol, lchnk )
       do k = 1,pver
          wrk(:,k) = vmr(:,k,id_o2)*(rxt_rates(:,k,rid_jeuv(17)) + rxt_rates(:,k,rid_jeuv(19)) &
               + rxt_rates(:,k,rid_jeuv(20)) + rxt_rates(:,k,rid_jeuv(21)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV3e', wrk, ncol, lchnk )
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = un2(:)*(rxt_rates(:,k,rid_jeuv(6)) + rxt_rates(:,k,rid_jeuv(10)) + rxt_rates(:,k,rid_jeuv(11)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV4', wrk, ncol, lchnk )
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = un2(:)*(rxt_rates(:,k,rid_jeuv(18)) + rxt_rates(:,k,rid_jeuv(22)) + rxt_rates(:,k,rid_jeuv(23)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUV4e', wrk, ncol, lchnk )
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = un2(:)*(rxt_rates(:,k,rid_jeuv(11)) + rxt_rates(:,k,rid_jeuv(13)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUVN2D', wrk, ncol, lchnk )
       do k = 1,pver
          un2(:)   = 1._r8 - (vmr(:,k,id_o) + vmr(:,k,id_o2) + vmr(:,k,id_h))
          wrk(:,k) = un2(:)*(rxt_rates(:,k,rid_jeuv(23)) + rxt_rates(:,k,rid_jeuv(25)))
          wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
       end do
       call outfld( 'PEUVN2De', wrk, ncol, lchnk )
       end if
    endif jeuvs

    if ( has_jno_i ) then
       if (.not. chm_diags_use_native_impl) then
          call chm_diags_euv_codon_wrap(11, ncol, rid_jno_i, vmr, rxt_rates, invariants, wrk)
       else
          do k = 1,pver
             wrk(:,k) = vmr(:,k,id_no)*rxt_rates(:,k,rid_jno_i)
             wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
          end do
       end if
       call outfld( 'PJNO_I', wrk, ncol, lchnk )
    endif
    if ( has_jno ) then
       if (.not. chm_diags_use_native_impl) then
          call chm_diags_euv_codon_wrap(11, ncol, rid_jno, vmr, rxt_rates, invariants, wrk)
       else
          do k = 1,pver
             wrk(:,k) = vmr(:,k,id_no)*rxt_rates(:,k,rid_jno)
             wrk(:,k) = wrk(:,k) * invariants(:,k,indexm)
          end do
       end if
       call outfld( 'PJNO', wrk, ncol, lchnk )
    endif

  end subroutine chm_diags

  subroutine het_diags( het_rates, mmr, pdel, lchnk, ncol )

    use cam_history,  only : outfld
    use phys_grid,    only : get_wght_all_p
    implicit none

    integer,  intent(in)  :: lchnk
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: het_rates(ncol,pver,max(1,gas_pcnst))
    real(r8), intent(in)  :: mmr(ncol,pver,gas_pcnst)
    real(r8), intent(in)  :: pdel(ncol,pver)

    real(r8), dimension(ncol) :: noy_wk, sox_wk, nhx_wk, wrk_wd
    integer :: m, k, j
    integer :: plat
    real(r8) :: wght(ncol)
    !
    ! output integrated wet deposition field
    !
    noy_wk(:) = 0._r8
    sox_wk(:) = 0._r8
    nhx_wk(:) = 0._r8

    call get_wght_all_p(lchnk, ncol, wght)

    do m = 1,gas_pcnst
       !
       ! compute vertical integral
       !
       wrk_wd(:ncol) = 0._r8
       do k = 1,pver
          wrk_wd(:ncol) = wrk_wd(:ncol) + het_rates(:ncol,k,m) * mmr(:ncol,k,m) * pdel(:ncol,k) 
       end do
       !
       wrk_wd(:ncol) = wrk_wd(:ncol) * rgrav * wght(:ncol) * rearth**2
       !

       call outfld( wetdep_name(m), wrk_wd(:ncol),               ncol, lchnk )
       call outfld( wtrate_name(m), het_rates(:ncol,:,m), ncol, lchnk )

       if ( any(noy_species == m ) ) then
          noy_wk(:ncol) = noy_wk(:ncol) + wrk_wd(:ncol)*N_molwgt/adv_mass(m)
       endif
       if ( m == id_n2o5 ) then  ! 2 NOy molecules in N2O5
          noy_wk(:ncol) = noy_wk(:ncol) + wrk_wd(:ncol)*N_molwgt/adv_mass(m)
       endif
       if ( any(sox_species == m ) ) then
          sox_wk(:ncol) = sox_wk(:ncol) + wrk_wd(:ncol)*S_molwgt/adv_mass(m)
       endif
       if ( any(nhx_species == m ) ) then
          nhx_wk(:ncol) = nhx_wk(:ncol) + wrk_wd(:ncol)*N_molwgt/adv_mass(m)
       endif

    end do
    
    call outfld( 'WD_NOY', noy_wk(:ncol), ncol, lchnk )
    call outfld( 'WD_SOX', sox_wk(:ncol), ncol, lchnk )
    call outfld( 'WD_NHX', nhx_wk(:ncol), ncol, lchnk )

  end subroutine het_diags

  subroutine chm_diags_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (chm_diags_impl_selected) return

    impl_name = 'native'
    call cam_codon_get_impl('CHM_DIAGS_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       chm_diags_use_native_impl = trim(adjustl(impl_name(:n))) /= 'codon'
    else
       chm_diags_use_native_impl = .true.
    end if

    chm_diags_impl_selected = .true.

    if (masterproc) then
       if (chm_diags_use_native_impl) then
          write(iulog,*) 'chm_diags implementation = native'
       else
          write(iulog,*) 'chm_diags implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine chm_diags_select_impl

  subroutine chm_diags_write_proof_line(line)

    character(len=*), intent(in) :: line
    character(len=512) :: proof_file
    integer :: status, n, proof_unit, ios

    if (.not. masterproc .or. chm_diags_proof_written) return

    write(iulog,*) trim(line)
    proof_file = ''
    call get_environment_variable('CHM_DIAGS_PROOF_FILE', value=proof_file, length=n, status=status)
    if (status == 0 .and. n > 0) then
       open(newunit=proof_unit, file=trim(proof_file(:n)), status='old', position='append', action='write', iostat=ios)
       if (ios /= 0) then
          open(newunit=proof_unit, file=trim(proof_file(:n)), status='replace', action='write', iostat=ios)
       end if
       if (ios == 0) then
          write(proof_unit,'(A)') trim(line)
          close(proof_unit)
       end if
    end if
    chm_diags_proof_written = .true.
    call flush(iulog)

  end subroutine chm_diags_write_proof_line

  subroutine chm_diags_mass_codon_wrap(ncol, area, pdel, mass)

    integer, intent(in) :: ncol
    real(r8), target, intent(inout) :: area(ncol)
    real(r8), target, intent(in) :: pdel(ncol,pver)
    real(r8), target, intent(out) :: mass(ncol,pver)

    interface
       subroutine chm_diags_mass_codon(ncol_c, pver_c, rgrav_c, rearth_c, area_p, pdel_p, mass_p) &
            bind(c, name="chm_diags_mass_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c
         real(c_double), value :: rgrav_c, rearth_c
         type(c_ptr), value :: area_p, pdel_p, mass_p
       end subroutine chm_diags_mass_codon
    end interface

    call chm_diags_mass_codon(int(ncol, c_int64_t), int(pver, c_int64_t), real(rgrav, c_double), &
         real(rearth, c_double), c_loc(area), c_loc(pdel), c_loc(mass))

  end subroutine chm_diags_mass_codon_wrap

  subroutine chm_diags_species_codon_wrap(m, ncol, do3chm, vmr, mmr, depflx, mmr_tend, mass, pmid, ltrop, &
       vmr_nox, vmr_noy, vmr_clox, vmr_cloy, vmr_tcly, vmr_brox, vmr_broy, vmr_toth, vmr_tbry, vmr_foy, &
       vmr_tfy, mmr_noy, mmr_sox, mmr_nhx, df_noy, df_sox, df_nhx, net_chem, do3chm_trp, do3chm_lms)

    integer, intent(in) :: m, ncol
    logical, intent(in) :: do3chm
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst), mmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: depflx(ncol,gas_pcnst), mmr_tend(ncol,pver,gas_pcnst)
    real(r8), target, intent(in) :: mass(ncol,pver), pmid(ncol,pver)
    integer, target, intent(in) :: ltrop(ncol)
    real(r8), target, intent(inout) :: vmr_nox(ncol,pver), vmr_noy(ncol,pver), vmr_clox(ncol,pver)
    real(r8), target, intent(inout) :: vmr_cloy(ncol,pver), vmr_tcly(ncol,pver), vmr_brox(ncol,pver)
    real(r8), target, intent(inout) :: vmr_broy(ncol,pver), vmr_toth(ncol,pver), vmr_tbry(ncol,pver)
    real(r8), target, intent(inout) :: vmr_foy(ncol,pver), vmr_tfy(ncol,pver), mmr_noy(ncol,pver)
    real(r8), target, intent(inout) :: mmr_sox(ncol,pver), mmr_nhx(ncol,pver), df_noy(ncol)
    real(r8), target, intent(inout) :: df_sox(ncol), df_nhx(ncol)
    real(r8), target, intent(out) :: net_chem(ncol,pver), do3chm_trp(ncol), do3chm_lms(ncol)

    integer(c_int64_t), target :: ids(21), flags(13)
    real(c_double), target :: adv_mass_c(gas_pcnst)

    interface
       subroutine chm_diags_species_packed_codon(ncol_c, pver_c, gas_pcnst_c, m_c, do3chm_c, fillvalue_c, &
            n_molwgt_c, s_molwgt_c, ids_p, flags_p, vmr_p, mmr_p, depflx_p, mmr_tend_p, mass_p, pmid_p, &
            adv_mass_p, ltrop_p, vmr_nox_p, vmr_noy_p, vmr_clox_p, vmr_cloy_p, vmr_tcly_p, vmr_brox_p, &
            vmr_broy_p, vmr_toth_p, vmr_tbry_p, vmr_foy_p, vmr_tfy_p, mmr_noy_p, mmr_sox_p, mmr_nhx_p, &
            df_noy_p, df_sox_p, df_nhx_p, net_chem_p, do3chm_trp_p, do3chm_lms_p) &
            bind(c, name="chm_diags_species_packed_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, gas_pcnst_c, m_c, do3chm_c
         real(c_double), value :: fillvalue_c, n_molwgt_c, s_molwgt_c
         type(c_ptr), value :: ids_p, flags_p, vmr_p, mmr_p, depflx_p, mmr_tend_p, mass_p, pmid_p
         type(c_ptr), value :: adv_mass_p, ltrop_p, vmr_nox_p, vmr_noy_p, vmr_clox_p, vmr_cloy_p
         type(c_ptr), value :: vmr_tcly_p, vmr_brox_p, vmr_broy_p, vmr_toth_p, vmr_tbry_p, vmr_foy_p
         type(c_ptr), value :: vmr_tfy_p, mmr_noy_p, mmr_sox_p, mmr_nhx_p, df_noy_p, df_sox_p, df_nhx_p
         type(c_ptr), value :: net_chem_p, do3chm_trp_p, do3chm_lms_p
       end subroutine chm_diags_species_packed_codon
    end interface

    ids(:) = int((/ id_ch4, id_n2o5, id_cfc12, id_cl2, id_cl2o2, id_cfc114, id_hcfc141b, id_h1202, &
         id_h2402, id_ch2br2, id_cfc11, id_cfc113, id_ch3ccl3, id_chbr3, id_ccl4, id_hcfc22, &
         id_cf2clbr, id_hcfc142b, id_cof2, id_cf3br, id_cfc115 /), c_int64_t)
    flags(:) = int((/ merge(1, 0, any(foy_species == m)), merge(1, 0, any(tfy_species == m)), &
         merge(1, 0, any(nox_species == m)), merge(1, 0, any(noy_species == m)), &
         merge(1, 0, any(sox_species == m)), merge(1, 0, any(nhx_species == m)), &
         merge(1, 0, any(clox_species == m)), merge(1, 0, any(cloy_species == m)), &
         merge(1, 0, any(tcly_species == m)), merge(1, 0, any(brox_species == m)), &
         merge(1, 0, any(broy_species == m)), merge(1, 0, any(tbry_species == m)), &
         merge(1, 0, any(toth_species == m)) /), c_int64_t)
    adv_mass_c(:) = real(adv_mass(:), c_double)

    call chm_diags_species_packed_codon(int(ncol, c_int64_t), int(pver, c_int64_t), &
         int(gas_pcnst, c_int64_t), int(m, c_int64_t), int(merge(1, 0, do3chm), c_int64_t), &
         real(fillvalue, c_double), real(N_molwgt, c_double), real(S_molwgt, c_double), c_loc(ids), &
         c_loc(flags), c_loc(vmr), c_loc(mmr), c_loc(depflx), c_loc(mmr_tend), c_loc(mass), c_loc(pmid), &
         c_loc(adv_mass_c), c_loc(ltrop), c_loc(vmr_nox), c_loc(vmr_noy), c_loc(vmr_clox), c_loc(vmr_cloy), &
         c_loc(vmr_tcly), c_loc(vmr_brox), c_loc(vmr_broy), c_loc(vmr_toth), c_loc(vmr_tbry), c_loc(vmr_foy), &
         c_loc(vmr_tfy), c_loc(mmr_noy), c_loc(mmr_sox), c_loc(mmr_nhx), c_loc(df_noy), c_loc(df_sox), &
         c_loc(df_nhx), c_loc(net_chem), c_loc(do3chm_trp), c_loc(do3chm_lms))

  end subroutine chm_diags_species_codon_wrap

  subroutine chm_diags_euv_codon_wrap(stage, ncol, rid_scalar, vmr, rxt_rates, invariants, wrk)

    integer, intent(in) :: stage, ncol, rid_scalar
    real(r8), target, intent(in) :: vmr(ncol,pver,gas_pcnst), rxt_rates(ncol,pver,rxntot)
    real(r8), target, intent(in) :: invariants(ncol,pver,max(1,nfs))
    real(r8), target, intent(out) :: wrk(ncol,pver)
    integer(c_int64_t), target :: rid_jeuv_c(NJEUV)

    interface
       subroutine chm_diags_euv_codon(stage_c, ncol_c, pver_c, indexm_c, id_o_c, id_o2_c, id_h_c, id_n_c, &
            id_no_c, rid_scalar_c, rid_jeuv_p, vmr_p, rxt_rates_p, invariants_p, wrk_p) &
            bind(c, name="chm_diags_euv_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pver_c, indexm_c, id_o_c, id_o2_c, id_h_c
         integer(c_int64_t), value :: id_n_c, id_no_c, rid_scalar_c
         type(c_ptr), value :: rid_jeuv_p, vmr_p, rxt_rates_p, invariants_p, wrk_p
       end subroutine chm_diags_euv_codon
    end interface

    rid_jeuv_c(:) = int(rid_jeuv(:), c_int64_t)
    call chm_diags_euv_codon(int(stage, c_int64_t), int(ncol, c_int64_t), int(pver, c_int64_t), &
         int(indexm, c_int64_t), int(id_o, c_int64_t), int(id_o2, c_int64_t), int(id_h, c_int64_t), &
         int(id_n, c_int64_t), int(id_no, c_int64_t), int(rid_scalar, c_int64_t), c_loc(rid_jeuv_c), &
         c_loc(vmr), c_loc(rxt_rates), c_loc(invariants), c_loc(wrk))

  end subroutine chm_diags_euv_codon_wrap

end module mo_chm_diags
