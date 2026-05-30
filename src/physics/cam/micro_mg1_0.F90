module micro_mg1_0

!---------------------------------------------------------------------------------
! Purpose:
!   MG microphysics
!
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!
! NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
! microphysics in cooperation with the MG liquid microphysics. This is
! controlled by the do_cldice variable.
!
! NOTE: If do_cldice is false, then MG microphysics should not update CLDICE
! or NUMICE; however, it is assumed that the other microphysics scheme will have
! updated CLDICE and NUMICE. The other microphysics should handle the following
! processes that would have been done by MG:
!   - Detrainment (liquid and ice)
!   - Homogeneous ice nucleation
!   - Heterogeneous ice nucleation
!   - Bergeron process
!   - Melting of ice
!   - Freezing of cloud drops
!   - Autoconversion (ice -> snow)
!   - Growth/Sublimation of ice
!   - Sedimentation of ice
!---------------------------------------------------------------------------------
! modification for sub-columns, HM, (orig 8/11/10)
! This is done using the logical 'microp_uniform' set to .true. = uniform for subcolumns
!---------------------------------------------------------------------------------

! Procedures required:
! 1) An implementation of the gamma function (if not intrinsic).
! 2) saturation vapor pressure to specific humidity formula
! 3) svp over water
! 4) svp over ice

#ifndef HAVE_GAMMA_INTRINSICS
use shr_spfn_mod, only: gamma => shr_spfn_gamma
#endif

  use wv_sat_methods, only: &
       svp_water => wv_sat_svp_water, &
       svp_ice => wv_sat_svp_ice, &
       svp_to_qsat => wv_sat_svp_to_qsat

  use phys_control, only: phys_getopts
  use spmd_utils, only: masterproc
  use cam_logfile, only: iulog

implicit none
private
save

interface
   function micro_mg_tend_codon(stage_c) result(stage_out) bind(c, name="micro_mg_tend_codon")
     use iso_c_binding, only: c_int64_t
     integer(c_int64_t), value :: stage_c
     integer(c_int64_t) :: stage_out
   end function micro_mg_tend_codon
end interface

! Note: The liu_in option has been removed, as there was a serious bug with this
! option being set to false. The code now behaves as if the default liu_in=.true.
! is always on. Addition/reinstatement of ice nucleation options will likely be
! done outside of this module.

public :: &
     micro_mg_init, &
     micro_mg_get_cols, &
     micro_mg_tend

integer, parameter :: r8 = selected_real_kind(12)      ! 8 byte real

real(r8) :: g              !gravity
real(r8) :: r              !Dry air Gas constant
real(r8) :: rv             !water vapor gas contstant
real(r8) :: cpp            !specific heat of dry air
real(r8) :: rhow           !density of liquid water
real(r8) :: tmelt          ! Freezing point of water (K)
real(r8) :: xxlv           ! latent heat of vaporization
real(r8) :: xlf            !latent heat of freezing
real(r8) :: xxls           !latent heat of sublimation

real(r8) :: rhosn  ! bulk density snow
real(r8) :: rhoi   ! bulk density ice

real(r8) :: ac,bc,as,bs,ai,bi,ar,br  !fall speed parameters 
real(r8) :: ci,di    !ice mass-diameter relation parameters
real(r8) :: cs,ds    !snow mass-diameter relation parameters
real(r8) :: cr,dr    !drop mass-diameter relation parameters
real(r8) :: f1s,f2s  !ventilation param for snow
real(r8) :: Eii      !collection efficiency aggregation of ice
real(r8) :: Ecr      !collection efficiency cloud droplets/rain
real(r8) :: f1r,f2r  !ventilation param for rain
real(r8) :: DCS      !autoconversion size threshold
real(r8) :: qsmall   !min mixing ratio 
real(r8) :: bimm,aimm !immersion freezing
real(r8) :: rhosu     !typical 850mn air density
real(r8) :: mi0       ! new crystal mass
real(r8) :: rin       ! radius of contact nuclei
real(r8) :: pi       ! pi

! Additional constants to help speed up code

real(r8) :: cons1
real(r8) :: cons4
real(r8) :: cons5
real(r8) :: cons6
real(r8) :: cons7
real(r8) :: cons8
real(r8) :: cons11
real(r8) :: cons13
real(r8) :: cons14
real(r8) :: cons16
real(r8) :: cons17
real(r8) :: cons22
real(r8) :: cons23
real(r8) :: cons24
real(r8) :: cons25
real(r8) :: cons27
real(r8) :: cons28

real(r8) :: lammini
real(r8) :: lammaxi
real(r8) :: lamminr
real(r8) :: lammaxr
real(r8) :: lammins
real(r8) :: lammaxs

! parameters for snow/rain fraction for convective clouds
real(r8) :: tmax_fsnow ! max temperature for transition to convective snow
real(r8) :: tmin_fsnow ! min temperature for transition to convective snow

!needed for findsp
real(r8) :: tt0       ! Freezing temperature

real(r8) :: csmin,csmax,minrefl,mindbz

real(r8) :: rhmini     ! Minimum rh for ice cloud fraction > 0.

logical :: use_hetfrz_classnuc ! option to use heterogeneous freezing

character(len=16)  :: micro_mg_precip_frac_method  ! type of precipitation fraction method
real(r8)           :: micro_mg_berg_eff_factor     ! berg efficiency factor

logical :: micro_mg1_0_init_use_native_impl = .false.
logical :: micro_mg1_0_init_impl_selected = .false.
logical :: micro_mg1_0_init_wrapper_logged = .false.
logical :: micro_mg1_0_init_scalars_logged = .false.
logical :: micro_mg1_0_colzero_use_native_impl = .false.
logical :: micro_mg1_0_colzero_impl_selected = .false.
logical :: micro_mg1_0_colzero_wrapper_logged = .false.
logical :: micro_mg1_0_rate1ord_logged = .false.
logical :: micro_mg1_0_substep_setup_logged = .false.
logical :: micro_mg1_0_substep_accum_logged = .false.
logical :: micro_mg1_0_incloud_activation_logged = .false.
logical :: micro_mg1_0_conservation_limiter_logged = .false.
logical :: micro_mg1_0_process_output_logged = .false.
logical :: micro_mg1_0_post_iter_avg_logged = .false.
logical :: micro_mg1_0_phase_change_logged = .false.
logical :: micro_mg1_0_number_cleanup_logged = .false.
logical :: micro_mg1_0_reflectivity_flags_logged = .false.
logical :: micro_mg1_0_tend_use_native_impl = .false.
logical :: micro_mg1_0_tend_impl_selected = .false.
logical :: micro_mg1_0_tend_logged = .false.
logical :: micro_mg1_0_tail_diag_use_native_impl = .false.
logical :: micro_mg1_0_tail_diag_impl_selected = .false.
logical :: micro_mg1_0_tail_diag_logged = .false.
logical :: micro_mg1_0_get_cols_use_native_impl = .false.
logical :: micro_mg1_0_get_cols_impl_selected = .false.
logical :: micro_mg1_0_get_cols_logged = .false.


!===============================================================================
contains
!===============================================================================

subroutine micro_mg_init( &
     kind, gravit, rair, rh2o, cpair,  &
     rhoh2o, tmelt_in, latvap, latice, &
     rhmini_in, micro_mg_dcs, use_hetfrz_classnuc_in, &
     micro_mg_precip_frac_method_in, micro_mg_berg_eff_factor_in, errstring)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! initialize constants for the morrison microphysics
! 
! Author: Andrew Gettelman Dec 2005
! 
!-----------------------------------------------------------------------

use iso_c_binding, only: c_loc

integer,          intent(in)  :: kind            ! Kind used for reals
real(r8),         intent(in)  :: gravit
real(r8),         intent(in)  :: rair
real(r8),         intent(in)  :: rh2o
real(r8),         intent(in)  :: cpair
real(r8),         intent(in)  :: rhoh2o
real(r8),         intent(in)  :: tmelt_in        ! Freezing point of water (K)
real(r8),         intent(in)  :: latvap
real(r8),         intent(in)  :: latice
real(r8),         intent(in)  :: rhmini_in       ! Minimum rh for ice cloud fraction > 0.
real(r8),         intent(in)  :: micro_mg_dcs
logical,          intent(in)  :: use_hetfrz_classnuc_in
character(len=16),intent(in)  :: micro_mg_precip_frac_method_in  ! type of precipitation fraction method
real(r8),         intent(in)  :: micro_mg_berg_eff_factor_in     ! berg efficiency factor

character(128),   intent(out) :: errstring       ! Output status (non-blank for error return)

integer k

integer l,m, iaer
real(r8) surften       ! surface tension of water w/respect to air (N/m)
real(r8) arg
real(r8), target :: init_scalars(55)

interface
   subroutine micro_mg1_0_init_scalars_codon(gravit_c, rair_c, rh2o_c, &
        cpair_c, rhoh2o_c, tmelt_c, rhmini_c, berg_eff_c, latvap_c, &
        latice_c, micro_mg_dcs_c, scalars_p) bind(c, name="micro_mg1_0_init_scalars_codon")
      use iso_c_binding, only: c_double, c_ptr
      real(c_double), value :: gravit_c, rair_c, rh2o_c, cpair_c
      real(c_double), value :: rhoh2o_c, tmelt_c, rhmini_c, berg_eff_c
      real(c_double), value :: latvap_c, latice_c, micro_mg_dcs_c
      type(c_ptr), value :: scalars_p
   end subroutine micro_mg1_0_init_scalars_codon
end interface
!-----------------------------------------------------------------------

errstring = ' '

if( kind .ne. r8 ) then
   errstring = 'micro_mg_init: KIND of reals does not match'
   return
end if

!declarations for morrison codes (transforms variable names)

call micro_mg1_0_select_init_impl()
if (micro_mg1_0_init_use_native_impl) then
   g= gravit                  !gravity
   r= rair                    !Dry air Gas constant: note units(phys_constants are in J/K/kmol)
   rv= rh2o                   !water vapor gas contstant
   cpp = cpair                !specific heat of dry air
   rhow = rhoh2o              !density of liquid water
   tmelt = tmelt_in
   rhmini = rhmini_in
   micro_mg_precip_frac_method = micro_mg_precip_frac_method_in
   micro_mg_berg_eff_factor    = micro_mg_berg_eff_factor_in

   ! latent heats

   xxlv = latvap         ! latent heat vaporization
   xlf = latice          ! latent heat freezing
   xxls = xxlv + xlf     ! latent heat of sublimation

   ! flags
   use_hetfrz_classnuc = use_hetfrz_classnuc_in

   ! parameters for snow/rain fraction for convective clouds

   tmax_fsnow = tmelt
   tmin_fsnow = tmelt-5._r8

   init_scalars(14) = 250._r8
   init_scalars(15) = 500._r8
   init_scalars(16) = 1000._r8
   init_scalars(17) = 3.e7_r8
   init_scalars(18) = 2._r8
   init_scalars(19) = 11.72_r8
   init_scalars(20) = 0.41_r8
   init_scalars(21) = 700._r8
   init_scalars(22) = 1._r8
   init_scalars(23) = 841.99667_r8
   init_scalars(24) = 0.8_r8
   init_scalars(25) = 3.1415927_r8
   init_scalars(26) = init_scalars(15)*init_scalars(25)/6._r8
   init_scalars(27) = 3._r8
   init_scalars(28) = init_scalars(14)*init_scalars(25)/6._r8
   init_scalars(29) = 3._r8
   init_scalars(30) = init_scalars(16)*init_scalars(25)/6._r8
   init_scalars(31) = 3._r8
   init_scalars(32) = 0.86_r8
   init_scalars(33) = 0.28_r8
   init_scalars(34) = 0.1_r8
   init_scalars(35) = 1.0_r8
   init_scalars(36) = 0.78_r8
   init_scalars(37) = 0.32_r8
   init_scalars(38) = micro_mg_dcs
   init_scalars(39) = 1.e-18_r8
   init_scalars(40) = 100._r8
   init_scalars(41) = 0.66_r8
   init_scalars(42) = 85000._r8/(rair * tmelt)
   init_scalars(43) = 0.1e-6_r8
   init_scalars(44) = 273.15_r8
   init_scalars(45) = -30._r8
   init_scalars(46) = 26._r8
   init_scalars(47) = -99._r8
   init_scalars(48) = 1.26e-10_r8
   init_scalars(49) = 1._r8/10.e-6_r8
   init_scalars(50) = 1._r8/(2._r8*micro_mg_dcs)
   init_scalars(51) = 1._r8/20.e-6_r8
   init_scalars(52) = 1._r8/500.e-6_r8
   init_scalars(53) = 1._r8/10.e-6_r8
   init_scalars(54) = 1._r8/2000.e-6_r8
   init_scalars(55) = 4._r8/3._r8*init_scalars(25)*init_scalars(15)* &
        (10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)
else
   call micro_mg1_0_log_init_scalars_entered()
   call micro_mg1_0_init_scalars_codon(gravit, rair, rh2o, cpair, rhoh2o, &
        tmelt_in, rhmini_in, micro_mg_berg_eff_factor_in, latvap, latice, &
        micro_mg_dcs, c_loc(init_scalars(1)))
   g= init_scalars(1)         !gravity
   r= init_scalars(2)         !Dry air Gas constant: note units(phys_constants are in J/K/kmol)
   rv= init_scalars(3)        !water vapor gas contstant
   cpp = init_scalars(4)      !specific heat of dry air
   rhow = init_scalars(5)     !density of liquid water
   tmelt = init_scalars(6)
   rhmini = init_scalars(7)
   micro_mg_precip_frac_method = micro_mg_precip_frac_method_in
   micro_mg_berg_eff_factor    = init_scalars(8)

   ! latent heats

   xxlv = init_scalars(9)     ! latent heat vaporization
   xlf = init_scalars(10)     ! latent heat freezing
   xxls = init_scalars(11)    ! latent heat of sublimation

   ! flags
   use_hetfrz_classnuc = use_hetfrz_classnuc_in

   ! parameters for snow/rain fraction for convective clouds

   tmax_fsnow = init_scalars(12)
   tmin_fsnow = init_scalars(13)
end if

! parameters below from Reisner et al. (1998)
! density parameters (kg/m3)

rhosn = init_scalars(14)    ! bulk density snow  (++ ceh)
rhoi = init_scalars(15)     ! bulk density ice
rhow = init_scalars(16)     ! bulk density liquid


! fall speed parameters, V = aD^b
! V is in m/s

! droplets
ac = init_scalars(17)
bc = init_scalars(18)

! snow
as = init_scalars(19)
bs = init_scalars(20)

! cloud ice
ai = init_scalars(21)
bi = init_scalars(22)

! rain
ar = init_scalars(23)
br = init_scalars(24)

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d

pi= init_scalars(25)

! cloud ice mass-diameter relationship

ci = init_scalars(26)
di = init_scalars(27)

! snow mass-diameter relationship

cs = init_scalars(28)
ds = init_scalars(29)

! drop mass-diameter relationship

cr = init_scalars(30)
dr = init_scalars(31)

! ventilation parameters for snow
! hall and prupacher

f1s = init_scalars(32)
f2s = init_scalars(33)

! collection efficiency, aggregation of cloud ice and snow

Eii = init_scalars(34)

! collection efficiency, accretion of cloud water by rain

Ecr = init_scalars(35)

! ventilation constants for rain

f1r = init_scalars(36)
f2r = init_scalars(37)

! autoconversion size threshold for cloud ice to snow (m)

Dcs = init_scalars(38)

! smallest mixing ratio considered in microphysics

qsmall = init_scalars(39)

! immersion freezing parameters, bigg 1953

bimm = init_scalars(40)
aimm = init_scalars(41)

! typical air density at 850 mb

rhosu = init_scalars(42)

! mass of new crystal due to aerosol freezing and growth (kg)

mi0 = init_scalars(55)

! radius of contact nuclei aerosol (m)

rin = init_scalars(43)

! freezing temperature
tt0=init_scalars(44)

pi=4._r8*atan(1.0_r8)

!Range of cloudsat reflectivities (dBz) for analytic simulator
csmin= init_scalars(45)
csmax= init_scalars(46)
mindbz = init_scalars(47)
!      minrefl = 10._r8**(mindbz/10._r8)
minrefl = init_scalars(48)

! Define constants to help speed up code (limit calls to gamma function)

cons1=gamma(1._r8+di)
cons4=gamma(1._r8+br)
cons5=gamma(4._r8+br)
cons6=gamma(1._r8+ds)
cons7=gamma(1._r8+bs)     
cons8=gamma(4._r8+bs)     
cons11=gamma(3._r8+bs)
cons13=gamma(5._r8/2._r8+br/2._r8)
cons14=gamma(5._r8/2._r8+bs/2._r8)
cons16=gamma(1._r8+bi)
cons17=gamma(4._r8+bi)
cons22=(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
cons23=dcs**3
cons24=dcs**2
cons25=dcs**bs
cons27=xxlv**2
cons28=xxls**2

lammaxi = init_scalars(49)
lammini = init_scalars(50)
lammaxr = init_scalars(51)
lamminr = init_scalars(52)
lammaxs = init_scalars(53)
lammins = init_scalars(54)

end subroutine micro_mg_init

!===============================================================================
!microphysics routine for each timestep goes here...

subroutine micro_mg_tend ( &
     microp_uniform, pcols, pver, ncol, top_lev, deltatin,&
     tn, qn, qc, qi, nc,                              &
     ni, p, pdel, cldn, liqcldf,                      &
     relvar, accre_enhan,                             &
     icecldf, rate1ord_cw2pr_st, naai, npccnin,       &
     rndst, nacon, tlat, qvlat, qctend,               &
     qitend, nctend, nitend, effc, effc_fn,           &
     effi, prect, preci, nevapr, evapsnow, am_evp_st, &
     prain, prodsnow, cmeout, deffi, pgamrad,         &
     lamcrad, qsout, dsout, rflx, sflx,               &
     qrout, reff_rain, reff_snow, qcsevap, qisevap,   &
     qvres, cmeiout, vtrmc, vtrmi, qcsedten,          &
     qisedten, prao, prco, mnuccco, mnuccto,          &
     msacwio, psacwso, bergso, bergo, melto,          &
     homoo, qcreso, prcio, praio, qireso,             &
     mnuccro, pracso, meltsdt, frzrdt, mnuccdo,       &
     nrout, nsout, refl, arefl, areflz,               &
     frefl, csrfl, acsrfl, fcsrfl, rercld,            &
     ncai, ncal, qrout2, qsout2, nrout2,              &
     nsout2, drout2, dsout2, freqs, freqr,            &
     nfice, prer_evap, do_cldice, errstring,          &
     tnd_qsnow, tnd_nsnow, re_ice,                    &
     frzimm, frzcnt, frzdep, preo, prdso,             &
     frzro, meltso, wtfc, wtfi, wtprelat, wtpostlat )

use iso_c_binding, only: c_int64_t

! input arguments
logical,  intent(in) :: microp_uniform  ! True = configure uniform for sub-columns  False = use w/o sub-columns (standard)
integer,  intent(in) :: pcols                ! size of column (first) index
integer,  intent(in) :: pver                 ! number of layers in columns
integer,  intent(in) :: ncol                 ! number of columns
integer,  intent(in) :: top_lev              ! top level microphys is applied
real(r8), intent(in) :: deltatin             ! time step (s)
real(r8), intent(in) :: tn(pcols,pver)       ! input temperature (K)
real(r8), intent(in) :: qn(pcols,pver)       ! input h20 vapor mixing ratio (kg/kg)
real(r8), intent(in) :: relvar(pcols,pver)   ! relative variance of cloud water (-)
real(r8), intent(in) :: accre_enhan(pcols,pver) ! optional accretion enhancement factor (-)

! note: all input cloud variables are grid-averaged
real(r8), intent(inout) :: qc(pcols,pver)    ! cloud water mixing ratio (kg/kg)
real(r8), intent(inout) :: qi(pcols,pver)    ! cloud ice mixing ratio (kg/kg)
real(r8), intent(inout) :: nc(pcols,pver)    ! cloud water number conc (1/kg)
real(r8), intent(inout) :: ni(pcols,pver)    ! cloud ice number conc (1/kg)
real(r8), intent(in) :: p(pcols,pver)        ! air pressure (pa)
real(r8), intent(in) :: pdel(pcols,pver)     ! pressure difference across level (pa)
real(r8), intent(in) :: cldn(pcols,pver)     ! cloud fraction
real(r8), intent(in) :: icecldf(pcols,pver)  ! ice cloud fraction   
real(r8), intent(in) :: liqcldf(pcols,pver)  ! liquid cloud fraction
          
real(r8), intent(out) :: rate1ord_cw2pr_st(pcols,pver) ! 1st order rate for direct cw to precip conversion 
! used for scavenging
! Inputs for aerosol activation
real(r8), intent(in) :: naai(pcols,pver)      ! ice nulceation number (from microp_aero_ts) 
real(r8), intent(in) :: npccnin(pcols,pver)   ! ccn activated number tendency (from microp_aero_ts)
real(r8), intent(in) :: rndst(pcols,pver,4)   ! radius of 4 dust bins for contact freezing (from microp_aero_ts)
real(r8), intent(in) :: nacon(pcols,pver,4)   ! number in 4 dust bins for contact freezing  (from microp_aero_ts)

! Used with CARMA cirrus microphysics
! (or similar external microphysics model)
logical,  intent(in) :: do_cldice             ! Prognosing cldice

! output arguments

real(r8), intent(out) :: tlat(pcols,pver)    ! latent heating rate       (W/kg)
real(r8), intent(out) :: qvlat(pcols,pver)   ! microphysical tendency qv (1/s)
real(r8), intent(out) :: qctend(pcols,pver)  ! microphysical tendency qc (1/s) 
real(r8), intent(out) :: qitend(pcols,pver)  ! microphysical tendency qi (1/s)
real(r8), intent(out) :: nctend(pcols,pver)  ! microphysical tendency nc (1/(kg*s))
real(r8), intent(out) :: nitend(pcols,pver)  ! microphysical tendency ni (1/(kg*s))
real(r8), intent(out) :: effc(pcols,pver)    ! droplet effective radius (micron)
real(r8), intent(out) :: effc_fn(pcols,pver) ! droplet effective radius, assuming nc = 1.e8 kg-1
real(r8), intent(out) :: effi(pcols,pver)    ! cloud ice effective radius (micron)
real(r8), intent(out) :: prect(pcols)        ! surface precip rate (m/s)
real(r8), intent(out) :: preci(pcols)        ! cloud ice/snow precip rate (m/s)
real(r8), intent(out) :: nevapr(pcols,pver)  ! evaporation rate of rain + snow
real(r8), intent(out) :: evapsnow(pcols,pver)! sublimation rate of snow
real(r8), intent(out) :: am_evp_st(pcols,pver)! stratiform evaporation area
real(r8), intent(out) :: prain(pcols,pver)   ! production of rain + snow
real(r8), intent(out) :: prodsnow(pcols,pver)! production of snow
real(r8), intent(out) :: cmeout(pcols,pver)  ! evap/sub of cloud
real(r8), intent(out) :: deffi(pcols,pver)   ! ice effective diameter for optics (radiation)
real(r8), intent(out) :: pgamrad(pcols,pver) ! ice gamma parameter for optics (radiation)
real(r8), intent(out) :: lamcrad(pcols,pver) ! slope of droplet distribution for optics (radiation)
real(r8), intent(out) :: qsout(pcols,pver)   ! snow mixing ratio (kg/kg)
real(r8), intent(out) :: dsout(pcols,pver)   ! snow diameter (m)
real(r8), intent(out) :: rflx(pcols,pver+1)  ! grid-box average rain flux (kg m^-2 s^-1)
real(r8), intent(out) :: sflx(pcols,pver+1)  ! grid-box average snow flux (kg m^-2 s^-1)
real(r8), intent(out) :: qrout(pcols,pver)     ! grid-box average rain mixing ratio (kg/kg)
real(r8), intent(inout) :: reff_rain(pcols,pver) ! rain effective radius (micron)
real(r8), intent(inout) :: reff_snow(pcols,pver) ! snow effective radius (micron)
real(r8), intent(out) :: qcsevap(pcols,pver) ! cloud water evaporation due to sedimentation
real(r8), intent(out) :: qisevap(pcols,pver) ! cloud ice sublimation due to sublimation
real(r8), intent(out) :: qvres(pcols,pver) ! residual condensation term to ensure RH < 100%
real(r8), intent(out) :: cmeiout(pcols,pver) ! grid-mean cloud ice sub/dep
real(r8), intent(out) :: vtrmc(pcols,pver) ! mass-weighted cloud water fallspeed
real(r8), intent(out) :: vtrmi(pcols,pver) ! mass-weighted cloud ice fallspeed
real(r8), intent(out) :: qcsedten(pcols,pver) ! qc sedimentation tendency
real(r8), intent(out) :: qisedten(pcols,pver) ! qi sedimentation tendency
! microphysical process rates for output (mixing ratio tendencies)
real(r8), intent(out) :: prao(pcols,pver) ! accretion of cloud by rain 
real(r8), intent(out) :: prco(pcols,pver) ! autoconversion of cloud to rain
real(r8), intent(out) :: mnuccco(pcols,pver) ! mixing rat tend due to immersion freezing
real(r8), intent(out) :: mnuccto(pcols,pver) ! mixing ratio tend due to contact freezing
real(r8), intent(out) :: msacwio(pcols,pver) ! mixing ratio tend due to H-M splintering
real(r8), intent(out) :: psacwso(pcols,pver) ! collection of cloud water by snow
real(r8), intent(out) :: bergso(pcols,pver) ! bergeron process on snow
real(r8), intent(out) :: bergo(pcols,pver) ! bergeron process on cloud ice
real(r8), intent(out) :: melto(pcols,pver) ! melting of cloud ice
real(r8), intent(out) :: homoo(pcols,pver) ! homogeneos freezign cloud water
real(r8), intent(out) :: qcreso(pcols,pver) ! residual cloud condensation due to removal of excess supersat
real(r8), intent(out) :: prcio(pcols,pver) ! autoconversion of cloud ice to snow
real(r8), intent(out) :: praio(pcols,pver) ! accretion of cloud ice by snow
real(r8), intent(out) :: qireso(pcols,pver) ! residual ice deposition due to removal of excess supersat
real(r8), intent(out) :: mnuccro(pcols,pver) ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
real(r8), intent(out) :: pracso (pcols,pver) ! mixing ratio tendency due to accretion of rain by snow (1/s)
real(r8), intent(out) :: meltsdt(pcols,pver) ! latent heating rate due to melting of snow  (W/kg)
real(r8), intent(out) :: frzrdt (pcols,pver) ! latent heating rate due to homogeneous freezing of rain (W/kg)
real(r8), intent(out) :: mnuccdo(pcols,pver) ! mass tendency from ice nucleation
real(r8), intent(out) :: nrout(pcols,pver) ! rain number concentration (1/m3)
real(r8), intent(out) :: nsout(pcols,pver) ! snow number concentration (1/m3)
real(r8), intent(out) :: refl(pcols,pver)    ! analytic radar reflectivity        
real(r8), intent(out) :: arefl(pcols,pver)  !average reflectivity will zero points outside valid range
real(r8), intent(out) :: areflz(pcols,pver)  !average reflectivity in z.
real(r8), intent(out) :: frefl(pcols,pver)
real(r8), intent(out) :: csrfl(pcols,pver)   !cloudsat reflectivity 
real(r8), intent(out) :: acsrfl(pcols,pver)  !cloudsat average
real(r8), intent(out) :: fcsrfl(pcols,pver)
real(r8), intent(out) :: rercld(pcols,pver) ! effective radius calculation for rain + cloud
real(r8), intent(out) :: ncai(pcols,pver) ! output number conc of ice nuclei available (1/m3)
real(r8), intent(out) :: ncal(pcols,pver) ! output number conc of CCN (1/m3)
real(r8), intent(out) :: qrout2(pcols,pver)
real(r8), intent(out) :: qsout2(pcols,pver)
real(r8), intent(out) :: nrout2(pcols,pver)
real(r8), intent(out) :: nsout2(pcols,pver)
real(r8), intent(out) :: drout2(pcols,pver) ! mean rain particle diameter (m)
real(r8), intent(out) :: dsout2(pcols,pver) ! mean snow particle diameter (m)
real(r8), intent(out) :: freqs(pcols,pver)
real(r8), intent(out) :: freqr(pcols,pver)
real(r8), intent(out) :: nfice(pcols,pver)
real(r8), intent(out) :: prer_evap(pcols,pver)

real(r8) :: nevapr2(pcols,pver)

character(128),   intent(out) :: errstring       ! Output status (non-blank for error return)

! Tendencies calculated by external schemes that can replace MG's native
! process tendencies.

! Used with CARMA cirrus microphysics
! (or similar external microphysics model)
real(r8), intent(in), pointer :: tnd_qsnow(:,:) ! snow mass tendency (kg/kg/s)
real(r8), intent(in), pointer :: tnd_nsnow(:,:) ! snow number tendency (#/kg/s)
real(r8), intent(in), pointer :: re_ice(:,:)    ! ice effective radius (m)

! From external ice nucleation.
real(r8), intent(in), pointer :: frzimm(:,:) ! Number tendency due to immersion freezing (1/cm3)
real(r8), intent(in), pointer :: frzcnt(:,:) ! Number tendency due to contact freezing (1/cm3)
real(r8), intent(in), pointer :: frzdep(:,:) ! Number tendency due to deposition nucleation (1/cm3)

!Water tracer/isotopes
real(r8), intent(out) :: preo(pcols,pver)      ! mass tendency from rain evaporation
real(r8), intent(out) :: prdso(pcols,pver)     ! mass tendency from snow sublimation
real(r8), intent(out) :: frzro(pcols,pver)     ! mass tendency from freezing of rain
real(r8), intent(out) :: meltso(pcols,pver)    ! mass tendency from melting of snow 
real(r8), intent(out) :: wtfc(pcols,pver)      ! initial fall velocity of cloud liquid (?)
real(r8), intent(out) :: wtfi(pcols,pver)      ! initial fall velocity of cloud ice (?)
real(r8), intent(out) :: wtprelat(pcols,pver)  ! change in temperature due to pre-sed processes
real(r8), intent(out) :: wtpostlat(pcols,pver) ! change in temperature due to post-sed processes

! local workspace
! all units mks unless otherwise stated

! Additional constants to help speed up code
real(r8) :: cons2
real(r8) :: cons3
real(r8) :: cons9
real(r8) :: cons10
real(r8) :: cons12
real(r8) :: cons15
real(r8) :: cons18
real(r8) :: cons19
real(r8) :: cons20

! temporary variables for sub-stepping 
real(r8) :: t1(pcols,pver)
real(r8) :: q1(pcols,pver)
real(r8) :: qc1(pcols,pver)
real(r8) :: qi1(pcols,pver)
real(r8) :: nc1(pcols,pver)
real(r8) :: ni1(pcols,pver)
real(r8) :: tlat1(pcols,pver)
real(r8) :: qvlat1(pcols,pver)
real(r8) :: qctend1(pcols,pver)
real(r8) :: qitend1(pcols,pver)
real(r8) :: nctend1(pcols,pver)
real(r8) :: nitend1(pcols,pver)
real(r8) :: prect1(pcols)
real(r8) :: preci1(pcols)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(r8) :: deltat        ! sub-time step (s)
real(r8) :: omsm    ! number near unity for round-off issues
real(r8) :: dto2    ! dt/2 (s)
real(r8) :: mincld  ! minimum allowed cloud fraction
real(r8) :: q(pcols,pver) ! water vapor mixing ratio (kg/kg)
real(r8) :: t(pcols,pver) ! temperature (K)
real(r8) :: rho(pcols,pver) ! air density (kg m-3)
real(r8) :: dv(pcols,pver)  ! diffusivity of water vapor in air
real(r8) :: mu(pcols,pver)  ! viscocity of air
real(r8) :: sc(pcols,pver)  ! schmidt number
real(r8) :: kap(pcols,pver) ! thermal conductivity of air
real(r8) :: rhof(pcols,pver) ! air density correction factor for fallspeed
real(r8) :: cldmax(pcols,pver) ! precip fraction assuming maximum overlap
real(r8) :: cldm(pcols,pver)   ! cloud fraction
real(r8) :: icldm(pcols,pver)   ! ice cloud fraction
real(r8) :: lcldm(pcols,pver)   ! liq cloud fraction
real(r8) :: icwc(pcols)    ! in cloud water content (liquid+ice)
real(r8) :: calpha(pcols)  ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cbeta(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cbetah(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cgamma(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cgamah(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: rcgama(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec1(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec2(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec3(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec4(pcols) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: qtmp ! dummy qv 
real(r8) :: dum  ! temporary dummy variable

real(r8) :: cme(pcols,pver)  ! total (liquid+ice) cond/evap rate of cloud

real(r8) :: cmei(pcols,pver) ! dep/sublimation rate of cloud ice
real(r8) :: cwml(pcols,pver) ! cloud water mixing ratio
real(r8) :: cwmi(pcols,pver) ! cloud ice mixing ratio
real(r8) :: nnuccd(pver)   ! ice nucleation rate from deposition/cond.-freezing
real(r8) :: mnuccd(pver)   ! mass tendency from ice nucleation
real(r8) :: qcld              ! total cloud water
real(r8) :: lcldn(pcols,pver) ! fractional coverage of new liquid cloud
real(r8) :: lcldo(pcols,pver) ! fractional coverage of old liquid cloud
real(r8) :: nctend_mixnuc(pcols,pver)
real(r8) :: arg ! argument of erfc

! for calculation of rate1ord_cw2pr_st
real(r8), target :: qcsinksum_rate1ord(pver)   ! sum over iterations of cw to precip sink
real(r8), target :: qcsum_rate1ord(pver)    ! sum over iterations of cloud water

real(r8) :: alpha

real(r8) :: dum1,dum2   !general dummy variables

real(r8) :: npccn(pver)     ! droplet activation rate
real(r8) :: qcic(pcols,pver) ! in-cloud cloud liquid mixing ratio
real(r8) :: qiic(pcols,pver) ! in-cloud cloud ice mixing ratio
real(r8) :: qniic(pcols,pver) ! in-precip snow mixing ratio
real(r8) :: qric(pcols,pver) ! in-precip rain mixing ratio
real(r8) :: ncic(pcols,pver) ! in-cloud droplet number conc
real(r8) :: niic(pcols,pver) ! in-cloud cloud ice number conc
real(r8) :: nsic(pcols,pver) ! in-precip snow number conc
real(r8) :: nric(pcols,pver) ! in-precip rain number conc
real(r8) :: lami(pver) ! slope of cloud ice size distr
real(r8) :: n0i(pver) ! intercept of cloud ice size distr
real(r8) :: lamc(pver) ! slope of cloud liquid size distr
real(r8) :: n0c(pver) ! intercept of cloud liquid size distr
real(r8) :: lams(pver) ! slope of snow size distr
real(r8) :: n0s(pver) ! intercept of snow size distr
real(r8) :: lamr(pver) ! slope of rain size distr
real(r8) :: n0r(pver) ! intercept of rain size distr
real(r8) :: cdist1(pver) ! size distr parameter to calculate droplet freezing
! combined size of precip & cloud drops
real(r8) :: arcld(pcols,pver) ! averaging control flag
real(r8) :: Actmp  !area cross section of drops
real(r8) :: Artmp  !area cross section of rain

real(r8) :: pgam(pver) ! spectral width parameter of droplet size distr
real(r8) :: lammax  ! maximum allowed slope of size distr
real(r8) :: lammin  ! minimum allowed slope of size distr
real(r8) :: nacnt   ! number conc of contact ice nuclei
real(r8) :: mnuccc(pver) ! mixing ratio tendency due to freezing of cloud water
real(r8) :: nnuccc(pver) ! number conc tendency due to freezing of cloud water

real(r8) :: mnucct(pver) ! mixing ratio tendency due to contact freezing of cloud water
real(r8) :: nnucct(pver) ! number conc tendency due to contact freezing of cloud water
real(r8) :: msacwi(pver) ! mixing ratio tendency due to HM ice multiplication
real(r8) :: nsacwi(pver) ! number conc tendency due to HM ice multiplication

real(r8) :: prc(pver) ! qc tendency due to autoconversion of cloud droplets
real(r8) :: nprc(pver) ! number conc tendency due to autoconversion of cloud droplets
real(r8) :: nprc1(pver) ! qr tendency due to autoconversion of cloud droplets
real(r8) :: nsagg(pver) ! ns tendency due to self-aggregation of snow
real(r8) :: dc0  ! mean size droplet size distr
real(r8) :: ds0  ! mean size snow size distr (area weighted)
real(r8) :: eci  ! collection efficiency for riming of snow by droplets
real(r8) :: psacws(pver) ! mixing rat tendency due to collection of droplets by snow
real(r8) :: npsacws(pver) ! number conc tendency due to collection of droplets by snow
real(r8) :: uni ! number-weighted cloud ice fallspeed
real(r8) :: umi ! mass-weighted cloud ice fallspeed
real(r8) :: uns(pver) ! number-weighted snow fallspeed
real(r8) :: ums(pver) ! mass-weighted snow fallspeed
real(r8) :: unr(pver) ! number-weighted rain fallspeed
real(r8) :: umr(pver) ! mass-weighted rain fallspeed
real(r8) :: unc ! number-weighted cloud droplet fallspeed
real(r8) :: umc ! mass-weighted cloud droplet fallspeed
real(r8) :: pracs(pver) ! mixing rat tendency due to collection of rain by snow
real(r8) :: npracs(pver) ! number conc tendency due to collection of rain by snow
real(r8) :: mnuccr(pver) ! mixing rat tendency due to freezing of rain
real(r8) :: nnuccr(pver) ! number conc tendency due to freezing of rain
real(r8) :: pra(pver) ! mixing rat tendnency due to accretion of droplets by rain
real(r8) :: npra(pver) ! nc tendnency due to accretion of droplets by rain
real(r8) :: nragg(pver) ! nr tendency due to self-collection of rain
real(r8) :: prci(pver) ! mixing rat tendency due to autoconversion of cloud ice to snow
real(r8) :: nprci(pver) ! number conc tendency due to autoconversion of cloud ice to snow
real(r8) :: prai(pver) ! mixing rat tendency due to accretion of cloud ice by snow
real(r8) :: nprai(pver) ! number conc tendency due to accretion of cloud ice by snow
real(r8) :: qvs ! liquid saturation vapor mixing ratio
real(r8) :: qvi ! ice saturation vapor mixing ratio
real(r8) :: dqsdt ! change of sat vapor mixing ratio with temperature
real(r8) :: dqsidt ! change of ice sat vapor mixing ratio with temperature
real(r8) :: ab ! correction factor for rain evap to account for latent heat
real(r8) :: qclr ! water vapor mixing ratio in clear air
real(r8) :: abi ! correction factor for snow sublimation to account for latent heat
real(r8) :: epss ! 1/ sat relaxation timescale for snow
real(r8) :: epsr ! 1/ sat relaxation timescale for rain
real(r8) :: pre(pver) ! rain mixing rat tendency due to evaporation
real(r8) :: prds(pver) ! snow mixing rat tendency due to sublimation
real(r8) :: qce ! dummy qc for conservation check
real(r8) :: qie ! dummy qi for conservation check
real(r8) :: nce ! dummy nc for conservation check
real(r8) :: nie ! dummy ni for conservation check
real(r8) :: ratio ! parameter for conservation check
real(r8) :: dumc(pcols,pver) ! dummy in-cloud qc
real(r8) :: dumnc(pcols,pver) ! dummy in-cloud nc
real(r8) :: dumi(pcols,pver) ! dummy in-cloud qi
real(r8) :: dumni(pcols,pver) ! dummy in-cloud ni
real(r8) :: dums(pcols,pver) ! dummy in-cloud snow mixing rat
real(r8) :: dumns(pcols,pver) ! dummy in-cloud snow number conc
real(r8) :: dumr(pcols,pver) ! dummy in-cloud rain mixing rat
real(r8) :: dumnr(pcols,pver) ! dummy in-cloud rain number conc
! below are parameters for cloud water and cloud ice sedimentation calculations
real(r8) :: fr(pver)
real(r8) :: fnr(pver)
real(r8) :: fc(pver)
real(r8) :: fnc(pver)
real(r8) :: fi(pver)
real(r8) :: fni(pver)
real(r8) :: fs(pver)
real(r8) :: fns(pver)
real(r8) :: faloutr(pver)
real(r8) :: faloutnr(pver)
real(r8) :: faloutc(pver)
real(r8) :: faloutnc(pver)
real(r8) :: falouti(pver)
real(r8) :: faloutni(pver)
real(r8) :: falouts(pver)
real(r8) :: faloutns(pver)
real(r8) :: faltndr
real(r8) :: faltndnr
real(r8) :: faltndc
real(r8) :: faltndnc
real(r8) :: faltndi
real(r8) :: faltndni
real(r8) :: faltnds
real(r8) :: faltndns
real(r8) :: faltndqie
real(r8) :: faltndqce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(r8) :: relhum(pcols,pver) ! relative humidity
real(r8) :: csigma(pcols) ! parameter for cond/evap of cloud water/ice
real(r8) :: rgvm ! max fallspeed for all species
real(r8) :: arn(pcols,pver) ! air density corrected rain fallspeed parameter
real(r8) :: asn(pcols,pver) ! air density corrected snow fallspeed parameter
real(r8) :: acn(pcols,pver) ! air density corrected cloud droplet fallspeed parameter
real(r8) :: ain(pcols,pver) ! air density corrected cloud ice fallspeed parameter
real(r8) :: nsubi(pver) ! evaporation of cloud ice number
real(r8) :: nsubc(pver) ! evaporation of droplet number
real(r8) :: nsubs(pver) ! evaporation of snow number
real(r8) :: nsubr(pver) ! evaporation of rain number
real(r8) :: mtime ! factor to account for droplet activation timescale
real(r8) :: dz(pcols,pver) ! height difference across model vertical level


!! add precip flux variables for sub-stepping
real(r8) :: rflx1(pcols,pver+1)
real(r8) :: sflx1(pcols,pver+1)

! returns from function/subroutine calls
real(r8) :: tsp(pcols,pver)      ! saturation temp (K)
real(r8) :: qsp(pcols,pver)      ! saturation mixing ratio (kg/kg)
real(r8) :: qsphy(pcols,pver)      ! saturation mixing ratio (kg/kg): hybrid rh
real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
real(r8) :: esl(pcols,pver)      ! liquid sat vapor pressure (pa)
real(r8) :: esi(pcols,pver)      ! ice sat vapor pressure (pa)

! sum of source/sink terms for diagnostic precip

real(r8) :: qnitend(pcols,pver) ! snow mixing ratio source/sink term
real(r8) :: nstend(pcols,pver)  ! snow number concentration source/sink term
real(r8) :: qrtend(pcols,pver) ! rain mixing ratio source/sink term
real(r8) :: nrtend(pcols,pver)  ! rain number concentration source/sink term
real(r8) :: qrtot ! vertically-integrated rain mixing rat source/sink term
real(r8) :: nrtot ! vertically-integrated rain number conc source/sink term
real(r8) :: qstot ! vertically-integrated snow mixing rat source/sink term
real(r8) :: nstot ! vertically-integrated snow number conc source/sink term

! new terms for Bergeron process

real(r8) :: dumnnuc ! provisional ice nucleation rate (for calculating bergeron)
real(r8) :: ninew  ! provisional cloud ice number conc (for calculating bergeron)
real(r8) :: qinew ! provisional cloud ice mixing ratio (for calculating bergeron)
real(r8) :: qvl  ! liquid sat mixing ratio   
real(r8) :: epsi ! 1/ sat relaxation timecale for cloud ice
real(r8) :: prd ! provisional deposition rate of cloud ice at water sat 
real(r8) :: berg(pcols,pver) ! mixing rat tendency due to bergeron process for cloud ice
real(r8) :: bergs(pver) ! mixing rat tendency due to bergeron process for snow

!bergeron terms
real(r8) :: bergtsf   !bergeron timescale to remove all liquid
real(r8) :: rhin      !modified RH for vapor deposition

! diagnostic rain/snow for output to history
! values are in-precip (local) !!!!

real(r8) :: drout(pcols,pver)     ! rain diameter (m)

!averageed rain/snow for history
real(r8) :: dumfice

!ice nucleation, droplet activation
real(r8) :: dum2i(pcols,pver) ! number conc of ice nuclei available (1/kg)
real(r8) :: dum2l(pcols,pver) ! number conc of CCN (1/kg)
real(r8) :: ncmax
real(r8) :: nimax

real(r8) :: qcvar     ! 1/relative variance of sub-grid qc

! loop array variables
integer i,k,nstep,n, l
integer ii,kk, m

! loop variables for sub-step solution
integer iter,it,ltrue(pcols)

! used in contact freezing via dust particles
real(r8)  tcnt, viscosity, mfp
real(r8)  slip1, slip2, slip3, slip4
!        real(r8)  dfaer1, dfaer2, dfaer3, dfaer4
!        real(r8)  nacon1,nacon2,nacon3,nacon4
real(r8)  ndfaer1, ndfaer2, ndfaer3, ndfaer4
real(r8)  nslip1, nslip2, nslip3, nslip4

! used in ice effective radius
real(r8)  bbi, cci, ak, iciwc, rvi

! used in Bergeron processe and water vapor deposition
real(r8)  Tk, deles, Aprpr, Bprpr, Cice, qi0, Crate, qidep

! mean cloud fraction over the time step
real(r8)  cldmw(pcols,pver)

! used in secondary ice production
real(r8) ni_secp

! variabels to check for RH after rain evap

real(r8) :: esn
real(r8) :: qsn
real(r8) :: ttmp

!water tracers:
real(r8) :: qrtend_copy(pcols,pver) !copy of qrtend.
real(r8) :: qnitend_copy(pcols,pver)!copy of qnitend.

real(r8) :: rainrt(pcols,pver)  ! rain rate for reflectivity calculation
real(r8) :: rainrt1(pcols,pver)
real(r8) :: tmp

integer(c_int64_t) :: touch_c

real(r8) dmc,ssmc,dstrn  ! variables for modal scheme.

real(r8), parameter :: cdnl    = 0.e6_r8    ! cloud droplet number limiter

! heterogeneous freezing
real(r8) :: mnudep(pver) ! mixing ratio tendency due to deposition of water vapor
real(r8) :: nnudep(pver) ! number conc tendency due to deposition of water vapor
real(r8) :: con1 ! work cnstant
real(r8) :: r3lx ! Mean volume radius (m)
real(r8) :: mi0l
real(r8) :: frztmp

logical  :: do_clubb_sgs

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

call micro_mg1_0_select_tend_impl()

if (.not. micro_mg1_0_tend_use_native_impl) then
   touch_c = micro_mg_tend_codon(1_c_int64_t)
   if (touch_c == 1_c_int64_t) then
      call micro_mg1_0_log_tend_entered()
   end if
end if

! Return error message
errstring = ' '

if (.not. (do_cldice .or. &
     (associated(tnd_qsnow) .and. associated(tnd_nsnow) .and. associated(re_ice)))) then
   errstring = "MG's native cloud ice processes are disabled, but &
        &no replacement values were passed in."
end if

if (use_hetfrz_classnuc .and. (.not. &
     (associated(frzimm) .and. associated(frzcnt) .and. associated(frzdep)))) then
   errstring = "Hoose heterogeneous freezing is enabled, but the &
        &required tendencies were not all passed in."
end if

call phys_getopts(do_clubb_sgs_out = do_clubb_sgs)

! assign variable deltat for sub-stepping...
deltat=deltatin

! parameters for scheme
omsm=0.99999_r8
dto2=0.5_r8*deltat
mincld=0.0001_r8

call micro_mg1_0_select_init_impl()
if (micro_mg1_0_init_use_native_impl) then

! initialize  output fields for number conc qand ice nucleation
ncai(1:ncol,1:pver)=0._r8 
ncal(1:ncol,1:pver)=0._r8  

!Initialize rain size
rercld(1:ncol,1:pver)=0._r8
arcld(1:ncol,1:pver)=0._r8

!initialize radiation output variables
pgamrad(1:ncol,1:pver)=0._r8 ! liquid gamma parameter for optics (radiation)
lamcrad(1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics (radiation)
deffi  (1:ncol,1:pver)=0._r8 ! slope of droplet distribution for optics (radiation)
!initialize radiation output variables
!initialize water vapor tendency term output
qcsevap(1:ncol,1:pver)=0._r8 
qisevap(1:ncol,1:pver)=0._r8 
qvres  (1:ncol,1:pver)=0._r8 
cmeiout (1:ncol,1:pver)=0._r8
vtrmc (1:ncol,1:pver)=0._r8
vtrmi (1:ncol,1:pver)=0._r8
qcsedten (1:ncol,1:pver)=0._r8
qisedten (1:ncol,1:pver)=0._r8    

prao(1:ncol,1:pver)=0._r8 
prco(1:ncol,1:pver)=0._r8 
mnuccco(1:ncol,1:pver)=0._r8 
mnuccto(1:ncol,1:pver)=0._r8 
msacwio(1:ncol,1:pver)=0._r8 
psacwso(1:ncol,1:pver)=0._r8 
bergso(1:ncol,1:pver)=0._r8 
bergo(1:ncol,1:pver)=0._r8 
melto(1:ncol,1:pver)=0._r8 
homoo(1:ncol,1:pver)=0._r8 
qcreso(1:ncol,1:pver)=0._r8 
prcio(1:ncol,1:pver)=0._r8 
praio(1:ncol,1:pver)=0._r8 
qireso(1:ncol,1:pver)=0._r8 
mnuccro(1:ncol,1:pver)=0._r8 
pracso (1:ncol,1:pver)=0._r8 
meltsdt(1:ncol,1:pver)=0._r8
frzrdt (1:ncol,1:pver)=0._r8
mnuccdo(1:ncol,1:pver)=0._r8

rflx(:,:)=0._r8
sflx(:,:)=0._r8
effc(:,:)=0._r8
effc_fn(:,:)=0._r8
effi(:,:)=0._r8

!water tracers/isotopes:
preo(1:ncol,1:pver)=0._r8
prdso(1:ncol,1:pver)=0._r8
frzro (1:ncol,1:pver)=0._r8
meltso(1:ncol,1:pver)=0._r8
wtfc(1:ncol,1:pver)=0._r8
wtfi(1:ncol,1:pver)=0._r8
wtprelat(1:ncol,1:pver)=0._r8
wtpostlat(1:ncol,1:pver)=0._r8

! initialize multi-level fields
q(1:ncol,1:pver)=qn(1:ncol,1:pver)
t(1:ncol,1:pver)=tn(1:ncol,1:pver)

! initialization
qc(1:ncol,1:top_lev-1) = 0._r8
qi(1:ncol,1:top_lev-1) = 0._r8
nc(1:ncol,1:top_lev-1) = 0._r8
ni(1:ncol,1:top_lev-1) = 0._r8
t1(1:ncol,1:pver) = t(1:ncol,1:pver)
q1(1:ncol,1:pver) = q(1:ncol,1:pver)
qc1(1:ncol,1:pver) = qc(1:ncol,1:pver)
qi1(1:ncol,1:pver) = qi(1:ncol,1:pver)
nc1(1:ncol,1:pver) = nc(1:ncol,1:pver)
ni1(1:ncol,1:pver) = ni(1:ncol,1:pver)

! initialize tendencies to zero
tlat1(1:ncol,1:pver)=0._r8
qvlat1(1:ncol,1:pver)=0._r8
qctend1(1:ncol,1:pver)=0._r8
qitend1(1:ncol,1:pver)=0._r8
nctend1(1:ncol,1:pver)=0._r8
nitend1(1:ncol,1:pver)=0._r8

! initialize precip output
qrout(1:ncol,1:pver)=0._r8
qsout(1:ncol,1:pver)=0._r8
nrout(1:ncol,1:pver)=0._r8
nsout(1:ncol,1:pver)=0._r8
dsout(1:ncol,1:pver)=0._r8

drout(1:ncol,1:pver)=0._r8

reff_rain(1:ncol,1:pver)=0._r8
reff_snow(1:ncol,1:pver)=0._r8

! initialize variables for trop_mozart
nevapr(1:ncol,1:pver) = 0._r8
nevapr2(1:ncol,1:pver) = 0._r8
evapsnow(1:ncol,1:pver) = 0._r8
prain(1:ncol,1:pver) = 0._r8
prodsnow(1:ncol,1:pver) = 0._r8
cmeout(1:ncol,1:pver) = 0._r8

am_evp_st(1:ncol,1:pver) = 0._r8

! for refl calc
rainrt1(1:ncol,1:pver) = 0._r8

! initialize precip fraction and output tendencies
cldmax(1:ncol,1:pver)=mincld

!initialize aerosol number
!        naer2(1:ncol,1:pver,:)=0._r8
dum2l(1:ncol,1:pver)=0._r8
dum2i(1:ncol,1:pver)=0._r8

! initialize avg precip rate
prect1(1:ncol)=0._r8
preci1(1:ncol)=0._r8

else

call micro_mg1_0_init_fields_codon_wrap(ncol, pcols, pver, top_lev, mincld, qn, tn, &
     qc, qi, nc, ni, ncai, ncal, rercld, arcld, pgamrad, lamcrad, deffi, &
     qcsevap, qisevap, qvres, cmeiout, vtrmc, vtrmi, qcsedten, qisedten, &
     prao, prco, mnuccco, mnuccto, msacwio, psacwso, bergso, bergo, melto, &
     homoo, qcreso, prcio, praio, qireso, mnuccro, pracso, meltsdt, frzrdt, &
     mnuccdo, rflx, sflx, effc, effc_fn, effi, preo, prdso, frzro, meltso, &
     wtfc, wtfi, wtprelat, wtpostlat, q, t, t1, q1, qc1, qi1, nc1, ni1, &
     tlat1, qvlat1, qctend1, qitend1, nctend1, nitend1, qrout, qsout, nrout, &
     nsout, dsout, drout, reff_rain, reff_snow, nevapr, nevapr2, evapsnow, &
     prain, prodsnow, cmeout, am_evp_st, rainrt1, cldmax, dum2l, dum2i, &
     prect1, preci1)

end if

! initialize time-varying parameters
! Keep the exponent-heavy expressions in the original Fortran context while
! the surrounding zero/copy initialization is selected between native/Codon.
do k=1,pver
   do i=1,ncol
      rho(i,k)=p(i,k)/(r*t(i,k))
      dv(i,k) = 8.794E-5_r8*t(i,k)**1.81_r8/p(i,k)
      mu(i,k) = 1.496E-6_r8*t(i,k)**1.5_r8/(t(i,k)+120._r8)
      sc(i,k) = mu(i,k)/(rho(i,k)*dv(i,k))
      kap(i,k) = 1.414e3_r8*1.496e-6_r8*t(i,k)**1.5_r8/(t(i,k)+120._r8)

      ! air density adjustment for fallspeed parameters
      ! includes air density correction factor to the
      ! power of 0.54 following Heymsfield and Bansemer 2007

      rhof(i,k)=(rhosu/rho(i,k))**0.54_r8

      arn(i,k)=ar*rhof(i,k)
      asn(i,k)=as*rhof(i,k)
      acn(i,k)=ac*rhof(i,k)
      ain(i,k)=ai*rhof(i,k)

      ! get dz from dp and hydrostatic approx
      ! keep dz positive (define as layer k-1 - layer k)

      dz(i,k)= pdel(i,k)/(rho(i,k)*g)

   end do
end do

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Get humidity and saturation vapor pressures

do k=top_lev,pver

   do i=1,ncol

      ! find wet bulk temperature and saturation value for provisional t and q without
      ! condensation
      
      es(i) = svp_water(t(i,k))
      qs(i) = svp_to_qsat(es(i), p(i,k))

      ! Prevents negative values.
      if (qs(i) < 0.0_r8) then
         qs(i) = 1.0_r8
         es(i) = p(i,k)
      end if

      esl(i,k)=svp_water(t(i,k))
      esi(i,k)=svp_ice(t(i,k))

      ! hm fix, make sure when above freezing that esi=esl, not active yet
      if (t(i,k).gt.tmelt)esi(i,k)=esl(i,k)

      relhum(i,k)=q(i,k)/qs(i)

      ! get cloud fraction, check for minimum

      cldm(i,k)=max(cldn(i,k),mincld)
      cldmw(i,k)=max(cldn(i,k),mincld)

      icldm(i,k)=max(icecldf(i,k),mincld)
      lcldm(i,k)=max(liqcldf(i,k),mincld)

      ! subcolumns, set cloud fraction variables to one
      ! if cloud water or ice is present, if not present
      ! set to mincld (mincld used instead of zero, to prevent
      ! possible division by zero errors

      if (microp_uniform) then

         cldm(i,k)=mincld
         cldmw(i,k)=mincld
         icldm(i,k)=mincld
         lcldm(i,k)=mincld

         if (qc(i,k).ge.qsmall) then
            lcldm(i,k)=1._r8           
            cldm(i,k)=1._r8
            cldmw(i,k)=1._r8
         end if

         if (qi(i,k).ge.qsmall) then             
            cldm(i,k)=1._r8
            icldm(i,k)=1._r8
         end if

      end if               ! sub-columns

      ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)

      nfice(i,k)=0._r8
      dumfice=qc(i,k)+qi(i,k)
      if (dumfice.gt.qsmall .and. qi(i,k).gt.qsmall) then
         nfice(i,k)=qi(i,k)/dumfice
      endif

      if (do_cldice .and. (t(i,k).lt.tmelt - 5._r8)) then

         ! if aerosols interact with ice set number of activated ice nuclei
         dum2=naai(i,k)

         dumnnuc=(dum2-ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
         dumnnuc=max(dumnnuc,0._r8)
         ! get provisional ni and qi after nucleation in order to calculate
         ! Bergeron process below
         ninew=ni(i,k)+dumnnuc*deltat
         qinew=qi(i,k)+dumnnuc*deltat*mi0

         !T>268
      else
         ninew=ni(i,k)
         qinew=qi(i,k)
      end if

      ! Initialize CME components

      cme(i,k) = 0._r8
      cmei(i,k)=0._r8


      !-------------------------------------------------------------------
      !Bergeron process

      ! make sure to initialize bergeron process to zero
      berg(i,k)=0._r8
      prd = 0._r8

      !condensation loop.

      ! get in-cloud qi and ni after nucleation
      if (icldm(i,k) .gt. 0._r8) then 
         qiic(i,k)=qinew/icldm(i,k)
         niic(i,k)=ninew/icldm(i,k)
      else
         qiic(i,k)=0._r8
         niic(i,k)=0._r8
      endif

      !if T < 0 C then bergeron.
      if (do_cldice .and. (t(i,k).lt.273.15_r8)) then

         !if ice exists
         if (qi(i,k).gt.qsmall) then

            bergtsf = 0._r8 ! bergeron time scale (fraction of timestep)

            qvi = svp_to_qsat(esi(i,k), p(i,k))
            qvl = svp_to_qsat(esl(i,k), p(i,k))

            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
            abi = 1._r8+dqsidt*xxls/cpp

            ! get ice size distribution parameters

            if (qiic(i,k).ge.qsmall) then
               lami(k) = (cons1*ci* &
                    niic(i,k)/qiic(i,k))**(1._r8/di)
               n0i(k) = niic(i,k)*lami(k)

               ! check for slope
               ! adjust vars
               if (lami(k).lt.lammini) then

                  lami(k) = lammini
                  n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
               else if (lami(k).gt.lammaxi) then
                  lami(k) = lammaxi
                  n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
               end if

               epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))

               !if liquid exists  
               if (qc(i,k).gt. qsmall) then 

                  !begin bergeron process
                  !     do bergeron (vapor deposition with RHw=1)
                  !     code to find berg (a rate) goes here

                  ! calculate Bergeron process

                  prd = epsi*(qvl-qvi)/abi

               else
                  prd = 0._r8
               end if

               ! multiply by cloud fraction

               prd = prd*min(icldm(i,k),lcldm(i,k))

               !     transfer of existing cloud liquid to ice

               berg(i,k)=max(0._r8,prd)

            end if  !end liquid exists bergeron

            if (berg(i,k).gt.0._r8) then
               bergtsf=max(0._r8,(qc(i,k)/berg(i,k))/deltat) 

               if(bergtsf.lt.1._r8) berg(i,k) = max(0._r8,qc(i,k)/deltat)

            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (bergtsf.lt.1._r8.or.icldm(i,k).gt.lcldm(i,k)) then

               if (qiic(i,k).ge.qsmall) then

                  ! first case is for case when liquid water is present, but is completely depleted 
                  ! in time step, i.e., bergrsf > 0 but < 1

                  if (qc(i,k).ge.qsmall) then
                     rhin  = (1.0_r8 + relhum(i,k)) / 2._r8
                     if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then
                        prd = epsi*(rhin*qvl-qvi)/abi

                        ! multiply by cloud fraction assuming liquid/ice maximum overlap
                        prd = prd*min(icldm(i,k),lcldm(i,k))

                        ! add to cmei
                        cmei(i,k) = cmei(i,k) + (prd * (1._r8- bergtsf))

                     end if ! rhin 
                  end if ! qc > qsmall

                  ! second case is for pure ice cloud, either no liquid, or icldm > lcldm

                  if (qc(i,k).lt.qsmall.or.icldm(i,k).gt.lcldm(i,k)) then

                     ! note: for case of no liquid, need to set liquid cloud fraction to zero
                     ! store liquid cloud fraction in 'dum'

                     if (qc(i,k).lt.qsmall) then 
                        dum=0._r8 
                     else
                        dum=lcldm(i,k)
                     end if

                     ! set RH to grid-mean value for pure ice cloud
                     rhin = relhum(i,k)

                     if ((rhin*esl(i,k)/esi(i,k)) > 1._r8) then

                        prd = epsi*(rhin*qvl-qvi)/abi

                        ! multiply by relevant cloud fraction for pure ice cloud
                        ! assuming maximum overlap of liquid/ice
                        prd = prd*max((icldm(i,k)-dum),0._r8)
                        cmei(i,k) = cmei(i,k) + prd

                     end if ! rhin
                  end if ! qc or icldm > lcldm
               end if ! qiic
            end if ! bergtsf or icldm > lcldm

            !     if deposition, it should not reduce grid mean rhi below 1.0
            if(cmei(i,k) > 0.0_r8 .and. (relhum(i,k)*esl(i,k)/esi(i,k)) > 1._r8 ) &
                 cmei(i,k)=min(cmei(i,k),(q(i,k)-qs(i)*esi(i,k)/esl(i,k))/abi/deltat)

         end if            !end ice exists loop
         !this ends temperature < 0. loop

         !-------------------------------------------------------------------
      end if  ! 
      !..............................................................

      ! evaporation should not exceed available water

      if ((-berg(i,k)).lt.-qc(i,k)/deltat) berg(i,k) = max(qc(i,k)/deltat,0._r8)

      !sublimation process...
      if (do_cldice .and. ((relhum(i,k)*esl(i,k)/esi(i,k)).lt.1._r8 .and. qiic(i,k).ge.qsmall )) then

         qvi = svp_to_qsat(esi(i,k), p(i,k))
         qvl = svp_to_qsat(esl(i,k), p(i,k))
         dqsidt =  xxls*qvi/(rv*t(i,k)**2)
         abi = 1._r8+dqsidt*xxls/cpp

         ! get ice size distribution parameters

         lami(k) = (cons1*ci* &
              niic(i,k)/qiic(i,k))**(1._r8/di)
         n0i(k) = niic(i,k)*lami(k)

         ! check for slope
         ! adjust vars
         if (lami(k).lt.lammini) then

            lami(k) = lammini
            n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
         else if (lami(k).gt.lammaxi) then
            lami(k) = lammaxi
            n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
         end if

         epsi = 2._r8*pi*n0i(k)*rho(i,k)*Dv(i,k)/(lami(k)*lami(k))

         ! modify for ice fraction below
         prd = epsi*(relhum(i,k)*qvl-qvi)/abi * icldm(i,k)
         cmei(i,k)=min(prd,0._r8)

      endif

      ! sublimation should not exceed available ice
      if (cmei(i,k).lt.-qi(i,k)/deltat) cmei(i,k)=-qi(i,k)/deltat

      ! sublimation should not increase grid mean rhi above 1.0 
      if(cmei(i,k) < 0.0_r8 .and. (relhum(i,k)*esl(i,k)/esi(i,k)) < 1._r8 ) &
           cmei(i,k)=min(0._r8,max(cmei(i,k),(q(i,k)-qs(i)*esi(i,k)/esl(i,k))/abi/deltat))

      ! limit cmei due for roundoff error

      cmei(i,k)=cmei(i,k)*omsm

      ! conditional for ice nucleation 
      if (do_cldice .and. (t(i,k).lt.(tmelt - 5._r8))) then 

         ! using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
         ! ice nucleation rate (dum2) has already been calculated and read in (naai)

         dum2i(i,k)=naai(i,k)
      else
         dum2i(i,k)=0._r8
      end if

   end do ! i loop
end do ! k loop


call micro_mg1_0_select_colzero_impl()
if (micro_mg1_0_colzero_use_native_impl) then
   !! initialize sub-step precip flux variables
   do i=1,ncol
      !! flux is zero at top interface, so these should stay as 0.
      rflx1(i,1)=0._r8
      sflx1(i,1)=0._r8
      do k=top_lev,pver

         ! initialize normal and sub-step precip flux variables
         rflx1(i,k+1)=0._r8
         sflx1(i,k+1)=0._r8
      end do ! i loop
   end do ! k loop
   !! initialize final precip flux variables.
   do i=1,ncol
      !! flux is zero at top interface, so these should stay as 0.
      rflx(i,1)=0._r8
      sflx(i,1)=0._r8
      do k=top_lev,pver
         ! initialize normal and sub-step precip flux variables
         rflx(i,k+1)=0._r8
         sflx(i,k+1)=0._r8
      end do ! i loop
   end do ! k loop

   do i=1,ncol
      ltrue(i)=0
      do k=top_lev,pver
         ! skip microphysical calculations if no cloud water

         if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall.or.cmei(i,k).ge.qsmall) ltrue(i)=1
      end do
   end do
else
   call micro_mg1_0_flux_ltrue_init_codon_wrap(ncol, pcols, pver, top_lev, qsmall, &
        rflx1, sflx1, rflx, sflx, qc, qi, cmei, ltrue)
end if

! assign number of sub-steps to iter
! use 2 sub-steps, following tests described in MG2008
iter = 2

! get sub-step time step
deltat=deltat/real(iter)

! since activation/nucleation processes are fast, need to take into account
! factor mtime = mixing timescale in cloud / model time step
! mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
! for now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk

!        note: mtime for bulk aerosols was set to: mtime=deltat/1200._r8

mtime=1._r8
if (micro_mg1_0_colzero_use_native_impl) then
   rate1ord_cw2pr_st(:,:)=0._r8 ! rce 2010/05/01
else
   call micro_mg1_0_rate1ord_zero_codon_wrap(ncol, pcols, pver, rate1ord_cw2pr_st)
end if
!!!! skip calculations if no cloud water
do i=1,ncol
   if (ltrue(i).eq.0) then
      if (micro_mg1_0_colzero_use_native_impl) then
      tlat(i,1:pver)=0._r8
      qvlat(i,1:pver)=0._r8
      qctend(i,1:pver)=0._r8
      qitend(i,1:pver)=0._r8
      qnitend(i,1:pver)=0._r8
      qrtend(i,1:pver)=0._r8
      nctend(i,1:pver)=0._r8
      nitend(i,1:pver)=0._r8
      nrtend(i,1:pver)=0._r8
      nstend(i,1:pver)=0._r8
      prect(i)=0._r8
      preci(i)=0._r8
      qniic(i,1:pver)=0._r8
      qric(i,1:pver)=0._r8
      nsic(i,1:pver)=0._r8
      nric(i,1:pver)=0._r8
      rainrt(i,1:pver)=0._r8
      !water tracers:
      !-------------
      qrtend_copy(i,1:pver)=0._r8
      qnitend_copy(i,1:pver)=0._r8
      !-------------
      else
         call micro_mg1_0_no_cloud_zero_column_codon_wrap(i, pcols, pver, tlat, qvlat, qctend, qitend, &
              qnitend, qrtend, nctend, nitend, nrtend, nstend, prect, preci, qniic, qric, nsic, nric, &
              rainrt, qrtend_copy, qnitend_copy)
      end if
      goto 300
   end if

   qcsinksum_rate1ord(1:pver)=0._r8 
   qcsum_rate1ord(1:pver)=0._r8 


!!!!!!!!! begin sub-step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !.....................................................................................................
   do it=1,iter

      ! initialize sub-step microphysical tendencies

      if (micro_mg1_0_colzero_use_native_impl) then
      tlat(i,1:pver)=0._r8
      qvlat(i,1:pver)=0._r8
      qctend(i,1:pver)=0._r8
      qitend(i,1:pver)=0._r8
      qnitend(i,1:pver)=0._r8
      qrtend(i,1:pver)=0._r8
      nctend(i,1:pver)=0._r8
      nitend(i,1:pver)=0._r8
      nrtend(i,1:pver)=0._r8
      nstend(i,1:pver)=0._r8

      !water tracers:
      qrtend_copy(i,1:pver)=0._r8
      qnitend_copy(i,1:pver)=0._r8

      ! initialize diagnostic precipitation to zero

      qniic(i,1:pver)=0._r8
      qric(i,1:pver)=0._r8
      nsic(i,1:pver)=0._r8
      nric(i,1:pver)=0._r8

      rainrt(i,1:pver)=0._r8

      else
         call micro_mg1_0_substep_zero_column_codon_wrap(i, pcols, pver, tlat, qvlat, qctend, qitend, &
              qnitend, qrtend, nctend, nitend, nrtend, nstend, qniic, qric, nsic, nric, rainrt, &
              qrtend_copy, qnitend_copy)
      end if

      ! begin new i,k loop, calculate new cldmax after adjustment to cldm above

      ! initialize vertically-integrated rain and snow tendencies

      qrtot = 0._r8
      nrtot = 0._r8
      qstot = 0._r8
      nstot = 0._r8

      ! initialize precip at surface

      prect(i)=0._r8
      preci(i)=0._r8

      if (.not. micro_mg1_0_colzero_use_native_impl) then
         call micro_mg1_0_substep_setup_column_codon_wrap(i, pcols, pver, top_lev, qsmall, mincld, &
              qc, qi, ni, cldm, cmei, cwml, cwmi, ums, uns, umr, unr, nsubi, nsubc)
      end if

      do k=top_lev,pver
      
         qcvar=relvar(i,k)
         cons2=gamma(qcvar+2.47_r8)
         cons3=gamma(qcvar)
         cons9=gamma(qcvar+2._r8)
         cons10=gamma(qcvar+1._r8)
         cons12=gamma(qcvar+1.15_r8) 
         cons15=gamma(qcvar+bc/3._r8)
         cons18=qcvar**2.47_r8
         cons19=qcvar**2
         cons20=qcvar**1.15_r8

         ! set cwml and cwmi to current qc and qi

         if (micro_mg1_0_colzero_use_native_impl) then
            cwml(i,k)=qc(i,k)
            cwmi(i,k)=qi(i,k)

            ! initialize precip fallspeeds to zero

            ums(k)=0._r8
            uns(k)=0._r8
            umr(k)=0._r8
            unr(k)=0._r8
         end if

         ! calculate precip fraction based on maximum overlap assumption

         ! for sub-columns cldm has already been set to 1 if cloud
         ! water or ice is present, so cldmax will be correctly set below
         ! and nothing extra needs to be done here

         if (k.eq.top_lev) then
            cldmax(i,k)=cldm(i,k)
         else
            ! if rain or snow mix ratio is smaller than
            ! threshold, then set cldmax to cloud fraction at current level

            if (do_clubb_sgs) then
               if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall) then
                  cldmax(i,k)=cldm(i,k)
               else
                  cldmax(i,k)=cldmax(i,k-1)
               end if
            else

               if (qric(i,k-1).ge.qsmall.or.qniic(i,k-1).ge.qsmall) then
                  cldmax(i,k)=max(cldmax(i,k-1),cldm(i,k))
               else
                  cldmax(i,k)=cldm(i,k)
               end if
            endif
         end if

         ! decrease in number concentration due to sublimation/evap
         ! divide by cloud fraction to get in-cloud decrease
         ! don't reduce Nc due to bergeron process

         if (micro_mg1_0_colzero_use_native_impl) then
            if (cmei(i,k) < 0._r8 .and. qi(i,k) > qsmall .and. cldm(i,k) > mincld) then
               nsubi(k)=cmei(i,k)/qi(i,k)*ni(i,k)/cldm(i,k)
            else
               nsubi(k)=0._r8
            end if
            nsubc(k)=0._r8
         end if


         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%

         if (do_cldice .and. dum2i(i,k).gt.0._r8.and.t(i,k).lt.(tmelt - 5._r8).and. &
              relhum(i,k)*esl(i,k)/esi(i,k).gt. rhmini+0.05_r8) then

            !if NCAI > 0. then set numice = ncai (as before)
            !note: this is gridbox averaged

            nnuccd(k)=(dum2i(i,k)-ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
            nnuccd(k)=max(nnuccd(k),0._r8)
            nimax = dum2i(i,k)*icldm(i,k)

            !Calc mass of new particles using new crystal mass...
            !also this will be multiplied by mtime as nnuccd is...

            mnuccd(k) = nnuccd(k) * mi0

            !  add mnuccd to cmei....
            cmei(i,k)= cmei(i,k) + mnuccd(k) * mtime

            !  limit cmei

            qvi = svp_to_qsat(esi(i,k), p(i,k))
            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
            abi = 1._r8+dqsidt*xxls/cpp
            cmei(i,k)=min(cmei(i,k),(q(i,k)-qvi)/abi/deltat)

            ! limit for roundoff error
            cmei(i,k)=cmei(i,k)*omsm

         else
            nnuccd(k)=0._r8
            nimax = 0._r8
            mnuccd(k) = 0._r8
         end if

         !c............................................................................
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
         ! for microphysical process calculations
         ! units are kg/kg for mixing ratio, 1/kg for number conc

         ! limit in-cloud values to 0.005 kg/kg

         if (micro_mg1_0_colzero_use_native_impl) then
            qcic(i,k)=min(cwml(i,k)/lcldm(i,k),5.e-3_r8)
            qiic(i,k)=min(cwmi(i,k)/icldm(i,k),5.e-3_r8)
            ncic(i,k)=max(nc(i,k)/lcldm(i,k),0._r8)
            niic(i,k)=max(ni(i,k)/icldm(i,k),0._r8)

            if (qc(i,k) - berg(i,k)*deltat.lt.qsmall) then
               qcic(i,k)=0._r8
               ncic(i,k)=0._r8
               if (qc(i,k)-berg(i,k)*deltat.lt.0._r8) then
                  berg(i,k)=qc(i,k)/deltat*omsm
               end if
            end if

            if (do_cldice .and. qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.qsmall) then
               qiic(i,k)=0._r8
               niic(i,k)=0._r8
               if (qi(i,k)+(cmei(i,k)+berg(i,k))*deltat.lt.0._r8) then
                  cmei(i,k)=(-qi(i,k)/deltat-berg(i,k))*omsm
               end if
            end if

            ! add to cme output

            cmeout(i,k) = cmeout(i,k)+cmei(i,k)

            !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ! droplet activation
            ! calculate potential for droplet activation if cloud water is present
            ! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
            ! number tendency (npccnin) is read in from companion routine

            ! assume aerosols already activated are equal to number of existing droplets for simplicity
            ! multiply by cloud fraction to obtain grid-average tendency

            if (qcic(i,k).ge.qsmall) then
               npccn(k) = max(0._r8,npccnin(i,k))
               dum2l(i,k)=(nc(i,k)+npccn(k)*deltat)/lcldm(i,k)
               dum2l(i,k)=max(dum2l(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3
               ncmax = dum2l(i,k)*lcldm(i,k)
            else
               npccn(k)=0._r8
               dum2l(i,k)=0._r8
               ncmax = 0._r8
            end if
         else
            ncmax = micro_mg1_0_incloud_activation_prep_codon_wrap(i, k, pcols, pver, deltat, qsmall, &
                 omsm, cdnl, do_cldice, cwml, cwmi, lcldm, icldm, nc, ni, qc, qi, berg, cmei, cmeout, &
                 npccnin, rho, qcic, qiic, ncic, niic, npccn, dum2l)
         end if

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! get size distribution parameters based on in-cloud cloud water/ice 
         ! these calculations also ensure consistency between number and mixing ratio
         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         !......................................................................
         ! cloud ice

         if (qiic(i,k).ge.qsmall) then

            ! add upper limit to in-cloud number concentration to prevent numerical error
            niic(i,k)=min(niic(i,k),qiic(i,k)*1.e20_r8)

            lami(k) = (cons1*ci*niic(i,k)/qiic(i,k))**(1._r8/di)
            n0i(k) = niic(i,k)*lami(k)

            ! check for slope
            ! adjust vars

            if (lami(k).lt.lammini) then

               lami(k) = lammini
               n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
               niic(i,k) = n0i(k)/lami(k)
            else if (lami(k).gt.lammaxi) then
               lami(k) = lammaxi
               n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*cons1)
               niic(i,k) = n0i(k)/lami(k)
            end if

         else
            lami(k) = 0._r8
            n0i(k) = 0._r8
         end if

         if (qcic(i,k).ge.qsmall) then

            ! add upper limit to in-cloud number concentration to prevent numerical error
            ncic(i,k)=min(ncic(i,k),qcic(i,k)*1.e20_r8)

            ncic(i,k)=max(ncic(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm  

            ! get pgam from fit to observations of martin et al. 1994

            pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
            pgam(k)=1._r8/(pgam(k)**2)-1._r8
            pgam(k)=max(pgam(k),2._r8)
            pgam(k)=min(pgam(k),15._r8)

            ! calculate lamc

            lamc(k) = (pi/6._r8*rhow*ncic(i,k)*gamma(pgam(k)+4._r8)/ &
                 (qcic(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)

            ! lammin, 50 micron diameter max mean size

            lammin = (pgam(k)+1._r8)/50.e-6_r8
            lammax = (pgam(k)+1._r8)/2.e-6_r8

            if (lamc(k).lt.lammin) then
               lamc(k) = lammin
               ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
            else if (lamc(k).gt.lammax) then
               lamc(k) = lammax
               ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
            end if

            ! parameter to calculate droplet freezing

            cdist1(k) = ncic(i,k)/gamma(pgam(k)+1._r8) 

         else
            lamc(k) = 0._r8
            cdist1(k) = 0._r8
         end if

         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! begin micropysical process calculations 
         !.................................................................
         ! autoconversion of cloud liquid water to rain
         ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
         ! minimum qc of 1 x 10^-8 prevents floating point error

         if (qcic(i,k).ge.1.e-8_r8) then

            ! nprc is increase in rain number conc due to autoconversion
            ! nprc1 is decrease in cloud droplet conc due to autoconversion

            ! assume exponential sub-grid distribution of qc, resulting in additional
            ! factor related to qcvar below

            ! hm switch for sub-columns, don't include sub-grid qc
            if (microp_uniform) then

               prc(k) = 1350._r8*qcic(i,k)**2.47_r8* &
                    (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
               nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
               nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))

            else

               prc(k) = cons2/(cons3*cons18)*1350._r8*qcic(i,k)**2.47_r8* &
                    (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
               nprc(k) = prc(k)/cons22
               nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))

            end if               ! sub-column switch

         else
            prc(k)=0._r8
            nprc(k)=0._r8
            nprc1(k)=0._r8
         end if

         ! add autoconversion to precip from above to get provisional rain mixing ratio
         ! and number concentration (qric and nric)

         ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)

         dum=0.45_r8
         dum1=0.45_r8

         if (k.eq.top_lev) then
            qric(i,k)=prc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
            nric(i,k)=nprc(k)*lcldm(i,k)*dz(i,k)/cldmax(i,k)/dum
         else
            if (qric(i,k-1).ge.qsmall) then
               dum=umr(k-1)
               dum1=unr(k-1)
            end if

            ! no autoconversion of rain number if rain/snow falling from above
            ! this assumes that new drizzle drops formed by autoconversion are rapidly collected
            ! by the existing rain/snow particles from above

            if (qric(i,k-1).ge.1.e-9_r8.or.qniic(i,k-1).ge.1.e-9_r8) then
               nprc(k)=0._r8
            end if

            qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
                 (rho(i,k)*dz(i,k)*((pra(k-1)+prc(k))*lcldm(i,k)+(pre(k-1)-pracs(k-1)-mnuccr(k-1))*cldmax(i,k))))&
                 /(dum*rho(i,k)*cldmax(i,k))
            nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
                 (rho(i,k)*dz(i,k)*(nprc(k)*lcldm(i,k)+(nsubr(k-1)-npracs(k-1)-nnuccr(k-1)+nragg(k-1))*cldmax(i,k))))&
                 /(dum1*rho(i,k)*cldmax(i,k))

         end if

         !.......................................................................
         ! Autoconversion of cloud ice to snow
         ! similar to Ferrier (1994)

         if (do_cldice) then
            if (t(i,k).le.273.15_r8.and.qiic(i,k).ge.qsmall) then

               ! note: assumes autoconversion timescale of 180 sec
               
               nprci(k) = n0i(k)/(lami(k)*180._r8)*exp(-lami(k)*dcs)

               prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
                    (cons23/lami(k)+3._r8*cons24/lami(k)**2+ &
                    6._r8*dcs/lami(k)**3+6._r8/lami(k)**4)*exp(-lami(k)*dcs)
            else
               prci(k)=0._r8
               nprci(k)=0._r8
            end if
         else
            ! Add in the particles that we have already converted to snow, and
            ! don't do any further autoconversion of ice.
            prci(k)  = tnd_qsnow(i, k) / cldm(i,k)
            nprci(k) = tnd_nsnow(i, k) / cldm(i,k)
         end if

         ! add autoconversion to flux from level above to get provisional snow mixing ratio
         ! and number concentration (qniic and nsic)

         dum=(asn(i,k)*cons25)
         dum1=(asn(i,k)*cons25)

         if (k.eq.top_lev) then
            qniic(i,k)=prci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
            nsic(i,k)=nprci(k)*icldm(i,k)*dz(i,k)/cldmax(i,k)/dum
         else
            if (qniic(i,k-1).ge.qsmall) then
               dum=ums(k-1)
               dum1=uns(k-1)
            end if

            qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
                 (rho(i,k)*dz(i,k)*((prci(k)+prai(k-1)+psacws(k-1)+bergs(k-1))*icldm(i,k)+(prds(k-1)+ &
                 pracs(k-1)+mnuccr(k-1))*cldmax(i,k))))&
                 /(dum*rho(i,k)*cldmax(i,k))

            nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
                 (rho(i,k)*dz(i,k)*(nprci(k)*icldm(i,k)+(nsubs(k-1)+nsagg(k-1)+nnuccr(k-1))*cldmax(i,k))))&
                 /(dum1*rho(i,k)*cldmax(i,k))

         end if

         ! if precip mix ratio is zero so should number concentration

         if (qniic(i,k).lt.qsmall) then
            qniic(i,k)=0._r8
            nsic(i,k)=0._r8
         end if

         if (qric(i,k).lt.qsmall) then
            qric(i,k)=0._r8
            nric(i,k)=0._r8
         end if

         ! make sure number concentration is a positive number to avoid 
         ! taking root of negative later

         nric(i,k)=max(nric(i,k),0._r8)
         nsic(i,k)=max(nsic(i,k),0._r8)

         !.......................................................................
         ! get size distribution parameters for precip
         !......................................................................
         ! rain

         if (qric(i,k).ge.qsmall) then
            lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
            n0r(k) = nric(i,k)*lamr(k)

            ! check for slope
            ! adjust vars

            if (lamr(k).lt.lamminr) then

               lamr(k) = lamminr

               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
            else if (lamr(k).gt.lammaxr) then
               lamr(k) = lammaxr
               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
            end if

            ! provisional rain number and mass weighted mean fallspeed (m/s)

            unr(k) = min(arn(i,k)*cons4/lamr(k)**br,9.1_r8*rhof(i,k))
            umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(i,k))

         else
            lamr(k) = 0._r8
            n0r(k) = 0._r8
            umr(k) = 0._r8
            unr(k) = 0._r8
         end if

         !......................................................................
         ! snow

         if (qniic(i,k).ge.qsmall) then
            lams(k) = (cons6*cs*nsic(i,k)/qniic(i,k))**(1._r8/ds)
            n0s(k) = nsic(i,k)*lams(k)

            ! check for slope
            ! adjust vars

            if (lams(k).lt.lammins) then
               lams(k) = lammins
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
               nsic(i,k) = n0s(k)/lams(k)

            else if (lams(k).gt.lammaxs) then
               lams(k) = lammaxs
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
               nsic(i,k) = n0s(k)/lams(k)
            end if

            ! provisional snow number and mass weighted mean fallspeed (m/s)

            ums(k) = min(asn(i,k)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof(i,k))
            uns(k) = min(asn(i,k)*cons7/lams(k)**bs,1.2_r8*rhof(i,k))

         else
            lams(k) = 0._r8
            n0s(k) = 0._r8
            ums(k) = 0._r8
            uns(k) = 0._r8
         end if

         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         ! heterogeneous freezing of cloud water

         if (.not. use_hetfrz_classnuc) then

            if (do_cldice .and. qcic(i,k).ge.qsmall .and. t(i,k).lt.269.15_r8) then

               ! immersion freezing (Bigg, 1953)


               ! subcolumns

               if (microp_uniform) then

                  mnuccc(k) = &
                     pi*pi/36._r8*rhow* &
                     cdist1(k)*gamma(7._r8+pgam(k))* &
                     bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
                     lamc(k)**3/lamc(k)**3

                  nnuccc(k) = &
                     pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
                     *bimm* &
                     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3

               else

                  mnuccc(k) = cons9/(cons3*cons19)* &
                     pi*pi/36._r8*rhow* &
                     cdist1(k)*gamma(7._r8+pgam(k))* &
                     bimm*(exp(aimm*(273.15_r8-t(i,k)))-1._r8)/ &
                     lamc(k)**3/lamc(k)**3

                  nnuccc(k) = cons10/(cons3*qcvar)* &
                     pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
                     *bimm* &
                     (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamc(k)**3
               end if           ! sub-columns


               ! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
               ! dust size and number in 4 bins are read in from companion routine

               tcnt=(270.16_r8-t(i,k))**1.3_r8
               viscosity=1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
               mfp=2.0_r8*viscosity/(p(i,k)  &                   ! Mean free path (m)
                  *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i,k))))           

               nslip1=1.0_r8+(mfp/rndst(i,k,1))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,1)/mfp))))! Slip correction factor
               nslip2=1.0_r8+(mfp/rndst(i,k,2))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,2)/mfp))))
               nslip3=1.0_r8+(mfp/rndst(i,k,3))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,3)/mfp))))
               nslip4=1.0_r8+(mfp/rndst(i,k,4))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(i,k,4)/mfp))))

               ndfaer1=1.381e-23_r8*t(i,k)*nslip1/(6._r8*pi*viscosity*rndst(i,k,1))  ! aerosol diffusivity (m2/s)
               ndfaer2=1.381e-23_r8*t(i,k)*nslip2/(6._r8*pi*viscosity*rndst(i,k,2))
               ndfaer3=1.381e-23_r8*t(i,k)*nslip3/(6._r8*pi*viscosity*rndst(i,k,3))
               ndfaer4=1.381e-23_r8*t(i,k)*nslip4/(6._r8*pi*viscosity*rndst(i,k,4))


               if (microp_uniform) then

                  mnucct(k) = &
                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
                     cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4

                  nnucct(k) = (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
                     cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)

               else

                  mnucct(k) = gamma(qcvar+4._r8/3._r8)/(cons3*qcvar**(4._r8/3._r8))*  &
                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*pi*pi/3._r8*rhow* &
                     cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4

                  nnucct(k) =  gamma(qcvar+1._r8/3._r8)/(cons3*qcvar**(1._r8/3._r8))*  &
                     (ndfaer1*(nacon(i,k,1)*tcnt)+ndfaer2*(nacon(i,k,2)*tcnt)+ &
                     ndfaer3*(nacon(i,k,3)*tcnt)+ndfaer4*(nacon(i,k,4)*tcnt))*2._r8*pi*  &
                     cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)

               end if      ! sub-column switch

               ! make sure number of droplets frozen does not exceed available ice nuclei concentration
               ! this prevents 'runaway' droplet freezing

               if (nnuccc(k)*lcldm(i,k).gt.nnuccd(k)) then
                  dum=(nnuccd(k)/(nnuccc(k)*lcldm(i,k)))
                  ! scale mixing ratio of droplet freezing with limit
                  mnuccc(k)=mnuccc(k)*dum
                  nnuccc(k)=nnuccd(k)/lcldm(i,k)
               end if

            else
               mnuccc(k)=0._r8
               nnuccc(k)=0._r8
               mnucct(k)=0._r8
               nnucct(k)=0._r8
            end if

         else
            if (do_cldice .and. qcic(i,k) >= qsmall) then
               con1 = 1._r8/(1.333_r8*pi)**0.333_r8
               r3lx = con1*(rho(i,k)*qcic(i,k)/(rhow*max(ncic(i,k)*rho(i,k), 1.0e6_r8)))**0.333_r8 ! in m
               r3lx = max(4.e-6_r8, r3lx)
               mi0l = 4._r8/3._r8*pi*rhow*r3lx**3_r8
                
               nnuccc(k) = frzimm(i,k)*1.0e6_r8/rho(i,k)
               mnuccc(k) = nnuccc(k)*mi0l 

               nnucct(k) = frzcnt(i,k)*1.0e6_r8/rho(i,k)
               mnucct(k) = nnucct(k)*mi0l 

               nnudep(k) = frzdep(i,k)*1.0e6_r8/rho(i,k)
               mnudep(k) = nnudep(k)*mi0
            else
               nnuccc(k) = 0._r8
               mnuccc(k) = 0._r8

               nnucct(k) = 0._r8
               mnucct(k) = 0._r8

               nnudep(k) = 0._r8
               mnudep(k) = 0._r8
            end if
         endif


         !.......................................................................
         ! snow self-aggregation from passarelli, 1978, used by reisner, 1998
         ! this is hard-wired for bs = 0.4 for now
         ! ignore self-collection of cloud ice

         if (qniic(i,k).ge.qsmall .and. t(i,k).le.273.15_r8) then
            nsagg(k) = -1108._r8*asn(i,k)*Eii* &
                 pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)*rho(i,k)** &
                 ((2._r8+bs)/3._r8)*qniic(i,k)**((2._r8+bs)/3._r8)* &
                 (nsic(i,k)*rho(i,k))**((4._r8-bs)/3._r8)/ &
                 (4._r8*720._r8*rho(i,k))
         else
            nsagg(k)=0._r8
         end if

         !.......................................................................
         ! accretion of cloud droplets onto snow/graupel
         ! here use continuous collection equation with
         ! simple gravitational collection kernel
         ! ignore collisions between droplets/cloud ice
         ! since minimum size ice particle for accretion is 50 - 150 micron

         ! ignore collision of snow with droplets above freezing

         if (qniic(i,k).ge.qsmall .and. t(i,k).le.tmelt .and. &
              qcic(i,k).ge.qsmall) then

            ! put in size dependent collection efficiency
            ! mean diameter of snow is area-weighted, since
            ! accretion is function of crystal geometric area
            ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)

            dc0 = (pgam(k)+1._r8)/lamc(k)
            ds0 = 1._r8/lams(k)
            dum = dc0*dc0*uns(k)*rhow/(9._r8*mu(i,k)*ds0)
            eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))

            eci = max(eci,0._r8)
            eci = min(eci,1._r8)


            ! no impact of sub-grid distribution of qc since psacws
            ! is linear in qc

            psacws(k) = pi/4._r8*asn(i,k)*qcic(i,k)*rho(i,k)* &
                 n0s(k)*Eci*cons11/ &
                 lams(k)**(bs+3._r8)
            npsacws(k) = pi/4._r8*asn(i,k)*ncic(i,k)*rho(i,k)* &
                 n0s(k)*Eci*cons11/ &
                 lams(k)**(bs+3._r8)
         else
            psacws(k)=0._r8
            npsacws(k)=0._r8
         end if

         ! add secondary ice production due to accretion of droplets by snow 
         ! (Hallet-Mossop process) (from Cotton et al., 1986)

         if (.not. do_cldice) then
            ni_secp   = 0.0_r8
            nsacwi(k) = 0.0_r8
            msacwi(k) = 0.0_r8
         else if((t(i,k).lt.270.16_r8) .and. (t(i,k).ge.268.16_r8)) then
            ni_secp   = 3.5e8_r8*(270.16_r8-t(i,k))/2.0_r8*psacws(k)
            nsacwi(k) = ni_secp
            msacwi(k) = min(ni_secp*mi0,psacws(k))
         else if((t(i,k).lt.268.16_r8) .and. (t(i,k).ge.265.16_r8)) then
            ni_secp   = 3.5e8_r8*(t(i,k)-265.16_r8)/3.0_r8*psacws(k)
            nsacwi(k) = ni_secp
            msacwi(k) = min(ni_secp*mi0,psacws(k))
         else
            ni_secp   = 0.0_r8
            nsacwi(k) = 0.0_r8
            msacwi(k) = 0.0_r8
         endif
         psacws(k) = max(0.0_r8,psacws(k)-ni_secp*mi0)

         !.......................................................................
         ! accretion of rain water by snow
         ! formula from ikawa and saito, 1991, used by reisner et al., 1998

         if (qric(i,k).ge.1.e-8_r8 .and. qniic(i,k).ge.1.e-8_r8 .and. & 
              t(i,k).le.273.15_r8) then

            pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k)-0.95_r8*ums(k))**2+ &
                 0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho(i,k)* &
                 n0r(k)*n0s(k)* &
                 (5._r8/(lamr(k)**6*lams(k))+ &
                 2._r8/(lamr(k)**5*lams(k)**2)+ &
                 0.5_r8/(lamr(k)**4*lams(k)**3)))

            npracs(k) = pi/2._r8*rho(i,k)*ecr*(1.7_r8*(unr(k)-uns(k))**2+ &
                 0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
                 (1._r8/(lamr(k)**3*lams(k))+ &
                 1._r8/(lamr(k)**2*lams(k)**2)+ &
                 1._r8/(lamr(k)*lams(k)**3))

         else
            pracs(k)=0._r8
            npracs(k)=0._r8
         end if

         !.......................................................................
         ! heterogeneous freezing of rain drops
         ! follows from Bigg (1953)

         if (t(i,k).lt.269.15_r8 .and. qric(i,k).ge.qsmall) then

            mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
                 (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamr(k)**3 &
                 /lamr(k)**3

            nnuccr(k) = pi*nric(i,k)*bimm* &
                 (exp(aimm*(273.15_r8-t(i,k)))-1._r8)/lamr(k)**3
         else
            mnuccr(k)=0._r8
            nnuccr(k)=0._r8
         end if

         !.......................................................................
         ! accretion of cloud liquid water by rain
         ! formula from Khrouditnov and Kogan (2000)
         ! gravitational collection kernel, droplet fall speed neglected

         if (qric(i,k).ge.qsmall .and. qcic(i,k).ge.qsmall) then

            ! include sub-grid distribution of cloud water

            ! add sub-column switch

            if (microp_uniform) then

               pra(k) = 67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
               npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))

            else

               pra(k) = accre_enhan(i,k)*(cons12/(cons3*cons20)*67._r8*(qcic(i,k)*qric(i,k))**1.15_r8)
               npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))

            end if               ! sub-column switch

         else
            pra(k)=0._r8
            npra(k)=0._r8
         end if

         !.......................................................................
         ! Self-collection of rain drops
         ! from Beheng(1994)

         if (qric(i,k).ge.qsmall) then
            nragg(k) = -8._r8*nric(i,k)*qric(i,k)*rho(i,k)
         else
            nragg(k)=0._r8
         end if

         !.......................................................................
         ! Accretion of cloud ice by snow
         ! For this calculation, it is assumed that the Vs >> Vi
         ! and Ds >> Di for continuous collection

         if (do_cldice .and. qniic(i,k).ge.qsmall.and.qiic(i,k).ge.qsmall &
              .and.t(i,k).le.273.15_r8) then

            prai(k) = pi/4._r8*asn(i,k)*qiic(i,k)*rho(i,k)* &
                 n0s(k)*Eii*cons11/ &
                 lams(k)**(bs+3._r8)
            nprai(k) = pi/4._r8*asn(i,k)*niic(i,k)* &
                 rho(i,k)*n0s(k)*Eii*cons11/ &
                 lams(k)**(bs+3._r8)
         else
            prai(k)=0._r8
            nprai(k)=0._r8
         end if

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! calculate evaporation/sublimation of rain and snow
         ! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
         ! in-cloud condensation/deposition of rain and snow is neglected
         ! except for transfer of cloud water to snow through bergeron process

         ! initialize evap/sub tendncies
         pre(k)=0._r8
         prds(k)=0._r8

         ! evaporation of rain
         ! only calculate if there is some precip fraction > cloud fraction

         if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8.or.cldmax(i,k).gt.lcldm(i,k)) then

            ! set temporary cloud fraction to zero if cloud water + ice is very small
            ! this will ensure that evaporation/sublimation of precip occurs over
            ! entire grid cell, since min cloud fraction is specified otherwise
            if (qcic(i,k)+qiic(i,k).lt.1.e-6_r8) then
               dum=0._r8
            else
               dum=lcldm(i,k)
            end if

            ! saturation vapor pressure
            esn=svp_water(t(i,k))
            qsn=svp_to_qsat(esn, p(i,k))

            ! recalculate saturation vapor pressure for liquid and ice
            esl(i,k)=esn
            esi(i,k)=svp_ice(t(i,k))
            ! hm fix, make sure when above freezing that esi=esl, not active yet
            if (t(i,k).gt.tmelt)esi(i,k)=esl(i,k)

            ! calculate q for out-of-cloud region
            qclr=(q(i,k)-dum*qsn)/(1._r8-dum)

            if (qric(i,k).ge.qsmall) then

               qvs=svp_to_qsat(esl(i,k), p(i,k))
               dqsdt = xxlv*qvs/(rv*t(i,k)**2)
               ab = 1._r8+dqsdt*xxlv/cpp
               epsr = 2._r8*pi*n0r(k)*rho(i,k)*Dv(i,k)* &
                    (f1r/(lamr(k)*lamr(k))+ &
                    f2r*(arn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                    sc(i,k)**(1._r8/3._r8)*cons13/ &
                    (lamr(k)**(5._r8/2._r8+br/2._r8)))

               pre(k) = epsr*(qclr-qvs)/ab

               ! only evaporate in out-of-cloud region
               ! and distribute across cldmax
               pre(k)=min(pre(k)*(cldmax(i,k)-dum),0._r8)
               pre(k)=pre(k)/cldmax(i,k)
               am_evp_st(i,k) = max(cldmax(i,k)-dum, 0._r8)
            end if

            ! sublimation of snow
            if (qniic(i,k).ge.qsmall) then
               qvi=svp_to_qsat(esi(i,k), p(i,k))
               dqsidt =  xxls*qvi/(rv*t(i,k)**2)
               abi = 1._r8+dqsidt*xxls/cpp
               epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                    (f1s/(lams(k)*lams(k))+ &
                    f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                    sc(i,k)**(1._r8/3._r8)*cons14/ &
                    (lams(k)**(5._r8/2._r8+bs/2._r8)))
               prds(k) = epss*(qclr-qvi)/abi

               ! only sublimate in out-of-cloud region and distribute over cldmax
               prds(k)=min(prds(k)*(cldmax(i,k)-dum),0._r8)
               prds(k)=prds(k)/cldmax(i,k)
               am_evp_st(i,k) = max(cldmax(i,k)-dum, 0._r8)
            end if

            ! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
            ! get updated RH at end of time step based on cloud water/ice condensation/evap

            qtmp=q(i,k)-(cmei(i,k)+(pre(k)+prds(k))*cldmax(i,k))*deltat
            ttmp=t(i,k)+((pre(k)*cldmax(i,k))*xxlv+ &
                 (cmei(i,k)+prds(k)*cldmax(i,k))*xxls)*deltat/cpp

            !limit range of temperatures!
            ttmp=max(180._r8,min(ttmp,323._r8))

            esn=svp_water(ttmp)  ! use rhw to allow ice supersaturation
            qsn=svp_to_qsat(esn, p(i,k))

            ! modify precip evaporation rate if q > qsat
            if (qtmp.gt.qsn) then
               if (pre(k)+prds(k).lt.-1.e-20_r8) then
                  dum1=pre(k)/(pre(k)+prds(k))
                  ! recalculate q and t after cloud water cond but without precip evap
                  qtmp=q(i,k)-(cmei(i,k))*deltat
                  ttmp=t(i,k)+(cmei(i,k)*xxls)*deltat/cpp
                  esn=svp_water(ttmp) ! use rhw to allow ice supersaturation
                  qsn=svp_to_qsat(esn, p(i,k))
                  dum=(qtmp-qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp**2))
                  dum=min(dum,0._r8)

                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
                  pre(k)=dum*dum1/deltat/cldmax(i,k)

                  ! do separately using RHI for prds....
                  esn=svp_ice(ttmp) ! use rhi to allow ice supersaturation
                  qsn=svp_to_qsat(esn, p(i,k))
                  dum=(qtmp-qsn)/(1._r8 + cons28*qsn/(cpp*rv*ttmp**2))
                  dum=min(dum,0._r8)

                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
                  prds(k)=dum*(1._r8-dum1)/deltat/cldmax(i,k)
               end if
            end if
         end if

         ! bergeron process - evaporation of droplets and deposition onto snow

         if (qniic(i,k).ge.qsmall.and.qcic(i,k).ge.qsmall.and.t(i,k).lt.tmelt) then
            qvi=svp_to_qsat(esi(i,k), p(i,k))
            qvs=svp_to_qsat(esl(i,k), p(i,k))
            dqsidt =  xxls*qvi/(rv*t(i,k)**2)
            abi = 1._r8+dqsidt*xxls/cpp
            epss = 2._r8*pi*n0s(k)*rho(i,k)*Dv(i,k)* &
                 (f1s/(lams(k)*lams(k))+ &
                 f2s*(asn(i,k)*rho(i,k)/mu(i,k))**0.5_r8* &
                 sc(i,k)**(1._r8/3._r8)*cons14/ &
                 (lams(k)**(5._r8/2._r8+bs/2._r8)))
            bergs(k)=epss*(qvs-qvi)/abi
         else
            bergs(k)=0._r8
         end if

         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! conservation to ensure no negative values of cloud water/precipitation
         ! in case microphysical process rates are large

         ! make sure and use end-of-time step values for cloud water, ice, due
         ! condensation/deposition

         ! note: for check on conservation, processes are multiplied by omsm
         ! to prevent problems due to round off error

         ! include mixing timescale  (mtime)

         qce=(qc(i,k) - berg(i,k)*deltat)
         nce=(nc(i,k)+npccn(k)*deltat*mtime)
         qie=(qi(i,k)+(cmei(i,k)+berg(i,k))*deltat)
         nie=(ni(i,k)+nnuccd(k)*deltat*mtime)

         ! conservation of qc

         if (micro_mg1_0_colzero_use_native_impl) then
            dum = (prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+ &
                 psacws(k)+bergs(k))*lcldm(i,k)*deltat

            if (dum.gt.qce) then
               ratio = qce/deltat/lcldm(i,k)/(prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+psacws(k)+bergs(k))*omsm

               prc(k) = prc(k)*ratio
               pra(k) = pra(k)*ratio
               mnuccc(k) = mnuccc(k)*ratio
               mnucct(k) = mnucct(k)*ratio
               msacwi(k) = msacwi(k)*ratio
               psacws(k) = psacws(k)*ratio
               bergs(k) = bergs(k)*ratio
            end if

            ! conservation of nc

            dum = (nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+ &
                 npsacws(k)-nsubc(k))*lcldm(i,k)*deltat

            if (dum.gt.nce) then
               ratio = nce/deltat/((nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+&
                    npsacws(k)-nsubc(k))*lcldm(i,k))*omsm

               nprc1(k) = nprc1(k)*ratio
               npra(k) = npra(k)*ratio
               nnuccc(k) = nnuccc(k)*ratio
               nnucct(k) = nnucct(k)*ratio
               npsacws(k) = npsacws(k)*ratio
               nsubc(k)=nsubc(k)*ratio
            end if

            ! conservation of qi

            if (do_cldice) then

               frztmp = -mnuccc(k) - mnucct(k) - msacwi(k)
               if (use_hetfrz_classnuc) frztmp = -mnuccc(k)-mnucct(k)-mnudep(k)-msacwi(k)
               dum = ( frztmp*lcldm(i,k) + (prci(k)+prai(k))*icldm(i,k) )*deltat

               if (dum.gt.qie) then

                  frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
                  if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
                  ratio = (qie/deltat + frztmp*lcldm(i,k))/((prci(k)+prai(k))*icldm(i,k))*omsm
                  prci(k) = prci(k)*ratio
                  prai(k) = prai(k)*ratio
               end if

               ! conservation of ni
               frztmp = -nnucct(k) - nsacwi(k)
               if (use_hetfrz_classnuc) frztmp = -nnucct(k) - nnuccc(k) - nnudep(k) - nsacwi(k)
               dum = ( frztmp*lcldm(i,k) + (nprci(k)+nprai(k)-nsubi(k))*icldm(i,k) )*deltat

               if (dum.gt.nie) then

                  frztmp = nnucct(k) + nsacwi(k)
                  if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
                  ratio = (nie/deltat + frztmp*lcldm(i,k))/ &
                        ((nprci(k)+nprai(k)-nsubi(k))*icldm(i,k))*omsm
                  nprci(k) = nprci(k)*ratio
                  nprai(k) = nprai(k)*ratio
                  nsubi(k) = nsubi(k)*ratio
               end if
            end if

            ! for precipitation conservation, use logic that vertical integral
            ! of tendency from current level to top of model (i.e., qrtot) cannot be negative

            ! conservation of rain mixing rat

            if (((prc(k)+pra(k))*lcldm(i,k)+(-mnuccr(k)+pre(k)-pracs(k))*&
                 cldmax(i,k))*dz(i,k)*rho(i,k)+qrtot.lt.0._r8) then

               if (-pre(k)+pracs(k)+mnuccr(k).ge.qsmall) then

                  ratio = (qrtot/(dz(i,k)*rho(i,k))+(prc(k)+pra(k))*lcldm(i,k))/&
                       ((-pre(k)+pracs(k)+mnuccr(k))*cldmax(i,k))*omsm

                  pre(k) = pre(k)*ratio
                  pracs(k) = pracs(k)*ratio
                  mnuccr(k) = mnuccr(k)*ratio
               end if
            end if

            ! conservation of nr
            ! for now neglect evaporation of nr
            nsubr(k)=0._r8

            if ((nprc(k)*lcldm(i,k)+(-nnuccr(k)+nsubr(k)-npracs(k)&
                 +nragg(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+nrtot.lt.0._r8) then

               if (-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k).ge.qsmall) then

                  ratio = (nrtot/(dz(i,k)*rho(i,k))+nprc(k)*lcldm(i,k))/&
                       ((-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k))*cldmax(i,k))*omsm

                  nsubr(k) = nsubr(k)*ratio
                  npracs(k) = npracs(k)*ratio
                  nnuccr(k) = nnuccr(k)*ratio
                  nragg(k) = nragg(k)*ratio
               end if
            end if

            ! conservation of snow mix ratio

            if (((bergs(k)+psacws(k))*lcldm(i,k)+(prai(k)+prci(k))*icldm(i,k)+(pracs(k)+&
                 mnuccr(k)+prds(k))*cldmax(i,k))*dz(i,k)*rho(i,k)+qstot.lt.0._r8) then

               if (-prds(k).ge.qsmall) then

                  ratio = (qstot/(dz(i,k)*rho(i,k))+(bergs(k)+psacws(k))*lcldm(i,k)+(prai(k)+prci(k))*icldm(i,k)+&
                       (pracs(k)+mnuccr(k))*cldmax(i,k))/(-prds(k)*cldmax(i,k))*omsm

                  prds(k) = prds(k)*ratio
               end if
            end if

            ! conservation of ns

            ! calculate loss of number due to sublimation
            ! for now neglect sublimation of ns
            nsubs(k)=0._r8

            if ((nprci(k)*icldm(i,k)+(nnuccr(k)+nsubs(k)+nsagg(k))*cldmax(i,k))*&
                 dz(i,k)*rho(i,k)+nstot.lt.0._r8) then

               if (-nsubs(k)-nsagg(k).ge.qsmall) then

                  ratio = (nstot/(dz(i,k)*rho(i,k))+nprci(k)*icldm(i,k)+&
                       nnuccr(k)*cldmax(i,k))/((-nsubs(k)-nsagg(k))*cldmax(i,k))*omsm

                  nsubs(k) = nsubs(k)*ratio
                  nsagg(k) = nsagg(k)*ratio
               end if
            end if
         else
            call micro_mg1_0_conservation_limiter_codon_wrap(i, k, pcols, pver, deltat, omsm, qsmall, &
                 qce, nce, qie, nie, qrtot, nrtot, qstot, nstot, do_cldice, use_hetfrz_classnuc, &
                 lcldm, icldm, cldmax, dz, rho, prc, pra, mnuccc, mnucct, msacwi, psacws, bergs, &
                 nprc1, npra, nnuccc, nnucct, npsacws, nsubc, prci, prai, mnudep, nprci, nprai, &
                 nsubi, nnudep, nsacwi, mnuccr, pre, pracs, nsubr, npracs, nnuccr, nragg, prds, &
                 nsubs, nsagg, nprc)
         end if

         ! get tendencies due to microphysical conversion processes
         ! note: tendencies are multiplied by appropaiate cloud/precip 
         ! fraction to get grid-scale values
         ! note: cmei is already grid-average values

         qvlat(i,k) = qvlat(i,k)-(pre(k)+prds(k))*cldmax(i,k)-cmei(i,k) 

         tlat(i,k) = tlat(i,k)+((pre(k)*cldmax(i,k)) &
              *xxlv+(prds(k)*cldmax(i,k)+cmei(i,k))*xxls+ &
              ((bergs(k)+psacws(k)+mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(i,k)+(mnuccr(k)+ &
              pracs(k))*cldmax(i,k)+berg(i,k))*xlf)

         qctend(i,k) = qctend(i,k)+ &
              (-pra(k)-prc(k)-mnuccc(k)-mnucct(k)-msacwi(k)- & 
              psacws(k)-bergs(k))*lcldm(i,k)-berg(i,k)

         if (do_cldice) then

            frztmp = mnuccc(k) + mnucct(k) + msacwi(k)
            if (use_hetfrz_classnuc) frztmp = mnuccc(k) + mnucct(k) + mnudep(k) + msacwi(k)
            qitend(i,k) = qitend(i,k) + frztmp*lcldm(i,k) + &
               (-prci(k)-prai(k))*icldm(i,k) + cmei(i,k) + berg(i,k)

         end if

         qrtend(i,k) = qrtend(i,k)+ &
              (pra(k)+prc(k))*lcldm(i,k)+(pre(k)-pracs(k)- &
              mnuccr(k))*cldmax(i,k)

         qnitend(i,k) = qnitend(i,k)+ &
              (prai(k)+prci(k))*icldm(i,k)+(psacws(k)+bergs(k))*lcldm(i,k)+(prds(k)+ &
              pracs(k)+mnuccr(k))*cldmax(i,k)

         if (micro_mg1_0_colzero_use_native_impl) then
            !water tracers/isotopes:
            qrtend_copy(i,k) = qrtend_copy(i,k) + qrtend(i,k)
            qnitend_copy(i,k) = qnitend_copy(i,k) + qnitend(i,k)

            ! add output for cmei (accumulate)
            cmeiout(i,k) = cmeiout(i,k) + cmei(i,k)

            ! assign variables for trop_mozart, these are grid-average
            ! evaporation/sublimation is stored here as positive term

            evapsnow(i,k) = evapsnow(i,k)-prds(k)*cldmax(i,k)
            nevapr(i,k) = nevapr(i,k)-pre(k)*cldmax(i,k)
            nevapr2(i,k) = nevapr2(i,k)-pre(k)*cldmax(i,k)

            ! change to make sure prain is positive: do not remove snow from
            ! prain used for wet deposition
            prain(i,k) = prain(i,k)+(pra(k)+prc(k))*lcldm(i,k)+(-pracs(k)- &
                 mnuccr(k))*cldmax(i,k)
            prodsnow(i,k) = prodsnow(i,k)+(prai(k)+prci(k))*icldm(i,k)+(psacws(k)+bergs(k))*lcldm(i,k)+(&
                 pracs(k)+mnuccr(k))*cldmax(i,k)

            ! following are used to calculate 1st order conversion rate of cloud water
            !    to rain and snow (1/s), for later use in aerosol wet removal routine
            ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
            !    used to calculate pra, prc, ... in this routine
            ! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
            !                      (no cloud ice or bergeron terms)
            ! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }

            qcsinksum_rate1ord(k) = qcsinksum_rate1ord(k) + (pra(k)+prc(k)+psacws(k))*lcldm(i,k)
            qcsum_rate1ord(k) = qcsum_rate1ord(k) + qc(i,k)

            ! microphysics output, note this is grid-averaged
            prao(i,k)=prao(i,k)+pra(k)*lcldm(i,k)
            prco(i,k)=prco(i,k)+prc(k)*lcldm(i,k)
            mnuccco(i,k)=mnuccco(i,k)+mnuccc(k)*lcldm(i,k)
            mnuccto(i,k)=mnuccto(i,k)+mnucct(k)*lcldm(i,k)
            mnuccdo(i,k)=mnuccdo(i,k)+mnuccd(k)*lcldm(i,k)
            msacwio(i,k)=msacwio(i,k)+msacwi(k)*lcldm(i,k)
            psacwso(i,k)=psacwso(i,k)+psacws(k)*lcldm(i,k)
            bergso(i,k)=bergso(i,k)+bergs(k)*lcldm(i,k)
            bergo(i,k)=bergo(i,k)+berg(i,k)
            prcio(i,k)=prcio(i,k)+prci(k)*icldm(i,k)
            praio(i,k)=praio(i,k)+prai(k)*icldm(i,k)
            mnuccro(i,k)=mnuccro(i,k)+mnuccr(k)*cldmax(i,k)
            pracso (i,k)=pracso (i,k)+pracs (k)*cldmax(i,k)
            !water tracers:
            preo(i,k)=preo(i,k)+pre(k)*cldmax(i,k)
            prdso(i,k)=prdso(i,k)+prds(k)*cldmax(i,k)
         else
            call micro_mg1_0_process_output_accum_codon_wrap(i, k, pcols, pver, qrtend, qnitend, &
                 qrtend_copy, qnitend_copy, cmei, cmeiout, prds, pre, cldmax, evapsnow, nevapr, &
                 nevapr2, pra, prc, lcldm, pracs, mnuccr, prain, prai, prci, icldm, psacws, bergs, &
                 prodsnow, qcsinksum_rate1ord, qcsum_rate1ord, qc, prao, prco, mnuccc, mnucct, &
                 mnuccd, msacwi, mnuccco, mnuccto, mnuccdo, msacwio, psacwso, bergso, berg, &
                 bergo, prcio, praio, mnuccro, pracso, preo, prdso)
         end if

         ! multiply activation/nucleation by mtime to account for fast timescale

         nctend(i,k) = nctend(i,k)+ npccn(k)*mtime+&
              (-nnuccc(k)-nnucct(k)-npsacws(k)+nsubc(k) & 
              -npra(k)-nprc1(k))*lcldm(i,k)      

         if (do_cldice) then

            frztmp = nnucct(k) + nsacwi(k)
            if (use_hetfrz_classnuc) frztmp = nnucct(k) + nnuccc(k) + nnudep(k) + nsacwi(k)
            nitend(i,k) = nitend(i,k) + nnuccd(k)*mtime + & 
                  frztmp*lcldm(i,k) + (nsubi(k)-nprci(k)-nprai(k))*icldm(i,k)

         end if

         nstend(i,k) = nstend(i,k)+(nsubs(k)+ &
              nsagg(k)+nnuccr(k))*cldmax(i,k)+nprci(k)*icldm(i,k)

         nrtend(i,k) = nrtend(i,k)+ &
              nprc(k)*lcldm(i,k)+(nsubr(k)-npracs(k)-nnuccr(k) &
              +nragg(k))*cldmax(i,k)

         ! make sure that nc and ni at advanced time step do not exceed
         ! maximum (existing N + source terms*dt), which is possible due to
         ! fast nucleation timescale

         if (nctend(i,k).gt.0._r8.and.nc(i,k)+nctend(i,k)*deltat.gt.ncmax) then
            nctend(i,k)=max(0._r8,(ncmax-nc(i,k))/deltat)
         end if

         if (do_cldice .and. nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax) then
            nitend(i,k)=max(0._r8,(nimax-ni(i,k))/deltat)
         end if

         ! get final values for precipitation q and N, based on
         ! flux of precip from above, source/sink term, and terminal fallspeed
         ! see eq. 15-16 in MG2008

         ! rain

         if (qric(i,k).ge.qsmall) then
            if (k.eq.top_lev) then
               qric(i,k)=qrtend(i,k)*dz(i,k)/cldmax(i,k)/umr(k)
               nric(i,k)=nrtend(i,k)*dz(i,k)/cldmax(i,k)/unr(k)
            else
               qric(i,k) = (rho(i,k-1)*umr(k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
                    (rho(i,k)*dz(i,k)*qrtend(i,k)))/(umr(k)*rho(i,k)*cldmax(i,k))
               nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
                    (rho(i,k)*dz(i,k)*nrtend(i,k)))/(unr(k)*rho(i,k)*cldmax(i,k))

            end if
         else
            qric(i,k)=0._r8
            nric(i,k)=0._r8
         end if

         ! snow

         if (qniic(i,k).ge.qsmall) then
            if (k.eq.top_lev) then
               qniic(i,k)=qnitend(i,k)*dz(i,k)/cldmax(i,k)/ums(k)
               nsic(i,k)=nstend(i,k)*dz(i,k)/cldmax(i,k)/uns(k)
            else
               qniic(i,k) = (rho(i,k-1)*ums(k-1)*qniic(i,k-1)*cldmax(i,k-1)+ &
                    (rho(i,k)*dz(i,k)*qnitend(i,k)))/(ums(k)*rho(i,k)*cldmax(i,k))
               nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
                    (rho(i,k)*dz(i,k)*nstend(i,k)))/(uns(k)*rho(i,k)*cldmax(i,k))
            end if
         else
            qniic(i,k)=0._r8
            nsic(i,k)=0._r8
         end if

         ! calculate precipitation flux at surface
         ! divide by density of water to get units of m/s

         prect(i) = prect(i)+(qrtend(i,k)*dz(i,k)*rho(i,k)+&
              qnitend(i,k)*dz(i,k)*rho(i,k))/rhow
         preci(i) = preci(i)+qnitend(i,k)*dz(i,k)*rho(i,k)/rhow

         ! convert rain rate from m/s to mm/hr

         rainrt(i,k)=qric(i,k)*rho(i,k)*umr(k)/rhow*3600._r8*1000._r8

         ! vertically-integrated precip source/sink terms (note: grid-averaged)

         qrtot = max(qrtot+qrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
         qstot = max(qstot+qnitend(i,k)*dz(i,k)*rho(i,k),0._r8)
         nrtot = max(nrtot+nrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
         nstot = max(nstot+nstend(i,k)*dz(i,k)*rho(i,k),0._r8)

         ! calculate melting and freezing of precip

         ! melt snow at +2 C

         if (t(i,k)+tlat(i,k)/cpp*deltat > 275.15_r8) then
            if (qstot > 0._r8) then

               ! make sure melting snow doesn't reduce temperature below threshold
               dum = -xlf/cpp*qstot/(dz(i,k)*rho(i,k))
               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.275.15_r8) then
                  dum = (t(i,k)+tlat(i,k)/cpp*deltat-275.15_r8)*cpp/xlf
                  dum = dum/(xlf/cpp*qstot/(dz(i,k)*rho(i,k)))
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qric(i,k)=qric(i,k)+dum*qniic(i,k)
               nric(i,k)=nric(i,k)+dum*nsic(i,k)
               qniic(i,k)=(1._r8-dum)*qniic(i,k)
               nsic(i,k)=(1._r8-dum)*nsic(i,k)

               !water tracers/isotopes:
               !original (Chuck Bardeen) version (ensures column integral matches final distribution) - JN
            !   meltso(i,:)  = meltso(i,:) + dum*qnitend_copy(i,:)  
            !   qnitend_copy(i,:) = qnitend_copy(i,:) - dum*qnitend_copy(i,:)
            !   qrtend_copy(i,:)  = qrtend_copy(i,:)  + dum*qnitend_copy(i,:) 

               !save total melted amount at this vertical level:
               meltso(i,k) = meltso(i,k) + dum*qstot*g  
 
               ! heating tendency 
               tmp=-xlf*dum*qstot/(dz(i,k)*rho(i,k))
               meltsdt(i,k)=meltsdt(i,k) + tmp

               tlat(i,k)=tlat(i,k)+tmp
               qrtot=qrtot+dum*qstot
               nrtot=nrtot+dum*nstot
               qstot=(1._r8-dum)*qstot
               nstot=(1._r8-dum)*nstot
               preci(i)=(1._r8-dum)*preci(i)
            end if
         end if

         ! freeze all rain at -5C for Arctic

         if (t(i,k)+tlat(i,k)/cpp*deltat < (tmelt - 5._r8)) then

            if (qrtot > 0._r8) then

               ! make sure freezing rain doesn't increase temperature above threshold
               dum = xlf/cpp*qrtot/(dz(i,k)*rho(i,k))
               if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.(tmelt - 5._r8)) then
                  dum = -(t(i,k)+tlat(i,k)/cpp*deltat-(tmelt-5._r8))*cpp/xlf
                  dum = dum/(xlf/cpp*qrtot/(dz(i,k)*rho(i,k)))
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qniic(i,k)=qniic(i,k)+dum*qric(i,k)
               nsic(i,k)=nsic(i,k)+dum*nric(i,k)
               qric(i,k)=(1._r8-dum)*qric(i,k)
               nric(i,k)=(1._r8-dum)*nric(i,k)

               !water tracers/isotopes:
               !original (Chuck Bardeen) version (ensures column integral matches final distribution) - JN
             !  frzro(i,:)   = frzro(i,:) + dum*qrtend_copy(i,:)
             !  qnitend_copy(i,:) = qnitend_copy(i,:) + dum*qrtend_copy(i,:)
             !  qrtend_copy(i,:)  = qrtend_copy(i,:) - dum*qrtend_copy(i,:)

               !save total frozen amount at this vertical level:
               frzro(i,k) = frzro(i,k) + dum*qrtot*g 

               ! heating tendency 
               tmp = xlf*dum*qrtot/(dz(i,k)*rho(i,k))
               frzrdt(i,k)=frzrdt(i,k) + tmp

               tlat(i,k)=tlat(i,k)+tmp
               qstot=qstot+dum*qrtot
               qrtot=(1._r8-dum)*qrtot
               nstot=nstot+dum*nrtot
               nrtot=(1._r8-dum)*nrtot
               preci(i)=preci(i)+dum*(prect(i)-preci(i))
            end if
         end if

         ! if rain/snow mix ratio is zero so should number concentration

         if (qniic(i,k).lt.qsmall) then
            qniic(i,k)=0._r8
            nsic(i,k)=0._r8
         end if

         if (qric(i,k).lt.qsmall) then
            qric(i,k)=0._r8
            nric(i,k)=0._r8
         end if

         ! make sure number concentration is a positive number to avoid 
         ! taking root of negative

         nric(i,k)=max(nric(i,k),0._r8)
         nsic(i,k)=max(nsic(i,k),0._r8)

         !.......................................................................
         ! get size distribution parameters for fallspeed calculations
         !......................................................................
         ! rain

         if (qric(i,k).ge.qsmall) then
            lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
            n0r(k) = nric(i,k)*lamr(k)

            ! check for slope
            ! change lammax and lammin for rain and snow
            ! adjust vars

            if (lamr(k).lt.lamminr) then

               lamr(k) = lamminr

               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
            else if (lamr(k).gt.lammaxr) then
               lamr(k) = lammaxr
               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
            end if


            ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

            unr(k) = min(arn(i,k)*cons4/lamr(k)**br,9.1_r8*rhof(i,k))
            umr(k) = min(arn(i,k)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(i,k))

         else
            lamr(k) = 0._r8
            n0r(k) = 0._r8
            umr(k)=0._r8
            unr(k)=0._r8
         end if

         !calculate mean size of combined rain and snow

         if (lamr(k).gt.0._r8) then
            Artmp = n0r(k) * pi / (2._r8 * lamr(k)**3._r8)
         else 
            Artmp = 0._r8
         endif

         if (lamc(k).gt.0._r8) then
            Actmp = cdist1(k) * pi * gamma(pgam(k)+3._r8)/(4._r8 * lamc(k)**2._r8)
         else 
            Actmp = 0._r8
         endif

         if (Actmp.gt.0_r8.or.Artmp.gt.0) then
            rercld(i,k)=rercld(i,k) + 3._r8 *(qric(i,k) + qcic(i,k)) / (4._r8 * rhow * (Actmp + Artmp))
            arcld(i,k)=arcld(i,k)+1._r8
         endif

         !......................................................................
         ! snow

         if (qniic(i,k).ge.qsmall) then
            lams(k) = (cons6*cs*nsic(i,k)/ &
                 qniic(i,k))**(1._r8/ds)
            n0s(k) = nsic(i,k)*lams(k)

            ! check for slope
            ! adjust vars

            if (lams(k).lt.lammins) then
               lams(k) = lammins
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
               nsic(i,k) = n0s(k)/lams(k)

            else if (lams(k).gt.lammaxs) then
               lams(k) = lammaxs
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*cons6)
               nsic(i,k) = n0s(k)/lams(k)
            end if

            ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)

            ums(k) = min(asn(i,k)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof(i,k))
            uns(k) = min(asn(i,k)*cons7/lams(k)**bs,1.2_r8*rhof(i,k))

         else
            lams(k) = 0._r8
            n0s(k) = 0._r8
            ums(k) = 0._r8
            uns(k) = 0._r8
         end if

         if (micro_mg1_0_colzero_use_native_impl) then
            !c........................................................................
            ! sum over sub-step for average process rates

            ! convert rain/snow q and N for output to history, note,
            ! output is for gridbox average

            qrout(i,k)=qrout(i,k)+qric(i,k)*cldmax(i,k)
            qsout(i,k)=qsout(i,k)+qniic(i,k)*cldmax(i,k)
            nrout(i,k)=nrout(i,k)+nric(i,k)*rho(i,k)*cldmax(i,k)
            nsout(i,k)=nsout(i,k)+nsic(i,k)*rho(i,k)*cldmax(i,k)

            tlat1(i,k)=tlat1(i,k)+tlat(i,k)
            qvlat1(i,k)=qvlat1(i,k)+qvlat(i,k)
            qctend1(i,k)=qctend1(i,k)+qctend(i,k)
            qitend1(i,k)=qitend1(i,k)+qitend(i,k)
            nctend1(i,k)=nctend1(i,k)+nctend(i,k)
            nitend1(i,k)=nitend1(i,k)+nitend(i,k)

            t(i,k)=t(i,k)+tlat(i,k)*deltat/cpp
            q(i,k)=q(i,k)+qvlat(i,k)*deltat
            qc(i,k)=qc(i,k)+qctend(i,k)*deltat
            qi(i,k)=qi(i,k)+qitend(i,k)*deltat
            nc(i,k)=nc(i,k)+nctend(i,k)*deltat
            ni(i,k)=ni(i,k)+nitend(i,k)*deltat

            rainrt1(i,k)=rainrt1(i,k)+rainrt(i,k)

            !divide rain radius over substeps for average
            if (arcld(i,k) .gt. 0._r8) then
               rercld(i,k)=rercld(i,k)/arcld(i,k)
            end if

            !calculate precip fluxes and adding them to summing sub-stepping variables
            !! flux is zero at top interface
            rflx(i,1)=0.0_r8
            sflx(i,1)=0.0_r8

            !! calculating the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
            rflx(i,k+1)=qrout(i,k)*rho(i,k)*umr(k)
            sflx(i,k+1)=qsout(i,k)*rho(i,k)*ums(k)

            !! add to summing sub-stepping variable
            rflx1(i,k+1)=rflx1(i,k+1)+rflx(i,k+1)
            sflx1(i,k+1)=sflx1(i,k+1)+sflx(i,k+1)
         end if

         !c........................................................................

      end do ! k loop

      if (.not. micro_mg1_0_colzero_use_native_impl) then
         call micro_mg1_0_substep_accum_column_codon_wrap(i, pcols, pver, top_lev, deltat, cpp, &
              qric, qniic, nric, nsic, rho, cldmax, qrout, qsout, nrout, nsout, tlat, qvlat, qctend, &
              qitend, nctend, nitend, tlat1, qvlat1, qctend1, qitend1, nctend1, nitend1, t, q, qc, qi, &
              nc, ni, rainrt, rainrt1, arcld, rercld, rflx, sflx, rflx1, sflx1, umr, ums)
      end if

      prect1(i)=prect1(i)+prect(i)
      preci1(i)=preci1(i)+preci(i)

   end do ! it loop, sub-step

   if (micro_mg1_0_colzero_use_native_impl) then
      do k = top_lev, pver
         rate1ord_cw2pr_st(i,k) = qcsinksum_rate1ord(k)/max(qcsum_rate1ord(k),1.0e-30_r8)
      end do
   else
      call micro_mg1_0_rate1ord_column_codon_wrap(i, pcols, pver, top_lev, &
           qcsinksum_rate1ord, qcsum_rate1ord, rate1ord_cw2pr_st)
   end if

300 continue  ! continue if no cloud water
end do ! i loop

! convert dt from sub-step back to full time step
deltat=deltat*real(iter)

!c.............................................................................

do i=1,ncol

   ! skip all calculations if no cloud water
   if (ltrue(i).eq.0) then

      do k=1,top_lev-1
         ! assign zero values for effective radius above 1 mbar
         effc(i,k)=0._r8
         effi(i,k)=0._r8
         effc_fn(i,k)=0._r8
         lamcrad(i,k)=0._r8
         pgamrad(i,k)=0._r8
         deffi(i,k)=0._r8
      end do

      do k=top_lev,pver
         ! assign default values for effective radius
         effc(i,k)=10._r8
         effi(i,k)=25._r8
         effc_fn(i,k)=10._r8
         lamcrad(i,k)=0._r8
         pgamrad(i,k)=0._r8
         deffi(i,k)=0._r8
      end do
      goto 500
   end if

   ! initialize nstep for sedimentation sub-steps
   nstep = 1

   call micro_mg1_0_post_iter_avg_codon_wrap(i, pcols, pver, top_lev, iter, prect1, preci1, &
        prect, preci, t1, q1, qc1, qi1, nc1, ni1, t, q, qc, qi, nc, ni, &
        tlat1, qvlat1, qctend1, qitend1, nctend1, nitend1, tlat, qvlat, qctend, qitend, &
        nctend, nitend, rainrt1, rainrt, rflx1, sflx1, rflx, sflx, qrout, qsout, nrout, nsout, &
        nevapr, nevapr2, evapsnow, prain, prodsnow, cmeout, cmeiout, meltsdt, frzrdt, &
        prao, prco, mnuccco, mnuccto, msacwio, psacwso, bergso, bergo, prcio, praio, &
        mnuccro, pracso, mnuccdo, preo, prdso, frzro, meltso, wtprelat, prer_evap)

   do k=top_lev,pver
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! calculate sedimentation for cloud water and ice
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! update in-cloud cloud mixing ratio and number concentration 
      ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
      ! note: these are in-cloud values***, hence we divide by cloud fraction

      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)
      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)/icldm(i,k)
      dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k),0._r8)
      dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)/icldm(i,k),0._r8)

      ! obtain new slope parameter to avoid possible singularity

      if (dumi(i,k).ge.qsmall) then
         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)

         lami(k) = (cons1*ci* &
              dumni(i,k)/dumi(i,k))**(1._r8/di)
         lami(k)=max(lami(k),lammini)
         lami(k)=min(lami(k),lammaxi)
      else
         lami(k)=0._r8
      end if

      if (dumc(i,k).ge.qsmall) then
         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)
         ! add lower limit to in-cloud number concentration
         dumnc(i,k)=max(dumnc(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3 
         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
         pgam(k)=1._r8/(pgam(k)**2)-1._r8
         pgam(k)=max(pgam(k),2._r8)
         pgam(k)=min(pgam(k),15._r8)

         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
         lammin = (pgam(k)+1._r8)/50.e-6_r8
         lammax = (pgam(k)+1._r8)/2.e-6_r8
         lamc(k)=max(lamc(k),lammin)
         lamc(k)=min(lamc(k),lammax)
      else
         lamc(k)=0._r8
      end if

      ! calculate number and mass weighted fall velocity for droplets
      ! include effects of sub-grid distribution of cloud water


      if (dumc(i,k).ge.qsmall) then
         unc = acn(i,k)*gamma(1._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+1._r8))
         umc = acn(i,k)*gamma(4._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+4._r8))
         ! fallspeed for output
         vtrmc(i,k)=umc
      else
         umc = 0._r8
         unc = 0._r8
      end if

      ! calculate number and mass weighted fall velocity for cloud ice

      if (dumi(i,k).ge.qsmall) then
         uni =  ain(i,k)*cons16/lami(k)**bi
         umi = ain(i,k)*cons17/(6._r8*lami(k)**bi)
         uni=min(uni,1.2_r8*rhof(i,k))
         umi=min(umi,1.2_r8*rhof(i,k))

         ! fallspeed
         vtrmi(i,k)=umi
      else
         umi = 0._r8
         uni = 0._r8
      end if

      fi(k) = g*rho(i,k)*umi
      fni(k) = g*rho(i,k)*uni
      fc(k) = g*rho(i,k)*umc
      fnc(k) = g*rho(i,k)*unc

      !water tracers:
      wtfc(i,k) = fc(k)
      wtfi(i,k) = fi(k)

      ! calculate number of split time steps to ensure courant stability criteria
      ! for sedimentation calculations

      rgvm = max(fi(k),fc(k),fni(k),fnc(k))
      nstep = max(int(rgvm*deltat/pdel(i,k)+1._r8),nstep)

      ! redefine dummy variables - sedimentation is calculated over grid-scale
      ! quantities to ensure conservation

      dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
      dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
      dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat),0._r8)
      dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat),0._r8)

      if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
      if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8

   end do       !!! vertical loop
   do n = 1,nstep  !! loop over sub-time step to ensure stability

      do k = top_lev,pver
         if (do_cldice) then
            falouti(k) = fi(k)*dumi(i,k)
            faloutni(k) = fni(k)*dumni(i,k)
         else
            falouti(k)  = 0._r8
            faloutni(k) = 0._r8
         end if

         faloutc(k) = fc(k)*dumc(i,k)
         faloutnc(k) = fnc(k)*dumnc(i,k)
      end do

      ! top of model

      k = top_lev
      faltndi = falouti(k)/pdel(i,k)
      faltndni = faloutni(k)/pdel(i,k)
      faltndc = faloutc(k)/pdel(i,k)
      faltndnc = faloutnc(k)/pdel(i,k)

      ! add fallout terms to microphysical tendencies

      qitend(i,k) = qitend(i,k)-faltndi/nstep
      nitend(i,k) = nitend(i,k)-faltndni/nstep
      qctend(i,k) = qctend(i,k)-faltndc/nstep
      nctend(i,k) = nctend(i,k)-faltndnc/nstep

      ! sedimentation tendencies for output
      qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
      qisedten(i,k)=qisedten(i,k)-faltndi/nstep

      dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
      dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
      dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
      dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

      do k = top_lev+1,pver

         ! for cloud liquid and ice, if cloud fraction increases with height
         ! then add flux from above to both vapor and cloud water of current level
         ! this means that flux entering clear portion of cell from above evaporates
         ! instantly

         dum=lcldm(i,k)/lcldm(i,k-1)
         dum=min(dum,1._r8)
         dum1=icldm(i,k)/icldm(i,k-1)
         dum1=min(dum1,1._r8)

         faltndqie=(falouti(k)-falouti(k-1))/pdel(i,k)
         faltndi=(falouti(k)-dum1*falouti(k-1))/pdel(i,k)
         faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel(i,k)
         faltndqce=(faloutc(k)-faloutc(k-1))/pdel(i,k)
         faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(i,k)
         faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(i,k)

         ! add fallout terms to eulerian tendencies

         qitend(i,k) = qitend(i,k)-faltndi/nstep
         nitend(i,k) = nitend(i,k)-faltndni/nstep
         qctend(i,k) = qctend(i,k)-faltndc/nstep
         nctend(i,k) = nctend(i,k)-faltndnc/nstep

         ! sedimentation tendencies for output
         qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
         qisedten(i,k)=qisedten(i,k)-faltndi/nstep

         ! add terms to to evap/sub of cloud water

         qvlat(i,k)=qvlat(i,k)-(faltndqie-faltndi)/nstep
         ! for output
         qisevap(i,k)=qisevap(i,k)-(faltndqie-faltndi)/nstep
         qvlat(i,k)=qvlat(i,k)-(faltndqce-faltndc)/nstep
         ! for output
         qcsevap(i,k)=qcsevap(i,k)-(faltndqce-faltndc)/nstep

         tlat(i,k)=tlat(i,k)+(faltndqie-faltndi)*xxls/nstep
         tlat(i,k)=tlat(i,k)+(faltndqce-faltndc)*xxlv/nstep

         dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
         dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
         dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
         dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

         Fni(K)=MAX(Fni(K)/pdel(i,K),Fni(K-1)/pdel(i,K-1))*pdel(i,K)
         FI(K)=MAX(FI(K)/pdel(i,K),FI(K-1)/pdel(i,K-1))*pdel(i,K)
         fnc(k)=max(fnc(k)/pdel(i,k),fnc(k-1)/pdel(i,k-1))*pdel(i,k)
         Fc(K)=MAX(Fc(K)/pdel(i,K),Fc(K-1)/pdel(i,K-1))*pdel(i,K)

      end do   !! k loop

      ! units below are m/s
      ! cloud water/ice sedimentation flux at surface 
      ! is added to precip flux at surface to get total precip (cloud + precip water)
      ! rate

      prect(i) = prect(i)+(faloutc(pver)+falouti(pver))/g/nstep/1000._r8  
      preci(i) = preci(i)+(falouti(pver))/g/nstep/1000._r8

   end do   !! nstep loop

   ! end sedimentation
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   ! get new update for variables that includes sedimentation tendency
   ! note : here dum variables are grid-average, NOT in-cloud

   do k=top_lev,pver

      call micro_mg1_0_phase_change_codon_wrap(i, k, pcols, pver, deltat, cpp, xlf, tmelt, qsmall, &
           pi, rhow, do_cldice, qc, qi, nc, ni, t, qctend, qitend, nctend, nitend, tlat, &
           dumc, dumi, dumnc, dumni, melto, homoo, wtpostlat)

      if (do_cldice) then

         ! remove any excess over-saturation, which is possible due to non-linearity when adding 
         ! together all microphysical processes
         ! follow code similar to old CAM scheme

         qtmp=q(i,k)+qvlat(i,k)*deltat
         ttmp=t(i,k)+tlat(i,k)/cpp*deltat

         esn = svp_water(ttmp)  ! use rhw to allow ice supersaturation
         qsn = svp_to_qsat(esn, p(i,k))

         if (qtmp > qsn .and. qsn > 0) then
            ! expression below is approximate since there may be ice deposition
            dum = (qtmp-qsn)/(1._r8+cons27*qsn/(cpp*rv*ttmp**2))/deltat
            ! add to output cme
            cmeout(i,k) = cmeout(i,k)+dum
            ! now add to tendencies, partition between liquid and ice based on temperature
            if (ttmp > 268.15_r8) then
               dum1=0.0_r8
               ! now add to tendencies, partition between liquid and ice based on te
            else if (ttmp < 238.15_r8) then
               dum1=1.0_r8
            else
               dum1=(268.15_r8-ttmp)/30._r8
            end if

            dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                 *qsn/(cpp*rv*ttmp**2))/deltat
            qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
            ! for output
            qcreso(i,k)=dum*(1._r8-dum1)
            qitend(i,k)=qitend(i,k)+dum*dum1
            qireso(i,k)=dum*dum1
            qvlat(i,k)=qvlat(i,k)-dum
            ! for output
            qvres(i,k)=-dum
            tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
            !water tracers:
            wtpostlat(i,k) = wtpostlat(i,k)+(dum*(1._r8-dum1)*xxlv+dum*dum1*xxls)
         end if
      end if

      !...............................................................................
      ! calculate effective radius for pass to radiation code
      ! if no cloud water, default value is 10 micron for droplets,
      ! 25 micron for cloud ice

      ! update cloud variables after instantaneous processes to get effective radius
      ! variables are in-cloud to calculate size dist parameters

      dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)/lcldm(i,k)
      dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)/icldm(i,k)
      dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)/lcldm(i,k)
      dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)/icldm(i,k)

      ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1

      dumc(i,k)=min(dumc(i,k),5.e-3_r8)
      dumi(i,k)=min(dumi(i,k),5.e-3_r8)

      !...................
      ! cloud ice effective radius

      if (dumi(i,k).ge.qsmall) then
         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumni(i,k)=min(dumni(i,k),dumi(i,k)*1.e20_r8)
         lami(k) = (cons1*ci*dumni(i,k)/dumi(i,k))**(1._r8/di)

         if (lami(k).lt.lammini) then
            lami(k) = lammini
            n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*cons1)
            niic(i,k) = n0i(k)/lami(k)
            ! adjust number conc if needed to keep mean size in reasonable range
            if (do_cldice) nitend(i,k)=(niic(i,k)*icldm(i,k)-ni(i,k))/deltat

         else if (lami(k).gt.lammaxi) then
            lami(k) = lammaxi
            n0i(k) = lami(k)**(di+1._r8)*dumi(i,k)/(ci*cons1)
            niic(i,k) = n0i(k)/lami(k)
            ! adjust number conc if needed to keep mean size in reasonable range
            if (do_cldice) nitend(i,k)=(niic(i,k)*icldm(i,k)-ni(i,k))/deltat
         end if
         effi(i,k) = 1.5_r8/lami(k)*1.e6_r8

      else
         effi(i,k) = 25._r8
      end if

      ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
      ! radius has already been determined from the size distribution.
      if (.not. do_cldice) then
         effi(i,k) = re_ice(i,k) * 1e6_r8      ! m -> um
      end if

      !...................
      ! cloud droplet effective radius

      if (dumc(i,k).ge.qsmall) then

         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumnc(i,k)=min(dumnc(i,k),dumc(i,k)*1.e20_r8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! set tendency to ensure minimum droplet concentration
         ! after update by microphysics, except when lambda exceeds bounds on mean drop
         ! size or if there is no cloud water
         if (dumnc(i,k).lt.cdnl/rho(i,k)) then   
            nctend(i,k)=(cdnl/rho(i,k)*lcldm(i,k)-nc(i,k))/deltat   
         end if
         dumnc(i,k)=max(dumnc(i,k),cdnl/rho(i,k)) ! sghan minimum in #/cm3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
         pgam(k)=1._r8/(pgam(k)**2)-1._r8
         pgam(k)=max(pgam(k),2._r8)
         pgam(k)=min(pgam(k),15._r8)

         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
         lammin = (pgam(k)+1._r8)/50.e-6_r8
         ! Multiply by omsm to fit within RRTMG's table.
         lammax = (pgam(k)+1._r8)*omsm/2.e-6_r8
         if (lamc(k).lt.lammin) then
            lamc(k) = lammin
            ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)* &
                 gamma(pgam(k)+1._r8)/ &
                 (pi*rhow*gamma(pgam(k)+4._r8))
            ! adjust number conc if needed to keep mean size in reasonable range
            nctend(i,k)=(ncic(i,k)*lcldm(i,k)-nc(i,k))/deltat

         else if (lamc(k).gt.lammax) then
            lamc(k) = lammax
            ncic(i,k) = 6._r8*lamc(k)**3*dumc(i,k)* &
                 gamma(pgam(k)+1._r8)/ &
                 (pi*rhow*gamma(pgam(k)+4._r8))
            ! adjust number conc if needed to keep mean size in reasonable range
            nctend(i,k)=(ncic(i,k)*lcldm(i,k)-nc(i,k))/deltat
         end if

         effc(i,k) = &
              gamma(pgam(k)+4._r8)/ &
              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
         !assign output fields for shape here
         lamcrad(i,k)=lamc(k)
         pgamrad(i,k)=pgam(k)

      else
         effc(i,k) = 10._r8
         lamcrad(i,k)=0._r8
         pgamrad(i,k)=0._r8
      end if

      ! ice effective diameter for david mitchell's optics
      if (do_cldice) then
         deffi(i,k)=effi(i,k)*rhoi/917._r8*2._r8
      else
         deffi(i,k)=effi(i,k) * 2._r8
      end if


!!! recalculate effective radius for constant number, in order to separate
      ! first and second indirect effects
      ! assume constant number of 10^8 kg-1

      dumnc(i,k)=1.e8_r8

      if (dumc(i,k).ge.qsmall) then
         pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
         pgam(k)=1._r8/(pgam(k)**2)-1._r8
         pgam(k)=max(pgam(k),2._r8)
         pgam(k)=min(pgam(k),15._r8)

         lamc(k) = (pi/6._r8*rhow*dumnc(i,k)*gamma(pgam(k)+4._r8)/ &
              (dumc(i,k)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
         lammin = (pgam(k)+1._r8)/50.e-6_r8
         lammax = (pgam(k)+1._r8)/2.e-6_r8
         if (lamc(k).lt.lammin) then
            lamc(k) = lammin
         else if (lamc(k).gt.lammax) then
            lamc(k) = lammax
         end if
         effc_fn(i,k) = &
              gamma(pgam(k)+4._r8)/ &
              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8

      else
         effc_fn(i,k) = 10._r8
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

   end do ! vertical k loop

500 continue

   if (micro_mg1_0_colzero_use_native_impl) then
      do k=top_lev,pver
         ! if updated q (after microphysics) is zero, then ensure updated n is also zero

         if (qc(i,k)+qctend(i,k)*deltat.lt.qsmall) nctend(i,k)=-nc(i,k)/deltat
         if (do_cldice .and. qi(i,k)+qitend(i,k)*deltat.lt.qsmall) nitend(i,k)=-ni(i,k)/deltat
      end do
   end if

end do ! i loop

if (.not. micro_mg1_0_colzero_use_native_impl) then
   call micro_mg1_0_number_cleanup_codon_wrap(ncol, pcols, pver, top_lev, deltat, qsmall, &
        do_cldice, qc, qi, nc, ni, qctend, qitend, nctend, nitend)
end if

! add snow ouptut
do i = 1,ncol
   do k=top_lev,pver
      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
         dsout(i,k)=3._r8*rhosn/917._r8*(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
      endif
   end do
end do

!calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual
do i = 1,ncol
   do k=top_lev,pver
      !! RAIN
      if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
         reff_rain(i,k)=1.5_r8*(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)*1.e6_r8
      endif
      !! SNOW
      if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
         reff_snow(i,k)=1.5_r8*(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)*1.e6_r8
      end if
   end do
end do

! analytic radar reflectivity
! formulas from Matthew Shupe, NOAA/CERES
! *****note: radar reflectivity is local (in-precip average)
! units of mm^6/m^3

do i = 1,ncol
   do k=top_lev,pver
      if (qc(i,k)+qctend(i,k)*deltat.ge.qsmall) then
         dum=((qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)*1000._r8)**2 &
              /(0.109_r8*(nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)/1.e6_r8)*lcldm(i,k)/cldmax(i,k)
      else
         dum=0._r8
      end if
      if (qi(i,k)+qitend(i,k)*deltat.ge.qsmall) then
         dum1=((qi(i,k)+qitend(i,k)*deltat)*rho(i,k)/icldm(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm(i,k)/cldmax(i,k)
      else 
         dum1=0._r8
      end if

      if (qsout(i,k).ge.qsmall) then
         dum1=dum1+(qsout(i,k)*rho(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
      end if

      refl(i,k)=dum+dum1

      ! add rain rate, but for 37 GHz formulation instead of 94 GHz
      ! formula approximated from data of Matrasov (2007)
      ! rainrt is the rain rate in mm/hr
      ! reflectivity (dum) is in DBz

      if (rainrt(i,k).ge.0.001_r8) then
         dum=log10(rainrt(i,k)**6._r8)+16._r8

         ! convert from DBz to mm^6/m^3

         dum = 10._r8**(dum/10._r8)
      else
         ! don't include rain rate in R calculation for values less than 0.001 mm/hr
         dum=0._r8
      end if

      ! add to refl

      refl(i,k)=refl(i,k)+dum

      !output reflectivity in Z.
      areflz(i,k)=refl(i,k)

      ! convert back to DBz 

      if (refl(i,k).gt.minrefl) then 
         refl(i,k)=10._r8*log10(refl(i,k))
      else
         refl(i,k)=-9999._r8
      end if

      if (micro_mg1_0_colzero_use_native_impl) then
         !set averaging flag
         if (refl(i,k).gt.mindbz) then
            arefl(i,k)=refl(i,k)
            frefl(i,k)=1.0_r8
         else
            arefl(i,k)=0._r8
            areflz(i,k)=0._r8
            frefl(i,k)=0._r8
         end if

         ! bound cloudsat reflectivity

         csrfl(i,k)=min(csmax,refl(i,k))

         !set averaging flag
         if (csrfl(i,k).gt.csmin) then
            acsrfl(i,k)=refl(i,k)
            fcsrfl(i,k)=1.0_r8
         else
            acsrfl(i,k)=0._r8
            fcsrfl(i,k)=0._r8
         end if
      end if

   end do
end do

if (.not. micro_mg1_0_colzero_use_native_impl) then
   call micro_mg1_0_reflectivity_flags_codon_wrap(ncol, pcols, pver, top_lev, mindbz, csmin, csmax, &
        refl, arefl, areflz, frefl, csrfl, acsrfl, fcsrfl)
end if


! averaging for snow and rain number and diameter

call micro_mg1_0_select_tail_diag_impl()
if (micro_mg1_0_tail_diag_use_native_impl) then
   qrout2(:,:)=0._r8
   qsout2(:,:)=0._r8
   nrout2(:,:)=0._r8
   nsout2(:,:)=0._r8
   drout2(:,:)=0._r8
   dsout2(:,:)=0._r8
   freqs(:,:)=0._r8
   freqr(:,:)=0._r8
   do i = 1,ncol
      do k=top_lev,pver
         if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
            qrout2(i,k)=qrout(i,k)
            nrout2(i,k)=nrout(i,k)
            drout2(i,k)=(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)
            freqr(i,k)=1._r8
         endif
         if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
            qsout2(i,k)=qsout(i,k)
            nsout2(i,k)=nsout(i,k)
            dsout2(i,k)=(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
            freqs(i,k)=1._r8
         endif
      end do
   end do
else
   call micro_mg1_0_tail_avg_codon_wrap(ncol, pcols, pver, top_lev, qrout, qsout, nrout, nsout, &
        qrout2, qsout2, nrout2, nsout2, drout2, dsout2, freqs, freqr)
   do i = 1,ncol
      do k=top_lev,pver
         if (qrout(i,k).gt.1.e-7_r8.and.nrout(i,k).gt.0._r8) then
            drout2(i,k)=(pi * rhow * nrout(i,k)/qrout(i,k))**(-1._r8/3._r8)
         endif
         if (qsout(i,k).gt.1.e-7_r8.and.nsout(i,k).gt.0._r8) then
            dsout2(i,k)=(pi * rhosn * nsout(i,k)/qsout(i,k))**(-1._r8/3._r8)
         endif
      end do
   end do
end if

! output activated liquid and ice (convert from #/kg -> #/m3)
if (micro_mg1_0_tail_diag_use_native_impl) then
   do i = 1,ncol
      do k=top_lev,pver
         ncai(i,k)=dum2i(i,k)*rho(i,k)
         ncal(i,k)=dum2l(i,k)*rho(i,k)
      end do
   end do
else
   call micro_mg1_0_tail_activation_codon_wrap(ncol, pcols, pver, top_lev, dum2i, dum2l, rho, ncai, ncal)
end if


!redefine fice here....
if (micro_mg1_0_tail_diag_use_native_impl) then
   nfice(:,:)=0._r8
   do k=top_lev,pver
      do i=1,ncol
         dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
         dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
         dumfice=qsout(i,k) + qrout(i,k) + dumc(i,k) + dumi(i,k)

         if (dumfice.gt.qsmall.and.(qsout(i,k)+dumi(i,k).gt.qsmall)) then
            nfice(i,k)=(qsout(i,k) + dumi(i,k))/dumfice
         endif

         if (nfice(i,k).gt.1._r8) then
            nfice(i,k)=1._r8
         endif

      enddo
   enddo
else
   call micro_mg1_0_tail_fice_codon_wrap(ncol, pcols, pver, top_lev, deltat, qsmall, qc, qi, qctend, &
        qitend, qsout, qrout, dumc, dumi, nfice)
end if


end subroutine micro_mg_tend

subroutine micro_mg1_0_select_tend_impl()
  character(len=32) :: impl_name
  integer :: n, status

  if (micro_mg1_0_tend_impl_selected) return

  call get_environment_variable('MICRO_MG1_0_TEND_IMPL', value=impl_name, length=n, status=status)
  if (status == 0 .and. n > 0) then
     micro_mg1_0_tend_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     micro_mg1_0_tend_use_native_impl = .false.
  end if

  if (micro_mg1_0_tend_use_native_impl) then
     micro_mg1_0_init_use_native_impl = .true.
     micro_mg1_0_init_impl_selected = .true.
     micro_mg1_0_colzero_use_native_impl = .true.
     micro_mg1_0_colzero_impl_selected = .true.
     micro_mg1_0_tail_diag_use_native_impl = .true.
     micro_mg1_0_tail_diag_impl_selected = .true.
  end if

  if (masterproc) then
     if (micro_mg1_0_tend_use_native_impl) then
        write(iulog,*) 'micro_mg_tend implementation = native'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_TEND_PROOF_FILE', &
             'micro_mg_tend implementation = native')
     else
        write(iulog,*) 'micro_mg_tend implementation = codon'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_TEND_PROOF_FILE', &
             'micro_mg_tend implementation = codon')
     end if
  end if

  micro_mg1_0_tend_impl_selected = .true.
end subroutine micro_mg1_0_select_tend_impl

subroutine micro_mg1_0_log_tend_entered()
  if (micro_mg1_0_tend_logged) return
  micro_mg1_0_tend_logged = .true.
  if (masterproc) then
     write(iulog,*) 'micro_mg_tend direct = codon helper stages; native MG process core and sensitive diagnostics'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_TEND_PROOF_FILE', &
          'micro_mg_tend direct = codon helper stages; native MG process core and sensitive diagnostics')
  end if
end subroutine micro_mg1_0_log_tend_entered

subroutine micro_mg1_0_select_tail_diag_impl()
  character(len=32) :: impl_name
  integer :: n, status

  if (micro_mg1_0_tail_diag_impl_selected) return

  call get_environment_variable('MICRO_MG1_0_TAIL_DIAG_IMPL', value=impl_name, length=n, status=status)
  if (status == 0 .and. n > 0) then
     micro_mg1_0_tail_diag_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     micro_mg1_0_tail_diag_use_native_impl = .false.
  end if

  if (masterproc) then
     if (micro_mg1_0_tail_diag_use_native_impl) then
        write(iulog,*) 'micro_mg1_0_tail_diag implementation = native'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_TAIL_DIAG_PROOF_FILE', &
             'micro_mg1_0_tail_diag implementation = native')
     else
        write(iulog,*) 'micro_mg1_0_tail_diag implementation = codon'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_TAIL_DIAG_PROOF_FILE', &
             'micro_mg1_0_tail_diag implementation = codon')
     end if
  end if

  micro_mg1_0_tail_diag_impl_selected = .true.
end subroutine micro_mg1_0_select_tail_diag_impl

subroutine micro_mg1_0_tail_diag_log_entry()
  if (masterproc .and. .not. micro_mg1_0_tail_diag_logged) then
     write(iulog,*) 'micro_mg1_0_tail_diag entered (activation and nfice tail diagnostics direct = codon; ' // &
          'native fractional-power reflectivity/diameter diagnostics)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_TAIL_DIAG_PROOF_FILE', &
          'micro_mg1_0_tail_diag entered (activation and nfice tail diagnostics direct = codon; ' // &
          'native fractional-power reflectivity/diameter diagnostics)')
     micro_mg1_0_tail_diag_logged = .true.
  end if
end subroutine micro_mg1_0_tail_diag_log_entry

subroutine micro_mg1_0_select_init_impl()
  character(len=32) :: impl_name
  integer :: n, status

  if (micro_mg1_0_init_impl_selected) return

  call get_environment_variable('MICRO_MG1_0_INIT_IMPL', value=impl_name, length=n, status=status)
  if (status == 0 .and. n > 0) then
     micro_mg1_0_init_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     micro_mg1_0_init_use_native_impl = .false.
  end if

  if (masterproc) then
     if (micro_mg1_0_init_use_native_impl) then
        write(iulog,*) 'micro_mg1_0_init implementation = native'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_INIT_PROOF_FILE', &
             'micro_mg1_0_init implementation = native')
     else
        write(iulog,*) 'micro_mg1_0_init implementation = codon'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_INIT_PROOF_FILE', &
             'micro_mg1_0_init implementation = codon')
     end if
  end if

  micro_mg1_0_init_impl_selected = .true.
end subroutine micro_mg1_0_select_init_impl

subroutine micro_mg1_0_append_impl_proof(env_name, proof_line)
  character(len=*), intent(in) :: env_name
  character(len=*), intent(in) :: proof_line
  character(len=512) :: proof_file
  integer :: n, status, unitno

  if (.not. masterproc) return

  call get_environment_variable(env_name, value=proof_file, length=n, status=status)
  if (status /= 0 .or. n <= 0) return

  open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
  write(unitno,'(A)') trim(proof_line)
  close(unitno)
end subroutine micro_mg1_0_append_impl_proof

subroutine micro_mg1_0_log_init_scalars_entered()
  if (micro_mg1_0_init_scalars_logged) return
  micro_mg1_0_init_scalars_logged = .true.
  if (masterproc) then
     write(iulog,*) 'micro_mg1_0_init scalars entered (module constants direct = codon; gamma/atan/power native island)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_INIT_PROOF_FILE', &
          'micro_mg1_0_init scalars entered (module constants direct = codon; gamma/atan/power native island)')
  end if
end subroutine micro_mg1_0_log_init_scalars_entered

subroutine micro_mg1_0_init_fields_codon_wrap(ncol_local, pcols_local, pver_local, top_lev_local, mincld_local, &
     qn_local, tn_local, qc_local, qi_local, nc_local, ni_local, ncai_local, ncal_local, rercld_local, arcld_local, &
     pgamrad_local, lamcrad_local, deffi_local, qcsevap_local, qisevap_local, qvres_local, cmeiout_local, &
     vtrmc_local, vtrmi_local, qcsedten_local, qisedten_local, prao_local, prco_local, mnuccco_local, &
     mnuccto_local, msacwio_local, psacwso_local, bergso_local, bergo_local, melto_local, homoo_local, &
     qcreso_local, prcio_local, praio_local, qireso_local, mnuccro_local, pracso_local, meltsdt_local, &
     frzrdt_local, mnuccdo_local, rflx_local, sflx_local, effc_local, effc_fn_local, effi_local, preo_local, &
     prdso_local, frzro_local, meltso_local, wtfc_local, wtfi_local, wtprelat_local, wtpostlat_local, q_local, &
     t_local, t1_local, q1_local, qc1_local, qi1_local, nc1_local, ni1_local, tlat1_local, qvlat1_local, &
     qctend1_local, qitend1_local, nctend1_local, nitend1_local, qrout_local, qsout_local, nrout_local, &
     nsout_local, dsout_local, drout_local, reff_rain_local, reff_snow_local, nevapr_local, nevapr2_local, &
     evapsnow_local, prain_local, prodsnow_local, cmeout_local, am_evp_st_local, rainrt1_local, cldmax_local, &
     dum2l_local, dum2i_local, prect1_local, preci1_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: ncol_local, pcols_local, pver_local, top_lev_local
  real(r8), intent(in) :: mincld_local
  real(r8), target, intent(in) :: qn_local(pcols_local,pver_local), tn_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qc_local(pcols_local,pver_local), qi_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nc_local(pcols_local,pver_local), ni_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: ncai_local(pcols_local,pver_local), ncal_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: rercld_local(pcols_local,pver_local), arcld_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: pgamrad_local(pcols_local,pver_local), lamcrad_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: deffi_local(pcols_local,pver_local), qcsevap_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qisevap_local(pcols_local,pver_local), qvres_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: cmeiout_local(pcols_local,pver_local), vtrmc_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: vtrmi_local(pcols_local,pver_local), qcsedten_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qisedten_local(pcols_local,pver_local), prao_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prco_local(pcols_local,pver_local), mnuccco_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: mnuccto_local(pcols_local,pver_local), msacwio_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: psacwso_local(pcols_local,pver_local), bergso_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: bergo_local(pcols_local,pver_local), melto_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: homoo_local(pcols_local,pver_local), qcreso_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prcio_local(pcols_local,pver_local), praio_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qireso_local(pcols_local,pver_local), mnuccro_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: pracso_local(pcols_local,pver_local), meltsdt_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: frzrdt_local(pcols_local,pver_local), mnuccdo_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: rflx_local(pcols_local,pver_local+1), sflx_local(pcols_local,pver_local+1)
  real(r8), target, intent(inout) :: effc_local(pcols_local,pver_local), effc_fn_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: effi_local(pcols_local,pver_local), preo_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prdso_local(pcols_local,pver_local), frzro_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: meltso_local(pcols_local,pver_local), wtfc_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: wtfi_local(pcols_local,pver_local), wtprelat_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: wtpostlat_local(pcols_local,pver_local), q_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: t_local(pcols_local,pver_local), t1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: q1_local(pcols_local,pver_local), qc1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qi1_local(pcols_local,pver_local), nc1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: ni1_local(pcols_local,pver_local), tlat1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qvlat1_local(pcols_local,pver_local), qctend1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qitend1_local(pcols_local,pver_local), nctend1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nitend1_local(pcols_local,pver_local), qrout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qsout_local(pcols_local,pver_local), nrout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nsout_local(pcols_local,pver_local), dsout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: drout_local(pcols_local,pver_local), reff_rain_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: reff_snow_local(pcols_local,pver_local), nevapr_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nevapr2_local(pcols_local,pver_local), evapsnow_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prain_local(pcols_local,pver_local), prodsnow_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: cmeout_local(pcols_local,pver_local), am_evp_st_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: rainrt1_local(pcols_local,pver_local), cldmax_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: dum2l_local(pcols_local,pver_local), dum2i_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prect1_local(pcols_local), preci1_local(pcols_local)

  interface
     subroutine micro_mg1_0_init_fields_codon(ncol_c, pcols_c, pver_c, top_lev_c, mincld_c, qn_p, tn_p, &
          qc_p, qi_p, nc_p, ni_p, ncai_p, ncal_p, rercld_p, arcld_p, pgamrad_p, lamcrad_p, deffi_p, &
          qcsevap_p, qisevap_p, qvres_p, cmeiout_p, vtrmc_p, vtrmi_p, qcsedten_p, qisedten_p, prao_p, &
          prco_p, mnuccco_p, mnuccto_p, msacwio_p, psacwso_p, bergso_p, bergo_p, melto_p, homoo_p, &
          qcreso_p, prcio_p, praio_p, qireso_p, mnuccro_p, pracso_p, meltsdt_p, frzrdt_p, mnuccdo_p, &
          rflx_p, sflx_p, effc_p, effc_fn_p, effi_p, preo_p, prdso_p, frzro_p, meltso_p, wtfc_p, &
          wtfi_p, wtprelat_p, wtpostlat_p, q_p, t_p, t1_p, q1_p, qc1_p, qi1_p, nc1_p, ni1_p, &
          tlat1_p, qvlat1_p, qctend1_p, qitend1_p, nctend1_p, nitend1_p, qrout_p, qsout_p, nrout_p, &
          nsout_p, dsout_p, drout_p, reff_rain_p, reff_snow_p, nevapr_p, nevapr2_p, evapsnow_p, &
          prain_p, prodsnow_p, cmeout_p, am_evp_st_p, rainrt1_p, cldmax_p, dum2l_p, dum2i_p, &
          prect1_p, preci1_p) bind(c, name="micro_mg1_0_init_fields_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       real(c_double), value :: mincld_c
       type(c_ptr), value :: qn_p, tn_p, qc_p, qi_p, nc_p, ni_p, ncai_p, ncal_p, rercld_p, arcld_p
       type(c_ptr), value :: pgamrad_p, lamcrad_p, deffi_p, qcsevap_p, qisevap_p, qvres_p, cmeiout_p
       type(c_ptr), value :: vtrmc_p, vtrmi_p, qcsedten_p, qisedten_p, prao_p, prco_p, mnuccco_p
       type(c_ptr), value :: mnuccto_p, msacwio_p, psacwso_p, bergso_p, bergo_p, melto_p, homoo_p
       type(c_ptr), value :: qcreso_p, prcio_p, praio_p, qireso_p, mnuccro_p, pracso_p, meltsdt_p
       type(c_ptr), value :: frzrdt_p, mnuccdo_p, rflx_p, sflx_p, effc_p, effc_fn_p, effi_p, preo_p
       type(c_ptr), value :: prdso_p, frzro_p, meltso_p, wtfc_p, wtfi_p, wtprelat_p, wtpostlat_p
       type(c_ptr), value :: q_p, t_p, t1_p, q1_p, qc1_p, qi1_p, nc1_p, ni1_p, tlat1_p, qvlat1_p
       type(c_ptr), value :: qctend1_p, qitend1_p, nctend1_p, nitend1_p, qrout_p, qsout_p, nrout_p
       type(c_ptr), value :: nsout_p, dsout_p, drout_p, reff_rain_p, reff_snow_p, nevapr_p, nevapr2_p
       type(c_ptr), value :: evapsnow_p, prain_p, prodsnow_p, cmeout_p, am_evp_st_p, rainrt1_p, cldmax_p
       type(c_ptr), value :: dum2l_p, dum2i_p, prect1_p, preci1_p
     end subroutine micro_mg1_0_init_fields_codon
  end interface

  call micro_mg1_0_append_impl_proof('MICRO_MG1_0_INIT_PROOF_FILE', &
       'micro_mg1_0_init_fields_codon_wrap entered (zero/copy initialization = codon)')

  if (masterproc .and. .not. micro_mg1_0_init_wrapper_logged) then
     write(iulog,*) 'micro_mg1_0_init_fields_codon_wrap entered (zero/copy initialization = codon)'
     micro_mg1_0_init_wrapper_logged = .true.
  end if

  call micro_mg1_0_init_fields_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), real(mincld_local, c_double), &
       c_loc(qn_local), c_loc(tn_local), c_loc(qc_local), c_loc(qi_local), c_loc(nc_local), c_loc(ni_local), &
       c_loc(ncai_local), c_loc(ncal_local), c_loc(rercld_local), c_loc(arcld_local), c_loc(pgamrad_local), &
       c_loc(lamcrad_local), c_loc(deffi_local), c_loc(qcsevap_local), c_loc(qisevap_local), c_loc(qvres_local), &
       c_loc(cmeiout_local), c_loc(vtrmc_local), c_loc(vtrmi_local), c_loc(qcsedten_local), c_loc(qisedten_local), &
       c_loc(prao_local), c_loc(prco_local), c_loc(mnuccco_local), c_loc(mnuccto_local), c_loc(msacwio_local), &
       c_loc(psacwso_local), c_loc(bergso_local), c_loc(bergo_local), c_loc(melto_local), c_loc(homoo_local), &
       c_loc(qcreso_local), c_loc(prcio_local), c_loc(praio_local), c_loc(qireso_local), c_loc(mnuccro_local), &
       c_loc(pracso_local), c_loc(meltsdt_local), c_loc(frzrdt_local), c_loc(mnuccdo_local), c_loc(rflx_local), &
       c_loc(sflx_local), c_loc(effc_local), c_loc(effc_fn_local), c_loc(effi_local), c_loc(preo_local), &
       c_loc(prdso_local), c_loc(frzro_local), c_loc(meltso_local), c_loc(wtfc_local), c_loc(wtfi_local), &
       c_loc(wtprelat_local), c_loc(wtpostlat_local), c_loc(q_local), c_loc(t_local), c_loc(t1_local), &
       c_loc(q1_local), c_loc(qc1_local), c_loc(qi1_local), c_loc(nc1_local), c_loc(ni1_local), &
       c_loc(tlat1_local), c_loc(qvlat1_local), c_loc(qctend1_local), c_loc(qitend1_local), &
       c_loc(nctend1_local), c_loc(nitend1_local), c_loc(qrout_local), c_loc(qsout_local), &
       c_loc(nrout_local), c_loc(nsout_local), c_loc(dsout_local), c_loc(drout_local), &
       c_loc(reff_rain_local), c_loc(reff_snow_local), c_loc(nevapr_local), c_loc(nevapr2_local), &
       c_loc(evapsnow_local), c_loc(prain_local), c_loc(prodsnow_local), c_loc(cmeout_local), &
       c_loc(am_evp_st_local), c_loc(rainrt1_local), c_loc(cldmax_local), c_loc(dum2l_local), &
       c_loc(dum2i_local), c_loc(prect1_local), c_loc(preci1_local))
end subroutine micro_mg1_0_init_fields_codon_wrap

subroutine micro_mg1_0_select_colzero_impl()
  character(len=32) :: impl_name
  integer :: n, status

  if (micro_mg1_0_colzero_impl_selected) return

  call get_environment_variable('MICRO_MG1_0_COLZERO_IMPL', value=impl_name, length=n, status=status)
  if (status == 0 .and. n > 0) then
     micro_mg1_0_colzero_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     micro_mg1_0_colzero_use_native_impl = .false.
  end if

  if (masterproc) then
     if (micro_mg1_0_colzero_use_native_impl) then
        write(iulog,*) 'micro_mg1_0_colzero implementation = native'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
             'micro_mg1_0_colzero implementation = native')
     else
        write(iulog,*) 'micro_mg1_0_colzero implementation = codon'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
             'micro_mg1_0_colzero implementation = codon')
     end if
  end if

  micro_mg1_0_colzero_impl_selected = .true.
end subroutine micro_mg1_0_select_colzero_impl

subroutine micro_mg1_0_colzero_log_entry()
  if (masterproc .and. .not. micro_mg1_0_colzero_wrapper_logged) then
     write(iulog,*) 'micro_mg1_0_colzero_wrap entered (column/flux/ltrue initialization = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_colzero_wrap entered (column/flux/ltrue initialization = codon)')
     micro_mg1_0_colzero_wrapper_logged = .true.
  end if
end subroutine micro_mg1_0_colzero_log_entry

subroutine micro_mg1_0_rate1ord_log_entry()
  if (masterproc .and. .not. micro_mg1_0_rate1ord_logged) then
     write(iulog,*) 'micro_mg1_0_rate1ord entered (cw-to-precip rate bookkeeping = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_rate1ord entered (cw-to-precip rate bookkeeping = codon)')
     micro_mg1_0_rate1ord_logged = .true.
  end if
end subroutine micro_mg1_0_rate1ord_log_entry

subroutine micro_mg1_0_substep_setup_log_entry()
  if (masterproc .and. .not. micro_mg1_0_substep_setup_logged) then
     write(iulog,*) 'micro_mg1_0_substep_setup entered (substep precip/cloud scratch setup = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_substep_setup entered (substep precip/cloud scratch setup = codon)')
     micro_mg1_0_substep_setup_logged = .true.
  end if
end subroutine micro_mg1_0_substep_setup_log_entry

subroutine micro_mg1_0_substep_accum_log_entry()
  if (masterproc .and. .not. micro_mg1_0_substep_accum_logged) then
     write(iulog,*) 'micro_mg1_0_substep_accum entered (substep tendency/state/flux accumulation = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_substep_accum entered (substep tendency/state/flux accumulation = codon)')
     micro_mg1_0_substep_accum_logged = .true.
  end if
end subroutine micro_mg1_0_substep_accum_log_entry

subroutine micro_mg1_0_incloud_activation_log_entry()
  if (masterproc .and. .not. micro_mg1_0_incloud_activation_logged) then
     write(iulog,*) 'micro_mg1_0_incloud_activation entered (in-cloud state and droplet activation prep = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_incloud_activation entered (in-cloud state and droplet activation prep = codon)')
     micro_mg1_0_incloud_activation_logged = .true.
  end if
end subroutine micro_mg1_0_incloud_activation_log_entry

subroutine micro_mg1_0_conservation_limiter_log_entry()
  if (masterproc .and. .not. micro_mg1_0_conservation_limiter_logged) then
     write(iulog,*) 'micro_mg1_0_conservation_limiter entered (process-rate conservation limiter = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_conservation_limiter entered (process-rate conservation limiter = codon)')
     micro_mg1_0_conservation_limiter_logged = .true.
  end if
end subroutine micro_mg1_0_conservation_limiter_log_entry

subroutine micro_mg1_0_process_output_log_entry()
  if (masterproc .and. .not. micro_mg1_0_process_output_logged) then
     write(iulog,*) 'micro_mg1_0_process_output entered (process/output accumulation = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_process_output entered (process/output accumulation = codon)')
     micro_mg1_0_process_output_logged = .true.
  end if
end subroutine micro_mg1_0_process_output_log_entry

subroutine micro_mg1_0_post_iter_avg_log_entry()
  if (masterproc .and. .not. micro_mg1_0_post_iter_avg_logged) then
     write(iulog,*) 'micro_mg1_0_post_iter_avg entered (post-iteration averaging/state restore = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_post_iter_avg entered (post-iteration averaging/state restore = codon)')
     micro_mg1_0_post_iter_avg_logged = .true.
  end if
end subroutine micro_mg1_0_post_iter_avg_log_entry

subroutine micro_mg1_0_phase_change_log_entry()
  if (masterproc .and. .not. micro_mg1_0_phase_change_logged) then
     write(iulog,*) 'micro_mg1_0_phase_change entered (cloud ice melting/freezing = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_phase_change entered (cloud ice melting/freezing = codon)')
     micro_mg1_0_phase_change_logged = .true.
  end if
end subroutine micro_mg1_0_phase_change_log_entry

subroutine micro_mg1_0_number_cleanup_log_entry()
  if (masterproc .and. .not. micro_mg1_0_number_cleanup_logged) then
     write(iulog,*) 'micro_mg1_0_number_cleanup entered (post-diagnostic number cleanup = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_number_cleanup entered (post-diagnostic number cleanup = codon)')
     micro_mg1_0_number_cleanup_logged = .true.
  end if
end subroutine micro_mg1_0_number_cleanup_log_entry

subroutine micro_mg1_0_reflectivity_flags_log_entry()
  if (masterproc .and. .not. micro_mg1_0_reflectivity_flags_logged) then
     write(iulog,*) 'micro_mg1_0_reflectivity_flags entered (reflectivity flags/bounds = codon)'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_COLZERO_PROOF_FILE', &
          'micro_mg1_0_reflectivity_flags entered (reflectivity flags/bounds = codon)')
     micro_mg1_0_reflectivity_flags_logged = .true.
  end if
end subroutine micro_mg1_0_reflectivity_flags_log_entry

subroutine micro_mg1_0_flux_ltrue_init_codon_wrap(ncol_local, pcols_local, pver_local, top_lev_local, &
     qsmall_local, rflx1_local, sflx1_local, rflx_local, sflx_local, qc_local, qi_local, cmei_local, &
     ltrue_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: ncol_local, pcols_local, pver_local, top_lev_local
  real(r8), intent(in) :: qsmall_local
  real(r8), target, intent(inout) :: rflx1_local(pcols_local,pver_local+1)
  real(r8), target, intent(inout) :: sflx1_local(pcols_local,pver_local+1)
  real(r8), target, intent(inout) :: rflx_local(pcols_local,pver_local+1)
  real(r8), target, intent(inout) :: sflx_local(pcols_local,pver_local+1)
  real(r8), target, intent(in) :: qc_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qi_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: cmei_local(pcols_local,pver_local)
  integer, target, intent(inout) :: ltrue_local(pcols_local)

  interface
     subroutine micro_mg1_0_flux_ltrue_init_codon(ncol_c, pcols_c, pver_c, top_lev_c, qsmall_c, &
          rflx1_p, sflx1_p, rflx_p, sflx_p, qc_p, qi_p, cmei_p, ltrue_p) &
          bind(c, name="micro_mg1_0_flux_ltrue_init_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       real(c_double), value :: qsmall_c
       type(c_ptr), value :: rflx1_p, sflx1_p, rflx_p, sflx_p
       type(c_ptr), value :: qc_p, qi_p, cmei_p, ltrue_p
     end subroutine micro_mg1_0_flux_ltrue_init_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_flux_ltrue_init_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), real(qsmall_local, c_double), &
       c_loc(rflx1_local), c_loc(sflx1_local), c_loc(rflx_local), c_loc(sflx_local), c_loc(qc_local), &
       c_loc(qi_local), c_loc(cmei_local), c_loc(ltrue_local))
end subroutine micro_mg1_0_flux_ltrue_init_codon_wrap

subroutine micro_mg1_0_rate1ord_zero_codon_wrap(ncol_local, pcols_local, pver_local, rate1ord_cw2pr_st_local)
  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  integer, intent(in) :: ncol_local, pcols_local, pver_local
  real(r8), target, intent(inout) :: rate1ord_cw2pr_st_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_rate1ord_zero_codon(ncol_c, pcols_c, pver_c, rate1ord_cw2pr_st_p) &
          bind(c, name="micro_mg1_0_rate1ord_zero_codon")
       import c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
       type(c_ptr), value :: rate1ord_cw2pr_st_p
     end subroutine micro_mg1_0_rate1ord_zero_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_rate1ord_log_entry()
  call micro_mg1_0_rate1ord_zero_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), c_loc(rate1ord_cw2pr_st_local))
end subroutine micro_mg1_0_rate1ord_zero_codon_wrap

subroutine micro_mg1_0_rate1ord_column_codon_wrap(i_local, pcols_local, pver_local, top_lev_local, &
     qcsinksum_rate1ord_local, qcsum_rate1ord_local, rate1ord_cw2pr_st_local)
  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, pcols_local, pver_local, top_lev_local
  real(r8), target, intent(in) :: qcsinksum_rate1ord_local(pver_local)
  real(r8), target, intent(in) :: qcsum_rate1ord_local(pver_local)
  real(r8), target, intent(inout) :: rate1ord_cw2pr_st_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_rate1ord_column_codon(i_c, pcols_c, pver_c, top_lev_c, &
          qcsinksum_rate1ord_p, qcsum_rate1ord_p, rate1ord_cw2pr_st_p) &
          bind(c, name="micro_mg1_0_rate1ord_column_codon")
       import c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: qcsinksum_rate1ord_p, qcsum_rate1ord_p, rate1ord_cw2pr_st_p
     end subroutine micro_mg1_0_rate1ord_column_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_rate1ord_log_entry()
  call micro_mg1_0_rate1ord_column_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), c_loc(qcsinksum_rate1ord_local), &
       c_loc(qcsum_rate1ord_local), c_loc(rate1ord_cw2pr_st_local))
end subroutine micro_mg1_0_rate1ord_column_codon_wrap

subroutine micro_mg1_0_no_cloud_zero_column_codon_wrap(i_local, pcols_local, pver_local, &
     tlat_local, qvlat_local, qctend_local, qitend_local, qnitend_local, qrtend_local, &
     nctend_local, nitend_local, nrtend_local, nstend_local, prect_local, preci_local, &
     qniic_local, qric_local, nsic_local, nric_local, rainrt_local, qrtend_copy_local, &
     qnitend_copy_local)
  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, pcols_local, pver_local
  real(r8), target, intent(inout) :: tlat_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qvlat_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qctend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qnitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qrtend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nctend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nrtend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nstend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prect_local(pcols_local), preci_local(pcols_local)
  real(r8), target, intent(inout) :: qniic_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qric_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nsic_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nric_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: rainrt_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qrtend_copy_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qnitend_copy_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_no_cloud_zero_column_codon(i_c, pcols_c, pver_c, &
          tlat_p, qvlat_p, qctend_p, qitend_p, qnitend_p, qrtend_p, &
          nctend_p, nitend_p, nrtend_p, nstend_p, prect_p, preci_p, &
          qniic_p, qric_p, nsic_p, nric_p, rainrt_p, qrtend_copy_p, &
          qnitend_copy_p) bind(c, name="micro_mg1_0_no_cloud_zero_column_codon")
       import c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, pcols_c, pver_c
       type(c_ptr), value :: tlat_p, qvlat_p, qctend_p, qitend_p, qnitend_p, qrtend_p
       type(c_ptr), value :: nctend_p, nitend_p, nrtend_p, nstend_p, prect_p, preci_p
       type(c_ptr), value :: qniic_p, qric_p, nsic_p, nric_p, rainrt_p, qrtend_copy_p
       type(c_ptr), value :: qnitend_copy_p
     end subroutine micro_mg1_0_no_cloud_zero_column_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_no_cloud_zero_column_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), c_loc(tlat_local), c_loc(qvlat_local), c_loc(qctend_local), &
       c_loc(qitend_local), c_loc(qnitend_local), c_loc(qrtend_local), c_loc(nctend_local), &
       c_loc(nitend_local), c_loc(nrtend_local), c_loc(nstend_local), c_loc(prect_local), &
       c_loc(preci_local), c_loc(qniic_local), c_loc(qric_local), c_loc(nsic_local), &
       c_loc(nric_local), c_loc(rainrt_local), c_loc(qrtend_copy_local), c_loc(qnitend_copy_local))
end subroutine micro_mg1_0_no_cloud_zero_column_codon_wrap

subroutine micro_mg1_0_substep_zero_column_codon_wrap(i_local, pcols_local, pver_local, &
     tlat_local, qvlat_local, qctend_local, qitend_local, qnitend_local, qrtend_local, &
     nctend_local, nitend_local, nrtend_local, nstend_local, qniic_local, qric_local, &
     nsic_local, nric_local, rainrt_local, qrtend_copy_local, qnitend_copy_local)
  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, pcols_local, pver_local
  real(r8), target, intent(inout) :: tlat_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qvlat_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qctend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qnitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qrtend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nctend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nrtend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nstend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qniic_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qric_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nsic_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nric_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: rainrt_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qrtend_copy_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qnitend_copy_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_substep_zero_column_codon(i_c, pcols_c, pver_c, &
          tlat_p, qvlat_p, qctend_p, qitend_p, qnitend_p, qrtend_p, &
          nctend_p, nitend_p, nrtend_p, nstend_p, qniic_p, qric_p, &
          nsic_p, nric_p, rainrt_p, qrtend_copy_p, qnitend_copy_p) &
          bind(c, name="micro_mg1_0_substep_zero_column_codon")
       import c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, pcols_c, pver_c
       type(c_ptr), value :: tlat_p, qvlat_p, qctend_p, qitend_p, qnitend_p, qrtend_p
       type(c_ptr), value :: nctend_p, nitend_p, nrtend_p, nstend_p, qniic_p, qric_p
       type(c_ptr), value :: nsic_p, nric_p, rainrt_p, qrtend_copy_p, qnitend_copy_p
     end subroutine micro_mg1_0_substep_zero_column_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_substep_zero_column_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), c_loc(tlat_local), c_loc(qvlat_local), c_loc(qctend_local), &
       c_loc(qitend_local), c_loc(qnitend_local), c_loc(qrtend_local), c_loc(nctend_local), &
       c_loc(nitend_local), c_loc(nrtend_local), c_loc(nstend_local), c_loc(qniic_local), &
       c_loc(qric_local), c_loc(nsic_local), c_loc(nric_local), c_loc(rainrt_local), &
       c_loc(qrtend_copy_local), c_loc(qnitend_copy_local))
end subroutine micro_mg1_0_substep_zero_column_codon_wrap

subroutine micro_mg1_0_substep_setup_column_codon_wrap(i_local, pcols_local, pver_local, &
     top_lev_local, qsmall_local, mincld_local, qc_local, qi_local, ni_local, cldm_local, &
     cmei_local, cwml_local, cwmi_local, ums_local, uns_local, umr_local, unr_local, &
     nsubi_local, nsubc_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, pcols_local, pver_local, top_lev_local
  real(r8), intent(in) :: qsmall_local, mincld_local
  real(r8), target, intent(in) :: qc_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qi_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: ni_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: cldm_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: cmei_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: cwml_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: cwmi_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: ums_local(pver_local)
  real(r8), target, intent(inout) :: uns_local(pver_local)
  real(r8), target, intent(inout) :: umr_local(pver_local)
  real(r8), target, intent(inout) :: unr_local(pver_local)
  real(r8), target, intent(inout) :: nsubi_local(pver_local)
  real(r8), target, intent(inout) :: nsubc_local(pver_local)

  interface
     subroutine micro_mg1_0_substep_setup_column_codon(i_c, pcols_c, pver_c, top_lev_c, &
          qsmall_c, mincld_c, qc_p, qi_p, ni_p, cldm_p, cmei_p, cwml_p, cwmi_p, &
          ums_p, uns_p, umr_p, unr_p, nsubi_p, nsubc_p) &
          bind(c, name="micro_mg1_0_substep_setup_column_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c
       real(c_double), value :: qsmall_c, mincld_c
       type(c_ptr), value :: qc_p, qi_p, ni_p, cldm_p, cmei_p, cwml_p, cwmi_p
       type(c_ptr), value :: ums_p, uns_p, umr_p, unr_p, nsubi_p, nsubc_p
     end subroutine micro_mg1_0_substep_setup_column_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_substep_setup_log_entry()
  call micro_mg1_0_substep_setup_column_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), &
       real(qsmall_local, c_double), real(mincld_local, c_double), &
       c_loc(qc_local), c_loc(qi_local), c_loc(ni_local), c_loc(cldm_local), c_loc(cmei_local), &
       c_loc(cwml_local), c_loc(cwmi_local), c_loc(ums_local), c_loc(uns_local), c_loc(umr_local), &
       c_loc(unr_local), c_loc(nsubi_local), c_loc(nsubc_local))
end subroutine micro_mg1_0_substep_setup_column_codon_wrap

real(r8) function micro_mg1_0_incloud_activation_prep_codon_wrap(i_local, k_local, pcols_local, pver_local, &
     deltat_local, qsmall_local, omsm_local, cdnl_local, do_cldice_local, cwml_local, cwmi_local, &
     lcldm_local, icldm_local, nc_local, ni_local, qc_local, qi_local, berg_local, cmei_local, cmeout_local, &
     npccnin_local, rho_local, qcic_local, qiic_local, ncic_local, niic_local, npccn_local, dum2l_local) &
     result(ncmax_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, k_local, pcols_local, pver_local
  real(r8), intent(in) :: deltat_local, qsmall_local, omsm_local, cdnl_local
  logical, intent(in) :: do_cldice_local
  real(r8), target, intent(in) :: cwml_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: cwmi_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: lcldm_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: icldm_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nc_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: ni_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qc_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qi_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: berg_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: cmei_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: cmeout_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: npccnin_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: rho_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qcic_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qiic_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: ncic_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: niic_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: npccn_local(pver_local)
  real(r8), target, intent(inout) :: dum2l_local(pcols_local,pver_local)

  interface
     function micro_mg1_0_incloud_activation_prep_codon(i_c, k_c, pcols_c, pver_c, deltat_c, qsmall_c, &
          omsm_c, cdnl_c, do_cldice_c, cwml_p, cwmi_p, lcldm_p, icldm_p, nc_p, ni_p, qc_p, qi_p, &
          berg_p, cmei_p, cmeout_p, npccnin_p, rho_p, qcic_p, qiic_p, ncic_p, niic_p, npccn_p, &
          dum2l_p) result(ncmax_c) bind(c, name="micro_mg1_0_incloud_activation_prep_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, do_cldice_c
       real(c_double), value :: deltat_c, qsmall_c, omsm_c, cdnl_c
       type(c_ptr), value :: cwml_p, cwmi_p, lcldm_p, icldm_p, nc_p, ni_p, qc_p, qi_p
       type(c_ptr), value :: berg_p, cmei_p, cmeout_p, npccnin_p, rho_p
       type(c_ptr), value :: qcic_p, qiic_p, ncic_p, niic_p, npccn_p, dum2l_p
       real(c_double) :: ncmax_c
     end function micro_mg1_0_incloud_activation_prep_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_incloud_activation_log_entry()
  ncmax_local = real(micro_mg1_0_incloud_activation_prep_codon(int(i_local, c_int64_t), &
       int(k_local, c_int64_t), int(pcols_local, c_int64_t), int(pver_local, c_int64_t), &
       real(deltat_local, c_double), real(qsmall_local, c_double), real(omsm_local, c_double), &
       real(cdnl_local, c_double), merge(1_c_int64_t, 0_c_int64_t, do_cldice_local), c_loc(cwml_local), &
       c_loc(cwmi_local), c_loc(lcldm_local), c_loc(icldm_local), c_loc(nc_local), c_loc(ni_local), &
       c_loc(qc_local), c_loc(qi_local), c_loc(berg_local), c_loc(cmei_local), c_loc(cmeout_local), &
       c_loc(npccnin_local), c_loc(rho_local), c_loc(qcic_local), c_loc(qiic_local), c_loc(ncic_local), &
       c_loc(niic_local), c_loc(npccn_local), c_loc(dum2l_local)), r8)
end function micro_mg1_0_incloud_activation_prep_codon_wrap

subroutine micro_mg1_0_conservation_limiter_codon_wrap(i_local, k_local, pcols_local, pver_local, &
     deltat_local, omsm_local, qsmall_local, qce_local, nce_local, qie_local, nie_local, qrtot_local, &
     nrtot_local, qstot_local, nstot_local, do_cldice_local, use_hetfrz_classnuc_local, lcldm_local, &
     icldm_local, cldmax_local, dz_local, rho_local, prc_local, pra_local, mnuccc_local, mnucct_local, &
     msacwi_local, psacws_local, bergs_local, nprc1_local, npra_local, nnuccc_local, nnucct_local, &
     npsacws_local, nsubc_local, prci_local, prai_local, mnudep_local, nprci_local, nprai_local, &
     nsubi_local, nnudep_local, nsacwi_local, mnuccr_local, pre_local, pracs_local, nsubr_local, &
     npracs_local, nnuccr_local, nragg_local, prds_local, nsubs_local, nsagg_local, nprc_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, k_local, pcols_local, pver_local
  real(r8), intent(in) :: deltat_local, omsm_local, qsmall_local
  real(r8), intent(in) :: qce_local, nce_local, qie_local, nie_local
  real(r8), intent(in) :: qrtot_local, nrtot_local, qstot_local, nstot_local
  logical, intent(in) :: do_cldice_local, use_hetfrz_classnuc_local
  real(r8), target, intent(in) :: lcldm_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: icldm_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: cldmax_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: dz_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: rho_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prc_local(pver_local)
  real(r8), target, intent(inout) :: pra_local(pver_local)
  real(r8), target, intent(inout) :: mnuccc_local(pver_local)
  real(r8), target, intent(inout) :: mnucct_local(pver_local)
  real(r8), target, intent(inout) :: msacwi_local(pver_local)
  real(r8), target, intent(inout) :: psacws_local(pver_local)
  real(r8), target, intent(inout) :: bergs_local(pver_local)
  real(r8), target, intent(inout) :: nprc1_local(pver_local)
  real(r8), target, intent(inout) :: npra_local(pver_local)
  real(r8), target, intent(inout) :: nnuccc_local(pver_local)
  real(r8), target, intent(inout) :: nnucct_local(pver_local)
  real(r8), target, intent(inout) :: npsacws_local(pver_local)
  real(r8), target, intent(inout) :: nsubc_local(pver_local)
  real(r8), target, intent(inout) :: prci_local(pver_local)
  real(r8), target, intent(inout) :: prai_local(pver_local)
  real(r8), target, intent(in) :: mnudep_local(pver_local)
  real(r8), target, intent(inout) :: nprci_local(pver_local)
  real(r8), target, intent(inout) :: nprai_local(pver_local)
  real(r8), target, intent(inout) :: nsubi_local(pver_local)
  real(r8), target, intent(in) :: nnudep_local(pver_local)
  real(r8), target, intent(in) :: nsacwi_local(pver_local)
  real(r8), target, intent(inout) :: mnuccr_local(pver_local)
  real(r8), target, intent(inout) :: pre_local(pver_local)
  real(r8), target, intent(inout) :: pracs_local(pver_local)
  real(r8), target, intent(inout) :: nsubr_local(pver_local)
  real(r8), target, intent(inout) :: npracs_local(pver_local)
  real(r8), target, intent(inout) :: nnuccr_local(pver_local)
  real(r8), target, intent(inout) :: nragg_local(pver_local)
  real(r8), target, intent(inout) :: prds_local(pver_local)
  real(r8), target, intent(inout) :: nsubs_local(pver_local)
  real(r8), target, intent(inout) :: nsagg_local(pver_local)
  real(r8), target, intent(in) :: nprc_local(pver_local)

  interface
     subroutine micro_mg1_0_conservation_limiter_codon(i_c, k_c, pcols_c, pver_c, deltat_c, &
          omsm_c, qsmall_c, qce_c, nce_c, qie_c, nie_c, qrtot_c, nrtot_c, qstot_c, nstot_c, &
          do_cldice_c, use_hetfrz_classnuc_c, lcldm_p, icldm_p, cldmax_p, dz_p, rho_p, prc_p, &
          pra_p, mnuccc_p, mnucct_p, msacwi_p, psacws_p, bergs_p, nprc1_p, npra_p, nnuccc_p, &
          nnucct_p, npsacws_p, nsubc_p, prci_p, prai_p, mnudep_p, nprci_p, nprai_p, nsubi_p, &
          nnudep_p, nsacwi_p, mnuccr_p, pre_p, pracs_p, nsubr_p, npracs_p, nnuccr_p, nragg_p, &
          prds_p, nsubs_p, nsagg_p, nprc_p) bind(c, name="micro_mg1_0_conservation_limiter_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, do_cldice_c, use_hetfrz_classnuc_c
       real(c_double), value :: deltat_c, omsm_c, qsmall_c, qce_c, nce_c, qie_c, nie_c
       real(c_double), value :: qrtot_c, nrtot_c, qstot_c, nstot_c
       type(c_ptr), value :: lcldm_p, icldm_p, cldmax_p, dz_p, rho_p
       type(c_ptr), value :: prc_p, pra_p, mnuccc_p, mnucct_p, msacwi_p, psacws_p, bergs_p
       type(c_ptr), value :: nprc1_p, npra_p, nnuccc_p, nnucct_p, npsacws_p, nsubc_p
       type(c_ptr), value :: prci_p, prai_p, mnudep_p, nprci_p, nprai_p, nsubi_p
       type(c_ptr), value :: nnudep_p, nsacwi_p, mnuccr_p, pre_p, pracs_p, nsubr_p, npracs_p
       type(c_ptr), value :: nnuccr_p, nragg_p, prds_p, nsubs_p, nsagg_p, nprc_p
     end subroutine micro_mg1_0_conservation_limiter_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_conservation_limiter_log_entry()
  call micro_mg1_0_conservation_limiter_codon(int(i_local, c_int64_t), int(k_local, c_int64_t), &
       int(pcols_local, c_int64_t), int(pver_local, c_int64_t), real(deltat_local, c_double), &
       real(omsm_local, c_double), real(qsmall_local, c_double), real(qce_local, c_double), &
       real(nce_local, c_double), real(qie_local, c_double), real(nie_local, c_double), &
       real(qrtot_local, c_double), real(nrtot_local, c_double), real(qstot_local, c_double), &
       real(nstot_local, c_double), merge(1_c_int64_t, 0_c_int64_t, do_cldice_local), &
       merge(1_c_int64_t, 0_c_int64_t, use_hetfrz_classnuc_local), c_loc(lcldm_local), &
       c_loc(icldm_local), c_loc(cldmax_local), c_loc(dz_local), c_loc(rho_local), &
       c_loc(prc_local), c_loc(pra_local), c_loc(mnuccc_local), c_loc(mnucct_local), &
       c_loc(msacwi_local), c_loc(psacws_local), c_loc(bergs_local), c_loc(nprc1_local), &
       c_loc(npra_local), c_loc(nnuccc_local), c_loc(nnucct_local), c_loc(npsacws_local), &
       c_loc(nsubc_local), c_loc(prci_local), c_loc(prai_local), c_loc(mnudep_local), &
       c_loc(nprci_local), c_loc(nprai_local), c_loc(nsubi_local), c_loc(nnudep_local), &
       c_loc(nsacwi_local), c_loc(mnuccr_local), c_loc(pre_local), c_loc(pracs_local), &
       c_loc(nsubr_local), c_loc(npracs_local), c_loc(nnuccr_local), c_loc(nragg_local), &
       c_loc(prds_local), c_loc(nsubs_local), c_loc(nsagg_local), c_loc(nprc_local))
end subroutine micro_mg1_0_conservation_limiter_codon_wrap

subroutine micro_mg1_0_process_output_accum_codon_wrap(i_local, k_local, pcols_local, pver_local, &
     qrtend_local, qnitend_local, qrtend_copy_local, qnitend_copy_local, cmei_local, cmeiout_local, &
     prds_local, pre_local, cldmax_local, evapsnow_local, nevapr_local, nevapr2_local, pra_local, &
     prc_local, lcldm_local, pracs_local, mnuccr_local, prain_local, prai_local, prci_local, &
     icldm_local, psacws_local, bergs_local, prodsnow_local, qcsinksum_rate1ord_local, &
     qcsum_rate1ord_local, qc_local, prao_local, prco_local, mnuccc_local, mnucct_local, &
     mnuccd_local, msacwi_local, mnuccco_local, mnuccto_local, mnuccdo_local, msacwio_local, &
     psacwso_local, bergso_local, berg_local, bergo_local, prcio_local, praio_local, &
     mnuccro_local, pracso_local, preo_local, prdso_local)
  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, k_local, pcols_local, pver_local
  real(r8), target, intent(in) :: qrtend_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qnitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qrtend_copy_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qnitend_copy_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: cmei_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: cmeiout_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: cldmax_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: evapsnow_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nevapr_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nevapr2_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: lcldm_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prain_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: icldm_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prodsnow_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qcsinksum_rate1ord_local(pver_local)
  real(r8), target, intent(inout) :: qcsum_rate1ord_local(pver_local)
  real(r8), target, intent(in) :: qc_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prao_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prco_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: mnuccco_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: mnuccto_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: mnuccdo_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: msacwio_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: psacwso_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: bergso_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: berg_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: bergo_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prcio_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: praio_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: mnuccro_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: pracso_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: preo_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prdso_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: prds_local(pver_local), pre_local(pver_local)
  real(r8), target, intent(in) :: pra_local(pver_local), prc_local(pver_local)
  real(r8), target, intent(in) :: pracs_local(pver_local), mnuccr_local(pver_local)
  real(r8), target, intent(in) :: prai_local(pver_local), prci_local(pver_local)
  real(r8), target, intent(in) :: psacws_local(pver_local), bergs_local(pver_local)
  real(r8), target, intent(in) :: mnuccc_local(pver_local), mnucct_local(pver_local)
  real(r8), target, intent(in) :: mnuccd_local(pver_local), msacwi_local(pver_local)

  interface
     subroutine micro_mg1_0_process_output_accum_codon(i_c, k_c, pcols_c, pver_c, qrtend_p, &
          qnitend_p, qrtend_copy_p, qnitend_copy_p, cmei_p, cmeiout_p, prds_p, pre_p, &
          cldmax_p, evapsnow_p, nevapr_p, nevapr2_p, pra_p, prc_p, lcldm_p, pracs_p, &
          mnuccr_p, prain_p, prai_p, prci_p, icldm_p, psacws_p, bergs_p, prodsnow_p, &
          qcsinksum_rate1ord_p, qcsum_rate1ord_p, qc_p, prao_p, prco_p, mnuccc_p, &
          mnucct_p, mnuccd_p, msacwi_p, mnuccco_p, mnuccto_p, mnuccdo_p, msacwio_p, &
          psacwso_p, bergso_p, berg_p, bergo_p, prcio_p, praio_p, mnuccro_p, pracso_p, &
          preo_p, prdso_p) bind(c, name="micro_mg1_0_process_output_accum_codon")
       import c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c
       type(c_ptr), value :: qrtend_p, qnitend_p, qrtend_copy_p, qnitend_copy_p, cmei_p, cmeiout_p
       type(c_ptr), value :: prds_p, pre_p, cldmax_p, evapsnow_p, nevapr_p, nevapr2_p
       type(c_ptr), value :: pra_p, prc_p, lcldm_p, pracs_p, mnuccr_p, prain_p
       type(c_ptr), value :: prai_p, prci_p, icldm_p, psacws_p, bergs_p, prodsnow_p
       type(c_ptr), value :: qcsinksum_rate1ord_p, qcsum_rate1ord_p, qc_p, prao_p, prco_p
       type(c_ptr), value :: mnuccc_p, mnucct_p, mnuccd_p, msacwi_p, mnuccco_p, mnuccto_p
       type(c_ptr), value :: mnuccdo_p, msacwio_p, psacwso_p, bergso_p, berg_p, bergo_p
       type(c_ptr), value :: prcio_p, praio_p, mnuccro_p, pracso_p, preo_p, prdso_p
     end subroutine micro_mg1_0_process_output_accum_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_process_output_log_entry()
  call micro_mg1_0_process_output_accum_codon(int(i_local, c_int64_t), int(k_local, c_int64_t), &
       int(pcols_local, c_int64_t), int(pver_local, c_int64_t), c_loc(qrtend_local), &
       c_loc(qnitend_local), c_loc(qrtend_copy_local), c_loc(qnitend_copy_local), c_loc(cmei_local), &
       c_loc(cmeiout_local), c_loc(prds_local), c_loc(pre_local), c_loc(cldmax_local), &
       c_loc(evapsnow_local), c_loc(nevapr_local), c_loc(nevapr2_local), c_loc(pra_local), &
       c_loc(prc_local), c_loc(lcldm_local), c_loc(pracs_local), c_loc(mnuccr_local), &
       c_loc(prain_local), c_loc(prai_local), c_loc(prci_local), c_loc(icldm_local), &
       c_loc(psacws_local), c_loc(bergs_local), c_loc(prodsnow_local), c_loc(qcsinksum_rate1ord_local), &
       c_loc(qcsum_rate1ord_local), c_loc(qc_local), c_loc(prao_local), c_loc(prco_local), &
       c_loc(mnuccc_local), c_loc(mnucct_local), c_loc(mnuccd_local), c_loc(msacwi_local), &
       c_loc(mnuccco_local), c_loc(mnuccto_local), c_loc(mnuccdo_local), c_loc(msacwio_local), &
       c_loc(psacwso_local), c_loc(bergso_local), c_loc(berg_local), c_loc(bergo_local), &
       c_loc(prcio_local), c_loc(praio_local), c_loc(mnuccro_local), c_loc(pracso_local), &
       c_loc(preo_local), c_loc(prdso_local))
end subroutine micro_mg1_0_process_output_accum_codon_wrap

subroutine micro_mg1_0_post_iter_avg_codon_wrap(i_local, pcols_local, pver_local, top_lev_local, iter_local, &
     prect1_local, preci1_local, prect_local, preci_local, t1_local, q1_local, qc1_local, qi1_local, &
     nc1_local, ni1_local, t_local, q_local, qc_local, qi_local, nc_local, ni_local, &
     tlat1_local, qvlat1_local, qctend1_local, qitend1_local, nctend1_local, nitend1_local, &
     tlat_local, qvlat_local, qctend_local, qitend_local, nctend_local, nitend_local, rainrt1_local, &
     rainrt_local, rflx1_local, sflx1_local, rflx_local, sflx_local, qrout_local, qsout_local, &
     nrout_local, nsout_local, nevapr_local, nevapr2_local, evapsnow_local, prain_local, &
     prodsnow_local, cmeout_local, cmeiout_local, meltsdt_local, frzrdt_local, prao_local, prco_local, &
     mnuccco_local, mnuccto_local, msacwio_local, psacwso_local, bergso_local, bergo_local, &
     prcio_local, praio_local, mnuccro_local, pracso_local, mnuccdo_local, preo_local, prdso_local, &
     frzro_local, meltso_local, wtprelat_local, prer_evap_local)
  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, pcols_local, pver_local, top_lev_local, iter_local
  real(r8), target, intent(in) :: prect1_local(pcols_local), preci1_local(pcols_local)
  real(r8), target, intent(inout) :: prect_local(pcols_local), preci_local(pcols_local)
  real(r8), target, intent(in) :: t1_local(pcols_local,pver_local), q1_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qc1_local(pcols_local,pver_local), qi1_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nc1_local(pcols_local,pver_local), ni1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: t_local(pcols_local,pver_local), q_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qc_local(pcols_local,pver_local), qi_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nc_local(pcols_local,pver_local), ni_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: tlat1_local(pcols_local,pver_local), qvlat1_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qctend1_local(pcols_local,pver_local), qitend1_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nctend1_local(pcols_local,pver_local), nitend1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: tlat_local(pcols_local,pver_local), qvlat_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qctend_local(pcols_local,pver_local), qitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nctend_local(pcols_local,pver_local), nitend_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: rainrt1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: rainrt_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: rflx1_local(pcols_local,pver_local+1), sflx1_local(pcols_local,pver_local+1)
  real(r8), target, intent(inout) :: rflx_local(pcols_local,pver_local+1), sflx_local(pcols_local,pver_local+1)
  real(r8), target, intent(inout) :: qrout_local(pcols_local,pver_local), qsout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nrout_local(pcols_local,pver_local), nsout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nevapr_local(pcols_local,pver_local), nevapr2_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: evapsnow_local(pcols_local,pver_local), prain_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prodsnow_local(pcols_local,pver_local), cmeout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: cmeiout_local(pcols_local,pver_local), meltsdt_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: frzrdt_local(pcols_local,pver_local), prao_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: prco_local(pcols_local,pver_local), mnuccco_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: mnuccto_local(pcols_local,pver_local), msacwio_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: psacwso_local(pcols_local,pver_local), bergso_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: bergo_local(pcols_local,pver_local), prcio_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: praio_local(pcols_local,pver_local), mnuccro_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: pracso_local(pcols_local,pver_local), mnuccdo_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: preo_local(pcols_local,pver_local), prdso_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: frzro_local(pcols_local,pver_local), meltso_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: wtprelat_local(pcols_local,pver_local), prer_evap_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_post_iter_avg_codon(i_c, pcols_c, pver_c, top_lev_c, iter_c, &
          prect1_p, preci1_p, prect_p, preci_p, t1_p, q1_p, qc1_p, qi1_p, nc1_p, ni1_p, &
          t_p, q_p, qc_p, qi_p, nc_p, ni_p, tlat1_p, qvlat1_p, qctend1_p, qitend1_p, &
          nctend1_p, nitend1_p, tlat_p, qvlat_p, qctend_p, qitend_p, nctend_p, nitend_p, &
          rainrt1_p, rainrt_p, rflx1_p, sflx1_p, rflx_p, sflx_p, qrout_p, qsout_p, nrout_p, &
          nsout_p, nevapr_p, nevapr2_p, evapsnow_p, prain_p, prodsnow_p, cmeout_p, cmeiout_p, &
          meltsdt_p, frzrdt_p, prao_p, prco_p, mnuccco_p, mnuccto_p, msacwio_p, psacwso_p, &
          bergso_p, bergo_p, prcio_p, praio_p, mnuccro_p, pracso_p, mnuccdo_p, preo_p, prdso_p, &
          frzro_p, meltso_p, wtprelat_p, prer_evap_p) bind(c, name="micro_mg1_0_post_iter_avg_codon")
       import c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c, iter_c
       type(c_ptr), value :: prect1_p, preci1_p, prect_p, preci_p, t1_p, q1_p, qc1_p, qi1_p
       type(c_ptr), value :: nc1_p, ni1_p, t_p, q_p, qc_p, qi_p, nc_p, ni_p
       type(c_ptr), value :: tlat1_p, qvlat1_p, qctend1_p, qitend1_p, nctend1_p, nitend1_p
       type(c_ptr), value :: tlat_p, qvlat_p, qctend_p, qitend_p, nctend_p, nitend_p
       type(c_ptr), value :: rainrt1_p, rainrt_p, rflx1_p, sflx1_p, rflx_p, sflx_p
       type(c_ptr), value :: qrout_p, qsout_p, nrout_p, nsout_p, nevapr_p, nevapr2_p
       type(c_ptr), value :: evapsnow_p, prain_p, prodsnow_p, cmeout_p, cmeiout_p, meltsdt_p, frzrdt_p
       type(c_ptr), value :: prao_p, prco_p, mnuccco_p, mnuccto_p, msacwio_p, psacwso_p
       type(c_ptr), value :: bergso_p, bergo_p, prcio_p, praio_p, mnuccro_p, pracso_p, mnuccdo_p
       type(c_ptr), value :: preo_p, prdso_p, frzro_p, meltso_p, wtprelat_p, prer_evap_p
     end subroutine micro_mg1_0_post_iter_avg_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_post_iter_avg_log_entry()
  call micro_mg1_0_post_iter_avg_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), int(iter_local, c_int64_t), &
       c_loc(prect1_local), c_loc(preci1_local), c_loc(prect_local), c_loc(preci_local), &
       c_loc(t1_local), c_loc(q1_local), c_loc(qc1_local), c_loc(qi1_local), c_loc(nc1_local), &
       c_loc(ni1_local), c_loc(t_local), c_loc(q_local), c_loc(qc_local), c_loc(qi_local), &
       c_loc(nc_local), c_loc(ni_local), c_loc(tlat1_local), c_loc(qvlat1_local), &
       c_loc(qctend1_local), c_loc(qitend1_local), c_loc(nctend1_local), c_loc(nitend1_local), &
       c_loc(tlat_local), c_loc(qvlat_local), c_loc(qctend_local), c_loc(qitend_local), &
       c_loc(nctend_local), c_loc(nitend_local), c_loc(rainrt1_local), c_loc(rainrt_local), &
       c_loc(rflx1_local), c_loc(sflx1_local), c_loc(rflx_local), c_loc(sflx_local), &
       c_loc(qrout_local), c_loc(qsout_local), c_loc(nrout_local), c_loc(nsout_local), &
       c_loc(nevapr_local), c_loc(nevapr2_local), c_loc(evapsnow_local), c_loc(prain_local), &
       c_loc(prodsnow_local), c_loc(cmeout_local), c_loc(cmeiout_local), c_loc(meltsdt_local), &
       c_loc(frzrdt_local), c_loc(prao_local), c_loc(prco_local), c_loc(mnuccco_local), &
       c_loc(mnuccto_local), c_loc(msacwio_local), c_loc(psacwso_local), c_loc(bergso_local), &
       c_loc(bergo_local), c_loc(prcio_local), c_loc(praio_local), c_loc(mnuccro_local), &
       c_loc(pracso_local), c_loc(mnuccdo_local), c_loc(preo_local), c_loc(prdso_local), &
       c_loc(frzro_local), c_loc(meltso_local), c_loc(wtprelat_local), c_loc(prer_evap_local))
end subroutine micro_mg1_0_post_iter_avg_codon_wrap

subroutine micro_mg1_0_phase_change_codon_wrap(i_local, k_local, pcols_local, pver_local, &
     deltat_local, cpp_local, xlf_local, tmelt_local, qsmall_local, pi_local, rhow_local, do_cldice_local, &
     qc_local, qi_local, nc_local, ni_local, t_local, qctend_local, qitend_local, nctend_local, nitend_local, &
     tlat_local, dumc_local, dumi_local, dumnc_local, dumni_local, melto_local, homoo_local, wtpostlat_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, k_local, pcols_local, pver_local
  real(r8), intent(in) :: deltat_local, cpp_local, xlf_local, tmelt_local, qsmall_local, pi_local, rhow_local
  logical, intent(in) :: do_cldice_local
  real(r8), target, intent(in) :: qc_local(pcols_local,pver_local), qi_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nc_local(pcols_local,pver_local), ni_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: t_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qctend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nctend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: tlat_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: dumc_local(pcols_local,pver_local), dumi_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: dumnc_local(pcols_local,pver_local), dumni_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: melto_local(pcols_local,pver_local), homoo_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: wtpostlat_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_phase_change_codon(i_c, k_c, pcols_c, pver_c, deltat_c, cpp_c, xlf_c, &
          tmelt_c, qsmall_c, pi_c, rhow_c, do_cldice_c, qc_p, qi_p, nc_p, ni_p, t_p, qctend_p, &
          qitend_p, nctend_p, nitend_p, tlat_p, dumc_p, dumi_p, dumnc_p, dumni_p, melto_p, &
          homoo_p, wtpostlat_p) bind(c, name="micro_mg1_0_phase_change_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, k_c, pcols_c, pver_c, do_cldice_c
       real(c_double), value :: deltat_c, cpp_c, xlf_c, tmelt_c, qsmall_c, pi_c, rhow_c
       type(c_ptr), value :: qc_p, qi_p, nc_p, ni_p, t_p, qctend_p, qitend_p, nctend_p, nitend_p
       type(c_ptr), value :: tlat_p, dumc_p, dumi_p, dumnc_p, dumni_p, melto_p, homoo_p, wtpostlat_p
     end subroutine micro_mg1_0_phase_change_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_phase_change_log_entry()
  call micro_mg1_0_phase_change_codon(int(i_local, c_int64_t), int(k_local, c_int64_t), &
       int(pcols_local, c_int64_t), int(pver_local, c_int64_t), real(deltat_local, c_double), &
       real(cpp_local, c_double), real(xlf_local, c_double), real(tmelt_local, c_double), &
       real(qsmall_local, c_double), real(pi_local, c_double), real(rhow_local, c_double), &
       merge(1_c_int64_t, 0_c_int64_t, do_cldice_local), c_loc(qc_local), c_loc(qi_local), &
       c_loc(nc_local), c_loc(ni_local), c_loc(t_local), c_loc(qctend_local), c_loc(qitend_local), &
       c_loc(nctend_local), c_loc(nitend_local), c_loc(tlat_local), c_loc(dumc_local), &
       c_loc(dumi_local), c_loc(dumnc_local), c_loc(dumni_local), c_loc(melto_local), &
       c_loc(homoo_local), c_loc(wtpostlat_local))
end subroutine micro_mg1_0_phase_change_codon_wrap

subroutine micro_mg1_0_number_cleanup_codon_wrap(ncol_local, pcols_local, pver_local, top_lev_local, &
     deltat_local, qsmall_local, do_cldice_local, qc_local, qi_local, nc_local, ni_local, &
     qctend_local, qitend_local, nctend_local, nitend_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: ncol_local, pcols_local, pver_local, top_lev_local
  real(r8), intent(in) :: deltat_local, qsmall_local
  logical, intent(in) :: do_cldice_local
  real(r8), target, intent(in) :: qc_local(pcols_local,pver_local), qi_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nc_local(pcols_local,pver_local), ni_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qctend_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qitend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nctend_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nitend_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_number_cleanup_codon(ncol_c, pcols_c, pver_c, top_lev_c, &
          deltat_c, qsmall_c, do_cldice_c, qc_p, qi_p, nc_p, ni_p, qctend_p, qitend_p, &
          nctend_p, nitend_p) bind(c, name="micro_mg1_0_number_cleanup_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, do_cldice_c
       real(c_double), value :: deltat_c, qsmall_c
       type(c_ptr), value :: qc_p, qi_p, nc_p, ni_p, qctend_p, qitend_p, nctend_p, nitend_p
     end subroutine micro_mg1_0_number_cleanup_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_number_cleanup_log_entry()
  call micro_mg1_0_number_cleanup_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), real(deltat_local, c_double), &
       real(qsmall_local, c_double), merge(1_c_int64_t, 0_c_int64_t, do_cldice_local), &
       c_loc(qc_local), c_loc(qi_local), c_loc(nc_local), c_loc(ni_local), c_loc(qctend_local), &
       c_loc(qitend_local), c_loc(nctend_local), c_loc(nitend_local))
end subroutine micro_mg1_0_number_cleanup_codon_wrap

subroutine micro_mg1_0_reflectivity_flags_codon_wrap(ncol_local, pcols_local, pver_local, top_lev_local, &
     mindbz_local, csmin_local, csmax_local, refl_local, arefl_local, areflz_local, frefl_local, &
     csrfl_local, acsrfl_local, fcsrfl_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: ncol_local, pcols_local, pver_local, top_lev_local
  real(r8), intent(in) :: mindbz_local, csmin_local, csmax_local
  real(r8), target, intent(inout) :: refl_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: arefl_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: areflz_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: frefl_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: csrfl_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: acsrfl_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: fcsrfl_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_reflectivity_flags_codon(ncol_c, pcols_c, pver_c, top_lev_c, &
          mindbz_c, csmin_c, csmax_c, refl_p, arefl_p, areflz_p, frefl_p, csrfl_p, acsrfl_p, &
          fcsrfl_p) bind(c, name="micro_mg1_0_reflectivity_flags_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       real(c_double), value :: mindbz_c, csmin_c, csmax_c
       type(c_ptr), value :: refl_p, arefl_p, areflz_p, frefl_p, csrfl_p, acsrfl_p, fcsrfl_p
     end subroutine micro_mg1_0_reflectivity_flags_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_reflectivity_flags_log_entry()
  call micro_mg1_0_reflectivity_flags_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), real(mindbz_local, c_double), &
       real(csmin_local, c_double), real(csmax_local, c_double), c_loc(refl_local), c_loc(arefl_local), &
       c_loc(areflz_local), c_loc(frefl_local), c_loc(csrfl_local), c_loc(acsrfl_local), &
       c_loc(fcsrfl_local))
end subroutine micro_mg1_0_reflectivity_flags_codon_wrap

subroutine micro_mg1_0_substep_accum_column_codon_wrap(i_local, pcols_local, pver_local, &
     top_lev_local, deltat_local, cpp_local, qric_local, qniic_local, nric_local, nsic_local, rho_local, &
     cldmax_local, qrout_local, qsout_local, nrout_local, nsout_local, tlat_local, qvlat_local, &
     qctend_local, qitend_local, nctend_local, nitend_local, tlat1_local, qvlat1_local, &
     qctend1_local, qitend1_local, nctend1_local, nitend1_local, t_local, q_local, qc_local, qi_local, &
     nc_local, ni_local, rainrt_local, rainrt1_local, arcld_local, rercld_local, rflx_local, sflx_local, &
     rflx1_local, sflx1_local, umr_local, ums_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: i_local, pcols_local, pver_local, top_lev_local
  real(r8), intent(in) :: deltat_local, cpp_local
  real(r8), target, intent(in) :: qric_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qniic_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nric_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nsic_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: rho_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: cldmax_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: tlat_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qvlat_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qctend_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qitend_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nctend_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nitend_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: rainrt_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: arcld_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: umr_local(pver_local), ums_local(pver_local)
  real(r8), target, intent(inout) :: qrout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qsout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nrout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nsout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: tlat1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qvlat1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qctend1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qitend1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nctend1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nitend1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: t_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: q_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qc_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qi_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nc_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: ni_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: rainrt1_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: rercld_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: rflx_local(pcols_local,pver_local+1)
  real(r8), target, intent(inout) :: sflx_local(pcols_local,pver_local+1)
  real(r8), target, intent(inout) :: rflx1_local(pcols_local,pver_local+1)
  real(r8), target, intent(inout) :: sflx1_local(pcols_local,pver_local+1)

  interface
     subroutine micro_mg1_0_substep_accum_column_codon(i_c, pcols_c, pver_c, top_lev_c, &
          deltat_c, cpp_c, qric_p, qniic_p, nric_p, nsic_p, rho_p, cldmax_p, qrout_p, qsout_p, &
          nrout_p, nsout_p, tlat_p, qvlat_p, qctend_p, qitend_p, nctend_p, nitend_p, tlat1_p, &
          qvlat1_p, qctend1_p, qitend1_p, nctend1_p, nitend1_p, t_p, q_p, qc_p, qi_p, nc_p, &
          ni_p, rainrt_p, rainrt1_p, arcld_p, rercld_p, rflx_p, sflx_p, rflx1_p, sflx1_p, &
          umr_p, ums_p) bind(c, name="micro_mg1_0_substep_accum_column_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: i_c, pcols_c, pver_c, top_lev_c
       real(c_double), value :: deltat_c, cpp_c
       type(c_ptr), value :: qric_p, qniic_p, nric_p, nsic_p, rho_p, cldmax_p, qrout_p, qsout_p
       type(c_ptr), value :: nrout_p, nsout_p, tlat_p, qvlat_p, qctend_p, qitend_p, nctend_p, nitend_p
       type(c_ptr), value :: tlat1_p, qvlat1_p, qctend1_p, qitend1_p, nctend1_p, nitend1_p
       type(c_ptr), value :: t_p, q_p, qc_p, qi_p, nc_p, ni_p, rainrt_p, rainrt1_p, arcld_p, rercld_p
       type(c_ptr), value :: rflx_p, sflx_p, rflx1_p, sflx1_p, umr_p, ums_p
     end subroutine micro_mg1_0_substep_accum_column_codon
  end interface

  call micro_mg1_0_colzero_log_entry()
  call micro_mg1_0_substep_accum_log_entry()
  call micro_mg1_0_substep_accum_column_codon(int(i_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), &
       real(deltat_local, c_double), real(cpp_local, c_double), c_loc(qric_local), c_loc(qniic_local), &
       c_loc(nric_local), c_loc(nsic_local), c_loc(rho_local), c_loc(cldmax_local), c_loc(qrout_local), &
       c_loc(qsout_local), c_loc(nrout_local), c_loc(nsout_local), c_loc(tlat_local), c_loc(qvlat_local), &
       c_loc(qctend_local), c_loc(qitend_local), c_loc(nctend_local), c_loc(nitend_local), &
       c_loc(tlat1_local), c_loc(qvlat1_local), c_loc(qctend1_local), c_loc(qitend1_local), &
       c_loc(nctend1_local), c_loc(nitend1_local), c_loc(t_local), c_loc(q_local), c_loc(qc_local), &
       c_loc(qi_local), c_loc(nc_local), c_loc(ni_local), c_loc(rainrt_local), c_loc(rainrt1_local), &
       c_loc(arcld_local), c_loc(rercld_local), c_loc(rflx_local), c_loc(sflx_local), c_loc(rflx1_local), &
       c_loc(sflx1_local), c_loc(umr_local), c_loc(ums_local))
end subroutine micro_mg1_0_substep_accum_column_codon_wrap

subroutine micro_mg1_0_tail_activation_codon_wrap(ncol_local, pcols_local, pver_local, top_lev_local, &
     dum2i_local, dum2l_local, rho_local, ncai_local, ncal_local)
  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  integer, intent(in) :: ncol_local, pcols_local, pver_local, top_lev_local
  real(r8), target, intent(in) :: dum2i_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: dum2l_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: rho_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: ncai_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: ncal_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_tail_activation_codon(ncol_c, pcols_c, pver_c, top_lev_c, &
          dum2i_p, dum2l_p, rho_p, ncai_p, ncal_p) bind(c, name="micro_mg1_0_tail_activation_codon")
       import c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: dum2i_p, dum2l_p, rho_p, ncai_p, ncal_p
     end subroutine micro_mg1_0_tail_activation_codon
  end interface

  call micro_mg1_0_tail_diag_log_entry()
  call micro_mg1_0_tail_activation_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), c_loc(dum2i_local), &
       c_loc(dum2l_local), c_loc(rho_local), c_loc(ncai_local), c_loc(ncal_local))
end subroutine micro_mg1_0_tail_activation_codon_wrap

subroutine micro_mg1_0_tail_avg_codon_wrap(ncol_local, pcols_local, pver_local, top_lev_local, &
     qrout_local, qsout_local, nrout_local, nsout_local, qrout2_local, qsout2_local, nrout2_local, &
     nsout2_local, drout2_local, dsout2_local, freqs_local, freqr_local)
  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  integer, intent(in) :: ncol_local, pcols_local, pver_local, top_lev_local
  real(r8), target, intent(in) :: qrout_local(pcols_local,pver_local), qsout_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: nrout_local(pcols_local,pver_local), nsout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: qrout2_local(pcols_local,pver_local), qsout2_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nrout2_local(pcols_local,pver_local), nsout2_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: drout2_local(pcols_local,pver_local), dsout2_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: freqs_local(pcols_local,pver_local), freqr_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_tail_avg_codon(ncol_c, pcols_c, pver_c, top_lev_c, &
          qrout_p, qsout_p, nrout_p, nsout_p, qrout2_p, qsout2_p, nrout2_p, nsout2_p, &
          drout2_p, dsout2_p, freqs_p, freqr_p) bind(c, name="micro_mg1_0_tail_avg_codon")
       import c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       type(c_ptr), value :: qrout_p, qsout_p, nrout_p, nsout_p, qrout2_p, qsout2_p
       type(c_ptr), value :: nrout2_p, nsout2_p, drout2_p, dsout2_p, freqs_p, freqr_p
     end subroutine micro_mg1_0_tail_avg_codon
  end interface

  call micro_mg1_0_tail_diag_log_entry()
  call micro_mg1_0_tail_avg_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), c_loc(qrout_local), &
       c_loc(qsout_local), c_loc(nrout_local), c_loc(nsout_local), c_loc(qrout2_local), &
       c_loc(qsout2_local), c_loc(nrout2_local), c_loc(nsout2_local), c_loc(drout2_local), &
       c_loc(dsout2_local), c_loc(freqs_local), c_loc(freqr_local))
end subroutine micro_mg1_0_tail_avg_codon_wrap

subroutine micro_mg1_0_tail_fice_codon_wrap(ncol_local, pcols_local, pver_local, top_lev_local, &
     deltat_local, qsmall_local, qc_local, qi_local, qctend_local, qitend_local, qsout_local, qrout_local, &
     dumc_local, dumi_local, nfice_local)
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  integer, intent(in) :: ncol_local, pcols_local, pver_local, top_lev_local
  real(r8), intent(in) :: deltat_local, qsmall_local
  real(r8), target, intent(in) :: qc_local(pcols_local,pver_local), qi_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qctend_local(pcols_local,pver_local), qitend_local(pcols_local,pver_local)
  real(r8), target, intent(in) :: qsout_local(pcols_local,pver_local), qrout_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: dumc_local(pcols_local,pver_local), dumi_local(pcols_local,pver_local)
  real(r8), target, intent(inout) :: nfice_local(pcols_local,pver_local)

  interface
     subroutine micro_mg1_0_tail_fice_codon(ncol_c, pcols_c, pver_c, top_lev_c, deltat_c, qsmall_c, &
          qc_p, qi_p, qctend_p, qitend_p, qsout_p, qrout_p, dumc_p, dumi_p, nfice_p) &
          bind(c, name="micro_mg1_0_tail_fice_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
       real(c_double), value :: deltat_c, qsmall_c
       type(c_ptr), value :: qc_p, qi_p, qctend_p, qitend_p, qsout_p, qrout_p, dumc_p, dumi_p, nfice_p
     end subroutine micro_mg1_0_tail_fice_codon
  end interface

  call micro_mg1_0_tail_diag_log_entry()
  call micro_mg1_0_tail_fice_codon(int(ncol_local, c_int64_t), int(pcols_local, c_int64_t), &
       int(pver_local, c_int64_t), int(top_lev_local, c_int64_t), real(deltat_local, c_double), &
       real(qsmall_local, c_double), c_loc(qc_local), c_loc(qi_local), c_loc(qctend_local), &
       c_loc(qitend_local), c_loc(qsout_local), c_loc(qrout_local), c_loc(dumc_local), &
       c_loc(dumi_local), c_loc(nfice_local))
end subroutine micro_mg1_0_tail_fice_codon_wrap

!========================================================================
!UTILITIES
!========================================================================

subroutine micro_mg1_0_select_get_cols_impl()
  character(len=32) :: impl_name
  integer :: n, status

  if (micro_mg1_0_get_cols_impl_selected) return

  call get_environment_variable('MICRO_MG1_0_GET_COLS_IMPL', value=impl_name, length=n, status=status)
  if (status == 0 .and. n > 0) then
     micro_mg1_0_get_cols_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     micro_mg1_0_get_cols_use_native_impl = .false.
  end if

  if (masterproc) then
     if (micro_mg1_0_get_cols_use_native_impl) then
        write(iulog,*) 'micro_mg1_0_get_cols implementation = native'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_GET_COLS_PROOF_FILE', &
             'micro_mg1_0_get_cols implementation = native')
     else
        write(iulog,*) 'micro_mg1_0_get_cols implementation = codon'
        call micro_mg1_0_append_impl_proof('MICRO_MG1_0_GET_COLS_PROOF_FILE', &
             'micro_mg1_0_get_cols implementation = codon')
     end if
  end if

  micro_mg1_0_get_cols_impl_selected = .true.
end subroutine micro_mg1_0_select_get_cols_impl

subroutine micro_mg1_0_get_cols_log_entry()
  if (masterproc .and. .not. micro_mg1_0_get_cols_logged) then
     write(iulog,*) 'micro_mg_get_cols direct = codon'
     call micro_mg1_0_append_impl_proof('MICRO_MG1_0_GET_COLS_PROOF_FILE', &
          'micro_mg_get_cols direct = codon')
     micro_mg1_0_get_cols_logged = .true.
  end if
end subroutine micro_mg1_0_get_cols_log_entry

subroutine micro_mg_get_cols(ncol, nlev, top_lev, qcn, qin, &
     mgncol, mgcols)

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr, c_null_ptr

  ! Determines which columns microphysics should operate over by
  ! checking for non-zero cloud water/ice.

  integer, intent(in) :: ncol      ! Number of columns with meaningful data
  integer, intent(in) :: nlev      ! Number of levels to use
  integer, intent(in) :: top_lev   ! Top level for microphysics

  real(r8), target, intent(in) :: qcn(:,:) ! cloud water mixing ratio (kg/kg)
  real(r8), target, intent(in) :: qin(:,:) ! cloud ice mixing ratio (kg/kg)

  integer, intent(out) :: mgncol   ! Number of columns MG will use
  integer, allocatable, target, intent(out) :: mgcols(:) ! column indices

  integer :: lev_offset  ! top_lev - 1 (defined here for consistency)
  logical :: ltrue(ncol) ! store tests for each column

  integer :: i, ii ! column indices
  integer(c_int64_t), target :: mgncol_c

  interface
     subroutine micro_mg_get_cols_codon(stage_c, ncol_c, ldq_c, nlev_c, top_lev_c, qsmall_c, &
          qcn_p, qin_p, mgncol_p, mgcols_p) bind(c, name="micro_mg_get_cols_codon")
       import c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: stage_c, ncol_c, ldq_c, nlev_c, top_lev_c
       real(c_double), value :: qsmall_c
       type(c_ptr), value :: qcn_p, qin_p, mgncol_p, mgcols_p
     end subroutine micro_mg_get_cols_codon
  end interface

  if (allocated(mgcols)) deallocate(mgcols)

  call micro_mg1_0_select_get_cols_impl()

  if (.not. micro_mg1_0_get_cols_use_native_impl) then
     mgncol_c = 0_c_int64_t
     call micro_mg1_0_get_cols_log_entry()
     call micro_mg_get_cols_codon(0_c_int64_t, int(ncol, c_int64_t), int(size(qcn, 1), c_int64_t), &
          int(nlev, c_int64_t), int(top_lev, c_int64_t), real(qsmall, c_double), &
          c_loc(qcn(1,1)), c_loc(qin(1,1)), c_loc(mgncol_c), c_null_ptr)
     mgncol = int(mgncol_c)
     allocate(mgcols(mgncol))
     if (mgncol > 0) then
        call micro_mg_get_cols_codon(1_c_int64_t, int(ncol, c_int64_t), int(size(qcn, 1), c_int64_t), &
             int(nlev, c_int64_t), int(top_lev, c_int64_t), real(qsmall, c_double), &
             c_loc(qcn(1,1)), c_loc(qin(1,1)), c_loc(mgncol_c), c_loc(mgcols(1)))
     end if
     return
  end if

  lev_offset = top_lev - 1

  ! Using "any" along dimension 2 collapses across levels, but
  ! not columns, so we know if water is present at any level
  ! in each column.

  ltrue = any(qcn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qin(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)

  ! Scan for true values to get a usable list of indices.

  mgncol = count(ltrue)
  allocate(mgcols(mgncol))
  i = 0
  do ii = 1,ncol
     if (ltrue(ii)) then
        i = i + 1
        mgcols(i) = ii
     end if
  end do

end subroutine micro_mg_get_cols

end module micro_mg1_0
