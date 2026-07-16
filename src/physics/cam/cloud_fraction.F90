module cloud_fraction

  ! Cloud fraction parameterization.


  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver, pverp
  use ref_pres,       only: pref_mid
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use ref_pres,       only: trop_cloud_top_lev
  use ap_cldfrc_scheme, only: cldfrc_run

  implicit none
  private
  save

  ! Public interfaces
  public &
     cldfrc_readnl,    &! read cldfrc_nl namelist
     cldfrc_register,  &! add fields to pbuf
     cldfrc_init,      &! Inititialization of cloud_fraction run-time parameters
     cldfrc_getparams, &! public access of tuning parameters
     cldfrc,           &! Computation of cloud fraction
     cldfrc_fice        ! Calculate fraction of condensate in ice phase (radiation partitioning)

  ! Private data
  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  ! Top level
  integer :: top_lev = 1

  ! Physics buffer indices
  integer :: sh_frac_idx   = 0
  integer :: dp_frac_idx   = 0

  ! Namelist variables
  logical  :: cldfrc_freeze_dry           ! switch for Vavrus correction
  logical  :: cldfrc_ice                  ! switch to compute ice cloud fraction
  real(r8) :: cldfrc_rhminl = unset_r8    ! minimum rh for low stable clouds
  real(r8) :: cldfrc_rhminl_adj_land = unset_r8   ! rhminl adjustment for snowfree land
  real(r8) :: cldfrc_rhminh = unset_r8    ! minimum rh for high stable clouds
  real(r8) :: cldfrc_rhminp = unset_r8    ! minimum rh for high stable clouds poleward of 60 degrees
  real(r8) :: cldfrc_rhminp_botmb = 300._r8 ! and pressures less than cldfrc_rhminp_botmb (hPa)
  real(r8) :: cldfrc_sh1    = unset_r8    ! parameter for shallow convection cloud fraction
  real(r8) :: cldfrc_sh2    = unset_r8    ! parameter for shallow convection cloud fraction
  real(r8) :: cldfrc_dp1    = unset_r8    ! parameter for deep convection cloud fraction
  real(r8) :: cldfrc_dp2    = unset_r8    ! parameter for deep convection cloud fraction
  real(r8) :: cldfrc_premit = unset_r8    ! top pressure bound for mid level cloud
  real(r8) :: cldfrc_premib  = unset_r8   ! bottom pressure bound for mid level cloud
  integer  :: cldfrc_iceopt               ! option for ice cloud closure
                                          ! 1=wang & sassen 2=schiller (iciwc)
                                          ! 3=wood & field, 4=Wilson (based on smith)
  real(r8) :: cldfrc_icecrit = unset_r8   ! Critical RH for ice clouds in Wilson & Ballard closure (smaller = more ice clouds)

  real(r8) :: rhminl             ! set from namelist input cldfrc_rhminl
  real(r8) :: rhminl_adj_land    ! set from namelist input cldfrc_rhminl_adj_land
  real(r8) :: rhminh             ! set from namelist input cldfrc_rhminh
  real(r8) :: rhminp             ! set from namelist input cldfrc_rhminp
  real(r8) :: sh1, sh2           ! set from namelist input cldfrc_sh1, cldfrc_sh2
  real(r8) :: dp1,dp2            ! set from namelist input cldfrc_dp1, cldfrc_dp2
  real(r8) :: premit             ! set from namelist input cldfrc_premit
  real(r8) :: premib             ! set from namelist input cldfrc_premib
  integer  :: iceopt             ! set from namelist input cldfrc_iceopt
  real(r8) :: icecrit            ! set from namelist input cldfrc_icecrit

  ! constants
  real(r8), parameter :: pnot = 1.e5_r8         ! reference pressure
  real(r8), parameter :: lapse = 6.5e-3_r8      ! U.S. Standard Atmosphere lapse rate
  real(r8), parameter :: pretop = 1.0e2_r8      ! pressure bounding high cloud

  integer count

  logical :: inversion_cld_off    ! Turns off stratification-based cld frc

  integer :: k700   ! model level nearest 700 mb

!================================================================================================
  contains
!================================================================================================

subroutine cldfrc_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cldfrc_readnl'

   namelist /cldfrc_nl/ cldfrc_freeze_dry,      cldfrc_ice,    cldfrc_rhminl, &
                        cldfrc_rhminl_adj_land, cldfrc_rhminh, cldfrc_sh1,    &
                        cldfrc_rhminp,          cldfrc_rhminp_botmb, &
                        cldfrc_sh2,             cldfrc_dp1,    cldfrc_dp2,    &
                        cldfrc_premit,          cldfrc_premib, cldfrc_iceopt, &
                        cldfrc_icecrit
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cldfrc_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cldfrc_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      rhminl = cldfrc_rhminl
      rhminl_adj_land = cldfrc_rhminl_adj_land
      rhminh = cldfrc_rhminh
      rhminp = cldfrc_rhminp
      sh1    = cldfrc_sh1
      sh2    = cldfrc_sh2
      dp1    = cldfrc_dp1
      dp2    = cldfrc_dp2
      premit = cldfrc_premit
      premib  = cldfrc_premib
      iceopt  = cldfrc_iceopt
      icecrit = cldfrc_icecrit

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(cldfrc_freeze_dry, 1, mpilog, 0, mpicom)
   call mpibcast(cldfrc_ice,        1, mpilog, 0, mpicom)
   call mpibcast(rhminl,            1, mpir8,  0, mpicom)
   call mpibcast(rhminl_adj_land,   1, mpir8,  0, mpicom)
   call mpibcast(rhminh,            1, mpir8,  0, mpicom)
   call mpibcast(rhminp,            1, mpir8,  0, mpicom)
   call mpibcast(sh1   ,            1, mpir8,  0, mpicom)
   call mpibcast(sh2   ,            1, mpir8,  0, mpicom)
   call mpibcast(dp1   ,            1, mpir8,  0, mpicom)
   call mpibcast(dp2   ,            1, mpir8,  0, mpicom)
   call mpibcast(premit,            1, mpir8,  0, mpicom)
   call mpibcast(premib,            1, mpir8,  0, mpicom)
   call mpibcast(iceopt,            1, mpiint, 0, mpicom)
   call mpibcast(icecrit,           1, mpir8,  0, mpicom)
#endif

end subroutine cldfrc_readnl

!================================================================================================

subroutine cldfrc_register

   ! Register fields in the physics buffer.

   use physics_buffer, only : pbuf_add_field, dtype_r8

   !-----------------------------------------------------------------------

   call pbuf_add_field('SH_FRAC', 'physpkg', dtype_r8, (/pcols,pver/), sh_frac_idx)
   call pbuf_add_field('DP_FRAC', 'physpkg', dtype_r8, (/pcols,pver/), dp_frac_idx)

end subroutine cldfrc_register

!================================================================================================

subroutine cldfrc_getparams(rhminl_out, rhminl_adj_land_out, rhminh_out,  premit_out, &
                            rhminp_out, premib_out, iceopt_out, icecrit_out)
!-----------------------------------------------------------------------
! Purpose: Return cldfrc tuning parameters
!-----------------------------------------------------------------------

   real(r8),          intent(out), optional :: rhminl_out
   real(r8),          intent(out), optional :: rhminl_adj_land_out
   real(r8),          intent(out), optional :: rhminh_out
   real(r8),          intent(out), optional :: rhminp_out
   real(r8),          intent(out), optional :: premit_out
   real(r8),          intent(out), optional :: premib_out
   integer,           intent(out), optional :: iceopt_out
   real(r8),          intent(out), optional :: icecrit_out

   if ( present(rhminl_out) )      rhminl_out = rhminl
   if ( present(rhminl_adj_land_out) ) rhminl_adj_land_out = rhminl_adj_land
   if ( present(rhminh_out) )      rhminh_out = rhminh
   if ( present(rhminp_out) )      rhminp_out = rhminp
   if ( present(premit_out) )      premit_out = premit
   if ( present(premib_out) )      premib_out  = premib
   if ( present(iceopt_out) )      iceopt_out  = iceopt
   if ( present(icecrit_out) )     icecrit_out = icecrit

end subroutine cldfrc_getparams

!===============================================================================

subroutine cldfrc_init

   ! Initialize cloud fraction run-time parameters

   use cam_history,   only:  phys_decomp, addfld
   use dycore,        only:  dycore_is, get_resolution
   use phys_control,  only:  phys_getopts

   ! horizontal grid specifier
   character(len=32) :: hgrid

   ! query interfaces for scheme settings
   character(len=16) :: shallow_scheme, eddy_scheme, macrop_scheme

   integer :: k
   !-----------------------------------------------------------------------------

   call phys_getopts(shallow_scheme_out = shallow_scheme ,&
                     eddy_scheme_out    = eddy_scheme    ,&
                     macrop_scheme_out  = macrop_scheme  )

   ! Limit CAM5 cloud physics to below top cloud level.
   if (macrop_scheme /= "rk") top_lev = trop_cloud_top_lev

   hgrid = get_resolution()

   ! Turn off inversion_cld if any UW PBL scheme is being used
   if ( (eddy_scheme .eq. 'diag_TKE' ) .or. (shallow_scheme .eq.  'UW' )) then
      inversion_cld_off = .true.
   else
      inversion_cld_off = .false.
   endif

   if ( masterproc ) then
      write(iulog,*)'tuning parameters cldfrc_init: inversion_cld_off',inversion_cld_off
      write(iulog,*)'tuning parameters cldfrc_init: dp1',dp1,'dp2',dp2,'sh1',sh1,'sh2',sh2
      if (shallow_scheme .ne. 'UW' ) then
         write(iulog,*)'tuning parameters cldfrc_init: rhminl',rhminl,'rhminl_adj_land',rhminl_adj_land, &
                       'rhminh',rhminh,'premit',premit,'premib',premib
         write(iulog,*)'tuning parameters cldfrc_init: iceopt',iceopt,'icecrit',icecrit
      endif
   endif

   if (pref_mid(top_lev) > 7.e4_r8) &
        call endrun ('cldfrc_init: model levels bracketing 700 mb not found')

   ! Find vertical level nearest 700 mb.
   k700 = minloc(abs(pref_mid(top_lev:pver) - 7.e4_r8), 1)

   if (masterproc) then
      write(iulog,*)'cldfrc_init: model level nearest 700 mb is',k700,'which is',pref_mid(k700),'pascals'
   end if

   call addfld ('SH_CLD   ', 'fraction', pver, 'A', 'Shallow convective cloud cover'                          ,phys_decomp)
   call addfld ('DP_CLD   ', 'fraction', pver, 'A', 'Deep convective cloud cover'                             ,phys_decomp)

end subroutine cldfrc_init

!===============================================================================

!> \section arg_table_cldfrc_run Argument Table
!! \htmlinclude cldfrc_run.html
subroutine cldfrc(lchnk, ncol, pbuf, &
       pmid, temp, q, omga, phis, shfrc, use_shfrc, &
       cloud, rhcloud, clc, pdel, cmfmc, cmfmc2, landfrac, snowh, &
       concld, cldst, ts, sst, ps, zdu, ocnfrac, rhu00, cldice, &
       icecldf, liqcldf, relhum, dindex)

    use cam_history, only: outfld
    use physconst, only: cappa, gravit, rair, pi
    use wv_saturation, only: qsat, qsat_water, svp_ice
    use phys_grid, only: get_rlat_all_p, get_rlon_all_p
    use perf_mod, only: t_startf, t_stopf
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field

    integer, intent(in) :: lchnk, ncol, dindex
    type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:)
    real(r8), intent(in) :: pmid(:,:), temp(:,:), q(:,:), omga(:,:)
    real(r8), intent(in) :: phis(:), shfrc(:,:), pdel(:,:)
    real(r8), intent(in) :: cmfmc(:,:), cmfmc2(:,:)
    real(r8), intent(in) :: landfrac(:), snowh(:), ts(:), sst(:), ps(:)
    real(r8), intent(in) :: zdu(:,:), ocnfrac(:), cldice(:,:)
    logical, intent(in) :: use_shfrc
    real(r8), intent(out) :: cloud(:,:), rhcloud(:,:), clc(:)
    real(r8), intent(out) :: concld(:,:), cldst(:,:), rhu00(:,:)
    real(r8), intent(out) :: icecldf(:,:), liqcldf(:,:), relhum(:,:)

    real(r8), pointer :: shallowcu(:,:), deepcu(:,:)
    real(r8) :: clat(pcols), clon(pcols)
    real(r8) :: es(pcols,pver), qs(pcols,pver)
    real(r8) :: esl(pcols,pver), esi(pcols,pver)
    integer :: cold_sst_index

    call t_startf('ap_cldfrc_run')

    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)
    call pbuf_get_field(pbuf, sh_frac_idx, shallowcu)
    call pbuf_get_field(pbuf, dp_frac_idx, deepcu)

    if (cldfrc_ice) then
       call qsat_water(temp(1:ncol,top_lev:pver), pmid(1:ncol,top_lev:pver), &
            esl(1:ncol,top_lev:pver), qs(1:ncol,top_lev:pver))
       esi(1:ncol,top_lev:pver) = svp_ice(temp(1:ncol,top_lev:pver))
    else
       call qsat(temp(1:ncol,top_lev:pver), pmid(1:ncol,top_lev:pver), &
            es(1:ncol,top_lev:pver), qs(1:ncol,top_lev:pver))
    end if

    call cldfrc_run( &
         pcols, pver, ncol, top_lev, k700, dindex, iceopt, &
         cldfrc_freeze_dry, cldfrc_ice, inversion_cld_off, use_shfrc, &
         rhminl, rhminl_adj_land, rhminh, rhminp, cldfrc_rhminp_botmb, &
         sh1, sh2, dp1, dp2, premit, premib, icecrit, unset_r8, &
         cappa, gravit, rair, pi, pref_mid, clat, pmid, temp, q, phis, &
         shfrc, cmfmc, cmfmc2, landfrac, snowh, sst, ps, ocnfrac, cldice, &
         qs, esl, esi, cloud, rhcloud, clc, concld, cldst, rhu00, &
         icecldf, liqcldf, relhum, shallowcu, deepcu, cold_sst_index)

    if (cold_sst_index > 0) then
       write(iulog,*) 'COLDSST: encountered in cldfrc:', lchnk, cold_sst_index, &
            ocnfrac(cold_sst_index), sst(cold_sst_index)
    end if

    call outfld('SH_CLD  ', shallowcu, pcols, lchnk)
    call outfld('DP_CLD  ', deepcu, pcols, lchnk)

    call t_stopf('ap_cldfrc_run')
  end subroutine cldfrc

!================================================================================================

  subroutine cldfrc_fice(ncol, t, fice, fsnow)
    use physconst, only: tmelt
    use ap_cloud_fraction_fice_scheme, only: cloud_fraction_fice_run
    use perf_mod, only: t_startf, t_stopf
    integer, intent(in) :: ncol
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(out) :: fice(pcols,pver)
    real(r8), intent(out) :: fsnow(pcols,pver)

    call t_startf('ap_cloud_fraction_fice_run')
    call cloud_fraction_fice_run(pcols, pver, ncol, t, tmelt, top_lev, &
         fice, fsnow)
    call t_stopf('ap_cloud_fraction_fice_run')
  end subroutine cldfrc_fice

  !-----------------------------------------------------------------------------
  ! Sets rhmin to a different value (rhminp) poleward of +/- 60 deg latitude and
  ! pressure levels less than cldfrc_rhminp_botmb (hPa) if cldfrc_rhminp is specified
  ! ** This is used only for special waccm/cam-chem cases with cam4 physics **
  !-----------------------------------------------------------------------------
  function relhum_min(press,lat) result(rh)
    use physconst, only: pi

    real(r8), intent(in) :: press, lat
    real(r8) :: rh

    rh = rhminh
    if (rhminp .eq. unset_r8 ) return

    if ((press .lt. cldfrc_rhminp_botmb*1.e2_r8) .and. &
        ( abs( lat*180._r8/pi ) .gt. 60._r8 ) ) then
       rh = rhminp
    endif

  end function relhum_min

end module cloud_fraction
