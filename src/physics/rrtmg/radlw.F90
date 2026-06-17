
module radlw
!----------------------------------------------------------------------- 
! 
! Purpose: Longwave radiation calculations.
!
!-----------------------------------------------------------------------
use shr_kind_mod,      only: r8 => shr_kind_r8
use ppgrid,            only: pcols, pver, pverp
use scamMod,           only: single_column, scm_crm_mode
use parrrtm,           only: nbndlw, ngptlw
use rrtmg_lw_init,     only: rrtmg_lw_ini
use rrtmg_lw_rad,      only: rrtmg_lw
use spmd_utils,        only: masterproc
use perf_mod,          only: t_startf, t_stopf
use cam_logfile,       only: iulog
use cam_abortutils,    only: endrun
use radconstants,      only: nlwbands
use iso_c_binding,     only: c_int64_t

implicit none

private
save

! Public methods

public ::&
   radlw_init,   &! initialize constants
   rad_rrtmg_lw   ! driver for longwave radiation code
   
! Private data
integer :: ntoplw    ! top level to solve for longwave cooling
logical :: use_native_rrtmg_lw_driver_impl = .false.
logical :: rrtmg_lw_driver_impl_selected = .false.
logical :: rrtmg_lw_driver_entered_logged = .false.

!===============================================================================
CONTAINS
!===============================================================================

subroutine rad_rrtmg_lw(lchnk   ,ncol      ,rrtmg_levs,r_state,       &
                        pmid    ,aer_lw_abs,cld       ,tauc_lw,       &
                        qrl     ,qrlc      ,                          &
                        flns    ,flnt      ,flnsc     ,flntc  ,flwds, &
                        flut    ,flutc     ,fnl       ,fcnl   ,fldsc, &
                        lu      ,ld        )

!-----------------------------------------------------------------------
   use cam_history,         only: outfld
   use mcica_subcol_gen_lw, only: mcica_subcol_lw
   use physconst,           only: cpair
   use rrtmg_state,         only: rrtmg_state_t
   use iso_c_binding,       only: c_double, c_int64_t, c_loc, c_ptr

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: rrtmg_levs            ! number of levels rad is applied

!
! Input arguments which are only passed to other routines
!
    type(rrtmg_state_t), intent(in), target :: r_state

   real(r8), intent(in) :: pmid(pcols,pver)     ! Level pressure (Pascals)

   real(r8), target, intent(in) :: aer_lw_abs (pcols,pver,nbndlw) ! aerosol absorption optics depth (LW)

   real(r8), intent(in) :: cld(pcols,pver)      ! Cloud cover
   real(r8), intent(in) :: tauc_lw(nbndlw,pcols,pver)   ! Cloud longwave optical depth by band

!
! Output arguments
!
   real(r8), target, intent(out) :: qrl (pcols,pver)     ! Longwave heating rate
   real(r8), target, intent(out) :: qrlc(pcols,pver)     ! Clearsky longwave heating rate
   real(r8), target, intent(out) :: flns(pcols)          ! Surface cooling flux
   real(r8), target, intent(out) :: flnt(pcols)          ! Net outgoing flux
   real(r8), target, intent(out) :: flut(pcols)          ! Upward flux at top of model
   real(r8), target, intent(out) :: flnsc(pcols)         ! Clear sky surface cooing
   real(r8), target, intent(out) :: flntc(pcols)         ! Net clear sky outgoing flux
   real(r8), target, intent(out) :: flutc(pcols)         ! Upward clear-sky flux at top of model
   real(r8), target, intent(out) :: flwds(pcols)         ! Down longwave flux at surface
   real(r8), target, intent(out) :: fldsc(pcols)         ! Down longwave clear flux at surface
   real(r8), target, intent(out) :: fcnl(pcols,pverp)    ! clear sky net flux at interfaces
   real(r8), target, intent(out) :: fnl(pcols,pverp)     ! net flux at interfaces

   real(r8), pointer, dimension(:,:,:) :: lu ! longwave spectral flux up
   real(r8), pointer, dimension(:,:,:) :: ld ! longwave spectral flux down
   
!
!---------------------------Local variables-----------------------------
!
   integer :: i, k, kk, nbnd         ! indices

   real(r8), target :: ful(pcols,pverp)     ! Total upwards longwave flux
   real(r8), target :: fsul(pcols,pverp)    ! Clear sky upwards longwave flux
   real(r8), target :: fdl(pcols,pverp)     ! Total downwards longwave flux
   real(r8), target :: fsdl(pcols,pverp)    ! Clear sky downwards longwv flux

   integer :: inflglw               ! Flag for cloud parameterization method
   integer :: iceflglw              ! Flag for ice cloud param method
   integer :: liqflglw              ! Flag for liquid cloud param method
   integer :: icld                  ! Flag for cloud overlap method
                                 ! 0=clear, 1=random, 2=maximum/random, 3=maximum

   real(r8), target :: tsfc(pcols)          ! surface temperature
   real(r8), target :: emis(pcols,nbndlw)   ! surface emissivity

   real(r8), target :: taua_lw(pcols,rrtmg_levs-1,nbndlw)     ! aerosol optical depth by band

   real(r8), parameter :: dps = 1._r8/86400._r8 ! Inverse of seconds per day

   ! Cloud arrays for McICA 
   integer, parameter :: nsubclw = ngptlw       ! rrtmg_lw g-point (quadrature point) dimension
   integer :: permuteseed                       ! permute seed for sub-column generator

   real(r8), target :: cicewp(pcols,rrtmg_levs-1)   ! in-cloud cloud ice water path
   real(r8), target :: cliqwp(pcols,rrtmg_levs-1)   ! in-cloud cloud liquid water path
   real(r8), target :: rei(pcols,rrtmg_levs-1)      ! ice particle effective radius (microns)
   real(r8), target :: rel(pcols,rrtmg_levs-1)      ! liquid particle radius (micron)

   real(r8) :: cld_stolw(nsubclw, pcols, rrtmg_levs-1)     ! cloud fraction (mcica)
   real(r8) :: cicewp_stolw(nsubclw, pcols, rrtmg_levs-1)  ! cloud ice water path (mcica)
   real(r8) :: cliqwp_stolw(nsubclw, pcols, rrtmg_levs-1)  ! cloud liquid water path (mcica)
   real(r8) :: rei_stolw(pcols,rrtmg_levs-1)               ! ice particle size (mcica)
   real(r8) :: rel_stolw(pcols,rrtmg_levs-1)               ! liquid particle size (mcica)
   real(r8) :: tauc_stolw(nsubclw, pcols, rrtmg_levs-1)    ! cloud optical depth (mcica - optional)

   ! Includes extra layer above model top
   real(r8), target :: uflx(pcols,rrtmg_levs+1)  ! Total upwards longwave flux
   real(r8), target :: uflxc(pcols,rrtmg_levs+1) ! Clear sky upwards longwave flux
   real(r8), target :: dflx(pcols,rrtmg_levs+1)  ! Total downwards longwave flux
   real(r8), target :: dflxc(pcols,rrtmg_levs+1) ! Clear sky downwards longwv flux
   real(r8), target :: hr(pcols,rrtmg_levs)      ! Longwave heating rate (K/d)
   real(r8), target :: hrc(pcols,rrtmg_levs)     ! Clear sky longwave heating rate (K/d)
   real(r8) lwuflxs(nbndlw,pcols,pverp+1)  ! Longwave spectral flux up
   real(r8) lwdflxs(nbndlw,pcols,pverp+1)  ! Longwave spectral flux down
   interface
      subroutine rrtmg_lw_zero_cloud_inputs_codon(ncol_c, pcols_c, nlay_c, &
           cicewp_p, cliqwp_p, rei_p, rel_p) bind(c, name="rrtmg_lw_zero_cloud_inputs_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, nlay_c
         type(c_ptr), value :: cicewp_p, cliqwp_p, rei_p, rel_p
      end subroutine rrtmg_lw_zero_cloud_inputs_codon
      subroutine rrtmg_lw_pre_codon(ncol_c, pcols_c, pver_c, pverp_c, rrtmg_levs_c, nbndlw_c, &
           aer_lw_abs_p, tlev_p, emis_p, tsfc_p, taua_lw_p) bind(c, name="rrtmg_lw_pre_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c, rrtmg_levs_c, nbndlw_c
         type(c_ptr), value :: aer_lw_abs_p, tlev_p, emis_p, tsfc_p, taua_lw_p
      end subroutine rrtmg_lw_pre_codon
      subroutine rrtmg_lw_post_codon(ncol_c, pcols_c, pver_c, pverp_c, rrtmg_levs_c, ntoplw_c, cpair_c, &
           uflx_p, dflx_p, hr_p, uflxc_p, dflxc_p, hrc_p, flwds_p, fldsc_p, flns_p, flnsc_p, flnt_p, flntc_p, &
           flut_p, flutc_p, ful_p, fdl_p, fsul_p, fsdl_p, fnl_p, fcnl_p, qrl_p, qrlc_p) &
           bind(c, name="rrtmg_lw_post_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c, rrtmg_levs_c, ntoplw_c
         real(c_double), value :: cpair_c
         type(c_ptr), value :: uflx_p, dflx_p, hr_p, uflxc_p, dflxc_p, hrc_p
         type(c_ptr), value :: flwds_p, fldsc_p, flns_p, flnsc_p, flnt_p, flntc_p, flut_p, flutc_p
         type(c_ptr), value :: ful_p, fdl_p, fsul_p, fsdl_p, fnl_p, fcnl_p, qrl_p, qrlc_p
      end subroutine rrtmg_lw_post_codon
   end interface
   !-----------------------------------------------------------------------

   ! mji/rrtmg

   ! Calculate cloud optical properties here if using CAM method, or if using one of the
   ! methods in RRTMG_LW, then pass in cloud physical properties and zero out cloud optical 
   ! properties here
   
   ! Zero optional cloud optical depth input array tauc_lw, 
   ! if inputting cloud physical properties into RRTMG_LW
   !          tauc_lw(:,:,:) = 0.
   ! Or, pass in CAM cloud longwave optical depth to RRTMG_LW
   ! do nbnd = 1, nbndlw
   !    tauc_lw(nbnd,:ncol,:pver) = cldtau(:ncol,:pver)
   ! end do

   ! Call mcica sub-column generator for RRTMG_LW

   ! Call sub-column generator for McICA in radiation
   call t_startf('mcica_subcol_lw')

   ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
   icld = 2
   ! Set permute seed (must be offset between LW and SW by at least 140 to insure 
   ! effective randomization)
   permuteseed = 150

   ! These fields are no longer supplied by CAM.
   call rrtmg_lw_driver_select_impl()
   if (use_native_rrtmg_lw_driver_impl) then
      cicewp = 0.0_r8
      cliqwp = 0.0_r8
      rei = 0.0_r8
      rel = 0.0_r8
   else
      call rrtmg_lw_driver_log_entered()
      call rrtmg_lw_zero_cloud_inputs_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(rrtmg_levs-1, c_int64_t), &
           c_loc(cicewp(1,1)), c_loc(cliqwp(1,1)), c_loc(rei(1,1)), c_loc(rel(1,1)) &
      )
   end if

   call mcica_subcol_lw(lchnk, ncol, rrtmg_levs-1, icld, permuteseed, pmid(:, pverp-rrtmg_levs+1:pverp-1), &
      cld(:, pverp-rrtmg_levs+1:pverp-1), cicewp, cliqwp, rei, rel, tauc_lw(:, :ncol, pverp-rrtmg_levs+1:pverp-1), &
      cld_stolw, cicewp_stolw, cliqwp_stolw, rei_stolw, rel_stolw, tauc_stolw)

   call t_stopf('mcica_subcol_lw')

   
   call t_startf('rrtmg_lw')

   !
   ! Call RRTMG_LW model
   !
   ! Set input flags for cloud parameterizations
   ! Use separate specification of ice and liquid cloud optical depth.
   ! Use either Ebert and Curry ice parameterization (iceflglw = 0 or 1), 
   ! or use Key (Streamer) approach (iceflglw = 2), or use Fu method
   ! (iceflglw = 3), and Hu/Stamnes for liquid (liqflglw = 1).
   ! For use in Fu method (iceflglw = 3), rei is converted in RRTMG_LW
   ! from effective radius to generalized effective size using the
   ! conversion of D. Mitchell, JAS, 2002.  For ice particles outside
   ! the effective range of either the Key or Fu approaches, the 
   ! Ebert and Curry method is applied. 

   ! Input CAM cloud optical depth directly
   inflglw = 0
   iceflglw = 0
   liqflglw = 0
   ! Use E&C approach for ice to mimic CAM3
   !   inflglw = 2
   !   iceflglw = 1
   !   liqflglw = 1
   ! Use merged Fu and E&C params for ice
   !   inflglw = 2
   !   iceflglw = 3
   !   liqflglw = 1

   ! Convert incoming water amounts from specific humidity to vmr as needed;
   ! Convert other incoming molecular amounts from mmr to vmr as needed;
   ! Convert pressures from Pa to hPa;
   ! Set surface emissivity to 1.0 here, this is treated in land surface model;
   ! Set surface temperature
   ! Set aerosol optical depth to zero for now

   if (associated(lu)) lu(1:ncol,:,:) = 0.0_r8
   if (associated(ld)) ld(1:ncol,:,:) = 0.0_r8

   if (use_native_rrtmg_lw_driver_impl) then
      emis(:ncol,:nbndlw) = 1._r8
      tsfc(:ncol) = r_state%tlev(:ncol,rrtmg_levs+1)
      taua_lw(:ncol, 1:rrtmg_levs-1, :nbndlw) = aer_lw_abs(:ncol,pverp-rrtmg_levs+1:pverp-1,:nbndlw)
   else
      call rrtmg_lw_driver_log_entered()
      call rrtmg_lw_pre_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), &
           int(rrtmg_levs, c_int64_t), int(nbndlw, c_int64_t), &
           c_loc(aer_lw_abs(1,1,1)), c_loc(r_state%tlev(1,1)), c_loc(emis(1,1)), c_loc(tsfc(1)), c_loc(taua_lw(1,1,1)) &
      )
   end if

   call rrtmg_lw(lchnk  ,ncol ,rrtmg_levs    ,icld    ,                 &
        r_state%pmidmb  ,r_state%pintmb  ,r_state%tlay    ,r_state%tlev    ,tsfc    ,r_state%h2ovmr, &
        r_state%o3vmr   ,r_state%co2vmr  ,r_state%ch4vmr  ,r_state%o2vmr   ,r_state%n2ovmr  ,r_state%cfc11vmr,r_state%cfc12vmr, &
        r_state%cfc22vmr,r_state%ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
        cld_stolw,tauc_stolw,cicewp_stolw,cliqwp_stolw ,rei, rel, &
        taua_lw, &
        uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, &
        lwuflxs, lwdflxs)

   if (use_native_rrtmg_lw_driver_impl) then
      !
      !----------------------------------------------------------------------
      ! All longitudes: store history tape quantities
      ! Flux units are in W/m2 on output from rrtmg_lw and contain output for
      ! extra layer above model top with vertical indexing from bottom to top.
      ! Heating units are in K/d on output from RRTMG and contain output for
      ! extra layer above model top with vertical indexing from bottom to top.
      ! Heating units are converted to J/kg/s below for use in CAM.

      flwds(:ncol) = dflx (:ncol,1)
      fldsc(:ncol) = dflxc(:ncol,1)
      flns(:ncol)  = uflx (:ncol,1) - dflx (:ncol,1)
      flnsc(:ncol) = uflxc(:ncol,1) - dflxc(:ncol,1)
      flnt(:ncol)  = uflx (:ncol,rrtmg_levs) - dflx (:ncol,rrtmg_levs)
      flntc(:ncol) = uflxc(:ncol,rrtmg_levs) - dflxc(:ncol,rrtmg_levs)
      flut(:ncol)  = uflx (:ncol,rrtmg_levs)
      flutc(:ncol) = uflxc(:ncol,rrtmg_levs)

      !
      ! Reverse vertical indexing here for CAM arrays to go from top to bottom.
      !
      ful = 0._r8
      fdl = 0._r8
      fsul = 0._r8
      fsdl = 0._r8
      ful (:ncol,pverp-rrtmg_levs+1:pverp)= uflx(:ncol,rrtmg_levs:1:-1)
      fdl (:ncol,pverp-rrtmg_levs+1:pverp)= dflx(:ncol,rrtmg_levs:1:-1)
      fsul(:ncol,pverp-rrtmg_levs+1:pverp)=uflxc(:ncol,rrtmg_levs:1:-1)
      fsdl(:ncol,pverp-rrtmg_levs+1:pverp)=dflxc(:ncol,rrtmg_levs:1:-1)

      fnl(:ncol,:) = ful(:ncol,:) - fdl(:ncol,:)
      ! mji/ cam excluded this?
      fcnl(:ncol,:) = fsul(:ncol,:) - fsdl(:ncol,:)

      ! Pass longwave heating to CAM arrays and convert from K/d to J/kg/s
      qrl = 0._r8
      qrlc = 0._r8
      qrl (:ncol,pverp-rrtmg_levs+1:pver)=hr (:ncol,rrtmg_levs-1:1:-1)*cpair*dps
      qrlc(:ncol,pverp-rrtmg_levs+1:pver)=hrc(:ncol,rrtmg_levs-1:1:-1)*cpair*dps

      ! Return 0 above solution domain
      if ( ntoplw > 1 )then
         qrl(:ncol,:ntoplw-1) = 0._r8
         qrlc(:ncol,:ntoplw-1) = 0._r8
      end if
   else
      call rrtmg_lw_post_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), &
           int(rrtmg_levs, c_int64_t), int(ntoplw, c_int64_t), real(cpair, c_double), &
           c_loc(uflx(1,1)), c_loc(dflx(1,1)), c_loc(hr(1,1)), c_loc(uflxc(1,1)), c_loc(dflxc(1,1)), c_loc(hrc(1,1)), &
           c_loc(flwds(1)), c_loc(fldsc(1)), c_loc(flns(1)), c_loc(flnsc(1)), c_loc(flnt(1)), c_loc(flntc(1)), &
           c_loc(flut(1)), c_loc(flutc(1)), c_loc(ful(1,1)), c_loc(fdl(1,1)), c_loc(fsul(1,1)), c_loc(fsdl(1,1)), &
           c_loc(fnl(1,1)), c_loc(fcnl(1,1)), c_loc(qrl(1,1)), c_loc(qrlc(1,1)) &
      )
   end if

   if (single_column.and.scm_crm_mode) then
      call outfld('FUL     ',ful,pcols,lchnk)
      call outfld('FDL     ',fdl,pcols,lchnk)
      call outfld('FULC    ',fsul,pcols,lchnk)
      call outfld('FDLC    ',fsdl,pcols,lchnk)
   endif

   ! Pass spectral fluxes, reverse layering
   ! order=(/3,1,2/) maps the first index of lwuflxs to the third index of lu.
   if (associated(lu)) then
      lu(:ncol,pverp-rrtmg_levs+1:pverp,:) = reshape(lwuflxs(:,:ncol,rrtmg_levs:1:-1), &
           (/ncol,rrtmg_levs,nbndlw/), order=(/3,1,2/))
   end if
   
   if (associated(ld)) then
      ld(:ncol,pverp-rrtmg_levs+1:pverp,:) = reshape(lwdflxs(:,:ncol,rrtmg_levs:1:-1), &
           (/ncol,rrtmg_levs,nbndlw/), order=(/3,1,2/))
   end if
   
   call t_stopf('rrtmg_lw')

end subroutine rad_rrtmg_lw

!-------------------------------------------------------------------------------

subroutine radlw_init()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize various constants for radiation scheme.
!
!-----------------------------------------------------------------------

   use ref_pres, only : pref_mid

   integer :: k

#define CAM_MISC_TAG 348
#define CAM_MISC_LABEL 'radlw_init'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG

   ! If the top model level is above ~90 km (0.1 Pa), set the top level to compute
   ! longwave cooling to about 80 km (1 Pa)
   if (pref_mid(1) .lt. 0.1_r8) then
      do k = 1, pver
         if (pref_mid(k) .lt. 1._r8) ntoplw  = k
      end do
   else
      ntoplw  = 1
   end if
   if (masterproc) then
      write(iulog,*) 'radlw_init: ntoplw =',ntoplw
   endif

   call rrtmg_lw_ini

end subroutine radlw_init

!-------------------------------------------------------------------------------

subroutine rrtmg_lw_driver_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (rrtmg_lw_driver_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('RRTMG_LW_DRIVER_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_rrtmg_lw_driver_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_rrtmg_lw_driver_impl = .false.
   end if

   rrtmg_lw_driver_impl_selected = .true.

   if (masterproc) then
      if (use_native_rrtmg_lw_driver_impl) then
         write(iulog,*) 'rrtmg_lw_driver implementation = native'
      else
         write(iulog,*) 'rrtmg_lw_driver implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine rrtmg_lw_driver_select_impl

!-------------------------------------------------------------------------------

subroutine rrtmg_lw_driver_log_entered()

   if (rrtmg_lw_driver_entered_logged) return
   rrtmg_lw_driver_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'rrtmg_lw_driver entered (cloud input zero/pre/post helpers = codon; ' // &
           'native pre/post blocks skipped; rrtmg_lw core = native)'
      call flush(iulog)
   end if

end subroutine rrtmg_lw_driver_log_entered

!-------------------------------------------------------------------------------

end module radlw
