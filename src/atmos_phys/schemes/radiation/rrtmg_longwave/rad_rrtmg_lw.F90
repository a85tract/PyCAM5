
module ap_rad_rrtmg_lw_scheme
!-----------------------------------------------------------------------
!
! Purpose: Longwave radiation calculations.
!
!-----------------------------------------------------------------------
use shr_kind_mod,      only: r8 => shr_kind_r8
use parrrtm,           only: nbndlw, ngptlw
use rrtmg_lw_rad,      only: rrtmg_lw
use radiation_host_hooks, only: radiation_timer_start, radiation_timer_stop, &
                                radiation_outfld_real2d

implicit none

private
save

! Public methods

public ::&
   rad_rrtmg_lw_run   ! driver for longwave radiation code

!===============================================================================
CONTAINS
!===============================================================================

!> \section arg_table_rad_rrtmg_lw_run Argument Table
!! \htmlinclude rad_rrtmg_lw_run.html
subroutine rad_rrtmg_lw_run(pcols, pver, pverp, lchnk, ncol, rrtmg_levs, ntoplw, cpair, &
                        single_column, scm_crm_mode, &
                        pmidmb, pintmb, tlay, tlev, h2ovmr, o3vmr, co2vmr, ch4vmr, &
                        o2vmr, n2ovmr, cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, &
                        pmid    ,aer_lw_abs,cld       ,tauc_lw,       &
                        qrl     ,qrlc      ,                          &
                        flns    ,flnt      ,flnsc     ,flntc  ,flwds, &
                        flut    ,flutc     ,fnl       ,fcnl   ,fldsc, &
                        lu      ,ld, ful, fsul, fdl, fsdl)

!-----------------------------------------------------------------------
   use mcica_subcol_gen_lw, only: mcica_subcol_lw

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: pcols, pver, pverp
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: rrtmg_levs            ! number of levels rad is applied
   integer, intent(in) :: ntoplw                 ! top level to solve for longwave cooling
   real(r8), intent(in) :: cpair                 ! dry-air heat capacity
   logical, intent(in) :: single_column, scm_crm_mode

!
! Input arguments which are only passed to other routines
!
   real(r8), intent(in) :: pmidmb(:,:), pintmb(:,:)
   real(r8), intent(in) :: tlay(:,:), tlev(:,:)
   real(r8), intent(in) :: h2ovmr(:,:), o3vmr(:,:), co2vmr(:,:), ch4vmr(:,:)
   real(r8), intent(in) :: o2vmr(:,:), n2ovmr(:,:), cfc11vmr(:,:), cfc12vmr(:,:)
   real(r8), intent(in) :: cfc22vmr(:,:), ccl4vmr(:,:)

   real(r8), intent(in) :: pmid(:,:)             ! Level pressure (Pascals)

   real(r8), intent(in) :: aer_lw_abs(:,:,:) ! aerosol absorption optics depth (LW)

   real(r8), intent(in) :: cld(:,:)              ! Cloud cover
   real(r8), intent(in) :: tauc_lw(:,:,:)        ! Cloud longwave optical depth by band

!
! Output arguments
!
   real(r8), intent(out) :: qrl(:,:)              ! Longwave heating rate
   real(r8), intent(out) :: qrlc(:,:)             ! Clearsky longwave heating rate
   real(r8), intent(out) :: flns(:)               ! Surface cooling flux
   real(r8), intent(out) :: flnt(:)               ! Net outgoing flux
   real(r8), intent(out) :: flut(:)               ! Upward flux at top of model
   real(r8), intent(out) :: flnsc(:)              ! Clear sky surface cooing
   real(r8), intent(out) :: flntc(:)              ! Net clear sky outgoing flux
   real(r8), intent(out) :: flutc(:)              ! Upward clear-sky flux at top of model
   real(r8), intent(out) :: flwds(:)              ! Down longwave flux at surface
   real(r8), intent(out) :: fldsc(:)              ! Down longwave clear flux at surface
   real(r8), intent(out) :: fcnl(:,:)             ! clear sky net flux at interfaces
   real(r8), intent(out) :: fnl(:,:)              ! net flux at interfaces

   real(r8), pointer, intent(inout), dimension(:,:,:) :: lu ! longwave spectral flux up
   real(r8), pointer, intent(inout), dimension(:,:,:) :: ld ! longwave spectral flux down

!
!---------------------------Local variables-----------------------------
!
   integer :: i, k, kk, nbnd         ! indices

   real(r8), intent(out) :: ful(pcols,pverp)     ! Total upwards longwave flux
   real(r8), intent(out) :: fsul(pcols,pverp)    ! Clear sky upwards longwave flux
   real(r8), intent(out) :: fdl(pcols,pverp)     ! Total downwards longwave flux
   real(r8), intent(out) :: fsdl(pcols,pverp)    ! Clear sky downwards longwv flux

   integer :: inflglw               ! Flag for cloud parameterization method
   integer :: iceflglw              ! Flag for ice cloud param method
   integer :: liqflglw              ! Flag for liquid cloud param method
   integer :: icld                  ! Flag for cloud overlap method
                                 ! 0=clear, 1=random, 2=maximum/random, 3=maximum

   real(r8) :: tsfc(pcols)          ! surface temperature
   real(r8) :: emis(pcols,nbndlw)   ! surface emissivity

   real(r8) :: taua_lw(pcols,rrtmg_levs-1,nbndlw)     ! aerosol optical depth by band

   real(r8), parameter :: dps = 1._r8/86400._r8 ! Inverse of seconds per day

   ! Cloud arrays for McICA
   integer, parameter :: nsubclw = ngptlw       ! rrtmg_lw g-point (quadrature point) dimension
   integer :: permuteseed                       ! permute seed for sub-column generator

   real(r8) :: cicewp(pcols,rrtmg_levs-1)   ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,rrtmg_levs-1)   ! in-cloud cloud liquid water path
   real(r8) :: rei(pcols,rrtmg_levs-1)      ! ice particle effective radius (microns)
   real(r8) :: rel(pcols,rrtmg_levs-1)      ! liquid particle radius (micron)

   real(r8) :: cld_stolw(nsubclw, pcols, rrtmg_levs-1)     ! cloud fraction (mcica)
   real(r8) :: cicewp_stolw(nsubclw, pcols, rrtmg_levs-1)  ! cloud ice water path (mcica)
   real(r8) :: cliqwp_stolw(nsubclw, pcols, rrtmg_levs-1)  ! cloud liquid water path (mcica)
   real(r8) :: rei_stolw(pcols,rrtmg_levs-1)               ! ice particle size (mcica)
   real(r8) :: rel_stolw(pcols,rrtmg_levs-1)               ! liquid particle size (mcica)
   real(r8) :: tauc_stolw(nsubclw, pcols, rrtmg_levs-1)    ! cloud optical depth (mcica - optional)

   ! Includes extra layer above model top
   real(r8) :: uflx(pcols,rrtmg_levs+1)  ! Total upwards longwave flux
   real(r8) :: uflxc(pcols,rrtmg_levs+1) ! Clear sky upwards longwave flux
   real(r8) :: dflx(pcols,rrtmg_levs+1)  ! Total downwards longwave flux
   real(r8) :: dflxc(pcols,rrtmg_levs+1) ! Clear sky downwards longwv flux
   real(r8) :: hr(pcols,rrtmg_levs)      ! Longwave heating rate (K/d)
   real(r8) :: hrc(pcols,rrtmg_levs)     ! Clear sky longwave heating rate (K/d)
   real(r8) lwuflxs(nbndlw,pcols,pverp+1)  ! Longwave spectral flux up
   real(r8) lwdflxs(nbndlw,pcols,pverp+1)  ! Longwave spectral flux down
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
   call radiation_timer_start('mcica_subcol_lw')

   ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
   icld = 2
   ! Set permute seed (must be offset between LW and SW by at least 140 to insure
   ! effective randomization)
   permuteseed = 150

   ! These fields are no longer supplied by CAM.
   cicewp = 0.0_r8
   cliqwp = 0.0_r8
   rei = 0.0_r8
   rel = 0.0_r8

   call mcica_subcol_lw(lchnk, ncol, rrtmg_levs-1, icld, permuteseed, pmid(:, pverp-rrtmg_levs+1:pverp-1), &
      cld(:, pverp-rrtmg_levs+1:pverp-1), cicewp, cliqwp, rei, rel, tauc_lw(:, :ncol, pverp-rrtmg_levs+1:pverp-1), &
      cld_stolw, cicewp_stolw, cliqwp_stolw, rei_stolw, rel_stolw, tauc_stolw)

   call radiation_timer_stop('mcica_subcol_lw')
   call radiation_timer_start('rrtmg_lw')

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

   emis(:ncol,:nbndlw) = 1._r8
   tsfc(:ncol) = tlev(:ncol,rrtmg_levs+1)
   taua_lw(:ncol, 1:rrtmg_levs-1, :nbndlw) = aer_lw_abs(:ncol,pverp-rrtmg_levs+1:pverp-1,:nbndlw)

   if (associated(lu)) lu(1:ncol,:,:) = 0.0_r8
   if (associated(ld)) ld(1:ncol,:,:) = 0.0_r8

   call rrtmg_lw(lchnk  ,ncol ,rrtmg_levs    ,icld    ,                 &
        pmidmb  ,pintmb  ,tlay    ,tlev    ,tsfc    ,h2ovmr, &
        o3vmr   ,co2vmr  ,ch4vmr  ,o2vmr   ,n2ovmr  ,cfc11vmr,cfc12vmr, &
        cfc22vmr,ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
        cld_stolw,tauc_stolw,cicewp_stolw,cliqwp_stolw ,rei, rel, &
        taua_lw, &
        uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, &
        lwuflxs, lwdflxs)

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

   if (single_column .and. scm_crm_mode) then
      call radiation_outfld_real2d('FUL     ', ful, pcols, lchnk)
      call radiation_outfld_real2d('FDL     ', fdl, pcols, lchnk)
      call radiation_outfld_real2d('FULC    ', fsul, pcols, lchnk)
      call radiation_outfld_real2d('FDLC    ', fsdl, pcols, lchnk)
   end if

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

   call radiation_timer_stop('rrtmg_lw')

end subroutine rad_rrtmg_lw_run

end module ap_rad_rrtmg_lw_scheme
