module modal_aero_wateruptake

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use physconst,        only: pi, rhoh2o
use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use wv_saturation,    only: qsat_water
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props, rad_cnst_get_mode_num
use cam_history,      only: addfld, add_default, phys_decomp, outfld
use cam_logfile,      only: iulog
use ref_pres,         only: top_lev => clim_modal_aero_top_lev
use phys_control,     only: phys_getopts
use cam_abortutils,   only: endrun

implicit none
private
save

public :: &
   modal_aero_wateruptake_init, &
   modal_aero_wateruptake_dr

public :: modal_aero_wateruptake_reg

real(r8), parameter :: third = 1._r8/3._r8
real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8


! Physics buffer indices
integer :: cld_idx        = 0
integer :: dgnum_idx      = 0
integer :: dgnumwet_idx   = 0
integer :: sulfeq_idx     = 0
integer :: wetdens_ap_idx = 0
integer :: qaerwat_idx    = 0

logical, public :: modal_strat_sulfate = .false.   ! If .true. then MAM sulfate surface area density used in stratospheric heterogeneous chemistry

!===============================================================================
contains
!===============================================================================

subroutine modal_aero_wateruptake_reg()

  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use rad_constituents, only: rad_cnst_get_info

   integer :: nmodes
   
   call rad_cnst_get_info(0, nmodes=nmodes)
   call pbuf_add_field('DGNUMWET',   'global',  dtype_r8, (/pcols, pver, nmodes/), dgnumwet_idx)
   call pbuf_add_field('WETDENS_AP', 'physpkg', dtype_r8, (/pcols, pver, nmodes/), wetdens_ap_idx)

   ! 1st order rate for direct conversion of strat. cloud water to precip (1/s)
   call pbuf_add_field('QAERWAT',    'physpkg', dtype_r8, (/pcols, pver, nmodes/), qaerwat_idx)  
   if (modal_strat_sulfate) then
     call pbuf_add_field('MAMH2SO4EQ', 'global',  dtype_r8, (/pcols, pver, nmodes/), sulfeq_idx)
   end if

end subroutine modal_aero_wateruptake_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_wateruptake_init(pbuf2d)
   use time_manager,  only: is_first_step
   use physics_buffer,only: pbuf_set_field
   use infnan,       only : nan, assignment(=)

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   real(r8) :: real_nan

   integer :: m, nmodes
   logical :: history_aerosol      ! Output the MAM aerosol variables and tendencies

   character(len=3) :: trnum       ! used to hold mode number (as characters)
   !----------------------------------------------------------------------------

   real_nan = nan
    
   cld_idx        = pbuf_get_index('CLD')    
   dgnum_idx      = pbuf_get_index('DGNUM')    

   ! assume for now that will compute wateruptake for climate list modes only

   call rad_cnst_get_info(0, nmodes=nmodes)

   do m = 1, nmodes
      write(trnum, '(i3.3)') m
      call addfld('dgnd_a'//trnum(2:3), 'm', pver, 'A', &
         'dry dgnum, interstitial, mode '//trnum(2:3), phys_decomp)
      call addfld('dgnw_a'//trnum(2:3), 'm', pver, 'A', &
         'wet dgnum, interstitial, mode '//trnum(2:3), phys_decomp)
      call addfld('wat_a'//trnum(3:3), 'm', pver, 'A', &
         'aerosol water, interstitial, mode '//trnum(2:3), phys_decomp)
      
      ! determine default variables
      call phys_getopts(history_aerosol_out = history_aerosol)

      if (history_aerosol) then  
         call add_default('dgnd_a'//trnum(2:3), 1, ' ')
         call add_default('dgnw_a'//trnum(2:3), 1, ' ')
         call add_default('wat_a'//trnum(3:3),  1, ' ')
      endif

   end do
   
   if (is_first_step()) then
      ! initialize fields in physics buffer
      call pbuf_set_field(pbuf2d, dgnumwet_idx, 0.0_r8)
      if (modal_strat_sulfate) then
      ! initialize fields in physics buffer to NaN (not a number) 
      ! so model will crash if used before initialization
         call pbuf_set_field(pbuf2d, sulfeq_idx, real_nan)
      endif
   endif

end subroutine modal_aero_wateruptake_init

!===============================================================================


subroutine modal_aero_wateruptake_dr(state, pbuf, list_idx_in, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m)

   use ap_modal_aero_wateruptake_dr_scheme, only : modal_aero_wateruptake_dr_run
   use perf_mod, only : t_startf, t_stopf
!-----------------------------------------------------------------------
!
! CAM specific driver for modal aerosol water uptake code.
!
! *** N.B. *** The calculation has been enabled for diagnostic mode lists
!              via optional arguments.  If the list_idx arg is present then
!              all the optional args must be present.
!
!-----------------------------------------------------------------------

   use time_manager,  only: is_first_step
   use cam_history,   only: outfld, fieldname_len
   use tropopause,    only: tropopause_find, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE
   ! Arguments
   type(physics_state), target, intent(in)    :: state          ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)        ! physics buffer

   integer,  optional,          intent(in)    :: list_idx_in
   real(r8), optional, target,  intent(in)    :: dgnumdry_m(:,:,:)
   real(r8), optional,          pointer       :: dgnumwet_m(:,:,:)
   real(r8), optional,          pointer       :: qaerwat_m(:,:,:)
   real(r8), optional,          pointer       :: wetdens_m(:,:,:)


   call t_startf('ap_modal_aero_wateruptake_dr_run')
   call modal_aero_wateruptake_dr_run(state, pbuf, cld_idx, dgnum_idx, dgnumwet_idx, &
        sulfeq_idx, wetdens_ap_idx, qaerwat_idx, modal_strat_sulfate, &
        list_idx_in, dgnumdry_m, dgnumwet_m, qaerwat_m, wetdens_m)
   call t_stopf('ap_modal_aero_wateruptake_dr_run')

end subroutine modal_aero_wateruptake_dr


!----------------------------------------------------------------------

   end module modal_aero_wateruptake

