module modal_aero_wateruptake

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use physconst,        only: pi, rhoh2o, epsilo, tmelt
use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props
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
   use time_manager, only : is_first_step
   use cam_history, only : fieldname_len
   use tropopause, only : tropopause_find, TROP_ALG_HYBSTOB, TROP_ALG_CLIMATE
   use wv_sat_methods, only : wv_sat_get_default_idx
!-----------------------------------------------------------------------
!
! CAM specific driver for modal aerosol water uptake code.
!
! *** N.B. *** The calculation has been enabled for diagnostic mode lists
!              via optional arguments.  If the list_idx arg is present then
!              all the optional args must be present.
!
!-----------------------------------------------------------------------

   ! Arguments
   type(physics_state), target, intent(in)    :: state          ! Physics state variables
   type(physics_buffer_desc),   pointer       :: pbuf(:)        ! physics buffer

   integer,  optional,          intent(in)    :: list_idx_in
   real(r8), optional, target,  intent(in)    :: dgnumdry_m(:,:,:)
   real(r8), optional,          pointer, intent(out) :: dgnumwet_m(:,:,:)
   real(r8), optional,          pointer, intent(out) :: qaerwat_m(:,:,:)
   real(r8), optional,          pointer, intent(out) :: wetdens_m(:,:,:)

   integer :: i, k, l, m
   integer :: itim_old
   integer :: lchnk
   integer :: list_idx
   integer :: maxspec
   integer :: ncol
   integer :: nmodes
   integer :: nspec
   integer :: stat
   integer :: trop_lev(pcols)
   integer, allocatable :: nspec_mode(:)

   logical :: first_step
   logical, allocatable :: species_is_sulfate(:,:)

   character(len=3) :: trnum
   character(len=32) :: spectype
   character(len=fieldname_len+3) :: fieldname
   character(len=512) :: errmsg
   integer :: errflg

   real(r8), parameter :: tboil = 373.16_r8
   real(r8), pointer :: cldn(:,:)
   real(r8), pointer :: dgncur_a(:,:,:)
   real(r8), pointer :: dgncur_awet(:,:,:)
   real(r8), pointer :: qaerwat(:,:,:)
   real(r8), pointer :: raer(:,:)
   real(r8), pointer :: sulfeq(:,:,:)
   real(r8), pointer :: wetdens(:,:,:)
   real(r8), allocatable, target :: sulfeq_work(:,:,:)
   real(r8), allocatable :: rhcrystal(:)
   real(r8), allocatable :: rhdeliques(:)
   real(r8), allocatable :: sigmag_mode(:)
   real(r8), allocatable :: species_density(:,:)
   real(r8), allocatable :: species_hygro(:,:)
   real(r8), allocatable :: species_mmr(:,:,:,:)
   real(r8), allocatable :: sulden(:,:,:)
   real(r8), allocatable :: wtpct(:,:,:)


   call t_startf('ap_modal_aero_wateruptake_dr_run')

   lchnk = state%lchnk
   ncol = state%ncol
   list_idx = 0
   if (present(list_idx_in)) then
      list_idx = list_idx_in
      if (.not. present(dgnumdry_m) .or. .not. present(dgnumwet_m) .or. &
          .not. present(qaerwat_m) .or. .not. present(wetdens_m)) then
         call endrun('modal_aero_wateruptake_dr called for diagnostic list but required args not present')
      end if
   end if

   call rad_cnst_get_info(list_idx, nmodes=nmodes)
   allocate(nspec_mode(nmodes), sigmag_mode(nmodes), rhcrystal(nmodes), &
        rhdeliques(nmodes), stat=stat)
   if (stat /= 0) call endrun('modal_aero_wateruptake_dr: mode property allocation failure')

   maxspec = 0
   do m = 1, nmodes
      call rad_cnst_get_info(list_idx, m, nspec=nspec_mode(m))
      maxspec = max(maxspec, nspec_mode(m))
   end do

   allocate(species_density(maxspec,nmodes), species_hygro(maxspec,nmodes), &
        species_is_sulfate(maxspec,nmodes), &
        species_mmr(pcols,pver,maxspec,nmodes), &
        wtpct(pcols,pver,nmodes), sulden(pcols,pver,nmodes), stat=stat)
   if (stat /= 0) call endrun('modal_aero_wateruptake_dr: constituent allocation failure')

   species_density(:,:) = 0._r8
   species_hygro(:,:) = 0._r8
   species_is_sulfate(:,:) = .false.
   species_mmr(:,:,:,:) = 0._r8

   do m = 1, nmodes
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigmag_mode(m), &
           rhcrystal=rhcrystal(m), rhdeliques=rhdeliques(m))
      nspec = nspec_mode(m)
      do l = 1, nspec
         call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, raer)
         call rad_cnst_get_aer_props(list_idx, m, l, &
              density_aer=species_density(l,m), &
              hygro_aer=species_hygro(l,m), spectype=spectype)
         species_is_sulfate(l,m) = trim(spectype) == 'sulfate'
         do k = 1, pver
            do i = 1, pcols
               species_mmr(i,k,l,m) = raer(i,k)
            end do
         end do
      end do
   end do

   if (list_idx == 0) then
      call pbuf_get_field(pbuf, dgnum_idx, dgncur_a)
      call pbuf_get_field(pbuf, dgnumwet_idx, dgncur_awet)
      call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens)
      call pbuf_get_field(pbuf, qaerwat_idx, qaerwat)
   else
      dgncur_a => dgnumdry_m
      allocate(dgnumwet_m(pcols,pver,nmodes), qaerwat_m(pcols,pver,nmodes), &
           wetdens_m(pcols,pver,nmodes), stat=stat)
      if (stat /= 0) call endrun('modal_aero_wateruptake_dr: diagnostic output allocation failure')
      dgncur_awet => dgnumwet_m
      qaerwat => qaerwat_m
      wetdens => wetdens_m
      dgncur_awet(:,:,:) = dgncur_a(:,:,:)
   end if

   if (modal_strat_sulfate) then
      call pbuf_get_field(pbuf, sulfeq_idx, sulfeq)
      call tropopause_find(state, trop_lev, primary=TROP_ALG_HYBSTOB, &
           backup=TROP_ALG_CLIMATE)
   else
      allocate(sulfeq_work(pcols,pver,nmodes), stat=stat)
      if (stat /= 0) call endrun('modal_aero_wateruptake_dr: sulfate workspace allocation failure')
      sulfeq_work(:,:,:) = 0._r8
      sulfeq => sulfeq_work
      trop_lev(:) = 0
   end if

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), &
        kount=(/pcols,pver,1/))
   first_step = is_first_step()

   call modal_aero_wateruptake_dr_run( &
        ncol, pcols, pver, nmodes, maxspec, top_lev, &
        modal_strat_sulfate, first_step, wv_sat_get_default_idx(), &
        pi, pi43, rhoh2o, epsilo, tmelt, tboil, nspec_mode, sigmag_mode, &
        rhcrystal, rhdeliques, species_density, species_hygro, &
        species_is_sulfate, species_mmr, state%q(:,:,1), state%t, &
        state%pmid, cldn, trop_lev, dgncur_a, dgncur_awet, sulfeq, &
        wetdens, qaerwat, wtpct, sulden, errmsg, errflg)

   if (errflg /= 0) then
      write(iulog,*) trim(errmsg)
      call endrun('modal_aero_wateruptake_dr: standalone kernel failure')
   end if

   if (modal_strat_sulfate) then
      do m = 1, nmodes
         fieldname = ' '
         write(fieldname,fmt='(a,i1)') 'wtpct_a',m
         call outfld(fieldname, wtpct(1:ncol,1:pver,m), ncol, lchnk)
         fieldname = ' '
         write(fieldname,fmt='(a,i1)') 'sulfeq_a',m
         call outfld(fieldname, sulfeq(1:ncol,1:pver,m), ncol, lchnk)
         fieldname = ' '
         write(fieldname,fmt='(a,i1)') 'sulden_a',m
         call outfld(fieldname, sulden(1:ncol,1:pver,m), ncol, lchnk)
      end do
   end if

   if (list_idx == 0) then
      do m = 1, nmodes
         write(trnum, '(i3.3)') m
         call outfld('wat_a'//trnum(3:3), qaerwat(:,:,m), pcols, lchnk)
         call outfld('dgnd_a'//trnum(2:3), dgncur_a(:,:,m), pcols, lchnk)
         call outfld('dgnw_a'//trnum(2:3), dgncur_awet(:,:,m), pcols, lchnk)
      end do
   end if

   deallocate(nspec_mode, sigmag_mode, rhcrystal, rhdeliques, &
        species_density, species_hygro, species_is_sulfate, species_mmr, &
        wtpct, sulden)
   if (allocated(sulfeq_work)) deallocate(sulfeq_work)

   call t_stopf('ap_modal_aero_wateruptake_dr_run')

end subroutine modal_aero_wateruptake_dr


!----------------------------------------------------------------------

   end module modal_aero_wateruptake
