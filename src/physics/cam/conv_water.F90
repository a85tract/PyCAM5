  module conv_water

   ! --------------------------------------------------------------------- !
   ! Purpose:                                                              !
   ! Computes grid-box average liquid (and ice) from stratus and cumulus   !
   ! Just for the purposes of radiation.                                   !
   !                                                                       !
   ! Method:                                                               !
   ! Extract information about deep+shallow liquid and cloud fraction from !
   ! the physics buffer.                                                   !
   !                                                                       !
   ! Author: Rich Neale, August 2006                                       !
   !         October 2006: Allow averaging of liquid to give a linear      !
   !                       average in emissivity.                          !
   !         Andrew Gettelman October 2010  Separate module                !
   !---------------------------------------------------------------------- !

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use spmd_utils,     only: masterproc
  use ppgrid,         only: pcols, pver, pverp
  use physconst,      only: gravit, latvap, latice
  use cam_abortutils, only: endrun

  use perf_mod
  use cam_logfile,    only: iulog
  use iso_c_binding,  only: c_double, c_int64_t

  implicit none
  private
  save

  public :: &
     conv_water_readnl,   &
     conv_water_register, &
     conv_water_init,     &
     conv_water_4rad,     &
     conv_water_in_rad

! pbuf indices

  integer :: icwmrsh_idx, icwmrdp_idx, fice_idx, sh_frac_idx, dp_frac_idx, &
             ast_idx, sh_cldliq1_idx, sh_cldice1_idx, rei_idx

  integer :: ixcldice, ixcldliq

! Namelist
integer, parameter :: unset_int = huge(1)

integer  :: conv_water_in_rad = unset_int  ! 0==> No; 1==> Yes-Arithmetic average;
                                           ! 2==> Yes-Average in emissivity.
integer  :: conv_water_mode
real(r8) :: frac_limit
logical :: use_native_impl = .false.
logical :: impl_selected = .false.
logical :: conv_water_readnl_logged = .false.
logical :: conv_water_init_logged = .false.
logical :: conv_water_4rad_logged = .false.

interface
   subroutine conv_water_4rad_codon(ncol_c, pcols_c, pver_c, conv_water_mode_c, microp_is_rk_c, frac_limit_c, gravit_c, &
        pdel_p, ls_liq_p, ls_ice_p, ast_p, sh_frac_p, dp_frac_p, rei_p, dp_icwmr_p, sh_icwmr_p, fice_p, &
        totg_liq_p, totg_ice_p, conv_ice_p, conv_liq_p, tot_ice_p, tot_liq_p, totg_ice_sh_p, totg_liq_sh_p, &
        totg_ice_dp_p, totg_liq_dp_p, fresh_p, fredp_p, frecu_p, fretot_p, sh_cldliq_p, sh_cldice_p) &
        bind(c, name="conv_water_4rad_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, conv_water_mode_c, microp_is_rk_c
      real(c_double), value :: frac_limit_c, gravit_c
      type(c_ptr), value :: pdel_p, ls_liq_p, ls_ice_p, ast_p, sh_frac_p, dp_frac_p, rei_p
      type(c_ptr), value :: dp_icwmr_p, sh_icwmr_p, fice_p
      type(c_ptr), value :: totg_liq_p, totg_ice_p, conv_ice_p, conv_liq_p, tot_ice_p, tot_liq_p
      type(c_ptr), value :: totg_ice_sh_p, totg_liq_sh_p, totg_ice_dp_p, totg_liq_dp_p
      type(c_ptr), value :: fresh_p, fredp_p, frecu_p, fretot_p, sh_cldliq_p, sh_cldice_p
   end subroutine conv_water_4rad_codon
   function conv_water_readnl_codon(value_c) result(out_c) bind(c, name="conv_water_readnl_codon")
      use iso_c_binding, only: c_double
      real(c_double), value :: value_c
      real(c_double) :: out_c
   end function conv_water_readnl_codon
   function conv_water_register_codon(flag_c) result(out_c) bind(c, name="conv_water_register_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function conv_water_register_codon
   function conv_water_init_codon(flag_c) result(out_c) bind(c, name="conv_water_init_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function conv_water_init_codon
end interface

!=============================================================================================
contains
!=============================================================================================

subroutine conv_water_log_direct(logged, proof_line)

   logical, intent(inout) :: logged
   character(len=*), intent(in) :: proof_line

   if (logged) return
   logged = .true.

   if (masterproc) then
      write(iulog,'(A)') trim(proof_line)
      call flush(iulog)
   end if

end subroutine conv_water_log_direct

!=============================================================================================

subroutine conv_water_4rad_log_direct()

   if (conv_water_4rad_logged) return
   conv_water_4rad_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'conv_water_4rad direct = codon; pbuf/outfld native CAM API island'
      call flush(iulog)
   end if

end subroutine conv_water_4rad_log_direct

subroutine conv_water_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CONV_WATER_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_impl = .false.
   end if

   impl_selected = .true.

   if (masterproc) then
      if (use_native_impl) then
         write(iulog,*) 'conv_water implementation = native'
      else
         write(iulog,*) 'conv_water implementation = codon'
      end if
   end if

end subroutine conv_water_select_impl

subroutine conv_water_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'conv_water_readnl'

   real(r8) :: conv_water_frac_limit

   namelist /conv_water_nl/ conv_water_in_rad, conv_water_frac_limit
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'conv_water_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, conv_water_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(conv_water_in_rad,     1, mpiint, 0, mpicom)
   call mpibcast(conv_water_frac_limit, 1, mpir8,  0, mpicom)
#endif

   conv_water_mode = conv_water_in_rad
   frac_limit      = conv_water_readnl_codon(real(conv_water_frac_limit, c_double))
   call conv_water_log_direct(conv_water_readnl_logged, &
        'conv_water_readnl direct = codon; namelist I/O and MPI broadcast native islands')

end subroutine conv_water_readnl

!=============================================================================================

  subroutine conv_water_register

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Register the fields in the physics buffer.                            !
  !                                                                       !
  !---------------------------------------------------------------------- !

    use constituents, only: cnst_add, pcnst
    use physconst,    only: mwdry, cpair

    use physics_buffer, only : pbuf_add_field, dtype_r8

  !-----------------------------------------------------------------------
    if (conv_water_register_codon(1_c_int64_t) == 0_c_int64_t) return

    ! these calls were already done in convect_shallow...so here I add the same fields to the physics buffer with a "1" at the end
! shallow gbm cloud liquid water (kg/kg)
    call pbuf_add_field('SH_CLDLIQ1','physpkg',dtype_r8,(/pcols,pver/),sh_cldliq1_idx)
! shallow gbm cloud ice water (kg/kg)
    call pbuf_add_field('SH_CLDICE1','physpkg',dtype_r8,(/pcols,pver/),sh_cldice1_idx)

  end subroutine conv_water_register


  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine conv_water_init()
   ! --------------------------------------------------------------------- !
   ! Purpose:                                                              !
   !   Initializes the pbuf indices required by conv_water                 !
   ! --------------------------------------------------------------------- !


   use physics_buffer, only : pbuf_get_index
   use cam_history,    only : phys_decomp, addfld

   use constituents,  only: cnst_get_ind

   implicit none
   if (conv_water_init_codon(1_c_int64_t) == 0_c_int64_t) return
   call conv_water_log_direct(conv_water_init_logged, &
        'conv_water_init direct = codon; pbuf/cnst/history native CAM API islands')

   call cnst_get_ind('CLDICE', ixcldice)
   call cnst_get_ind('CLDLIQ', ixcldliq)

   icwmrsh_idx  = pbuf_get_index('ICWMRSH')
   icwmrdp_idx  = pbuf_get_index('ICWMRDP')
   fice_idx     = pbuf_get_index('FICE')
   sh_frac_idx  = pbuf_get_index('SH_FRAC')
   dp_frac_idx  = pbuf_get_index('DP_FRAC')
   ast_idx      = pbuf_get_index('AST')
   rei_idx      = pbuf_get_index('REI')

   ! Convective cloud water variables.
   call addfld ('ICIMRCU  ', 'kg/kg   ', pver, 'A', 'Convection in-cloud ice mixing ratio '   , phys_decomp)
   call addfld ('ICLMRCU  ', 'kg/kg   ', pver, 'A', 'Convection in-cloud liquid mixing ratio ', phys_decomp)
   call addfld ('ICIMRTOT ', 'kg/kg   ', pver, 'A', 'Total in-cloud ice mixing ratio '        , phys_decomp)
   call addfld ('ICLMRTOT ', 'kg/kg   ', pver, 'A', 'Total in-cloud liquid mixing ratio '     , phys_decomp)

   call addfld ('GCLMRDP  ', 'kg/kg   ', pver, 'A', 'Grid-mean deep convective LWC'           , phys_decomp)
   call addfld ('GCIMRDP  ', 'kg/kg   ', pver, 'A', 'Grid-mean deep convective IWC'           , phys_decomp)
   call addfld ('GCLMRSH  ', 'kg/kg   ', pver, 'A', 'Grid-mean shallow convective LWC'        , phys_decomp)
   call addfld ('GCIMRSH  ', 'kg/kg   ', pver, 'A', 'Grid-mean shallow convective IWC'        , phys_decomp)
   call addfld ('FRESH  ', '1', pver, 'A', 'Fractional occurrence of shallow cumulus with condensate', phys_decomp)
   call addfld ('FREDP  ', '1', pver, 'A', 'Fractional occurrence of deep cumulus with condensate', phys_decomp)
   call addfld ('FRECU  ', '1', pver, 'A', 'Fractional occurrence of cumulus with condensate', phys_decomp)
   call addfld ('FRETOT ', '1', pver, 'A', 'Fractional occurrence of cloud with condensate', phys_decomp)

   end subroutine conv_water_init

   subroutine conv_water_4rad(state, pbuf, totg_liq, totg_ice)

   ! --------------------------------------------------------------------- !
   ! Purpose:                                                              !
   ! Computes grid-box average liquid (and ice) from stratus and cumulus   !
   ! Just for the purposes of radiation.                                   !
   !                                                                       !
   ! Method:                                                               !
   ! Extract information about deep+shallow liquid and cloud fraction from !
   ! the physics buffer.                                                   !
   !                                                                       !
   !---------------------------------------------------------------------- !

   use iso_c_binding,  only: c_double, c_int64_t, c_loc, c_ptr
   use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

   use physics_types,   only: physics_state
   use cam_history,     only: outfld
   use phys_control,    only: phys_getopts

   implicit none

   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)

   real(r8), target, intent(out) :: totg_liq(pcols,pver)
   real(r8), target, intent(out) :: totg_ice(pcols,pver)

   real(r8), pointer, dimension(:,:) :: pdel
   real(r8), pointer, dimension(:,:) :: ls_liq
   real(r8), pointer, dimension(:,:) :: ls_ice
   real(r8), pointer, dimension(:,:) :: ast
   real(r8), pointer, dimension(:,:) :: sh_frac
   real(r8), pointer, dimension(:,:) :: dp_frac
   real(r8), pointer, dimension(:,:) :: rei
   real(r8), pointer, dimension(:,:) :: dp_icwmr
   real(r8), pointer, dimension(:,:) :: sh_icwmr
   real(r8), pointer, dimension(:,:) :: fice
   real(r8), pointer, dimension(:,:) :: sh_cldliq
   real(r8), pointer, dimension(:,:) :: sh_cldice

   real(r8), target :: conv_ice(pcols,pver)
   real(r8), target :: conv_liq(pcols,pver)
   real(r8), target :: tot_ice(pcols,pver)
   real(r8), target :: tot_liq(pcols,pver)
   real(r8), target :: totg_ice_sh(pcols,pver)
   real(r8), target :: totg_liq_sh(pcols,pver)
   real(r8), target :: totg_ice_dp(pcols,pver)
   real(r8), target :: totg_liq_dp(pcols,pver)
   real(r8), target :: fresh(pcols,pver)
   real(r8), target :: fredp(pcols,pver)
   real(r8), target :: frecu(pcols,pver)
   real(r8), target :: fretot(pcols,pver)

   integer :: itim_old
   integer :: lchnk
   integer :: ncol
   integer :: microp_is_rk
   character(len=16) :: microp_scheme

   call conv_water_select_impl()

   if (use_native_impl) then
      call conv_water_4rad_native(state, pbuf, totg_liq, totg_ice)
      return
   end if

   ncol  = state%ncol
   lchnk = state%lchnk
   pdel   => state%pdel
   ls_liq => state%q(:,:,ixcldliq)
   ls_ice => state%q(:,:,ixcldice)

   call phys_getopts(microp_scheme_out = microp_scheme)
   microp_is_rk = 0
   if (microp_scheme == 'RK') microp_is_rk = 1

   call pbuf_get_field(pbuf, icwmrsh_idx, sh_icwmr)
   call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr)
   call pbuf_get_field(pbuf, fice_idx,    fice)

   call pbuf_get_field(pbuf, sh_frac_idx, sh_frac)
   call pbuf_get_field(pbuf, dp_frac_idx, dp_frac)
   call pbuf_get_field(pbuf, rei_idx,     rei)

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, sh_cldliq1_idx, sh_cldliq)
   call pbuf_get_field(pbuf, sh_cldice1_idx, sh_cldice)

   call conv_water_4rad_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(conv_water_mode, c_int64_t), &
        int(microp_is_rk, c_int64_t), real(frac_limit, c_double), real(gravit, c_double), &
        c_loc(pdel), c_loc(ls_liq), c_loc(ls_ice), c_loc(ast), c_loc(sh_frac), c_loc(dp_frac), c_loc(rei), &
        c_loc(dp_icwmr), c_loc(sh_icwmr), c_loc(fice), c_loc(totg_liq), c_loc(totg_ice), c_loc(conv_ice), &
        c_loc(conv_liq), c_loc(tot_ice), c_loc(tot_liq), c_loc(totg_ice_sh), c_loc(totg_liq_sh), c_loc(totg_ice_dp), &
        c_loc(totg_liq_dp), c_loc(fresh), c_loc(fredp), c_loc(frecu), c_loc(fretot), c_loc(sh_cldliq), c_loc(sh_cldice) &
   )
   call conv_water_4rad_log_direct()

   call outfld( 'ICLMRCU ', conv_liq  , pcols, lchnk )
   call outfld( 'ICIMRCU ', conv_ice  , pcols, lchnk )
   call outfld( 'ICLMRTOT', tot_liq   , pcols, lchnk )
   call outfld( 'ICIMRTOT', tot_ice   , pcols, lchnk )

   call outfld('GCLMRDP', totg_liq_dp, pcols, lchnk)
   call outfld('GCIMRDP', totg_ice_dp, pcols, lchnk)
   call outfld('GCLMRSH', totg_liq_sh, pcols, lchnk)
   call outfld('GCIMRSH', totg_ice_sh, pcols, lchnk)
   call outfld('FRESH',   fresh,       pcols, lchnk)
   call outfld('FREDP',   fredp,       pcols, lchnk)
   call outfld('FRECU',   frecu,       pcols, lchnk)
   call outfld('FRETOT',  fretot,      pcols, lchnk)

   end subroutine conv_water_4rad

   subroutine conv_water_4rad_native(state, pbuf, totg_liq, totg_ice)

   use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx

   use physics_types,   only: physics_state
   use cam_history,     only: outfld
   use phys_control,    only: phys_getopts

   implicit none

   type(physics_state), target, intent(in) :: state        ! state variables
   type(physics_buffer_desc),   pointer    :: pbuf(:)

   real(r8), intent(out):: totg_liq(pcols,pver)   ! Total GBA in-cloud liquid
   real(r8), intent(out):: totg_ice(pcols,pver)   ! Total GBA in-cloud ice

   real(r8), pointer, dimension(:,:) ::  pdel     ! Moist pressure difference across layer
   real(r8), pointer, dimension(:,:) ::  ls_liq   ! Large-scale contributions to GBA cloud liq
   real(r8), pointer, dimension(:,:) ::  ls_ice   ! Large-scale contributions to GBA cloud ice

   real(r8), pointer, dimension(:,:) ::  ast      ! Physical liquid+ice stratus cloud fraction
   real(r8), pointer, dimension(:,:) ::  sh_frac  ! Shallow convective cloud fraction
   real(r8), pointer, dimension(:,:) ::  dp_frac  ! Deep convective cloud fraction
   real(r8), pointer, dimension(:,:) ::  rei      ! Ice effective drop size (microns)

   real(r8), pointer, dimension(:,:) ::  dp_icwmr ! Deep conv. cloud water
   real(r8), pointer, dimension(:,:) ::  sh_icwmr ! Shallow conv. cloud water
   real(r8), pointer, dimension(:,:) ::  fice     ! Ice partitioning ratio
   real(r8), pointer, dimension(:,:) ::  sh_cldliq ! shallow convection gbx liq cld mixing ratio for COSP
   real(r8), pointer, dimension(:,:) ::  sh_cldice ! shallow convection gbx ice cld mixing ratio for COSP

   real(r8) :: conv_ice(pcols,pver)               ! Convective contributions to IC cloud ice
   real(r8) :: conv_liq(pcols,pver)               ! Convective contributions to IC cloud liquid
   real(r8) :: tot_ice(pcols,pver)                ! Total IC ice
   real(r8) :: tot_liq(pcols,pver)                ! Total IC liquid

   integer  :: i,k,itim_old                       ! Lon, lev indices buff stuff.
   real(r8) :: cu_icwmr                           ! Convective  water for this grid-box.
   real(r8) :: ls_icwmr                           ! Large-scale water for this grid-box.
   real(r8) :: tot_icwmr                          ! Large-scale water for this grid-box.
   real(r8) :: ls_frac                            ! Large-scale cloud frac for this grid-box.
   real(r8) :: tot0_frac, cu0_frac, dp0_frac, sh0_frac
   real(r8) :: kabs, kabsi, kabsl, alpha, dp0, sh0, ic_limit
   real(r8) :: wrk1

   real(r8) :: totg_ice_sh(pcols,pver)   ! Grid-mean IWP from shallow convective cloud
   real(r8) :: totg_liq_sh(pcols,pver)   ! Grid-mean LWP from shallow convective cloud
   real(r8) :: totg_ice_dp(pcols,pver)   ! Grid-mean IWP from deep convective cloud
   real(r8) :: totg_liq_dp(pcols,pver)   ! Grid-mean LWP from deep convective cloud
   real(r8) :: fresh(pcols,pver)         ! Fractional occurrence of shallow cumulus
   real(r8) :: fredp(pcols,pver)         ! Fractional occurrence of deep cumulus
   real(r8) :: frecu(pcols,pver)         ! Fractional occurrence of cumulus
   real(r8) :: fretot(pcols,pver)        ! Fractional occurrence of cloud

   integer :: lchnk
   integer :: ncol

   parameter( kabsl = 0.090361_r8, ic_limit = 1.e-12_r8 )
   character(len=16) :: microp_scheme

   ncol  = state%ncol
   lchnk = state%lchnk
   pdel   => state%pdel
   ls_liq => state%q(:,:,ixcldliq)
   ls_ice => state%q(:,:,ixcldice)

   call phys_getopts( microp_scheme_out = microp_scheme )

   call pbuf_get_field(pbuf, icwmrsh_idx, sh_icwmr )
   call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr )
   call pbuf_get_field(pbuf, fice_idx,    fice )

   call pbuf_get_field(pbuf, sh_frac_idx,  sh_frac )
   call pbuf_get_field(pbuf, dp_frac_idx,  dp_frac )
   call pbuf_get_field(pbuf, rei_idx,      rei )

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,  ast,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   fresh(:,:)  = 0._r8
   fredp(:,:)  = 0._r8
   frecu(:,:)  = 0._r8
   fretot(:,:) = 0._r8

   do k = 1, pver
   do i = 1, ncol

      if( sh_frac(i,k) <= frac_limit .or. sh_icwmr(i,k) <= ic_limit ) then
          sh0_frac = 0._r8
      else
          sh0_frac = sh_frac(i,k)
      endif
      if( dp_frac(i,k) <= frac_limit .or. dp_icwmr(i,k) <= ic_limit ) then
          dp0_frac = 0._r8
      else
          dp0_frac = dp_frac(i,k)
      endif
      cu0_frac = sh0_frac + dp0_frac

      wrk1 = min(1._r8,max(0._r8, ls_ice(i,k)/(ls_ice(i,k)+ls_liq(i,k)+1.e-36_r8)))

      if( ( cu0_frac < frac_limit ) .or. ( ( sh_icwmr(i,k) + dp_icwmr(i,k) ) < ic_limit ) ) then

            cu0_frac = 0._r8
            cu_icwmr = 0._r8

            ls_frac = ast(i,k)
            if( ls_frac < frac_limit ) then
                ls_frac  = 0._r8
                ls_icwmr = 0._r8
            else
                ls_icwmr = ( ls_liq(i,k) + ls_ice(i,k) )/max(frac_limit,ls_frac)
            end if

            tot0_frac = ls_frac
            tot_icwmr = ls_icwmr

      else

            if( microp_scheme == 'RK' ) then
               kabsi = 0.005_r8 + 1._r8/rei(i,k)
            else
               kabsi = 0.005_r8 + 1._r8/min(max(13._r8,rei(i,k)),130._r8)
            endif
            kabs  = kabsl * ( 1._r8 - wrk1 ) + kabsi * wrk1
            alpha = -1.66_r8*kabs*pdel(i,k)/gravit*1000.0_r8

            select case (conv_water_mode)
            case (1)
               cu_icwmr = ( sh0_frac * sh_icwmr(i,k) + dp0_frac*dp_icwmr(i,k))/max(frac_limit,cu0_frac)
            case (2)
               sh0 = exp(alpha*sh_icwmr(i,k))
               dp0 = exp(alpha*dp_icwmr(i,k))
               cu_icwmr = log((sh0_frac*sh0+dp0_frac*dp0)/max(frac_limit,cu0_frac))
               cu_icwmr = cu_icwmr/alpha
            case default
            end select

            ls_frac   = ast(i,k)
            ls_icwmr  = (ls_liq(i,k) + ls_ice(i,k))/max(frac_limit,ls_frac)
            tot0_frac = (ls_frac + cu0_frac)

            select case (conv_water_mode)
            case (1)
               tot_icwmr = (ls_frac*ls_icwmr + cu0_frac*cu_icwmr)/max(frac_limit,tot0_frac)
            case (2)
               tot_icwmr = log((ls_frac*exp(alpha*ls_icwmr)+cu0_frac*exp(alpha*cu_icwmr))/max(frac_limit,tot0_frac))
               tot_icwmr = tot_icwmr/alpha
            case default
            end select

      end if

      conv_ice(i,k) = cu_icwmr * wrk1
      conv_liq(i,k) = cu_icwmr * (1._r8-wrk1)

      tot_ice(i,k)  = tot_icwmr * wrk1
      tot_liq(i,k)  = tot_icwmr * (1._r8-wrk1)

      totg_ice(i,k) = tot0_frac * tot_icwmr * wrk1
      totg_liq(i,k) = tot0_frac * tot_icwmr * (1._r8-wrk1)

      totg_ice_sh(i,k)  = sh0_frac * sh_icwmr(i,k) * wrk1
      totg_ice_dp(i,k)  = dp0_frac * dp_icwmr(i,k) * wrk1
      totg_liq_sh(i,k)  = sh0_frac * sh_icwmr(i,k) * (1._r8-wrk1)
      totg_liq_dp(i,k)  = dp0_frac * dp_icwmr(i,k) * (1._r8-wrk1)
      if( sh0_frac > frac_limit ) then
          fresh(i,k) = 1._r8
      endif
      if( dp0_frac > frac_limit ) then
          fredp(i,k) = 1._r8
      endif
      if( cu0_frac > frac_limit ) then
          frecu(i,k) = 1._r8
      endif
      if( tot0_frac > frac_limit ) then
          fretot(i,k) = 1._r8
      endif

   end do
   end do

   call pbuf_get_field(pbuf, sh_cldliq1_idx, sh_cldliq  )
   call pbuf_get_field(pbuf, sh_cldice1_idx, sh_cldice  )

   sh_cldliq(:ncol,:pver)=sh_icwmr(:ncol,:pver)*(1-fice(:ncol,:pver))*sh_frac(:ncol,:pver)
   sh_cldice(:ncol,:pver)=sh_icwmr(:ncol,:pver)*fice(:ncol,:pver)*sh_frac(:ncol,:pver)

   call outfld( 'ICLMRCU ', conv_liq  , pcols, lchnk )
   call outfld( 'ICIMRCU ', conv_ice  , pcols, lchnk )
   call outfld( 'ICLMRTOT', tot_liq   , pcols, lchnk )
   call outfld( 'ICIMRTOT', tot_ice   , pcols, lchnk )

   call outfld('GCLMRDP', totg_liq_dp, pcols, lchnk)
   call outfld('GCIMRDP', totg_ice_dp, pcols, lchnk)
   call outfld('GCLMRSH', totg_liq_sh, pcols, lchnk)
   call outfld('GCIMRSH', totg_ice_sh, pcols, lchnk)
   call outfld('FRESH',   fresh,       pcols, lchnk)
   call outfld('FREDP',   fredp,       pcols, lchnk)
   call outfld('FRECU',   frecu,       pcols, lchnk)
   call outfld('FRETOT',  fretot,      pcols, lchnk)

   end subroutine conv_water_4rad_native

end module conv_water
