module tidal_diag 

  !---------------------------------------------------------------------------------
  ! Module to compute fourier coefficients for the diurnal and semidiurnal tide 
  !
  ! Created by: Dan Marsh
  ! Date: 12 May 2008
  !---------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use ppgrid,        only: pcols, pver
  use spmd_utils,    only: masterproc
  use cam_logfile,   only: iulog
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

  implicit none

  private

  ! Public interfaces

  public :: tidal_diag_init   ! create coefficient history file variables
  public :: tidal_diag_write  ! calculate and output dignostics
  public :: get_tidal_coeffs

  logical :: use_native_tidal_diag_impl = .false.
  logical :: tidal_diag_impl_selected = .false.
  logical :: tidal_diag_proof_written = .false.
  logical :: tidal_diag_init_logged = .false.
  logical :: tidal_diag_write_logged = .false.
  logical :: tidal_diag_write_scale_logged = .false.
  logical :: get_tidal_coeffs_logged = .false.

  interface
    function tidal_diag_int_codon(value_c, force_one_c) result(out_c) bind(c, name="tidal_diag_int_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: value_c, force_one_c
      integer(c_int64_t) :: out_c
    end function tidal_diag_int_codon
    subroutine get_tidal_coeffs_codon(tod_c, pi_c, cday_c, dcoef_p) bind(c, name="get_tidal_coeffs_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: tod_c
      real(c_double), value :: pi_c, cday_c
      type(c_ptr), value :: dcoef_p
    end subroutine get_tidal_coeffs_codon
    subroutine tidal_diag_scale_2d_codon(ncol_c, pcols_c, pver_c, coef_c, src_p, dst_p) &
         bind(c, name="tidal_diag_scale_2d_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
      real(c_double), value :: coef_c
      type(c_ptr), value :: src_p, dst_p
    end subroutine tidal_diag_scale_2d_codon
    subroutine tidal_diag_scale_1d_codon(ncol_c, coef_c, src_p, dst_p) &
         bind(c, name="tidal_diag_scale_1d_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c
      real(c_double), value :: coef_c
      type(c_ptr), value :: src_p, dst_p
    end subroutine tidal_diag_scale_1d_codon
  end interface

contains

  !===============================================================================

  subroutine tidal_diag_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (tidal_diag_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('TIDAL_DIAG_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_tidal_diag_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_tidal_diag_impl = .false.
    end if

    tidal_diag_impl_selected = .true.

    if (masterproc) then
       if (use_native_tidal_diag_impl) then
          write(iulog,*) 'tidal_diag implementation = native'
       else
          write(iulog,*) 'tidal_diag implementation = codon'
       end if
    end if

  end subroutine tidal_diag_select_impl

  !===============================================================================

  subroutine tidal_diag_proof_once()

    if (tidal_diag_proof_written) return
    tidal_diag_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'tidal_diag entered (history scalar helpers = codon)'
    end if

  end subroutine tidal_diag_proof_once

  !===============================================================================

  integer function tidal_diag_int(value_in, force_one)

    integer, intent(in) :: value_in
    logical, intent(in) :: force_one
    integer(c_int64_t) :: force_one_c, out_c

    call tidal_diag_select_impl()

    if (use_native_tidal_diag_impl) then
       if (force_one) then
          tidal_diag_int = 1
       else
          tidal_diag_int = value_in
       end if
       return
    end if

    call tidal_diag_proof_once()
    if (force_one) then
       force_one_c = 1_c_int64_t
    else
       force_one_c = 0_c_int64_t
    end if
    out_c = tidal_diag_int_codon(int(value_in, c_int64_t), force_one_c)
    tidal_diag_int = int(out_c)

  end function tidal_diag_int

  !===============================================================================

  subroutine tidal_diag_log_direct(logged, proof_line)

    logical, intent(inout) :: logged
    character(len=*), intent(in) :: proof_line

    if (logged) return
    logged = .true.

    if (masterproc) then
       write(iulog,'(A)') trim(proof_line)
    end if

  end subroutine tidal_diag_log_direct

  !===============================================================================

  subroutine  tidal_diag_init()
    !----------------------------------------------------------------------- 
    ! Purpose: create fourier coefficient history file variables
    !-----------------------------------------------------------------------

    use cam_history,        only: addfld, phys_decomp

    integer :: lev3d, levs
    integer(c_int64_t) :: out_c

    call tidal_diag_select_impl()
    if (use_native_tidal_diag_impl) then
       lev3d = pver
       levs = 1
    else
       call tidal_diag_proof_once()
       out_c = tidal_diag_int_codon(int(pver, c_int64_t), 0_c_int64_t)
       lev3d = int(out_c)
       out_c = tidal_diag_int_codon(int(pver, c_int64_t), 1_c_int64_t)
       levs = int(out_c)
       call tidal_diag_log_direct(tidal_diag_init_logged, 'tidal_diag_init direct = codon')
    end if

    call addfld ('T_24_COS','K       ',lev3d, 'A','Temperature 24hr. cos coeff.',phys_decomp)
    call addfld ('T_24_SIN','K       ',lev3d, 'A','Temperature 24hr. sin coeff.',phys_decomp)
    call addfld ('T_12_COS','K       ',lev3d, 'A','Temperature 12hr. cos coeff.',phys_decomp)
    call addfld ('T_12_SIN','K       ',lev3d, 'A','Temperature 12hr. sin coeff.',phys_decomp)

    call addfld ('U_24_COS','m/s     ',lev3d, 'A','Zonal wind 24hr. cos coeff.',phys_decomp)
    call addfld ('U_24_SIN','m/s     ',lev3d, 'A','Zonal wind 24hr. sin coeff.',phys_decomp)
    call addfld ('U_12_COS','m/s     ',lev3d, 'A','Zonal wind 12hr. cos coeff.',phys_decomp)
    call addfld ('U_12_SIN','m/s     ',lev3d, 'A','Zonal wind 12hr. sin coeff.',phys_decomp)

    call addfld ('V_24_COS','m/s     ',lev3d, 'A','Meridional wind 24hr. cos coeff.',phys_decomp)
    call addfld ('V_24_SIN','m/s     ',lev3d, 'A','Meridional wind 24hr. sin coeff.',phys_decomp)
    call addfld ('V_12_COS','m/s     ',lev3d, 'A','Meridional wind 12hr. cos coeff.',phys_decomp)
    call addfld ('V_12_SIN','m/s     ',lev3d, 'A','Meridional wind 12hr. sin coeff.',phys_decomp)

    call addfld ('PS_24_COS','Pa     ',levs,  'A','surface pressure 24hr. cos coeff.',phys_decomp)
    call addfld ('PS_24_SIN','Pa     ',levs,  'A','surface pressure 24hr. sin coeff.',phys_decomp)
    call addfld ('PS_12_COS','Pa     ',levs,  'A','surface pressure 12hr. cos coeff.',phys_decomp)
    call addfld ('PS_12_SIN','Pa     ',levs,  'A','surface pressure 12hr. sin coeff.',phys_decomp)

    call addfld ('OMEGA_24_COS','Pa/s',lev3d, 'A','vertical pressure velocity 24hr. cos coeff.',phys_decomp)
    call addfld ('OMEGA_24_SIN','Pa/s',lev3d, 'A','vertical pressure velocity 24hr. sin coeff.',phys_decomp)
    call addfld ('OMEGA_12_COS','Pa/s',lev3d, 'A','vertical pressure velocity 12hr. cos coeff.',phys_decomp)
    call addfld ('OMEGA_12_SIN','Pa/s',lev3d, 'A','vertical pressure velocity 12hr. sin coeff.',phys_decomp)

    return

  end subroutine tidal_diag_init

  !===============================================================================

  subroutine  tidal_diag_write(state)

    !----------------------------------------------------------------------- 
    ! Purpose: calculate fourier coefficients and save to history files 
    !-----------------------------------------------------------------------
    use cam_history,   only: outfld, hist_fld_active
    use physics_types, only: physics_state

    implicit none

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(physics_state), intent(in) :: state
    !
    !---------------------------Local workspace-----------------------------

    integer  :: lchnk
    integer(c_int64_t) :: out_c

    real(r8) :: dcoef(4) 
    real(r8), target :: field2d(pcols,pver)
    real(r8), target :: field1d(pcols)
    integer :: ncol

    !-----------------------------------------------------------------------

    call tidal_diag_select_impl()
    if (use_native_tidal_diag_impl) then
       lchnk = state%lchnk
       ncol = state%ncol
    else
       call tidal_diag_proof_once()
       out_c = tidal_diag_int_codon(int(state%lchnk, c_int64_t), 0_c_int64_t)
       lchnk = int(out_c)
       out_c = tidal_diag_int_codon(int(state%ncol, c_int64_t), 0_c_int64_t)
       ncol = int(out_c)
       call tidal_diag_log_direct(tidal_diag_write_logged, 'tidal_diag_write direct = codon')
    end if

    call get_tidal_coeffs( dcoef )
    if (.not. use_native_tidal_diag_impl) then
       ! Fixed-case history settings can leave all tidal output fields inactive.
       call tidal_diag_scale_2d(state%t, dcoef(1), ncol, field2d)
       call tidal_diag_scale_1d(state%ps, dcoef(1), ncol, field1d)
    end if

    if ( hist_fld_active('T_24_COS') .or. hist_fld_active('T_24_SIN') ) then
       call tidal_diag_scale_2d(state%t, dcoef(1), ncol, field2d)
       call outfld( 'T_24_SIN', field2d(:ncol,:), ncol, lchnk )
       call tidal_diag_scale_2d(state%t, dcoef(2), ncol, field2d)
       call outfld( 'T_24_COS', field2d(:ncol,:), ncol, lchnk )
    endif
    if ( hist_fld_active('T_12_COS') .or. hist_fld_active('T_12_SIN') ) then
       call tidal_diag_scale_2d(state%t, dcoef(3), ncol, field2d)
       call outfld( 'T_12_SIN', field2d(:ncol,:), ncol, lchnk )
       call tidal_diag_scale_2d(state%t, dcoef(4), ncol, field2d)
       call outfld( 'T_12_COS', field2d(:ncol,:), ncol, lchnk )
    endif

    if ( hist_fld_active('U_24_COS') .or. hist_fld_active('U_24_SIN') ) then
       call tidal_diag_scale_2d(state%u, dcoef(1), ncol, field2d)
       call outfld( 'U_24_SIN', field2d(:ncol,:), ncol, lchnk )
       call tidal_diag_scale_2d(state%u, dcoef(2), ncol, field2d)
       call outfld( 'U_24_COS', field2d(:ncol,:), ncol, lchnk )
    endif
    if ( hist_fld_active('U_12_COS') .or. hist_fld_active('U_12_SIN') ) then
       call tidal_diag_scale_2d(state%u, dcoef(3), ncol, field2d)
       call outfld( 'U_12_SIN', field2d(:ncol,:), ncol, lchnk )
       call tidal_diag_scale_2d(state%u, dcoef(4), ncol, field2d)
       call outfld( 'U_12_COS', field2d(:ncol,:), ncol, lchnk )
    endif

    if ( hist_fld_active('V_24_COS') .or. hist_fld_active('V_24_SIN') ) then
       call tidal_diag_scale_2d(state%v, dcoef(1), ncol, field2d)
       call outfld( 'V_24_SIN', field2d(:ncol,:), ncol, lchnk )
       call tidal_diag_scale_2d(state%v, dcoef(2), ncol, field2d)
       call outfld( 'V_24_COS', field2d(:ncol,:), ncol, lchnk )
    endif
    if ( hist_fld_active('V_12_COS') .or. hist_fld_active('V_12_SIN') ) then
       call tidal_diag_scale_2d(state%v, dcoef(3), ncol, field2d)
       call outfld( 'V_12_SIN', field2d(:ncol,:), ncol, lchnk )
       call tidal_diag_scale_2d(state%v, dcoef(4), ncol, field2d)
       call outfld( 'V_12_COS', field2d(:ncol,:), ncol, lchnk )
    endif

    if ( hist_fld_active('PS_24_COS') .or. hist_fld_active('PS_24_SIN') ) then
       call tidal_diag_scale_1d(state%ps, dcoef(1), ncol, field1d)
       call outfld( 'PS_24_SIN', field1d(:ncol), ncol, lchnk )
       call tidal_diag_scale_1d(state%ps, dcoef(2), ncol, field1d)
       call outfld( 'PS_24_COS', field1d(:ncol), ncol, lchnk )
    endif
    if ( hist_fld_active('PS_12_COS') .or. hist_fld_active('PS_12_SIN') ) then
       call tidal_diag_scale_1d(state%ps, dcoef(3), ncol, field1d)
       call outfld( 'PS_12_SIN', field1d(:ncol), ncol, lchnk )
       call tidal_diag_scale_1d(state%ps, dcoef(4), ncol, field1d)
       call outfld( 'PS_12_COS', field1d(:ncol), ncol, lchnk )
    endif

    if ( hist_fld_active('OMEGA_24_COS') .or. hist_fld_active('OMEGA_24_SIN') ) then
       call tidal_diag_scale_2d(state%omega, dcoef(1), ncol, field2d)
       call outfld( 'OMEGA_24_SIN', field2d(:ncol,:), ncol, lchnk )
       call tidal_diag_scale_2d(state%omega, dcoef(2), ncol, field2d)
       call outfld( 'OMEGA_24_COS', field2d(:ncol,:), ncol, lchnk )
    endif
    if ( hist_fld_active('OMEGA_12_COS') .or. hist_fld_active('OMEGA_12_SIN') ) then
       call tidal_diag_scale_2d(state%omega, dcoef(3), ncol, field2d)
       call outfld( 'OMEGA_12_SIN', field2d(:ncol,:), ncol, lchnk )
       call tidal_diag_scale_2d(state%omega, dcoef(4), ncol, field2d)
       call outfld( 'OMEGA_12_COS', field2d(:ncol,:), ncol, lchnk )
    endif

    return

  end subroutine tidal_diag_write

  !===============================================================================

  subroutine tidal_diag_scale_2d(src, coef, ncol, dst)

    real(r8), target, intent(in)  :: src(pcols,pver)
    real(r8),         intent(in)  :: coef
    integer,          intent(in)  :: ncol
    real(r8), target, intent(out) :: dst(pcols,pver)

    if (use_native_tidal_diag_impl) then
       dst(:ncol,:) = src(:ncol,:) * coef
       return
    end if

    call tidal_diag_log_direct(tidal_diag_write_scale_logged, 'tidal_diag_write field scaling direct = codon')
    call tidal_diag_scale_2d_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         real(coef, c_double), c_loc(src(1,1)), c_loc(dst(1,1)))

  end subroutine tidal_diag_scale_2d

  !===============================================================================

  subroutine tidal_diag_scale_1d(src, coef, ncol, dst)

    real(r8), target, intent(in)  :: src(pcols)
    real(r8),         intent(in)  :: coef
    integer,          intent(in)  :: ncol
    real(r8), target, intent(out) :: dst(pcols)

    if (use_native_tidal_diag_impl) then
       dst(:ncol) = src(:ncol) * coef
       return
    end if

    call tidal_diag_log_direct(tidal_diag_write_scale_logged, 'tidal_diag_write field scaling direct = codon')
    call tidal_diag_scale_1d_codon(int(ncol, c_int64_t), real(coef, c_double), c_loc(src(1)), c_loc(dst(1)))

  end subroutine tidal_diag_scale_1d

  !===============================================================================

  subroutine get_tidal_coeffs( dcoef )

    !----------------------------------------------------------------------- 
    ! Purpose: calculate fourier coefficients
    !-----------------------------------------------------------------------

    use time_manager,  only: get_curr_date               
    use physconst, only: pi, cday

    real(r8), target, intent(out) :: dcoef(4)

 !  variables to calculate tidal coeffs
    real(r8), parameter :: pi_x_2 = 2._r8*pi
    real(r8), parameter :: pi_x_4 = 4._r8*pi
    integer  :: year, month
    integer  :: day              ! day of month
    integer  :: tod              ! time of day (seconds past 0Z) 
    real(r8) :: gmtfrac 

 !  calculate multipliers for Fourier transform in time (tidal analysis)
    call get_curr_date(year, month, day, tod)
    call tidal_diag_select_impl()
    if (use_native_tidal_diag_impl) then
       gmtfrac = tod / cday

       dcoef(1) = 2._r8*sin(pi_x_2*gmtfrac)
       dcoef(2) = 2._r8*cos(pi_x_2*gmtfrac)
       dcoef(3) = 2._r8*sin(pi_x_4*gmtfrac)
       dcoef(4) = 2._r8*cos(pi_x_4*gmtfrac)
       return
    end if

    call tidal_diag_proof_once()
    call get_tidal_coeffs_codon(int(tod, c_int64_t), real(pi, c_double), real(cday, c_double), c_loc(dcoef(1)))
    call tidal_diag_log_direct(get_tidal_coeffs_logged, 'get_tidal_coeffs direct = codon')

  end subroutine get_tidal_coeffs

end module tidal_diag
