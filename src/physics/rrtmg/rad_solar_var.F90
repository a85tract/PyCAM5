!-------------------------------------------------------------------------------
! This module uses the Lean solar irradiance data to provide a solar cycle
! scaling factor used in heating rate calculations 
!-------------------------------------------------------------------------------
module rad_solar_var

  use shr_kind_mod ,     only : r8 => shr_kind_r8
  use solar_data,        only : sol_irrad, we, nbins, has_spectrum, sol_tsi
  use solar_data,        only : do_spctrl_scaling
  use cam_abortutils,    only : endrun
  use cam_logfile,       only : iulog
  use spmd_utils,        only : masterproc
  use iso_c_binding,     only : c_int64_t, c_loc, c_ptr

  implicit none
  save

  private
  public :: rad_solar_var_init
  public :: get_variability

  real(r8), allocatable :: ref_band_irrad(:)  ! scaling will be relative to ref_band_irrad in each band
  real(r8), allocatable :: irrad(:)           ! solar irradiance at model timestep in each band
  real(r8)              :: tsi_ref            ! total solar irradiance assumed by rrtmg                                                 

  real(r8), allocatable, target :: radbinmax(:)
  real(r8), allocatable :: radbinmin(:)
  integer :: nradbins
  logical :: use_native_rad_solar_var_init_impl = .false.
  logical :: rad_solar_var_init_impl_selected = .false.
  logical :: rad_solar_var_init_logged = .false.
  logical :: use_native_rrtmg_solar_variability_impl = .false.
  logical :: rrtmg_solar_variability_impl_selected = .false.
  logical :: rrtmg_solar_variability_entered_logged = .false.
  interface
     function rad_solar_var_init_far_ir_codon(nradbins_c, radbinmax_p) result(radmax_loc_c) &
          bind(c, name="rad_solar_var_init_far_ir_codon")
       import :: c_int64_t, c_ptr
       integer(c_int64_t), value :: nradbins_c
       type(c_ptr), value :: radbinmax_p
       integer(c_int64_t) :: radmax_loc_c
     end function rad_solar_var_init_far_ir_codon
  end interface
contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine rad_solar_var_init( )
    use radconstants,  only : get_number_sw_bands
    use radconstants,  only : get_sw_spectral_boundaries
    use radconstants,  only : get_ref_solar_band_irrad
    use radconstants,  only : get_ref_total_solar_irrad

    integer :: i
    integer :: ierr
    integer :: yr, mon, tod
    integer :: radmax_loc


    call get_number_sw_bands(nradbins)

    if ( do_spctrl_scaling ) then

       if ( .not.has_spectrum ) then
          call endrun('rad_solar_var_init: solar input file must have irradiance spectrum')
       endif

       allocate (radbinmax(nradbins),stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for radbinmax')
       end if

       allocate (radbinmin(nradbins),stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for radbinmin')
       end if

       allocate (ref_band_irrad(nradbins), stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for ref_band_irrad')
       end if

       allocate (irrad(nradbins), stat=ierr)
       if (ierr /= 0) then
          call endrun('rad_solar_var_init: Error allocating space for irrad')
       end if

       call get_sw_spectral_boundaries(radbinmin, radbinmax, 'nm')

       ! Make sure that the far-IR is included, even if RRTMG does not
       ! extend that far down. 10^5 nm corresponds to a wavenumber of
       ! 100 cm^-1.
      call rad_solar_var_init_select_impl()
      if (use_native_rad_solar_var_init_impl) then
         radmax_loc = maxloc(radbinmax,1)
         radbinmax(radmax_loc) = max(100000._r8,radbinmax(radmax_loc))
      else
         call rad_solar_var_init_log()
         radmax_loc = int(rad_solar_var_init_far_ir_codon(int(nradbins, c_int64_t), c_loc(radbinmax(1))))
      end if

       ! for rrtmg, reference spectrum from rrtmg
       call get_ref_solar_band_irrad( ref_band_irrad )

    else

       call get_ref_total_solar_irrad(tsi_ref)

    endif

  endsubroutine rad_solar_var_init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine rad_solar_var_init_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (rad_solar_var_init_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('RRTMG_RAD_SOLAR_VAR_INIT_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_rad_solar_var_init_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_rad_solar_var_init_impl = .false.
    end if

    rad_solar_var_init_impl_selected = .true.

  end subroutine rad_solar_var_init_select_impl

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine rad_solar_var_init_log()

    if (rad_solar_var_init_logged) return
    rad_solar_var_init_logged = .true.

    if (masterproc) then
       write(iulog,*) 'rad_solar_var_init implementation = codon'
       call flush(iulog)
    end if

  end subroutine rad_solar_var_init_log

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine get_variability( sfac )

    real(r8), intent(out) :: sfac(nradbins)       ! scaling factors for CAM heating

    integer :: yr, mon, day, tod

    if ( do_spctrl_scaling ) then

      call rrtmg_solar_variability_select_impl()
      if (use_native_rrtmg_solar_variability_impl) then
         call integrate_spectrum( nbins, nradbins, we, radbinmin, radbinmax, sol_irrad, irrad)
         sfac(:nradbins) = irrad(:nradbins)/ref_band_irrad(:nradbins)
      else
         call rrtmg_solar_variability_log_entered()
         call rrtmg_solar_variability_codon_wrap(nbins, nradbins, we, radbinmin, radbinmax, &
              sol_irrad, irrad, ref_band_irrad, sfac)
      end if

    else

       sfac(:nradbins) = sol_tsi/tsi_ref

    endif

  endsubroutine get_variability

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine rrtmg_solar_variability_codon_wrap(nsrc, ntrg, src_x, min_trg, max_trg, src, trg, ref_irrad, sfac)

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    implicit none

    integer,  intent(in)  :: nsrc
    integer,  intent(in)  :: ntrg
    real(r8), target, intent(in)  :: src_x(nsrc+1)
    real(r8), target, intent(in)  :: max_trg(ntrg)
    real(r8), target, intent(in)  :: min_trg(ntrg)
    real(r8), target, intent(in)  :: src(nsrc)
    real(r8), target, intent(out) :: trg(ntrg)
    real(r8), target, intent(in)  :: ref_irrad(ntrg)
    real(r8), target, intent(out) :: sfac(ntrg)

    interface
       subroutine rrtmg_solar_variability_codon(nsrc_c, ntrg_c, src_x_p, min_trg_p, max_trg_p, &
            src_p, trg_p, ref_irrad_p, sfac_p) bind(c, name="rrtmg_solar_variability_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: nsrc_c, ntrg_c
         type(c_ptr), value :: src_x_p, min_trg_p, max_trg_p, src_p, trg_p, ref_irrad_p, sfac_p
       end subroutine rrtmg_solar_variability_codon
    end interface

    call rrtmg_solar_variability_codon(int(nsrc, c_int64_t), int(ntrg, c_int64_t), &
         c_loc(src_x(1)), c_loc(min_trg(1)), c_loc(max_trg(1)), c_loc(src(1)), &
         c_loc(trg(1)), c_loc(ref_irrad(1)), c_loc(sfac(1)))

  end subroutine rrtmg_solar_variability_codon_wrap

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine rrtmg_solar_variability_select_impl()

    implicit none

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (rrtmg_solar_variability_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('RRTMG_SOLAR_VARIABILITY_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_rrtmg_solar_variability_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_rrtmg_solar_variability_impl = .false.
    end if

    rrtmg_solar_variability_impl_selected = .true.

    if (masterproc) then
       if (use_native_rrtmg_solar_variability_impl) then
          write(iulog,*) 'rrtmg_solar_variability implementation = native'
       else
          write(iulog,*) 'rrtmg_solar_variability implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine rrtmg_solar_variability_select_impl

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine rrtmg_solar_variability_log_entered()

    implicit none

    if (rrtmg_solar_variability_entered_logged) return
    rrtmg_solar_variability_entered_logged = .true.

    if (masterproc) then
       write(iulog,*) 'rrtmg_solar_variability entered (solar spectral integrate/scale = codon)'
       call flush(iulog)
    end if

  end subroutine rrtmg_solar_variability_log_entered

!-------------------------------------------------------------------------------
! private method.........
!-------------------------------------------------------------------------------

  subroutine integrate_spectrum( nsrc, ntrg, src_x, min_trg, max_trg, src, trg )

    use mo_util, only : rebin

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer,  intent(in)  :: nsrc                  ! dimension source array
    integer,  intent(in)  :: ntrg                  ! dimension target array
    real(r8), intent(in)  :: src_x(nsrc+1)         ! source coordinates
    real(r8), intent(in)  :: max_trg(ntrg)         ! target coordinates
    real(r8), intent(in)  :: min_trg(ntrg)         ! target coordinates
    real(r8), intent(in)  :: src(nsrc)             ! source array
    real(r8), intent(out) :: trg(ntrg)             ! target array
 
    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    real(r8) :: trg_x(2), targ(1)         ! target coordinates
    integer  :: i

    do i = 1, ntrg

       trg_x(1) = min_trg(i)
       trg_x(2) = max_trg(i)

       call rebin( nsrc, 1, src_x, trg_x, src(1:nsrc), targ(:) )
       ! W/m2/nm --> W/m2
       trg( i ) = targ(1)*(trg_x(2)-trg_x(1))

    enddo


  end subroutine integrate_spectrum

endmodule rad_solar_var
