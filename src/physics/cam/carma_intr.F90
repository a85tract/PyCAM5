!! This module is stub for a coupler between the CAM model and the Community Aerosol
!! and Radiation Model for Atmospheres (CARMA) microphysics model. It is used when
!! CARMA is not being used, so that the CAM code that calls CARMA does not need to
!! be changed. The real version of this routine exists in the directory
!! physics/carma/cam. A CARMA model can be activated by using configure with the
!! option:
!!
!!  -carma <carma_pkg>
!!
!! where carma_pkg is the name for a particular microphysical model.
!!
!! @author  Chuck Bardeen
!! @version May 2009
module carma_intr

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use pmgrid,         only: plat, plev, plevp, plon
  use ppgrid,         only: pcols, pver, pverp
  use constituents,   only: pcnst
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use iso_c_binding,  only: c_int64_t, c_loc, c_null_ptr, c_ptr


  implicit none
  
  private
  save

  ! Public interfaces
  
  ! CAM Physics Interface
  public carma_register                 ! register consituents
  public carma_is_active                ! retrns true if this package is active (microphysics = .true.)
  public carma_implements_cnst          ! returns true if consituent is implemented by this package
  public carma_init_cnst                ! initialize constituent mixing ratios, if not read from initial file
  public carma_init                     ! initialize timestep independent variables
  public carma_final                    ! finalize the CARMA module
  public carma_timestep_init            ! initialize timestep dependent variables
  public carma_timestep_tend            ! interface to tendency computation
  public carma_accumulate_stats         ! collect stats from all MPI tasks
  
  ! Other Microphysics
  public carma_emission_tend            ! calculate tendency from emission source function
  public carma_wetdep_tend              ! calculate tendency from wet deposition
  
  logical :: use_native_carma_intr_impl = .false.
  logical :: carma_intr_impl_selected = .false.
  logical :: carma_intr_proof_written = .false.
  logical :: carma_register_logged = .false.
  logical :: carma_implements_cnst_logged = .false.
  logical :: carma_init_logged = .false.
  logical :: carma_final_logged = .false.
  logical :: carma_timestep_init_logged = .false.
  logical :: carma_timestep_tend_logged = .false.
  logical :: carma_accumulate_stats_logged = .false.

  interface
    function carma_intr_false_codon() result(out_c) bind(c, name="carma_intr_false_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_intr_false_codon

    function carma_intr_touch_codon() result(out_c) bind(c, name="carma_intr_touch_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_intr_touch_codon

    function carma_register_codon() result(out_c) bind(c, name="carma_register_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_register_codon

    function carma_implements_cnst_codon() result(out_c) bind(c, name="carma_implements_cnst_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_implements_cnst_codon

    function carma_init_codon() result(out_c) bind(c, name="carma_init_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_init_codon

    function carma_final_codon() result(out_c) bind(c, name="carma_final_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_final_codon

    function carma_timestep_init_codon() result(out_c) bind(c, name="carma_timestep_init_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_timestep_init_codon

    subroutine carma_timestep_tend_codon(pcols_c, prec_str_present_c, snow_str_present_c, &
         prec_sed_present_c, snow_sed_present_c, prec_str_p, snow_str_p, prec_sed_p, snow_sed_p) &
         bind(c, name="carma_timestep_tend_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: pcols_c, prec_str_present_c, snow_str_present_c
      integer(c_int64_t), value :: prec_sed_present_c, snow_sed_present_c
      type(c_ptr), value :: prec_str_p, snow_str_p, prec_sed_p, snow_sed_p
    end subroutine carma_timestep_tend_codon

    function carma_accumulate_stats_codon() result(out_c) bind(c, name="carma_accumulate_stats_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_accumulate_stats_codon
  end interface

contains

  subroutine carma_intr_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (carma_intr_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('CARMA_INTR_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_carma_intr_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_carma_intr_impl = .false.
    end if

    carma_intr_impl_selected = .true.

    if (masterproc) then
       if (use_native_carma_intr_impl) then
          write(iulog,*) 'carma_intr implementation = native'
       else
          write(iulog,*) 'carma_intr implementation = codon'
       end if
    end if

  end subroutine carma_intr_select_impl

  !================================================================================================

  subroutine carma_intr_proof_once()

    if (carma_intr_proof_written) return
    carma_intr_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'carma_intr entered (stub false/no-op helpers = codon)'
    end if

  end subroutine carma_intr_proof_once

  !================================================================================================

  subroutine carma_intr_log_direct(logged, proof_line)

    logical, intent(inout) :: logged
    character(len=*), intent(in) :: proof_line

    if (logged) return
    logged = .true.

    if (masterproc) then
       write(iulog,'(A)') trim(proof_line)
    end if

  end subroutine carma_intr_log_direct

  !================================================================================================

  subroutine carma_intr_touch()

    integer(c_int64_t) :: out_c

    call carma_intr_select_impl()

    if (use_native_carma_intr_impl) then
       return
    end if

    call carma_intr_proof_once()
    out_c = carma_intr_touch_codon()

  end subroutine carma_intr_touch

  !================================================================================================

  logical function carma_intr_false()

    integer(c_int64_t) :: out_c

    call carma_intr_select_impl()

    if (use_native_carma_intr_impl) then
       carma_intr_false = .false.
       return
    end if

    call carma_intr_proof_once()
    out_c = carma_intr_false_codon()
    carma_intr_false = out_c /= 0_c_int64_t

  end function carma_intr_false

  !================================================================================================


  subroutine carma_register
    integer(c_int64_t) :: out_c

    call carma_intr_select_impl()

    if (use_native_carma_intr_impl) then
       call carma_register_native()
       return
    end if

    call carma_intr_log_direct(carma_register_logged, 'carma_register direct = codon')
    out_c = carma_register_codon()

    return
  end subroutine carma_register

  subroutine carma_register_native
    implicit none

    return
  end subroutine carma_register_native


  function carma_is_active()
    implicit none
  
    logical :: carma_is_active
  
    carma_is_active = carma_intr_false()
    
    return
  end function carma_is_active


  function carma_implements_cnst(name)
    implicit none
    
    character(len=*), intent(in) :: name   !! constituent name
    logical :: carma_implements_cnst       ! return value
    integer(c_int64_t) :: out_c

    call carma_intr_select_impl()

    if (use_native_carma_intr_impl) then
       carma_implements_cnst = carma_implements_cnst_native(name)
       return
    end if

    call carma_intr_log_direct(carma_implements_cnst_logged, 'carma_implements_cnst direct = codon')
    out_c = carma_implements_cnst_codon()
    carma_implements_cnst = out_c /= 0_c_int64_t

    return
  end function carma_implements_cnst

  function carma_implements_cnst_native(name)
    implicit none

    character(len=*), intent(in) :: name
    logical :: carma_implements_cnst_native

    carma_implements_cnst_native = .false.

    return
  end function carma_implements_cnst_native
  

  subroutine carma_init
    integer(c_int64_t) :: out_c

    call carma_intr_select_impl()

    if (use_native_carma_intr_impl) then
       call carma_init_native()
       return
    end if

    call carma_intr_log_direct(carma_init_logged, 'carma_init direct = codon')
    out_c = carma_init_codon()
    
    return
  end subroutine carma_init

  subroutine carma_init_native
    implicit none

    return
  end subroutine carma_init_native


  subroutine carma_final
    integer(c_int64_t) :: out_c

    call carma_intr_select_impl()

    if (use_native_carma_intr_impl) then
       call carma_final_native()
       return
    end if

    call carma_intr_log_direct(carma_final_logged, 'carma_final direct = codon')
    out_c = carma_final_codon()
        
    return
  end subroutine carma_final

  subroutine carma_final_native
    implicit none

    return
  end subroutine carma_final_native
  

  subroutine carma_timestep_init
    integer(c_int64_t) :: out_c

    call carma_intr_select_impl()

    if (use_native_carma_intr_impl) then
       call carma_timestep_init_native()
       return
    end if

    call carma_intr_log_direct(carma_timestep_init_logged, 'carma_timestep_init direct = codon')
    out_c = carma_timestep_init_codon()

    return
  end subroutine carma_timestep_init

  subroutine carma_timestep_init_native
    implicit none

    return
  end subroutine carma_timestep_init_native


  subroutine carma_timestep_tend(state, cam_in, cam_out, ptend, dt, pbuf, dlf, rliq, prec_str, snow_str, &
    prec_sed, snow_sed, ustar, obklen)
    use hycoef,           only: hyai, hybi, hyam, hybm
    use time_manager,     only: get_nstep, get_step_size, is_first_step
    use camsrfexch,       only: cam_in_t, cam_out_t
    use scamMod,          only: single_column
 
    implicit none

    type(physics_state), intent(inout) :: state                 !! physics state variables
    type(cam_in_t), intent(in)         :: cam_in                !! surface inputs
    type(cam_out_t), intent(inout)     :: cam_out               !! cam output to surface models
    type(physics_ptend), intent(out)   :: ptend                 !! constituent tendencies
    real(r8), intent(in)               :: dt                    !! time step (s)
    type(physics_buffer_desc), pointer :: pbuf(:)               !! physics buffer
    real(r8), intent(in), optional     :: dlf(pcols,pver)       !! Detraining cld H20 from convection (kg/kg/s)
    real(r8), intent(inout), optional  :: rliq(pcols)           !! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), target, intent(out), optional    :: prec_str(pcols)       !! [Total] sfc flux of precip from stratiform (m/s)
    real(r8), target, intent(out), optional    :: snow_str(pcols)       !! [Total] sfc flux of snow from stratiform   (m/s)
    real(r8), target, intent(out), optional    :: prec_sed(pcols)       !! total precip from cloud sedimentation (m/s)
    real(r8), target, intent(out), optional    :: snow_sed(pcols)       !! snow from cloud ice sedimentation (m/s)
    real(r8), intent(in), optional     :: ustar(pcols)          !! friction velocity (m/s)
    real(r8), intent(in), optional     :: obklen(pcols)         !! Obukhov length [ m ]
    type(c_ptr) :: prec_str_p, snow_str_p, prec_sed_p, snow_sed_p
    integer(c_int64_t) :: prec_str_present, snow_str_present, prec_sed_present, snow_sed_present

    call carma_intr_select_impl()
    
    call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update

    if (use_native_carma_intr_impl) then
       if (present(prec_str))  prec_str(:)    = 0._r8
       if (present(snow_str))  snow_str(:)    = 0._r8
       if (present(prec_sed))  prec_sed(:)    = 0._r8
       if (present(snow_sed))  snow_sed(:)    = 0._r8
       return
    end if

    call carma_intr_log_direct(carma_timestep_tend_logged, 'carma_timestep_tend direct = codon')
    prec_str_present = 0_c_int64_t
    snow_str_present = 0_c_int64_t
    prec_sed_present = 0_c_int64_t
    snow_sed_present = 0_c_int64_t
    prec_str_p = c_null_ptr
    snow_str_p = c_null_ptr
    prec_sed_p = c_null_ptr
    snow_sed_p = c_null_ptr
    if (present(prec_str)) then
       prec_str_present = 1_c_int64_t
       prec_str_p = c_loc(prec_str(1))
    end if
    if (present(snow_str)) then
       snow_str_present = 1_c_int64_t
       snow_str_p = c_loc(snow_str(1))
    end if
    if (present(prec_sed)) then
       prec_sed_present = 1_c_int64_t
       prec_sed_p = c_loc(prec_sed(1))
    end if
    if (present(snow_sed)) then
       snow_sed_present = 1_c_int64_t
       snow_sed_p = c_loc(snow_sed(1))
    end if
    call carma_timestep_tend_codon(int(pcols, c_int64_t), prec_str_present, snow_str_present, &
         prec_sed_present, snow_sed_present, prec_str_p, snow_str_p, prec_sed_p, snow_sed_p)

    return
  end subroutine carma_timestep_tend


  subroutine carma_init_cnst(name, q, gcid)
    implicit none

    character(len=*), intent(in) :: name               !! constituent name
    real(r8), intent(out)        :: q(plon,plev,plat)  !! mass mixing ratio
    integer, intent(in)          :: gcid(:)            !! global column id

    call carma_intr_touch()
    
    if (name == "carma") then
      q = 0._r8
    end if 
    
    return
  end subroutine carma_init_cnst


  subroutine carma_emission_tend(state, ptend, cam_in, dt)
    use camsrfexch,       only: cam_in_t

    implicit none
    
    type(physics_state), intent(in )    :: state                !! physics state
    type(physics_ptend), intent(inout)  :: ptend                !! physics state tendencies
    type(cam_in_t),      intent(inout)  :: cam_in               !! surface inputs
    real(r8),            intent(in)     :: dt                   !! time step (s)

    call carma_intr_touch()

    return
  end subroutine carma_emission_tend 


  subroutine carma_wetdep_tend(state, ptend, dt,  pbuf, dlf, cam_out)
    use camsrfexch,       only: cam_out_t

    implicit none

    real(r8),             intent(in)    :: dt             !! time step (s)
    type(physics_state),  intent(in )   :: state          !! physics state
    type(physics_ptend),  intent(inout) :: ptend          !! physics state tendencies
    type(physics_buffer_desc), pointer  :: pbuf(:)        !! physics buffer
    real(r8), intent(in)                :: dlf(pcols,pver)       !! Detraining cld H20 from convection (kg/kg/s)
    type(cam_out_t),      intent(inout) :: cam_out        !! cam output to surface models

    call carma_intr_touch()

    return
  end subroutine carma_wetdep_tend


  subroutine carma_accumulate_stats()
    integer(c_int64_t) :: out_c

    call carma_intr_select_impl()

    if (use_native_carma_intr_impl) then
       call carma_accumulate_stats_native()
       return
    end if

    call carma_intr_log_direct(carma_accumulate_stats_logged, 'carma_accumulate_stats direct = codon')
    out_c = carma_accumulate_stats_codon()

  end subroutine carma_accumulate_stats

  subroutine carma_accumulate_stats_native()
    implicit none

  end subroutine carma_accumulate_stats_native
end module carma_intr
