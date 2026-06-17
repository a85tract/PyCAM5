!! This module handles reading the namelist and provides access to some other flags
!! that control a specific CARMA model's behavior.
!!
!! By default the specific CARMA model does not have any unique namelist values. If
!! a CARMA model wishes to have its own namelist, then this file needs to be copied
!! from physics/cam to physics/model/<model_name> and the code needed to read in the
!! namelist values added there. This file will take the place of the one in
!! physics/cam. 
!!
!! It needs to be in its own file to resolve some circular dependencies.
!!
!! @author  Chuck Bardeen
!! @version Mar-2011
module carma_model_flags_mod

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use iso_c_binding,  only: c_int64_t

  ! Flags for integration with CAM Microphysics
  public carma_model_readnl                   ! read the carma model namelist
  

  ! Namelist flags
  !
  ! Create a public definition of any new namelist variables that you wish to have,
  ! and default them to an inital value.
  logical, public                :: carma_flag        = .false.   ! If .true. then turn on CARMA microphysics in CAM
  real(r8), public               :: carma_vf_const    = 0.0_r8    ! If specified and non-zero, constant fall velocity for all particles [cm/s]
  character(len=256), public     :: carma_reftfile    = 'carma_reft.nc'  ! path to the file containing the reference temperature profile

  logical :: use_native_carma_model_flags_impl = .false.
  logical :: carma_model_flags_impl_selected = .false.
  logical :: carma_model_flags_proof_written = .false.
  logical :: carma_model_readnl_logged = .false.

  private :: carma_model_flags_select_impl, carma_model_flags_proof_once, carma_model_flags_touch

  interface
    function carma_flags_touch_codon() result(out_c) bind(c, name="carma_flags_touch_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_flags_touch_codon
    function carma_model_readnl_codon() result(out_c) bind(c, name="carma_model_readnl_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
    end function carma_model_readnl_codon
  end interface

contains

  subroutine carma_model_flags_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (carma_model_flags_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('CARMA_FLAGS_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_carma_model_flags_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_carma_model_flags_impl = .false.
    end if

    carma_model_flags_impl_selected = .true.

    if (masterproc) then
       if (use_native_carma_model_flags_impl) then
          write(iulog,*) 'carma_model_flags implementation = native'
       else
          write(iulog,*) 'carma_model_flags implementation = codon'
       end if
    end if

  end subroutine carma_model_flags_select_impl

  !================================================================================================
  !================================================================================================
  subroutine carma_model_flags_proof_once()

    if (carma_model_flags_proof_written) return
    carma_model_flags_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'carma_model_flags entered (model runtime flags = codon)'
    end if

  end subroutine carma_model_flags_proof_once

  !================================================================================================
  !================================================================================================
  subroutine carma_model_flags_log_direct(logged, proof_line)

    logical, intent(inout) :: logged
    character(len=*), intent(in) :: proof_line

    if (logged) return
    logged = .true.

    if (masterproc) then
       write(iulog,'(A)') trim(proof_line)
    end if

  end subroutine carma_model_flags_log_direct

  !================================================================================================
  !================================================================================================
  subroutine carma_model_flags_touch()

    integer(c_int64_t) :: out_c

    call carma_model_flags_select_impl()

    if (use_native_carma_model_flags_impl) then
       return
    end if

    call carma_model_flags_proof_once()
    out_c = carma_flags_touch_codon()

  end subroutine carma_model_flags_touch

  !================================================================================================
  !================================================================================================

  !! Read the CARMA model runtime options from the namelist
  !!
  !! @author  Chuck Bardeen
  !! @version Mar-2011
  subroutine carma_model_readnl(nlfile)
  
    ! Read carma namelist group.
  
    use cam_abortutils,  only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
  
    ! args
  
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
  
    ! local vars
  
    integer :: unitn, ierr
    integer(c_int64_t) :: out_c
  
    ! read namelist for CARMA
!    namelist /carma_model_nl/ &
!      carma_flag, &
!      carma_maxretries, &
!      carma_conmax, &
!      carma_reftfile
  
!    if (masterproc) then
!       unitn = getunit()
!       open( unitn, file=trim(nlfile), status='old' )
!       call find_group_name(unitn, 'carma_model_nl', status=ierr)
!       if (ierr == 0) then
!          read(unitn, carma_model_nl, iostat=ierr)
!          if (ierr /= 0) then
!             call endrun('carma_model_readnl: ERROR reading namelist')
!          end if
!       end if
!       close(unitn)
!       call freeunit(unitn)
!    end if
  
#ifdef SPMD
!    call mpibcast (carma_flag,            1 ,mpilog, 0,mpicom)
!    call mpibcast (carma_maxretries,      1 ,mpiint, 0,mpicom)
!    call mpibcast (carma_conmax,          1 ,mpir8,  0,mpicom)
!    call mpibcast (carma_reftfile, len(carma_reftfile), mpichar, 0, mpicom)
#endif

    call carma_model_flags_select_impl()
    if (use_native_carma_model_flags_impl) then
       call carma_model_flags_touch()
    else
       call carma_model_flags_proof_once()
       out_c = carma_model_readnl_codon()
       call carma_model_flags_log_direct(carma_model_readnl_logged, 'carma_model_readnl direct = codon')
    end if
  
  end subroutine carma_model_readnl

end module carma_model_flags_mod
