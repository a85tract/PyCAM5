module iondrag
  !-------------------------------------------------------------------------------
  !  Dummy interface for waccm/iondrag module
  !-------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use iso_c_binding, only: c_int64_t
  use ppgrid       ,only: pver
  use physics_types,only: physics_state, physics_ptend
  use physics_buffer ,only: physics_buffer_desc
  use spmd_utils, only: masterproc
  use cam_logfile, only: iulog

  implicit none

  save

  private                         ! Make default type private to the module

  !-------------------------------------------------------------------------------
  ! Public interfaces:
  !-------------------------------------------------------------------------------
  public :: iondrag_register         ! Register variables in pbuf physics buffer
  public :: iondrag_init             ! Initialization
  public :: iondrag_calc             ! ion drag tensors lxx,lyy,lxy,lyx
  public :: iondrag_readnl
  public :: do_waccm_ions

  interface iondrag_calc
     module procedure iondrag_calc_ions
     module procedure iondrag_calc_ghg
  end interface

  logical, parameter :: do_waccm_ions = .false.

  logical :: use_native_iondrag_impl = .false.
  logical :: iondrag_impl_selected = .false.
  logical :: iondrag_proof_written = .false.

  interface
     function iondrag_touch_codon() result(out_c) bind(c, name="iondrag_touch_codon")
        use iso_c_binding, only: c_int64_t
        integer(c_int64_t) :: out_c
     end function iondrag_touch_codon
  end interface

contains

  !================================================================================================

  subroutine iondrag_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (iondrag_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('IONDRAG_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_iondrag_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_iondrag_impl = .false.
    end if

    iondrag_impl_selected = .true.

    if (masterproc) then
       if (use_native_iondrag_impl) then
          write(iulog,*) 'iondrag implementation = native'
       else
          write(iulog,*) 'iondrag implementation = codon'
       end if
    end if

  end subroutine iondrag_select_impl

  !================================================================================================

  subroutine iondrag_proof_once()

    if (iondrag_proof_written) return
    iondrag_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'iondrag entered (dummy no-op helpers = codon)'
    end if

  end subroutine iondrag_proof_once

  !================================================================================================

  subroutine iondrag_touch()

    integer(c_int64_t) :: out_c

    if (iondrag_proof_written) return

    call iondrag_select_impl()

    if (use_native_iondrag_impl) return

    call iondrag_proof_once()
    out_c = iondrag_touch_codon()
    if (out_c /= 0_c_int64_t) then
       return
    end if

  end subroutine iondrag_touch

  !================================================================================================

  subroutine iondrag_readnl(nlfile)

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    call iondrag_touch()

  end subroutine iondrag_readnl

  !==============================================================================     

  subroutine iondrag_register

    call iondrag_touch()

  end subroutine iondrag_register

  !================================================================================================

  subroutine iondrag_init( pref_mid )
   
    !-------------------------------------------------------------------------------
    ! dummy arguments
    !-------------------------------------------------------------------------------
    real(r8), intent(in) :: pref_mid(pver)

    call iondrag_touch()

  end subroutine iondrag_init

  !================================================================================================
  subroutine iondrag_calc_ions( lchnk, ncol, state, ptend, pbuf, delt )

    !-------------------------------------------------------------------------------
    ! dummy arguments
    !-------------------------------------------------------------------------------
    integer,intent(in)   :: lchnk               ! current chunk index
    integer,intent(in)   :: ncol                ! number of atmospheric columns
    real(r8), intent(in) :: delt                ! time step (s)
    type(physics_state), intent(in), target    :: state ! Physics state variables
    type(physics_ptend), intent(out)   :: ptend   ! Physics tendencies
    type(physics_buffer_desc), pointer :: pbuf(:) ! physics buffer

    call iondrag_touch()

  end subroutine iondrag_calc_ions

  !=========================================================================

  subroutine iondrag_calc_ghg (lchnk,ncol,state,ptend)

    !--------------------Input arguments------------------------------------

    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns

    type(physics_state), intent(in) :: state
    type(physics_ptend), intent(out):: ptend

    call iondrag_touch()

  end subroutine iondrag_calc_ghg

  !===================================================================================

end module iondrag
