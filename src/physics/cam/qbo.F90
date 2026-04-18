module qbo

! Stub version of qbo module

use cam_logfile, only: iulog
use spmd_utils,  only: masterproc

implicit none
private
save

logical :: use_native_impl = .false.
logical :: impl_selected = .false.
logical :: use_native_tstep_impl = .false.
logical :: tstep_impl_selected = .false.

!---------------------------------------------------------------------
! Public methods
!---------------------------------------------------------------------
public               :: qbo_readnl             ! read namelist
public               :: qbo_init               ! initialize qbo package
public               :: qbo_timestep_init      ! interpolate to current time
public               :: qbo_relax              ! relax zonal mean wind

logical, public, parameter :: qbo_use_forcing  = .FALSE.

contains

subroutine qbo_relax_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('QBO_RELAX_IMPL', value=impl_name, length=n, status=status)

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
        write(iulog,*) 'qbo_relax implementation = native'
     else
        write(iulog,*) 'qbo_relax implementation = codon'
     end if
  end if

end subroutine qbo_relax_select_impl

subroutine qbo_timestep_init_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tstep_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('QBO_TSTEP_INIT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tstep_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tstep_impl = .false.
  end if

  tstep_impl_selected = .true.

  if (masterproc) then
     if (use_native_tstep_impl) then
        write(iulog,*) 'qbo_timestep_init implementation = native'
     else
        write(iulog,*) 'qbo_timestep_init implementation = codon'
     end if
  end if

end subroutine qbo_timestep_init_select_impl

subroutine qbo_readnl(nlfile)

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Stub; do nothing.

end subroutine qbo_readnl

subroutine qbo_init

  ! Stub; do nothing.

end subroutine qbo_init

subroutine qbo_timestep_init
  interface
     subroutine qbo_tstep_init_codon() bind(c, name="qbo_tstep_init_codon")
     end subroutine qbo_tstep_init_codon
  end interface

  call qbo_timestep_init_select_impl()

  if (use_native_tstep_impl) then
     call qbo_timestep_init_native()
     return
  end if

  call qbo_tstep_init_codon()

end subroutine qbo_timestep_init

subroutine qbo_timestep_init_native

  ! Stub; do nothing.

end subroutine qbo_timestep_init_native

subroutine qbo_relax( state, pbuf, ptend )

  use iso_c_binding, only: c_int64_t
  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc

  interface
     subroutine qbo_relax_codon() bind(c, name="qbo_relax_codon")
     end subroutine qbo_relax_codon
  end interface

!--------------------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------------------
  type(physics_state), intent(in)    :: state                ! Physics state variables
  type(physics_buffer_desc), pointer :: pbuf(:)              ! Physics buffer
  type(physics_ptend), intent(out)   :: ptend                ! individual parameterization tendencies

  call qbo_relax_select_impl()

  if (use_native_impl) then
     call qbo_relax_native(state, pbuf, ptend)
     return
  end if

  call physics_ptend_init(ptend, state%psetcols, 'qbo (stub)')
  call qbo_relax_codon()

end subroutine qbo_relax

subroutine qbo_relax_native( state, pbuf, ptend )

  use physics_types,  only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only: physics_buffer_desc
!--------------------------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------------------------
  type(physics_state), intent(in)    :: state                ! Physics state variables
  type(physics_buffer_desc), pointer :: pbuf(:)              ! Physics buffer
  type(physics_ptend), intent(out)   :: ptend                ! individual parameterization tendencies

  ! Stub; do nothing except init unused ptend.
  call physics_ptend_init(ptend, state%psetcols, 'qbo (stub)')

end subroutine qbo_relax_native

end module qbo
