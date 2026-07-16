module cloud_microphysics_host_hooks

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  abstract interface
     subroutine timer_hook(timer_name)
       character(len=*), intent(in) :: timer_name
     end subroutine timer_hook

     subroutine outfld_real2d_hook(field_name, field, dim1, lchnk)
       import :: r8
       character(len=*), intent(in) :: field_name
       real(r8), intent(in) :: field(:,:)
       integer, intent(in) :: dim1, lchnk
     end subroutine outfld_real2d_hook

     logical function history_active_hook(field_name)
       character(len=*), intent(in) :: field_name
     end function history_active_hook

     subroutine endrun_hook(message)
       character(len=*), intent(in) :: message
     end subroutine endrun_hook
  end interface

  procedure(timer_hook), pointer :: timer_start_impl => null()
  procedure(timer_hook), pointer :: timer_stop_impl => null()
  procedure(outfld_real2d_hook), pointer :: outfld_real2d_impl => null()
  procedure(history_active_hook), pointer :: history_active_impl => null()
  procedure(endrun_hook), pointer :: endrun_impl => null()

  public :: register_cloud_microphysics_host_hooks
  public :: cloud_microphysics_timer_start, cloud_microphysics_timer_stop
  public :: cloud_microphysics_outfld_real2d
  public :: cloud_microphysics_history_active
  public :: cloud_microphysics_endrun

contains

  subroutine register_cloud_microphysics_host_hooks(timer_start, timer_stop, &
       outfld_real2d, history_active, endrun_handler)
    procedure(timer_hook) :: timer_start, timer_stop
    procedure(outfld_real2d_hook) :: outfld_real2d
    procedure(history_active_hook) :: history_active
    procedure(endrun_hook) :: endrun_handler

    timer_start_impl => timer_start
    timer_stop_impl => timer_stop
    outfld_real2d_impl => outfld_real2d
    history_active_impl => history_active
    endrun_impl => endrun_handler
  end subroutine register_cloud_microphysics_host_hooks

  subroutine cloud_microphysics_timer_start(timer_name)
    character(len=*), intent(in) :: timer_name

    if (associated(timer_start_impl)) call timer_start_impl(timer_name)
  end subroutine cloud_microphysics_timer_start

  subroutine cloud_microphysics_timer_stop(timer_name)
    character(len=*), intent(in) :: timer_name

    if (associated(timer_stop_impl)) call timer_stop_impl(timer_name)
  end subroutine cloud_microphysics_timer_stop

  subroutine cloud_microphysics_outfld_real2d(field_name, field, dim1, lchnk)
    character(len=*), intent(in) :: field_name
    real(r8), intent(in) :: field(:,:)
    integer, intent(in) :: dim1, lchnk

    if (associated(outfld_real2d_impl)) then
       call outfld_real2d_impl(field_name, field, dim1, lchnk)
    end if
  end subroutine cloud_microphysics_outfld_real2d

  logical function cloud_microphysics_history_active(field_name)
    character(len=*), intent(in) :: field_name

    cloud_microphysics_history_active = .false.
    if (associated(history_active_impl)) then
       cloud_microphysics_history_active = history_active_impl(field_name)
    end if
  end function cloud_microphysics_history_active

  subroutine cloud_microphysics_endrun(message)
    character(len=*), intent(in) :: message

    if (associated(endrun_impl)) then
       call endrun_impl(message)
    else
       error stop 'cloud microphysics fatal error'
    end if
  end subroutine cloud_microphysics_endrun

end module cloud_microphysics_host_hooks
