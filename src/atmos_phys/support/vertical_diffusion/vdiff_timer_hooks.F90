module vdiff_timer_hooks

  implicit none
  private

  abstract interface
     subroutine vdiff_timer_hook(event, handle)
       character(len=*), intent(in) :: event
       integer, optional :: handle
     end subroutine vdiff_timer_hook
  end interface

  procedure(vdiff_timer_hook), pointer :: start_hook => null()
  procedure(vdiff_timer_hook), pointer :: stop_hook => null()

  public :: register_vdiff_timer_hooks
  public :: vdiff_timer_start
  public :: vdiff_timer_stop

contains

  subroutine register_vdiff_timer_hooks(start_proc, stop_proc)
    procedure(vdiff_timer_hook) :: start_proc
    procedure(vdiff_timer_hook) :: stop_proc

    start_hook => start_proc
    stop_hook => stop_proc

  end subroutine register_vdiff_timer_hooks

  subroutine vdiff_timer_start(event)
    character(len=*), intent(in) :: event

    if (associated(start_hook)) call start_hook(event)

  end subroutine vdiff_timer_start

  subroutine vdiff_timer_stop(event)
    character(len=*), intent(in) :: event

    if (associated(stop_hook)) call stop_hook(event)

  end subroutine vdiff_timer_stop

end module vdiff_timer_hooks
