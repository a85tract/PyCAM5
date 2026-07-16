module ap_modal_aerosol_timer_hooks

  implicit none
  private

  abstract interface
     subroutine modal_aerosol_timer_hook(event, handle)
       character(len=*), intent(in) :: event
       integer, optional :: handle
     end subroutine modal_aerosol_timer_hook
  end interface

  procedure(modal_aerosol_timer_hook), pointer :: start_hook => null()
  procedure(modal_aerosol_timer_hook), pointer :: stop_hook => null()

  public :: register_modal_aerosol_timer_hooks
  public :: modal_aerosol_timer_start
  public :: modal_aerosol_timer_stop

contains

  subroutine register_modal_aerosol_timer_hooks(start_proc, stop_proc)
    procedure(modal_aerosol_timer_hook) :: start_proc
    procedure(modal_aerosol_timer_hook) :: stop_proc

    start_hook => start_proc
    stop_hook => stop_proc

  end subroutine register_modal_aerosol_timer_hooks

  subroutine modal_aerosol_timer_start(event)
    character(len=*), intent(in) :: event

    if (associated(start_hook)) call start_hook(event)

  end subroutine modal_aerosol_timer_start

  subroutine modal_aerosol_timer_stop(event)
    character(len=*), intent(in) :: event

    if (associated(stop_hook)) call stop_hook(event)

  end subroutine modal_aerosol_timer_stop

end module ap_modal_aerosol_timer_hooks
