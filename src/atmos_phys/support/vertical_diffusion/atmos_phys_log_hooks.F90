module atmos_phys_log_hooks

  implicit none
  private

  abstract interface
    subroutine log_message_hook(message)
      character(len=*), intent(in) :: message
    end subroutine log_message_hook
  end interface

  procedure(log_message_hook), pointer :: production_log => null()

  public :: register_atmos_phys_log_hook
  public :: atmos_phys_log

contains

  subroutine register_atmos_phys_log_hook(log_proc)
    procedure(log_message_hook) :: log_proc

    production_log => log_proc
  end subroutine register_atmos_phys_log_hook

  subroutine atmos_phys_log(message)
    character(len=*), intent(in) :: message

    if (associated(production_log)) call production_log(message)
  end subroutine atmos_phys_log

end module atmos_phys_log_hooks
