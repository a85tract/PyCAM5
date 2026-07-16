module dropmixnuc_host_hooks

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  abstract interface
     subroutine timer_hook(timer_name)
       character(len=*), intent(in) :: timer_name
     end subroutine timer_hook

     subroutine outfld_real1d_hook(field_name, field, dim1, lchnk)
       import :: r8
       character(len=*), intent(in) :: field_name
       real(r8), intent(in) :: field(:)
       integer, intent(in) :: dim1, lchnk
     end subroutine outfld_real1d_hook

     subroutine outfld_real2d_hook(field_name, field, dim1, lchnk)
       import :: r8
       character(len=*), intent(in) :: field_name
       real(r8), intent(in) :: field(:,:)
       integer, intent(in) :: dim1, lchnk
     end subroutine outfld_real2d_hook

     subroutine endrun_hook(message)
       character(len=*), intent(in) :: message
     end subroutine endrun_hook
  end interface

  procedure(timer_hook), pointer :: timer_start_impl => null()
  procedure(timer_hook), pointer :: timer_stop_impl => null()
  procedure(outfld_real1d_hook), pointer :: outfld_real1d_impl => null()
  procedure(outfld_real2d_hook), pointer :: outfld_real2d_impl => null()
  procedure(endrun_hook), pointer :: endrun_impl => null()

  public :: register_dropmixnuc_host_hooks
  public :: dropmixnuc_timer_start, dropmixnuc_timer_stop
  public :: dropmixnuc_outfld_real1d, dropmixnuc_outfld_real2d
  public :: dropmixnuc_endrun

contains

  subroutine register_dropmixnuc_host_hooks(timer_start, timer_stop, &
       outfld_real1d, outfld_real2d, endrun_handler)
    procedure(timer_hook) :: timer_start, timer_stop
    procedure(outfld_real1d_hook) :: outfld_real1d
    procedure(outfld_real2d_hook) :: outfld_real2d
    procedure(endrun_hook) :: endrun_handler

    timer_start_impl => timer_start
    timer_stop_impl => timer_stop
    outfld_real1d_impl => outfld_real1d
    outfld_real2d_impl => outfld_real2d
    endrun_impl => endrun_handler
  end subroutine register_dropmixnuc_host_hooks

  subroutine dropmixnuc_timer_start(timer_name)
    character(len=*), intent(in) :: timer_name
    if (associated(timer_start_impl)) call timer_start_impl(timer_name)
  end subroutine dropmixnuc_timer_start

  subroutine dropmixnuc_timer_stop(timer_name)
    character(len=*), intent(in) :: timer_name
    if (associated(timer_stop_impl)) call timer_stop_impl(timer_name)
  end subroutine dropmixnuc_timer_stop

  subroutine dropmixnuc_outfld_real1d(field_name, field, dim1, lchnk)
    character(len=*), intent(in) :: field_name
    real(r8), intent(in) :: field(:)
    integer, intent(in) :: dim1, lchnk
    if (associated(outfld_real1d_impl)) then
       call outfld_real1d_impl(field_name, field, dim1, lchnk)
    end if
  end subroutine dropmixnuc_outfld_real1d

  subroutine dropmixnuc_outfld_real2d(field_name, field, dim1, lchnk)
    character(len=*), intent(in) :: field_name
    real(r8), intent(in) :: field(:,:)
    integer, intent(in) :: dim1, lchnk
    if (associated(outfld_real2d_impl)) then
       call outfld_real2d_impl(field_name, field, dim1, lchnk)
    end if
  end subroutine dropmixnuc_outfld_real2d

  subroutine dropmixnuc_endrun(message)
    character(len=*), intent(in) :: message
    if (associated(endrun_impl)) then
       call endrun_impl(message)
    else
       error stop 'dropmixnuc fatal error'
    end if
  end subroutine dropmixnuc_endrun

end module dropmixnuc_host_hooks
