module aero_model_wetdep_host_hooks

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  abstract interface
     subroutine timer_hook(timer_name)
       character(len=*), intent(in) :: timer_name
     end subroutine timer_hook

     subroutine outfld_real1d_hook(field_index, cloudborne, suffix, field, &
          dim1, lchnk)
       import :: r8
       integer, intent(in) :: field_index
       logical, intent(in) :: cloudborne
       character(len=*), intent(in) :: suffix
       real(r8), intent(in) :: field(:)
       integer, intent(in) :: dim1, lchnk
     end subroutine outfld_real1d_hook

     subroutine outfld_real2d_hook(field_index, cloudborne, suffix, field, &
          dim1, lchnk)
       import :: r8
       integer, intent(in) :: field_index
       logical, intent(in) :: cloudborne
       character(len=*), intent(in) :: suffix
       real(r8), intent(in) :: field(:,:)
       integer, intent(in) :: dim1, lchnk
     end subroutine outfld_real2d_hook
  end interface

  procedure(timer_hook), pointer :: timer_start_impl => null()
  procedure(timer_hook), pointer :: timer_stop_impl => null()
  procedure(outfld_real1d_hook), pointer :: outfld_real1d_impl => null()
  procedure(outfld_real2d_hook), pointer :: outfld_real2d_impl => null()

  public :: register_aero_model_wetdep_host_hooks
  public :: aero_model_wetdep_timer_start, aero_model_wetdep_timer_stop
  public :: aero_model_wetdep_outfld_real1d, aero_model_wetdep_outfld_real2d

contains

  subroutine register_aero_model_wetdep_host_hooks(timer_start, timer_stop, &
       outfld_real1d, outfld_real2d)
    procedure(timer_hook) :: timer_start, timer_stop
    procedure(outfld_real1d_hook) :: outfld_real1d
    procedure(outfld_real2d_hook) :: outfld_real2d

    timer_start_impl => timer_start
    timer_stop_impl => timer_stop
    outfld_real1d_impl => outfld_real1d
    outfld_real2d_impl => outfld_real2d
  end subroutine register_aero_model_wetdep_host_hooks

  subroutine aero_model_wetdep_timer_start(timer_name)
    character(len=*), intent(in) :: timer_name
    if (associated(timer_start_impl)) call timer_start_impl(timer_name)
  end subroutine aero_model_wetdep_timer_start

  subroutine aero_model_wetdep_timer_stop(timer_name)
    character(len=*), intent(in) :: timer_name
    if (associated(timer_stop_impl)) call timer_stop_impl(timer_name)
  end subroutine aero_model_wetdep_timer_stop

  subroutine aero_model_wetdep_outfld_real1d(field_index, cloudborne, &
       suffix, field, dim1, lchnk)
    integer, intent(in) :: field_index
    logical, intent(in) :: cloudborne
    character(len=*), intent(in) :: suffix
    real(r8), intent(in) :: field(:)
    integer, intent(in) :: dim1, lchnk

    if (associated(outfld_real1d_impl)) then
       call outfld_real1d_impl(field_index, cloudborne, suffix, field, &
            dim1, lchnk)
    end if
  end subroutine aero_model_wetdep_outfld_real1d

  subroutine aero_model_wetdep_outfld_real2d(field_index, cloudborne, &
       suffix, field, dim1, lchnk)
    integer, intent(in) :: field_index
    logical, intent(in) :: cloudborne
    character(len=*), intent(in) :: suffix
    real(r8), intent(in) :: field(:,:)
    integer, intent(in) :: dim1, lchnk

    if (associated(outfld_real2d_impl)) then
       call outfld_real2d_impl(field_index, cloudborne, suffix, field, &
            dim1, lchnk)
    end if
  end subroutine aero_model_wetdep_outfld_real2d

end module aero_model_wetdep_host_hooks
