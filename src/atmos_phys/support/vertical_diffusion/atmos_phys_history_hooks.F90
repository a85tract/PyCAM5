module atmos_phys_history_hooks

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  abstract interface
     subroutine history_real_1d_hook(name, field, horizontal_size, chunk)
       import :: r8
       character(len=*), intent(in) :: name
       real(r8), intent(in) :: field(:)
       integer, intent(in) :: horizontal_size, chunk
     end subroutine history_real_1d_hook

     subroutine history_real_2d_hook(name, field, horizontal_size, chunk)
       import :: r8
       character(len=*), intent(in) :: name
       real(r8), intent(in) :: field(:,:)
       integer, intent(in) :: horizontal_size, chunk
     end subroutine history_real_2d_hook
  end interface

  procedure(history_real_1d_hook), pointer :: real_1d_hook => null()
  procedure(history_real_2d_hook), pointer :: real_2d_hook => null()

  interface atmos_phys_outfld
     module procedure atmos_phys_outfld_real_1d
     module procedure atmos_phys_outfld_real_2d
  end interface atmos_phys_outfld

  public :: register_atmos_phys_history_hooks
  public :: atmos_phys_outfld

contains

  subroutine register_atmos_phys_history_hooks(real_1d_proc, real_2d_proc)
    procedure(history_real_1d_hook) :: real_1d_proc
    procedure(history_real_2d_hook) :: real_2d_proc

    real_1d_hook => real_1d_proc
    real_2d_hook => real_2d_proc
  end subroutine register_atmos_phys_history_hooks

  subroutine atmos_phys_outfld_real_1d(name, field, horizontal_size, chunk)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: field(:)
    integer, intent(in) :: horizontal_size, chunk

    if (associated(real_1d_hook)) then
       call real_1d_hook(name, field, horizontal_size, chunk)
    end if
  end subroutine atmos_phys_outfld_real_1d

  subroutine atmos_phys_outfld_real_2d(name, field, horizontal_size, chunk)
    character(len=*), intent(in) :: name
    real(r8), intent(in) :: field(:,:)
    integer, intent(in) :: horizontal_size, chunk

    if (associated(real_2d_hook)) then
       call real_2d_hook(name, field, horizontal_size, chunk)
    end if
  end subroutine atmos_phys_outfld_real_2d

end module atmos_phys_history_hooks
