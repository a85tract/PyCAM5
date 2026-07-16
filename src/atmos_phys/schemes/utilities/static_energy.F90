module static_energy

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: update_dry_static_energy_run

contains

  !> \section arg_table_update_dry_static_energy_run Argument Table
  !! \htmlinclude update_dry_static_energy_run.html
  !!
  subroutine update_dry_static_energy_run(ncol, num_layers, gravit, temperature, &
       geopotential_height, surface_geopotential, cpair, include_geopotential, &
       dry_static_energy, errcode, errmsg)
    integer, intent(in) :: ncol, num_layers
    real(r8), intent(in) :: gravit, temperature(:,:), geopotential_height(:,:)
    real(r8), intent(in) :: surface_geopotential(:), cpair
    logical, intent(in) :: include_geopotential
    real(r8), intent(inout) :: dry_static_energy(:,:)
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    integer :: i, k

    errcode = 0
    errmsg = ''
    do k = 1, num_layers
       do i = 1, ncol
          if (include_geopotential) then
             dry_static_energy(i,k) = cpair*temperature(i,k) + &
                  gravit*geopotential_height(i,k) + surface_geopotential(i)
          else
             dry_static_energy(i,k) = cpair*temperature(i,k) + &
                  surface_geopotential(i)
          end if
       end do
    end do
  end subroutine update_dry_static_energy_run

end module static_energy
