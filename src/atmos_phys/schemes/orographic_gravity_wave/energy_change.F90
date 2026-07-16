module ap_energy_change_scheme

  use gw_utils, only: r8
  use coords_1d, only: Coords1D

  implicit none
  private

  public :: energy_change_run

contains

  !> \section arg_table_energy_change_run Argument Table
  !! \htmlinclude energy_change_run.html
  !!
subroutine energy_change_run(pver, ncol, gravit, dt, p, u, v, dudt, dvdt, dsdt, de)

  integer, intent(in) :: pver, ncol
  real(r8), intent(in) :: gravit

  ! Time step.
  real(r8), intent(in) :: dt
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Winds at start of time step.
  real(r8), intent(in) :: u(:,:), v(:,:)
  ! Wind tendencies.
  real(r8), intent(in) :: dudt(:,:), dvdt(:,:)
  ! Heating tendency.
  real(r8), intent(in) :: dsdt(:,:)
  ! Change in energy.
  real(r8), intent(out) :: de(:)

  ! Level index.
  integer :: k

  ! Net gain/loss of total energy in the column.
  de = 0.0_r8
  do k = 1, pver
     de = de + p%del(:,k)/gravit * (dsdt(:,k) + &
          dudt(:,k)*(u(:,k)+dudt(:,k)*0.5_r8*dt) + &
          dvdt(:,k)*(v(:,k)+dvdt(:,k)*0.5_r8*dt) )
  end do

end subroutine energy_change_run

end module ap_energy_change_scheme
