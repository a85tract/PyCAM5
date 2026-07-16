module physics_tendency_updaters

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: apply_tendency_of_eastward_wind_run
  public :: apply_tendency_of_northward_wind_run
  public :: apply_heating_rate_run
  public :: apply_constituent_tendencies_run

contains

  !> \section arg_table_apply_tendency_of_eastward_wind_run Argument Table
  !! \htmlinclude apply_tendency_of_eastward_wind_run.html
  !!
  subroutine apply_tendency_of_eastward_wind_run(ncol, top_level, bot_level, &
       dudt, u, dudt_total, dt, accumulate_total, errcode, errmsg)
    integer, intent(in) :: ncol, top_level, bot_level
    real(r8), intent(in) :: dudt(:,:), dt
    real(r8), intent(inout) :: u(:,:), dudt_total(:,:)
    logical, intent(in) :: accumulate_total
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    integer :: k

    errcode = 0
    errmsg = ''
    do k = top_level, bot_level
       u(:ncol,k) = u(:ncol,k) + dudt(:ncol,k) * dt
       if (accumulate_total) &
            dudt_total(:ncol,k) = dudt_total(:ncol,k) + dudt(:ncol,k)
    end do
  end subroutine apply_tendency_of_eastward_wind_run

  !> \section arg_table_apply_tendency_of_northward_wind_run Argument Table
  !! \htmlinclude apply_tendency_of_northward_wind_run.html
  !!
  subroutine apply_tendency_of_northward_wind_run(ncol, top_level, bot_level, &
       dvdt, v, dvdt_total, dt, accumulate_total, errcode, errmsg)
    integer, intent(in) :: ncol, top_level, bot_level
    real(r8), intent(in) :: dvdt(:,:), dt
    real(r8), intent(inout) :: v(:,:), dvdt_total(:,:)
    logical, intent(in) :: accumulate_total
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    integer :: k

    errcode = 0
    errmsg = ''
    do k = top_level, bot_level
       v(:ncol,k) = v(:ncol,k) + dvdt(:ncol,k) * dt
       if (accumulate_total) &
            dvdt_total(:ncol,k) = dvdt_total(:ncol,k) + dvdt(:ncol,k)
    end do
  end subroutine apply_tendency_of_northward_wind_run

  !> \section arg_table_apply_heating_rate_run Argument Table
  !! \htmlinclude apply_heating_rate_run.html
  !!
  subroutine apply_heating_rate_run(ncol, top_level, bot_level, heating_rate, &
       dry_static_energy, dtdt_total, dt, cpair, accumulate_total, errcode, errmsg)
    integer, intent(in) :: ncol, top_level, bot_level
    real(r8), intent(in) :: heating_rate(:,:), dt, cpair(:,:)
    real(r8), intent(inout) :: dry_static_energy(:,:), dtdt_total(:,:)
    logical, intent(in) :: accumulate_total
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    integer :: k

    errcode = 0
    errmsg = ''
    do k = top_level, bot_level
       dry_static_energy(:ncol,k) = dry_static_energy(:ncol,k) + &
            heating_rate(:ncol,k) * dt
       if (accumulate_total) &
            dtdt_total(:ncol,k) = dtdt_total(:ncol,k) + &
            heating_rate(:ncol,k) / cpair(:ncol,k)
    end do
  end subroutine apply_heating_rate_run

  !> \section arg_table_apply_constituent_tendencies_run Argument Table
  !! \htmlinclude apply_constituent_tendencies_run.html
  !!
  subroutine apply_constituent_tendencies_run(ncol, top_level, bot_level, &
       constituent_tendency, constituent, dt, errcode, errmsg)
    integer, intent(in) :: ncol, top_level, bot_level
    real(r8), intent(in) :: constituent_tendency(:,:), dt
    real(r8), intent(inout) :: constituent(:,:)
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    integer :: k

    errcode = 0
    errmsg = ''
    do k = top_level, bot_level
       constituent(:ncol,k) = constituent(:ncol,k) + &
            constituent_tendency(:ncol,k) * dt
    end do
  end subroutine apply_constituent_tendencies_run

end module physics_tendency_updaters
