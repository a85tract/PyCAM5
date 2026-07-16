module ap_cloud_fraction_fice_scheme

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: cloud_fraction_fice_run

contains

  !> \section arg_table_cloud_fraction_fice_run Argument Table
  !! \htmlinclude cloud_fraction_fice_run.html
  subroutine cloud_fraction_fice_run(pcols, pver, ncol, t, tmelt, top_lev, &
       fice, fsnow)
    integer, intent(in) :: pcols, pver, ncol, top_lev
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: tmelt
    real(r8), intent(out) :: fice(pcols,pver)
    real(r8), intent(out) :: fsnow(pcols,pver)

    real(r8) :: tmax_fice
    real(r8) :: tmin_fice
    real(r8) :: tmax_fsnow
    real(r8) :: tmin_fsnow
    integer :: i, k

    tmax_fice = tmelt - 10._r8
    tmin_fice = tmax_fice - 30._r8
    tmax_fsnow = tmelt
    tmin_fsnow = tmelt - 5._r8

    fice(:,:top_lev-1) = 0._r8
    fsnow(:,:top_lev-1) = 0._r8

    do k=top_lev,pver
       do i=1,ncol
          if (t(i,k) > tmax_fice) then
             fice(i,k) = 0.0_r8
          else if (t(i,k) < tmin_fice) then
             fice(i,k) = 1.0_r8
          else
             fice(i,k) =(tmax_fice - t(i,k)) / (tmax_fice - tmin_fice)
          end if

          if (t(i,k) > tmax_fsnow) then
             fsnow(i,k) = 0.0_r8
          else if (t(i,k) < tmin_fsnow) then
             fsnow(i,k) = 1.0_r8
          else
             fsnow(i,k) =(tmax_fsnow - t(i,k)) / (tmax_fsnow - tmin_fsnow)
          end if
       end do
    end do
  end subroutine cloud_fraction_fice_run

end module ap_cloud_fraction_fice_scheme
