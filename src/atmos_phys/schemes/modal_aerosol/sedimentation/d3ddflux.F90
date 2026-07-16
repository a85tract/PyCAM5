module ap_d3ddflux_scheme

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: d3ddflux_run

contains

  !> \section arg_table_d3ddflux_run Argument Table
  !! \htmlinclude d3ddflux_run.html
  subroutine d3ddflux_run(pcols, pver, ncol, vlc_dry, q, pmid, pdel, tv, &
       dep_dry, dep_dry_tend, dt, rair, gravit)
    integer, intent(in) :: pcols, pver, ncol
    real(r8), intent(in) :: vlc_dry(pcols,pver)
    real(r8), intent(in) :: q(pcols,pver)
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: tv(pcols,pver)
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: rair
    real(r8), intent(in) :: gravit
    real(r8), intent(out) :: dep_dry(pcols)
    real(r8), intent(out) :: dep_dry_tend(pcols,pver)

    real(r8) :: flux(pcols,0:pver)
    integer :: i, k

    do i=1,ncol
       flux(i,0)=0._r8
    enddo
    do k=1,pver
       do i = 1, ncol
          flux(i,k) = -min(vlc_dry(i,k) * q(i,k) * pmid(i,k) /(tv(i,k) * rair), &
                    q(i,k)*pdel(i,k)/gravit/dt)
          dep_dry_tend(i,k)=(flux(i,k)-flux(i,k-1))/pdel(i,k)*gravit
       end do
    enddo
    do i=1,ncol
       dep_dry(i)=flux(i,pver)
    enddo
    return
  end subroutine d3ddflux_run

end module ap_d3ddflux_scheme
