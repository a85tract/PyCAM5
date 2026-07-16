module ap_vertical_diffusion_positive_moisture_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: vertical_diffusion_positive_moisture_run

contains

  !> \section arg_table_vertical_diffusion_positive_moisture_run Argument Table
  !! \htmlinclude vertical_diffusion_positive_moisture_run.html
  !!
  subroutine vertical_diffusion_positive_moisture_run(cp, xlv, xls, ncol, &
       mkx, dt, qvmin, qlmin, qimin, dp, qv, ql, qi, t, s, qvten, qlten, &
       qiten, sten, impossible_count, errmsg, errflg)

    integer,  intent(in)    :: ncol, mkx
    real(r8), intent(in)    :: cp, xlv, xls
    real(r8), intent(in)    :: dt, qvmin, qlmin, qimin
    real(r8), intent(in)    :: dp(ncol,mkx)
    real(r8), intent(inout) :: qv(ncol,mkx)
    real(r8), intent(inout) :: ql(ncol,mkx)
    real(r8), intent(inout) :: qi(ncol,mkx)
    real(r8), intent(inout) :: t(ncol,mkx)
    real(r8), intent(inout) :: s(ncol,mkx)
    real(r8), intent(inout) :: qvten(ncol,mkx)
    real(r8), intent(inout) :: qlten(ncol,mkx)
    real(r8), intent(inout) :: qiten(ncol,mkx)
    real(r8), intent(inout) :: sten(ncol,mkx)
    integer, intent(out) :: impossible_count
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer  :: i, k
    real(r8) :: dql, dqi, dqv, sum, aa, dum

    impossible_count = 0
    errmsg = ''
    errflg = 0

    do i = 1, ncol
       do k = mkx, 1, -1
          dql        = max(0._r8,1._r8*qlmin-ql(i,k))
          dqi        = max(0._r8,1._r8*qimin-qi(i,k))
          qlten(i,k) = qlten(i,k) +  dql/dt
          qiten(i,k) = qiten(i,k) +  dqi/dt
          qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
          sten(i,k)  = sten(i,k)  + xlv * (dql/dt) + xls * (dqi/dt)
          ql(i,k)    = ql(i,k) +  dql
          qi(i,k)    = qi(i,k) +  dqi
          qv(i,k)    = qv(i,k) -  dql - dqi
          s(i,k)     = s(i,k)  +  xlv * dql + xls * dqi
          t(i,k)     = t(i,k)  + (xlv * dql + xls * dqi)/cp
          dqv        = max(0._r8,1._r8*qvmin-qv(i,k))
          qvten(i,k) = qvten(i,k) + dqv/dt
          qv(i,k)    = qv(i,k)    + dqv
          if (k .ne. 1) then
             qv(i,k-1)    = qv(i,k-1)    - dqv*dp(i,k)/dp(i,k-1)
             qvten(i,k-1) = qvten(i,k-1) - dqv*dp(i,k)/dp(i,k-1)/dt
          endif
          qv(i,k) = max(qv(i,k),qvmin)
          ql(i,k) = max(ql(i,k),qlmin)
          qi(i,k) = max(qi(i,k),qimin)
       end do

       if (dqv .gt. 1.e-20_r8) then
          sum = 0._r8
          do k = 1, mkx
             if (qv(i,k) .gt. 2._r8*qvmin) sum = sum + qv(i,k)*dp(i,k)
          enddo
          aa = dqv*dp(i,1)/max(1.e-20_r8,sum)
          if (aa .lt. 0.5_r8) then
             do k = 1, mkx
                if (qv(i,k) .gt. 2._r8*qvmin) then
                   dum        = aa*qv(i,k)
                   qv(i,k)    = qv(i,k) - dum
                   qvten(i,k) = qvten(i,k) - dum/dt
                endif
             enddo
          else
             impossible_count = impossible_count + 1
          endif
       endif
    end do

  end subroutine vertical_diffusion_positive_moisture_run

end module ap_vertical_diffusion_positive_moisture_scheme
