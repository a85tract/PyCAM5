module qneg4_scheme

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: qneg4_run

contains

  !> \section arg_table_qneg4_run Argument Table
  !! \htmlinclude qneg4_run.html
  !!
  subroutine qneg4_run(ncol, nconst, dt, qbot, srfrpdel, qminimum, gravit, &
       latvap, shflx, lhflx, qflx, nptsexc, excess_indices, excess, worst, &
       worst_column, errcode, errmsg)
    integer, intent(in) :: ncol, nconst
    real(r8), intent(in) :: dt, qbot(:,:), srfrpdel(:), qminimum
    real(r8), intent(in) :: gravit, latvap
    real(r8), intent(inout) :: shflx(:), lhflx(:), qflx(:,:)
    integer, intent(out) :: nptsexc, excess_indices(:), worst_column, errcode
    real(r8), intent(out) :: excess(:), worst
    character(len=512), intent(out) :: errmsg

    integer :: i, ii

    errcode = 0
    errmsg = ''

    nptsexc = 0
    do i = 1, ncol
       excess(i) = qflx(i,1) - &
            (qminimum - qbot(i,1))/(dt*gravit*srfrpdel(i))
       if (excess(i) < 0._r8) then
          nptsexc = nptsexc + 1
          excess_indices(nptsexc) = i
          qflx(i,1) = qflx(i,1) - excess(i)
          lhflx(i) = lhflx(i) - excess(i)*latvap
          shflx(i) = shflx(i) + excess(i)*latvap
       end if
    end do

    worst = 0._r8
    worst_column = -1
    if (nptsexc > 10) then
       do ii = 1, nptsexc
          i = excess_indices(ii)
          if (excess(i) < worst) then
             worst = excess(i)
             worst_column = i
          end if
       end do
    end if

  end subroutine qneg4_run

end module qneg4_scheme
