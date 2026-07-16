module qneg

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: qneg_run

contains

  !> \section arg_table_qneg_run Argument Table
  !! \htmlinclude qneg_run.html
  !!
  subroutine qneg_run(num_columns, declared_columns, num_levels, qminimum, &
       constituent, nvals, worst, worst_column, worst_level, errcode, errmsg)
    integer, intent(in) :: num_columns, declared_columns, num_levels
    real(r8), intent(in) :: qminimum
    real(r8), intent(inout) :: constituent(declared_columns,num_levels)
    integer, intent(out) :: nvals, worst_column, worst_level, errcode
    real(r8), intent(out) :: worst
    character(len=512), intent(out) :: errmsg

    integer :: indx(num_columns,num_levels)
    integer :: nval(num_levels)
    integer :: nn, iwtmp, i, ii, k
    logical :: found

    errcode = 0
    errmsg = ''
    nvals = 0
    found = .false.
    worst = 1.e35_r8
    worst_column = -1
    worst_level = -1

!DIR$ preferstream
    do k = 1, num_levels
       nval(k) = 0
!DIR$ prefervector
       nn = 0
       do i = 1, num_columns
          if (constituent(i,k) < qminimum) then
             nn = nn + 1
             indx(nn,k) = i
          end if
       end do
       nval(k) = nn
    end do

    do k = 1, num_levels
       if (nval(k) > 0) then
          found = .true.
          nvals = nvals + nval(k)
          iwtmp = -1
!cdir nodep,altcode=loopcnt
          do ii = 1, nval(k)
             i = indx(ii,k)
             if (constituent(i,k) < worst) then
                worst = constituent(i,k)
                iwtmp = ii
             end if
          end do
          if (iwtmp /= -1) worst_level = k
          if (iwtmp /= -1) worst_column = indx(iwtmp,k)
!cdir nodep,altcode=loopcnt
          do ii = 1, nval(k)
             i = indx(ii,k)
             constituent(i,k) = qminimum
          end do
       end if
    end do

    ! Keep the original qneg3 interface responsible for warnings.
    if (.not. found) nvals = 0
  end subroutine qneg_run

end module qneg
