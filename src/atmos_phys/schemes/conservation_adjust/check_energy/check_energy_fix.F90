module ap_check_energy_fix_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private
  public :: check_energy_fix_run

contains

  !> \section arg_table_check_energy_fix_run Argument Table
  !! \htmlinclude check_energy_fix_run.html
  !!
  subroutine check_energy_fix_run(ncol, pver, pint, gravit, heat_glob, &
       ptend_s, eshflx, errmsg, errflg)
    integer,  intent(in)    :: ncol, pver
    real(r8), intent(in)    :: pint(:,:), gravit, heat_glob
    real(r8), intent(inout) :: ptend_s(:,:)
    real(r8), intent(out)   :: eshflx(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: i

    errmsg = ''
    errflg = 0

    ptend_s(:ncol,:pver) = heat_glob
    do i = 1, ncol
       eshflx(i) = heat_glob * (pint(i,pver+1) - pint(i,1)) / gravit
    end do
  end subroutine check_energy_fix_run

end module ap_check_energy_fix_scheme
