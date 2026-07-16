module ap_dycore_energy_consistency_adjust_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private
  public :: dycore_energy_consistency_adjust_run

contains

  !> \section arg_table_dycore_energy_consistency_adjust_run Argument Table
  !! \htmlinclude dycore_energy_consistency_adjust_run.html
  !!
  subroutine dycore_energy_consistency_adjust_run(ncol, pver, do_adjust, &
       scaling_dycore, tend_dtdt, tend_dtdt_local, errmsg, errflg)
    integer,  intent(in)  :: ncol, pver
    logical,  intent(in)  :: do_adjust
    real(r8), intent(in)  :: scaling_dycore(:,:), tend_dtdt(:,:)
    real(r8), intent(out) :: tend_dtdt_local(:,:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0
    tend_dtdt_local(:ncol,:pver) = 0._r8

    ! This branch is intentionally disabled for the old PI-atm energy formula.
    ! It exists as an explicit lifecycle boundary, not as an import of the new
    ! CESM cp/cv correction.  If enabled in a future case, it uses the standard
    ! tendency form without changing today's PI result.
    if (do_adjust) then
       tend_dtdt_local(:ncol,:pver) = &
            (scaling_dycore(:ncol,:pver) - 1._r8) * tend_dtdt(:ncol,:pver)
    end if
  end subroutine dycore_energy_consistency_adjust_run

end module ap_dycore_energy_consistency_adjust_scheme
