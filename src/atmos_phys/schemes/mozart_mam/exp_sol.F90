module ap_exp_sol_scheme

  implicit none
  private

  public :: exp_sol_run

contains

  !> \section arg_table_exp_sol_run Argument Table
  !! \htmlinclude exp_sol_run.html
  !!
  subroutine exp_sol_run(explicit_species_count, entry_count, errmsg, errflg)
    integer, intent(in)    :: explicit_species_count
    integer, intent(inout) :: entry_count
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    entry_count = entry_count + 1
    errmsg = ''
    errflg = 0

    ! The generated pp_trop_mam3 mechanism has clscnt1=0.  Its active
    ! exp_sol call is therefore a numerical no-op, but remains a real process
    ! boundary.  Fail closed if this mechanism-specific scheme is ever wired
    ! to a chemistry package with explicit species.
    if (explicit_species_count /= 0) then
       errmsg = 'exp_sol_run is locked to pp_trop_mam3 (clscnt1=0)'
       errflg = 1
    end if
  end subroutine exp_sol_run

end module ap_exp_sol_scheme
