module ap_drydep_update_scheme

  implicit none
  private

  public :: drydep_update_run

contains

  !> \section arg_table_drydep_update_run Argument Table
  !! \htmlinclude drydep_update_run.html
  !!
  subroutine drydep_update_run(number_drydep_species, drydep_method, &
       land_drydep_method, do_update, errmsg, errflg)
    integer, intent(in) :: number_drydep_species
    character(len=*), intent(in) :: drydep_method
    character(len=*), intent(in) :: land_drydep_method
    logical, intent(out) :: do_update
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    errmsg = ''
    errflg = 0
    do_update = number_drydep_species >= 1 .and. drydep_method == land_drydep_method
  end subroutine drydep_update_run

end module ap_drydep_update_scheme
