module ap_set_dry_to_wet_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private
  public :: set_dry_to_wet_run

contains

  !> \section arg_table_set_dry_to_wet_run Argument Table
  !! \htmlinclude set_dry_to_wet_run.html
  !!
  subroutine set_dry_to_wet_run(ncol, pver, pcnst, constituent_is_dry, &
       pdel, pdeldry, q, errmsg, errflg)
    integer,  intent(in)    :: ncol, pver, pcnst
    logical,  intent(in)    :: constituent_is_dry(:)
    real(r8), intent(in)    :: pdel(:,:), pdeldry(:,:)
    real(r8), intent(inout) :: q(:,:,:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: m

    errmsg = ''
    errflg = 0

    do m = 1, pcnst
       if (constituent_is_dry(m)) then
          q(:ncol,:pver,m) = q(:ncol,:pver,m)*pdeldry(:ncol,:pver)/pdel(:ncol,:pver)
       end if
    end do
  end subroutine set_dry_to_wet_run

end module ap_set_dry_to_wet_scheme
