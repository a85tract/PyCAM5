module ap_calc_obklen_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: calc_obklen_run

contains

  !> \section arg_table_calc_obklen_run Argument Table
  !! \htmlinclude calc_obklen_run.html
  !!
  subroutine calc_obklen_run(ncol, g, vk, cpair, zvir, ths, thvs, qflx, &
       shflx, rrho, ustar, khfs, kqfs, kbfs, obklen, errmsg, errflg)

    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: g
    real(r8), intent(in)  :: vk
    real(r8), intent(in)  :: cpair
    real(r8), intent(in)  :: zvir
    real(r8), intent(in)  :: ths(ncol)
    real(r8), intent(in)  :: thvs(ncol)
    real(r8), intent(in)  :: qflx(ncol)
    real(r8), intent(in)  :: shflx(ncol)
    real(r8), intent(in)  :: rrho(ncol)
    real(r8), intent(in)  :: ustar(ncol)
    real(r8), intent(out) :: khfs(ncol)
    real(r8), intent(out) :: kqfs(ncol)
    real(r8), intent(out) :: kbfs(ncol)
    real(r8), intent(out) :: obklen(ncol)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    call calc_obklen_point(g, vk, cpair, zvir, ths, thvs, qflx, shflx, &
         rrho, ustar, khfs, kqfs, kbfs, obklen)

  end subroutine calc_obklen_run

  elemental subroutine calc_obklen_point(g, vk, cpair, zvir, ths, thvs, &
       qflx, shflx, rrho, ustar, khfs, kqfs, kbfs, obklen)

    real(r8), intent(in)  :: g, vk, cpair, zvir
    real(r8), intent(in)  :: ths, thvs, qflx, shflx, rrho, ustar
    real(r8), intent(out) :: khfs, kqfs, kbfs, obklen

    khfs = shflx*rrho/cpair
    kqfs = qflx*rrho
    kbfs = khfs + zvir*ths*kqfs
    obklen = -thvs * ustar**3 / &
         (g*vk*(kbfs + sign(1.e-10_r8,kbfs)))

  end subroutine calc_obklen_point

end module ap_calc_obklen_scheme
