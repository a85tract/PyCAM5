module ap_pbl_radiative_feedback

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  public :: compute_radf

contains

subroutine compute_radf(choice_radf, i, pcols, pver, ncvmax, ncvfin, ktop, qmin, &
                        ql, pi, qrlw, g, cldeff, zi, chs, lwp_cl, opt_depth_cl, &
                        radinvfrac_cl, radf_cl)
  character(len=6), intent(in) :: choice_radf
  integer, intent(in) :: i
  integer, intent(in) :: pcols
  integer, intent(in) :: pver
  integer, intent(in) :: ncvmax
  integer, intent(in) :: ncvfin(pcols)
  integer, intent(in) :: ktop(pcols,ncvmax)
  real(r8), intent(in) :: qmin
  real(r8), intent(in) :: ql(pcols,pver)
  real(r8), intent(in) :: pi(pcols,pver+1)
  real(r8), intent(in) :: qrlw(pcols,pver)
  real(r8), intent(in) :: g
  real(r8), intent(in) :: cldeff(pcols,pver)
  real(r8), intent(in) :: zi(pcols,pver+1)
  real(r8), intent(in) :: chs(pcols,pver+1)
  real(r8), intent(out) :: lwp_cl(ncvmax)
  real(r8), intent(out) :: opt_depth_cl(ncvmax)
  real(r8), intent(out) :: radinvfrac_cl(ncvmax)
  real(r8), intent(out) :: radf_cl(ncvmax)

  integer :: kt, ncv
  real(r8) :: lwp, opt_depth, radinvfrac, radf

  lwp_cl = 0._r8
  opt_depth_cl = 0._r8
  radinvfrac_cl = 0._r8
  radf_cl = 0._r8

  do ncv = 1, ncvfin(i)
    kt = ktop(i,ncv)
    if (choice_radf == 'orig') then
      if (ql(i,kt) > qmin .and. ql(i,kt-1) < qmin) then
        lwp = ql(i,kt)*(pi(i,kt+1)-pi(i,kt))/g
        opt_depth = 156._r8*lwp
        radinvfrac = opt_depth*(4._r8+opt_depth)/ &
             (6._r8*(4._r8+opt_depth)+opt_depth**2)
        radf = qrlw(i,kt)/(pi(i,kt)-pi(i,kt+1))
        radf = max(radinvfrac*radf*(zi(i,kt)-zi(i,kt+1)),0._r8)*chs(i,kt)
      end if
    else if (choice_radf == 'ramp') then
      lwp = ql(i,kt)*(pi(i,kt+1)-pi(i,kt))/g
      opt_depth = 156._r8*lwp
      radinvfrac = opt_depth*(4._r8+opt_depth)/ &
           (6._r8*(4._r8+opt_depth)+opt_depth**2)
      radinvfrac = max(cldeff(i,kt)-cldeff(i,kt-1),0._r8)*radinvfrac
      radf = qrlw(i,kt)/(pi(i,kt)-pi(i,kt+1))
      radf = max(radinvfrac*radf*(zi(i,kt)-zi(i,kt+1)),0._r8)*chs(i,kt)
    else if (choice_radf == 'maxi') then
      lwp = ql(i,kt)*(pi(i,kt+1)-pi(i,kt))/g
      opt_depth = 156._r8*lwp
      radinvfrac = opt_depth*(4._r8+opt_depth)/ &
           (6._r8*(4._r8+opt_depth)+opt_depth**2)
      radf = max(radinvfrac*qrlw(i,kt)/(pi(i,kt)-pi(i,kt+1))* &
           (zi(i,kt)-zi(i,kt+1)),0._r8)
      lwp = ql(i,kt-1)*(pi(i,kt)-pi(i,kt-1))/g
      opt_depth = 156._r8*lwp
      radinvfrac = opt_depth*(4._r8+opt_depth)/ &
           (6._r8*(4._r8+opt_depth)+opt_depth**2)
      radf = radf + max(radinvfrac*qrlw(i,kt-1)/(pi(i,kt-1)-pi(i,kt))* &
           (zi(i,kt-1)-zi(i,kt)),0._r8)
      radf = max(radf,0._r8)*chs(i,kt)
    end if

    lwp_cl(ncv) = lwp
    opt_depth_cl(ncv) = opt_depth
    radinvfrac_cl(ncv) = radinvfrac
    radf_cl(ncv) = radf
  end do
end subroutine compute_radf

end module ap_pbl_radiative_feedback
