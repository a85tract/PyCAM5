module ap_wtrc_mass_fixer_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private
  public :: wtrc_mass_fixer_run

contains

  !> \section arg_table_wtrc_mass_fixer_run Argument Table
  !! \htmlinclude wtrc_mass_fixer_run.html
  !!
  subroutine wtrc_mass_fixer_run( &
       ncol, pver, pwtype, nwset, isphdo, wisotope, &
       iatype, bulk_indices, rstd_by_wset, wtrc_qmin, &
       limiter_hdo_hgh, limiter_hdo_low, limiter_18o_hgh, limiter_18o_low, &
       limiter_phis_crit, radtodeg, lat, phis, q, errmsg, errflg)

    integer, intent(in) :: ncol, pver, pwtype, nwset, isphdo
    logical, intent(in) :: wisotope
    integer, intent(in) :: iatype(:,:), bulk_indices(:)
    real(r8), intent(in) :: rstd_by_wset(:,:)
    real(r8), intent(in) :: wtrc_qmin
    real(r8), intent(in) :: limiter_hdo_hgh, limiter_hdo_low
    real(r8), intent(in) :: limiter_18o_hgh, limiter_18o_low
    real(r8), intent(in) :: limiter_phis_crit, radtodeg
    real(r8), intent(in) :: lat(:), phis(:)
    real(r8), intent(inout) :: q(:,:,:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer  :: i, k, p, m
    real(r8) :: oval
    real(r8) :: diff
    real(r8) :: ratio
    real(r8) :: wtlat, wtphis
    real(r8) :: diff_limit

    errmsg = ''
    errflg = 0

    ! PI-atm wtrc_mass_fixer, retaining the original i-k-p-m nesting and
    ! statement/expression order.
    do i = 1, ncol
       do k = 1, pver
          do p = 1, pwtype
             if (q(i,k,iatype(1,p)) .ne. q(i,k,bulk_indices(p))) then
                diff = q(i,k,iatype(1,p)) - q(i,k,bulk_indices(p))
                oval = q(i,k,iatype(1,p))
                do m = 1, nwset
                   ratio = ratio_from_mass(q(i,k,iatype(m,p)), oval, &
                        rstd_by_wset(m,p), wtrc_qmin)
                   q(i,k,iatype(m,p)) = q(i,k,iatype(m,p)) - ratio*diff
                end do
             end if
          end do
       end do
    end do

    if (wisotope) then
       do i = 1, ncol
          do k = 1, pver
             do p = 1, pwtype
                if (q(i,k,iatype(2,p)) .ne. q(i,k,bulk_indices(p))) then
                   diff = q(i,k,iatype(2,p)) - q(i,k,bulk_indices(p))
                   oval = q(i,k,iatype(2,p))
                   do m = 2, nwset
                      ratio = ratio_from_mass(q(i,k,iatype(m,p)), oval, &
                           rstd_by_wset(m,p), wtrc_qmin)
                      q(i,k,iatype(m,p)) = q(i,k,iatype(m,p)) - ratio*diff
                   end do
                end if
             end do
          end do
       end do
    end if

    do i = 1, ncol
       wtlat = lat(i)*radtodeg
       wtphis = phis(i)

       do k = 1, pver
          do p = 1, pwtype
             do m = 2, nwset
                if (m == isphdo) then
                   diff_limit = limiter_hdo_hgh
                else
                   diff_limit = limiter_18o_hgh
                end if
                if (q(i,k,iatype(m,p)) .gt. diff_limit*q(i,k,iatype(1,p))) then
                   q(i,k,iatype(m,p)) = q(i,k,iatype(1,p))
                else
                   if (abs(wtlat) < 60._r8 .and. wtphis > limiter_phis_crit) then
                      if (m == isphdo) then
                         diff_limit = limiter_hdo_low + 0.005_r8*(k - pver)
                      else
                         diff_limit = limiter_18o_low + 0.005_r8*(k - pver)
                      end if
                      if (q(i,k,iatype(m,p)) < diff_limit*q(i,k,iatype(1,p))) then
                         q(i,k,iatype(m,p)) = q(i,k,iatype(1,p))
                      end if
                   end if
                end if
             end do
          end do
       end do
    end do

  end subroutine wtrc_mass_fixer_run


  function ratio_from_mass(qtrc, qtot, rstd, qmin) result(ratio)
    real(r8), intent(in) :: qtrc, qtot, rstd, qmin
    real(r8) :: ratio

    if (abs(qtot) < qmin) then
       ratio = rstd
    else
       ratio = qtrc / qtot
    end if
  end function ratio_from_mass

end module ap_wtrc_mass_fixer_scheme
