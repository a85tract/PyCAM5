module ap_check_energy_chng_scheme
  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: check_energy_chng_timestep_init
  public :: check_energy_chng_run

contains

  !> \section arg_table_check_energy_chng_timestep_init Argument Table
  !! \htmlinclude check_energy_chng_timestep_init.html
  !!
  subroutine check_energy_chng_timestep_init( &
       ncol, pver, ixcldice, ixcldliq, ixrain, ixsnow, &
       gravit, latvap, latice, u, v, s, q, pdel, &
       te_ini, tw_ini, te_cur, tw_cur, tend_te_tnd, tend_tw_tnd, count, &
       errmsg, errflg)

    integer,  intent(in)    :: ncol, pver
    integer,  intent(in)    :: ixcldice, ixcldliq, ixrain, ixsnow
    real(r8), intent(in)    :: gravit, latvap, latice
    real(r8), intent(in)    :: u(:,:), v(:,:), s(:,:), q(:,:,:), pdel(:,:)
    real(r8), intent(out)   :: te_ini(:), tw_ini(:)
    real(r8), intent(out)   :: te_cur(:), tw_cur(:)
    real(r8), intent(out)   :: tend_te_tnd(:), tend_tw_tnd(:)
    integer,  intent(out)   :: count
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(r8) :: ke(ncol)
    real(r8) :: se(ncol)
    real(r8) :: wv(ncol)
    real(r8) :: wl(ncol)
    real(r8) :: wi(ncol)
    integer  :: i, k

    errmsg = ''
    errflg = 0

    ! This is the PI-atm algorithm from check_energy_timestep_init.  Preserve
    ! its loop nesting and expression order; do not replace it with the newer
    ! CAM hydrostatic-energy formulation.
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(u(i,k)**2 + v(i,k)**2)*pdel(i,k)/gravit
          se(i) = se(i) + s(i,k)*pdel(i,k)/gravit
          wv(i) = wv(i) + q(i,k,1)*pdel(i,k)/gravit
       end do
    end do

    if (ixcldliq > 1 .and. ixcldice > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + q(i,k,ixcldliq)*pdel(i,k)/gravit
             wi(i) = wi(i) + q(i,k,ixcldice)*pdel(i,k)/gravit
          end do
       end do
    end if

    if (ixrain > 1 .and. ixsnow > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + q(i,k,ixrain)*pdel(i,k)/gravit
             wi(i) = wi(i) + q(i,k,ixsnow)*pdel(i,k)/gravit
          end do
       end do
    end if

    do i = 1, ncol
       te_ini(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       tw_ini(i) = wv(i) + wl(i) + wi(i)

       te_cur(i) = te_ini(i)
       tw_cur(i) = tw_ini(i)
    end do

    tend_te_tnd(:ncol) = 0._r8
    tend_tw_tnd(:ncol) = 0._r8
    count = 0

  end subroutine check_energy_chng_timestep_init


  !> \section arg_table_check_energy_chng_run Argument Table
  !! \htmlinclude check_energy_chng_run.html
  !!
  subroutine check_energy_chng_run( &
       ncol, pver, ixcldice, ixcldliq, ixrain, ixsnow, &
       gravit, latvap, latice, u, v, s, q, pdel, &
       te_cur, tw_cur, tend_te_tnd, tend_tw_tnd, ztodt, &
       flx_vap, flx_cnd, flx_ice, flx_sen, &
       te, tw, te_xpd, tw_xpd, te_dif, tw_dif, te_tnd, tw_tnd, &
       te_rer, tw_rer, errmsg, errflg)

    integer,  intent(in)    :: ncol, pver
    integer,  intent(in)    :: ixcldice, ixcldliq, ixrain, ixsnow
    real(r8), intent(in)    :: gravit, latvap, latice
    real(r8), intent(in)    :: u(:,:), v(:,:), s(:,:), q(:,:,:), pdel(:,:)
    real(r8), intent(in)    :: te_cur(:), tw_cur(:)
    real(r8), intent(inout) :: tend_te_tnd(:), tend_tw_tnd(:)
    real(r8), intent(in)    :: ztodt
    real(r8), intent(in)    :: flx_vap(:), flx_cnd(:), flx_ice(:), flx_sen(:)
    real(r8), intent(out)   :: te(:), tw(:)
    real(r8), intent(out)   :: te_xpd(:), tw_xpd(:)
    real(r8), intent(out)   :: te_dif(:), tw_dif(:)
    real(r8), intent(out)   :: te_tnd(:), tw_tnd(:)
    real(r8), intent(out)   :: te_rer(:), tw_rer(:)
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    real(r8) :: ke(ncol)
    real(r8) :: se(ncol)
    real(r8) :: wv(ncol)
    real(r8) :: wl(ncol)
    real(r8) :: wi(ncol)
    integer  :: i, k

    errmsg = ''
    errflg = 0

    ! This is the PI-atm algorithm from check_energy_chng.  Keep the original
    ! vertical/column loop order and the original frozen-static-energy formula.
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(u(i,k)**2 + v(i,k)**2)*pdel(i,k)/gravit
          se(i) = se(i) + s(i,k)*pdel(i,k)/gravit
          wv(i) = wv(i) + q(i,k,1)*pdel(i,k)/gravit
       end do
    end do

    if (ixcldliq > 1 .and. ixcldice > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + q(i,k,ixcldliq)*pdel(i,k)/gravit
             wi(i) = wi(i) + q(i,k,ixcldice)*pdel(i,k)/gravit
          end do
       end do
    end if

    if (ixrain > 1 .and. ixsnow > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + q(i,k,ixrain)*pdel(i,k)/gravit
             wi(i) = wi(i) + q(i,k,ixsnow)*pdel(i,k)/gravit
          end do
       end do
    end if

    do i = 1, ncol
       te(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       tw(i) = wv(i) + wl(i) + wi(i)
    end do

    do i = 1, ncol
       te_dif(i) = te(i) - te_cur(i)
       tw_dif(i) = tw(i) - tw_cur(i)

       te_tnd(i) = flx_vap(i)*(latvap+latice) - &
            (flx_cnd(i) - flx_ice(i))*1000._r8*latice + flx_sen(i)
       tw_tnd(i) = flx_vap(i) - flx_cnd(i)*1000._r8

       tend_te_tnd(i) = tend_te_tnd(i) + te_tnd(i)
       tend_tw_tnd(i) = tend_tw_tnd(i) + tw_tnd(i)

       te_xpd(i) = te_cur(i) + te_tnd(i)*ztodt
       tw_xpd(i) = tw_cur(i) + tw_tnd(i)*ztodt

       te_rer(i) = (te_xpd(i) - te(i)) / te_cur(i)
    end do

    tw_rer = 0._r8
    where (tw_cur(:ncol) > 0._r8)
       tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / tw_cur(:ncol)
    end where

  end subroutine check_energy_chng_run

end module ap_check_energy_chng_scheme
