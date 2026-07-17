module wv_saturation_portable

  ! Thin portable adapter around wv_sat_methods.  The saturation pressure
  ! formulae live only in wv_sat_methods; this module supplies the derivative
  ! outputs and water-only wet-bulb iteration used by Park macrophysics.

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: cpair  => shr_const_cpdair, &
       latvap => shr_const_latvap, latice => shr_const_latice
  use wv_sat_methods, only: wv_sat_methods_init, &
       wv_sat_svp_water, wv_sat_svp_ice, &
       wv_sat_qsat_water, wv_sat_qsat_ice

  implicit none
  private
  save

  real(r8), parameter :: tboil = 373.16_r8
  real(r8), parameter :: ttrice = 20.00_r8
  real(r8), parameter :: tmin = 127.16_r8
  real(r8), parameter :: tmax = 375.16_r8

  real(r8) :: epsilo
  real(r8) :: omeps
  real(r8) :: rh2o
  real(r8) :: tmelt
  real(r8) :: c3

  public :: wv_saturation_portable_init
  public :: svp_water, svp_ice
  public :: qsat_water, qsat_ice
  public :: findsp_vc

contains

  subroutine wv_saturation_portable_init(epsilo_in, latvap_in, latice_in, &
       rh2o_in, cpair_in, tmelt_in, h2otrip_in, errmsg, errflg)
    real(r8), intent(in) :: epsilo_in, latvap_in, latice_in
    real(r8), intent(in) :: rh2o_in, cpair_in, tmelt_in, h2otrip_in
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    epsilo = epsilo_in
    omeps = 1._r8 - epsilo
    rh2o = rh2o_in
    tmelt = tmelt_in
    c3 = 287.04_r8*(7.5_r8*log(10._r8))/cpair

    if (cpair_in /= cpair .or. latvap_in /= latvap .or. &
         latice_in /= latice) then
       errmsg = 'wv_saturation_portable_init: immutable host constants differ'
       errflg = 1
       return
    end if

    call wv_sat_methods_init(r8, tmelt_in, h2otrip_in, tboil, ttrice, &
         epsilo, errmsg)
    if (len_trim(errmsg) == 0) then
       errflg = 0
    else
       errflg = 1
    end if
  end subroutine wv_saturation_portable_init

  elemental function svp_water(t) result(es)
    real(r8), intent(in) :: t
    real(r8) :: es

    es = wv_sat_svp_water(t)
  end function svp_water

  elemental function svp_ice(t) result(es)
    real(r8), intent(in) :: t
    real(r8) :: es

    es = wv_sat_svp_ice(t)
  end function svp_ice

  elemental function tq_enthalpy(t, q, hltalt) result(enthalpy)
    real(r8), intent(in) :: t
    real(r8), intent(in) :: q
    real(r8), intent(in) :: hltalt
    real(r8) :: enthalpy

    enthalpy = cpair * t + hltalt * q
  end function tq_enthalpy

  elemental subroutine no_ip_hltalt(t, hltalt)
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hltalt

    hltalt = latvap
    if (t >= tmelt) then
       hltalt = hltalt - 2369.0_r8*(t-tmelt)
    end if
  end subroutine no_ip_hltalt

  elemental subroutine deriv_outputs(t, p, es, qs, hltalt, tterm, &
       gam, dqsdt)
    real(r8), intent(in) :: t
    real(r8), intent(in) :: p
    real(r8), intent(in) :: es
    real(r8), intent(in) :: qs
    real(r8), intent(in) :: hltalt
    real(r8), intent(in) :: tterm
    real(r8), intent(out), optional :: gam
    real(r8), intent(out), optional :: dqsdt

    real(r8) :: desdt
    real(r8) :: dqsdt_loc

    if (qs == 1.0_r8) then
       dqsdt_loc = 0._r8
    else
       desdt = hltalt*es/(rh2o*t*t) + tterm
       dqsdt_loc = qs*p*desdt/(es*(p-omeps*es))
    end if

    if (present(dqsdt)) dqsdt = dqsdt_loc
    if (present(gam))   gam   = dqsdt_loc * (hltalt/cpair)
  end subroutine deriv_outputs

  elemental subroutine qsat_water(t, p, es, qs, gam, dqsdt, enthalpy)
    real(r8), intent(in) :: t, p
    real(r8), intent(out) :: es, qs
    real(r8), intent(out), optional :: gam, dqsdt, enthalpy

    real(r8) :: hltalt

    call wv_sat_qsat_water(t, p, es, qs)

    if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then
       call no_ip_hltalt(t, hltalt)

       if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

       call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
            gam=gam, dqsdt=dqsdt)
    end if
  end subroutine qsat_water

  elemental subroutine qsat_ice(t, p, es, qs, gam, dqsdt, enthalpy)
    real(r8), intent(in) :: t, p
    real(r8), intent(out) :: es, qs
    real(r8), intent(out), optional :: gam, dqsdt, enthalpy

    real(r8) :: desdt, dqsdt_loc, hltalt

    call wv_sat_qsat_ice(t, p, es, qs)

    if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then
       hltalt = latvap + latice
       if (present(enthalpy)) enthalpy = cpair*t + hltalt*qs

       if (qs == 1._r8) then
          dqsdt_loc = 0._r8
       else
          desdt = hltalt*es/(rh2o*t*t)
          dqsdt_loc = qs*p*desdt/(es*(p-omeps*es))
       end if
       if (present(dqsdt)) dqsdt = dqsdt_loc
       if (present(gam)) gam = dqsdt_loc*(hltalt/cpair)
    end if
  end subroutine qsat_ice

  subroutine findsp_vc(q, t, p, use_ice, tsp, qsp)
    use cloud_microphysics_host_hooks, only: cloud_microphysics_endrun

    real(r8), intent(in) :: q(:), t(:), p(:)
    logical, intent(in) :: use_ice
    real(r8), intent(out) :: tsp(:), qsp(:)

    integer :: status(size(q))
    integer :: i

    call findsp(q, t, p, use_ice, tsp, qsp, status)

    do i = 1, size(q)
       if (status(i) == 2) then
          call cloud_microphysics_endrun( &
               'wv_saturation::FINDSP -- not converging')
       else if (status(i) == 8) then
          call cloud_microphysics_endrun( &
               'wv_saturation::FINDSP -- enthalpy is not conserved')
       end if
    end do
  end subroutine findsp_vc

  elemental subroutine findsp(q, t, p, use_ice, tsp, qsp, status)
    real(r8), intent(in) :: q, t, p
    logical, intent(in) :: use_ice
    real(r8), intent(out) :: tsp, qsp
    integer, intent(out) :: status

    integer, parameter :: iter = 8
    real(r8), parameter :: dttol = 1.e-4_r8
    real(r8), parameter :: dqtol = 1.e-4_r8
    integer :: l
    real(r8) :: es, gam, dgdt, g, hltalt, qs
    real(r8) :: t1, q1, dt, dq, qvd, r1b, c1, c2
    real(r8) :: enin, enout

    if (use_ice) then
       call qsat_ice(t, p, es, qs)
    else
       call qsat_water(t, p, es, qs)
    end if

    if (p <= 5._r8*es .or. qs <= 0._r8 .or. qs >= 0.5_r8 .or. &
         t < tmin .or. t > tmax) then
       status = 1
       tsp = t
       qsp = q
       enin = 1._r8
       enout = 1._r8
       return
    end if

    status = 2
    if (use_ice) then
       hltalt = latvap + latice
    else
       call no_ip_hltalt(t, hltalt)
    end if
    enin = tq_enthalpy(t, q, hltalt)

    c1 = hltalt*c3
    c2 = (t + 36._r8)**2
    r1b = c2/(c2 + c1*qs)
    qvd = r1b * (q - qs)
    tsp = t + ((hltalt/cpair)*qvd)

    if (use_ice) then
       call qsat_ice(tsp, p, es, qsp, gam=gam, enthalpy=enout)
    else
       call qsat_water(tsp, p, es, qsp, gam=gam, enthalpy=enout)
    end if

    do l = 1, iter
       g = enin - enout
       dgdt = -cpair * (1 + gam)
       t1 = tsp - g/dgdt
       dt = abs(t1 - tsp)/t1
       tsp = t1

       if (tsp < tmin) then
          tsp = tmin
          if (use_ice) then
             hltalt = latvap + latice
          else
             call no_ip_hltalt(tsp, hltalt)
          end if
          qsp = (enin - cpair*tsp)/hltalt
          enout = tq_enthalpy(tsp, qsp, hltalt)
          status = 4
          exit
       end if

       if (use_ice) then
          call qsat_ice(tsp, p, es, q1, gam=gam, enthalpy=enout)
       else
          call qsat_water(tsp, p, es, q1, gam=gam, enthalpy=enout)
       end if
       dq = abs(q1 - qsp)/max(q1,1.e-12_r8)
       qsp = q1

       if (dt < dttol .and. dq < dqtol) then
          status = 0
          exit
       end if
    end do

    if (abs((enin-enout)/(enin+enout)) > 1.e-4_r8) status = 8
  end subroutine findsp

end module wv_saturation_portable
