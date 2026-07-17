module wv_saturation_portable

  ! Portable adapter around wv_sat_methods.  The saturation pressure formulae
  ! live only in the standalone wv_sat_methods copy; this module supplies the
  ! mixed-phase table, derivative outputs, and wet-bulb iteration used by Park
  ! macrophysics.

  use shr_kind_mod, only: r8 => shr_kind_r8
  use shr_const_mod, only: cpair  => shr_const_cpdair, &
       latvap => shr_const_latvap, latice => shr_const_latice
  use wv_sat_methods, only: wv_sat_methods_init, &
       wv_sat_svp_water, wv_sat_svp_ice, &
       wv_sat_svp_trans, wv_sat_svp_to_qsat, &
       wv_sat_qsat_water, wv_sat_qsat_ice

  implicit none
  private
  save

  real(r8), parameter :: tboil = 373.16_r8
  real(r8), parameter :: ttrice = 20.00_r8
  real(r8), parameter :: tmin = 127.16_r8
  real(r8), parameter :: tmax = 375.16_r8

  integer :: plenest
  real(r8), allocatable :: estbl(:)

  real(r8) :: epsilo
  real(r8) :: omeps
  real(r8) :: rh2o
  real(r8) :: tmelt
  real(r8) :: c3
  integer :: iulog

  real(r8) :: pcf(5) = (/ &
       5.04469588506e-01_r8, &
       -5.47288442819e+00_r8, &
       -3.67471858735e-01_r8, &
       -8.95963532403e-03_r8, &
       -7.78053686625e-05_r8 /)

  public :: wv_saturation_portable_init
  public :: svp_water, svp_ice
  public :: qsat, qsat_water, qsat_ice
  public :: findsp_vc

contains

  subroutine wv_saturation_portable_init(epsilo_in, latvap_in, latice_in, &
       rh2o_in, cpair_in, tmelt_in, h2otrip_in, iulog_in, errmsg, errflg)
    real(r8), intent(in) :: epsilo_in, latvap_in, latice_in
    real(r8), intent(in) :: rh2o_in, cpair_in, tmelt_in, h2otrip_in
    integer, intent(in) :: iulog_in
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: status
    real(r8) :: t
    integer :: i

    epsilo = epsilo_in
    omeps = 1._r8 - epsilo
    rh2o = rh2o_in
    tmelt = tmelt_in
    iulog = iulog_in
    c3 = 287.04_r8*(7.5_r8*log(10._r8))/cpair

    if (cpair_in /= cpair .or. latvap_in /= latvap .or. &
         latice_in /= latice) then
       errmsg = 'wv_saturation_portable_init: immutable host constants differ'
       errflg = 1
       return
    end if

    call wv_sat_methods_init(r8, tmelt_in, h2otrip_in, tboil, ttrice, &
         epsilo, errmsg)
    if (len_trim(errmsg) /= 0) then
       errflg = 1
       return
    end if

    plenest = ceiling(tmax-tmin) + 2

    allocate(estbl(plenest), stat=status)
    if (status /= 0) then
       errmsg = 'wv_saturation_portable_init: ERROR allocating saturation vapor pressure table'
       errflg = 1
       return
    end if

    do i = 1, plenest
       estbl(i) = svp_trans(tmin + real(i-1,r8))
    end do

    errflg = 0
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

  elemental function svp_trans(t) result(es)
    real(r8), intent(in) :: t
    real(r8) :: es

    es = wv_sat_svp_trans(t)
  end function svp_trans

  elemental function estblf(t) result(es)
    real(r8), intent(in) :: t
    real(r8) :: es

    integer :: i
    real(r8) :: t_tmp
    real(r8) :: weight

    t_tmp = max(min(t,tmax)-tmin, 0._r8)
    i = int(t_tmp) + 1
    weight = t_tmp - aint(t_tmp, r8)
    es = (1._r8 - weight)*estbl(i) + weight*estbl(i+1)
  end function estblf

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

  elemental subroutine calc_hltalt(t, hltalt, tterm)
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hltalt
    real(r8), intent(out), optional :: tterm

    real(r8) :: tc
    real(r8) :: weight
    integer :: i

    if (present(tterm)) tterm = 0.0_r8

    call no_ip_hltalt(t,hltalt)
    if (t < tmelt) then
       tc = t - tmelt

       if (tc >= -ttrice) then
          weight = -tc/ttrice

          if (present(tterm)) then
             do i = size(pcf), 1, -1
                tterm = pcf(i) + tc*tterm
             end do
             tterm = tterm/ttrice
          end if

       else
          weight = 1.0_r8
       end if

       hltalt = hltalt + weight*latice

    end if
  end subroutine calc_hltalt

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

  elemental subroutine qsat(t, p, es, qs, gam, dqsdt, enthalpy)
    real(r8), intent(in) :: t
    real(r8), intent(in) :: p
    real(r8), intent(out) :: es
    real(r8), intent(out) :: qs
    real(r8), intent(out), optional :: gam
    real(r8), intent(out), optional :: dqsdt
    real(r8), intent(out), optional :: enthalpy

    real(r8) :: hltalt
    real(r8) :: tterm

    es = estblf(t)

    qs = wv_sat_svp_to_qsat(es, p)

    es = min(es, p)

    if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

       call calc_hltalt(t, hltalt, tterm)

       if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

       call deriv_outputs(t, p, es, qs, hltalt, tterm, &
            gam=gam, dqsdt=dqsdt)

    end if
  end subroutine qsat

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

    real(r8) :: hltalt

    call wv_sat_qsat_ice(t, p, es, qs)

    if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then
       hltalt = latvap + latice

       if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

       call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
            gam=gam, dqsdt=dqsdt)
    end if
  end subroutine qsat_ice

  subroutine findsp_vc(q, t, p, use_ice, tsp, qsp)
    use cloud_microphysics_host_hooks, only: endrun => cloud_microphysics_endrun

    real(r8), intent(in) :: q(:), t(:), p(:)
    logical, intent(in) :: use_ice
    real(r8), intent(out) :: tsp(:), qsp(:)

    integer :: status(size(q))
    integer :: n, i

    n = size(q)

    call findsp(q, t, p, use_ice, tsp, qsp, status)

    do i = 1,n
       if (status(i) == 2) then
          write(iulog,*) ' findsp not converging at i = ', i
          write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
          write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
          call endrun ('wv_saturation::FINDSP -- not converging')
       else if (status(i) == 8) then
          write(iulog,*) ' the enthalpy is not conserved at i = ', i
          write(iulog,*) ' t, q, p ', t(i), q(i), p(i)
          write(iulog,*) ' tsp, qsp ', tsp(i), qsp(i)
          call endrun ('wv_saturation::FINDSP -- enthalpy is not conserved')
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
       call qsat(t, p, es, qs)
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
       call calc_hltalt(t,hltalt)
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
       call qsat(tsp, p, es, qsp, gam=gam, enthalpy=enout)
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
             call calc_hltalt(tsp,hltalt)
          else
             call no_ip_hltalt(tsp, hltalt)
          end if
          qsp = (enin - cpair*tsp)/hltalt
          enout = tq_enthalpy(tsp, qsp, hltalt)
          status = 4
          exit
       end if

       if (use_ice) then
          call qsat(tsp, p, es, q1, gam=gam, enthalpy=enout)
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
