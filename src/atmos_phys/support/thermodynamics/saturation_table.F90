module ap_saturation_table

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

  public :: saturation_table_init
  public :: qsat_table
  public :: findsp_table_vc
  public :: register_qsat_table_hook
  public :: register_findsp_table_hook

  abstract interface
    subroutine qsat_scalar_hook(t, p, es, qs, gam, dqsdt, enthalpy)
      import :: r8
      real(r8), intent(in) :: t, p
      real(r8), intent(out) :: es, qs
      real(r8), intent(out), optional :: gam, dqsdt, enthalpy
    end subroutine qsat_scalar_hook

    subroutine findsp_vector_hook(q, t, p, use_ice, tsp, qsp)
      import :: r8
      real(r8), intent(in) :: q(:), t(:), p(:)
      logical, intent(in) :: use_ice
      real(r8), intent(out) :: tsp(:), qsp(:)
    end subroutine findsp_vector_hook
  end interface

  procedure(qsat_scalar_hook), pointer :: production_qsat => null()
  procedure(findsp_vector_hook), pointer :: production_findsp => null()

  interface qsat_table
    module procedure qsat_table_scalar
    module procedure qsat_table_1d
  end interface qsat_table

  real(r8), parameter :: tmin = 127.16_r8
  real(r8), parameter :: tmax = 375.16_r8
  real(r8), parameter :: ttrice = 20._r8
  real(r8), parameter :: pcf(5) = (/ &
       5.04469588506e-01_r8, &
      -5.47288442819e+00_r8, &
      -3.67471858735e-01_r8, &
      -8.95963532403e-03_r8, &
      -7.78053686625e-05_r8 /)

  real(r8), allocatable :: estbl(:)
  real(r8) :: epsilo
  real(r8) :: omeps
  real(r8) :: latvap
  real(r8) :: latice
  real(r8) :: rh2o
  real(r8) :: cpair
  real(r8) :: tmelt
  real(r8) :: c3

contains

  subroutine register_qsat_table_hook(qsat_proc)
    procedure(qsat_scalar_hook) :: qsat_proc

    production_qsat => qsat_proc
  end subroutine register_qsat_table_hook

  subroutine register_findsp_table_hook(findsp_proc)
    procedure(findsp_vector_hook) :: findsp_proc

    production_findsp => findsp_proc
  end subroutine register_findsp_table_hook

  subroutine saturation_table_init(table, epsilo_in, latvap_in, latice_in, &
                                   rh2o_in, cpair_in, tmelt_in)
    real(r8), intent(in) :: table(:)
    real(r8), intent(in) :: epsilo_in
    real(r8), intent(in) :: latvap_in
    real(r8), intent(in) :: latice_in
    real(r8), intent(in) :: rh2o_in
    real(r8), intent(in) :: cpair_in
    real(r8), intent(in) :: tmelt_in

    if (allocated(estbl)) deallocate(estbl)
    allocate(estbl(size(table)))
    estbl = table
    epsilo = epsilo_in
    omeps = 1._r8 - epsilo
    latvap = latvap_in
    latice = latice_in
    rh2o = rh2o_in
    cpair = cpair_in
    tmelt = tmelt_in
    c3 = 287.04_r8*(7.5_r8*log(10._r8))/cpair
  end subroutine saturation_table_init

  subroutine findsp_table_vc(q, t, p, use_ice, tsp, qsp, errflg)
    real(r8), intent(in) :: q(:)
    real(r8), intent(in) :: t(:)
    real(r8), intent(in) :: p(:)
    logical, intent(in) :: use_ice
    real(r8), intent(out) :: tsp(:)
    real(r8), intent(out) :: qsp(:)
    integer, intent(out), optional :: errflg

    integer :: i
    integer :: status
    integer :: local_errflg

    if (associated(production_findsp)) then
      call production_findsp(q, t, p, use_ice, tsp, qsp)
      local_errflg = 0
      if (present(errflg)) errflg = local_errflg
      return
    end if

    local_errflg = 0
    do i = 1, size(q)
      call findsp_table(q(i), t(i), p(i), use_ice, tsp(i), qsp(i), status)
      if (status == 2 .or. status == 8) local_errflg = 1
    end do
    if (present(errflg)) errflg = local_errflg
  end subroutine findsp_table_vc

  subroutine findsp_table(q, t, p, use_ice, tsp, qsp, status)
    real(r8), intent(in) :: q
    real(r8), intent(in) :: t
    real(r8), intent(in) :: p
    logical, intent(in) :: use_ice
    real(r8), intent(out) :: tsp
    real(r8), intent(out) :: qsp
    integer, intent(out) :: status

    integer, parameter :: iter = 8
    real(r8), parameter :: dttol = 1.e-4_r8
    real(r8), parameter :: dqtol = 1.e-4_r8
    integer :: l
    real(r8) :: es, gam, dgdt, gwork, hltalt, qs
    real(r8) :: t1, q1, dt, dq, qvd, r1b, c1, c2
    real(r8) :: enin, enout

    ! UW shallow convection always requests the mixed ice-water path.
    ! Retain the argument so this routine mirrors CAM's findsp interface.
    call qsat_table(t, p, es, qs)
    if (p <= 5._r8*es .or. qs <= 0._r8 .or. qs >= 0.5_r8 .or. &
        t < tmin .or. t > tmax) then
      status = 1
      tsp = t
      qsp = q
      return
    end if

    status = 2
    call calc_hltalt(t, hltalt, enout)
    enin = cpair*t + hltalt*q
    c1 = hltalt*c3
    c2 = (t+36._r8)**2
    r1b = c2/(c2+c1*qs)
    qvd = r1b*(q-qs)
    tsp = t+(hltalt/cpair)*qvd
    call qsat_table(tsp, p, es, qsp, gam=gam, enthalpy=enout)

    do l = 1, iter
      gwork = enin-enout
      dgdt = -cpair*(1._r8+gam)
      t1 = tsp-gwork/dgdt
      dt = abs(t1-tsp)/t1
      tsp = t1
      if (tsp < tmin) then
        tsp = tmin
        call calc_hltalt(tsp, hltalt, enout)
        qsp = (enin-cpair*tsp)/hltalt
        enout = cpair*tsp+hltalt*qsp
        status = 4
        exit
      end if
      call qsat_table(tsp, p, es, q1, gam=gam, enthalpy=enout)
      dq = abs(q1-qsp)/max(q1,1.e-12_r8)
      qsp = q1
      if (dt < dttol .and. dq < dqtol) then
        status = 0
        exit
      end if
    end do
    if (abs((enin-enout)/(enin+enout)) > 1.e-4_r8) status = 8
  end subroutine findsp_table

  subroutine qsat_table_scalar(t, p, es, qs, gam, dqsdt, enthalpy)
    real(r8), intent(in) :: t
    real(r8), intent(in) :: p
    real(r8), intent(out) :: es
    real(r8), intent(out) :: qs
    real(r8), intent(out), optional :: gam
    real(r8), intent(out), optional :: dqsdt
    real(r8), intent(out), optional :: enthalpy

    real(r8) :: hltalt
    real(r8) :: tterm

    if (associated(production_qsat)) then
      call production_qsat(t, p, es, qs, gam, dqsdt, enthalpy)
      return
    end if

    es = estbl_lookup(t)
    qs = svp_to_qsat(es, p)
    es = min(es, p)

    if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then
      call calc_hltalt(t, hltalt, tterm)
      if (present(enthalpy)) enthalpy = cpair*t + hltalt*qs
      call deriv_outputs(t, p, es, qs, hltalt, tterm, gam, dqsdt)
    end if
  end subroutine qsat_table_scalar

  subroutine qsat_table_1d(t, p, es, qs, gam, dqsdt, enthalpy)
    real(r8), intent(in) :: t(:), p(:)
    real(r8), intent(out) :: es(:), qs(:)
    real(r8), intent(out), optional :: gam(:), dqsdt(:), enthalpy(:)
    integer :: i

    do i = 1, size(t)
      if (present(gam)) then
        if (present(dqsdt)) then
          if (present(enthalpy)) then
            call qsat_table_scalar(t(i), p(i), es(i), qs(i), &
                 gam(i), dqsdt(i), enthalpy(i))
          else
            call qsat_table_scalar(t(i), p(i), es(i), qs(i), &
                 gam=gam(i), dqsdt=dqsdt(i))
          end if
        else if (present(enthalpy)) then
          call qsat_table_scalar(t(i), p(i), es(i), qs(i), &
               gam=gam(i), enthalpy=enthalpy(i))
        else
          call qsat_table_scalar(t(i), p(i), es(i), qs(i), gam=gam(i))
        end if
      else if (present(dqsdt)) then
        if (present(enthalpy)) then
          call qsat_table_scalar(t(i), p(i), es(i), qs(i), &
               dqsdt=dqsdt(i), enthalpy=enthalpy(i))
        else
          call qsat_table_scalar(t(i), p(i), es(i), qs(i), dqsdt=dqsdt(i))
        end if
      else if (present(enthalpy)) then
        call qsat_table_scalar(t(i), p(i), es(i), qs(i), &
             enthalpy=enthalpy(i))
      else
        call qsat_table_scalar(t(i), p(i), es(i), qs(i))
      end if
    end do
  end subroutine qsat_table_1d

  elemental function estbl_lookup(t) result(es)
    real(r8), intent(in) :: t
    real(r8) :: es
    real(r8) :: t_tmp
    real(r8) :: weight
    integer :: i

    t_tmp = max(min(t,tmax)-tmin, 0._r8)
    i = int(t_tmp) + 1
    weight = t_tmp - aint(t_tmp, r8)
    es = (1._r8-weight)*estbl(i) + weight*estbl(i+1)
  end function estbl_lookup

  elemental function svp_to_qsat(es, p) result(qs)
    real(r8), intent(in) :: es
    real(r8), intent(in) :: p
    real(r8) :: qs

    if ((p-es) <= 0._r8) then
      qs = 1._r8
    else
      qs = epsilo*es/(p-omeps*es)
    end if
  end function svp_to_qsat

  elemental subroutine no_ip_hltalt(t, hltalt)
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hltalt

    hltalt = latvap
    if (t >= tmelt) hltalt = hltalt - 2369._r8*(t-tmelt)
  end subroutine no_ip_hltalt

  elemental subroutine calc_hltalt(t, hltalt, tterm)
    real(r8), intent(in) :: t
    real(r8), intent(out) :: hltalt
    real(r8), intent(out) :: tterm

    real(r8) :: tc
    real(r8) :: weight
    integer :: i

    tterm = 0._r8
    call no_ip_hltalt(t, hltalt)
    if (t < tmelt) then
      tc = t-tmelt
      if (tc >= -ttrice) then
        weight = -tc/ttrice
        do i = size(pcf), 1, -1
          tterm = pcf(i) + tc*tterm
        end do
        tterm = tterm/ttrice
      else
        weight = 1._r8
      end if
      hltalt = hltalt + weight*latice
    end if
  end subroutine calc_hltalt

  elemental subroutine deriv_outputs(t, p, es, qs, hltalt, tterm, gam, dqsdt)
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

    if (qs == 1._r8) then
      dqsdt_loc = 0._r8
    else
      desdt = hltalt*es/(rh2o*t*t) + tterm
      dqsdt_loc = qs*p*desdt/(es*(p-omeps*es))
    end if
    if (present(dqsdt)) dqsdt = dqsdt_loc
    if (present(gam)) gam = dqsdt_loc*(hltalt/cpair)
  end subroutine deriv_outputs

end module ap_saturation_table
