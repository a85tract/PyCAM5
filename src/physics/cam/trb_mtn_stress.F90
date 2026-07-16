module trb_mtn_stress

  use ap_compute_tms_scheme, only : compute_tms_run
  use perf_mod,              only : t_startf, t_stopf

  implicit none
  private      
  save

  public init_tms                             ! Initialization
  public compute_tms                          ! Full routine

  ! ------------ !
  ! Private data !
  ! ------------ !

  integer,  parameter :: r8 = selected_real_kind(12) ! 8 byte real

  real(r8)            :: orocnst              ! Converts from standard deviation to height [ no unit ]
  real(r8)            :: z0fac                ! Factor determining z_0 from orographic standard deviation [ no unit ] 
  real(r8)            :: karman               ! von Karman constant
  real(r8)            :: gravit               ! Acceleration due to gravity
  real(r8)            :: rair                 ! Gas constant for dry air

contains

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_tms( kind, oro_in, z0fac_in, karman_in, gravit_in, rair_in, &
       errstring)

    integer, intent(in) :: kind

    real(r8), intent(in) :: oro_in, z0fac_in, karman_in, gravit_in, rair_in

    character(len=*), intent(out) :: errstring

    errstring = ' '

    if ( kind /= r8 ) then
       errstring = 'inconsistent KIND of reals passed to init_tms'
       return
    endif

    orocnst  = oro_in
    z0fac    = z0fac_in
    karman   = karman_in
    gravit   = gravit_in
    rair     = rair_in
    
  end subroutine init_tms

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine compute_tms( pcols    , pver    , ncol    ,                     &
                          u        , v       , t       , pmid    , exner   , &
                          zm       , sgh     , ksrf    , taux    , tauy    , & 
                          landfrac )

    !------------------------------------------------------------------------------ !
    ! Turbulent mountain stress parameterization                                    !  
    !                                                                               !
    ! Returns surface drag coefficient and stress associated with subgrid mountains !
    ! For points where the orographic variance is small ( including ocean ),        !
    ! the returned surface drag coefficient and stress is zero.                     !
    !                                                                               !
    ! Lastly arranged : Sungsu Park. Jan. 2010.                                     !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments ! 
    ! ---------------------- !

    integer,  intent(in)  :: pcols                 ! Number of columns dimensioned
    integer,  intent(in)  :: pver                  ! Number of model layers
    integer,  intent(in)  :: ncol                  ! Number of columns actually used

    real(r8), intent(in)  :: u(pcols,pver)         ! Layer mid-point zonal wind [ m/s ]
    real(r8), intent(in)  :: v(pcols,pver)         ! Layer mid-point meridional wind [ m/s ]
    real(r8), intent(in)  :: t(pcols,pver)         ! Layer mid-point temperature [ K ]
    real(r8), intent(in)  :: pmid(pcols,pver)      ! Layer mid-point pressure [ Pa ]
    real(r8), intent(in)  :: exner(pcols,pver)     ! Layer mid-point exner function [ no unit ]
    real(r8), intent(in)  :: zm(pcols,pver)        ! Layer mid-point height [ m ]
    real(r8), intent(in)  :: sgh(pcols)            ! Standard deviation of orography [ m ]
    real(r8), intent(in)  :: landfrac(pcols)       ! Land fraction [ fraction ]
    
    real(r8), intent(out) :: ksrf(pcols)           ! Surface drag coefficient [ kg/s/m2 ]
    real(r8), intent(out) :: taux(pcols)           ! Surface zonal      wind stress [ N/m2 ]
    real(r8), intent(out) :: tauy(pcols)           ! Surface meridional wind stress [ N/m2 ]

    character(len=512) :: errmsg
    integer            :: errflg

    call t_startf('ap_compute_tms_run')
    call compute_tms_run(pcols, pver, ncol, orocnst, z0fac, karman, &
         gravit, rair, u, v, t, pmid, exner, zm, sgh, landfrac, &
         ksrf, taux, tauy, errmsg, errflg)
    call t_stopf('ap_compute_tms_run')

  end subroutine compute_tms

end module trb_mtn_stress
