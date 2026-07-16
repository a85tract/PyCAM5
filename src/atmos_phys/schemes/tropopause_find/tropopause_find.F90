! This module is used to diagnose the location of the tropopause. Multiple
! algorithms are provided, some of which may not be able to identify a
! tropopause in all situations. To handle these cases, an analytic
! definition and a climatology are provided that can be used to fill in
! when the original algorithm fails. The tropopause temperature and
! pressure are determined and can be output to the history file.
!
! Author: Charles Bardeen
! Created: April, 2009
!
! CCPP-ized: Haipeng Lin, August 2024
! Adapted to preserve the PI-atm single-request API: July 2026
module tropopause_find_scheme
  !---------------------------------------------------------------
  ! ... variables for the tropopause module
  !---------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none

  private

  ! CCPP-compliant subroutines
  public :: tropopause_find_init
  public :: tropopause_find_run

  ! Switches for tropopause method
  public :: TROP_ALG_NONE, TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  public :: TROP_ALG_STOBIE, TROP_ALG_HYBSTOB, TROP_ALG_TWMO, TROP_ALG_WMO
  public :: TROP_ALG_CPP
  public :: NOTFOUND

  save

  ! These parameters define an enumeration to be used to define the primary
  ! and backup algorithms to be used with the tropopause_find() method. The
  ! backup algorithm is meant to provide a solution when the primary algorithm
  ! fails. The algorithms that can't fail (i.e., always find a tropopause) are:
  ! TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE, and TROP_ALG_STOBIE.
  integer, parameter    :: TROP_ALG_NONE      = 1    ! Don't evaluate
  integer, parameter    :: TROP_ALG_ANALYTIC  = 2    ! Analytic Expression
  integer, parameter    :: TROP_ALG_CLIMATE   = 3    ! Climatology
  integer, parameter    :: TROP_ALG_STOBIE    = 4    ! Stobie Algorithm
  integer, parameter    :: TROP_ALG_TWMO      = 5    ! WMO Definition, Reichler et al. [2003]
  integer, parameter    :: TROP_ALG_WMO       = 6    ! WMO Definition
  integer, parameter    :: TROP_ALG_HYBSTOB   = 7    ! Hybrid Stobie Algorithm
  integer, parameter    :: TROP_ALG_CPP       = 8    ! Cold Point Parabolic
  integer, parameter    :: TROP_ALG_CHEMTROP  = 9    ! Chemical tropopause

  integer, parameter    :: default_primary    = TROP_ALG_TWMO        ! default primary algorithm
  integer, parameter    :: default_backup     = TROP_ALG_CLIMATE     ! default backup algorithm

  integer, parameter    :: NOTFOUND = -1

  real(r8), parameter :: ALPHA  = 0.03_r8

  ! physical constants
  ! These constants are set in module variables rather than as parameters
  ! to support the aquaplanet mode in which the constants have values determined
  ! by the experiment protocol
  real(r8) :: cnst_kap     ! = cappa
  real(r8) :: cnst_faktor  ! = -gravit/rair
  real(r8) :: cnst_rga     ! = 1/gravit
  real(r8) :: cnst_ka1     ! = cnst_kap - 1._r8
  real(r8) :: cnst_rad2deg ! = 180/pi

!================================================================================================
contains
!================================================================================================

!> \section arg_table_tropopause_find_init Argument Table
!! \htmlinclude tropopause_find_init.html
  subroutine tropopause_find_init(cappa, rair, gravit, pi, errmsg, errflg)

    real(r8), intent(in)         :: cappa   ! R/Cp
    real(r8), intent(in)         :: rair    ! Dry air gas constant (J K-1 kg-1)
    real(r8), intent(in)         :: gravit  ! Gravitational acceleration (m s-2)
    real(r8), intent(in)         :: pi      ! Pi

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! define physical constants
    cnst_kap     = cappa
    cnst_faktor  = -gravit/rair
    cnst_rga     = 1._r8/gravit              ! Reciprocal of gravit (s2 m-1)
    cnst_ka1     = cnst_kap - 1._r8
    cnst_rad2deg = 180._r8/pi               ! radians to degrees conversion factor

  end subroutine tropopause_find_init


!> \section arg_table_tropopause_find_run Argument Table
!! \htmlinclude tropopause_find_run.html
  ! Searches all the columns and attempts to identify the tropopause.
  ! Two routines can be specifed, a primary routine which is tried first and a
  ! backup routine which will be tried only if the first routine fails. If the
  ! tropopause can not be identified by either routine, then a NOTFOUND is returned
  ! for the tropopause level, temperature and pressure.
  subroutine tropopause_find_run(ncol, pver, fillvalue, lat, pint, pmid, t, zi, zm, phis, &
                                       calday, tropp_p_loc, tropp_days, &
                                       tropLev, tropP, tropT, tropZ, &
                                       hstobie_trop, hstobie_linoz, hstobie_tropop, &
                                       primary, backup, &
                                       errmsg, errflg)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: fillvalue     ! Fill value for diagnostic outputs
    real(r8), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    real(r8), intent(in)         :: calday        ! Day of year including fraction of day

    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    real(r8), intent(in)         :: tropp_p_loc(:,:)
    real(r8), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12

    integer,                   intent(out)     :: tropLev(:)          ! tropopause level index
    real(r8), optional, intent(out)     :: tropP(:)            ! tropopause pressure (Pa)
    real(r8), optional, intent(out)     :: tropT(:)            ! tropopause temperature (K)
    real(r8), optional, intent(out)     :: tropZ(:)            ! tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(r8), optional, intent(inout)   :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(r8), optional, intent(inout)   :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(r8), optional, intent(inout)   :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    ! primary and backup are no longer optional arguments for CCPP-compliance.
    ! specify defaults when calling (TWMO, CLIMO)
    integer, intent(in)       :: primary                   ! primary detection algorithm
    integer, intent(in)       :: backup                    ! backup detection algorithm

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! Initialize the results to a missing value, so that the algorithms will
    ! attempt to find the tropopause for all of them.
    tropLev(:) = NOTFOUND
    if(present(tropP)) tropP(:) = fillvalue
    if(present(tropT)) tropT(:) = fillvalue
    if(present(tropZ)) tropZ(:) = fillvalue

    ! Try to find the tropopause using the primary algorithm.
    if (primary /= TROP_ALG_NONE) then
      call tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                calday, tropp_p_loc, tropp_days, &
                                primary, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, &
                                hstobie_trop=hstobie_trop, hstobie_linoz=hstobie_linoz, hstobie_tropop=hstobie_tropop, & ! only for HYBSTOB
                                errmsg=errmsg, errflg=errflg)
    end if

    if ((backup /= TROP_ALG_NONE) .and. any(tropLev(:) == NOTFOUND)) then
      call tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                calday, tropp_p_loc, tropp_days, &
                                backup, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, &
                                hstobie_trop=hstobie_trop, hstobie_linoz=hstobie_linoz, hstobie_tropop=hstobie_tropop, &
                                errmsg=errmsg, errflg=errflg)
    end if

  end subroutine tropopause_find_run

  ! Call the appropriate tropopause detection routine based upon the algorithm
  ! specifed.
  !
  ! NOTE: It is assumed that the output fields have been initialized by the
  ! caller, and only output values set to fillvalue will be detected.
  subroutine tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                  calday, tropp_p_loc, tropp_days, &
                                  algorithm, tropLev, tropP, tropT, tropZ, &
                                  hstobie_trop, hstobie_linoz, hstobie_tropop, &
                                  errmsg, errflg)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    real(r8), intent(in)         :: calday        ! Day of year including fraction of day

    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    real(r8), intent(in)         :: tropp_p_loc(:,:)
    real(r8), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12

    integer,                   intent(in)      :: algorithm             ! detection algorithm
    integer,                   intent(inout)   :: tropLev(:)            ! tropopause level index
    real(r8), optional, intent(inout)   :: tropP(:)              ! tropopause pressure (Pa)
    real(r8), optional, intent(inout)   :: tropT(:)              ! tropopause temperature (K)
    real(r8), optional, intent(inout)   :: tropZ(:)              ! tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(r8), optional, intent(inout)   :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(r8), optional, intent(inout)   :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(r8), optional, intent(inout)   :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! Dispatch the request to the appropriate routine.
    select case(algorithm)
      case(TROP_ALG_ANALYTIC)
        call tropopause_analytic(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_CLIMATE)
        call tropopause_climate(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                calday, tropp_p_loc, tropp_days, &
                                tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_STOBIE)
        call tropopause_stobie(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_HYBSTOB)
        if(present(hstobie_trop) .and. present(hstobie_linoz) .and. present(hstobie_tropop)) then
          call tropopause_hybridstobie(ncol, pver, pmid, t, zm, &
                                       tropLev, tropP, tropT, tropZ, &
                                       hstobie_trop, hstobie_linoz, hstobie_tropop)
        else
          call tropopause_hybridstobie(ncol, pver, pmid, t, zm, &
                                       tropLev, tropP, tropT, tropZ)
        endif

      case(TROP_ALG_TWMO)
        call tropopause_twmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_WMO)
        call tropopause_wmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_CPP)
        call tropopause_cpp(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_CHEMTROP)
        ! hplin: needs climatological arguments as calling tropopause_climate from within findChemTrop
        call tropopause_findChemTrop(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                     calday, tropp_p_loc, tropp_days, &
                                     tropLev, tropP, tropT, tropZ, &
                                     errmsg, errflg)

      case default
        errflg = 1
        write(errmsg,*) 'tropopause: Invalid detection algorithm (',  algorithm, ') specified.'
    end select

  end subroutine tropopause_findUsing

  ! This analytic expression closely matches the mean tropopause determined
  ! by the NCEP reanalysis and has been used by the radiation code.
  subroutine tropopause_analytic(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                 tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,                   intent(inout)   :: tropLev(:)             ! tropopause level index
    real(r8), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(r8), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(r8), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)

    ! Local Variables
    integer       :: i
    integer       :: k
    real(r8)      :: tP                       ! tropopause pressure (Pa)

    ! Iterate over all of the columns.
    do i = 1, ncol

      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then

        ! Calculate the pressure of the tropopause.
        tP = (25000.0_r8 - 15000.0_r8 * (cos(lat(i)))**2)

        ! Find the level that contains the tropopause.
        do k = pver, 2, -1
          if (tP >= pint(i, k)) then
            tropLev(i) = k
            exit
          end if
        end do

        ! Return the optional outputs
        if (present(tropP)) tropP(i) = tP

        if (present(tropT)) then
          tropT(i) = tropopause_interpolateT(pver, pmid, t, i, tropLev(i), tP)
        end if

        if (present(tropZ)) then
          tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
        end if
      end if
    end do

  end subroutine tropopause_analytic

  ! Determine the tropopause pressure from a climatology,
  ! interpolated to the current day of year and latitude.
  subroutine tropopause_climate(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                calday, tropp_p_loc, tropp_days, tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    real(r8), intent(in)         :: calday        ! Day of year including fraction of day

    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    real(r8), intent(in)         :: tropp_p_loc(:,:)
    real(r8), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12

    integer,                   intent(inout)  :: tropLev(:)            ! tropopause level index
    real(r8), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(r8), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(r8), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables
    integer       :: i
    integer       :: k
    integer       :: m
    real(r8)      :: tP                       ! tropopause pressure (Pa)
    real(r8)      :: dels
    integer       :: last
    integer       :: next

    ! If any columns remain to be indentified, then get the current
    ! day from the calendar.

    if (any(tropLev == NOTFOUND)) then

      !--------------------------------------------------------
      ! ... setup the time interpolation
      !--------------------------------------------------------
      if( calday < tropp_days(1) ) then
        next = 1
        last = 12
        dels = (365._r8 + calday - tropp_days(12)) / (365._r8 + tropp_days(1) - tropp_days(12))
      else if( calday >= tropp_days(12) ) then
        next = 1
        last = 12
        dels = (calday - tropp_days(12)) / (365._r8 + tropp_days(1) - tropp_days(12))
      else
        do m = 11,1,-1
           if( calday >= tropp_days(m) ) then
              exit
           end if
        end do
        last = m
        next = m + 1
        dels = (calday - tropp_days(m)) / (tropp_days(m+1) - tropp_days(m))
      end if

      dels = max( min( 1._r8,dels ),0._r8 )


      ! Iterate over all of the columns.
      do i = 1, ncol

        ! Skip column in which the tropopause has already been found.
        if (tropLev(i) == NOTFOUND) then

        !--------------------------------------------------------
        ! ... get tropopause level from climatology
        !--------------------------------------------------------
          ! Interpolate the tropopause pressure.
          tP = tropp_p_loc(i,last) &
            + dels * (tropp_p_loc(i,next) - tropp_p_loc(i,last))

          ! Find the associated level.
          do k = pver, 2, -1
            if (tP >= pint(i, k)) then
              tropLev(i) = k
              exit
            end if
          end do

          ! Return the optional outputs
          if (present(tropP)) tropP(i) = tP

          if (present(tropT)) then
            tropT(i) = tropopause_interpolateT(pver, pmid, t, i, tropLev(i), tP)
          end if

          if (present(tropZ)) then
            tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
          end if
        end if
      end do
    end if

  end subroutine tropopause_climate

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine tropopause_hybridstobie(ncol, pver, pmid, t, zm, &
                                     tropLev, tropP, tropT, tropZ, &
                                     hstobie_trop, hstobie_linoz, hstobie_tropop)

    !-----------------------------------------------------------------------
    ! Originally written by Philip Cameron-Smith, LLNL
    !
    !   Stobie-Linoz hybrid: the highest altitude of
    !          a) Stobie algorithm, or
    !          b) minimum Linoz pressure.
    !
    ! NOTE: the ltrop(i) gridbox itself is assumed to be a STRATOSPHERIC gridbox.
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !        ... Local variables
    !-----------------------------------------------------------------------
    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levelserp
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)

    integer,            intent(inout)   :: tropLev(:)             ! tropopause level index
    real(r8), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(r8), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(r8), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(r8), optional, intent(inout)   :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(r8), optional, intent(inout)   :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(r8), optional, intent(inout)   :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    real(r8),parameter  ::  min_Stobie_Pressure= 40.E2_r8 !For case 2 & 4.  [Pa]
    real(r8),parameter  ::  max_Linoz_Pressure =208.E2_r8 !For case     4.  [Pa]

    integer      :: i, k
    real(r8)     :: stobie_min, shybrid_temp      !temporary variable for case 2 & 3.
    integer      :: ltrop_linoz(ncol)            !Lowest possible Linoz vertical level
    integer      :: ltrop_trop(ncol)             !Tropopause level for hybrid case.
    logical      :: ltrop_linoz_set               !Flag that lowest linoz level already found.
    real(r8)     :: trop_output(ncol,pver)        !For output purposes only.
    real(r8)     :: trop_linoz_output(ncol,pver)  !For output purposes only.
    real(r8)     :: trop_trop_output(ncol,pver)   !For output purposes only.

    ltrop_linoz(:) = 1  ! Initialize to default value.
    ltrop_trop(:) = 1   ! Initialize to default value.

    LOOP_COL4: do i=1,ncol

       ! Skip column in which the tropopause has already been found.
       not_found: if (tropLev(i) == NOTFOUND) then

          stobie_min = 1.e10_r8    ! An impossibly large number
          ltrop_linoz_set = .FALSE.
          LOOP_LEV: do k=pver,1,-1
             IF (pmid(i,k) < min_stobie_pressure) cycle
             shybrid_temp = ALPHA * t(i,k) - Log10(pmid(i,k))
             !PJC_NOTE: the units of pmid won't matter, because it is just an additive offset.
             IF (shybrid_temp<stobie_min) then
                ltrop_trop(i)=k
                stobie_min = shybrid_temp
             ENDIF
             IF (pmid(i,k) < max_Linoz_pressure .AND. .NOT. ltrop_linoz_set) THEN
                ltrop_linoz(i) = k
                ltrop_linoz_set = .TRUE.
             ENDIF
          enddo LOOP_LEV

          tropLev(i) = MIN(ltrop_trop(i),ltrop_linoz(i))

          if (present(tropP)) then
             tropP(i) = pmid(i,tropLev(i))
          endif
          if (present(tropT)) then
             tropT(i) = t(i,tropLev(i))
          endif
          if (present(tropZ)) then
             tropZ(i) = zm(i,tropLev(i))
          endif

       endif not_found

    enddo LOOP_COL4

    trop_output(:,:)=0._r8
    trop_linoz_output(:,:)=0._r8
    trop_trop_output(:,:)=0._r8
    do i=1,ncol
       if (tropLev(i)>0) then
          trop_output(i,tropLev(i))=1._r8
          trop_linoz_output(i,ltrop_linoz(i))=1._r8
          trop_trop_output(i,ltrop_trop(i))=1._r8
       endif
    enddo

    if(present(hstobie_trop)) then
      hstobie_trop(:,:) = trop_output(:,:)
    endif

    if(present(hstobie_linoz)) then
      hstobie_linoz(:,:) = trop_linoz_output(:,:)
    endif

    if(present(hstobie_tropop)) then
      hstobie_tropop(:,:) = trop_trop_output(:,:)
    endif

  end subroutine tropopause_hybridstobie

  ! This routine originates with Stobie at NASA Goddard, but does not have a
  ! known reference. It was supplied by Philip Cameron-Smith of LLNL.
  !
  subroutine tropopause_stobie(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                               tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    integer,            intent(inout)   :: tropLev(:)             ! tropopause level index
    real(r8), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(r8), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(r8), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)

    ! Local Variables
    integer       :: i
    integer       :: k
    integer       :: tLev                     ! tropopause level
    real(r8)      :: tP                       ! tropopause pressure (Pa)
    real(r8)      :: stobie(pver)             ! stobie weighted temperature
    real(r8)      :: sTrop                    ! stobie value at the tropopause

    ! Iterate over all of the columns.
    do i = 1, ncol

      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then

        ! Caclulate a pressure weighted temperature.
        stobie(:) = ALPHA * t(i,:) - log10(pmid(i, :))

        ! Search from the bottom up, looking for the first minimum.
        tLev  = -1

        do k = pver-1, 1, -1

          if (pmid(i, k) <= 4000._r8) then
            exit
          end if

          if (pmid(i, k) >= 55000._r8) then
            cycle
          end if

          if ((tLev == -1) .or. (stobie(k) < sTrop)) then
            tLev  = k
            tP    = pmid(i, k)
            sTrop = stobie(k)
          end if
        end do

        if (tLev /= -1) then
          tropLev(i) = tLev

          ! Return the optional outputs
          if (present(tropP)) tropP(i) = tP

          if (present(tropT)) then
            tropT(i) = tropopause_interpolateT(pver, pmid, t, i, tropLev(i), tP)
          end if

          if (present(tropZ)) then
            tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
          end if
        end if
      end if
    end do

  end subroutine tropopause_stobie


  ! This routine is an implementation of Reichler et al. [2003] done by
  ! Reichler and downloaded from his web site. Minimal modifications were
  ! made to have the routine work within the CAM framework (i.e. using
  ! CAM constants and types).
  !
  ! NOTE: I am not a big fan of the goto's and multiple returns in this
  ! code, but for the moment I have left them to preserve as much of the
  ! original and presumably well tested code as possible.
  ! UPDATE: The most "obvious" substitutions have been made to replace
  ! goto/return statements with cycle/exit. The structure is still
  ! somewhat tangled.
  ! UPDATE 2: "gamma" renamed to "gam" in order to avoid confusion
  ! with the Fortran 2008 intrinsic. "level" argument removed because
  ! a physics column is not contiguous, so using explicit dimensions
  ! will cause the data to be needlessly copied.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! determination of tropopause height from gridded temperature data
  !
  ! Reichler, T., M. Dameris, and R. Sausen (2003),
  ! Determining the tropopause height from gridded data,
  ! Geophys. Res. Lett., 30, 2042, doi:10.1029/2003GL018240, 20.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine twmo(t, p, plimu, pliml, gam, trp)

    real(r8), intent(in), dimension(:)      :: t, p
    real(r8), intent(in)                    :: plimu, pliml, gam
    real(r8), intent(out)                   :: trp

    real(r8), parameter                     :: deltaz = 2000.0_r8

    real(r8)                                :: pmk, pm, a, b, tm, dtdp, dtdz
    real(r8)                                :: ag, bg, ptph
    real(r8)                                :: pm0, pmk0, dtdz0
    real(r8)                                :: p2km, asum, aquer
    real(r8)                                :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
    integer                                 :: level
    integer                                 :: icount, jj
    integer                                 :: j


    trp=-99.0_r8                           ! negative means not valid

    ! initialize start level
    ! dt/dz
    level = size(t)
    pmk= .5_r8 * (p(level-1)**cnst_kap+p(level)**cnst_kap)
    pm = pmk**(1/cnst_kap)
    a = (t(level-1)-t(level))/(p(level-1)**cnst_kap-p(level)**cnst_kap)
    b = t(level)-(a*p(level)**cnst_kap)
    tm = a * pmk + b
    dtdp = a * cnst_kap * (pm**cnst_ka1)
    dtdz = cnst_faktor*dtdp*pm/tm

    main_loop: do j=level-1,2,-1
      pm0 = pm
      pmk0 = pmk
      dtdz0  = dtdz

      ! dt/dz
      pmk= .5_r8 * (p(j-1)**cnst_kap+p(j)**cnst_kap)
      pm = pmk**(1/cnst_kap)
      a = (t(j-1)-t(j))/(p(j-1)**cnst_kap-p(j)**cnst_kap)
      b = t(j)-(a*p(j)**cnst_kap)
      tm = a * pmk + b
      dtdp = a * cnst_kap * (pm**cnst_ka1)
      dtdz = cnst_faktor*dtdp*pm/tm
      ! dt/dz valid?
      if (dtdz.le.gam) cycle main_loop    ! no, dt/dz < -2 K/km
      if (pm.gt.plimu)   cycle main_loop    ! no, too low

      ! dtdz is valid, calculate tropopause pressure
      if (dtdz0.lt.gam) then
        ag = (dtdz-dtdz0) / (pmk-pmk0)
        bg = dtdz0 - (ag * pmk0)
        ptph = exp(log((gam-bg)/ag)/cnst_kap)
      else
        ptph = pm
      endif

      if (ptph.lt.pliml) cycle main_loop
      if (ptph.gt.plimu) cycle main_loop

      ! 2nd test: dtdz above 2 km must not exceed gam
      p2km = ptph + deltaz*(pm/tm)*cnst_faktor     ! p at ptph + 2km
      asum = 0.0_r8                                ! dtdz above
      icount = 0                                   ! number of levels above

      ! test until apm < p2km
      in_loop: do jj=j,2,-1

        pmk2 = .5_r8 * (p(jj-1)**cnst_kap+p(jj)**cnst_kap) ! p mean ^kappa
        pm2 = pmk2**(1/cnst_kap)                           ! p mean
        if(pm2.gt.ptph) cycle in_loop            ! doesn't happen
        if(pm2.lt.p2km) exit in_loop             ! ptropo is valid

        a2 = (t(jj-1)-t(jj))                     ! a
        a2 = a2/(p(jj-1)**cnst_kap-p(jj)**cnst_kap)
        b2 = t(jj)-(a2*p(jj)**cnst_kap)          ! b
        tm2 = a2 * pmk2 + b2                     ! T mean
        dtdp2 = a2 * cnst_kap * (pm2**(cnst_kap-1))  ! dt/dp
        dtdz2 = cnst_faktor*dtdp2*pm2/tm2
        asum = asum+dtdz2
        icount = icount+1
        aquer = asum/float(icount)               ! dt/dz mean

        ! discard ptropo ?
        if (aquer.le.gam) cycle main_loop      ! dt/dz above < gam

      enddo in_loop  ! test next level

      trp = ptph
      exit main_loop
    enddo main_loop

  end subroutine twmo


  ! This routine uses an implementation of Reichler et al. [2003] done by
  ! Reichler and downloaded from his web site. This is similar to the WMO
  ! routines, but is designed for GCMs with a coarse vertical grid.
  subroutine tropopause_twmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                             tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(r8), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(r8), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(r8), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables
    real(r8), parameter     :: gam    = -0.002_r8         ! K/m
    real(r8), parameter     :: plimu    = 45000._r8         ! Pa
    real(r8), parameter     :: pliml    = 7500._r8          ! Pa

    integer                 :: i
    integer                 :: k
    real(r8)                :: tP                       ! tropopause pressure (Pa)

    ! Iterate over all of the columns.
    do i = 1, ncol

      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then

        ! Use the routine from Reichler.
        call twmo(t(i, :), pmid(i, :), plimu, pliml, gam, tP)

        ! if successful, store of the results and find the level and temperature.
        if (tP > 0) then

          ! Find the associated level.
          do k = pver, 2, -1
            if (tP >= pint(i, k)) then
              tropLev(i) = k
              exit
            end if
          end do

          ! Return the optional outputs
          if (present(tropP)) tropP(i) = tP

          if (present(tropT)) then
            tropT(i) = tropopause_interpolateT(pver, pmid, t, i, tropLev(i), tP)
          end if

          if (present(tropZ)) then
            tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
          end if
        end if
      end if
    end do

  end subroutine tropopause_twmo

  ! This routine implements the WMO definition of the tropopause (WMO, 1957; Seidel and Randel, 2006).
  ! Seidel, D. J., and W. J. Randel (2006),
  ! Variability and trends in the global tropopause estimated from radiosonde data,
  ! J. Geophys. Res., 111, D21101, doi:10.1029/2006JD007363.
  !
  ! This requires that the lapse rate be less than 2 K/km for an altitude range
  ! of 2 km. The search starts at the surface and stops the first time this
  ! criteria is met.
  subroutine tropopause_wmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                            tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(r8), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(r8), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(r8), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables
    real(r8), parameter    :: ztrop_low   = 5000._r8        ! lowest tropopause level allowed (m)
    real(r8), parameter    :: ztrop_high  = 20000._r8       ! highest tropopause level allowed (m)
    real(r8), parameter    :: max_dtdz    = 0.002_r8        ! max dt/dz for tropopause level (K/m)
    real(r8), parameter    :: min_trop_dz = 2000._r8        ! min tropopause thickness (m)

    integer                 :: i
    integer                 :: k
    integer                 :: k2
    real(r8)                :: tP                           ! tropopause pressure (Pa)
    real(r8)                :: dt

    ! Iterate over all of the columns.
    do i = 1, ncol

      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then

        kloop: do k = pver-1, 2, -1

          ! Skip levels below the minimum and stop if nothing is found
          ! before the maximum.
          if (zm(i, k) < ztrop_low) then
            cycle kloop
          else if (zm(i, k) > ztrop_high) then
            exit kloop
          end if

          ! Compare the actual lapse rate to the threshold
          dt = t(i, k) - t(i, k-1)

          if (dt <= (max_dtdz * (zm(i, k-1) - zm(i, k)))) then

            ! Make sure that the lapse rate stays below the threshold for the
            ! specified range.
            k2loop: do k2 = k-1, 2, -1
              if ((zm(i, k2) - zm(i, k)) >= min_trop_dz) then
                tP = pmid(i, k)
                tropLev(i) = k
                exit k2loop
              end if

              dt = t(i, k) - t(i, k2)
              if (dt > (max_dtdz * (zm(i, k2) - zm(i, k)))) then
                exit k2loop
              end if
           end do k2loop

           if (tropLev(i) == NOTFOUND) then
              cycle kloop
           else

              ! Return the optional outputs
              if (present(tropP)) tropP(i) = tP

              if (present(tropT)) then
                tropT(i) = tropopause_interpolateT(pver, pmid, t, i, tropLev(i), tP)
              end if

              if (present(tropZ)) then
                tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
              end if

              exit kloop
            end if
          end if
        end do kloop
      end if
    end do

  end subroutine tropopause_wmo


  ! This routine searches for the cold point tropopause, and uses a parabolic
  ! fit of the coldest point and two adjacent points to interpolate the cold point
  ! between model levels.
  subroutine tropopause_cpp(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                            tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(r8), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(r8), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(r8), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables
    real(r8), parameter    :: ztrop_low   = 5000._r8        ! lowest tropopause level allowed (m)
    real(r8), parameter    :: ztrop_high  = 25000._r8       ! highest tropopause level allowed (m)

    integer                 :: i
    integer                 :: k, firstk, lastk
    integer                 :: k2
    real(r8)                :: tZ                           ! tropopause height (m)
    real(r8)                :: tmin
    real(r8)                :: f0, f1, f2
    real(r8)                :: x0, x1, x2
    real(r8)                :: c0, c1, c2
    real(r8)                :: a, b, c

    ! Iterate over all of the columns.
    do i = 1, ncol

      firstk = 0
      lastk  = pver+1

      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then
        tmin = 1e6_r8

        kloop: do k = pver-1, 2, -1

          ! Skip levels below the minimum and stop if nothing is found
          ! before the maximum.
          if (zm(i, k) < ztrop_low) then
            firstk = k
            cycle kloop
          else if (zm(i, k) > ztrop_high) then
            lastk = k
            exit kloop
          end if

          ! Find the coldest point
          if (t(i, k) < tmin) then
            tropLev(i) = k
            tmin = t(i,k)
          end if
        end do kloop


        ! If the minimum is at the edge of the search range, then don't
        ! consider this to be a minima
        if ((tropLev(i) >= (firstk-1)) .or. (tropLev(i) <= (lastk+1))) then
          tropLev(i) = NOTFOUND
        else

          ! If returning P, Z, or T, then do a parabolic fit using the
          ! cold point and it its 2 surrounding points to interpolate
          ! between model levels.
          if (present(tropP) .or. present(tropZ) .or. present(tropT)) then
            f0 = t(i, tropLev(i)-1)
            f1 = t(i, tropLev(i))
            f2 = t(i, tropLev(i)+1)

            x0 = zm(i, tropLev(i)-1)
            x1 = zm(i, tropLev(i))
            x2 = zm(i, tropLev(i)+1)

            c0 = (x0-x1)*(x0-x2)
            c1 = (x1-x0)*(x1-x2)
            c2 = (x2-x0)*(x2-x1)

            ! Determine the quadratic coefficients of:
            !   T = a * z^2 - b*z + c
            a = (f0/c0 + f1/c1 + f2/c2)
            b = (f0/c0*(x1+x2) + f1/c1*(x0+x2) + f2/c2*(x0+x1))
            c = f0/c0*x1*x2 + f1/c1*x0*x2 + f2/c2*x0*x1

            ! Find the altitude of the minimum temperature
            tZ = 0.5_r8 * b / a

            ! The fit should be between the upper and lower points,
            ! so skip the point if the fit fails.
            if ((tZ >= x0) .or. (tZ <= x2)) then
              tropLev(i) = NOTFOUND
            else

              ! Return the optional outputs
              if (present(tropP)) then
                tropP(i) = tropopause_interpolateP(pver, pmid, zm, i, tropLev(i), tZ)
              end if

              if (present(tropT)) then
                tropT(i) = a * tZ*tZ - b*tZ + c
              end if

              if (present(tropZ)) then
                tropZ(i) = tZ
              end if
            end if
          end if
        end if
      end if
    end do

  end subroutine tropopause_cpp

  ! Searches all the columns and attempts to identify the "chemical"
  ! tropopause. This is the lapse rate tropopause, backed up by the climatology
  ! if the lapse rate fails to find the tropopause at pressures higher than a certain
  ! threshold. This pressure threshold depends on latitude. Between 50S and 50N,
  ! the climatology is used if the lapse rate tropopause is not found at P > 75 hPa.
  ! At high latitude (poleward of 50), the threshold is increased to 125 hPa to
  ! eliminate false events that are sometimes detected in the cold polar stratosphere.
  !
  ! NOTE: This routine was adapted from code in chemistry.F90 and mo_gasphase_chemdr.F90.
  ! During the CCPP-ization, findChemTrop is now called from tropopause_find_run using method CHEMTROP
  ! and now also returns the standard tropLev, tropP, tropT, tropZ outputs (optional).
  ! The "backup" option is dropped as it is not used anywhere in current CAM.
  subroutine tropopause_findChemTrop(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                     calday, tropp_p_loc, tropp_days, &
                                     tropLev, tropP, tropT, tropZ, errmsg, errflg)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    real(r8), intent(in)         :: calday        ! Day of year including fraction of day

    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    real(r8), intent(in)         :: tropp_p_loc(:,:)
    real(r8), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12

    integer,                   intent(out)     :: tropLev(:)            ! tropopause level index
    real(r8), optional, intent(inout)   :: tropP(:)              ! tropopause pressure (Pa)
    real(r8), optional, intent(inout)   :: tropT(:)              ! tropopause temperature (K)
    real(r8), optional, intent(inout)   :: tropZ(:)              ! tropopause height (m)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local Variable
    real(r8)     :: dlats(ncol)
    integer             :: i

    errmsg = ' '
    errflg = 0

    ! First use the lapse rate tropopause.
    call tropopause_twmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

    ! Now check high latitudes (poleward of 50) and set the level to the
    ! climatology if the level was not found or is at P <= 125 hPa.
    dlats(:ncol) = lat(:ncol) * cnst_rad2deg ! convert to degrees

    do i = 1, ncol
      if (abs(dlats(i)) > 50._r8) then
        if (tropLev(i) .ne. NOTFOUND) then
          if (pmid(i, tropLev(i)) <= 12500._r8) then
            tropLev(i) = NOTFOUND
          end if
        end if
      end if
    end do

    ! Now use the climatology backup
    if (any(tropLev(:) == NOTFOUND)) then
      call tropopause_climate(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                              calday, tropp_p_loc, tropp_days, &
                              tropLev, tropP, tropT, tropZ)
    end if

  end subroutine tropopause_findChemTrop

  ! This routine interpolates the pressures in the physics state to
  ! find the pressure at the specified tropopause altitude.
  function tropopause_interpolateP(pver, pmid, zm, icol, tropLev, tropZ)

    implicit none

    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    integer, intent(in)                 :: icol               ! column being processed
    integer, intent(in)                 :: tropLev            ! tropopause level index
    real(r8), optional, intent(in)      :: tropZ              ! tropopause pressure (m)
    real(r8)                            :: tropopause_interpolateP

    ! Local Variables
    real(r8)   :: tropP              ! tropopause pressure (Pa)
    real(r8)   :: dlogPdZ            ! dlog(p)/dZ

    ! Interpolate the temperature linearly against log(P)

    ! Is the tropopause at the midpoint?
    if (tropZ == zm(icol, tropLev)) then
      tropP = pmid(icol, tropLev)

    else if (tropZ > zm(icol, tropLev)) then

      ! It is above the midpoint? Make sure we aren't at the top.
      if (tropLev > 1) then
        dlogPdZ = (log(pmid(icol, tropLev)) - log(pmid(icol, tropLev - 1))) / &
          (zm(icol, tropLev) - zm(icol, tropLev - 1))
        tropP = pmid(icol, tropLev) + exp((tropZ - zm(icol, tropLev)) * dlogPdZ)
      end if
    else

      ! It is below the midpoint. Make sure we aren't at the bottom.
      if (tropLev < pver) then
        dlogPdZ =  (log(pmid(icol, tropLev + 1)) - log(pmid(icol, tropLev))) / &
          (zm(icol, tropLev + 1) - zm(icol, tropLev))
        tropP = pmid(icol, tropLev) + exp((tropZ - zm(icol, tropLev)) * dlogPdZ)
      end if
    end if

    tropopause_interpolateP = tropP
  end function tropopause_interpolateP


  ! This routine interpolates the temperatures in the physics state to
  ! find the temperature at the specified tropopause pressure.
  function tropopause_interpolateT(pver, pmid, t, icol, tropLev, tropP)

    implicit none

    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: t(:,:)        ! Temperature (K)
    integer, intent(in)                 :: icol               ! column being processed
    integer, intent(in)                 :: tropLev            ! tropopause level index
    real(r8), optional, intent(in)      :: tropP              ! tropopause pressure (Pa)
    real(r8)                            :: tropopause_interpolateT

    ! Local Variables
    real(r8)   :: tropT              ! tropopause temperature (K)
    real(r8)   :: dTdlogP            ! dT/dlog(P)

    ! Interpolate the temperature linearly against log(P)

    ! Is the tropopause at the midpoint?
    if (tropP == pmid(icol, tropLev)) then
      tropT = t(icol, tropLev)

    else if (tropP < pmid(icol, tropLev)) then

      ! It is above the midpoint? Make sure we aren't at the top.
      if (tropLev > 1) then
        dTdlogP = (t(icol, tropLev) - t(icol, tropLev - 1)) / &
          (log(pmid(icol, tropLev)) - log(pmid(icol, tropLev - 1)))
        tropT = t(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dTdlogP
      end if
    else

      ! It is below the midpoint. Make sure we aren't at the bottom.
      if (tropLev < pver) then
        dTdlogP = (t(icol, tropLev + 1) - t(icol, tropLev)) / &
          (log(pmid(icol, tropLev + 1)) - log(pmid(icol, tropLev)))
        tropT = t(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dTdlogP
      end if
    end if

    tropopause_interpolateT = tropT
  end function tropopause_interpolateT


  ! This routine interpolates the geopotential height in the physics state to
  ! find the geopotential height at the specified tropopause pressure.
  function tropopause_interpolateZ(pint, pmid, zi, zm, phis, icol, tropLev, tropP)

    real(r8), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(r8), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(r8), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(r8), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(r8), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer, intent(in)                 :: icol               ! column being processed
    integer, intent(in)                 :: tropLev            ! tropopause level index
    real(r8), optional, intent(in)      :: tropP              ! tropopause pressure (Pa)
    real(r8)                            :: tropopause_interpolateZ

    ! Local Variables
    real(r8)   :: tropZ              ! tropopause geopotential height (m)
    real(r8)   :: dZdlogP            ! dZ/dlog(P)

    ! Interpolate the geopotential height linearly against log(P)

    ! Is the tropopause at the midpoint?
    if (tropP == pmid(icol, tropLev)) then
      tropZ = zm(icol, tropLev)

    else if (tropP < pmid(icol, tropLev)) then

      ! It is above the midpoint? Make sure we aren't at the top.
      dZdlogP = (zm(icol, tropLev) - zi(icol, tropLev)) / &
        (log(pmid(icol, tropLev)) - log(pint(icol, tropLev)))
      tropZ = zm(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dZdlogP
    else

      ! It is below the midpoint. Make sure we aren't at the bottom.
      dZdlogP = (zm(icol, tropLev) - zi(icol, tropLev+1)) / &
        (log(pmid(icol, tropLev)) - log(pint(icol, tropLev+1)))
      tropZ = zm(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dZdlogP
    end if

    ! PI-atm zi/zm already use the host physics-state reference height.
    tropopause_interpolateZ = tropZ
  end function tropopause_interpolateZ
end module tropopause_find_scheme
