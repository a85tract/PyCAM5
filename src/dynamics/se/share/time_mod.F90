#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module time_mod
  !------------------
  use kinds, only : real_kind
  use iso_c_binding, only : c_int64_t, c_loc, c_ptr
  !------------------
  implicit none
  integer,public                :: nsplit=1
  integer,public                :: nmax          ! Max number of timesteps
  integer,public                :: nEndStep      ! Number of End Step
  integer,public                :: ndays         ! Max number of days
  real (kind=real_kind), public :: tstep         ! Dynamics timestep
  real (kind=real_kind), public :: phys_tscale=0 ! Physics time scale

  real (kind=real_kind), public, parameter :: secphr = 3600.0D0 ! Timestep filter
  real (kind=real_kind), public, parameter :: secpday = 86400.0D0 ! Timestep filter

  ! smooth now in namelist
  real (kind=real_kind), public :: smooth  = 0.05D0    ! Timestep filter
  integer, parameter :: ptimelevels = 3                           ! number of time levels in the dycore

  type, public :: TimeLevel_t
     integer nm1      ! relative time level n-1
     integer n0       ! relative time level n
     integer np1      ! relative time level n+1
     integer nstep    ! time level since simulation start
     integer nstep0   ! timelevel of first complete leapfrog timestep
  end type TimeLevel_t

  ! Methods
  public :: Time_at
  public :: TimeLevel_update
  public :: TimeLevel_init
  public :: TimeLevel_Qdp

  interface TimeLevel_init
     module procedure TimeLevel_init_default
     module procedure TimeLevel_init_specific
     module procedure TimeLevel_init_copy
  end interface

contains

  function Time_at(nstep) result(tat)
    integer, intent(in) :: nstep
    real (kind=real_kind) :: tat
    tat = nstep*tstep
  end function Time_at

  subroutine TimeLevel_init_default(tl)
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use cam_logfile, only : iulog
    type (TimeLevel_t), target, intent(out) :: tl
    interface
       subroutine se_timelevel_init_default_codon(nm1_p, n0_p, np1_p, nstep_p, nstep0_p) &
            bind(c, name='se_timelevel_init_default_codon')
         use iso_c_binding, only : c_ptr
         type(c_ptr), value :: nm1_p, n0_p, np1_p, nstep_p, nstep0_p
       end subroutine se_timelevel_init_default_codon
    end interface

#define SE_MISC_TAG 29
#define SE_MISC_LABEL 'time_mod'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

    call se_timelevel_init_default_codon(c_loc(tl%nm1), c_loc(tl%n0), c_loc(tl%np1), &
         c_loc(tl%nstep), c_loc(tl%nstep0))
    write(iulog,*) 'timelevel_init_default implementation = codon'
  end subroutine TimeLevel_init_default

  subroutine TimeLevel_init_copy(tl, tin)
    type (TimeLevel_t), intent(in) :: tin
    type (TimeLevel_t), intent(out) :: tl
    tl%nm1   = tin%nm1
    tl%n0    = tin%n0
    tl%np1   = tin%np1
    tl%nstep = tin%nstep
    tl%nstep0= tin%nstep0
  end subroutine TimeLevel_init_copy

  subroutine TimeLevel_init_specific(tl,n0,n1,n2,nstep)
    type (TimeLevel_t) :: tl
    integer, intent(in) :: n0,n1,n2,nstep
    tl%nm1= n0
    tl%n0 = n1
    tl%np1= n2
    tl%nstep= nstep
  end subroutine TimeLevel_init_specific


  !this subroutine returns the proper
  !locations for nm1 and n0 for Qdp - because
  !it only has 2 levels for storage
  subroutine TimeLevel_Qdp(tl, qsplit, n0, np1)
    use cam_logfile, only : iulog
    interface
      subroutine timelevel_qdp_codon(nstep_c, qsplit_c, has_np1_c, n0_p, np1_p) &
           bind(c, name='timelevel_qdp_codon')
        import :: c_int64_t, c_ptr
        integer(c_int64_t), value :: nstep_c, qsplit_c, has_np1_c
        type(c_ptr), value :: n0_p, np1_p
      end subroutine timelevel_qdp_codon
    end interface
    type (TimeLevel_t) :: tl
    integer, intent(in) :: qsplit
    integer, intent(inout), target :: n0
    integer, intent(inout), optional, target :: np1
    integer, target :: np1_dummy
    logical, save :: proof_seen = .false.

    if (present(np1)) then
       call timelevel_qdp_codon(int(tl%nstep, c_int64_t), int(qsplit, c_int64_t), &
            1_c_int64_t, c_loc(n0), c_loc(np1))
    else
       np1_dummy = 0
       call timelevel_qdp_codon(int(tl%nstep, c_int64_t), int(qsplit, c_int64_t), &
            0_c_int64_t, c_loc(n0), c_loc(np1_dummy))
    endif
    if (.not. proof_seen) then
       write(iulog,*) 'timelevel_qdp implementation = codon'
       proof_seen = .true.
    endif

    !print * ,'nstep = ', tl%nstep, 'qsplit= ', qsplit, 'i_temp = ', i_temp, 'n0 = ', n0

  end subroutine TimeLevel_Qdp

  subroutine TimeLevel_update(tl,uptype)
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use cam_logfile, only : iulog
    type (TimeLevel_t), target :: tl
    character(len=*)   :: uptype

    ! Local Variable

    integer :: ierr, uptype_code
    interface
       function se_timelevel_update_codon(nm1_p, n0_p, np1_p, nstep_p, uptype_code_c) result(ierr_c) &
            bind(c, name='se_timelevel_update_codon')
         use iso_c_binding, only : c_int64_t, c_ptr
         type(c_ptr), value :: nm1_p, n0_p, np1_p, nstep_p
         integer(c_int64_t), value :: uptype_code_c
         integer(c_int64_t) :: ierr_c
       end function se_timelevel_update_codon
    end interface
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
!$OMP MASTER
#endif
    if (uptype == "leapfrog") then
       uptype_code = 1
    else if (uptype == "forward") then
       uptype_code = 2
    else 
       print *,'WARNING: TimeLevel_update called wint invalid uptype=',uptype
       uptype_code = 0
    end if
    ierr = int(se_timelevel_update_codon(c_loc(tl%nm1), c_loc(tl%n0), c_loc(tl%np1), &
         c_loc(tl%nstep), int(uptype_code, c_int64_t)))
    if (ierr == 0) write(iulog,*) 'timelevel_update implementation = codon'
#if (defined HORIZ_OPENMP)
!$OMP END MASTER
!$OMP BARRIER    
#endif
  end subroutine TimeLevel_update

end module time_mod
