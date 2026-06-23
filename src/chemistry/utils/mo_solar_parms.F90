module mo_solar_parms

  use shr_kind_mod,     only : r8 => shr_kind_r8, shr_kind_cl
  use cam_abortutils,   only : endrun
  use cam_logfile,      only : iulog
  use time_utils,       only : flt_date
  use spmd_utils,       only : masterproc
  use mo_util,          only : chemistry_misc_codon_touch

  implicit none

  private

  public :: solar_parms_readnl
  public :: solar_parms_init
  public :: solar_parms_timestep_init
  public :: solar_parms_get
  public :: solar_parms_on

  save

  integer               :: ntimes
  integer               :: tim_ndx
  integer,  allocatable :: dates(:)
  real(r8)              :: dels
  real(r8), allocatable :: times(:)
  real(r8), allocatable :: f107(:)
  real(r8), allocatable :: f107a(:)
  real(r8), allocatable :: kp(:)
  real(r8), allocatable :: ap(:)

  logical,protected :: solar_parms_on = .false.

  character(len=shr_kind_cl) :: solar_parms_file = ' '     ! solar variability parameters

contains
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  subroutine solar_parms_readnl(nlfile)
    use iso_c_binding, only : c_int64_t

    use namelist_utils, only: find_group_name
    use units,          only: getunit, freeunit
#ifdef SPMD
    use mpishorthand,   only: mpichar, mpicom
#endif

    interface
       function solar_parms_readnl_codon(tag) result(tag_out) bind(c, name='solar_parms_readnl_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function solar_parms_readnl_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.

    ! arguments
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    namelist /solar_parms_nl/ solar_parms_file

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'solar_parms_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, solar_parms_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun('solar_parms_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! broadcast the options to all MPI tasks
    call mpibcast(solar_parms_file, len(solar_parms_file),   mpichar, 0, mpicom)
#endif

    solar_parms_on = len_trim(solar_parms_file)>0
    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('SOLAR_PARMS_READNL_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = solar_parms_readnl_codon(int(139, c_int64_t))
       if (rt_codon_tag_out /= int(139, c_int64_t)) then
          write(iulog,*) 'solar_parms_readnl_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'solar_parms_readnl implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

  end subroutine solar_parms_readnl
  subroutine solar_parms_init ()
    !---------------------------------------------------------------
    !	... initialize solar parmaters
    !---------------------------------------------------------------

    use ioFileMod
    use time_manager,   only: get_curr_date
    use error_messages, only: alloc_err
    use cam_pio_utils,  only: cam_pio_openfile
    use pio,            only: file_desc_t, var_desc_t, pio_get_var, pio_inq_dimid, &
                              pio_inq_varid, pio_closefile, pio_inq_dimlen, pio_nowrite
    use iso_c_binding, only : c_int64_t

    interface
       function solar_parms_init_codon(tag) result(tag_out) bind(c, name='solar_parms_init_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function solar_parms_init_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.


    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    type(file_desc_t)  :: ncid
    integer  :: n
    integer  :: dimid
    type(var_desc_t)  :: varid
    integer  :: astat
    integer  :: wrk_date
    integer  :: yr, mon, day, ncsec
    real(r8) :: wrk_time
    real(r8), allocatable :: bz(:)
    integer  :: ndx(1)
    character(len=256) :: locfn
    integer :: ierr

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('SOLAR_PARMS_INIT_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = solar_parms_init_codon(int(160, c_int64_t))
       if (rt_codon_tag_out /= int(160, c_int64_t)) then
          write(iulog,*) 'solar_parms_init_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'solar_parms_init implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

    if (.not.solar_parms_on) return

    !-----------------------------------------------------------------------
    !	... readin the solar parms dataset
    !-----------------------------------------------------------------------

    if(masterproc) write(iulog,*) 'SOLAR_PARMS: getting file ', trim(solar_parms_file)
    call getfil(solar_parms_file,  locfn, 0)
    if(masterproc) write(iulog,*) 'SOLAR_PARMS: opening file ', trim(locfn)
    call cam_pio_openfile ( ncid, locfn, PIO_NOWRITE)
    ierr = pio_inq_dimid( ncid, 'time', dimid )
    ierr = pio_inq_dimlen( ncid, dimid, ntimes )
    allocate( dates(ntimes), times(ntimes),stat=astat )
    if( astat /= 0 ) then
       call alloc_err( astat, 'solar_parms_init', 'dates,times', ntimes )
    end if
    ierr = pio_inq_varid( ncid, 'date', varid )
    ierr = pio_get_var( ncid, varid, dates )

    do n = 1,ntimes
       times(n) = flt_date( dates(n), 0 )
    end do
    call get_curr_date( yr, mon, day, ncsec )
    wrk_date = 10000*yr + 100*mon + day
    if(masterproc) write(iulog,*) ' '
    if(masterproc) write(iulog,*) '--------------------------------------------------'
    if(masterproc) write(iulog,*) 'solar_parms_init: values for date = ',wrk_date
    wrk_time = flt_date( wrk_date, 0 )
    if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
       write(iulog,*) 'solar_parms_init: initial time is out of range of solar parm times'
       call endrun
    end if
    do n = 2,ntimes
       if( wrk_time <= times(n) ) then
          exit
       end if
    end do
    tim_ndx = n - 1
    dels    = (wrk_time - times(tim_ndx))/(times(tim_ndx+1) - times(tim_ndx))
    if(masterproc) write(iulog,*) 'solar_parms_init: tim_ndx, dels, times(tim_ndx:tim_ndx+1) = ', &
                                                     tim_ndx, dels, dates(tim_ndx:tim_ndx+1)
    if(masterproc) write(iulog,*) '--------------------------------------------------'
    if(masterproc) write(iulog,*) ' '
    !---------------------------------------------------------------
    !	... allocate and read solar parms
    !---------------------------------------------------------------
    allocate( f107(ntimes), f107a(ntimes), &
         kp(ntimes), ap(ntimes), stat=astat )
    if( astat /= 0 ) then
       call alloc_err( astat, 'solar_parms_init', 'f107 ... ap ', ntimes )
    end if
    ierr = pio_inq_varid( ncid, 'f107', varid )
    ierr = pio_get_var( ncid, varid, f107 )
    ierr = pio_inq_varid( ncid, 'f107a', varid )
    ierr = pio_get_var( ncid, varid, f107a )
    ierr = pio_inq_varid( ncid, 'kp', varid )
    ierr = pio_get_var( ncid, varid, kp )
    ierr = pio_inq_varid( ncid, 'ap', varid )
    ierr = pio_get_var( ncid, varid, ap )

    call pio_closefile( ncid )

end subroutine solar_parms_init
subroutine solar_parms_timestep_init
  !---------------------------------------------------------------
  !	... set solar parameters timing
  !---------------------------------------------------------------

 use time_manager,   only : get_curr_date, is_end_curr_day
 use spmd_utils,     only : masterproc
    use iso_c_binding, only : c_int64_t

 implicit none

    interface
       function solar_parms_timestep_init_codon(tag) result(tag_out) bind(c, name='solar_parms_timestep_init_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function solar_parms_timestep_init_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.


 !---------------------------------------------------------------
 !	... local variables
 !---------------------------------------------------------------
 integer  :: n
 integer  :: wrk_date
 integer  :: yr, mon, day, ncsec
 real(r8) :: wrk_time

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('SOLAR_PARMS_TIMESTEP_INIT_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = solar_parms_timestep_init_codon(int(161, c_int64_t))
       if (rt_codon_tag_out /= int(161, c_int64_t)) then
          write(iulog,*) 'solar_parms_timestep_init_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'solar_parms_timestep_init implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

 if (.not.solar_parms_on) return

 if( is_end_curr_day() ) then
    call get_curr_date( yr, mon, day, ncsec )
    wrk_date = 10000*yr + 100*mon + day
    if (masterproc) &
         write(iulog,*) 'solar_parms_timestep_init: values for date = ',wrk_date
    wrk_time = flt_date( wrk_date, 0 )
    if( wrk_time < times(1) .or. wrk_time > times(ntimes) ) then
       write(iulog,*) 'solar_parms_timestep_init: time is out of range of solar parm times'
       call endrun('solar_parms_timestep_init: time is out of range of solar parm times')
    end if
    do n = 2,ntimes
       if( wrk_time <= times(n) ) then
          exit
       end if
    end do
    tim_ndx = n - 1
    dels    = (wrk_time - times(tim_ndx))/(times(tim_ndx+1) - times(tim_ndx))
 end if


end subroutine solar_parms_timestep_init
subroutine solar_parms_get( f107_s, f107a_s, ap_s, kp_s, hp_s )
    use iso_c_binding, only : c_int64_t
  !---------------------------------------------------------------
  !	... set,retrieve solar parmaters
  !---------------------------------------------------------------

 implicit none

    interface
       function solar_parms_get_codon(tag) result(tag_out) bind(c, name='solar_parms_get_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function solar_parms_get_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.


 !---------------------------------------------------------------
 !	... dummy arguments
 !---------------------------------------------------------------
 real(r8), optional, intent(out) :: f107_s                   ! solar euv factor
 real(r8), optional, intent(out) :: f107a_s                  ! averaged solar euv factor
 real(r8), optional, intent(out) :: ap_s                     ! solar mag factor
 real(r8), optional, intent(out) :: kp_s                     ! solar mag factor
 real(r8), optional, intent(out) :: hp_s                     ! hemispheric power

 !---------------------------------------------------------------
 !	... local variables
 !---------------------------------------------------------------
 integer  :: tnp
 real(r8) :: wkp                                             ! wrk solar mag factor

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('SOLAR_PARMS_GET_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = solar_parms_get_codon(int(162, c_int64_t))
       if (rt_codon_tag_out /= int(162, c_int64_t)) then
          write(iulog,*) 'solar_parms_get_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'solar_parms_get implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

 if (.not.solar_parms_on) return

 tnp = tim_ndx + 1
 if( present( f107_s ) ) then
    f107_s  =  f107(tim_ndx) + dels*(f107(tnp) - f107(tim_ndx))
 end if
 if( present( f107a_s ) ) then
    f107a_s  =  f107a(tim_ndx) + dels*(f107a(tnp) - f107a(tim_ndx))
 end if
 if( present( kp_s ) ) then
    kp_s  =  kp(tim_ndx) + dels*(kp(tnp) - kp(tim_ndx))
 end if
 if( present( ap_s ) ) then
    ap_s  =  ap(tim_ndx) + dels*(ap(tnp) - ap(tim_ndx))
 end if
 if( present( hp_s ) ) then
    wkp  =  kp(tim_ndx) + dels*(kp(tnp) - kp(tim_ndx))
    hp_s = max( 3._r8,-2.78_r8 + 9.39_r8*wkp )
 end if

end subroutine solar_parms_get

end module mo_solar_parms
