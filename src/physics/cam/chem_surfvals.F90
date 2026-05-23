
module chem_surfvals

!-----------------------------------------------------------------------------------
! Purpose: Provides greenhouse gas (ghg) values at the Earth's surface.
!          These values may be time dependent.
!
! Author: Brian Eaton (assembled module from existing scattered code pieces)
!-----------------------------------------------------------------------------------

   use shr_kind_mod,   only: r8=>shr_kind_r8
   use spmd_utils,     only: masterproc
   use time_manager,   only: get_curr_date, get_start_date, is_end_curr_day, &
                             timemgr_datediff, get_curr_calday
   use cam_abortutils, only: endrun
   use netcdf
   use error_messages, only: handle_ncerr  
   use cam_logfile,    only: iulog
   use m_types,        only: time_ramp
   use constituents,   only: pcnst
   use iso_c_binding,  only: c_double, c_int64_t

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make default access private
   save

! Public methods
   public ::&
      chem_surfvals_readnl,  &! read namelist input
      chem_surfvals_init,    &! initialize options that depend on namelist input
      chem_surfvals_set,     &! set ghg surface values when scenario_ghg is 'RAMPED' or 'CHEM_LBC_FILE'
      chem_surfvals_get,     &! return surface values for: CO2VMR, CO2MMR, CH4VMR
                              ! N2OVMR, F11VMR, and F12VMR
      chem_surfvals_co2_rad   ! return co2 for radiation

   public :: flbc_list

! Private module data

   ! Default values for namelist variables -- now set by build-namelist
   real(r8) :: o2mmr = .23143_r8               ! o2 mass mixing ratio
   real(r8) :: co2vmr_rad = -1.0_r8            ! co2 vmr override for radiation
   real(r8) :: co2vmr = -1.0_r8                ! co2   volume mixing ratio 
   real(r8) :: n2ovmr = -1.0_r8                ! n2o   volume mixing ratio 
   real(r8) :: ch4vmr = -1.0_r8                ! ch4   volume mixing ratio 
   real(r8) :: f11vmr = -1.0_r8                ! cfc11 volume mixing ratio 
   real(r8) :: f12vmr = -1.0_r8                ! cfc12 volume mixing ratio 
   character(len=16) :: scenario_ghg = 'FIXED' ! 'FIXED','RAMPED' or 'RAMP_CO2_ONLY'
   integer  :: rampYear_ghg = 0                ! ramped gases fixed at this year (if > 0)
   character(len=256) :: bndtvghg = ' '        ! filename for ramped data
   integer  :: ramp_co2_start_ymd = 0          ! start date for co2 ramping (yyyymmdd)
   real(r8) :: ramp_co2_annual_rate = 1.0_r8      ! % amount of co2 ramping per yr; default is 1% 
   real(r8) :: ramp_co2_cap = -9999.0_r8          ! co2 ramp cap if rate>0, floor otherwise 
                                               ! as multiple or fraction of inital value
                                               ! ex. 4.0 => cap at 4x initial co2 setting 
   integer  :: ghg_yearStart_model = 0         ! model start year
   integer  :: ghg_yearStart_data  = 0         ! data  start year   

   logical  :: ghg_use_calendar                ! true => data year = model year
   logical :: doRamp_ghg    ! true => turn on ramping for ghg
   logical :: ramp_just_co2 ! true => ramping to be done just for co2 and not other ghg's
   integer :: fixYear_ghg   ! year at which Ramped gases are fixed
   integer :: co2_start     ! date at which co2 begins ramping
   real(r8) :: co2_daily_factor    ! daily multiplier to achieve annual rate of co2 ramp
   real(r8) :: co2_limit    ! value of co2vmr where ramping ends
   real(r8) :: co2_base     ! initial co2 volume mixing ratio, before any ramping
   integer :: ntim = -1               ! number of yearly data values
   integer,  allocatable, target :: yrdata(:) ! yearly data values
   real(r8), allocatable, target :: co2(:)    ! co2 mixing ratios in ppmv 
   real(r8), allocatable, target :: ch4(:)    ! ppbv
   real(r8), allocatable, target :: n2o(:)    ! ppbv
   real(r8), allocatable, target :: f11(:)    ! pptv
   real(r8), allocatable, target :: f12(:)    ! pptv
   real(r8), allocatable, target :: adj(:)    ! unitless adjustment factor for f11 & f12
   
   ! fixed lower boundary 
   
   character(len=256) :: flbc_file = ' '
   character(len=16)  :: flbc_list(pcnst) = ''
   type(time_ramp)    :: flbc_timing     != time_ramp( "CYCLICAL",  19970101, 0 )
   logical :: use_native_impl = .false.
   logical :: impl_selected = .false.

   interface
      function chem_surfvals_readnl_codon(flag_c) result(out_c) bind(c, name="chem_surfvals_readnl_codon")
         use iso_c_binding, only: c_int64_t
         integer(c_int64_t), value :: flag_c
         integer(c_int64_t) :: out_c
      end function chem_surfvals_readnl_codon
      function chem_surfvals_init_codon(flag_c) result(out_c) bind(c, name="chem_surfvals_init_codon")
         use iso_c_binding, only: c_int64_t
         integer(c_int64_t), value :: flag_c
         integer(c_int64_t) :: out_c
      end function chem_surfvals_init_codon
      function chem_surfvals_get_codon(value_c) result(out_c) bind(c, name="chem_surfvals_get_codon")
         use iso_c_binding, only: c_double
         real(c_double), value :: value_c
         real(c_double) :: out_c
      end function chem_surfvals_get_codon
      function chem_surfvals_co2_rad_codon(value_c) result(out_c) bind(c, name="chem_surfvals_co2_rad_codon")
         use iso_c_binding, only: c_double
         real(c_double), value :: value_c
         real(c_double) :: out_c
      end function chem_surfvals_co2_rad_codon
   end interface

!=========================================================================================
contains
!=========================================================================================

subroutine chem_surfvals_readnl(nlfile)

   ! Read chem_surfvals_nl namelist group.

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr, i
   character(len=*), parameter :: subname = 'chem_surfvals_readnl'
   
   character(len=8)   :: flbc_type = 'CYCLICAL'     ! 'CYCLICAL' | 'SERIAL' | 'FIXED'
   integer            :: flbc_cycle_yr = 0
   integer            :: flbc_fixed_ymd = 0
   integer            :: flbc_fixed_tod = 0
   integer(c_int64_t) :: active_c

   namelist /chem_surfvals_nl/ co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, &
                               co2vmr_rad, scenario_ghg, rampyear_ghg, bndtvghg, &
                               ramp_co2_start_ymd, ramp_co2_annual_rate, ramp_co2_cap, &
                               ghg_yearStart_model, ghg_yearStart_data
   ! waccm/cam-chem naemlist
   namelist /chem_surfvals_nl/ flbc_type, flbc_cycle_yr, flbc_fixed_ymd, flbc_fixed_tod, flbc_list, flbc_file

   !-----------------------------------------------------------------------------
   active_c = chem_surfvals_readnl_codon(1_c_int64_t)
   if (active_c == 0_c_int64_t) return

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'chem_surfvals_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, chem_surfvals_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast (co2vmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (n2ovmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (ch4vmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (f11vmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (f12vmr,                          1,   mpir8, 0, mpicom)
   call mpibcast (co2vmr_rad,                      1,   mpir8, 0, mpicom)
   call mpibcast (scenario_ghg,    len(scenario_ghg), mpichar, 0, mpicom)
   call mpibcast (rampyear_ghg,                    1,  mpiint, 0, mpicom)
   call mpibcast (bndtvghg,            len(bndtvghg), mpichar, 0, mpicom)
   call mpibcast (ramp_co2_start_ymd,              1,  mpiint, 0, mpicom)
   call mpibcast (ramp_co2_annual_rate,            1,   mpir8, 0, mpicom)
   call mpibcast (ramp_co2_cap,                    1,   mpir8, 0, mpicom)
   call mpibcast (ghg_yearstart_model,             1,  mpiint, 0, mpicom)
   call mpibcast (ghg_yearstart_data,              1,  mpiint, 0, mpicom)
   
   ! waccm/cam-chem fixed lower boundary 
   
   call mpibcast (flbc_type,         len(flbc_type),                  mpichar, 0, mpicom)
   call mpibcast (flbc_cycle_yr,     1,                               mpiint,  0, mpicom)
   call mpibcast (flbc_fixed_ymd,    1,                               mpiint,  0, mpicom)
   call mpibcast (flbc_fixed_tod,    1,                               mpiint,  0, mpicom)
   call mpibcast (flbc_list,         len(flbc_list(1))*pcnst,         mpichar, 0, mpicom)
   call mpibcast (flbc_file,         len(flbc_file),                  mpichar, 0, mpicom)

#endif

   flbc_timing%type      = flbc_type
   flbc_timing%cycle_yr  = flbc_cycle_yr
   flbc_timing%fixed_ymd = flbc_fixed_ymd
   flbc_timing%fixed_tod = flbc_fixed_tod

   if ( len_trim(bndtvghg) > 0 .and. len_trim(flbc_file) > 0 ) then
      call endrun('chem_surfvals_readnl: Cannot specify both bndtvghg and flbc_file ')
   endif

   if (co2vmr_rad > 0._r8) then
      if (masterproc) &
         write(iulog,*) trim(subname)//': co2vmr_rad override is set to ', co2vmr_rad
   end if

end subroutine chem_surfvals_readnl

!================================================================================================

subroutine chem_surfvals_init()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize the ramp options that are controlled by namelist input.
! Set surface values at initial time.
! N.B. This routine must be called after the time manager has been initialized
!      since chem_surfvals_set calls time manager methods.
! 
! Author: B. Eaton - merged code from parse_namelist and rampnl_ghg.
! 
!-----------------------------------------------------------------------

   use infnan,  only : posinf, assignment(=)
   use mo_flbc, only : flbc_inti

   !---------------------------Local variables-----------------------------
   integer :: yr, mon, day, ncsec
   integer(c_int64_t) :: active_c
   !-----------------------------------------------------------------------
   active_c = chem_surfvals_init_codon(1_c_int64_t)
   if (active_c == 0_c_int64_t) return

   if (scenario_ghg == 'FIXED') then
      doRamp_ghg = .false.
      ramp_just_co2 = .false.
      if (masterproc) &
         write(iulog,*)'chem_surfvals_init: ghg surface values are fixed as follows'

   else if (scenario_ghg == 'RAMPED') then
      doRamp_ghg = .true.
      ramp_just_co2 = .false.
      call ghg_ramp_read

      fixYear_ghg = rampYear_ghg     ! set private member to namelist var
      if (masterproc) then
         if ( fixYear_ghg > 0 ) then
            write(iulog,*) '  FIXED values from year ',fixYear_ghg
         else
            write(iulog,*) '  RAMPED values initialized to'
         end if
      end if
      call chem_surfvals_set()

   else if (scenario_ghg == 'RAMP_CO2_ONLY') then
      if(ramp_co2_start_ymd == 0) then
         ! by default start the ramp at the initial run time
         call get_start_date(yr, mon, day, ncsec)
         ramp_co2_start_ymd = yr*10000 + mon*100 + day
      end if
      co2_start = ramp_co2_start_ymd

      if(ramp_co2_annual_rate <= -100.0_r8) then
         write(iulog,*) 'RAMP_CO2:  invalid ramp_co2_annual_rate= ',ramp_co2_annual_rate
         call endrun ('chem_surfvals_init: RAMP_CO2_ANNUAL_RATE must be greater than -100.0')
      end if

      doRamp_ghg = .true.
      ramp_just_co2 = .true.
      co2_base = co2vmr        ! save initial setting 
      if (masterproc) &
           write(iulog,*) '  RAMPED values initialized to'

      co2_daily_factor = (ramp_co2_annual_rate*0.01_r8+1.0_r8)**(1.0_r8/365.0_r8)

      if(ramp_co2_cap > 0.0_r8) then  
         co2_limit = ramp_co2_cap * co2_base
      else                                  ! if no cap/floor specified, provide default
         if(ramp_co2_annual_rate < 0.0_r8) then
            co2_limit = 0.0_r8
         else
            co2_limit = posinf
         end if
      end if
      if((ramp_co2_annual_rate<0.0_r8 .and. co2_limit>co2_base) .or. &
         (ramp_co2_annual_rate>0.0_r8 .and. co2_limit<co2_base)) then
         write(iulog,*) 'RAMP_CO2: ramp_co2_cap is unreachable'
         write(iulog,*) 'RAMP_CO2: ramp_co2_annual_rate= ',ramp_co2_annual_rate,' ramp_co2_cap= ',ramp_co2_cap
         call endrun('chem_surfvals_init:  ramp_co2_annual_rate and ramp_co2_cap incompatible')
      end if

      call chem_surfvals_set()
   else if (scenario_ghg == 'CHEM_LBC_FILE') then
      ! set by lower boundary conditions file
      call flbc_inti( flbc_file, flbc_list, flbc_timing, co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr )
      call chem_surfvals_set()
   else
      call endrun ('chem_surfvals_init: input namelist SCENARIO_GHG must be set to either FIXED, RAMPED, RAMP_CO2_ONLY, &
                   & or CHEM_LBC_FILE')
   endif

   if (masterproc) then
      write(iulog,*) '  co2 volume mixing ratio = ',co2vmr
      write(iulog,*) '  ch4 volume mixing ratio = ',ch4vmr
      write(iulog,*) '  n2o volume mixing ratio = ',n2ovmr
      write(iulog,*) '  f11 volume mixing ratio = ',f11vmr
      write(iulog,*) '  f12 volume mixing ratio = ',f12vmr
   end if

end subroutine chem_surfvals_init

!=========================================================================================

subroutine ghg_ramp_read()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read ramped greenhouse gas surface data.  
! 
! Author: T. Henderson
! 
!-----------------------------------------------------------------------

   use ioFileMod, only: getfil
#if ( defined SPMD )
   use mpishorthand, only: mpicom, mpiint, mpir8
#endif
   character(len=*), parameter :: subname = 'ghg_ramp_read'

!---------------------------Local variables-----------------------------
   integer :: ncid
   integer :: co2_id
   integer :: ch4_id
   integer :: n2o_id
   integer :: f11_id
   integer :: f12_id
   integer :: adj_id
   integer :: date_id
   integer :: time_id
   integer :: ierror
   character(len=256) :: locfn          ! netcdf local filename to open

   if (masterproc) then
     call getfil (bndtvghg, locfn, 0)
     call handle_ncerr( nf90_open (trim(locfn), NF90_NOWRITE, ncid),subname,__LINE__)

     write(iulog,*)'GHG_RAMP_READ:  reading ramped greenhouse gas surface data from file ',trim(locfn)

     call handle_ncerr( nf90_inq_varid( ncid, 'date', date_id ),subname,__LINE__)
     call handle_ncerr( nf90_inq_varid( ncid, 'CO2', co2_id ),subname,__LINE__)
     call handle_ncerr( nf90_inq_varid( ncid, 'CH4', ch4_id ),subname,__LINE__)
     call handle_ncerr( nf90_inq_varid( ncid, 'N2O', n2o_id ),subname,__LINE__)
     call handle_ncerr( nf90_inq_varid( ncid, 'f11', f11_id ),subname,__LINE__)
     call handle_ncerr( nf90_inq_varid( ncid, 'f12', f12_id ),subname,__LINE__)
     call handle_ncerr( nf90_inq_varid( ncid, 'adj', adj_id ),subname,__LINE__)
     call handle_ncerr( nf90_inq_dimid( ncid, 'time', time_id ),subname,__LINE__)
     call handle_ncerr( nf90_inquire_dimension( ncid, time_id, len=ntim ),subname,__LINE__)

   endif
#if (defined SPMD )
   call mpibcast (ntim, 1, mpiint, 0, mpicom)
#endif
   ! these arrays are never deallocated
   allocate ( yrdata(ntim), co2(ntim), ch4(ntim), n2o(ntim),    &
                 f11(ntim), f12(ntim), adj(ntim), stat=ierror )
   if (ierror /= 0) then
     write(iulog,*)'GHG_RAMP_READ:  ERROR, allocate() failed!'
     call endrun
   endif
   if (masterproc) then
     call handle_ncerr( nf90_get_var (ncid, date_id, yrdata ),subname,__LINE__)
     yrdata = yrdata / 10000
     call handle_ncerr( nf90_get_var (ncid, co2_id, co2 ),subname,__LINE__)
     call handle_ncerr( nf90_get_var (ncid, ch4_id, ch4 ),subname,__LINE__)
     call handle_ncerr( nf90_get_var (ncid, n2o_id, n2o ),subname,__LINE__)
     call handle_ncerr( nf90_get_var (ncid, f11_id, f11 ),subname,__LINE__)
     call handle_ncerr( nf90_get_var (ncid, f12_id, f12 ),subname,__LINE__)
     call handle_ncerr( nf90_get_var (ncid, adj_id, adj ),subname,__LINE__)
     call handle_ncerr( nf90_close (ncid),subname,__LINE__)
     write(iulog,*)'GHG_RAMP_READ:  successfully read ramped greenhouse gas surface data from years ',&
	yrdata(1),' through ',yrdata(ntim)
   endif
#if (defined SPMD )
   call mpibcast (co2, ntim, mpir8, 0, mpicom)
   call mpibcast (ch4, ntim, mpir8, 0, mpicom)
   call mpibcast (n2o, ntim, mpir8, 0, mpicom)
   call mpibcast (f11, ntim, mpir8, 0, mpicom)
   call mpibcast (f12, ntim, mpir8, 0, mpicom)
   call mpibcast (adj, ntim, mpir8, 0, mpicom)
   call mpibcast (yrdata, ntim, mpiint, 0, mpicom)
#endif

   return

end subroutine ghg_ramp_read

!=========================================================================================

function chem_surfvals_get(name)
  use physconst,    only: mwdry, mwco2

  character(len=*), intent(in) :: name

  real(r8) :: rmwco2 
  real(r8) :: chem_surfvals_get

  rmwco2 = mwco2/mwdry    ! ratio of molecular weights of co2 to dry air
  select case (name)
  case ('CO2VMR')
     chem_surfvals_get = co2vmr
  case ('CO2MMR')
     chem_surfvals_get = rmwco2 * co2vmr
  case ('N2OVMR')
     chem_surfvals_get = n2ovmr
  case ('CH4VMR')
     chem_surfvals_get = ch4vmr
  case ('F11VMR')
     chem_surfvals_get = f11vmr
  case ('F12VMR')
     chem_surfvals_get = f12vmr
  case ('O2MMR')
     chem_surfvals_get = o2mmr
  case default
     call endrun('chem_surfvals_get does not know name')
  end select
  chem_surfvals_get = chem_surfvals_get_codon(real(chem_surfvals_get, c_double))

end function chem_surfvals_get


!=========================================================================================

function chem_surfvals_co2_rad(vmr_in)
 
   ! Return the value of CO2 (as mmr) that is radiatively active.

   ! This method is used by ghg_data to set the prescribed value of CO2 in
   ! the physics buffer.  If the user has set the co2vmr_rad namelist
   ! variable then that value will override either the value set by the
   ! co2vmr namelist variable, or the values time interpolated from a
   ! dataset.
   
   ! This method is also used by cam_history to write the radiatively active
   ! CO2 to the history file.  The optional argument allows returning the
   ! value as vmr.

   use physconst,    only: mwdry, mwco2

   ! Arguments
   logical, intent(in), optional :: vmr_in  ! return CO2 as vmr

   ! Return value
   real(r8) :: chem_surfvals_co2_rad

   ! Local variables
   real(r8) :: convert_vmr      ! convert vmr to desired output
   !-----------------------------------------------------------------------

   ! by default convert vmr to mmr
   convert_vmr = mwco2/mwdry    ! ratio of molecular weights of co2 to dry air
   if (present(vmr_in)) then
      ! if request return vmr
      if (vmr_in) convert_vmr = 1.0_r8
   end if

   if (co2vmr_rad > 0._r8) then
      chem_surfvals_co2_rad = convert_vmr * co2vmr_rad
   else                           
      chem_surfvals_co2_rad = convert_vmr * co2vmr     
   end if
   chem_surfvals_co2_rad = chem_surfvals_co2_rad_codon(real(chem_surfvals_co2_rad, c_double))

end function chem_surfvals_co2_rad

!=========================================================================================

subroutine chem_surfvals_set()

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   use ppgrid,         only: begchunk, endchunk
   use mo_flbc,        only: flbc_gmean_vmr, flbc_chk

!---------------------------Local variables-----------------------------

   integer  :: yr, mon, day, ncsec ! components of a date
   integer  :: ncdate              ! current date in integer format [yyyymmdd]
   real(r8) :: calday
   real(r8) :: daydiff
   integer(c_int64_t), target :: status_code
   integer(c_int64_t), allocatable, target :: yrdata64(:)
   real(c_double), target :: co2vmr_work
   real(c_double), target :: ch4vmr_work
   real(c_double), target :: n2ovmr_work
   real(c_double), target :: f11vmr_work
   real(c_double), target :: f12vmr_work

   interface
      subroutine chem_surfvals_set_all_codon(fixYear_ghg_c, ghg_yearStart_model_c, ghg_yearStart_data_c, &
           yr_c, calday_c, ntim_c, yrdata_p, co2_p, ch4_p, n2o_p, f11_p, f12_p, adj_p, &
           co2vmr_p, ch4vmr_p, n2ovmr_p, f11vmr_p, f12vmr_p, status_p) bind(c, name="chem_surfvals_set_all_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: fixYear_ghg_c, ghg_yearStart_model_c, ghg_yearStart_data_c
         integer(c_int64_t), value :: yr_c, ntim_c
         real(c_double), value :: calday_c
         type(c_ptr), value :: yrdata_p, co2_p, ch4_p, n2o_p, f11_p, f12_p, adj_p
         type(c_ptr), value :: co2vmr_p, ch4vmr_p, n2ovmr_p, f11vmr_p, f12vmr_p, status_p
      end subroutine chem_surfvals_set_all_codon
      subroutine chem_surfvals_set_co2_codon(daydiff_c, co2_base_c, co2_daily_factor_c, &
           co2_limit_c, co2vmr_p) bind(c, name="chem_surfvals_set_co2_codon")
         use iso_c_binding, only: c_double, c_ptr
         real(c_double), value :: daydiff_c, co2_base_c, co2_daily_factor_c, co2_limit_c
         type(c_ptr), value :: co2vmr_p
      end subroutine chem_surfvals_set_co2_codon
   end interface
   
   call chem_surfvals_select_impl()

   if (use_native_impl) then
      call chem_surfvals_set_native()
      return
   end if

   if ( doRamp_ghg ) then
      if(ramp_just_co2) then
         call get_curr_date(yr, mon, day, ncsec)
         ncdate = yr*10000 + mon*100 + day
         call timemgr_datediff(co2_start, 0, ncdate, ncsec, daydiff)
         co2vmr_work = real(co2vmr, c_double)
         call chem_surfvals_set_co2_codon(real(daydiff, c_double), real(co2_base, c_double), &
              real(co2_daily_factor, c_double), real(co2_limit, c_double), c_loc(co2vmr_work))
         co2vmr = real(co2vmr_work, r8)
      else
         if (ntim <= 0 .or. .not. allocated(yrdata) .or. .not. allocated(co2) .or. .not. allocated(ch4) .or. &
             .not. allocated(n2o) .or. .not. allocated(f11) .or. .not. allocated(f12) .or. .not. allocated(adj)) then
            call chem_surfvals_set_native()
            return
         end if

         call get_curr_date(yr, mon, day, ncsec)
         calday = get_curr_calday()

         allocate(yrdata64(ntim))
         yrdata64(:) = int(yrdata(:), c_int64_t)

         co2vmr_work = real(co2vmr, c_double)
         ch4vmr_work = real(ch4vmr, c_double)
         n2ovmr_work = real(n2ovmr, c_double)
         f11vmr_work = real(f11vmr, c_double)
         f12vmr_work = real(f12vmr, c_double)
         status_code = 0_c_int64_t

         call chem_surfvals_set_all_codon( &
              int(fixYear_ghg, c_int64_t), int(ghg_yearStart_model, c_int64_t), int(ghg_yearStart_data, c_int64_t), &
              int(yr, c_int64_t), real(calday, c_double), int(ntim, c_int64_t), &
              c_loc(yrdata64), c_loc(co2), c_loc(ch4), c_loc(n2o), c_loc(f11), c_loc(f12), c_loc(adj), &
              c_loc(co2vmr_work), c_loc(ch4vmr_work), c_loc(n2ovmr_work), c_loc(f11vmr_work), c_loc(f12vmr_work), &
              c_loc(status_code) &
         )

         if (status_code /= 0_c_int64_t) then
            deallocate(yrdata64)
            call chem_surfvals_set_native()
            return
         end if

         co2vmr = real(co2vmr_work, r8)
         ch4vmr = real(ch4vmr_work, r8)
         n2ovmr = real(n2ovmr_work, r8)
         f11vmr = real(f11vmr_work, r8)
         f12vmr = real(f12vmr_work, r8)
         deallocate(yrdata64)
      end if
   elseif (scenario_ghg == 'CHEM_LBC_FILE') then
      call chem_surfvals_set_native()
      return
   endif

   if (masterproc .and. is_end_curr_day()) then
      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day
      write(iulog,*) 'chem_surfvals_set: ncdate= ',ncdate,' co2vmr=',co2vmr

      if (.not. ramp_just_co2 .and. mon==1 .and. day==1) then
         write(iulog,*) 'chem_surfvals_set: ch4vmr=', ch4vmr, ' n2ovmr=', n2ovmr, &
                        ' f11vmr=', f11vmr, ' f12vmr=', f12vmr
      end if

   end if

   return
end subroutine chem_surfvals_set

!=========================================================================================

subroutine chem_surfvals_set_native()

   use ppgrid,         only: begchunk, endchunk
   use mo_flbc,        only: flbc_gmean_vmr, flbc_chk

!---------------------------Local variables-----------------------------

   integer  :: yr, mon, day, ncsec ! components of a date
   integer  :: ncdate              ! current date in integer format [yyyymmdd]
   
   if ( doRamp_ghg ) then
      if(ramp_just_co2) then
         call chem_surfvals_set_co2()
      else
         call chem_surfvals_set_all()
      end if
   elseif (scenario_ghg == 'CHEM_LBC_FILE') then
      ! set mixing ratios from cam-chem/waccm lbc file 
      call flbc_chk()
      call flbc_gmean_vmr(co2vmr,ch4vmr,n2ovmr,f11vmr,f12vmr)
   endif

   if (masterproc .and. is_end_curr_day()) then
      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day
      write(iulog,*) 'chem_surfvals_set: ncdate= ',ncdate,' co2vmr=',co2vmr

      if (.not. ramp_just_co2 .and. mon==1 .and. day==1) then
         write(iulog,*) 'chem_surfvals_set: ch4vmr=', ch4vmr, ' n2ovmr=', n2ovmr, &
                        ' f11vmr=', f11vmr, ' f12vmr=', f12vmr
      end if

   end if

   return
end subroutine chem_surfvals_set_native

!=========================================================================================

subroutine chem_surfvals_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CHEM_SURFVALS_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_impl = .false.
   end if

   impl_selected = .true.

   if (masterproc) then
      if (use_native_impl) then
         write(iulog,*) 'chem_surfvals implementation = native'
      else
         write(iulog,*) 'chem_surfvals implementation = codon'
      end if
   end if

end subroutine chem_surfvals_select_impl

!=========================================================================================

subroutine chem_surfvals_set_all()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes greenhouse gas volume mixing ratios via interpolation of
! yearly input data.
! 
! Author: B. Eaton - updated ramp_ghg for use in chem_surfvals module
! 
!-----------------------------------------------------------------------
   use interpolate_data, only: get_timeinterp_factors

!---------------------------Local variables-----------------------------

   integer yrmodel           ! model year
   integer nyrm              ! year index
   integer nyrp              ! year index
   integer :: yr, mon, day   ! components of a date
   integer :: ncdate         ! current date in integer format [yyyymmdd]
   integer :: ncsec          ! current time of day [seconds]

   real(r8) :: calday            ! current calendar day
   real(r8) doymodel             ! model day of year
   real(r8) doydatam             ! day of year for input data yrdata(nyrm)
   real(r8) doydatap             ! day or year for input data yrdata(nyrp)
   real(r8) deltat               ! delta time
   real(r8) fact1, fact2         ! time interpolation factors
   real(r8) cfcscl               ! cfc scale factor for f11

   integer yearRan_model         ! model ran year
!
! ---------------------------------------------------------------------
!
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day
!
! determine ghg_use_calendar      
!
   if ( ghg_yearStart_model > 0 .and. ghg_yearStart_data > 0 ) then
      ghg_use_calendar = .false.
   else
      ghg_use_calendar = .true.
   end if
!
! determine index into input data
!
   if ( fixYear_ghg > 0) then
      yrmodel  = fixYear_ghg
      nyrm = fixYear_ghg - yrdata(1) + 1
   else
      if ( ghg_use_calendar) then
         yrmodel  = yr          
         nyrm = yr - yrdata(1) + 1
      else 
         yearRan_model = yr - ghg_yearStart_model
         if ( yearRan_model < 0 ) then
            call endrun('chem_surfvals_set_all: incorrect ghg_yearStart_model')
         endif
         yrmodel  = yearRan_model + ghg_yearStart_data
 
         nyrm = ghg_yearStart_data + yearRan_model - yrdata(1) + 1
      end if
   end if

   nyrp       = nyrm + 1
!
! if current date is before yrdata(1), quit
!
   if (nyrm < 1) then
      write(iulog,*)'chem_surfvals_set_all: data time index is out of bounds'
      write(iulog,*)'nyrm = ',nyrm,' nyrp= ',nyrp, ' ncdate= ', ncdate
      call endrun
   endif
!
! if current date later than yrdata(ntim), call endrun.
! if want to use ntim values - uncomment the following lines
! below and comment the call to endrun and previous write
!
   if (nyrp > ntim) then
      call endrun ('chem_surfvals_set_all: error - current date is past the end of valid data')
!         write(iulog,*)'chem_surfvals_set_all: using ghg data for ',yrdata(ntim)
!         co2vmr = co2(ntim)*1.e-06
!         ch4vmr = ch4(ntim)*1.e-09
!         n2ovmr = n2o(ntim)*1.e-09
!         f11vmr = f11(ntim)*1.e-12*(1.+cfcscl)
!         f12vmr = f12(ntim)*1.e-12
!         co2mmr = rmwco2 * co2vmr
!         return
   endif
!
! determine time interpolation factors, check sanity
! of interpolation factors to within 32-bit roundoff
! assume that day of year is 1 for all input data
!
   doymodel = yrmodel*365._r8    + calday
   doydatam = yrdata(nyrm)*365._r8 + 1._r8
   doydatap = yrdata(nyrp)*365._r8 + 1._r8

   call get_timeinterp_factors(.false.,2,doydatam,doydatap, doymodel, &
        fact1, fact2,'chem_surfvals')

!
! do time interpolation:
!   co2     in ppmv
!   n2o,ch4 in ppbv
!   f11,f12 in pptv
!
   co2vmr = (co2(nyrm)*fact1 + co2(nyrp)*fact2)*1.e-06_r8
   ch4vmr = (ch4(nyrm)*fact1 + ch4(nyrp)*fact2)*1.e-09_r8
   n2ovmr = (n2o(nyrm)*fact1 + n2o(nyrp)*fact2)*1.e-09_r8

   cfcscl = (adj(nyrm)*fact1 + adj(nyrp)*fact2)
   f11vmr = (f11(nyrm)*fact1 + f11(nyrp)*fact2)*1.e-12_r8*(1._r8+cfcscl)
   f12vmr = (f12(nyrm)*fact1 + f12(nyrp)*fact2)*1.e-12_r8

   return
end subroutine chem_surfvals_set_all

!=========================================================================================

subroutine chem_surfvals_set_co2()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes co2 greenhouse gas volume mixing ratio via ramping info 
! provided in namelist var's
! 
! Author: B. Eaton - updated ramp_ghg for use in chem_surfvals module
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8

!---------------------------Local variables-----------------------------

   real(r8) :: daydiff             ! number of days of co2 ramping
   integer  :: yr, mon, day, ncsec ! components of a date
   integer  :: ncdate              ! current date in integer format [yyyymmdd]
!-----------------------------------------------------------------------

   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day

   call timemgr_datediff(co2_start, 0, ncdate, ncsec, daydiff)

   if (daydiff > 0.0_r8) then

      co2vmr = co2_base*(co2_daily_factor)**daydiff

      if(co2_daily_factor < 1.0_r8) then
         co2vmr = max(co2vmr,co2_limit)
      else
         co2vmr = min(co2vmr,co2_limit)
      end if
   end if

   return
end subroutine chem_surfvals_set_co2


!=========================================================================================

end module chem_surfvals
