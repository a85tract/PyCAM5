! This module is used to diagnose the location of the tropopause. Multiple
! algorithms are provided, some of which may not be able to identify a
! tropopause in all situations. To handle these cases, an analytic
! definition and a climatology are provided that can be used to fill in
! when the original algorithm fails. The tropopause temperature and
! pressure are determined and can be output to the history file.
!
! These routines are based upon code in the WACCM chemistry module
! including mo_tropoause.F90 and llnl_set_chem_trop.F90. The code
! for the Reichler et al. [2003] algorithm is from:
!
!   http://www.gfdl.noaa.gov/~tjr/TROPO/tropocode.htm
!
! Author: Charles Bardeen
! Created: April, 2009

module tropopause
  !---------------------------------------------------------------
  ! ... variables for the tropopause module
  !---------------------------------------------------------------

  use shr_kind_mod,         only : r8 => shr_kind_r8
  use ppgrid,               only : pcols, pver, begchunk, endchunk
  use cam_abortutils,       only : endrun
  use cam_logfile,          only : iulog
  use cam_history_support,  only : fillvalue
  use physics_types,        only : physics_state
  use physconst,            only : cappa, rair, gravit, pi
  use spmd_utils,           only : masterproc
  use tropopause_find_scheme, only : tropopause_find_init, tropopause_find_run

  implicit none

  private
  
  public  :: tropopause_readnl, tropopause_init, tropopause_find, tropopause_output
  public  :: TROP_ALG_NONE, TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  public  :: TROP_ALG_STOBIE, TROP_ALG_HYBSTOB, TROP_ALG_TWMO, TROP_ALG_WMO

  save

  ! These parameters define and enumeration to be used to define the primary
  ! and backup algorithms to be used with the tropopause_find() method. The
  ! backup algorithm is meant to provide a solution when the primary algorithm
  ! fail. The algorithms that can't fail are: TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  ! and TROP_ALG_STOBIE.
  integer, parameter    :: TROP_ALG_NONE      = 1    ! Don't evaluate
  integer, parameter    :: TROP_ALG_ANALYTIC  = 2    ! Analytic Expression
  integer, parameter    :: TROP_ALG_CLIMATE   = 3    ! Climatology
  integer, parameter    :: TROP_ALG_STOBIE    = 4    ! Stobie Algorithm
  integer, parameter    :: TROP_ALG_TWMO      = 5    ! WMO Definition, Reichler et al. [2003]
  integer, parameter    :: TROP_ALG_WMO       = 6    ! WMO Definition
  integer, parameter    :: TROP_ALG_HYBSTOB   = 7    ! Hybrid Stobie Algorithm
  
  integer, parameter    :: TROP_NALG          = 7    ! Number of Algorithms  
  character,parameter   :: TROP_LETTER(TROP_NALG) = (/ ' ', 'A', 'C', 'S', 'T', 'W', 'H' /)
                                                     ! unique identifier for output, don't use P

  ! These variables should probably be controlled by namelist entries.
  logical ,parameter    :: output_all         = .False.              ! output tropopause info from all algorithms
  integer ,parameter    :: default_primary    = TROP_ALG_TWMO        ! default primary algorithm
  integer ,parameter    :: default_backup     = TROP_ALG_CLIMATE     ! default backup algorithm

  ! Namelist variables
  character(len=256)    :: tropopause_climo_file = 'trop_climo'      ! absolute filepath of climatology file

  ! These variables are used to store the climatology data.
  real(r8)              :: days(12)                                  ! days in the climatology
  real(r8), pointer     :: tropp_p_loc(:,:,:)                        ! climatological tropopause pressures

  integer, parameter :: NOTFOUND = -1

  real(r8),parameter :: ALPHA  = 0.03_r8
    
  ! physical constants
  ! These constants are set in module variables rather than as parameters 
  ! to support the aquaplanet mode in which the constants have values determined
  ! by the experiment protocol
  real(r8) :: cnst_kap     ! = cappa
  real(r8) :: cnst_faktor  ! = -gravit/rair
  real(r8) :: cnst_ka1     ! = cnst_kap - 1._r8

!================================================================================================
contains
!================================================================================================

   ! Read namelist variables.
   subroutine tropopause_readnl(nlfile)

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use mpishorthand

      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer :: unitn, ierr
      character(len=*), parameter :: subname = 'tropopause_readnl'

      namelist /tropopause_nl/ tropopause_climo_file
      !-----------------------------------------------------------------------------

      if (masterproc) then
         unitn = getunit()
         open( unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'tropopause_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, tropopause_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(tropopause_climo_file, len(tropopause_climo_file), mpichar, 0, mpicom)
#endif

   end subroutine tropopause_readnl


  ! This routine is called during intialization and must be called before the
  ! other methods in this module can be used. Its main tasks are to read in the
  ! climatology from a file and to define the output fields. Much of this code
  ! is taken from mo_tropopause.
  subroutine tropopause_init()
  

    use ppgrid,        only: pver
    use cam_pio_utils, only: phys_decomp
    use cam_history,   only: addfld


    implicit none

    character(len=512) :: errmsg
    integer :: errflg

    ! define physical constants
    cnst_kap    = cappa
    cnst_faktor = -gravit/rair
    cnst_ka1    = cnst_kap - 1._r8

    call tropopause_find_init(cappa, rair, gravit, pi, errmsg, errflg)
    if (errflg /= 0) call endrun(trim(errmsg))

    ! Define the output fields.
    call addfld('TROP_P',  'Pa',          1,    'A', 'Tropopause Pressure',    phys_decomp, flag_xyfill=.True.)
    call addfld('TROP_T',  'K',           1,    'A', 'Tropopause Temperature', phys_decomp, flag_xyfill=.True.)
    call addfld('TROP_Z',  'm',           1,    'A', 'Tropopause Height',      phys_decomp, flag_xyfill=.True.)
    call addfld('TROP_DZ', 'm',           pver, 'A', 'Relative Tropopause Height',  phys_decomp)
    call addfld('TROP_PD', 'probability', pver, 'A', 'Tropopause Probabilty',  phys_decomp)
    call addfld('TROP_FD', 'probability', 1,    'A', 'Tropopause Found',       phys_decomp)
    
    call addfld('TROPP_P',  'Pa',          1,    'A', 'Tropopause Pressure (primary)',     phys_decomp, flag_xyfill=.True.)
    call addfld('TROPP_T',  'K',           1,    'A', 'Tropopause Temperature (primary)',  phys_decomp, flag_xyfill=.True.)
    call addfld('TROPP_Z',  'm',           1,    'A', 'Tropopause Height (primary)',       phys_decomp, flag_xyfill=.True.)
    call addfld('TROPP_DZ', 'm',         pver,   'A', 'Relalive Tropopause Height (primary)',  phys_decomp)
    call addfld('TROPP_PD', 'probability', pver, 'A', 'Tropopause Distribution (primary)', phys_decomp)
    call addfld('TROPP_FD', 'probability', 1,    'A', 'Tropopause Found (primary)',        phys_decomp)
    
    call addfld( 'hstobie_trop',   'fraction of model time', pver, 'I', 'Lowest level with stratospheric chemsitry', phys_decomp )
    call addfld( 'hstobie_linoz',  'fraction of model time', pver, 'I', 'Lowest possible Linoz level', phys_decomp )
    call addfld( 'hstobie_tropop', 'fraction of model time', pver, 'I', &
         'Troposphere boundary calculated in chemistry', phys_decomp )

    ! If requested, be prepared to output results from all of the methods.
    if (output_all) then
      call addfld('TROPA_P',  'Pa',          1,    'A', 'Tropopause Pressure (analytic)',     phys_decomp, flag_xyfill=.True.)
      call addfld('TROPA_T',  'K',           1,    'A', 'Tropopause Temperature (analytic)',  phys_decomp, flag_xyfill=.True.)
      call addfld('TROPA_Z',  'm',           1,    'A', 'Tropopause Height (analytic)',       phys_decomp, flag_xyfill=.True.)
      call addfld('TROPA_PD', 'probability', pver, 'A', 'Tropopause Distribution (analytic)', phys_decomp)
      call addfld('TROPA_FD', 'probability', 1,    'A', 'Tropopause Found (analytic)',        phys_decomp)

      call addfld('TROPC_P',  'Pa',          1,    'A', 'Tropopause Pressure (climatology)',     phys_decomp, flag_xyfill=.True.)
      call addfld('TROPC_T',  'K',           1,    'A', 'Tropopause Temperature (climatology)',  phys_decomp, flag_xyfill=.True.)
      call addfld('TROPC_Z',  'm',           1,    'A', 'Tropopause Height (climatology)',       phys_decomp, flag_xyfill=.True.)
      call addfld('TROPC_PD', 'probability', pver, 'A', 'Tropopause Distribution (climatology)', phys_decomp)
      call addfld('TROPC_FD', 'probability', 1,    'A', 'Tropopause Found (climatology)',        phys_decomp)

      call addfld('TROPS_P',  'Pa',          1,    'A', 'Tropopause Pressure (stobie)',     phys_decomp, flag_xyfill=.True.)
      call addfld('TROPS_T',  'K',           1,    'A', 'Tropopause Temperature (stobie)',  phys_decomp, flag_xyfill=.True.)
      call addfld('TROPS_Z',  'm',           1,    'A', 'Tropopause Height (stobie)',       phys_decomp, flag_xyfill=.True.)
      call addfld('TROPS_PD', 'probability', pver, 'A', 'Tropopause Distribution (stobie)', phys_decomp)
      call addfld('TROPS_FD', 'probability', 1,    'A', 'Tropopause Found (stobie)',        phys_decomp)

      call addfld('TROPT_P',  'Pa',          1,    'A', 'Tropopause Pressure (twmo)',     phys_decomp, flag_xyfill=.True.)
      call addfld('TROPT_T',  'K',           1,    'A', 'Tropopause Temperature (twmo)',  phys_decomp, flag_xyfill=.True.)
      call addfld('TROPT_Z',  'm',           1,    'A', 'Tropopause Height (twmo)',       phys_decomp, flag_xyfill=.True.)
      call addfld('TROPT_PD', 'probability', pver, 'A', 'Tropopause Distribution (twmo)', phys_decomp)
      call addfld('TROPT_FD', 'probability', 1,    'A', 'Tropopause Found (twmo)',        phys_decomp)

      call addfld('TROPW_P',  'Pa',          1,    'A', 'Tropopause Pressure (WMO)',     phys_decomp, flag_xyfill=.True.)
      call addfld('TROPW_T',  'K',           1,    'A', 'Tropopause Temperature (WMO)',  phys_decomp, flag_xyfill=.True.)
      call addfld('TROPW_Z',  'm',           1,    'A', 'Tropopause Height (WMO)',       phys_decomp, flag_xyfill=.True.)
      call addfld('TROPW_PD', 'probability', pver, 'A', 'Tropopause Distribution (WMO)', phys_decomp)
      call addfld('TROPW_FD', 'probability', 1,    'A', 'Tropopause Found (WMO)',        phys_decomp)

      call addfld('TROPH_P',  'Pa',          1,    'A', 'Tropopause Pressure (Hybrid Stobie)',     phys_decomp, flag_xyfill=.True.)
      call addfld('TROPH_T',  'K',           1,    'A', 'Tropopause Temperature (Hybrid Stobie)',  phys_decomp, flag_xyfill=.True.)
      call addfld('TROPH_Z',  'm',           1,    'A', 'Tropopause Height (Hybrid Stobie)',       phys_decomp, flag_xyfill=.True.)
      call addfld('TROPH_PD', 'probability', pver, 'A', 'Tropopause Distribution (Hybrid Stobie)', phys_decomp)
      call addfld('TROPH_FD', 'probability', 1,    'A', 'Tropopause Found (Hybrid Stobie)',        phys_decomp)
     end if


    call tropopause_read_file()


  end subroutine tropopause_init
  

  subroutine tropopause_read_file
    !------------------------------------------------------------------
    ! ... initialize upper boundary values
    !------------------------------------------------------------------
    use interpolate_data,  only : lininterp_init, lininterp, interp_type, lininterp_finish
    use dyn_grid,     only : get_dyn_grid_parm
    use phys_grid,    only : get_ncols_p, get_rlat_all_p, get_rlon_all_p	
    use ioFileMod,    only : getfil
    use time_manager, only : get_calday
    use physconst,    only : pi
    use cam_pio_utils, only: cam_pio_openfile
    use pio,          only : file_desc_t, var_desc_t, pio_inq_dimid, pio_inq_dimlen, &
         pio_inq_varid, pio_get_var, pio_closefile, pio_nowrite

    !------------------------------------------------------------------
    ! ... local variables
    !------------------------------------------------------------------
    integer :: i, j, n
    integer :: ierr
    type(file_desc_t) :: pio_id
    integer :: dimid
    type(var_desc_t) :: vid
    integer :: nlon, nlat, ntimes
    integer :: start(3)
    integer :: count(3)
    integer, parameter :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
         716, 816, 915, 1016, 1115, 1216 /)
    integer :: plon, plat
    type(interp_type) :: lon_wgts, lat_wgts
    real(r8), allocatable :: tropp_p_in(:,:,:)
    real(r8), allocatable :: lat(:)
    real(r8), allocatable :: lon(:)
    real(r8) :: to_lats(pcols), to_lons(pcols)
    real(r8), parameter :: d2r=pi/180._r8, zero=0._r8, twopi=pi*2._r8
    character(len=256) :: locfn
    integer  :: c, ncols


    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')


    !-----------------------------------------------------------------------
    !       ... open netcdf file
    !-----------------------------------------------------------------------
    call getfil (tropopause_climo_file, locfn, 0)
    call cam_pio_openfile(pio_id, trim(locfn), PIO_NOWRITE)

    !-----------------------------------------------------------------------
    !       ... get time dimension
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'time', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, ntimes )
    if( ntimes /= 12 )then
       write(iulog,*) 'tropopause_init: number of months = ',ntimes,'; expecting 12'
       call endrun
    end if
    !-----------------------------------------------------------------------
    !       ... get latitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'lat', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, nlat )
    allocate( lat(nlat), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_init: lat allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( pio_id, 'lat', vid )
    ierr = pio_get_var( pio_id, vid, lat )
    lat(:nlat) = lat(:nlat) * d2r
    !-----------------------------------------------------------------------
    !       ... get longitudes
    !-----------------------------------------------------------------------
    ierr = pio_inq_dimid( pio_id, 'lon', dimid )
    ierr = pio_inq_dimlen( pio_id, dimid, nlon )
    allocate( lon(nlon), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_init: lon allocation error = ',ierr
       call endrun
    end if
    ierr = pio_inq_varid( pio_id, 'lon', vid )
    ierr = pio_get_var( pio_id, vid, lon )
    lon(:nlon) = lon(:nlon) * d2r

    !------------------------------------------------------------------
    !  ... allocate arrays
    !------------------------------------------------------------------
    allocate( tropp_p_in(nlon,nlat,ntimes), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'tropopause_init: tropp_p_in allocation error = ',ierr
       call endrun
    end if
    !------------------------------------------------------------------
    !  ... read in the tropopause pressure
    !------------------------------------------------------------------
    ierr = pio_inq_varid( pio_id, 'trop_p', vid )
    start = (/ 1, 1, 1 /)
    count = (/ nlon, nlat, ntimes /)
    ierr = pio_get_var( pio_id, vid, start, count, tropp_p_in )

    !------------------------------------------------------------------
    !  ... close the netcdf file
    !------------------------------------------------------------------
    call pio_closefile( pio_id )

    !--------------------------------------------------------------------
    !  ... regrid
    !--------------------------------------------------------------------

    allocate( tropp_p_loc(pcols,begchunk:endchunk,ntimes), stat=ierr )

    if( ierr /= 0 ) then
      write(iulog,*) 'tropopause_init: tropp_p_loc allocation error = ',ierr
      call endrun
    end if

    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)
       call lininterp_init(lon, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(lat, nlat, to_lats, ncols, 1, lat_wgts)
       do n=1,ntimes
          call lininterp(tropp_p_in(:,:,n), nlon, nlat, tropp_p_loc(1:ncols,c,n), ncols, lon_wgts, lat_wgts)    
       end do
       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)
    end do
    deallocate(lon)
    deallocate(lat)
    deallocate(tropp_p_in)

    !--------------------------------------------------------
    ! ... initialize the monthly day of year times
    !--------------------------------------------------------

    do n = 1,12
       days(n) = get_calday( dates(n), 0 )
    end do
    if (masterproc) then
       write(iulog,*) 'tropopause_init : days'
       write(iulog,'(1p,5g15.8)') days(:)
    endif

  end subroutine tropopause_read_file
  
  ! for the tropopause level, temperature and pressure.
  subroutine tropopause_find(pstate, tropLev, tropP, tropT, tropZ, primary, backup)

    use cam_history,  only : outfld
    use perf_mod,     only : t_startf, t_stopf
    use time_manager, only : get_curr_calday

    implicit none

    type(physics_state), intent(in)     :: pstate 
    integer, optional, intent(in)       :: primary                   ! primary detection algorithm
    integer, optional, intent(in)       :: backup                    ! backup detection algorithm
    integer,            intent(out)     :: tropLev(pcols)            ! tropopause level index   
    real(r8), optional, intent(out)     :: tropP(pcols)              ! tropopause pressure (Pa)  
    real(r8), optional, intent(out)     :: tropT(pcols)              ! tropopause temperature (K)
    real(r8), optional, intent(out)     :: tropZ(pcols)              ! tropopause height (m)
    
    ! Local Variables
    integer       :: primAlg            ! Primary algorithm
    integer       :: backAlg            ! Backup algorithm
    integer       :: ncol
    integer       :: lchnk
    integer       :: errflg
    real(r8)      :: calday
    real(r8)      :: hstobie_trop(pcols,pver)
    real(r8)      :: hstobie_linoz(pcols,pver)
    real(r8)      :: hstobie_tropop(pcols,pver)
    character(len=512) :: errmsg
    
    ! Set the algorithms to be used, either the ones provided or the defaults.
    if (present(primary)) then
      primAlg = primary
    else
      primAlg = default_primary
    end if
    
    if (present(backup)) then
      backAlg = backup
    else
      backAlg = default_backup
    end if
    
    ncol = pstate%ncol
    lchnk = pstate%lchnk
    calday = get_curr_calday()

    call t_startf('ap_tropopause_find_run')
    call tropopause_find_run(ncol, pver, fillvalue, pstate%lat, pstate%pint, &
         pstate%pmid, pstate%t, pstate%zi, pstate%zm, pstate%phis, calday, &
         tropp_p_loc(:,lchnk,:), days, tropLev, tropP=tropP, tropT=tropT, &
         tropZ=tropZ, hstobie_trop=hstobie_trop(:ncol,:), &
         hstobie_linoz=hstobie_linoz(:ncol,:), &
         hstobie_tropop=hstobie_tropop(:ncol,:), primary=primAlg, &
         backup=backAlg, errmsg=errmsg, errflg=errflg)
    call t_stopf('ap_tropopause_find_run')

    if (errflg /= 0) call endrun(trim(errmsg))

    ! CAM history ownership remains in the adapter.
    if ((primAlg == TROP_ALG_HYBSTOB) .or. (backAlg == TROP_ALG_HYBSTOB)) then
       call outfld('hstobie_trop',   hstobie_trop(:ncol,:),   ncol, lchnk)
       call outfld('hstobie_linoz',  hstobie_linoz(:ncol,:),  ncol, lchnk)
       call outfld('hstobie_tropop', hstobie_tropop(:ncol,:), ncol, lchnk)
    end if
  end subroutine tropopause_find
  
  
  ! Output the tropopause pressure and temperature to the history files. Two sets
  ! of output will be generated, one for the default algorithm and another one
  ! using the default routine, but backed by a climatology when the default
  ! algorithm fails.
  subroutine tropopause_output(pstate)
    use cam_history,  only : outfld
    
    implicit none

    type(physics_state), intent(in)     :: pstate
  
    ! Local Variables
    integer       :: i
    integer       :: alg
    integer       :: ncol                     ! number of cloumns in the chunk
    integer       :: lchnk                    ! chunk identifier
    integer       :: tropLev(pcols)           ! tropopause level index   
    real(r8)      :: tropP(pcols)             ! tropopause pressure (Pa)  
    real(r8)      :: tropT(pcols)             ! tropopause temperature (K) 
    real(r8)      :: tropZ(pcols)             ! tropopause height (m) 
    real(r8)      :: tropFound(pcols)         ! tropopause found  
    real(r8)      :: tropDZ(pcols, pver)      ! relative tropopause height (m) 
    real(r8)      :: tropPdf(pcols, pver)     ! tropopause probability distribution  

    ! Information about the chunk.  
    lchnk = pstate%lchnk
    ncol  = pstate%ncol

    ! Find the tropopause using the default algorithm backed by the climatology.
    call tropopause_find(pstate, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ)
    
    tropPdf(:,:) = 0._r8
    tropFound(:) = 0._r8
    tropDZ(:,:) = fillvalue 
    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._r8
        tropFound(i) = 1._r8
        tropDZ(i,:) = pstate%zm(i,:) - tropZ(i) 
      end if
    end do

    call outfld('TROP_P',   tropP(:ncol),      ncol, lchnk)
    call outfld('TROP_T',   tropT(:ncol),      ncol, lchnk)
    call outfld('TROP_Z',   tropZ(:ncol),      ncol, lchnk)
    call outfld('TROP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    call outfld('TROP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    call outfld('TROP_FD',  tropFound(:ncol),  ncol, lchnk)
    
    
    ! Find the tropopause using just the primary algorithm.
    call tropopause_find(pstate, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, backup=TROP_ALG_NONE)

    tropPdf(:,:) = 0._r8
    tropFound(:) = 0._r8
    tropDZ(:,:) = fillvalue 
    
    do i = 1, ncol
      if (tropLev(i) /= NOTFOUND) then
        tropPdf(i, tropLev(i)) = 1._r8
        tropFound(i) = 1._r8
        tropDZ(i,:) = pstate%zm(i,:) - tropZ(i) 
      end if
    end do

    call outfld('TROPP_P',   tropP(:ncol),      ncol, lchnk)
    call outfld('TROPP_T',   tropT(:ncol),      ncol, lchnk)
    call outfld('TROPP_Z',   tropZ(:ncol),      ncol, lchnk)
    call outfld('TROPP_DZ',  tropDZ(:ncol, :), ncol, lchnk)
    call outfld('TROPP_PD',  tropPdf(:ncol, :), ncol, lchnk)
    call outfld('TROPP_FD',  tropFound(:ncol),  ncol, lchnk)
    
    
    ! If requested, do all of the algorithms.
    if (output_all) then
    
      do alg = 2, TROP_NALG
    
        ! Find the tropopause using just the analytic algorithm.
        call tropopause_find(pstate, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, primary=alg, backup=TROP_ALG_NONE)
  
        tropPdf(:,:) = 0._r8
        tropFound(:) = 0._r8
      
        do i = 1, ncol
          if (tropLev(i) /= NOTFOUND) then
            tropPdf(i, tropLev(i)) = 1._r8
            tropFound(i) = 1._r8
          end if
        end do
  
        call outfld('TROP' // TROP_LETTER(alg) // '_P',   tropP(:ncol),      ncol, lchnk)
        call outfld('TROP' // TROP_LETTER(alg) // '_T',   tropT(:ncol),      ncol, lchnk)
        call outfld('TROP' // TROP_LETTER(alg) // '_Z',   tropZ(:ncol),      ncol, lchnk)
        call outfld('TROP' // TROP_LETTER(alg) // '_PD',  tropPdf(:ncol, :), ncol, lchnk)
        call outfld('TROP' // TROP_LETTER(alg) // '_FD',  tropFound(:ncol),  ncol, lchnk)
      end do
    end if
    
    return
  end subroutine tropopause_output
end module tropopause
