module mo_srf_emissions
  !---------------------------------------------------------------
  ! 	... surface emissions module
  !---------------------------------------------------------------

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use chem_mods,     only : gas_pcnst
  use spmd_utils,    only : masterproc,iam
  use mo_tracname,   only : solsym
  use cam_abortutils,only : endrun
  use ioFileMod,     only : getfil
  use ppgrid,        only : pcols, begchunk, endchunk
  use cam_logfile,   only : iulog
  use tracer_data,   only : trfld,trfile
  use perf_mod,      only : t_startf, t_stopf
  use ap_set_srf_emissions_scheme, only : set_srf_emissions_run

  implicit none

  type :: emission
     integer           :: spc_ndx
     real(r8)          :: mw
     real(r8)          :: scalefactor
     character(len=256):: filename
     character(len=16) :: species
     character(len=8)  :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type emission

  private

  public  :: srf_emissions_inti, set_srf_emissions, set_srf_emissions_time 

  save

  real(r8), parameter :: amufac = 1.65979e-23_r8         ! 1.e4* kg / amu
  logical :: has_emis(gas_pcnst)
  type(emission), allocatable :: emissions(:)
  integer                     :: n_emis_files 
  integer :: c10h16_ndx, isop_ndx

contains

  subroutine srf_emissions_inti( srf_emis_specifier, emis_type_in, emis_cycle_yr, emis_fixed_ymd, emis_fixed_tod )

    !-----------------------------------------------------------------------
    ! 	... initialize the surface emissions
    !-----------------------------------------------------------------------

    use chem_mods,        only : adv_mass
    use mo_constants,     only : d2r, pi, rearth
    use string_utils,     only : to_upper
    use mo_chem_utls,     only : get_spc_ndx 
    use tracer_data,      only : trcdata_init
    use cam_pio_utils,    only : cam_pio_openfile
    use pio,              only : pio_inquire, pio_nowrite, pio_closefile, pio_inq_varndims
    use pio,              only : pio_inq_varname, file_desc_t, pio_get_att, PIO_NOERR, PIO_GLOBAL
    use pio,              only : pio_seterrorhandling, PIO_BCAST_ERROR,PIO_INTERNAL_ERROR
    use chem_surfvals,    only : flbc_list
    use string_utils,     only : GLC

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: srf_emis_specifier(:)
    character(len=*), intent(in) :: emis_type_in
    integer,          intent(in) :: emis_cycle_yr
    integer,          intent(in) :: emis_fixed_ymd
    integer,          intent(in) :: emis_fixed_tod

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer  :: astat
    integer  :: j, l, m, n, i, nn                     ! Indices
    character(len=16)  :: spc_name
    character(len=256) :: filename

    character(len=16)  :: emis_species(gas_pcnst)
    character(len=256) :: emis_filenam(gas_pcnst)
    integer  :: emis_indexes(gas_pcnst)
    real(r8) :: emis_scalefactor(gas_pcnst)

    integer :: vid, nvars, isec
    integer, allocatable :: vndims(:)
    type(file_desc_t) :: ncid
    character(len=32)  :: varname
    character(len=256) :: locfn
    integer :: ierr
    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    logical         , parameter :: rmv_file = .false.

    character(len=32) :: emis_type = ' '
    character(len=80) :: file_interp_type = ' '
    character(len=256) :: tmp_string = ' '
    character(len=32) :: xchr = ' '
    real(r8) :: xdbl

    has_emis(:) = .false.
    nn = 0

    count_emis: do n=1,gas_pcnst
       if ( len_trim(srf_emis_specifier(n) ) == 0 ) then
          exit count_emis
       endif

       i = scan(srf_emis_specifier(n),'->')
       spc_name = trim(adjustl(srf_emis_specifier(n)(:i-1)))
       
       ! need to parse out scalefactor ...
       tmp_string = adjustl(srf_emis_specifier(n)(i+2:))
       j = scan( tmp_string, '*' )
       if (j>0) then
          xchr = tmp_string(1:j-1) ! get the multipler (left of the '*')
          read( xchr, * ) xdbl   ! convert the string to a real
          tmp_string = adjustl(tmp_string(j+1:)) ! get the filepath name (right of the '*')
       else
          xdbl = 1._r8
       endif
       filename = trim(tmp_string)

       m = get_spc_ndx(spc_name)

       if (m > 0) then
          has_emis(m) = .true.
          has_emis(m) = has_emis(m) .and. ( .not. any( flbc_list == spc_name ) )
       else 
          write(iulog,*) 'srf_emis_inti: spc_name ',spc_name,' is not included in the simulation'
          call endrun('srf_emis_inti: invalid surface emission specification')
       endif

       if ( has_emis(m) ) then
          nn = nn+1
          emis_species(nn) = spc_name
          emis_filenam(nn) = filename
          emis_indexes(nn) = m
          emis_scalefactor(nn) = xdbl
       endif

    enddo count_emis

    n_emis_files = nn

    if (masterproc) write(iulog,*) 'srf_emis_inti: n_emis_files = ',n_emis_files

    allocate( emissions(n_emis_files), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'srf_emis_inti: failed to allocate emissions array; error = ',astat
       call endrun
    end if

    !-----------------------------------------------------------------------
    ! 	... setup the emission type array
    !-----------------------------------------------------------------------
    do m=1,n_emis_files 
       emissions(m)%spc_ndx          = emis_indexes(m)
       emissions(m)%units            = 'Tg/y'
       emissions(m)%species          = emis_species(m)
       emissions(m)%mw               = adv_mass(emis_indexes(m))                     ! g / mole
       emissions(m)%filename         = emis_filenam(m)
       emissions(m)%scalefactor      = emis_scalefactor(m)
    enddo

    !-----------------------------------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------------------------------
    spc_loop: do m = 1, n_emis_files

       emissions(m)%nsectors = 0
       
       call getfil (emissions(m)%filename, locfn, 0)
       call cam_pio_openfile ( ncid, trim(locfn), PIO_NOWRITE)
       ierr = pio_inquire (ncid, nvariables=nvars)

       allocate(vndims(nvars))

       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, vndims(vid))

          if( vndims(vid) < 3 ) then
             cycle
          elseif( vndims(vid) > 3 ) then
             ierr = pio_inq_varname (ncid, vid, varname)
             write(iulog,*) 'srf_emis_inti: Skipping variable ', trim(varname),', ndims = ',vndims(vid), &
                  ' , species=',trim(emissions(m)%species)
             cycle
          end if

          emissions(m)%nsectors = emissions(m)%nsectors+1

       enddo

       allocate( emissions(m)%sectors(emissions(m)%nsectors), stat=astat )
       if( astat/= 0 ) then
         write(iulog,*) 'srf_emis_inti: failed to allocate emissions(m)%sectors array; error = ',astat
         call endrun
       end if

       isec = 1

       do vid = 1,nvars
          if( vndims(vid) == 3 ) then
             ierr = pio_inq_varname(ncid, vid, emissions(m)%sectors(isec))
             isec = isec+1
          endif

       enddo
       deallocate(vndims)

       ! Global attribute 'input_method' overrides the srf_emis_type namelist setting on
       ! a file-by-file basis.  If the emis file does not contain the 'input_method' 
       ! attribute then the srf_emis_type namelist setting is used.
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ierr = pio_get_att(ncid, PIO_GLOBAL, 'input_method', file_interp_type)
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
       if ( ierr == PIO_NOERR) then
          l = GLC(file_interp_type)
          emis_type(1:l) = file_interp_type(1:l)
          emis_type(l+1:) = ' '
       else
          emis_type = trim(emis_type_in)
       endif

       call pio_closefile (ncid)

       allocate(emissions(m)%file%in_pbuf(size(emissions(m)%sectors)))
       emissions(m)%file%in_pbuf(:) = .false.

       call trcdata_init( emissions(m)%sectors, &
                          emissions(m)%filename, filelist, datapath, &
                          emissions(m)%fields,  &
                          emissions(m)%file, &
                          rmv_file, emis_cycle_yr, emis_fixed_ymd, emis_fixed_tod, trim(emis_type) )

    enddo spc_loop

    c10h16_ndx = get_spc_ndx('C10H16')
    isop_ndx = get_spc_ndx('ISOP')

  end subroutine srf_emissions_inti

  subroutine set_srf_emissions_time( pbuf2d, state )
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use tracer_data,  only : advance_trcdata
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    do m = 1,n_emis_files
       call advance_trcdata( emissions(m)%fields, emissions(m)%file, state, pbuf2d  )
    end do

  end subroutine set_srf_emissions_time

  ! adds surf flux specified in file to sflx
  subroutine set_srf_emissions( lchnk, ncol, sflx )
    !--------------------------------------------------------
    !	... form the surface fluxes for this latitude slice
    !--------------------------------------------------------

    use mo_constants, only : pi
    use time_manager, only : get_curr_calday
    use string_utils, only : to_lower, GLC
    use phys_grid,    only : get_rlat_all_p, get_rlon_all_p

    implicit none

    !--------------------------------------------------------
    !	... Dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)  :: ncol                  ! columns in chunk
    integer,  intent(in)  :: lchnk                 ! chunk index
    real(r8), intent(out) :: sflx(:,:) ! surface emissions ( kg/m^2/s )

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer :: m, isec
    integer :: max_sector_count, allocated_emission_files
    integer :: errflg
    real(r8) :: calday
    real(r8), parameter :: twopi = 2.0_r8 * pi
    real(r8), parameter :: pid2 = 0.5_r8 * pi
    real(r8), parameter :: dec_max = 23.45_r8 * pi/180._r8
    real(r8), allocatable :: sector_flux(:,:,:)
    real(r8), allocatable :: emission_scale(:), molecular_weight(:)
    integer, allocatable :: species_index(:), sector_count(:)
    logical, allocatable :: units_are_mks(:)
    logical :: has_c10h16_emission, has_isoprene_emission
    character(len=512) :: errmsg

    character(len=12),parameter :: mks_units(4) = (/ "kg/m2/s     ", &
                                                     "kg/m2/sec   ", &
                                                     "kg/m^2/s    ", &
                                                     "kg/m^2/sec  " /)
    character(len=12) :: units

    real(r8), dimension(ncol) :: rlats, rlons

    allocated_emission_files = max(1,n_emis_files)
    max_sector_count = 1
    do m = 1,n_emis_files
       max_sector_count = max(max_sector_count,emissions(m)%nsectors)
    end do

    allocate(sector_flux(size(sflx,1),max_sector_count,allocated_emission_files))
    allocate(emission_scale(allocated_emission_files))
    allocate(molecular_weight(allocated_emission_files))
    allocate(species_index(allocated_emission_files))
    allocate(sector_count(allocated_emission_files))
    allocate(units_are_mks(allocated_emission_files))

    sector_flux = 0._r8
    emission_scale = 0._r8
    molecular_weight = 0._r8
    species_index = 1
    sector_count = 0
    units_are_mks = .false.

    ! Keep tracer-data pointers and unit-string handling on the CAM side.
    do m = 1,n_emis_files
       species_index(m) = emissions(m)%spc_ndx
       sector_count(m) = emissions(m)%nsectors
       emission_scale(m) = emissions(m)%scalefactor
       molecular_weight(m) = emissions(m)%mw
       do isec = 1,emissions(m)%nsectors
          sector_flux(:ncol,isec,m) = emissions(m)%fields(isec)%data(:ncol,1,lchnk)
       end do

       units = to_lower(trim(emissions(m)%fields(1)%units(:GLC(emissions(m)%fields(1)%units))))
       units_are_mks(m) = any(mks_units(:) == units)
    end do

    call get_rlat_all_p( lchnk, ncol, rlats )
    call get_rlon_all_p( lchnk, ncol, rlons )

    calday = get_curr_calday()

    has_c10h16_emission = .false.
    if (c10h16_ndx > 0) has_c10h16_emission = has_emis(c10h16_ndx)
    has_isoprene_emission = .false.
    if (isop_ndx > 0) has_isoprene_emission = has_emis(isop_ndx)

    call t_startf('ap_set_srf_emissions_run')
    call set_srf_emissions_run(ncol, size(sflx,1), size(sflx,2), &
         n_emis_files, allocated_emission_files, max_sector_count, &
         species_index, sector_count, sector_flux, emission_scale, &
         molecular_weight, units_are_mks, amufac, rlats, rlons, calday, &
         pi, twopi, pid2, dec_max, c10h16_ndx, has_c10h16_emission, isop_ndx, &
         has_isoprene_emission, sflx, errmsg, errflg)
    call t_stopf('ap_set_srf_emissions_run')

    deallocate(sector_flux, emission_scale, molecular_weight)
    deallocate(species_index, sector_count, units_are_mks)

  end subroutine set_srf_emissions

end module mo_srf_emissions
