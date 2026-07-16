!===============================================================================
! Dust for Modal Aerosol Model
!===============================================================================
module dust_model 
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun
  use perf_mod,         only: t_startf, t_stopf
  use ap_dust_emis_scheme, only: dust_emis_run

  implicit none
  private

  public :: dust_names
  public :: dust_nbin
  public :: dust_nnum
  public :: dust_indices
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init
  public :: dust_active

#if  ( defined MODAL_AERO_3MODE )
  integer, parameter :: dust_nbin = 2
  integer, parameter :: dust_nnum = 2
  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = (/ 'dst_a1', 'dst_a3', 'num_a1', 'num_a3' /)
  real(r8),         parameter :: dust_dmt_grd(dust_nbin+1) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.032_r8,0.968_r8 /)
#elif ( defined MODAL_AERO_4MODE )
  integer, parameter :: dust_nbin = 3
  integer, parameter :: dust_nnum = 3
  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = &
                                 (/ 'dst_a2', 'dst_a1', 'dst_a3', 'num_a2', 'num_a1', 'num_a3' /) ! Aitken dust
  real(r8),         parameter :: dust_dmt_grd(dust_nbin+1) = &
                                 (/ 0.01e-6_r8, 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8 /) ! Aitken dust
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = &
                                 (/ 1.65E-05_r8, 0.011_r8, 0.999_r8 /) ! Aitken dust
#elif ( defined MODAL_AERO_7MODE )
  integer, parameter :: dust_nbin = 2
  integer, parameter :: dust_nnum = 2
  character(len=6), parameter :: dust_names(dust_nbin+dust_nnum) = (/ 'dst_a5', 'dst_a7', 'num_a5', 'num_a7' /)
  real(r8),         parameter :: dust_dmt_grd(dust_nbin+1) = (/ 0.1e-6_r8, 2.0e-6_r8, 10.0e-6_r8/)
  real(r8),         parameter :: dust_emis_sclfctr(dust_nbin) = (/ 0.13_r8, 0.87_r8 /)
#endif

  integer  :: dust_indices(dust_nbin+dust_nnum)
  real(r8) :: dust_dmt_vwr(dust_nbin)
  real(r8) :: dust_stk_crc(dust_nbin)

  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  logical :: dust_active = .false.

 contains

  !=============================================================================
  ! reads dust namelist options
  !=============================================================================
  subroutine dust_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(dust_emis_fact, 1,                   mpir8,   0, mpicom)
    call mpibcast(soil_erod_file, len(soil_erod_file), mpichar, 0, mpicom)
#endif

  end subroutine dust_readnl

  !=============================================================================
  !=============================================================================
  subroutine dust_init()
    use soil_erod_mod, only: soil_erod_init
    use constituents,  only: cnst_get_ind
    use dust_common,   only: dust_set_params

    integer :: n

    do n = 1, dust_nbin
       call cnst_get_ind(dust_names(n), dust_indices(n),abort=.false.)
    end do
    do n = 1, dust_nnum
       call cnst_get_ind(dust_names(dust_nbin+n), dust_indices(dust_nbin+n),abort=.false.)
    enddo 
    dust_active = any(dust_indices(:) > 0)
    if (.not.dust_active) return
   
    call  soil_erod_init( dust_emis_fact, soil_erod_file )

    call dust_set_params( dust_nbin, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )

  end subroutine dust_init

  !===============================================================================
  !===============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility
    use mo_constants,  only : dust_density
    use physconst,     only : pi

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

  ! local vars
    integer :: i
    integer :: errflg
    character(len=512) :: errmsg

    ! set dust emissions

    do i = 1,ncol
       soil_erod(i) = soil_erodibility( i, lchnk )
    end do

    call t_startf('ap_dust_emis_run')
    call dust_emis_run(ncol, size(cflx,1), size(dust_flux_in,2), size(cflx,2), &
         dust_nbin, size(dust_indices), dust_flux_in, cflx, soil_erod, &
         dust_indices, dust_emis_sclfctr, dust_dmt_vwr, soil_erod_fact, &
         pi, dust_density, errmsg, errflg)
    call t_stopf('ap_dust_emis_run')

  end subroutine dust_emis

end module dust_model
