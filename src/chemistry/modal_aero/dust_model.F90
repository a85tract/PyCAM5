!===============================================================================
! Dust for Modal Aerosol Model
!===============================================================================
module dust_model 
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun
  use mo_util,          only: chemistry_misc_codon_touch

  implicit none
  private

  public :: dust_names
  public :: dust_nbin
  public :: dust_nnum
  public :: dust_indices
  public :: dust_emis_sclfctr
  public :: dust_dmt_vwr
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
  real(r8), target :: dust_dmt_vwr(dust_nbin)
  real(r8) :: dust_stk_crc(dust_nbin)

  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  logical :: dust_active = .false.
  logical :: dust_emis_codon_logged = .false.

  interface
     subroutine dust_emis_codon(ncol_c, pcols_c, ndstflx_c, dust_nbin_c, soil_erod_fact_c, pi_c, &
                                dust_density_c, soil_erod_threshold_c, dust_flux_in_p, cflx_p, soil_erod_p, &
                                soil_erodibility_p, dust_indices_p, dust_emis_sclfctr_p, dust_dmt_vwr_p) &
                                bind(c, name="dust_emis_codon")
       use iso_c_binding, only : c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, ndstflx_c, dust_nbin_c
       real(c_double), value :: soil_erod_fact_c, pi_c, dust_density_c, soil_erod_threshold_c
       type(c_ptr), value :: dust_flux_in_p, cflx_p, soil_erod_p, soil_erodibility_p
       type(c_ptr), value :: dust_indices_p, dust_emis_sclfctr_p, dust_dmt_vwr_p
     end subroutine dust_emis_codon
  end interface

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

    call chemistry_misc_codon_touch('dust_model', 127)

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
    use ppgrid,        only : pcols
    use cam_logfile,   only : iulog
    use iso_c_binding, only : c_double, c_int64_t, c_loc

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in), target    :: dust_flux_in(:,:)
    real(r8), intent(inout), target :: cflx(:,:)
    real(r8), intent(out), target   :: soil_erod(:)

  ! local vars
    integer :: i, m, idst, inum
    integer :: status, n, code, ndstflx
    character(len=32) :: impl_name
    logical :: use_native_impl
    integer(c_int64_t), target :: dust_indices_c(dust_nbin+dust_nnum)
    real(c_double), target :: dust_emis_sclfctr_c(dust_nbin)
    real(c_double), target :: dust_dmt_vwr_c(dust_nbin)
    real(r8) :: x_mton
    real(r8),parameter :: soil_erod_threshold = 0.1_r8

    impl_name = 'codon'
    call cam_codon_get_impl('DUST_EMIS_IMPL', impl_name, n, status)
    use_native_impl = .false.
    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    end if

    if (.not. use_native_impl) then
       ndstflx = size(dust_flux_in, 2)
       dust_indices_c(:) = int(dust_indices(:), c_int64_t)
       dust_emis_sclfctr_c(:) = real(dust_emis_sclfctr(:), c_double)
       dust_dmt_vwr_c(:) = real(dust_dmt_vwr(:), c_double)

       if (masterproc .and. .not. dust_emis_codon_logged) then
          write(iulog,'(A)') 'dust_emis implementation = codon'
          call flush(iulog)
          dust_emis_codon_logged = .true.
       end if

       call dust_emis_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(ndstflx, c_int64_t), int(dust_nbin, c_int64_t), &
            real(soil_erod_fact, c_double), real(pi, c_double), real(dust_density, c_double), &
            real(soil_erod_threshold, c_double), c_loc(dust_flux_in(1,1)), c_loc(cflx(1,1)), c_loc(soil_erod(1)), &
            c_loc(soil_erodibility(1,lchnk)), c_loc(dust_indices_c(1)), c_loc(dust_emis_sclfctr_c(1)), &
            c_loc(dust_dmt_vwr_c(1)) &
       )
       return
    end if

    ! set dust emissions

    col_loop: do i =1,ncol

       soil_erod(i) = soil_erodibility( i, lchnk )

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       ! rebin and adjust dust emissons..
       do m = 1,dust_nbin

          idst = dust_indices(m)

          cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
               * dust_emis_sclfctr(m)*soil_erod(i)/soil_erod_fact*1.15_r8

          x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))                

          inum = dust_indices(m+dust_nbin)

          cflx(i,inum) = cflx(i,idst)*x_mton

       enddo

    end do col_loop

  end subroutine dust_emis

end module dust_model
