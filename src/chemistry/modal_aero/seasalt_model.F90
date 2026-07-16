!===============================================================================
! Seasalt for Modal Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,         only: pcols, pver
  use modal_aero_data,only: ntot_amode
  use perf_mod,       only: t_startf, t_stopf
  use ap_seasalt_emis_scheme, only: seasalt_emis_run

  implicit none
  private

  public :: seasalt_nbin
  public :: seasalt_nnum
  public :: seasalt_names
  public :: seasalt_indices
  public :: seasalt_init
  public :: seasalt_emis
  public :: seasalt_active

  integer, parameter :: nslt = max(3,ntot_amode-3)
  integer, parameter :: nnum = nslt
  integer, parameter :: seasalt_nbin = nslt
  integer, parameter :: seasalt_nnum = nnum

#if  ( defined MODAL_AERO_7MODE )
  character(len=6),parameter :: seasalt_names(nslt+nnum) = &
       (/ 'ncl_a1', 'ncl_a2', 'ncl_a4', 'ncl_a6', 'num_a1', 'num_a2', 'num_a4', 'num_a6' /)
#elif( defined MODAL_AERO_3MODE || defined MODAL_AERO_4MODE )
  character(len=6),parameter :: seasalt_names(nslt+nnum) = &
       (/ 'ncl_a1', 'ncl_a2', 'ncl_a3', 'num_a1', 'num_a2', 'num_a3'/)
#endif

  integer :: seasalt_indices(nslt+nnum)

  logical :: seasalt_active = .false.

contains
  
  !=============================================================================
  !=============================================================================
  subroutine seasalt_init
    use ap_sslt_sections_private, only: ap_sslt_sections_init
    use constituents,  only: cnst_get_ind

    integer :: m

    do m = 1, seasalt_nbin
       call cnst_get_ind(seasalt_names(m), seasalt_indices(m),abort=.false.)
    enddo
    do m = 1, seasalt_nnum
       call cnst_get_ind(seasalt_names(seasalt_nbin+m), seasalt_indices(seasalt_nbin+m),abort=.false.)
    enddo

    seasalt_active = any(seasalt_indices(:) > 0)

    if (.not.seasalt_active) return

    call ap_sslt_sections_init()

  end subroutine seasalt_init

  !=============================================================================
  !=============================================================================
  subroutine seasalt_emis( u10cubed,  srf_temp, ocnfrc, ncol, cflx )

    use mo_constants,  only: dns_aer_sst=>seasalt_density, pi

    ! dummy arguments
    real(r8), intent(in) :: u10cubed(:)
    real(r8), intent(in) :: srf_temp(:)
    real(r8), intent(in) :: ocnfrc(:)
    integer,  intent(in) :: ncol
    real(r8), intent(inout) :: cflx(:,:)
    integer :: errflg
    character(len=512) :: errmsg

#if  ( defined MODAL_AERO_7MODE )
    real(r8), parameter :: emis_scale = 1.62_r8
    real(r8), parameter :: sst_sz_range_lo (nslt) = (/ 0.08e-6_r8, 0.02e-6_r8, 0.3e-6_r8,  1.0e-6_r8 /)  ! accu, aitken, fine, coarse
    real(r8), parameter :: sst_sz_range_hi (nslt) = (/ 0.3e-6_r8,  0.08e-6_r8, 1.0e-6_r8, 10.0e-6_r8 /)
#elif( defined MODAL_AERO_3MODE || defined MODAL_AERO_4MODE )
    real(r8), parameter :: emis_scale = 1.35_r8 
    real(r8), parameter :: sst_sz_range_lo (nslt) =  (/ 0.08e-6_r8,  0.02e-6_r8,  1.0e-6_r8 /)  ! accu, aitken, coarse
    real(r8), parameter :: sst_sz_range_hi (nslt) =  (/ 1.0e-6_r8,   0.08e-6_r8, 10.0e-6_r8 /)
#endif

    call t_startf('ap_seasalt_emis_run')
    call seasalt_emis_run(ncol, size(cflx,1), size(cflx,2), seasalt_nbin, &
         size(seasalt_indices), u10cubed, srf_temp, ocnfrc, cflx, &
         seasalt_indices, sst_sz_range_lo, sst_sz_range_hi, emis_scale, &
         pi, dns_aer_sst, errmsg, errflg)
    call t_stopf('ap_seasalt_emis_run')

  end subroutine seasalt_emis

end module seasalt_model
