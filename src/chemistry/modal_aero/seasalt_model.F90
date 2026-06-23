!===============================================================================
! Seasalt for Modal Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,         only: pcols, pver
  use modal_aero_data,only: ntot_amode
  use mo_util,        only: chemistry_misc_codon_touch

  implicit none
  private

  public :: seasalt_nbin
  public :: seasalt_nnum
  public :: seasalt_names
  public :: seasalt_indices
  public :: seasalt_emis_scale
  public :: seasalt_sz_range_lo
  public :: seasalt_sz_range_hi
  public :: seasalt_init
  public :: seasalt_emis
  public :: seasalt_active

  integer, parameter :: nslt = max(3,ntot_amode-3)
  integer, parameter :: nnum = nslt
  integer, parameter :: seasalt_nbin = nslt
  integer, parameter :: seasalt_nnum = nnum

#if  ( defined MODAL_AERO_7MODE )
  real(r8), parameter :: seasalt_emis_scale = 1.62_r8
  real(r8), parameter :: seasalt_sz_range_lo(nslt) = (/ 0.08e-6_r8, 0.02e-6_r8, 0.3e-6_r8,  1.0e-6_r8 /)
  real(r8), parameter :: seasalt_sz_range_hi(nslt) = (/ 0.3e-6_r8,  0.08e-6_r8, 1.0e-6_r8, 10.0e-6_r8 /)
  character(len=6),parameter :: seasalt_names(nslt+nnum) = &
       (/ 'ncl_a1', 'ncl_a2', 'ncl_a4', 'ncl_a6', 'num_a1', 'num_a2', 'num_a4', 'num_a6' /)
#elif( defined MODAL_AERO_3MODE || defined MODAL_AERO_4MODE )
  real(r8), parameter :: seasalt_emis_scale = 1.35_r8
  real(r8), parameter :: seasalt_sz_range_lo(nslt) = (/ 0.08e-6_r8,  0.02e-6_r8,  1.0e-6_r8 /)
  real(r8), parameter :: seasalt_sz_range_hi(nslt) = (/ 1.0e-6_r8,   0.08e-6_r8, 10.0e-6_r8 /)
  character(len=6),parameter :: seasalt_names(nslt+nnum) = &
       (/ 'ncl_a1', 'ncl_a2', 'ncl_a3', 'num_a1', 'num_a2', 'num_a3'/)
#endif

  integer :: seasalt_indices(nslt+nnum)

  logical :: seasalt_active = .false.

contains
  
  !=============================================================================
  !=============================================================================
  subroutine seasalt_init
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog
    use sslt_sections, only: sslt_sections_init
    use constituents,  only: cnst_get_ind

    interface
       function seasalt_init_codon(tag) result(tag_out) bind(c, name='seasalt_init_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function seasalt_init_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.


    integer :: m

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('SEASALT_INIT_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = seasalt_init_codon(int(135, c_int64_t))
       if (rt_codon_tag_out /= int(135, c_int64_t)) then
          write(iulog,*) 'seasalt_init_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'seasalt_init implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

    do m = 1, seasalt_nbin
       call cnst_get_ind(seasalt_names(m), seasalt_indices(m),abort=.false.)
    enddo
    do m = 1, seasalt_nnum
       call cnst_get_ind(seasalt_names(seasalt_nbin+m), seasalt_indices(seasalt_nbin+m),abort=.false.)
    enddo

    seasalt_active = any(seasalt_indices(:) > 0)

    if (.not.seasalt_active) return

    call sslt_sections_init()

  end subroutine seasalt_init

  !=============================================================================
  !=============================================================================
  subroutine seasalt_emis( u10cubed,  srf_temp, ocnfrc, ncol, cflx )

    use sslt_sections, only: nsections, fluxes, Dg, rdry
    use mo_constants,  only: dns_aer_sst=>seasalt_density, pi

    ! dummy arguments
    real(r8), intent(in) :: u10cubed(:)
    real(r8), intent(in) :: srf_temp(:)
    real(r8), intent(in) :: ocnfrc(:)
    integer,  intent(in) :: ncol
    real(r8), intent(inout) :: cflx(:,:)

    ! local vars
    integer  :: mn, mm, ibin, isec, i
    real(r8) :: fi(ncol,nsections)

    fi(:ncol,:nsections) = fluxes( srf_temp, u10cubed, ncol )

    do ibin = 1,nslt
       mm = seasalt_indices(ibin)
       mn = seasalt_indices(nslt+ibin)
       
       if (mn>0) then
          do i=1, nsections
             if (Dg(i).ge.seasalt_sz_range_lo(ibin) .and. Dg(i).lt.seasalt_sz_range_hi(ibin)) then
                cflx(:ncol,mn)=cflx(:ncol,mn)+fi(:ncol,i)*ocnfrc(:ncol)*seasalt_emis_scale  !++ ag: scale sea-salt
             endif
          enddo
       endif

       cflx(:ncol,mm)=0.0_r8
       do i=1, nsections
          if (Dg(i).ge.seasalt_sz_range_lo(ibin) .and. Dg(i).lt.seasalt_sz_range_hi(ibin)) then
             cflx(:ncol,mm)=cflx(:ncol,mm)+fi(:ncol,i)*ocnfrc(:ncol)*seasalt_emis_scale  &   !++ ag: scale sea-salt
                  *4._r8/3._r8*pi*rdry(i)**3*dns_aer_sst  ! should use dry size, convert from number to mass flux (kg/m2/s)
          endif
       enddo

    enddo

  end subroutine seasalt_emis

end module seasalt_model
