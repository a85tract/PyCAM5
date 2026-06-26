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
  logical :: seasalt_emis_codon_logged = .false.

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

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use cam_logfile, only: iulog
    use spmd_utils, only: masterproc
    use sslt_sections, only: nsections, fluxes, Dg, rdry, consta, constb
    use mo_constants,  only: dns_aer_sst=>seasalt_density, pi

    ! dummy arguments
    real(r8), target, intent(in) :: u10cubed(:)
    real(r8), target, intent(in) :: srf_temp(:)
    real(r8), target, intent(in) :: ocnfrc(:)
    integer,  intent(in) :: ncol
    real(r8), target, intent(inout) :: cflx(:,:)

    ! local vars
    integer  :: mn, mm, ibin, isec, i
    integer :: rt_codon_n, rt_codon_status, rt_codon_i, rt_codon_code
    integer(c_int64_t), target :: seasalt_indices_c(nslt+nnum)
    character(len=32) :: rt_codon_impl_name
    logical :: rt_codon_use_native
    real(c_double), target :: seasalt_sz_range_lo_c(nslt)
    real(c_double), target :: seasalt_sz_range_hi_c(nslt)
    real(c_double), target :: dg_c(nsections)
    real(c_double), target :: rdry_c(nsections)
    real(r8), target :: fi(pcols,nsections)
    real(r8), target :: sflx(pcols)
    real(r8), target :: whitecap(pcols)

    interface
       subroutine seasalt_emis_codon(ncol_c, pcols_c, seasalt_nbin_c, nsections_c, &
            seasalt_emis_scale_c, pi_c, seasalt_density_c, u10cubed_p, srf_temp_p, ocnfrc_p, &
            cflx_p, sflx_p, seasalt_indices_p, seasalt_sz_range_lo_p, seasalt_sz_range_hi_p, &
            dg_p, rdry_p, fi_p, whitecap_p, consta_p, constb_p) bind(c, name="seasalt_emis_codon")
         import :: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, seasalt_nbin_c, nsections_c
         real(c_double), value :: seasalt_emis_scale_c, pi_c, seasalt_density_c
         type(c_ptr), value :: u10cubed_p, srf_temp_p, ocnfrc_p, cflx_p, sflx_p
         type(c_ptr), value :: seasalt_indices_p, seasalt_sz_range_lo_p, seasalt_sz_range_hi_p
         type(c_ptr), value :: dg_p, rdry_p, fi_p, whitecap_p, consta_p, constb_p
       end subroutine seasalt_emis_codon
    end interface

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('SEASALT_EMIS_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    rt_codon_use_native = .false.
    if (rt_codon_status == 0 .and. rt_codon_n > 0) then
       do rt_codon_i = 1, rt_codon_n
          rt_codon_code = iachar(rt_codon_impl_name(rt_codon_i:rt_codon_i))
          if (rt_codon_code >= iachar('A') .and. rt_codon_code <= iachar('Z')) then
             rt_codon_impl_name(rt_codon_i:rt_codon_i) = achar(rt_codon_code + iachar('a') - iachar('A'))
          end if
       end do
       rt_codon_use_native = trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native'
    end if

    if (.not. rt_codon_use_native) then
       seasalt_indices_c(:) = int(seasalt_indices(:), c_int64_t)
       seasalt_sz_range_lo_c(:) = real(seasalt_sz_range_lo(:), c_double)
       seasalt_sz_range_hi_c(:) = real(seasalt_sz_range_hi(:), c_double)
       do isec = 1, nsections
          dg_c(isec) = real(Dg(isec), c_double)
          rdry_c(isec) = real(rdry(isec), c_double)
       end do
       call seasalt_emis_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(nslt, c_int64_t), &
            int(nsections, c_int64_t), real(seasalt_emis_scale, c_double), real(pi, c_double), &
            real(dns_aer_sst, c_double), c_loc(u10cubed(1)), c_loc(srf_temp(1)), c_loc(ocnfrc(1)), &
            c_loc(cflx(1,1)), c_loc(sflx(1)), c_loc(seasalt_indices_c(1)), c_loc(seasalt_sz_range_lo_c(1)), &
            c_loc(seasalt_sz_range_hi_c(1)), c_loc(dg_c(1)), c_loc(rdry_c(1)), c_loc(fi(1,1)), &
            c_loc(whitecap(1)), c_loc(consta(1,1)), c_loc(constb(1,1)))
       if (masterproc .and. .not. seasalt_emis_codon_logged) then
          seasalt_emis_codon_logged = .true.
          write(iulog,'(A)') 'seasalt_emis implementation = codon'
          call flush(iulog)
       end if
       return
    end if

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
