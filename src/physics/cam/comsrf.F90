
!-----------------------------------------------------------------------
!
! !MODULE: comsrf
!
! !DESCRIPTION:	Module to handle surface fluxes for the subcomponents of cam/csm
!                    Currently this is a hodge-podge of 2D arrays without a lot
!                    of thought or design. We are under the process of removing
!                    this completely and moving the relevent arrays to the modules
!                    that actually use the data.
!
!			See: http://swiki.ucar.edu/start/66
!
!-----------------------------------------------------------------------
module comsrf
!
! USES:
!
  use shr_kind_mod,   only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use ppgrid,         only: pcols, begchunk, endchunk
  use infnan,         only: nan, assignment(=)
  use cam_abortutils, only: endrun
  use cam_logfile,    only: iulog
  use spmd_utils,     only: masterproc

  implicit none

!----------------------------------------------------------------------- 
! PRIVATE: Make default data and interfaces private
!----------------------------------------------------------------------- 
  private     ! By default all data is private to this module
!
! ! PUBLIC MEMBER FUNCTIONS:
!
  public initialize_comsrf          ! Set the surface temperature and sea-ice fraction
!
! Public data
!
  public landm, sgh, sgh30, fv, ram1, soilw, fsns, fsds
  public fsnt, flns, flnt, srfrpdel, psm1, prcsnw
  public trefmxav, trefmnav

  real(r8), allocatable, target :: landm(:,:)     ! land/ocean/sea ice flag
  real(r8), allocatable, target :: sgh(:,:)       ! land/ocean/sea ice flag
  real(r8), allocatable, target :: sgh30(:,:)     ! land/ocean/sea ice flag
  real(r8), allocatable:: fv(:,:)        ! needed for dry dep velocities (over land)
  real(r8), allocatable:: ram1(:,:)      ! needed for dry dep velocities (over land)
  real(r8), allocatable:: soilw(:,:)     ! needed for dust emission (over land)
  real(r8), allocatable, target :: fsns(:,:)      ! surface absorbed solar flux
  real(r8), allocatable, target :: fsds(:,:)      ! downward solar flux
  real(r8), allocatable, target :: fsnt(:,:)      ! Net column abs solar flux at model top
  real(r8), allocatable, target :: flns(:,:)      ! Srf longwave cooling (up-down) flux
  real(r8), allocatable, target :: flnt(:,:)      ! Net outgoing lw flux at model top
  real(r8), allocatable, target :: srfrpdel(:,:)  ! 1./(pint(k+1)-pint(k))
  real(r8), allocatable, target :: psm1(:,:)      ! surface pressure
  real(r8), allocatable, target :: prcsnw(:,:)    ! cam tot snow precip
  real(r8), allocatable, target :: trefmxav(:,:)  ! diagnostic: tref max over the day
  real(r8), allocatable, target :: trefmnav(:,:)  ! diagnostic: tref min over the day

! Private module data
  logical :: use_native_comsrf_init_impl = .false.
  logical :: comsrf_init_impl_selected = .false.
  logical :: comsrf_init_proof_written = .false.

  interface
     subroutine initialize_comsrf_codon(total_len_c, nan_value_c, landm_p, sgh_p, sgh30_p, &
          fsns_p, fsds_p, fsnt_p, flns_p, flnt_p, srfrpdel_p, psm1_p, prcsnw_p, &
          trefmxav_p, trefmnav_p) bind(c, name="initialize_comsrf_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: total_len_c
       real(c_double), value :: nan_value_c
       type(c_ptr), value :: landm_p, sgh_p, sgh30_p, fsns_p, fsds_p, fsnt_p
       type(c_ptr), value :: flns_p, flnt_p, srfrpdel_p, psm1_p, prcsnw_p
       type(c_ptr), value :: trefmxav_p, trefmnav_p
     end subroutine initialize_comsrf_codon
  end interface

!===============================================================================
CONTAINS
!===============================================================================
  subroutine comsrf_init_select_impl()
    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (comsrf_init_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('COMSRF_INIT_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_comsrf_init_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_comsrf_init_impl = .false.
    end if

    comsrf_init_impl_selected = .true.

    if (masterproc) then
       if (use_native_comsrf_init_impl) then
          write(iulog,*) 'comsrf_init implementation = native'
       else
          write(iulog,*) 'comsrf_init implementation = codon'
       end if
    end if
  end subroutine comsrf_init_select_impl

!===============================================================================
  subroutine comsrf_init_proof_once()
    if (comsrf_init_proof_written) return
    comsrf_init_proof_written = .true.
    if (masterproc) then
       write(iulog,'(A)') 'comsrf_init entered (surface field initialization = codon)'
    end if
  end subroutine comsrf_init_proof_once

!===============================================================================

!======================================================================
! PUBLIC ROUTINES: Following routines are publically accessable
!======================================================================
!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: initialize_comsrf
!
! !DESCRIPTION:
!
! Initialize the procedure for specifying sea surface temperatures
! Do initial read of time-varying ice boundary dataset, reading two
! consecutive months on either side of the current model date.
!
! Method: 
! 
! Author: 
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine initialize_comsrf
  use cam_control_mod,  only: ideal_phys, adiabatic
  use iso_c_binding,    only: c_double, c_int64_t, c_loc
!-----------------------------------------------------------------------
!
! Purpose:
! Initialize surface data
!
! Method:
!
! Author: Mariana Vertenstein
!
!-----------------------------------------------------------------------
    integer k,c      ! level, constituent indices
    real(r8) :: real_nan

    if(.not. (adiabatic .or. ideal_phys)) then
       allocate (landm   (pcols,begchunk:endchunk))
       allocate (sgh     (pcols,begchunk:endchunk))
       allocate (sgh30   (pcols,begchunk:endchunk))

       allocate (fv      (pcols,begchunk:endchunk))
       allocate (ram1    (pcols,begchunk:endchunk))
       allocate (soilw   (pcols,begchunk:endchunk))
       allocate (fsns    (pcols,begchunk:endchunk))         
       allocate (fsds    (pcols,begchunk:endchunk))         
       allocate (fsnt    (pcols,begchunk:endchunk))         
       allocate (flns    (pcols,begchunk:endchunk))         
       allocate (flnt    (pcols,begchunk:endchunk))         
       allocate (srfrpdel(pcols,begchunk:endchunk))
       allocate (psm1    (pcols,begchunk:endchunk))
       allocate (prcsnw  (pcols,begchunk:endchunk))
       allocate (trefmxav(pcols,begchunk:endchunk))
       allocate (trefmnav(pcols,begchunk:endchunk))
       !
       ! Initialize to NaN or Inf
       ! elements of the array outside valid surface points must be set to
       ! zero if these fields are to be written to netcdf history files.
       !
       call comsrf_init_select_impl()
       if (use_native_comsrf_init_impl) then
          landm    (:,:) = nan
          sgh      (:,:) = nan
          sgh30    (:,:) = nan
          fsns     (:,:) = nan
          fsds     (:,:) = nan
          fsnt     (:,:) = nan
          flns     (:,:) = nan
          flnt     (:,:) = nan
          srfrpdel (:,:) = nan
          psm1     (:,:) = nan
          prcsnw   (:,:) = nan
          trefmxav (:,:) = -1.0e36_r8
          trefmnav (:,:) =  1.0e36_r8
       else
          real_nan = nan
          call comsrf_init_proof_once()
          call initialize_comsrf_codon(int(pcols*(endchunk-begchunk+1), c_int64_t), &
               real(real_nan, c_double), c_loc(landm(1,begchunk)), c_loc(sgh(1,begchunk)), &
               c_loc(sgh30(1,begchunk)), c_loc(fsns(1,begchunk)), c_loc(fsds(1,begchunk)), &
               c_loc(fsnt(1,begchunk)), c_loc(flns(1,begchunk)), c_loc(flnt(1,begchunk)), &
               c_loc(srfrpdel(1,begchunk)), c_loc(psm1(1,begchunk)), c_loc(prcsnw(1,begchunk)), &
               c_loc(trefmxav(1,begchunk)), c_loc(trefmnav(1,begchunk)))
       end if
    end if
  end subroutine initialize_comsrf

end module comsrf
