module camsrfexch
!-----------------------------------------------------------------------
!
! Module to handle data that is exchanged between the CAM atmosphere
! model and the surface models (land, sea-ice, and ocean).
!
!-----------------------------------------------------------------------
!
! USES:
!
  use shr_kind_mod,  only: r8 => shr_kind_r8, r4 => shr_kind_r4
  use constituents,  only: pcnst
  use ppgrid,        only: pcols, begchunk, endchunk
  use phys_grid,     only: get_ncols_p, phys_grid_initialized
  use infnan,        only: posinf, assignment(=)
  use cam_abortutils,only: endrun
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc

  implicit none

  logical :: cam_export_use_native_impl = .false.
  logical :: cam_export_impl_selected = .false.
  logical :: cam_export_entered_logged = .false.
  logical :: alloc_init_use_native_impl = .false.
  logical :: alloc_init_impl_selected = .false.
  logical :: hub2atm_alloc_init_entered_logged = .false.
  logical :: atm2hub_alloc_init_entered_logged = .false.

!----------------------------------------------------------------------- 
! PRIVATE: Make default data and interfaces private
!----------------------------------------------------------------------- 
  private     ! By default all data is private to this module
!
! Public interfaces
!
  public atm2hub_alloc              ! Atmosphere to surface data allocation method
  public hub2atm_alloc              ! Merged hub surface to atmosphere data allocation method
  public atm2hub_deallocate
  public hub2atm_deallocate
  public cam_export
!
! Public data types
!
  public cam_out_t                  ! Data from atmosphere
  public cam_in_t                   ! Merged surface data

!---------------------------------------------------------------------------
! This is the data that is sent from the atmosphere to the surface models
!---------------------------------------------------------------------------

  type cam_out_t 
     integer  :: lchnk               ! chunk index
     integer  :: ncol                ! number of columns in chunk
     real(r8) :: tbot(pcols)         ! bot level temperature
     real(r8) :: zbot(pcols)         ! bot level height above surface
     real(r8) :: ubot(pcols)         ! bot level u wind
     real(r8) :: vbot(pcols)         ! bot level v wind
     real(r8) :: qbot(pcols,pcnst)   ! bot level specific humidity
     real(r8) :: pbot(pcols)         ! bot level pressure
     real(r8) :: rho(pcols)          ! bot level density	
     real(r8) :: netsw(pcols)        !	
     real(r8) :: flwds(pcols)        ! 
     real(r8) :: precsc(pcols)       !
     real(r8) :: precsl(pcols)       !
     real(r8) :: precc(pcols)        ! 
     real(r8) :: precl(pcols)        ! 
     real(r8) :: soll(pcols)         ! 
     real(r8) :: sols(pcols)         ! 
     real(r8) :: solld(pcols)        !
     real(r8) :: solsd(pcols)        !
     real(r8) :: thbot(pcols)        ! 
     real(r8) :: co2prog(pcols)      ! prognostic co2
     real(r8) :: co2diag(pcols)      ! diagnostic co2
     real(r8) :: psl(pcols)
     real(r8) :: bcphiwet(pcols)     ! wet deposition of hydrophilic black carbon
     real(r8) :: bcphidry(pcols)     ! dry deposition of hydrophilic black carbon
     real(r8) :: bcphodry(pcols)     ! dry deposition of hydrophobic black carbon
     real(r8) :: ocphiwet(pcols)     ! wet deposition of hydrophilic organic carbon
     real(r8) :: ocphidry(pcols)     ! dry deposition of hydrophilic organic carbon
     real(r8) :: ocphodry(pcols)     ! dry deposition of hydrophobic organic carbon
     real(r8) :: dstwet1(pcols)      ! wet deposition of dust (bin1)
     real(r8) :: dstdry1(pcols)      ! dry deposition of dust (bin1)
     real(r8) :: dstwet2(pcols)      ! wet deposition of dust (bin2)
     real(r8) :: dstdry2(pcols)      ! dry deposition of dust (bin2)
     real(r8) :: dstwet3(pcols)      ! wet deposition of dust (bin3)
     real(r8) :: dstdry3(pcols)      ! dry deposition of dust (bin3)
     real(r8) :: dstwet4(pcols)      ! wet deposition of dust (bin4)
     real(r8) :: dstdry4(pcols)      ! dry deposition of dust (bin4)
     !water tracers/isotopes:
     real(r8) :: precrl_16O(pcols)   !Large-scale rain
     real(r8) :: precrl_HDO(pcols)
     real(r8) :: precrl_18O(pcols)
     real(r8) :: precsl_16O(pcols)   !Large-scale snow
     real(r8) :: precsl_HDO(pcols)
     real(r8) :: precsl_18O(pcols)
     real(r8) :: precrc_16O(pcols)   !Convective rain
     real(r8) :: precrc_HDO(pcols)
     real(r8) :: precrc_18O(pcols)   
     real(r8) :: precsc_16O(pcols)   !Convective snow
     real(r8) :: precsc_HDO(pcols)
     real(r8) :: precsc_18O(pcols)
  end type cam_out_t 

!---------------------------------------------------------------------------
! This is the merged state of sea-ice, land and ocean surface parameterizations
!---------------------------------------------------------------------------

  type cam_in_t    
     integer  :: lchnk                   ! chunk index
     integer  :: ncol                    ! number of active columns
     real(r8) :: asdir(pcols)            ! albedo: shortwave, direct
     real(r8) :: asdif(pcols)            ! albedo: shortwave, diffuse
     real(r8) :: aldir(pcols)            ! albedo: longwave, direct
     real(r8) :: aldif(pcols)            ! albedo: longwave, diffuse
     real(r8) :: lwup(pcols)             ! longwave up radiative flux
     real(r8) :: lhf(pcols)              ! latent heat flux
     real(r8) :: shf(pcols)              ! sensible heat flux
     real(r8) :: wsx(pcols)              ! surface u-stress (N)
     real(r8) :: wsy(pcols)              ! surface v-stress (N)
     real(r8) :: tref(pcols)             ! ref height surface air temp
     real(r8) :: qref(pcols)             ! ref height specific humidity 
     real(r8) :: u10(pcols)              ! 10m wind speed
     real(r8) :: ts(pcols)               ! merged surface temp 
     real(r8) :: sst(pcols)              ! sea surface temp
     real(r8) :: snowhland(pcols)        ! snow depth (liquid water equivalent) over land 
     real(r8) :: snowhice(pcols)         ! snow depth over ice
     real(r8) :: fco2_lnd(pcols)         ! co2 flux from lnd
     real(r8) :: fco2_ocn(pcols)         ! co2 flux from ocn
     real(r8) :: fdms(pcols)             ! dms flux
     real(r8) :: landfrac(pcols)         ! land area fraction
     real(r8) :: icefrac(pcols)          ! sea-ice areal fraction
     real(r8) :: ocnfrac(pcols)          ! ocean areal fraction
     real(r8), pointer, dimension(:) :: ram1  !aerodynamical resistance (s/m) (pcols)
     real(r8), pointer, dimension(:) :: fv    !friction velocity (m/s) (pcols)
     real(r8), pointer, dimension(:) :: soilw !volumetric soil water (m3/m3)
     real(r8) :: cflx(pcols,pcnst)       ! constituent flux (emissions)
     real(r8) :: ustar(pcols)            ! atm/ocn saved version of ustar
     real(r8) :: re(pcols)               ! atm/ocn saved version of re
     real(r8) :: ssq(pcols)              ! atm/ocn saved version of ssq
     real(r8), pointer, dimension(:,:) :: depvel ! deposition velocities
     real(r8), pointer, dimension(:,:) :: dstflx ! dust fluxes
     real(r8), pointer, dimension(:,:) :: meganflx ! MEGAN fluxes
     !simple land model:
     real(r8) :: buckH(pcols)        !bulk water bucket mass
     real(r8) :: buck16(pcols)       !isotopic water bucket masses
     real(r8) :: buckD(pcols)
     real(r8) :: buck18(pcols)
  end type cam_in_t    

!===============================================================================
CONTAINS
!===============================================================================

subroutine camsrfexch_alloc_init_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line
   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('CAMSRFEXCH_ALLOC_INIT_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine camsrfexch_alloc_init_append_proof

subroutine camsrfexch_alloc_init_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (alloc_init_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CAMSRFEXCH_ALLOC_INIT_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      alloc_init_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      alloc_init_use_native_impl = .false.
   end if

   alloc_init_impl_selected = .true.

   if (masterproc) then
      if (alloc_init_use_native_impl) then
         write(iulog,*) 'camsrfexch_alloc_init implementation = native'
         call camsrfexch_alloc_init_append_proof('camsrfexch_alloc_init selector entered implementation = native')
      else
         write(iulog,*) 'camsrfexch_alloc_init implementation = codon'
         call camsrfexch_alloc_init_append_proof('camsrfexch_alloc_init selector entered implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine camsrfexch_alloc_init_select_impl

subroutine hub2atm_alloc_init_log_entered()

   if (hub2atm_alloc_init_entered_logged) return
   hub2atm_alloc_init_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'hub2atm_alloc_init entered (unified surface input init stage dispatch = codon)'
      call camsrfexch_alloc_init_append_proof('hub2atm_alloc_init entered (unified surface input init stage dispatch = codon)')
      call flush(iulog)
   end if

end subroutine hub2atm_alloc_init_log_entered

subroutine atm2hub_alloc_init_log_entered()

   if (atm2hub_alloc_init_entered_logged) return
   atm2hub_alloc_init_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'atm2hub_alloc_init entered (unified surface output init stage dispatch = codon)'
      call camsrfexch_alloc_init_append_proof('atm2hub_alloc_init entered (unified surface output init stage dispatch = codon)')
      call flush(iulog)
   end if

end subroutine atm2hub_alloc_init_log_entered

!===============================================================================

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: hub2atm_alloc
!
! !DESCRIPTION:
!
!   Allocate space for the surface to atmosphere data type. And initialize
!   the values.
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine hub2atm_alloc( cam_in )
    use seq_drydep_mod,  only: lnd_drydep, n_drydep
    use cam_cpl_indices, only: index_x2a_Sl_ram1, index_x2a_Sl_fv, index_x2a_Sl_soilw, index_x2a_Fall_flxdst1
    use cam_cpl_indices, only: index_x2a_Fall_flxvoc
    use shr_megan_mod,   only: shr_megan_mechcomps_n
    use iso_c_binding,   only: c_double, c_int64_t, c_loc, c_ptr

   use water_tracer_vars, only: wtrc_srf_bucket_mode
!
!!ARGUMENTS:
!
   type(cam_in_t), pointer ::  cam_in(:)     ! Merged surface state
!
!!LOCAL VARIABLES:
!
    integer :: c        ! chunk index
    integer :: ierror   ! Error code
    real(r8) :: posinf_r8
    real(r8), target :: dummy_1d(1)
    real(r8), target :: dummy_2d(1,1)
    type(c_ptr) :: ram1_p, fv_p, soilw_p, dstflx_p, meganflx_p, depvel_p
    interface
       subroutine hub2atm_alloc_init_stage_dispatch_codon(pcols_c, pcnst_c, n_drydep_c, n_megan_c, posinf_c, &
            init_bucket_c, has_ram1_c, has_fv_c, has_soilw_c, has_dstflx_c, has_meganflx_c, has_depvel_c, &
            asdir_p, asdif_p, aldir_p, aldif_p, lwup_p, lhf_p, shf_p, wsx_p, wsy_p, tref_p, qref_p, &
            u10_p, ts_p, sst_p, snowhland_p, snowhice_p, fco2_lnd_p, fco2_ocn_p, fdms_p, &
            landfrac_p, icefrac_p, ocnfrac_p, ram1_p, fv_p, soilw_p, cflx_p, ustar_p, re_p, ssq_p, &
            dstflx_p, meganflx_p, depvel_p, buckH_p, buck16_p, buckD_p, buck18_p) &
            bind(c, name="hub2atm_alloc_init_stage_dispatch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: pcols_c, pcnst_c, n_drydep_c, n_megan_c
         real(c_double), value :: posinf_c
         integer(c_int64_t), value :: init_bucket_c, has_ram1_c, has_fv_c, has_soilw_c
         integer(c_int64_t), value :: has_dstflx_c, has_meganflx_c, has_depvel_c
         type(c_ptr), value :: asdir_p, asdif_p, aldir_p, aldif_p, lwup_p, lhf_p, shf_p, wsx_p, wsy_p
         type(c_ptr), value :: tref_p, qref_p, u10_p, ts_p, sst_p, snowhland_p, snowhice_p
         type(c_ptr), value :: fco2_lnd_p, fco2_ocn_p, fdms_p, landfrac_p, icefrac_p, ocnfrac_p
         type(c_ptr), value :: ram1_p, fv_p, soilw_p, cflx_p, ustar_p, re_p, ssq_p
         type(c_ptr), value :: dstflx_p, meganflx_p, depvel_p, buckH_p, buck16_p, buckD_p, buck18_p
       end subroutine hub2atm_alloc_init_stage_dispatch_codon
    end interface
!----------------------------------------------------------------------- 
! 
! EOP
!
    if ( .not. phys_grid_initialized() ) call endrun( "HUB2ATM_ALLOC error: phys_grid not called yet" )
    allocate (cam_in(begchunk:endchunk), stat=ierror)
    if ( ierror /= 0 )then
      write(iulog,*) 'Allocation error: ', ierror
      call endrun('HUB2ATM_ALLOC error: allocation error')
    end if

    do c = begchunk,endchunk
       nullify(cam_in(c)%ram1)
       nullify(cam_in(c)%fv)
       nullify(cam_in(c)%soilw)
       nullify(cam_in(c)%depvel)
       nullify(cam_in(c)%dstflx)
       nullify(cam_in(c)%meganflx)
    enddo  
    do c = begchunk,endchunk 
       if (index_x2a_Sl_ram1>0) then
          allocate (cam_in(c)%ram1(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error ram1')
       endif
       if (index_x2a_Sl_fv>0) then
          allocate (cam_in(c)%fv(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error fv')
       endif
       if (index_x2a_Sl_soilw /= 0) then
          allocate (cam_in(c)%soilw(pcols), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error soilw')
       end if
       if (index_x2a_Fall_flxdst1>0) then
          ! Assume 4 bins from surface model ....
          allocate (cam_in(c)%dstflx(pcols,4), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error dstflx')
       endif
       if ( index_x2a_Fall_flxvoc>0 .and. shr_megan_mechcomps_n>0 ) then
          allocate (cam_in(c)%meganflx(pcols,shr_megan_mechcomps_n), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error meganflx')
       endif
    end do

    if (lnd_drydep .and. n_drydep>0) then
       do c = begchunk,endchunk 
          allocate (cam_in(c)%depvel(pcols,n_drydep), stat=ierror)
          if ( ierror /= 0 ) call endrun('HUB2ATM_ALLOC error: allocation error depvel')
       end do
    endif

    call camsrfexch_alloc_init_select_impl()

    posinf_r8 = posinf
    dummy_1d(:) = 0._r8
    dummy_2d(:,:) = 0._r8
    do c = begchunk,endchunk
       cam_in(c)%lchnk = c
       cam_in(c)%ncol  = get_ncols_p(c)
       if (alloc_init_use_native_impl) then
          cam_in(c)%asdir    (:) = 0._r8
          cam_in(c)%asdif    (:) = 0._r8
          cam_in(c)%aldir    (:) = 0._r8
          cam_in(c)%aldif    (:) = 0._r8
          cam_in(c)%lwup     (:) = 0._r8
          cam_in(c)%lhf      (:) = 0._r8
          cam_in(c)%shf      (:) = 0._r8
          cam_in(c)%wsx      (:) = 0._r8
          cam_in(c)%wsy      (:) = 0._r8
          cam_in(c)%tref     (:) = 0._r8
          cam_in(c)%qref     (:) = 0._r8
          cam_in(c)%u10      (:) = 0._r8
          cam_in(c)%ts       (:) = 0._r8
          cam_in(c)%sst      (:) = 0._r8
          cam_in(c)%snowhland(:) = 0._r8
          cam_in(c)%snowhice (:) = 0._r8
          cam_in(c)%fco2_lnd (:) = 0._r8
          cam_in(c)%fco2_ocn (:) = 0._r8
          cam_in(c)%fdms     (:) = 0._r8
          cam_in(c)%landfrac (:) = posinf
          cam_in(c)%icefrac  (:) = posinf
          cam_in(c)%ocnfrac  (:) = posinf

          if (associated(cam_in(c)%ram1)) &
               cam_in(c)%ram1  (:) = 0.1_r8
          if (associated(cam_in(c)%fv)) &
               cam_in(c)%fv    (:) = 0.1_r8
          if (associated(cam_in(c)%soilw)) &
               cam_in(c)%soilw (:) = 0.0_r8
          if (associated(cam_in(c)%dstflx)) &
               cam_in(c)%dstflx(:,:) = 0.0_r8
          if (associated(cam_in(c)%meganflx)) &
               cam_in(c)%meganflx(:,:) = 0.0_r8

          cam_in(c)%cflx   (:,:) = 0._r8
          cam_in(c)%ustar    (:) = 0._r8
          cam_in(c)%re       (:) = 0._r8
          cam_in(c)%ssq      (:) = 0._r8
          if (lnd_drydep .and. n_drydep>0) then
             cam_in(c)%depvel (:,:) = 0._r8
          endif
          !For simple land model:
          if ( wtrc_srf_bucket_mode ) then
             cam_in(c)%buckH(:)    = 0.021_r8 !initalize as half-full bucket
             cam_in(c)%buck16(:)   = 0.021_r8
             cam_in(c)%buckD(:)    = 0.021_r8
             cam_in(c)%buck18(:)   = 0.021_r8
          end if
       else
          if (associated(cam_in(c)%ram1)) then
             ram1_p = c_loc(cam_in(c)%ram1(1))
          else
             ram1_p = c_loc(dummy_1d(1))
          end if
          if (associated(cam_in(c)%fv)) then
             fv_p = c_loc(cam_in(c)%fv(1))
          else
             fv_p = c_loc(dummy_1d(1))
          end if
          if (associated(cam_in(c)%soilw)) then
             soilw_p = c_loc(cam_in(c)%soilw(1))
          else
             soilw_p = c_loc(dummy_1d(1))
          end if
          if (associated(cam_in(c)%dstflx)) then
             dstflx_p = c_loc(cam_in(c)%dstflx(1,1))
          else
             dstflx_p = c_loc(dummy_2d(1,1))
          end if
          if (associated(cam_in(c)%meganflx)) then
             meganflx_p = c_loc(cam_in(c)%meganflx(1,1))
          else
             meganflx_p = c_loc(dummy_2d(1,1))
          end if
          if (associated(cam_in(c)%depvel)) then
             depvel_p = c_loc(cam_in(c)%depvel(1,1))
          else
             depvel_p = c_loc(dummy_2d(1,1))
          end if
          call hub2atm_alloc_init_log_entered()
          call hub2atm_alloc_init_stage_dispatch_codon( &
               int(pcols, c_int64_t), int(pcnst, c_int64_t), int(n_drydep, c_int64_t), &
               int(shr_megan_mechcomps_n, c_int64_t), real(posinf_r8, c_double), &
               merge(1_c_int64_t, 0_c_int64_t, wtrc_srf_bucket_mode), &
               merge(1_c_int64_t, 0_c_int64_t, associated(cam_in(c)%ram1)), &
               merge(1_c_int64_t, 0_c_int64_t, associated(cam_in(c)%fv)), &
               merge(1_c_int64_t, 0_c_int64_t, associated(cam_in(c)%soilw)), &
               merge(1_c_int64_t, 0_c_int64_t, associated(cam_in(c)%dstflx)), &
               merge(1_c_int64_t, 0_c_int64_t, associated(cam_in(c)%meganflx)), &
               merge(1_c_int64_t, 0_c_int64_t, associated(cam_in(c)%depvel)), &
               c_loc(cam_in(c)%asdir(1)), c_loc(cam_in(c)%asdif(1)), c_loc(cam_in(c)%aldir(1)), &
               c_loc(cam_in(c)%aldif(1)), c_loc(cam_in(c)%lwup(1)), c_loc(cam_in(c)%lhf(1)), &
               c_loc(cam_in(c)%shf(1)), c_loc(cam_in(c)%wsx(1)), c_loc(cam_in(c)%wsy(1)), &
               c_loc(cam_in(c)%tref(1)), c_loc(cam_in(c)%qref(1)), c_loc(cam_in(c)%u10(1)), &
               c_loc(cam_in(c)%ts(1)), c_loc(cam_in(c)%sst(1)), c_loc(cam_in(c)%snowhland(1)), &
               c_loc(cam_in(c)%snowhice(1)), c_loc(cam_in(c)%fco2_lnd(1)), c_loc(cam_in(c)%fco2_ocn(1)), &
               c_loc(cam_in(c)%fdms(1)), c_loc(cam_in(c)%landfrac(1)), c_loc(cam_in(c)%icefrac(1)), &
               c_loc(cam_in(c)%ocnfrac(1)), ram1_p, fv_p, soilw_p, c_loc(cam_in(c)%cflx(1,1)), &
               c_loc(cam_in(c)%ustar(1)), c_loc(cam_in(c)%re(1)), c_loc(cam_in(c)%ssq(1)), &
               dstflx_p, meganflx_p, depvel_p, c_loc(cam_in(c)%buckH(1)), c_loc(cam_in(c)%buck16(1)), &
               c_loc(cam_in(c)%buckD(1)), c_loc(cam_in(c)%buck18(1)) &
          )
       end if
    end do

  end subroutine hub2atm_alloc

!
!===============================================================================
!

!----------------------------------------------------------------------- 
! 
! BOP
!
! !IROUTINE: atm2hub_alloc
!
! !DESCRIPTION:
!
!   Allocate space for the atmosphere to surface data type. And initialize
!   the values.
! 
!-----------------------------------------------------------------------
!
! !INTERFACE
!
  subroutine atm2hub_alloc( cam_out )
!
!!USES:
   use water_tracer_vars, only: wtrc_srf_bucket_mode
   use iso_c_binding,     only: c_int64_t, c_loc, c_ptr
!
!
!!ARGUMENTS:
!
   type(cam_out_t), pointer :: cam_out(:)    ! Atmosphere to surface input
!
!!LOCAL VARIABLES:
!
    integer :: c            ! chunk index
    integer :: ierror       ! Error code
    interface
       subroutine atm2hub_alloc_init_stage_dispatch_codon(pcols_c, pcnst_c, tbot_p, zbot_p, ubot_p, vbot_p, qbot_p, &
            pbot_p, rho_p, netsw_p, flwds_p, precsc_p, precsl_p, precc_p, precl_p, soll_p, sols_p, &
            solld_p, solsd_p, thbot_p, co2prog_p, co2diag_p, psl_p, bcphidry_p, bcphodry_p, bcphiwet_p, &
            ocphidry_p, ocphodry_p, ocphiwet_p, dstdry1_p, dstwet1_p, dstdry2_p, dstwet2_p, &
            dstdry3_p, dstwet3_p, dstdry4_p, dstwet4_p, precrl_16O_p, precrl_HDO_p, precrl_18O_p, &
            precsl_16O_p, precsl_HDO_p, precsl_18O_p, precrc_16O_p, precrc_HDO_p, precrc_18O_p, &
            precsc_16O_p, precsc_HDO_p, precsc_18O_p) bind(c, name="atm2hub_alloc_init_stage_dispatch_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: pcols_c, pcnst_c
         type(c_ptr), value :: tbot_p, zbot_p, ubot_p, vbot_p, qbot_p, pbot_p, rho_p
         type(c_ptr), value :: netsw_p, flwds_p, precsc_p, precsl_p, precc_p, precl_p
         type(c_ptr), value :: soll_p, sols_p, solld_p, solsd_p, thbot_p, co2prog_p, co2diag_p, psl_p
         type(c_ptr), value :: bcphidry_p, bcphodry_p, bcphiwet_p, ocphidry_p, ocphodry_p, ocphiwet_p
         type(c_ptr), value :: dstdry1_p, dstwet1_p, dstdry2_p, dstwet2_p, dstdry3_p, dstwet3_p, dstdry4_p, dstwet4_p
         type(c_ptr), value :: precrl_16O_p, precrl_HDO_p, precrl_18O_p, precsl_16O_p, precsl_HDO_p, precsl_18O_p
         type(c_ptr), value :: precrc_16O_p, precrc_HDO_p, precrc_18O_p, precsc_16O_p, precsc_HDO_p, precsc_18O_p
       end subroutine atm2hub_alloc_init_stage_dispatch_codon
    end interface
    !----------------------------------------------------------------------- 

    if ( .not. phys_grid_initialized() ) call endrun( "ATM2HUB_ALLOC error: phys_grid not called yet" )
    allocate (cam_out(begchunk:endchunk), stat=ierror)
    if ( ierror /= 0 )then
      write(iulog,*) 'Allocation error: ', ierror
      call endrun('ATM2HUB_ALLOC error: allocation error')
    end if

    call camsrfexch_alloc_init_select_impl()

    do c = begchunk,endchunk
       cam_out(c)%lchnk       = c
       cam_out(c)%ncol        = get_ncols_p(c)
       if (alloc_init_use_native_impl) then
          cam_out(c)%tbot(:)     = 0._r8
          cam_out(c)%zbot(:)     = 0._r8
          cam_out(c)%ubot(:)     = 0._r8
          cam_out(c)%vbot(:)     = 0._r8
          cam_out(c)%qbot(:,:)   = 0._r8
          cam_out(c)%pbot(:)     = 0._r8
          cam_out(c)%rho(:)      = 0._r8
          cam_out(c)%netsw(:)    = 0._r8
          cam_out(c)%flwds(:)    = 0._r8
          cam_out(c)%precsc(:)   = 0._r8
          cam_out(c)%precsl(:)   = 0._r8
          cam_out(c)%precc(:)    = 0._r8
          cam_out(c)%precl(:)    = 0._r8
          cam_out(c)%soll(:)     = 0._r8
          cam_out(c)%sols(:)     = 0._r8
          cam_out(c)%solld(:)    = 0._r8
          cam_out(c)%solsd(:)    = 0._r8
          cam_out(c)%thbot(:)    = 0._r8
          cam_out(c)%co2prog(:)  = 0._r8
          cam_out(c)%co2diag(:)  = 0._r8
          cam_out(c)%psl(:)      = 0._r8
          cam_out(c)%bcphidry(:) = 0._r8
          cam_out(c)%bcphodry(:) = 0._r8
          cam_out(c)%bcphiwet(:) = 0._r8
          cam_out(c)%ocphidry(:) = 0._r8
          cam_out(c)%ocphodry(:) = 0._r8
          cam_out(c)%ocphiwet(:) = 0._r8
          cam_out(c)%dstdry1(:)  = 0._r8
          cam_out(c)%dstwet1(:)  = 0._r8
          cam_out(c)%dstdry2(:)  = 0._r8
          cam_out(c)%dstwet2(:)  = 0._r8
          cam_out(c)%dstdry3(:)  = 0._r8
          cam_out(c)%dstwet3(:)  = 0._r8
          cam_out(c)%dstdry4(:)  = 0._r8
          cam_out(c)%dstwet4(:)  = 0._r8
          !water tracers/isotopes:
          cam_out(c)%precrl_16O  = 0._r8
          cam_out(c)%precrl_HDO  = 0._r8
          cam_out(c)%precrl_18O  = 0._r8
          cam_out(c)%precsl_16O  = 0._r8
          cam_out(c)%precsl_HDO  = 0._r8
          cam_out(c)%precsl_18O  = 0._r8
          cam_out(c)%precrc_16O  = 0._r8
          cam_out(c)%precrc_HDO  = 0._r8
          cam_out(c)%precrc_18O  = 0._r8
          cam_out(c)%precsc_16O  = 0._r8
          cam_out(c)%precsc_HDO  = 0._r8
          cam_out(c)%precsc_18O  = 0._r8
       else
          call atm2hub_alloc_init_log_entered()
          call atm2hub_alloc_init_stage_dispatch_codon( &
               int(pcols, c_int64_t), int(pcnst, c_int64_t), c_loc(cam_out(c)%tbot(1)), &
               c_loc(cam_out(c)%zbot(1)), c_loc(cam_out(c)%ubot(1)), c_loc(cam_out(c)%vbot(1)), &
               c_loc(cam_out(c)%qbot(1,1)), c_loc(cam_out(c)%pbot(1)), c_loc(cam_out(c)%rho(1)), &
               c_loc(cam_out(c)%netsw(1)), c_loc(cam_out(c)%flwds(1)), c_loc(cam_out(c)%precsc(1)), &
               c_loc(cam_out(c)%precsl(1)), c_loc(cam_out(c)%precc(1)), c_loc(cam_out(c)%precl(1)), &
               c_loc(cam_out(c)%soll(1)), c_loc(cam_out(c)%sols(1)), c_loc(cam_out(c)%solld(1)), &
               c_loc(cam_out(c)%solsd(1)), c_loc(cam_out(c)%thbot(1)), c_loc(cam_out(c)%co2prog(1)), &
               c_loc(cam_out(c)%co2diag(1)), c_loc(cam_out(c)%psl(1)), c_loc(cam_out(c)%bcphidry(1)), &
               c_loc(cam_out(c)%bcphodry(1)), c_loc(cam_out(c)%bcphiwet(1)), c_loc(cam_out(c)%ocphidry(1)), &
               c_loc(cam_out(c)%ocphodry(1)), c_loc(cam_out(c)%ocphiwet(1)), c_loc(cam_out(c)%dstdry1(1)), &
               c_loc(cam_out(c)%dstwet1(1)), c_loc(cam_out(c)%dstdry2(1)), c_loc(cam_out(c)%dstwet2(1)), &
               c_loc(cam_out(c)%dstdry3(1)), c_loc(cam_out(c)%dstwet3(1)), c_loc(cam_out(c)%dstdry4(1)), &
               c_loc(cam_out(c)%dstwet4(1)), c_loc(cam_out(c)%precrl_16O(1)), c_loc(cam_out(c)%precrl_HDO(1)), &
               c_loc(cam_out(c)%precrl_18O(1)), c_loc(cam_out(c)%precsl_16O(1)), c_loc(cam_out(c)%precsl_HDO(1)), &
               c_loc(cam_out(c)%precsl_18O(1)), c_loc(cam_out(c)%precrc_16O(1)), c_loc(cam_out(c)%precrc_HDO(1)), &
               c_loc(cam_out(c)%precrc_18O(1)), c_loc(cam_out(c)%precsc_16O(1)), c_loc(cam_out(c)%precsc_HDO(1)), &
               c_loc(cam_out(c)%precsc_18O(1)) &
          )
       end if
    end do

  end subroutine atm2hub_alloc

  subroutine atm2hub_deallocate(cam_out)
    type(cam_out_t), pointer :: cam_out(:)    ! Atmosphere to surface input
    if(associated(cam_out)) then
       deallocate(cam_out)
    end if
    nullify(cam_out)

  end subroutine atm2hub_deallocate
  subroutine hub2atm_deallocate(cam_in)
    type(cam_in_t), pointer :: cam_in(:)    ! Atmosphere to surface input
    integer :: c

    if(associated(cam_in)) then
       do c=begchunk,endchunk
          if(associated(cam_in(c)%ram1)) then
             deallocate(cam_in(c)%ram1)
             nullify(cam_in(c)%ram1)
          end if
          if(associated(cam_in(c)%fv)) then
             deallocate(cam_in(c)%fv)
             nullify(cam_in(c)%fv)
          end if
          if(associated(cam_in(c)%soilw)) then
             deallocate(cam_in(c)%soilw)
             nullify(cam_in(c)%soilw)
          end if
          if(associated(cam_in(c)%dstflx)) then
             deallocate(cam_in(c)%dstflx)
             nullify(cam_in(c)%dstflx)
          end if
          if(associated(cam_in(c)%meganflx)) then
             deallocate(cam_in(c)%meganflx)
             nullify(cam_in(c)%meganflx)
          end if
          if(associated(cam_in(c)%depvel)) then
             deallocate(cam_in(c)%depvel)
             nullify(cam_in(c)%depvel)
          end if
          
       enddo

       deallocate(cam_in)
    end if
    nullify(cam_in)

  end subroutine hub2atm_deallocate


!======================================================================

subroutine cam_export_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line
   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('CAM_EXPORT_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine cam_export_append_proof

subroutine cam_export_log_entered()

   if (cam_export_entered_logged) return
   cam_export_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'cam_export entered (unified state/qbot/precip/surface/co2/water transfer stage dispatch = codon)'
      call cam_export_append_proof('cam_export entered (unified state/qbot/precip/surface/co2/water transfer stage dispatch = codon)')
      call flush(iulog)
   end if

end subroutine cam_export_log_entered

subroutine cam_export_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cam_export_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CAM_EXPORT_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      cam_export_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      cam_export_use_native_impl = .false.
   end if

   cam_export_impl_selected = .true.

   if (masterproc) then
      if (cam_export_use_native_impl) then
         write(iulog,*) 'cam_export implementation = native'
         call cam_export_append_proof('cam_export selector entered implementation = native')
      else
         write(iulog,*) 'cam_export implementation = codon'
         call cam_export_append_proof('cam_export selector entered implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine cam_export_select_impl

!======================================================================

subroutine cam_export(state,cam_out,pbuf)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Transfer atmospheric fields into necessary surface data structures
! 
! Author: L. Bath  CMS Contact: M. Vertenstein
! 
!-----------------------------------------------------------------------
   use physics_types,    only: physics_state
   use ppgrid,           only: pver
   use cam_history,      only: outfld
   use comsrf,           only: psm1, srfrpdel, prcsnw
   use chem_surfvals,    only: chem_surfvals_get
   use co2_cycle,        only: co2_transport, c_i
   use physconst,        only: mwdry, mwco2
   use constituents,     only: pcnst
   use cam_control_mod,  only: rair
   use physics_buffer,   only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
   !water tracers:
   use water_tracer_vars,only: trace_water, wtrc_srfpcp_indices, wtrc_nwset, wtrc_iatype,&
                               iwspec
   use water_types,      only: iwtstrain, iwtstsnow, iwtcvrain, iwtcvsnow
   use water_isotopes,   only: isph2o, isph216o, isphdo, isph218o
   use iso_c_binding,    only: c_double, c_int64_t, c_loc, c_ptr
   implicit none

   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
   !
   type(physics_state),  target, intent(in)    :: state
   type (cam_out_t),     target, intent(inout) :: cam_out
   type(physics_buffer_desc), pointer  :: pbuf(:)

   !
   !---------------------------Local variables-----------------------------
   !
   integer :: i              ! Longitude index
   integer :: m              ! constituent index
   integer :: lchnk          ! Chunk index
   integer :: ncol
   integer :: prec_dp_idx, snow_dp_idx, prec_sh_idx, snow_sh_idx
   integer :: prec_sed_idx,snow_sed_idx,prec_pcw_idx,snow_pcw_idx
   real(r8), target :: dummy_1d(1)
   real(r8) :: co2diag_val
   type(c_ptr) :: precrl_16O_p, precsl_16O_p, precrc_16O_p, precsc_16O_p
   type(c_ptr) :: precrl_HDO_p, precsl_HDO_p, precrc_HDO_p, precsc_HDO_p
   type(c_ptr) :: precrl_18O_p, precsl_18O_p, precrc_18O_p, precsc_18O_p

   real(r8), pointer :: prec_dp(:)                 ! total precipitation   from ZM convection
   real(r8), pointer :: snow_dp(:)                 ! snow from ZM   convection
   real(r8), pointer :: prec_sh(:)                 ! total precipitation   from Hack convection
   real(r8), pointer :: snow_sh(:)                 ! snow from   Hack   convection
   real(r8), pointer :: prec_sed(:)                ! total precipitation   from ZM convection
   real(r8), pointer :: snow_sed(:)                ! snow from ZM   convection
   real(r8), pointer :: prec_pcw(:)                ! total precipitation   from Hack convection
   real(r8), pointer :: snow_pcw(:)                ! snow from Hack   convection
   !water tracers/isotopes:
   logical           :: exist16, existD, exist18   !logicals that determine whether or not species exists in run.
   logical           :: pass16, passD, pass18      !logicals that prevent the passing of water tag infromation to iCLM4.
   real(r8), pointer :: precrl_16O(:)
   real(r8), pointer :: precrl_HDO(:)
   real(r8), pointer :: precrl_18O(:)
   real(r8), pointer :: precsl_16O(:)
   real(r8), pointer :: precsl_HDO(:)
   real(r8), pointer :: precsl_18O(:)
   real(r8), pointer :: precrc_16O(:)
   real(r8), pointer :: precrc_HDO(:)
   real(r8), pointer :: precrc_18O(:)
   real(r8), pointer :: precsc_16O(:)
   real(r8), pointer :: precsc_HDO(:)
   real(r8), pointer :: precsc_18O(:)

   interface
      subroutine cam_export_core_stage_dispatch_codon(ncol_c, pcols_c, pver_c, pcnst_c, rair_c, &
           mwdry_c, mwco2_c, co2diag_val_c, co2_transport_c, co2_idx_c, trace_water_c, &
           exist16_c, existD_c, exist18_c, &
           state_t_p, state_exner_p, state_zm_p, state_u_p, state_v_p, state_pmid_p, state_q_p, &
           state_ps_p, state_rpdel_p, psm1_p, srfrpdel_p, co2diag_p, co2prog_p, prcsnw_p, &
           prec_dp_p, snow_dp_p, prec_sh_p, snow_sh_p, prec_sed_p, snow_sed_p, prec_pcw_p, snow_pcw_p, &
           tbot_p, thbot_p, zbot_p, ubot_p, vbot_p, pbot_p, rho_p, qbot_p, &
           precc_p, precl_p, precsc_p, precsl_p, &
           precrl_16O_in_p, precsl_16O_in_p, precrc_16O_in_p, precsc_16O_in_p, &
           precrl_HDO_in_p, precsl_HDO_in_p, precrc_HDO_in_p, precsc_HDO_in_p, &
           precrl_18O_in_p, precsl_18O_in_p, precrc_18O_in_p, precsc_18O_in_p, &
           precrl_16O_out_p, precsl_16O_out_p, precrc_16O_out_p, precsc_16O_out_p, &
           precrl_HDO_out_p, precsl_HDO_out_p, precrc_HDO_out_p, precsc_HDO_out_p, &
           precrl_18O_out_p, precsl_18O_out_p, precrc_18O_out_p, precsc_18O_out_p) &
           bind(c, name="cam_export_core_stage_dispatch_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c
         integer(c_int64_t), value :: co2_transport_c, co2_idx_c, trace_water_c
         integer(c_int64_t), value :: exist16_c, existD_c, exist18_c
         real(c_double), value :: rair_c, mwdry_c, mwco2_c, co2diag_val_c
         type(c_ptr), value :: state_t_p, state_exner_p, state_zm_p, state_u_p, state_v_p, state_pmid_p, state_q_p
         type(c_ptr), value :: state_ps_p, state_rpdel_p, psm1_p, srfrpdel_p, co2diag_p, co2prog_p, prcsnw_p
         type(c_ptr), value :: prec_dp_p, snow_dp_p, prec_sh_p, snow_sh_p, prec_sed_p, snow_sed_p, prec_pcw_p, snow_pcw_p
         type(c_ptr), value :: tbot_p, thbot_p, zbot_p, ubot_p, vbot_p, pbot_p, rho_p, qbot_p
         type(c_ptr), value :: precc_p, precl_p, precsc_p, precsl_p
         type(c_ptr), value :: precrl_16O_in_p, precsl_16O_in_p, precrc_16O_in_p, precsc_16O_in_p
         type(c_ptr), value :: precrl_HDO_in_p, precsl_HDO_in_p, precrc_HDO_in_p, precsc_HDO_in_p
         type(c_ptr), value :: precrl_18O_in_p, precsl_18O_in_p, precrc_18O_in_p, precsc_18O_in_p
         type(c_ptr), value :: precrl_16O_out_p, precsl_16O_out_p, precrc_16O_out_p, precsc_16O_out_p
         type(c_ptr), value :: precrl_HDO_out_p, precsl_HDO_out_p, precrc_HDO_out_p, precsc_HDO_out_p
         type(c_ptr), value :: precrl_18O_out_p, precsl_18O_out_p, precrc_18O_out_p, precsc_18O_out_p
      end subroutine cam_export_core_stage_dispatch_codon
   end interface

   !NOTE:  This water tracer code can currently only handle three water tracers/isotopes.
   !It might be beneficial in the future to make this setup more flexible, with
   !a few variables containing all of the needed data [e.g., prec(n), where n is the
   !number of water tracers or isotopes]. - JN.

   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   dummy_1d(:) = 0._r8
   call cam_export_select_impl()

   prec_dp_idx = pbuf_get_index('PREC_DP')
   snow_dp_idx = pbuf_get_index('SNOW_DP')
   prec_sh_idx = pbuf_get_index('PREC_SH')
   snow_sh_idx = pbuf_get_index('SNOW_SH')
   prec_sed_idx = pbuf_get_index('PREC_SED')
   snow_sed_idx = pbuf_get_index('SNOW_SED')
   prec_pcw_idx = pbuf_get_index('PREC_PCW')
   snow_pcw_idx = pbuf_get_index('SNOW_PCW')

   call pbuf_get_field(pbuf, prec_dp_idx, prec_dp)
   call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
   call pbuf_get_field(pbuf, prec_sh_idx, prec_sh)
   call pbuf_get_field(pbuf, snow_sh_idx, snow_sh)
   call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
   call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
   call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
   call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)

  !water tracers/isotopes:
  !----------------------
   exist16 = .false.
   existD  = .false.
   exist18 = .false.
   pass16  = .true.
   passD   = .true.
   pass18  = .true.
   if(trace_water) then
     do m=1, wtrc_nwset !loop over water tracer precip.
       select case(iwspec(wtrc_iatype(m,iwtstrain))) !determine water species
         case(isph2o) !H2O
         !Do nothing. This call may not be needed... -JN
         case(isph216o) !H216O 
           if(pass16) then !first call for H216O?
             exist16 = .true.
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtstrain,m), precrl_16O)
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtstsnow,m), precsl_16O)   
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtcvrain,m), precrc_16O)
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtcvsnow,m), precsc_16O)
             pass16 = .false.
           end if
         case(isphdo) !HD16O
           if(passD) then !first call for HDO?
             existD = .true.
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtstrain,m), precrl_HDO)
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtstsnow,m), precsl_HDO)
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtcvrain,m), precrc_HDO)
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtcvsnow,m), precsc_HDO)
             passD = .false.
           end if
         case(isph218o) !H218O
           if(pass18) then !first call for H218O?
             exist18 = .true.
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtstrain,m), precrl_18O)
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtstsnow,m), precsl_18O)
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtcvrain,m), precrc_18O)
             call pbuf_get_field(pbuf,wtrc_srfpcp_indices(iwtcvsnow,m), precsc_18O)
             pass18 = .false.
           end if
       end select
     end do
   end if
  !-------------------------

   precrl_16O_p = c_loc(dummy_1d(1))
   precsl_16O_p = c_loc(dummy_1d(1))
   precrc_16O_p = c_loc(dummy_1d(1))
   precsc_16O_p = c_loc(dummy_1d(1))
   precrl_HDO_p = c_loc(dummy_1d(1))
   precsl_HDO_p = c_loc(dummy_1d(1))
   precrc_HDO_p = c_loc(dummy_1d(1))
   precsc_HDO_p = c_loc(dummy_1d(1))
   precrl_18O_p = c_loc(dummy_1d(1))
   precsl_18O_p = c_loc(dummy_1d(1))
   precrc_18O_p = c_loc(dummy_1d(1))
   precsc_18O_p = c_loc(dummy_1d(1))
   if (trace_water .and. exist16) then
      precrl_16O_p = c_loc(precrl_16O(1))
      precsl_16O_p = c_loc(precsl_16O(1))
      precrc_16O_p = c_loc(precrc_16O(1))
      precsc_16O_p = c_loc(precsc_16O(1))
   end if
   if (trace_water .and. existD) then
      precrl_HDO_p = c_loc(precrl_HDO(1))
      precsl_HDO_p = c_loc(precsl_HDO(1))
      precrc_HDO_p = c_loc(precrc_HDO(1))
      precsc_HDO_p = c_loc(precsc_HDO(1))
   end if
   if (trace_water .and. exist18) then
      precrl_18O_p = c_loc(precrl_18O(1))
      precsl_18O_p = c_loc(precsl_18O(1))
      precrc_18O_p = c_loc(precrc_18O(1))
      precsc_18O_p = c_loc(precsc_18O(1))
   end if
   co2diag_val = chem_surfvals_get('CO2VMR') * 1.0e+6_r8

   if (cam_export_use_native_impl) then
      do i=1,ncol
         cam_out%tbot(i)  = state%t(i,pver)
         cam_out%thbot(i) = state%t(i,pver) * state%exner(i,pver)
         cam_out%zbot(i)  = state%zm(i,pver)
         cam_out%ubot(i)  = state%u(i,pver)
         cam_out%vbot(i)  = state%v(i,pver)
         cam_out%pbot(i)  = state%pmid(i,pver)
         cam_out%rho(i)   = cam_out%pbot(i)/(rair*cam_out%tbot(i))
         psm1(i,lchnk)    = state%ps(i)
         srfrpdel(i,lchnk)= state%rpdel(i,pver)
      end do
      do m = 1, pcnst
        do i = 1, ncol
           cam_out%qbot(i,m) = state%q(i,pver,m)
        end do
      end do
   else
      call cam_export_log_entered()
      call cam_export_core_stage_dispatch_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnst, c_int64_t), &
           real(rair, c_double), real(mwdry, c_double), real(mwco2, c_double), real(co2diag_val, c_double), &
           merge(1_c_int64_t, 0_c_int64_t, co2_transport()), int(c_i(4), c_int64_t), &
           merge(1_c_int64_t, 0_c_int64_t, trace_water), merge(1_c_int64_t, 0_c_int64_t, exist16), &
           merge(1_c_int64_t, 0_c_int64_t, existD), merge(1_c_int64_t, 0_c_int64_t, exist18), &
           c_loc(state%t), c_loc(state%exner), c_loc(state%zm), c_loc(state%u), &
           c_loc(state%v), c_loc(state%pmid), c_loc(state%q), c_loc(state%ps), c_loc(state%rpdel), &
           c_loc(psm1(1,lchnk)), c_loc(srfrpdel(1,lchnk)), c_loc(cam_out%co2diag(1)), &
           c_loc(cam_out%co2prog(1)), c_loc(prcsnw(1,lchnk)), c_loc(prec_dp(1)), c_loc(snow_dp(1)), &
           c_loc(prec_sh(1)), c_loc(snow_sh(1)), c_loc(prec_sed(1)), c_loc(snow_sed(1)), &
           c_loc(prec_pcw(1)), c_loc(snow_pcw(1)), c_loc(cam_out%tbot(1)), c_loc(cam_out%thbot(1)), &
           c_loc(cam_out%zbot(1)), c_loc(cam_out%ubot(1)), c_loc(cam_out%vbot(1)), c_loc(cam_out%pbot(1)), &
           c_loc(cam_out%rho(1)), c_loc(cam_out%qbot(1,1)), c_loc(cam_out%precc(1)), c_loc(cam_out%precl(1)), &
           c_loc(cam_out%precsc(1)), c_loc(cam_out%precsl(1)), &
           precrl_16O_p, precsl_16O_p, precrc_16O_p, precsc_16O_p, &
           precrl_HDO_p, precsl_HDO_p, precrc_HDO_p, precsc_HDO_p, &
           precrl_18O_p, precsl_18O_p, precrc_18O_p, precsc_18O_p, &
           c_loc(cam_out%precrl_16O(1)), c_loc(cam_out%precsl_16O(1)), &
           c_loc(cam_out%precrc_16O(1)), c_loc(cam_out%precsc_16O(1)), &
           c_loc(cam_out%precrl_HDO(1)), c_loc(cam_out%precsl_HDO(1)), &
           c_loc(cam_out%precrc_HDO(1)), c_loc(cam_out%precsc_HDO(1)), &
           c_loc(cam_out%precrl_18O(1)), c_loc(cam_out%precsl_18O(1)), &
           c_loc(cam_out%precrc_18O(1)), c_loc(cam_out%precsc_18O(1)) &
      )
   end if

   if (cam_export_use_native_impl) then
      cam_out%co2diag(:ncol) = co2diag_val
   end if
   if (cam_export_use_native_impl .and. co2_transport()) then
      do i=1,ncol
         cam_out%co2prog(i) = state%q(i,pver,c_i(4)) * 1.0e+6_r8 *mwdry/mwco2
      end do
   end if
   !
   ! Precipation and snow rates from shallow convection, deep convection and stratiform processes.
   ! Compute total convective and stratiform precipitation and snow rates
   !
   if (cam_export_use_native_impl) then
      do i=1,ncol
         cam_out%precc (i) = prec_dp(i)  + prec_sh(i)
         cam_out%precl (i) = prec_sed(i) + prec_pcw(i)
         cam_out%precsc(i) = snow_dp(i)  + snow_sh(i)
         cam_out%precsl(i) = snow_sed(i) + snow_pcw(i)

         ! jrm These checks should not be necessary if they exist in the parameterizations
         if (cam_out%precc(i) .lt.0._r8) cam_out%precc(i)=0._r8
         if (cam_out%precl(i) .lt.0._r8) cam_out%precl(i)=0._r8
         if (cam_out%precsc(i).lt.0._r8) cam_out%precsc(i)=0._r8
         if (cam_out%precsl(i).lt.0._r8) cam_out%precsl(i)=0._r8
         if (cam_out%precsc(i).gt.cam_out%precc(i)) cam_out%precsc(i)=cam_out%precc(i)
         if (cam_out%precsl(i).gt.cam_out%precl(i)) cam_out%precsl(i)=cam_out%precl(i)
         ! end jrm
         !water tracers/isotopes:
         !----------------------
         if(trace_water) then
           if(exist16) then
             cam_out%precrl_16O(i)  = precrl_16O(i)
             cam_out%precsl_16O(i)  = precsl_16O(i)
             cam_out%precrc_16O(i)  = precrc_16O(i)
             cam_out%precsc_16O(i)  = precsc_16O(i)
           end if
           if(existD) then
             cam_out%precrl_HDO(i)  = precrl_HDO(i)
             cam_out%precsl_HDO(i)  = precsl_HDO(i)
             cam_out%precrc_HDO(i)  = precrc_HDO(i)
             cam_out%precsc_HDO(i)  = precsc_HDO(i)
           end if
           if(exist18) then
             cam_out%precrl_18O(i)  = precrl_18O(i)
             cam_out%precsl_18O(i)  = precsl_18O(i)
             cam_out%precrc_18O(i)  = precrc_18O(i)
             cam_out%precsc_18O(i)  = precsc_18O(i)
           end if
          !negative value prevention:
           if (cam_out%precrl_16O(i) .lt. 0._r8) cam_out%precrl_16O(i)=0._r8
           if (cam_out%precrl_HDO(i) .lt. 0._r8) cam_out%precrl_HDO(i)=0._r8
           if (cam_out%precrl_18O(i) .lt. 0._r8) cam_out%precrl_18O(i)=0._r8
           if (cam_out%precsl_16O(i) .lt. 0._r8) cam_out%precsl_16O(i)=0._r8
           if (cam_out%precsl_HDO(i) .lt. 0._r8) cam_out%precsl_HDO(i)=0._r8
           if (cam_out%precsl_18O(i) .lt. 0._r8) cam_out%precsl_18O(i)=0._r8
           if (cam_out%precrc_16O(i) .lt. 0._r8) cam_out%precrc_16O(i)=0._r8
           if (cam_out%precrc_HDO(i) .lt. 0._r8) cam_out%precrc_HDO(i)=0._r8
           if (cam_out%precrc_18O(i) .lt. 0._r8) cam_out%precrc_18O(i)=0._r8
           if (cam_out%precsc_16O(i) .lt. 0._r8) cam_out%precsc_16O(i)=0._r8
           if (cam_out%precsc_HDO(i) .lt. 0._r8) cam_out%precsc_HDO(i)=0._r8
           if (cam_out%precsc_18O(i) .lt. 0._r8) cam_out%precsc_18O(i)=0._r8
         end if
         !----------------------
      end do
   end if

   ! total snowfall rate: needed by slab ocean model
   if (cam_export_use_native_impl) then
      prcsnw(:ncol,lchnk) = cam_out%precsc(:ncol) + cam_out%precsl(:ncol)
   end if

end subroutine cam_export

end module camsrfexch
