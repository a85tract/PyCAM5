module cloud_fraction

  ! Cloud fraction parameterization.


  use shr_kind_mod,   only: r8 => shr_kind_r8
  use ppgrid,         only: pcols, pver, pverp
  use ref_pres,       only: pref_mid 
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun
  use ref_pres,       only: trop_cloud_top_lev
  use iso_c_binding,  only: c_int64_t

  implicit none
  private
  save

  ! Public interfaces
  public &
     cldfrc_readnl,    &! read cldfrc_nl namelist
     cldfrc_register,  &! add fields to pbuf
     cldfrc_init,      &! Inititialization of cloud_fraction run-time parameters
     cldfrc_getparams, &! public access of tuning parameters
     cldfrc,           &! Computation of cloud fraction
     cldfrc_fice        ! Calculate fraction of condensate in ice phase (radiation partitioning)

  ! Private data
  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  ! Top level
  integer :: top_lev = 1

  ! Physics buffer indices 
  integer :: sh_frac_idx   = 0  
  integer :: dp_frac_idx   = 0 

  ! Namelist variables
  logical  :: cldfrc_freeze_dry           ! switch for Vavrus correction
  logical  :: cldfrc_ice                  ! switch to compute ice cloud fraction
  real(r8) :: cldfrc_rhminl = unset_r8    ! minimum rh for low stable clouds
  real(r8) :: cldfrc_rhminl_adj_land = unset_r8   ! rhminl adjustment for snowfree land
  real(r8) :: cldfrc_rhminh = unset_r8    ! minimum rh for high stable clouds
  real(r8) :: cldfrc_rhminp = unset_r8    ! minimum rh for high stable clouds poleward of 60 degrees
  real(r8) :: cldfrc_rhminp_botmb = 300._r8 ! and pressures less than cldfrc_rhminp_botmb (hPa)
  real(r8) :: cldfrc_sh1    = unset_r8    ! parameter for shallow convection cloud fraction
  real(r8) :: cldfrc_sh2    = unset_r8    ! parameter for shallow convection cloud fraction
  real(r8) :: cldfrc_dp1    = unset_r8    ! parameter for deep convection cloud fraction
  real(r8) :: cldfrc_dp2    = unset_r8    ! parameter for deep convection cloud fraction
  real(r8) :: cldfrc_premit = unset_r8    ! top pressure bound for mid level cloud
  real(r8) :: cldfrc_premib  = unset_r8   ! bottom pressure bound for mid level cloud
  integer  :: cldfrc_iceopt               ! option for ice cloud closure
                                          ! 1=wang & sassen 2=schiller (iciwc)
                                          ! 3=wood & field, 4=Wilson (based on smith)
  real(r8) :: cldfrc_icecrit = unset_r8   ! Critical RH for ice clouds in Wilson & Ballard closure (smaller = more ice clouds)

  real(r8) :: rhminl             ! set from namelist input cldfrc_rhminl
  real(r8) :: rhminl_adj_land    ! set from namelist input cldfrc_rhminl_adj_land
  real(r8) :: rhminh             ! set from namelist input cldfrc_rhminh
  real(r8) :: rhminp             ! set from namelist input cldfrc_rhminp
  real(r8) :: sh1, sh2           ! set from namelist input cldfrc_sh1, cldfrc_sh2
  real(r8) :: dp1,dp2            ! set from namelist input cldfrc_dp1, cldfrc_dp2
  real(r8) :: premit             ! set from namelist input cldfrc_premit
  real(r8) :: premib             ! set from namelist input cldfrc_premib
  integer  :: iceopt             ! set from namelist input cldfrc_iceopt
  real(r8) :: icecrit            ! set from namelist input cldfrc_icecrit

  ! constants
  real(r8), parameter :: pnot = 1.e5_r8         ! reference pressure
  real(r8), parameter :: lapse = 6.5e-3_r8      ! U.S. Standard Atmosphere lapse rate
  real(r8), parameter :: pretop = 1.0e2_r8      ! pressure bounding high cloud

  integer count

  logical :: inversion_cld_off    ! Turns off stratification-based cld frc

  integer :: k700   ! model level nearest 700 mb
  logical :: use_native_cldfrc_fice_impl = .false.
  logical :: cldfrc_fice_impl_selected = .false.
  logical :: use_native_cldfrc_convective_cover_impl = .false.
  logical :: cldfrc_convective_cover_impl_selected = .false.
  logical :: use_native_cldfrc_state_init_impl = .false.
  logical :: cldfrc_state_init_impl_selected = .false.
  logical :: use_native_cldfrc_layer_rh_impl = .false.
  logical :: cldfrc_layer_rh_impl_selected = .false.
  logical :: use_native_cldfrc_ice_wilson_impl = .false.
  logical :: cldfrc_ice_wilson_impl_selected = .false.
  logical :: use_native_cldfrc_total_cloud_impl = .false.
  logical :: cldfrc_total_cloud_impl_selected = .false.
  logical :: use_native_cldfrc_batch_impl = .false.
  logical :: cldfrc_batch_impl_selected = .false.
  logical :: use_native_cldfrc_getparams_impl = .false.
  logical :: cldfrc_getparams_impl_selected = .false.
  logical :: cldfrc_batch_entered_logged = .false.

  interface
     subroutine cldfrc_batch_dispatch_codon(stage_c, ncol_c, pcols_c, pver_c, top_lev_c, flag_c, &
          scalar1_c, scalar2_c, scalar3_c, scalar4_c, scalar5_c, scalar6_c, scalar7_c, scalar8_c, scalar9_c, &
          p1_p, p2_p, p3_p, p4_p, p5_p, p6_p, p7_p, p8_p, p9_p, p10_p) &
          bind(c, name="cldfrc_batch_dispatch_codon")
        use iso_c_binding, only: c_double, c_int64_t, c_ptr
        integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, top_lev_c, flag_c
        real(c_double), value :: scalar1_c, scalar2_c, scalar3_c, scalar4_c, scalar5_c, scalar6_c
        real(c_double), value :: scalar7_c, scalar8_c, scalar9_c
        type(c_ptr), value :: p1_p, p2_p, p3_p, p4_p, p5_p, p6_p, p7_p, p8_p, p9_p, p10_p
     end subroutine cldfrc_batch_dispatch_codon

     subroutine cldfrc_getparams_codon(flags_c, rhminl_c, rhminl_adj_land_c, rhminh_c, rhminp_c, &
          premit_c, premib_c, iceopt_c, icecrit_c, rhminl_p, rhminl_adj_land_p, rhminh_p, rhminp_p, &
          premit_p, premib_p, iceopt_p, icecrit_p) bind(c, name="cldfrc_getparams_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: flags_c, iceopt_c
       real(c_double), value :: rhminl_c, rhminl_adj_land_c, rhminh_c, rhminp_c
       real(c_double), value :: premit_c, premib_c, icecrit_c
       type(c_ptr), value :: rhminl_p, rhminl_adj_land_p, rhminh_p, rhminp_p
       type(c_ptr), value :: premit_p, premib_p, iceopt_p, icecrit_p
     end subroutine cldfrc_getparams_codon
     function cldfrc_register_codon(flag_c) result(out_c) bind(c, name="cldfrc_register_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: out_c
     end function cldfrc_register_codon
  end interface

!================================================================================================
  contains
!================================================================================================

subroutine cldfrc_fice_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cldfrc_fice_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CLDFRC_FICE_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_cldfrc_fice_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_cldfrc_fice_impl = .false.
   end if

   cldfrc_fice_impl_selected = .true.

   if (masterproc) then
      if (use_native_cldfrc_fice_impl) then
         write(iulog,*) 'cldfrc_fice implementation = native'
      else
         write(iulog,*) 'cldfrc_fice implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine cldfrc_fice_select_impl

!================================================================================================

subroutine cldfrc_convective_cover_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cldfrc_convective_cover_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CLDFRC_CONVECTIVE_COVER_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_cldfrc_convective_cover_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_cldfrc_convective_cover_impl = .false.
   end if

   cldfrc_convective_cover_impl_selected = .true.

   if (masterproc) then
      if (use_native_cldfrc_convective_cover_impl) then
         write(iulog,*) 'cldfrc_convective_cover implementation = native'
      else
         write(iulog,*) 'cldfrc_convective_cover implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine cldfrc_convective_cover_select_impl

!================================================================================================

subroutine cldfrc_state_init_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cldfrc_state_init_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CLDFRC_STATE_INIT_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_cldfrc_state_init_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_cldfrc_state_init_impl = .false.
   end if

   cldfrc_state_init_impl_selected = .true.

   if (masterproc) then
      if (use_native_cldfrc_state_init_impl) then
         write(iulog,*) 'cldfrc_state_init implementation = native'
      else
         write(iulog,*) 'cldfrc_state_init implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine cldfrc_state_init_select_impl

!================================================================================================

subroutine cldfrc_layer_rh_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cldfrc_layer_rh_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CLDFRC_LAYER_RH_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_cldfrc_layer_rh_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_cldfrc_layer_rh_impl = .false.
   end if

   cldfrc_layer_rh_impl_selected = .true.

   if (masterproc) then
      if (use_native_cldfrc_layer_rh_impl) then
         write(iulog,*) 'cldfrc_layer_rh implementation = native'
      else
         write(iulog,*) 'cldfrc_layer_rh implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine cldfrc_layer_rh_select_impl

!================================================================================================

subroutine cldfrc_ice_wilson_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cldfrc_ice_wilson_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CLDFRC_ICE_WILSON_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_cldfrc_ice_wilson_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_cldfrc_ice_wilson_impl = .false.
   end if

   cldfrc_ice_wilson_impl_selected = .true.

   if (masterproc) then
      if (use_native_cldfrc_ice_wilson_impl) then
         write(iulog,*) 'cldfrc_ice_wilson implementation = native'
      else
         write(iulog,*) 'cldfrc_ice_wilson implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine cldfrc_ice_wilson_select_impl

!================================================================================================

subroutine cldfrc_total_cloud_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cldfrc_total_cloud_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CLDFRC_TOTAL_CLOUD_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_cldfrc_total_cloud_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_cldfrc_total_cloud_impl = .false.
   end if

   cldfrc_total_cloud_impl_selected = .true.

   if (masterproc) then
      if (use_native_cldfrc_total_cloud_impl) then
         write(iulog,*) 'cldfrc_total_cloud implementation = native'
      else
         write(iulog,*) 'cldfrc_total_cloud implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine cldfrc_total_cloud_select_impl

!================================================================================================

subroutine cldfrc_batch_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line

   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('CLDFRC_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine cldfrc_batch_append_proof

!================================================================================================

subroutine cldfrc_batch_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cldfrc_batch_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CLDFRC_BATCH_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_cldfrc_batch_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_cldfrc_batch_impl = .false.
   end if

   cldfrc_batch_impl_selected = .true.

   if (masterproc) then
      if (use_native_cldfrc_batch_impl) then
         write(iulog,*) 'cldfrc_batch implementation = native'
         call cldfrc_batch_append_proof('cldfrc_batch selector entered implementation = native')
      else
         write(iulog,*) 'cldfrc_batch implementation = codon'
         call cldfrc_batch_append_proof('cldfrc_batch selector entered implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine cldfrc_batch_select_impl

!================================================================================================

subroutine cldfrc_batch_log_entered()

   if (cldfrc_batch_entered_logged) return
   cldfrc_batch_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'cldfrc_batch entered (unified stage-dispatch state/layer/ice/convective/total/fice direct = codon)'
      call cldfrc_batch_append_proof( &
           'cldfrc_batch entered (unified stage-dispatch state/layer/ice/convective/total/fice direct = codon)')
      call flush(iulog)
   end if

end subroutine cldfrc_batch_log_entered

!================================================================================================

subroutine cldfrc_getparams_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cldfrc_getparams_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CLDFRC_GETPARAMS_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_cldfrc_getparams_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_cldfrc_getparams_impl = .false.
   end if

   cldfrc_getparams_impl_selected = .true.

   if (masterproc) then
      if (use_native_cldfrc_getparams_impl) then
         write(iulog,*) 'cldfrc_getparams implementation = native'
      else
         write(iulog,*) 'cldfrc_getparams implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine cldfrc_getparams_select_impl

!================================================================================================

subroutine cldfrc_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cldfrc_readnl'

   namelist /cldfrc_nl/ cldfrc_freeze_dry,      cldfrc_ice,    cldfrc_rhminl, &
                        cldfrc_rhminl_adj_land, cldfrc_rhminh, cldfrc_sh1,    &
                        cldfrc_rhminp,          cldfrc_rhminp_botmb, &
                        cldfrc_sh2,             cldfrc_dp1,    cldfrc_dp2,    &
                        cldfrc_premit,          cldfrc_premib, cldfrc_iceopt, &
                        cldfrc_icecrit
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cldfrc_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cldfrc_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      rhminl = cldfrc_rhminl
      rhminl_adj_land = cldfrc_rhminl_adj_land
      rhminh = cldfrc_rhminh
      rhminp = cldfrc_rhminp
      sh1    = cldfrc_sh1
      sh2    = cldfrc_sh2
      dp1    = cldfrc_dp1
      dp2    = cldfrc_dp2
      premit = cldfrc_premit
      premib  = cldfrc_premib
      iceopt  = cldfrc_iceopt
      icecrit = cldfrc_icecrit

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(cldfrc_freeze_dry, 1, mpilog, 0, mpicom)
   call mpibcast(cldfrc_ice,        1, mpilog, 0, mpicom)
   call mpibcast(rhminl,            1, mpir8,  0, mpicom)
   call mpibcast(rhminl_adj_land,   1, mpir8,  0, mpicom)
   call mpibcast(rhminh,            1, mpir8,  0, mpicom)
   call mpibcast(rhminp,            1, mpir8,  0, mpicom)
   call mpibcast(sh1   ,            1, mpir8,  0, mpicom)
   call mpibcast(sh2   ,            1, mpir8,  0, mpicom)
   call mpibcast(dp1   ,            1, mpir8,  0, mpicom)
   call mpibcast(dp2   ,            1, mpir8,  0, mpicom)
   call mpibcast(premit,            1, mpir8,  0, mpicom)
   call mpibcast(premib,            1, mpir8,  0, mpicom)
   call mpibcast(iceopt,            1, mpiint, 0, mpicom)
   call mpibcast(icecrit,           1, mpir8,  0, mpicom)
#endif

end subroutine cldfrc_readnl

!================================================================================================

subroutine cldfrc_register

   ! Register fields in the physics buffer.

   use physics_buffer, only : pbuf_add_field, dtype_r8

   !-----------------------------------------------------------------------
   if (cldfrc_register_codon(1_c_int64_t) == 0_c_int64_t) return

   call pbuf_add_field('SH_FRAC', 'physpkg', dtype_r8, (/pcols,pver/), sh_frac_idx) 
   call pbuf_add_field('DP_FRAC', 'physpkg', dtype_r8, (/pcols,pver/), dp_frac_idx) 

end subroutine cldfrc_register

!================================================================================================

subroutine cldfrc_getparams(rhminl_out, rhminl_adj_land_out, rhminh_out,  premit_out, &
                            rhminp_out, premib_out, iceopt_out, icecrit_out)
   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
!-----------------------------------------------------------------------
! Purpose: Return cldfrc tuning parameters
!-----------------------------------------------------------------------

   real(r8), target, intent(out), optional :: rhminl_out
   real(r8), target, intent(out), optional :: rhminl_adj_land_out
   real(r8), target, intent(out), optional :: rhminh_out
   real(r8), target, intent(out), optional :: rhminp_out
   real(r8), target, intent(out), optional :: premit_out
   real(r8), target, intent(out), optional :: premib_out
   integer,  target, intent(out), optional :: iceopt_out
   real(r8), target, intent(out), optional :: icecrit_out

   integer(c_int64_t) :: flags
   type(c_ptr) :: rhminl_p, rhminl_adj_land_p, rhminh_p, rhminp_p
   type(c_ptr) :: premit_p, premib_p, iceopt_p, icecrit_p

   call cldfrc_getparams_select_impl()

   if (use_native_cldfrc_getparams_impl) then
      if ( present(rhminl_out) )      rhminl_out = rhminl
      if ( present(rhminl_adj_land_out) ) rhminl_adj_land_out = rhminl_adj_land
      if ( present(rhminh_out) )      rhminh_out = rhminh
      if ( present(rhminp_out) )      rhminp_out = rhminp
      if ( present(premit_out) )      premit_out = premit
      if ( present(premib_out) )      premib_out  = premib
      if ( present(iceopt_out) )      iceopt_out  = iceopt
      if ( present(icecrit_out) )     icecrit_out = icecrit
      return
   end if

   flags = 0_c_int64_t
   rhminl_p = c_null_ptr
   rhminl_adj_land_p = c_null_ptr
   rhminh_p = c_null_ptr
   rhminp_p = c_null_ptr
   premit_p = c_null_ptr
   premib_p = c_null_ptr
   iceopt_p = c_null_ptr
   icecrit_p = c_null_ptr

   if ( present(rhminl_out) ) then
      flags = flags + 1_c_int64_t
      rhminl_p = c_loc(rhminl_out)
   end if
   if ( present(rhminl_adj_land_out) ) then
      flags = flags + 2_c_int64_t
      rhminl_adj_land_p = c_loc(rhminl_adj_land_out)
   end if
   if ( present(rhminh_out) ) then
      flags = flags + 4_c_int64_t
      rhminh_p = c_loc(rhminh_out)
   end if
   if ( present(rhminp_out) ) then
      flags = flags + 8_c_int64_t
      rhminp_p = c_loc(rhminp_out)
   end if
   if ( present(premit_out) ) then
      flags = flags + 16_c_int64_t
      premit_p = c_loc(premit_out)
   end if
   if ( present(premib_out) ) then
      flags = flags + 32_c_int64_t
      premib_p = c_loc(premib_out)
   end if
   if ( present(iceopt_out) ) then
      flags = flags + 64_c_int64_t
      iceopt_p = c_loc(iceopt_out)
   end if
   if ( present(icecrit_out) ) then
      flags = flags + 128_c_int64_t
      icecrit_p = c_loc(icecrit_out)
   end if

   call cldfrc_getparams_codon(flags, real(rhminl, c_double), real(rhminl_adj_land, c_double), &
        real(rhminh, c_double), real(rhminp, c_double), real(premit, c_double), real(premib, c_double), &
        int(iceopt, c_int64_t), real(icecrit, c_double), rhminl_p, rhminl_adj_land_p, rhminh_p, &
        rhminp_p, premit_p, premib_p, iceopt_p, icecrit_p)

end subroutine cldfrc_getparams

!===============================================================================

subroutine cldfrc_init

   ! Initialize cloud fraction run-time parameters

   use cam_history,   only:  phys_decomp, addfld
   use dycore,        only:  dycore_is, get_resolution
   use phys_control,  only:  phys_getopts

   ! horizontal grid specifier
   character(len=32) :: hgrid

   ! query interfaces for scheme settings
   character(len=16) :: shallow_scheme, eddy_scheme, macrop_scheme

   integer :: k
   !-----------------------------------------------------------------------------

   call phys_getopts(shallow_scheme_out = shallow_scheme ,&
                     eddy_scheme_out    = eddy_scheme    ,&
                     macrop_scheme_out  = macrop_scheme  )

   ! Limit CAM5 cloud physics to below top cloud level.
   if (macrop_scheme /= "rk") top_lev = trop_cloud_top_lev

   hgrid = get_resolution()

   ! Turn off inversion_cld if any UW PBL scheme is being used
   if ( (eddy_scheme .eq. 'diag_TKE' ) .or. (shallow_scheme .eq.  'UW' )) then
      inversion_cld_off = .true.
   else
      inversion_cld_off = .false.
   endif

   if ( masterproc ) then 
      write(iulog,*)'tuning parameters cldfrc_init: inversion_cld_off',inversion_cld_off
      write(iulog,*)'tuning parameters cldfrc_init: dp1',dp1,'dp2',dp2,'sh1',sh1,'sh2',sh2
      if (shallow_scheme .ne. 'UW' ) then
         write(iulog,*)'tuning parameters cldfrc_init: rhminl',rhminl,'rhminl_adj_land',rhminl_adj_land, &
                       'rhminh',rhminh,'premit',premit,'premib',premib
         write(iulog,*)'tuning parameters cldfrc_init: iceopt',iceopt,'icecrit',icecrit
      endif
   endif

   if (pref_mid(top_lev) > 7.e4_r8) &
        call endrun ('cldfrc_init: model levels bracketing 700 mb not found')

   ! Find vertical level nearest 700 mb.
   k700 = minloc(abs(pref_mid(top_lev:pver) - 7.e4_r8), 1)

   if (masterproc) then
      write(iulog,*)'cldfrc_init: model level nearest 700 mb is',k700,'which is',pref_mid(k700),'pascals'
   end if

   call addfld ('SH_CLD   ', 'fraction', pver, 'A', 'Shallow convective cloud cover'                          ,phys_decomp)
   call addfld ('DP_CLD   ', 'fraction', pver, 'A', 'Deep convective cloud cover'                             ,phys_decomp)

end subroutine cldfrc_init

!===============================================================================

subroutine cldfrc(lchnk   ,ncol    , pbuf,  &
       pmid    ,temp    ,q       ,omga    , phis, &
       shfrc   ,use_shfrc, &
       cloud   ,rhcloud, clc     ,pdel    , &
       cmfmc   ,cmfmc2  ,landfrac,snowh   ,concld  ,cldst   , &
       ts      ,sst     ,ps      ,zdu     ,ocnfrac ,&
       rhu00   ,cldice  ,icecldf ,liqcldf ,relhum  ,dindex )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute cloud fraction 
    ! 
    ! 
    ! Method: 
    ! This calculate cloud fraction using a relative humidity threshold
    ! The threshold depends upon pressure, and upon the presence or absence 
    ! of convection as defined by a reasonably large vertical mass flux 
    ! entering that layer from below.
    ! 
    ! Author: Many. Last modified by Jim McCaa
    ! 
    !-----------------------------------------------------------------------
    use cam_history,   only: outfld
    use physconst,     only: cappa, gravit, rair, tmelt
    use wv_saturation, only: qsat, qsat_water, svp_ice
    use phys_grid,     only: get_rlat_all_p, get_rlon_all_p
    use dycore,        only: dycore_is, get_resolution

   
!RBN - Need this to write shallow,deep fraction to phys buffer.
!PJR - we should probably make seperate modules for determining convective
!      clouds and make this one just responsible for relative humidity clouds
    
    use physics_buffer, only: physics_buffer_desc, pbuf_get_field

    ! Arguments
    integer, intent(in) :: lchnk                  ! chunk identifier
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    integer, intent(in) :: dindex                 ! 0 or 1 to perturb rh
    
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(in) :: pmid(pcols,pver)      ! midpoint pressures
    real(r8), intent(in) :: temp(pcols,pver)      ! temperature
    real(r8), intent(in) :: q(pcols,pver)         ! specific humidity
    real(r8), intent(in) :: omga(pcols,pver)      ! vertical pressure velocity
    real(r8), intent(in) :: cmfmc(pcols,pverp)    ! convective mass flux--m sub c
    real(r8), intent(in) :: cmfmc2(pcols,pverp)   ! shallow convective mass flux--m sub c
    real(r8), intent(in) :: snowh(pcols)          ! snow depth (liquid water equivalent)
    real(r8), intent(in) :: pdel(pcols,pver)      ! pressure depth of layer
    real(r8), intent(in) :: landfrac(pcols)       ! Land fraction
    real(r8), intent(in) :: ocnfrac(pcols)        ! Ocean fraction
    real(r8), intent(in) :: ts(pcols)             ! surface temperature
    real(r8), intent(in) :: sst(pcols)            ! sea surface temperature
    real(r8), intent(in) :: ps(pcols)             ! surface pressure
    real(r8), intent(in) :: zdu(pcols,pver)       ! detrainment rate from deep convection
    real(r8), intent(in) :: phis(pcols)           ! surface geopotential
    real(r8), intent(in) :: shfrc(pcols,pver)     ! cloud fraction from convect_shallow
    real(r8), intent(in) :: cldice(pcols,pver)    ! cloud ice mixing ratio
    logical,  intent(in)  :: use_shfrc

    ! Output arguments
    real(r8), intent(out) :: cloud(pcols,pver)     ! cloud fraction
    real(r8), intent(out) :: rhcloud(pcols,pver)   ! cloud fraction
    real(r8), intent(out) :: clc(pcols)            ! column convective cloud amount
    real(r8), intent(out) :: cldst(pcols,pver)     ! cloud fraction
    real(r8), intent(out) :: rhu00(pcols,pver)     ! RH threshold for cloud
    real(r8), intent(out) :: relhum(pcols,pver)    ! RH 
    real(r8), intent(out) :: icecldf(pcols,pver)   ! ice cloud fraction
    real(r8), intent(out) :: liqcldf(pcols,pver)   ! liquid cloud fraction (combined into cloud)

    !---------------------------Local workspace-----------------------------
    !
    real(r8) concld(pcols,pver)    ! convective cloud cover
    real(r8) cld                   ! intermediate scratch variable (low cld)
    real(r8) dthdpmn(pcols)         ! most stable lapse rate below 750 mb
    real(r8) dthdp                 ! lapse rate (intermediate variable)
    real(r8) es(pcols,pver)        ! saturation vapor pressure
    real(r8) qs(pcols,pver)        ! saturation specific humidity
    real(r8) rhwght                ! weighting function for rhlim transition
    real(r8) rh(pcols,pver)        ! relative humidity
    real(r8) rhdif                 ! intermediate scratch variable
    real(r8) strat                 ! intermediate scratch variable
    real(r8) theta(pcols,pver)     ! potential temperature
    real(r8) rhlim                 ! local rel. humidity threshold estimate
    real(r8) coef1                 ! coefficient to convert mass flux to mb/d
    real(r8) clrsky(pcols)         ! temporary used in random overlap calc
    real(r8) rpdeli(pcols,pver-1) ! 1./(pmid(k+1)-pmid(k))
    real(r8) rhpert                !the specified perturbation to rh

    real(r8), pointer, dimension(:,:) :: deepcu      ! deep convection cloud fraction
    real(r8), pointer, dimension(:,:) :: shallowcu   ! shallow convection cloud fraction

    integer i, ierror, k           ! column, level indices
    integer kp1, ifld
    integer kdthdp(pcols)
    integer numkcld                ! number of levels in which to allow clouds

    !  In Cloud Ice Content variables
    real(r8) :: a,b,c,as,bs,cs        !fit parameters
    real(r8) :: Kc                    !constant for ice cloud calc (wood & field)
    real(r8) :: ttmp                  !limited temperature
    real(r8) :: icicval               !empirical iwc value
    real(r8) :: rho                   !local air density
    real(r8) :: esl(pcols,pver)       !liq sat vapor pressure
    real(r8) :: esi(pcols,pver)       !ice sat vapor pressure
    real(r8) :: ncf,phi               !Wilson and Ballard parameters

    real(r8) thetas(pcols)                    ! ocean surface potential temperature
    real(r8) :: clat(pcols)                   ! current latitudes(radians)
    real(r8) :: clon(pcols)                   ! current longitudes(radians)

    ! Statement functions
    logical land
    land(i) = nint(landfrac(i)) == 1

    call get_rlat_all_p(lchnk, ncol, clat)
    call get_rlon_all_p(lchnk, ncol, clon)

    call pbuf_get_field(pbuf, sh_frac_idx, shallowcu )
    call pbuf_get_field(pbuf, dp_frac_idx, deepcu )

    ! Initialise cloud fraction
    shallowcu = 0._r8
    deepcu    = 0._r8

    !==================================================================================
    ! PHILOSOPHY OF PRESENT IMPLEMENTATION
    !++ag ice3
    ! Modification to philosophy for ice supersaturation
    ! philosophy below is based on RH water only. This is 'liquid condensation'
    ! or liquid cloud (even though it will freeze immediately to ice)
    ! The idea is that the RH limits for condensation are strict only for
    ! water saturation
    !
    ! Ice clouds are formed by explicit parameterization of ice nucleation. 
    ! Closure for ice cloud fraction is done on available cloud ice, such that
    ! the in-cloud ice content matches an empirical fit
    ! thus, icecldf = min(cldice/icicval,1) where icicval = f(temp,cldice,numice)
    ! for a first cut, icicval=f(temp) only.
    ! Combined cloud fraction is maximum overlap  cloud=max(1,max(icecldf,liqcldf))
    ! No dA/dt term for ice?
    !--ag
    !
    ! There are three co-existing cloud types: convective, inversion related low-level
    ! stratocumulus, and layered cloud (based on relative humidity).  Layered and 
    ! stratocumulus clouds do not compete with convective cloud for which one creates 
    ! the most cloud.  They contribute collectively to the total grid-box average cloud 
    ! amount.  This is reflected in the way in which the total cloud amount is evaluated 
    ! (a sum as opposed to a logical "or" operation)
    !
    !==================================================================================
    ! set defaults for rhu00
    rhu00(:,:) = 2.0_r8
    ! define rh perturbation in order to estimate rhdfda
    rhpert = 0.01_r8 

    !set Wang and Sassen IWC paramters
    a=26.87_r8
    b=0.569_r8
    c=0.002892_r8
    !set schiller parameters
    as=-68.4202_r8
    bs=0.983917_r8
    cs=2.81795_r8
    !set wood and field paramters...
    Kc=75._r8

    ! Evaluate potential temperature and relative humidity
    ! If not computing ice cloud fraction then hybrid RH, if MG then water RH
    if ( cldfrc_ice ) then
       call qsat_water(temp(1:ncol,top_lev:pver), pmid(1:ncol,top_lev:pver), &
            esl(1:ncol,top_lev:pver), qs(1:ncol,top_lev:pver))

       esi(1:ncol,top_lev:pver) = svp_ice(temp(1:ncol,top_lev:pver))
    else
       call qsat(temp(1:ncol,top_lev:pver), pmid(1:ncol,top_lev:pver), &
            es(1:ncol,top_lev:pver), qs(1:ncol,top_lev:pver))
    endif

    do k=top_lev,pver
       theta(:ncol,k)    = temp(:ncol,k)*(pnot/pmid(:ncol,k))**cappa
    end do

    call cldfrc_state_init(ncol, q, qs, relhum, rh, cloud, icecldf, liqcldf, rhcloud, cldst, concld, dindex, rhpert)

    ! Initialize other temporary variables
    ierror = 0
    do i=1,ncol
       ! Adjust thetas(i) in the presence of non-zero ocean heights.
       ! This reduces the temperature for positive heights according to a standard lapse rate.
       if(ocnfrac(i).gt.0.01_r8) thetas(i)  = &
            ( sst(i) - lapse * phis(i) / gravit) * (pnot/ps(i))**cappa
       if(ocnfrac(i).gt.0.01_r8.and.sst(i).lt.260._r8) ierror = i
       clc(i) = 0.0_r8
    end do
    coef1 = gravit*864.0_r8    ! conversion to millibars/day

    if (ierror > 0) then
       write(iulog,*) 'COLDSST: encountered in cldfrc:', lchnk,ierror,ocnfrac(ierror),sst(ierror)
    endif

    do k=top_lev,pver-1
       rpdeli(:ncol,k) = 1._r8/(pmid(:ncol,k+1) - pmid(:ncol,k))
    end do

    !
    ! Estimate of local convective cloud cover based on convective mass flux
    ! Modify local large-scale relative humidity to account for presence of 
    ! convective cloud when evaluating relative humidity based layered cloud amount
    !
    concld(:ncol,top_lev:pver) = 0.0_r8
    !
    ! cloud mass flux in SI units of kg/m2/s; should produce typical numbers of 20%
    ! shallow and deep convective cloudiness are evaluated separately (since processes
    ! are evaluated separately) and summed
    !   
#ifndef PERGRO
    call cldfrc_convective_cover(ncol, use_shfrc, shfrc, cmfmc, cmfmc2, shallowcu, deepcu, concld, rh)
#endif
    !==================================================================================
    !
    !          ****** Compute layer cloudiness ******
    !
    !====================================================================
    ! Begin the evaluation of layered cloud amount based on (modified) RH 
    !====================================================================
    !
    numkcld = pver
    call cldfrc_layer_rh(ncol, landfrac, snowh, clat, pmid, pref_mid, q, rh, rhcloud, rhu00)

    if (cldfrc_ice .and. iceopt > 3) then
       call cldfrc_ice_wilson(ncol, cldice, qs, rhcloud, icecldf, liqcldf, cloud)
    else
       do k=top_lev+1,numkcld
          do i=1,ncol
             if (cldfrc_ice) then

                ! Evaluate ice cloud fraction based on in-cloud ice content

                !--------ICE CLOUD OPTION 1--------Wang & Sassen 2002
                !         Evaluate desired in-cloud water content
                !               icicval = f(temp,cldice,numice)
                !         Start with a function of temperature.
                !         Wang & Sassen 2002 (JAS), based on ARM site MMCR (midlat cirrus)
                !           parameterization valid for 203-253K
                !           icival > 0 for t>195K
                if (iceopt.lt.3) then
                   if (iceopt.eq.1) then
                      ttmp=max(195._r8,min(temp(i,k),253._r8)) - 273.16_r8
                      icicval=a + b * ttmp + c * ttmp**2._r8
                      !convert units
                      rho=pmid(i,k)/(rair*temp(i,k))
                      icicval= icicval * 1.e-6_r8 / rho
                   else
                      !--------ICE CLOUD OPTION 2--------Schiller 2008 (JGR)
                      !          Use a curve based on FISH measurements in
                      !          tropics, mid-lats and arctic. Curve is for 180-250K (raise to 273K?)
                      !          use median all flights

                      ttmp=max(190._r8,min(temp(i,k),273.16_r8))
                      icicval = 10._r8 **(as * bs**ttmp + cs)
                      !convert units from ppmv to kg/kg
                      icicval= icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
                   endif
                   !set icecldfraction  for OPTION 1 or OPTION2
                   icecldf(i,k) =  max(0._r8,min(cldice(i,k)/icicval,1._r8))

                else if (iceopt.eq.3) then

                   !--------ICE CLOUD OPTION 3--------Wood & Field 2000 (JAS)
                   ! eq 6: cloud fraction = 1 - exp (-K * qc/qsati)
        
                   icecldf(i,k)=1._r8 - exp(-Kc*cldice(i,k)/(qs(i,k)*(esi(i,k)/esl(i,k))))
                   icecldf(i,k)=max(0._r8,min(icecldf(i,k),1._r8))
                else
                   !--------ICE CLOUD OPTION 4--------Wilson and ballard 1999
                   ! inversion of smith....
                   !       ncf = cldice / ((1-RHcrit)*qs)
                   ! then a function of ncf....
                   ncf =cldice(i,k)/((1._r8 - icecrit)*qs(i,k))
                   if (ncf.le.0._r8) then
                      icecldf(i,k)=0._r8
                   else if (ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8) then
                      icecldf(i,k)=0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
                   else if (ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8) then
                      phi=(acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
                      icecldf(i,k)=(1._r8 - 4._r8 * cos(phi) * cos(phi))
                   else
                      icecldf(i,k)=1._r8
                   endif
                   icecldf(i,k)=max(0._r8,min(icecldf(i,k),1._r8))
                endif

                !TEST: if ice present, icecldf=1.
                !          if (cldice(i,k).ge.1.e-8_r8) then
                !             icecldf(i,k) = 0.99_r8
                !          endif

                !!          if ((cldice(i,k) .gt. icicval) .or. ((cldice(i,k) .gt. 0._r8) .and. (icecldf(i,k) .eq. 0._r8))) then
                !          if (cldice(i,k) .gt. 1.e-8_r8) then
                !             write(iulog,*) 'i,k,pmid,rho,t,cldice,icicval,icecldf,rhcloud: ', &
                !                i,k,pmid(i,k),rho,temp(i,k),cldice(i,k),icicval,icecldf(i,k),rhcloud(i,k)
                !          endif

                !         Combine ice and liquid cloud fraction assuming maximum overlap.
                ! Combined cloud fraction is maximum overlap
                !          cloud(i,k)=min(1._r8,max(icecldf(i,k),rhcloud(i,k)))

                liqcldf(i,k)=(1._r8 - icecldf(i,k))* rhcloud(i,k)
                cloud(i,k)=liqcldf(i,k) + icecldf(i,k)
             else
                ! For RK microphysics
                cloud(i,k) = rhcloud(i,k)
             end if
          end do
       end do
    end if
    !
    ! Add in the marine strat
    ! MARINE STRATUS SHOULD BE A SPECIAL CASE OF LAYERED CLOUD
    ! CLOUD CURRENTLY CONTAINS LAYERED CLOUD DETERMINED BY RH CRITERIA
    ! TAKE THE MAXIMUM OF THE DIAGNOSED LAYERED CLOUD OR STRATOCUMULUS
    !
    !===================================================================================
    !
    !  SOME OBSERVATIONS ABOUT THE FOLLOWING SECTION OF CODE (missed in earlier look)
    !  K700 IS SET AS A CONSTANT BASED ON HYBRID COORDINATE: IT DOES NOT DEPEND ON 
    !  LOCAL PRESSURE; THERE IS NO PRESSURE RAMP => LOOKS LEVEL DEPENDENT AND 
    !  DISCONTINUOUS IN SPACE (I.E., STRATUS WILL END SUDDENLY WITH NO TRANSITION)
    !
    !  IT APPEARS THAT STRAT IS EVALUATED ACCORDING TO KLEIN AND HARTMANN; HOWEVER,
    !  THE ACTUAL STRATUS AMOUNT (CLDST) APPEARS TO DEPEND DIRECTLY ON THE RH BELOW
    !  THE STRONGEST PART OF THE LOW LEVEL INVERSION.  
    !PJR answers: 1) the rh limitation is a physical/mathematical limitation
    !             cant have more cloud than there is RH
    !             allowed the cloud to exist two layers below the inversion
    !             because the numerics frequently make 50% relative humidity
    !             in level below the inversion which would allow no cloud
    !             2) since  the cloud is only allowed over ocean, it should
    !             be very insensitive to surface pressure (except due to 
    !             spectral ringing, which also causes so many other problems
    !             I didnt worry about it.
    !
    !==================================================================================
    if (.not.inversion_cld_off) then
    !
    ! Find most stable level below 750 mb for evaluating stratus regimes
    !
    do i=1,ncol
       ! Nothing triggers unless a stability greater than this minimum threshold is found
       dthdpmn(i) = -0.125_r8
       kdthdp(i) = 0
    end do
    !
    do k=top_lev+1,pver
       do i=1,ncol
          if (pmid(i,k) >= premib .and. ocnfrac(i).gt. 0.01_r8) then
             ! I think this is done so that dtheta/dp is in units of dg/mb (JJH)
             dthdp = 100.0_r8*(theta(i,k) - theta(i,k-1))*rpdeli(i,k-1)
             if (dthdp < dthdpmn(i)) then
                dthdpmn(i) = dthdp
                kdthdp(i) = k     ! index of interface of max inversion
             end if
          end if
       end do
    end do

    ! Also check between the bottom layer and the surface
    ! Only perform this check if the criteria were not met above

    do i = 1,ncol
       if ( kdthdp(i) .eq. 0 .and. ocnfrac(i).gt.0.01_r8) then
          dthdp = 100.0_r8 * (thetas(i) - theta(i,pver)) / (ps(i)-pmid(i,pver))
          if (dthdp < dthdpmn(i)) then
             dthdpmn(i) = dthdp
             kdthdp(i) = pver     ! index of interface of max inversion
          endif
       endif
    enddo

    do i=1,ncol
       if (kdthdp(i) /= 0) then
          k = kdthdp(i)
          kp1 = min(k+1,pver)
          ! Note: strat will be zero unless ocnfrac > 0.01
          strat = min(1._r8,max(0._r8, ocnfrac(i) * ((theta(i,k700)-thetas(i))*.057_r8-.5573_r8) ) )
          !
          ! assign the stratus to the layer just below max inversion
          ! the relative humidity changes so rapidly across the inversion
          ! that it is not safe to just look immediately below the inversion
          ! so limit the stratus cloud by rh in both layers below the inversion
          !
          cldst(i,k) = min(strat,max(rh(i,k),rh(i,kp1)))
       end if
    end do
    end if  ! .not.inversion_cld_off

    call cldfrc_total_cloud(ncol, rhcloud, cldst, concld, cloud)

    call outfld( 'SH_CLD  ', shallowcu   , pcols, lchnk )
    call outfld( 'DP_CLD  ', deepcu      , pcols, lchnk )

    !
    return
  end subroutine cldfrc

!================================================================================================

  subroutine cldfrc_layer_rh(ncol, landfrac_local, snowh_local, clat_local, pmid_local, pref_mid_local, q_local, &
       rh_local, rhcloud_local, rhu00_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr, c_double
    use physconst, only: pi

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: landfrac_local(pcols)
    real(r8), target, intent(in) :: snowh_local(pcols)
    real(r8), target, intent(in) :: clat_local(pcols)
    real(r8), target, intent(in) :: pmid_local(pcols,pver)
    real(r8), target, intent(in) :: pref_mid_local(pver)
    real(r8), target, intent(in) :: q_local(pcols,pver)
    real(r8), target, intent(in) :: rh_local(pcols,pver)
    real(r8), target, intent(inout) :: rhcloud_local(pcols,pver)
    real(r8), target, intent(inout) :: rhu00_local(pcols,pver)

    integer(c_int64_t) :: cldfrc_freeze_dry_i

    interface
       subroutine cldfrc_batch_layer_rh_codon(ncol_c, pcols_c, pver_c, top_lev_c, cldfrc_freeze_dry_c, premib_c, &
            premit_c, rhminl_c, rhminl_adj_land_c, rhminh_c, rhminp_c, cldfrc_rhminp_botmb_c, unset_r8_c, pi_c, &
            landfrac_p, snowh_p, clat_p, pmid_p, pref_mid_p, q_p, rh_p, rhcloud_p, rhu00_p) &
            bind(c, name="cldfrc_batch_layer_rh_codon")
         use iso_c_binding, only: c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, cldfrc_freeze_dry_c
         real(c_double), value :: premib_c, premit_c, rhminl_c, rhminl_adj_land_c, rhminh_c, rhminp_c
         real(c_double), value :: cldfrc_rhminp_botmb_c, unset_r8_c, pi_c
         type(c_ptr), value :: landfrac_p, snowh_p, clat_p, pmid_p, pref_mid_p, q_p, rh_p, rhcloud_p, rhu00_p
       end subroutine cldfrc_batch_layer_rh_codon
    end interface

    call cldfrc_batch_select_impl()

    if (use_native_cldfrc_batch_impl) then
       call cldfrc_layer_rh_native(ncol, landfrac_local, snowh_local, clat_local, pmid_local, pref_mid_local, &
            q_local, rh_local, rhcloud_local, rhu00_local)
       return
    end if

    if (cldfrc_freeze_dry) then
       cldfrc_freeze_dry_i = 1_c_int64_t
    else
       cldfrc_freeze_dry_i = 0_c_int64_t
    end if

    call cldfrc_batch_log_entered()

    call cldfrc_batch_dispatch_codon( &
         2_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(top_lev, c_int64_t), cldfrc_freeze_dry_i, real(premib, c_double), real(premit, c_double), &
         real(rhminl, c_double), real(rhminl_adj_land, c_double), real(rhminh, c_double), real(rhminp, c_double), &
         real(cldfrc_rhminp_botmb, c_double), real(unset_r8, c_double), real(pi, c_double), c_loc(landfrac_local), &
         c_loc(snowh_local), c_loc(clat_local), c_loc(pmid_local), c_loc(pref_mid_local), c_loc(q_local), &
         c_loc(rh_local), c_loc(rhcloud_local), c_loc(rhu00_local), c_loc(rhu00_local))

  end subroutine cldfrc_layer_rh

!================================================================================================

  subroutine cldfrc_ice_wilson(ncol, cldice_local, qs_local, rhcloud_local, icecldf_local, liqcldf_local, cloud_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr, c_double

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: cldice_local(pcols,pver)
    real(r8), target, intent(in) :: qs_local(pcols,pver)
    real(r8), target, intent(in) :: rhcloud_local(pcols,pver)
    real(r8), target, intent(inout) :: icecldf_local(pcols,pver)
    real(r8), target, intent(inout) :: liqcldf_local(pcols,pver)
    real(r8), target, intent(inout) :: cloud_local(pcols,pver)

    interface
       subroutine cldfrc_batch_ice_wilson_codon(ncol_c, pcols_c, pver_c, top_lev_c, icecrit_c, one_sixth_c, &
            two_thirds_c, two_pow_three_halves_c, phi_offset_c, cldice_p, qs_p, rhcloud_p, icecldf_p, liqcldf_p, &
            cloud_p) bind(c, name="cldfrc_batch_ice_wilson_codon")
         use iso_c_binding, only: c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         real(c_double), value :: icecrit_c, one_sixth_c, two_thirds_c, two_pow_three_halves_c, phi_offset_c
         type(c_ptr), value :: cldice_p, qs_p, rhcloud_p, icecldf_p, liqcldf_p, cloud_p
       end subroutine cldfrc_batch_ice_wilson_codon
    end interface

    call cldfrc_batch_select_impl()

    if (use_native_cldfrc_batch_impl) then
       call cldfrc_ice_wilson_native(ncol, cldice_local, qs_local, rhcloud_local, icecldf_local, liqcldf_local, &
            cloud_local)
       return
    end if

    call cldfrc_batch_log_entered()

    call cldfrc_batch_dispatch_codon( &
         3_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(top_lev, c_int64_t), 0_c_int64_t, real(icecrit, c_double), real(1._r8/6._r8, c_double), &
         real(2._r8/3._r8, c_double), real(2._r8**(3._r8/2._r8), c_double), real(4._r8*3.1415927_r8, c_double), &
         0._c_double, 0._c_double, 0._c_double, 0._c_double, c_loc(cldice_local), c_loc(qs_local), &
         c_loc(rhcloud_local), c_loc(icecldf_local), c_loc(liqcldf_local), c_loc(cloud_local), c_loc(cloud_local), &
         c_loc(cloud_local), c_loc(cloud_local), c_loc(cloud_local))

  end subroutine cldfrc_ice_wilson

!================================================================================================

  subroutine cldfrc_state_init(ncol, q_local, qs_local, relhum_local, rh_local, cloud_local, icecldf_local, &
       liqcldf_local, rhcloud_local, cldst_local, concld_local, dindex_local, rhpert_local)

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr, c_double

    integer, intent(in) :: ncol
    integer, intent(in) :: dindex_local
    real(r8), intent(in) :: rhpert_local
    real(r8), target, intent(in) :: q_local(pcols,pver)
    real(r8), target, intent(in) :: qs_local(pcols,pver)
    real(r8), target, intent(inout) :: relhum_local(pcols,pver)
    real(r8), target, intent(inout) :: rh_local(pcols,pver)
    real(r8), target, intent(inout) :: cloud_local(pcols,pver)
    real(r8), target, intent(inout) :: icecldf_local(pcols,pver)
    real(r8), target, intent(inout) :: liqcldf_local(pcols,pver)
    real(r8), target, intent(inout) :: rhcloud_local(pcols,pver)
    real(r8), target, intent(inout) :: cldst_local(pcols,pver)
    real(r8), target, intent(inout) :: concld_local(pcols,pver)

    interface
       subroutine cldfrc_batch_state_init_codon(ncol_c, pcols_c, pver_c, top_lev_c, dindex_c, rhpert_c, q_p, qs_p, &
            relhum_p, rh_p, cloud_p, icecldf_p, liqcldf_p, rhcloud_p, cldst_p, concld_p) &
            bind(c, name="cldfrc_batch_state_init_codon")
         use iso_c_binding, only: c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, dindex_c
         real(c_double), value :: rhpert_c
         type(c_ptr), value :: q_p, qs_p, relhum_p, rh_p, cloud_p, icecldf_p, liqcldf_p, rhcloud_p, cldst_p, concld_p
       end subroutine cldfrc_batch_state_init_codon
    end interface

    call cldfrc_batch_select_impl()

    if (use_native_cldfrc_batch_impl) then
       call cldfrc_state_init_native(ncol, q_local, qs_local, relhum_local, rh_local, cloud_local, icecldf_local, &
            liqcldf_local, rhcloud_local, cldst_local, concld_local, dindex_local, rhpert_local)
       return
    end if

    call cldfrc_batch_log_entered()

    call cldfrc_batch_dispatch_codon( &
         1_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(top_lev, c_int64_t), int(dindex_local, c_int64_t), real(rhpert_local, c_double), &
         0._c_double, 0._c_double, 0._c_double, 0._c_double, 0._c_double, 0._c_double, 0._c_double, 0._c_double, &
         c_loc(q_local), c_loc(qs_local), c_loc(relhum_local), c_loc(rh_local), c_loc(cloud_local), &
         c_loc(icecldf_local), c_loc(liqcldf_local), c_loc(rhcloud_local), c_loc(cldst_local), c_loc(concld_local))

  end subroutine cldfrc_state_init

!================================================================================================

  subroutine cldfrc_state_init_native(ncol, q_local, qs_local, relhum_local, rh_local, cloud_local, icecldf_local, &
       liqcldf_local, rhcloud_local, cldst_local, concld_local, dindex_local, rhpert_local)

    integer, intent(in) :: ncol
    integer, intent(in) :: dindex_local
    real(r8), intent(in) :: rhpert_local
    real(r8), intent(in) :: q_local(pcols,pver)
    real(r8), intent(in) :: qs_local(pcols,pver)
    real(r8), intent(inout) :: relhum_local(pcols,pver)
    real(r8), intent(inout) :: rh_local(pcols,pver)
    real(r8), intent(inout) :: cloud_local(pcols,pver)
    real(r8), intent(inout) :: icecldf_local(pcols,pver)
    real(r8), intent(inout) :: liqcldf_local(pcols,pver)
    real(r8), intent(inout) :: rhcloud_local(pcols,pver)
    real(r8), intent(inout) :: cldst_local(pcols,pver)
    real(r8), intent(inout) :: concld_local(pcols,pver)

    integer :: i, k

    cloud_local   = 0._r8
    icecldf_local = 0._r8
    liqcldf_local = 0._r8
    rhcloud_local = 0._r8
    cldst_local   = 0._r8
    concld_local  = 0._r8

    do k=top_lev,pver
       do i=1,ncol
          rh_local(i,k)     = q_local(i,k)/qs_local(i,k)*(1.0_r8+real(dindex_local,r8)*rhpert_local)
          relhum_local(i,k) = rh_local(i,k)
       end do
    end do

  end subroutine cldfrc_state_init_native

!================================================================================================

  subroutine cldfrc_layer_rh_native(ncol, landfrac_local, snowh_local, clat_local, pmid_local, pref_mid_local, &
       q_local, rh_local, rhcloud_local, rhu00_local)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: landfrac_local(pcols)
    real(r8), intent(in) :: snowh_local(pcols)
    real(r8), intent(in) :: clat_local(pcols)
    real(r8), intent(in) :: pmid_local(pcols,pver)
    real(r8), intent(in) :: pref_mid_local(pver)
    real(r8), intent(in) :: q_local(pcols,pver)
    real(r8), intent(in) :: rh_local(pcols,pver)
    real(r8), intent(inout) :: rhcloud_local(pcols,pver)
    real(r8), intent(inout) :: rhu00_local(pcols,pver)

    integer :: i, k
    real(r8) :: rhdif, rhlim, rhwght

    do k=top_lev+1,pver
       do i=1,ncol
          if ( pmid_local(i,k).ge.premib ) then
             if (nint(landfrac_local(i)) == 1 .and. (snowh_local(i) <= 0.000001_r8)) then
                rhlim = rhminl - rhminl_adj_land
             else
                rhlim = rhminl
             endif

             rhdif = (rh_local(i,k) - rhlim)/(1.0_r8-rhlim)
             rhcloud_local(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)

             if (cldfrc_freeze_dry) then
                rhcloud_local(i,k) = rhcloud_local(i,k)*max(0.15_r8,min(1.0_r8,q_local(i,k)/0.0030_r8))
             endif

          else if ( pmid_local(i,k).lt.premit ) then
             rhlim = relhum_min(pref_mid_local(k),clat_local(i))
             rhdif = (rh_local(i,k) - rhlim)/(1.0_r8-rhlim)
             rhcloud_local(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          else
             rhwght = (premib-(max(pmid_local(i,k),premit)))/(premib-premit)

             if (nint(landfrac_local(i)) == 1 .and. (snowh_local(i) <= 0.000001_r8)) then
                rhlim = relhum_min(pref_mid_local(k),clat_local(i))*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
             else
                rhlim = relhum_min(pref_mid_local(k),clat_local(i))*rhwght + rhminl*(1.0_r8-rhwght)
             endif
             rhdif = (rh_local(i,k) - rhlim)/(1.0_r8-rhlim)
             rhcloud_local(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          end if

          rhu00_local(i,k)=rhlim
       end do
    end do

  end subroutine cldfrc_layer_rh_native

!================================================================================================

  subroutine cldfrc_ice_wilson_native(ncol, cldice_local, qs_local, rhcloud_local, icecldf_local, liqcldf_local, &
       cloud_local)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: cldice_local(pcols,pver)
    real(r8), intent(in) :: qs_local(pcols,pver)
    real(r8), intent(in) :: rhcloud_local(pcols,pver)
    real(r8), intent(inout) :: icecldf_local(pcols,pver)
    real(r8), intent(inout) :: liqcldf_local(pcols,pver)
    real(r8), intent(inout) :: cloud_local(pcols,pver)

    integer :: i, k
    real(r8) :: ncf_local, phi_local

    do k=top_lev+1,pver
       do i=1,ncol
          ncf_local = cldice_local(i,k)/((1._r8 - icecrit)*qs_local(i,k))
          if (ncf_local.le.0._r8) then
             icecldf_local(i,k)=0._r8
          else if (ncf_local.gt.0._r8 .and. ncf_local.le.1._r8/6._r8) then
             icecldf_local(i,k)=0.5_r8*(6._r8 * ncf_local)**(2._r8/3._r8)
          else if (ncf_local.gt.1._r8/6._r8 .and. ncf_local.lt.1._r8) then
             phi_local=(acos(3._r8*(1._r8-ncf_local)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
             icecldf_local(i,k)=(1._r8 - 4._r8 * cos(phi_local) * cos(phi_local))
          else
             icecldf_local(i,k)=1._r8
          endif
          icecldf_local(i,k)=max(0._r8,min(icecldf_local(i,k),1._r8))

          liqcldf_local(i,k)=(1._r8 - icecldf_local(i,k))* rhcloud_local(i,k)
          cloud_local(i,k)=liqcldf_local(i,k) + icecldf_local(i,k)
       end do
    end do

  end subroutine cldfrc_ice_wilson_native

!================================================================================================

  subroutine cldfrc_total_cloud(ncol, rhcloud_local, cldst_local, concld_local, cloud_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: rhcloud_local(pcols,pver)
    real(r8), target, intent(in) :: cldst_local(pcols,pver)
    real(r8), target, intent(in) :: concld_local(pcols,pver)
    real(r8), target, intent(inout) :: cloud_local(pcols,pver)

    interface
       subroutine cldfrc_batch_total_cloud_codon(ncol_c, pcols_c, pver_c, top_lev_c, rhcloud_p, cldst_p, concld_p, &
            cloud_p) bind(c, name="cldfrc_batch_total_cloud_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         type(c_ptr), value :: rhcloud_p, cldst_p, concld_p, cloud_p
       end subroutine cldfrc_batch_total_cloud_codon
    end interface

    call cldfrc_batch_select_impl()

    if (use_native_cldfrc_batch_impl) then
       call cldfrc_total_cloud_native(ncol, rhcloud_local, cldst_local, concld_local, cloud_local)
       return
    end if

    call cldfrc_batch_log_entered()

    call cldfrc_batch_dispatch_codon( &
         5_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
         0_c_int64_t, 0._c_double, 0._c_double, 0._c_double, 0._c_double, 0._c_double, 0._c_double, &
         0._c_double, 0._c_double, 0._c_double, c_loc(rhcloud_local), c_loc(cldst_local), c_loc(concld_local), &
         c_loc(cloud_local), c_loc(cloud_local), c_loc(cloud_local), c_loc(cloud_local), c_loc(cloud_local), &
         c_loc(cloud_local), c_loc(cloud_local))

  end subroutine cldfrc_total_cloud

!================================================================================================

  subroutine cldfrc_total_cloud_native(ncol, rhcloud_local, cldst_local, concld_local, cloud_local)

    integer, intent(in) :: ncol
    real(r8), intent(in) :: rhcloud_local(pcols,pver)
    real(r8), intent(in) :: cldst_local(pcols,pver)
    real(r8), intent(in) :: concld_local(pcols,pver)
    real(r8), intent(inout) :: cloud_local(pcols,pver)

    integer :: i, k

    do k=top_lev,pver
       do i=1,ncol
          cloud_local(i,k) = max(rhcloud_local(i,k),cldst_local(i,k))
          cloud_local(i,k) = min(cloud_local(i,k)+concld_local(i,k), 1.0_r8)
       end do
    end do

  end subroutine cldfrc_total_cloud_native

!================================================================================================

  subroutine cldfrc_convective_cover(ncol, use_shfrc_local, shfrc_local, cmfmc_local, cmfmc2_local, &
       shallowcu_local, deepcu_local, concld_local, rh_local)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    logical, intent(in) :: use_shfrc_local
    real(r8), target, intent(in) :: shfrc_local(pcols,pver)
    real(r8), target, intent(in) :: cmfmc_local(pcols,pverp)
    real(r8), target, intent(in) :: cmfmc2_local(pcols,pverp)
    real(r8), target, intent(inout) :: shallowcu_local(pcols,pver)
    real(r8), target, intent(inout) :: deepcu_local(pcols,pver)
    real(r8), target, intent(inout) :: concld_local(pcols,pver)
    real(r8), target, intent(inout) :: rh_local(pcols,pver)

    integer(c_int64_t) :: use_shfrc_c

    interface
       subroutine cldfrc_batch_convective_cover_codon(ncol_c, pcols_c, pver_c, top_lev_c, use_shfrc_c, &
            sh1_c, sh2_c, dp1_c, dp2_c, shfrc_p, cmfmc_p, cmfmc2_p, shallowcu_p, deepcu_p, concld_p, rh_p) &
            bind(c, name="cldfrc_batch_convective_cover_codon")
         use iso_c_binding, only: c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, use_shfrc_c
         real(c_double), value :: sh1_c, sh2_c, dp1_c, dp2_c
         type(c_ptr), value :: shfrc_p, cmfmc_p, cmfmc2_p, shallowcu_p, deepcu_p, concld_p, rh_p
       end subroutine cldfrc_batch_convective_cover_codon
    end interface

    call cldfrc_batch_select_impl()

    if (use_native_cldfrc_batch_impl) then
       call cldfrc_convective_cover_native(ncol, use_shfrc_local, shfrc_local, cmfmc_local, cmfmc2_local, &
            shallowcu_local, deepcu_local, concld_local, rh_local)
       return
    end if

    use_shfrc_c = 0_c_int64_t
    if (use_shfrc_local) use_shfrc_c = 1_c_int64_t

    call cldfrc_batch_log_entered()

    call cldfrc_batch_dispatch_codon( &
         4_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
         use_shfrc_c, real(sh1, c_double), real(sh2, c_double), real(dp1, c_double), real(dp2, c_double), &
         0._c_double, 0._c_double, 0._c_double, 0._c_double, 0._c_double, c_loc(shfrc_local), c_loc(cmfmc_local), &
         c_loc(cmfmc2_local), c_loc(shallowcu_local), c_loc(deepcu_local), c_loc(concld_local), c_loc(rh_local), &
         c_loc(rh_local), c_loc(rh_local), c_loc(rh_local))

  end subroutine cldfrc_convective_cover

!================================================================================================

  subroutine cldfrc_convective_cover_native(ncol, use_shfrc_local, shfrc_local, cmfmc_local, cmfmc2_local, &
       shallowcu_local, deepcu_local, concld_local, rh_local)

    integer, intent(in) :: ncol
    logical, intent(in) :: use_shfrc_local
    real(r8), intent(in) :: shfrc_local(pcols,pver)
    real(r8), intent(in) :: cmfmc_local(pcols,pverp)
    real(r8), intent(in) :: cmfmc2_local(pcols,pverp)
    real(r8), intent(inout) :: shallowcu_local(pcols,pver)
    real(r8), intent(inout) :: deepcu_local(pcols,pver)
    real(r8), intent(inout) :: concld_local(pcols,pver)
    real(r8), intent(inout) :: rh_local(pcols,pver)

    integer :: i, k

    do k=top_lev,pver
       do i=1,ncol
          if ( .not. use_shfrc_local ) then
             shallowcu_local(i,k) = max(0.0_r8,min(sh1*log(1.0_r8+sh2*cmfmc2_local(i,k+1)),0.30_r8))
          else
             shallowcu_local(i,k) = shfrc_local(i,k)
          endif
          deepcu_local(i,k) = max(0.0_r8,min(dp1*log(1.0_r8+dp2*(cmfmc_local(i,k+1)-cmfmc2_local(i,k+1))),0.60_r8))
          concld_local(i,k) = min(shallowcu_local(i,k) + deepcu_local(i,k),0.80_r8)
          rh_local(i,k) = (rh_local(i,k) - concld_local(i,k))/(1.0_r8 - concld_local(i,k))
       end do
    end do

  end subroutine cldfrc_convective_cover_native

!================================================================================================

  subroutine cldfrc_fice(ncol, t, fice, fsnow)
!
! Compute the fraction of the total cloud water which is in ice phase.
! The fraction depends on temperature only. 
! This is the form that was used for radiation, the code came from cldefr originally
! 
! Author: B. A. Boville Sept 10, 2002
!  modified: PJR 3/13/03 (added fsnow to ascribe snow production for convection )
!-----------------------------------------------------------------------
    use physconst, only: tmelt
    use iso_c_binding, only: c_int64_t, c_loc, c_ptr, c_double

! Arguments
    integer,  intent(in)  :: ncol                 ! number of active columns
    real(r8), target, intent(in)  :: t(pcols,pver)        ! temperature

    real(r8), target, intent(out) :: fice(pcols,pver)     ! Fractional ice content within cloud
    real(r8), target, intent(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection

! Local variables
    real(r8) :: tmax_fice                         ! max temperature for cloud ice formation
    real(r8) :: tmin_fice                         ! min temperature for cloud ice formation
    real(r8) :: tmax_fsnow                        ! max temperature for transition to convective snow
    real(r8) :: tmin_fsnow                        ! min temperature for transition to convective snow

    interface
       subroutine cldfrc_batch_fice_codon(ncol_c, pcols_c, pver_c, top_lev_c, t_p, fice_p, fsnow_p, &
            tmax_fice_c, tmin_fice_c, tmax_fsnow_c, tmin_fsnow_c) bind(c, name="cldfrc_batch_fice_codon")
         use iso_c_binding, only: c_int64_t, c_ptr, c_double
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c
         type(c_ptr), value :: t_p, fice_p, fsnow_p
         real(c_double), value :: tmax_fice_c, tmin_fice_c, tmax_fsnow_c, tmin_fsnow_c
       end subroutine cldfrc_batch_fice_codon
    end interface

!-----------------------------------------------------------------------

    call cldfrc_batch_select_impl()

    tmax_fice = tmelt - 10._r8        ! max temperature for cloud ice formation
    tmin_fice = tmax_fice - 30._r8    ! min temperature for cloud ice formation
    tmax_fsnow = tmelt                ! max temperature for transition to convective snow
    tmin_fsnow = tmelt - 5._r8        ! min temperature for transition to convective snow

    if (use_native_cldfrc_batch_impl) then
       call cldfrc_fice_native(ncol, t, fice, fsnow)
       return
    end if

    call cldfrc_batch_log_entered()

    call cldfrc_batch_dispatch_codon( &
         6_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev, c_int64_t), &
         0_c_int64_t, real(tmax_fice, c_double), real(tmin_fice, c_double), real(tmax_fsnow, c_double), &
         real(tmin_fsnow, c_double), 0._c_double, 0._c_double, 0._c_double, 0._c_double, 0._c_double, &
         c_loc(t), c_loc(fice), c_loc(fsnow), c_loc(t), c_loc(t), c_loc(t), c_loc(t), c_loc(t), c_loc(t), c_loc(t))

  end subroutine cldfrc_fice

!================================================================================================

  subroutine cldfrc_fice_native(ncol, t, fice, fsnow)
!
! Compute the fraction of the total cloud water which is in ice phase.
! The fraction depends on temperature only. 
! This is the form that was used for radiation, the code came from cldefr originally
! 
! Author: B. A. Boville Sept 10, 2002
!  modified: PJR 3/13/03 (added fsnow to ascribe snow production for convection )
!-----------------------------------------------------------------------
    use physconst, only: tmelt

! Arguments
    integer,  intent(in)  :: ncol                 ! number of active columns
    real(r8), intent(in)  :: t(pcols,pver)        ! temperature

    real(r8), intent(out) :: fice(pcols,pver)     ! Fractional ice content within cloud
    real(r8), intent(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection

! Local variables
    real(r8) :: tmax_fice                         ! max temperature for cloud ice formation
    real(r8) :: tmin_fice                         ! min temperature for cloud ice formation
    real(r8) :: tmax_fsnow                        ! max temperature for transition to convective snow
    real(r8) :: tmin_fsnow                        ! min temperature for transition to convective snow

    integer :: i,k                                ! loop indexes

!-----------------------------------------------------------------------

    tmax_fice = tmelt - 10._r8        ! max temperature for cloud ice formation
    tmin_fice = tmax_fice - 30._r8    ! min temperature for cloud ice formation
    tmax_fsnow = tmelt                ! max temperature for transition to convective snow
    tmin_fsnow = tmelt - 5._r8        ! min temperature for transition to convective snow

    fice(:,:top_lev-1) = 0._r8
    fsnow(:,:top_lev-1) = 0._r8

! Define fractional amount of cloud that is ice
    do k=top_lev,pver
       do i=1,ncol

! If warmer than tmax then water phase
          if (t(i,k) > tmax_fice) then
             fice(i,k) = 0.0_r8

! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fice) then
             fice(i,k) = 1.0_r8

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fice(i,k) =(tmax_fice - t(i,k)) / (tmax_fice - tmin_fice)
          end if

! snow fraction partitioning

! If warmer than tmax then water phase
          if (t(i,k) > tmax_fsnow) then
             fsnow(i,k) = 0.0_r8

! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fsnow) then
             fsnow(i,k) = 1.0_r8

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fsnow(i,k) =(tmax_fsnow - t(i,k)) / (tmax_fsnow - tmin_fsnow)
          end if

       end do
    end do

  end subroutine cldfrc_fice_native

  !-----------------------------------------------------------------------------
  ! Sets rhmin to a different value (rhminp) poleward of +/- 60 deg latitude and 
  ! pressure levels less than cldfrc_rhminp_botmb (hPa) if cldfrc_rhminp is specified
  ! ** This is used only for special waccm/cam-chem cases with cam4 physics **
  !-----------------------------------------------------------------------------
  function relhum_min(press,lat) result(rh)
    use physconst, only: pi

    real(r8), intent(in) :: press, lat
    real(r8) :: rh

    rh = rhminh
    if (rhminp .eq. unset_r8 ) return

    if ((press .lt. cldfrc_rhminp_botmb*1.e2_r8) .and. &
        ( abs( lat*180._r8/pi ) .gt. 60._r8 ) ) then
       rh = rhminp
    endif

  end function relhum_min

end module cloud_fraction
