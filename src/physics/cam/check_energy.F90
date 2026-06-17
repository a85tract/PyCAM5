
module check_energy

!---------------------------------------------------------------------------------
! Purpose:
!
! Module to check 
!   1. vertically integrated total energy and water conservation for each
!      column within the physical parameterizations
!
!   2. global mean total energy conservation between the physics output state
!      and the input state on the next time step.
!
!   3. add a globally uniform heating term to account for any change of total energy in 2.
!
! Author: Byron Boville  Oct 31, 2002
!         
! Modifications:
!   03.03.29  Boville  Add global energy check and fixer.        
!
!---------------------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, begchunk, endchunk
  use spmd_utils,      only: masterproc
  
  use phys_gmean,      only: gmean
  use physconst,       only: gravit, latvap, latice
  use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init
  use constituents,    only: cnst_get_ind, pcnst, cnst_name, cnst_get_type_byind
  use time_manager,    only: is_first_step
  use cam_logfile,     only: iulog
  use cam_abortutils,  only: endrun
  use iso_c_binding,   only: c_int64_t

  implicit none
  private

! Public types:
  public check_tracers_data

! Public methods
  public :: check_energy_defaultopts ! set default namelist values
  public :: check_energy_setopts   ! set namelist values
  public :: check_energy_register  ! register fields in physics buffer
  public :: check_energy_get_integrals ! get energy integrals computed in check_energy_gmean
  public :: check_energy_init      ! initialization of module
  public :: check_energy_timestep_init  ! timestep initialization of energy integrals and cumulative boundary fluxes
  public :: check_energy_chng      ! check changes in integrals against cumulative boundary fluxes
  public :: check_energy_gmean     ! global means of physics input and output total energy
  public :: check_energy_fix       ! add global mean energy difference as a heating
  public :: check_energy_fix_native
  public :: check_tracers_init      ! initialize tracer integrals and cumulative boundary fluxes
  public :: check_tracers_chng      ! check changes in integrals against cumulative boundary fluxes


! Private module data

  logical  :: print_energy_errors = .false.
  logical  :: use_native_fix_impl = .false.
  logical  :: fix_impl_selected = .false.
  logical  :: use_native_timestep_init_impl = .false.
  logical  :: timestep_init_impl_selected = .false.
  logical  :: use_native_chng_impl = .false.
  logical  :: chng_impl_selected = .false.
  logical  :: use_native_gmean_impl = .false.
  logical  :: gmean_impl_selected = .false.
  logical  :: use_native_energy_batch_impl = .false.
  logical  :: energy_batch_impl_selected = .false.
  logical  :: energy_batch_entered_logged = .false.
  logical  :: use_native_tracers_init_impl = .false.
  logical  :: tracers_init_impl_selected = .false.
  logical  :: use_native_tracers_chng_impl = .false.
  logical  :: tracers_chng_impl_selected = .false.
  logical  :: use_native_tracers_batch_impl = .false.
  logical  :: tracers_batch_impl_selected = .false.
  logical  :: tracers_batch_entered_logged = .false.
  logical  :: check_energy_defaultopts_logged = .false.
  logical  :: check_energy_setopts_logged = .false.
  logical  :: check_energy_init_logged = .false.
  logical  :: check_energy_timestep_init_logged = .false.
  logical  :: check_energy_chng_logged = .false.

  real(r8) :: teout_glob           ! global mean energy of output state
  real(r8) :: teinp_glob           ! global mean energy of input state
  real(r8) :: tedif_glob           ! global mean energy difference
  real(r8) :: psurf_glob           ! global mean surface pressure
  real(r8) :: ptopb_glob           ! global mean top boundary pressure
  real(r8) :: heat_glob            ! global mean heating rate

! Physics buffer indices
  
  integer  :: teout_idx  = 0       ! teout index in physics buffer 
  integer  :: dtcore_idx = 0       ! dtcore index in physics buffer 

  type check_tracers_data
     real(r8) :: tracer(pcols,pcnst)       ! initial vertically integrated total (kinetic + static) energy
     real(r8) :: tracer_tnd(pcols,pcnst)   ! cumulative boundary flux of total energy
     integer :: count(pcnst)               ! count of values with significant imbalances
  end type check_tracers_data

  interface
     subroutine check_energy_batch_dispatch_codon(stage_c, ncol_c, pver_c, pcols_c, psetcols_c, pcnst_c, &
          ixcldliq_c, ixcldice_c, ixrain_c, ixsnow_c, scalar1_c, scalar2_c, scalar3_c, &
          p1_p, p2_p, p3_p, p4_p, p5_p, p6_p, p7_p, p8_p, p9_p, p10_p, p11_p, p12_p, p13_p, p14_p, &
          p15_p, p16_p, p17_p, p18_p) bind(c, name="check_energy_batch_dispatch_codon")
        use iso_c_binding, only: c_double, c_int64_t, c_ptr
        integer(c_int64_t), value :: stage_c, ncol_c, pver_c, pcols_c, psetcols_c, pcnst_c
        integer(c_int64_t), value :: ixcldliq_c, ixcldice_c, ixrain_c, ixsnow_c
        real(c_double), value :: scalar1_c, scalar2_c, scalar3_c
        type(c_ptr), value :: p1_p, p2_p, p3_p, p4_p, p5_p, p6_p, p7_p, p8_p, p9_p, p10_p
        type(c_ptr), value :: p11_p, p12_p, p13_p, p14_p, p15_p, p16_p, p17_p, p18_p
     end subroutine check_energy_batch_dispatch_codon
     subroutine check_energy_timestep_init_codon(ncol_c, pver_c, psetcols_c, pcnst_c, latvap_c, latice_c, gravit_c, &
          ixcldliq_c, ixcldice_c, ixrain_c, ixsnow_c, state_u_p, state_v_p, state_s_p, state_q_p, state_pdel_p, &
          ke_p, se_p, wv_p, wl_p, wi_p, state_te_ini_p, state_tw_ini_p) bind(c, name="check_energy_timestep_init_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pver_c, psetcols_c, pcnst_c
       integer(c_int64_t), value :: ixcldliq_c, ixcldice_c, ixrain_c, ixsnow_c
       real(c_double), value :: latvap_c, latice_c, gravit_c
       type(c_ptr), value :: state_u_p, state_v_p, state_s_p, state_q_p, state_pdel_p
       type(c_ptr), value :: ke_p, se_p, wv_p, wl_p, wi_p, state_te_ini_p, state_tw_ini_p
     end subroutine check_energy_timestep_init_codon
     subroutine check_energy_chng_codon(ncol_c, pver_c, psetcols_c, latvap_c, latice_c, gravit_c, &
          ixcldliq_c, ixcldice_c, ixrain_c, ixsnow_c, state_u_p, state_v_p, state_s_p, state_q_p, state_pdel_p, &
          flx_vap_p, flx_cnd_p, flx_ice_p, flx_sen_p, ke_p, se_p, wv_p, wl_p, wi_p, tend_te_tnd_p, &
          tend_tw_tnd_p, state_te_cur_p, state_tw_cur_p) bind(c, name="check_energy_chng_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pver_c, psetcols_c
       integer(c_int64_t), value :: ixcldliq_c, ixcldice_c, ixrain_c, ixsnow_c
       real(c_double), value :: latvap_c, latice_c, gravit_c
       type(c_ptr), value :: state_u_p, state_v_p, state_s_p, state_q_p, state_pdel_p
       type(c_ptr), value :: flx_vap_p, flx_cnd_p, flx_ice_p, flx_sen_p
       type(c_ptr), value :: ke_p, se_p, wv_p, wl_p, wi_p, tend_te_tnd_p, tend_tw_tnd_p
       type(c_ptr), value :: state_te_cur_p, state_tw_cur_p
     end subroutine check_energy_chng_codon
     subroutine check_energy_gmean_codon(ncol_c, state_te_ini_p, teout_p, pint_surf_p, te1_p, te2_p, te3_p) &
          bind(c, name="check_energy_gmean_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c
       type(c_ptr), value :: state_te_ini_p, teout_p, pint_surf_p, te1_p, te2_p, te3_p
     end subroutine check_energy_gmean_codon
     subroutine check_energy_fix_codon(ncol_c, pcols_c, pver_c, psetcols_c, heat_glob_c, gravit_c, &
          state_pint_p, ptend_s_p, eshflx_p) bind(c, name="check_energy_fix_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, psetcols_c
       real(c_double), value :: heat_glob_c, gravit_c
       type(c_ptr), value :: state_pint_p, ptend_s_p, eshflx_p
     end subroutine check_energy_fix_codon
     subroutine check_tracers_init_codon() bind(c, name="check_tracers_init_codon")
     end subroutine check_tracers_init_codon
     subroutine check_tracers_chng_codon() bind(c, name="check_tracers_chng_codon")
     end subroutine check_tracers_chng_codon
     function check_energy_defaultopts_codon(flag_c) result(out_c) bind(c, name="check_energy_defaultopts_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: out_c
     end function check_energy_defaultopts_codon
     function check_energy_setopts_codon(flag_c) result(out_c) bind(c, name="check_energy_setopts_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: out_c
     end function check_energy_setopts_codon
     function check_energy_register_codon(flag_c) result(out_c) bind(c, name="check_energy_register_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: out_c
     end function check_energy_register_codon
     function check_energy_init_codon(flag_c) result(out_c) bind(c, name="check_energy_init_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: out_c
     end function check_energy_init_codon
  end interface


!===============================================================================
contains
!===============================================================================

subroutine check_energy_fix_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (fix_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CHECK_ENERGY_FIX_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_fix_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_fix_impl = .false.
   end if

   fix_impl_selected = .true.

   if (masterproc) then
      if (use_native_fix_impl) then
         write(iulog,*) 'check_energy_fix implementation = native'
      else
         write(iulog,*) 'check_energy_fix implementation = codon'
      end if
   end if

end subroutine check_energy_fix_select_impl

subroutine check_energy_timestep_init_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (timestep_init_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CHECK_ENERGY_TIMESTEP_INIT_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_timestep_init_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_timestep_init_impl = .false.
   end if

   timestep_init_impl_selected = .true.

   if (masterproc) then
      if (use_native_timestep_init_impl) then
         write(iulog,*) 'check_energy_timestep_init implementation = native'
      else
         write(iulog,*) 'check_energy_timestep_init implementation = codon'
      end if
   end if

end subroutine check_energy_timestep_init_select_impl

subroutine check_energy_chng_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (chng_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CHECK_ENERGY_CHNG_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_chng_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_chng_impl = .false.
   end if

   chng_impl_selected = .true.

   if (masterproc) then
      if (use_native_chng_impl) then
         write(iulog,*) 'check_energy_chng implementation = native'
      else
         write(iulog,*) 'check_energy_chng implementation = codon'
      end if
   end if

end subroutine check_energy_chng_select_impl

subroutine check_energy_gmean_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (gmean_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CHECK_ENERGY_GMEAN_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_gmean_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_gmean_impl = .false.
   end if

   gmean_impl_selected = .true.

   if (masterproc) then
      if (use_native_gmean_impl) then
         write(iulog,*) 'check_energy_gmean implementation = native'
      else
         write(iulog,*) 'check_energy_gmean implementation = codon'
      end if
   end if

end subroutine check_energy_gmean_select_impl

subroutine check_energy_batch_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line

   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('CHECK_ENERGY_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine check_energy_batch_append_proof

subroutine check_energy_batch_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (energy_batch_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CHECK_ENERGY_BATCH_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_energy_batch_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_energy_batch_impl = .false.
   end if

   energy_batch_impl_selected = .true.

   if (masterproc) then
      if (use_native_energy_batch_impl) then
         write(iulog,*) 'check_energy_batch implementation = native'
         call check_energy_batch_append_proof('check_energy_batch selector entered implementation = native')
      else
         write(iulog,*) 'check_energy_batch implementation = codon'
         call check_energy_batch_append_proof('check_energy_batch selector entered implementation = codon')
      end if
   end if

end subroutine check_energy_batch_select_impl

subroutine check_energy_batch_log_entered()

   if (energy_batch_entered_logged) return
   energy_batch_entered_logged = .true.

   if (masterproc) then
      write(iulog,*) 'check_energy_batch entered (unified stage-dispatch timestep_init/chng/gmean/fix direct = codon)'
      call check_energy_batch_append_proof( &
           'check_energy_batch entered (unified stage-dispatch timestep_init/chng/gmean/fix direct = codon)')
   end if

end subroutine check_energy_batch_log_entered

subroutine check_tracers_init_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (tracers_init_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CHECK_TRACERS_INIT_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_tracers_init_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_tracers_init_impl = .false.
   end if

   tracers_init_impl_selected = .true.

   if (masterproc) then
      if (use_native_tracers_init_impl) then
         write(iulog,*) 'check_tracers_init implementation = native'
      else
         write(iulog,*) 'check_tracers_init implementation = codon'
      end if
   end if

end subroutine check_tracers_init_select_impl

subroutine check_tracers_chng_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (tracers_chng_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CHECK_TRACERS_CHNG_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_tracers_chng_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_tracers_chng_impl = .false.
   end if

   tracers_chng_impl_selected = .true.

   if (masterproc) then
      if (use_native_tracers_chng_impl) then
         write(iulog,*) 'check_tracers_chng implementation = native'
      else
         write(iulog,*) 'check_tracers_chng implementation = codon'
      end if
   end if

end subroutine check_tracers_chng_select_impl

subroutine check_tracers_batch_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line

   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('CHECK_TRACERS_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine check_tracers_batch_append_proof

subroutine check_tracers_batch_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (tracers_batch_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CHECK_TRACERS_BATCH_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_tracers_batch_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_tracers_batch_impl = .false.
   end if

   tracers_batch_impl_selected = .true.

   if (masterproc) then
      if (use_native_tracers_batch_impl) then
         write(iulog,*) 'check_tracers_batch implementation = native'
         call check_tracers_batch_append_proof('check_tracers_batch selector entered implementation = native')
      else
         write(iulog,*) 'check_tracers_batch implementation = codon'
         call check_tracers_batch_append_proof('check_tracers_batch selector entered implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine check_tracers_batch_select_impl

subroutine check_tracers_batch_log_entered()

   if (tracers_batch_entered_logged) return
   tracers_batch_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'check_tracers_batch entered (unified stage-dispatch init/chng direct = codon)'
      call check_tracers_batch_append_proof('check_tracers_batch entered (unified stage-dispatch init/chng direct = codon)')
      call flush(iulog)
   end if

end subroutine check_tracers_batch_log_entered

subroutine check_energy_log_direct(logged, proof_line)

   logical, intent(inout) :: logged
   character(len=*), intent(in) :: proof_line

   if (logged) return
   logged = .true.

   if (masterproc) then
      write(iulog,'(A)') trim(proof_line)
      call flush(iulog)
   end if

end subroutine check_energy_log_direct

subroutine check_energy_defaultopts( &
   print_energy_errors_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   logical,          intent(out), optional :: print_energy_errors_out
   integer(c_int64_t) :: out_c
!-----------------------------------------------------------------------

   if ( present(print_energy_errors_out) ) then
      out_c = check_energy_defaultopts_codon(merge(1_c_int64_t, 0_c_int64_t, print_energy_errors))
      call check_energy_log_direct(check_energy_defaultopts_logged, 'check_energy_defaultopts direct = codon')
      print_energy_errors_out = out_c /= 0_c_int64_t
   endif

end subroutine check_energy_defaultopts

!================================================================================================

subroutine check_energy_setopts( &
   print_energy_errors_in)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   logical,          intent(in), optional :: print_energy_errors_in
   integer(c_int64_t) :: out_c
!-----------------------------------------------------------------------

   if ( present(print_energy_errors_in) ) then
      out_c = check_energy_setopts_codon(merge(1_c_int64_t, 0_c_int64_t, print_energy_errors_in))
      call check_energy_log_direct(check_energy_setopts_logged, 'check_energy_setopts direct = codon')
      print_energy_errors = out_c /= 0_c_int64_t
   endif

end subroutine check_energy_setopts

!================================================================================================

  subroutine check_energy_register()
!
! Register fields in the physics buffer.
! 
!-----------------------------------------------------------------------
    
    use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls
    use physics_buffer, only : pbuf_register_subcol
    use subcol_utils,   only : is_subcol_on

!-----------------------------------------------------------------------
    if (check_energy_register_codon(1_c_int64_t) == 0_c_int64_t) return

! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add_field('TEOUT', 'global',dtype_r8 , (/pcols,dyn_time_lvls/),      teout_idx)
    call pbuf_add_field('DTCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dtcore_idx)
    if(is_subcol_on()) then
      call pbuf_register_subcol('TEOUT', 'phys_register', teout_idx)
      call pbuf_register_subcol('DTCORE', 'phys_register', dtcore_idx)
    end if

  end subroutine check_energy_register

!===============================================================================

subroutine check_energy_get_integrals( tedif_glob_out, heat_glob_out )

!----------------------------------------------------------------------- 
! Purpose: Return energy integrals
!-----------------------------------------------------------------------

     real(r8), intent(out), optional :: tedif_glob_out
     real(r8), intent(out), optional :: heat_glob_out

!-----------------------------------------------------------------------

   if ( present(tedif_glob_out) ) then
      tedif_glob_out = tedif_glob
   endif
   if ( present(heat_glob_out) ) then
      heat_glob_out = heat_glob
   endif

end subroutine check_energy_get_integrals

!================================================================================================

  subroutine check_energy_init()
!
! Initialize the energy conservation module
! 
!-----------------------------------------------------------------------
    use cam_history,       only: addfld, add_default, phys_decomp
    use phys_control,      only: phys_getopts

    implicit none

    logical          :: history_budget
    integer          :: history_budget_histfile_num ! output history file number for budget fields

!-----------------------------------------------------------------------
    if (check_energy_init_codon(1_c_int64_t) == 0_c_int64_t) return
    call check_energy_log_direct(check_energy_init_logged, &
         'check_energy_init direct = codon; phys_getopts/history native CAM API islands')

    call phys_getopts( history_budget_out = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num)

! register history variables
    call addfld('TEINP   ', 'W/m2', 1,    'A', 'Total energy of physics input',    phys_decomp)
    call addfld('TEOUT   ', 'W/m2', 1,    'A', 'Total energy of physics output',   phys_decomp)
    call addfld('TEFIX   ', 'W/m2', 1,    'A', 'Total energy after fixer',         phys_decomp)
    call addfld('DTCORE'  , 'K/s' , pver, 'A', 'T tendency due to dynamical core', phys_decomp)

    if (masterproc) then
       write (iulog,*) ' print_energy_errors is set', print_energy_errors
    endif

    if ( history_budget ) then
       call add_default ('DTCORE   '  , history_budget_histfile_num, ' ')
    end if

  end subroutine check_energy_init

!===============================================================================

  subroutine check_energy_timestep_init(state, tend, pbuf, col_type)
    use iso_c_binding, only: c_double, c_int64_t, c_ptr, c_loc
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field
!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state),   target, intent(inout)    :: state
    type(physics_tend ),   intent(inout)    :: tend
    type(physics_buffer_desc), pointer      :: pbuf(:)
    integer, optional                       :: col_type  ! Flag inidicating whether using grid or subcolumns
!---------------------------Local storage-------------------------------

    real(r8), target :: ke(state%ncol)             ! vertical integral of kinetic energy
    real(r8), target :: se(state%ncol)             ! vertical integral of static energy
    real(r8), target :: wv(state%ncol)             ! vertical integral of water (vapor)
    real(r8), target :: wl(state%ncol)             ! vertical integral of water (liquid)
    real(r8), target :: wi(state%ncol)             ! vertical integral of water (ice)

    integer ncol                                   ! number of atmospheric columns
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
!-----------------------------------------------------------------------

    call check_energy_batch_select_impl()
    if (use_native_energy_batch_impl) then
       call check_energy_timestep_init_native(state, tend, pbuf, col_type)
       return
    end if

    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)

    call check_energy_batch_log_entered()
    call check_energy_log_direct(check_energy_timestep_init_logged, &
         'check_energy_timestep_init direct = codon; native cnst lookup/pbuf boundary')
    call check_energy_timestep_init_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), &
         int(state%psetcols, c_int64_t), int(pcnst, c_int64_t), &
         real(latvap, c_double), real(latice, c_double), real(gravit, c_double), &
         int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), int(ixrain, c_int64_t), int(ixsnow, c_int64_t), &
         c_loc(state%u), c_loc(state%v), c_loc(state%s), c_loc(state%q), c_loc(state%pdel), &
         c_loc(ke), c_loc(se), c_loc(wv), c_loc(wl), c_loc(wi), c_loc(state%te_ini), c_loc(state%tw_ini) &
    )

    state%te_cur(:ncol) = state%te_ini(:ncol)
    state%tw_cur(:ncol) = state%tw_ini(:ncol)

! zero cummulative boundary fluxes 
    tend%te_tnd(:ncol) = 0._r8
    tend%tw_tnd(:ncol) = 0._r8

    state%count = 0

! initialize physics buffer
    if (is_first_step()) then
       call pbuf_set_field(pbuf, teout_idx, state%te_ini, col_type=col_type)
    end if

  end subroutine check_energy_timestep_init

!===============================================================================

  subroutine check_energy_timestep_init_native(state, tend, pbuf, col_type)
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field
!-----------------------------------------------------------------------
! Compute initial values of energy and water integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state),   intent(inout)    :: state
    type(physics_tend ),   intent(inout)    :: tend
    type(physics_buffer_desc), pointer      :: pbuf(:)
    integer, optional                       :: col_type  ! Flag inidicating whether using grid or subcolumns
!---------------------------Local storage-------------------------------

    real(r8) :: ke(state%ncol)                     ! vertical integral of kinetic energy
    real(r8) :: se(state%ncol)                     ! vertical integral of static energy
    real(r8) :: wv(state%ncol)                     ! vertical integral of water (vapor)
    real(r8) :: wl(state%ncol)                     ! vertical integral of water (liquid)
    real(r8) :: wi(state%ncol)                     ! vertical integral of water (ice)

    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
!-----------------------------------------------------------------------

    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)

! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
          wv(i) = wv(i) + state%q(i,k,1       )*state%pdel(i,k)/gravit
       end do
    end do

    ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
    if (ixcldliq > 1  .and.  ixcldice > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
             wi(i) = wi(i) + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    ! Don't require precip either, if microphysics doesn't add it.
    if (ixrain > 1  .and.  ixsnow > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + state%q(i,k,ixrain)*state%pdel(i,k)/gravit
             wi(i) = wi(i) + state%q(i,k,ixsnow)*state%pdel(i,k)/gravit
          end do
       end do
    end if

! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
       state%te_ini(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       state%tw_ini(i) = wv(i) + wl(i) + wi(i)

       state%te_cur(i) = state%te_ini(i)
       state%tw_cur(i) = state%tw_ini(i)
    end do

! zero cummulative boundary fluxes 
    tend%te_tnd(:ncol) = 0._r8
    tend%tw_tnd(:ncol) = 0._r8

    state%count = 0

! initialize physics buffer
    if (is_first_step()) then
       call pbuf_set_field(pbuf, teout_idx, state%te_ini, col_type=col_type)
    end if

  end subroutine check_energy_timestep_init_native

!===============================================================================

  subroutine check_energy_chng(state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)
    use iso_c_binding, only: c_double, c_int64_t, c_ptr, c_loc

!-----------------------------------------------------------------------
! Check that the energy and water change matches the boundary fluxes
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state)    , target, intent(inout) :: state
    type(physics_tend )    , target, intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), target, intent(in   ) :: flx_vap(pcols)          ! (pcols) - boundary flux of vapor         (kg/m2/s)
    real(r8), target, intent(in   ) :: flx_cnd(pcols)          ! (pcols) -boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), target, intent(in   ) :: flx_ice(pcols)          ! (pcols) -boundary flux of ice           (m/s) (snow?)
    real(r8), target, intent(in   ) :: flx_sen(pcols)          ! (pcols) -boundary flux of sensible heat (w/m2)

!******************** BAB ******************************************************
!******* Note that the precip and ice fluxes are in precip units (m/s). ********
!******* I would prefer to have kg/m2/s.                                ********
!******* I would also prefer liquid (not total) and ice fluxes          ********
!*******************************************************************************

!---------------------------Local storage-------------------------------

    real(r8), target :: ke(state%ncol)             ! vertical integral of kinetic energy
    real(r8), target :: se(state%ncol)             ! vertical integral of static energy
    real(r8), target :: wv(state%ncol)             ! vertical integral of water (vapor)
    real(r8), target :: wl(state%ncol)             ! vertical integral of water (liquid)
    real(r8), target :: wi(state%ncol)             ! vertical integral of water (ice)

    integer ncol                                   ! number of atmospheric columns
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
!-----------------------------------------------------------------------

    call check_energy_batch_select_impl()
    if (use_native_energy_batch_impl .or. print_energy_errors) then
       call check_energy_chng_native(state, tend, name, nstep, ztodt, flx_vap, flx_cnd, flx_ice, flx_sen)
       return
    end if

    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)

    call check_energy_batch_log_entered()
    call check_energy_log_direct(check_energy_chng_logged, &
         'check_energy_chng direct = codon; native print_energy_errors debug fallback')
    call check_energy_chng_codon( &
         int(ncol, c_int64_t), int(pver, c_int64_t), int(state%psetcols, c_int64_t), &
         real(latvap, c_double), real(latice, c_double), real(gravit, c_double), &
         int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), int(ixrain, c_int64_t), int(ixsnow, c_int64_t), &
         c_loc(state%u), c_loc(state%v), c_loc(state%s), c_loc(state%q), c_loc(state%pdel), &
         c_loc(flx_vap), c_loc(flx_cnd), c_loc(flx_ice), c_loc(flx_sen), &
         c_loc(ke), c_loc(se), c_loc(wv), c_loc(wl), c_loc(wi), &
         c_loc(tend%te_tnd), c_loc(tend%tw_tnd), c_loc(state%te_cur), c_loc(state%tw_cur) &
    )

  end subroutine check_energy_chng


!===============================================================================

  subroutine check_energy_chng_native(state, tend, name, nstep, ztodt,        &
       flx_vap, flx_cnd, flx_ice, flx_sen)

!-----------------------------------------------------------------------
! Check that the energy and water change matches the boundary fluxes
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(inout) :: state
    type(physics_tend )    , intent(inout) :: tend
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: flx_vap(pcols)      ! (pcols) - boundary flux of vapor         (kg/m2/s)
    real(r8), intent(in   ) :: flx_cnd(pcols)      ! (pcols) -boundary flux of liquid+ice    (m/s) (precip?)
    real(r8), intent(in   ) :: flx_ice(pcols)      ! (pcols) -boundary flux of ice           (m/s) (snow?)
    real(r8), intent(in   ) :: flx_sen(pcols)      ! (pcols) -boundary flux of sensible heat (w/m2)

!******************** BAB ******************************************************
!******* Note that the precip and ice fluxes are in precip units (m/s). ********
!******* I would prefer to have kg/m2/s.                                ********
!******* I would also prefer liquid (not total) and ice fluxes          ********
!*******************************************************************************

!---------------------------Local storage-------------------------------

    real(r8) :: te_xpd(state%ncol)                 ! expected value (f0 + dt*boundary_flux)
    real(r8) :: te_dif(state%ncol)                 ! energy of input state - original energy
    real(r8) :: te_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: te_rer(state%ncol)                 ! relative error in energy column

    real(r8) :: tw_xpd(state%ncol)                 ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tw_dif(state%ncol)                 ! tw_inp - original water
    real(r8) :: tw_tnd(state%ncol)                 ! tendency from last process
    real(r8) :: tw_rer(state%ncol)                 ! relative error in water column

    real(r8) :: ke(state%ncol)                     ! vertical integral of kinetic energy
    real(r8) :: se(state%ncol)                     ! vertical integral of static energy
    real(r8) :: wv(state%ncol)                     ! vertical integral of water (vapor)
    real(r8) :: wl(state%ncol)                     ! vertical integral of water (liquid)
    real(r8) :: wi(state%ncol)                     ! vertical integral of water (ice)

    real(r8) :: te(state%ncol)                     ! vertical integral of total energy
    real(r8) :: tw(state%ncol)                     ! vertical integral of total water

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
!-----------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)

    ! Compute vertical integrals of dry static energy and water (vapor, liquid, ice)
    ke = 0._r8
    se = 0._r8
    wv = 0._r8
    wl = 0._r8
    wi = 0._r8
    do k = 1, pver
       do i = 1, ncol
          ke(i) = ke(i) + 0.5_r8*(state%u(i,k)**2 + state%v(i,k)**2)*state%pdel(i,k)/gravit
          se(i) = se(i) + state%s(i,k         )*state%pdel(i,k)/gravit
          wv(i) = wv(i) + state%q(i,k,1       )*state%pdel(i,k)/gravit
       end do
    end do

    ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
    if (ixcldliq > 1  .and.  ixcldice > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + state%q(i,k,ixcldliq)*state%pdel(i,k)/gravit
             wi(i) = wi(i) + state%q(i,k,ixcldice)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    ! Don't require precip either, if microphysics doesn't add it.
    if (ixrain > 1  .and.  ixsnow > 1) then
       do k = 1, pver
          do i = 1, ncol
             wl(i) = wl(i) + state%q(i,k,ixrain)*state%pdel(i,k)/gravit
             wi(i) = wi(i) + state%q(i,k,ixsnow)*state%pdel(i,k)/gravit
          end do
       end do
    end if

    ! Compute vertical integrals of frozen static energy and total water.
    do i = 1, ncol
       te(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
       tw(i) = wv(i) + wl(i) + wi(i)
    end do

    ! compute expected values and tendencies
    do i = 1, ncol
       ! change in static energy and total water
       te_dif(i) = te(i) - state%te_cur(i)
       tw_dif(i) = tw(i) - state%tw_cur(i)

       ! expected tendencies from boundary fluxes for last process
       te_tnd(i) = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._r8*latice + flx_sen(i)
       tw_tnd(i) = flx_vap(i) - flx_cnd(i) *1000._r8

       ! cummulative tendencies from boundary fluxes
       tend%te_tnd(i) = tend%te_tnd(i) + te_tnd(i)
       tend%tw_tnd(i) = tend%tw_tnd(i) + tw_tnd(i)

       ! expected new values from previous state plus boundary fluxes
       te_xpd(i) = state%te_cur(i) + te_tnd(i)*ztodt
       tw_xpd(i) = state%tw_cur(i) + tw_tnd(i)*ztodt

       ! relative error, expected value - input state / previous state 
       te_rer(i) = (te_xpd(i) - te(i)) / state%te_cur(i)
    end do

    ! relative error for total water (allow for dry atmosphere)
    tw_rer = 0._r8
    where (state%tw_cur(:ncol) > 0._r8) 
       tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / state%tw_cur(:ncol)
    end where

    ! error checking
    if (print_energy_errors) then
       if (any(abs(te_rer(1:ncol)) > 1.E-14_r8 .or. abs(tw_rer(1:ncol)) > 1.E-10_r8)) then
          do i = 1, ncol
             ! the relative error threshold for the water budget has been reduced to 1.e-10
             ! to avoid messages generated by QNEG3 calls
             ! PJR- change to identify if error in energy or water 
             if (abs(te_rer(i)) > 1.E-14_r8 ) then 
                state%count = state%count + 1
                write(iulog,*) "significant energy conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) te(i),te_xpd(i),te_dif(i),tend%te_tnd(i)*ztodt,  &
                      te_tnd(i)*ztodt,te_rer(i)
             endif
             if ( abs(tw_rer(i)) > 1.E-10_r8) then
                state%count = state%count + 1
                write(iulog,*) "significant water conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
                write(iulog,*) tw(i),tw_xpd(i),tw_dif(i),tend%tw_tnd(i)*ztodt,  &
                      tw_tnd(i)*ztodt,tw_rer(i)
             end if
          end do
       end if
    end if

    ! copy new value to state
    do i = 1, ncol
       state%te_cur(i) = te(i)
       state%tw_cur(i) = tw(i)
    end do

  end subroutine check_energy_chng_native


!===============================================================================
  subroutine check_energy_gmean(state, pbuf2d, dtime, nstep)
    use iso_c_binding, only: c_double, c_int64_t, c_ptr, c_loc

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    
!-----------------------------------------------------------------------
! Compute global mean total energy of physics input and output states
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), target, intent(in   ), dimension(begchunk:endchunk) :: state
    type(physics_buffer_desc),    pointer    :: pbuf2d(:,:)

    real(r8), intent(in) :: dtime        ! physics time step
    integer , intent(in) :: nstep        ! current timestep number

!---------------------------Local storage-------------------------------
    integer :: ncol                      ! number of active columns
    integer :: lchnk                     ! chunk index

    real(r8), target :: te(pcols,begchunk:endchunk,3)
                                         ! total energy of input/output states (copy)
    real(r8) :: te_glob(3)               ! global means of total energy
    real(r8), pointer :: teout(:)
!-----------------------------------------------------------------------

    call check_energy_batch_select_impl()
    if (use_native_energy_batch_impl) then
       call check_energy_gmean_native(state, pbuf2d, dtime, nstep)
       return
    end if

    ! Copy total energy out of input and output states
!DIR$ CONCURRENT
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk),teout_idx, teout)
      if (ncol > 0) then
         call check_energy_batch_log_entered()
         call check_energy_gmean_codon( &
              int(ncol, c_int64_t), &
              c_loc(state(lchnk)%te_ini(1)), c_loc(teout(1)), c_loc(state(lchnk)%pint(1,pver+1)), &
              c_loc(te(1,lchnk,1)), c_loc(te(1,lchnk,2)), c_loc(te(1,lchnk,3)) &
         )
      end if
    end do

    ! Compute global means of input and output energies and of
    ! surface pressure for heating rate (assume uniform ptop)
    call gmean(te, te_glob, 3)

    if (begchunk .le. endchunk) then
       teinp_glob = te_glob(1)
       teout_glob = te_glob(2)
       psurf_glob = te_glob(3)
       ptopb_glob = state(begchunk)%pint(1,1)

       ! Global mean total energy difference
       tedif_glob =  teinp_glob - teout_glob
       heat_glob  = -tedif_glob/dtime * gravit / (psurf_glob - ptopb_glob)

       if (masterproc) then
          write(iulog,'(1x,a9,1x,i8,4(1x,e25.17))') "nstep, te", nstep, teinp_glob, teout_glob, heat_glob, psurf_glob
       end if
    else
       heat_glob = 0._r8
    end if  !  (begchunk .le. endchunk)
    
  end subroutine check_energy_gmean

!===============================================================================
  subroutine check_energy_gmean_native(state, pbuf2d, dtime, nstep)

    use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
    
!-----------------------------------------------------------------------
! Compute global mean total energy of physics input and output states
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(in   ), dimension(begchunk:endchunk) :: state
    type(physics_buffer_desc),    pointer    :: pbuf2d(:,:)

    real(r8), intent(in) :: dtime        ! physics time step
    integer , intent(in) :: nstep        ! current timestep number

!---------------------------Local storage-------------------------------
    integer :: ncol                      ! number of active columns
    integer :: lchnk                     ! chunk index

    real(r8) :: te(pcols,begchunk:endchunk,3)   
                                         ! total energy of input/output states (copy)
    real(r8) :: te_glob(3)               ! global means of total energy
    real(r8), pointer :: teout(:)
!-----------------------------------------------------------------------

    ! Copy total energy out of input and output states
!DIR$ CONCURRENT
    do lchnk = begchunk, endchunk
       ncol = state(lchnk)%ncol
       ! input energy
       te(:ncol,lchnk,1) = state(lchnk)%te_ini(:ncol)
       ! output energy
       call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk),teout_idx, teout)

       te(:ncol,lchnk,2) = teout(1:ncol)
       ! surface pressure for heating rate
       te(:ncol,lchnk,3) = state(lchnk)%pint(:ncol,pver+1)
    end do

    ! Compute global means of input and output energies and of
    ! surface pressure for heating rate (assume uniform ptop)
    call gmean(te, te_glob, 3)

    if (begchunk .le. endchunk) then
       teinp_glob = te_glob(1)
       teout_glob = te_glob(2)
       psurf_glob = te_glob(3)
       ptopb_glob = state(begchunk)%pint(1,1)

       ! Global mean total energy difference
       tedif_glob =  teinp_glob - teout_glob
       heat_glob  = -tedif_glob/dtime * gravit / (psurf_glob - ptopb_glob)

       if (masterproc) then
          write(iulog,'(1x,a9,1x,i8,4(1x,e25.17))') "nstep, te", nstep, teinp_glob, teout_glob, heat_glob, psurf_glob
       end if
    else
       heat_glob = 0._r8
    end if  !  (begchunk .le. endchunk)
    
  end subroutine check_energy_gmean_native

!===============================================================================
  subroutine check_energy_fix(state, ptend, nstep, eshflx)
    use iso_c_binding, only: c_double, c_int64_t, c_ptr, c_loc

!-----------------------------------------------------------------------
! Add heating rate required for global mean total energy conservation
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), target, intent(in   ) :: state
    type(physics_ptend), target, intent(inout) :: ptend

    integer , intent(in   ) :: nstep          ! time step number
    real(r8), target, intent(inout  ) :: eshflx(pcols)  ! effective sensible heat flux

!---------------------------Local storage-------------------------------
    integer  :: ncol                     ! number of atmospheric columns in chunk
    integer  :: status, n, i, k
    integer  :: first_i_ptend, first_k_ptend, first_i_esh, first_i_pint, first_k_pint
    integer  :: first_i_hflux_srf, first_i_hflux_top
    integer  :: first_i_state_s, first_k_state_s, first_i_te_cur
    logical  :: debug_compare, ptend_mismatch, esh_mismatch, pint_mismatch
    logical  :: hflux_srf_mismatch, hflux_top_mismatch
    logical  :: state_s_mismatch, te_cur_mismatch
    logical, save :: debug_announced = .false.
    character(len=32) :: debug_env
    real(r8), allocatable :: ptend_s_dbg(:,:), pint_dbg(:,:), state_s_dbg(:,:), te_cur_dbg(:)
    real(r8) :: eshflx_dbg(pcols)
    real(r8) :: max_abs_s, max_abs_esh, max_abs_pint, diff
    real(r8) :: max_abs_hflux_srf, max_abs_hflux_top, max_abs_state_s, max_abs_te_cur
!-----------------------------------------------------------------------
    ncol = state%ncol
    debug_compare = .false.
    debug_env = ''

    call check_energy_batch_select_impl()
    call get_environment_variable('CHECK_ENERGY_FIX_DEBUG', value=debug_env, length=n, status=status)
    if (status == 0 .and. n > 0) then
       debug_compare = trim(adjustl(debug_env(:n))) /= '0'
    end if
    if (debug_compare .and. .not. debug_announced .and. masterproc) then
       write(iulog,*) 'check_energy_fix debug compare = enabled'
       debug_announced = .true.
    end if

    call physics_ptend_init(ptend, state%psetcols, 'chkenergyfix', ls=.true.)
    eshflx(:) = 0._r8

#if ( defined OFFLINE_DYN )
    ! disable the energy fix for offline driver
    heat_glob = 0._r8
#endif

    if (use_native_energy_batch_impl) then
       call check_energy_fix_native(state, ptend, nstep, eshflx)
       return
    end if

    if (debug_compare) then
       allocate(ptend_s_dbg(state%psetcols,pver), stat=status)
       if (status /= 0) call endrun('check_energy_fix debug allocate failed: ptend_s_dbg')
       allocate(pint_dbg(state%psetcols,pver+1), stat=status)
       if (status /= 0) call endrun('check_energy_fix debug allocate failed: pint_dbg')
       allocate(state_s_dbg(state%psetcols,pver), stat=status)
       if (status /= 0) call endrun('check_energy_fix debug allocate failed: state_s_dbg')
       allocate(te_cur_dbg(state%psetcols), stat=status)
       if (status /= 0) call endrun('check_energy_fix debug allocate failed: te_cur_dbg')
       ptend_s_dbg(:,:) = 0._r8
       pint_dbg(:,:) = state%pint(:,:)
       state_s_dbg(:,:) = state%s(:,:)
       te_cur_dbg(:) = state%te_cur(:)
       eshflx(:) = 0._r8
    end if

    call check_energy_batch_log_entered()
    call check_energy_fix_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(state%psetcols, c_int64_t), real(heat_glob, c_double), real(gravit, c_double), &
         c_loc(state%pint), c_loc(ptend%s), c_loc(eshflx) &
    )

    if (debug_compare) then
       eshflx_dbg(:) = 0._r8
       do k = 1, pver
          do i = 1, ncol
             ptend_s_dbg(i,k) = heat_glob
          end do
       end do
       do i = 1, ncol
          eshflx_dbg(i) = heat_glob * (state%pint(i,pver+1) - state%pint(i,1)) / gravit
       end do

       max_abs_s = 0._r8
       max_abs_esh = 0._r8
       max_abs_pint = 0._r8
       max_abs_hflux_srf = 0._r8
       max_abs_hflux_top = 0._r8
       max_abs_state_s = 0._r8
       max_abs_te_cur = 0._r8
       first_i_ptend = 0
       first_k_ptend = 0
       first_i_esh = 0
       first_i_pint = 0
       first_k_pint = 0
       first_i_hflux_srf = 0
       first_i_hflux_top = 0
       first_i_state_s = 0
       first_k_state_s = 0
       first_i_te_cur = 0
       ptend_mismatch = .false.
       esh_mismatch = .false.
       pint_mismatch = .false.
       hflux_srf_mismatch = .false.
       hflux_top_mismatch = .false.
       state_s_mismatch = .false.
       te_cur_mismatch = .false.
       do k = 1, pver
          do i = 1, state%psetcols
             diff = abs(ptend%s(i,k) - ptend_s_dbg(i,k))
             if (diff > max_abs_s) max_abs_s = diff
             if (.not. ptend_mismatch .and. ptend%s(i,k) /= ptend_s_dbg(i,k)) then
                first_i_ptend = i
                first_k_ptend = k
                ptend_mismatch = .true.
             end if
          end do
       end do
       do i = 1, pcols
          diff = abs(eshflx(i) - eshflx_dbg(i))
          if (diff > max_abs_esh) max_abs_esh = diff
          if (.not. esh_mismatch .and. eshflx(i) /= eshflx_dbg(i)) then
             first_i_esh = i
             esh_mismatch = .true.
          end if
       end do
       do i = 1, state%psetcols
          diff = abs(ptend%hflux_srf(i))
          if (diff > max_abs_hflux_srf) max_abs_hflux_srf = diff
          if (.not. hflux_srf_mismatch .and. ptend%hflux_srf(i) /= 0._r8) then
             first_i_hflux_srf = i
             hflux_srf_mismatch = .true.
          end if
          diff = abs(ptend%hflux_top(i))
          if (diff > max_abs_hflux_top) max_abs_hflux_top = diff
          if (.not. hflux_top_mismatch .and. ptend%hflux_top(i) /= 0._r8) then
             first_i_hflux_top = i
             hflux_top_mismatch = .true.
          end if
       end do
       do k = 1, pver + 1
          do i = 1, state%psetcols
             diff = abs(state%pint(i,k) - pint_dbg(i,k))
             if (diff > max_abs_pint) max_abs_pint = diff
             if (.not. pint_mismatch .and. state%pint(i,k) /= pint_dbg(i,k)) then
                first_i_pint = i
                first_k_pint = k
                pint_mismatch = .true.
             end if
          end do
       end do
       do k = 1, pver
          do i = 1, state%psetcols
             diff = abs(state%s(i,k) - state_s_dbg(i,k))
             if (diff > max_abs_state_s) max_abs_state_s = diff
             if (.not. state_s_mismatch .and. state%s(i,k) /= state_s_dbg(i,k)) then
                first_i_state_s = i
                first_k_state_s = k
                state_s_mismatch = .true.
             end if
          end do
       end do
       do i = 1, state%psetcols
          diff = abs(state%te_cur(i) - te_cur_dbg(i))
          if (diff > max_abs_te_cur) max_abs_te_cur = diff
          if (.not. te_cur_mismatch .and. state%te_cur(i) /= te_cur_dbg(i)) then
             first_i_te_cur = i
             te_cur_mismatch = .true.
          end if
       end do

       if (ptend_mismatch .or. esh_mismatch .or. pint_mismatch .or. hflux_srf_mismatch .or. hflux_top_mismatch .or. &
            state_s_mismatch .or. te_cur_mismatch) then

          write(iulog,*) 'check_energy_fix debug mismatch at nstep=', nstep, ' lchnk=', state%lchnk
          write(iulog,*) '  max_abs_ptend_s=', max_abs_s, ' max_abs_eshflx=', max_abs_esh, &
               ' max_abs_pint=', max_abs_pint, ' max_abs_hflux_srf=', max_abs_hflux_srf, &
               ' max_abs_hflux_top=', max_abs_hflux_top, ' max_abs_state_s=', max_abs_state_s, &
               ' max_abs_te_cur=', max_abs_te_cur
          if (ptend_mismatch) then
             write(iulog,*) '  first ptend%s mismatch i,k=', first_i_ptend, first_k_ptend, &
                  ptend%s(first_i_ptend,first_k_ptend), ptend_s_dbg(first_i_ptend,first_k_ptend)
          end if
          if (esh_mismatch) then
             write(iulog,*) '  first eshflx mismatch i=', first_i_esh, eshflx(first_i_esh), eshflx_dbg(first_i_esh)
          end if
          if (pint_mismatch) then
             write(iulog,*) '  first pint mismatch i,k=', first_i_pint, first_k_pint, &
                  state%pint(first_i_pint,first_k_pint), pint_dbg(first_i_pint,first_k_pint)
          end if
          if (hflux_srf_mismatch) then
             write(iulog,*) '  first hflux_srf mismatch i=', first_i_hflux_srf, ptend%hflux_srf(first_i_hflux_srf)
          end if
          if (hflux_top_mismatch) then
             write(iulog,*) '  first hflux_top mismatch i=', first_i_hflux_top, ptend%hflux_top(first_i_hflux_top)
          end if
          if (state_s_mismatch) then
             write(iulog,*) '  first state%s mismatch i,k=', first_i_state_s, first_k_state_s, &
                  state%s(first_i_state_s,first_k_state_s), state_s_dbg(first_i_state_s,first_k_state_s)
          end if
          if (te_cur_mismatch) then
             write(iulog,*) '  first te_cur mismatch i=', first_i_te_cur, state%te_cur(first_i_te_cur), &
                  te_cur_dbg(first_i_te_cur)
          end if
          flush(iulog)
          deallocate(ptend_s_dbg, pint_dbg, state_s_dbg, te_cur_dbg)
          call endrun('check_energy_fix debug mismatch')
       end if
       deallocate(ptend_s_dbg, pint_dbg, state_s_dbg, te_cur_dbg)
    end if

    return
  end subroutine check_energy_fix

!===============================================================================
  subroutine check_energy_fix_native(state, ptend, nstep, eshflx)

!-----------------------------------------------------------------------
! Add heating rate required for global mean total energy conservation
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

    type(physics_state), intent(in   ) :: state
    type(physics_ptend), intent(inout) :: ptend

    integer , intent(in   ) :: nstep          ! time step number
    real(r8), intent(inout  ) :: eshflx(pcols)  ! effective sensible heat flux

!---------------------------Local storage-------------------------------
    integer  :: i                        ! column
    integer  :: ncol                     ! number of atmospheric columns in chunk
!-----------------------------------------------------------------------
    ncol = state%ncol

! add (-) global mean total energy difference as heating
    ptend%s(:ncol,:pver) = heat_glob

! compute effective sensible heat flux
    eshflx(:) = 0._r8
    do i = 1, ncol
       eshflx(i) = heat_glob * (state%pint(i,pver+1) - state%pint(i,1)) / gravit
    end do

    return
  end subroutine check_energy_fix_native


!===============================================================================
  subroutine check_tracers_init(state, tracerint)
    use iso_c_binding, only: c_double, c_int64_t

!-----------------------------------------------------------------------
! Compute initial values of tracers integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    type(physics_state),   intent(in)    :: state
    type(check_tracers_data), intent(out)   :: tracerint

!-----------------------------------------------------------------------

    call check_tracers_batch_select_impl()
    if (use_native_tracers_batch_impl) then
       call check_tracers_init_native(state, tracerint)
       return
    end if

    call check_tracers_batch_log_entered()
    call check_tracers_init_codon()

    return
  end subroutine check_tracers_init

!===============================================================================
  subroutine check_tracers_init_native(state, tracerint)

!-----------------------------------------------------------------------
! Compute initial values of tracers integrals, 
! zero cumulative tendencies
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------

    type(physics_state),   intent(in)    :: state
    type(check_tracers_data), intent(out)   :: tracerint

!---------------------------Local storage-------------------------------

    real(r8) :: tr(pcols)                          ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                ! pdel for tracer

    integer ncol                                   ! number of atmospheric columns
    integer  i,k,m                                 ! column, level,constituent indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices

!-----------------------------------------------------------------------

    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)

    do m = 1,pcnst

       if ( any(m == (/ 1, ixcldliq, ixcldice, &
                           ixrain,   ixsnow    /)) ) exit   ! dont process water substances
                                                            ! they are checked in check_energy
       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals of tracer
       tr = 0._r8
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
          end do
       end do

       ! Compute vertical integrals of frozen static tracers and total water.
       do i = 1, ncol
          tracerint%tracer(i,m) = tr(i)
       end do

       ! zero cummulative boundary fluxes 
       tracerint%tracer_tnd(:ncol,m) = 0._r8

       tracerint%count(m) = 0

    end do

    return
  end subroutine check_tracers_init_native

!===============================================================================
  subroutine check_tracers_chng(state, tracerint, name, nstep, ztodt, cflx)
    use iso_c_binding, only: c_double, c_int64_t

!-----------------------------------------------------------------------
! Check that the tracers and water change matches the boundary fluxes
! these checks are not save when there are tracers transformations, as 
! they only check to see whether a mass change in the column is
! associated with a flux
!-----------------------------------------------------------------------

    use cam_abortutils, only: endrun 


    implicit none

!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(in   ) :: state
    type(check_tracers_data), intent(inout) :: tracerint! tracers integrals and boundary fluxes
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: cflx(pcols,pcnst)       ! boundary flux of tracers       (kg/m2/s)

    call check_tracers_batch_select_impl()
    if (use_native_tracers_batch_impl) then
       call check_tracers_chng_native(state, tracerint, name, nstep, ztodt, cflx)
       return
    end if

    call check_tracers_batch_log_entered()
    call check_tracers_chng_codon()

    return
  end subroutine check_tracers_chng

!===============================================================================
  subroutine check_tracers_chng_native(state, tracerint, name, nstep, ztodt, cflx)

!-----------------------------------------------------------------------
! Check that the tracers and water change matches the boundary fluxes
! these checks are not save when there are tracers transformations, as 
! they only check to see whether a mass change in the column is
! associated with a flux
!-----------------------------------------------------------------------

    use cam_abortutils, only: endrun 


    implicit none

!------------------------------Arguments--------------------------------

    type(physics_state)    , intent(in   ) :: state
    type(check_tracers_data), intent(inout) :: tracerint! tracers integrals and boundary fluxes
    character*(*),intent(in) :: name               ! parameterization name for fluxes
    integer , intent(in   ) :: nstep               ! current timestep number
    real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
    real(r8), intent(in   ) :: cflx(pcols,pcnst)       ! boundary flux of tracers       (kg/m2/s)

!---------------------------Local storage-------------------------------

    real(r8) :: tracer_inp(pcols,pcnst)                   ! total tracer of new (input) state
    real(r8) :: tracer_xpd(pcols,pcnst)                   ! expected value (w0 + dt*boundary_flux)
    real(r8) :: tracer_dif(pcols,pcnst)                   ! tracer_inp - original tracer
    real(r8) :: tracer_tnd(pcols,pcnst)                   ! tendency from last process
    real(r8) :: tracer_rer(pcols,pcnst)                   ! relative error in tracer column

    real(r8) :: tr(pcols)                           ! vertical integral of tracer
    real(r8) :: trpdel(pcols, pver)                       ! pdel for tracer

    integer lchnk                                  ! chunk identifier
    integer ncol                                   ! number of atmospheric columns
    integer  i,k                                   ! column, level indices
    integer :: ixcldice, ixcldliq                  ! CLDICE and CLDLIQ indices
    integer :: ixrain, ixsnow                      ! RAINQM and SNOWQM indices
    integer :: m                            ! tracer index
    character(len=8) :: tracname   ! tracername
!-----------------------------------------------------------------------
!!$    if (.true.) return

    lchnk = state%lchnk
    ncol  = state%ncol
    call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
    call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
    call cnst_get_ind('RAINQM', ixrain,   abort=.false.)
    call cnst_get_ind('SNOWQM', ixsnow,   abort=.false.)

    do m = 1,pcnst

       if ( any(m == (/ 1, ixcldliq, ixcldice, &
                           ixrain,   ixsnow    /)) ) exit   ! dont process water substances
                                                            ! they are checked in check_energy

       tracname = cnst_name(m)
       if (cnst_get_type_byind(m).eq.'dry') then
          trpdel(:ncol,:) = state%pdeldry(:ncol,:)
       else
          trpdel(:ncol,:) = state%pdel(:ncol,:)
       endif

       ! Compute vertical integrals tracers
       tr = 0._r8
       do k = 1, pver
          do i = 1, ncol
             tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
          end do
       end do

       ! Compute vertical integrals of tracer
       do i = 1, ncol
          tracer_inp(i,m) = tr(i)
       end do

       ! compute expected values and tendencies
       do i = 1, ncol
          ! change in tracers 
          tracer_dif(i,m) = tracer_inp(i,m) - tracerint%tracer(i,m)

          ! expected tendencies from boundary fluxes for last process
          tracer_tnd(i,m) = cflx(i,m)

          ! cummulative tendencies from boundary fluxes
          tracerint%tracer_tnd(i,m) = tracerint%tracer_tnd(i,m) + tracer_tnd(i,m)

          ! expected new values from original values plus boundary fluxes
          tracer_xpd(i,m) = tracerint%tracer(i,m) + tracerint%tracer_tnd(i,m)*ztodt

          ! relative error, expected value - input value / original 
          tracer_rer(i,m) = (tracer_xpd(i,m) - tracer_inp(i,m)) / tracerint%tracer(i,m)
       end do

!! final loop for error checking
!    do i = 1, ncol

!! error messages
!       if (abs(enrgy_rer(i)) > 1.E-14 .or. abs(water_rer(i)) > 1.E-14) then
!          tracerint%count = tracerint%count + 1
!          write(iulog,*) "significant conservations error after ", name,        &
!               " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col", i
!          write(iulog,*) enrgy_inp(i),enrgy_xpd(i),enrgy_dif(i),tracerint%enrgy_tnd(i)*ztodt,  &
!               enrgy_tnd(i)*ztodt,enrgy_rer(i)
!          write(iulog,*) water_inp(i),water_xpd(i),water_dif(i),tracerint%water_tnd(i)*ztodt,  &
!               water_tnd(i)*ztodt,water_rer(i)
!       end if
!    end do


       ! final loop for error checking
       if ( maxval(tracer_rer) > 1.E-14_r8 ) then
          write(iulog,*) "CHECK_TRACERS TRACER large rel error"
          write(iulog,*) tracer_rer
       endif

       do i = 1, ncol
          ! error messages
          if (abs(tracer_rer(i,m)) > 1.E-14_r8 ) then
             tracerint%count = tracerint%count + 1
             write(iulog,*) "CHECK_TRACERS TRACER significant conservation error after ", name,        &
                  " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col",i
             write(iulog,*)' process name, tracname, index ',  name, tracname, m
             write(iulog,*)" input integral              ",tracer_inp(i,m)
             write(iulog,*)" expected integral           ", tracer_xpd(i,m)
             write(iulog,*)" input - inital integral     ",tracer_dif(i,m)
             write(iulog,*)" cumulative tend      ",tracerint%tracer_tnd(i,m)*ztodt
             write(iulog,*)" process tend         ",tracer_tnd(i,m)*ztodt
             write(iulog,*)" relative error       ",tracer_rer(i,m)
             call endrun()
          end if
       end do
    end do

    return
  end subroutine check_tracers_chng_native


end module check_energy
