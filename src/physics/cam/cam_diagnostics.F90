module cam_diagnostics

!---------------------------------------------------------------------------------
! Module to compute a variety of diagnostics quantities for history files
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use camsrfexch,     only: cam_in_t, cam_out_t
use physics_types,  only: physics_state, physics_tend
use ppgrid,         only: pcols, pver, pverp, begchunk, endchunk
use physics_buffer, only: physics_buffer_desc, pbuf_add_field, dtype_r8, dyn_time_lvls, &
                          pbuf_get_field, pbuf_get_index, pbuf_old_tim_idx



use cam_history,    only: outfld, write_inithist, hist_fld_active
use cam_logfile,    only: iulog
use constituents,   only: pcnst, cnst_name, cnst_longname, cnst_cam_outfld, ptendnam, dmetendnam, apcnst, bpcnst, &
                          cnst_get_ind
use dycore,         only: dycore_is
use phys_control,   only: phys_getopts
use wv_saturation,  only: qsat, qsat_water, svp_ice
use time_manager,   only: is_first_step

use scamMod,        only: single_column, wfld
use cam_abortutils, only: endrun
use spmd_utils,     only: masterproc

use water_tracer_vars, only: trace_water, wtrc_nwset, wtrc_iatype, wtrc_srfvap_names, wtrc_srfpcp_indices, &
                             wtrc_out_names
use water_types,       only: iwtvap, iwtcvrain, iwtcvsnow, iwtstrain, iwtstsnow

implicit none
private
save

! Public interfaces

public :: &
   diag_register,      &! register pbuf space
   diag_init,          &! initialization
   diag_allocate,      &! allocate memory for module variables
   diag_deallocate,    &! deallocate memory for module variables
   diag_conv_tend_ini, &! initialize convective tendency calcs
   diag_phys_writeout, &! output diagnostics of the dynamics
   diag_phys_tend_writeout, & ! output physics tendencies
   diag_state_b4_phys_write,& ! output state before physics execution
   diag_conv,          &! output diagnostics of convective processes
   diag_surf,          &! output diagnostics of the surface
   diag_export,        &! output export state
   diag_physvar_ic,    &
   diag_readnl          ! read namelist options

logical, public :: inithist_all = .false. ! Flag to indicate set of fields to be
                                          ! included on IC file
                                          !  .false.  include only required fields
                                          !  .true.   include required *and* optional fields

! Private data

integer :: dqcond_num                     ! number of constituents to compute convective
character(len=16) :: dcconnam(pcnst)      ! names of convection tendencies
                                          ! tendencies for
real(r8), allocatable, target :: dtcond(:,:,:)    ! temperature tendency due to convection
type dqcond_t
   real(r8), allocatable :: cnst(:,:,:)   ! constituent tendency due to convection
end type dqcond_t
type(dqcond_t), allocatable :: dqcond(:)

character(len=8) :: diag_cnst_conv_tend = 'q_only' ! output constituent tendencies due to convection
                                                   ! 'none', 'q_only' or 'all'

logical          :: history_amwg                   ! output the variables used by the AMWG diag package
logical          :: history_vdiag                  ! output the variables used by the AMWG variability diag package
logical          :: history_eddy                   ! output the eddy variables
logical          :: history_budget                 ! output tendencies and state variables for CAM4
                                                   ! temperature, water vapor, cloud ice and cloud
                                                   ! liquid budgets.
integer          :: history_budget_histfile_num    ! output history file number for budget fields
logical          :: history_waccm                  ! outputs typically used for WACCM

!Physics buffer indices
integer  ::      qcwat_idx  = 0
integer  ::      tcwat_idx  = 0
integer  ::      lcwat_idx  = 0
integer  ::      cld_idx    = 0
integer  ::      concld_idx = 0
integer  ::      tke_idx    = 0
integer  ::      kvm_idx    = 0
integer  ::      kvh_idx    = 0
integer  ::      cush_idx   = 0
integer  ::      t_ttend_idx = 0

integer  ::      prec_dp_idx  = 0
integer  ::      snow_dp_idx  = 0
integer  ::      prec_sh_idx  = 0
integer  ::      snow_sh_idx  = 0
integer  ::      prec_sed_idx = 0
integer  ::      snow_sed_idx = 0
integer  ::      prec_pcw_idx = 0
integer  ::      snow_pcw_idx = 0


integer :: tpert_idx=-1, qpert_idx=-1, pblh_idx=-1
logical :: diag_surf_use_native_impl = .false.
logical :: diag_surf_impl_selected = .false.
logical :: diag_phys_writeout_use_native_impl = .false.
logical :: diag_phys_writeout_impl_selected = .false.
logical :: diag_phys_writeout_batch_entered_logged = .false.
logical :: diag_physvar_ic_use_native_impl = .false.
logical :: diag_physvar_ic_impl_selected = .false.
logical :: diag_physvar_ic_logged = .false.
logical :: diag_phys_tend_use_native_impl = .false.
logical :: diag_phys_tend_impl_selected = .false.
logical :: diag_phys_tend_entered_logged = .false.
logical :: cam_diag_conv_batch_use_native_impl = .false.
logical :: cam_diag_conv_batch_impl_selected = .false.
logical :: cam_diag_conv_tend_ini_entered_logged = .false.
logical :: cam_diag_conv_entered_logged = .false.
logical :: cam_diag_conv_precip_dtcond_entered_logged = .false.
logical :: cam_diag_init_helpers_use_native_impl = .false.
logical :: cam_diag_init_helpers_impl_selected = .false.
logical :: cam_diag_init_helpers_entered_logged = .false.
logical :: diag_readnl_use_native_impl = .false.
logical :: diag_readnl_impl_selected = .false.
logical :: diag_readnl_logged = .false.
logical :: cam_diag_parent_use_native_impl = .false.
logical :: cam_diag_parent_impl_selected = .false.
logical :: diag_register_logged = .false.
logical :: diag_init_logged = .false.
logical :: diag_allocate_logged = .false.
logical :: diag_deallocate_logged = .false.
logical :: diag_phys_writeout_logged = .false.
logical :: diag_surf_logged = .false.
logical :: diag_export_logged = .false.
logical :: diag_state_b4_phys_write_logged = .false.
logical :: diag_conv_logged = .false.

interface
   function diag_readnl_codon() result(out_c) bind(c, name="diag_readnl_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t) :: out_c
   end function diag_readnl_codon
   function cam_diagnostics_touch_codon(stage_c) result(stage_out) bind(c, name="cam_diagnostics_touch_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: stage_c
      integer(c_int64_t) :: stage_out
   end function cam_diagnostics_touch_codon
   function diag_register_codon(stage_c) result(stage_out) bind(c, name="diag_register_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: stage_c
      integer(c_int64_t) :: stage_out
   end function diag_register_codon
   function diag_allocate_codon(stage_c) result(stage_out) bind(c, name="diag_allocate_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: stage_c
      integer(c_int64_t) :: stage_out
   end function diag_allocate_codon
   function diag_deallocate_codon(stage_c) result(stage_out) bind(c, name="diag_deallocate_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: stage_c
      integer(c_int64_t) :: stage_out
   end function diag_deallocate_codon
   function diag_export_codon(stage_c) result(stage_out) bind(c, name="diag_export_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: stage_c
      integer(c_int64_t) :: stage_out
   end function diag_export_codon
   function diag_state_b4_phys_write_codon(stage_c) result(stage_out) bind(c, name="diag_state_b4_phys_write_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: stage_c
      integer(c_int64_t) :: stage_out
   end function diag_state_b4_phys_write_codon
   function diag_conv_codon(stage_c) result(stage_out) bind(c, name="diag_conv_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: stage_c
      integer(c_int64_t) :: stage_out
   end function diag_conv_codon
   function diag_phys_writeout_codon(stage_c) result(stage_out) bind(c, name="diag_phys_writeout_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: stage_c
      integer(c_int64_t) :: stage_out
   end function diag_phys_writeout_codon
   subroutine diag_phys_writeout_batch_dispatch_codon(group_c, mode_c, submode_c, ncol_c, pcols_c, pver_c, &
        scalar1_c, scalar2_c, scalar3_c, a_p, b_p, c_p, d_p, e_p, f_p, out1_p, out2_p, out3_p) &
        bind(c, name="diag_phys_writeout_batch_dispatch_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: group_c, mode_c, submode_c, ncol_c, pcols_c, pver_c
      real(c_double), value :: scalar1_c, scalar2_c, scalar3_c
      type(c_ptr), value :: a_p, b_p, c_p, d_p, e_p, f_p, out1_p, out2_p, out3_p
   end subroutine diag_phys_writeout_batch_dispatch_codon
   subroutine diag_conv_batch_dispatch_codon(group_c, mode_c, ncol_c, pcols_c, pver_c, pcnst_c, m_c, &
        scalar1_c, scalar2_c, a_p, b_p, c_p, d_p, e_p, f_p, g_p, h_p, out1_p, out2_p, out3_p, out4_p, &
        out5_p, aux1_p, aux2_p) bind(c, name="diag_conv_batch_dispatch_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: group_c, mode_c, ncol_c, pcols_c, pver_c, pcnst_c, m_c
      real(c_double), value :: scalar1_c, scalar2_c
      type(c_ptr), value :: a_p, b_p, c_p, d_p, e_p, f_p, g_p, h_p
      type(c_ptr), value :: out1_p, out2_p, out3_p, out4_p, out5_p, aux1_p, aux2_p
   end subroutine diag_conv_batch_dispatch_codon
   subroutine diag_conv_scale_2d_codon(ncol_c, pcols_c, pver_c, scale_c, src_p, out_p) &
        bind(c, name="diag_conv_scale_2d_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
      real(c_double), value :: scale_c
      type(c_ptr), value :: src_p, out_p
   end subroutine diag_conv_scale_2d_codon
   subroutine diag_conv_tend_ini_codon(mode_c, ncol_c, pcols_c, pver_c, pcnst_c, m_c, src_p, dst_p) &
        bind(c, name="diag_conv_tend_ini_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c, pcnst_c, m_c
      type(c_ptr), value :: src_p, dst_p
   end subroutine diag_conv_tend_ini_codon
end interface

contains

! ===============================================================================

subroutine diag_register

   use iso_c_binding, only: c_int64_t

   call cam_diag_touch_and_log(1_c_int64_t, diag_register_logged, &
        'diag_register direct = codon; register selector/touch direct = codon; pbuf_add_field native CAM API island')

   ! Request physics buffer space for fields that persist across timesteps.
   call pbuf_add_field('T_TTEND', 'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), t_ttend_idx)

end subroutine diag_register

!===============================================================================

subroutine diag_readnl(nlfile)
  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand
  use spmd_utils,      only: masterproc
  use iso_c_binding,   only: c_int64_t

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  integer(c_int64_t) :: out_c
  logical :: group_found
  character(len=*), parameter :: subname = 'diag_readnl'

  namelist /cam_diag_opts/ diag_cnst_conv_tend
  !-----------------------------------------------------------------------------

  call diag_readnl_select_impl()
  if (.not. diag_readnl_use_native_impl) then
     group_found = .false.
     if (masterproc) then
        unitn = getunit()
        open( unitn, file=trim(nlfile), status='old' )
        call find_group_name(unitn, 'cam_diag_opts', status=ierr)
        group_found = ierr == 0
        close(unitn)
        call freeunit(unitn)
     end if
#ifdef SPMD
     call mpibcast(group_found, 1, mpilog,  0, mpicom)
#endif
     if (.not. group_found) then
        out_c = diag_readnl_codon()
        if (out_c == 0_c_int64_t) then
           call diag_readnl_log_direct('diag_readnl direct = codon')
           return
        end if
     end if
  end if

  if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cam_diag_opts', status=ierr)
      if (ierr == 0) then
         read(unitn, cam_diag_opts, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(diag_cnst_conv_tend, len(diag_cnst_conv_tend), mpichar,  0, mpicom)
#endif

end subroutine diag_readnl

!================================================================================================

subroutine diag_readnl_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (diag_readnl_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CAM_DIAGNOSTICS_READNL_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      diag_readnl_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      diag_readnl_use_native_impl = .false.
   end if

   diag_readnl_impl_selected = .true.

   if (masterproc) then
      if (diag_readnl_use_native_impl) then
         write(iulog,*) 'diag_readnl implementation = native'
      else
         write(iulog,*) 'diag_readnl implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine diag_readnl_select_impl

!================================================================================================

subroutine diag_readnl_log_direct(proof_line)

   character(len=*), intent(in) :: proof_line

   if (diag_readnl_logged) return
   diag_readnl_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') trim(proof_line)
      call flush(iulog)
   end if

end subroutine diag_readnl_log_direct

!================================================================================================

subroutine cam_diag_parent_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cam_diag_parent_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CAM_DIAG_PARENT_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      cam_diag_parent_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      cam_diag_parent_use_native_impl = .false.
   end if

   cam_diag_parent_impl_selected = .true.

   if (masterproc) then
      if (cam_diag_parent_use_native_impl) then
         write(iulog,*) 'cam_diag_parent implementation = native'
      else
         write(iulog,*) 'cam_diag_parent implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine cam_diag_parent_select_impl

!================================================================================================

subroutine cam_diag_log_direct(logged, proof_line)

   logical, intent(inout) :: logged
   character(len=*), intent(in) :: proof_line

   if (logged) return
   logged = .true.

   if (masterproc) then
      write(iulog,'(A)') trim(proof_line)
      call flush(iulog)
   end if

end subroutine cam_diag_log_direct

!================================================================================================

subroutine cam_diag_touch_and_log(stage_c, logged, proof_line)

   use iso_c_binding, only: c_int64_t

   integer(c_int64_t), intent(in) :: stage_c
   logical, intent(inout) :: logged
   character(len=*), intent(in) :: proof_line
   integer(c_int64_t) :: exact_stage_c

   call cam_diag_parent_select_impl()
   if (cam_diag_parent_use_native_impl) return

   select case (stage_c)
   case (1_c_int64_t)
      exact_stage_c = diag_register_codon(stage_c)
   case (3_c_int64_t)
      exact_stage_c = diag_allocate_codon(stage_c)
   case (4_c_int64_t)
      exact_stage_c = diag_deallocate_codon(stage_c)
   case (5_c_int64_t)
      exact_stage_c = diag_export_codon(stage_c)
   case (6_c_int64_t)
      exact_stage_c = diag_state_b4_phys_write_codon(stage_c)
   case default
      exact_stage_c = cam_diagnostics_touch_codon(stage_c)
   end select

   if (exact_stage_c == stage_c) then
      call cam_diag_log_direct(logged, proof_line)
   end if

end subroutine cam_diag_touch_and_log

!================================================================================================

subroutine diag_init()

  ! Declare the history fields for which this module contains outfld calls.

   use cam_history,        only: addfld, add_default, phys_decomp
   use constituent_burden, only: constituent_burden_init
   use cam_control_mod,    only: moist_physics, ideal_phys
   use tidal_diag,         only: tidal_diag_init
   use iso_c_binding,      only: c_int64_t

   integer :: k, m
   integer :: ixcldice, ixcldliq ! constituent indices for cloud liquid and ice water.
   integer :: ierr
   integer :: conv_tend_code

   interface
      function diag_init_codon(history_budget_c, conv_tend_code_c, pcnst_c) result(result_c) &
           bind(c, name="diag_init_codon")
         use iso_c_binding, only: c_int64_t
         integer(c_int64_t), value :: history_budget_c, conv_tend_code_c, pcnst_c
         integer(c_int64_t) :: result_c
      end function diag_init_codon
   end interface

   ! outfld calls in diag_phys_writeout

   call addfld ('NSTEP   ','timestep',1,    'A','Model timestep',phys_decomp)
   call addfld ('PHIS    ','m2/s2   ',1,    'I','Surface geopotential',phys_decomp)

   call addfld ('PS      ','Pa      ',1,    'A','Surface pressure',phys_decomp)
   call addfld ('T       ','K       ',pver, 'A','Temperature',phys_decomp)
   call addfld ('U       ','m/s     ',pver, 'A','Zonal wind',phys_decomp)
   call addfld ('V       ','m/s     ',pver, 'A','Meridional wind',phys_decomp)
   call addfld (cnst_name(1),'kg/kg ',pver, 'A',cnst_longname(1),phys_decomp)

   ! State before physics
   call addfld ('TBP     ','K       ',pver, 'A','Temperature (before physics)'       ,phys_decomp)
   call addfld (bpcnst(1) ,'kg/kg   ',pver, 'A',cnst_longname(1)//' (before physics)',phys_decomp)
   ! State after physics
   call addfld ('TAP     ','K       ',pver, 'A','Temperature (after physics)'       ,phys_decomp)
   call addfld ('UAP     ','m/s     ',pver, 'A','Zonal wind (after physics)'        ,phys_decomp)
   call addfld ('VAP     ','m/s     ',pver, 'A','Meridional wind (after physics)'   ,phys_decomp)
   call addfld (apcnst(1) ,'kg/kg   ',pver, 'A',cnst_longname(1)//' (after physics)',phys_decomp)
   if ( dycore_is('LR') ) then
      call addfld ('TFIX    ','K/s     ',1,    'A'     ,'T fixer (T equivalent of Energy correction)',phys_decomp)
      call addfld ('PTTEND_RESID','K/s ',pver, 'A'     ,&
                   'T-tendency due to BAB kluge at end of tphysac (diagnostic not part of T-budget)' ,phys_decomp)
   end if
   call addfld ('TTEND_TOT   ','K/s' ,pver, 'A','Total temperature tendency'   ,phys_decomp)

   ! column burdens for all constituents except water vapor
   call constituent_burden_init

   call addfld ('Z3      ','m       ',pver, 'A','Geopotential Height (above sea level)',phys_decomp)
   call addfld ('Z1000   ','m       ',1,    'A','Geopotential Z at 1000 mbar pressure surface',phys_decomp)
   call addfld ('Z700    ','m       ',1,    'A','Geopotential Z at 700 mbar pressure surface',phys_decomp)
   ! nanr
   call addfld ('Z850    ','m       ',1,    'A','Geopotential Z at 850 mbar pressure surface',phys_decomp)
   call addfld ('Z500    ','m       ',1,    'A','Geopotential Z at 500 mbar pressure surface',phys_decomp)
   call addfld ('Z300    ','m       ',1,    'A','Geopotential Z at 300 mbar pressure surface',phys_decomp)
   call addfld ('Z250    ','m       ',1,    'A','Geopotential Z at 250 mbar pressure surface',phys_decomp)
   call addfld ('Z200    ','m       ',1,    'A','Geopotential Z at 200 mbar pressure surface',phys_decomp)
   call addfld ('Z100    ','m       ',1,    'A','Geopotential Z at 100 mbar pressure surface',phys_decomp)
   call addfld ('Z050    ','m       ',1,    'A','Geopotential Z at 50 mbar pressure surface',phys_decomp)

   call addfld ('ZZ      ','m2      ',pver, 'A','Eddy height variance' ,phys_decomp)
   call addfld ('VZ      ','m2/s    ',pver, 'A','Meridional transport of geopotential energy',phys_decomp)
   call addfld ('VT      ','K m/s   ',pver, 'A','Meridional heat transport',phys_decomp)
   ! nanr
   call addfld ('UT      ','K m/s   ',pver, 'A','Zonal heat transport',phys_decomp)
   call addfld ('VU      ','m2/s2   ',pver, 'A','Meridional flux of zonal momentum' ,phys_decomp)
   call addfld ('VV      ','m2/s2   ',pver, 'A','Meridional velocity squared' ,phys_decomp)
   ! nanr
   call addfld ('UQ      ','m/skg/kg',pver, 'A','Zonal water transport',phys_decomp)
   call addfld ('VQ      ','m/skg/kg',pver, 'A','Meridional water transport',phys_decomp)
   call addfld ('QQ      ','kg2/kg2 ',pver, 'A','Eddy moisture variance',phys_decomp)
   call addfld ('OMEGAV  ','m Pa/s2 ',pver ,'A','Vertical flux of meridional momentum' ,phys_decomp)
   call addfld ('OMGAOMGA','Pa2/s2  ',pver ,'A','Vertical flux of vertical momentum' ,phys_decomp)
   call addfld ('OMEGAQ  ','kgPa/kgs',pver ,'A','Vertical water transport' ,phys_decomp)

   call addfld ('UU      ','m2/s2   ',pver, 'A','Zonal velocity squared' ,phys_decomp)
   call addfld ('WSPEED  ','m/s     ',pver, 'X','Horizontal total wind speed maximum' ,phys_decomp)
   call addfld ('WSPDSRFMX','m/s    ',1,    'X','Horizontal total wind speed maximum at the surface' ,phys_decomp)
   call addfld ('WSPDSRFAV','m/s    ',1,    'A','Horizontal total wind speed average at the surface' ,phys_decomp)

   call addfld ('OMEGA   ','Pa/s    ',pver, 'A','Vertical velocity (pressure)',phys_decomp)
   call addfld ('OMEGAT  ','K Pa/s  ',pver, 'A','Vertical heat flux' ,phys_decomp)
   call addfld ('OMEGAU  ','m Pa/s2 ',pver, 'A','Vertical flux of zonal momentum' ,phys_decomp)
   call addfld ('OMEGA850','Pa/s    ',1,    'A','Vertical velocity at 850 mbar pressure surface',phys_decomp)
   call addfld ('OMEGA500','Pa/s    ',1,    'A','Vertical velocity at 500 mbar pressure surface',phys_decomp)
   ! nanr
   call addfld ('OMEGA200','Pa/s    ',1,    'A','Vertical velocity at 200 mbar pressure surface',phys_decomp)

   call addfld ('MQ      ','kg/m2   ',pver, 'A','Water vapor mass in layer',phys_decomp)
   call addfld ('TMQ     ','kg/m2   ',1,    'A','Total (vertically integrated) precipitable water',phys_decomp)
   call addfld ('IVT     ','kg/m/s  ',1,    'A','Total (vertically integrated) vapor transport',phys_decomp)
   call addfld ('uIVT    ','kg/m/s  ',1,    'A','u-component (vertically integrated) vapor transport',phys_decomp)
   call addfld ('vIVT    ','kg/m/s  ',1,    'A','v-component (vertically integrated) vapor transport',phys_decomp)
   call addfld ('RELHUM  ','percent ',pver, 'A','Relative humidity',phys_decomp)
   call addfld ('RHW  ','percent '   ,pver, 'A','Relative humidity with respect to liquid',phys_decomp)
   call addfld ('RHI  ','percent '   ,pver, 'A','Relative humidity with respect to ice',phys_decomp)
   call addfld ('RHCFMIP','percent ' ,pver, 'A','Relative humidity with respect to water above 273 K, ice below 273 K',phys_decomp)
   call addfld ('PSL     ','Pa      ',1,    'A','Sea level pressure',phys_decomp)

   !**********************
   !Water tracers/isotopes
   !**********************
   if(trace_water) then
     do m=1,wtrc_nwset !loop over water tracers
       call addfld ('TMQ_'//trim(wtrc_out_names(m)),'kg/m2   ',1,'A','Total (vertically integrated) precipitable water for '//trim(wtrc_out_names(m)),phys_decomp)
       call addfld ('TVQ_'//trim(wtrc_out_names(m)),'kg/m/s   ',1,'A','Total (vertically integrated) meridional flux for '//trim(wtrc_out_names(m)),phys_decomp)
       call addfld ('TUQ_'//trim(wtrc_out_names(m)),'kg/m/s   ',1,'A','Total (vertically integrated) zonal flux for '//trim(wtrc_out_names(m)),phys_decomp)
     end do
   end if
   !**********************

   call addfld ('T850    ','K       ',1,    'A','Temperature at 850 mbar pressure surface',phys_decomp)
   call addfld ('T500    ','K       ',1,    'A','Temperature at 500 mbar pressure surface',phys_decomp)
   call addfld ('T300    ','K       ',1,    'A','Temperature at 300 mbar pressure surface',phys_decomp)
   call addfld ('T200    ','K       ',1,    'A','Temperature at 200 mbar pressure surface',phys_decomp)
   call addfld ('Q850    ','kg/kg   ',1,    'A','Specific Humidity at 850 mbar pressure surface',phys_decomp)
   call addfld ('Q500    ','kg/kg   ',1,    'A','Specific Humidity at 500 mbar pressure surface',phys_decomp)
   call addfld ('Q200    ','kg/kg   ',1,    'A','Specific Humidity at 200 mbar pressure surface',phys_decomp)
   call addfld ('U925    ','m/s     ',1,    'A','Zonal wind at 925 mbar pressure surface',phys_decomp)
   call addfld ('U850    ','m/s     ',1,    'A','Zonal wind at 850 mbar pressure surface',phys_decomp)
   ! nanr
   call addfld ('U500    ','m/s     ',1,    'A','Zonal wind at 500 mbar pressure surface',phys_decomp)
   call addfld ('U600    ','m/s     ',1,    'A','Zonal wind at 600 mbar pressure surface',phys_decomp)
   call addfld ('U700    ','m/s     ',1,    'A','Zonal wind at 700 mbar pressure surface',phys_decomp)
   call addfld ('U250    ','m/s     ',1,    'A','Zonal wind at 250 mbar pressure surface',phys_decomp)
   call addfld ('U200    ','m/s     ',1,    'A','Zonal wind at 200 mbar pressure surface',phys_decomp)
   call addfld ('U010    ','m/s     ',1,    'A','Zonal wind at  10 mbar pressure surface',phys_decomp)
   call addfld ('V925    ','m/s     ',1,    'A','Meridional wind at 925 mbar pressure surface',phys_decomp)
   call addfld ('V850    ','m/s     ',1,    'A','Meridional wind at 850 mbar pressure surface',phys_decomp)
   call addfld ('V200    ','m/s     ',1,    'A','Meridional wind at 200 mbar pressure surface',phys_decomp)
   call addfld ('V250    ','m/s     ',1,    'A','Meridional wind at 250 mbar pressure surface',phys_decomp)
   ! nanr
   call addfld ('V500    ','m/s     ',1,    'A','Meridional wind at 500 mbar pressure surface',phys_decomp)
   call addfld ('V600    ','m/s     ',1,    'A','Meridional wind at 600 mbar pressure surface',phys_decomp)
   call addfld ('V700    ','m/s     ',1,    'A','Meridional wind at 700 mbar pressure surface',phys_decomp)

   call addfld ('TT      ','K2      ',pver, 'A','Eddy temperature variance' ,phys_decomp)

   call addfld ('UBOT    ','m/s     ',1,    'A','Lowest model level zonal wind',phys_decomp)
   call addfld ('VBOT    ','m/s     ',1,    'A','Lowest model level meridional wind',phys_decomp)
   call addfld ('QBOT    ','kg/kg   ',1,    'A','Lowest model level water vapor mixing ratio',phys_decomp)
   call addfld ('ZBOT    ','m       ',1,    'A','Lowest model level height', phys_decomp)

   ! Water tracers:
   ! NOTE:  may need better method for handling multiple water tracers versus
   ! just hard-coding them in, but for now, this will do, at least for water
   ! isotopes. -JN
   if(trace_water) then !Are water tracers on?
     do m=1,wtrc_nwset
       call addfld (trim(wtrc_srfvap_names(m))//'BT','kg/kg   ',1,    'A','Lowest model level mixing ratio for '//trim(wtrc_srfvap_names(m)),phys_decomp)
     end do

!     call addfld ('H2OVBT  ','kg/kg   ',1,    'A','Lowest model level H2O vapor mixing ratio',phys_decomp)
!     call addfld ('H216OVBT','kg/kg   ',1,    'A','Lowest model level H216O vapor mixing ratio',phys_decomp)
!     call addfld ('HDOVBT  ','kg/kg   ',1,    'A','Lowest model level HD16O vapor mixing ratio',phys_decomp)
!     call addfld ('H218OVBT','kg/kg   ',1,    'A','Lowest model level H218O vapor mixing ratio',phys_decomp)
   end if

   call addfld ('ATMEINT  ','J/m2    ',1, 'A','Vertically integrated total atmospheric energy ',phys_decomp)

   call addfld ('T1000      ','K     ',1,   'A','Temperature at 1000 mbar pressure surface',phys_decomp)
   call addfld ('T925       ','K     ',1,   'A','Temperature at 925 mbar pressure surface',phys_decomp)
   call addfld ('T700       ','K     ',1,   'A','Temperature at 700 mbar pressure surface',phys_decomp)
   call addfld ('T010       ','K     ',1,   'A','Temperature at 10 mbar pressure surface',phys_decomp)
   call addfld ('Q1000      ','kg/kg ',1,   'A','Specific Humidity at 1000 mbar pressure surface',phys_decomp)
   call addfld ('Q925       ','kg/kg ',1,   'A','Specific Humidity at 925 mbar pressure surface',phys_decomp)

   call addfld ('T7001000   ','K     ',1,   'A','Temperature difference 700 mb - 1000 mb',phys_decomp)
   call addfld ('TH7001000  ','K     ',1,   'A','Theta difference 700 mb - 1000 mb',phys_decomp)
   call addfld ('THE7001000 ','K     ',1,   'A','ThetaE difference 700 mb - 1000 mb',phys_decomp)

   call addfld ('T8501000   ','K     ',1,   'A','Temperature difference 850 mb - 1000 mb',phys_decomp)
   call addfld ('TH8501000  ','K     ',1,   'A','Theta difference 850 mb - 1000 mb',phys_decomp)
   call addfld ('THE8501000 ','K     ',1,   'A','ThetaE difference 850 mb - 1000 mb',phys_decomp)
   call addfld ('T9251000   ','K     ',1,   'A','Temperature difference 925 mb - 1000 mb',phys_decomp)
   call addfld ('TH9251000  ','K     ',1,   'A','Theta difference 925 mb - 1000 mb',phys_decomp)
   call addfld ('THE9251000 ','K     ',1,   'A','ThetaE difference 925 mb - 1000 mb',phys_decomp)

   ! This field is added by radiation when full physics is used
   if ( ideal_phys )then
      call addfld('QRS     ', 'K/s     ', pver, 'A', 'Solar heating rate', phys_decomp)
   end if

   ! ----------------------------
   ! determine default variables
   ! ----------------------------
   call phys_getopts(history_amwg_out   = history_amwg    , &
                     history_vdiag_out  = history_vdiag   , &
                     history_eddy_out   = history_eddy    , &
                     history_budget_out = history_budget  , &
                     history_budget_histfile_num_out = history_budget_histfile_num, &
                     history_waccm_out  = history_waccm)

   if (history_amwg) then
      call add_default ('PHIS    '  , 1, ' ')
      call add_default ('PS      '  , 1, ' ')
      call add_default ('T       '  , 1, ' ')
      call add_default ('U       '  , 1, ' ')
      call add_default ('V       '  , 1, ' ')
      call add_default (cnst_name(1), 1, ' ')
      call add_default ('Z3      '  , 1, ' ')
      call add_default ('OMEGA   '  , 1, ' ')
      call add_default ('VT      ', 1, ' ')
      call add_default ('VU      ', 1, ' ')
      call add_default ('VV      ', 1, ' ')
      call add_default ('VQ      ', 1, ' ')
      call add_default ('UU      ', 1, ' ')
      call add_default ('OMEGAT  ', 1, ' ')
      call add_default ('TMQ     ', 1, ' ')
      call add_default ('IVT     ', 1, ' ')
      call add_default ('PSL     ', 1, ' ')
      if (moist_physics) then
         call add_default ('RELHUM  ', 1, ' ')
      end if
      ! This field is added by radiation when full physics is used
      if ( ideal_phys )then
         call add_default('QRS     ', 1, ' ')
      end if
   end if

   if (history_vdiag) then
     call add_default ('U200', 2, ' ')
     call add_default ('V200', 2, ' ')
     call add_default ('U925', 2, ' ')
     call add_default ('U850', 2, ' ')
     call add_default ('U200', 3, ' ')
     call add_default ('U850', 3, ' ')
     call add_default ('OMEGA500', 3, ' ')
   end if

   if (history_eddy) then
      call add_default ('VT      ', 1, ' ')
      call add_default ('VU      ', 1, ' ')
      call add_default ('VV      ', 1, ' ')
      call add_default ('VQ      ', 1, ' ')
      call add_default ('UU      ', 1, ' ')
      call add_default ('OMEGAT  ', 1, ' ')
      call add_default ('OMEGAQ  ', 1, ' ')
      call add_default ('OMEGAU  ', 1, ' ')
      call add_default ('OMEGAV  ', 1, ' ')
   endif

   if ( history_budget ) then
      call add_default ('PHIS    '  , history_budget_histfile_num, ' ')
      call add_default ('PS      '  , history_budget_histfile_num, ' ')
      call add_default ('T       '  , history_budget_histfile_num, ' ')
      call add_default ('U       '  , history_budget_histfile_num, ' ')
      call add_default ('V       '  , history_budget_histfile_num, ' ')
      call add_default (cnst_name(1), history_budget_histfile_num, ' ')
      call add_default ('TTEND_TOT' , history_budget_histfile_num, ' ')

      ! State before physics (FV)
      call add_default ('TBP     '  , history_budget_histfile_num, ' ')
      call add_default (bpcnst(1)   , history_budget_histfile_num, ' ')
      ! State after physics (FV)
      call add_default ('TAP     '  , history_budget_histfile_num, ' ')
      call add_default ('UAP     '  , history_budget_histfile_num, ' ')
      call add_default ('VAP     '  , history_budget_histfile_num, ' ')
      call add_default (apcnst(1)   , history_budget_histfile_num, ' ')
      if ( dycore_is('LR') ) then
         call add_default ('TFIX    '    , history_budget_histfile_num, ' ')
         call add_default ('PTTEND_RESID', history_budget_histfile_num, ' ')
      end if
   end if

   ! create history variables for fourier coefficients of the diurnal
   ! and semidiurnal tide in T, U, V, and Z3
   call tidal_diag_init()

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Exit here for adiabatic/ideal physics cases !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (.not. moist_physics) return


   call addfld ('PDELDRY ','Pa      ',pver, 'A','Dry pressure difference between levels',phys_decomp)
   call addfld ('PSDRY   ','Pa      ',1,    'A','Surface pressure',phys_decomp)

   if (history_waccm) then
      call add_default ('PS      ', 2, ' ')
      call add_default ('T       ', 2, ' ')
   end if

   ! outfld calls in diag_conv

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call addfld ('DTCOND  ','K/s     ',pver, 'A','T tendency - moist processes',phys_decomp)
   call addfld ('DTCOND_24_COS','K/s',pver, 'A','T tendency - moist processes 24hr. cos coeff.',phys_decomp)
   call addfld ('DTCOND_24_SIN','K/s',pver, 'A','T tendency - moist processes 24hr. sin coeff.',phys_decomp)
   call addfld ('DTCOND_12_COS','K/s',pver, 'A','T tendency - moist processes 12hr. cos coeff.',phys_decomp)
   call addfld ('DTCOND_12_SIN','K/s',pver, 'A','T tendency - moist processes 12hr. sin coeff.',phys_decomp)

   ! determine number of constituents for which convective tendencies must be computed
   call cam_diag_init_helpers_select_impl()
   if (.not. cam_diag_init_helpers_use_native_impl) then
      conv_tend_code = 0
      if (diag_cnst_conv_tend == 'q_only') conv_tend_code = 1
      if (diag_cnst_conv_tend == 'all')    conv_tend_code = 2
      call cam_diag_init_helpers_log_entered()
      call cam_diag_log_direct(diag_init_logged, &
           'diag_init direct = codon; same-routine dqcond policy direct; addfld/add_default/pbuf/history native CAM API islands')
      dqcond_num = int(diag_init_codon( &
           merge(1_c_int64_t, 0_c_int64_t, history_budget), int(conv_tend_code, c_int64_t), int(pcnst, c_int64_t)))
   else
      if (history_budget) then
         dqcond_num = pcnst
      else
         if (diag_cnst_conv_tend == 'none')   dqcond_num = 0
         if (diag_cnst_conv_tend == 'q_only') dqcond_num = 1
         if (diag_cnst_conv_tend == 'all')    dqcond_num = pcnst
      end if
   end if

   do m = 1, dqcond_num
      dcconnam(m) = 'DC'//cnst_name(m)
   end do

   if (diag_cnst_conv_tend == 'q_only' .or. diag_cnst_conv_tend == 'all' .or. history_budget) then
      call addfld (dcconnam(1), 'kg/kg/s',pver,'A',trim(cnst_name(1))//' tendency due to moist processes',phys_decomp)
      if ( diag_cnst_conv_tend == 'q_only' .or. diag_cnst_conv_tend == 'all' ) then
         call add_default (dcconnam(1),                           1, ' ')
      end if
      if( history_budget ) then
         call add_default (dcconnam(1), history_budget_histfile_num, ' ')
      end if
      if (diag_cnst_conv_tend == 'all' .or. history_budget) then
         do m = 2, pcnst
            call addfld (dcconnam(m), 'kg/kg/s',pver,'A',trim(cnst_name(m))//' tendency due to moist processes',phys_decomp)
            if( diag_cnst_conv_tend == 'all' ) then
               call add_default (dcconnam(m),                           1, ' ')
            end if
            if( history_budget .and. (m == ixcldliq .or. m == ixcldice) ) then
               call add_default (dcconnam(m), history_budget_histfile_num, ' ')
            end if
         end do
      end if
   end if

   call addfld ('PRECL   ','m/s     ',1,    'A','Large-scale (stable) precipitation rate (liq + ice)'                ,phys_decomp)
   call addfld ('PRECC   ','m/s     ',1,    'A','Convective precipitation rate (liq + ice)'                          ,phys_decomp)
   call addfld ('PRECT   ','m/s     ',1,    'A','Total (convective and large-scale) precipitation rate (liq + ice)'  ,phys_decomp)
   call addfld ('PREC_PCW','m/s     ',1,    'A','LS_pcw precipitation rate',phys_decomp)
   call addfld ('PREC_zmc','m/s     ',1,    'A','CV_zmc precipitation rate',phys_decomp)
   call addfld ('PRECTMX ','m/s     ',1,    'X','Maximum (convective and large-scale) precipitation rate (liq+ice)'  ,phys_decomp)
   call addfld ('PRECSL  ','m/s     ',1,    'A','Large-scale (stable) snow rate (water equivalent)'                  ,phys_decomp)
   call addfld ('PRECSC  ','m/s     ',1,    'A','Convective snow rate (water equivalent)'                            ,phys_decomp)
   call addfld ('PRECCav ','m/s     ',1,    'A','Average large-scale precipitation (liq + ice)'                      ,phys_decomp)
   call addfld ('PRECLav ','m/s     ',1,    'A','Average convective precipitation  (liq + ice)'                      ,phys_decomp)

   !**********************
   !water tracers/isotopes
   !**********************
    if(trace_water) then
      do m=1,wtrc_nwset
        call addfld ('PRECT_'//trim(wtrc_out_names(m)),'m/s     ',1,'A', &
                     'Total (convective and large-scale) precipitation rate (liq + ice) for '//trim(wtrc_out_names(m)) ,phys_decomp)
      end do
    end if
   !**********************

   ! outfld calls in diag_surf

   call addfld ('SHFLX   ','W/m2    ',1,    'A','Surface sensible heat flux',phys_decomp)
   call addfld ('LHFLX   ','W/m2    ',1,    'A','Surface latent heat flux',phys_decomp)
   call addfld ('QFLX    ','kg/m2/s ',1,    'A','Surface water flux',phys_decomp)

   !**********************
   !water tracers/isotopes
   !**********************
   if(trace_water) then !Are water tracers on?
     do m=1,wtrc_nwset
       call addfld ('QFLX_'//trim(wtrc_out_names(m)),'kg/m2/s ',1,   'A','Surface water flux for '//trim(wtrc_out_names(m)),phys_decomp)
     end do

     call addfld ('buckH', 'm ',1,    'A','Bucket depth for bulk water',phys_decomp)
     call addfld ('buckD', 'm ',1,    'A','Bucket depth for HDO',phys_decomp)
     call addfld ('buck16','m ',1,    'A','Bucket depth for H2O16',phys_decomp)
     call addfld ('buck18','m ',1,    'A','Bucket depth for H2O18',phys_decomp)
   end if
   !**********************

   call addfld ('TAUX    ','N/m2    ',1,    'A','Zonal surface stress',phys_decomp)
   call addfld ('TAUY    ','N/m2    ',1,    'A','Meridional surface stress',phys_decomp)
   call addfld ('TREFHT  ','K       ',1,    'A','Reference height temperature',phys_decomp)
   call addfld ('TREFHTMN','K       ',1,    'M','Minimum reference height temperature over output period',phys_decomp)
   call addfld ('TREFHTMX','K       ',1,    'X','Maximum reference height temperature over output period',phys_decomp)
   call addfld ('QREFHT  ','kg/kg   ',1,    'A','Reference height humidity',phys_decomp)
   call addfld ('U10     ','m/s     ',1,    'A','10m wind speed',phys_decomp)
   call addfld ('RHREFHT ','fraction',1,    'A','Reference height relative humidity',phys_decomp)

   call addfld ('LANDFRAC','fraction',1,    'A','Fraction of sfc area covered by land',phys_decomp)
   call addfld ('ICEFRAC ','fraction',1,    'A','Fraction of sfc area covered by sea-ice',phys_decomp)
   call addfld ('OCNFRAC ','fraction',1,    'A','Fraction of sfc area covered by ocean',phys_decomp)

   call addfld ('TREFMNAV','K       ',1,    'A','Average of TREFHT daily minimum',phys_decomp)
   call addfld ('TREFMXAV','K       ',1,    'A','Average of TREFHT daily maximum',phys_decomp)

   call addfld ('TS      ','K       ',1,    'A','Surface temperature (radiative)',phys_decomp)
   call addfld ('TSMN    ','K       ',1,    'M','Minimum surface temperature over output period',phys_decomp)
   call addfld ('TSMX    ','K       ',1,    'X','Maximum surface temperature over output period',phys_decomp)
   call addfld ('SNOWHLND','m       ',1,    'A','Water equivalent snow depth',phys_decomp)
   call addfld ('SNOWHICE','m       ',1,    'A','Snow depth over ice',phys_decomp, fill_value = 1.e30_r8)
   ! nanr
   call addfld ('TBOT    ','K       ',1,    'A','Lowest model level temperature', phys_decomp)

   call addfld ('ASDIR',   '1',       1,    'A','albedo: shortwave, direct', phys_decomp)
   call addfld ('ASDIF',   '1',       1,    'A','albedo: shortwave, diffuse', phys_decomp)
   call addfld ('ALDIR',   '1',       1,    'A','albedo: longwave, direct', phys_decomp)
   call addfld ('ALDIF',   '1',       1,    'A','albedo: longwave, diffuse', phys_decomp)
   call addfld ('SST',     'K',       1,    'A','sea surface temperature', phys_decomp)

   ! defaults
   if (history_amwg) then
       call add_default ('DTCOND  ', 1, ' ')
       call add_default ('PRECL   ', 1, ' ')
       call add_default ('PRECC   ', 1, ' ')
       call add_default ('PRECSL  ', 1, ' ')
       call add_default ('PRECSC  ', 1, ' ')
       call add_default ('SHFLX   ', 1, ' ')
       call add_default ('LHFLX   ', 1, ' ')
       call add_default ('QFLX    ', 1, ' ')
       call add_default ('TAUX    ', 1, ' ')
       call add_default ('TAUY    ', 1, ' ')
       call add_default ('TREFHT  ', 1, ' ')
       call add_default ('LANDFRAC', 1, ' ')
       call add_default ('OCNFRAC ', 1, ' ')
       call add_default ('QREFHT  ', 1, ' ')
       call add_default ('U10     ', 1, ' ')
       call add_default ('ICEFRAC ', 1, ' ')
       call add_default ('TS      ', 1, ' ')
       call add_default ('TSMN    ', 1, ' ')
       call add_default ('TSMX    ', 1, ' ')
       call add_default ('SNOWHLND', 1, ' ')
       call add_default ('SNOWHICE', 1, ' ')
    endif

    if (history_vdiag) then
        call add_default ('PRECT   ', 2, ' ')
        call add_default ('PRECT   ', 3, ' ')
        call add_default ('PRECT   ', 4, ' ')
    end if

   ! outfld calls in diag_phys_tend_writeout

   call addfld ('PTTEND  '   ,'K/s     ',pver, 'A','T total physics tendency'                             ,phys_decomp)
   call addfld (ptendnam(       1),  'kg/kg/s ',pver, 'A',trim(cnst_name(       1))//' total physics tendency '      ,phys_decomp)
   call addfld (ptendnam(ixcldliq),  'kg/kg/s ',pver, 'A',trim(cnst_name(ixcldliq))//' total physics tendency '      ,phys_decomp)
   call addfld (ptendnam(ixcldice),  'kg/kg/s ',pver, 'A',trim(cnst_name(ixcldice))//' total physics tendency '      ,phys_decomp)
   if ( dycore_is('LR') )then
      call addfld (dmetendnam(       1),'kg/kg/s ',pver, 'A', &
           trim(cnst_name(       1))//' dme adjustment tendency (FV) ',phys_decomp)
      call addfld (dmetendnam(ixcldliq),'kg/kg/s ',pver, 'A', &
           trim(cnst_name(ixcldliq))//' dme adjustment tendency (FV) ',phys_decomp)
      call addfld (dmetendnam(ixcldice),'kg/kg/s ',pver, 'A', &
           trim(cnst_name(ixcldice))//' dme adjustment tendency (FV) ',phys_decomp)
   end if

   if ( history_budget ) then
      call add_default ('PTTEND'          , history_budget_histfile_num, ' ')
      call add_default (ptendnam(       1), history_budget_histfile_num, ' ')
      call add_default (ptendnam(ixcldliq), history_budget_histfile_num, ' ')
      call add_default (ptendnam(ixcldice), history_budget_histfile_num, ' ')
      if ( dycore_is('LR') )then
         call add_default(dmetendnam(1)       , history_budget_histfile_num, ' ')
         call add_default(dmetendnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(dmetendnam(ixcldice), history_budget_histfile_num, ' ')
      end if
      if( history_budget_histfile_num > 1 ) then
         call add_default ('DTCOND  '         , history_budget_histfile_num, ' ')
      end if
   end if

   ! outfld calls in diag_physvar_ic

   call addfld ('QCWAT&IC   ','kg/kg   ',pver, 'I','q associated with cloud water'                   ,phys_decomp)
   call addfld ('TCWAT&IC   ','kg/kg   ',pver, 'I','T associated with cloud water'                   ,phys_decomp)
   call addfld ('LCWAT&IC   ','kg/kg   ',pver, 'I','Cloud water (ice + liq'                          ,phys_decomp)
   call addfld ('CLOUD&IC   ','fraction',pver, 'I','Cloud fraction'                                  ,phys_decomp)
   call addfld ('CONCLD&IC   ','fraction',pver, 'I','Convective cloud fraction'                      ,phys_decomp)
   call addfld ('TKE&IC     ','m2/s2   ',pverp,'I','Turbulent Kinetic Energy'                        ,phys_decomp)
   call addfld ('CUSH&IC    ','m       ',1,    'I','Convective Scale Height'                         ,phys_decomp)
   call addfld ('KVH&IC     ','m2/s    ',pverp,'I','Vertical diffusion diffusivities (heat/moisture)',phys_decomp)
   call addfld ('KVM&IC     ','m2/s    ',pverp,'I','Vertical diffusion diffusivities (momentum)'     ,phys_decomp)
   call addfld ('PBLH&IC    ','m       ',1,    'I','PBL height'                                      ,phys_decomp)
   call addfld ('TPERT&IC   ','K       ',1,    'I','Perturbation temperature (eddies in PBL)'        ,phys_decomp)
   call addfld ('QPERT&IC   ','kg/kg   ',1,    'I','Perturbation specific humidity (eddies in PBL)'  ,phys_decomp)
   ! nanr
   call addfld ('TBOT&IC    ','K       ',1,    'I','Lowest model level temperature '                 ,phys_decomp)


   ! Initial file - Optional fields

   if (inithist_all) then
      call add_default ('CONCLD&IC  ',0, 'I')
      call add_default ('QCWAT&IC   ',0, 'I')
      call add_default ('TCWAT&IC   ',0, 'I')
      call add_default ('LCWAT&IC   ',0, 'I')
      call add_default ('PBLH&IC    ',0, 'I')
      call add_default ('TPERT&IC   ',0, 'I')
      call add_default ('QPERT&IC   ',0, 'I')
      call add_default ('CLOUD&IC   ',0, 'I')
      call add_default ('TKE&IC     ',0, 'I')
      call add_default ('CUSH&IC    ',0, 'I')
      call add_default ('KVH&IC     ',0, 'I')
      call add_default ('KVM&IC     ',0, 'I')
      ! nanr
      call add_default ('TBOT&IC    ',0, 'I')
   end if

   ! CAM export state
   call addfld('a2x_BCPHIWET', 'kg/m2/s', 1, 'A', 'wetdep of hydrophilic black carbon',   phys_decomp)
   call addfld('a2x_BCPHIDRY', 'kg/m2/s', 1, 'A', 'drydep of hydrophilic black carbon',   phys_decomp)
   call addfld('a2x_BCPHODRY', 'kg/m2/s', 1, 'A', 'drydep of hydrophobic black carbon',   phys_decomp)
   call addfld('a2x_OCPHIWET', 'kg/m2/s', 1, 'A', 'wetdep of hydrophilic organic carbon', phys_decomp)
   call addfld('a2x_OCPHIDRY', 'kg/m2/s', 1, 'A', 'drydep of hydrophilic organic carbon', phys_decomp)
   call addfld('a2x_OCPHODRY', 'kg/m2/s', 1, 'A', 'drydep of hydrophobic organic carbon', phys_decomp)
   call addfld('a2x_DSTWET1',  'kg/m2/s', 1, 'A', 'wetdep of dust (bin1)',                phys_decomp)
   call addfld('a2x_DSTDRY1',  'kg/m2/s', 1, 'A', 'drydep of dust (bin1)',                phys_decomp)
   call addfld('a2x_DSTWET2',  'kg/m2/s', 1, 'A', 'wetdep of dust (bin2)',                phys_decomp)
   call addfld('a2x_DSTDRY2',  'kg/m2/s', 1, 'A', 'drydep of dust (bin2)',                phys_decomp)
   call addfld('a2x_DSTWET3',  'kg/m2/s', 1, 'A', 'wetdep of dust (bin3)',                phys_decomp)
   call addfld('a2x_DSTDRY3',  'kg/m2/s', 1, 'A', 'drydep of dust (bin3)',                phys_decomp)
   call addfld('a2x_DSTWET4',  'kg/m2/s', 1, 'A', 'wetdep of dust (bin4)',                phys_decomp)
   call addfld('a2x_DSTDRY4',  'kg/m2/s', 1, 'A', 'drydep of dust (bin4)',                phys_decomp)

  qcwat_idx  = pbuf_get_index('QCWAT',ierr)
  tcwat_idx  = pbuf_get_index('TCWAT',ierr)
  lcwat_idx  = pbuf_get_index('LCWAT',ierr)
  cld_idx    = pbuf_get_index('CLD')
  concld_idx = pbuf_get_index('CONCLD')

  tke_idx  = pbuf_get_index('tke')
  kvm_idx  = pbuf_get_index('kvm')
  kvh_idx  = pbuf_get_index('kvh')
  cush_idx = pbuf_get_index('cush')

  pblh_idx  = pbuf_get_index('pblh')
  tpert_idx = pbuf_get_index('tpert')
  qpert_idx = pbuf_get_index('qpert',ierr)

  prec_dp_idx  = pbuf_get_index('PREC_DP')
  snow_dp_idx  = pbuf_get_index('SNOW_DP')
  prec_sh_idx  = pbuf_get_index('PREC_SH')
  snow_sh_idx  = pbuf_get_index('SNOW_SH')
  prec_sed_idx = pbuf_get_index('PREC_SED')
  snow_sed_idx = pbuf_get_index('SNOW_SED')
  prec_pcw_idx = pbuf_get_index('PREC_PCW')
  snow_pcw_idx = pbuf_get_index('SNOW_PCW')

end subroutine diag_init

!===============================================================================

subroutine diag_allocate()
   use iso_c_binding, only: c_int64_t
   use infnan, only: nan, assignment(=)

   ! Allocate memory for module variables.
   ! Done at the begining of a physics step at same point as the pbuf allocate for
   ! variables with "physpkg" scope.

   ! Local variables
   character(len=*), parameter :: sub = 'diag_allocate'
   integer :: i, istat

   call cam_diag_touch_and_log(3_c_int64_t, diag_allocate_logged, &
        'diag_allocate direct = codon; allocation selector/touch direct = codon; native allocate/nan-fill island')

   allocate(dtcond(pcols,pver,begchunk:endchunk), stat=istat)
   if ( istat /= 0 ) call endrun (sub//': ERROR: allocate failed')
   dtcond = nan

   if (dqcond_num > 0) then
      allocate(dqcond(dqcond_num))
      do i = 1, dqcond_num
         allocate(dqcond(i)%cnst(pcols,pver,begchunk:endchunk), stat=istat)
         if ( istat /= 0 ) call endrun (sub//': ERROR: allocate failed')
         dqcond(i)%cnst = nan
      end do
   end if

end subroutine diag_allocate

!===============================================================================

subroutine diag_deallocate()

   use iso_c_binding, only: c_int64_t

! Deallocate memory for module variables.
! Done at the end of a physics step at same point as the pbuf deallocate for
! variables with "physpkg" scope.

! Local variables
   character(len=*), parameter :: sub = 'diag_deallocate'
   integer :: i, istat

   call cam_diag_touch_and_log(4_c_int64_t, diag_deallocate_logged, &
        'diag_deallocate direct = codon; deallocation selector/touch direct = codon; native deallocate/error island')

   deallocate(dtcond, stat=istat)
   if ( istat /= 0 ) call endrun (sub//': ERROR: deallocate failed')

   if (dqcond_num > 0) then
      do i = 1, dqcond_num
         deallocate(dqcond(i)%cnst, stat=istat)
         if ( istat /= 0 ) call endrun (sub//': ERROR: deallocate failed')
      end do
      deallocate(dqcond, stat=istat)
      if ( istat /= 0 ) call endrun (sub//': ERROR: deallocate failed')
   end if

end subroutine diag_deallocate
!===============================================================================

subroutine cam_diag_init_helpers_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cam_diag_init_helpers_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CAM_DIAG_INIT_HELPERS_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      cam_diag_init_helpers_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      cam_diag_init_helpers_use_native_impl = .false.
   end if

   cam_diag_init_helpers_impl_selected = .true.

   if (masterproc) then
      if (cam_diag_init_helpers_use_native_impl) then
         write(iulog,*) 'cam_diag_init_helpers implementation = native'
      else
         write(iulog,*) 'cam_diag_init_helpers implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine cam_diag_init_helpers_select_impl

!===============================================================================

subroutine cam_diag_init_helpers_log_entered()

   if (cam_diag_init_helpers_entered_logged) return
   cam_diag_init_helpers_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'cam_diag_init_helpers entered (dqcond policy direct = codon)'
      call flush(iulog)
   end if

end subroutine cam_diag_init_helpers_log_entered

!===============================================================================

subroutine cam_diag_conv_batch_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line
   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('CAM_DIAG_CONV_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine cam_diag_conv_batch_append_proof

subroutine cam_diag_conv_batch_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (cam_diag_conv_batch_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('CAM_DIAG_CONV_BATCH_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      cam_diag_conv_batch_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      cam_diag_conv_batch_use_native_impl = .false.
   end if

   cam_diag_conv_batch_impl_selected = .true.

   if (masterproc) then
      if (cam_diag_conv_batch_use_native_impl) then
         write(iulog,*) 'cam_diag_conv_batch implementation = native'
         call cam_diag_conv_batch_append_proof('cam_diag_conv_batch selector entered implementation = native')
      else
         write(iulog,*) 'cam_diag_conv_batch implementation = codon'
         call cam_diag_conv_batch_append_proof('cam_diag_conv_batch selector entered implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine cam_diag_conv_batch_select_impl

subroutine cam_diag_conv_batch_log_tend_ini_entered()

   if (cam_diag_conv_tend_ini_entered_logged) return
   cam_diag_conv_tend_ini_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'diag_conv_tend_ini direct = codon; state/dqcond/T_TTEND copy body direct; pbuf_get_field native CAM API island'
      call cam_diag_conv_batch_append_proof( &
           'diag_conv_tend_ini direct = codon; state/dqcond/T_TTEND copy body direct; pbuf_get_field native CAM API island')
      call flush(iulog)
   end if

end subroutine cam_diag_conv_batch_log_tend_ini_entered

subroutine cam_diag_conv_batch_log_diag_conv_entered()

   if (cam_diag_conv_entered_logged) return
   cam_diag_conv_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'cam_diag_conv_batch entered (unified stage-dispatch update direct = codon)'
      call cam_diag_conv_batch_append_proof('cam_diag_conv_batch entered (unified stage-dispatch update direct = codon)')
      call flush(iulog)
   end if

end subroutine cam_diag_conv_batch_log_diag_conv_entered

!===============================================================================

subroutine cam_diag_conv_batch_log_precip_dtcond_entered()

   if (cam_diag_conv_precip_dtcond_entered_logged) return
   cam_diag_conv_precip_dtcond_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'cam_diag_conv_batch entered (unified stage-dispatch precip totals + dtcond direct = codon)'
      call cam_diag_conv_batch_append_proof('cam_diag_conv_batch entered (unified stage-dispatch precip totals + dtcond direct = codon)')
      call flush(iulog)
   end if

end subroutine cam_diag_conv_batch_log_precip_dtcond_entered

!===============================================================================

subroutine diag_phys_tend_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line
   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('DIAG_PHYS_TEND_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine diag_phys_tend_append_proof

subroutine diag_phys_tend_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (diag_phys_tend_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('DIAG_PHYS_TEND_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      diag_phys_tend_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      diag_phys_tend_use_native_impl = .false.
   end if

   diag_phys_tend_impl_selected = .true.

   if (masterproc) then
      if (diag_phys_tend_use_native_impl) then
         write(iulog,*) 'diag_phys_tend implementation = native'
         call diag_phys_tend_append_proof('diag_phys_tend selector entered implementation = native')
      else
         write(iulog,*) 'diag_phys_tend implementation = codon'
         call diag_phys_tend_append_proof('diag_phys_tend selector entered implementation = codon')
      end if
      call flush(iulog)
   end if

end subroutine diag_phys_tend_select_impl

subroutine diag_phys_tend_log_entered()

   if (diag_phys_tend_entered_logged) return
   diag_phys_tend_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'diag_phys_tend_writeout direct = codon; numeric tendency-update stages direct; cnst_get_ind/outfld/pbuf/check_energy native CAM API islands'
      call diag_phys_tend_append_proof( &
           'diag_phys_tend_writeout direct = codon; numeric tendency-update stages direct; cnst_get_ind/outfld/pbuf/check_energy native CAM API islands')
      call flush(iulog)
   end if

end subroutine diag_phys_tend_log_entered

!===============================================================================

subroutine diag_conv_tend_ini(state,pbuf)

! Initialize convective tendency calcs.

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

! Argument:

   type(physics_state), target, intent(in) :: state

   type(physics_buffer_desc), pointer :: pbuf(:)

! Local variables:

   integer :: i, k, m, lchnk, ncol
   real(r8), pointer, dimension(:,:) :: t_ttend
   real(r8), target :: dqcond_work(pcols,pver)

   lchnk = state%lchnk
   ncol  = state%ncol

   call cam_diag_conv_batch_select_impl()

   if (cam_diag_conv_batch_use_native_impl) then
      do k = 1, pver
         do i = 1, ncol
            dtcond(i,k,lchnk) = state%s(i,k)
         end do
      end do

      do m = 1, dqcond_num
         do k = 1, pver
            do i = 1, ncol
               dqcond(m)%cnst(i,k,lchnk) = state%q(i,k,m)
            end do
         end do
      end do
   else
      call cam_diag_conv_batch_log_tend_ini_entered()

      call diag_conv_tend_ini_codon( &
           1_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           int(pcnst, c_int64_t), 0_c_int64_t, c_loc(state%s), c_loc(dtcond(1,1,lchnk)) &
      )

      do m = 1, dqcond_num
         call diag_conv_tend_ini_codon( &
              2_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
              int(pcnst, c_int64_t), int(m, c_int64_t), c_loc(state%q), c_loc(dqcond_work) &
         )
         dqcond(m)%cnst(:,:,lchnk) = dqcond_work(:,:)
      end do
   end if

   !! initialize to pbuf T_TTEND to temperature at first timestep
   if (is_first_step()) then
      do m = 1, dyn_time_lvls
         call pbuf_get_field(pbuf, t_ttend_idx, t_ttend, start=(/1,1,m/), kount=(/pcols,pver,1/))
         if (cam_diag_conv_batch_use_native_impl) then
            t_ttend(:ncol,:) = state%t(:ncol,:)
         else
            call diag_conv_tend_ini_codon( &
                 3_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
                 int(pcnst, c_int64_t), 0_c_int64_t, c_loc(state%t), c_loc(t_ttend) &
            )
         end if
      end do
   end if

end subroutine diag_conv_tend_ini
!===============================================================================

  subroutine diag_phys_writeout(state, psl)

!-----------------------------------------------------------------------
!
! Purpose: record dynamics variables on physics grid
!
!-----------------------------------------------------------------------
    use physconst,          only: gravit, rga, rair, cpair, latvap, rearth, pi, cappa
    use time_manager,       only: get_nstep
    use interpolate_data,   only: vertinterp
    use constituent_burden, only: constituent_burden_comp
    use cam_control_mod,    only: moist_physics
    use co2_cycle,          only: c_i, co2_transport
    use iso_c_binding,      only: c_int64_t

    use tidal_diag,         only: tidal_diag_write
!-----------------------------------------------------------------------
!
! Arguments
!
   type(physics_state), intent(inout) :: state
   real(r8), optional , intent(out)   :: psl(pcols)
!
!---------------------------Local workspace-----------------------------
!
    real(r8) ftem(pcols,pver) ! temporary workspace
    real(r8) ftem1(pcols,pver) ! another temporary workspace
    real(r8) ftem2(pcols,pver) ! another temporary workspace
    real(r8) ftem4(pcols,pver) ! another temporary workspace
    real(r8) ftem5(pcols,pver) ! another temporary workspace
    real(r8) psl_tmp(pcols)   ! Sea Level Pressure
    real(r8) z3(pcols,pver)   ! geo-potential height
    real(r8) p_surf(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_t1(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_t2(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_q1(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_q2(pcols)    ! data interpolated to a pressure surface
    real(r8) tem2(pcols,pver) ! temporary workspace
    real(r8) timestep(pcols)  ! used for outfld call
    real(r8) esl(pcols,pver)   ! saturation vapor pressures
    real(r8) esi(pcols,pver)   !
    real(r8) dlon(pcols)      ! width of grid cell (meters)
    integer  plon             ! number of longitudes

    integer i, k, m, lchnk, ncol, nstep
    integer(c_int64_t) :: direct_stage_c
!
!-----------------------------------------------------------------------
!
    lchnk = state%lchnk
    ncol  = state%ncol

    call diag_phys_writeout_select_impl()
    if (.not. diag_phys_writeout_use_native_impl) then
       direct_stage_c = diag_phys_writeout_codon(7_c_int64_t)
       if (direct_stage_c == 7_c_int64_t) then
          call diag_phys_writeout_batch_log_entered()
          call cam_diag_log_direct(diag_phys_writeout_logged, &
               'diag_phys_writeout same-routine direct = codon; numeric field helpers direct = codon; ' // &
               'outfld/vertinterp/tidal/history native CAM API islands')
       end if
    end if

    ! Output NSTEP for debugging
    nstep = get_nstep()
    timestep(:ncol) = nstep
    call outfld ('NSTEP   ',timestep, pcols, lchnk)

    call outfld('T       ',state%t , pcols   ,lchnk   )
    call outfld('PS      ',state%ps, pcols   ,lchnk   )
    call outfld('U       ',state%u , pcols   ,lchnk   )
    call outfld('V       ',state%v , pcols   ,lchnk   )
    do m=1,pcnst
       if ( cnst_cam_outfld(m) ) then
          call outfld(cnst_name(m),state%q(1,1,m),pcols ,lchnk )
       end if
    end do

    if (co2_transport()) then
       do m = 1,4
          call outfld(trim(cnst_name(c_i(m)))//'_BOT', state%q(1,pver,c_i(m)), pcols, lchnk)
       end do
    end if

    ! column burdens of all constituents except water vapor
    call constituent_burden_comp(state)

    if ( moist_physics) then
       call outfld('PDELDRY ',state%pdeldry, pcols, lchnk)
       call outfld('PSDRY',   state%psdry,   pcols, lchnk)
    end if

    call outfld('PHIS    ',state%phis,    pcols,   lchnk     )



#if (defined BFB_CAM_SCAM_IOP )
    call outfld('phis    ',state%phis,    pcols,   lchnk     )
#endif

!
! Add height of surface to midpoint height above surface
!
    call diag_phys_writeout_fill_z3(ncol, state%zm, state%phis, z3)
    call outfld('Z3      ',z3,pcols,lchnk)
!
! Output Z3 on pressure surfaces
!
    if (hist_fld_active('Z1000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 100000._r8, z3, p_surf)
       call outfld('Z1000    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z700')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, z3, p_surf)
       call outfld('Z700    ', p_surf, pcols, lchnk)
    end if
    ! -- nanr
    if (hist_fld_active('Z850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, z3, p_surf)
       call outfld('Z850    ', p_surf, pcols, lchnk)
    end if
    ! -- end
    if (hist_fld_active('Z500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, z3, p_surf)
       call outfld('Z500    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z300')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, z3, p_surf)
       call outfld('Z300    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z250')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 25000._r8, z3, p_surf)
       call outfld('Z250    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, z3, p_surf)
       call outfld('Z200    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z100')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 10000._r8, z3, p_surf)
       call outfld('Z100    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z050')) then
       call vertinterp(ncol, pcols, pver, state%pmid,  5000._r8, z3, p_surf)
       call outfld('Z050    ', p_surf, pcols, lchnk)
    end if
!
! Quadratic height fiels Z3*Z3
!
    call diag_phys_writeout_square_field(ncol, z3, ftem)
    call outfld('ZZ      ',ftem,pcols,lchnk)

    call diag_phys_writeout_mul_scalar_field(ncol, z3, state%v, gravit, ftem)
    call outfld('VZ      ',ftem,  pcols,lchnk)
!
! Meridional advection fields
!
    call diag_phys_writeout_mul_field(ncol, state%v, state%t, ftem)
    call outfld ('VT      ',ftem    ,pcols   ,lchnk     )

! -- nanr
    call diag_phys_writeout_mul_field(ncol, state%u, state%t, ftem)
    call outfld ('UT      ',ftem    ,pcols   ,lchnk     )

    call diag_phys_writeout_mul_field(ncol, state%u, state%q(:,:,1), ftem)
    call outfld ('UQ      ',ftem    ,pcols   ,lchnk     )
! -- end

    call diag_phys_writeout_mul_field(ncol, state%v, state%q(:,:,1), ftem)
    call outfld ('VQ      ',ftem    ,pcols   ,lchnk     )

    call diag_phys_writeout_square_field(ncol, state%q(:,:,1), ftem)
    call outfld ('QQ      ',ftem    ,pcols   ,lchnk     )

    call diag_phys_writeout_square_field(ncol, state%v, ftem)
    call outfld ('VV      ',ftem    ,pcols   ,lchnk     )

    call diag_phys_writeout_mul_field(ncol, state%v, state%u, ftem)
    call outfld ('VU      ',ftem    ,pcols   ,lchnk     )

! zonal advection

    call diag_phys_writeout_square_field(ncol, state%u, ftem)
    call outfld ('UU      ',ftem    ,pcols   ,lchnk     )

! Wind speed
    call diag_phys_writeout_wspeed_field(ncol, state%u, state%v, ftem)
    call outfld ('WSPEED  ',ftem    ,pcols   ,lchnk     )
    call outfld ('WSPDSRFMX',ftem(:,pver)   ,pcols   ,lchnk     )
    call outfld ('WSPDSRFAV',ftem(:,pver)   ,pcols   ,lchnk     )

! Vertical velocity and advection

    if (single_column) then
       call outfld('OMEGA   ',wfld,    pcols,   lchnk     )
    else
       call outfld('OMEGA   ',state%omega,    pcols,   lchnk     )
    endif

#if (defined BFB_CAM_SCAM_IOP )
    call outfld('omega   ',state%omega,    pcols,   lchnk     )
#endif

    call diag_phys_writeout_mul_field(ncol, state%omega, state%t, ftem)
    call outfld('OMEGAT  ',ftem,    pcols,   lchnk     )
    call diag_phys_writeout_mul_field(ncol, state%omega, state%u, ftem)
    call outfld('OMEGAU  ',ftem,    pcols,   lchnk     )
    call diag_phys_writeout_mul_field(ncol, state%omega, state%v, ftem)
    call outfld('OMEGAV  ',ftem,    pcols,   lchnk     )
    call diag_phys_writeout_mul_field(ncol, state%omega, state%q(:,:,1), ftem)
    call outfld('OMEGAQ  ',ftem,    pcols,   lchnk     )
    call diag_phys_writeout_square_field(ncol, state%omega, ftem)
    call outfld('OMGAOMGA',ftem,    pcols,   lchnk     )
!
! Output omega at 850 and 500 mb pressure levels
!
    if (hist_fld_active('OMEGA850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%omega, p_surf)
       call outfld('OMEGA850', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('OMEGA500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%omega, p_surf)
       call outfld('OMEGA500', p_surf, pcols, lchnk)
    end if
! -- nanr
    if (hist_fld_active('OMEGA200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%omega, p_surf)
       call outfld('OMEGA200', p_surf, pcols, lchnk)
    end if
! -- end
!
! Mass of q, by layer and vertically integrated
!
    call diag_phys_writeout_mass_and_tmq(ncol, state%q(:,:,1), state%pdel, ftem, p_surf)
    call outfld ('MQ      ',ftem    ,pcols   ,lchnk     )
    call outfld ('TMQ     ',p_surf, pcols   ,lchnk     )

    !**********************
    !Water tracers/isotopes
    !**********************
     if(trace_water) then
       do m=1,wtrc_nwset
       !-----------
         call diag_phys_writeout_wtrc_column(ncol, state%q(:,:,wtrc_iatype(m,iwtvap)), state%v, state%pdel, &
              rga, 1, ftem)
         call outfld ('TMQ_'//trim(wtrc_out_names(m)), ftem, pcols, lchnk)
       !-----------
         call diag_phys_writeout_wtrc_column(ncol, state%q(:,:,wtrc_iatype(m,iwtvap)), state%v, state%pdel, &
              rga, 2, ftem)
         call outfld ('TVQ_'//trim(wtrc_out_names(m)), ftem, pcols, lchnk)
       !-----------
         call diag_phys_writeout_wtrc_column(ncol, state%q(:,:,wtrc_iatype(m,iwtvap)), state%u, state%pdel, &
              rga, 2, ftem)
         call outfld ('TUQ_'//trim(wtrc_out_names(m)), ftem, pcols, lchnk)
       end do
       !-----------
     end if
    !**********************

!CAS integrated vapor transport calculation

    call diag_phys_writeout_ivt(ncol, state%q(:,:,1), state%u, state%v, state%pdel, rga, ftem4, ftem5, ftem)

    call outfld ('IVT     ',ftem, pcols   ,lchnk     )

    !just output uq*dp/g
    call diag_phys_writeout_copy_col1(ncol, ftem4, ftem)
    call outfld ('uIVT     ',ftem, pcols   ,lchnk     )

    !just output vq*dp/g
    call diag_phys_writeout_copy_col1(ncol, ftem5, ftem)
    call outfld ('vIVT     ',ftem, pcols   ,lchnk     )

!CAS


    if (moist_physics) then

       ! Relative humidity
       if (hist_fld_active('RELHUM')) then
          call qsat(state%t(:ncol,:), state%pmid(:ncol,:), &
               tem2(:ncol,:), ftem(:ncol,:))
          call diag_phys_writeout_scale_relhum(ncol, state%q(:,:,1), ftem)
          call outfld ('RELHUM  ',ftem    ,pcols   ,lchnk     )
       end if

       if (hist_fld_active('RHW') .or. hist_fld_active('RHI') .or. hist_fld_active('RHCFMIP') ) then

          ! RH w.r.t liquid (water)
          call qsat_water (state%t(:ncol,:), state%pmid(:ncol,:), &
               esl(:ncol,:), ftem(:ncol,:))
          call diag_phys_writeout_scale_relhum(ncol, state%q(:,:,1), ftem)
          call outfld ('RHW  ',ftem    ,pcols   ,lchnk     )

          ! Convert to RHI (ice)
          do i=1,ncol
             do k=1,pver
                esi(i,k)=svp_ice(state%t(i,k))
             end do
          end do
          call diag_phys_writeout_rhi_rhcfmip(ncol, state%t, esl, esi, ftem, ftem1, ftem2)
          call outfld ('RHI  ',ftem1    ,pcols   ,lchnk     )

          call outfld ('RHCFMIP  ',ftem2    ,pcols   ,lchnk     )

       end if

    end if
!
! Sea level pressure
!
    if (present(psl) .or. hist_fld_active('PSL')) then
       call cpslec (ncol, state%pmid, state%phis, state%ps, state%t,psl_tmp, gravit, rair)
       call outfld ('PSL     ',psl_tmp  ,pcols, lchnk     )
       if (present(psl)) then
          psl(:ncol) = psl_tmp(:ncol)
       end if
    end if
!
! Output T,q,u,v fields on pressure surfaces
!
    if (hist_fld_active('T850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%t, p_surf)
       call outfld('T850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%t, p_surf)
       call outfld('T500    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T300')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, state%t, p_surf)
       call outfld('T300    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%t, p_surf)
       call outfld('T200    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('Q850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%q(1,1,1), p_surf)
       call outfld('Q850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('Q500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%q(1,1,1), p_surf)
       call outfld('Q500    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('Q200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%q(1,1,1), p_surf)
       call outfld('Q200    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U925')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 92500._r8, state%u, p_surf)
       call outfld('U925    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%u, p_surf)
       call outfld('U850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U250')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 25000._r8, state%u, p_surf)
       call outfld('U250    ', p_surf, pcols, lchnk )
    end if
! -- nanr
    if (hist_fld_active('U500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%u, p_surf)
       call outfld('U500    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U600')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 60000._r8, state%u, p_surf)
       call outfld('U600    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U700')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, state%u, p_surf)
       call outfld('U700    ', p_surf, pcols, lchnk )
    end if
! -- end
    if (hist_fld_active('U200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%u, p_surf)
       call outfld('U200    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U010')) then
       call vertinterp(ncol, pcols, pver, state%pmid,  1000._r8, state%u, p_surf)
       call outfld('U010    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V925')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 92500._r8, state%v, p_surf)
       call outfld('V925    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%v, p_surf)
       call outfld('V850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V250')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 25000._r8, state%v, p_surf)
       call outfld('V250    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%v, p_surf)
       call outfld('V200    ', p_surf, pcols, lchnk )
    end if
! -- nanr
    if (hist_fld_active('V500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%v, p_surf)
       call outfld('V500    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V600')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 60000._r8, state%v, p_surf)
       call outfld('V600    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V700')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, state%v, p_surf)
       call outfld('V700    ', p_surf, pcols, lchnk )
    end if
! -- end

    call diag_phys_writeout_square_field(ncol, state%t, ftem)
    call outfld('TT      ',ftem    ,pcols   ,lchnk   )
!
! Output U, V, T, Q, P and Z at bottom level
!
    call outfld ('UBOT    ', state%u(1,pver)  ,  pcols, lchnk)
    call outfld ('VBOT    ', state%v(1,pver)  ,  pcols, lchnk)
    call outfld ('QBOT    ', state%q(1,pver,1),  pcols, lchnk)
    call outfld ('ZBOT    ', state%zm(1,pver) , pcols, lchnk)

    if(trace_water) then !using water tracers or isotopes?
      do m=1, wtrc_nwset
        call outfld (trim(wtrc_srfvap_names(m))//'BT', state%q(1,pver,wtrc_iatype(m,iwtvap)), pcols, lchnk)
      end do

!      call outfld ('H2OVBT  ', state%q(1,pver,wtrc_iatype(1,iwtvap)),  pcols, lchnk)
!      call outfld ('H216OVBT', state%q(1,pver,wtrc_iatype(2,iwtvap)),  pcols, lchnk)
!      call outfld ('HDOVBT  ', state%q(1,pver,wtrc_iatype(3,iwtvap)),  pcols, lchnk)
!      call outfld ('H218OVBT', state%q(1,pver,wtrc_iatype(4,iwtvap)),  pcols, lchnk)
    end if

! Total energy of the atmospheric column for atmospheric heat storage calculations

    call diag_phys_writeout_atmeint(ncol, state%t, state%q(:,:,1), state%u, state%v, state%pdel, state%phis, timestep)
    call outfld ('ATMEINT   ',timestep  ,pcols   ,lchnk     )

!! Boundary layer atmospheric stability, temperature, water vapor diagnostics

    if (hist_fld_active('T1000')      .or. &
        hist_fld_active('T9251000')   .or. &
        hist_fld_active('TH9251000')  .or. &
        hist_fld_active('THE9251000') .or. &
        hist_fld_active('T8501000')   .or. &
        hist_fld_active('TH8501000')  .or. &
        hist_fld_active('THE8501000') .or. &
        hist_fld_active('T7001000')   .or. &
        hist_fld_active('TH7001000')  .or. &
        hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 100000._r8, state%t, p_surf_t1)
    end if

    if (hist_fld_active('T925')       .or. &
        hist_fld_active('T9251000')   .or. &
        hist_fld_active('TH9251000')  .or. &
        hist_fld_active('THE9251000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 92500._r8, state%t, p_surf_t2)
    end if

    if (hist_fld_active('Q1000')      .or. &
        hist_fld_active('THE9251000') .or. &
        hist_fld_active('THE8501000') .or. &
        hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 100000._r8, state%q(1,1,1), p_surf_q1)
    end if

    if (hist_fld_active('Q925')       .or. &
        hist_fld_active('THE9251000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 92500._r8, state%q(1,1,1), p_surf_q2)
    end if

    !!! at 1000 mb and 925 mb
    if (hist_fld_active('T1000')) then
       call outfld('T1000    ', p_surf_t1, pcols, lchnk )
    end if

    if (hist_fld_active('T925')) then
       call outfld('T925    ', p_surf_t2, pcols, lchnk )
    end if

    if (hist_fld_active('Q1000')) then
       call outfld('Q1000    ', p_surf_q1, pcols, lchnk )
    end if

    if (hist_fld_active('Q925')) then
       call outfld('Q925    ', p_surf_q2, pcols, lchnk )
    end if

    if (hist_fld_active('T9251000')) then
       p_surf = p_surf_t2-p_surf_t1
       call outfld('T9251000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('TH9251000')) then
       p_surf = (p_surf_t2*(1000.0_r8/925.0_r8)**cappa)-(p_surf_t1*(1.0_r8)**cappa)
       call outfld('TH9251000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('THE9251000')) then
       p_surf = (p_surf_t2*(1000.0_r8/925.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q2)/(1004.0_r8*p_surf_t2))- &
            (p_surf_t1*(1.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q1)/(1004.0_r8*p_surf_t1))
       call outfld('THE9251000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('T8501000')  .or. &
        hist_fld_active('TH8501000') .or. &
        hist_fld_active('THE8501000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%t, p_surf_t2)
    end if

    !!! at 1000 mb and 850 mb
    if (hist_fld_active('T8501000')) then
       p_surf = p_surf_t2-p_surf_t1
       call outfld('T8501000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('TH8501000')) then
       p_surf = (p_surf_t2*(1000.0_r8/850.0_r8)**cappa)-(p_surf_t1*(1.0_r8)**cappa)
       call outfld('TH8501000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('THE8501000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%q(1,1,1), p_surf_q2)
       p_surf = (p_surf_t2*(1000.0_r8/850.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q2)/(1004.0_r8*p_surf_t2))- &
            (p_surf_t1*(1.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q1)/(1004.0_r8*p_surf_t1))
       call outfld('THE8501000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('T7001000')  .or. &
        hist_fld_active('TH7001000') .or. &
        hist_fld_active('T700') .or. &
        hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, state%t, p_surf_t2)
    end if

   !!! at 700 mb
    if (hist_fld_active('T700')) then
       call outfld('T700    ', p_surf_t2, pcols, lchnk )
    end if

    !!! at 1000 mb and 700 mb
    if (hist_fld_active('T7001000')) then
       p_surf = p_surf_t2-p_surf_t1
       call outfld('T7001000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('TH7001000')) then
       p_surf = (p_surf_t2*(1000.0_r8/700.0_r8)**cappa)-(p_surf_t1*(1.0_r8)**cappa)
       call outfld('TH7001000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, state%q(1,1,1), p_surf_q2)
       p_surf = (p_surf_t2*(1000.0_r8/700.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q2)/(1004.0_r8*p_surf_t2))- &
            (p_surf_t1*(1.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q1)/(1004.0_r8*p_surf_t1))
       call outfld('THE7001000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('T010')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 1000._r8, state%t, p_surf)
       call outfld('T010           ', p_surf, pcols, lchnk )
    end if


    !---------------------------------------------------------
    ! tidal diagnostics
    !---------------------------------------------------------
    call tidal_diag_write(state)

    return
  end subroutine diag_phys_writeout
!===============================================================================

  subroutine diag_phys_writeout_fill_z3(ncol, zm, phis, z3)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use physconst,     only: rga

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: zm(pcols,pver)
    real(r8), target, intent(in) :: phis(pcols)
    real(r8), target, intent(out) :: z3(pcols,pver)

    integer :: k

    interface
       subroutine diag_phys_writeout_basic_2d_codon(mode_c, ncol_c, pcols_c, pver_c, scale_c, a_p, b_p, out_p) &
            bind(c, name="diag_phys_writeout_basic_2d_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scale_c
         type(c_ptr), value :: a_p, b_p, out_p
       end subroutine diag_phys_writeout_basic_2d_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       do k = 1, pver
          z3(:ncol,k) = zm(:ncol,k) + phis(:ncol)*rga
       end do
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         1_c_int64_t, 1_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), real(rga, c_double), 0._c_double, 0._c_double, &
         c_loc(zm), c_loc(phis), c_loc(zm), c_loc(zm), c_loc(zm), c_loc(zm), c_loc(z3), c_loc(z3), c_loc(z3) &
    )

  end subroutine diag_phys_writeout_fill_z3

!===============================================================================

  subroutine diag_phys_writeout_mul_field(ncol, a, b, out)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: a(pcols,pver), b(pcols,pver)
    real(r8), target, intent(out) :: out(pcols,pver)

    integer :: k

    interface
       subroutine diag_phys_writeout_basic_2d_codon(mode_c, ncol_c, pcols_c, pver_c, scale_c, a_p, b_p, out_p) &
            bind(c, name="diag_phys_writeout_basic_2d_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scale_c
         type(c_ptr), value :: a_p, b_p, out_p
       end subroutine diag_phys_writeout_basic_2d_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       do k = 1, pver
          out(:ncol,k) = a(:ncol,k) * b(:ncol,k)
       end do
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         1_c_int64_t, 2_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), 0._c_double, 0._c_double, 0._c_double, &
         c_loc(a), c_loc(b), c_loc(a), c_loc(a), c_loc(a), c_loc(a), c_loc(out), c_loc(out), c_loc(out) &
    )

  end subroutine diag_phys_writeout_mul_field

!===============================================================================

  subroutine diag_phys_writeout_mul_scalar_field(ncol, a, b, scale, out)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: a(pcols,pver), b(pcols,pver)
    real(r8), intent(in) :: scale
    real(r8), target, intent(out) :: out(pcols,pver)

    integer :: k

    interface
       subroutine diag_phys_writeout_basic_2d_codon(mode_c, ncol_c, pcols_c, pver_c, scale_c, a_p, b_p, out_p) &
            bind(c, name="diag_phys_writeout_basic_2d_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scale_c
         type(c_ptr), value :: a_p, b_p, out_p
       end subroutine diag_phys_writeout_basic_2d_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       do k = 1, pver
          out(:ncol,k) = a(:ncol,k) * b(:ncol,k) * scale
       end do
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         1_c_int64_t, 3_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), real(scale, c_double), 0._c_double, 0._c_double, &
         c_loc(a), c_loc(b), c_loc(a), c_loc(a), c_loc(a), c_loc(a), c_loc(out), c_loc(out), c_loc(out) &
    )

  end subroutine diag_phys_writeout_mul_scalar_field

!===============================================================================

  subroutine diag_phys_writeout_square_field(ncol, a, out)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: a(pcols,pver)
    real(r8), target, intent(out) :: out(pcols,pver)

    integer :: k

    interface
       subroutine diag_phys_writeout_basic_2d_codon(mode_c, ncol_c, pcols_c, pver_c, scale_c, a_p, b_p, out_p) &
            bind(c, name="diag_phys_writeout_basic_2d_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scale_c
         type(c_ptr), value :: a_p, b_p, out_p
       end subroutine diag_phys_writeout_basic_2d_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       do k = 1, pver
          out(:ncol,k) = a(:ncol,k) * a(:ncol,k)
       end do
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         1_c_int64_t, 4_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), 0._c_double, 0._c_double, 0._c_double, &
         c_loc(a), c_loc(a), c_loc(a), c_loc(a), c_loc(a), c_loc(a), c_loc(out), c_loc(out), c_loc(out) &
    )

  end subroutine diag_phys_writeout_square_field

!===============================================================================

  subroutine diag_phys_writeout_wspeed_field(ncol, u, v, out)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: u(pcols,pver), v(pcols,pver)
    real(r8), target, intent(out) :: out(pcols,pver)

    interface
       subroutine diag_phys_writeout_basic_2d_codon(mode_c, ncol_c, pcols_c, pver_c, scale_c, a_p, b_p, out_p) &
            bind(c, name="diag_phys_writeout_basic_2d_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scale_c
         type(c_ptr), value :: a_p, b_p, out_p
       end subroutine diag_phys_writeout_basic_2d_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       out(:ncol,:) = sqrt(u(:ncol,:)*u(:ncol,:) + v(:ncol,:)*v(:ncol,:))
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         1_c_int64_t, 5_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), 0._c_double, 0._c_double, 0._c_double, &
         c_loc(u), c_loc(v), c_loc(u), c_loc(u), c_loc(u), c_loc(u), c_loc(out), c_loc(out), c_loc(out) &
    )

  end subroutine diag_phys_writeout_wspeed_field

!===============================================================================

  subroutine diag_phys_writeout_mass_and_tmq(ncol, q, pdel, mq, tmq)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use physconst,     only: rga

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: q(pcols,pver), pdel(pcols,pver)
    real(r8), target, intent(out) :: mq(pcols,pver)
    real(r8), target, intent(out) :: tmq(pcols)

    integer :: k

    interface
       subroutine diag_phys_writeout_column_reduce_codon(mode_c, ncol_c, pcols_c, pver_c, scalar1_c, scalar2_c, &
            scalar3_c, a_p, b_p, c_p, d_p, e_p, f_p, out2d_p, out1d_p) bind(c, name="diag_phys_writeout_column_reduce_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scalar1_c, scalar2_c, scalar3_c
         type(c_ptr), value :: a_p, b_p, c_p, d_p, e_p, f_p, out2d_p, out1d_p
       end subroutine diag_phys_writeout_column_reduce_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       mq(:ncol,:) = q(:ncol,:) * pdel(:ncol,:) * rga
       tmq(:ncol) = mq(:ncol,1)
       do k = 2, pver
          tmq(:ncol) = tmq(:ncol) + mq(:ncol,k)
       end do
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         2_c_int64_t, 1_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), real(rga, c_double), 0._c_double, 0._c_double, &
         c_loc(q), c_loc(pdel), c_loc(q), c_loc(q), c_loc(q), c_loc(q), c_loc(mq), c_loc(tmq), c_loc(tmq) &
    )

  end subroutine diag_phys_writeout_mass_and_tmq

!===============================================================================

  subroutine diag_phys_writeout_atmeint(ncol, t, q, u, v, pdel, phis, atmeint)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use physconst,     only: cpair, latvap, gravit

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: t(pcols,pver), q(pcols,pver), u(pcols,pver), v(pcols,pver), pdel(pcols,pver)
    real(r8), target, intent(in) :: phis(pcols)
    real(r8), target, intent(out) :: atmeint(pcols)

    integer :: k
    real(r8) :: work(pcols,pver)

    interface
       subroutine diag_phys_writeout_column_reduce_codon(mode_c, ncol_c, pcols_c, pver_c, scalar1_c, scalar2_c, &
            scalar3_c, a_p, b_p, c_p, d_p, e_p, f_p, out2d_p, out1d_p) bind(c, name="diag_phys_writeout_column_reduce_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scalar1_c, scalar2_c, scalar3_c
         type(c_ptr), value :: a_p, b_p, c_p, d_p, e_p, f_p, out2d_p, out1d_p
       end subroutine diag_phys_writeout_column_reduce_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       do k = 1, pver
          work(:ncol,k) = (cpair*t(:ncol,k) + phis(:ncol) + latvap*q(:ncol,k) + &
               0.5_r8*(u(:ncol,k)*u(:ncol,k) + v(:ncol,k)*v(:ncol,k))) * (pdel(:ncol,k)/gravit)
       end do
       atmeint(:ncol) = work(:ncol,1)
       do k = 2, pver
          atmeint(:ncol) = atmeint(:ncol) + work(:ncol,k)
       end do
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         2_c_int64_t, 2_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), real(cpair, c_double), real(latvap, c_double), real(gravit, c_double), &
         c_loc(t), c_loc(q), c_loc(u), c_loc(v), c_loc(pdel), c_loc(phis), c_loc(t), c_loc(atmeint), c_loc(atmeint) &
    )

  end subroutine diag_phys_writeout_atmeint

!===============================================================================

  subroutine diag_phys_writeout_wtrc_column(ncol, qtr, wind, pdel, rga_in, mode, out)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol, mode
    real(r8), intent(in) :: rga_in
    real(r8), target, intent(in) :: qtr(pcols,pver), wind(pcols,pver), pdel(pcols,pver)
    real(r8), target, intent(out) :: out(pcols,pver)

    integer :: k

    interface
       subroutine diag_phys_writeout_transport_moisture_codon(mode_c, submode_c, ncol_c, pcols_c, pver_c, &
            scalar_c, a_p, b_p, c_p, d_p, out1_p, out2_p, out3_p) bind(c, name="diag_phys_writeout_transport_moisture_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, submode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scalar_c
         type(c_ptr), value :: a_p, b_p, c_p, d_p, out1_p, out2_p, out3_p
       end subroutine diag_phys_writeout_transport_moisture_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       if (mode == 1) then
          out(:ncol,:) = qtr(:ncol,:) * pdel(:ncol,:) * rga_in
       else
          out(:ncol,:) = wind(:ncol,:)*qtr(:ncol,:) * pdel(:ncol,:) * rga_in
       end if
       do k=2,pver
          out(:ncol,1) = out(:ncol,1) + out(:ncol,k)
       end do
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         3_c_int64_t, 1_c_int64_t, int(mode, c_int64_t), int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), real(rga_in, c_double), 0._c_double, 0._c_double, &
         c_loc(qtr), c_loc(wind), c_loc(pdel), c_loc(qtr), c_loc(qtr), c_loc(qtr), c_loc(out), &
         c_loc(out), c_loc(out) &
    )

  end subroutine diag_phys_writeout_wtrc_column

!===============================================================================

  subroutine diag_phys_writeout_ivt(ncol, q, u, v, pdel, rga_in, uqdp, vqdp, ivt)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), intent(in) :: rga_in
    real(r8), target, intent(in) :: q(pcols,pver), u(pcols,pver), v(pcols,pver), pdel(pcols,pver)
    real(r8), target, intent(out) :: uqdp(pcols,pver), vqdp(pcols,pver), ivt(pcols,pver)

    integer :: k

    interface
       subroutine diag_phys_writeout_transport_moisture_codon(mode_c, submode_c, ncol_c, pcols_c, pver_c, &
            scalar_c, a_p, b_p, c_p, d_p, out1_p, out2_p, out3_p) bind(c, name="diag_phys_writeout_transport_moisture_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, submode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scalar_c
         type(c_ptr), value :: a_p, b_p, c_p, d_p, out1_p, out2_p, out3_p
       end subroutine diag_phys_writeout_transport_moisture_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       uqdp(:ncol,:) = q(:ncol,:) * u(:ncol,:) *pdel(:ncol,:) * rga_in
       vqdp(:ncol,:) = q(:ncol,:) * v(:ncol,:) *pdel(:ncol,:) * rga_in
       do k=2,pver
          uqdp(:ncol,1) = uqdp(:ncol,1) + uqdp(:ncol,k)
          vqdp(:ncol,1) = vqdp(:ncol,1) + vqdp(:ncol,k)
       end do
       ivt(:ncol,1) = sqrt( uqdp(:ncol,1)**2 + vqdp(:ncol,1)**2)
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         3_c_int64_t, 2_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), real(rga_in, c_double), 0._c_double, 0._c_double, &
         c_loc(q), c_loc(u), c_loc(v), c_loc(pdel), c_loc(q), c_loc(q), c_loc(uqdp), c_loc(vqdp), c_loc(ivt) &
    )

  end subroutine diag_phys_writeout_ivt

!===============================================================================

  subroutine diag_phys_writeout_copy_col1(ncol, src, dst)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: src(pcols,pver)
    real(r8), target, intent(inout) :: dst(pcols,pver)

    interface
       subroutine diag_phys_writeout_basic_2d_codon(mode_c, ncol_c, pcols_c, pver_c, scale_c, a_p, b_p, out_p) &
            bind(c, name="diag_phys_writeout_basic_2d_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scale_c
         type(c_ptr), value :: a_p, b_p, out_p
       end subroutine diag_phys_writeout_basic_2d_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       dst(:ncol,1) = src(:ncol,1)
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         1_c_int64_t, 6_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), 0._c_double, 0._c_double, 0._c_double, &
         c_loc(src), c_loc(src), c_loc(src), c_loc(src), c_loc(src), c_loc(src), c_loc(dst), c_loc(dst), c_loc(dst) &
    )

  end subroutine diag_phys_writeout_copy_col1

!===============================================================================

  subroutine diag_phys_writeout_scale_relhum(ncol, q, rh)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: q(pcols,pver)
    real(r8), target, intent(inout) :: rh(pcols,pver)

    interface
       subroutine diag_phys_writeout_basic_2d_codon(mode_c, ncol_c, pcols_c, pver_c, scale_c, a_p, b_p, out_p) &
            bind(c, name="diag_phys_writeout_basic_2d_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scale_c
         type(c_ptr), value :: a_p, b_p, out_p
       end subroutine diag_phys_writeout_basic_2d_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       rh(:ncol,:) = q(:ncol,:)/rh(:ncol,:)*100._r8
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         1_c_int64_t, 7_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), 100._c_double, 0._c_double, 0._c_double, &
         c_loc(q), c_loc(rh), c_loc(q), c_loc(q), c_loc(q), c_loc(q), c_loc(rh), c_loc(rh), c_loc(rh) &
    )

  end subroutine diag_phys_writeout_scale_relhum

!===============================================================================

  subroutine diag_phys_writeout_rhi_rhcfmip(ncol, t, esl, esi, rhw, rhi, rhcfmip)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    integer, intent(in) :: ncol
    real(r8), target, intent(in) :: t(pcols,pver), esl(pcols,pver), esi(pcols,pver), rhw(pcols,pver)
    real(r8), target, intent(out) :: rhi(pcols,pver), rhcfmip(pcols,pver)

    integer :: i, k

    interface
       subroutine diag_phys_writeout_transport_moisture_codon(mode_c, submode_c, ncol_c, pcols_c, pver_c, &
            scalar_c, a_p, b_p, c_p, d_p, out1_p, out2_p, out3_p) bind(c, name="diag_phys_writeout_transport_moisture_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: mode_c, submode_c, ncol_c, pcols_c, pver_c
         real(c_double), value :: scalar_c
         type(c_ptr), value :: a_p, b_p, c_p, d_p, out1_p, out2_p, out3_p
       end subroutine diag_phys_writeout_transport_moisture_codon
    end interface

    if (diag_phys_writeout_use_native_impl) then
       do i=1,ncol
          do k=1,pver
             rhi(i,k)=rhw(i,k)*esl(i,k)/esi(i,k)
          end do
       end do

       rhcfmip(:ncol,:)=rhw(:ncol,:)

       do i=1,ncol
          do k=1,pver
             if (t(i,k) .gt. 273) then
                rhcfmip(i,k)=rhw(i,k)
             else
                rhcfmip(i,k)=rhi(i,k)
             end if
          end do
       end do
       return
    end if

    call diag_phys_writeout_batch_dispatch_codon( &
         3_c_int64_t, 3_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), &
         int(pver, c_int64_t), 0._c_double, 0._c_double, 0._c_double, &
         c_loc(t), c_loc(esl), c_loc(esi), c_loc(rhw), c_loc(t), c_loc(t), c_loc(rhi), &
         c_loc(rhcfmip), c_loc(rhcfmip) &
    )

  end subroutine diag_phys_writeout_rhi_rhcfmip

!===============================================================================

subroutine diag_phys_writeout_batch_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line
   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('DIAG_PHYS_WRITEOUT_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine diag_phys_writeout_batch_append_proof

subroutine diag_phys_writeout_batch_log_entered()

   if (diag_phys_writeout_batch_entered_logged) return
   diag_phys_writeout_batch_entered_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'diag_phys_writeout_batch entered (unified stage-dispatch basic/column/transport/relhum/tt direct = codon)'
      call diag_phys_writeout_batch_append_proof( &
           'diag_phys_writeout_batch entered (unified stage-dispatch basic/column/transport/relhum/tt ' // &
           'direct = codon)')
      call flush(iulog)
   end if

end subroutine diag_phys_writeout_batch_log_entered

!===============================================================================

subroutine diag_phys_writeout_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (diag_phys_writeout_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('DIAG_PHYS_WRITEOUT_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      diag_phys_writeout_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      diag_phys_writeout_use_native_impl = .false.
   end if

   diag_phys_writeout_impl_selected = .true.

   if (masterproc) then
      if (diag_phys_writeout_use_native_impl) then
         write(iulog,*) 'diag_phys_writeout implementation = native'
      else
         write(iulog,*) 'diag_phys_writeout implementation = codon'
      end if
   end if

end subroutine diag_phys_writeout_select_impl

!===============================================================================

subroutine diag_conv(state, ztodt, pbuf)

!-----------------------------------------------------------------------
!
! Output diagnostics associated with all convective processes.
!
!-----------------------------------------------------------------------
   use physconst,     only: cpair
   use tidal_diag,    only: get_tidal_coeffs
   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

! Arguments:

   real(r8),            intent(in) :: ztodt   ! timestep for computing physics tendencies
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)

! convective precipitation variables
   real(r8), pointer :: prec_dp(:)                 ! total precipitation   from ZM convection
   real(r8), pointer :: snow_dp(:)                 ! snow from ZM   convection
   real(r8), pointer :: prec_sh(:)                 ! total precipitation   from Hack convection
   real(r8), pointer :: snow_sh(:)                 ! snow from   Hack   convection
   real(r8), pointer :: prec_sed(:)                ! total precipitation   from ZM convection
   real(r8), pointer :: snow_sed(:)                ! snow from ZM   convection
   real(r8), pointer :: prec_pcw(:)                ! total precipitation   from Hack convection
   real(r8), pointer :: snow_pcw(:)                ! snow from Hack   convection

   !water tracers/isotopes:
   real(r8), pointer :: wtprec(:)                  !water tracer precipitation
   real(r8), pointer :: wtprec_cvrain(:), wtprec_cvsnow(:), wtprec_strain(:), wtprec_stsnow(:)
   real(r8), target  :: wtprect(pcols)             !total water tracer precipitation

! Local variables:

   integer :: i, k, m, lchnk, ncol

   real(r8) :: rtdt

   real(r8), target :: precc(pcols)       ! convective precip rate
   real(r8), target :: precl(pcols)       ! stratiform precip rate
   real(r8), target :: snowc(pcols)       ! convective snow rate
   real(r8), target :: snowl(pcols)       ! stratiform snow rate
   real(r8), target :: prect(pcols)       ! total (conv+large scale) precip rate
   real(r8), target :: dqcond_work(pcols,pver)
   real(r8), target :: tidal_work(pcols,pver)
   real(r8) :: dcoef(4)                   ! for tidal component of T tend

   lchnk = state%lchnk
   ncol  = state%ncol

   rtdt = 1._r8/ztodt

   call cam_diag_conv_batch_select_impl()
   if (.not. cam_diag_conv_batch_use_native_impl) then
      if (diag_conv_codon(7_c_int64_t) == 7_c_int64_t) then
         call cam_diag_log_direct(diag_conv_logged, &
              'diag_conv direct = codon; convective diagnostic arithmetic direct = codon; outfld/pbuf native CAM API islands')
      end if
   end if

   call pbuf_get_field(pbuf, prec_dp_idx, prec_dp)
   call pbuf_get_field(pbuf, snow_dp_idx, snow_dp)
   call pbuf_get_field(pbuf, prec_sh_idx, prec_sh)
   call pbuf_get_field(pbuf, snow_sh_idx, snow_sh)
   call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
   call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
   call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
   call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)

! Precipitation rates (multi-process)
   if (cam_diag_conv_batch_use_native_impl) then
      precc(:ncol) = prec_dp(:ncol)  + prec_sh(:ncol)
      precl(:ncol) = prec_sed(:ncol) + prec_pcw(:ncol)
      snowc(:ncol) = snow_dp(:ncol)  + snow_sh(:ncol)
      snowl(:ncol) = snow_sed(:ncol) + snow_pcw(:ncol)
      prect(:ncol) = precc(:ncol)    + precl(:ncol)
   else
      call cam_diag_conv_batch_log_diag_conv_entered()
      call cam_diag_conv_batch_log_precip_dtcond_entered()
      call diag_conv_batch_dispatch_codon( &
           3_c_int64_t, 0_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           int(pcnst, c_int64_t), 0_c_int64_t, &
           real(rtdt, c_double), real(cpair, c_double), &
           c_loc(prec_dp), c_loc(snow_dp), c_loc(prec_sh), c_loc(snow_sh), &
           c_loc(prec_sed), c_loc(snow_sed), c_loc(prec_pcw), c_loc(snow_pcw), &
           c_loc(precc), c_loc(precl), c_loc(snowc), c_loc(snowl), c_loc(prect), &
           c_loc(state%s), c_loc(dtcond(1,1,lchnk)) &
      )
   end if

   call outfld('PRECC   ', precc, pcols, lchnk )
   call outfld('PRECL   ', precl, pcols, lchnk )
   call outfld('PREC_PCW', prec_pcw,pcols   ,lchnk )
   call outfld('PREC_zmc', prec_dp ,pcols   ,lchnk )
   call outfld('PRECSC  ', snowc, pcols, lchnk )
   call outfld('PRECSL  ', snowl, pcols, lchnk )
   call outfld('PRECT   ', prect, pcols, lchnk )
   call outfld('PRECTMX ', prect, pcols, lchnk )

   !**********************
   !Water tracers/isotopes
   !**********************
    if(trace_water) then
      do m=1,wtrc_nwset
        !convective rain:
        call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvrain,m), wtprec_cvrain)
        !convective snow:
        call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvsnow,m), wtprec_cvsnow)
        !stratiform rain:
        call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstrain,m), wtprec_strain)
        !stratiform snow:
        call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtstsnow,m), wtprec_stsnow)
        if (cam_diag_conv_batch_use_native_impl) then
           wtprect(:ncol) = wtprec_cvrain(:ncol) !add to sum
           wtprect(:ncol) = wtprect(:ncol) + wtprec_cvsnow(:ncol)
           wtprect(:ncol) = wtprect(:ncol) + wtprec_strain(:ncol)
           wtprect(:ncol) = wtprect(:ncol) + wtprec_stsnow(:ncol)
        else
           call diag_conv_batch_dispatch_codon( &
                2_c_int64_t, 2_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
                int(pcnst, c_int64_t), 0_c_int64_t, 0._c_double, 0._c_double, &
                c_loc(wtprec_cvrain), c_loc(wtprec_cvsnow), c_loc(wtprec_strain), c_loc(wtprec_stsnow), &
                c_loc(wtprec_cvrain), c_loc(wtprec_cvsnow), c_loc(wtprec_strain), c_loc(wtprec_stsnow), &
                c_loc(wtprect), c_loc(wtprect), c_loc(wtprect), c_loc(wtprect), c_loc(wtprect), &
                c_loc(wtprec_cvrain), c_loc(wtprect) &
           )
        end if
        !add to output variable:
        call outfld('PRECT_'//trim(wtrc_out_names(m)), wtprect, pcols, lchnk)
      end do
    end  if
   !**********************

   call outfld('PRECLav ', precl, pcols, lchnk )
   call outfld('PRECCav ', precc, pcols, lchnk )

#if ( defined BFB_CAM_SCAM_IOP )
   call outfld('Prec   ' , prect, pcols, lchnk )
#endif

   ! Total convection tendencies.

   if (cam_diag_conv_batch_use_native_impl) then
      do k = 1, pver
         do i = 1, ncol
            dtcond(i,k,lchnk) = (state%s(i,k) - dtcond(i,k,lchnk))*rtdt / cpair
         end do
      end do
   end if
   call outfld('DTCOND  ', dtcond(:,:,lchnk), pcols, lchnk)

   ! output tidal coefficients
   call get_tidal_coeffs( dcoef )
   if (cam_diag_conv_batch_use_native_impl) then
      call outfld( 'DTCOND_24_SIN', dtcond(:ncol,:,lchnk)*dcoef(1), ncol, lchnk )
      call outfld( 'DTCOND_24_COS', dtcond(:ncol,:,lchnk)*dcoef(2), ncol, lchnk )
      call outfld( 'DTCOND_12_SIN', dtcond(:ncol,:,lchnk)*dcoef(3), ncol, lchnk )
      call outfld( 'DTCOND_12_COS', dtcond(:ncol,:,lchnk)*dcoef(4), ncol, lchnk )
   else
      call diag_conv_scale_2d_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           real(dcoef(1), c_double), c_loc(dtcond(1,1,lchnk)), c_loc(tidal_work(1,1)))
      call outfld( 'DTCOND_24_SIN', tidal_work, ncol, lchnk )
      call diag_conv_scale_2d_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           real(dcoef(2), c_double), c_loc(dtcond(1,1,lchnk)), c_loc(tidal_work(1,1)))
      call outfld( 'DTCOND_24_COS', tidal_work, ncol, lchnk )
      call diag_conv_scale_2d_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           real(dcoef(3), c_double), c_loc(dtcond(1,1,lchnk)), c_loc(tidal_work(1,1)))
      call outfld( 'DTCOND_12_SIN', tidal_work, ncol, lchnk )
      call diag_conv_scale_2d_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           real(dcoef(4), c_double), c_loc(dtcond(1,1,lchnk)), c_loc(tidal_work(1,1)))
      call outfld( 'DTCOND_12_COS', tidal_work, ncol, lchnk )
   end if

   do m = 1, dqcond_num
      if ( cnst_cam_outfld(m) ) then
         if (cam_diag_conv_batch_use_native_impl) then
            do k = 1, pver
               do i = 1, ncol
                  dqcond(m)%cnst(i,k,lchnk) = (state%q(i,k,m) - dqcond(m)%cnst(i,k,lchnk))*rtdt
               end do
            end do
         else
            dqcond_work(:,:) = dqcond(m)%cnst(:,:,lchnk)
            call diag_conv_batch_dispatch_codon( &
                 2_c_int64_t, 4_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
                 int(pcnst, c_int64_t), int(m, c_int64_t), real(rtdt, c_double), 0._c_double, &
                 c_loc(state%q), c_loc(state%q), c_loc(state%q), c_loc(state%q), &
                 c_loc(state%q), c_loc(state%q), c_loc(state%q), c_loc(state%q), &
                 c_loc(dqcond_work), c_loc(dqcond_work), c_loc(dqcond_work), c_loc(dqcond_work), c_loc(dqcond_work), &
                 c_loc(state%q), c_loc(dqcond_work) &
            )
            dqcond(m)%cnst(:,:,lchnk) = dqcond_work(:,:)
         end if
         call outfld(dcconnam(m), dqcond(m)%cnst(:,:,lchnk), pcols, lchnk)
      end if
   end do

end subroutine diag_conv

!===============================================================================

subroutine diag_surf_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (diag_surf_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('DIAG_SURF_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      diag_surf_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      diag_surf_use_native_impl = .false.
   end if

   diag_surf_impl_selected = .true.

   if (masterproc) then
      if (diag_surf_use_native_impl) then
         write(iulog,*) 'diag_surf implementation = native'
      else
         write(iulog,*) 'diag_surf implementation = codon'
      end if
   end if

end subroutine diag_surf_select_impl

subroutine diag_physvar_ic_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (diag_physvar_ic_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('DIAG_PHYSVAR_IC_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      diag_physvar_ic_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      diag_physvar_ic_use_native_impl = .false.
   end if

   diag_physvar_ic_impl_selected = .true.

   if (masterproc) then
      if (diag_physvar_ic_use_native_impl) then
         write(iulog,*) 'diag_physvar_ic implementation = native'
      else
         write(iulog,*) 'diag_physvar_ic implementation = codon'
      end if
   end if

end subroutine diag_physvar_ic_select_impl

subroutine diag_surf (cam_in, cam_out, ps, trefmxav, trefmnav )

!-----------------------------------------------------------------------
!
! Purpose: record surface diagnostics
!
!-----------------------------------------------------------------------

   use iso_c_binding,   only: c_double, c_int64_t, c_loc, c_ptr
   use time_manager,     only: is_end_curr_day
   use co2_cycle,        only: c_i, co2_transport
   use constituents,     only: sflxnam

!-----------------------------------------------------------------------
!
! Input arguments
!
    type(cam_in_t),  target, intent(in) :: cam_in
    type(cam_out_t), intent(in) :: cam_out

    real(r8), target, intent(inout) :: trefmnav(pcols) ! daily minimum tref
    real(r8), target, intent(inout) :: trefmxav(pcols) ! daily maximum tref

    real(r8), intent(in)    :: ps(pcols)       ! Surface pressure.
!
!---------------------------Local workspace-----------------------------
!
    integer :: i, k, m      ! indexes
    integer :: lchnk        ! chunk identifier
    integer :: ncol         ! longitude dimension
    real(r8), target :: tem2(pcols)    ! temporary workspace
    real(r8), target :: ftem(pcols)    ! temporary workspace
    real(r8), target :: trefmx_day(pcols)
    real(r8), target :: trefmn_day(pcols)
    interface
       subroutine diag_surf_codon(ncol_c, pcols_c, end_day_c, qref_p, rhref_p, tref_p, &
            trefmxav_p, trefmnav_p, trefmx_day_p, trefmn_day_p) bind(c, name="diag_surf_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, end_day_c
         type(c_ptr), value :: qref_p, rhref_p, tref_p, trefmxav_p, trefmnav_p, trefmx_day_p, trefmn_day_p
       end subroutine diag_surf_codon
    end interface
!
!-----------------------------------------------------------------------
!
    call diag_surf_select_impl()

    if (diag_surf_use_native_impl) then
       call diag_surf_native(cam_in, cam_out, ps, trefmxav, trefmnav)
       return
    end if

    lchnk = cam_in%lchnk
    ncol  = cam_in%ncol

    call outfld('SHFLX',    cam_in%shf,       pcols, lchnk)
    call outfld('LHFLX',    cam_in%lhf,       pcols, lchnk)
    call outfld('QFLX',     cam_in%cflx(1,1), pcols, lchnk)

    call outfld('TAUX',     cam_in%wsx,       pcols, lchnk)
    call outfld('TAUY',     cam_in%wsy,       pcols, lchnk)
    call outfld('TREFHT  ', cam_in%tref,      pcols, lchnk)
    call outfld('TREFHTMX', cam_in%tref,      pcols, lchnk)
    call outfld('TREFHTMN', cam_in%tref,      pcols, lchnk)
    call outfld('QREFHT',   cam_in%qref,      pcols, lchnk)
    call outfld('U10',      cam_in%u10,       pcols, lchnk)

    ! Water tracers:
    if(trace_water) then
      do m=1,wtrc_nwset
        call outfld ('QFLX_'//trim(wtrc_out_names(m)), cam_in%cflx(1,wtrc_iatype(m,iwtvap)), pcols, lchnk)
      end do

      call outfld('buckH',  cam_in%buckH, pcols, lchnk)
      call outfld('buckD',  cam_in%buckD, pcols, lchnk)
      call outfld('buck16', cam_in%buck16, pcols, lchnk)
      call outfld('buck18', cam_in%buck18, pcols, lchnk)
    end if

!
! Calculate and output reference height RH (RHREFHT)

   call qsat(cam_in%tref(:ncol), ps(:ncol), tem2(:ncol), ftem(:ncol))
   trefmx_day(:) = trefmxav(:)
   trefmn_day(:) = trefmnav(:)
   call diag_surf_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), merge(1_c_int64_t, 0_c_int64_t, is_end_curr_day()), &
        c_loc(cam_in%qref), c_loc(ftem), c_loc(cam_in%tref), &
        c_loc(trefmxav), c_loc(trefmnav), c_loc(trefmx_day), c_loc(trefmn_day) &
   )
   call cam_diag_log_direct(diag_surf_logged, &
        'diag_surf direct = codon; RHREFHT/daily TREF minmax numeric body direct = codon; outfld/qsat/CO2 native CAM API islands')


    call outfld('RHREFHT',   ftem,      pcols, lchnk)


#if (defined BFB_CAM_SCAM_IOP )
    call outfld('shflx   ',cam_in%shf,   pcols,   lchnk)
    call outfld('lhflx   ',cam_in%lhf,   pcols,   lchnk)
    call outfld('trefht  ',cam_in%tref,  pcols,   lchnk)
#endif
!
! Ouput ocn and ice fractions
!
    call outfld('LANDFRAC', cam_in%landfrac, pcols, lchnk)
    call outfld('ICEFRAC',  cam_in%icefrac,  pcols, lchnk)
    call outfld('OCNFRAC',  cam_in%ocnfrac,  pcols, lchnk)
!
! Compute daily minimum and maximum of TREF
!
    if (is_end_curr_day()) then
       call outfld('TREFMXAV', trefmx_day,pcols,   lchnk     )
       call outfld('TREFMNAV', trefmn_day,pcols,   lchnk     )
    endif

    call outfld('TBOT',     cam_out%tbot,     pcols, lchnk)
    call outfld('TS',       cam_in%ts,        pcols, lchnk)
    call outfld('TSMN',     cam_in%ts,        pcols, lchnk)
    call outfld('TSMX',     cam_in%ts,        pcols, lchnk)
    call outfld('SNOWHLND', cam_in%snowhland, pcols, lchnk)
    call outfld('SNOWHICE', cam_in%snowhice,  pcols, lchnk)
    call outfld('ASDIR',    cam_in%asdir,     pcols, lchnk)
    call outfld('ASDIF',    cam_in%asdif,     pcols, lchnk)
    call outfld('ALDIR',    cam_in%aldir,     pcols, lchnk)
    call outfld('ALDIF',    cam_in%aldif,     pcols, lchnk)
    call outfld('SST',      cam_in%sst,       pcols, lchnk)

    if (co2_transport()) then
       do m = 1,4
          call outfld(sflxnam(c_i(m)), cam_in%cflx(:,c_i(m)), pcols, lchnk)
       end do
    end if

end subroutine diag_surf

subroutine diag_surf_native (cam_in, cam_out, ps, trefmxav, trefmnav )

!-----------------------------------------------------------------------
!
! Purpose: record surface diagnostics
!
!-----------------------------------------------------------------------

   use time_manager,     only: is_end_curr_day
   use co2_cycle,        only: c_i, co2_transport
   use constituents,     only: sflxnam

!-----------------------------------------------------------------------
!
! Input arguments
!
    type(cam_in_t),  intent(in) :: cam_in
    type(cam_out_t), intent(in) :: cam_out

    real(r8), intent(inout) :: trefmnav(pcols) ! daily minimum tref
    real(r8), intent(inout) :: trefmxav(pcols) ! daily maximum tref

    real(r8), intent(in)    :: ps(pcols)       ! Surface pressure.
!
!---------------------------Local workspace-----------------------------
!
    integer :: i, k, m      ! indexes
    integer :: lchnk        ! chunk identifier
    integer :: ncol         ! longitude dimension
    real(r8) tem2(pcols)    ! temporary workspace
    real(r8) ftem(pcols)    ! temporary workspace
!
!-----------------------------------------------------------------------
!
    lchnk = cam_in%lchnk
    ncol  = cam_in%ncol

    call outfld('SHFLX',    cam_in%shf,       pcols, lchnk)
    call outfld('LHFLX',    cam_in%lhf,       pcols, lchnk)
    call outfld('QFLX',     cam_in%cflx(1,1), pcols, lchnk)

    call outfld('TAUX',     cam_in%wsx,       pcols, lchnk)
    call outfld('TAUY',     cam_in%wsy,       pcols, lchnk)
    call outfld('TREFHT  ', cam_in%tref,      pcols, lchnk)
    call outfld('TREFHTMX', cam_in%tref,      pcols, lchnk)
    call outfld('TREFHTMN', cam_in%tref,      pcols, lchnk)
    call outfld('QREFHT',   cam_in%qref,      pcols, lchnk)
    call outfld('U10',      cam_in%u10,       pcols, lchnk)

    ! Water tracers:
    if(trace_water) then
      do m=1,wtrc_nwset
        call outfld ('QFLX_'//trim(wtrc_out_names(m)), cam_in%cflx(1,wtrc_iatype(m,iwtvap)), pcols, lchnk)
      end do

      call outfld('buckH',  cam_in%buckH, pcols, lchnk)
      call outfld('buckD',  cam_in%buckD, pcols, lchnk)
      call outfld('buck16', cam_in%buck16, pcols, lchnk)
      call outfld('buck18', cam_in%buck18, pcols, lchnk)
    end if

!
! Calculate and output reference height RH (RHREFHT)

   call qsat(cam_in%tref(:ncol), ps(:ncol), tem2(:ncol), ftem(:ncol))
       ftem(:ncol) = cam_in%qref(:ncol)/ftem(:ncol)*100._r8


    call outfld('RHREFHT',   ftem,      pcols, lchnk)


#if (defined BFB_CAM_SCAM_IOP )
    call outfld('shflx   ',cam_in%shf,   pcols,   lchnk)
    call outfld('lhflx   ',cam_in%lhf,   pcols,   lchnk)
    call outfld('trefht  ',cam_in%tref,  pcols,   lchnk)
#endif
!
! Ouput ocn and ice fractions
!
    call outfld('LANDFRAC', cam_in%landfrac, pcols, lchnk)
    call outfld('ICEFRAC',  cam_in%icefrac,  pcols, lchnk)
    call outfld('OCNFRAC',  cam_in%ocnfrac,  pcols, lchnk)
!
! Compute daily minimum and maximum of TREF
!
    do i = 1,ncol
       trefmxav(i) = max(cam_in%tref(i),trefmxav(i))
       trefmnav(i) = min(cam_in%tref(i),trefmnav(i))
    end do
    if (is_end_curr_day()) then
       call outfld('TREFMXAV', trefmxav,pcols,   lchnk     )
       call outfld('TREFMNAV', trefmnav,pcols,   lchnk     )
       trefmxav(:ncol) = -1.0e36_r8
       trefmnav(:ncol) =  1.0e36_r8
    endif

    call outfld('TBOT',     cam_out%tbot,     pcols, lchnk)
    call outfld('TS',       cam_in%ts,        pcols, lchnk)
    call outfld('TSMN',     cam_in%ts,        pcols, lchnk)
    call outfld('TSMX',     cam_in%ts,        pcols, lchnk)
    call outfld('SNOWHLND', cam_in%snowhland, pcols, lchnk)
    call outfld('SNOWHICE', cam_in%snowhice,  pcols, lchnk)
    call outfld('ASDIR',    cam_in%asdir,     pcols, lchnk)
    call outfld('ASDIF',    cam_in%asdif,     pcols, lchnk)
    call outfld('ALDIR',    cam_in%aldir,     pcols, lchnk)
    call outfld('ALDIF',    cam_in%aldif,     pcols, lchnk)
    call outfld('SST',      cam_in%sst,       pcols, lchnk)

    if (co2_transport()) then
       do m = 1,4
          call outfld(sflxnam(c_i(m)), cam_in%cflx(:,c_i(m)), pcols, lchnk)
       end do
    end if

end subroutine diag_surf_native

!===============================================================================

subroutine diag_export(cam_out)

!-----------------------------------------------------------------------
!
! Purpose: Write export state to history file
!
!-----------------------------------------------------------------------

   use iso_c_binding, only: c_int64_t

   ! arguments
   type(cam_out_t), intent(inout) :: cam_out

   ! Local variables:
   integer :: lchnk        ! chunk identifier
   logical :: atm_dep_flux ! true ==> sending deposition fluxes to coupler.
                           ! Otherwise, set them to zero.
   !-----------------------------------------------------------------------

   lchnk = cam_out%lchnk

   call cam_diag_touch_and_log(5_c_int64_t, diag_export_logged, &
        'diag_export direct = codon; export selector/touch direct = codon; phys_getopts/outfld/coupler native CAM API island')

   call phys_getopts(atm_dep_flux_out=atm_dep_flux)

   if (.not. atm_dep_flux) then
      ! set the fluxes to zero before outfld and sending them to the
      ! coupler
      cam_out%bcphiwet = 0.0_r8
      cam_out%bcphidry = 0.0_r8
      cam_out%bcphodry = 0.0_r8
      cam_out%ocphiwet = 0.0_r8
      cam_out%ocphidry = 0.0_r8
      cam_out%ocphodry = 0.0_r8
      cam_out%dstwet1  = 0.0_r8
      cam_out%dstdry1  = 0.0_r8
      cam_out%dstwet2  = 0.0_r8
      cam_out%dstdry2  = 0.0_r8
      cam_out%dstwet3  = 0.0_r8
      cam_out%dstdry3  = 0.0_r8
      cam_out%dstwet4  = 0.0_r8
      cam_out%dstdry4  = 0.0_r8
   end if

   call outfld('a2x_BCPHIWET', cam_out%bcphiwet, pcols, lchnk)
   call outfld('a2x_BCPHIDRY', cam_out%bcphidry, pcols, lchnk)
   call outfld('a2x_BCPHODRY', cam_out%bcphodry, pcols, lchnk)
   call outfld('a2x_OCPHIWET', cam_out%ocphiwet, pcols, lchnk)
   call outfld('a2x_OCPHIDRY', cam_out%ocphidry, pcols, lchnk)
   call outfld('a2x_OCPHODRY', cam_out%ocphodry, pcols, lchnk)
   call outfld('a2x_DSTWET1',  cam_out%dstwet1,  pcols, lchnk)
   call outfld('a2x_DSTDRY1',  cam_out%dstdry1,  pcols, lchnk)
   call outfld('a2x_DSTWET2',  cam_out%dstwet2,  pcols, lchnk)
   call outfld('a2x_DSTDRY2',  cam_out%dstdry2,  pcols, lchnk)
   call outfld('a2x_DSTWET3',  cam_out%dstwet3,  pcols, lchnk)
   call outfld('a2x_DSTDRY3',  cam_out%dstdry3,  pcols, lchnk)
   call outfld('a2x_DSTWET4',  cam_out%dstwet4,  pcols, lchnk)
   call outfld('a2x_DSTDRY4',  cam_out%dstdry4,  pcols, lchnk)

end subroutine diag_export

!#######################################################################

   subroutine diag_physvar_ic (lchnk,  pbuf, cam_out, cam_in)
   use iso_c_binding, only: c_int64_t
!
!---------------------------------------------
!
! Purpose: record physics variables on IC file
!
!---------------------------------------------
!

!
! Arguments
!
   integer       , intent(in) :: lchnk  ! chunk identifier
   type(physics_buffer_desc), pointer :: pbuf(:)

   type(cam_out_t), intent(inout) :: cam_out
   type(cam_in_t),  intent(inout) :: cam_in
   interface
      function diag_physvar_ic_codon() result(out_c) bind(c, name="diag_physvar_ic_codon")
         use iso_c_binding, only: c_int64_t
         integer(c_int64_t) :: out_c
      end function diag_physvar_ic_codon
   end interface
   integer(c_int64_t) :: diag_physvar_ic_touch_c
!
!-----------------------------------------------------------------------
!
   call diag_physvar_ic_select_impl()

   if (diag_physvar_ic_use_native_impl) then
      call diag_physvar_ic_native(lchnk, pbuf, cam_out, cam_in)
      return
   end if

   diag_physvar_ic_touch_c = diag_physvar_ic_codon()
   if (diag_physvar_ic_touch_c == 1_c_int64_t) then
      call cam_diag_log_direct(diag_physvar_ic_logged, &
           'diag_physvar_ic direct = codon; inithist branch selection direct; pbuf_get_field/outfld native CAM API island')
   end if

   if (write_inithist()) then
      call diag_physvar_ic_native(lchnk, pbuf, cam_out, cam_in)
   end if

   end subroutine diag_physvar_ic

!#######################################################################

   subroutine diag_physvar_ic_native (lchnk,  pbuf, cam_out, cam_in)
!
!---------------------------------------------
!
! Purpose: record physics variables on IC file
!
!---------------------------------------------
!
!
! Arguments
!
   integer       , intent(in) :: lchnk  ! chunk identifier
   type(physics_buffer_desc), pointer :: pbuf(:)

   type(cam_out_t), intent(inout) :: cam_out
   type(cam_in_t),  intent(inout) :: cam_in
!
!---------------------------Local workspace-----------------------------
!
   integer  :: k                 ! indices
   integer  :: itim_old          ! indices

   real(r8), pointer, dimension(:,:) :: cwat_var
   real(r8), pointer, dimension(:,:) :: conv_var_3d
   real(r8), pointer, dimension(:  ) :: conv_var_2d
   real(r8), pointer :: tpert(:), pblh(:), qpert(:)
!
!-----------------------------------------------------------------------
!
   if( write_inithist() ) then

      !
      ! Associate pointers with physics buffer fields
      !
      itim_old = pbuf_old_tim_idx()

      if (qcwat_idx > 0) then
         call pbuf_get_field(pbuf, qcwat_idx, cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
         call outfld('QCWAT&IC   ',cwat_var, pcols,lchnk)
      end if

      if (tcwat_idx > 0) then
         call pbuf_get_field(pbuf, tcwat_idx,  cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
         call outfld('TCWAT&IC   ',cwat_var, pcols,lchnk)
      end if

      if (lcwat_idx > 0) then
         call pbuf_get_field(pbuf, lcwat_idx,  cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
         call outfld('LCWAT&IC   ',cwat_var, pcols,lchnk)
      end if

      call pbuf_get_field(pbuf, cld_idx,    cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      call outfld('CLOUD&IC   ',cwat_var, pcols,lchnk)

      call pbuf_get_field(pbuf, concld_idx, cwat_var, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      call outfld('CONCLD&IC   ',cwat_var, pcols,lchnk)

      call pbuf_get_field(pbuf, tke_idx, conv_var_3d)
      call outfld('TKE&IC    ',conv_var_3d, pcols,lchnk)

      call pbuf_get_field(pbuf, kvm_idx,  conv_var_3d)
      call outfld('KVM&IC    ',conv_var_3d, pcols,lchnk)

      call pbuf_get_field(pbuf, kvh_idx,  conv_var_3d)
      call outfld('KVH&IC    ',conv_var_3d, pcols,lchnk)

      call pbuf_get_field(pbuf, cush_idx, conv_var_2d ,(/1,itim_old/),  (/pcols,1/))
      call outfld('CUSH&IC   ',conv_var_2d, pcols,lchnk)

      if (qpert_idx > 0) then
         call pbuf_get_field(pbuf, qpert_idx, qpert)
         call outfld('QPERT&IC   ', qpert, pcols, lchnk)
      end if

      call pbuf_get_field(pbuf, pblh_idx,  pblh)
      call outfld('PBLH&IC    ', pblh,  pcols, lchnk)

      call pbuf_get_field(pbuf, tpert_idx, tpert)
      call outfld('TPERT&IC   ', tpert, pcols, lchnk)

   end if

   end subroutine diag_physvar_ic_native


!#######################################################################

subroutine diag_phys_tend_writeout(state, pbuf,  tend, ztodt, tmp_q, tmp_cldliq, tmp_cldice, &
                                   tmp_t, qini, cldliqini, cldiceini)

   !---------------------------------------------------------------
   !
   ! Purpose:  Dump physics tendencies for moisture and temperature
   !
   !---------------------------------------------------------------

   use check_energy,    only: check_energy_get_integrals
   use physconst,       only: cpair
   use iso_c_binding,   only: c_double, c_int64_t, c_loc, c_ptr, c_null_ptr

   ! Arguments

   type(physics_state), target, intent(in   ) :: state

   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_tend ), target, intent(in   ) :: tend
   real(r8)           , intent(in   ) :: ztodt                  ! physics timestep
   real(r8) , target   , intent(inout) :: tmp_q     (pcols,pver) ! As input, holds pre-adjusted tracers (FV)
   real(r8) , target   , intent(inout) :: tmp_cldliq(pcols,pver) ! As input, holds pre-adjusted tracers (FV)
   real(r8) , target   , intent(inout) :: tmp_cldice(pcols,pver) ! As input, holds pre-adjusted tracers (FV)
   real(r8) , target   , intent(inout) :: tmp_t     (pcols,pver) ! holds last physics_updated T (FV)
   real(r8) , target   , intent(in   ) :: qini      (pcols,pver) ! tracer fields at beginning of physics
   real(r8) , target   , intent(in   ) :: cldliqini (pcols,pver) ! tracer fields at beginning of physics
   real(r8) , target   , intent(in   ) :: cldiceini (pcols,pver) ! tracer fields at beginning of physics

   !---------------------------Local workspace-----------------------------

   integer  :: m      ! constituent index
   integer  :: lchnk  ! chunk index
   integer  :: ncol   ! number of columns in chunk
   real(r8), target :: ftem2(pcols     ) ! Temporary workspace for outfld variables
   real(r8), target :: ftem3(pcols,pver) ! Temporary workspace for outfld variables
   real(r8) :: rtdt
   real(r8) :: heat_glob         ! global energy integral (FV only)
   integer  :: ixcldice, ixcldliq! constituent indices for cloud liquid and ice water.
   ! CAM pointers to get variables from the physics buffer
   real(r8), pointer, dimension(:,:) :: t_ttend
   integer  :: itim_old
   interface
      subroutine diag_phys_tend_writeout_codon(stage_c, ncol_c, pcols_c, pver_c, pcnst_c, &
           ixcldliq_c, ixcldice_c, ztodt_c, rtdt_c, heat_glob_c, cpair_c, state_t_p, state_q_p, &
           tend_dtdt_p, tmp_t_p, tmp_q_p, tmp_cldliq_p, tmp_cldice_p, qini_p, cldliqini_p, &
           cldiceini_p, ftem2_p, ftem3_p, t_ttend_p) bind(c, name="diag_phys_tend_writeout_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: stage_c, ncol_c, pcols_c, pver_c, pcnst_c
         integer(c_int64_t), value :: ixcldliq_c, ixcldice_c
         real(c_double), value :: ztodt_c, rtdt_c, heat_glob_c, cpair_c
         type(c_ptr), value :: state_t_p, state_q_p, tend_dtdt_p, tmp_t_p, tmp_q_p
         type(c_ptr), value :: tmp_cldliq_p, tmp_cldice_p, qini_p, cldliqini_p, cldiceini_p
         type(c_ptr), value :: ftem2_p, ftem3_p, t_ttend_p
      end subroutine diag_phys_tend_writeout_codon
   end interface

   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   rtdt  = 1._r8/ztodt
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call diag_phys_tend_select_impl()
   if (.not. diag_phys_tend_use_native_impl) then
      call diag_phys_tend_log_entered()
   end if

   ! Dump out post-physics state (FV only)

   if (dycore_is('LR')) then
      if (diag_phys_tend_use_native_impl) then
         tmp_t(:ncol,:pver) = (tmp_t(:ncol,:pver) - state%t(:ncol,:pver))/ztodt
      else
         call diag_phys_tend_writeout_codon( &
              1_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
              int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
              real(ztodt, c_double), real(rtdt, c_double), 0._c_double, real(cpair, c_double), &
              c_loc(state%t), c_null_ptr, c_null_ptr, c_loc(tmp_t), c_loc(tmp_q), c_loc(tmp_cldliq), &
              c_loc(tmp_cldice), c_loc(qini), c_loc(cldliqini), c_loc(cldiceini), c_loc(ftem2), &
              c_loc(ftem3), c_null_ptr &
         )
      end if
      call outfld('PTTEND_RESID', tmp_t, pcols, lchnk   )
   end if
   call outfld('TAP', state%t, pcols, lchnk   )
   call outfld('UAP', state%u, pcols, lchnk   )
   call outfld('VAP', state%v, pcols, lchnk   )

   if ( cnst_cam_outfld(       1) ) call outfld (apcnst(       1), state%q(1,1,       1), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldliq) ) call outfld (apcnst(ixcldliq), state%q(1,1,ixcldliq), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldice) ) call outfld (apcnst(ixcldice), state%q(1,1,ixcldice), pcols, lchnk)

   ! T-tendency due to FV Energy fixer (remove from total physics tendency diagnostic)

   if (dycore_is('LR')) then
      call check_energy_get_integrals( heat_glob_out=heat_glob )
      if (diag_phys_tend_use_native_impl) then
         ftem2(:ncol)  = heat_glob/cpair
      else
         call diag_phys_tend_writeout_codon( &
              2_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
              int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
              real(ztodt, c_double), real(rtdt, c_double), real(heat_glob, c_double), &
              real(cpair, c_double), c_null_ptr, c_null_ptr, c_loc(tend%dtdt), c_loc(tmp_t), &
              c_loc(tmp_q), c_loc(tmp_cldliq), c_loc(tmp_cldice), c_loc(qini), c_loc(cldliqini), &
              c_loc(cldiceini), c_loc(ftem2), c_loc(ftem3), c_null_ptr &
         )
      end if
      call outfld('TFIX', ftem2, pcols, lchnk   )
      if (diag_phys_tend_use_native_impl) then
         ftem3(:ncol,:pver)  = tend%dtdt(:ncol,:pver) - heat_glob/cpair
      end if
   else
      if (diag_phys_tend_use_native_impl) then
         ftem3(:ncol,:pver)  = tend%dtdt(:ncol,:pver)
      else
         call diag_phys_tend_writeout_codon( &
              3_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
              int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
              real(ztodt, c_double), real(rtdt, c_double), 0._c_double, real(cpair, c_double), &
              c_null_ptr, c_null_ptr, c_loc(tend%dtdt), c_loc(tmp_t), c_loc(tmp_q), c_loc(tmp_cldliq), &
              c_loc(tmp_cldice), c_loc(qini), c_loc(cldliqini), c_loc(cldiceini), c_loc(ftem2), &
              c_loc(ftem3), c_null_ptr &
         )
      end if
   end if

   ! Total physics tendency for Temperature

   call outfld('PTTEND',ftem3, pcols, lchnk )

   ! Tendency for dry mass adjustment of q (valid for FV only)

   if (dycore_is('LR')) then
      if (diag_phys_tend_use_native_impl) then
         tmp_q     (:ncol,:pver) = (state%q(:ncol,:pver,       1) - tmp_q     (:ncol,:pver))*rtdt
         tmp_cldliq(:ncol,:pver) = (state%q(:ncol,:pver,ixcldliq) - tmp_cldliq(:ncol,:pver))*rtdt
         tmp_cldice(:ncol,:pver) = (state%q(:ncol,:pver,ixcldice) - tmp_cldice(:ncol,:pver))*rtdt
      else
         call diag_phys_tend_writeout_codon( &
              4_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
              int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
              real(ztodt, c_double), real(rtdt, c_double), 0._c_double, real(cpair, c_double), &
              c_null_ptr, c_loc(state%q), c_null_ptr, c_loc(tmp_t), c_loc(tmp_q), c_loc(tmp_cldliq), &
              c_loc(tmp_cldice), c_loc(qini), c_loc(cldliqini), c_loc(cldiceini), c_loc(ftem2), &
              c_loc(ftem3), c_null_ptr &
         )
      end if
      if ( cnst_cam_outfld(       1) ) call outfld (dmetendnam(       1), tmp_q     , pcols, lchnk)
      if ( cnst_cam_outfld(ixcldliq) ) call outfld (dmetendnam(ixcldliq), tmp_cldliq, pcols, lchnk)
      if ( cnst_cam_outfld(ixcldice) ) call outfld (dmetendnam(ixcldice), tmp_cldice, pcols, lchnk)
   end if

   ! Total physics tendency for moisture and other tracers

   if ( cnst_cam_outfld(       1) ) then
      if (diag_phys_tend_use_native_impl) then
         ftem3(:ncol,:pver) = (state%q(:ncol,:pver,       1) - qini     (:ncol,:pver) )*rtdt
      else
         call diag_phys_tend_writeout_codon( &
              5_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
              int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
              real(ztodt, c_double), real(rtdt, c_double), 0._c_double, real(cpair, c_double), &
              c_null_ptr, c_loc(state%q), c_null_ptr, c_loc(tmp_t), c_loc(tmp_q), c_loc(tmp_cldliq), &
              c_loc(tmp_cldice), c_loc(qini), c_loc(cldliqini), c_loc(cldiceini), c_loc(ftem2), &
              c_loc(ftem3), c_null_ptr &
         )
      end if
      call outfld (ptendnam(       1), ftem3, pcols, lchnk)
   end if
   if ( cnst_cam_outfld(ixcldliq) ) then
      if (diag_phys_tend_use_native_impl) then
         ftem3(:ncol,:pver) = (state%q(:ncol,:pver,ixcldliq) - cldliqini(:ncol,:pver) )*rtdt
      else
         call diag_phys_tend_writeout_codon( &
              6_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
              int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
              real(ztodt, c_double), real(rtdt, c_double), 0._c_double, real(cpair, c_double), &
              c_null_ptr, c_loc(state%q), c_null_ptr, c_loc(tmp_t), c_loc(tmp_q), c_loc(tmp_cldliq), &
              c_loc(tmp_cldice), c_loc(qini), c_loc(cldliqini), c_loc(cldiceini), c_loc(ftem2), &
              c_loc(ftem3), c_null_ptr &
         )
      end if
      call outfld (ptendnam(ixcldliq), ftem3, pcols, lchnk)
   end if
   if ( cnst_cam_outfld(ixcldice) ) then
      if (diag_phys_tend_use_native_impl) then
         ftem3(:ncol,:pver) = (state%q(:ncol,:pver,ixcldice) - cldiceini(:ncol,:pver) )*rtdt
      else
         call diag_phys_tend_writeout_codon( &
              7_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
              int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
              real(ztodt, c_double), real(rtdt, c_double), 0._c_double, real(cpair, c_double), &
              c_null_ptr, c_loc(state%q), c_null_ptr, c_loc(tmp_t), c_loc(tmp_q), c_loc(tmp_cldliq), &
              c_loc(tmp_cldice), c_loc(qini), c_loc(cldliqini), c_loc(cldiceini), c_loc(ftem2), &
              c_loc(ftem3), c_null_ptr &
         )
      end if
      call outfld (ptendnam(ixcldice), ftem3, pcols, lchnk)
   end if

   ! Total (physics+dynamics, everything!) tendency for Temperature

   !! get temperature stored in physics buffer
   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, t_ttend_idx, t_ttend, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   !! calculate and outfld the total temperature tendency
   if (diag_phys_tend_use_native_impl) then
      ftem3(:ncol,:) = (state%t(:ncol,:) - t_ttend(:ncol,:))/ztodt
   else
      call diag_phys_tend_writeout_codon( &
           8_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
           real(ztodt, c_double), real(rtdt, c_double), 0._c_double, real(cpair, c_double), &
           c_loc(state%t), c_loc(state%q), c_loc(tend%dtdt), c_loc(tmp_t), c_loc(tmp_q), &
           c_loc(tmp_cldliq), c_loc(tmp_cldice), c_loc(qini), c_loc(cldliqini), c_loc(cldiceini), &
           c_loc(ftem2), c_loc(ftem3), c_loc(t_ttend) &
      )
   end if
   call outfld('TTEND_TOT', ftem3, pcols, lchnk)

   !! update physics buffer with this time-step's temperature
   if (diag_phys_tend_use_native_impl) then
      t_ttend(:ncol,:) = state%t(:ncol,:)
   else
      call diag_phys_tend_writeout_codon( &
           9_c_int64_t, int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
           int(pcnst, c_int64_t), int(ixcldliq, c_int64_t), int(ixcldice, c_int64_t), &
           real(ztodt, c_double), real(rtdt, c_double), 0._c_double, real(cpair, c_double), &
           c_loc(state%t), c_loc(state%q), c_loc(tend%dtdt), c_loc(tmp_t), c_loc(tmp_q), &
           c_loc(tmp_cldliq), c_loc(tmp_cldice), c_loc(qini), c_loc(cldliqini), c_loc(cldiceini), &
           c_loc(ftem2), c_loc(ftem3), c_loc(t_ttend) &
      )
   end if

end subroutine diag_phys_tend_writeout

!#######################################################################

   subroutine diag_state_b4_phys_write (state)
!
!---------------------------------------------------------------
!
! Purpose:  Dump state just prior to executing physics
!
!---------------------------------------------------------------
!
! Arguments
!
   use iso_c_binding, only: c_int64_t

   type(physics_state), intent(in) :: state
!
!---------------------------Local workspace-----------------------------
!
   integer :: ixcldice, ixcldliq ! constituent indices for cloud liquid and ice water.
   integer :: lchnk              ! chunk index
!
!-----------------------------------------------------------------------
!
   lchnk = state%lchnk

   call cam_diag_touch_and_log(6_c_int64_t, diag_state_b4_phys_write_logged, &
        'diag_state_b4_phys_write direct = codon; state snapshot selector/touch direct = codon; ' // &
        'cnst_get_ind/outfld native CAM API island')

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call outfld('TBP', state%t, pcols, lchnk   )
   if ( cnst_cam_outfld(       1) ) call outfld (bpcnst(       1), state%q(1,1,       1), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldliq) ) call outfld (bpcnst(ixcldliq), state%q(1,1,ixcldliq), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldice) ) call outfld (bpcnst(ixcldice), state%q(1,1,ixcldice), pcols, lchnk)

   end subroutine diag_state_b4_phys_write

end module cam_diagnostics
