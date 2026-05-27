module phys_control
!-----------------------------------------------------------------------
! Purpose:
!
! Provides a control interface to CAM physics packages
!
! Revision history:
! 2006-05-01  D. B. Coleman,  Creation of module
! 2009-02-13  Eaton           Replace *_{default,set}opts methods with module namelist.
!                             Add vars to indicate physics version and chemistry type.
!-----------------------------------------------------------------------

use spmd_utils,     only: masterproc
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use shr_kind_mod,   only: r8 => shr_kind_r8
use iso_c_binding,  only: c_int64_t, c_loc

implicit none
private
save

public :: &
   phys_ctl_readnl,   &! read namelist from file
   phys_getopts,      &! generic query method
   phys_deepconv_pbl, &! return true if deep convection is allowed in the PBL
   phys_do_flux_avg,  &! return true to average surface fluxes
   cam_physpkg_is,    &! query for the name of the physics package
   cam_chempkg_is,    &! query for the name of the chemistry package
   waccmx_is

! Private module data

character(len=16), parameter :: unset_str = 'UNSET'
integer,           parameter :: unset_int = huge(1)

! Namelist variables:
character(len=16) :: cam_physpkg          = unset_str  ! CAM physics package [cam3 | cam4 | cam5 |
                                                       !   ideal | adiabatic].
character(len=32) :: cam_chempkg          = unset_str  ! CAM chemistry package 
character(len=16) :: waccmx_opt           = unset_str  ! WACCMX run option [ionosphere | neutral | off
character(len=16) :: deep_scheme          = unset_str  ! deep convection package
character(len=16) :: shallow_scheme       = unset_str  ! shallow convection package
character(len=16) :: eddy_scheme          = unset_str  ! vertical diffusion package
character(len=16) :: microp_scheme        = unset_str  ! microphysics package
character(len=16) :: macrop_scheme        = unset_str  ! macrophysics package
character(len=16) :: radiation_scheme     = unset_str  ! radiation package
integer           :: srf_flux_avg         = unset_int  ! 1 => smooth surface fluxes, 0 otherwise

logical           :: use_subcol_microp    = .false.    ! if .true. then use sub-columns in microphysics

logical           :: atm_dep_flux         = .true.     ! true => deposition fluxes will be provided
                                                       ! to the coupler
logical           :: history_amwg         = .true.     ! output the variables used by the AMWG diag package
logical           :: history_vdiag        = .false.    ! output the variables used by the AMWG variability diag package
logical           :: history_aerosol      = .false.    ! output the MAM aerosol variables and tendencies
logical           :: history_aero_optics  = .false.    ! output the aerosol
logical           :: history_eddy         = .false.    ! output the eddy variables
logical           :: history_budget       = .false.    ! output tendencies and state variables for CAM4
                                                       ! temperature, water vapor, cloud ice and cloud
                                                       ! liquid budgets.
integer           :: history_budget_histfile_num = 1   ! output history file number for budget fields
logical           :: history_waccm        = .false.    ! output variables of interest for WACCM runs
logical           :: history_waccmx       = .false.    ! output variables of interest for WACCM-X runs
logical           :: history_chemistry    = .true.     ! output default chemistry-related variables
logical           :: history_carma        = .true.     ! output default CARMA-related variables
logical           :: history_clubb        = .true.     ! output default CLUBB-related variables
logical           :: do_clubb_sgs
logical           :: do_tms
logical           :: micro_do_icesupersat
! Check validity of physics_state objects in physics_update.
logical           :: state_debug_checks   = .false.

! Macro/micro-physics co-substeps
integer           :: cld_macmic_num_steps = 1

logical           :: offline_driver       = .false.    ! true => offline driver is being used

logical :: prog_modal_aero ! determines whether prognostic modal aerosols are present in the run.

! Option to use heterogeneous freezing
logical, public, protected :: use_hetfrz_classnuc = .false.

! Which gravity wave sources are used?
! Orography.
logical, public, protected :: use_gw_oro = .true.
! Frontogenesis.
logical, public, protected :: use_gw_front = .false.
! Frontogenesis to inertial spectrum.
logical, public, protected :: use_gw_front_igw = .false.
! Deep convection.
logical, public, protected :: use_gw_convect_dp = .false.
! Shallow convection.
logical, public, protected :: use_gw_convect_sh = .false.

logical :: use_native_phys_control_bool_helpers_impl = .false.
logical :: phys_control_bool_helpers_impl_selected = .false.
logical :: phys_control_bool_helpers_proof_written = .false.
logical :: cam_physpkg_is_logged = .false.
logical :: cam_chempkg_is_logged = .false.
logical :: waccmx_is_logged = .false.
logical :: phys_ctl_readnl_logged = .false.
logical :: phys_getopts_logged = .false.
logical :: phys_deepconv_pbl_logged = .false.
logical :: phys_do_flux_avg_logged = .false.

interface
   function phys_control_deepconv_pbl_codon(eddy_diag_tke_c, shallow_uw_c) &
        result(deepconv_pbl_c) bind(c, name="phys_control_deepconv_pbl_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: eddy_diag_tke_c, shallow_uw_c
      integer(c_int64_t) :: deepconv_pbl_c
   end function phys_control_deepconv_pbl_codon

   function phys_control_do_flux_avg_codon(srf_flux_avg_c) &
        result(do_flux_avg_c) bind(c, name="phys_control_do_flux_avg_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: srf_flux_avg_c
      integer(c_int64_t) :: do_flux_avg_c
   end function phys_control_do_flux_avg_codon

   function phys_control_bool_flag_codon(flag_c) &
        result(flag_out_c) bind(c, name="phys_control_bool_flag_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: flag_out_c
   end function phys_control_bool_flag_codon

   function phys_control_index_positive_codon(index_c) &
        result(flag_out_c) bind(c, name="phys_control_index_positive_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: index_c
      integer(c_int64_t) :: flag_out_c
   end function phys_control_index_positive_codon

   function phys_control_int_value_codon(value_c) &
        result(value_out_c) bind(c, name="phys_control_int_value_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: value_c
      integer(c_int64_t) :: value_out_c
   end function phys_control_int_value_codon

   function cam_physpkg_is_codon(name_len_c, name_ascii_p, pkg_len_c, pkg_ascii_p) &
        result(match_c) bind(c, name="cam_physpkg_is_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: name_len_c, pkg_len_c
      type(c_ptr), value :: name_ascii_p, pkg_ascii_p
      integer(c_int64_t) :: match_c
   end function cam_physpkg_is_codon

   function cam_chempkg_is_codon(name_len_c, name_ascii_p, pkg_len_c, pkg_ascii_p) &
        result(match_c) bind(c, name="cam_chempkg_is_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: name_len_c, pkg_len_c
      type(c_ptr), value :: name_ascii_p, pkg_ascii_p
      integer(c_int64_t) :: match_c
   end function cam_chempkg_is_codon

   function waccmx_is_codon(name_len_c, name_ascii_p, opt_len_c, opt_ascii_p) &
        result(match_c) bind(c, name="waccmx_is_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: name_len_c, opt_len_c
      type(c_ptr), value :: name_ascii_p, opt_ascii_p
      integer(c_int64_t) :: match_c
   end function waccmx_is_codon
end interface

!======================================================================= 
contains
!======================================================================= 

subroutine phys_control_bool_helpers_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (phys_control_bool_helpers_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('PHYS_CONTROL_BOOL_HELPERS_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_phys_control_bool_helpers_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_phys_control_bool_helpers_impl = .false.
   end if

   phys_control_bool_helpers_impl_selected = .true.

   if (masterproc) then
      if (use_native_phys_control_bool_helpers_impl) then
         write(iulog,*) 'phys_control_bool_helpers implementation = native'
      else
         write(iulog,*) 'phys_control_bool_helpers implementation = codon'
      end if
   end if

end subroutine phys_control_bool_helpers_select_impl

!===============================================================================

subroutine phys_control_bool_helpers_proof_once()

   if (phys_control_bool_helpers_proof_written) return
   phys_control_bool_helpers_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'phys_control_bool_helpers entered (runtime option boolean helpers = codon)'
   end if

end subroutine phys_control_bool_helpers_proof_once

!===============================================================================

subroutine phys_control_log_direct(logged, proof_line)

   logical, intent(inout) :: logged
   character(len=*), intent(in) :: proof_line

   if (logged) return
   logged = .true.

   if (masterproc) then
      write(iulog,'(A)') trim(proof_line)
      call flush(iulog)
   end if

end subroutine phys_control_log_direct

!===============================================================================

subroutine phys_control_ascii_pack(text, ascii)

   character(len=*), intent(in) :: text
   integer(c_int64_t), intent(out) :: ascii(:)
   integer :: i

   do i = 1, len(text)
      ascii(i) = int(iachar(text(i:i)), c_int64_t)
   end do

end subroutine phys_control_ascii_pack

subroutine phys_ctl_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'phys_ctl_readnl'

   namelist /phys_ctl_nl/ cam_physpkg, cam_chempkg, waccmx_opt, deep_scheme, shallow_scheme, &
      eddy_scheme, microp_scheme,  macrop_scheme, radiation_scheme, srf_flux_avg, &
      use_subcol_microp, atm_dep_flux, history_amwg, history_vdiag, history_aerosol, history_aero_optics, &
      history_eddy, history_budget,  history_budget_histfile_num, history_waccm, &
      history_waccmx, history_chemistry, history_carma, history_clubb, &
      do_clubb_sgs, do_tms, state_debug_checks, use_hetfrz_classnuc, use_gw_oro, use_gw_front, &
      use_gw_front_igw, use_gw_convect_dp, use_gw_convect_sh, cld_macmic_num_steps, &
      offline_driver, micro_do_icesupersat
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'phys_ctl_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, phys_ctl_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(deep_scheme,      len(deep_scheme)      , mpichar, 0, mpicom)
   call mpibcast(cam_physpkg,      len(cam_physpkg)      , mpichar, 0, mpicom)
   call mpibcast(cam_chempkg,      len(cam_chempkg)      , mpichar, 0, mpicom)
   call mpibcast(waccmx_opt,       len(waccmx_opt)       , mpichar, 0, mpicom)
   call mpibcast(shallow_scheme,   len(shallow_scheme)   , mpichar, 0, mpicom)
   call mpibcast(eddy_scheme,      len(eddy_scheme)      , mpichar, 0, mpicom)
   call mpibcast(microp_scheme,    len(microp_scheme)    , mpichar, 0, mpicom)
   call mpibcast(radiation_scheme, len(radiation_scheme) , mpichar, 0, mpicom)
   call mpibcast(macrop_scheme,    len(macrop_scheme)    , mpichar, 0, mpicom)
   call mpibcast(srf_flux_avg,                    1 , mpiint,  0, mpicom)
   call mpibcast(use_subcol_microp,               1 , mpilog,  0, mpicom)
   call mpibcast(atm_dep_flux,                    1 , mpilog,  0, mpicom)
   call mpibcast(history_amwg,                    1 , mpilog,  0, mpicom)
   call mpibcast(history_vdiag,                   1 , mpilog,  0, mpicom)
   call mpibcast(history_eddy,                    1 , mpilog,  0, mpicom)
   call mpibcast(history_aerosol,                 1 , mpilog,  0, mpicom)
   call mpibcast(history_aero_optics,             1 , mpilog,  0, mpicom)
   call mpibcast(history_budget,                  1 , mpilog,  0, mpicom)
   call mpibcast(history_budget_histfile_num,     1 , mpiint,  0, mpicom)
   call mpibcast(history_waccm,                   1 , mpilog,  0, mpicom)
   call mpibcast(history_waccmx,                  1 , mpilog,  0, mpicom)
   call mpibcast(history_chemistry,               1 , mpilog,  0, mpicom)
   call mpibcast(history_carma,                   1 , mpilog,  0, mpicom)
   call mpibcast(history_clubb,                   1 , mpilog,  0, mpicom)
   call mpibcast(do_clubb_sgs,                    1 , mpilog,  0, mpicom)
   call mpibcast(do_tms,                          1 , mpilog,  0, mpicom)
   call mpibcast(micro_do_icesupersat,            1 , mpilog,  0, mpicom)
   call mpibcast(state_debug_checks,              1 , mpilog,  0, mpicom)
   call mpibcast(use_hetfrz_classnuc,             1 , mpilog,  0, mpicom)
   call mpibcast(use_gw_oro,                      1 , mpilog,  0, mpicom)
   call mpibcast(use_gw_front,                    1 , mpilog,  0, mpicom)
   call mpibcast(use_gw_front_igw,                1 , mpilog,  0, mpicom)
   call mpibcast(use_gw_convect_dp,               1 , mpilog,  0, mpicom)
   call mpibcast(use_gw_convect_sh,               1 , mpilog,  0, mpicom)
   call mpibcast(cld_macmic_num_steps,            1 , mpiint,  0, mpicom)
   call mpibcast(offline_driver,                  1 , mpilog,  0, mpicom)
#endif

   ! Error checking:

   ! Check compatibility of eddy & shallow schemes
   if (( shallow_scheme .eq. 'UW' ) .and. ( eddy_scheme .ne. 'diag_TKE' )) then
      write(iulog,*)'Do you really want to run UW shallow scheme without diagnostic TKE eddy scheme? Quiting'
      call endrun('shallow convection and eddy scheme may be incompatible')
   endif

   if (( shallow_scheme .eq. 'Hack' ) .and. ( ( eddy_scheme .ne. 'HB' ) .and. ( eddy_scheme .ne. 'HBR' ))) then
      write(iulog,*)'Do you really want to run Hack shallow scheme with a non-standard eddy scheme? Quiting.'
      call endrun('shallow convection and eddy scheme may be incompatible')
   endif

   ! Check compatibility of PBL and Microphysics schemes
   if (( eddy_scheme .eq. 'diag_TKE' ) .and. ( microp_scheme .eq. 'RK' )) then
      write(iulog,*)'UW PBL is not compatible with RK microphysics.  Quiting'
      call endrun('PBL and Microphysics schemes incompatible')
   endif
   
   ! Add a check to make sure CLUBB and MG are used together
   if ( do_clubb_sgs .and. ( microp_scheme .ne. 'MG')) then
      write(iulog,*)'CLUBB is only compatible with MG microphysics.  Quiting'
      call endrun('CLUBB and microphysics schemes incompatible')
   endif

   ! Check that eddy_scheme, macrop_scheme, shallow_scheme are all set to CLUBB_SGS if do_clubb_sgs is true
   if (do_clubb_sgs) then
      if (eddy_scheme .ne. 'CLUBB_SGS' .or. macrop_scheme .ne. 'CLUBB_SGS' .or. shallow_scheme .ne. 'CLUBB_SGS') then
         write(iulog,*)'eddy_scheme, macrop_scheme and shallow_scheme must all be CLUBB_SGS.  Quiting'
         call endrun('CLUBB and eddy, macrop or shallow schemes incompatible')
      endif
   endif
      
   ! Macro/micro co-substepping support.
   if (cld_macmic_num_steps > 1) then
      if (microp_scheme /= "MG" .or. (macrop_scheme /= "park" .and. macrop_scheme /= "CLUBB_SGS")) then
         call endrun ("Setting cld_macmic_num_steps > 1 is only &
              &supported with Park or CLUBB macrophysics and MG microphysics.")
      end if
   end if

   ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
   call phys_control_bool_helpers_select_impl()
   if (use_native_phys_control_bool_helpers_impl) then
      prog_modal_aero = index(cam_chempkg,'_mam')>0
   else
      call phys_control_bool_helpers_proof_once()
      prog_modal_aero = (phys_control_index_positive_codon(int(index(cam_chempkg,'_mam'), c_int64_t)) /= 0_c_int64_t)
      call phys_control_log_direct(phys_ctl_readnl_logged, &
           'phys_ctl_readnl direct = codon; namelist/MPI/compatibility native islands')
   end if

end subroutine phys_ctl_readnl

!===============================================================================

logical function cam_physpkg_is(name)

   ! query for the name of the physics package

   character(len=*) :: name
   integer(c_int64_t), target :: name_ascii(len(name))
   integer(c_int64_t), target :: pkg_ascii(len(cam_physpkg))
   
   call phys_control_bool_helpers_select_impl()
   if (use_native_phys_control_bool_helpers_impl) then
      cam_physpkg_is = (trim(name) == trim(cam_physpkg))
      return
   end if

   call phys_control_bool_helpers_proof_once()
   call phys_control_ascii_pack(name, name_ascii)
   call phys_control_ascii_pack(cam_physpkg, pkg_ascii)
   cam_physpkg_is = (cam_physpkg_is_codon(int(len(name), c_int64_t), c_loc(name_ascii(1)), &
        int(len(cam_physpkg), c_int64_t), c_loc(pkg_ascii(1))) /= 0_c_int64_t)
   call phys_control_log_direct(cam_physpkg_is_logged, 'cam_physpkg_is direct = codon')
end function cam_physpkg_is

!===============================================================================

logical function cam_chempkg_is(name)

   ! query for the name of the chemics package

   character(len=*) :: name
   integer(c_int64_t), target :: name_ascii(len(name))
   integer(c_int64_t), target :: pkg_ascii(len(cam_chempkg))
   
   call phys_control_bool_helpers_select_impl()
   if (use_native_phys_control_bool_helpers_impl) then
      cam_chempkg_is = (trim(name) == trim(cam_chempkg))
      return
   end if

   call phys_control_bool_helpers_proof_once()
   call phys_control_ascii_pack(name, name_ascii)
   call phys_control_ascii_pack(cam_chempkg, pkg_ascii)
   cam_chempkg_is = (cam_chempkg_is_codon(int(len(name), c_int64_t), c_loc(name_ascii(1)), &
        int(len(cam_chempkg), c_int64_t), c_loc(pkg_ascii(1))) /= 0_c_int64_t)
   call phys_control_log_direct(cam_chempkg_is_logged, 'cam_chempkg_is direct = codon')
end function cam_chempkg_is

!===============================================================================

logical function waccmx_is(name)

   ! query for the name of the waccmx run option

   character(len=*) :: name
   integer(c_int64_t), target :: name_ascii(len(name))
   integer(c_int64_t), target :: opt_ascii(len(waccmx_opt))
   
   call phys_control_bool_helpers_select_impl()
   if (use_native_phys_control_bool_helpers_impl) then
      waccmx_is = (trim(name) == trim(waccmx_opt))
      return
   end if

   call phys_control_bool_helpers_proof_once()
   call phys_control_ascii_pack(name, name_ascii)
   call phys_control_ascii_pack(waccmx_opt, opt_ascii)
   waccmx_is = (waccmx_is_codon(int(len(name), c_int64_t), c_loc(name_ascii(1)), &
        int(len(waccmx_opt), c_int64_t), c_loc(opt_ascii(1))) /= 0_c_int64_t)
   call phys_control_log_direct(waccmx_is_logged, 'waccmx_is direct = codon')
end function waccmx_is

!===============================================================================

subroutine phys_getopts(deep_scheme_out, shallow_scheme_out, eddy_scheme_out, microp_scheme_out, &
                        radiation_scheme_out, use_subcol_microp_out, atm_dep_flux_out, &
                         history_amwg_out, history_vdiag_out, history_aerosol_out, history_aero_optics_out, history_eddy_out, &
                        history_budget_out, history_budget_histfile_num_out, &
                        history_waccm_out, history_waccmx_out, history_chemistry_out, &
                        history_carma_out, history_clubb_out, &
                        cam_chempkg_out, prog_modal_aero_out, macrop_scheme_out, &
                        do_clubb_sgs_out, do_tms_out, state_debug_checks_out, cld_macmic_num_steps_out, &
                        offline_driver_out, micro_do_icesupersat_out)
!-----------------------------------------------------------------------
! Purpose: Return runtime settings
!          deep_scheme_out   : deep convection scheme
!          shallow_scheme_out: shallow convection scheme
!          eddy_scheme_out   : vertical diffusion scheme
!          microp_scheme_out : microphysics scheme
!          radiation_scheme_out : radiation_scheme
!-----------------------------------------------------------------------

   character(len=16), intent(out), optional :: deep_scheme_out
   character(len=16), intent(out), optional :: shallow_scheme_out
   character(len=16), intent(out), optional :: eddy_scheme_out
   character(len=16), intent(out), optional :: microp_scheme_out
   character(len=16), intent(out), optional :: radiation_scheme_out
   character(len=16), intent(out), optional :: macrop_scheme_out
   logical,           intent(out), optional :: use_subcol_microp_out
   logical,           intent(out), optional :: atm_dep_flux_out
   logical,           intent(out), optional :: history_amwg_out
   logical,           intent(out), optional :: history_vdiag_out
   logical,           intent(out), optional :: history_eddy_out
   logical,           intent(out), optional :: history_aerosol_out
   logical,           intent(out), optional :: history_aero_optics_out
   logical,           intent(out), optional :: history_budget_out
   integer,           intent(out), optional :: history_budget_histfile_num_out
   logical,           intent(out), optional :: history_waccm_out
   logical,           intent(out), optional :: history_waccmx_out
   logical,           intent(out), optional :: history_chemistry_out
   logical,           intent(out), optional :: history_carma_out
   logical,           intent(out), optional :: history_clubb_out
   logical,           intent(out), optional :: do_clubb_sgs_out
   logical,           intent(out), optional :: micro_do_icesupersat_out        
   character(len=32), intent(out), optional :: cam_chempkg_out
   logical,           intent(out), optional :: prog_modal_aero_out
   logical,           intent(out), optional :: do_tms_out
   logical,           intent(out), optional :: state_debug_checks_out
   integer,           intent(out), optional :: cld_macmic_num_steps_out
   logical,           intent(out), optional :: offline_driver_out

   call phys_control_bool_helpers_select_impl()

   if ( present(deep_scheme_out         ) ) deep_scheme_out          = deep_scheme
   if ( present(shallow_scheme_out      ) ) shallow_scheme_out       = shallow_scheme
   if ( present(eddy_scheme_out         ) ) eddy_scheme_out          = eddy_scheme
   if ( present(microp_scheme_out       ) ) microp_scheme_out        = microp_scheme
   if ( present(radiation_scheme_out    ) ) radiation_scheme_out     = radiation_scheme

   if ( present(use_subcol_microp_out   ) ) use_subcol_microp_out    = use_subcol_microp
   if ( present(macrop_scheme_out       ) ) macrop_scheme_out        = macrop_scheme
   if ( present(atm_dep_flux_out        ) ) atm_dep_flux_out         = atm_dep_flux
   if ( present(history_aerosol_out     ) ) history_aerosol_out      = history_aerosol
   if ( present(history_aero_optics_out ) ) history_aero_optics_out  = history_aero_optics
   if ( present(history_budget_out      ) ) history_budget_out       = history_budget
   if ( present(history_amwg_out        ) ) history_amwg_out         = history_amwg
   if ( present(history_vdiag_out       ) ) history_vdiag_out        = history_vdiag
   if ( present(history_eddy_out        ) ) history_eddy_out         = history_eddy
   if ( present(history_budget_histfile_num_out ) ) history_budget_histfile_num_out = history_budget_histfile_num
   if ( present(history_waccm_out       ) ) history_waccm_out        = history_waccm
   if ( present(history_waccmx_out      ) ) history_waccmx_out       = history_waccmx
   if ( present(history_chemistry_out   ) ) history_chemistry_out    = history_chemistry
   if ( present(history_carma_out       ) ) history_carma_out        = history_carma
   if ( present(history_clubb_out       ) ) history_clubb_out        = history_clubb
   if ( present(do_clubb_sgs_out        ) ) do_clubb_sgs_out         = do_clubb_sgs
   if ( present(micro_do_icesupersat_out )) micro_do_icesupersat_out = micro_do_icesupersat
   if ( present(cam_chempkg_out         ) ) cam_chempkg_out          = cam_chempkg
   if ( present(prog_modal_aero_out     ) ) prog_modal_aero_out      = prog_modal_aero
   if ( present(do_tms_out              ) ) do_tms_out               = do_tms
   if ( present(state_debug_checks_out  ) ) state_debug_checks_out   = state_debug_checks
   if ( present(cld_macmic_num_steps_out) ) cld_macmic_num_steps_out = cld_macmic_num_steps
   if ( present(offline_driver_out      ) ) offline_driver_out       = offline_driver

   if (.not. use_native_phys_control_bool_helpers_impl) then
      call phys_control_bool_helpers_proof_once()
      call phys_control_log_direct(phys_getopts_logged, &
           'phys_getopts direct = codon; optional output assignment native boundary')
      if (present(use_subcol_microp_out)) &
         use_subcol_microp_out = (phys_control_bool_flag_codon( &
         merge(1_c_int64_t, 0_c_int64_t, use_subcol_microp_out)) /= 0_c_int64_t)
      if (present(atm_dep_flux_out)) &
         atm_dep_flux_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, atm_dep_flux_out)) /= 0_c_int64_t)
      if (present(history_amwg_out)) &
         history_amwg_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, history_amwg_out)) /= 0_c_int64_t)
      if (present(history_vdiag_out)) &
         history_vdiag_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, history_vdiag_out)) /= 0_c_int64_t)
      if (present(history_aerosol_out)) &
         history_aerosol_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, history_aerosol_out)) /= 0_c_int64_t)
      if (present(history_aero_optics_out)) &
         history_aero_optics_out = (phys_control_bool_flag_codon( &
         merge(1_c_int64_t, 0_c_int64_t, history_aero_optics_out)) /= 0_c_int64_t)
      if (present(history_eddy_out)) &
         history_eddy_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, history_eddy_out)) /= 0_c_int64_t)
      if (present(history_budget_out)) &
         history_budget_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, history_budget_out)) /= 0_c_int64_t)
      if (present(history_waccm_out)) &
         history_waccm_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, history_waccm_out)) /= 0_c_int64_t)
      if (present(history_waccmx_out)) &
         history_waccmx_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, history_waccmx_out)) /= 0_c_int64_t)
      if (present(history_chemistry_out)) &
         history_chemistry_out = (phys_control_bool_flag_codon( &
         merge(1_c_int64_t, 0_c_int64_t, history_chemistry_out)) /= 0_c_int64_t)
      if (present(history_carma_out)) &
         history_carma_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, history_carma_out)) /= 0_c_int64_t)
      if (present(history_clubb_out)) &
         history_clubb_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, history_clubb_out)) /= 0_c_int64_t)
      if (present(do_clubb_sgs_out)) &
         do_clubb_sgs_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, do_clubb_sgs_out)) /= 0_c_int64_t)
      if (present(micro_do_icesupersat_out)) &
         micro_do_icesupersat_out = (phys_control_bool_flag_codon( &
         merge(1_c_int64_t, 0_c_int64_t, micro_do_icesupersat_out)) /= 0_c_int64_t)
      if (present(prog_modal_aero_out)) &
         prog_modal_aero_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, prog_modal_aero_out)) /= 0_c_int64_t)
      if (present(do_tms_out)) &
         do_tms_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, do_tms_out)) /= 0_c_int64_t)
      if (present(state_debug_checks_out)) &
         state_debug_checks_out = (phys_control_bool_flag_codon( &
         merge(1_c_int64_t, 0_c_int64_t, state_debug_checks_out)) /= 0_c_int64_t)
      if (present(offline_driver_out)) &
         offline_driver_out = (phys_control_bool_flag_codon(merge(1_c_int64_t, 0_c_int64_t, offline_driver_out)) /= 0_c_int64_t)
      if (present(history_budget_histfile_num_out)) &
         history_budget_histfile_num_out = int(phys_control_int_value_codon(int(history_budget_histfile_num_out, c_int64_t)))
      if (present(cld_macmic_num_steps_out)) &
         cld_macmic_num_steps_out = int(phys_control_int_value_codon(int(cld_macmic_num_steps_out, c_int64_t)))
   end if

end subroutine phys_getopts

!===============================================================================

function phys_deepconv_pbl()

  logical phys_deepconv_pbl
  integer(c_int64_t) :: eddy_diag_tke_c, shallow_uw_c, result_c

   call phys_control_bool_helpers_select_impl()
   if (use_native_phys_control_bool_helpers_impl) then
      phys_deepconv_pbl = phys_deepconv_pbl_native()
      return
   end if

   call phys_control_bool_helpers_proof_once()
   eddy_diag_tke_c = 0_c_int64_t
   shallow_uw_c = 0_c_int64_t
   if (eddy_scheme .eq. 'diag_TKE') eddy_diag_tke_c = 1_c_int64_t
   if (shallow_scheme .eq. 'UW') shallow_uw_c = 1_c_int64_t
   result_c = phys_control_deepconv_pbl_codon(eddy_diag_tke_c, shallow_uw_c)
   phys_deepconv_pbl = (result_c /= 0_c_int64_t)
   call phys_control_log_direct(phys_deepconv_pbl_logged, 'phys_deepconv_pbl direct = codon')

   return

end function phys_deepconv_pbl

!===============================================================================

function phys_deepconv_pbl_native()

  logical phys_deepconv_pbl_native

   ! Don't allow deep convection in PBL if running UW PBL scheme
   if ( (eddy_scheme .eq. 'diag_TKE' ) .or. (shallow_scheme .eq. 'UW' ) ) then
      phys_deepconv_pbl_native = .true.
   else
      phys_deepconv_pbl_native = .false.
   endif

   return

end function phys_deepconv_pbl_native

!===============================================================================

function phys_do_flux_avg()

   logical :: phys_do_flux_avg
   integer(c_int64_t) :: result_c
   !----------------------------------------------------------------------

   call phys_control_bool_helpers_select_impl()
   if (use_native_phys_control_bool_helpers_impl) then
      phys_do_flux_avg = phys_do_flux_avg_native()
      return
   end if

   call phys_control_bool_helpers_proof_once()
   result_c = phys_control_do_flux_avg_codon(int(srf_flux_avg, c_int64_t))
   phys_do_flux_avg = (result_c /= 0_c_int64_t)
   call phys_control_log_direct(phys_do_flux_avg_logged, 'phys_do_flux_avg direct = codon')

end function phys_do_flux_avg

!===============================================================================

function phys_do_flux_avg_native()

   logical :: phys_do_flux_avg_native
   !----------------------------------------------------------------------

   phys_do_flux_avg_native = .false.
   if (srf_flux_avg == 1) phys_do_flux_avg_native = .true.

end function phys_do_flux_avg_native

!===============================================================================
end module phys_control
