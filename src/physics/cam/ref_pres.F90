module ref_pres
!--------------------------------------------------------------------------
! 
! Provides access to reference pressures for use by the physics
! parameterizations.  The pressures are provided by the dynamical core
! since it currently determines the grid used by the physics.
! 
! Note that the init method for this module is called before the init
! method in physpkg; therefore, most physics modules can use these
! reference pressures during their init phases.
! 
!--------------------------------------------------------------------------

use shr_kind_mod, only: r8=>shr_kind_r8
use ppgrid,       only: pver, pverp

implicit none
public
save

! Reference pressures (Pa)
real(r8), protected, target :: pref_edge(pverp)     ! Layer edges
real(r8), protected, target :: pref_mid(pver)       ! Layer midpoints
real(r8), protected, target :: pref_mid_norm(pver)  ! Layer midpoints normalized by
                                            ! surface pressure ('eta' coordinate)

real(r8), protected :: ptop_ref             ! Top of model
real(r8), protected :: psurf_ref            ! Surface pressure

! Number of top levels using pure pressure representation
integer, protected :: num_pr_lev

! Pressure used to set troposphere cloud physics top (Pa)
real(r8), protected :: trop_cloud_top_press = 0._r8
! Top level for troposphere cloud physics
integer, protected :: trop_cloud_top_lev

! Pressure used to set MAM process top (Pa)
real(r8), protected :: clim_modal_aero_top_press = 0._r8
! Top level for MAM processes that impact climate
integer, protected :: clim_modal_aero_top_lev

! Molecular diffusion is calculated only if the model top is below this
! pressure (Pa).
real(r8), protected :: do_molec_press = 0.1_r8
! Pressure used to set bottom of molecular diffusion region (Pa).
real(r8), protected :: molec_diff_bot_press = 50._r8
! Flag for molecular diffusion, and molecular diffusion level indices.
logical, protected :: do_molec_diff = .false.
integer, protected :: ntop_molec = 1
integer, protected :: nbot_molec = 0

logical :: use_native_ref_pres_init_impl = .false.
logical :: ref_pres_init_impl_selected = .false.
logical :: ref_pres_init_proof_written = .false.

interface
   subroutine ref_pres_init_finalize_codon(pver_c, pverp_c, trop_cloud_top_press_c, &
        clim_modal_aero_top_press_c, do_molec_press_c, molec_diff_bot_press_c, pref_edge_p, &
        pref_mid_p, pref_mid_norm_p, scalar_out_p, int_out_p, flag_out_p) &
        bind(c, name="ref_pres_init_finalize_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: pver_c, pverp_c
      real(c_double), value :: trop_cloud_top_press_c, clim_modal_aero_top_press_c
      real(c_double), value :: do_molec_press_c, molec_diff_bot_press_c
      type(c_ptr), value :: pref_edge_p, pref_mid_p, pref_mid_norm_p
      type(c_ptr), value :: scalar_out_p, int_out_p, flag_out_p
   end subroutine ref_pres_init_finalize_codon
   function press_lim_idx_codon(p_c, top_c, pver_c, pref_mid_p) result(k_lim_c) &
        bind(c, name="press_lim_idx_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      real(c_double), value :: p_c
      integer(c_int64_t), value :: top_c, pver_c
      type(c_ptr), value :: pref_mid_p
      integer(c_int64_t) :: k_lim_c
   end function press_lim_idx_codon
end interface

!====================================================================================
contains
!====================================================================================

subroutine ref_pres_init_select_impl()

   use cam_logfile,  only: iulog
   use spmd_utils,   only: masterproc

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (ref_pres_init_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('REF_PRES_INIT_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_ref_pres_init_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_ref_pres_init_impl = .false.
   end if

   ref_pres_init_impl_selected = .true.

   if (masterproc) then
      if (use_native_ref_pres_init_impl) then
         write(iulog,*) 'ref_pres_init implementation = native'
      else
         write(iulog,*) 'ref_pres_init implementation = codon'
      end if
   end if

end subroutine ref_pres_init_select_impl

!====================================================================================

subroutine ref_pres_init_proof_once()

   use cam_logfile,  only: iulog
   use spmd_utils,   only: masterproc

   if (ref_pres_init_proof_written) return
   ref_pres_init_proof_written = .true.

   if (masterproc) then
      write(iulog,'(A)') 'ref_pres_init entered (pressure normalization/limit helpers = codon)'
   end if

end subroutine ref_pres_init_proof_once

!====================================================================================

subroutine ref_pres_readnl(nlfile)

   use spmd_utils,      only: masterproc
   use cam_abortutils,  only: endrun
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'ref_pres_readnl'

   namelist /ref_pres_nl/ trop_cloud_top_press, clim_modal_aero_top_press,&
        do_molec_press, molec_diff_bot_press
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'ref_pres_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, ref_pres_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! Check that top for modal aerosols is not lower than
      ! top for clouds.
      if (clim_modal_aero_top_press > trop_cloud_top_press) &
           call endrun("ERROR: clim_modal_aero_top press must be less &
           &than or equal to trop_cloud_top_press.")
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(trop_cloud_top_press,            1 , mpir8,   0, mpicom)
   call mpibcast(clim_modal_aero_top_press,       1 , mpir8,   0, mpicom)
   call mpibcast(do_molec_press,                  1 , mpir8,   0, mpicom)
   call mpibcast(molec_diff_bot_press,            1 , mpir8,   0, mpicom)
#endif

end subroutine ref_pres_readnl

!====================================================================================

subroutine ref_pres_init

  use iso_c_binding, only: c_double, c_int64_t, c_loc
  use dyn_grid,     only: dyn_grid_get_pref

  real(c_double), target :: scalar_out(2)
  integer(c_int64_t), target :: int_out(3), flag_out(1)

  ! Get reference pressures from the dynamical core.
  call dyn_grid_get_pref(pref_edge, pref_mid, num_pr_lev)

  call ref_pres_init_select_impl()
  if (use_native_ref_pres_init_impl) then
     ptop_ref = pref_edge(1)
     psurf_ref = pref_edge(pverp)

     pref_mid_norm = pref_mid/psurf_ref

     ! Find level corresponding to the top of troposphere clouds.
     trop_cloud_top_lev = press_lim_idx(trop_cloud_top_press, &
          top=.true.)

     ! Find level corresponding to the top for MAM processes.
     clim_modal_aero_top_lev = press_lim_idx(clim_modal_aero_top_press, &
          top=.true.)

     ! Find level corresponding to the molecular diffusion bottom.
     do_molec_diff = (ptop_ref < do_molec_press)
     if (do_molec_diff) then
        nbot_molec = press_lim_idx(molec_diff_bot_press, &
             top=.false.)
     end if
  else
     call ref_pres_init_proof_once()
     call ref_pres_init_finalize_codon(int(pver, c_int64_t), int(pverp, c_int64_t), &
          real(trop_cloud_top_press, c_double), real(clim_modal_aero_top_press, c_double), &
          real(do_molec_press, c_double), real(molec_diff_bot_press, c_double), &
          c_loc(pref_edge(1)), c_loc(pref_mid(1)), c_loc(pref_mid_norm(1)), &
          c_loc(scalar_out(1)), c_loc(int_out(1)), c_loc(flag_out(1)))
     ptop_ref = scalar_out(1)
     psurf_ref = scalar_out(2)
     trop_cloud_top_lev = int(int_out(1))
     clim_modal_aero_top_lev = int(int_out(2))
     nbot_molec = int(int_out(3))
     do_molec_diff = flag_out(1) /= 0_c_int64_t
  end if

end subroutine ref_pres_init

!====================================================================================

! Convert pressure limiters to the appropriate level.
function press_lim_idx(p, top) result(k_lim)

  use iso_c_binding, only: c_double, c_int64_t, c_loc

  ! Pressure
  real(r8), intent(in) :: p
  ! Is this a top or bottom limit?
  logical,  intent(in) :: top
  integer :: k_lim

  k_lim = int(press_lim_idx_codon(real(p, c_double), merge(1_c_int64_t, 0_c_int64_t, top), &
       int(pver, c_int64_t), c_loc(pref_mid(1))))

end function press_lim_idx

!====================================================================================

end module ref_pres
