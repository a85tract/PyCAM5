module microp_driver

!-------------------------------------------------------------------------------------------------------
!
! Driver for CAM microphysics parameterizations
!
!-------------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use ppgrid,         only: pver
use physics_types,  only: physics_state, physics_ptend, physics_tend,  &
                          physics_ptend_copy, physics_ptend_sum
use physics_buffer, only: pbuf_get_index, pbuf_get_field, physics_buffer_desc
use phys_control,   only: phys_getopts

use micro_mg_cam,   only: micro_mg_cam_readnl, micro_mg_cam_register, &
                          micro_mg_cam_implements_cnst, micro_mg_cam_init_cnst, &
                          micro_mg_cam_init, micro_mg_cam_tend
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use perf_mod,       only: t_startf, t_stopf
use spmd_utils,     only: masterproc

implicit none
private
save

public :: &
   microp_driver_readnl,          &
   microp_driver_register,        &
   microp_driver_init_cnst,       &
   microp_driver_implements_cnst, &
   microp_driver_init,            &
   microp_driver_tend

character(len=16)  :: microp_scheme   ! Microphysics scheme
logical :: use_native_impl = .false.
logical :: impl_selected = .false.
integer :: codon_scheme_code = 0
logical :: codon_scheme_selected = .false.
logical :: microp_driver_implements_cnst_logged = .false.
logical :: microp_driver_readnl_logged = .false.
logical :: microp_driver_register_logged = .false.
logical :: microp_driver_init_logged = .false.
logical :: microp_driver_tend_logged = .false.

interface
   function microp_driver_readnl_codon(flag_c) result(out_c) bind(c, name="microp_driver_readnl_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function microp_driver_readnl_codon
   function microp_driver_register_codon(flag_c) result(out_c) bind(c, name="microp_driver_register_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function microp_driver_register_codon
   function microp_driver_init_codon(flag_c) result(out_c) bind(c, name="microp_driver_init_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function microp_driver_init_codon
   function microp_driver_tend_codon(flag_c) result(out_c) bind(c, name="microp_driver_tend_codon")
      use iso_c_binding, only: c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function microp_driver_tend_codon
   function microp_driver_implements_cnst_codon(scheme_len_c, scheme_ascii_p, name_len_c, name_ascii_p) result(out_c) &
        bind(c, name="microp_driver_implements_cnst_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: scheme_len_c, name_len_c
      type(c_ptr), value :: scheme_ascii_p, name_ascii_p
      integer(c_int64_t) :: out_c
   end function microp_driver_implements_cnst_codon
end interface

!===============================================================================
contains
!===============================================================================

subroutine microp_driver_readnl(nlfile)

   use iso_c_binding, only: c_int64_t

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
   integer(c_int64_t) :: active_c

   ! Read in namelist for microphysics scheme
   !-----------------------------------------------------------------------

   call phys_getopts(microp_scheme_out=microp_scheme)
   active_c = microp_driver_readnl_codon(1_c_int64_t)
   if (active_c == 0_c_int64_t) return
   call microp_driver_log_direct(microp_driver_readnl_logged, &
        'microp_driver_readnl direct = codon; scheme dispatch direct = codon; micro_mg_cam_readnl native/Codon child boundary')

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_readnl(nlfile)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_readnl:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_readnl

subroutine microp_driver_register

   ! Register microphysics constituents and fields in the physics buffer.
   !-----------------------------------------------------------------------

   use iso_c_binding, only: c_int64_t

   integer(c_int64_t) :: active_c

   active_c = microp_driver_register_codon(1_c_int64_t)
   if (active_c == 0_c_int64_t) return
   call microp_driver_log_direct(microp_driver_register_logged, &
        'microp_driver_register direct = codon; scheme dispatch direct = codon; micro_mg_cam_register child boundary')

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_register()
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_register:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_register

!===============================================================================

function microp_driver_implements_cnst(name)

   use iso_c_binding, only: c_int64_t, c_loc

   ! Return true if specified constituent is implemented by the
   ! microphysics package

   character(len=*), intent(in) :: name        ! constituent name
   logical :: microp_driver_implements_cnst    ! return value

   ! Local workspace
   integer :: i
   integer(c_int64_t) :: out_c
   integer(c_int64_t), target :: name_ascii(max(1, len(name)))
   integer(c_int64_t), target :: scheme_ascii(len(microp_scheme))
   !-----------------------------------------------------------------------

   microp_driver_implements_cnst = .false.

   call microp_driver_select_impl()
   if (.not. use_native_impl) then
      select case (microp_scheme)
      case ('MG', 'RK')
         continue
      case default
         call endrun('microp_driver_implements_cnst:: unrecognized microp_scheme')
      end select
      do i = 1, len(name)
         name_ascii(i) = int(iachar(name(i:i)), c_int64_t)
      end do
      do i = 1, len(microp_scheme)
         scheme_ascii(i) = int(iachar(microp_scheme(i:i)), c_int64_t)
      end do
      out_c = microp_driver_implements_cnst_codon(int(len(microp_scheme), c_int64_t), &
           c_loc(scheme_ascii(1)), int(len(name), c_int64_t), c_loc(name_ascii(1)))
      microp_driver_implements_cnst = out_c /= 0_c_int64_t
      call microp_driver_log_direct(microp_driver_implements_cnst_logged, &
           'microp_driver_implements_cnst direct = codon')
      return
   end if

   select case (microp_scheme)
   case ('MG')
      microp_driver_implements_cnst = micro_mg_cam_implements_cnst(name)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_implements_cnst:: unrecognized microp_scheme')
   end select

end function microp_driver_implements_cnst

!===============================================================================

subroutine microp_driver_init_cnst(name, q, gcid)

   ! Initialize the microphysics constituents, if they are
   ! not read from the initial file.

   character(len=*), intent(in)  :: name     ! constituent name
   real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,          intent(in)  :: gcid(:)  ! global column id
   !-----------------------------------------------------------------------

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_init_cnst(name, q, gcid)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_init_cnst:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_init_cnst

!===============================================================================

subroutine microp_driver_init(pbuf2d)

   use iso_c_binding, only: c_int64_t

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   integer(c_int64_t) :: active_c

   ! Initialize the microphysics parameterizations
   !-----------------------------------------------------------------------
   active_c = microp_driver_init_codon(1_c_int64_t)
   if (active_c == 0_c_int64_t) return
   call microp_driver_log_direct(microp_driver_init_logged, &
        'microp_driver_init direct = codon; scheme dispatch direct = codon; micro_mg_cam_init child boundary')

   select case (microp_scheme)
   case ('MG')
      call micro_mg_cam_init(pbuf2d)
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_init:: unrecognized microp_scheme')
   end select


end subroutine microp_driver_init

!===============================================================================

subroutine microp_driver_tend(state, ptend, dtime, pbuf)

   ! Call the microphysics parameterization run methods.

   use iso_c_binding, only: c_int64_t

   ! Input arguments

   type(physics_state), intent(in)    :: state       ! State variables
   type(physics_ptend), intent(out)   :: ptend       ! Package tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in)  :: dtime                    ! Timestep
   integer(c_int64_t) :: tend_touch_c

   !======================================================================

   call microp_driver_select_impl()

   if (use_native_impl) then
      call microp_driver_tend_native(state, ptend, dtime, pbuf)
      return
   end if

   call microp_driver_select_codon_scheme()
   tend_touch_c = microp_driver_tend_codon(1_c_int64_t)
   if (tend_touch_c /= 0_c_int64_t) then
      call microp_driver_log_direct(microp_driver_tend_logged, &
           'microp_driver_tend direct = codon; scheme dispatch direct = codon; micro_mg_cam_tend child boundary')
   end if

   ! Call MG Microphysics

   select case (codon_scheme_code)
   case (1)
      call t_startf('microp_mg_tend')
      call micro_mg_cam_tend(state, ptend, dtime, pbuf)
      call t_stopf('microp_mg_tend')
   case (2)
      ! microp_driver doesn't handle this one
      continue
   case default
      call microp_driver_tend_native(state, ptend, dtime, pbuf)
   end select

end subroutine microp_driver_tend

!===============================================================================

subroutine microp_driver_select_codon_scheme()

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   integer :: i
   integer(c_int64_t), target :: scheme_ascii(len(microp_scheme))
   integer(c_int64_t), target :: scheme_code
   integer(c_int64_t), target :: status_code

   interface
      subroutine microp_driver_select_scheme_codon(scheme_len_c, scheme_ascii_p, scheme_code_p, status_p) &
           bind(c, name="microp_driver_select_scheme_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: scheme_len_c
         type(c_ptr), value :: scheme_ascii_p, scheme_code_p, status_p
      end subroutine microp_driver_select_scheme_codon
   end interface

   if (codon_scheme_selected) return

   do i = 1, len(microp_scheme)
      scheme_ascii(i) = int(iachar(microp_scheme(i:i)), c_int64_t)
   end do

   scheme_code = 0_c_int64_t
   status_code = 0_c_int64_t
   call microp_driver_select_scheme_codon( &
        int(len(microp_scheme), c_int64_t), c_loc(scheme_ascii(1)), c_loc(scheme_code), c_loc(status_code) &
   )

   if (status_code /= 0_c_int64_t) then
      codon_scheme_code = -1
   else
      codon_scheme_code = int(scheme_code)
   end if

   codon_scheme_selected = .true.

end subroutine microp_driver_select_codon_scheme

!===============================================================================

subroutine microp_driver_tend_native(state, ptend, dtime, pbuf)

   ! Call the native microphysics parameterization run methods.

   type(physics_state), intent(in)    :: state
   type(physics_ptend), intent(out)   :: ptend
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in)               :: dtime

   select case (microp_scheme)
   case ('MG')
      call t_startf('microp_mg_tend')
      call micro_mg_cam_tend(state, ptend, dtime, pbuf)
      call t_stopf('microp_mg_tend')
   case ('RK')
      ! microp_driver doesn't handle this one
      continue
   case default
      call endrun('microp_driver_tend:: unrecognized microp_scheme')
   end select

end subroutine microp_driver_tend_native

!===============================================================================

subroutine microp_driver_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('MICROP_DRIVER_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_impl = .false.
   end if

   impl_selected = .true.

   if (masterproc) then
      if (use_native_impl) then
         write(iulog,*) 'microp_driver_tend implementation = native'
      else
         write(iulog,*) 'microp_driver_tend implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine microp_driver_select_impl

!===============================================================================

subroutine microp_driver_log_direct(logged, proof_line)

   logical, intent(inout) :: logged
   character(len=*), intent(in) :: proof_line

   if (logged) return
   logged = .true.

   if (masterproc) then
      write(iulog,'(A)') trim(proof_line)
      call flush(iulog)
   end if

end subroutine microp_driver_log_direct

end module microp_driver
