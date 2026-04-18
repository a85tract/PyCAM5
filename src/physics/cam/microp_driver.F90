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

!===============================================================================
contains
!===============================================================================

subroutine microp_driver_readnl(nlfile)

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Read in namelist for microphysics scheme
   !-----------------------------------------------------------------------

   call phys_getopts(microp_scheme_out=microp_scheme)

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

   ! Return true if specified constituent is implemented by the
   ! microphysics package

   character(len=*), intent(in) :: name        ! constituent name
   logical :: microp_driver_implements_cnst    ! return value

   ! Local workspace
   integer :: m
   !-----------------------------------------------------------------------

   microp_driver_implements_cnst = .false.

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

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! Initialize the microphysics parameterizations
   !-----------------------------------------------------------------------

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

   ! Input arguments

   type(physics_state), intent(in)    :: state       ! State variables
   type(physics_ptend), intent(out)   :: ptend       ! Package tendencies
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in)  :: dtime                    ! Timestep

   !======================================================================

   call microp_driver_select_impl()

   if (use_native_impl) then
      call microp_driver_tend_native(state, ptend, dtime, pbuf)
      return
   end if

   call microp_driver_select_codon_scheme()

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
   call get_environment_variable('MICROP_DRIVER_IMPL', value=impl_name, length=n, status=status)

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

end module microp_driver
