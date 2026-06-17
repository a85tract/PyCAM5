module micro_mg_data

!
! Packing and time averaging for the MG interface.
!
! Use is as follows:
!
! 1) Figure out which columns will do averaging (mgncol) and the number of
!    levels where the microphysics will run (nlev).
!
! 2) Create an MGPacker object and assign it as follows:
!
!      packer = MGPacker(pcols, pver, mgcols, top_lev)
!
!    Where [pcols, pver] is the shape of the ultimate input/output arrays
!    that are defined at level midpoints.
!
! 3) Create a post-processing array of type MGPostProc:
!
!      post_proc = MGPostProc(packer)
!
! 4) Add pairs of pointers for packed and unpacked representations, already
!    associated with buffers of the correct dimensions:
!
!      call post_proc%add_field(unpacked_pointer, packed_pointer, &
!             fillvalue, accum_mean)
!
!    The third value is the default value used to "unpack" for points with
!    no "packed" part, and the fourth value is the method used to
!    accumulate values over time steps. These two arguments can be omitted,
!    in which case the default value will be 0 and the accumulation method
!    will take the mean.
!
! 5) Use the packed fields in MG, and for each MG iteration, do:
!
!      call post_proc%accumulate()
!
! 6) Perform final accumulation and scatter values into the unpacked arrays:
!
!      call post_proc%process_and_unpack()
!
! 7) Destroy the object when complete:
!
!      call post_proc%finalize()
!
! Caveat: MGFieldPostProc will hit a divide-by-zero error if you try to
!         take the mean over 0 steps.
!

! This include header defines CPP macros that only have an effect for debug
! builds.
#include "shr_assert.h"

use shr_kind_mod, only: r8 => shr_kind_r8
use iso_c_binding, only: c_double, c_int64_t, c_loc, c_null_ptr, c_ptr
use shr_log_mod, only: &
     errMsg => shr_log_errMsg, &
     OOBMsg => shr_log_OOBMsg
use shr_sys_mod, only: shr_sys_abort
use spmd_utils, only: masterproc, iam
use cam_logfile, only: iulog

implicit none
private

public :: MGPacker
public :: MGFieldPostProc
public :: accum_null
public :: accum_mean
public :: MGPostProc

type :: MGPacker
   ! Unpacked array dimensions.
   integer :: pcols
   integer :: pver
   ! Calculated packed dimensions, stored for convenience.
   integer :: mgncol
   integer :: nlev
   ! Which columns are packed.
   integer, allocatable :: mgcols(:)
   ! Topmost level to copy into the packed array.
   integer :: top_lev
 contains
   procedure, private :: pack_1D
   procedure, private :: pack_2D
   procedure, private :: pack_3D
   generic :: pack => pack_1D, pack_2D, pack_3D
   procedure :: pack_interface
   procedure, private :: unpack_1D
   procedure, private :: unpack_1D_array_fill
   procedure, private :: unpack_2D
   procedure, private :: unpack_2D_array_fill
   procedure, private :: unpack_3D
   procedure, private :: unpack_3D_array_fill
   generic :: unpack => unpack_1D, unpack_1D_array_fill, &
        unpack_2D, unpack_2D_array_fill, unpack_3D, unpack_3D_array_fill
   procedure :: finalize => MGPacker_finalize
end type MGPacker

interface MGPacker
   module procedure new_MGPacker
end interface

! Enum for time accumulation/averaging methods.
integer, parameter :: accum_null = 0
integer, parameter :: accum_mean = 1

integer, parameter :: pack_mode_pack_1D = 1
integer, parameter :: pack_mode_pack_2D = 2
integer, parameter :: pack_mode_pack_interface = 3
integer, parameter :: pack_mode_pack_3D = 4
integer, parameter :: pack_mode_unpack_1D = 5
integer, parameter :: pack_mode_unpack_1D_array_fill = 6
integer, parameter :: pack_mode_unpack_2D = 7
integer, parameter :: pack_mode_unpack_2D_array_fill = 8
integer, parameter :: pack_mode_unpack_3D = 9
integer, parameter :: pack_mode_unpack_3D_array_fill = 10
integer, parameter :: pack_mode_accumulate_1D = 11
integer, parameter :: pack_mode_accumulate_2D = 12
integer, parameter :: pack_mode_mean_1D = 13
integer, parameter :: pack_mode_mean_2D = 14

logical :: micro_mg_data_packer_use_native_impl = .false.
logical :: micro_mg_data_packer_impl_selected = .false.
logical :: micro_mg_data_packer_entered_logged = .false.
logical :: micro_mg_data_pack_2d_direct_logged = .false.
logical :: micro_mg_data_pack_interface_direct_logged = .false.
logical :: micro_mg_data_pack_3d_direct_logged = .false.
logical :: micro_mg_data_unpack_1d_direct_logged = .false.
logical :: micro_mg_data_unpack_2d_direct_logged = .false.
logical :: micro_mg_data_unpack_2d_array_fill_direct_logged = .false.
logical :: micro_mg_data_accumulate_1d_direct_logged = .false.
logical :: micro_mg_data_accumulate_2d_direct_logged = .false.
logical :: micro_mg_data_new_mgpacker_direct_logged = .false.
logical :: micro_mg_data_mgpacker_finalize_direct_logged = .false.
logical :: micro_mg_data_mgfield_finalize_direct_logged = .false.
logical :: micro_mg_data_new_mgpostproc_direct_logged = .false.
logical :: micro_mg_data_mgpostproc_finalize_direct_logged = .false.
logical :: micro_mg_data_add_field_1d_direct_logged = .false.
logical :: micro_mg_data_add_field_2d_direct_logged = .false.

interface
   subroutine micro_mg_data_pack_unpack_codon(mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c, &
        extent2_c, extent3_c, count1_c, num_steps_c, fillvalue_c, mgcols_p, src_p, fill_p, dst_p) &
        bind(c, name="micro_mg_data_pack_unpack_codon")
      import c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c
      integer(c_int64_t), value :: extent2_c, extent3_c, count1_c, num_steps_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: mgcols_p, src_p, fill_p, dst_p
   end subroutine micro_mg_data_pack_unpack_codon
   subroutine pack_2d_codon(mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c, &
        extent2_c, extent3_c, count1_c, num_steps_c, fillvalue_c, mgcols_p, src_p, fill_p, dst_p) &
        bind(c, name="pack_2d_codon")
      import c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c
      integer(c_int64_t), value :: extent2_c, extent3_c, count1_c, num_steps_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: mgcols_p, src_p, fill_p, dst_p
   end subroutine pack_2d_codon
   subroutine pack_interface_codon(mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c, &
        extent2_c, extent3_c, count1_c, num_steps_c, fillvalue_c, mgcols_p, src_p, fill_p, dst_p) &
        bind(c, name="pack_interface_codon")
      import c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c
      integer(c_int64_t), value :: extent2_c, extent3_c, count1_c, num_steps_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: mgcols_p, src_p, fill_p, dst_p
   end subroutine pack_interface_codon
   subroutine pack_3d_codon(mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c, &
        extent2_c, extent3_c, count1_c, num_steps_c, fillvalue_c, mgcols_p, src_p, fill_p, dst_p) &
        bind(c, name="pack_3d_codon")
      import c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c
      integer(c_int64_t), value :: extent2_c, extent3_c, count1_c, num_steps_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: mgcols_p, src_p, fill_p, dst_p
   end subroutine pack_3d_codon
   subroutine unpack_1d_codon(mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c, &
        extent2_c, extent3_c, count1_c, num_steps_c, fillvalue_c, mgcols_p, src_p, fill_p, dst_p) &
        bind(c, name="unpack_1d_codon")
      import c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c
      integer(c_int64_t), value :: extent2_c, extent3_c, count1_c, num_steps_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: mgcols_p, src_p, fill_p, dst_p
   end subroutine unpack_1d_codon
   subroutine unpack_2d_codon(mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c, &
        extent2_c, extent3_c, count1_c, num_steps_c, fillvalue_c, mgcols_p, src_p, fill_p, dst_p) &
        bind(c, name="unpack_2d_codon")
      import c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c
      integer(c_int64_t), value :: extent2_c, extent3_c, count1_c, num_steps_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: mgcols_p, src_p, fill_p, dst_p
   end subroutine unpack_2d_codon
   subroutine unpack_2d_array_fill_codon(mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c, &
        extent2_c, extent3_c, count1_c, num_steps_c, fillvalue_c, mgcols_p, src_p, fill_p, dst_p) &
        bind(c, name="unpack_2d_array_fill_codon")
      import c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: mode_c, pcols_c, pver_c, mgncol_c, nlev_c, top_lev_c
      integer(c_int64_t), value :: extent2_c, extent3_c, count1_c, num_steps_c
      real(c_double), value :: fillvalue_c
      type(c_ptr), value :: mgcols_p, src_p, fill_p, dst_p
   end subroutine unpack_2d_array_fill_codon
   subroutine new_mgpacker_codon(pcols_c, pver_c, mgcols_count_c, top_lev_c, plan_p) &
        bind(c, name="new_mgpacker_codon")
      import c_int64_t, c_ptr
      integer(c_int64_t), value :: pcols_c, pver_c, mgcols_count_c, top_lev_c
      type(c_ptr), value :: plan_p
   end subroutine new_mgpacker_codon
   subroutine mgpacker_finalize_codon(plan_p) bind(c, name="mgpacker_finalize_codon")
      import c_ptr
      type(c_ptr), value :: plan_p
   end subroutine mgpacker_finalize_codon
   function mgfieldpostproc_1d_codon(flag_c) result(out_c) bind(c, name="mgfieldpostproc_1d_codon")
      import c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function mgfieldpostproc_1d_codon
   function mgfieldpostproc_2d_codon(flag_c) result(out_c) bind(c, name="mgfieldpostproc_2d_codon")
      import c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function mgfieldpostproc_2d_codon
   subroutine mgfieldpostproc_finalize_codon(plan_p) bind(c, name="mgfieldpostproc_finalize_codon")
      import c_ptr
      type(c_ptr), value :: plan_p
   end subroutine mgfieldpostproc_finalize_codon
   function mgfieldpostproc_process_and_unpack_codon(flag_c) result(out_c) &
        bind(c, name="mgfieldpostproc_process_and_unpack_codon")
      import c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function mgfieldpostproc_process_and_unpack_codon
   function mgfieldpostproc_unpack_only_codon(flag_c) result(out_c) bind(c, name="mgfieldpostproc_unpack_only_codon")
      import c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function mgfieldpostproc_unpack_only_codon
   subroutine new_mgpostproc_codon(plan_p) bind(c, name="new_mgpostproc_codon")
      import c_ptr
      type(c_ptr), value :: plan_p
   end subroutine new_mgpostproc_codon
   subroutine mgpostproc_finalize_codon(field_count_c, plan_p) bind(c, name="mgpostproc_finalize_codon")
      import c_int64_t, c_ptr
      integer(c_int64_t), value :: field_count_c
      type(c_ptr), value :: plan_p
   end subroutine mgpostproc_finalize_codon
   subroutine add_field_1d_codon(fill_present_c, accum_present_c, fill_in_c, accum_in_c, &
        plan_p, fill_out_p) bind(c, name="add_field_1d_codon")
      import c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: fill_present_c, accum_present_c, accum_in_c
      real(c_double), value :: fill_in_c
      type(c_ptr), value :: plan_p, fill_out_p
   end subroutine add_field_1d_codon
   subroutine add_field_2d_codon(fill_present_c, accum_present_c, fill_in_c, accum_in_c, &
        plan_p, fill_out_p) bind(c, name="add_field_2d_codon")
      import c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: fill_present_c, accum_present_c, accum_in_c
      real(c_double), value :: fill_in_c
      type(c_ptr), value :: plan_p, fill_out_p
   end subroutine add_field_2d_codon
   function mgpostproc_accumulate_codon(flag_c) result(out_c) bind(c, name="mgpostproc_accumulate_codon")
      import c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function mgpostproc_accumulate_codon
   function mgpostproc_process_and_unpack_codon(flag_c) result(out_c) bind(c, name="mgpostproc_process_and_unpack_codon")
      import c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function mgpostproc_process_and_unpack_codon
   function mgpostproc_unpack_only_codon(flag_c) result(out_c) bind(c, name="mgpostproc_unpack_only_codon")
      import c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function mgpostproc_unpack_only_codon
   function mgpostproc_copy_codon(flag_c) result(out_c) bind(c, name="mgpostproc_copy_codon")
      import c_int64_t
      integer(c_int64_t), value :: flag_c
      integer(c_int64_t) :: out_c
   end function mgpostproc_copy_codon
end interface

type :: MGFieldPostProc
   integer :: accum_method = -1
   integer :: rank = -1
   integer :: num_steps = 0
   real(r8) :: fillvalue = 0._r8
   real(r8), pointer :: unpacked_1D(:) => null()
   real(r8), pointer :: packed_1D(:) => null()
   real(r8), allocatable :: buffer_1D(:)
   real(r8), pointer :: unpacked_2D(:,:) => null()
   real(r8), pointer :: packed_2D(:,:) => null()
   real(r8), allocatable :: buffer_2D(:,:)
 contains
   procedure :: accumulate => MGFieldPostProc_accumulate
   procedure :: process_and_unpack => MGFieldPostProc_process_and_unpack
   procedure :: unpack_only => MGFieldPostProc_unpack_only
   procedure :: finalize => MGFieldPostProc_finalize
end type MGFieldPostProc

interface MGFieldPostProc
   module procedure MGFieldPostProc_1D
   module procedure MGFieldPostProc_2D
end interface MGFieldPostProc

#define VECTOR_NAME MGFieldPostProcVec
#define TYPE_NAME type(MGFieldPostProc)
#define THROW(string) call shr_sys_abort(string)

public :: VECTOR_NAME

#include "dynamic_vector_typedef.inc"

type MGPostProc
   type(MGPacker) :: packer
   type(MGFieldPostProcVec) :: field_procs
 contains
   procedure, private :: add_field_1D
   procedure, private :: add_field_2D
   generic :: add_field => add_field_1D, add_field_2D
   procedure :: accumulate => MGPostProc_accumulate
   procedure :: process_and_unpack => MGPostProc_process_and_unpack
   procedure :: unpack_only => MGPostProc_unpack_only
   procedure :: finalize => MGPostProc_finalize
   procedure, private :: MGPostProc_copy
   generic :: assignment(=) => MGPostProc_copy
end type MGPostProc

interface MGPostProc
   module procedure new_MGPostProc
end interface MGPostProc

contains

function new_MGPacker(pcols, pver, mgcols, top_lev)
  integer, intent(in) :: pcols, pver
  integer, intent(in) :: mgcols(:)
  integer, intent(in) :: top_lev

  type(MGPacker) :: new_MGPacker
  integer(c_int64_t), target :: plan(4)

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     new_MGPacker%pcols = pcols
     new_MGPacker%pver = pver
     new_MGPacker%mgncol = size(mgcols)
     new_MGPacker%nlev = pver - top_lev + 1
     allocate(new_MGPacker%mgcols(new_MGPacker%mgncol))
     new_MGPacker%mgcols = mgcols
     new_MGPacker%top_lev = top_lev
     return
  end if

  call new_mgpacker_codon(int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(size(mgcols), c_int64_t), int(top_lev, c_int64_t), c_loc(plan))
  call micro_mg_data_log_new_mgpacker_direct()

  new_MGPacker%pcols = int(plan(1))
  new_MGPacker%pver = int(plan(2))
  new_MGPacker%mgncol = int(plan(3))
  new_MGPacker%nlev = int(plan(4))

  allocate(new_MGPacker%mgcols(new_MGPacker%mgncol))
  new_MGPacker%mgcols = mgcols
  new_MGPacker%top_lev = top_lev

end function new_MGPacker

! Rely on the fact that intent(out) forces the compiler to deallocate all
! allocatable components and restart the type from scratch. Although
! compiler support for finalization varies, this seems to be one of the few
! cases where all major compilers are reliable, and humans are not.
subroutine MGPacker_finalize(self)
  class(MGPacker), intent(out) :: self
  integer(c_int64_t), target :: plan(1)

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) return

  call mgpacker_finalize_codon(c_loc(plan))
  if (plan(1) /= 0_c_int64_t) call micro_mg_data_log_mgpacker_finalize_direct()
end subroutine MGPacker_finalize

subroutine micro_mg_data_select_packer_impl()
  character(len=32) :: impl_name
  integer :: n, status

  if (micro_mg_data_packer_impl_selected) return

  call cam_codon_get_impl('MICRO_MG_DATA_PACKER_IMPL', impl_name, n, status)
  if (status == 0 .and. n > 0) then
     micro_mg_data_packer_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     micro_mg_data_packer_use_native_impl = .false.
  end if

  if (masterproc) then
     if (micro_mg_data_packer_use_native_impl) then
        write(iulog,*) 'micro_mg_data_packer implementation = native'
     else
        write(iulog,*) 'micro_mg_data_packer implementation = codon'
     end if
  end if

  micro_mg_data_packer_impl_selected = .true.
end subroutine micro_mg_data_select_packer_impl

subroutine micro_mg_data_log_packer_entry()
  if (masterproc .and. .not. micro_mg_data_packer_entered_logged) then
     write(iulog,*) 'micro_mg_data_packer entered (pack/unpack helpers = codon)'
     micro_mg_data_packer_entered_logged = .true.
  end if
end subroutine micro_mg_data_log_packer_entry

subroutine micro_mg_data_append_direct_proof(proof_line)
  character(len=*), intent(in) :: proof_line
  character(len=512) :: proof_file
  character(len=640) :: rank_proof_file
  integer :: n, status, unitno

  proof_file = ''
  call get_environment_variable('MICRO_MG_DATA_PACKER_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     write(rank_proof_file,'(A,".",I0,".proof")') trim(proof_file(:n)), iam
     open(newunit=unitno, file=trim(rank_proof_file), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if
end subroutine micro_mg_data_append_direct_proof

subroutine micro_mg_data_log_pack_2d_direct()
  if (.not. micro_mg_data_pack_2d_direct_logged) then
     call micro_mg_data_append_direct_proof('pack_2D direct = codon')
     if (masterproc) write(iulog,*) 'pack_2D direct = codon'
     micro_mg_data_pack_2d_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_pack_2d_direct

subroutine micro_mg_data_log_pack_interface_direct()
  if (.not. micro_mg_data_pack_interface_direct_logged) then
     call micro_mg_data_append_direct_proof('pack_interface direct = codon')
     if (masterproc) write(iulog,*) 'pack_interface direct = codon'
     micro_mg_data_pack_interface_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_pack_interface_direct

subroutine micro_mg_data_log_pack_3d_direct()
  if (.not. micro_mg_data_pack_3d_direct_logged) then
     call micro_mg_data_append_direct_proof('pack_3D direct = codon')
     if (masterproc) write(iulog,*) 'pack_3D direct = codon'
     micro_mg_data_pack_3d_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_pack_3d_direct

subroutine micro_mg_data_log_unpack_1d_direct()
  if (.not. micro_mg_data_unpack_1d_direct_logged) then
     call micro_mg_data_append_direct_proof('unpack_1D direct = codon')
     if (masterproc) write(iulog,*) 'unpack_1D direct = codon'
     micro_mg_data_unpack_1d_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_unpack_1d_direct

subroutine micro_mg_data_log_unpack_2d_direct()
  if (.not. micro_mg_data_unpack_2d_direct_logged) then
     call micro_mg_data_append_direct_proof('unpack_2D direct = codon')
     if (masterproc) write(iulog,*) 'unpack_2D direct = codon'
     micro_mg_data_unpack_2d_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_unpack_2d_direct

subroutine micro_mg_data_log_unpack_2d_array_fill_direct()
  if (.not. micro_mg_data_unpack_2d_array_fill_direct_logged) then
     call micro_mg_data_append_direct_proof('unpack_2D_array_fill direct = codon')
     if (masterproc) write(iulog,*) 'unpack_2D_array_fill direct = codon'
     micro_mg_data_unpack_2d_array_fill_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_unpack_2d_array_fill_direct

subroutine micro_mg_data_log_accumulate_1d_direct()
  if (.not. micro_mg_data_accumulate_1d_direct_logged) then
     call micro_mg_data_append_direct_proof('MGFieldPostProc_accumulate rank1 direct = codon; allocation native island')
     if (masterproc) write(iulog,*) 'MGFieldPostProc_accumulate rank1 direct = codon; allocation native island'
     micro_mg_data_accumulate_1d_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_accumulate_1d_direct

subroutine micro_mg_data_log_accumulate_2d_direct()
  if (.not. micro_mg_data_accumulate_2d_direct_logged) then
     call micro_mg_data_append_direct_proof('MGFieldPostProc_accumulate rank2 direct = codon; allocation native island')
     if (masterproc) write(iulog,*) 'MGFieldPostProc_accumulate rank2 direct = codon; allocation native island'
     micro_mg_data_accumulate_2d_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_accumulate_2d_direct

subroutine micro_mg_data_log_new_mgpacker_direct()
  if (.not. micro_mg_data_new_mgpacker_direct_logged) then
     call micro_mg_data_append_direct_proof( &
          'new_MGPacker direct = codon; dimension plan direct = codon; mgcols allocation native boundary')
     if (masterproc) write(iulog,*) &
          'new_MGPacker direct = codon; dimension plan direct = codon; mgcols allocation native boundary'
     micro_mg_data_new_mgpacker_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_new_mgpacker_direct

subroutine micro_mg_data_log_mgpacker_finalize_direct()
  if (.not. micro_mg_data_mgpacker_finalize_direct_logged) then
     call micro_mg_data_append_direct_proof( &
          'MGPacker_finalize direct = codon; finalization action plan direct = codon; allocatable cleanup native boundary')
     if (masterproc) write(iulog,*) &
          'MGPacker_finalize direct = codon; finalization action plan direct = codon; allocatable cleanup native boundary'
     micro_mg_data_mgpacker_finalize_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_mgpacker_finalize_direct

subroutine micro_mg_data_log_mgfield_finalize_direct()
  if (.not. micro_mg_data_mgfield_finalize_direct_logged) then
     call micro_mg_data_append_direct_proof( &
          'MGFieldPostProc_finalize direct = codon; finalization action plan direct = codon; pointer cleanup native boundary')
     if (masterproc) write(iulog,*) &
          'MGFieldPostProc_finalize direct = codon; finalization action plan direct = codon; pointer cleanup native boundary'
     micro_mg_data_mgfield_finalize_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_mgfield_finalize_direct

subroutine micro_mg_data_log_new_mgpostproc_direct()
  if (.not. micro_mg_data_new_mgpostproc_direct_logged) then
     call micro_mg_data_append_direct_proof( &
          'new_MGPostProc direct = codon; constructor action plan direct = codon; vector clear native boundary')
     if (masterproc) write(iulog,*) &
          'new_MGPostProc direct = codon; constructor action plan direct = codon; vector clear native boundary'
     micro_mg_data_new_mgpostproc_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_new_mgpostproc_direct

subroutine micro_mg_data_log_mgpostproc_finalize_direct()
  if (.not. micro_mg_data_mgpostproc_finalize_direct_logged) then
     call micro_mg_data_append_direct_proof( &
          'MGPostProc_finalize direct = codon; recursive finalization plan direct = codon; vector cleanup native boundary')
     if (masterproc) write(iulog,*) &
          'MGPostProc_finalize direct = codon; recursive finalization plan direct = codon; vector cleanup native boundary'
     micro_mg_data_mgpostproc_finalize_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_mgpostproc_finalize_direct

subroutine micro_mg_data_log_add_field_1d_direct()
  if (.not. micro_mg_data_add_field_1d_direct_logged) then
     call micro_mg_data_append_direct_proof( &
          'add_field_1D direct = codon; field registration defaults direct = codon; vector push native boundary')
     if (masterproc) write(iulog,*) &
          'add_field_1D direct = codon; field registration defaults direct = codon; vector push native boundary'
     micro_mg_data_add_field_1d_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_add_field_1d_direct

subroutine micro_mg_data_log_add_field_2d_direct()
  if (.not. micro_mg_data_add_field_2d_direct_logged) then
     call micro_mg_data_append_direct_proof( &
          'add_field_2D direct = codon; field registration defaults direct = codon; vector push native boundary')
     if (masterproc) write(iulog,*) &
          'add_field_2D direct = codon; field registration defaults direct = codon; vector push native boundary'
     micro_mg_data_add_field_2d_direct_logged = .true.
  end if
end subroutine micro_mg_data_log_add_field_2d_direct

subroutine micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  class(MGPacker), intent(in) :: self
  integer(c_int64_t), intent(out) :: mgcols_c(:)

  if (self%mgncol > 0) then
     mgcols_c(1:self%mgncol) = int(self%mgcols(1:self%mgncol), c_int64_t)
  end if
end subroutine micro_mg_data_copy_mgcols_i64

subroutine micro_mg_data_pack_1D_codon_wrap(self, unpacked, packed)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: unpacked(:)
  real(r8), target, intent(out) :: packed(:)
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  src_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(unpacked) > 0) src_p = c_loc(unpacked)
  if (size(packed) > 0) dst_p = c_loc(packed)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_pack_1D, c_int64_t), &
       int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
       int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), 1_c_int64_t, 1_c_int64_t, &
       int(self%mgncol, c_int64_t), 1_c_int64_t, 0.0_c_double, c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
end subroutine micro_mg_data_pack_1D_codon_wrap

subroutine micro_mg_data_pack_2D_codon_wrap(self, unpacked, packed, extent2, mode)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: unpacked(:,:)
  real(r8), target, intent(out) :: packed(:,:)
  integer, intent(in) :: extent2, mode
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  src_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(unpacked) > 0) src_p = c_loc(unpacked)
  if (size(packed) > 0) dst_p = c_loc(packed)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_pack_unpack_codon(int(mode, c_int64_t), &
       int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
       int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(extent2, c_int64_t), &
       1_c_int64_t, int(self%mgncol, c_int64_t), 1_c_int64_t, 0.0_c_double, &
       c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
end subroutine micro_mg_data_pack_2D_codon_wrap

subroutine micro_mg_data_pack_3D_codon_wrap(self, unpacked, packed, extent3)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: unpacked(:,:,:)
  real(r8), target, intent(out) :: packed(:,:,:)
  integer, intent(in) :: extent3
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  src_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(unpacked) > 0) src_p = c_loc(unpacked)
  if (size(packed) > 0) dst_p = c_loc(packed)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_pack_3D, c_int64_t), &
       int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
       int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(self%nlev, c_int64_t), &
       int(extent3, c_int64_t), int(self%mgncol, c_int64_t), 1_c_int64_t, 0.0_c_double, &
       c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
end subroutine micro_mg_data_pack_3D_codon_wrap

subroutine micro_mg_data_unpack_1D_codon_wrap(self, packed, unpacked, fill)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:)
  real(r8), target, intent(out) :: unpacked(:)
  real(r8), intent(in) :: fill
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  src_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(packed) > 0) src_p = c_loc(packed)
  if (size(unpacked) > 0) dst_p = c_loc(unpacked)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_unpack_1D, c_int64_t), &
       int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
       int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), 1_c_int64_t, 1_c_int64_t, &
       int(self%pcols, c_int64_t), 1_c_int64_t, real(fill, c_double), c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
end subroutine micro_mg_data_unpack_1D_codon_wrap

subroutine micro_mg_data_unpack_1D_array_fill_codon_wrap(self, packed, fill, unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:)
  real(r8), target, intent(in) :: fill(:)
  real(r8), target, intent(out) :: unpacked(:)
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, fill_p, dst_p

  call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  src_p = c_null_ptr
  fill_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(packed) > 0) src_p = c_loc(packed)
  if (size(fill) > 0) fill_p = c_loc(fill)
  if (size(unpacked) > 0) dst_p = c_loc(unpacked)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_unpack_1D_array_fill, c_int64_t), &
       int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
       int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), 1_c_int64_t, 1_c_int64_t, &
       int(self%pcols, c_int64_t), 1_c_int64_t, 0.0_c_double, c_loc(mgcols_c), src_p, fill_p, dst_p)
end subroutine micro_mg_data_unpack_1D_array_fill_codon_wrap

subroutine micro_mg_data_unpack_2D_codon_wrap(self, packed, unpacked, fill, extent2)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:,:)
  real(r8), target, intent(out) :: unpacked(:,:)
  real(r8), intent(in) :: fill
  integer, intent(in) :: extent2
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  src_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(packed) > 0) src_p = c_loc(packed)
  if (size(unpacked) > 0) dst_p = c_loc(unpacked)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_unpack_2D, c_int64_t), &
       int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
       int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(extent2, c_int64_t), &
       1_c_int64_t, int(self%pcols, c_int64_t), 1_c_int64_t, real(fill, c_double), &
       c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
end subroutine micro_mg_data_unpack_2D_codon_wrap

subroutine micro_mg_data_unpack_2D_array_fill_codon_wrap(self, packed, fill, unpacked, extent2)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:,:)
  real(r8), target, intent(in) :: fill(:,:)
  real(r8), target, intent(out) :: unpacked(:,:)
  integer, intent(in) :: extent2
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, fill_p, dst_p

  call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  src_p = c_null_ptr
  fill_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(packed) > 0) src_p = c_loc(packed)
  if (size(fill) > 0) fill_p = c_loc(fill)
  if (size(unpacked) > 0) dst_p = c_loc(unpacked)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_unpack_2D_array_fill, c_int64_t), &
       int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
       int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(extent2, c_int64_t), &
       1_c_int64_t, int(self%pcols, c_int64_t), 1_c_int64_t, 0.0_c_double, &
       c_loc(mgcols_c), src_p, fill_p, dst_p)
end subroutine micro_mg_data_unpack_2D_array_fill_codon_wrap

subroutine micro_mg_data_unpack_3D_codon_wrap(self, packed, unpacked, fill, extent3)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:,:,:)
  real(r8), target, intent(out) :: unpacked(:,:,:)
  real(r8), intent(in) :: fill
  integer, intent(in) :: extent3
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  src_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(packed) > 0) src_p = c_loc(packed)
  if (size(unpacked) > 0) dst_p = c_loc(unpacked)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_unpack_3D, c_int64_t), &
       int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
       int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(self%nlev, c_int64_t), &
       int(extent3, c_int64_t), int(self%pcols, c_int64_t), 1_c_int64_t, real(fill, c_double), &
       c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
end subroutine micro_mg_data_unpack_3D_codon_wrap

subroutine micro_mg_data_unpack_3D_array_fill_codon_wrap(self, packed, fill, unpacked, extent3)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:,:,:)
  real(r8), target, intent(in) :: fill(:,:,:)
  real(r8), target, intent(out) :: unpacked(:,:,:)
  integer, intent(in) :: extent3
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, fill_p, dst_p

  call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
  src_p = c_null_ptr
  fill_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(packed) > 0) src_p = c_loc(packed)
  if (size(fill) > 0) fill_p = c_loc(fill)
  if (size(unpacked) > 0) dst_p = c_loc(unpacked)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_unpack_3D_array_fill, c_int64_t), &
       int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
       int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(self%nlev, c_int64_t), &
       int(extent3, c_int64_t), int(self%pcols, c_int64_t), 1_c_int64_t, 0.0_c_double, &
       c_loc(mgcols_c), src_p, fill_p, dst_p)
end subroutine micro_mg_data_unpack_3D_array_fill_codon_wrap

subroutine micro_mg_data_accumulate_1D_codon_wrap(packed, buffer)
  real(r8), target, intent(in) :: packed(:)
  real(r8), target, intent(inout) :: buffer(:)
  type(c_ptr) :: src_p, dst_p

  src_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(packed) > 0) src_p = c_loc(packed)
  if (size(buffer) > 0) dst_p = c_loc(buffer)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_log_accumulate_1d_direct()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_accumulate_1D, c_int64_t), &
       0_c_int64_t, 0_c_int64_t, 0_c_int64_t, 0_c_int64_t, 0_c_int64_t, &
       1_c_int64_t, 1_c_int64_t, int(size(packed), c_int64_t), 1_c_int64_t, &
       0.0_c_double, c_null_ptr, src_p, c_null_ptr, dst_p)
end subroutine micro_mg_data_accumulate_1D_codon_wrap

subroutine micro_mg_data_accumulate_2D_codon_wrap(packed, buffer)
  real(r8), target, intent(in) :: packed(:,:)
  real(r8), target, intent(inout) :: buffer(:,:)
  type(c_ptr) :: src_p, dst_p

  src_p = c_null_ptr
  dst_p = c_null_ptr
  if (size(packed) > 0) src_p = c_loc(packed)
  if (size(buffer) > 0) dst_p = c_loc(buffer)

  call micro_mg_data_log_packer_entry()
  call micro_mg_data_log_accumulate_2d_direct()
  call micro_mg_data_pack_unpack_codon(int(pack_mode_accumulate_2D, c_int64_t), &
       0_c_int64_t, 0_c_int64_t, 0_c_int64_t, 0_c_int64_t, 0_c_int64_t, &
       int(size(packed, 2), c_int64_t), 1_c_int64_t, int(size(packed, 1), c_int64_t), &
       1_c_int64_t, 0.0_c_double, c_null_ptr, src_p, c_null_ptr, dst_p)
end subroutine micro_mg_data_accumulate_2D_codon_wrap

function pack_1D(self, unpacked) result(packed)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: unpacked(:)

  real(r8), target :: packed(self%mgncol)

  SHR_ASSERT(size(unpacked) == self%pcols, errMsg(__FILE__, __LINE__))

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     packed = unpacked(self%mgcols)
  else
     call micro_mg_data_pack_1D_codon_wrap(self, unpacked, packed)
  end if

end function pack_1D

! Separation of pack and pack_interface is to workaround a PGI bug.
function pack_2D(self, unpacked) result(packed)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: unpacked(:,:)

  real(r8), target :: packed(self%mgncol,self%nlev)
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  SHR_ASSERT(size(unpacked, 1) == self%pcols, errMsg(__FILE__, __LINE__))

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     packed = unpacked(self%mgcols,self%top_lev:)
  else
     call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
     src_p = c_null_ptr
     dst_p = c_null_ptr
     if (size(unpacked) > 0) src_p = c_loc(unpacked)
     if (size(packed) > 0) dst_p = c_loc(packed)
     call micro_mg_data_log_packer_entry()
     call micro_mg_data_log_pack_2d_direct()
     call pack_2d_codon(int(pack_mode_pack_2D, c_int64_t), &
          int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
          int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(self%nlev, c_int64_t), &
          1_c_int64_t, int(self%mgncol, c_int64_t), 1_c_int64_t, 0.0_c_double, &
          c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
  end if

end function pack_2D

function pack_interface(self, unpacked) result(packed)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: unpacked(:,:)

  real(r8), target :: packed(self%mgncol,self%nlev+1)
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     packed = unpacked(self%mgcols,self%top_lev:)
  else
     call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
     src_p = c_null_ptr
     dst_p = c_null_ptr
     if (size(unpacked) > 0) src_p = c_loc(unpacked)
     if (size(packed) > 0) dst_p = c_loc(packed)
     call micro_mg_data_log_packer_entry()
     call micro_mg_data_log_pack_interface_direct()
     call pack_interface_codon(int(pack_mode_pack_interface, c_int64_t), &
          int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
          int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(self%nlev+1, c_int64_t), &
          1_c_int64_t, int(self%mgncol, c_int64_t), 1_c_int64_t, 0.0_c_double, &
          c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
  end if

end function pack_interface

function pack_3D(self, unpacked) result(packed)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: unpacked(:,:,:)

  real(r8), target :: packed(self%mgncol,self%nlev,size(unpacked, 3))
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  SHR_ASSERT(size(unpacked,1) == self%pcols, errMsg(__FILE__, __LINE__))

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     packed = unpacked(self%mgcols,self%top_lev:,:)
  else
     call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
     src_p = c_null_ptr
     dst_p = c_null_ptr
     if (size(unpacked) > 0) src_p = c_loc(unpacked)
     if (size(packed) > 0) dst_p = c_loc(packed)
     call micro_mg_data_log_packer_entry()
     call micro_mg_data_log_pack_3d_direct()
     call pack_3d_codon(int(pack_mode_pack_3D, c_int64_t), &
          int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
          int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(self%nlev, c_int64_t), &
          int(size(unpacked, 3), c_int64_t), int(self%mgncol, c_int64_t), 1_c_int64_t, 0.0_c_double, &
          c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
  end if

end function pack_3D

function unpack_1D(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:)
  real(r8), intent(in) :: fill

  real(r8), target :: unpacked(self%pcols)
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  SHR_ASSERT(size(packed) == self%mgncol, errMsg(__FILE__, __LINE__))

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     unpacked = fill
     unpacked(self%mgcols) = packed
  else
     call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
     src_p = c_null_ptr
     dst_p = c_null_ptr
     if (size(packed) > 0) src_p = c_loc(packed)
     if (size(unpacked) > 0) dst_p = c_loc(unpacked)
     call micro_mg_data_log_packer_entry()
     call micro_mg_data_log_unpack_1d_direct()
     call unpack_1d_codon(int(pack_mode_unpack_1D, c_int64_t), &
          int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
          int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), 1_c_int64_t, 1_c_int64_t, &
          int(self%pcols, c_int64_t), 1_c_int64_t, real(fill, c_double), c_loc(mgcols_c), &
          src_p, c_null_ptr, dst_p)
  end if

end function unpack_1D

function unpack_1D_array_fill(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:)
  real(r8), target, intent(in) :: fill(:)

  real(r8), target :: unpacked(self%pcols)

  SHR_ASSERT(size(packed) == self%mgncol, errMsg(__FILE__, __LINE__))

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     unpacked = fill
     unpacked(self%mgcols) = packed
  else
     call micro_mg_data_unpack_1D_array_fill_codon_wrap(self, packed, fill, unpacked)
  end if

end function unpack_1D_array_fill

function unpack_2D(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:,:)
  real(r8), intent(in) :: fill

  real(r8), target :: unpacked(self%pcols,self%pver+size(packed, 2)-self%nlev)
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, dst_p

  SHR_ASSERT(size(packed, 1) == self%mgncol, errMsg(__FILE__, __LINE__))

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     unpacked = fill
     unpacked(self%mgcols,self%top_lev:) = packed
  else
     call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
     src_p = c_null_ptr
     dst_p = c_null_ptr
     if (size(packed) > 0) src_p = c_loc(packed)
     if (size(unpacked) > 0) dst_p = c_loc(unpacked)
     call micro_mg_data_log_packer_entry()
     call micro_mg_data_log_unpack_2d_direct()
     call unpack_2d_codon(int(pack_mode_unpack_2D, c_int64_t), &
          int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
          int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(size(packed, 2), c_int64_t), &
          1_c_int64_t, int(self%pcols, c_int64_t), 1_c_int64_t, real(fill, c_double), &
          c_loc(mgcols_c), src_p, c_null_ptr, dst_p)
  end if

end function unpack_2D

function unpack_2D_array_fill(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:,:)
  real(r8), target, intent(in) :: fill(:,:)

  real(r8), target :: unpacked(self%pcols,self%pver+size(packed, 2)-self%nlev)
  integer(c_int64_t), target :: mgcols_c(max(1,self%mgncol))
  type(c_ptr) :: src_p, fill_p, dst_p

  SHR_ASSERT(size(packed, 1) == self%mgncol, errMsg(__FILE__, __LINE__))

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     unpacked = fill
     unpacked(self%mgcols,self%top_lev:) = packed
  else
     call micro_mg_data_copy_mgcols_i64(self, mgcols_c)
     src_p = c_null_ptr
     fill_p = c_null_ptr
     dst_p = c_null_ptr
     if (size(packed) > 0) src_p = c_loc(packed)
     if (size(fill) > 0) fill_p = c_loc(fill)
     if (size(unpacked) > 0) dst_p = c_loc(unpacked)
     call micro_mg_data_log_packer_entry()
     call micro_mg_data_log_unpack_2d_array_fill_direct()
     call unpack_2d_array_fill_codon(int(pack_mode_unpack_2D_array_fill, c_int64_t), &
          int(self%pcols, c_int64_t), int(self%pver, c_int64_t), int(self%mgncol, c_int64_t), &
          int(self%nlev, c_int64_t), int(self%top_lev, c_int64_t), int(size(packed, 2), c_int64_t), &
          1_c_int64_t, int(self%pcols, c_int64_t), 1_c_int64_t, 0.0_c_double, &
          c_loc(mgcols_c), src_p, fill_p, dst_p)
  end if

end function unpack_2D_array_fill

function unpack_3D(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:,:,:)
  real(r8), intent(in) :: fill

  real(r8), target :: unpacked(self%pcols,self%pver,size(packed, 3))

  SHR_ASSERT(size(packed, 1) == self%mgncol, errMsg(__FILE__, __LINE__))

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     unpacked = fill
     unpacked(self%mgcols,self%top_lev:,:) = packed
  else
     call micro_mg_data_unpack_3D_codon_wrap(self, packed, unpacked, fill, size(packed, 3))
  end if

end function unpack_3D

function unpack_3D_array_fill(self, packed, fill) result(unpacked)
  class(MGPacker), intent(in) :: self
  real(r8), target, intent(in) :: packed(:,:,:)
  real(r8), target, intent(in) :: fill(:,:,:)

  real(r8), target :: unpacked(self%pcols,self%pver,size(packed, 3))

  SHR_ASSERT(size(packed, 1) == self%mgncol, errMsg(__FILE__, __LINE__))

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     unpacked = fill
     unpacked(self%mgcols,self%top_lev:,:) = packed
  else
     call micro_mg_data_unpack_3D_array_fill_codon_wrap(self, packed, fill, unpacked, size(packed, 3))
  end if

end function unpack_3D_array_fill

function MGFieldPostProc_1D(unpacked_ptr, packed_ptr, fillvalue, &
     accum_method) result(field_proc)
  real(r8), pointer, intent(in) :: unpacked_ptr(:)
  real(r8), pointer, intent(in) :: packed_ptr(:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  type(MGFieldPostProc) :: field_proc

  call micro_mg_data_select_packer_impl()
  if (.not. micro_mg_data_packer_use_native_impl) then
     if (mgfieldpostproc_1d_codon(1_c_int64_t) == 0_c_int64_t) return
  end if

  field_proc%rank = 1
  field_proc%unpacked_1D => unpacked_ptr
  field_proc%packed_1D => packed_ptr
  if (present(fillvalue)) then
     field_proc%fillvalue = fillvalue
  else
     field_proc%fillvalue = 0._r8
  end if
  if (present(accum_method)) then
     field_proc%accum_method = accum_method
  else
     field_proc%accum_method = accum_mean
  end if

end function MGFieldPostProc_1D

function MGFieldPostProc_2D(unpacked_ptr, packed_ptr, fillvalue, &
     accum_method) result(field_proc)
  real(r8), pointer, intent(in) :: unpacked_ptr(:,:)
  real(r8), pointer, intent(in) :: packed_ptr(:,:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  type(MGFieldPostProc) :: field_proc

  call micro_mg_data_select_packer_impl()
  if (.not. micro_mg_data_packer_use_native_impl) then
     if (mgfieldpostproc_2d_codon(1_c_int64_t) == 0_c_int64_t) return
  end if

  field_proc%rank = 2
  field_proc%unpacked_2D => unpacked_ptr
  field_proc%packed_2D => packed_ptr
  if (present(fillvalue)) then
     field_proc%fillvalue = fillvalue
  else
     field_proc%fillvalue = 0._r8
  end if
  if (present(accum_method)) then
     field_proc%accum_method = accum_method
  else
     field_proc%accum_method = accum_mean
  end if

end function MGFieldPostProc_2D

! Use the same intent(out) trick as for MGPacker, which is actually more
! useful here.
subroutine MGFieldPostProc_finalize(self)
  class(MGFieldPostProc), intent(out) :: self
  integer(c_int64_t), target :: plan(1)

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) return

  call mgfieldpostproc_finalize_codon(c_loc(plan))
  if (plan(1) /= 0_c_int64_t) call micro_mg_data_log_mgfield_finalize_direct()
end subroutine MGFieldPostProc_finalize

subroutine MGFieldPostProc_accumulate(self)
  class(MGFieldPostProc), intent(inout) :: self

  call micro_mg_data_select_packer_impl()

  select case (self%accum_method)
  case (accum_null)
     ! "Null" method does nothing.
  case (accum_mean)
     ! Allocation is done on the first accumulation step to allow the
     ! MGFieldPostProc to be copied after construction without copying the
     ! allocated array (until this function is first called).
     self%num_steps = self%num_steps + 1
     select case (self%rank)
     case (1)
        SHR_ASSERT(associated(self%packed_1D), errMsg(__FILE__, __LINE__))
        if (.not. allocated(self%buffer_1D)) then
           allocate(self%buffer_1D(size(self%packed_1D)))
           self%buffer_1D = 0._r8
        end if
        if (micro_mg_data_packer_use_native_impl) then
           self%buffer_1D = self%buffer_1D + self%packed_1D
        else
           call micro_mg_data_accumulate_1D_codon_wrap(self%packed_1D, self%buffer_1D)
        end if
     case (2)
        SHR_ASSERT(associated(self%packed_2D), errMsg(__FILE__, __LINE__))
        if (.not. allocated(self%buffer_2D)) then
           ! Awkward; in F2008 can be replaced by source/mold.
           allocate(self%buffer_2D(&
                size(self%packed_2D, 1),size(self%packed_2D, 2)))
           self%buffer_2D = 0._r8
        end if
        if (micro_mg_data_packer_use_native_impl) then
           self%buffer_2D = self%buffer_2D + self%packed_2D
        else
           call micro_mg_data_accumulate_2D_codon_wrap(self%packed_2D, self%buffer_2D)
        end if
     case default
        call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
             " Unsupported rank for MGFieldPostProc accumulation.")
     end select
  case default
     call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
          " Unrecognized MGFieldPostProc accumulation method.")
  end select

end subroutine MGFieldPostProc_accumulate

subroutine MGFieldPostProc_process_and_unpack(self, packer)
  class(MGFieldPostProc), intent(inout) :: self
  class(MGPacker), intent(in) :: packer

  call micro_mg_data_select_packer_impl()
  if (.not. micro_mg_data_packer_use_native_impl) then
     if (mgfieldpostproc_process_and_unpack_codon(1_c_int64_t) == 0_c_int64_t) return
  end if

  select case (self%accum_method)
  case (accum_null)
     ! "Null" method just leaves the value as the last time step, so don't
     ! actually need to do anything.
  case (accum_mean)
     select case (self%rank)
     case (1)
        SHR_ASSERT(associated(self%packed_1D), errMsg(__FILE__, __LINE__))
        self%packed_1D = self%buffer_1D/self%num_steps
     case (2)
        SHR_ASSERT(associated(self%packed_2D), errMsg(__FILE__, __LINE__))
        self%packed_2D = self%buffer_2D/self%num_steps
     case default
        call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
             " Unsupported rank for MGFieldPostProc accumulation.")
     end select
  case default
     call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
          " Unrecognized MGFieldPostProc accumulation method.")
  end select

  call self%unpack_only(packer)

end subroutine MGFieldPostProc_process_and_unpack

subroutine MGFieldPostProc_unpack_only(self, packer)
  class(MGFieldPostProc), intent(inout) :: self
  class(MGPacker), intent(in) :: packer

  call micro_mg_data_select_packer_impl()
  if (.not. micro_mg_data_packer_use_native_impl) then
     if (mgfieldpostproc_unpack_only_codon(1_c_int64_t) == 0_c_int64_t) return
  end if

  select case (self%rank)
  case (1)
     SHR_ASSERT(associated(self%unpacked_1D), errMsg(__FILE__, __LINE__))
     self%unpacked_1D = packer%unpack(self%packed_1D, self%fillvalue)
  case (2)
     SHR_ASSERT(associated(self%unpacked_2D), errMsg(__FILE__, __LINE__))
     self%unpacked_2D = packer%unpack(self%packed_2D, self%fillvalue)
  case default
     call shr_sys_abort(errMsg(__FILE__, __LINE__) // &
          " Unsupported rank for MGFieldPostProc unpacking.")
  end select

end subroutine MGFieldPostProc_unpack_only

#include "dynamic_vector_procdef.inc"

function new_MGPostProc(packer) result(post_proc)
  type(MGPacker), intent(in) :: packer

  type(MGPostProc) :: post_proc
  integer(c_int64_t), target :: plan(1)

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     post_proc%packer = packer
     call post_proc%field_procs%clear()
     return
  end if

  call new_mgpostproc_codon(c_loc(plan))
  if (plan(1) == 0_c_int64_t) return
  call micro_mg_data_log_new_mgpostproc_direct()

  post_proc%packer = packer
  call post_proc%field_procs%clear()

end function new_MGPostProc

! Can't use the same intent(out) trick, because PGI doesn't get the
! recursive deallocation right.
subroutine MGPostProc_finalize(self)
  class(MGPostProc), intent(inout) :: self

  integer :: i
  integer(c_int64_t), target :: plan(3)

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     call self%packer%finalize()
     do i = 1, self%field_procs%vsize()
        call self%field_procs%data(i)%finalize()
     end do
     call self%field_procs%clear()
     call self%field_procs%shrink_to_fit()
     return
  end if

  call mgpostproc_finalize_codon(int(self%field_procs%vsize(), c_int64_t), c_loc(plan))
  if (plan(1) == 0_c_int64_t) return
  call micro_mg_data_log_mgpostproc_finalize_direct()

  if (plan(2) /= 0_c_int64_t) call self%packer%finalize()
  do i = 1, int(plan(3))
     call self%field_procs%data(i)%finalize()
  end do
  call self%field_procs%clear()
  call self%field_procs%shrink_to_fit()

end subroutine MGPostProc_finalize

subroutine add_field_1D(self, unpacked_ptr, packed_ptr, fillvalue, &
     accum_method)
  class(MGPostProc), intent(inout) :: self
  real(r8), pointer, intent(in) :: unpacked_ptr(:)
  real(r8), pointer, intent(in) :: packed_ptr(:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  integer(c_int64_t), target :: plan(2)
  real(r8), target :: fill_plan(1)
  integer(c_int64_t) :: fill_present_c, accum_present_c, accum_in_c
  real(c_double) :: fill_in_c

  fill_present_c = merge(1_c_int64_t, 0_c_int64_t, present(fillvalue))
  accum_present_c = merge(1_c_int64_t, 0_c_int64_t, present(accum_method))
  fill_in_c = 0.0_c_double
  if (present(fillvalue)) fill_in_c = real(fillvalue, c_double)
  accum_in_c = 0_c_int64_t
  if (present(accum_method)) accum_in_c = int(accum_method, c_int64_t)

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     fill_plan(1) = real(fill_in_c, r8)
     plan(2) = int(accum_mean, c_int64_t)
     if (present(accum_method)) plan(2) = int(accum_method, c_int64_t)
     call self%field_procs%push_back(MGFieldPostProc(unpacked_ptr, &
          packed_ptr, fill_plan(1), int(plan(2))))
     return
  end if

  call add_field_1d_codon(fill_present_c, accum_present_c, fill_in_c, accum_in_c, &
       c_loc(plan), c_loc(fill_plan))
  if (plan(1) == 0_c_int64_t) return
  call micro_mg_data_log_add_field_1d_direct()

  call self%field_procs%push_back(MGFieldPostProc(unpacked_ptr, &
       packed_ptr, fill_plan(1), int(plan(2))))

end subroutine add_field_1D

subroutine add_field_2D(self, unpacked_ptr, packed_ptr, fillvalue, &
     accum_method)
  class(MGPostProc), intent(inout) :: self
  real(r8), pointer, intent(in) :: unpacked_ptr(:,:)
  real(r8), pointer, intent(in) :: packed_ptr(:,:)
  real(r8), intent(in), optional :: fillvalue
  integer, intent(in), optional :: accum_method
  integer(c_int64_t), target :: plan(2)
  real(r8), target :: fill_plan(1)
  integer(c_int64_t) :: fill_present_c, accum_present_c, accum_in_c
  real(c_double) :: fill_in_c

  fill_present_c = merge(1_c_int64_t, 0_c_int64_t, present(fillvalue))
  accum_present_c = merge(1_c_int64_t, 0_c_int64_t, present(accum_method))
  fill_in_c = 0.0_c_double
  if (present(fillvalue)) fill_in_c = real(fillvalue, c_double)
  accum_in_c = 0_c_int64_t
  if (present(accum_method)) accum_in_c = int(accum_method, c_int64_t)

  call micro_mg_data_select_packer_impl()
  if (micro_mg_data_packer_use_native_impl) then
     fill_plan(1) = real(fill_in_c, r8)
     plan(2) = int(accum_mean, c_int64_t)
     if (present(accum_method)) plan(2) = int(accum_method, c_int64_t)
     call self%field_procs%push_back(MGFieldPostProc(unpacked_ptr, &
          packed_ptr, fill_plan(1), int(plan(2))))
     return
  end if

  call add_field_2d_codon(fill_present_c, accum_present_c, fill_in_c, accum_in_c, &
       c_loc(plan), c_loc(fill_plan))
  if (plan(1) == 0_c_int64_t) return
  call micro_mg_data_log_add_field_2d_direct()

  call self%field_procs%push_back(MGFieldPostProc(unpacked_ptr, &
       packed_ptr, fill_plan(1), int(plan(2))))

end subroutine add_field_2D

subroutine MGPostProc_accumulate(self)
  class(MGPostProc), intent(inout) :: self

  integer :: i

  call micro_mg_data_select_packer_impl()
  if (.not. micro_mg_data_packer_use_native_impl) then
     if (mgpostproc_accumulate_codon(1_c_int64_t) == 0_c_int64_t) return
  end if

  do i = 1, self%field_procs%vsize()
     call self%field_procs%data(i)%accumulate()
  end do

end subroutine MGPostProc_accumulate

subroutine MGPostProc_process_and_unpack(self)
  class(MGPostProc), intent(inout) :: self

  integer :: i

  call micro_mg_data_select_packer_impl()
  if (.not. micro_mg_data_packer_use_native_impl) then
     if (mgpostproc_process_and_unpack_codon(1_c_int64_t) == 0_c_int64_t) return
  end if

  do i = 1, self%field_procs%vsize()
     call self%field_procs%data(i)%process_and_unpack(self%packer)
  end do

end subroutine MGPostProc_process_and_unpack

subroutine MGPostProc_unpack_only(self)
  class(MGPostProc), intent(inout) :: self

  integer :: i

  call micro_mg_data_select_packer_impl()
  if (.not. micro_mg_data_packer_use_native_impl) then
     if (mgpostproc_unpack_only_codon(1_c_int64_t) == 0_c_int64_t) return
  end if

  do i = 1, self%field_procs%vsize()
     call self%field_procs%data(i)%unpack_only(self%packer)
  end do

end subroutine MGPostProc_unpack_only

! This is necessary only to work around Intel/PGI bugs.
subroutine MGPostProc_copy(lhs, rhs)
  class(MGPostProc), intent(out) :: lhs
  type(MGPostProc), intent(in) :: rhs

  call micro_mg_data_select_packer_impl()
  if (.not. micro_mg_data_packer_use_native_impl) then
     if (mgpostproc_copy_codon(1_c_int64_t) == 0_c_int64_t) return
  end if

  lhs%packer = rhs%packer
  lhs%field_procs = rhs%field_procs
end subroutine MGPostProc_copy

end module micro_mg_data
