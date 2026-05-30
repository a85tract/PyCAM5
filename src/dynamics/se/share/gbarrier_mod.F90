module gbarrier_mod
  use gbarriertype_mod, only: gbarrier_t
  implicit none

  integer, parameter :: LOG2MAX = 6 
  integer, parameter :: MAXTHREADS = 64
  logical :: gbarrier_use_native_impl = .false.
  logical :: gbarrier_impl_selected = .false.
  logical :: gbarrier_init_proof_seen = .false.
  logical :: gbarrier_delete_proof_seen = .false.
  logical :: gbarrier_proof_seen = .false.

  public :: gbarrier_init
  public :: gbarrier_info
  public :: gbarrier

  contains

    subroutine gbarrier_init(barrier, nthreads)
      use iso_c_binding, only : c_int64_t
      use cam_logfile, only : iulog
      type (gbarrier_t), intent(out) :: barrier
      integer, intent(in) :: nthreads

      interface 
        subroutine gbarrier_initialize(c_barrier, nthreads) bind(C)
          use, intrinsic :: ISO_C_Binding, only: C_ptr, C_int
          implicit none

          type (C_ptr), intent(out) :: c_barrier
          integer (C_int), intent(in), value :: nthreads
        end subroutine gbarrier_initialize
      end interface
      interface
        subroutine gbarrier_init_codon(c_barrier, nthreads) bind(C, name="gbarrier_init_codon")
          use, intrinsic :: ISO_C_Binding, only: C_ptr, c_int64_t
          implicit none

          type (C_ptr), intent(out) :: c_barrier
          integer (c_int64_t), intent(in), value :: nthreads
        end subroutine gbarrier_init_codon
      end interface

#define SE_MISC_TAG 35
#define SE_MISC_LABEL 'gbarrier_mod'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

      call gbarrier_select_impl()
      if (gbarrier_use_native_impl) then
        call gbarrier_initialize(barrier%c_barrier, nthreads)
      else
        call gbarrier_init_codon(barrier%c_barrier, int(nthreads, c_int64_t))
        call gbarrier_write_proof('gbarrier_init', gbarrier_init_proof_seen)
      endif
    end subroutine gbarrier_init

    subroutine gbarrier_delete(barrier)
      type (gbarrier_t), intent(in) :: barrier

      interface
        subroutine gbarrier_free(c_barrier) bind(C)
          use, intrinsic :: ISO_C_Binding, only: C_ptr
          implicit none

          type (C_ptr), intent(in) :: c_barrier
        end subroutine gbarrier_free
      end interface
      interface
        subroutine gbarrier_delete_codon(c_barrier) bind(C, name="gbarrier_delete_codon")
          use, intrinsic :: ISO_C_Binding, only: C_ptr
          implicit none

          type (C_ptr), intent(in) :: c_barrier
        end subroutine gbarrier_delete_codon
      end interface

      call gbarrier_select_impl()
      if (gbarrier_use_native_impl) then
        call gbarrier_free(barrier%c_barrier)
      else
        call gbarrier_delete_codon(barrier%c_barrier)
        call gbarrier_write_proof('gbarrier_delete', gbarrier_delete_proof_seen)
      endif
    end subroutine gbarrier_delete

    subroutine gbarrier_info(barrier)
      type (gbarrier_t), intent(in) :: barrier

      interface
        subroutine gbarrier_print(c_barrier) bind(C)
          use, intrinsic :: ISO_C_Binding, only: C_ptr
          implicit none
          type (C_ptr), value :: c_barrier
        end subroutine gbarrier_print
      end interface

      call gbarrier_print(barrier%c_barrier)
    end subroutine gbarrier_info


    subroutine gbarrier(barrier, threadID)
      use iso_c_binding, only : c_int64_t
      type (gbarrier_t), intent(in) :: barrier
      integer, intent(in) :: threadID

      interface 
        subroutine gbarrier_synchronize(c_barrier, thread) bind(C)
          use, intrinsic :: ISO_C_Binding, only: C_ptr, C_int
          implicit none

          type (C_ptr), intent(in), value :: c_barrier
          integer (C_int), intent(in), value :: thread
        end subroutine gbarrier_synchronize
      end interface
      interface
        subroutine gbarrier_synchronize_codon(c_barrier, thread) bind(C, name="gbarrier_synchronize_codon")
          use, intrinsic :: ISO_C_Binding, only: C_ptr, c_int64_t
          implicit none

          type (C_ptr), intent(in), value :: c_barrier
          integer (c_int64_t), intent(in), value :: thread
        end subroutine gbarrier_synchronize_codon
      end interface

      call gbarrier_select_impl()
      if (gbarrier_use_native_impl) then
        call gbarrier_synchronize(barrier%c_barrier, threadID)
      else
        call gbarrier_synchronize_codon(barrier%c_barrier, int(threadID, c_int64_t))
        if (threadID == 0) then
          call gbarrier_write_proof('gbarrier', gbarrier_proof_seen)
        endif
      endif
    end subroutine gbarrier

    subroutine gbarrier_select_impl()
      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (gbarrier_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('GBARRIER_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
        do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          endif
        enddo
        gbarrier_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
        gbarrier_use_native_impl = .false.
      endif

      gbarrier_impl_selected = .true.
    end subroutine gbarrier_select_impl

    subroutine gbarrier_write_proof(routine_name, proof_seen)
      use cam_logfile, only : iulog
      character(len=*), intent(in) :: routine_name
      logical, intent(inout) :: proof_seen

      if (.not. proof_seen) then
        write(iulog,*) trim(routine_name)//' implementation = codon'
        proof_seen = .true.
      endif
    end subroutine gbarrier_write_proof

end module gbarrier_mod
