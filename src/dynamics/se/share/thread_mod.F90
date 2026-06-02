#ifdef HAVE_CONFIG_H  
#include "config.h"
#endif

module thread_mod

#ifdef _OPENMP
  use omp_lib, only: omp_get_thread_num, &
       omp_in_parallel, &
       omp_set_num_threads, &
       omp_get_max_threads, &
       omp_get_num_threads, &
       omp_get_nested,      &
       omp_set_nested
#endif

  implicit none
  private

  integer, public :: max_num_threads=1    ! maximum number of OpenMP threads
  integer, public :: horz_num_threads=1   ! number of OpenMP threads in horizontal
  integer, public :: vert_num_threads=1   ! number of OpenMP threads in vertical
  integer, public :: tracer_num_threads=1 ! number of OpenMP threads in tracers

  public :: omp_get_thread_num
  public :: omp_in_parallel
  public :: omp_set_num_threads
  public :: omp_get_max_threads
  public :: omp_get_num_threads
  public :: omp_get_nested
  public :: omp_set_nested
  public :: initomp
contains

#ifndef _OPENMP
  function omp_get_thread_num() result(ithr)
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog
    interface
      function se_omp_get_thread_num_codon() result(ithr_c) bind(c, name='se_omp_get_thread_num_codon')
        import :: c_int64_t
        integer(c_int64_t) :: ithr_c
      end function se_omp_get_thread_num_codon
    end interface
    integer ithr
    logical, save :: proof_seen = .false.
    ithr=int(se_omp_get_thread_num_codon())
    if (.not. proof_seen) then
       write(iulog,*) 'omp_get_thread_num implementation = codon'
       proof_seen = .true.
    endif
  end function omp_get_thread_num

  function omp_get_num_threads() result(ithr)
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog
    interface
      function se_omp_get_num_threads_codon() result(ithr_c) bind(c, name='se_omp_get_num_threads_codon')
        import :: c_int64_t
        integer(c_int64_t) :: ithr_c
      end function se_omp_get_num_threads_codon
    end interface
    integer ithr
    logical, save :: proof_seen = .false.
    ithr=int(se_omp_get_num_threads_codon())
    if (.not. proof_seen) then
       write(iulog,*) 'omp_get_num_threads implementation = codon'
       proof_seen = .true.
    endif
  end function omp_get_num_threads

  function omp_in_parallel() result(ans)
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog
    interface
      function se_omp_in_parallel_codon() result(ans_c) bind(c, name='se_omp_in_parallel_codon')
        import :: c_int64_t
        integer(c_int64_t) :: ans_c
      end function se_omp_in_parallel_codon
    end interface
    logical ans
    logical, save :: proof_seen = .false.
    ans=(se_omp_in_parallel_codon() /= 0_c_int64_t)
    if (.not. proof_seen) then
       write(iulog,*) 'omp_in_parallel implementation = codon'
       proof_seen = .true.
    endif
  end function omp_in_parallel

  subroutine omp_set_num_threads(NThreads)
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog
    interface
      function se_omp_set_num_threads_codon(nthreads_c) result(nthreads_out_c) &
           bind(c, name='se_omp_set_num_threads_codon')
        import :: c_int64_t
        integer(c_int64_t), value :: nthreads_c
        integer(c_int64_t) :: nthreads_out_c
      end function se_omp_set_num_threads_codon
    end interface
    integer Nthreads
    logical, save :: proof_seen = .false.

#define SE_MISC_TAG 33
#define SE_MISC_LABEL 'thread_mod'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

    NThreads=int(se_omp_set_num_threads_codon(int(NThreads, c_int64_t)))
    if (.not. proof_seen) then
       write(iulog,*) 'omp_set_num_threads implementation = codon'
       proof_seen = .true.
    endif
  end subroutine omp_set_num_threads

  integer function omp_get_max_threads()
    omp_get_max_threads=1
  end function omp_get_max_threads

  integer function omp_get_nested()
    omp_get_nested=0
  end function omp_get_nested

  integer function omp_set_nested()
    omp_set_nested=0
  end function omp_set_nested

  subroutine initomp
    max_num_threads = 1
  end subroutine initomp

#else
  subroutine initomp
    max_num_threads = omp_get_max_threads()
  end subroutine initomp
#endif

end module thread_mod
