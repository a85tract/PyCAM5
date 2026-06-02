#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module reduction_mod
  use kinds, only : real_kind
  implicit none
  private

  type, public :: ReductionBuffer_int_1d_t
     integer, dimension(:), pointer :: buf
     integer :: len=0
     integer :: ctr
  end type ReductionBuffer_int_1d_t

  type, public :: ReductionBuffer_r_1d_t
     real (kind=real_kind), dimension(:), pointer :: buf
     integer :: len=0
     integer :: ctr
  end type ReductionBuffer_r_1d_t

  type, public :: ReductionBuffer_ordered_1d_t
     real (kind=real_kind), dimension(:,:),pointer :: buf
     integer :: len=0
     integer :: ctr
  end type ReductionBuffer_ordered_1d_t

  public :: ParallelMin,ParallelMax

  !type (ReductionBuffer_ordered_1d_t), public :: red_sum
  type (ReductionBuffer_int_1d_t),       public :: red_max_int
  type (ReductionBuffer_int_1d_t),       public :: red_sum_int
  type (ReductionBuffer_r_1d_t),       public :: red_sum
  type (ReductionBuffer_r_1d_t),       public :: red_max,red_min
  type (ReductionBuffer_r_1d_t),       public :: red_flops,red_timer

  !JMD new addition
#ifndef Darwin
  SAVE red_sum,red_max,red_min,red_flops,red_timer,red_max_int,red_sum_int
#endif
  interface ParallelMin
     module procedure ParallelMin1d
     module procedure ParallelMin0d
  end interface
  interface ParallelMax
     module procedure ParallelMax1d_int
     module procedure ParallelMax2d_int
     module procedure ParallelMax1d
     module procedure ParallelMax0d
     module procedure ParallelMax0d_int
  end interface

  interface pmax_mt
     module procedure pmax_mt_int_1d
     module procedure pmax_mt_r_1d
  end interface

  interface pmin_mt
     module procedure pmin_mt_r_1d
  end interface

  interface InitReductionBuffer
     module procedure InitReductionBuffer_int_1d
     module procedure InitReductionBuffer_r_1d
     module procedure InitReductionBuffer_ordered_1d
  end interface

  public :: InitReductionBuffer
  public :: pmax_mt, pmin_mt
  public :: ElementSum_1d

contains

  function ParallelMin1d(data,hybrid) result(pmin)
    use hybrid_mod, only : hybrid_t
    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    use cam_logfile, only : iulog
    implicit none
    interface
      function parallelmin1d_local_codon(data_p, len_c) result(pmin_c) &
           bind(c, name='parallelmin1d_local_codon')
        import :: c_double, c_int64_t, c_ptr
        type(c_ptr), value :: data_p
        integer(c_int64_t), value :: len_c
        real(c_double) :: pmin_c
      end function parallelmin1d_local_codon
    end interface
    real(kind=real_kind), intent(in), target :: data(:)
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: pmin

    real(kind=real_kind)                :: tmp(1)
    logical, save :: proof_seen = .false.


    tmp(1) = parallelmin1d_local_codon(c_loc(data(1)), int(size(data), c_int64_t))
    call pmin_mt(red_min,tmp,1,hybrid)
    pmin = red_min%buf(1)
    if (.not. proof_seen) then
       write(iulog,*) 'parallelmin1d implementation = codon'
       proof_seen = .true.
    endif

  end function ParallelMin1d

  function ParallelMin0d(data,hybrid) result(pmin)
    use hybrid_mod, only : hybrid_t
    use iso_c_binding, only : c_double
    use cam_logfile, only : iulog
    implicit none
    interface
      function parallelmin0d_local_codon(data_c) result(pmin_c) &
           bind(c, name='parallelmin0d_local_codon')
        import :: c_double
        real(c_double), value :: data_c
        real(c_double) :: pmin_c
      end function parallelmin0d_local_codon
    end interface
    real(kind=real_kind), intent(in)    :: data
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: pmin
    real(kind=real_kind)                :: tmp(1)
    logical, save :: proof_seen = .false.
    tmp(1) = parallelmin0d_local_codon(data)
    call pmin_mt(red_min,tmp,1,hybrid)
    pmin = red_min%buf(1)
    if (.not. proof_seen) then
       write(iulog,*) 'parallelmin0d implementation = codon'
       proof_seen = .true.
    endif

  end function ParallelMin0d
  !==================================================
  function ParallelMax2d_int(data, n, m, hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    integer, intent(in)                 :: n,m
    integer, intent(in), dimension(n,m) :: data
    type (hybrid_t), intent(in)         :: hybrid
    integer, dimension(n,m)             :: pmax
    integer, dimension(n*m)             :: tmp
    integer :: ierr,i,j
    do i=1,n 
      do j=1,m
        tmp(i+(j-1)*n) = data(i,j)
      enddo 
    enddo 
    call pmax_mt(red_max_int,tmp,n*m,hybrid)
    do i=1,n 
      do j=1,m
        pmax(i,j) = red_max_int%buf(i+(j-1)*n) 
      enddo 
    enddo 
  end function ParallelMax2d_int

  function ParallelMax1d_int(data, len, hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    integer, intent(in)                 :: len
    integer, intent(in), dimension(len) :: data
    type (hybrid_t), intent(in)         :: hybrid
    integer, dimension(len)             :: pmax, tmp
    integer :: ierr

    tmp = data(:)
    call pmax_mt(red_max_int,tmp,len,hybrid)
    pmax(:) = red_max_int%buf(1:len)

  end function ParallelMax1d_int
  function ParallelMax1d(data,hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
    use cam_logfile, only : iulog
    implicit none
    interface
      function parallelmax1d_local_codon(data_p, len_c) result(pmax_c) &
           bind(c, name='parallelmax1d_local_codon')
        import :: c_double, c_int64_t, c_ptr
        type(c_ptr), value :: data_p
        integer(c_int64_t), value :: len_c
        real(c_double) :: pmax_c
      end function parallelmax1d_local_codon
    end interface
    real(kind=real_kind), intent(in), target :: data(:)
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: pmax

    real(kind=real_kind)                :: tmp(1)
    logical, save :: proof_seen = .false.


    tmp(1) = parallelmax1d_local_codon(c_loc(data(1)), int(size(data), c_int64_t))
    call pmax_mt(red_max,tmp,1,hybrid)
    pmax = red_max%buf(1)
    if (.not. proof_seen) then
       write(iulog,*) 'parallelmax1d implementation = codon'
       proof_seen = .true.
    endif

  end function ParallelMax1d
  function ParallelMax0d(data,hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    use iso_c_binding, only : c_double
    use cam_logfile, only : iulog
    implicit none
    interface
      function parallelmax0d_local_codon(data_c) result(pmax_c) &
           bind(c, name='parallelmax0d_local_codon')
        import :: c_double
        real(c_double), value :: data_c
        real(c_double) :: pmax_c
      end function parallelmax0d_local_codon
    end interface
    real(kind=real_kind), intent(in)    :: data
    type (hybrid_t),      intent(in)    :: hybrid
    real(kind=real_kind)                :: pmax
    real(kind=real_kind)                :: tmp(1)
    logical, save :: proof_seen = .false.

    tmp(1)=parallelmax0d_local_codon(data)

    call pmax_mt(red_max,tmp,1,hybrid)
    pmax = red_max%buf(1)
    if (.not. proof_seen) then
       write(iulog,*) 'parallelmax0d implementation = codon'
       proof_seen = .true.
    endif

  end function ParallelMax0d
  function ParallelMax0d_int(data,hybrid) result(pmax)
    use hybrid_mod, only : hybrid_t
    implicit none
    integer             , intent(in)    :: data
    type (hybrid_t),      intent(in)    :: hybrid
    integer                             :: pmax
    integer                             :: tmp(1)

    tmp(1)=data

    call pmax_mt(red_max_int,tmp,1,hybrid)
    pmax = red_max_int%buf(1)

  end function ParallelMax0d_int
  !==================================================
  subroutine InitReductionBuffer_int_1d(red,len)
    use parallel_mod, only: abortmp
    use thread_mod, only: omp_get_num_threads
    use iso_c_binding, only : c_int64_t, c_int, c_loc, c_ptr
    use cam_logfile, only : iulog
    integer, intent(in)           :: len
    type (ReductionBuffer_int_1d_t), intent(inout), target :: red
    integer(c_int), target :: new_len_c, new_ctr_c
    integer(c_int64_t) :: realloc_c
    logical, save :: proof_seen = .false.

    interface
       function initreductionbuffer_int_1d_codon(current_len_c, requested_len_c, len_p, ctr_p) result(realloc_c) &
            bind(c, name='initreductionbuffer_int_1d_codon')
         import :: c_int64_t, c_ptr
         integer(c_int64_t), value :: current_len_c
         integer(c_int64_t), value :: requested_len_c
         type(c_ptr), value :: len_p
         type(c_ptr), value :: ctr_p
         integer(c_int64_t) :: realloc_c
       end function initreductionbuffer_int_1d_codon
    end interface

#define SE_MISC_TAG 21
#define SE_MISC_LABEL 'reduction_mod'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

    if (omp_get_num_threads()>1) then
       call abortmp("Error: attempt to allocate reduction buffer in threaded region")
    endif

    realloc_c = initreductionbuffer_int_1d_codon(int(red%len, c_int64_t), int(len, c_int64_t), &
         c_loc(new_len_c), c_loc(new_ctr_c))
    if (realloc_c /= 0_c_int64_t) then
       if (red%len>0) deallocate(red%buf)
       red%len  = int(new_len_c)
       allocate(red%buf(red%len))
       red%buf  = 0
    endif
    red%ctr = int(new_ctr_c)
    if (.not. proof_seen) then
       write(iulog,*) 'initreductionbuffer_int_1d implementation = codon'
       proof_seen = .true.
    endif

  end subroutine InitReductionBuffer_int_1d
  !****************************************************************
  subroutine InitReductionBuffer_r_1d(red,len)
    use parallel_mod, only: abortmp
    use thread_mod, only: omp_get_num_threads
    use iso_c_binding, only : c_int64_t, c_int, c_loc, c_ptr
    use cam_logfile, only : iulog
    integer, intent(in)           :: len
    type (ReductionBuffer_r_1d_t), intent(inout), target :: red
    integer(c_int), target :: new_len_c, new_ctr_c
    integer(c_int64_t) :: realloc_c
    logical, save :: proof_seen = .false.

    interface
       function initreductionbuffer_r_1d_codon(current_len_c, requested_len_c, len_p, ctr_p) result(realloc_c) &
            bind(c, name='initreductionbuffer_r_1d_codon')
         import :: c_int64_t, c_ptr
         integer(c_int64_t), value :: current_len_c
         integer(c_int64_t), value :: requested_len_c
         type(c_ptr), value :: len_p
         type(c_ptr), value :: ctr_p
         integer(c_int64_t) :: realloc_c
       end function initreductionbuffer_r_1d_codon
    end interface

    if (omp_get_num_threads()>1) then
       call abortmp("Error: attempt to allocate reduction buffer in threaded region")
    endif

    realloc_c = initreductionbuffer_r_1d_codon(int(red%len, c_int64_t), int(len, c_int64_t), &
         c_loc(new_len_c), c_loc(new_ctr_c))
    if (realloc_c /= 0_c_int64_t) then
       if (red%len>0) deallocate(red%buf)
       red%len  = int(new_len_c)
       allocate(red%buf(red%len))
       red%buf  = 0.0D0
    endif
    red%ctr = int(new_ctr_c)
    if (.not. proof_seen) then
       write(iulog,*) 'initreductionbuffer_r_1d implementation = codon'
       proof_seen = .true.
    endif
  end subroutine InitReductionBuffer_r_1d
  !****************************************************************
  subroutine InitReductionBuffer_ordered_1d(red,len,nthread)
    use parallel_mod, only: abortmp
    use thread_mod, only: omp_get_num_threads
    use iso_c_binding, only : c_int64_t, c_int, c_loc, c_ptr
    use cam_logfile, only : iulog
    integer, intent(in)           :: len
    integer, intent(in)           :: nthread
    type (ReductionBuffer_ordered_1d_t), intent(inout), target :: red
    integer(c_int), target :: new_len_c, new_ctr_c
    integer(c_int64_t) :: realloc_c
    logical, save :: proof_seen = .false.

    interface
       function initreductionbuffer_ordered_1d_codon(current_len_c, requested_len_c, nthread_c, len_p, ctr_p) &
            result(realloc_c) bind(c, name='initreductionbuffer_ordered_1d_codon')
         import :: c_int64_t, c_ptr
         integer(c_int64_t), value :: current_len_c
         integer(c_int64_t), value :: requested_len_c
         integer(c_int64_t), value :: nthread_c
         type(c_ptr), value :: len_p
         type(c_ptr), value :: ctr_p
         integer(c_int64_t) :: realloc_c
       end function initreductionbuffer_ordered_1d_codon
    end interface

    if (omp_get_num_threads()>1) then
       call abortmp("Error: attempt to allocate reduction buffer in threaded region")
    endif

    realloc_c = initreductionbuffer_ordered_1d_codon(int(red%len, c_int64_t), int(len, c_int64_t), &
         int(nthread, c_int64_t), c_loc(new_len_c), c_loc(new_ctr_c))
    if (realloc_c /= 0_c_int64_t) then
       if (red%len>0) deallocate(red%buf)
       red%len  = int(new_len_c)
       allocate(red%buf(len,nthread+1))
       red%buf  = 0.0D0
    endif
    red%ctr = int(new_ctr_c)
    if (.not. proof_seen) then
       write(iulog,*) 'initreductionbuffer_ordered_1d implementation = codon'
       proof_seen = .true.
    endif
  end subroutine InitReductionBuffer_ordered_1d

  ! =======================================
  ! pmax_mt:
  !
  ! thread safe, parallel reduce maximum
  ! of a one dimensional reduction vector
  ! =======================================

  subroutine pmax_mt_int_1d(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t
#ifdef _MPI
    use parallel_mod, only: mpi_min, mpi_max, mpiinteger_t,abortmp
#else
    use parallel_mod, only: abortmp
#endif

    type (ReductionBuffer_int_1d_t)   :: red       ! shared memory reduction buffer struct
    integer,               intent(in) :: len       ! buffer length
    integer, intent(inout)            :: redp(len) ! thread private vector of partial sum
    type (hybrid_t),       intent(in) :: hybrid    ! parallel handle

    ! Local variables
#ifdef _MPI
    integer ierr
#endif

    integer  :: k
    if (len>red%len) call abortmp('ERROR: threadsafe reduction buffer too small')


#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
    !$OMP CRITICAL (CRITMAX)
#endif
    if (red%ctr == 0) red%buf(1:len)= -9999
    if (red%ctr < hybrid%NThreads) then
       do k=1,len
          red%buf(k)=MAX(red%buf(k),redp(k))
       enddo
       red%ctr=red%ctr+1
    end if
    if (red%ctr == hybrid%NThreads) red%ctr=0
#if (defined HORIZ_OPENMP)
    !$OMP END CRITICAL (CRITMAX)
#endif
#ifdef _MPI
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1),redp,len,MPIinteger_t, &
            MPI_MAX,hybrid%par%comm,ierr)

       red%buf(1:len)=redp(1:len)
    end if
#endif
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif


  end subroutine pmax_mt_int_1d
  
  subroutine pmax_mt_r_1d(red,redp,len,hybrid)
    use hybrid_mod, only : hybrid_t
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
#ifdef _MPI
    use parallel_mod, only: mpi_min, mpi_max, mpireal_t,abortmp
#else
    use parallel_mod, only: abortmp
#endif
    use cam_logfile, only : iulog
    interface
      subroutine reduction_max_r_local_codon(buf_p, ctr_p, redp_p, len_c, nthreads_c) &
           bind(c, name='reduction_max_r_local_codon')
        import :: c_int64_t, c_ptr
        type(c_ptr), value :: buf_p, ctr_p, redp_p
        integer(c_int64_t), value :: len_c, nthreads_c
      end subroutine reduction_max_r_local_codon
    end interface

    type (ReductionBuffer_r_1d_t), target :: red ! shared memory reduction buffer struct
    real (kind=real_kind), intent(inout), target :: redp(:) ! thread private vector of partial sum
    integer,               intent(in) :: len     ! buffer length
    type (hybrid_t),       intent(in) :: hybrid  ! parallel handle

    ! Local variables
#ifdef _MPI
    integer ierr
#endif

    logical, save :: proof_seen = .false.
    if (len>red%len) call abortmp('ERROR: threadsafe reduction buffer too small')

#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
    !$OMP CRITICAL (CRITMAX)
#endif
    call reduction_max_r_local_codon(c_loc(red%buf(1)), c_loc(red%ctr), c_loc(redp(1)), &
         int(len, c_int64_t), int(hybrid%NThreads, c_int64_t))
    if (.not. proof_seen) then
       write(iulog,*) 'pmax_mt_r_1d implementation = codon'
       proof_seen = .true.
    endif
#if (defined HORIZ_OPENMP)
    !$OMP END CRITICAL (CRITMAX)
#endif
#ifdef _MPI
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1),redp,len,MPIreal_t, &
            MPI_MAX,hybrid%par%comm,ierr)

       red%buf(1:len)=redp(1:len)
    end if
#endif
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif


  end subroutine pmax_mt_r_1d

  ! =======================================
  ! pmin_mt:
  !
  ! thread safe, parallel reduce maximum
  ! of a one dimensional reduction vector
  ! =======================================

  subroutine pmin_mt_r_1d(red,redp,len,hybrid)
    use kinds, only : int_kind
    use hybrid_mod, only : hybrid_t
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
#ifdef _MPI
    use parallel_mod, only: mpi_min, mpireal_t,abortmp
#else
    use parallel_mod, only: abortmp
#endif
    use cam_logfile, only : iulog
    interface
      subroutine reduction_min_r_local_codon(buf_p, ctr_p, redp_p, len_c, nthreads_c) &
           bind(c, name='reduction_min_r_local_codon')
        import :: c_int64_t, c_ptr
        type(c_ptr), value :: buf_p, ctr_p, redp_p
        integer(c_int64_t), value :: len_c, nthreads_c
      end subroutine reduction_min_r_local_codon
    end interface

    type (ReductionBuffer_r_1d_t), target :: red ! shared memory reduction buffer struct
    real (kind=real_kind), intent(inout), target :: redp(:) ! thread private vector of partial sum
    integer,               intent(in) :: len     ! buffer length
    type (hybrid_t),       intent(in) :: hybrid  ! parallel handle

    ! Local variables

#ifdef _MPI
    integer ierr
#endif
    logical, save :: proof_seen = .false.

    if (len>red%len) call abortmp('ERROR: threadsafe reduction buffer too small')

#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
    !$OMP CRITICAL (CRITMAX)
#endif
    call reduction_min_r_local_codon(c_loc(red%buf(1)), c_loc(red%ctr), c_loc(redp(1)), &
         int(len, c_int64_t), int(hybrid%NThreads, c_int64_t))
    if (.not. proof_seen) then
       write(iulog,*) 'pmin_mt_r_1d implementation = codon'
       proof_seen = .true.
    endif
#if (defined HORIZ_OPENMP)
    !$OMP END CRITICAL (CRITMAX)
#endif
#ifdef _MPI
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if (hybrid%ithr==0) then

       call MPI_Allreduce(red%buf(1),redp,len,MPIreal_t, &
            MPI_MIN,hybrid%par%comm,ierr)

       red%buf(1:len)=redp(1:len)
    end if
#endif
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif


  end subroutine pmin_mt_r_1d

  subroutine ElementSum_1d(res,variable,type,hybrid)
    use hybrid_mod, only : hybrid_t
    use dimensions_mod, only : nelem
#ifdef _MPI
  use parallel_mod, only : ORDERED, mpireal_t, mpi_min, mpi_max, mpi_sum, mpi_success
#else
  use parallel_mod, only : ORDERED
#endif
    implicit none

    ! ==========================
    !     Arguments
    ! ==========================
    real(kind=real_kind),intent(out) :: res
    real(kind=real_kind),intent(in)  :: variable(:)
    integer,intent(in)               :: type
    type (hybrid_t), intent(in)      :: hybrid 

    ! ==========================
    !       Local Variables
    ! ==========================

    !
    ! Note this is a real kludge here since it may be used for
    !  arrays of size other then nelem
    !

    integer                          :: i
#if 0
    real(kind=real_kind),allocatable :: Global(:)
    real(kind=real_kind),allocatable :: buffer(:)
#endif

#ifdef _MPI
    integer                           :: errorcode,errorlen
    character*(80) errorstring

    real(kind=real_kind)             :: local_sum
    integer                          :: ierr
#endif

#ifdef _MPI
    if(hybrid%ithr == 0) then 
#if 0
       if(type == ORDERED) then
          allocate(buffer(nelem))
          call MPI_Gatherv(variable,nelemd,MPIreal_t,buffer, &
               recvcount,displs,MPIreal_t,hybrid%par%root, &
               hybrid%par%comm,ierr)
          if(ierr .ne. MPI_SUCCESS) then 
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'ElementSum_1d: Error after call to MPI_Gatherv: ',errorstring
          endif

          if(hybrid%par%masterproc) then
             allocate(Global(nelem))
             do ip=1,hybrid%par%nprocs
                nelemr = recvcount(ip)
                disp   = displs(ip)
                do ie=1,nelemr
                   ig = Schedule(ip)%Local2Global(ie)
                   Global(ig) = buffer(disp+ie)
                enddo
             enddo
             ! ===========================
             !  Perform the ordererd sum
             ! ===========================
             res = 0.0d0
             do i=1,nelem
                res = res + Global(i)
             enddo
             deallocate(Global)
          endif
          ! =============================================
          !  Broadcast the results back everybody
          ! =============================================
          call MPI_Bcast(res,1,MPIreal_t,hybrid%par%root, &
               hybrid%par%comm,ierr)
          if(ierr .ne. MPI_SUCCESS) then 
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'ElementSum_1d: Error after call to MPI_Bcast: ',errorstring
          endif

          deallocate(buffer)
       else
#endif
          local_sum=SUM(variable)
          call MPI_Barrier(hybrid%par%comm,ierr)

          call MPI_Allreduce(local_sum,res,1,MPIreal_t, &
               MPI_SUM,hybrid%par%comm,ierr)
          if(ierr .ne. MPI_SUCCESS) then 
             errorcode=ierr
             call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
             print *,'ElementSum_1d: Error after call to MPI_Allreduce: ',errorstring
          endif
#if 0
       endif
#endif
    endif
#else
    if(hybrid%ithr == 0) then 
       if(type == ORDERED) then
          ! ===========================
          !  Perform the ordererd sum
          ! ===========================
          res = 0.0d0
          do i=1,nelem
             res = res + variable(i)
          enddo
       else
          res=SUM(variable)
       endif
    endif
#endif

  end subroutine ElementSum_1d

end module reduction_mod
