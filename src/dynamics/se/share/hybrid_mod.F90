#ifdef HAVE_CONFIG_H  
#include "config.h"
#endif

! ===========================================
! Module to support hybrid programming model
! hybrid_t is assumed to be a private struct
! ===========================================

module hybrid_mod

use parallel_mod  , only : parallel_t, copy_par
use thread_mod    , only : omp_set_num_threads, omp_get_thread_num 
use thread_mod    , only : horz_num_threads, vert_num_threads, tracer_num_threads
use dimensions_mod, only : nlev, qsize, ntrac
use iso_c_binding, only : c_int64_t, c_loc, c_ptr

implicit none
private

  type, private :: hybridold_t
     type (parallel_t)    :: par
     integer              :: ithr
     integer              :: nthreads
     logical              :: masterthread
  end type

  type, private :: hybrid_p
     integer :: ibeg, iend
     integer :: kbeg, kend
     integer :: qbeg, qend
  end type

  type, public :: hybrid_t
     type (parallel_t) :: par
     integer           :: ithr
     integer           :: localsense
     integer           :: nthreads
     integer           :: ibeg, iend
     integer           :: kbeg, kend
     integer           :: qbeg, qend
     logical           :: masterthread
  end type 

  integer, allocatable, target :: work_pool_horz(:,:)
  integer, allocatable, target :: work_pool_vert(:,:)
  integer, allocatable, target :: work_pool_trac(:,:)

  integer :: nelemd_save
  logical :: init_ranges = .true.
  integer :: region_num_threads
  character(len=64) :: region_name

  private :: hybrid_create
  public :: PrintHybrid
  public :: set_region_num_threads
  private :: set_loop_ranges
  public :: get_loop_ranges
  public :: init_loop_ranges
  public :: threadOwnsTracer, threadOwnsVertlevel
  public :: config_thread_region
  public :: get_number_threads

  interface config_thread_region 
      module procedure config_thread_region_par
      module procedure config_thread_region_hybrid
  end interface
  interface PrintHybrid 
!      module procedure PrintHybridold
      module procedure PrintHybridnew
  end interface

contains

  subroutine PrintHybridold(hybt,hybp,vname)
    type (hybridold_t) :: hybt
    type (hybrid_p) :: hybp
    character(len=*) :: vname
     
    write(*,21) vname, hybt%par%rank, hybt%ithr, hybt%nthreads, &
                hybp%ibeg, hybp%iend,hybp%kbeg,hybp%kend, &
                hybp%qbeg, hybp%qend
21  format('PrintHybrid: (',a, ', rank: ',i8, ', ithrd: ',i4,',  nthreads: ',i4, &
           ', i{beg,end}: ',2(i4),', k{beg,end}: ',2(i4),', q{beg,end}: ',2(i4),')')

  end subroutine PrintHybridold

  subroutine PrintHybridnew(hybt,vname)
    type (hybrid_t) :: hybt
    character(len=*) :: vname

    write(*,21) vname, hybt%par%rank, hybt%ithr, hybt%nthreads, &
                hybt%ibeg, hybt%iend,hybt%kbeg,hybt%kend, &
                hybt%qbeg, hybt%qend
21  format('PrintHybrid: (',a, ', rank: ',i8, ', ithrd: ',i4,',  nthreads: ',i4, &
           ', i{beg,end}: ',2(i4),', k{beg,end}: ',2(i4),', q{beg,end}: ',2(i4),')')

  end subroutine PrintHybridnew

  
  function hybrid_create(par,ithr,nthreads) result(hybrid)
      type (parallel_t) , intent(in) :: par
      integer           , intent(in) :: ithr
      integer, optional , intent(in) :: nthreads
      type (hybridold_t)                :: hybrid

      hybrid%par      = par      ! relies on parallel_mod copy constructor
      hybrid%ithr     = ithr     
      if ( present(nthreads) ) then
        hybrid%nthreads = nthreads
      else
        hybrid%nthreads = region_num_threads
      endif

      hybrid%masterthread = (par%masterproc .and. ithr==0)

  end function hybrid_create 

  function config_thread_region_hybrid(old,region_name) result(new)
     type (hybrid_t), intent(in) :: old
     character(len=*), intent(in) :: region_name
     type (hybrid_t) :: new

     integer :: ithr
     integer :: kbeg_range, kend_range, qbeg_range, qend_range
     

     ithr = omp_get_thread_num()

     if ( TRIM(region_name) == 'serial') then 
         region_num_threads = 1
         new%ibeg = old%ibeg;      new%iend = old%iend
         new%kbeg = old%kbeg;      new%kend = old%kend
         new%qbeg = old%qbeg;      new%qend = old%qend
     endif
     if ( TRIM(region_name) == 'vertical') then
         region_num_threads = vert_num_threads
         call set_thread_ranges_1D ( work_pool_vert, kbeg_range, kend_range, ithr )
         new%ibeg = old%ibeg;      new%iend = old%iend
         new%kbeg = kbeg_range;    new%kend = kend_range
         new%qbeg = old%qbeg;      new%qend = old%qend
      endif

      if ( TRIM(region_name) == 'tracer' ) then
         region_num_threads = tracer_num_threads
         call set_thread_ranges_1D ( work_pool_trac, qbeg_range, qend_range, ithr)
         new%ibeg = old%ibeg;      new%iend = old%iend
         new%kbeg = old%kbeg;      new%kend = old%kend
         new%qbeg = qbeg_range;    new%qend = qend_range
      endif


      if ( TRIM(region_name) == 'vertical_and_tracer' ) then
         region_num_threads = vert_num_threads*tracer_num_threads
         call set_thread_ranges_2D ( work_pool_vert, work_pool_trac, kbeg_range, kend_range, &
                                                       qbeg_range, qend_range, ithr )
         new%ibeg = old%ibeg;      new%iend = old%iend
         new%kbeg = kbeg_range;    new%kend = kend_range
         new%qbeg = qbeg_range;    new%qend = qend_range
      endif


      new%par          = old%par      ! relies on parallel_mod copy constructor
      new%nthreads     = old%nthreads * region_num_threads
      if( region_num_threads .ne. 1 ) then 
          new%ithr         = old%ithr * region_num_threads + ithr
      else 
          new%ithr         = old%ithr
      endif
      new%masterthread = old%masterthread
      new%localsense   = old%localsense
!  Do we want to make this following call?
!      call omp_set_num_threads(new%nthreads)

  end function config_thread_region_hybrid

  function config_thread_region_par(par,region_name) result(hybrid)
      type (parallel_t) , intent(in) :: par
      character(len=*), intent(in) :: region_name
      type (hybrid_t)                :: hybrid
      ! local 
      integer    :: ithr
      integer    :: ibeg_range, iend_range
      integer    :: kbeg_range, kend_range
      integer    :: qbeg_range, qend_range
      integer    :: nthreads

      ithr            = omp_get_thread_num()

      if ( TRIM(region_name) == 'serial') then
         region_num_threads = 1
         call set_thread_ranges_1D ( work_pool_horz, ibeg_range, iend_range, ithr )
         hybrid%ibeg = 1;          hybrid%iend = nelemd_save
         hybrid%kbeg = 1;          hybrid%kend = nlev
         hybrid%qbeg = 1;          hybrid%qend = qsize
      endif

      if ( TRIM(region_name) == 'horizontal') then
         region_num_threads = horz_num_threads 
         call set_thread_ranges_1D ( work_pool_horz, ibeg_range, iend_range, ithr )
         hybrid%ibeg = ibeg_range; hybrid%iend = iend_range
         hybrid%kbeg = 1;          hybrid%kend = nlev
         hybrid%qbeg = 1;          hybrid%qend = qsize
      endif

      if ( TRIM(region_name) == 'vertical') then
         region_num_threads = vert_num_threads 
         call set_thread_ranges_1D ( work_pool_vert, kbeg_range, kend_range, ithr )
         hybrid%ibeg = 1;          hybrid%iend = nelemd_save
         hybrid%kbeg = kbeg_range; hybrid%kend = kend_range
         hybrid%qbeg = 1;          hybrid%qend = qsize
      endif
  
      if ( TRIM(region_name) == 'tracer' ) then
         region_num_threads = tracer_num_threads
         call set_thread_ranges_1D ( work_pool_trac, qbeg_range, qend_range, ithr)
         hybrid%ibeg = 1;          hybrid%iend = nelemd_save
         hybrid%kbeg = 1;          hybrid%kend = nlev
         hybrid%qbeg = qbeg_range; hybrid%qend = qend_range
      endif

    
      if ( TRIM(region_name) == 'vertical_and_tracer' ) then
         region_num_threads = vert_num_threads*tracer_num_threads
         call set_thread_ranges_2D ( work_pool_vert, work_pool_trac, kbeg_range, kend_range, &
                                                       qbeg_range, qend_range, ithr )
         hybrid%ibeg = 1;          hybrid%iend = nelemd_save
         hybrid%kbeg = kbeg_range; hybrid%kend = kend_range
         hybrid%qbeg = qbeg_range; hybrid%qend = qend_range
      endif
      call omp_set_num_threads(region_num_threads)

!      hybrid%par          = par      ! relies on parallel_mod copy constructor
      call copy_par(hybrid%par,par)
      hybrid%nthreads     = region_num_threads
      hybrid%ithr         = ithr
      hybrid%localsense   = 0
      hybrid%masterthread = (par%masterproc .and. ithr==0)

  end function config_thread_region_par

  subroutine init_loop_ranges(nelemd)
      use iso_c_binding, only : c_int64_t, c_loc, c_ptr
      use cam_logfile, only : iulog

      integer, intent(in) :: nelemd
      logical, save :: proof_seen = .false.

      interface
         subroutine init_loop_ranges_codon(nelemd_c, nlev_c, qsize_c, horz_num_threads_c, &
              vert_num_threads_c, tracer_num_threads_c, work_pool_horz_p, work_pool_vert_p, &
              work_pool_trac_p) bind(c, name='init_loop_ranges_codon')
           import :: c_int64_t, c_ptr
           integer(c_int64_t), value :: nelemd_c, nlev_c, qsize_c
           integer(c_int64_t), value :: horz_num_threads_c, vert_num_threads_c, tracer_num_threads_c
           type(c_ptr), value :: work_pool_horz_p, work_pool_vert_p, work_pool_trac_p
         end subroutine init_loop_ranges_codon
      end interface

#define SE_MISC_TAG 22
#define SE_MISC_LABEL 'init_loop_ranges'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

      
      if ( init_ranges ) then
        nelemd_save=nelemd
!JMD#ifdef _OPENMP
        if ( .NOT. allocated(work_pool_horz) ) allocate(work_pool_horz(horz_num_threads,2))
        if ( .NOT. allocated(work_pool_vert) ) allocate(work_pool_vert(vert_num_threads,2))
        if ( .NOT. allocated(work_pool_trac) ) allocate(work_pool_trac(tracer_num_threads,2))
        call init_loop_ranges_codon(int(nelemd, c_int64_t), int(nlev, c_int64_t), &
             int(qsize, c_int64_t), int(horz_num_threads, c_int64_t), &
             int(vert_num_threads, c_int64_t), int(tracer_num_threads, c_int64_t), &
             c_loc(work_pool_horz(1,1)), c_loc(work_pool_vert(1,1)), c_loc(work_pool_trac(1,1)))
        if (.not. proof_seen) then
          write(iulog,*) 'init_loop_ranges implementation = codon'
          proof_seen = .true.
        endif

!JMD#endif
        init_ranges = .false.
      endif

  end subroutine init_loop_ranges

 subroutine set_region_num_threads( local_name )

  character(len=*), intent(in) :: local_name

  region_name = local_name

#ifdef _OPENMP

  if ( TRIM(region_name) == 'horizontal') then
    region_num_threads = horz_num_threads 
    call omp_set_num_threads(region_num_threads)
    return
  endif

  if ( TRIM(region_name) == 'vertical') then
    region_num_threads = vert_num_threads 
    call omp_set_num_threads(region_num_threads)
    return
  endif
  
  if ( TRIM(region_name) == 'tracer' ) then
    region_num_threads = tracer_num_threads
    call omp_set_num_threads(region_num_threads)
    return
  endif
    
  if ( TRIM(region_name) == 'vertical_and_tracer' ) then
    region_num_threads = vert_num_threads*tracer_num_threads
    call omp_set_num_threads(region_num_threads)
    return
  endif
 
#endif
    
  end subroutine set_region_num_threads

  subroutine set_loop_ranges (pybrid)

  type (hybrid_p) :: pybrid

  integer :: ibeg_range, iend_range
  integer :: kbeg_range, kend_range
  integer :: qbeg_range, qend_range
  integer :: idthread

#ifdef _OPENMP
  idthread = omp_get_thread_num()

  if ( TRIM(region_name) == 'horizontal' ) then
    call set_thread_ranges_1D ( work_pool_horz, ibeg_range, iend_range, idthread )
    pybrid%ibeg = ibeg_range; pybrid%iend = iend_range
    pybrid%kbeg = 1;          pybrid%kend = nlev
    pybrid%qbeg = 1;          pybrid%qend = qsize
  endif

  if ( TRIM(region_name) == 'vertical' ) then
    call set_thread_ranges_1D ( work_pool_vert, kbeg_range, kend_range, idthread )
    !FIXME: need to set ibeg, iend as well
    pybrid%kbeg = kbeg_range; pybrid%kend = kend_range
    pybrid%qbeg = 1;          pybrid%qend = qsize
  endif

  if ( TRIM(region_name) == 'tracer' ) then
    call set_thread_ranges_1D ( work_pool_trac, qbeg_range, qend_range, idthread )
    !FIXME: need to set ibeg, iend as well
    pybrid%kbeg = 1;          pybrid%kend = nlev
    pybrid%qbeg = qbeg_range; pybrid%qend = qend_range
  endif

  if ( TRIM(region_name) == 'vertical_and_tracer' ) then
    call set_thread_ranges_2D ( work_pool_vert, work_pool_trac, kbeg_range, kend_range, &
                                                       qbeg_range, qend_range, idthread )
    !FIXME: need to set ibeg, iend as well
    pybrid%kbeg = kbeg_range; pybrid%kend = kend_range
    pybrid%qbeg = qbeg_range; pybrid%qend = qend_range
  endif

#else
  call reset_loop_ranges(pybrid, region_name)
#endif 

  end subroutine set_loop_ranges

  subroutine get_loop_ranges (pybrid, ibeg, iend, kbeg, kend, qbeg, qend)
  use cam_logfile, only : iulog
  interface
    subroutine get_loop_ranges_codon(ibeg_in_c, iend_in_c, kbeg_in_c, kend_in_c, &
         qbeg_in_c, qend_in_c, mask_c, ibeg_p, iend_p, kbeg_p, kend_p, qbeg_p, qend_p) &
         bind(c, name='get_loop_ranges_codon')
      import :: c_int64_t, c_ptr
      integer(c_int64_t), value :: ibeg_in_c, iend_in_c, kbeg_in_c, kend_in_c
      integer(c_int64_t), value :: qbeg_in_c, qend_in_c, mask_c
      type(c_ptr), value :: ibeg_p, iend_p, kbeg_p, kend_p, qbeg_p, qend_p
    end subroutine get_loop_ranges_codon
  end interface

  type (hybrid_t), intent(in) :: pybrid
  integer, optional, intent(out), target :: ibeg, iend, kbeg, kend, qbeg, qend
  integer, target :: ibeg_dummy, iend_dummy, kbeg_dummy, kend_dummy, qbeg_dummy, qend_dummy
  type(c_ptr) :: ibeg_ptr, iend_ptr, kbeg_ptr, kend_ptr, qbeg_ptr, qend_ptr
  integer :: mask
  logical, save :: proof_seen = .false.

  mask = 0
  ibeg_dummy = 0; iend_dummy = 0; kbeg_dummy = 0
  kend_dummy = 0; qbeg_dummy = 0; qend_dummy = 0
  if (.not. proof_seen) then
    write(iulog,*) 'get_loop_ranges implementation = codon'
    proof_seen = .true.
  endif
  ibeg_ptr = c_loc(ibeg_dummy)
  iend_ptr = c_loc(iend_dummy)
  kbeg_ptr = c_loc(kbeg_dummy)
  kend_ptr = c_loc(kend_dummy)
  qbeg_ptr = c_loc(qbeg_dummy)
  qend_ptr = c_loc(qend_dummy)
  if ( present(ibeg) ) then
    mask = mask + 1
    ibeg_ptr = c_loc(ibeg)
  endif
  if ( present(iend) ) then
    mask = mask + 2
    iend_ptr = c_loc(iend)
  endif
  if ( present(kbeg) ) then
    mask = mask + 4
    kbeg_ptr = c_loc(kbeg)
  endif
  if ( present(kend) ) then
    mask = mask + 8
    kend_ptr = c_loc(kend)
  endif
  if ( present(qbeg) ) then
    mask = mask + 16
    qbeg_ptr = c_loc(qbeg)
  endif
  if ( present(qend) ) then
    mask = mask + 32
    qend_ptr = c_loc(qend)
  endif
  call get_loop_ranges_codon(int(pybrid%ibeg, c_int64_t), int(pybrid%iend, c_int64_t), &
       int(pybrid%kbeg, c_int64_t), int(pybrid%kend, c_int64_t), &
       int(pybrid%qbeg, c_int64_t), int(pybrid%qend, c_int64_t), int(mask, c_int64_t), &
       ibeg_ptr, iend_ptr, kbeg_ptr, kend_ptr, qbeg_ptr, qend_ptr)

  end subroutine get_loop_ranges

  function threadOwnsVertlevel(hybrid,value) result(found) 

   type (hybrid_t), intent(in) :: hybrid
   integer, intent(in) :: value
   logical :: found

   found = .false.
   if ((value >= hybrid%kbeg) .and. (value <= hybrid%kend)) then 
      found = .true. 
   endif
 
  end function threadOwnsVertlevel

  function threadOwnsTracer(hybrid,value) result(found) 

   type (hybrid_t), intent(in) :: hybrid
   integer, intent(in) :: value
   logical :: found

   found = .false.
   if ((value >= hybrid%qbeg) .and. (value <= hybrid%qend)) then 
      found = .true. 
   endif
 
  end function threadOwnsTracer

  subroutine reset_loop_ranges (pybrid, region_name)

  type (hybrid_p)              :: pybrid
  character(len=*), intent(in) :: region_name

  if ( TRIM(region_name) == 'vertical' ) then
    pybrid%kbeg = 1; pybrid%kend = nlev
  endif

  if ( TRIM(region_name) == 'tracer' ) then
    pybrid%qbeg = 1; pybrid%qend = qsize
  endif

  if ( TRIM(region_name) == 'vertical_and_tracer' ) then
    pybrid%kbeg = 1; pybrid%kend = nlev
    pybrid%qbeg = 1; pybrid%qend = qsize
  endif

  end subroutine reset_loop_ranges 

  subroutine set_thread_ranges_3D ( work_pool_x, work_pool_y, work_pool_z, &
                       beg_range_1, end_range_1, beg_range_2, end_range_2, &
                                    beg_range_3, end_range_3, idthread )

  integer, intent (in   ) :: work_pool_x(:,:)
  integer, intent (in   ) :: work_pool_y(:,:)
  integer, intent (in   ) :: work_pool_z(:,:)
  integer, intent (inout) :: beg_range_1 
  integer, intent (inout) :: end_range_1 
  integer, intent (inout) :: beg_range_2 
  integer, intent (inout) :: end_range_2 
  integer, intent (inout) :: beg_range_3 
  integer, intent (inout) :: end_range_3 
  integer, intent (inout) :: idthread

  integer :: index(3)
  integer :: i, j, k, ind, irange, jrange, krange

  ind = 0

  krange = SIZE(work_pool_z,1)
  jrange = SIZE(work_pool_y,1)
  irange = SIZE(work_pool_x,1)
  do k = 1, krange
    do j = 1, jrange
      do i = 1, irange
        if( ind == idthread ) then
          index(1) = i
          index(2) = j
          index(3) = k
        endif
        ind = ind + 1
      enddo
    enddo
  enddo
  beg_range_1 = work_pool_x(index(1),1)
  end_range_1 = work_pool_x(index(1),2)
  beg_range_2 = work_pool_y(index(2),1)
  end_range_2 = work_pool_y(index(2),2)
  beg_range_3 = work_pool_z(index(3),1)
  end_range_3 = work_pool_z(index(3),2)

!  write(6,1000) idthread, beg_range_1, end_range_1, &
!                          beg_range_2, end_range_2, &
!                          beg_range_3, end_range_3
!  call flush(6)
1000 format( 'set_thread_ranges_3D', 7(i4) )

  end subroutine set_thread_ranges_3D

  subroutine set_thread_ranges_2D( work_pool_x, work_pool_y, beg_range_1, end_range_1, &
                                                    beg_range_2, end_range_2, idthread )

  integer, intent (in   ) :: work_pool_x(:,:)
  integer, intent (in   ) :: work_pool_y(:,:)
  integer, intent (inout) :: beg_range_1 
  integer, intent (inout) :: end_range_1 
  integer, intent (inout) :: beg_range_2 
  integer, intent (inout) :: end_range_2 
  integer, intent (inout) :: idthread

  integer :: index(2)
  integer :: i, j, ind, irange, jrange

  ind = 0

  jrange = SIZE(work_pool_y,1)
  irange = SIZE(work_pool_x,1)
  do j = 1, jrange
    do i = 1, irange
      if( ind == idthread ) then
        index(1) = i
        index(2) = j
      endif
      ind = ind + 1
    enddo
  enddo
  beg_range_1 = work_pool_x(index(1),1)
  end_range_1 = work_pool_x(index(1),2)
  beg_range_2 = work_pool_y(index(2),1)
  end_range_2 = work_pool_y(index(2),2)

! write(6,1000) idthread, beg_range_1, end_range_1, &
!                         beg_range_2, end_range_2
! call flush(6)

1000 format( 'set_thread_ranges_2D', 7(i4) )

  end subroutine set_thread_ranges_2D

  subroutine set_thread_ranges_1D( work_pool, beg_range, end_range, idthread )
  use cam_logfile, only : iulog

  interface
    subroutine set_thread_ranges_1d_codon(work_pool_p, nrows_c, idthread_c, &
         beg_range_p, end_range_p) bind(c, name='set_thread_ranges_1d_codon')
      import :: c_int64_t, c_ptr
      type(c_ptr), value :: work_pool_p
      integer(c_int64_t), value :: nrows_c, idthread_c
      type(c_ptr), value :: beg_range_p, end_range_p
    end subroutine set_thread_ranges_1d_codon
  end interface

  integer, intent (in), target :: work_pool(:,:)
  integer, intent (inout), target :: beg_range
  integer, intent (inout), target :: end_range
  integer, intent (inout) :: idthread
  logical, save :: proof_seen = .false.

  call set_thread_ranges_1d_codon(c_loc(work_pool(1,1)), int(size(work_pool, 1), c_int64_t), &
       int(idthread, c_int64_t), c_loc(beg_range), c_loc(end_range))
  if (.not. proof_seen) then
    write(iulog,*) 'set_thread_ranges_1d implementation = codon'
    proof_seen = .true.
  endif

! write(6,1000) idthread, beg_range, end_range
! call flush(6)
1000 format( 'set_thread_ranges_1D', 7(i4) )

  end subroutine set_thread_ranges_1D

  subroutine create_work_pool( start_domain, end_domain, ndomains, ipe, beg_index, end_index )
  use cam_logfile, only : iulog

  interface
    subroutine create_work_pool_codon(start_domain_c, end_domain_c, ndomains_c, ipe_c, &
         beg_index_p, end_index_p) bind(c, name='create_work_pool_codon')
      import :: c_int64_t, c_ptr
      integer(c_int64_t), value :: start_domain_c, end_domain_c, ndomains_c, ipe_c
      type(c_ptr), value :: beg_index_p, end_index_p
    end subroutine create_work_pool_codon
  end interface

  integer, intent(in) :: start_domain, end_domain
  integer, intent(in) :: ndomains, ipe
  integer, intent(out), target :: beg_index, end_index
  logical, save :: proof_seen = .false.

  call create_work_pool_codon(int(start_domain, c_int64_t), int(end_domain, c_int64_t), &
       int(ndomains, c_int64_t), int(ipe, c_int64_t), c_loc(beg_index), c_loc(end_index))
  if (.not. proof_seen) then
    write(iulog,*) 'create_work_pool implementation = codon'
    proof_seen = .true.
  endif

  end subroutine create_work_pool

  subroutine get_number_threads(maxthreads)
  use thread_mod, only : omp_get_nested

    integer, INTENT(OUT) :: maxthreads

    character(len=32) :: instring, string1, string2
    integer :: thr1=1, thr2=1

    call GET_ENVIRONMENT_VARIABLE('OMP_NUM_THREADS', instring)
    ! omp_get_nested is deprecated
!    if ( omp_get_nested() ) then
!      call split_string(instring, string1, string2, ',')
!      read (string1,'(I4)') thr1; read (string2,'(I4)') thr2
!    else
      read (instring,'(I4)') thr1
!    endif
    maxthreads = thr1 * thr2

  end subroutine get_number_threads

  subroutine split_string(instring, string1, string2, delim)
    character(*), INTENT(INOUT) :: instring
    character(*), INTENT(IN)    :: delim
    character(*), INTENT(OUT)   :: string1, string2

    integer :: index

    instring = TRIM(instring)

    index = SCAN(instring,delim)
    string1 = instring(1:index-1)
    string2 = instring(index+1:)

  end subroutine split_string

end module hybrid_mod
