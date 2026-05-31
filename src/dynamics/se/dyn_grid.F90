module dyn_grid
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: Definition of dynamics computational grid.
  !
  ! Method: Variables are private; interface routines used to extract
  !         information for use in user code. Global column index range
  !         defined using full (unreduced) grid. 
  ! 
  ! Entry points:
  !      get_block_bounds_d       get first and last indices in global 
  !                               block ordering
  !      get_block_gcol_d         get column indices for given block
  !      get_block_gcol_cnt_d     get number of columns in given block
  !      get_block_lvl_cnt_d      get number of vertical levels in column
  !      get_block_levels_d       get vertical levels in column
  !      get_gcol_block_d         get global block indices and local columns 
  !                               index for given global column index
  !      get_gcol_block_cnt_d     get number of blocks containing data
  !                               from a given global column index
  !      get_block_owner_d        get process "owning" given block
  !      get_horiz_grid_d         get horizontal grid coordinates
  !      get_horiz_grid_dim_d     get horizontal dimensions of dynamics grid
  !      dyn_grid_get_pref        get reference pressures for the dynamics grid
  !      dyn_grid_get_elem_coords get coordinates of a specified block element 
  !                               of the dynamics grid
  !      dyn_grid_find_gcols      finds nearest column for given lat/lon
  !      dyn_grid_get_colndx      get element block/column and MPI process indices 
  !                               corresponding to a specified global column index
  !
  ! Author: Jim Edwards and Patrick Worley
  ! 
  !-----------------------------------------------------------------------
  use element_mod, only : element_t
  use cam_logfile, only : iulog
  use shr_kind_mod, only: r8 => shr_kind_r8
  use fvm_control_volume_mod, only : fvm_struct
  use iso_c_binding, only : c_int64_t

  implicit none
  private
  save

  integer           :: ngcols_d = 0     ! number of dynamics columns
  integer, public, parameter :: ptimelevels=2
  logical :: gblocks_need_initialized=.true.


  type block_global_data
     integer :: UniquePtOffset
     integer :: NumUniqueP
     integer :: LocalID
     integer :: Owner
  end type block_global_data


  type(block_global_data), allocatable :: gblocks(:)
  type(element_t), public, pointer :: elem(:) => null()
  type(fvm_struct), public, pointer :: fvm(:) => null()

  public :: dyn_grid_init, get_block_owner_d, get_gcol_block_d, get_gcol_block_cnt_d, &
       get_block_gcol_cnt_d, get_horiz_grid_dim_d, get_block_levels_d, &
       get_block_gcol_d, get_block_bounds_d, get_horiz_grid_d, &
       get_block_lvl_cnt_d, set_horiz_grid_cnt_d, get_dyn_grid_parm, &
       get_dyn_grid_parm_real2d, get_dyn_grid_parm_real1d, get_block_ldof_d, &
       dyn_grid_get_pref

  public :: dyn_grid_get_elem_coords
  public :: dyn_grid_find_gcols
  public :: dyn_grid_get_colndx

  integer, pointer :: fdofP_local(:,:,:) => null()

  real(r8), public, pointer :: w(:) => null()        ! weights
  real(r8), pointer :: clat(:) => null()     ! model latitudes (radians)
  real(r8), pointer :: clon(:,:) => null()   ! model longitudes (radians)
  real(r8), pointer :: latdeg(:) => null()   ! model latitudes (degrees)
  real(r8), pointer :: londeg(:,:) => null() ! model longitudes (degrees)


!========================================================================
contains
!========================================================================

subroutine dyn_grid_init()
   use iso_c_binding, only : c_int64_t

#define SE_MISC_TAG 7
#define SE_MISC_LABEL 'dyn_grid_init'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

   nullify(clat, clon)

end subroutine dyn_grid_init

  !========================================================================
  !
  Subroutine get_block_bounds_d(block_first,block_last)
    use dimensions_mod, only: nelem
    !----------------------------------------------------------------------- 
    ! 
    !                          
    ! Purpose: Return first and last indices used in global block ordering
    ! 
    ! Method: 
    ! 
    ! Author: Jim Edwards
    ! 
    !-----------------------------------------------------------------------
    !------------------------------Arguments--------------------------------
    integer, intent(out) :: block_first  ! first (global) index used for blocks
    integer, intent(out) :: block_last   ! last (global) index used for blocks

    !-----------------------------------------------------------------------
    block_first = 1
    block_last = nelem

    return
  end subroutine get_block_bounds_d
  !
  !========================================================================
  !
  subroutine get_block_gcol_d(blockid,size,cdex)
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use cam_logfile, only : iulog
    !----------------------------------------------------------------------- 
    ! 
    !                          
    ! Purpose: Return list of dynamics column indices in given block
    ! 
    ! Method: 
    ! 
    ! Author: Jim Edwards
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !------------------------------Arguments--------------------------------
    integer, intent(in) :: blockid      ! global block id
    integer, intent(in) :: size         ! array size

    integer, intent(out), target :: cdex(size)   ! global column indices
    !
    logical, save :: proof_seen = .false.
    interface
       subroutine get_block_gcol_d_codon(size_c, unique_pt_offset_c, cdex_p) &
            bind(c, name='get_block_gcol_d_codon')
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: size_c, unique_pt_offset_c
         type(c_ptr), value :: cdex_p
       end subroutine get_block_gcol_d_codon
    end interface

    if(gblocks_need_initialized) call gblocks_init()
    call get_block_gcol_d_codon(int(size, c_int64_t), &
         int(gblocks(blockid)%UniquePtOffset, c_int64_t), c_loc(cdex))
    if (.not. proof_seen) then
       write(iulog,*) 'get_block_gcol_d implementation = codon'
       proof_seen = .true.
    endif

    return
  end subroutine get_block_gcol_d

  subroutine get_block_ldof_d(nlev, ldof)
    use dof_mod, only : CreateMetaData
    use dimensions_mod, only : npsq, np, nelemd
    use parallel_mod,   only : par
    use element_mod, only : index_t
    use pio, only : pio_offset_kind

    integer, intent(in) :: nlev
    integer(kind=pio_offset_kind), pointer    :: ldof(:)

    integer :: i, j, k, ie, ii
    type (index_t), pointer  :: idx     
    integer :: mydof(npsq)

    logical, save :: firstcall = .true.


    if(firstcall) then
       allocate(fdofp_local(np,np,nelemd))
       fdofp_local=0
       do ie=1,nelemd
          idx => elem(ie)%idxP
          do ii=1,idx%NumUniquePts
             i=idx%ia(ii)
             j=idx%ja(ii)          
             fdofp_local(i,j,ie) = (idx%UniquePtoffset+ii-1)
          end do
       end do
       firstcall=.false.
!    print *,__FILE__,__LINE__,maxval(fdofp_local),count(fdofp_local/=0)
    end if

    allocate(ldof(nelemd*npsq*nlev))       
    ldof = 0
    do ie=1,nelemd
       do j=1,np
          do i=1,np
             if(fdofp_local(i,j,ie)>0) then
                do k=0,nlev-1
                   ldof(i+(j-1)*np+k*npsq+(ie-1)*npsq*nlev ) = fdofP_local(i,j,ie) + ngcols_d*k
                end do
             end if
          end do
       end do
    end do
!    print *,__FILE__,__LINE__,maxval(ldof),count(ldof/=0)

  end subroutine get_block_ldof_d




  !
  !========================================================================
  !
  integer function get_block_gcol_cnt_d(blockid)
    !----------------------------------------------------------------------- 
    ! 
    !                          
    ! Purpose: Return number of dynamics columns in indicated block
    ! 
    ! Method: 
    ! 
    ! Author: Jim Edwards
    ! 
    !-----------------------------------------------------------------------
    integer, intent(in) :: blockid
    integer :: ie

    if(gblocks_need_initialized) call gblocks_init()
    get_block_gcol_cnt_d=gblocks(blockid)%NumUniqueP

    return
  end function get_block_gcol_cnt_d

  !
  !========================================================================
  !
  integer function get_block_lvl_cnt_d(blockid,bcid)
    use iso_c_binding, only : c_int64_t
    use pmgrid, only: plevp

    !-----------------------------------------------------------------------
    !
    !
    ! Purpose: Return number of levels in indicated column. If column
    !          includes surface fields, then it is defined to also
    !          include level 0.
    !
    ! Method:
    !
    ! Author: Patrick Worley
    !
    !-----------------------------------------------------------------------

    implicit none
    interface
       function get_block_lvl_cnt_d_codon(plevp_c) result(level_count_c) &
            bind(c, name='get_block_lvl_cnt_d_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: plevp_c
         integer(c_int64_t) :: level_count_c
       end function get_block_lvl_cnt_d_codon
    end interface
    !------------------------------Arguments--------------------------------
    integer, intent(in) :: blockid  ! global block id
    integer, intent(in) :: bcid    ! column index within block
    character(len=32) :: impl_name
    integer :: impl_n, impl_status
    integer(c_int64_t) :: level_count_c
    logical, save :: proof_seen = .false.

    !-----------------------------------------------------------------------

    impl_name = 'codon'
    call get_environment_variable('DYN_GRID_IMPL', value=impl_name, &
         length=impl_n, status=impl_status)
    if (impl_status == 0 .and. impl_n > 0 .and. &
         trim(adjustl(impl_name(:impl_n))) == 'native') then
       get_block_lvl_cnt_d = plevp
       return
    endif

    level_count_c = get_block_lvl_cnt_d_codon(int(plevp, c_int64_t))
    get_block_lvl_cnt_d = int(level_count_c)
    if (.not. proof_seen) then
       write(iulog,*) 'get_block_lvl_cnt_d implementation = codon'
       proof_seen = .true.
    endif

    return
  end function get_block_lvl_cnt_d
  !
  !========================================================================
  !
  subroutine get_block_levels_d(blockid, bcid, lvlsiz, levels)
    use pmgrid,         only: plev
    use cam_abortutils, only: endrun

    !-----------------------------------------------------------------------
    !
    !
    ! Purpose: Return level indices in indicated column. If column
    !          includes surface fields, then it is defined to also
    !          include level 0.
    !
    ! Method:
    !
    ! Author: Patrick Worley
    !
    !-----------------------------------------------------------------------

    implicit none
    !------------------------------Arguments--------------------------------
    integer, intent(in) :: blockid  ! global block id
    integer, intent(in) :: bcid    ! column index within block
    integer, intent(in) :: lvlsiz   ! dimension of levels array

    integer, intent(out) :: levels(lvlsiz) ! levels indices for block

    !---------------------------Local workspace-----------------------------
    !

    !---------------------------Local workspace-----------------------------
    !
    integer k                      ! loop index
    !-----------------------------------------------------------------------
    !  latitude slice block
    if (lvlsiz < plev + 1) then
       write(iulog,*)'GET_BLOCK_LEVELS_D: levels array not large enough (', &
            lvlsiz,' < ',plev + 1,' ) '
       call endrun
    else
       do k=0,plev
          levels(k+1) = k
       enddo
       do k=plev+2,lvlsiz
          levels(k) = -1
       enddo
    endif

    return
  end subroutine get_block_levels_d

  !
  !========================================================================
  !
  subroutine get_gcol_block_d(gcol,cnt,blockid,bcid,localblockid)
    use dimensions_mod, only: nelemd, nelem
    use kinds,          only: int_kind
    use cam_abortutils, only: endrun
    !----------------------------------------------------------------------- 
    ! 
    !                          
    ! Purpose: Return global block index and local column index
    !          for global column index
    !           element array is naturally ordered
    ! 
    ! Method: 
    ! 
    ! Author: Jim Edwards
    ! Modified: replaced linear search with best guess and binary search   
    !           (Pat Worley, 7/2/09)
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !------------------------------Arguments--------------------------------
    integer, intent(in) :: gcol     ! global column index
    integer, intent(in) :: cnt      ! size of blockid and bcid arrays

    integer, intent(out) :: blockid(cnt) ! block index
    integer, intent(out) :: bcid(cnt)    ! column index within block
    integer, intent(out), optional :: localblockid(cnt)
    integer :: sb, eb, ie, high, low
    logical :: found
    integer, save :: iedex_save = 1
    !---------------------------Local workspace-----------------------------
    !
    if(gblocks_need_initialized) call gblocks_init()
    found = .false.
    low = 1
    high = nelem

    ! check whether previous found element is the same here
    if (.not. found) then
       ie = iedex_save
       sb = gblocks(ie)%UniquePtOffset
       if (gcol >= sb) then
          eb = sb+gblocks(ie)%NumUniqueP
          if (gcol < eb) then
             found = .true.
          else
             low = ie
          endif
       else
          high = ie
       endif
    endif

    ! check whether next element  is the one wanted
    if ((.not. found) .and. &
        ((low .eq. iedex_save) .or. (iedex_save .eq. nelem))) then
       ie = iedex_save + 1
       if (ie > nelem) ie = 1

       sb = gblocks(ie)%UniquePtOffset
       if (gcol >= sb) then
          eb = sb+gblocks(ie)%NumUniqueP
          if (gcol < eb) then
             found = .true.
          else
             low = ie
          endif
       else
          high = ie
       endif
    endif

    ! otherwise, use a binary search to find element
    if (.not. found) then
       ! (start with a sanity check)
       ie = low
       sb = gblocks(ie)%UniquePtOffset

       ie = high
       eb = gblocks(ie)%UniquePtOffset + gblocks(ie)%NumUniqueP

       if ((gcol < sb) .or.  (gcol >= eb)) then
          do ie=1,nelemd
             write(iulog,*) __LINE__,ie,elem(ie)%idxP%UniquePtOffset,elem(ie)%idxP%NumUniquePts
          end do
          call endrun()
       end if

       do while (.not. found)
          ie = low + (high-low)/2;

          sb = gblocks(ie)%UniquePtOffset
          if (gcol >= sb) then
             eb = sb+gblocks(ie)%NumUniqueP
             if (gcol < eb) then
                found = .true.
             else
                low = ie+1
             endif
          else
             high = ie-1
          endif
       enddo
    endif

    blockid(1)=ie
    if(present(localblockid)) localblockid(1)=gblocks(ie)%LocalID
    bcid(1)=gcol-sb+1
    iedex_save = ie

 return
end subroutine get_gcol_block_d

!
!========================================================================
!
integer function get_gcol_block_cnt_d(gcol)
 use iso_c_binding, only : c_int64_t

  !----------------------------------------------------------------------- 
  ! 
  !                          
  ! Purpose: Return number of blocks contain data for the vertical column
  !          with the given global column index
  ! 
  ! Method: 
  ! 
  ! Author: Patrick Worley
  ! 
  !-----------------------------------------------------------------------

 implicit none
 interface
    function get_gcol_block_cnt_d_codon() result(block_count_c) &
         bind(c, name='get_gcol_block_cnt_d_codon')
      import :: c_int64_t
      integer(c_int64_t) :: block_count_c
    end function get_gcol_block_cnt_d_codon
 end interface
 !------------------------------Arguments--------------------------------
 integer, intent(in) :: gcol     ! global column index
 character(len=32) :: impl_name
 integer :: impl_n, impl_status
 integer(c_int64_t) :: block_count_c
 logical, save :: proof_seen = .false.
 !-----------------------------------------------------------------------
 impl_name = 'codon'
 call get_environment_variable('DYN_GRID_IMPL', value=impl_name, &
      length=impl_n, status=impl_status)
 if (impl_status == 0 .and. impl_n > 0 .and. &
      trim(adjustl(impl_name(:impl_n))) == 'native') then
    get_gcol_block_cnt_d = 1
    return
 endif

 block_count_c = get_gcol_block_cnt_d_codon()
 get_gcol_block_cnt_d = int(block_count_c)
 if (.not. proof_seen) then
    write(iulog,*) 'get_gcol_block_cnt_d implementation = codon'
    proof_seen = .true.
 endif

 return
end function get_gcol_block_cnt_d

!
!========================================================================
!
integer function get_block_owner_d(blockid)
 use cam_abortutils, only: endrun
 use iso_c_binding, only : c_int64_t
 use cam_logfile, only : iulog
 !----------------------------------------------------------------------- 
 ! 
 !                          
 ! Purpose: Return id of processor that "owns" the indicated block
 ! 
 ! Method: 
 ! 
 ! Author: Jim Edwards
 ! 
 !-----------------------------------------------------------------------

 implicit none
 !------------------------------Arguments--------------------------------
 integer, intent(in) :: blockid  ! global block id
 logical, save :: proof_seen = .false.
 interface
    function get_block_owner_d_codon(owner_c) result(owner_out) bind(c, name='get_block_owner_d_codon')
      use iso_c_binding, only : c_int64_t
      integer(c_int64_t), value :: owner_c
      integer(c_int64_t) :: owner_out
    end function get_block_owner_d_codon
 end interface

 !-----------------------------------------------------------------------
 if(gblocks_need_initialized) call gblocks_init()

 get_block_owner_d = int(get_block_owner_d_codon(int(gblocks(blockid)%Owner, c_int64_t)))
 if(get_block_owner_d>-1) then
    if (.not. proof_seen) then
       write(iulog,*) 'get_block_owner_d implementation = codon'
       proof_seen = .true.
    endif
 else
    call endrun('Block owner not assigned in gblocks_init')
 end if

 return
end function get_block_owner_d

!
!========================================================================
!
subroutine get_horiz_grid_dim_d(hdim1_d,hdim2_d)

  !----------------------------------------------------------------------- 
  ! 
  !                          
  ! Purpose: Returns declared horizontal dimensions of computational grid.
  !          For non-lon/lat grids, declare grid to be one-dimensional,
  !          i.e., (ngcols_d x 1)
  ! 
  ! Method: 
  ! 
  ! Author: Patrick Worley
  ! 
  !-----------------------------------------------------------------------

  !------------------------------Arguments--------------------------------
 integer, intent(out) :: hdim1_d           ! first horizontal dimension
 integer, intent(out), optional :: hdim2_d           ! second horizontal dimension
 !-----------------------------------------------------------------------
 hdim1_d = ngcols_d
 if(present(hdim2_d)) hdim2_d = 1

 return
end subroutine get_horiz_grid_dim_d
!
!========================================================================
!
subroutine set_horiz_grid_cnt_d(NumUniqueCols)
 integer, intent(in) :: NumUniqueCols
 !----------------------------------------------------------------------- 
 ! 
 !                          
 ! Purpose: Set number of columns in the dynamics computational grid
 ! 
 ! Method: 
 ! 
 ! Author: Jim Edwards
 ! 
 !-----------------------------------------------------------------------
 ngcols_d = NumUniqueCols


 return
end subroutine set_horiz_grid_cnt_d
!
!========================================================================
!
subroutine get_horiz_grid_d(nxy,clat_d_out,clon_d_out,area_d_out, &
    wght_d_out)

 use shr_kind_mod,   only: r8 => shr_kind_r8
 use physconst,      only: pi
 use dof_mod,        only: UniqueCoords, UniquePoints
 use dimensions_mod, only: nelem, nelemd, nelemdmax, np
 use spmd_utils,     only: iam, npes, mpi_integer, mpir8, mpicom

 use cam_abortutils, only: endrun
 !----------------------------------------------------------------------- 
 ! 
 !                          
 ! Purpose: Return latitude and longitude (in radians), column surface
 !          area (in radians squared) and surface integration weights
 !          for global column indices that will be passed to/from physics
 ! 
 ! Method: 
 ! 
 ! Author: Jim Edwards
 ! 
 !-----------------------------------------------------------------------
 integer, intent(in)   :: nxy                     ! array sizes

 real(r8), intent(out), optional :: clat_d_out(:) ! column latitudes
 real(r8), intent(out), optional :: clon_d_out(:) ! column longitudes
 real(r8), intent(out), target, optional :: area_d_out(:) ! column surface 
 !  area
 real(r8), intent(out), target, optional :: wght_d_out(:) ! column integration
 !  weight

 integer :: ie, i, sb, eb, p, j

 real (r8), dimension(np,np) :: areaw
 real(r8), pointer :: area_d(:)
 real(r8) :: scale_factor
 integer :: ierr
 integer :: rdispls(npes), recvcounts(npes), gid(npes)
 integer :: ibuf
 real(r8) :: rbuf(nxy)
 logical :: compute_area 

 if (nxy < ngcols_d) then
    write(iulog,*)'GET_HORIZ_GRID_D: arrays not large enough (', &
         nxy,' < ',ngcols_d,' ) '
    call endrun
 endif

 compute_area = .false.
 if ( present(area_d_out) ) then
    if (size(area_d_out) .ne. nxy) then
       call endrun('bad area_d_out array size in dyn_grid')
    end if
    area_d => area_d_out
    compute_area = .true.
 else if ( present(wght_d_out) ) then
    if (size(wght_d_out) .ne. nxy) then
       call endrun('bad wght_d_out array size in dyn_grid')
    end if
    area_d => wght_d_out
    compute_area = .true.
 endif

 if(gblocks_need_initialized) call gblocks_init()
 ! could improve this
 if(.not. associated(clat)) then	
    allocate(clat(nxy), clon(nxy,1))
 end if

 clon(:,1) = -iam 
 do ie=1,nelemdmax        
    if(ie<=nelemd) then
       rdispls(iam+1)=elem(ie)%idxp%UniquePtOffset-1
       eb = rdispls(iam+1) + elem(ie)%idxp%NumUniquePts
       recvcounts(iam+1)=elem(ie)%idxP%NumUniquePts
       call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep, &
            clat(rdispls(iam+1)+1:eb), &
            clon(rdispls(iam+1)+1:eb,1)) 

       if (compute_area) then
          areaw = 1.0_r8/elem(ie)%rspheremp(:,:)         
          call UniquePoints(elem(ie)%idxP, areaw, area_d(rdispls(iam+1)+1:eb))
       endif
    else     
       rdispls(iam+1) = 0
       recvcounts(iam+1) = 0
    endif
    ibuf=rdispls(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, rdispls, &
         1, mpi_integer, mpicom, ierr)

    ibuf=recvcounts(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, recvcounts, &
         1, mpi_integer, mpicom, ierr)

    sb=rdispls(iam+1)+1
    eb =rdispls(iam+1) + recvcounts(iam+1)

    rbuf(1:recvcounts(iam+1)) = clat(sb:eb)  ! whats going to happen if end=0?
    call mpi_allgatherv(rbuf,recvcounts(iam+1),mpir8, clat, &
         recvcounts(:), rdispls(:), mpir8, mpicom, ierr)

    rbuf(1:recvcounts(iam+1)) = clon(sb:eb,1)
    call mpi_allgatherv(rbuf,recvcounts(iam+1),mpir8, clon, &
         recvcounts(:), rdispls(:), mpir8, mpicom, ierr)

    if (compute_area) then
       rbuf(1:recvcounts(iam+1)) = area_d(sb:eb)
       call mpi_allgatherv(rbuf,recvcounts(iam+1),mpir8, area_d, &
         recvcounts(:), rdispls(:), mpir8, mpicom, ierr)
    endif
 end do
 call mpi_barrier(mpicom,ierr)

 if ( present(clat_d_out) ) then
    if (size(clat_d_out) .ne. nxy) then
       call endrun('bad clat_d_out array size in dyn_grid')
    end if
    clat_d_out(:) = clat(:)
 endif

 if ( present(clon_d_out) ) then
    if (size(clon_d_out) .ne. nxy) then
       call endrun('bad clon_d_out array size in dyn_grid')
    end if
    clon_d_out(:) = clon(:,1)
 endif

 !if one of area_d_out  or wght_d_out was present, than it was computed
 !above.  if they were *both* present, then do this: 
 if ( present(area_d_out) .and. present(wght_d_out) ) then
    wght_d_out(:) = area_d_out(:)
 endif


 return
end subroutine get_horiz_grid_d


subroutine get_gcol_lat(gcid, lat)
  use shr_kind_mod, only: r8 => shr_kind_r8

  integer, intent(in) :: gcid(:)
  real(r8), intent(out) :: lat(:)
  integer :: k, glen

  glen = size(gcid)

  do k=1,glen
     lat(k) = clat(gcid(k))
  end do

end subroutine get_gcol_lat

subroutine get_gcol_lon(gcid, lon)

  integer, intent(in) :: gcid(:)
  real(r8), intent(out) :: lon(:)
  integer :: k, glen

  glen = size(gcid)

  do k=1,glen
     lon(k) = clon(gcid(k),1)
  end do

end subroutine get_gcol_lon



subroutine gblocks_init()
 use dimensions_mod, only: nelem, nelemd, nelemdmax
 use spmd_utils,     only: iam, npes, mpi_integer, mpicom
 use cam_abortutils, only: endrun
 integer :: ie, p
 integer :: ibuf
 integer :: ierr
 integer :: rdispls(npes), recvcounts(npes), gid(npes), lid(npes)


 if(.not.allocated(gblocks)) then
    allocate(gblocks(nelem))
    do ie=1,nelem
       gblocks(ie)%Owner=-1
       gblocks(ie)%UniquePtOffset=-1
       gblocks(ie)%NumUniqueP    =-1
       gblocks(ie)%LocalID       =-1
    end do
 end if

 do ie=1,nelemdmax        
    if(ie<=nelemd) then
       rdispls(iam+1)=elem(ie)%idxP%UniquePtOffset-1
       gid(iam+1)=elem(ie)%GlobalID
       lid(iam+1)=ie
       recvcounts(iam+1)=elem(ie)%idxP%NumUniquePts
    else     
       rdispls(iam+1) = 0
       recvcounts(iam+1) = 0
       gid(iam+1)=0
    endif
    ibuf=lid(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, lid, &
         1, mpi_integer, mpicom, ierr)

    ibuf=gid(iam+1)
    call mpi_allgather(ibuf,1,mpi_integer, gid, &
         1, mpi_integer, mpicom, ierr)

    ibuf=rdispls(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, rdispls, &
         1, mpi_integer, mpicom, ierr)

    ibuf=recvcounts(iam+1)
    call mpi_allgather(ibuf, 1, mpi_integer, recvcounts, &
         1, mpi_integer, mpicom, ierr)
    do p=1,npes
       if(gid(p)>0) then
          gblocks(gid(p))%UniquePtOffset=rdispls(p)+1
          gblocks(gid(p))%NumUniqueP=recvcounts(p)
          gblocks(gid(p))%LocalID=lid(p)
          gblocks(gid(p))%Owner=p-1
       end if
    end do
 end do
 gblocks_need_initialized=.false.

end subroutine gblocks_init



!#######################################################################
   function get_dyn_grid_parm_real2d(name) result(rval)
     use iso_c_binding, only : c_int64_t
     use cam_logfile, only : iulog

     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:,:)
     integer(c_int64_t) :: selector
     logical, save :: proof_seen = .false.
     interface
        function get_dyn_grid_parm_real2d_codon(name_c) result(selector_c) &
             bind(c, name='get_dyn_grid_parm_real2d_codon')
          use iso_c_binding, only : c_int64_t
          integer(c_int64_t), value :: name_c
          integer(c_int64_t) :: selector_c
        end function get_dyn_grid_parm_real2d_codon
     end interface

     selector = get_dyn_grid_parm_real2d_codon(dyn_grid_parm_name_code(name))
     if(selector == 1_c_int64_t) then
        rval => clon
     else if(selector == 2_c_int64_t) then
        rval => londeg
     else
        nullify(rval)
     end if
     if (.not. proof_seen) then
        write(iulog,*) 'get_dyn_grid_parm_real2d implementation = codon'
        proof_seen = .true.
     endif
   end function get_dyn_grid_parm_real2d

!#######################################################################
   function get_dyn_grid_parm_real1d(name) result(rval)
     use iso_c_binding, only : c_int64_t
     use cam_logfile, only : iulog

     character(len=*), intent(in) :: name
     real(r8), pointer :: rval(:)
     integer(c_int64_t) :: selector
     logical, save :: proof_seen = .false.
     interface
        function get_dyn_grid_parm_real1d_codon(name_c) result(selector_c) &
             bind(c, name='get_dyn_grid_parm_real1d_codon')
          use iso_c_binding, only : c_int64_t
          integer(c_int64_t), value :: name_c
          integer(c_int64_t) :: selector_c
        end function get_dyn_grid_parm_real1d_codon
     end interface

     selector = get_dyn_grid_parm_real1d_codon(dyn_grid_parm_name_code(name))
     if(selector == 1_c_int64_t) then
        rval => clat
     else if(selector == 2_c_int64_t) then
        rval => latdeg
     else if(selector == 3_c_int64_t) then
        rval => w
     else
        nullify(rval)
     end if
     if (.not. proof_seen) then
        write(iulog,*) 'get_dyn_grid_parm_real1d implementation = codon'
        proof_seen = .true.
     endif
   end function get_dyn_grid_parm_real1d

   integer(c_int64_t) function dyn_grid_parm_name_code(name) result(code)
     use iso_c_binding, only : c_int64_t
     character(len=*), intent(in) :: name

     select case (trim(name))
     case ('clon')
        code = 1_c_int64_t
     case ('londeg')
        code = 2_c_int64_t
     case ('clat')
        code = 3_c_int64_t
     case ('latdeg')
        code = 4_c_int64_t
     case ('w')
        code = 5_c_int64_t
     case ('ne')
        code = 6_c_int64_t
     case ('np')
        code = 7_c_int64_t
     case ('npsq')
        code = 8_c_int64_t
     case ('nelemd')
        code = 9_c_int64_t
     case ('beglat')
        code = 10_c_int64_t
     case ('endlat')
        code = 11_c_int64_t
     case ('beglonxy')
        code = 12_c_int64_t
     case ('endlonxy')
        code = 13_c_int64_t
     case ('beglatxy')
        code = 14_c_int64_t
     case ('endlatxy')
        code = 15_c_int64_t
     case ('plat')
        code = 16_c_int64_t
     case ('plon')
        code = 17_c_int64_t
     case ('plev')
        code = 18_c_int64_t
     case ('plevp')
        code = 19_c_int64_t
     case ('nlon')
        code = 20_c_int64_t
     case ('nlat')
        code = 21_c_int64_t
     case default
        code = 0_c_int64_t
     end select
   end function dyn_grid_parm_name_code




integer function get_dyn_grid_parm(name) result(ival)
 use iso_c_binding, only : c_int64_t
 use cam_logfile, only : iulog
 use pmgrid, only : beglat, endlat, plat, plon, plev, plevp
 use dimensions_mod, only: ne, np, nelemd, npsq
 use interpolate_mod, only : get_interp_parameter
 character(len=*), intent(in) :: name
 integer(c_int64_t) :: name_code
 integer(c_int64_t) :: nlon_c, nlat_c
 integer(c_int64_t) :: ival_c
 logical, save :: proof_seen = .false.
 interface
    function get_dyn_grid_parm_codon(name_c, ne_c, np_c, npsq_c, nelemd_c, &
         beglat_c, endlat_c, ngcols_d_c, plat_c, plev_c, plevp_c, nlon_c, &
         nlat_c) result(ival_c) bind(c, name='get_dyn_grid_parm_codon')
      import :: c_int64_t
      integer(c_int64_t), value :: name_c, ne_c, np_c, npsq_c, nelemd_c
      integer(c_int64_t), value :: beglat_c, endlat_c, ngcols_d_c, plat_c
      integer(c_int64_t), value :: plev_c, plevp_c, nlon_c, nlat_c
      integer(c_int64_t) :: ival_c
    end function get_dyn_grid_parm_codon
 end interface

 name_code = dyn_grid_parm_name_code(name)
 nlon_c = -1_c_int64_t
 nlat_c = -1_c_int64_t
 if (name_code == 20_c_int64_t) then
    nlon_c = int(get_interp_parameter('nlon'), c_int64_t)
 else if (name_code == 21_c_int64_t) then
    nlat_c = int(get_interp_parameter('nlat'), c_int64_t)
 endif

 ival_c = get_dyn_grid_parm_codon(name_code, &
      int(ne, c_int64_t), int(np, c_int64_t), int(npsq, c_int64_t), &
      int(nelemd, c_int64_t), int(beglat, c_int64_t), &
      int(endlat, c_int64_t), int(ngcols_d, c_int64_t), &
      int(plat, c_int64_t), int(plev, c_int64_t), int(plevp, c_int64_t), &
      nlon_c, nlat_c)
 ival = int(ival_c)
 if (.not. proof_seen) then
    write(iulog,*) 'get_dyn_grid_parm implementation = codon'
    proof_seen = .true.
 endif
end function get_dyn_grid_parm

!#######################################################################

subroutine dyn_grid_get_pref(pref_edge, pref_mid, num_pr_lev)

   ! return reference pressures for the dynamics grid

   use pmgrid, only: plev
   use hycoef, only: hypi, hypm, nprlev

   ! arguments
   real(r8), intent(out) :: pref_edge(:) ! reference pressure at layer edges (Pa)
   real(r8), intent(out) :: pref_mid(:)  ! reference pressure at layer midpoints (Pa)
   integer,  intent(out) :: num_pr_lev   ! number of top levels using pure pressure representation

   integer :: k
   !-----------------------------------------------------------------------

   do k = 1, plev
      pref_edge(k) = hypi(k)
      pref_mid(k)  = hypm(k)
   end do
   pref_edge(plev+1) = hypi(plev+1)
   
   num_pr_lev = nprlev

end subroutine dyn_grid_get_pref

!#######################################################################

!-------------------------------------------------------------------------------
! This returns the lat/lon information (and corresponding MPI task numbers (owners)) 
! of the global model grid columns nearest to the input satellite coordinate (lat,lon)
!-------------------------------------------------------------------------------
subroutine dyn_grid_find_gcols( lat, lon, nclosest, owners, col, lbk, rlat, rlon, idyn_dists ) 
  use spmd_utils,     only: masterproc, iam
  use dimensions_mod, only: np
  use cam_abortutils, only: endrun
  use shr_const_mod,  only: SHR_CONST_PI, SHR_CONST_REARTH

  real(r8), intent(in) :: lat
  real(r8), intent(in) :: lon
  integer, intent(in)  :: nclosest
  integer, intent(out) :: owners(nclosest)
  integer, intent(out) :: col(nclosest)
  integer, intent(out) :: lbk(nclosest)

  real(r8),optional, intent(out) :: rlon(nclosest)
  real(r8),optional, intent(out) :: rlat(nclosest)
  real(r8),optional, intent(out) :: idyn_dists(nclosest)

  real(r8) :: dist            ! the distance (in radians**2 from lat, lon)
  real(r8) :: latr, lonr      ! lat, lon inputs converted to radians
  integer  :: blockid(1), bcid(1), lclblockid(1)
  integer  :: i,j,k, ii
  real(r8), allocatable :: clat_d(:), clon_d(:), distmin(:)
  integer, allocatable :: igcol(:)
  real(r8), parameter :: rad2deg = 180._r8/SHR_CONST_PI

  latr = lat/rad2deg
  lonr = lon/rad2deg

  allocate( clat_d(1:ngcols_d) )
  allocate( clon_d(1:ngcols_d) )
  allocate( igcol(nclosest) )
  allocate( distmin(nclosest) )

  call get_horiz_grid_d(ngcols_d, clat_d_out=clat_d, clon_d_out=clon_d)
  
  igcol(:)    = -999
  distmin(:) = 1.e10_r8

  do i = 1,ngcols_d
     
     ! Use the Spherical Law of Cosines to find the great-circle distance.
     dist = acos(sin(latr) * sin(clat_d(i)) + cos(latr) * cos(clat_d(i)) * cos(clon_d(i) - lonr)) * SHR_CONST_REARTH
     do j = nclosest, 1, -1
        if (dist < distmin(j)) then
           
           if (j < nclosest) then
              distmin(j+1) = distmin(j)
              igcol(j+1)    = igcol(j)
           end if
           
           distmin(j) = dist
           igcol(j)    = i
        else
           exit
        end if
     enddo
     
  enddo

  do i = 1,nclosest

     call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
     owners(i) = get_block_owner_d(blockid(1))
     
     if (owners(i)==iam) then
        lbk(i) = lclblockid(1)
        ii = igcol(i)-elem(lbk(i))%idxp%UniquePtoffset+1
        k=elem(lbk(i))%idxp%ia(ii)
        j=elem(lbk(i))%idxp%ja(ii)
        col(i) = k+(j-1)*np
     else
        lbk(i) = -1
        col(i) = -1
     endif
     if ( present(rlat) ) rlat(i) = clat_d(igcol(i)) * rad2deg
     if ( present(rlon) ) rlon(i) = clon_d(igcol(i)) * rad2deg
     
     if (present(idyn_dists)) then
        idyn_dists(i) = distmin(i)
     end if
     
  end do

  deallocate( clat_d )
  deallocate( clon_d )
  deallocate( igcol )
  deallocate( distmin )

 endsubroutine dyn_grid_find_gcols

!#######################################################################
subroutine dyn_grid_get_colndx( igcol, nclosest, owners, col, lbk ) 
  use spmd_utils, only: iam
  use dimensions_mod, only: np

  integer, intent(in)  :: nclosest
  integer, intent(in)  :: igcol(nclosest)
  integer, intent(out) :: owners(nclosest)
  integer, intent(out) :: col(nclosest)
  integer, intent(out) :: lbk(nclosest)

  integer  :: i,j,k, ii
  integer  :: blockid(1), bcid(1), lclblockid(1)

  do i = 1,nclosest

     call  get_gcol_block_d( igcol(i), 1, blockid, bcid, lclblockid )
     owners(i) = get_block_owner_d(blockid(1))
     
     if (owners(i)==iam) then
        lbk(i) = lclblockid(1)
        ii = igcol(i)-elem(lbk(i))%idxp%UniquePtoffset+1
        k=elem(lbk(i))%idxp%ia(ii)
        j=elem(lbk(i))%idxp%ja(ii)
        col(i) = k+(j-1)*np
     else
        lbk(i) = -1
        col(i) = -1
     endif
     
  end do

end subroutine dyn_grid_get_colndx
!#######################################################################
! this returns coordinates of a specified block element of the dyn grid
subroutine dyn_grid_get_elem_coords( ie, rlon, rlat, cdex )
  use dof_mod, only : UniqueCoords
  use dimensions_mod, only : np

  integer, intent(in) :: ie ! block element index

  real(r8),optional, intent(out) :: rlon(:) ! longitudes of the columns in the element
  real(r8),optional, intent(out) :: rlat(:) ! latitudes of the columns in the element
  integer, optional, intent(out) :: cdex(:) ! global column index

  integer :: sb,eb, ii, i,j, icol, igcol
  real(r8), allocatable :: clat(:), clon(:)

  sb = elem(ie)%idxp%UniquePtOffset
  eb = sb + elem(ie)%idxp%NumUniquePts-1

  allocate( clat(sb:eb), clon(sb:eb) )
  call UniqueCoords( elem(ie)%idxP, elem(ie)%spherep, clat(sb:eb), clon(sb:eb) )

  if (present(cdex)) cdex(:) = -1
  if (present(rlat)) rlat(:) = -999._r8
  if (present(rlon)) rlon(:) = -999._r8

  do ii=1,elem(ie)%idxp%NumUniquePts
     i=elem(ie)%idxp%ia(ii)
     j=elem(ie)%idxp%ja(ii)
     icol = i+(j-1)*np
     igcol = elem(ie)%idxp%UniquePtoffset+ii-1
     if (present(cdex)) cdex(icol) = igcol
     if (present(rlat)) rlat(icol) = clat( igcol )
     if (present(rlon)) rlon(icol) = clon( igcol )
  end do

  deallocate( clat, clon )

endsubroutine dyn_grid_get_elem_coords

!#######################################################################

end module dyn_grid
