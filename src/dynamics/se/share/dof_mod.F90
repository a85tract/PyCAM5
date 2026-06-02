#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dof_mod
  use kinds, only : real_kind,int_kind,long_kind
  use dimensions_mod, only : np, npsq, nelem, nelemd
  use quadrature_mod, only : quadrature_t
  use element_mod, only : element_t,index_t
  use parallel_mod, only : parallel_t, mpiinteger_t
  use edgetype_mod, only : longedgebuffer_t
  use edge_mod, only : initEdgebuffer,freeEdgebuffer, &
		       longedgevpack, longedgevunpackmin
  use bndry_mod, only : bndry_exchangev
implicit none
private
  ! public data
  ! public subroutines
  logical, save :: createuniqueindex_codon_logged = .false.
  public :: global_dof
  public :: genLocalDof
  public :: PrintDofP
  public :: UniquePoints
  public :: PutUniquePoints
  public :: UniqueNcolsP
  public :: UniqueCoords
  public :: CreateUniqueIndex
  public :: SetElemOffset
  public :: CreateMetaData

  interface UniquePoints
     module procedure UniquePoints2D
     module procedure UniquePoints3D
     module procedure UniquePoints4D
  end interface
  interface PutUniquePoints
     module procedure PutUniquePoints2D
     module procedure PutUniquePoints3D
     module procedure PutUniquePoints4D
  end interface


contains

  subroutine genLocalDof(ig,npts,ldof)
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use cam_logfile, only : iulog

    integer(kind=int_kind), intent(in) :: ig
    integer(kind=int_kind), intent(in) :: npts
    integer(kind=int_kind), intent(inout), target :: ldof(:,:)

    logical, save :: proof_seen = .false.
    interface
       subroutine genlocaldof_codon(ig_c, npts_c, ldof_p) bind(c, name='genlocaldof_codon')
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ig_c, npts_c
         type(c_ptr), value :: ldof_p
       end subroutine genlocaldof_codon
    end interface
  
    call genlocaldof_codon(int(ig, c_int64_t), int(npts, c_int64_t), c_loc(ldof))
    if (.not. proof_seen) then
       write(iulog,*) 'genlocaldof implementation = codon'
       proof_seen = .true.
    endif

  end subroutine genLocalDOF

! ===========================================
! global_dof
!
! Compute the global degree of freedom for each element...
! ===========================================

  subroutine global_dof(par,elem)
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog

    type (parallel_t),intent(in) :: par
    type (element_t)             :: elem(:)

    type (LongEdgeBuffer_t)    :: edge

    real(kind=real_kind)  da                     ! area element

    type (quadrature_t) :: gp

    integer (kind=int_kind) :: ldofP(np,np,nelemd)

    integer ii
    integer i,j,ig,ie
    integer kptr
    integer iptr

    ! ===================
    ! begin code
    ! ===================
#define SE_MISC_TAG 11
#define SE_MISC_LABEL 'dof_mod'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

    call initEdgeBuffer(edge,1)

    ! =================================================
    ! mass matrix on the velocity grid
    ! =================================================    

 
    do ie=1,nelemd
       ig = elem(ie)%vertex%number
       call genLocalDOF(ig,np,ldofP(:,:,ie))
	 
       kptr=0
       call LongEdgeVpack(edge,ldofP(:,:,ie),1,kptr,elem(ie)%desc)
    end do

    ! ==============================
    ! Insert boundary exchange here
    ! ==============================

    call bndry_exchangeV(par,edge)

    do ie=1,nelemd
       ! we should unpack directly into elem(ie)%gdofV, but we dont have
       ! a VunpackMIN that takes integer*8.  gdofV integer*8 means  
       ! more than 2G grid points.
       kptr=0
       call LongEdgeVunpackMIN(edge,ldofP(:,:,ie),1,kptr,elem(ie)%desc)
       elem(ie)%gdofP(:,:)=ldofP(:,:,ie)
    end do
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
    call FreeEdgeBuffer(edge)
       
  end subroutine global_dof


  subroutine UniquePoints2D(idxUnique,src,dest)
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use cam_logfile, only : iulog
    type (index_t), target :: idxUnique
    real (kind=real_kind), target :: src(:,:)
    real (kind=real_kind), target :: dest(:)

    logical, save :: proof_seen = .false.
    interface
       subroutine uniquepoints2d_codon(num_unique_pts_c, ia_p, ja_p, ni_c, src_p, dest_p) &
            bind(c, name='uniquepoints2d_codon')
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: num_unique_pts_c, ni_c
         type(c_ptr), value :: ia_p, ja_p, src_p, dest_p
       end subroutine uniquepoints2d_codon
    end interface
    
    call uniquepoints2d_codon(int(idxUnique%NumUniquePts, c_int64_t), &
         c_loc(idxUnique%ia), c_loc(idxUnique%ja), int(size(src, 1), c_int64_t), &
         c_loc(src), c_loc(dest))
    if (.not. proof_seen) then
       write(iulog,*) 'uniquepoints2d implementation = codon'
       proof_seen = .true.
    endif

  end subroutine UniquePoints2D

! putUniquePoints first zeros out the destination array, then fills the unique points of the 
! array with values from src.  A boundary communication should then be called to fill in the 
! redundent points of the array

  subroutine putUniquePoints2D(idxUnique,src,dest)
    type (index_t) :: idxUnique
    real (kind=real_kind),intent(in) :: src(:)
    real (kind=real_kind),intent(out) :: dest(:,:)

    integer(kind=int_kind) :: i,j,ii
    
    dest=0.0D0
    do ii=1,idxUnique%NumUniquePts
       i=idxUnique%ia(ii)
       j=idxUnique%ja(ii)
       dest(i,j)=src(ii)
    enddo

  end subroutine putUniquePoints2D

  subroutine UniqueNcolsP(elem,idxUnique,cid)    
    use element_mod, only : GetColumnIdP, element_t
    type (element_t), intent(in) :: elem
    type (index_t), intent(in) :: idxUnique
    integer,intent(out) :: cid(:)
    integer(kind=int_kind) :: i,j,ii


    do ii=1,idxUnique%NumUniquePts
       i=idxUnique%ia(ii)
       j=idxUnique%ja(ii)
       cid(ii)=GetColumnIdP(elem,i,j)
    enddo
    
  end subroutine UniqueNcolsP


  subroutine UniqueCoords(idxUnique,src,lat,lon)

    use coordinate_systems_mod, only  : spherical_polar_t
    type (index_t), intent(in) :: idxUnique

    type (spherical_polar_t) :: src(:,:)
    real (kind=real_kind), intent(out) :: lat(:)
    real (kind=real_kind), intent(out) :: lon(:)

    integer(kind=int_kind) :: i,j,ii

    do ii=1,idxUnique%NumUniquePts
       i=idxUnique%ia(ii)
       j=idxUnique%ja(ii)
       lat(ii)=src(i,j)%lat
       lon(ii)=src(i,j)%lon
    enddo

  end subroutine UniqueCoords

  subroutine UniquePoints3D(idxUnique,nlyr,src,dest)
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use cam_logfile, only : iulog
    type (index_t), target :: idxUnique
    integer(kind=int_kind) :: nlyr
    real (kind=real_kind), target, contiguous :: src(:,:,:)
    real (kind=real_kind), target, contiguous :: dest(:,:)

    logical, save :: proof_seen = .false.
    interface
       subroutine uniquepoints3d_codon(num_unique_pts_c, nlyr_c, ia_p, ja_p, ni_c, nj_c, &
            src_p, dest_p) bind(c, name='uniquepoints3d_codon')
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: num_unique_pts_c, nlyr_c, ni_c, nj_c
         type(c_ptr), value :: ia_p, ja_p, src_p, dest_p
       end subroutine uniquepoints3d_codon
    end interface

    call uniquepoints3d_codon(int(idxUnique%NumUniquePts, c_int64_t), int(nlyr, c_int64_t), &
         c_loc(idxUnique%ia(1)), c_loc(idxUnique%ja(1)), int(size(src, 1), c_int64_t), &
         int(size(src, 2), c_int64_t), c_loc(src(1,1,1)), c_loc(dest(1,1)))
    if (.not. proof_seen) then
       write(iulog,*) 'uniquepoints3d implementation = codon'
       proof_seen = .true.
    endif

  end subroutine UniquePoints3D
  subroutine UniquePoints4D(idxUnique,d3,d4,src,dest)
    type (index_t) :: idxUnique
    integer(kind=int_kind) :: d3,d4
    real (kind=real_kind) :: src(:,:,:,:)
    real (kind=real_kind) :: dest(:,:,:)
    
    integer(kind=int_kind) :: i,j,k,n,ii

    do n=1,d4
       do k=1,d3
          do ii=1,idxUnique%NumUniquePts
             i=idxUnique%ia(ii)
             j=idxUnique%ja(ii)
             dest(ii,k,n)=src(i,j,k,n)
          enddo
       end do
    enddo

  end subroutine UniquePoints4D

! putUniquePoints first zeros out the destination array, then fills the unique points of the 
! array with values from src.  A boundary communication should then be called to fill in the 
! redundent points of the array

  subroutine putUniquePoints3D(idxUnique,nlyr,src,dest)
    type (index_t) :: idxUnique
    integer(kind=int_kind) :: nlyr
    real (kind=real_kind),intent(in) :: src(:,:)
    real (kind=real_kind),intent(out) :: dest(:,:,:)
    
    integer(kind=int_kind) :: i,j,k,ii

    dest=0.0D0
    do k=1,nlyr
       do ii=1,idxUnique%NumUniquePts
          i=idxUnique%ia(ii)
          j=idxUnique%ja(ii)
          dest(i,j,k)=src(ii,k)
       enddo
    enddo

  end subroutine putUniquePoints3D

  subroutine putUniquePoints4D(idxUnique,d3,d4,src,dest)
    type (index_t) :: idxUnique
    integer(kind=int_kind) :: d3,d4
    real (kind=real_kind),intent(in) :: src(:,:,:)
    real (kind=real_kind),intent(out) :: dest(:,:,:,:)
    
    integer(kind=int_kind) :: i,j,k,n,ii

    dest=0.0D0
    do n=1,d4
       do k=1,d3
          do ii=1,idxunique%NumUniquePts
             i=idxUnique%ia(ii)
             j=idxUnique%ja(ii)
             dest(i,j,k,n)=src(ii,k,n)
          enddo
       enddo
    end do
  end subroutine putUniquePoints4D

  subroutine SetElemOffset(par,elem,GlobalUniqueColsP)
#ifdef _MPI
     use parallel_mod, only : mpi_sum
#endif
     type (parallel_t) :: par
     type (element_t) :: elem(:)
     integer, intent(out) :: GlobalUniqueColsP

     integer(kind=int_kind), allocatable :: numElemP(:),numElem2P(:)
     integer(kind=int_kind), allocatable :: numElemV(:),numElem2V(:)
     integer(kind=int_kind), allocatable :: gOffset(:)
    
     integer(kind=int_kind) :: ie,ig,nprocs,ierr

     logical,parameter :: Debug = .FALSE.

     nprocs = par%nprocs
     allocate(numElemP(nelem))
     allocate(numElem2P(nelem))
     allocate(gOffset(nelem))
     numElemP=0;numElem2P=0;gOffset=0

     do ie=1,nelemd
	ig = elem(ie)%GlobalId
	numElemP(ig) = elem(ie)%idxP%NumUniquePts
     enddo
#ifdef _MPI
     call MPI_Allreduce(numElemP,numElem2P,nelem,MPIinteger_t,MPI_SUM,par%comm,ierr) 
#else
     numElem2P=numElemP
#endif

     gOffset(1)=1
     do ig=2,nelem
	gOffset(ig) = gOffset(ig-1)+numElem2P(ig-1)
     enddo
     do ie=1,nelemd
        ig = elem(ie)%GlobalId
        elem(ie)%idxP%UniquePtOffset=gOffset(ig)
     enddo
     GlobalUniqueColsP = gOffset(nelem)+numElem2P(nelem)-1

     deallocate(numElemP)
     deallocate(numElem2P)
     deallocate(gOffset)
  end subroutine SetElemOffset

  subroutine CreateUniqueIndex(ig,gdof,idx)
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use cam_logfile, only : iulog
    use spmd_utils, only : masterproc

    integer(kind=int_kind) :: ig
    type (index_t), target :: idx
    integer(kind=long_kind), target :: gdof(:,:)
    
    integer :: npts

    interface
       function createuniqueindex_codon(ig_c, npts_c, gdof_p, ia_p, ja_p) result(num_unique_c) &
            bind(c, name="createuniqueindex_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ig_c, npts_c
         type(c_ptr), value :: gdof_p, ia_p, ja_p
         integer(c_int64_t) :: num_unique_c
       end function createuniqueindex_codon
    end interface


    npts = size(gdof,dim=1)
    idx%NumUniquePts = int(createuniqueindex_codon( &
         int(ig, c_int64_t), int(npts, c_int64_t), &
         c_loc(gdof(1,1)), c_loc(idx%ia(1)), c_loc(idx%ja(1))), int_kind)

    if (masterproc .and. .not. createuniqueindex_codon_logged) then
       write(iulog,*) 'createuniqueindex implementation = codon'
       createuniqueindex_codon_logged = .true.
       call flush(iulog)
    end if

  end subroutine CreateUniqueIndex


  subroutine CreateMetaData(par,elem,subelement_corners, fdofp)
    type (parallel_t),intent(in) :: par
    type (element_t), target    :: elem(:)

    integer, intent(out),optional         :: subelement_corners((np-1)*(np-1)*nelemd,4)
    integer(kind=int_kind), optional :: fdofp(np,np,nelemd)

    type (index_t), pointer  :: idx 
    type (LongEdgeBuffer_t)    :: edge
    integer :: i, j, ii, ie, base
    integer(kind=long_kind), pointer :: gdof(:,:)
    integer :: fdofp_local(np,np,nelemd)

    call initEdgeBuffer(edge,1)
    fdofp_local=0
    
    do ie=1,nelemd
       idx => elem(ie)%idxP
       do ii=1,idx%NumUniquePts
          i=idx%ia(ii)
          j=idx%ja(ii)
          
          fdofp_local(i,j,ie) = -(idx%UniquePtoffset+ii-1)
       end do
       call LongEdgeVpack(edge,fdofp_local(:,:,ie),1,0,elem(ie)%desc)
    end do
    call bndry_exchangeV(par,edge)
    do ie=1,nelemd
       base = (ie-1)*(np-1)*(np-1)
       call LongEdgeVunpackMIN(edge,fdofp_local(:,:,ie),1,0,elem(ie)%desc)
       if(present(subelement_corners)) then
          ii=0       
          do j=1,np-1
             do i=1,np-1
                ii=ii+1
                subelement_corners(base+ii,1) = -fdofp_local(i,j,ie)
                subelement_corners(base+ii,2) = -fdofp_local(i,j+1,ie)
                subelement_corners(base+ii,3) = -fdofp_local(i+1,j+1,ie)
                subelement_corners(base+ii,4) = -fdofp_local(i+1,j,ie)
             end do
          end do
       end if
    end do
    if(present(fdofp)) then
       fdofp=-fdofp_local
    end if
    


  end subroutine CreateMetaData


! ==========================================
!  PrintDofP
!
!   Prints the degree of freedom 
! ==========================================
  subroutine PrintDofP(elem)

   implicit none
   type (element_t), intent(in) :: elem(:)

   integer :: ie,nse,i,j
   

   nse = SIZE(elem)
 
   do ie=1,nse
      print *,'Element # ',elem(ie)%vertex%number
      do j=np,1,-1
         write(6,*) (elem(ie)%gdofP(i,j), i=1,np)
      enddo
   enddo
 10 format('I5')

 end subroutine PrintDofP

end module dof_mod
