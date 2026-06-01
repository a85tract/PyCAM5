#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module ll_mod
  implicit none
  private
  type :: node_t
     integer :: id
     integer :: Src,Dest
     logical :: valid
     type(node_t), pointer :: prev,next
  end type node_t

  type :: root_t
     integer     :: number
     type(node_t), pointer :: first
  end type root_t
  public :: node_t, root_t
  integer, public :: NumEdges 

  public :: PrintEdgeList
  public :: LLAddEdge,LLFindEdge, LLInsertEdge
  public :: LLSetEdgeCount,LLGetEdgeCount
  public :: LLFree

contains 

  subroutine LLSetEdgeCount(value)
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog
    implicit none
    integer,intent(in)   :: value
    integer(c_int64_t) :: value_c
    logical, save :: edge_count_proof_seen = .false.

    interface
       function llsetedgecount_codon(value_in_c) result(value_out_c) bind(c, name='llsetedgecount_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: value_in_c
         integer(c_int64_t) :: value_out_c
       end function llsetedgecount_codon
    end interface

#define SE_MISC_TAG 24
#define SE_MISC_LABEL 'll_mod'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

    value_c = llsetedgecount_codon(int(value, c_int64_t))
    NumEdges = int(value_c)
    if (.not. edge_count_proof_seen) then
       write(iulog,*) 'llsetedgecount implementation = codon'
       edge_count_proof_seen = .true.
    endif
  end subroutine LLSetEdgeCount

  subroutine LLGetEdgeCount(value)
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog
    implicit none
    integer,intent(out)  :: value
    integer(c_int64_t) :: value_c
    logical, save :: proof_seen = .false.

    interface
       function llgetedgecount_codon(value_in_c) result(value_out_c) bind(c, name='llgetedgecount_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: value_in_c
         integer(c_int64_t) :: value_out_c
       end function llgetedgecount_codon
    end interface

    value_c = llgetedgecount_codon(int(NumEdges, c_int64_t))
    value = int(value_c)
    if (.not. proof_seen) then
       write(iulog,*) 'llgetedgecount implementation = codon'
       proof_seen = .true.
    endif
  end subroutine LLGetEdgeCount

  subroutine PrintEdgeList(EdgeList)

    type(root_t) :: EdgeList(:)
    type(node_t), pointer :: temp_node
    integer :: nlist, i
    nlist = SIZE(EdgeList)

    do i=1,nlist
       temp_node => EdgeList(i)%first
       do while(associated(temp_node)) 
          print *,'Vertex: ',EdgeList(i)%number ,temp_node%Src,'->' ,temp_node%dest, '(',temp_node%id,')'
          temp_node => temp_node%next
       enddo
    enddo

  end subroutine PrintEdgeList

  subroutine LLFree(List)

    implicit none
    type(root_t) :: List
    type(node_t), pointer :: temp_node
    integer :: nlist,i


    temp_node => List%first
    ! Find the end of the list
    do while(associated(temp_node%next))
       temp_node => temp_node%next
    enddo

    temp_node => temp_node%prev
    !Now step back and deallocate all entries  
    do while(associated(temp_node))
       deallocate(temp_node%next)
       temp_node => temp_node%prev
    enddo

  end subroutine LLFree

  subroutine LLInsertEdge(EdgeList,src,dest,eNum)
    type (root_t), intent(inout) :: EdgeList
    integer, intent(in) :: src,dest
    integer, intent(out) :: eNum
    logical :: found

    call LLFindEdge(EdgeList,src,dest,eNum,found) 
    if(.not. found) then 
       call LLAddEdge(EdgeList,src,dest,eNum) 
    endif

  end subroutine LLInsertEdge

  subroutine LLFindEdge(Edge,src,dest,id,found)

    type (root_t), intent(in) :: Edge
    integer, intent(in)  :: src,dest
    integer, intent(out) :: id
    logical, intent(out) :: found

    type (node_t), pointer :: temp_node

    found =.FALSE.

    temp_node => Edge%first
    do while(associated(temp_node) .and. (.not. found))
       if((dest .eq. temp_node%dest) .and. (src .eq. temp_node%Src) ) then 
          found = .TRUE. 
          id=temp_node%id
       else
          temp_node => temp_node%next
       endif
    enddo
  end subroutine LLFindEdge

  subroutine LLAddEdge(EdgeList,src,dest,id)
    type (root_t), intent(inout) :: EdgeList
    integer, intent(in)  :: src
    integer, intent(in)  :: dest
    integer, intent(out)  :: id

    type(node_t), pointer :: temp_node
    type(node_t), pointer  :: new_node
    type(node_t), pointer :: parent

    temp_node => EdgeList%first
    parent    => EdgeList%first

    do while(associated(temp_node))
       parent => temp_node
       temp_node => parent%next
    enddo
    allocate(new_node)
    NumEdges = NumEdges + 1

    new_node%src=src
    new_node%dest=dest
    new_node%id=NumEdges
    NULLIFY(new_node%next)
    new_node%prev => parent

    if(associated(EdgeList%first)) then
       parent%next => new_node 
    else
       EdgeList%first => new_node 
    endif
    id = NumEdges

  end subroutine LLAddEdge

end module ll_mod
