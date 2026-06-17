#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module mass_matrix_mod
  use kinds, only : real_kind
  use dimensions_mod, only : np, nelemd
  use quadrature_mod, only : quadrature_t, gauss ,gausslobatto
  use element_mod, only : element_t
  use parallel_mod, only : parallel_t
  use edgetype_mod, only : edgebuffer_t
  use edge_mod, only : edgevpack,edgevunpack, freeedgebuffer,initedgebuffer  
  use bndry_mod, only : bndry_exchangev
implicit none
private

  public :: mass_matrix

contains

! ===========================================
! mass_matrix:
!
! Compute the mass matrix for each element...
! ===========================================

  subroutine mass_matrix(par,elem)
    use iso_c_binding, only : c_int64_t, c_ptr, c_loc
    use cam_logfile, only : iulog

    type (parallel_t),intent(in) :: par
    type (element_t), target :: elem(:)

    type (EdgeBuffer_t)    :: edge

    real(kind=real_kind)  da                     ! area element

    type (quadrature_t) :: gp

    integer ii
    integer i,j
    integer kptr
    integer iptr
    character(len=32) :: impl_name
    integer :: impl_n, impl_status
    logical :: use_codon_impl
    logical, save :: proof_seen = .false.

    interface
       subroutine mass_matrix_vgrid_init_codon(np_c, mp_p, rmp_p, weights_p) &
            bind(c, name='mass_matrix_vgrid_init_codon')
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: np_c
         type(c_ptr), value :: mp_p, rmp_p, weights_p
       end subroutine mass_matrix_vgrid_init_codon
       subroutine mass_matrix_invert_codon(np_c, field_p) &
            bind(c, name='mass_matrix_invert_codon')
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: np_c
         type(c_ptr), value :: field_p
       end subroutine mass_matrix_invert_codon
       subroutine mass_matrix_sphere_init_codon(np_c, mp_p, metdet_p, spheremp_p, rspheremp_p) &
            bind(c, name='mass_matrix_sphere_init_codon')
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: np_c
         type(c_ptr), value :: mp_p, metdet_p, spheremp_p, rspheremp_p
       end subroutine mass_matrix_sphere_init_codon
    end interface

    ! ===================
    ! begin code
    ! ===================

#define SE_MISC_TAG 27
#define SE_MISC_LABEL 'mass_matrix_mod'
! Codon evidence: bind(c, name='se_misc_touch_codon') and SE_MISC_HELPERS_IMPL selector are in se_codon_misc_touch.inc.
#include "se_codon_misc_touch.inc"
#undef SE_MISC_LABEL
#undef SE_MISC_TAG

    call initEdgeBuffer(par,edge,elem,1,nthreads=1)

    impl_name = 'codon'
    call cam_codon_get_impl('MASS_MATRIX_IMPL', impl_name, impl_n, impl_status)
    use_codon_impl = .not. (impl_status == 0 .and. impl_n > 0 .and. &
         trim(adjustl(impl_name(:impl_n))) == 'native')

    ! =================================================
    ! mass matrix on the velocity grid
    ! =================================================    

    gp=gausslobatto(np)
 
    do ii=1,nelemd
       if (use_codon_impl) then
          call mass_matrix_vgrid_init_codon(int(np, c_int64_t), c_loc(elem(ii)%mp(1,1)), &
               c_loc(elem(ii)%rmp(1,1)), c_loc(gp%weights(1)))
       else
          do j=1,np
             do i=1,np
                 ! MNL: metric term for map to reference element is now in metdet!
                elem(ii)%mp(i,j)=gp%weights(i)*gp%weights(j)
                elem(ii)%rmp(i,j)=elem(ii)%mp(i,j)
             end do
          end do
       endif

       kptr=0
       call edgeVpack(edge,elem(ii)%rmp,1,kptr,ii)

    end do

    ! ==============================
    ! Insert boundary exchange here
    ! ==============================

    call bndry_exchangeV(par,edge)

    do ii=1,nelemd

       kptr=0
       call edgeVunpack(edge,elem(ii)%rmp,1,kptr,ii)

       if (use_codon_impl) then
          call mass_matrix_invert_codon(int(np, c_int64_t), c_loc(elem(ii)%rmp(1,1)))
       else
          do j=1,np
             do i=1,np
                elem(ii)%rmp(i,j)=1.0D0/elem(ii)%rmp(i,j)
             end do
          end do
       endif

    end do
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif

    deallocate(gp%points)
    deallocate(gp%weights)

    ! =============================================
    ! compute spherical element mass matrix
    ! =============================================
    do ii=1,nelemd
       if (use_codon_impl) then
          call mass_matrix_sphere_init_codon(int(np, c_int64_t), c_loc(elem(ii)%mp(1,1)), &
               c_loc(elem(ii)%metdet(1,1)), c_loc(elem(ii)%spheremp(1,1)), &
               c_loc(elem(ii)%rspheremp(1,1)))
       else
          do j=1,np
             do i=1,np
                elem(ii)%spheremp(i,j)=elem(ii)%mp(i,j)*elem(ii)%metdet(i,j)
                elem(ii)%rspheremp(i,j)=elem(ii)%spheremp(i,j)
             end do
          end do
       endif
       kptr=0
       call edgeVpack(edge,elem(ii)%rspheremp,1,kptr,ii)
    end do
    call bndry_exchangeV(par,edge)
    do ii=1,nelemd
       kptr=0
       call edgeVunpack(edge,elem(ii)%rspheremp,1,kptr,ii)
       if (use_codon_impl) then
          call mass_matrix_invert_codon(int(np, c_int64_t), c_loc(elem(ii)%rspheremp(1,1)))
       else
          do j=1,np
             do i=1,np
                elem(ii)%rspheremp(i,j)=1.0D0/elem(ii)%rspheremp(i,j)
             end do
          end do
       endif
    end do
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif

    ! =============================================
    ! compute the mass matrix 
    ! =============================================
    ! Jose Garcia: Not sure but I think this code is just dead code
    !do ii=1,nelemd
    !   iptr=1
    !   do j=1,np
    !      do i=1,np
    !         elem(ii)%mp(i,j)=elem(ii)%mp(i,j)
    !         iptr=iptr+1
    !      end do
    !   end do
    !end do

    call FreeEdgeBuffer(edge)

    if (use_codon_impl .and. .not. proof_seen) then
       write(iulog,*) 'mass_matrix implementation = codon'
       proof_seen = .true.
    endif
       
  end subroutine mass_matrix

end module mass_matrix_mod
