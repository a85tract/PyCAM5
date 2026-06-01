!===========================================================================
! Combine several reactions into one pseudo reaction to correct the 
! photolysis rate J(O1D) to incorporate the effect of the other reactions. 
!
! Creator: Philip Cameron-Smith
!===========================================================================

module llnl_O1D_to_2OH_adj

  use shr_kind_mod, only : r8 => shr_kind_r8
  use iso_c_binding, only : c_int64_t

  implicit none

  private
  public :: O1D_to_2OH_adj, O1D_to_2OH_adj_init

  logical :: o1d_to_2oh_adj_use_native_impl = .false.
  logical :: o1d_to_2oh_adj_impl_selected = .false.
  logical :: o1d_to_2oh_adj_init_proof_written = .false.

  integer :: jo1d_ndx

  interface
     function o1d_to_2oh_adj_init_active_codon(active) result(out_c) bind(c, name="o1d_to_2oh_adj_init_active_codon")
       use iso_c_binding, only : c_int64_t
       integer(c_int64_t), value :: active
       integer(c_int64_t) :: out_c
     end function o1d_to_2oh_adj_init_active_codon
  end interface

contains
!===========================================================================

!===========================================================================
!===========================================================================
  subroutine O1D_to_2OH_adj_init
    use mo_chem_utls, only : get_rxt_ndx
    use cam_logfile,  only : iulog
    use spmd_utils,       only : masterproc

    implicit none
    integer(c_int64_t) :: active_c

    jo1d_ndx  = get_rxt_ndx( 'j2oh' )
    active_c = o1d_to_2oh_adj_init_active_codon(merge(1_c_int64_t, 0_c_int64_t, jo1d_ndx > 0))
    if (masterproc) then
       write (iulog,*) 'O1D_to_2OH_adj_init: Found j2oh index in O1D_to_2OH_adj_init of   ', jo1d_ndx
       if (active_c == 0_c_int64_t) then
          write (iulog,'(A)') 'o1d_to_2oh_adj_init direct = codon missing-j2oh no-op'
       else
          write (iulog,'(A)') 'o1d_to_2oh_adj_init selector = codon; active j2oh setup body = native'
       end if
       call flush(iulog)
    endif

  end subroutine O1D_to_2OH_adj_init

!===========================================================================
!===========================================================================
  subroutine O1D_to_2OH_adj( p_rate, inv, m, ncol, tfld )

    use chem_mods,    only : nfs, phtcnt, rxntot, nfs !PJC added rxntot, nfs
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr
    use ppgrid,       only : pcols, pver              !PJC added pcols
    use mo_setinv,    only : n2_ndx, o2_ndx, h2o_ndx  !PJC

    implicit none

    !--------------------------------------------------------------------
    ! ... dummy arguments
    !--------------------------------------------------------------------
    integer,  intent(in) :: ncol
    real(r8), target, intent(in) :: inv(ncol,pver,nfs)
    real(r8), intent(in) :: m(ncol,pver)
    real(r8), target, intent(inout) :: p_rate(ncol,pver,rxntot)
    real(r8), target, intent(in)    :: tfld(pcols,pver)               ! midpoint temperature (K)

    interface
       subroutine O1D_to_2OH_adj_codon(ncol_c, pcols_c, pver_c, rxntot_c, nfs_c, jo1d_ndx_c, n2_ndx_c, o2_ndx_c, &
            h2o_ndx_c, p_rate_p, inv_p, tfld_p) bind(c, name="O1D_to_2OH_adj_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, rxntot_c, nfs_c
         integer(c_int64_t), value :: jo1d_ndx_c, n2_ndx_c, o2_ndx_c, h2o_ndx_c
         type(c_ptr), value :: p_rate_p, inv_p, tfld_p
       end subroutine O1D_to_2OH_adj_codon
    end interface

    call O1D_to_2OH_adj_select_impl()

    if (.not. o1d_to_2oh_adj_use_native_impl) then
       call O1D_to_2OH_adj_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
            int(rxntot, c_int64_t), int(nfs, c_int64_t), int(jo1d_ndx, c_int64_t), int(n2_ndx, c_int64_t), &
            int(o2_ndx, c_int64_t), int(h2o_ndx, c_int64_t), c_loc(p_rate), c_loc(inv), c_loc(tfld))
       return
    end if

    call O1D_to_2OH_adj_native( p_rate, inv, m, ncol, tfld )

  end subroutine O1D_to_2OH_adj

  subroutine O1D_to_2OH_adj_native( p_rate, inv, m, ncol, tfld )

    use chem_mods,    only : nfs, phtcnt, rxntot, nfs !PJC added rxntot, nfs
    use ppgrid,       only : pcols, pver              !PJC added pcols
    use mo_setinv,    only : n2_ndx, o2_ndx, h2o_ndx  !PJC

    implicit none

    !--------------------------------------------------------------------
    ! ... dummy arguments
    !--------------------------------------------------------------------
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: inv(ncol,pver,nfs)
    real(r8), intent(in) :: m(ncol,pver)
    real(r8), intent(inout) :: p_rate(ncol,pver,rxntot)
    real(r8), intent(in)    :: tfld(pcols,pver)               ! midpoint temperature (K)

    !--------------------------------------------------------------------
    ! ... local variables
    !--------------------------------------------------------------------
    integer :: k
    real(r8) :: im(ncol)
    real(r8) :: n2_rate(ncol,pver)
    real(r8) :: o2_rate(ncol,pver)
    real(r8) :: h2o_rate(ncol,pver)

    real(r8), parameter :: x1 = 2.15e-11_r8
    real(r8), parameter :: x2 = 3.30e-11_r8
    real(r8), parameter :: x3 = 1.63e-10_r8
    real(r8), parameter :: y1 = 110.0_r8
    real(r8), parameter :: y2 =  55.0_r8
    real(r8), parameter :: y3 =  60.0_r8

    if (jo1d_ndx<1) return

    n2_rate(:,:)  = x1 * Exp ( y1 / tfld(:ncol,:)) * inv(:,:,n2_ndx)
    o2_rate(:,:)  = x2 * Exp ( y2 / tfld(:ncol,:)) * inv(:,:,o2_ndx)
    h2o_rate(:,:) = x3 * Exp ( y3 / tfld(:ncol,:)) * inv(:,:,h2o_ndx)

    p_rate(:,:,jo1d_ndx) = p_rate(:,:,jo1d_ndx) *   &
                          (h2o_rate(:,:) / (h2o_rate(:,:) + n2_rate(:,:) + o2_rate(:,:)))

  end subroutine O1D_to_2OH_adj_native

  subroutine O1D_to_2OH_adj_select_impl()

    use cam_logfile, only : iulog
    use spmd_utils, only : masterproc

    implicit none

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (o1d_to_2oh_adj_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('O1D_TO_2OH_ADJ_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       o1d_to_2oh_adj_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       o1d_to_2oh_adj_use_native_impl = .false.
    end if

    o1d_to_2oh_adj_impl_selected = .true.

    if (masterproc) then
       if (o1d_to_2oh_adj_use_native_impl) then
          write(iulog,*) 'O1D_to_2OH_adj implementation = native'
       else
          write(iulog,*) 'O1D_to_2OH_adj implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine O1D_to_2OH_adj_select_impl

end module llnl_O1D_to_2OH_adj
