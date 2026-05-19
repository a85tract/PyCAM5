!-------------------------------------------------------------------
! rebins the 4 sea salt bins into 2 bins for the radiation
!
!  N.B. This code looks for the constituents of SSLTA and SSLTC
!       in the physics buffer first, and uses those if found.
!       Consequently, it is not possible to have prognostic sea
!       salt be radiatively active if the prescribed sea salt is
!       also present.  The current (cam3_5_52) chemistry configurations
!       don't allow both prescribed and prognostic to be present
!       simultaneously, but a more flexible chemistry package that
!       allows this would break this code.
!
! Created by: Francis Vitt
! Date: 9 May 2008
!-------------------------------------------------------------------
module sslt_rebin

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_logfile,    only: iulog
  use iso_c_binding,  only: c_int64_t
  use spmd_utils,     only: masterproc

  implicit none

  integer :: indices(4)
  integer :: sslta_idx, ssltc_idx

  logical :: has_sslt = .false.
  character(len=1) :: source
  character(len=1), parameter :: DATA = 'D'
  character(len=1), parameter :: PROG = 'P'
  logical :: use_native_sslt_rebin_impl = .false.
  logical :: sslt_rebin_impl_selected = .false.
  logical :: sslt_rebin_proof_written = .false.

  interface
     function sslt_rebin_has_four_codon(i1_c, i2_c, i3_c, i4_c) result(has_c) &
          bind(c, name="sslt_rebin_has_four_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: i1_c, i2_c, i3_c, i4_c
       integer(c_int64_t) :: has_c
     end function sslt_rebin_has_four_codon
  end interface

  private
  public :: sslt_rebin_init, sslt_rebin_adv, sslt_rebin_register
contains

  subroutine sslt_rebin_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (sslt_rebin_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('SSLT_REBIN_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_sslt_rebin_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_sslt_rebin_impl = .false.
    end if

    sslt_rebin_impl_selected = .true.

    if (masterproc) then
       if (use_native_sslt_rebin_impl) then
          write(iulog,*) 'sslt_rebin implementation = native'
       else
          write(iulog,*) 'sslt_rebin implementation = codon'
       end if
    end if

  end subroutine sslt_rebin_select_impl

  subroutine sslt_rebin_proof_once()

    if (sslt_rebin_proof_written) return
    sslt_rebin_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'sslt_rebin entered (source-index helper = codon)'
    end if

  end subroutine sslt_rebin_proof_once

  logical function sslt_rebin_has_four(i1, i2, i3, i4) result(has_four)
    integer, intent(in) :: i1, i2, i3, i4

    call sslt_rebin_select_impl()

    if (use_native_sslt_rebin_impl) then
       has_four = i1 > 0 .and. i2 > 0 .and. i3 > 0 .and. i4 > 0
       return
    end if

    call sslt_rebin_proof_once()
    has_four = sslt_rebin_has_four_codon(int(i1, c_int64_t), &
         int(i2, c_int64_t), int(i3, c_int64_t), int(i4, c_int64_t)) /= 0_c_int64_t

  end function sslt_rebin_has_four


!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_register
    use ppgrid,       only : pver,pcols
    
    use physics_buffer, only : pbuf_add_field, dtype_r8

    ! add SSLTA and SSLTC to physics buffer
    call pbuf_add_field('SSLTA','physpkg',dtype_r8,(/pcols,pver/),sslta_idx)
    call pbuf_add_field('SSLTC','physpkg',dtype_r8,(/pcols,pver/),ssltc_idx)

  endsubroutine sslt_rebin_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_init()

    use constituents, only : cnst_get_ind
    
    use physics_buffer, only : pbuf_get_index, pbuf_set_field, physics_buffer_desc
    use ppgrid,       only : pver
    use cam_history,  only : addfld, phys_decomp

    implicit none

    integer :: errcode


    indices(1) = pbuf_get_index('sslt1',errcode)
    indices(2) = pbuf_get_index('sslt2',errcode)
    indices(3) = pbuf_get_index('sslt3',errcode)
    indices(4) = pbuf_get_index('sslt4',errcode)

    has_sslt = sslt_rebin_has_four(indices(1), indices(2), indices(3), indices(4))
    if ( has_sslt ) source = DATA

    if ( .not. has_sslt ) then
       call cnst_get_ind ('SSLT01', indices(1), abort=.false.)
       call cnst_get_ind ('SSLT02', indices(2), abort=.false.)
       call cnst_get_ind ('SSLT03', indices(3), abort=.false.)
       call cnst_get_ind ('SSLT04', indices(4), abort=.false.)
       has_sslt = sslt_rebin_has_four(indices(1), indices(2), indices(3), indices(4))
       if ( has_sslt ) source = PROG
    endif

    if ( has_sslt ) then
       call addfld('SSLTA','kg/kg', pver, 'A', 'sea salt', phys_decomp )
       call addfld('SSLTC','kg/kg', pver, 'A', 'sea salt', phys_decomp )
    endif

  end subroutine sslt_rebin_init
  
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_adv(pbuf,  phys_state)

    use physics_types,only : physics_state
    
    use ppgrid,       only : pver, pcols
    use cam_history,  only : outfld
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    implicit none

    
    type(physics_state), target, intent(in) :: phys_state
    type(physics_buffer_desc), pointer :: pbuf(:)

!++ changed wgt_sscm declaration for roundoff validation with earlier code
!    real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode 
    real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode 

    real(r8), dimension(:,:), pointer :: sslt1, sslt2, sslt3, sslt4
    real(r8), dimension(:,:), pointer :: sslta, ssltc
    integer :: lchnk, ncol
    real(r8) :: sslt_sum(pcols,pver)

    lchnk = phys_state%lchnk
    ncol = phys_state%ncol

    if (.not. has_sslt) return

    select case( source )
    case (PROG)
       sslt1 => phys_state%q(:,:,indices(1))
       sslt2 => phys_state%q(:,:,indices(2))
       sslt3 => phys_state%q(:,:,indices(3))
       sslt4 => phys_state%q(:,:,indices(4))
    case (DATA)
       call pbuf_get_field(pbuf, indices(1), sslt1)
       call pbuf_get_field(pbuf, indices(2), sslt2)
       call pbuf_get_field(pbuf, indices(3), sslt3)
       call pbuf_get_field(pbuf, indices(4), sslt4)
    end select

    call pbuf_get_field(pbuf, sslta_idx, sslta )
    call pbuf_get_field(pbuf, ssltc_idx, ssltc )

    sslt_sum(:ncol,:) = sslt1(:ncol,:) + sslt2(:ncol,:) + sslt3(:ncol,:) + sslt4(:ncol,:)
    sslta(:ncol,:) = (1._r8-wgt_sscm)*sslt_sum(:ncol,:) ! fraction of seasalt mass in accumulation mode
    ssltc(:ncol,:) = wgt_sscm*sslt_sum(:ncol,:) ! fraction of seasalt mass in coagulation mode

    call outfld( 'SSLTA', sslta(:ncol,:), ncol, lchnk )
    call outfld( 'SSLTC', ssltc(:ncol,:), ncol, lchnk )

  end subroutine sslt_rebin_adv

end module sslt_rebin
