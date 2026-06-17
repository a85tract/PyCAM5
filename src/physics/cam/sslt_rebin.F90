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
  logical :: sslt_rebin_register_logged = .false.
  logical :: sslt_rebin_init_logged = .false.
  logical :: sslt_rebin_adv_logged = .false.

  interface
     function sslt_rebin_has_four_codon(i1_c, i2_c, i3_c, i4_c) result(has_c) &
          bind(c, name="sslt_rebin_has_four_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: i1_c, i2_c, i3_c, i4_c
       integer(c_int64_t) :: has_c
     end function sslt_rebin_has_four_codon
     function sslt_rebin_active_codon(has_sslt_c) result(active_c) &
          bind(c, name="sslt_rebin_active_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: has_sslt_c
       integer(c_int64_t) :: active_c
     end function sslt_rebin_active_codon
     function sslt_rebin_register_codon(pcols_c, pver_c) result(mask_c) &
          bind(c, name="sslt_rebin_register_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: pcols_c, pver_c
       integer(c_int64_t) :: mask_c
     end function sslt_rebin_register_codon
     subroutine sslt_rebin_adv_codon(ncol_c, pver_c, pcols_c, wgt_sscm_c, &
          sslt1_p, sslt2_p, sslt3_p, sslt4_p, sslta_p, ssltc_p) &
          bind(c, name="sslt_rebin_adv_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pver_c, pcols_c
       real(c_double), value :: wgt_sscm_c
       type(c_ptr), value :: sslt1_p, sslt2_p, sslt3_p, sslt4_p
       type(c_ptr), value :: sslta_p, ssltc_p
     end subroutine sslt_rebin_adv_codon
  end interface

  private
  public :: sslt_rebin_init, sslt_rebin_adv, sslt_rebin_register
contains

  subroutine sslt_rebin_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (sslt_rebin_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('SSLT_REBIN_IMPL', impl_name, n, status)

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

  subroutine sslt_rebin_log_direct(logged, proof_line)

    logical, intent(inout) :: logged
    character(len=*), intent(in) :: proof_line

    if (logged) return
    logged = .true.

    if (masterproc) then
       write(iulog,'(A)') trim(proof_line)
    end if

  end subroutine sslt_rebin_log_direct


!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine sslt_rebin_register
    use ppgrid,       only : pver,pcols
    
    use physics_buffer, only : pbuf_add_field, dtype_r8

    integer(c_int64_t) :: register_mask_c

    call sslt_rebin_select_impl()
    if (.not. use_native_sslt_rebin_impl) then
       call sslt_rebin_proof_once()
       register_mask_c = sslt_rebin_register_codon(int(pcols, c_int64_t), int(pver, c_int64_t))
       if (register_mask_c > 0_c_int64_t) then
          call sslt_rebin_log_direct(sslt_rebin_register_logged, &
               'sslt_rebin_register direct = codon; register shape mask direct = codon; pbuf_add_field native CAM API island')
       end if
    end if

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
    integer(c_int64_t) :: has_c


    indices(1) = pbuf_get_index('sslt1',errcode)
    indices(2) = pbuf_get_index('sslt2',errcode)
    indices(3) = pbuf_get_index('sslt3',errcode)
    indices(4) = pbuf_get_index('sslt4',errcode)

    call sslt_rebin_select_impl()
    if (use_native_sslt_rebin_impl) then
       has_sslt = indices(1) > 0 .and. indices(2) > 0 .and. indices(3) > 0 .and. indices(4) > 0
    else
       call sslt_rebin_proof_once()
       has_c = sslt_rebin_has_four_codon(int(indices(1), c_int64_t), int(indices(2), c_int64_t), &
            int(indices(3), c_int64_t), int(indices(4), c_int64_t))
       has_sslt = has_c /= 0_c_int64_t
       call sslt_rebin_log_direct(sslt_rebin_init_logged, 'sslt_rebin_init direct = codon')
    end if
    if ( has_sslt ) source = DATA

    if ( .not. has_sslt ) then
       call cnst_get_ind ('SSLT01', indices(1), abort=.false.)
       call cnst_get_ind ('SSLT02', indices(2), abort=.false.)
       call cnst_get_ind ('SSLT03', indices(3), abort=.false.)
       call cnst_get_ind ('SSLT04', indices(4), abort=.false.)
       if (use_native_sslt_rebin_impl) then
          has_sslt = indices(1) > 0 .and. indices(2) > 0 .and. indices(3) > 0 .and. indices(4) > 0
       else
          has_c = sslt_rebin_has_four_codon(int(indices(1), c_int64_t), int(indices(2), c_int64_t), &
               int(indices(3), c_int64_t), int(indices(4), c_int64_t))
          has_sslt = has_c /= 0_c_int64_t
       end if
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
    use iso_c_binding, only : c_double, c_loc

    implicit none

    
    type(physics_state), target, intent(in) :: phys_state
    type(physics_buffer_desc), pointer :: pbuf(:)

!++ changed wgt_sscm declaration for roundoff validation with earlier code
!    real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode 
    real(r8), parameter :: wgt_sscm = 6.0_r8 / 7.0_r8 ! Fraction of total seasalt mass in coarse mode 

    real(r8), dimension(:,:), pointer, contiguous :: sslt1, sslt2, sslt3, sslt4
    real(r8), dimension(:,:), pointer, contiguous :: sslta, ssltc
    integer :: lchnk, ncol
    integer(c_int64_t) :: active_c
    real(r8) :: sslt_sum(pcols,pver)

    lchnk = phys_state%lchnk
    ncol = phys_state%ncol

    call sslt_rebin_select_impl()
    if (use_native_sslt_rebin_impl) then
       if (.not. has_sslt) return
    else
       call sslt_rebin_proof_once()
       active_c = sslt_rebin_active_codon(merge(1_c_int64_t, 0_c_int64_t, has_sslt))
       call sslt_rebin_log_direct(sslt_rebin_adv_logged, 'sslt_rebin_adv direct = codon')
       if (active_c == 0_c_int64_t) return
    end if

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

    if (use_native_sslt_rebin_impl) then
       sslt_sum(:ncol,:) = sslt1(:ncol,:) + sslt2(:ncol,:) + sslt3(:ncol,:) + sslt4(:ncol,:)
       sslta(:ncol,:) = (1._r8-wgt_sscm)*sslt_sum(:ncol,:) ! fraction of seasalt mass in accumulation mode
       ssltc(:ncol,:) = wgt_sscm*sslt_sum(:ncol,:) ! fraction of seasalt mass in coagulation mode
    else
       call sslt_rebin_adv_codon(int(ncol, c_int64_t), int(pver, c_int64_t), int(pcols, c_int64_t), &
            real(wgt_sscm, c_double), c_loc(sslt1(1,1)), c_loc(sslt2(1,1)), &
            c_loc(sslt3(1,1)), c_loc(sslt4(1,1)), c_loc(sslta(1,1)), c_loc(ssltc(1,1)))
    end if

    call outfld( 'SSLTA', sslta(:ncol,:), ncol, lchnk )
    call outfld( 'SSLTC', ssltc(:ncol,:), ncol, lchnk )

  end subroutine sslt_rebin_adv

end module sslt_rebin
