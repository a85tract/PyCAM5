
module constituent_burden

!-----------------------------------------------------------------------------------------
! Purpose: subroutines to generate constituent burden history variables
!
! Revision history:
! 2005-12-21  K. Lindsay       Original version
!-----------------------------------------------------------------------------------------

  use constituents, only: pcnst
  use cam_logfile, only: iulog
  use iso_c_binding, only: c_double, c_int64_t, c_loc
  use spmd_utils, only: masterproc

  implicit none

! Public interfaces

  public constituent_burden_init
  public constituent_burden_comp

  private

  character(len=18) :: burdennam(pcnst)     ! name of burden history variables
  logical :: use_native_constituent_burden_impl = .false.
  logical :: constituent_burden_impl_selected = .false.
  logical :: constituent_burden_proof_written = .false.
  logical :: constituent_burden_comp_logged = .false.
  logical :: constituent_burden_init_logged = .false.

  interface
     function constituent_burden_flag_codon(flag_c) result(active_c) &
          bind(c, name="constituent_burden_flag_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: active_c
     end function constituent_burden_flag_codon
     function constituent_burden_init_codon(flag_c) result(active_c) &
          bind(c, name="constituent_burden_init_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: active_c
     end function constituent_burden_init_codon
     function constituent_burden_comp_codon(flag_c) result(active_c) &
          bind(c, name="constituent_burden_comp_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: active_c
     end function constituent_burden_comp_codon
     subroutine constituent_burden_comp_integral_codon(ncol_c, psetcols_c, pver_c, pcnst_c, m_c, dry_c, &
          rga_c, q_p, pdel_p, pdeldry_p, ftem_p) bind(c, name="constituent_burden_comp_integral_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, psetcols_c, pver_c, pcnst_c, m_c, dry_c
       real(c_double), value :: rga_c
       type(c_ptr), value :: q_p, pdel_p, pdeldry_p, ftem_p
     end subroutine constituent_burden_comp_integral_codon
  end interface

  save

!=========================================================================================

contains

!=========================================================================================

subroutine constituent_burden_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (constituent_burden_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('CONSTITUENT_BURDEN_IMPL', value=impl_name, &
       length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_constituent_burden_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_constituent_burden_impl = .false.
  end if

  constituent_burden_impl_selected = .true.

  if (masterproc) then
     if (use_native_constituent_burden_impl) then
        write(iulog,*) 'constituent_burden implementation = native'
     else
        write(iulog,*) 'constituent_burden implementation = codon'
     end if
  end if

end subroutine constituent_burden_select_impl

subroutine constituent_burden_proof_once()

  if (constituent_burden_proof_written) return
  constituent_burden_proof_written = .true.

  if (masterproc) then
     write(iulog,'(A)') 'constituent_burden entered (history active gate = codon)'
  end if

end subroutine constituent_burden_proof_once

subroutine constituent_burden_log_comp_direct()

  if (constituent_burden_comp_logged) return
  constituent_burden_comp_logged = .true.

  if (masterproc) then
     write(iulog,'(A)') 'constituent_burden_comp direct = codon; column burden integral direct = codon; hist_fld_active/outfld native CAM API boundary'
  end if

end subroutine constituent_burden_log_comp_direct

subroutine constituent_burden_log_direct(logged, proof_line)

  logical, intent(inout) :: logged
  character(len=*), intent(in) :: proof_line

  if (logged) return
  logged = .true.

  if (masterproc) then
     write(iulog,'(A)') trim(proof_line)
     call flush(iulog)
  end if

end subroutine constituent_burden_log_direct

logical function constituent_burden_flag(flag) result(active)
  logical, intent(in) :: flag
  integer(c_int64_t) :: flag_c

  call constituent_burden_select_impl()

  if (use_native_constituent_burden_impl) then
     active = flag
     return
  end if

  call constituent_burden_proof_once()
  flag_c = 0_c_int64_t
  if (flag) flag_c = 1_c_int64_t
  active = constituent_burden_flag_codon(flag_c) /= 0_c_int64_t

end function constituent_burden_flag

!=========================================================================================

subroutine constituent_burden_init

  use cam_history,   only: addfld, phys_decomp
  use constituents,  only: cnst_name

  integer :: m
  if (constituent_burden_init_codon(1_c_int64_t) == 0_c_int64_t) return
  call constituent_burden_log_direct(constituent_burden_init_logged, &
       'constituent_burden_init direct = codon; history registration native CAM API island')

  do m = 2, pcnst
    burdennam(m) = 'TM'//cnst_name(m)
    call addfld (burdennam(m), 'kg/m2', 1, 'A', &
                 trim(cnst_name(m)) // ' column burden', phys_decomp)
  end do

end subroutine constituent_burden_init

!=========================================================================================

subroutine constituent_burden_comp(state)

  use physics_types, only: physics_state
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use constituents,  only: cnst_type
  use ppgrid,        only: pcols
  use physconst,     only: rga
  use cam_history,   only: outfld, hist_fld_active

!-----------------------------------------------------------------------
!
! Arguments
!
   type(physics_state), target, intent(inout) :: state
!
!---------------------------Local workspace-----------------------------

  real(r8), target :: ftem(pcols)      ! temporary workspace

  integer :: m, lchnk, ncol
  integer(c_int64_t) :: active_c

  call constituent_burden_select_impl()
  if (use_native_constituent_burden_impl) then
     call constituent_burden_comp_native(state)
     return
  end if

  call constituent_burden_proof_once()
  active_c = constituent_burden_comp_codon(1_c_int64_t)
  if (active_c == 0_c_int64_t) return

  lchnk = state%lchnk
  ncol  = state%ncol
  call constituent_burden_log_comp_direct()

  do m = 2, pcnst
     if (.not. constituent_burden_flag(hist_fld_active(burdennam(m)))) cycle
     call constituent_burden_comp_integral_codon(int(ncol, c_int64_t), int(state%psetcols, c_int64_t), &
          int(size(state%q, 2), c_int64_t), int(size(state%q, 3), c_int64_t), int(m, c_int64_t), &
          merge(1_c_int64_t, 0_c_int64_t, cnst_type(m) .eq. 'dry'), real(rga, c_double), &
          c_loc(state%q), c_loc(state%pdel), c_loc(state%pdeldry), c_loc(ftem))
     call outfld (burdennam(m), ftem, pcols, lchnk)
  end do

end subroutine constituent_burden_comp

!=========================================================================================

subroutine constituent_burden_comp_native(state)

  use physics_types, only: physics_state
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use constituents,  only: cnst_type
  use ppgrid,        only: pcols
  use physconst,     only: rga
  use cam_history,   only: outfld, hist_fld_active

   type(physics_state), intent(inout) :: state

  real(r8) :: ftem(pcols)

  integer :: m, lchnk, ncol

  lchnk = state%lchnk
  ncol  = state%ncol

  do m = 2, pcnst
     if (.not. hist_fld_active(burdennam(m))) cycle
     if (cnst_type(m) .eq. 'dry') then
        ftem(:ncol) = sum(state%q(:ncol,:,m) * state%pdeldry(:ncol,:), dim=2) * rga
     else
        ftem(:ncol) = sum(state%q(:ncol,:,m) * state%pdel(:ncol,:), dim=2) * rga
     endif
     call outfld (burdennam(m), ftem, pcols, lchnk)
  end do

end subroutine constituent_burden_comp_native

!=========================================================================================

end module constituent_burden
