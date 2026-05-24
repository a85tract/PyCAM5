
module radheat
!-----------------------------------------------------------------------
!
! Purpose:  Provide an interface to convert shortwave and longwave
!           radiative heating terms into net heating.
!
!           This module provides a hook to allow incorporating additional
!           radiative terms (eUV heating and nonLTE longwave cooling).
! 
! Original version: B.A. Boville
!-----------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use ppgrid,        only: pcols, pver
use physics_types, only: physics_state, physics_ptend, physics_ptend_init
use cam_logfile,   only: iulog
use spmd_utils,    only: masterproc

use physics_buffer, only : physics_buffer_desc

implicit none
private
save

logical :: use_native_impl = .false.
logical :: impl_selected = .false.
logical :: use_native_tstep_init_impl = .false.
logical :: tstep_init_impl_selected = .false.
logical :: radheat_batch_use_native_impl = .false.
logical :: radheat_batch_impl_selected = .false.
logical :: radheat_batch_entered_logged = .false.
logical :: radheat_readnl_logged = .false.
logical :: radheat_init_logged = .false.
logical :: radheat_timestep_init_logged = .false.
logical :: radheat_tend_logged = .false.

interface
   subroutine radheat_readnl_codon() bind(c, name="radheat_readnl_codon")
   end subroutine radheat_readnl_codon

   subroutine radheat_init_codon() bind(c, name="radheat_init_codon")
   end subroutine radheat_init_codon

   subroutine radheat_timestep_init_codon() bind(c, name="radheat_timestep_init_codon")
   end subroutine radheat_timestep_init_codon

   subroutine radheat_batch_timestep_init_stage_dispatch_codon() &
        bind(c, name="radheat_batch_timestep_init_stage_dispatch_codon")
   end subroutine radheat_batch_timestep_init_stage_dispatch_codon

   subroutine radheat_tend_codon(ncol_c, pcols_c, pver_c, psetcols_c, &
        qrl_p, qrs_p, ptend_s_p, fsns_p, fsnt_p, flns_p, flnt_p, net_flx_p) &
        bind(c, name="radheat_tend_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, psetcols_c
      type(c_ptr), value :: qrl_p, qrs_p, ptend_s_p, fsns_p, fsnt_p, flns_p, flnt_p, net_flx_p
   end subroutine radheat_tend_codon

   subroutine radheat_batch_tend_stage_dispatch_codon(ncol_c, pcols_c, pver_c, psetcols_c, &
        qrl_p, qrs_p, ptend_s_p, fsns_p, fsnt_p, flns_p, flnt_p, net_flx_p) &
        bind(c, name="radheat_batch_tend_stage_dispatch_codon")
      use iso_c_binding, only: c_int64_t, c_ptr
      integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, psetcols_c
      type(c_ptr), value :: qrl_p, qrs_p, ptend_s_p, fsns_p, fsnt_p, flns_p, flnt_p, net_flx_p
   end subroutine radheat_batch_tend_stage_dispatch_codon
end interface

! Public interfaces
public  &
   radheat_readnl,        &!
   radheat_init,          &!
   radheat_timestep_init, &!
   radheat_tend            ! return net radiative heating

public :: radheat_disable_waccm ! disable waccm heating in the upper atm

!===============================================================================
contains
!===============================================================================

subroutine radheat_batch_append_proof(proof_line)

  character(len=*), intent(in) :: proof_line

  character(len=512) :: proof_file
  integer :: status, n, unitno

  proof_file = ''
  call get_environment_variable('RADHEAT_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
  if (status == 0 .and. n > 0) then
     open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
     write(unitno,'(A)') trim(proof_line)
     close(unitno)
  end if

end subroutine radheat_batch_append_proof

subroutine radheat_batch_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (radheat_batch_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('RADHEAT_BATCH_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     radheat_batch_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     radheat_batch_use_native_impl = .false.
  end if

  radheat_batch_impl_selected = .true.

  if (masterproc) then
     if (radheat_batch_use_native_impl) then
        write(iulog,*) 'radheat_batch implementation = native'
        call radheat_batch_append_proof('radheat_batch selector entered implementation = native')
     else
        write(iulog,*) 'radheat_batch implementation = codon'
        call radheat_batch_append_proof('radheat_batch selector entered implementation = codon')
     end if
     call flush(iulog)
  end if

end subroutine radheat_batch_select_impl

subroutine radheat_batch_log_entered()

  if (radheat_batch_entered_logged) return
  radheat_batch_entered_logged = .true.

  if (masterproc) then
     write(iulog,'(A)') 'radheat_batch entered (unified timestep/tend stage dispatch = codon)'
     call radheat_batch_append_proof('radheat_batch entered (unified timestep/tend stage dispatch = codon)')
     call flush(iulog)
  end if

end subroutine radheat_batch_log_entered

subroutine radheat_log_direct(logged, proof_line)

  logical, intent(inout) :: logged
  character(len=*), intent(in) :: proof_line

  if (logged) return
  logged = .true.

  if (masterproc) then
     write(iulog,'(A)') trim(proof_line)
     call radheat_batch_append_proof(proof_line)
     call flush(iulog)
  end if

end subroutine radheat_log_direct

subroutine radheat_readnl(nlfile)

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  call radheat_batch_select_impl()

  if (radheat_batch_use_native_impl) return

  call radheat_readnl_codon()
  call radheat_log_direct(radheat_readnl_logged, 'radheat_readnl direct = codon')

end subroutine radheat_readnl

!================================================================================================

subroutine radheat_init(pref_mid)

   use pmgrid, only: plev
   use physics_buffer, only : physics_buffer_desc

   real(r8), intent(in) :: pref_mid(plev)

   call radheat_batch_select_impl()

   if (radheat_batch_use_native_impl) return

   call radheat_init_codon()
   call radheat_log_direct(radheat_init_logged, 'radheat_init direct = codon')

end subroutine radheat_init

!================================================================================================

subroutine radheat_timestep_init (state, pbuf2d)
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    call radheat_batch_select_impl()

    if (radheat_batch_use_native_impl) then
       call radheat_timestep_init_native(state, pbuf2d)
       return
    end if

    call radheat_batch_log_entered()
    call radheat_timestep_init_codon()
    call radheat_log_direct(radheat_timestep_init_logged, 'radheat_timestep_init direct = codon')


end subroutine radheat_timestep_init

!================================================================================================

subroutine radheat_timestep_init_native (state, pbuf2d)
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc

    type(physics_state), intent(in):: state(begchunk:endchunk)                 
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


end subroutine radheat_timestep_init_native

!================================================================================================

subroutine radheat_tend(state, pbuf,  ptend, qrl, qrs, fsns, &
                        fsnt, flns, flnt, asdir, net_flx)
#if ( defined OFFLINE_DYN )
   use metdata, only: met_rlx, met_srf_feedback
#endif
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!-----------------------------------------------------------------------

! Arguments
   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   type(physics_state), intent(in)  :: state             ! Physics state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_ptend), intent(out) :: ptend             ! indivdual parameterization tendencie
   real(r8),            target, intent(in)  :: qrl(pcols,pver)   ! longwave heating
   real(r8),            target, intent(in)  :: qrs(pcols,pver)   ! shortwave heating
   real(r8),            target, intent(in)  :: fsns(pcols)       ! Surface solar absorbed flux
   real(r8),            target, intent(in)  :: fsnt(pcols)       ! Net column abs solar flux at model top
   real(r8),            target, intent(in)  :: flns(pcols)       ! Srf longwave cooling (up-down) flux
   real(r8),            target, intent(in)  :: flnt(pcols)       ! Net outgoing lw flux at model top
   real(r8),            intent(in)  :: asdir(pcols)      ! shortwave, direct albedo
   real(r8),            target, intent(out) :: net_flx(pcols)  


! Local variables
   integer :: ncol
   real(r8), target :: ptend_s_work(state%psetcols,pver)
!-----------------------------------------------------------------------

   call radheat_batch_select_impl()

#if ( defined OFFLINE_DYN )
   call radheat_tend_native(state, pbuf, ptend, qrl, qrs, fsns, fsnt, flns, flnt, asdir, net_flx)
   return
#endif

   if (radheat_batch_use_native_impl) then
      call radheat_tend_native(state, pbuf, ptend, qrl, qrs, fsns, fsnt, flns, flnt, asdir, net_flx)
      return
   end if

   ncol = state%ncol

   call physics_ptend_init(ptend,state%psetcols, 'radheat', ls=.true.)
   ptend_s_work = 0._r8

   call radheat_batch_log_entered()
   call radheat_tend_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(state%psetcols, c_int64_t), &
        c_loc(qrl), c_loc(qrs), c_loc(ptend_s_work), c_loc(fsns), c_loc(fsnt), c_loc(flns), c_loc(flnt), c_loc(net_flx) &
   )
   call radheat_log_direct(radheat_tend_logged, 'radheat_tend direct = codon')
   ptend%s = ptend_s_work

end subroutine radheat_tend

!================================================================================================
subroutine radheat_tend_native(state, pbuf,  ptend, qrl, qrs, fsns, &
                               fsnt, flns, flnt, asdir, net_flx)
#if ( defined OFFLINE_DYN )
   use metdata, only: met_rlx, met_srf_feedback
#endif
!-----------------------------------------------------------------------
! Compute net radiative heating from qrs and qrl, and the associated net
! boundary flux.
!-----------------------------------------------------------------------

! Arguments
   type(physics_state), intent(in)  :: state             ! Physics state variables
   
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_ptend), intent(out) :: ptend             ! indivdual parameterization tendencie
   real(r8),            intent(in)  :: qrl(pcols,pver)   ! longwave heating
   real(r8),            intent(in)  :: qrs(pcols,pver)   ! shortwave heating
   real(r8),            intent(in)  :: fsns(pcols)       ! Surface solar absorbed flux
   real(r8),            intent(in)  :: fsnt(pcols)       ! Net column abs solar flux at model top
   real(r8),            intent(in)  :: flns(pcols)       ! Srf longwave cooling (up-down) flux
   real(r8),            intent(in)  :: flnt(pcols)       ! Net outgoing lw flux at model top
   real(r8),            intent(in)  :: asdir(pcols)      ! shortwave, direct albedo
   real(r8),            intent(out) :: net_flx(pcols)


! Local variables
   integer :: i, k
   integer :: ncol
!-----------------------------------------------------------------------

   ncol = state%ncol

   call physics_ptend_init(ptend,state%psetcols, 'radheat', ls=.true.)

#if ( defined OFFLINE_DYN )
   ptend%s(:ncol,:) = 0._r8
   do k = 1,pver
     if (met_rlx(k) < 1._r8 .or. met_srf_feedback) then
       ptend%s(:ncol,k) = (qrs(:ncol,k) + qrl(:ncol,k))
     endif
   enddo 
#else
   ptend%s(:ncol,:) = (qrs(:ncol,:) + qrl(:ncol,:))
#endif

   do i = 1, ncol
      net_flx(i) = fsnt(i) - fsns(i) - flnt(i) + flns(i)
   end do

end subroutine radheat_tend_native

!================================================================================================
subroutine radheat_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('RADHEAT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_impl = .false.
  end if

  impl_selected = .true.

  if (masterproc) then
     if (use_native_impl) then
        write(iulog,*) 'radheat implementation = native'
     else
        write(iulog,*) 'radheat implementation = codon'
     end if
  end if

end subroutine radheat_select_impl

!================================================================================================
subroutine radheat_timestep_init_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tstep_init_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('RADHEAT_TSTEP_INIT_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_tstep_init_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_tstep_init_impl = .false.
  end if

  tstep_init_impl_selected = .true.

  if (masterproc) then
     if (use_native_tstep_init_impl) then
        write(iulog,*) 'radheat_timestep_init implementation = native'
     else
        write(iulog,*) 'radheat_timestep_init implementation = codon'
     end if
  end if

end subroutine radheat_timestep_init_select_impl

!================================================================================================
  subroutine radheat_disable_waccm()
  end subroutine radheat_disable_waccm
end module radheat
