! The routine under test, extracted byte-for-byte from
! src/physics/cam/vertical_diffusion.F90 so that it compiles without the
! seventeen framework modules that file use-imports. run_threeway.py
! re-derives this extraction and fails if it no longer matches, so the
! copy cannot drift from the original.
!
! Only the environment is stood in for: pcols, pver and pcnst are the
! dimensions CAM use-imports from ppgrid and constituents. The routine
! text below is untouched.
module vd_ptend_core
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
  integer, parameter :: pcols = 8
  integer, parameter :: pver  = 30
  integer, parameter :: pcnst = 4
contains

  subroutine vertical_diffusion_ptend_core_native(ncol, psetcols_local, q_tmp_local, s_tmp_local, u_tmp_local, &
       v_tmp_local, state_q_local, state_s_local, state_u_local, state_v_local, sl_local, qt_local, &
       sl_prePBL_local, qt_prePBL_local, rztodt_local, ptend_q_local, ptend_s_local, ptend_u_local, ptend_v_local, &
       slten_local, qtten_local)

    integer, intent(in) :: ncol
    integer, intent(in) :: psetcols_local
    real(r8), intent(in) :: q_tmp_local(pcols,pver,pcnst)
    real(r8), intent(in) :: s_tmp_local(pcols,pver)
    real(r8), intent(in) :: u_tmp_local(pcols,pver)
    real(r8), intent(in) :: v_tmp_local(pcols,pver)
    real(r8), intent(in) :: state_q_local(pcols,pver,pcnst)
    real(r8), intent(in) :: state_s_local(pcols,pver)
    real(r8), intent(in) :: state_u_local(pcols,pver)
    real(r8), intent(in) :: state_v_local(pcols,pver)
    real(r8), intent(in) :: sl_local(pcols,pver)
    real(r8), intent(in) :: qt_local(pcols,pver)
    real(r8), intent(in) :: sl_prePBL_local(pcols,pver)
    real(r8), intent(in) :: qt_prePBL_local(pcols,pver)
    real(r8), intent(in) :: rztodt_local
    real(r8), intent(inout) :: ptend_q_local(psetcols_local,pver,pcnst)
    real(r8), intent(inout) :: ptend_s_local(psetcols_local,pver)
    real(r8), intent(inout) :: ptend_u_local(psetcols_local,pver)
    real(r8), intent(inout) :: ptend_v_local(psetcols_local,pver)
    real(r8), intent(inout) :: slten_local(pcols,pver)
    real(r8), intent(inout) :: qtten_local(pcols,pver)

    ptend_s_local(:ncol,:)       = ( s_tmp_local(:ncol,:) - state_s_local(:ncol,:) ) * rztodt_local
    ptend_u_local(:ncol,:)       = ( u_tmp_local(:ncol,:) - state_u_local(:ncol,:) ) * rztodt_local
    ptend_v_local(:ncol,:)       = ( v_tmp_local(:ncol,:) - state_v_local(:ncol,:) ) * rztodt_local
    ptend_q_local(:ncol,:pver,:) = ( q_tmp_local(:ncol,:pver,:) - state_q_local(:ncol,:pver,:) ) * rztodt_local
    slten_local(:ncol,:)         = ( sl_local(:ncol,:) - sl_prePBL_local(:ncol,:) ) * rztodt_local
    qtten_local(:ncol,:)         = ( qt_local(:ncol,:) - qt_prePBL_local(:ncol,:) ) * rztodt_local

  end subroutine vertical_diffusion_ptend_core_native

end module vd_ptend_core
