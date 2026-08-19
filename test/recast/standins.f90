! Stand-ins for the framework modules dadadj.F90 use-imports, so the file
! compiles standalone for the module-level gate. Values here must match the
! shim's (test/recast/dadadj_shim.py): cappa is spelled as the same
! expression on both sides so the float64 arithmetic agrees to the bit.
module shr_kind_mod
  implicit none
  integer, parameter :: shr_kind_r8 = selected_real_kind(12)
end module shr_kind_mod

module ppgrid
  implicit none
  integer, parameter :: pcols = 8
  integer, parameter :: pver  = 30
  integer, parameter :: pverp = 31
end module ppgrid

module physconst
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  ! (SHR_CONST_RGAS / MWDAIR) / CPDAIR, spelled as its definition.
  real(r8), parameter :: cappa = (6.02214e26_r8 * 1.38065e-23_r8 / 28.966_r8) / 1.00464e3_r8
end module physconst

module cam_control_mod
  implicit none
  integer :: nlvdry = 3
end module cam_control_mod

module cam_logfile
  implicit none
  integer :: iulog = 6
end module cam_logfile

module spmd_utils
  implicit none
  logical :: masterproc = .false.
end module spmd_utils

module phys_grid
  implicit none
contains
  integer function get_lat_p(lchnk, icol)
    integer, intent(in) :: lchnk, icol
    get_lat_p = 0
  end function get_lat_p
  integer function get_lon_p(lchnk, icol)
    integer, intent(in) :: lchnk, icol
    get_lon_p = 0
  end function get_lon_p
end module phys_grid

module cam_abortutils
  implicit none
contains
  subroutine endrun(msg)
    character(len=*), intent(in), optional :: msg
    ! The gate calls dadadj_native directly and compares state, so a
    ! non-convergence "abort" becomes a return: both sides completed the
    ! same fixed iteration count and their state is still comparable.
  end subroutine endrun
end module cam_abortutils

! Link-time stubs for the dispatcher path the gate never takes.
subroutine cam_codon_get_impl(name, value, n, status)
  character(len=*), intent(in) :: name
  character(len=*), intent(inout) :: value
  integer, intent(out) :: n, status
  n = 0
  status = 1
end subroutine cam_codon_get_impl

subroutine dadadj_codon_stub_guard()
end subroutine dadadj_codon_stub_guard

subroutine dadadj_codon(ncol_c, pcols_c, nlvdry_c, cappa_c, &
     pmid_p, pint_p, pdel_p, t_p, q_p, &
     c1dad_p, c2dad_p, c3dad_p, c4dad_p, dodad_p, &
     status_p, zeps_fail_p, fail_i_p) bind(c, name="dadadj_codon")
  use iso_c_binding, only: c_double, c_int64_t, c_ptr
  integer(c_int64_t), value :: ncol_c, pcols_c, nlvdry_c
  real(c_double), value :: cappa_c
  type(c_ptr), value :: pmid_p, pint_p, pdel_p, t_p, q_p
  type(c_ptr), value :: c1dad_p, c2dad_p, c3dad_p, c4dad_p, dodad_p
  type(c_ptr), value :: status_p, zeps_fail_p, fail_i_p
  ! Never called: the gate compares dadadj_native (Fortran) against the
  ! Codon library invoked directly from Python. This satisfies the linker
  ! for the dispatcher object that rides along in dadadj.F90.
end subroutine dadadj_codon
