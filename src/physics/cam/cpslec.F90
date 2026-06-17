
subroutine cpslec (ncol, pmid, phis, ps, t, psl, gravit, rair)

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver
  use cam_logfile, only: iulog
  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
  use spmd_utils, only: masterproc

  implicit none

  integer , intent(in) :: ncol
  real(r8), intent(in), target :: pmid(pcols,pver)
  real(r8), intent(in), target :: phis(pcols)
  real(r8), intent(in), target :: ps(pcols)
  real(r8), intent(in), target :: t(pcols,pver)
  real(r8), intent(out), target :: psl(pcols)
  real(r8), intent(in) :: gravit
  real(r8), intent(in) :: rair

  logical, save :: use_native_cpslec_impl = .false.
  logical, save :: cpslec_impl_selected = .false.
  logical, save :: cpslec_proof_written = .false.

  interface
     subroutine cpslec_codon(ncol_c, pmid_p, phis_p, ps_p, t_p, psl_p, gravit_c, &
          rair_c, pcols_c, pver_c) bind(c, name="cpslec_codon")
        use iso_c_binding, only: c_double, c_int64_t, c_ptr
        integer(c_int64_t), value :: ncol_c
        type(c_ptr), value :: pmid_p, phis_p, ps_p, t_p, psl_p
        real(c_double), value :: gravit_c
        real(c_double), value :: rair_c
        integer(c_int64_t), value :: pcols_c
        integer(c_int64_t), value :: pver_c
     end subroutine cpslec_codon
  end interface

  call cpslec_select_impl()

  if (use_native_cpslec_impl) then
     call cpslec_native(ncol, pmid, phis, ps, t, psl, gravit, rair)
     return
  end if

  call cpslec_proof_once()
  call cpslec_codon(int(ncol, c_int64_t), c_loc(pmid(1,1)), c_loc(phis(1)), &
       c_loc(ps(1)), c_loc(t(1,1)), c_loc(psl(1)), real(gravit, c_double), &
       real(rair, c_double), int(pcols, c_int64_t), int(pver, c_int64_t))

  return

contains

subroutine cpslec_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (cpslec_impl_selected) return

  impl_name = 'codon'
  call cam_codon_get_impl('CPSLEC_IMPL', impl_name, n, status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_cpslec_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_cpslec_impl = .false.
  end if

  cpslec_impl_selected = .true.

  if (masterproc) then
     if (use_native_cpslec_impl) then
        write(iulog,*) 'cpslec implementation = native'
     else
        write(iulog,*) 'cpslec implementation = codon'
     end if
  end if

end subroutine cpslec_select_impl

subroutine cpslec_proof_once()

  if (cpslec_proof_written) return
  cpslec_proof_written = .true.

  if (masterproc) then
     write(iulog,'(A)') 'cpslec entered (sea-level pressure loop = codon)'
  end if

end subroutine cpslec_proof_once

end subroutine cpslec

subroutine cpslec_native (ncol, pmid, phis, ps, t, psl, gravit, rair)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Hybrid coord version:  Compute sea level pressure for a latitude line
! 
! Method: 
! CCM2 hybrid coord version using ECMWF formulation
! Algorithm: See section 3.1.b in NCAR NT-396 "Vertical 
! Interpolation and Truncation of Model-Coordinate Data
!
! Author: Stolen from the Processor by Erik Kluzek
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid, only: pcols, pver

  implicit none

!-----------------------------Arguments---------------------------------
  integer , intent(in) :: ncol             ! longitude dimension

  real(r8), intent(in) :: pmid(pcols,pver) ! Atmospheric pressure (pascals)
  real(r8), intent(in) :: phis(pcols)      ! Surface geopotential (m**2/sec**2)
  real(r8), intent(in) :: ps(pcols)        ! Surface pressure (pascals)
  real(r8), intent(in) :: T(pcols,pver)    ! Vertical slice of temperature (top to bot)
  real(r8), intent(in) :: gravit           ! Gravitational acceleration
  real(r8), intent(in) :: rair             ! gas constant for dry air

  real(r8), intent(out):: psl(pcols)       ! Sea level pressures (pascals)
!-----------------------------------------------------------------------

!-----------------------------Parameters--------------------------------
  real(r8), parameter :: xlapse = 6.5e-3_r8   ! Temperature lapse rate (K/m)
!-----------------------------------------------------------------------

!-----------------------------Local Variables---------------------------
  integer i              ! Loop index
  real(r8) alpha         ! Temperature lapse rate in terms of pressure ratio (unitless)
  real(r8) Tstar         ! Computed surface temperature
  real(r8) TT0           ! Computed temperature at sea-level
  real(r8) alph          ! Power to raise P/Ps to get rate of increase of T with pressure
  real(r8) beta          ! alpha*phis/(R*T) term used in approximation of PSL
!-----------------------------------------------------------------------
!
  alpha = rair*xlapse/gravit
  do i=1,ncol
     if ( abs(phis(i)/gravit) < 1.e-4_r8 )then
        psl(i)=ps(i)
     else
        Tstar=T(i,pver)*(1._r8+alpha*(ps(i)/pmid(i,pver)-1._r8)) ! pg 7 eq 5

        TT0=Tstar + xlapse*phis(i)/gravit                  ! pg 8 eq 13

        if ( Tstar<=290.5_r8 .and. TT0>290.5_r8 ) then           ! pg 8 eq 14.1
           alph=rair/phis(i)*(290.5_r8-Tstar)  
        else if (Tstar>290.5_r8  .and. TT0>290.5_r8) then        ! pg 8 eq 14.2
           alph=0._r8
           Tstar= 0.5_r8 * (290.5_r8 + Tstar)  
        else  
           alph=alpha  
           if (Tstar<255._r8) then  
              Tstar= 0.5_r8 * (255._r8 + Tstar)                  ! pg 8 eq 14.3
           endif
        endif

        beta = phis(i)/(rair*Tstar)
        psl(i)=ps(i)*exp( beta*(1._r8-alph*beta/2._r8+((alph*beta)**2)/3._r8))
     end if
  enddo

  return
end subroutine cpslec_native
