
subroutine dadadj (lchnk   ,ncol    , &
                   pmid    ,pint    ,pdel    ,t       , &
                   q       )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! GFDL style dry adiabatic adjustment
! 
! Method: 
! if stratification is unstable, adjustment to the dry adiabatic lapse
! rate is forced subject to the condition that enthalpy is conserved.
! 
! Author: CMS Contact J.Hack
! 
!-----------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use perf_mod,        only: t_startf, t_stopf
   use ap_dadadj_scheme, only: dadadj_run
   implicit none

!
! Arguments
!
   integer, intent(in) :: lchnk               ! chunk identifier
   integer, intent(in) :: ncol                ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)   ! pressure at model levels
   real(r8), intent(in) :: pint(pcols,pverp)  ! pressure at model interfaces
   real(r8), intent(in) :: pdel(pcols,pver)   ! vertical delta-p

!
! Input/output arguments
!
   real(r8), intent(inout) :: t(pcols,pver)      ! temperature (K)
   real(r8), intent(inout) :: q(pcols,pver)      ! specific humidity

   call t_startf('ap_dadadj_run')
   call dadadj_run(lchnk, ncol, pcols, pver, pverp, &
                   pmid, pint, pdel, t, q)
   call t_stopf('ap_dadadj_run')
end subroutine dadadj
