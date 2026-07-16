
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
   use phys_grid,       only: get_lat_p, get_lon_p
   use physconst,       only: cappa
   use cam_control_mod, only: nlvdry
   use cam_abortutils,  only: endrun
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

   real(r8) :: lat(pcols), lon(pcols)
   character(len=512) :: errmsg
   integer :: errflg, i

   do i = 1, ncol
      lat(i) = real(get_lat_p(lchnk, i), r8)
      lon(i) = real(get_lon_p(lchnk, i), r8)
   end do

   call t_startf('ap_dadadj_run')
   call dadadj_run(lchnk, ncol, pcols, pver, pverp, &
                   pmid, pint, pdel, t, q, cappa, nlvdry, lat, lon, &
                   errmsg, errflg)
   call t_stopf('ap_dadadj_run')
   if (errflg /= 0) call endrun(trim(errmsg))
end subroutine dadadj
