
subroutine qneg4 (subnam  ,lchnk   ,ncol    ,ztodt   ,        &
                  qbot    ,srfrpdel,shflx   ,lhflx   ,qflx    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check if moisture flux into the ground is exceeding the total
! moisture content of the lowest model layer (creating negative moisture
! values).  If so, then subtract the excess from the moisture and
! latent heat fluxes and add it to the sensible heat flux.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Olson
! 
! Water isotopes added by J. Nusbaumer - Mar 2011
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,    only: get_lat_p, get_lon_p
   use physconst,    only: gravit, latvap
   use constituents, only: qmin, pcnst
   use cam_logfile,  only: iulog
   use cam_abortutils, only: endrun
   use qneg4_scheme, only: qneg4_run
   use perf_mod, only: t_startf, t_stopf

   !water isotopes:
   use water_types,   only: iwtvap
   use water_tracer_vars, only: trace_water, wtrc_iatype, wtrc_ntype, iwspec
   use water_tracers, only: wtrc_ratio

   implicit none

!
! Input arguments
!
   character*8, intent(in) :: subnam         ! name of calling routine
!
   integer, intent(in) :: lchnk              ! chunk index
   integer, intent(in) :: ncol               ! number of atmospheric columns
!
   real(r8), intent(in) :: ztodt             ! two times model timestep (2 delta-t)
   real(r8), intent(in) :: qbot(pcols,pcnst) ! moisture at lowest model level
   real(r8), intent(in) :: srfrpdel(pcols)   ! 1./(pint(K+1)-pint(K))
!
! Input/Output arguments
!
   real(r8), intent(inout) :: shflx(pcols)   ! Surface sensible heat flux (J/m2/s)
   real(r8), intent(inout) :: lhflx(pcols)   ! Surface latent   heat flux (J/m2/s)
   real(r8), intent(inout) :: qflx (pcols,pcnst)   ! surface water flux (kg/m^2/s)
!
!---------------------------Local workspace-----------------------------
!
   integer :: iw                ! i index of worst violator
   integer :: nptsexc           ! number of points with excess flux
   integer :: i, ii, ivap       ! column and isotope indices
   integer :: m, errcode
   integer :: indxexc(pcols)
   real(r8) :: excess(pcols), qfxo(pcols,pcnst), rat, worst
   character(len=512) :: errmsg

!
!-----------------------------------------------------------------------
!

   do m = 1, pcnst
      do i = 1, ncol
         qfxo(i,m) = qflx(i,m)
      end do
   end do

   call t_startf('ap_qneg4_run')
   call qneg4_run(ncol, pcnst, ztodt, qbot, srfrpdel, qmin(1), &
        gravit, latvap, shflx, lhflx, qflx, nptsexc, &
        indxexc(:ncol), excess(:ncol), worst, iw, errcode, errmsg)
   if (errcode /= 0) call endrun(trim(errmsg))

   if (nptsexc.gt.10) then
      write(iulog,9000) subnam,nptsexc,worst, lchnk, iw, &
           get_lat_p(lchnk,iw),get_lon_p(lchnk,iw)
   end if

   if (trace_water) then
      do ivap = 1, wtrc_ntype(iwtvap)
         m = wtrc_iatype(ivap, iwtvap)
         do ii = 1, nptsexc
            i = indxexc(ii)
            rat = wtrc_ratio(iwspec(m), qfxo(i,m), &
                 qfxo(i,wtrc_iatype(1,iwtvap)))
            qflx(i,m) = qflx(i,m) - rat*excess(i)
         end do
      end do
   end if
   call t_stopf('ap_qneg4_run')

   return
9000 format(' QNEG4 WARNING from ',a8 &
            ,' Max possible LH flx exceeded at ',i6,' points. ' &
            ,', Worst excess = ',1pe12.4 &
            ,', lchnk = ',i6 &
            ,', i = ',i6 &
            ,', same as indices lat =', i6 &
            ,', lon =', i6 &
           )
end subroutine qneg4
