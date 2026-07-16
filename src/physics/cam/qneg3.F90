subroutine qneg3 (subnam  ,idx     ,ncol    ,ncold   ,lver    ,lconst_beg  , &
                  lconst_end       ,qmin    ,q       )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check moisture and tracers for minimum value, reset any below
! minimum value to minimum value and return information to allow
! warning message to be printed. The global average is NOT preserved.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_logfile,  only: iulog
   use cam_abortutils, only: endrun
   use qneg,         only: qneg_run
   use perf_mod,     only: t_startf, t_stopf
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   character*(*), intent(in) :: subnam ! name of calling routine

   integer, intent(in) :: idx          ! chunk/latitude index
   integer, intent(in) :: ncol         ! number of atmospheric columns
   integer, intent(in) :: ncold        ! declared number of atmospheric columns
   integer, intent(in) :: lver         ! number of vertical levels in column
   integer, intent(in) :: lconst_beg   ! beginning constituent
   integer, intent(in) :: lconst_end   ! ending    constituent

   real(r8), intent(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

!
! Input/Output arguments
!
   real(r8), intent(inout) :: q(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field
!
!---------------------------Local workspace-----------------------------
!
   integer nvals            ! number of values found < qmin
   integer m                ! constituent index
   integer iw,kw            ! i,k indices of worst violator
   integer errcode
   real(r8) worst           ! biggest violator
   character(len=512) :: errmsg

!
!-----------------------------------------------------------------------
!

   do m=lconst_beg,lconst_end
      call t_startf('ap_qneg_run')
      call qneg_run(ncol, ncold, lver, qmin(m), q(:,:,m), nvals, worst, &
           iw, kw, errcode, errmsg)
      call t_stopf('ap_qneg_run')
      if (errcode /= 0) call endrun(trim(errmsg))
      if (nvals>100 .and. abs(worst)>max(qmin(m),1.e-12_r8)) then
         write(iulog,9000)subnam,m,idx,nvals,qmin(m),worst,iw,kw
      end if
   end do
!
   return
9000 format(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i7, &
            ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
            1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)
end subroutine qneg3

