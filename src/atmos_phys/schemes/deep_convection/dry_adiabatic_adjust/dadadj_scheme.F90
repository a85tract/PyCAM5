module ap_dadadj_scheme

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: dadadj_run

contains

!> \section arg_table_dadadj_run Argument Table
!! \htmlinclude dadadj_run.html
subroutine dadadj_run (lchnk   ,ncol    ,pcols   ,pver    ,pverp   , &
                   pmid    ,pint    ,pdel    ,t       , &
                   q       ,cappa   ,nlvdry  ,lat     ,lon     , &
                   errmsg  ,errflg  )
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
   implicit none

   integer niter           ! number of iterations for convergence
   parameter (niter = 15)

!
! Arguments
!
   integer, intent(in) :: lchnk               ! chunk identifier
   integer, intent(in) :: ncol                ! number of atmospheric columns
   integer, intent(in) :: pcols               ! declared horizontal dimension
   integer, intent(in) :: pver                ! vertical layer dimension
   integer, intent(in) :: pverp               ! vertical interface dimension
   integer, intent(in) :: nlvdry              ! deepest interface considered

   real(r8), intent(in) :: pmid(pcols,pver)   ! pressure at model levels
   real(r8), intent(in) :: pint(pcols,pverp)  ! pressure at model interfaces
   real(r8), intent(in) :: pdel(pcols,pver)   ! vertical delta-p
   real(r8), intent(in) :: cappa               ! R/cp for dry air
   real(r8), intent(in) :: lat(pcols)          ! host latitude index/value for diagnostics
   real(r8), intent(in) :: lon(pcols)          ! host longitude index/value for diagnostics

!
! Input/output arguments
!
   real(r8), intent(inout) :: t(pcols,pver)      ! temperature (K)
   real(r8), intent(inout) :: q(pcols,pver)      ! specific humidity
   character(len=512), intent(out) :: errmsg
   integer, intent(out) :: errflg
!
!---------------------------Local workspace-----------------------------
!
   integer i,k             ! longitude, level indices
   integer jiter           ! iteration index

   real(r8) c1dad(pver)        ! intermediate constant
   real(r8) c2dad(pver)        ! intermediate constant
   real(r8) c3dad(pver)        ! intermediate constant
   real(r8) c4dad(pver)        ! intermediate constant
   real(r8) gammad             ! dry adiabatic lapse rate (deg/Pa)
   real(r8) zeps               ! convergence criterion (deg/Pa)
   real(r8) rdenom             ! reciprocal of denominator of expression
   real(r8) dtdp               ! delta-t/delta-p
   real(r8) zepsdp             ! zeps*delta-p
   real(r8) zgamma             ! intermediate constant
   real(r8) qave               ! mean q between levels

   logical ilconv          ! .TRUE. ==> convergence was attained
   logical dodad(pcols)    ! .TRUE. ==> do dry adjustment
!
!-----------------------------------------------------------------------
!
   errmsg = ''
   errflg = 0

   if (pver < 2 .or. pverp < pver + 1 .or. nlvdry < 1 .or. nlvdry >= pver) then
      errmsg = 'dadadj_run: invalid vertical dimensions or nlvdry'
      errflg = 1
      return
   end if

   zeps = 2.0e-5_r8           ! set convergence criteria
!
! Find gridpoints with unstable stratification
!
   do i=1,ncol
      gammad = cappa*0.5_r8*(t(i,2) + t(i,1))/pint(i,2)
      dtdp = (t(i,2) - t(i,1))/(pmid(i,2) - pmid(i,1))
      dodad(i) = (dtdp + zeps) .gt. gammad
   end do
   do k=2,nlvdry
      do i=1,ncol
         gammad = cappa*0.5_r8*(t(i,k+1) + t(i,k))/pint(i,k+1)
         dtdp = (t(i,k+1) - t(i,k))/(pmid(i,k+1) - pmid(i,k))
         dodad(i) = dodad(i) .or. (dtdp + zeps).gt.gammad
      end do
   end do
!
! Make a dry adiabatic adjustment
! Note: nlvdry ****MUST**** be < pver
!
   do 80 i=1,ncol
      if (dodad(i)) then
         zeps = 2.0e-5_r8
         do k=1,nlvdry
            c1dad(k) = cappa*0.5_r8*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
            c2dad(k) = (1._r8 - c1dad(k))/(1._r8 + c1dad(k))
            rdenom = 1._r8/(pdel(i,k)*c2dad(k) + pdel(i,k+1))
            c3dad(k) = rdenom*pdel(i,k)
            c4dad(k) = rdenom*pdel(i,k+1)
         end do
50       do jiter=1,niter
            ilconv = .true.
            do k=1,nlvdry
               zepsdp = zeps*(pmid(i,k+1) - pmid(i,k))
               zgamma = c1dad(k)*(t(i,k) + t(i,k+1))
               if ((t(i,k+1)-t(i,k)) >= (zgamma+zepsdp)) then
                  ilconv = .false.
                  t(i,k+1) = t(i,k)*c3dad(k) + t(i,k+1)*c4dad(k)
                  t(i,k) = c2dad(k)*t(i,k+1)
                  qave = (pdel(i,k+1)*q(i,k+1) + pdel(i,k)*q(i,k))/(pdel(i,k+1)+ pdel(i,k))
                  q(i,k+1) = qave
                  q(i,k) = qave
               end if
            end do
            if (ilconv) go to 80 ! convergence => next longitude
         end do
!
! Double convergence criterion if no convergence in niter iterations
!
         zeps = zeps + zeps
         if (zeps > 1.e-4_r8) then
            write(errmsg,800) lat(i), lon(i), zeps
            errflg = 1
            return
         else
            go to 50
         end if
      end if
80    continue
      return
!
! Formats
!
800   format('dadadj_run: no convergence at lat/lon=',2f12.5,', zeps=',e9.4)
end subroutine dadadj_run

end module ap_dadadj_scheme
