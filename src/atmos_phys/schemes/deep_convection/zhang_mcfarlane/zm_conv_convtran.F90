module ap_zm_conv_convtran_scheme

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: zm_conv_convtran_run

contains

!> \section arg_table_zm_conv_convtran_run Argument Table
!! \htmlinclude zm_conv_convtran_run.html
subroutine zm_conv_convtran_run(lchnk   ,pcols  ,pver    ,nwtrc   ,nphase  , &
                    doconvtran,q       ,ncnst   ,constituent_is_dry, &
                    trace_water,nwater ,liq_indices,ice_indices, &
                    liq_rstd,ice_rstd,wtrc_qmin,mu      ,md      , &
                    du      ,eu      ,ed      ,dp      ,dsubcld , &
                    jt      ,mx      ,ideep   ,il1g    ,il2g    , &
                    nstep   ,fracis  ,dqdt    ,dpdry, Rwt )
!-----------------------------------------------------------------------
!
! Purpose:
! Convective transport of trace species
!
! Mixing ratios may be with respect to either dry or moist air
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: P. Rasch
!
! Rwt added by J. Nusbaumer to fix mystery variable passing error
!
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: pcols                 ! declared horizontal dimension
   integer, intent(in) :: pver                  ! vertical layer dimension
   integer, intent(in) :: nwtrc                 ! water tracer dimension
   integer, intent(in) :: nphase                ! liquid/ice water phase dimension
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   integer, intent(in) :: nwater                ! number of liquid/ice water tracers
   logical, intent(in) :: doconvtran(ncnst)     ! flag for doing convective transport
   logical, intent(in) :: constituent_is_dry(ncnst) ! true for dry-mass mixing ratios
   logical, intent(in) :: trace_water           ! true when water tracers are active
   integer, intent(in) :: liq_indices(nwtrc)    ! constituent index for each liquid tracer
   integer, intent(in) :: ice_indices(nwtrc)    ! constituent index for each ice tracer
   real(r8), intent(in) :: liq_rstd(nwtrc)      ! fallback liquid tracer ratios
   real(r8), intent(in) :: ice_rstd(nwtrc)      ! fallback ice tracer ratios
   real(r8), intent(in) :: wtrc_qmin            ! minimum base-water tendency
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Tracer array including moisture
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: dsubcld(pcols)       ! Delta pressure from cloud base to sfc
   real(r8), intent(in) :: fracis(pcols,pver,ncnst) ! fraction of tracer that is insoluble

   integer, intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: mx(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index

   real(r8), intent(in) :: dpdry(pcols,pver)       ! Delta pressure between interfaces


! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

!water tracers:
   !NOTE:  the 2 at the end is for liquid and ice - JN
   real(r8), intent(out) :: Rwt(pcols,pver,nwtrc,nphase) !water tracer ratio

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index

   real(r8) cabv                 ! Mix ratio of constituent above
   real(r8) cbel                 ! Mix ratio of constituent below
   real(r8) cdifr                ! Normalized diff between cabv and cbel
   real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
   real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
   real(r8) const(pcols,pver)    ! Gathered tracer array
   real(r8) fisg(pcols,pver)     ! gathered insoluble fraction of tracer
   real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
   real(r8) dcondt(pcols,pver)   ! Gathered tend array
   real(r8) small                ! A small number
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) mupdudp              ! A work variable
   real(r8) minc                 ! A work variable
   real(r8) maxc                 ! A work variable
   real(r8) fluxin               ! A work variable
   real(r8) fluxout              ! A work variable
   real(r8) netflux              ! A work variable

   real(r8) dutmp(pcols,pver)       ! Mass detraining from updraft
   real(r8) eutmp(pcols,pver)       ! Mass entraining from updraft
   real(r8) edtmp(pcols,pver)       ! Mass entraining from downdraft
   real(r8) dptmp(pcols,pver)    ! Delta pressure between interfaces

   !Water tracers?
   real(r8) chtmp(pcols,pver)       !Use for debugging...

!NOTE:  ixcldice = 3  - JN

!-----------------------------------------------------------------------
!
   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! Loop ever each constituent
   do m = 2, ncnst
      if (doconvtran(m)) then

         if (constituent_is_dry(m)) then
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dpdry(i,k)
                  dutmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eutmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  edtmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
               end do
            end do
         else
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dp(i,k)
                  dutmp(i,k) = du(i,k)
                  eutmp(i,k) = eu(i,k)
                  edtmp(i,k) = ed(i,k)
               end do
            end do
         endif
!        dptmp = dp

! Gather up the constituent and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc < 0) then
                  cdifr = 0._r8
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr > 1.E-6_r8) then
                  cabv = max(const(i,km1),maxc*1.e-12_r8)
                  cbel = max(const(i,k),maxc*1.e-12_r8)
                  chat(i,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))
               end if

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._r8

            end do
         end do

!Delete:
!         if(m .eq. 3) then
!           chtmp(:,:) = chat(i,k)
!         end if

!         if(m .eq. wtrc_iatype(3,iwtice)) then
!           write(iulog,*) 'chat ice 3 diff',sum(4._r8*chat(:il2g,:)-chtmp(:il2g,:))
!         end if

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         do i = il1g,il2g
            mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
            if (mupdudp > mbsth) then
               conu(i,kk) = (+eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1))/md(i,k)
            endif
         end do

! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
               if (mupdudp > mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eutmp(i,kk)*fisg(i,kk)* &
                                  const(i,kk)*dptmp(i,kk) )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-edtmp(i,km1)*fisg(i,km1)*const(i,km1) &
                                  *dptmp(i,km1) )/md(i,k)
               endif
            end do
         end do


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g

! version 1 hard to check for roundoff errors
!               dcondt(i,k) =
!     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
!     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
!     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
!     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
!     $                   )/dp(i,k)

! version 2 hard to limit fluxes
!               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)
!     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
!               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
!     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(i,kp1)*conu(i,kp1)+ mu(i,k)*min(chat(i,k),const(i,km1)) &
                         -(md(i,k)  *cond(i,k) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout = mu(i,k)*conu(i,k) + mu(i,kp1)*min(chat(i,kp1),const(i,k)) &
                         -(md(i,kp1)*cond(i,kp1) + md(i,k)*min(chat(i,k),const(i,k)))

               netflux = fluxin - fluxout
               if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                  netflux = 0._r8
               endif
               dcondt(i,k) = netflux/dptmp(i,k)
            end do
         end do
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!DIR$ NOINTERCHANGE
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k == mx(i)) then

! version 1
!                  dcondt(i,k) = (1./dsubcld(i))*
!     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
!     $               -md(i,k)*(cond(i,k)-chat(i,k))
!     $              )

! version 2
!                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
!                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
! version 3
                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1)) - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*min(chat(i,k),const(i,k))

                  netflux = fluxin - fluxout
                  if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                     netflux = 0._r8
                  endif
!                  dcondt(i,k) = netflux/dsubcld(i)
                  dcondt(i,k) = netflux/dptmp(i,k)
               else if (k > mx(i)) then
!                  dcondt(i,k) = dcondt(i,k-1)
                  dcondt(i,k) = 0._r8
               end if
            end do
         end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8
         do k = 1,pver
            kp1 = min(pver,k+1)
!DIR$ CONCURRENT
            do i = il1g,il2g
               dqdt(ideep(i),k,m) = dcondt(i,k)
            end do
         end do

      end if      ! for doconvtran

   end do

   if ( trace_water )then
!Calculate the water tracer ratio:
Rwt(:,:,:,:) = 1._r8 !initalize ratio
if(doconvtran(liq_indices(1))) then !are water tracers being transported?
  do m=2,nwater !loop over water tracers
    do k=1,pver
      do i=il1g,il2g
        Rwt(ideep(i),k,m,1) = tracer_ratio( &
                              dqdt(ideep(i),k,liq_indices(m)),   &
                              dqdt(ideep(i),k,liq_indices(1)),   &
                              liq_rstd(m),wtrc_qmin)
        Rwt(ideep(i),k,m,2) = tracer_ratio( &
                              dqdt(ideep(i),k,ice_indices(m)),   &
                              dqdt(ideep(i),k,ice_indices(1)),   &
                              ice_rstd(m),wtrc_qmin)
      end do
    end do
  end do
end if
   end if

   return
end subroutine zm_conv_convtran_run

pure real(r8) function tracer_ratio(qtrc, qtot, rstd, qmin)
   real(r8), intent(in) :: qtrc
   real(r8), intent(in) :: qtot
   real(r8), intent(in) :: rstd
   real(r8), intent(in) :: qmin

   if (abs(qtot) < qmin) then
      tracer_ratio = rstd
   else
      tracer_ratio = qtrc / qtot
   end if
end function tracer_ratio

end module ap_zm_conv_convtran_scheme
