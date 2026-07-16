module ap_cldfrc_scheme

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  real(r8), parameter :: pnot = 1.e5_r8
  real(r8), parameter :: lapse = 6.5e-3_r8
  real(r8), parameter :: pretop = 1.e2_r8

  public :: cldfrc_run

contains

  !> \section arg_table_cldfrc_run Argument Table
  !! \htmlinclude cldfrc_run.html
  subroutine cldfrc_run( &
       pcols, pver, ncol, top_lev, k700, dindex, iceopt, &
       cldfrc_freeze_dry, cldfrc_ice, inversion_cld_off, &
       shallow_fraction_provided, &
       rhminl, rhminl_adj_land, rhminh, rhminp, rhminp_botmb, &
       sh1, sh2, dp1, dp2, &
       premit, premib, icecrit, unset_r8, cappa, gravit, rair, pi, &
       pref_mid, clat, pmid, temp, q, phis, shfrc, cmfmc, cmfmc2, &
       landfrac, snowh, sst, ps, ocnfrac, cldice, qs, esl, esi, &
       cloud, rhcloud, clc, concld, cldst, rhu00, icecldf, liqcldf, &
       relhum, shallowcu, deepcu, cold_sst_index)

    integer, intent(in) :: pcols, pver, ncol, top_lev, k700, dindex, iceopt
    logical, intent(in) :: cldfrc_freeze_dry, cldfrc_ice
    logical, intent(in) :: inversion_cld_off, shallow_fraction_provided
    real(r8), intent(in) :: rhminl, rhminl_adj_land, rhminh, rhminp
    real(r8), intent(in) :: rhminp_botmb, sh1, sh2, dp1, dp2
    real(r8), intent(in) :: premit, premib, icecrit, unset_r8
    real(r8), intent(in) :: cappa, gravit, rair, pi
    real(r8), intent(in) :: pref_mid(:), clat(:)
    real(r8), intent(in) :: pmid(:,:), temp(:,:), q(:,:), phis(:)
    real(r8), intent(in) :: shfrc(:,:), cmfmc(:,:), cmfmc2(:,:)
    real(r8), intent(in) :: landfrac(:), snowh(:), sst(:), ps(:)
    real(r8), intent(in) :: ocnfrac(:), cldice(:,:)
    real(r8), intent(in) :: qs(:,:), esl(:,:), esi(:,:)
    real(r8), intent(out) :: cloud(:,:), rhcloud(:,:), clc(:)
    real(r8), intent(out) :: concld(:,:), cldst(:,:), rhu00(:,:)
    real(r8), intent(out) :: icecldf(:,:), liqcldf(:,:), relhum(:,:)
    real(r8), intent(out) :: shallowcu(:,:), deepcu(:,:)
    integer, intent(out) :: cold_sst_index

    real(r8) :: dthdpmn(pcols)
    real(r8) :: dthdp
    real(r8) :: rhwght
    real(r8) :: rh(pcols,pver)
    real(r8) :: rhdif
    real(r8) :: strat
    real(r8) :: theta(pcols,pver)
    real(r8) :: rhlim
    real(r8) :: coef1
    real(r8) :: rpdeli(pcols,pver-1)
    real(r8) :: rhpert
    logical :: cldbnd(pcols)
    integer :: i, k, kp1, numkcld
    integer :: kdthdp(pcols)
    real(r8) :: a, b, c, as, bs, cs
    real(r8) :: Kc, ttmp, icicval, rho, ncf, phi
    real(r8) :: thetas(pcols)

    logical :: land
    land(i) = nint(landfrac(i)) == 1

    ! Initialise cloud fraction
    shallowcu = 0._r8
    deepcu    = 0._r8

    !==================================================================================
    ! PHILOSOPHY OF PRESENT IMPLEMENTATION
    !++ag ice3
    ! Modification to philosophy for ice supersaturation
    ! philosophy below is based on RH water only. This is 'liquid condensation'
    ! or liquid cloud (even though it will freeze immediately to ice)
    ! The idea is that the RH limits for condensation are strict only for
    ! water saturation
    !
    ! Ice clouds are formed by explicit parameterization of ice nucleation.
    ! Closure for ice cloud fraction is done on available cloud ice, such that
    ! the in-cloud ice content matches an empirical fit
    ! thus, icecldf = min(cldice/icicval,1) where icicval = f(temp,cldice,numice)
    ! for a first cut, icicval=f(temp) only.
    ! Combined cloud fraction is maximum overlap  cloud=max(1,max(icecldf,liqcldf))
    ! No dA/dt term for ice?
    !--ag
    !
    ! There are three co-existing cloud types: convective, inversion related low-level
    ! stratocumulus, and layered cloud (based on relative humidity).  Layered and
    ! stratocumulus clouds do not compete with convective cloud for which one creates
    ! the most cloud.  They contribute collectively to the total grid-box average cloud
    ! amount.  This is reflected in the way in which the total cloud amount is evaluated
    ! (a sum as opposed to a logical "or" operation)
    !
    !==================================================================================
    ! set defaults for rhu00
    rhu00(:,:) = 2.0_r8
    ! define rh perturbation in order to estimate rhdfda
    rhpert = 0.01_r8

    !set Wang and Sassen IWC paramters
    a=26.87_r8
    b=0.569_r8
    c=0.002892_r8
    !set schiller parameters
    as=-68.4202_r8
    bs=0.983917_r8
    cs=2.81795_r8
    !set wood and field paramters...
    Kc=75._r8

    cloud    = 0._r8
    icecldf  = 0._r8
    liqcldf  = 0._r8
    rhcloud  = 0._r8
    cldst    = 0._r8
    concld   = 0._r8

    do k=top_lev,pver
       theta(:ncol,k)    = temp(:ncol,k)*(pnot/pmid(:ncol,k))**cappa

       do i=1,ncol
          rh(i,k)     = q(i,k)/qs(i,k)*(1.0_r8+real(dindex,r8)*rhpert)
          !  record relhum, rh itself will later be modified related with concld
          relhum(i,k) = rh(i,k)
       end do
    end do

    ! Initialize other temporary variables
    cold_sst_index = 0
    do i=1,ncol
       ! Adjust thetas(i) in the presence of non-zero ocean heights.
       ! This reduces the temperature for positive heights according to a standard lapse rate.
       if(ocnfrac(i).gt.0.01_r8) thetas(i)  = &
            ( sst(i) - lapse * phis(i) / gravit) * (pnot/ps(i))**cappa
       if(ocnfrac(i).gt.0.01_r8.and.sst(i).lt.260._r8) cold_sst_index = i
       clc(i) = 0.0_r8
    end do
    coef1 = gravit*864.0_r8    ! conversion to millibars/day

    do k=top_lev,pver-1
       rpdeli(:ncol,k) = 1._r8/(pmid(:ncol,k+1) - pmid(:ncol,k))
    end do

    !
    ! Estimate of local convective cloud cover based on convective mass flux
    ! Modify local large-scale relative humidity to account for presence of
    ! convective cloud when evaluating relative humidity based layered cloud amount
    !
    concld(:ncol,top_lev:pver) = 0.0_r8
    !
    ! cloud mass flux in SI units of kg/m2/s; should produce typical numbers of 20%
    ! shallow and deep convective cloudiness are evaluated separately (since processes
    ! are evaluated separately) and summed
    !
#ifndef PERGRO
    do k=top_lev,pver
       do i=1,ncol
          if ( .not. shallow_fraction_provided ) then
             shallowcu(i,k) = max(0.0_r8,min(sh1*log(1.0_r8+sh2*cmfmc2(i,k+1)),0.30_r8))
          else
             shallowcu(i,k) = shfrc(i,k)
          endif
          deepcu(i,k) = max(0.0_r8,min(dp1*log(1.0_r8+dp2*(cmfmc(i,k+1)-cmfmc2(i,k+1))),0.60_r8))
          concld(i,k) = min(shallowcu(i,k) + deepcu(i,k),0.80_r8)
          rh(i,k) = (rh(i,k) - concld(i,k))/(1.0_r8 - concld(i,k))
       end do
    end do
#endif
    !==================================================================================
    !
    !          ****** Compute layer cloudiness ******
    !
    !====================================================================
    ! Begin the evaluation of layered cloud amount based on (modified) RH
    !====================================================================
    !
    numkcld = pver
    do k=top_lev+1,numkcld
       kp1 = min(k + 1,pver)
       do i=1,ncol

          !++ag   This is now designed to apply FOR LIQUID CLOUDS (condensation > RH water)

          cldbnd(i) = pmid(i,k).ge.pretop

          if ( pmid(i,k).ge.premib ) then
             !==============================================================
             ! This is the low cloud (below premib) block
             !==============================================================
             ! enhance low cloud activation over land with no snow cover
             if (land(i) .and. (snowh(i) <= 0.000001_r8)) then
                rhlim = rhminl - rhminl_adj_land
             else
                rhlim = rhminl
             endif

             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             rhcloud(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)

             ! SJV: decrease cloud amount if very low water vapor content
             ! (thus very cold): "freeze dry"
             if (cldfrc_freeze_dry) then
                rhcloud(i,k) = rhcloud(i,k)*max(0.15_r8,min(1.0_r8,q(i,k)/0.0030_r8))
             endif

          else if ( pmid(i,k).lt.premit ) then
             !==============================================================
             ! This is the high cloud (above premit) block
             !==============================================================
             !
             rhlim = relhum_min(pref_mid(k), clat(i), rhminh, rhminp, unset_r8, rhminp_botmb, pi)
             !
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             rhcloud(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          else
             !==============================================================
             ! This is the middle cloud block
             !==============================================================
             !
             !       linear rh threshold transition between thresholds for low & high cloud
             !
             rhwght = (premib-(max(pmid(i,k),premit)))/(premib-premit)

             if (land(i) .and. (snowh(i) <= 0.000001_r8)) then
                rhlim = relhum_min(pref_mid(k), clat(i), rhminh, rhminp, unset_r8, rhminp_botmb, pi)*rhwght + (rhminl - rhminl_adj_land)*(1.0_r8-rhwght)
             else
                rhlim = relhum_min(pref_mid(k), clat(i), rhminh, rhminp, unset_r8, rhminp_botmb, pi)*rhwght + rhminl*(1.0_r8-rhwght)
             endif
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             rhcloud(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          end if
          !==================================================================================
          ! WE NEED TO DOCUMENT THE PURPOSE OF THIS TYPE OF CODE (ASSOCIATED WITH 2ND CALL)
          !==================================================================================
          !      !
          !      ! save rhlim to rhu00, it handles well by itself for low/high cloud
          !      !
          rhu00(i,k)=rhlim
          !==================================================================================

          if (cldfrc_ice) then

             ! Evaluate ice cloud fraction based on in-cloud ice content

             !--------ICE CLOUD OPTION 1--------Wang & Sassen 2002
             !         Evaluate desired in-cloud water content
             !               icicval = f(temp,cldice,numice)
             !         Start with a function of temperature.
             !         Wang & Sassen 2002 (JAS), based on ARM site MMCR (midlat cirrus)
             !           parameterization valid for 203-253K
             !           icival > 0 for t>195K
             if (iceopt.lt.3) then
                if (iceopt.eq.1) then
                   ttmp=max(195._r8,min(temp(i,k),253._r8)) - 273.16_r8
                   icicval=a + b * ttmp + c * ttmp**2._r8
                   !convert units
                   rho=pmid(i,k)/(rair*temp(i,k))
                   icicval= icicval * 1.e-6_r8 / rho
                else
                   !--------ICE CLOUD OPTION 2--------Schiller 2008 (JGR)
                   !          Use a curve based on FISH measurements in
                   !          tropics, mid-lats and arctic. Curve is for 180-250K (raise to 273K?)
                   !          use median all flights

                   ttmp=max(190._r8,min(temp(i,k),273.16_r8))
                   icicval = 10._r8 **(as * bs**ttmp + cs)
                   !convert units from ppmv to kg/kg
                   icicval= icicval * 1.e-6_r8 * 18._r8 / 28.97_r8
                endif
                !set icecldfraction  for OPTION 1 or OPTION2
                icecldf(i,k) =  max(0._r8,min(cldice(i,k)/icicval,1._r8))

             else if (iceopt.eq.3) then

                !--------ICE CLOUD OPTION 3--------Wood & Field 2000 (JAS)
                ! eq 6: cloud fraction = 1 - exp (-K * qc/qsati)

                icecldf(i,k)=1._r8 - exp(-Kc*cldice(i,k)/(qs(i,k)*(esi(i,k)/esl(i,k))))
                icecldf(i,k)=max(0._r8,min(icecldf(i,k),1._r8))
             else
                !--------ICE CLOUD OPTION 4--------Wilson and ballard 1999
                ! inversion of smith....
                !       ncf = cldice / ((1-RHcrit)*qs)
                ! then a function of ncf....
                ncf =cldice(i,k)/((1._r8 - icecrit)*qs(i,k))
                if (ncf.le.0._r8) then
                   icecldf(i,k)=0._r8
                else if (ncf.gt.0._r8 .and. ncf.le.1._r8/6._r8) then
                   icecldf(i,k)=0.5_r8*(6._r8 * ncf)**(2._r8/3._r8)
                else if (ncf.gt.1._r8/6._r8 .and. ncf.lt.1._r8) then
                   phi=(acos(3._r8*(1._r8-ncf)/2._r8**(3._r8/2._r8))+4._r8*3.1415927_r8)/3._r8
                   icecldf(i,k)=(1._r8 - 4._r8 * cos(phi) * cos(phi))
                else
                   icecldf(i,k)=1._r8
                endif
                icecldf(i,k)=max(0._r8,min(icecldf(i,k),1._r8))
             endif
             !TEST: if ice present, icecldf=1.
             !          if (cldice(i,k).ge.1.e-8_r8) then
             !             icecldf(i,k) = 0.99_r8
             !          endif

             !!          if ((cldice(i,k) .gt. icicval) .or. ((cldice(i,k) .gt. 0._r8) .and. (icecldf(i,k) .eq. 0._r8))) then
             !          if (cldice(i,k) .gt. 1.e-8_r8) then
             !             write(iulog,*) 'i,k,pmid,rho,t,cldice,icicval,icecldf,rhcloud: ', &
             !                i,k,pmid(i,k),rho,temp(i,k),cldice(i,k),icicval,icecldf(i,k),rhcloud(i,k)
             !          endif

             !         Combine ice and liquid cloud fraction assuming maximum overlap.
             ! Combined cloud fraction is maximum overlap
             !          cloud(i,k)=min(1._r8,max(icecldf(i,k),rhcloud(i,k)))

             liqcldf(i,k)=(1._r8 - icecldf(i,k))* rhcloud(i,k)
             cloud(i,k)=liqcldf(i,k) + icecldf(i,k)
          else
             ! For RK microphysics
             cloud(i,k) = rhcloud(i,k)
          end if
       end do
    end do
    !
    ! Add in the marine strat
    ! MARINE STRATUS SHOULD BE A SPECIAL CASE OF LAYERED CLOUD
    ! CLOUD CURRENTLY CONTAINS LAYERED CLOUD DETERMINED BY RH CRITERIA
    ! TAKE THE MAXIMUM OF THE DIAGNOSED LAYERED CLOUD OR STRATOCUMULUS
    !
    !===================================================================================
    !
    !  SOME OBSERVATIONS ABOUT THE FOLLOWING SECTION OF CODE (missed in earlier look)
    !  K700 IS SET AS A CONSTANT BASED ON HYBRID COORDINATE: IT DOES NOT DEPEND ON
    !  LOCAL PRESSURE; THERE IS NO PRESSURE RAMP => LOOKS LEVEL DEPENDENT AND
    !  DISCONTINUOUS IN SPACE (I.E., STRATUS WILL END SUDDENLY WITH NO TRANSITION)
    !
    !  IT APPEARS THAT STRAT IS EVALUATED ACCORDING TO KLEIN AND HARTMANN; HOWEVER,
    !  THE ACTUAL STRATUS AMOUNT (CLDST) APPEARS TO DEPEND DIRECTLY ON THE RH BELOW
    !  THE STRONGEST PART OF THE LOW LEVEL INVERSION.
    !PJR answers: 1) the rh limitation is a physical/mathematical limitation
    !             cant have more cloud than there is RH
    !             allowed the cloud to exist two layers below the inversion
    !             because the numerics frequently make 50% relative humidity
    !             in level below the inversion which would allow no cloud
    !             2) since  the cloud is only allowed over ocean, it should
    !             be very insensitive to surface pressure (except due to
    !             spectral ringing, which also causes so many other problems
    !             I didnt worry about it.
    !
    !==================================================================================
    if (.not.inversion_cld_off) then
    !
    ! Find most stable level below 750 mb for evaluating stratus regimes
    !
    do i=1,ncol
       ! Nothing triggers unless a stability greater than this minimum threshold is found
       dthdpmn(i) = -0.125_r8
       kdthdp(i) = 0
    end do
    !
    do k=top_lev+1,pver
       do i=1,ncol
          if (pmid(i,k) >= premib .and. ocnfrac(i).gt. 0.01_r8) then
             ! I think this is done so that dtheta/dp is in units of dg/mb (JJH)
             dthdp = 100.0_r8*(theta(i,k) - theta(i,k-1))*rpdeli(i,k-1)
             if (dthdp < dthdpmn(i)) then
                dthdpmn(i) = dthdp
                kdthdp(i) = k     ! index of interface of max inversion
             end if
          end if
       end do
    end do

    ! Also check between the bottom layer and the surface
    ! Only perform this check if the criteria were not met above

    do i = 1,ncol
       if ( kdthdp(i) .eq. 0 .and. ocnfrac(i).gt.0.01_r8) then
          dthdp = 100.0_r8 * (thetas(i) - theta(i,pver)) / (ps(i)-pmid(i,pver))
          if (dthdp < dthdpmn(i)) then
             dthdpmn(i) = dthdp
             kdthdp(i) = pver     ! index of interface of max inversion
          endif
       endif
    enddo

    do i=1,ncol
       if (kdthdp(i) /= 0) then
          k = kdthdp(i)
          kp1 = min(k+1,pver)
          ! Note: strat will be zero unless ocnfrac > 0.01
          strat = min(1._r8,max(0._r8, ocnfrac(i) * ((theta(i,k700)-thetas(i))*.057_r8-.5573_r8) ) )
          !
          ! assign the stratus to the layer just below max inversion
          ! the relative humidity changes so rapidly across the inversion
          ! that it is not safe to just look immediately below the inversion
          ! so limit the stratus cloud by rh in both layers below the inversion
          !
          cldst(i,k) = min(strat,max(rh(i,k),rh(i,kp1)))
       end if
    end do
    end if  ! .not.inversion_cld_off

    do k=top_lev,pver
       do i=1,ncol
          !
          !       which is greater; standard layered cloud amount or stratocumulus diagnosis
          !
          cloud(i,k) = max(rhcloud(i,k),cldst(i,k))
          !
          !       add in the contributions of convective cloud (determined separately and accounted
          !       for by modifications to the large-scale relative humidity.
          !
          cloud(i,k) = min(cloud(i,k)+concld(i,k), 1.0_r8)
       end do
    end do

  end subroutine cldfrc_run

  function relhum_min(press, lat, rhminh, rhminp, unset_r8, rhminp_botmb, pi) result(rh)
    real(r8), intent(in) :: press, lat, rhminh, rhminp
    real(r8), intent(in) :: unset_r8, rhminp_botmb, pi
    real(r8) :: rh

    rh = rhminh
    if (rhminp == unset_r8) return

    if (press < rhminp_botmb * 1.e2_r8 .and. &
         abs(lat * 180._r8 / pi) > 60._r8) then
       rh = rhminp
    end if
  end function relhum_min

end module ap_cldfrc_scheme
