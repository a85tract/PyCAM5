module ap_wetdepa_v2_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: wetdepa_v2_run

contains

  !> \section arg_table_wetdepa_v2_run Argument Table
  !! \htmlinclude wetdepa_v2_run.html
  !!
subroutine wetdepa_v2_run(                                  &
   pcols, pver, gravit, iulog, p, q, pdel, cldt, cldc,  &
   cmfdqr, evapc, conicw, precs, conds,                 &
   evaps, cwat, tracer, deltat, scavt,                  &
   iscavt, cldvcu, cldvst, dlf, fracis,                 &
   sol_fact, ncol, scavcoef, is_strat_cloudborne, qqcw, &
   f_act_conv, icscavt, isscavt, bcscavt, bsscavt,      &
   sol_facti_in, sol_factic_in, negative_dblchek,         &
   negative_srct, negative_rat, negative_fracev )

   !-----------------------------------------------------------------------
   !
   ! scavenging code for very soluble aerosols
   !
   !-----------------------------------------------------------------------

   integer, intent(in) :: pcols, pver, iulog
   real(r8), intent(in) :: gravit

   real(r8), intent(in) ::&
      p(pcols,pver),        &! pressure
      q(pcols,pver),        &! moisture
      pdel(pcols,pver),     &! pressure thikness
      cldt(pcols,pver),     &! total cloud fraction
      cldc(pcols,pver),     &! convective cloud fraction
      cmfdqr(pcols,pver),   &! rate of production of convective precip
      evapc(pcols,pver),    &! Evaporation rate of convective precipitation
      conicw(pcols,pver),   &! convective cloud water
      cwat(pcols,pver),     &! cloud water amount
      precs(pcols,pver),    &! rate of production of stratiform precip
      conds(pcols,pver),    &! rate of production of condensate
      evaps(pcols,pver),    &! rate of evaporation of precip
      cldvcu(pcols,pver),   &! Convective precipitation area at the top interface of each layer
      cldvst(pcols,pver),   &! Stratiform precipitation area at the top interface of each layer
      dlf(pcols,pver),      &! Detrainment of convective condensate [kg/kg/s]
      deltat,               &! time step
      tracer(pcols,pver)     ! trace species

   ! If subroutine is called with just sol_fact:
   !    sol_fact is used for both in- and below-cloud scavenging
   ! If subroutine is called with optional argument sol_facti_in:
   !    sol_fact  is used for below cloud scavenging
   !    sol_facti is used for in cloud scavenging

   real(r8), intent(in)  :: sol_fact
   integer,  intent(in)  :: ncol
   real(r8), intent(in)  :: scavcoef(pcols,pver) ! Dana and Hales coefficient (/mm) (0.1 if not MODAL_AERO)
   real(r8), intent(out) ::&
      scavt(pcols,pver),   &! scavenging tend
      iscavt(pcols,pver),  &! incloud scavenging tends
      fracis(pcols,pver)    ! fraction of species not scavenged

   ! Setting is_strat_cloudborne=.true. indicates that tracer is stratiform-cloudborne aerosol.
   !   This is only used by MAM code.  The optional args qqcw and f_act_conv are not referenced
   !   in this case.
   ! Setting is_strat_cloudborne=.false. is being used to indicate that the tracers are the
   !   interstitial modal aerosols.  In this case the optional qqcw (the cloud borne mixing ratio
   !   corresponding to the interstitial aerosol) must be provided, as well as the optional f_act_conv.
   logical,  intent(in), optional :: is_strat_cloudborne
   real(r8), intent(in), optional :: qqcw(pcols,pver)
   real(r8), intent(in), optional :: f_act_conv(pcols,pver)

   real(r8), intent(in), optional :: sol_facti_in   ! solubility factor (frac of aerosol scavenged in cloud)
   real(r8), intent(in), optional :: sol_factic_in(pcols,pver)  ! sol_facti_in for convective clouds


   real(r8), intent(out), optional :: icscavt(pcols,pver)     ! incloud, convective
   real(r8), intent(out), optional :: isscavt(pcols,pver)     ! incloud, stratiform
   real(r8), intent(out), optional :: bcscavt(pcols,pver)     ! below cloud, convective
   real(r8), intent(out), optional :: bsscavt(pcols,pver)     ! below cloud, stratiform

   ! Diagnostic values are returned to the CAM host facade so the standalone
   ! science kernel does not depend on CAM logging infrastructure.
   real(r8), intent(out) :: negative_dblchek(pcols,pver)
   real(r8), intent(out) :: negative_srct(pcols,pver)
   real(r8), intent(out) :: negative_rat(pcols,pver)
   real(r8), intent(out) :: negative_fracev(pcols,pver)

   ! local variables

   integer :: i, k

   real(r8) :: omsm                 ! 1 - (a small number)
   real(r8) :: clds(pcols)          ! stratiform cloud fraction
   real(r8) :: fracev(pcols)        ! fraction of precip from above that is evaporating
   real(r8) :: fracev_cu(pcols)     ! Fraction of convective precip from above that is evaporating
   real(r8) :: fracp(pcols)         ! fraction of cloud water converted to precip
   real(r8) :: pdog(pcols)          ! work variable (pdel/gravit)
   real(r8) :: rpdog(pcols)         ! work variable (gravit/pdel)
   real(r8) :: precabc(pcols)       ! conv precip from above (work array)
   real(r8) :: precabs(pcols)       ! strat precip from above (work array)
   real(r8) :: rat(pcols)           ! ratio of amount available to amount removed
   real(r8) :: scavab(pcols)        ! scavenged tracer flux from above (work array)
   real(r8) :: scavabc(pcols)       ! scavenged tracer flux from above (work array)
   real(r8) :: srcc(pcols)          ! tend for convective rain
   real(r8) :: srcs(pcols)          ! tend for stratiform rain
   real(r8) :: srct(pcols)          ! work variable

   real(r8) :: fins(pcols)          ! fraction of rem. rate by strat rain
   real(r8) :: finc(pcols)          ! fraction of rem. rate by conv. rain
   real(r8) :: conv_scav_ic(pcols)  ! convective scavenging incloud
   real(r8) :: conv_scav_bc(pcols)  ! convective scavenging below cloud
   real(r8) :: st_scav_ic(pcols)    ! stratiform scavenging incloud
   real(r8) :: st_scav_bc(pcols)    ! stratiform scavenging below cloud

   real(r8) :: odds(pcols)          ! limit on removal rate (proportional to prec)
   real(r8) :: dblchek(pcols)
   real(r8) :: trac_qqcw(pcols)
   real(r8) :: tracer_incu(pcols)
   real(r8) :: tracer_mean(pcols)

   ! For stratiform cloud, cloudborne aerosol is treated explicitly,
   !    and sol_facti is 1.0 for cloudborne, 0.0 for interstitial.
   ! For convective cloud, cloudborne aerosol is not treated explicitly,
   !    and sol_factic is 1.0 for both cloudborne and interstitial.

   real(r8) :: sol_facti              ! in cloud fraction of aerosol scavenged
   real(r8) :: sol_factb              ! below cloud fraction of aerosol scavenged
   real(r8) :: sol_factic(pcols,pver) ! in cloud fraction of aerosol scavenged for convective clouds

   real(r8) :: rdeltat
   ! ------------------------------------------------------------------------

   omsm = 1._r8-2*epsilon(1._r8) ! used to prevent roundoff errors below zero

   ! default (if other sol_facts aren't in call, set all to required sol_fact)
   sol_facti = sol_fact
   sol_factb = sol_fact

   if ( present(sol_facti_in) )  sol_facti = sol_facti_in

   sol_factic  = sol_facti
   if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in

   ! this section of code is for highly soluble aerosols,
   ! the assumption is that within the cloud that
   ! all the tracer is in the cloud water
   !
   ! for both convective and stratiform clouds,
   ! the fraction of cloud water converted to precip defines
   ! the amount of tracer which is pulled out.

   precabs(:ncol) = 0.0_r8
   precabc(:ncol) = 0.0_r8
   scavab(:ncol)  = 0.0_r8
   scavabc(:ncol) = 0.0_r8

   do k = 1, pver
      do i = 1, ncol

         clds(i)  = cldt(i,k) - cldc(i,k)
         pdog(i)  = pdel(i,k)/gravity
         rpdog(i) = gravity/pdel(i,k)
         rdeltat  = 1.0_r8/deltat

         ! ****************** Evaporation **************************
         ! calculate the fraction of strat precip from above
         !                 which evaporates within this layer
         fracev(i) = evaps(i,k)*pdog(i) &
                     /max(1.e-12_r8,precabs(i))

         ! trap to ensure reasonable ratio bounds
         fracev(i) = max(0._r8,min(1._r8,fracev(i)))

         ! Same as above but convective precipitation part
         fracev_cu(i) = evapc(i,k)*pdog(i)/max(1.e-12_r8,precabc(i))
         fracev_cu(i) = max(0._r8,min(1._r8,fracev_cu(i)))

         ! ****************** Convection ***************************
         !
         ! set odds proportional to fraction of the grid box that is swept by the
         ! precipitation =precabc/rhoh20*(area of sphere projected on plane
         !                                /volume of sphere)*deltat
         ! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
         ! unless the fraction of the area that is cloud is less than odds, in which
         ! case use the cloud fraction (assumes precabs is in kg/m2/s)
         ! is really: precabs*3/4/1000./1e-3*deltat
         ! here I use .1 from Balkanski
         !
         ! use a local rate of convective rain production for incloud scav
         !
         ! Fraction of convective cloud water converted to rain.  This version is used
         ! in 2 of the 3 branches below before fracp is reused in the stratiform calc.
         ! NB: In below formula for fracp conicw is a LWC/IWC that has already
         !     precipitated out, i.e., conicw does not contain precipitation

         fracp(i) = cmfdqr(i,k)*deltat / &
                    max( 1.e-12_r8, cldc(i,k)*conicw(i,k) + (cmfdqr(i,k)+dlf(i,k))*deltat )
         fracp(i) = max( min( 1._r8, fracp(i)), 0._r8 )

         if ( present(is_strat_cloudborne) ) then

            if ( is_strat_cloudborne ) then

               ! convective scavenging

               conv_scav_ic(i) = 0._r8

               conv_scav_bc(i) = 0._r8

               ! stratiform scavenging

               fracp(i) = precs(i,k)*deltat / &
                          max( 1.e-12_r8, cwat(i,k) + precs(i,k)*deltat )
               fracp(i) = max( 0._r8, min(1._r8, fracp(i)) )
               st_scav_ic(i) = sol_facti *fracp(i)*tracer(i,k)*rdeltat

               st_scav_bc(i) = 0._r8

            else

               ! convective scavenging

               trac_qqcw(i) = min(qqcw(i,k), &
                                  tracer(i,k)*( clds(i)/max( 0.01_r8, 1._r8-clds(i) ) ) )

               tracer_incu(i) = f_act_conv(i,k)*(tracer(i,k) + trac_qqcw(i))

               conv_scav_ic(i) = sol_factic(i,k)*cldc(i,k)*fracp(i)*tracer_incu(i)*rdeltat

               tracer_mean(i) = tracer(i,k)*(1._r8 - cldc(i,k)*f_act_conv(i,k)) - &
                                cldc(i,k)*f_act_conv(i,k)*trac_qqcw(i)
               tracer_mean(i) = max(0._r8,tracer_mean(i))

               odds(i) = precabc(i)/max(cldvcu(i,k),1.e-5_r8)*scavcoef(i,k)*deltat
               odds(i) = max(min(1._r8,odds(i)),0._r8)
               conv_scav_bc(i) = sol_factb *cldvcu(i,k)*odds(i)*tracer_mean(i)*rdeltat


               ! stratiform scavenging

               st_scav_ic(i) = 0._r8

               odds(i) = precabs(i)/max(cldvst(i,k),1.e-5_r8)*scavcoef(i,k)*deltat
               odds(i) = max(min(1._r8,odds(i)),0._r8)
               st_scav_bc(i) = sol_factb *cldvst(i,k)*odds(i)*tracer_mean(i)*rdeltat

            end if

         else

            ! convective scavenging

            conv_scav_ic(i) = sol_factic(i,k)*cldc(i,k)*fracp(i)*tracer(i,k)*rdeltat

            odds(i) = precabc(i)/max(cldvcu(i,k), 1.e-5_r8)*scavcoef(i,k)*deltat
            odds(i) = max( min(1._r8, odds(i)), 0._r8)
            conv_scav_bc(i) = sol_factb*cldvcu(i,k)*odds(i)*tracer(i,k)*rdeltat

            ! stratiform scavenging

            ! fracp is the fraction of cloud water converted to precip
            ! NB: In below formula for fracp cwat is a LWC/IWC that has already
            !     precipitated out, i.e., cwat does not contain precipitation
            fracp(i) = precs(i,k)*deltat / &
                       max( 1.e-12_r8, cwat(i,k) + precs(i,k)*deltat )
            fracp(i) = max( 0._r8, min( 1._r8, fracp(i) ) )

            ! assume the corresponding amnt of tracer is removed
            st_scav_ic(i) = sol_facti*clds(i)*fracp(i)*tracer(i,k)*rdeltat

            odds(i) = precabs(i)/max(cldvst(i,k),1.e-5_r8)*scavcoef(i,k)*deltat
            odds(i) = max(min(1._r8,odds(i)),0._r8)
            st_scav_bc(i) =sol_factb*(cldvst(i,k)*odds(i)) *tracer(i,k)*rdeltat

         end if

         ! total convective scavenging
         srcc(i) = conv_scav_ic(i) + conv_scav_bc(i)
         finc(i) = conv_scav_ic(i)/(srcc(i) + 1.e-36_r8)

         ! total stratiform scavenging
         srcs(i) = st_scav_ic(i) + st_scav_bc(i)
         fins(i) = st_scav_ic(i)/(srcs(i) + 1.e-36_r8)

         ! make sure we dont take out more than is there
         ! ratio of amount available to amount removed
         rat(i) = tracer(i,k)/max(deltat*(srcc(i)+srcs(i)),1.e-36_r8)
         if (rat(i).lt.1._r8) then
            srcs(i) = srcs(i)*rat(i)
            srcc(i) = srcc(i)*rat(i)
         endif
         srct(i) = (srcc(i)+srcs(i))*omsm


         ! fraction that is not removed within the cloud
         ! (assumed to be interstitial, and subject to convective transport)
         fracp(i) = deltat*srct(i)/max(cldvst(i,k)*tracer(i,k),1.e-36_r8)  ! amount removed
         fracp(i) = max(0._r8,min(1._r8,fracp(i)))
         fracis(i,k) = 1._r8 - fracp(i)

         ! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
         ! Sungsu added cumulus contribution in the below 3 blocks
         scavt(i,k) = -srct(i) + (fracev(i)*scavab(i)+fracev_cu(i)*scavabc(i))*rpdog(i)
         iscavt(i,k) = -(srcc(i)*finc(i) + srcs(i)*fins(i))*omsm

         if ( present(icscavt) ) icscavt(i,k) = -(srcc(i)*finc(i)) * omsm
         if ( present(isscavt) ) isscavt(i,k) = -(srcs(i)*fins(i)) * omsm
         if ( present(bcscavt) ) bcscavt(i,k) = -(srcc(i) * (1-finc(i))) * omsm +  &
            fracev_cu(i)*scavabc(i)*rpdog(i)

         if ( present(bsscavt) ) bsscavt(i,k) = -(srcs(i) * (1-fins(i))) * omsm +  &
            fracev(i)*scavab(i)*rpdog(i)
         dblchek(i) = tracer(i,k) + deltat*scavt(i,k)
         negative_dblchek(i,k) = dblchek(i)
         negative_srct(i,k) = srct(i)
         negative_rat(i,k) = rat(i)
         negative_fracev(i,k) = fracev(i)

         ! now keep track of scavenged mass and precip
         scavab(i) = scavab(i)*(1-fracev(i)) + srcs(i)*pdog(i)
         precabs(i) = precabs(i) + (precs(i,k) - evaps(i,k))*pdog(i)
         scavabc(i) = scavabc(i)*(1-fracev_cu(i)) + srcc(i)*pdog(i)
         precabc(i) = precabc(i) + (cmfdqr(i,k) - evapc(i,k))*pdog(i)

      end do ! End of i = 1, ncol

   end do ! End of k = 1, pver

end subroutine wetdepa_v2_run

end module ap_wetdepa_v2_scheme
