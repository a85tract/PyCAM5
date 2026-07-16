module ap_aero_model_wetdep_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ap_wetdepa_v2_scheme, only : wetdepa_v2_run
  use aero_model_wetdep_host_hooks, only : &
       aero_model_wetdep_timer_start, aero_model_wetdep_timer_stop, &
       aero_model_wetdep_outfld_real1d, aero_model_wetdep_outfld_real2d

  implicit none
  private
  save

  public :: aero_model_wetdep_init, aero_model_wetdep_run

  integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12

  integer :: pcols, pver, pcnst, ntot_amode, iulog
  integer :: modeptr_coarse, modeptr_pcarbon, modeptr_finedust
  integer :: modeptr_coardust
  real(r8) :: gravit
  integer, allocatable :: nspec_amode(:), numptr_amode(:)
  integer, allocatable :: numptrcw_amode(:), lmassptr_amode(:,:)
  integer, allocatable :: lmassptrcw_amode(:,:), lspectype_amode(:,:)
  integer, allocatable :: lptr_dust_a_amode(:), lptr_nacl_a_amode(:)
  real(r8), allocatable :: dgnum_amode(:), spechygro(:), specdens_amode(:)

contains

  subroutine aero_model_wetdep_init(pcols_in, pver_in, pcnst_in, &
       gravit_in, iulog_in, nspec_amode_in, numptr_amode_in, &
       numptrcw_amode_in, lmassptr_amode_in, lmassptrcw_amode_in, &
       lspectype_amode_in, modeptr_coarse_in, modeptr_pcarbon_in, &
       modeptr_finedust_in, modeptr_coardust_in, lptr_dust_a_amode_in, &
       lptr_nacl_a_amode_in, dgnum_amode_in, spechygro_in, &
       specdens_amode_in)

    integer, intent(in) :: pcols_in, pver_in, pcnst_in, iulog_in
    real(r8), intent(in) :: gravit_in
    integer, intent(in) :: nspec_amode_in(:), numptr_amode_in(:)
    integer, intent(in) :: numptrcw_amode_in(:)
    integer, intent(in) :: lmassptr_amode_in(:,:), lmassptrcw_amode_in(:,:)
    integer, intent(in) :: lspectype_amode_in(:,:)
    integer, intent(in) :: modeptr_coarse_in, modeptr_pcarbon_in
    integer, intent(in) :: modeptr_finedust_in, modeptr_coardust_in
    integer, intent(in) :: lptr_dust_a_amode_in(:), lptr_nacl_a_amode_in(:)
    real(r8), intent(in) :: dgnum_amode_in(:), spechygro_in(:)
    real(r8), intent(in) :: specdens_amode_in(:)

    pcols = pcols_in
    pver = pver_in
    pcnst = pcnst_in
    gravit = gravit_in
    iulog = iulog_in
    ntot_amode = size(nspec_amode_in)
    modeptr_coarse = modeptr_coarse_in
    modeptr_pcarbon = modeptr_pcarbon_in
    modeptr_finedust = modeptr_finedust_in
    modeptr_coardust = modeptr_coardust_in

    allocate(nspec_amode(ntot_amode), numptr_amode(ntot_amode), &
         numptrcw_amode(ntot_amode), &
         lmassptr_amode(size(lmassptr_amode_in,1),ntot_amode), &
         lmassptrcw_amode(size(lmassptrcw_amode_in,1),ntot_amode), &
         lspectype_amode(size(lspectype_amode_in,1),ntot_amode), &
         lptr_dust_a_amode(ntot_amode), lptr_nacl_a_amode(ntot_amode), &
         dgnum_amode(ntot_amode), spechygro(size(spechygro_in)), &
         specdens_amode(size(specdens_amode_in)))

    nspec_amode = nspec_amode_in
    numptr_amode = numptr_amode_in
    numptrcw_amode = numptrcw_amode_in
    lmassptr_amode = lmassptr_amode_in
    lmassptrcw_amode = lmassptrcw_amode_in
    lspectype_amode = lspectype_amode_in
    lptr_dust_a_amode = lptr_dust_a_amode_in
    lptr_nacl_a_amode = lptr_nacl_a_amode_in
    dgnum_amode = dgnum_amode_in
    spechygro = spechygro_in
    specdens_amode = specdens_amode_in
  end subroutine aero_model_wetdep_init

  !> \section arg_table_aero_model_wetdep_run Argument Table
  !! \htmlinclude aero_model_wetdep_run.html
  !!
  subroutine aero_model_wetdep_run(lchnk, ncol, dt, dlf, q, pmid, pdel, &
       ptend_q, ptend_lq, cldt, cldcu, cmfdqr, evapc, conicw, prain, &
       qme, evapr, totcond, cldvcu, cldvst, dgnumwet, fracis, qqcw_data, &
       sol_facti_cloud_borne_in, sol_factb_interstitial_in, &
       sol_factic_interstitial_in, &
       dlndg_nimptblgrow_in, scavimptblnum_in, scavimptblvol_in, &
       aerdepwetis, aerdepwetcw)

    ! args
    integer, intent(in) :: lchnk, ncol
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: dlf(:,:), q(:,:,:), pmid(:,:), pdel(:,:)
    real(r8), intent(inout) :: ptend_q(:,:,:), fracis(:,:,:)
    logical, intent(inout) :: ptend_lq(:)
    real(r8), intent(in) :: cldt(:,:), cldcu(:,:), cmfdqr(:,:), evapc(:,:)
    real(r8), intent(in) :: conicw(:,:), prain(:,:), qme(:,:), evapr(:,:)
    real(r8), intent(in) :: totcond(:,:), cldvcu(:,:), cldvst(:,:)
    real(r8), intent(in) :: dgnumwet(:,:,:)
    real(r8), intent(inout) :: qqcw_data(:,:,:)
    real(r8), intent(in) :: sol_facti_cloud_borne_in
    real(r8), intent(in) :: sol_factb_interstitial_in
    real(r8), intent(in) :: sol_factic_interstitial_in
    real(r8), intent(in) :: dlndg_nimptblgrow_in
    real(r8), intent(in) :: scavimptblnum_in(-7:,:)
    real(r8), intent(in) :: scavimptblvol_in(-7:,:)
    real(r8), intent(out) :: aerdepwetis(:,:), aerdepwetcw(:,:)

    ! local vars

    integer :: m ! tracer index

    real(r8) :: iscavt(pcols, pver)

    integer :: mm
    integer :: i,k

    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)
    real(r8) :: sol_factb, sol_facti
    real(r8) :: sol_factic(pcols,pver)

    real(r8) :: sflx(pcols) ! deposition flux

    integer :: jnv ! index for scavcoefnv 3rd dimension
    integer :: lphase ! index for interstitial / cloudborne aerosol
    integer :: lspec ! index for aerosol number / chem-mass / water-mass
    integer :: lcoardust, lcoarnacl ! indices for coarse mode dust and seasalt masses
    real(r8) :: dqdt_tmp(pcols,pver) ! temporary array to hold tendency for 1 species
    real(r8) :: f_act_conv(pcols,pver) ! prescribed aerosol activation fraction for convective cloud ! rce 2010/05/01
    real(r8) :: f_act_conv_coarse(pcols,pver) ! similar but for coarse mode ! rce 2010/05/02
    real(r8) :: f_act_conv_coarse_dust, f_act_conv_coarse_nacl ! rce 2010/05/02
    real(r8) :: fracis_cw(pcols,pver)
    real(r8) :: hygro_sum_old(pcols,pver) ! before removal [sum of (mass*hydro/dens)]
    real(r8) :: hygro_sum_del(pcols,pver) ! removal change to [sum of (mass*hydro/dens)]
    real(r8) :: hygro_sum_old_ik, hygro_sum_new_ik
    real(r8) :: prec(pcols) ! precipitation rate
    real(r8) :: q_tmp(pcols,pver) ! temporary array to hold "most current" mixing ratio for 1 species
    real(r8) :: scavcoefnv(pcols,pver,0:2) ! Dana and Hales coefficient (/mm) for
                                           ! cloud-borne num & vol (0),
                                           ! interstitial num (1), interstitial vol (2)
    real(r8) :: tmpa, tmpb
    real(r8) :: tmpdust, tmpnacl
    ! The standalone wetdepa kernel returns the values that the legacy CAM
    ! facade uses for negative-value diagnostics.  This caller does not emit
    ! that diagnostic, but it must still provide storage for those outputs.
    real(r8) :: negative_dblchek(pcols,pver)
    real(r8) :: negative_srct(pcols,pver)
    real(r8) :: negative_rat(pcols,pver)
    real(r8) :: negative_fracev(pcols,pver)
    logical  :: isprx(pcols,pver) ! true if precipation

    prec(:ncol)=0._r8
    do k=1,pver
       where (prec(:ncol) >= 1.e-7_r8)
          isprx(:ncol,k) = .true.
       elsewhere
          isprx(:ncol,k) = .false.
       endwhere
       prec(:ncol) = prec(:ncol) + (prain(:ncol,k) + cmfdqr(:ncol,k) - &
            evapr(:ncol,k))*pdel(:ncol,k)/gravit
    end do

    ! calculate the mass-weighted sol_factic for coarse mode species
    ! sol_factic_coarse(:,:) = 0.30_r8 ! tuned 1/4
    f_act_conv_coarse(:,:) = 0.60_r8 ! rce 2010/05/02
    f_act_conv_coarse_dust = 0.40_r8 ! rce 2010/05/02
    f_act_conv_coarse_nacl = 0.80_r8 ! rce 2010/05/02
    if (modeptr_coarse > 0) then
       lcoardust = lptr_dust_a_amode(modeptr_coarse)
       lcoarnacl = lptr_nacl_a_amode(modeptr_coarse)
       if ((lcoardust > 0) .and. (lcoarnacl > 0)) then
          do k = 1, pver
             do i = 1, ncol
                tmpdust = max(0.0_r8, q(i,k,lcoardust) + &
                     ptend_q(i,k,lcoardust)*dt)
                tmpnacl = max(0.0_r8, q(i,k,lcoarnacl) + &
                     ptend_q(i,k,lcoarnacl)*dt)
                if ((tmpdust+tmpnacl) > 1.0e-30_r8) then
                   ! sol_factic_coarse(i,k) = (0.2_r8*tmpdust + 0.4_r8*tmpnacl)/(tmpdust+tmpnacl) ! tuned 1/6
                   f_act_conv_coarse(i,k) = (f_act_conv_coarse_dust*tmpdust &
                        + f_act_conv_coarse_nacl*tmpnacl)/(tmpdust+tmpnacl) ! rce 2010/05/02
                end if
             end do
          end do
       end if
    end if

    scavcoefnv(:,:,0) = 0.0_r8 ! below-cloud scavcoef = 0.0 for cloud-borne species

    do m = 1, ntot_amode ! main loop over aerosol modes

       do lphase = 1, 2 ! loop over interstitial (1) and cloud-borne (2) forms

          ! sol_factb and sol_facti values
          ! sol_factb - currently this is basically a tuning factor
          ! sol_facti & sol_factic - currently has a physical basis, and reflects activation fraction
          !
          ! 2008-mar-07 rce - sol_factb (interstitial) changed from 0.3 to 0.1
          ! - sol_factic (interstitial, dust modes) changed from 1.0 to 0.5
          ! - sol_factic (cloud-borne, pcarb modes) no need to set it to 0.0
          ! because the cloud-borne pcarbon == 0 (no activation)
          !
          ! rce 2010/05/02
          ! prior to this date, sol_factic was used for convective in-cloud wet removal,
          ! and its value reflected a combination of an activation fraction (which varied between modes)
          ! and a tuning factor
          ! from this date forward, two parameters are used for convective in-cloud wet removal
          ! f_act_conv is the activation fraction
          ! note that "non-activation" of aerosol in air entrained into updrafts should
          ! be included here
          ! eventually we might use the activate routine (with w ~= 1 m/s) to calculate
          ! this, but there is still the entrainment issue
          ! sol_factic is strictly a tuning factor
          !
          if (lphase == 1) then ! interstial aerosol
             hygro_sum_old(:,:) = 0.0_r8
             hygro_sum_del(:,:) = 0.0_r8
             call modal_aero_bcscavcoef_get( m, ncol, isprx, dgnumwet, &
                  scavcoefnv(:,:,1), scavcoefnv(:,:,2), dlndg_nimptblgrow_in, &
                  scavimptblnum_in, scavimptblvol_in )

             sol_factb = sol_factb_interstitial_in ! all below-cloud scav ON (0.1 "tuning factor")

             sol_facti = 0.0_r8 ! strat in-cloud scav totally OFF for institial

             sol_factic = sol_factic_interstitial_in

             if (m == modeptr_pcarbon) then
                ! sol_factic = 0.0_r8 ! conv in-cloud scav OFF (0.0 activation fraction)
                f_act_conv = 0.0_r8 ! rce 2010/05/02
             else if ((m == modeptr_finedust) .or. (m == modeptr_coardust)) then
                ! sol_factic = 0.2_r8 ! conv in-cloud scav ON (0.5 activation fraction) ! tuned 1/4
                f_act_conv = 0.4_r8 ! rce 2010/05/02
             else
                ! sol_factic = 0.4_r8 ! conv in-cloud scav ON (1.0 activation fraction) ! tuned 1/4
                f_act_conv = 0.8_r8 ! rce 2010/05/02
             end if

          else ! cloud-borne aerosol (borne by stratiform cloud drops)

             sol_factb  = 0.0_r8   ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
             sol_facti  = sol_facti_cloud_borne_in   ! strat  in-cloud scav cloud-borne tuning factor
             sol_factic = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean
                                   !        that conv precip collects strat droplets)
             f_act_conv = 0.0_r8   ! conv   in-cloud scav OFF (having this on would mean

          end if
          !
          ! rce 2010/05/03
          ! wetdepa has "sol_fact" parameters:
          ! sol_facti, sol_factic, sol_factb for liquid cloud

          do lspec = 0, nspec_amode(m)+1 ! loop over number + chem constituents + water

             if (lspec == 0) then ! number
                if (lphase == 1) then
                   mm = numptr_amode(m)
                   jnv = 1
                else
                   mm = numptrcw_amode(m)
                   jnv = 0
                endif
             else if (lspec <= nspec_amode(m)) then ! non-water mass
                if (lphase == 1) then
                   mm = lmassptr_amode(lspec,m)
                   jnv = 2
                else
                   mm = lmassptrcw_amode(lspec,m)
                   jnv = 0
                endif
             else ! water mass
                ! bypass wet removal of aerosol water
                cycle
                if (lphase == 1) then
                   mm = 0
                   ! mm = lwaterptr_amode(m)
                   jnv = 2
                else
                   mm = 0
                   jnv = 0
                endif
             endif

             if (mm <= 0) cycle


             ! set f_act_conv for interstitial (lphase=1) coarse mode species
             ! for the convective in-cloud, we conceptually treat the coarse dust and seasalt
             ! as being externally mixed, and apply f_act_conv = f_act_conv_coarse_dust/nacl to dust/seasalt
             ! number and sulfate are conceptually partitioned to the dust and seasalt
             ! on a mass basis, so the f_act_conv for number and sulfate are
             ! mass-weighted averages of the values used for dust/seasalt
             if ((lphase == 1) .and. (m == modeptr_coarse)) then
                ! sol_factic = sol_factic_coarse
                f_act_conv = f_act_conv_coarse ! rce 2010/05/02
                if (lspec > 0) then
                   if (lmassptr_amode(lspec,m) == lptr_dust_a_amode(m)) then
                      ! sol_factic = 0.2_r8 ! tuned 1/4
                      f_act_conv = f_act_conv_coarse_dust ! rce 2010/05/02
                   else if (lmassptr_amode(lspec,m) == lptr_nacl_a_amode(m)) then
                      ! sol_factic = 0.4_r8 ! tuned 1/6
                      f_act_conv = f_act_conv_coarse_nacl ! rce 2010/05/02
                   end if
                end if
             end if


             if ((lphase == 1) .and. (lspec <= nspec_amode(m))) then
                ptend_lq(mm) = .true.
                dqdt_tmp(:,:) = 0.0_r8
                ! q_tmp reflects changes from modal_aero_calcsize and is the "most current" q
                q_tmp(1:ncol,:) = q(1:ncol,:,mm) + ptend_q(1:ncol,:,mm)*dt

                call aero_model_wetdep_timer_start('ap_wetdepa_v2_run')
                call wetdepa_v2_run(pcols, pver, gravit, iulog, pmid, q(:,:,1), &
                     pdel, cldt, cldcu, cmfdqr, evapc, conicw, prain, qme, &
                     evapr, totcond, q_tmp, dt, dqdt_tmp, iscavt, cldvcu, cldvst, &
                     dlf, fracis(:,:,mm), sol_factb, ncol, &
                     scavcoefnv(:,:,jnv), &
                     is_strat_cloudborne=.false.,  &
                     qqcw=qqcw_data(:,:,mm),  &
                     f_act_conv=f_act_conv, &
                     icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                     sol_facti_in=sol_facti, sol_factic_in=sol_factic, &
                     negative_dblchek=negative_dblchek, negative_srct=negative_srct, &
                     negative_rat=negative_rat, negative_fracev=negative_fracev)
                call aero_model_wetdep_timer_stop('ap_wetdepa_v2_run')

                ptend_q(1:ncol,:,mm) = ptend_q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

                call aero_model_wetdep_outfld_real2d(mm, .false., 'WET', &
                     dqdt_tmp, pcols, lchnk)
                call aero_model_wetdep_outfld_real2d(mm, .false., 'SIC', &
                     icscavt, pcols, lchnk)
                call aero_model_wetdep_outfld_real2d(mm, .false., 'SIS', &
                     isscavt, pcols, lchnk)
                call aero_model_wetdep_outfld_real2d(mm, .false., 'SBC', &
                     bcscavt, pcols, lchnk)
                call aero_model_wetdep_outfld_real2d(mm, .false., 'SBS', &
                     bsscavt, pcols, lchnk)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+dqdt_tmp(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .false., 'SFWET', &
                     sflx, pcols, lchnk)
                aerdepwetis(:ncol,mm) = sflx(:ncol)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+icscavt(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .false., 'SFSIC', &
                     sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+isscavt(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .false., 'SFSIS', &
                     sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bcscavt(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .false., 'SFSBC', &
                     sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bsscavt(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .false., 'SFSBS', &
                     sflx, pcols, lchnk)

                if (lspec > 0) then
                   tmpa = spechygro(lspectype_amode(lspec,m))/ &
                        specdens_amode(lspectype_amode(lspec,m))
                   tmpb = tmpa*dt
                   hygro_sum_old(1:ncol,:) = hygro_sum_old(1:ncol,:) &
                        + tmpa*q_tmp(1:ncol,:)
                   hygro_sum_del(1:ncol,:) = hygro_sum_del(1:ncol,:) &
                        + tmpb*dqdt_tmp(1:ncol,:)
                end if

             else ! lphase == 2
                dqdt_tmp(:,:) = 0.0_r8

                call aero_model_wetdep_timer_start('ap_wetdepa_v2_run')
                call wetdepa_v2_run(pcols, pver, gravit, iulog, pmid, q(:,:,1), &
                     pdel, cldt, cldcu, cmfdqr, evapc, conicw, prain, qme, &
                     evapr, totcond, qqcw_data(:,:,mm), dt, &
                     dqdt_tmp, iscavt, cldvcu, cldvst, &
                     dlf, fracis_cw, sol_factb, ncol, &
                     scavcoefnv(:,:,jnv), &
                     is_strat_cloudborne=.true.,  &
                     icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                     sol_facti_in=sol_facti, sol_factic_in=sol_factic, &
                     negative_dblchek=negative_dblchek, negative_srct=negative_srct, &
                     negative_rat=negative_rat, negative_fracev=negative_fracev)
                call aero_model_wetdep_timer_stop('ap_wetdepa_v2_run')

                qqcw_data(1:ncol,:,mm) = qqcw_data(1:ncol,:,mm) + &
                     dqdt_tmp(1:ncol,:) * dt

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+dqdt_tmp(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .true., 'SFWET', &
                     sflx, pcols, lchnk)
                aerdepwetcw(:ncol,mm) = sflx(:ncol)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+icscavt(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .true., 'SFSIC', &
                     sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+isscavt(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .true., 'SFSIS', &
                     sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bcscavt(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .true., 'SFSBC', &
                     sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bsscavt(i,k)*pdel(i,k)/gravit
                   enddo
                enddo
                call aero_model_wetdep_outfld_real1d(mm, .true., 'SFSBS', &
                     sflx, pcols, lchnk)

             endif

          enddo ! lspec = 0, nspec_amode(m)+1
       enddo ! lphase = 1, 2
    enddo ! m = 1, ntot_amode

  end subroutine aero_model_wetdep_run

  subroutine modal_aero_bcscavcoef_get( m, ncol, isprx, dgn_awet, scavcoefnum, scavcoefvol, &
       dlndg_nimptblgrow_in, scavimptblnum_in, scavimptblvol_in )

    !-----------------------------------------------------------------------
    implicit none

    integer,intent(in) :: m, ncol
    logical,intent(in):: isprx(pcols,pver)
    real(r8), intent(in) :: dgn_awet(pcols,pver,ntot_amode)
    real(r8), intent(out) :: scavcoefnum(pcols,pver), scavcoefvol(pcols,pver)
    real(r8), intent(in) :: dlndg_nimptblgrow_in
    real(r8), intent(in) :: scavimptblnum_in(-7:,:)
    real(r8), intent(in) :: scavimptblvol_in(-7:,:)

    integer i, k, jgrow
    real(r8) dumdgratio, xgrow, dumfhi, dumflo, scavimpvol, scavimpnum


    do k = 1, pver
       do i = 1, ncol

          ! do only if no precip
          if ( isprx(i,k) ) then
             !
             ! interpolate table values using log of (actual-wet-size)/(base-dry-size)

             dumdgratio = dgn_awet(i,k,m)/dgnum_amode(m)

             if ((dumdgratio .ge. 0.99_r8) .and. (dumdgratio .le. 1.01_r8)) then
                scavimpvol = scavimptblvol_in(0,m)
                scavimpnum = scavimptblnum_in(0,m)
             else
                xgrow = log( dumdgratio ) / dlndg_nimptblgrow_in
                jgrow = int( xgrow )
                if (xgrow .lt. 0._r8) jgrow = jgrow - 1
                if (jgrow .lt. nimptblgrow_mind) then
                   jgrow = nimptblgrow_mind
                   xgrow = jgrow
                else
                   jgrow = min( jgrow, nimptblgrow_maxd-1 )
                end if

                dumfhi = xgrow - jgrow
                dumflo = 1._r8 - dumfhi

                scavimpvol = dumflo*scavimptblvol_in(jgrow,m) + &
                     dumfhi*scavimptblvol_in(jgrow+1,m)
                scavimpnum = dumflo*scavimptblnum_in(jgrow,m) + &
                     dumfhi*scavimptblnum_in(jgrow+1,m)

             end if

             ! impaction scavenging removal amount for volume
             scavcoefvol(i,k) = exp( scavimpvol )
             ! impaction scavenging removal amount to number
             scavcoefnum(i,k) = exp( scavimpnum )

             ! scavcoef = impaction scav rate (1/h) for precip = 1 mm/h
             ! scavcoef = impaction scav rate (1/s) for precip = pfx_inrain
             ! (scavcoef/3600) = impaction scav rate (1/s) for precip = 1 mm/h
             ! (pfx_inrain*3600) = in-rain-area precip rate (mm/h)
             ! impactrate = (scavcoef/3600) * (pfx_inrain*3600)
          else
             scavcoefvol(i,k) = 0._r8
             scavcoefnum(i,k) = 0._r8
          end if

       end do
    end do

    return
  end subroutine modal_aero_bcscavcoef_get

end module ap_aero_model_wetdep_scheme
