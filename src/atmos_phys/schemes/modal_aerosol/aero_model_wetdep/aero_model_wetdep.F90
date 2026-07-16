module ap_aero_model_wetdep_scheme

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use constituents,   only : pcnst, cnst_name
  use ppgrid,         only : pcols, pver
  use camsrfexch,     only : cam_out_t
  use aerodep_flx,    only : aerodep_flx_prescribed
  use physics_types,  only : physics_state, physics_ptend, physics_ptend_init
  use physics_buffer, only : physics_buffer_desc, pbuf_get_field
  use physconst,      only : gravit
  use cam_history,    only : outfld
  use perf_mod,       only : t_startf, t_stopf
  use modal_aero_data
  use modal_aero_deposition, only : set_srf_wetdep
  use wetdep,                only : wetdepa_v2, wetdep_inputs_set, wetdep_inputs_t
  use modal_aero_calcsize,   only : modal_aero_calcsize_sub
  use modal_aero_wateruptake,only : modal_aero_wateruptake_dr

  implicit none
  private

  public :: aero_model_wetdep_run

  integer, parameter :: nimptblgrow_mind=-7, nimptblgrow_maxd=12

contains

  !> \section arg_table_aero_model_wetdep_run Argument Table
  !! \htmlinclude aero_model_wetdep_run.html
  !!
  subroutine aero_model_wetdep_run(state, dt, dlf, cam_out, ptend, pbuf, &
       nmodes_in, dgnumwet_idx_in, qaerwat_idx_in, fracis_idx_in, &
       sol_facti_cloud_borne_in, sol_factb_interstitial_in, &
       sol_factic_interstitial_in, nwetdep_in, wetdep_lq_in, &
       dlndg_nimptblgrow_in, scavimptblnum_in, scavimptblvol_in)

    use modal_aero_deposition, only: set_srf_wetdep
    use wetdep,                only: wetdepa_v2, wetdep_inputs_set, wetdep_inputs_t
    use modal_aero_data
    use modal_aero_calcsize,   only: modal_aero_calcsize_sub
    use modal_aero_wateruptake,only: modal_aero_wateruptake_dr


    ! args

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:)
    integer,  intent(in) :: nmodes_in
    integer,  intent(in) :: dgnumwet_idx_in, qaerwat_idx_in, fracis_idx_in
    real(r8), intent(in) :: sol_facti_cloud_borne_in
    real(r8), intent(in) :: sol_factb_interstitial_in
    real(r8), intent(in) :: sol_factic_interstitial_in
    integer,  intent(in) :: nwetdep_in
    logical,  intent(in) :: wetdep_lq_in(:)
    real(r8), intent(in) :: dlndg_nimptblgrow_in
    real(r8), intent(in) :: scavimptblnum_in(-7:,:)
    real(r8), intent(in) :: scavimptblvol_in(-7:,:)

    ! local vars

    integer :: m ! tracer index

    integer :: lchnk ! chunk identifier
    integer :: ncol ! number of atmospheric columns

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
    real(r8) :: water_old, water_new ! temporary old/new aerosol water mix-rat
    logical  :: isprx(pcols,pver) ! true if precipation
    real(r8) :: aerdepwetis(pcols,pcnst) ! aerosol wet deposition (interstitial)
    real(r8) :: aerdepwetcw(pcols,pcnst) ! aerosol wet deposition (cloud water)
    real(r8), pointer :: fldcw(:,:)

    real(r8), pointer :: dgnumwet(:,:,:)
    real(r8), pointer :: qaerwat(:,:,:)  ! aerosol water

    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble

    type(wetdep_inputs_t) :: dep_inputs

    lchnk = state%lchnk
    ncol  = state%ncol

    call physics_ptend_init(ptend, state%psetcols, 'aero_model_wetdep', lq=wetdep_lq_in)

    ! Do calculations of mode radius and water uptake if:
    ! 1) modal aerosols are affecting the climate, or
    ! 2) prognostic modal aerosols are enabled

    call t_startf('calcsize')
    ! for prognostic modal aerosols the transfer of mass between aitken and accumulation
    ! modes is done in conjunction with the dry radius calculation
    call modal_aero_calcsize_sub(state, ptend, dt, pbuf)
    call t_stopf('calcsize')

    call t_startf('wateruptake')
    call modal_aero_wateruptake_dr(state, pbuf)
    call t_stopf('wateruptake')

    if (nwetdep_in<1) return

    call wetdep_inputs_set( state, pbuf, dep_inputs )

    call pbuf_get_field(pbuf, dgnumwet_idx_in,       dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes_in/) )
    call pbuf_get_field(pbuf, qaerwat_idx_in,        qaerwat,  start=(/1,1,1/), kount=(/pcols,pver,nmodes_in/) )
    call pbuf_get_field(pbuf, fracis_idx_in,         fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

    prec(:ncol)=0._r8
    do k=1,pver
       where (prec(:ncol) >= 1.e-7_r8)
          isprx(:ncol,k) = .true.
       elsewhere
          isprx(:ncol,k) = .false.
       endwhere
       prec(:ncol) = prec(:ncol) + (dep_inputs%prain(:ncol,k) + dep_inputs%cmfdqr(:ncol,k) - dep_inputs%evapr(:ncol,k)) &
            *state%pdel(:ncol,k)/gravit
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
                tmpdust = max( 0.0_r8, state%q(i,k,lcoardust) + ptend%q(i,k,lcoardust)*dt )
                tmpnacl = max( 0.0_r8, state%q(i,k,lcoarnacl) + ptend%q(i,k,lcoarnacl)*dt )
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
                ptend%lq(mm) = .TRUE.
                dqdt_tmp(:,:) = 0.0_r8
                ! q_tmp reflects changes from modal_aero_calcsize and is the "most current" q
                q_tmp(1:ncol,:) = state%q(1:ncol,:,mm) + ptend%q(1:ncol,:,mm)*dt
                fldcw => qqcw_get_field(pbuf, mm,lchnk)

                call wetdepa_v2( state%pmid, state%q(:,:,1), state%pdel, &
                     dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                     dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                     dep_inputs%evapr, dep_inputs%totcond, q_tmp, dt, &
                     dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                     dlf, fracis(:,:,mm), sol_factb, ncol, &
                     scavcoefnv(:,:,jnv), &
                     is_strat_cloudborne=.false.,  &
                     qqcw=fldcw,  &
                     f_act_conv=f_act_conv, &
                     icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                     sol_facti_in=sol_facti, sol_factic_in=sol_factic )

                ptend%q(1:ncol,:,mm) = ptend%q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

                call outfld( trim(cnst_name(mm))//'WET', dqdt_tmp(:,:), pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIC', icscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SIS', isscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBC', bcscavt, pcols, lchnk)
                call outfld( trim(cnst_name(mm))//'SBS', bsscavt, pcols, lchnk)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFWET', sflx, pcols, lchnk)
                aerdepwetis(:ncol,mm) = sflx(:ncol)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSIC', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSIS', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSBC', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(mm))//'SFSBS', sflx, pcols, lchnk)

                if (lspec > 0) then
                   tmpa = spechygro(lspectype_amode(lspec,m))/ &
                        specdens_amode(lspectype_amode(lspec,m))
                   tmpb = tmpa*dt
                   hygro_sum_old(1:ncol,:) = hygro_sum_old(1:ncol,:) &
                        + tmpa*q_tmp(1:ncol,:)
                   hygro_sum_del(1:ncol,:) = hygro_sum_del(1:ncol,:) &
                        + tmpb*dqdt_tmp(1:ncol,:)
                end if

             else if ((lphase == 1) .and. (lspec == nspec_amode(m)+1)) then
                ! aerosol water -- because of how wetdepa treats evaporation of stratiform
                ! precip, it is not appropriate to apply wetdepa to aerosol water
                ! instead, "hygro_sum" = [sum of (mass*hygro/dens)] is calculated before and
                ! after wet removal, and new water is calculated using
                ! new_water = old_water*min(10,(hygro_sum_new/hygro_sum_old))
                ! the "min(10,...)" is to avoid potential problems when hygro_sum_old ~= 0
                ! also, individual wet removal terms (ic,is,bc,bs) are not output to history
                ! ptend%lq(mm) = .TRUE.
                ! dqdt_tmp(:,:) = 0.0_r8
                do k = 1, pver
                   do i = 1, ncol
                      ! water_old = max( 0.0_r8, state%q(i,k,mm)+ptend%q(i,k,mm)*dt )
                      water_old = max( 0.0_r8, qaerwat(i,k,mm) )
                      hygro_sum_old_ik = max( 0.0_r8, hygro_sum_old(i,k) )
                      hygro_sum_new_ik = max( 0.0_r8, hygro_sum_old_ik+hygro_sum_del(i,k) )
                      if (hygro_sum_new_ik >= 10.0_r8*hygro_sum_old_ik) then
                         water_new = 10.0_r8*water_old
                      else
                         water_new = water_old*(hygro_sum_new_ik/hygro_sum_old_ik)
                      end if
                      ! dqdt_tmp(i,k) = (water_new - water_old)/dt
                      qaerwat(i,k,mm) = water_new
                   end do
                end do

                ! ptend%q(1:ncol,:,mm) = ptend%q(1:ncol,:,mm) + dqdt_tmp(1:ncol,:)

                ! call outfld( trim(cnst_name(mm))

                ! sflx(:)=0._r8
                ! do k=1,pver
                ! do i=1,ncol
                ! sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                ! enddo
                ! enddo
                ! call outfld( trim(cnst_name(mm))

             else ! lphase == 2
                dqdt_tmp(:,:) = 0.0_r8
                fldcw => qqcw_get_field(pbuf, mm,lchnk)

                call wetdepa_v2(state%pmid, state%q(:,:,1), state%pdel, &
                     dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                     dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                     dep_inputs%evapr, dep_inputs%totcond, fldcw, dt, &
                     dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                     dlf, fracis_cw, sol_factb, ncol, &
                     scavcoefnv(:,:,jnv), &
                     is_strat_cloudborne=.true.,  &
                     icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                     sol_facti_in=sol_facti, sol_factic_in=sol_factic )

                fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+dqdt_tmp(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFWET', sflx, pcols, lchnk)
                aerdepwetcw(:ncol,mm) = sflx(:ncol)

                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+icscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFSIC', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+isscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFSIS', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bcscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFSBC', sflx, pcols, lchnk)
                sflx(:)=0._r8
                do k=1,pver
                   do i=1,ncol
                      sflx(i)=sflx(i)+bsscavt(i,k)*state%pdel(i,k)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name_cw(mm))//'SFSBS', sflx, pcols, lchnk)

             endif

          enddo ! lspec = 0, nspec_amode(m)+1
       enddo ! lphase = 1, 2
    enddo ! m = 1, ntot_amode

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not.aerodep_flx_prescribed()) then
       call set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)
    endif

  end subroutine aero_model_wetdep_run

  subroutine modal_aero_bcscavcoef_get( m, ncol, isprx, dgn_awet, scavcoefnum, scavcoefvol, &
       dlndg_nimptblgrow_in, scavimptblnum_in, scavimptblvol_in )

    use modal_aero_data
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
