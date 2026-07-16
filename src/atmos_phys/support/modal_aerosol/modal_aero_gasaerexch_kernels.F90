! Standalone scientific closure for CAM modal aerosol gas/aerosol exchange.
module ap_modal_aero_gasaerexch_kernels

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_const_mod, only : r_universal => shr_const_rgas
  use ap_modal_aero_rename_kernels, only : ap_modal_aero_rename_kernel

  implicit none
  private

  public :: modal_aero_gasaerexch_kernel

contains


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_gasaerexch_sub --- ...
!
! !INTERFACE:
!> \section arg_table_modal_aero_gasaerexch_sub_run Argument Table
!! \htmlinclude modal_aero_gasaerexch_sub_run.html
subroutine modal_aero_gasaerexch_kernel(                         &
                        lchnk, ncol, nstep,                      &
                        pcols_in, pver_in, pcnstxx_in,           &
                        ntot_amode_in, pcnst_in, top_lev_in,     &
                        maxspec_in, ntot_aspectype_in,           &
                        loffset, deltat,                         &
                        l_so4g_in, l_nh4g_in,                   &
                        l_msag_in, l_soag_in,                   &
                        modefrm_pcage_in, modetoo_pcage_in,     &
                        nspecfrm_pcage_in,                      &
                        lspecfrm_pcage_in, lspectoo_pcage_in,   &
                        modeptr_pcarbon_in,                     &
                        numptr_amode_in, nspec_amode_in,        &
                        lmassptr_amode_in, lspectype_amode_in,  &
                        lptr_so4_a_amode_in,                    &
                        lptr_nh4_a_amode_in,                    &
                        lptr_soa_a_amode_in,                    &
                        lptr_pom_a_amode_in,                    &
                        sigmag_amode_in, alnsg_amode_in,        &
                        specmw_amode_in, specdens_amode_in,     &
                        specmw_so4_amode_in,                    &
                        specmw_nh4_amode_in,                    &
                        specmw_soa_amode_in,                    &
                        specdens_so4_amode_in,                  &
                        specdens_nh4_amode_in,                  &
                        specdens_soa_amode_in,                  &
                        dr_so4_monolayers_pcage_in,             &
                        soa_equivso4_factor_in,                 &
                        gravity_in, mwdry_in, rair_in,          &
                        adv_mass_in,                            &
                        t, pmid, pdel, qh2o, troplev,           &
                        q, qqcw, dqdt_other, dqqcwdt_other,     &
                        dgncur_a, dgncur_awet,                  &
                        has_sulfeq, sulfeq,                     &
                        qsrflx, qqcwsrflx,                     &
                        dotend, dotendqqcw,                    &
                        dotendrn, dotendqqcwrn,                &
                        ierr_out, rename_error_out             )


implicit none

! !PARAMETERS:
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: pcols_in             ! declared horizontal dimension
   integer,  intent(in)    :: pver_in              ! vertical layer dimension
   integer,  intent(in)    :: pcnstxx_in           ! gas/aerosol species dimension
   integer,  intent(in)    :: ntot_amode_in        ! aerosol mode dimension
   integer,  intent(in)    :: pcnst_in, top_lev_in
   integer,  intent(in)    :: maxspec_in, ntot_aspectype_in
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   integer,  intent(in)    :: troplev(pcols_in)    ! tropopause vertical index
   integer,  intent(in)    :: l_so4g_in, l_nh4g_in, l_msag_in, l_soag_in
   integer,  intent(in)    :: modefrm_pcage_in, modetoo_pcage_in
   integer,  intent(in)    :: nspecfrm_pcage_in, modeptr_pcarbon_in
   integer,  intent(in)    :: lspecfrm_pcage_in(maxspec_in)
   integer,  intent(in)    :: lspectoo_pcage_in(maxspec_in)
   integer,  intent(in)    :: numptr_amode_in(ntot_amode_in)
   integer,  intent(in)    :: nspec_amode_in(ntot_amode_in)
   integer,  intent(in)    :: lmassptr_amode_in(maxspec_in,ntot_amode_in)
   integer,  intent(in)    :: lspectype_amode_in(maxspec_in,ntot_amode_in)
   integer,  intent(in)    :: lptr_so4_a_amode_in(ntot_amode_in)
   integer,  intent(in)    :: lptr_nh4_a_amode_in(ntot_amode_in)
   integer,  intent(in)    :: lptr_soa_a_amode_in(ntot_amode_in)
   integer,  intent(in)    :: lptr_pom_a_amode_in(ntot_amode_in)
   real(r8), intent(in)    :: deltat               ! time step (s)
   real(r8), intent(in)    :: sigmag_amode_in(ntot_amode_in)
   real(r8), intent(in)    :: alnsg_amode_in(ntot_amode_in)
   real(r8), intent(in)    :: specmw_amode_in(ntot_aspectype_in)
   real(r8), intent(in)    :: specdens_amode_in(ntot_aspectype_in)
   real(r8), intent(in)    :: specmw_so4_amode_in, specmw_nh4_amode_in
   real(r8), intent(in)    :: specmw_soa_amode_in
   real(r8), intent(in)    :: specdens_so4_amode_in
   real(r8), intent(in)    :: specdens_nh4_amode_in
   real(r8), intent(in)    :: specdens_soa_amode_in
   real(r8), intent(in)    :: dr_so4_monolayers_pcage_in
   real(r8), intent(in)    :: soa_equivso4_factor_in
   real(r8), intent(in)    :: gravity_in, mwdry_in, rair_in
   real(r8), intent(in)    :: adv_mass_in(pcnstxx_in)

   real(r8), intent(inout) :: q(ncol,pver_in,pcnstxx_in) ! tracer mixing ratio (TMR) array
                                                   ! *** MUST BE  #/kmol-air for number
                                                   ! *** MUST BE mol/mol-air for mass
                                                   ! *** NOTE ncol dimension
   real(r8), intent(inout) :: qqcw(ncol,pver_in,pcnstxx_in)
                                                   ! like q but for cloud-borner tracers
   real(r8), intent(in)    :: dqdt_other(ncol,pver_in,pcnstxx_in)
                                                   ! TMR tendency from other continuous
                                                   ! growth processes (aqchem, soa??)
                                                   ! *** NOTE ncol dimension
   real(r8), intent(in)    :: dqqcwdt_other(ncol,pver_in,pcnstxx_in)
                                                   ! like dqdt_other but for cloud-borner tracers
   real(r8), intent(in)    :: t(pcols_in,pver_in)        ! temperature at model levels (K)
   real(r8), intent(in)    :: pmid(pcols_in,pver_in)     ! pressure at model levels (Pa)
   real(r8), intent(in)    :: pdel(pcols_in,pver_in)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: qh2o(pcols_in,pver_in)     ! water vapor mixing ratio (kg/kg)
   real(r8), intent(in)    :: dgncur_a(pcols_in,pver_in,ntot_amode_in)
   real(r8), intent(in)    :: dgncur_awet(pcols_in,pver_in,ntot_amode_in)
   logical, intent(in) :: has_sulfeq
   real(r8), intent(in) :: sulfeq(pcols_in,pver_in,ntot_amode_in)
   real(r8), intent(out) :: qsrflx(pcols_in,pcnstxx_in,2)
   real(r8), intent(out) :: qqcwsrflx(pcols_in,pcnstxx_in,2)
   logical, intent(out) :: dotend(pcnstxx_in), dotendqqcw(pcnstxx_in)
   logical, intent(out) :: dotendrn(pcnstxx_in), dotendqqcwrn(pcnstxx_in)
   integer, intent(out) :: ierr_out, rename_error_out

                                 ! dry & wet geo. mean dia. (m) of number distrib.

! !DESCRIPTION:
! computes TMR (tracer mixing ratio) tendencies for gas condensation
!    onto aerosol particles
!
! this version does condensation of H2SO4, NH3, and MSA, both treated as
! completely non-volatile (gas --> aerosol, but no aerosol --> gas)
!    gas H2SO4 goes to aerosol SO4
!    gas MSA (if present) goes to aerosol SO4
!       aerosol MSA is not distinguished from aerosol SO4
!    gas NH3 (if present) goes to aerosol NH4
!       if gas NH3 is not present, then ????
!
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: jsrflx_gaexch = 1
   integer, parameter :: jsrflx_rename = 2
   integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1, ldiag4=-1
   integer, parameter :: method_soa = 2
!     method_soa=0 is no uptake
!     method_soa=1 is irreversible uptake done like h2so4 uptake
!     method_soa=2 is reversible uptake using subr modal_aero_soaexch

   integer :: i, icol_diag, iq
   integer :: idiagss
   integer :: ido_so4a(ntot_amode_in), ido_nh4a(ntot_amode_in)
   integer :: ido_soaa(ntot_amode_in)
   integer :: jac, jsrf
   integer :: k
   integer :: l, l2, lb, lsfrm, lstoo
   integer :: l_so4g, l_nh4g, l_msag, l_soag
   integer :: n, niter, niter_max, ntot_soamode

   logical :: is_dorename_atik, dorename_atik(ncol,pver_in)

   real (r8) :: avg_uprt_nh4, avg_uprt_so4, avg_uprt_soa
   real (r8) :: deltatxx
   real (r8) :: dqdt_nh4(ntot_amode_in), dqdt_so4(ntot_amode_in)
   real (r8) :: dqdt_soa(ntot_amode_in)
   real (r8) :: fac_m2v_nh4, fac_m2v_so4, fac_m2v_soa
   real (r8) :: fac_m2v_pcarbon(maxspec_in)
   real (r8) :: fac_volsfc_pcarbon
   real (r8) :: fgain_nh4(ntot_amode_in), fgain_so4(ntot_amode_in)
   real (r8) :: fgain_soa(ntot_amode_in)
   real (r8) :: g0_soa
   real (r8) :: pdel_fac
   real (r8) :: qmax_nh4, qnew_nh4, qnew_so4
   real (r8) :: qold_nh4(ntot_amode_in), qold_so4(ntot_amode_in)
   real (r8) :: qold_soa(ntot_amode_in), qold_poa(ntot_amode_in)
   real (r8) :: sum_dqdt_msa, sum_dqdt_so4, sum_dqdt_soa
   real (r8) :: sum_dqdt_nh4, sum_dqdt_nh4_b
   real (r8) :: sum_uprt_msa, sum_uprt_nh4, sum_uprt_so4, sum_uprt_soa
   real (r8) :: tmp1, tmp2
   real (r8) :: tmp_kxt, tmp_pxt
   real (r8) :: tmp_so4a_bgn, tmp_so4a_end
   real (r8) :: tmp_so4g_avg, tmp_so4g_bgn, tmp_so4g_equ
   real (r8) :: uptkrate(ntot_amode_in,pcols_in,pver_in)
   real (r8) :: uptkratebb(ntot_amode_in), uptkrate_soa(ntot_amode_in)
                ! gas-to-aerosol mass transfer rates (1/s)
   real (r8) :: vol_core, vol_shell
   real (r8) :: xferfrac_pcage, xferfrac_max
   real (r8) :: xferrate

   logical  :: do_msag         ! true if msa gas is a species
   logical  :: do_nh4g         ! true if nh3 gas is a species
   logical  :: do_soag         ! true if soa gas is a species

   integer, parameter :: nsrflx = 2     ! last dimension of qsrflx
   real(r8) :: dqdt(ncol,pver_in,pcnstxx_in)
   real(r8) :: dqqcwdt(ncol,pver_in,pcnstxx_in)

!  following only needed for diagnostics
   real(r8) :: qold(ncol,pver_in,pcnstxx_in)  ! NOTE dims
   real(r8) :: qnew(ncol,pver_in,pcnstxx_in)  ! NOTE dims
   real(r8) :: qdel(ncol,pver_in,pcnstxx_in)  ! NOTE dims
   real(r8) :: dumavec(1000), dumbvec(1000), dumcvec(1000)
   real(r8) :: qqcwold(ncol,pver_in,pcnstxx_in)
   real(r8) :: dqdtsv1(ncol,pver_in,pcnstxx_in)
   real(r8) :: dqqcwdtsv1(ncol,pver_in,pcnstxx_in)


!----------------------------------------------------------------------

   ierr_out = 0
   rename_error_out = 0
   icol_diag = -1
!!$   if (ldiag1 > 0) then
!!$   if (nstep < 3) then
!!$      do i = 1, ncol
!!$         if ((latndx(i) == 23) .and. (lonndx(i) == 37)) icol_diag = i
!!$      end do
!!$   end if
!!$   end if


! set gas species indices
   l_so4g = l_so4g_in
   l_nh4g = l_nh4g_in
   l_msag = l_msag_in
   l_soag = l_soag_in
   if ((l_so4g <= 0) .or. (l_so4g > pcnstxx_in)) then
      ierr_out = 1
      return
   end if
   do_nh4g = .false.
   do_msag = .false.
   if ((l_nh4g > 0) .and. (l_nh4g <= pcnstxx_in)) do_nh4g = .true.
   if ((l_msag > 0) .and. (l_msag <= pcnstxx_in)) do_msag = .true.
   do_soag = .false.
   if ((method_soa == 1) .or. (method_soa == 2)) then
      if ((l_soag > 0) .and. (l_soag <= pcnstxx_in)) do_soag = .true.
   else if (method_soa /= 0) then
      ierr_out = 2
      return
   end if

! set tendency flags
   dotend(:) = .false.
   dotendqqcw(:) = .false.
   ido_so4a(:) = 0
   ido_nh4a(:) = 0
   ido_soaa(:) = 0

   dotend(l_so4g) = .true.
   if ( do_nh4g ) dotend(l_nh4g) = .true.
   if ( do_msag ) dotend(l_msag) = .true.
   if ( do_soag ) dotend(l_soag) = .true.
   ntot_soamode = 0
   do n = 1, ntot_amode_in
      l = lptr_so4_a_amode_in(n)-loffset
      if ((l > 0) .and. (l <= pcnstxx_in)) then
         dotend(l) = .true.
         ido_so4a(n) = 1
         if ( do_nh4g ) then
            l = lptr_nh4_a_amode_in(n)-loffset
            if ((l > 0) .and. (l <= pcnstxx_in)) then
               dotend(l) = .true.
               ido_nh4a(n) = 1
            end if
         end if
      end if
      if ( do_soag ) then
         l = lptr_soa_a_amode_in(n)-loffset
         if ((l > 0) .and. (l <= pcnstxx_in)) then
            dotend(l) = .true.
            ido_soaa(n) = 1
            ntot_soamode = n
         end if
      end if
   end do
   if ( do_soag ) ntot_soamode = max( ntot_soamode, modefrm_pcage_in )

   if (modefrm_pcage_in > 0) then
      ido_so4a(modefrm_pcage_in) = 2
      if (ido_nh4a(modetoo_pcage_in) == 1) ido_nh4a(modefrm_pcage_in) = 2
      if (ido_soaa(modetoo_pcage_in) == 1) ido_soaa(modefrm_pcage_in) = 2
      do iq = 1, nspecfrm_pcage_in
         lsfrm = lspecfrm_pcage_in(iq)-loffset
         lstoo = lspectoo_pcage_in(iq)-loffset
         if ((lsfrm > 0) .and. (lsfrm <= pcnst_in)) then
            dotend(lsfrm) = .true.
            if ((lstoo > 0) .and. (lstoo <= pcnst_in)) then
               dotend(lstoo) = .true.
            end if
         end if
      end do


      fac_m2v_so4 = specmw_so4_amode_in / specdens_so4_amode_in
      fac_m2v_nh4 = specmw_nh4_amode_in / specdens_nh4_amode_in
      fac_m2v_soa = specmw_soa_amode_in / specdens_soa_amode_in
      fac_m2v_pcarbon(:) = 0.0_r8
      n = modeptr_pcarbon_in
      do l = 1, nspec_amode_in(n)
         l2 = lspectype_amode_in(l,n)
!   fac_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!        [m3-AP/kmol-AP]    = [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
         fac_m2v_pcarbon(l) = specmw_amode_in(l2) / specdens_amode_in(l2)
      end do
      fac_volsfc_pcarbon = exp( 2.5_r8*(alnsg_amode_in(n)**2) )
      xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps
   end if


! zero out tendencies and other
   dqdt(:,:,:) = 0.0_r8
   dqqcwdt(:,:,:) = 0.0_r8
   qsrflx(:,:,:) = 0.0_r8
   qqcwsrflx(:,:,:) = 0.0_r8

! compute gas-to-aerosol mass transfer rates
   call gas_aer_uptkrates( ncol, loffset, pcols_in, pver_in,   &
                           pcnstxx_in, ntot_amode_in,          &
                           top_lev_in, numptr_amode_in,        &
                           sigmag_amode_in, mwdry_in, rair_in, &
                           q, t, pmid, dgncur_awet, uptkrate   )


! use this for tendency calcs to avoid generating very small negative values
   deltatxx = deltat * (1.0_r8 + 1.0e-15_r8)


   jsrf = jsrflx_gaexch
   do k=top_lev_in,pver_in
     do i=1,ncol

!   fgain_so4(n) = fraction of total h2so4 uptake going to mode n
!   fgain_nh4(n) = fraction of total  nh3  uptake going to mode n
        sum_uprt_so4 = 0.0_r8
        sum_uprt_nh4 = 0.0_r8
        sum_uprt_soa = 0.0_r8
        do n = 1, ntot_amode_in
            uptkratebb(n) = uptkrate(n,i,k)
            if (ido_so4a(n) > 0) then
                fgain_so4(n) = uptkratebb(n)
                sum_uprt_so4 = sum_uprt_so4 + fgain_so4(n)
                if (ido_so4a(n) == 1) then
                    qold_so4(n) = q(i,k,lptr_so4_a_amode_in(n)-loffset)
                else
                    qold_so4(n) = 0.0_r8
                end if
            else
                fgain_so4(n) = 0.0_r8
                qold_so4(n) = 0.0_r8
            end if

            if (ido_nh4a(n) > 0) then
!   2.08 factor is for gas diffusivity (nh3/h2so4)
!   differences in fuch-sutugin and accom coef ignored
                fgain_nh4(n) = uptkratebb(n)*2.08_r8
                sum_uprt_nh4 = sum_uprt_nh4 + fgain_nh4(n)
                if (ido_nh4a(n) == 1) then
                    qold_nh4(n) = q(i,k,lptr_nh4_a_amode_in(n)-loffset)
                else
                    qold_nh4(n) = 0.0_r8
                end if
            else
                fgain_nh4(n) = 0.0_r8
                qold_nh4(n) = 0.0_r8
            end if

            if (ido_soaa(n) > 0) then
!   0.81 factor is for gas diffusivity (soa/h2so4)
!   (differences in fuch-sutugin and accom coef ignored)
                fgain_soa(n) = uptkratebb(n)*0.81_r8
                sum_uprt_soa = sum_uprt_soa + fgain_soa(n)
                if (ido_soaa(n) == 1) then
                    qold_soa(n) = q(i,k,lptr_soa_a_amode_in(n)-loffset)
                    l = lptr_pom_a_amode_in(n)-loffset
                    if (l > 0) then
                       qold_poa(n) = q(i,k,l)
                    else
                       qold_poa(n) = 0.0_r8
                    end if
                else
                    qold_soa(n) = 0.0_r8
                    qold_poa(n) = 0.0_r8
                end if
            else
                fgain_soa(n) = 0.0_r8
                qold_soa(n) = 0.0_r8
                qold_poa(n) = 0.0_r8
            end if
            uptkrate_soa(n) = fgain_soa(n)
        end do

        if (sum_uprt_so4 > 0.0_r8) then
            do n = 1, ntot_amode_in
                fgain_so4(n) = fgain_so4(n) / sum_uprt_so4
            end do
        end if
!       at this point (sum_uprt_so4 <= 0.0) only when all the fgain_so4 are zero
        if (sum_uprt_nh4 > 0.0_r8) then
            do n = 1, ntot_amode_in
                fgain_nh4(n) = fgain_nh4(n) / sum_uprt_nh4
            end do
        end if
        if (sum_uprt_soa > 0.0_r8) then
            do n = 1, ntot_amode_in
                fgain_soa(n) = fgain_soa(n) / sum_uprt_soa
            end do
        end if

!   uptake amount (fraction of gas uptaken) over deltat
        avg_uprt_so4 = (1.0_r8 - exp(-deltatxx*sum_uprt_so4))/deltatxx
        avg_uprt_nh4 = (1.0_r8 - exp(-deltatxx*sum_uprt_nh4))/deltatxx
        avg_uprt_soa = (1.0_r8 - exp(-deltatxx*sum_uprt_soa))/deltatxx

!   sum_dqdt_so4 = so4_a tendency from h2so4 gas uptake (mol/mol/s)
!   sum_dqdt_msa = msa_a tendency from msa   gas uptake (mol/mol/s)
!   sum_dqdt_nh4 = nh4_a tendency from nh3   gas uptake (mol/mol/s)
!   sum_dqdt_soa = soa_a tendency from soa   gas uptake (mol/mol/s)
        sum_dqdt_so4 = q(i,k,l_so4g) * avg_uprt_so4
        if ( do_msag ) then
            sum_dqdt_msa = q(i,k,l_msag) * avg_uprt_so4
        else
            sum_dqdt_msa = 0.0_r8
        end if
        if ( do_nh4g ) then
            sum_dqdt_nh4 = q(i,k,l_nh4g) * avg_uprt_nh4
        else
            sum_dqdt_nh4 = 0.0_r8
        end if
        if ( do_soag ) then
            sum_dqdt_soa = q(i,k,l_soag) * avg_uprt_soa
        else
            sum_dqdt_soa = 0.0_r8
        end if

        if ( has_sulfeq .and. (k <= troplev(i)) ) then
           !   compute TMR tendencies for so4 interstial aerosol due to reversible gas uptake
           !   only above the tropopause

           tmp_kxt = deltatxx*sum_uprt_so4  ! sum over modes of uptake_rate*deltat
           tmp_pxt = 0.0_r8
           do n = 1, ntot_amode_in
              if (ido_so4a(n) <= 0) cycle
              tmp_pxt = tmp_pxt + uptkratebb(n)*sulfeq(i,k,n)
           end do
           tmp_pxt = max( 0.0_r8, tmp_pxt*deltatxx )  ! sum over modes of uptake_rate*sulfeq*deltat
           tmp_so4g_bgn = q(i,k,l_so4g)
           ! calc avg h2so4(g) over deltat
           if (tmp_kxt >= 1.0e-5_r8) then
              ! exponential decay towards equilibrium value solution
              tmp_so4g_equ = tmp_pxt/tmp_kxt
              tmp_so4g_avg = tmp_so4g_equ + (tmp_so4g_bgn-tmp_so4g_equ)*(1.0_r8-exp(-tmp_kxt))/tmp_kxt
           else
              ! first order approx for tmp_kxt small
              tmp_so4g_avg = tmp_so4g_bgn*(1.0_r8-0.5_r8*tmp_kxt) + 0.5_r8*tmp_pxt
           end if
           sum_dqdt_so4 = 0.0_r8
           do n = 1, ntot_amode_in
              if (ido_so4a(n) <= 0) cycle
              ! calc change to so4(a) in mode n
              if (ido_so4a(n) == 1) then
                 l = lptr_so4_a_amode_in(n)-loffset
                 tmp_so4a_bgn = q(i,k,l)
              else
                 tmp_so4a_bgn = 0.0_r8
              end if
              tmp_so4a_end = tmp_so4a_bgn + deltatxx*uptkratebb(n)*(tmp_so4g_avg-sulfeq(i,k,n))
              tmp_so4a_end = max( 0.0_r8, tmp_so4a_end )
              dqdt_so4(n) = (tmp_so4a_end - tmp_so4a_bgn)/deltatxx
              sum_dqdt_so4 = sum_dqdt_so4 + dqdt_so4(n)
           end do
           ! do not allow msa condensation in stratosphere
           ! ( Note that the code for msa has never been used.
           !   The plan was to simulate msa(g), treat it as non-volatile (like h2so4(g)),
           !   and treat condensed msa as sulfate, so just one additional tracer. )
           if ( do_msag ) sum_dqdt_msa = 0.0_r8

        else
           !   compute TMR tendencies for so4 interstial aerosol due to simple gas uptake
           do n = 1, ntot_amode_in
              dqdt_so4(n) = fgain_so4(n)*(sum_dqdt_so4 + sum_dqdt_msa)
           end do
        end if

        !   compute TMR tendencies for nh4 interstial aerosol due to simple gas uptake
        !   but force nh4/so4 molar ratio <= 2
        sum_dqdt_nh4_b = 0.0_r8
        dqdt_nh4(:) = 0._r8
        if ( do_nh4g ) then
           do n = 1, ntot_amode_in
              dqdt_nh4(n) = fgain_nh4(n)*sum_dqdt_nh4
              qnew_nh4 = qold_nh4(n) + dqdt_nh4(n)*deltat
              qnew_so4 = qold_so4(n) + dqdt_so4(n)*deltat
              qmax_nh4 = 2.0_r8*qnew_so4
              if (qnew_nh4 > qmax_nh4) then
                 dqdt_nh4(n) = (qmax_nh4 - qold_nh4(n))/deltatxx
              end if
              sum_dqdt_nh4_b = sum_dqdt_nh4_b + dqdt_nh4(n)
           end do
        end if

        if (( do_soag ) .and. (method_soa > 1)) then
!   compute TMR tendencies for soag and soa interstial aerosol
!   using soa parameterization
           niter_max = 1000
           dqdt_soa(:) = 0.0_r8

           call modal_aero_soaexch( deltat, t(i,k), pmid(i,k), &
              niter, niter_max, ntot_soamode, &
              q(i,k,l_soag), qold_soa, qold_poa, uptkrate_soa, &
              tmp1, dqdt_soa )
           sum_dqdt_soa = -tmp1

        else if ( do_soag ) then
!   compute TMR tendencies for soa interstial aerosol
!   due to simple gas uptake
           do n = 1, ntot_amode_in
                dqdt_soa(n) = fgain_soa(n)*sum_dqdt_soa
           end do
        else
           dqdt_soa(:) = 0.0_r8
        end if

        pdel_fac = pdel(i,k)/gravity_in
        do n = 1, ntot_amode_in
            if (ido_so4a(n) == 1) then
                l = lptr_so4_a_amode_in(n)-loffset
                dqdt(i,k,l) = dqdt_so4(n)
                qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(n)*pdel_fac
            end if

            if ( do_nh4g ) then
                if (ido_nh4a(n) == 1) then
                    l = lptr_nh4_a_amode_in(n)-loffset
                    dqdt(i,k,l) = dqdt_nh4(n)
                    qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(n)*pdel_fac
                end if
            end if

            if ( do_soag ) then
                if (ido_soaa(n) == 1) then
                    l = lptr_soa_a_amode_in(n)-loffset
                    dqdt(i,k,l) = dqdt_soa(n)
                    qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(n)*pdel_fac
                end if
            end if
        end do

!   compute TMR tendencies for h2so4, nh3, and msa gas
!   due to simple gas uptake
        l = l_so4g
        dqdt(i,k,l) = -sum_dqdt_so4
        qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac

        if ( do_msag ) then
           l = l_msag
           dqdt(i,k,l) = -sum_dqdt_msa
           qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
        end if

        if ( do_nh4g ) then
           l = l_nh4g
           dqdt(i,k,l) = -sum_dqdt_nh4_b
           qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
        end if

        if ( do_soag ) then
           l = l_soag
           dqdt(i,k,l) = -sum_dqdt_soa
           qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt(i,k,l)*pdel_fac
        end if

!   compute TMR tendencies associated with primary carbon aging
        if (modefrm_pcage_in > 0) then
           n = modeptr_pcarbon_in
           vol_shell = deltat *   &
                     ( dqdt_so4(n)*fac_m2v_so4 + dqdt_nh4(n)*fac_m2v_nh4 +   &
                       dqdt_soa(n)*fac_m2v_soa*soa_equivso4_factor_in )
           vol_core = 0.0_r8
           do l = 1, nspec_amode_in(n)
              vol_core = vol_core + &
                    q(i,k,lmassptr_amode_in(l,n)-loffset)*fac_m2v_pcarbon(l)
           end do
!   ratio1 = vol_shell/vol_core =
!      actual hygroscopic-shell-volume/carbon-core-volume after gas uptake
!   ratio2 = 6.0_r8*dr_so4_monolayers_pcage_in/(dgncur_a*fac_volsfc_pcarbon)
!      = (shell-volume corresponding to n_so4_monolayers_pcage)/core-volume
!      The 6.0/(dgncur_a*fac_volsfc_pcarbon) = (mode-surface-area/mode-volume)
!   Note that vol_shell includes both so4+nh4 AND soa as "equivalent so4",
!      The soa_equivso4_factor_in accounts for the lower hygroscopicity of soa.
!
!   Define xferfrac_pcage = min( 1.0, ratio1/ratio2)
!   But ratio1/ratio2 == tmp1/tmp2, and coding below avoids possible overflow
!
           tmp1 = vol_shell*dgncur_a(i,k,n)*fac_volsfc_pcarbon
           tmp2 = max( 6.0_r8*dr_so4_monolayers_pcage_in*vol_core, 0.0_r8 )
           if (tmp1 >= tmp2) then
              xferfrac_pcage = xferfrac_max
           else
              xferfrac_pcage = min( tmp1/tmp2, xferfrac_max )
           end if

           if (xferfrac_pcage > 0.0_r8) then
              do iq = 1, nspecfrm_pcage_in
                 lsfrm = lspecfrm_pcage_in(iq)-loffset
                 lstoo = lspectoo_pcage_in(iq)-loffset
                 xferrate = (xferfrac_pcage/deltat)*q(i,k,lsfrm)
                 dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xferrate
                 qsrflx(i,lsfrm,jsrf) = qsrflx(i,lsfrm,jsrf) - xferrate*pdel_fac
                 if ((lstoo > 0) .and. (lstoo <= pcnst_in)) then
                     dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xferrate
                     qsrflx(i,lstoo,jsrf) = qsrflx(i,lstoo,jsrf) + xferrate*pdel_fac
                 end if
              end do

              if (ido_so4a(modetoo_pcage_in) > 0) then
                 l = lptr_so4_a_amode_in(modetoo_pcage_in)-loffset
                 dqdt(i,k,l) = dqdt(i,k,l) + dqdt_so4(modefrm_pcage_in)
                 qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_so4(modefrm_pcage_in)*pdel_fac
              end if

              if (ido_nh4a(modetoo_pcage_in) > 0) then
                 l = lptr_nh4_a_amode_in(modetoo_pcage_in)-loffset
                 dqdt(i,k,l) = dqdt(i,k,l) + dqdt_nh4(modefrm_pcage_in)
                 qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_nh4(modefrm_pcage_in)*pdel_fac
              end if

              if (ido_soaa(modetoo_pcage_in) > 0) then
                 l = lptr_soa_a_amode_in(modetoo_pcage_in)-loffset
                 dqdt(i,k,l) = dqdt(i,k,l) + dqdt_soa(modefrm_pcage_in)
                 qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) + dqdt_soa(modefrm_pcage_in)*pdel_fac
              end if
           end if

        end if


! diagnostics start -------------------------------------------------------
!!$   if (ldiag2 > 0) then
!!$   if (i == icol_diag) then
!!$   if (mod(k-1,5) == 0) then
!!$      write(*,'(a,43i5)') 'gasaerexch aaa nstep,lat,lon,k', nstep, latndx(i), lonndx(i), k
!!$      write(*,'(a,1p,10e12.4)') 'uptkratebb   ', uptkratebb(:)
!!$      write(*,'(a,1p,10e12.4)') 'sum_uprt_so4 ', sum_uprt_so4
!!$      write(*,'(a,1p,10e12.4)') 'fgain_so4    ', fgain_so4(:)
!!$      write(*,'(a,1p,10e12.4)') 'sum_uprt_nh4 ', sum_uprt_nh4
!!$      write(*,'(a,1p,10e12.4)') 'fgain_nh4    ', fgain_nh4(:)
!!$      write(*,'(a,1p,10e12.4)') 'sum_uprt_soa ', sum_uprt_soa
!!$      write(*,'(a,1p,10e12.4)') 'fgain_soa    ', fgain_soa(:)
!!$      write(*,'(a,1p,10e12.4)') 'so4g o,dqdt,n', q(i,k,l_so4g), sum_dqdt_so4, &
!!$                                                (q(i,k,l_so4g)-deltat*sum_dqdt_so4)
!!$      write(*,'(a,1p,10e12.4)') 'nh3g o,dqdt,n', q(i,k,l_nh4g), sum_dqdt_nh4, sum_dqdt_nh4_b, &
!!$                                                (q(i,k,l_nh4g)-deltat*sum_dqdt_nh4_b)
!!$      write(*,'(a,1p,10e12.4)') 'soag o,dqdt,n', q(i,k,l_soag), sum_dqdt_soa, &
!!$                                                (q(i,k,l_soag)-deltat*sum_dqdt_soa)
!!$      write(*,'(a,i12,1p,10e12.4)') &
!!$                                'method,g0,t,p', method_soa, g0_soa, t(i,k), pmid(i,k)
!!$      write(*,'(a,1p,10e12.4)') 'so4 old      ', qold_so4(:)
!!$      write(*,'(a,1p,10e12.4)') 'so4 dqdt     ', dqdt_so4(:)
!!$      write(*,'(a,1p,10e12.4)') 'so4 new      ', (qold_so4(:)+deltat*dqdt_so4(:))
!!$      write(*,'(a,1p,10e12.4)') 'nh4 old      ', qold_nh4(:)
!!$      write(*,'(a,1p,10e12.4)') 'nh4 dqdt     ', dqdt_nh4(:)
!!$      write(*,'(a,1p,10e12.4)') 'nh4 new      ', (qold_nh4(:)+deltat*dqdt_nh4(:))
!!$      write(*,'(a,1p,10e12.4)') 'soa old      ', qold_soa(:)
!!$      write(*,'(a,1p,10e12.4)') 'soa dqdt     ', dqdt_soa(:)
!!$      write(*,'(a,1p,10e12.4)') 'soa new      ', (qold_soa(:)+deltat*dqdt_soa(:))
!!$      write(*,'(a,1p,10e12.4)') 'vshell, core ', vol_shell, vol_core
!!$      write(*,'(a,1p,10e12.4)') 'dr_mono, ... ', dr_so4_monolayers_pcage_in,   &
!!$                                 soa_equivso4_factor_in
!!$      write(*,'(a,1p,10e12.4)') 'dgn, ...     ', dgncur_a(i,k,modefrm_pcage_in),   &
!!$                                 fac_volsfc_pcarbon
!!$      write(*,'(a,1p,10e12.4)') 'tmp1, tmp2   ', tmp1, tmp2
!!$      write(*,'(a,1p,10e12.4)') 'xferfrac_age ', xferfrac_pcage
!!$   end if
!!$   end if
!!$   end if
! diagnostics end ---------------------------------------------------------


     end do   ! "i = 1, ncol"
   end do     ! "k = top_lev_in, pver_in"


! set "temporary testing arrays"
   qold(:,:,:) = q(:,:,:)
   qqcwold(:,:,:) = qqcw(:,:,:)
   dqdtsv1(:,:,:) = dqdt(:,:,:)
   dqqcwdtsv1(:,:,:) = dqqcwdt(:,:,:)


!
! do renaming calcs
!
   dotendrn(:) = .false.
   dotendqqcwrn(:) = .false.
   dorename_atik(1:ncol,:) = .true.
   is_dorename_atik = .true.
   call ap_modal_aero_rename_kernel(                        &
        'modal_aero_gasaerexch_sub',            &
        lchnk,             ncol,      nstep,    &
        loffset,           deltat,              &
        gravity_in,                            &
        pdel,              troplev,             &
        dotendrn,          q,                   &
        dqdt,              dqdt_other,          &
        dotendqqcwrn,      qqcw,                &
        dqqcwdt,           dqqcwdt_other,       &
        is_dorename_atik,  dorename_atik,       &
        jsrflx_rename,     nsrflx,              &
        qsrflx,            qqcwsrflx,           &
        rename_error_out                       )
   if (rename_error_out /= 0) return


!
!  apply the dqdt to update q (and same for qqcw)
!
   do l = 1, pcnstxx_in
      if ( dotend(l) .or. dotendrn(l) ) then
         do k = top_lev_in, pver_in
         do i = 1, ncol
            q(i,k,l) = q(i,k,l) + dqdt(i,k,l)*deltat
         end do
         end do
      end if
      if ( dotendqqcw(l) .or. dotendqqcwrn(l) ) then
         do k = top_lev_in, pver_in
         do i = 1, ncol
            qqcw(i,k,l) = qqcw(i,k,l) + dqqcwdt(i,k,l)*deltat
         end do
         end do
      end if
   end do


! diagnostics start -------------------------------------------------------
!!$   if (ldiag3 > 0) then
!!$   if (icol_diag > 0) then
!!$      i = icol_diag
!!$      write(*,'(a,3i5)') 'gasaerexch ppp nstep,lat,lon', nstep, latndx(i), lonndx(i)
!!$      write(*,'(2i5,3(2x,a))') 0, 0, 'ppp', 'pdel for all k'
!!$      write(*,'(1p,7e12.4)') (pdel(i,k), k=top_lev_in,pver_in)
!!$
!!$      write(*,'(a,3i5)') 'gasaerexch ddd nstep,lat,lon', nstep, latndx(i), lonndx(i)
!!$      do l = 1, pcnstxx_in
!!$         lb = l + loffset
!!$
!!$         if ( dotend(l) .or. dotendrn(l) ) then
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd1', cnst_name(lb),    'qold for all k'
!!$            write(*,'(1p,7e12.4)') (qold(i,k,l), k=top_lev_in,pver_in)
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd2', cnst_name(lb),    'qnew for all k'
!!$            write(*,'(1p,7e12.4)') (q(i,k,l), k=top_lev_in,pver_in)
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd3', cnst_name(lb),    'dqdt from conden for all k'
!!$            write(*,'(1p,7e12.4)') (dqdtsv1(i,k,l), k=top_lev_in,pver_in)
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd4', cnst_name(lb),    'dqdt from rename for all k'
!!$            write(*,'(1p,7e12.4)') ((dqdt(i,k,l)-dqdtsv1(i,k,l)), k=top_lev_in,pver_in)
!!$            write(*,'(2i5,3(2x,a))') 1, l, 'ddd5', cnst_name(lb),    'dqdt other for all k'
!!$            write(*,'(1p,7e12.4)') (dqdt_other(i,k,l), k=top_lev_in,pver_in)
!!$         end if
!!$
!!$         if ( dotendqqcw(l) .or. dotendqqcwrn(l) ) then
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd1', cnst_name_cw(lb), 'qold for all k'
!!$            write(*,'(1p,7e12.4)') (qqcwold(i,k,l), k=top_lev_in,pver_in)
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd2', cnst_name_cw(lb), 'qnew for all k'
!!$            write(*,'(1p,7e12.4)') (qqcw(i,k,l), k=top_lev_in,pver_in)
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd3', cnst_name_cw(lb), 'dqdt from conden for all k'
!!$            write(*,'(1p,7e12.4)') (dqqcwdtsv1(i,k,l), k=top_lev_in,pver_in)
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd4', cnst_name_cw(lb), 'dqdt from rename for all k'
!!$            write(*,'(1p,7e12.4)') ((dqqcwdt(i,k,l)-dqqcwdtsv1(i,k,l)), k=top_lev_in,pver_in)
!!$            write(*,'(2i5,3(2x,a))') 2, l, 'ddd5', cnst_name_cw(lb), 'dqdt other for all k'
!!$            write(*,'(1p,7e12.4)') (dqqcwdt_other(i,k,l), k=top_lev_in,pver_in)
!!$         end if
!!$
!!$      end do
!!$
!!$      write(*,'(a,3i5)') 'gasaerexch fff nstep,lat,lon', nstep, latndx(i), lonndx(i)
!!$      do l = 1, pcnstxx_in
!!$         lb = l + loffset
!!$         if ( dotend(l) .or. dotendrn(l) .or. dotendqqcw(l) .or. dotendqqcwrn(l) ) then
!!$            write(*,'(i5,2(2x,a,2l3))') l, &
!!$               cnst_name(lb), dotend(l), dotendrn(l), &
!!$               cnst_name_cw(lb), dotendqqcw(l), dotendqqcwrn(l)
!!$         end if
!!$      end do
!!$
!!$   end if
!!$   end if
! diagnostics end ---------------------------------------------------------


! Scale the returned column tendencies exactly where the historical routine
! performed the conversion; the CAM facade owns only field naming/output.
        do l = 1, pcnstxx_in
           do jsrf = 1, 2
              if ((jsrf == jsrflx_gaexch .and. dotend(l)) .or. &
                  (jsrf == jsrflx_rename .and. dotendrn(l))) then
                 do i = 1, ncol
                    qsrflx(i,l,jsrf) = qsrflx(i,l,jsrf) * &
                         (adv_mass_in(l)/mwdry_in)
                 end do
              end if
              if (jsrf == jsrflx_rename .and. dotendqqcwrn(l)) then
                 do i = 1, ncol
                    qqcwsrflx(i,l,jsrf) = qqcwsrflx(i,l,jsrf) * &
                         (adv_mass_in(l)/mwdry_in)
                 end do
              end if
           end do
        end do


        return
!EOC
   end subroutine modal_aero_gasaerexch_kernel


!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine gas_aer_uptkrates( ncol, loffset, pcols_in, pver_in,   &
                              pcnstxx_in, ntot_amode_in,          &
                              top_lev_in, numptr_amode_in,        &
                              sigmag_amode_in, mwdry_in, rair_in, &
                              q, t, pmid, dgncur_awet, uptkrate   )

!
!                         /
!   computes   uptkrate = | dx  dN/dx  gas_conden_rate(Dp(x))
!                         /
!   using Gauss-Hermite quadrature of order nghq=2
!
!       Dp = particle diameter (cm)
!       x = ln(Dp)
!       dN/dx = log-normal particle number density distribution
!       gas_conden_rate(Dp) = 2 * pi * gasdiffus * Dp * F(Kn,ac)
!           F(Kn,ac) = Fuchs-Sutugin correction factor
!           Kn = Knudsen number
!           ac = accomodation coefficient
!

implicit none


   integer,  intent(in) :: ncol                 ! number of atmospheric column
   integer,  intent(in) :: loffset
   integer,  intent(in) :: pcols_in, pver_in, pcnstxx_in
   integer,  intent(in) :: ntot_amode_in, top_lev_in
   integer,  intent(in) :: numptr_amode_in(ntot_amode_in)
   real(r8), intent(in) :: sigmag_amode_in(ntot_amode_in)
   real(r8), intent(in) :: mwdry_in, rair_in
   real(r8), intent(in) :: q(ncol,pver_in,pcnstxx_in) ! Tracer array (mol,#/mol-air)
   real(r8), intent(in) :: t(pcols_in,pver_in)        ! Temperature in Kelvin
   real(r8), intent(in) :: pmid(pcols_in,pver_in)     ! Air pressure in Pa
   real(r8), intent(in) :: dgncur_awet(pcols_in,pver_in,ntot_amode_in)

   real(r8), intent(out) :: uptkrate(ntot_amode_in,pcols_in,pver_in)
                            ! gas-to-aerosol mass transfer rates (1/s)


! local
   integer, parameter :: nghq = 2
   integer :: i, iq, k, l1, l2, la, n

   ! Can use sqrt here once Lahey is gone.
   real(r8), parameter :: tworootpi = 3.5449077_r8
   real(r8), parameter :: root2 = 1.4142135_r8
   real(r8), parameter :: beta = 2.0_r8

   real(r8) :: aircon
   real(r8) :: const
   real(r8) :: dp, dum_m2v
   real(r8) :: dryvol_a(pcols_in,pver_in)
   real(r8) :: gasdiffus, gasspeed
   real(r8) :: freepathx2, fuchs_sutugin
   real(r8) :: knudsen
   real(r8) :: lndp, lndpgn, lnsg
   real(r8) :: num_a
   real(r8) :: rhoair
   real(r8) :: sumghq
   real(r8), save :: xghq(nghq), wghq(nghq) ! quadrature abscissae and weights

   data xghq / 0.70710678_r8, -0.70710678_r8 /
   data wghq / 0.88622693_r8,  0.88622693_r8 /


! outermost loop over all modes
   do n = 1, ntot_amode_in

! 22-aug-2007 rc easter - get number from q array rather
!    than computing a "bounded" number conc.
!! compute dry volume = sum_over_components{ component_mass / density }
!!    (m3-AP/mol-air)
!! compute it for all i,k to improve accessing q array
!      dryvol_a(1:ncol,:) = 0.0_r8
!      do l1 = 1, nspec_amode_in(n)
!         l2 = lspectype_amode_in(l1,n)
!! dum_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!! [m3-AP/kmol-AP]= [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
!         dum_m2v = specmw_amode_in(l2) / specdens_amode_in(l2)
!         la = lmassptr_amode_in(l1,n)
!         dryvol_a(1:ncol,:) = dryvol_a(1:ncol,:)    &
!                            + max(0.0_r8,q(1:ncol,:,la))*dum_m2v
!      end do

! loops k and i
      do k=top_lev_in,pver_in
      do i=1,ncol

         rhoair = pmid(i,k)/(rair_in*t(i,k))   ! (kg-air/m3)
!        aircon = 1.0e3*rhoair/mwdry_in        ! (mol-air/m3)

!!   "bounded" number conc. (#/m3)
!        num_a = dryvol_a(i,k)*v2ncur_a(i,k,n)*aircon

!   number conc. (#/m3) -- note q(i,k,numptr) is (#/kmol-air)
!   so need aircon in (kmol-air/m3)
         aircon = rhoair/mwdry_in              ! (kmol-air/m3)
         num_a = q(i,k,numptr_amode_in(n)-loffset)*aircon

!   gasdiffus = h2so4 gas diffusivity from mosaic code (m^2/s)
!               (pmid must be Pa)
         gasdiffus = 0.557e-4_r8 * (t(i,k)**1.75_r8) / pmid(i,k)
!   gasspeed = h2so4 gas mean molecular speed from mosaic code (m/s)
         gasspeed  = 1.470e1_r8 * sqrt(t(i,k))
!   freepathx2 = 2 * (h2so4 mean free path)  (m)
         freepathx2 = 6.0_r8*gasdiffus/gasspeed

         lnsg   = log( sigmag_amode_in(n) )
         lndpgn = log( dgncur_awet(i,k,n) )   ! (m)
         const  = tworootpi * num_a * exp(beta*lndpgn + 0.5_r8*(beta*lnsg)**2)

!   sum over gauss-hermite quadrature points
         sumghq = 0.0_r8
         do iq = 1, nghq
            lndp = lndpgn + beta*lnsg**2 + root2*lnsg*xghq(iq)
            dp = exp(lndp)

!   knudsen number
            knudsen = freepathx2/dp
!   following assumes accomodation coefficient = ac = 0.65
!   (Adams & Seinfeld, 2002, JGR, and references therein)
!           fuchs_sutugin = (0.75*ac*(1. + knudsen)) /
!                           (knudsen*(1.0 + knudsen + 0.283*ac) + 0.75*ac)
            fuchs_sutugin = (0.4875_r8*(1._r8 + knudsen)) /   &
                            (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)

            sumghq = sumghq + wghq(iq)*dp*fuchs_sutugin/(dp**beta)
         end do
         uptkrate(n,i,k) = const * gasdiffus * sumghq

      end do   ! "do i = 1, ncol"
      end do   ! "do k = 1, pver_in"

   end do   ! "do n = 1, ntot_soamode"


   return
   end subroutine gas_aer_uptkrates

!----------------------------------------------------------------------


      subroutine modal_aero_soaexch( dtfull, temp, pres, &
          niter, niter_max, ntot_soamode, &
          g_soa_in, a_soa_in, a_poa_in, xferrate, &
          g_soa_tend, a_soa_tend )
!         g_soa_tend, a_soa_tend, g0_soa, idiagss )

!-----------------------------------------------------------------------
!
! Purpose:
!
! calculates condensation/evaporation of "soa gas"
! to/from multiple aerosol modes in 1 grid cell
!
! key assumptions
! (1) ambient equilibrium vapor pressure of soa gas
!     is given by p0_soa_298 and delh_vap_soa
! (2) equilibrium vapor pressure of soa gas at aerosol
!     particle surface is given by raoults law in the form
!     g_star = g0_soa*[a_soa/(a_soa + a_opoa)]
! (3) (oxidized poa)/(total poa) is equal to frac_opoa (constant)
!
!
! Author: R. Easter and R. Zaveri
!
!-----------------------------------------------------------------------
      implicit none

      real(r8), intent(in)  :: dtfull     ! full integration time step (s)
      real(r8), intent(in)  :: temp       ! air temperature (K)
      real(r8), intent(in)  :: pres       ! air pressure (Pa)
      integer,  intent(out) :: niter      ! number of iterations performed
      integer,  intent(in)  :: niter_max  ! max allowed number of iterations
      integer,  intent(in)  :: ntot_soamode             ! number of modes having soa
      real(r8), intent(in)  :: g_soa_in                 ! initial soa gas mixrat (mol/mol)
      real(r8), intent(in)  :: a_soa_in(ntot_soamode)   ! initial soa aerosol mixrat (mol/mol)
      real(r8), intent(in)  :: a_poa_in(ntot_soamode)   ! initial poa aerosol mixrat (mol/mol)
      real(r8), intent(in)  :: xferrate(ntot_soamode)   ! gas-aerosol mass transfer rate (1/s)
      real(r8), intent(out) :: g_soa_tend               ! soa gas mixrat tendency (mol/mol/s)
      real(r8), intent(out) :: a_soa_tend(ntot_soamode) ! soa aerosol mixrat tendency (mol/mol/s)
!     real(r8), intent(out) :: g0_soa   ! ambient soa gas equilib mixrat (mol/mol)
!     integer,  intent(in)  :: idiagss

      integer :: luna=6
      integer :: m

      real(r8), parameter :: alpha = 0.05_r8  ! parameter used in calc of time step
      real(r8), parameter :: g_min1 = 1.0e-20_r8
      real(r8), parameter :: opoa_frac = 0.1_r8  ! fraction of poa that is opoa
      real(r8), parameter :: delh_vap_soa = 156.0e3_r8
      ! delh_vap_soa = heat of vaporization for gas soa (J/mol)
      real(r8), parameter :: p0_soa_298 = 1.0e-10_r8
      real(r8), parameter :: rgas = r_universal*1.e-3_r8

      real(r8) :: a_opoa(ntot_soamode)    ! oxidized-poa aerosol mixrat (mol/mol)
      real(r8) :: a_soa(ntot_soamode)     ! soa aerosol mixrat (mol/mol)
      real(r8) :: a_soa_tmp               ! temporary soa aerosol mixrat (mol/mol)
      real(r8) :: beta(ntot_soamode)      ! dtcur*xferrate
      real(r8) :: dtcur    ! current time step (s)
      real(r8) :: dtmax    ! = (dtfull-tcur)
      real(r8) :: g_soa    ! soa gas mixrat (mol/mol)
      real(r8) :: g0_soa   ! ambient soa gas equilib mixrat (mol/mol)
      real(r8) :: g_star(ntot_soamode)    ! soa gas mixrat that is in equilib
                                          ! with each aerosol mode (mol/mol)
      real(r8) :: phi(ntot_soamode)       ! "relative driving force"
      real(r8) :: p0_soa   ! soa gas equilib vapor presssure (atm)
      real(r8) :: sat(ntot_soamode)
      real(r8) :: tcur     ! current integration time (from 0 s)
      real(r8) :: tmpa, tmpb
      real(r8) :: tot_soa  ! g_soa + sum( a_soa(:) )


! force things to be non-negative and calc tot_soa
! calc a_opoa (always slightly >0)
      g_soa = max( g_soa_in, 0.0_r8 )
      tot_soa = g_soa
      do m = 1, ntot_soamode
         a_soa(m) = max( a_soa_in(m), 0.0_r8 )
         tot_soa = tot_soa + a_soa(m)
         a_opoa(m) = opoa_frac*a_poa_in(m)
         a_opoa(m) = max( a_opoa(m), 1.0e-20_r8 )  ! force to small non-zero value
      end do

! calc ambient equilibrium soa gas
      p0_soa = p0_soa_298 * &
               exp( -(delh_vap_soa/rgas)*((1.0_r8/temp)-(1.0_r8/298.0_r8)) )
      g0_soa = 1.01325e5_r8*p0_soa/pres
! molecular weight adjustment
! the soa parameterization assumes that real gsoa and asoa have mw=150
! currently in cam3,
!    mw=12 is used (this has to do with the mozart preprocessor)
!    soag emission files (molec/cm2/s units) are set to give the desired
!       mass emissions (kg/m2/s) and mass mixing ratios (kg/kg)
!       when mw=12 is applied
!    as a result, the molar mixing ratios for both gsoa and asoa
!       are artificially scaled up by (150/12)
!       and g0_soa must be similarly scaled up
      g0_soa = g0_soa*(150.0_r8/12.0_r8)
!     g0_soa = 0.0   ! force irreversible uptake

      niter = 0
      tcur = 0.0_r8
      dtcur = 0.0_r8
      phi(:) = 0.0_r8
      g_star(:) = 0.0_r8

!     if (idiagss > 0) then
!        write(luna,'(a,1p,10e11.3)') 'p0, g0_soa', p0_soa, g0_soa
!        write(luna,'(3a)') &
!           'niter, tcur,   dtcur,    phi(:),                       ', &
!           'g_star(:),                    ', &
!           'a_soa(:),                     g_soa'
!        write(luna,'(3a)') &
!           '                         sat(:),                       ', &
!           'sat(:)*a_soa(:)               ', &
!           'a_opoa(:)'
!        write(luna,'(i3,1p,20e10.2)') niter, tcur, dtcur, &
!           phi(:), g_star(:), a_soa(:), g_soa
!     end if


! integration loop -- does multiple substeps to reach dtfull
timeloop: do while (tcur < dtfull-1.0e-3_r8 )

      niter = niter + 1
      if (niter > niter_max) exit

      tmpa = 0.0_r8
      do m = 1, ntot_soamode
         sat(m) = g0_soa/(a_soa(m) + a_opoa(m))
         g_star(m) = sat(m)*a_soa(m)
         phi(m) = (g_soa - g_star(m))/max(g_soa,g_star(m),g_min1)
         tmpa = tmpa + xferrate(m)*abs(phi(m))
      end do

      dtmax = dtfull-tcur
      if (dtmax*tmpa <= alpha) then
! here alpha/tmpa >= dtmax, so this is final substep
         dtcur = dtmax
         tcur = dtfull
      else
         dtcur = alpha/tmpa
         tcur = tcur + dtcur
      end if

! step 1 - for modes where soa is condensing, estimate "new" a_soa(m)
!    using an explicit calculation with "old" g_soa
!    and g_star(m) calculated using "old" a_soa(m)
! do this to get better estimate of "new" a_soa(m) and sat(m)
      do m = 1, ntot_soamode
         beta(m) = dtcur*xferrate(m)
         tmpa = g_soa - g_star(m)
         if (tmpa > 0.0_r8) then
            a_soa_tmp = a_soa(m) + beta(m)*tmpa
            sat(m) = g0_soa/(a_soa_tmp + a_opoa(m))
            g_star(m) = sat(m)*a_soa_tmp   ! this just needed for diagnostics
         end if
      end do

! step 2 - implicit in g_soa and semi-implicit in a_soa,
!    with g_star(m) calculated semi-implicitly
      tmpa = 0.0_r8
      tmpb = 0.0_r8
      do m = 1, ntot_soamode
         tmpa = tmpa + a_soa(m)/(1.0_r8 + beta(m)*sat(m))
         tmpb = tmpb + beta(m)/(1.0_r8 + beta(m)*sat(m))
      end do
      g_soa = (tot_soa - tmpa)/(1.0_r8 + tmpb)
      g_soa = max( 0.0_r8, g_soa )
      do m = 1, ntot_soamode
         a_soa(m) = (a_soa(m) + beta(m)*g_soa)/   &
                    (1.0_r8 + beta(m)*sat(m))
      end do

!     if (idiagss > 0) then
!        write(luna,'(i3,1p,20e10.2)') niter, tcur, dtcur, &
!           phi(:), g_star(:), a_soa(:), g_soa
!        write(luna,'(23x,1p,20e10.2)') &
!           sat(:), sat(:)*a_soa(:), a_opoa(:)
!     end if

!     if (niter > 9992000) then
!        write(luna,'(a)') '*** to many iterations'
!        exit
!     end if

      end do timeloop


      g_soa_tend = (g_soa - g_soa_in)/dtfull
      do m = 1, ntot_soamode
         a_soa_tend(m) = (a_soa(m) - a_soa_in(m))/dtfull
      end do


      return
      end subroutine modal_aero_soaexch

!----------------------------------------------------------------------




!----------------------------------------------------------------------

end module ap_modal_aero_gasaerexch_kernels
