! Standalone scientific closure for CAM modal-aerosol mode renaming.
module ap_modal_aero_rename_kernels

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use shr_const_mod, only : pi => shr_const_pi

  implicit none
  private
  save

  public :: ap_modal_aero_rename_configure
  public :: ap_modal_aero_rename_kernel

  integer, parameter :: method_optbb_renamexf = 2

  logical :: modal_accum_coarse_exch = .false.
  logical :: rename_configured = .false.
  integer :: npair_renamexf = 0
  integer :: modeptr_coarse = 0
  integer :: modeptr_accum = 0

  integer, allocatable :: modefrm_renamexf(:)
  integer, allocatable :: modetoo_renamexf(:)
  integer, allocatable :: nspecfrm_renamexf(:)
  integer, allocatable :: lspecfrma_renamexf(:,:)
  integer, allocatable :: lspecfrmc_renamexf(:,:)
  integer, allocatable :: lspectooa_renamexf(:,:)
  integer, allocatable :: lspectooc_renamexf(:,:)
  integer, allocatable :: igrow_shrink_renamexf(:)
  integer, allocatable :: ixferable_all_renamexf(:)
  integer, allocatable :: ixferable_all_needed_renamexf(:)
  integer, allocatable :: ixferable_a_renamexf(:,:)
  integer, allocatable :: ixferable_c_renamexf(:,:)
  integer, allocatable :: ido_mode_calcaa(:)

  integer, allocatable :: nspec_amode(:)
  integer, allocatable :: lspectype_amode(:,:)
  integer, allocatable :: lmassptr_amode(:,:)
  integer, allocatable :: lmassptrcw_amode(:,:)
  integer, allocatable :: numptr_amode(:)
  integer, allocatable :: numptrcw_amode(:)

  real(r8), allocatable :: dp_belowcut(:)
  real(r8), allocatable :: dp_cut(:)
  real(r8), allocatable :: dp_xferall_thresh(:)
  real(r8), allocatable :: dp_xfernone_threshaa(:)
  real(r8), allocatable :: dryvol_smallest(:)
  real(r8), allocatable :: factoraa(:)
  real(r8), allocatable :: factoryy(:)
  real(r8), allocatable :: lndp_cut(:)
  real(r8), allocatable :: factor_3alnsg2(:)
  real(r8), allocatable :: v2nhirlx(:)
  real(r8), allocatable :: v2nlorlx(:)

  real(r8), allocatable :: alnsg_amode(:)
  real(r8), allocatable :: voltonumblo_amode(:)
  real(r8), allocatable :: voltonumbhi_amode(:)
  real(r8), allocatable :: dgnum_amode(:)
  real(r8), allocatable :: specmw_amode(:)
  real(r8), allocatable :: specdens_amode(:)

contains

  subroutine ap_modal_aero_rename_configure(                         &
       modal_accum_coarse_exch_in, npair_in, modeptr_coarse_in,      &
       modeptr_accum_in, modefrm_in, modetoo_in, nspecfrm_in,       &
       lspecfrma_in, lspecfrmc_in, lspectooa_in, lspectooc_in,      &
       igrow_shrink_in, ixferable_all_in, ixferable_all_needed_in,  &
       ixferable_a_in, ixferable_c_in, ido_mode_calcaa_in,          &
       dp_belowcut_in, dp_cut_in, dp_xferall_thresh_in,             &
       dp_xfernone_thresh_in, dryvol_smallest_in, factoraa_in,      &
       factoryy_in, lndp_cut_in, factor_3alnsg2_in, v2nhirlx_in,    &
       v2nlorlx_in, alnsg_in, voltonumblo_in, voltonumbhi_in,       &
       dgnum_in, nspec_in, lspectype_in, specmw_in, specdens_in,    &
       lmassptr_in, lmassptrcw_in, numptr_in, numptrcw_in)

    logical, intent(in) :: modal_accum_coarse_exch_in
    integer, intent(in) :: npair_in, modeptr_coarse_in, modeptr_accum_in
    integer, intent(in) :: modefrm_in(:), modetoo_in(:), nspecfrm_in(:)
    integer, intent(in) :: lspecfrma_in(:,:), lspecfrmc_in(:,:)
    integer, intent(in) :: lspectooa_in(:,:), lspectooc_in(:,:)
    integer, intent(in) :: igrow_shrink_in(:), ixferable_all_in(:)
    integer, intent(in) :: ixferable_all_needed_in(:)
    integer, intent(in) :: ixferable_a_in(:,:), ixferable_c_in(:,:)
    integer, intent(in) :: ido_mode_calcaa_in(:)
    integer, intent(in) :: nspec_in(:), lspectype_in(:,:)
    integer, intent(in) :: lmassptr_in(:,:), lmassptrcw_in(:,:)
    integer, intent(in) :: numptr_in(:), numptrcw_in(:)
    real(r8), intent(in) :: dp_belowcut_in(:), dp_cut_in(:)
    real(r8), intent(in) :: dp_xferall_thresh_in(:)
    real(r8), intent(in) :: dp_xfernone_thresh_in(:)
    real(r8), intent(in) :: dryvol_smallest_in(:), factoraa_in(:)
    real(r8), intent(in) :: factoryy_in(:), lndp_cut_in(:)
    real(r8), intent(in) :: factor_3alnsg2_in(:)
    real(r8), intent(in) :: v2nhirlx_in(:), v2nlorlx_in(:)
    real(r8), intent(in) :: alnsg_in(:), voltonumblo_in(:)
    real(r8), intent(in) :: voltonumbhi_in(:), dgnum_in(:)
    real(r8), intent(in) :: specmw_in(:), specdens_in(:)

    modal_accum_coarse_exch = modal_accum_coarse_exch_in
    npair_renamexf = npair_in
    modeptr_coarse = modeptr_coarse_in
    modeptr_accum = modeptr_accum_in

    modefrm_renamexf = modefrm_in
    modetoo_renamexf = modetoo_in
    nspecfrm_renamexf = nspecfrm_in
    lspecfrma_renamexf = lspecfrma_in
    lspecfrmc_renamexf = lspecfrmc_in
    lspectooa_renamexf = lspectooa_in
    lspectooc_renamexf = lspectooc_in
    igrow_shrink_renamexf = igrow_shrink_in
    ixferable_all_renamexf = ixferable_all_in
    ixferable_all_needed_renamexf = ixferable_all_needed_in
    ixferable_a_renamexf = ixferable_a_in
    ixferable_c_renamexf = ixferable_c_in
    ido_mode_calcaa = ido_mode_calcaa_in

    dp_belowcut = dp_belowcut_in
    dp_cut = dp_cut_in
    dp_xferall_thresh = dp_xferall_thresh_in
    dp_xfernone_threshaa = dp_xfernone_thresh_in
    dryvol_smallest = dryvol_smallest_in
    factoraa = factoraa_in
    factoryy = factoryy_in
    lndp_cut = lndp_cut_in
    factor_3alnsg2 = factor_3alnsg2_in
    v2nhirlx = v2nhirlx_in
    v2nlorlx = v2nlorlx_in

    alnsg_amode = alnsg_in
    voltonumblo_amode = voltonumblo_in
    voltonumbhi_amode = voltonumbhi_in
    dgnum_amode = dgnum_in
    nspec_amode = nspec_in
    lspectype_amode = lspectype_in
    specmw_amode = specmw_in
    specdens_amode = specdens_in
    lmassptr_amode = lmassptr_in
    lmassptrcw_amode = lmassptrcw_in
    numptr_amode = numptr_in
    numptrcw_amode = numptrcw_in
    rename_configured = .true.
  end subroutine ap_modal_aero_rename_configure

  !------------------------------------------------------------------
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine ap_modal_aero_rename_kernel(                 &
       fromwhere,         lchnk,               &
       ncol,              nstep,               &
       loffset,           deltat,              &
       gravity,                                &
       pdel,              troplev,             &
       dotendrn,          q,                   &
       dqdt,              dqdt_other,          &
       dotendqqcwrn,      qqcw,                &
       dqqcwdt,           dqqcwdt_other,       &
       is_dorename_atik,  dorename_atik,       &
       jsrflx_rename,     nsrflx,              &
       qsrflx,            qqcwsrflx,           &
       rename_error,      dqdt_rnpos           )


    ! !PARAMETERS:
    character(len=*), intent(in) :: fromwhere    ! identifies which module
    ! is making the call
    integer,  intent(in)    :: lchnk                ! chunk identifier
    integer,  intent(in)    :: ncol                 ! number of atmospheric column
    integer,  intent(in)    :: nstep                ! model time-step number
    integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
    real(r8), intent(in)    :: deltat               ! time step (s)
    real(r8), intent(in)    :: gravity
    integer,  intent(in)    :: troplev(:)

    real(r8), intent(in)    :: pdel(:,:)             ! pressure thickness of levels (Pa)
    real(r8), intent(in)    :: q(:,:,:)              ! tracer mixing ratio array
    ! *** MUST BE mol/mol-air or #/mol-air
    ! *** NOTE ncol and pcnstxx dimensions
    real(r8), intent(in)    :: qqcw(:,:,:) ! like q but for cloud-borne species

    real(r8), intent(inout) :: dqdt(:,:,:)  ! TMR tendency array;
    ! incoming dqdt = tendencies for the 
    !     "fromwhere" continuous growth process 
    ! the renaming tendencies are added on
    ! *** NOTE ncol and pcnstxx dimensions
    real(r8), intent(inout) :: dqqcwdt(:,:,:)
    real(r8), intent(in)    :: dqdt_other(:,:,:)  
    ! tendencies for "other" continuous growth process 
    ! currently in cam3
    !     dqdt is from gas (h2so4, nh3) condensation
    !     dqdt_other is from aqchem and soa
    ! *** NOTE ncol and pcnstxx dimensions
    real(r8), intent(in)    :: dqqcwdt_other(:,:,:)  
    logical,  intent(inout) :: dotendrn(:) ! identifies the species for which
    !     renaming dqdt is computed
    logical,  intent(inout) :: dotendqqcwrn(:)

    logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
    logical,  intent(in)    :: dorename_atik(:,:) ! true if renaming should
    ! be done at i,k
    integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
    integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

    real(r8), intent(inout) :: qsrflx(:,:,:)
    ! process-specific column tracer tendencies 
    real(r8), intent(inout) :: qqcwsrflx(:,:,:)
    integer, intent(out) :: rename_error
    real(r8), optional, intent(out) &
         :: dqdt_rnpos(:,:,:)
    ! the positive (production) part of the renaming tendency

    rename_error = 0
    if (.not. rename_configured) then
       rename_error = -1
       return
    end if

    if (modal_accum_coarse_exch) then
       call modal_aero_rename_acc_crs_sub(        &
            fromwhere,         lchnk,               &
            ncol,              nstep,               &
            loffset,           deltat,              &
            gravity,                                &
            pdel,              troplev,             &
            dotendrn,          q,                   &
            dqdt,              dqdt_other,          &
            dotendqqcwrn,      qqcw,                &
            dqqcwdt,           dqqcwdt_other,       &
            is_dorename_atik,  dorename_atik,       &
            jsrflx_rename,     nsrflx,              &
            qsrflx,            qqcwsrflx,           &
            rename_error,      dqdt_rnpos           )
    else
       call modal_aero_rename_no_acc_crs_sub(             &
            fromwhere,         lchnk,               &
            ncol,              nstep,               &
            loffset,           deltat,              &
            gravity,                                &
            pdel,                                   &
            dotendrn,          q,                   &
            dqdt,              dqdt_other,          &
            dotendqqcwrn,      qqcw,                &
            dqqcwdt,           dqqcwdt_other,       &
            is_dorename_atik,  dorename_atik,       &
            jsrflx_rename,     nsrflx,              &
            qsrflx,            qqcwsrflx,           &
            rename_error                           )
    endif
  end subroutine ap_modal_aero_rename_kernel

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! private methods
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_rename_no_acc_crs_sub --- ...
!
! !INTERFACE:
	subroutine modal_aero_rename_no_acc_crs_sub(                       &
                        fromwhere,         lchnk,               &
                        ncol,              nstep,               &
                        loffset,           deltat,              &
                        gravity,                                &
                        pdel,                                   &
                        dotendrn,          q,                   &
                        dqdt,              dqdt_other,          &
                        dotendqqcwrn,      qqcw,                &
                        dqqcwdt,           dqqcwdt_other,       &
                        is_dorename_atik,  dorename_atik,       &
                        jsrflx_rename,     nsrflx,              &
                        qsrflx,            qqcwsrflx,           &
                        rename_error                           )

! !USES:
   use shr_spfn_mod, only: erfc => shr_spfn_erfc

   implicit none


! !PARAMETERS:
   character(len=*), intent(in) :: fromwhere    ! identifies which module
                                                ! is making the call
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   real(r8), intent(in)    :: deltat               ! time step (s)
   real(r8), intent(in)    :: gravity

   real(r8), intent(in)    :: pdel(:,:)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: q(:,:,:) ! tracer mixing ratio array
                                                   ! *** MUST BE mol/mol-air or #/mol-air
                                                   ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: qqcw(:,:,:) ! like q but for cloud-borne species

   real(r8), intent(inout) :: dqdt(:,:,:)  ! TMR tendency array;
                              ! incoming dqdt = tendencies for the 
                              !     "fromwhere" continuous growth process 
                              ! the renaming tendencies are added on
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(inout) :: dqqcwdt(:,:,:)
   real(r8), intent(in)    :: dqdt_other(:,:,:)  
                              ! tendencies for "other" continuous growth process 
                              ! currently in cam3
                              !     dqdt is from gas (h2so4, nh3) condensation
                              !     dqdt_other is from aqchem and soa
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: dqqcwdt_other(:,:,:)  
   logical,  intent(inout) :: dotendrn(:) ! identifies the species for which
                              !     renaming dqdt is computed
   logical,  intent(inout) :: dotendqqcwrn(:)

   logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
   logical,  intent(in)    :: dorename_atik(:,:) ! true if renaming should
                                                        ! be done at i,k
   integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
   integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

   real(r8), intent(inout) :: qsrflx(:,:,:)
                              ! process-specific column tracer tendencies 
   real(r8), intent(inout) :: qqcwsrflx(:,:,:)
   integer, intent(out) :: rename_error

! !DESCRIPTION: 
! computes TMR (tracer mixing ratio) tendencies for "mode renaming"
!    during a continuous growth process
! currently this transfers number and mass (and surface) from the aitken
!    to accumulation mode after gas condensation or stratiform-cloud
!    aqueous chemistry
! (convective cloud aqueous chemistry not yet implemented)
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: ldiag1=-1
   integer :: i, icol_diag, ipair, iq, j, k, l, l1, l2, la, lc
   integer :: lsfrma, lsfrmc, lstooa, lstooc
   integer :: mfrm, mtoo, n, n1, n2, ntot_msa_a
   integer :: idomode(size(nspec_amode))


   real (r8) :: deldryvol_a(ncol,size(q,2),size(nspec_amode))
   real (r8) :: deldryvol_c(ncol,size(q,2),size(nspec_amode))
   real (r8) :: deltatinv
   real (r8) :: dp_belowcut(size(modefrm_renamexf))
   real (r8) :: dp_cut(size(modefrm_renamexf))
   real (r8) :: dgn_aftr, dgn_xfer
   real (r8) :: dgn_t_new, dgn_t_old
   real (r8) :: dryvol_t_del, dryvol_t_new
   real (r8) :: dryvol_t_old, dryvol_t_oldbnd
   real (r8) :: dryvol_a(ncol,size(q,2),size(nspec_amode))
   real (r8) :: dryvol_c(ncol,size(q,2),size(nspec_amode))
   real (r8) :: dryvol_smallest(size(nspec_amode))
   real (r8) :: dum
   real (r8) :: dum3alnsg2(size(modefrm_renamexf))
   real (r8) :: dum_m2v, dum_m2vdt
   real (r8) :: factoraa(size(nspec_amode))
   real (r8) :: factoryy(size(nspec_amode))
   real (r8) :: frelax
   real (r8) :: lndp_cut(size(modefrm_renamexf))
   real (r8) :: lndgn_new, lndgn_old
   real (r8) :: lndgv_new, lndgv_old
   real (r8) :: num_t_old, num_t_oldbnd
   real (r8) :: onethird
   real (r8) :: pdel_fac
   real (r8) :: tailfr_volnew, tailfr_volold
   real (r8) :: tailfr_numnew, tailfr_numold
   real (r8) :: v2nhirlx(size(nspec_amode)), v2nlorlx(size(nspec_amode))
   real (r8) :: xfercoef, xfertend
   real (r8) :: xferfrac_vol, xferfrac_num, xferfrac_max

   real (r8) :: yn_tail, yv_tail

! begin
	rename_error = 0


!
!   calculations done once on initial entry
!
!   "init" is now done through chem_init (and things under it)
!	if (npair_renamexf .eq. -123456789) then
!	    npair_renamexf = 0
!	    call modal_aero_rename_init
!	end if

!
!   check if any renaming pairs exist
!
	if (npair_renamexf .le. 0) return
! 	if (ncol .ne. -123456789) return
!	if (fromwhere .eq. 'aqchem') return

!
!   compute aerosol dry-volume for the "from mode" of each renaming pair
!   also compute dry-volume change during the continuous growth process
!	using the incoming dqdt*deltat
!
	deltatinv = 1.0_r8/(deltat*(1.0_r8 + 1.0e-15_r8))
	onethird = 1.0_r8/3.0_r8
	frelax = 27.0_r8
	xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

	do n = 1, size(nspec_amode)
	    idomode(n) = 0
	end do

	do ipair = 1, npair_renamexf
	    if (ipair .gt. 1) goto 8100
	    idomode(modefrm_renamexf(ipair)) = 1

	    mfrm = modefrm_renamexf(ipair)
	    mtoo = modetoo_renamexf(ipair)
	    factoraa(mfrm) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mfrm)**2))
	    factoraa(mtoo) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mtoo)**2))
	    factoryy(mfrm) = sqrt( 0.5_r8 )/alnsg_amode(mfrm)
!   dryvol_smallest is a very small volume mixing ratio (m3-AP/kmol-air)
!   used for avoiding overflow.  it corresponds to dp = 1 nm
!   and number = 1e-5 #/mg-air ~= 1e-5 #/cm3-air
	    dryvol_smallest(mfrm) = 1.0e-25_r8
	    v2nlorlx(mfrm) = voltonumblo_amode(mfrm)*frelax
	    v2nhirlx(mfrm) = voltonumbhi_amode(mfrm)/frelax

	    dum3alnsg2(ipair) = 3.0_r8 * (alnsg_amode(mfrm)**2)
	    dp_cut(ipair) = sqrt(   &
		dgnum_amode(mfrm)*exp(1.5_r8*(alnsg_amode(mfrm)**2)) *   &
		dgnum_amode(mtoo)*exp(1.5_r8*(alnsg_amode(mtoo)**2)) )
	    lndp_cut(ipair) = log( dp_cut(ipair) )
	    dp_belowcut(ipair) = 0.99_r8*dp_cut(ipair)
	end do

	do n = 1, size(nspec_amode)
	    if (idomode(n) .gt. 0) then
		dryvol_a(1:ncol,:,n) = 0.0_r8
		dryvol_c(1:ncol,:,n) = 0.0_r8
		deldryvol_a(1:ncol,:,n) = 0.0_r8
		deldryvol_c(1:ncol,:,n) = 0.0_r8
		do l1 = 1, nspec_amode(n)
		    l2 = lspectype_amode(l1,n)
!   dum_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!            [m3-AP/kmol-AP]= [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
		    dum_m2v = specmw_amode(l2) / specdens_amode(l2)
		    dum_m2vdt = dum_m2v*deltat
		    la = lmassptr_amode(l1,n)-loffset
		    if (la > 0) then
		    dryvol_a(1:ncol,:,n) = dryvol_a(1:ncol,:,n)    &
			+ dum_m2v*max( 0.0_r8,   &
                          q(1:ncol,:,la)-deltat*dqdt_other(1:ncol,:,la) )
		    deldryvol_a(1:ncol,:,n) = deldryvol_a(1:ncol,:,n)    &
			+ (dqdt_other(1:ncol,:,la) + dqdt(1:ncol,:,la))*dum_m2vdt
		    end if

		    lc = lmassptrcw_amode(l1,n)-loffset
		    if (lc > 0) then
		    dryvol_c(1:ncol,:,n) = dryvol_c(1:ncol,:,n)    &
			+ dum_m2v*max( 0.0_r8,   &
                          qqcw(1:ncol,:,lc)-deltat*dqqcwdt_other(1:ncol,:,lc) )
		    deldryvol_c(1:ncol,:,n) = deldryvol_c(1:ncol,:,n)    &
			+ (dqqcwdt_other(1:ncol,:,lc) +   &
			         dqqcwdt(1:ncol,:,lc))*dum_m2vdt
		    end if
		end do
	    end if
	end do



!
!   loop over levels and columns to calc the renaming
!
mainloop1_k:  do k = 1, size(q,2)
mainloop1_i:  do i = 1, ncol

!   if dorename_atik is provided, then check if renaming needed at this i,k
	if (is_dorename_atik) then
	    if (.not. dorename_atik(i,k)) cycle mainloop1_i
	end if
	pdel_fac = pdel(i,k)/gravity

!
!   loop over renameing pairs
!
mainloop1_ipair:  do ipair = 1, npair_renamexf

	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)

!   dryvol_t_old is the old total (a+c) dry-volume for the "from" mode 
!	in m^3-AP/kmol-air
!   dryvol_t_new is the new total dry-volume
!	(old/new = before/after the continuous growth)
	dryvol_t_old = dryvol_a(i,k,mfrm) + dryvol_c(i,k,mfrm)
	dryvol_t_del = deldryvol_a(i,k,mfrm) + deldryvol_c(i,k,mfrm)
	dryvol_t_new = dryvol_t_old + dryvol_t_del
	dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )

!   no renaming if dryvol_t_new ~ 0 or dryvol_t_del ~ 0
	if (dryvol_t_new .le. dryvol_smallest(mfrm)) cycle mainloop1_ipair
	if (dryvol_t_del .le. 1.0e-6_r8*dryvol_t_oldbnd) cycle mainloop1_ipair

!   num_t_old is total number in particles/kmol-air
	num_t_old = q(i,k,numptr_amode(mfrm)-loffset)
	num_t_old = num_t_old + qqcw(i,k,numptrcw_amode(mfrm)-loffset)
	num_t_old = max( 0.0_r8, num_t_old )
	dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
	num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx(mfrm), num_t_old )
	num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx(mfrm), num_t_oldbnd )

!   no renaming if dgnum < "base" dgnum, 
	dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa(mfrm)))**onethird
	if (dgn_t_new .le. dgnum_amode(mfrm)) cycle mainloop1_ipair

!   compute new fraction of number and mass in the tail (dp > dp_cut)
	lndgn_new = log( dgn_t_new )
	lndgv_new = lndgn_new + dum3alnsg2(ipair)
	yn_tail = (lndp_cut(ipair) - lndgn_new)*factoryy(mfrm)
	yv_tail = (lndp_cut(ipair) - lndgv_new)*factoryy(mfrm)
	tailfr_numnew = 0.5_r8*erfc( yn_tail )
	tailfr_volnew = 0.5_r8*erfc( yv_tail )

!   compute old fraction of number and mass in the tail (dp > dp_cut)
	dgn_t_old =   &
		(dryvol_t_oldbnd/(num_t_oldbnd*factoraa(mfrm)))**onethird
!   if dgn_t_new exceeds dp_cut, use the minimum of dgn_t_old and 
!   dp_belowcut to guarantee some transfer
	if (dgn_t_new .ge. dp_cut(ipair)) then
	    dgn_t_old = min( dgn_t_old, dp_belowcut(ipair) )
	end if
	lndgn_old = log( dgn_t_old )
	lndgv_old = lndgn_old + dum3alnsg2(ipair)
	yn_tail = (lndp_cut(ipair) - lndgn_old)*factoryy(mfrm)
	yv_tail = (lndp_cut(ipair) - lndgv_old)*factoryy(mfrm)
	tailfr_numold = 0.5_r8*erfc( yn_tail )
	tailfr_volold = 0.5_r8*erfc( yv_tail )

!   transfer fraction is difference between new and old tail-fractions
!   transfer fraction for number cannot exceed that of mass
	dum = tailfr_volnew*dryvol_t_new - tailfr_volold*dryvol_t_old
	if (dum .le. 0.0_r8) cycle mainloop1_ipair

	xferfrac_vol = min( dum, dryvol_t_new )/dryvol_t_new
	xferfrac_vol = min( xferfrac_vol, xferfrac_max ) 
	xferfrac_num = tailfr_numnew - tailfr_numold
	xferfrac_num = max( 0.0_r8, min( xferfrac_num, xferfrac_vol ) )

!   diagnostic output start ----------------------------------------
!!$ 	if (ldiag1 > 0) then
!!$ 	icol_diag = -1
!!$ 	if ((lonndx(i) == 37) .and. (latndx(i) == 23)) icol_diag = i
!!$ 	if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
!!$ !	write(lun,97010) fromwhere, nstep, lchnk, i, k, ipair
!!$ 	write(lun,97010) fromwhere, nstep, latndx(i), lonndx(i), k, ipair
!!$ 	write(lun,97020) 'drv old/oldbnd/new/del     ',   &
!!$ 		dryvol_t_old, dryvol_t_oldbnd, dryvol_t_new, dryvol_t_del
!!$ 	write(lun,97020) 'num old/oldbnd, dgnold/new ',   &
!!$ 		num_t_old, num_t_oldbnd, dgn_t_old, dgn_t_new
!!$ 	write(lun,97020) 'tailfr v_old/new, n_old/new',   &
!!$ 		tailfr_volold, tailfr_volnew, tailfr_numold, tailfr_numnew
!!$ 	dum = max(1.0e-10_r8,xferfrac_vol) / max(1.0e-10_r8,xferfrac_num)
!!$ 	dgn_xfer = dgn_t_new * dum**onethird
!!$ 	dum = max(1.0e-10_r8,(1.0_r8-xferfrac_vol)) /   &
!!$               max(1.0e-10_r8,(1.0_r8-xferfrac_num))
!!$ 	dgn_aftr = dgn_t_new * dum**onethird
!!$ 	write(lun,97020) 'xferfrac_v/n; dgn_xfer/aftr',   &
!!$ 		xferfrac_vol, xferfrac_num, dgn_xfer, dgn_aftr
!!$ !97010	format( / 'RENAME ', a, '  nx,lc,i,k,ip', i8, 4i4 )
!!$ 97010	format( / 'RENAME ', a, '  nx,lat,lon,k,ip', i8, 4i4 )
!!$ 97020	format( a, 6(1pe15.7) )
!!$ 	end if
!!$ 	end if
!   diagnostic output end   ------------------------------------------


!
!   compute tendencies for the renaming transfer
!
	j = jsrflx_rename
	do iq = 1, nspecfrm_renamexf(ipair)
	    xfercoef = xferfrac_vol*deltatinv
	    if (iq .eq. 1) xfercoef = xferfrac_num*deltatinv

	    lsfrma = lspecfrma_renamexf(iq,ipair)-loffset
	    lsfrmc = lspecfrmc_renamexf(iq,ipair)-loffset
	    lstooa = lspectooa_renamexf(iq,ipair)-loffset
	    lstooc = lspectooc_renamexf(iq,ipair)-loffset

	    if (lsfrma .gt. 0) then
		xfertend = xfercoef*max( 0.0_r8,   &
			    (q(i,k,lsfrma)+dqdt(i,k,lsfrma)*deltat) )

			dqdt(i,k,lsfrma) = dqdt(i,k,lsfrma) - xfertend
		qsrflx(i,lsfrma,j) = qsrflx(i,lsfrma,j) - xfertend*pdel_fac
		if (lstooa .gt. 0) then
		    dqdt(i,k,lstooa) = dqdt(i,k,lstooa) + xfertend
		    qsrflx(i,lstooa,j) = qsrflx(i,lstooa,j) + xfertend*pdel_fac
		end if
	    end if

	    if (lsfrmc .gt. 0) then
		xfertend = xfercoef*max( 0.0_r8,   &
			    (qqcw(i,k,lsfrmc)+dqqcwdt(i,k,lsfrmc)*deltat) )
		dqqcwdt(i,k,lsfrmc) = dqqcwdt(i,k,lsfrmc) - xfertend
		qqcwsrflx(i,lsfrmc,j) = qqcwsrflx(i,lsfrmc,j) - xfertend*pdel_fac
		if (lstooc .gt. 0) then
		    dqqcwdt(i,k,lstooc) = dqqcwdt(i,k,lstooc) + xfertend
		    qqcwsrflx(i,lstooc,j) = qqcwsrflx(i,lstooc,j) + xfertend*pdel_fac
		end if
	    end if

	end do   ! "iq = 1, nspecfrm_renamexf(ipair)"


	end do mainloop1_ipair


	end do mainloop1_i
	end do mainloop1_k

!
!   set dotend's
!
	dotendrn(:) = .false.
	dotendqqcwrn(:) = .false.
	do ipair = 1, npair_renamexf
	do iq = 1, nspecfrm_renamexf(ipair)
	    lsfrma = lspecfrma_renamexf(iq,ipair) - loffset
	    lsfrmc = lspecfrmc_renamexf(iq,ipair) - loffset
	    lstooa = lspectooa_renamexf(iq,ipair) - loffset
	    lstooc = lspectooc_renamexf(iq,ipair) - loffset
	    if (lsfrma .gt. 0) then
		dotendrn(lsfrma) = .true.
		if (lstooa .gt. 0) dotendrn(lstooa) = .true.
	    end if
	    if (lsfrmc .gt. 0) then
		dotendqqcwrn(lsfrmc) = .true.
		if (lstooc .gt. 0) dotendqqcwrn(lstooc) = .true.
	    end if
	end do
	end do


	return


!
!   error -- renaming currently just works for 1 pair
!
8100	rename_error = ipair
        return

!EOC
	end subroutine modal_aero_rename_no_acc_crs_sub



!-------------------------------------------------------------------------

!----------------------------------------------------------------------
! code for troposphere and stratosphere
! -- allows accumulation to coarse mode exchange
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_rename_acc_crs_sub --- ...
!
! !INTERFACE:
	subroutine modal_aero_rename_acc_crs_sub(                       &
                        fromwhere,         lchnk,               &
                        ncol,              nstep,               &
                        loffset,           deltat,              &
                        gravity,                                &
                        pdel,              troplev,             &
                        dotendrn,          q,                   &
                        dqdt,              dqdt_other,          &
                        dotendqqcwrn,      qqcw,                &
                        dqqcwdt,           dqqcwdt_other,       &
                        is_dorename_atik,  dorename_atik,       &
                        jsrflx_rename,     nsrflx,              &
                        qsrflx,            qqcwsrflx,           &
                        rename_error,      dqdt_rnpos           )

! !USES:

   use shr_spfn_mod, only: erfc => shr_spfn_erfc

   implicit none


! !PARAMETERS:
   character(len=*), intent(in) :: fromwhere    ! identifies which module
                                                ! is making the call
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   real(r8), intent(in)    :: deltat               ! time step (s)
   real(r8), intent(in)    :: gravity
   integer,  intent(in)    :: troplev(:)

   real(r8), intent(in)    :: pdel(:,:)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: q(:,:,:) ! tracer mixing ratio array
                                                   ! *** MUST BE mol/mol-air or #/mol-air
                                                   ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: qqcw(:,:,:) ! like q but for cloud-borne species

   real(r8), intent(inout) :: dqdt(:,:,:)  ! TMR tendency array;
                              ! incoming dqdt = tendencies for the 
                              !     "fromwhere" continuous growth process 
                              ! the renaming tendencies are added on
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(inout) :: dqqcwdt(:,:,:)
   real(r8), intent(in)    :: dqdt_other(:,:,:)  
                              ! tendencies for "other" continuous growth process 
                              ! currently in cam3
                              !     dqdt is from gas (h2so4, nh3) condensation
                              !     dqdt_other is from aqchem and soa
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: dqqcwdt_other(:,:,:)  
   logical,  intent(inout) :: dotendrn(:) ! identifies the species for which
                              !     renaming dqdt is computed
   logical,  intent(inout) :: dotendqqcwrn(:)

   logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
   logical,  intent(in)    :: dorename_atik(:,:) ! true if renaming should
                                                        ! be done at i,k
   integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
   integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

   real(r8), intent(inout) :: qsrflx(:,:,:)
                              ! process-specific column tracer tendencies 
   real(r8), intent(inout) :: qqcwsrflx(:,:,:)
   integer, intent(out) :: rename_error
   real(r8), optional, intent(out) &
                           :: dqdt_rnpos(:,:,:)
                              ! the positive (production) part of the renaming tendency

! !DESCRIPTION: 
! computes TMR (tracer mixing ratio) tendencies for "mode renaming"
!    during a continuous growth process
! currently this transfers number and mass (and surface) from the aitken
!    to accumulation mode after gas condensation or stratiform-cloud
!    aqueous chemistry
! (convective cloud aqueous chemistry not yet implemented)
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: ldiag1 = -1
   integer :: i, icol_diag, ipair, iq
   integer :: j, k
   integer :: l, l1, l2, la, lc
   integer :: lsfrma, lsfrmc, lstooa, lstooc
   integer :: mfrm, mtoo, n, n1, n2, ntot_msa_a
   logical :: l_dqdt_rnpos
   logical :: flagaa_shrink, flagbb_shrink

   real (r8) :: deldryvol_a(ncol,size(q,2))
   real (r8) :: deldryvol_c(ncol,size(q,2))
   real (r8) :: deltatinv
   real (r8) :: dgn_aftr, dgn_xfer
   real (r8) :: dgn_t_new, dgn_t_old, dgn_t_oldb
   real (r8) :: dryvol_t_del, dryvol_t_new, dryvol_t_new_xfab
   real (r8) :: dryvol_t_old, dryvol_t_oldb, dryvol_t_oldbnd
   real (r8) :: dryvol_a(ncol,size(q,2))
   real (r8) :: dryvol_c(ncol,size(q,2))
   real (r8) :: dryvol_a_xfab(ncol,size(q,2))
   real (r8) :: dryvol_c_xfab(ncol,size(q,2))
   real (r8) :: dryvol_xferamt
   real (r8) :: lndgn_new, lndgn_old
   real (r8) :: lndgv_new, lndgv_old
   real (r8) :: num_t_old, num_t_oldbnd
   real (r8) :: onethird
   real (r8) :: pdel_fac
   real (r8) :: tailfr_volnew, tailfr_volold
   real (r8) :: tailfr_numnew, tailfr_numold
   real (r8) :: tmpa, tmpf
   real (r8) :: tmp_m2v, tmp_m2vdt
   real (r8) :: xfercoef, xfertend
   real (r8) :: xferfrac_vol, xferfrac_num, xferfrac_max

   real (r8) :: yn_tail, yv_tail

! begin
	rename_error = 0


!
!   calculations done once on initial entry
!
!   "init" is now done through chem_init (and things under it)
!	if (npair_renamexf .eq. -123456789) then
!	    npair_renamexf = 0
!	    call modal_aero_rename_init
!	end if

!
!   check if any renaming pairs exist
!
	if (npair_renamexf .le. 0) return
! 	if (ncol .ne. -123456789) return
!	if (fromwhere .eq. 'aqchem') return


	deltatinv = 1.0_r8/(deltat*(1.0_r8 + 1.0e-15_r8))
	onethird = 1.0_r8/3.0_r8
	xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

	if ( present( dqdt_rnpos ) ) then
	    l_dqdt_rnpos = .true.
	    dqdt_rnpos(:,:,:) = 0.0_r8
	else
	    l_dqdt_rnpos = .false.
	end if



!
!   loop over renaming pairs
!
mainloop1_ipair:  do ipair = 1, npair_renamexf

	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)

	flagaa_shrink = .false.
	if ((mfrm==modeptr_coarse) .and. (mtoo==modeptr_accum)) &
	    flagaa_shrink = .true.

!
!   compute aerosol dry-volume for the "from mode" of each renaming pair
!   also compute dry-volume change during the continuous growth process
!	using the incoming dqdt*deltat
!
	dryvol_a(:,:) = 0.0_r8
	dryvol_c(:,:) = 0.0_r8
	deldryvol_a(:,:) = 0.0_r8
	deldryvol_c(:,:) = 0.0_r8
	if (ixferable_all_renamexf(ipair) <= 0) then
	    dryvol_a_xfab(:,:) = 0.0_r8
	    dryvol_c_xfab(:,:) = 0.0_r8
	end if

	n = mfrm
	do l1 = 1, nspec_amode(n)
	    l2 = lspectype_amode(l1,n)
!   tmp_m2v converts (kmol-AP/kmol-air) to (m3-AP/kmol-air)
!            [m3-AP/kmol-AP]= [kg-AP/kmol-AP]  / [kg-AP/m3-AP]
	    tmp_m2v = specmw_amode(l2) / specdens_amode(l2)
	    tmp_m2vdt = tmp_m2v*deltat
	    la = lmassptr_amode(l1,n)-loffset
	    if (la > 0) then
		dryvol_a(1:ncol,:) = dryvol_a(1:ncol,:)    &
		    + tmp_m2v*max( 0.0_r8,   &
		      q(1:ncol,:,la)-deltat*dqdt_other(1:ncol,:,la) )
		deldryvol_a(1:ncol,:) = deldryvol_a(1:ncol,:)    &
		    + (dqdt_other(1:ncol,:,la) + dqdt(1:ncol,:,la))*tmp_m2vdt
		if ( (ixferable_all_renamexf(ipair) <= 0) .and. &
		     (ixferable_a_renamexf(l1,ipair) > 0) ) then
		    dryvol_a_xfab(1:ncol,:) = dryvol_a_xfab(1:ncol,:)    &
			+ tmp_m2v*max( 0.0_r8,   &
			q(1:ncol,:,la)+deltat*dqdt(1:ncol,:,la) )
		end if
	    end if

	    lc = lmassptrcw_amode(l1,n)-loffset
	    if (lc > 0) then
		dryvol_c(1:ncol,:) = dryvol_c(1:ncol,:)    &
		    + tmp_m2v*max( 0.0_r8,   &
		      qqcw(1:ncol,:,lc)-deltat*dqqcwdt_other(1:ncol,:,lc) )
		deldryvol_c(1:ncol,:) = deldryvol_c(1:ncol,:)    &
		    + (dqqcwdt_other(1:ncol,:,lc) +   &
		             dqqcwdt(1:ncol,:,lc))*tmp_m2vdt
		if ( (ixferable_all_renamexf(ipair) <= 0) .and. &
		     (ixferable_c_renamexf(l1,ipair) > 0) ) then
		    dryvol_c_xfab(1:ncol,:) = dryvol_c_xfab(1:ncol,:)    &
			+ tmp_m2v*max( 0.0_r8,   &
			  qqcw(1:ncol,:,lc)+deltat*dqqcwdt(1:ncol,:,lc) )
		end if
	    end if
	end do

!
!
!   loop over levels and columns to calc the renaming
!
!
mainloop1_k:  do k = 1, size(q,2)
mainloop1_i:  do i = 1, ncol

!   if dorename_atik is provided, then check if renaming needed at this i,k
	if (is_dorename_atik) then
	    if (.not. dorename_atik(i,k)) cycle mainloop1_i
	end if


!   dryvol_t_old is the old total (a+c) dry-volume for the "from" mode 
!	in m^3-AP/kmol-air
!   dryvol_t_new is the new total dry-volume
!	(old/new = before/after the continuous growth)
	dryvol_t_old = dryvol_a(i,k) + dryvol_c(i,k)
	dryvol_t_del = deldryvol_a(i,k) + deldryvol_c(i,k)
	dryvol_t_new = dryvol_t_old + dryvol_t_del
	dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )

grow_shrink_conditional1: &
	if (igrow_shrink_renamexf(ipair) > 0) then
!   do renaming for growing particles

!   no renaming if dryvol_t_new ~ 0
	if (dryvol_t_new .le. dryvol_smallest(mfrm)) cycle mainloop1_i
!   no renaming if delta_dryvol is very small or negative
	if ( (method_optbb_renamexf /= 2) .and. &
	     (dryvol_t_del .le. 1.0e-6_r8*dryvol_t_oldbnd) ) cycle mainloop1_i

!   num_t_old is total number in particles/kmol-air
	num_t_old = q(i,k,numptr_amode(mfrm)-loffset)
	num_t_old = num_t_old + qqcw(i,k,numptrcw_amode(mfrm)-loffset)
	num_t_old = max( 0.0_r8, num_t_old )
	dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
	num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx(mfrm), num_t_old )
	num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx(mfrm), num_t_oldbnd )

!   compute new dgnum
	dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa(mfrm)))**onethird
!   no renaming if dgn_t_new < threshold value
	if (dgn_t_new .le. dp_xfernone_threshaa(ipair)) cycle mainloop1_i

!   compute old dgnum and possibly a smaller value to get more renaming transfer
	dgn_t_old =   &
		(dryvol_t_oldbnd/(num_t_oldbnd*factoraa(mfrm)))**onethird
	dgn_t_oldb = dgn_t_old
	dryvol_t_oldb = dryvol_t_old
	if ( method_optbb_renamexf == 2) then
	    if (dgn_t_old .ge. dp_cut(ipair)) then
		! this revised volume corresponds to dgn_t_old == dp_belowcut, and same number conc
		dryvol_t_oldb = dryvol_t_old * (dp_belowcut(ipair)/dgn_t_old)**3
		dgn_t_oldb = dp_belowcut(ipair)
	    end if
	    if (dgn_t_new .lt. dp_xferall_thresh(ipair)) then
		!   no renaming if delta_dryvol is very small or negative
		if ((dryvol_t_new-dryvol_t_oldb) .le. 1.0e-6_r8*dryvol_t_oldbnd) cycle mainloop1_i
	    end if

	else if (dgn_t_new .ge. dp_cut(ipair)) then
!   if dgn_t_new exceeds dp_cut, use the minimum of dgn_t_oldb and 
!   dp_belowcut to guarantee some transfer
	    dgn_t_oldb = min( dgn_t_oldb, dp_belowcut(ipair) )
	end if

!   compute new fraction of number and mass in the tail (dp > dp_cut)
	lndgn_new = log( dgn_t_new )
	lndgv_new = lndgn_new + factor_3alnsg2(ipair)
	yn_tail = (lndp_cut(ipair) - lndgn_new)*factoryy(mfrm)
	yv_tail = (lndp_cut(ipair) - lndgv_new)*factoryy(mfrm)
	tailfr_numnew = 0.5_r8*erfc( yn_tail )
	tailfr_volnew = 0.5_r8*erfc( yv_tail )

!   compute old fraction of number and mass in the tail (dp > dp_cut)
	lndgn_old = log( dgn_t_oldb )
	lndgv_old = lndgn_old + factor_3alnsg2(ipair)
	yn_tail = (lndp_cut(ipair) - lndgn_old)*factoryy(mfrm)
	yv_tail = (lndp_cut(ipair) - lndgv_old)*factoryy(mfrm)
	tailfr_numold = 0.5_r8*erfc( yn_tail )
	tailfr_volold = 0.5_r8*erfc( yv_tail )

!   transfer fraction is difference between new and old tail-fractions
!   transfer fraction for number cannot exceed that of mass
	if ( (method_optbb_renamexf == 2) .and. &
	     (dgn_t_new .ge. dp_xferall_thresh(ipair)) ) then
	    dryvol_xferamt = dryvol_t_new
	else
	    dryvol_xferamt = tailfr_volnew*dryvol_t_new - tailfr_volold*dryvol_t_oldb
	end if
	if (dryvol_xferamt .le. 0.0_r8) cycle mainloop1_i

	xferfrac_vol = max( 0.0_r8, (dryvol_xferamt/dryvol_t_new) )
	if ( method_optbb_renamexf == 2 .and. &
	     (xferfrac_vol >= xferfrac_max) ) then
	    ! transfer entire contents of mode
	    xferfrac_vol = 1.0_r8
	    xferfrac_num = 1.0_r8
	else
	    xferfrac_vol = min( xferfrac_vol, xferfrac_max ) 
	    xferfrac_num = tailfr_numnew - tailfr_numold
	    xferfrac_num = max( 0.0_r8, min( xferfrac_num, xferfrac_vol ) )
	end if

	if (ixferable_all_renamexf(ipair) <= 0) then
	    ! not all species are xferable
	    dryvol_t_new_xfab = max( 0.0_r8, (dryvol_a_xfab(i,k) + dryvol_c_xfab(i,k)) )
	    dryvol_xferamt = xferfrac_vol*dryvol_t_new
	    if (dryvol_t_new_xfab >= 0.999999_r8*dryvol_xferamt) then
		! xferable dryvol can supply the needed dryvol_xferamt
		! but xferfrac_vol must be increased
		xferfrac_vol = min( 1.0_r8, (dryvol_xferamt/dryvol_t_new_xfab) )
	    else if (dryvol_t_new_xfab >= 1.0e-7_r8*dryvol_xferamt) then
		! xferable dryvol cannot supply the needed dryvol_xferamt
		! so transfer all of it, and reduce the number transfer
		xferfrac_vol = 1.0_r8
		xferfrac_num = xferfrac_num*(dryvol_t_new_xfab/dryvol_xferamt)
	    else
		! xferable dryvol << needed dryvol_xferamt
		cycle mainloop1_i
	    end if
	end if

	else grow_shrink_conditional1
!   do renaming for shrinking particles

!   no renaming if (dryvol_t_old ~ 0)
	if (dryvol_t_old .le. dryvol_smallest(mfrm)) cycle mainloop1_i

!   when (delta_dryvol is very small or positive), 
!      which means particles are not evaporating,
!      only do renaming if [(flagaa_shrink true) and (in stratosphere)]],
!   and set flagbb_shrink true to identify this special case
	if (dryvol_t_del .ge. -1.0e-6_r8*dryvol_t_oldbnd) then
	    if ( ( flagaa_shrink ) .and. ( k < troplev(i) ) ) then
		flagbb_shrink = .true.
	    else
		cycle mainloop1_i
	    end if
	else
	    flagbb_shrink = .false.
	end if

!   num_t_old is total number in particles/kmol-air
	num_t_old = q(i,k,numptr_amode(mfrm)-loffset)
	num_t_old = num_t_old + qqcw(i,k,numptrcw_amode(mfrm)-loffset)
	num_t_old = max( 0.0_r8, num_t_old )
	dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest(mfrm) )
	num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx(mfrm), num_t_old )
	num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx(mfrm), num_t_oldbnd )

!   compute new dgnum
	dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa(mfrm)))**onethird
!   no renaming if (dgn_t_new > xfernone threshold value)
	if (dgn_t_new .ge. dp_xfernone_threshaa(ipair)) cycle mainloop1_i
!   if (flagbb_shrink true), renaming only when (dgn_t_new <= dp_cut value)
	if ( flagbb_shrink ) then
	    if (dgn_t_new .gt. dp_cut(ipair)) cycle mainloop1_i
	end if

	if ( dgn_t_new .le. dp_xferall_thresh(ipair) ) then
!   special case of (dgn_t_new <= xferall threshold value)
	    tailfr_numnew = 1.0_r8
	    tailfr_volnew = 1.0_r8
	else
!   compute new fraction of number and mass in the tail (dp < dp_cut)
	    lndgn_new = log( dgn_t_new )
	    lndgv_new = lndgn_new + factor_3alnsg2(ipair)
	    yn_tail = (lndp_cut(ipair) - lndgn_new)*factoryy(mfrm)
	    yv_tail = (lndp_cut(ipair) - lndgv_new)*factoryy(mfrm)
	    tailfr_numnew = 1.0_r8 - 0.5_r8*erfc( yn_tail )
	    tailfr_volnew = 1.0_r8 - 0.5_r8*erfc( yv_tail )
	end if

!   compute old dgnum
	dgn_t_old =   &
		(dryvol_t_oldbnd/(num_t_oldbnd*factoraa(mfrm)))**onethird
	dgn_t_oldb = dgn_t_old
	dryvol_t_oldb = dryvol_t_old

!   no need to compute old fraction of number and mass in the tail
	tailfr_numold = 0.0_r8
	tailfr_volold = 0.0_r8

!   transfer fraction is new tail-fraction
	xferfrac_vol = tailfr_volnew
	if (xferfrac_vol .le. 0.0_r8) cycle mainloop1_i
	xferfrac_num = tailfr_numnew

	if (xferfrac_vol >= xferfrac_max) then
	    ! transfer entire contents of mode
	    xferfrac_vol = 1.0_r8
	    xferfrac_num = 1.0_r8
	else
	    xferfrac_vol = min( xferfrac_vol, xferfrac_max ) 
!   transfer fraction for number cannot be less than that of volume
	    xferfrac_num = max( xferfrac_num, xferfrac_vol )
	    xferfrac_num = min( xferfrac_max, xferfrac_num )
	end if

	if (ixferable_all_renamexf(ipair) <= 0) then
	    ! not all species are xferable
	    dryvol_t_new_xfab = max( 0.0_r8, (dryvol_a_xfab(i,k) + dryvol_c_xfab(i,k)) )
	    dryvol_xferamt = xferfrac_vol*dryvol_t_new
	    if (dryvol_t_new_xfab >= 0.999999_r8*dryvol_xferamt) then
		! xferable dryvol can supply the needed dryvol_xferamt
		! but xferfrac_vol must be increased
		xferfrac_vol = min( 1.0_r8, (dryvol_xferamt/dryvol_t_new_xfab) )
	    else if (dryvol_t_new_xfab >= 1.0e-7_r8*dryvol_xferamt) then
		! xferable dryvol cannot supply the needed dryvol_xferamt
		! so transfer all of it, and reduce the number transfer
		xferfrac_vol = 1.0_r8
		xferfrac_num = xferfrac_num*(dryvol_t_new_xfab/dryvol_xferamt)
	    else
		! xferable dryvol << needed dryvol_xferamt
		cycle mainloop1_i
	    end if
	end if

	endif grow_shrink_conditional1


!!   diagnostic output start ----------------------------------------
!!	if (ldiag1 > 0) then
! 	icol_diag = -1
! 	if ((lonndx(i) == 37) .and. (latndx(i) == 23)) icol_diag = i
!!	if ((i == icol_diag) .and. (mod(k-1,5) == 0)) then
!! qak
! 	if (ldiag1 <= 0) then
! 	if ((i == 1) .and. (k == 1)) then
!! qak
! !	write(lun,97010) fromwhere, nstep, lchnk, i, k, ipair
!! 	write(lun,97010) fromwhere, nstep, latndx(i), lonndx(i), k, ipair
!! 	write(lun,97020) 'drv old/oldb/oldbnd/new/del    ',   &
!! 		dryvol_t_old, dryvol_t_oldb, dryvol_t_oldbnd, &
!! 		dryvol_t_new, dryvol_t_del
!! 	write(lun,97020) 'num old/oldbnd, dgnold/oldb/new',   &
!! 		num_t_old, num_t_oldbnd, dgn_t_old, dgn_t_oldb, dgn_t_new
!! 	write(lun,97020) 'tailfr v_old/new, n_old/new    ',   &
!! 		tailfr_volold, tailfr_volnew, tailfr_numold, tailfr_numnew
! 	tmpa = max(1.0e-10_r8,xferfrac_vol) / max(1.0e-10_r8,xferfrac_num)
! 	dgn_xfer = dgn_t_new * tmpa**onethird
! 	tmpa = max(1.0e-10_r8,(1.0_r8-xferfrac_vol)) /   &
!               max(1.0e-10_r8,(1.0_r8-xferfrac_num))
!! 	dgn_aftr = dgn_t_new * tmpa**onethird
!! 	write(lun,97020) 'xferfrac_v/n; dgn_xfer/aftr    ',   &
!! 		xferfrac_vol, xferfrac_num, dgn_xfer, dgn_aftr
! !97010	format( / 'RENAME ', a, '  nx,lc,i,k,ip', i8, 4i4 )
! 97010	format( / 'RENAME ', a, '  nx,lat,lon,k,ip', i8, 4i4 )
! 97020	format( a, 6(1pe15.7) )
! 	end if
! 	end if
!   diagnostic output end   ------------------------------------------


!
!   compute tendencies for the renaming transfer
!
	pdel_fac = pdel(i,k)/gravity
	j = jsrflx_rename
	do iq = 1, nspecfrm_renamexf(ipair)
	    xfercoef = xferfrac_vol*deltatinv
	    if (iq .eq. 1) xfercoef = xferfrac_num*deltatinv

	    lsfrma = lspecfrma_renamexf(iq,ipair)-loffset
	    lsfrmc = lspecfrmc_renamexf(iq,ipair)-loffset
	    lstooa = lspectooa_renamexf(iq,ipair)-loffset
	    lstooc = lspectooc_renamexf(iq,ipair)-loffset

	    if (lsfrma .gt. 0) then
		xfertend = xfercoef*max( 0.0_r8,   &
			    (q(i,k,lsfrma)+dqdt(i,k,lsfrma)*deltat) )

		dqdt(i,k,lsfrma) = dqdt(i,k,lsfrma) - xfertend
		qsrflx(i,lsfrma,j) = qsrflx(i,lsfrma,j) - xfertend*pdel_fac
		if (lstooa .gt. 0) then
		    dqdt(i,k,lstooa) = dqdt(i,k,lstooa) + xfertend
		    qsrflx(i,lstooa,j) = qsrflx(i,lstooa,j) + xfertend*pdel_fac
		    if ( l_dqdt_rnpos ) &
			dqdt_rnpos(i,k,lstooa) = dqdt_rnpos(i,k,lstooa) + xfertend
		end if
	    end if

	    if (lsfrmc .gt. 0) then
		xfertend = xfercoef*max( 0.0_r8,   &
			    (qqcw(i,k,lsfrmc)+dqqcwdt(i,k,lsfrmc)*deltat) )
		dqqcwdt(i,k,lsfrmc) = dqqcwdt(i,k,lsfrmc) - xfertend
		qqcwsrflx(i,lsfrmc,j) = qqcwsrflx(i,lsfrmc,j) - xfertend*pdel_fac
		if (lstooc .gt. 0) then
		    dqqcwdt(i,k,lstooc) = dqqcwdt(i,k,lstooc) + xfertend
		    qqcwsrflx(i,lstooc,j) = qqcwsrflx(i,lstooc,j) + xfertend*pdel_fac
		end if
	    end if

	end do   ! "iq = 1, nspecfrm_renamexf(ipair)"


	end do mainloop1_i
	end do mainloop1_k


	end do mainloop1_ipair

!
!   set dotend's
!
	dotendrn(:) = .false.
	dotendqqcwrn(:) = .false.
	do ipair = 1, npair_renamexf
	do iq = 1, nspecfrm_renamexf(ipair)
	    lsfrma = lspecfrma_renamexf(iq,ipair) - loffset
	    lsfrmc = lspecfrmc_renamexf(iq,ipair) - loffset
	    lstooa = lspectooa_renamexf(iq,ipair) - loffset
	    lstooc = lspectooc_renamexf(iq,ipair) - loffset
	    if (lsfrma .gt. 0) then
		dotendrn(lsfrma) = .true.
		if (lstooa .gt. 0) dotendrn(lstooa) = .true.
	    end if
	    if (lsfrmc .gt. 0) then
		dotendqqcwrn(lsfrmc) = .true.
		if (lstooc .gt. 0) dotendqqcwrn(lstooc) = .true.
	    end if
	end do
	end do


	return


!
!   error -- renaming currently just works for 1 pair
!
8100	rename_error = ipair
        return

!EOC
	end subroutine modal_aero_rename_acc_crs_sub



!-------------------------------------------------------------------------
! for modal aerosols in the troposphere and stratophere
! -- allows accumulation to coarse mode exchange
!-------------------------------------------------------------------------

!----------------------------------------------------------------------

end module ap_modal_aero_rename_kernels
