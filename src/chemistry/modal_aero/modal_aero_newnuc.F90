! modal_aero_newnuc.F90


!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_newnuc --- modal aerosol new-particle nucleation
!
! !INTERFACE:
module modal_aero_newnuc

! !USES:
   use shr_kind_mod,  only:  r8 => shr_kind_r8
   use shr_kind_mod,  only:  r4 => shr_kind_r4
   use mo_constants,  only:  pi
   use chem_mods,     only:  gas_pcnst

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public modal_aero_newnuc_sub, modal_aero_newnuc_init

! !PUBLIC DATA MEMBERS:
  integer, parameter  :: pcnstxx = gas_pcnst
  integer  :: l_h2so4_sv, l_nh3_sv, lnumait_sv, lnh4ait_sv, lso4ait_sv
  logical :: modal_aero_newnuc_zero_tendencies_use_native_impl = .false.
  logical :: modal_aero_newnuc_zero_tendencies_impl_selected = .false.
  logical :: modal_aero_newnuc_prepare_box_inputs_use_native_impl = .false.
  logical :: modal_aero_newnuc_prepare_box_inputs_impl_selected = .false.
  logical :: pbl_nuc_wang2008_use_native_impl = .false.
  logical :: pbl_nuc_wang2008_impl_selected = .false.
  logical :: binary_nuc_vehk2002_use_native_impl = .false.
  logical :: binary_nuc_vehk2002_impl_selected = .false.
  logical :: mer07_veh02_nuc_mosaic_prepare_rates_use_native_impl = .false.
  logical :: mer07_veh02_nuc_mosaic_prepare_rates_impl_selected = .false.
  logical :: mer07_veh02_nuc_mosaic_finalize_use_native_impl = .false.
  logical :: mer07_veh02_nuc_mosaic_finalize_impl_selected = .false.

! min h2so4 vapor for nuc calcs = 4.0e-16 mol/mol-air ~= 1.0e4 molecules/cm3, 
  real(r8), parameter :: qh2so4_cutoff = 4.0e-16_r8

! !DESCRIPTION: This module implements ...
!
! !REVISION HISTORY:
!
!   R.Easter 2007.09.14:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! list private module data here

!EOC
!----------------------------------------------------------------------


  contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_newnuc_zero_tendencies_select_impl()

   use cam_logfile, only: iulog
   use spmd_utils, only: masterproc

   implicit none

   character(len=48) :: impl_name
   integer :: status, n, i, code

   if (modal_aero_newnuc_zero_tendencies_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('MODAL_AERO_NEWNUC_ZERO_TENDENCIES_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      modal_aero_newnuc_zero_tendencies_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      modal_aero_newnuc_zero_tendencies_use_native_impl = .false.
   end if

   modal_aero_newnuc_zero_tendencies_impl_selected = .true.

   if (masterproc) then
      if (modal_aero_newnuc_zero_tendencies_use_native_impl) then
         write(iulog,*) 'modal_aero_newnuc_zero_tendencies implementation = native'
      else
         write(iulog,*) 'modal_aero_newnuc_zero_tendencies implementation = codon'
      end if
   end if

  end subroutine modal_aero_newnuc_zero_tendencies_select_impl

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_newnuc_zero_tendencies(ncol, pcnst_in, nsrflx_in, dqdt, qsrflx)

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
   use ppgrid, only: pcols, pver

   implicit none

   integer, intent(in) :: ncol, pcnst_in, nsrflx_in
   real(r8), target, intent(out) :: dqdt(ncol,pver,pcnstxx)
   real(r8), target, intent(out) :: qsrflx(pcols,pcnst_in,nsrflx_in)

   interface
      subroutine modal_aero_newnuc_zero_tendencies_codon(ncol_c, pcols_c, pver_c, pcnstxx_c, pcnst_c, nsrflx_c, dqdt_p, qsrflx_p) &
           bind(c, name="modal_aero_newnuc_zero_tendencies_codon")
        use iso_c_binding, only: c_int64_t, c_ptr
        integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnstxx_c, pcnst_c, nsrflx_c
        type(c_ptr), value :: dqdt_p, qsrflx_p
      end subroutine modal_aero_newnuc_zero_tendencies_codon
   end interface

   call modal_aero_newnuc_zero_tendencies_select_impl()

   if (modal_aero_newnuc_zero_tendencies_use_native_impl) then
      dqdt(1:ncol,:,:) = 0.0_r8
      qsrflx(1:ncol,:,:) = 0.0_r8
      return
   end if

   call modal_aero_newnuc_zero_tendencies_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnstxx, c_int64_t), int(pcnst_in, c_int64_t), &
        int(nsrflx_in, c_int64_t), c_loc(dqdt(1,1,1)), c_loc(qsrflx(1,1,1)) &
   )

  end subroutine modal_aero_newnuc_zero_tendencies

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_newnuc_prepare_box_inputs_select_impl()

   use cam_logfile, only: iulog
   use spmd_utils, only: masterproc

   implicit none

   character(len=48) :: impl_name
   integer :: status, n, i, code

   if (modal_aero_newnuc_prepare_box_inputs_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('MODAL_AERO_NEWNUC_PREPARE_BOX_INPUTS_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      modal_aero_newnuc_prepare_box_inputs_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      modal_aero_newnuc_prepare_box_inputs_use_native_impl = .false.
   end if

   modal_aero_newnuc_prepare_box_inputs_impl_selected = .true.

   if (masterproc) then
      if (modal_aero_newnuc_prepare_box_inputs_use_native_impl) then
         write(iulog,*) 'modal_aero_newnuc_prepare_box_inputs implementation = native'
      else
         write(iulog,*) 'modal_aero_newnuc_prepare_box_inputs implementation = codon'
      end if
   end if

  end subroutine modal_aero_newnuc_prepare_box_inputs_select_impl

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_newnuc_prepare_box_inputs( &
       ncol, top_lev_in, l_h2so4_in, l_nh3_in, do_nh3_in, deltat, q, qv, cld, qv_sat, del_h2so4_gasprod, del_h2so4_aeruptk, &
       active_mask, cldx_out, qh2so4_cur_out, qh2so4_avg_out, qnh3_cur_out, tmp_uptkrate_out, relhumnn_out)

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr
   use ppgrid, only: pcols, pver

   implicit none

   integer, intent(in) :: ncol, top_lev_in, l_h2so4_in, l_nh3_in
   logical, intent(in) :: do_nh3_in
   real(r8), intent(in) :: deltat
   real(r8), target, intent(in) :: q(ncol,pver,pcnstxx)
   real(r8), target, intent(in) :: qv(pcols,pver)
   real(r8), target, intent(in) :: cld(ncol,pver)
   real(r8), target, intent(in) :: qv_sat(pcols,pver)
   real(r8), target, intent(in) :: del_h2so4_gasprod(ncol,pver)
   real(r8), target, intent(in) :: del_h2so4_aeruptk(ncol,pver)
   integer(c_int64_t), target, intent(out) :: active_mask(ncol,pver)
   real(r8), target, intent(out) :: cldx_out(ncol,pver)
   real(r8), target, intent(out) :: qh2so4_cur_out(ncol,pver)
   real(r8), target, intent(out) :: qh2so4_avg_out(ncol,pver)
   real(r8), target, intent(out) :: qnh3_cur_out(ncol,pver)
   real(r8), target, intent(out) :: tmp_uptkrate_out(ncol,pver)
   real(r8), target, intent(out) :: relhumnn_out(ncol,pver)

   interface
      subroutine modal_aero_newnuc_prepare_box_inputs_codon( &
           ncol_c, pcols_c, pver_c, top_lev_c, l_h2so4_c, l_nh3_c, do_nh3_c, deltat_c, qh2so4_cutoff_c, &
           q_p, qv_p, cld_p, qv_sat_p, del_h2so4_gasprod_p, del_h2so4_aeruptk_p, active_mask_p, cldx_p, &
           qh2so4_cur_p, qh2so4_avg_p, qnh3_cur_p, tmp_uptkrate_p, relhumnn_p) &
           bind(c, name="modal_aero_newnuc_prepare_box_inputs_codon")
        use iso_c_binding, only: c_double, c_int64_t, c_ptr
        integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, top_lev_c, l_h2so4_c, l_nh3_c, do_nh3_c
        real(c_double), value :: deltat_c, qh2so4_cutoff_c
        type(c_ptr), value :: q_p, qv_p, cld_p, qv_sat_p, del_h2so4_gasprod_p, del_h2so4_aeruptk_p
        type(c_ptr), value :: active_mask_p, cldx_p, qh2so4_cur_p, qh2so4_avg_p, qnh3_cur_p, tmp_uptkrate_p, relhumnn_p
      end subroutine modal_aero_newnuc_prepare_box_inputs_codon
   end interface

   integer :: i, k
   integer(c_int64_t) :: do_nh3_c
   real(r8) :: cldx
   real(r8) :: qh2so4_cur, qh2so4_avg, qnh3_cur
   real(r8) :: qvswtr, relhum, relhumav
   real(r8) :: tmpa, tmpb, tmpc
   real(r8) :: tmp_q2, tmp_q3
   real(r8) :: tmp_uptkrate

   call modal_aero_newnuc_prepare_box_inputs_select_impl()

   if (do_nh3_in) then
      do_nh3_c = 1_c_int64_t
   else
      do_nh3_c = 0_c_int64_t
   end if

   if (modal_aero_newnuc_prepare_box_inputs_use_native_impl) then
      active_mask(:,:) = 0_c_int64_t
      cldx_out(:,:) = 0.0_r8
      qh2so4_cur_out(:,:) = 0.0_r8
      qh2so4_avg_out(:,:) = 0.0_r8
      qnh3_cur_out(:,:) = 0.0_r8
      tmp_uptkrate_out(:,:) = 0.0_r8
      relhumnn_out(:,:) = 0.0_r8

      do k = top_lev_in, pver
         do i = 1, ncol
            if (cld(i,k) >= 0.99_r8) cycle

            qh2so4_cur = q(i,k,l_h2so4_in)
            if (qh2so4_cur <= qh2so4_cutoff) cycle

            tmpa = max(0.0_r8, del_h2so4_gasprod(i,k))
            tmp_q3 = qh2so4_cur
            tmp_q2 = tmp_q3 + max(0.0_r8, -del_h2so4_aeruptk(i,k))

            if (tmp_q2 <= tmp_q3) then
               tmpb = 0.0_r8
            else
               tmpc = tmp_q2 * exp(-20.0_r8)
               if (tmp_q3 <= tmpc) then
                  tmp_q3 = tmpc
                  tmpb = 20.0_r8
               else
                  tmpb = log(tmp_q2/tmp_q3)
               end if
            end if

            tmp_uptkrate = tmpb/deltat

            if (tmpb <= 0.1_r8) then
               qh2so4_avg = tmp_q3*(1.0_r8 + 0.5_r8*tmpb) - 0.5_r8*tmpa
            else
               tmpc = tmpa/tmpb
               qh2so4_avg = (tmp_q3 - tmpc)*((exp(tmpb)-1.0_r8)/tmpb) + tmpc
            end if
            if (qh2so4_avg <= qh2so4_cutoff) cycle

            if (do_nh3_in) then
               qnh3_cur = max(0.0_r8, q(i,k,l_nh3_in))
            else
               qnh3_cur = 0.0_r8
            end if

            qvswtr = qv_sat(i,k)
            qvswtr = max(qvswtr, 1.0e-20_r8)
            relhumav = qv(i,k) / qvswtr
            relhumav = max(0.0_r8, min(1.0_r8, relhumav))

            cldx = max(0.0_r8, cld(i,k))
            relhum = (relhumav - cldx) / (1.0_r8 - cldx)
            relhum = max(0.0_r8, min(1.0_r8, relhum))

            active_mask(i,k) = 1_c_int64_t
            cldx_out(i,k) = cldx
            qh2so4_cur_out(i,k) = qh2so4_cur
            qh2so4_avg_out(i,k) = qh2so4_avg
            qnh3_cur_out(i,k) = qnh3_cur
            tmp_uptkrate_out(i,k) = tmp_uptkrate
            relhumnn_out(i,k) = max(0.01_r8, min(0.99_r8, relhum))
         end do
      end do
      return
   end if

   call modal_aero_newnuc_prepare_box_inputs_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(top_lev_in, c_int64_t), &
        int(l_h2so4_in, c_int64_t), int(l_nh3_in, c_int64_t), do_nh3_c, deltat, qh2so4_cutoff, &
        c_loc(q(1,1,1)), c_loc(qv(1,1)), c_loc(cld(1,1)), c_loc(qv_sat(1,1)), c_loc(del_h2so4_gasprod(1,1)), &
        c_loc(del_h2so4_aeruptk(1,1)), c_loc(active_mask(1,1)), c_loc(cldx_out(1,1)), c_loc(qh2so4_cur_out(1,1)), &
        c_loc(qh2so4_avg_out(1,1)), c_loc(qnh3_cur_out(1,1)), c_loc(tmp_uptkrate_out(1,1)), c_loc(relhumnn_out(1,1)) &
   )

  end subroutine modal_aero_newnuc_prepare_box_inputs

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_newnuc_sub --- ...
!
! !INTERFACE:
   subroutine modal_aero_newnuc_sub(                             &
                        lchnk,    ncol,     nstep,               &
                        loffset,  deltat,                        &
                        t,        pmid,     pdel,                &
                        zm,       pblh,                          &
                        qv,       cld,                           &
                        q,                                       &
                        del_h2so4_gasprod,  del_h2so4_aeruptk    )


! !USES:
   use iso_c_binding,     only: c_int64_t
   use modal_aero_data
   use cam_abortutils,    only: endrun
   use cam_history,       only: outfld, fieldname_len
   use chem_mods,         only: adv_mass
   use constituents,      only: pcnst, cnst_name
   use physconst,         only: gravit, mwdry, r_universal
   use ppgrid,            only: pcols, pver
   use spmd_utils,        only: iam, masterproc
   use wv_saturation,     only: qsat
   use ref_pres,          only: top_lev=>clim_modal_aero_top_lev

   implicit none

! !PARAMETERS:
   integer, intent(in)  :: lchnk            ! chunk identifier
   integer, intent(in)  :: ncol             ! number of columns in chunk
   integer, intent(in)  :: nstep            ! model step
   integer, intent(in)  :: loffset          ! offset applied to modal aero "pointers"
   real(r8), intent(in) :: deltat           ! model timestep (s)

   real(r8), intent(in) :: t(pcols,pver)    ! temperature (K)
   real(r8), intent(in) :: pmid(pcols,pver) ! pressure at model levels (Pa)
   real(r8), intent(in) :: pdel(pcols,pver) ! pressure thickness of levels (Pa)
   real(r8), intent(in) :: zm(pcols,pver)   ! midpoint height above surface (m)
   real(r8), intent(in) :: pblh(pcols)      ! pbl height (m)
   real(r8), intent(in) :: qv(pcols,pver)   ! specific humidity (kg/kg)
   real(r8), intent(in) :: cld(ncol,pver)   ! stratiform cloud fraction
                                            ! *** NOTE ncol dimension
   real(r8), intent(inout) :: q(ncol,pver,pcnstxx) 
                                            ! tracer mixing ratio (TMR) array
                                            ! *** MUST BE mol/mol-air or #/mol-air
                                            ! *** NOTE ncol & pcnstxx dimensions
   real(r8), intent(in) :: del_h2so4_gasprod(ncol,pver) 
                                            ! h2so4 gas-phase production
                                            ! change over deltat (mol/mol)
   real(r8), intent(in) :: del_h2so4_aeruptk(ncol,pver) 
                                            ! h2so4 gas-phase loss to
                                            ! aerosol over deltat (mol/mol)

! !DESCRIPTION: 
!   computes changes due to aerosol nucleation (new particle formation)
!       treats both nucleation and subsequent growth of new particles
!	    to aitken mode size
!   uses the following parameterizations
!       vehkamaki et al. (2002) parameterization for binary
!           homogeneous nucleation (h2so4-h2o) plus
!       kerminen and kulmala (2002) parameterization for
!           new particle loss during growth to aitken size
!
! !REVISION HISTORY:
!   R.Easter 2007.09.14:  Adapted from MIRAGE2 code and CMAQ V4.6 code
!
!EOP
!----------------------------------------------------------------------
!BOC

!   local variables
	integer :: i, itmp, k, l, lmz, lun, m, mait
	integer :: lnumait, lso4ait, lnh4ait
	integer :: l_h2so4, l_nh3
	integer :: ldiagveh02
	integer, parameter :: ldiag1=-1, ldiag2=-1, ldiag3=-1, ldiag4=-1
        integer, parameter :: newnuc_method_flagaa = 11
!       integer, parameter :: newnuc_method_flagaa = 12
        !  1=merikanto et al (2007) ternary   2=vehkamaki et al (2002) binary
        ! 11=merikanto ternary + first-order boundary layer
        ! 12=merikanto ternary + second-order boundary layer

	real(r8) :: adjust_factor
	real(r8) :: aircon
	real(r8) :: cldx 
	real(r8) :: dens_nh4so4a
	real(r8) :: dmdt_ait, dmdt_aitsv1, dmdt_aitsv2, dmdt_aitsv3
	real(r8) :: dndt_ait, dndt_aitsv1, dndt_aitsv2, dndt_aitsv3
	real(r8) :: dnh4dt_ait, dso4dt_ait
	real(r8) :: dpnuc
	real(r8) :: dplom_mode(1), dphim_mode(1)
	real(r8) :: ev_sat(pcols,pver)
	real(r8) :: mass1p
	real(r8) :: mass1p_aithi, mass1p_aitlo 
	real(r8) :: mw_so4a_host
	real(r8) :: pdel_fac
	real(r8) :: qh2so4_cur, qh2so4_avg, qh2so4_del
	real(r8) :: qnh3_cur, qnh3_del, qnh4a_del
	real(r8) :: qnuma_del
	real(r8) :: qso4a_del
		real(r8) :: qv_sat(pcols,pver)
		real(r8) :: qvswtr
		real(r8) :: relhum, relhumav, relhumnn
		real(r8) :: tmpa, tmpb, tmpc
		real(r8) :: tmp_q1, tmp_q2, tmp_q3
		real(r8) :: tmp_frso4, tmp_uptkrate
		real(r8) :: cldx_work(ncol,pver)
		real(r8) :: qh2so4_cur_work(ncol,pver), qh2so4_avg_work(ncol,pver)
		real(r8) :: qnh3_cur_work(ncol,pver), tmp_uptkrate_work(ncol,pver)
		real(r8) :: relhumnn_work(ncol,pver)

		integer, parameter :: nsrflx = 1     ! last dimension of qsrflx
		integer(c_int64_t) :: active_mask(ncol,pver)
		real(r8) :: qsrflx(pcols,pcnst,nsrflx)
	                              ! process-specific column tracer tendencies
	                              ! 1 = nucleation (for aerocom)
	real(r8) :: dqdt(ncol,pver,pcnstxx)  ! TMR tendency array -- NOTE dims
	logical  :: dotend(pcnst)            ! flag for doing tendency
	logical  :: do_nh3                   ! flag for doing nh3/nh4


	character(len=1) :: tmpch1, tmpch2, tmpch3
        character(len=fieldname_len+3) :: fieldname


! begin
	lun = 6

!--------------------------------------------------------------------------------
!!$   if (ldiag1 > 0) then
!!$   do i = 1, ncol
!!$   if (lonndx(i) /= 37) cycle
!!$   if (latndx(i) /= 23) cycle
!!$   if (nstep > 3)       cycle
!!$   write( lun, '(/a,i7,3i5,f10.2)' )   &
!!$         '*** modal_aero_newnuc_sub -- nstep, iam, lat, lon =',   &
!!$         nstep, iam, latndx(i), lonndx(i)
!!$   end do
!!$   if (nstep > 3) call endrun( '*** modal_aero_newnuc_sub -- testing halt after step 3' )
!!$!  if (ncol /= -999888777) return
!!$   end if
!--------------------------------------------------------------------------------

!-----------------------------------------------------------------------
	l_h2so4 = l_h2so4_sv - loffset
	l_nh3   = l_nh3_sv   - loffset
	lnumait = lnumait_sv - loffset
	lnh4ait = lnh4ait_sv - loffset
	lso4ait = lso4ait_sv - loffset

!   skip if no aitken mode OR if no h2so4 species
	if ((l_h2so4 <= 0) .or. (lso4ait <= 0) .or. (lnumait <= 0)) return

	dotend(:) = .false.
	call modal_aero_newnuc_zero_tendencies(ncol, pcnst, nsrflx, dqdt, qsrflx)

!   set dotend
	mait = modeptr_aitken
	dotend(lnumait) = .true.
	dotend(lso4ait) = .true.
	dotend(l_h2so4) = .true.

	lnh4ait = lptr_nh4_a_amode(mait) - loffset
	if ((l_nh3   > 0) .and. (l_nh3   <= pcnst) .and. &
	    (lnh4ait > 0) .and. (lnh4ait <= pcnst)) then
	    do_nh3 = .true.
	    dotend(lnh4ait) = .true.
	    dotend(l_nh3) = .true.
	else
	    do_nh3 = .false.
	end if


!   dry-diameter limits for "grown" new particles
	dplom_mode(1) = exp( 0.67_r8*log(dgnumlo_amode(mait))   &
	                   + 0.33_r8*log(dgnum_amode(mait)) )
	dphim_mode(1) = dgnumhi_amode(mait)

!   mass1p_... = mass (kg) of so4 & nh4 in a single particle of diameter ...
!                (assuming same dry density for so4 & nh4)
!	mass1p_aitlo - dp = dplom_mode(1)
!	mass1p_aithi - dp = dphim_mode(1)
	tmpa = specdens_so4_amode*pi/6.0_r8
	mass1p_aitlo = tmpa*(dplom_mode(1)**3)
	mass1p_aithi = tmpa*(dphim_mode(1)**3)

!   compute qv_sat = saturation specific humidity
	call qsat(t(1:ncol, 1:pver), pmid(1:ncol, 1:pver), &
	            ev_sat(1:ncol, 1:pver), qv_sat(1:ncol, 1:pver))

	call modal_aero_newnuc_prepare_box_inputs( &
	     ncol, top_lev, l_h2so4, l_nh3, do_nh3, deltat, q, qv, cld, qv_sat, del_h2so4_gasprod, del_h2so4_aeruptk, &
	     active_mask, cldx_work, qh2so4_cur_work, qh2so4_avg_work, qnh3_cur_work, tmp_uptkrate_work, relhumnn_work)

!   mw_so4a_host is molec-wght of sulfate aerosol in host code
!      96 when nh3/nh4 are simulated
!      something else when nh3/nh4 are not simulated
	mw_so4a_host = specmw_so4_amode


!
!   loop over levels and columns to calc the renaming
!
main_k:	do k = top_lev, pver
main_i:	do i = 1, ncol

	if (active_mask(i,k) == 0_c_int64_t) cycle main_i

	cldx = cldx_work(i,k)
	qh2so4_cur = qh2so4_cur_work(i,k)
	qh2so4_avg = qh2so4_avg_work(i,k)
	qnh3_cur = qnh3_cur_work(i,k)
	tmp_uptkrate = tmp_uptkrate_work(i,k)
	relhumnn = relhumnn_work(i,k)


!   call ... routine to get nucleation rates
 	ldiagveh02 = -1
!!$ 	if (ldiag2 > 0) then
!!$ 	if ((lonndx(i) == 37) .and. (latndx(i) == 23)) then
!!$ 	if ((k >= 24) .or. (mod(k,4) == 0)) then
!!$ 	    ldiagveh02 = +1
!!$            write(lun,'(/a,i8,3i4,f8.2,1p,4e10.2)')   &
!!$ 		'veh02 call - nstep,lat,lon,k; tk,rh,p,cair',   &
!!$ 		nstep, latndx(i), lonndx(i), k,   &
!!$ 		t(i,k), relhumnn, pmid(k,k), aircon
!!$ 	end if
!!$ 	end if
!!$ 	end if
        call mer07_veh02_nuc_mosaic_1box(   &
           newnuc_method_flagaa,   &
           deltat, t(i,k), relhumnn, pmid(i,k),   &
           zm(i,k), pblh(i),   &
           qh2so4_cur, qh2so4_avg, qnh3_cur, tmp_uptkrate,   &
           mw_so4a_host,   &
           1, 1, dplom_mode, dphim_mode,   &
           itmp, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a, ldiagveh02 )
!          qh2so4_del, qnh3_del, dens_nh4so4a )
!----------------------------------------------------------------------
!       subr mer07_veh02_nuc_mosaic_1box(   &
!          newnuc_method_flagaa,   &
!          dtnuc, temp_in, rh_in, press_in,   &
!          qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate,   &
!          nsize, maxd_asize, dplom_sect, dphim_sect,   &
!          isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
!          qh2so4_del, qnh3_del, dens_nh4so4a )
!
!! subr arguments (in)
!        real(r8), intent(in) :: dtnuc             ! nucleation time step (s)
!        real(r8), intent(in) :: temp_in           ! temperature, in k
!        real(r8), intent(in) :: rh_in             ! relative humidity, as fraction
!        real(r8), intent(in) :: press_in          ! air pressure (pa)
!
!        real(r8), intent(in) :: qh2so4_cur, qh2so4_avg
!                                                  ! gas h2so4 mixing ratios (mol/mol-air)
!        real(r8), intent(in) :: qnh3_cur          ! gas nh3 mixing ratios (mol/mol-air)
!             ! qxxx_cur = current value (after gas chem and condensation)
!             ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
!        real(r8), intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)

!
!        integer, intent(in) :: nsize                    ! number of aerosol size bins
!        integer, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
!        real(r8), intent(in) :: dplom_sect(maxd_asize)  ! dry diameter at lower bnd of bin (m)
!        real(r8), intent(in) :: dphim_sect(maxd_asize)  ! dry diameter at upper bnd of bin (m)
!
!! subr arguments (out)
!        integer, intent(out) :: isize_nuc         ! size bin into which new particles go
!        real(r8), intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mol-air)
!        real(r8), intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mol/mol-air)
!        real(r8), intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mol/mol-air)
!        real(r8), intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mol/mol-air)
!        real(r8), intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mol/mol-air)
!                                                  ! aerosol changes are > 0; gas changes are < 0
!        real(r8), intent(out) :: dens_nh4so4a     ! dry-density of the new nh4-so4 aerosol mass (kg/m3)
!----------------------------------------------------------------------


!   convert qnuma_del from (#/mol-air) to (#/kmol-air)
        qnuma_del = qnuma_del*1.0e3_r8
!   number nuc rate (#/kmol-air/s) from number nuc amt
        dndt_ait = qnuma_del/deltat
!   fraction of mass nuc going to so4
        tmpa = qso4a_del*specmw_so4_amode
        tmpb = tmpa + qnh4a_del*specmw_nh4_amode
        tmp_frso4 = max( tmpa, 1.0e-35_r8 )/max( tmpb, 1.0e-35_r8 )
!   mass nuc rate (kg/kmol-air/s or g/mol...) hhfrom mass nuc amts
        dmdt_ait = max( 0.0_r8, (tmpb/deltat) ) 

	dndt_aitsv1 = dndt_ait
	dmdt_aitsv1 = dmdt_ait
	dndt_aitsv2 = 0.0_r8
	dmdt_aitsv2 = 0.0_r8
	dndt_aitsv3 = 0.0_r8
	dmdt_aitsv3 = 0.0_r8
        tmpch1 = ' '
        tmpch2 = ' '

	if (dndt_ait < 1.0e2_r8) then
!   ignore newnuc if number rate < 100 #/kmol-air/s ~= 0.3 #/mg-air/d
            dndt_ait = 0.0_r8
            dmdt_ait = 0.0_r8
            tmpch1 = 'A'

	else
	    dndt_aitsv2 = dndt_ait
	    dmdt_aitsv2 = dmdt_ait
            tmpch1 = 'B'

!   mirage2 code checked for complete h2so4 depletion here,
!   but this is now done in mer07_veh02_nuc_mosaic_1box
	    mass1p = dmdt_ait/dndt_ait
	    dndt_aitsv3 = dndt_ait
	    dmdt_aitsv3 = dmdt_ait

!   apply particle size constraints
	    if (mass1p < mass1p_aitlo) then
!   reduce dndt to increase new particle size
		dndt_ait = dmdt_ait/mass1p_aitlo
                tmpch1 = 'C'
	    else if (mass1p > mass1p_aithi) then
!   reduce dmdt to decrease new particle size
		dmdt_ait = dndt_ait*mass1p_aithi
                tmpch1 = 'E'
	    end if
	end if

! *** apply adjustment factor to avoid unrealistically high
!     aitken number concentrations in mid and upper troposphere
!	adjust_factor = 0.5
!	dndt_ait = dndt_ait * adjust_factor
!	dmdt_ait = dmdt_ait * adjust_factor

!   set tendencies
	pdel_fac = pdel(i,k)/gravit

!   dso4dt_ait, dnh4dt_ait are (kmol/kmol-air/s)
        dso4dt_ait = dmdt_ait*tmp_frso4/specmw_so4_amode
        dnh4dt_ait = dmdt_ait*(1.0_r8 - tmp_frso4)/specmw_nh4_amode

	dqdt(i,k,l_h2so4) = -dso4dt_ait*(1.0_r8-cldx)
	qsrflx(i,l_h2so4,1) = qsrflx(i,l_h2so4,1) + dqdt(i,k,l_h2so4)*pdel_fac
	q(i,k,l_h2so4) = q(i,k,l_h2so4) + dqdt(i,k,l_h2so4)*deltat

	dqdt(i,k,lso4ait) = dso4dt_ait*(1.0_r8-cldx)
	qsrflx(i,lso4ait,1) = qsrflx(i,lso4ait,1) + dqdt(i,k,lso4ait)*pdel_fac
	q(i,k,lso4ait) = q(i,k,lso4ait) + dqdt(i,k,lso4ait)*deltat
	if (lnumait > 0) then
	    dqdt(i,k,lnumait) = dndt_ait*(1.0_r8-cldx)
	    qsrflx(i,lnumait,1) = qsrflx(i,lnumait,1)   &
	                        + dqdt(i,k,lnumait)*pdel_fac
	    q(i,k,lnumait) = q(i,k,lnumait) + dqdt(i,k,lnumait)*deltat
	end if

	if (( do_nh3 ) .and. (dnh4dt_ait > 0.0_r8)) then
	    dqdt(i,k,l_nh3) = -dnh4dt_ait*(1.0_r8-cldx)
	    qsrflx(i,l_nh3,1) = qsrflx(i,l_nh3,1) + dqdt(i,k,l_nh3)*pdel_fac
	    q(i,k,l_nh3) = q(i,k,l_nh3) + dqdt(i,k,l_nh3)*deltat

	    dqdt(i,k,lnh4ait) = dnh4dt_ait*(1.0_r8-cldx)
	    qsrflx(i,lnh4ait,1) = qsrflx(i,lnh4ait,1) + dqdt(i,k,lnh4ait)*pdel_fac
	    q(i,k,lnh4ait) = q(i,k,lnh4ait) + dqdt(i,k,lnh4ait)*deltat
	end if

!!   temporary diagnostic
!        if (ldiag3 > 0) then
!        if ((dndt_ait /= 0.0_r8) .or. (dmdt_ait /= 0.0_r8)) then
!           write(lun,'(3a,1x,i7,3i5,1p,5e12.4)')   &
!              'newnucxx', tmpch1, tmpch2, nstep, lchnk, i, k,   &
!              dndt_ait, dmdt_ait, cldx
!!          call endrun( 'modal_aero_newnuc_sub' )
!        end if
!        end if


!   diagnostic output start ----------------------------------------
!!$ 	if (ldiag4 > 0) then
!!$ 	if ((lonndx(i) == 37) .and. (latndx(i) == 23)) then
!!$ 	if ((k >= 24) .or. (mod(k,4) == 0)) then
!!$        write(lun,97010) nstep, latndx(i), lonndx(i), k, t(i,k), aircon
!!$        write(lun,97020) 'pmid, pdel                   ',   &
!!$                pmid(i,k), pdel(i,k)
!!$        write(lun,97030) 'qv,qvsw, cld, rh_av, rh_clr  ',   &
!!$                qv(i,k), qvswtr, cldx, relhumav, relhum
!!$        write(lun,97020) 'h2so4_cur, _pre, _av, nh3_cur',   &
!!$ 		qh2so4_cur, tmp_q2, qh2so4_avg, qnh3_cur
!!$        write(lun,97020) 'del_h2so4_gasprod, _aeruptk  ',   &
!!$ 		del_h2so4_gasprod(i,k), del_h2so4_aeruptk(i,k),   &
!!$ 		tmp_uptkrate*3600.0_r8
!!$        write(lun,97020) ' '
!!$        write(lun,97050) 'tmpch1, tmpch2               ', tmpch1, tmpch2
!!$        write(lun,97020) 'dndt_, dmdt_aitsv1           ',   &
!!$ 				 dndt_aitsv1, dmdt_aitsv1
!!$        write(lun,97020) 'dndt_, dmdt_aitsv2           ',   &
!!$ 				 dndt_aitsv2, dmdt_aitsv2
!!$        write(lun,97020) 'dndt_, dmdt_aitsv3           ',   &
!!$ 				 dndt_aitsv3, dmdt_aitsv3
!!$        write(lun,97020) 'dndt_, dmdt_ait              ',   &
!!$ 				 dndt_ait, dmdt_ait
!!$        write(lun,97020) 'dso4dt_, dnh4dt_ait          ',   &
!!$ 				 dso4dt_ait, dnh4dt_ait
!!$        write(lun,97020) 'qso4a_del, qh2so4_del        ',   &
!!$ 				 qso4a_del, qh2so4_del
!!$        write(lun,97020) 'qnh4a_del, qnh3_del          ',   &
!!$ 				 qnh4a_del, qnh3_del
!!$        write(lun,97020) 'dqdt(h2so4), (nh3)           ',   &
!!$ 		 dqdt(i,k,l_h2so4), dqdt(i,k,l_nh3) 
!!$        write(lun,97020) 'dqdt(so4a), (nh4a), (numa)   ',   &
!!$ 		 dqdt(i,k,lso4ait), dqdt(i,k,lnh4ait), dqdt(i,k,lnumait)
!!$ 
!!$ 	dpnuc = 0.0_r8
!!$ 	if (dndt_aitsv1 > 1.0e-5_r8) dpnuc = (6.0_r8*dmdt_aitsv1/   &
!!$ 			(pi*specdens_so4_amode*dndt_aitsv1))**0.3333333_r8
!!$        if (dpnuc > 0.0_r8) then
!!$        write(lun,97020) 'dpnuc,      dp_aitlo, _aithi ',   &
!!$ 			 dpnuc, dplom_mode(1), dphim_mode(1)
!!$        write(lun,97020) 'mass1p, mass1p_aitlo, _aithi ',   &
!!$ 			 mass1p, mass1p_aitlo, mass1p_aithi
!!$        end if
!!$ 
!!$ 97010  format( / 'NEWNUC nstep,lat,lon,k,tk,cair', i8, 3i4, f8.2, 1pe12.4 )
!!$ 97020  format( a, 1p, 6e12.4 )
!!$ 97030  format( a, 1p, 2e12.4, 0p, 5f10.6 )
!!$ 97040  format( 29x, 1p, 6e12.4 )
!!$ 97050  format( a, 2(3x,a) )
!!$        end if
!!$        end if
!!$        end if
!   diagnostic output end   ------------------------------------------


	end do main_i
	end do main_k


!   do history file column-tendency fields
	do l = loffset+1, pcnst
	    lmz = l - loffset
	    if ( .not. dotend(lmz) ) cycle

	    do i = 1, ncol
		qsrflx(i,lmz,1) = qsrflx(i,lmz,1)*(adv_mass(lmz)/mwdry)
	    end do
	    fieldname = trim(cnst_name(l)) // '_sfnnuc1'
	    call outfld( fieldname, qsrflx(:,lmz,1), pcols, lchnk )

!	    if (( masterproc ) .and. (nstep < 1)) &
!		write(lun,'(2(a,2x),1p,e11.3)') &
!		'modal_aero_newnuc_sub outfld', fieldname, adv_mass(lmz)
	end do ! l = ...


	return
!EOC
	end subroutine modal_aero_newnuc_sub



!----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine mer07_veh02_nuc_mosaic_1box(   &
           newnuc_method_flagaa, dtnuc, temp_in, rh_in, press_in,   &
           zm_in, pblh_in,   &
           qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate,   &
           mw_so4a_host,   &
           nsize, maxd_asize, dplom_sect, dphim_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a, ldiagaa )
!          qh2so4_del, qnh3_del, dens_nh4so4a )
          use mo_constants, only: rgas, &               ! Gas constant (J/K/kmol)
                                  avogad => avogadro    ! Avogadro's number (1/kmol)
          use physconst,    only: mw_so4a => mwso4, &   ! Molecular weight of sulfate
                                  mw_nh4a => mwnh4      ! Molecular weight of ammonium
!.......................................................................
!
! calculates new particle production from homogeneous nucleation
!    over timestep dtnuc, using nucleation rates from either
!    merikanto et al. (2007) h2so4-nh3-h2o ternary parameterization
!    vehkamaki et al. (2002) h2so4-h2o binary parameterization
!
! the new particles are "grown" to the lower-bound size of the host code's 
!    smallest size bin.  (this "growth" is somewhat ad hoc, and would not be
!    necessary if the host code's size bins extend down to ~1 nm.)
!
!    if the h2so4 and nh3 mass mixing ratios (mixrats) of the grown new 
!    particles exceed the current gas mixrats, the new particle production
!    is reduced so that the new particle mass mixrats match the gas mixrats.
!
!    the correction of kerminen and kulmala (2002) is applied to account
!    for loss of the new particles by coagulation as they are
!    growing to the "host code mininum size"
!
! revision history
!    coded by rc easter, pnnl, xx-apr-2007
!
! key routines called: subr ternary_nuc_napari
!
! references:
!    merikanto, j., i. napari, h. vehkamaki, t. anttila,
!     and m. kulmala, 2007, new parameterization of
!     sulfuric acid-ammonia-water ternary nucleation
!     rates at tropospheric conditions,
!       j. geophys. res., 112, d15207, doi:10.1029/2006jd0027977
!
!    vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!       c. timmreck, m. noppel and a. laaksonen, 2002,
!       an improved parameterization for sulfuric acid-water nucleation
!       rates for tropospheric and stratospheric conditions,
!       j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
!
!    kerminen, v., and m. kulmala, 2002,
!	analytical formulae connecting the "real" and the "apparent"
!	nucleation rate and the nuclei number concentration
!	for atmospheric nucleation events
!
!.......................................................................
      implicit none

! subr arguments (in)
        real(r8), intent(in) :: dtnuc             ! nucleation time step (s)
        real(r8), intent(in) :: temp_in           ! temperature, in k
        real(r8), intent(in) :: rh_in             ! relative humidity, as fraction
        real(r8), intent(in) :: press_in          ! air pressure (pa)
        real(r8), intent(in) :: zm_in             ! layer midpoint height (m)
        real(r8), intent(in) :: pblh_in           ! pbl height (m)

        real(r8), intent(in) :: qh2so4_cur, qh2so4_avg
                                                  ! gas h2so4 mixing ratios (mol/mol-air)
        real(r8), intent(in) :: qnh3_cur          ! gas nh3 mixing ratios (mol/mol-air)
             ! qxxx_cur = current value (after gas chem and condensation)
             ! qxxx_avg = estimated average value (for simultaneous source/sink calcs)
        real(r8), intent(in) :: h2so4_uptkrate    ! h2so4 uptake rate to aerosol (1/s)
        real(r8), intent(in) :: mw_so4a_host      ! mw of so4 aerosol in host code (g/mol)

        integer, intent(in) :: newnuc_method_flagaa     ! 1=merikanto et al (2007) ternary
                                                        ! 2=vehkamaki et al (2002) binary
        integer, intent(in) :: nsize                    ! number of aerosol size bins
        integer, intent(in) :: maxd_asize               ! dimension for dplom_sect, ...
        real(r8), intent(in) :: dplom_sect(maxd_asize)  ! dry diameter at lower bnd of bin (m)
        real(r8), intent(in) :: dphim_sect(maxd_asize)  ! dry diameter at upper bnd of bin (m)
        integer, intent(in) :: ldiagaa

! subr arguments (out)
        integer, intent(out) :: isize_nuc         ! size bin into which new particles go
        real(r8), intent(out) :: qnuma_del        ! change to aerosol number mixing ratio (#/mol-air)
        real(r8), intent(out) :: qso4a_del        ! change to aerosol so4 mixing ratio (mol/mol-air)
        real(r8), intent(out) :: qnh4a_del        ! change to aerosol nh4 mixing ratio (mol/mol-air)
        real(r8), intent(out) :: qh2so4_del       ! change to gas h2so4 mixing ratio (mol/mol-air)
        real(r8), intent(out) :: qnh3_del         ! change to gas nh3 mixing ratio (mol/mol-air)
                                                  ! aerosol changes are > 0; gas changes are < 0
        real(r8), intent(out) :: dens_nh4so4a     ! dry-density of the new nh4-so4 aerosol mass (kg/m3)

! subr arguments (out) passed via common block  
!    these are used to duplicate the outputs of yang zhang's original test driver
!    they are not really needed in wrf-chem
        real(r8) :: ratenuclt        ! j = ternary nucleation rate from napari param. (cm-3 s-1)
        real(r8) :: rateloge         ! ln (j)
        real(r8) :: cnum_h2so4       ! number of h2so4 molecules in the critical nucleus
        real(r8) :: cnum_nh3         ! number of nh3   molecules in the critical nucleus
        real(r8) :: cnum_tot         ! total number of molecules in the critical nucleus
        real(r8) :: radius_cluster   ! the radius of cluster (nm)


! local variables
        integer :: i
        integer :: igrow
        integer, save :: icase = 0, icase_reldiffmax = 0
!       integer, parameter :: ldiagaa = -1
        integer :: lun
        integer :: newnuc_method_flagaa2
        integer :: use_ternary_rate, use_binary_rate, do_pbl_rate

        real(r8), parameter :: onethird = 1.0_r8/3.0_r8

        real(r8), parameter :: accom_coef_h2so4 = 0.65_r8   ! accomodation coef for h2so4 conden

! dry densities (kg/m3) molecular weights of aerosol 
! ammsulf, ammbisulf, and sulfacid (from mosaic  dens_electrolyte values)
!       real(r8), parameter :: dens_ammsulf   = 1.769e3
!       real(r8), parameter :: dens_ammbisulf = 1.78e3
!       real(r8), parameter :: dens_sulfacid  = 1.841e3
! use following to match cam3 modal_aero densities
        real(r8), parameter :: dens_ammsulf   = 1.770e3_r8
        real(r8), parameter :: dens_ammbisulf = 1.770e3_r8
        real(r8), parameter :: dens_sulfacid  = 1.770e3_r8

! molecular weights (g/mol) of aerosol ammsulf, ammbisulf, and sulfacid
!    for ammbisulf and sulfacid, use 114 & 96 here rather than 115 & 98
!    because we don't keep track of aerosol hion mass
        real(r8), parameter :: mw_ammsulf   = 132.0_r8
        real(r8), parameter :: mw_ammbisulf = 114.0_r8
        real(r8), parameter :: mw_sulfacid  =  96.0_r8

        real(r8), save :: reldiffmax = 0.0_r8

        real(r8) cair                     ! dry-air molar density (mol/m3)
        real(r8) cs_prime_kk              ! kk2002 "cs_prime" parameter (1/m2)
        real(r8) cs_kk                    ! kk2002 "cs" parameter (1/s)
        real(r8) dens_part                ! "grown" single-particle dry density (kg/m3)
        real(r8) dfin_kk, dnuc_kk         ! kk2002 final/initial new particle wet diameter (nm)
        real(r8) dpdry_clus               ! critical cluster diameter (m)
        real(r8) dpdry_part               ! "grown" single-particle dry diameter (m)
        real(r8) tmpa, tmpb, tmpc, tmpe, tmpq
        real(r8) tmpa1, tmpb1
        real(r8) tmp_m1, tmp_m2, tmp_m3, tmp_n1, tmp_n2, tmp_n3
        real(r8) tmp_spd                  ! h2so4 vapor molecular speed (m/s)
        real(r8) factor_kk
        real(r8) fogas, foso4a, fonh4a, fonuma
        real(r8) freduce                  ! reduction factor applied to nucleation rate
                                          ! due to limited availability of h2so4 & nh3 gases
        real(r8) freducea, freduceb
        real(r8) gamma_kk                 ! kk2002 "gamma" parameter (nm2*m2/h)
        real(r8) gr_kk                    ! kk2002 "gr" parameter (nm/h)
        real(r8) kgaero_per_moleso4a      ! (kg dry aerosol)/(mol aerosol so4)
        real(r8) mass_part                ! "grown" single-particle dry mass (kg)
        real(r8) molenh4a_per_moleso4a    ! (mol aerosol nh4)/(mol aerosol so4)
        real(r8) nh3ppt, nh3ppt_bb        ! actual and bounded nh3 (ppt)
        real(r8) nu_kk                    ! kk2002 "nu" parameter (nm)
        real(r8) qmolnh4a_del_max         ! max production of aerosol nh4 over dtnuc (mol/mol-air)
        real(r8) qmolso4a_del_max         ! max production of aerosol so4 over dtnuc (mol/mol-air)
        real(r8) ratenuclt_bb             ! nucleation rate (#/m3/s)
        real(r8) ratenuclt_kk             ! nucleation rate after kk2002 adjustment (#/m3/s)
        real(r8) rh_bb                    ! bounded value of rh_in
        real(r8) so4vol_in                ! concentration of h2so4 for nucl. calc., molecules cm-3
        real(r8) so4vol_bb                ! bounded value of so4vol_in
        real(r8) temp_bb                  ! bounded value of temp_in
        real(r8) voldry_clus              ! critical-cluster dry volume (m3)
        real(r8) voldry_part              ! "grown" single-particle dry volume (m3)
        real(r8) wetvol_dryvol            ! grown particle (wet-volume)/(dry-volume)
        real(r8) wet_volfrac_so4a         ! grown particle (dry-volume-from-so4)/(wet-volume)



!
! if h2so4 vapor < qh2so4_cutoff
! exit with new particle formation = 0
!
        isize_nuc = 1
        qnuma_del = 0.0_r8
        qso4a_del = 0.0_r8
        qnh4a_del = 0.0_r8
        qh2so4_del = 0.0_r8
        qnh3_del = 0.0_r8
!       if (qh2so4_avg .le. qh2so4_cutoff) return   ! this no longer needed
!       if (qh2so4_cur .le. qh2so4_cutoff) return   ! this no longer needed

        if ((newnuc_method_flagaa /=  1) .and. &
            (newnuc_method_flagaa /=  2) .and. &
            (newnuc_method_flagaa /= 11) .and. &
            (newnuc_method_flagaa /= 12)) return


!
! make call to parameterization routine
!

! calc h2so4 in molecules/cm3 and nh3 in ppt, and prepare bounded inputs
        call mer07_veh02_nuc_mosaic_prepare_rates( &
           newnuc_method_flagaa, temp_in, rh_in, press_in, zm_in, pblh_in, qh2so4_avg, qnh3_cur, &
           cair, so4vol_in, nh3ppt, ratenuclt, rateloge, temp_bb, rh_bb, so4vol_bb, nh3ppt_bb, &
           newnuc_method_flagaa2, use_ternary_rate, use_binary_rate, do_pbl_rate )

        if (use_ternary_rate /= 0) then
! make call to merikanto ternary parameterization routine
! (when nh3ppt < 0.1, use binary param instead)
           call ternary_nuc_merik2007(   &
              temp_bb, rh_bb, so4vol_bb, nh3ppt_bb,   &
              rateloge,   &
              cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )

        else
! make call to vehkamaki binary parameterization routine
           if (use_binary_rate /= 0) then
              call binary_nuc_vehk2002(   &
                 temp_bb, rh_bb, so4vol_bb,   &
                 ratenuclt, rateloge,   &
                 cnum_h2so4, cnum_tot, radius_cluster )
           end if
           cnum_nh3 = 0.0_r8
        end if


! do boundary layer nuc
        if (do_pbl_rate /= 0) then
           call pbl_nuc_wang2008( so4vol_in,   &
              newnuc_method_flagaa, newnuc_method_flagaa2,   &
              ratenuclt, rateloge,   &
              cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )
        end if


! if nucleation rate is less than 1e-6 #/m3/s ~= 0.1 #/cm3/day,
! exit with new particle formation = 0
        if (rateloge  .le. -13.82_r8) return
!       if (ratenuclt .le. 1.0e-6) return
        ratenuclt = exp( rateloge )
        ratenuclt_bb = ratenuclt*1.0e6_r8

        if (ldiagaa <= 0) then
           call mer07_veh02_nuc_mosaic_finalize( &
              dtnuc, temp_in, rh_in, press_in, qh2so4_cur, qnh3_cur, h2so4_uptkrate, mw_so4a_host, &
              nsize, maxd_asize, dplom_sect, dphim_sect, cnum_h2so4, cnum_nh3, radius_cluster, so4vol_in, ratenuclt_bb, &
              isize_nuc, qnuma_del, qso4a_del, qnh4a_del, qh2so4_del, qnh3_del, dens_nh4so4a )
           return
        end if


! wet/dry volume ratio - use simple kohler approx for ammsulf/ammbisulf
        tmpa = max( 0.10_r8, min( 0.95_r8, rh_in ) )
        wetvol_dryvol = 1.0_r8 - 0.56_r8/log(tmpa)


! determine size bin into which the new particles go
! (probably it will always be bin #1, but ...)
        voldry_clus = ( max(cnum_h2so4,1.0_r8)*mw_so4a + cnum_nh3*mw_nh4a ) /   &
                      (1.0e3_r8*dens_sulfacid*avogad)
! correction when host code sulfate is really ammonium bisulfate/sulfate
        voldry_clus = voldry_clus * (mw_so4a_host/mw_so4a)
        dpdry_clus = (voldry_clus*6.0_r8/pi)**onethird

        isize_nuc = 1
        dpdry_part = dplom_sect(1)
        if (dpdry_clus <= dplom_sect(1)) then
           igrow = 1   ! need to clusters to larger size
        else if (dpdry_clus >= dphim_sect(nsize)) then
           igrow = 0
           isize_nuc = nsize
           dpdry_part = dphim_sect(nsize)
        else
           igrow = 0
           do i = 1, nsize
              if (dpdry_clus < dphim_sect(i)) then
                 isize_nuc = i
                 dpdry_part = dpdry_clus
                 dpdry_part = min( dpdry_part, dphim_sect(i) )
                 dpdry_part = max( dpdry_part, dplom_sect(i) )
                 exit
              end if
           end do
        end if
        voldry_part = (pi/6.0_r8)*(dpdry_part**3)


!
! determine composition and density of the "grown particles"
! the grown particles are assumed to be liquid
!    (since critical clusters contain water)
!    so any (nh4/so4) molar ratio between 0 and 2 is allowed
! assume that the grown particles will have 
!    (nh4/so4 molar ratio) = min( 2, (nh3/h2so4 gas molar ratio) )
!
        if (igrow .le. 0) then
! no "growing" so pure sulfuric acid
           tmp_n1 = 0.0_r8
           tmp_n2 = 0.0_r8
           tmp_n3 = 1.0_r8
        else if (qnh3_cur .ge. qh2so4_cur) then
! combination of ammonium sulfate and ammonium bisulfate
! tmp_n1 & tmp_n2 = mole fractions of the ammsulf & ammbisulf
           tmp_n1 = (qnh3_cur/qh2so4_cur) - 1.0_r8
           tmp_n1 = max( 0.0_r8, min( 1.0_r8, tmp_n1 ) )
           tmp_n2 = 1.0_r8 - tmp_n1
           tmp_n3 = 0.0_r8
        else
! combination of ammonium bisulfate and sulfuric acid
! tmp_n2 & tmp_n3 = mole fractions of the ammbisulf & sulfacid
           tmp_n1 = 0.0_r8
           tmp_n2 = (qnh3_cur/qh2so4_cur)
           tmp_n2 = max( 0.0_r8, min( 1.0_r8, tmp_n2 ) )
           tmp_n3 = 1.0_r8 - tmp_n2
	end if

        tmp_m1 = tmp_n1*mw_ammsulf
        tmp_m2 = tmp_n2*mw_ammbisulf
        tmp_m3 = tmp_n3*mw_sulfacid
        dens_part = (tmp_m1 + tmp_m2 + tmp_m3)/   &
           ((tmp_m1/dens_ammsulf) + (tmp_m2/dens_ammbisulf)   &
                                  + (tmp_m3/dens_sulfacid))
        dens_nh4so4a = dens_part
        mass_part  = voldry_part*dens_part 
! (mol aerosol nh4)/(mol aerosol so4)
        molenh4a_per_moleso4a = 2.0_r8*tmp_n1 + tmp_n2  
! (kg dry aerosol)/(mol aerosol so4)
        kgaero_per_moleso4a = 1.0e-3_r8*(tmp_m1 + tmp_m2 + tmp_m3)  
! correction when host code sulfate is really ammonium bisulfate/sulfate
        kgaero_per_moleso4a = kgaero_per_moleso4a * (mw_so4a_host/mw_so4a)

! fraction of wet volume due to so4a
        tmpb = 1.0_r8 + molenh4a_per_moleso4a*17.0_r8/98.0_r8
        wet_volfrac_so4a = 1.0_r8 / ( wetvol_dryvol * tmpb )


!
! calc kerminen & kulmala (2002) correction
!
        if (igrow <=  0) then
            factor_kk = 1.0_r8

        else
! "gr" parameter (nm/h) = condensation growth rate of new particles
! use kk2002 eqn 21 for h2so4 uptake, and correct for nh3 & h2o uptake
            tmp_spd = 14.7_r8*sqrt(temp_in)   ! h2so4 molecular speed (m/s)
            gr_kk = 3.0e-9_r8*tmp_spd*mw_sulfacid*so4vol_in/   &
                    (dens_part*wet_volfrac_so4a)

! "gamma" parameter (nm2/m2/h)
! use kk2002 eqn 22
!
! dfin_kk = wet diam (nm) of grown particle having dry dia = dpdry_part (m)
            dfin_kk = 1.0e9_r8 * dpdry_part * (wetvol_dryvol**onethird)
! dnuc_kk = wet diam (nm) of cluster
            dnuc_kk = 2.0_r8*radius_cluster
            dnuc_kk = max( dnuc_kk, 1.0_r8 )
! neglect (dmean/150)**0.048 factor, 
! which should be very close to 1.0 because of small exponent
            gamma_kk = 0.23_r8 * (dnuc_kk)**0.2_r8   &
                     * (dfin_kk/3.0_r8)**0.075_r8   &
                     * (dens_part*1.0e-3_r8)**(-0.33_r8)   &
                     * (temp_in/293.0_r8)**(-0.75_r8)

! "cs_prime parameter" (1/m2) 
! instead kk2002 eqn 3, use
!     cs_prime ~= tmpa / (4*pi*tmpb * h2so4_accom_coef)
! where
!     tmpa = -d(ln(h2so4))/dt by conden to particles   (1/h units)
!     tmpb = h2so4 vapor diffusivity (m2/h units)
! this approx is generally within a few percent of the cs_prime
!     calculated directly from eqn 2, 
!     which is acceptable, given overall uncertainties
! tmpa = -d(ln(h2so4))/dt by conden to particles   (1/h units)
            tmpa = h2so4_uptkrate * 3600.0_r8
            tmpa1 = tmpa
            tmpa = max( tmpa, 0.0_r8 )
! tmpb = h2so4 gas diffusivity (m2/s, then m2/h)
            tmpb = 6.7037e-6_r8 * (temp_in**0.75_r8) / cair
            tmpb1 = tmpb         ! m2/s
            tmpb = tmpb*3600.0_r8   ! m2/h
            cs_prime_kk = tmpa/(4.0_r8*pi*tmpb*accom_coef_h2so4)
            cs_kk = cs_prime_kk*4.0_r8*pi*tmpb1

! "nu" parameter (nm) -- kk2002 eqn 11
            nu_kk = gamma_kk*cs_prime_kk/gr_kk
! nucleation rate adjustment factor (--) -- kk2002 eqn 13
            factor_kk = exp( (nu_kk/dfin_kk) - (nu_kk/dnuc_kk) )

        end if
        ratenuclt_kk = ratenuclt_bb*factor_kk


! max production of aerosol dry mass (kg-aero/m3-air)
        tmpa = max( 0.0_r8, (ratenuclt_kk*dtnuc*mass_part) )
! max production of aerosol so4 (mol-so4a/mol-air)
        tmpe = tmpa/(kgaero_per_moleso4a*cair)
! max production of aerosol so4 (mol/mol-air)
! based on ratenuclt_kk and mass_part
        qmolso4a_del_max = tmpe

! check if max production exceeds available h2so4 vapor
        freducea = 1.0_r8
        if (qmolso4a_del_max .gt. qh2so4_cur) then
           freducea = qh2so4_cur/qmolso4a_del_max
        end if

! check if max production exceeds available nh3 vapor
        freduceb = 1.0_r8
        if (molenh4a_per_moleso4a .ge. 1.0e-10_r8) then
! max production of aerosol nh4 (ppm) based on ratenuclt_kk and mass_part
           qmolnh4a_del_max = qmolso4a_del_max*molenh4a_per_moleso4a
           if (qmolnh4a_del_max .gt. qnh3_cur) then
              freduceb = qnh3_cur/qmolnh4a_del_max
           end if
        end if
        freduce = min( freducea, freduceb )

! if adjusted nucleation rate is less than 1e-12 #/m3/s ~= 0.1 #/cm3/day,
! exit with new particle formation = 0
        if (freduce*ratenuclt_kk .le. 1.0e-12_r8) return


! note:  suppose that at this point, freduce < 1.0 (no gas-available 
!    constraints) and molenh4a_per_moleso4a < 2.0
! if the gas-available constraints is do to h2so4 availability,
!    then it would be possible to condense "additional" nh3 and have
!    (nh3/h2so4 gas molar ratio) < (nh4/so4 aerosol molar ratio) <= 2 
! one could do some additional calculations of 
!    dens_part & molenh4a_per_moleso4a to realize this
! however, the particle "growing" is a crude approximate way to get
!    the new particles to the host code's minimum particle size,
! are such refinements worth the effort?


! changes to h2so4 & nh3 gas (in mol/mol-air), limited by amounts available
        tmpa = 0.9999_r8
        qh2so4_del = min( tmpa*qh2so4_cur, freduce*qmolso4a_del_max )
        qnh3_del   = min( tmpa*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del

! changes to so4 & nh4 aerosol (in mol/mol-air)
        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del
! change to aerosol number (in #/mol-air)
        qnuma_del = 1.0e-3_r8*(qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part

! do the following (tmpa, tmpb, tmpc) calculations as a check
! max production of aerosol number (#/mol-air)
        tmpa = max( 0.0_r8, (ratenuclt_kk*dtnuc/cair) )
! adjusted production of aerosol number (#/mol-air)
        tmpb = tmpa*freduce
! relative difference from qnuma_del
        tmpc = (tmpb - qnuma_del)/max(tmpb, qnuma_del, 1.0e-35_r8)


!
! diagnostic output to fort.41
! (this should be commented-out or deleted in the wrf-chem version)
!
        if (ldiagaa <= 0) return

        icase = icase + 1
        if (abs(tmpc) .gt. abs(reldiffmax)) then
           reldiffmax = tmpc
           icase_reldiffmax = icase
        end if
!       do lun = 41, 51, 10
        do lun = 6, 6
!          write(lun,'(/)')
           write(lun,'(a,2i9,1p,e10.2)')   &
               'vehkam bin-nuc icase, icase_rdmax =',   &
               icase, icase_reldiffmax, reldiffmax
           if (freduceb .lt. freducea) then
              if (abs(freducea-freduceb) .gt.   &
                   3.0e-7_r8*max(freduceb,freducea)) write(lun,'(a,1p,2e15.7)')   &
                 'freducea, b =', freducea, freduceb
           end if
        end do

! output factors so that output matches that of ternucl03
!       fogas  = 1.0e6                     ! convert mol/mol-air to ppm
!       foso4a = 1.0e9*mw_so4a/mw_air      ! convert mol-so4a/mol-air to ug/kg-air
!       fonh4a = 1.0e9*mw_nh4a/mw_air      ! convert mol-nh4a/mol-air to ug/kg-air
!       fonuma = 1.0e3/mw_air              ! convert #/mol-air to #/kg-air
        fogas  = 1.0_r8
        foso4a = 1.0_r8
        fonh4a = 1.0_r8
        fonuma = 1.0_r8

!       do lun = 41, 51, 10
        do lun = 6, 6

        write(lun,'(a,2i5)') 'newnuc_method_flagaa/aa2',   &
           newnuc_method_flagaa, newnuc_method_flagaa2

        write(lun,9210)
        write(lun,9201) temp_in, rh_in,   &
           ratenuclt, 2.0_r8*radius_cluster*1.0e-7_r8, dpdry_part*1.0e2_r8,   &
           voldry_part*1.0e6_r8, float(igrow)
        write(lun,9215)
        write(lun,9201)   &
           qh2so4_avg*fogas, 0.0_r8,  &
           qh2so4_cur*fogas, qnh3_cur*fogas,  &
           qh2so4_del*fogas, qnh3_del*fogas,  &
           qso4a_del*foso4a, qnh4a_del*fonh4a

        write(lun,9220)
        write(lun,9201)   &
           dtnuc, dens_nh4so4a*1.0e-3_r8,   &
           (qnh3_cur/qh2so4_cur), molenh4a_per_moleso4a,   &
           qnuma_del*fonuma, tmpb*fonuma, tmpc, freduce

        end do

!       lun = 51
        lun = 6
        write(lun,9230)
        write(lun,9201)   &
           press_in, cair*1.0e-6_r8, so4vol_in,   &
           wet_volfrac_so4a, wetvol_dryvol, dens_part*1.0e-3_r8

        if (igrow > 0) then
        write(lun,9240)
        write(lun,9201)   &
           tmp_spd, gr_kk, dnuc_kk, dfin_kk,   &
           gamma_kk, tmpa1, tmpb1, cs_kk

        write(lun,9250)
        write(lun,9201)   &
           cs_prime_kk, nu_kk, factor_kk, ratenuclt,   &
           ratenuclt_kk*1.0e-6_r8
        end if

9201    format ( 1p, 40e10.2  )
9210    format (   &
        '      temp        rh',   &
        '   ratenuc  dia_clus ddry_part',   &
        ' vdry_part     igrow' )
9215    format (   &
        '  h2so4avg  h2so4pre',   &
        '  h2so4cur   nh3_cur',   &
        '  h2so4del   nh3_del',   &
        '  so4a_del  nh4a_del' )
9220    format (    &
        '     dtnuc    dens_a   nh/so g   nh/so a',   &
        '  numa_del  numa_dl2   reldiff   freduce' )
9230    format (   &
        '  press_in      cair so4_volin',   &
        ' wet_volfr wetv_dryv dens_part' )
9240    format (   &
        '   tmp_spd     gr_kk   dnuc_kk   dfin_kk',   &
        '  gamma_kk     tmpa1     tmpb1     cs_kk' )
9250    format (   &
        ' cs_pri_kk     nu_kk factor_kk ratenuclt',   &
        ' ratenu_kk' )


        return
        end subroutine mer07_veh02_nuc_mosaic_1box

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine mer07_veh02_nuc_mosaic_prepare_rates_select_impl()

   use cam_logfile, only: iulog
   use spmd_utils, only: masterproc

   implicit none

   character(len=48) :: impl_name
   integer :: status, n, i, code

   if (mer07_veh02_nuc_mosaic_prepare_rates_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('MER07_VEH02_NUC_MOSAIC_PREPARE_RATES_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      mer07_veh02_nuc_mosaic_prepare_rates_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      mer07_veh02_nuc_mosaic_prepare_rates_use_native_impl = .false.
   end if

   mer07_veh02_nuc_mosaic_prepare_rates_impl_selected = .true.

   if (masterproc) then
      if (mer07_veh02_nuc_mosaic_prepare_rates_use_native_impl) then
         write(iulog,*) 'mer07_veh02_nuc_mosaic_prepare_rates implementation = native'
      else
         write(iulog,*) 'mer07_veh02_nuc_mosaic_prepare_rates implementation = codon'
      end if
   end if

  end subroutine mer07_veh02_nuc_mosaic_prepare_rates_select_impl

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine mer07_veh02_nuc_mosaic_prepare_rates( newnuc_method_flagaa, temp_in, rh_in, press_in, zm_in, pblh_in, &
       qh2so4_avg, qnh3_cur, cair, so4vol_in, nh3ppt, ratenuclt, rateloge, temp_bb, rh_bb, so4vol_bb, nh3ppt_bb, &
       newnuc_method_flagaa2, use_ternary_rate, use_binary_rate, do_pbl_rate )

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
   use mo_constants, only: rgas, avogad => avogadro

   implicit none

   integer, intent(in) :: newnuc_method_flagaa
   real(r8), intent(in) :: temp_in, rh_in, press_in, zm_in, pblh_in
   real(r8), intent(in) :: qh2so4_avg, qnh3_cur
   real(r8), intent(out) :: cair, so4vol_in, nh3ppt, ratenuclt, rateloge
   real(r8), intent(out) :: temp_bb, rh_bb, so4vol_bb, nh3ppt_bb
   integer, intent(out) :: newnuc_method_flagaa2, use_ternary_rate, use_binary_rate, do_pbl_rate

   integer(c_int64_t), target :: newnuc_method_flagaa2_work
   integer(c_int64_t), target :: use_ternary_rate_work, use_binary_rate_work, do_pbl_rate_work
   real(c_double), target :: cair_work, so4vol_in_work, nh3ppt_work
   real(c_double), target :: ratenuclt_work, rateloge_work
   real(c_double), target :: temp_bb_work, rh_bb_work, so4vol_bb_work, nh3ppt_bb_work

   interface
      subroutine mer07_veh02_nuc_mosaic_prepare_rates_codon( newnuc_method_flagaa_c, temp_in_c, rh_in_c, press_in_c, zm_in_c, &
           pblh_in_c, qh2so4_avg_c, qnh3_cur_c, rgas_c, avogad_c, cair_p, so4vol_in_p, nh3ppt_p, ratenuclt_p, rateloge_p, &
           temp_bb_p, rh_bb_p, so4vol_bb_p, nh3ppt_bb_p, newnuc_method_flagaa2_p, use_ternary_rate_p, use_binary_rate_p, &
           do_pbl_rate_p ) bind(c, name="mer07_veh02_nuc_mosaic_prepare_rates_codon")
        use iso_c_binding, only: c_double, c_int64_t, c_ptr
        integer(c_int64_t), value :: newnuc_method_flagaa_c
        real(c_double), value :: temp_in_c, rh_in_c, press_in_c, zm_in_c, pblh_in_c, qh2so4_avg_c, qnh3_cur_c
        real(c_double), value :: rgas_c, avogad_c
        type(c_ptr), value :: cair_p, so4vol_in_p, nh3ppt_p, ratenuclt_p, rateloge_p
        type(c_ptr), value :: temp_bb_p, rh_bb_p, so4vol_bb_p, nh3ppt_bb_p
        type(c_ptr), value :: newnuc_method_flagaa2_p, use_ternary_rate_p, use_binary_rate_p, do_pbl_rate_p
      end subroutine mer07_veh02_nuc_mosaic_prepare_rates_codon
   end interface

   call mer07_veh02_nuc_mosaic_prepare_rates_select_impl()

   if (mer07_veh02_nuc_mosaic_prepare_rates_use_native_impl) then
      call mer07_veh02_nuc_mosaic_prepare_rates_native( newnuc_method_flagaa, temp_in, rh_in, press_in, zm_in, pblh_in, &
           qh2so4_avg, qnh3_cur, cair, so4vol_in, nh3ppt, ratenuclt, rateloge, temp_bb, rh_bb, so4vol_bb, nh3ppt_bb, &
           newnuc_method_flagaa2, use_ternary_rate, use_binary_rate, do_pbl_rate )
      return
   end if

   call mer07_veh02_nuc_mosaic_prepare_rates_codon( &
        int(newnuc_method_flagaa, c_int64_t), real(temp_in, c_double), real(rh_in, c_double), real(press_in, c_double), &
        real(zm_in, c_double), real(pblh_in, c_double), real(qh2so4_avg, c_double), real(qnh3_cur, c_double), &
        real(rgas, c_double), real(avogad, c_double), c_loc(cair_work), c_loc(so4vol_in_work), c_loc(nh3ppt_work), &
        c_loc(ratenuclt_work), c_loc(rateloge_work), c_loc(temp_bb_work), c_loc(rh_bb_work), c_loc(so4vol_bb_work), &
        c_loc(nh3ppt_bb_work), c_loc(newnuc_method_flagaa2_work), c_loc(use_ternary_rate_work), c_loc(use_binary_rate_work), &
        c_loc(do_pbl_rate_work) )

   cair = real(cair_work, r8)
   so4vol_in = real(so4vol_in_work, r8)
   nh3ppt = real(nh3ppt_work, r8)
   ratenuclt = real(ratenuclt_work, r8)
   rateloge = real(rateloge_work, r8)
   temp_bb = real(temp_bb_work, r8)
   rh_bb = real(rh_bb_work, r8)
   so4vol_bb = real(so4vol_bb_work, r8)
   nh3ppt_bb = real(nh3ppt_bb_work, r8)
   newnuc_method_flagaa2 = int(newnuc_method_flagaa2_work, kind(newnuc_method_flagaa2))
   use_ternary_rate = int(use_ternary_rate_work, kind(use_ternary_rate))
   use_binary_rate = int(use_binary_rate_work, kind(use_binary_rate))
   do_pbl_rate = int(do_pbl_rate_work, kind(do_pbl_rate))

  end subroutine mer07_veh02_nuc_mosaic_prepare_rates

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine mer07_veh02_nuc_mosaic_prepare_rates_native( newnuc_method_flagaa, temp_in, rh_in, press_in, zm_in, pblh_in, &
       qh2so4_avg, qnh3_cur, cair, so4vol_in, nh3ppt, ratenuclt, rateloge, temp_bb, rh_bb, so4vol_bb, nh3ppt_bb, &
       newnuc_method_flagaa2, use_ternary_rate, use_binary_rate, do_pbl_rate )

   use mo_constants, only: rgas, avogad => avogadro

   implicit none

   integer, intent(in) :: newnuc_method_flagaa
   real(r8), intent(in) :: temp_in, rh_in, press_in, zm_in, pblh_in
   real(r8), intent(in) :: qh2so4_avg, qnh3_cur
   real(r8), intent(out) :: cair, so4vol_in, nh3ppt, ratenuclt, rateloge
   real(r8), intent(out) :: temp_bb, rh_bb, so4vol_bb, nh3ppt_bb
   integer, intent(out) :: newnuc_method_flagaa2, use_ternary_rate, use_binary_rate, do_pbl_rate

   cair = press_in/(temp_in*rgas)
   so4vol_in = qh2so4_avg * cair * avogad * 1.0e-6_r8
   nh3ppt = qnh3_cur * 1.0e12_r8
   ratenuclt = 1.0e-38_r8
   rateloge = log( ratenuclt )
   temp_bb = 0.0_r8
   rh_bb = 0.0_r8
   so4vol_bb = 0.0_r8
   nh3ppt_bb = 0.0_r8
   use_ternary_rate = 0
   use_binary_rate = 0
   do_pbl_rate = 0

   if ( (newnuc_method_flagaa /= 2) .and. (nh3ppt >= 0.1_r8) ) then
      if (so4vol_in >= 5.0e4_r8) then
         temp_bb = max( 235.0_r8, min( 295.0_r8, temp_in ) )
         rh_bb = max( 0.05_r8, min( 0.95_r8, rh_in ) )
         so4vol_bb = max( 5.0e4_r8, min( 1.0e9_r8, so4vol_in ) )
         nh3ppt_bb = max( 0.1_r8, min( 1.0e3_r8, nh3ppt ) )
         use_ternary_rate = 1
      end if
      newnuc_method_flagaa2 = 1
   else
      if (so4vol_in >= 1.0e4_r8) then
         temp_bb = max( 230.15_r8, min( 305.15_r8, temp_in ) )
         rh_bb = max( 1.0e-4_r8, min( 1.0_r8, rh_in ) )
         so4vol_bb = max( 1.0e4_r8, min( 1.0e11_r8, so4vol_in ) )
         use_binary_rate = 1
      end if
      newnuc_method_flagaa2 = 2
   end if

   if ((newnuc_method_flagaa == 11) .or. (newnuc_method_flagaa == 12)) then
      if (zm_in <= max(pblh_in, 100.0_r8)) then
         do_pbl_rate = 1
      end if
   end if

  end subroutine mer07_veh02_nuc_mosaic_prepare_rates_native

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine mer07_veh02_nuc_mosaic_finalize_select_impl()

   use cam_logfile, only: iulog
   use spmd_utils, only: masterproc

   implicit none

   character(len=48) :: impl_name
   integer :: status, n, i, code

   if (mer07_veh02_nuc_mosaic_finalize_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('MER07_VEH02_NUC_MOSAIC_FINALIZE_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      mer07_veh02_nuc_mosaic_finalize_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      mer07_veh02_nuc_mosaic_finalize_use_native_impl = .false.
   end if

   mer07_veh02_nuc_mosaic_finalize_impl_selected = .true.

   if (masterproc) then
      if (mer07_veh02_nuc_mosaic_finalize_use_native_impl) then
         write(iulog,*) 'mer07_veh02_nuc_mosaic_finalize implementation = native'
      else
         write(iulog,*) 'mer07_veh02_nuc_mosaic_finalize implementation = codon'
      end if
   end if

  end subroutine mer07_veh02_nuc_mosaic_finalize_select_impl

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine mer07_veh02_nuc_mosaic_finalize( dtnuc, temp_in, rh_in, press_in, qh2so4_cur, qnh3_cur, h2so4_uptkrate, &
       mw_so4a_host, nsize, maxd_asize, dplom_sect, dphim_sect, cnum_h2so4, cnum_nh3, radius_cluster, so4vol_in, &
       ratenuclt_bb, isize_nuc, qnuma_del, qso4a_del, qnh4a_del, qh2so4_del, qnh3_del, dens_nh4so4a )

   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
   use mo_constants, only: pi, rgas, avogad => avogadro
   use physconst, only: mw_so4a => mwso4, mw_nh4a => mwnh4

   implicit none

   real(r8), intent(in) :: dtnuc, temp_in, rh_in, press_in
   real(r8), intent(in) :: qh2so4_cur, qnh3_cur, h2so4_uptkrate, mw_so4a_host
   integer, intent(in) :: nsize, maxd_asize
   real(r8), intent(in), target :: dplom_sect(maxd_asize), dphim_sect(maxd_asize)
   real(r8), intent(in) :: cnum_h2so4, cnum_nh3, radius_cluster, so4vol_in, ratenuclt_bb
   integer, intent(inout) :: isize_nuc
   real(r8), intent(inout) :: qnuma_del, qso4a_del, qnh4a_del, qh2so4_del, qnh3_del, dens_nh4so4a

   integer(c_int64_t), target :: isize_nuc_work
   real(c_double), target :: qnuma_del_work, qso4a_del_work, qnh4a_del_work
   real(c_double), target :: qh2so4_del_work, qnh3_del_work, dens_nh4so4a_work

   interface
      subroutine mer07_veh02_nuc_mosaic_finalize_codon( dtnuc_c, temp_in_c, rh_in_c, press_in_c, qh2so4_cur_c, qnh3_cur_c, &
           h2so4_uptkrate_c, mw_so4a_host_c, nsize_c, dplom_sect_p, dphim_sect_p, cnum_h2so4_c, cnum_nh3_c, radius_cluster_c, &
           so4vol_in_c, ratenuclt_bb_c, pi_c, rgas_c, avogad_c, mw_so4a_c, mw_nh4a_c, isize_nuc_p, qnuma_del_p, qso4a_del_p, &
           qnh4a_del_p, qh2so4_del_p, qnh3_del_p, dens_nh4so4a_p ) bind(c, name="mer07_veh02_nuc_mosaic_finalize_codon")
        use iso_c_binding, only: c_double, c_int64_t, c_ptr
        real(c_double), value :: dtnuc_c, temp_in_c, rh_in_c, press_in_c, qh2so4_cur_c, qnh3_cur_c
        real(c_double), value :: h2so4_uptkrate_c, mw_so4a_host_c, cnum_h2so4_c, cnum_nh3_c, radius_cluster_c
        real(c_double), value :: so4vol_in_c, ratenuclt_bb_c, pi_c, rgas_c, avogad_c, mw_so4a_c, mw_nh4a_c
        integer(c_int64_t), value :: nsize_c
        type(c_ptr), value :: dplom_sect_p, dphim_sect_p, isize_nuc_p, qnuma_del_p, qso4a_del_p, qnh4a_del_p
        type(c_ptr), value :: qh2so4_del_p, qnh3_del_p, dens_nh4so4a_p
      end subroutine mer07_veh02_nuc_mosaic_finalize_codon
   end interface

   call mer07_veh02_nuc_mosaic_finalize_select_impl()

   if (mer07_veh02_nuc_mosaic_finalize_use_native_impl) then
      call mer07_veh02_nuc_mosaic_finalize_native( dtnuc, temp_in, rh_in, press_in, qh2so4_cur, qnh3_cur, h2so4_uptkrate, &
           mw_so4a_host, nsize, maxd_asize, dplom_sect, dphim_sect, cnum_h2so4, cnum_nh3, radius_cluster, so4vol_in, ratenuclt_bb, &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del, qh2so4_del, qnh3_del, dens_nh4so4a )
      return
   end if

   isize_nuc_work = int(isize_nuc, c_int64_t)
   qnuma_del_work = real(qnuma_del, c_double)
   qso4a_del_work = real(qso4a_del, c_double)
   qnh4a_del_work = real(qnh4a_del, c_double)
   qh2so4_del_work = real(qh2so4_del, c_double)
   qnh3_del_work = real(qnh3_del, c_double)

   call mer07_veh02_nuc_mosaic_finalize_codon( &
        real(dtnuc, c_double), real(temp_in, c_double), real(rh_in, c_double), &
        real(press_in, c_double), real(qh2so4_cur, c_double), real(qnh3_cur, c_double), &
        real(h2so4_uptkrate, c_double), real(mw_so4a_host, c_double), int(nsize, c_int64_t), &
        c_loc(dplom_sect(1)), c_loc(dphim_sect(1)), real(cnum_h2so4, c_double), &
        real(cnum_nh3, c_double), real(radius_cluster, c_double), real(so4vol_in, c_double), &
        real(ratenuclt_bb, c_double), real(pi, c_double), real(rgas, c_double), &
        real(avogad, c_double), real(mw_so4a, c_double), real(mw_nh4a, c_double), &
        c_loc(isize_nuc_work), c_loc(qnuma_del_work), c_loc(qso4a_del_work), &
        c_loc(qnh4a_del_work), c_loc(qh2so4_del_work), c_loc(qnh3_del_work), &
        c_loc(dens_nh4so4a_work) )

   isize_nuc = int(isize_nuc_work, kind(isize_nuc))
   qnuma_del = real(qnuma_del_work, r8)
   qso4a_del = real(qso4a_del_work, r8)
   qnh4a_del = real(qnh4a_del_work, r8)
   qh2so4_del = real(qh2so4_del_work, r8)
   qnh3_del = real(qnh3_del_work, r8)
   dens_nh4so4a = real(dens_nh4so4a_work, r8)

  end subroutine mer07_veh02_nuc_mosaic_finalize

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine mer07_veh02_nuc_mosaic_finalize_native( dtnuc, temp_in, rh_in, press_in, qh2so4_cur, qnh3_cur, h2so4_uptkrate, &
       mw_so4a_host, nsize, maxd_asize, dplom_sect, dphim_sect, cnum_h2so4, cnum_nh3, radius_cluster, so4vol_in, ratenuclt_bb, &
       isize_nuc, qnuma_del, qso4a_del, qnh4a_del, qh2so4_del, qnh3_del, dens_nh4so4a )

   use mo_constants, only: pi, rgas, avogad => avogadro
   use physconst, only: mw_so4a => mwso4, mw_nh4a => mwnh4

   implicit none

   real(r8), intent(in) :: dtnuc, temp_in, rh_in, press_in
   real(r8), intent(in) :: qh2so4_cur, qnh3_cur, h2so4_uptkrate, mw_so4a_host
   integer, intent(in) :: nsize, maxd_asize
   real(r8), intent(in) :: dplom_sect(maxd_asize), dphim_sect(maxd_asize)
   real(r8), intent(in) :: cnum_h2so4, cnum_nh3, radius_cluster, so4vol_in, ratenuclt_bb
   integer, intent(inout) :: isize_nuc
   real(r8), intent(inout) :: qnuma_del, qso4a_del, qnh4a_del, qh2so4_del, qnh3_del, dens_nh4so4a

   integer :: i
   integer :: igrow
   real(r8), parameter :: onethird = 1.0_r8/3.0_r8
   real(r8), parameter :: accom_coef_h2so4 = 0.65_r8
   real(r8), parameter :: dens_ammsulf   = 1.770e3_r8
   real(r8), parameter :: dens_ammbisulf = 1.770e3_r8
   real(r8), parameter :: dens_sulfacid  = 1.770e3_r8
   real(r8), parameter :: mw_ammsulf   = 132.0_r8
   real(r8), parameter :: mw_ammbisulf = 114.0_r8
   real(r8), parameter :: mw_sulfacid  = 96.0_r8
   real(r8) :: cair, cs_prime_kk, dens_part, dfin_kk, dnuc_kk
   real(r8) :: dpdry_clus, dpdry_part, factor_kk, freduce, freducea, freduceb
   real(r8) :: gamma_kk, gr_kk, kgaero_per_moleso4a, mass_part, molenh4a_per_moleso4a
   real(r8) :: nu_kk, qmolnh4a_del_max, qmolso4a_del_max, ratenuclt_kk
   real(r8) :: tmp_m1, tmp_m2, tmp_m3, tmp_n1, tmp_n2, tmp_n3
   real(r8) :: tmp_spd, tmpa, tmpb, tmpe, voldry_clus, voldry_part
   real(r8) :: wet_volfrac_so4a, wetvol_dryvol

   cair = press_in/(temp_in*rgas)

   tmpa = max( 0.10_r8, min( 0.95_r8, rh_in ) )
   wetvol_dryvol = 1.0_r8 - 0.56_r8/log(tmpa)

   voldry_clus = ( max(cnum_h2so4,1.0_r8)*mw_so4a + cnum_nh3*mw_nh4a ) / &
                 (1.0e3_r8*dens_sulfacid*avogad)
   voldry_clus = voldry_clus * (mw_so4a_host/mw_so4a)
   dpdry_clus = (voldry_clus*6.0_r8/pi)**onethird

   dpdry_part = dplom_sect(1)
   if (dpdry_clus <= dplom_sect(1)) then
      igrow = 1
   else if (dpdry_clus >= dphim_sect(nsize)) then
      igrow = 0
      isize_nuc = nsize
      dpdry_part = dphim_sect(nsize)
   else
      igrow = 0
      do i = 1, nsize
         if (dpdry_clus < dphim_sect(i)) then
            isize_nuc = i
            dpdry_part = dpdry_clus
            dpdry_part = min( dpdry_part, dphim_sect(i) )
            dpdry_part = max( dpdry_part, dplom_sect(i) )
            exit
         end if
      end do
   end if
   voldry_part = (pi/6.0_r8)*(dpdry_part**3)

   if (igrow .le. 0) then
      tmp_n1 = 0.0_r8
      tmp_n2 = 0.0_r8
      tmp_n3 = 1.0_r8
   else if (qnh3_cur .ge. qh2so4_cur) then
      tmp_n1 = (qnh3_cur/qh2so4_cur) - 1.0_r8
      tmp_n1 = max( 0.0_r8, min( 1.0_r8, tmp_n1 ) )
      tmp_n2 = 1.0_r8 - tmp_n1
      tmp_n3 = 0.0_r8
   else
      tmp_n1 = 0.0_r8
      tmp_n2 = (qnh3_cur/qh2so4_cur)
      tmp_n2 = max( 0.0_r8, min( 1.0_r8, tmp_n2 ) )
      tmp_n3 = 1.0_r8 - tmp_n2
   end if

   tmp_m1 = tmp_n1*mw_ammsulf
   tmp_m2 = tmp_n2*mw_ammbisulf
   tmp_m3 = tmp_n3*mw_sulfacid
   dens_part = (tmp_m1 + tmp_m2 + tmp_m3)/((tmp_m1/dens_ammsulf) + (tmp_m2/dens_ammbisulf) + (tmp_m3/dens_sulfacid))
   dens_nh4so4a = dens_part
   mass_part = voldry_part*dens_part
   molenh4a_per_moleso4a = 2.0_r8*tmp_n1 + tmp_n2
   kgaero_per_moleso4a = 1.0e-3_r8*(tmp_m1 + tmp_m2 + tmp_m3)
   kgaero_per_moleso4a = kgaero_per_moleso4a * (mw_so4a_host/mw_so4a)

   tmpb = 1.0_r8 + molenh4a_per_moleso4a*17.0_r8/98.0_r8
   wet_volfrac_so4a = 1.0_r8 / ( wetvol_dryvol * tmpb )

   if (igrow <= 0) then
      factor_kk = 1.0_r8
   else
      tmp_spd = 14.7_r8*sqrt(temp_in)
      gr_kk = 3.0e-9_r8*tmp_spd*mw_sulfacid*so4vol_in/(dens_part*wet_volfrac_so4a)
      dfin_kk = 1.0e9_r8 * dpdry_part * (wetvol_dryvol**onethird)
      dnuc_kk = 2.0_r8*radius_cluster
      dnuc_kk = max( dnuc_kk, 1.0_r8 )
      gamma_kk = 0.23_r8 * (dnuc_kk)**0.2_r8 &
               * (dfin_kk/3.0_r8)**0.075_r8 &
               * (dens_part*1.0e-3_r8)**(-0.33_r8) &
               * (temp_in/293.0_r8)**(-0.75_r8)
      tmpa = h2so4_uptkrate * 3600.0_r8
      tmpa = max( tmpa, 0.0_r8 )
      tmpb = 6.7037e-6_r8 * (temp_in**0.75_r8) / cair
      tmpb = tmpb*3600.0_r8
      cs_prime_kk = tmpa/(4.0_r8*pi*tmpb*accom_coef_h2so4)
      nu_kk = gamma_kk*cs_prime_kk/gr_kk
      factor_kk = exp( (nu_kk/dfin_kk) - (nu_kk/dnuc_kk) )
   end if
   ratenuclt_kk = ratenuclt_bb*factor_kk

   tmpa = max( 0.0_r8, (ratenuclt_kk*dtnuc*mass_part) )
   tmpe = tmpa/(kgaero_per_moleso4a*cair)
   qmolso4a_del_max = tmpe

   freducea = 1.0_r8
   if (qmolso4a_del_max .gt. qh2so4_cur) then
      freducea = qh2so4_cur/qmolso4a_del_max
   end if

   freduceb = 1.0_r8
   if (molenh4a_per_moleso4a .ge. 1.0e-10_r8) then
      qmolnh4a_del_max = qmolso4a_del_max*molenh4a_per_moleso4a
      if (qmolnh4a_del_max .gt. qnh3_cur) then
         freduceb = qnh3_cur/qmolnh4a_del_max
      end if
   end if
   freduce = min( freducea, freduceb )

   if (freduce*ratenuclt_kk .le. 1.0e-12_r8) return

   tmpa = 0.9999_r8
   qh2so4_del = min( tmpa*qh2so4_cur, freduce*qmolso4a_del_max )
   qnh3_del = min( tmpa*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
   qh2so4_del = -qh2so4_del
   qnh3_del = -qnh3_del
   qso4a_del = -qh2so4_del
   qnh4a_del = -qnh3_del
   qnuma_del = 1.0e-3_r8*(qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part

  end subroutine mer07_veh02_nuc_mosaic_finalize_native



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine pbl_nuc_wang2008_select_impl()

   use cam_logfile, only: iulog
   use spmd_utils, only: masterproc

   implicit none

   character(len=48) :: impl_name
   integer :: status, n, i, code

   if (pbl_nuc_wang2008_impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('PBL_NUC_WANG2008_IMPL', value=impl_name, length=n, status=status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      pbl_nuc_wang2008_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      pbl_nuc_wang2008_use_native_impl = .false.
   end if

   pbl_nuc_wang2008_impl_selected = .true.

   if (masterproc) then
      if (pbl_nuc_wang2008_use_native_impl) then
         write(iulog,*) 'pbl_nuc_wang2008 implementation = native'
      else
         write(iulog,*) 'pbl_nuc_wang2008 implementation = codon'
      end if
   end if

  end subroutine pbl_nuc_wang2008_select_impl

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine pbl_nuc_wang2008( so4vol,   &
            newnuc_method_flagaa, newnuc_method_flagaa2,   &
            ratenucl, rateloge,   &
            cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )

        use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

        implicit none

! subr arguments (in)
        real(r8), intent(in) :: so4vol
        integer, intent(in)  :: newnuc_method_flagaa

! subr arguments (inout)
        integer, intent(inout)  :: newnuc_method_flagaa2
        real(r8), intent(inout) :: ratenucl
        real(r8), intent(inout) :: rateloge
        real(r8), intent(inout) :: cnum_tot
        real(r8), intent(inout) :: cnum_h2so4
        real(r8), intent(inout) :: cnum_nh3
        real(r8), intent(inout) :: radius_cluster

        integer(c_int64_t), target :: newnuc_method_flagaa2_work
        real(c_double), target :: ratenucl_work, rateloge_work
        real(c_double), target :: cnum_tot_work, cnum_h2so4_work
        real(c_double), target :: cnum_nh3_work, radius_cluster_work

        interface
           subroutine pbl_nuc_wang2008_codon(so4vol_c, newnuc_method_flagaa_c, newnuc_method_flagaa2_p, ratenucl_p, rateloge_p, &
                cnum_tot_p, cnum_h2so4_p, cnum_nh3_p, radius_cluster_p) bind(c, name="pbl_nuc_wang2008_codon")
             use iso_c_binding, only: c_double, c_int64_t, c_ptr
             real(c_double), value :: so4vol_c
             integer(c_int64_t), value :: newnuc_method_flagaa_c
             type(c_ptr), value :: newnuc_method_flagaa2_p, ratenucl_p, rateloge_p
             type(c_ptr), value :: cnum_tot_p, cnum_h2so4_p, cnum_nh3_p, radius_cluster_p
           end subroutine pbl_nuc_wang2008_codon
        end interface

        call pbl_nuc_wang2008_select_impl()

        if (pbl_nuc_wang2008_use_native_impl) then
           call pbl_nuc_wang2008_native( &
                so4vol, newnuc_method_flagaa, newnuc_method_flagaa2, ratenucl, rateloge, cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster)
           return
        end if

        newnuc_method_flagaa2_work = int(newnuc_method_flagaa2, c_int64_t)
        ratenucl_work = real(ratenucl, c_double)
        rateloge_work = real(rateloge, c_double)
        cnum_tot_work = real(cnum_tot, c_double)
        cnum_h2so4_work = real(cnum_h2so4, c_double)
        cnum_nh3_work = real(cnum_nh3, c_double)
        radius_cluster_work = real(radius_cluster, c_double)

        call pbl_nuc_wang2008_codon( &
             real(so4vol, c_double), int(newnuc_method_flagaa, c_int64_t), c_loc(newnuc_method_flagaa2_work), c_loc(ratenucl_work), &
             c_loc(rateloge_work), c_loc(cnum_tot_work), c_loc(cnum_h2so4_work), c_loc(cnum_nh3_work), c_loc(radius_cluster_work) &
        )

        newnuc_method_flagaa2 = int(newnuc_method_flagaa2_work, kind(newnuc_method_flagaa2))
        ratenucl = real(ratenucl_work, r8)
        rateloge = real(rateloge_work, r8)
        cnum_tot = real(cnum_tot_work, r8)
        cnum_h2so4 = real(cnum_h2so4_work, r8)
        cnum_nh3 = real(cnum_nh3_work, r8)
        radius_cluster = real(radius_cluster_work, r8)

        end subroutine pbl_nuc_wang2008

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine pbl_nuc_wang2008_native( so4vol,   &
            newnuc_method_flagaa, newnuc_method_flagaa2,   &
            ratenucl, rateloge,   &
            cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster )
!
! calculates boundary nucleation nucleation rate
! using the first or second-order parameterization in  
!     wang, m., and j.e. penner, 2008,
!        aerosol indirect forcing in a global model with particle nucleation,
!        atmos. chem. phys. discuss., 8, 13943-13998
!
        implicit none

! subr arguments (in)
        real(r8), intent(in) :: so4vol            ! concentration of h2so4 (molecules cm-3)
        integer, intent(in)  :: newnuc_method_flagaa  
                                ! [11,12] value selects [first,second]-order parameterization

! subr arguments (inout)
        integer, intent(inout)  :: newnuc_method_flagaa2
        real(r8), intent(inout) :: ratenucl         ! binary nucleation rate, j (# cm-3 s-1)
        real(r8), intent(inout) :: rateloge         ! log( ratenucl )

        real(r8), intent(inout) :: cnum_tot         ! total number of molecules
                                                    ! in the critical nucleus
        real(r8), intent(inout) :: cnum_h2so4       ! number of h2so4 molecules
        real(r8), intent(inout) :: cnum_nh3         ! number of nh3 molecules
        real(r8), intent(inout) :: radius_cluster   ! the radius of cluster (nm)


! local variables
        real(r8) :: tmp_diam, tmp_mass, tmp_volu
        real(r8) :: tmp_rateloge, tmp_ratenucl

! executable


! nucleation rate
        if (newnuc_method_flagaa == 11) then
           tmp_ratenucl = 1.0e-6_r8 * so4vol
        else if (newnuc_method_flagaa == 12) then
           tmp_ratenucl = 1.0e-12_r8 * (so4vol**2)
        else
           return
        end if
        tmp_rateloge = log( tmp_ratenucl )

! exit if pbl nuc rate is lower than (incoming) ternary/binary rate
        if (tmp_rateloge <= rateloge) return

        rateloge = tmp_rateloge
        ratenucl = tmp_ratenucl
        newnuc_method_flagaa2 = newnuc_method_flagaa

! following wang 2002, assume fresh nuclei are 1 nm diameter
!    subsequent code will "grow" them to aitken mode size
        radius_cluster = 0.5_r8

! assume fresh nuclei are pure h2so4
!    since aitken size >> initial size, the initial composition 
!    has very little impact on the results
        tmp_diam = radius_cluster * 2.0e-7_r8   ! diameter in cm
        tmp_volu = (tmp_diam**3) * (pi/6.0_r8)  ! volume in cm^3
        tmp_mass = tmp_volu * 1.8_r8            ! mass in g
        cnum_h2so4 = (tmp_mass / 98.0_r8) * 6.023e23_r8   ! no. of h2so4 molec assuming pure h2so4
        cnum_tot = cnum_h2so4
        cnum_nh3 = 0.0_r8


        return
        end subroutine pbl_nuc_wang2008_native



!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine binary_nuc_vehk2002_select_impl()

        use cam_logfile, only: iulog
        use spmd_utils, only: masterproc

        implicit none

        character(len=48) :: impl_name
        integer :: status, n, i, code

        if (binary_nuc_vehk2002_impl_selected) return

        impl_name = 'codon'
        call get_environment_variable('BINARY_NUC_VEHK2002_IMPL', value=impl_name, length=n, status=status)

        if (status == 0 .and. n > 0) then
           do i = 1, n
              code = iachar(impl_name(i:i))
              if (code >= iachar('A') .and. code <= iachar('Z')) then
                 impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
              end if
           end do
           binary_nuc_vehk2002_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
        else
           binary_nuc_vehk2002_use_native_impl = .false.
        end if

        binary_nuc_vehk2002_impl_selected = .true.

        if (masterproc) then
           if (binary_nuc_vehk2002_use_native_impl) then
              write(iulog,*) 'binary_nuc_vehk2002 implementation = native'
           else
              write(iulog,*) 'binary_nuc_vehk2002 implementation = codon'
           end if
        end if

        end subroutine binary_nuc_vehk2002_select_impl

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine binary_nuc_vehk2002( temp, rh, so4vol,   &
            ratenucl, rateloge,   &
            cnum_h2so4, cnum_tot, radius_cluster )

        use iso_c_binding, only: c_double, c_loc, c_ptr

        implicit none

! subr arguments (in)
        real(r8), intent(in) :: temp
        real(r8), intent(in) :: rh
        real(r8), intent(in) :: so4vol

! subr arguments (out)
        real(r8), intent(out) :: ratenucl
        real(r8), intent(out) :: rateloge
        real(r8), intent(out) :: cnum_h2so4
        real(r8), intent(out) :: cnum_tot
        real(r8), intent(out) :: radius_cluster

        real(c_double), target :: ratenucl_work, rateloge_work
        real(c_double), target :: cnum_h2so4_work, cnum_tot_work, radius_cluster_work

        interface
           subroutine binary_nuc_vehk2002_codon(temp_c, rh_c, so4vol_c, ratenucl_p, rateloge_p, cnum_h2so4_p, cnum_tot_p, radius_cluster_p) &
                bind(c, name="binary_nuc_vehk2002_codon")
             use iso_c_binding, only: c_double, c_ptr
             real(c_double), value :: temp_c, rh_c, so4vol_c
             type(c_ptr), value :: ratenucl_p, rateloge_p, cnum_h2so4_p, cnum_tot_p, radius_cluster_p
           end subroutine binary_nuc_vehk2002_codon
        end interface

        call binary_nuc_vehk2002_select_impl()

        if (binary_nuc_vehk2002_use_native_impl) then
           call binary_nuc_vehk2002_native(temp, rh, so4vol, ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster)
           return
        end if

        ratenucl_work = 0.0_c_double
        rateloge_work = 0.0_c_double
        cnum_h2so4_work = 0.0_c_double
        cnum_tot_work = 0.0_c_double
        radius_cluster_work = 0.0_c_double

        call binary_nuc_vehk2002_codon( &
             real(temp, c_double), real(rh, c_double), real(so4vol, c_double), c_loc(ratenucl_work), c_loc(rateloge_work), &
             c_loc(cnum_h2so4_work), c_loc(cnum_tot_work), c_loc(radius_cluster_work) &
        )

        ratenucl = real(ratenucl_work, r8)
        rateloge = real(rateloge_work, r8)
        cnum_h2so4 = real(cnum_h2so4_work, r8)
        cnum_tot = real(cnum_tot_work, r8)
        radius_cluster = real(radius_cluster_work, r8)

        return

        end subroutine binary_nuc_vehk2002

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine binary_nuc_vehk2002_native( temp, rh, so4vol,   &
            ratenucl, rateloge,   &
            cnum_h2so4, cnum_tot, radius_cluster )
!
! calculates binary nucleation rate and critical cluster size
! using the parameterization in  
!     vehkamäki, h., m. kulmala, i. napari, k.e.j. lehtinen,
!        c. timmreck, m. noppel and a. laaksonen, 2002,
!        an improved parameterization for sulfuric acid-water nucleation
!        rates for tropospheric and stratospheric conditions,
!        j. geophys. res., 107, 4622, doi:10.1029/2002jd002184
!
        implicit none

! subr arguments (in)
        real(r8), intent(in) :: temp              ! temperature (k)  
        real(r8), intent(in) :: rh                ! relative humidity (0-1)
        real(r8), intent(in) :: so4vol            ! concentration of h2so4 (molecules cm-3)

! subr arguments (out)
        real(r8), intent(out) :: ratenucl         ! binary nucleation rate, j (# cm-3 s-1)
        real(r8), intent(out) :: rateloge         ! log( ratenucl )

        real(r8), intent(out) :: cnum_h2so4       ! number of h2so4 molecules
                                                  ! in the critical nucleus
        real(r8), intent(out) :: cnum_tot         ! total number of molecules
                                                  ! in the critical nucleus
        real(r8), intent(out) :: radius_cluster   ! the radius of cluster (nm)


! local variables
        real(r8) :: crit_x
        real(r8) :: acoe, bcoe, ccoe, dcoe, ecoe, fcoe, gcoe, hcoe, icoe, jcoe
        real(r8) :: tmpa, tmpb

! executable


! calc sulfuric acid mole fraction in critical cluster
        crit_x = 0.740997_r8 - 0.00266379_r8 * temp   &
               - 0.00349998_r8 * log (so4vol)   &
               + 0.0000504022_r8 * temp * log (so4vol)   &
               + 0.00201048_r8 * log (rh)   &
               - 0.000183289_r8 * temp * log (rh)   &
               + 0.00157407_r8 * (log (rh)) ** 2.0_r8   &
               - 0.0000179059_r8 * temp * (log (rh)) ** 2.0_r8   &
               + 0.000184403_r8 * (log (rh)) ** 3.0_r8   &
               - 1.50345e-6_r8 * temp * (log (rh)) ** 3.0_r8


! calc nucleation rate
        acoe    = 0.14309_r8+2.21956_r8*temp   &
                - 0.0273911_r8 * temp**2.0_r8   &
                + 0.0000722811_r8 * temp**3.0_r8 + 5.91822_r8/crit_x

        bcoe    = 0.117489_r8 + 0.462532_r8 *temp   &
                - 0.0118059_r8 * temp**2.0_r8   &
                + 0.0000404196_r8 * temp**3.0_r8 + 15.7963_r8/crit_x

        ccoe    = -0.215554_r8-0.0810269_r8 * temp   &
                + 0.00143581_r8 * temp**2.0_r8   &
                - 4.7758e-6_r8 * temp**3.0_r8   &
                - 2.91297_r8/crit_x

        dcoe    = -3.58856_r8+0.049508_r8 * temp   &
                - 0.00021382_r8 * temp**2.0_r8   &
                + 3.10801e-7_r8 * temp**3.0_r8   &
                - 0.0293333_r8/crit_x

        ecoe    = 1.14598_r8 - 0.600796_r8 * temp   &
                + 0.00864245_r8 * temp**2.0_r8   &
                - 0.0000228947_r8 * temp**3.0_r8   &
                - 8.44985_r8/crit_x

        fcoe    = 2.15855_r8 + 0.0808121_r8 * temp   &
                -0.000407382_r8 * temp**2.0_r8   &
                -4.01957e-7_r8 * temp**3.0_r8   &
                + 0.721326_r8/crit_x

        gcoe    = 1.6241_r8 - 0.0160106_r8 * temp   &
                + 0.0000377124_r8 * temp**2.0_r8   &
                + 3.21794e-8_r8 * temp**3.0_r8   &
                - 0.0113255_r8/crit_x

        hcoe    = 9.71682_r8 - 0.115048_r8 * temp   &
                + 0.000157098_r8 * temp**2.0_r8   &
                + 4.00914e-7_r8 * temp**3.0_r8   &
                + 0.71186_r8/crit_x

        icoe    = -1.05611_r8 + 0.00903378_r8 * temp   &
                - 0.0000198417_r8 * temp**2.0_r8   &
                + 2.46048e-8_r8  * temp**3.0_r8   &
                - 0.0579087_r8/crit_x

        jcoe    = -0.148712_r8 + 0.00283508_r8 * temp   &
                - 9.24619e-6_r8  * temp**2.0_r8   &
                + 5.00427e-9_r8 * temp**3.0_r8   &
                - 0.0127081_r8/crit_x

        tmpa     =     (   &
                  acoe   &
                + bcoe * log (rh)   &
                + ccoe * ( log (rh))**2.0_r8   &
                + dcoe * ( log (rh))**3.0_r8   &
                + ecoe * log (so4vol)   &
                + fcoe * (log (rh)) * (log (so4vol))   &
                + gcoe * ((log (rh) ) **2.0_r8)   &
                       * (log (so4vol))   &
                + hcoe * (log (so4vol)) **2.0_r8   &
                + icoe * log (rh)   &
                       * ((log (so4vol)) **2.0_r8)   &
                + jcoe * (log (so4vol)) **3.0_r8   &
                )
        rateloge = tmpa
        tmpa = min( tmpa, log(1.0e38_r8) )
        ratenucl = exp ( tmpa )
!       write(*,*) 'tmpa, ratenucl =', tmpa, ratenucl



! calc number of molecules in critical cluster
        acoe    = -0.00295413_r8 - 0.0976834_r8*temp   &
                + 0.00102485_r8 * temp**2.0_r8   &
                - 2.18646e-6_r8 * temp**3.0_r8 - 0.101717_r8/crit_x

        bcoe    = -0.00205064_r8 - 0.00758504_r8*temp   &
                + 0.000192654_r8 * temp**2.0_r8   &
                - 6.7043e-7_r8 * temp**3.0_r8 - 0.255774_r8/crit_x

        ccoe    = +0.00322308_r8 + 0.000852637_r8 * temp   &
                - 0.0000154757_r8 * temp**2.0_r8   &
                + 5.66661e-8_r8 * temp**3.0_r8   &
                + 0.0338444_r8/crit_x

        dcoe    = +0.0474323_r8 - 0.000625104_r8 * temp   &
                + 2.65066e-6_r8 * temp**2.0_r8   &
                - 3.67471e-9_r8 * temp**3.0_r8   &
                - 0.000267251_r8/crit_x

        ecoe    = -0.0125211_r8 + 0.00580655_r8 * temp   &
                - 0.000101674_r8 * temp**2.0_r8   &
                + 2.88195e-7_r8 * temp**3.0_r8   &
                + 0.0942243_r8/crit_x

        fcoe    = -0.038546_r8 - 0.000672316_r8 * temp   &
                + 2.60288e-6_r8 * temp**2.0_r8   &
                + 1.19416e-8_r8 * temp**3.0_r8   &
                - 0.00851515_r8/crit_x

        gcoe    = -0.0183749_r8 + 0.000172072_r8 * temp   &
                - 3.71766e-7_r8 * temp**2.0_r8   &
                - 5.14875e-10_r8 * temp**3.0_r8   &
                + 0.00026866_r8/crit_x

        hcoe    = -0.0619974_r8 + 0.000906958_r8 * temp   &
                - 9.11728e-7_r8 * temp**2.0_r8   &
                - 5.36796e-9_r8 * temp**3.0_r8   &
                - 0.00774234_r8/crit_x

        icoe    = +0.0121827_r8 - 0.00010665_r8 * temp   &
                + 2.5346e-7_r8 * temp**2.0_r8   &
                - 3.63519e-10_r8 * temp**3.0_r8   &
                + 0.000610065_r8/crit_x

        jcoe    = +0.000320184_r8 - 0.0000174762_r8 * temp   &
                + 6.06504e-8_r8 * temp**2.0_r8   &
                - 1.4177e-11_r8 * temp**3.0_r8   &
                + 0.000135751_r8/crit_x

        cnum_tot = exp (   &
                  acoe   &
                + bcoe * log (rh)   &
                + ccoe * ( log (rh))**2.0_r8   &
                + dcoe * ( log (rh))**3.0_r8   &
                + ecoe * log (so4vol)   &
                + fcoe * (log (rh)) * (log (so4vol))   &
                + gcoe * ((log (rh) ) **2.0_r8)   &
                       * (log (so4vol))   &
                + hcoe * (log (so4vol)) **2.0_r8   &
                + icoe * log (rh)   &
                       * ((log (so4vol)) **2.0_r8)   &
                + jcoe * (log (so4vol)) **3.0_r8   &
                )

        cnum_h2so4 = cnum_tot * crit_x

!   calc radius (nm) of critical cluster
        radius_cluster = exp( -1.6524245_r8 + 0.42316402_r8*crit_x   &
                              + 0.3346648_r8*log(cnum_tot) )
      

      return
      end subroutine binary_nuc_vehk2002_native



!----------------------------------------------------------------------
!----------------------------------------------------------------------
subroutine modal_aero_newnuc_init

!-----------------------------------------------------------------------
!
! Purpose:
!    set do_adjust and do_aitken flags
!    create history fields for column tendencies associated with
!       modal_aero_calcsize
!
! Author: R. Easter
!
!-----------------------------------------------------------------------

use modal_aero_data
use modal_aero_rename

use cam_abortutils,   only:  endrun
use cam_history,      only:  addfld, add_default, fieldname_len, phys_decomp
use constituents,     only:  pcnst, cnst_get_ind, cnst_name
use spmd_utils,       only:  masterproc
use phys_control,     only: phys_getopts


implicit none

!-----------------------------------------------------------------------
! arguments

!-----------------------------------------------------------------------
! local
   integer  :: l_h2so4, l_nh3
   integer  :: lnumait, lnh4ait, lso4ait
   integer  :: l
   integer  :: m, mait

   character(len=fieldname_len)   :: tmpname
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit

   logical                        :: dotend(pcnst)
   logical                        :: history_aerosol      ! Output the MAM aerosol tendencies

   !-----------------------------------------------------------------------     
   
        call phys_getopts( history_aerosol_out        = history_aerosol   )


!   set these indices
!   skip if no h2so4 species
!   skip if no aitken mode so4 or num species
	l_h2so4_sv = 0
	l_nh3_sv = 0
	lnumait_sv = 0
	lnh4ait_sv = 0
	lso4ait_sv = 0

	call cnst_get_ind( 'H2SO4', l_h2so4, .false. )
	call cnst_get_ind( 'NH3', l_nh3, .false. )

	mait = modeptr_aitken
	if (mait > 0) then
	    lnumait = numptr_amode(mait)
	    lso4ait = lptr_so4_a_amode(mait)
	    lnh4ait = lptr_nh4_a_amode(mait)
	end if
	if ((l_h2so4  <= 0) .or. (l_h2so4 > pcnst)) then
	    write(*,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- l_h2so4 <= 0'
	    return
	else if ((lso4ait <= 0) .or. (lso4ait > pcnst)) then
	    write(*,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- lso4ait <= 0'
	    return
	else if ((lnumait <= 0) .or. (lnumait > pcnst)) then
	    write(*,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- lnumait <= 0'
	    return
	else if ((mait <= 0) .or. (mait > ntot_amode)) then
	    write(*,'(/a/)')   &
		'*** modal_aero_newnuc bypass -- modeptr_aitken <= 0'
	    return
	end if

	l_h2so4_sv = l_h2so4
	l_nh3_sv   = l_nh3
	lnumait_sv = lnumait
	lnh4ait_sv = lnh4ait
	lso4ait_sv = lso4ait

!
!   create history file column-tendency fields
!
	dotend(:) = .false.
	dotend(lnumait) = .true.
	dotend(lso4ait) = .true.
	dotend(l_h2so4) = .true.
	if ((l_nh3   > 0) .and. (l_nh3   <= pcnst) .and. &
	    (lnh4ait > 0) .and. (lnh4ait <= pcnst)) then
	    dotend(lnh4ait) = .true.
	    dotend(l_nh3) = .true.
	end if

	do l = 1, pcnst
	    if ( .not. dotend(l) ) cycle
	    tmpname = cnst_name(l)
	    unit = 'kg/m2/s'
	    do m = 1, ntot_amode
	        if (l == numptr_amode(m)) unit = '#/m2/s'
	    end do
	    fieldname = trim(tmpname) // '_sfnnuc1'
	    long_name = trim(tmpname) // ' modal_aero new particle nucleation column tendency'
	    call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
            if ( history_aerosol ) then 
               call add_default( fieldname, 1, ' ' )
            endif
	    if ( masterproc ) write(*,'(3(a,2x))') &
		'modal_aero_newnuc_init addfld', fieldname, unit
	end do ! l = ...


      return
      end subroutine modal_aero_newnuc_init



subroutine ternary_nuc_merik2007( t, rh, c2, c3, j_log, ntot, nacid, namm, r )
!subroutine ternary_fit(          t, rh, c2, c3, j_log, ntot, nacid, namm, r )
! *************************** ternary_fit.f90 ********************************
! joonas merikanto, 2006
!
! fortran 90 subroutine that calculates the parameterized composition 
! and nucleation rate of critical clusters in h2o-h2so4-nh3 vapor
!
! warning: the fit should not be used outside its limits of validity
! (limits indicated below)
!
! in:
! t:     temperature (k), limits 235-295 k
! rh:    relative humidity as fraction (eg. 0.5=50%) limits 0.05-0.95
! c2:    sulfuric acid concentration (molecules/cm3) limits 5x10^4 - 10^9 molecules/cm3
! c3:    ammonia mixing ratio (ppt) limits 0.1 - 1000 ppt
!
! out:
! j_log: logarithm of nucleation rate (1/(s cm3))
! ntot:  total number of molecules in the critical cluster
! nacid: number of sulfuric acid molecules in the critical cluster
! namm:  number of ammonia molecules in the critical cluster
! r:     radius of the critical cluster (nm)
!  ****************************************************************************
implicit none

real(r8), intent(in) :: t, rh, c2, c3
real(r8), intent(out) :: j_log, ntot, nacid, namm, r
real(r8) :: j, t_onset

t_onset=143.6002929064716_r8 + 1.0178856665693992_r8*rh + &
   10.196398812974294_r8*log(c2) - &
   0.1849879416839113_r8*log(c2)**2 - 17.161783213150173_r8*log(c3) + &
   (109.92469248546053_r8*log(c3))/log(c2) + &
   0.7734119613144357_r8*log(c2)*log(c3) - 0.15576469879527022_r8*log(c3)**2

if(t_onset.gt.t) then 

   j_log=-12.861848898625231_r8 + 4.905527742256349_r8*c3 - 358.2337705052991_r8*rh -& 
   0.05463019231872484_r8*c3*t + 4.8630382337426985_r8*rh*t + &
   0.00020258394697064567_r8*c3*t**2 - 0.02175548069741675_r8*rh*t**2 - &
   2.502406532869512e-7_r8*c3*t**3 + 0.00003212869941055865_r8*rh*t**3 - &
   4.39129415725234e6_r8/log(c2)**2 + (56383.93843154586_r8*t)/log(c2)**2 -& 
   (239.835990963361_r8*t**2)/log(c2)**2 + &
   (0.33765136625580167_r8*t**3)/log(c2)**2 - &
   (629.7882041830943_r8*rh)/(c3**3*log(c2)) + &
   (7.772806552631709_r8*rh*t)/(c3**3*log(c2)) - &
   (0.031974053936299256_r8*rh*t**2)/(c3**3*log(c2)) + &
   (0.00004383764128775082_r8*rh*t**3)/(c3**3*log(c2)) + &
   1200.472096232311_r8*log(c2) - 17.37107890065621_r8*t*log(c2) + &
   0.08170681335921742_r8*t**2*log(c2) - 0.00012534476159729881_r8*t**3*log(c2) - &
   14.833042158178936_r8*log(c2)**2 + 0.2932631303555295_r8*t*log(c2)**2 - &
   0.0016497524241142845_r8*t**2*log(c2)**2 + &
   2.844074805239367e-6_r8*t**3*log(c2)**2 - 231375.56676032578_r8*log(c3) - &
   100.21645273730675_r8*rh*log(c3) + 2919.2852552424706_r8*t*log(c3) + &
   0.977886555834732_r8*rh*t*log(c3) - 12.286497122264588_r8*t**2*log(c3) - &
   0.0030511783284506377_r8*rh*t**2*log(c3) + &
   0.017249301826661612_r8*t**3*log(c3) + 2.967320346100855e-6_r8*rh*t**3*log(c3) + &
   (2.360931724951942e6_r8*log(c3))/log(c2) - &
   (29752.130254319443_r8*t*log(c3))/log(c2) + &
   (125.04965118142027_r8*t**2*log(c3))/log(c2) - &
   (0.1752996881934318_r8*t**3*log(c3))/log(c2) + &
   5599.912337254629_r8*log(c2)*log(c3) - 70.70896612937771_r8*t*log(c2)*log(c3) + &
   0.2978801613269466_r8*t**2*log(c2)*log(c3) - &
   0.00041866525019504_r8*t**3*log(c2)*log(c3) + 75061.15281456841_r8*log(c3)**2 - &
   931.8802278173565_r8*t*log(c3)**2 + 3.863266220840964_r8*t**2*log(c3)**2 - &
   0.005349472062284983_r8*t**3*log(c3)**2 - &
   (732006.8180571689_r8*log(c3)**2)/log(c2) + &
   (9100.06398573816_r8*t*log(c3)**2)/log(c2) - &
   (37.771091915932004_r8*t**2*log(c3)**2)/log(c2) + &
   (0.05235455395566905_r8*t**3*log(c3)**2)/log(c2) - &
   1911.0303773001353_r8*log(c2)*log(c3)**2 + &
   23.6903969622286_r8*t*log(c2)*log(c3)**2 - &
   0.09807872005428583_r8*t**2*log(c2)*log(c3)**2 + &
   0.00013564560238552576_r8*t**3*log(c2)*log(c3)**2 - &
   3180.5610833308_r8*log(c3)**3 + 39.08268568672095_r8*t*log(c3)**3 - &
   0.16048521066690752_r8*t**2*log(c3)**3 + &
   0.00022031380023793877_r8*t**3*log(c3)**3 + &
   (40751.075322248245_r8*log(c3)**3)/log(c2) - &
   (501.66977622013934_r8*t*log(c3)**3)/log(c2) + &
   (2.063469732254135_r8*t**2*log(c3)**3)/log(c2) - &
   (0.002836873785758324_r8*t**3*log(c3)**3)/log(c2) + &
   2.792313345723013_r8*log(c2)**2*log(c3)**3 - &
   0.03422552111802899_r8*t*log(c2)**2*log(c3)**3 + &
   0.00014019195277521142_r8*t**2*log(c2)**2*log(c3)**3 - &
   1.9201227328396297e-7_r8*t**3*log(c2)**2*log(c3)**3 - &
   980.923146020468_r8*log(rh) + 10.054155220444462_r8*t*log(rh) - &
   0.03306644502023841_r8*t**2*log(rh) + 0.000034274041225891804_r8*t**3*log(rh) + &
   (16597.75554295064_r8*log(rh))/log(c2) - &
   (175.2365504237746_r8*t*log(rh))/log(c2) + &
   (0.6033215603167458_r8*t**2*log(rh))/log(c2) - &
   (0.0006731787599587544_r8*t**3*log(rh))/log(c2) - &
   89.38961120336789_r8*log(c3)*log(rh) + 1.153344219304926_r8*t*log(c3)*log(rh) - &
   0.004954549700267233_r8*t**2*log(c3)*log(rh) + &
   7.096309866238719e-6_r8*t**3*log(c3)*log(rh) + &
   3.1712136610383244_r8*log(c3)**3*log(rh) - &
   0.037822330602328806_r8*t*log(c3)**3*log(rh) + &
   0.0001500555743561457_r8*t**2*log(c3)**3*log(rh) - &
   1.9828365865570703e-7_r8*t**3*log(c3)**3*log(rh)

   j=exp(j_log)

   ntot=57.40091052369212_r8 - 0.2996341884645408_r8*t + &
   0.0007395477768531926_r8*t**2 - &
   5.090604835032423_r8*log(c2) + 0.011016634044531128_r8*t*log(c2) + &
   0.06750032251225707_r8*log(c2)**2 - 0.8102831333223962_r8*log(c3) + &
   0.015905081275952426_r8*t*log(c3) - 0.2044174683159531_r8*log(c2)*log(c3) + &
   0.08918159167625832_r8*log(c3)**2 - 0.0004969033586666147_r8*t*log(c3)**2 + &
   0.005704394549007816_r8*log(c3)**3 + 3.4098703903474368_r8*log(j) - &
   0.014916956508210809_r8*t*log(j) + 0.08459090011666293_r8*log(c3)*log(j) - &
   0.00014800625143907616_r8*t*log(c3)*log(j) + 0.00503804694656905_r8*log(j)**2
 
   r=3.2888553966535506e-10_r8 - 3.374171768439839e-12_r8*t + &
   1.8347359507774313e-14_r8*t**2 + 2.5419844298881856e-12_r8*log(c2) - &
   9.498107643050827e-14_r8*t*log(c2) + 7.446266520834559e-13_r8*log(c2)**2 + &
   2.4303397746137294e-11_r8*log(c3) + 1.589324325956633e-14_r8*t*log(c3) - &
   2.034596219775266e-12_r8*log(c2)*log(c3) - 5.59303954457172e-13_r8*log(c3)**2 - &
   4.889507104645867e-16_r8*t*log(c3)**2 + 1.3847024107506764e-13_r8*log(c3)**3 + &
   4.141077193427042e-15_r8*log(j) - 2.6813110884009767e-14_r8*t*log(j) + &
   1.2879071621313094e-12_r8*log(c3)*log(j) - &
   3.80352446061867e-15_r8*t*log(c3)*log(j) - 1.8790172502456827e-14_r8*log(j)**2
 
   nacid=-4.7154180661803595_r8 + 0.13436423483953885_r8*t - & 
   0.00047184686478816176_r8*t**2 - & 
   2.564010713640308_r8*log(c2) + 0.011353312899114723_r8*t*log(c2) + &
   0.0010801941974317014_r8*log(c2)**2 + 0.5171368624197119_r8*log(c3) - &
   0.0027882479896204665_r8*t*log(c3) + 0.8066971907026886_r8*log(c3)**2 - & 
   0.0031849094214409335_r8*t*log(c3)**2 - 0.09951184152927882_r8*log(c3)**3 + &
   0.00040072788891745513_r8*t*log(c3)**3 + 1.3276469271073974_r8*log(j) - &
   0.006167654171986281_r8*t*log(j) - 0.11061390967822708_r8*log(c3)*log(j) + &
   0.0004367575329273496_r8*t*log(c3)*log(j) + 0.000916366357266258_r8*log(j)**2
 
   namm=71.20073903979772_r8 - 0.8409600103431923_r8*t + &
   0.0024803006590334922_r8*t**2 + &
   2.7798606841602607_r8*log(c2) - 0.01475023348171676_r8*t*log(c2) + &
   0.012264508212031405_r8*log(c2)**2 - 2.009926050440182_r8*log(c3) + &
   0.008689123511431527_r8*t*log(c3) - 0.009141180198955415_r8*log(c2)*log(c3) + &
   0.1374122553905617_r8*log(c3)**2 - 0.0006253227821679215_r8*t*log(c3)**2 + &
   0.00009377332742098946_r8*log(c3)**3 + 0.5202974341687757_r8*log(j) - &
   0.002419872323052805_r8*t*log(j) + 0.07916392322884074_r8*log(c3)*log(j) - &
   0.0003021586030317366_r8*t*log(c3)*log(j) + 0.0046977006608603395_r8*log(j)**2

else
! nucleation rate less that 5e-6, setting j_log arbitrary small
   j_log=-300._r8
end if

return

end  subroutine ternary_nuc_merik2007

!----------------------------------------------------------------------
end module modal_aero_newnuc
