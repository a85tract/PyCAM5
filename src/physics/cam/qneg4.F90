
subroutine qneg4 (subnam  ,lchnk   ,ncol    ,ztodt   ,        &
                  qbot    ,srfrpdel,shflx   ,lhflx   ,qflx    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check if moisture flux into the ground is exceeding the total
! moisture content of the lowest model layer (creating negative moisture
! values).  If so, then subtract the excess from the moisture and
! latent heat fluxes and add it to the sensible heat flux.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Olson
! 
! Water isotopes added by J. Nusbaumer - Mar 2011
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,    only: get_lat_p, get_lon_p
   use physconst,    only: gravit, latvap
   use constituents, only: qmin, pcnst
   use cam_logfile,  only: iulog
   use spmd_utils,   only: masterproc

   !water isotopes:
   use water_types,       only: iwtvap
   use water_tracer_vars, only: trace_water, wisotope, WTRC_MAX_CNST, wtrc_iatype, wtrc_ntype, iwspec, wtrc_qmin, wtrc_fixed_rstd
   use water_tracers,     only: wtrc_get_rstd
   use water_isotopes,    only: pwtspec
   use iso_c_binding,     only: c_double, c_int64_t, c_loc, c_ptr

   implicit none

!
! Input arguments
!
   character*8, intent(in) :: subnam         ! name of calling routine
!
   integer, intent(in) :: lchnk              ! chunk index
   integer, intent(in) :: ncol               ! number of atmospheric columns
!
   real(r8), intent(in) :: ztodt             ! two times model timestep (2 delta-t)
   real(r8), target, intent(in) :: qbot(pcols,pcnst) ! moisture at lowest model level
   real(r8), target, intent(in) :: srfrpdel(pcols)   ! 1./(pint(K+1)-pint(K))
!
! Input/Output arguments
!
   real(r8), target, intent(inout) :: shflx(pcols)   ! Surface sensible heat flux (J/m2/s)
   real(r8), target, intent(inout) :: lhflx(pcols)   ! Surface latent   heat flux (J/m2/s)
   real(r8), target, intent(inout) :: qflx (pcols,pcnst)   ! surface water flux (kg/m^2/s)

   logical, save :: use_native_impl = .false.
   logical, save :: impl_selected = .false.

   integer :: ivap
   integer :: m
   integer(c_int64_t), target :: indxexc(pcols)    ! index array of points with excess flux
   integer(c_int64_t), target :: nptsexc           ! number of points with excess flux
   integer(c_int64_t), target :: iw                ! i index of worst violator
   integer(c_int64_t), target :: wtrc_iatype_vap(WTRC_MAX_CNST)
   integer(c_int64_t), target :: iwspec64(pcnst)
   real(c_double), target :: worst                 ! biggest violator
   real(c_double), target :: excess(pcols)         ! Excess downward sfc latent heat flux
   real(c_double), target :: qfxo(pcols,pcnst)     ! initial tracer flux
   real(c_double), target :: rstd(pwtspec)

   interface
      subroutine qneg4_codon(ncol_c, pcols_c, pcnst_c, trace_water_on_c, iwtvap_n_c, &
           qmin1_c, ztodt_c, gravit_c, latvap_c, wtrc_qmin_c, &
           qbot_p, srfrpdel_p, shflx_p, lhflx_p, qflx_p, indxexc_p, excess_p, qfxo_p, &
           nptsexc_p, worst_p, iw_p, wtrc_iatype_vap_p, iwspec_p, rstd_p) bind(c, name="qneg4_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pcnst_c, trace_water_on_c, iwtvap_n_c
         real(c_double), value :: qmin1_c, ztodt_c, gravit_c, latvap_c, wtrc_qmin_c
         type(c_ptr), value :: qbot_p, srfrpdel_p, shflx_p, lhflx_p, qflx_p, indxexc_p, excess_p, qfxo_p
         type(c_ptr), value :: nptsexc_p, worst_p, iw_p, wtrc_iatype_vap_p, iwspec_p, rstd_p
      end subroutine qneg4_codon
   end interface

!
!-----------------------------------------------------------------------
!

   call qneg4_batch_select_impl()

   if (use_native_impl) then
      call qneg4_native(subnam, lchnk, ncol, ztodt, qbot, srfrpdel, shflx, lhflx, qflx)
      return
   end if

   if (trace_water) then
      do ivap = 1, WTRC_MAX_CNST
         wtrc_iatype_vap(ivap) = 0_c_int64_t
      end do
      do ivap = 1, wtrc_ntype(iwtvap)
         wtrc_iatype_vap(ivap) = int(wtrc_iatype(ivap, iwtvap), c_int64_t)
      end do
      do m = 1, pcnst
         iwspec64(m) = int(iwspec(m), c_int64_t)
      end do
      do m = 1, pwtspec
         rstd(m) = real(wtrc_get_rstd(m), c_double)
      end do
   else
      wtrc_iatype_vap = 0_c_int64_t
      iwspec64 = 0_c_int64_t
      rstd = 0._c_double
   end if

   call qneg4_batch_log_entered()
   call qneg4_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), int(pcnst, c_int64_t), &
        merge(1_c_int64_t, 0_c_int64_t, trace_water), int(wtrc_ntype(iwtvap), c_int64_t), &
        real(qmin(1), c_double), real(ztodt, c_double), real(gravit, c_double), real(latvap, c_double), real(wtrc_qmin, c_double), &
        c_loc(qbot), c_loc(srfrpdel), c_loc(shflx), c_loc(lhflx), c_loc(qflx), c_loc(indxexc), c_loc(excess), c_loc(qfxo), &
        c_loc(nptsexc), c_loc(worst), c_loc(iw), c_loc(wtrc_iatype_vap), c_loc(iwspec64), c_loc(rstd) &
   )

   if (nptsexc > 10_c_int64_t) then
      write(iulog,9000) subnam, int(nptsexc), real(worst, r8), lchnk, int(iw), get_lat_p(lchnk, int(iw)), get_lon_p(lchnk, int(iw))
   end if

   return
9000 format(' QNEG4 WARNING from ',a8 &
            ,' Max possible LH flx exceeded at ',i6,' points. ' &
            ,', Worst excess = ',1pe12.4 &
            ,', lchnk = ',i6 &
            ,', i = ',i6 &
            ,', same as indices lat =', i6 &
            ,', lon =', i6 &
           )

contains

   subroutine qneg4_batch_append_proof(proof_line)
      character(len=*), intent(in) :: proof_line
      character(len=512) :: proof_file
      integer :: status, n, unitno

      proof_file = ''
      call get_environment_variable('QNEG_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
      if (status == 0 .and. n > 0) then
         open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
         write(unitno,'(A)') trim(proof_line)
         close(unitno)
      end if
   end subroutine qneg4_batch_append_proof

   subroutine qneg4_batch_select_impl()
      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('QNEG_BATCH_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_impl = .false.
      end if

      impl_selected = .true.

      if (masterproc) then
         if (use_native_impl) then
            write(iulog,*) 'qneg_batch qneg4 implementation = native'
            call qneg4_batch_append_proof('qneg_batch selector entered implementation = native (qneg4)')
         else
            write(iulog,*) 'qneg_batch qneg4 implementation = codon'
            call qneg4_batch_append_proof('qneg_batch selector entered implementation = codon (qneg4)')
         end if
         call flush(iulog)
      end if
   end subroutine qneg4_batch_select_impl

   subroutine qneg4_batch_log_entered()
      logical, save :: entered_logged = .false.

      if (entered_logged) return
      entered_logged = .true.

      if (masterproc) then
         write(iulog,'(A)') 'qneg_batch entered (qneg4 unified stage dispatch = codon)'
         call qneg4_batch_append_proof('qneg_batch entered (qneg4 unified stage dispatch = codon)')
         call flush(iulog)
      end if
   end subroutine qneg4_batch_log_entered

end subroutine qneg4

subroutine qneg4_native (subnam  ,lchnk   ,ncol    ,ztodt   ,        &
                         qbot    ,srfrpdel,shflx   ,lhflx   ,qflx    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check if moisture flux into the ground is exceeding the total
! moisture content of the lowest model layer (creating negative moisture
! values).  If so, then subtract the excess from the moisture and
! latent heat fluxes and add it to the sensible heat flux.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Olson
! 
! Water isotopes added by J. Nusbaumer - Mar 2011
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,    only: get_lat_p, get_lon_p
   use physconst,    only: gravit, latvap
   use constituents, only: qmin, pcnst
   use cam_logfile,  only: iulog

   !water isotopes:
   use water_types,   only: iwtvap
   use water_tracer_vars, only: trace_water, wtrc_iatype, wtrc_ntype, iwspec
   use water_tracers, only: wtrc_ratio

   implicit none

!
! Input arguments
!
   character*8, intent(in) :: subnam         ! name of calling routine
!
   integer, intent(in) :: lchnk              ! chunk index
   integer, intent(in) :: ncol               ! number of atmospheric columns
!
   real(r8), intent(in) :: ztodt             ! two times model timestep (2 delta-t)
   real(r8), intent(in) :: qbot(pcols,pcnst) ! moisture at lowest model level
   real(r8), intent(in) :: srfrpdel(pcols)   ! 1./(pint(K+1)-pint(K))
!
! Input/Output arguments
!
   real(r8), intent(inout) :: shflx(pcols)   ! Surface sensible heat flux (J/m2/s)
   real(r8), intent(inout) :: lhflx(pcols)   ! Surface latent   heat flux (J/m2/s)
   real(r8), intent(inout) :: qflx (pcols,pcnst)   ! surface water flux (kg/m^2/s)
!
!---------------------------Local workspace-----------------------------
!
   integer :: i,ii              ! longitude indices
   integer :: iw                ! i index of worst violator
   integer :: indxexc(pcols)    ! index array of points with excess flux
   integer :: nptsexc           ! number of points with excess flux
   integer :: m                 ! loop control variable for water isotopes
   integer :: ivap              ! isotope index
!
   real(r8):: worst             ! biggest violator
   real(r8):: excess(pcols)     ! Excess downward sfc latent heat flux

!water isotopes:
   real(r8):: qfxo(pcols,pcnst)   ! initial tracer flux
   real(r8):: rat                 ! tracer ratio

!
!-----------------------------------------------------------------------
!

! Store old value to input for water tracers

   do m = 1, pcnst
     do i = 1, ncol
       qfxo(i,m) = qflx(i,m)
     end do
   end do

! Compute excess downward (negative) q flux compared to a theoretical
! maximum downward q flux.  The theoretical max is based upon the
! given moisture content of lowest level of the model atmosphere.
!
   nptsexc = 0
   do i = 1,ncol
      excess(i) = qflx(i,1) - (qmin(1) - qbot(i,1))/(ztodt*gravit*srfrpdel(i))
!
! If there is an excess downward (negative) q flux, then subtract
! excess from "qflx" and "lhflx" and add to "shflx".
!
      if (excess(i) < 0._r8) then
         nptsexc = nptsexc + 1
         indxexc(nptsexc) = i
         qflx (i,1) = qflx (i,1) - excess(i)
         lhflx(i) = lhflx(i) - excess(i)*latvap
         shflx(i) = shflx(i) + excess(i)*latvap
      end if
   end do
!
! Write out worst value if excess
!
   if (nptsexc.gt.10) then
      worst = 0._r8
      do ii=1,nptsexc
         i = indxexc(ii)
         if (excess(i) < worst) then
            worst = excess(i)
            iw = i
         end if
      end do
      write(iulog,9000) subnam,nptsexc,worst, lchnk, iw, get_lat_p(lchnk,iw),get_lon_p(lchnk,iw)
   end if
!
! Water tracers: where total has change, modify tracers to conserve ratios
!

   if (trace_water) then
     !NOTE:  qfxo may not be needed, as ratio is against H2O tracer, not q. - JN
     do ivap = 1, wtrc_ntype(iwtvap)
       m = wtrc_iatype(ivap, iwtvap)
      
       do ii = 1, nptsexc
         i = indxexc(ii)
!         rat = wtrc_ratio(iwspec(m), qfxo(i,m),qfxo(i,1))
         rat = wtrc_ratio(iwspec(m),qfxo(i,m),qfxo(i,wtrc_iatype(1,iwtvap)))
         qflx(i,m) = qflx(i,m) - rat*excess(i)
       end do
     end do
   end if
!

   return
9000 format(' QNEG4 WARNING from ',a8 &
            ,' Max possible LH flx exceeded at ',i6,' points. ' &
            ,', Worst excess = ',1pe12.4 &
            ,', lchnk = ',i6 &
            ,', i = ',i6 &
            ,', same as indices lat =', i6 &
            ,', lon =', i6 &
           )
end subroutine qneg4_native
