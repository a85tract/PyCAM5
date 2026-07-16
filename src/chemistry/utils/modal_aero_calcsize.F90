module modal_aero_calcsize

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: masterproc
use physconst,        only: pi, rhoh2o, gravit

use ppgrid,           only: pcols, pver
use physics_types,    only: physics_state, physics_ptend
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field

use phys_control,     only: phys_getopts
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_props, rad_cnst_get_mode_num

use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun
use cam_history,      only: addfld, add_default, fieldname_len, phys_decomp, outfld
use constituents,     only: pcnst, cnst_name

use ref_pres,         only: top_lev => clim_modal_aero_top_lev

#ifdef MODAL_AERO

! these are the variables needed for the diagnostic calculation of dry radius
use modal_aero_data, only: ntot_amode, nspec_amode, &
                           numptr_amode, &
                           alnsg_amode, &
                           voltonumbhi_amode, voltonumblo_amode, &
                           dgnum_amode, dgnumhi_amode, dgnumlo_amode


! these variables are needed for the prognostic calculations to exchange mass
! between modes
use modal_aero_data,  only: numptrcw_amode, mprognum_amode, qqcw_get_field, lmassptrcw_amode, &
           lmassptr_amode, modeptr_accum, modeptr_aitken, ntot_aspectype, &
           lspectype_amode, specmw_amode, specdens_amode, voltonumb_amode, &
           cnst_name_cw

use modal_aero_rename, only: lspectooa_renamexf, lspecfrma_renamexf, lspectooc_renamexf, lspecfrmc_renamexf, &
           modetoo_renamexf, nspecfrm_renamexf, npair_renamexf, modefrm_renamexf


#endif


implicit none
private
save

public modal_aero_calcsize_init, modal_aero_calcsize_sub, modal_aero_calcsize_diag
public :: modal_aero_calcsize_reg

logical :: do_adjust_default
logical :: do_aitacc_transfer_default

integer :: dgnum_idx = -1

!===============================================================================
contains
!===============================================================================

subroutine modal_aero_calcsize_reg()
  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use rad_constituents, only: rad_cnst_get_info

  integer :: nmodes
  
  call rad_cnst_get_info(0, nmodes=nmodes)

  call pbuf_add_field('DGNUM', 'global',  dtype_r8, (/pcols, pver, nmodes/), dgnum_idx)    

end subroutine modal_aero_calcsize_reg

!===============================================================================
!===============================================================================

subroutine modal_aero_calcsize_init(pbuf2d)
   use time_manager,  only: is_first_step
   use physics_buffer,only: pbuf_set_field

   !-----------------------------------------------------------------------
   !
   ! Purpose:
   !    set do_adjust_default and do_aitacc_transfer_default flags
   !    create history fields for column tendencies associated with
   !       modal_aero_calcsize
   !
   ! Author: R. Easter
   !
   !-----------------------------------------------------------------------

   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local
   integer  :: ipair, iq
   integer  :: jac
   integer  :: lsfrm, lstoo
   integer  :: n, nacc, nait
   logical  :: history_aerosol

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit
   !-----------------------------------------------------------------------

   call phys_getopts(history_aerosol_out=history_aerosol)

   ! init entities required for both prescribed and prognostic modes

   if (is_first_step()) then
      ! initialize fields in physics buffer
      call pbuf_set_field(pbuf2d, dgnum_idx, 0.0_r8)
   endif

#ifndef MODAL_AERO
   do_adjust_default          = .false.
   do_aitacc_transfer_default = .false.
#else
   !  do_adjust_default allows adjustment to be turned on/off
   do_adjust_default = .true.

   !  do_aitacc_transfer_default allows aitken <--> accum mode transfer to be turned on/off
   !  *** it can only be true when aitken & accum modes are both present
   !      and have prognosed number and diagnosed surface/sigmag
   nait = modeptr_aitken
   nacc = modeptr_accum
   do_aitacc_transfer_default = .false.
   if ((modeptr_aitken > 0) .and.   &
      (modeptr_accum  > 0) .and.   &
      (modeptr_aitken /= modeptr_accum)) then
      do_aitacc_transfer_default = .true.
      if (mprognum_amode(nait) <= 0) do_aitacc_transfer_default = .false.
      if (mprognum_amode(nacc) <= 0) do_aitacc_transfer_default = .false.
   end if

   if ( .not. do_adjust_default ) return

   !  define history fields for number-adjust source-sink for all modes
   do n = 1, ntot_amode 
      if (mprognum_amode(n) <= 0) cycle

      do jac = 1, 2
         if (jac == 1) then
            tmpnamea = cnst_name(numptr_amode(n))
         else
            tmpnamea = cnst_name_cw(numptrcw_amode(n))
         end if
         unit = '#/m2/s'
         fieldname = trim(tmpnamea) // '_sfcsiz1'
         long_name = trim(tmpnamea) // ' calcsize number-adjust column source'
         call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

         fieldname = trim(tmpnamea) // '_sfcsiz2'
         long_name = trim(tmpnamea) // ' calcsize number-adjust column sink'
         call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname
      end do   ! jac = ...
   end do   ! n = ...

   if ( .not. do_aitacc_transfer_default ) return

   ! check that renaming ipair=1 is aitken-->accum
   ipair = 1
   if ((modefrm_renamexf(ipair) .ne. nait) .or.   &
      (modetoo_renamexf(ipair) .ne. nacc)) then
      write( 6, '(//2a//)' )   &
         '*** modal_aero_calcaersize_init error -- ',   &
         'modefrm/too_renamexf(1) are wrong'
      call endrun( 'modal_aero_calcaersize_init error' )
   end if

   ! define history fields for aitken-accum transfer
   do iq = 1, nspecfrm_renamexf(ipair)

      ! jac=1 does interstitial ("_a"); jac=2 does activated ("_c"); 
      do jac = 1, 2

         ! the lspecfrma_renamexf (and lspecfrmc_renamexf) are aitken species
         ! the lspectooa_renamexf (and lspectooc_renamexf) are accum  species
         if (jac .eq. 1) then
            lsfrm = lspecfrma_renamexf(iq,ipair)
            lstoo = lspectooa_renamexf(iq,ipair)
         else
            lsfrm = lspecfrmc_renamexf(iq,ipair)
            lstoo = lspectooc_renamexf(iq,ipair)
         end if
         if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

         if (jac .eq. 1) then
            tmpnamea = cnst_name(lsfrm)
            tmpnameb = cnst_name(lstoo)
         else
            tmpnamea = cnst_name_cw(lsfrm)
            tmpnameb = cnst_name_cw(lstoo)
         end if

         unit = 'kg/m2/s'
         if ((tmpnamea(1:3) == 'num') .or. &
            (tmpnamea(1:3) == 'NUM')) unit = '#/m2/s'
         fieldname = trim(tmpnamea) // '_sfcsiz3'
         long_name = trim(tmpnamea) // ' calcsize aitken-to-accum adjust column tendency'
         call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

         fieldname = trim(tmpnameb) // '_sfcsiz3'
         long_name = trim(tmpnameb) // ' calcsize aitken-to-accum adjust column tendency'
         call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

         fieldname = trim(tmpnamea) // '_sfcsiz4'
         long_name = trim(tmpnamea) // ' calcsize accum-to-aitken adjust column tendency'
         call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

         fieldname = trim(tmpnameb) // '_sfcsiz4'
         long_name = trim(tmpnameb) // ' calcsize accum-to-aitken adjust column tendency'
         call addfld( fieldname, unit, 1, 'A', long_name, phys_decomp )
         if (history_aerosol) then
            call add_default(fieldname, 1, ' ')
         end if
         if ( masterproc ) write(*,'(2a)') 'calcsize addfld - ', fieldname

      end do   ! jac = ...
   end do   ! iq = ...

#endif

end subroutine modal_aero_calcsize_init

!===============================================================================

subroutine modal_aero_calcsize_sub(state, ptend, deltat, pbuf, do_adjust_in, &
   do_aitacc_transfer_in)

   use ap_modal_aero_calcsize_sub_scheme, only : modal_aero_calcsize_sub_run
   use perf_mod, only : t_startf, t_stopf

   !-----------------------------------------------------------------------
   !
   ! Calculates aerosol size distribution parameters 
   !    mprognum_amode >  0
   !       calculate Dgnum from mass, number, and fixed sigmag
   !    mprognum_amode <= 0
   !       calculate number from mass, fixed Dgnum, and fixed sigmag
   !
   ! Also (optionally) adjusts prognostic number to
   !    be within bounds determined by mass, Dgnum bounds, and sigma bounds
   !
   ! Author: R. Easter
   !
   !-----------------------------------------------------------------------

   ! arguments
   type(physics_state), target, intent(in)    :: state       ! Physics state variables
   type(physics_ptend), target, intent(inout) :: ptend       ! indivdual parameterization tendencies
   real(r8),                    intent(in)    :: deltat      ! model time-step size (s)
   type(physics_buffer_desc),   pointer       :: pbuf(:)     ! physics buffer

   logical, optional :: do_adjust_in
   logical, optional :: do_aitacc_transfer_in


   call t_startf('ap_modal_aero_calcsize_sub_run')
   call modal_aero_calcsize_sub_run(state, ptend, deltat, pbuf, &
        do_adjust_default, do_aitacc_transfer_default, dgnum_idx, &
        do_adjust_in, do_aitacc_transfer_in)
   call t_stopf('ap_modal_aero_calcsize_sub_run')

end subroutine modal_aero_calcsize_sub
 

!----------------------------------------------------------------------


subroutine modal_aero_calcsize_diag(state, pbuf, list_idx_in, dgnum_m)

   !-----------------------------------------------------------------------
   !
   ! Calculate aerosol size distribution parameters 
   !
   ! ***N.B.*** DGNUM for the modes in the climate list are put directly into
   !            the physics buffer.  For diagnostic list calculations use the
   !            optional list_idx and dgnum args.
   !-----------------------------------------------------------------------

   ! arguments
   type(physics_state), intent(in), target :: state   ! Physics state variables
   type(physics_buffer_desc), pointer :: pbuf(:)      ! physics buffer

   integer,  optional, intent(in)   :: list_idx_in    ! diagnostic list index
   real(r8), optional, pointer      :: dgnum_m(:,:,:) ! interstital aerosol dry number mode radius (m)

   ! local
   integer  :: i, k, l1, n
   integer  :: lchnk, ncol
   integer  :: list_idx, stat
   integer  :: nmodes
   integer  :: nspec

   real(r8), pointer :: dgncur_a(:,:) ! (pcols,pver)


   real(r8), parameter :: third = 1.0_r8/3.0_r8

   real(r8), pointer :: mode_num(:,:) ! mode number mixing ratio
   real(r8), pointer :: specmmr(:,:)  ! specie mmr
   real(r8)          :: specdens      ! specie density

   real(r8) :: dryvol_a(pcols,pver)   ! interstital aerosol dry volume (cm^3/mol_air)

   real(r8) :: dgnum, dgnumhi, dgnumlo
   real(r8) :: dgnyy, dgnxx           ! dgnumlo/hi of current mode
   real(r8) :: drv_a                  ! dry volume (cm3/mol_air)
   real(r8) :: dumfac, dummwdens      ! work variables
   real(r8) :: num_a0                 ! initial number (#/mol_air)
   real(r8) :: num_a                  ! final number (#/mol_air)
   real(r8) :: voltonumbhi, voltonumblo
   real(r8) :: v2nyy, v2nxx           ! voltonumblo/hi of current mode
   real(r8) :: sigmag, alnsg
   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   list_idx = 0  ! climate list by default
   if (present(list_idx_in)) list_idx = list_idx_in

   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx /= 0) then
      if (.not. present(dgnum_m)) then
         call endrun('modal_aero_calcsize_diag called for'// &
                     'diagnostic list but dgnum_m pointer not present')
      end if
      allocate(dgnum_m(pcols,pver,nmodes), stat=stat)
      if (stat > 0) then
         call endrun('modal_aero_calcsize_diag: allocation FAILURE: dgnum_m')
      end if
   end if

   do n = 1, nmodes

      if (list_idx == 0) then
         call pbuf_get_field(pbuf, dgnum_idx, dgncur_a, start=(/1,1,n/), kount=(/pcols,pver,1/))
      else
         dgncur_a => dgnum_m(:,:,n)
      end if

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, n, dgnum=dgnum, dgnumhi=dgnumhi, dgnumlo=dgnumlo, &
                                   sigmag=sigmag)

      ! get mode number mixing ratio
      call rad_cnst_get_mode_num(list_idx, n, 'a', state, pbuf, mode_num)

      dgncur_a(:,:) = dgnum
      dryvol_a(:,:) = 0.0_r8

      ! compute dry volume mixrats = 
      !      sum_over_components{ component_mass mixrat / density }
      call rad_cnst_get_info(list_idx, n, nspec=nspec)
      do l1 = 1, nspec

         call rad_cnst_get_aer_mmr(list_idx, n, l1, 'a', state, pbuf, specmmr)
         call rad_cnst_get_aer_props(list_idx, n, l1, density_aer=specdens)

         ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
         dummwdens = 1.0_r8 / specdens

         do k=top_lev,pver
            do i=1,ncol
               dryvol_a(i,k) = dryvol_a(i,k)    &
                  + max(0.0_r8, specmmr(i,k))*dummwdens
            end do
         end do
      end do

      alnsg  = log( sigmag )
      dumfac = exp(4.5_r8*alnsg**2)*pi/6.0_r8
      voltonumblo = 1._r8 / ( (pi/6._r8)*(dgnumlo**3)*exp(4.5_r8*alnsg**2) )
      voltonumbhi = 1._r8 / ( (pi/6._r8)*(dgnumhi**3)*exp(4.5_r8*alnsg**2) )
      v2nxx = voltonumbhi
      v2nyy = voltonumblo
      dgnxx = dgnumhi
      dgnyy = dgnumlo

      do k = top_lev, pver
         do i = 1, ncol

            drv_a = dryvol_a(i,k)
            num_a0 = mode_num(i,k)
            num_a = max( 0.0_r8, num_a0 )

            if (drv_a > 0.0_r8) then
               if (num_a <= drv_a*v2nxx) then
                  dgncur_a(i,k) = dgnxx
               else if (num_a >= drv_a*v2nyy) then
                  dgncur_a(i,k) = dgnyy
               else
                  dgncur_a(i,k) = (drv_a/(dumfac*num_a))**third
               end if
            end if

         end do
      end do

   end do ! nmodes

end subroutine modal_aero_calcsize_diag

!----------------------------------------------------------------------

end module modal_aero_calcsize
