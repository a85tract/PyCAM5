module modal_aero_calcsize

!   RCE 07.04.13:  Adapted from MIRAGE2 code

use shr_kind_mod,     only: r8 => shr_kind_r8
use iso_c_binding,    only: c_int64_t
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
use modal_aero_data, only: ntot_amode, nspec_amode, maxd_aspectype, &
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
           modetoo_renamexf, nspecfrm_renamexf, npair_renamexf, modefrm_renamexf, &
           maxpair_renamexf, maxspec_renamexf


#endif


implicit none
private
save

public modal_aero_calcsize_init, modal_aero_calcsize_sub, modal_aero_calcsize_diag
public :: modal_aero_calcsize_reg

logical :: do_adjust_default
logical :: do_aitacc_transfer_default
logical :: modal_aero_calcsize_sub_use_native_impl = .false.
logical :: modal_aero_calcsize_sub_impl_selected = .false.
logical :: modal_aero_calcsize_sub_proof_written = .false.
logical :: modal_aero_calcsize_sub_wrap_proof_written = .false.
logical :: modal_aero_calcsize_reg_logged = .false.

integer :: dgnum_idx = -1

interface
   function modal_aero_calcsize_reg_codon(stage_c) result(out_c) bind(c, name="modal_aero_calcsize_reg_codon")
      import :: c_int64_t
      integer(c_int64_t), value :: stage_c
      integer(c_int64_t) :: out_c
   end function modal_aero_calcsize_reg_codon
end interface

!===============================================================================
contains
!===============================================================================

subroutine modal_aero_calcsize_reg()
  use physics_buffer,   only: pbuf_add_field, dtype_r8
  use rad_constituents, only: rad_cnst_get_info

  integer :: nmodes
  integer(c_int64_t) :: active_c

  active_c = modal_aero_calcsize_reg_codon(1_c_int64_t)
  if (.not. modal_aero_calcsize_reg_logged) then
     modal_aero_calcsize_reg_logged = .true.
     if (masterproc) then
        write(iulog,'(A)') &
             'modal_aero_calcsize_reg direct = codon; rad constituent query/pbuf registration native CAM API island'
        call flush(iulog)
     end if
  end if
  if (active_c == 0_c_int64_t) return
  
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

subroutine modal_aero_calcsize_sub_append_impl_proof(env_name, proof_line)

  character(len=*), intent(in) :: env_name, proof_line

  character(len=512) :: proof_path
  integer :: status, n, unit_id

  call get_environment_variable(env_name, value=proof_path, length=n, status=status)
  if (status /= 0 .or. n <= 0) return

  open(newunit=unit_id, file=trim(adjustl(proof_path(:n))), status='unknown', action='write', &
       position='append', iostat=status)
  if (status /= 0) return

  write(unit_id,'(A)') trim(proof_line)
  close(unit_id)

end subroutine modal_aero_calcsize_sub_append_impl_proof

!===============================================================================
!===============================================================================

subroutine modal_aero_calcsize_sub_select_impl()

  character(len=48) :: impl_name
  integer :: status, n, i, code

  if (modal_aero_calcsize_sub_impl_selected) return

  impl_name = 'codon'
  call cam_codon_get_impl('MODAL_AERO_CALCSIZE_SUB_IMPL', impl_name, n, status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     modal_aero_calcsize_sub_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     modal_aero_calcsize_sub_use_native_impl = .false.
  end if

  modal_aero_calcsize_sub_impl_selected = .true.

  if (masterproc) then
     if (modal_aero_calcsize_sub_use_native_impl) then
        write(iulog,*) 'modal_aero_calcsize_sub implementation = native'
     else
        write(iulog,*) 'modal_aero_calcsize_sub implementation = codon'
        if (.not. modal_aero_calcsize_sub_proof_written) then
           call modal_aero_calcsize_sub_append_impl_proof('MODAL_AERO_CALCSIZE_SUB_PROOF_FILE', &
                'modal_aero_calcsize_sub selector entered implementation = codon')
           modal_aero_calcsize_sub_proof_written = .true.
        end if
     end if
     call flush(iulog)
  end if

end subroutine modal_aero_calcsize_sub_select_impl

!===============================================================================
!===============================================================================

subroutine modal_aero_calcsize_sub(state, ptend, deltat, pbuf, do_adjust_in, &
   do_aitacc_transfer_in)

   use iso_c_binding, only: c_int64_t

   ! arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend), target, intent(inout) :: ptend
   real(r8),                    intent(in)    :: deltat
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   logical, optional :: do_adjust_in
   logical, optional :: do_aitacc_transfer_in

#ifdef MODAL_AERO

   integer, parameter :: nsrflx = 4

   logical :: do_adjust
   logical :: do_aitacc_transfer

   integer :: iq, jac, l, lchnk, lsfrm, lstoo, n, ncol

   logical,  pointer :: dotend(:)
   real(r8), pointer :: dqdt(:,:,:)
   real(r8), pointer :: dgncur_a(:,:,:)
   real(r8), pointer :: fldcw(:,:)
   real(r8), pointer :: pdel(:,:)
   real(r8), pointer :: q(:,:,:)

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname

   real(r8), target :: qqcw_work(pcols,pver,pcnst)
   real(r8), target :: dqqcwdt_work(pcols,pver,pcnst)
   real(r8), target :: qsrflx_work(pcols,pcnst,nsrflx,2)
   real(r8), target :: dryvol_a_work(pcols,pver)
   real(r8), target :: dryvol_c_work(pcols,pver)
   real(r8), target :: drv_a_aitsv_work(pcols,pver), num_a_aitsv_work(pcols,pver)
   real(r8), target :: drv_c_aitsv_work(pcols,pver), num_c_aitsv_work(pcols,pver)
   real(r8), target :: drv_a_accsv_work(pcols,pver), num_a_accsv_work(pcols,pver)
   real(r8), target :: drv_c_accsv_work(pcols,pver), num_c_accsv_work(pcols,pver)
   real(r8), target :: dgnum_amode_work(ntot_amode)
   real(r8), target :: dgnumhi_amode_work(ntot_amode)
   real(r8), target :: dgnumlo_amode_work(ntot_amode)
   real(r8), target :: alnsg_amode_work(ntot_amode)
   real(r8), target :: voltonumb_amode_work(ntot_amode)
   real(r8), target :: voltonumblo_amode_work(ntot_amode)
   real(r8), target :: voltonumbhi_amode_work(ntot_amode)
   real(r8), target :: specdens_amode_work(maxd_aspectype)

   integer(c_int64_t) :: do_adjust_c, do_aitacc_transfer_c, nspecfrm_pair1_c
   integer(c_int64_t), target :: dotend_mask(pcnst), dotendqqcw_mask(pcnst)
   integer(c_int64_t), target :: mprognum_amode_c(ntot_amode)
   integer(c_int64_t), target :: numptr_amode_c(ntot_amode), numptrcw_amode_c(ntot_amode)
   integer(c_int64_t), target :: nspec_amode_c(ntot_amode)
   integer(c_int64_t), target :: lspectype_amode_c(maxd_aspectype,ntot_amode)
   integer(c_int64_t), target :: lmassptr_amode_c(maxd_aspectype,ntot_amode)
   integer(c_int64_t), target :: lmassptrcw_amode_c(maxd_aspectype,ntot_amode)
   integer(c_int64_t), target :: lspecfrma_pair1_c(maxspec_renamexf)
   integer(c_int64_t), target :: lspecfrmc_pair1_c(maxspec_renamexf)
   integer(c_int64_t), target :: lspectooa_pair1_c(maxspec_renamexf)
   integer(c_int64_t), target :: lspectooc_pair1_c(maxspec_renamexf)

   if (present(do_adjust_in)) then
      do_adjust = do_adjust_in
   else
      do_adjust = do_adjust_default
   end if

   if (present(do_aitacc_transfer_in)) then
      do_aitacc_transfer = do_aitacc_transfer_in
   else
      do_aitacc_transfer = do_aitacc_transfer_default
   end if

   call modal_aero_calcsize_sub_select_impl()
   if (modal_aero_calcsize_sub_use_native_impl) then
      call modal_aero_calcsize_sub_native(state, ptend, deltat, pbuf, do_adjust, do_aitacc_transfer)
      return
   end if

   if (do_aitacc_transfer) then
      if (npair_renamexf .le. 0) then
         write( 6, '(//a//)' ) '*** modal_aero_calcaersize_sub error -- npair_renamexf <= 0'
         call endrun( 'modal_aero_calcaersize_sub error' )
      end if
      if ((modefrm_renamexf(1) .ne. modeptr_aitken) .or. (modetoo_renamexf(1) .ne. modeptr_accum)) then
         write( 6, '(//2a//)' ) '*** modal_aero_calcaersize_sub error -- ', 'modefrm/too_renamexf(1) are wrong'
         call endrun( 'modal_aero_calcaersize_sub error' )
      end if
   end if

   lchnk = state%lchnk
   ncol  = state%ncol

   pdel => state%pdel
   q    => state%q
   dotend => ptend%lq
   dqdt   => ptend%q

   call pbuf_get_field(pbuf, dgnum_idx, dgncur_a)

   qqcw_work(:,:,:) = 0.0_r8
   do l = 1, pcnst
      nullify(fldcw)
      fldcw => qqcw_get_field(pbuf, l, lchnk, .true.)
      if (associated(fldcw)) then
         qqcw_work(1:ncol,top_lev:pver,l) = fldcw(1:ncol,top_lev:pver)
      end if
   end do

   dgnum_amode_work(:) = dgnum_amode(:)
   dgnumhi_amode_work(:) = dgnumhi_amode(:)
   dgnumlo_amode_work(:) = dgnumlo_amode(:)
   alnsg_amode_work(:) = alnsg_amode(:)
   voltonumb_amode_work(:) = voltonumb_amode(:)
   voltonumblo_amode_work(:) = voltonumblo_amode(:)
   voltonumbhi_amode_work(:) = voltonumbhi_amode(:)
   specdens_amode_work(:) = specdens_amode(:)

   do n = 1, ntot_amode
      mprognum_amode_c(n) = int(mprognum_amode(n), c_int64_t)
      numptr_amode_c(n) = int(numptr_amode(n), c_int64_t)
      numptrcw_amode_c(n) = int(numptrcw_amode(n), c_int64_t)
      nspec_amode_c(n) = int(nspec_amode(n), c_int64_t)
      do l = 1, maxd_aspectype
         lspectype_amode_c(l,n) = int(lspectype_amode(l,n), c_int64_t)
         lmassptr_amode_c(l,n) = int(lmassptr_amode(l,n), c_int64_t)
         lmassptrcw_amode_c(l,n) = int(lmassptrcw_amode(l,n), c_int64_t)
      end do
   end do

   nspecfrm_pair1_c = 0_c_int64_t
   lspecfrma_pair1_c(:) = 0_c_int64_t
   lspecfrmc_pair1_c(:) = 0_c_int64_t
   lspectooa_pair1_c(:) = 0_c_int64_t
   lspectooc_pair1_c(:) = 0_c_int64_t
   if (do_aitacc_transfer) then
      nspecfrm_pair1_c = int(nspecfrm_renamexf(1), c_int64_t)
      do iq = 1, nspecfrm_renamexf(1)
         lspecfrma_pair1_c(iq) = int(lspecfrma_renamexf(iq,1), c_int64_t)
         lspecfrmc_pair1_c(iq) = int(lspecfrmc_renamexf(iq,1), c_int64_t)
         lspectooa_pair1_c(iq) = int(lspectooa_renamexf(iq,1), c_int64_t)
         lspectooc_pair1_c(iq) = int(lspectooc_renamexf(iq,1), c_int64_t)
      end do
   end if

   do_adjust_c = merge(1_c_int64_t, 0_c_int64_t, do_adjust)
   do_aitacc_transfer_c = merge(1_c_int64_t, 0_c_int64_t, do_aitacc_transfer)
   dotend_mask(:) = 0_c_int64_t
   dotendqqcw_mask(:) = 0_c_int64_t

   call modal_aero_calcsize_sub_codon_wrap( &
        ncol, deltat, do_adjust_c, do_aitacc_transfer_c, q, qqcw_work, pdel, dqdt, dqqcwdt_work, qsrflx_work, &
        dgncur_a, dryvol_a_work, dryvol_c_work, drv_a_aitsv_work, num_a_aitsv_work, drv_c_aitsv_work, &
        num_c_aitsv_work, drv_a_accsv_work, num_a_accsv_work, drv_c_accsv_work, num_c_accsv_work, dotend_mask, &
        dotendqqcw_mask, mprognum_amode_c, numptr_amode_c, numptrcw_amode_c, nspec_amode_c, lspectype_amode_c, &
        lmassptr_amode_c, lmassptrcw_amode_c, dgnum_amode_work, dgnumhi_amode_work, dgnumlo_amode_work, &
        alnsg_amode_work, voltonumb_amode_work, voltonumblo_amode_work, voltonumbhi_amode_work, &
        specdens_amode_work, nspecfrm_pair1_c, lspecfrma_pair1_c, lspecfrmc_pair1_c, lspectooa_pair1_c, &
        lspectooc_pair1_c )

   do l = 1, pcnst
      if (dotend_mask(l) /= 0_c_int64_t) dotend(l) = .true.
      if (dotendqqcw_mask(l) == 0_c_int64_t) cycle
      nullify(fldcw)
      fldcw => qqcw_get_field(pbuf, l, lchnk, .true.)
      if (.not. associated(fldcw)) cycle
      fldcw(1:ncol,top_lev:pver) = qqcw_work(1:ncol,top_lev:pver,l)
   end do

   if ( .not. do_adjust ) return

   do n = 1, ntot_amode
      if (mprognum_amode(n) <= 0) cycle

      do jac = 1, 2
         if (jac == 1) then
            l = numptr_amode(n)
            tmpnamea = cnst_name(l)
         else
            l = numptrcw_amode(n)
            tmpnamea = cnst_name_cw(l)
         end if
         fieldname = trim(tmpnamea) // '_sfcsiz1'
         call outfld( fieldname, qsrflx_work(:,l,1,jac), pcols, lchnk)

         fieldname = trim(tmpnamea) // '_sfcsiz2'
         call outfld( fieldname, qsrflx_work(:,l,2,jac), pcols, lchnk)
      end do
   end do

   if ( .not. do_aitacc_transfer ) return

   do iq = 1, nspecfrm_renamexf(1)
      do jac = 1, 2
         if (jac .eq. 1) then
            lsfrm = lspecfrma_renamexf(iq,1)
            lstoo = lspectooa_renamexf(iq,1)
         else
            lsfrm = lspecfrmc_renamexf(iq,1)
            lstoo = lspectooc_renamexf(iq,1)
         end if
         if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

         if (jac .eq. 1) then
            tmpnamea = cnst_name(lsfrm)
            tmpnameb = cnst_name(lstoo)
         else
            tmpnamea = cnst_name_cw(lsfrm)
            tmpnameb = cnst_name_cw(lstoo)
         end if

         fieldname = trim(tmpnamea) // '_sfcsiz3'
         call outfld( fieldname, qsrflx_work(:,lsfrm,3,jac), pcols, lchnk)

         fieldname = trim(tmpnameb) // '_sfcsiz3'
         call outfld( fieldname, qsrflx_work(:,lstoo,3,jac), pcols, lchnk)

         fieldname = trim(tmpnamea) // '_sfcsiz4'
         call outfld( fieldname, qsrflx_work(:,lsfrm,4,jac), pcols, lchnk)

         fieldname = trim(tmpnameb) // '_sfcsiz4'
         call outfld( fieldname, qsrflx_work(:,lstoo,4,jac), pcols, lchnk)
      end do
   end do

#else

   call modal_aero_calcsize_sub_native(state, ptend, deltat, pbuf, do_adjust_in, do_aitacc_transfer_in)

#endif

end subroutine modal_aero_calcsize_sub

!===============================================================================

subroutine modal_aero_calcsize_sub_native(state, ptend, deltat, pbuf, do_adjust_in, &
   do_aitacc_transfer_in)

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

#ifdef MODAL_AERO

   ! local

   logical :: do_adjust
   logical :: do_aitacc_transfer

   integer  :: lchnk                ! chunk identifier
   integer  :: ncol                 ! number of columns

   real(r8), pointer :: t(:,:)      ! Temperature in Kelvin
   real(r8), pointer :: pmid(:,:)   ! pressure at model levels (Pa)
   real(r8), pointer :: pdel(:,:)   ! pressure thickness of levels
   real(r8), pointer :: q(:,:,:)    ! Tracer MR array 

   logical,  pointer :: dotend(:)   ! flag for doing tendency
   real(r8), pointer :: dqdt(:,:,:) ! TMR tendency array

   real(r8), pointer :: dgncur_a(:,:,:)

   integer  :: i, icol_diag, iduma, ipair, iq
   integer  :: ixfer_acc2ait, ixfer_ait2acc
   integer  :: ixfer_acc2ait_sv(pcols,pver), ixfer_ait2acc_sv(pcols,pver)
   integer  :: j, jac, jsrflx, k 
   integer  :: l, l1, la, lc, lna, lnc, lsfrm, lstoo
   integer  :: n, nacc, nait

   integer, save  :: idiagaa = 1

   logical  :: dotendqqcw(pcnst)
   logical  :: noxf_acc2ait(ntot_aspectype)

   character(len=fieldname_len)   :: tmpnamea, tmpnameb
   character(len=fieldname_len+3) :: fieldname

   real(r8), parameter :: third = 1.0_r8/3.0_r8
   real(r8), pointer :: fldcw(:,:)
   real(r8) :: delnum_a2, delnum_c2            !  work variables
   real(r8) :: delnum_a3, delnum_c3, delnum_t3 !  work variables
   real(r8) :: deltatinv                     ! 1/deltat
   real(r8) :: dgncur_c(pcols,pver,ntot_amode)
   real(r8) :: dgnyy, dgnxx                  ! dgnumlo/hi of current mode
   real(r8) :: dqqcwdt(pcols,pver,pcnst)     ! cloudborne TMR tendency array
   real(r8) :: drv_a, drv_c, drv_t           ! dry volume (cm3/mol_air)
   real(r8) :: drv_t0
   real(r8) :: drv_a_noxf, drv_c_noxf, drv_t_noxf 
   real(r8) :: drv_a_acc, drv_c_acc
   real(r8) :: drv_a_accsv(pcols,pver), drv_c_accsv(pcols,pver)
   real(r8) :: drv_a_aitsv(pcols,pver), drv_c_aitsv(pcols,pver)
   real(r8) :: drv_a_sv(pcols,pver,ntot_amode), drv_c_sv(pcols,pver,ntot_amode)
   real(r8) :: dryvol_a(pcols,pver)          ! interstital aerosol dry 
   ! volume (cm^3/mol_air)
   real(r8) :: dryvol_c(pcols,pver)          ! activated aerosol dry volume
   real(r8) :: duma, dumb, dumc, dumd        ! work variables
   real(r8) :: dumfac, dummwdens             ! work variables
   real(r8) :: frelaxadj                     ! relaxation factor applied
   ! to size bounds
   real(r8) :: fracadj                       ! deltat/tadj
   real(r8) :: num_a0, num_c0, num_t0        ! initial number (#/mol_air)
   real(r8) :: num_a1, num_c1                ! working number (#/mol_air)
   real(r8) :: num_a2, num_c2, num_t2        ! working number (#/mol_air)
   real(r8) :: num_a, num_c, num_t           ! final number (#/mol_air)
   real(r8) :: num_t_noxf
   real(r8) :: numbnd                        ! bounded number
   real(r8) :: num_a_acc, num_c_acc
   real(r8) :: num_a_accsv(pcols,pver), num_c_accsv(pcols,pver)
   real(r8) :: num_a_aitsv(pcols,pver), num_c_aitsv(pcols,pver)
   real(r8) :: num_a_sv(pcols,pver,ntot_amode), num_c_sv(pcols,pver,ntot_amode)
   real(r8) :: pdel_fac                      ! 
   real(r8) :: tadj                          ! adjustment time scale
   real(r8) :: tadjinv                       ! 1/tadj
   real(r8) :: v2ncur_a(pcols,pver,ntot_amode)
   real(r8) :: v2ncur_c(pcols,pver,ntot_amode)
   real(r8) :: v2nyy, v2nxx, v2nzz           ! voltonumblo/hi of current mode
   real(r8) :: v2nyyrl, v2nxxrl              ! relaxed voltonumblo/hi 
   real(r8) :: xfercoef
   real(r8) :: xfercoef_num_acc2ait, xfercoef_vol_acc2ait
   real(r8) :: xfercoef_num_ait2acc, xfercoef_vol_ait2acc
   real(r8) :: xferfrac_num_acc2ait, xferfrac_vol_acc2ait
   real(r8) :: xferfrac_num_ait2acc, xferfrac_vol_ait2acc
   real(r8) :: xfertend, xfertend_num(2,2)

   integer, parameter :: nsrflx = 4    ! last dimension of qsrflx
   real(r8) :: qsrflx(pcols,pcnst,nsrflx,2)
   ! process-specific column tracer tendencies
   ! 3rd index -- 
   !    1="standard" number adjust gain;
   !    2="standard" number adjust loss;
   !    3=aitken-->accum renaming; 4=accum-->aitken)
   ! 4th index -- 
   !    1="a" species; 2="c" species
   !-----------------------------------------------------------------------

   if (present(do_adjust_in)) then
      do_adjust = do_adjust_in
   else
      do_adjust = do_adjust_default
   end if

   if (present(do_aitacc_transfer_in)) then
      do_aitacc_transfer = do_aitacc_transfer_in
   else
      do_aitacc_transfer = do_aitacc_transfer_default
   end if

   lchnk = state%lchnk
   ncol  = state%ncol

   t    => state%t
   pmid => state%pmid
   pdel => state%pdel
   q    => state%q
      
   dotend => ptend%lq
   dqdt   => ptend%q

   call pbuf_get_field(pbuf, dgnum_idx, dgncur_a)

   dotendqqcw(:) = .false.
   dqqcwdt(:,:,:) = 0.0_r8
   qsrflx(:,:,:,:) = 0.0_r8

   nait = modeptr_aitken
   nacc = modeptr_accum

   deltatinv = 1.0_r8/(deltat*(1.0_r8 + 1.0e-15_r8))
   ! tadj = adjustment time scale for number, surface when they are prognosed
   !           currently set to deltat
   tadj = deltat
   tadj = 86400
   tadj = max( tadj, deltat )
   tadjinv = 1.0_r8/(tadj*(1.0_r8 + 1.0e-15_r8))
   fracadj = deltat*tadjinv
   fracadj = max( 0.0_r8, min( 1.0_r8, fracadj ) )

   
   !
   !
   ! the "do 40000" loop does the original (pre jan-2006)
   !   number adjustment, one mode at a time
   ! this artificially adjusts number when mean particle size is too large
   !   or too small
   !
   !
   do n = 1, ntot_amode


      ! initialize all parameters to the default values for the mode
      do k=top_lev,pver
         do i=1,ncol
            !    sgcur_a(i,k,n) = sigmag_amode(n)
            !    sgcur_c(i,k,n) = sigmag_amode(n)
            dgncur_a(i,k,n) = dgnum_amode(n)
            dgncur_c(i,k,n) = dgnum_amode(n)
            v2ncur_a(i,k,n) = voltonumb_amode(n)
            v2ncur_c(i,k,n) = voltonumb_amode(n)
            dryvol_a(i,k) = 0.0_r8
            dryvol_c(i,k) = 0.0_r8
         end do
      end do

      ! compute dry volume mixrats = 
      !      sum_over_components{ component_mass mixrat / density }
      do l1 = 1, nspec_amode(n)
         ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
         dummwdens = 1.0_r8 / specdens_amode(lspectype_amode(l1,n))
         la = lmassptr_amode(l1,n)
         do k=top_lev,pver
            do i=1,ncol
               dryvol_a(i,k) = dryvol_a(i,k)    &
                  + max(0.0_r8,q(i,k,la))*dummwdens
            end do
         end do

         fldcw => qqcw_get_field(pbuf,lmassptrcw_amode(l1,n),lchnk)
         do k=top_lev,pver
            do i=1,ncol
               dryvol_c(i,k) = dryvol_c(i,k)    &
                  + max(0.0_r8,fldcw(i,k))*dummwdens
            end do
         end do
      end do

      ! set "short-hand" number pointers
      lna = numptr_amode(n)
      lnc = numptrcw_amode(n)
      fldcw => qqcw_get_field(pbuf,numptrcw_amode(n),lchnk,.true.)


      ! go to section for appropriate number/surface diagnosed/prognosed options
      if (mprognum_amode(n) <= 0) then

         ! option 1 -- number diagnosed (fixed dgnum and sigmag)
         !    compute number tendencies that will bring numbers to their
         !    current diagnosed values
         !
         if (lna > 0) then
            dotend(lna) = .true.
            do k=top_lev,pver
               do i=1,ncol
                  dqdt(i,k,lna) = (dryvol_a(i,k)*voltonumb_amode(n)   &
                     - q(i,k,lna)) * deltatinv
               end do
            end do
         end if
         if (lnc > 0) then
            dotendqqcw(lnc) = .true.
            do k=top_lev,pver
               do i=1,ncol
                  dqqcwdt(i,k,lnc) = (dryvol_c(i,k)*voltonumb_amode(n)   &
                     - fldcw(i,k)) * deltatinv
               end do
            end do
         end if
      else


         !
         ! option 2 -- number prognosed (variable dgnum, fixed sigmag)
         !       Compute number tendencies to adjust numbers if they are outside
         !    the limits determined by current volume and dgnumlo/hi
         !       The interstitial and activated aerosol fractions can, at times,
         !    be the lower or upper tail of the "total" distribution.  Thus they
         !    can be expected to have a greater range of size parameters than
         !    what is specified for the total distribution (via dgnumlo/hi)
         !       When both the interstitial and activated dry volumes are positive,
         !    the adjustment strategy is to (1) adjust the interstitial and activated
         !    numbers towards relaxed bounds, then (2) adjust the total/combined
         !    number towards the primary bounds.
         !
         ! note
         !    v2nyy = voltonumblo_amode is proportional to dgnumlo**(-3), 
         !            and produces the maximum allowed number for a given volume
         !    v2nxx = voltonumbhi_amode is proportional to dgnumhi**(-3), 
         !            and produces the minimum allowed number for a given volume
         !    v2nxxrl and v2nyyrl are their "relaxed" equivalents.  
         !            Setting frelaxadj=27=3**3 means that 
         !            dgnumlo_relaxed = dgnumlo/3 and dgnumhi_relaxed = dgnumhi*3
         !
         ! if do_aitacc_transfer is .true., then
         !     for n=nacc, multiply v2nyy by 1.0e6 to effectively turn off the
         !         adjustment when number is too big (size is too small)
         !     for n=nait, divide   v2nxx by 1.0e6 to effectively turn off the
         !         adjustment when number is too small (size is too big)
         !OLD  however, do not change the v2nyyrl/v2nxxrl so that
         !OLD      the interstitial<-->activated adjustment is not changed
         !NEW  also change the v2nyyrl/v2nxxrl so that
         !NEW      the interstitial<-->activated adjustment is turned off 
         !
      end if
      frelaxadj = 27.0_r8
      dumfac = exp(4.5_r8*alnsg_amode(n)**2)*pi/6.0_r8
      v2nxx = voltonumbhi_amode(n)
      v2nyy = voltonumblo_amode(n)
      v2nxxrl = v2nxx/frelaxadj
      v2nyyrl = v2nyy*frelaxadj
      dgnxx = dgnumhi_amode(n)
      dgnyy = dgnumlo_amode(n)
      if ( do_aitacc_transfer ) then
         if (n == nait) v2nxx = v2nxx/1.0e6_r8
         if (n == nacc) v2nyy = v2nyy*1.0e6_r8
         v2nxxrl = v2nxx/frelaxadj   ! NEW
         v2nyyrl = v2nyy*frelaxadj   ! NEW
      end if

      if (do_adjust) then
         dotend(lna) = .true.
         dotendqqcw(lnc) = .true.
      end if

      do  k = top_lev, pver
         do  i = 1, ncol

            drv_a = dryvol_a(i,k)
            num_a0 = q(i,k,lna)
            num_a = max( 0.0_r8, num_a0 )
            drv_c = dryvol_c(i,k)
            num_c0 = fldcw(i,k)
            num_c = max( 0.0_r8, num_c0 )

            if ( do_adjust) then

               !
               ! do number adjustment for interstitial and activated particles
               !    adjustments that (1) make numbers non-negative or (2) make numbers
               !       zero when volume is zero are applied over time-scale deltat
               !    adjustments that bring numbers to within specified bounds are
               !       applied over time-scale tadj
               !
               if ((drv_a <= 0.0_r8) .and. (drv_c <= 0.0_r8)) then
                  ! both interstitial and activated volumes are zero
                  ! adjust both numbers to zero
                  num_a = 0.0_r8
                  dqdt(i,k,lna) = -num_a0*deltatinv
                  num_c = 0.0_r8
                  dqqcwdt(i,k,lnc) = -num_c0*deltatinv
               else if (drv_c <= 0.0_r8) then
                  ! activated volume is zero, so interstitial number/volume == total/combined
                  ! apply step 1 and 3, but skip the relaxed adjustment (step 2, see below)
                  num_c = 0.0_r8
                  dqqcwdt(i,k,lnc) = -num_c0*deltatinv
                  num_a1 = num_a
                  numbnd = max( drv_a*v2nxx, min( drv_a*v2nyy, num_a1 ) )
                  num_a  = num_a1 + (numbnd - num_a1)*fracadj
                  dqdt(i,k,lna) = (num_a - num_a0)*deltatinv

               else if (drv_a <= 0.0_r8) then
                  ! interstitial volume is zero, treat similar to above
                  num_a = 0.0_r8
                  dqdt(i,k,lna) = -num_a0*deltatinv
                  num_c1 = num_c
                  numbnd = max( drv_c*v2nxx, min( drv_c*v2nyy, num_c1 ) )
                  num_c  = num_c1 + (numbnd - num_c1)*fracadj
                  dqqcwdt(i,k,lnc) = (num_c - num_c0)*deltatinv
               else
                  ! both volumes are positive
                  ! apply 3 adjustment steps
                  ! step1:  num_a,c0 --> num_a,c1 forces non-negative values
                  num_a1 = num_a
                  num_c1 = num_c
                  ! step2:  num_a,c1 --> num_a,c2 applies relaxed bounds to the interstitial
                  !    and activated number (individually)
                  !    if only only a or c changes, adjust the other in the opposite direction
                  !    as much as possible to conserve a+c
                  numbnd = max( drv_a*v2nxxrl, min( drv_a*v2nyyrl, num_a1 ) )
                  delnum_a2 = (numbnd - num_a1)*fracadj
                  num_a2 = num_a1 + delnum_a2
                  numbnd = max( drv_c*v2nxxrl, min( drv_c*v2nyyrl, num_c1 ) )
                  delnum_c2 = (numbnd - num_c1)*fracadj
                  num_c2 = num_c1 + delnum_c2
                  if ((delnum_a2 == 0.0_r8) .and. (delnum_c2 /= 0.0_r8)) then
                     num_a2 = max( drv_a*v2nxxrl, min( drv_a*v2nyyrl,   &
                        num_a1-delnum_c2 ) )
                  else if ((delnum_a2 /= 0.0_r8) .and. (delnum_c2 == 0.0_r8)) then
                     num_c2 = max( drv_c*v2nxxrl, min( drv_c*v2nyyrl,   &
                        num_c1-delnum_a2 ) )
                  end if
                  ! step3:  num_a,c2 --> num_a,c3 applies stricter bounds to the 
                  !    combined/total number
                  drv_t = drv_a + drv_c
                  num_t2 = num_a2 + num_c2
                  delnum_a3 = 0.0_r8
                  delnum_c3 = 0.0_r8
                  if (num_t2 < drv_t*v2nxx) then
                     delnum_t3 = (drv_t*v2nxx - num_t2)*fracadj
                     ! if you are here then (num_a2 < drv_a*v2nxx) and/or
                     !                      (num_c2 < drv_c*v2nxx) must be true
                     if ((num_a2 < drv_a*v2nxx) .and. (num_c2 < drv_c*v2nxx)) then
                        delnum_a3 = delnum_t3*(num_a2/num_t2)
                        delnum_c3 = delnum_t3*(num_c2/num_t2)
                     else if (num_c2 < drv_c*v2nxx) then
                        delnum_c3 = delnum_t3
                     else if (num_a2 < drv_a*v2nxx) then
                        delnum_a3 = delnum_t3
                     end if
                  else if (num_t2 > drv_t*v2nyy) then
                     delnum_t3 = (drv_t*v2nyy - num_t2)*fracadj
                     ! if you are here then (num_a2 > drv_a*v2nyy) and/or
                     !                      (num_c2 > drv_c*v2nyy) must be true
                     if ((num_a2 > drv_a*v2nyy) .and. (num_c2 > drv_c*v2nyy)) then
                        delnum_a3 = delnum_t3*(num_a2/num_t2)
                        delnum_c3 = delnum_t3*(num_c2/num_t2)
                     else if (num_c2 > drv_c*v2nyy) then
                        delnum_c3 = delnum_t3
                     else if (num_a2 > drv_a*v2nyy) then
                        delnum_a3 = delnum_t3
                     end if
                  end if
                  num_a = num_a2 + delnum_a3
                  dqdt(i,k,lna) = (num_a - num_a0)*deltatinv
                  num_c = num_c2 + delnum_c3
                  dqqcwdt(i,k,lnc) = (num_c - num_c0)*deltatinv
               end if

            end if ! do_adjust

            !
            ! now compute current dgn and v2n
            !
            if (drv_a > 0.0_r8) then
               if (num_a <= drv_a*v2nxx) then
                  dgncur_a(i,k,n) = dgnxx
                  v2ncur_a(i,k,n) = v2nxx
               else if (num_a >= drv_a*v2nyy) then
                  dgncur_a(i,k,n) = dgnyy
                  v2ncur_a(i,k,n) = v2nyy
               else
                  dgncur_a(i,k,n) = (drv_a/(dumfac*num_a))**third
                  v2ncur_a(i,k,n) = num_a/drv_a
               end if
            end if
            pdel_fac = pdel(i,k)/gravit   ! = rho*dz
            jac = 1
            qsrflx(i,lna,1,jac) = qsrflx(i,lna,1,jac) + max(0.0_r8,dqdt(i,k,lna))*pdel_fac
            qsrflx(i,lna,2,jac) = qsrflx(i,lna,2,jac) + min(0.0_r8,dqdt(i,k,lna))*pdel_fac

            if (drv_c > 0.0_r8) then
               if (num_c <= drv_c*v2nxx) then
                  dgncur_c(i,k,n) = dgnumhi_amode(n)
                  v2ncur_c(i,k,n) = v2nxx
               else if (num_c >= drv_c*v2nyy) then
                  dgncur_c(i,k,n) = dgnumlo_amode(n)
                  v2ncur_c(i,k,n) = v2nyy
               else
                  dgncur_c(i,k,n) = (drv_c/(dumfac*num_c))**third
                  v2ncur_c(i,k,n) = num_c/drv_c
               end if
            end if
            jac = 2
            qsrflx(i,lnc,1,jac) = qsrflx(i,lnc,1,jac) + max(0.0_r8,dqqcwdt(i,k,lnc))*pdel_fac
            qsrflx(i,lnc,2,jac) = qsrflx(i,lnc,2,jac) + min(0.0_r8,dqqcwdt(i,k,lnc))*pdel_fac


            ! save number and dryvol for aitken <--> accum renaming
            if ( do_aitacc_transfer ) then
               if (n == nait) then
                  drv_a_aitsv(i,k) = drv_a
                  num_a_aitsv(i,k) = num_a
                  drv_c_aitsv(i,k) = drv_c
                  num_c_aitsv(i,k) = num_c
               else if (n == nacc) then
                  drv_a_accsv(i,k) = drv_a
                  num_a_accsv(i,k) = num_a
                  drv_c_accsv(i,k) = drv_c
                  num_c_accsv(i,k) = num_c
               end if
            end if
            drv_a_sv(i,k,n) = drv_a
            num_a_sv(i,k,n) = num_a
            drv_c_sv(i,k,n) = drv_c
            num_c_sv(i,k,n) = num_c

         end do
      end do


      !
      ! option 3 -- number and surface prognosed (variable dgnum and sigmag)
      !             this is not implemented
      !
   end do  ! do n = 1, ntot_amode


   !
   !
   ! the following section (from here to label 49000) 
   !    does aitken <--> accum mode transfer 
   !
   ! when the aitken mode mean size is too big, the largest
   !    aitken particles are transferred into the accum mode
   !    to reduce the aitken mode mean size
   ! when the accum mode mean size is too small, the smallest
   !    accum particles are transferred into the aitken mode
   !    to increase the accum mode mean size
   !
   !
   ixfer_ait2acc_sv(:,:) = 0
   ixfer_acc2ait_sv(:,:) = 0
   if ( do_aitacc_transfer ) then

      ! old - on time first step, npair_renamexf will be <= 0,
      !       in which case need to do modal_aero_rename_init
      ! new - init is now done through chem_init and things below it
      if (npair_renamexf .le. 0) then
         npair_renamexf = 0
         !        call modal_aero_rename_init
         if (npair_renamexf .le. 0) then
            write( 6, '(//a//)' )   &
               '*** modal_aero_calcaersize_sub error -- npair_renamexf <= 0'
            call endrun( 'modal_aero_calcaersize_sub error' )
         end if
      end if

      ! check that renaming ipair=1 is aitken-->accum
      ipair = 1
      if ((modefrm_renamexf(ipair) .ne. nait) .or.   &
         (modetoo_renamexf(ipair) .ne. nacc)) then
         write( 6, '(//2a//)' )   &
            '*** modal_aero_calcaersize_sub error -- ',   &
            'modefrm/too_renamexf(1) are wrong'
         call endrun( 'modal_aero_calcaersize_sub error' )
      end if

      ! set dotend() for species that will be transferred
      do iq = 1, nspecfrm_renamexf(ipair)
         lsfrm = lspecfrma_renamexf(iq,ipair)
         lstoo = lspectooa_renamexf(iq,ipair)
         if ((lsfrm > 0) .and. (lstoo > 0)) then
            dotend(lsfrm) = .true.
            dotend(lstoo) = .true.
         end if
         lsfrm = lspecfrmc_renamexf(iq,ipair)
         lstoo = lspectooc_renamexf(iq,ipair)
         if ((lsfrm > 0) .and. (lstoo > 0)) then
            dotendqqcw(lsfrm) = .true.
            dotendqqcw(lstoo) = .true.
         end if
      end do

      ! identify accum species cannot be transferred to aitken mode
      noxf_acc2ait(:) = .true.
      do l1 = 1, nspec_amode(nacc)
         la = lmassptr_amode(l1,nacc)
         do iq = 1, nspecfrm_renamexf(ipair)
            if (lspectooa_renamexf(iq,ipair) == la) then
               noxf_acc2ait(l1) = .false.
            end if
         end do
      end do

      ! v2nzz is voltonumb at the "geometrically-defined" mid-point
      ! between the aitken and accum modes
      v2nzz = sqrt(voltonumb_amode(nait)*voltonumb_amode(nacc))

      ! loop over columns and levels
      do  k = top_lev, pver
         do  i = 1, ncol

            pdel_fac = pdel(i,k)/gravit   ! = rho*dz
            xfertend_num(:,:) = 0.0_r8

            ! compute aitken --> accum transfer rates
            ixfer_ait2acc = 0
            xfercoef_num_ait2acc = 0.0_r8
            xfercoef_vol_ait2acc = 0.0_r8

            drv_t = drv_a_aitsv(i,k) + drv_c_aitsv(i,k)
            num_t = num_a_aitsv(i,k) + num_c_aitsv(i,k)
            if (drv_t > 0.0_r8) then
               if (num_t < drv_t*v2nzz) then
                  ixfer_ait2acc = 1
                  if (num_t < drv_t*voltonumb_amode(nacc)) then
                     xferfrac_num_ait2acc = 1.0_r8
                     xferfrac_vol_ait2acc = 1.0_r8
                  else
                     xferfrac_vol_ait2acc = ((num_t/drv_t) - v2nzz)/   &
                        (voltonumb_amode(nacc) - v2nzz)
                     xferfrac_num_ait2acc = xferfrac_vol_ait2acc*   &
                        (drv_t*voltonumb_amode(nacc)/num_t)
                     if ((xferfrac_num_ait2acc <= 0.0_r8) .or.   &
                        (xferfrac_vol_ait2acc <= 0.0_r8)) then
                        xferfrac_num_ait2acc = 0.0_r8
                        xferfrac_vol_ait2acc = 0.0_r8
                     else if ((xferfrac_num_ait2acc >= 1.0_r8) .or.   &
                        (xferfrac_vol_ait2acc >= 1.0_r8)) then
                        xferfrac_num_ait2acc = 1.0_r8
                        xferfrac_vol_ait2acc = 1.0_r8
                     end if
                  end if
                  xfercoef_num_ait2acc = xferfrac_num_ait2acc*tadjinv
                  xfercoef_vol_ait2acc = xferfrac_vol_ait2acc*tadjinv
                  xfertend_num(1,1) = num_a_aitsv(i,k)*xfercoef_num_ait2acc
                  xfertend_num(1,2) = num_c_aitsv(i,k)*xfercoef_num_ait2acc
               end if
            end if

            ! compute accum --> aitken transfer rates
            ! accum may have some species (seasalt, dust, poa, lll) that are
            !    not in aitken mode
            ! so first divide the accum drv & num into not-transferred (noxf) species 
            !    and transferred species, and use the transferred-species 
            !    portion in what follows
            ixfer_acc2ait = 0
            xfercoef_num_acc2ait = 0.0_r8
            xfercoef_vol_acc2ait = 0.0_r8

            drv_t = drv_a_accsv(i,k) + drv_c_accsv(i,k)
            num_t = num_a_accsv(i,k) + num_c_accsv(i,k)
            drv_a_noxf = 0.0_r8
            drv_c_noxf = 0.0_r8
            if (drv_t > 0.0_r8) then
               if (num_t > drv_t*v2nzz) then
                  do l1 = 1, nspec_amode(nacc)

                     if ( noxf_acc2ait(l1) ) then
                        ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
                        dummwdens = 1.0_r8 / specdens_amode(lspectype_amode(l1,nacc))
                        la = lmassptr_amode(l1,nacc)
                        drv_a_noxf = drv_a_noxf    &
                           + max(0.0_r8,q(i,k,la))*dummwdens
                        lc = lmassptrcw_amode(l1,nacc)
                        
                        fldcw => qqcw_get_field(pbuf,lmassptrcw_amode(l1,nacc),lchnk)
                        drv_c_noxf = drv_c_noxf    &
                           + max(0.0_r8,fldcw(i,k))*dummwdens
                     end if
                  end do
                  drv_t_noxf = drv_a_noxf + drv_c_noxf
                  num_t_noxf = drv_t_noxf*voltonumblo_amode(nacc)
                  num_t0 = num_t
                  drv_t0 = drv_t
                  num_t = max( 0.0_r8, num_t - num_t_noxf )
                  drv_t = max( 0.0_r8, drv_t - drv_t_noxf )
               end if
            end if

            if (drv_t > 0.0_r8) then
               if (num_t > drv_t*v2nzz) then
                  ixfer_acc2ait = 1
                  if (num_t > drv_t*voltonumb_amode(nait)) then
                     xferfrac_num_acc2ait = 1.0_r8
                     xferfrac_vol_acc2ait = 1.0_r8
                  else
                     xferfrac_vol_acc2ait = ((num_t/drv_t) - v2nzz)/   &
                        (voltonumb_amode(nait) - v2nzz)
                     xferfrac_num_acc2ait = xferfrac_vol_acc2ait*   &
                        (drv_t*voltonumb_amode(nait)/num_t)
                     if ((xferfrac_num_acc2ait <= 0.0_r8) .or.   &
                        (xferfrac_vol_acc2ait <= 0.0_r8)) then
                        xferfrac_num_acc2ait = 0.0_r8
                        xferfrac_vol_acc2ait = 0.0_r8
                     else if ((xferfrac_num_acc2ait >= 1.0_r8) .or.   &
                        (xferfrac_vol_acc2ait >= 1.0_r8)) then
                        xferfrac_num_acc2ait = 1.0_r8
                        xferfrac_vol_acc2ait = 1.0_r8
                     end if
                  end if
                  duma = 1.0e-37_r8
                  xferfrac_num_acc2ait = xferfrac_num_acc2ait*   &
                     num_t/max( duma, num_t0 )
                  xfercoef_num_acc2ait = xferfrac_num_acc2ait*tadjinv
                  xfercoef_vol_acc2ait = xferfrac_vol_acc2ait*tadjinv
                  xfertend_num(2,1) = num_a_accsv(i,k)*xfercoef_num_acc2ait
                  xfertend_num(2,2) = num_c_accsv(i,k)*xfercoef_num_acc2ait
               end if
            end if

            ! jump to end-of-loop if no transfer is needed at current i,k
            if (ixfer_ait2acc+ixfer_acc2ait > 0) then
               ixfer_ait2acc_sv(i,k) = ixfer_ait2acc
               ixfer_acc2ait_sv(i,k) = ixfer_acc2ait

               !
               ! compute new dgncur & v2ncur for aitken & accum modes
               !
               ! currently inactive
               do n = nait, nacc, (nacc-nait)
                  if (n .eq. nait) then
                     duma = (xfertend_num(1,1) - xfertend_num(2,1))*deltat
                     num_a     = max( 0.0_r8, num_a_aitsv(i,k) - duma )
                     num_a_acc = max( 0.0_r8, num_a_accsv(i,k) + duma )
                     duma = (drv_a_aitsv(i,k)*xfercoef_vol_ait2acc -   &
                        (drv_a_accsv(i,k)-drv_a_noxf)*xfercoef_vol_acc2ait)*deltat
                     drv_a     = max( 0.0_r8, drv_a_aitsv(i,k) - duma )
                     drv_a_acc = max( 0.0_r8, drv_a_accsv(i,k) + duma )
                     duma = (xfertend_num(1,2) - xfertend_num(2,2))*deltat
                     num_c     = max( 0.0_r8, num_c_aitsv(i,k) - duma )
                     num_c_acc = max( 0.0_r8, num_c_accsv(i,k) + duma )
                     duma = (drv_c_aitsv(i,k)*xfercoef_vol_ait2acc -   &
                        (drv_c_accsv(i,k)-drv_c_noxf)*xfercoef_vol_acc2ait)*deltat
                     drv_c     = max( 0.0_r8, drv_c_aitsv(i,k) - duma )
                     drv_c_acc = max( 0.0_r8, drv_c_accsv(i,k) + duma )
                  else
                     num_a = num_a_acc
                     drv_a = drv_a_acc
                     num_c = num_c_acc
                     drv_c = drv_c_acc
                  end if

                  if (drv_a > 0.0_r8) then
                     if (num_a <= drv_a*voltonumbhi_amode(n)) then
                        dgncur_a(i,k,n) = dgnumhi_amode(n)
                        v2ncur_a(i,k,n) = voltonumbhi_amode(n)
                     else if (num_a >= drv_a*voltonumblo_amode(n)) then
                        dgncur_a(i,k,n) = dgnumlo_amode(n)
                        v2ncur_a(i,k,n) = voltonumblo_amode(n)
                     else
                        dgncur_a(i,k,n) = (drv_a/(dumfac*num_a))**third
                        v2ncur_a(i,k,n) = num_a/drv_a
                     end if
                  else
                     dgncur_a(i,k,n) = dgnum_amode(n)
                     v2ncur_a(i,k,n) = voltonumb_amode(n)
                  end if
                  
                  if (drv_c > 0.0_r8) then
                     if (num_c <= drv_c*voltonumbhi_amode(n)) then
                        dgncur_c(i,k,n) = dgnumhi_amode(n)
                        v2ncur_c(i,k,n) = voltonumbhi_amode(n)
                     else if (num_c >= drv_c*voltonumblo_amode(n)) then
                        dgncur_c(i,k,n) = dgnumlo_amode(n)
                        v2ncur_c(i,k,n) = voltonumblo_amode(n)
                     else
                        dgncur_c(i,k,n) = (drv_c/(dumfac*num_c))**third
                        v2ncur_c(i,k,n) = num_c/drv_c
                     end if
                  else
                     dgncur_c(i,k,n) = dgnum_amode(n)
                     v2ncur_c(i,k,n) = voltonumb_amode(n)
                  end if

               end do


               !
               ! compute tendency amounts for aitken <--> accum transfer
               !
               
               if ( masterproc ) then
                  if (idiagaa > 0) then
                     do j = 1, 2
                        do iq = 1, nspecfrm_renamexf(ipair)
                           do jac = 1, 2
                              if (j .eq. 1) then
                                 if (jac .eq. 1) then
                                    lsfrm = lspecfrma_renamexf(iq,ipair)
                                    lstoo = lspectooa_renamexf(iq,ipair)
                                 else
                                    lsfrm = lspecfrmc_renamexf(iq,ipair)
                                    lstoo = lspectooc_renamexf(iq,ipair)
                                 end if
                              else
                                 if (jac .eq. 1) then
                                    lsfrm = lspectooa_renamexf(iq,ipair)
                                    lstoo = lspecfrma_renamexf(iq,ipair)
                                 else
                                    lsfrm = lspectooc_renamexf(iq,ipair)
                                    lstoo = lspecfrmc_renamexf(iq,ipair)
                                 end if
                              end if
                              write( 6, '(a,3i3,2i4)' ) 'calcsize j,iq,jac, lsfrm,lstoo',   &
                                 j,iq,jac, lsfrm,lstoo
                           end do
                        end do
                     end do
                  end if
               end if
               idiagaa = -1


               ! j=1 does aitken-->accum; j=2 does accum-->aitken 
               do  j = 1, 2

                  if ((j .eq. 1 .and. ixfer_ait2acc > 0) .or. &
                     (j .eq. 2 .and. ixfer_acc2ait > 0)) then

                     jsrflx = j+2
                     if (j .eq. 1) then
                        xfercoef = xfercoef_vol_ait2acc
                     else
                        xfercoef = xfercoef_vol_acc2ait
                     end if

                     do  iq = 1, nspecfrm_renamexf(ipair)

                        ! jac=1 does interstitial ("_a"); jac=2 does activated ("_c"); 
                        do  jac = 1, 2

                           ! the lspecfrma_renamexf (and lspecfrmc_renamexf) are aitken species
                           ! the lspectooa_renamexf (and lspectooc_renamexf) are accum  species
                           ! for j=1, want lsfrm=aitken species, lstoo=accum  species
                           ! for j=2, want lsfrm=accum  species,  lstoo=aitken species
                           if (j .eq. 1) then
                              if (jac .eq. 1) then
                                 lsfrm = lspecfrma_renamexf(iq,ipair)
                                 lstoo = lspectooa_renamexf(iq,ipair)
                              else
                                 lsfrm = lspecfrmc_renamexf(iq,ipair)
                                 lstoo = lspectooc_renamexf(iq,ipair)
                              end if
                           else
                              if (jac .eq. 1) then
                                 lsfrm = lspectooa_renamexf(iq,ipair)
                                 lstoo = lspecfrma_renamexf(iq,ipair)
                              else
                                 lsfrm = lspectooc_renamexf(iq,ipair)
                                 lstoo = lspecfrmc_renamexf(iq,ipair)
                              end if
                           end if

                           if ((lsfrm > 0) .and. (lstoo > 0)) then
                              if (jac .eq. 1) then
                                 if (iq .eq. 1) then
                                    xfertend = xfertend_num(j,jac)
                                 else
                                    xfertend = max(0.0_r8,q(i,k,lsfrm))*xfercoef
                                 end if
                                 dqdt(i,k,lsfrm) = dqdt(i,k,lsfrm) - xfertend
                                 dqdt(i,k,lstoo) = dqdt(i,k,lstoo) + xfertend
                              else
                                 if (iq .eq. 1) then
                                    xfertend = xfertend_num(j,jac)
                                 else
                                    fldcw => qqcw_get_field(pbuf,lsfrm,lchnk)
                                    xfertend = max(0.0_r8,fldcw(i,k))*xfercoef
                                 end if
                                 dqqcwdt(i,k,lsfrm) = dqqcwdt(i,k,lsfrm) - xfertend
                                 dqqcwdt(i,k,lstoo) = dqqcwdt(i,k,lstoo) + xfertend
                              end if
                              qsrflx(i,lsfrm,jsrflx,jac) = qsrflx(i,lsfrm,jsrflx,jac) - xfertend*pdel_fac
                              qsrflx(i,lstoo,jsrflx,jac) = qsrflx(i,lstoo,jsrflx,jac) + xfertend*pdel_fac
                           end if

                        end do
                     end do
                  end if
               end do

            end if
         end do
      end do


   end if  !  do_aitacc_transfer 
   lsfrm = -123456789   ! executable statement for debugging


   !
   ! apply tendencies to cloud-borne species MRs
   !
   do l = 1, pcnst
      lc = l
      if ( lc>0 .and. dotendqqcw(lc) ) then
         fldcw=> qqcw_get_field(pbuf,l,lchnk)
         do k = top_lev, pver
            do i = 1, ncol
               fldcw(i,k) = max( 0.0_r8,   &
                  (fldcw(i,k) + dqqcwdt(i,k,lc)*deltat) )
            end do
         end do
      end if
   end do

   !
   ! do outfld calls
   !

   ! history fields for number-adjust source-sink for all modes
   if ( .not. do_adjust ) return
   
   do n = 1, ntot_amode 
      if (mprognum_amode(n) <= 0) cycle

      do jac = 1, 2
         if (jac == 1) then
            l = numptr_amode(n)
            tmpnamea = cnst_name(l)
         else
            l = numptrcw_amode(n)
            tmpnamea = cnst_name_cw(l)
         end if
         fieldname = trim(tmpnamea) // '_sfcsiz1'
         call outfld( fieldname, qsrflx(:,l,1,jac), pcols, lchnk)
         
         fieldname = trim(tmpnamea) // '_sfcsiz2'
         call outfld( fieldname, qsrflx(:,l,2,jac), pcols, lchnk)
      end do   ! jac = ...

   end do   ! n = ...


   ! history fields for aitken-accum transfer
   if ( .not. do_aitacc_transfer ) return

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
         if ((lsfrm <= 0) .or. (lstoo <= 0)) cycle

         fieldname = trim(tmpnamea) // '_sfcsiz3'
         call outfld( fieldname, qsrflx(:,lsfrm,3,jac), pcols, lchnk)

         fieldname = trim(tmpnameb) // '_sfcsiz3'
         call outfld( fieldname, qsrflx(:,lstoo,3,jac), pcols, lchnk)

         fieldname = trim(tmpnamea) // '_sfcsiz4'
         call outfld( fieldname, qsrflx(:,lsfrm,4,jac), pcols, lchnk)

         fieldname = trim(tmpnameb) // '_sfcsiz4'
         call outfld( fieldname, qsrflx(:,lstoo,4,jac), pcols, lchnk)

      end do   ! jac = ...
   end do   ! iq = ...

#endif

end subroutine modal_aero_calcsize_sub_native

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine modal_aero_calcsize_sub_codon_wrap( &
     ncol, deltat, do_adjust_c, do_aitacc_transfer_c, q, qqcw, pdel, dqdt, dqqcwdt, qsrflx, dgncur_a, &
     dryvol_a, dryvol_c, drv_a_aitsv, num_a_aitsv, drv_c_aitsv, num_c_aitsv, drv_a_accsv, num_a_accsv, &
     drv_c_accsv, num_c_accsv, dotend_mask, dotendqqcw_mask, mprognum_amode_c, numptr_amode_c, &
     numptrcw_amode_c, nspec_amode_c, lspectype_amode_c, lmassptr_amode_c, lmassptrcw_amode_c, &
     dgnum_amode_work, dgnumhi_amode_work, dgnumlo_amode_work, alnsg_amode_work, voltonumb_amode_work, &
     voltonumblo_amode_work, voltonumbhi_amode_work, specdens_amode_work, nspecfrm_pair1_c, &
     lspecfrma_pair1_c, lspecfrmc_pair1_c, lspectooa_pair1_c, lspectooc_pair1_c )

  use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

  integer, intent(in) :: ncol
  real(r8), intent(in) :: deltat
  integer(c_int64_t), intent(in) :: do_adjust_c, do_aitacc_transfer_c, nspecfrm_pair1_c
  real(r8), target, intent(inout) :: q(pcols,pver,pcnst), qqcw(pcols,pver,pcnst)
  real(r8), target, intent(in) :: pdel(pcols,pver)
  real(r8), target, intent(inout) :: dqdt(pcols,pver,pcnst), dqqcwdt(pcols,pver,pcnst)
  real(r8), target, intent(inout) :: qsrflx(pcols,pcnst,4,2), dgncur_a(pcols,pver,ntot_amode)
  real(r8), target, intent(inout) :: dryvol_a(pcols,pver), dryvol_c(pcols,pver)
  real(r8), target, intent(inout) :: drv_a_aitsv(pcols,pver), num_a_aitsv(pcols,pver)
  real(r8), target, intent(inout) :: drv_c_aitsv(pcols,pver), num_c_aitsv(pcols,pver)
  real(r8), target, intent(inout) :: drv_a_accsv(pcols,pver), num_a_accsv(pcols,pver)
  real(r8), target, intent(inout) :: drv_c_accsv(pcols,pver), num_c_accsv(pcols,pver)
  integer(c_int64_t), target, intent(inout) :: dotend_mask(pcnst), dotendqqcw_mask(pcnst)
  integer(c_int64_t), target, intent(in) :: mprognum_amode_c(ntot_amode)
  integer(c_int64_t), target, intent(in) :: numptr_amode_c(ntot_amode), numptrcw_amode_c(ntot_amode)
  integer(c_int64_t), target, intent(in) :: nspec_amode_c(ntot_amode)
  integer(c_int64_t), target, intent(in) :: lspectype_amode_c(maxd_aspectype,ntot_amode)
  integer(c_int64_t), target, intent(in) :: lmassptr_amode_c(maxd_aspectype,ntot_amode)
  integer(c_int64_t), target, intent(in) :: lmassptrcw_amode_c(maxd_aspectype,ntot_amode)
  real(r8), target, intent(in) :: dgnum_amode_work(ntot_amode), dgnumhi_amode_work(ntot_amode)
  real(r8), target, intent(in) :: dgnumlo_amode_work(ntot_amode), alnsg_amode_work(ntot_amode)
  real(r8), target, intent(in) :: voltonumb_amode_work(ntot_amode), voltonumblo_amode_work(ntot_amode)
  real(r8), target, intent(in) :: voltonumbhi_amode_work(ntot_amode), specdens_amode_work(maxd_aspectype)
  integer(c_int64_t), target, intent(in) :: lspecfrma_pair1_c(maxspec_renamexf)
  integer(c_int64_t), target, intent(in) :: lspecfrmc_pair1_c(maxspec_renamexf)
  integer(c_int64_t), target, intent(in) :: lspectooa_pair1_c(maxspec_renamexf)
  integer(c_int64_t), target, intent(in) :: lspectooc_pair1_c(maxspec_renamexf)

  character(len=96) :: wrap_proof_line

  interface
     subroutine modal_aero_calcsize_sub_codon( &
          ncol_c, pcols_c, pver_c, pcnst_c, top_lev_c, ntot_amode_c, maxd_aspectype_c, maxspec_renamexf_c, &
          nait_c, nacc_c, do_adjust_c, do_aitacc_transfer_c, nspecfrm_pair1_c, deltat_c, gravit_c, q_p, qqcw_p, &
          pdel_p, dqdt_p, dqqcwdt_p, qsrflx_p, dgncur_a_p, dryvol_a_p, dryvol_c_p, drv_a_aitsv_p, num_a_aitsv_p, &
          drv_c_aitsv_p, num_c_aitsv_p, drv_a_accsv_p, num_a_accsv_p, drv_c_accsv_p, num_c_accsv_p, dotend_p, &
          dotendqqcw_p, mprognum_amode_p, numptr_amode_p, numptrcw_amode_p, nspec_amode_p, lspectype_amode_p, &
          lmassptr_amode_p, lmassptrcw_amode_p, dgnum_amode_p, dgnumhi_amode_p, dgnumlo_amode_p, alnsg_amode_p, &
          voltonumb_amode_p, voltonumblo_amode_p, voltonumbhi_amode_p, specdens_amode_p, lspecfrma_pair1_p, &
          lspecfrmc_pair1_p, lspectooa_pair1_p, lspectooc_pair1_p) bind(c, name="modal_aero_calcsize_sub_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, top_lev_c, ntot_amode_c
       integer(c_int64_t), value :: maxd_aspectype_c, maxspec_renamexf_c, nait_c, nacc_c
       integer(c_int64_t), value :: do_adjust_c, do_aitacc_transfer_c, nspecfrm_pair1_c
       real(c_double), value :: deltat_c, gravit_c
       type(c_ptr), value :: q_p, qqcw_p, pdel_p, dqdt_p, dqqcwdt_p, qsrflx_p, dgncur_a_p
       type(c_ptr), value :: dryvol_a_p, dryvol_c_p, drv_a_aitsv_p, num_a_aitsv_p, drv_c_aitsv_p, num_c_aitsv_p
       type(c_ptr), value :: drv_a_accsv_p, num_a_accsv_p, drv_c_accsv_p, num_c_accsv_p, dotend_p, dotendqqcw_p
       type(c_ptr), value :: mprognum_amode_p, numptr_amode_p, numptrcw_amode_p, nspec_amode_p
       type(c_ptr), value :: lspectype_amode_p, lmassptr_amode_p, lmassptrcw_amode_p
       type(c_ptr), value :: dgnum_amode_p, dgnumhi_amode_p, dgnumlo_amode_p, alnsg_amode_p
       type(c_ptr), value :: voltonumb_amode_p, voltonumblo_amode_p, voltonumbhi_amode_p, specdens_amode_p
       type(c_ptr), value :: lspecfrma_pair1_p, lspecfrmc_pair1_p, lspectooa_pair1_p, lspectooc_pair1_p
     end subroutine modal_aero_calcsize_sub_codon
  end interface

  if (masterproc .and. .not. modal_aero_calcsize_sub_wrap_proof_written) then
     wrap_proof_line = 'modal_aero_calcsize_sub_codon_wrap entered'
     write(iulog,'(A)') trim(wrap_proof_line)
     call modal_aero_calcsize_sub_append_impl_proof('MODAL_AERO_CALCSIZE_SUB_PROOF_FILE', trim(wrap_proof_line))
     modal_aero_calcsize_sub_wrap_proof_written = .true.
     call flush(iulog)
  end if

  call modal_aero_calcsize_sub_codon( &
       int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnst, c_int64_t), &
       int(top_lev, c_int64_t), int(ntot_amode, c_int64_t), int(maxd_aspectype, c_int64_t), &
       int(maxspec_renamexf, c_int64_t), int(modeptr_aitken, c_int64_t), int(modeptr_accum, c_int64_t), &
       do_adjust_c, do_aitacc_transfer_c, nspecfrm_pair1_c, real(deltat, c_double), real(gravit, c_double), &
       c_loc(q(1,1,1)), c_loc(qqcw(1,1,1)), c_loc(pdel(1,1)), c_loc(dqdt(1,1,1)), c_loc(dqqcwdt(1,1,1)), &
       c_loc(qsrflx(1,1,1,1)), c_loc(dgncur_a(1,1,1)), c_loc(dryvol_a(1,1)), c_loc(dryvol_c(1,1)), &
       c_loc(drv_a_aitsv(1,1)), c_loc(num_a_aitsv(1,1)), c_loc(drv_c_aitsv(1,1)), c_loc(num_c_aitsv(1,1)), &
       c_loc(drv_a_accsv(1,1)), c_loc(num_a_accsv(1,1)), c_loc(drv_c_accsv(1,1)), c_loc(num_c_accsv(1,1)), &
       c_loc(dotend_mask(1)), c_loc(dotendqqcw_mask(1)), c_loc(mprognum_amode_c(1)), c_loc(numptr_amode_c(1)), &
       c_loc(numptrcw_amode_c(1)), c_loc(nspec_amode_c(1)), c_loc(lspectype_amode_c(1,1)), &
       c_loc(lmassptr_amode_c(1,1)), c_loc(lmassptrcw_amode_c(1,1)), c_loc(dgnum_amode_work(1)), &
       c_loc(dgnumhi_amode_work(1)), c_loc(dgnumlo_amode_work(1)), c_loc(alnsg_amode_work(1)), &
       c_loc(voltonumb_amode_work(1)), c_loc(voltonumblo_amode_work(1)), c_loc(voltonumbhi_amode_work(1)), &
       c_loc(specdens_amode_work(1)), c_loc(lspecfrma_pair1_c(1)), c_loc(lspecfrmc_pair1_c(1)), &
       c_loc(lspectooa_pair1_c(1)), c_loc(lspectooc_pair1_c(1)) )

end subroutine modal_aero_calcsize_sub_codon_wrap


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
