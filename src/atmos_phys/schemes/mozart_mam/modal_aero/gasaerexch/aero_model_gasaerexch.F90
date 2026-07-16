! Decoupled CAM modal-aerosol gas/aerosol process driver.
module ap_aero_model_gasaerexch_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8
  use constituents, only : pcnst
  use ppgrid, only : pcols, pver
  use physconst, only : gravit
  use chem_mods, only : gas_pcnst, adv_mass
  use mo_tracname, only : solsym
  use cam_history, only : outfld
  use perf_mod, only : t_startf, t_stopf
  use physics_buffer, only : physics_buffer_desc, pbuf_get_field
  use modal_aero_data, only : ntot_amode, cnst_name_cw, qqcw_get_field
  use time_manager, only : get_nstep
  use modal_aero_coag, only : modal_aero_coag_sub
  use modal_aero_gasaerexch, only : modal_aero_gasaerexch_sub
  use modal_aero_newnuc, only : modal_aero_newnuc_sub
  use mo_setsox, only : setsox, has_sox

  implicit none
  private

  public :: aero_model_gasaerexch_run

contains

  !> \section arg_table_aero_model_gasaerexch_run Argument Table
  !! \htmlinclude aero_model_gasaerexch_run.html
  subroutine aero_model_gasaerexch_run(                         &
       loffset, ncol, lchnk, troplev, delt, tfld, pmid, pdel,    &
       mbar, zm, qh2o, cwat, cldfr, cldnum, airdens, invariants, &
       del_h2so4_gasprod, vmr0, vmr, pbuf, dgnum_idx,            &
       dgnumwet_idx, wetdens_ap_idx, pblh_idx, sulfeq_idx,       &
       ndx_h2so4, dgnum_name, dgnumwet_name)

    integer, intent(in) :: loffset
    integer, intent(in) :: ncol
    integer, intent(in) :: lchnk
    integer, intent(in) :: troplev(:)
    real(r8), intent(in) :: delt
    real(r8), intent(in) :: tfld(:,:)
    real(r8), intent(in) :: pmid(:,:)
    real(r8), intent(in) :: pdel(:,:)
    real(r8), intent(in) :: mbar(:,:)
    real(r8), intent(in) :: airdens(:,:)
    real(r8), intent(in) :: invariants(:,:,:)
    real(r8), intent(in) :: del_h2so4_gasprod(:,:)
    real(r8), intent(in) :: zm(:,:)
    real(r8), intent(in) :: qh2o(:,:)
    real(r8), intent(in) :: cwat(:,:)
    real(r8), intent(in) :: cldfr(:,:)
    real(r8), intent(in) :: cldnum(:,:)
    real(r8), intent(in) :: vmr0(:,:,:)
    real(r8), intent(inout) :: vmr(:,:,:)
    type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:)
    integer, intent(in) :: dgnum_idx
    integer, intent(in) :: dgnumwet_idx
    integer, intent(in) :: wetdens_ap_idx
    integer, intent(in) :: pblh_idx
    integer, intent(in) :: sulfeq_idx
    integer, intent(in) :: ndx_h2so4
    character(len=*), intent(in) :: dgnum_name(:)
    character(len=*), intent(in) :: dgnumwet_name(:)

    integer :: n, m
    integer :: i, k
    integer :: nstep

    real(r8) :: del_h2so4_aeruptk(ncol,pver)

    real(r8), pointer :: dgnum(:,:,:), dgnumwet(:,:,:), wetdens(:,:,:)
    real(r8), pointer :: pblh(:)

    real(r8), dimension(ncol) :: wrk
    character(len=32) :: name
    real(r8) :: dvmrcwdt(ncol,pver,gas_pcnst)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw(ncol,pver,gas_pcnst)

    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: sulfeq(:,:,:)

    call pbuf_get_field(pbuf, dgnum_idx,      dgnum )
    call pbuf_get_field(pbuf, dgnumwet_idx,   dgnumwet )
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens )
    call pbuf_get_field(pbuf, pblh_idx,       pblh)

    do n=1,ntot_amode
       call outfld(dgnum_name(n), dgnum(1:ncol,1:pver,n), ncol, lchnk )
       call outfld(dgnumwet_name(n), dgnumwet(1:ncol,1:pver,n), ncol, lchnk )
    end do

    nstep = get_nstep()

    dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
    do m = 1, gas_pcnst
      wrk(:) = 0._r8
      do k = 1,pver
        wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m)*adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
      end do
      name = 'GS_'//trim(solsym(m))
      call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

    call qqcw2vmr( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)
    dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)

    if( has_sox ) then
       call setsox(   &
            ncol,     &
            lchnk,    &
            loffset,  &
            delt,     &
            pmid,     &
            pdel,     &
            tfld,     &
            mbar,     &
            cwat,     &
            cldfr,    &
            cldnum,   &
            airdens,  &
            invariants, &
            vmrcw,    &
            vmr       &
            )
    endif

    dvmrdt = (vmr - dvmrdt) / delt
    dvmrcwdt = (vmrcw - dvmrcwdt) / delt
    do m = 1, gas_pcnst
      wrk(:) = 0._r8
      do k = 1,pver
        wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,k,m) * adv_mass(m)/mbar(:ncol,k)*pdel(:ncol,k)/gravit
      end do
      name = 'AQ_'//trim(solsym(m))
      call outfld( name, wrk(:ncol), ncol, lchnk )
    enddo

    if (ndx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4)
    else
       del_h2so4_aeruptk(:,:) = 0.0_r8
    endif

    call t_startf('modal_gas-aer_exchng')

    if ( sulfeq_idx>0 ) then
       call pbuf_get_field( pbuf, sulfeq_idx, sulfeq )
    else
       nullify( sulfeq )
    endif

    call modal_aero_gasaerexch_sub(            &
         lchnk,    ncol,     nstep,            &
         loffset,            delt,             &
         tfld,     pmid,     pdel,             &
         qh2o,               troplev,          &
         vmr,                vmrcw,            &
         dvmrdt,             dvmrcwdt,         &
         dgnum,              dgnumwet,         &
         sulfeq     )

    if (ndx_h2so4 > 0) then
       del_h2so4_aeruptk(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - del_h2so4_aeruptk(1:ncol,:)
    endif

    call t_stopf('modal_gas-aer_exchng')

    call t_startf('modal_nucl')

    call modal_aero_newnuc_sub(                             &
         lchnk,    ncol,     nstep,            &
         loffset,            delt,             &
         tfld,     pmid,     pdel,             &
         zm,       pblh,                       &
         qh2o,     cldfr,                      &
         vmr,                                  &
         del_h2so4_gasprod,  del_h2so4_aeruptk )

    call t_stopf('modal_nucl')

    call t_startf('modal_coag')

    call modal_aero_coag_sub(                               &
         lchnk,    ncol,     nstep,            &
         loffset,            delt,             &
         tfld,     pmid,     pdel,             &
         vmr,                                  &
         dgnum,              dgnumwet,         &
         wetdens                          )

    call t_stopf('modal_coag')

    call vmr2qqcw( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    do n = 1,pcnst
       fldcw => qqcw_get_field(pbuf,n,lchnk,errorhandle=.true.)
       if(associated(fldcw)) then
          call outfld( cnst_name_cw(n), fldcw(:,:), pcols, lchnk )
       endif
    end do

  end subroutine aero_model_gasaerexch_run

  subroutine qqcw2vmr(lchnk, vmr, mbar, ncol, im, pbuf)

    integer, intent(in) :: lchnk, ncol, im
    real(r8), intent(in) :: mbar(ncol,pver)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: k, m
    real(r8), pointer :: fldcw(:,:)

    do m=1,gas_pcnst
       if( adv_mass(m) /= 0._r8 ) then
          fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
          if(associated(fldcw)) then
             do k=1,pver
                vmr(:ncol,k,m) = mbar(:ncol,k) * fldcw(:ncol,k) / adv_mass(m)
             end do
          else
             vmr(:,:,m) = 0.0_r8
          end if
       end if
    end do
  end subroutine qqcw2vmr

  subroutine vmr2qqcw( lchnk, vmr, mbar, ncol, im, pbuf )

    use m_spc_id

    integer, intent(in) :: lchnk, ncol, im
    real(r8), intent(in) :: mbar(ncol,pver)
    real(r8), intent(in) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: k, m
    real(r8), pointer :: fldcw(:,:)

    do m = 1,gas_pcnst
       fldcw => qqcw_get_field(pbuf, m+im,lchnk,errorhandle=.true.)
       if( adv_mass(m) /= 0._r8 .and. associated(fldcw)) then
          do k = 1,pver
             fldcw(:ncol,k) = adv_mass(m) * vmr(:ncol,k,m) / mbar(:ncol,k)
          end do
       end if
    end do

  end subroutine vmr2qqcw

end module ap_aero_model_gasaerexch_scheme
