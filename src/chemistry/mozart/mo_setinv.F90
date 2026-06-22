
module mo_setinv

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile,  only : iulog
  use chem_mods,    only : inv_lst, nfs, gas_pcnst
  use cam_history,  only : addfld, phys_decomp, outfld
  use ppgrid,       only : pcols, pver

  implicit none

  save

  integer :: id_o, id_o2, id_h
  integer :: m_ndx, o2_ndx, n2_ndx, h2o_ndx, o3_ndx
  logical :: has_o2, has_n2, has_h2o, has_o3, has_var_o2
  logical :: setinv_use_native_impl = .false.
  logical :: setinv_impl_selected = .false.
  logical :: setinv_postprocess_entered_logged = .false.

  private
  public :: setinv_inti, setinv, has_h2o, o2_ndx, h2o_ndx, n2_ndx

contains

  subroutine setinv_append_impl_proof(proof_line)

    character(len=*), intent(in) :: proof_line
    character(len=512) :: proof_file
    integer :: status, n, unitno

    proof_file = ''
    call get_environment_variable('SETINV_PROOF_FILE', value=proof_file, length=n, status=status)
    if (status == 0 .and. n > 0) then
       open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
       write(unitno,'(A)') trim(proof_line)
       close(unitno)
    end if

  end subroutine setinv_append_impl_proof

  subroutine setinv_log_postprocess_entered()

    use spmd_utils, only : masterproc

    if (setinv_postprocess_entered_logged) return
    setinv_postprocess_entered_logged = .true.

    if (masterproc) then
       write(iulog,*) 'setinv postprocess entered (unified tracer/output stage dispatch = codon)'
       call setinv_append_impl_proof('setinv postprocess entered (unified tracer/output stage dispatch = codon)')
       call flush(iulog)
    end if

  end subroutine setinv_log_postprocess_entered

  subroutine setinv_inti
    !-----------------------------------------------------------------
    !        ... initialize the module
    !-----------------------------------------------------------------

    use mo_chem_utls, only : get_inv_ndx, get_spc_ndx
    use spmd_utils,   only : masterproc
    use iso_c_binding, only : c_int64_t, c_loc

    implicit none

    integer :: i
    integer(c_int64_t), target :: lookup_ids(8)
    integer(c_int64_t), target :: ids_c(8)
    integer(c_int64_t), target :: flags_c(5)

    interface
       subroutine setinv_inti_codon(lookup_ids_p, ids_p, flags_p) bind(c, name="setinv_inti_codon")
         use iso_c_binding, only : c_ptr
         type(c_ptr), value :: lookup_ids_p, ids_p, flags_p
       end subroutine setinv_inti_codon
    end interface

    lookup_ids(1) = int(get_inv_ndx( 'M' ), c_int64_t)
    lookup_ids(2) = int(get_inv_ndx( 'N2' ), c_int64_t)
    lookup_ids(3) = int(get_inv_ndx( 'O2' ), c_int64_t)
    lookup_ids(4) = int(get_inv_ndx( 'H2O' ), c_int64_t)
    lookup_ids(5) = int(get_inv_ndx( 'O3' ), c_int64_t)

    lookup_ids(6) = int(get_spc_ndx('O'), c_int64_t)
    lookup_ids(7) = int(get_spc_ndx('O2'), c_int64_t)
    lookup_ids(8) = int(get_spc_ndx('H'), c_int64_t)
    ids_c(:) = 0_c_int64_t
    flags_c(:) = 0_c_int64_t

    call setinv_inti_codon(c_loc(lookup_ids), c_loc(ids_c), c_loc(flags_c))
    call setinv_inti_log_codon()

    m_ndx   = int(ids_c(1))
    n2_ndx  = int(ids_c(2))
    o2_ndx  = int(ids_c(3))
    h2o_ndx = int(ids_c(4))
    o3_ndx  = int(ids_c(5))

    id_o  = int(ids_c(6))
    id_o2 = int(ids_c(7))
    id_h  = int(ids_c(8))

    has_var_o2 = flags_c(1) /= 0_c_int64_t
    has_n2     = flags_c(2) /= 0_c_int64_t
    has_o2     = flags_c(3) /= 0_c_int64_t
    has_h2o    = flags_c(4) /= 0_c_int64_t
    has_o3     = flags_c(5) /= 0_c_int64_t

    if (masterproc) write(iulog,*) 'setinv_inti: m,n2,o2,h2o ndx = ',m_ndx,n2_ndx,o2_ndx,h2o_ndx

    do i = 1,nfs
      call addfld( trim(inv_lst(i))//'_dens', 'molecules/cm3', pver,'A', 'invariant density', phys_decomp )
      !call addfld( trim(inv_lst(i))//'_mmr', 'kg/kg', pver,'A', 'invariant density', phys_decomp )
      call addfld( trim(inv_lst(i))//'_vmr', 'mole/mole', pver,'A', 'invariant density', phys_decomp )
    enddo
      
  end subroutine setinv_inti

  subroutine setinv_inti_log_codon()

    use spmd_utils, only : masterproc

    implicit none

    if (masterproc) then
       write(iulog,*) 'setinv_inti implementation = codon'
       call flush(iulog)
    end if

  end subroutine setinv_inti_log_codon

  subroutine setinv( invariants, tfld, h2ovmr, vmr, pmid, ncol, lchnk, pbuf )
    !-----------------------------------------------------------------
    !        ... set the invariant densities (molecules/cm**3)
    !-----------------------------------------------------------------

    use mo_constants,  only : boltz_cgs
    use tracer_cnst,   only : num_tracer_cnst, tracer_cnst_flds, get_cnst_data
    use mo_chem_utls,  only : get_inv_ndx
    use physics_buffer, only : physics_buffer_desc
    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr

    implicit none

    !-----------------------------------------------------------------
    !        ... dummy arguments
    !-----------------------------------------------------------------
    integer,  intent(in)  ::      ncol                      ! chunk column count
    real(r8), target, intent(in)  ::      tfld(pcols,pver)          ! temperature
    real(r8), target, intent(in)  ::      h2ovmr(ncol,pver)         ! water vapor vmr
    real(r8), target, intent(in)  ::      pmid(pcols,pver)          ! pressure
    integer,  intent(in)  ::      lchnk                     ! chunk number
    real(r8), target, intent(in)  ::      vmr(ncol,pver,gas_pcnst)  ! vmr
    real(r8), target, intent(out) ::      invariants(ncol,pver,nfs) ! invariant array
    type(physics_buffer_desc), pointer :: pbuf(:)


    real(r8), target :: cnst_offline( ncol, pver )

    !-----------------------------------------------------------------
    !        .. local variables
    !-----------------------------------------------------------------
    integer :: k, i, ndx
    real(r8), parameter ::  Pa_xfac = 10._r8                 ! Pascals to dyne/cm^2
    real(r8) :: sum1(ncol)
    real(r8), target :: tmp_out(ncol,pver)
    real(r8), target :: tmp_vmr_out(ncol,pver)

    interface
       subroutine setinv_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, nfs_c, m_ndx_c, n2_ndx_c, o2_ndx_c, h2o_ndx_c, &
            id_o_c, id_o2_c, id_h_c, has_n2_c, has_o2_c, has_h2o_c, has_var_o2_c, pa_xfac_c, boltz_cgs_c, tfld_p, &
            h2ovmr_p, vmr_p, pmid_p, invariants_p) bind(c, name="setinv_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c, nfs_c
         integer(c_int64_t), value :: m_ndx_c, n2_ndx_c, o2_ndx_c, h2o_ndx_c
         integer(c_int64_t), value :: id_o_c, id_o2_c, id_h_c
         integer(c_int64_t), value :: has_n2_c, has_o2_c, has_h2o_c, has_var_o2_c
         real(c_double), value :: pa_xfac_c, boltz_cgs_c
         type(c_ptr), value :: tfld_p, h2ovmr_p, vmr_p, pmid_p, invariants_p
       end subroutine setinv_codon
       subroutine setinv_apply_tracer_cnst_stage_dispatch_codon(ncol_c, pver_c, nfs_c, ndx_c, m_ndx_c, cnst_offline_p, &
            invariants_p) bind(c, name="setinv_apply_tracer_cnst_stage_dispatch_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, nfs_c, ndx_c, m_ndx_c
         type(c_ptr), value :: cnst_offline_p, invariants_p
       end subroutine setinv_apply_tracer_cnst_stage_dispatch_codon
       subroutine setinv_copy_invariant_codon(ncol_c, pver_c, nfs_c, inv_ndx_c, invariants_p, tmp_out_p) &
            bind(c, name="setinv_copy_invariant_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, nfs_c, inv_ndx_c
         type(c_ptr), value :: invariants_p, tmp_out_p
       end subroutine setinv_copy_invariant_codon
       subroutine setinv_vmr_output_codon(ncol_c, pver_c, nfs_c, inv_ndx_c, m_ndx_c, invariants_p, tmp_out_p) &
            bind(c, name="setinv_vmr_output_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, nfs_c, inv_ndx_c, m_ndx_c
         type(c_ptr), value :: invariants_p, tmp_out_p
       end subroutine setinv_vmr_output_codon
       subroutine setinv_output_pair_stage_dispatch_codon(ncol_c, pver_c, nfs_c, inv_ndx_c, m_ndx_c, invariants_p, tmp_dens_p, tmp_vmr_p) &
            bind(c, name="setinv_output_pair_stage_dispatch_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, nfs_c, inv_ndx_c, m_ndx_c
         type(c_ptr), value :: invariants_p, tmp_dens_p, tmp_vmr_p
       end subroutine setinv_output_pair_stage_dispatch_codon
    end interface

    call setinv_select_impl()

    !-----------------------------------------------------------------
    !        note: invariants are in cgs density units.
    !              the pmid array is in pascals and must be
    !	       mutiplied by 10. to yield dynes/cm**2.
    !-----------------------------------------------------------------
    if (setinv_use_native_impl) then
       invariants(:,:,:) = 0._r8
       !-----------------------------------------------------------------
       !	... set m, n2, o2, and h2o densities
       !-----------------------------------------------------------------
       do k = 1,pver
          invariants(:ncol,k,m_ndx) = Pa_xfac * pmid(:ncol,k) / (boltz_cgs*tfld(:ncol,k))
       end do

       if( has_n2 ) then
          if ( has_var_o2 ) then
             do k = 1,pver
                sum1(:ncol) = (vmr(:ncol,k,id_o) + vmr(:ncol,k,id_o2) + vmr(:ncol,k,id_h))
                invariants(:ncol,k,n2_ndx) = (1._r8 - sum1(:)) * invariants(:ncol,k,m_ndx)
             end do
          else
             do k = 1,pver
                invariants(:ncol,k,n2_ndx) = .79_r8 * invariants(:ncol,k,m_ndx)
             end do
          endif
       end if
       if( has_o2 ) then
          do k = 1,pver
             invariants(:ncol,k,o2_ndx) = .21_r8 * invariants(:ncol,k,m_ndx)
          end do
       end if
       if( has_h2o ) then
          do k = 1,pver
             invariants(:ncol,k,h2o_ndx) = h2ovmr(:ncol,k) * invariants(:ncol,k,m_ndx)
          end do
       end if
    else
       call setinv_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
            int(nfs, c_int64_t), int(m_ndx, c_int64_t), int(n2_ndx, c_int64_t), int(o2_ndx, c_int64_t), &
            int(h2o_ndx, c_int64_t), int(id_o, c_int64_t), int(id_o2, c_int64_t), int(id_h, c_int64_t), &
            merge(1_c_int64_t, 0_c_int64_t, has_n2), merge(1_c_int64_t, 0_c_int64_t, has_o2), &
            merge(1_c_int64_t, 0_c_int64_t, has_h2o), merge(1_c_int64_t, 0_c_int64_t, has_var_o2), &
            real(Pa_xfac, c_double), real(boltz_cgs, c_double), c_loc(tfld), c_loc(h2ovmr), c_loc(vmr), c_loc(pmid), &
            c_loc(invariants) &
       )
    end if

    do i = 1,num_tracer_cnst

       call get_cnst_data( tracer_cnst_flds(i), cnst_offline,  ncol, lchnk, pbuf )
       ndx =  get_inv_ndx( tracer_cnst_flds(i) )

       if (setinv_use_native_impl) then
          do k = 1,pver
             invariants(:ncol,k,ndx) = cnst_offline(:ncol,k)*invariants(:ncol,k,m_ndx)
          enddo
       else
          call setinv_log_postprocess_entered()
          call setinv_apply_tracer_cnst_stage_dispatch_codon( &
               int(ncol, c_int64_t), int(pver, c_int64_t), int(nfs, c_int64_t), int(ndx, c_int64_t), &
               int(m_ndx, c_int64_t), c_loc(cnst_offline), c_loc(invariants) &
          )
       end if

    enddo

    do i = 1,nfs
      if (setinv_use_native_impl) then
        tmp_out(:ncol,:) =  invariants(:ncol,:,i)
      else
        call setinv_log_postprocess_entered()
        call setinv_output_pair_stage_dispatch_codon( &
             int(ncol, c_int64_t), int(pver, c_int64_t), int(nfs, c_int64_t), int(i, c_int64_t), int(m_ndx, c_int64_t), &
             c_loc(invariants), c_loc(tmp_out), c_loc(tmp_vmr_out) &
        )
      end if
      call outfld( trim(inv_lst(i))//'_dens', tmp_out(:ncol,:), ncol, lchnk )
      if (setinv_use_native_impl) then
        tmp_out(:ncol,:) =  invariants(:ncol,:,i) / invariants(:ncol,:,m_ndx)
        call outfld( trim(inv_lst(i))//'_vmr',  tmp_out(:ncol,:), ncol, lchnk )
      else
        call outfld( trim(inv_lst(i))//'_vmr',  tmp_vmr_out(:ncol,:), ncol, lchnk )
      end if
    enddo

  end subroutine setinv

  subroutine setinv_select_impl()

    use spmd_utils, only : masterproc

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (setinv_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('SETINV_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       setinv_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       setinv_use_native_impl = .false.
    end if

    setinv_impl_selected = .true.

    if (masterproc) then
       if (setinv_use_native_impl) then
          write(iulog,*) 'setinv implementation = native'
          call setinv_append_impl_proof('setinv selector entered implementation = native')
       else
          write(iulog,*) 'setinv implementation = codon'
          call setinv_append_impl_proof('setinv selector entered implementation = codon')
       end if
       call flush(iulog)
    end if

  end subroutine setinv_select_impl

end module mo_setinv
