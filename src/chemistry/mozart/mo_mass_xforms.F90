module mo_mass_xforms

  use ppgrid,       only : pcols, pver
  use shr_kind_mod, only : r8 => shr_kind_r8


  private
  public :: mmr2vmr, mmr2vmri, vmr2mmr, vmr2mmri, h2o_to_vmr, h2o_to_mmr, init_mass_xforms
  save

  real(r8) :: adv_mass_h2o = 18._r8
  logical :: mmr2vmr_use_native_impl = .false.
  logical :: mmr2vmr_impl_selected = .false.
  logical :: vmr2mmr_use_native_impl = .false.
  logical :: vmr2mmr_impl_selected = .false.
  logical :: h2o_to_vmr_use_native_impl = .false.
  logical :: h2o_to_vmr_impl_selected = .false.
  logical :: mass_xforms_batch_use_native_impl = .false.
  logical :: mass_xforms_batch_impl_selected = .false.
  logical :: mass_xforms_batch_entered_logged = .false.

contains

  subroutine mass_xforms_batch_append_proof(proof_line)

    character(len=*), intent(in) :: proof_line

    character(len=512) :: proof_file
    integer :: status, n, unitno

    proof_file = ''
    call get_environment_variable('MASS_XFORMS_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
    if (status == 0 .and. n > 0) then
       open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
       write(unitno,'(A)') trim(proof_line)
       close(unitno)
    end if

  end subroutine mass_xforms_batch_append_proof

  subroutine mass_xforms_batch_select_impl()

    use cam_logfile, only : iulog
    use spmd_utils,  only : masterproc

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (mass_xforms_batch_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('MASS_XFORMS_BATCH_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       mass_xforms_batch_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       mass_xforms_batch_use_native_impl = .false.
    end if

    mass_xforms_batch_impl_selected = .true.

    if (masterproc) then
       if (mass_xforms_batch_use_native_impl) then
          write(iulog,*) 'mass_xforms_batch implementation = native'
          call mass_xforms_batch_append_proof('mass_xforms_batch selector entered implementation = native')
       else
          write(iulog,*) 'mass_xforms_batch implementation = codon'
          call mass_xforms_batch_append_proof('mass_xforms_batch selector entered implementation = codon')
       end if
       call flush(iulog)
    end if

  end subroutine mass_xforms_batch_select_impl

  subroutine mass_xforms_batch_log_entered()

    use cam_logfile, only : iulog
    use spmd_utils,  only : masterproc

    if (mass_xforms_batch_entered_logged) return
    mass_xforms_batch_entered_logged = .true.

    if (masterproc) then
       write(iulog,'(A)') 'mass_xforms_batch entered (mmr2vmr/vmr2mmr/h2o_to_vmr direct = codon)'
       call mass_xforms_batch_append_proof('mass_xforms_batch entered (mmr2vmr/vmr2mmr/h2o_to_vmr direct = codon)')
       call flush(iulog)
    end if

  end subroutine mass_xforms_batch_log_entered

  subroutine init_mass_xforms
    use mo_chem_utls, only : get_spc_ndx
    use chem_mods,    only : adv_mass
    use mo_util,      only : chemistry_misc_codon_touch

    implicit none

    integer  :: id_h2o

    call chemistry_misc_codon_touch('init_mass_xforms', 173)

    id_h2o = get_spc_ndx('H2O')

    if ( id_h2o > 0 ) then
       adv_mass_h2o = adv_mass(id_h2o)
    else
       adv_mass_h2o = 18._r8
    endif

  endsubroutine init_mass_xforms

  subroutine mmr2vmr( mmr, vmr, mbar, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    use chem_mods, only : adv_mass, gas_pcnst
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    real(r8), target, intent(in)    :: mbar(ncol,pver)
    real(r8), target, intent(in)    :: mmr(pcols,pver,gas_pcnst)
    real(r8), target, intent(inout) :: vmr(ncol,pver,gas_pcnst)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), target :: adv_mass_local(gas_pcnst)

    interface
       subroutine mass_xforms_batch_mmr2vmr_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, mbar_p, mmr_p, adv_mass_p, vmr_p) &
            bind(c, name="mass_xforms_batch_mmr2vmr_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c
         type(c_ptr), value :: mbar_p, mmr_p, adv_mass_p, vmr_p
       end subroutine mass_xforms_batch_mmr2vmr_codon
    end interface

    call mass_xforms_batch_select_impl()

    if (mass_xforms_batch_use_native_impl) then
       do m = 1,gas_pcnst
          if( adv_mass(m) /= 0._r8 ) then
             do k = 1,pver
                vmr(:ncol,k,m) = mbar(:ncol,k) * mmr(:ncol,k,m) / adv_mass(m)
             end do
          end if
       end do
       return
    end if

    adv_mass_local(:) = adv_mass(:)

    call mass_xforms_batch_log_entered()
    call mass_xforms_batch_mmr2vmr_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         c_loc(mbar), c_loc(mmr), c_loc(adv_mass_local), c_loc(vmr) &
    )

  end subroutine mmr2vmr

  subroutine mmr2vmri( mmr, vmr, mbar, mi, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    real(r8), intent(in)    :: mi
    real(r8), intent(in)    :: mbar(:,:)
    real(r8), intent(in)    :: mmr(:,:)
    real(r8), intent(inout) :: vmr(:,:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer  :: k
    real(r8) :: rmi

    rmi = 1._r8/mi
    do k = 1,pver
       vmr(:ncol,k) = mbar(:ncol,k) * mmr(:ncol,k) * rmi
    end do

  end subroutine mmr2vmri

  subroutine vmr2mmr( vmr, mmr, mbar, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    use m_spc_id
    use chem_mods, only : adv_mass, gas_pcnst
    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    implicit none

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    real(r8), target, intent(in)    :: mbar(ncol,pver)
    real(r8), target, intent(in)    :: vmr(ncol,pver,gas_pcnst)
    real(r8), target, intent(inout) :: mmr(pcols,pver,gas_pcnst)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: k, m
    real(r8), target :: adv_mass_local(gas_pcnst)

    interface
       subroutine mass_xforms_batch_vmr2mmr_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, mbar_p, vmr_p, adv_mass_p, mmr_p) &
            bind(c, name="mass_xforms_batch_vmr2mmr_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c
         type(c_ptr), value :: mbar_p, vmr_p, adv_mass_p, mmr_p
       end subroutine mass_xforms_batch_vmr2mmr_codon
    end interface

    call mass_xforms_batch_select_impl()

    !-----------------------------------------------------------------
    !	... The non-group species
    !-----------------------------------------------------------------
    if (mass_xforms_batch_use_native_impl) then
       do m = 1,gas_pcnst
          if( adv_mass(m) /= 0._r8 ) then
             do k = 1,pver
                mmr(:ncol,k,m) = adv_mass(m) * vmr(:ncol,k,m) / mbar(:ncol,k)
             end do
          end if
       end do
       return
    end if

    adv_mass_local(:) = adv_mass(:)

    call mass_xforms_batch_log_entered()
    call mass_xforms_batch_vmr2mmr_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         c_loc(mbar), c_loc(vmr), c_loc(adv_mass_local), c_loc(mmr) &
    )

  end subroutine vmr2mmr

  subroutine vmr2mmri( vmr, mmr, mbar, mi, ncol )
    !-----------------------------------------------------------------
    !	... Xfrom from volume to mass mixing ratio
    !-----------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------
    !	... dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: ncol
    real(r8), intent(in)    :: mi
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(in)    :: vmr(ncol,pver)
    real(r8), intent(inout) :: mmr(pcols,pver)

    !-----------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------
    integer :: k, m

    !-----------------------------------------------------------------
    !	... mass to volume mixing for individual species
    !-----------------------------------------------------------------
    do k = 1,pver
       mmr(:ncol,k) = mi * vmr(:ncol,k) / mbar(:ncol,k)
    end do

  end subroutine vmr2mmri

  subroutine h2o_to_vmr( h2o_mmr, h2o_vmr, mbar, ncol )
    !-----------------------------------------------------------------------
    !     ... Transform water vapor from mass to volumetric mixing ratio
    !-----------------------------------------------------------------------

    use chem_mods, only : adv_mass
    use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr

    implicit none

    !-----------------------------------------------------------------------
    !	... Dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) ::    ncol
    real(r8), target, dimension(pcols,pver), intent(in) :: &
         h2o_mmr                ! specific humidity ( mmr )
    real(r8), target, dimension(ncol,pver), intent(in)  :: &
         mbar                   ! atmos mean mass
    real(r8), target, dimension(ncol,pver), intent(out) :: &
         h2o_vmr                ! water vapor vmr

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    integer ::   k

    interface
       subroutine mass_xforms_batch_h2o_to_vmr_codon(ncol_c, pcols_c, pver_c, adv_mass_h2o_c, h2o_mmr_p, mbar_p, h2o_vmr_p) &
            bind(c, name="mass_xforms_batch_h2o_to_vmr_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c
         real(c_double), value :: adv_mass_h2o_c
         type(c_ptr), value :: h2o_mmr_p, mbar_p, h2o_vmr_p
       end subroutine mass_xforms_batch_h2o_to_vmr_codon
    end interface

    call mass_xforms_batch_select_impl()

    if (mass_xforms_batch_use_native_impl) then
       do k = 1,pver
          h2o_vmr(:ncol,k) = mbar(:ncol,k) * h2o_mmr(:ncol,k) / adv_mass_h2o
       end do
       return
    end if

    call mass_xforms_batch_log_entered()
    call mass_xforms_batch_h2o_to_vmr_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), real(adv_mass_h2o, c_double), &
         c_loc(h2o_mmr), c_loc(mbar), c_loc(h2o_vmr) &
    )

  end subroutine h2o_to_vmr

  subroutine h2o_to_mmr( h2o_vmr, h2o_mmr, mbar, ncol )
    !-----------------------------------------------------------------------
    !     ... Transform water vapor from volumetric to mass mixing ratio
    !-----------------------------------------------------------------------

    use chem_mods, only : adv_mass

    implicit none

    !-----------------------------------------------------------------------
    !	... Dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) ::    ncol
    real(r8), dimension(ncol,pver), intent(in)  :: &
         mbar                   ! atmos mean mass
    real(r8), dimension(ncol,pver), intent(in)  :: &
         h2o_vmr               ! water vapor vmr
    real(r8), dimension(pcols,pver), intent(out) :: &
         h2o_mmr                ! specific humidity ( mmr )

    !-----------------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------------
    integer ::   k

    do k = 1,pver
       h2o_mmr(:ncol,k) = h2o_vmr(:ncol,k) * adv_mass_h2o / mbar(:ncol,k)
    end do

  end subroutine h2o_to_mmr

  subroutine mmr2vmr_select_impl()

    use cam_logfile, only : iulog
    use spmd_utils,  only : masterproc

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (mmr2vmr_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('MMR2VMR_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       mmr2vmr_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       mmr2vmr_use_native_impl = .false.
    end if

    mmr2vmr_impl_selected = .true.

    if (masterproc) then
       if (mmr2vmr_use_native_impl) then
          write(iulog,*) 'mmr2vmr implementation = native'
       else
          write(iulog,*) 'mmr2vmr implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine mmr2vmr_select_impl

  subroutine vmr2mmr_select_impl()

    use cam_logfile, only : iulog
    use spmd_utils,  only : masterproc

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (vmr2mmr_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('VMR2MMR_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       vmr2mmr_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       vmr2mmr_use_native_impl = .false.
    end if

    vmr2mmr_impl_selected = .true.

    if (masterproc) then
       if (vmr2mmr_use_native_impl) then
          write(iulog,*) 'vmr2mmr implementation = native'
       else
          write(iulog,*) 'vmr2mmr implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine vmr2mmr_select_impl

  subroutine h2o_to_vmr_select_impl()

    use cam_logfile, only : iulog
    use spmd_utils,  only : masterproc

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (h2o_to_vmr_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('H2O_TO_VMR_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       h2o_to_vmr_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       h2o_to_vmr_use_native_impl = .false.
    end if

    h2o_to_vmr_impl_selected = .true.

    if (masterproc) then
       if (h2o_to_vmr_use_native_impl) then
          write(iulog,*) 'h2o_to_vmr implementation = native'
       else
          write(iulog,*) 'h2o_to_vmr implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine h2o_to_vmr_select_impl

end module mo_mass_xforms
