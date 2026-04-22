
module mo_mean_mass

  implicit none

  private
  public :: set_mean_mass, init_mean_mass

  integer :: id_o2, id_o, id_h, id_n
  logical :: set_mean_mass_use_native_impl = .false.
  logical :: set_mean_mass_impl_selected = .false.

contains

  subroutine init_mean_mass
    use mo_chem_utls, only : get_spc_ndx

    implicit none

    id_o2 = get_spc_ndx('O2')
    id_o  = get_spc_ndx('O')
    id_h  = get_spc_ndx('H')
    id_n  = get_spc_ndx('N')

  endsubroutine init_mean_mass

  subroutine set_mean_mass( ncol, mmr, mbar )
    !-----------------------------------------------------------------
    !        ... Set the invariant densities (molecules/cm**3)
    !-----------------------------------------------------------------

    use shr_kind_mod,     only : r8 => shr_kind_r8
    use ppgrid,           only : pver, pcols
    use chem_mods,        only : adv_mass, gas_pcnst
    use physconst,        only : mwdry                   ! molecular weight of dry air
    use cam_abortutils,   only : endrun
    use phys_control,     only : waccmx_is               !WACCM-X runtime switch
    use iso_c_binding,    only : c_double, c_int64_t, c_loc, c_ptr

    implicit none

    !-----------------------------------------------------------------
    !        ... Dummy arguments
    !-----------------------------------------------------------------
    integer, intent(in)   ::      ncol
    real(r8), target, intent(in)  ::      mmr(pcols,pver,gas_pcnst) ! species concentrations (kg/kg)
    real(r8), target, intent(out) ::      mbar(ncol,pver)           ! mean mass (g/mole)

    !-----------------------------------------------------------------
    !        ... Local variables
    !-----------------------------------------------------------------
    integer  :: k
    real(r8) :: xn2(ncol)                                  ! n2 mmr
    real(r8) :: fn2(ncol)                                  ! n2 vmr
    real(r8) :: fo(ncol)                                   ! o  vmr
    real(r8) :: fo2(ncol)                                  ! o2 vmr
    real(r8) :: fh(ncol)                                   ! h vmr
    real(r8) :: ftot(ncol)                                 ! total vmr
    real(r8) :: mean_mass(ncol)                            ! wrk variable

    logical  :: fixed_mbar                                 ! Fixed mean mass flag
    real(r8), target :: adv_mass_local(gas_pcnst)

    interface
       subroutine set_mean_mass_codon(ncol_c, pcols_c, pver_c, gas_pcnst_c, id_o2_c, id_o_c, id_h_c, id_n_c, &
            fixed_mbar_c, mwdry_c, mmr_p, adv_mass_p, mbar_p) bind(c, name="set_mean_mass_codon")
         use iso_c_binding, only : c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, gas_pcnst_c
         integer(c_int64_t), value :: id_o2_c, id_o_c, id_h_c, id_n_c, fixed_mbar_c
         real(c_double), value :: mwdry_c
         type(c_ptr), value :: mmr_p, adv_mass_p, mbar_p
       end subroutine set_mean_mass_codon
    end interface

    call set_mean_mass_select_impl()

    !-------------------------------------------
    !  Mean mass not fixed for WACCM-X
    !-------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
      fixed_mbar = .false.
    else
      fixed_mbar = .true.
    endif

    if (set_mean_mass_use_native_impl) then
       if( fixed_mbar ) then
          !-----------------------------------------------------------------
          !	... use CAM meam molecular weight 
          !-----------------------------------------------------------------
          mbar(:ncol,:pver) = mwdry  
       else
          if ( id_o2 > 0 .and. id_o > 0 .and. id_h > 0 .and. id_n > 0 ) then
             !-----------------------------------------------------------------
             !	... set the mean mass
             !-----------------------------------------------------------------
             do k = 1,pver
                xn2(:)    = 1._r8 - (mmr(:ncol,k,id_o2) + mmr(:ncol,k,id_o) + mmr(:ncol,k,id_h))
                fn2(:)    = .5_r8 * xn2(:) / adv_mass(id_n)
                fo2(:)    = mmr(:ncol,k,id_o2) / adv_mass(id_o2)
                fo(:)     = mmr(:ncol,k,id_o) / adv_mass(id_o)
                fh(:)     = mmr(:ncol,k,id_h) / adv_mass(id_h)
                mbar(:ncol,k) = 1._r8 / (fn2(:) + fo2(:) + fo(:) + fh(:))
             end do
          else
             call endrun('set_mean_mass: not able to compute mean mass')
          endif
       endif
       return
    end if

    if (.not. fixed_mbar) then
       if (.not. (id_o2 > 0 .and. id_o > 0 .and. id_h > 0 .and. id_n > 0)) then
          call endrun('set_mean_mass: not able to compute mean mass')
       end if
    end if

    adv_mass_local(:) = adv_mass(:)

    call set_mean_mass_codon( &
         int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(gas_pcnst, c_int64_t), &
         int(id_o2, c_int64_t), int(id_o, c_int64_t), int(id_h, c_int64_t), int(id_n, c_int64_t), &
         merge(1_c_int64_t, 0_c_int64_t, fixed_mbar), real(mwdry, c_double), &
         c_loc(mmr), c_loc(adv_mass_local), c_loc(mbar) &
    )

  end subroutine set_mean_mass

  subroutine set_mean_mass_select_impl()

    use cam_logfile, only : iulog
    use spmd_utils,  only : masterproc

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (set_mean_mass_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('SET_MEAN_MASS_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       set_mean_mass_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       set_mean_mass_use_native_impl = .false.
    end if

    set_mean_mass_impl_selected = .true.

    if (masterproc) then
       if (set_mean_mass_use_native_impl) then
          write(iulog,*) 'set_mean_mass implementation = native'
       else
          write(iulog,*) 'set_mean_mass implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine set_mean_mass_select_impl

end module mo_mean_mass
