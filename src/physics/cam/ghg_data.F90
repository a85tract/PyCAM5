
module ghg_data

!------------------------------------------------------------------------------------------------
! Purpose:
! Provide default distributions of CH4, N2O, CFC11 and CFC12 to the radiation routines.
! **NOTE** CO2 is assumed by the radiation to a be constant value.  This value is
!          currently supplied directly by the chem_surfvals module.
!
! Revision history:
! 2004-08-29  B. Eaton        Create CAM interface to trcmix.
!------------------------------------------------------------------------------------------------

use shr_kind_mod,   only: r8 => shr_kind_r8
use ppgrid,         only: pcols, pver, begchunk, endchunk
use physics_types,  only: physics_state
use physconst,      only: mwdry, mwch4, mwn2o, mwf11, mwf12, mwco2
use chem_surfvals,  only: chem_surfvals_get, chem_surfvals_co2_rad
use cam_abortutils, only: endrun
use error_messages, only: handle_err
use iso_c_binding, only: c_int64_t, c_ptr
use spmd_utils, only: masterproc
use cam_logfile, only: iulog


implicit none
private
save

! Public interfaces
public ::&
   ghg_data_register, &! register ghg's with pbuf2d
   ghg_data_timestep_init    ! place data model of ghg's in pbuf2d

! Private variables

real(r8), target :: rmwn2o ! = mwn2o/mwdry ! ratio of molecular weight n2o   to dry air
real(r8), target :: rmwch4 ! = mwch4/mwdry ! ratio of molecular weight ch4   to dry air
real(r8), target :: rmwf11 ! = mwf11/mwdry ! ratio of molecular weight cfc11 to dry air
real(r8), target :: rmwf12 ! = mwf12/mwdry ! ratio of molecular weight cfc12 to dry air
real(r8), target :: rmwco2 ! = mwco2/mwdry ! ratio of molecular weights of co2 to dry air

integer, parameter :: ncnst = 6                        ! number of constituents
character(len=8), dimension(ncnst), parameter :: &
   cnst_names = (/'N2O  ', 'CH4  ', 'CFC11', 'CFC12', 'CO2  ', 'O2   '/) ! constituent names
integer  :: pbuf_idx(ncnst)

logical :: use_native_ghg_data_mw_ratios_impl = .false.
logical :: ghg_data_mw_ratios_impl_selected = .false.
logical :: ghg_data_mw_ratios_proof_written = .false.
logical :: use_native_ghg_data_trcmix_scale_impl = .false.
logical :: ghg_data_trcmix_scale_impl_selected = .false.
logical :: ghg_data_trcmix_scale_proof_written = .false.
logical :: ghg_data_trcmix_direct_logged = .false.
logical :: ghg_data_timestep_init_direct_logged = .false.

interface
  subroutine ghg_data_mw_ratios_codon(mwdry_c, mwn2o_c, mwch4_c, mwf11_c, mwf12_c, mwco2_c, &
       rmwn2o_p, rmwch4_p, rmwf11_p, rmwf12_p, rmwco2_p) bind(c, name="ghg_data_mw_ratios_codon")
    use iso_c_binding, only: c_double, c_ptr
    real(c_double), value :: mwdry_c, mwn2o_c, mwch4_c, mwf11_c, mwf12_c, mwco2_c
    type(c_ptr), value :: rmwn2o_p, rmwch4_p, rmwf11_p, rmwf12_p, rmwco2_p
  end subroutine ghg_data_mw_ratios_codon

	  function ghg_data_trcmix_scale_codon(gas_id_c, dlat_c) result(scale_c) &
	       bind(c, name="ghg_data_trcmix_scale_codon")
	    use iso_c_binding, only: c_int64_t, c_double
	    integer(c_int64_t), value :: gas_id_c
	    real(c_double), value :: dlat_c
	    real(c_double) :: scale_c
	  end function ghg_data_trcmix_scale_codon
	  function ghg_data_register_codon(flag_c) result(out_c) bind(c, name="ghg_data_register_codon")
	    use iso_c_binding, only: c_int64_t
	    integer(c_int64_t), value :: flag_c
	    integer(c_int64_t) :: out_c
	  end function ghg_data_register_codon
	  function ghg_data_timestep_init_codon(flag_c) result(out_c) bind(c, name="ghg_data_timestep_init_codon")
	    use iso_c_binding, only: c_int64_t
	    integer(c_int64_t), value :: flag_c
	    integer(c_int64_t) :: out_c
	  end function ghg_data_timestep_init_codon
	  subroutine ghg_data_trcmix_codon(gas_id_c, ncol_c, pcols_c, pver_c, trop_mmr_c, constant_mmr_c, &
	       clat_p, pmid_p, q_p) bind(c, name="ghg_data_trcmix_codon")
	    use iso_c_binding, only: c_double, c_int64_t, c_ptr
	    integer(c_int64_t), value :: gas_id_c, ncol_c, pcols_c, pver_c
	    real(c_double), value :: trop_mmr_c, constant_mmr_c
	    type(c_ptr), value :: clat_p, pmid_p, q_p
	  end subroutine ghg_data_trcmix_codon
	end interface

!================================================================================================
contains
!================================================================================================

logical function ghg_data_use_native(selector)
  character(len=*), intent(in) :: selector
  character(len=32) :: impl_name
  integer :: status, n, i, code

  impl_name = 'codon'
  call get_environment_variable(selector, value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     ghg_data_use_native = trim(adjustl(impl_name(:n))) == 'native'
  else
     ghg_data_use_native = .false.
  end if
end function ghg_data_use_native

subroutine ghg_data_mw_ratios_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (ghg_data_mw_ratios_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GHG_DATA_MW_RATIOS_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
    do i = 1, n
      code = iachar(impl_name(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
    use_native_ghg_data_mw_ratios_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
    use_native_ghg_data_mw_ratios_impl = .false.
  end if

  ghg_data_mw_ratios_impl_selected = .true.

  if (masterproc) then
    if (use_native_ghg_data_mw_ratios_impl) then
      write(iulog,*) 'ghg_data_mw_ratios implementation = native'
    else
      write(iulog,*) 'ghg_data_mw_ratios implementation = codon'
    end if
  end if

end subroutine ghg_data_mw_ratios_select_impl

subroutine ghg_data_trcmix_scale_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (ghg_data_trcmix_scale_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('GHG_DATA_TRCMIX_SCALE_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
    do i = 1, n
      code = iachar(impl_name(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
    use_native_ghg_data_trcmix_scale_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
    use_native_ghg_data_trcmix_scale_impl = .false.
  end if

  ghg_data_trcmix_scale_impl_selected = .true.

  if (masterproc) then
    if (use_native_ghg_data_trcmix_scale_impl) then
      write(iulog,*) 'ghg_data_trcmix_scale implementation = native'
    else
      write(iulog,*) 'ghg_data_trcmix_scale implementation = codon'
    end if
  end if

end subroutine ghg_data_trcmix_scale_select_impl

subroutine ghg_data_mw_ratios_proof_once()

  if (ghg_data_mw_ratios_proof_written) return
  ghg_data_mw_ratios_proof_written = .true.

  if (masterproc) then
    write(iulog,'(A)') 'ghg_data_mw_ratios entered (greenhouse gas molecular-weight ratios = codon)'
  end if

end subroutine ghg_data_mw_ratios_proof_once

subroutine ghg_data_trcmix_scale_proof_once()

  if (ghg_data_trcmix_scale_proof_written) return
  ghg_data_trcmix_scale_proof_written = .true.

  if (masterproc) then
    write(iulog,'(A)') 'ghg_data_trcmix_scale entered (trace-gas latitude scale helper = codon)'
  end if

end subroutine ghg_data_trcmix_scale_proof_once

subroutine ghg_data_trcmix_direct_proof_once()

  if (ghg_data_trcmix_direct_logged) return
  ghg_data_trcmix_direct_logged = .true.

  if (masterproc) then
    write(iulog,'(A)') 'trcmix direct = codon'
    call flush(iulog)
  end if

end subroutine ghg_data_trcmix_direct_proof_once

subroutine ghg_data_timestep_init_direct_proof_once()

  if (ghg_data_timestep_init_direct_logged) return
  ghg_data_timestep_init_direct_logged = .true.

  if (masterproc) then
    write(iulog,'(A)') 'ghg_data_timestep_init direct = codon control shell; ' // &
         'native pbuf chunk API boundary; trcmix body direct = codon'
    call flush(iulog)
  end if

end subroutine ghg_data_timestep_init_direct_proof_once

subroutine ghg_data_update_mw_ratios()
!-------------------------------------------------------------------------------
! update molecular-weight ratios used by the greenhouse-gas profile helper
!-------------------------------------------------------------------------------
  use iso_c_binding, only: c_double, c_loc

  call ghg_data_mw_ratios_select_impl()

  if (use_native_ghg_data_mw_ratios_impl) then
    call ghg_data_update_mw_ratios_native()
    return
  end if

  call ghg_data_mw_ratios_proof_once()
  call ghg_data_mw_ratios_codon(real(mwdry, c_double), real(mwn2o, c_double), real(mwch4, c_double), &
       real(mwf11, c_double), real(mwf12, c_double), real(mwco2, c_double), c_loc(rmwn2o), c_loc(rmwch4), &
       c_loc(rmwf11), c_loc(rmwf12), c_loc(rmwco2))

end subroutine ghg_data_update_mw_ratios

subroutine ghg_data_update_mw_ratios_native()
!-------------------------------------------------------------------------------
! native fallback for molecular-weight ratios
!-------------------------------------------------------------------------------

  rmwn2o = mwn2o/mwdry      ! ratio of molecular weight n2o   to dry air
  rmwch4 = mwch4/mwdry      ! ratio of molecular weight ch4   to dry air
  rmwf11 = mwf11/mwdry      ! ratio of molecular weight cfc11 to dry air
  rmwf12 = mwf12/mwdry      ! ratio of molecular weight cfc12 to dry air
  rmwco2 = mwco2/mwdry      ! ratio of molecular weights of co2 to dry air

end subroutine ghg_data_update_mw_ratios_native

subroutine ghg_data_register()
!-------------------------------------------------------------------------------
! register ghg's with pbuf2d
!-------------------------------------------------------------------------------
  use physics_buffer, only : pbuf_add_field, dtype_r8

  integer iconst
  if (ghg_data_register_codon(1_c_int64_t) == 0_c_int64_t) return

 
  do iconst = 1,ncnst
     call pbuf_add_field(cnst_names(iconst),'physpkg',dtype_r8,(/pcols,pver/),pbuf_idx(iconst))
  enddo

end subroutine ghg_data_register

subroutine ghg_data_timestep_init(pbuf2d, state)
!-------------------------------------------------------------------------------
! place data model of ghg's in pbuf2d at each timestep
!-------------------------------------------------------------------------------
  use ppgrid,              only: begchunk, endchunk, pcols, pver
  use physics_types,       only: physics_state
  use physics_buffer,      only: physics_buffer_desc, pbuf_get_field, pbuf_get_chunk

  
  type(physics_state), intent(in), dimension(begchunk:endchunk) :: state
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
 
  type(physics_buffer_desc), pointer :: pbuf_chnk(:)
  real(r8), pointer :: tmpptr(:,:)

  integer iconst
  integer lchnk
  if (ghg_data_use_native('GHG_DATA_TIMESTEP_INIT_IMPL')) then
     if (masterproc .and. .not. ghg_data_timestep_init_direct_logged) then
        ghg_data_timestep_init_direct_logged = .true.
        write(iulog,'(A)') 'ghg_data_timestep_init direct = native'
        call flush(iulog)
     end if
  else
     if (ghg_data_timestep_init_codon(1_c_int64_t) == 0_c_int64_t) return
     call ghg_data_timestep_init_direct_proof_once()
  end if

  call ghg_data_update_mw_ratios()

   do iconst = 1,ncnst
!$OMP PARALLEL DO PRIVATE (LCHNK,tmpptr,pbuf_chnk)
     do lchnk = begchunk, endchunk
       pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
       call pbuf_get_field(pbuf_chnk, pbuf_idx(iconst), tmpptr) 
       call trcmix(cnst_names(iconst), state(lchnk)%ncol, &
                   state(lchnk)%lat, state(lchnk)%pmid, &
                   tmpptr)
     enddo
  enddo

end subroutine ghg_data_timestep_init


!================================================================================================

subroutine trcmix(name, ncol, clat, pmid, q)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Specify zonal mean mass mixing ratios of CH4, N2O, CFC11 and
! CFC12
! 
! Method: 
! Distributions assume constant mixing ratio in the troposphere
! and a decrease of mixing ratio in the stratosphere. Tropopause
! defined by ptrop. The scale height of the particular trace gas
! depends on latitude. This assumption produces a more realistic
! stratospheric distribution of the various trace gases.
! 
! Author: J. Kiehl
! 
!-----------------------------------------------------------------------

   use iso_c_binding, only: c_int64_t, c_double, c_loc

   ! Arguments
   character(len=*), intent(in)  :: name              ! constituent name
   integer,          intent(in)  :: ncol              ! number of columns
   real(r8), target, intent(in)  :: clat(pcols)       ! latitude in radians for columns
   real(r8), target, intent(in)  :: pmid(pcols,pver)  ! model pressures
   real(r8), target, intent(out) :: q(pcols,pver)     ! constituent mass mixing ratio

   integer i                ! longitude loop index
   integer k                ! level index

   real(r8) coslat(pcols)   ! cosine of latitude
   real(r8) dlat            ! latitude in degrees
   real(r8) ptrop           ! pressure level of tropopause
   real(r8) pratio          ! pressure divided by ptrop
   real(r8) trop_mmr        ! tropospheric mass mixing ratio
   real(r8) scale           ! pressure scale height
   real(r8) constant_mmr    ! spatially constant mass mixing ratio
   integer(c_int64_t) :: gas_id
   logical :: use_codon_scale
!-----------------------------------------------------------------------

   call ghg_data_trcmix_scale_select_impl()
   if (.not. use_native_ghg_data_trcmix_scale_impl) then
      gas_id = 0_c_int64_t
      trop_mmr = 0.0_r8
      constant_mmr = 0.0_r8

      if (name == 'O2') then
         gas_id = 5_c_int64_t
         constant_mmr = chem_surfvals_get('O2MMR')
      else if (name == 'CO2') then
         gas_id = 6_c_int64_t
         constant_mmr = chem_surfvals_co2_rad()
      else if (name == 'CH4') then
         gas_id = 1_c_int64_t
         trop_mmr = rmwch4 * chem_surfvals_get('CH4VMR')
      else if (name == 'N2O') then
         gas_id = 2_c_int64_t
         trop_mmr = rmwn2o * chem_surfvals_get('N2OVMR')
      else if (name == 'CFC11') then
         gas_id = 3_c_int64_t
         trop_mmr = rmwf11 * chem_surfvals_get('F11VMR')
      else if (name == 'CFC12') then
         gas_id = 4_c_int64_t
         trop_mmr = rmwf12 * chem_surfvals_get('F12VMR')
      end if

      if (gas_id /= 0_c_int64_t) then
         call ghg_data_trcmix_codon(gas_id, int(ncol, c_int64_t), int(pcols, c_int64_t), &
              int(pver, c_int64_t), real(trop_mmr, c_double), real(constant_mmr, c_double), &
              c_loc(clat(1)), c_loc(pmid(1,1)), c_loc(q(1,1)))
         call ghg_data_trcmix_direct_proof_once()
         return
      end if
   end if

   use_codon_scale = .false.

   do i = 1, ncol
      coslat(i) = cos(clat(i))
   end do

   if (name == 'O2') then

      q = chem_surfvals_get('O2MMR')

   else if (name == 'CO2') then

      q = chem_surfvals_co2_rad()

   else if (name == 'CH4') then

      ! set tropospheric mass mixing ratios
      trop_mmr = rmwch4 * chem_surfvals_get('CH4VMR')

      do k = 1,pver
         do i = 1,ncol
            ! set stratospheric scale height factor for gases
            dlat = abs(57.2958_r8 * clat(i))
            if (use_codon_scale) then
               call ghg_data_trcmix_scale_proof_once()
               scale = real(ghg_data_trcmix_scale_codon(1_c_int64_t, real(dlat, c_double)), r8)
            else
               if(dlat.le.45.0_r8) then
                  scale = 0.2353_r8
               else
                  scale = 0.2353_r8 + 0.0225489_r8 * (dlat - 45)
               end if
            end if

            ! pressure of tropopause
            ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

            ! determine output mass mixing ratios
            if (pmid(i,k) >= ptrop) then
               q(i,k) = trop_mmr
            else
               pratio = pmid(i,k)/ptrop
               q(i,k) = trop_mmr * (pratio)**scale
            end if
         end do
      end do

   else if (name == 'N2O') then

      ! set tropospheric mass mixing ratios
      trop_mmr = rmwn2o * chem_surfvals_get('N2OVMR')

      do k = 1,pver
         do i = 1,ncol
            ! set stratospheric scale height factor for gases
            dlat = abs(57.2958_r8 * clat(i))
            if (use_codon_scale) then
               call ghg_data_trcmix_scale_proof_once()
               scale = real(ghg_data_trcmix_scale_codon(2_c_int64_t, real(dlat, c_double)), r8)
            else
               if(dlat.le.45.0_r8) then
                  scale = 0.3478_r8 + 0.00116_r8 * dlat
               else
                  scale = 0.4000_r8 + 0.013333_r8 * (dlat - 45)
               end if
            end if

            ! pressure of tropopause
            ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

            ! determine output mass mixing ratios
            if (pmid(i,k) >= ptrop) then
               q(i,k) = trop_mmr
            else
               pratio = pmid(i,k)/ptrop
               q(i,k) = trop_mmr * (pratio)**scale
            end if
         end do
      end do

   else if (name == 'CFC11') then

      ! set tropospheric mass mixing ratios
      trop_mmr = rmwf11 * chem_surfvals_get('F11VMR')

      do k = 1,pver
         do i = 1,ncol
            ! set stratospheric scale height factor for gases
            dlat = abs(57.2958_r8 * clat(i))
            if (use_codon_scale) then
               call ghg_data_trcmix_scale_proof_once()
               scale = real(ghg_data_trcmix_scale_codon(3_c_int64_t, real(dlat, c_double)), r8)
            else
               if(dlat.le.45.0_r8) then
                  scale = 0.7273_r8 + 0.00606_r8 * dlat
               else
                  scale = 1.00_r8 + 0.013333_r8 * (dlat - 45)
               end if
            end if

            ! pressure of tropopause
            ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

            ! determine output mass mixing ratios
            if (pmid(i,k) >= ptrop) then
               q(i,k) = trop_mmr
            else
               pratio = pmid(i,k)/ptrop
               q(i,k) = trop_mmr * (pratio)**scale
            end if
         end do
      end do

   else if (name == 'CFC12') then

      ! set tropospheric mass mixing ratios
      trop_mmr = rmwf12 * chem_surfvals_get('F12VMR')

      do k = 1,pver
         do i = 1,ncol
            ! set stratospheric scale height factor for gases
            dlat = abs(57.2958_r8 * clat(i))
            if (use_codon_scale) then
               call ghg_data_trcmix_scale_proof_once()
               scale = real(ghg_data_trcmix_scale_codon(4_c_int64_t, real(dlat, c_double)), r8)
            else
               if(dlat.le.45.0_r8) then
                  scale = 0.4000_r8 + 0.00222_r8 * dlat
               else
                  scale = 0.50_r8 + 0.024444_r8 * (dlat - 45)
               end if
            end if

            ! pressure of tropopause
            ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

            ! determine output mass mixing ratios
            if (pmid(i,k) >= ptrop) then
               q(i,k) = trop_mmr
            else
               pratio = pmid(i,k)/ptrop
               q(i,k) = trop_mmr * (pratio)**scale
            end if
         end do
      end do

   end if

end subroutine trcmix

end module ghg_data
