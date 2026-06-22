!-----------------------------------------------------------------------
! Reads namelist options for gas-phase wet deposition
!
! Created by Francis Vitt -- 22 Apr 2011
!-----------------------------------------------------------------------
module gas_wetdep_opts

  use constituents,     only : pcnst
  use cam_logfile,      only : iulog
  use constituents,     only : pcnst
  use spmd_utils,       only : masterproc
  use cam_abortutils,   only : endrun
  use mo_util,          only : chemistry_misc_codon_touch
  use iso_c_binding,    only : c_int64_t, c_loc

  implicit none

  character(len=8), target :: gas_wetdep_list(pcnst) = ' '
  character(len=3) :: gas_wetdep_method = 'MOZ'
  integer :: gas_wetdep_cnt = 0

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  subroutine gas_wetdep_readnl(nlfile)

    use cam_abortutils,  only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
#ifdef SPMD
    use mpishorthand,    only: mpichar, mpicom
#endif

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    integer :: unitn, ierr
    integer(c_int64_t), target :: status_c(2)
    character(len=3), target :: gas_wetdep_method_c

    interface
       subroutine gas_wetdep_readnl_codon(pcnst_c, list_p, method_p, status_p) &
            bind(c, name="gas_wetdep_readnl_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: pcnst_c
         type(c_ptr), value :: list_p, method_p, status_p
       end subroutine gas_wetdep_readnl_codon
    end interface

    namelist /wetdep_inparm/ gas_wetdep_list
    namelist /wetdep_inparm/ gas_wetdep_method

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'wetdep_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, wetdep_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun('mo_neu_wetdep->wetdep_readnl: ERROR reading wetdep_inparm namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast (gas_wetdep_list, len(gas_wetdep_list(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast (gas_wetdep_method, len(gas_wetdep_method), mpichar, 0, mpicom)
#endif

    gas_wetdep_method_c = gas_wetdep_method
    status_c(:) = 0_c_int64_t
    call gas_wetdep_readnl_codon(int(pcnst, c_int64_t), c_loc(gas_wetdep_list), c_loc(gas_wetdep_method_c), &
         c_loc(status_c))
    call gas_wetdep_readnl_log_codon()

    gas_wetdep_cnt = int(status_c(1))

    if (status_c(2) /= 0_c_int64_t) then
       call endrun('gas_wetdep_readnl; gas_wetdep_method must be set to either MOZ or NEU')
    endif
    call chemistry_misc_codon_touch('gas_wetdep_readnl', 137)

  end subroutine gas_wetdep_readnl

  subroutine gas_wetdep_readnl_log_codon()

    implicit none

    if (masterproc) then
       write(iulog,*) 'gas_wetdep_readnl implementation = codon'
       call flush(iulog)
    end if

  end subroutine gas_wetdep_readnl_log_codon

end module gas_wetdep_opts
