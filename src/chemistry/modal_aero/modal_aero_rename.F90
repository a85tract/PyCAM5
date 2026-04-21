! modal_aero_rename.F90
!----------------------------------------------------------------------
!BOP
!
! !MODULE: modal_aero_rename --- modal aerosol mode merging (renaming)
!
! !INTERFACE:
  module modal_aero_rename

! !USES:
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use cam_abortutils,  only: endrun
  use cam_logfile,     only: iulog
  use mo_constants,    only: pi
  use chem_mods,       only: gas_pcnst
  use ppgrid,          only: pcols, pver
  use constituents,    only: pcnst, cnst_name
  use spmd_utils,      only: masterproc
  use modal_aero_data, only: maxd_aspectype, ntot_amode
  use modal_aero_data, only: alnsg_amode, voltonumblo_amode, voltonumbhi_amode, dgnum_amode, nspec_amode
  use modal_aero_data, only: lspectype_amode, specmw_amode, specdens_amode, lmassptr_amode, lmassptrcw_amode
  use modal_aero_data, only: numptr_amode, numptrcw_amode, modeptr_coarse, modeptr_accum, lspectype_amode
  use modal_aero_data, only: specmw_amode, specdens_amode, lmassptr_amode, lmassptrcw_amode, numptr_amode, numptrcw_amode
  use modal_aero_data, only: dgnumhi_amode, dgnumlo_amode, cnst_name_cw, modeptr_aitken

  implicit none
  private
  save

! !PUBLIC MEMBER FUNCTIONS:
  public modal_aero_rename_sub, modal_aero_rename_init

! !PUBLIC DATA MEMBERS:
  integer, parameter :: pcnstxx = gas_pcnst
  integer, parameter, public :: maxspec_renamexf = maxd_aspectype

! *** select one of the 3 following options
! *** for maxpair_renamexf = 2 or 3, use mode definition files with
!     dgnumhi_amode(modeptr_accum)  = 1.1e-6 m
!     dgnumlo_amode(modeptr_coarse) = 0.9e-6 m

! integer, parameter, public :: maxpair_renamexf = 1
! integer, parameter, public :: ipair_select_renamexf(maxpair_renamexf) = (/ 2001 /)

! integer, parameter, public :: maxpair_renamexf = 2
! integer, parameter, public :: ipair_select_renamexf(maxpair_renamexf) = (/ 2001, 1003 /)

  integer, parameter, public :: maxpair_renamexf = 3
  integer, parameter, public :: ipair_select_renamexf(maxpair_renamexf) = (/ 2001, 1003, 3001 /)
! ipair_select_renamexf defines the mode_from and mode_too for each renaming pair
! 2001 = aitken --> accum
! 1003 = accum  --> coarse
! 3001 = coarse --> accum

  integer, parameter, public :: method_optbb_renamexf = 2

  integer, public :: npair_renamexf = -123456789
  integer, public :: modefrm_renamexf(maxpair_renamexf)
  integer, public :: modetoo_renamexf(maxpair_renamexf)
  integer, public :: nspecfrm_renamexf(maxpair_renamexf)

  integer, public :: lspecfrma_renamexf(maxspec_renamexf,maxpair_renamexf)
  integer, public :: lspecfrmc_renamexf(maxspec_renamexf,maxpair_renamexf)
  integer, public :: lspectooa_renamexf(maxspec_renamexf,maxpair_renamexf)
  integer, public :: lspectooc_renamexf(maxspec_renamexf,maxpair_renamexf)

  integer, public :: igrow_shrink_renamexf(maxpair_renamexf)
  integer, public :: ixferable_all_renamexf(maxpair_renamexf)
  integer, public :: ixferable_all_needed_renamexf(maxpair_renamexf)
  integer, public :: ixferable_a_renamexf(maxspec_renamexf,maxpair_renamexf)
  integer, public :: ixferable_c_renamexf(maxspec_renamexf,maxpair_renamexf)

! !PRIVATE DATA MEMBERS:
  integer :: ido_mode_calcaa(ntot_amode)
  real (r8) :: dp_belowcut(maxpair_renamexf)
  real (r8) :: dp_cut(maxpair_renamexf)
  real (r8) :: dp_xferall_thresh(maxpair_renamexf)
  real (r8) :: dp_xfernone_threshaa(maxpair_renamexf)
  real (r8) :: dryvol_smallest(ntot_amode)
  real (r8) :: factoraa(ntot_amode)
  real (r8) :: factoryy(ntot_amode)
  real (r8) :: lndp_cut(maxpair_renamexf)
  real (r8) :: factor_3alnsg2(maxpair_renamexf)
  real (r8) :: v2nhirlx(ntot_amode), v2nlorlx(ntot_amode)

  logical :: modal_accum_coarse_exch = .false.
  logical :: modal_aero_rename_no_acc_crs_dryvols_use_native_impl = .false.
  logical :: modal_aero_rename_no_acc_crs_dryvols_impl_selected = .false.
  logical :: modal_aero_rename_no_acc_crs_xferfracs_use_native_impl = .false.
  logical :: modal_aero_rename_no_acc_crs_xferfracs_impl_selected = .false.
  logical :: modal_aero_rename_no_acc_crs_tendencies_use_native_impl = .false.
  logical :: modal_aero_rename_no_acc_crs_tendencies_impl_selected = .false.
  logical :: modal_aero_rename_acc_crs_dryvols_use_native_impl = .false.
  logical :: modal_aero_rename_acc_crs_dryvols_impl_selected = .false.
  logical :: modal_aero_rename_acc_crs_xferfracs_use_native_impl = .false.
  logical :: modal_aero_rename_acc_crs_xferfracs_impl_selected = .false.
  logical :: modal_aero_rename_acc_crs_tendencies_use_native_impl = .false.
  logical :: modal_aero_rename_acc_crs_tendencies_impl_selected = .false.
  logical :: modal_aero_rename_acc_crs_pair_use_native_impl = .false.
  logical :: modal_aero_rename_acc_crs_pair_impl_selected = .false.
  logical :: modal_aero_rename_acc_crs_sub_use_native_impl = .false.
  logical :: modal_aero_rename_acc_crs_sub_impl_selected = .false.
  logical :: modal_aero_rename_set_dotend_flags_use_native_impl = .false.
  logical :: modal_aero_rename_set_dotend_flags_impl_selected = .false.

! !DESCRIPTION: This module implements ...
!
! !REVISION HISTORY:
!
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! list private module data here

!EOC
!----------------------------------------------------------------------
contains

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_init(modal_accum_coarse_exch_in)
    logical, optional, intent(in) :: modal_accum_coarse_exch_in
    
    if (present(modal_accum_coarse_exch_in)) then
       modal_accum_coarse_exch = modal_accum_coarse_exch_in
    endif

    if (modal_accum_coarse_exch) then
       call modal_aero_rename_acc_crs_init()
    else
       call modal_aero_rename_no_acc_crs_init()
    endif

  end subroutine modal_aero_rename_init

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_no_acc_crs_dryvols_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_rename_no_acc_crs_dryvols_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('MODAL_AERO_RENAME_NO_ACC_CRS_DRYVOLS_IMPL', &
         value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       modal_aero_rename_no_acc_crs_dryvols_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       modal_aero_rename_no_acc_crs_dryvols_use_native_impl = .false.
    end if

    modal_aero_rename_no_acc_crs_dryvols_impl_selected = .true.

    if (masterproc) then
       if (modal_aero_rename_no_acc_crs_dryvols_use_native_impl) then
          write(iulog,*) 'modal_aero_rename_no_acc_crs_dryvols implementation = native'
       else
          write(iulog,*) 'modal_aero_rename_no_acc_crs_dryvols implementation = codon'
       end if
    end if

  end subroutine modal_aero_rename_no_acc_crs_dryvols_select_impl

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_no_acc_crs_xferfracs_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_rename_no_acc_crs_xferfracs_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('MODAL_AERO_RENAME_NO_ACC_CRS_XFERFRACS_IMPL', &
         value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       modal_aero_rename_no_acc_crs_xferfracs_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       modal_aero_rename_no_acc_crs_xferfracs_use_native_impl = .false.
    end if

    modal_aero_rename_no_acc_crs_xferfracs_impl_selected = .true.

    if (masterproc) then
       if (modal_aero_rename_no_acc_crs_xferfracs_use_native_impl) then
          write(iulog,*) 'modal_aero_rename_no_acc_crs_xferfracs implementation = native'
       else
          write(iulog,*) 'modal_aero_rename_no_acc_crs_xferfracs implementation = codon'
       end if
    end if

  end subroutine modal_aero_rename_no_acc_crs_xferfracs_select_impl

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_no_acc_crs_tendencies_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_rename_no_acc_crs_tendencies_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('MODAL_AERO_RENAME_NO_ACC_CRS_TENDENCIES_IMPL', &
         value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       modal_aero_rename_no_acc_crs_tendencies_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       modal_aero_rename_no_acc_crs_tendencies_use_native_impl = .false.
    end if

    modal_aero_rename_no_acc_crs_tendencies_impl_selected = .true.

    if (masterproc) then
       if (modal_aero_rename_no_acc_crs_tendencies_use_native_impl) then
          write(iulog,*) 'modal_aero_rename_no_acc_crs_tendencies implementation = native'
       else
          write(iulog,*) 'modal_aero_rename_no_acc_crs_tendencies implementation = codon'
       end if
    end if

  end subroutine modal_aero_rename_no_acc_crs_tendencies_select_impl

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_dryvols_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_rename_acc_crs_dryvols_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('MODAL_AERO_RENAME_ACC_CRS_DRYVOLS_IMPL', &
         value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       modal_aero_rename_acc_crs_dryvols_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       modal_aero_rename_acc_crs_dryvols_use_native_impl = .false.
    end if

    modal_aero_rename_acc_crs_dryvols_impl_selected = .true.

    if (masterproc) then
       if (modal_aero_rename_acc_crs_dryvols_use_native_impl) then
          write(iulog,*) 'modal_aero_rename_acc_crs_dryvols implementation = native'
       else
          write(iulog,*) 'modal_aero_rename_acc_crs_dryvols implementation = codon'
       end if
    end if

  end subroutine modal_aero_rename_acc_crs_dryvols_select_impl

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_xferfracs_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_rename_acc_crs_xferfracs_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('MODAL_AERO_RENAME_ACC_CRS_XFERFRACS_IMPL', &
         value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       modal_aero_rename_acc_crs_xferfracs_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       modal_aero_rename_acc_crs_xferfracs_use_native_impl = .false.
    end if

    modal_aero_rename_acc_crs_xferfracs_impl_selected = .true.

    if (masterproc) then
       if (modal_aero_rename_acc_crs_xferfracs_use_native_impl) then
          write(iulog,*) 'modal_aero_rename_acc_crs_xferfracs implementation = native'
       else
          write(iulog,*) 'modal_aero_rename_acc_crs_xferfracs implementation = codon'
       end if
    end if

  end subroutine modal_aero_rename_acc_crs_xferfracs_select_impl

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_tendencies_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_rename_acc_crs_tendencies_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('MODAL_AERO_RENAME_ACC_CRS_TENDENCIES_IMPL', &
         value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       modal_aero_rename_acc_crs_tendencies_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       modal_aero_rename_acc_crs_tendencies_use_native_impl = .false.
    end if

    modal_aero_rename_acc_crs_tendencies_impl_selected = .true.

    if (masterproc) then
       if (modal_aero_rename_acc_crs_tendencies_use_native_impl) then
          write(iulog,*) 'modal_aero_rename_acc_crs_tendencies implementation = native'
       else
          write(iulog,*) 'modal_aero_rename_acc_crs_tendencies implementation = codon'
       end if
    end if

  end subroutine modal_aero_rename_acc_crs_tendencies_select_impl

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_set_dotend_flags_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_rename_set_dotend_flags_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('MODAL_AERO_RENAME_SET_DOTEND_FLAGS_IMPL', &
         value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       modal_aero_rename_set_dotend_flags_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       modal_aero_rename_set_dotend_flags_use_native_impl = .false.
    end if

    modal_aero_rename_set_dotend_flags_impl_selected = .true.

    if (masterproc) then
       if (modal_aero_rename_set_dotend_flags_use_native_impl) then
          write(iulog,*) 'modal_aero_rename_set_dotend_flags implementation = native'
       else
          write(iulog,*) 'modal_aero_rename_set_dotend_flags implementation = codon'
       end if
    end if

  end subroutine modal_aero_rename_set_dotend_flags_select_impl

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_pair_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_rename_acc_crs_pair_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('MODAL_AERO_RENAME_ACC_CRS_PAIR_IMPL', &
         value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       modal_aero_rename_acc_crs_pair_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       modal_aero_rename_acc_crs_pair_use_native_impl = .false.
    end if

    modal_aero_rename_acc_crs_pair_impl_selected = .true.

    if (masterproc) then
       if (modal_aero_rename_acc_crs_pair_use_native_impl) then
          write(iulog,*) 'modal_aero_rename_acc_crs_pair implementation = native'
       else
          write(iulog,*) 'modal_aero_rename_acc_crs_pair implementation = codon'
       end if
    end if

  end subroutine modal_aero_rename_acc_crs_pair_select_impl

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_sub_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (modal_aero_rename_acc_crs_sub_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('MODAL_AERO_RENAME_ACC_CRS_SUB_IMPL', &
         value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       modal_aero_rename_acc_crs_sub_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       modal_aero_rename_acc_crs_sub_use_native_impl = .false.
    end if

    modal_aero_rename_acc_crs_sub_impl_selected = .true.

    if (masterproc) then
       if (modal_aero_rename_acc_crs_sub_use_native_impl) then
          write(iulog,*) 'modal_aero_rename_acc_crs_sub implementation = native'
       else
          write(iulog,*) 'modal_aero_rename_acc_crs_sub implementation = codon'
       end if
    end if

  end subroutine modal_aero_rename_acc_crs_sub_select_impl

  !------------------------------------------------------------------
  !------------------------------------------------------------------
  subroutine modal_aero_rename_sub(                       &
       fromwhere,         lchnk,               &
       ncol,              nstep,               &
       loffset,           deltat,              &
       pdel,              troplev,             &
       dotendrn,          q,                   &
       dqdt,              dqdt_other,          &
       dotendqqcwrn,      qqcw,                &
       dqqcwdt,           dqqcwdt_other,       &
       is_dorename_atik,  dorename_atik,       &
       jsrflx_rename,     nsrflx,              &
       qsrflx,            qqcwsrflx,           &
       dqdt_rnpos                              )


    ! !PARAMETERS:
    character(len=*), intent(in) :: fromwhere    ! identifies which module
    ! is making the call
    integer,  intent(in)    :: lchnk                ! chunk identifier
    integer,  intent(in)    :: ncol                 ! number of atmospheric column
    integer,  intent(in)    :: nstep                ! model time-step number
    integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
    real(r8), intent(in)    :: deltat               ! time step (s)
    integer,  intent(in)    :: troplev(pcols)

    real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
    real(r8), intent(in)    :: q(ncol,pver,pcnstxx) ! tracer mixing ratio array
    ! *** MUST BE mol/mol-air or #/mol-air
    ! *** NOTE ncol and pcnstxx dimensions
    real(r8), intent(in)    :: qqcw(ncol,pver,pcnstxx) ! like q but for cloud-borne species

    real(r8), intent(inout) :: dqdt(ncol,pver,pcnstxx)  ! TMR tendency array;
    ! incoming dqdt = tendencies for the 
    !     "fromwhere" continuous growth process 
    ! the renaming tendencies are added on
    ! *** NOTE ncol and pcnstxx dimensions
    real(r8), intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), intent(in)    :: dqdt_other(ncol,pver,pcnstxx)  
    ! tendencies for "other" continuous growth process 
    ! currently in cam3
    !     dqdt is from gas (h2so4, nh3) condensation
    !     dqdt_other is from aqchem and soa
    ! *** NOTE ncol and pcnstxx dimensions
    real(r8), intent(in)    :: dqqcwdt_other(ncol,pver,pcnstxx)  
    logical,  intent(inout) :: dotendrn(pcnstxx) ! identifies the species for which
    !     renaming dqdt is computed
    logical,  intent(inout) :: dotendqqcwrn(pcnstxx)

    logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
    logical,  intent(in)    :: dorename_atik(ncol,pver) ! true if renaming should
    ! be done at i,k
    integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
    integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

    real(r8), intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
    ! process-specific column tracer tendencies 
    real(r8), intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
    real(r8), optional, intent(out) &
         :: dqdt_rnpos(ncol,pver,pcnstxx)
    ! the positive (production) part of the renaming tendency

    if (modal_accum_coarse_exch) then
       call modal_aero_rename_acc_crs_sub(        &
            fromwhere,         lchnk,               &
            ncol,              nstep,               &
            loffset,           deltat,              &
            pdel,              troplev,             &
            dotendrn,          q,                   &
            dqdt,              dqdt_other,          &
            dotendqqcwrn,      qqcw,                &
            dqqcwdt,           dqqcwdt_other,       &
            is_dorename_atik,  dorename_atik,       &
            jsrflx_rename,     nsrflx,              &
            qsrflx,            qqcwsrflx,           &
            dqdt_rnpos                              )
    else
       call modal_aero_rename_no_acc_crs_sub(             &
            fromwhere,         lchnk,               &
            ncol,              nstep,               &
            loffset,           deltat,              &
            pdel,                                   &
            dotendrn,          q,                   &
            dqdt,              dqdt_other,          &
            dotendqqcwrn,      qqcw,                &
            dqqcwdt,           dqqcwdt_other,       &
            is_dorename_atik,  dorename_atik,       &
            jsrflx_rename,     nsrflx,              &
            qsrflx,            qqcwsrflx            )
    endif
  end subroutine modal_aero_rename_sub

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! private methods
!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_no_acc_crs_dryvols( ncol, loffset, deltat, &
       idomode, q, qqcw, dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, &
       nspec_amode_in, lspectype_amode_in, specmw_amode_in, specdens_amode_in, &
       lmassptr_amode_in, lmassptrcw_amode_in, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c )

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    real(r8), intent(in) :: deltat
    integer, target, intent(in) :: idomode(ntot_amode)
    real(r8), target, intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqdt(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqdt_other(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqqcwdt_other(ncol,pver,pcnstxx)
    integer, target, intent(in) :: nspec_amode_in(ntot_amode)
    integer, target, intent(in) :: lspectype_amode_in(maxspec_renamexf,ntot_amode)
    real(r8), target, intent(in) :: specmw_amode_in(*)
    real(r8), target, intent(in) :: specdens_amode_in(*)
    integer, target, intent(in) :: lmassptr_amode_in(maxspec_renamexf,ntot_amode)
    integer, target, intent(in) :: lmassptrcw_amode_in(maxspec_renamexf,ntot_amode)
    real(r8), target, intent(out) :: dryvol_a(ncol,pver,ntot_amode)
    real(r8), target, intent(out) :: dryvol_c(ncol,pver,ntot_amode)
    real(r8), target, intent(out) :: deldryvol_a(ncol,pver,ntot_amode)
    real(r8), target, intent(out) :: deldryvol_c(ncol,pver,ntot_amode)
    integer(c_int64_t), target :: idomode_c(ntot_amode)
    integer(c_int64_t), target :: nspec_amode_c(ntot_amode)
    integer(c_int64_t), target :: lspectype_amode_c(maxspec_renamexf,ntot_amode)
    integer(c_int64_t), target :: lmassptr_amode_c(maxspec_renamexf,ntot_amode)
    integer(c_int64_t), target :: lmassptrcw_amode_c(maxspec_renamexf,ntot_amode)
    integer :: n, l1

    interface
       subroutine modal_aero_rename_no_acc_crs_dryvols_codon( &
            ncol_c, pver_c, pcnstxx_c, ntot_amode_c, maxspec_renamexf_c, loffset_c, deltat_c, &
            idomode_p, q_p, qqcw_p, dqdt_p, dqdt_other_p, dqqcwdt_p, dqqcwdt_other_p, &
            nspec_amode_p, lspectype_amode_p, specmw_amode_p, specdens_amode_p, &
            lmassptr_amode_p, lmassptrcw_amode_p, dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p ) &
            bind(c, name="modal_aero_rename_no_acc_crs_dryvols_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, pcnstxx_c, ntot_amode_c, maxspec_renamexf_c, loffset_c
         real(c_double), value :: deltat_c
         type(c_ptr), value :: idomode_p, q_p, qqcw_p, dqdt_p, dqdt_other_p, dqqcwdt_p, dqqcwdt_other_p
         type(c_ptr), value :: nspec_amode_p, lspectype_amode_p, specmw_amode_p, specdens_amode_p
         type(c_ptr), value :: lmassptr_amode_p, lmassptrcw_amode_p
         type(c_ptr), value :: dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p
       end subroutine modal_aero_rename_no_acc_crs_dryvols_codon
    end interface

    call modal_aero_rename_no_acc_crs_dryvols_select_impl()

    if (.not. modal_aero_rename_no_acc_crs_dryvols_use_native_impl) then
       do n = 1, ntot_amode
          idomode_c(n) = int(idomode(n), c_int64_t)
          nspec_amode_c(n) = int(nspec_amode_in(n), c_int64_t)
          do l1 = 1, maxspec_renamexf
             lspectype_amode_c(l1,n) = int(lspectype_amode_in(l1,n), c_int64_t)
             lmassptr_amode_c(l1,n) = int(lmassptr_amode_in(l1,n), c_int64_t)
             lmassptrcw_amode_c(l1,n) = int(lmassptrcw_amode_in(l1,n), c_int64_t)
          end do
       end do

       call modal_aero_rename_no_acc_crs_dryvols_codon( &
            int(ncol, c_int64_t), int(pver, c_int64_t), int(pcnstxx, c_int64_t), &
            int(ntot_amode, c_int64_t), int(maxspec_renamexf, c_int64_t), int(loffset, c_int64_t), &
            real(deltat, c_double), c_loc(idomode_c(1)), c_loc(q(1,1,1)), c_loc(qqcw(1,1,1)), &
            c_loc(dqdt(1,1,1)), c_loc(dqdt_other(1,1,1)), c_loc(dqqcwdt(1,1,1)), c_loc(dqqcwdt_other(1,1,1)), &
            c_loc(nspec_amode_c(1)), c_loc(lspectype_amode_c(1,1)), c_loc(specmw_amode_in(1)), c_loc(specdens_amode_in(1)), &
            c_loc(lmassptr_amode_c(1,1)), c_loc(lmassptrcw_amode_c(1,1)), c_loc(dryvol_a(1,1,1)), c_loc(dryvol_c(1,1,1)), &
            c_loc(deldryvol_a(1,1,1)), c_loc(deldryvol_c(1,1,1)) &
       )
       return
    end if

    call modal_aero_rename_no_acc_crs_dryvols_native( ncol, loffset, deltat, &
         idomode, q, qqcw, dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, &
         nspec_amode_in, lspectype_amode_in, specmw_amode_in, specdens_amode_in, &
         lmassptr_amode_in, lmassptrcw_amode_in, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c )

  end subroutine modal_aero_rename_no_acc_crs_dryvols

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_no_acc_crs_dryvols_native( ncol, loffset, deltat, &
       idomode, q, qqcw, dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, &
       nspec_amode_in, lspectype_amode_in, specmw_amode_in, specdens_amode_in, &
       lmassptr_amode_in, lmassptrcw_amode_in, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c )

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    real(r8), intent(in) :: deltat
    integer, intent(in) :: idomode(ntot_amode)
    real(r8), intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqdt(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqdt_other(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqqcwdt_other(ncol,pver,pcnstxx)
    integer, intent(in) :: nspec_amode_in(ntot_amode)
    integer, intent(in) :: lspectype_amode_in(maxspec_renamexf,ntot_amode)
    real(r8), intent(in) :: specmw_amode_in(*)
    real(r8), intent(in) :: specdens_amode_in(*)
    integer, intent(in) :: lmassptr_amode_in(maxspec_renamexf,ntot_amode)
    integer, intent(in) :: lmassptrcw_amode_in(maxspec_renamexf,ntot_amode)
    real(r8), intent(out) :: dryvol_a(ncol,pver,ntot_amode)
    real(r8), intent(out) :: dryvol_c(ncol,pver,ntot_amode)
    real(r8), intent(out) :: deldryvol_a(ncol,pver,ntot_amode)
    real(r8), intent(out) :: deldryvol_c(ncol,pver,ntot_amode)

    integer :: n, l1, l2, la, lc
    real(r8) :: dum_m2v, dum_m2vdt

    do n = 1, ntot_amode
       if (idomode(n) .gt. 0) then
          dryvol_a(1:ncol,:,n) = 0.0_r8
          dryvol_c(1:ncol,:,n) = 0.0_r8
          deldryvol_a(1:ncol,:,n) = 0.0_r8
          deldryvol_c(1:ncol,:,n) = 0.0_r8
          do l1 = 1, nspec_amode_in(n)
             l2 = lspectype_amode_in(l1,n)
             dum_m2v = specmw_amode_in(l2) / specdens_amode_in(l2)
             dum_m2vdt = dum_m2v*deltat
             la = lmassptr_amode_in(l1,n)-loffset
             if (la > 0) then
                dryvol_a(1:ncol,:,n) = dryvol_a(1:ncol,:,n)    &
                   + dum_m2v*max( 0.0_r8,   &
                     q(1:ncol,:,la)-deltat*dqdt_other(1:ncol,:,la) )
                deldryvol_a(1:ncol,:,n) = deldryvol_a(1:ncol,:,n)    &
                   + (dqdt_other(1:ncol,:,la) + dqdt(1:ncol,:,la))*dum_m2vdt
             end if

             lc = lmassptrcw_amode_in(l1,n)-loffset
             if (lc > 0) then
                dryvol_c(1:ncol,:,n) = dryvol_c(1:ncol,:,n)    &
                   + dum_m2v*max( 0.0_r8,   &
                     qqcw(1:ncol,:,lc)-deltat*dqqcwdt_other(1:ncol,:,lc) )
                deldryvol_c(1:ncol,:,n) = deldryvol_c(1:ncol,:,n)    &
                   + (dqqcwdt_other(1:ncol,:,lc) + dqqcwdt(1:ncol,:,lc))*dum_m2vdt
             end if
          end do
       end if
    end do

  end subroutine modal_aero_rename_no_acc_crs_dryvols_native

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_no_acc_crs_xferfracs( ncol, loffset, npair_renamexf_in, &
       q, qqcw, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, &
       modefrm_renamexf_in, modetoo_renamexf_in, numptr_amode_in, numptrcw_amode_in, &
       dgnum_amode_in, factoraa_in, factoryy_in, dryvol_smallest_in, v2nlorlx_in, v2nhirlx_in, &
       dum3alnsg2_in, dp_cut_in, lndp_cut_in, dp_belowcut_in, onethird_in, xferfrac_max_in, &
       xferfrac_vol_out, xferfrac_num_out )

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    integer, intent(in) :: npair_renamexf_in
    real(r8), target, intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dryvol_a(ncol,pver,ntot_amode)
    real(r8), target, intent(in) :: dryvol_c(ncol,pver,ntot_amode)
    real(r8), target, intent(in) :: deldryvol_a(ncol,pver,ntot_amode)
    real(r8), target, intent(in) :: deldryvol_c(ncol,pver,ntot_amode)
    integer, target, intent(in) :: modefrm_renamexf_in(maxpair_renamexf)
    integer, target, intent(in) :: modetoo_renamexf_in(maxpair_renamexf)
    integer, target, intent(in) :: numptr_amode_in(ntot_amode)
    integer, target, intent(in) :: numptrcw_amode_in(ntot_amode)
    real(r8), target, intent(in) :: dgnum_amode_in(ntot_amode)
    real(r8), target, intent(in) :: factoraa_in(ntot_amode)
    real(r8), target, intent(in) :: factoryy_in(ntot_amode)
    real(r8), target, intent(in) :: dryvol_smallest_in(ntot_amode)
    real(r8), target, intent(in) :: v2nlorlx_in(ntot_amode)
    real(r8), target, intent(in) :: v2nhirlx_in(ntot_amode)
    real(r8), target, intent(in) :: dum3alnsg2_in(maxpair_renamexf)
    real(r8), target, intent(in) :: dp_cut_in(maxpair_renamexf)
    real(r8), target, intent(in) :: lndp_cut_in(maxpair_renamexf)
    real(r8), target, intent(in) :: dp_belowcut_in(maxpair_renamexf)
    real(r8), intent(in) :: onethird_in
    real(r8), intent(in) :: xferfrac_max_in
    real(r8), target, intent(out) :: xferfrac_vol_out(ncol,pver,maxpair_renamexf)
    real(r8), target, intent(out) :: xferfrac_num_out(ncol,pver,maxpair_renamexf)
    integer(c_int64_t), target :: modefrm_renamexf_c(maxpair_renamexf)
    integer(c_int64_t), target :: modetoo_renamexf_c(maxpair_renamexf)
    integer(c_int64_t), target :: numptr_amode_c(ntot_amode)
    integer(c_int64_t), target :: numptrcw_amode_c(ntot_amode)
    integer :: n, ipair

    interface
       subroutine modal_aero_rename_no_acc_crs_xferfracs_codon( &
            ncol_c, pver_c, pcnstxx_c, ntot_amode_c, maxpair_renamexf_c, loffset_c, npair_renamexf_c, &
            q_p, qqcw_p, dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p, &
            modefrm_renamexf_p, modetoo_renamexf_p, numptr_amode_p, numptrcw_amode_p, &
            dgnum_amode_p, factoraa_p, factoryy_p, dryvol_smallest_p, v2nlorlx_p, v2nhirlx_p, &
            dum3alnsg2_p, dp_cut_p, lndp_cut_p, dp_belowcut_p, onethird_c, xferfrac_max_c, &
            xferfrac_vol_p, xferfrac_num_p ) bind(c, name="modal_aero_rename_no_acc_crs_xferfracs_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, pcnstxx_c, ntot_amode_c, maxpair_renamexf_c
         integer(c_int64_t), value :: loffset_c, npair_renamexf_c
         real(c_double), value :: onethird_c, xferfrac_max_c
         type(c_ptr), value :: q_p, qqcw_p, dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p
         type(c_ptr), value :: modefrm_renamexf_p, modetoo_renamexf_p, numptr_amode_p, numptrcw_amode_p
         type(c_ptr), value :: dgnum_amode_p, factoraa_p, factoryy_p, dryvol_smallest_p, v2nlorlx_p, v2nhirlx_p
         type(c_ptr), value :: dum3alnsg2_p, dp_cut_p, lndp_cut_p, dp_belowcut_p
         type(c_ptr), value :: xferfrac_vol_p, xferfrac_num_p
       end subroutine modal_aero_rename_no_acc_crs_xferfracs_codon
    end interface

    call modal_aero_rename_no_acc_crs_xferfracs_select_impl()

    if (.not. modal_aero_rename_no_acc_crs_xferfracs_use_native_impl) then
       do ipair = 1, maxpair_renamexf
          modefrm_renamexf_c(ipair) = int(modefrm_renamexf_in(ipair), c_int64_t)
          modetoo_renamexf_c(ipair) = int(modetoo_renamexf_in(ipair), c_int64_t)
       end do
       do n = 1, ntot_amode
          numptr_amode_c(n) = int(numptr_amode_in(n), c_int64_t)
          numptrcw_amode_c(n) = int(numptrcw_amode_in(n), c_int64_t)
       end do

       call modal_aero_rename_no_acc_crs_xferfracs_codon( &
            int(ncol, c_int64_t), int(pver, c_int64_t), int(pcnstxx, c_int64_t), int(ntot_amode, c_int64_t), &
            int(maxpair_renamexf, c_int64_t), int(loffset, c_int64_t), int(npair_renamexf_in, c_int64_t), &
            c_loc(q(1,1,1)), c_loc(qqcw(1,1,1)), c_loc(dryvol_a(1,1,1)), c_loc(dryvol_c(1,1,1)), &
            c_loc(deldryvol_a(1,1,1)), c_loc(deldryvol_c(1,1,1)), c_loc(modefrm_renamexf_c(1)), c_loc(modetoo_renamexf_c(1)), &
            c_loc(numptr_amode_c(1)), c_loc(numptrcw_amode_c(1)), c_loc(dgnum_amode_in(1)), c_loc(factoraa_in(1)), &
            c_loc(factoryy_in(1)), c_loc(dryvol_smallest_in(1)), c_loc(v2nlorlx_in(1)), c_loc(v2nhirlx_in(1)), &
            c_loc(dum3alnsg2_in(1)), c_loc(dp_cut_in(1)), c_loc(lndp_cut_in(1)), c_loc(dp_belowcut_in(1)), &
            real(onethird_in, c_double), real(xferfrac_max_in, c_double), c_loc(xferfrac_vol_out(1,1,1)), &
            c_loc(xferfrac_num_out(1,1,1)) &
       )
       return
    end if

    call modal_aero_rename_no_acc_crs_xferfracs_native( ncol, loffset, npair_renamexf_in, &
         q, qqcw, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, modefrm_renamexf_in, modetoo_renamexf_in, &
         numptr_amode_in, numptrcw_amode_in, dgnum_amode_in, factoraa_in, factoryy_in, dryvol_smallest_in, &
         v2nlorlx_in, v2nhirlx_in, dum3alnsg2_in, dp_cut_in, lndp_cut_in, dp_belowcut_in, onethird_in, &
         xferfrac_max_in, xferfrac_vol_out, xferfrac_num_out )

  end subroutine modal_aero_rename_no_acc_crs_xferfracs

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_no_acc_crs_xferfracs_native( ncol, loffset, npair_renamexf_in, &
       q, qqcw, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, &
       modefrm_renamexf_in, modetoo_renamexf_in, numptr_amode_in, numptrcw_amode_in, &
       dgnum_amode_in, factoraa_in, factoryy_in, dryvol_smallest_in, v2nlorlx_in, v2nhirlx_in, &
       dum3alnsg2_in, dp_cut_in, lndp_cut_in, dp_belowcut_in, onethird_in, xferfrac_max_in, &
       xferfrac_vol_out, xferfrac_num_out )

    use shr_spfn_mod, only: erfc => shr_spfn_erfc

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    integer, intent(in) :: npair_renamexf_in
    real(r8), intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dryvol_a(ncol,pver,ntot_amode)
    real(r8), intent(in) :: dryvol_c(ncol,pver,ntot_amode)
    real(r8), intent(in) :: deldryvol_a(ncol,pver,ntot_amode)
    real(r8), intent(in) :: deldryvol_c(ncol,pver,ntot_amode)
    integer, intent(in) :: modefrm_renamexf_in(maxpair_renamexf)
    integer, intent(in) :: modetoo_renamexf_in(maxpair_renamexf)
    integer, intent(in) :: numptr_amode_in(ntot_amode)
    integer, intent(in) :: numptrcw_amode_in(ntot_amode)
    real(r8), intent(in) :: dgnum_amode_in(ntot_amode)
    real(r8), intent(in) :: factoraa_in(ntot_amode)
    real(r8), intent(in) :: factoryy_in(ntot_amode)
    real(r8), intent(in) :: dryvol_smallest_in(ntot_amode)
    real(r8), intent(in) :: v2nlorlx_in(ntot_amode)
    real(r8), intent(in) :: v2nhirlx_in(ntot_amode)
    real(r8), intent(in) :: dum3alnsg2_in(maxpair_renamexf)
    real(r8), intent(in) :: dp_cut_in(maxpair_renamexf)
    real(r8), intent(in) :: lndp_cut_in(maxpair_renamexf)
    real(r8), intent(in) :: dp_belowcut_in(maxpair_renamexf)
    real(r8), intent(in) :: onethird_in
    real(r8), intent(in) :: xferfrac_max_in
    real(r8), intent(out) :: xferfrac_vol_out(ncol,pver,maxpair_renamexf)
    real(r8), intent(out) :: xferfrac_num_out(ncol,pver,maxpair_renamexf)

    integer :: i, k, ipair, mfrm, mtoo
    real(r8) :: dgn_t_new, dgn_t_old
    real(r8) :: dryvol_t_del, dryvol_t_new
    real(r8) :: dryvol_t_old, dryvol_t_oldbnd
    real(r8) :: dum
    real(r8) :: lndgn_new, lndgn_old
    real(r8) :: lndgv_new, lndgv_old
    real(r8) :: num_t_old, num_t_oldbnd
    real(r8) :: tailfr_volnew, tailfr_volold
    real(r8) :: tailfr_numnew, tailfr_numold
    real(r8) :: yn_tail, yv_tail

    xferfrac_vol_out(:,:,:) = 0.0_r8
    xferfrac_num_out(:,:,:) = 0.0_r8

    do ipair = 1, npair_renamexf_in
       mfrm = modefrm_renamexf_in(ipair)
       mtoo = modetoo_renamexf_in(ipair)

       do k = 1, pver
       do i = 1, ncol
          dryvol_t_old = dryvol_a(i,k,mfrm) + dryvol_c(i,k,mfrm)
          dryvol_t_del = deldryvol_a(i,k,mfrm) + deldryvol_c(i,k,mfrm)
          dryvol_t_new = dryvol_t_old + dryvol_t_del
          dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest_in(mfrm) )

          if (dryvol_t_new .le. dryvol_smallest_in(mfrm)) cycle
          if (dryvol_t_del .le. 1.0e-6_r8*dryvol_t_oldbnd) cycle

          num_t_old = q(i,k,numptr_amode_in(mfrm)-loffset)
          num_t_old = num_t_old + qqcw(i,k,numptrcw_amode_in(mfrm)-loffset)
          num_t_old = max( 0.0_r8, num_t_old )
          dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest_in(mfrm) )
          num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx_in(mfrm), num_t_old )
          num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx_in(mfrm), num_t_oldbnd )

          dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa_in(mfrm)))**onethird_in
          if (dgn_t_new .le. dgnum_amode_in(mfrm)) cycle

          lndgn_new = log( dgn_t_new )
          lndgv_new = lndgn_new + dum3alnsg2_in(ipair)
          yn_tail = (lndp_cut_in(ipair) - lndgn_new)*factoryy_in(mfrm)
          yv_tail = (lndp_cut_in(ipair) - lndgv_new)*factoryy_in(mfrm)
          tailfr_numnew = 0.5_r8*erfc( yn_tail )
          tailfr_volnew = 0.5_r8*erfc( yv_tail )

          dgn_t_old = (dryvol_t_oldbnd/(num_t_oldbnd*factoraa_in(mfrm)))**onethird_in
          if (dgn_t_new .ge. dp_cut_in(ipair)) then
             dgn_t_old = min( dgn_t_old, dp_belowcut_in(ipair) )
          end if
          lndgn_old = log( dgn_t_old )
          lndgv_old = lndgn_old + dum3alnsg2_in(ipair)
          yn_tail = (lndp_cut_in(ipair) - lndgn_old)*factoryy_in(mfrm)
          yv_tail = (lndp_cut_in(ipair) - lndgv_old)*factoryy_in(mfrm)
          tailfr_numold = 0.5_r8*erfc( yn_tail )
          tailfr_volold = 0.5_r8*erfc( yv_tail )

          dum = tailfr_volnew*dryvol_t_new - tailfr_volold*dryvol_t_old
          if (dum .le. 0.0_r8) cycle

          xferfrac_vol_out(i,k,ipair) = min( dum, dryvol_t_new )/dryvol_t_new
          xferfrac_vol_out(i,k,ipair) = min( xferfrac_vol_out(i,k,ipair), xferfrac_max_in )
          xferfrac_num_out(i,k,ipair) = tailfr_numnew - tailfr_numold
          xferfrac_num_out(i,k,ipair) = max( 0.0_r8, min( xferfrac_num_out(i,k,ipair), xferfrac_vol_out(i,k,ipair) ) )
       end do
       end do
    end do

  end subroutine modal_aero_rename_no_acc_crs_xferfracs_native

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_no_acc_crs_tendencies( ncol, loffset, npair_renamexf_in, &
       deltat, deltatinv, is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, &
       pdel, q, qqcw, dqdt, dqqcwdt, qsrflx, qqcwsrflx, xferfrac_vol_ik, xferfrac_num_ik, &
       nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, &
       lspectooa_renamexf_in, lspectooc_renamexf_in )

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use physconst, only: gravit

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    integer, intent(in) :: npair_renamexf_in
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: deltatinv
    logical, intent(in) :: is_dorename_atik
    logical, intent(in) :: dorename_atik(ncol,pver)
    integer, intent(in) :: jsrflx_rename
    integer, intent(in) :: nsrflx
    real(r8), target, intent(in) :: pdel(pcols,pver)
    real(r8), target, intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), target, intent(inout) :: dqdt(ncol,pver,pcnstxx)
    real(r8), target, intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), target, intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
    real(r8), target, intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
    real(r8), target, intent(in) :: xferfrac_vol_ik(ncol,pver,maxpair_renamexf)
    real(r8), target, intent(in) :: xferfrac_num_ik(ncol,pver,maxpair_renamexf)
    integer, target, intent(in) :: nspecfrm_renamexf_in(maxpair_renamexf)
    integer, target, intent(in) :: lspecfrma_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, target, intent(in) :: lspecfrmc_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, target, intent(in) :: lspectooa_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, target, intent(in) :: lspectooc_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer(c_int64_t), target :: dorename_atik_c(ncol,pver)
    integer(c_int64_t), target :: nspecfrm_renamexf_c(maxpair_renamexf)
    integer(c_int64_t), target :: lspecfrma_renamexf_c(maxspec_renamexf,maxpair_renamexf)
    integer(c_int64_t), target :: lspecfrmc_renamexf_c(maxspec_renamexf,maxpair_renamexf)
    integer(c_int64_t), target :: lspectooa_renamexf_c(maxspec_renamexf,maxpair_renamexf)
    integer(c_int64_t), target :: lspectooc_renamexf_c(maxspec_renamexf,maxpair_renamexf)
    integer(c_int64_t) :: is_dorename_atik_c
    integer :: i, k, ipair, iq

    interface
       subroutine modal_aero_rename_no_acc_crs_tendencies_codon( &
            ncol_c, pcols_c, pver_c, pcnstxx_c, maxpair_renamexf_c, maxspec_renamexf_c, &
            loffset_c, npair_renamexf_c, jsrflx_rename_c, nsrflx_c, is_dorename_atik_c, &
            deltat_c, deltatinv_c, gravit_c, pdel_p, dorename_atik_p, q_p, qqcw_p, dqdt_p, &
            dqqcwdt_p, qsrflx_p, qqcwsrflx_p, xferfrac_vol_p, xferfrac_num_p, nspecfrm_renamexf_p, &
            lspecfrma_renamexf_p, lspecfrmc_renamexf_p, lspectooa_renamexf_p, lspectooc_renamexf_p ) &
            bind(c, name="modal_aero_rename_no_acc_crs_tendencies_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnstxx_c, maxpair_renamexf_c
         integer(c_int64_t), value :: maxspec_renamexf_c, loffset_c, npair_renamexf_c
         integer(c_int64_t), value :: jsrflx_rename_c, nsrflx_c, is_dorename_atik_c
         real(c_double), value :: deltat_c, deltatinv_c, gravit_c
         type(c_ptr), value :: pdel_p, dorename_atik_p, q_p, qqcw_p, dqdt_p, dqqcwdt_p
         type(c_ptr), value :: qsrflx_p, qqcwsrflx_p, xferfrac_vol_p, xferfrac_num_p
         type(c_ptr), value :: nspecfrm_renamexf_p, lspecfrma_renamexf_p, lspecfrmc_renamexf_p
         type(c_ptr), value :: lspectooa_renamexf_p, lspectooc_renamexf_p
       end subroutine modal_aero_rename_no_acc_crs_tendencies_codon
    end interface

    call modal_aero_rename_no_acc_crs_tendencies_select_impl()

    if (.not. modal_aero_rename_no_acc_crs_tendencies_use_native_impl) then
       if (is_dorename_atik) then
          is_dorename_atik_c = 1_c_int64_t
          do k = 1, pver
             do i = 1, ncol
                if (dorename_atik(i,k)) then
                   dorename_atik_c(i,k) = 1_c_int64_t
                else
                   dorename_atik_c(i,k) = 0_c_int64_t
                end if
             end do
          end do
       else
          is_dorename_atik_c = 0_c_int64_t
          dorename_atik_c(:,:) = 0_c_int64_t
       end if

       do ipair = 1, maxpair_renamexf
          nspecfrm_renamexf_c(ipair) = int(nspecfrm_renamexf_in(ipair), c_int64_t)
          do iq = 1, maxspec_renamexf
             lspecfrma_renamexf_c(iq,ipair) = int(lspecfrma_renamexf_in(iq,ipair), c_int64_t)
             lspecfrmc_renamexf_c(iq,ipair) = int(lspecfrmc_renamexf_in(iq,ipair), c_int64_t)
             lspectooa_renamexf_c(iq,ipair) = int(lspectooa_renamexf_in(iq,ipair), c_int64_t)
             lspectooc_renamexf_c(iq,ipair) = int(lspectooc_renamexf_in(iq,ipair), c_int64_t)
          end do
       end do

       call modal_aero_rename_no_acc_crs_tendencies_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnstxx, c_int64_t), &
            int(maxpair_renamexf, c_int64_t), int(maxspec_renamexf, c_int64_t), int(loffset, c_int64_t), &
            int(npair_renamexf_in, c_int64_t), int(jsrflx_rename, c_int64_t), int(nsrflx, c_int64_t), &
            is_dorename_atik_c, real(deltat, c_double), real(deltatinv, c_double), real(gravit, c_double), &
            c_loc(pdel(1,1)), c_loc(dorename_atik_c(1,1)), c_loc(q(1,1,1)), c_loc(qqcw(1,1,1)), c_loc(dqdt(1,1,1)), &
            c_loc(dqqcwdt(1,1,1)), c_loc(qsrflx(1,1,1)), c_loc(qqcwsrflx(1,1,1)), c_loc(xferfrac_vol_ik(1,1,1)), &
            c_loc(xferfrac_num_ik(1,1,1)), c_loc(nspecfrm_renamexf_c(1)), c_loc(lspecfrma_renamexf_c(1,1)), &
            c_loc(lspecfrmc_renamexf_c(1,1)), c_loc(lspectooa_renamexf_c(1,1)), c_loc(lspectooc_renamexf_c(1,1)) &
       )
       return
    end if

    call modal_aero_rename_no_acc_crs_tendencies_native( ncol, loffset, npair_renamexf_in, &
         deltat, deltatinv, is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, &
         pdel, q, qqcw, dqdt, dqqcwdt, qsrflx, qqcwsrflx, xferfrac_vol_ik, xferfrac_num_ik, &
         nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, &
         lspectooa_renamexf_in, lspectooc_renamexf_in )

  end subroutine modal_aero_rename_no_acc_crs_tendencies

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_no_acc_crs_tendencies_native( ncol, loffset, npair_renamexf_in, &
       deltat, deltatinv, is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, &
       pdel, q, qqcw, dqdt, dqqcwdt, qsrflx, qqcwsrflx, xferfrac_vol_ik, xferfrac_num_ik, &
       nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, &
       lspectooa_renamexf_in, lspectooc_renamexf_in )

    use physconst, only: gravit

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    integer, intent(in) :: npair_renamexf_in
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: deltatinv
    logical, intent(in) :: is_dorename_atik
    logical, intent(in) :: dorename_atik(ncol,pver)
    integer, intent(in) :: jsrflx_rename
    integer, intent(in) :: nsrflx
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: dqdt(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
    real(r8), intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
    real(r8), intent(in) :: xferfrac_vol_ik(ncol,pver,maxpair_renamexf)
    real(r8), intent(in) :: xferfrac_num_ik(ncol,pver,maxpair_renamexf)
    integer, intent(in) :: nspecfrm_renamexf_in(maxpair_renamexf)
    integer, intent(in) :: lspecfrma_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, intent(in) :: lspecfrmc_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, intent(in) :: lspectooa_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, intent(in) :: lspectooc_renamexf_in(maxspec_renamexf,maxpair_renamexf)

    integer :: i, j, k, ipair, iq
    integer :: lsfrma, lsfrmc, lstooa, lstooc
    real(r8) :: pdel_fac, xfercoef, xfertend, xferfrac_vol, xferfrac_num

mainloop1_k:  do k = 1, pver
mainloop1_i:  do i = 1, ncol

       if (is_dorename_atik) then
          if (.not. dorename_atik(i,k)) cycle mainloop1_i
       end if
       pdel_fac = pdel(i,k)/gravit

mainloop1_ipair:  do ipair = 1, npair_renamexf_in

          xferfrac_vol = xferfrac_vol_ik(i,k,ipair)
          xferfrac_num = xferfrac_num_ik(i,k,ipair)
          if (xferfrac_vol .le. 0.0_r8) cycle mainloop1_ipair

          j = jsrflx_rename
          do iq = 1, nspecfrm_renamexf_in(ipair)
             xfercoef = xferfrac_vol*deltatinv
             if (iq .eq. 1) xfercoef = xferfrac_num*deltatinv

             lsfrma = lspecfrma_renamexf_in(iq,ipair)-loffset
             lsfrmc = lspecfrmc_renamexf_in(iq,ipair)-loffset
             lstooa = lspectooa_renamexf_in(iq,ipair)-loffset
             lstooc = lspectooc_renamexf_in(iq,ipair)-loffset

             if (lsfrma .gt. 0) then
                xfertend = xfercoef*max( 0.0_r8, (q(i,k,lsfrma)+dqdt(i,k,lsfrma)*deltat) )
                dqdt(i,k,lsfrma) = dqdt(i,k,lsfrma) - xfertend
                qsrflx(i,lsfrma,j) = qsrflx(i,lsfrma,j) - xfertend*pdel_fac
                if (lstooa .gt. 0) then
                   dqdt(i,k,lstooa) = dqdt(i,k,lstooa) + xfertend
                   qsrflx(i,lstooa,j) = qsrflx(i,lstooa,j) + xfertend*pdel_fac
                end if
             end if

             if (lsfrmc .gt. 0) then
                xfertend = xfercoef*max( 0.0_r8, (qqcw(i,k,lsfrmc)+dqqcwdt(i,k,lsfrmc)*deltat) )
                dqqcwdt(i,k,lsfrmc) = dqqcwdt(i,k,lsfrmc) - xfertend
                qqcwsrflx(i,lsfrmc,j) = qqcwsrflx(i,lsfrmc,j) - xfertend*pdel_fac
                if (lstooc .gt. 0) then
                   dqqcwdt(i,k,lstooc) = dqqcwdt(i,k,lstooc) + xfertend
                   qqcwsrflx(i,lstooc,j) = qqcwsrflx(i,lstooc,j) + xfertend*pdel_fac
                end if
             end if
          end do

       end do mainloop1_ipair

    end do mainloop1_i
    end do mainloop1_k

  end subroutine modal_aero_rename_no_acc_crs_tendencies_native

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_rename_no_acc_crs_sub --- ...
!
! !INTERFACE:
	subroutine modal_aero_rename_no_acc_crs_sub(                       &
                        fromwhere,         lchnk,               &
                        ncol,              nstep,               &
                        loffset,           deltat,              &
                        pdel,                                   &
                        dotendrn,          q,                   &
                        dqdt,              dqdt_other,          &
                        dotendqqcwrn,      qqcw,                &
                        dqqcwdt,           dqqcwdt_other,       &
                        is_dorename_atik,  dorename_atik,       &
                        jsrflx_rename,     nsrflx,              &
                        qsrflx,            qqcwsrflx            )

! !USES:
   use physconst, only: gravit, mwdry
   use units, only: getunit
   use shr_spfn_mod, only: erfc => shr_spfn_erfc

   implicit none


! !PARAMETERS:
   character(len=*), intent(in) :: fromwhere    ! identifies which module
                                                ! is making the call
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   real(r8), intent(in)    :: deltat               ! time step (s)

   real(r8), intent(in)    :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
   real(r8), intent(in)    :: q(ncol,pver,pcnstxx) ! tracer mixing ratio array
                                                   ! *** MUST BE mol/mol-air or #/mol-air
                                                   ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: qqcw(ncol,pver,pcnstxx) ! like q but for cloud-borne species

   real(r8), intent(inout) :: dqdt(ncol,pver,pcnstxx)  ! TMR tendency array;
                              ! incoming dqdt = tendencies for the 
                              !     "fromwhere" continuous growth process 
                              ! the renaming tendencies are added on
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
   real(r8), intent(in)    :: dqdt_other(ncol,pver,pcnstxx)  
                              ! tendencies for "other" continuous growth process 
                              ! currently in cam3
                              !     dqdt is from gas (h2so4, nh3) condensation
                              !     dqdt_other is from aqchem and soa
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), intent(in)    :: dqqcwdt_other(ncol,pver,pcnstxx)  
   logical,  intent(inout) :: dotendrn(pcnstxx) ! identifies the species for which
                              !     renaming dqdt is computed
   logical,  intent(inout) :: dotendqqcwrn(pcnstxx)

   logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
   logical,  intent(in)    :: dorename_atik(ncol,pver) ! true if renaming should
                                                        ! be done at i,k
   integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
   integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

   real(r8), intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
                              ! process-specific column tracer tendencies 
   real(r8), intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)

! !DESCRIPTION: 
! computes TMR (tracer mixing ratio) tendencies for "mode renaming"
!    during a continuous growth process
! currently this transfers number and mass (and surface) from the aitken
!    to accumulation mode after gas condensation or stratiform-cloud
!    aqueous chemistry
! (convective cloud aqueous chemistry not yet implemented)
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: ldiag1=-1
   integer :: i, icol_diag, ipair, iq, j, k, l, l1, l2, la, lc, lunout
   integer :: lsfrma, lsfrmc, lstooa, lstooc
   integer :: mfrm, mtoo, n, n1, n2, ntot_msa_a
   integer :: idomode(ntot_amode)
   integer, save :: lun = -1  ! logical unit for diagnostics (6, or other
                              ! if a special diagnostics file is opened)


   real (r8) :: deldryvol_a(ncol,pver,ntot_amode)
   real (r8) :: deldryvol_c(ncol,pver,ntot_amode)
   real (r8) :: deltatinv
   real (r8) :: dp_belowcut(maxpair_renamexf)
   real (r8) :: dp_cut(maxpair_renamexf)
   real (r8) :: dgn_aftr, dgn_xfer
   real (r8) :: dgn_t_new, dgn_t_old
   real (r8) :: dryvol_t_del, dryvol_t_new
   real (r8) :: dryvol_t_old, dryvol_t_oldbnd
   real (r8) :: dryvol_a(ncol,pver,ntot_amode)
   real (r8) :: dryvol_c(ncol,pver,ntot_amode)
   real (r8) :: dryvol_smallest(ntot_amode)
   real (r8) :: dum
   real (r8) :: dum3alnsg2(maxpair_renamexf)
   real (r8) :: dum_m2v, dum_m2vdt
   real (r8) :: factoraa(ntot_amode)
   real (r8) :: factoryy(ntot_amode)
   real (r8) :: frelax
   real (r8) :: lndp_cut(maxpair_renamexf)
   real (r8) :: lndgn_new, lndgn_old
   real (r8) :: lndgv_new, lndgv_old
   real (r8) :: num_t_old, num_t_oldbnd
   real (r8) :: onethird
   real (r8) :: pdel_fac
   real (r8) :: tailfr_volnew, tailfr_volold
   real (r8) :: tailfr_numnew, tailfr_numold
   real (r8) :: v2nhirlx(ntot_amode), v2nlorlx(ntot_amode)
   real (r8) :: xfercoef, xfertend
   real (r8) :: xferfrac_vol, xferfrac_num, xferfrac_max
   real (r8) :: xferfrac_vol_ik(ncol,pver,maxpair_renamexf)
   real (r8) :: xferfrac_num_ik(ncol,pver,maxpair_renamexf)

   real (r8) :: yn_tail, yv_tail

! begin
	lunout = iulog

!   get logical unit (for output to dumpconv, deactivate the "lun = 6")
 	lun = iulog
	if (lun < 1) then
	   lun = getunit()
 	   open( unit=lun, file='dump.rename',   &
 			status='unknown', form='formatted' )
	end if


!
!   calculations done once on initial entry
!
!   "init" is now done through chem_init (and things under it)
!	if (npair_renamexf .eq. -123456789) then
!	    npair_renamexf = 0
!	    call modal_aero_rename_init
!	end if

!
!   check if any renaming pairs exist
!
	if (npair_renamexf .le. 0) return
! 	if (ncol .ne. -123456789) return
!	if (fromwhere .eq. 'aqchem') return

!
!   compute aerosol dry-volume for the "from mode" of each renaming pair
!   also compute dry-volume change during the continuous growth process
!	using the incoming dqdt*deltat
!
	deltatinv = 1.0_r8/(deltat*(1.0_r8 + 1.0e-15_r8))
	onethird = 1.0_r8/3.0_r8
	frelax = 27.0_r8
	xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

	do n = 1, ntot_amode
	    idomode(n) = 0
	end do

	do ipair = 1, npair_renamexf
	    if (ipair .gt. 1) goto 8100
	    idomode(modefrm_renamexf(ipair)) = 1

	    mfrm = modefrm_renamexf(ipair)
	    mtoo = modetoo_renamexf(ipair)
	    factoraa(mfrm) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mfrm)**2))
	    factoraa(mtoo) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mtoo)**2))
	    factoryy(mfrm) = sqrt( 0.5_r8 )/alnsg_amode(mfrm)
!   dryvol_smallest is a very small volume mixing ratio (m3-AP/kmol-air)
!   used for avoiding overflow.  it corresponds to dp = 1 nm
!   and number = 1e-5 #/mg-air ~= 1e-5 #/cm3-air
	    dryvol_smallest(mfrm) = 1.0e-25_r8
	    v2nlorlx(mfrm) = voltonumblo_amode(mfrm)*frelax
	    v2nhirlx(mfrm) = voltonumbhi_amode(mfrm)/frelax

	    dum3alnsg2(ipair) = 3.0_r8 * (alnsg_amode(mfrm)**2)
	    dp_cut(ipair) = sqrt(   &
		dgnum_amode(mfrm)*exp(1.5_r8*(alnsg_amode(mfrm)**2)) *   &
		dgnum_amode(mtoo)*exp(1.5_r8*(alnsg_amode(mtoo)**2)) )
	    lndp_cut(ipair) = log( dp_cut(ipair) )
	    dp_belowcut(ipair) = 0.99_r8*dp_cut(ipair)
	end do

	call modal_aero_rename_no_acc_crs_dryvols( ncol, loffset, deltat, &
	     idomode, q, qqcw, dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, &
	     nspec_amode, lspectype_amode, specmw_amode, specdens_amode, &
	     lmassptr_amode, lmassptrcw_amode, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c )

	call modal_aero_rename_no_acc_crs_xferfracs( ncol, loffset, npair_renamexf, &
	     q, qqcw, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, &
	     modefrm_renamexf, modetoo_renamexf, numptr_amode, numptrcw_amode, dgnum_amode, &
	     factoraa, factoryy, dryvol_smallest, v2nlorlx, v2nhirlx, dum3alnsg2, dp_cut, lndp_cut, dp_belowcut, &
	     onethird, xferfrac_max, xferfrac_vol_ik, xferfrac_num_ik )

	call modal_aero_rename_no_acc_crs_tendencies( ncol, loffset, npair_renamexf, &
	     deltat, deltatinv, is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, &
	     pdel, q, qqcw, dqdt, dqqcwdt, qsrflx, qqcwsrflx, xferfrac_vol_ik, xferfrac_num_ik, &
	     nspecfrm_renamexf, lspecfrma_renamexf, lspecfrmc_renamexf, lspectooa_renamexf, lspectooc_renamexf )

!
!   set dotend's
!
	call modal_aero_rename_set_dotend_flags( loffset, npair_renamexf, nspecfrm_renamexf, &
	     lspecfrma_renamexf, lspecfrmc_renamexf, lspectooa_renamexf, lspectooc_renamexf, &
	     dotendrn, dotendqqcwrn )


	return


!
!   error -- renaming currently just works for 1 pair
!
8100	write(lunout,9050) ipair
	call endrun( 'modal_aero_rename_no_acc_crs_sub error' )
9050	format( / '*** subr. modal_aero_rename_no_acc_crs_sub ***' /   &
      	    4x, 'aerosol renaming not implemented for ipair =', i5 )

!EOC
	end subroutine modal_aero_rename_no_acc_crs_sub

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_dryvols( ncol, loffset, deltat, ixferable_all_in, &
       q, qqcw, dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, nspec_mfrm_in, &
       lspectype_mfrm_in, specmw_amode_in, specdens_amode_in, lmassptr_mfrm_in, &
       lmassptrcw_mfrm_in, ixferable_a_in, ixferable_c_in, dryvol_a, dryvol_c, &
       deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab )

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    integer, intent(in) :: ixferable_all_in
    real(r8), intent(in) :: deltat
    real(r8), target, intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqdt(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqdt_other(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqqcwdt_other(ncol,pver,pcnstxx)
    integer, intent(in) :: nspec_mfrm_in
    integer, target, intent(in) :: lspectype_mfrm_in(maxspec_renamexf)
    real(r8), target, intent(in) :: specmw_amode_in(*)
    real(r8), target, intent(in) :: specdens_amode_in(*)
    integer, target, intent(in) :: lmassptr_mfrm_in(maxspec_renamexf)
    integer, target, intent(in) :: lmassptrcw_mfrm_in(maxspec_renamexf)
    integer, target, intent(in) :: ixferable_a_in(maxspec_renamexf)
    integer, target, intent(in) :: ixferable_c_in(maxspec_renamexf)
    real(r8), target, intent(out) :: dryvol_a(ncol,pver)
    real(r8), target, intent(out) :: dryvol_c(ncol,pver)
    real(r8), target, intent(out) :: deldryvol_a(ncol,pver)
    real(r8), target, intent(out) :: deldryvol_c(ncol,pver)
    real(r8), target, intent(out) :: dryvol_a_xfab(ncol,pver)
    real(r8), target, intent(out) :: dryvol_c_xfab(ncol,pver)
    integer(c_int64_t) :: ixferable_all_c
    integer(c_int64_t), target :: lspectype_mfrm_c(maxspec_renamexf)
    integer(c_int64_t), target :: lmassptr_mfrm_c(maxspec_renamexf)
    integer(c_int64_t), target :: lmassptrcw_mfrm_c(maxspec_renamexf)
    integer(c_int64_t), target :: ixferable_a_c(maxspec_renamexf)
    integer(c_int64_t), target :: ixferable_c_c(maxspec_renamexf)
    integer :: l1

    interface
       subroutine modal_aero_rename_acc_crs_dryvols_codon( &
            ncol_c, pver_c, pcnstxx_c, maxspec_renamexf_c, loffset_c, ixferable_all_c, &
            nspec_mfrm_c, deltat_c, q_p, qqcw_p, dqdt_p, dqdt_other_p, dqqcwdt_p, dqqcwdt_other_p, &
            lspectype_mfrm_p, specmw_amode_p, specdens_amode_p, lmassptr_mfrm_p, lmassptrcw_mfrm_p, &
            ixferable_a_p, ixferable_c_p, dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p, &
            dryvol_a_xfab_p, dryvol_c_xfab_p ) bind(c, name="modal_aero_rename_acc_crs_dryvols_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pver_c, pcnstxx_c, maxspec_renamexf_c
         integer(c_int64_t), value :: loffset_c, ixferable_all_c, nspec_mfrm_c
         real(c_double), value :: deltat_c
         type(c_ptr), value :: q_p, qqcw_p, dqdt_p, dqdt_other_p, dqqcwdt_p, dqqcwdt_other_p
         type(c_ptr), value :: lspectype_mfrm_p, specmw_amode_p, specdens_amode_p
         type(c_ptr), value :: lmassptr_mfrm_p, lmassptrcw_mfrm_p, ixferable_a_p, ixferable_c_p
         type(c_ptr), value :: dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p
         type(c_ptr), value :: dryvol_a_xfab_p, dryvol_c_xfab_p
       end subroutine modal_aero_rename_acc_crs_dryvols_codon
    end interface

    call modal_aero_rename_acc_crs_dryvols_select_impl()

    if (.not. modal_aero_rename_acc_crs_dryvols_use_native_impl) then
       ixferable_all_c = int(ixferable_all_in, c_int64_t)
       do l1 = 1, maxspec_renamexf
          lspectype_mfrm_c(l1) = int(lspectype_mfrm_in(l1), c_int64_t)
          lmassptr_mfrm_c(l1) = int(lmassptr_mfrm_in(l1), c_int64_t)
          lmassptrcw_mfrm_c(l1) = int(lmassptrcw_mfrm_in(l1), c_int64_t)
          ixferable_a_c(l1) = int(ixferable_a_in(l1), c_int64_t)
          ixferable_c_c(l1) = int(ixferable_c_in(l1), c_int64_t)
       end do

       call modal_aero_rename_acc_crs_dryvols_codon( &
            int(ncol, c_int64_t), int(pver, c_int64_t), int(pcnstxx, c_int64_t), int(maxspec_renamexf, c_int64_t), &
            int(loffset, c_int64_t), ixferable_all_c, int(nspec_mfrm_in, c_int64_t), real(deltat, c_double), &
            c_loc(q(1,1,1)), c_loc(qqcw(1,1,1)), c_loc(dqdt(1,1,1)), c_loc(dqdt_other(1,1,1)), &
            c_loc(dqqcwdt(1,1,1)), c_loc(dqqcwdt_other(1,1,1)), c_loc(lspectype_mfrm_c(1)), c_loc(specmw_amode_in(1)), &
            c_loc(specdens_amode_in(1)), c_loc(lmassptr_mfrm_c(1)), c_loc(lmassptrcw_mfrm_c(1)), c_loc(ixferable_a_c(1)), &
            c_loc(ixferable_c_c(1)), c_loc(dryvol_a(1,1)), c_loc(dryvol_c(1,1)), c_loc(deldryvol_a(1,1)), &
            c_loc(deldryvol_c(1,1)), c_loc(dryvol_a_xfab(1,1)), c_loc(dryvol_c_xfab(1,1)) &
       )
       return
    end if

    call modal_aero_rename_acc_crs_dryvols_native( ncol, loffset, deltat, ixferable_all_in, &
         q, qqcw, dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, nspec_mfrm_in, lspectype_mfrm_in, &
         specmw_amode_in, specdens_amode_in, lmassptr_mfrm_in, lmassptrcw_mfrm_in, ixferable_a_in, &
         ixferable_c_in, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab )

  end subroutine modal_aero_rename_acc_crs_dryvols

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_dryvols_native( ncol, loffset, deltat, ixferable_all_in, &
       q, qqcw, dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, nspec_mfrm_in, lspectype_mfrm_in, &
       specmw_amode_in, specdens_amode_in, lmassptr_mfrm_in, lmassptrcw_mfrm_in, ixferable_a_in, &
       ixferable_c_in, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab )

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    integer, intent(in) :: ixferable_all_in
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqdt(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqdt_other(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqqcwdt_other(ncol,pver,pcnstxx)
    integer, intent(in) :: nspec_mfrm_in
    integer, intent(in) :: lspectype_mfrm_in(maxspec_renamexf)
    real(r8), intent(in) :: specmw_amode_in(*)
    real(r8), intent(in) :: specdens_amode_in(*)
    integer, intent(in) :: lmassptr_mfrm_in(maxspec_renamexf)
    integer, intent(in) :: lmassptrcw_mfrm_in(maxspec_renamexf)
    integer, intent(in) :: ixferable_a_in(maxspec_renamexf)
    integer, intent(in) :: ixferable_c_in(maxspec_renamexf)
    real(r8), intent(out) :: dryvol_a(ncol,pver)
    real(r8), intent(out) :: dryvol_c(ncol,pver)
    real(r8), intent(out) :: deldryvol_a(ncol,pver)
    real(r8), intent(out) :: deldryvol_c(ncol,pver)
    real(r8), intent(out) :: dryvol_a_xfab(ncol,pver)
    real(r8), intent(out) :: dryvol_c_xfab(ncol,pver)

    integer :: l1, l2, la, lc
    real(r8) :: tmp_m2v, tmp_m2vdt

    dryvol_a(:,:) = 0.0_r8
    dryvol_c(:,:) = 0.0_r8
    deldryvol_a(:,:) = 0.0_r8
    deldryvol_c(:,:) = 0.0_r8
    dryvol_a_xfab(:,:) = 0.0_r8
    dryvol_c_xfab(:,:) = 0.0_r8

    do l1 = 1, nspec_mfrm_in
       l2 = lspectype_mfrm_in(l1)
       tmp_m2v = specmw_amode_in(l2) / specdens_amode_in(l2)
       tmp_m2vdt = tmp_m2v*deltat
       la = lmassptr_mfrm_in(l1)-loffset
       if (la > 0) then
          dryvol_a(1:ncol,:) = dryvol_a(1:ncol,:)    &
               + tmp_m2v*max( 0.0_r8, q(1:ncol,:,la)-deltat*dqdt_other(1:ncol,:,la) )
          deldryvol_a(1:ncol,:) = deldryvol_a(1:ncol,:)    &
               + (dqdt_other(1:ncol,:,la) + dqdt(1:ncol,:,la))*tmp_m2vdt
          if ( (ixferable_all_in <= 0) .and. (ixferable_a_in(l1) > 0) ) then
             dryvol_a_xfab(1:ncol,:) = dryvol_a_xfab(1:ncol,:)    &
                  + tmp_m2v*max( 0.0_r8, q(1:ncol,:,la)+deltat*dqdt(1:ncol,:,la) )
          end if
       end if

       lc = lmassptrcw_mfrm_in(l1)-loffset
       if (lc > 0) then
          dryvol_c(1:ncol,:) = dryvol_c(1:ncol,:)    &
               + tmp_m2v*max( 0.0_r8, qqcw(1:ncol,:,lc)-deltat*dqqcwdt_other(1:ncol,:,lc) )
          deldryvol_c(1:ncol,:) = deldryvol_c(1:ncol,:)    &
               + (dqqcwdt_other(1:ncol,:,lc) + dqqcwdt(1:ncol,:,lc))*tmp_m2vdt
          if ( (ixferable_all_in <= 0) .and. (ixferable_c_in(l1) > 0) ) then
             dryvol_c_xfab(1:ncol,:) = dryvol_c_xfab(1:ncol,:)    &
                  + tmp_m2v*max( 0.0_r8, qqcw(1:ncol,:,lc)+deltat*dqqcwdt(1:ncol,:,lc) )
          end if
       end if
    end do

  end subroutine modal_aero_rename_acc_crs_dryvols_native

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_xferfracs( ncol, loffset, troplev, q, qqcw, &
       dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab, &
       mfrm_in, numptr_amode_mfrm_in, numptrcw_amode_mfrm_in, dgnum_amode_mfrm_in, &
       factoraa_in, factoryy_in, dryvol_smallest_in, v2nlorlx_in, v2nhirlx_in, &
       factor_3alnsg2_in, dp_cut_in, lndp_cut_in, dp_belowcut_in, dp_xfernone_thresh_in, &
       dp_xferall_thresh_in, flagaa_shrink_in, igrow_shrink_in, ixferable_all_in, &
       method_optbb_in, onethird_in, xferfrac_max_in, xferfrac_vol_out, xferfrac_num_out )

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    integer, target, intent(in) :: troplev(pcols)
    real(r8), target, intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dryvol_a(ncol,pver)
    real(r8), target, intent(in) :: dryvol_c(ncol,pver)
    real(r8), target, intent(in) :: deldryvol_a(ncol,pver)
    real(r8), target, intent(in) :: deldryvol_c(ncol,pver)
    real(r8), target, intent(in) :: dryvol_a_xfab(ncol,pver)
    real(r8), target, intent(in) :: dryvol_c_xfab(ncol,pver)
    integer, intent(in) :: mfrm_in
    integer, intent(in) :: numptr_amode_mfrm_in
    integer, intent(in) :: numptrcw_amode_mfrm_in
    real(r8), intent(in) :: dgnum_amode_mfrm_in
    real(r8), intent(in) :: factoraa_in
    real(r8), intent(in) :: factoryy_in
    real(r8), intent(in) :: dryvol_smallest_in
    real(r8), intent(in) :: v2nlorlx_in
    real(r8), intent(in) :: v2nhirlx_in
    real(r8), intent(in) :: factor_3alnsg2_in
    real(r8), intent(in) :: dp_cut_in
    real(r8), intent(in) :: lndp_cut_in
    real(r8), intent(in) :: dp_belowcut_in
    real(r8), intent(in) :: dp_xfernone_thresh_in
    real(r8), intent(in) :: dp_xferall_thresh_in
    logical, intent(in) :: flagaa_shrink_in
    integer, intent(in) :: igrow_shrink_in
    integer, intent(in) :: ixferable_all_in
    integer, intent(in) :: method_optbb_in
    real(r8), intent(in) :: onethird_in
    real(r8), intent(in) :: xferfrac_max_in
    real(r8), target, intent(out) :: xferfrac_vol_out(ncol,pver)
    real(r8), target, intent(out) :: xferfrac_num_out(ncol,pver)
    integer(c_int64_t), target :: troplev_c(pcols)
    integer(c_int64_t) :: flagaa_shrink_c
    integer :: i

    interface
       subroutine modal_aero_rename_acc_crs_xferfracs_codon( &
            ncol_c, pcols_c, pver_c, pcnstxx_c, loffset_c, mfrm_c, numptr_amode_mfrm_c, &
            numptrcw_amode_mfrm_c, igrow_shrink_c, ixferable_all_c, method_optbb_c, flagaa_shrink_c, &
            dgnum_amode_mfrm_c, factoraa_c, factoryy_c, dryvol_smallest_c, v2nlorlx_c, v2nhirlx_c, &
            factor_3alnsg2_c, dp_cut_c, lndp_cut_c, dp_belowcut_c, dp_xfernone_thresh_c, dp_xferall_thresh_c, &
            onethird_c, xferfrac_max_c, troplev_p, q_p, qqcw_p, dryvol_a_p, dryvol_c_p, deldryvol_a_p, &
            deldryvol_c_p, dryvol_a_xfab_p, dryvol_c_xfab_p, xferfrac_vol_p, xferfrac_num_p ) &
            bind(c, name="modal_aero_rename_acc_crs_xferfracs_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnstxx_c, loffset_c, mfrm_c
         integer(c_int64_t), value :: numptr_amode_mfrm_c, numptrcw_amode_mfrm_c
         integer(c_int64_t), value :: igrow_shrink_c, ixferable_all_c, method_optbb_c, flagaa_shrink_c
         real(c_double), value :: dgnum_amode_mfrm_c, factoraa_c, factoryy_c, dryvol_smallest_c
         real(c_double), value :: v2nlorlx_c, v2nhirlx_c, factor_3alnsg2_c, dp_cut_c, lndp_cut_c
         real(c_double), value :: dp_belowcut_c, dp_xfernone_thresh_c, dp_xferall_thresh_c
         real(c_double), value :: onethird_c, xferfrac_max_c
         type(c_ptr), value :: troplev_p, q_p, qqcw_p, dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p
         type(c_ptr), value :: dryvol_a_xfab_p, dryvol_c_xfab_p, xferfrac_vol_p, xferfrac_num_p
       end subroutine modal_aero_rename_acc_crs_xferfracs_codon
    end interface

    call modal_aero_rename_acc_crs_xferfracs_select_impl()

    if (.not. modal_aero_rename_acc_crs_xferfracs_use_native_impl) then
       do i = 1, pcols
          troplev_c(i) = int(troplev(i), c_int64_t)
       end do
       if (flagaa_shrink_in) then
          flagaa_shrink_c = 1_c_int64_t
       else
          flagaa_shrink_c = 0_c_int64_t
       end if

       call modal_aero_rename_acc_crs_xferfracs_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnstxx, c_int64_t), &
            int(loffset, c_int64_t), int(mfrm_in, c_int64_t), int(numptr_amode_mfrm_in, c_int64_t), &
            int(numptrcw_amode_mfrm_in, c_int64_t), int(igrow_shrink_in, c_int64_t), int(ixferable_all_in, c_int64_t), &
            int(method_optbb_in, c_int64_t), flagaa_shrink_c, real(dgnum_amode_mfrm_in, c_double), real(factoraa_in, c_double), &
            real(factoryy_in, c_double), real(dryvol_smallest_in, c_double), real(v2nlorlx_in, c_double), real(v2nhirlx_in, c_double), &
            real(factor_3alnsg2_in, c_double), real(dp_cut_in, c_double), real(lndp_cut_in, c_double), real(dp_belowcut_in, c_double), &
            real(dp_xfernone_thresh_in, c_double), real(dp_xferall_thresh_in, c_double), real(onethird_in, c_double), &
            real(xferfrac_max_in, c_double), c_loc(troplev_c(1)), c_loc(q(1,1,1)), c_loc(qqcw(1,1,1)), c_loc(dryvol_a(1,1)), &
            c_loc(dryvol_c(1,1)), c_loc(deldryvol_a(1,1)), c_loc(deldryvol_c(1,1)), c_loc(dryvol_a_xfab(1,1)), &
            c_loc(dryvol_c_xfab(1,1)), c_loc(xferfrac_vol_out(1,1)), c_loc(xferfrac_num_out(1,1)) &
       )
       return
    end if

    call modal_aero_rename_acc_crs_xferfracs_native( ncol, loffset, troplev, q, qqcw, dryvol_a, dryvol_c, &
         deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab, mfrm_in, numptr_amode_mfrm_in, &
         numptrcw_amode_mfrm_in, dgnum_amode_mfrm_in, factoraa_in, factoryy_in, dryvol_smallest_in, &
         v2nlorlx_in, v2nhirlx_in, factor_3alnsg2_in, dp_cut_in, lndp_cut_in, dp_belowcut_in, &
         dp_xfernone_thresh_in, dp_xferall_thresh_in, flagaa_shrink_in, igrow_shrink_in, ixferable_all_in, &
         method_optbb_in, onethird_in, xferfrac_max_in, xferfrac_vol_out, xferfrac_num_out )

  end subroutine modal_aero_rename_acc_crs_xferfracs

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_xferfracs_native( ncol, loffset, troplev, q, qqcw, dryvol_a, dryvol_c, &
       deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab, mfrm_in, numptr_amode_mfrm_in, &
       numptrcw_amode_mfrm_in, dgnum_amode_mfrm_in, factoraa_in, factoryy_in, dryvol_smallest_in, &
       v2nlorlx_in, v2nhirlx_in, factor_3alnsg2_in, dp_cut_in, lndp_cut_in, dp_belowcut_in, &
       dp_xfernone_thresh_in, dp_xferall_thresh_in, flagaa_shrink_in, igrow_shrink_in, ixferable_all_in, &
       method_optbb_in, onethird_in, xferfrac_max_in, xferfrac_vol_out, xferfrac_num_out )

    use shr_spfn_mod, only: erfc => shr_spfn_erfc

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    integer, intent(in) :: troplev(pcols)
    real(r8), intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dryvol_a(ncol,pver)
    real(r8), intent(in) :: dryvol_c(ncol,pver)
    real(r8), intent(in) :: deldryvol_a(ncol,pver)
    real(r8), intent(in) :: deldryvol_c(ncol,pver)
    real(r8), intent(in) :: dryvol_a_xfab(ncol,pver)
    real(r8), intent(in) :: dryvol_c_xfab(ncol,pver)
    integer, intent(in) :: mfrm_in
    integer, intent(in) :: numptr_amode_mfrm_in
    integer, intent(in) :: numptrcw_amode_mfrm_in
    real(r8), intent(in) :: dgnum_amode_mfrm_in
    real(r8), intent(in) :: factoraa_in
    real(r8), intent(in) :: factoryy_in
    real(r8), intent(in) :: dryvol_smallest_in
    real(r8), intent(in) :: v2nlorlx_in
    real(r8), intent(in) :: v2nhirlx_in
    real(r8), intent(in) :: factor_3alnsg2_in
    real(r8), intent(in) :: dp_cut_in
    real(r8), intent(in) :: lndp_cut_in
    real(r8), intent(in) :: dp_belowcut_in
    real(r8), intent(in) :: dp_xfernone_thresh_in
    real(r8), intent(in) :: dp_xferall_thresh_in
    logical, intent(in) :: flagaa_shrink_in
    integer, intent(in) :: igrow_shrink_in
    integer, intent(in) :: ixferable_all_in
    integer, intent(in) :: method_optbb_in
    real(r8), intent(in) :: onethird_in
    real(r8), intent(in) :: xferfrac_max_in
    real(r8), intent(out) :: xferfrac_vol_out(ncol,pver)
    real(r8), intent(out) :: xferfrac_num_out(ncol,pver)

    integer :: i, k
    logical :: flagbb_shrink
    real(r8) :: dgn_t_new, dgn_t_old, dgn_t_oldb
    real(r8) :: dryvol_t_del, dryvol_t_new, dryvol_t_new_xfab
    real(r8) :: dryvol_t_old, dryvol_t_oldb, dryvol_t_oldbnd
    real(r8) :: dryvol_xferamt
    real(r8) :: lndgn_new, lndgn_old
    real(r8) :: lndgv_new, lndgv_old
    real(r8) :: num_t_old, num_t_oldbnd
    real(r8) :: tailfr_volnew, tailfr_volold
    real(r8) :: tailfr_numnew, tailfr_numold
    real(r8) :: xferfrac_vol, xferfrac_num
    real(r8) :: yn_tail, yv_tail

    xferfrac_vol_out(:,:) = 0.0_r8
    xferfrac_num_out(:,:) = 0.0_r8

    do k = 1, pver
       do i = 1, ncol

          dryvol_t_old = dryvol_a(i,k) + dryvol_c(i,k)
          dryvol_t_del = deldryvol_a(i,k) + deldryvol_c(i,k)
          dryvol_t_new = dryvol_t_old + dryvol_t_del
          dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest_in )

grow_shrink_conditional1: &
          if (igrow_shrink_in > 0) then
             if (dryvol_t_new .le. dryvol_smallest_in) cycle
             if ( (method_optbb_in /= 2) .and. (dryvol_t_del .le. 1.0e-6_r8*dryvol_t_oldbnd) ) cycle

             num_t_old = q(i,k,numptr_amode_mfrm_in-loffset)
             num_t_old = num_t_old + qqcw(i,k,numptrcw_amode_mfrm_in-loffset)
             num_t_old = max( 0.0_r8, num_t_old )
             dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest_in )
             num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx_in, num_t_old )
             num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx_in, num_t_oldbnd )

             dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa_in))**onethird_in
             if (dgn_t_new .le. dp_xfernone_thresh_in) cycle

             dgn_t_old = (dryvol_t_oldbnd/(num_t_oldbnd*factoraa_in))**onethird_in
             dgn_t_oldb = dgn_t_old
             dryvol_t_oldb = dryvol_t_old
             if (method_optbb_in == 2) then
                if (dgn_t_old .ge. dp_cut_in) then
                   dryvol_t_oldb = dryvol_t_old * (dp_belowcut_in/dgn_t_old)**3
                   dgn_t_oldb = dp_belowcut_in
                end if
                if (dgn_t_new .lt. dp_xferall_thresh_in) then
                   if ((dryvol_t_new-dryvol_t_oldb) .le. 1.0e-6_r8*dryvol_t_oldbnd) cycle
                end if
             else if (dgn_t_new .ge. dp_cut_in) then
                dgn_t_oldb = min( dgn_t_oldb, dp_belowcut_in )
             end if

             lndgn_new = log( dgn_t_new )
             lndgv_new = lndgn_new + factor_3alnsg2_in
             yn_tail = (lndp_cut_in - lndgn_new)*factoryy_in
             yv_tail = (lndp_cut_in - lndgv_new)*factoryy_in
             tailfr_numnew = 0.5_r8*erfc( yn_tail )
             tailfr_volnew = 0.5_r8*erfc( yv_tail )

             lndgn_old = log( dgn_t_oldb )
             lndgv_old = lndgn_old + factor_3alnsg2_in
             yn_tail = (lndp_cut_in - lndgn_old)*factoryy_in
             yv_tail = (lndp_cut_in - lndgv_old)*factoryy_in
             tailfr_numold = 0.5_r8*erfc( yn_tail )
             tailfr_volold = 0.5_r8*erfc( yv_tail )

             if ( (method_optbb_in == 2) .and. (dgn_t_new .ge. dp_xferall_thresh_in) ) then
                dryvol_xferamt = dryvol_t_new
             else
                dryvol_xferamt = tailfr_volnew*dryvol_t_new - tailfr_volold*dryvol_t_oldb
             end if
             if (dryvol_xferamt .le. 0.0_r8) cycle

             xferfrac_vol = max( 0.0_r8, (dryvol_xferamt/dryvol_t_new) )
             if ( method_optbb_in == 2 .and. (xferfrac_vol >= xferfrac_max_in) ) then
                xferfrac_vol = 1.0_r8
                xferfrac_num = 1.0_r8
             else
                xferfrac_vol = min( xferfrac_vol, xferfrac_max_in )
                xferfrac_num = tailfr_numnew - tailfr_numold
                xferfrac_num = max( 0.0_r8, min( xferfrac_num, xferfrac_vol ) )
             end if

             if (ixferable_all_in <= 0) then
                dryvol_t_new_xfab = max( 0.0_r8, (dryvol_a_xfab(i,k) + dryvol_c_xfab(i,k)) )
                dryvol_xferamt = xferfrac_vol*dryvol_t_new
                if (dryvol_t_new_xfab >= 0.999999_r8*dryvol_xferamt) then
                   xferfrac_vol = min( 1.0_r8, (dryvol_xferamt/dryvol_t_new_xfab) )
                else if (dryvol_t_new_xfab >= 1.0e-7_r8*dryvol_xferamt) then
                   xferfrac_vol = 1.0_r8
                   xferfrac_num = xferfrac_num*(dryvol_t_new_xfab/dryvol_xferamt)
                else
                   cycle
                end if
             end if

          else grow_shrink_conditional1
             if (dryvol_t_old .le. dryvol_smallest_in) cycle

             if (dryvol_t_del .ge. -1.0e-6_r8*dryvol_t_oldbnd) then
                if ( flagaa_shrink_in .and. (k < troplev(i)) ) then
                   flagbb_shrink = .true.
                else
                   cycle
                end if
             else
                flagbb_shrink = .false.
             end if

             num_t_old = q(i,k,numptr_amode_mfrm_in-loffset)
             num_t_old = num_t_old + qqcw(i,k,numptrcw_amode_mfrm_in-loffset)
             num_t_old = max( 0.0_r8, num_t_old )
             dryvol_t_oldbnd = max( dryvol_t_old, dryvol_smallest_in )
             num_t_oldbnd = min( dryvol_t_oldbnd*v2nlorlx_in, num_t_old )
             num_t_oldbnd = max( dryvol_t_oldbnd*v2nhirlx_in, num_t_oldbnd )

             dgn_t_new = (dryvol_t_new/(num_t_oldbnd*factoraa_in))**onethird_in
             if (dgn_t_new .ge. dp_xfernone_thresh_in) cycle
             if (flagbb_shrink) then
                if (dgn_t_new .gt. dp_cut_in) cycle
             end if

             if ( dgn_t_new .le. dp_xferall_thresh_in ) then
                tailfr_numnew = 1.0_r8
                tailfr_volnew = 1.0_r8
             else
                lndgn_new = log( dgn_t_new )
                lndgv_new = lndgn_new + factor_3alnsg2_in
                yn_tail = (lndp_cut_in - lndgn_new)*factoryy_in
                yv_tail = (lndp_cut_in - lndgv_new)*factoryy_in
                tailfr_numnew = 1.0_r8 - 0.5_r8*erfc( yn_tail )
                tailfr_volnew = 1.0_r8 - 0.5_r8*erfc( yv_tail )
             end if

             dgn_t_old = (dryvol_t_oldbnd/(num_t_oldbnd*factoraa_in))**onethird_in
             dgn_t_oldb = dgn_t_old
             dryvol_t_oldb = dryvol_t_old

             tailfr_numold = 0.0_r8
             tailfr_volold = 0.0_r8

             xferfrac_vol = tailfr_volnew
             if (xferfrac_vol .le. 0.0_r8) cycle
             xferfrac_num = tailfr_numnew

             if (xferfrac_vol >= xferfrac_max_in) then
                xferfrac_vol = 1.0_r8
                xferfrac_num = 1.0_r8
             else
                xferfrac_vol = min( xferfrac_vol, xferfrac_max_in )
                xferfrac_num = max( xferfrac_num, xferfrac_vol )
                xferfrac_num = min( xferfrac_max_in, xferfrac_num )
             end if

             if (ixferable_all_in <= 0) then
                dryvol_t_new_xfab = max( 0.0_r8, (dryvol_a_xfab(i,k) + dryvol_c_xfab(i,k)) )
                dryvol_xferamt = xferfrac_vol*dryvol_t_new
                if (dryvol_t_new_xfab >= 0.999999_r8*dryvol_xferamt) then
                   xferfrac_vol = min( 1.0_r8, (dryvol_xferamt/dryvol_t_new_xfab) )
                else if (dryvol_t_new_xfab >= 1.0e-7_r8*dryvol_xferamt) then
                   xferfrac_vol = 1.0_r8
                   xferfrac_num = xferfrac_num*(dryvol_t_new_xfab/dryvol_xferamt)
                else
                   cycle
                end if
             end if
          end if grow_shrink_conditional1

          xferfrac_vol_out(i,k) = xferfrac_vol
          xferfrac_num_out(i,k) = xferfrac_num
       end do
    end do

  end subroutine modal_aero_rename_acc_crs_xferfracs_native

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_tendencies( ncol, loffset, deltat, deltatinv, &
       is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, pdel, q, qqcw, dqdt, dqqcwdt, &
       qsrflx, qqcwsrflx, xferfrac_vol_ik, xferfrac_num_ik, nspecfrm_renamexf_in, &
       lspecfrma_renamexf_in, lspecfrmc_renamexf_in, lspectooa_renamexf_in, lspectooc_renamexf_in, &
       l_dqdt_rnpos, dqdt_rnpos )

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use physconst, only: gravit

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: deltatinv
    logical, intent(in) :: is_dorename_atik
    logical, intent(in) :: dorename_atik(ncol,pver)
    integer, intent(in) :: jsrflx_rename
    integer, intent(in) :: nsrflx
    real(r8), target, intent(in) :: pdel(pcols,pver)
    real(r8), target, intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), target, intent(inout) :: dqdt(ncol,pver,pcnstxx)
    real(r8), target, intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), target, intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
    real(r8), target, intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
    real(r8), target, intent(in) :: xferfrac_vol_ik(ncol,pver)
    real(r8), target, intent(in) :: xferfrac_num_ik(ncol,pver)
    integer, target, intent(in) :: nspecfrm_renamexf_in
    integer, target, intent(in) :: lspecfrma_renamexf_in(maxspec_renamexf)
    integer, target, intent(in) :: lspecfrmc_renamexf_in(maxspec_renamexf)
    integer, target, intent(in) :: lspectooa_renamexf_in(maxspec_renamexf)
    integer, target, intent(in) :: lspectooc_renamexf_in(maxspec_renamexf)
    logical, intent(in) :: l_dqdt_rnpos
    real(r8), optional, target, intent(inout) :: dqdt_rnpos(ncol,pver,pcnstxx)
    integer(c_int64_t), target :: dorename_atik_c(ncol,pver)
    integer(c_int64_t), target :: lspecfrma_renamexf_c(maxspec_renamexf)
    integer(c_int64_t), target :: lspecfrmc_renamexf_c(maxspec_renamexf)
    integer(c_int64_t), target :: lspectooa_renamexf_c(maxspec_renamexf)
    integer(c_int64_t), target :: lspectooc_renamexf_c(maxspec_renamexf)
    integer(c_int64_t) :: is_dorename_atik_c
    integer(c_int64_t) :: l_dqdt_rnpos_c
    integer(c_int64_t) :: nspecfrm_renamexf_c
    real(r8), target :: dqdt_rnpos_dummy(1,1,1)
    type(c_ptr) :: dqdt_rnpos_p
    integer :: i, k, iq

    interface
       subroutine modal_aero_rename_acc_crs_tendencies_codon( &
            ncol_c, pcols_c, pver_c, pcnstxx_c, maxspec_renamexf_c, loffset_c, &
            nspecfrm_renamexf_c, jsrflx_rename_c, nsrflx_c, is_dorename_atik_c, l_dqdt_rnpos_c, &
            deltat_c, deltatinv_c, gravit_c, pdel_p, dorename_atik_p, q_p, qqcw_p, dqdt_p, dqqcwdt_p, &
            qsrflx_p, qqcwsrflx_p, xferfrac_vol_p, xferfrac_num_p, lspecfrma_renamexf_p, &
            lspecfrmc_renamexf_p, lspectooa_renamexf_p, lspectooc_renamexf_p, dqdt_rnpos_p ) &
            bind(c, name="modal_aero_rename_acc_crs_tendencies_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnstxx_c, maxspec_renamexf_c
         integer(c_int64_t), value :: loffset_c, nspecfrm_renamexf_c, jsrflx_rename_c, nsrflx_c
         integer(c_int64_t), value :: is_dorename_atik_c, l_dqdt_rnpos_c
         real(c_double), value :: deltat_c, deltatinv_c, gravit_c
         type(c_ptr), value :: pdel_p, dorename_atik_p, q_p, qqcw_p, dqdt_p, dqqcwdt_p
         type(c_ptr), value :: qsrflx_p, qqcwsrflx_p, xferfrac_vol_p, xferfrac_num_p
         type(c_ptr), value :: lspecfrma_renamexf_p, lspecfrmc_renamexf_p
         type(c_ptr), value :: lspectooa_renamexf_p, lspectooc_renamexf_p, dqdt_rnpos_p
       end subroutine modal_aero_rename_acc_crs_tendencies_codon
    end interface

    call modal_aero_rename_acc_crs_tendencies_select_impl()

    if (.not. modal_aero_rename_acc_crs_tendencies_use_native_impl) then
       if (is_dorename_atik) then
          is_dorename_atik_c = 1_c_int64_t
          do k = 1, pver
             do i = 1, ncol
                if (dorename_atik(i,k)) then
                   dorename_atik_c(i,k) = 1_c_int64_t
                else
                   dorename_atik_c(i,k) = 0_c_int64_t
                end if
             end do
          end do
       else
          is_dorename_atik_c = 0_c_int64_t
          dorename_atik_c(:,:) = 0_c_int64_t
       end if

       if (l_dqdt_rnpos) then
          l_dqdt_rnpos_c = 1_c_int64_t
       else
          l_dqdt_rnpos_c = 0_c_int64_t
       end if

       nspecfrm_renamexf_c = int(nspecfrm_renamexf_in, c_int64_t)
       do iq = 1, maxspec_renamexf
          lspecfrma_renamexf_c(iq) = int(lspecfrma_renamexf_in(iq), c_int64_t)
          lspecfrmc_renamexf_c(iq) = int(lspecfrmc_renamexf_in(iq), c_int64_t)
          lspectooa_renamexf_c(iq) = int(lspectooa_renamexf_in(iq), c_int64_t)
          lspectooc_renamexf_c(iq) = int(lspectooc_renamexf_in(iq), c_int64_t)
       end do

       if (present(dqdt_rnpos)) then
          dqdt_rnpos_p = c_loc(dqdt_rnpos(1,1,1))
       else
          dqdt_rnpos_dummy(1,1,1) = 0.0_r8
          dqdt_rnpos_p = c_loc(dqdt_rnpos_dummy(1,1,1))
       end if

       call modal_aero_rename_acc_crs_tendencies_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnstxx, c_int64_t), &
            int(maxspec_renamexf, c_int64_t), int(loffset, c_int64_t), nspecfrm_renamexf_c, int(jsrflx_rename, c_int64_t), &
            int(nsrflx, c_int64_t), is_dorename_atik_c, l_dqdt_rnpos_c, real(deltat, c_double), real(deltatinv, c_double), &
            real(gravit, c_double), c_loc(pdel(1,1)), c_loc(dorename_atik_c(1,1)), c_loc(q(1,1,1)), c_loc(qqcw(1,1,1)), &
            c_loc(dqdt(1,1,1)), c_loc(dqqcwdt(1,1,1)), c_loc(qsrflx(1,1,1)), c_loc(qqcwsrflx(1,1,1)), c_loc(xferfrac_vol_ik(1,1)), &
            c_loc(xferfrac_num_ik(1,1)), c_loc(lspecfrma_renamexf_c(1)), c_loc(lspecfrmc_renamexf_c(1)), c_loc(lspectooa_renamexf_c(1)), &
            c_loc(lspectooc_renamexf_c(1)), dqdt_rnpos_p &
       )
       return
    end if

    call modal_aero_rename_acc_crs_tendencies_native( ncol, loffset, deltat, deltatinv, &
         is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, pdel, q, qqcw, dqdt, dqqcwdt, &
         qsrflx, qqcwsrflx, xferfrac_vol_ik, xferfrac_num_ik, nspecfrm_renamexf_in, lspecfrma_renamexf_in, &
         lspecfrmc_renamexf_in, lspectooa_renamexf_in, lspectooc_renamexf_in, l_dqdt_rnpos, dqdt_rnpos )

  end subroutine modal_aero_rename_acc_crs_tendencies

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_tendencies_native( ncol, loffset, deltat, deltatinv, &
       is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, pdel, q, qqcw, dqdt, dqqcwdt, &
       qsrflx, qqcwsrflx, xferfrac_vol_ik, xferfrac_num_ik, nspecfrm_renamexf_in, &
       lspecfrma_renamexf_in, lspecfrmc_renamexf_in, lspectooa_renamexf_in, lspectooc_renamexf_in, &
       l_dqdt_rnpos, dqdt_rnpos )

    use physconst, only: gravit

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: deltatinv
    logical, intent(in) :: is_dorename_atik
    logical, intent(in) :: dorename_atik(ncol,pver)
    integer, intent(in) :: jsrflx_rename
    integer, intent(in) :: nsrflx
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: dqdt(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
    real(r8), intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
    real(r8), intent(in) :: xferfrac_vol_ik(ncol,pver)
    real(r8), intent(in) :: xferfrac_num_ik(ncol,pver)
    integer, intent(in) :: nspecfrm_renamexf_in
    integer, intent(in) :: lspecfrma_renamexf_in(maxspec_renamexf)
    integer, intent(in) :: lspecfrmc_renamexf_in(maxspec_renamexf)
    integer, intent(in) :: lspectooa_renamexf_in(maxspec_renamexf)
    integer, intent(in) :: lspectooc_renamexf_in(maxspec_renamexf)
    logical, intent(in) :: l_dqdt_rnpos
    real(r8), optional, intent(inout) :: dqdt_rnpos(ncol,pver,pcnstxx)

    integer :: i, j, k, iq
    integer :: lsfrma, lsfrmc, lstooa, lstooc
    real(r8) :: pdel_fac, xfercoef, xfertend, xferfrac_vol, xferfrac_num

mainloop1_k:  do k = 1, pver
mainloop1_i:  do i = 1, ncol

       if (is_dorename_atik) then
          if (.not. dorename_atik(i,k)) cycle mainloop1_i
       end if

       xferfrac_vol = xferfrac_vol_ik(i,k)
       xferfrac_num = xferfrac_num_ik(i,k)
       if (xferfrac_vol .le. 0.0_r8) cycle mainloop1_i

       pdel_fac = pdel(i,k)/gravit
       j = jsrflx_rename
       do iq = 1, nspecfrm_renamexf_in
          xfercoef = xferfrac_vol*deltatinv
          if (iq .eq. 1) xfercoef = xferfrac_num*deltatinv

          lsfrma = lspecfrma_renamexf_in(iq)-loffset
          lsfrmc = lspecfrmc_renamexf_in(iq)-loffset
          lstooa = lspectooa_renamexf_in(iq)-loffset
          lstooc = lspectooc_renamexf_in(iq)-loffset

          if (lsfrma .gt. 0) then
             xfertend = xfercoef*max( 0.0_r8, (q(i,k,lsfrma)+dqdt(i,k,lsfrma)*deltat) )
             dqdt(i,k,lsfrma) = dqdt(i,k,lsfrma) - xfertend
             qsrflx(i,lsfrma,j) = qsrflx(i,lsfrma,j) - xfertend*pdel_fac
             if (lstooa .gt. 0) then
                dqdt(i,k,lstooa) = dqdt(i,k,lstooa) + xfertend
                qsrflx(i,lstooa,j) = qsrflx(i,lstooa,j) + xfertend*pdel_fac
                if ( l_dqdt_rnpos .and. present(dqdt_rnpos) ) then
                   dqdt_rnpos(i,k,lstooa) = dqdt_rnpos(i,k,lstooa) + xfertend
                end if
             end if
          end if

          if (lsfrmc .gt. 0) then
             xfertend = xfercoef*max( 0.0_r8, (qqcw(i,k,lsfrmc)+dqqcwdt(i,k,lsfrmc)*deltat) )
             dqqcwdt(i,k,lsfrmc) = dqqcwdt(i,k,lsfrmc) - xfertend
             qqcwsrflx(i,lsfrmc,j) = qqcwsrflx(i,lsfrmc,j) - xfertend*pdel_fac
             if (lstooc .gt. 0) then
                dqqcwdt(i,k,lstooc) = dqqcwdt(i,k,lstooc) + xfertend
                qqcwsrflx(i,lstooc,j) = qqcwsrflx(i,lstooc,j) + xfertend*pdel_fac
             end if
          end if
       end do

    end do mainloop1_i
    end do mainloop1_k

  end subroutine modal_aero_rename_acc_crs_tendencies_native


!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_pair( ncol, loffset, deltat, deltatinv, troplev, pdel, q, qqcw, &
       dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, &
       qsrflx, qqcwsrflx, ixferable_all_in, nspec_mfrm_in, lspectype_mfrm_in, specmw_amode_in, specdens_amode_in, &
       lmassptr_mfrm_in, lmassptrcw_mfrm_in, ixferable_a_in, ixferable_c_in, mfrm_in, numptr_amode_mfrm_in, &
       numptrcw_amode_mfrm_in, dgnum_amode_mfrm_in, factoraa_in, factoryy_in, dryvol_smallest_in, v2nlorlx_in, &
       v2nhirlx_in, factor_3alnsg2_in, dp_cut_in, lndp_cut_in, dp_belowcut_in, dp_xfernone_thresh_in, &
       dp_xferall_thresh_in, flagaa_shrink_in, igrow_shrink_in, method_optbb_in, onethird_in, xferfrac_max_in, &
       nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, lspectooa_renamexf_in, lspectooc_renamexf_in, &
       dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab, xferfrac_vol_ik, xferfrac_num_ik, &
       l_dqdt_rnpos, dqdt_rnpos )

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use physconst, only: gravit

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: deltatinv
    integer, target, intent(in) :: troplev(pcols)
    real(r8), target, intent(in) :: pdel(pcols,pver)
    real(r8), target, intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), target, intent(inout) :: dqdt(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqdt_other(ncol,pver,pcnstxx)
    real(r8), target, intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), target, intent(in) :: dqqcwdt_other(ncol,pver,pcnstxx)
    logical, intent(in) :: is_dorename_atik
    logical, target, intent(in) :: dorename_atik(ncol,pver)
    integer, intent(in) :: jsrflx_rename
    integer, intent(in) :: nsrflx
    real(r8), target, intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
    real(r8), target, intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
    integer, intent(in) :: ixferable_all_in
    integer, intent(in) :: nspec_mfrm_in
    integer, target, intent(in) :: lspectype_mfrm_in(maxspec_renamexf)
    real(r8), target, intent(in) :: specmw_amode_in(*)
    real(r8), target, intent(in) :: specdens_amode_in(*)
    integer, target, intent(in) :: lmassptr_mfrm_in(maxspec_renamexf)
    integer, target, intent(in) :: lmassptrcw_mfrm_in(maxspec_renamexf)
    integer, target, intent(in) :: ixferable_a_in(maxspec_renamexf)
    integer, target, intent(in) :: ixferable_c_in(maxspec_renamexf)
    integer, intent(in) :: mfrm_in
    integer, intent(in) :: numptr_amode_mfrm_in
    integer, intent(in) :: numptrcw_amode_mfrm_in
    real(r8), intent(in) :: dgnum_amode_mfrm_in
    real(r8), intent(in) :: factoraa_in
    real(r8), intent(in) :: factoryy_in
    real(r8), intent(in) :: dryvol_smallest_in
    real(r8), intent(in) :: v2nlorlx_in
    real(r8), intent(in) :: v2nhirlx_in
    real(r8), intent(in) :: factor_3alnsg2_in
    real(r8), intent(in) :: dp_cut_in
    real(r8), intent(in) :: lndp_cut_in
    real(r8), intent(in) :: dp_belowcut_in
    real(r8), intent(in) :: dp_xfernone_thresh_in
    real(r8), intent(in) :: dp_xferall_thresh_in
    logical, intent(in) :: flagaa_shrink_in
    integer, intent(in) :: igrow_shrink_in
    integer, intent(in) :: method_optbb_in
    real(r8), intent(in) :: onethird_in
    real(r8), intent(in) :: xferfrac_max_in
    integer, intent(in) :: nspecfrm_renamexf_in
    integer, target, intent(in) :: lspecfrma_renamexf_in(maxspec_renamexf)
    integer, target, intent(in) :: lspecfrmc_renamexf_in(maxspec_renamexf)
    integer, target, intent(in) :: lspectooa_renamexf_in(maxspec_renamexf)
    integer, target, intent(in) :: lspectooc_renamexf_in(maxspec_renamexf)
    real(r8), target, intent(out) :: dryvol_a(ncol,pver)
    real(r8), target, intent(out) :: dryvol_c(ncol,pver)
    real(r8), target, intent(out) :: deldryvol_a(ncol,pver)
    real(r8), target, intent(out) :: deldryvol_c(ncol,pver)
    real(r8), target, intent(out) :: dryvol_a_xfab(ncol,pver)
    real(r8), target, intent(out) :: dryvol_c_xfab(ncol,pver)
    real(r8), target, intent(out) :: xferfrac_vol_ik(ncol,pver)
    real(r8), target, intent(out) :: xferfrac_num_ik(ncol,pver)
    logical, intent(in) :: l_dqdt_rnpos
    real(r8), optional, target, intent(inout) :: dqdt_rnpos(ncol,pver,pcnstxx)
    integer(c_int64_t), target :: troplev_c(pcols)
    integer(c_int64_t), target :: dorename_atik_c(ncol,pver)
    integer(c_int64_t), target :: lspectype_mfrm_c(maxspec_renamexf)
    integer(c_int64_t), target :: lmassptr_mfrm_c(maxspec_renamexf)
    integer(c_int64_t), target :: lmassptrcw_mfrm_c(maxspec_renamexf)
    integer(c_int64_t), target :: ixferable_a_c(maxspec_renamexf)
    integer(c_int64_t), target :: ixferable_c_c(maxspec_renamexf)
    integer(c_int64_t), target :: lspecfrma_renamexf_c(maxspec_renamexf)
    integer(c_int64_t), target :: lspecfrmc_renamexf_c(maxspec_renamexf)
    integer(c_int64_t), target :: lspectooa_renamexf_c(maxspec_renamexf)
    integer(c_int64_t), target :: lspectooc_renamexf_c(maxspec_renamexf)
    integer(c_int64_t) :: is_dorename_atik_c
    integer(c_int64_t) :: l_dqdt_rnpos_c
    integer(c_int64_t) :: flagaa_shrink_c
    real(r8), target :: dqdt_rnpos_dummy(1,1,1)
    type(c_ptr) :: dqdt_rnpos_p
    integer :: i, k, iq

    interface
       subroutine modal_aero_rename_acc_crs_pair_codon( &
            ncol_c, pcols_c, pver_c, pcnstxx_c, maxspec_renamexf_c, loffset_c, is_dorename_atik_c, &
            l_dqdt_rnpos_c, jsrflx_rename_c, nsrflx_c, ixferable_all_c, nspec_mfrm_c, mfrm_c, &
            numptr_amode_mfrm_c, numptrcw_amode_mfrm_c, igrow_shrink_c, method_optbb_c, flagaa_shrink_c, &
            nspecfrm_renamexf_c, deltat_c, deltatinv_c, gravit_c, dgnum_amode_mfrm_c, factoraa_c, factoryy_c, &
            dryvol_smallest_c, v2nlorlx_c, v2nhirlx_c, factor_3alnsg2_c, dp_cut_c, lndp_cut_c, dp_belowcut_c, &
            dp_xfernone_thresh_c, dp_xferall_thresh_c, onethird_c, xferfrac_max_c, troplev_p, pdel_p, &
            dorename_atik_p, q_p, qqcw_p, dqdt_p, dqdt_other_p, dqqcwdt_p, dqqcwdt_other_p, qsrflx_p, qqcwsrflx_p, &
            lspectype_mfrm_p, specmw_amode_p, specdens_amode_p, lmassptr_mfrm_p, lmassptrcw_mfrm_p, ixferable_a_p, &
            ixferable_c_p, lspecfrma_renamexf_p, lspecfrmc_renamexf_p, lspectooa_renamexf_p, lspectooc_renamexf_p, &
            dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p, dryvol_a_xfab_p, dryvol_c_xfab_p, &
            xferfrac_vol_p, xferfrac_num_p, dqdt_rnpos_p ) bind(c, name="modal_aero_rename_acc_crs_pair_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnstxx_c, maxspec_renamexf_c
         integer(c_int64_t), value :: loffset_c, is_dorename_atik_c, l_dqdt_rnpos_c
         integer(c_int64_t), value :: jsrflx_rename_c, nsrflx_c, ixferable_all_c, nspec_mfrm_c, mfrm_c
         integer(c_int64_t), value :: numptr_amode_mfrm_c, numptrcw_amode_mfrm_c
         integer(c_int64_t), value :: igrow_shrink_c, method_optbb_c, flagaa_shrink_c, nspecfrm_renamexf_c
         real(c_double), value :: deltat_c, deltatinv_c, gravit_c, dgnum_amode_mfrm_c, factoraa_c
         real(c_double), value :: factoryy_c, dryvol_smallest_c, v2nlorlx_c, v2nhirlx_c
         real(c_double), value :: factor_3alnsg2_c, dp_cut_c, lndp_cut_c, dp_belowcut_c
         real(c_double), value :: dp_xfernone_thresh_c, dp_xferall_thresh_c, onethird_c, xferfrac_max_c
         type(c_ptr), value :: troplev_p, pdel_p, dorename_atik_p, q_p, qqcw_p, dqdt_p, dqdt_other_p
         type(c_ptr), value :: dqqcwdt_p, dqqcwdt_other_p, qsrflx_p, qqcwsrflx_p
         type(c_ptr), value :: lspectype_mfrm_p, specmw_amode_p, specdens_amode_p, lmassptr_mfrm_p
         type(c_ptr), value :: lmassptrcw_mfrm_p, ixferable_a_p, ixferable_c_p
         type(c_ptr), value :: lspecfrma_renamexf_p, lspecfrmc_renamexf_p, lspectooa_renamexf_p
         type(c_ptr), value :: lspectooc_renamexf_p, dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p
         type(c_ptr), value :: dryvol_a_xfab_p, dryvol_c_xfab_p, xferfrac_vol_p, xferfrac_num_p, dqdt_rnpos_p
       end subroutine modal_aero_rename_acc_crs_pair_codon
    end interface

    call modal_aero_rename_acc_crs_pair_select_impl()

    if (.not. modal_aero_rename_acc_crs_pair_use_native_impl) then
       if (is_dorename_atik) then
          is_dorename_atik_c = 1_c_int64_t
          do k = 1, pver
             do i = 1, ncol
                if (dorename_atik(i,k)) then
                   dorename_atik_c(i,k) = 1_c_int64_t
                else
                   dorename_atik_c(i,k) = 0_c_int64_t
                end if
             end do
          end do
       else
          is_dorename_atik_c = 0_c_int64_t
          dorename_atik_c(:,:) = 0_c_int64_t
       end if

       if (l_dqdt_rnpos) then
          l_dqdt_rnpos_c = 1_c_int64_t
       else
          l_dqdt_rnpos_c = 0_c_int64_t
       end if

       if (flagaa_shrink_in) then
          flagaa_shrink_c = 1_c_int64_t
       else
          flagaa_shrink_c = 0_c_int64_t
       end if

       do iq = 1, maxspec_renamexf
          lspectype_mfrm_c(iq) = int(lspectype_mfrm_in(iq), c_int64_t)
          lmassptr_mfrm_c(iq) = int(lmassptr_mfrm_in(iq), c_int64_t)
          lmassptrcw_mfrm_c(iq) = int(lmassptrcw_mfrm_in(iq), c_int64_t)
          ixferable_a_c(iq) = int(ixferable_a_in(iq), c_int64_t)
          ixferable_c_c(iq) = int(ixferable_c_in(iq), c_int64_t)
          lspecfrma_renamexf_c(iq) = int(lspecfrma_renamexf_in(iq), c_int64_t)
          lspecfrmc_renamexf_c(iq) = int(lspecfrmc_renamexf_in(iq), c_int64_t)
          lspectooa_renamexf_c(iq) = int(lspectooa_renamexf_in(iq), c_int64_t)
          lspectooc_renamexf_c(iq) = int(lspectooc_renamexf_in(iq), c_int64_t)
       end do
       do i = 1, pcols
          troplev_c(i) = int(troplev(i), c_int64_t)
       end do

       if (present(dqdt_rnpos)) then
          dqdt_rnpos_p = c_loc(dqdt_rnpos(1,1,1))
       else
          dqdt_rnpos_dummy(1,1,1) = 0.0_r8
          dqdt_rnpos_p = c_loc(dqdt_rnpos_dummy(1,1,1))
       end if

       call modal_aero_rename_acc_crs_pair_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnstxx, c_int64_t), &
            int(maxspec_renamexf, c_int64_t), int(loffset, c_int64_t), is_dorename_atik_c, l_dqdt_rnpos_c, &
            int(jsrflx_rename, c_int64_t), int(nsrflx, c_int64_t), int(ixferable_all_in, c_int64_t), int(nspec_mfrm_in, c_int64_t), &
            int(mfrm_in, c_int64_t), int(numptr_amode_mfrm_in, c_int64_t), int(numptrcw_amode_mfrm_in, c_int64_t), &
            int(igrow_shrink_in, c_int64_t), int(method_optbb_in, c_int64_t), flagaa_shrink_c, int(nspecfrm_renamexf_in, c_int64_t), &
            real(deltat, c_double), real(deltatinv, c_double), real(gravit, c_double), real(dgnum_amode_mfrm_in, c_double), &
            real(factoraa_in, c_double), real(factoryy_in, c_double), real(dryvol_smallest_in, c_double), real(v2nlorlx_in, c_double), &
            real(v2nhirlx_in, c_double), real(factor_3alnsg2_in, c_double), real(dp_cut_in, c_double), real(lndp_cut_in, c_double), &
            real(dp_belowcut_in, c_double), real(dp_xfernone_thresh_in, c_double), real(dp_xferall_thresh_in, c_double), &
            real(onethird_in, c_double), real(xferfrac_max_in, c_double), c_loc(troplev_c(1)), c_loc(pdel(1,1)), c_loc(dorename_atik_c(1,1)), &
            c_loc(q(1,1,1)), c_loc(qqcw(1,1,1)), c_loc(dqdt(1,1,1)), c_loc(dqdt_other(1,1,1)), c_loc(dqqcwdt(1,1,1)), &
            c_loc(dqqcwdt_other(1,1,1)), c_loc(qsrflx(1,1,1)), c_loc(qqcwsrflx(1,1,1)), c_loc(lspectype_mfrm_c(1)), &
            c_loc(specmw_amode_in(1)), c_loc(specdens_amode_in(1)), c_loc(lmassptr_mfrm_c(1)), c_loc(lmassptrcw_mfrm_c(1)), &
            c_loc(ixferable_a_c(1)), c_loc(ixferable_c_c(1)), c_loc(lspecfrma_renamexf_c(1)), c_loc(lspecfrmc_renamexf_c(1)), &
            c_loc(lspectooa_renamexf_c(1)), c_loc(lspectooc_renamexf_c(1)), c_loc(dryvol_a(1,1)), c_loc(dryvol_c(1,1)), &
            c_loc(deldryvol_a(1,1)), c_loc(deldryvol_c(1,1)), c_loc(dryvol_a_xfab(1,1)), c_loc(dryvol_c_xfab(1,1)), &
            c_loc(xferfrac_vol_ik(1,1)), c_loc(xferfrac_num_ik(1,1)), dqdt_rnpos_p &
       )
       return
    end if

    call modal_aero_rename_acc_crs_pair_native( ncol, loffset, deltat, deltatinv, troplev, pdel, q, qqcw, &
         dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, &
         qsrflx, qqcwsrflx, ixferable_all_in, nspec_mfrm_in, lspectype_mfrm_in, specmw_amode_in, specdens_amode_in, &
         lmassptr_mfrm_in, lmassptrcw_mfrm_in, ixferable_a_in, ixferable_c_in, mfrm_in, numptr_amode_mfrm_in, &
         numptrcw_amode_mfrm_in, dgnum_amode_mfrm_in, factoraa_in, factoryy_in, dryvol_smallest_in, v2nlorlx_in, &
         v2nhirlx_in, factor_3alnsg2_in, dp_cut_in, lndp_cut_in, dp_belowcut_in, dp_xfernone_thresh_in, &
         dp_xferall_thresh_in, flagaa_shrink_in, igrow_shrink_in, method_optbb_in, onethird_in, xferfrac_max_in, &
         nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, lspectooa_renamexf_in, lspectooc_renamexf_in, &
         dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab, xferfrac_vol_ik, xferfrac_num_ik, &
         l_dqdt_rnpos, dqdt_rnpos )

  end subroutine modal_aero_rename_acc_crs_pair

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_acc_crs_pair_native( ncol, loffset, deltat, deltatinv, troplev, pdel, q, qqcw, &
       dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, is_dorename_atik, dorename_atik, jsrflx_rename, nsrflx, &
       qsrflx, qqcwsrflx, ixferable_all_in, nspec_mfrm_in, lspectype_mfrm_in, specmw_amode_in, specdens_amode_in, &
       lmassptr_mfrm_in, lmassptrcw_mfrm_in, ixferable_a_in, ixferable_c_in, mfrm_in, numptr_amode_mfrm_in, &
       numptrcw_amode_mfrm_in, dgnum_amode_mfrm_in, factoraa_in, factoryy_in, dryvol_smallest_in, v2nlorlx_in, &
       v2nhirlx_in, factor_3alnsg2_in, dp_cut_in, lndp_cut_in, dp_belowcut_in, dp_xfernone_thresh_in, &
       dp_xferall_thresh_in, flagaa_shrink_in, igrow_shrink_in, method_optbb_in, onethird_in, xferfrac_max_in, &
       nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, lspectooa_renamexf_in, lspectooc_renamexf_in, &
       dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab, xferfrac_vol_ik, xferfrac_num_ik, &
       l_dqdt_rnpos, dqdt_rnpos )

    implicit none

    integer, intent(in) :: ncol
    integer, intent(in) :: loffset
    real(r8), intent(in) :: deltat
    real(r8), intent(in) :: deltatinv
    integer, intent(in) :: troplev(pcols)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: q(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: dqdt(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqdt_other(ncol,pver,pcnstxx)
    real(r8), intent(in) :: qqcw(ncol,pver,pcnstxx)
    real(r8), intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
    real(r8), intent(in) :: dqqcwdt_other(ncol,pver,pcnstxx)
    logical, intent(in) :: is_dorename_atik
    logical, intent(in) :: dorename_atik(ncol,pver)
    integer, intent(in) :: jsrflx_rename
    integer, intent(in) :: nsrflx
    real(r8), intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
    real(r8), intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
    integer, intent(in) :: ixferable_all_in
    integer, intent(in) :: nspec_mfrm_in
    integer, intent(in) :: lspectype_mfrm_in(maxspec_renamexf)
    real(r8), intent(in) :: specmw_amode_in(*)
    real(r8), intent(in) :: specdens_amode_in(*)
    integer, intent(in) :: lmassptr_mfrm_in(maxspec_renamexf)
    integer, intent(in) :: lmassptrcw_mfrm_in(maxspec_renamexf)
    integer, intent(in) :: ixferable_a_in(maxspec_renamexf)
    integer, intent(in) :: ixferable_c_in(maxspec_renamexf)
    integer, intent(in) :: mfrm_in
    integer, intent(in) :: numptr_amode_mfrm_in
    integer, intent(in) :: numptrcw_amode_mfrm_in
    real(r8), intent(in) :: dgnum_amode_mfrm_in
    real(r8), intent(in) :: factoraa_in
    real(r8), intent(in) :: factoryy_in
    real(r8), intent(in) :: dryvol_smallest_in
    real(r8), intent(in) :: v2nlorlx_in
    real(r8), intent(in) :: v2nhirlx_in
    real(r8), intent(in) :: factor_3alnsg2_in
    real(r8), intent(in) :: dp_cut_in
    real(r8), intent(in) :: lndp_cut_in
    real(r8), intent(in) :: dp_belowcut_in
    real(r8), intent(in) :: dp_xfernone_thresh_in
    real(r8), intent(in) :: dp_xferall_thresh_in
    logical, intent(in) :: flagaa_shrink_in
    integer, intent(in) :: igrow_shrink_in
    integer, intent(in) :: method_optbb_in
    real(r8), intent(in) :: onethird_in
    real(r8), intent(in) :: xferfrac_max_in
    integer, intent(in) :: nspecfrm_renamexf_in
    integer, intent(in) :: lspecfrma_renamexf_in(maxspec_renamexf)
    integer, intent(in) :: lspecfrmc_renamexf_in(maxspec_renamexf)
    integer, intent(in) :: lspectooa_renamexf_in(maxspec_renamexf)
    integer, intent(in) :: lspectooc_renamexf_in(maxspec_renamexf)
    real(r8), intent(out) :: dryvol_a(ncol,pver)
    real(r8), intent(out) :: dryvol_c(ncol,pver)
    real(r8), intent(out) :: deldryvol_a(ncol,pver)
    real(r8), intent(out) :: deldryvol_c(ncol,pver)
    real(r8), intent(out) :: dryvol_a_xfab(ncol,pver)
    real(r8), intent(out) :: dryvol_c_xfab(ncol,pver)
    real(r8), intent(out) :: xferfrac_vol_ik(ncol,pver)
    real(r8), intent(out) :: xferfrac_num_ik(ncol,pver)
    logical, intent(in) :: l_dqdt_rnpos
    real(r8), optional, intent(inout) :: dqdt_rnpos(ncol,pver,pcnstxx)

    call modal_aero_rename_acc_crs_dryvols( ncol, loffset, deltat, ixferable_all_in, q, qqcw, dqdt, dqdt_other, &
         dqqcwdt, dqqcwdt_other, nspec_mfrm_in, lspectype_mfrm_in, specmw_amode_in, specdens_amode_in, lmassptr_mfrm_in, &
         lmassptrcw_mfrm_in, ixferable_a_in, ixferable_c_in, dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, &
         dryvol_a_xfab, dryvol_c_xfab )

    call modal_aero_rename_acc_crs_xferfracs( ncol, loffset, troplev, q, qqcw, dryvol_a, dryvol_c, deldryvol_a, &
         deldryvol_c, dryvol_a_xfab, dryvol_c_xfab, mfrm_in, numptr_amode_mfrm_in, numptrcw_amode_mfrm_in, &
         dgnum_amode_mfrm_in, factoraa_in, factoryy_in, dryvol_smallest_in, v2nlorlx_in, v2nhirlx_in, factor_3alnsg2_in, &
         dp_cut_in, lndp_cut_in, dp_belowcut_in, dp_xfernone_thresh_in, dp_xferall_thresh_in, flagaa_shrink_in, &
         igrow_shrink_in, ixferable_all_in, method_optbb_in, onethird_in, xferfrac_max_in, xferfrac_vol_ik, xferfrac_num_ik )

    if (l_dqdt_rnpos) then
       call modal_aero_rename_acc_crs_tendencies( ncol, loffset, deltat, deltatinv, is_dorename_atik, dorename_atik, &
            jsrflx_rename, nsrflx, pdel, q, qqcw, dqdt, dqqcwdt, qsrflx, qqcwsrflx, xferfrac_vol_ik, xferfrac_num_ik, &
            nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, lspectooa_renamexf_in, lspectooc_renamexf_in, &
            l_dqdt_rnpos, dqdt_rnpos )
    else
       call modal_aero_rename_acc_crs_tendencies( ncol, loffset, deltat, deltatinv, is_dorename_atik, dorename_atik, &
            jsrflx_rename, nsrflx, pdel, q, qqcw, dqdt, dqqcwdt, qsrflx, qqcwsrflx, xferfrac_vol_ik, xferfrac_num_ik, &
            nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, lspectooa_renamexf_in, lspectooc_renamexf_in, &
            l_dqdt_rnpos )
    end if

  end subroutine modal_aero_rename_acc_crs_pair_native


!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_set_dotend_flags( loffset, npair_renamexf_in, &
       nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, &
       lspectooa_renamexf_in, lspectooc_renamexf_in, dotendrn, dotendqqcwrn )

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr

    implicit none

    integer, intent(in) :: loffset
    integer, intent(in) :: npair_renamexf_in
    integer, target, intent(in) :: nspecfrm_renamexf_in(maxpair_renamexf)
    integer, target, intent(in) :: lspecfrma_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, target, intent(in) :: lspecfrmc_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, target, intent(in) :: lspectooa_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, target, intent(in) :: lspectooc_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    logical, intent(out) :: dotendrn(pcnstxx)
    logical, intent(out) :: dotendqqcwrn(pcnstxx)
    integer(c_int64_t), target :: nspecfrm_renamexf_c(maxpair_renamexf)
    integer(c_int64_t), target :: lspecfrma_renamexf_c(maxspec_renamexf,maxpair_renamexf)
    integer(c_int64_t), target :: lspecfrmc_renamexf_c(maxspec_renamexf,maxpair_renamexf)
    integer(c_int64_t), target :: lspectooa_renamexf_c(maxspec_renamexf,maxpair_renamexf)
    integer(c_int64_t), target :: lspectooc_renamexf_c(maxspec_renamexf,maxpair_renamexf)
    integer(c_int64_t), target :: dotendrn_c(pcnstxx)
    integer(c_int64_t), target :: dotendqqcwrn_c(pcnstxx)
    integer :: ipair, iq, l

    interface
       subroutine modal_aero_rename_set_dotend_flags_codon( &
            pcnstxx_c, maxpair_renamexf_c, maxspec_renamexf_c, loffset_c, npair_renamexf_c, &
            nspecfrm_renamexf_p, lspecfrma_renamexf_p, lspecfrmc_renamexf_p, lspectooa_renamexf_p, &
            lspectooc_renamexf_p, dotendrn_p, dotendqqcwrn_p ) &
            bind(c, name="modal_aero_rename_set_dotend_flags_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: pcnstxx_c, maxpair_renamexf_c, maxspec_renamexf_c
         integer(c_int64_t), value :: loffset_c, npair_renamexf_c
         type(c_ptr), value :: nspecfrm_renamexf_p, lspecfrma_renamexf_p, lspecfrmc_renamexf_p
         type(c_ptr), value :: lspectooa_renamexf_p, lspectooc_renamexf_p, dotendrn_p, dotendqqcwrn_p
       end subroutine modal_aero_rename_set_dotend_flags_codon
    end interface

    call modal_aero_rename_set_dotend_flags_select_impl()

    if (.not. modal_aero_rename_set_dotend_flags_use_native_impl) then
       do ipair = 1, maxpair_renamexf
          nspecfrm_renamexf_c(ipair) = int(nspecfrm_renamexf_in(ipair), c_int64_t)
          do iq = 1, maxspec_renamexf
             lspecfrma_renamexf_c(iq,ipair) = int(lspecfrma_renamexf_in(iq,ipair), c_int64_t)
             lspecfrmc_renamexf_c(iq,ipair) = int(lspecfrmc_renamexf_in(iq,ipair), c_int64_t)
             lspectooa_renamexf_c(iq,ipair) = int(lspectooa_renamexf_in(iq,ipair), c_int64_t)
             lspectooc_renamexf_c(iq,ipair) = int(lspectooc_renamexf_in(iq,ipair), c_int64_t)
          end do
       end do

       call modal_aero_rename_set_dotend_flags_codon( &
            int(pcnstxx, c_int64_t), int(maxpair_renamexf, c_int64_t), int(maxspec_renamexf, c_int64_t), &
            int(loffset, c_int64_t), int(npair_renamexf_in, c_int64_t), c_loc(nspecfrm_renamexf_c(1)), &
            c_loc(lspecfrma_renamexf_c(1,1)), c_loc(lspecfrmc_renamexf_c(1,1)), c_loc(lspectooa_renamexf_c(1,1)), &
            c_loc(lspectooc_renamexf_c(1,1)), c_loc(dotendrn_c(1)), c_loc(dotendqqcwrn_c(1)) &
       )
       do l = 1, pcnstxx
          dotendrn(l) = dotendrn_c(l) /= 0_c_int64_t
          dotendqqcwrn(l) = dotendqqcwrn_c(l) /= 0_c_int64_t
       end do
       return
    end if

    call modal_aero_rename_set_dotend_flags_native( loffset, npair_renamexf_in, &
         nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, &
         lspectooa_renamexf_in, lspectooc_renamexf_in, dotendrn, dotendqqcwrn )

  end subroutine modal_aero_rename_set_dotend_flags

!----------------------------------------------------------------------
!----------------------------------------------------------------------
  subroutine modal_aero_rename_set_dotend_flags_native( loffset, npair_renamexf_in, &
       nspecfrm_renamexf_in, lspecfrma_renamexf_in, lspecfrmc_renamexf_in, &
       lspectooa_renamexf_in, lspectooc_renamexf_in, dotendrn, dotendqqcwrn )

    implicit none

    integer, intent(in) :: loffset
    integer, intent(in) :: npair_renamexf_in
    integer, intent(in) :: nspecfrm_renamexf_in(maxpair_renamexf)
    integer, intent(in) :: lspecfrma_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, intent(in) :: lspecfrmc_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, intent(in) :: lspectooa_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    integer, intent(in) :: lspectooc_renamexf_in(maxspec_renamexf,maxpair_renamexf)
    logical, intent(out) :: dotendrn(pcnstxx)
    logical, intent(out) :: dotendqqcwrn(pcnstxx)

    integer :: ipair, iq
    integer :: lsfrma, lsfrmc, lstooa, lstooc

    dotendrn(:) = .false.
    dotendqqcwrn(:) = .false.
    do ipair = 1, npair_renamexf_in
       do iq = 1, nspecfrm_renamexf_in(ipair)
          lsfrma = lspecfrma_renamexf_in(iq,ipair) - loffset
          lsfrmc = lspecfrmc_renamexf_in(iq,ipair) - loffset
          lstooa = lspectooa_renamexf_in(iq,ipair) - loffset
          lstooc = lspectooc_renamexf_in(iq,ipair) - loffset
          if (lsfrma .gt. 0) then
             dotendrn(lsfrma) = .true.
             if (lstooa .gt. 0) dotendrn(lstooa) = .true.
          end if
          if (lsfrmc .gt. 0) then
             dotendqqcwrn(lsfrmc) = .true.
             if (lstooc .gt. 0) dotendqqcwrn(lstooc) = .true.
          end if
       end do
    end do

  end subroutine modal_aero_rename_set_dotend_flags_native



!-------------------------------------------------------------------------
	subroutine modal_aero_rename_no_acc_crs_init
!
!   computes pointers for species transfer during aerosol renaming
!	(a2 --> a1 transfer)
!   transfers include number_a, number_c, mass_a, mass_c and
!	water_a
!

	implicit none

!   local variables
	integer ipair, iq, iqfrm, iqfrm_aa, iqtoo, iqtoo_aa,   &
      	  lsfrma, lsfrmc, lstooa, lstooc, lunout,   &
      	  mfrm, mtoo, n1, n2, nsamefrm, nsametoo, nspec


	lunout = iulog
!
!   define "from mode" and "to mode" for each tail-xfer pairing
!	currently just a2-->a1
!
	n1 = modeptr_accum
	n2 = modeptr_aitken
	if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
	    npair_renamexf = 1
	    modefrm_renamexf(1) = n2
	    modetoo_renamexf(1) = n1
	else
	    npair_renamexf = 0
	    return
	end if

!
!   define species involved in each tail-xfer pairing
!	(include aerosol water)
!
	do 1900 ipair = 1, npair_renamexf
	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)

	nspec = 0
	do 1490 iqfrm = -1, nspec_amode(mfrm)
	    iqtoo = iqfrm
	    if (iqfrm .eq. -1) then
		lsfrma = numptr_amode(mfrm)
		lstooa = numptr_amode(mtoo)
		lsfrmc = numptrcw_amode(mfrm)
		lstooc = numptrcw_amode(mtoo)
	    else if (iqfrm .eq. 0) then
!   bypass transfer of aerosol water due to renaming
                goto 1490
!               lsfrma = lwaterptr_amode(mfrm)
!               lsfrmc = 0
!               lstooa = lwaterptr_amode(mtoo)
!               lstooc = 0
	    else
		lsfrma = lmassptr_amode(iqfrm,mfrm)
		lsfrmc = lmassptrcw_amode(iqfrm,mfrm)
		lstooa = 0
		lstooc = 0
	    end if

	    if ((lsfrma .lt. 1) .or. (lsfrma .gt. pcnst)) then
		write(lunout,9100) mfrm, iqfrm, lsfrma
		call endrun( 'modal_aero_rename_no_acc_crs_init error' )
	    end if
	    if (iqfrm .le. 0) goto 1430

	    if ((lsfrmc .lt. 1) .or. (lsfrmc .gt. pcnst)) then
		write(lunout,9102) mfrm, iqfrm, lsfrmc
		call endrun( 'modal_aero_rename_no_acc_crs_init error' )
	    end if

! find "too" species having same lspectype_amode as the "frm" species
! several species in a mode may have the same lspectype_amode, so also
!    use the ordering as a criterion (e.g., 1st <--> 1st, 2nd <--> 2nd)
	    iqfrm_aa = 1
	    iqtoo_aa = 1
	    if (iqfrm .gt. nspec_amode(mfrm)) then
		iqfrm_aa = nspec_amode(mfrm) + 1
		iqtoo_aa = nspec_amode(mtoo) + 1
	    end if
	    nsamefrm = 0
	    do iq = iqfrm_aa, iqfrm
		if ( lspectype_amode(iq   ,mfrm) .eq.   &
      		     lspectype_amode(iqfrm,mfrm) ) then
		    nsamefrm = nsamefrm + 1
		end if
	    end do
	    nsametoo = 0
	    do iqtoo = iqtoo_aa, nspec_amode(mtoo)
		if ( lspectype_amode(iqtoo,mtoo) .eq.   &
      		     lspectype_amode(iqfrm,mfrm) ) then
		    nsametoo = nsametoo + 1
		    if (nsametoo .eq. nsamefrm) then
			lstooc = lmassptrcw_amode(iqtoo,mtoo)
			lstooa = lmassptr_amode(iqtoo,mtoo)
			goto 1430
		    end if
		end if
	    end do

1430	    nspec = nspec + 1
	    if ((lstooc .lt. 1) .or. (lstooc .gt. pcnst)) lstooc = 0
	    if ((lstooa .lt. 1) .or. (lstooa .gt. pcnst)) lstooa = 0
	    if (lstooa .eq. 0) then
		write(lunout,9104) mfrm, iqfrm, lsfrma, iqtoo, lstooa
		call endrun( 'modal_aero_rename_no_acc_crs_init error' )
	    end if
	    if ((lstooc .eq. 0) .and. (iqfrm .ne. 0)) then
		write(lunout,9104) mfrm, iqfrm, lsfrmc, iqtoo, lstooc
		call endrun( 'modal_aero_rename_no_acc_crs_init error' )
	    end if
	    lspecfrma_renamexf(nspec,ipair) = lsfrma
	    lspectooa_renamexf(nspec,ipair) = lstooa
	    lspecfrmc_renamexf(nspec,ipair) = lsfrmc
	    lspectooc_renamexf(nspec,ipair) = lstooc
1490	continue

	nspecfrm_renamexf(ipair) = nspec
1900	continue

9100	format( / '*** subr. modal_aero_rename_no_acc_crs_init' /   &
      	'lspecfrma out of range' /   &
      	'modefrm, ispecfrm, lspecfrma =', 3i6 / )
9102	format( / '*** subr. modal_aero_rename_no_acc_crs_init' /   &
      	'lspecfrmc out of range' /   &
      	'modefrm, ispecfrm, lspecfrmc =', 3i6 / )
9104	format( / '*** subr. modal_aero_rename_no_acc_crs_init' /   &
      	'lspectooa out of range' /   &
      	'modefrm, ispecfrm, lspecfrma, ispectoo, lspectooa =', 5i6 / )
9106	format( / '*** subr. modal_aero_rename_no_acc_crs_init' /   &
      	'lspectooc out of range' /   &
      	'modefrm, ispecfrm, lspecfrmc, ispectoo, lspectooc =', 5i6 / )

!
!   output results
!
	if ( masterproc ) then

	write(lunout,9310)

	do 2900 ipair = 1, npair_renamexf
	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)
	write(lunout,9320) ipair, mfrm, mtoo

	do iq = 1, nspecfrm_renamexf(ipair)
	    lsfrma = lspecfrma_renamexf(iq,ipair)
	    lstooa = lspectooa_renamexf(iq,ipair)
	    lsfrmc = lspecfrmc_renamexf(iq,ipair)
	    lstooc = lspectooc_renamexf(iq,ipair)
	    if (lstooa .gt. 0) then
		write(lunout,9330) lsfrma, cnst_name(lsfrma),   &
				   lstooa, cnst_name(lstooa)
	    else
		write(lunout,9340) lsfrma, cnst_name(lsfrma)
	    end if
	    if (lstooc .gt. 0) then
		write(lunout,9330) lsfrmc, cnst_name_cw(lsfrmc),   &
				   lstooc, cnst_name_cw(lstooc)
	    else if (lsfrmc .gt. 0) then
		write(lunout,9340) lsfrmc, cnst_name_cw(lsfrmc)
	    else
		write(lunout,9350)
	    end if
	end do

2900	continue
	write(lunout,*)

	end if ! ( masterproc )

9310	format( / 'subr. modal_aero_rename_no_acc_crs_init' )
9320	format( 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3 )
9330	format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340	format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )
9350	format( 5x, 'no corresponding activated species' )

	return
	end subroutine modal_aero_rename_no_acc_crs_init

!----------------------------------------------------------------------
! code for troposphere and stratosphere
! -- allows accumulation to coarse mode exchange
!----------------------------------------------------------------------
!BOP
! !ROUTINE:  modal_aero_rename_acc_crs_sub --- ...
!
! !INTERFACE:
	subroutine modal_aero_rename_acc_crs_sub(                       &
                        fromwhere,         lchnk,               &
                        ncol,              nstep,               &
                        loffset,           deltat,              &
                        pdel,              troplev,             &
                        dotendrn,          q,                   &
                        dqdt,              dqdt_other,          &
                        dotendqqcwrn,      qqcw,                &
                        dqqcwdt,           dqqcwdt_other,       &
                        is_dorename_atik,  dorename_atik,       &
                        jsrflx_rename,     nsrflx,              &
                        qsrflx,            qqcwsrflx,           &
                        dqdt_rnpos                              )

! !USES:

   use physconst, only: gravit, mwdry
   use units, only: getunit
   use shr_spfn_mod, only: erfc => shr_spfn_erfc
   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

   implicit none


! !PARAMETERS:
   character(len=*), intent(in) :: fromwhere    ! identifies which module
                                                ! is making the call
   integer,  intent(in)    :: lchnk                ! chunk identifier
   integer,  intent(in)    :: ncol                 ! number of atmospheric column
   integer,  intent(in)    :: nstep                ! model time-step number
   integer,  intent(in)    :: loffset              ! offset applied to modal aero "ptrs"
   real(r8), intent(in)    :: deltat               ! time step (s)
   integer,  intent(in)    :: troplev(pcols)

   real(r8), target, intent(in) :: pdel(pcols,pver)     ! pressure thickness of levels (Pa)
   real(r8), target, intent(in) :: q(ncol,pver,pcnstxx) ! tracer mixing ratio array
                                                   ! *** MUST BE mol/mol-air or #/mol-air
                                                   ! *** NOTE ncol and pcnstxx dimensions
   real(r8), target, intent(in) :: qqcw(ncol,pver,pcnstxx) ! like q but for cloud-borne species

   real(r8), target, intent(inout) :: dqdt(ncol,pver,pcnstxx)  ! TMR tendency array;
                              ! incoming dqdt = tendencies for the 
                              !     "fromwhere" continuous growth process 
                              ! the renaming tendencies are added on
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), target, intent(inout) :: dqqcwdt(ncol,pver,pcnstxx)
   real(r8), target, intent(in) :: dqdt_other(ncol,pver,pcnstxx)  
                              ! tendencies for "other" continuous growth process 
                              ! currently in cam3
                              !     dqdt is from gas (h2so4, nh3) condensation
                              !     dqdt_other is from aqchem and soa
                              ! *** NOTE ncol and pcnstxx dimensions
   real(r8), target, intent(in) :: dqqcwdt_other(ncol,pver,pcnstxx)  
   logical,  intent(inout) :: dotendrn(pcnstxx) ! identifies the species for which
                              !     renaming dqdt is computed
   logical,  intent(inout) :: dotendqqcwrn(pcnstxx)

   logical,  intent(in)    :: is_dorename_atik          ! true if dorename_atik is provided
   logical,  intent(in)    :: dorename_atik(ncol,pver) ! true if renaming should
                                                        ! be done at i,k
   integer,  intent(in)    :: jsrflx_rename        ! qsrflx index for renaming
   integer,  intent(in)    :: nsrflx               ! last dimension of qsrflx

   real(r8), target, intent(inout) :: qsrflx(pcols,pcnstxx,nsrflx)
                              ! process-specific column tracer tendencies 
   real(r8), target, intent(inout) :: qqcwsrflx(pcols,pcnstxx,nsrflx)
   real(r8), optional, target, intent(out) &
                           :: dqdt_rnpos(ncol,pver,pcnstxx)
                              ! the positive (production) part of the renaming tendency

! !DESCRIPTION: 
! computes TMR (tracer mixing ratio) tendencies for "mode renaming"
!    during a continuous growth process
! currently this transfers number and mass (and surface) from the aitken
!    to accumulation mode after gas condensation or stratiform-cloud
!    aqueous chemistry
! (convective cloud aqueous chemistry not yet implemented)
!
! !REVISION HISTORY:
!   RCE 07.04.13:  Adapted from MIRAGE2 code
!
!EOP
!----------------------------------------------------------------------
!BOC

! local variables
   integer, parameter :: ldiag1 = -1
   integer :: i, icol_diag, ipair, iq
   integer :: j, k
   integer :: l, l1, l2, la, lc, lunout
   integer :: lsfrma, lsfrmc, lstooa, lstooc
   integer :: mfrm, mtoo, n, n1, n2, ntot_msa_a
   integer, save :: lun = -1  ! logical unit for diagnostics (6, or other
                              ! if a special diagnostics file is opened)

   logical :: l_dqdt_rnpos
   logical :: flagaa_shrink, flagbb_shrink

   real (r8), target :: deldryvol_a(ncol,pver)
   real (r8), target :: deldryvol_c(ncol,pver)
   real (r8) :: deltatinv
   real (r8) :: dgn_aftr, dgn_xfer
   real (r8) :: dgn_t_new, dgn_t_old, dgn_t_oldb
   real (r8) :: dryvol_t_del, dryvol_t_new, dryvol_t_new_xfab
   real (r8) :: dryvol_t_old, dryvol_t_oldb, dryvol_t_oldbnd
   real (r8), target :: dryvol_a(ncol,pver)
   real (r8), target :: dryvol_c(ncol,pver)
   real (r8), target :: dryvol_a_xfab(ncol,pver)
   real (r8), target :: dryvol_c_xfab(ncol,pver)
   real (r8) :: dryvol_xferamt
   real (r8) :: lndgn_new, lndgn_old
   real (r8) :: lndgv_new, lndgv_old
   real (r8) :: num_t_old, num_t_oldbnd
   real (r8) :: onethird
   real (r8) :: pdel_fac
   real (r8) :: tailfr_volnew, tailfr_volold
   real (r8) :: tailfr_numnew, tailfr_numold
   real (r8) :: tmpa, tmpf
   real (r8) :: tmp_m2v, tmp_m2vdt
   real (r8) :: xfercoef, xfertend
   real (r8) :: xferfrac_vol, xferfrac_num, xferfrac_max
   real (r8), target :: xferfrac_vol_ik(ncol,pver)
   real (r8), target :: xferfrac_num_ik(ncol,pver)

   real (r8) :: yn_tail, yv_tail
   integer(c_int64_t), target :: troplev_c(pcols)
   integer(c_int64_t), target :: dorename_atik_c(ncol,pver)
   integer(c_int64_t), target :: modefrm_renamexf_c(maxpair_renamexf)
   integer(c_int64_t), target :: modetoo_renamexf_c(maxpair_renamexf)
   integer(c_int64_t), target :: nspec_amode_c(ntot_amode)
   integer(c_int64_t), target :: lspectype_amode_c(maxspec_renamexf,ntot_amode)
   integer(c_int64_t), target :: lmassptr_amode_c(maxspec_renamexf,ntot_amode)
   integer(c_int64_t), target :: lmassptrcw_amode_c(maxspec_renamexf,ntot_amode)
   integer(c_int64_t), target :: numptr_amode_c(ntot_amode)
   integer(c_int64_t), target :: numptrcw_amode_c(ntot_amode)
   integer(c_int64_t), target :: igrow_shrink_renamexf_c(maxpair_renamexf)
   integer(c_int64_t), target :: ixferable_all_renamexf_c(maxpair_renamexf)
   integer(c_int64_t), target :: ixferable_a_renamexf_c(maxspec_renamexf,maxpair_renamexf)
   integer(c_int64_t), target :: ixferable_c_renamexf_c(maxspec_renamexf,maxpair_renamexf)
   integer(c_int64_t), target :: nspecfrm_renamexf_c(maxpair_renamexf)
   integer(c_int64_t), target :: lspecfrma_renamexf_c(maxspec_renamexf,maxpair_renamexf)
   integer(c_int64_t), target :: lspecfrmc_renamexf_c(maxspec_renamexf,maxpair_renamexf)
   integer(c_int64_t), target :: lspectooa_renamexf_c(maxspec_renamexf,maxpair_renamexf)
   integer(c_int64_t), target :: lspectooc_renamexf_c(maxspec_renamexf,maxpair_renamexf)
   integer(c_int64_t), target :: dotendrn_c(pcnstxx)
   integer(c_int64_t), target :: dotendqqcwrn_c(pcnstxx)
   integer(c_int64_t) :: is_dorename_atik_c
   integer(c_int64_t) :: l_dqdt_rnpos_c
   real(r8), target :: specmw_amode_c(size(specmw_amode))
   real(r8), target :: specdens_amode_c(size(specdens_amode))
   real(r8), target :: dgnum_amode_c(ntot_amode)
   real(r8), target :: factoraa_c(ntot_amode)
   real(r8), target :: factoryy_c(ntot_amode)
   real(r8), target :: dryvol_smallest_c(ntot_amode)
   real(r8), target :: v2nlorlx_c(ntot_amode)
   real(r8), target :: v2nhirlx_c(ntot_amode)
   real(r8), target :: factor_3alnsg2_c(maxpair_renamexf)
   real(r8), target :: dp_cut_c(maxpair_renamexf)
   real(r8), target :: lndp_cut_c(maxpair_renamexf)
   real(r8), target :: dp_belowcut_c(maxpair_renamexf)
   real(r8), target :: dp_xfernone_threshaa_c(maxpair_renamexf)
   real(r8), target :: dp_xferall_thresh_c(maxpair_renamexf)
   real(r8), target :: dqdt_rnpos_dummy(1,1,1)
   type(c_ptr) :: dqdt_rnpos_p

   interface
      subroutine modal_aero_rename_acc_crs_sub_codon( &
           ncol_c, pcols_c, pver_c, pcnstxx_c, maxpair_renamexf_c, maxspec_renamexf_c, loffset_c, &
           npair_renamexf_c, is_dorename_atik_c, l_dqdt_rnpos_c, jsrflx_rename_c, nsrflx_c, &
           modeptr_coarse_c, modeptr_accum_c, method_optbb_c, deltat_c, deltatinv_c, onethird_c, &
           xferfrac_max_c, gravit_c, troplev_p, pdel_p, dorename_atik_p, q_p, qqcw_p, dqdt_p, &
           dqdt_other_p, dqqcwdt_p, dqqcwdt_other_p, qsrflx_p, qqcwsrflx_p, modefrm_renamexf_p, &
           modetoo_renamexf_p, nspec_amode_p, lspectype_amode_p, specmw_amode_p, specdens_amode_p, &
           lmassptr_amode_p, lmassptrcw_amode_p, numptr_amode_p, numptrcw_amode_p, dgnum_amode_p, &
           factoraa_p, factoryy_p, dryvol_smallest_p, v2nlorlx_p, v2nhirlx_p, factor_3alnsg2_p, &
           dp_cut_p, lndp_cut_p, dp_belowcut_p, dp_xfernone_threshaa_p, dp_xferall_thresh_p, &
           igrow_shrink_renamexf_p, ixferable_all_renamexf_p, ixferable_a_renamexf_p, ixferable_c_renamexf_p, &
           nspecfrm_renamexf_p, lspecfrma_renamexf_p, lspecfrmc_renamexf_p, lspectooa_renamexf_p, &
           lspectooc_renamexf_p, dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p, dryvol_a_xfab_p, &
           dryvol_c_xfab_p, xferfrac_vol_p, xferfrac_num_p, dotendrn_p, dotendqqcwrn_p, dqdt_rnpos_p ) &
           bind(c, name="modal_aero_rename_acc_crs_sub_codon")
        use iso_c_binding, only: c_double, c_int64_t, c_ptr
        integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnstxx_c, maxpair_renamexf_c, maxspec_renamexf_c
        integer(c_int64_t), value :: loffset_c, npair_renamexf_c, is_dorename_atik_c, l_dqdt_rnpos_c
        integer(c_int64_t), value :: jsrflx_rename_c, nsrflx_c, modeptr_coarse_c, modeptr_accum_c, method_optbb_c
        real(c_double), value :: deltat_c, deltatinv_c, onethird_c, xferfrac_max_c, gravit_c
        type(c_ptr), value :: troplev_p, pdel_p, dorename_atik_p, q_p, qqcw_p, dqdt_p, dqdt_other_p
        type(c_ptr), value :: dqqcwdt_p, dqqcwdt_other_p, qsrflx_p, qqcwsrflx_p
        type(c_ptr), value :: modefrm_renamexf_p, modetoo_renamexf_p, nspec_amode_p, lspectype_amode_p
        type(c_ptr), value :: specmw_amode_p, specdens_amode_p, lmassptr_amode_p, lmassptrcw_amode_p
        type(c_ptr), value :: numptr_amode_p, numptrcw_amode_p, dgnum_amode_p, factoraa_p, factoryy_p
        type(c_ptr), value :: dryvol_smallest_p, v2nlorlx_p, v2nhirlx_p, factor_3alnsg2_p, dp_cut_p
        type(c_ptr), value :: lndp_cut_p, dp_belowcut_p, dp_xfernone_threshaa_p, dp_xferall_thresh_p
        type(c_ptr), value :: igrow_shrink_renamexf_p, ixferable_all_renamexf_p, ixferable_a_renamexf_p
        type(c_ptr), value :: ixferable_c_renamexf_p, nspecfrm_renamexf_p, lspecfrma_renamexf_p
        type(c_ptr), value :: lspecfrmc_renamexf_p, lspectooa_renamexf_p, lspectooc_renamexf_p
        type(c_ptr), value :: dryvol_a_p, dryvol_c_p, deldryvol_a_p, deldryvol_c_p, dryvol_a_xfab_p
        type(c_ptr), value :: dryvol_c_xfab_p, xferfrac_vol_p, xferfrac_num_p, dotendrn_p, dotendqqcwrn_p
        type(c_ptr), value :: dqdt_rnpos_p
      end subroutine modal_aero_rename_acc_crs_sub_codon
   end interface

! begin
	lunout = iulog

!   get logical unit (for output to dumpconv, deactivate the "lun = 6")
 	lun = iulog
	if (lun < 1) then
	   lun = getunit()
 	   open( unit=lun, file='dump.rename',   &
 			status='unknown', form='formatted' )
	end if


!
!   calculations done once on initial entry
!
!   "init" is now done through chem_init (and things under it)
!	if (npair_renamexf .eq. -123456789) then
!	    npair_renamexf = 0
!	    call modal_aero_rename_init
!	end if

!
!   check if any renaming pairs exist
!
	if (npair_renamexf .le. 0) return
! 	if (ncol .ne. -123456789) return
!	if (fromwhere .eq. 'aqchem') return


	deltatinv = 1.0_r8/(deltat*(1.0_r8 + 1.0e-15_r8))
	onethird = 1.0_r8/3.0_r8
	xferfrac_max = 1.0_r8 - 10.0_r8*epsilon(1.0_r8)   ! 1-eps

	if ( present( dqdt_rnpos ) ) then
	    l_dqdt_rnpos = .true.
	    dqdt_rnpos(:,:,:) = 0.0_r8
	else
	    l_dqdt_rnpos = .false.
	end if

        call modal_aero_rename_acc_crs_sub_select_impl()

        if (.not. modal_aero_rename_acc_crs_sub_use_native_impl) then
           if (is_dorename_atik) then
              is_dorename_atik_c = 1_c_int64_t
              do k = 1, pver
                 do i = 1, ncol
                    if (dorename_atik(i,k)) then
                       dorename_atik_c(i,k) = 1_c_int64_t
                    else
                       dorename_atik_c(i,k) = 0_c_int64_t
                    end if
                 end do
              end do
           else
              is_dorename_atik_c = 0_c_int64_t
              dorename_atik_c(:,:) = 0_c_int64_t
           end if

           if (l_dqdt_rnpos) then
              l_dqdt_rnpos_c = 1_c_int64_t
           else
              l_dqdt_rnpos_c = 0_c_int64_t
           end if

           do i = 1, pcols
              troplev_c(i) = int(troplev(i), c_int64_t)
           end do
           do ipair = 1, maxpair_renamexf
              modefrm_renamexf_c(ipair) = int(modefrm_renamexf(ipair), c_int64_t)
              modetoo_renamexf_c(ipair) = int(modetoo_renamexf(ipair), c_int64_t)
              igrow_shrink_renamexf_c(ipair) = int(igrow_shrink_renamexf(ipair), c_int64_t)
              ixferable_all_renamexf_c(ipair) = int(ixferable_all_renamexf(ipair), c_int64_t)
              nspecfrm_renamexf_c(ipair) = int(nspecfrm_renamexf(ipair), c_int64_t)
              factor_3alnsg2_c(ipair) = factor_3alnsg2(ipair)
              dp_cut_c(ipair) = dp_cut(ipair)
              lndp_cut_c(ipair) = lndp_cut(ipair)
              dp_belowcut_c(ipair) = dp_belowcut(ipair)
              dp_xfernone_threshaa_c(ipair) = dp_xfernone_threshaa(ipair)
              dp_xferall_thresh_c(ipair) = dp_xferall_thresh(ipair)
              do iq = 1, maxspec_renamexf
                 ixferable_a_renamexf_c(iq,ipair) = int(ixferable_a_renamexf(iq,ipair), c_int64_t)
                 ixferable_c_renamexf_c(iq,ipair) = int(ixferable_c_renamexf(iq,ipair), c_int64_t)
                 lspecfrma_renamexf_c(iq,ipair) = int(lspecfrma_renamexf(iq,ipair), c_int64_t)
                 lspecfrmc_renamexf_c(iq,ipair) = int(lspecfrmc_renamexf(iq,ipair), c_int64_t)
                 lspectooa_renamexf_c(iq,ipair) = int(lspectooa_renamexf(iq,ipair), c_int64_t)
                 lspectooc_renamexf_c(iq,ipair) = int(lspectooc_renamexf(iq,ipair), c_int64_t)
              end do
           end do
           do n = 1, ntot_amode
              nspec_amode_c(n) = int(nspec_amode(n), c_int64_t)
              numptr_amode_c(n) = int(numptr_amode(n), c_int64_t)
              numptrcw_amode_c(n) = int(numptrcw_amode(n), c_int64_t)
              dgnum_amode_c(n) = dgnum_amode(n)
              factoraa_c(n) = factoraa(n)
              factoryy_c(n) = factoryy(n)
              dryvol_smallest_c(n) = dryvol_smallest(n)
              v2nlorlx_c(n) = v2nlorlx(n)
              v2nhirlx_c(n) = v2nhirlx(n)
              do iq = 1, maxspec_renamexf
                 lspectype_amode_c(iq,n) = int(lspectype_amode(iq,n), c_int64_t)
                 lmassptr_amode_c(iq,n) = int(lmassptr_amode(iq,n), c_int64_t)
                 lmassptrcw_amode_c(iq,n) = int(lmassptrcw_amode(iq,n), c_int64_t)
              end do
           end do
           specmw_amode_c(:) = specmw_amode(:)
           specdens_amode_c(:) = specdens_amode(:)

           if (present(dqdt_rnpos)) then
              dqdt_rnpos_p = c_loc(dqdt_rnpos(1,1,1))
           else
              dqdt_rnpos_dummy(1,1,1) = 0.0_r8
              dqdt_rnpos_p = c_loc(dqdt_rnpos_dummy(1,1,1))
           end if

           call modal_aero_rename_acc_crs_sub_codon( &
                int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pcnstxx, c_int64_t), &
                int(maxpair_renamexf, c_int64_t), int(maxspec_renamexf, c_int64_t), int(loffset, c_int64_t), &
                int(npair_renamexf, c_int64_t), is_dorename_atik_c, l_dqdt_rnpos_c, int(jsrflx_rename, c_int64_t), &
                int(nsrflx, c_int64_t), int(modeptr_coarse, c_int64_t), int(modeptr_accum, c_int64_t), &
                int(method_optbb_renamexf, c_int64_t), real(deltat, c_double), real(deltatinv, c_double), &
                real(onethird, c_double), real(xferfrac_max, c_double), real(gravit, c_double), c_loc(troplev_c(1)), &
                c_loc(pdel(1,1)), c_loc(dorename_atik_c(1,1)), c_loc(q(1,1,1)), c_loc(qqcw(1,1,1)), c_loc(dqdt(1,1,1)), &
                c_loc(dqdt_other(1,1,1)), c_loc(dqqcwdt(1,1,1)), c_loc(dqqcwdt_other(1,1,1)), c_loc(qsrflx(1,1,1)), &
                c_loc(qqcwsrflx(1,1,1)), c_loc(modefrm_renamexf_c(1)), c_loc(modetoo_renamexf_c(1)), c_loc(nspec_amode_c(1)), &
                c_loc(lspectype_amode_c(1,1)), c_loc(specmw_amode_c(1)), c_loc(specdens_amode_c(1)), c_loc(lmassptr_amode_c(1,1)), &
                c_loc(lmassptrcw_amode_c(1,1)), c_loc(numptr_amode_c(1)), c_loc(numptrcw_amode_c(1)), c_loc(dgnum_amode_c(1)), &
                c_loc(factoraa_c(1)), c_loc(factoryy_c(1)), c_loc(dryvol_smallest_c(1)), c_loc(v2nlorlx_c(1)), c_loc(v2nhirlx_c(1)), &
                c_loc(factor_3alnsg2_c(1)), c_loc(dp_cut_c(1)), c_loc(lndp_cut_c(1)), c_loc(dp_belowcut_c(1)), &
                c_loc(dp_xfernone_threshaa_c(1)), c_loc(dp_xferall_thresh_c(1)), c_loc(igrow_shrink_renamexf_c(1)), &
                c_loc(ixferable_all_renamexf_c(1)), c_loc(ixferable_a_renamexf_c(1,1)), c_loc(ixferable_c_renamexf_c(1,1)), &
                c_loc(nspecfrm_renamexf_c(1)), c_loc(lspecfrma_renamexf_c(1,1)), c_loc(lspecfrmc_renamexf_c(1,1)), &
                c_loc(lspectooa_renamexf_c(1,1)), c_loc(lspectooc_renamexf_c(1,1)), c_loc(dryvol_a(1,1)), c_loc(dryvol_c(1,1)), &
                c_loc(deldryvol_a(1,1)), c_loc(deldryvol_c(1,1)), c_loc(dryvol_a_xfab(1,1)), c_loc(dryvol_c_xfab(1,1)), &
                c_loc(xferfrac_vol_ik(1,1)), c_loc(xferfrac_num_ik(1,1)), c_loc(dotendrn_c(1)), c_loc(dotendqqcwrn_c(1)), &
                dqdt_rnpos_p )

           do l = 1, pcnstxx
              dotendrn(l) = dotendrn_c(l) /= 0_c_int64_t
              dotendqqcwrn(l) = dotendqqcwrn_c(l) /= 0_c_int64_t
           end do
           return
        end if



!
!   loop over renaming pairs
!
mainloop1_ipair:  do ipair = 1, npair_renamexf

	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)

	flagaa_shrink = .false.
	if ((mfrm==modeptr_coarse) .and. (mtoo==modeptr_accum)) &
	    flagaa_shrink = .true.

        if (l_dqdt_rnpos) then
           call modal_aero_rename_acc_crs_pair( ncol, loffset, deltat, deltatinv, troplev, pdel, q, qqcw, &
                dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, is_dorename_atik, dorename_atik, jsrflx_rename, &
                nsrflx, qsrflx, qqcwsrflx, ixferable_all_renamexf(ipair), nspec_amode(mfrm), lspectype_amode(:,mfrm), &
                specmw_amode, specdens_amode, lmassptr_amode(:,mfrm), lmassptrcw_amode(:,mfrm), ixferable_a_renamexf(:,ipair), &
                ixferable_c_renamexf(:,ipair), mfrm, numptr_amode(mfrm), numptrcw_amode(mfrm), dgnum_amode(mfrm), &
                factoraa(mfrm), factoryy(mfrm), dryvol_smallest(mfrm), v2nlorlx(mfrm), v2nhirlx(mfrm), factor_3alnsg2(ipair), &
                dp_cut(ipair), lndp_cut(ipair), dp_belowcut(ipair), dp_xfernone_threshaa(ipair), dp_xferall_thresh(ipair), &
                flagaa_shrink, igrow_shrink_renamexf(ipair), method_optbb_renamexf, onethird, xferfrac_max, &
                nspecfrm_renamexf(ipair), lspecfrma_renamexf(:,ipair), lspecfrmc_renamexf(:,ipair), lspectooa_renamexf(:,ipair), &
                lspectooc_renamexf(:,ipair), dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab, &
                xferfrac_vol_ik, xferfrac_num_ik, l_dqdt_rnpos, dqdt_rnpos )
        else
           call modal_aero_rename_acc_crs_pair( ncol, loffset, deltat, deltatinv, troplev, pdel, q, qqcw, &
                dqdt, dqdt_other, dqqcwdt, dqqcwdt_other, is_dorename_atik, dorename_atik, jsrflx_rename, &
                nsrflx, qsrflx, qqcwsrflx, ixferable_all_renamexf(ipair), nspec_amode(mfrm), lspectype_amode(:,mfrm), &
                specmw_amode, specdens_amode, lmassptr_amode(:,mfrm), lmassptrcw_amode(:,mfrm), ixferable_a_renamexf(:,ipair), &
                ixferable_c_renamexf(:,ipair), mfrm, numptr_amode(mfrm), numptrcw_amode(mfrm), dgnum_amode(mfrm), &
                factoraa(mfrm), factoryy(mfrm), dryvol_smallest(mfrm), v2nlorlx(mfrm), v2nhirlx(mfrm), factor_3alnsg2(ipair), &
                dp_cut(ipair), lndp_cut(ipair), dp_belowcut(ipair), dp_xfernone_threshaa(ipair), dp_xferall_thresh(ipair), &
                flagaa_shrink, igrow_shrink_renamexf(ipair), method_optbb_renamexf, onethird, xferfrac_max, &
                nspecfrm_renamexf(ipair), lspecfrma_renamexf(:,ipair), lspecfrmc_renamexf(:,ipair), lspectooa_renamexf(:,ipair), &
                lspectooc_renamexf(:,ipair), dryvol_a, dryvol_c, deldryvol_a, deldryvol_c, dryvol_a_xfab, dryvol_c_xfab, &
                xferfrac_vol_ik, xferfrac_num_ik, l_dqdt_rnpos )
        end if


	end do mainloop1_ipair

!
!   set dotend's
!
	call modal_aero_rename_set_dotend_flags( loffset, npair_renamexf, nspecfrm_renamexf, &
	     lspecfrma_renamexf, lspecfrmc_renamexf, lspectooa_renamexf, lspectooc_renamexf, &
	     dotendrn, dotendqqcwrn )


	return


!
!   error -- renaming currently just works for 1 pair
!
8100	write(lunout,9050) ipair
	call endrun( 'modal_aero_rename_acc_crs_sub error' )
9050	format( / '*** subr. modal_aero_rename_acc_crs_sub ***' /   &
      	    4x, 'aerosol renaming not implemented for ipair =', i5 )

!EOC
	end subroutine modal_aero_rename_acc_crs_sub



!-------------------------------------------------------------------------
! for modal aerosols in the troposphere and stratophere
! -- allows accumulation to coarse mode exchange
!-------------------------------------------------------------------------
	subroutine modal_aero_rename_acc_crs_init
!
!   computes pointers for species transfer during aerosol renaming
!	(a2 --> a1 transfer)
!   transfers include number_a, number_c, mass_a, mass_c and
!	water_a
!

	implicit none

!   local variables
	integer :: i, ipair, iq, iqfrm, iqtooa, iqtooc, itmpa
	integer :: l, lsfrma, lsfrmc, lstooa, lstooc, lunout
	integer :: mfrm, mtoo
	integer :: n1, n2, nspec
	integer :: nch_lfrm, nch_ltoo, nch_mfrmid, nch_mtooid

	real (r8) :: frelax

	lunout = iulog

!
!   define "from mode" and "to mode" for each tail-xfer pairing
!	using the values in ipair_select_renamexf(:)
!
	npair_renamexf = 0
	do ipair = 1, maxpair_renamexf
	    itmpa = ipair_select_renamexf(ipair)
	    if (itmpa == 0) then
		exit
	    else if (itmpa == 2001) then
		mfrm = modeptr_aitken
		mtoo = modeptr_accum
		igrow_shrink_renamexf(ipair) = 1
		ixferable_all_needed_renamexf(ipair) = 1
	    else if (itmpa == 1003) then
		mfrm = modeptr_accum
		mtoo = modeptr_coarse
		igrow_shrink_renamexf(ipair) = 1
		ixferable_all_needed_renamexf(ipair) = 0
	    else if (itmpa == 3001) then
		mfrm = modeptr_coarse
		mtoo = modeptr_accum
		igrow_shrink_renamexf(ipair) = -1
		ixferable_all_needed_renamexf(ipair) = 0
	    else
		write(lunout,'(/2a,3(1x,i12))') &
		    '*** subr. modal_aero_rename_acc_crs_init', &
		    'bad ipair_select_renamexf', ipair, itmpa
		call endrun( 'modal_aero_rename_acc_crs_init error' )
	    end if

	    do i = 1, ipair-1
		if (itmpa .eq. ipair_select_renamexf(i)) then
		    write(lunout,'(/2a/10(1x,i12))') &
			'*** subr. modal_aero_rename_acc_crs_init', &
			'duplicates in ipair_select_renamexf', &
			ipair_select_renamexf(1:ipair)
		    call endrun( 'modal_aero_rename_acc_crs_init error' )
		end if
	    end do

	    if ( (mfrm .ge. 1) .and. (mfrm .le. ntot_amode) .and. &
	         (mtoo .ge. 1) .and. (mtoo .le. ntot_amode) ) then
		npair_renamexf = ipair
		modefrm_renamexf(ipair) = mfrm
		modetoo_renamexf(ipair) = mtoo
	    else
		write(lunout,'(/2a,3(1x,i12))') &
		    '*** subr. modal_aero_rename_acc_crs_init', &
		    'bad mfrm or mtoo', ipair, mfrm, mtoo
		call endrun( 'modal_aero_rename_acc_crs_init error' )
	    end if
	end do ! ipair

	if (npair_renamexf .le. 0) then
	    write(lunout,'(/a/a,3(1x,i12))') &
		'*** subr. modal_aero_rename_acc_crs_init -- npair_renamexf = 0'
	    return
	end if


!
!   define species involved in each tail-xfer pairing
!	(include aerosol water)
!
	do 1900 ipair = 1, npair_renamexf
	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)
	ixferable_all_renamexf(ipair) = 1

	if (mfrm < 10) then
	    nch_mfrmid = 1
	else if (mfrm < 100) then
	    nch_mfrmid = 2
	else
	    nch_mfrmid = 3
	end if
	if (mtoo < 10) then
	    nch_mtooid = 1
	else if (mtoo < 100) then
	    nch_mtooid = 2
	else
	    nch_mtooid = 3
	end if

	nspec = 0
	do 1490 iqfrm = -1, nspec_amode(mfrm)
	    if (iqfrm .eq. -1) then
		lsfrma = numptr_amode(mfrm)
		lstooa = numptr_amode(mtoo)
		lsfrmc = numptrcw_amode(mfrm)
		lstooc = numptrcw_amode(mtoo)
	    else if (iqfrm .eq. 0) then
!   bypass transfer of aerosol water due to renaming
                goto 1490
!               lsfrma = lwaterptr_amode(mfrm)
!               lsfrmc = 0
!               lstooa = lwaterptr_amode(mtoo)
!               lstooc = 0
	    else
		lsfrma = lmassptr_amode(iqfrm,mfrm)
		lsfrmc = lmassptrcw_amode(iqfrm,mfrm)
		lstooa = 0
		lstooc = 0
	    end if

	    if ((lsfrma .lt. 1) .or. (lsfrma .gt. pcnst)) then
		write(lunout,9100) ipair, mfrm, iqfrm, lsfrma
		call endrun( 'modal_aero_rename_acc_crs_init error' )
	    end if
	    if (iqfrm .le. 0) goto 1430

	    if ((lsfrmc .lt. 1) .or. (lsfrmc .gt. pcnst)) then
		write(lunout,9102) ipair, mfrm, iqfrm, lsfrmc
		call endrun( 'modal_aero_rename_acc_crs_init error' )
	    end if

! find "too" species having same name (except for mode number) as the "frm" species
	    nch_lfrm = len(trim(cnst_name(lsfrma))) - nch_mfrmid
	    iqtooa = -99
	    do iq = 1, nspec_amode(mtoo)
		l = lmassptr_amode(iq,mtoo)
		if ((l .lt. 1) .or. (l .gt. pcnst)) cycle
		nch_ltoo = len(trim(cnst_name(l))) - nch_mtooid
		if ( cnst_name(lsfrma)(1:nch_lfrm) == &
		     cnst_name(l     )(1:nch_ltoo) ) then
		    lstooa = l
		    iqtooa = iq
		    exit
		end if
	    end do

	    nch_lfrm = len(trim(cnst_name_cw(lsfrmc))) - nch_mfrmid
	    iqtooc = -99
	    do iq = 1, nspec_amode(mtoo)
		l = lmassptrcw_amode(iq,mtoo)
		if ((l .lt. 1) .or. (l .gt. pcnst)) cycle
		nch_ltoo = len(trim(cnst_name_cw(l))) - nch_mtooid
		if ( cnst_name_cw(lsfrmc)(1:nch_lfrm) == &
		     cnst_name_cw(l     )(1:nch_ltoo) ) then
		    lstooc = l
		    iqtooc = iq
		    exit
		end if
	    end do

1430	    if ((lstooc .lt. 1) .or. (lstooc .gt. pcnst)) lstooc = 0
	    if ((lstooa .lt. 1) .or. (lstooa .gt. pcnst)) lstooa = 0

	    if ((lstooa .eq. 0) .or. (lstooc .eq. 0)) then
		if ( ( masterproc                                  ) .or. &
		     ( (lstooa .ne. 0) .or. (lstooc .ne. 0)        ) .or. &
		     ( ixferable_all_needed_renamexf(ipair) .gt. 0 ) ) then
		    if (lstooa .eq. 0) &
			write(lunout,9104) trim(cnst_name(lsfrma)), &
			    ipair, mfrm, iqfrm, lsfrma, iqtooa, lstooa
		    if (lstooc .eq. 0) &
			write(lunout,9106) trim(cnst_name_cw(lsfrmc)), &
			    ipair, mfrm, iqfrm, lsfrmc, iqtooc, lstooc
		end if
		if ((lstooa .ne. 0) .or. (lstooc .ne. 0)) then
		    write(lunout,9108)
		    call endrun( 'modal_aero_rename_acc_crs_init error' )
		end if
		if (ixferable_all_needed_renamexf(ipair) .gt. 0) then
		    write(lunout,9109)
		    call endrun( 'modal_aero_rename_acc_crs_init error' )
		end if
		ixferable_all_renamexf(ipair) = 0
		if (iqfrm .gt. 0) then
		    ixferable_a_renamexf(iqfrm,ipair) = 0
		    ixferable_c_renamexf(iqfrm,ipair) = 0
		end if
	    else
		nspec = nspec + 1
		lspecfrma_renamexf(nspec,ipair) = lsfrma
		lspectooa_renamexf(nspec,ipair) = lstooa
		lspecfrmc_renamexf(nspec,ipair) = lsfrmc
		lspectooc_renamexf(nspec,ipair) = lstooc
		if (iqfrm .gt. 0) then
		    ixferable_a_renamexf(iqfrm,ipair) = 1
		    ixferable_c_renamexf(iqfrm,ipair) = 1
		end if
	    end if
1490	continue

	nspecfrm_renamexf(ipair) = nspec
1900	continue

9100	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'lspecfrma out of range' /   &
      	'ipair, modefrm, ispecfrm, lspecfrma =', 4i6 )
9102	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'lspecfrmc out of range' /   &
      	'ipair, modefrm, ispecfrm, lspecfrmc =', 4i6 )
9104	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'lspectooa out of range for', 2x, a /   &
      	'ipair, modefrm, ispecfrm, lspecfrma, ispectoo, lspectooa =', 6i6 )
9106	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'lspectooc out of range for', 2x, a /   &
      	'ipair, modefrm, ispecfrm, lspecfrmc, ispectoo, lspectooc =', 6i6 )
9108	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'only one of lspectooa and lspectooc is out of range' )
9109	format( / '*** subr. modal_aero_rename_acc_crs_init' /   &
      	'all species must be xferable for this pair' )


!
!
!   initialize some working variables
!
!
	ido_mode_calcaa(:) = 0
	frelax = 27.0_r8

	do ipair = 1, npair_renamexf
	    mfrm = modefrm_renamexf(ipair)
	    mtoo = modetoo_renamexf(ipair)
	    ido_mode_calcaa(mfrm) = 1

	    factoraa(mfrm) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mfrm)**2))
	    factoraa(mtoo) = (pi/6._r8)*exp(4.5_r8*(alnsg_amode(mtoo)**2))
	    factoryy(mfrm) = sqrt( 0.5_r8 )/alnsg_amode(mfrm)

!   dryvol_smallest is a very small volume mixing ratio (m3-AP/kmol-air)
!   used for avoiding overflow.  it corresponds to dp = 1 nm
!   and number = 1e-5 #/mg-air ~= 1e-5 #/cm3-air
	    dryvol_smallest(mfrm) = 1.0e-25_r8
	    v2nlorlx(mfrm) = voltonumblo_amode(mfrm)*frelax
	    v2nhirlx(mfrm) = voltonumbhi_amode(mfrm)/frelax

	    factor_3alnsg2(ipair) = 3.0_r8 * (alnsg_amode(mfrm)**2)

	    dp_cut(ipair) = sqrt(   &
		dgnum_amode(mfrm)*exp(1.5_r8*(alnsg_amode(mfrm)**2)) *   &
		dgnum_amode(mtoo)*exp(1.5_r8*(alnsg_amode(mtoo)**2)) )
	    dp_xferall_thresh(ipair) = dgnum_amode(mtoo)
	    dp_xfernone_threshaa(ipair) = dgnum_amode(mfrm)

	    if ((mfrm == modeptr_accum) .and. (mtoo == modeptr_coarse)) then
		dp_cut(ipair)               = 1.0e-6_r8
		dp_xfernone_threshaa(ipair) = 0.9e-6_r8
		dp_xferall_thresh(ipair)    = 1.1e-6_r8
	    else if ((mfrm == modeptr_coarse) .and. (mtoo == modeptr_accum)) then
		dp_cut(ipair)               = 1.0e-6_r8
		dp_xfernone_threshaa(ipair) = 1.0e-6_r8
		dp_xferall_thresh(ipair)    = 0.9e-6_r8
	    end if

	    lndp_cut(ipair) = log( dp_cut(ipair) )
	    dp_belowcut(ipair) = 0.99_r8*dp_cut(ipair)
	end do


!
!   output results
!
	if ( masterproc ) then

	write(lunout,9310)
	write(lunout,'(a,1x,i12)') 'method_optbb_renamexf', method_optbb_renamexf

	do 2900 ipair = 1, npair_renamexf
	mfrm = modefrm_renamexf(ipair)
	mtoo = modetoo_renamexf(ipair)
	write(lunout,9320) ipair, mfrm, mtoo, &
	    igrow_shrink_renamexf(ipair), ixferable_all_renamexf(ipair)

	do iq = 1, nspecfrm_renamexf(ipair)
	    lsfrma = lspecfrma_renamexf(iq,ipair)
	    lstooa = lspectooa_renamexf(iq,ipair)
	    lsfrmc = lspecfrmc_renamexf(iq,ipair)
	    lstooc = lspectooc_renamexf(iq,ipair)
	    if (lstooa .gt. 0) then
		write(lunout,9330) lsfrma, cnst_name(lsfrma),   &
				   lstooa, cnst_name(lstooa)
	    else
		write(lunout,9340) lsfrma, cnst_name(lsfrma)
	    end if
	    if (lstooc .gt. 0) then
		write(lunout,9330) lsfrmc, cnst_name_cw(lsfrmc),   &
				   lstooc, cnst_name_cw(lstooc)
	    else if (lsfrmc .gt. 0) then
		write(lunout,9340) lsfrmc, cnst_name_cw(lsfrmc)
	    else
		write(lunout,9350)
	    end if
	end do

	if (igrow_shrink_renamexf(ipair) > 0) then
	write(lunout,'(5x,a,1p,2e12.3)') 'mfrm dgnum, dgnumhi ', &
		dgnum_amode(mfrm), dgnumhi_amode(mfrm)
	write(lunout,'(5x,a,1p,2e12.3)') 'mtoo dgnum, dgnumlo ', &
		dgnum_amode(mtoo), dgnumlo_amode(mtoo)
	else
	write(lunout,'(5x,a,1p,2e12.3)') 'mfrm dgnum, dgnumlo ', &
		dgnum_amode(mfrm), dgnumlo_amode(mfrm)
	write(lunout,'(5x,a,1p,2e12.3)') 'mtoo dgnum, dgnumhi ', &
		dgnum_amode(mtoo), dgnumhi_amode(mtoo)
	end if

	write(lunout,'(5x,a,1p,2e12.3)') 'dp_cut              ', &
		dp_cut(ipair)
	write(lunout,'(5x,a,1p,2e12.3)') 'dp_xfernone_threshaa', &
		dp_xfernone_threshaa(ipair)
	write(lunout,'(5x,a,1p,2e12.3)') 'dp_xferall_thresh   ', &
		dp_xferall_thresh(ipair)

2900	continue
	write(lunout,*)

	end if ! ( masterproc )

9310	format( / 'subr. modal_aero_rename_acc_crs_init' )
9320	format( / 'pair', i3, 5x, 'mode', i3, ' ---> mode', i3, &
	        5x, 'igrow_shrink', i3, 5x, 'ixferable_all', i3 )
9330	format( 5x, 'spec', i3, '=', a, ' ---> spec', i3, '=', a )
9340	format( 5x, 'spec', i3, '=', a, ' ---> LOSS' )
9350	format( 5x, 'no corresponding activated species' )


	return
	end subroutine modal_aero_rename_acc_crs_init

!----------------------------------------------------------------------

   end module modal_aero_rename
