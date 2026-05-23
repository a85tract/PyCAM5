!===============================================================================
! Age of air test tracers
! provides dissipation rate and surface fluxes for diagnostic constituents
!===============================================================================

module aoa_tracers

  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver
  use constituents, only: pcnst, cnst_add, cnst_name, cnst_longname
  use cam_logfile,  only: iulog
  use ref_pres,     only: pref_mid_norm

  implicit none
  private
  save

  ! Public interfaces
  public :: aoa_tracers_register         ! register constituents
  public :: aoa_tracers_implements_cnst  ! true if named constituent is implemented by this package
  public :: aoa_tracers_init_cnst        ! initialize constituent field
  public :: aoa_tracers_init             ! initialize history fields, datasets
  public :: aoa_tracers_timestep_init    ! place to perform per timestep initialization
  public :: aoa_tracers_timestep_tend    ! calculate tendencies
  public :: aoa_tracers_readnl           ! read namelist options

  ! Private module data

  integer, parameter :: ncnst=4  ! number of constituents implemented by this module

  ! constituent names
  character(len=8), parameter :: c_names(ncnst) = (/'AOA1', 'AOA2', 'HORZ', 'VERT'/)

  ! constituent source/sink names
  character(len=8), parameter :: src_names(ncnst) = (/'AOA1SRC', 'AOA2SRC', 'HORZSRC', 'VERTSRC'/)

  integer :: ifirst ! global index of first constituent
  integer :: ixaoa1 ! global index for AOA1 tracer
  integer :: ixaoa2 ! global index for AOA2 tracer
  integer :: ixht   ! global index for HORZ tracer
  integer :: ixvt   ! global index for VERT tracer

  ! Data from namelist variables
  logical :: aoa_tracers_flag  = .false.    ! true => turn on test tracer code, namelist variable
  logical :: aoa_read_from_ic_file = .true. ! true => tracers initialized from IC file
  logical :: use_native_impl = .false.
  logical :: impl_selected = .false.
  logical :: use_native_tstep_init_impl = .false.
  logical :: tstep_init_impl_selected = .false.
  logical :: aoa_tracers_register_logged = .false.
  logical :: aoa_tracers_implements_cnst_logged = .false.
  logical :: aoa_tracers_init_logged = .false.

  interface
     function aoa_tracers_flag_codon(flag_c) result(out_c) bind(c, name="aoa_tracers_flag_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: out_c
     end function aoa_tracers_flag_codon
     function aoa_tracers_implements_cnst_codon(flag_c) result(out_c) &
          bind(c, name="aoa_tracers_implements_cnst_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: out_c
     end function aoa_tracers_implements_cnst_codon
  end interface
  
  real(r8),  parameter ::  treldays = 15._r8
  real(r8),  parameter ::  vert_offset = 10._r8

  ! 15-days used for diagnostic of transport circulation and K-tensors  
  ! relaxation (in the original papers PM-1987 and YSGD-2000) => Zonal Mean 
  ! to evaluate eddy-fluxes for 2D-diagnostics, here relaxation to the GLOBAL MEAN  IC
  ! it may help to keep gradients but will rule-out 2D-transport diagnostics
  ! in km  to avoid negative values of  vertical tracers
  ! VERT(k) = -7._r8*alog(hyam(k)+hybm(k)) + vert_offset
  
  ! PM-1987:
  ! Plumb, R. A., and J. D. Mahlman (1987), The zonally averaged transport
  ! characteristics of the GFDL general circulation/transport model,
  ! J. Atmos.Sci.,44, 298–327

  ! YSGD-2000:
  ! Yudin, Valery A., Sergey P. Smyshlyaev, Marvin A. Geller, Victor L. Dvortsov, 2000: 
  ! Transport Diagnostics of GCMs and Implications for 2D Chemistry-Transport Model of 
  ! Troposphere and Stratosphere. J. Atmos. Sci., 57, 673–699.
  ! doi: http://dx.doi.org/10.1175/1520-0469(2000)057<0673:TDOGAI>2.0.CO;2

  real(r8), target :: qrel_vert(pver)  ! = -7._r8*log(pref_mid_norm(k)) + vert_offset

!===============================================================================
contains
!===============================================================================

!================================================================================
  subroutine aoa_tracers_readnl(nlfile)

    use namelist_utils,     only: find_group_name
    use units,              only: getunit, freeunit
    use mpishorthand
    use cam_abortutils,     only: endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aoa_tracers_readnl'


    namelist /aoa_tracers_nl/ aoa_tracers_flag, aoa_read_from_ic_file

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aoa_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aoa_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(aoa_tracers_flag, 1, mpilog,  0, mpicom)
    call mpibcast(aoa_read_from_ic_file, 1, mpilog,  0, mpicom)
#endif

  endsubroutine aoa_tracers_readnl

!================================================================================

  subroutine aoa_tracers_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents
    ! 
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    use iso_c_binding, only: c_int64_t
    !-----------------------------------------------------------------------
    integer(c_int64_t) :: active_c

    call aoa_tracers_select_impl()
    if (use_native_impl) then
       if (.not. aoa_tracers_flag) return
    else
       active_c = aoa_tracers_flag_codon(merge(1_c_int64_t, 0_c_int64_t, aoa_tracers_flag))
       call aoa_tracers_log_direct(aoa_tracers_register_logged, 'aoa_tracers_register direct = codon')
       if (active_c == 0_c_int64_t) return
    end if

    call cnst_add(c_names(1), mwdry, cpair, 0._r8, ixaoa1, readiv=aoa_read_from_ic_file, &
                  longname='Age-of_air tracer 1')
    ifirst = ixaoa1
    call cnst_add(c_names(2), mwdry, cpair, 0._r8, ixaoa2, readiv=aoa_read_from_ic_file, &
                  longname='Age-of_air tracer 2')
    call cnst_add(c_names(3), mwdry, cpair, 1._r8, ixht,   readiv=aoa_read_from_ic_file, &
                  longname='horizontal tracer')
    call cnst_add(c_names(4), mwdry, cpair, 0._r8, ixvt,   readiv=aoa_read_from_ic_file, &
                  longname='vertical tracer')

  end subroutine aoa_tracers_register

!===============================================================================

  function aoa_tracers_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    !-----------------------------------------------------------------------

    use iso_c_binding, only: c_int64_t

    character(len=*), intent(in) :: name   ! constituent name
    logical :: aoa_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    integer(c_int64_t) :: active_c, out_c
    !-----------------------------------------------------------------------

    aoa_tracers_implements_cnst = .false.

    call aoa_tracers_select_impl()
    if (use_native_impl) then
       if (.not. aoa_tracers_flag) return
    else
       active_c = aoa_tracers_flag_codon(merge(1_c_int64_t, 0_c_int64_t, aoa_tracers_flag))
       if (active_c == 0_c_int64_t) return
    end if

    do m = 1, ncnst
       if (name == c_names(m)) then
          if (use_native_impl) then
             aoa_tracers_implements_cnst = .true.
          else
             out_c = aoa_tracers_implements_cnst_codon(1_c_int64_t)
             aoa_tracers_implements_cnst = out_c /= 0_c_int64_t
             call aoa_tracers_log_direct(aoa_tracers_implements_cnst_logged, &
                  'aoa_tracers_implements_cnst direct = codon')
          end if
          return
       end if
    end do

  end function aoa_tracers_implements_cnst

!===============================================================================

  subroutine aoa_tracers_init_cnst(name, q, gcid)

    !----------------------------------------------------------------------- 
    !
    ! Purpose: initialize test tracers mixing ratio fields 
    !  This subroutine is called at the beginning of an initial run ONLY
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in)  :: name
    real(r8),         intent(out) :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id

    integer :: m
    !-----------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, q, gcid)
       endif
    end do

  end subroutine aoa_tracers_init_cnst

!===============================================================================

  subroutine aoa_tracers_init

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize age of air constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default, phys_decomp
    use iso_c_binding, only: c_int64_t

    integer :: m, mm, k
    integer(c_int64_t) :: active_c
    !-----------------------------------------------------------------------

    call aoa_tracers_select_impl()
    if (use_native_impl) then
       if (.not. aoa_tracers_flag) return
    else
       active_c = aoa_tracers_flag_codon(merge(1_c_int64_t, 0_c_int64_t, aoa_tracers_flag))
       call aoa_tracers_log_direct(aoa_tracers_init_logged, 'aoa_tracers_init direct = codon')
       if (active_c == 0_c_int64_t) return
    end if

    ! Set names of tendencies and declare them as history variables

    do m = 1, ncnst
       mm = ifirst+m-1
       call addfld (cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm), phys_decomp)
       call addfld (src_names(m),  'kg/kg/s ', pver, 'A', trim(cnst_name(mm))//' source/sink', phys_decomp)

       call add_default (cnst_name(mm), 1, ' ')
       call add_default (src_names(m),  1, ' ')
    end do

    do k = 1,pver
       qrel_vert(k) = -7._r8*log(pref_mid_norm(k)) + vert_offset
    enddo

  end subroutine aoa_tracers_init

!===============================================================================

  subroutine aoa_tracers_timestep_init( phys_state )
    !-----------------------------------------------------------------------
    ! Provides a place to reinitialize diagnostic constituents HORZ and VERT
    !-----------------------------------------------------------------------

    use iso_c_binding, only: c_int64_t, c_loc, c_ptr
    use time_manager,   only: get_curr_date
    use ppgrid,         only: begchunk, endchunk
    use physics_types,  only: physics_state

    type(physics_state), target, intent(inout), dimension(begchunk:endchunk), optional :: phys_state    


    integer c, ncol
    integer yr, mon, day, tod
    interface
       subroutine aoa_tracers_tstep_init_codon(ncol_c, pcols_c, pver_c, ixht_c, ixvt_c, &
            qrel_vert_p, state_lat_p, state_q_p) bind(c, name="aoa_tracers_tstep_init_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, ixht_c, ixvt_c
         type(c_ptr), value :: qrel_vert_p, state_lat_p, state_q_p
       end subroutine aoa_tracers_tstep_init_codon
    end interface
    !--------------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    call aoa_tracers_tstep_init_select_impl()

    if (use_native_tstep_init_impl) then
       call aoa_tracers_timestep_init_native(phys_state)
       return
    end if

    call get_curr_date (yr,mon,day,tod)

    if ( day == 1 .and. tod == 0) then
       if (masterproc) then
         write(iulog,*) 'AGE_OF_AIR_CONSTITUENTS: RE-INITIALIZING HORZ/VERT CONSTITUENTS'
       endif

       do c = begchunk, endchunk
          ncol = phys_state(c)%ncol
          call aoa_tracers_tstep_init_codon( &
               int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
               int(ixht, c_int64_t), int(ixvt, c_int64_t), c_loc(qrel_vert), &
               c_loc(phys_state(c)%lat), c_loc(phys_state(c)%q) &
          )
       end do

    end if

  end subroutine aoa_tracers_timestep_init

!===============================================================================

  subroutine aoa_tracers_tstep_init_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (tstep_init_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('AOA_TRACERS_TSTEP_INIT_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_tstep_init_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_tstep_init_impl = .false.
    end if

    tstep_init_impl_selected = .true.

    if (masterproc) then
       if (use_native_tstep_init_impl) then
          write(iulog,*) 'aoa_tracers_timestep_init implementation = native'
       else
          write(iulog,*) 'aoa_tracers_timestep_init implementation = codon'
       end if
    end if

  end subroutine aoa_tracers_tstep_init_select_impl

!===============================================================================

  subroutine aoa_tracers_timestep_init_native( phys_state )
    !-----------------------------------------------------------------------
    ! Provides a place to reinitialize diagnostic constituents HORZ and VERT
    !-----------------------------------------------------------------------

    use time_manager,   only: get_curr_date
    use ppgrid,         only: begchunk, endchunk
    use physics_types,  only: physics_state

    type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state

    integer c, i, k, ncol
    integer yr, mon, day, tod
    !--------------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    call get_curr_date (yr,mon,day,tod)

    if ( day == 1 .and. tod == 0) then
       if (masterproc) then
         write(iulog,*) 'AGE_OF_AIR_CONSTITUENTS: RE-INITIALIZING HORZ/VERT CONSTITUENTS'
       endif

       do c = begchunk, endchunk
          ncol = phys_state(c)%ncol
          do k = 1, pver
             do i = 1, ncol
                phys_state(c)%q(i,k,ixht) = 2._r8 + sin(phys_state(c)%lat(i))
                phys_state(c)%q(i,k,ixvt) = qrel_vert(k)
             end do
          end do
       end do

    end if

  end subroutine aoa_tracers_timestep_init_native

!===============================================================================

  subroutine aoa_tracers_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('AOA_TRACERS_IMPL', value=impl_name, length=n, status=status)

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
          write(iulog,*) 'aoa_tracers implementation = native'
       else
          write(iulog,*) 'aoa_tracers implementation = codon'
       end if
    end if

  end subroutine aoa_tracers_select_impl

!===============================================================================

  subroutine aoa_tracers_log_direct(logged, proof_line)

    logical, intent(inout) :: logged
    character(len=*), intent(in) :: proof_line

    if (logged) return
    logged = .true.

    if (masterproc) then
       write(iulog,'(A)') trim(proof_line)
       call flush(iulog)
    end if

  end subroutine aoa_tracers_log_direct

!===============================================================================

  subroutine aoa_tracers_timestep_tend(state, ptend, cflx, landfrac, dt)

    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use cam_history,   only: outfld
    use time_manager,  only: get_nstep

    ! Arguments
    type(physics_state), target, intent(in)    :: state              ! state variables
    type(physics_ptend), target, intent(out)   :: ptend              ! package tendencies
    real(r8), target, intent(inout) :: cflx(pcols,pcnst)  ! Surface constituent flux (kg/m^2/s)
    real(r8), target, intent(in)    :: landfrac(pcols)    ! Land fraction
    real(r8),            intent(in)    :: dt                 ! timestep

    integer :: nstep                          ! current timestep number
    integer :: lchnk                          ! chunk identifier
    logical  :: lq(pcnst)
    interface
       subroutine aoa_tracers_timestep_tend_codon(ncol_c, pcols_c, pver_c, pcnst_c, psetcols_c, &
            ixaoa1_c, ixaoa2_c, ixht_c, ixvt_c, nstep_c, dt_c, qrel_vert_p, state_lat_p, &
            state_q_p, ptend_q_p, cflx_p, landfrac_p) bind(c, name="aoa_tracers_timestep_tend_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pcnst_c, psetcols_c
         integer(c_int64_t), value :: ixaoa1_c, ixaoa2_c, ixht_c, ixvt_c, nstep_c
         real(c_double), value :: dt_c
         type(c_ptr), value :: qrel_vert_p, state_lat_p, state_q_p, ptend_q_p, cflx_p, landfrac_p
       end subroutine aoa_tracers_timestep_tend_codon
    end interface

    if (.not. aoa_tracers_flag) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
    end if

    call aoa_tracers_select_impl()

    if (use_native_impl) then
       call aoa_tracers_timestep_tend_native(state, ptend, cflx, landfrac, dt)
       return
    end if

    lq(:)      = .FALSE.
    lq(ixaoa1) = .TRUE.
    lq(ixaoa2) = .TRUE.
    lq(ixht)   = .TRUE.
    lq(ixvt)   = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, 'aoa_tracers', lq=lq)

    nstep = get_nstep()
    lchnk = state%lchnk

    call aoa_tracers_timestep_tend_codon( &
         int(state%ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
         int(pcnst, c_int64_t), int(state%psetcols, c_int64_t), &
         int(ixaoa1, c_int64_t), int(ixaoa2, c_int64_t), int(ixht, c_int64_t), int(ixvt, c_int64_t), &
         int(nstep, c_int64_t), real(dt, c_double), c_loc(qrel_vert), c_loc(state%lat), &
         c_loc(state%q), c_loc(ptend%q), c_loc(cflx), c_loc(landfrac) &
    )

    ! record tendencies on history files
    call outfld (src_names(1), ptend%q(:,:,ixaoa1), pcols, lchnk)
    call outfld (src_names(2), ptend%q(:,:,ixaoa2), pcols, lchnk)
    call outfld (src_names(3), ptend%q(:,:,ixht),   pcols, lchnk)
    call outfld (src_names(4), ptend%q(:,:,ixvt),   pcols, lchnk)

  end subroutine aoa_tracers_timestep_tend

!===============================================================================

  subroutine aoa_tracers_timestep_tend_native(state, ptend, cflx, landfrac, dt)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use phys_grid,     only: get_rlat_all_p , get_lat_all_p
    use cam_history,   only: outfld
    use time_manager,  only: get_nstep

    ! Arguments
    type(physics_state), intent(in)    :: state              ! state variables
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(inout) :: cflx(pcols,pcnst)  ! Surface constituent flux (kg/m^2/s)
    real(r8),            intent(in)    :: landfrac(pcols)    ! Land fraction
    real(r8),            intent(in)    :: dt                 ! timestep

    !----------------- Local workspace-------------------------------

    integer :: i, k
    integer :: lchnk                          ! chunk identifier
    integer :: ncol                           ! no. of column in chunk
    integer :: nstep                          ! current timestep number
    real(r8) :: qrel                          ! value to be relaxed to
    real(r8) :: xhorz                         ! updated value of HORZ
    real(r8) :: xvert                         ! updated value of VERT
    logical  :: lq(pcnst)
    real(r8) :: teul                          ! relaxation in  1/sec*dt/2 = k*dt/2
    real(r8) :: wimp                          !     1./(1.+ k*dt/2)
    real(r8) :: wsrc                          !  teul*wimp    
    !------------------------------------------------------------------

    teul = .5_r8*dt/(86400._r8 * treldays)   ! 1/2 for the semi-implicit scheme if dt=time step
    wimp = 1._r8/(1._r8 +teul)
    wsrc = teul*wimp

    if (.not. aoa_tracers_flag) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
    end if

    lq(:)      = .FALSE.
    lq(ixaoa1) = .TRUE.
    lq(ixaoa2) = .TRUE.
    lq(ixht)   = .TRUE.
    lq(ixvt)   = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, 'aoa_tracers', lq=lq)

    nstep = get_nstep()
    lchnk = state%lchnk
    ncol  = state%ncol

    do k = 1, pver
       do i = 1, ncol

          ! AOA1
          ptend%q(i,k,ixaoa1) = 0.0_r8

          ! AOA2
          ptend%q(i,k,ixaoa2) = 0.0_r8

          ! HORZ
          qrel              = 2._r8 + sin(state%lat(i))          ! qrel  should zonal mean
          xhorz             = state%q(i,k,ixht)*wimp + wsrc*qrel ! Xnew = weight*3D-tracer + (1.-weight)*1D-tracer
          ptend%q(i,k,ixht) = (xhorz - state%q(i,k,ixht)) / dt   ! Xnew = weight*3D-tracer + (1.-weight)*2D-tracer  zonal mean
                                                                 !  Can be still used .... to diagnose fluxes OT-tracers
          ! VERT
          qrel              = qrel_vert(k)                       ! qrel  should zonal mean
          xvert             = wimp*state%q(i,k,ixvt) + wsrc*qrel
          ptend%q(i,k,ixvt) = (xvert - state%q(i,k,ixvt)) / dt

       end do
    end do

    ! record tendencies on history files
    call outfld (src_names(1), ptend%q(:,:,ixaoa1), pcols, lchnk)
    call outfld (src_names(2), ptend%q(:,:,ixaoa2), pcols, lchnk)
    call outfld (src_names(3), ptend%q(:,:,ixht),   pcols, lchnk)
    call outfld (src_names(4), ptend%q(:,:,ixvt),   pcols, lchnk)

    ! Set tracer fluxes
    do i = 1, ncol

       ! AOA1
       cflx(i,ixaoa1) = 1.e-6_r8

       ! AOA2
       if (landfrac(i) .eq. 1._r8  .and.  state%lat(i) .gt. 0.35_r8) then
          cflx(i,ixaoa2) = 1.e-6_r8 + 1e-6_r8*0.0434_r8*real(nstep,r8)*dt/(86400._r8*365._r8)
       else
          cflx(i,ixaoa2) = 0._r8
       endif

       ! HORZ
       cflx(i,ixht) = 0._r8

       ! VERT
       cflx(i,ixvt) = 0._r8

    end do

  end subroutine aoa_tracers_timestep_tend_native

!===========================================================================

  subroutine init_cnst_3d(m, q, gcid)

    use dyn_grid, only : get_horiz_grid_d, get_horiz_grid_dim_d
    use dycore,   only : dycore_is

    integer,  intent(in)  :: m       ! global constituent index
    real(r8), intent(out) :: q(:,:)  ! kg tracer/kg dry air (gcol,plev)
    integer,  intent(in)  :: gcid(:) ! global column id

    real(r8), allocatable :: lat(:)
    integer :: plon, plat, ngcols
    integer :: j, k, gsize
    !-----------------------------------------------------------------------

    if (masterproc) write(iulog,*) 'AGE-OF-AIR CONSTITUENTS: INITIALIZING ',cnst_name(m),m

    if (m == ixaoa1) then

       q(:,:) = 0.0_r8

    else if (m == ixaoa2) then

       q(:,:) = 0.0_r8

    else if (m == ixht) then

       call get_horiz_grid_dim_d( plon, plat )
       ngcols = plon*plat
       gsize = size(gcid)
       allocate(lat(ngcols))
       call get_horiz_grid_d(ngcols,clat_d_out=lat)
       do j = 1, gsize
          q(j,:) = 2._r8 + sin(lat(gcid(j)))
       end do
       deallocate(lat)

    else if (m == ixvt) then

       do k = 1, pver
          do j = 1, size(q,1)
             q(j,k) = qrel_vert(k)
          end do
       end do

    end if

  end subroutine init_cnst_3d

!=====================================================================


end module aoa_tracers
