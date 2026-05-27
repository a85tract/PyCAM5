!======================================================================
! This is an interface for 3 test tracers: tr1,tr2,tr3
!
! This uses the tracers_suite module to initialize 
!   mixing ratios, fluxes and calculate tendencies.
!   Details of calling tree below. All of the detailed information about the
!  tracers should be store in the suite file, including the number & names of tracers. 
!
! Author B. Eaton
! History  D. Bundy, June 2003 modified to the format of physics interface
!        
!
!---------------------------------------------------------------
!
!  ------------  calling tree --------------
!  Register the tracers as advected fields, pass names to model
!  initindx.F90:			call tracers_register()
!
!  Initialize the tracer mixing ratio field
!  inidat.F90:read_inidat
!  	-> tracers.F90: tracers_init_cnst 
!  		-> tracers_suite.F90:init_cnst_tr
!
!  Initialize data set, things that need to be done at the beginning of a 
!  run (whether restart or initial)
!  inti.F90
!  	-> tracers.F90: tracers_init
!  		-> tracers_suite.F90:init_tr
!  		-> addfld/add default for surface flux (SF)
!
!  Timestepping:
!  advnce.F90
!  	-> tracers_timestep_init
!  		-> tracers_suite.F90:timestep_init_tr
!
!  tphysac.F90
!  	-> tracers_timestep_tend
!  		-> tracers_suite.F90:flux_tr
!  		-> tracers_suite.F90:tend_tr
!
!======================================================================


module tracers

  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile,  only: iulog
  use spmd_utils,   only: masterproc

  implicit none
  private
  save

! Public interfaces
  public tracers_register                  ! register constituent
  public tracers_implements_cnst           ! true if named constituent is implemented by this package
  public tracers_init_cnst                 ! initialize constituent field
  public tracers_init                      ! initialize history fields, datasets
  public tracers_timestep_tend             ! calculate tendencies
  public tracers_timestep_init             ! interpolate dataset for constituent each timestep

! Data from namelist variables
  logical, public :: tracers_flag  = .false.     ! true => turn on test tracer code, namelist variable

! Private module data

  integer :: trac_ncnst                    ! total number of test tracers
  integer :: ixtrct=-999                   ! index of 1st constituent
  logical :: debug = .false.
  logical :: use_native_impl = .false.
  logical :: impl_selected = .false.
  logical :: use_native_tstep_init_impl = .false.
  logical :: tstep_init_impl_selected = .false.
  logical :: tracers_register_logged = .false.
  logical :: tracers_implements_cnst_logged = .false.
  logical :: tracers_init_logged = .false.
  logical :: tracers_timestep_init_logged = .false.
  logical :: tracers_timestep_tend_logged = .false.

  interface
     subroutine tracers_timestep_init_codon() bind(c, name="tracers_timestep_init_codon")
     end subroutine tracers_timestep_init_codon
     function tracers_flag_codon(flag_c) result(out_c) bind(c, name="tracers_flag_codon")
        use iso_c_binding, only: c_int64_t
        integer(c_int64_t), value :: flag_c
        integer(c_int64_t) :: out_c
     end function tracers_flag_codon
     function tracers_implements_cnst_codon(flag_c, name_len_c, name_ascii_p, ncnst_c) result(out_c) &
          bind(c, name="tracers_implements_cnst_codon")
        use iso_c_binding, only: c_int64_t, c_ptr
        integer(c_int64_t), value :: flag_c, name_len_c, ncnst_c
        type(c_ptr), value :: name_ascii_p
        integer(c_int64_t) :: out_c
     end function tracers_implements_cnst_codon
     function tracers_register_codon(flag_c) result(out_c) bind(c, name="tracers_register_codon")
        use iso_c_binding, only: c_int64_t
        integer(c_int64_t), value :: flag_c
        integer(c_int64_t) :: out_c
     end function tracers_register_codon
     function tracers_init_codon(flag_c) result(out_c) bind(c, name="tracers_init_codon")
        use iso_c_binding, only: c_int64_t
        integer(c_int64_t), value :: flag_c
        integer(c_int64_t) :: out_c
     end function tracers_init_codon
  end interface
  
contains
!======================================================================
subroutine tracers_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TRACERS_IMPL', value=impl_name, length=n, status=status)

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
        write(iulog,*) 'tracers implementation = native'
     else
        write(iulog,*) 'tracers implementation = codon'
     end if
  end if

end subroutine tracers_select_impl
!======================================================================
subroutine tracers_tstep_init_select_impl()

  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (tstep_init_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('TRACERS_TSTEP_INIT_IMPL', value=impl_name, length=n, status=status)

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
        write(iulog,*) 'tracers_timestep_init implementation = native'
     else
        write(iulog,*) 'tracers_timestep_init implementation = codon'
     end if
  end if

end subroutine tracers_tstep_init_select_impl
!======================================================================
subroutine tracers_log_direct(logged, proof_line)

  logical, intent(inout) :: logged
  character(len=*), intent(in) :: proof_line

  if (logged) return
  logged = .true.

  if (masterproc) then
     write(iulog,'(A)') trim(proof_line)
     call flush(iulog)
  end if

end subroutine tracers_log_direct
!======================================================================
subroutine tracers_register
!----------------------------------------------------------------------- 
!
! Purpose: register advected tracers. Called by initindx.F90
!  The registration lets the model know what the tracer names are
!  and returns the index number ixtrct for the constituent array
! 
! Author: D. Bundy
!-----------------------------------------------------------------------

   use physconst,    only: mwdry, cpair
   use constituents, only: cnst_add, cnst_num_avail
   use tracers_suite, only: get_tracer_name
   use iso_c_binding, only: c_int64_t
   
   implicit none
!---------------------------Local workspace-----------------------------
   integer :: mm,m                                 ! dummy
   character(len=8) :: name   ! constituent name
   real(r8) minc
   integer(c_int64_t) :: active_c

!-----------------------------------------------------------------------
   call tracers_select_impl()
   if (use_native_impl) then
      if (.not. tracers_flag) return
   else
      active_c = tracers_register_codon(merge(1_c_int64_t, 0_c_int64_t, tracers_flag))
      call tracers_log_direct(tracers_register_logged, 'tracers_register direct = codon')
      if (active_c == 0_c_int64_t) return
   end if

      minc = 0        ! min mixing ratio (normal setting)
      minc = -1.e36_r8   ! min mixing ratio (disable qneg3)
      
      ! Set the number of test tracers equal to the number of slots available
      ! in the constituent array
      trac_ncnst = cnst_num_avail()

      do m = 1,trac_ncnst 
         name = get_tracer_name(m)  ! get name from suite file
         
         ! add constituent name to list of advected, save index number ixtrct
         call cnst_add(name, mwdry, cpair, minc, mm, &  
              readiv=.false.,mixtype='dry')
         if ( m .eq. 1 ) ixtrct = mm  ! save index number of first tracer
         
      end do

end subroutine tracers_register
!======================================================================

function tracers_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------

  use tracers_suite, only: get_tracer_name
  use iso_c_binding, only: c_int64_t, c_loc
  
  implicit none
!-----------------------------Arguments---------------------------------
  
  character(len=*), intent(in) :: name   ! constituent name
  logical :: tracers_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
   integer :: i, m
   integer(c_int64_t) :: active_c, out_c
   integer(c_int64_t), target :: name_ascii(max(1, len(name)))
!-----------------------------------------------------------------------

   tracers_implements_cnst = .false.
   call tracers_select_impl()
   if (use_native_impl) then
      if (.not. tracers_flag) return
   else
      active_c = tracers_flag_codon(merge(1_c_int64_t, 0_c_int64_t, tracers_flag))
      do i = 1, len(name)
         name_ascii(i) = int(iachar(name(i:i)), c_int64_t)
      end do
      out_c = tracers_implements_cnst_codon(active_c, int(len(name), c_int64_t), &
           c_loc(name_ascii(1)), int(trac_ncnst, c_int64_t))
      tracers_implements_cnst = out_c /= 0_c_int64_t
      call tracers_log_direct(tracers_implements_cnst_logged, 'tracers_implements_cnst direct = codon')
      return
   end if

   do m = 1, trac_ncnst
      if (name == get_tracer_name(m)) then
         tracers_implements_cnst = .true.
         return
      end if
   end do
end function tracers_implements_cnst

!===============================================================================
subroutine tracers_init_cnst(name, q, gcid)

!----------------------------------------------------------------------- 
!
! Purpose: initialize test tracers mixing ratio fields 
!  This subroutine is called at the beginning of an initial run ONLY
!
!-----------------------------------------------------------------------

  use tracers_suite,   only: init_cnst_tr, get_tracer_name

  implicit none

  character(len=*), intent(in) :: name
  real(r8), intent(out), dimension(:,:) :: q    ! kg tracer/kg dry air (gcol,plev)
  integer,  intent(in)                  :: gcid(:)  ! global column id
! Local
  integer m
  if ( tracers_flag ) then 
     do m = 1, trac_ncnst
        if (name ==  get_tracer_name(m))  then
           call init_cnst_tr(m,q, gcid)
        endif
     end do
  end if

end subroutine tracers_init_cnst

!===============================================================================
subroutine tracers_init

!----------------------------------------------------------------------- 
!
! Purpose: declare history variables, initialize data sets
!  This subroutine is called at the beginning of an initial or restart run
!
!-----------------------------------------------------------------------

   use tracers_suite,   only: init_tr, get_tracer_name
   use cam_history,     only: addfld, add_default, phys_decomp
   use ppgrid,          only: pver
   use constituents,    only: cnst_get_ind, cnst_name, cnst_longname, sflxnam
   use iso_c_binding,   only: c_int64_t

   ! Local
   integer m, mm
   character(len=8) :: name   ! constituent name
   integer(c_int64_t) :: active_c

   call tracers_select_impl()
   if (use_native_impl) then
      if (.not. tracers_flag) return
   else
      active_c = tracers_init_codon(merge(1_c_int64_t, 0_c_int64_t, tracers_flag))
      call tracers_log_direct(tracers_init_logged, 'tracers_init direct = codon')
      if (active_c == 0_c_int64_t) return
   end if
     
      do m = 1,trac_ncnst 
         name = get_tracer_name(m)
         call cnst_get_ind(name, mm)
         call addfld (cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm), phys_decomp)
         call addfld (sflxnam(mm),   'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)

         call add_default (cnst_name(mm), 1, ' ')
         call add_default (sflxnam(mm),   1, ' ')
      end do
     
      ! initialize datasets, etc, needed for constituents.
      call init_tr  

end subroutine tracers_init

!======================================================================

subroutine tracers_timestep_init( phys_state )
!----------------------------------------------------------------------- 
!
! Purpose: At the beginning of a timestep, there are some things to do
! that just the masterproc should do. This currently just interpolates
! the emissions boundary data set to the current time step.
!
!-----------------------------------------------------------------------

  use tracers_suite, only: timestep_init_tr

  ! phys_state argument is unused in this version
  use ppgrid,         only: begchunk, endchunk
  use physics_types,  only: physics_state
  type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state    
!-----------------------------------------------------------------------

  call tracers_tstep_init_select_impl()

  if (use_native_tstep_init_impl) then
     call tracers_timestep_init_native(phys_state)
     return
  end if

  if (tracers_flag) then
     call tracers_log_direct(tracers_timestep_init_logged, &
          'tracers_timestep_init direct = codon control shell; native enabled timestep_init_tr island')
     call tracers_timestep_init_native(phys_state)
     return
  end if

  call tracers_timestep_init_codon()
  call tracers_log_direct(tracers_timestep_init_logged, &
       'tracers_timestep_init direct = codon flag-off no-op')

end subroutine tracers_timestep_init

!======================================================================

subroutine tracers_timestep_init_native( phys_state )
!----------------------------------------------------------------------- 
!
! Purpose: At the beginning of a timestep, there are some things to do
! that just the masterproc should do. This currently just interpolates
! the emissions boundary data set to the current time step.
!
!-----------------------------------------------------------------------

  use tracers_suite, only: timestep_init_tr

  ! phys_state argument is unused in this version
  use ppgrid,         only: begchunk, endchunk
  use physics_types,  only: physics_state
  type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state    
!-----------------------------------------------------------------------

  if ( tracers_flag ) then 
     
     call timestep_init_tr
     
     if (debug) write(iulog,*)'tracers_timestep_init done'
  endif

end subroutine tracers_timestep_init_native

!======================================================================

subroutine tracers_timestep_tend(state, ptend, cflx, landfrac, deltat)

!----------------------------------------------------------------------- 
!
! Purpose: During the timestep, compute test tracer mixing ratio 
! tendencies and surface fluxes.
! 
! Author: D. Bundy
!-----------------------------------------------------------------------

  use iso_c_binding, only: c_int64_t, c_loc, c_ptr
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use ppgrid,        only: pcols, pver
  use constituents,  only: pcnst, sflxnam, cnst_cam_outfld
  use cam_history,   only: outfld

  implicit none

  ! Arguments
   type(physics_state), intent(in)  :: state          ! state variables
   type(physics_ptend), target, intent(out) :: ptend  ! package tendencies
   real(r8),            intent(in)  :: deltat         ! timestep
   real(r8),            intent(in)  :: landfrac(pcols) ! Land fraction
   real(r8), target, intent(inout) :: cflx(pcols,pcnst) ! Surface constituent flux (kg/m^2/s)

! Local variables
   integer  :: m               ! tracer number (internal)

   logical  :: lq(pcnst)
   interface
      subroutine tracers_timestep_tend_codon(ncol_c, pcols_c, pver_c, psetcols_c, &
           ixtrct_c, trac_ncnst_c, ptend_q_p, cflx_p) bind(c, name="tracers_timestep_tend_codon")
        use iso_c_binding, only: c_int64_t, c_ptr
        integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, psetcols_c
        integer(c_int64_t), value :: ixtrct_c, trac_ncnst_c
        type(c_ptr), value :: ptend_q_p, cflx_p
      end subroutine tracers_timestep_tend_codon
   end interface
!-----------------------------------------------------------------------

  call tracers_select_impl()

  if (.not. tracers_flag) then
       if (.not. use_native_impl) then
          if (tracers_flag_codon(0_c_int64_t) /= 0_c_int64_t) then
             call tracers_timestep_tend_native(state, ptend, cflx, landfrac, deltat)
             return
          end if
          call tracers_log_direct(tracers_timestep_tend_logged, &
               'tracers_timestep_tend direct = codon flag-off control shell; native empty ptend allocation boundary')
       end if
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
  endif

  if (use_native_impl) then
     call tracers_timestep_tend_native(state, ptend, cflx, landfrac, deltat)
     return
  end if

  lq(:)      = .FALSE.
  lq(ixtrct:ixtrct+trac_ncnst-1) = .TRUE.
  call physics_ptend_init(ptend, state%psetcols, 'tracers', lq=lq)

  call tracers_timestep_tend_codon( &
       int(state%ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), &
       int(state%psetcols, c_int64_t), int(ixtrct, c_int64_t), int(trac_ncnst, c_int64_t), &
       c_loc(ptend%q), c_loc(cflx) &
  )
  call tracers_log_direct(tracers_timestep_tend_logged, &
       'tracers_timestep_tend direct = codon tendency/flux body; native ptend/history/debug boundaries')

  do  m = 1,trac_ncnst
     if (debug) write(iulog,*)'tracers.F90 calling for tracer ',m

     if ( cnst_cam_outfld(ixtrct+m-1) ) then
        call outfld (sflxnam(ixtrct+m-1),cflx(:,ixtrct+m-1),pcols,state%lchnk)
     end if
  end do

  if ( debug ) then
     do  m = 1,trac_ncnst
        write(iulog,*)'tracers_timestep_tend ixtrct,m,ixtrct+m-1',ixtrct,m,ixtrct+m-1
        write(iulog,*)'tracers_timestep_tend min max flux',minval(cflx(:,ixtrct+m-1)),maxval(cflx(:,ixtrct+m-1))
        write(iulog,*)'tracers_timestep_tend min max tend',minval(ptend%q(:,:,ixtrct+m-1)),maxval(ptend%q(:,:,ixtrct+m-1))
     end do
     write(iulog,*)'tracers_timestep_tend end'
  endif

end subroutine tracers_timestep_tend

!======================================================================

subroutine tracers_timestep_tend_native(state, ptend, cflx, landfrac, deltat)

  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use ppgrid,        only: pcols, pver
  use constituents,  only: pcnst, sflxnam, cnst_cam_outfld
  use tracers_suite, only: flux_tr, tend_tr
  use cam_history,   only: outfld

  implicit none

  ! Arguments
   type(physics_state), intent(in)  :: state          ! state variables
   type(physics_ptend), intent(out) :: ptend          ! package tendencies
   real(r8),            intent(in)  :: deltat         ! timestep
   real(r8),            intent(in)  :: landfrac(pcols) ! Land fraction
   real(r8),            intent(inout) :: cflx(pcols,pcnst) ! Surface constituent flux (kg/m^2/s)

 ! Local variables
   integer  :: m               ! tracer number (internal)

   logical  :: lq(pcnst)

 !-----------------------------------------------------------------------

   if (.not. tracers_flag) then
      call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
      return
   else
      lq(:)      = .FALSE.
      lq(ixtrct:ixtrct+trac_ncnst-1) = .TRUE.
      call physics_ptend_init(ptend, state%psetcols, 'tracers', lq=lq)

      do  m = 1,trac_ncnst
         if (debug) write(iulog,*)'tracers.F90 calling for tracer ',m

         !calculate flux
         call flux_tr(m,state%ncol,state%lchnk, landfrac, cflx(:,ixtrct+m-1))

         !calculate tendency
         call tend_tr(m,state%ncol, state%q(:,:,ixtrct+m-1), deltat, ptend%q(:,:,ixtrct+m-1))

         !outfld calls could go here
         if ( cnst_cam_outfld(ixtrct+m-1) ) then
            call outfld (sflxnam(ixtrct+m-1),cflx(:,ixtrct+m-1),pcols,state%lchnk)
         end if
      end do

      if ( debug ) then
         do  m = 1,trac_ncnst
            write(iulog,*)'tracers_timestep_tend ixtrct,m,ixtrct+m-1',ixtrct,m,ixtrct+m-1
            write(iulog,*)'tracers_timestep_tend min max flux',minval(cflx(:,ixtrct+m-1)),maxval(cflx(:,ixtrct+m-1))
            write(iulog,*)'tracers_timestep_tend min max tend',minval(ptend%q(:,:,ixtrct+m-1)),maxval(ptend%q(:,:,ixtrct+m-1))
         end do
         write(iulog,*)'tracers_timestep_tend end'
      endif
   endif

end subroutine tracers_timestep_tend_native

end module tracers
