
module rayleigh_friction

  !---------------------------------------------------------------------------------
  ! Module to apply rayleigh friction in region of model top.
  ! We specify a decay rate profile that is largest at the model top and
  ! drops off vertically using a hyperbolic tangent profile.
  ! We compute the tendencies in u and v using an Euler backward scheme.
  ! We then apply the negative of the kinetic energy tendency to "s", the dry
  ! static energy.
  !
  ! calling sequence:
  !
  !  rayleigh_friction_init          initializes rayleigh friction constants
  !  rayleigh_friction_tend          computes rayleigh friction tendencies
  !
  !---------------------------Code history--------------------------------
  ! This is a new routine written by Art Mirin in collaboration with Phil Rasch.
  ! Initial coding for this version:  Art Mirin, May 2007.
  !---------------------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use ppgrid,           only: pver
  use spmd_utils,       only: masterproc
  use cam_logfile,      only: iulog

  implicit none
  private          ! Make default type private to the module
  save
  
  !
  ! Public interfaces
  !
  public rayleigh_friction_init          ! Initialization
  public rayleigh_friction_tend          ! Computation of tendencies

  !
  ! Public data
  !
  integer, public   :: rayk0 = 2           ! vertical level at which rayleigh friction term is centered
  real (r8), public :: raykrange = 0._r8   ! range of rayleigh friction profile 
                                           ! if 0, range is set to satisfy x=2 (see below)
  real (r8), public :: raytau0 = 0._r8     ! approximate value of decay time at model top (days)
                                           ! if 0., no rayleigh friction is applied

  ! 
  ! Private data
  !
  logical :: use_native_impl = .false.
  logical :: impl_selected = .false.
  logical :: rayleigh_friction_init_logged = .false.
  logical :: rayleigh_friction_tend_logged = .false.
  real (r8), target :: krange         ! range of rayleigh friction profile
  real (r8), target :: tau0           ! approximate value of decay time at model top
  real (r8), target :: otau0          ! inverse of tau0
  real (r8), target :: otau(pver)     ! inverse decay time versus vertical level

  interface
    subroutine rayleigh_friction_init_codon(rayk0_c, raykrange_c, raytau0_c, pver_c, &
         krange_p, tau0_p, otau0_p, otau_p) bind(c, name="rayleigh_friction_init_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      integer(c_int64_t), value :: rayk0_c
      real(c_double), value :: raykrange_c, raytau0_c
      integer(c_int64_t), value :: pver_c
      type(c_ptr), value :: krange_p, tau0_p, otau0_p, otau_p
    end subroutine rayleigh_friction_init_codon

    subroutine rayleigh_friction_tend_codon(ztodt_c, pver_c, psetcols_c, ncol_c, &
         otau_p, state_u_p, state_v_p, ptend_u_p, ptend_v_p, ptend_s_p) &
         bind(c, name="rayleigh_friction_tend_codon")
      use iso_c_binding, only: c_double, c_int64_t, c_ptr
      real(c_double), value :: ztodt_c
      integer(c_int64_t), value :: pver_c, psetcols_c, ncol_c
      type(c_ptr), value :: otau_p, state_u_p, state_v_p, ptend_u_p, ptend_v_p, ptend_s_p
    end subroutine rayleigh_friction_tend_codon
  end interface

  ! We apply a profile of the form otau0 * [1 + tanh (x)] / 2 , where
  ! x = (k0 - k) / krange. The default is for x to equal 2 at k=1, meaning
  ! krange = (k0 - 1) / 2. The default is applied when raykrange is set to 0.
  ! If otau0 = 0, no term is applied.

contains

  !===============================================================================
  subroutine rayleigh_friction_log_direct(logged, proof_line)
    logical, intent(inout) :: logged
    character(len=*), intent(in) :: proof_line

    if (logged) return
    logged = .true.

    if (masterproc) then
       write(iulog,'(A)') trim(proof_line)
       call flush(iulog)
    end if
  end subroutine rayleigh_friction_log_direct

  !===============================================================================
  subroutine rayleigh_friction_select_impl()
    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('RAYLEIGH_FRICTION_IMPL', impl_name, n, status)

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
       impl_name = 'codon'
    end if

    impl_selected = .true.

    if (masterproc) then
       if (use_native_impl) then
          write(iulog,*) 'Rayleigh friction implementation = native'
       else
          write(iulog,*) 'Rayleigh friction implementation = codon'
       end if
    end if
  end subroutine rayleigh_friction_select_impl

  !===============================================================================
  subroutine rayleigh_friction_init()
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr

    !---------------------------Local storage-------------------------------
    integer k

    !-----------------------------------------------------------------------
    ! Compute tau array
    !-----------------------------------------------------------------------
    call rayleigh_friction_select_impl()

    if (use_native_impl) then
       call rayleigh_friction_init_native()
       return
    end if

    call rayleigh_friction_init_codon( &
         int(rayk0, c_int64_t), real(raykrange, c_double), real(raytau0, c_double), int(pver, c_int64_t), &
         c_loc(krange), c_loc(tau0), c_loc(otau0), c_loc(otau) &
    )
    call rayleigh_friction_log_direct(rayleigh_friction_init_logged, 'rayleigh_friction_init direct = codon')

    if (masterproc) then
       write (iulog,*) 'Rayleigh friction - rayk0 = ', rayk0
       write (iulog,*) 'Rayleigh friction - raykrange = ', raykrange
       write (iulog,*) 'Rayleigh friction - raytau0 = ', raytau0
       write (iulog,*) 'Rayleigh friction - krange = ', krange
       write (iulog,*) 'Rayleigh friction - otau0 = ', otau0
       write (iulog,*) 'Rayleigh friction decay rate profile'
       do k = 1, pver
          write (iulog,*) '   k = ', k, '   otau = ', otau(k)
       enddo
    end if

    return

  end subroutine rayleigh_friction_init

  !===============================================================================
  subroutine rayleigh_friction_init_native()
    !------------------------------Arguments--------------------------------

    !---------------------------Local storage-------------------------------
    real (r8) x
    integer k

    !-----------------------------------------------------------------------
    ! Compute tau array
    !-----------------------------------------------------------------------

    krange = raykrange
    if (raykrange .eq. 0._r8) krange = (rayk0 - 1) / 2._r8

    tau0 = (86400._r8) * raytau0   ! convert to seconds
    otau0 = 0._r8
    if (tau0 .ne. 0._r8) otau0 = 1._r8/tau0

    do k = 1, pver
       x = (rayk0 - k) / krange
       otau(k) = otau0 * (1 + tanh(x)) / (2._r8)
    enddo

    if (masterproc) then
       write (iulog,*) 'Rayleigh friction - rayk0 = ', rayk0
       write (iulog,*) 'Rayleigh friction - raykrange = ', raykrange
       write (iulog,*) 'Rayleigh friction - raytau0 = ', raytau0
       write (iulog,*) 'Rayleigh friction - krange = ', krange
       write (iulog,*) 'Rayleigh friction - otau0 = ', otau0
       write (iulog,*) 'Rayleigh friction decay rate profile'
       do k = 1, pver
          write (iulog,*) '   k = ', k, '   otau = ', otau(k)
       enddo
    end if

    return

  end subroutine rayleigh_friction_init_native
  
!=========================================================================================
  subroutine rayleigh_friction_tend(                                     &
       ztodt    ,state    ,ptend    )
    !-----------------------------------------------------------------------
    ! interface routine for rayleigh friction
    !-----------------------------------------------------------------------
    use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
    use physics_types, only: physics_state, physics_ptend, physics_ptend_init


    !------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ztodt                  ! physics timestep
    type(physics_state), target, intent(in)  :: state      ! physics state variables
    
    type(physics_ptend), target, intent(out) :: ptend      ! individual parameterization tendencies
    !
    !---------------------------Local storage-------------------------------
    !-----------------------------------------------------------------------

    call rayleigh_friction_select_impl()

    if (use_native_impl) then
       call rayleigh_friction_tend_native(ztodt, state, ptend)
       return
    end if

    call physics_ptend_init(ptend, state%psetcols, 'rayleigh friction', ls=.true., lu=.true., lv=.true.)

    if (otau0 .eq. 0._r8) then
       call rayleigh_friction_tend_codon( &
            real(ztodt, c_double), int(pver, c_int64_t), int(state%psetcols, c_int64_t), int(0, c_int64_t), &
            c_loc(otau), c_loc(state%u), c_loc(state%v), c_loc(ptend%u), c_loc(ptend%v), c_loc(ptend%s) &
       )
       call rayleigh_friction_log_direct(rayleigh_friction_tend_logged, 'rayleigh_friction_tend direct = codon')
       return
    end if

    call rayleigh_friction_tend_codon( &
         real(ztodt, c_double), int(pver, c_int64_t), int(state%psetcols, c_int64_t), int(state%ncol, c_int64_t), &
         c_loc(otau), c_loc(state%u), c_loc(state%v), c_loc(ptend%u), c_loc(ptend%v), c_loc(ptend%s) &
    )
    call rayleigh_friction_log_direct(rayleigh_friction_tend_logged, 'rayleigh_friction_tend direct = codon')

    return
  end subroutine rayleigh_friction_tend

!=========================================================================================
  subroutine rayleigh_friction_tend_native(                                     &
       ztodt    ,state    ,ptend    )
    !-----------------------------------------------------------------------
    ! interface routine for rayleigh friction
    !-----------------------------------------------------------------------
    use physics_types, only: physics_state, physics_ptend, physics_ptend_init


    !------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ztodt                  ! physics timestep
    type(physics_state), intent(in)  :: state      ! physics state variables
    
    type(physics_ptend), intent(out) :: ptend      ! individual parameterization tendencies
    !
    !---------------------------Local storage-------------------------------
    integer :: ncol                                ! number of atmospheric columns
    integer :: k                                   ! level
    real(r8) :: rztodt                             ! 1./ztodt
    real(r8) :: c1, c2, c3                         ! temporary variables
    !-----------------------------------------------------------------------

    call physics_ptend_init(ptend, state%psetcols, 'rayleigh friction', ls=.true., lu=.true., lv=.true.)

    if (otau0 .eq. 0._r8) return

    rztodt = 1._r8/ztodt
    ncol  = state%ncol

    ! u, v and s are modified by rayleigh friction

    do k = 1, pver
       c2 = 1._r8 / (1._r8 + otau(k)*ztodt)
       c1 = -otau(k) * c2
       c3 = 0.5_r8 * (1._r8 - c2*c2) * rztodt
       ptend%u(:ncol,k) = c1 * state%u(:ncol,k)
       ptend%v(:ncol,k) = c1 * state%v(:ncol,k)
       ptend%s(:ncol,k) = c3 * (state%u(:ncol,k)**2 + state%v(:ncol,k)**2)
    enddo

    return
  end subroutine rayleigh_friction_tend_native

end module rayleigh_friction
