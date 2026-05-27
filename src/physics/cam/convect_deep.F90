
module convect_deep
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to several deep convection interfaces. Currently includes:
!    Zhang-McFarlane (default)
!    Kerry Emanuel 
!
!
! Author: D.B. Coleman, Sep 2004
!
!---------------------------------------------------------------------------------
   use shr_kind_mod, only: r8=>shr_kind_r8
   use ppgrid,       only: pver, pcols, pverp, begchunk, endchunk
   use cam_logfile,  only: iulog
   use spmd_utils,   only: masterproc
   use iso_c_binding, only: c_int64_t

   implicit none

   save
   private                         ! Make default type private to the module

! Public methods

   public ::&
      convect_deep_register,           &! register fields in physics buffer
      convect_deep_init,               &! initialize donner_deep module
      convect_deep_tend,               &! return tendencies
      convect_deep_tend_2,             &! return tendencies
      deep_scheme_does_scav_trans             ! = .t. if scheme does scavenging and conv. transport

! Private module data
   character(len=16) :: deep_scheme    ! default set in phys_control.F90, use namelist to change
   logical :: use_native_impl = .false.
   logical :: impl_selected = .false.
   integer :: codon_scheme_code = 0
   logical :: codon_scheme_selected = .false.
   logical :: convect_deep_tend_logged = .false.
   logical :: convect_deep_tend_2_logged = .false.
   logical :: convect_deep_register_logged = .false.
   logical :: convect_deep_init_logged = .false.
! Physics buffer indices 
   integer     ::  icwmrdp_idx      = 0 
   integer     ::  rprddp_idx       = 0 
   integer     ::  nevapr_dpcu_idx  = 0 
   integer     ::  cldtop_idx       = 0 
   integer     ::  cldbot_idx       = 0 
   integer     ::  cld_idx          = 0 
   integer     ::  fracis_idx       = 0 

   integer     ::  pblh_idx        = 0 
   integer     ::  tpert_idx       = 0 
   integer     ::  prec_dp_idx     = 0
   integer     ::  snow_dp_idx     = 0

   integer     ::  ttend_dp_idx        = 0

   interface
      function deep_scheme_does_scav_trans_codon(scheme_code_c) result(flag_c) &
           bind(c, name="deep_scheme_does_scav_trans_codon")
        use iso_c_binding, only: c_int64_t
        integer(c_int64_t), value :: scheme_code_c
        integer(c_int64_t) :: flag_c
      end function deep_scheme_does_scav_trans_codon
      function convect_deep_register_codon(flag_c) result(out_c) bind(c, name="convect_deep_register_codon")
        use iso_c_binding, only: c_int64_t
        integer(c_int64_t), value :: flag_c
        integer(c_int64_t) :: out_c
      end function convect_deep_register_codon
      function convect_deep_init_codon(flag_c) result(out_c) bind(c, name="convect_deep_init_codon")
        use iso_c_binding, only: c_int64_t
        integer(c_int64_t), value :: flag_c
        integer(c_int64_t) :: out_c
      end function convect_deep_init_codon
      function convect_deep_tend_2_action_codon(scheme_code_c) result(action_c) &
           bind(c, name="convect_deep_tend_2_action_codon")
        use iso_c_binding, only: c_int64_t
        integer(c_int64_t), value :: scheme_code_c
        integer(c_int64_t) :: action_c
      end function convect_deep_tend_2_action_codon
      function convect_deep_tend_action_codon(scheme_code_c) result(action_c) &
           bind(c, name="convect_deep_tend_action_codon")
        use iso_c_binding, only: c_int64_t
        integer(c_int64_t), value :: scheme_code_c
        integer(c_int64_t) :: action_c
      end function convect_deep_tend_action_codon
   end interface

!=========================================================================================
  contains 

!=========================================================================================
function deep_scheme_does_scav_trans()
  use iso_c_binding, only: c_int64_t
!
! Function called by tphysbc to determine if it needs to do scavenging and convective transport
! or if those have been done by the deep convection scheme. Each scheme could have its own
! identical query function for a less-knowledgable interface but for now, we know that KE 
! does scavenging & transport, and ZM doesn't
!

  logical deep_scheme_does_scav_trans

  call convect_deep_select_impl()

  if (.not. use_native_impl) then
     call convect_deep_select_codon_scheme()
     deep_scheme_does_scav_trans = deep_scheme_does_scav_trans_codon(int(codon_scheme_code, c_int64_t)) /= 0_c_int64_t
     return
  end if

  deep_scheme_does_scav_trans = .false.

  if ( deep_scheme .eq. 'KE' ) deep_scheme_does_scav_trans = .true.

  return

end function deep_scheme_does_scav_trans

!=========================================================================================
subroutine convect_deep_register

!----------------------------------------
! Purpose: register fields with the physics buffer
!----------------------------------------

  
  use physics_buffer, only : pbuf_add_field, dtype_r8
  use zm_conv_intr, only: zm_conv_register
  use phys_control, only: phys_getopts, use_gw_convect_dp

  implicit none

  integer idx
  integer(c_int64_t) :: active_c

  ! get deep_scheme setting from phys_control
  call phys_getopts(deep_scheme_out = deep_scheme)
  active_c = convect_deep_register_codon(1_c_int64_t)
  if (active_c == 0_c_int64_t) return
  call convect_deep_log_direct(convect_deep_register_logged, &
       'convect_deep_register direct = codon; zm_conv_register/pbuf_add_field native CAM API islands')

  select case ( deep_scheme )
  case('ZM') !    Zhang-McFarlane (default)
     call zm_conv_register
  end select

  call pbuf_add_field('ICWMRDP',    'physpkg',dtype_r8,(/pcols,pver/),icwmrdp_idx)
  call pbuf_add_field('RPRDDP',     'physpkg',dtype_r8,(/pcols,pver/),rprddp_idx)
  call pbuf_add_field('NEVAPR_DPCU','physpkg',dtype_r8,(/pcols,pver/),nevapr_dpcu_idx)
  call pbuf_add_field('PREC_DP',    'physpkg',dtype_r8,(/pcols/),     prec_dp_idx)
  call pbuf_add_field('SNOW_DP',   'physpkg',dtype_r8,(/pcols/),      snow_dp_idx)

  ! If gravity waves from deep convection are on, output this field.
  if (use_gw_convect_dp) then
     call pbuf_add_field('TTEND_DP','physpkg',dtype_r8,(/pcols,pver/),ttend_dp_idx)
  end if

end subroutine convect_deep_register

!=========================================================================================



subroutine convect_deep_init(pref_edge)

!----------------------------------------
! Purpose:  declare output fields, initialize variables needed by convection
!----------------------------------------

  use cam_history,    only: phys_decomp, addfld                          
  use pmgrid,         only: plevp
  use spmd_utils,     only: masterproc
  use zm_conv_intr,   only: zm_conv_init
  use cam_abortutils, only: endrun
  
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index

  implicit none

  real(r8),intent(in) :: pref_edge(plevp)        ! reference pressures at interfaces
  integer(c_int64_t) :: active_c

  active_c = convect_deep_init_codon(1_c_int64_t)
  if (active_c == 0_c_int64_t) return
  call convect_deep_log_direct(convect_deep_init_logged, &
       'convect_deep_init direct = codon; scheme/action shell direct = codon; zm_conv_init/pbuf/history native CAM API islands')

  select case ( deep_scheme )
  case('off')
     if (masterproc) write(iulog,*)'convect_deep: no deep convection selected'
  case('CLUBB_SGS')
     if (masterproc) write(iulog,*)'convect_deep: CLUBB_SGS selected'
  case('ZM')
     if (masterproc) write(iulog,*)'convect_deep initializing Zhang-McFarlane convection'
     call zm_conv_init(pref_edge)
  case('UNICON')
     if (masterproc) write(iulog,*)'convect_deep: deep convection done by UNICON'
  case default
     if (masterproc) write(iulog,*)'WARNING: convect_deep: no deep convection scheme. May fail.'
  end select

  cldtop_idx = pbuf_get_index('CLDTOP')
  cldbot_idx = pbuf_get_index('CLDBOT')
  cld_idx    = pbuf_get_index('CLD')
  fracis_idx = pbuf_get_index('FRACIS')

  pblh_idx   = pbuf_get_index('pblh')
  tpert_idx  = pbuf_get_index('tpert')

  call addfld ('ICWMRDP  ', 'kg/kg   ', pver, 'A', 'Deep Convection in-cloud water mixing ratio '            ,phys_decomp)

end subroutine convect_deep_init
!=========================================================================================
!subroutine convect_deep_tend(state, ptend, tdt, pbuf)

subroutine convect_deep_tend( &
     mcon    ,cme     ,          &
     dlf     ,pflx    ,zdu      , &
     rliq    , &
     ztodt   , &
     state   ,ptend   ,landfrac ,&
     pbuf    ,wtdlf )


   use physics_types, only: physics_state, physics_ptend, physics_tend, physics_ptend_init
   
   use cam_history,    only: outfld
   use constituents,   only: pcnst, cnst_name
   use zm_conv_intr,   only: zm_conv_tend
   use cam_history,    only: outfld
   use physconst,      only: cpair
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index

 !Water tracers:
    use water_tracer_vars, only: trace_water, wtrc_ntype, wtrc_srfpcp_indices
    use water_types,    only: iwtvap, iwtcvrain, iwtcvsnow


! Arguments
   type(physics_state), intent(in ) :: state   ! Physics state variables
   type(physics_ptend), intent(out) :: ptend   ! individual parameterization tendencies
   

   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in) :: ztodt               ! 2 delta t (model time increment)
   real(r8), intent(in) :: landfrac(pcols)     ! Land fraction
      

   real(r8), intent(out) :: mcon(pcols,pverp)  ! Convective mass flux--m sub c
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)    ! cmf condensation - evaporation
   real(r8), intent(out) :: zdu(pcols,pver)    ! detraining mass flux

   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals

   !Water tracers:
   real(r8), intent(out) :: wtdlf(pcols,pver,wtrc_ntype(iwtvap))   !Detraining tracer liquid
   integer                         :: wtpcidx                      !Physics Buffer index
   integer                         :: wtsnidx                      !Physics Buffer index
   real(r8), pointer, dimension(:) :: wtprec                       !Tracer surface rain 
   real(r8), pointer, dimension(:) :: wtsnow                       !Tracer surface snow

   real(r8), pointer :: prec(:)   ! total precipitation
   real(r8), pointer :: snow(:)   ! snow from ZM convection 

   real(r8), pointer, dimension(:) :: jctop
   real(r8), pointer, dimension(:) :: jcbot
   real(r8), pointer, dimension(:,:,:) :: cld        
   real(r8), pointer, dimension(:,:) :: ql        ! wg grid slice of cloud liquid water.
   real(r8), pointer, dimension(:,:) :: rprd      ! rain production rate
   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

   real(r8), pointer, dimension(:,:) :: evapcdp   ! Evaporation of deep convective precipitation

   real(r8), pointer :: pblh(:)                ! Planetary boundary layer height
   real(r8), pointer :: tpert(:)               ! Thermal temperature excess

   ! Temperature tendency from deep convection (pbuf pointer).
   real(r8), pointer, dimension(:,:) :: ttend_dp

   real(r8) zero(pcols, pver)

   integer i, k
   integer :: scheme_code
   integer(c_int64_t) :: tend_action_c

   call convect_deep_select_impl()
   if (.not. use_native_impl) call convect_deep_select_codon_scheme()

   call pbuf_get_field(pbuf, cldtop_idx,  jctop )
   call pbuf_get_field(pbuf, cldbot_idx,  jcbot )
   call pbuf_get_field(pbuf, icwmrdp_idx, ql    )

   if (use_native_impl) then
      select case ( deep_scheme )
      case('ZM')
         scheme_code = 1
      case('off')
         scheme_code = 2
      case('UNICON')
         scheme_code = 3
      case('CLUBB_SGS')
         scheme_code = 4
      case default
         scheme_code = -1
      end select
   else
      scheme_code = codon_scheme_code
   end if

   if (.not. use_native_impl) then
      tend_action_c = convect_deep_tend_action_codon(int(codon_scheme_code, c_int64_t))
      call convect_deep_log_tend_direct()
      if (tend_action_c == 1_c_int64_t) then
         scheme_code = 1
      else if (tend_action_c == 2_c_int64_t) then
         scheme_code = 2
      else
         scheme_code = -1
      end if
   end if

  select case ( scheme_code )
  case(2, 3, 4) ! in UNICON case the run method is called from convect_shallow_tend
    zero = 0     
    mcon = 0
    dlf = 0
    pflx = 0
    cme = 0
    zdu = 0
    rliq = 0

    !water tracers: 
    wtdlf(:,:,:) = 0

    call physics_ptend_init(ptend, state%psetcols, 'convect_deep')

!
! Associate pointers with physics buffer fields
!

    call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1/),   kount=(/pcols,pver/) ) 
    call pbuf_get_field(pbuf, rprddp_idx,      rprd )
    call pbuf_get_field(pbuf, fracis_idx,      fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )
    call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )
    call pbuf_get_field(pbuf, prec_dp_idx,     prec )
    call pbuf_get_field(pbuf, snow_dp_idx,     snow )

    prec=0
    snow=0

    jctop = pver
    jcbot = 1._r8
    cld = 0
    ql = 0
    rprd = 0
    fracis = 0
    evapcdp = 0

    !Water tracers
    do i=1,wtrc_ntype(iwtcvrain)
      call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvrain,i), wtprec)
      call pbuf_get_field(pbuf, wtrc_srfpcp_indices(iwtcvsnow,i), wtsnow)
      wtprec(:) = 0._r8 !assign values
      wtsnow(:) = 0._r8
    end do

  case(1) !    1 ==> Zhang-McFarlane (default)
     call pbuf_get_field(pbuf, pblh_idx,  pblh)
     call pbuf_get_field(pbuf, tpert_idx, tpert)

     call zm_conv_tend( pblh    ,mcon    ,cme     , &
          tpert   ,dlf     ,pflx    ,zdu      , &
          rliq    , &
          ztodt   , &
          jctop, jcbot , &
          state   ,ptend   ,landfrac ,pbuf , &
          wtdlf )

  end select

  ! If we added temperature tendency to pbuf, set it now.

  if (ttend_dp_idx > 0) then
     call pbuf_get_field(pbuf, ttend_dp_idx, ttend_dp)
     ttend_dp(:state%ncol,:pver) = ptend%s(:state%ncol,:pver)/cpair
  end if

  call outfld( 'ICWMRDP ', ql  , pcols, state%lchnk )


end subroutine convect_deep_tend
!=========================================================================================


subroutine convect_deep_tend_2( state,  ptend,  ztodt, pbuf)

   use physics_types, only: physics_state, physics_ptend, physics_ptend_init
   
   use physics_buffer,  only: physics_buffer_desc
   use constituents, only: pcnst
   use zm_conv_intr, only: zm_conv_tend_2

! Arguments
   type(physics_state), intent(in ) :: state          ! Physics state variables
   type(physics_ptend), intent(out) :: ptend          ! indivdual parameterization tendencies
   
   type(physics_buffer_desc), pointer :: pbuf(:)

   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   integer(c_int64_t) :: action_c


   call convect_deep_select_impl()
   if (.not. use_native_impl) call convect_deep_select_codon_scheme()

   if (.not. use_native_impl) then
      action_c = convect_deep_tend_2_action_codon(int(codon_scheme_code, c_int64_t))
      call convect_deep_log_tend_2_direct()
      if (action_c == 1_c_int64_t) then
         call zm_conv_tend_2( state,   ptend,  ztodt,  pbuf)
      else
         call physics_ptend_init(ptend, state%psetcols, 'convect_deep')
      end if
      return
   end if

   if (deep_scheme .eq. 'ZM') then
      call zm_conv_tend_2( state,   ptend,  ztodt,  pbuf)
   else
      call physics_ptend_init(ptend, state%psetcols, 'convect_deep')
   end if

end subroutine convect_deep_tend_2

!=========================================================================================

subroutine convect_deep_append_proof(proof_line)

   character(len=*), intent(in) :: proof_line
   character(len=512) :: proof_file
   integer :: status, n, unitno

   proof_file = ''
   call get_environment_variable('CONVECT_DEEP_PROOF_FILE', value=proof_file, length=n, status=status)
   if (status == 0 .and. n > 0) then
      open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
      write(unitno,'(A)') trim(proof_line)
      close(unitno)
   end if

end subroutine convect_deep_append_proof

!=========================================================================================

subroutine convect_deep_log_tend_direct()

   if (convect_deep_tend_logged) return
   convect_deep_tend_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') &
           'convect_deep_tend direct = codon; scheme/action dispatch direct = codon; zm_conv_tend/pbuf/outfld native islands'
      call convect_deep_append_proof( &
           'convect_deep_tend direct = codon; scheme/action dispatch direct = codon; zm_conv_tend/pbuf/outfld native islands')
      call flush(iulog)
   end if

end subroutine convect_deep_log_tend_direct

!=========================================================================================

subroutine convect_deep_log_tend_2_direct()

   if (convect_deep_tend_2_logged) return
   convect_deep_tend_2_logged = .true.

   if (masterproc) then
      write(iulog,'(A)') 'convect_deep_tend_2 direct = codon; zm_conv_tend_2/physics_ptend_init native CAM API islands'
      call convect_deep_append_proof( &
           'convect_deep_tend_2 direct = codon; zm_conv_tend_2/physics_ptend_init native CAM API islands')
      call flush(iulog)
   end if

end subroutine convect_deep_log_tend_2_direct

!=========================================================================================

subroutine convect_deep_log_direct(logged, proof_line)

   logical, intent(inout) :: logged
   character(len=*), intent(in) :: proof_line

   if (logged) return
   logged = .true.

   if (masterproc) then
      write(iulog,'(A)') trim(proof_line)
      call convect_deep_append_proof(trim(proof_line))
      call flush(iulog)
   end if

end subroutine convect_deep_log_direct

!=========================================================================================

subroutine convect_deep_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (impl_selected) return

   impl_name = 'codon'
   call get_environment_variable('CONVECT_DEEP_IMPL', value=impl_name, length=n, status=status)

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
         write(iulog,*) 'convect_deep_tend implementation = native'
      else
         write(iulog,*) 'convect_deep_tend implementation = codon'
      end if
      call flush(iulog)
   end if

end subroutine convect_deep_select_impl

!=========================================================================================

subroutine convect_deep_select_codon_scheme()

   use iso_c_binding, only: c_int64_t, c_loc, c_ptr

   integer :: i
   integer(c_int64_t), target :: scheme_ascii(len(deep_scheme))
   integer(c_int64_t), target :: scheme_code
   integer(c_int64_t), target :: status_code

   interface
      subroutine convect_deep_select_scheme_codon(scheme_len_c, scheme_ascii_p, scheme_code_p, status_p) &
           bind(c, name="convect_deep_select_scheme_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: scheme_len_c
         type(c_ptr), value :: scheme_ascii_p, scheme_code_p, status_p
      end subroutine convect_deep_select_scheme_codon
   end interface

   if (codon_scheme_selected) return

   do i = 1, len(deep_scheme)
      scheme_ascii(i) = int(iachar(deep_scheme(i:i)), c_int64_t)
   end do

   scheme_code = 0_c_int64_t
   status_code = 0_c_int64_t
   call convect_deep_select_scheme_codon( &
        int(len(deep_scheme), c_int64_t), c_loc(scheme_ascii(1)), c_loc(scheme_code), c_loc(status_code) &
   )

   if (status_code /= 0_c_int64_t) then
      codon_scheme_code = -1
   else
      codon_scheme_code = int(scheme_code)
   end if

   codon_scheme_selected = .true.

end subroutine convect_deep_select_codon_scheme


end module convect_deep
