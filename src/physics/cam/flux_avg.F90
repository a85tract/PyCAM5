module flux_avg

!---------------------------------------------------------------------------------
! Purpose: Contains code to smooth the surface fluxes to reduce
!          instabilities in the surface layer.
!---------------------------------------------------------------------------------

  use shr_kind_mod,     only: r8=>shr_kind_r8
  use ppgrid,           only: begchunk, endchunk, pcols
  use cam_logfile,      only: iulog
  use spmd_utils,       only: masterproc
  
  use physics_types,    only: physics_state
  use camsrfexch,       only: cam_in_t    
  use phys_grid,        only: get_ncols_p
  use physics_buffer, only : pbuf_add_field, dtype_r8
  implicit none
  private
  save

  ! Public interfaces

  public :: flux_avg_register
  public :: flux_avg_init
  public :: flux_avg_run
  
  ! Private module data

  integer :: lhflx_idx      ! lhflx index in physics buffer
  integer :: shflx_idx      ! shflx index in physics buffer
  integer :: qflx_idx       ! qflx index in physics buffer
  integer :: taux_idx       ! taux index in physics buffer
  integer :: tauy_idx       ! tauy index in physics buffer
  integer :: lhflx_res_idx  ! lhflx_res index in physics buffer
  integer :: shflx_res_idx  ! shflx_res index in physics buffer
  integer :: qflx_res_idx   ! qflx_res index in physics buffer
  integer :: taux_res_idx   ! taux_res index in physics buffer
  integer :: tauy_res_idx   ! tauy_res index in physics buffer
  logical :: use_native_impl = .false.
  logical :: impl_selected = .false.
  logical :: use_native_init_impl = .false.
  logical :: init_impl_selected = .false.

!===============================================================================
contains
!===============================================================================

subroutine flux_avg_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('FLUX_AVG_IMPL', impl_name, n, status)

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
         write(iulog,*) 'flux_avg implementation = native'
      else
         write(iulog,*) 'flux_avg implementation = codon'
      end if
   end if

end subroutine flux_avg_select_impl

subroutine flux_avg_init_select_impl()

   character(len=32) :: impl_name
   integer :: status, n, i, code

   if (init_impl_selected) return

   impl_name = 'codon'
   call cam_codon_get_impl('FLUX_AVG_INIT_IMPL', impl_name, n, status)

   if (status == 0 .and. n > 0) then
      do i = 1, n
         code = iachar(impl_name(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
         end if
      end do
      use_native_init_impl = trim(adjustl(impl_name(:n))) == 'native'
   else
      use_native_init_impl = .false.
   end if

   init_impl_selected = .true.

   if (masterproc) then
      if (use_native_init_impl) then
         write(iulog,*) 'flux_avg_init implementation = native'
      else
         write(iulog,*) 'flux_avg_init implementation = codon'
      end if
   end if

end subroutine flux_avg_init_select_impl

subroutine flux_avg_register()

   !----------------------------------------------------------------------
   !
   ! Register the fluxes in the physics buffer.
   ! 
   !-----------------------------------------------------------------------

   ! Request physics buffer space for fields that persist across timesteps.
   call pbuf_add_field('LHFLX',    'global',dtype_r8,(/pcols,1/),lhflx_idx)
   call pbuf_add_field('SHFLX',    'global',dtype_r8,(/pcols,1/),shflx_idx)
   call pbuf_add_field('TAUX',     'global',dtype_r8,(/pcols,1/),taux_idx)
   call pbuf_add_field('TAUY',     'global',dtype_r8,(/pcols,1/),tauy_idx)
   call pbuf_add_field('QFLX',     'global',dtype_r8,(/pcols,1/),qflx_idx)
   call pbuf_add_field('LHFLX_RES','global',dtype_r8,(/pcols,1/),lhflx_res_idx)
   call pbuf_add_field('SHFLX_RES','global',dtype_r8,(/pcols,1/),shflx_res_idx)
   call pbuf_add_field('TAUX_RES', 'global',dtype_r8,(/pcols,1/),taux_res_idx)
   call pbuf_add_field('TAUY_RES', 'global',dtype_r8,(/pcols,1/),tauy_res_idx)
   call pbuf_add_field('QFLX_RES', 'global',dtype_r8,(/pcols,1/),qflx_res_idx)

end subroutine flux_avg_register

!===============================================================================

subroutine flux_avg_init(cam_in,  pbuf2d)
  use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
  use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk, pbuf_get_field
   ! Initialize the surface fluxes in the physics buffer using the cam import state

   type(cam_in_t), target, intent(in)    :: cam_in(begchunk:endchunk)
   
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   integer :: lchnk
   integer :: ncol
   type(physics_buffer_desc), pointer :: pbuf2d_chunk(:)
   real(r8), pointer, dimension(:) :: lhflx, shflx, qflx, taux, tauy
   real(r8), pointer, dimension(:) :: lhflx_res, shflx_res, qflx_res, taux_res, tauy_res
   interface
      subroutine flux_avg_init_codon(ncol_c, pcols_c, cam_lhf_p, cam_shf_p, cam_cflx1_p, cam_wsx_p, cam_wsy_p, &
           lhflx_p, shflx_p, qflx_p, taux_p, tauy_p, &
           lhflx_res_p, shflx_res_p, qflx_res_p, taux_res_p, tauy_res_p) bind(c, name="flux_avg_init_codon")
        use iso_c_binding, only : c_int64_t, c_ptr
        integer(c_int64_t), value :: ncol_c, pcols_c
        type(c_ptr), value :: cam_lhf_p, cam_shf_p, cam_cflx1_p, cam_wsx_p, cam_wsy_p
        type(c_ptr), value :: lhflx_p, shflx_p, qflx_p, taux_p, tauy_p
        type(c_ptr), value :: lhflx_res_p, shflx_res_p, qflx_res_p, taux_res_p, tauy_res_p
      end subroutine flux_avg_init_codon
   end interface

   !----------------------------------------------------------------------- 

   call flux_avg_init_select_impl()

   if (use_native_init_impl) then
      call flux_avg_init_native(cam_in, pbuf2d)
      return
   end if

   do lchnk = begchunk, endchunk
      ncol = get_ncols_p(lchnk)
      pbuf2d_chunk => pbuf_get_chunk(pbuf2d, lchnk)

      call pbuf_get_field(pbuf2d_chunk, lhflx_idx,     lhflx )
      call pbuf_get_field(pbuf2d_chunk, shflx_idx,     shflx )
      call pbuf_get_field(pbuf2d_chunk, qflx_idx,      qflx  )
      call pbuf_get_field(pbuf2d_chunk, taux_idx,      taux  )
      call pbuf_get_field(pbuf2d_chunk, tauy_idx,      tauy  )
      call pbuf_get_field(pbuf2d_chunk, lhflx_res_idx, lhflx_res )
      call pbuf_get_field(pbuf2d_chunk, shflx_res_idx, shflx_res )
      call pbuf_get_field(pbuf2d_chunk, qflx_res_idx,  qflx_res  )
      call pbuf_get_field(pbuf2d_chunk, taux_res_idx,  taux_res  )
      call pbuf_get_field(pbuf2d_chunk, tauy_res_idx,  tauy_res  )

      call flux_avg_init_codon( &
           int(ncol, c_int64_t), int(pcols, c_int64_t), &
           c_loc(cam_in(lchnk)%lhf), c_loc(cam_in(lchnk)%shf), c_loc(cam_in(lchnk)%cflx(1,1)), &
           c_loc(cam_in(lchnk)%wsx), c_loc(cam_in(lchnk)%wsy), &
           c_loc(lhflx), c_loc(shflx), c_loc(qflx), c_loc(taux), c_loc(tauy), &
           c_loc(lhflx_res), c_loc(shflx_res), c_loc(qflx_res), c_loc(taux_res), c_loc(tauy_res) &
      )
   end do


end subroutine flux_avg_init

!===============================================================================

subroutine flux_avg_init_native(cam_in,  pbuf2d)
  use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_get_chunk
   ! Initialize the surface fluxes in the physics buffer using the cam import state

   type(cam_in_t),      intent(in)    :: cam_in(begchunk:endchunk)
   
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   integer :: lchnk
   integer :: ncol
   type(physics_buffer_desc), pointer :: pbuf2d_chunk(:)

   !----------------------------------------------------------------------- 

   do lchnk = begchunk, endchunk
      ncol = get_ncols_p(lchnk)
      pbuf2d_chunk => pbuf_get_chunk(pbuf2d, lchnk)
      call pbuf_set_field(pbuf2d_chunk, lhflx_idx,  cam_in(lchnk)%lhf(:ncol))
      call pbuf_set_field(pbuf2d_chunk, shflx_idx,  cam_in(lchnk)%shf(:ncol))
      call pbuf_set_field(pbuf2d_chunk, qflx_idx,   cam_in(lchnk)%cflx(:ncol,1))
      call pbuf_set_field(pbuf2d_chunk, taux_idx,   cam_in(lchnk)%wsx(:ncol))
      call pbuf_set_field(pbuf2d_chunk, tauy_idx,   cam_in(lchnk)%wsy(:ncol))

      call pbuf_set_field(pbuf2d,       shflx_res_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d_chunk, lhflx_res_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d_chunk, qflx_res_idx,  0.0_r8)
      call pbuf_set_field(pbuf2d_chunk, taux_res_idx,  0.0_r8)
      call pbuf_set_field(pbuf2d_chunk, tauy_res_idx,  0.0_r8)
   end do


end subroutine flux_avg_init_native

!===============================================================================

subroutine flux_avg_run(state, cam_in,  pbuf, nstep, deltat)
  use iso_c_binding, only : c_double, c_int64_t, c_loc, c_ptr
  use physics_buffer, only : physics_buffer_desc, pbuf_get_field
   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   !
   !----------------------------------------------------------------------- 
!++ debug code to be removed after PBL code validated
   use phys_debug,       only: phys_debug_flux1, phys_debug_flux2
!-- debug code to be removed after PBL code validated

   ! Input arguments

   type(physics_state), intent(in)    :: state
   type(cam_in_t), target, intent(inout) :: cam_in
   type(physics_buffer_desc), pointer :: pbuf(:)
   
   integer,             intent(in)    :: nstep
   real(r8),            intent(in)    :: deltat

   ! Local variables
   integer :: lchnk                  ! chunk identifier
   integer :: ncol                   ! number of atmospheric columns

   ! physics buffer fields
   real(r8), pointer, dimension(:) :: lhflx   ! latent heat flux
   real(r8), pointer, dimension(:) :: shflx   ! sensible heat flux
   real(r8), pointer, dimension(:) :: qflx    ! water vapor heat flux
   real(r8), pointer, dimension(:) :: taux    ! x momentum flux
   real(r8), pointer, dimension(:) :: tauy    ! y momentum flux
   real(r8), pointer, dimension(:) :: lhflx_res   ! latent heat flux
   real(r8), pointer, dimension(:) :: shflx_res   ! sensible heat flux
   real(r8), pointer, dimension(:) :: qflx_res    ! water vapor heat flux
   real(r8), pointer, dimension(:) :: taux_res    ! x momentum flux
   real(r8), pointer, dimension(:) :: tauy_res    ! y momentum flux
   real(r8), target :: temp(pcols)
   interface
      subroutine flux_avg_run_codon(ncol_c, pcols_c, nstep_c, deltat_c, &
           cam_lhf_p, cam_shf_p, cam_wsx_p, cam_wsy_p, cam_cflx1_p, &
           lhflx_p, shflx_p, qflx_p, taux_p, tauy_p, &
           lhflx_res_p, shflx_res_p, qflx_res_p, taux_res_p, tauy_res_p, temp_p) &
           bind(c, name="flux_avg_run_codon")
        use iso_c_binding, only : c_double, c_int64_t, c_ptr
        integer(c_int64_t), value :: ncol_c, pcols_c, nstep_c
        real(c_double), value :: deltat_c
        type(c_ptr), value :: cam_lhf_p, cam_shf_p, cam_wsx_p, cam_wsy_p, cam_cflx1_p
        type(c_ptr), value :: lhflx_p, shflx_p, qflx_p, taux_p, tauy_p
        type(c_ptr), value :: lhflx_res_p, shflx_res_p, qflx_res_p, taux_res_p, tauy_res_p, temp_p
      end subroutine flux_avg_run_codon
   end interface
   !----------------------------------------------------------------------- 

   call flux_avg_select_impl()

   if (use_native_impl) then
      call flux_avg_run_native(state, cam_in, pbuf, nstep, deltat)
      return
   end if

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Associate pointers with physics buffer fields
   call pbuf_get_field(pbuf, lhflx_idx,     lhflx )
   call pbuf_get_field(pbuf, shflx_idx,     shflx )
   call pbuf_get_field(pbuf, qflx_idx,      qflx  )
   call pbuf_get_field(pbuf, taux_idx,      taux  )
   call pbuf_get_field(pbuf, tauy_idx,      tauy  )

   call pbuf_get_field(pbuf, lhflx_res_idx, lhflx_res )
   call pbuf_get_field(pbuf, shflx_res_idx, shflx_res )
   call pbuf_get_field(pbuf, qflx_res_idx,  qflx_res  )
   call pbuf_get_field(pbuf, taux_res_idx,  taux_res  )
   call pbuf_get_field(pbuf, tauy_res_idx,  tauy_res  )

!++ debug code to be removed after PBL code validated
   call phys_debug_flux1(lchnk, cam_in, lhflx, shflx, taux, tauy, qflx, &
                         lhflx_res, shflx_res, taux_res, tauy_res, qflx_res)
!-- debug code to be removed after PBL code validated

   call flux_avg_run_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), int(nstep, c_int64_t), real(deltat, c_double), &
        c_loc(cam_in%lhf), c_loc(cam_in%shf), c_loc(cam_in%wsx), c_loc(cam_in%wsy), c_loc(cam_in%cflx(1,1)), &
        c_loc(lhflx), c_loc(shflx), c_loc(qflx), c_loc(taux), c_loc(tauy), &
        c_loc(lhflx_res), c_loc(shflx_res), c_loc(qflx_res), c_loc(taux_res), c_loc(tauy_res), c_loc(temp) &
   )

!++ debug code to be removed after PBL code validated
   call phys_debug_flux2(lchnk, cam_in, lhflx, &
                         lhflx_res, shflx_res, taux_res, tauy_res, qflx_res)
!-- debug code to be removed after PBL code validated

end subroutine flux_avg_run

!===============================================================================

subroutine flux_avg_run_native(state, cam_in,  pbuf, nstep, deltat)
  use physics_buffer, only : physics_buffer_desc, pbuf_get_field
   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   !
   !----------------------------------------------------------------------- 
!++ debug code to be removed after PBL code validated
   use phys_debug,       only: phys_debug_flux1, phys_debug_flux2
!-- debug code to be removed after PBL code validated

   ! Input arguments

   type(physics_state), intent(in)    :: state
   type(cam_in_t),      intent(inout) :: cam_in
   type(physics_buffer_desc), pointer :: pbuf(:)
   
   integer,             intent(in)    :: nstep
   real(r8),            intent(in)    :: deltat

   ! Local variables
   integer :: lchnk                  ! chunk identifier
   integer :: ncol                   ! number of atmospheric columns

   ! physics buffer fields
   real(r8), pointer, dimension(:) :: lhflx   ! latent heat flux
   real(r8), pointer, dimension(:) :: shflx   ! sensible heat flux
   real(r8), pointer, dimension(:) :: qflx    ! water vapor heat flux
   real(r8), pointer, dimension(:) :: taux    ! x momentum flux
   real(r8), pointer, dimension(:) :: tauy    ! y momentum flux
   real(r8), pointer, dimension(:) :: lhflx_res   ! latent heat flux
   real(r8), pointer, dimension(:) :: shflx_res   ! sensible heat flux
   real(r8), pointer, dimension(:) :: qflx_res    ! water vapor heat flux
   real(r8), pointer, dimension(:) :: taux_res    ! x momentum flux
   real(r8), pointer, dimension(:) :: tauy_res    ! y momentum flux
   !----------------------------------------------------------------------- 

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Associate pointers with physics buffer fields
   call pbuf_get_field(pbuf, lhflx_idx,     lhflx )
   call pbuf_get_field(pbuf, shflx_idx,     shflx )
   call pbuf_get_field(pbuf, qflx_idx,      qflx  )
   call pbuf_get_field(pbuf, taux_idx,      taux  )
   call pbuf_get_field(pbuf, tauy_idx,      tauy  )

   call pbuf_get_field(pbuf, lhflx_res_idx, lhflx_res )
   call pbuf_get_field(pbuf, shflx_res_idx, shflx_res )
   call pbuf_get_field(pbuf, qflx_res_idx,  qflx_res  )
   call pbuf_get_field(pbuf, taux_res_idx,  taux_res  )
   call pbuf_get_field(pbuf, tauy_res_idx,  tauy_res  )

!++ debug code to be removed after PBL code validated
   call phys_debug_flux1(lchnk, cam_in, lhflx, shflx, taux, tauy, qflx, &
                         lhflx_res, shflx_res, taux_res, tauy_res, qflx_res)
!-- debug code to be removed after PBL code validated

   call smooth (cam_in%lhf, lhflx, lhflx_res, nstep, deltat, ncol)
   call smooth (cam_in%shf, shflx, shflx_res, nstep, deltat, ncol)
   call smooth (cam_in%wsx, taux, taux_res, nstep, deltat, ncol)
   call smooth (cam_in%wsy, tauy, tauy_res, nstep, deltat, ncol)
   call smooth (cam_in%cflx(:pcols,1), qflx, qflx_res, nstep, deltat, ncol)

!++ debug code to be removed after PBL code validated
   call phys_debug_flux2(lchnk, cam_in, lhflx, &
                         lhflx_res, shflx_res, taux_res, tauy_res, qflx_res)
!-- debug code to be removed after PBL code validated

end subroutine flux_avg_run_native

!===============================================================================

subroutine smooth(new, old, res, nstep, deltat, ncol)

   real(r8), intent(inout) :: new(pcols)
   real(r8), intent(inout) :: old(pcols)
   real(r8), intent(inout) :: res(pcols)
   real(r8), intent(in)    :: deltat
   integer,  intent(in)    :: nstep
   integer,  intent(in)    :: ncol

   real(r8) :: temp(pcols)
   integer i

   temp(1:ncol) = new(1:ncol)
   if (nstep > 0) then
      new(1:ncol) = 0.5_r8*(new(1:ncol)+old(1:ncol))
   else
      old(1:ncol) = new(1:ncol)
      res(1:ncol) = 0._r8
   endif

   ! storing the old value for smoothing on the next step
   ! doesnt seem to be stable
   ! old(1:ncol) = temp(1:ncol)

   ! storing the smoothed value for the next step

   ! first add the flux that the surface model wanted to provide less
   ! the flux the atmosphere will actually see to the residual
   res(1:ncol) = res(1:ncol) + temp(1:ncol)-new(1:ncol)

   ! now calculate the amount that we might increment the new flux
   ! to include some of the residual
   ! If the residual is small we will just add it all, 
   ! but if it is large we will add it at the rate required to put
   ! the residual back into the flux over a 2 hour period
   do i = 1,ncol
      if (abs(res(i)).lt.max(abs(new(i)),abs(old(i)))*0.05_r8) then
         temp(i) = res(i)
         res(i) = 0._r8
      else
         temp(i) = res(i)*deltat/7200._r8
         !     temp(i) = res(i)*deltat*0.5/7200.
         res(i) = res(i)-temp(i)
      endif
   end do

   ! dont do conservative smoothing for first 12 hours
   if (nstep*deltat/86400._r8 < 0.5_r8) then
      ! use this line if your dont want to use the residual
      !if (.true.) then
      temp = 0._r8
      res = 0._r8
   endif

   ! make the new flux the average of the sfc model and last timestep
   ! plus some of the residual
   new(1:ncol) = new(1:ncol) + temp(1:ncol)
   old(1:ncol) = new(1:ncol)

end subroutine smooth

!===============================================================================

end module flux_avg
