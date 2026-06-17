module trb_mtn_stress

  implicit none
  private      
  save

  public init_tms                             ! Initialization
  public compute_tms                          ! Full routine

  ! ------------ !
  ! Private data !
  ! ------------ !

  integer,  parameter :: r8 = selected_real_kind(12) ! 8 byte real

  real(r8), parameter :: horomin= 1._r8       ! Minimum value of subgrid orographic height for mountain stress [ m ]
  real(r8), parameter :: z0max  = 100._r8     ! Maximum value of z_0 for orography [ m ]
  real(r8), parameter :: dv2min = 0.01_r8     ! Minimum shear squared [ m2/s2 ]
  real(r8), target    :: orocnst              ! Converts from standard deviation to height [ no unit ]
  real(r8), target    :: z0fac                ! Factor determining z_0 from orographic standard deviation [ no unit ]
  real(r8), target    :: karman               ! von Karman constant
  real(r8), target    :: gravit               ! Acceleration due to gravity
  real(r8), target    :: rair                 ! Gas constant for dry air

  logical             :: use_native_tms_impl = .false.
  logical             :: tms_impl_selected = .false.
  logical             :: tms_proof_written = .false.
  logical             :: use_native_init_tms_impl = .false.
  logical             :: init_tms_impl_selected = .false.
  logical             :: init_tms_proof_written = .false.

  interface
     subroutine trb_mtn_stress_init_codon(oro_in_c, z0fac_in_c, karman_in_c, gravit_in_c, rair_in_c, &
          orocnst_p, z0fac_p, karman_p, gravit_p, rair_p) bind(c, name="trb_mtn_stress_init_codon")
       use iso_c_binding, only: c_double, c_ptr
       real(c_double), value :: oro_in_c, z0fac_in_c, karman_in_c, gravit_in_c, rair_in_c
       type(c_ptr), value :: orocnst_p, z0fac_p, karman_p, gravit_p, rair_p
     end subroutine trb_mtn_stress_init_codon

     subroutine init_tms_codon(oro_in_c, z0fac_in_c, karman_in_c, gravit_in_c, rair_in_c, &
          orocnst_p, z0fac_p, karman_p, gravit_p, rair_p) bind(c, name="init_tms_codon")
       use iso_c_binding, only: c_double, c_ptr
       real(c_double), value :: oro_in_c, z0fac_in_c, karman_in_c, gravit_in_c, rair_in_c
       type(c_ptr), value :: orocnst_p, z0fac_p, karman_p, gravit_p, rair_p
     end subroutine init_tms_codon

     subroutine trb_mtn_stress_compute_codon(pcols_c, pver_c, ncol_c, orocnst_c, z0fac_c, karman_c, gravit_c, rair_c, &
          u_p, v_p, t_p, pmid_p, exner_p, zm_p, sgh_p, landfrac_p, ksrf_p, taux_p, tauy_p) &
          bind(c, name="trb_mtn_stress_compute_codon")
       use iso_c_binding, only: c_double, c_int64_t, c_ptr
       integer(c_int64_t), value :: pcols_c, pver_c, ncol_c
       real(c_double), value :: orocnst_c, z0fac_c, karman_c, gravit_c, rair_c
       type(c_ptr), value :: u_p, v_p, t_p, pmid_p, exner_p, zm_p
       type(c_ptr), value :: sgh_p, landfrac_p, ksrf_p, taux_p, tauy_p
     end subroutine trb_mtn_stress_compute_codon
  end interface

contains

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_tms_select_impl()

    use cam_logfile, only: iulog
    use spmd_utils,  only: masterproc

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (init_tms_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('TRB_MTN_STRESS_INIT_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_init_tms_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_init_tms_impl = .false.
    end if

    init_tms_impl_selected = .true.

    if (masterproc) then
       if (use_native_init_tms_impl) then
          write(iulog,*) 'trb_mtn_stress_init implementation = native'
       else
          write(iulog,*) 'trb_mtn_stress_init implementation = codon'
       end if
    end if

  end subroutine init_tms_select_impl

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_tms_proof_once()

    use cam_logfile, only: iulog
    use spmd_utils,  only: masterproc

    if (init_tms_proof_written) return
    init_tms_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'trb_mtn_stress_init direct = codon'
    end if

  end subroutine init_tms_proof_once

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine tms_select_impl()

    use cam_logfile, only: iulog
    use spmd_utils,  only: masterproc

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (tms_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('TRB_MTN_STRESS_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_tms_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_tms_impl = .false.
    end if

    tms_impl_selected = .true.

    if (masterproc) then
       if (use_native_tms_impl) then
          write(iulog,*) 'trb_mtn_stress implementation = native'
       else
          write(iulog,*) 'trb_mtn_stress implementation = codon'
       end if
    end if

  end subroutine tms_select_impl

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine tms_proof_once()

    use cam_logfile, only: iulog
    use spmd_utils,  only: masterproc

    if (tms_proof_written) return
    tms_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'trb_mtn_stress compute_tms entered (surface drag helper = codon)'
    end if

  end subroutine tms_proof_once

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine init_tms( kind, oro_in, z0fac_in, karman_in, gravit_in, rair_in, &
       errstring)

    use iso_c_binding, only: c_double, c_loc

    integer, intent(in) :: kind

    real(r8), intent(in) :: oro_in, z0fac_in, karman_in, gravit_in, rair_in

    character(len=*), intent(out) :: errstring

    errstring = ' '

    if ( kind /= r8 ) then
       errstring = 'inconsistent KIND of reals passed to init_tms'
       return
    endif

    call init_tms_select_impl()

    if (use_native_init_tms_impl) then
       orocnst  = oro_in
       z0fac    = z0fac_in
       karman   = karman_in
       gravit   = gravit_in
       rair     = rair_in
    else
       call init_tms_proof_once()
       call init_tms_codon(real(oro_in, c_double), real(z0fac_in, c_double), real(karman_in, c_double), &
            real(gravit_in, c_double), real(rair_in, c_double), c_loc(orocnst), c_loc(z0fac), c_loc(karman), &
            c_loc(gravit), c_loc(rair))
    end if
    
  end subroutine init_tms

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine compute_tms( pcols    , pver    , ncol    ,                     &
                          u        , v       , t       , pmid    , exner   , &
                          zm       , sgh     , ksrf    , taux    , tauy    , &
                          landfrac )

    use iso_c_binding, only: c_double, c_int64_t, c_loc

    integer,  intent(in)  :: pcols
    integer,  intent(in)  :: pver
    integer,  intent(in)  :: ncol
    real(r8), target, intent(in)  :: u(pcols,pver), v(pcols,pver), t(pcols,pver)
    real(r8), target, intent(in)  :: pmid(pcols,pver), exner(pcols,pver), zm(pcols,pver)
    real(r8), target, intent(in)  :: sgh(pcols), landfrac(pcols)
    real(r8), target, intent(out) :: ksrf(pcols), taux(pcols), tauy(pcols)

    call tms_select_impl()

    if (use_native_tms_impl) then
       call compute_tms_native(pcols, pver, ncol, u, v, t, pmid, exner, zm, sgh, ksrf, taux, tauy, landfrac)
       return
    end if

    call tms_proof_once()
    call trb_mtn_stress_compute_codon(int(pcols, c_int64_t), int(pver, c_int64_t), int(ncol, c_int64_t), &
         real(orocnst, c_double), real(z0fac, c_double), real(karman, c_double), real(gravit, c_double), &
         real(rair, c_double), c_loc(u(1,1)), c_loc(v(1,1)), c_loc(t(1,1)), c_loc(pmid(1,1)), c_loc(exner(1,1)), &
         c_loc(zm(1,1)), c_loc(sgh(1)), c_loc(landfrac(1)), c_loc(ksrf(1)), c_loc(taux(1)), c_loc(tauy(1)))

  end subroutine compute_tms

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine compute_tms_native( pcols    , pver    , ncol    ,                     &
                          u        , v       , t       , pmid    , exner   , &
                          zm       , sgh     , ksrf    , taux    , tauy    , & 
                          landfrac )

    !------------------------------------------------------------------------------ !
    ! Turbulent mountain stress parameterization                                    !  
    !                                                                               !
    ! Returns surface drag coefficient and stress associated with subgrid mountains !
    ! For points where the orographic variance is small ( including ocean ),        !
    ! the returned surface drag coefficient and stress is zero.                     !
    !                                                                               !
    ! Lastly arranged : Sungsu Park. Jan. 2010.                                     !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments ! 
    ! ---------------------- !

    integer,  intent(in)  :: pcols                 ! Number of columns dimensioned
    integer,  intent(in)  :: pver                  ! Number of model layers
    integer,  intent(in)  :: ncol                  ! Number of columns actually used

    real(r8), intent(in)  :: u(pcols,pver)         ! Layer mid-point zonal wind [ m/s ]
    real(r8), intent(in)  :: v(pcols,pver)         ! Layer mid-point meridional wind [ m/s ]
    real(r8), intent(in)  :: t(pcols,pver)         ! Layer mid-point temperature [ K ]
    real(r8), intent(in)  :: pmid(pcols,pver)      ! Layer mid-point pressure [ Pa ]
    real(r8), intent(in)  :: exner(pcols,pver)     ! Layer mid-point exner function [ no unit ]
    real(r8), intent(in)  :: zm(pcols,pver)        ! Layer mid-point height [ m ]
    real(r8), intent(in)  :: sgh(pcols)            ! Standard deviation of orography [ m ]
    real(r8), intent(in)  :: landfrac(pcols)       ! Land fraction [ fraction ]
    
    real(r8), intent(out) :: ksrf(pcols)           ! Surface drag coefficient [ kg/s/m2 ]
    real(r8), intent(out) :: taux(pcols)           ! Surface zonal      wind stress [ N/m2 ]
    real(r8), intent(out) :: tauy(pcols)           ! Surface meridional wind stress [ N/m2 ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer  :: i                                  ! Loop index
    integer  :: kb, kt                             ! Bottom and top of source region
    
    real(r8) :: horo                               ! Orographic height [ m ]
    real(r8) :: z0oro                              ! Orographic z0 for momentum [ m ]
    real(r8) :: dv2                                ! (delta v)**2 [ m2/s2 ]
    real(r8) :: ri                                 ! Richardson number [ no unit ]
    real(r8) :: stabfri                            ! Instability function of Richardson number [ no unit ]
    real(r8) :: rho                                ! Density [ kg/m3 ]
    real(r8) :: cd                                 ! Drag coefficient [ no unit ]
    real(r8) :: vmag                               ! Velocity magnitude [ m /s ]

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
       
    do i = 1, ncol

     ! determine subgrid orgraphic height ( mean to peak )

       horo = orocnst * sgh(i)

     ! No mountain stress if horo is too small

       if( horo < horomin ) then

           ksrf(i) = 0._r8
           taux(i) = 0._r8
           tauy(i) = 0._r8

       else

         ! Determine z0m for orography

           z0oro = min( z0fac * horo, z0max )

         ! Calculate neutral drag coefficient

           cd = ( karman / log( ( zm(i,pver) + z0oro ) / z0oro) )**2

         ! Calculate the Richardson number over the lowest 2 layers

           kt  = pver - 1
           kb  = pver
           dv2 = max( ( u(i,kt) - u(i,kb) )**2 + ( v(i,kt) - v(i,kb) )**2, dv2min )

         ! Modification : Below computation of Ri is wrong. Note that 'Exner' function here is
         !                inverse exner function. Here, exner function is not multiplied in
         !                the denominator. Also, we should use moist Ri not dry Ri.
         !                Also, this approach using the two lowest model layers can be potentially
         !                sensitive to the vertical resolution.  
         ! OK. I only modified the part associated with exner function.

           ri  = 2._r8 * gravit * ( t(i,kt) * exner(i,kt) - t(i,kb) * exner(i,kb) ) * ( zm(i,kt) - zm(i,kb) ) &
                                / ( ( t(i,kt) * exner(i,kt) + t(i,kb) * exner(i,kb) ) * dv2 )

         ! ri  = 2._r8 * gravit * ( t(i,kt) * exner(i,kt) - t(i,kb) * exner(i,kb) ) * ( zm(i,kt) - zm(i,kb) ) &
         !                      / ( ( t(i,kt) + t(i,kb) ) * dv2 )

         ! Calculate the instability function and modify the neutral drag cofficient.
         ! We should probably follow more elegant approach like Louis et al (1982) or Bretherton and Park (2009) 
         ! but for now we use very crude approach : just 1 for ri < 0, 0 for ri > 1, and linear ramping.

           stabfri = max( 0._r8, min( 1._r8, 1._r8 - ri ) )
           cd      = cd * stabfri

         ! Compute density, velocity magnitude and stress using bottom level properties

           rho     = pmid(i,pver) / ( rair * t(i,pver) ) 
           vmag    = sqrt( u(i,pver)**2 + v(i,pver)**2 )
           ksrf(i) = rho * cd * vmag * landfrac(i)
           taux(i) = -ksrf(i) * u(i,pver)
           tauy(i) = -ksrf(i) * v(i,pver)

       end if

    end do
    
    return
  end subroutine compute_tms_native

end module trb_mtn_stress
