!--------------------------------------------------------------------------------
! Manages the absorber concentrations in the layers RRTMG operates 
! including an extra layer over the model if needed.
!
! Creator: Francis Vitt 
! 9 May 2011
!--------------------------------------------------------------------------------
module rrtmg_state

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, pverp
  use cam_logfile,     only: iulog
  use spmd_utils,      only: masterproc

  implicit none
  private
  save
  
  public :: rrtmg_state_t
  public :: rrtmg_state_init
  public :: rrtmg_state_create
  public :: rrtmg_state_update
  public :: rrtmg_state_destroy
  public :: num_rrtmg_levs

  type rrtmg_state_t

     real(r8), allocatable :: h2ovmr(:,:)   ! h2o volume mixing ratio
     real(r8), allocatable :: o3vmr(:,:)    ! o3 volume mixing ratio
     real(r8), allocatable :: co2vmr(:,:)   ! co2 volume mixing ratio 
     real(r8), allocatable :: ch4vmr(:,:)   ! ch4 volume mixing ratio 
     real(r8), allocatable :: o2vmr(:,:)    ! o2  volume mixing ratio 
     real(r8), allocatable :: n2ovmr(:,:)   ! n2o volume mixing ratio 
     real(r8), allocatable :: cfc11vmr(:,:) ! cfc11 volume mixing ratio
     real(r8), allocatable :: cfc12vmr(:,:) ! cfc12 volume mixing ratio
     real(r8), allocatable :: cfc22vmr(:,:) ! cfc22 volume mixing ratio
     real(r8), allocatable :: ccl4vmr(:,:)  ! ccl4 volume mixing ratio

     real(r8), allocatable :: pmidmb(:,:)   ! Level pressure (hPa)
     real(r8), allocatable :: pintmb(:,:)   ! Model interface pressure (hPa)
     real(r8), allocatable :: tlay(:,:)     ! mid point temperature
     real(r8), allocatable :: tlev(:,:)     ! interface temperature

  end type rrtmg_state_t

  integer :: num_rrtmg_levs ! number of pressure levels greate than 1.e-4_r8 mbar

  real(r8), parameter :: amdw = 1.607793_r8    ! Molecular weight of dry air / water vapor
  real(r8), parameter :: amdc = 0.658114_r8    ! Molecular weight of dry air / carbon dioxide
  real(r8), parameter :: amdo = 0.603428_r8    ! Molecular weight of dry air / ozone
  real(r8), parameter :: amdm = 1.805423_r8    ! Molecular weight of dry air / methane
  real(r8), parameter :: amdn = 0.658090_r8    ! Molecular weight of dry air / nitrous oxide
  real(r8), parameter :: amdo2 = 0.905140_r8   ! Molecular weight of dry air / oxygen
  real(r8), parameter :: amdc1 = 0.210852_r8   ! Molecular weight of dry air / CFC11
  real(r8), parameter :: amdc2 = 0.239546_r8   ! Molecular weight of dry air / CFC12

  logical :: use_native_rrtmg_state_impl = .false.
  logical :: rrtmg_state_impl_selected = .false.
  logical :: rrtmg_state_entered_logged = .false.
  logical :: rrtmg_state_create_logged = .false.

contains

!--------------------------------------------------------------------------------
! sets the number of model levels RRTMG operates
!--------------------------------------------------------------------------------
  subroutine rrtmg_state_init

    use ref_pres,       only : pref_edge
    use iso_c_binding,  only : c_int64_t, c_loc, c_ptr
    implicit none

    interface
       function rrtmg_state_init_codon(pverp_c, pref_edge_p) result(num_rrtmg_levs_c) &
            bind(c, name="rrtmg_state_init_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: pverp_c
         type(c_ptr), value :: pref_edge_p
         integer(c_int64_t) :: num_rrtmg_levs_c
       end function rrtmg_state_init_codon
    end interface

    call rrtmg_state_select_impl()

    ! The following cuts off RRTMG at roughly the point where it becomes
    ! invalid due to low pressure.
    if (use_native_rrtmg_state_impl) then
       num_rrtmg_levs = count( pref_edge(:) > 1._r8 ) ! pascals (1.e-2 mbar)
    else
       call rrtmg_state_log_init()
       num_rrtmg_levs = int(rrtmg_state_init_codon(int(pverp, c_int64_t), c_loc(pref_edge(1))))
    end if

  end subroutine rrtmg_state_init
  
!--------------------------------------------------------------------------------
! creates (alloacates) an rrtmg_state object
!--------------------------------------------------------------------------------

  function rrtmg_state_create( pstate, cam_in ) result( rstate )
    use physics_types,    only: physics_state
    use camsrfexch,       only: cam_in_t
    use physconst,        only: stebol
    use iso_c_binding,    only: c_double, c_int64_t, c_loc, c_ptr

    implicit none

    type(physics_state), intent(in), target :: pstate
    type(cam_in_t),      intent(in), target :: cam_in

    type(rrtmg_state_t), pointer  :: rstate

    real(r8) dy                   ! Temporary layer pressure thickness
    real(r8) :: tint(pcols,pverp)    ! Model interface temperature
    integer  :: ncol, i, kk, k

    interface
       subroutine rrtmg_state_create_codon(ncol_c, pcols_c, pver_c, pverp_c, num_rrtmg_levs_c, stebol_c, &
            t_p, lnpint_p, lnpmid_p, pmid_p, pint_p, lwup_p, pmidmb_p, pintmb_p, tlay_p, tlev_p) &
            bind(c, name="rrtmg_state_create_codon")
         use iso_c_binding, only: c_double, c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c, num_rrtmg_levs_c
         real(c_double), value :: stebol_c
         type(c_ptr), value :: t_p, lnpint_p, lnpmid_p, pmid_p, pint_p, lwup_p
         type(c_ptr), value :: pmidmb_p, pintmb_p, tlay_p, tlev_p
       end subroutine rrtmg_state_create_codon
    end interface

    allocate( rstate )

    allocate( rstate%h2ovmr(pcols,num_rrtmg_levs) )
    allocate( rstate%o3vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%co2vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%ch4vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%o2vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%n2ovmr(pcols,num_rrtmg_levs) )
    allocate( rstate%cfc11vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%cfc12vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%cfc22vmr(pcols,num_rrtmg_levs) )
    allocate( rstate%ccl4vmr(pcols,num_rrtmg_levs) )

    allocate( rstate%pmidmb(pcols,num_rrtmg_levs) )
    allocate( rstate%pintmb(pcols,num_rrtmg_levs+1) )
    allocate( rstate%tlay(pcols,num_rrtmg_levs) )
    allocate( rstate%tlev(pcols,num_rrtmg_levs+1) )

    ncol = pstate%ncol

    call rrtmg_state_select_impl()

    if (.not. use_native_rrtmg_state_impl) then
       call rrtmg_state_create_log()
       call rrtmg_state_log_entered()
       call rrtmg_state_create_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pver, c_int64_t), int(pverp, c_int64_t), &
            int(num_rrtmg_levs, c_int64_t), real(stebol, c_double), &
            c_loc(pstate%t(1,1)), c_loc(pstate%lnpint(1,1)), c_loc(pstate%lnpmid(1,1)), &
            c_loc(pstate%pmid(1,1)), c_loc(pstate%pint(1,1)), c_loc(cam_in%lwup(1)), &
            c_loc(rstate%pmidmb(1,1)), c_loc(rstate%pintmb(1,1)), c_loc(rstate%tlay(1,1)), c_loc(rstate%tlev(1,1)) &
       )
       return
    end if

    ! Calculate interface temperatures (following method
    ! used in radtpl for the longwave), using surface upward flux and
    ! stebol constant in mks units
    do i = 1,ncol
       tint(i,1) = pstate%t(i,1)
       tint(i,pverp) = sqrt(sqrt(cam_in%lwup(i)/stebol))
       do k = 2,pver
          dy = (pstate%lnpint(i,k) - pstate%lnpmid(i,k)) / (pstate%lnpmid(i,k-1) - pstate%lnpmid(i,k))
          tint(i,k) = pstate%t(i,k) - dy * (pstate%t(i,k) - pstate%t(i,k-1))
       end do
    end do

    do k = 1, num_rrtmg_levs

       kk = max(k + (pverp-num_rrtmg_levs)-1,1)

       rstate%pmidmb(:ncol,k) = pstate%pmid(:ncol,kk) * 1.e-2_r8
       rstate%pintmb(:ncol,k) = pstate%pint(:ncol,kk) * 1.e-2_r8

       rstate%tlay(:ncol,k) = pstate%t(:ncol,kk)
       rstate%tlev(:ncol,k) = tint(:ncol,kk)

    enddo

    ! bottom interface
    rstate%pintmb(:ncol,num_rrtmg_levs+1) = pstate%pint(:ncol,pverp) * 1.e-2_r8 ! mbar
    rstate%tlev(:ncol,num_rrtmg_levs+1) = tint(:ncol,pverp)

    ! top layer thickness
    if (num_rrtmg_levs==pverp) then
       rstate%pmidmb(:ncol,1) = 0.5_r8 * rstate%pintmb(:ncol,2) 
       rstate%pintmb(:ncol,1) = 1.e-4_r8 ! mbar
    endif

  endfunction rrtmg_state_create

!--------------------------------------------------------------------------------
! updates the concentration fields
!--------------------------------------------------------------------------------
  subroutine rrtmg_state_update(pstate,pbuf,icall,rstate)
    use physics_types,    only: physics_state
    use physics_buffer,   only: physics_buffer_desc
    use rad_constituents, only: rad_cnst_get_gas
    use iso_c_binding,    only: c_int64_t, c_loc, c_ptr

    implicit none

    type(physics_state), intent(in), target :: pstate
    type(physics_buffer_desc),  pointer :: pbuf(:)
    integer,             intent(in) :: icall                     ! index through climate/diagnostic radiation calls
    type(rrtmg_state_t), pointer    :: rstate

    real(r8), pointer, dimension(:,:) :: sp_hum ! specific humidity
    real(r8), pointer, dimension(:,:) :: n2o    ! nitrous oxide mass mixing ratio
    real(r8), pointer, dimension(:,:) :: ch4    ! methane mass mixing ratio
    real(r8), pointer, dimension(:,:) :: o2     ! O2 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc11  ! cfc11 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: cfc12  ! cfc12 mass mixing ratio
    real(r8), pointer, dimension(:,:) :: o3     ! Ozone mass mixing ratio
    real(r8), pointer, dimension(:,:) :: co2    ! co2   mass mixing ratio
    
    integer  :: ncol, i, kk, k

    interface
       subroutine rrtmg_state_update_codon(ncol_c, pcols_c, pverp_c, num_rrtmg_levs_c, &
            sp_hum_p, o2_p, o3_p, co2_p, n2o_p, ch4_p, cfc11_p, cfc12_p, &
            ch4vmr_p, h2ovmr_p, o3vmr_p, co2vmr_p, o2vmr_p, n2ovmr_p, cfc11vmr_p, cfc12vmr_p, &
            cfc22vmr_p, ccl4vmr_p) bind(c, name="rrtmg_state_update_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, pcols_c, pverp_c, num_rrtmg_levs_c
         type(c_ptr), value :: sp_hum_p, o2_p, o3_p, co2_p, n2o_p, ch4_p, cfc11_p, cfc12_p
         type(c_ptr), value :: ch4vmr_p, h2ovmr_p, o3vmr_p, co2vmr_p, o2vmr_p, n2ovmr_p
         type(c_ptr), value :: cfc11vmr_p, cfc12vmr_p, cfc22vmr_p, ccl4vmr_p
       end subroutine rrtmg_state_update_codon
    end interface

    ncol = pstate%ncol

    ! Get specific humidity
    call rad_cnst_get_gas(icall,'H2O', pstate, pbuf, sp_hum)
    ! Get oxygen mass mixing ratio.
    call rad_cnst_get_gas(icall,'O2',  pstate, pbuf, o2)
    ! Get ozone mass mixing ratio.
    call rad_cnst_get_gas(icall,'O3',  pstate, pbuf, o3)
    ! Get CO2 mass mixing ratio
    call rad_cnst_get_gas(icall,'CO2', pstate, pbuf, co2)
    ! Get N2O mass mixing ratio
    call rad_cnst_get_gas(icall,'N2O', pstate, pbuf, n2o)
    ! Get CH4 mass mixing ratio
    call rad_cnst_get_gas(icall,'CH4', pstate, pbuf, ch4)
    ! Get CFC mass mixing ratios
    call rad_cnst_get_gas(icall,'CFC11', pstate, pbuf, cfc11)
    call rad_cnst_get_gas(icall,'CFC12', pstate, pbuf, cfc12)

    call rrtmg_state_select_impl()

    if (.not. use_native_rrtmg_state_impl) then
       call rrtmg_state_log_entered()
       call rrtmg_state_update_codon( &
            int(ncol, c_int64_t), int(pcols, c_int64_t), int(pverp, c_int64_t), int(num_rrtmg_levs, c_int64_t), &
            c_loc(sp_hum(1,1)), c_loc(o2(1,1)), c_loc(o3(1,1)), c_loc(co2(1,1)), &
            c_loc(n2o(1,1)), c_loc(ch4(1,1)), c_loc(cfc11(1,1)), c_loc(cfc12(1,1)), &
            c_loc(rstate%ch4vmr(1,1)), c_loc(rstate%h2ovmr(1,1)), c_loc(rstate%o3vmr(1,1)), &
            c_loc(rstate%co2vmr(1,1)), c_loc(rstate%o2vmr(1,1)), c_loc(rstate%n2ovmr(1,1)), &
            c_loc(rstate%cfc11vmr(1,1)), c_loc(rstate%cfc12vmr(1,1)), c_loc(rstate%cfc22vmr(1,1)), &
            c_loc(rstate%ccl4vmr(1,1)) &
       )
       return
    end if

    do k = 1, num_rrtmg_levs

       kk = max(k + (pverp-num_rrtmg_levs)-1,1)

       rstate%ch4vmr(:ncol,k)   = ch4(:ncol,kk) * amdm
       rstate%h2ovmr(:ncol,k)   = (sp_hum(:ncol,kk) / (1._r8 - sp_hum(:ncol,kk))) * amdw
       rstate%o3vmr(:ncol,k)    = o3(:ncol,kk) * amdo
       rstate%co2vmr(:ncol,k)   = co2(:ncol,kk) * amdc
       rstate%ch4vmr(:ncol,k)   = ch4(:ncol,kk) * amdm
       rstate%o2vmr(:ncol,k)    = o2(:ncol,kk) * amdo2
       rstate%n2ovmr(:ncol,k)   = n2o(:ncol,kk) * amdn
       rstate%cfc11vmr(:ncol,k) = cfc11(:ncol,kk) * amdc1
       rstate%cfc12vmr(:ncol,k) = cfc12(:ncol,kk) * amdc2
       rstate%cfc22vmr(:ncol,k) = 0._r8
       rstate%ccl4vmr(:ncol,k)  = 0._r8

    enddo

  end subroutine rrtmg_state_update

!--------------------------------------------------------------------------------
! de-allocates an rrtmg_state object
!--------------------------------------------------------------------------------
  subroutine rrtmg_state_destroy(rstate)

    use iso_c_binding, only : c_int64_t

    implicit none

    type(rrtmg_state_t), pointer   :: rstate
    integer(c_int64_t) :: destroy_token

    interface
       function rrtmg_state_destroy_codon_touch(value_c) result(result_c) &
            bind(c, name="rrtmg_state_destroy_codon_touch")
         use iso_c_binding, only: c_int64_t
         integer(c_int64_t), value :: value_c
         integer(c_int64_t) :: result_c
       end function rrtmg_state_destroy_codon_touch
    end interface

    call rrtmg_state_select_impl()
    if (.not. use_native_rrtmg_state_impl) then
       destroy_token = rrtmg_state_destroy_codon_touch(int(num_rrtmg_levs, c_int64_t))
       if (masterproc .and. destroy_token >= 0_c_int64_t) then
          write(iulog,*) 'rrtmg_state_destroy implementation = codon'
          call flush(iulog)
       endif
    else if (masterproc) then
       write(iulog,*) 'rrtmg_state_destroy implementation = native'
       call flush(iulog)
    end if

    deallocate(rstate%h2ovmr)
    deallocate(rstate%o3vmr)
    deallocate(rstate%co2vmr)
    deallocate(rstate%ch4vmr)
    deallocate(rstate%o2vmr)
    deallocate(rstate%n2ovmr)
    deallocate(rstate%cfc11vmr)
    deallocate(rstate%cfc12vmr)
    deallocate(rstate%cfc22vmr)
    deallocate(rstate%ccl4vmr)

    deallocate(rstate%pmidmb)
    deallocate(rstate%pintmb)
    deallocate(rstate%tlay)
    deallocate(rstate%tlev)

    deallocate( rstate )
    nullify(rstate)

  endsubroutine rrtmg_state_destroy

!--------------------------------------------------------------------------------

  subroutine rrtmg_state_select_impl()

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (rrtmg_state_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('RRTMG_STATE_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_rrtmg_state_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_rrtmg_state_impl = .false.
    end if

    rrtmg_state_impl_selected = .true.

    if (masterproc) then
       if (use_native_rrtmg_state_impl) then
          write(iulog,*) 'rrtmg_state implementation = native'
       else
          write(iulog,*) 'rrtmg_state implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine rrtmg_state_select_impl

!--------------------------------------------------------------------------------

  subroutine rrtmg_state_log_entered()

    if (rrtmg_state_entered_logged) return
    rrtmg_state_entered_logged = .true.

    if (masterproc) then
       write(iulog,*) 'rrtmg_state entered (create/update helpers = codon)'
       call flush(iulog)
    end if

  end subroutine rrtmg_state_log_entered

!--------------------------------------------------------------------------------

  subroutine rrtmg_state_create_log()

    if (rrtmg_state_create_logged) return
    rrtmg_state_create_logged = .true.

    if (masterproc) then
       write(iulog,*) 'rrtmg_state_create implementation = codon'
       call flush(iulog)
    end if

  end subroutine rrtmg_state_create_log

!--------------------------------------------------------------------------------

  subroutine rrtmg_state_log_init()

    if (masterproc) then
       write(iulog,*) 'rrtmg_state_init implementation = codon'
       call flush(iulog)
    end if

  end subroutine rrtmg_state_log_init

end module rrtmg_state
