module geopotential_temp

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: geopotential_temp_run

contains

  !> \section arg_table_geopotential_temp_run Argument Table
  !! \htmlinclude geopotential_temp_run.html
  !!
  subroutine geopotential_temp_run(num_layers, num_interfaces, finite_volume_dynamics, &
       piln, pmln, pint, pmid, pdel, rpdel, dry_static_energy, specific_humidity, &
       surface_geopotential, rair, gravit, cpair, zvir, temperature, zi, zm, ncol, &
       errcode, errmsg)
    integer, intent(in) :: num_layers, num_interfaces, ncol
    logical, intent(in) :: finite_volume_dynamics
    real(r8), intent(in) :: piln(:,:), pmln(:,:), pint(:,:), pmid(:,:)
    real(r8), intent(in) :: pdel(:,:), rpdel(:,:), dry_static_energy(:,:)
    real(r8), intent(in) :: specific_humidity(:,:), surface_geopotential(:)
    real(r8), intent(in) :: rair(:,:), gravit, cpair(:,:), zvir(:,:)
    real(r8), intent(out) :: temperature(:,:), zi(:,:), zm(:,:)
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    integer :: i, k
    real(r8) :: hkk(ncol), hkl(ncol), rog(ncol,num_layers)
    real(r8) :: tv, tvfac

    errcode = 0
    errmsg = ''
    if (num_interfaces /= num_layers + 1) then
       errcode = 1
       errmsg = 'geopotential_temp_run: interface dimension is not layer dimension + 1'
       return
    end if

    rog(:ncol,:) = rair(:ncol,:) / gravit

    do i = 1, ncol
       zi(i,num_interfaces) = 0.0_r8
    end do

    do k = num_layers, 1, -1
       if (finite_volume_dynamics) then
          do i = 1, ncol
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1._r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          end do
       else
          do i = 1, ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          end do
       end if

       do i = 1, ncol
          tvfac = 1._r8 + zvir(i,k) * specific_humidity(i,k)
          tv = (dry_static_energy(i,k) - surface_geopotential(i) - &
               gravit*zi(i,k+1)) / ((cpair(i,k) / tvfac) + rair(i,k)*hkk(i))
          temperature(i,k) = tv / tvfac
          zm(i,k) = zi(i,k+1) + rog(i,k) * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog(i,k) * tv * hkl(i)
       end do
    end do

    ! pmln is retained in the interface because it is part of the CAM
    ! geopotential_dse public contract, although the PI-atm algorithm does not use it.
  end subroutine geopotential_temp_run

end module geopotential_temp
