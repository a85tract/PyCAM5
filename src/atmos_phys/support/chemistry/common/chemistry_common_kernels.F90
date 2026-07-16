! Intrinsic-only state transforms shared by the standalone chemistry schemes.
module ap_chemistry_common_kernels

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: chemistry_map_constituents
  public :: chemistry_set_mean_mass
  public :: chemistry_mmr_to_vmr
  public :: chemistry_vmr_to_mmr
  public :: chemistry_mmr_to_vmri
  public :: chemistry_h2o_to_vmr
  public :: chemistry_form_tendencies
  public :: chemistry_map_surface_flux
  public :: chemistry_clip_negative

contains

  subroutine chemistry_map_constituents(ncol, pcols, pver, pcnst,       &
       gas_pcnst, map2chm, q, mmr)
    integer, intent(in) :: ncol, pcols, pver, pcnst, gas_pcnst
    integer, intent(in) :: map2chm(pcnst)
    real(r8), intent(in) :: q(pcols,pver,pcnst)
    real(r8), intent(inout) :: mmr(pcols,pver,gas_pcnst)
    integer :: m, n

    do m = 1, pcnst
       n = map2chm(m)
       if (n > 0) mmr(:ncol,:,n) = q(:ncol,:,m)
    end do
  end subroutine chemistry_map_constituents

  subroutine chemistry_set_mean_mass(ncol, pver, gas_pcnst, fixed_mbar, &
       mwdry, id_o2, id_o, id_h, id_n, adv_mass, mmr, mbar, ierr)
    integer, intent(in) :: ncol, pver, gas_pcnst
    logical, intent(in) :: fixed_mbar
    real(r8), intent(in) :: mwdry
    integer, intent(in) :: id_o2, id_o, id_h, id_n
    real(r8), intent(in) :: adv_mass(gas_pcnst)
    real(r8), intent(in) :: mmr(ncol,pver,gas_pcnst)
    real(r8), intent(out) :: mbar(ncol,pver)
    integer, intent(out) :: ierr
    integer :: k
    real(r8) :: xn2(ncol), fn2(ncol), fo(ncol), fo2(ncol), fh(ncol)

    ierr = 0
    if (fixed_mbar) then
       mbar = mwdry
    else if (id_o2 > 0 .and. id_o > 0 .and. id_h > 0 .and. id_n > 0) then
       do k = 1, pver
          xn2 = 1._r8 - (mmr(:,k,id_o2) + mmr(:,k,id_o) + mmr(:,k,id_h))
          fn2 = .5_r8 * xn2 / adv_mass(id_n)
          fo2 = mmr(:,k,id_o2) / adv_mass(id_o2)
          fo = mmr(:,k,id_o) / adv_mass(id_o)
          fh = mmr(:,k,id_h) / adv_mass(id_h)
          mbar(:,k) = 1._r8 / (fn2 + fo2 + fo + fh)
       end do
    else
       mbar = mwdry
       ierr = 1
    end if
  end subroutine chemistry_set_mean_mass

  subroutine chemistry_mmr_to_vmr(ncol, pver, gas_pcnst, adv_mass,     &
       mmr, mbar, vmr)
    integer, intent(in) :: ncol, pver, gas_pcnst
    real(r8), intent(in) :: adv_mass(gas_pcnst)
    real(r8), intent(in) :: mmr(ncol,pver,gas_pcnst)
    real(r8), intent(in) :: mbar(ncol,pver)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)
    integer :: k, m

    do m = 1, gas_pcnst
       if (adv_mass(m) /= 0._r8) then
          do k = 1, pver
             vmr(:,k,m) = mbar(:,k) * mmr(:,k,m) / adv_mass(m)
          end do
       end if
    end do
  end subroutine chemistry_mmr_to_vmr

  subroutine chemistry_vmr_to_mmr(ncol, pver, gas_pcnst, adv_mass,     &
       vmr, mbar, mmr)
    integer, intent(in) :: ncol, pver, gas_pcnst
    real(r8), intent(in) :: adv_mass(gas_pcnst)
    real(r8), intent(in) :: vmr(ncol,pver,gas_pcnst)
    real(r8), intent(in) :: mbar(ncol,pver)
    real(r8), intent(inout) :: mmr(ncol,pver,gas_pcnst)
    integer :: k, m

    do m = 1, gas_pcnst
       if (adv_mass(m) /= 0._r8) then
          do k = 1, pver
             mmr(:,k,m) = adv_mass(m) * vmr(:,k,m) / mbar(:,k)
          end do
       end if
    end do
  end subroutine chemistry_vmr_to_mmr

  subroutine chemistry_mmr_to_vmri(ncol, pver, molecular_weight, mmr, &
       mbar, vmr)
    integer, intent(in) :: ncol, pver
    real(r8), intent(in) :: molecular_weight
    real(r8), intent(in) :: mmr(ncol,pver), mbar(ncol,pver)
    real(r8), intent(out) :: vmr(ncol,pver)
    integer :: k
    real(r8) :: reciprocal_weight

    reciprocal_weight = 1._r8 / molecular_weight
    do k = 1, pver
       vmr(:,k) = mbar(:,k) * mmr(:,k) * reciprocal_weight
    end do
  end subroutine chemistry_mmr_to_vmri

  subroutine chemistry_h2o_to_vmr(ncol, pver, h2o_molecular_weight,   &
       h2o_mmr, mbar, h2o_vmr)
    integer, intent(in) :: ncol, pver
    real(r8), intent(in) :: h2o_molecular_weight
    real(r8), intent(in) :: h2o_mmr(ncol,pver), mbar(ncol,pver)
    real(r8), intent(out) :: h2o_vmr(ncol,pver)
    integer :: k

    do k = 1, pver
       h2o_vmr(:,k) = mbar(:,k) * h2o_mmr(:,k) / h2o_molecular_weight
    end do
  end subroutine chemistry_h2o_to_vmr

  subroutine chemistry_form_tendencies(ncol, pcols, pver, pcnst,      &
       gas_pcnst, map2chm, reciprocal_timestep, mmr_old, mmr_new,    &
       mmr_tendency, qtend)
    integer, intent(in) :: ncol, pcols, pver, pcnst, gas_pcnst
    integer, intent(in) :: map2chm(pcnst)
    real(r8), intent(in) :: reciprocal_timestep
    real(r8), intent(in) :: mmr_old(pcols,pver,gas_pcnst)
    real(r8), intent(in) :: mmr_new(pcols,pver,gas_pcnst)
    real(r8), intent(out) :: mmr_tendency(pcols,pver,gas_pcnst)
    real(r8), intent(inout) :: qtend(pcols,pver,pcnst)
    integer :: m, n

    do m = 1, gas_pcnst
       mmr_tendency(:ncol,:,m) = (mmr_new(:ncol,:,m) -             &
            mmr_old(:ncol,:,m)) * reciprocal_timestep
    end do
    do m = 1, pcnst
       n = map2chm(m)
       if (n > 0) qtend(:ncol,:,m) = qtend(:ncol,:,m) +            &
            mmr_tendency(:ncol,:,n)
    end do
  end subroutine chemistry_form_tendencies

  subroutine chemistry_map_surface_flux(ncol, pcols, pcnst, gas_pcnst, &
       map2chm, species_flux, cflx, drydepflx)
    integer, intent(in) :: ncol, pcols, pcnst, gas_pcnst
    integer, intent(in) :: map2chm(pcnst)
    real(r8), intent(in) :: species_flux(pcols,gas_pcnst)
    real(r8), intent(inout) :: cflx(pcols,pcnst)
    real(r8), intent(out) :: drydepflx(pcols,pcnst)
    integer :: m, n

    drydepflx = 0._r8
    do m = 1, pcnst
       n = map2chm(m)
       if (n > 0) then
          cflx(:ncol,m) = cflx(:ncol,m) - species_flux(:ncol,n)
          drydepflx(:ncol,m) = species_flux(:ncol,n)
       end if
    end do
  end subroutine chemistry_map_surface_flux

  subroutine chemistry_clip_negative(ncol, pver, gas_pcnst, field)
    integer, intent(in) :: ncol, pver, gas_pcnst
    real(r8), intent(inout) :: field(ncol,pver,gas_pcnst)
    integer :: m

    do m = 1, gas_pcnst
       where (field(:,:,m) < 0._r8) field(:,:,m) = 0._r8
    end do
  end subroutine chemistry_clip_negative

end module ap_chemistry_common_kernels
