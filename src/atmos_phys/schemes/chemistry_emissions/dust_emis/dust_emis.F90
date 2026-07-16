module ap_dust_emis_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: dust_emis_run

contains


  !> \section arg_table_dust_emis_run Argument Table
  !! \htmlinclude dust_emis_run.html
  !!
  subroutine dust_emis_run(num_columns, declared_columns, num_dust_flux_components, &
       num_constituents, num_dust_bins, num_dust_emission_tracers, dust_flux_in, &
       constituent_flux, soil_erodibility, dust_indices, dust_emission_scale, &
       dust_volume_mean_diameter, soil_erodibility_factor, pi_value, dust_density, &
       errmsg, errflg)
    integer, intent(in) :: num_columns
    integer, intent(in) :: declared_columns
    integer, intent(in) :: num_dust_flux_components
    integer, intent(in) :: num_constituents
    integer, intent(in) :: num_dust_bins
    integer, intent(in) :: num_dust_emission_tracers
    real(r8), intent(in) :: dust_flux_in(declared_columns,num_dust_flux_components)
    real(r8), intent(inout) :: constituent_flux(declared_columns,num_constituents)
    real(r8), intent(inout) :: soil_erodibility(declared_columns)
    integer, intent(in) :: dust_indices(num_dust_emission_tracers)
    real(r8), intent(in) :: dust_emission_scale(num_dust_bins)
    real(r8), intent(in) :: dust_volume_mean_diameter(num_dust_bins)
    real(r8), intent(in) :: soil_erodibility_factor
    real(r8), intent(in) :: pi_value
    real(r8), intent(in) :: dust_density
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: i, m, idst, inum
    real(r8) :: x_mton
    real(r8), parameter :: soil_erod_threshold = 0.1_r8

    errmsg = ''
    errflg = 0

    col_loop: do i = 1,num_columns
       if (soil_erodibility(i) .lt. soil_erod_threshold) soil_erodibility(i) = 0._r8

       do m = 1,num_dust_bins
          idst = dust_indices(m)

          constituent_flux(i,idst) = sum(-dust_flux_in(i,:)) &
               * dust_emission_scale(m)*soil_erodibility(i)/soil_erodibility_factor*1.15_r8

          x_mton = 6._r8 / (pi_value * dust_density &
               * (dust_volume_mean_diameter(m)**3._r8))

          inum = dust_indices(m+num_dust_bins)
          constituent_flux(i,inum) = constituent_flux(i,idst)*x_mton
       end do
    end do col_loop
  end subroutine dust_emis_run

end module ap_dust_emis_scheme
