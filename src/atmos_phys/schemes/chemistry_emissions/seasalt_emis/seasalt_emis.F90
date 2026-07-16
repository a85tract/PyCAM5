module ap_seasalt_emis_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8
  use ap_sslt_sections_private, only : nsections, ap_sslt_section_fluxes, Dg, rdry

  implicit none
  private

  public :: seasalt_emis_run

contains


  !> \section arg_table_seasalt_emis_run Argument Table
  !! \htmlinclude seasalt_emis_run.html
  !!
  subroutine seasalt_emis_run(num_columns, declared_columns, num_constituents, &
       num_seasalt_bins, num_seasalt_tracers, u10_power, surface_temperature, &
       ocean_fraction, constituent_flux, seasalt_indices, size_range_low, &
       size_range_high, emission_scale, pi_value, seasalt_density, errmsg, errflg)
    integer, intent(in) :: num_columns
    integer, intent(in) :: declared_columns
    integer, intent(in) :: num_constituents
    integer, intent(in) :: num_seasalt_bins
    integer, intent(in) :: num_seasalt_tracers
    real(r8), intent(in) :: u10_power(declared_columns)
    real(r8), intent(in) :: surface_temperature(declared_columns)
    real(r8), intent(in) :: ocean_fraction(declared_columns)
    real(r8), intent(inout) :: constituent_flux(declared_columns,num_constituents)
    integer, intent(in) :: seasalt_indices(num_seasalt_tracers)
    real(r8), intent(in) :: size_range_low(num_seasalt_bins)
    real(r8), intent(in) :: size_range_high(num_seasalt_bins)
    real(r8), intent(in) :: emission_scale
    real(r8), intent(in) :: pi_value
    real(r8), intent(in) :: seasalt_density
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: mn, mm, ibin, i
    real(r8) :: fi(num_columns,nsections)

    errmsg = ''
    errflg = 0

    fi(:num_columns,:nsections) = ap_sslt_section_fluxes( &
         surface_temperature, u10_power, num_columns)

    do ibin = 1,num_seasalt_bins
       mm = seasalt_indices(ibin)
       mn = seasalt_indices(num_seasalt_bins+ibin)

       if (mn > 0) then
          do i = 1,nsections
             if (Dg(i) .ge. size_range_low(ibin) .and. &
                  Dg(i) .lt. size_range_high(ibin)) then
                constituent_flux(:num_columns,mn) = &
                     constituent_flux(:num_columns,mn) &
                     + fi(:num_columns,i)*ocean_fraction(:num_columns)*emission_scale
             endif
          end do
       endif

       constituent_flux(:num_columns,mm) = 0.0_r8
       do i = 1,nsections
          if (Dg(i) .ge. size_range_low(ibin) .and. &
               Dg(i) .lt. size_range_high(ibin)) then
             constituent_flux(:num_columns,mm) = constituent_flux(:num_columns,mm) &
                  + fi(:num_columns,i)*ocean_fraction(:num_columns)*emission_scale &
                  *4._r8/3._r8*pi_value*rdry(i)**3*seasalt_density
          endif
       end do
    end do
  end subroutine seasalt_emis_run

end module ap_seasalt_emis_scheme
