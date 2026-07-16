module ap_set_srf_emissions_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: set_srf_emissions_run

contains


  !> \section arg_table_set_srf_emissions_run Argument Table
  !! \htmlinclude set_srf_emissions_run.html
  !!
  subroutine set_srf_emissions_run(num_columns, declared_columns, num_species, &
       num_emission_files, allocated_emission_files, max_sector_count, &
       species_index, sector_count, sector_flux, emission_scale, &
       molecular_weight, units_are_mks, mass_conversion_factor, latitude, &
       longitude, calendar_day, pi_value, twopi_value, pid2_value, dec_max_value, &
       c10h16_index, has_c10h16_emission, &
       isoprene_index, has_isoprene_emission, surface_flux, errmsg, errflg)
    integer, intent(in) :: num_columns
    integer, intent(in) :: declared_columns
    integer, intent(in) :: num_species
    integer, intent(in) :: num_emission_files
    integer, intent(in) :: allocated_emission_files
    integer, intent(in) :: max_sector_count
    integer, intent(in) :: species_index(allocated_emission_files)
    integer, intent(in) :: sector_count(allocated_emission_files)
    real(r8), intent(in) :: sector_flux(declared_columns,max_sector_count, &
         allocated_emission_files)
    real(r8), intent(in) :: emission_scale(allocated_emission_files)
    real(r8), intent(in) :: molecular_weight(allocated_emission_files)
    logical, intent(in) :: units_are_mks(allocated_emission_files)
    real(r8), intent(in) :: mass_conversion_factor
    real(r8), intent(in) :: latitude(num_columns)
    real(r8), intent(in) :: longitude(num_columns)
    real(r8), intent(in) :: calendar_day
    real(r8), intent(in) :: pi_value
    real(r8), intent(in) :: twopi_value
    real(r8), intent(in) :: pid2_value
    real(r8), intent(in) :: dec_max_value
    integer, intent(in) :: c10h16_index
    logical, intent(in) :: has_c10h16_emission
    integer, intent(in) :: isoprene_index
    logical, intent(in) :: has_isoprene_emission
    real(r8), intent(out) :: surface_flux(declared_columns,num_species)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: i, m, n, isec
    real(r8) :: factor
    real(r8) :: dayfrac
    real(r8) :: iso_off, iso_on
    logical :: polar_day, polar_night
    real(r8) :: doy_loc
    real(r8) :: sunon, sunoff
    real(r8) :: loc_angle
    real(r8) :: column_latitude
    real(r8) :: declination
    real(r8) :: tod
    real(r8) :: flux(num_columns)
    real(r8) :: mfactor
    real(r8), parameter :: dayspy = 365._r8

    errmsg = ''
    errflg = 0

    surface_flux(:,:) = 0._r8

    emis_loop: do m = 1,num_emission_files
       n = species_index(m)

       flux(:) = 0._r8
       do isec = 1,sector_count(m)
          flux(:num_columns) = flux(:num_columns) &
               + sector_flux(:num_columns,isec,m)
       end do

       flux(:num_columns) = emission_scale(m)*flux(:num_columns)

       if (units_are_mks(m)) then
          surface_flux(:num_columns,n) = surface_flux(:num_columns,n) &
               + flux(:num_columns)
       else
          mfactor = mass_conversion_factor * molecular_weight(m)
          surface_flux(:num_columns,n) = surface_flux(:num_columns,n) &
               + flux(:num_columns) * mfactor
       endif
    end do emis_loop

    doy_loc = aint(calendar_day)
    declination = dec_max_value * cos((doy_loc - 172._r8)*twopi_value/dayspy)
    tod = (calendar_day - doy_loc) + .5_r8

    do i = 1,num_columns
       polar_day = .false.
       polar_night = .false.

       loc_angle = tod * twopi_value + longitude(i)
       loc_angle = mod(loc_angle,twopi_value)
       column_latitude = latitude(i)

       if (abs(column_latitude) >= (pid2_value - abs(declination))) then
          if (sign(1._r8,declination) == sign(1._r8,column_latitude)) then
             polar_day = .true.
             sunoff = 2._r8*twopi_value
             sunon = -twopi_value
          else
             polar_night = .true.
          end if
       else
          sunoff = acos(-tan(declination)*tan(column_latitude))
          sunon = twopi_value - sunoff
       end if

       if (c10h16_index > 0) then
          if (has_c10h16_emission) then
             if (.not. polar_night .and. .not. polar_day) then
                dayfrac = sunoff / pi_value
                surface_flux(i,c10h16_index) = surface_flux(i,c10h16_index) &
                     / (.7_r8 + .3_r8*dayfrac)
                if (loc_angle >= sunoff .and. loc_angle <= sunon) then
                   surface_flux(i,c10h16_index) = surface_flux(i,c10h16_index) * .7_r8
                endif
             end if
          end if
       end if

       if (isoprene_index > 0) then
          if (has_isoprene_emission) then
             if (.not. polar_night) then
                if (polar_day) then
                   iso_off = .8_r8 * pi_value
                   iso_on = 1.2_r8 * pi_value
                else
                   iso_off = .8_r8 * sunoff
                   iso_on = 2._r8 * pi_value - iso_off
                end if
                if (loc_angle >= iso_off .and. loc_angle <= iso_on) then
                   surface_flux(i,isoprene_index) = 0._r8
                else
                   factor = loc_angle - iso_on
                   if (factor <= 0._r8) then
                      factor = factor + 2._r8*pi_value
                   end if
                   factor = factor / (2._r8*iso_off + 1.e-6_r8)
                   surface_flux(i,isoprene_index) = surface_flux(i,isoprene_index) &
                        * 2._r8 / iso_off * pi_value * (sin(pi_value*factor))**2
                end if
             else
                surface_flux(i,isoprene_index) = 0._r8
             end if
          end if
       end if
    end do
  end subroutine set_srf_emissions_run

end module ap_set_srf_emissions_scheme
