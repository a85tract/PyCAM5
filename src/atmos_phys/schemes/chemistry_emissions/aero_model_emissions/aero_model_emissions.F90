module ap_aero_model_emissions_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: aero_model_emissions_run

contains


  !> \section arg_table_aero_model_emissions_run Argument Table
  !! \htmlinclude aero_model_emissions_run.html
  !!
  subroutine aero_model_emissions_run(num_columns, declared_columns, &
       eastward_wind, northward_wind, midpoint_height, ocean_roughness_length, &
       u10_emission_factor, errmsg, errflg)
    integer, intent(in) :: num_columns
    integer, intent(in) :: declared_columns
    real(r8), intent(in) :: eastward_wind(declared_columns)
    real(r8), intent(in) :: northward_wind(declared_columns)
    real(r8), intent(in) :: midpoint_height(declared_columns)
    real(r8), intent(in) :: ocean_roughness_length
    real(r8), intent(out) :: u10_emission_factor(declared_columns)
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    errmsg = ''
    errflg = 0

    u10_emission_factor(:num_columns) = sqrt( &
         eastward_wind(:num_columns)**2+northward_wind(:num_columns)**2)

    ! Move the winds from the lowest-layer midpoint to 10 m.
    u10_emission_factor(:num_columns) = u10_emission_factor(:num_columns) &
         *log(10._r8/ocean_roughness_length) &
         /log(midpoint_height(:num_columns)/ocean_roughness_length)

    ! Gong et al. (1997) sea-salt emission wind factor.
    u10_emission_factor(:num_columns) = u10_emission_factor(:num_columns)**3.41_r8
  end subroutine aero_model_emissions_run

end module ap_aero_model_emissions_scheme
