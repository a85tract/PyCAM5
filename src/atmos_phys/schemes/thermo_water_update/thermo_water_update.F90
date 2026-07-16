module thermo_water_update

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: thermo_water_update_run

contains

  !> \section arg_table_thermo_water_update_run Argument Table
  !! \htmlinclude thermo_water_update_run.html
  !!
  subroutine thermo_water_update_run(psetcols, standard_columns, num_layers, &
       first_chunk, last_chunk, cpair, rair, cpairv, rairv, cpairv_local, &
       rairv_local, errcode, errmsg)
    integer, intent(in) :: psetcols, standard_columns, num_layers
    integer, intent(in) :: first_chunk, last_chunk
    real(r8), intent(in) :: cpair, rair
    real(r8), intent(in) :: cpairv(:,:,:), rairv(:,:,:)
    real(r8), intent(out) :: cpairv_local(psetcols,num_layers,first_chunk:last_chunk)
    real(r8), intent(out) :: rairv_local(psetcols,num_layers,first_chunk:last_chunk)
    integer, intent(out) :: errcode
    character(len=512), intent(out) :: errmsg

    errcode = 0
    errmsg = ''

    if (psetcols == standard_columns) then
       cpairv_local(:,:,:) = cpairv(:,:,:)
    else if (psetcols > standard_columns .and. all(cpairv(:,:,:) == cpair)) then
       cpairv_local(:,:,:) = cpair
    else
       errcode = 1
       errmsg = 'thermo_water_update_run: variable cpair is not supported with subcolumns'
       return
    end if

    if (psetcols == standard_columns) then
       rairv_local(:,:,:) = rairv(:,:,:)
    else if (psetcols > standard_columns .and. all(rairv(:,:,:) == rair)) then
       rairv_local(:,:,:) = rair
    else
       errcode = 1
       errmsg = 'thermo_water_update_run: variable rair is not supported with subcolumns'
       return
    end if
  end subroutine thermo_water_update_run

end module thermo_water_update
