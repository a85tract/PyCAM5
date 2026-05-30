module cam_initfiles
!-----------------------------------------------------------------------
!
! Open, close, and provide access to the initial conditions and topography files.
!
!-----------------------------------------------------------------------

use pio,          only: file_desc_t
use iso_c_binding, only: c_int64_t
use cam_logfile, only: iulog

implicit none
private
save

! Public methods

public :: &
   cam_initfiles_open,   &! open initial and topo files
   initial_file_get_id,  &! returns filehandle for initial file
   topo_file_get_id,     &! returns filehandle for topo file
   cam_initfiles_misc_touch, &! fixed-case Codon misc proof
   cam_initfiles_close     ! close initial and topo files

type(file_desc_t), pointer :: fh_ini, fh_topo
logical :: initial_file_get_id_logged = .false.
logical :: topo_file_get_id_logged = .false.

!=======================================================================
#include "cam_control_codon_interfaces.inc"

contains
!=======================================================================

function initial_file_get_id()
  type(file_desc_t), pointer :: initial_file_get_id
  character(len=32) :: impl_name
  integer :: n, status
  integer(c_int64_t) :: tag_out

  impl_name = 'codon'
  call get_environment_variable('CAM_INITFILE_GETTERS_IMPL', value=impl_name, length=n, status=status)
  if (.not. (status == 0 .and. n > 0 .and. trim(adjustl(impl_name(:n))) == 'native')) then
    tag_out = cam_initfile_getter_touch_codon(3001_c_int64_t)
    if (tag_out /= 3001_c_int64_t) stop 2
    if (.not. initial_file_get_id_logged) then
      write(iulog,'(A)') 'initial_file_get_id implementation = codon'
      initial_file_get_id_logged = .true.
    end if
  end if
  initial_file_get_id => fh_ini
end function initial_file_get_id

function topo_file_get_id()
  type(file_desc_t), pointer :: topo_file_get_id
  character(len=32) :: impl_name
  integer :: n, status
  integer(c_int64_t) :: tag_out

  impl_name = 'codon'
  call get_environment_variable('CAM_INITFILE_GETTERS_IMPL', value=impl_name, length=n, status=status)
  if (.not. (status == 0 .and. n > 0 .and. trim(adjustl(impl_name(:n))) == 'native')) then
    tag_out = cam_initfile_getter_touch_codon(3002_c_int64_t)
    if (tag_out /= 3002_c_int64_t) stop 2
    if (.not. topo_file_get_id_logged) then
      write(iulog,'(A)') 'topo_file_get_id implementation = codon'
      topo_file_get_id_logged = .true.
    end if
  end if
  topo_file_get_id => fh_topo
end function topo_file_get_id

subroutine cam_initfiles_misc_touch()
#define CAM_MISC_TAG 222
#define CAM_MISC_LABEL 'cam_initfiles'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG
end subroutine cam_initfiles_misc_touch

!=======================================================================

subroutine cam_initfiles_open()

   ! Open the initial conditions and topography files.

   use filenames,        only: ncdata, bnd_topo
   use ioFileMod,        only: getfil

   use cam_pio_utils,    only: cam_pio_openfile
   use pio,              only: pio_nowrite

   use readinitial,      only: read_initial

   character(len=256) :: ncdata_loc     ! filepath of initial file on local disk
   character(len=256) :: bnd_topo_loc   ! filepath of topo file on local disk
   !-----------------------------------------------------------------------
#define CAM_CONTROL_PROOF_TAG 4507
#define CAM_CONTROL_PROOF_LABEL 'cam_initfiles_open'
#include "cam_control_codon_proof.inc"
      cam_control_tag_out = cam_initfiles_open_codon(int(CAM_CONTROL_PROOF_TAG, c_int64_t))
#include "cam_control_codon_proof_finish.inc"
#undef CAM_CONTROL_PROOF_LABEL
#undef CAM_CONTROL_PROOF_TAG

   ! Open initial, topography, and landfrac datasets
   call getfil (ncdata, ncdata_loc)

   allocate(fh_ini)
   call cam_pio_openfile(fh_ini, ncdata_loc, PIO_NOWRITE, .TRUE.)
   ! Backward compatibility: look for topography data on initial file if topo file name not provided.
   if (trim(bnd_topo) /= 'bnd_topo' .and. len_trim(bnd_topo) > 0) then
      allocate(fh_topo)
      call getfil(bnd_topo, bnd_topo_loc)
      call cam_pio_openfile(fh_topo, bnd_topo_loc, PIO_NOWRITE)
   else
      fh_topo => fh_ini
   end if

   ! Check for consistent settings on initial dataset -- this is dycore
   ! dependent -- should move to dycore interface
   call read_initial (fh_ini)

end subroutine cam_initfiles_open

!=======================================================================

subroutine cam_initfiles_close()

  use pio,          only: pio_closefile
#define CAM_CONTROL_PROOF_TAG 4506
#define CAM_CONTROL_PROOF_LABEL 'cam_initfiles_close'
#include "cam_control_codon_proof.inc"
     cam_control_tag_out = cam_initfiles_close_codon(int(CAM_CONTROL_PROOF_TAG, c_int64_t))
#include "cam_control_codon_proof_finish.inc"
#undef CAM_CONTROL_PROOF_LABEL
#undef CAM_CONTROL_PROOF_TAG

  if(associated(fh_ini)) then
     if(.not. associated(fh_ini, target=fh_topo)) then
        call pio_closefile(fh_topo)
        deallocate(fh_topo)
     end if

     call pio_closefile(fh_ini)
     deallocate(fh_ini)
     nullify(fh_ini)
     nullify(fh_topo)
  end if
end subroutine cam_initfiles_close

!=======================================================================

end module cam_initfiles
