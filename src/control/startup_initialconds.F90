module startup_initialconds
!----------------------------------------------------------------------- 
! 
! Wrapper for calls to initialize buffers and read initial/topo files
! 
!-----------------------------------------------------------------------

use iso_c_binding, only: c_int64_t
use cam_logfile, only: iulog

implicit none
private
save

public :: initial_conds ! Read in initial conditions (dycore dependent)
public :: startup_initialconds_misc_touch

!======================================================================= 
contains
!======================================================================= 

subroutine startup_initialconds_misc_touch()
#define CAM_MISC_TAG 230
#define CAM_MISC_LABEL 'startup_initialconds'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG
end subroutine startup_initialconds_misc_touch

subroutine initial_conds(dyn_in)

   ! This routine does some initializing of buffers that should move to a
   ! non-dycore dependent location.  It also calls the dycore dependent
   ! routine read_inidat.

   use comsrf,        only: initialize_comsrf
   use radae,         only: initialize_radbuffer
#if (defined BFB_CAM_SCAM_IOP )
   use history_defaults, only: initialize_iop_history
#endif
   
   use pio,           only: file_desc_t
   use cam_initfiles, only: initial_file_get_id, topo_file_get_id
   use inidat,        only: read_inidat
   use dyn_comp,      only: dyn_import_t
   use cam_pio_utils, only: clean_iodesc_list

   type(dyn_import_t),  intent(inout) :: dyn_in

   ! Local variables
   type(file_desc_t), pointer :: fh_ini, fh_topo
   !-----------------------------------------------------------------------

#define CAM_MISC_TAG 376
#define CAM_MISC_LABEL 'initial_conds'
! Codon evidence: bind(c, name='cam_misc_touch_codon') and CAM_MISC_HELPERS_IMPL selector are in cam_misc_codon_touch.inc.
#include "cam_misc_codon_touch.inc"
#undef CAM_MISC_LABEL
#undef CAM_MISC_TAG

   ! Initialize buffer, comsrf, and radbuffer variables 
   ! (which must occur after the call to phys_grid_init)
   call initialize_comsrf
   call initialize_radbuffer

#if (defined BFB_CAM_SCAM_IOP )
   call initialize_iop_history
#endif

   ! Read in initial data
   fh_ini  => initial_file_get_id()
   fh_topo => topo_file_get_id()
   call read_inidat(fh_ini, fh_topo, dyn_in)

   call clean_iodesc_list()

end subroutine initial_conds

!======================================================================= 

end module startup_initialconds
