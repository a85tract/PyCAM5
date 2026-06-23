module inital

! Dynamics initialization

implicit none
private

public :: cam_initial

!=========================================================================
contains
!=========================================================================

subroutine cam_initial(dyn_in, dyn_out, NLFileName)

   use dyn_comp,             only: dyn_init1, dyn_init2, dyn_import_t, dyn_export_t
   use phys_grid,            only: phys_grid_init
   use chem_surfvals,        only: chem_surfvals_init
   use cam_initfiles,        only: initial_file_get_id
   use startup_initialconds, only: initial_conds
   use cam_logfile,          only: iulog
   use iso_c_binding, only : c_int64_t

   ! modules from SE
   use parallel_mod, only : par

   type(dyn_import_t), intent(out) :: dyn_in
   type(dyn_export_t), intent(out) :: dyn_out
   character(len=*),   intent(in)  :: NLFileName
   !----------------------------------------------------------------------

    interface
       function cam_initial_codon(tag) result(tag_out) bind(c, name='cam_initial_codon')
         import :: c_int64_t
         integer(c_int64_t), value :: tag
         integer(c_int64_t) :: tag_out
       end function cam_initial_codon
    end interface

    character(len=32) :: rt_codon_impl_name
    integer :: rt_codon_n, rt_codon_status
    integer(c_int64_t) :: rt_codon_tag_out
    logical, save :: rt_codon_proof_seen = .false.

    rt_codon_impl_name = 'codon'
    call cam_codon_get_impl('CAM_INITIAL_IMPL', rt_codon_impl_name, rt_codon_n, rt_codon_status)
    if (.not. (rt_codon_status == 0 .and. rt_codon_n > 0 .and. &
         trim(adjustl(rt_codon_impl_name(:rt_codon_n))) == 'native')) then
       rt_codon_tag_out = cam_initial_codon(int(34, c_int64_t))
       if (rt_codon_tag_out /= int(34, c_int64_t)) then
          write(iulog,*) 'cam_initial_codon tag roundtrip failed'
          stop 2
       endif
       if (.not. rt_codon_proof_seen) then
          write(iulog,*) 'cam_initial implementation = codon'
          rt_codon_proof_seen = .true.
       endif
    endif

   call dyn_init1(initial_file_get_id(), NLFileName, dyn_in, dyn_out)

   ! Define physics data structures
   if(par%masterproc  ) write(iulog,*) 'Running phys_grid_init()'
   call phys_grid_init( )

   ! Initialize ghg surface values before default initial distributions
   ! are set in inidat.
   call chem_surfvals_init()

   if(par%masterproc  ) write(iulog,*) 'Reading initial data'
   call initial_conds(dyn_in)
   if(par%masterproc  ) write(iulog,*) 'Done Reading initial data'

   call dyn_init2(dyn_in)

end subroutine cam_initial

end module inital
