
module mo_chem_utls

  use iso_c_binding, only : c_int64_t, c_loc, c_ptr
  use cam_logfile,   only : iulog
  use spmd_utils,    only : masterproc

  private
  public :: get_spc_ndx, get_het_ndx, get_extfrc_ndx, get_rxt_ndx, get_inv_ndx

  save

  logical :: get_spc_ndx_codon_logged = .false.
  logical :: get_inv_ndx_codon_logged = .false.
  logical :: get_extfrc_ndx_codon_logged = .false.
  logical :: get_rxt_ndx_codon_logged = .false.

  interface
     function chem_lookup_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, list_count) result(idx_c) &
          bind(c, name="chem_lookup_name_codon")
       import :: c_int64_t, c_ptr
       integer(c_int64_t), value :: name_len, list_len, list_count
       type(c_ptr), value :: name_ascii_p, list_ascii_p
       integer(c_int64_t) :: idx_c
     end function chem_lookup_name_codon

     function chem_lookup_mapped_name_codon(name_len, name_ascii_p, list_len, list_ascii_p, map_p, list_count) result(idx_c) &
          bind(c, name="chem_lookup_mapped_name_codon")
       import :: c_int64_t, c_ptr
       integer(c_int64_t), value :: name_len, list_len, list_count
       type(c_ptr), value :: name_ascii_p, list_ascii_p, map_p
       integer(c_int64_t) :: idx_c
     end function chem_lookup_mapped_name_codon
  end interface

contains

  integer function get_spc_ndx( spc_name )
    !-----------------------------------------------------------------------
    !     ... return overall species index associated with spc_name
    !-----------------------------------------------------------------------

    use chem_mods,     only : gas_pcnst
    use mo_tracname,   only : tracnam => solsym

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: spc_name

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m
    integer :: ichar
    integer(c_int64_t), allocatable, target :: name_ascii(:)
    integer(c_int64_t), allocatable, target :: tracnam_ascii(:,:)
    integer(c_int64_t) :: idx_c

    get_spc_ndx = -1
    allocate(name_ascii(max(1, len(spc_name))))
    allocate(tracnam_ascii(len(tracnam), max(1, gas_pcnst)))
    do ichar = 1,len(spc_name)
       name_ascii(ichar) = int(iachar(spc_name(ichar:ichar)), c_int64_t)
    end do
    if (len(spc_name) == 0) name_ascii(1) = 32_c_int64_t
    do m = 1,gas_pcnst
       do ichar = 1,len(tracnam(1))
          tracnam_ascii(ichar,m) = int(iachar(tracnam(m)(ichar:ichar)), c_int64_t)
       end do
    end do
    idx_c = chem_lookup_name_codon(int(len(spc_name), c_int64_t), c_loc(name_ascii(1)), &
         int(len(tracnam), c_int64_t), c_loc(tracnam_ascii(1,1)), int(gas_pcnst, c_int64_t))
    get_spc_ndx = int(idx_c)
    call mo_chem_utls_log_direct(get_spc_ndx_codon_logged, 'get_spc_ndx direct = codon')

  end function get_spc_ndx

  integer function get_inv_ndx( invariant )
    !-----------------------------------------------------------------------
    !     ... return overall external frcing index associated with spc_name
    !-----------------------------------------------------------------------

    use chem_mods,  only : nfs, inv_lst

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: invariant

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m
    integer :: ichar
    integer(c_int64_t), allocatable, target :: name_ascii(:)
    integer(c_int64_t), allocatable, target :: inv_ascii(:,:)
    integer(c_int64_t) :: idx_c

    get_inv_ndx = -1
    if (nfs <= 0) then
       call mo_chem_utls_log_direct(get_inv_ndx_codon_logged, 'get_inv_ndx direct = codon')
       return
    end if
    allocate(name_ascii(max(1, len(invariant))))
    allocate(inv_ascii(len(inv_lst), nfs))
    do ichar = 1,len(invariant)
       name_ascii(ichar) = int(iachar(invariant(ichar:ichar)), c_int64_t)
    end do
    if (len(invariant) == 0) name_ascii(1) = 32_c_int64_t
    do m = 1,nfs
       do ichar = 1,len(inv_lst)
          inv_ascii(ichar,m) = int(iachar(inv_lst(m)(ichar:ichar)), c_int64_t)
       end do
    end do
    idx_c = chem_lookup_name_codon(int(len(invariant), c_int64_t), c_loc(name_ascii(1)), &
         int(len(inv_lst), c_int64_t), c_loc(inv_ascii(1,1)), int(nfs, c_int64_t))
    get_inv_ndx = int(idx_c)
    call mo_chem_utls_log_direct(get_inv_ndx_codon_logged, 'get_inv_ndx direct = codon')

  end function get_inv_ndx

  integer function get_het_ndx( het_name )
    !-----------------------------------------------------------------------
    !     ... return overall het process index associated with spc_name
    !-----------------------------------------------------------------------

    use gas_wetdep_opts,only : gas_wetdep_method, gas_wetdep_list, gas_wetdep_cnt

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: het_name

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    get_het_ndx=-1

    do m=1,gas_wetdep_cnt

       if( trim( het_name ) == trim( gas_wetdep_list(m) ) ) then
          get_het_ndx = get_spc_ndx( gas_wetdep_list(m) )
          return
       endif
  
    enddo

  end function get_het_ndx

  integer function get_extfrc_ndx( frc_name )
    !-----------------------------------------------------------------------
    !     ... return overall external frcing index associated with spc_name
    !-----------------------------------------------------------------------

    use chem_mods,  only : extcnt, extfrc_lst

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: frc_name

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m
    integer :: ichar
    integer(c_int64_t), allocatable, target :: name_ascii(:)
    integer(c_int64_t), allocatable, target :: extfrc_ascii(:,:)
    integer(c_int64_t) :: idx_c

    get_extfrc_ndx = -1
    if (extcnt <= 0) then
       call mo_chem_utls_log_direct(get_extfrc_ndx_codon_logged, 'get_extfrc_ndx direct = codon')
       return
    end if
    allocate(name_ascii(max(1, len(frc_name))))
    allocate(extfrc_ascii(len(extfrc_lst), extcnt))
    do ichar = 1,len(frc_name)
       name_ascii(ichar) = int(iachar(frc_name(ichar:ichar)), c_int64_t)
    end do
    if (len(frc_name) == 0) name_ascii(1) = 32_c_int64_t
    do m = 1,extcnt
       do ichar = 1,len(extfrc_lst)
          extfrc_ascii(ichar,m) = int(iachar(extfrc_lst(m)(ichar:ichar)), c_int64_t)
       end do
    end do
    idx_c = chem_lookup_name_codon(int(len(frc_name), c_int64_t), c_loc(name_ascii(1)), &
         int(len(extfrc_lst), c_int64_t), c_loc(extfrc_ascii(1,1)), int(extcnt, c_int64_t))
    get_extfrc_ndx = int(idx_c)
    call mo_chem_utls_log_direct(get_extfrc_ndx_codon_logged, 'get_extfrc_ndx direct = codon')

  end function get_extfrc_ndx

  integer function get_rxt_ndx( rxt_tag )
    !-----------------------------------------------------------------------
    !     ... return overall external frcing index associated with spc_name
    !-----------------------------------------------------------------------

    use chem_mods,  only : rxt_tag_cnt, rxt_tag_lst, rxt_tag_map

    implicit none

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: rxt_tag

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m
    integer :: ichar
    integer(c_int64_t), allocatable, target :: name_ascii(:)
    integer(c_int64_t), allocatable, target :: rxt_ascii(:,:)
    integer(c_int64_t), allocatable, target :: rxt_map(:)
    integer(c_int64_t) :: idx_c

    get_rxt_ndx = -1
    if (rxt_tag_cnt <= 0) then
       call mo_chem_utls_log_direct(get_rxt_ndx_codon_logged, 'get_rxt_ndx direct = codon')
       return
    end if
    allocate(name_ascii(max(1, len(rxt_tag))))
    allocate(rxt_ascii(len(rxt_tag_lst), rxt_tag_cnt))
    allocate(rxt_map(rxt_tag_cnt))
    do ichar = 1,len(rxt_tag)
       name_ascii(ichar) = int(iachar(rxt_tag(ichar:ichar)), c_int64_t)
    end do
    if (len(rxt_tag) == 0) name_ascii(1) = 32_c_int64_t
    do m = 1,rxt_tag_cnt
       do ichar = 1,len(rxt_tag_lst)
          rxt_ascii(ichar,m) = int(iachar(rxt_tag_lst(m)(ichar:ichar)), c_int64_t)
       end do
       rxt_map(m) = int(rxt_tag_map(m), c_int64_t)
    end do
    idx_c = chem_lookup_mapped_name_codon(int(len(rxt_tag), c_int64_t), c_loc(name_ascii(1)), &
         int(len(rxt_tag_lst), c_int64_t), c_loc(rxt_ascii(1,1)), c_loc(rxt_map(1)), &
         int(rxt_tag_cnt, c_int64_t))
    get_rxt_ndx = int(idx_c)
    call mo_chem_utls_log_direct(get_rxt_ndx_codon_logged, 'get_rxt_ndx direct = codon')

  end function get_rxt_ndx

  subroutine mo_chem_utls_log_direct(logged, line)
    implicit none

    logical, intent(inout) :: logged
    character(len=*), intent(in) :: line

    if (.not. logged) then
       if (masterproc) then
          write(iulog,'(A)') trim(line)
          call flush(iulog)
       end if
       logged = .true.
    end if
  end subroutine mo_chem_utls_log_direct

end module mo_chem_utls
