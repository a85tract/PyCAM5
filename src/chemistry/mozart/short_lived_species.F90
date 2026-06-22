!---------------------------------------------------------------------
! Manages the storage of non-transported short-lived chemical species
! in the physics buffer.
!
! Created by: Francis Vitt -- 20 Aug 2008
!---------------------------------------------------------------------
module short_lived_species

  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : slvd_lst, nslvd, gas_pcnst
  use cam_logfile,  only : iulog
  use iso_c_binding, only : c_int64_t, c_loc, c_ptr
  use ppgrid,       only : pcols, pver, begchunk, endchunk
  use spmd_utils,   only : masterproc
  use mo_util,      only : chemistry_misc_codon_touch


  implicit none

  save
  private
  public :: map
  public :: register_short_lived_species
  public :: initialize_short_lived_species
  public :: set_short_lived_species
  public :: get_short_lived_species
  public :: slvd_index
  public :: pbf_idx

  integer :: pbf_idx
  integer :: map(nslvd)

  character(len=16), parameter :: pbufname = 'ShortLivedSpecies'

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine register_short_lived_species
    use physics_buffer, only : pbuf_add_field, dtype_r8

    implicit none

    integer(c_int64_t) :: active_c

    interface
       function register_short_lived_species_codon(active) result(out_c) bind(c, name="register_short_lived_species_codon")
         import :: c_int64_t
         integer(c_int64_t), value :: active
         integer(c_int64_t) :: out_c
       end function register_short_lived_species_codon
    end interface

    active_c = register_short_lived_species_codon(merge(1_c_int64_t, 0_c_int64_t, nslvd >= 1))
    call chemistry_misc_codon_touch('register_short_lived_species', 140)
    if (active_c == 0_c_int64_t) then
       if (masterproc) then
          write(iulog,'(A)') 'register_short_lived_species implementation = codon no-short-lived no-op'
          call flush(iulog)
       end if
       return
    end if
    if (masterproc) then
       write(iulog,'(A)') 'register_short_lived_species implementation = codon; pbuf registration native island'
       call flush(iulog)
    end if

    call pbuf_add_field(pbufname,'global',dtype_r8,(/pcols,pver,nslvd/),pbf_idx)

  end subroutine register_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine initialize_short_lived_species(ncid_ini, pbuf2d)
    use ioFileMod,      only : getfil
    use error_messages, only : handle_ncerr
    use dycore,         only : dycore_is
    use mo_tracname,    only : solsym
    use ncdio_atm,      only : infld
    use pio,            only : file_desc_t
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_get_chunk, pbuf_get_field

    implicit none

    type(file_desc_t), intent(inout) :: ncid_ini
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)


    integer          :: m,n,lchnk
    character(len=8) :: fieldname
    character(len=4) :: dim1name
    logical          :: found
    real(r8),pointer :: tmpptr(:,:,:)   ! temporary pointer
    real(r8),pointer :: tmpptr2(:,:,:)   ! temporary pointer

    call chemistry_misc_codon_touch('initialize_short_lived_species', 163)
    if ( nslvd < 1 ) return

    found = .false.

    if(dycore_is('se')) then  
       dim1name='ncol'
    else
       dim1name='lon'
    end if

    call pbuf_set_field(pbuf2d, pbf_idx, 0._r8)

    allocate(tmpptr(pcols,pver,begchunk:endchunk))

    do m=1,nslvd
       n = map(m)
       fieldname = solsym(n)
       call infld( fieldname,ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
                   tmpptr, found, grid_map='PHYS')

       if (.not.found) then
          tmpptr(:,:,:) = 1.e-36_r8
       endif

       call pbuf_set_field(pbuf2d, pbf_idx, tmpptr, start=(/1,1,m/),kount=(/pcols,pver,1/))
       
       if (masterproc) write(iulog,*)  fieldname, ' is set to short-lived'
  
    enddo

    deallocate(tmpptr)

  end subroutine initialize_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine set_short_lived_species( q, lchnk, ncol, pbuf )

    use physics_buffer, only : physics_buffer_desc, pbuf_set_field

    implicit none 

    real(r8), intent(in)               :: q(pcols,pver,gas_pcnst)
    integer,  intent(in)               :: lchnk, ncol
    type(physics_buffer_desc), pointer :: pbuf(:)

    integer :: m,n

    call chemistry_misc_codon_touch('set_short_lived_species', 164)
    if ( nslvd < 1 ) return

    do m=1,nslvd
       n = map(m)
       call pbuf_set_field(pbuf, pbf_idx, q(:,:,n), start=(/1,1,m/),kount=(/pcols,pver,1/))
    enddo

  end subroutine set_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine get_short_lived_species( q, lchnk, ncol, pbuf )
    use physics_buffer, only : physics_buffer_desc, pbuf_get_field

    implicit none 

    real(r8), intent(inout)            :: q(pcols,pver,gas_pcnst)
    integer,  intent(in)               :: lchnk, ncol
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8),pointer                   :: tmpptr(:,:)


    integer :: m,n 

    if ( nslvd < 1 ) return

    do m=1,nslvd
       n = map(m)
       call pbuf_get_field(pbuf, pbf_idx, tmpptr, start=(/1,1,m/), kount=(/ pcols,pver,1 /))
       q(:ncol,:,n) = tmpptr(:ncol,:)
    enddo

  endsubroutine get_short_lived_species

!---------------------------------------------------------------------
!---------------------------------------------------------------------
  function slvd_index( name )
    implicit none

    interface
       function slvd_index_codon(name_len_c, name_ascii_p, list_len_c, list_ascii_p, list_count_c) result(idx_c) &
            bind(c, name="slvd_index_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: name_len_c, list_len_c, list_count_c
         type(c_ptr), value :: name_ascii_p, list_ascii_p
         integer(c_int64_t) :: idx_c
       end function slvd_index_codon
    end interface

    character(len=*) :: name
    integer :: slvd_index

    integer :: m, ichar, i, code, status, n
    character(len=32) :: impl_name
    integer(c_int64_t), allocatable, target :: name_ascii(:)
    integer(c_int64_t), allocatable, target :: slvd_ascii(:,:)
    logical, save :: logged = .false.

    slvd_index = -1
    impl_name = 'codon'
    call cam_codon_get_impl('SLVD_INDEX_IMPL', impl_name, n, status)
    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
    end if

    if (.not. (status == 0 .and. n > 0 .and. trim(adjustl(impl_name(:n))) == 'native')) then
       allocate(name_ascii(max(1, len(name))))
       allocate(slvd_ascii(len(slvd_lst), max(1, nslvd)))

       name_ascii(:) = 32_c_int64_t
       do ichar = 1, len(name)
          name_ascii(ichar) = int(iachar(name(ichar:ichar)), c_int64_t)
       end do

       slvd_ascii(:,:) = 32_c_int64_t
       do m = 1, nslvd
          do ichar = 1, len(slvd_lst)
             slvd_ascii(ichar,m) = int(iachar(slvd_lst(m)(ichar:ichar)), c_int64_t)
          end do
       end do

       slvd_index = int(slvd_index_codon(int(len(name), c_int64_t), c_loc(name_ascii(1)), &
            int(len(slvd_lst), c_int64_t), c_loc(slvd_ascii(1,1)), int(nslvd, c_int64_t)))
       if (masterproc .and. .not. logged) then
          write(iulog,'(A)') 'slvd_index direct = codon'
          call flush(iulog)
          logged = .true.
       end if
       return
    end if

    if ( nslvd < 1 ) return

    do m=1,nslvd
       if ( name == slvd_lst(m) ) then
          slvd_index = m
          return 
       endif
    enddo

  endfunction slvd_index

end module short_lived_species
