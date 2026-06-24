!-----------------------------------------------------------------------
!
! Manages the adjustment of ClOy and BrOy family components in response
! to conservation issues resulting from advection.
!
! Created by: Francis Vitt
! Date: 21 May 2008
! Modified by Stacy Walters
! Date: 13 August 2008
!-----------------------------------------------------------------------

module clybry_fam

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use ppgrid,        only : pcols, pver
  use chem_mods,     only : gas_pcnst, adv_mass
  use constituents,  only : pcnst
  use short_lived_species,only: set_short_lived_species,get_short_lived_species
  use mo_util,       only : chemistry_misc_codon_touch

  implicit none

  save

  private
  public :: clybry_fam_set
  public :: clybry_fam_adj
  public :: clybry_fam_init

  integer :: id_cly,id_bry

  integer :: id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_clono2
  integer :: id_br,id_bro,id_hbr,id_brono2,id_brcl,id_hobr

  logical :: has_clybry

contains

  !------------------------------------------
  !------------------------------------------
  subroutine clybry_fam_init

    use mo_chem_utls, only : get_spc_ndx
    use iso_c_binding, only : c_int64_t, c_loc
    implicit none

    integer(c_int64_t), target :: lookup_ids(16)
    integer(c_int64_t), target :: ids_c(16)
    integer(c_int64_t), target :: has_clybry_c

    interface
       subroutine clybry_fam_init_codon(lookup_ids_p, ids_p, has_clybry_p) &
            bind(c, name="clybry_fam_init_codon")
          use iso_c_binding, only : c_ptr
          type(c_ptr), value :: lookup_ids_p, ids_p, has_clybry_p
       end subroutine clybry_fam_init_codon
    end interface

    call chemistry_misc_codon_touch('clybry_fam_init', 132)
    lookup_ids(1) = int(get_spc_ndx('CLY'), c_int64_t)
    lookup_ids(2) = int(get_spc_ndx('BRY'), c_int64_t)

    lookup_ids(3) = int(get_spc_ndx('CL'), c_int64_t)
    lookup_ids(4) = int(get_spc_ndx('CLO'), c_int64_t)
    lookup_ids(5) = int(get_spc_ndx('HOCL'), c_int64_t)
    lookup_ids(6) = int(get_spc_ndx('CL2'), c_int64_t)
    lookup_ids(7) = int(get_spc_ndx('CL2O2'), c_int64_t)
    lookup_ids(8) = int(get_spc_ndx('OCLO'), c_int64_t)
    lookup_ids(9) = int(get_spc_ndx('HCL'), c_int64_t)
    lookup_ids(10) = int(get_spc_ndx('CLONO2'), c_int64_t)

    lookup_ids(11) = int(get_spc_ndx('BR'), c_int64_t)
    lookup_ids(12) = int(get_spc_ndx('BRO'), c_int64_t)
    lookup_ids(13) = int(get_spc_ndx('HBR'), c_int64_t)
    lookup_ids(14) = int(get_spc_ndx('BRONO2'), c_int64_t)
    lookup_ids(15) = int(get_spc_ndx('BRCL'), c_int64_t)
    lookup_ids(16) = int(get_spc_ndx('HOBR'), c_int64_t)
    ids_c(:) = 0_c_int64_t
    has_clybry_c = 0_c_int64_t

    call clybry_fam_init_codon(c_loc(lookup_ids), c_loc(ids_c), c_loc(has_clybry_c))
    call clybry_fam_init_log_codon()

    id_cly = int(ids_c(1))
    id_bry = int(ids_c(2))

    id_cl = int(ids_c(3))
    id_clo = int(ids_c(4))
    id_hocl = int(ids_c(5))
    id_cl2 = int(ids_c(6))
    id_cl2o2 = int(ids_c(7))
    id_oclo = int(ids_c(8))
    id_hcl = int(ids_c(9))
    id_clono2 = int(ids_c(10))

    id_br = int(ids_c(11))
    id_bro = int(ids_c(12))
    id_hbr = int(ids_c(13))
    id_brono2 = int(ids_c(14))
    id_brcl = int(ids_c(15))
    id_hobr = int(ids_c(16))

    has_clybry = has_clybry_c /= 0_c_int64_t

  endsubroutine clybry_fam_init

  subroutine clybry_fam_init_log_codon()

    use cam_logfile, only : iulog
    use spmd_utils, only : masterproc

    implicit none

    if (masterproc) then
       write(iulog,*) 'clybry_fam_init implementation = codon'
       call flush(iulog)
    end if

  end subroutine clybry_fam_init_log_codon

!--------------------------------------------------------------
! set the ClOy and BrOy mass mixing ratios
!  - this is call before advection
!--------------------------------------------------------------
  subroutine clybry_fam_set( ncol, lchnk, map2chm, q, pbuf )

    use time_manager,  only : get_nstep
    use physics_buffer, only : physics_buffer_desc
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog
    use spmd_utils, only : masterproc

    implicit none

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in)    :: ncol, lchnk
    integer,  intent(in)    :: map2chm(pcnst)
    real(r8), intent(inout) :: q(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8) :: wrk(ncol,pver,2)
    real(r8) :: mmr(pcols,pver,gas_pcnst)
    integer  :: n, m, status, i, code
    integer(c_int64_t) :: active_c
    character(len=32) :: impl_name
    logical :: use_native_impl

    interface
       function clybry_fam_set_codon(active) result(out_c) bind(c, name="clybry_fam_set_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t), value :: active
         integer(c_int64_t) :: out_c
       end function clybry_fam_set_codon
    end interface

    impl_name = 'codon'
    call cam_codon_get_impl('CLYBRY_FAM_SET_IMPL', impl_name, n, status)
    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
       end do
       use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_impl = .false.
    end if

    if (.not. use_native_impl) then
       active_c = clybry_fam_set_codon(merge(1_c_int64_t, 0_c_int64_t, has_clybry))
       if (active_c == 0_c_int64_t) then
          if (masterproc) then
             write(iulog,'(A)') 'clybry_fam_set direct = codon no-clybry no-op'
             call flush(iulog)
          end if
          return
       end if
    end if

    if (.not. has_clybry) return

    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          mmr(:ncol,:,m) = q(:ncol,:, n)
       endif
    enddo
    call get_short_lived_species( mmr, lchnk, ncol, pbuf )

!--------------------------------------------------------------
!       ... form updated chlorine, bromine atom mass mixing ratios
!--------------------------------------------------------------
    wrk(:,:,1) = cloy( mmr, pcols, ncol )
    wrk(:,:,2) = broy( mmr, pcols, ncol )

    mmr(:ncol,:,id_cly) = wrk(:,:,1)
    mmr(:ncol,:,id_bry) = wrk(:,:,2)

    call set_short_lived_species( mmr, lchnk, ncol, pbuf )
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          q(:ncol,:, n) = mmr(:ncol,:,m)
       endif
    enddo

  end subroutine clybry_fam_set

!--------------------------------------------------------------
! adjust the ClOy and BrOy individual family members 
!  - this is call after advection
!--------------------------------------------------------------
  subroutine clybry_fam_adj( ncol, lchnk, map2chm, q, pbuf )

    use time_manager,  only : is_first_step
    use physics_buffer, only : physics_buffer_desc
    use iso_c_binding, only : c_int64_t
    use cam_logfile, only : iulog
    use spmd_utils, only : masterproc

    implicit none

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in)    :: ncol, lchnk
    integer,  intent(in)    :: map2chm(pcnst)
    real(r8), intent(inout) :: q(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: factor(ncol,pver)
    real(r8) :: wrk(ncol,pver)
    real(r8) :: mmr(pcols,pver,gas_pcnst)

    integer  :: n, m, status, i, code
    integer(c_int64_t) :: active_c
    character(len=32) :: impl_name
    logical :: use_native_impl

    interface
       function clybry_fam_adj_codon(active) result(out_c) bind(c, name="clybry_fam_adj_codon")
         use iso_c_binding, only : c_int64_t
         integer(c_int64_t), value :: active
         integer(c_int64_t) :: out_c
       end function clybry_fam_adj_codon
    end interface

    impl_name = 'codon'
    call cam_codon_get_impl('CLYBRY_FAM_ADJ_IMPL', impl_name, n, status)
    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
       end do
       use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_impl = .false.
    end if

    if (.not. use_native_impl) then
       active_c = clybry_fam_adj_codon(merge(1_c_int64_t, 0_c_int64_t, has_clybry .and. .not. is_first_step()))
       if (active_c == 0_c_int64_t) then
          if (masterproc) then
             write(iulog,'(A)') 'clybry_fam_adj direct = codon no-clybry/no-adjust no-op'
             call flush(iulog)
          end if
          return
       end if
    end if

    if (.not. has_clybry) return

!--------------------------------------------------------------
!       ... CLY,BRY are not adjusted until the end of the first timestep
!--------------------------------------------------------------
    if (is_first_step()) return

    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          mmr(:ncol,:,m) = q(:ncol,:, n)
       endif
    enddo
    call get_short_lived_species( mmr, lchnk, ncol, pbuf )

!--------------------------------------------------------------
!       ... form updated chlorine atom mass mixing ratio
!--------------------------------------------------------------
    wrk(:,:) = cloy( mmr, pcols, ncol )

    factor(:ncol,:) = mmr(:ncol,:,id_cly) / wrk(:ncol,:)
!--------------------------------------------------------------
!       ... adjust "group" members
!--------------------------------------------------------------
    mmr(:ncol,:,id_cl)     = factor(:ncol,:)*mmr(:ncol,:,id_cl)
    mmr(:ncol,:,id_clo)    = factor(:ncol,:)*mmr(:ncol,:,id_clo)
    mmr(:ncol,:,id_hocl)   = factor(:ncol,:)*mmr(:ncol,:,id_hocl)
    mmr(:ncol,:,id_cl2)    = factor(:ncol,:)*mmr(:ncol,:,id_cl2)
    mmr(:ncol,:,id_cl2o2)  = factor(:ncol,:)*mmr(:ncol,:,id_cl2o2)
    mmr(:ncol,:,id_oclo)   = factor(:ncol,:)*mmr(:ncol,:,id_oclo)
    mmr(:ncol,:,id_hcl)    = factor(:ncol,:)*mmr(:ncol,:,id_hcl)
    mmr(:ncol,:,id_clono2) = factor(:ncol,:)*mmr(:ncol,:,id_clono2)

!--------------------------------------------------------------
!        ... form updated bromine atom mass mixing ratio
!--------------------------------------------------------------
    wrk(:,:) = broy( mmr, pcols, ncol )

    factor(:ncol,:) = mmr(:ncol,:,id_bry) / wrk(:ncol,:)
!--------------------------------------------------------------
!       ... adjust "group" members
!--------------------------------------------------------------
    mmr(:ncol,:,id_br)     = factor(:ncol,:)*mmr(:ncol,:,id_br)
    mmr(:ncol,:,id_bro)    = factor(:ncol,:)*mmr(:ncol,:,id_bro)
    mmr(:ncol,:,id_hbr)    = factor(:ncol,:)*mmr(:ncol,:,id_hbr)
    mmr(:ncol,:,id_brono2) = factor(:ncol,:)*mmr(:ncol,:,id_brono2)
    mmr(:ncol,:,id_brcl)   = factor(:ncol,:)*mmr(:ncol,:,id_brcl)
    mmr(:ncol,:,id_hobr)   = factor(:ncol,:)*mmr(:ncol,:,id_hobr)

    call set_short_lived_species( mmr, lchnk, ncol, pbuf )
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          q(:ncol,:, n) = mmr(:ncol,:,m)
       endif
    enddo

  end subroutine clybry_fam_adj

!--------------------------------------------------------------
! private methods
!--------------------------------------------------------------

!--------------------------------------------------------------
! compute the mass mixing retio of ClOy
!--------------------------------------------------------------
  function cloy( q, pcols, ncol )

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in) :: pcols
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)

!--------------------------------------------------------------
!       ... function declaration
!--------------------------------------------------------------
    real(r8) :: cloy(ncol,pver)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: wrk(ncol)
    integer  :: k

    do k = 1,pver
       wrk(:) = q(:ncol,k,id_cl)           /adv_mass(id_cl) &
              + q(:ncol,k,id_clo)          /adv_mass(id_clo) &
              + q(:ncol,k,id_hocl)         /adv_mass(id_hocl) &
              + 2._r8*( q(:ncol,k,id_cl2)  /adv_mass(id_cl2) &
                      + q(:ncol,k,id_cl2o2)/adv_mass(id_cl2o2) ) &
              + q(:ncol,k,id_oclo)         /adv_mass(id_oclo) &
              + q(:ncol,k,id_hcl)          /adv_mass(id_hcl) &
              + q(:ncol,k,id_clono2)       /adv_mass(id_clono2) 
       cloy(:,k) = adv_mass(id_cl) * wrk(:)
    end do

  end function cloy

!--------------------------------------------------------------
! compute the mass mixing retio of BrOy
!--------------------------------------------------------------
  function broy( q, pcols, ncol )

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in) :: pcols
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)

!--------------------------------------------------------------
!       ... function declaration
!--------------------------------------------------------------
    real(r8) :: broy(ncol,pver)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: wrk(ncol)
    integer  :: k

    do k = 1,pver
       wrk(:) = q(:ncol,k,id_br)    /adv_mass(id_br) &
              + q(:ncol,k,id_bro)   /adv_mass(id_bro) &
              + q(:ncol,k,id_hbr)   /adv_mass(id_hbr) &
              + q(:ncol,k,id_brono2)/adv_mass(id_brono2) &
              + q(:ncol,k,id_brcl)  /adv_mass(id_brcl) &
              + q(:ncol,k,id_hobr)  /adv_mass(id_hobr)
       broy(:,k) = adv_mass(id_br) * wrk(:)
    end do

  end function broy

end module clybry_fam
