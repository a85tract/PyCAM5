
module mo_solarproton

	  use shr_kind_mod,  only: r8 => shr_kind_r8
	  use physconst,      only: pi
	  use mo_util,        only: chemistry_misc_codon_touch
	  use cam_logfile,    only: iulog
	  use spmd_utils,     only: masterproc
	  use iso_c_binding,  only: c_int64_t, c_loc, c_ptr

	  implicit none

  save
  logical :: spe_prod_codon_logged = .false.

contains

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine spe_init

    use spedata, only : spedata_init

    implicit none

    !-----------------------------------------------------------------------
    !      ... read in SPE ionization rates
    !-----------------------------------------------------------------------

	    call chemistry_misc_codon_touch('spe_init', 148)
	    call spedata_init()

  end subroutine spe_init

  !-----------------------------------------------------------------------
  !
  !     ... calculates NO production on calday (output in molec/cm3/s)
  !
  !-----------------------------------------------------------------------
  subroutine spe_prod( noxprod, hoxprod, pmid, zmid, lchnk, ncol)

    use mo_apex, only : alatm              ! magnetic latitude grid (radians)
    use ppgrid,  only : pcols, pver
    use spedata, only : get_ionpairs_profile, spe_run
    use spehox,  only : hox_prod_factor

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in) ::  &
         ncol, &                           ! column count
         lchnk                             ! chunk index
    real(r8), intent(in) :: &
         pmid(pcols,pver)                 ! midpoint pressure (Pa)
    real(r8), intent(in) :: &
         zmid(ncol,pver)                  ! midpoint altitude (km)
    real(r8), intent(out), target :: &
         noxprod(ncol,pver)                ! NO production
    real(r8), intent(out), target :: &
         hoxprod(ncol,pver)               ! HOx production

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------

    integer  :: i
    real(r8) :: dlat_aur
    logical  :: do_spe(ncol)
    integer(c_int64_t) :: active_c

    real(r8) :: ion_pairs(pver)
    real(r8), parameter :: noxprod_factor = 1._r8 
    real(r8) :: hoxprod_factor(pver)

    interface
       function spe_prod_codon(active, ncol_c, pver_c, noxprod_p, hoxprod_p) result(out_c) &
            bind(c, name="spe_prod_codon")
         import :: c_int64_t, c_ptr
         integer(c_int64_t), value :: active, ncol_c, pver_c
         type(c_ptr), value :: noxprod_p, hoxprod_p
         integer(c_int64_t) :: out_c
       end function spe_prod_codon
    end interface

    !-----------------------------------------------------------------------
    ! 	... intialize NO production
    !-----------------------------------------------------------------------

    active_c = spe_prod_codon(merge(1_c_int64_t, 0_c_int64_t, spe_run), &
         int(ncol, c_int64_t), int(pver, c_int64_t), c_loc(noxprod), c_loc(hoxprod))
    if (.not. spe_prod_codon_logged) then
       spe_prod_codon_logged = .true.
       if (masterproc) then
          if (active_c == 0_c_int64_t) then
             write(iulog,'(A)') 'spe_prod implementation = codon flag-off zero no-op'
          else
             write(iulog,'(A)') 'spe_prod implementation = codon; active ion-pair body = native island'
          end if
          call flush(iulog)
       end if
    end if

    if (active_c == 0_c_int64_t) return

    !-----------------------------------------------------------------------
    ! 	... check magnetic latitudes, and return if all below 60 deg
    !-----------------------------------------------------------------------
    do i = 1,ncol
       dlat_aur = alatm(i,lchnk)
       do_spe(i) = abs( dlat_aur ) > pi/3._r8
    enddo

    if( all( .not. do_spe(:) ) ) then
       return
    end if

    do i = 1,ncol
       if( do_spe(i) ) then
          call get_ionpairs_profile( pmid(i,:pver), ion_pairs(:pver) )
          noxprod(i,:pver) = noxprod_factor * ion_pairs(:pver)
          hoxprod_factor(:pver) = hox_prod_factor( ion_pairs(:pver), zmid(i,:pver) )
          hoxprod(i,:pver) = hoxprod_factor(:pver)* ion_pairs(:pver)
       end if
    end do

  end subroutine spe_prod

end module mo_solarproton
