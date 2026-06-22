      module euvac

	      use shr_kind_mod,     only : r8 => shr_kind_r8
	      use cam_abortutils,   only : endrun
	      use cam_logfile,      only : iulog
	      use mo_util,          only : chemistry_misc_codon_touch
	      implicit none

      private
      public :: euvac_init
      public :: euvac_set_etf
      public :: euvac_etf

      save

      integer               :: nstruct
      integer               :: nbins
      real(r8), allocatable :: wc(:)                ! wave interval center (nm)
      real(r8), allocatable :: we(:)                ! wave interval edges (nm)
      real(r8), allocatable :: wlintv(:)            ! wave interval (nm)
      real(r8), allocatable :: wlintvi(:)           ! inverse wave interval (nm)
      real(r8), allocatable :: refmin(:)
      real(r8), allocatable :: afac(:)
      real(r8), protected, allocatable :: euvac_etf(:)

      logical :: euvac_on 

      contains

      subroutine euvac_init (euvac_file)
!---------------------------------------------------------------
!	... initialize euvac etf module
!---------------------------------------------------------------

      use cam_pio_utils,  only : cam_pio_openfile
      use pio,            only : pio_nowrite, pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, &
                                 pio_get_var, file_desc_t, pio_closefile
      use spmd_utils,     only : masterproc
      use error_messages, only : alloc_err
      use ioFileMod,      only : getfil
      use iso_c_binding,  only : c_int64_t
      implicit none

      character(len=*), intent(in) :: euvac_file

!---------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------
      type(file_desc_t)  :: ncid
      integer  :: n, ierr
      integer  :: dimid
      integer  :: varid
	      integer  :: astat
	      character(len=256) :: locfn
	      integer(c_int64_t) :: active_c

      interface
         function euvac_init_codon(active) result(out_c) bind(c, name="euvac_init_codon")
            use iso_c_binding, only : c_int64_t
            integer(c_int64_t), value :: active
            integer(c_int64_t) :: out_c
         end function euvac_init_codon
      end interface

	      call chemistry_misc_codon_touch('euvac_init', 151)
	      euvac_on = len_trim(euvac_file)>0
      active_c = euvac_init_codon(merge(1_c_int64_t, 0_c_int64_t, euvac_on))
      if (active_c == 0_c_int64_t) return

!-----------------------------------------------------------------------
!	... readin the etf data
!-----------------------------------------------------------------------
      call getfil( euvac_file, locfn, 0 )
      call cam_pio_openfile (ncid, trim(locfn), PIO_NOWRITE)
!-----------------------------------------------------------------------
!	... check primary dimension consistency
!-----------------------------------------------------------------------
      ierr = pio_inq_dimid( ncid, 'dim1_WC', dimid )
      ierr = pio_inq_dimlen( ncid, dimid, nbins )
      ierr = pio_inq_dimid( ncid, 'dim1_WLINT', dimid )
      ierr = pio_inq_dimlen( ncid, dimid, n )
      if( n /= nbins ) then
         write(iulog,*) 'euvac_init: WLINT dimension(',n,') does not match bin count ',nbins
         call endrun
      end if
      ierr = pio_inq_dimid( ncid, 'dim1_REFMIN', dimid )
      ierr = pio_inq_dimlen( ncid, dimid, n )
      if( n /= nbins ) then
         write(iulog,*) 'euvac_init: REFMIN dimension(',n,') does not match bin count ',nbins
         call endrun
      end if
      ierr = pio_inq_dimid( ncid, 'dim1_AFAC', dimid )
      ierr = pio_inq_dimlen( ncid, dimid, n )
      if( n /= nbins ) then
         write(iulog,*) 'euvac_init: AFAC dimension(',n,') does not match bin count ',nbins
         call endrun
      end if

!-----------------------------------------------------------------------
!	... allocate primary arrays
!-----------------------------------------------------------------------
      allocate( wc(nbins), we(nbins+1), wlintv(nbins), wlintvi(nbins), &
           refmin(nbins), afac(nbins), euvac_etf(nbins), stat=astat )
      if( astat /= 0 ) then
         call alloc_err( astat, 'euvac_init', 'wc ... euvac_etf', nbins )
      end if
!-----------------------------------------------------------------------
!	... read primary arrays
!-----------------------------------------------------------------------
      ierr = pio_inq_varid( ncid, 'WC', varid )
      ierr = pio_get_var( ncid, varid, wc )
      ierr = pio_inq_varid( ncid, 'WLINT', varid )
      ierr = pio_get_var( ncid, varid, wlintv )
      ierr = pio_inq_varid( ncid, 'REFMIN', varid )
      ierr = pio_get_var( ncid, varid, refmin )
      ierr = pio_inq_varid( ncid, 'AFAC', varid )
      ierr = pio_get_var( ncid, varid, afac )
      
      call pio_closefile( ncid )


      wlintvi(:)   = 1._r8/wlintv(:)
      we(:nbins)   = wc(:nbins) - .5_r8*wlintv(:nbins)
      we(nbins+1)  = wc(nbins) + .5_r8*wlintv(nbins)

      end subroutine euvac_init

      subroutine euvac_set_etf( f107, f107a )
!---------------------------------------------------------------
!	... set euvac etf
!---------------------------------------------------------------

      use iso_c_binding,  only : c_int64_t
      use spmd_utils,     only : masterproc

      implicit none

!---------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------
      real(r8), intent(in) :: f107
      real(r8), intent(in) :: f107a

!---------------------------------------------------------------
!	... local variables
!---------------------------------------------------------------
      real(r8), parameter :: factor = 80._r8
      integer  :: w
      real(r8) :: pindex
      integer(c_int64_t) :: active_c

      interface
         function euvac_set_etf_codon(active_in) result(active_out) &
              bind(c, name="euvac_set_etf_codon")
            use iso_c_binding, only : c_int64_t
            integer(c_int64_t), value :: active_in
            integer(c_int64_t) :: active_out
         end function euvac_set_etf_codon
      end interface

      active_c = euvac_set_etf_codon(merge(1_c_int64_t, 0_c_int64_t, euvac_on))
      if (active_c /= 0_c_int64_t .and. active_c /= 1_c_int64_t) then
         call endrun('euvac_set_etf_codon: unexpected return value')
      end if
      if (masterproc) then
         write(iulog,*) 'euvac_set_etf implementation = codon'
         call flush(iulog)
      end if

      if (active_c == 0_c_int64_t) return

      pindex = .5_r8*(f107 + f107a) - factor
      euvac_etf(:) = refmin(:) * max( .8_r8,(1._r8 + afac(:)*pindex) )

      if( masterproc ) then
         write(iulog,*) ' '
         write(iulog,*) '--------------------------------------------------------'
         write(iulog,*) 'euvac_set_etf: f107,f107a = ',f107,f107a
#ifdef EUVAC_DIAGS
         write(iulog,*) 'euvac_set_etf:  wc, etf'
         do w = 1,nbins
            write(iulog,'(1p,2g15.7)') wc(w),euvac_etf(w)
         end do
#endif
         write(iulog,*) '--------------------------------------------------------'
         write(iulog,*) ' '
      end if

      end subroutine euvac_set_etf

      end module euvac
