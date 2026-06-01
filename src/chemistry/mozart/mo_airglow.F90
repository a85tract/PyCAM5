
      module mo_airglow

	      use shr_kind_mod,  only : r8 => shr_kind_r8
	      use physconst,     only : avogad
	      use cam_abortutils,    only : endrun
	      use mo_util,       only : chemistry_misc_codon_touch
	      use cam_logfile,   only : iulog
	      use spmd_utils,    only : masterproc
	      use iso_c_binding, only : c_int64_t

	      implicit none

      save

      integer , parameter :: nag      = 3
      real(r8), parameter :: secpday  = 86400._r8
      real(r8), parameter :: daypsec  = 1._r8/secpday
      real(r8), parameter :: hc       = 6.62608e-34_r8*2.9979e8_r8/1.e-9_r8
      real(r8), parameter :: wc_o2_1s = 1._r8/762._r8
      real(r8), parameter :: wc_o2_1d = 1._r8/1270._r8
      real(r8), parameter :: wc_o1d   = 1._r8/630._r8

      integer :: rid_ag1, rid_ag2, rid_ag3
      logical :: has_airglow
      logical :: init_airglow_proof_written = .false.

      private
      public :: airglow, init_airglow

      interface
        function init_airglow_active_codon(active) result(out_c) bind(c, name="init_airglow_active_codon")
          use iso_c_binding, only : c_int64_t
          integer(c_int64_t), value :: active
          integer(c_int64_t) :: out_c
        end function init_airglow_active_codon
      end interface

      contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
        subroutine init_airglow
          use mo_chem_utls, only : get_rxt_ndx
          use cam_history,  only : addfld, phys_decomp
          use ppgrid,       only : pver

	          implicit none
          integer(c_int64_t) :: active_c

	          call chemistry_misc_codon_touch('mo_airglow', 152)
	          rid_ag1 = get_rxt_ndx( 'ag1' )
          rid_ag2 = get_rxt_ndx( 'ag2' )
          rid_ag3 = get_rxt_ndx( 'ag3' )

          has_airglow = rid_ag1 > 0 .and. rid_ag2 > 0 .and. rid_ag3 > 0

          active_c = init_airglow_active_codon(merge(1_c_int64_t, 0_c_int64_t, has_airglow))
          if (.not. init_airglow_proof_written) then
             init_airglow_proof_written = .true.
             if (masterproc) then
                if (active_c == 0_c_int64_t) then
                   write(iulog,'(A)') 'init_airglow direct = codon missing-airglow-reactions no-op'
                else
                   write(iulog,'(A)') 'init_airglow selector = codon; active airglow addfld body = native'
                end if
                call flush(iulog)
             end if
          end if

          if (active_c == 0_c_int64_t) return

          call addfld( 'AIRGLW1',   'K/s ', pver, 'I', 'O2_1D -> O2 + 1.27 micron airglow loss', phys_decomp )
          call addfld( 'AIRGLW2',   'K/s ', pver, 'I', 'O2_1S -> O2 + 762nm airglow loss', phys_decomp )
          call addfld( 'AIRGLW3',   'K/s ', pver, 'I', 'O1D -> O + 630 nm airglow loss', phys_decomp )
          call addfld( 'AIRGLWTOT', 'K/s ', pver, 'I', 'airglow total loss', phys_decomp )

        endsubroutine init_airglow

      subroutine airglow( ag_tot, o2_1s, o2_1d, o1d, rxt, cp, &
                          ncol, lchnk )
!-----------------------------------------------------------------------
!      	... forms the airglow heating rates
!-----------------------------------------------------------------------

      use chem_mods,     only : rxntot
      use ppgrid,        only : pver
      use cam_history,   only : outfld
      use mo_constants,  only : avo => avogadro
      
      implicit none

!-----------------------------------------------------------------------
!     	... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)   ::  ncol                                ! columns in chunck
      integer, intent(in)   ::  lchnk                               ! chunk index
      real(r8), intent(in)  ::  rxt(ncol,pver,rxntot)               ! rxt rates (1/cm^3/s)
      real(r8), intent(in)  ::  o2_1s(ncol,pver)                    ! concentration (mol/mol)
      real(r8), intent(in)  ::  o2_1d(ncol,pver)                    ! concentration (mol/mol)
      real(r8), intent(in)  ::  o1d(ncol,pver)                      ! concentration (mol/mol)
      real(r8), intent(in)  ::  cp(ncol,pver)                       ! specific heat capacity
      real(r8), intent(out) ::  ag_tot(ncol,pver)                   ! airglow total heating rate (K/s)

!-----------------------------------------------------------------------
!     	... local variables
!-----------------------------------------------------------------------
      integer  ::  k
      real(r8) ::  tmp(ncol)
      real(r8) ::  ag_rate(ncol,pver,nag)

      if (.not. has_airglow) return

      do k = 1,pver
         tmp(:)          = hc * avo / cp(:,k)
         ag_rate(:,k,1)  = tmp(:)*rxt(:,k,rid_ag1)*o2_1d(:,k)*wc_o2_1d
         ag_rate(:,k,2)  = tmp(:)*rxt(:,k,rid_ag2)*o2_1s(:,k)*wc_o2_1s
         ag_rate(:,k,3)  = tmp(:)*rxt(:,k,rid_ag3)*o1d(:,k)*wc_o1d
         ag_tot(:,k)     = ag_rate(:,k,1) + ag_rate(:,k,2) + ag_rate(:,k,3)
      end do

!-----------------------------------------------------------------------
!     	... output the rates
!-----------------------------------------------------------------------
      call outfld( 'AIRGLW1', ag_rate(:,:,1), ncol, lchnk )
      call outfld( 'AIRGLW2', ag_rate(:,:,2), ncol, lchnk )
      call outfld( 'AIRGLW3', ag_rate(:,:,3), ncol, lchnk )
      call outfld( 'AIRGLWTOT', ag_tot, ncol, lchnk )

      end subroutine airglow

      end module mo_airglow
