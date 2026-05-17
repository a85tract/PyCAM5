!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_vrtqdr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:15 $
!
      module rrtmg_sw_vrtqdr

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use shr_kind_mod, only: r8 => shr_kind_r8
      use cam_logfile, only: iulog
      use spmd_utils, only: masterproc

!      use parkind, only: jpim, jprb
!      use parrrsw, only: ngptsw

      implicit none
      save

      logical :: use_native_vrtqdr_sw_impl = .false.
      logical :: vrtqdr_sw_impl_selected = .false.
      logical :: vrtqdr_sw_entered_logged = .false.

      contains

! --------------------------------------------------------------------------
      subroutine vrtqdr_sw(klev, kw, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu)
! --------------------------------------------------------------------------
      use iso_c_binding, only: c_int64_t, c_loc, c_ptr

      integer, intent (in) :: klev
      integer, intent (in) :: kw

      real(kind=r8), target, intent(in) :: pref(:)
      real(kind=r8), target, intent(in) :: prefd(:)
      real(kind=r8), target, intent(in) :: ptra(:)
      real(kind=r8), target, intent(in) :: ptrad(:)

      real(kind=r8), target, intent(in) :: pdbt(:)
      real(kind=r8), target, intent(in) :: ptdbt(:)

      real(kind=r8), target, intent(inout) :: prdnd(:)
      real(kind=r8), target, intent(inout) :: prup(:)
      real(kind=r8), target, intent(inout) :: prupd(:)

      real(kind=r8), target, intent(out) :: pfd(:,:)
      real(kind=r8), target, intent(out) :: pfu(:,:)

      real(kind=r8), target :: ztdn(klev+1)

      interface
         subroutine rrtmg_sw_vrtqdr_codon(klev_c, kw_c, pref_p, prefd_p, ptra_p, ptrad_p, pdbt_p, &
              prdnd_p, prup_p, prupd_p, ptdbt_p, pfd_p, pfu_p, ztdn_p) bind(c, name="rrtmg_sw_vrtqdr_codon")
            use iso_c_binding, only: c_int64_t, c_ptr
            integer(c_int64_t), value :: klev_c, kw_c
            type(c_ptr), value :: pref_p, prefd_p, ptra_p, ptrad_p, pdbt_p
            type(c_ptr), value :: prdnd_p, prup_p, prupd_p, ptdbt_p, pfd_p, pfu_p, ztdn_p
         end subroutine rrtmg_sw_vrtqdr_codon
      end interface

      call vrtqdr_sw_select_impl()
      if (use_native_vrtqdr_sw_impl) then
         call vrtqdr_sw_native(klev, kw, pref, prefd, ptra, ptrad, pdbt, prdnd, prup, prupd, ptdbt, pfd, pfu)
      else
         call vrtqdr_sw_log_entered()
         call rrtmg_sw_vrtqdr_codon( &
              int(klev, c_int64_t), int(kw, c_int64_t), &
              c_loc(pref(1)), c_loc(prefd(1)), c_loc(ptra(1)), c_loc(ptrad(1)), c_loc(pdbt(1)), &
              c_loc(prdnd(1)), c_loc(prup(1)), c_loc(prupd(1)), c_loc(ptdbt(1)), &
              c_loc(pfd(1,1)), c_loc(pfu(1,1)), c_loc(ztdn(1)) &
         )
      end if

      end subroutine vrtqdr_sw

! --------------------------------------------------------------------------
      subroutine vrtqdr_sw_native(klev, kw, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu)
! --------------------------------------------------------------------------
 
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input

      integer, intent (in) :: klev                   ! number of model layers
      integer, intent (in) :: kw                     ! g-point index

      real(kind=r8), intent(in) :: pref(:)                    ! direct beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: prefd(:)                   ! diffuse beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptra(:)                    ! direct beam transmissivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptrad(:)                   ! diffuse beam transmissivity
                                                                 !   Dimensions: (nlayers+1)

      real(kind=r8), intent(in) :: pdbt(:)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptdbt(:)
                                                                 !   Dimensions: (nlayers+1)

      real(kind=r8), intent(inout) :: prdnd(:)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prup(:)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prupd(:)
                                                                 !   Dimensions: (nlayers+1)

! Output
      real(kind=r8), intent(out) :: pfd(:,:)                   ! downwelling flux (W/m2)
                                                                 !   Dimensions: (nlayers+1,ngptsw)
                                                                 ! unadjusted for earth/sun distance or zenith angle
      real(kind=r8), intent(out) :: pfu(:,:)                   ! upwelling flux (W/m2)
                                                                 !   Dimensions: (nlayers+1,ngptsw)
                                                                 ! unadjusted for earth/sun distance or zenith angle

! Local

      integer :: ikp, ikx, jk

      real(kind=r8) :: zreflect
      real(kind=r8) :: ztdn(klev+1)  

! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------
                   
! Link lowest layer with surface
             
      zreflect = 1._r8 / (1._r8 - prefd(klev+1) * prefd(klev))
      prup(klev) = pref(klev) + (ptrad(klev) * &
                 ((ptra(klev) - pdbt(klev)) * prefd(klev+1) + &
                   pdbt(klev) * pref(klev+1))) * zreflect
      prupd(klev) = prefd(klev) + ptrad(klev) * ptrad(klev) * &
                    prefd(klev+1) * zreflect

! Pass from bottom to top 

      do jk = 1,klev-1
         ikp = klev+1-jk                       
         ikx = ikp-1
         zreflect = 1._r8 / (1._r8 -prupd(ikp) * prefd(ikx))
         prup(ikx) = pref(ikx) + (ptrad(ikx) * &
                   ((ptra(ikx) - pdbt(ikx)) * prupd(ikp) + &
                     pdbt(ikx) * prup(ikp))) * zreflect
         prupd(ikx) = prefd(ikx) + ptrad(ikx) * ptrad(ikx) * &
                      prupd(ikp) * zreflect
      enddo
    
! Upper boundary conditions

      ztdn(1) = 1._r8
      prdnd(1) = 0._r8
      ztdn(2) = ptra(1)
      prdnd(2) = prefd(1)

! Pass from top to bottom

      do jk = 2,klev
         ikp = jk+1
         zreflect = 1._r8 / (1._r8 - prefd(jk) * prdnd(jk))
         ztdn(ikp) = ptdbt(jk) * ptra(jk) + &
                    (ptrad(jk) * ((ztdn(jk) - ptdbt(jk)) + &
                     ptdbt(jk) * pref(jk) * prdnd(jk))) * zreflect
         prdnd(ikp) = prefd(jk) + ptrad(jk) * ptrad(jk) * &
                      prdnd(jk) * zreflect
      enddo
    
! Up and down-welling fluxes at levels

      do jk = 1,klev+1
         zreflect = 1._r8 / (1._r8 - prdnd(jk) * prupd(jk))
         pfu(jk,kw) = (ptdbt(jk) * prup(jk) + &
                      (ztdn(jk) - ptdbt(jk)) * prupd(jk)) * zreflect
         pfd(jk,kw) = ptdbt(jk) + (ztdn(jk) - ptdbt(jk)+ &
                      ptdbt(jk) * prup(jk) * prdnd(jk)) * zreflect
      enddo

      end subroutine vrtqdr_sw_native

! --------------------------------------------------------------------------
      subroutine vrtqdr_sw_select_impl()

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (vrtqdr_sw_impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('RRTMG_SW_VRTQDR_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_vrtqdr_sw_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_vrtqdr_sw_impl = .false.
      end if

      vrtqdr_sw_impl_selected = .true.

      if (masterproc) then
         if (use_native_vrtqdr_sw_impl) then
            write(iulog,*) 'rrtmg_sw_vrtqdr implementation = native'
         else
            write(iulog,*) 'rrtmg_sw_vrtqdr implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine vrtqdr_sw_select_impl

! --------------------------------------------------------------------------
      subroutine vrtqdr_sw_log_entered()

      if (vrtqdr_sw_entered_logged) return
      vrtqdr_sw_entered_logged = .true.

      if (masterproc) then
         write(iulog,*) 'rrtmg_sw_vrtqdr entered (vertical quadrature = codon)'
         call flush(iulog)
      end if

      end subroutine vrtqdr_sw_log_entered

      end module rrtmg_sw_vrtqdr
