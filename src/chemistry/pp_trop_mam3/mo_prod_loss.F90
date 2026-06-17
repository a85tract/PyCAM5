



      module mo_prod_loss

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: exp_prod_loss
      public :: imp_prod_loss

      logical :: exp_prod_loss_use_native_impl = .false.
      logical :: exp_prod_loss_impl_selected = .false.
      logical :: exp_prod_loss_impl_logged = .false.
      logical :: imp_prod_loss_use_native_impl = .false.
      logical :: imp_prod_loss_impl_selected = .false.

      contains

      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:,:,:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:,:,:)
      real(r8), intent(in) :: rxt(:,:,:)
      real(r8), intent(in) :: het_rates(:,:,:)

      interface
         subroutine exp_prod_loss_codon() bind(c, name="exp_prod_loss_codon")
         end subroutine exp_prod_loss_codon
      end interface

      call exp_prod_loss_select_impl()

      if (.not. exp_prod_loss_use_native_impl) then
         call exp_prod_loss_codon()
      end if

      end subroutine exp_prod_loss

      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )

      use ppgrid, only : pver
      use iso_c_binding, only : c_loc, c_ptr

      implicit none

!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:), target, intent(out) :: &
            prod, &
            loss
      real(r8), target, intent(in) :: y(:)
      real(r8), target, intent(in) :: rxt(:)
      real(r8), target, intent(in) :: het_rates(:)

      interface
         subroutine imp_prod_loss_codon(prod_p, loss_p, y_p, rxt_p, het_rates_p) bind(c, name="imp_prod_loss_codon")
            use iso_c_binding, only : c_ptr
            type(c_ptr), value :: prod_p, loss_p, y_p, rxt_p, het_rates_p
         end subroutine imp_prod_loss_codon
      end interface

      call imp_prod_loss_select_impl()

      if (.not. imp_prod_loss_use_native_impl) then
         call imp_prod_loss_codon(c_loc(prod), c_loc(loss), c_loc(y), c_loc(rxt), c_loc(het_rates))
         return
      end if



!--------------------------------------------------------------------
! ... loss and production for Implicit method
!--------------------------------------------------------------------


         loss(1) = ( + rxt(1) + rxt(3) + het_rates(1))* y(1)
         prod(1) = 0._r8
         loss(2) = ( + het_rates(2))* y(2)
         prod(2) =rxt(4)*y(3)
         loss(3) = ( + rxt(4) + het_rates(3))* y(3)
         prod(3) = (rxt(5) +.500_r8*rxt(6) +rxt(7))*y(4)
         loss(4) = ( + rxt(5) + rxt(6) + rxt(7) + het_rates(4))* y(4)
         prod(4) = 0._r8
         loss(5) = ( + het_rates(5))* y(5)
         prod(5) = 0._r8
         loss(6) = ( + het_rates(6))* y(6)
         prod(6) = 0._r8
         loss(7) = ( + het_rates(7))* y(7)
         prod(7) = 0._r8
         loss(8) = ( + het_rates(8))* y(8)
         prod(8) = 0._r8
         loss(9) = ( + het_rates(9))* y(9)
         prod(9) = 0._r8
         loss(10) = ( + het_rates(10))* y(10)
         prod(10) = 0._r8
         loss(11) = ( + het_rates(11))* y(11)
         prod(11) = 0._r8
         loss(12) = ( + het_rates(12))* y(12)
         prod(12) = 0._r8
         loss(13) = ( + het_rates(13))* y(13)
         prod(13) = 0._r8
         loss(14) = ( + het_rates(14))* y(14)
         prod(14) = 0._r8
         loss(15) = ( + het_rates(15))* y(15)
         prod(15) = 0._r8
         loss(16) = ( + het_rates(16))* y(16)
         prod(16) = 0._r8
         loss(17) = ( + het_rates(17))* y(17)
         prod(17) = 0._r8
         loss(18) = ( + het_rates(18))* y(18)
         prod(18) = 0._r8
         loss(19) = ( + het_rates(19))* y(19)
         prod(19) = 0._r8
         loss(20) = ( + het_rates(20))* y(20)
         prod(20) = 0._r8

      end subroutine imp_prod_loss

      subroutine imp_prod_loss_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (imp_prod_loss_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('IMP_PROD_LOSS_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         imp_prod_loss_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         imp_prod_loss_use_native_impl = .false.
      end if

      imp_prod_loss_impl_selected = .true.

      if (masterproc) then
         if (imp_prod_loss_use_native_impl) then
            write(iulog,*) 'imp_prod_loss implementation = native'
         else
            write(iulog,*) 'imp_prod_loss implementation = codon'
         end if
         call flush(iulog)
      end if

      end subroutine imp_prod_loss_select_impl

      subroutine exp_prod_loss_select_impl()

      use cam_logfile, only : iulog
      use spmd_utils, only : masterproc

      implicit none

      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (exp_prod_loss_impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('EXP_PROD_LOSS_IMPL', impl_name, n, status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         exp_prod_loss_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         exp_prod_loss_use_native_impl = .false.
      end if

      exp_prod_loss_impl_selected = .true.

      if (masterproc .and. .not. exp_prod_loss_impl_logged) then
         if (exp_prod_loss_use_native_impl) then
            write(iulog,*) 'exp_prod_loss implementation = native'
         else
            write(iulog,*) 'exp_prod_loss implementation = codon'
         end if
         call flush(iulog)
         exp_prod_loss_impl_logged = .true.
      end if

      end subroutine exp_prod_loss_select_impl

      end module mo_prod_loss
