module mo_util

  use shr_kind_mod, only : r8 => shr_kind_r8
  use cam_logfile, only : iulog
  use spmd_utils, only : masterproc

  implicit none

  private
  public :: rebin

  logical :: rebin_use_native_impl = .false.
  logical :: rebin_impl_selected = .false.

contains

  subroutine rebin( nsrc, ntrg, src_x, trg_x, src, trg )
    !---------------------------------------------------------------
    !	... rebin src to trg
    !---------------------------------------------------------------

    use iso_c_binding, only : c_int64_t, c_loc, c_ptr

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer, intent(in)   :: nsrc                  ! dimension source array
    integer, intent(in)   :: ntrg                  ! dimension target array
    real(r8), target, intent(in)  :: src_x(nsrc+1) ! source coordinates
    real(r8), target, intent(in)  :: trg_x(ntrg+1) ! target coordinates
    real(r8), target, intent(in)  :: src(nsrc)     ! source array
    real(r8), target, intent(out) :: trg(ntrg)     ! target array

    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    integer  :: i, l
    integer  :: si, si1
    integer  :: sil, siu
    real(r8)     :: y
    real(r8)     :: sl, su
    real(r8)     :: tl, tu

    interface
       subroutine rebin_codon(nsrc_c, ntrg_c, src_x_p, trg_x_p, src_p, trg_p) bind(c, name="rebin_codon")
         use iso_c_binding, only : c_int64_t, c_ptr
         integer(c_int64_t), value :: nsrc_c, ntrg_c
         type(c_ptr), value :: src_x_p, trg_x_p, src_p, trg_p
       end subroutine rebin_codon
    end interface

    call rebin_select_impl()

    if (.not. rebin_use_native_impl) then
       call rebin_codon( int(nsrc, c_int64_t), int(ntrg, c_int64_t), c_loc(src_x), c_loc(trg_x), c_loc(src), c_loc(trg) )
       return
    end if

    !---------------------------------------------------------------
    !	... check interval overlap
    !---------------------------------------------------------------
    !     if( trg_x(1) < src_x(1) .or. trg_x(ntrg+1) > src_x(nsrc+1) ) then
    !        write(iulog,*) 'rebin: target grid is outside source grid'
    !        write(iulog,*) '       target grid from ',trg_x(1),' to ',trg_x(ntrg+1)
    !        write(iulog,*) '       source grid from ',src_x(1),' to ',src_x(nsrc+1)
    !        call endrun
    !     end if

    do i = 1,ntrg
       tl = trg_x(i)
       if( tl < src_x(nsrc+1) ) then
          do sil = 1,nsrc+1
             if( tl <= src_x(sil) ) then
                exit
             end if
          end do
          tu = trg_x(i+1)
          do siu = 1,nsrc+1
             if( tu <= src_x(siu) ) then
                exit
             end if
          end do
          y   = 0._r8
          sil = max( sil,2 )
          siu = min( siu,nsrc+1 )
          do si = sil,siu
             si1 = si - 1
             sl  = max( tl,src_x(si1) )
             su  = min( tu,src_x(si) )
             y   = y + (su - sl)*src(si1)
          end do
          trg(i) = y/(trg_x(i+1) - trg_x(i))
       else
          trg(i) = 0._r8
       end if
    end do

  end subroutine rebin

  subroutine rebin_select_impl()

    implicit none

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (rebin_impl_selected) return

    impl_name = 'codon'
    call get_environment_variable('REBIN_IMPL', value=impl_name, length=n, status=status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       rebin_use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       rebin_use_native_impl = .false.
    end if

    rebin_impl_selected = .true.

    if (masterproc) then
       if (rebin_use_native_impl) then
          write(iulog,*) 'rebin implementation = native'
       else
          write(iulog,*) 'rebin implementation = codon'
       end if
       call flush(iulog)
    end if

  end subroutine rebin_select_impl


end module mo_util
