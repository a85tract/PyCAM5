module cmparray_mod

  use shr_kind_mod, only : r8 => shr_kind_r8
  use iso_c_binding, only : c_int64_t, c_loc
  
  implicit none
  private
  save
  
  public expdaynite, cmpdaynite

  logical :: use_native_cmparray_impl = .false.
  logical :: cmparray_impl_selected = .false.
  logical :: cmparray_proof_written = .false.

  interface
     subroutine cmpdaynite_3d_r_copy_codon(in_array_p, out_array_p, idxday_p, idxnite_p, &
          nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c) &
          bind(c, name="cmpdaynite_3d_r_copy_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       type(c_ptr), value :: in_array_p, out_array_p, idxday_p, idxnite_p
       integer(c_int64_t), value :: nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c
     end subroutine cmpdaynite_3d_r_copy_codon
     subroutine cmpdaynite_1d_r_copy_codon(in_array_p, out_array_p, idxday_p, idxnite_p, &
          nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c) &
          bind(c, name="cmpdaynite_1d_r_copy_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       type(c_ptr), value :: in_array_p, out_array_p, idxday_p, idxnite_p
       integer(c_int64_t), value :: nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c
     end subroutine cmpdaynite_1d_r_copy_codon
     subroutine cmpdaynite_2d_r_copy_codon(in_array_p, out_array_p, idxday_p, idxnite_p, &
          nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c) &
          bind(c, name="cmpdaynite_2d_r_copy_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       type(c_ptr), value :: in_array_p, out_array_p, idxday_p, idxnite_p
       integer(c_int64_t), value :: nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c
     end subroutine cmpdaynite_2d_r_copy_codon

     subroutine expdaynite_3d_r_codon(array_p, tmp_p, idxday_p, idxnite_p, &
          nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c) &
          bind(c, name="expdaynite_3d_r_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       type(c_ptr), value :: array_p, tmp_p, idxday_p, idxnite_p
       integer(c_int64_t), value :: nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c
     end subroutine expdaynite_3d_r_codon
     subroutine expdaynite_1d_r_codon(array_p, tmp_p, idxday_p, idxnite_p, &
          nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c) &
          bind(c, name="expdaynite_1d_r_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       type(c_ptr), value :: array_p, tmp_p, idxday_p, idxnite_p
       integer(c_int64_t), value :: nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c
     end subroutine expdaynite_1d_r_codon
     subroutine expdaynite_2d_r_codon(array_p, tmp_p, idxday_p, idxnite_p, &
          nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c) &
          bind(c, name="expdaynite_2d_r_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       type(c_ptr), value :: array_p, tmp_p, idxday_p, idxnite_p
       integer(c_int64_t), value :: nday_c, nnite_c, il1_c, iu1_c, il2_c, iu2_c, il3_c, iu3_c
     end subroutine expdaynite_2d_r_codon
  end interface

  interface CmpDayNite
    module procedure CmpDayNite_1d_R
    module procedure CmpDayNite_2d_R
    module procedure CmpDayNite_3d_R
    module procedure CmpDayNite_1d_R_Copy
    module procedure CmpDayNite_2d_R_Copy
    module procedure CmpDayNite_3d_R_Copy
    module procedure CmpDayNite_1d_I
    module procedure CmpDayNite_2d_I
    module procedure CmpDayNite_3d_I
  end interface ! CmpDayNite

  interface ExpDayNite
    module procedure ExpDayNite_1d_R
    module procedure ExpDayNite_2d_R
    module procedure ExpDayNite_3d_R
    module procedure ExpDayNite_1d_I
    module procedure ExpDayNite_2d_I
    module procedure ExpDayNite_3d_I
  end interface ! ExpDayNite

  interface cmparray
    module procedure cmparray_1d_R
    module procedure cmparray_2d_R
    module procedure cmparray_3d_R
  end interface ! cmparray

  interface chksum
    module procedure chksum_1d_R
    module procedure chksum_2d_R
    module procedure chksum_3d_R
    module procedure chksum_1d_I
    module procedure chksum_2d_I
    module procedure chksum_3d_I
  end interface ! chksum

  contains

  subroutine cmparray_select_impl()

    use cam_logfile, only: iulog
    use spmd_utils, only: masterproc

    character(len=32) :: impl_name
    integer :: status, n, i, code

    if (cmparray_impl_selected) return

    impl_name = 'codon'
    call cam_codon_get_impl('CMPARRAY_IMPL', impl_name, n, status)

    if (status == 0 .and. n > 0) then
       do i = 1, n
          code = iachar(impl_name(i:i))
          if (code >= iachar('A') .and. code <= iachar('Z')) then
             impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
          end if
       end do
       use_native_cmparray_impl = trim(adjustl(impl_name(:n))) == 'native'
    else
       use_native_cmparray_impl = .false.
    end if

    cmparray_impl_selected = .true.

    if (masterproc) then
       if (use_native_cmparray_impl) then
          write(iulog,*) 'cmparray implementation = native'
       else
          write(iulog,*) 'cmparray implementation = codon'
       end if
    end if

  end subroutine cmparray_select_impl

  subroutine cmparray_proof_once()

    use cam_logfile, only: iulog
    use spmd_utils, only: masterproc

    if (cmparray_proof_written) return
    cmparray_proof_written = .true.

    if (masterproc) then
       write(iulog,'(A)') 'cmparray entered (real day-night reorder = codon)'
    end if

  end subroutine cmparray_proof_once

  subroutine CmpDayNite_1d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1) :: Array

    call CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine CmpDayNite_1d_R

  subroutine CmpDayNite_2d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1,il2:iu2) :: Array

    call CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine CmpDayNite_2d_R

  subroutine CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array

    real(r8), dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nnite) = Array(IdxNite(1:Nnite),j,k)
        Array(il1:il1+Nday-1,j,k) = Array(IdxDay(1:Nday),j,k)
        Array(il1+Nday:il1+Nday+Nnite-1,j,k) = tmp(1:Nnite)

      end do
    end do

    return
  end subroutine CmpDayNite_3d_R

  subroutine CmpDayNite_1d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, target, intent(in), dimension(Nday) :: IdxDay
    integer, target, intent(in), dimension(Nnite) :: IdxNite
    real(r8), target, intent(in), dimension(il1:iu1) :: InArray
    real(r8), target, intent(out), dimension(il1:iu1) :: OutArray

    call cmparray_select_impl()

    if (.not. use_native_cmparray_impl) then
       call cmparray_proof_once()
       call cmpdaynite_1d_r_copy_codon(c_loc(InArray), c_loc(OutArray), &
            c_loc(IdxDay), c_loc(IdxNite), int(Nday, c_int64_t), int(Nnite, c_int64_t), &
            int(il1, c_int64_t), int(iu1, c_int64_t), 1_c_int64_t, 1_c_int64_t, &
            1_c_int64_t, 1_c_int64_t)
       return
    end if

    call CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine CmpDayNite_1d_R_Copy

  subroutine CmpDayNite_2d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, target, intent(in), dimension(Nday) :: IdxDay
    integer, target, intent(in), dimension(Nnite) :: IdxNite
    real(r8), target, intent(in), dimension(il1:iu1,il2:iu2) :: InArray
    real(r8), target, intent(out), dimension(il1:iu1,il2:iu2) :: OutArray

    call cmparray_select_impl()

    if (.not. use_native_cmparray_impl) then
       call cmparray_proof_once()
       call cmpdaynite_2d_r_copy_codon(c_loc(InArray), c_loc(OutArray), &
            c_loc(IdxDay), c_loc(IdxNite), int(Nday, c_int64_t), int(Nnite, c_int64_t), &
            int(il1, c_int64_t), int(iu1, c_int64_t), int(il2, c_int64_t), int(iu2, c_int64_t), &
            1_c_int64_t, 1_c_int64_t)
       return
    end if

    call CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine CmpDayNite_2d_R_Copy

  subroutine CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, target, intent(in), dimension(Nday) :: IdxDay
    integer, target, intent(in), dimension(Nnite) :: IdxNite
    real(r8), target, intent(in), dimension(il1:iu1,il2:iu2,il3:iu3) :: InArray
    real(r8), target, intent(out), dimension(il1:iu1,il2:iu2,il3:iu3) :: OutArray

    integer :: i, j, k

    call cmparray_select_impl()

    if (.not. use_native_cmparray_impl) then
       call cmparray_proof_once()
       call cmpdaynite_3d_r_copy_codon(c_loc(InArray), c_loc(OutArray), &
            c_loc(IdxDay), c_loc(IdxNite), int(Nday, c_int64_t), int(Nnite, c_int64_t), &
            int(il1, c_int64_t), int(iu1, c_int64_t), int(il2, c_int64_t), int(iu2, c_int64_t), &
            int(il3, c_int64_t), int(iu3, c_int64_t))
       return
    end if


    do k = il3, iu3
      do j = il2, iu2

         do i=il1,il1+Nday-1
            OutArray(i,j,k) = InArray(IdxDay(i-il1+1),j,k)
         enddo
         do i=il1+Nday,il1+Nday+Nnite-1
            OutArray(i,j,k) = InArray(IdxNite(i-(il1+Nday)+1),j,k)
         enddo
        

      end do
    end do

    return
  end subroutine CmpDayNite_3d_R_Copy

  subroutine CmpDayNite_1d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1) :: Array

    call CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine CmpDayNite_1d_I

  subroutine CmpDayNite_2d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1,il2:iu2) :: Array

    call CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine CmpDayNite_2d_I

  subroutine CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array

    integer, dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nnite) = Array(IdxNite(1:Nnite),j,k)
        Array(il1:il1+Nday-1,j,k) = Array(IdxDay(1:Nday),j,k)
        Array(il1+Nday:il1+Nday+Nnite-1,j,k) = tmp(1:Nnite)

      end do
    end do

    return
  end subroutine CmpDayNite_3d_I

  subroutine ExpDayNite_1d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, target, intent(in), dimension(Nday) :: IdxDay
    integer, target, intent(in), dimension(Nnite) :: IdxNite
    real(r8), target, intent(inout), dimension(il1:iu1) :: Array

    real(r8), target, dimension(il1:iu1) :: tmp

    call cmparray_select_impl()

    if (.not. use_native_cmparray_impl) then
       call cmparray_proof_once()
       call expdaynite_1d_r_codon(c_loc(Array), c_loc(tmp), &
            c_loc(IdxDay), c_loc(IdxNite), int(Nday, c_int64_t), int(Nnite, c_int64_t), &
            int(il1, c_int64_t), int(iu1, c_int64_t), 1_c_int64_t, 1_c_int64_t, &
            1_c_int64_t, 1_c_int64_t)
       return
    end if

    call ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine ExpDayNite_1d_R

  subroutine ExpDayNite_2d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, target, intent(in), dimension(Nday) :: IdxDay
    integer, target, intent(in), dimension(Nnite) :: IdxNite
    real(r8), target, intent(inout), dimension(il1:iu1,il2:iu2) :: Array

    real(r8), target, dimension(il1:iu1) :: tmp

    call cmparray_select_impl()

    if (.not. use_native_cmparray_impl) then
       call cmparray_proof_once()
       call expdaynite_2d_r_codon(c_loc(Array), c_loc(tmp), &
            c_loc(IdxDay), c_loc(IdxNite), int(Nday, c_int64_t), int(Nnite, c_int64_t), &
            int(il1, c_int64_t), int(iu1, c_int64_t), int(il2, c_int64_t), int(iu2, c_int64_t), &
            1_c_int64_t, 1_c_int64_t)
       return
    end if

    call ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine ExpDayNite_2d_R

  subroutine ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, target, intent(in), dimension(Nday) :: IdxDay
    integer, target, intent(in), dimension(Nnite) :: IdxNite
    real(r8), target, intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array

    real(r8), target, dimension(il1:iu1) :: tmp
    integer :: i, j, k

    call cmparray_select_impl()

    if (.not. use_native_cmparray_impl) then
       call cmparray_proof_once()
       call expdaynite_3d_r_codon(c_loc(Array), c_loc(tmp), &
            c_loc(IdxDay), c_loc(IdxNite), int(Nday, c_int64_t), int(Nnite, c_int64_t), &
            int(il1, c_int64_t), int(iu1, c_int64_t), int(il2, c_int64_t), int(iu2, c_int64_t), &
            int(il3, c_int64_t), int(iu3, c_int64_t))
       return
    end if


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nday) = Array(1:Nday,j,k)
        Array(IdxNite(1:Nnite),j,k) = Array(il1+Nday:il1+Nday+Nnite-1,j,k)
        Array(IdxDay(1:Nday),j,k) = tmp(1:Nday)

      end do
    end do

    return
  end subroutine ExpDayNite_3d_R

  subroutine ExpDayNite_1d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1) :: Array

    call ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, 1, 1, 1, 1)

    return
  end subroutine ExpDayNite_1d_I

  subroutine ExpDayNite_2d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1,il2:iu2) :: Array

    call ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2, iu2, 1, 1)

    return
  end subroutine ExpDayNite_2d_I

  subroutine ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1, il2,iu2, il3, iu3)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1,il2:iu2,il3:iu3) :: Array

    integer, dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nday) = Array(1:Nday,j,k)
        Array(IdxNite(1:Nnite),j,k) = Array(il1+Nday:il1+Nday+Nnite-1,j,k)
        Array(IdxDay(1:Nday),j,k) = tmp(1:Nday)

      end do
    end do

    return
  end subroutine ExpDayNite_3d_I

!******************************************************************************!
!                                                                              !
!                                 DEBUG                                        !
!                                                                              !
!******************************************************************************!

  subroutine cmparray_1d_R(name, Ref, New, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    real(r8), intent(in), dimension(id1) :: Ref
    real(r8), intent(in), dimension(id1) :: New

    call cmparray_3d_R(name, Ref, New, id1, is1, ie1, 1, 1, 1, 1, 1, 1)
  end subroutine cmparray_1d_R

  subroutine cmparray_2d_R(name, Ref, New, id1, is1, ie1, id2, is2, ie2)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    real(r8), intent(in), dimension(id1, id2) :: Ref
    real(r8), intent(in), dimension(id1, id2) :: New

    call cmparray_3d_R(name, Ref, New, id1, is1, ie1, id2, is2, ie2, 1, 1, 1)
  end subroutine cmparray_2d_R

  subroutine cmparray_3d_R(name, Ref, New, id1, is1, ie1, id2, is2, ie2, id3, is3, ie3)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in) :: id3, is3, ie3
    real(r8), intent(in), dimension(id1, id2, id3) :: Ref
    real(r8), intent(in), dimension(id1, id2, id3) :: New

    integer :: i, j, k
    integer :: nerr
    logical :: found
    real(r8):: rdiff
    real(r8), parameter :: rtol = 1.0e-13_r8

    nerr = 0

    do k = is3, ie3
      do j = is2, ie2

        found = .false.
        do i = is1, ie1
          rdiff = abs(New(i,j,k)-Ref(i,j,k))
          rdiff = rdiff / merge(abs(Ref(i,j,k)), 1.0_r8, Ref(i,j,k) /= 0.0_r8)
          if ( rdiff > rtol ) then
            found = .true.
            exit
          end if
        end do

        if ( found ) then
          do i = is1, ie1
            rdiff = abs(New(i,j,k)-Ref(i,j,k))
            rdiff = rdiff / merge(abs(Ref(i,j,k)), 1.0_r8, Ref(i,j,k) /= 0.0_r8)
            if ( rdiff > rtol ) then
              print 666, name, i, j, k, Ref(i, j, k), New(i, j, k), rdiff
              nerr = nerr + 1
              if ( nerr > 10 ) stop
            end if
          end do
        end if

      end do
    end do

    return
666 format('cmp3d: ', a10, 3(1x, i4), 3(1x, e20.14))

  end subroutine cmparray_3d_R

  subroutine chksum_1d_R(name, Ref, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    real(r8), intent(in), dimension(id1) :: Ref

    call chksum_3d_R(name, Ref, id1, is1, ie1, 1, 1, 1, 1, 1, 1)
  end subroutine chksum_1d_R

  subroutine chksum_1d_I(name, Ref, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in), dimension(id1) :: Ref

    call chksum_3d_I(name, Ref, id1, is1, ie1, 1, 1, 1, 1, 1, 1)
  end subroutine chksum_1d_I

  subroutine chksum_2d_R(name, Ref, id1, is1, ie1, id2, is2, ie2)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    real(r8), intent(in), dimension(id1, id2) :: Ref

    call chksum_3d_R(name, Ref, id1, is1, ie1, id2, is2, ie2, 1, 1, 1)
  end subroutine chksum_2d_R

  subroutine chksum_2d_I(name, Ref, id1, is1, ie1, id2, is2, ie2)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in), dimension(id1, id2) :: Ref

    call chksum_3d_I(name, Ref, id1, is1, ie1, id2, is2, ie2, 1, 1, 1)
  end subroutine chksum_2d_I

  subroutine chksum_3d_R(name, Ref, id1, is1, ie1, id2, is2, ie2, id3, is3, ie3)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in) :: id3, is3, ie3
!orig    real(r8), intent(in), dimension(id1, id2, id3) :: Ref
    real(r8), intent(in), dimension(is1:ie1, is2:ie2, is3:ie3) :: Ref

    real(r8) :: chksum
    real(r8) :: rmin, rmax
    integer :: i, j, k
    integer :: imin, jmin, kmin
    integer :: imax, jmax, kmax

    imin = is1 ; jmin = is2 ; kmin = is3
    imax = is1 ; jmax = is2 ; kmax = is3
    rmin = Ref(is1, is2, is3) ; rmax = rmin

    chksum = 0.0_r8

    do k = is3, ie3
      do j = is2, ie2
        do i = is1, ie1
          chksum = chksum + abs(Ref(i,j,k))
          if ( Ref(i,j,k) < rmin ) then
            rmin = Ref(i,j,k)
            imin = i ; jmin = j ; kmin = k
          end if
          if ( Ref(i,j,k) > rmax ) then
            rmax = Ref(i,j,k)
            imax = i ; jmax = j ; kmax = k
          end if
        end do
      end do
    end do

    print 666, name, chksum, imin, jmin, kmin, imax, jmax, kmax
666 format('chksum: ', a8, 1x, e20.14, 6(1x, i4))

  end subroutine chksum_3d_R

  subroutine chksum_3d_I(name, Ref, id1, is1, ie1, id2, is2, ie2, id3, is3, ie3)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in) :: id3, is3, ie3
    integer,  intent(in), dimension(id1, id2, id3) :: Ref

    integer :: i, j, k
    integer :: chksum
    chksum = 0

    do k = is3, ie3
      do j = is2, ie2
        do i = is1, ie1
          chksum = chksum + abs(Ref(i,j,k))
        end do
      end do
    end do

    print 666, name, chksum
666 format('chksum: ', a8, 1x, i8)

  end subroutine chksum_3d_I

end module cmparray_mod
