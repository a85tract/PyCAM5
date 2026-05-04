subroutine qneg3 (subnam  ,idx     ,ncol    ,ncold   ,lver    ,lconst_beg  , &
                  lconst_end       ,qmin    ,q       )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check moisture and tracers for minimum value, reset any below
! minimum value to minimum value and return information to allow
! warning message to be printed. The global average is NOT preserved.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_logfile,  only: iulog
   use spmd_utils,   only: masterproc
   use iso_c_binding, only: c_double, c_int64_t, c_loc, c_ptr
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   character*(*), intent(in) :: subnam ! name of calling routine

   integer, intent(in) :: idx          ! chunk/latitude index
   integer, intent(in) :: ncol         ! number of atmospheric columns
   integer, intent(in) :: ncold        ! declared number of atmospheric columns
   integer, intent(in) :: lver         ! number of vertical levels in column
   integer, intent(in) :: lconst_beg   ! beginning constituent
   integer, intent(in) :: lconst_end   ! ending    constituent

   real(r8), target, intent(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

!
! Input/Output arguments
!
   real(r8), target, intent(inout) :: q(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field

   logical, save :: use_native_impl = .false.
   logical, save :: impl_selected = .false.

   integer(c_int64_t), target :: indx(ncol,lver)  ! array of indices of points < qmin
   integer(c_int64_t), target :: nval(lver)       ! number of points < qmin for 1 level
   integer(c_int64_t), target :: nvals(lconst_beg:lconst_end)
   integer(c_int64_t), target :: iw(lconst_beg:lconst_end)
   integer(c_int64_t), target :: kw(lconst_beg:lconst_end)
   integer :: m
   integer :: mloc
   real(c_double), target :: worst(lconst_beg:lconst_end)

   interface
      subroutine qneg_batch_3_codon(ncol_c, ncold_c, lver_c, nconst_c, &
           qmin_p, q_p, indx_p, nval_p, nvals_p, worst_p, iw_p, kw_p) bind(c, name="qneg_batch_3_codon")
         use iso_c_binding, only: c_int64_t, c_ptr
         integer(c_int64_t), value :: ncol_c, ncold_c, lver_c, nconst_c
         type(c_ptr), value :: qmin_p, q_p, indx_p, nval_p, nvals_p, worst_p, iw_p, kw_p
      end subroutine qneg_batch_3_codon
   end interface

   call qneg3_batch_select_impl()

   if (use_native_impl) then
      call qneg3_native(subnam, idx, ncol, ncold, lver, lconst_beg, lconst_end, qmin, q)
      return
   end if

   call qneg3_batch_log_entered()
   call qneg_batch_3_codon( &
        int(ncol, c_int64_t), int(ncold, c_int64_t), int(lver, c_int64_t), &
        int(lconst_end-lconst_beg+1, c_int64_t), &
        c_loc(qmin), c_loc(q), c_loc(indx), c_loc(nval), c_loc(nvals), c_loc(worst), c_loc(iw), c_loc(kw) &
   )

   do m=lconst_beg,lconst_end
      mloc = m - lconst_beg + 1
      if (nvals(m) > 100_c_int64_t .and. abs(real(worst(m), r8)) > max(qmin(m), 1.e-12_r8)) then
         write(iulog,9000)subnam,m,idx,int(nvals(m)),qmin(m),real(worst(m), r8),int(iw(m)),int(kw(m))
      end if
   end do

   return
9000 format(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i7, &
            ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
            1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)

contains

   subroutine qneg3_batch_append_proof(proof_line)
      character(len=*), intent(in) :: proof_line
      character(len=512) :: proof_file
      integer :: status, n, unitno

      proof_file = ''
      call get_environment_variable('QNEG_BATCH_PROOF_FILE', value=proof_file, length=n, status=status)
      if (status == 0 .and. n > 0) then
         open(newunit=unitno, file=trim(proof_file(:n)), status='unknown', position='append', action='write')
         write(unitno,'(A)') trim(proof_line)
         close(unitno)
      end if
   end subroutine qneg3_batch_append_proof

   subroutine qneg3_batch_select_impl()
      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (impl_selected) return

      impl_name = 'codon'
      call get_environment_variable('QNEG_BATCH_IMPL', value=impl_name, length=n, status=status)

      if (status == 0 .and. n > 0) then
         do i = 1, n
            code = iachar(impl_name(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
               impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
         end do
         use_native_impl = trim(adjustl(impl_name(:n))) == 'native'
      else
         use_native_impl = .false.
      end if

      impl_selected = .true.

      if (masterproc) then
         if (use_native_impl) then
            write(iulog,*) 'qneg_batch qneg3 implementation = native'
            call qneg3_batch_append_proof('qneg_batch selector entered implementation = native (qneg3)')
         else
            write(iulog,*) 'qneg_batch qneg3 implementation = codon'
            call qneg3_batch_append_proof('qneg_batch selector entered implementation = codon (qneg3)')
         end if
         call flush(iulog)
      end if
   end subroutine qneg3_batch_select_impl

   subroutine qneg3_batch_log_entered()
      logical, save :: entered_logged = .false.

      if (entered_logged) return
      entered_logged = .true.

      if (masterproc) then
         write(iulog,'(A)') 'qneg_batch entered (qneg3 direct = codon)'
         call qneg3_batch_append_proof('qneg_batch entered (qneg3 direct = codon)')
         call flush(iulog)
      end if
   end subroutine qneg3_batch_log_entered

end subroutine qneg3

subroutine qneg3_native (subnam  ,idx     ,ncol    ,ncold   ,lver    ,lconst_beg  , &
                         lconst_end       ,qmin    ,q       )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check moisture and tracers for minimum value, reset any below
! minimum value to minimum value and return information to allow
! warning message to be printed. The global average is NOT preserved.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Rosinski
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use cam_logfile,  only: iulog
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   character*(*), intent(in) :: subnam ! name of calling routine

   integer, intent(in) :: idx          ! chunk/latitude index
   integer, intent(in) :: ncol         ! number of atmospheric columns
   integer, intent(in) :: ncold        ! declared number of atmospheric columns
   integer, intent(in) :: lver         ! number of vertical levels in column
   integer, intent(in) :: lconst_beg   ! beginning constituent
   integer, intent(in) :: lconst_end   ! ending    constituent

   real(r8), intent(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

!
! Input/Output arguments
!
   real(r8), intent(inout) :: q(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field
!
!---------------------------Local workspace-----------------------------
!
   integer indx(ncol,lver)  ! array of indices of points < qmin
   integer nval(lver)       ! number of points < qmin for 1 level
   integer nvals            ! number of values found < qmin
   integer nn
   integer iwtmp
   integer i,ii,k           ! longitude, level indices
   integer m                ! constituent index
   integer iw,kw            ! i,k indices of worst violator

   logical found            ! true => at least 1 minimum violator found

   real(r8) worst           ! biggest violator
 
!
!-----------------------------------------------------------------------
!

   do m=lconst_beg,lconst_end
      nvals = 0
      found = .false.
      worst = 1.e35_r8
      iw = -1
!
! Test all field values for being less than minimum value. Set q = qmin
! for all such points. Trace offenders and identify worst one.
!
!DIR$ preferstream
      do k=1,lver
         nval(k) = 0
!DIR$ prefervector
         nn = 0
         do i=1,ncol
            if (q(i,k,m) < qmin(m)) then
               nn = nn + 1
               indx(nn,k) = i
            end if
         end do
         nval(k) = nn
      end do

      do k=1,lver
         if (nval(k) > 0) then
            found = .true.
            nvals = nvals + nval(k)
            iwtmp = -1
!cdir nodep,altcode=loopcnt
            do ii=1,nval(k)
               i = indx(ii,k)
               if (q(i,k,m) < worst) then
                  worst = q(i,k,m)
                  iwtmp = ii
               end if
            end do
            if (iwtmp /= -1 ) kw = k
            if (iwtmp /= -1 ) iw = indx(iwtmp,k)
!cdir nodep,altcode=loopcnt
            do ii=1,nval(k)
               i = indx(ii,k)
               q(i,k,m) = qmin(m)
            end do
         end if
      end do

      if (found .and. nvals>100 .and.  abs(worst)>max(qmin(m),1.e-12_r8)) then
         write(iulog,9000)subnam,m,idx,nvals,qmin(m),worst,iw,kw
      end if
   end do
!
   return
9000 format(' QNEG3 from ',a,':m=',i3,' lat/lchnk=',i7, &
            ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
            1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)
end subroutine qneg3_native
