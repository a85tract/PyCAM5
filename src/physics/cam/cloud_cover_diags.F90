!===============================================================================
! cloud cover output
!===============================================================================
module cloud_cover_diags

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver,pverp
  use cam_history,   only: addfld, add_default, phys_decomp, outfld
  use phys_control,  only: phys_getopts
  use cam_logfile,   only: iulog
  use spmd_utils,    only: masterproc
  use iso_c_binding, only: c_int64_t, c_loc

  implicit none

  private

  public :: cloud_cover_diags_init
  public :: cloud_cover_diags_out

  logical :: use_native_cloud_cover_diags_impl = .false.
  logical :: cloud_cover_diags_impl_selected = .false.
  logical :: cloud_cover_diags_proof_written = .false.
  logical :: cloud_cover_diags_init_logged = .false.

  interface
     function cloud_cover_diags_init_codon(flag_c) result(out_c) bind(c, name="cloud_cover_diags_init_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: out_c
     end function cloud_cover_diags_init_codon
     function cloud_cover_diags_out_codon(flag_c) result(out_c) bind(c, name="cloud_cover_diags_out_codon")
       use iso_c_binding, only: c_int64_t
       integer(c_int64_t), value :: flag_c
       integer(c_int64_t) :: out_c
     end function cloud_cover_diags_out_codon
     subroutine cloud_cover_cldsav_codon(ncol_c, pcols_c, pver_c, pverp_c, cld_p, pmid_p, &
          pmxrgn_p, nmxrgn_p, cldtot_p, cldlow_p, cldmed_p, cldhgh_p, irgn_p, &
          clrsky_p, clrskymax_p) bind(c, name="cloud_cover_cldsav_codon")
       use iso_c_binding, only: c_int64_t, c_ptr
       integer(c_int64_t), value :: ncol_c, pcols_c, pver_c, pverp_c
       type(c_ptr), value :: cld_p, pmid_p, pmxrgn_p, nmxrgn_p
       type(c_ptr), value :: cldtot_p, cldlow_p, cldmed_p, cldhgh_p
       type(c_ptr), value :: irgn_p, clrsky_p, clrskymax_p
     end subroutine cloud_cover_cldsav_codon
  end interface

contains

!===============================================================================
subroutine cloud_cover_diags_select_impl()
  character(len=32) :: impl_name
  integer :: status, n, i, code

  if (cloud_cover_diags_impl_selected) return

  impl_name = 'codon'
  call get_environment_variable('CLOUD_COVER_DIAGS_IMPL', value=impl_name, length=n, status=status)

  if (status == 0 .and. n > 0) then
     do i = 1, n
        code = iachar(impl_name(i:i))
        if (code >= iachar('A') .and. code <= iachar('Z')) then
           impl_name(i:i) = achar(code + iachar('a') - iachar('A'))
        end if
     end do
     use_native_cloud_cover_diags_impl = trim(adjustl(impl_name(:n))) == 'native'
  else
     use_native_cloud_cover_diags_impl = .false.
  end if

  cloud_cover_diags_impl_selected = .true.

  if (masterproc) then
     if (use_native_cloud_cover_diags_impl) then
        write(iulog,*) 'cloud_cover_diags implementation = native'
     else
        write(iulog,*) 'cloud_cover_diags implementation = codon'
     end if
  end if
end subroutine cloud_cover_diags_select_impl

!===============================================================================
subroutine cloud_cover_diags_proof_once()
  if (cloud_cover_diags_proof_written) return
  cloud_cover_diags_proof_written = .true.
  if (masterproc) then
     write(iulog,'(A)') 'cloud_cover_diags entered (cldsav max-random overlap helper = codon)'
  end if
end subroutine cloud_cover_diags_proof_once

!===============================================================================
subroutine cloud_cover_diags_log_direct(logged, proof_line)
  logical, intent(inout) :: logged
  character(len=*), intent(in) :: proof_line

  if (logged) return
  logged = .true.

  if (masterproc) then
     write(iulog,'(A)') trim(proof_line)
     call flush(iulog)
  end if
end subroutine cloud_cover_diags_log_direct

!===============================================================================
!===============================================================================
subroutine cloud_cover_diags_init(sampling_seq)

  character(len=*), intent(in) :: sampling_seq
  logical :: history_amwg         ! output the variables used by the AMWG diag package
  integer(c_int64_t) :: active_c

  active_c = cloud_cover_diags_init_codon(1_c_int64_t)
  if (active_c == 0_c_int64_t) return
  call cloud_cover_diags_log_direct(cloud_cover_diags_init_logged, &
       'cloud_cover_diags_init direct = codon; history registration native CAM API island')

  call addfld ('CLOUD   ','fraction',pver, 'A','Cloud fraction'                        ,phys_decomp, sampling_seq=sampling_seq)
  call addfld ('CLDTOT  ','fraction',1,    'A','Vertically-integrated total cloud'     ,phys_decomp, sampling_seq=sampling_seq)
  call addfld ('CLDLOW  ','fraction',1,    'A','Vertically-integrated low cloud'       ,phys_decomp, sampling_seq=sampling_seq)
  call addfld ('CLDMED  ','fraction',1,    'A','Vertically-integrated mid-level cloud' ,phys_decomp, sampling_seq=sampling_seq)
  call addfld ('CLDHGH  ','fraction',1,    'A','Vertically-integrated high cloud'      ,phys_decomp, sampling_seq=sampling_seq)

  ! determine the add_default fields
  call phys_getopts(history_amwg_out           = history_amwg  )
 
  if (history_amwg) then
      call add_default ('CLOUD   ', 1, ' ')
      call add_default ('CLDTOT  ', 1, ' ')
      call add_default ('CLDLOW  ', 1, ' ')
      call add_default ('CLDMED  ', 1, ' ')
      call add_default ('CLDHGH  ', 1, ' ')
  endif
    

end subroutine cloud_cover_diags_init

!===============================================================================
!===============================================================================
subroutine cloud_cover_diags_out(lchnk, ncol, cld, pmid, nmxrgn, pmxrgn )

  integer,  intent(in) :: lchnk, ncol
  real(r8), target, intent(in) :: cld(pcols,pver)
  real(r8), target, intent(in) :: pmid(pcols,pver)
  integer,  target, intent(in) :: nmxrgn(pcols)
  real(r8), target, intent(in) :: pmxrgn(pcols,pverp)

  real(r8), target :: cltot(pcols)            ! Diagnostic total cloud cover
  real(r8), target :: cllow(pcols)            !       "     low  cloud cover
  real(r8), target :: clmed(pcols)            !       "     mid  cloud cover
  real(r8), target :: clhgh(pcols)            !       "     hgh  cloud cover
  integer, target :: irgn(pcols)
  real(r8), target :: clrsky(pcols)
  real(r8), target :: clrskymax(pcols)
  integer(c_int64_t) :: active_c

  active_c = cloud_cover_diags_out_codon(1_c_int64_t)
  if (active_c == 0_c_int64_t) return

  call cloud_cover_diags_select_impl()
  if (use_native_cloud_cover_diags_impl) then
     call cldsav (lchnk, ncol, cld, pmid, cltot, cllow, clmed, clhgh, nmxrgn, pmxrgn)
  else
     call cloud_cover_diags_proof_once()
     call cloud_cover_cldsav_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
          int(pver, c_int64_t), int(pverp, c_int64_t), c_loc(cld(1,1)), c_loc(pmid(1,1)), &
          c_loc(pmxrgn(1,1)), c_loc(nmxrgn(1)), c_loc(cltot(1)), c_loc(cllow(1)), &
          c_loc(clmed(1)), c_loc(clhgh(1)), c_loc(irgn(1)), c_loc(clrsky(1)), &
          c_loc(clrskymax(1)))
     if (masterproc) then
        write(iulog,'(A)') 'cloud_cover_diags_out direct = codon; history outfld native'
     end if
  end if

  !
  ! Dump cloud field information to history tape buffer (diagnostics)
  !
  call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
  call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
  call outfld('CLDMED  ',clmed  ,pcols,lchnk)
  call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)

  call outfld('CLOUD   ',cld    ,pcols,lchnk) 

end subroutine cloud_cover_diags_out

!===============================================================================
!===============================================================================
subroutine cldsav(lchnk   ,ncol    , &
                  cld     ,pmid    ,cldtot  ,cldlow  ,cldmed  , &
                  cldhgh  ,nmxrgn  ,pmxrgn  )
!-----------------------------------------------------------------------
   use iso_c_binding, only: c_int64_t, c_loc
   integer, intent(in) :: lchnk
   integer, intent(in) :: ncol
   real(r8), target, intent(in) :: cld(pcols,pver)
   real(r8), target, intent(in) :: pmid(pcols,pver)
   real(r8), target, intent(in) :: pmxrgn(pcols,pverp)
   integer, target, intent(in) :: nmxrgn(pcols)
   real(r8), target, intent(out) :: cldtot(pcols)
   real(r8), target, intent(out) :: cldlow(pcols)
   real(r8), target, intent(out) :: cldmed(pcols)
   real(r8), target, intent(out) :: cldhgh(pcols)

   integer, target :: irgn(pcols)
   real(r8), target :: clrsky(pcols)
   real(r8), target :: clrskymax(pcols)

   call cloud_cover_diags_select_impl()
   if (use_native_cloud_cover_diags_impl) then
      call cldsav_native(lchnk, ncol, cld, pmid, cldtot, cldlow, cldmed, cldhgh, nmxrgn, pmxrgn)
   else
      call cloud_cover_diags_proof_once()
      call cloud_cover_cldsav_codon(int(ncol, c_int64_t), int(pcols, c_int64_t), &
           int(pver, c_int64_t), int(pverp, c_int64_t), c_loc(cld(1,1)), c_loc(pmid(1,1)), &
           c_loc(pmxrgn(1,1)), c_loc(nmxrgn(1)), c_loc(cldtot(1)), c_loc(cldlow(1)), &
           c_loc(cldmed(1)), c_loc(cldhgh(1)), c_loc(irgn(1)), c_loc(clrsky(1)), &
           c_loc(clrskymax(1)))
   end if

   return
end subroutine cldsav

!===============================================================================
subroutine cldsav_native(lchnk   ,ncol    , &
                  cld     ,pmid    ,cldtot  ,cldlow  ,cldmed  , &
                  cldhgh  ,nmxrgn  ,pmxrgn  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute total & 3 levels of cloud fraction assuming maximum-random overlap.
! Pressure ranges for the 3 cloud levels are specified.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------

   implicit none
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: cld(pcols,pver)     ! Cloud fraction
   real(r8), intent(in) :: pmid(pcols,pver)    ! Level pressures
   real(r8), intent(in) :: pmxrgn(pcols,pverp) ! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc

   integer, intent(in) :: nmxrgn(pcols)        ! Number of maximally overlapped regions
!
! Output arguments
!
   real(r8), intent(out) :: cldtot(pcols)       ! Total random overlap cloud cover
   real(r8), intent(out) :: cldlow(pcols)       ! Low random overlap cloud cover
   real(r8), intent(out) :: cldmed(pcols)       ! Middle random overlap cloud cover
   real(r8), intent(out) :: cldhgh(pcols)       ! High random overlap cloud cover

!
!---------------------------Local workspace-----------------------------
!
   integer i,k                  ! Longitude,level indices
   integer irgn(pcols)          ! Max-overlap region index
   integer max_nmxrgn           ! maximum value of nmxrgn over columns
   integer ityp                 ! Type counter
   real(r8) clrsky(pcols)       ! Max-random clear sky fraction
   real(r8) clrskymax(pcols)    ! Maximum overlap clear sky fraction
!------------------------------Parameters-------------------------------
   real(r8) plowmax             ! Max prs for low cloud cover range
   real(r8) plowmin             ! Min prs for low cloud cover range
   real(r8) pmedmax             ! Max prs for mid cloud cover range
   real(r8) pmedmin             ! Min prs for mid cloud cover range
   real(r8) phghmax             ! Max prs for hgh cloud cover range
   real(r8) phghmin             ! Min prs for hgh cloud cover range
!
   parameter (plowmax = 120000._r8,plowmin = 70000._r8, &
              pmedmax =  70000._r8,pmedmin = 40000._r8, &
              phghmax =  40000._r8,phghmin =  5000._r8)

   real(r8) ptypmin(4)
   real(r8) ptypmax(4)

   data ptypmin /phghmin, plowmin, pmedmin, phghmin/
   data ptypmax /plowmax, plowmax, pmedmax, phghmax/
!
!-----------------------------------------------------------------------
!
! Initialize region number
!
   max_nmxrgn = -1
   do i=1,ncol
      max_nmxrgn = max(max_nmxrgn,nmxrgn(i))
   end do

   do ityp = 1, 4
      irgn(1:ncol) = 1
      do k =1,max_nmxrgn-1
         do i=1,ncol
            if (pmxrgn(i,irgn(i)) < ptypmin(ityp) .and. irgn(i) < nmxrgn(i)) then
               irgn(i) = irgn(i) + 1
            end if
         end do
      end do
!
! Compute cloud amount by estimating clear-sky amounts
!
      clrsky(1:ncol)    = 1.0_r8
      clrskymax(1:ncol) = 1.0_r8
      do k = 1, pver
         do i=1,ncol
            if (pmid(i,k) >= ptypmin(ityp) .and. pmid(i,k) <= ptypmax(ityp)) then
               if (pmxrgn(i,irgn(i)) < pmid(i,k) .and. irgn(i) < nmxrgn(i)) then
                  irgn(i) = irgn(i) + 1
                  clrsky(i) = clrsky(i) * clrskymax(i)
                  clrskymax(i) = 1.0_r8
               endif
               clrskymax(i) = min(clrskymax(i),1.0_r8-cld(i,k))
            endif
         end do
      end do
      if (ityp == 1) cldtot(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
      if (ityp == 2) cldlow(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
      if (ityp == 3) cldmed(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
      if (ityp == 4) cldhgh(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
   end do

   return
end subroutine cldsav_native

end module cloud_cover_diags
