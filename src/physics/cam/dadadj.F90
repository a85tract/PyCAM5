subroutine dadadj (lchnk   ,ncol    , &
                   pmid    ,pint    ,pdel    ,t       , &
                   q       )
!-----------------------------------------------------------------------
!
! Purpose:
! GFDL style dry adiabatic adjustment
!
! Method:
! if stratification is unstable, adjustment to the dry adiabatic lapse
! rate is forced subject to the condition that enthalpy is conserved.
!
! Author: CMS Contact J.Hack
!
!-----------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,       only: get_lat_p, get_lon_p
   use physconst,       only: cappa
   use cam_abortutils,  only: endrun
   use cam_control_mod, only: nlvdry
   use cam_logfile,     only: iulog
   use spmd_utils,      only: masterproc
   use iso_c_binding,   only: c_double, c_int64_t, c_loc, c_ptr
   implicit none

   integer, intent(in) :: lchnk               ! chunk identifier
   integer, intent(in) :: ncol                ! number of atmospheric columns

   real(r8), target, intent(in)    :: pmid(pcols,pver)   ! pressure at model levels
   real(r8), target, intent(in)    :: pint(pcols,pverp)  ! pressure at model interfaces
   real(r8), target, intent(in)    :: pdel(pcols,pver)   ! vertical delta-p
   real(r8), target, intent(inout) :: t(pcols,pver)      ! temperature (K)
   real(r8), target, intent(inout) :: q(pcols,pver)      ! specific humidity

   logical, save :: use_native_impl = .false.
   logical, save :: impl_selected = .false.

   real(r8), target :: c1dad(pver)
   real(r8), target :: c2dad(pver)
   real(r8), target :: c3dad(pver)
   real(r8), target :: c4dad(pver)
   integer(c_int64_t), target :: dodad(pcols)
   integer(c_int64_t), target :: status_code
   integer(c_int64_t), target :: fail_i
   real(c_double), target :: zeps_fail

   interface
      subroutine dadadj_codon(ncol_c, pcols_c, nlvdry_c, cappa_c, &
           pmid_p, pint_p, pdel_p, t_p, q_p, &
           c1dad_p, c2dad_p, c3dad_p, c4dad_p, dodad_p, &
           status_p, zeps_fail_p, fail_i_p) bind(c, name="dadadj_codon")
        use iso_c_binding, only: c_double, c_int64_t, c_ptr
        integer(c_int64_t), value :: ncol_c, pcols_c, nlvdry_c
        real(c_double), value :: cappa_c
        type(c_ptr), value :: pmid_p, pint_p, pdel_p, t_p, q_p
        type(c_ptr), value :: c1dad_p, c2dad_p, c3dad_p, c4dad_p, dodad_p
        type(c_ptr), value :: status_p, zeps_fail_p, fail_i_p
      end subroutine dadadj_codon
   end interface

   call dadadj_select_impl()

   if (use_native_impl) then
      call dadadj_native(lchnk, ncol, pmid, pint, pdel, t, q)
      return
   end if

   call dadadj_codon( &
        int(ncol, c_int64_t), int(pcols, c_int64_t), int(nlvdry, c_int64_t), real(cappa, c_double), &
        c_loc(pmid), c_loc(pint), c_loc(pdel), c_loc(t), c_loc(q), &
        c_loc(c1dad), c_loc(c2dad), c_loc(c3dad), c_loc(c4dad), c_loc(dodad), &
        c_loc(status_code), c_loc(zeps_fail), c_loc(fail_i) &
   )

   if (status_code /= 0_c_int64_t) then
      write(iulog,*)'DADADJ: No convergence in dry adiabatic adjustment'
      write(iulog,800) get_lat_p(lchnk, int(fail_i)), get_lon_p(lchnk, int(fail_i)), real(zeps_fail, r8)
      call endrun
   end if

   return

800 format(' lat,lon = ',2i5,', zeps= ',e9.4)

contains

   subroutine dadadj_select_impl()
      character(len=32) :: impl_name
      integer :: status, n, i, code

      if (impl_selected) return

      impl_name = 'codon'
      call cam_codon_get_impl('DADADJ_IMPL', impl_name, n, status)

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
            write(iulog,*) 'dadadj implementation = native'
         else
            write(iulog,*) 'dadadj implementation = codon'
         end if
      end if
   end subroutine dadadj_select_impl

end subroutine dadadj

subroutine dadadj_native (lchnk   ,ncol    , &
                          pmid    ,pint    ,pdel    ,t       , &
                          q       )
!-----------------------------------------------------------------------
!
! Purpose:
! GFDL style dry adiabatic adjustment
!
! Method:
! if stratification is unstable, adjustment to the dry adiabatic lapse
! rate is forced subject to the condition that enthalpy is conserved.
!
! Author: CMS Contact J.Hack
!
!-----------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,       only: get_lat_p, get_lon_p
   use physconst,       only: cappa
   use cam_abortutils,  only: endrun
   use cam_control_mod, only: nlvdry
   use cam_logfile,     only: iulog
   implicit none

   integer niter           ! number of iterations for convergence
   parameter (niter = 15)

!
! Arguments
!
   integer, intent(in) :: lchnk               ! chunk identifier
   integer, intent(in) :: ncol                ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)   ! pressure at model levels
   real(r8), intent(in) :: pint(pcols,pverp)  ! pressure at model interfaces
   real(r8), intent(in) :: pdel(pcols,pver)   ! vertical delta-p

!
! Input/output arguments
!
   real(r8), intent(inout) :: t(pcols,pver)      ! temperature (K)
   real(r8), intent(inout) :: q(pcols,pver)      ! specific humidity
!
!---------------------------Local workspace-----------------------------
!
   integer i,k             ! longitude, level indices
   integer jiter           ! iteration index

   real(r8) c1dad(pver)        ! intermediate constant
   real(r8) c2dad(pver)        ! intermediate constant
   real(r8) c3dad(pver)        ! intermediate constant
   real(r8) c4dad(pver)        ! intermediate constant
   real(r8) gammad             ! dry adiabatic lapse rate (deg/Pa)
   real(r8) zeps               ! convergence criterion (deg/Pa)
   real(r8) rdenom             ! reciprocal of denominator of expression
   real(r8) dtdp               ! delta-t/delta-p
   real(r8) zepsdp             ! zeps*delta-p
   real(r8) zgamma             ! intermediate constant
   real(r8) qave               ! mean q between levels

   logical ilconv          ! .TRUE. ==> convergence was attained
   logical dodad(pcols)    ! .TRUE. ==> do dry adjustment
!
!-----------------------------------------------------------------------
!
   zeps = 2.0e-5_r8           ! set convergence criteria
!
! Find gridpoints with unstable stratification
!
   do i=1,ncol
      gammad = cappa*0.5_r8*(t(i,2) + t(i,1))/pint(i,2)
      dtdp = (t(i,2) - t(i,1))/(pmid(i,2) - pmid(i,1))
      dodad(i) = (dtdp + zeps) .gt. gammad
   end do
   do k=2,nlvdry
      do i=1,ncol
         gammad = cappa*0.5_r8*(t(i,k+1) + t(i,k))/pint(i,k+1)
         dtdp = (t(i,k+1) - t(i,k))/(pmid(i,k+1) - pmid(i,k))
         dodad(i) = dodad(i) .or. (dtdp + zeps).gt.gammad
      end do
   end do
!
! Make a dry adiabatic adjustment
! Note: nlvdry ****MUST**** be < pver
!
   do 80 i=1,ncol
      if (dodad(i)) then
         zeps = 2.0e-5_r8
         do k=1,nlvdry
            c1dad(k) = cappa*0.5_r8*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
            c2dad(k) = (1._r8 - c1dad(k))/(1._r8 + c1dad(k))
            rdenom = 1._r8/(pdel(i,k)*c2dad(k) + pdel(i,k+1))
            c3dad(k) = rdenom*pdel(i,k)
            c4dad(k) = rdenom*pdel(i,k+1)
         end do
50       do jiter=1,niter
            ilconv = .true.
            do k=1,nlvdry
               zepsdp = zeps*(pmid(i,k+1) - pmid(i,k))
               zgamma = c1dad(k)*(t(i,k) + t(i,k+1))
               if ((t(i,k+1)-t(i,k)) >= (zgamma+zepsdp)) then
                  ilconv = .false.
                  t(i,k+1) = t(i,k)*c3dad(k) + t(i,k+1)*c4dad(k)
                  t(i,k) = c2dad(k)*t(i,k+1)
                  qave = (pdel(i,k+1)*q(i,k+1) + pdel(i,k)*q(i,k))/(pdel(i,k+1)+ pdel(i,k))
                  q(i,k+1) = qave
                  q(i,k) = qave
               end if
            end do
            if (ilconv) go to 80 ! convergence => next longitude
         end do
!
! Double convergence criterion if no convergence in niter iterations
!
         zeps = zeps + zeps
         if (zeps > 1.e-4_r8) then
            write(iulog,*)'DADADJ: No convergence in dry adiabatic adjustment'
            write(iulog,800) get_lat_p(lchnk,i),get_lon_p(lchnk,i),zeps
            call endrun
         else
            write(iulog,810) zeps,get_lat_p(lchnk,i),get_lon_p(lchnk,i)
            go to 50
         end if
      end if
80    continue
      return
!
! Formats
!
800   format(' lat,lon = ',2i5,', zeps= ',e9.4)
810   format(//,'DADADJ: Convergence criterion doubled to EPS=',E9.4, &
             ' for'/'        DRY CONVECTIVE ADJUSTMENT at Lat,Lon=', &
             2i5)
end subroutine dadadj_native
