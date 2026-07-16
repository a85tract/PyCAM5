module ap_calcram_scheme

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private

  public :: calcram_run

contains

  !> \section arg_table_calcram_run Argument Table
  !! \htmlinclude calcram_run.html
  !!
      subroutine calcram_run(pcols, rair, gravit, ncol,landfrac,icefrac,ocnfrac,obklen,&
           ustar,ram1in,ram1,t,pmid,&
           pdel,fvin,fv)
        !
        ! !DESCRIPTION:
        !
        ! Calc aerodynamic resistance over oceans and sea ice (comes in from land model)
        ! from Seinfeld and Pandis, p.963.
        !
        ! Author: Natalie Mahowald
        !
        implicit none
        integer, intent(in) :: pcols
        real(r8), intent(in) :: rair, gravit
        integer, intent(in) :: ncol
        real(r8),intent(in) :: ram1in(pcols)         !aerodynamical resistance (s/m)
        real(r8),intent(in) :: fvin(pcols)                 ! sfc frc vel from land
        real(r8),intent(out) :: ram1(pcols)         !aerodynamical resistance (s/m)
        real(r8),intent(out) :: fv(pcols)                 ! sfc frc vel from land
        real(r8), intent(in) :: obklen(pcols)                 ! obklen
        real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
        real(r8), intent(in) :: landfrac(pcols)               ! land fraction
        real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
        real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
        real(r8), intent(in) :: t(pcols)       !atm temperature (K)
        real(r8), intent(in) :: pmid(pcols)    !atm pressure (Pa)
        real(r8), intent(in) :: pdel(pcols)    !atm pressure (Pa)
        real(r8), parameter :: zzocen = 0.0001_r8   ! Ocean aerodynamic roughness length
        real(r8), parameter :: zzsice = 0.0400_r8   ! Sea ice aerodynamic roughness length
        real(r8), parameter :: xkar   = 0.4_r8      ! Von Karman constant

        ! local variables
        real(r8) :: z,psi,psi0,nu,nu0,temp,ram
        integer :: i
        !    write(iulog,*) rair,zzsice,zzocen,gravit,xkar


        do i=1,ncol
           z=pdel(i)*rair*t(i)/pmid(i)/gravit/2.0_r8   !use half the layer height like Ganzefeld and Lelieveld, 1995
           if(obklen(i).eq.0) then
              psi=0._r8
              psi0=0._r8
           else
              psi=min(max(z/obklen(i),-1.0_r8),1.0_r8)
              psi0=min(max(zzocen/obklen(i),-1.0_r8),1.0_r8)
           endif
           temp=z/zzocen
           if(icefrac(i) > 0.5_r8) then
              if(obklen(i).gt.0) then
                 psi0=min(max(zzsice/obklen(i),-1.0_r8),1.0_r8)
              else
                 psi0=0.0_r8
              endif
              temp=z/zzsice
	   endif
           if(psi> 0._r8) then
              ram=1/xkar/ustar(i)*(log(temp)+4.7_r8*(psi-psi0))
           else
              nu=(1.00_r8-15.000_r8*psi)**(.25_r8)
              nu0=(1.000_r8-15.000_r8*psi0)**(.25_r8)
              if(ustar(i).ne.0._r8) then
                 ram=1/xkar/ustar(i)*(log(temp) &
                      +log(((nu0**2+1.00_r8)*(nu0+1.0_r8)**2)/((nu**2+1.0_r8)*(nu+1.00_r8)**2)) &
                      +2.0_r8*(atan(nu)-atan(nu0)))
              else
	         ram=0._r8
              endif
           endif
           if(landfrac(i) < 0.000000001_r8) then
              fv(i)=ustar(i)
              ram1(i)=ram
           else
              fv(i)=fvin(i)
              ram1(i)=ram1in(i)
           endif
           !          write(iulog,*) i,pdel(i),t(i),pmid(i),gravit,obklen(i),psi,psi0,icefrac(i),nu,nu0,ram,ustar(i),&
           !             log(((nu0**2+1.00)*(nu0+1.0)**2)/((nu**2+1.0)*(nu+1.00)**2)),2.0*(atan(nu)-atan(nu0))

        enddo

        ! fvitt -- fv == 0 causes a floating point exception in
        ! dry dep of sea salts and dust
        where ( fv(:ncol) == 0._r8 )
           fv(:ncol) = 1.e-12_r8
        endwhere

        return
      end subroutine calcram_run

end module ap_calcram_scheme
