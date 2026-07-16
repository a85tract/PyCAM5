
module mo_exp_sol

  private
  public :: exp_sol
  public :: exp_sol_inti

  integer, save :: exp_sol_scheme_entry_count = 0
!$omp threadprivate(exp_sol_scheme_entry_count)

contains

  subroutine exp_sol_inti

    use mo_tracname, only : solsym
    use chem_mods,   only : clscnt1, clsmap
    use ppgrid,      only : pver
    use cam_history, only : addfld, phys_decomp

    implicit none

    integer :: i,j

    do i = 1,clscnt1

       j = clsmap(i,1)
       call addfld( trim(solsym(j))//'_CHMP', '/cm3/s ', pver, 'I', 'chemical production rate', phys_decomp )
       call addfld( trim(solsym(j))//'_CHML', '/cm3/s ', pver, 'I', 'chemical loss rate',       phys_decomp )

    enddo
  end subroutine exp_sol_inti


  subroutine exp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, xhnm, ncol, lchnk, ltrop )
    !-----------------------------------------------------------------------
    !      	... Exp_sol advances the volumetric mixing ratio
    !           forward one time step via the fully explicit
    !           Euler scheme
    !-----------------------------------------------------------------------

    use chem_mods,     only : clscnt1, extcnt, gas_pcnst, clsmap, rxntot
    use ppgrid,        only : pcols, pver
    use shr_kind_mod,  only : r8 => shr_kind_r8
    use cam_history,   only : outfld
    use mo_tracname,   only : solsym
    use perf_mod,      only : t_startf, t_stopf
    use cam_logfile,   only : iulog
    use cam_abortutils, only : endrun
    use ap_exp_sol_scheme, only : exp_sol_run

    implicit none
    !-----------------------------------------------------------------------
    !     	... Dummy arguments
    !-----------------------------------------------------------------------
    integer,  intent(in)    ::  ncol                                ! columns in chunck
    integer,  intent(in)    ::  lchnk                               ! chunk id
    real(r8), intent(in)    ::  delt                                ! time step (s)
    real(r8), intent(in)    ::  het_rates(ncol,pver,max(1,gas_pcnst))  ! het rates (1/cm^3/s)
    real(r8), intent(in)    ::  reaction_rates(ncol,pver,rxntot)    ! rxt rates (1/cm^3/s)
    real(r8), intent(in)    ::  extfrc(ncol,pver,extcnt)            ! "external insitu forcing" (1/cm^3/s)
    real(r8), intent(in)    ::  xhnm(ncol,pver)
    integer,  intent(in)    ::  ltrop(pcols)                        ! chemistry troposphere boundary (index)
    real(r8), intent(inout) ::  base_sol(ncol,pver,gas_pcnst)       ! working mixing ratios (vmr)

    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    integer  ::  i, k, l, m
    integer  ::  errflg
    character(len=512) :: errmsg
    real(r8), dimension(ncol,pver,clscnt1) :: &
         prod, &
         loss, &
         ind_prd

    real(r8), dimension(ncol,pver) :: wrk

    call t_startf('ap_exp_sol_run')
    call exp_sol_run(clscnt1, exp_sol_scheme_entry_count, errmsg, errflg)
    call t_stopf('ap_exp_sol_run')
    if (errflg /= 0) then
       write(iulog,*) trim(errmsg)
       call endrun
    end if

    ! Names and history remain CAM adapter responsibilities.  For the active
    ! pp_trop_mam3 mechanism clscnt1=0, so this loop is intentionally empty.
    do m = 1,clscnt1
       l = clsmap(m,1)
       do i = 1,ncol
          do k = ltrop(i)+1,pver
             base_sol(i,k,l)  = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
          end do
       end do

       wrk(:,:) = (prod(:,:,m) + ind_prd(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHMP', wrk(:,:), ncol, lchnk )
       wrk(:,:) = (loss(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHML', wrk(:,:), ncol, lchnk )
       
    end do

  end subroutine exp_sol

end module mo_exp_sol
