module ap_imp_sol_trop_mam3_scheme

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  ! These are generated-mechanism constants, not general CAM dimensions.
  integer, parameter :: mechanism_gas_pcnst = 20
  integer, parameter :: mechanism_rxntot = 7
  integer, parameter :: mechanism_extcnt = 7
  integer, parameter :: mechanism_clscnt4 = 20
  integer, parameter :: mechanism_nzcnt = 22
  integer, parameter :: mechanism_itermax = 11
  integer, parameter :: mechanism_cut_limit = 5

  public :: imp_sol_run

contains

  !> \section arg_table_imp_sol_run Argument Table
  !! \htmlinclude imp_sol_run.html
  !!
  subroutine imp_sol_run(ncol, num_levels, number_gas_species, number_reactions, &
       number_external_forcings, number_implicit_species, number_jacobian_entries, &
       maximum_iterations, base_sol, reaction_rates, het_rates, extfrc, delt, &
       xhnm, ltrop, clsmap4, permute4, class_independent_reaction_count, &
       epsilon, factor, small, prod_out, loss_out, solver_failure_count, &
       last_failure_dt, entry_count, errmsg, errflg)

    integer, intent(in) :: ncol, num_levels
    integer, intent(in) :: number_gas_species, number_reactions
    integer, intent(in) :: number_external_forcings, number_implicit_species
    integer, intent(in) :: number_jacobian_entries, maximum_iterations
    real(r8), intent(inout) :: base_sol(ncol,num_levels,number_gas_species)
    real(r8), intent(in) :: reaction_rates(ncol,num_levels,number_reactions)
    real(r8), intent(in) :: het_rates(ncol,num_levels,number_gas_species)
    real(r8), intent(in) :: extfrc(ncol,num_levels,number_external_forcings)
    real(r8), intent(in) :: delt
    real(r8), intent(in) :: xhnm(ncol,num_levels)
    integer, intent(in) :: ltrop(ncol)
    integer, intent(in) :: clsmap4(number_implicit_species)
    integer, intent(in) :: permute4(number_implicit_species)
    integer, intent(in) :: class_independent_reaction_count
    real(r8), intent(in) :: epsilon(number_implicit_species)
    logical, intent(in) :: factor(maximum_iterations)
    real(r8), intent(in) :: small
    real(r8), intent(out) :: prod_out(ncol,num_levels,number_implicit_species)
    real(r8), intent(out) :: loss_out(ncol,num_levels,number_implicit_species)
    integer, intent(out) :: solver_failure_count
    real(r8), intent(out) :: last_failure_dt
    integer, intent(inout) :: entry_count
    character(len=512), intent(out) :: errmsg
    integer, intent(out) :: errflg

    integer :: nr_iter, lev, i, j, k, m
    integer :: fail_cnt, cut_cnt, stp_con_cnt
    real(r8) :: interval_done, dt, dti
    real(r8) :: max_delta(mechanism_clscnt4)
    real(r8) :: sys_jac(mechanism_nzcnt)
    real(r8) :: lin_jac(mechanism_nzcnt)
    real(r8) :: solution(mechanism_clscnt4)
    real(r8) :: forcing(mechanism_clscnt4)
    real(r8) :: iter_invariant(mechanism_clscnt4)
    real(r8) :: prod(mechanism_clscnt4)
    real(r8) :: loss(mechanism_clscnt4)
    real(r8) :: lrxt(mechanism_rxntot)
    real(r8) :: lsol(mechanism_gas_pcnst)
    real(r8) :: lhet(mechanism_gas_pcnst)
    real(r8) :: ind_prd(ncol,num_levels,mechanism_clscnt4)
    logical :: convergence
    logical :: frc_mask
    logical :: converged(mechanism_clscnt4)

    entry_count = entry_count + 1
    errmsg = ''
    errflg = 0
    solver_failure_count = 0
    last_failure_dt = 0._r8

    if (number_gas_species /= mechanism_gas_pcnst .or. &
        number_reactions /= mechanism_rxntot .or. &
        number_external_forcings /= mechanism_extcnt .or. &
        number_implicit_species /= mechanism_clscnt4 .or. &
        number_jacobian_entries /= mechanism_nzcnt .or. &
        maximum_iterations /= mechanism_itermax) then
       errmsg = 'imp_sol_run dimension mismatch: expected pp_trop_mam3 mechanism'
       errflg = 1
       prod_out = 0._r8
       loss_out = 0._r8
       return
    end if

    prod_out(:,:,:) = 0._r8
    loss_out(:,:,:) = 0._r8
    solution(:) = 0._r8

    ! Class-independent forcing.  Keep the generated helper call before the
    ! level/column solve, as in mo_imp_sol.
    if (class_independent_reaction_count > 0 .or. &
        number_external_forcings > 0) then
       call indprd(4, ind_prd, mechanism_clscnt4, base_sol, extfrc, &
            reaction_rates, ncol, num_levels)
    else
       do m = 1,mechanism_clscnt4
          ind_prd(:,:,m) = 0._r8
       end do
    end if

    level_loop : do lev = 1,num_levels
       column_loop : do i = 1,ncol
          if (lev <= ltrop(i)) cycle column_loop

          do m = 1,mechanism_rxntot
             lrxt(m) = reaction_rates(i,lev,m)
          end do
          do m = 1,mechanism_gas_pcnst
             lhet(m) = het_rates(i,lev,m)
          end do

          dt = delt
          cut_cnt = 0
          fail_cnt = 0
          stp_con_cnt = 0
          interval_done = 0._r8

          time_step_loop : do
             dti = 1._r8 / dt

             do m = 1,mechanism_gas_pcnst
                lsol(m) = base_sol(i,lev,m)
             end do

             do k = 1,mechanism_clscnt4
                j = clsmap4(k)
                m = permute4(k)
                solution(m) = lsol(j)
             end do

             if (class_independent_reaction_count > 0 .or. &
                 number_external_forcings > 0) then
                do m = 1,mechanism_clscnt4
                   iter_invariant(m) = dti * solution(m) + ind_prd(i,lev,m)
                end do
             else
                do m = 1,mechanism_clscnt4
                   iter_invariant(m) = dti * solution(m)
                end do
             end if

             call linmat(lin_jac, lsol, lrxt, lhet)

             iter_loop : do nr_iter = 1,mechanism_itermax
                if (factor(nr_iter)) then
                   call nlnmat(sys_jac, lsol, lrxt, lin_jac, dti)
                   call lu_fac(sys_jac)
                end if

                call imp_prod_loss(prod, loss, lsol, lrxt, lhet)
                do m = 1,mechanism_clscnt4
                   forcing(m) = solution(m)*dti - &
                        (iter_invariant(m) + prod(m) - loss(m))
                end do

                call lu_slv(sys_jac, forcing)
                do m = 1,mechanism_clscnt4
                   solution(m) = solution(m) + forcing(m)
                end do

                if (nr_iter > 1) then
                   do k = 1,mechanism_clscnt4
                      m = permute4(k)
                      if (abs(solution(m)) > 1.e-20_r8) then
                         max_delta(k) = abs(forcing(m)/solution(m))
                      else
                         max_delta(k) = 0._r8
                      end if
                   end do
                end if

                where (solution(:) < 0._r8)
                   solution(:) = 0._r8
                end where

                do k = 1,mechanism_clscnt4
                   j = clsmap4(k)
                   m = permute4(k)
                   lsol(j) = solution(m)
                end do

                converged(:) = .true.
                if (nr_iter > 1) then
                   do k = 1,mechanism_clscnt4
                      m = permute4(k)
                      frc_mask = abs(forcing(m)) > small
                      if (frc_mask) then
                         converged(k) = abs(forcing(m)) <= &
                              epsilon(k)*abs(solution(m))
                      else
                         converged(k) = .true.
                      end if
                   end do
                   convergence = all(converged(:))
                   if (convergence) exit
                end if
             end do iter_loop

             if (.not. convergence) then
                fail_cnt = fail_cnt + 1
                solver_failure_count = solver_failure_count + 1
                last_failure_dt = dt
                stp_con_cnt = 0
                if (cut_cnt < mechanism_cut_limit) then
                   cut_cnt = cut_cnt + 1
                   if (cut_cnt < mechanism_cut_limit) then
                      dt = .5_r8 * dt
                   else
                      dt = .1_r8 * dt
                   end if
                   cycle time_step_loop
                end if
             end if

             interval_done = interval_done + dt
             if (abs(delt - interval_done) <= .0001_r8) then
                exit time_step_loop
             else
                if (convergence) then
                   stp_con_cnt = stp_con_cnt + 1
                end if
                do m = 1,mechanism_gas_pcnst
                   base_sol(i,lev,m) = lsol(m)
                end do
                if (stp_con_cnt >= 2) then
                   dt = 2._r8*dt
                   stp_con_cnt = 0
                end if
                dt = min(dt,delt-interval_done)
             end if
          end do time_step_loop

          do k = 1,mechanism_clscnt4
             j = clsmap4(k)
             m = permute4(k)
             base_sol(i,lev,j) = solution(m)
          end do

          ! pp_trop_mam3 has no OX/O3 implicit diagnostic branch.  Therefore
          ! the generated mechanism always takes the original else branch.
          do k = 1,mechanism_clscnt4
             m = permute4(k)
             prod_out(i,lev,k) = prod(m) + ind_prd(i,lev,m)
             loss_out(i,lev,k) = loss(m)
          end do
       end do column_loop
    end do level_loop

    do i = 1,mechanism_clscnt4
       prod_out(:,:,i) = prod_out(:,:,i)*xhnm
       loss_out(:,:,i) = loss_out(:,:,i)*xhnm
    end do
  end subroutine imp_sol_run


  ! The following private helpers are the generated pp_trop_mam3 solver
  ! kernels.  Their statement and call ordering intentionally matches the
  ! generated chemistry source; they are not a public process API.
  subroutine indprd(class, prod, nprod, y, extfrc, rxt, ncol, num_levels)
    integer, intent(in) :: class, nprod, ncol, num_levels
    real(r8), intent(in) :: y(ncol,num_levels,mechanism_gas_pcnst)
    real(r8), intent(in) :: rxt(ncol,num_levels,mechanism_rxntot)
    real(r8), intent(in) :: extfrc(ncol,num_levels,mechanism_extcnt)
    real(r8), intent(inout) :: prod(ncol,num_levels,nprod)

    if (class == 4) then
       prod(:,:,1) = rxt(:,:,2)
       prod(:,:,2) = 0._r8
       prod(:,:,3) = + extfrc(:,:,1)
       prod(:,:,4) = 0._r8
       prod(:,:,5) = 0._r8
       prod(:,:,6) = + extfrc(:,:,2)
       prod(:,:,7) = + extfrc(:,:,4)
       prod(:,:,8) = 0._r8
       prod(:,:,9) = + extfrc(:,:,5)
       prod(:,:,10) = 0._r8
       prod(:,:,11) = 0._r8
       prod(:,:,12) = + extfrc(:,:,6)
       prod(:,:,13) = + extfrc(:,:,3)
       prod(:,:,14) = 0._r8
       prod(:,:,15) = 0._r8
       prod(:,:,16) = + extfrc(:,:,7)
       prod(:,:,17) = 0._r8
       prod(:,:,18) = 0._r8
       prod(:,:,19) = 0._r8
       prod(:,:,20) = 0._r8
    end if
  end subroutine indprd


  subroutine linmat(mat, y, rxt, het_rates)
    real(r8), intent(in) :: y(mechanism_gas_pcnst)
    real(r8), intent(in) :: rxt(mechanism_rxntot)
    real(r8), intent(in) :: het_rates(mechanism_gas_pcnst)
    real(r8), intent(inout) :: mat(mechanism_nzcnt)

    call linmat01(mat, y, rxt, het_rates)
  end subroutine linmat


  subroutine linmat01(mat, y, rxt, het_rates)
    real(r8), intent(in) :: y(mechanism_gas_pcnst)
    real(r8), intent(in) :: rxt(mechanism_rxntot)
    real(r8), intent(in) :: het_rates(mechanism_gas_pcnst)
    real(r8), intent(inout) :: mat(mechanism_nzcnt)

    mat(1) = -(rxt(1) + rxt(3) + het_rates(1))
    mat(2) = -(het_rates(2))
    mat(3) = rxt(4)
    mat(4) = -(rxt(4) + het_rates(3))
    mat(5) = rxt(5) + .500_r8*rxt(6) + rxt(7)
    mat(6) = -(rxt(5) + rxt(6) + rxt(7) + het_rates(4))
    mat(7) = -(het_rates(5))
    mat(8) = -(het_rates(6))
    mat(9) = -(het_rates(7))
    mat(10) = -(het_rates(8))
    mat(11) = -(het_rates(9))
    mat(12) = -(het_rates(10))
    mat(13) = -(het_rates(11))
    mat(14) = -(het_rates(12))
    mat(15) = -(het_rates(13))
    mat(16) = -(het_rates(14))
    mat(17) = -(het_rates(15))
    mat(18) = -(het_rates(16))
    mat(19) = -(het_rates(17))
    mat(20) = -(het_rates(18))
    mat(21) = -(het_rates(19))
    mat(22) = -(het_rates(20))
  end subroutine linmat01


  subroutine nlnmat(mat, y, rxt, lmat, dti)
    real(r8), intent(in) :: dti
    real(r8), intent(in) :: lmat(mechanism_nzcnt)
    real(r8), intent(in) :: y(mechanism_gas_pcnst)
    real(r8), intent(in) :: rxt(mechanism_rxntot)
    real(r8), intent(inout) :: mat(mechanism_nzcnt)

    call nlnmat_finit(mat, lmat, dti)
  end subroutine nlnmat


  subroutine nlnmat_finit(mat, lmat, dti)
    real(r8), intent(in) :: dti
    real(r8), intent(in) :: lmat(mechanism_nzcnt)
    real(r8), intent(inout) :: mat(mechanism_nzcnt)
    mat(1) = lmat(1)
    mat(2) = lmat(2)
    mat(3) = lmat(3)
    mat(4) = lmat(4)
    mat(5) = lmat(5)
    mat(6) = lmat(6)
    mat(7) = lmat(7)
    mat(8) = lmat(8)
    mat(9) = lmat(9)
    mat(10) = lmat(10)
    mat(11) = lmat(11)
    mat(12) = lmat(12)
    mat(13) = lmat(13)
    mat(14) = lmat(14)
    mat(15) = lmat(15)
    mat(16) = lmat(16)
    mat(17) = lmat(17)
    mat(18) = lmat(18)
    mat(19) = lmat(19)
    mat(20) = lmat(20)
    mat(21) = lmat(21)
    mat(22) = lmat(22)
    mat(1) = mat(1) - dti
    mat(2) = mat(2) - dti
    mat(4) = mat(4) - dti
    mat(6) = mat(6) - dti
    mat(7) = mat(7) - dti
    mat(8) = mat(8) - dti
    mat(9) = mat(9) - dti
    mat(10) = mat(10) - dti
    mat(11) = mat(11) - dti
    mat(12) = mat(12) - dti
    mat(13) = mat(13) - dti
    mat(14) = mat(14) - dti
    mat(15) = mat(15) - dti
    mat(16) = mat(16) - dti
    mat(17) = mat(17) - dti
    mat(18) = mat(18) - dti
    mat(19) = mat(19) - dti
    mat(20) = mat(20) - dti
    mat(21) = mat(21) - dti
    mat(22) = mat(22) - dti
  end subroutine nlnmat_finit


  subroutine lu_fac(lu)
    real(r8), intent(inout) :: lu(:)

    call lu_fac01(lu)
  end subroutine lu_fac


  subroutine lu_fac01(lu)
    real(r8), intent(inout) :: lu(:)

    lu(1) = 1._r8 / lu(1)
    lu(2) = 1._r8 / lu(2)
    lu(4) = 1._r8 / lu(4)
    lu(6) = 1._r8 / lu(6)
    lu(7) = 1._r8 / lu(7)
    lu(8) = 1._r8 / lu(8)
    lu(9) = 1._r8 / lu(9)
    lu(10) = 1._r8 / lu(10)
    lu(11) = 1._r8 / lu(11)
    lu(12) = 1._r8 / lu(12)
    lu(13) = 1._r8 / lu(13)
    lu(14) = 1._r8 / lu(14)
    lu(15) = 1._r8 / lu(15)
    lu(16) = 1._r8 / lu(16)
    lu(17) = 1._r8 / lu(17)
    lu(18) = 1._r8 / lu(18)
    lu(19) = 1._r8 / lu(19)
    lu(20) = 1._r8 / lu(20)
    lu(21) = 1._r8 / lu(21)
    lu(22) = 1._r8 / lu(22)
  end subroutine lu_fac01


  subroutine lu_slv(lu, b)
    real(r8), intent(in) :: lu(:)
    real(r8), intent(inout) :: b(:)

    call lu_slv01(lu, b)
  end subroutine lu_slv


  subroutine lu_slv01(lu, b)
    real(r8), intent(in) :: lu(:)
    real(r8), intent(inout) :: b(:)

    b(20) = b(20) * lu(22)
    b(19) = b(19) * lu(21)
    b(18) = b(18) * lu(20)
    b(17) = b(17) * lu(19)
    b(16) = b(16) * lu(18)
    b(15) = b(15) * lu(17)
    b(14) = b(14) * lu(16)
    b(13) = b(13) * lu(15)
    b(12) = b(12) * lu(14)
    b(11) = b(11) * lu(13)
    b(10) = b(10) * lu(12)
    b(9) = b(9) * lu(11)
    b(8) = b(8) * lu(10)
    b(7) = b(7) * lu(9)
    b(6) = b(6) * lu(8)
    b(5) = b(5) * lu(7)
    b(4) = b(4) * lu(6)
    b(3) = b(3) - lu(5) * b(4)
    b(3) = b(3) * lu(4)
    b(2) = b(2) - lu(3) * b(3)
    b(2) = b(2) * lu(2)
    b(1) = b(1) * lu(1)
  end subroutine lu_slv01


  subroutine imp_prod_loss(prod, loss, y, rxt, het_rates)
    real(r8), intent(out) :: prod(:), loss(:)
    real(r8), intent(in) :: y(:), rxt(:), het_rates(:)

    loss(1) = (+ rxt(1) + rxt(3) + het_rates(1))*y(1)
    prod(1) = 0._r8
    loss(2) = (+ het_rates(2))*y(2)
    prod(2) = rxt(4)*y(3)
    loss(3) = (+ rxt(4) + het_rates(3))*y(3)
    prod(3) = (rxt(5) +.500_r8*rxt(6) +rxt(7))*y(4)
    loss(4) = (+ rxt(5) + rxt(6) + rxt(7) + het_rates(4))*y(4)
    prod(4) = 0._r8
    loss(5) = (+ het_rates(5))*y(5)
    prod(5) = 0._r8
    loss(6) = (+ het_rates(6))*y(6)
    prod(6) = 0._r8
    loss(7) = (+ het_rates(7))*y(7)
    prod(7) = 0._r8
    loss(8) = (+ het_rates(8))*y(8)
    prod(8) = 0._r8
    loss(9) = (+ het_rates(9))*y(9)
    prod(9) = 0._r8
    loss(10) = (+ het_rates(10))*y(10)
    prod(10) = 0._r8
    loss(11) = (+ het_rates(11))*y(11)
    prod(11) = 0._r8
    loss(12) = (+ het_rates(12))*y(12)
    prod(12) = 0._r8
    loss(13) = (+ het_rates(13))*y(13)
    prod(13) = 0._r8
    loss(14) = (+ het_rates(14))*y(14)
    prod(14) = 0._r8
    loss(15) = (+ het_rates(15))*y(15)
    prod(15) = 0._r8
    loss(16) = (+ het_rates(16))*y(16)
    prod(16) = 0._r8
    loss(17) = (+ het_rates(17))*y(17)
    prod(17) = 0._r8
    loss(18) = (+ het_rates(18))*y(18)
    prod(18) = 0._r8
    loss(19) = (+ het_rates(19))*y(19)
    prod(19) = 0._r8
    loss(20) = (+ het_rates(20))*y(20)
    prod(20) = 0._r8
  end subroutine imp_prod_loss

end module ap_imp_sol_trop_mam3_scheme
