module ap_gw_drag_prof_scheme

  use gw_utils, only: r8
  use coords_1d, only: Coords1D

  implicit none
  private

  public :: gw_drag_prof_run

contains

  !> \section arg_table_gw_drag_prof_run Argument Table
  !! \htmlinclude gw_drag_prof_run.html
  !!
subroutine gw_drag_prof_run(pver, pverp, nconst, ngwv, nwave, ktop, &
     tau_0_ubc, dback, rog, alpha, taumin, tndmax, umcfac, ubmc2mn, &
     gravit, kwv, effkwv, ncol, p, src_level, tend_level, dt, &
     t,    &
     piln, rhoi,    nm,   ni,  ubm,  ubi,  xv,    yv,   &
     effgw,      c, kvtt, q,   dse,  tau,  utgw,  vtgw, &
     ttgw, qtgw, egwdffi,   gwut, dttdf, dttke, ro_adjust)

  !-----------------------------------------------------------------------
  ! Solve for the drag profile from the multiple gravity wave drag
  ! parameterization.
  ! 1. scan up from the wave source to determine the stress profile
  ! 2. scan down the stress profile to determine the tendencies
  !     => apply bounds to the tendency
  !          a. from wkb solution
  !          b. from computational stability constraints
  !     => adjust stress on interface below to reflect actual bounded
  !        tendency
  !-----------------------------------------------------------------------

  use gw_diffusion, only: gw_ediff, gw_diff_tend
  use linear_1d_operators, only: TriDiagDecomp

  !------------------------------Arguments--------------------------------
  integer, intent(in) :: pver, pverp, nconst, ngwv, nwave, ktop
  logical, intent(in) :: tau_0_ubc
  real(r8), intent(in) :: dback, rog, alpha(pverp)
  real(r8), intent(in) :: taumin, tndmax, umcfac, ubmc2mn
  real(r8), intent(in) :: gravit, kwv, effkwv
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Level from which gravity waves are propagated upward.
  integer, intent(in) :: src_level(ncol)
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(ncol)
  ! Using tend_level > src_level allows the orographic waves to prescribe
  ! wave propagation up to a certain level, but then allow wind tendencies
  ! and adjustments to tau below that level.

  ! Time step.
  real(r8), intent(in) :: dt

  ! Midpoint and interface temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Log of interface pressures.
  real(r8), intent(in) :: piln(ncol,pverp)
  ! Interface densities.
  real(r8), intent(in) :: rhoi(ncol,pverp)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver), ni(ncol,pverp)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(in) :: ubm(ncol,pver), ubi(ncol,pverp)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(in) :: xv(ncol), yv(ncol)
  ! Tendency efficiency.
  real(r8), intent(in) :: effgw(ncol)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,nwave)
  ! Molecular thermal diffusivity.
  real(r8), intent(in) :: kvtt(ncol,pverp)
  ! Constituent array.
  real(r8), intent(in) :: q(ncol,pver,nconst)
  ! Dry static energy.
  real(r8), intent(in) :: dse(ncol,pver)

  ! Wave Reynolds stress.
  real(r8), intent(inout) :: tau(ncol,nwave,pverp)
  ! Zonal/meridional wind tendencies.
  real(r8), intent(out) :: utgw(ncol,pver), vtgw(ncol,pver)
  ! Gravity wave heating tendency.
  real(r8), intent(out) :: ttgw(ncol,pver)
  ! Gravity wave constituent tendency.
  real(r8), intent(out) :: qtgw(ncol,pver,nconst)

  ! Effective gravity wave diffusivity at interfaces.
  real(r8), intent(out) :: egwdffi(ncol,pverp)

  ! Gravity wave wind tendency for each wave.
  real(r8), intent(out) :: gwut(ncol,pver,nwave)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8), intent(out) :: dttdf(ncol,pver)
  real(r8), intent(out) :: dttke(ncol,pver)

  ! Adjustment parameter for IGWs.
  real(r8), intent(in), optional :: &
       ro_adjust(ncol,nwave,pverp)

  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, and constituent loop indices.
  integer :: k, l, m

  ! Lowest tendency and source levels.
  integer :: kbot_tend, kbot_src

  ! "Total" and saturation diffusivity.
  real(r8) :: d(ncol)
  ! Imaginary part of vertical wavenumber.
  real(r8) :: mi(ncol)
  ! Stress after damping.
  real(r8) :: taudmp(ncol)
  ! Saturation stress.
  real(r8) :: tausat(ncol)
  ! (ub-c) and (ub-c)**2
  real(r8) :: ubmc(ncol), ubmc2(ncol)
  ! Temporary ubar tendencies (overall, and at wave l).
  real(r8) :: ubt(ncol,pver), ubtl(ncol)
  real(r8) :: wrk(ncol)
  ! Ratio used for ubt tndmax limiting.
  real(r8) :: ubt_lim_ratio(ncol)

  ! LU decomposition.
  type(TriDiagDecomp) :: decomp

  !------------------------------------------------------------------------

  ! Lowest levels that loops need to iterate over.
  kbot_tend = maxval(tend_level)
  kbot_src = maxval(src_level)

  ! Initialize gravity wave drag tendencies to zero.

  utgw = 0._r8
  vtgw = 0._r8

  gwut = 0._r8

  dttke = 0._r8
  ttgw = 0._r8

  ! Workaround floating point exception issues on Intel by initializing
  ! everything that's first set in a where block.
  mi = 0._r8
  taudmp = 0._r8
  tausat = 0._r8
  ubmc = 0._r8
  ubmc2 = 0._r8
  wrk = 0._r8

  !------------------------------------------------------------------------
  ! Compute the stress profiles and diffusivities
  !------------------------------------------------------------------------

  ! Loop from bottom to top to get stress profiles.
  do k = kbot_src, ktop, -1

     ! Determine the diffusivity for each column.

     d = dback + kvtt(:,k)

     do l = 1, nwave

        ! Determine the absolute value of the saturation stress.
        ! Define critical levels where the sign of (u-c) changes between
        ! interfaces.
        ubmc = ubi(:,k) - c(:,l)

        tausat = 0.0_r8
        where (src_level >= k)
           ! Test to see if u-c has the same sign here as the level below.
           where (ubmc > 0.0_r8 .eqv. ubi(:,k+1) > c(:,l))
              tausat = abs(effkwv * rhoi(:,k) * ubmc**3 / &
                   (2._r8*ni(:,k)))
           end where
        end where

        if (present(ro_adjust)) then
           where (src_level >= k)
              tausat = tausat * sqrt(ro_adjust(:,l,k))
           end where
        end if

        where (src_level >= k)

           ! Compute stress for each wave. The stress at this level is the
           ! min of the saturation stress and the stress at the level below
           ! reduced by damping. The sign of the stress must be the same as
           ! at the level below.

           ubmc2 = max(ubmc**2, ubmc2mn)
           mi = ni(:,k) / (2._r8 * kwv * ubmc2) * &
                (alpha(k) + ni(:,k)**2/ubmc2 * d)
           wrk = -2._r8*mi*rog*t(:,k)*(piln(:,k+1) - piln(:,k))

           taudmp = tau(:,l,k+1) * exp(wrk)

           ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
           ! we limit tau, so instead limit the arrays used to set it.
           where (tausat <= taumin) tausat = 0._r8
           where (taudmp <= taumin) taudmp = 0._r8

           tau(:,l,k) = min(taudmp, tausat)

        end where
     end do

  end do

  ! Force tau at the top of the model to zero, if requested.
  if (tau_0_ubc) tau(:,:,ktop) = 0._r8

  ! Apply efficiency to completed stress profile.
  do k = ktop, kbot_tend+1
     do l = 1, nwave
        where (k-1 <= tend_level)
           tau(:,l,k) = tau(:,l,k) * effgw
        end where
     end do
  end do

  !------------------------------------------------------------------------
  ! Compute the tendencies from the stress divergence.
  !------------------------------------------------------------------------

  ! Loop over levels from top to bottom
  do k = ktop, kbot_tend

     ! Accumulate the mean wind tendency over wavenumber.
     ubt(:,k) = 0.0_r8

     do l = 1, nwave    ! loop over wave

        ! Determine the wind tendency, including excess stress carried down
        ! from above.
        ubtl = gravit * (tau(:,l,k+1)-tau(:,l,k)) * p%rdel(:,k)

        ! Apply first tendency limit to maintain numerical stability.
        ! Enforce du/dt < |c-u|/dt  so u-c cannot change sign
        !    (u^n+1 = u^n + du/dt * dt)
        ! The limiter is somewhat stricter, so that we don't come anywhere
        ! near reversing c-u.
        ubtl = min(ubtl, umcfac * abs(c(:,l)-ubm(:,k)) / dt)

        where (k <= tend_level)

           ! Save tendency for each wave (for later computation of kzz):
           gwut(:,k,l) = sign(ubtl, c(:,l)-ubm(:,k))
           ubt(:,k) = ubt(:,k) + gwut(:,k,l)

        end where

     end do

     ! Apply second tendency limit to maintain numerical stability.
     ! Enforce du/dt < tndmax so that ridicuously large tendencies are not
     ! permitted.
     ! This can only happen above tend_level, so don't bother checking the
     ! level explicitly.
     where (abs(ubt(:,k)) > tndmax)
        ubt_lim_ratio = tndmax/abs(ubt(:,k))
        ubt(:,k) = ubt_lim_ratio * ubt(:,k)
     elsewhere
        ubt_lim_ratio = 1._r8
     end where

     do l = 1, nwave
        gwut(:,k,l) = ubt_lim_ratio*gwut(:,k,l)
        ! Redetermine the effective stress on the interface below from the
        ! wind tendency. If the wind tendency was limited above, then the
        ! new stress will be smaller than the old stress, causing stress
        ! divergence in the next layer down. This smoothes large stress
        ! divergences downward while conserving total stress.
        where (k <= tend_level)
           tau(:,l,k+1) = tau(:,l,k) + &
                abs(gwut(:,k,l)) * p%del(:,k) / gravit
        end where
     end do

     ! Project the mean wind tendency onto the components.
     where (k <= tend_level)
        utgw(:,k) = ubt(:,k) * xv
        vtgw(:,k) = ubt(:,k) * yv
     end where

     ! End of level loop.
  end do

  ! Calculate effective diffusivity and LU decomposition for the
  ! vertical diffusion solver.
  call gw_ediff (ncol, pver, ngwv, kbot_tend, ktop, tend_level, &
       gwut, ubm, nm, rhoi, dt, gravit, p, c, &
       egwdffi, decomp, ro_adjust=ro_adjust)

  ! Calculate tendency on each constituent.
  do m = 1, size(q,3)

     call gw_diff_tend(ncol, pver, kbot_tend, ktop, q(:,:,m), &
          dt, decomp, qtgw(:,:,m))

  enddo

  ! Calculate tendency from diffusing dry static energy (dttdf).
  call gw_diff_tend(ncol, pver, kbot_tend, ktop, dse, dt, decomp, dttdf)

  ! Evaluate second temperature tendency term: Conversion of kinetic
  ! energy into thermal.
  do l = 1, nwave
     do k = ktop, kbot_tend
        dttke(:,k) = dttke(:,k) - (ubm(:,k) - c(:,l)) * gwut(:,k,l)
     end do
  end do

  ttgw = dttke + dttdf

  ! Deallocate decomp.
  call decomp%finalize()

end subroutine gw_drag_prof_run

end module ap_gw_drag_prof_scheme
