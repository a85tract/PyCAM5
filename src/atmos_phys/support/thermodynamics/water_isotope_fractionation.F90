module ap_water_isotope_fractionation

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ap_water_isotope_coefficients, only: wiso_alpl, wiso_alpi, &
       wiso_akel, wiso_akci

  implicit none
  private

  abstract interface
     function water_isotope_alpha_hook(ispec, isrctype, idsttype, tk, rh, &
          do_kinetic) result(alpha)
       import :: r8
       integer, intent(in) :: ispec, isrctype, idsttype
       real(r8), intent(in) :: tk, rh
       logical, intent(in) :: do_kinetic
       real(r8) :: alpha
     end function water_isotope_alpha_hook
  end interface

  procedure(water_isotope_alpha_hook), pointer :: alpha_hook => null()

  public :: water_isotope_alpha
  public :: register_water_isotope_alpha_hook

contains

  subroutine register_water_isotope_alpha_hook(alpha_proc)
    procedure(water_isotope_alpha_hook) :: alpha_proc

    alpha_hook => alpha_proc
  end subroutine register_water_isotope_alpha_hook

  function water_isotope_alpha(ispec, isrctype, idsttype, tk, rh, &
       do_kinetic) result(alpha)
    integer, intent(in) :: ispec, isrctype, idsttype
    real(r8), intent(in) :: tk, rh
    logical, intent(in) :: do_kinetic
    real(r8) :: alpha

    if (associated(alpha_hook)) then
       alpha = alpha_hook(ispec, isrctype, idsttype, tk, rh, do_kinetic)
       return
    end if

    alpha = fallback_water_isotope_alpha(ispec, isrctype, idsttype, tk, &
         rh, do_kinetic)
  end function water_isotope_alpha

  function fallback_water_isotope_alpha(ispec, isrctype, idsttype, tk, &
       rh, do_kinetic) result(alpha)
    integer, intent(in) :: ispec, isrctype, idsttype
    real(r8), intent(in) :: tk, rh
    logical, intent(in) :: do_kinetic
    real(r8) :: alpha
    integer, parameter :: iwtvap = 1
    integer, parameter :: iwtliq = 2
    integer, parameter :: iwtstrain = 4
    integer, parameter :: iwtcvrain = 6

    alpha = 1._r8

    if (isrctype /= idsttype) then
       if (isrctype == iwtvap) then
          if (idsttype == iwtliq .or. idsttype == iwtstrain .or. &
              idsttype == iwtcvrain) then
             alpha = wiso_alpl(ispec, tk)
             if (do_kinetic) alpha = wiso_akel(ispec, tk, rh, alpha)
          else
             alpha = wiso_alpi(ispec, tk)
             if (do_kinetic) alpha = wiso_akci(ispec, tk, alpha)
          end if
       else if (idsttype == iwtvap) then
          if (isrctype == iwtliq .or. isrctype == iwtstrain .or. &
              isrctype == iwtcvrain) then
             alpha = wiso_alpl(ispec, tk)
             if (do_kinetic) alpha = wiso_akel(ispec, tk, rh, alpha)
             alpha = 1._r8/alpha
          else
             alpha = 1._r8
          end if
       end if
    end if

  end function fallback_water_isotope_alpha

end module ap_water_isotope_fractionation
