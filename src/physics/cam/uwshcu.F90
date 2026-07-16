module uwshcu

  ! Preserve CAM's original module API while the numerical process bodies and
  ! their shared initialized state live in the grouped scheme module.
  use ap_uwshcu_processes_scheme, only: uwshcu_readnl, init_uwshcu, &
       compute_uwshcu => compute_uwshcu_run, &
       compute_uwshcu_inv => compute_uwshcu_inv_run

  implicit none
  private

  public :: uwshcu_readnl
  public :: init_uwshcu
  public :: compute_uwshcu
  public :: compute_uwshcu_inv

end module uwshcu
