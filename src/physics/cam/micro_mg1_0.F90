module micro_mg1_0

  use ap_micro_mg_tend_scheme, only: &
       micro_mg_init,                 &
       micro_mg_get_cols,             &
       micro_mg_tend => micro_mg_tend_run

  implicit none
  private

  public :: micro_mg_init
  public :: micro_mg_get_cols
  public :: micro_mg_tend

end module micro_mg1_0
