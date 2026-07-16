module ndrop

  use ap_dropmixnuc_scheme, only: &
       ndrop_init,                 &
       dropmixnuc => dropmixnuc_run

  implicit none
  private

  public :: ndrop_init
  public :: dropmixnuc

end module ndrop
