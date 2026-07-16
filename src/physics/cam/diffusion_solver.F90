module diffusion_solver

  use ap_compute_vdiff_scheme, only : scheme_init_vdiff => init_vdiff, &
       compute_vdiff_run, compute_vdiff_cam_adapter, &
       vdiff_selector, new_fieldlist_vdiff, &
       vdiff_select, vdiff_selector_to_flags, operator(.not.), any
  use perf_mod, only : t_startf, t_stopf

  implicit none
  private

  integer, parameter :: r8 = selected_real_kind(12)

  public :: init_vdiff
  public :: new_fieldlist_vdiff
  public :: compute_vdiff
  public :: vdiff_selector
  public :: vdiff_select
  public :: operator(.not.)
  public :: any

contains

  subroutine init_vdiff(kind, iulog_in, rair_in, gravit_in, do_iss_in, &
       errstring)

    integer,        intent(in)  :: kind
    integer,        intent(in)  :: iulog_in
    real(r8),       intent(in)  :: rair_in
    real(r8),       intent(in)  :: gravit_in
    logical,        intent(in)  :: do_iss_in
    character(128), intent(out) :: errstring

    call scheme_init_vdiff(kind, iulog_in, rair_in, gravit_in, do_iss_in, &
         errstring)

  end subroutine init_vdiff

  subroutine compute_vdiff( lchnk, pcols, pver, ncnst, ncol, pmid, pint, &
       pdel, rpdel, t, ztodt, taux, tauy, shflx, cflx, ntop, nbot, kvh, &
       kvm, kvq, cgs, cgh, zi, ksrftms, qmincg, fieldlist, fieldlistm, &
       u, v, q, dse, tautmsx, tautmsy, dtk, topflx, errstring, tauresx, &
       tauresy, itaures, cpairv, rairi, do_molec_diff, &
       compute_molec_diff, vd_lu_qdecomp, kvt )

    use molec_diff, only : host_compute_molec_diff => compute_molec_diff, &
         host_vd_lu_qdecomp => vd_lu_qdecomp

    integer,  intent(in) :: lchnk
    integer,  intent(in) :: pcols
    integer,  intent(in) :: pver
    integer,  intent(in) :: ncnst
    integer,  intent(in) :: ncol
    integer,  intent(in) :: ntop
    integer,  intent(in) :: nbot
    integer,  intent(in) :: itaures
    real(r8), intent(in) :: pmid(pcols,pver)
    real(r8), intent(in) :: pint(pcols,pver+1)
    real(r8), intent(in) :: pdel(pcols,pver)
    real(r8), intent(in) :: rpdel(pcols,pver)
    real(r8), intent(in) :: t(pcols,pver)
    real(r8), intent(in) :: ztodt
    real(r8), intent(in) :: taux(pcols)
    real(r8), intent(in) :: tauy(pcols)
    real(r8), intent(in) :: shflx(pcols)
    real(r8), intent(in) :: cflx(pcols,ncnst)
    real(r8), intent(in) :: kvh(pcols,pver+1)
    real(r8), intent(in) :: zi(pcols,pver+1)
    real(r8), intent(in) :: ksrftms(pcols)
    real(r8), intent(in) :: qmincg(ncnst)
    real(r8), intent(in) :: cpairv(pcols,pver)
    real(r8), intent(in) :: rairi(pcols,pver+1)
    logical,  intent(in) :: do_molec_diff
    type(vdiff_selector), intent(in) :: fieldlist
    type(vdiff_selector), intent(in) :: fieldlistm

    real(r8), intent(inout) :: kvm(pcols,pver+1)
    real(r8), intent(inout) :: kvq(pcols,pver+1)
    real(r8), intent(inout) :: cgs(pcols,pver+1)
    real(r8), intent(inout) :: cgh(pcols,pver+1)
    real(r8), intent(inout) :: u(pcols,pver)
    real(r8), intent(inout) :: v(pcols,pver)
    real(r8), intent(inout) :: q(pcols,pver,ncnst)
    real(r8), intent(inout) :: dse(pcols,pver)
    real(r8), intent(inout) :: tauresx(pcols)
    real(r8), intent(inout) :: tauresy(pcols)

    real(r8), intent(out) :: tautmsx(pcols)
    real(r8), intent(out) :: tautmsy(pcols)
    real(r8), intent(out) :: dtk(pcols,pver)
    real(r8), intent(out) :: topflx(pcols)
    character(128), intent(out) :: errstring
    real(r8), intent(out), optional :: kvt(pcols,pver+1)
    procedure(), optional :: compute_molec_diff
    procedure(), optional :: vd_lu_qdecomp

    character(len=512) :: errmsg
    character(len=128) :: core_errstring
    integer :: errflg
    logical :: field_flags(3+ncnst)
    logical :: molecular_field_flags(3+ncnst)

    call vdiff_selector_to_flags(fieldlist, field_flags)
    call vdiff_selector_to_flags(fieldlistm, molecular_field_flags)

    errmsg = ''
    core_errstring = ''
    errflg = 0

    call t_startf('ap_compute_vdiff_run')
    if (do_molec_diff) then
       if (present(kvt)) then
          call compute_vdiff_cam_adapter(lchnk, pcols, pver, ncnst, ncol, &
               pmid, pint, pdel, rpdel, t, ztodt, taux, tauy, shflx, cflx, &
               ntop, nbot, kvh, kvm, kvq, cgs, cgh, zi, ksrftms, qmincg, &
               field_flags, molecular_field_flags, u, v, q, dse, tautmsx, &
               tautmsy, dtk, topflx, core_errstring, tauresx, tauresy, &
               itaures, cpairv, rairi, do_molec_diff, &
               host_compute_molec_diff, host_vd_lu_qdecomp, kvt)
       else
          call compute_vdiff_cam_adapter(lchnk, pcols, pver, ncnst, ncol, &
               pmid, pint, pdel, rpdel, t, ztodt, taux, tauy, shflx, cflx, &
               ntop, nbot, kvh, kvm, kvq, cgs, cgh, zi, ksrftms, qmincg, &
               field_flags, molecular_field_flags, u, v, q, dse, tautmsx, &
               tautmsy, dtk, topflx, core_errstring, tauresx, tauresy, &
               itaures, cpairv, rairi, do_molec_diff, &
               host_compute_molec_diff, host_vd_lu_qdecomp)
       end if
       errmsg = trim(core_errstring)
       if (len_trim(core_errstring) > 0) errflg = 1
    else if (present(kvt)) then
       call compute_vdiff_run(lchnk, pcols, pver, pver+1, ncnst, 3+ncnst, ncol, pmid, &
            pint, pdel, rpdel, t, ztodt, taux, tauy, shflx, cflx, ntop, &
            nbot, kvh, kvm, kvq, cgs, cgh, zi, ksrftms, qmincg, &
            field_flags, molecular_field_flags, u, v, q, dse, tautmsx, tautmsy, dtk, &
            topflx, tauresx, tauresy, itaures, cpairv, rairi, &
            do_molec_diff, kvt, errmsg, errflg)
    else
       call compute_vdiff_run(lchnk, pcols, pver, pver+1, ncnst, 3+ncnst, ncol, pmid, &
            pint, pdel, rpdel, t, ztodt, taux, tauy, shflx, cflx, ntop, &
            nbot, kvh, kvm, kvq, cgs, cgh, zi, ksrftms, qmincg, &
            field_flags, molecular_field_flags, u, v, q, dse, tautmsx, tautmsy, dtk, &
            topflx, tauresx, tauresy, itaures, cpairv, rairi, &
            do_molec_diff, errmsg=errmsg, errflg=errflg)
    end if
    call t_stopf('ap_compute_vdiff_run')

    errstring = trim(errmsg)

  end subroutine compute_vdiff

end module diffusion_solver
