module ndrop

  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst
  use cam_history, only: fieldname_len
  use physics_types, only: physics_state, physics_ptend
  use physics_buffer, only: physics_buffer_desc
  use ap_dropmixnuc_scheme, only: dropmixnuc_init, dropmixnuc_run
  use dropmixnuc_host_hooks, only: register_dropmixnuc_host_hooks

  implicit none
  private
  save

  integer :: numliq_idx = -1
  integer :: kvh_idx = -1
  integer :: ntot_amode, ncnst_tot
  integer, allocatable :: nspec_amode(:)
  integer, allocatable :: mam_idx(:,:), mam_cnst_idx(:,:)
  real(r8), allocatable :: sigmag_amode(:), dgnumlo_amode(:), dgnumhi_amode(:)
  real(r8), allocatable :: specdens_amode(:,:), spechygro_amode(:,:)
  character(len=fieldname_len), allocatable :: fieldname(:), fieldname_cw(:)
  logical :: prog_modal_aero
  logical :: lq(pcnst) = .false.

  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  public :: ndrop_init
  public :: dropmixnuc

contains

  subroutine ndrop_init
    use cam_history, only: addfld, add_default, phys_decomp
    use cam_logfile, only: iulog
    use constituents, only: cnst_get_ind
    use phys_control, only: phys_getopts
    use physconst, only: pi, rhoh2o, mwh2o, r_universal, rh2o, &
         gravit, latvap, cpair, rair
    use physics_buffer, only: pbuf_get_index
    use ppgrid, only: pcols, pver, pverp
    use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_props, &
         rad_cnst_get_aer_props, rad_cnst_get_mam_mmr_idx, &
         rad_cnst_get_mode_num_idx
    use ref_pres, only: top_lev => trop_cloud_top_lev

    integer :: ii, l, lptr, m, mm, nspec_max
    character(len=32) :: tmpname, tmpname_cw
    character(len=128) :: long_name
    character(len=8) :: unit
    logical :: history_amwg, history_aerosol

    call cnst_get_ind('NUMLIQ', numliq_idx)
    kvh_idx = pbuf_get_index('kvh')

    call rad_cnst_get_info(0, nmodes=ntot_amode)
    allocate(nspec_amode(ntot_amode), sigmag_amode(ntot_amode), &
         dgnumlo_amode(ntot_amode), dgnumhi_amode(ntot_amode))

    do m = 1, ntot_amode
       call rad_cnst_get_info(0, m, nspec=nspec_amode(m))
       call rad_cnst_get_mode_props(0, m, sigmag=sigmag_amode(m), &
            dgnumhi=dgnumhi_amode(m), dgnumlo=dgnumlo_amode(m))
    end do

    nspec_max = nspec_amode(1)
    ncnst_tot = nspec_amode(1) + 1
    do m = 2, ntot_amode
       nspec_max = max(nspec_max, nspec_amode(m))
       ncnst_tot = ncnst_tot + nspec_amode(m) + 1
    end do

    allocate(mam_idx(ntot_amode,0:nspec_max), &
         mam_cnst_idx(ntot_amode,0:nspec_max), &
         specdens_amode(ntot_amode,nspec_max), &
         spechygro_amode(ntot_amode,nspec_max), &
         fieldname(ncnst_tot), fieldname_cw(ncnst_tot))
    mam_idx = 0
    mam_cnst_idx = 0
    specdens_amode = 0._r8
    spechygro_amode = 0._r8

    ii = 0
    do m = 1, ntot_amode
       do l = 0, nspec_amode(m)
          ii = ii + 1
          mam_idx(m,l) = ii
       end do
    end do

    call phys_getopts(history_amwg_out=history_amwg, &
         history_aerosol_out=history_aerosol, &
         prog_modal_aero_out=prog_modal_aero)

    do m = 1, ntot_amode
       do l = 0, nspec_amode(m)
          mm = mam_idx(m,l)
          unit = 'kg/m2/s'
          if (l == 0) unit = '#/m2/s'

          if (l == 0) then
             call rad_cnst_get_info(0, m, num_name=tmpname, &
                  num_name_cw=tmpname_cw)
          else
             call rad_cnst_get_info(0, m, l, spec_name=tmpname, &
                  spec_name_cw=tmpname_cw)
             call rad_cnst_get_aer_props(0, m, l, &
                  density_aer=specdens_amode(m,l), &
                  hygro_aer=spechygro_amode(m,l))
          end if

          fieldname(mm) = trim(tmpname)//'_mixnuc1'
          fieldname_cw(mm) = trim(tmpname_cw)//'_mixnuc1'

          if (prog_modal_aero) then
             if (l == 0) then
                call rad_cnst_get_mode_num_idx(m, lptr)
             else
                call rad_cnst_get_mam_mmr_idx(m, l, lptr)
             end if
             mam_cnst_idx(m,l) = lptr
             lq(lptr) = .true.

             long_name = trim(tmpname)//' dropmixnuc mixnuc column tendency'
             call addfld(fieldname(mm), unit, 1, 'A', long_name, phys_decomp)
             long_name = trim(tmpname_cw)//' dropmixnuc mixnuc column tendency'
             call addfld(fieldname_cw(mm), unit, 1, 'A', long_name, phys_decomp)
             if (history_aerosol) then
                call add_default(fieldname(mm), 1, ' ')
                call add_default(fieldname_cw(mm), 1, ' ')
             end if
          end if
       end do
    end do

    call addfld('CCN1    ','#/cm3   ',pver,'A', &
         'CCN concentration at S=0.02%',phys_decomp)
    call addfld('CCN2    ','#/cm3   ',pver,'A', &
         'CCN concentration at S=0.05%',phys_decomp)
    call addfld('CCN3    ','#/cm3   ',pver,'A', &
         'CCN concentration at S=0.1%',phys_decomp)
    call addfld('CCN4    ','#/cm3   ',pver,'A', &
         'CCN concentration at S=0.2%',phys_decomp)
    call addfld('CCN5    ','#/cm3   ',pver,'A', &
         'CCN concentration at S=0.5%',phys_decomp)
    call addfld('CCN6    ','#/cm3   ',pver,'A', &
         'CCN concentration at S=1.0%',phys_decomp)
    call addfld('WTKE     ','m/s     ',pver,'A', &
         'Standard deviation of updraft velocity',phys_decomp)
    call addfld('NDROPMIX ','#/kg/s  ',pver,'A', &
         'Droplet number mixing',phys_decomp)
    call addfld('NDROPSRC ','#/kg/s  ',pver,'A', &
         'Droplet number source',phys_decomp)
    call addfld('NDROPSNK ','#/kg/s  ',pver,'A', &
         'Droplet number loss by microphysics',phys_decomp)
    call addfld('NDROPCOL ','#/m2    ',1,'A', &
         'Column droplet number',phys_decomp)

    if (history_amwg) call add_default('CCN3', 1, ' ')

    call register_dropmixnuc_host_hooks(cam_timer_start, cam_timer_stop, &
         cam_outfld_real1d, cam_outfld_real2d, cam_endrun)
    call dropmixnuc_init(pcols, pver, pverp, top_lev, iulog, pi, rhoh2o, &
         mwh2o, r_universal, rh2o, gravit, latvap, cpair, rair, &
         nspec_amode, sigmag_amode, dgnumlo_amode, dgnumhi_amode, &
         mam_idx, mam_cnst_idx, specdens_amode, spechygro_amode, &
         prog_modal_aero, fieldname, fieldname_cw)
  end subroutine ndrop_init

  subroutine dropmixnuc(state, ptend, dtmicro, pbuf, wsub, cldn, cldo, &
       tendnd, factnum)
    use physics_buffer, only: pbuf_get_field
    use physics_types, only: physics_ptend_init
    use ppgrid, only: pcols, pver
    use rad_constituents, only: rad_cnst_get_mode_num, rad_cnst_get_aer_mmr
    use ref_pres, only: top_lev => trop_cloud_top_lev
    use physconst, only: rair
    use wv_saturation, only: qsat

    type(physics_state), target, intent(in) :: state
    type(physics_ptend), intent(out) :: ptend
    real(r8), intent(in) :: dtmicro
    type(physics_buffer_desc), pointer, intent(inout) :: pbuf(:)
    real(r8), intent(in) :: wsub(:,:), cldn(:,:), cldo(:,:)
    real(r8), intent(out) :: tendnd(:,:), factnum(:,:,:)

    type(ptr2d_t), allocatable :: qqcw_ptr(:)
    real(r8), allocatable :: raer_data(:,:,:), qqcw_data(:,:,:)
    real(r8), allocatable :: qsat_es(:,:), qsat_qs(:,:)
    real(r8), pointer :: field(:,:), kvh(:,:)
    real(r8) :: pres, rhoair
    integer :: i, k, l, m, mm

    allocate(raer_data(pcols,pver,ncnst_tot), &
         qqcw_data(pcols,pver,ncnst_tot), qqcw_ptr(ncnst_tot), &
         qsat_es(pcols,pver), qsat_qs(pcols,pver))

    do m = 1, ntot_amode
       mm = mam_idx(m,0)
       call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, field)
       raer_data(:,:,mm) = field
       call rad_cnst_get_mode_num(0, m, 'c', state, pbuf, qqcw_ptr(mm)%fld)
       qqcw_data(:,:,mm) = qqcw_ptr(mm)%fld
       do l = 1, nspec_amode(m)
          mm = mam_idx(m,l)
          call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, field)
          raer_data(:,:,mm) = field
          call rad_cnst_get_aer_mmr(0, m, l, 'c', state, pbuf, qqcw_ptr(mm)%fld)
          qqcw_data(:,:,mm) = qqcw_ptr(mm)%fld
       end do
    end do

    qsat_es = 0._r8
    qsat_qs = 0._r8
    do k = top_lev, pver
       do i = 1, state%ncol
          rhoair = state%pmid(i,k)/(rair*state%t(i,k))
          pres = rair*rhoair*state%t(i,k)
          call qsat(state%t(i,k), pres, qsat_es(i,k), qsat_qs(i,k))
       end do
    end do

    call pbuf_get_field(pbuf, kvh_idx, kvh)
    if (prog_modal_aero) then
       call physics_ptend_init(ptend, state%psetcols, 'ndrop', lq=lq)
    else
       call physics_ptend_init(ptend, state%psetcols, 'ndrop')
    end if

    call dropmixnuc_run(state%lchnk, state%ncol, dtmicro, &
         state%q(:,:,numliq_idx), state%t, state%pmid, state%pint, &
         state%pdel, state%rpdel, state%zm, kvh, qsat_es, qsat_qs, &
         raer_data, qqcw_data, ptend%q, wsub, cldn, cldo, tendnd, factnum)

    if (prog_modal_aero) then
       do m = 1, ntot_amode
          do l = 0, nspec_amode(m)
             mm = mam_idx(m,l)
             qqcw_ptr(mm)%fld = qqcw_data(:,:,mm)
          end do
       end do
    end if

    deallocate(raer_data, qqcw_data, qqcw_ptr, qsat_es, qsat_qs)
  end subroutine dropmixnuc

  subroutine cam_timer_start(timer_name)
    use perf_mod, only: t_startf
    character(len=*), intent(in) :: timer_name
    call t_startf(timer_name)
  end subroutine cam_timer_start

  subroutine cam_timer_stop(timer_name)
    use perf_mod, only: t_stopf
    character(len=*), intent(in) :: timer_name
    call t_stopf(timer_name)
  end subroutine cam_timer_stop

  subroutine cam_outfld_real1d(field_name, field, dim1, lchnk)
    use cam_history, only: outfld
    character(len=*), intent(in) :: field_name
    real(r8), intent(in) :: field(:)
    integer, intent(in) :: dim1, lchnk
    call outfld(field_name, field, dim1, lchnk)
  end subroutine cam_outfld_real1d

  subroutine cam_outfld_real2d(field_name, field, dim1, lchnk)
    use cam_history, only: outfld
    character(len=*), intent(in) :: field_name
    real(r8), intent(in) :: field(:,:)
    integer, intent(in) :: dim1, lchnk
    call outfld(field_name, field, dim1, lchnk)
  end subroutine cam_outfld_real2d

  subroutine cam_endrun(message)
    use cam_abortutils, only: endrun
    character(len=*), intent(in) :: message
    call endrun(message)
  end subroutine cam_endrun

end module ndrop
