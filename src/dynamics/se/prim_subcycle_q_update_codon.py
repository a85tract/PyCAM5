import se_dynamics_prim_subcycle_codon as _prim
import se_dynamics_euler_hypervis_codon as _euler
import se_dynamics_sphere_ops_codon as _sphere
import se_dynamics_vertical_remap_codon as _vertical
import se_dynamics_remap_q_ppm_codon as _remap
import se_dynamics_misc_codon as _misc

@export
def get_block_lvl_cnt_d_codon(plevp: int) -> int:
    return plevp

@export
def get_gcol_block_cnt_d_codon() -> int:
    return 1

@export
def set_interp_parameter_codon(
    parm_code: int,
    value: int,
    gridtype_in: int,
    itype_in: int,
    nlon_in: int,
    nlat_in: int,
    auto_grid_in: int,
    itype_out_p: cobj,
    nlon_out_p: cobj,
    nlat_out_p: cobj,
    gridtype_out_p: cobj,
    auto_grid_out_p: cobj,
) -> int:
    return _misc.set_interp_parameter_codon(
        parm_code,
        value,
        gridtype_in,
        itype_in,
        nlon_in,
        nlat_in,
        auto_grid_in,
        itype_out_p,
        nlon_out_p,
        nlat_out_p,
        gridtype_out_p,
        auto_grid_out_p,
    )

@export
def prim_subcycle_dp3d_init_codon(
    np: int,
    nlev: int,
    ps0: float,
    hyai_p: cobj,
    hybi_p: cobj,
    ps_v_p: cobj,
    dp3d_p: cobj,
):
    return _prim.prim_subcycle_dp3d_init_codon(
        np,
        nlev,
        ps0,
        hyai_p,
        hybi_p,
        ps_v_p,
        dp3d_p,
    )

@export
def prim_subcycle_q_update_codon(
    np: int,
    nlev: int,
    qsize: int,
    ps0: float,
    hyai_p: cobj,
    hybi_p: cobj,
    ps_v_p: cobj,
    dp_np1_p: cobj,
    qdp_p: cobj,
    q_p: cobj,
):
    return _prim.prim_subcycle_q_update_codon(
        np,
        nlev,
        qsize,
        ps0,
        hyai_p,
        hybi_p,
        ps_v_p,
        dp_np1_p,
        qdp_p,
        q_p,
    )

@export
def se_log2_codon(n: int) -> int:
    return _misc.se_log2_codon(n)

@export
def se_factor_fill_codon(num: int, factors_p: cobj, numfact_p: cobj):
    return _misc.se_factor_fill_codon(num, factors_p, numfact_p)

@export
def se_calcsegmentlength_codon(lenp: int, lens: int, mpattern: int, nlyr: int, hme_mpattern_s: int, hme_mpattern_p: int) -> int:
    return _misc.se_calcsegmentlength_codon(lenp, lens, mpattern, nlyr, hme_mpattern_s, hme_mpattern_p)

@export
def se_timelevel_init_default_codon(nm1_p: cobj, n0_p: cobj, np1_p: cobj, nstep_p: cobj, nstep0_p: cobj):
    return _misc.se_timelevel_init_default_codon(nm1_p, n0_p, np1_p, nstep_p, nstep0_p)

@export
def se_timelevel_update_codon(nm1_p: cobj, n0_p: cobj, np1_p: cobj, nstep_p: cobj, uptype_code: int) -> int:
    return _misc.se_timelevel_update_codon(nm1_p, n0_p, np1_p, nstep_p, uptype_code)

@export
def qdp_time_avg_codon(
    np: int,
    nlev: int,
    qsize: int,
    rkstage: int,
    n0_qdp: int,
    np1_qdp: int,
    qdp_p: cobj,
):
    return _prim.qdp_time_avg_codon(
        np,
        nlev,
        qsize,
        rkstage,
        n0_qdp,
        np1_qdp,
        qdp_p,
    )

@export
def euler_step_vstar_prepare_codon(
    np: int,
    nlev: int,
    dt: float,
    rhs_multiplier: int,
    dp_in_p: cobj,
    divdp_proj_p: cobj,
    vn0_p: cobj,
    dp_out_p: cobj,
    vstar_p: cobj,
):
    return _euler.euler_step_vstar_prepare_codon(
        np,
        nlev,
        dt,
        rhs_multiplier,
        dp_in_p,
        divdp_proj_p,
        vn0_p,
        dp_out_p,
        vstar_p,
    )

@export
def euler_step_limiter_dpstar_codon(
    np: int,
    nlev: int,
    dt: float,
    rhs_viss: int,
    nu_q: float,
    nu_p: float,
    dp_in_p: cobj,
    divdp_p: cobj,
    dpdiss_biharmonic_p: cobj,
    spheremp_p: cobj,
    dp_star_p: cobj,
):
    return _euler.euler_step_limiter_dpstar_codon(
        np,
        nlev,
        dt,
        rhs_viss,
        nu_q,
        nu_p,
        dp_in_p,
        divdp_p,
        dpdiss_biharmonic_p,
        spheremp_p,
        dp_star_p,
    )

@export
def euler_step_qdp_writeback_codon(
    np: int,
    nlev: int,
    qsize: int,
    qidx: int,
    np1_qdp: int,
    qdp_p: cobj,
    spheremp_p: cobj,
    qtens_p: cobj,
):
    return _euler.euler_step_qdp_writeback_codon(
        np,
        nlev,
        qsize,
        qidx,
        np1_qdp,
        qdp_p,
        spheremp_p,
        qtens_p,
    )

@export
def euler_step_qdp_restore_codon(
    np: int,
    nlev: int,
    qdp_p: cobj,
    rspheremp_p: cobj,
):
    return _euler.euler_step_qdp_restore_codon(
        np,
        nlev,
        qdp_p,
        rspheremp_p,
    )

@export
def advance_hypervis_qtens_prepare_codon(
    np: int,
    nlev: int,
    qsize: int,
    ps0: float,
    dt2: float,
    nu_p: float,
    hyai_p: cobj,
    hybi_p: cobj,
    dp_in_p: cobj,
    divdp_proj_p: cobj,
    dpdiss_ave_p: cobj,
    qdp_p: cobj,
    dp_out_p: cobj,
    qtens_p: cobj,
):
    return _euler.advance_hypervis_qtens_prepare_codon(
        np,
        nlev,
        qsize,
        ps0,
        dt2,
        nu_p,
        hyai_p,
        hybi_p,
        dp_in_p,
        divdp_proj_p,
        dpdiss_ave_p,
        qdp_p,
        dp_out_p,
        qtens_p,
    )

@export
def advance_hypervis_qdp_update_codon(
    np: int,
    nlev: int,
    dt: float,
    nu_q: float,
    qdp_p: cobj,
    spheremp_p: cobj,
    qtens_p: cobj,
):
    return _euler.advance_hypervis_qdp_update_codon(
        np,
        nlev,
        dt,
        nu_q,
        qdp_p,
        spheremp_p,
        qtens_p,
    )

@export
def euler_step_dssvar_restore_codon(
    np: int,
    nlev: int,
    dssvar_p: cobj,
    rspheremp_p: cobj,
):
    return _euler.euler_step_dssvar_restore_codon(
        np,
        nlev,
        dssvar_p,
        rspheremp_p,
    )

@export
def euler_step_dssvar_pack_codon(
    np: int,
    nlev: int,
    dssvar_p: cobj,
    spheremp_p: cobj,
):
    return _euler.euler_step_dssvar_pack_codon(
        np,
        nlev,
        dssvar_p,
        spheremp_p,
    )

@export
def euler_step_qtens_base_codon(
    np: int,
    dt: float,
    qdp_p: cobj,
    dp_star_p: cobj,
    qtens_p: cobj,
):
    return _euler.euler_step_qtens_base_codon(
        np,
        dt,
        qdp_p,
        dp_star_p,
        qtens_p,
    )

@export
def euler_step_gradq_prepare_codon(
    np: int,
    vstar1_p: cobj,
    vstar2_p: cobj,
    qdp_p: cobj,
    gradq1_p: cobj,
    gradq2_p: cobj,
):
    return _euler.euler_step_gradq_prepare_codon(
        np,
        vstar1_p,
        vstar2_p,
        qdp_p,
        gradq1_p,
        gradq2_p,
    )

@export
def euler_step_qtens_biharmonic_add_codon(
    np: int,
    qtens_p: cobj,
    qtens_biharmonic_p: cobj,
):
    return _euler.euler_step_qtens_biharmonic_add_codon(
        np,
        qtens_p,
        qtens_biharmonic_p,
    )

@export
def euler_step_qtens_biharmonic_init_codon(
    np: int,
    nlev: int,
    qsize: int,
    dt: float,
    rhs_multiplier: int,
    dp_in_p: cobj,
    divdp_proj_p: cobj,
    qdp_p: cobj,
    dp_out_p: cobj,
    qtens_biharmonic_p: cobj,
):
    return _euler.euler_step_qtens_biharmonic_init_codon(
        np,
        nlev,
        qsize,
        dt,
        rhs_multiplier,
        dp_in_p,
        divdp_proj_p,
        qdp_p,
        dp_out_p,
        qtens_biharmonic_p,
    )

@export
def euler_step_qtens_biharmonic_scale_codon(
    np: int,
    nlev: int,
    qsize: int,
    qtens_biharmonic_p: cobj,
    dpdiss_ave_p: cobj,
    dp0_p: cobj,
):
    return _euler.euler_step_qtens_biharmonic_scale_codon(
        np,
        nlev,
        qsize,
        qtens_biharmonic_p,
        dpdiss_ave_p,
        dp0_p,
    )

@export
def euler_step_qtens_biharmonic_unapply_codon(
    np: int,
    nlev: int,
    rhs_viss: int,
    dt: float,
    nu_q: float,
    qtens_biharmonic_p: cobj,
    spheremp_p: cobj,
    dp0_p: cobj,
):
    return _euler.euler_step_qtens_biharmonic_unapply_codon(
        np,
        nlev,
        rhs_viss,
        dt,
        nu_q,
        qtens_biharmonic_p,
        spheremp_p,
        dp0_p,
    )

@export
def euler_step_qminmax_update_codon(
    np: int,
    nlev: int,
    qsize: int,
    rhs_multiplier: int,
    qtens_biharmonic_p: cobj,
    qmin_p: cobj,
    qmax_p: cobj,
):
    return _euler.euler_step_qminmax_update_codon(
        np,
        nlev,
        qsize,
        rhs_multiplier,
        qtens_biharmonic_p,
        qmin_p,
        qmax_p,
    )

@export
def limiter2d_zero_codon(
    np: int,
    nlev: int,
    q_p: cobj,
):
    return _euler.limiter2d_zero_codon(
        np,
        nlev,
        q_p,
    )

@export
def limiter_optim_iter_full_codon(
    np: int,
    nlev: int,
    ptens_p: cobj,
    sphweights_p: cobj,
    minp_p: cobj,
    maxp_p: cobj,
    dpmass_p: cobj,
    weights_p: cobj,
    whois_neg_p: cobj,
    whois_pos_p: cobj,
    x_p: cobj,
    c_p: cobj,
    al_neg_p: cobj,
    al_pos_p: cobj,
):
    return _euler.limiter_optim_iter_full_codon(
        np,
        nlev,
        ptens_p,
        sphweights_p,
        minp_p,
        maxp_p,
        dpmass_p,
        weights_p,
        whois_neg_p,
        whois_pos_p,
        x_p,
        c_p,
        al_neg_p,
        al_pos_p,
    )

@export
def divergence_sphere_codon(
    np: int,
    rrearth: float,
    v_p: cobj,
    dvv_p: cobj,
    metdet_p: cobj,
    dinv_p: cobj,
    rmetdet_p: cobj,
    gv_p: cobj,
    vvtemp_p: cobj,
    div_p: cobj,
):
    return _sphere.divergence_sphere_codon(
        np,
        rrearth,
        v_p,
        dvv_p,
        metdet_p,
        dinv_p,
        rmetdet_p,
        gv_p,
        vvtemp_p,
        div_p,
    )

@export
def divergence_sphere_wk_codon(
    np: int,
    rrearth: float,
    v_p: cobj,
    dvv_p: cobj,
    spheremp_p: cobj,
    dinv_p: cobj,
    vtemp_p: cobj,
    div_p: cobj,
):
    return _sphere.divergence_sphere_wk_codon(
        np,
        rrearth,
        v_p,
        dvv_p,
        spheremp_p,
        dinv_p,
        vtemp_p,
        div_p,
    )

@export
def laplace_sphere_wk_codon(
    np: int,
    rrearth: float,
    hypervis_power: int,
    hypervis_scaling: int,
    var_coef: int,
    s_p: cobj,
    dvv_p: cobj,
    spheremp_p: cobj,
    dinv_p: cobj,
    variable_hyperviscosity_p: cobj,
    tensorvisc_p: cobj,
    grads_p: cobj,
    oldgrads_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    laplace_p: cobj,
):
    return _sphere.laplace_sphere_wk_codon(
        np,
        rrearth,
        hypervis_power,
        hypervis_scaling,
        var_coef,
        s_p,
        dvv_p,
        spheremp_p,
        dinv_p,
        variable_hyperviscosity_p,
        tensorvisc_p,
        grads_p,
        oldgrads_p,
        v1_p,
        v2_p,
        laplace_p,
    )

@export
def vlaplace_sphere_wk_codon(
    np: int,
    rrearth: float,
    hypervis_power: int,
    hypervis_scaling: int,
    var_coef: int,
    has_nu_ratio: int,
    nu_ratio: float,
    v_p: cobj,
    dvv_p: cobj,
    mp_p: cobj,
    spheremp_p: cobj,
    metinv_p: cobj,
    metdet_p: cobj,
    rmetdet_p: cobj,
    d_p: cobj,
    dinv_p: cobj,
    variable_hyperviscosity_p: cobj,
    tensorvisc_p: cobj,
    vec_sphere2cart_p: cobj,
    dum_cart_p: cobj,
    dum_tmp_p: cobj,
    div_p: cobj,
    vor_p: cobj,
    lap_tmp_p: cobj,
    lap_tmp2_p: cobj,
    work1_p: cobj,
    work2_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    laplace_p: cobj,
):
    return _sphere.vlaplace_sphere_wk_codon(
        np,
        rrearth,
        hypervis_power,
        hypervis_scaling,
        var_coef,
        has_nu_ratio,
        nu_ratio,
        v_p,
        dvv_p,
        mp_p,
        spheremp_p,
        metinv_p,
        metdet_p,
        rmetdet_p,
        d_p,
        dinv_p,
        variable_hyperviscosity_p,
        tensorvisc_p,
        vec_sphere2cart_p,
        dum_cart_p,
        dum_tmp_p,
        div_p,
        vor_p,
        lap_tmp_p,
        lap_tmp2_p,
        work1_p,
        work2_p,
        v1_p,
        v2_p,
        laplace_p,
    )

@export
def vorticity_sphere_codon(
    np: int,
    rrearth: float,
    v_p: cobj,
    dvv_p: cobj,
    d_p: cobj,
    rmetdet_p: cobj,
    vco_p: cobj,
    vtemp_p: cobj,
    vort_p: cobj,
):
    return _sphere.vorticity_sphere_codon(
        np,
        rrearth,
        v_p,
        dvv_p,
        d_p,
        rmetdet_p,
        vco_p,
        vtemp_p,
        vort_p,
    )

@export
def gradient_sphere_codon(
    np: int,
    rrearth: float,
    s_p: cobj,
    dvv_p: cobj,
    dinv_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    ds_p: cobj,
):
    return _sphere.gradient_sphere_codon(
        np,
        rrearth,
        s_p,
        dvv_p,
        dinv_p,
        v1_p,
        v2_p,
        ds_p,
    )

@export
def curl_sphere_codon(
    np: int,
    rrearth: float,
    s_p: cobj,
    dvv_p: cobj,
    d_p: cobj,
    metdet_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    ds_p: cobj,
):
    return _sphere.curl_sphere_codon(
        np,
        rrearth,
        s_p,
        dvv_p,
        d_p,
        metdet_p,
        v1_p,
        v2_p,
        ds_p,
    )

@export
def ugradv_sphere_codon(
    np: int,
    rrearth: float,
    u_p: cobj,
    v_p: cobj,
    dvv_p: cobj,
    dinv_p: cobj,
    vec_sphere2cart_p: cobj,
    dum_cart_p: cobj,
    tmp_p: cobj,
    v1_p: cobj,
    v2_p: cobj,
    ugradv_p: cobj,
):
    return _sphere.ugradv_sphere_codon(
        np,
        rrearth,
        u_p,
        v_p,
        dvv_p,
        dinv_p,
        vec_sphere2cart_p,
        dum_cart_p,
        tmp_p,
        v1_p,
        v2_p,
        ugradv_p,
    )

@export
def vertical_remap_rsplit_prepare_codon(
    np: int,
    nlev: int,
    ps0: float,
    hyai_p: cobj,
    hybi_p: cobj,
    ps_v_p: cobj,
    dp3d_p: cobj,
    dp_p: cobj,
    dp_star_p: cobj,
):
    return _vertical.vertical_remap_rsplit_prepare_codon(
        np,
        nlev,
        ps0,
        hyai_p,
        hybi_p,
        ps_v_p,
        dp3d_p,
        dp_p,
        dp_star_p,
    )

@export
def vertical_remap_t_scale_codon(
    np: int,
    nlev: int,
    t_p: cobj,
    dp_star_p: cobj,
    ttmp_p: cobj,
):
    return _vertical.vertical_remap_t_scale_codon(
        np,
        nlev,
        t_p,
        dp_star_p,
        ttmp_p,
    )

@export
def vertical_remap_t_unscale_codon(
    np: int,
    nlev: int,
    ttmp_p: cobj,
    dp_p: cobj,
    t_p: cobj,
):
    return _vertical.vertical_remap_t_unscale_codon(
        np,
        nlev,
        ttmp_p,
        dp_p,
        t_p,
    )

@export
def vertical_remap_v_scale_codon(
    np: int,
    nlev: int,
    v_p: cobj,
    dp_star_p: cobj,
    ttmp_p: cobj,
):
    return _vertical.vertical_remap_v_scale_codon(
        np,
        nlev,
        v_p,
        dp_star_p,
        ttmp_p,
    )

@export
def vertical_remap_v_unscale_codon(
    np: int,
    nlev: int,
    ttmp_p: cobj,
    dp_p: cobj,
    v_p: cobj,
):
    return _vertical.vertical_remap_v_unscale_codon(
        np,
        nlev,
        ttmp_p,
        dp_p,
        v_p,
    )

@export
def vertical_remap_ps_v_update_codon(
    np: int,
    nlev: int,
    hyai1: float,
    ps0: float,
    dp3d_p: cobj,
    ps_v_p: cobj,
):
    return _vertical.vertical_remap_ps_v_update_codon(
        np,
        nlev,
        hyai1,
        ps0,
        dp3d_p,
        ps_v_p,
    )

@export
def remap_q_ppm_interval_setup_codon(
    nlev: int,
    pio_p: cobj,
    pin_p: cobj,
    dpo_p: cobj,
    kid_p: cobj,
    z1_p: cobj,
    z2_p: cobj,
):
    return _remap.remap_q_ppm_interval_setup_codon(
        nlev,
        pio_p,
        pin_p,
        dpo_p,
        kid_p,
        z1_p,
        z2_p,
    )

@export
def remap_q_ppm_mass_prep_codon(
    nx: int,
    nlev: int,
    qsize: int,
    iidx: int,
    jidx: int,
    qidx: int,
    qdp_p: cobj,
    dpo_p: cobj,
    masso_p: cobj,
    ao_p: cobj,
):
    return _remap.remap_q_ppm_mass_prep_codon(
        nx,
        nlev,
        qsize,
        iidx,
        jidx,
        qidx,
        qdp_p,
        dpo_p,
        masso_p,
        ao_p,
    )

@export
def remap_q_ppm_compute_ppm_grids_codon(
    nlev: int,
    vert_remap_q_alg: int,
    dx_p: cobj,
    rslt_p: cobj,
):
    return _remap.remap_q_ppm_compute_ppm_grids_codon(
        nlev,
        vert_remap_q_alg,
        dx_p,
        rslt_p,
    )

@export
def remap_q_ppm_compute_ppm_codon(
    nlev: int,
    vert_remap_q_alg: int,
    a_p: cobj,
    dx_p: cobj,
    ai_p: cobj,
    dma_p: cobj,
    coefs_p: cobj,
):
    return _remap.remap_q_ppm_compute_ppm_codon(
        nlev,
        vert_remap_q_alg,
        a_p,
        dx_p,
        ai_p,
        dma_p,
        coefs_p,
    )

@export
def remap_q_ppm_mass_apply_codon(
    nx: int,
    nlev: int,
    iidx: int,
    jidx: int,
    qidx: int,
    kid_p: cobj,
    masso_p: cobj,
    coefs_p: cobj,
    z1_p: cobj,
    z2_p: cobj,
    dpo_p: cobj,
    qdp_p: cobj,
):
    return _remap.remap_q_ppm_mass_apply_codon(
        nx,
        nlev,
        iidx,
        jidx,
        qidx,
        kid_p,
        masso_p,
        coefs_p,
        z1_p,
        z2_p,
        dpo_p,
        qdp_p,
    )

@export
def se_misc_touch_codon(
    tag: int,
):
    return _misc.se_misc_touch_codon(
        tag,
    )

@export
def get_block_gcol_d_codon(
    size: int,
    unique_pt_offset: int,
    cdex_p: cobj,
):
    return _misc.get_block_gcol_d_codon(size, unique_pt_offset, cdex_p)

@export
def get_block_owner_d_codon(
    owner: int,
) -> int:
    return _misc.get_block_owner_d_codon(owner)

@export
def get_dyn_grid_parm_real2d_codon(
    name_code: int,
) -> int:
    return _misc.get_dyn_grid_parm_real2d_codon(name_code)

@export
def get_dyn_grid_parm_real1d_codon(
    name_code: int,
) -> int:
    return _misc.get_dyn_grid_parm_real1d_codon(name_code)

@export
def get_dyn_grid_parm_codon(
    name_code: int,
    ne: int,
    np: int,
    npsq: int,
    nelemd: int,
    beglat: int,
    endlat: int,
    ngcols_d: int,
    plat: int,
    plev: int,
    plevp: int,
    nlon: int,
    nlat: int,
) -> int:
    return _misc.get_dyn_grid_parm_codon(
        name_code,
        ne,
        np,
        npsq,
        nelemd,
        beglat,
        endlat,
        ngcols_d,
        plat,
        plev,
        plevp,
        nlon,
        nlat,
    )

@export
def get_ldof_fill_codon(
    nlev: int,
    nelemd: int,
    hdim: int,
    num_unique_pts_p: cobj,
    unique_pt_offsets_p: cobj,
    ldof_p: cobj,
) -> int:
    return _misc.get_ldof_fill_codon(
        nlev,
        nelemd,
        hdim,
        num_unique_pts_p,
        unique_pt_offsets_p,
        ldof_p,
    )

@export
def latlon_interpolation_codon(
    t: int,
    n: int,
    value: int,
) -> int:
    return _misc.latlon_interpolation_codon(t, n, value)

@export
def dycore_is_codon(
    is_match: int,
) -> int:
    return _misc.dycore_is_codon(is_match)

@export
def isfactorable_codon(
    n: int,
) -> int:
    return _misc.isfactorable_codon(n)

@export
def genlocaldof_codon(
    ig: int,
    npts: int,
    ldof_p: cobj,
):
    return _misc.genlocaldof_codon(ig, npts, ldof_p)

@export
def uniquepoints2d_codon(
    num_unique_pts: int,
    ia_p: cobj,
    ja_p: cobj,
    ni: int,
    src_p: cobj,
    dest_p: cobj,
):
    return _misc.uniquepoints2d_codon(num_unique_pts, ia_p, ja_p, ni, src_p, dest_p)

@export
def convert_gbl_index_codon(
    number: int,
    ne: int,
    ie_p: cobj,
    je_p: cobj,
    face_no_p: cobj,
):
    return _misc.convert_gbl_index_codon(number, ne, ie_p, je_p, face_no_p)

@export
def gridedge_type_codon(
    head_processor: int,
    tail_processor: int,
    internal_edge: int,
    external_edge: int,
) -> int:
    return _misc.gridedge_type_codon(head_processor, tail_processor, internal_edge, external_edge)

@export
def copy_buffer_codon(
    nthreads: int,
    ithr: int,
    len_move_ptr: int,
    buf_p: cobj,
    receive_p: cobj,
    move_ptr_p: cobj,
    move_length_p: cobj,
):
    return _misc.copy_buffer_codon(
        nthreads,
        ithr,
        len_move_ptr,
        buf_p,
        receive_p,
        move_ptr_p,
        move_length_p,
    )

@export
def var_is_vector_codon(
    name_len: int,
    name_ascii_p: cobj,
    entry_len: int,
    entries_ascii_p: cobj,
    nentries: int,
) -> int:
    return _misc.var_is_vector_codon(
        name_len,
        name_ascii_p,
        entry_len,
        entries_ascii_p,
        nentries,
    )

@export
def reduction_max_r_local_codon(
    buf_p: cobj,
    ctr_p: cobj,
    redp_p: cobj,
    length: int,
    nthreads: int,
):
    return _misc.reduction_max_r_local_codon(
        buf_p,
        ctr_p,
        redp_p,
        length,
        nthreads,
    )

@export
def reduction_min_r_local_codon(
    buf_p: cobj,
    ctr_p: cobj,
    redp_p: cobj,
    length: int,
    nthreads: int,
):
    return _misc.reduction_min_r_local_codon(
        buf_p,
        ctr_p,
        redp_p,
        length,
        nthreads,
    )

@export
def copy_par_codon(
    rank2_p: cobj,
    root2_p: cobj,
    nprocs2_p: cobj,
    comm2_p: cobj,
    intercomm2_p: cobj,
    intracomm2_p: cobj,
    intracommsize2_p: cobj,
    intracommrank2_p: cobj,
    comm_graph_full2_p: cobj,
    comm_graph_inter2_p: cobj,
    comm_graph_intra2_p: cobj,
    group_graph_full2_p: cobj,
    masterproc2_p: cobj,
    rank1_p: cobj,
    root1_p: cobj,
    nprocs1_p: cobj,
    comm1_p: cobj,
    intercomm1_p: cobj,
    intracomm1_p: cobj,
    intracommsize1_p: cobj,
    intracommrank1_p: cobj,
    comm_graph_full1_p: cobj,
    comm_graph_inter1_p: cobj,
    comm_graph_intra1_p: cobj,
    group_graph_full1_p: cobj,
    masterproc1_p: cobj,
):
    return _misc.copy_par_codon(
        rank2_p,
        root2_p,
        nprocs2_p,
        comm2_p,
        intercomm2_p,
        intracomm2_p,
        intracommsize2_p,
        intracommrank2_p,
        comm_graph_full2_p,
        comm_graph_inter2_p,
        comm_graph_intra2_p,
        group_graph_full2_p,
        masterproc2_p,
        rank1_p,
        root1_p,
        nprocs1_p,
        comm1_p,
        intercomm1_p,
        intracomm1_p,
        intracommsize1_p,
        intracommrank1_p,
        comm_graph_full1_p,
        comm_graph_inter1_p,
        comm_graph_intra1_p,
        group_graph_full1_p,
        masterproc1_p,
    )

@export
def init_edge_buffer_i8_header_codon(
    np: int,
    max_corner_elem: int,
    nelemd: int,
    nlyr: int,
    nlyr_p: cobj,
    nbuf_p: cobj,
):
    return _misc.init_edge_buffer_i8_header_codon(
        np,
        max_corner_elem,
        nelemd,
        nlyr,
        nlyr_p,
        nbuf_p,
    )

@export
def zero_i32_buffer_codon(n: int, buf_p: cobj):
    return _misc.zero_i32_buffer_codon(n, buf_p)

@export
def projectpoint_codon(
    cart_x: float,
    cart_y: float,
    face_no: int,
    r_p: cobj,
    lon_p: cobj,
    lat_p: cobj,
):
    return _misc.projectpoint_codon(cart_x, cart_y, face_no, r_p, lon_p, lat_p)

@export
def ref2sphere_double_codon(
    a: float,
    b: float,
    face_no: int,
    c1x: float,
    c1y: float,
    c2x: float,
    c2y: float,
    c3x: float,
    c3y: float,
    c4x: float,
    c4y: float,
    r_p: cobj,
    lon_p: cobj,
    lat_p: cobj,
):
    return _misc.ref2sphere_double_codon(
        a, b, face_no, c1x, c1y, c2x, c2y, c3x, c3y, c4x, c4y, r_p, lon_p, lat_p
    )

@export
def dmap_equiangular_codon(
    a: float,
    b: float,
    face_no: int,
    c1x: float,
    c1y: float,
    c2x: float,
    c2y: float,
    c3x: float,
    c3y: float,
    c4x: float,
    c4y: float,
    u11: float,
    u12: float,
    u21: float,
    u22: float,
    u31: float,
    u32: float,
    u41: float,
    u42: float,
    d_p: cobj,
):
    return _misc.dmap_equiangular_codon(
        a, b, face_no, c1x, c1y, c2x, c2y, c3x, c3y, c4x, c4y,
        u11, u12, u21, u22, u31, u32, u41, u42, d_p
    )

@export
def create_work_pool_codon(
    start_domain: int,
    end_domain: int,
    ndomains: int,
    ipe: int,
    beg_index_p: cobj,
    end_index_p: cobj,
):
    return _misc.create_work_pool_codon(
        start_domain,
        end_domain,
        ndomains,
        ipe,
        beg_index_p,
        end_index_p,
    )

@export
def set_thread_ranges_1d_codon(
    work_pool_p: cobj,
    nrows: int,
    idthread: int,
    beg_range_p: cobj,
    end_range_p: cobj,
):
    return _misc.set_thread_ranges_1d_codon(
        work_pool_p,
        nrows,
        idthread,
        beg_range_p,
        end_range_p,
    )

@export
def get_loop_ranges_codon(
    ibeg_in: int,
    iend_in: int,
    kbeg_in: int,
    kend_in: int,
    qbeg_in: int,
    qend_in: int,
    mask: int,
    ibeg_p: cobj,
    iend_p: cobj,
    kbeg_p: cobj,
    kend_p: cobj,
    qbeg_p: cobj,
    qend_p: cobj,
):
    return _misc.get_loop_ranges_codon(
        ibeg_in,
        iend_in,
        kbeg_in,
        kend_in,
        qbeg_in,
        qend_in,
        mask,
        ibeg_p,
        iend_p,
        kbeg_p,
        kend_p,
        qbeg_p,
        qend_p,
    )

@export
def timelevel_qdp_codon(
    nstep: int,
    qsplit: int,
    has_np1: int,
    n0_p: cobj,
    np1_p: cobj,
):
    return _misc.timelevel_qdp_codon(
        nstep,
        qsplit,
        has_np1,
        n0_p,
        np1_p,
    )

@export
def elem_jacobians_codon(
    coords_xy_p: cobj,
    unif2quadmap_p: cobj,
):
    return _misc.elem_jacobians_codon(
        coords_xy_p,
        unif2quadmap_p,
    )

@export
def element_var_coordinates_codon(
    npts: int,
    corners_xy_p: cobj,
    points_p: cobj,
    cart_xy_p: cobj,
):
    return _misc.element_var_coordinates_codon(
        npts,
        corners_xy_p,
        points_p,
        cart_xy_p,
    )

@export
def gausslobatto_wts_codon(
    np1: int,
    glpts_p: cobj,
    wts_p: cobj,
):
    return _misc.gausslobatto_wts_codon(
        np1,
        glpts_p,
        wts_p,
    )

@export
def find_buffer_slot_codon(
    inbr: int,
    length: int,
    tmp_p: cobj,
    n: int,
    ptr_p: cobj,
):
    return _misc.find_buffer_slot_codon(
        inbr,
        length,
        tmp_p,
        n,
        ptr_p,
    )

@export
def cubesetupedgeindex_codon(
    s_face: int,
    d_face: int,
    south: int,
    east: int,
    north: int,
    west: int,
    reverse_p: cobj,
):
    return _misc.cubesetupedgeindex_codon(
        s_face,
        d_face,
        south,
        east,
        north,
        west,
        reverse_p,
    )

@export
def copy_buffer_codon(
    nthreads: int,
    ithr: int,
    len_move_ptr: int,
    buf_p: cobj,
    receive_p: cobj,
    move_ptr_p: cobj,
    move_length_p: cobj,
):
    return _misc.copy_buffer_codon(
        nthreads,
        ithr,
        len_move_ptr,
        buf_p,
        receive_p,
        move_ptr_p,
        move_length_p,
    )

@export
def var_is_vector_codon(
    name_len: int,
    name_ascii_p: cobj,
    entry_len: int,
    entries_ascii_p: cobj,
    nentries: int,
) -> int:
    return _misc.var_is_vector_codon(
        name_len,
        name_ascii_p,
        entry_len,
        entries_ascii_p,
        nentries,
    )

@export
def reduction_max_r_local_codon(
    buf_p: cobj,
    ctr_p: cobj,
    redp_p: cobj,
    length: int,
    nthreads: int,
):
    return _misc.reduction_max_r_local_codon(
        buf_p,
        ctr_p,
        redp_p,
        length,
        nthreads,
    )

@export
def reduction_min_r_local_codon(
    buf_p: cobj,
    ctr_p: cobj,
    redp_p: cobj,
    length: int,
    nthreads: int,
):
    return _misc.reduction_min_r_local_codon(
        buf_p,
        ctr_p,
        redp_p,
        length,
        nthreads,
    )

@export
def copy_par_codon(
    rank2_p: cobj,
    root2_p: cobj,
    nprocs2_p: cobj,
    comm2_p: cobj,
    intercomm2_p: cobj,
    intracomm2_p: cobj,
    intracommsize2_p: cobj,
    intracommrank2_p: cobj,
    comm_graph_full2_p: cobj,
    comm_graph_inter2_p: cobj,
    comm_graph_intra2_p: cobj,
    group_graph_full2_p: cobj,
    masterproc2_p: cobj,
    rank1_p: cobj,
    root1_p: cobj,
    nprocs1_p: cobj,
    comm1_p: cobj,
    intercomm1_p: cobj,
    intracomm1_p: cobj,
    intracommsize1_p: cobj,
    intracommrank1_p: cobj,
    comm_graph_full1_p: cobj,
    comm_graph_inter1_p: cobj,
    comm_graph_intra1_p: cobj,
    group_graph_full1_p: cobj,
    masterproc1_p: cobj,
):
    return _misc.copy_par_codon(
        rank2_p,
        root2_p,
        nprocs2_p,
        comm2_p,
        intercomm2_p,
        intracomm2_p,
        intracommsize2_p,
        intracommrank2_p,
        comm_graph_full2_p,
        comm_graph_inter2_p,
        comm_graph_intra2_p,
        group_graph_full2_p,
        masterproc2_p,
        rank1_p,
        root1_p,
        nprocs1_p,
        comm1_p,
        intercomm1_p,
        intracomm1_p,
        intracommsize1_p,
        intracommrank1_p,
        comm_graph_full1_p,
        comm_graph_inter1_p,
        comm_graph_intra1_p,
        group_graph_full1_p,
        masterproc1_p,
    )

@export
def init_edge_buffer_i8_header_codon(
    np: int,
    max_corner_elem: int,
    nelemd: int,
    nlyr: int,
    nlyr_p: cobj,
    nbuf_p: cobj,
):
    return _misc.init_edge_buffer_i8_header_codon(
        np,
        max_corner_elem,
        nelemd,
        nlyr,
        nlyr_p,
        nbuf_p,
    )

@export
def zero_i32_buffer_codon(n: int, buf_p: cobj):
    return _misc.zero_i32_buffer_codon(n, buf_p)

@export
def gbarrier_init_codon(
    c_barrier_p: cobj,
    nthreads: int,
):
    return _misc.gbarrier_init_codon(
        c_barrier_p,
        nthreads,
    )

@export
def gbarrier_delete_codon(
    c_barrier_p: cobj,
):
    return _misc.gbarrier_delete_codon(
        c_barrier_p,
    )

@export
def gbarrier_synchronize_codon(
    c_barrier: cobj,
    thread: int,
):
    return _misc.gbarrier_synchronize_codon(
        c_barrier,
        thread,
    )

@export
def legendre_codon(
    x: float,
    n: int,
    leg_p: cobj,
):
    return _misc.legendre_codon(
        x,
        n,
        leg_p,
    )

@export
def se_gausslobatto_fill_codon(
    npts: int,
    points_p: cobj,
    weights_p: cobj,
) -> int:
    return _misc.se_gausslobatto_fill_codon(
        npts,
        points_p,
        weights_p,
    )

@export
def allocate_gridvertex_nbrs_select_dim_codon(
    has_dim: int,
    dim: int,
    default_dim: int,
) -> int:
    return _misc.allocate_gridvertex_nbrs_select_dim_codon(
        has_dim,
        dim,
        default_dim,
    )

@export
def deallocate_gridvertex_nbrs_touch_codon(
    tag: int,
) -> int:
    return _misc.deallocate_gridvertex_nbrs_touch_codon(tag)

@export
def applycamforcing_dynamics_codon(
    np: int,
    nlev: int,
    dt_q: float,
    t_p: cobj,
    ft_p: cobj,
    v_p: cobj,
    fm_p: cobj,
):
    return _misc.applycamforcing_dynamics_codon(
        np,
        nlev,
        dt_q,
        t_p,
        ft_p,
        v_p,
        fm_p,
    )

@export
def createuniqueindex_codon(
    ig: int,
    npts: int,
    gdof_p: cobj,
    ia_p: cobj,
    ja_p: cobj,
) -> int:
    return _misc.createuniqueindex_codon(
        ig,
        npts,
        gdof_p,
        ia_p,
        ja_p,
    )
