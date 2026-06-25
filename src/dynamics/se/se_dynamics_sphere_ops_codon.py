import se_dynamics_common_codon as _common

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
    v = Ptr[float](v_p)
    dvv = Ptr[float](dvv_p)
    metdet = Ptr[float](metdet_p)
    dinv = Ptr[float](dinv_p)
    rmetdet = Ptr[float](rmetdet_p)
    gv = Ptr[float](gv_p)
    vvtemp = Ptr[float](vvtemp_p)
    div = Ptr[float](div_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _common._plane_idx(i, j, np)
            v1 = v[_common._vec2_idx(i, j, 1, np)]
            v2 = v[_common._vec2_idx(i, j, 2, np)]
            gv[_common._vec2_idx(i, j, 1, np)] = metdet[plane_idx] * (
                dinv[_common._mat22_idx(i, j, 1, 1, np)] * v1 + dinv[_common._mat22_idx(i, j, 1, 2, np)] * v2
            )
            gv[_common._vec2_idx(i, j, 2, np)] = metdet[plane_idx] * (
                dinv[_common._mat22_idx(i, j, 2, 1, np)] * v1 + dinv[_common._mat22_idx(i, j, 2, 2, np)] * v2
            )

    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dudx00 = 0.0
            dvdy00 = 0.0
            for i in range(1, np + 1):
                dudx00 = dudx00 + dvv[_common._plane_idx(i, l, np)] * gv[_common._vec2_idx(i, j, 1, np)]
                dvdy00 = dvdy00 + dvv[_common._plane_idx(i, l, np)] * gv[_common._vec2_idx(j, i, 2, np)]
            div[_common._plane_idx(l, j, np)] = dudx00
            vvtemp[_common._plane_idx(j, l, np)] = dvdy00

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _common._plane_idx(i, j, np)
            div[plane_idx] = (div[plane_idx] + vvtemp[plane_idx]) * (rmetdet[plane_idx] * rrearth)


def _divergence_sphere_wk_vec(
    np: int,
    rrearth: float,
    v: Ptr[float],
    dvv: Ptr[float],
    spheremp: Ptr[float],
    dinv: Ptr[float],
    vtemp: Ptr[float],
    div: Ptr[float],
):
    for j in range(1, np + 1):
        for i in range(1, np + 1):
            v1 = v[_common._vec2_idx(i, j, 1, np)]
            v2 = v[_common._vec2_idx(i, j, 2, np)]
            vtemp[_common._vec2_idx(i, j, 1, np)] = (
                dinv[_common._mat22_idx(i, j, 1, 1, np)] * v1 + dinv[_common._mat22_idx(i, j, 1, 2, np)] * v2
            )
            vtemp[_common._vec2_idx(i, j, 2, np)] = (
                dinv[_common._mat22_idx(i, j, 2, 1, np)] * v1 + dinv[_common._mat22_idx(i, j, 2, 2, np)] * v2
            )

    for n in range(1, np + 1):
        for m in range(1, np + 1):
            div_idx = _common._plane_idx(m, n, np)
            div[div_idx] = 0.0
            for j in range(1, np + 1):
                div[div_idx] = div[div_idx] - (
                    spheremp[_common._plane_idx(j, n, np)] * vtemp[_common._vec2_idx(j, n, 1, np)] * dvv[_common._plane_idx(m, j, np)]
                    + spheremp[_common._plane_idx(m, j, np)] * vtemp[_common._vec2_idx(m, j, 2, np)] * dvv[_common._plane_idx(n, j, np)]
                ) * rrearth


def _gradient_sphere_field(
    np: int,
    rrearth: float,
    s: Ptr[float],
    s_offset: int,
    dvv: Ptr[float],
    dinv: Ptr[float],
    v1: Ptr[float],
    v2: Ptr[float],
    ds: Ptr[float],
):
    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dsdx00 = 0.0
            dsdy00 = 0.0
            for i in range(1, np + 1):
                dsdx00 = dsdx00 + dvv[_common._plane_idx(i, l, np)] * s[s_offset + _common._plane_idx(i, j, np)]
                dsdy00 = dsdy00 + dvv[_common._plane_idx(i, l, np)] * s[s_offset + _common._plane_idx(j, i, np)]
            v1[_common._plane_idx(l, j, np)] = dsdx00 * rrearth
            v2[_common._plane_idx(j, l, np)] = dsdy00 * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _common._plane_idx(i, j, np)
            v1_val = v1[plane_idx]
            v2_val = v2[plane_idx]
            ds[_common._vec2_idx(i, j, 1, np)] = (
                dinv[_common._mat22_idx(i, j, 1, 1, np)] * v1_val + dinv[_common._mat22_idx(i, j, 2, 1, np)] * v2_val
            )
            ds[_common._vec2_idx(i, j, 2, np)] = (
                dinv[_common._mat22_idx(i, j, 1, 2, np)] * v1_val + dinv[_common._mat22_idx(i, j, 2, 2, np)] * v2_val
            )


def _laplace_sphere_wk_field(
    np: int,
    rrearth: float,
    hypervis_power: int,
    hypervis_scaling: int,
    var_coef: int,
    s: Ptr[float],
    s_offset: int,
    dvv: Ptr[float],
    spheremp: Ptr[float],
    dinv: Ptr[float],
    variable_hyperviscosity: Ptr[float],
    tensorvisc: Ptr[float],
    grads: Ptr[float],
    oldgrads: Ptr[float],
    v1: Ptr[float],
    v2: Ptr[float],
    laplace: Ptr[float],
):
    _gradient_sphere_field(np, rrearth, s, s_offset, dvv, dinv, v1, v2, grads)

    if var_coef != 0:
        if hypervis_power != 0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _common._plane_idx(i, j, np)
                    scale = variable_hyperviscosity[plane_idx]
                    grads[_common._vec2_idx(i, j, 1, np)] = grads[_common._vec2_idx(i, j, 1, np)] * scale
                    grads[_common._vec2_idx(i, j, 2, np)] = grads[_common._vec2_idx(i, j, 2, np)] * scale
        elif hypervis_scaling != 0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    oldgrads[_common._vec2_idx(i, j, 1, np)] = grads[_common._vec2_idx(i, j, 1, np)]
                    oldgrads[_common._vec2_idx(i, j, 2, np)] = grads[_common._vec2_idx(i, j, 2, np)]

            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    oldgrad1 = oldgrads[_common._vec2_idx(i, j, 1, np)]
                    oldgrad2 = oldgrads[_common._vec2_idx(i, j, 2, np)]
                    grads[_common._vec2_idx(i, j, 1, np)] = (
                        oldgrad1 * tensorvisc[_common._mat22_idx(i, j, 1, 1, np)]
                        + oldgrad2 * tensorvisc[_common._mat22_idx(i, j, 1, 2, np)]
                    )
                    grads[_common._vec2_idx(i, j, 2, np)] = (
                        oldgrad1 * tensorvisc[_common._mat22_idx(i, j, 2, 1, np)]
                        + oldgrad2 * tensorvisc[_common._mat22_idx(i, j, 2, 2, np)]
                    )

    _divergence_sphere_wk_vec(np, rrearth, grads, dvv, spheremp, dinv, oldgrads, laplace)


def _curl_sphere_wk_testcov(
    np: int,
    rrearth: float,
    s: Ptr[float],
    dvv: Ptr[float],
    mp: Ptr[float],
    d: Ptr[float],
    dscontra: Ptr[float],
    ds: Ptr[float],
):
    for n in range(1, np + 1):
        for m in range(1, np + 1):
            dscontra[_common._vec2_idx(m, n, 1, np)] = 0.0
            dscontra[_common._vec2_idx(m, n, 2, np)] = 0.0
            for j in range(1, np + 1):
                dscontra[_common._vec2_idx(m, n, 1, np)] = dscontra[_common._vec2_idx(m, n, 1, np)] - (
                    mp[_common._plane_idx(m, j, np)] * s[_common._plane_idx(m, j, np)] * dvv[_common._plane_idx(n, j, np)]
                ) * rrearth
                dscontra[_common._vec2_idx(m, n, 2, np)] = dscontra[_common._vec2_idx(m, n, 2, np)] + (
                    mp[_common._plane_idx(j, n, np)] * s[_common._plane_idx(j, n, np)] * dvv[_common._plane_idx(m, j, np)]
                ) * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            ds[_common._vec2_idx(i, j, 1, np)] = (
                d[_common._mat22_idx(i, j, 1, 1, np)] * dscontra[_common._vec2_idx(i, j, 1, np)]
                + d[_common._mat22_idx(i, j, 1, 2, np)] * dscontra[_common._vec2_idx(i, j, 2, np)]
            )
            ds[_common._vec2_idx(i, j, 2, np)] = (
                d[_common._mat22_idx(i, j, 2, 1, np)] * dscontra[_common._vec2_idx(i, j, 1, np)]
                + d[_common._mat22_idx(i, j, 2, 2, np)] * dscontra[_common._vec2_idx(i, j, 2, np)]
            )


def _gradient_sphere_wk_testcov(
    np: int,
    rrearth: float,
    s: Ptr[float],
    dvv: Ptr[float],
    mp: Ptr[float],
    metinv: Ptr[float],
    metdet: Ptr[float],
    d: Ptr[float],
    dscontra: Ptr[float],
    ds: Ptr[float],
):
    for n in range(1, np + 1):
        for m in range(1, np + 1):
            dscontra[_common._vec2_idx(m, n, 1, np)] = 0.0
            dscontra[_common._vec2_idx(m, n, 2, np)] = 0.0
            for j in range(1, np + 1):
                plane_idx = _common._plane_idx(m, n, np)
                dscontra[_common._vec2_idx(m, n, 1, np)] = dscontra[_common._vec2_idx(m, n, 1, np)] - (
                    (
                        mp[_common._plane_idx(j, n, np)]
                        * metinv[_common._mat22_idx(m, n, 1, 1, np)]
                        * metdet[plane_idx]
                        * s[_common._plane_idx(j, n, np)]
                        * dvv[_common._plane_idx(m, j, np)]
                    )
                    + (
                        mp[_common._plane_idx(m, j, np)]
                        * metinv[_common._mat22_idx(m, n, 2, 1, np)]
                        * metdet[plane_idx]
                        * s[_common._plane_idx(m, j, np)]
                        * dvv[_common._plane_idx(n, j, np)]
                    )
                ) * rrearth
                dscontra[_common._vec2_idx(m, n, 2, np)] = dscontra[_common._vec2_idx(m, n, 2, np)] - (
                    (
                        mp[_common._plane_idx(j, n, np)]
                        * metinv[_common._mat22_idx(m, n, 1, 2, np)]
                        * metdet[plane_idx]
                        * s[_common._plane_idx(j, n, np)]
                        * dvv[_common._plane_idx(m, j, np)]
                    )
                    + (
                        mp[_common._plane_idx(m, j, np)]
                        * metinv[_common._mat22_idx(m, n, 2, 2, np)]
                        * metdet[plane_idx]
                        * s[_common._plane_idx(m, j, np)]
                        * dvv[_common._plane_idx(n, j, np)]
                    )
                ) * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            ds[_common._vec2_idx(i, j, 1, np)] = (
                d[_common._mat22_idx(i, j, 1, 1, np)] * dscontra[_common._vec2_idx(i, j, 1, np)]
                + d[_common._mat22_idx(i, j, 1, 2, np)] * dscontra[_common._vec2_idx(i, j, 2, np)]
            )
            ds[_common._vec2_idx(i, j, 2, np)] = (
                d[_common._mat22_idx(i, j, 2, 1, np)] * dscontra[_common._vec2_idx(i, j, 1, np)]
                + d[_common._mat22_idx(i, j, 2, 2, np)] * dscontra[_common._vec2_idx(i, j, 2, np)]
            )


def curl_sphere_wk_testcov_codon(
    np: int,
    rrearth: float,
    s_p: cobj,
    dvv_p: cobj,
    mp_p: cobj,
    d_p: cobj,
    dscontra_p: cobj,
    ds_p: cobj,
):
    _curl_sphere_wk_testcov(
        np,
        rrearth,
        Ptr[float](s_p),
        Ptr[float](dvv_p),
        Ptr[float](mp_p),
        Ptr[float](d_p),
        Ptr[float](dscontra_p),
        Ptr[float](ds_p),
    )


def gradient_sphere_wk_testcov_codon(
    np: int,
    rrearth: float,
    s_p: cobj,
    dvv_p: cobj,
    mp_p: cobj,
    metinv_p: cobj,
    metdet_p: cobj,
    d_p: cobj,
    dscontra_p: cobj,
    ds_p: cobj,
):
    _gradient_sphere_wk_testcov(
        np,
        rrearth,
        Ptr[float](s_p),
        Ptr[float](dvv_p),
        Ptr[float](mp_p),
        Ptr[float](metinv_p),
        Ptr[float](metdet_p),
        Ptr[float](d_p),
        Ptr[float](dscontra_p),
        Ptr[float](ds_p),
    )


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
    _divergence_sphere_wk_vec(
        np,
        rrearth,
        Ptr[float](v_p),
        Ptr[float](dvv_p),
        Ptr[float](spheremp_p),
        Ptr[float](dinv_p),
        Ptr[float](vtemp_p),
        Ptr[float](div_p),
    )


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
    _laplace_sphere_wk_field(
        np,
        rrearth,
        hypervis_power,
        hypervis_scaling,
        var_coef,
        Ptr[float](s_p),
        0,
        Ptr[float](dvv_p),
        Ptr[float](spheremp_p),
        Ptr[float](dinv_p),
        Ptr[float](variable_hyperviscosity_p),
        Ptr[float](tensorvisc_p),
        Ptr[float](grads_p),
        Ptr[float](oldgrads_p),
        Ptr[float](v1_p),
        Ptr[float](v2_p),
        Ptr[float](laplace_p),
    )


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
    v = Ptr[float](v_p)
    mp = Ptr[float](mp_p)
    spheremp = Ptr[float](spheremp_p)
    metinv = Ptr[float](metinv_p)
    metdet = Ptr[float](metdet_p)
    d = Ptr[float](d_p)
    variable_hyperviscosity = Ptr[float](variable_hyperviscosity_p)
    vec_sphere2cart = Ptr[float](vec_sphere2cart_p)
    dum_cart = Ptr[float](dum_cart_p)
    dum_tmp = Ptr[float](dum_tmp_p)
    div = Ptr[float](div_p)
    vor = Ptr[float](vor_p)
    lap_tmp = Ptr[float](lap_tmp_p)
    lap_tmp2 = Ptr[float](lap_tmp2_p)
    work1 = Ptr[float](work1_p)
    work2 = Ptr[float](work2_p)
    laplace = Ptr[float](laplace_p)

    if hypervis_scaling != 0 and var_coef != 0:
        for component in range(1, 4):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    dum_cart[_common._vec3_idx(i, j, component, np)] = (
                        vec_sphere2cart[_common._mat32_idx(i, j, component, 1, np)] * v[_common._vec2_idx(i, j, 1, np)]
                        + vec_sphere2cart[_common._mat32_idx(i, j, component, 2, np)] * v[_common._vec2_idx(i, j, 2, np)]
                    )

        for component in range(1, 4):
            _laplace_sphere_wk_field(
                np,
                rrearth,
                hypervis_power,
                hypervis_scaling,
                var_coef,
                dum_cart,
                (component - 1) * np * np,
                Ptr[float](dvv_p),
                spheremp,
                Ptr[float](dinv_p),
                variable_hyperviscosity,
                Ptr[float](tensorvisc_p),
                work1,
                work2,
                Ptr[float](v1_p),
                Ptr[float](v2_p),
                dum_tmp,
            )
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    dum_cart[_common._vec3_idx(i, j, component, np)] = dum_tmp[_common._plane_idx(i, j, np)]

        for component in range(1, 3):
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    laplace[_common._vec2_idx(i, j, component, np)] = (
                        dum_cart[_common._vec3_idx(i, j, 1, np)] * vec_sphere2cart[_common._mat32_idx(i, j, 1, component, np)]
                        + dum_cart[_common._vec3_idx(i, j, 2, np)] * vec_sphere2cart[_common._mat32_idx(i, j, 2, component, np)]
                        + dum_cart[_common._vec3_idx(i, j, 3, np)] * vec_sphere2cart[_common._mat32_idx(i, j, 3, component, np)]
                    )
    else:
        divergence_sphere_codon(np, rrearth, v_p, dvv_p, metdet_p, dinv_p, rmetdet_p, lap_tmp_p, dum_tmp_p, div_p)
        vorticity_sphere_codon(np, rrearth, v_p, dvv_p, d_p, rmetdet_p, lap_tmp2_p, dum_tmp_p, vor_p)

        if var_coef != 0 and hypervis_power != 0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _common._plane_idx(i, j, np)
                    scale = variable_hyperviscosity[plane_idx]
                    div[plane_idx] = div[plane_idx] * scale
                    vor[plane_idx] = vor[plane_idx] * scale

        if has_nu_ratio != 0:
            for j in range(1, np + 1):
                for i in range(1, np + 1):
                    plane_idx = _common._plane_idx(i, j, np)
                    div[plane_idx] = nu_ratio * div[plane_idx]

        _gradient_sphere_wk_testcov(np, rrearth, div, Ptr[float](dvv_p), mp, metinv, metdet, d, work1, lap_tmp)
        _curl_sphere_wk_testcov(np, rrearth, vor, Ptr[float](dvv_p), mp, d, work2, lap_tmp2)

        rrearth_sq = rrearth * rrearth
        for n in range(1, np + 1):
            for m in range(1, np + 1):
                plane_idx = _common._plane_idx(m, n, np)
                laplace[_common._vec2_idx(m, n, 1, np)] = (
                    lap_tmp[_common._vec2_idx(m, n, 1, np)]
                    - lap_tmp2[_common._vec2_idx(m, n, 1, np)]
                    + 2.0 * spheremp[plane_idx] * v[_common._vec2_idx(m, n, 1, np)] * rrearth_sq
                )
                laplace[_common._vec2_idx(m, n, 2, np)] = (
                    lap_tmp[_common._vec2_idx(m, n, 2, np)]
                    - lap_tmp2[_common._vec2_idx(m, n, 2, np)]
                    + 2.0 * spheremp[plane_idx] * v[_common._vec2_idx(m, n, 2, np)] * rrearth_sq
                )


def vlaplace_sphere_wk_contra_codon(
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
    vlaplace_sphere_wk_codon(
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
    v = Ptr[float](v_p)
    dvv = Ptr[float](dvv_p)
    d = Ptr[float](d_p)
    rmetdet = Ptr[float](rmetdet_p)
    vco = Ptr[float](vco_p)
    vtemp = Ptr[float](vtemp_p)
    vort = Ptr[float](vort_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            v1 = v[_common._vec2_idx(i, j, 1, np)]
            v2 = v[_common._vec2_idx(i, j, 2, np)]
            vco[_common._vec2_idx(i, j, 1, np)] = (
                d[_common._mat22_idx(i, j, 1, 1, np)] * v1 + d[_common._mat22_idx(i, j, 2, 1, np)] * v2
            )
            vco[_common._vec2_idx(i, j, 2, np)] = (
                d[_common._mat22_idx(i, j, 1, 2, np)] * v1 + d[_common._mat22_idx(i, j, 2, 2, np)] * v2
            )

    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dudy00 = 0.0
            dvdx00 = 0.0
            for i in range(1, np + 1):
                dvdx00 = dvdx00 + dvv[_common._plane_idx(i, l, np)] * vco[_common._vec2_idx(i, j, 2, np)]
                dudy00 = dudy00 + dvv[_common._plane_idx(i, l, np)] * vco[_common._vec2_idx(j, i, 1, np)]
            vort[_common._plane_idx(l, j, np)] = dvdx00
            vtemp[_common._plane_idx(j, l, np)] = dudy00

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _common._plane_idx(i, j, np)
            vort[plane_idx] = (vort[plane_idx] - vtemp[plane_idx]) * (rmetdet[plane_idx] * rrearth)


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
    s = Ptr[float](s_p)
    dvv = Ptr[float](dvv_p)
    dinv = Ptr[float](dinv_p)
    v1 = Ptr[float](v1_p)
    v2 = Ptr[float](v2_p)
    ds = Ptr[float](ds_p)

    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dsdx00 = 0.0
            dsdy00 = 0.0
            for i in range(1, np + 1):
                dsdx00 = dsdx00 + dvv[_common._plane_idx(i, l, np)] * s[_common._plane_idx(i, j, np)]
                dsdy00 = dsdy00 + dvv[_common._plane_idx(i, l, np)] * s[_common._plane_idx(j, i, np)]
            v1[_common._plane_idx(l, j, np)] = dsdx00 * rrearth
            v2[_common._plane_idx(j, l, np)] = dsdy00 * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _common._plane_idx(i, j, np)
            v1_val = v1[plane_idx]
            v2_val = v2[plane_idx]
            ds[_common._vec2_idx(i, j, 1, np)] = (
                dinv[_common._mat22_idx(i, j, 1, 1, np)] * v1_val + dinv[_common._mat22_idx(i, j, 2, 1, np)] * v2_val
            )
            ds[_common._vec2_idx(i, j, 2, np)] = (
                dinv[_common._mat22_idx(i, j, 1, 2, np)] * v1_val + dinv[_common._mat22_idx(i, j, 2, 2, np)] * v2_val
            )


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
    s = Ptr[float](s_p)
    dvv = Ptr[float](dvv_p)
    d = Ptr[float](d_p)
    metdet = Ptr[float](metdet_p)
    v1 = Ptr[float](v1_p)
    v2 = Ptr[float](v2_p)
    ds = Ptr[float](ds_p)

    for j in range(1, np + 1):
        for l in range(1, np + 1):
            dsdx00 = 0.0
            dsdy00 = 0.0
            for i in range(1, np + 1):
                dsdx00 = dsdx00 + dvv[_common._plane_idx(i, l, np)] * s[_common._plane_idx(i, j, np)]
                dsdy00 = dsdy00 + dvv[_common._plane_idx(i, l, np)] * s[_common._plane_idx(j, i, np)]
            v2[_common._plane_idx(l, j, np)] = -dsdx00 * rrearth
            v1[_common._plane_idx(j, l, np)] = dsdy00 * rrearth

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            plane_idx = _common._plane_idx(i, j, np)
            v1_val = v1[plane_idx]
            v2_val = v2[plane_idx]
            ds[_common._vec2_idx(i, j, 1, np)] = (
                d[_common._mat22_idx(i, j, 1, 1, np)] * v1_val + d[_common._mat22_idx(i, j, 1, 2, np)] * v2_val
            ) / metdet[plane_idx]
            ds[_common._vec2_idx(i, j, 2, np)] = (
                d[_common._mat22_idx(i, j, 2, 1, np)] * v1_val + d[_common._mat22_idx(i, j, 2, 2, np)] * v2_val
            ) / metdet[plane_idx]


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
    u = Ptr[float](u_p)
    v = Ptr[float](v_p)
    dvv = Ptr[float](dvv_p)
    dinv = Ptr[float](dinv_p)
    vec_sphere2cart = Ptr[float](vec_sphere2cart_p)
    dum_cart = Ptr[float](dum_cart_p)
    tmp = Ptr[float](tmp_p)
    v1 = Ptr[float](v1_p)
    v2 = Ptr[float](v2_p)
    ugradv = Ptr[float](ugradv_p)

    for component in range(1, 4):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                dum_cart[_common._vec3_idx(i, j, component, np)] = (
                    vec_sphere2cart[_common._mat32_idx(i, j, component, 1, np)] * v[_common._vec2_idx(i, j, 1, np)]
                    + vec_sphere2cart[_common._mat32_idx(i, j, component, 2, np)] * v[_common._vec2_idx(i, j, 2, np)]
                )

    for component in range(1, 4):
        for j in range(1, np + 1):
            for l in range(1, np + 1):
                dsdx00 = 0.0
                dsdy00 = 0.0
                for i in range(1, np + 1):
                    dsdx00 = dsdx00 + dvv[_common._plane_idx(i, l, np)] * dum_cart[_common._vec3_idx(i, j, component, np)]
                    dsdy00 = dsdy00 + dvv[_common._plane_idx(i, l, np)] * dum_cart[_common._vec3_idx(j, i, component, np)]
                v1[_common._plane_idx(l, j, np)] = dsdx00 * rrearth
                v2[_common._plane_idx(j, l, np)] = dsdy00 * rrearth

        for j in range(1, np + 1):
            for i in range(1, np + 1):
                plane_idx = _common._plane_idx(i, j, np)
                v1_val = v1[plane_idx]
                v2_val = v2[plane_idx]
                tmp[_common._vec2_idx(i, j, 1, np)] = (
                    dinv[_common._mat22_idx(i, j, 1, 1, np)] * v1_val + dinv[_common._mat22_idx(i, j, 2, 1, np)] * v2_val
                )
                tmp[_common._vec2_idx(i, j, 2, np)] = (
                    dinv[_common._mat22_idx(i, j, 1, 2, np)] * v1_val + dinv[_common._mat22_idx(i, j, 2, 2, np)] * v2_val
                )
                dum_cart[_common._vec3_idx(i, j, component, np)] = (
                    u[_common._vec2_idx(i, j, 1, np)] * tmp[_common._vec2_idx(i, j, 1, np)]
                    + u[_common._vec2_idx(i, j, 2, np)] * tmp[_common._vec2_idx(i, j, 2, np)]
                )

    for component in range(1, 3):
        for j in range(1, np + 1):
            for i in range(1, np + 1):
                ugradv[_common._vec2_idx(i, j, component, np)] = (
                    dum_cart[_common._vec3_idx(i, j, 1, np)] * vec_sphere2cart[_common._mat32_idx(i, j, 1, component, np)]
                    + dum_cart[_common._vec3_idx(i, j, 2, np)] * vec_sphere2cart[_common._mat32_idx(i, j, 2, component, np)]
                    + dum_cart[_common._vec3_idx(i, j, 3, np)] * vec_sphere2cart[_common._mat32_idx(i, j, 3, component, np)]
                )
