@inline
def _idx(i: int) -> int:
    return i - 1


@inline
def _smooth(
    new: Ptr[float],
    old: Ptr[float],
    res: Ptr[float],
    temp: Ptr[float],
    nstep: int,
    deltat: float,
    ncol: int,
    pcols: int,
):
    for i in range(1, ncol + 1):
        temp[_idx(i)] = new[_idx(i)]

    if nstep > 0:
        for i in range(1, ncol + 1):
            new[_idx(i)] = 0.5 * (new[_idx(i)] + old[_idx(i)])
    else:
        for i in range(1, ncol + 1):
            old[_idx(i)] = new[_idx(i)]
            res[_idx(i)] = 0.0

    for i in range(1, ncol + 1):
        res[_idx(i)] = res[_idx(i)] + temp[_idx(i)] - new[_idx(i)]

    for i in range(1, ncol + 1):
        if abs(res[_idx(i)]) < max(abs(new[_idx(i)]), abs(old[_idx(i)])) * 0.05:
            temp[_idx(i)] = res[_idx(i)]
            res[_idx(i)] = 0.0
        else:
            temp[_idx(i)] = res[_idx(i)] * deltat / 7200.0
            res[_idx(i)] = res[_idx(i)] - temp[_idx(i)]

    if float(nstep) * deltat / 86400.0 < 0.5:
        for i in range(1, pcols + 1):
            temp[_idx(i)] = 0.0
            res[_idx(i)] = 0.0

    for i in range(1, ncol + 1):
        new[_idx(i)] = new[_idx(i)] + temp[_idx(i)]
        old[_idx(i)] = new[_idx(i)]


@export
def flux_avg_init_codon(
    ncol: int,
    pcols: int,
    cam_lhf_p: cobj,
    cam_shf_p: cobj,
    cam_cflx1_p: cobj,
    cam_wsx_p: cobj,
    cam_wsy_p: cobj,
    lhflx_p: cobj,
    shflx_p: cobj,
    qflx_p: cobj,
    taux_p: cobj,
    tauy_p: cobj,
    lhflx_res_p: cobj,
    shflx_res_p: cobj,
    qflx_res_p: cobj,
    taux_res_p: cobj,
    tauy_res_p: cobj,
):
    cam_lhf = Ptr[float](cam_lhf_p)
    cam_shf = Ptr[float](cam_shf_p)
    cam_cflx1 = Ptr[float](cam_cflx1_p)
    cam_wsx = Ptr[float](cam_wsx_p)
    cam_wsy = Ptr[float](cam_wsy_p)

    lhflx = Ptr[float](lhflx_p)
    shflx = Ptr[float](shflx_p)
    qflx = Ptr[float](qflx_p)
    taux = Ptr[float](taux_p)
    tauy = Ptr[float](tauy_p)

    lhflx_res = Ptr[float](lhflx_res_p)
    shflx_res = Ptr[float](shflx_res_p)
    qflx_res = Ptr[float](qflx_res_p)
    taux_res = Ptr[float](taux_res_p)
    tauy_res = Ptr[float](tauy_res_p)

    for i in range(1, ncol + 1):
        lhflx[_idx(i)] = cam_lhf[_idx(i)]
        shflx[_idx(i)] = cam_shf[_idx(i)]
        qflx[_idx(i)] = cam_cflx1[_idx(i)]
        taux[_idx(i)] = cam_wsx[_idx(i)]
        tauy[_idx(i)] = cam_wsy[_idx(i)]

    for i in range(1, pcols + 1):
        lhflx_res[_idx(i)] = 0.0
        shflx_res[_idx(i)] = 0.0
        qflx_res[_idx(i)] = 0.0
        taux_res[_idx(i)] = 0.0
        tauy_res[_idx(i)] = 0.0


@export
def flux_avg_run_codon(
    ncol: int,
    pcols: int,
    nstep: int,
    deltat: float,
    cam_lhf_p: cobj,
    cam_shf_p: cobj,
    cam_wsx_p: cobj,
    cam_wsy_p: cobj,
    cam_cflx1_p: cobj,
    lhflx_p: cobj,
    shflx_p: cobj,
    qflx_p: cobj,
    taux_p: cobj,
    tauy_p: cobj,
    lhflx_res_p: cobj,
    shflx_res_p: cobj,
    qflx_res_p: cobj,
    taux_res_p: cobj,
    tauy_res_p: cobj,
    temp_p: cobj,
):
    cam_lhf = Ptr[float](cam_lhf_p)
    cam_shf = Ptr[float](cam_shf_p)
    cam_wsx = Ptr[float](cam_wsx_p)
    cam_wsy = Ptr[float](cam_wsy_p)
    cam_cflx1 = Ptr[float](cam_cflx1_p)

    lhflx = Ptr[float](lhflx_p)
    shflx = Ptr[float](shflx_p)
    qflx = Ptr[float](qflx_p)
    taux = Ptr[float](taux_p)
    tauy = Ptr[float](tauy_p)

    lhflx_res = Ptr[float](lhflx_res_p)
    shflx_res = Ptr[float](shflx_res_p)
    qflx_res = Ptr[float](qflx_res_p)
    taux_res = Ptr[float](taux_res_p)
    tauy_res = Ptr[float](tauy_res_p)

    temp = Ptr[float](temp_p)

    _smooth(cam_lhf, lhflx, lhflx_res, temp, nstep, deltat, ncol, pcols)
    _smooth(cam_shf, shflx, shflx_res, temp, nstep, deltat, ncol, pcols)
    _smooth(cam_wsx, taux, taux_res, temp, nstep, deltat, ncol, pcols)
    _smooth(cam_wsy, tauy, tauy_res, temp, nstep, deltat, ncol, pcols)
    _smooth(cam_cflx1, qflx, qflx_res, temp, nstep, deltat, ncol, pcols)
