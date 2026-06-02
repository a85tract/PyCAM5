import se_dynamics_common_codon as _common


def preq_hydrostatic_codon(
    np: int,
    nlev: int,
    rgas: float,
    phi_p: cobj,
    phis_p: cobj,
    tv_p: cobj,
    p_p: cobj,
    dp_p: cobj,
    phii_p: cobj,
):
    phi = Ptr[float](phi_p)
    phis = Ptr[float](phis_p)
    tv = Ptr[float](tv_p)
    p = Ptr[float](p_p)
    dp = Ptr[float](dp_p)
    phii = Ptr[float](phii_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            idx = _common._vol_idx(i, j, nlev, np)
            plane = _common._plane_idx(i, j, np)
            hkk = dp[idx] * 0.5 / p[idx]
            hkl = 2.0 * hkk
            phii[idx] = rgas * tv[idx] * hkl
            phi[idx] = phis[plane] + rgas * tv[idx] * hkk

        k = nlev - 1
        while k >= 2:
            for i in range(1, np + 1):
                idx = _common._vol_idx(i, j, k, np)
                idxp1 = _common._vol_idx(i, j, k + 1, np)
                plane = _common._plane_idx(i, j, np)
                hkk = dp[idx] * 0.5 / p[idx]
                hkl = 2.0 * hkk
                phii[idx] = phii[idxp1] + rgas * tv[idx] * hkl
                phi[idx] = phis[plane] + phii[idxp1] + rgas * tv[idx] * hkk
            k -= 1

        for i in range(1, np + 1):
            idx = _common._vol_idx(i, j, 1, np)
            idx2 = _common._vol_idx(i, j, 2, np)
            plane = _common._plane_idx(i, j, np)
            hkk = 0.5 * dp[idx] / p[idx]
            phi[idx] = phis[plane] + phii[idx2] + rgas * tv[idx] * hkk


def preq_omega_ps_codon(np: int, nlev: int, omega_p_p: cobj, p_p: cobj, vgrad_p_p: cobj, divdp_p: cobj, suml_p: cobj):
    omega_p = Ptr[float](omega_p_p)
    p = Ptr[float](p_p)
    vgrad_p = Ptr[float](vgrad_p_p)
    divdp = Ptr[float](divdp_p)
    suml = Ptr[float](suml_p)

    for j in range(1, np + 1):
        for i in range(1, np + 1):
            idx = _common._vol_idx(i, j, 1, np)
            plane = _common._plane_idx(i, j, np)
            ckk = 0.5 / p[idx]
            term = divdp[idx]
            omega_p[idx] = vgrad_p[idx] / p[idx]
            omega_p[idx] = omega_p[idx] - ckk * term
            suml[plane] = term

        for k in range(2, nlev):
            for i in range(1, np + 1):
                idx = _common._vol_idx(i, j, k, np)
                plane = _common._plane_idx(i, j, np)
                ckk = 0.5 / p[idx]
                ckl = 2.0 * ckk
                term = divdp[idx]
                omega_p[idx] = vgrad_p[idx] / p[idx]
                omega_p[idx] = omega_p[idx] - ckl * suml[plane] - ckk * term
                suml[plane] = suml[plane] + term

        for i in range(1, np + 1):
            idx = _common._vol_idx(i, j, nlev, np)
            plane = _common._plane_idx(i, j, np)
            ckk = 0.5 / p[idx]
            ckl = 2.0 * ckk
            term = divdp[idx]
            omega_p[idx] = vgrad_p[idx] / p[idx]
            omega_p[idx] = omega_p[idx] - ckl * suml[plane] - ckk * term
