from chemistry_common_codon import _idx2, _idx3

@inline
def _aero_model_gasaerexch_column_flux(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass: float,
    gravit: float,
    field: Ptr[float],
    mbar: Ptr[float],
    pdel: Ptr[float],
    wrk: Ptr[float],
):
    for i in range(1, ncol + 1):
        wrk[i - 1] = 0.0

    for k in range(1, pver + 1):
        for i in range(1, ncol + 1):
            wrk[i - 1] += field[_idx2(i, k, ncol)] * adv_mass / mbar[_idx2(i, k, pcols)] * pdel[
                _idx2(i, k, pcols)
            ] / gravit

@inline
def _aero_model_gasaerexch_all_column_fluxes(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    gravit: float,
    field: Ptr[float],
    mbar: Ptr[float],
    pdel: Ptr[float],
    adv_mass: Ptr[float],
    wrk: Ptr[float],
):
    for m in range(1, gas_pcnst + 1):
        mass = adv_mass[m - 1]
        for i in range(1, ncol + 1):
            wrk[_idx2(i, m, ncol)] = 0.0

        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                wrk[_idx2(i, m, ncol)] += (
                    field[_idx3(i, k, m, ncol, pver)]
                    * mass
                    / mbar[_idx2(i, k, pcols)]
                    * pdel[_idx2(i, k, pcols)]
                    / gravit
                )

@inline
def _aero_model_gasaerexch_h2so4_save_or_delta(
    ncol: int,
    pver: int,
    ndx_h2so4: int,
    stage3_mode: int,
    vmr: Ptr[float],
    del_h2so4_aeruptk: Ptr[float],
):
    if stage3_mode == 0:
        if ndx_h2so4 > 0:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    del_h2so4_aeruptk[_idx2(i, k, ncol)] = vmr[_idx3(i, k, ndx_h2so4, ncol, pver)]
        else:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    del_h2so4_aeruptk[_idx2(i, k, ncol)] = 0.0
        return

    if ndx_h2so4 > 0:
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx2(i, k, ncol)
                del_h2so4_aeruptk[idx] = vmr[_idx3(i, k, ndx_h2so4, ncol, pver)] - del_h2so4_aeruptk[idx]

@inline
def _aero_model_gasaerexch_gas_tend(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr0: Ptr[float],
    vmr: Ptr[float],
    dvmrdt: Ptr[float],
):
    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dvmrdt[idx] = (vmr[idx] - vmr0[idx]) / delt

@inline
def _aero_model_gasaerexch_aq_tend(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr: Ptr[float],
    vmrcw: Ptr[float],
    dvmrdt: Ptr[float],
    dvmrcwdt: Ptr[float],
):
    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dvmrdt[idx] = (vmr[idx] - dvmrdt[idx]) / delt

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dvmrcwdt[idx] = (vmrcw[idx] - dvmrcwdt[idx]) / delt

def _aero_model_gasaerexch_vmrcw_batch(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    qqcw_offset: int,
    mbar_ld1: int,
    qqcw_ptrs: Ptr[cobj],
    qqcw_present: Ptr[int],
    mbar: Ptr[float],
    adv_mass: Ptr[float],
    vmr: Ptr[float],
):
    for m in range(1, gas_pcnst + 1):
        if adv_mass[m - 1] == 0.0:
            continue

        qqcw_index = m + qqcw_offset
        if qqcw_present[qqcw_index - 1] != 0:
            fldcw = Ptr[float](qqcw_ptrs[qqcw_index - 1])
            if mode == 1:
                for k in range(1, pver + 1):
                    for i in range(1, ncol + 1):
                        vmr_idx = _idx3(i, k, m, ncol, pver)
                        vmr[vmr_idx] = (
                            mbar[_idx2(i, k, mbar_ld1)] * fldcw[_idx2(i, k, pcols)] / adv_mass[m - 1]
                        )
            else:
                for k in range(1, pver + 1):
                    for i in range(1, ncol + 1):
                        vmr_idx = _idx3(i, k, m, ncol, pver)
                        fldcw[_idx2(i, k, pcols)] = (
                            adv_mass[m - 1] * vmr[vmr_idx] / mbar[_idx2(i, k, mbar_ld1)]
                        )
        elif mode == 1:
            for k in range(1, pver + 1):
                for i in range(1, ncol + 1):
                    vmr[_idx3(i, k, m, ncol, pver)] = 0.0

def _aero_model_gasaerexch_snapshot_state(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    vmr: Ptr[float],
    vmrcw: Ptr[float],
    dvmrdt: Ptr[float],
    dvmrcwdt: Ptr[float],
):
    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dvmrdt[idx] = vmr[idx]

    for m in range(1, gas_pcnst + 1):
        for k in range(1, pver + 1):
            for i in range(1, ncol + 1):
                idx = _idx3(i, k, m, ncol, pver)
                dvmrcwdt[idx] = vmrcw[idx]

def aero_model_gasaerexch_codon(
    stage: int,
    stage3_mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    delt: float,
    gravit: float,
    vmr0_p: cobj,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    adv_mass_p: cobj,
    wrk_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    dvmrdt = Ptr[float](dvmrdt_p)

    if stage == 1:
        vmr0 = Ptr[float](vmr0_p)
        mbar = Ptr[float](mbar_p)
        pdel = Ptr[float](pdel_p)
        adv_mass = Ptr[float](adv_mass_p)
        wrk = Ptr[float](wrk_p)

        _aero_model_gasaerexch_gas_tend(ncol, pver, gas_pcnst, delt, vmr0, vmr, dvmrdt)
        _aero_model_gasaerexch_all_column_fluxes(ncol, pcols, pver, gas_pcnst, gravit, dvmrdt, mbar, pdel, adv_mass, wrk)
        return

    if stage == 2:
        vmrcw = Ptr[float](vmrcw_p)
        dvmrcwdt = Ptr[float](dvmrcwdt_p)
        mbar = Ptr[float](mbar_p)
        pdel = Ptr[float](pdel_p)
        adv_mass = Ptr[float](adv_mass_p)
        wrk = Ptr[float](wrk_p)

        _aero_model_gasaerexch_aq_tend(ncol, pver, gas_pcnst, delt, vmr, vmrcw, dvmrdt, dvmrcwdt)
        _aero_model_gasaerexch_all_column_fluxes(ncol, pcols, pver, gas_pcnst, gravit, dvmrdt, mbar, pdel, adv_mass, wrk)
        return

    if stage == 3:
        del_h2so4_aeruptk = Ptr[float](del_h2so4_aeruptk_p)
        _aero_model_gasaerexch_h2so4_save_or_delta(
            ncol, pver, ndx_h2so4, stage3_mode, vmr, del_h2so4_aeruptk
        )
        return

    if stage == 4:
        vmrcw = Ptr[float](vmrcw_p)
        dvmrcwdt = Ptr[float](dvmrcwdt_p)
        mbar = Ptr[float](mbar_p)
        pdel = Ptr[float](pdel_p)
        adv_mass = Ptr[float](adv_mass_p)
        wrk = Ptr[float](wrk_p)
        del_h2so4_aeruptk = Ptr[float](del_h2so4_aeruptk_p)

        _aero_model_gasaerexch_aq_tend(ncol, pver, gas_pcnst, delt, vmr, vmrcw, dvmrdt, dvmrcwdt)
        _aero_model_gasaerexch_all_column_fluxes(ncol, pcols, pver, gas_pcnst, gravit, dvmrdt, mbar, pdel, adv_mass, wrk)
        _aero_model_gasaerexch_h2so4_save_or_delta(ncol, pver, ndx_h2so4, stage3_mode, vmr, del_h2so4_aeruptk)

def aero_model_gasaerexch_presetsox_shell_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    gravit: float,
    vmr0_p: cobj,
    vmr_p: cobj,
    dvmrdt_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    adv_mass_p: cobj,
    wrk_p: cobj,
):
    vmr0 = Ptr[float](vmr0_p)
    vmr = Ptr[float](vmr_p)
    dvmrdt = Ptr[float](dvmrdt_p)
    mbar = Ptr[float](mbar_p)
    pdel = Ptr[float](pdel_p)
    adv_mass = Ptr[float](adv_mass_p)
    wrk = Ptr[float](wrk_p)

    _aero_model_gasaerexch_gas_tend(ncol, pver, gas_pcnst, delt, vmr0, vmr, dvmrdt)
    _aero_model_gasaerexch_all_column_fluxes(ncol, pcols, pver, gas_pcnst, gravit, dvmrdt, mbar, pdel, adv_mass, wrk)

def aero_model_gasaerexch_column_flux_codon(
    ncol: int,
    pcols: int,
    pver: int,
    adv_mass: float,
    gravit: float,
    field_p: cobj,
    mbar_p: cobj,
    pdel_p: cobj,
    wrk_p: cobj,
):
    field = Ptr[float](field_p)
    mbar = Ptr[float](mbar_p)
    pdel = Ptr[float](pdel_p)
    wrk = Ptr[float](wrk_p)

    _aero_model_gasaerexch_column_flux(ncol, pcols, pver, adv_mass, gravit, field, mbar, pdel, wrk)

def aero_model_gasaerexch_h2so4_save_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    del_h2so4_aeruptk = Ptr[float](del_h2so4_aeruptk_p)

    _aero_model_gasaerexch_h2so4_save_or_delta(ncol, pver, ndx_h2so4, 0, vmr, del_h2so4_aeruptk)

def aero_model_gasaerexch_h2so4_delta_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    ndx_h2so4: int,
    vmr_p: cobj,
    del_h2so4_aeruptk_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    del_h2so4_aeruptk = Ptr[float](del_h2so4_aeruptk_p)

    _aero_model_gasaerexch_h2so4_save_or_delta(ncol, pver, ndx_h2so4, 1, vmr, del_h2so4_aeruptk)

def aero_model_gasaerexch_gas_tend_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr0_p: cobj,
    vmr_p: cobj,
    dvmrdt_p: cobj,
):
    vmr0 = Ptr[float](vmr0_p)
    vmr = Ptr[float](vmr_p)
    dvmrdt = Ptr[float](dvmrdt_p)

    _aero_model_gasaerexch_gas_tend(ncol, pver, gas_pcnst, delt, vmr0, vmr, dvmrdt)

def aero_model_gasaerexch_aq_tend_codon(
    ncol: int,
    pver: int,
    gas_pcnst: int,
    delt: float,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
):
    vmr = Ptr[float](vmr_p)
    vmrcw = Ptr[float](vmrcw_p)
    dvmrdt = Ptr[float](dvmrdt_p)
    dvmrcwdt = Ptr[float](dvmrcwdt_p)

    _aero_model_gasaerexch_aq_tend(ncol, pver, gas_pcnst, delt, vmr, vmrcw, dvmrdt, dvmrcwdt)

def aero_model_gasaerexch_vmrcw_batch_codon(
    mode: int,
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    qqcw_offset: int,
    mbar_ld1: int,
    qqcw_ptrs_p: cobj,
    qqcw_present_p: cobj,
    mbar_p: cobj,
    adv_mass_p: cobj,
    vmr_p: cobj,
):
    qqcw_ptrs = Ptr[cobj](qqcw_ptrs_p)
    qqcw_present = Ptr[int](qqcw_present_p)
    mbar = Ptr[float](mbar_p)
    adv_mass = Ptr[float](adv_mass_p)
    vmr = Ptr[float](vmr_p)
    _aero_model_gasaerexch_vmrcw_batch(
        mode, ncol, pcols, pver, gas_pcnst, qqcw_offset, mbar_ld1, qqcw_ptrs, qqcw_present, mbar, adv_mass, vmr
    )

def aero_model_gasaerexch_load_snapshot_codon(
    ncol: int,
    pcols: int,
    pver: int,
    gas_pcnst: int,
    qqcw_offset: int,
    mbar_ld1: int,
    qqcw_ptrs_p: cobj,
    qqcw_present_p: cobj,
    mbar_p: cobj,
    adv_mass_p: cobj,
    vmr_p: cobj,
    vmrcw_p: cobj,
    dvmrdt_p: cobj,
    dvmrcwdt_p: cobj,
):
    qqcw_ptrs = Ptr[cobj](qqcw_ptrs_p)
    qqcw_present = Ptr[int](qqcw_present_p)
    mbar = Ptr[float](mbar_p)
    adv_mass = Ptr[float](adv_mass_p)
    vmr = Ptr[float](vmr_p)
    vmrcw = Ptr[float](vmrcw_p)
    dvmrdt = Ptr[float](dvmrdt_p)
    dvmrcwdt = Ptr[float](dvmrcwdt_p)

    _aero_model_gasaerexch_vmrcw_batch(
        1, ncol, pcols, pver, gas_pcnst, qqcw_offset, mbar_ld1, qqcw_ptrs, qqcw_present, mbar, adv_mass, vmrcw
    )
    _aero_model_gasaerexch_snapshot_state(ncol, pver, gas_pcnst, vmr, vmrcw, dvmrdt, dvmrcwdt)
