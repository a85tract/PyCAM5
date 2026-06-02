@inline
def _idx2(i: int, j: int, ld1: int) -> int:
    return (i - 1) + (j - 1) * ld1


@export
def modal_aero_deposition_init_codon(
    active: int,
    bc1_present: int,
    pom1_present: int,
    soa1_present: int,
    soa2_present: int,
    dst1_present: int,
    dst3_present: int,
    ncl3_present: int,
    so43_present: int,
    bc4_present: int,
    pom4_present: int,
    bc1_opt: int,
    pom1_opt: int,
    soa1_opt: int,
    soa2_opt: int,
    dst1_opt: int,
    dst3_opt: int,
    ncl3_opt: int,
    so43_opt: int,
    bc4_opt: int,
    pom4_opt: int,
    bc1_lookup: int,
    pom1_lookup: int,
    soa1_lookup: int,
    soa2_lookup: int,
    dst1_lookup: int,
    dst3_lookup: int,
    ncl3_lookup: int,
    so43_lookup: int,
    bc4_lookup: int,
    pom4_lookup: int,
    idx_bc1_p: cobj,
    idx_pom1_p: cobj,
    idx_soa1_p: cobj,
    idx_soa2_p: cobj,
    idx_dst1_p: cobj,
    idx_dst3_p: cobj,
    idx_ncl3_p: cobj,
    idx_so43_p: cobj,
    idx_bc4_p: cobj,
    idx_pom4_p: cobj,
    bin_fluxes_p: cobj,
    initialized_p: cobj,
) -> int:
    idx_bc1 = Ptr[int](idx_bc1_p)
    idx_pom1 = Ptr[int](idx_pom1_p)
    idx_soa1 = Ptr[int](idx_soa1_p)
    idx_soa2 = Ptr[int](idx_soa2_p)
    idx_dst1 = Ptr[int](idx_dst1_p)
    idx_dst3 = Ptr[int](idx_dst3_p)
    idx_ncl3 = Ptr[int](idx_ncl3_p)
    idx_so43 = Ptr[int](idx_so43_p)
    idx_bc4 = Ptr[int](idx_bc4_p)
    idx_pom4 = Ptr[int](idx_pom4_p)
    bin_fluxes = Ptr[int](bin_fluxes_p)
    initialized = Ptr[int](initialized_p)

    if active == 0:
        return active

    idx_bc1[0] = bc1_opt if bc1_present != 0 else bc1_lookup
    idx_pom1[0] = pom1_opt if pom1_present != 0 else pom1_lookup
    idx_soa1[0] = soa1_opt if soa1_present != 0 else soa1_lookup
    idx_soa2[0] = soa2_opt if soa2_present != 0 else soa2_lookup
    idx_dst1[0] = dst1_opt if dst1_present != 0 else dst1_lookup
    idx_dst3[0] = dst3_opt if dst3_present != 0 else dst3_lookup
    idx_ncl3[0] = ncl3_opt if ncl3_present != 0 else ncl3_lookup
    idx_so43[0] = so43_opt if so43_present != 0 else so43_lookup
    idx_bc4[0] = bc4_opt if bc4_present != 0 else bc4_lookup
    idx_pom4[0] = pom4_opt if pom4_present != 0 else pom4_lookup
    bin_fluxes[0] = (
        1
        if idx_dst1[0] > 0 and idx_dst3[0] > 0 and idx_ncl3[0] > 0 and idx_so43[0] > 0
        else 0
    )
    initialized[0] = 1
    return active


@export
def set_srf_wetdep_codon(
    ncol: int,
    pcols: int,
    idx_bc1: int,
    idx_bc4: int,
    idx_pom1: int,
    idx_pom4: int,
    idx_soa1: int,
    idx_soa2: int,
    idx_dst1: int,
    idx_dst3: int,
    aerdepwetis_p: cobj,
    aerdepwetcw_p: cobj,
    bcphiwet_p: cobj,
    ocphiwet_p: cobj,
    dstwet1_p: cobj,
    dstwet2_p: cobj,
    dstwet3_p: cobj,
    dstwet4_p: cobj,
):
    aerdepwetis = Ptr[float](aerdepwetis_p)
    aerdepwetcw = Ptr[float](aerdepwetcw_p)
    bcphiwet = Ptr[float](bcphiwet_p)
    ocphiwet = Ptr[float](ocphiwet_p)
    dstwet1 = Ptr[float](dstwet1_p)
    dstwet2 = Ptr[float](dstwet2_p)
    dstwet3 = Ptr[float](dstwet3_p)
    dstwet4 = Ptr[float](dstwet4_p)

    for i in range(1, pcols + 1):
        bcphiwet[i - 1] = 0.0
        ocphiwet[i - 1] = 0.0

    for i in range(1, ncol + 1):
        if idx_bc1 > 0:
            bcphiwet[i - 1] = (
                bcphiwet[i - 1]
                - (aerdepwetis[_idx2(i, idx_bc1, pcols)] + aerdepwetcw[_idx2(i, idx_bc1, pcols)])
            )
        if idx_bc4 > 0:
            bcphiwet[i - 1] = (
                bcphiwet[i - 1]
                - (aerdepwetis[_idx2(i, idx_bc4, pcols)] + aerdepwetcw[_idx2(i, idx_bc4, pcols)])
            )

        if idx_soa1 > 0:
            ocphiwet[i - 1] = (
                ocphiwet[i - 1]
                - (aerdepwetis[_idx2(i, idx_soa1, pcols)] + aerdepwetcw[_idx2(i, idx_soa1, pcols)])
            )
        if idx_soa2 > 0:
            ocphiwet[i - 1] = (
                ocphiwet[i - 1]
                - (aerdepwetis[_idx2(i, idx_soa2, pcols)] + aerdepwetcw[_idx2(i, idx_soa2, pcols)])
            )
        if idx_pom1 > 0:
            ocphiwet[i - 1] = (
                ocphiwet[i - 1]
                - (aerdepwetis[_idx2(i, idx_pom1, pcols)] + aerdepwetcw[_idx2(i, idx_pom1, pcols)])
            )
        if idx_pom4 > 0:
            ocphiwet[i - 1] = (
                ocphiwet[i - 1]
                - (aerdepwetis[_idx2(i, idx_pom4, pcols)] + aerdepwetcw[_idx2(i, idx_pom4, pcols)])
            )

        dstwet1[i - 1] = -(
            aerdepwetis[_idx2(i, idx_dst1, pcols)] + aerdepwetcw[_idx2(i, idx_dst1, pcols)]
        )
        dstwet2[i - 1] = 0.0
        dstwet3[i - 1] = -(
            aerdepwetis[_idx2(i, idx_dst3, pcols)] + aerdepwetcw[_idx2(i, idx_dst3, pcols)]
        )
        dstwet4[i - 1] = 0.0

        if bcphiwet[i - 1] < 0.0:
            bcphiwet[i - 1] = 0.0
        if ocphiwet[i - 1] < 0.0:
            ocphiwet[i - 1] = 0.0
        if dstwet1[i - 1] < 0.0:
            dstwet1[i - 1] = 0.0
        if dstwet3[i - 1] < 0.0:
            dstwet3[i - 1] = 0.0


@export
def set_srf_drydep_codon(
    ncol: int,
    pcols: int,
    idx_bc1: int,
    idx_bc4: int,
    idx_pom1: int,
    idx_pom4: int,
    idx_soa1: int,
    idx_soa2: int,
    idx_dst1: int,
    idx_dst3: int,
    aerdepdryis_p: cobj,
    aerdepdrycw_p: cobj,
    bcphidry_p: cobj,
    bcphodry_p: cobj,
    ocphidry_p: cobj,
    ocphodry_p: cobj,
    dstdry1_p: cobj,
    dstdry2_p: cobj,
    dstdry3_p: cobj,
    dstdry4_p: cobj,
):
    aerdepdryis = Ptr[float](aerdepdryis_p)
    aerdepdrycw = Ptr[float](aerdepdrycw_p)
    bcphidry = Ptr[float](bcphidry_p)
    bcphodry = Ptr[float](bcphodry_p)
    ocphidry = Ptr[float](ocphidry_p)
    ocphodry = Ptr[float](ocphodry_p)
    dstdry1 = Ptr[float](dstdry1_p)
    dstdry2 = Ptr[float](dstdry2_p)
    dstdry3 = Ptr[float](dstdry3_p)
    dstdry4 = Ptr[float](dstdry4_p)

    for i in range(1, pcols + 1):
        bcphidry[i - 1] = 0.0
        bcphodry[i - 1] = 0.0
        ocphidry[i - 1] = 0.0
        ocphodry[i - 1] = 0.0

    for i in range(1, ncol + 1):
        if idx_bc1 > 0:
            bcphidry[i - 1] = (
                bcphidry[i - 1]
                + aerdepdryis[_idx2(i, idx_bc1, pcols)]
                + aerdepdrycw[_idx2(i, idx_bc1, pcols)]
            )
        if idx_bc4 > 0:
            bcphodry[i - 1] = (
                bcphodry[i - 1]
                + aerdepdryis[_idx2(i, idx_bc4, pcols)]
                + aerdepdrycw[_idx2(i, idx_bc4, pcols)]
            )

        if idx_pom1 > 0:
            ocphidry[i - 1] = (
                ocphidry[i - 1]
                + aerdepdryis[_idx2(i, idx_pom1, pcols)]
                + aerdepdrycw[_idx2(i, idx_pom1, pcols)]
            )
        if idx_pom4 > 0:
            ocphodry[i - 1] = (
                ocphodry[i - 1]
                + aerdepdryis[_idx2(i, idx_pom4, pcols)]
                + aerdepdrycw[_idx2(i, idx_pom4, pcols)]
            )
        if idx_soa1 > 0:
            ocphidry[i - 1] = (
                ocphidry[i - 1]
                + aerdepdryis[_idx2(i, idx_soa1, pcols)]
                + aerdepdrycw[_idx2(i, idx_soa1, pcols)]
            )
        if idx_soa2 > 0:
            ocphodry[i - 1] = (
                ocphodry[i - 1]
                + aerdepdryis[_idx2(i, idx_soa2, pcols)]
                + aerdepdrycw[_idx2(i, idx_soa2, pcols)]
            )

        dstdry1[i - 1] = aerdepdryis[_idx2(i, idx_dst1, pcols)] + aerdepdrycw[
            _idx2(i, idx_dst1, pcols)
        ]
        dstdry2[i - 1] = 0.0
        dstdry3[i - 1] = aerdepdryis[_idx2(i, idx_dst3, pcols)] + aerdepdrycw[
            _idx2(i, idx_dst3, pcols)
        ]
        dstdry4[i - 1] = 0.0

        if bcphidry[i - 1] < 0.0:
            bcphidry[i - 1] = 0.0
        if bcphodry[i - 1] < 0.0:
            bcphodry[i - 1] = 0.0
        if ocphidry[i - 1] < 0.0:
            ocphidry[i - 1] = 0.0
        if ocphodry[i - 1] < 0.0:
            ocphodry[i - 1] = 0.0
        if dstdry1[i - 1] < 0.0:
            dstdry1[i - 1] = 0.0
        if dstdry3[i - 1] < 0.0:
            dstdry3[i - 1] = 0.0
