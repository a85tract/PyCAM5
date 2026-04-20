@inline
def _idx2(i: int, j: int, ld1: int) -> int:
    return (i - 1) + (j - 1) * ld1


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
