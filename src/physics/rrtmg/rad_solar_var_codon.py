@export
def integrate_spectrum_codon(
    nsrc: int,
    ntrg: int,
    src_x_p: cobj,
    min_trg_p: cobj,
    max_trg_p: cobj,
    src_p: cobj,
    trg_p: cobj,
):
    src_x = Ptr[float](src_x_p)
    min_trg = Ptr[float](min_trg_p)
    max_trg = Ptr[float](max_trg_p)
    src = Ptr[float](src_p)
    trg = Ptr[float](trg_p)

    for i in range(1, ntrg + 1):
        tl = min_trg[i - 1]
        tu = max_trg[i - 1]
        if tl < src_x[nsrc]:
            sil = 1
            for l in range(1, nsrc + 2):
                if tl <= src_x[l - 1]:
                    sil = l
                    break

            siu = 1
            for l in range(1, nsrc + 2):
                if tu <= src_x[l - 1]:
                    siu = l
                    break

            y = 0.0
            if sil < 2:
                sil = 2
            if siu > nsrc + 1:
                siu = nsrc + 1

            for si in range(sil, siu + 1):
                si1 = si - 1
                src_l = src_x[si1 - 1]
                if tl > src_l:
                    sl = tl
                else:
                    sl = src_l
                src_u = src_x[si - 1]
                if tu < src_u:
                    su = tu
                else:
                    su = src_u
                y = y + (su - sl) * src[si1 - 1]

            targ = y / (tu - tl)
        else:
            targ = 0.0

        trg[i - 1] = targ * (tu - tl)
