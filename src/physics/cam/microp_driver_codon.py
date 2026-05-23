@export
def microp_driver_readnl_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def microp_driver_register_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def microp_driver_init_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def microp_driver_implements_cnst_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def microp_driver_select_scheme_codon(
    scheme_len: int,
    scheme_ascii_p: cobj,
    scheme_code_p: cobj,
    status_p: cobj,
):
    scheme_ascii = Ptr[int](scheme_ascii_p)
    scheme_code = Ptr[int](scheme_code_p)
    status = Ptr[int](status_p)

    status[0] = 0
    scheme_code[0] = 0

    n = scheme_len
    while n > 0 and scheme_ascii[n - 1] == 32:
        n -= 1

    if n == 2:
        c1 = scheme_ascii[0]
        c2 = scheme_ascii[1]

        if c1 >= 65 and c1 <= 90:
            c1 += 32
        if c2 >= 65 and c2 <= 90:
            c2 += 32

        if c1 == 109 and c2 == 103:
            scheme_code[0] = 1
            return

        if c1 == 114 and c2 == 107:
            scheme_code[0] = 2
            return

    status[0] = 1
