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
def microp_driver_tend_codon(flag: int) -> int:
    if flag != 0:
        return 1
    return 0


@export
def microp_driver_implements_cnst_codon(
    scheme_len: int,
    scheme_ascii_p: cobj,
    name_len: int,
    name_ascii_p: cobj,
) -> int:
    scheme_ascii = Ptr[int](scheme_ascii_p)
    name_ascii = Ptr[int](name_ascii_p)
    if _scheme_is_mg(scheme_len, scheme_ascii) == 0:
        return 0
    return _micro_mg_name_match(name_len, name_ascii)


@inline
def _scheme_is_mg(scheme_len: int, scheme_ascii: Ptr[int]) -> int:
    n = scheme_len
    while n > 0 and scheme_ascii[n - 1] == 32:
        n -= 1
    if n != 2:
        return 0
    c1 = scheme_ascii[0]
    c2 = scheme_ascii[1]
    if c1 >= 65 and c1 <= 90:
        c1 += 32
    if c2 >= 65 and c2 <= 90:
        c2 += 32
    if c1 == 109 and c2 == 103:
        return 1
    return 0


@inline
def _micro_mg_name_match(name_len: int, name_ascii: Ptr[int]) -> int:
    if _name_eq8(name_len, name_ascii, 67, 76, 68, 76, 73, 81, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 67, 76, 68, 73, 67, 69, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 78, 85, 77, 76, 73, 81, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 78, 85, 77, 73, 67, 69, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 82, 65, 73, 78, 81, 77, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 83, 78, 79, 87, 81, 77, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 78, 85, 77, 82, 65, 73, 32, 32) != 0:
        return 1
    if _name_eq8(name_len, name_ascii, 78, 85, 77, 83, 78, 79, 32, 32) != 0:
        return 1
    return 0


@inline
def _name_eq8(
    name_len: int,
    name_ascii: Ptr[int],
    c0: int,
    c1: int,
    c2: int,
    c3: int,
    c4: int,
    c5: int,
    c6: int,
    c7: int,
) -> int:
    if name_len > 8:
        i = 8
        while i < name_len:
            if name_ascii[i] != 32:
                return 0
            i += 1
    values = (c0, c1, c2, c3, c4, c5, c6, c7)
    i = 0
    while i < 8:
        left = 32
        if i < name_len:
            left = name_ascii[i]
        if left != values[i]:
            return 0
        i += 1
    return 1


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
