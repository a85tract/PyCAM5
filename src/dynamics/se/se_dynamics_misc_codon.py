def se_misc_touch_codon(tag: int) -> int:
    return tag


from C import gbarrier_initialize(cobj, i32)
from C import gbarrier_free(cobj)
from C import gbarrier_synchronize(cobj, i32)


def gbarrier_init_codon(c_barrier_p: cobj, nthreads: int):
    gbarrier_initialize(c_barrier_p, i32(nthreads))


def gbarrier_delete_codon(c_barrier_p: cobj):
    gbarrier_free(c_barrier_p)


def gbarrier_synchronize_codon(c_barrier: cobj, thread: int):
    gbarrier_synchronize(c_barrier, i32(thread))
