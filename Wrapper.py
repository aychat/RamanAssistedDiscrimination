import os
from ctypes import cdll, c_int, c_double, POINTER, Structure


__doc__ = """
Python wrapper for RamanControl.c
Compile with:
gcc -O3 -shared -o RamanControl.so RamanControl.c -lnlopt -lm -fopenmp -fPIC
"""


class c_complex(Structure):
    """
    Complex double ctypes
    """
    _fields_ = [
        ('real', c_double),
        ('imag', c_double)
    ]


class Parameters(Structure):
    """
    Parameters structure ctypes
    """
    _fields_ = [
        ('time', POINTER(c_double)),
        ('rho_0', POINTER(c_complex)),
        ('mu', POINTER(c_complex)),

        ('A_R', c_double),
        ('width_R', c_double),
        ('t0_R', c_double),

        ('A_EE', c_double),
        ('width_EE', c_double),
        ('t0_EE', c_double),

        ('w_R', c_double),
        ('w_v', c_double),
        ('w_EE', c_double),

        ('nDIM', c_int),
        ('timeDIM', c_int),

        ('field_out', POINTER(c_complex)),
    ]


class Molecule(Structure):
    """
    Parameters structure ctypes
    """
    _fields_ = [
        ('energies', POINTER(c_double)),
        ('gamma_decay', POINTER(c_double)),
        ('gamma_pcd', POINTER(c_double)),
        ('rho', POINTER(c_complex)),
        ('dyn_rho', POINTER(c_complex))
    ]

try:
    # Load the shared library assuming that it is in the same directory
    lib = cdll.LoadLibrary(os.getcwd() + "/RamanControl.so")
except OSError:
    raise NotImplementedError(
        """
        The library is absent. You must compile the C shared library using the commands:
        gcc -O3 -shared -o RamanControl.so RamanControl.c -lnlopt -lm -fopenmp -fPIC
        """
    )

#####################################################
#                                                   #
#          Wrapper call for C functions             #
#                                                   #
#####################################################


lib.RamanControlFunction.argtypes = (
    POINTER(Molecule),      # molecule molA
    POINTER(Molecule),      # molecule molB
    POINTER(Parameters),    # parameter field_params
)
lib.RamanControlFunction.restype = None


def RamanControlFunction(molA, molB, func_params):
    return lib.RamanControlFunction(
        molA,
        molB,
        func_params
    )

