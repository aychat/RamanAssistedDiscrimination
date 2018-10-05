import os
import ctypes
from ctypes import c_int, c_double, c_char_p, POINTER, Structure


__doc__ = """
Python wrapper for RamanAssistedControl.c
Compile with:
gcc -O3 -shared -o RamanAssistedControl.so RamanAssistedControl.c -lm -fopenmp -fPIC
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

        ('A_exc', c_double),
        ('width_exc', c_double),
        ('w_exc', c_double),

        ('nDIM', c_int),
        ('timeDIM', c_int),

        ('field_out', POINTER(c_complex))
    ]


class Molecule(Structure):
    """
    Parameters structure ctypes
    """
    _fields_ = [
        ('nDIM', c_int),
        ('energies', POINTER(c_double)),
        ('gamma_decay', POINTER(c_double)),
        ('gamma_pure_dephasing', POINTER(c_double)),
        ('mu', POINTER(c_complex)),
        ('rho', POINTER(c_complex)),
        ('dyn_rho', POINTER(c_complex)),
        ('rho_0', POINTER(c_complex))
    ]


try:
    # Load the shared library assuming that it is in the same directory
    lib1 = ctypes.cdll.LoadLibrary(os.getcwd() + "/Dynamics.so")
except OSError:
    raise NotImplementedError(
        """
        The library is absent. You must compile the C shared library using the commands:
        gcc -O3 -shared -o Dynamics.so Dynamics.c -lnlopt -lm -fopenmp -fPIC
        """
    )

#####################################################
#                                                   #
#          Declaring CalculateSpectra function      #
#                                                   #
#####################################################

lib1.Dynamics.argtypes = (
    POINTER(Molecule),                  # molecule molA
    POINTER(Parameters),     # parameter field_params
)
lib1.Dynamics.restype = POINTER(c_complex)


def Dynamics(mol, func_params):
    return lib1.Dynamics(
        mol,
        func_params
    )
