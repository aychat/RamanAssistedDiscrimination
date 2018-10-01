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

        ('lower_bounds', POINTER(c_double)),
        ('upper_bounds', POINTER(c_double)),
        ('guess', POINTER(c_double)),
        ('MAX_EVAL', c_int)
    ]


class Parameters_AbsSpectra(Structure):
    """
    Parameters structure ctypes
    """
    _fields_ = [
        ('time_spectra_abs', POINTER(c_double)),
        ('frequency_abs', POINTER(c_double)),

        ('A_abs', c_double),
        ('width_abs', c_double),

        ('nDIM', c_int),
        ('timeDIM_abs', c_int),
        ('freqDIM_abs', c_int),

        ('field_abs', POINTER(c_complex))
    ]


class Parameters_VibSpectra(Structure):
    """
    Parameters structure ctypes
    """
    _fields_ = [
        ('time_spectra_vib', POINTER(c_double)),
        ('frequency_vib', POINTER(c_double)),

        ('A_vib', c_double),
        ('width_vib', c_double),

        ('nDIM', c_int),
        ('timeDIM_vib', c_int),
        ('freqDIM_vib', c_int),

        ('field_vib', POINTER(c_complex))
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
        ('rho_0', POINTER(c_complex)),
        ('abs_spectra', POINTER(c_double)),
        ('vib_spectra', POINTER(c_double)),
    ]


try:
    # Load the shared library assuming that it is in the same directory
    lib = ctypes.cdll.LoadLibrary(os.getcwd() + "/RamanAssistedControl.so")
except OSError:
    raise NotImplementedError(
        """
        The library is absent. You must compile the C shared library using the commands:
        gcc -O3 -shared -o RamanAssistedControl.so RamanAssistedControl.c -lnlopt -lm -fopenmp -fPIC
        """
    )

#####################################################
#                                                   #
#          Declaring RamanControlFunction           #
#                                                   #
#####################################################

lib.RamanControlFunction.argtypes = (
    POINTER(Molecule),      # molecule molA
    POINTER(Molecule),      # molecule molB
    POINTER(Parameters),    # parameter field_params
)
lib.RamanControlFunction.restype = POINTER(c_complex)


def RamanControlFunction(molA, molB, func_params):
    return lib.RamanControlFunction(
        molA,
        molB,
        func_params
    )


try:
    # Load the shared library assuming that it is in the same directory
    lib1 = ctypes.cdll.LoadLibrary(os.getcwd() + "/Spectra.so")
except OSError:
    raise NotImplementedError(
        """
        The library is absent. You must compile the C shared library using the commands:
        gcc -O3 -shared -o Spectra.so Spectra.c -lnlopt -lm -fopenmp -fPIC
        """
    )

#####################################################
#                                                   #
#          Declaring CalculateSpectra function      #
#                                                   #
#####################################################

lib1.CalculateSpectra.argtypes = (
    POINTER(Molecule),                  # molecule molA
    POINTER(Molecule),                  # molecule molB
    POINTER(Parameters_AbsSpectra),     # parameter field_params
    POINTER(Parameters_VibSpectra),     # parameter field_params
)
lib1.CalculateSpectra.restype = POINTER(c_complex)


def CalculateSpectra(molA, molB, abs_spec_params, vib_spec_params):
    return lib1.CalculateSpectra(
        molA,
        molB,
        abs_spec_params,
        vib_spec_params
    )
