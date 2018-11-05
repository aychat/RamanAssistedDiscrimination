import numpy as np
from types import MethodType, FunctionType
from RamanControl_wrapper import *
from ctypes import c_int, c_double, c_char_p, POINTER, Structure


class ADict(dict):
    """
    Dictionary where you can access keys as attributes
    """
    def __getattr__(self, item):
        try:
            return self[item]
        except KeyError:
            dict.__getattribute__(self, item)


class RamanControl:
    """
    """

    def __init__(self, params, **kwargs):
        """
        __init__ function call to initialize variables from the
        parameters for the class instance provided in __main__ and
        add new variables for use in other functions in this class.
        """

        for name, value in kwargs.items():
            if isinstance(value, FunctionType):
                setattr(self, name, MethodType(value, self, self.__class__))
            else:
                setattr(self, name, value)

        self.time_AE = np.linspace(-params.timeAMP_AE, params.timeAMP_AE, params.timeDIM_AE)
        self.time_VR = np.linspace(-params.timeAMP_VR, params.timeAMP_VR, params.timeDIM_VR)

        self.frequency_A = 1./np.linspace(1./params.frequencyMAX_A, 1./params.frequencyMIN_A, params.frequencyDIM_AE)
        self.frequency_E = 1./np.linspace(1./params.frequencyMAX_E, 1./params.frequencyMIN_E, params.frequencyDIM_AE)
        self.frequency_VR = 1./np.linspace(1./params.frequencyMAX_VR, 1./params.frequencyMIN_VR, params.frequencyDIM_VR)

        self.field_A = np.empty(params.timeDIM_AE, dtype=np.complex)
        self.field_E = np.empty(params.timeDIM_AE, dtype=np.complex)
        self.field_V = np.empty(params.timeDIM_VR, dtype=np.complex)
        self.field_R = np.empty(params.timeDIM_VR, dtype=np.complex)

        self.gamma_decay = np.ascontiguousarray(self.gamma_decay)
        self.gamma_pure_dephasingA = np.ascontiguousarray(self.gamma_pure_dephasingA)
        self.gamma_pure_dephasingB = np.ascontiguousarray(self.gamma_pure_dephasingB)
        self.mu = np.ascontiguousarray(self.mu)
        self.rho_0 = np.ascontiguousarray(params.rho_0_A)
        self.rhoA = np.ascontiguousarray(params.rho_0_A.copy())
        self.rhoB = np.ascontiguousarray(params.rho_0_A.copy())
        self.energies_A = np.ascontiguousarray(self.energies_A)
        self.energies_B = np.ascontiguousarray(self.energies_B)

        N = len(self.energies_A)

        self.abs_spectraA = np.ascontiguousarray(np.zeros(len(self.frequency_A)))
        self.abs_spectraB = np.ascontiguousarray(np.zeros(len(self.frequency_A)))
        self.ems_spectraA = np.ascontiguousarray(np.zeros(len(self.frequency_E)))
        self.ems_spectraB = np.ascontiguousarray(np.zeros(len(self.frequency_E)))
        self.vib_spectraA = np.ascontiguousarray(np.zeros(len(self.frequency_VR)))
        self.vib_spectraB = np.ascontiguousarray(np.zeros(len(self.frequency_VR)))
        self.Raman_spectraA = np.ascontiguousarray(np.zeros(len(self.frequency_VR)))
        self.Raman_spectraB = np.ascontiguousarray(np.zeros(len(self.frequency_VR)))

        self.dyn_rho_A = np.ascontiguousarray(np.zeros((N, params.timeDIM_VR)), dtype=np.complex)
        self.dyn_rho_B = np.ascontiguousarray(np.zeros((N, params.timeDIM_VR)), dtype=np.complex)

    def create_molecules(self, molA, molB):
        molA.nDIM = len(self.energies_A)
        molA.energies = self.energies_A.ctypes.data_as(POINTER(c_double))
        molA.gamma_decay = self.gamma_decay.ctypes.data_as(POINTER(c_double))
        molA.gamma_pure_dephasing = self.gamma_pure_dephasingA.ctypes.data_as(POINTER(c_double))
        molA.mu = self.mu.ctypes.data_as(POINTER(c_complex))
        molA.rho = self.rhoA.ctypes.data_as(POINTER(c_complex))
        molA.rho_0 = self.rho_0.ctypes.data_as(POINTER(c_complex))
        molA.abs_spectra = self.abs_spectraA.ctypes.data_as(POINTER(c_double))
        molA.ems_spectra = self.ems_spectraA.ctypes.data_as(POINTER(c_double))
        molA.vib_spectra = self.vib_spectraA.ctypes.data_as(POINTER(c_double))
        molA.Raman_spectra = self.Raman_spectraA.ctypes.data_as(POINTER(c_double))
        molA.dyn_rho = self.dyn_rho_A.ctypes.data_as(POINTER(c_complex))

        molB.nDIM = len(self.energies_B)
        molB.energies = self.energies_B.ctypes.data_as(POINTER(c_double))
        molB.gamma_decay = self.gamma_decay.ctypes.data_as(POINTER(c_double))
        molB.gamma_pure_dephasing = self.gamma_pure_dephasingB.ctypes.data_as(POINTER(c_double))
        molB.mu = self.mu.ctypes.data_as(POINTER(c_complex))
        molB.rho = self.rhoB.ctypes.data_as(POINTER(c_complex))
        molB.rho_0 = self.rho_0.ctypes.data_as(POINTER(c_complex))
        molB.abs_spectra = self.abs_spectraB.ctypes.data_as(POINTER(c_double))
        molB.ems_spectra = self.ems_spectraB.ctypes.data_as(POINTER(c_double))
        molB.vib_spectra = self.vib_spectraB.ctypes.data_as(POINTER(c_double))
        molB.Raman_spectra = self.Raman_spectraB.ctypes.data_as(POINTER(c_double))
        molB.dyn_rho = self.dyn_rho_B.ctypes.data_as(POINTER(c_complex))

    def create_parameters_spectra(self, spectra_params, params):
        spectra_params.rho_0_A = params.rho_0_A.ctypes.data_as(POINTER(c_complex))
        spectra_params.rho_0_E = params.rho_0_E.ctypes.data_as(POINTER(c_complex))
        spectra_params.time_AE = self.time_AE.ctypes.data_as(POINTER(c_double))
        spectra_params.time_VR = self.time_VR.ctypes.data_as(POINTER(c_double))
        spectra_params.frequency_A = self.frequency_A.ctypes.data_as(POINTER(c_double))
        spectra_params.frequency_E = self.frequency_E.ctypes.data_as(POINTER(c_double))
        spectra_params.frequency_VR = self.frequency_VR.ctypes.data_as(POINTER(c_double))

        spectra_params.field_amp_AE = params.field_amp_AE
        spectra_params.field_amp_VR = params.field_amp_VR
        spectra_params.omega_R = params.omega_R

        spectra_params.nDIM = len(self.energies_A)
        spectra_params.nEXC = params.nEXC

        spectra_params.timeDIM_AE = len(self.time_AE)
        spectra_params.timeDIM_VR = len(self.time_VR)

        spectra_params.freqDIM_A = len(self.frequency_A)
        spectra_params.freqDIM_E = len(self.frequency_E)
        spectra_params.freqDIM_VR = len(self.frequency_VR)

        spectra_params.field_A = self.field_A.ctypes.data_as(POINTER(c_complex))
        spectra_params.field_E = self.field_E.ctypes.data_as(POINTER(c_complex))
        spectra_params.field_V = self.field_V.ctypes.data_as(POINTER(c_complex))
        spectra_params.field_R = self.field_R.ctypes.data_as(POINTER(c_complex))

        spectra_params.omega_v1 = params.omega_v1
        spectra_params.omega_v2 = params.omega_v2
        spectra_params.omega_v3 = params.omega_v3
        spectra_params.omega_v4 = params.omega_v4
        spectra_params.omega_e1 = params.omega_e1

    def calculate_spectra(self, params):
        molA = Molecule()
        molB = Molecule()
        self.create_molecules(molA, molB)
        params_spectra = Parameters_Spectra()
        self.create_parameters_spectra(params_spectra, params)
        CalculateSpectra(molA, molB, params_spectra)
        return


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import time

    np.set_printoptions(precision=4)
    energy_factor = 1. / 27.211385
    time_factor = .02418884 / 1000

    energies_A = np.array((0.000, 0.08233, 0.09832, 0.16304, 0.20209, 1.7679256, 1.85871, 1.87855, 1.96783, 2.02991)) * energy_factor
    energies_B = np.array((0.000, 0.08313, 0.09931, 0.15907, 0.19924, 1.5879256, 1.66871, 1.67546, 1.76783, 1.82991)) * energy_factor
    N = len(energies_A)
    N_vib = N - 5
    N_exc = N - N_vib
    rho_0_ems = np.zeros((N, N), dtype=np.complex)
    rho_0_ems[N_vib, N_vib] = 1. + 0j
    rho_0_abs = np.zeros((N, N), dtype=np.complex)
    rho_0_abs[0, 0] = 1. + 0j

    mu = 4.97738 * np.ones_like(rho_0_abs)
    np.fill_diagonal(mu, 0j)
    population_decay = 2.418884e-8
    electronic_dephasingA = 2.7 * 2.418884e-4
    electronic_dephasingB = 5.0 * 2.418884e-4
    vibrational_dephasing = 0.1 * 2.418884e-5

    gamma_decay = np.ones((N, N)) * population_decay
    np.fill_diagonal(gamma_decay, 0.0)
    gamma_decay = np.tril(gamma_decay)

    gamma_pure_dephasingA = np.ones_like(gamma_decay) * vibrational_dephasing
    np.fill_diagonal(gamma_pure_dephasingA, 0.0)

    for i in range(N_vib):
        for j in range(N_vib, N):
            gamma_pure_dephasingA[i, j] = electronic_dephasingA
            gamma_pure_dephasingA[j, i] = electronic_dephasingA

    gamma_pure_dephasingA[5, 4] = electronic_dephasingA * 0.65
    gamma_pure_dephasingA[4, 5] = electronic_dephasingA * 0.65
    gamma_pure_dephasingA[0, 9] = electronic_dephasingA * 0.65
    gamma_pure_dephasingA[9, 0] = electronic_dephasingA * 0.65

    gamma_pure_dephasingA[5, 3] = electronic_dephasingA * 0.70
    gamma_pure_dephasingA[3, 5] = electronic_dephasingA * 0.70
    gamma_pure_dephasingA[0, 8] = electronic_dephasingA * 0.70
    gamma_pure_dephasingA[8, 0] = electronic_dephasingA * 0.70

    gamma_pure_dephasingA[5, 2] = electronic_dephasingA * 0.20
    gamma_pure_dephasingA[2, 5] = electronic_dephasingA * 0.20
    gamma_pure_dephasingA[0, 7] = electronic_dephasingA * 0.20
    gamma_pure_dephasingA[7, 0] = electronic_dephasingA * 0.20

    gamma_pure_dephasingA[5, 1] = electronic_dephasingA * 0.18
    gamma_pure_dephasingA[1, 5] = electronic_dephasingA * 0.18
    gamma_pure_dephasingA[0, 6] = electronic_dephasingA * 0.18
    gamma_pure_dephasingA[6, 0] = electronic_dephasingA * 0.18

    gamma_pure_dephasingA[5, 0] = electronic_dephasingA * 0.60
    gamma_pure_dephasingA[0, 5] = electronic_dephasingA * 0.60
    mu[5, 0] *= 0.10
    mu[0, 5] *= 0.10

    gamma_pure_dephasingB = np.ones_like(gamma_decay) * vibrational_dephasing
    np.fill_diagonal(gamma_pure_dephasingB, 0.0)
    for i in range(N_vib):
        for j in range(N_vib, N):
            gamma_pure_dephasingB[i, j] = electronic_dephasingB
            gamma_pure_dephasingB[j, i] = electronic_dephasingB

    gamma_pure_dephasingB[5, 4] = electronic_dephasingB * 0.65
    gamma_pure_dephasingB[4, 5] = electronic_dephasingB * 0.65
    gamma_pure_dephasingB[0, 9] = electronic_dephasingB * 0.65
    gamma_pure_dephasingB[9, 0] = electronic_dephasingB * 0.65

    gamma_pure_dephasingB[5, 3] = electronic_dephasingB * 0.70
    gamma_pure_dephasingB[3, 5] = electronic_dephasingB * 0.70
    gamma_pure_dephasingB[0, 8] = electronic_dephasingB * 0.70
    gamma_pure_dephasingB[8, 0] = electronic_dephasingB * 0.70

    gamma_pure_dephasingB[5, 2] = electronic_dephasingB * 0.20
    gamma_pure_dephasingB[2, 5] = electronic_dephasingB * 0.20
    gamma_pure_dephasingB[0, 7] = electronic_dephasingB * 0.20
    gamma_pure_dephasingB[7, 0] = electronic_dephasingB * 0.20

    gamma_pure_dephasingB[5, 1] = electronic_dephasingB * 0.18
    gamma_pure_dephasingB[1, 5] = electronic_dephasingB * 0.18
    gamma_pure_dephasingB[0, 6] = electronic_dephasingB * 0.18
    gamma_pure_dephasingB[6, 0] = electronic_dephasingB * 0.18

    gamma_pure_dephasingB[5, 0] = electronic_dephasingB * 0.60
    gamma_pure_dephasingB[0, 5] = electronic_dephasingB * 0.60

    np.set_printoptions(precision=2)

    print(gamma_pure_dephasingA)

    params = ADict(
        energy_factor=energy_factor,
        time_factor=time_factor,
        rho_0_A=rho_0_abs,
        rho_0_E=rho_0_ems,

        timeDIM_AE=1000,
        timeAMP_AE=2000,

        timeDIM_VR=5000,
        timeAMP_VR=50000,

        frequencyDIM_AE=250,
        frequencyMIN_A=1.5*energy_factor,
        frequencyMAX_A=2.7*energy_factor,

        frequencyMIN_E=1.2 * energy_factor,
        frequencyMAX_E=2.3 * energy_factor,

        frequencyDIM_VR=250,
        frequencyMIN_VR=0.075*energy_factor,
        frequencyMAX_VR=0.21*energy_factor,

        field_amp_AE=0.000003,
        field_amp_VR=0.000014,

        omega_R=0.5*energy_factor,
        nEXC=N_exc,

        omega_v1=energies_A[1],
        omega_v2=energies_A[2],
        omega_v3=energies_A[3],
        omega_v4=energies_A[4],

        omega_e1=energies_A[5]-energies_A[4]
    )

    FourLevels = dict(
        energies_A=energies_A,
        energies_B=energies_B,
        gamma_decay=gamma_decay,
        gamma_pure_dephasingA=gamma_pure_dephasingA,
        gamma_pure_dephasingB=gamma_pure_dephasingB,
        mu=mu,
    )

    def render_ticks(axes):
        axes.get_xaxis().set_tick_params(which='both', direction='in', width=1, labelrotation=0, labelsize='large')
        axes.get_yaxis().set_tick_params(which='both', direction='in', width=1, labelcolor='r', labelsize='large')
        axes.get_xaxis().set_ticks_position('both')
        axes.get_yaxis().set_ticks_position('both')
        axes.grid()


    molecules = RamanControl(params, **FourLevels)

    start = time.time()
    molecules.calculate_spectra(params)

    print(time.time() - start)

    # fig, axes = plt.subplots(nrows=3, ncols=1)
    # axes[0].plot(molecules.time_VR, molecules.field_R.real)
    # axes[0].plot(molecules.time_AE, molecules.field_A.real)
    # axes[1].plot(energy_factor * 1239.84 / molecules.frequency_A, molecules.abs_spectraA)
    # axes[1].plot(energy_factor * 1239.84 / molecules.frequency_A, molecules.abs_spectraB)
    # axes[2].plot(energy_factor * 1239.84 / molecules.frequency_VR, molecules.Raman_spectraA)
    # axes[2].plot(energy_factor * 1239.84 / molecules.frequency_VR, molecules.Raman_spectraB)
    # plt.show()

    fig, axes = plt.subplots(nrows=3, ncols=1)
    axes[0].plot(molecules.time_VR, molecules.field_R.real)

    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[0], label='g1')
    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[1], label='g2')
    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[2], label='g3')
    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[3], label='g4')
    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[4], label='g5')
    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[5])
    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[6])
    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[7])
    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[8])
    axes[1].plot(molecules.time_VR, molecules.dyn_rho_A[9])

    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[0], label='g1')
    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[1], label='g2')
    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[2], label='g3')
    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[3], label='g4')
    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[4], label='g5')
    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[5])
    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[6])
    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[7])
    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[8])
    axes[2].plot(molecules.time_VR, molecules.dyn_rho_B[9])

    render_ticks(axes[0])
    render_ticks(axes[1])
    render_ticks(axes[2])

    axes[1].legend()
    axes[2].legend()

    print(molecules.rhoA.diagonal())
    print()
    print(molecules.rhoB.diagonal())

    plt.show()
