import numpy as np
from types import MethodType, FunctionType
from RamanAssistedControl_wrapper import *
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

        self.time = np.linspace(-params.timeAMP, params.timeAMP, params.timeDIM)
        self.time_spectra_abs_ems = np.linspace(-params.timeAMP_spectra_abs_ems, params.timeAMP_spectra_abs_ems, params.timeDIM_spectra_abs_ems)
        self.time_spectra_vib = np.linspace(-params.timeAMP_spectra_vib, params.timeAMP_spectra_vib, params.timeDIM_spectra_vib)

        self.frequency_abs = 1./np.linspace(1./params.frequencyMAX_abs, 1./params.frequencyMIN_abs, params.frequencyDIM_abs)
        self.frequency_ems = 1./np.linspace(1./params.frequencyMAX_ems, 1./params.frequencyMIN_ems, params.frequencyDIM_ems)
        self.frequency_vib = 1./np.linspace(1./params.frequencyMAX_vib, 1./params.frequencyMIN_vib, params.frequencyDIM_vib)

        self.field_t = np.empty(params.timeDIM, dtype=np.complex)
        self.field_abs = np.empty(params.timeDIM_spectra_abs_ems, dtype=np.complex)
        self.field_ems = np.empty(params.timeDIM_spectra_abs_ems, dtype=np.complex)
        self.field_vib = np.empty(params.timeDIM_spectra_vib, dtype=np.complex)

        self.gamma_decay = np.ascontiguousarray(self.gamma_decay)
        self.gamma_pure_dephasingA = np.ascontiguousarray(self.gamma_pure_dephasingA)
        self.gamma_pure_dephasingB = np.ascontiguousarray(self.gamma_pure_dephasingB)
        self.mu = np.ascontiguousarray(self.mu)
        self.rho_0 = np.ascontiguousarray(params.rho_0_abs)
        self.rhoA = np.ascontiguousarray(params.rho_0_abs.copy())
        self.rhoB = np.ascontiguousarray(params.rho_0_abs.copy())
        self.energies_A = np.ascontiguousarray(self.energies_A)
        self.energies_B = np.ascontiguousarray(self.energies_B)

        N = len(self.energies_A)
        self.N_exc = N - params.N_vib
        self.dyn_rhoA = np.ascontiguousarray(np.zeros((N+1, params.timeDIM), dtype=np.complex))
        self.dyn_rhoB = np.ascontiguousarray(np.zeros((N+1, params.timeDIM), dtype=np.complex))

        self.abs_spectraA = np.ascontiguousarray(np.zeros(len(self.frequency_abs)))
        self.abs_spectraB = np.ascontiguousarray(np.zeros(len(self.frequency_abs)))
        self.ems_spectraA = np.ascontiguousarray(np.zeros(len(self.frequency_ems)))
        self.ems_spectraB = np.ascontiguousarray(np.zeros(len(self.frequency_ems)))
        self.vib_spectraA = np.ascontiguousarray(np.zeros(len(self.frequency_vib)))
        self.vib_spectraB = np.ascontiguousarray(np.zeros(len(self.frequency_vib)))

    def create_molecules(self, molA, molB):
        molA.nDIM = len(self.energies_A)
        molA.energies = self.energies_A.ctypes.data_as(POINTER(c_double))
        molA.gamma_decay = self.gamma_decay.ctypes.data_as(POINTER(c_double))
        molA.gamma_pure_dephasing = self.gamma_pure_dephasingA.ctypes.data_as(POINTER(c_double))
        molA.mu = self.mu.ctypes.data_as(POINTER(c_complex))
        molA.rho = self.rhoA.ctypes.data_as(POINTER(c_complex))
        molA.dyn_rho = self.dyn_rhoA.ctypes.data_as(POINTER(c_complex))
        molA.rho_0 = self.rho_0.ctypes.data_as(POINTER(c_complex))
        molA.abs_spectra = self.abs_spectraA.ctypes.data_as(POINTER(c_double))
        molA.ems_spectra = self.ems_spectraA.ctypes.data_as(POINTER(c_double))
        molA.vib_spectra = self.vib_spectraA.ctypes.data_as(POINTER(c_double))

        molB.nDIM = len(self.energies_A)
        molB.energies = self.energies_B.ctypes.data_as(POINTER(c_double))
        molB.gamma_decay = self.gamma_decay.ctypes.data_as(POINTER(c_double))
        molB.gamma_pure_dephasing = self.gamma_pure_dephasingB.ctypes.data_as(POINTER(c_double))
        molB.mu = self.mu.ctypes.data_as(POINTER(c_complex))
        molB.rho = self.rhoB.ctypes.data_as(POINTER(c_complex))
        molB.dyn_rho = self.dyn_rhoB.ctypes.data_as(POINTER(c_complex))
        molB.rho_0 = self.rho_0.ctypes.data_as(POINTER(c_complex))
        molB.abs_spectra = self.abs_spectraB.ctypes.data_as(POINTER(c_double))
        molB.ems_spectra = self.ems_spectraB.ctypes.data_as(POINTER(c_double))
        molB.vib_spectra = self.vib_spectraB.ctypes.data_as(POINTER(c_double))

    def create_parameters(self, func_params, params):
        func_params.time = self.time.ctypes.data_as(POINTER(c_double))
        func_params.rho_0 = self.rho_0.ctypes.data_as(POINTER(c_complex))

        t0_R = params.t0_R * params.timeAMP
        t0_EE = params.t0_EE * params.timeAMP
        width_R = params.timeAMP / params.width_R
        width_EE = params.timeAMP / params.width_EE

        func_params.A_R = params.A_R
        func_params.width_R = width_R
        func_params.t0_R = t0_R

        func_params.A_EE = params.A_EE
        func_params.width_EE = width_EE
        func_params.t0_EE = t0_EE

        func_params.w_R = params.w_R
        func_params.w_v1 = params.w_v1
        func_params.w_v2 = params.w_v2
        func_params.w_EE = params.w_EE

        func_params.nDIM = len(self.energies_A)
        func_params.timeDIM = len(self.time)

        func_params.field_out = self.field_t.ctypes.data_as(POINTER(c_complex))

        func_params.lower_bounds = params.lower_bounds.ctypes.data_as(POINTER(c_double))
        func_params.upper_bounds = params.upper_bounds.ctypes.data_as(POINTER(c_double))
        func_params.guess = guess.ctypes.data_as(POINTER(c_double))

        func_params.MAX_EVAL = params.MAX_EVAL

    def create_parameters_abs_spectra(self, spectra_params, params):
        spectra_params.rho_0_abs = params.rho_0_abs.ctypes.data_as(POINTER(c_complex))
        spectra_params.rho_0_ems = params.rho_0_ems.ctypes.data_as(POINTER(c_complex))
        spectra_params.time_spectra_abs_ems = self.time_spectra_abs_ems.ctypes.data_as(POINTER(c_double))
        spectra_params.frequency_abs = self.frequency_abs.ctypes.data_as(POINTER(c_double))
        spectra_params.frequency_ems = self.frequency_ems.ctypes.data_as(POINTER(c_double))

        spectra_params.A_abs_ems = params.A_abs_ems

        spectra_params.nDIM = len(self.energies_A)
        spectra_params.nEXC = self.N_exc
        spectra_params.timeDIM_abs_ems = len(self.time_spectra_abs_ems)
        spectra_params.freqDIM_abs = len(self.frequency_abs)
        spectra_params.freqDIM_ems = len(self.frequency_ems)

        spectra_params.field_abs = self.field_abs.ctypes.data_as(POINTER(c_complex))
        spectra_params.field_ems = self.field_ems.ctypes.data_as(POINTER(c_complex))

    def create_parameters_vib_spectra(self, spectra_params, params):
        spectra_params.time_spectra_vib = self.time_spectra_vib.ctypes.data_as(POINTER(c_double))
        spectra_params.frequency_vib = self.frequency_vib.ctypes.data_as(POINTER(c_double))

        spectra_params.A_vib = params.A_vib
        spectra_params.w_R = params.vib_R

        spectra_params.nDIM = len(self.energies_A)
        spectra_params.timeDIM_vib = len(self.time_spectra_vib)
        spectra_params.freqDIM_vib = len(self.frequency_vib)

        spectra_params.field_vib = self.field_vib.ctypes.data_as(POINTER(c_complex))

        print(self.field_vib.shape)

    def call_raman_control_func(self, params):
        molA = Molecule()
        molB = Molecule()
        self.create_molecules(molA, molB)
        func_params = Parameters()
        self.create_parameters(func_params, params)
        RamanControlFunction(molA, molB, func_params)
        return

    def calculate_spectra(self, params):
        molA = Molecule()
        molB = Molecule()
        self.create_molecules(molA, molB)
        abs_spec_params = Parameters_AbsEmsSpectra()
        vib_spec_params = Parameters_VibSpectra()
        self.create_parameters_abs_spectra(abs_spec_params, params)
        self.create_parameters_vib_spectra(vib_spec_params, params)
        CalculateSpectra(molA, molB, abs_spec_params, vib_spec_params)
        return


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import time
    # from itertools import *

    np.set_printoptions(precision=4)
    energy_factor = 1. / 27.211385
    time_factor = .02418884 / 1000

    energies_A = np.array((0.000, 0.08233, 0.09832, 0.16304, 0.20209, 1.7679256, 1.85871, 1.87855, 1.96783, 2.02991)) * energy_factor
    energies_B = np.array((0.000, 0.08313, 0.09931, 0.15907, 0.19924, 1.7712000, 1.80871, 1.82855, 1.86783, 1.92991)) * energy_factor

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
    electronic_dephasingB = 4.0 * 2.418884e-4
    vibrational_dephasing = 0.05 * 2.418884e-5

    gamma_decay = np.ones((N, N)) * population_decay
    np.fill_diagonal(gamma_decay, 0.0)
    gamma_decay = np.tril(gamma_decay)

    gamma_pure_dephasingA = np.ones_like(gamma_decay) * vibrational_dephasing
    np.fill_diagonal(gamma_pure_dephasingA, 0.0)

    for i in range(N_vib):
        for j in range(N_vib, N):
            gamma_pure_dephasingA[i, j] = electronic_dephasingA
            gamma_pure_dephasingA[j, i] = electronic_dephasingA

    gamma_pure_dephasingA[5, 4] = electronic_dephasingA*0.65
    gamma_pure_dephasingA[4, 5] = electronic_dephasingA*0.65
    gamma_pure_dephasingA[0, 9] = electronic_dephasingA*0.65
    gamma_pure_dephasingA[9, 0] = electronic_dephasingA*0.65

    gamma_pure_dephasingA[5, 3] = electronic_dephasingA*0.70
    gamma_pure_dephasingA[3, 5] = electronic_dephasingA*0.70
    gamma_pure_dephasingA[0, 8] = electronic_dephasingA*0.70
    gamma_pure_dephasingA[8, 0] = electronic_dephasingA*0.70

    gamma_pure_dephasingA[5, 2] = electronic_dephasingA*0.20
    gamma_pure_dephasingA[2, 5] = electronic_dephasingA*0.20
    gamma_pure_dephasingA[0, 7] = electronic_dephasingA*0.20
    gamma_pure_dephasingA[7, 0] = electronic_dephasingA*0.20

    gamma_pure_dephasingA[5, 1] = electronic_dephasingA*0.18
    gamma_pure_dephasingA[1, 5] = electronic_dephasingA*0.18
    gamma_pure_dephasingA[0, 6] = electronic_dephasingA*0.18
    gamma_pure_dephasingA[6, 0] = electronic_dephasingA*0.18

    gamma_pure_dephasingA[5, 0] = electronic_dephasingA*0.60
    gamma_pure_dephasingA[0, 5] = electronic_dephasingA*0.60
    mu[5, 0] *= 0.10
    mu[0, 5] *= 0.10

    gamma_pure_dephasingB = np.ones_like(gamma_decay) * vibrational_dephasing
    np.fill_diagonal(gamma_pure_dephasingB, 0.0)
    for i in range(N_vib):
        for j in range(N_vib, N):
            if i == 1 or i == 2:
                gamma_pure_dephasingB[i, j] = electronic_dephasingB*0.35
                gamma_pure_dephasingB[j, i] = electronic_dephasingB*0.35
            else:
                gamma_pure_dephasingB[i, j] = electronic_dephasingB*0.95
                gamma_pure_dephasingB[j, i] = electronic_dephasingB*0.95

    np.set_printoptions(precision=2)

    # print(gamma_decay)
    print(gamma_pure_dephasingA)

    lower_bounds = np.asarray([0.000005, 3.5, 0.9 * energies_A[2], 0.9 * energies_A[3], 0.35 * energy_factor])
    upper_bounds = np.asarray([0.000045, 10.0, 1.1 * energies_A[2], 1.1 * energies_A[3], 0.5 * energy_factor])
    guess = np.asarray([0.000709642, 3.52037, 0.00547927, 0.0069986, 0.0179549])

    params = ADict(
        energy_factor=energy_factor,
        time_factor=time_factor,
        timeDIM=120000,
        timeAMP=60000,

        A_R=0.000576595,
        width_R=5.017,
        t0_R=0.0,

        A_EE=0.000366972,
        width_EE=19.32,
        t0_EE=0.55,

        w_R=0.6 * energy_factor,
        w_v1=energies_A[2],
        w_v2=energies_A[3],
        w_EE=(energies_A[4] - energies_A[1]),
        rho_0_abs=rho_0_abs,
        rho_0_ems=rho_0_ems,

        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
        guess=guess,

        MAX_EVAL=20,

        timeDIM_spectra_abs_ems=1000,
        timeAMP_spectra_abs_ems=2000,

        timeDIM_spectra_vib=10000,
        timeAMP_spectra_vib=200000,

        frequencyDIM_abs=250,
        frequencyMIN_abs=1.5*energy_factor,
        frequencyMAX_abs=2.7*energy_factor,

        frequencyDIM_ems=250,
        frequencyMIN_ems=1.2 * energy_factor,
        frequencyMAX_ems=2.3 * energy_factor,

        frequencyDIM_vib=250,
        frequencyMIN_vib=0.05*energy_factor,
        frequencyMAX_vib=0.24*energy_factor,

        A_abs_ems=0.000003,
        A_vib=0.000005,

        vib_R=0.5*energy_factor,
        N_vib=N_vib
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

    start = time.time()
    molecules = RamanControl(params, **FourLevels)
    molecules.calculate_spectra(params)
    end_spectra = time.time()
    print("Time to calculate spectra: ", end_spectra - start)
    print()

    fig, axes = plt.subplots(nrows=4, ncols=1)
    # gs1 = gridspec.GridSpec(2, 1)
    # gs2 = gridspec.GridSpec(2, 1)
    # for ax in axes[1, :]:
    #     ax.remove()
    # for ax in axes[0, :]:
    #     ax.remove()
    # axes[0, 0] = fig.add_subplot(gs1[0, 0:])
    # axes_spectra = fig.add_subplot(gs2[1, 0:])

    axes[0].set_title(
        'Field for Electronic Spectra \n $\\tau_E$= {} fs'
            .format(int(1e3*time_factor/electronic_dephasingA)))
    axes[0].plot(molecules.time_spectra_abs_ems * time_factor, molecules.field_abs.real, 'r')
    axes[0].plot(molecules.time_spectra_abs_ems * time_factor, molecules.field_ems.real, 'b')
    render_ticks(axes[0])
    axes[0].set_xlabel("Time (in ps)")

    axes[2].set_title(
        'Field for Vibrational Spectra \n $\\tau_E$= {} fs'
            .format(int(1e3 * time_factor / electronic_dephasingA)))
    axes[2].plot(molecules.time_spectra_vib * time_factor, molecules.field_vib.real, 'k')
    render_ticks(axes[2])
    axes[2].set_xlabel("Time (in ps)")

    axes[1].set_title("Absorption & Emission spectra")
    axes[1].plot(1239.84 / (molecules.frequency_abs / energy_factor), molecules.abs_spectraA, 'r')
    axes[1].plot(1239.84 / (molecules.frequency_ems / energy_factor), molecules.ems_spectraA, 'b')
    axes[1].set_xlim(450., 850.)

    axes[3].set_title("Vibrational spectra")
    axes[3].plot((molecules.frequency_vib / energy_factor), molecules.vib_spectraA, 'k')

    axes[0].set_xlabel("Time (in ps)")
    axes[0].yaxis.tick_left()
    axes[0].yaxis.set_label_position("left")
    axes[1].yaxis.tick_right()
    axes[1].yaxis.set_label_position("right")
    axes[2].yaxis.tick_right()
    axes[2].yaxis.set_label_position("right")

    axes[0].ticklabel_format(scilimits=(-2, 2))
    axes[2].ticklabel_format(scilimits=(-4, 4))
    render_ticks(axes[1])
    render_ticks(axes[3])

    fig.subplots_adjust(left=0.32, bottom=None, right=0.68, top=0.8, wspace=0.025, hspace=0.55)

    plt.show()

