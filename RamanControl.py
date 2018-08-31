import numpy as np
from types import MethodType, FunctionType
from Wrapper import *


class RamanControl:

    def __init__(self, **kwargs):
        """
        __init__ function call to initialize variables from the
        parameters for the class instance provided in __main__ and
        add new variables for use in other functions in this class.
        """

        for name, value in kwargs.items():
            if isinstance(value, FunctionType):
                setattr(self, name, MethodType(value, self))
            else:
                setattr(self, name, value)

        self.dyn_rho_A = np.empty((int(self.nDIM * (self.nDIM + 3) / 2), self.timeDIM), dtype=np.complex)
        self.dyn_rho_B = np.empty_like(self.dyn_rho_A, dtype=np.complex)

        self.rho_A = np.empty((self.nDIM, self.nDIM), dtype=np.complex)
        self.rho_B = np.empty_like(self.rho_A, dtype=np.complex)

        self.time = np.linspace(-self.timeAMP, self.timeAMP, self.timeDIM, endpoint=False)
        self.time_spectra = np.linspace(-self.timeAMP_spectra, self.timeAMP_spectra, self.timeDIM_spectra, endpoint=False)
        self.field = np.empty(self.timeDIM)
        self.field_spectra = np.empty(self.timeDIM_spectra, dtype=np.complex)
        self.freq_spectra = np.linspace(self.freq_min, self.freq_max, self.freqDIM)
        self.spectra_A = np.empty(self.freqDIM, dtype=np.complex)
        self.spectra_B = np.empty(self.freqDIM, dtype=np.complex)

    def create_molecule(self, MolA, MolB):
        """
        Initializes two instance of the Structure "Molecule"
        :param MolA: Molecule A of type Structure(Molecule)
        :param MolB: Molecule B of type Structure(Molecule)
        :return: None
        """

        MolA.energies = self.energies_A.ctypes.data_as(POINTER(c_double))
        MolB.energies = self.energies_B.ctypes.data_as(POINTER(c_double))
        MolA.gamma_decay = self.gamma_decay_A.ctypes.data_as(POINTER(c_double))
        MolB.gamma_decay = self.gamma_decay_B.ctypes.data_as(POINTER(c_double))
        MolA.gamma_pcd = self.gamma_pcd_A.ctypes.data_as(POINTER(c_double))
        MolB.gamma_pcd = self.gamma_pcd_B.ctypes.data_as(POINTER(c_double))
        MolA.rho = self.rho_A.ctypes.data_as(POINTER(c_complex))
        MolB.rho = self.rho_B.ctypes.data_as(POINTER(c_complex))
        MolA.dyn_rho = self.dyn_rho_A.ctypes.data_as(POINTER(c_complex))
        MolB.dyn_rho = self.dyn_rho_B.ctypes.data_as(POINTER(c_complex))
        MolA.spectra = self.spectra_A.ctypes.data_as(POINTER(c_complex))
        MolB.spectra = self.spectra_B.ctypes.data_as(POINTER(c_complex))

        return

    def create_parameters(self, func_params):
        """
        Initializes an instance of the Structure "Parameters"
        :param func_params: Parameter Structure instance to be updated
        :return: None
        """

        func_params.time = self.time.ctypes.data_as(POINTER(c_double))
        func_params.rho_0 = self.rho_0.ctypes.data_as(POINTER(c_complex))
        func_params.mu = self.mu.ctypes.data_as(POINTER(c_complex))
        func_params.A_R = self.A_R
        func_params.width_R = self.timeAMP / self.width_R
        func_params.t0_R = self.timeAMP * self.t0_R
        func_params.A_EE = self.A_EE
        func_params.width_EE = self.timeAMP / self.width_EE
        func_params.t0_EE = self.timeAMP * self.t0_EE
        func_params.w_R = self.w_R
        func_params.w_v = self.w_v
        func_params.w_EE = self.w_EE
        func_params.nDIM = self.nDIM
        func_params.timeDIM = self.timeDIM
        func_params.field_out = self.field.ctypes.data_as(POINTER(c_double))

        return

    def create_spectra_parameters(self, spec_params):
        """
        Initializes an instance of the Structure "Parameters"
        :param func_params: Parameter Structure instance to be updated
        :return: None
        """

        spec_params.field_spectra = self.field_spectra.ctypes.data_as(POINTER(c_complex))
        spec_params.freq_spectra = self.freq_spectra.ctypes.data_as(POINTER(c_double))
        spec_params.freqDIM = self.freqDIM
        spec_params.w_R = self.w_R
        spec_params.A_S = self.A_S
        spec_params.width_S = self.timeAMP_spectra / self.width_S
        spec_params.time_spectra = self.time_spectra.ctypes.data_as(POINTER(c_double))
        spec_params.timeDIM_spectra = self.timeDIM_spectra
        spec_params.nDIM = self.nDIM

        return

    def call_raman_optimization(self):
        """

        :return:
        """
        MolA = Molecule()
        MolB = Molecule()
        func_params = Parameters()
        spec_params = SpectraParams()
        self.create_molecule(MolA, MolB)
        self.create_parameters(func_params)
        self.create_spectra_parameters(spec_params)
        RamanControlFunction(MolA, MolB, func_params, spec_params)


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    np.set_printoptions(precision=4)

    energy_factor = 1. / 27.211385
    time_factor = .02418
    timeAMP = 20000.
                                                                                ##################################################################
                                                                                #                                                                #
    energies_A = np.array((0.000, 0.07439, 1.94655, 2.02094)) * energy_factor   #   Energies in a.u. for molecule A                              #
    energies_B = np.array((0.000, 0.09639, 1.92655, 2.00094)) * energy_factor   #   Energies in a.u. for molecule B                              #
                                                                                #                                                                #
    N = len(energies_A)                                                         #                                                                #
    rho_0 = np.zeros((N, N), dtype=np.complex)                                  #                                                                #
    rho_0[0, 0] = 1. + 0j                                                       #   Only ground state populated                                  #
    mu = 4.97738 * np.ones_like(rho_0)                                          #   Dipole moment in a.u. (same for all allowed transition)      #
    np.fill_diagonal(mu, 0j)                                                    #   \mu_{i, i} = 0.                                              #
                                                                                #                                                                #
    population_decay = 2.418884e-8                                              #   All population relaxation lifetimes equal 1 ns               #
    electronic_dephasing = 3. * 2.418884e-4                                     #   All electronic dephasing 50 fs                               #
    vibrational_dephasing = 1.422872e-5                                         #   All vibrational dephasing 1 ps                               #
                                                                                #                                                                #
    gamma_decay = np.ones((N, N)) * population_decay                            #                                                                #
    np.fill_diagonal(gamma_decay, 0.0)                                          #   No decay for diagonal elements                               #
    gamma_decay = np.tril(gamma_decay)                                          #   Decay non-zero only from higher to lower energy levels       #
                                                                                #                                                                #
    gamma_pure_dephasing = np.ones_like(gamma_decay) * electronic_dephasing     #                                                                #
    np.fill_diagonal(gamma_pure_dephasing, 0.0)                                 #   No dephasing for population terms                            #
    gamma_pure_dephasing[1, 0] = vibrational_dephasing                          #                                                                #
    gamma_pure_dephasing[3, 2] = vibrational_dephasing                          ##################################################################
    gamma_pure_dephasing[0, 1] = vibrational_dephasing
    gamma_pure_dephasing[2, 3] = vibrational_dephasing
    gamma_pure_dephasing[2, 0] *= 1.
    gamma_pure_dephasing[0, 2] *= 1.
    gamma_pure_dephasing[3, 0] *= 2.
    gamma_pure_dephasing[0, 3] *= 2.
    gamma_pure_dephasing[1, 3] *= 1.
    gamma_pure_dephasing[3, 1] *= 1.

    lower_bounds = np.asarray([0.00020, timeAMP / 25.0, 0.9 * energies_A[1], 1.2 * energy_factor])
    upper_bounds = np.asarray([0.00070, timeAMP / 6.0, 1.1 * energies_A[1], 1.6 * energy_factor])

    guess = np.asarray([np.random.uniform(lower_bounds[i], upper_bounds[i]) for i in range(len(lower_bounds))])

    params = dict(
        energy_factor=energy_factor,
        time_factor=time_factor,
        timeDIM=40000,
        timeAMP=timeAMP,
        nDIM=N,

        A_R=0.000576595,
        width_R=5.017,
        t0_R=-0.35,

        A_EE=0.000366972,
        width_EE=19.32,
        t0_EE=0.55,

        w_R=1.4 * energy_factor,
        w_v=energies_A[1],
        w_EE=(energies_A[2] - energies_A[1]),
        rho_0=rho_0,

        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
        guess=guess,

        energies_A=energies_A,
        energies_B=energies_B,
        gamma_decay_A=gamma_decay,
        gamma_decay_B=gamma_decay,
        gamma_pcd_A=gamma_pure_dephasing,
        gamma_pcd_B=gamma_pure_dephasing,
        mu=mu,

        freqDIM=100,
        freq_min=.05 * energy_factor,
        freq_max=.12 * energy_factor,
        A_S=0.00005,
        width_S=4.,
        timeDIM_spectra=40000,
        timeAMP_spectra=20000
    )

    gamma = 9. * gamma_pure_dephasing.copy()
    for i in range(4):
        for j in range(4):
            for k in range(4):
                if k > i:
                    gamma[i, j] += gamma_decay[k, i]
                if k > j:
                    gamma[i, j] += gamma_decay[k, j]
            gamma[i, j] *= 0.5

    print(gamma)
    N = 100
    omega = np.linspace(params['freq_min'], params['freq_max'], params['freqDIM'])
    spectraA = np.zeros(N)
    spectraB = np.zeros(N)

    for i in range(N):
        spectraA[i] *= 0.
        spectraB[i] *= 0.
        for j in range(1, 4):
            # spectraA[i] += energies_A[j] * np.abs(mu[j, 0]) ** 2 * gamma[j, 0] / (
            #         (energies_A[j] - omega[i]) ** 2 + gamma[j, 0] ** 2)
            # spectraB[i] += energies_B[j] * np.abs(mu[j, 0]) ** 2 * gamma[j, 0] / (
            #         (energies_B[j] - omega[i]) ** 2 + gamma[j, 0] ** 2)
            spectraA[i] += (energies_A[j] / (energies_A[j]**2 - (omega[i] - 1j*gamma[j, 0])**2)).imag * omega[i]
            spectraB[i] += (energies_B[j] / (energies_B[j]**2 - (omega[i] - 1j*gamma[j, 0])**2)).imag * omega[i]

    system = RamanControl(**params)
    system.call_raman_optimization()
    fig = plt.figure()
    plt.plot(system.time_spectra, system.field_spectra)

    fig1 = plt.figure()
    plt.plot(1239.84193 / system.freq_spectra * energy_factor, system.spectra_A.real / np.abs(system.spectra_A).max(), 'r')
    plt.plot(1239.84193 / system.freq_spectra * energy_factor, system.spectra_B.real / np.abs(system.spectra_A).max(), 'b')
    # plt.plot(1239.84193 / omega * energy_factor, -spectraA / np.abs(spectraA).max(), 'r--')
    # plt.plot(1239.84193 / omega * energy_factor, -spectraB / np.abs(spectraA).max(), 'b--')
    plt.show()

