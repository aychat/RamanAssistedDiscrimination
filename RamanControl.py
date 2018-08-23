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
        self.field = np.empty(self.timeDIM, dtype=np.complex)

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
        func_params.field_out = self.field.ctypes.data_as(POINTER(c_complex))

        return

    def call_raman_optimization(self):
        """

        :return:
        """
        MolA = Molecule()
        MolB = Molecule()
        func_params = Parameters()
        self.create_molecule(MolA, MolB)
        self.create_parameters(func_params)
        RamanControlFunction(MolA, MolB, func_params)


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    np.set_printoptions(precision=4)

    energy_factor = 1. / 27.211385
    time_factor = .02418
    timeAMP = 20000.
                                                                                ##################################################################
                                                                                #                                                                #
    energies_A = np.array((0.000, 0.07439, 1.94655, 2.02094)) * energy_factor   #   Energies in a.u. for molecule A                              #
    energies_B = np.array((0.000, 0.07639, 1.96655, 2.04094)) * energy_factor   #   Energies in a.u. for molecule B                              #
                                                                                #                                                                #
    N = len(energies_A)                                                         #                                                                #
    rho_0 = np.zeros((N, N), dtype=np.complex)                                  #                                                                #
    rho_0[0, 0] = 1. + 0j                                                       #   Only ground state populated                                  #
    mu = 4.97738 * np.ones_like(rho_0)                                          #   Dipole moment in a.u. (same for all allowed transition)      #
    np.fill_diagonal(mu, 0j)                                                    #   \mu_{i, i} = 0.                                              #
                                                                                #                                                                #
    population_decay = 2.418884e-8                                              #   All population relaxation lifetimes equal 1 ns               #
    electronic_dephasing = 2.418884e-4                                          #   All electronic dephasing 50 fs                               #
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
    )

    for i in range(10):
        system = RamanControl(**params)
        system.A_R *= (i / 10.) + 1.
        system.call_raman_optimization()
        plt.plot(system.time, system.dyn_rho_A.real[0], label=str(system.A_R))
        del system
    plt.legend()
    plt.show()