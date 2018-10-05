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
        self.field_t = np.empty(params.timeDIM, dtype=np.complex)

        self.gamma_decay = np.ascontiguousarray(self.gamma_decay)
        self.gamma_pure_dephasing = np.ascontiguousarray(self.gamma_pure_dephasing)
        self.mu = np.ascontiguousarray(self.mu)
        self.rho_0 = np.ascontiguousarray(params.rho_0)
        self.rho = np.ascontiguousarray(params.rho_0.copy())
        self.energies = np.ascontiguousarray(self.energies)
        N = len(self.energies)
        self.dyn_rho = np.ascontiguousarray(np.zeros((N+1, params.timeDIM), dtype=np.complex))

    def create_molecule(self, mol):
        mol.nDIM = len(self.energies)
        mol.energies = self.energies.ctypes.data_as(POINTER(c_double))
        mol.gamma_decay = self.gamma_decay.ctypes.data_as(POINTER(c_double))
        mol.gamma_pure_dephasing = self.gamma_pure_dephasing.ctypes.data_as(POINTER(c_double))
        mol.mu = self.mu.ctypes.data_as(POINTER(c_complex))
        mol.rho = self.rho.ctypes.data_as(POINTER(c_complex))
        mol.dyn_rho = self.dyn_rho.ctypes.data_as(POINTER(c_complex))
        mol.rho_0 = self.rho_0.ctypes.data_as(POINTER(c_complex))

    def create_parameters(self, func_params, params):
        func_params.time = self.time.ctypes.data_as(POINTER(c_double))
        func_params.rho_0 = self.rho_0.ctypes.data_as(POINTER(c_complex))

        width_exc = params.timeAMP / params.width_exc
        func_params.A_exc = params.A_exc
        func_params.width_exc = width_exc
        func_params.w_exc = params.w_exc
        func_params.nDIM = len(self.energies)
        func_params.timeDIM = len(self.time)
        func_params.field_out = self.field_t.ctypes.data_as(POINTER(c_complex))

    def call_dynamics(self, params):
        mol = Molecule()
        self.create_molecule(mol)
        func_params = Parameters()
        self.create_parameters(func_params, params)
        Dynamics(mol, func_params)
        return


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    # from itertools import *

    np.set_printoptions(precision=4)
    energy_factor = 1. / 27.211385
    time_factor = .02418884 / 1000

    energies = np.array((0.000, 1.87855)) * energy_factor
    N = len(energies)
    rho_0 = np.zeros((N, N), dtype=np.complex)
    rho_0[0, 0] = 1. + 0j

    mu = 4.97738 * np.ones_like(rho_0)
    np.fill_diagonal(mu, 0j)
    population_decay = 2.418884e-8
    dephasing_rate = 5. * 2.418884e-4
    gamma_decay = np.ones((N, N)) * population_decay
    np.fill_diagonal(gamma_decay, 0.0)
    gamma_decay = np.tril(gamma_decay)
    gamma_pure_dephasing = np.ones_like(gamma_decay) * dephasing_rate
    np.fill_diagonal(gamma_pure_dephasing, 0.0)

    gamma_pure_dephasing *= 0.0
    gamma_decay *= 0.0

    params = ADict(
        energy_factor=energy_factor,
        time_factor=time_factor,
        timeDIM=12000,
        timeAMP=2500,

        A_exc=0.0025,
        width_exc=5.,
        w_exc=energies[1],
        rho_0=rho_0,
    )

    FourLevels = dict(
        energies=energies,
        gamma_decay=gamma_decay,
        gamma_pure_dephasing=gamma_pure_dephasing,
        mu=mu,
    )

    def render_ticks(axes):
        axes.get_xaxis().set_tick_params(which='both', direction='in', width=1, labelrotation=0, labelsize='large')
        axes.get_yaxis().set_tick_params(which='both', direction='in', width=1, labelcolor='r', labelsize='large')
        axes.get_xaxis().set_ticks_position('both')
        axes.get_yaxis().set_ticks_position('both')
        axes.grid()

    molecules = RamanControl(params, **FourLevels)
    molecules.call_dynamics(params)

    fig, axes = plt.subplots(nrows=1, ncols=1)
    axes.plot(molecules.time*time_factor, molecules.field_t.real, 'r')

    fig1, axes = plt.subplots(nrows=2, ncols=1, sharex=True)

    axes[0].plot(molecules.time*time_factor, molecules.dyn_rho[0, :], label='11', linewidth=2.)
    # axes[0].plot(molecules.time*time_factor, molecules.dyn_rho[1, :], label='22', linewidth=2.)
    # axes[0].plot(molecules.time*time_factor, molecules.dyn_rho[2, :], 'k', label='Tr[$\\rho^2$]', linewidth=2.)
    axes[0].plot(molecules.time*time_factor, np.cos(np.abs(mu[1, 0]*params.A_exc)*molecules.time)**2, 'k--', label='11 - 22', linewidth=2.)
    # axes[0].plot(molecules.time*time_factor, np.sin(np.abs(mu[1, 0]*params.A_exc)*molecules.time)**2, 'm--', label='11 - 22', linewidth=2.)

    axes[0].legend(loc=4)

    axes[1].plot(molecules.time*time_factor, molecules.dyn_rho[0, :] - molecules.dyn_rho[1, :], 'k', label='11 - 22', linewidth=2.)

    # T1 = 1. / population_decay
    # Tc = 1. / dephasing_rate
    # T2 = 1. / ((0.5 / T1) + (1. / Tc))
    # ref = (1 + (params.w_exc - energies[1])**2 * T2**2)/(1 + (params.w_exc - energies[1])**2 * T2**2 + 4*mu[1, 0]**2 * params.A_exc**2 * T1 * T2)
    # print(ref)
    axes[1].legend(loc=4)

    render_ticks(axes[0])
    render_ticks(axes[1])

    print(molecules.rho.real, "\n")

    plt.show()

