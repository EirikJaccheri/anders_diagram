import matplotlib.pyplot as plt
import numpy as np
from thermopack.multiparameter import multiparam


# initial state
T0 = 295.3  # (K)
P0 = 2000000  # (Pa)
Z0 = [1.]
T_C = 278.15  # (K) temperature of heat exchanger
T_TANK = 273.15 - 53  # (K) temperature inside final tank

# after pumping
P1 = 45.0e5  # (Pa)


def get_phase_path(eos, t0, p0, z0):
    # plot two phase envelope
    t, p, v = eos.get_envelope_twophase(
        5.0e3, z0, maximum_pressure=P0, minimum_temperature=273.3, calc_v=True)
    s = [eos.entropy_tv(t[i], v[i], z0) for i in range(len(t))]
    fig, ax = plt.subplots()
    ax.plot(s, t, label="Two phase envelope")
    ax.set_xlabel("Entropy (J/(K*kg)")
    ax.set_ylabel("Temperature (K)")

    # initial preasure increase to 45 bar (starting in pure gas)
    v0, = eos.specific_volume(t0, p0, z0, eos.VAPPH)
    s0, = eos.entropy_tv(t0, v0, z0)
    s1 = s0
    t1 = eos.two_phase_psflash(P1, z0, s1).T
    ax.plot([s0, s0], [t0, t1], label="isentropic compression")

    # isobar as the gas is cooled
    t_bub, xbub = eos.bubble_temperature(P1, z0)
    T_iso_p_vap = np.linspace(t1, t_bub, 100)
    T_iso_p_liq = np.linspace(t_bub, T_C, 100)
    T_iso_p = np.concatenate((T_iso_p_vap, T_iso_p_liq))
    s_iso_p = [eos.entropy(t, P1, z0, eos.VAPPH) for t in T_iso_p_vap] + [
        eos.entropy(t, P1, z0, eos.LIQPH) for t in T_iso_p_liq]
    ax.plot(s_iso_p, T_iso_p, label="isobaric cooling")

    # isentalpic cooling to get gas liquid mixture
    # p_c = eos.pressure_ts(T_C, z0)
    # h, = eos.enthalpy(T_C, P1, z0, eos.LIQPH)
    # T_START = T_iso_p[-1]
    # T_iso_h, p_iso_h, v_iso_h, s_iso_h = eos.get_isenthalp(h, z0, minimum_pressure=1e5, maximum_pressure=P1,
    #                                                        minimum_temperature=T_TANK, maximum_temperature=T_START,
    #                                                        nmax=100)
    # ax.plot(s_iso_h, T_iso_h, label="isenthalpic cooling")

    h, = eos.enthalpy(T_C, P1, z0, eos.LIQPH)

    T_iso_h, p_iso_h, v_iso_h, s_iso_h = eos.get_isenthalp(h, z0, minimum_pressure=1e5,
                                                           maximum_pressure=P1,
                                                           minimum_temperature=T_TANK,
                                                           maximum_temperature=500.0,
                                                           nmax=100)

    # getting liquid and gas fraction
    p_end = p_iso_h[0]
    flsh = eos.two_phase_phflash(p_end, z0, h)

    print(flsh)

    ax.plot(s_iso_h, T_iso_h, label="Isenthalpic throttling")

    ax.set_xlim(-100, -20)
    ax.set_ylim(T_TANK, 380)
    ax.legend()
    plt.savefig("plots/co2_phase_path.png")
    plt.show()


if __name__ == "__main__":
    GERGCO2 = multiparam("CO2", "GERG2008")
    get_phase_path(GERGCO2, T0, P0, Z0)
