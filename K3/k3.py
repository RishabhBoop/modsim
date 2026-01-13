import numpy as np

# import scipy
import matplotlib.pyplot as plt
from numpy import pi as pi
from scipy.integrate import solve_ivp
from scipy import signal
from OMPython import ModelicaSystem
import os
import math

u_dach = 2
f = 1
R = 20
L = 9e-3
Cap = 1000e-6
C = Cap


def aufgabe1():
    print("-------------- Aufgabe 1 --------------")
    t_max = 2
    Delta_t = 1e-3
    t_aequidistant = np.linspace(0, t_max, int(t_max / Delta_t) + 1)

    A = np.array([[-1 / (R * Cap), -1 / Cap], [1 / L, 0]])
    B = np.array([[1 / (R * Cap)], [0]])
    C = np.array([[1, 0], [0, 1], [-1 / R, 0]])
    D = np.array([[0], [0], [1 / R]])

    def xdot_fkt(t, x):
        x = np.reshape(x, (2, 1))
        # Ergänzen Sie hier die Gleichung zur Berechnung der Ableitungen:
        u = u_dach * np.sin(2 * pi * f * t)
        xdot = A.dot(x) + B * u
        xdot = np.reshape(xdot, (2,))
        return xdot

    u_a_0 = 0
    i_l_0 = 0  # Anfangswerte
    x0 = np.array([u_a_0, i_l_0])

    sol = solve_ivp(xdot_fkt, [0, t_max], x0, t_eval=t_aequidistant, method="LSODA")
    # sol = solve_ivp(xdot_fkt, [0, t_max], x0, t_eval=t_aequidistant)
    x = sol.y
    # Ergänzen Sie hier die Gleichung zur Berechnung der Ausgangsgrößen
    y = C.dot(x) + D * (u_dach * np.sin(2 * pi * f * t_aequidistant))

    u_a = y[0]
    i_c = y[2] - y[1]
    i_e = y[2]

    print(
        "Ohne LSODA verwendet es zur Lösung der DGLs den Runge-Kutta-Fehlberg-Algorithmus (RK45) als Standardmethode. Dies ist nicht so gut geeignet für steife Gleichungssysteme, wie sie in diesem Fall vorliegen."
    )
    print("-" * 30)
    print("Eigenfrequenz der Schaltung:", 1 / (2 * pi * np.sqrt(L * Cap)), "Hz")
    print(
        "Annährendes ohmsches verhalten wird durch die Parallelschaltung von L und C bei Anregung mit niedriger Frequenz erreicht. "
        "Die Induktiviität wirkt wie ein Kurzschluss und der Kondensator wie eine offene Schaltung, so dass der Widerstand dominiert."
    )

    fig = plt.figure(1, figsize=(10, 6))
    fig.clf()
    ax = fig.add_subplot(211)
    ax.plot(t_aequidistant, u_a, color="black", label=r"$u_a$ (V)")
    ax.plot(t_aequidistant, i_c, color="red", label=r"$i_c$ (A)")
    # ax.plot(t_aequidistant, 0.001*u_dach*np.sin(2*pi*53*t_aequidistant), color='darkgreen', label=r'$u_e$ (V)')
    ax.grid()
    ax.set_ylabel(r"$u_a,\; i_c$")
    ax.legend(loc="upper right", frameon=True)

    ax = fig.add_subplot(212)
    ax.plot(t_aequidistant, 10 * i_e, color="blue", label=r"$10\cdot i_e$ (A)")
    ax.plot(
        t_aequidistant,
        u_dach * np.sin(2 * pi * f * t_aequidistant),
        color="darkgreen",
        label=r"$u_e$ (V)",
    )
    ax.set_ylabel(r"$u_e,\;10\cdot i_e$")
    ax.set_xlabel("t (s)")
    ax.grid()
    ax.legend(loc="upper right", frameon=True)

    plt.tight_layout()
    plt.show()


def aufgabe2():
    print("-------------- Aufgabe 2 --------------")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    modelname = "A2"
    model_path = os.path.join(script_dir, modelname + ".mo")
    mod = ModelicaSystem(model_path, modelname)
    mod.setSimulationOptions("stopTime=5.0")
    mod.simulate()

    # Get solution arrays
    [t] = mod.getSolutions("time")
    [u_a] = mod.getSolutions("u_a")
    [i_L] = mod.getSolutions("i_L")

    # plot u_a and i_L over time and show legend/grid
    plt.plot(t, u_a, label="u_a (V)")
    plt.plot(t, i_L, label="i_L (A)")
    plt.xlabel("Zeit (s)")
    plt.ylabel("Wert")
    plt.title("RLC-Schaltkreis: Spannung und Strom (Modelica)")
    plt.legend()
    plt.grid()
    plt.show()


def aufgabe3():
    print("-------------- Aufgabe 3 --------------")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    modelname = "A3"
    model_path = os.path.join(script_dir, modelname + ".mo")
    mod = ModelicaSystem(model_path, modelname)
    mod.setParameters([f"R_top = {R}", f"L_top = {L}", f"C_top = {C}"])
    mod.setSimulationOptions("stopTime=5.0")
    mod.simulate()

    [t] = mod.getSolutions("time")
    [u_a] = mod.getSolutions("u_a")
    [i_e] = mod.getSolutions("i_e")
    [i_c] = mod.getSolutions("i_c")
    [i_L] = mod.getSolutions("i_L")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    # Oben: u_a und i_c
    ax1.plot(t, u_a, "k", label=r"$u_a$ (V)")
    ax1.plot(t, i_c, "r", label=r"$i_c$ (A)")
    ax1.set_ylabel(r"$u_a$, $i_c$")
    ax1.grid(True)
    ax1.legend()

    # Unten: 10*i_e und u_e
    u_e = 2 * np.sin(2 * np.pi * 1 * t)  # rekonstruiert
    ax2.plot(t, 10 * i_e, "b", label=r"$10 \cdot i_e$ (A)")
    ax2.plot(t, u_e, "darkgreen", label=r"$u_e$ (V)")
    ax2.set_xlabel("Zeit $t$ (s)")
    ax2.set_ylabel(r"$u_e$, $10 \cdot i_e$")
    ax2.grid(True)
    ax2.legend()

    plt.tight_layout()
    plt.show()


def aufgabe4():
    pass


def aufgabe5():

    print("-------------- Aufgabe 5 --------------")

    num = [L, 0]
    denum = [R * L * C, L, R]

    sys = signal.TransferFunction(num, denum)

    w, mag, phase = signal.bode(sys)

    p = sys.poles
    print("Pol 1:", p[0])
    print("Pol 2:", p[1])
    print("Grenzfrequenz (|pol|)", abs(p[0]))
    print("Grenzfrequenz (1/sqrt(LC)):", 1/np.sqrt(L*C), "Hz")

    fig, axis = plt.subplots(2, 1)

    axis[0].semilogx(w, mag)
    axis[1].semilogx(w, phase)

    #plt.title("Bode Diagramm (Bandpass)")
    axis[0].set_title('Bode Diagramm')
    axis[0].set_ylabel('Amplitude [dB]')
    axis[0].set_xlabel("Frequenz [rad/s]")
    axis[0].grid(True, which="both")
    axis[1].grid(True, which='both')
    # axis[1].set_title('Phase')
    axis[1].set_ylabel('Phase [°]')
    axis[1].set_xlabel("Frequenz [rad/s]")

    plt.show()


aufgabe1()
aufgabe2()
aufgabe3()
aufgabe4()
aufgabe5()
