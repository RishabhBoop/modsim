from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from OMPython import ModelicaSystem
import os

# Check if EXPORT_PLOT is defined, otherwise default to False
try:
    EXPORT_PLOT
except NameError:
    EXPORT_PLOT = False

script_dir = os.path.dirname(os.path.abspath(__file__))
g = 9.81
l = 1.0
d = 0.3 # Dämpfung; relevant für a5

def om_non_lin_dgl(stoptime=5.0):
    modelname = "K2_3"
    model_path = os.path.join(script_dir, modelname + ".mo")
    mod = ModelicaSystem(model_path, modelname)
    mod.setParameters(
        [
            f"g = {g}",
            f"l = {l}",
            f"d = {d}"
        ]
    )
    mod.setSimulationOptions(f"stopTime={stoptime}")
    mod.simulate()
    [t] = mod.getSolutions("time")
    [phi] = mod.getSolutions("phi")
    [phi_der] = mod.getSolutions("phi_der")

    return t, phi, phi_der

def om_lin_dgl(stoptime=5.0):
    modelname = "K2_2"
    model_path = os.path.join(script_dir, modelname + ".mo")
    mod = ModelicaSystem(model_path, modelname)
    mod.setSimulationOptions(f"stopTime={stoptime}")
    mod.simulate()
    [t] = mod.getSolutions("time")
    [phi] = mod.getSolutions("phi")
    [phi_der] = mod.getSolutions("phi_der")

    return t, phi, phi_der

def xdot_fkt(t, x, *args):
    xdot = [x[1], -g / l * x[0]]
    return xdot

def lin_fkt(t, x):
        return [x[1], -d * x[1] - (g / l) *  x[0]]

def scipy_lin_dgl(x0, t_max):
    t = np.linspace(0, t_max, 101)
    sol = solve_ivp(xdot_fkt, [0, t_max], x0, t_eval=t, method="LSODA", args=(g, l))
    t = sol.t
    x = sol.y
    return t, x

def aufgabe1():
    # define global variables
    t_max = 5
    t = np.linspace(0, t_max, 101)

    # ------------------- Variante 1 -------------------
    print("VARIANTE 1 mit scipy.integrate.solve_ivp")

    x0 = [np.pi / 4, 0]

    # solve ODE
    t, x = scipy_lin_dgl(x0, t_max)

    plt.plot(t, x[0, :], label="Winkel (rad)")
    plt.xlabel("Zeit (s)")
    plt.title("Schwingung eines einfachen Pendels")
    plt.legend()
    plt.grid()
    if EXPORT_PLOT:
        plt.savefig('k2_a1_scipy.png', format='png', bbox_inches="tight", dpi=600, transparent=True)
        print("Plot saved to k2_a1_scipy.png")
    else:
        plt.show()

    # ------------------- Variante 2 -------------------
    print("VARIANTE 2 mit Modelica")
    # Get the directory where this script is located
    t, phi, phi_der = om_lin_dgl()
    plt.plot(t, phi, label="Winkel (rad)")
    plt.xlabel("Zeit (s)")
    plt.title("Kapitel 2 - Aufgabe 1: Schwingung eines einfachen Pendels (Modelica)")
    plt.legend()
    plt.grid()
    if EXPORT_PLOT:
        plt.savefig('k2_a1_modelica.png', format='png', bbox_inches="tight", dpi=600, transparent=True)
        print("Plot saved to k2_a1_modelica.png")
    else:
        plt.show()


def aufgabe2():
    # define global variables
    t_max = 5
    t = np.linspace(0, t_max, 101)

    # ------------------- Variante 1 -------------------
    print("VARIANTE 1 mit scipy.integrate.solve_ivp")

    # define the ODE system
    def xdot_fkt(t, x, *args):
        xdot = [x[1], -g / l * x[0] - d * x[1]]
        return xdot

    x0 = [np.pi / 4, 0]

    # solve ODE
    sol = solve_ivp(xdot_fkt, [0, t_max], x0, t_eval=t, method="LSODA", args=(g, l))
    t = sol.t
    x = sol.y

    plt.plot(t, x[0, :], label="Winkel (rad)")
    plt.xlabel("Zeit (s)")
    plt.title("Kapitel 2 - Aufgabe 2: Schwingung eines einfachen Pendels")
    plt.legend()
    plt.grid()
    if EXPORT_PLOT:
        plt.savefig('k2_a2_scipy.png', format='png', bbox_inches="tight", dpi=600, transparent=True)
        print("Plot saved to k2_a2_scipy.png")
    else:
        plt.show()

    # ------------------- Variante 2 -------------------
    print("VARIANTE 2 mit Modelica")
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    modelname = "K2_2"
    model_path = os.path.join(script_dir, modelname + ".mo")
    mod = ModelicaSystem(model_path, modelname)
    mod.setSimulationOptions("stopTime=5.0")
    mod.simulate()
    [t] = mod.getSolutions("time")
    [phi] = mod.getSolutions("phi")
    plt.plot(t, phi, label="Winkel (rad)")
    plt.xlabel("Zeit (s)")
    plt.title("Kapitel 2 - Aufgabe 2: Schwingung eines einfachen Pendels (Modelica)")
    plt.legend()
    plt.grid()
    if EXPORT_PLOT:
        plt.savefig('k2_a2_modelica.png', format='png', bbox_inches="tight", dpi=600, transparent=True)
        print("Plot saved to k2_a2_modelica.png")
    else:
        plt.show()


def aufgabe3():
    t, phi, phi_der = om_non_lin_dgl()
    plt.plot(t, phi, label="Nichtlineare DGL")
    # plt.text("Nichtlineare DGL")

    # A2 for comparison
    t, phi, phi_der = om_lin_dgl()
    plt.plot(t, phi, label="linearisierte DGL")
    # plt.text("Vereinfachte DGL")

    plt.xlabel("Zeit (s)")
    plt.title("Schwingung eines einfachen Pendels (Modelica)")
    plt.legend()
    plt.grid()
    if EXPORT_PLOT:
        plt.savefig('k2_a3.png', format='png', bbox_inches="tight", dpi=600, transparent=True)
        print("Plot saved to k2_a3.png")
    else:
        plt.show()
    pass


def aufgabe4():
    # Zusatzfrage Antwort:
    print(
        "\nAntwort zur Zusatzfrage:\n"
        "Bei der Berechnung des Trägheitsmomentes J = J_Steiner + J_Eigen\n"
        "wird das Eighenträgheitsmoment der Kugel (2/5 * m * r^2) vernachlässigt,\n"
        "da der Radius der Kugel (0.05m) im Vergleich zur Pendellänge (1m) klein ist."
    )


def aufgabe5():
    print("\n--- AUFGABE 5: Energetische Analyse ---")
    m = 1.0
    t_max = 10.0  # Increased time to see damping effect clearly
    
    # 1. Get simulation data
    t, phi, phi_der = om_non_lin_dgl(stoptime=t_max)

    # 2. Calculate Energies using NumPy arrays
    # Potential Energy: m*g*h where h = l - l*cos(phi)
    E_pot = m * g * l * (1 - np.cos(phi))
    
    # Kinetic Energy: 1/2 * m * v^2 where v = l * phi_der
    E_kin = 0.5 * m * (l * phi_der)**2
    
    # Total Energy
    E_total = E_pot + E_kin
    
    # Dissipated Energy (Energy lost since start)
    # E_total[0] is the initial energy
    E_diss = E_total[0] - E_total

    # 3. Plotting
    plt.figure()
    
    plt.plot(t, E_pot, label=r"$E_{pot}$ (Lage)")
    plt.plot(t, E_kin, label=r"$E_{kin}$ (Bewegung)")
    plt.plot(t, E_diss, label=r"$E_{diss}$ (Verlust/Reibung)")
    plt.plot(t, E_total, "k--", linewidth=2, label=r"$E_{ges}$ (Gesamt)")

    plt.title("Kapitel 2 - Aufgabe 5: Energetische Analyse (Nichtlineares Pendel)")
    plt.xlabel("Zeit (s)")
    plt.ylabel("Energie (J)")
    plt.legend(loc="right")
    plt.grid(True)
    if EXPORT_PLOT:
        plt.savefig('k2_a5_nichtlinear.png', format='png', bbox_inches="tight", dpi=600, transparent=True)
        print("Plot saved to k2_a5_nichtlinear.png")

    # --------------------------
    # Linear DGL
    t_lin, phi_lin, phi_der_lin = om_lin_dgl(t_max)

    E_pot_lin = m * g * l * (1 - np.cos(phi_lin))

    E_kin_lin = 0.5 * m * (l * phi_der_lin)**2

    E_total_lin = E_pot_lin + E_kin_lin
    
    E_diss_lin = E_total_lin[0] - E_total_lin


    plt.figure()
    
    plt.plot(t_lin, E_pot_lin, label=r"$E_{pot}$ (Lage)")
    plt.plot(t_lin, E_kin_lin, label=r"$E_{kin}$ (Bewegung)")
    plt.plot(t_lin, E_diss_lin, label=r"$E_{diss}$ (Verlust/Reibung)")
    plt.plot(t_lin, E_total_lin, "k--", linewidth=2, label=r"$E_{ges}$ (Gesamt)")

    plt.title("Kapitel 2 - Aufgabe 5: Energetische Analyse (Linearisiertes Pendel)")
    plt.xlabel("Zeit (s)")
    plt.ylabel("Energie (J)")
    plt.legend(loc="right")
    plt.grid(True)
    if EXPORT_PLOT:
        plt.savefig('k2_a5_linearisiert.png', format='png', bbox_inches="tight", dpi=600, transparent=True)
        print("Plot saved to k2_a5_linearisiert.png")



    if not EXPORT_PLOT:
        plt.show()

    # Wenn Sie mit der linearisierten DGL rechnen, aber die reale Energieformel (1−cosϕ) verwenden, verletzen Sie den Energieerhaltungssatz, da die Lösung ϕ(t) nicht exakt zu diesem Potenzial gehört. 
    # Die linearisierte DGL gilt nur für kleine Winkel, bei denen (1/2​)ϕ^2≈1−cos(ϕ) ist. Bei π/4 ist dieser Fehler nicht mehr vernachlässigbar.
    print("Antwort auf Frage:")
    print("Durch die Linearisierung wird die Simulation ungenau und die dissapated Energy rutscht teilweise ins Negative, was physikalisch unmöglich ist.")

if __name__ == "__main__":
    aufgabe1()
    aufgabe2()
    aufgabe3()
    aufgabe4()
    aufgabe5()
