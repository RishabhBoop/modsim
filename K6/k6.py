import numpy as np
from OMPython import ModelicaSystem
import matplotlib.pyplot as plt
import os


def aufgabe1(
    filename="A1.mo",
    modelname="ModSimBib.RLmitFreilauf",
    outputTitle="RL Circuit with Freewheeling Diode",
):
    """
    This Aufgabe is implementing an ohm-inductor circuit simulation using our own Modellica Block.
    """
    # fetch the modelica file
    current_path = os.path.realpath(__file__).strip(__file__.split("/")[-1])
    package_path = current_path + "ModSimBib/package.mo"
    model_path = current_path + filename
    mod = ModelicaSystem(model_path, modelname, [package_path])

    # set simulation options and simulate
    mod.setSimulationOptions(simOptions="stopTime=2.0")
    mod.simulate()
    [t] = mod.getSolutions("time")
    [i_L] = mod.getSolutions("inductor.i")
    [u_L] = mod.getSolutions("inductor.v")
    [i_q] = mod.getSolutions("switch.i")

    # plot u_L and i over time
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Top plot: Currents
    ax1.plot(t, i_L, label="i_L (Inductor current)", color="red")
    ax1.plot(t, i_q, label="i_q (Source current)", color="black")
    ax1.set_ylabel("Ströme (A)")
    ax1.set_title(outputTitle)
    ax1.legend()
    ax1.grid()

    # Bottom plot: Inductor voltage
    ax2.plot(t, u_L, label="u_L (Inductor voltage)", color="green")
    ax2.set_xlabel("Zeit (s)")
    ax2.set_ylabel("u_L (V)")
    ax2.legend()
    ax2.grid()

    plt.tight_layout()
    plt.show()


def aufgabe2():
    """
    This Aufgabe is implementing an ohm-inductor circuit simulation using Modellica's built-in blocks.
    """
    aufgabe1(
        "A2.mo",
        "ModSimBib.RLmitFreilauf_ownModels",
        "RL-Schaltkreis mit Freilaufdiode (eigene Modelle)",
    )


def aufgabe4():
    # Load OM Model
    current_path = os.path.realpath(__file__).strip(__file__.split("/")[-1])
    package_path = current_path + "ModSimBib/package.mo"
    model_path = current_path + "A4.mo"
    mod = ModelicaSystem(model_path, "myPWM", [package_path])
    
    # Set parameters
    R = 1 # Ohm
    L = 10e-3  # H
    U_DC = 300  # V
    U_dach = 300*0.7  # V
    I_start = 18.4  # A
    print("U_dach:", U_dach)
    mod.setParameters([f"R = {R}", f"L = {L}", f"U_dc = {U_DC}", f"U_dach = {U_dach}", f"i_start = {I_start}"])
    mod.setSimulationOptions("stopTime=0.04")

    # Simulate and get results
    mod.simulate()
    [t] = mod.getSolutions("time")
    [i_L] = mod.getSolutions("inductor1.i")
    [u] = mod.getSolutions("vierqst.u_out")

    # Plot results
    fig, ax1 = plt.subplots(1, 1, figsize=(12, 6))

    # Plot Voltage on primary y-axis
    ax1.plot(t * 1000, u, "g-", label="Spannung u", linewidth=1)
    ax1.set_xlabel("Zeit t / ms")
    ax1.set_ylabel("Spannung u / V", color="g")
    ax1.tick_params(axis="y", labelcolor="g")
    ax1.grid(True)
    ax1.set_title("Aufgabe 4: PWM Vierquadrantensteller")

    # Plot Current on secondary y-axis
    ax2 = ax1.twinx()
    ax2.plot(t * 1000, i_L, "r-", label="Strom i_L", linewidth=1)
    ax2.set_ylabel("Strom i / A", color="r")
    ax2.tick_params(axis="y", labelcolor="r")

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

    fig.tight_layout()
    plt.show()


# aufgabe1()
# aufgabe2()
aufgabe4()
