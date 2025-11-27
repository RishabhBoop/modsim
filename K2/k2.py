from scipy.integrate import solve_ivp
import numpy as np    
import matplotlib.pyplot as plt
from OMPython import ModelicaSystem
import time
import os


def aufgabe1():
    # define global variables
    t_max=5;
    t=np.linspace(0, t_max,101)
    g = 9.81
    l = 1.0

    # ------------------- Variante 1 -------------------
    print("VARIANTE 1 mit scipy.integrate.solve_ivp")
    # define the ODE system
    def xdot_fkt(t, x, *args):
        xdot= [x[1], -g/l*x[0]]
        return xdot
    
    x0=[np.pi/4, 0]

    # solve ODE
    sol = solve_ivp(xdot_fkt, [0, t_max], x0, t_eval=t, method='LSODA', args=(g, l))
    t=sol.t
    x=sol.y

    plt.plot(t, x[0,:], label='Winkel (rad)')
    plt.xlabel('Zeit (s)')
    plt.title('Schwingung eines einfachen Pendels')
    plt.legend()
    plt.grid()
    plt.show()


    # ------------------- Variante 2 -------------------
    print("VARIANTE 2 mit Modelica")
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    modelname = "K2_1"
    model_path = os.path.join(script_dir, modelname + ".mo")
    mod = ModelicaSystem(model_path, modelname)
    mod.setSimulationOptions('stopTime=5.0')
    mod.simulate()
    [t] = mod.getSolutions('time')
    [phi] = mod.getSolutions('phi')
    plt.plot(t, phi, label='Winkel (rad)')
    plt.xlabel('Zeit (s)')
    plt.title('Schwingung eines einfachen Pendels (Modelica)')
    plt.legend()
    plt.grid()
    plt.show()


def aufgabe2():
    # define global variables
    t_max=5;
    t=np.linspace(0, t_max,101)
    g = 9.81
    l = 1.0
    d = 0.3

    # ------------------- Variante 1 -------------------
    print("VARIANTE 1 mit scipy.integrate.solve_ivp")
    # define the ODE system
    def xdot_fkt(t, x, *args):
        xdot= [x[1], -g/l*x[0] - d*x[1]]
        return xdot
    
    x0=[np.pi/4, 0]

    # solve ODE
    sol = solve_ivp(xdot_fkt, [0, t_max], x0, t_eval=t, method='LSODA', args=(g, l))
    t=sol.t
    x=sol.y

    plt.plot(t, x[0,:], label='Winkel (rad)')
    plt.xlabel('Zeit (s)')
    plt.title('Schwingung eines einfachen Pendels')
    plt.legend()
    plt.grid()
    plt.show()


    # ------------------- Variante 2 -------------------
    print("VARIANTE 2 mit Modelica")
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    modelname = "K2_2"
    model_path = os.path.join(script_dir, modelname + ".mo")
    mod = ModelicaSystem(model_path, modelname)
    mod.setSimulationOptions('stopTime=5.0')
    mod.simulate()
    [t] = mod.getSolutions('time')
    [phi] = mod.getSolutions('phi')
    plt.plot(t, phi, label='Winkel (rad)')
    plt.xlabel('Zeit (s)')
    plt.title('Schwingung eines einfachen Pendels (Modelica)')
    plt.legend()
    plt.grid()
    plt.show()



def aufgabe3():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    print(script_dir)
    modelname = "K2_3"
    model_path = os.path.join(script_dir, modelname + ".mo")
    mod = ModelicaSystem(model_path, modelname)
    mod.setSimulationOptions('stopTime=5.0')
    mod.simulate()
    [t] = mod.getSolutions('time')
    [phi] = mod.getSolutions('phi')
    plt.plot(t, phi, label='Nichtlineare DGL')
    # plt.text("Nichtlineare DGL")


    # A2 for comparison
    modelname = "K2_2"
    model_path = os.path.join(script_dir, modelname + ".mo")
    mod = ModelicaSystem(model_path, modelname)
    mod.setSimulationOptions('stopTime=5.0')
    mod.simulate()
    [t] = mod.getSolutions('time')
    [phi] = mod.getSolutions('phi')
    plt.plot(t, phi, label='linearisierte DGL')
    # plt.text("Vereinfachte DGL")

    plt.xlabel('Zeit (s)')
    plt.title('Schwingung eines einfachen Pendels (Modelica)')
    plt.legend()
    plt.grid()
    plt.show()
    pass


aufgabe_choice = input("Welche Aufgabe möchten Sie ausführen? ")
if aufgabe_choice == "1":
    aufgabe1()
elif aufgabe_choice == "2":
    aufgabe2()
elif aufgabe_choice == "3":
    aufgabe3()
# elif aufgabe_choice == "4":
#     aufgabe4()
else:
    print("Ungültige Eingabe. Bitte gültige Zahl eingeben.")