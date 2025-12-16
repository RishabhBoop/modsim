import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import solve_ivp
from OMPython import ModelicaSystem
import os


def x_dot_a1(t, x):
    """
    Define the differential equation dx/dt = f(t, x).
    """
    return (-2 * t * x**2) / (t**2 + 1)


def exact_x_a1(t):
    """
    Provide the exact solution for comparison.
    """
    return 1 / (np.log(t**2 + 1) + 0.5)


def x_dot_a2(t, x):
    return x + t**3


def exact_x_a2(t):
    return -(t**3) - 3 * t**2 - 6 * t - 6 + (19 / 8) * math.exp(1.5 + t)

x_dot = x_dot_a1
exact_x = exact_x_a1


# --- Euler Method Implementation ---
def euler(x_k, t_k, h):
    """
    Implement the Euler method for numerically solving DGLs.
    x_k+1 = x_k + h * f(t_k, x_k)
    """
    return x_k + h * x_dot(t_k, x_k)


def solve_euler(x0, t0, h, steps):
    """
    Solve the DGL using the Euler method over a specified number of steps.
    """
    t_values = [t0 + i * h for i in range(steps + 1)]
    x_values = [x0]

    for k in range(steps):
        x_next = euler(x_values[-1], t_values[k], h)
        x_values.append(x_next)

    return t_values, x_values


# --- Heun Method Implementation ---
def heun(x_k, t_k, h):
    """
    Implement the Heun method for numerically solving DGLs.
    x_k+1 = x_k + (h/2) * [f(t_k, x_k) + f(t_k + h, x_k + h * f(t_k, x_k))]
    """
    pre1 = x_k + h * x_dot(t_k, x_k)
    return x_k + (h / 2) * (x_dot(t_k, x_k) + x_dot(t_k + h, pre1))


def solve_heun(x0, t0, h, steps):
    """
    Solve DGL using Heun Method
    """
    t_values = [t0 + i * h for i in range(steps + 1)]
    x_values = [x0]

    for k in range(steps):
        x_next = heun(x_values[-1], t_values[k], h)
        x_values.append(x_next)

    return t_values, x_values


# --- Runge-Kutta Method Implementation ---
def runge_kutta(x_k, t_k, h):
    """
    Implement the classical Runge-Kutta method for numerically solving DGLs.
    """
    k1 = x_dot(t_k, x_k)
    c = h / 2
    k2 = x_dot(t_k + c, x_k + c * k1)
    k3 = x_dot(t_k + c, x_k + c * k2)
    k4 = x_dot(t_k + h, x_k + h * k3)
    return x_k + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def solve_runge_kutta(x0, t0, h, steps):
    """
    Solve DGL using Runge-Kutta Method
    """
    t_values = [t0 + i * h for i in range(steps + 1)]
    x_values = [x0]

    for k in range(steps):
        x_next = runge_kutta(x_values[-1], t_values[k], h)
        x_values.append(x_next)

    return t_values, x_values


# --- Globaler Fehler ---
def global_error(numerical_values, exact_values):
    """
    Calculate the global error between numerical and exact solutions.
    """
    error = [num - exa for num, exa in zip(numerical_values, exact_values)]
    return error


def aufgabe1():
    # plot euler
    x0 = 2.0
    t0 = 0.0
    h = 0.02
    steps = 100

    t_values_euler, x_values_euler = solve_euler(x0, t0, h, steps)
    t_values_heun, x_values_heun = solve_heun(x0, t0, h, steps)
    t_values_rk, x_values_rk = solve_runge_kutta(x0, t0, h, steps)

    # global_errors_euler = global_error(x_values_euler, [exact_x(t) for t in t_values_euler])
    # global_errors_heun = global_error(x_values_heun, [exact_x(t) for t in t_values_heun])
    # global_errors_rk = global_error(x_values_rk, [exact_x(t) for t in t_values_rk])

    # plot approximation
    plt.plot(t_values_euler, x_values_euler, label="Euler Approximation")
    plt.plot(t_values_euler, x_values_heun, label="Heun Approximation")
    plt.plot(t_values_rk, x_values_rk, label="Runge-Kutta Approximation")
    plt.plot(
        t_values_euler, [exact_x(t) for t in t_values_euler], label="Exact Solution"
    )
    plt.xlabel("t")
    plt.ylabel("x")
    plt.legend()

    # # plot global error
    # plt.figure()
    # plt.plot(t_values_euler, global_errors_euler, label="Euler Global Error")
    # plt.ylabel(ylabel="Global Error Euler")
    # plt.xlabel(xlabel="t")
    # plt.legend()

    # plt.figure()
    # plt.plot(t_values_heun, global_errors_heun, label="Heun Global Error")
    # plt.ylabel(ylabel="Global Error Heun")
    # plt.xlabel("t")
    # plt.legend()

    # plt.figure()
    # plt.plot(t_values_rk, global_errors_rk, label="Runge-Kutta Global Error")
    # plt.ylabel(ylabel="Global Error Runge-Kutta")
    # plt.xlabel("t")
    # plt.legend()

    plt.show()

x_dot = x_dot_a2
exact_x = exact_x_a2

def aufgabe2():
    x0 = 2
    t0 = -1.5
    t_end = 2.0
    amount_steps = 5

    for n in range(amount_steps+1):
        steps = 2**n
        h = (t_end - t0) / steps
        t_val, x_val = solve_euler(x0, t0, h, steps)
        plt.plot(t_val, x_val, label=f"Steps: {steps}")

    plt.plot(
        np.linspace(t0, t_end, 100),
        [exact_x(t) for t in np.linspace(t0, t_end, 100)],
        "k--",
        label="Exact Solution",
    )
    plt.xlabel("t")
    plt.ylabel("x")
    plt.legend()
    plt.show()


def aufgabe3():
    """
    Berechnen Sie  x(t=2) für die Differentialgleichung aus Aufgabenteil 2 und 
    stellen Sie den globalen Fehler in Abhängigkeit der Zeitschrittweite für das 
    Euler- und das Runge-Kutta-Verfahren dar.
    Wie kann man aus dieser Abbildung die Ordnung der Verfahren ablesen?
    """
    x0 = 2
    t0 = -1.5
    t_end = 2.0
    amount_steps = 12

    exact_at_2 = exact_x(2)

    globerror_euler_x_at_2_vals = []
    globerror_rk_x_at_2_vals = []

    for n in range(amount_steps+1):
        steps = 2**n
        h = (t_end - t0) / steps
        t_val, x_val = solve_euler(x0, t0, h, steps)
        t_val, x_val_rk = solve_runge_kutta(x0, t0, h, steps)
        x_at_2 = x_val[-1]
        rk_x_at_2 = x_val_rk[-1]
        glob_error = global_error([x_at_2], [exact_at_2])[0]
        globerror_euler_x_at_2_vals.append(abs(glob_error))
        glob_error = global_error([rk_x_at_2], [exact_at_2])[0]
        globerror_rk_x_at_2_vals.append(abs(glob_error))
    
    # print("Global Error Euler at x(2):", globerror_euler_x_at_2_vals)
    # print("Global Error Runge-Kutta at x(2):", globerror_rk_x_at_2_vals)

    #plot errors vs h in log scale
    h_vals = [(t_end - t0) / (2**n) for n in range(amount_steps+1)]
    # print("h values:", h_vals)
    plt.loglog(h_vals, globerror_euler_x_at_2_vals, label="Euler Global Error at x(2)")
    plt.loglog(h_vals, globerror_rk_x_at_2_vals, label="Runge-Kutta Global Error at x(2)")
    plt.xlabel("Step size h")
    plt.ylabel("Global Error at x(2)")
    plt.legend()
    plt.show()

# --------- Selbsterarbeitung A4 ---------
t_max = 3
v_0 = 5
h_0 = 10
g = 9.81
alpha = 0.001


def xdot_fkt(t, x, *args):
    theta, omega = x
    xdot = [x[1], -g - alpha * x[1] ** 3]
    return xdot

def hit_ground(t, y, *args):
    return y[0]


hit_ground.terminal = True
hit_ground.direction = -1

def aufgabe4():
    h_0_neu = h_0
    v_0_neu = v_0
    t_0_neu = 0
    h = np.array([])
    v = np.array([])
    t_ges = np.array([])
    for k in range(11):
        # hier ist es ungeschickt, t_eval vorzugeben, weil dann der Zero Crossing Punkt nicht drin ist
        # besser: max_step
        sol = solve_ivp(
            xdot_fkt,
            [t_0_neu, 100],
            [h_0_neu, v_0_neu],
            max_step=0.1,
            events=hit_ground,
            args=(g, alpha),
        )
        t = sol.t
        x = sol.y
        t_0_neu = t[-1]
        h_0_neu = 0
        v_0_neu = -x[1, -1]
        t_ges = np.append(t_ges, t)
        h = np.append(h, x[0])
        v = np.append(v, x[1])

    # -------------- mycode --------------
    model_path = os.path.join(os.path.dirname(__file__), "bounce.mo")
    model = ModelicaSystem(model_path, "bounce")
    model.setSimulationOptions("stopTime=25")
    model.simulate()
    t_vec_mo = model.getSolutions("time")[0]
    h_mo = model.getSolutions("h")[0]
    v_mo = model.getSolutions("v")[0]
    # ------------------------------------

    fig = plt.figure(1, figsize=(10, 6))
    fig.clf()
    ax = fig.add_subplot(211)
    ax.plot(t_ges, h, "b", label="Berechnung mit solve_ivp")
    ax.plot(t_vec_mo, h_mo, "k+", label="Berechnung mit Modelica")
    ax.set_ylabel("Höhe")
    ax.legend(loc="best")
    ax.grid()
    ax = fig.add_subplot(212)
    ax.plot(t_ges, v, "g", label="Berechnung mit solve_ivp")
    ax.plot(t_vec_mo, v_mo, "k+", label="Berechnung mit Modelica")
    ax.set_ylabel("Geschwindigkeit")
    ax.legend(loc="best")
    ax.set_xlabel("t")
    ax.grid()
    plt.show()

def main():
    aufgabe4()


if __name__ == "__main__":
    main()
