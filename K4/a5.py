import numpy as np
import matplotlib.pyplot as plt

## Aufgabe 5


def i_t(t, U=300, R=10, L=10e-3, i0=0, t_start=0):
    """Calculate current in RL circuit with step voltage input."""
    # Calculate I_inf (long term current)
    I_inf = U / R
    # Calculate time relative to the start of the step
    dt = t - t_start
    return I_inf + (i0 - I_inf) * np.exp((-R / L) * dt)


def u_verlauf(u_o, duty, T_p):
    """Generate voltage profile for PWM cycle."""
    T_on = duty * T_p
    T_halbe = (T_p - T_on) / 2
    t_u = np.array([0, T_halbe, T_halbe + T_on, T_p])
    u = np.array([0, u_o, 0, 0])
    return u, t_u


def u_sin(t, U_hat=200, f=50):
    """Sinusoidal voltage u(t) = U_hat * cos(2*pi*f*t)"""
    return U_hat * np.cos(2 * np.pi * f * t)


def i_continuous(t, U_hat, R, L, omega):
    """
    Analytical steady-state current for continuous sinusoidal voltage.
    u(t) = U_hat * cos(omega*t)
    i(t) = I_hat * cos(omega*t - phi)
    where I_hat = U_hat / |Z| and phi = arctan(omega*L/R)
    """
    Z = np.sqrt(R**2 + (omega * L) ** 2)  # Impedance magnitude
    phi = np.arctan(omega * L / R)  # Phase angle
    I_hat = U_hat / Z  # Current amplitude
    return I_hat * np.cos(omega * t - phi)


def calculate_i0_no_transient(U_hat, R, L, omega):
    """
    Calculate the initial current i(0) for no transient.
    For u(t) = U_hat * cos(omega*t), steady-state: i(t) = I_hat * cos(omega*t - phi)
    At t=0: i(0) = I_hat * cos(-phi) = I_hat * cos(phi)
    """
    Z = np.sqrt(R**2 + (omega * L) ** 2)
    phi = np.arctan(omega * L / R)
    I_hat = U_hat / Z
    i_0 = I_hat * np.cos(-phi)  # = I_hat * cos(phi)
    return i_0


def RL_circuit(duty, T_P, N, U_0, i_0, R, L, t_start=0):
    """
    Simulate an RL circuit with PWM switching.

    :param duty: Duty cycle (0 to 1)
    :param T_p: Period duration in seconds
    :param N: Number of discretization points
    :param u_0: Supply voltage in volts
    :param i_0: Initial current in amperes
    :param R: Resistance in ohms
    :param L: Inductance in Henry
    :param t_start: Start time in seconds
    :return: Tuple of (time array, voltage array, current array)
    """
    i0_start = i_0

    T_1 = duty * T_P  # Dauer für Ein
    T_0 = (1 - duty) * T_P  # Dauer für aus
    T_0_halbe = T_0 / 2

    N_0 = np.round((1 - duty) / 2 * (N - 1))  # Abschnitte im ausschaltbereich
    N_1 = (
        N - 2 * N_0
    ) - 1  # Abschnitte im anschaltbereich (gesamte anzahl - 2*punkte im ausschaltbereich)

    if N_0 < 1:
        N_0 = 1

    if N_1 < 2:
        if N % 2 == 0:
            N_1 = 2
        else:
            N_1 = 1
        N_1 = 1  # Abschnitte im anschaltbereich (gesamte anzahl - 2*punkte im ausschaltbereich)

    delta_t1 = T_0 / (N_0 * 2)  # im aus
    delta_t2 = (T_1) / N_1  # im an

    # vectors
    t_a_vec = np.arange(0, int(N_0) + 1) * delta_t1
    t_b_vec = T_0_halbe + np.arange(0, int(N_1) + 1) * delta_t2
    t_c_vec = T_0 / 2 + T_1 + np.arange(0, int(N_0) + 1) * delta_t1

    # Calculate currents using numpy vectorization
    i_a_vec = i_t(t_a_vec, U=0, R=R, L=L, i0=i0_start, t_start=0)
    i_b_vec = i_t(t_b_vec, U=U_0, R=R, L=L, i0=i_a_vec[-1], t_start=t_b_vec[0])
    i_c_vec = i_t(t_c_vec, U=0, R=R, L=L, i0=i_b_vec[-1], t_start=t_c_vec[0])

    # Merge the vectors
    # We take all points from A except the last (overlap with B start)
    # All points from B except the last (overlap with C start)
    # All points from C
    t_merged = np.concatenate([t_a_vec[:-1], t_b_vec[:-1], t_c_vec])
    i_merged = np.concatenate([i_a_vec[:-1], i_b_vec[:-1], i_c_vec])

    # move each t vector by t_start
    t_merged += t_start

    return t_merged, i_merged


def plot_RL_circuit(t_pwm, i_pwm, t_u_total, u_total, t_cont, i_cont):
    """Plot PWM and continuous results for comparison."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # Plot Voltage
    ax1.plot(
        t_u_total * 1000,
        u_total,
        "g-",
        drawstyle="steps-post",
        label="PWM Spannung",
        linewidth=1,
    )
    ax1.plot(
        t_cont * 1000,
        200 * np.cos(2 * np.pi * 50 * t_cont),
        "b--",
        label="Kontinuierliche Spannung",
        linewidth=1.5,
    )
    ax1.set_ylabel("Spannung u / V")
    ax1.tick_params(axis="y", labelcolor="g")
    ax1.grid(True)
    ax1.set_title("Aufgabe 5: Vierquadrantensteller mit sinusförmiger Spannung")

    # Plot Current
    ax2 = ax1.twinx()
    ax2.plot(t_pwm * 1000, i_pwm, "r-", label="Kontinuierlicher Strom", linewidth=1)
    # ax2.set_xlabel('Zeit t / ms')
    ax2.set_ylabel("Strom i / A")
    ax2.tick_params(axis="y", labelcolor="r")
    ax2.grid(True)

    ax3 = ax2.twinx()
    ax3.plot(
        t_cont * 1000,
        i_cont,
        "m--",
        label="Kontinuierlicher Strom (analytisch)",
        linewidth=1,
    )
    # ax3.tick_params(axis='y', labelcolor='m')
    ax2.grid(True)

    fig.tight_layout()
    plt.show()


# ============ Parameters for Aufgabe 5 ============
U_0 = 300  # DC supply voltage in V for PWM
U_hat = 200  # Amplitude of sinusoidal voltage in V
f = 50  # Frequency in Hz
omega = 2 * np.pi * f
R = 1  # Resistance in Ohm
L = 10e-3  # Inductance in H
T_p = 1e-3  # PWM period in s
N = 6  # Discretization points per PWM cycle

# Simulate over one full period of the 50 Hz sine (20 ms)
T_sin = 1 / f  # = 20 ms
num_cycles = int(T_sin / T_p)  # Number of PWM cycles

# Calculate analytical starting value for no transient
i_0_analytical = calculate_i0_no_transient(U_hat, R, L, omega)
print(f"Analytischer Startwert für i(0) ohne Einschwingvorgang: {i_0_analytical:.4f} A")
print(f"Impedanz |Z| = {np.sqrt(R**2 + (omega*L)**2):.4f} Ohm")
print(f"Phasenwinkel phi = {np.degrees(np.arctan(omega*L/R)):.2f}°")

# Arrays to collect results
t_total = np.array([])
i_total = np.array([])
t_u_total = np.array([])
u_total = np.array([])

i_0_current = i_0_analytical  # Start with analytical value for no transient

for k in range(num_cycles):
    # Time at the middle of the current PWM cycle
    t_mid = k * T_p + T_p / 2

    # Mean voltage for this cycle (sampled at middle of cycle)
    u_mittel = u_sin(t_mid, U_hat, f)

    # Duty cycle is always positive
    duty = abs(u_mittel / U_0)

    # Determine the DC voltage sign based on the mean voltage sign
    if u_mittel >= 0:
        u_const = U_0
    else:
        u_const = -U_0

    # Simulate RL circuit for this PWM cycle
    t_merged, i_merged = RL_circuit(
        duty=duty, T_P=T_p, N=N, U_0=u_const, i_0=i_0_current, R=R, L=L, t_start=k * T_p
    )

    # Get voltage profile for this cycle
    u_merged, t_u = u_verlauf(u_const, duty, T_p)
    t_u += k * T_p

    # Append to total arrays (avoiding duplicate points)
    if k == 0:
        t_total = t_merged
        i_total = i_merged
        t_u_total = t_u
        u_total = u_merged
    else:
        t_total = np.append(t_total[:-1], t_merged)
        i_total = np.append(i_total[:-1], i_merged)
        t_u_total = np.append(t_u_total[:-1], t_u)
        u_total = np.append(u_total[:-1], u_merged)

    # Update initial current for next cycle
    i_0_current = i_merged[-1]

# Calculate continuous (analytical) current for comparison
t_continuous = np.linspace(0, T_sin, 1000)
i_continuous_result = i_continuous(t_continuous, U_hat, R, L, omega)

# Plot the results
plot_RL_circuit(t_total, i_total, t_u_total, u_total, t_continuous, i_continuous_result)
