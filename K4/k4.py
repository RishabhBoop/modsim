import numpy as np    
import matplotlib.pyplot as plt
import math
import time

## Aufgabe 1-4

def i_t(t, U=300, R=10, L=10e-3, i0=0, t_start=0):
    # Calculate I_inf (long term current)
    I_inf = U / R
    # Calculate time relative to the start of the step
    dt = t - t_start
    return I_inf + (i0 - I_inf) * np.exp((-R / L) * dt)


def u_t(u_o, duty, T_p):
    T_on = duty * T_p
    T_halbe = (T_p - T_on) / 2
    t_u = np.array([0, T_halbe, T_halbe + T_on, T_p])
    u = np.array([0, u_o, 0, 0])
    print("U(t) =", u)
    print("t_u =", t_u)
    return u, t_u



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

    T_1 = duty*T_P # Dauer für Ein
    T_0 = (1-duty)*T_P # Dauer für aus
    T_0_halbe = T_0/2

    N_0 = np.round((1-duty)/2*(N-1)) # Abschnitte im ausschaltbereich
    N_1 = (N - 2*N_0)-1 # Abschnitte im anschaltbereich (gesamte anzahl - 2*punkte im ausschaltbereich)

    if N_0 < 1:
        N_0 = 1

    if N_1 < 2:
        if N%2 == 0:
            N_1 = 2
        else:
            N_1 = 1
        N_1 = 1 # Abschnitte im anschaltbereich (gesamte anzahl - 2*punkte im ausschaltbereich)

    delta_t1 = T_0/(N_0*2) # im aus
    delta_t2 = (T_1)/N_1 # im an

    # vectors
    t_a_vec = np.arange(0, int(N_0)+1) * delta_t1
    t_b_vec = T_0_halbe + np.arange(0, int(N_1)+1) * delta_t2
    t_c_vec = T_0/2 + T_1 + np.arange(0, int(N_0)+1) * delta_t1

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

def plot_RL_circuit(t_merged, u_merged, i_merged, t_u_total):
    # Plot the merged vectors
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot voltage
    ax1.plot(t_u_total, u_merged, 'g-', drawstyle='steps-post', label='Voltage (V)')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Voltage (V)', color='g')
    ax1.tick_params(axis='y', labelcolor='g')
    ax1.grid(True)

    # Create a second y-axis for current
    ax2 = ax1.twinx()
    ax2.plot(t_merged, i_merged, 'r-o', label='Current (A)')
    ax2.set_ylabel('Current (A)', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    # Title
    plt.title('Lösung zu Aufgabe 1')
    fig.tight_layout()
    plt.show()

u_0 = 300
R = 10
L = 10e-3
T_p = 5e-3
N=12
u_mittel_vec=np.array([90, 150, 210])
t_merged, u_merged, i_merged = np.array([]), np.array([]), np.array([])
t_total, u_total, i_total, t_u_total = np.array([]), np.array([]), np.array([]), np.array([])
for k in range(u_mittel_vec.size):
    duty = abs(u_mittel_vec[k]/u_0)
    i_0 = i_merged[-1] if i_merged.size > 0 else 10
    t_merged, i_merged = RL_circuit(duty=duty, T_P=T_p, N=N, U_0=u_0, i_0=i_0, R=R, L=L, t_start=k*T_p)
    u_merged, t_u = u_t(u_0, duty, T_p)
    t_u += k*T_p
    
    t_u_total = np.append(t_u_total, t_u)
    t_total = np.append(t_total, t_merged)
    u_total = np.append(u_total, u_merged)
    i_total = np.append(i_total, i_merged)
    #plot_RL_circuit(t_merged, u_merged, i_merged)

plot_RL_circuit(t_total, u_total, i_total, t_u_total)