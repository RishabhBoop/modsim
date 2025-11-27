import numpy as np 
#import scipy 
import matplotlib.pyplot as plt
from numpy import pi as pi
from scipy.integrate import solve_ivp

u_dach=2
f=1
R=20
L=9e-3
Cap=1000e-6

t_max=2
Delta_t=1e-3
t_aequidistant=np.linspace(0, t_max, int(t_max/Delta_t)+1)

A=np.array([[-1/(R*Cap), -1/Cap], 
            [1/L,         0    ]])
B=np.array([[1/(R*Cap)],
               [0]])
C=np.array([[   1,  0],
            [   0,  1],
            [ -1/R, 0]])
D=np.array([[  0],
            [  0],
            [1/R]])

def xdot_fkt(t, x):
    x=np.reshape(x, (2,1))
# Ergänzen Sie hier die Gleichung zur Berechnung der Ableitungen:
    u = u_dach * np.sin(2*pi*f*t)
    xdot = A.dot(x) + B * u
    xdot = np.reshape(xdot, (2,))
    return xdot

u_a_0=0; i_l_0=0 # Anfangswerte
x0=np.array([u_a_0, i_l_0])

sol = solve_ivp(xdot_fkt, [0, t_max], x0, t_eval=t_aequidistant, method='LSODA')
#sol = solve_ivp(xdot_fkt, [0, t_max], x0, t_eval=t_aequidistant)
x=sol.y;
# Ergänzen Sie hier die Gleichung zur Berechnung der Ausgangsgrößen
y = C.dot(x) + D * (u_dach * np.sin(2*pi*f*t_aequidistant))


u_a=y[0]; i_c=y[2]-y[1]; i_e=y[2]



print("Ohne LSODA verwendet es zur Lösung der DGLs den Runge-Kutta-Fehlberg-Algorithmus (RK45) als Standardmethode. Dies ist nicht so gut geeignet für steife Gleichungssysteme, wie sie in diesem Fall vorliegen.")
print('-'*30)
print("Eigenfrequenz der Schaltung:", 1/(2*pi*np.sqrt(L*Cap)), "Hz")
print("Annährendes ohmsches verhalten wird durch die Parallelschaltung von L und C bei Anregung mit niedriger Frequenz erreicht. " \
"Die Induktiviität wirkt wie ein Kurzschluss und der Kondensator wie eine offene Schaltung, so dass der Widerstand dominiert.")

fig=plt.figure(1, figsize=(10,6)); fig.clf()
ax = fig.add_subplot(211)
ax.plot(t_aequidistant, u_a, color='black', label=r'$u_a$ (V)')
ax.plot(t_aequidistant, i_c, color='red',   label=r'$i_c$ (A)')
# ax.plot(t_aequidistant, 0.001*u_dach*np.sin(2*pi*53*t_aequidistant), color='darkgreen', label=r'$u_e$ (V)')
ax.grid()
ax.set_ylabel(r'$u_a,\; i_c$')
ax.legend(loc='upper right', frameon=True)

ax = fig.add_subplot(212)
ax.plot(t_aequidistant, 10*i_e, color='blue',      label=r'$10\cdot i_e$ (A)')
ax.plot(t_aequidistant, u_dach*np.sin(2*pi*f*t_aequidistant), color='darkgreen', label=r'$u_e$ (V)')
ax.set_ylabel(r'$u_e,\;10\cdot i_e$')
ax.set_xlabel('t (s)')
ax.grid()
ax.legend(loc='upper right', frameon=True)

plt.tight_layout()
plt.show()
