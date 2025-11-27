import numpy as np    
import matplotlib.pyplot as plt
from OMPython import ModelicaSystem
import os


R=20
C=1000e-6
L=9e-3

script_dir = os.path.dirname(os.path.abspath(__file__))
modelname = "A3"
model_path = os.path.join(script_dir, modelname + ".mo")
mod = ModelicaSystem(model_path, modelname)
mod.setParameters([
    f"R = {R}",
    f"L = {L}",
    f"C = {C}"
])
mod.setSimulationOptions('stopTime=5.0')
mod.simulate()


[t] = mod.getSolutions("time")
[u_a]= mod.getSolutions("u_a")
[i_e] = mod.getSolutions("i_e")
[i_c] = mod.getSolutions("i_c")
[i_L] = mod.getSolutions("i_L")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

# Oben: u_a und i_c
ax1.plot(t, u_a, 'k', label=r'$u_a$ (V)')
ax1.plot(t, i_c, 'r', label=r'$i_c$ (A)')
ax1.set_ylabel(r'$u_a$, $i_c$')
ax1.grid(True)
ax1.legend()

# Unten: 10*i_e und u_e
u_e = 2 * np.sin(2*np.pi*1*t)  # rekonstruiert
ax2.plot(t, 10*i_e, 'b', label=r'$10 \cdot i_e$ (A)')
ax2.plot(t, u_e, 'darkgreen', label=r'$u_e$ (V)')
ax2.set_xlabel('Zeit $t$ (s)')
ax2.set_ylabel(r'$u_e$, $10 \cdot i_e$')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()