import numpy as np    
import matplotlib.pyplot as plt
from OMPython import ModelicaSystem
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
modelname = "A2"
model_path = os.path.join(script_dir, modelname + ".mo")
mod = ModelicaSystem(model_path, modelname)
mod.setSimulationOptions('stopTime=5.0')
mod.simulate()

# Get solution arrays
[t] = mod.getSolutions('time')
[u_a] = mod.getSolutions('u_a')
[i_L] = mod.getSolutions('i_L')

# plot u_a and i_L over time and show legend/grid
plt.plot(t, u_a, label='u_a (V)')
plt.plot(t, i_L, label='i_L (A)')
plt.xlabel('Zeit (s)')
plt.ylabel('Wert')
plt.title('RLC-Schaltkreis: Spannung und Strom (Modelica)')
plt.legend()
plt.grid()
plt.show()