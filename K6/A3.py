import numpy as np 
from OMPython import ModelicaSystem
import matplotlib.pyplot as plt
import os

modelname='A3'
# von Köller, nicht auf mein System anwendbar
# mod=ModelicaSystem(os.getcwd()+'/'+modelname+'.mo',modelname, [os.getcwd().replace('\\','/')+'/ModSimBib/package.mo'])

# Angepasste Version:
current_path = os.path.realpath(__file__).strip(__file__.split('/')[-1])
package_path =  current_path + 'ModSimBib/package.mo'
model_path = current_path + 'A3.mo'
mod = ModelicaSystem(model_path, modelname, [package_path])

mod.setSimulationOptions('stopTime=0.015')
mod.simulate()
[t]=mod.getSolutions('time')
[i_L]=mod.getSolutions('inductor1.i')
[u]=mod.getSolutions('vierqst1.u_out')

fig=plt.figure(1, figsize=(10,6)); fig.clf()
ax = fig.add_subplot(111)
ax.plot(t, 10*i_L, 'r')
ax.plot(t, u, 'darkgreen')
ax.grid()
ax.set_xlabel('t/s')
ax.set_ylabel('u/v bzw. 10*i/A')
ax.set_xticks(np.arange(0, 0.015+0.005, 0.005))
fig.tight_layout()

plt.show()