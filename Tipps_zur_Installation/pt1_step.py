from OMPython import ModelicaSystem
import matplotlib.pyplot as plt
import os

modelname='pt1_step'
mod=ModelicaSystem(os.getcwd()+'/'+modelname+'.mo', modelname)
mod.setSimulationOptions('stopTime=5.0')
#mod.setParameters('firstOrder.k=2.0')
mod.simulate()

[t]=mod.getSolutions('time')
[phi]=mod.getSolutions('firstOrder.y')

fig=plt.figure(1, figsize=(10,6)); fig.clf()
plt.plot(t, phi, 'b')
plt.xlabel('t')
plt.grid()
plt.show()