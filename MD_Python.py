import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
from matplotlib.animation import FuncAnimation

with open("Data_forces.txt","r") as f1:
    lines1 = f1.readlines()
   
fx = []
fy = []
fz = []

for line in lines1:
    fx_ = (((line.split('\t'))[0]))
    fy_ = (((line.split('\t')))[1])
    fz_ = (((line.split('\t'))[2]))
    fx.append(float(fx_))
    fy.append(float(fy_))
    fz.append(float(fz_))

N = []
for i in range(len(fy)):
    y = i
    N.append(1+y)

plt.scatter(N, fx)
plt.xlabel('particle number')
plt.ylabel('fx/fx_lammps relation')
plt.show()

plt.scatter(N, fy)
plt.xlabel('particle number')
plt.ylabel('fy/fy_lammps relation')
plt.show()

plt.scatter(N, fz)
plt.xlabel('particle number')
plt.ylabel('fz/fz_lammps relation')
plt.show()


counter2 = 0
y2 = []
with open("Energy_K 10648.txt","r") as f:
    for line in f.readlines():
        y2.append(float(line.split(' ')[0]))
        counter2 = counter2 + 1
x2 = np.arange( 0, counter2,  1)
y2 = np.array(y2)
plt.scatter(x2,y2)
plt.xlabel('iterition step')
plt.ylabel('Kinetic Energy ()')
plt.show()

counter1 = 0
y = []
with open("Energy 10648.txt","r") as f:
    for line in f.readlines():
       y.append(float(line.split(' ')[0]))
       counter1 = counter1 + 1
x = np.arange( 0, counter1,  1) 
y = np.array(y)
plt.scatter(x,y)
plt.xlabel('iterition step')
plt.ylabel('Energy ()')
plt.show()



fig = plt.figure()
with open("List 10648.txt","r") as f:
    st9= f.readlines()
x = np.arange( 0, 50,  1) 
def top(i): 
    y = []
    vm_mid = float(st9[i].split(' ')[51])
    E_k = float(st9[i].split(' ')[52])
    N = float(st9[i].split(' ')[53])
    for j in range(50):
        h = j+1
        y.append((float(st9[i].split(' ')[h]) / 100 / N))
    y = np.array(y)
    ax = plt.axes(ylim=(0,0.1))
    plt.bar(x, y)
    te = np.arange(0, 50, 1)
    y3 = []   
    for i in range(50):
        v =  float(vm_mid / 50*(i+1))
        kT = float(2/3 * E_k / N)
        w =  float(4 * 3.14 * v * v  * (2 * 3.14 * kT )**(-1.5) * math.exp(-v*v/(2*kT)) /(50/vm_mid))
        p = (w)
        y3.append(p)  
    plt.plot(te, y3, 'red' , linewidth = 2)

hist_animation = animation.FuncAnimation(fig, top, interval = 1, cache_frame_data=False)
plt.show()


