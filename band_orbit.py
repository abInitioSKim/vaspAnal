#!/usr/bin/python
'''
Created on 2014. 03. 21.

@author: nwan
@author: shj
'''
import numpy as np
import matplotlib.pyplot as plt
from vaspData import readPROCAR

def orbitColor(proj, opt='spd'):
    if opt == 'spd':
        s = proj[0]
        if len(proj)>4:
            p = np.sum(proj[1:4])
            d = np.sum(proj[4:9])
        else :
            p = proj[1]
            d = proj[2]  
        c = np.array([s,p,d])/(np.max([s,p,d]))
        return c
    elif opt == 'dz2':
        dz2 = proj[6]/np.sum(proj[:-1])
        return (1-dz2,1,1)


#PROCAR = readPROCAR('PROCAR_MoS2')
PROCAR = readPROCAR('PROCAR_Si')
#PROCAR = readPROCAR('/home/users/zeneco/00_Si_direct/74_new5/01_368/05_Band/PROCAR')

kpts = PROCAR[3]
eigs = PROCAR[4]
proj = PROCAR[5]

kpts =np.append([kpts[0]],kpts,axis=0)
x = [ np.linalg.norm(kpts[i] - kpts[i-1]) for i,k in enumerate(kpts)]
x = x[:-1]
x[0]=0
#x = [dist if dist != 0  else 0 for dist in x]
x = np.cumsum(x)


#Proj = np.zeros((nKpt,nBands,nIons,nOrbits))
proj = np.sum(proj,axis=2)
proj = proj.reshape(np.prod(proj.shape)/proj.shape[-1],proj.shape[-1])

#proj = np.array([orbitColor(p, 'dz2' ) for p in proj])
color = np.array([orbitColor(p, 'spd' ) for p in proj])
print np.sum(color,axis=0)

X = np.tile(x,(1,eigs.shape[1])).flatten()
X = np.repeat(x,eigs.shape[1])
Y = eigs.flatten()




fig = plt.figure()
plt.axis([0,x.max(),-06,11])
for i in range(len(X)):
    plt.plot(X[i],Y[i],'o',c= color[i], alpha = 1.0)
plt.show()
