#!/usr/bin/python
'''
Created on 2014. 03. 21.

@author: nwan
@author: shj
'''
import numpy as np
import matplotlib.pyplot as plt
from vaspData import readPROCAR
import argparse

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

parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="PROCAR file path is needed. default is ./PROCAR",default = "./PROCAR")
parser.add_argument("-o", "--orbital", choices=['spd','dz2'], default='spd', help="how to colorize which orbital")

args = parser.parse_args()
procar_path = args.file
opt = args.orbital


PROCAR = readPROCAR(procar_path)

kpts = PROCAR[3]
eigs = PROCAR[4]
proj = PROCAR[5]

kpts =np.append([kpts[0]],kpts,axis=0)
x = [ np.linalg.norm(kpts[i] - kpts[i-1]) for i,k in enumerate(kpts)]
x = x[:-1]
x[0]=0
#x = [dist if dist != 0  else 0 for dist in x]
x = np.cumsum(x)


proj = np.sum(proj,axis=2)
proj = proj.reshape(np.prod(proj.shape)/proj.shape[-1],proj.shape[-1])

color = np.array([orbitColor(p, opt ) for p in proj])

X = np.tile(x,(1,eigs.shape[1])).flatten()
X = np.repeat(x,eigs.shape[1])
Y = eigs.flatten()




fig = plt.figure()
plt.axis([0,x.max(),-06,11])
for i in range(len(X)):
    plt.plot(X[i],Y[i],'o',c= color[i], alpha = 1.0)
plt.show()
