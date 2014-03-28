#!/usr/bin/python
from vaspData import *
import matplotlib.pyplot as plt
import sys

if len(sys.argv)>2:
    EIGENVAL = sys.argv[1]
    
else: 
    EIGENVAL = "EIGENVAL"

kpts, eigvals = readEIGENVAL(EIGENVAL)

kpts.append(kpts[0])
x = [ np.linalg.norm(kpts[i] - kpts[i-1]) for i,k in enumerate(kpts)]
x = x[:-1]
x = [dist if dist<0.3 else 0 for dist in x]
x = np.cumsum(x)

X = np.tile(x,eigvals.shape[1])
Y = eigvals.T.flatten()
fig = plt.figure()

plt.plot(X,Y,'ro', alpha = 0.3)
#fig.savefig(sys.argv[2])
plt.show()
