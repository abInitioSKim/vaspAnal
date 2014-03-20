#!/usr/bin/python
from vaspData import readLOCPOT
import numpy as np
import matplotlib.pyplot as plt

a,[scLatVecx,scLatVecy,scLatVecz],[gridx,gridy,gridz] , LOCPOT = readLOCPOT()

im = plt.imshow(np.average(LOCPOT,axis=2))

plt.colorbar(im)

plt.show()

