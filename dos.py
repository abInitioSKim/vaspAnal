#!/usr/bin/python
from vaspData import readDOSCAR
import matplotlib.pyplot as plt
import sys
if len(sys.argv)>1:
    fileName = sys.argv[1]
else:
    fileName = 'DOSCAR'
eSet,tDOSSet,sDOSSet,pDOSSet,dDOSSet = readDOSCAR(fileName,[1,2,3])

plt.plot(eSet,tDOSSet,'r-',label='DOS_total')

plt.plot(eSet,sDOSSet,'g-',label='Si_1_s')
plt.plot(eSet,pDOSSet,'b-',label='Si_1_p')

plt.legend()

plt.show()
