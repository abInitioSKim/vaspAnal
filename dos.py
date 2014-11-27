#!/usr//bin/python2.7
from vaspData import readDOSCAR, writeData
import matplotlib.pyplot as plt
import sys
import argparse

if __name__=='__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-f","--file", help="DOSCAR file path is needed. default is ./DOSCAR",default = "./DOSCAR")
    parser.add_argument('-a', '--atoms', nargs='+', type=int, help="atom index for pdos. ex) 0: all atoms, averaged. [1, 2, 3] ")
    parser.add_argument("-e", "--eRange", nargs=2, type=float, help = 'energy Range (eV)')
    parser.add_argument("-w","--write", help="write data",default = None)

    args = parser.parse_args()
    doscar = args.file
    atoms = args.atoms
    l_write = args.write

    eSet,tDOSSet,sDOSSet,pDOSSet,dDOSSet = readDOSCAR(doscar,atoms)

    plt.plot(eSet,tDOSSet,'r-',label='DOS_total')

    if atoms !=None:
        plt.plot(eSet,sDOSSet,'g-',label=str(atoms)+'_s')
        plt.plot(eSet,pDOSSet,'b-',label=str(atoms)+'_p')
        plt.plot(eSet,dDOSSet,'y-',label=str(atoms)+'_d')

    plt.legend()
    
    if args.eRange != None:
        eMin, eMax = args.eRange
        plt.xlim(eMin, eMax)    
    
    if l_write:
        writeData('./energy.dat', eSet, eSet.shape)
        writeData('./dos_total.dat', tDOSSet, tDOSSet.shape)
        writeData('./pdos_s.dat', sDOSSet, sDOSSet.shape)
        writeData('./pdos_p.dat', pDOSSet, pDOSSet.shape)
        writeData('./pdos_d.dat', dDOSSet, dDOSSet.shape)
    plt.show()
