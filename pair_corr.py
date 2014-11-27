#!/usr/bin/python2.7
import matplotlib.pyplot as plt
import numpy as np
from vaspData import readCONTCAR
import sys
import argparse

def rdf(contcar,atom1Symbol='Si',atom2Symbol='Si',binSpacing=0.1,rMax=0,nx=1,ny=1,nz=1):

    multi = np.array([nx,ny,nz])

    latConst,latticeVecs,atomSetDirect = readCONTCAR(contcar)

    if rMax == 0:
        x = np.linalg.norm(latticeVecs[0])
        y = np.linalg.norm(latticeVecs[1])
        z = np.linalg.norm(latticeVecs[2])
        rMax = latConst*np.min([x*nx,y*ny,z*nz])/2.

    rMax = (int(rMax/binSpacing))*binSpacing

    atomPos1 =[]
    atomPos2 =[]

    for atom in atomSetDirect:
        s,pos = atom
        if s==atom1Symbol:
            atomPos1.append(pos)
        if s==atom2Symbol: 
            atomPos2.append(pos)


    volume = np.dot(latticeVecs[0],np.cross(latticeVecs[1],latticeVecs[2]))*latConst**3

    N1 = len(atomPos1)*nx*ny*nz
    N2 = len(atomPos2)*nx*ny*nz
    volume = volume*nx*ny*nz


    # multiplying atoms nx, ny, nz

    temp = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for pos in atomPos1:
                    temp.append(pos+np.array([i,j,k]))
    atomPos1 = temp

    temp = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for pos in atomPos2:
                    temp.append(pos+np.array([i,j,k]))
    atomPos2 = temp


    latticeVecs=np.array(latticeVecs)

    r = []
    for pos1 in atomPos1:
        for pos2 in atomPos2:
            diff = pos2-pos1
            if np.linalg.norm(diff) ==  0:
                continue
            diff = diff - np.round(diff/multi)*multi
            #check check check this out
            diff = np.dot(latticeVecs.T,diff)
            diff = latConst*np.linalg.norm(diff)
            r.append(diff)


    bins = np.linspace(0, rMax, num = rMax/binSpacing+1, endpoint=True)

    hist, bin_edges = np.histogram(r, bins)
    rdf = hist
    rdf = hist/4./np.pi/(bin_edges+binSpacing)[:-1]**2/N1/N2*volume/binSpacing

    return bin_edges[:-1],rdf

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("atom1Symbol", help="atom1Symbol")
    parser.add_argument("atom2Symbol", help="atom2Symbol")
    parser.add_argument("-f","--file", help="CONTCAR file path is needed. default is ./CONTCAR",default = "./CONTCAR")
    parser.add_argument("-b", "--bins", type=float, default = 0.1, help = 'bin spacing (angstrom)')
    parser.add_argument("-R", "--rMax", type=float, default = 0, help = 'maximum r, I expect you know what I mean')

    parser.add_argument("-x", "--nx",  type=int, default = 1, help = 'multiplying supercell along x')
    parser.add_argument("-y", "--ny",  type=int, default = 1, help = 'multiplying supercell along y')
    parser.add_argument("-z", "--nz",  type=int, default = 1, help = 'multiplying supercell along z')


    args = parser.parse_args()
    
    atom1Symbol = args.atom1Symbol
    atom2Symbol = args.atom2Symbol

    contcar = args.file
    binSpacing =  args.bins
    rMax = args.rMax

    nx = args.nx
    ny = args.ny
    nz = args.nz

    binSpacing, rMax, nx, ny, nz = float(binSpacing), float(rMax), int(nx), int(ny), int(nz) 
    r,rdf = rdf(contcar, atom1Symbol, atom2Symbol, binSpacing, rMax, nx, ny, nz)
    plt.plot(r,rdf)
    plt.show()   
