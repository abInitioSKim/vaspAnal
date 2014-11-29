#!/usr/bin/python2.7
'''
Created on 2014. 03. 21.

@author: nwan
@author: shj
'''
import numpy as np
import matplotlib.pyplot as plt
from vasp_io import readPROCAR, readKPOINTS_linemode
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

def get_efermi(PROCAR):
    occs = PROCAR[6].flatten()
    eigs = PROCAR[4].flatten()
    eigs = [eig for idx, eig in enumerate(eigs)
            if occs[idx] > 0]
    return max(eigs)

def label_ticks(kpts, specialKPName, nkpt_line):
    assert len(kpts) % nkpt_line ==0, 'check KPOINTS and PROCAR files'

    ticks =[]
    for indx in range(len(kpts) / nkpt_line):
        # print 'kpt', kpts[indx * nkpt_line]
        ticks.append(kpts[indx * nkpt_line])
    ticks.append(kpts[-1])
    return ticks

def get_kpt_length(kpt_vec, nkpt_line):
    kpts =np.append([kpt_vec[0]], kpt_vec, axis=0)
    x = [np.linalg.norm(kpts[i] - kpts[i-1]) for i, k in enumerate(kpts)]
    x = x[1:]
    for indx in range(len(x)):
        if indx % nkpt_line == 0:
            print indx
            x[indx] = 0
    x = np.cumsum(x)
    return x

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="PROCAR file path is needed. default is ./PROCAR",default = "./PROCAR")
    parser.add_argument("-o", "--orbital", choices=['spd','dz2'], default='spd', help="how to colorize which orbital")

    args = parser.parse_args()
    procar_path = args.file
    opt = args.orbital

    # procar_path = '/home/users/nwan/02Project/16_MX2HETERO/slab_single/MoS2WS2/MoS2/01_BAND/PROCAR'
    # kpoints_path = '/home/users/nwan/02Project/16_MX2HETERO/slab_single/MoS2WS2/MoS2/01_BAND/KPOINTS'
    procar_path = './PROCAR'
    kpoints_path = './KPOINTS'

    PROCAR = readPROCAR(procar_path)
    e_fermi = get_efermi(PROCAR)

    nkpt_line, specialKPName = readKPOINTS_linemode(fileName=kpoints_path)

    kpt_vec = PROCAR[3]
    eigs = PROCAR[4] - e_fermi
    proj = PROCAR[5]

    # kpts =np.append([kpts[0]],kpts, axis=0)
    # x = [np.linalg.norm(kpts[i] - kpts[i-1]) for i, k in enumerate(kpts)]
    # x = x[:-1]
    # for indx in range(len(x)):
    #     if indx % nkpt_line == 0:
    #         x[indx] = 0
    # x = np.cumsum(x)
    x = get_kpt_length(kpt_vec, nkpt_line)

    ticks = label_ticks(x, specialKPName, nkpt_line)

    proj = np.sum(proj,axis=2)
    proj = proj.reshape(np.prod(proj.shape)/proj.shape[-1],proj.shape[-1])

    color = np.array([orbitColor(p, opt ) for p in proj])

    X = np.tile(x,(1,eigs.shape[1])).flatten()
    X = np.repeat(x,eigs.shape[1])
    Y = eigs.flatten()

    fig = plt.figure()
    plt.axis([0,x.max(),-06,11])
    for i in range(len(X)):
        plt.plot(X[i],Y[i],'o', c=color[i], markeredgecolor=color[i])

    # print label_ticks(X, specialKPName, nkpt_line)
    # print X
    plt.xticks(ticks, specialKPName)

    plt.show()
