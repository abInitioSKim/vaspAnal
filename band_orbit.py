#!/usr/bin/python2.7
'''
Created on 2014. 03. 21.

@author: nwan
@author: shj
'''
import numpy as np
import matplotlib.pyplot as plt
from vasp_io import readPROCAR, readKPOINTS_linemode, readCONTCAR, get_reciprocal_lattice
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
        c = np.array([s,p,d]) #/(np.sum([s,p,d]))
        return c
    elif opt == 'dz2':
        dz2 = proj[6]/np.sum(proj[:-1])
        return (1-dz2,1,1)

def get_efermi(PROCAR):
    CUTOFF_OCC = 0.5
    occs = PROCAR[6].flatten()
    eigs = PROCAR[4].flatten()
    eigs = [eig for idx, eig in enumerate(eigs)
            if occs[idx] > CUTOFF_OCC]
    return max(eigs)

def label_ticks(kpts, specialKPName, nkpt_line):
    assert len(kpts) % nkpt_line ==0, 'check KPOINTS and PROCAR files'

    ticks =[]
    for indx in range(len(kpts) / nkpt_line):
        # print 'kpt', kpts[indx * nkpt_line]
        ticks.append(kpts[indx * nkpt_line])
    ticks.append(kpts[-1])
    return ticks

def get_kpt_length(kpt_vec, nkpt_line, rec_mat):
    for index, kpt in enumerate(kpt_vec):
        kpt_vec[index] = np.dot(kpt, rec_mat)

    kpts =np.append([kpt_vec[0]], kpt_vec, axis=0)
    x = [np.linalg.norm(kpts[i] - kpts[i-1]) for i, k in enumerate(kpts)]
    x = x[1:]
    for indx in range(len(x)):
        if indx % nkpt_line == 0:
            x[indx] = 0
    x = np.cumsum(x)
    return x

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--path",
                        help="dir path is needed. default is ./",default = "./")
    parser.add_argument("-o", "--orbital", choices=['spd','dz2'], default='spd', help="how to colorize which orbital")
    parser.add_argument("-e","--energy", nargs=2,
        help="energy range with respect to Fermi level. default is -3 3",
        default = [-3, 3])

    args = parser.parse_args()
    path = args.path
    opt = args.orbital
    e_min, e_max = [float(energy) for energy in args.energy]

    # procar_path = '/home/users/nwan/02Project/16_MX2HETERO/slab_single/MoS2WS2/MoS2/01_BAND/PROCAR'
    # kpoints_path = '/home/users/nwan/02Project/16_MX2HETERO/slab_single/MoS2WS2/MoS2/01_BAND/KPOINTS'
    # contcar_path = '/home/users/nwan/02Project/16_MX2HETERO/slab_single/MoS2WS2/MoS2/01_BAND/CONTCAR'
    procar_path = '{}./PROCAR'.format(path)
    kpoints_path = '{}./KPOINTS'.format(path)
    contcar_path = '{}./CONTCAR'.format(path)
    # read files
    rec_mat = get_reciprocal_lattice(contcar_path)
    latConst, latticeVecs, atomSetDirect = readCONTCAR(contcar_path)
    elements = list(set([atom[0] for atom in atomSetDirect]))

    n_elements = []
    for element in elements:
        n_ele = len([atom for atom in atomSetDirect if atom[0] == element])
        n_elements.append(n_ele)
    print elements
    print n_elements

    PROCAR = readPROCAR(procar_path)
    e_fermi = get_efermi(PROCAR)

    nkpt_line, specialKPName = readKPOINTS_linemode(fileName=kpoints_path)

    kpt_vec = PROCAR[3]
    eigs = PROCAR[4] - e_fermi
    proj = PROCAR[5]

    x = get_kpt_length(kpt_vec, nkpt_line, rec_mat)
    
    n_ax = len(elements)
    fig = plt.figure()

    a_i = np.cumsum(n_elements)
    a_i = np.append([0], a_i)

    ax = []

    for ele_index in range(n_ax):
        # coloring depends on the projections
        st_index = int(a_i[ele_index])
        end_index = int(a_i[ele_index + 1])

        # separate element wise
        # Proj[kpt,band,ion,:] = orbital_proj
        proj_element = np.sum(proj[:, :, st_index: end_index, :], axis=2)
        proj_all =  np.sum(proj[:, :, : , :], axis=2)

        proj_all = proj_all.reshape(
            np.prod(proj_all.shape) / proj_all.shape[-1],
            proj_all.shape[-1])

        proj_element = proj_element.reshape(
            np.prod(proj_element.shape) / proj_element.shape[-1],
            proj_element.shape[-1])

        color = np.array([orbitColor(p / np.sum(proj_all[i][:-1]), opt, ) 
            for i, p in enumerate(proj_element)])

        if color.max() > 1.:
            print color.max()
            color /= color.max()

        # draw graph
        X = np.tile(x,(1,eigs.shape[1])).flatten()
        X = np.repeat(x,eigs.shape[1])
        Y = eigs.flatten()

        # make axis
        ax.append(fig.add_subplot(1, n_ax, ele_index))
        ax[ele_index].axis([0, x.max(), e_min, e_max])
        ticks = label_ticks(x, specialKPName, nkpt_line)
        ax[ele_index].set_xticks(ticks)
        ax[ele_index].set_xticklabels(specialKPName)
        ax[ele_index].set_title(elements[ele_index])

        # draw band
        for i in range(len(X)):
            if e_min < Y[i] < e_max:
                ax[ele_index].plot(X[i],Y[i],'o',
                    c=color[i], markeredgecolor=color[i])
        # draw fermi
        ax[ele_index].plot([0, x.max()], [0, 0],
            'k-', lw=2, alpha=0.5)
        

    plt.show()
