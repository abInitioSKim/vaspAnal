#!/usr/local/bin/python
from vasp_io import readCHGCAR, readEIGENVAL
def get_EINT(path):
    # Selected energy range (EINT)   -4.1000   4.1000
    pass

def gauss(x, a, mean, sigma):
    from scipy.stats import norm
    return a * norm.pdf(x, mean, sigma)

def get_ldos(path, extent, sigma=0.1, e_grid=100, axis=0):
    '''
    read PARCHG file in the path
    return ldos along x
    Arguments:
        extent: (x_min, x_max, energy_min, energy_max)
    '''
    ldos = []
    energy_list = []
    chg_list = []
    x_min, x_max, energy_min, energy_max = extent
    energy_space = (energy_max - energy_min) / float(e_grid)
    
    energy_bin = np.linspace(e_min, e_max, e_grid)

    # read EIGENVAL to get eigenvalue
    eigenval, n_elect, weight = readEIGENVAL('{}/EIGENVAL'.format(path),
                                             NELECT=True, lweight=True)
    # energy w.r.t e_fermi
    e_fermi = eigenval[0, n_elect / 2 -1]

    # number of eigenvalues within the range of energy
    num_eig = np.sum([1 for energy in eigenval.flat
                      if e_min <= energy - e_fermi <= e_max])

    print num_eig, 'PARCHG will be read'
    num = 0
    axis_sum = [0, 1, 2]
    axis_sum.pop(axis)
    print axis_sum
    for kpt_index, eigenval_kpt in enumerate(eigenval):
        # if kpt_index > 0:
        #     continue
        for band_index, energy in enumerate(eigenval_kpt):
            energy -= e_fermi
            if e_min <= energy <= e_max:
                num += 1
                parchg_file = '/PARCHG.{:0>4d}.{:0>4d}'.format(
                                band_index + 1, kpt_index + 1)
                print num, num_eig, parchg_file
                if os.path.isfile(path + '/' + parchg_file):
                    lat_const, lattice_matrix, CHGCAR = \
                        readCHGCAR(path + parchg_file)


                    chg_sum = np.sum(CHGCAR, axis=tuple(axis_sum))
                    x_grid = len(chg_sum)
                    x_space = (x_max - x_min) / float(x_grid)
                    chg_sum *= weight[kpt_index] / 2. / np.prod(CHGCAR.shape) \
                               / energy_space / x_space

                    ldos.append([energy, chg_sum])
                else:
                    print parchg_file, "does not exists. check energy range",\
                          "or your parchg calculations"

    x_max = np.linalg.norm(lattice_matrix[axis]) * lat_const

    # smoothing charge densities along energy axis
    ldos_interpol = np.zeros((x_grid, e_grid))
    for x_index in xrange(x_grid):
        temp_x = np.zeros(e_grid)
        for energy_index in xrange(len(ldos)):
            energy = ldos[energy_index][0]
            density = ldos[energy_index][1][x_index]
            energy_list.append(energy)
            temp_x += gauss(energy_bin, density, energy, sigma)

        ldos_interpol[x_index, :] = temp_x[::-1]
    return ldos_interpol

def draw_ldos(ldos, extent_calc, extent_show, 
              vmin=0., vmax=1., ticks=None, log_scale=False):
    from matplotlib import rcParams, gridspec
    from matplotlib.colors import LogNorm
    import mpl_toolkits.axisartist as axisartist
    aspect = 20

    vmax = ldos.max() * vmax
    vmin = ldos.max() * vmin
    # eMin , eMax = extent[2:]

    fig = plt.figure(figsize=(7, 5))
    gs = gridspec.GridSpec(1, 2,
        width_ratios=[1, 0.1])

    # ax_ldos = fig.add_subplot(axisartist.Subplot(fig,gs[:,:1]))
    ax_ldos = plt.subplot2grid((1,aspect), (0,0), colspan=aspect-1)#(111)
    ax_cbar = plt.subplot2grid((1,aspect), (0,aspect-1))#(111)
    ax_ldos.set_ylabel('Energy (eV)')
    ax_ldos.set_xlabel('x ($\AA$)') 
    print vmin, vmax
    if log_scale:
        im = ax_ldos.imshow(ldos.T, aspect='auto', 
               extent=extent_calc, norm=LogNorm(vmin=vmin, vmax=vmax))
    else:
        im = ax_ldos.imshow(ldos.T, aspect='auto', 
               extent=extent_calc, vmin=vmin+0.01, vmax=vmax)

    ax_ldos.axis(extent_show)

    # draw color bar
    cbar = plt.colorbar(im, cax=ax_cbar)
    ax_cbar.yaxis.tick_right()
    ax_cbar.set_ylabel('Local Density of states (States/eV/$\AA$)')
    # colorbar ticks
    if ticks is not None:
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(ticks)

    plt.show()

if __name__ == '__main__':
    import sys, os
    import argparse
    import os
    import numpy as np
    from pylab import *
    print 'Draw local density of states'

    # argument handle
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--path",
        help="PARCHG and EIGENVAL file path is needed.",
        required=True)
    parser.add_argument("-o","--output",
        help="Output file name also read it if exists. default is ldos.npy",
        default = "ldos.npy")
    parser.add_argument("-g","--egrid",
        help="number of grids for energy. default is 200",
        default = 200)
    parser.add_argument("-e","--energy", nargs=2,
        help="energy range with respect to Fermi level. default is -3 3",
        default = [-3, 3])
    parser.add_argument("-x","--xrange", nargs=2,
        help="x range with respect to Fermi level. default is 0 40",
        default = [0, 40])
    parser.add_argument("-v","--vrange", nargs=2,
        help="x range with respect to Fermi level. default is 0 1.0",
        default = [0.1, 1.0])
    parser.add_argument("-s","--sigma",
        help="Gaussian with for energy smearing. default is 0.1 eV",
        default = 0.2)
    parser.add_argument("-a","--axis",
        help="summation axis (0, 1, 2). default is 0",
        default = 0)
    parser.add_argument("--force-calc", dest='f_calc',
        help="force to recalculation. default is False",
        default = False, action='store_true')
    parser.add_argument("--log-scale", dest='log_scale',
        help="scale log. default is True",
        default = False, action='store_true')

    args = parser.parse_args()
    path = args.path
    e_min, e_max = [float(energy) for energy in args.energy]
    x_min, x_max = [float(x) for x in args.xrange]
    vmin, vmax = [float(v) for v in args.vrange]
    e_grid = int(args.egrid)
    sigma = float(args.sigma)
    f_calc = args.f_calc
    save_name = args.output
    axis = int(args.axis) % 3
    log_scale = args.log_scale

    extent_calc = (x_min, x_max, e_min - sigma, e_max + sigma)
    extent_show = (x_min, x_max, e_min, e_max)

    print f_calc
    if f_calc:
        ldos_interpol = get_ldos(path, extent_calc, sigma=sigma,
                                 e_grid=e_grid, axis=axis)
        np.save(path + '/' + save_name, ldos_interpol)
    else:
        try:
            print save_name, ' exist. read it'
            ldos_interpol = np.load(path + '/' + save_name)
        except IOError:
            print 'no save file. calculation starts'
            ldos_interpol = get_ldos(path, extent_calc, sigma=sigma,
                                     e_grid=e_grid, axis=axis)
        np.save(path + '/' + save_name, ldos_interpol)

    draw_ldos(ldos_interpol, extent_calc, extent_show,
              vmin, vmax, np.arange(10) * 10,  log_scale=log_scale)
