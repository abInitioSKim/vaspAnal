#!/usr/local/Python-2.7.3/bin/python2.7

import matplotlib.pyplot as plt
import numpy as np
import sys, itertools
from matplotlib import rcParams, gridspec
from matplotlib import rc
import mpl_toolkits.axisartist as axisartist
from vaspData import readKPOINTS_linemode, readDOSCAR, readPROCAR, readCONTCAR, getVBM, getEgap
if __name__ == '__main__':
    #rc('font',**{'family':'sans-serif','sans-serif':['Times New Roman']})
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Times New Roman']})
    # use usetex=True for publish quality
    # rc('text', usetex=False)
    #rc('text', usetex=True)
    params = {'backend': 'ps',
          'text.latex.preamble': [r"\usepackage{upgreek}",
                                  r"\usepackage{siunitx}",
                                  r"\usepackage{amsmath}",
                                  r"\usepackage{amstext}",],
          'axes.labelsize': 16,
          #'axes.linewidth': 1,
          #'text.fontsize':17,
          'legend.fontsize': 14,
          'xtick.labelsize': 16,
          # 'xtick.major.pad'      : 6,      # distance to major tick label in points
          # 'xtick.minor.pad'      : 6,      # distance to the minor tick label in points
          'ytick.major.pad'      : 6,     # distance to major tick label in points
          'ytick.minor.pad'      : 6,     # distance to the minor tick label in points
          #'xtick.major.width' : 0.75,
          'ytick.labelsize': 16,
          # 'figure.figsize': [8.8,6.8],
          # 'text.usetex': True,
          'axes.unicode_minus': True,
          'ps.usedistiller' : 'xpdf'}          
    rcParams.update(params)
    # rcParams.update({'figure.autolayout':True})


def readBand_KPOINTS(fileName ='EIGENVAL', NELECT = 0, kpoints_file = 'KPOINTS'):
    # read EIGENVAL
    f=open(fileName)
    buffer=f.readlines()
    f.close()
    [nElect,nKpt,nBand]=[int(i) for i in buffer[5].split()][:]
    #print [nElect, nBand,nKpt]
    # VBM = getVBM(nElect/2, fileName)
    
    bandInfo = []
    kpoints =[]
    eigenvals =np.zeros((nBand,nKpt))
    
    for j in range(nKpt):
        kpoint =np.array(buffer[-1 + 8 + (nBand+2)*j].split())[:3]
        kpoint = np.array([float(k) for k in kpoint])
        kpoints.append(kpoint)


        for i in range(nBand):
            eigenval = buffer[i + 8 + (nBand+2)*j].split()
            eigenval = float(eigenval[1])
            eigenvals[i,j] = eigenval
    
    nkpt_line, specialKPName = readKPOINTS_linemode(fileName = kpoints_file)
    specialKPPos = [i for i, k in enumerate(kpoints) if i == 0 or (i+1) % nkpt_line == 0]
    #print specialKPName 
    VBM = eigenvals[nElect/2-1,:].max()
    VBM = getVBM(fileName)

    return kpoints,eigenvals - VBM , specialKPPos, specialKPName


def drawBand(eigenval='EIGENVAL', isLine=True, ax=''):
    if ax=='':ax=plt
    
    kpts, eigvals, specialKPPos, specialKPName = readBand_KPOINTS(eigenval)
    # print specialKPPos, specialKPName
    kpts.append(kpts[0])
    x = [ np.linalg.norm(kpts[i] - kpts[i-1]) for i,k in enumerate(kpts)]
    x = x[:-1]
    x = [dist if dist<0.3 else 0 for dist in x]
    x = np.cumsum(x)

    for i in specialKPPos:
        ax.axvline(x=x[i],color = 'k')

    if isLine == True:
        for band in eigvals:
            ax.plot(x,band,'k-')
    else:
        X = np.tile(x,eigvals.shape[0],)
        Y = eigvals.flatten()
        plt.plot(X,Y,'ro', alpha = 0.3, )

    ax.set_xticks(x[specialKPPos])
    ax.set_xticklabels(specialKPName)
    ax.set_xlim((x.min(),x.max()))


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

def kpt_length(kpts):
    latConst, latticeVecs, atomSetDirect = readCONTCAR(fileName='CONTCAR')
    a, b, c = ([lat for lat in latticeVecs])
    G_a = 2 * np.pi * np.cross(b, c) / np.dot(c, np.cross(a, b))
    G_b = 2 * np.pi * np.cross(c, a) / np.dot(c, np.cross(a, b))
    G_c = 2 * np.pi * np.cross(a, b) / np.dot(c, np.cross(a, b))
    scale = np.array([np.linalg.norm(v) for v in [G_a, G_b, G_c]])
    kpts = [kpt * scale for kpt in kpts]
    kpts.append(kpts[0])
    x = [ np.linalg.norm(kpts[i] - kpts[i-1]) for i,k in enumerate(kpts)]
    x = x[:-1]
    x = [dist if dist < 0.3 * np.min(scale) else 0 for dist in x]
    x = np.cumsum(x)  
    return x


def draw_band_orbit(procar = 'PROCAR', eigenval = 'EIGENVAL', isLine = True, ax=''):
    if ax=='':ax=plt
    
    kpts, eigvals, specialKPPos, specialKPName = readBand_KPOINTS(eigenval)
    # kpts.append(kpts[0])
    # x = [ np.linalg.norm(kpts[i] - kpts[i-1]) for i,k in enumerate(kpts)]
    # x = x[:-1]
    # x = [dist if dist<0.3 else 0 for dist in x]
    # x = np.cumsum(x)
    x = kpt_length(kpts)

    opt = 'spd'
    PROCAR = readPROCAR(procar)
    proj = PROCAR[5]
    proj = np.sum(proj,axis=2)
    proj = proj.reshape(np.prod(proj.shape)/proj.shape[-1],proj.shape[-1])

    color = np.array([orbitColor(p, opt ) for p in proj])

    Egap = getEgap(eigenval)
    if Egap > 0:
        ax.axhspan(0, Egap, lw = 0, facecolor = '0.5', alpha = 0.2)
    for i in specialKPPos:
        ax.axvline(x=x[i],color = 'k')

    if isLine:
        for band in eigvals:
            ax.plot(x,band,'k-')
    elif True:
        # print color.shape
        # print eigvals.shape
        for i, band in enumerate(eigvals):
            ax.plot(x,band,'o', c = color[i], markeredgewidth = 0, alpha = .5)
        for i, band in enumerate(eigvals):
            ax.plot(x,band,'o', c = color[i], markeredgewidth = 0, alpha = .5)
        # for i, band in enumerate(eigvals[::-1]):
        #     ax.plot(x,band,'o', c = color[-i -1], markeredgewidth = 0, alpha = 1)
        # for i, band in enumerate(eigvals[::-1]):
        #     ax.plot(x,band,'o', c = color[-i - 1], markeredgewidth = 0, alpha = 0.5)
    else:
        X = np.tile(x,eigvals.shape[0],)
        Y = eigvals.flatten()
        plt.plot(X,Y,'ro', alpha = 0.3, )


    ax.set_xticks(x[specialKPPos])
    ax.set_xticklabels(specialKPName)
    ax.set_xlim((x.min(),x.max()))



def drawDOS(doscar="",atomNum=0,ax='',symbol='-',dos='tspd',legend=[],factor=1.):
    if ax=='':ax=plt
    #eSet,tDOSSet,sDOSSet,pDOSSet,dDOSSet = readPDOS(doscar,atomNum)
    eSet,tDOSSet,sDOSSet,pDOSSet,dDOSSet = readDOSCAR(doscar,atomNum)
    # plt.plot(eSet,sDOSSet)
    # print eSet.shape, sDOSSet.shape
    mozala = len(dos) - len(legend)
    if mozala>0:
        # legend += [''] *mozala
        legend = [s for s in dos]
    legend = legend[::-1]
    if 't' in dos:
        ax.plot(tDOSSet*factor,eSet,symbol,label=legend.pop())
    if 's' in dos:
        ax.plot(sDOSSet*factor,eSet,symbol,label=legend.pop()) #Si$_{}$
    if 'p' in dos:
        ax.plot(pDOSSet*factor,eSet,symbol,label=legend.pop()) #Si$_{}$
    if 'd' in dos:
        ax.plot(dDOSSet*factor,eSet,symbol,label=legend.pop()) #Si$_{}$

    handles, labels = ax.get_legend_handles_labels()
    
    max_dos = (np.average(tDOSSet) * 15 / 10.) * factor 

    ax.set_xlim((0, max_dos ))

    leg=ax.legend(handles, labels, \
    frameon=False ,loc='best',numpoints=1,handletextpad=0.5,borderpad=0.2, labelspacing=0.2, handlelength=1 )
    # for t in leg.get_texts():
    #     t.set_fontsize(legendFontSize)

if __name__ == '__main__':
# if False:
    ''' ================================================ setting ================================================ '''
    fig = plt.figure(figsize=(10,7))
    plt.subplots_adjust(top=0.92,bottom=0.10, left =0.12 ,right =0.97, wspace=0.0,hspace=0.0)
    gs = gridspec.GridSpec(2, 2,
        width_ratios=[2,1],
        height_ratios=[0.05,1])
    eMin , eMax = -15, 10

    EigenvalFile = "EIGENVAL"
    DoscarFile = "DOSCAR"
    
    ''' ================================================ draw band  ==================================== '''
    ax_band = fig.add_subplot(axisartist.Subplot(fig,gs[:,:1]))
    ax_band.set_ylim((eMin , eMax))
    ax_band.set_ylabel('Energy (eV)') 

    # drawBand(EigenvalFile,ax=ax_band)
    draw_band_orbit(procar = 'PROCAR', eigenval = EigenvalFile, isLine = False, ax = ax_band)

    ''' ================================================ draw dos  ===================================='''
    ax_dos = fig.add_subplot(gs[:,1:2])
    ax_dos.tick_params(direction='in', labelleft='off',labelright='off', labelbottom = 'off' )
    ax_dos.set_ylim((eMin , eMax))
    ax_dos.set_title('DOS')

    drawDOS(DoscarFile,0, dos = 'tspd', ax= ax_dos)

    # fig.savefig('Estruct16.eps')
    # fig.savefig('Estruct16.pdf')
    # fig.savefig('Estruct130.png')
    plt.show()
