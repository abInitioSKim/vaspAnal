#!/usr/local/Python-2.7.3/bin/python2.7
from vaspData import *
import matplotlib.pyplot as plt
import sys
import argparse

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", help="eigenval file path is needed. default is ./eigenval",default = "./EIGENVAL")
    args = parser.parse_args()
    eigenval_path = args.file

    kpts, eigvals = readEIGENVAL(eigenval_path)

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
