import sys
sys.path.append('./..')
from CFG.settings import *
import matplotlib
matplotlib.use('PDF')
from include.utils import smartAppend
from include.utils import smartDumpDictHdf5
import limix.stats.fdr as FDR
from include.data import QtlData

import scipy as sp
import h5py
import os
import pdb
import glob
import cPickle
from optparse import OptionParser
import pylab as pl
pl.ion()

def import_hdf5(in_file):
    f = h5py.File(in_file, 'r')
    R = {}
    for key in f.keys():
        R[key] = f[key][:]
    f.close()
    return R

def qq_plot(plt, pv, color='k', label=None):
    """ qq plot """
    if type(pv) is not list:
        pv = [pv]; color = [color];
    for i, _pv in enumerate(pv):
        opy = _pv.copy()
        opy = opy[opy>0]
        opv = -sp.log10(sp.sort(opy))
        tpv = -sp.log10(sp.linspace(0,1,opv.shape[0]+1)[1:])
        pl.plot(tpv, opv, '.', color=color[i], label=label)
    xlim1, xlim2 = plt.get_xlim()
    pl.plot([0,xlim2],[0,xlim2],'r')
    pl.xlabel('Expected -log10(pv)')
    pl.ylabel('Observed -log10(pv)')

if __name__=='__main__':

    parser = OptionParser()
    parser.add_option("--peer", action="store_true", dest='peer', default=False)
    (opt, args) = parser.parse_args()
    opt_dict = vars(opt)

    data = QtlData()

    in_dir = os.path.join(CFG['out'], 'eqtl', 'eqtl_cis')
    fname = 'summary'
    if opt.peer:    fname += '_peer'
    fname+= '.hdf5'
    in_file = os.path.join(in_dir, fname)
    Rtest = import_hdf5(in_file)

    fname = 'summary'
    if opt.peer:    fname += '_peer'
    fname+= '_perm.hdf5'
    in_file = os.path.join(in_dir, fname)
    Rperm = import_hdf5(in_file)

    pdb.set_trace()

    figdir = './../figures/eqtl/cis_eqtl'
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    if 1:
        """ relative position  plkot """
        pl.figure(1, figsize=(6, 4))
        qv = Rtest['pv_bonf']
        rpos = abs(Rtest['pos']-Rtest['gene_start'])
        plt = pl.subplot(111)
        pl.hist(rpos[qv<0.10], 30)
        pl.xlabel('eQTL relative position')
        pl.tight_layout()
        pl.savefig(os.path.join(figdir, 'relpos.pdf'))
        pdb.set_trace()
        pl.close()

    if 1:
        """ relative position  plkot """
        pl.figure(1, figsize=(10, 4))
        plt = pl.subplot(121)
        qq_plot(plt, Rtest['pv_bonf'], color='b', label='real data')
        qq_plot(plt, Rperm['pv_bonf'], color='k', label='perm data')
        pl.legend(loc=2)
        plt = pl.subplot(122)
        qq_plot(plt, Rperm['pv_bonf'], color='k')
        pl.tight_layout()
        pl.savefig(os.path.join(figdir, 'qq_plot.png'))
        pdb.set_trace()
        pl.close()
