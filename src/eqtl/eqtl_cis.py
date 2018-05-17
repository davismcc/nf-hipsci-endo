import sys
sys.path.append('./../')
from CFG.settings import *
import limix
import limix.modules.qtl as QTL
import limix.stats.fdr as FDR
from include.data import QtlData
from CFG.settings import *
from include.utils import smartDumpDictHdf5
from include.utils import getLambda
import scipy as sp
import pylab as pl
import h5py
import pdb
import copy
import warnings
import os
from optparse import OptionParser

out_dir = os.path.join(CFG['out'], 'eqtl', 'eqtl_cis')

if __name__=='__main__':

    parser = OptionParser()
    parser.add_option("--n_jobs", dest='n_jobs', type=int, default=1)
    parser.add_option("--job_i", dest='job_i', type=int, default=0)
    parser.add_option("--peer", action="store_true", dest='peer', default=False)
    parser.add_option("--perm", action="store_true", dest='perm', default=False)
    parser.add_option("--seed", dest='seed', type=int, default=None)
    parser.add_option("--debug", action="store_true", dest='debug',
                      default=False)
    (opt, args) = parser.parse_args()
    opt_dict = vars(opt)

    if opt.debug:
        fname = 'debug.hdf5'
        pdb.set_trace()
    else:
        runs_folder = 'runs'
        if opt.peer:    runs_folder += '_peer'
        if opt.perm:    runs_folder += '_perm'
        split_folder = '%.4d' % int(sp.ceil(opt.job_i/1000))
        out_dir = os.path.join(out_dir, runs_folder, split_folder)
        fname = '%.3d_%.3d.hdf5' % (opt.n_jobs, opt.job_i)

    # read data and split dataset into jobs
    data = QtlData()
    K = data.get_K()
    K/= K.diagonal().mean()
    K+= 1e-4 * sp.eye(K.shape[0])
    all_genes = data.geneID.copy() 
    n_genes = all_genes.shape[0]
    Icv = sp.floor(opt.n_jobs * sp.arange(n_genes) / n_genes)
    genes = all_genes[Icv==opt.job_i]

    # creates out file
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_file = os.path.join(out_dir, fname)
    fout = h5py.File(out_file,'w')

    # loops over genes
    for gene in genes:

        print '.. Analyzing gene %s' % gene
        gene_group = fout.create_group(gene)

        print '   .. Importing data'
        try:
            Xc, info = data.getGenotypes(gene, return_info=True)
        except:
            print 'Error: no SNPs found in cis'
            continue
        Y = data.getPhenotypes(gene, peer=opt.peer, gauss=True)
        o = gene_group.create_group('snp_info')
        smartDumpDictHdf5(info, o)

        if opt.perm:
            if opt.seed is not None:
                sp.random.seed(opt.seed)
            idxs = sp.random.permutation(Xc.shape[0])
            Xc = Xc[idxs, :]

        if 1:
            print "   .. single trait analysis"
            lmm = QTL.test_lmm(Xc, Y, K=K)
            pv = lmm.getPv()
            RV = {}
            RV['pv'] = pv
            RV['qv'] = FDR.qvalues(pv) 
            RV['beta'] = lmm.getBetaSNP() 
            RV['lambda'] = getLambda(pv)
            o = gene_group.create_group('st')
            smartDumpDictHdf5(RV, o)

    fout.close()

