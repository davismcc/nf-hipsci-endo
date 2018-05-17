import scipy as SP
import io
import sys
import os
import pdb
import re
import time
import glob
import string
import time
import h5py
import re
import cPickle
from sklearn import metrics



ct = 'iPS'
base_dir = '/homes/buettner/research/users/buettner/hipsci-singlecell/data/pilot3/'+ct+'/variantsAll'
counts_file = base_dir+'/counts.hdf5'
gene_file = h5py.File(counts_file,'r')['genes']
genes_c = gene_file['gene_ids'][:]

fvarTrue = h5py.File(os.path.join(base_dir,'variants.hdf5'),'r')
geno_id_map_true = fvarTrue['geno_id_map'][:]
geno_id = fvarTrue['geno_id'][:]
sample_ids = fvarTrue['sample_ids'][:]

Ngenes = len(genes_c)
topGenes = SP.logspace(1,SP.log2(10000),num=10,base=2.0, dtype='int')

#loop through all samples
geno_map_list = []
for cell in sample_ids:
    fvar = h5py.File(os.path.join(base_dir,cell+'.bam.hdf5'),'r')
    genes = fvar['gene_ids'][:]
    error = fvar['error']
    mapped = fvar['mapped']
    Itop = fvar['Itop'][genes[0]]
    count = fvar['count']

    gene_idx = SP.hstack(SP.array([SP.where(genes[i]==genes_c[Itop])[0] for i in range(len(genes))]))
    
    geno_map = []
    error_list = []
    for nGene in topGenes:
        if min(gene_idx)>nGene:
            error_list.append(SP.ones(4))
            geno_map.append("None") 
            continue
        idx_top = SP.where(gene_idx<=nGene)[0][-1]
        _genes = SP.array(genes)[0:idx_top+1]
        _error = SP.vstack(SP.array([error[g][:] for g in _genes]))
        _mapped = SP.concatenate(SP.array([mapped[g][:] for g in _genes]))
        _count = SP.concatenate(SP.array([count[g][:] for g in _genes]))
        _Iok = (_mapped>0.5) & (_count>100)        
        _error = _error[_Iok].sum(axis=0)/_Iok.sum()
        error_list.append(_error)
        _geno_id_map = _error.argmin()
        geno_map.append(geno_id[_geno_id_map])
    geno_map_list.append(geno_map)        
    print "processed sample ", cell            

geno_mapArray = SP.vstack(SP.array(geno_map_list))
f1_score = []
for nGene in topGenes:
   f1_score.append(metrics.f1_score(geno_mapTrue, geno_mapArray[, average="macro"))

#f1_scores = metrics.f1_score(l, pred, average=None)
#f1_score_av = metrics.f1_score(labs, pred, average="macro")

