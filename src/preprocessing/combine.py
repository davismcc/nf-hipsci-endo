"""combine all the steps of the pipeline in a single hdf5 file"""
import scipy as SP
import h5py 
import os
import glob
import sys
import re
import pdb

def calc_norm(group):
    """calc normalized counts (library size standardized) and store
    Iok: set of libraries to consider for normalization  
    """
    counts = SP.array(group['counts'][:],dtype='float')
    #library size in million reads
    L      = counts.sum(axis=0)
    LM = SP.median(L)
    C  = LM/L
    counts*=C
    #simple normaliation by library size
    group.create_dataset(name='counts_norm',data=counts)

if __name__ == '__main__':
    #1. in directory
    in_dir = sys.argv[1]
    fout = h5py.File(os.path.join(in_dir,'data_raw.hdf5'),'w')
    fcount = h5py.File(os.path.join(in_dir,'counts','counts.hdf5'),'r')
    fvar = h5py.File(os.path.join(in_dir,'variants','variants.hdf5'),'r')
    quality = SP.loadtxt(os.path.join(in_dir,'data_quality.csv'),delimiter=',',dtype='str')
    
    #reference sample ID from gene expression data
    SAMPLE_IDS = SP.sort(fcount['genes/sample_ids'][:])

    #1. summarize gene expression
    sample_ids = fcount['genes/sample_ids'][:]
    counts    = fcount['genes/counts'][:]
    gene_ids  = fcount['genes/gene_ids'][:]
    annotation = fcount['genes/annotation'][:]
    #remove samples with no counts or library sizes below threshold:
    L = counts.sum(axis=0)

    #fix sample_ids
    Imatch = SP.nonzero(SAMPLE_IDS[:,SP.newaxis]==sample_ids[:,SP.newaxis].T)[1]
    assert (sample_ids[Imatch]==SAMPLE_IDS).all(), 'outch'
    #filter gene_ids, only keeping ENSEMBL genes
    Iok = SP.array([(re.match('ENS.*',gene_id) is not None) for gene_id in gene_ids])
    counts = counts[Iok,:][:,Imatch]
    #lib size of end. RNA
    __library_size_end = counts.sum(axis=0)
    gene_ids = gene_ids[Iok]
    annotation = annotation[Iok]

    fgenes=fout.create_group('genes')
    fgenes.create_dataset(name='counts',data=counts)
    fgenes.create_dataset(name='annotation',data=annotation)
    fgenes.create_dataset(name='ids',data=gene_ids)
    calc_norm(fgenes)

    #2. summarize exon expression
    sample_ids = fcount['exons/sample_ids'][:]
    counts    = fcount['exons/counts'][:]
    exon_ids  = fcount['exons/exon_ids'][:]
    annotation = fcount['exons/annotation'][:]
    #fix sample_ids
    Imatch = SP.nonzero(SAMPLE_IDS[:,SP.newaxis]==sample_ids[:,SP.newaxis].T)[1]
    assert (sample_ids[Imatch]==SAMPLE_IDS).all(), 'outch'
    counts = counts[:,Imatch]

    #filter exon_ids, only keeping ENSEMBL genes
    Iens = SP.array([(re.match('ENS.*',exon_id) is not None) for exon_id in exon_ids])
    counts_ens = counts[Iens,:]
    exon_ids_ens = exon_ids[Iens]
    annotation_ens = annotation[Iens]

    fexons=fout.create_group('exons')
    fexons.create_dataset(name='counts',data=counts_ens)
    fexons.create_dataset(name='annotation',data=annotation)
    fexons.create_dataset(name='ids',data=exon_ids)
    calc_norm(fexons)

    #parse out ERCCs
    Iercc = SP.array([(re.match('ERCC.*',exon_id) is not None) for exon_id in exon_ids])
    counts_ercc = counts[Iercc,:]
    exon_ids_ercc = exon_ids[Iercc]
    fercc = fout.create_group('ERCC')
    fercc.create_dataset(name='counts',data=counts_ercc)
    fercc.create_dataset(name='ids',data=exon_ids_ercc)
    
    #diagnostic
    Idiag = SP.array([(re.match('__.*',exon_id) is not None) for exon_id in exon_ids])
    counts_diag = counts[Idiag,:]
    __library_size = counts.sum(axis=0)
    #convert to relative fractions
    counts_diag = SP.concatenate((__library_size[SP.newaxis,:],__library_size_end[SP.newaxis,:],counts_diag),axis=0)
    fractions_diag = SP.array(counts_diag,dtype='float')/counts_diag[0]
    exon_ids_diag = SP.concatenate((['__library_size','__library_size_endogenous'],exon_ids[Idiag]))
    fqual = fout.create_group('qual')
    fqual.create_dataset(name='counts',data=counts_diag)
    fqual.create_dataset(name='fractions',data=fractions_diag)
    fqual.create_dataset(name='ids',data=exon_ids_diag)
 
    #variant assignment
    sample_ids = fvar['sample_ids'][:]
    geno_id    = fvar['geno_id'][:]
    error      = fvar['error'][:]
    geno_id_map = fvar['geno_id_map'][:]

    Imatch = SP.nonzero(SAMPLE_IDS[:,SP.newaxis]==sample_ids[:,SP.newaxis].T)[1]
    assert (sample_ids[Imatch]==SAMPLE_IDS).all(), 'outch'
    geno_id_map = geno_id_map[Imatch]
    error       = error[Imatch]
    fvariant = fout.create_group('variants')
    fvariant.create_dataset(name='geno_ids',data=geno_id)
    fvariant.create_dataset(name='error',data=error)
    fvariant.create_dataset(name='geno_id_map',data=geno_id_map)

    fout.create_dataset(name='sample_ids',data=SAMPLE_IDS)
    fout.create_dataset(name='qualityOK',data=qualityOK)

    Iread=(counts_diag[0]-counts_diag[4]>1E6)
    LU = SP.unique(geno_id_map[qualityOK&Iread])
    for line in LU:
        print '%s:%d' %(line,(geno_id_map[qualityOK&Iread]==line).sum())
