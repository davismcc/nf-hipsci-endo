"""Run cis QTL analyses for endodiff single-cell data

Usage:
  eqtl_cis_endodiff.py --chrom=CHROM --genot_file=FILE_GENOT --exprs_file=FILE_EXPRS --kpop_file=FILE_KPOP --output_file=FILE_OUT [-d] [-p] [-h | -v]
  eqtl_cis_endodiff.py --version 
  eqtl_cis_endodiff.py --help

Process raw expression data from an HDF5 file and save processed data to another
HDF5 file.

Options:
  -c --chrom=CHROM              chromosome to run
  -i --genot_file=FILE_GENOT    input genotype HDF5 file
  -e --exprs_file=FILE_EXPRS    input expression HDF5 file
  -k --kpop_file=FILE_KPOP      input Kpop (GRM) HDF5 file
  -o --output_file=FILE_OUT     output file
  -d --design                   flag to fit covariates in 'design' slot of expression HDF5 file
  -p --permute                  permute genotypes to generate permutation results for calibration
  --version                     show version and exit
  -h --help                     show this help message and exit

Example usage command:
python src/eqtl/eqtl_cis_endodiff.py --input_file=data_raw/diff_1/run_19776/bam/19776_2#160.bam --output_file=data_raw/diff_1/run_19776/variants/19776_2#160_id_donor.hdf5 

Written by Davis McCarthy, January 2017
"""

################################################################################

import sys
import re
sys.path.append('./../')
import limix.modules.qtl as QTL
import limix.stats.fdr as FDR
#from include.utils import getLambda
import scipy as sp
import numpy as np
import h5py
import pdb
import warnings
import os
from docopt import docopt
from datetime import datetime
import pandas as pd
sys.path.append('./../include')
#from normalization import gaussianize


################################################################################

def getGenotypes(geno_file, indiv_idx=None, return_info=False, dosage=False, 
                     min_maf=0.05, min_mac=2):
    """ return all genotypes (hard calls or dosage) on the chromosome, 
    as appropriate for trans-eQTL mapping 
    """
    fgeno = h5py.File(geno_file)
    if dosage:
        X = fgeno['genotype']['matrix_dosage'][:]
    else:
        X = fgeno['genotype']['matrix'][:]
    ## compute MAF and filter variants based on min_maf threshold
    print "        ...filtering on MAF"
    sample_names = fgeno['sample_info']['sampleID'][:]
    donor = [re.sub("_[0-9]*", "", re.sub("HP.*-", "", x)) for x in sample_names]
    dat = pd.DataFrame({'donor': donor})
    dup = sp.array(dat.duplicated())
    X_for_maf = X[~dup,:]
    num_nan = sp.isnan(X_for_maf).sum(axis = 0)
    nvars = X_for_maf.shape[0] - num_nan
    maf = np.nansum(X_for_maf, axis = 0) / 2 / nvars
    maf_ok = ((maf > min_maf) & (maf < 1 - min_maf))
    ## compute MAC for individuals used and filter based on min_mac threshold
    print "        ...filtering on MAC"
    if indiv_idx is not None:
        sample_uniq = np.unique(sample_names[indiv_idx])
        mac_idx = match(sample_names, sample_uniq)
        mac = np.nansum(X[mac_idx,:], axis = 0)
        mac_ok = ((mac >= min_mac) & (mac <= 2 * len(sample_uniq) - min_mac))
        ## ensure at least two samples have minor allele
        nhomref = np.nansum(X[mac_idx,:] == 0, axis = 0)
        nhomalt = np.nansum(X[mac_idx,:] == 2, axis = 0)
        mac_ok = (mac_ok & (nhomref < len(sample_uniq) - min_mac) & (nhomalt < len(sample_uniq) - min_mac))
    else:
        mac = np.nansum(X, axis = 0)
        mac_ok = ((mac >= min_mac) & (mac <= 2 * len(sample_names) - min_mac))
        nhomref = np.nansum(X == 0, axis = 0)
        nhomalt = np.nansum(X == 2, axis = 0)
        mac_ok = (mac_ok & (nhomref < len(sample_uniq) - min_mac) & (nhomalt < len(sample_uniq) - min_mac))
    ## filter variants
    X = X[:,(maf_ok & mac_ok)]
    ## select individual lines by index
    if indiv_idx is not None:
        X = X[indiv_idx,:]
    if return_info:
        info = {}
        info['maf'] = maf[(maf_ok & mac_ok)]
        for key in fgeno['snp_info']:
            if key not in ['sampleID']:
                info[key] = fgeno['snp_info'][key][(maf_ok & mac_ok)]
        fgeno.close()
        return X, info
    else:
        fgeno.close()
        return X


def getLambda(pv):
    """
    return lambda genomic control given the pvs
    """
    rv = sp.array([sp.median(sp.log10(pv),1)/sp.log10(0.5)])
    return rv


def match(a, b):
    """
    return indices of a that match the elements of b
    """
    mlist =  [ (a == x).nonzero()[0][0] if x in a else None for x in b ]
    mlist = [ x for x in mlist if x is not None ]
    m = np.array(mlist)
    return m


if __name__=='__main__':
    ## Main function that can run an analysis
    ## parse command-line arguments
    arguments = docopt(__doc__, version = "0.0.1")
    chrom = str(arguments['--chrom'])
    chrom = chrom.replace('chr', '')
    genot_file = arguments['--genot_file']
    exprs_file = arguments['--exprs_file']
    kpop_file = arguments['--kpop_file']
    out_file = arguments['--output_file']
    if arguments['--design']:
        fit_design = True
    else:
        fit_design = False
    if arguments['--permute']:
        permute = True
    else:
        permute = False
    print 'HDF5 Genotype Filename: ' + genot_file
        
    ## get sample IDs for genotype and Kpop 
    fgeno = h5py.File(genot_file)
    sampleID_geno = fgeno['sample_info']['sampleID'][:]

    # read expression data and gene info
    fexprs = h5py.File(exprs_file)
    sample_gene_id = fexprs['sampleID'][:]
    sample_gene_idx = np.array([x in fgeno['sample_info']['sampleID'][:] for x in sample_gene_id])
    Y = fexprs['exprs'][:]
    # read in design if appropriate
    if fit_design:
        design = fexprs['design'][:]
        design = design[sample_gene_idx,:]
    # get sample indices for genotypes and Kpop
    sample_geno_idx = match(sampleID_geno, sample_gene_id)
    # get indices for genes on this chromosome
    gene_info = {}
    for key in fexprs['gene_info']:
        gene_info[key] = fexprs['gene_info'][key][:]
    gene_info = pd.DataFrame(gene_info)
    gene_idx = np.array(gene_info['chromosome_name'] == chrom)
    # subset expression matrix and gene info to genes on this chromosome
    gene_info = gene_info.iloc[gene_idx]
    Y = Y[sample_gene_idx, :][:, gene_idx]
    genes = gene_info['hgnc_symbol'].values.astype(str)
    ## gaussianize expression values
    fout = h5py.File(out_file, 'a')
    fout.create_dataset('geneID', data=genes)
    fout.create_dataset('sampleID', data=sample_gene_id[sample_gene_idx])
    fout.create_dataset('exprs', data=Y)
    # Y = gaussianize(Y)
    # fout.create_dataset('exprs_gauss', data=Y)
    
    # get Kpop, genotypic covariance matrix between individuals
    fK = h5py.File(kpop_file)
    sample_K_idx = match(fK['sampleID'][:], fexprs['sampleID'][:])
    K = fK['Kpop'][:]
    K = K[sample_K_idx,:][:,sample_K_idx]
    K /= K.diagonal().mean()
    K += 1e-4 * sp.eye(K.shape[0])

    # Try to get genotypes
    print '    ... Importing genotype data'
    try:
        X, info = getGenotypes(genot_file, indiv_idx=sample_geno_idx, 
                                   return_info=True, dosage=False, 
                                   min_maf=0.05, min_mac=5)
    except:
        print 'Error: no SNPs found'
    info_df = pd.DataFrame(info)
    fout.create_dataset('genotypes', data=X)
    fout.create_dataset('gdid', data=info['gdid'])
    fout.close()

    # open output file for saving pandas tables
    fout = pd.HDFStore(out_file, mode='a', complevel=9, complib='blosc')
    ## set options so that tables can be filtered on file
    pd.set_option('io.hdf.default_format','table')

    # loops over genes
    # take only unique genes
    all_results = None
    for gene in genes:
        #pdb.set_trace()
        # Define gene group
        print '.. Analyzing gene %s' % gene
        # gene_group = fout.create_group(gene)
        gene_idx = np.array(gene==genes)
        Y_gene = Y[:, gene_idx]
        # Apply permutation if desired
        out_dict = {}
        ## select SNPs to test for this gene - within 500 MB of gene start
        gene_start = gene_info['start_position'][gene_idx].values[0]
        min_pos = np.array([0, gene_start - 250000]).max()
        max_pos = gene_start + 250000
        snp_pos = info_df['pos'].values
        snp_idx = np.logical_and((snp_pos > min_pos), (snp_pos < max_pos))
        if not snp_idx.any():
            continue
        X_gene = X[:, snp_idx]
        if permute:
            X_gene = np.random.permutation(X_gene)
        info_gene = info_df.iloc[snp_idx]
        for item in info_gene.columns:
            out_dict[item] = info_gene[item].values
        assoc_gene = gene.repeat(snp_idx.sum())
        out_dict['assoc_gene'] = assoc_gene
        # Run the LMM analysis
        print "   .. single trait analysis"
        if fit_design:
            lmm = QTL.test_lmm(X_gene, Y_gene, K=K, covs = design) 
        else:
            lmm = QTL.test_lmm(X_gene, Y_gene, K=K) 
        pv = lmm.getPv()
        pv[0][np.isnan(pv[0])] = 1.0 # set any NaN p-values to 1
        out_dict['pv'] = pv[0]
        out_dict['qv'] = FDR.qvalues(pv)[0]
        out_dict['beta'] = lmm.getBetaSNP()[0]
        lambda_val = getLambda(pv)
        lambda_val = lambda_val.repeat(len(out_dict['pv']))
        out_dict['lambda_val'] = lambda_val
        out_df = pd.DataFrame(out_dict, index=out_dict['gdid'])
        ## convert full stops in gene name to underscore
        gene = gene.replace(".", "_")
        ## append results for chunk to gene's results df to HDF5 file
        print "    ....appending results..."
        fout.put(gene, out_df, table=True, data_columns=True)
        if all_results is None:
            all_results = out_df
        else:
            all_results = all_results.append(out_df, ignore_index=True)
        print "    ....appending done."
            
    ## put all results in one big dataframe
    print "Adding results in one big dataframe to output"
    chrname = 'chr' + chrom
    fout.put(chrname, all_results, table=True, data_columns=True)
    ## add gene info to output
    fout.put('gene_info', gene_info, table=True, data_columns=True)
    fout.close()
    ## add metadata to attributes
    fout = h5py.File(out_file)
    fout.attrs['time_stamp'] = datetime.now().isoformat()
    fout.attrs['author'] = 'Davis McCarthy'
    fout.attrs['institution'] = 'EMBL-EBI'
    fout.close()
    
    print "cis eQTL mapping completed"
    print "Processed %s genes for chromosome %s" %(len(genes), chrom)
    print "Results saved to %s" %out_file

# End
    

