"""Compute Kpop (kinship) matrices for HipSci cell lines for scRNA-seq eQTL mapping

Usage: 
  compute_Kpop.py [-h] [--debug]
  compute_Kpop.py --version

Combine imputed genotype data across chromosomes from HDF5 files, compute the
kinship matrices as XX' on normalised genotypes and save processed data to an 
HDF5 file.

Options:
  -h --help               show this help message and exit
  --version               show version and exit
  --debug                 run in debugging mode

For defaults for the input and output data, see the configuration settting in 
../CFG/settings.py

Written by Davis McCarthy, November 2016
"""

################################################################################

## import packages
import scipy as sp
import h5py
import os
from docopt import docopt
from datetime import datetime
import sys
import numpy as np
sys.path.append('./../')
sys.path.append('./')
sys.path.append('./src')
from CFG.settings import CFG
from include.data_rnaseq import QtlData
import pdb

################################################################################

def getGenotypes(geno_file, indiv_idx=None, dosage=False, return_info=False):
    """ return all genotypes (hard calls or dosage) on the chromosome, 
    as appropriate for trans-eQTL mapping 
    """
    fgeno = h5py.File(geno_file)
    if dosage:
        X = fgeno['genotype']['matrix_dosage'][:]
    else:
        X = fgeno['genotype']['matrix'][:]
    if indiv_idx is not None:
        X = X[indiv_idx,:]
    if return_info:
        info = {}
        for key in fgeno['snp_info']:
            if key not in ['sampleID']:
                info[key] = fgeno['snp_info'][key][:]
        fgeno.close()
        return X, info
    else:
        fgeno.close()
        return X


if __name__ == '__main__':
    ## parse command line arguments and define parameters
    arguments = docopt(__doc__, version = "0.0.1")

    ## get sample IDs
    fgeno = h5py.File(CFG['eqtl']['rnaseq_data']['1'])
    sampleID = fgeno['sample_info']['sampleID'][:]
    
    ## define output file
    outFile = CFG['eqtl']['Kpop_file']
    fOut = h5py.File(outFile, 'a')
    if 'sampleID' not in fOut:
        fOut.create_dataset('sampleID', data=sampleID)
    else:
        del fOut['sampleID']
        fOut.create_dataset('sampleID', data=sampleID)

    ## combine genotype data across all chromosomes
    chrom_list = [str(x) for x in range(1, 23)]
    #chrom_list.append('X')
    Kchrom_list = {}
    nvars_list = {}
    for dsg in ['no']:
        if dsg == 'yes':
            dosage = True
        else: 
            dosage = False
        print "Using dosages? ...%s" %dsg
        genoData_by_chrom = {}
        genoInfo_by_chrom = {}
        X = None
        for chrom in chrom_list:
            print "...chr%s" % chrom
            # read data
            geno_file = CFG['eqtl']['rnaseq_data'][chrom]
            # Try to get genotypes
            print '    ... Importing genotype data'
            try:
                Xchrom, info = getGenotypes(geno_file, return_info=True, dosage=dosage)
                genoData_by_chrom[chrom] = Xchrom
                genoInfo_by_chrom[chrom] = info
            except:
                print 'Error: no SNPs found'
            ## genotypes are arranged samples x variants
            ## normalise genotypes
            zero_sd = (Xchrom.std(0) == 0)
            Xchrom = Xchrom[:, ~zero_sd]
            Xchrom -= Xchrom.mean(0)
            Xchrom /= Xchrom.std(0)
            ## compute Kpop
            nvars = X.shape[1]
            Kchrom_list[chrom] = np.dot(X, X.transpose()) / nvars
            
        for chrom in chrom_list:
            if chrom == '1':
                Kpop_raw = Kchrom_list[chrom]
            else: 
                Kpop_raw += Kchrom_list[chrom]

        Kpop_raw /= len(chrom_list)
        Kpop = {}
        Kpop['Kpop_raw'] = Kpop_raw
        ## normalise Kpop values
        Kpop_norm = Kpop_raw
        Kpop_norm /= Kpop_norm.diagonal().mean()
        # numerical fix to matrix:
        Kpop_norm += 1e-3 * sp.eye(Kpop_norm.shape[0])
        Kpop['Kpop_norm'] = Kpop_norm

        ### Export Kpop into chromosome-wise datasets
        print "..exporting Kpop"
        ## add/replace Kpop data
        for key1 in ['Kpop_norm', 'Kpop_raw']:
            if key1 not in fOut:
                fOut.create_group(key1)
            if dosage:
                print "Kpop computed with genotype dosages"
                keys2 = ['Kpop_std_dosage', 'Kpop_dosage']
            else:
                print "Kpop computed with standard genotype"
                keys2 = ['Kpop_std', 'Kpop']
            for key2 in keys2:
                print "Key1: %s - Key2: %s" % (key1, key2)
                if key2 in fOut[key1]:
                    del fOut[key1][key2]
                fOut[key1].create_dataset(key2, data=Kpop[key1])
    
    ## add metadata about run
    fOut.attrs['time_stamp'] = datetime.now().isoformat()
    fOut.attrs['author'] = 'Davis McCarthy'
    fOut.attrs['institution'] = 'EMBL-EBI'
    ## close connection to HDF5 file
    fOut.close()
    print "Done!"



    
