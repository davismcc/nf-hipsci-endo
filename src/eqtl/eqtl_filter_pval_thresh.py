"""Filter aggregated chromosome-level trans-eQTL results and save to HDF5 file

Usage: 
  eqtl_filter_pval_thresh.py <output_file> [-hc CHROM] [-p THRESH]

Options:
  -h --help           show this help message and exit
  --version           show version and exit
  -c --chrom=CHROM    chromosome name (1-22, X)
  -p --pval=THRESH    p-value threshold: associations with p-value below this 
                      threshold will be output (default = 1e-06)

Written by Davis McCarthy, February 2016
"""

################################################################################

import sys
sys.path.append('./../')
import eqtl_analysis_tools as eat
from CFG.settings import CFG
import os
from docopt import docopt
import pandas as pd
from datetime import datetime
import h5py

################################################################################

if __name__ == '__main__':
    ## parse command line arguments and define parameters
    arguments = docopt(__doc__, version = "0.0.1")
    ## sort out chromosome
    if arguments['--chrom'] is None:
        chrom = '22'
    else:
        chrom = str(arguments['--chrom'])
    ## sort out p-value threshold
    if arguments['--pval'] is None:
        pval_thresh = 1e-06
    else:
        pval_thresh = float(arguments['--pval'])
    ## define input file
    input_dir = CFG['out_rnaseq_aggregated']
    input_file = 'aggregated_results_chr' + chrom + '.hdf5'
    input_file = os.path.join(input_dir, input_file)
    ## sort out output file
    if arguments['<output_file>'] is None:
        raise ValueError('Please provide an output file name.')
    else:
        output_file = arguments['<output_file>']
    ## open hdf5 input file
    ## read in gene info
    gene_info = pd.read_hdf(input_file, 'gene_info')
    gene_info = eat.add_numeric_gene_pos(gene_info)
    ## add plurinet gene status to gene_info
    plurinet = pd.read_csv('../../resources/plurinet_geneID_list.txt', 
                               header=None)
    gene_info['plurinet'] = gene_info.index.isin(plurinet[0])
    ## read in snp info
    snp_info = pd.read_hdf(input_file, 'snp_info')
    snp_info = eat.add_numeric_var_pos(snp_info)
    ## get results for all genes
    where_arg = 'pv<' + str(pval_thresh)
    print "Reading in filtered results...this can take some time..."
    ## get keys
    fIN = h5py.File(input_file)
    chr_keys = fIN.keys()
    gene_keys = {}
    for chrom in chr_keys:
        gene_keys[chrom] = fIN[chrom].keys()
    fIN.close()
    fIN = pd.HDFStore(input_file)
    res_all_dict = {}
    for chr_key in chr_keys:
        for gene_key in gene_keys[chr_key]:
            print ".",
            qry = os.path.join(chr_key, gene_key)
            try: 
                gene_results = pd.read_hdf(fIN, qry, where=where_arg)
            except KeyError:
                print "KeyError: group %s/p not found in HDF5 input file" %qry
            if not gene_results.empty:
                res_all_dict[gene_key.replace("_", ".")] = gene_results
    #res_all_dict = eat.get_all_results_one_chrom(fIN, gene_info, where=where_arg)
    res_all_df = eat.flatten_results_dict_to_df(
        res_all_dict, gene_info, snp_info)
    ## output results
    print "Saving results for %d associations passing the threshold" %(
        res_all_df.shape[0])
    key = 'chr' + chrom
    res_all_df.to_hdf(output_file, key, format='table', complevel=9, 
                          complib='blosc', mode='w', data_columns=True)
    ## add metadata about run
    print "Adding metadata about run."
    hdf_out = h5py.File(output_file, 'a')
    hdf_out.attrs['time_stamp'] = datetime.now().isoformat()
    hdf_out.attrs['author'] = 'Davis McCarthy'
    hdf_out.attrs['institution'] = 'EMBL-EBI'
    ## close HDF5 output file
    hdf_out.close()
    print "All done!"
    
# Fin.

