"""Identify HipSci donor line for a cell RNA-seq sample

Usage:
  eqtl_cis_endodiff_aggregate.py --run=RUN_TYPE [-h | -v]
  eqtl_cis_endodiff_aggregate.py --version 
  eqtl_cis_endodiff_aggregate.py --help

Process raw expression data from an HDF5 file and save processed data to another
HDF5 file.

Options:
  -r --run=RUN_TYPE             string name for the type of run/analysis
  --version                     show version and exit
  -h --help                     show this help message and exit

Example usage command:
python src/eqtl/eqtl_cis_endodiff_aggregate.py --run='eqtl_results_design_day3'

Written by Davis McCarthy, February 2017
"""

import numpy as np
import pandas as pd
import feather
from docopt import docopt

if __name__=='__main__':
    arguments = docopt(__doc__, version = "0.0.1")
    run_type = str(arguments['--run'])
    output_h5='eqtl/endodiff_cis.%s.aggregated.h5' %run_type
    output_df_feather='eqtl/endodiff_cis.%s.aggregated.feather' %run_type
    output_info_feather='eqtl/endodiff_cis.%s.agg_gene_info.feather' %run_type
    ## read in results
    df = pd.DataFrame()
    gene_info = pd.DataFrame()
    for chrom in np.arange(1, 23).astype(str):
        res_file = 'eqtl/endodiff_cis_%s.chr%s.h5' %(run_type, chrom)
        chrom_name = 'chr%s' %chrom
        chrom_res = pd.read_hdf(res_file, chrom_name)
        chrom_gene_info = pd.read_hdf(res_file, 'gene_info')
        # gb = chrom_res.groupby(['assoc_gene'])
        # min_qv = gb['qv'].min()
        # keep_gene = (min_qv < 0.1)
        # sig_genes = min_qv.index[np.array(keep_gene)].values.astype(str)
        # chrom_res_filt = chrom_res.loc[chrom_res['assoc_gene'].isin(sig_genes)]
        # chrom_gene_info_filt = chrom_gene_info.loc[chrom_gene_info['hgnc_symbol'].isin(sig_genes)]
        df = df.append(chrom_res)
        gene_info = gene_info.append(chrom_gene_info)
        gene_info.index = gene_info['hgnc_symbol']
    ## save results to file
    df.to_hdf(key='cis_results', path_or_buf=output_h5)
    gene_info.to_hdf(key='gene_info', path_or_buf=output_h5)
    feather.write_dataframe(df, output_df_feather)
    feather.write_dataframe(gene_info, output_info_feather)
