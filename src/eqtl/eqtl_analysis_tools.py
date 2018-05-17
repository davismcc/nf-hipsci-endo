"""A collection of tools for analysing aggregated results of trans-eQTL mapping 
in HipSci RNAseq data

Usage: import these functions into a Python or IPython session.

Written by Davis McCarthy, February 2016
"""

################################################################################

import sys
import h5py
import pdb
import copy
import warnings
import os
from datetime import datetime
import pandas as pd
import scipy as sp
import numpy as np
## local code
sys.path.append('./../')
from include.data_rnaseq import QtlData
from CFG.settings import CFG
## plotting
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from matplotlib import gridspec
# plotting and visualization utilties
import limix.utils.plot as lmxplt
# genotype summary stats
import limix.stats.geno_summary as lmxgeno


################################################################################
## ------------------ functions for wranging aggegated results  ----------------

def get_gene_group(gene_info, gene):
    """Convenience function to get a correct gene group to query HDF5 table
    from a Pandas data-frame containing the gene information.
    """
    chrom = gene_info.ix[gene.replace('_', '.')]['chr']
    gene_group = os.path.join(chrom, gene)
    ## make sure gene ID will work with the HDF5 table keys
    gene_group = gene_group.replace('.', '_')
    return gene_group

def get_gene_results_one_chrom(hdfstore, gene_info, gene_list, where='pv<5e-08'):
    """Access trans-eQTL results from a Pandas HDFStore for multiple genes
    """
    results = {}
    for gene in gene_list:
        qry = get_gene_group(gene_info, gene)
        gene_results = pd.read_hdf(hdfstore, qry, where=where)
        if not gene_results.empty:
            results[gene] = gene_results
    return results

def get_all_results_one_chrom(hdfstore, gene_info, where='pv<5e-08'):
    """Get results across all genes that satisfy a threshold for one chromosome
    """
    results = {}
    for chr_key in hdfstore.keys():
        for gene_key in hdfstore[chr_key]:
            print ".",
            qry = os.path.join(chr_key, gene_key)
            try: 
                gene_results = pd.read_hdf(hdfstore, qry, where=where)
            except KeyError:
                print "KeyError: group %s/p not found in HDF5 input file" %qry
            if not gene_results.empty:
                results[gene_key.replace("_", ".")] = gene_results
    return results

def add_numeric_gene_pos(gene_info):
    """
    Add numeric gene (start) genomic position to a gene_info dataframe
    """
    gene_chr_numeric = gene_info['chr']
    gene_chr_numeric = ['23' if x == 'X' else x for x in gene_chr_numeric]
    gene_chr_numeric = ['24' if x == 'Y' else x for x in gene_chr_numeric]
    gene_start_vec = gene_info['start']
    gene_start_vec = [str(x).zfill(10) for x in gene_start_vec]
    gene_pos_numeric = [x + '.' + y for x, y in zip(gene_chr_numeric, gene_start_vec)]
    gene_pos_numeric = np.array([float(x) for x in gene_pos_numeric])
    gene_info['genome_pos_numeric'] = gene_pos_numeric
    return gene_info

def add_numeric_var_pos(snp_info):
    """
    Add numeric variant genomic position to a snp_info dataframe
    """
    var_chr_vec = [s.split('_')[1] for s in snp_info['gdid']]
    var_chr_vec = ['23' if x == 'X' else x for x in var_chr_vec]
    var_chr_vec = ['24' if x == 'Y' else x for x in var_chr_vec]
    var_pos_vec = [s.split('_')[2].zfill(10) for s in snp_info['gdid']]
    var_pos_numeric = [x + '.' + y for x, y in zip(var_chr_vec, var_pos_vec)]
    var_pos_numeric = np.array([float(x) for x in var_pos_numeric])
    snp_info['genome_pos_numeric'] = var_pos_numeric
    return snp_info

def order_genes_by_pos(gene_info):
    """
    Return indices to order genes by genomic position
    """
    gene_chr_numeric = gene_info['chr']
    gene_chr_numeric = ['23' if x == 'X' else x for x in gene_chr_numeric]
    gene_chr_numeric = ['24' if x == 'Y' else x for x in gene_chr_numeric]
    gene_start_vec = gene_info['start']
    gene_start_vec = [str(x).zfill(10) for x in gene_start_vec]
    gene_pos_numeric = [x + '.' + y for x, y in zip(gene_chr_numeric, gene_start_vec)]
    gene_pos_numeric = np.array([float(x) for x in gene_pos_numeric])
    genes_order = gene_pos_numeric.ravel().argsort()
    return(genes_order)

def order_vars_by_pos(snp_info):
    """
    Return indices to order variants by genomic position given a list of gdid
    identifiers for variants
    """
    var_chr_vec = [s.split('_')[1] for s in snp_info['gdid']]
    var_pos_vec = [s.split('_')[2].zfill(10) for s in snp_info['gdid']]
    var_pos_numeric = [x + '.' + y for x, y in zip(var_chr_vec, var_pos_vec)]
    var_pos_numeric = np.array([float(x) for x in var_pos_numeric])
    var_order = var_pos_numeric.ravel().argsort()
    return(var_order)


def flatten_results_dict_to_df(results, gene_info, snp_info):
    """From a results dictionary where elements are gene IDs and contain results
    for a set of variants (SNPs)
    """
    new_df = pd.DataFrame()
    for gene in results.keys():
        this_df = results[gene]
        nhits = this_df.shape[0]
        ## add snp info
        snp_idx = this_df.index
        this_snp_info = snp_info.ix[snp_idx]
        this_df['snp_chrom'] = this_snp_info['chrom']
        this_df['snp_pos'] = this_snp_info['pos']
        if 'genome_pos_numeric' in this_snp_info:
            this_df['snp_pos_numeric'] = this_snp_info['genome_pos_numeric']
        this_df['info'] = this_snp_info['info']
        this_df['gdid'] = this_snp_info['gdid']
        this_df['maf'] = this_snp_info['maf']
        ## add gene info
        this_gene_info = gene_info.ix[gene]
        this_df['gene_chrom'] = [this_gene_info['chr'] for i in range(nhits)]
        this_df['gene_start'] = [this_gene_info['start'] for i in range(nhits)]
        this_df['gene_end'] = [this_gene_info['end'] for i in range(nhits)]
        this_df['gene_length'] = [this_gene_info['length'] for i in range(nhits)]
        this_df['gene_name'] = [this_gene_info['name'] for i in range(nhits)]
        this_df['geneID'] = [this_gene_info['geneID'] for i in range(nhits)]
        this_df['gene_strand'] = [this_gene_info['strand'] for i in range(nhits)]
        this_df['gene_biotype'] = [this_gene_info['type'] for i in range(nhits)]
        if 'genome_pos_numeric' in this_gene_info:
            gene_pos_numeric = this_gene_info['genome_pos_numeric'].repeat(nhits)
            this_df['gene_pos_numeric'] = gene_pos_numeric
        ## use an index that combines geneID and gdid
        #pdb.set_trace()
        new_idx = [x+'_'+y for x, y in zip(this_df['geneID'], this_df['gdid'])]
        this_df.index = new_idx
        ## append this df to the new_df object
        new_df = new_df.append(this_df)
    ## return new flat dataframe with all results
    return new_df
        
def neg_log10_p(x):
    """
    Return -log10(p-values) given a list of p-values
    """
    return -np.log10(x)

################################################################################
## -------------------------- plotting functions -------------------------------

def dotplot(df, filename = '../../figures/dotplot.png', 
            title='Significant associations', show=False):
    """Make a dotplot showing associations by plotting gene position against
    variant position.
    """
    print datetime.now().isoformat()
    print "Making dotplot for given dataframe ..."
    plt.style.use('grayscale')
    #sns.set_style("white")
    sns.set_context("talk")
    fig = plt.figure(figsize=(15, 15)) 
    max_snp_chrom_pos = range(1, 24)
    ## Setup the axes
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 7]) 
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex = ax0)
    ## make the plots
    ax0.hist(df['snp_pos_numeric'], bins=1000, color='b')
    ax0.grid(b=True, which='both', color='0.5')
    ax1.scatter(df['snp_pos_numeric'], df['gene_pos_numeric'], alpha=0.1)
    ## define titles, ticks, labels etc
    ### hist plot
    ax0.set_title(title)
    ax0.set_ylabel('count')
    ax0.set_xticks(max_snp_chrom_pos)
    ax0.set_xticklabels(np.arange(22) + 1)
    ### scatter plot
    ax1.set_xlabel('variant position (by chromosome)')
    ax1.set_xticks(max_snp_chrom_pos)
    ax1.set_xticklabels(np.arange(22) + 1)
    ax1.set_yticks(max_gene_chrom_pos[1:])
    ax1.set_yticklabels(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                         '11', '12', '13', '14', '15', '16', '17', '18', '19',
                         '20', '21', '22', 'X', 'Y'])
    ax1.set_ylabel('gene position (by chromosome)')
    ax1.grid(b=True, which='both', color='0.5')
    plt.xlim([0, df['snp_pos_numeric'].max(0)])
    plt.ylim([0, df['gene_pos_numeric'].max(0)])
    plt.subplots_adjust(wspace=0, hspace=0.1)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    if show:
        plt.show()

def plot_manhattan(df, chromBounds = range(1, 23), title='Manhattan plot', 
                       gene=None, filename='../../figures/manhattan_plot.png'):
    """Produce a manhattan plot for a set of associations
    """
    sns.set_style("white")
    sns.set_context("talk")
    plt.figure(figsize=[15,5])
    lmxplt.plot_manhattan(df['snp_pos_numeric'], df['pv'], chromBounds)
    plt.title(title)
    if gene is not None:
        plt.plot(df[df['gene_name']==gene]['gene_pos_numeric'][0], 1, '^r')
    plt.savefig(filename, dpi=300, bbox_inches='tight')

################################################################################


