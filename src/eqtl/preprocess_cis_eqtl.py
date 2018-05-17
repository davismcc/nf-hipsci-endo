"""Preprocessing of HipSci RNAseq data (REL-2016-01) for trans-eQTL mapping

Usage: 
  preprocess.py [-hp NAME] [--max_na_pheno=MAX_NA] [--input_data_dir=DIR] 
                [--output_data_dir=DIR] [--pheno_file=FILE] [--Kpop_file=FILE]
                [--output_file=FILE]
  preprocess.py --version

Process raw expression data from an HDF5 file and save processed data to another
HDF5 file.

Options:
  -h --help               show this help message and exit
  --version               show version and exit
  -p DIR --plot=DIR       plot results and save to filename
  --max_na_pheno=MAX_NA   maximum proportion of missing (NA, i.e zero here) 
                          values acceptable to retain a gene for analysis (default=0.1)
  --input_data_dir=DIR    directory in which to find input data 
                          (default: ../../../hipsci_eqtl_data/eQTL_rnaseq/)
  --output_data_dir=DIR   directory in which to save output data
  --pheno_file=FILE       filename for the input phenotype (expression) data
  --Kpop_file=FILE        filename for the Kpop (population kinship) data
  --output_file=FILE      filename for the output processed  phenotype data

For defaults for the input data, see the configuration settting in 
../CFG/settings.py

Written by Davis McCarthy, February 2016
Based on earlier code written by Helena Kilpinen, October/November 2015
"""

################################################################################

import sys
sys.path.append('./../')
from CFG.settings import CFG
import os
import re
import h5py
import scipy as SP
import scipy.linalg as LA
import numpy as np
import numpy.linalg as linalg
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend.
mpl.use('Agg')
import seaborn as sns
import pylab as PL
sys.path.append('./../include')
from normalization import gaussianize
from docopt import docopt
from datetime import date
from datetime import datetime
import pandas as pd
import pdb

################################################################################

def PC_varExplained(Y, standardize=True):
	"""Run PCA and calculate the cumulative fraction of variance

	Args:
	    Y (dbl): phenotype values
	    standardize (logical): if True, phenotypes are standardized

	Returns:
        var (dbl): cumulative distribution of variance explained
	"""
	# figuring out the number of latent factors
	if standardize:
		Y -= Y.mean(0)
		Y /= Y.std(0)
	covY = SP.cov(Y)
	S,U = linalg.eigh(covY + 1e-6 * SP.eye(covY.shape[0]))
	S = S[::-1]
	rv = np.array([S[0:i].sum() for i in range(1, S.shape[0])])
	rv /= S.sum()
	return rv

################################################################################

def export_results(outFile, sampleInfo, geneAnn, Kpop, Y_counts,
                   Y_exprs, Y_tpm, Y_residuals_raw,
                   Y_residuals_gaussianized, Y_residuals_standardized):
    """Export results to an HDF5 file

    Args:
        outFile (str): filename for the output HDF5 file
        sampleInfo (dict): dict containing sample information
        geneAnn (dict): dict containing gene annotation information
        Kpop (dict): dict of K (kinship) matrices for the population
        Y_counts (dbl): raw gene-level counts 
        Y_exprs (dbl): log2(tpm + 1) gene-level expression data
        Y_residuals_gaussianized (dbl): Gaussianized gene-level expression 
                           residuals after regressing out effects of gender and 
                           growing conditions
                           
    Returns:
        null

    This function is designed only to save gene-level expression data as exon 
    and probe level expression data is not as relevant for trans-eQTL mapping.
    """
    ## add/replace metadata information
    ### Export metadata (same for all):
    print "..exporting metadata:"
    fOut = h5py.File(outFile, 'w')       
    # pdb.set_trace()
    # fOut = pd.HDFStore(outFile, mode='w', complevel=9, complib='blosc')
    ## set options so that tables can be filtered on file
    # pd.set_option('io.hdf.default_format','table')
    if "/meta" in fOut:
        del fOut['meta']
    if "/sample_info" in fOut:
        del fOut['sample_info']
    gSampleInfo = fOut.create_group('sample_info')
    sample_info = sampleInfo.to_dict('list')
    for key in sample_info.keys():
        gSampleInfo.create_dataset(key, data=np.array(sample_info[key]))
    # sampleInfo.to_hdf(fOut, "sample_info", format = 'table', complevel = 9, 
    #                       complib = 'blosc')
    # fOut.put('sample_info', sampleInfo, table=True, data_columns=True)
    ## add/replace gene metadata information
    if "/gene_info" in fOut:
        del fOut['gene_info']
    gGeneInfo = fOut.create_group('gene_info')
    gene_info = geneAnn.to_dict('list')
    for key in gene_info.keys():
        gGeneInfo.create_dataset(key, data=np.array(gene_info[key]))
    # fOut.put('gene_info', geneAnn, table=True, data_columns=True)
    # geneAnn.to_hdf(fOut, "gene_info", format = 'table', complevel = 9, 
    #                    complib = 'blosc')
    # fOut.close()
    ## add gene expression data
    # fOut = h5py.File(outFile, 'a')
    ## add/replace gene expression data
    if "/phenotype_gene" in fOut:
        del fOut['phenotype_gene']
    gGENE = fOut.create_group('phenotype_gene')
    gGENE = fOut['phenotype_gene']
    gGENE.create_dataset('counts', data = Y_counts)
    gGENE.create_dataset('tpm', data = Y_tpm)
    gGENE.create_dataset('exprs', data = Y_exprs)
    gGENE.create_dataset('residuals_raw', data = Y_residuals_raw)
    gGENE.create_dataset('residuals_gaussianized', data = Y_residuals_gaussianized)
    gGENE.create_dataset('residuals_standardized', data = Y_residuals_standardized)
    ## add/replace Kpop data 
    if "/Kpop" in fOut:
        del fOut['Kpop']
    gKPOP = fOut.create_group('Kpop')
    for key1 in Kpop.keys():
        gKPOP.create_group(key1)
        for key2 in Kpop[key1].keys():
            gKPOP[key1].create_dataset(key2, data=Kpop[key1][key2])
    ## add metadata about run
    fOut.attrs['time_stamp'] = datetime.now().isoformat()
    fOut.attrs['author'] = 'Davis McCarthy'
    fOut.attrs['institution'] = 'EMBL-EBI'
    ## close connection to HDF5 file
    fOut.close()

################################################################################

def read_pheno_data(phenoFile):
    """Read in gene-based phenotype data from an HDF5 file
    
    Args:
        phenoFile (str): filename for the relevant HDF5 file

    Returns:
        dictionary with gene-level phenotype data
    """
    fpheno = h5py.File(phenoFile,'r')
    pheno_gene = {}
    pheno_gene['exprs'] = fpheno['gene_level/exprs'][:]
    pheno_gene['counts'] = fpheno['gene_level/counts'][:]
    pheno_gene['tpm'] = fpheno['gene_level/tpm'][:]
    pheno_gene['is_exprs'] = fpheno['gene_level/is_exprs'][:]
    fpheno.close()
    return pheno_gene

################################################################################

def read_Kpop_data(KpopFile):
    """Read in gene-based phenotype data from an HDF5 file
    
    Args:
        phenoFile (str): filename for the relevant HDF5 file

    Returns:
        dictionary with gene-level phenotype data
    """
    fKpop = h5py.File(KpopFile,'r')
    Kpop = {}
    for key1 in fKpop:
        Kpop[key1] = {}
        for key2 in fKpop[key1]:
            Kpop[key1][key2] = fKpop[key1][key2][:]
    fKpop.close()
    return Kpop

################################################################################

def read_sample_meta(phenoFile):
    """Read in gene-based metadata from an HDF5 file
    
    Args:
        phenoFile (str): filename for the relevant HDF5 file
    
    Returns:
        dictionary with sample metadata
    """
    sampleInfo = pd.read_hdf(phenoFile, 'gene_level/sample_metadata')
    sampleInfo.columns = [re.sub("\.", "_", x) for x in sampleInfo.columns]
    # fpheno = h5py.File(phenoFile,'r')
    # Sample meta (same for all):
    # sampleInfo = {}
    # sampleInfo['sampleID'] = fpheno['phenotype_gene']['colNames'][:]
    # sampleInfo['donor'] = fpheno['meta']['derived_from'][:]
    # sampleInfo['celltype'] = fpheno['meta']['derived_from_cell_type'][:]
    # sampleInfo['disease'] = fpheno['meta']['disease'][:]
    # sampleInfo['gender'] = fpheno['meta']['gender'][:]
    # sampleInfo['passage'] = fpheno['meta']['passage'][:]
    # sampleInfo['growing_conditions'] = fpheno['meta']['growing_conditions'][:]
    # sampleInfo['reprogramming'] = fpheno['meta']['reprogramming'][:]
    # sampleInfo['age'] = fpheno['meta']['age'][:]
    # sampleInfo['ethnicity'] = fpheno['meta']['ethnicity'][:]
    # sampleInfo['pluritest_novelty'] = fpheno['meta']['pluritest_novelty'][:]
    # sampleInfo['pluritest_raw'] = fpheno['meta']['pluritest_raw'][:]
    # sampleInfo['rnaseq_sendai_reads'] = fpheno['meta']['rnaseq_sendai_reads'][:]
    # sampleInfo['assaytime_rnaseq'] = fpheno['meta']['assaytime_rnaseq'][:]
    # sampleInfo['assayuser_rnaseq'] = fpheno['meta']['assayuser_rnaseq'][:]
    # fpheno.close()
    return sampleInfo

################################################################################

def read_feature_annos(phenoFile):
    """Read in gene-based metadata from an HDF5 file
    
    Args:
        phenoFile (str): filename for the relevant HDF5 file

    Returns:
        dictionary with feature annotations
    """
    geneAnn = pd.read_hdf(phenoFile, 'gene_level/feature_metadata')
    geneAnn.columns = [re.sub("\.", "_", x) for x in geneAnn.columns]
    # fpheno = h5py.File(phenoFile,'r')
    # # Feature annotations:
    # geneAnn = {}
    # geneAnn['geneID'] = fpheno['phenotype_gene']['rowNames'][:]
    # geneAnn['chr'] = fpheno['phenotype_gene']['gencode_v19']['chr'][:]
    # geneAnn['start'] = fpheno['phenotype_gene']['gencode_v19']['start'][:]
    # geneAnn['end'] = fpheno['phenotype_gene']['gencode_v19']['end'][:]
    # geneAnn['name'] = fpheno['phenotype_gene']['gencode_v19']['gene_name'][:]
    # geneAnn['type'] = fpheno['phenotype_gene']['gencode_v19']['gene_type'][:]
    # geneAnn['length'] = fpheno['phenotype_gene']['gencode_v19']['target_length'][:]
    # fpheno.close()
    return geneAnn

################################################################################

def regressOut(Y, X):
    """Regresses out X from Y

    Args: 
        Y (dbl): vector or matrix of dependent variable
        X (dbl): vector or matrix of independent variables to regress out

    Returns:
        vector or matrix of residuals from linear regression

    Assumes X is a p x n matrix where p is the number of variables and n the 
    number of samples.

    """
    Xd = LA.pinv(X)
    Y_out = Y - X.dot(Xd.dot(Y))
    return Y_out


################################################################################

if __name__ == '__main__':
    ## parse command line arguments and define parameters
    arguments = docopt(__doc__, version = "0.0.1")
    ## sort out max number of missing values allowed
    if arguments['--max_na_pheno']:
        max_na_pheno = float(arguments['--max_na_pheno'])
    else:
        max_na_pheno = 0.1
    print("Maximum acceptable number of missing values is %d" %max_na_pheno)
    ## sort out plotting
    plotDir = arguments['--plot']
    if plotDir:
        plot = True
        print("Plots will be saved in " + plotDir)
    else:
        plot = False
    ## sort out input directory
    if arguments['--input_data_dir'] is None:
        input_dir = CFG['REL-2016-01']['base']
        # input_dir = CFG['REL-2016-01']['input_dir']
        #input_dir = "../../data/eQTL_rnaseq/"
    else:
        input_dir = arguments['--input_data_dir']
    ## sort out output directory
    if arguments['--output_data_dir'] is None:
        output_dir = CFG['REL-2016-01']['base']
    else:
        output_dir = arguments['--output_data_dir']
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        print "Creating output directory: " + output_dir
    ## sort out input phenotype file
    if arguments['--pheno_file'] is None:
        pheno_file = CFG['REL-2016-01']['raw_pheno_file']
        #pheno_file = 'REL-2015-05.RNAseq.STAR.quantifications.HTSeq.meta_20151030.h5'
    else:
        pheno_file = arguments['--pheno_file']
    ## sort out Kpop file
    if arguments['--Kpop_file'] is None:
        Kpop_file = CFG['REL-2016-01']['Kpop_file']
        #Kpop_file = '/nfs/research/stegle/projects/hipsci/data/eQTL/Nov-15/eQTL_rnaseq/data/REL-2014-11.Kpop_maf01.ipsc.hdf5'
    else:
        Kpop_file = arguments['--Kpop_file']
    ## sort out output phenotype file
    if arguments['--output_file'] is None:
        output_file = CFG['REL-2016-01']['processed_pheno_file']
        #output_file = "REL-2016-02.RNAseq.trans.ipsc.meta_20160214.hdf5"
    else:
        output_file = arguments['--output_file']
    ## define file variables
    # phenoFile = os.path.join(input_dir, pheno_file)   #
    phenoFile = pheno_file
    outFile = os.path.join(output_dir, output_file)
    print "\nReading from: %s" %(phenoFile)
    print "\nReading Kpop from: %s" %(Kpop_file)
    print "\nOutputting data to: %s\n" %(outFile)
    sampleMetaFile = os.path.join(CFG['REL-2016-01']['base'], 
                                      "sample_metadata.tsv")
    geneMetaFile = os.path.join(CFG['REL-2016-01']['base'], 
                                      "gene_metadata.tsv")

    # get data from input file
    # sampleInfo = read_sample_meta(phenoFile)
    sampleInfo = pd.read_csv(sampleMetaFile, sep = "\t")
    phenoGene = read_pheno_data(phenoFile)
    # geneAnn = read_feature_annos(phenoFile)
    geneAnn = pd.read_csv(geneMetaFile, sep = "\t")
    Kpop = read_Kpop_data(Kpop_file)
    
    ## Filter out some samples
    ## Note. some missing cell type annotations .. all correspond to
    ## Filipa's lines -> exclude
    ## Yields 169 lines = 129 unique donors
    print ".. filtering samples for normal, fibroblast-only, sendai-only"
    flt = SP.ones(sampleInfo['name'].shape[0], dtype=bool)
    df_idx = pd.DataFrame({'disease': (sampleInfo['disease'] == "Normal"), 
             'derived': (sampleInfo['derived_from_cell_type'] == "Fibroblast"), 
             'repr_retro': (sampleInfo['reprogramming'] != "Retrovirus"), 
             'repr_episo': (sampleInfo['reprogramming'] != "Episomal")})
    flt_samp = df_idx.all(axis = 1)
    flt_array = np.array(flt_samp, dtype=bool)
    for key in phenoGene.keys():
        phenoGene[key] = phenoGene[key][flt_array,:]
    sampleInfo = sampleInfo[flt_array]
	## unique donor
    udonor = SP.unique(sampleInfo['donor'])
    N = sampleInfo['donor'].shape[0]
    print("After filtering, retain %s lines (corresponding to %s unique donors)"
          %(N, len(udonor)))
    ## Note. a couple of Filipa's lines still remain after this filter
    ## -> keep in for now, probably need to exclude later

    # Deal with zero counts
    Y_gene = phenoGene['exprs'] #  (249, 39293)
    ## expression values are log2(tpm + 1)
    Y = Y_gene
    ## sum up zero counts per gene
    counts = phenoGene['counts']
    ## take 5 as minimum counts to say observation is expressed
    zeros = SP.sum(counts < 5, axis = 0)
    prop_zero = zeros / float(Y.shape[0]) # zero_count_per_gene / n_donors
    flt = (prop_zero < max_na_pheno)
    print ".. filtering for zero gene counts (max_na: %s)" %(max_na_pheno)
    print ".. retaining %s genes" %(sum(flt))
    Y_gene = Y[:,flt] # (249, 17076)
    print ".. updating annotations accordingly\n"
    geneAnn = geneAnn[flt]

    # log and gaussianize:
    print "... log scaled counts and gaussianize expression gene-wise"
    ## converts the columns of Y (here: genes) to quantiles of standard normal
    # Y2_gene = SP.log10(Y_gene)     # log10(scaled_counts + 0.5)
    Y2_gene = Y_gene
    Y3_gene = gaussianize(Y_gene) # (249, 17076)
    # Y3_gene -= Y3_gene.mean(0)
    # Y3_gene /= Y3_gene.std(0)

    # plot
    if plot:
        print "... producing histograms of expression distributions"
        sns.set_style("white")
        plotFile = os.path.join(plotDir, 'preprocess_hist_gauss_counts.pdf')
        # gene
        mean1 = SP.mean(Y2_gene, axis=0)
        mean2 = SP.mean(Y3_gene, axis=0)
        f, (ax1, ax2) = PL.subplots(2, sharex=False, sharey=False)
        sns.distplot(mean1, ax=ax1, bins=20, kde=False)
        sns.despine()
        ax1.set_title('log10 scaled counts (mean)')
        sns.distplot(mean2, ax=ax2, bins=20, kde=False)
        sns.despine()
        ax2.set_title('log10 scaled gene counts + gaussianized (mean)')
        # tweak the title
        ttl1 = ax1.title
        ttl1.set_weight('bold')
        ttl2 = ax2.title
        ttl2.set_weight('bold')
        PL.figtext(0.01, 0.01, date.today().isoformat())
        PL.tight_layout()
        PL.savefig(plotFile)
        PL.close()

        # Produce a PCA plot of the samples
        covY = SP.cov(Y3_gene)
        eigenvals, eigenvecs = linalg.eigh(covY + 1e-6 * SP.eye(covY.shape[0]))
        eigenvals = eigenvals[::-1]
        eigenvecs = eigenvecs[::-1]
        df_pcs = pd.DataFrame({'PC1': eigenvecs[0], 'PC2': eigenvecs[1], 'PC3': eigenvecs[2], 'PC4': eigenvecs[3], 'PC5': eigenvecs[4]})
        ## PC pairs plot coloured by assay time
        ## just use date and time for assaytime
        print "... producing PC pairs plot coloured by assay time"
        assaytime = [str(dt).split(" ")[0][:-3] for dt in sampleInfo['assaytime_rnaseq']]
        df1 = pd.DataFrame({'assay_time': assaytime})
        df_to_plot = pd.concat([df_pcs, df1], axis = 1)
        ax = sns.pairplot(df_to_plot, hue = "assay_time")
        sns.despine()
        ax.savefig(os.path.join(plotDir, 'preprocess_rnaseq_gauss_counts_pairplot_pcs1-5_col_by_assaytime.pdf'))
        ## PC pairs plot coloured by assay user
        print "... producing PC pairs plot coloured by assay user"
        df1 = pd.DataFrame({'assay_user': sampleInfo['assayuser_rnaseq']})
        df_to_plot = pd.concat([df_pcs, df1], axis = 1)
        ax = sns.pairplot(df_to_plot, hue = "assay_user")
        sns.despine()
        ax.savefig(os.path.join(plotDir, 'preprocess_rnaseq_gauss_counts_pairplot_pcs1-5_col_by_assayuser.pdf'))
        ## PC pairs plot coloured by gender
        # print "... producing PC pairs plot coloured by gender"
        # df1 = pd.DataFrame({'gender': sampleInfo['gender']})
        # df_to_plot = pd.concat([df_pcs, df1], axis = 1)
        # ax = sns.pairplot(df_to_plot, hue = "gender")
        # sns.despine()
        # ax.savefig(os.path.join(plotDir, 'preprocess_rnaseq_gauss_counts_pairplot_pcs1-5_col_by_gender.pdf'))
        ## PC pairs plot coloured by growing condition
        # print "... producing PC pairs plot coloured by growing condition"
        # df1 = pd.DataFrame({'growing_conditions': sampleInfo['growing_conditions_rnaseq']})
        # df_to_plot = pd.concat([df_pcs, df1], axis = 1)
        # ax = sns.pairplot(df_to_plot, hue = "growing_conditions")
        # sns.despine()
        # ax.savefig(os.path.join(plotDir, 'preprocess_rnaseq_gauss_counts_pairplot_pcs1-5_col_by_growing_conditions.pdf'))
                
        # plot proportion of variance explained by PCs
        print "... plotting percentage of variance explained by PCs"
        sns.set_style("ticks")
        df = pd.DataFrame({'PC': np.arange(1, 50),'Pct_Variance_Explained': eigenvals[0:49] / eigenvals.sum() * 100})
        ax = df.plot(kind = "scatter", x = "PC", y = "Pct_Variance_Explained")  # df is an instance of DataFrame
        ax.set_ylim([0, 20])
        ax.set_ylabel("% variance explained")
        ax.set_xlim([0, 50])
        ax.set_xlabel("principal component")
        fig = ax.get_figure()
        sns.despine()
        fig.savefig(os.path.join(plotDir, 'preprocess_rnaseq_gauss_counts_var_expl_pcs.pdf'))

    # print out proportion of variance explained by PCs:
    ## no need to standardize as already normalized
    cum_var = PC_varExplained(Y3_gene, standardize=False) 
    print "Gene-based counts: "
    print "PCA (RNAseq log10 + gauss + std) - cumulative variance explained:"
    for i in [0, 1, 4, 9, 14, 19, 29, 39, 49]:
        print "PC%s: %s " %(i + 1, cum_var[i])

    if plot:
        # plot cumulative percentage of variance explained by PCs
        print "... plotting cumulative percentage of variance explained by PCs"
        sns.set_style("ticks")
        df = pd.DataFrame({'PC': np.arange(1, 50), 'Cumulative_Variance_Explained': cum_var[:49] * 100})
        ax = df.plot(kind = "scatter", x = "PC", y = "Cumulative_Variance_Explained")  # df is an instance of DataFrame
        ax.set_ylim([0, 100])
        ax.set_ylabel("cumulative % variance explained")
        ax.set_xlim([0, 50])
        ax.set_xlabel("principal component")
        fig = ax.get_figure()
        sns.despine()
        fig.savefig(os.path.join(plotDir, 'preprocess_rnaseq_gauss_counts_cum_var_expl_pcs.pdf'))

    # Define covariates:
    print "\n... defining covariates per line"
    gender_cov = 1. * (sampleInfo['gender']=='female')[:, SP.newaxis]
    ff_cov = (1. *
              (sampleInfo['growing_conditions_rnaseq'] == 'Feeder-free')[:,SP.newaxis])
    Cov = SP.concatenate([gender_cov, ff_cov], 1)
    Cov -= Cov.mean(0)
    Cov /= Cov.std(0)

    # Regress out covariates and obtain residuals
    ## regress out covar from Y:
    resids_raw = regressOut(Y3_gene, Cov)
    ## Gaussianize:
    resids_gauss = gaussianize(resids_raw)
    ## Standardize:
    resids_stand = (resids_gauss - resids_gauss.mean(0)) / resids_gauss.std(0)

    if plot:
        print "... producing histograms of residuals after covariate regression"
        sns.set_style("white")
        # mean residuals by gene
        plotFile = os.path.join(plotDir, 'preprocess_hist_expression_residuals_means.pdf')
        mean1 = SP.mean(resids_raw, axis=0)
        mean2 = SP.mean(resids_gauss, axis=0)
        mean3 = SP.mean(resids_stand, axis=0)   
        f, (ax1, ax2, ax3) = PL.subplots(3, sharex=False, sharey=False)
        sns.distplot(mean1, ax=ax1, bins=20, kde=False)
        sns.despine()
        ax1.set_title('mean raw residuals by gene')
        sns.distplot(mean2, ax=ax2, bins=20, kde=False)
        sns.despine()
        ax2.set_title('mean gaussianized residuals by gene')
        sns.distplot(mean3, ax=ax3, bins=20, kde=False)
        sns.despine()
        ax3.set_title('mean standardized gaussianized residuals by gene')        
        # tweak the title
        ttl1 = ax1.title
        ttl1.set_weight('bold')
        ttl2 = ax2.title
        ttl2.set_weight('bold')
        ttl3 = ax3.title
        ttl3.set_weight('bold')
        PL.figtext(0.01, 0.01, date.today().isoformat())
        PL.tight_layout()
        PL.savefig(plotFile)
        PL.close()

        # residual distributions for NANOG
        nanog_idx = (geneAnn['hgnc_symbol'] == 'NANOG').ravel().nonzero()[0][0]
        plotFile = os.path.join(plotDir, 'preprocess_hist_expression_residuals_nanog.pdf')
        f, (ax1, ax2, ax3) = PL.subplots(3, sharex=False, sharey=False)
        sns.distplot(resids_raw[:, nanog_idx], ax=ax1, bins=20, kde=False)
        sns.despine()
        ax1.set_title('raw residuals for NANOG')
        sns.distplot(resids_gauss[:, nanog_idx], ax=ax2, bins=20, kde=False)
        sns.despine()
        ax2.set_title('gaussianized residuals for NANOG')
        sns.distplot(resids_stand[:, nanog_idx], ax=ax3, bins=20, kde=False)
        sns.despine()
        ax3.set_title('standardized gaussianized residuals for NANOG')        
        # tweak the title
        ttl1 = ax1.title
        ttl1.set_weight('bold')
        ttl2 = ax2.title
        ttl2.set_weight('bold')
        ttl3 = ax3.title
        ttl3.set_weight('bold')
        PL.figtext(0.01, 0.01, date.today().isoformat())
        PL.tight_layout()
        PL.savefig(plotFile)
        PL.close()

    # export results
    print "\n exporting results to " + outFile
    counts = phenoGene['counts'][:,flt]
    exprs = phenoGene['exprs'][:,flt]
    tpm = phenoGene['tpm'][:,flt]
    ## drop problematic columns
    sampleInfo.drop('checks_gex_fail', axis=1, inplace=True)
    sampleInfo.drop('checks_qc2_swap', axis=1, inplace=True)
    sampleInfo.drop('study_blueprint', axis=1, inplace=True)
    sampleInfo.drop('study_reprogramming_comparison', axis=1, inplace=True)
    sampleInfo.drop('study_media_comparison', axis=1, inplace=True)

    export_results(outFile, sampleInfo, geneAnn, Kpop, counts,
                   exprs, tpm, resids_raw, resids_gauss, resids_stand)

    print "Done!"
    
################################################################################

# Fin.
