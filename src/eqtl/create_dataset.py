"""Create dataset(s) of HipSci RNA-seq and genotype data for trans-eQTL mapping

Usage: 
  create_dataset.py [-hc NAME] [--release=REL] [--nan_max=INT] [--maf_min=DBL] 
                    [--imp2_info_min==DBL] [--input_data_dir=DIR] 
                    [--output_data_dir=DIR] [--pheno_file=FILE] 
                    [--geno_file=FILE] [--geno_dosage_file=FILE] 
                    [--output_file=FILE] [--debug]
  create_dataset.py --version

Combine processed RNA-seq expression data and genotype data for a given 
chromosome from HDF5 files and save processed data to an HDF5 file for that 
chromosome.

Options:
  -h --help               show this help message and exit
  --version               show version and exit
  -c NAME --chrom=NAME    chromosome name
  -r REL --release=REL    name of the release of the data to be used (e.g. REL-2016-01)
  --nan_max=INT           maximum number of missing (NA) values acceptable 
                          to retain a sample for analysis (default=10)
  --maf_min=DBL           minimum minor allele frequency for genotyped variants
                          (default=0.01)
  --imp2_info_min=DBL     minimum IMPUTE2 info score required to retain 
                          genotyped variants (default=0.4)
  --input_data_dir=DIR    directory in which to find input data 
                          (default: ../../data/eQTL_rnaseq/)
  --output_data_dir=DIR   directory in which to find input data
  --pheno_file=FILE       filename for the input phenotype (expression) data
  --geno_file=FILE        filename for the genotypes
  --geno_dosage_file=FILE filename for the genotype dosages
  --output_file=FILE      filename for the output combined dataset file
  --debug                 run in debugging mode

For defaults for the input and output data, see the configuration settting in 
../CFG/settings.py

# HipSci REL-2014-11 eQTL map, imputed genotypes
# Step1: match expression data to genotypes, impute missing genotypes
# chromosome defined from command line
# 03/02/15 HK: adapted for genotype dosage
# 05/03/15 HK: filtering for IMP2_info score
# Updated 02/07/15 HK: for final versions
# Updated 13/10/15 for RNA-seq data and old/new genotypes
# 02/11/15 new genotypes only

Written by Davis McCarthy, February 2016
Based on earlier code written by Helena Kilpinen, February-November 2015
"""

################################################################################

## import packages
import scipy as SP
import h5py
import os
from docopt import docopt
from datetime import datetime
import sys
sys.path.append('./../')
from CFG.settings import CFG
sys.path.append('./../include')
from normalization import gaussianize
import pdb

################################################################################

def read_geno_data(genoFile, genoFile_dosage):
    """Read in gene-based phenotype data from an HDF5 file
    
    Args:
        genoFile (str): filename for the HDF5 file with standard genotypes
        genoFile_dosage (str): filename for the HDF5 file with dosage genotypes

    Returns:
        dictionary with genotype data
    """
    # import standard genotype info for old ond new genotypes
    geno_data = {}
    fGeno = h5py.File(genoFile) # single-sample imputation
    std_ids = fGeno['genotype']['row_header']['sample_ID'][:]
    std_pos = fGeno['genotype']['col_header']['pos'][:]
    geno_data['X'] = fGeno['genotype']['matrix'][:]
    geno_data['sampleID'] = fGeno['genotype']['row_header']['sample_ID'][:]
    geno_data['chrom'] = fGeno['genotype']['col_header']['chrom'][:]
    geno_data['pos'] = fGeno['genotype']['col_header']['pos'][:]
    fGeno.close()
    # deal with 255 values (if any)
    geno_data['X'] = geno_data['X'].astype(float)
    geno_data['X'][geno_data['X']==255] = -1
    # import genotype dosage info
    fGeno = h5py.File(genoFile_dosage) # single-sample imputation (original)
    dosage_ids = fGeno['genotype']['row_header']['sample_ID'][:]
    dosage_pos = fGeno['genotype']['col_header']['pos'][:]
    geno_data['Xd'] = fGeno['genotype']['matrix'][:]
    geno_data['gdid'] = fGeno['genotype']['col_header']['gdid'][:]
    geno_data['info'] = fGeno['genotype']['col_header']['IMP2_info'][0]
    fGeno.close()
    # deal with 255 values (if any)
    geno_data['Xd'] = geno_data['Xd'].astype(float)
    geno_data['Xd'][geno_data['Xd']==255] = -1
    geno_data['Xd'] = geno_data['Xd'].T # to match dimensions of X
    ## check that the two genotype sets match!
    if set(std_ids)==set(dosage_ids):
        print "Sample order in the two genotype matrices (std/dosage) match"
    else:
        raise Exception('Sample order in the two genotype matrices (std/dosage) do not match')
    if set(std_pos)==set(dosage_pos):
        print "Variant positions in the two genotype matrices (std/dosage) match"
    else:
        raise Exception('Variant positions in the two genotype matrices (std/dosage) do not match')
    ## return genotype data
    return geno_data

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
    pheno_gene['residuals_raw'] = fpheno['phenotype_gene']['residuals_raw'][:]
    pheno_gene['residuals_gaussianized'] = fpheno['phenotype_gene']['residuals_gaussianized'][:]
    pheno_gene['residuals_standardized'] = fpheno['phenotype_gene']['residuals_standardized'][:]
    pheno_gene['counts'] = fpheno['phenotype_gene']['counts'][:]
    if '/phenotype_gene/processed_counts' in fpheno:
        pheno_gene['processed_counts'] = fpheno['phenotype_gene']['processed_counts'][:]
    if '/phenotype_gene/exprs' in fpheno:
        pheno_gene['exprs'] = fpheno['phenotype_gene']['exprs'][:]
    if '/phenotype_gene/tpm' in fpheno:
        pheno_gene['tpm'] = fpheno['phenotype_gene']['tpm'][:]
    if '/phenotype_gene/peer_residuals_n10' in fpheno:
        pheno_gene['peer_residuals_n10'] = fpheno['phenotype_gene']['peer_residuals_n10'][:]
    if '/phenotype_gene/peer_residuals_n15' in fpheno:
        pheno_gene['peer_residuals_n15'] = fpheno['phenotype_gene']['peer_residuals_n15'][:]
    if '/phenotype_gene/peer_residuals_n20' in fpheno:
        pheno_gene['peer_residuals_n20'] = fpheno['phenotype_gene']['peer_residuals_n20'][:]
    fpheno.close()
    return pheno_gene

################################################################################

def read_Kpop_data(phenoFile):
    """Read in gene-based phenotype data from an HDF5 file
    
    Args:
        phenoFile (str): filename for the relevant HDF5 file

    Returns:
        dictionary with gene-level phenotype data
    """
    fKpop = h5py.File(phenoFile,'r')
    Kpop = {}
    for key1 in fKpop['Kpop']:
        Kpop[key1] = {}
        for key2 in fKpop['Kpop'][key1]:
            Kpop[key1][key2] = fKpop['Kpop'][key1][key2][:]
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
    fpheno = h5py.File(phenoFile,'r')
    # Sample meta (same for all):
    sampleInfo = {}
    for key in fpheno['sample_info'].keys():
    	sampleInfo[key] = fpheno['sample_info'][key][:]
    fpheno.close()
    return sampleInfo

################################################################################

def read_gene_annos(phenoFile):
    """Read in gene-based metadata from an HDF5 file
    
    Args:
        phenoFile (str): filename for the relevant HDF5 file

    Returns:
        dictionary with feature annotations
    """
    fpheno = h5py.File(phenoFile,'r')
    # Feature annotations:
    geneAnn = {}
    for key in fpheno['gene_info'].keys():
    	geneAnn[key] = fpheno['gene_info'][key][:]
    fpheno.close()
    return geneAnn

################################################################################

def matchID(id1, id2):
	""" match id1 and id2 """
	idx1 = []
	idx2 = []
	for i in range(id1.shape[0]):
		#print 1.*i/float(id1.shape[0])
		if id1[i] in id2.tolist():
			idx1.append(i)
			idx2.append(SP.where(id2==id1[i])[0][0])
	idx1 = SP.array(idx1)
	idx2 = SP.array(idx2)
	print (id1[idx1] == id2[idx2]).all()  # assert that's right
	return idx1, idx2

################################################################################

def export_results(outFile, sampleInfo, geneAnn, phenoGene, genoData, Kpop):
    """Export results to an HDF5 file

    Args:
        outFile (str): filename for the output HDF5 file
        sampleInfo (hdf5): hdf5 connection to and HDF5 group containing sample 
                           info
        geneAnn (hdf5): hdf5 connection to and HDF5 group containing gene
                           annotation info
        phenoGene (dbl): dict containing matrices of expression values
        genoData (dbl): dict containing genotype matrices and SNP information
                           
    Returns:
        null

    This function is designed only to save gene-level expression data as exon 
    and probe level expression data is not as relevant for trans-eQTL mapping.
    """
    ### Export metadata (same for all):
    print "..exporting metadata:"
    fOut = h5py.File(outFile, 'a')
    ## add/replace sample metadata information
    if "/sample_info" in fOut:
        del fOut['sample_info']
    gMETA = fOut.create_group('sample_info')
    for key in sampleInfo.keys():
        gMETA.create_dataset(key, data = sampleInfo[key])
    ## add/replace gene metadata information
    if "/gene_info" in fOut:
        del fOut['gene_info']
    gFEATURE = fOut.create_group('gene_info')
    for key in geneAnn.keys():
        gFEATURE.create_dataset(key, data=geneAnn[key])
    ## add/replace gene expression data
    if "/phenotype_gene" in fOut:
        del fOut['phenotype_gene']
    gGENE = fOut.create_group('phenotype_gene')
    for key in phenoGene.keys():
        gGENE.create_dataset(key, data=phenoGene[key])
    ## add/replace SNP info
    if "/snp_info" in fOut:
        del fOut['snp_info']
    gSNP = fOut.create_group('snp_info')
    for key in genoData.keys():
        if key not in ['X', 'Xd']:
            gSNP.create_dataset(key, data=genoData[key])
    ## add/replace genotype matrices
    if "/genotype" in fOut:
        del fOut['genotype']
    gGeno = fOut.create_group('genotype')
    gGeno.create_dataset('matrix', data=genoData['X'])
    gGeno.create_dataset('matrix_dosage', data=genoData['Xd'])
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

if __name__ == '__main__':
    arguments = docopt(__doc__, version = "0.0.1")
    # ----------------------------------------------
    # Parse command line arguments and define parameters
    ## sort out chromosome
    chrom = arguments['--chrom']
    if chrom:
        chrom = str(chrom)
    else:
        chrom = '1'
    print("Creating dataset for chromosome " + chrom)
    ## sort out release
    if arguments['--release'] == 'REL-2016-01':
        release = 'REL-2016-01'
    else:
        release = 'rnaseq_Jan-16'
    ## sort out max number of missing values allowed
    if arguments['--nan_max']:
        nan_max = float(arguments['--nan_max'])
    else:
        nan_max = 10
    print("Maximum acceptable number of missing values is %2f" %nan_max)
    ## sort out minimum MAF
    if arguments['--maf_min']:
        maf_min = float(arguments['--maf_min'])
    else:
        maf_min = 0.01
    print("Minimum minor allele frequency is %d" %maf_min)
    ## sort out minimum IMPUTE2 info score
    if arguments['--imp2_info_min']:
        imp2_info_min = float(arguments['--imp2_info_min'])
    else:
        imp2_info_min = 0.4
    print("Minimum IMPUTE2 info score thrshold is %2f" %imp2_info_min)
    ## sort out input directory
    if arguments['--input_data_dir'] is None:
        input_dir = CFG[release]['base']
    else:
        input_dir = arguments['--input_data_dir']
    print "Input data directory is " + input_dir
    ## sort out output directory
    if arguments['--output_data_dir'] is None:
        output_dir = CFG[release]['base']
    else:
        output_dir = arguments['--output_data_dir']
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    ## sort out input phenotype file
    if arguments['--pheno_file'] is None:
        #pheno_file = 'REL-2016-02.RNAseq.trans.ipsc.meta_20160214.hdf5'
        pheno_file = CFG[release]['processed_pheno_file']
    else:
        pheno_file = arguments['--pheno_file']
    ## sort out input genotype file
    if arguments['--geno_file'] is None:
        # Genotypes: single-sample imputation
        # geno_dir = "/nfs/research/stegle/projects/hipsci/data/genotypes/imputed/REL-2014-11_SS/hdf5"
        # geno_file = "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.20151005.genotypes.chr" + chrom + ".gdid.mac1.recode.GL.hdf5"
        # genoFile = os.path.join(geno_dir, geno_file)
        genoFile = CFG['genotypes']['standard'][chrom]
    else:
        genoFile = arguments['--geno_file']
    ## sort out input genotype dosage file
    if arguments['--geno_dosage_file'] is None:
        # Genotypes: single-sample imputation
        genoFile_dosage = CFG['genotypes']['dosage'][chrom]
    else:
        genoFile_dosage = arguments['--geno_dosage_file']
    ## sort out output phenotype file
    if arguments['--output_file'] is None:
        output_file = release + ".RNAseq.trans.ipsc.data.chr" + chrom + ".hdf5"
    else:
        output_file = arguments['--output_file']
    ## define file variables
    phenoFile = os.path.join(input_dir, pheno_file)
    outFile = os.path.join(output_dir, output_file)
    print "\nReading expression data and metadata from: %s" %(phenoFile)
    print "\nReading standard genotype data from: %s" %(genoFile)
    print "\nReading genotype dosage data from: %s" %(genoFile_dosage)
    print "\nOutputting data to: %s\n" %(outFile)
    ## check if we should be in debug mode
    if arguments['--debug'] is None:
        debug = False
    else:
        debug = True

    # ----------------------------------------------
    # get expression, gene and sample data from input pheno file
    sampleInfo = read_sample_meta(phenoFile)
    geneAnn = read_gene_annos(phenoFile)
    phenoGene = read_pheno_data(phenoFile)
    Kpop = read_Kpop_data(phenoFile)

    # gene strand info missing - get from original gencode
    # print ".. extracting gene strand from original genocode v19 annotations"
    # GC = SP.loadtxt(CFG['gencode_genes'], delimiter = '\t', dtype='str')
    # id1 = GC[1::, 5]
    # strand = GC[1::, 3]
    # id2 = SP.array([geneAnn['geneID'][x].split('.')[0] for x in range(geneAnn['geneID'].shape[0])])
    # target_strand = []
    # for p in range(id2.shape[0]):
    #     idx = SP.where(id1==id2[p])
    #     target_strand.append(strand[idx][0])
    # geneAnn['strand'] = SP.array(target_strand)

    # get standard genotype info
    genoData = read_geno_data(genoFile, genoFile_dosage)

    # match genotypes with phenotype, both datasets
    print "\nMatching genotypes with phenotypes..."
    if release == "REL-2016-01":
        sampleInfo['sampleID'] = sampleInfo['name']
        geneAnn['geneID'] = geneAnn['ensembl_gene_id']
    idx_p, idx_g = matchID(sampleInfo['sampleID'], genoData['sampleID'])
    genoData['X'] = genoData['X'][idx_g,:] # yields 167 lines; n_variants still differ
    ## yields 237 lines for REL-2016-01 data
    genoData['Xd'] = genoData['Xd'][idx_g,:]
    genoData['sampleID'] = genoData['sampleID'][idx_g]
    for key in phenoGene.keys():
        phenoGene[key] = phenoGene[key][idx_p,:]
    for key in sampleInfo.keys():
        sampleInfo[key] = sampleInfo[key][idx_p]
    print ".. dataset reduced to %s lines" %(sampleInfo['sampleID'].shape[0])
    print ".. finished matching sample ids between genotype and phenotype"

    ## standardize PEER residuals
    for nfactors in ['10', '15', '20']:
        peer_name = 'peer_residuals_n' + nfactors
        if peer_name in phenoGene:
            peer_stand_name = 'peer_residuals_n' + nfactors + 'std'
            ## Gaussianize:
            peer_gauss = gaussianize(phenoGene[peer_name])
            ## Standardize:
            peer_stand = (peer_gauss - peer_gauss.mean(0)) / peer_gauss.std(0)
            ## add to phenoGene
            phenoGene[peer_stand_name] = peer_stand

    # order all datasets by donorID
    print ".. ordering datasets by donorID"
    idx = SP.argsort(sampleInfo['donor'])
    for key in phenoGene.keys():
        phenoGene[key] = phenoGene[key][idx,:]
    for key in sampleInfo.keys():
        sampleInfo[key] = sampleInfo[key][idx]
    genoData['X'] = genoData['X'][idx,:]
    genoData['Xd'] = genoData['Xd'][idx,:]
    genoData['sampleID'] = genoData['sampleID'][idx]

    # -------------------------------------------------------
    # Assumes that genotype and genotype dosage information is matching:
    # MAF/NA based on the non-dosage genotype set
    # read_geno_data() should raise an exception if this is not the case
    # --------------------------------------------------------

    # ----------------------------------------------
    # filter and impute SNPs
    print ".. starting variant processing"
    print ".. filtering by MAF (%.2f) and percentage of missing values (%.2f)" %(maf_min,nan_max)
    # replace neg entries with nans
    genoData['X'][genoData['X'] < 0] = float('nan')
    Nvals = (~SP.isnan(genoData['X'])).sum(0)
    MAF = ((2 * SP.sum(genoData['X']==2, axis=0) +
           SP.sum(genoData['X']==1, axis=0)) / (2.*Nvals))
    MAF[SP.isnan(MAF)] = 0
    # invert minor and major if MAF > 0.5
    is_to_invert = MAF > 0.5
    Xi = genoData['X']
    for s in range(Xi.shape[1]):
        if is_to_invert[s]:
            Xi[Xi[:,s]==0, s] = 2
            Xi[Xi[:,s]==2, s] = 0
    # recalculate MAF after correction
    MAF = (2 * SP.sum(Xi==2, axis=0) + SP.sum(Xi==1, axis=0)) / (2. * Nvals)
    MAF[SP.isnan(MAF)] = 0
    genoData['maf'] = MAF
    # filter based on number of missing genotypes
    perc_miss = SP.isnan(genoData['X']).sum(0) / float(genoData['X'].shape[0])
    filt_var = (perc_miss < nan_max) * (MAF > maf_min) 
    genoData['X'] = genoData['X'][:,filt_var]
    genoData['Xd'] = genoData['Xd'][:,filt_var]
    genoData['info'] = genoData['info'][filt_var]
    genoData['pos'] = genoData['pos'][filt_var]
    genoData['chrom'] = genoData['chrom'][filt_var]
    genoData['maf'] = genoData['maf'][filt_var]
    genoData['gdid'] = genoData['gdid'][filt_var]
    print '%d out of %d SNPs have survived after MAF/NaN filter' %(filt_var.sum(), filt_var.shape[0])
    print 'allele encoding is 0 for REF allele, 1 for ALT allele'

    # filter out sites below defined imputation score
    print ".. filtering for imputation score"
    genoData['info'] = genoData['info'].astype(float)
    filt_var2 = SP.ones(genoData['info'].shape[0], dtype=bool)
    for i in range(genoData['info'].shape[0]):
        I = genoData['info'][i] > imp2_info_min
        filt_var2[i] = I
    genoData['X'] = genoData['X'][:,filt_var2]
    genoData['Xd'] = genoData['Xd'][:,filt_var2]
    genoData['info'] = genoData['info'][filt_var2]
    genoData['pos'] = genoData['pos'][filt_var2]
    genoData['chrom'] = genoData['chrom'][filt_var2]
    genoData['maf'] = genoData['maf'][filt_var2]
    genoData['gdid'] = genoData['gdid'][filt_var2]
    print '%d out of %d SNPs have survived after IMP2_info filter' %(filt_var2.sum(), filt_var2.shape[0])

    # Impute missing values for standard genotypes and dosages
    print ".. mean imputing SNPs"
    for s in range(genoData['X'].shape[1]):
        is_nan = SP.isnan(genoData['X'][:,s])
        genoData['X'][is_nan,s] = SP.mean(genoData['X'][~is_nan,s])
    print ".. mean imputing SNPs - dosage"
    for sd in range(genoData['Xd'].shape[1]):
        is_nan = SP.isnan(genoData['Xd'][:,sd])
        genoData['Xd'][is_nan,sd] = SP.mean(genoData['Xd'][~is_nan,sd])
	
    # order SNPs by chromosome and position
    print ".. ordering chrom and pos"
    pos_max = genoData['pos'].max()
    tmp_chrom = genoData['chrom']
    tmp_chrom[tmp_chrom == 'X'] = 23
    tmp_chrom = SP.array(tmp_chrom, dtype=int)
    _pos = genoData['pos']+ pos_max * tmp_chrom
    idx = SP.argsort(_pos)
    genoData['X'] = genoData['X'][:,idx]
    genoData['Xd'] = genoData['Xd'][:,idx]
    genoData['info'] = genoData['info'][idx]
    genoData['pos'] = genoData['pos'][idx]
    genoData['chrom'] = genoData['chrom'][idx]
    genoData['maf'] = genoData['maf'][idx]
    genoData['gdid'] = genoData['gdid'][idx]

    # ----------------------------------------------
    # Write out final hdf5
    print ".. exporting data to file:"
    print outFile
    export_results(outFile, sampleInfo, geneAnn, phenoGene, genoData, Kpop)
    print "Done!"
    
################################################################################

# Fin.
