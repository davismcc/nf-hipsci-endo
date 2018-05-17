import os
import vcf
import h5py
import numpy as np
import pandas as pd
from datetime import datetime
from cyvcf2 import VCF

# We can use the basic structure/format of the HDF5 genotype files produced by Helena (e.g. `/nfs/research2/stegle/projects/hipsci/data/eQTL/Dec-16/OUT/REL-2016-09.STAR_HTSeq.data.151216.chr1.hdf5`).

## define files and samples
geno_vcf_dir = "/hps/nobackup/stegle/projects/hipsci/hipsci_genotypes/REL-2016-09"
vcf_file = "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2016-09.imputed_phased.INFO_0.4_filtered.20160912.genotypes.chr22.vcf.gz"
h5_file = "hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2016-09.imputed_phased.INFO_0.4_filtered.20160912.genotypes.chr22.h5"
vcf_reader = vcf.Reader(filename=os.path.join(geno_vcf_dir, vcf_file), compressed=True)
samples = np.array([a.encode('utf8') for a in vcf_reader.samples])

## read in data from vcf
refs = {}
alts = {}
chroms = {}
poss = {}
starts = {}
ends = {}
ids = {}
gt_types = {}
for variant in VCF(os.path.join(geno_vcf_dir, vcf_file)): # or VCF('some.bcf')
    if variant.is_snp and variant.ID is not None:
        var_id = str(variant.ID, "utf-8")
        refs[var_id] = variant.REF
        alts[var_id] = variant.ALT[0] # e.g. REF='A', ALT=['C', 'T']
        chroms[var_id] = variant.CHROM
        poss[var_id] = variant.POS
        starts[var_id] = variant.start
        ends[var_id] = variant.end
        ids[var_id] = var_id
        #     # numpy arrays of specific things we pull from the sample fields.
        #     # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
        gts = np.array(variant.gt_types)
        # To make things more understandable change the UNKNOWN value from 2 to nan
        # and the HOM_ALT value from 3 to 2
        gts = gts.astype(float)
        gts[gts == 2] = np.nan
        gts[gts == 3] = 2.0
        gt_types[var_id] = gts

## put data into dataframes for convenience
snp_info_df = pd.DataFrame({'ID': ids, 'CHROM': chroms, 'POS': poss, 'START': starts, 'END': ends, 'REF': refs, 'ALT': alts})
gt_df = pd.DataFrame(gt_types)

# # Write to HDF5 file
fout = h5py.File(os.path.join(geno_vcf_dir, h5_file))
fout.attrs['time_stamp'] = datetime.now().isoformat()
fout.attrs['author'] = 'Davis McCarthy'
fout.attrs['institution'] = 'EMBL-EBI'
# Add sampleIDs to HDF5 file.
fout.create_dataset('sampleID', data=samples)
fout.create_group('sample_info')
fout.create_dataset('sample_info/sampleID', data=samples)
# Add genotype matrix.
fout.create_group('genotype')
fout.create_dataset('/genotype/matrix', data=gt_df.values)
# Add SNP Info to HDF5 file.
fout.create_group(name='snp_info')
fout.create_dataset('/snp_info/alt', data=np.array([a.encode('utf-8') for a in snp_info_df['ALT'].values]))
fout.create_dataset('/snp_info/ref', data=np.array([a.encode('utf-8') for a in snp_info_df['REF'].values]))
fout.create_dataset('/snp_info/chrom', data=np.array([a.encode('utf-8') for a in snp_info_df['CHROM'].values]))
fout.create_dataset('/snp_info/pos', data=snp_info_df['POS'].values)
fout.create_dataset('/snp_info/start', data=snp_info_df['START'].values)
fout.create_dataset('/snp_info/end', data=snp_info_df['END'].values)
fout.create_dataset('/snp_info/gdid', data=np.array([a.encode('utf-8') for a in snp_info_df['ID'].values]))
fout.close()



