"""merge single cell RNA-Seq fastq files"""

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
import pysam

python_cmd = 'python'
# settings
nthreads = 1
sleep = 0
mem_thread = 20000

genotype_file = '/nfs/research/stegle/projects/hipsci/data/genotypes/gtarray/REL-2014-05/genotype.hdf5'
cluster_cmd = 'bsub -n %d -q research-rh6 -R "rusage[mem=%d,tmp=10000]" -M %d -o ./cluster_out' % (nthreads,mem_thread,mem_thread)


if __name__ == '__main__':
    if 'cluster' in sys.argv:
        #run merge jobs
        expr_file = sys.argv[4]
        dir_name_out = sys.argv[3]
        dir_name_in = sys.argv[2]
        if not os.path.exists(dir_name_out):
            os.makedirs(dir_name_out)
        FL = glob.glob(os.path.join(dir_name_in,'*.bam'))
        #match sample numbers
        #FL = FL[0:1]
        for fn in FL:
            out_file = os.path.join(dir_name_out,'%s' % os.path.basename(fn)+'.hdf5')
            if os.path.exists(out_file):
                print "skipping file: %s" % (out_file)
                continue
            cmd = '%s %s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],fn,out_file,expr_file)
            print cmd
            #cmd = '%s %s %s %s' % (python_cmd, sys.argv[0],fn,out_file)
            os.system(cmd)
            time.sleep(sleep)
    elif 'collect' in sys.argv:
        dir_name_out = sys.argv[2]
        FL = glob.glob(os.path.join(dir_name_out,'sc*.hdf5'))
        geno_id = None
        sample_ids = []
        geno_id_map = []
        error = []
        for fn in FL:
            f = h5py.File(fn,'r')
            _sample_id = f.keys()[0]
            sample_ids.append(_sample_id)
            geno_id = f[_sample_id]['geno_id'][:]
            geno_id_map.append(f[_sample_id]['geno_id_map'][0])
            _error = f[_sample_id]['error'][:]
            if (len(_error.shape)==1):
                _error = SP.zeros(len(geno_id))
            else:
                _error = _error.sum(axis=0)
            error.append(_error)
            pass
        error = SP.array(error)
        geno_id_map = SP.array(geno_id_map)
        geno_id     = SP.array(geno_id)
        geno_id_map = geno_id[geno_id_map]
        sample_ids = SP.array(sample_ids,dtype='str')
        fout = h5py.File(os.path.join(dir_name_out,'variants.hdf5'),'w')
        #fout.create_dataset(name='error',data=error)
        fout.create_dataset(name='sample_ids',data=sample_ids)
        fout.create_dataset(name='geno_id_map',data=geno_id_map)
        fout.create_dataset(name='geno_id',data=geno_id)
        fout.create_dataset(name='error',data=error)
    else:
        #carry out computation
        in_file  = sys.argv[1]
        expr_file = sys.argv[3]
        out_file = sys.argv[2]
        sample_id = os.path.basename(in_file).split('.')[0]

        #load expression
        genes = h5py.File(expr_file,'r')['genes']


        #gnotypes present in the experiment
        geno_id     = 'coxy,xavk,ffdm_2,ffdj_1'.split(',')
        #geno_id     = 'coxy,xavk,veve'.split(',')
        #geno_id     = 'fsfe10a,fsps13b,coxy,xavk'.split(',')

        base_name = os.path.splitext(os.path.basename(in_file))[0]
        base_name=string.replace(base_name,'#','')

        fgeno = h5py.File(genotype_file,'r') 
        sample_ids = fgeno['genotype/row_header/sample_ID'][:]
        alleles = fgeno['genotype/col_header']['alleles'][:]
        chrom = fgeno['genotype/col_header']['chrom'][:]
        pos = fgeno['genotype/col_header']['pos'][:]
        #remove header 
        sample_ids=SP.array([re.match('HP.*-(.*)',s).group(1) for s in sample_ids])
        #mattch
        if 0:
            geno_id = sample_ids 
        Igeno = SP.array([SP.nonzero(sample_ids==g)[0][0] for g in geno_id])
        Mgeno = fgeno['genotype']['matrix'][:][Igeno,:]
        #get variable positions
        Ivar = Mgeno.std(axis=0)>0
        Mgeno = Mgeno[:,Ivar]
        pos = pos[Ivar]
        chrom = chrom[Ivar]
        alleles = alleles[Ivar]
       
        #determine most highly covered 100 genes and use these for sample identification
        counts = genes['counts'][:]
        annotation = genes['annotation'][:]
        LS     = SP.array(counts.sum(axis=0),dtype='float')
        counts = counts/LS
        Mcounts = counts.mean(axis=1)
        Mcounts = counts[:,SP.nonzero(genes['sample_ids'][:]==sample_id)[0][0]]

        #top expressed genes
        #Itop = Mcounts.argsort()[::-1][0:500]
        Itop = Mcounts.argsort()[::-1][0:5000]

        #open file
        f = pysam.Samfile( in_file, "rb" )


        #iterate over for 1,000 positions
        #total count
        count = [] 
        #percent that could be mathced to a variant
        mapped = []
        #squared error to expectation of either line
        error = []
        
        print "start"
        for i_gene in Itop:
            gene_id = genes['gene_ids'][i_gene]
            if not re.match('ENSG.*',gene_id):
                continue
            anno = annotation[i_gene]
            _chrom = re.match('chr(.*)',anno[0]).group(1)
            _pos0  = int(anno[1])
            _pos1  = int(anno[2])
            _strand = anno[3]
            try:
                Ivar = chrom==int(_chrom)
                Ivar = Ivar & (pos>=_pos0) & (pos<_pos1)
            except Exception:
                continue
            #select region 
            if Ivar.sum()==0:
                continue
            for iv in SP.nonzero(Ivar)[0]:
                _pos = int(pos[iv])
                _alleles = alleles[iv]
                _geno    = Mgeno[:,iv]
                fetch = f.fetch('chr%s'% _chrom,_pos,_pos+1)
                L = []
                #create pileup representation of reads
                for read in fetch:
                    _positions = SP.array(read.positions)+1-_pos
                    _seq       = SP.asarray(list(read.seq))
                    if (read.is_qcfail or read.is_duplicate or (not read.is_proper_pair) or read.is_secondary):
                        #ski QC failtures or secondary alignments, and all taht
                        continue
                    #matching positions?
                    Iok = SP.nonzero(_positions==0)[0]
                    if len(Iok)>0:
                        L.append(_seq[Iok])
                pass
                L = SP.array(L)
                #pdb.set_trace()
                if len(L)==0:
                    continue
                #total count
                _C = L.shape[0]
                count.append(_C)
                P = SP.array([(L==_alleles[0]).sum(),(L==_alleles[1]).sum()],dtype='float')
                mapped.append(P.sum()/_C)
                P /= P.sum()
                G = SP.zeros([len(geno_id),2])
                for ig in xrange(len(_geno)):
                    #reference per genotype: 
                    if _geno[ig]==0:
                        G[ig,0] = 1
                    elif _geno[ig]==2:
                        G[ig,1] = 1
                    else:
                        G[ig,:] = 1
                _error = ((G-P)**2).sum(axis=1)
                error.append(_error)
                pass
        error = SP.array(error)
        mapped = SP.array(mapped)
        count  = SP.array(count)
        Iok = (mapped>0.5) & (count>100)
        #Iok = (mapped>0.2) & (count>10)
        count = count[Iok]
        error = error[Iok]
        mapped=mapped[Iok]

        geno_id_map = error.sum(axis=0).argmin()
        geno_map    = geno_id[geno_id_map]
        
        #write to output file
        fout = h5py.File(out_file,'w')
        group = fout.create_group(os.path.basename(in_file).split('.')[0])
        group.create_dataset(name='error',data=error)
        group.create_dataset(name='count',data=count)
        group.create_dataset(name='mapped',data=mapped)
        group.create_dataset(name='geno_id',data=geno_id)
        group.create_dataset(name='geno_id_map',data=[geno_id_map])
