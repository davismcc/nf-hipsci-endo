"""count expression matrices using htseq"""

import scipy as SP
import io
import sys
import os
import pdb
import re
import time
import glob
import string
import h5py
import cPickle

python_cmd = 'python'
nthreads = 1
sleep = 20
mem_thread = 10000
reference_gtf = '/homes/stegle/research/datasets/references/human/GRCh37/gencode.v19.annotation_ERCC.gtf'

cluster_cmd = 'bsub -n %d -q research-rh6 -R "rusage[mem=%d,tmp=10000]" -M %d -o ./cluster_out' % (nthreads,mem_thread,mem_thread)

htseq_count_cmd = 'htseq-count --format=bam --order=pos --mode=intersection-strict --stranded=no'
htseq_qa_cmd = 'htseq-qa --type=bam'
samtools_cmd = 'samtools'


if __name__ == '__main__':
    if 'cluster' in sys.argv:
        #run merge jobs
        dir_name_out = sys.argv[3]
        dir_name_in = sys.argv[2]
        if not os.path.exists(dir_name_out):
            os.makedirs(dir_name_out)
        FL = glob.glob(os.path.join(dir_name_in,'*.bam'))
        for fn in FL:
            out_file = os.path.join(dir_name_out,'%s' % os.path.basename(fn))
            if os.path.exists(out_file+'.txt'):
                print "skipping file: %s" % (out_file)
                continue
            cmd = '%s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],fn,out_file)
            print cmd
            os.system(cmd)
            time.sleep(sleep)
    elif 'collect' in sys.argv:
        #create extendd hdf5 output 
        dir_name_in = sys.argv[2]
        fout = h5py.File(os.path.join(dir_name_in,'counts.hdf5'),'w')
        fout.create_group(name='genes')
        fout.create_group(name='exons')
 
        # add annotation
        # Load gtf file
        exon_id_pat = re.compile(r'.*exon_id\s+"(.*?)"')
        gene_id_pat = re.compile(r'.*gene_id\s+"(.*?)"')
        transcript_id_pat = re.compile(r'.*transcript_id\s+"(.*?)"')
        pickl_file = reference_gtf + '.pickle'

        if not os.path.exists(pickl_file) or 'recalc' in sys.argv:
            exons = {}
            genes = {}
            with open(reference_gtf) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.split('\t')
                    if fields[2]=='gene':
                        gene_id = gene_id_pat.match(fields[-1]).group(1)
                        chrom   = fields[0]
                        start   = int(fields[3])
                        stop    = int(fields[4])
                        strand   = fields[6]
                        genes[gene_id] = [chrom,start,stop,strand]

                    if fields[2]=='exon':
                        if exon_id_pat.search(fields[-1]) is not None:
                            exon_id = exon_id_pat.search(fields[-1]).group(1)
                        else:
                            exon_id = None
                        gene_id = gene_id_pat.match(fields[-1]).group(1)
                        transcript_id = transcript_id_pat.match(fields[-1]).group(1)
                        chrom   = fields[0]
                        start   = int(fields[3])
                        stop    = int(fields[4])
                        strand   = fields[6]
                        exons[exon_id] = [chrom,start,stop,strand,gene_id,transcript_id]
            cPickle.dump([genes,exons],open(pickl_file,'wb'))
        else:
            [genes,exons] = cPickle.load(open(pickl_file,'rb'))

        #collect all dataset and summarizied in hdf5 files
        FL = glob.glob(os.path.join(dir_name_in,'*.bam.quant.gene.txt'))
        sample_ids = []
        counts     = []
        gene_ids = None
        for fn in FL:
            print fn
            bn = os.path.basename(fn)
            sample_id = bn.split('.')[0]
            sample_ids.append(sample_id)
            M = SP.loadtxt(fn,dtype='str')
            gene_ids = M[:,0]
            _counts = SP.array(M[:,1],dtype='int')
            counts.append(_counts)
            pass
        counts = SP.array(counts).T
        sample_ids = SP.array(sample_ids)
        
        #create tabular output as matrix:
        fgenes = fout['genes']
        fgenes.create_dataset(name='counts',data=counts)
        fgenes.create_dataset(name='gene_ids',data=gene_ids)
 
        #gene annotation
        gene_annotation = []
        for gene_id in gene_ids:
            if re.match('ENS.*',gene_id) is not None:
                gene_annotation.append(genes[gene_id]) 
            else:
                gene_annotation.append([0,0,0,0])
        fgenes.create_dataset(name='annotation',data=SP.array(gene_annotation))
        fgenes.create_dataset(name='sample_ids',data=sample_ids)
      

        #exons
        #collect all dataset and summarizied in hdf5 files
        FL = glob.glob(os.path.join(dir_name_in,'*.bam.quant.exon.txt'))
        sample_ids = []
        counts     = []
        exon_ids = None
        for fn in FL:
            print fn
            bn = os.path.basename(fn)
            sample_id = bn.split('.')[0]
            sample_ids.append(sample_id)
            M = SP.loadtxt(fn,dtype='str')
            exon_ids = M[:,0]
            _counts = SP.array(M[:,1],dtype='int')
            counts.append(_counts)
            pass
        counts = SP.array(counts).T
        sample_ids = SP.array(sample_ids)
        
        #create tabular output as matrix:
        fexons = fout['exons']
        fexons.create_dataset(name='counts',data=counts)
        fexons.create_dataset(name='exon_ids',data=exon_ids)
    
        #exon annotation
        exon_annotation = []
        for exon_id in exon_ids:
            if (re.match('ENS.*',exon_id) is not None):
                exon_annotation.append(genes[exon_id]) 
            else:
                exon_annotation.append([0,0,0,0])
        fexons.create_dataset(name='annotation',data=SP.array(exon_annotation))
        fexons.create_dataset(name='sample_ids',data=sample_ids)

    else:
        #carry out computation
        out_file = sys.argv[2]
        in_file  = sys.argv[1]
        #0. count total reads and mapped reads
        cmd = "%s view -c %s > %s.qa.txt" % (samtools_cmd,in_file,out_file)
        print cmd
        os.system(cmd)
        cmd = "%s view -c -F 4 %s >> %s.qa.txt" % (samtools_cmd,in_file,out_file)
        print cmd
        os.system(cmd)
        #1. run htseq qual
        cmd = "%s -o %s.qa.pdf %s" % (htseq_qa_cmd,out_file,in_file)
        print cmd
        os.system(cmd)
        #2. run htseq count (gene)
        cmd = "%s --type=gene %s %s > %s.quant.gene.txt" % (htseq_count_cmd,in_file,reference_gtf,out_file)
        print cmd
        os.system(cmd)
        #3. run htseq count (exon)
        cmd = "%s --type=exon %s %s > %s.quant.exon.txt" % (htseq_count_cmd,in_file,reference_gtf,out_file)
        print cmd
        os.system(cmd)
