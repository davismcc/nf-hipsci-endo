"""
quantify transcript expression with Salmon using lightweight alignment

USAGE:
To farm jobs to the cluster - 
python run_salmon.py cluster dir_name_in dir_name_out

To run a single job - 
python run_salmon.py in_file out_file

ARGS:
cluster - flag to run Salmon on all files in dir_name_in
dir_name_in - directory in which to find BAM files for which to quantify expression
dir_name_out - directory in which to output transcript expression values
in_file - a BAM file for which to quantify transcript expression
out_file - a directory name for the output from Salmon

BAM files from Sanger in (dir_name_in):
../data/run1/

Salmon quantities out (dir_name_out):
../data/run1/quantSalmon_lwt

e.g. python run_salmon_lightweight.py cluster ../data/run1/raw ../data/run1/quantSalmon_lwt

Davis McCarthy
November 2015
"""

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
from datetime import datetime

# settings
#---------
python_cmd = 'python'
nthreads = 8
sleep = 1
mem_thread = 16000
bam2fast_cmd = 'bam2fastq --aligned --unaligned '
aligner='STAR'
star_reference = '/nfs/research/stegle/datasets/references/human/STAR_GRCh37.75_ERCC'
transcriptome_ref = '/nfs/research/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.rel75.cdna.all.ERCC92.fa'
transcriptome_idx = '/nfs/research/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.rel75.cdna.all.ERCC92.salmon_idx'
salmon_libtype = 'IU'
salmon_cmd = '/homes/davis/davis/src/SalmonBeta-0.5.1_DebianSqueeze/bin/salmon '
salmon_cmd += 'quant -i %s -l %s ' % (transcriptome_idx, salmon_libtype)
cluster_cmd = 'bsub -n %d -q research-rh6 -R "rusage[mem=%d,tmp=50000]" -M %d -o ./cluster_out' % (nthreads, mem_thread, mem_thread)
#----------

if __name__ == '__main__':
    if 'cluster' in sys.argv:
        #run merge jobs
        dir_name_out = sys.argv[3]
        dir_name_in = sys.argv[2]
        if not os.path.exists(dir_name_out):
            os.makedirs(dir_name_out)
        FL1 = [fn for fn in glob.glob(os.path.join(dir_name_in, '*_1.fastq'))]
        FL1.sort()
        FL2 = [fn for fn in glob.glob(os.path.join(dir_name_in,'*_2.fastq'))]
        FL2.sort()
        #match sample numbers
        for fn1, fn2 in zip(FL1, FL2):
            out_file_base = os.path.splitext(os.path.basename(fn1))[0]
            out_file_base = string.replace(out_file_base, '#','_')
            out_file = os.path.join(dir_name_out, out_file_base, 'quant.sf')
            out_dir = os.path.join(dir_name_out, out_file_base)
            # if os.path.exists(out_file):
            #     print "skipping file: %s" % (out_file)
            #     continue
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            cmd = '%s %s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],
                                         fn1, fn2, out_dir)
            out_file = os.path.join(dir_name_out, '%s' % os.path.basename(fn))
            if os.path.exists(out_file):
                print "skipping file: %s" % (out_file)
                continue
            print cmd
            os.system(cmd)
            time.sleep(sleep)
    else:
        #carry out computation
        in_file_1  = sys.argv[1]
        in_file_2  = sys.argv[2]
        out_dir = sys.argv[3]
        base_name = os.path.splitext(os.path.basename(in_file_1))[0]
        base_name = os.path.splitext(base_name)[0]
        base_name = string.replace(base_name, '#','')
        tmp_dir = os.path.join('/tmp', base_name)
        salmon_tmp = os.path.join('/tmp', 'salmonQuant_light', base_name)
        if not os.path.exists(salmon_tmp):
            os.makedirs(salmon_tmp)
        #2. quantify expression
        cmd = '%s -1 %s -2 %s -o %s -p %d' % (salmon_cmd, in_file_1, in_file_2,
                                           salmon_tmp, nthreads)
        datetime.now().isoformat()
        print cmd
        os.system(cmd)
        #3. copy result
        cmd = 'mv %s/* %s' % (salmon_tmp, out_dir)
        datetime.now().isoformat()
        print cmd
        os.system(cmd)
        #3. delete tmp files
        # pdb.set_trace()
        # for root, dirs, files in os.walk(salmon_tmp, topdown=False):
        #     for name in files:
        #         os.remove(os.path.join(root, name))
        #     for name in dirs:
        #         os.rmdir(os.path.join(root, name))
        # os.rmdir(tmp_dir)
        # pass



