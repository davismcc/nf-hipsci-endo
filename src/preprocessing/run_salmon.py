"""
quantify transcript expression with Salmon using STAR alignments

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

STAR alignments in (dir_name_in):
../data/run1/alignmentsSTAR

Salmon quantities out (dir_name_out):
../data/run1/quantSalmon

e.g. python run_salmon.py cluster ../data/run1/alignmentsSTAR ../data/run1/quantSalmon

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

# settings
#---------
python_cmd = 'python'
nthreads = 16
sleep = 120
mem_thread = 65000
aligner='STAR'
star_reference = '/nfs/research/stegle/datasets/references/human/STAR_GRCh37.75_ERCC'
transcriptome_ref = '/nfs/research/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.rel75.cdna.all.ERCC92.fa'
salmon_libtype = 'ISF'
salmon_cmd = '/homes/davis/davis/src/SalmonBeta-0.5.1_DebianSqueeze/bin/salmon '
salmon_cmd += 'quant -t %s -l %s ' % (transcriptome_ref, salmon_libtype)
cluster_cmd = 'bsub -n %d -q research-rh6 -R "rusage[mem=%d,tmp=50000]" -M %d -o ./cluster_out' % (nthreads, mem_thread, mem_thread)
#----------

if __name__ == '__main__':
    if 'cluster' in sys.argv:
        #run merge jobs
        dir_name_out = sys.argv[3]
        dir_name_in = sys.argv[2]
        if not os.path.exists(dir_name_out):
            os.makedirs(dir_name_out)
        FL = glob.glob(os.path.join(dir_name_in,'*.bam'))
        FL = [fn for fn in glob.glob(os.path.join(dir_name_in,'*.bam'))
                          if not re.search('#', fn)]
        #match sample numbers
        #FL = FL[0:1]
        for fn in FL:
            out_file_base = os.path.splitext(os.path.basename(fn))[0]
            out_file = os.path.join(dir_name_out,'%s.sf' % out_file_base)
            if os.path.exists(out_file):
                print "skipping file: %s" % (out_file)
                continue
            cmd = '%s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],
                                      fn, out_file)
            print cmd
            #cmd = '%s %s %s %s' % (python_cmd, sys.argv[0],fn,out_file)
            #os.system(cmd)
            time.sleep(sleep)
    else:
        #carry out computation
        out_file = sys.argv[2]
        in_file  = sys.argv[1]
        base_name = os.path.splitext(os.path.basename(in_file))[0]
        base_name = string.replace(base_name, '#','')
        tmp_dir = os.path.join('/tmp', base_name)
        salmon_tmp = os.path.join('/tmp', 'salmonQuant', base_name)
        if not os.path.exists(salmon_tmp):
            os.makedirs(salmon_tmp)
        #1. quantify expression
        cmd = '%s -a %s -o %s -p %d' % (salmon_cmd, in_file,
                                                  salmon_tmp, nthreads)
        print cmd
        os.system(cmd)
        #2. copy result
        cmd = 'mv %s/quant.sf %s' % (salmon_tmp, out_file)
        os.system(cmd)
        print cmd
        #3. delete tmp files
        # pdb.set_trace()
        # for root, dirs, files in os.walk(salmon_tmp, topdown=False):
        #     for name in files:
        #         os.remove(os.path.join(root, name))
        #     for name in dirs:
        #         os.rmdir(os.path.join(root, name))
        # os.rmdir(tmp_dir)
        # pass



