"""convert cram single cell RNA-Seq to bam and fastq files

USAGE:
To farm jobs to the cluster - 
python run_salmon.py cluster dir_name_in dir_name_out

To run a single job - 
python run_salmon.py in_file out_file

ARGS:
cluster - flag to run conversion on all files in dir_name_in
dir_name_in - directory in which to find BAM files for which to quantify expression
dir_name_out - directory in which to output transcript expression values
in_file - a BAM file for which to quantify transcript expression
out_file - a directory name for the output from Salmon
reconvert - optional flag as the final argument to force new conversion of existing files

BAM files from Sanger in (dir_name_in):
../data/run1/raw

Salmon quantities out (dir_name_out):
../data/run1/raw

e.g. python convert_cram.py cluster ../data/run1/raw ../data/run1/raw
e.g. python convert_cram.py cluster ../data/run1/raw ../data/run1/raw reconvert

Davis McCarthy & Florian Buettner
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

# Settings
#-----------
python_cmd = 'python'
nthreads = 4
sleep = 3
mem_thread = 65000
bam2fast_cmd = 'bam2fastq' # extract aligned and unaligned reads filtered s.t.
# reads failing QC will not be extracted(?)
cram2bam_cmd = 'java -jar /homes/buettner/research/users/buettner/hipsci-singlecell/preprocessing/CRAM/cramtools-3.0.jar bam'
cram2fast_cmd = 'java -jar /homes/buettner/research/users/buettner/hipsci-singlecell/preprocessing/CRAM/cramtools-3.0.jar fastq'
samtools_cmd = 'samtools'
cluster_cmd = 'bsub -n %d -q research-rh6 -R "rusage[mem=%d,tmp=50000]" -M %d -o ./cluster_out' % (nthreads, mem_thread, mem_thread)
#-----------

if __name__ == '__main__':
    if 'cluster' in sys.argv:
        #run merge jobs
        dir_name_out = sys.argv[3]
        dir_name_in = sys.argv[2]
        if not os.path.exists(dir_name_out):
            os.makedirs(dir_name_out)
        FL = glob.glob(os.path.join(dir_name_in, '*.cram'))
        #match sample numbers
        #FL = FL[0:1]
        for fn in FL:
            out_file_base = os.path.join(dir_name_out,
                                         '%s' % os.path.basename(fn))
            out_file_base = os.path.splitext(out_file_base)[0]
            skip_existing = (os.path.exists('%s.bam' % out_file_base) and not
                         'reconvert' in sys.argv)
            if skip_existing:
                print "skipping file conversion for: %s" % (fn)
                continue
            cmd = '%s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0], fn,
                                      out_file_base)
            print cmd
            #cmd = '%s %s %s %s' % (python_cmd, sys.argv[0], fn, out_file)
            os.system(cmd)
            time.sleep(sleep)
    else:
        #carry out computation
        out_file_base = sys.argv[2]
        in_file  = sys.argv[1]
        base_name = os.path.splitext(os.path.basename(in_file))[0]
        base_name = string.replace(base_name,'#','_')
        tmp_dir = os.path.join('/tmp', base_name)
        print 'temp dir is: %s' % tmp_dir
        print datetime.now().isoformat()
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        out_file_base_tmp = os.path.join(tmp_dir, base_name)
        fastq_base_tmp = string.replace(out_file_base_tmp, '#', '_')
        #1. convert to bam
        cmd = "%s -I %s -O %s.bam" % (cram2bam_cmd, in_file, out_file_base_tmp)
        print cmd
        print datetime.now().isoformat()
        os.system(cmd)
        #2. convert to fastq
        #cmd = "%s -I %s -F %s" % (cram2fast_cmd, in_file, out_file_base_tmp)
        cmd = "%s %s --output %s --overwrite" % (bam2fast_cmd,
                                     out_file_base_tmp + '.bam',
                                     fastq_base_tmp + '#.fastq')
        print cmd
        print datetime.now().isoformat()
        os.system(cmd)
        #3. copy result
        cmd = 'mv %s.bam %s.bam' % (out_file_base_tmp, out_file_base)
        print cmd
        print datetime.now().isoformat()
        os.system(cmd)
        #pdb.set_trace()
        cmd = 'mv %s_1.fastq %s_1.fastq' % (fastq_base_tmp, out_file_base)
        print cmd
        print datetime.now().isoformat()
        os.system(cmd)
        cmd = 'mv %s_2.fastq %s_2.fastq' % (fastq_base_tmp, out_file_base)
        print cmd
        print datetime.now().isoformat()
        os.system(cmd)
        #4. delete tmp files
        for root, dirs, files in os.walk(tmp_dir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(tmp_dir)
        pass

