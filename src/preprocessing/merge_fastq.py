"""merge single cell RNA-Seq fastq files"""

import scipy as SP
import io
import sys
import os
import pdb
import re
import time
import glob

cluster_cmd = 'bsub -n 1 -q research-rh6 -R "rusage[mem=3000]" -M 3000 -o ./cluster_out'
#cluster_cmd = ''
python_cmd = 'python'
bam2fastq_cmd  = 'bam2fastq'

if __name__ == '__main__':
    if 'cluster' in sys.argv:
        #run merge jobs
        dir_name_out = sys.argv[3]
        dir_name_in = sys.argv[2]
        if not os.path.exists(dir_name_out):
            os.makedirs(dir_name_out)
        FL = glob.glob(os.path.join(dir_name_in,'*.bam'))
        #match sample numbers
        SN = [int(re.match('.*#(\d*)\.bam',fl).group(1)) for fl in FL]
        R = {}
        for i in xrange(len(FL)):
            ln = SN[i]
            if ln not in R.keys():
                R[ln] = []
            R[ln].append(FL[i])
        #submit 
        for ln in R.keys():
            files = R[ln]
            out_file = os.path.join(dir_name_out,'sc_%d.bam' % ln)
            cmd = '%s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],out_file)
            for fn in files:
                cmd += ' %s' % fn
            os.system(cmd)
            #print cmd
    else:
        #carry out computation
        out_file = sys.argv[1]
        merge_files = sys.argv[2::]
        #1. create individual fastq files
        cmd = 'samtools merge %s' % (out_file)
        for mf in merge_files:
            cmd += ' %s' % mf
        os.system(cmd)
        pass

