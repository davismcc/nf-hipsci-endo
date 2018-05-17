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

python_cmd = 'python'
# settings
nthreads = 4
sleep = 120
mem_thread = 65000
aligner='STAR'
#aligner='tophat'
#aligner='gsnap'
outFilterMultimapNmax = 10
gsnap_reference = '/nfs/research/stegle/datasets/references/human/gsnap-GRHC37_ERCC'
gsnap_splice = '/nfs/research/stegle/datasets/references/human/gsnap-GRHC37_ERCC/GRCh37_splitesites'

star_reference = '/nfs/research/stegle/datasets/references/human/STAR_GRCh37.75_ERCC'
tophat_gtf = '/homes/stegle/research/datasets/references/human/tophat2-GRHC37_ERCC/Homo_sapiens.GRCh37.75_ERCC.gtf'
tophat_index = '/homes/stegle/research/datasets/references/human/tophat2-GRHC37_ERCC/trans_index/'
tophat_reference = '/homes/stegle/research/datasets/references/human/tophat2-GRHC37_ERCC/GRCh37.p13.genome.chr_only_ERCC'


cluster_cmd = 'bsub -n %d -q research-rh6 -R "rusage[mem=%d,tmp=50000]" -M %d -o ./cluster_out' % (nthreads,mem_thread,mem_thread)
#cluster_cmd = ''
bam2fast_cmd = 'bam2fastq'
samtools_cmd = 'samtools'
gsnap_cmd = 'gsnap -A sam -B 5 -t %d -n 1 -Q --nofails -d human -D %s -s %s' % (nthreads,gsnap_reference,gsnap_splice)
#STAR_cmd  = 'STAR --genomeDir %s --outFilterMultimapNmax %d --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterScoreMinOverLread 0.33 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM Unsorted --runThreadN %d' % (star_reference,outFilterMultimapNmax,nthreads)
STAR_cmd  = 'STAR --genomeDir %s --outFilterMultimapNmax %d --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterScoreMinOverLread 0.33 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --runThreadN %d' % (star_reference,outFilterMultimapNmax,nthreads)

tophat_cmd ='tophat2  -p %d --segment-length 25 --library-type=fr-unstranded --no-coverage-search --min-intron-length 6 --bowtie-n --max-multihits 20 --mate-inner-dist 350 --mate-std-dev 300 -G %s --transcriptome-index %s --rg-id sample --rg-sample sample' % (nthreads,tophat_gtf,tophat_index)

#deduplication using picard
#java -jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics#see http://gatkforums.broadinstitute.org/discussion/3891/calling-variants-in-rnaseq

if __name__ == '__main__':
    if 'cluster' in sys.argv:
        #run merge jobs
        dir_name_out = sys.argv[3]
        dir_name_in = sys.argv[2]
        if not os.path.exists(dir_name_out):
            os.makedirs(dir_name_out)
        FL = glob.glob(os.path.join(dir_name_in,'*.bam'))
        #match sample numbers
        #FL = FL[0:1]
        for fn in FL:
            out_file = os.path.join(dir_name_out,'%s' % os.path.basename(fn))
            if os.path.exists(out_file):
                print "skipping file: %s" % (out_file)
                continue
            cmd = '%s %s %s %s %s' % (cluster_cmd, python_cmd, sys.argv[0],fn,out_file)
            print cmd
            #cmd = '%s %s %s %s' % (python_cmd, sys.argv[0],fn,out_file)
            os.system(cmd)
            time.sleep(sleep)
    else:
        #carry out computation
        out_file = sys.argv[2]
        in_file  = sys.argv[1]
        base_name = os.path.splitext(os.path.basename(in_file))[0]
        base_name=string.replace(base_name,'#','')
        tmp_dir = os.path.join('/tmp',base_name)
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        fastq_file = os.path.join(tmp_dir,base_name+'#.fastq')
        fastq_file_1 = string.replace(fastq_file,'#','_1')
        fastq_file_2 = string.replace(fastq_file,'#','_2')
        sam_file_tmp = os.path.join(tmp_dir,base_name+'_tmp.sam')
        bam_file_tmp = os.path.join(tmp_dir,base_name+'_tmp.bam')
        bam_sorted_tmp = os.path.join(tmp_dir,base_name)
        #1. convert to fastq
        cmd = "%s %s --output %s" % (bam2fast_cmd, in_file, fastq_file)
        os.system(cmd)
        #2. run alignment
        if aligner=='gsnap':
            cmd = "%s %s %s > %s" % (gsnap_cmd,fastq_file_1,fastq_file_2,sam_file_tmp)
            os.system(cmd)
        elif aligner=='STAR':
            cmd = "%s --readFilesIn %s %s --outFileNamePrefix %s/" % (STAR_cmd,fastq_file_1,fastq_file_2,tmp_dir)
            os.system(cmd)
            #rename file
            cmd = 'mv %s %s' % (os.path.join(tmp_dir,'Aligned.out.sam'),sam_file_tmp)
            os.system(cmd)
        elif aligner=='tophat':
            cmd = '%s -o %s %s %s %s' % (tophat_cmd,os.path.join(tmp_dir,'Aligned.out'),tophat_reference,fastq_file_1,fastq_file_2)
            os.system(cmd)
            #merge mapped and unmapped resads
            cmd = 'mv %s %s' % (os.path.join(tmp_dir,'Aligned.out/accepted_hits.bam'),sam_file_tmp)
            os.system(cmd)
        #3. sam2bam
        cmd = "%s view -bS %s > %s" % (samtools_cmd,sam_file_tmp,bam_file_tmp)
        os.system(cmd)
        #4. sort
        cmd = "%s sort %s %s" % (samtools_cmd,bam_file_tmp,bam_sorted_tmp)
        os.system(cmd)
        #5. index
        cmd = "%s index %s.bam" % (samtools_cmd,bam_sorted_tmp)
        os.system(cmd)
        #6. copy result
        cmd = 'mv %s.bam %s' % (bam_sorted_tmp,out_file)
        os.system(cmd)
        print cmd
        cmd = 'mv %s.bam.bai %s.bai' % (bam_sorted_tmp,out_file)
        os.system(cmd)
        print cmd
     
        #3. delete tmp files
        pdb.set_trace()
        for root, dirs, files in os.walk(tmp_dir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        os.rmdir(tmp_dir)
        pass

