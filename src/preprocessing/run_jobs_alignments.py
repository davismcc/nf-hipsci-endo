import scipy as SP
import io
import sys
import os
import pdb
import re
import time

submit_cmd = 'bsub -n 4 -q research-rh6 -R "rusage[mem=3000]" -o /homes/stegle/research/projects/1000GenomesRNASeq/scripts/cluster_out /homes/stegle/projects/1000GenomesRNASeq/scripts/mRNA/align_reads.sh'
match = re.compile('(.*) extract')
prefix = 'alignments/v1'

if __name__ == '__main__':
    #bulk downlaod files from arrayExpress
    fn = sys.argv[1]
    #get dirname of of dir
    base_dir = os.path.dirname(fn)
    fasta_dir = os.path.join(base_dir,'fasta')
    out_dir = os.path.join(base_dir,prefix)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    description = SP.loadtxt(fn,dtype='str',delimiter=',')
    L=description[1::,8]
    Lok = []
    #match unqiue library names
    library = SP.array([match.match(l).group(1) for l in L])
    library = SP.unique(library)
    library = library[50::]
    for l in library:
        print "processing: %s" % l
        bam_file = os.path.join(out_dir,'%s.bam' % l)
        if os.path.exists(bam_file):
            print "bam file exists, skipping"
            continue
        else:
            fn_1 = os.path.join(fasta_dir,'%s_1.fastq.gz' % l)
            fn_2 = os.path.join(fasta_dir,'%s_2.fastq.gz' % l)
            if not (os.path.exists(fn_1) and os.path.exists(fn_2)):
                print "fasta files missing, skipping"
		continue
                pass
        Lok.append(l)
        os.system('%s %s' %(submit_cmd,l))
	time.sleep(120)
