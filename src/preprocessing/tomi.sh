
#Create index out of reference with GSNAP version 2013-10-12 called with args: /homes/ti1/tools/gmap-2013-10-12/bin/gsnap
bsub -M 40000 /homes/ti1/tools/gmap-2012-07-20/bin/gmap_build -d mouse /nfs/research2/teichmann/tomi/references_data/genomes/Mus_musculus.GRCm38.73.dna.primary_assembly.fa -D /nfs/research2/teichmann/tomi/references_data/genomes/gmap-2012-07-20


bsub "cat /nfs/research2/teichmann/tomi/references_data/genomes/Mus_musculus.GRCm38.73.gtf | ~/tools/gmap-2012-07-20/bin/gtf_splicesites > /nfs/research2/teichmann/tomi/references_data/genomes/gmap-2012-07-20/Mus_musculus.GRCm38.73.splicesites"


bsub -q research-rh6 -n 10 -M 32000 -J "join[1-2]" "~/tools/gmap-2012-07-20/bin/gsnap -A sam -B 5 -t 10 -n 1 -Q --nofails -d mouse -D /nfs/research2/teichmann/tomi/references_data/genomes/gmap-2012-07-20/mouse -s /nfs/research2/teichmann/tomi/references_data/genomes/gmap-2012-07-20/mouse.splicesites /nfs/research2/teichmann/tomi/data/mES/bulk/10671_7_\${LSB_JOBINDEX}_1.fastq /nfs/research2/teichmann/tomi/data/mES/bulk/10671_7_\${LSB_JOBINDEX}_2.fastq > /nfs/research2/teichmann/tomi/analysis/mappings/mES/bulk/10671_7_\${LSB_JOBINDEX}.sam"


bsub -q research-rh6 -J "join[1-2]" "/homes/ti1/tools/samtools-0.1.19/samtools view -bS /nfs/research2/teichmann/tomi/analysis/mappings/mES/bulk/10671_7_\${LSB_JOBINDEX}.sam | /homes/ti1/tools/samtools-0.1.19/samtools sort - /nfs/research2/teichmann/tomi/analysis/mappings/mES/bulk/ordered/10671_7_\${LSB_JOBINDEX}.sorted"

bsub -q research-rh6 -J "join[1-2]" "/homes/ti1/tools/samtools-0.1.19/samtools view -bS /nfs/research2/teichmann/tomi/analysis/mappings/mES/bulk/10671_7_\${LSB_JOBINDEX}.sam | /homes/ti1/tools/samtools-0.1.19/samtools sort - /nfs/research2/teichmann/tomi/analysis/mappings/mES/bulk/ordered/10671_7_\${LSB_JOBINDEX}.sorted"


